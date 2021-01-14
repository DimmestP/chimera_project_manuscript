"""
    Encodes fitting and interpolating using Gaussian processes following Rasmussen and Williams (2006).

    The Gaussian process algorithms come from chapter 3 and the Jacobian of the negative log likelihood from chapter 5 of Rasmussen and Williams.

    Covariance functions can either be linear, squared exponential, neural network-like, or squared exponential with a linear trend. Bounds for hyperparameters are specified in log10 space. Hyperparameters are given in log space.

    A typical workflow is:

    g= gp.maternGP({0: (-4, 4), 1: (-4, 4), 2: (-4, -2)}, x, y)
    g.findhyperparameters()
    g.results()
    g.predict(x, derivs= 1)
    plt.figure()
    plt.subplot(2,1,1)
    g.sketch('.')
    plt.subplot(2,1,2)
    g.sketch('.', derivs= 1)
    plt.show()

    The mean of the fitted GP is g.f and its variance is g.fvar.

    Running g.predict(x, derivs= 2) returns g.df, g.dfvar, g.ddf, and g.ddfvar.

    Prior functions can also be sampled. For example,

    g= gp.sqexplinGP({0: (-2,2), 1: (-2,2), 2: (-2,2), 3: (-2,2), 4: (-2,2)}, x, y)
    plot(x, g.sampleprior(3, th= [1.0, 0.1, 3.1, 1.3]))

    will plot three samples of the prior latent functions with hyperparameters 1.0, 0.1, 3.1, and 1.3. There is no need to specify the hyperparameter for measurement error: it is not used to generate prior functions.

    N.B. small (close to zero) values of the estimated measurement error can lead to instabilities in finding the hyperparameters.
"""

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt


class gaussianprocess:
    def __init__(self, lthbounds, x, y, merrors= False, warnings= True):
        '''
        Creates a Gaussian process.

        Arguments
        --
        lthbounds: a dictionary of pairs of the bounds on the hyperparameters in log10 space,
        such as {0: [0,6], 1: [-3,4], 2: [-6,-4]}
        x: a 1-d array of the abscissa data
        y: a multi-dimensional array of the ordinate data
        merrors: if specified, a 1-d array of the measurement errors (as variances)
        v'''
        x, y= np.asarray(x), np.asarray(y)
        if x.size == 0 or y.size == 0:
            raise gaussianprocessException('x or y is incorrectly initialised')
        else:
            if warnings and len(x) > 1.0e3:
                print('GP Warning: large data set -', len(x), 'datapoints - may slow optimisation and give memory issues')
            self.b= [lthbounds[a] for a in lthbounds.keys()]
            self.x, self.y, self.xnew= x, y, x
            self.merrors= merrors
            self.version= '1.06'





    def covfn(self):
        raise NotImplementedError(' No covariance function specified in class %s'
                                  % self.__class__.__name__)

    def d1covfn(self):
        raise NotImplementedError(' No first derivative of the covariance function specified in class %s' % self.__class__.__name__)

    def d1d2covfn(self):
        raise NotImplementedError(' No second derivative of the covariance function specified in class %s' % self.__class__.__name__)


    def kernelmatrix(self, lth, x):
        """
        Returns kernel matrix K(X,X) supplemented with measurement noise and its Cholesky decomposition.

        Arguments
        --
        lth: log of the hyperparameters
        x: abscissa values
        merrors: if specified, rescales the fitted measurement error
        """
        k= np.empty((len(x), len(x)))
        for i in range(len(x)):
            k[i,:]= self.covfn(x[i], x, lth)[0]
        if np.any(self.merrors):
            kn= k + np.exp(lth[-1])*np.diag(self.merrors)
        else:
            kn= k + np.exp(lth[-1])*np.identity(len(x))
        try:
            L= linalg.cho_factor(kn)
        except np.linalg.LinAlgError:
            raise gaussianprocessException('Kernel matrix is not positive definite')
        return k, L


    def nlml(self, lth):
        """
        Returns negative of log marginal likelihood.

        Arguments
        --
        lth: log of the hyperparameters
        """
        x, y= self.x, self.y
        k, L= self.kernelmatrix(lth, x)
        try:
            al= linalg.cho_solve(L, y)
        except:
            raise gaussianprocessException('Evaluating nlml failed')
        halfdetK= np.sum(np.log(np.diagonal(L[0])))
        return 0.5*np.dot(y, al) + halfdetK + 0.5*len(y)*np.log(2*np.pi)


    def jacnlml(self, lth):
        """
        Returns the Jacobian of negative log marginal likelihood with respect to the hyperparameters with deriviatives being taken assuming the hyperparmaters are in log space.

        Arguments
        --
        lth: log of the hyperparameters
        """
        x, y= self.x, self.y
        k, L= self.kernelmatrix(lth, x)
        # find derivatives of kernel matrix wrt hyperparameters
        kjac= np.empty((len(x), len(x), len(lth)))
        for i in range(len(x)):
            kjac[i, :, :-1]= self.covfn(x[i], x, lth)[1]
        if np.any(self.merrors):
            kjac[:, :, -1]= np.diag(self.merrors)*np.exp(lth[-1])
        else:
            kjac[:, :, -1]= np.identity(len(x))*np.exp(lth[-1])
        # calculate jacobian
        al= linalg.cho_solve(L, y)
        alal= np.outer(al, al)
        Kinv= linalg.cho_solve(L, np.identity(len(x)))
        return np.asarray([-0.5*np.trace(np.dot(alal - Kinv, kjac[:,:,i])) for i in range(len(lth))])


    def findhyperparameters(self, noruns= 1, noinits= 100, exitearly= False, stvals= False, optmethod= 'l_bfgs_b',
                            optmessages= False, quiet= True, linalgmax= 3):
        """
        Finds the best fit hyperparameters (.lth_opt) and the optimum value of negative log marginal likelihood (.nlml_opt).

        Arguments
        --
        noruns: number of attempts to find optimum hyperparameters (the best of all the runs is chosen)
        noinits: number of attempts to find a good initial condition for the optimization
        exitearly: if True, fitting stops at the first successful attempt
        stvals: an (optional) initial guess for the log hyperparameters
        optmethod: the optimization routine to be used, either 'l_bfgs_b' (default) or 'tnc'
        optmessages: if True, display messages from the optimization routine
        quiet: if True, print warning that if an optimum hyperparameter is at a bound
        linalgmax: number of attempts (default is 3) if a linear algebra (numerical) error is generated
        """
        b= self.b
        self.hparamerr= []
        lmlml= np.empty(noruns)
        lthf= np.empty((noruns, len(b)))
        success= np.empty(noruns)
        # convert b into exponential base
        b= np.array(b)*np.log(10)
        # run optimization
        for i in range(noruns):
            # find initial conditions
            bestnlml= np.inf
            for k in range(noinits):
                lth= [np.random.uniform(b[j][0], b[j][1]) for j in range(len(b))]
                try:
                    trialnlml= self.nlml(lth)
                    if trialnlml < bestnlml:
                        bestnlml= trialnlml
                        stvals= lth
                except gaussianprocessException:
                    pass
            # minimize
            linalgerror= 0
            while linalgerror < linalgmax:
                try:
                    if np.any(stvals):
                        # initial values given for hyperparameters
                        lth= stvals
                    else:
                        # choose random initial values for hyperparameters
                        lth= [np.random.uniform(b[j][0], b[j][1]) for j in range(len(b))]
                    # run Gaussian process
                    if optmethod == 'tnc':
                        from scipy.optimize import fmin_tnc
                        lthf[i,:], nf, success[i]= fmin_tnc(self.nlml, lth, fprime= self.jacnlml,
                                                          bounds= b, maxfun= 1000, messages= optmessages)
                        lmlml[i]= self.nlml(lthf[i,:])
                        break
                    elif optmethod == 'l_bfgs_b':
                        from scipy.optimize import fmin_l_bfgs_b
                        lthf[i,:], lmlml[i], dout= fmin_l_bfgs_b(self.nlml, lth, fprime= self.jacnlml,
                                                                bounds= b, disp= optmessages)
                        success[i]= dout['warnflag'] + 1
                        break
                    else:
                        print(optmethod + ' unrecognized.')
                except gaussianprocessException:
                    print(' Warning: linear algebra error - trying a different initial condition')
                    stvals= False
                    linalgerror += 1
            if success[i] != 1 or np.any(np.isnan(lthf)):
                print(' Warning: optimization of hyperparameters failed at run ' + str(i+1))
                # stop optimising initial condition
                stvals= False
                noinits= 0
            else:
                if exitearly: break
        # only process runs that did not converge
        if np.any(success == 1):
            lmlml= lmlml[success == 1]
            lthf= lthf[success == 1]
            # find best choice
            lthb= lthf[lmlml.argmin()]
            self.nlml_opt= lmlml.min()
            # print warning
            for i in range(len(b)):
                if (lthb[i] == b[i][1] or lthb[i] == b[i][0]):
                    if not quiet:
                        print( ' Warning: hparam[' + str(i) + '] is at a boundary')
                        print('\thparam[' + str(i) + ']= {:e}'.format(np.exp(lthb[i]))
                          + ' [{:e}'.format(np.exp(b[i][0])) + ', {:e}]'.format(np.exp(b[i][1])))
                    if lthb[i] == b[i][1]:
                        self.hparamerr.append([i, 'u'])
                    else:
                        self.hparamerr.append([i, 'l'])
            self.lth_opt= lthb
        else:
            raise gaussianprocessException('Optimization of hyperparameters failed')




    def results(self, warning= True):
        '''
        Displays results from optimizing hyperparameters.

        Arguments
        --
        warning: if True, warn when a hyperparameter hits a boundary
        '''
        print('log(max likelihood)= %e' % (-self.nlml_opt))
        for j, pv in enumerate(np.exp(self.lth_opt)):
            print('hparam[' + str(j) + ']= {:e}'.format(pv) + ' [{:e}'.format(10**(self.b[j][0]))
                + ', {:e}]'.format(10**(self.b[j][1])))
        if warning:
            for el in self.hparamerr:
                if el[1] == 'l':
                    print('Warning: hyperparameter ' + str(el[0]) + ' is at a lower bound')
                else:
                    print('Warning: hyperparameter ' + str(el[0]) + ' is at an upper bound')



    def sample(self, nosamples= 1, derivs= 0, xnew= False):
        '''
        Generate samples from the Gaussian process as an array.

        Arguments
        --
        nosamples: number of samples
        derivs: if 0, only the latent function is sampled; if 1, the latent function and the first derivative are sampled; if 2, the latent function and the first and second derivatives are sampled
        xnew: (default: False) potential new abscissa values for sampled ordinate values – g.predict() will be called
        '''
        if np.any(xnew):
            self.predict(xnew, derivs= derivs)
        else:
            try:
                xnew= self.xnew
                ss= np.transpose(np.random.multivariate_normal(self.mnp, self.covp, nosamples))
                if derivs == 0:
                    return ss
                elif derivs == 1:
                    return ss[:len(xnew),:], ss[len(xnew):2*len(xnew),:]
                elif derivs == 2:
                    return ss[:len(xnew),:], ss[len(xnew):2*len(xnew),:], ss[2*len(xnew):,:]
            except AttributeError:
                print( ' Run gp.predict() first before sampling.')




    def sampleprior(self, size= 1, lth= False):
        '''
        Generate samples from the prior of the Gaussian process as an array.

        Arguments
        --
        size: number of samples
        lth: log hyperparameters to use (if not specified, the hyperparameters are chosen at random)
        '''
        x, y, b= self.x, self.y, self.b
        if np.any(lth):
            # hyperparameters are given (measurement error is not necessary)
            if len(lth) == self.noparams:
                lth= np.concatenate((lth, [1.0]))
        else:
            # sample random hyperparameters
            lth= np.log(np.power(10, [np.random.uniform(b[i][0], b[i][1]) for i in range(len(b))]))
        cov= self.kernelmatrix(lth, x)[0]
        return np.transpose(np.random.multivariate_normal(np.zeros(len(x)), cov, size))




    def predict(self, xnew, merrorsnew= False, derivs= 0, addnoise= False):
        """
        Determines the predicted mean latent function (.f) and its variance (.fvar) and potentially the predicted mean first derivative (.df) and its variance (.dfvar) and the predicted mean second derivative (.ddf) and its variance (.ddfvar) . Also .mnp is the predicted combined array of the mean latent function and its mean derivatives and .covp is the corresponding covariance matrix.

        Arguments
        --
        xnew: abscissa values for which predicted ordinate values are desired
        merrorsnew: if specified, the expected measurements errors at xnew (need not be specified if xnew= x)
        derivs: if 0, only the latent function is inferred; if 1, the latent function and the first derivative are inferred; if 2, the latent function and the first and second derivatives are inferred
        addnoise: if True, add measuremnet noise to the predicted variance
        """
        xnew= np.asarray(xnew)
        if len(self.x) == len(xnew) and (self.x == xnew).all():
            xold= True
        else:
            xold= False
        if np.any(self.merrors) and not np.any(merrorsnew) and not xold:
            print('Length of xnew is different from x.')
            raise gaussianprocessException('Measurement errors were used to find the hyperparameters and measurement errors are therefore required for any predictions.')
        elif not hasattr(self, 'lth_opt'):
            raise gaussianprocessException(' Run gp.findhyperparameters() first before making predictions.')
        else:
            # set up
            self.xnew= xnew
            lth, x, y= self.lth_opt, self.x, self.y
            # work with an array of length 3*N: the first N values being the function,
            # the second N values being the first derivative, and the last N values being the second derivative
            Knewold= np.empty((len(xnew), len(x)))
            Knewnew= np.empty((len(xnew), len(xnew)))
            if derivs > 0:
                d1Knewold= np.empty((len(xnew), len(x)))
                d1Knewnew= np.empty((len(xnew), len(xnew)))
                d1d2Knewnew= np.empty((len(xnew), len(xnew)))
            if derivs > 1:
                d12Knewold= np.empty((len(xnew), len(x)))
                d12Knewnew= np.empty((len(xnew), len(xnew)))
                d12d2Knewnew= np.empty((len(xnew), len(xnew)))
                d12d22Knewnew= np.empty((len(xnew), len(xnew)))
            for i in range(len(xnew)):
                Knewold[i,:]= self.covfn(xnew[i], x, lth)[0]
                Knewnew[i,:]= self.covfn(xnew[i], xnew, lth)[0]
                if derivs > 0:
                    d1Knewold[i,:]= self.d1covfn(xnew[i], x, lth)[0]
                    d1Knewnew[i,:]= self.d1covfn(xnew[i], xnew, lth)[0]
                    d1d2Knewnew[i,:]= self.d1d2covfn(xnew[i], xnew, lth)[0]
                if derivs > 1:
                    d12Knewold[i,:]= self.d12covfn(xnew[i], x, lth)[0]
                    d12Knewnew[i,:]= self.d12covfn(xnew[i], xnew, lth)[0]
                    d12d2Knewnew[i,:]= self.d12d2covfn(xnew[i], xnew, lth)[0]
                    d12d22Knewnew[i,:]= self.d12d22covfn(xnew[i], xnew, lth)[0]
            if derivs == 0:
                kv= Knewold
                km= Knewnew
            elif derivs == 1:
                kv= np.bmat([[Knewold], [d1Knewold]])
                km= np.bmat([[Knewnew, np.transpose(d1Knewnew)],
                             [d1Knewnew, d1d2Knewnew]])
            elif derivs == 2:
                kv= np.bmat([[Knewold], [d1Knewold], [d12Knewold]])
                km= np.bmat([[Knewnew, np.transpose(d1Knewnew), np.transpose(d12Knewnew)],
                              [d1Knewnew, d1d2Knewnew, np.transpose(d12d2Knewnew)],
                              [d12Knewnew, d12d2Knewnew, d12d22Knewnew]])
            # find mean prediction
            k, L= self.kernelmatrix(lth, x)
            m= np.dot(kv, linalg.cho_solve(L, y))
            mnp= np.reshape(np.array(m), np.size(m))
            self.mnp= mnp
            # find variance of prediction
            covp= km - np.dot(kv, linalg.cho_solve(L, np.transpose(kv)))
            self.covp= covp
            varp= np.diag(covp)
            # for user
            self.f= mnp[:len(xnew)]
            self.fvar= varp[:len(xnew)]
            fvar= varp[:len(xnew)]
            if addnoise:
                # add measurement error to the variance of the latent function
                if np.any(self.merrors):
                    if xold:
                        self.fvar= fvar + np.exp(lth[-1])*np.diag(self.merrors)
                    else:
                        self.fvar= fvar + merrorsnew
                else:
                    self.fvar= fvar + np.exp(lth[-1])*np.identity(len(xnew))
            else:
                # just take the variance of the latent function
                self.fvar= fvar
            if derivs > 0:
                self.df= mnp[len(xnew):2*len(xnew)]
                self.dfvar= varp[len(xnew):2*len(xnew)]
            if derivs > 1:
                self.ddf= mnp[2*len(xnew):]
                self.ddfvar= varp[2*len(xnew):]




    def batchpredict(self, xnew, merrorsnew= False, derivs= 0, addnoise= False, maxlen= 1000):
        '''
        Run batch predictions of size maxlen to increase speed and reduce memory issues.

        Note that you must run predict not batchpredict if you wish to generate samples from the Gaussian process.

        Arguments
        --
        xnew: abscissa values for which predicted ordinate values are desired
        merrorsnew: if specified, the expected measurements errors at xnew (need not be specified if xnew= x)
        derivs: if 0, only the latent function is inferred; if 1, the latent function and the first derivative are inferred; if 2, the latent function and the first and second derivatives are inferred
        addnoise: if True, add measuremnet noise to the predicted variance
        maxlen: number of data points to process in batches (default: 1000)
        '''
        if len(xnew) < maxlen:
            self.predict(xnew, merrorsnew, derivs, addnoise)
        else:
            f= np.zeros(len(xnew))
            fvar= np.zeros(len(xnew))
            if derivs > 0:
                df= np.zeros(len(xnew))
                dfvar= np.zeros(len(xnew))
            if derivs > 1:
                ddf= np.zeros(len(xnew))
                ddfvar= np.zeros(len(xnew))
            bs= np.append(np.arange(0, len(xnew), maxlen), len(xnew) + 1)
            for i in range(1, len(bs)):
                self.predict(xnew[bs[i-1]:bs[i]], merrorsnew, derivs, addnoise)
                f[bs[i-1]:bs[i]]= self.f
                fvar[bs[i-1]:bs[i]]= self.fvar
                if derivs > 0:
                    df[bs[i-1]:bs[i]]= self.df
                    dfvar[bs[i-1]:bs[i]]= self.dfvar
                if derivs > 1:
                    ddf[bs[i-1]:bs[i]]= self.ddf
                    ddfvar[bs[i-1]:bs[i]]= self.ddfvar
            self.f= f
            self.fvar= fvar
            if derivs > 0:
                self.df= df
                self.dfvar= dfvar
            if derivs > 1:
                self.ddf= ddf
                self.ddfvar= ddfvar





    def sketch(self, datasymbol= 'o', derivs= 0, GPcolor= 'blue', nostds= 2):
        """
        Plots data with mean prediction plus band of twice the standard deviation.

        Arguments
        --
        datasymbol: the symbol used to mark the data points (if False do not plot data points)
        derivs: if 0, plot data and mean prediction; if 1, plot first derivative with respect to x; if 2, plot second derivative
        GPcolor: color to draw mean and standard deviation of Gaussian process
        nostds: number of standard deviations to use as errorbars
        """
        x, y, xnew= self.x, self.y, self.xnew
        if derivs == 0:
            f= self.f
            sd= np.sqrt(self.fvar)
            if datasymbol: plt.plot(x, y, 'r' + datasymbol)
        elif derivs == 1:
            f= self.df
            sd= np.sqrt(self.dfvar)
        elif derivs == 2:
            f= self.ddf
            sd= np.sqrt(self.ddfvar)
        else:
            print('sketch: error in derivs')
        plt.plot(xnew, f, color= GPcolor)
        plt.fill_between(xnew, f-nostds*sd, f+nostds*sd, facecolor= GPcolor, alpha=0.2)



####

class lnGP(gaussianprocess):
    '''
    Gaussian process with a linear covariance function
    '''
    noparams= 2
    description= 'linear Gaussian process'

    def covfn(self, x, xp, lth):
        '''
        Linear covariance function.
        Returns the kernel function and Jacobian of the kernel function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        k= th[0] + th[1]*x*xp
        jk= np.empty((len(xp), self.noparams))
        jk[:,0]= th[0]*np.ones(len(xp))
        jk[:,1]= th[1]*x*xp
        return k, jk

    def gradcovfn(self, x, xp, lth):
        '''
        Returns the gradient of the covariance function with respect to x.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        return th[1]*xp, False

    def hesscovfn(self, x, xp, lth):
        '''
        Returns the Hessian of the covariance function with respect to x and xp.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        return th[1]*np.ones(len(xp)), False


####

class nnGP(gaussianprocess):
    '''
    Gaussian process with a neural network covariance function.
    '''
    noparams= 2
    description= 'neural network Gaussian process'

    def info(self):
        print('hparam[0] determines the initial value')
        print('hparam[1] determines the flexibility')
        print('hparam[2] determines the variance of the measurement error')

    def covfn(self, x, xp, lth):
        """
        Neural network covariance function from Williams.
        Returns the kernel function and Jacobian of the kernel function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        """
        th= np.exp(lth)
        k= (np.arcsin(2*(th[0] + x*xp*th[1])/np.sqrt(1+2*(th[0]+x**2*th[1]))
                      /np.sqrt(1+2*(th[0]+xp**2*th[1]))))*2/np.pi
        jk= np.empty((len(xp), self.noparams))
        den= np.pi*(1+2*th[0]+2*th[1]*x**2)*(1+2*th[0]+2*th[1]*xp**2) \
          *np.sqrt(1+4*th[0]*(1+th[1]*(x-xp)**2)+2*th[1]*(x**2+xp**2))
        jk[:,0]= (4*(1+2*th[0]*(1+th[1]*(x-xp)**2) - 2*th[1]**2*x*(x-xp)**2*xp \
                     + 2*th[1]*(x**2-x*xp+xp**2)))/den*th[0]
        jk[:,1]= -(4*(2*th[0]**2*(x-xp)**2 - x*xp*(1+th[1]*(x**2+xp**2)) \
                      + th[0]*(-2*th[1]*x**3*xp+xp**2-2*x*xp*(2+th[1]*xp**2) \
                               +x**2*(1+4*th[1]*xp**2))))/den*th[1]
        return k, jk

    def d1covfn(self, x, xp, lth):
        '''
        Returns the d/dx of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        return 4*th[1]*(-2*th[0]*x+xp+2*th[0]*xp)/(1+2*th[0]+2*th[1]*x**2) \
          /np.sqrt(1+4*th[0]*(1+th[1]*(x-xp)**2)+2*th[1]*(x**2+xp**2))/np.pi, False

    def d1d2covfn(self, x, xp, lth):
        '''
        Returns the d/dx d/dxp of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        return 4*(th[1] + 4*th[0]*th[1])/(1 + 4*th[0]*(1+th[1]*(x-xp)**2) \
                                          + 2*th[1]*(x**2+xp**2))**1.5/np.pi, False

    def d12covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        return -8*th[1]*(8*th[0]**3 + th[0]**2*(6 - 8*th[1]*x*(x - 3*xp) - 16*th[1]**2*x*(x - xp)**3) + th[1]*x*xp*(3 + 6*th[1]*x**2 + 4*th[1]*xp**2) + th[0]*(1 - 2*th[1]*x*(x - 9*xp) - 8*th[1]**2*x*(x**3 - 3*x**2*xp + 3*x*xp**2 - 2*xp**3)))/(np.pi*(1 + 2*th[0] + 2*th[1]*x**2)**2*(1 + 4*th[0]*(1 + th[1]*(x - xp)**2) + 2*th[1]*(x**2 + xp**2))**1.5), False

    def d12d2covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 d/dxp of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        return -24*(1 + 4*th[0])*th[1]**2*(x + 2*th[0]*x - 2*th[0]*xp)/(np.pi*(1 + 4*th[0]*(1 + th[1]*(x - xp)**2) + 2*th[1]*(x**2 + xp**2))**2.5), False

    def d12d22covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 d^2/dxp^2 of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        return -48*(1 + 4*th[0])*th[1]**2*(4*th[0]**2*(-1 + 4*th[1]*(x - xp)**2) - 5*th[1]*x*xp + th[0]*(-1 + 4*th[1]*(2*x**2 - 5*x*xp + 2*xp**2)))/(np.pi*(1 + 4*th[0]*(1 + th[1]*(x - xp)**2) + 2*th[1]*(x**2 + xp**2))**3.5), False



####


class sqexpGP(gaussianprocess):
    '''
    Gaussian process with a squared exponential covariance function.
    '''
    noparams= 2
    description= 'squared exponential Gaussian process'

    def info(self):
        print('hparam[0] determines the amplitude of variation')
        print('hparam[1] determines the flexibility')
        print('hparam[2] determines the variance of the measurement error')

    def covfn(self, x, xp, lth):
        '''
        Squared exponential covariance function.
        Returns the kernel function and Jacobian of the kernel function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        xp= np.array(xp)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        k= th[0]*e
        jk= np.empty((len(xp), self.noparams))
        jk[:,0]= e*th[0]
        jk[:,1]= -th[0]*th[1]*e/2.0*(x-xp)**2
        return k, jk

    def d1covfn(self, x, xp, lth):
        '''
        Returns d/dx of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return -e*th[0]*th[1]*(x-xp), False

    def d1d2covfn(self, x, xp, lth):
        '''
        Returns d/dx d/dxp of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return -e*th[0]*th[1]*(-1 + th[1]*(x-xp)**2), False

    def d12covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return e*th[0]*th[1]*(-1 + th[1]*(x - xp)**2), False

    def d12d2covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 d/dxp of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return e*th[0]*th[1]**2*(-3 + th[1]*(x - xp)**2)*(x - xp), False

    def d12d22covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 d^2/dxp^2 of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return e*th[0]*th[1]**2*(3 - 6*th[1]*(x - xp)**2 + th[1]**2*(x - xp)**4), False

####


class sqexplinGP(gaussianprocess):
    '''
    Gaussian process with a squared exponential covariance function with a linear trend.
    Returns the kernel function and Jacobian of the kernel function.

    Arguments
    --
    x: a 1-d array of abscissa
    xp: a 1-d array of alternative abscissa
    lth: the log of the hyperparameters
    '''
    noparams= 3
    description= 'squared exponential Gaussian process with a linear trend'

    def info(self):
        print('hparam[0] determines the amplitude of variation')
        print('hparam[1] determines the flexibility')
        print('hparam[2] determines the linear trend with increasing input')
        print('hparam[3] determines the variance of the measurement error')

    def covfn(self, x, xp, lth):
        '''
        Squared exponential covariance function with a linear trend.
        Returns the kernel function and Jacobian of the kernel function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        xp= np.array(xp)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        k= th[0]*e + th[2]*x*xp
        jk= np.empty((len(xp), self.noparams))
        jk[:,0]= e*th[0]
        jk[:,1]= -th[0]*th[1]*e/2.0*(x-xp)**2
        jk[:,2]= x*xp*th[2]
        return k, jk


    def d1covfn(self, x, xp, lth):
        '''
        Returns d/dx of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return -e*th[0]*th[1]*(x-xp) + th[2]*xp, False


    def d1d2covfn(self, x, xp, lth):
        '''
        Returns d/dx d/xp of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return th[2] - e*th[0]*th[1]*(-1 + th[1]*(x-xp)**2), False

    def d12covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return e*th[0]*th[1]*(-1 + th[1]*(x - xp)**2), False

    def d12d2covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 d/dxp of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return e*th[0]*th[1]**2*(-3 + th[1]*(x - xp)**2)*(x - xp), False

    def d12d22covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 d^2/dxp^2 of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return e*th[0]*th[1]**2*(3 - 6*th[1]*(x - xp)**2 + th[1]**2*(x - xp)**4), False

####


class maternGP(gaussianprocess):
    '''
    Gaussian process with a Matern covariance function that is twice differentiable.
    Returns the kernel function and Jacobian of the kernel function.

    Arguments
    --
    x: a 1-d array of abscissa
    xp: a 1-d array of alternative abscissa
    lth: the log of the hyperparameters
    '''
    noparams= 2
    description= '(twice differentiable) Matern covariance function'

    def info(self):
        print('hparam[0] determines the amplitude of variation')
        print('hparam[1] determines the stiffness')
        print('hparam[2] determines the variance of the measurement error')

    def covfn(self, x, xp, lth):
        '''
        Squared exponential covariance function.
        Returns the kernel function and Jacobian of the kernel function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        xp= np.asarray(xp)
        r= np.abs(x - xp)
        s5= np.sqrt(5)
        e= np.exp(-s5*r/th[1])
        poly= 1 + 5*r**2/3/th[1]**2 + s5*r/th[1]
        k= th[0]*e*poly
        jk= np.empty((len(xp), self.noparams))
        jk[:,0]= e*poly
        jk[:,1]= 5*e*th[0]*r**2*(th[1] + s5*r)/3/th[1]**4
        return k, jk

    def d1covfn(self, x, xp, lth):
        '''
        Returns d/dx of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        r= np.abs(x - xp)
        s5= np.sqrt(5)
        e= np.exp(-s5*r/th[1])
        df= 5*e*th[0]*r*(th[1] + s5*r)/3/th[1]**3
        sns= np.ones(np.size(xp))
        sns[x > xp]= -1
        return sns*df, False

    def d1d2covfn(self, x, xp, lth):
        '''
        Returns d/dx d/dxp of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        r= np.abs(x - xp)
        s5= np.sqrt(5)
        e= np.exp(-s5*r/th[1])
        return 5*e*th[0]*(th[1]**2 + s5*th[1]*r - 5*r**2)/3/th[1]**4, False

    def d12covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        r= np.abs(x - xp)
        s5= np.sqrt(5)
        e= np.exp(-s5*r/th[1])
        return -5*e*th[0]*(th[1]**2 + s5*th[1]*r - 5*r**2)/3/th[1]**4, False

    def d12d2covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 d/dxp of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        r= np.abs(x - xp)
        s5= np.sqrt(5)
        e= np.exp(-s5*r/th[1])
        df= 25*e*th[0]*r*(3*th[1] - s5*r)/3/th[1]**5
        sns= np.ones(np.size(xp))
        sns[x > xp]= -1
        return sns*df, False

    def d12d22covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 d^2/dxp^2 of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        r= np.abs(x - xp)
        s5= np.sqrt(5)
        e= np.exp(-s5*r/th[1])
        return 25*e*th[0]*(3*th[1]**2 + 5*r**2 - 5*s5*th[1]*r)/3/th[1]**6, False


####
class periodicGP(gaussianprocess):
    '''
    Gaussian process with a periodic covariance function.
    Returns the kernel function and Jacobian of the kernel function.

    Arguments
    --
    x: a 1-d array of abscissa
    xp: a 1-d array of alternative abscissa
    lth: the log of the hyperparameters
    '''
    noparams= 3
    description= 'periodic covariance function'

    def info(self):
        print('hparam[0] determines the amplitude of variation')
        print('hparam[1] determines the stiffness')
        print('hparam[2] determines the period')

    def covfn(self, x, xp, lth):
        '''
        Periodic covariance function.
        Returns the kernel function and Jacobian of the kernel function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        xp= np.asarray(xp)
        pr= np.exp(-2*th[1]*np.sin((x-xp)/th[2])**2)
        k= th[0]*pr
        jk= np.empty((len(xp), self.noparams))
        jk[:,0]= pr
        jk[:,1]= -2*pr*th[0]*np.sin((x-xp)/th[2])**2
        jk[:,2]= 2*pr*th[0]*th[1]*(x-xp)*np.sin(2*(x-xp)/th[2])/th[2]**2
        return k, jk

    def d1covfn(self, x, xp, lth):
        '''
        Returns d/dx of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        pr= np.exp(-2*th[1]*np.sin((x-xp)/th[2])**2)
        return -2*pr*th[0]*th[1]*np.sin(2*(x-xp)/th[2])/th[2], False

    def d1d2covfn(self, x, xp, lth):
        '''
        Returns d/dx d/xp of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        pr= np.exp(-2*th[1]*np.sin((x-xp)/th[2])**2)
        return 2*pr*th[0]*th[1]/th[2]**2*(2*np.cos(2*(x-xp)/th[2]) + th[1]*(-1 + np.cos(4*(x-xp)/th[2]))), False

    def d12covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        pr= np.exp(-2*th[1]*np.sin((x-xp)/th[2])**2)
        return -2*pr*th[0]*th[1]/th[2]**2*(2*np.cos(2*(x-xp)/th[2]) + th[1]*(-1 + np.cos(4*(x-xp)/th[2]))), False

    def d12d2covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 d/dxp of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        pr= np.exp(-2*th[1]*np.sin((x-xp)/th[2])**2)
        return -4*pr*th[0]*th[1]/th[2]**3*(2 - th[1]**2 + 6*th[1]*np.cos(2*(x-xp)/th[2])
                                           + th[1]**2*np.cos(4*(x-xp)/th[2]))*np.sin(2*(x-xp)/th[2]), False

    def d12d22covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 d^2/dxp^2 of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        pr= np.exp(-2*th[1]*np.sin((x-xp)/th[2])**2)
        return 2*pr*th[0]*th[1]/th[2]**4*((8-12*th[1]**2)*np.cos(2*(x-xp)/th[2])
                                          + th[1]*(-4 + 3*th[1]**2 - 4*(-7 + th[1]**2)*np.cos(4*(x-xp)/th[2])
                                                   + 12*th[1]*np.cos(6*(x-xp)/th[2])
                                                   + th[1]**2*np.cos(8*(x-xp)/th[2]))), False

####
class localperiodicGP(gaussianprocess):
    '''
    Gaussian process with a locally periodic covariance function.
    Returns the kernel function and Jacobian of the kernel function.

    Arguments
    --
    x: a 1-d array of abscissa
    xp: a 1-d array of alternative abscissa
    lth: the log of the hyperparameters
    '''
    noparams= 4
    description= 'locally periodic covariance function'

    def info(self):
        print('hparam[0] determines the amplitude of variation')
        print('hparam[1] determines the stiffness')
        print('hparam[2] determines the period')
        print('hparam[3] determines the local length scale')

    def covfn(self, x, xp, lth):
        '''
        Periodic covariance function.
        Returns the kernel function and Jacobian of the kernel function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        xp= np.asarray(xp)
        lpr= np.exp(-2*th[1]*np.sin((x-xp)/th[2])**2 - (x-xp)**2*th[3]/2)
        k= th[0]*lpr
        jk= np.empty((len(xp), self.noparams))
        jk[:,0]= lpr
        jk[:,1]= -2*lpr*th[0]*np.sin((x-xp)/th[2])**2
        jk[:,2]= 2/th[2]**2*lpr*th[0]*th[1]*(x-xp)*np.sin(2*(x-xp)/th[2])
        jk[:,3]= -lpr*th[0]*(x-xp)**2/2
        return k, jk

    def d1covfn(self, x, xp, lth):
        '''
        Returns d/dx of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        lpr= np.exp(-2*th[1]*np.sin((x-xp)/th[2])**2 - (x-xp)**2*th[3]/2)
        return lpr*th[0]*(th[3]*(-x + xp) - (2*th[1]*np.sin((2*(x - xp))/th[2]))/th[2]), False

    def d1d2covfn(self, x, xp, lth):
        '''
        Returns d/dx d/xp of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        lpr= np.exp(-2*th[1]*np.sin((x-xp)/th[2])**2 - (x-xp)**2*th[3]/2)
        return lpr*th[0]*((-2*th[1]**2 + th[2]**2*th[3] - th[2]**2*th[3]**2*x**2 + 2*th[2]**2*th[3]**2*x*xp
                    - th[2]**2*th[3]**2*xp**2
                    + 4*th[1]*np.cos((2*(x - xp))/th[2])
                    + 2*th[1]**2*np.cos((4*(x - xp))/th[2])
                    - 4*th[1]*th[2]*th[3]*x*np.sin((2*(x - xp))/th[2])
                    + 4*th[1]*th[2]*th[3]*xp*np.sin((2*(x - xp))/th[2]))/th[2]**2), False

    def d12covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        lpr= np.exp(-2*th[1]*np.sin((x-xp)/th[2])**2 - (x-xp)**2*th[3]/2)
        return lpr*th[0]*(2*th[1]**2 - th[2]**2*th[3] + th[2]**2*th[3]**2*x**2 - 2*th[2]**2*th[3]**2*x*xp
                           + th[2]**2*th[3]**2*xp**2 - 4*th[1]*np.cos((2*(x - xp))/th[2])
                           - 2*th[1]**2*np.cos((4*(x - xp))/th[2]) + 4*th[1]*th[2]*th[3]*x*np.sin((2*(x - xp))/th[2])
                           - 4*th[1]*th[2]*th[3]*xp*np.sin((2*(x - xp))/th[2]))/th[2]**2, False

    def d12d2covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 d/dxp of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        lpr= np.exp(-2*th[1]*np.sin((x-xp)/th[2])**2 - (x-xp)**2*th[3]/2)
        return -lpr*th[0]*((-6*th[1]**2*th[2]*th[3]*x + 3*th[2]**3*th[3]**2*x - th[2]**3*th[3]**3*x**3
                             + 6*th[1]**2*th[2]*th[3]*xp - 3*th[2]**3*th[3]**2*xp + 3*th[2]**3*th[3]**3*x**2*xp
                             - 3*th[2]**3*th[3]**3*x*xp**2 + th[2]**3*th[3]**3*xp**3
                             + 12*th[1]*th[2]*th[3]*(x - xp)*np.cos((2*(x - xp))/th[2])
                             + 6*th[1]**2*th[2]*th[3]*(x - xp)*np.cos((4*(x - xp))/th[2])
                             + 8*th[1]*np.sin((2*(x - xp))/th[2])
                             - 6*th[1]**3*np.sin((2*(x - xp))/th[2])
                             + 6*th[1]*th[2]**2*th[3]*np.sin((2*(x - xp))/th[2])
                             - 6*th[1]*th[2]**2*th[3]**2*x**2*np.sin((2*(x - xp))/th[2])
                             + 12*th[1]*th[2]**2*th[3]**2*x*xp*np.sin((2*(x - xp))/th[2])
                             - 6*th[1]*th[2]**2*th[3]**2*xp**2*np.sin((2*(x - xp))/th[2])
                             + 12*th[1]**2*np.sin((4*(x - xp))/th[2])
                             + 2*th[1]**3*np.sin((6*(x - xp))/th[2]))/th[2]**3), False

    def d12d22covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 d^2/dxp^2 of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        lpr= np.exp(-2*th[1]*np.sin((x-xp)/th[2])**2 - (x-xp)**2*th[3]/2)
        return lpr*th[0]*(-8*th[1]**2 + 6*th[1]**4 - 12*th[1]**2*th[2]**2*th[3] + 3*th[2]**4*th[3]**2
                           + 12*th[1]**2*th[2]**2*th[3]**2*x**2 - 6*th[2]**4*th[3]**3*x**2 + th[2]**4*th[3]**4*x**4
                           - 24*th[1]**2*th[2]**2*th[3]**2*x*xp + 12*th[2]**4*th[3]**3*x*xp
                           - 4*th[2]**4*th[3]**4*x**3*xp + 12*th[1]**2*th[2]**2*th[3]**2*xp**2
                           - 6*th[2]**4*th[3]**3*xp**2 + 6*th[2]**4*th[3]**4*x**2*xp**2 - 4*th[2]**4*th[3]**4*x*xp**3
                           + th[2]**4*th[3]**4*xp**4
                           - 8*th[1]*(-2 + 3*th[1]**2 + 3*th[2]**2*th[3]*(-1 + th[3]*(x - xp)**2))*np.cos((2*(x - xp))/th[2])
                           - 4*th[1]**2*(-14 + 2*th[1]**2 + 3*th[2]**2*th[3]*(-1 + th[3]*(x - xp)**2))*np.cos((4*(x - xp))/th[2])
                           + 24*th[1]**3*np.cos((6*(x - xp))/th[2]) + 2*th[1]**4*np.cos((8*(x - xp))/th[2])
                           - 32*th[1]*th[2]*th[3]*x*np.sin((2*(x - xp))/th[2])
                           + 24*th[1]**3*th[2]*th[3]*x*np.sin((2*(x - xp))/th[2])
                           - 24*th[1]*th[2]**3*th[3]**2*x*np.sin((2*(x - xp))/th[2])
                           + 8*th[1]*th[2]**3*th[3]**3*x**3*np.sin((2*(x - xp))/th[2])
                           + 32*th[1]*th[2]*th[3]*xp*np.sin((2*(x - xp))/th[2])
                           - 24*th[1]**3*th[2]*th[3]*xp*np.sin((2*(x - xp))/th[2])
                           + 24*th[1]*th[2]**3*th[3]**2*xp*np.sin((2*(x - xp))/th[2])
                           - 24*th[1]*th[2]**3*th[3]**3*x**2*xp*np.sin((2*(x - xp))/th[2])
                           + 24*th[1]*th[2]**3*th[3]**3*x*xp**2*np.sin((2*(x - xp))/th[2])
                           - 8*th[1]*th[2]**3*th[3]**3*xp**3*np.sin((2*(x - xp))/th[2])
                           - 48*th[1]**2*th[2]*th[3]*x*np.sin((4*(x - xp))/th[2])
                           + 48*th[1]**2*th[2]*th[3]*xp*np.sin((4*(x - xp))/th[2])
                           - 8*th[1]**3*th[2]*th[3]*x*np.sin((6*(x - xp))/th[2])
                           + 8*th[1]**3*th[2]*th[3]*xp*np.sin((6*(x - xp))/th[2]))/th[2]**4, False




####
class spectralmixtureGP(gaussianprocess):
    '''
    Gaussian process with a spectral mixture covariance function (Wilson & Adams, 2013).
    Returns the kernel function and Jacobian of the kernel function.

    UNFINISHED.

    Arguments
    --
    x: a 1-d array of abscissa
    xp: a 1-d array of alternative abscissa
    lth: the log of the hyperparameters
    '''
    description= 'spectral mixture covariance function'

    def info(self):
        print('hparams are grouped in threes for each Gaussian in the spectral density')
        print('hparam[0] determines the variance')
        print('hparam[1] determines the mean')
        print('hparam[2] determines the weight')


    def covfn(self, x, xp, lth):
        '''
        Periodic covariance function.
        Returns the kernel function and Jacobian of the kernel function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        xp= np.asarray(xp)
        noparams= len(self.b) - 1
        k= 0
        for i in np.arange(0, noparams, 3):
            k += (np.exp(-(x-xp)**2*th[i]/2)*np.sqrt(th[i])*th[i+2]
                  *np.cos(2*np.pi*(x-xp)*th[i+1]))
        jk= np.empty((len(xp), noparams))
        for i in np.arange(0, noparams, 3):
            smexp= np.exp(-(x-xp)**2*th[i]/2)
            jk[:,i]= -th[i+2]*(-1 + th[i]*(x-xp)**2)*np.cos(2*np.pi*th[i+1]*(x-xp))/2/np.sqrt(th[i])
            jk[:,i+1]= -2*np.pi*np.sqrt(th[i])*th[i+2]*(x-xp)*np.sin(2*np.pi*th[i+1]*(x-xp))
            jk[:,i+2]= np.sqrt(th[i])*np.cos(2*np.pi*th[i+1]*(x-xp))
        return k, jk


    def plotspectrum(self, npts= 100):
        th= np.exp(self.lth_opt)
        noparams= len(self.b) - 1
        # need to choose points around each mu
        s= np.append(np.linspace(th[1] - 5*np.sqrt(th[0]), th[1], npts),
                     np.linspace(th[1], th[1] + 5*np.sqrt(th[0]), npts))
        for i in np.arange(1, noparams, 3):
            s= np.append(s, np.linspace(th[i+1] - 5*np.sqrt(th[i]), th[i+1], npts))
            s= np.append(s, np.linspace(th[i+1], th[i+1] + 5*np.sqrt(th[i]), npts))
        s= np.sort(np.unique(s))
        sp= np.sum([th[i+2]*np.exp(-(s - th[i+1])**2/2/th[i])/np.sqrt(th[i]) for i in np.arange(0, noparams, 3)], axis= 0)
        plt.figure()
        plt.plot(s, sp)
        plt.yscale('log')
        plt.show()

####
class gaussianprocessException(Exception):
    __doc__= ''
    pass
