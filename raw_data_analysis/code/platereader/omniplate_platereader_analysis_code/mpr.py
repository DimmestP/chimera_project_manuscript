import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import genutils as gu
import pandas as pd

class processbatch:

    def __init__(self, datadirs, igwells= {}, igstrains= [], pkldir= False, logfile= False):
        '''
        Runs slpr on all .xlsx files in a list of directories to create pkl files.

        Arguments
        --
        datadirs: list of directories containing raw .xlsx data and _contents files
        strains: strains whose statistics should be processed (default is 'all')
        igwells: a dictionary listing the wells to be ignored for particualr data sets
        igstrains: a dictionary listing the strains to be ignored for particular data sets
        pkldir: directory to which the pkl files should be moved
        logfile: output from slpr is sent to a log file for each directory
        '''
        import os
        from platereader import slpr
        cwd= os.getcwd()
        # create pkldir if necessary
        if pkldir and not os.path.isdir(cwd + '/' + pkldir):
            os.makedirs(cwd + '/' + pkldir)
        # process data
        for wdir in datadirs:
            wdir += '/'
            for f in os.listdir(wdir):
                if f.split('.')[-1] == 'xlsx' and 'contents' not in f:
                    fname= f.split('.')[0]
                    print('processing', fname)
                    igw= igwells[fname] if fname in igwells else []
                    igs= igstrains[fname] if fname in igstrains else []
                    if logfile:
                        import sys
                        oldstdout= sys.stdout
                        sys.stdout= open(wdir + wdir.split('/')[-2] + '.batchlog', 'a')
                    self.p= slpr(f, wdir= wdir, ignorestrains= igs)
                    if igw: self.p.ignorewells(igw)
                    self.process()
                    self.p.save()
                    if logfile: sys.stdout= oldstdout
                    if pkldir:
                        os.rename(cwd + '/' + wdir + fname + '.pkl', cwd + '/' + pkldir + '/' + fname + '.pkl')
            if pkldir and logfile:
                os.rename(cwd + '/' + wdir + wdir.split('/')[-2] + '.batchlog',
                          cwd + '/' + pkldir + '/' + wdir.split('/')[-2] + '.batchlog')

    def process(self):
        '''
        Runs bespoke platereader analysis
        '''
        p= self.p
        p.correctmedia(figs= False)
        p.correctOD(figs= False)
        p.getstats(figs= False, keysmessage= False)



#############

class mpr:

    def __init__(self, direc, sortby= False, excludekeys= [], includekeys= [], quiet= False):
        '''
        Arguments
        --
        direc : directory containing a list of pickle files
        sortby : if False, create one dataframe; if 'strain' or 'condition', create dataframes for each strain or condition
        excludekeys : if specified as a string, such as 'GFP', all keys containing this keyword wil be excluded, such as 'GFPmn'
        includekeys : if specified, only include strains that have at least one key that contains these keywords
        quiet : if True, don't print any messages
        '''

        if direc[-1] != '/': direc += '/'
        excludekeys, includekeys= gu.makelist(excludekeys), gu.makelist(includekeys)
        data= {}
        tmax, difft, keys= [], [], []
        allconditions, allstrains, alldata= [], [], []
        self.notincluded= []
        # load processed files
        from os import listdir
        dcontents= listdir(direc)
        for fname in dcontents:
            if '.pkl' in fname:
                p= np.load(direc + fname)
                data[p.name]= p
                # store time information
                tmax.append(p.t[-1])
                difft.append(np.median(np.diff(p.t)))
                # store conditions and strains
                allconditions.append(p.getconditions(nomedia= True))
                allstrains.append([item for item in p.allstrains if 'null' not in item])
                alldata.append([item for item in p.alldata if ('null' not in item) and ('media' not in item)])
                if not quiet:
                    # print info
                    print(p.name + '\n---')
                    print('conditions:', '; '.join(p.getconditions(nomedia= True)))
                    print('max time: {:.2f}'.format(p.t[-1]))
                    print('strains:', '; '.join([item for item in p.allstrains if 'null' not in item]))
                    print('data types:', '; '.join(p.datatypes))
                    if hasattr(p, 'description'): print(p.description)
                    print('\n')
                # find keys that don't contain replicate data
                for c, s in p.getconditionsandstrains(nomedia= True, nonull= True):
                    lkeys= [key for key in p.d[c][s] if np.asarray(p.d[c][s][key]).ndim < 2]
                    keys.append(lkeys)
        keys= np.unique(np.concatenate(keys).flatten())
        if not quiet: print(len(data), 'pkl files loaded')

        # remove uninteresting keys
        ignorekeys= ['ignoredwells', 'originalplateloc', 'plateloc', 'time', 'OD gp']
        for key in keys:
            for ex in excludekeys:
                if ex in key: ignorekeys.append(key)
        nkeys= np.setdiff1d(keys, ignorekeys)

        # create attributes
        self.data= data
        self.allconditions= list(np.sort(np.unique(np.hstack(allconditions))))
        self.allstrains= list(np.sort(np.unique(np.hstack(allstrains))))
        # check for replicates
        ualld, repalld= np.unique(np.sort(np.hstack(alldata)), return_counts= True)
        self.alldata, self.replicates= list(ualld), list(repalld)

        # define common time vector
        self.t= np.arange(0, np.min(np.array(tmax)), np.median(np.array(difft)))

        ####

        # amalgamate data
        self.dd= {}
        from scipy.interpolate import interp1d
        for key in nkeys:
            # first make dictionary
            tdic= {}
            repdic= {}
            # process each plate
            for pname in data:
                p= data[pname]
                for c, s in p.getconditionsandstrains(nomedia= True, nonull= True):
                    if includekeys and np.any([key not in p.d[c][s].keys() for key in includekeys]):
                        if s + ' in ' + c not in self.notincluded:
                            self.notincluded.append(s + ' in ' + c)
                    else:
                        d= p.d[c][s][key]
                        timevar= True if np.size(d) > 1 else False
                        # label replicate conditions
                        if s + ' in ' + c in repdic:
                            repdic[s + ' in ' + c] += 1
                        else:
                            repdic[s + ' in ' + c]= 1
                        if self.replicates[self.alldata.index(s + ' in ' + c)] > 1:
                            rlabel = ' replicate ' + str(repdic[s + ' in ' + c])
                        else:
                            rlabel= ''
                        # add data to dictionary
                        if timevar:
                            # time-dependent variable
                            tdic[s + ' in ' + c + rlabel]= interp1d(p.t, d)(self.t)
                        else:
                            # summary statistic
                            tdic[s + ' in ' + c + rlabel]= [d]
            # then make dataframe from dictionary
            if timevar:
                if sortby:
                    # make dataframes for each strain or each condition
                    tds= {}
                    for sc in tdic:
                        s, c= sc.split(' in ')
                        if sortby == 'strain':
                            if s in tds:
                                tds[s][c]= tdic[sc]
                            else:
                                tds[s]= {c : tdic[sc]}
                        elif sortby == 'condition':
                            if c in tds:
                                tds[c][s]= tdic[sc]
                            else:
                                tds[c]= {s : tdic[sc]}
                        else:
                            print('Error:', sortby, 'option is not recognized')
                            return
                    df= {}
                    for keytds in tds:
                        dft= pd.DataFrame(tds[keytds])
                        dft.insert(0, 'time', self.t)
                        df[keytds]= dft
                else:
                    # make one dataframe
                    df= pd.DataFrame(tdic)
                    df.insert(0, 'time', self.t)
            else:
                # summary statistics
                df= self.processsummarystats(key, tdic)
            self.dd[key]= df


    def processsummarystats(self, key, tdic):
        '''
        Customizable routine to convert a dictionary (tdic) of a summary statistic (key) to a dataframe
        '''
        import pandas as pd
        rows= []
        for sc in tdic:
            s, c= sc.split(' in ')
            if 'replicate' in sc:
                r= int(c.split('replicate')[-1])
                c= c.split('replicate')[0]
            else:
                r= 1
            for da in tdic[sc]:
                rows.append({'strain' : s, 'condition' : c, 'replicate' : r, key: da})
        df= pd.DataFrame(rows)
        df= df[[key for key in rows[0]]]
        return df


    ####
    def plotdf(self, key, sorc= 'condition', ysize= 15):
        '''
        key: key for a scalar, such as 'max growth rate'
        sorc: key to sort by (default: 'condition')

        or

        key: key for a time-dependent variable
        sorc: unused
        '''
        if self.dd[key].shape[1] == 3:
            f, ax= plt.subplots(figsize= (6, ysize))
            sns.set_color_codes('pastel')
            sns.barplot(x= key, y= sorc, data= self.dd[key].sort_values(key, ascending= False), color= 'b')
            ax.set(ylabel= '', xlabel= key)
            plt.show()
        else:
            self.dd[key].plot(x= 'time').legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.show()
