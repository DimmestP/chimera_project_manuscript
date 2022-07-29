#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import gaussianprocess as gp
from fitderiv import fitderiv
import genutils as gu
import pandas as pd
import datetime
import logging
from io import StringIO
import inspect
try:
    import seaborn as sns
except ImportError:
    pass
plt.rcParams['figure.max_open_warning']= 0
idx= pd.IndexSlice


class platereader:
    '''
    for analyzing platereader data, correcting for autofluorescence, and determining growth rates.

    All data is stored used Panda's dataframes and plotted using Seaborn.

    Three dataframes are created. If p is an instance of the platereader class:
        p.r contains the raw data for each well in the plate;
        p.s contains the processed time-series using the data from all relevant wells;
        p.sc constains any summary statistics, such as 'max gr'

    A typical work flow is:

    import omniplate as om

    EITHER

    p= om.platereader('GALdata.xls', 'GALcontents.xlsx', wdir= 'data/')

    OR

    p= om.platereader()
    p.load('GALdata.xls', 'GALcontents.xlsx')

    THEN

    p.plot('OD', plate= True)
    p.correctOD()
    p.plot(y= 'OD')
    p.plot(y= 'OD', by= 'strain', conditionincludes= ['Gal', 'Glu'], strainexcludes= 'HXT7')
    p.correctmedia()
    p.correctauto()
    p.plot('c-GFPperod', by= 'condition')
    p.getstats('OD')
    p.savefigs()
    p.exportdf()

    See also http://swainlab.bio.ed.ac.uk/software/platereader/omniplate.html

    Wells can be removed using, for example,

    p.ignorewells(['B1', 'B2'])

    but the analysis routines should then be re-run.

    All wells can be reinstated with:

    p.ignorewells()

    General properties of the data and of previous processing are shown with:

    p.info()
    p.attributes()
    p.corrections()
    p.log()

    The data are stored in three Pandas dataframes:

    p.r: contains the data stored by well
    p.s: contains time-series of processed data, which are created from all the relevant wells
    p.sc: contains summary statistics, such as the maximum growth rate

    Bounds on the hyperparameters for all Gaussian processes can be specified. For example,

    p.correctauto(bd= {1: (-6, -4)})

    changes the bounds on the second hyperparameter for the Gaussian process used to fit data for the reference strain to be [10^-6, 10^-4]. This hyperparameter describes the flexibility of the fitted curve.

    As a citation for this software, please use:

    Lichten CA, White R, Clark IBN, Swain PS. Unmixing of fluorescence spectra to resolve quantitative time-series measurements of gene expression in plate readers. BMC Biotechnol 14 (2014) 11

    The technique for estimating growth rates is given in:

    Swain PS, Stevenson K, Leary A, Montano-Gutierrez LF, Clark IB, Vogel J, Pilizota T. Inferring time derivatives including cell growth rates using Gaussian processes. Nat Commun 7 (2016) 13766

    '''

    #####
    def __init__(self, dnames= False, anames= False,  wdir= '', platereadertype= 'Tecan', dsheetnumbers= False,
                 asheetnumbers= False, ODfname= 'ODcorrection_Glucose_Haploid.txt', info= True):
        '''
        Initiate, and potentially immediately load data for processing.

        Arguments
        --
        dnames: file name of data from the plate reader or a list of file names
        anames: file name or list of files names for the annotation (if unspecified assumed to be the root of the corresponding dname + 'Contents.xlsx)
        wdir: the working directory - name of directory where the data files are stored
        platereadertype: type of plate reader ('Tecan' or 'Sunrise' or 'old Tecan')
        dsheetnumbers: specifies the sheets of the data Excel files
        asheetnumbers: specifies the sheets of the annotation Excel files
        ODfname: file name for dilution data for corrected OD (default is 'ODcorrection_Glucose_Haploid.txt')
        info: if True (default), display information on the data after loading
        '''
        self.version= '0.15'
        print('\nomniplate version=', self.version)
        self.wdir= wdir

        # enable logging
        self.time= '{:%Y-%b-%d %H:%M:%S}'.format(datetime.datetime.now())
        self.logger= logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.logstream= StringIO()
        loghandler= logging.StreamHandler(self.logstream)
        self.logger.addHandler(loghandler)
        # startlog
        self.logger.info('om.platereader version ' +  f'{self.version}')
        self.logger.info(f'{self.time}\n')
        self.logmethod(self.logger)

        # general parameters
        self.gamma= 0.114   # ratio of 585 to 525 for eGFP
        self.nosamples= 100  # for estimating error through sampling
        self.consist= 2   # number of stds of corrected reference strain data that is considered measurement noise
        self.overflow= -999.99

        if True:
            # warning generated occasionally when sampling from the Gaussian process likely because of numerical errors
            import warnings
            warnings.simplefilter("ignore", RuntimeWarning)

        # dictionary recording extent of analysis
        self.progress= {'ODcorrected' : {}, 'ODcorrectedformedia' : {}, 'ignoredwells' : {},
                        'negativevalues' : {}, 'mediaGP' : {}, 'getstatsGP' : {},
                        'refGP' : {}, 'ODfname' : {}, 'gc' : {}}
        self.allexperiments= []
        self.alldata= {}
        self.allconditions= {}
        self.allstrains= {}
        self.datatypes= {}

        if dnames == False:
            # list all files in current directory
            self.ls()
        else:
            # immediately load data
            self.load(dnames, anames, platereadertype, dsheetnumbers, asheetnumbers, ODfname, info)




    #####
    def ls(self):
        '''
        List all files in the working directory.

        Arguments
        --
        None
        '''
        import os
        wdir= os.getcwd() + '/' + self.wdir
        print('Working directory is', wdir)
        print('Files available are:', '\n---')
        for f in os.listdir(wdir):
            if os.path.isfile(wdir + '/' + f): print(f)


    def changewdir(self, wdir):
        '''
        Change working directory

        Arguments
        ---
        wdir: new working directory specified from the current directory, e.g. 'data/'
        '''
        self.wdir= wdir
        self.ls()


    #####
    def load(self, dnames, anames= False, platereadertype= 'Tecan', dsheetnumbers= False, asheetnumbers= False,
             ODfname= 'ODcorrection_Glucose_Haploid.txt', info= True):
        '''
        Loads raw data files generated by the plate reader and the corresponding annotation files.

        Arguments
        --
        dnames: file name of data from the plate reader or a list of file names
        anames: file name or list of files names for the annotation (if unspecified assumed to be the root of the corresponding dname + 'Contents.xlsx)
        platereadertype: type of plate reader ('Tecan' or 'Sunrise' or 'old Tecan')
        dsheetnumbers: specifies the sheets of the data Excel files
        asheetnumbers: specifies the sheets of the annotation Excel files
        ODfname: file name for dilution data for corrected OD (default is 'ODcorrection_Glucose_Haploid.txt')
        info: if True (default), display information on the data after loading
        '''
        dnames= gu.makelist(dnames)
        if not anames:
            anames= [dname.split('.')[0] + 'Contents.xlsx' for dname in dnames]
        else:
            anames= gu.makelist(anames)
        if not dsheetnumbers: dsheetnumbers= [0 for dname in dnames]
        if not asheetnumbers: asheetnumbers= [0 for dname in dnames]

        for i, dname in enumerate(dnames):
            # dataframe for raw data
            # defines self.allexperiments, self.allconditions, self.allstrains, self.alldata
            rdf= self.importdatafromplate(platereadertype, dname, dsheetnumbers[i], anames[i], asheetnumbers[i])
            self.r= pd.merge(self.r, rdf, how= 'outer') if hasattr(self, 'r') else rdf
            # update progress dictionary
            experiment= dname.split('.')[0]
            self.progress['ODcorrected'][experiment]= False
            self.progress['ODcorrectedformedia'][experiment]= False
            self.progress['ignoredwells'][experiment]= []
            self.progress['negativevalues'][experiment]= False
            self.progress['mediaGP'][experiment]= {dn : {} for dn in self.datatypes[experiment]}
            self.progress['getstatsGP'][experiment]= {c : {} for c in self.allconditions[experiment]}
            self.progress['refGP'][experiment]= {c : {} for c in self.allconditions[experiment]}
            self.progress['gc'][experiment]= None

        # define ODfname in progress dictionary
        if isinstance(ODfname, str):
            self.progress['ODfname']= {exp: ODfname for exp in self.allexperiments}
        else:
            self.progress['ODfname']= {exp: ODfname[i] for i, exp in enumerate(self.allexperiments)}

        # dataframe for summary stats and corrections
        alldfs= []
        for exp in self.allexperiments:
            strs, cons= [], []
            for cs in self.alldata[exp]:
                strs.append(cs.split(' in ')[0])
                cons.append(cs.split(' in ')[1])
            corrdict= {'experiment' : exp, 'strain' : strs, 'condition' : cons, 'OD corrected' : False}
            corrdict.update({dtype + ' corrected for media' : False for dtype in self.datatypes[exp]})
            corrdict.update({dtype + ' corrected for autofluorescence' : False for dtype in self.datatypes[exp]
                             if dtype not in ['AutoFL', 'OD']})
            alldfs.append(pd.DataFrame(corrdict))
        self.sc= pd.concat(alldfs)

        # dataframe of original data
        self.origr= self.r.copy()
        # dataframe for well content
        self.wellsdf= self.makewellsdf()
        # dataframe for summary data
        self.s= self.make_s()

        # display info on experiment, conditions and strains
        if info: self.info()



    #####
    def importdatafromplate(self, platereadertype, dname, dsheetnumber, aname, asheetnumber):
        '''
        Internal function: Creates dataframe from input files created by the plate reader.
        '''
        ###
        # import annotation
        ###
        blank= 'Blank'
        exp= dname.split('.')[0]
        try:
            anno= pd.read_excel(self.wdir + aname, index_col= 0, sheet_name= asheetnumber)
            self.alldata[exp]= []
            rcontents= {}
            # run through all wells
            for x in np.arange(1,13):
                for y in 'ABCDEFGH':
                    well= y + str(x)
                    if isinstance(anno[x][y], str):
                        s, c= anno[x][y].split(' in ')
                        rcontents[well]= [c.strip(), s.strip()]
                        self.alldata[exp].append(rcontents[well][1] + ' in ' + rcontents[well][0])
                    else:
                        rcontents[well]= [blank, blank]
            # summaries of data
            self.alldata[exp]= list(np.unique(self.alldata[exp]))
            self.allconditions[exp]= list(np.unique([rcontents[well][0] for well in rcontents
                                                if rcontents[well][0] != blank]))
            self.allstrains[exp]= list(np.unique([rcontents[well][1] for well in rcontents
                                            if rcontents[well][0] != blank]))
        except FileNotFoundError:
            raise(SystemExit("\nError: Can't find " + self.wdir + aname))

        ###
        # import data created by plate reader
        ###
        try:
            print('loading', dname)
            dfd= pd.read_excel(self.wdir + dname, sheet_name= dsheetnumber)
            experiment= dname.split('.')[0]
            self.allexperiments.append(experiment)
        except FileNotFoundError:
            raise(SystemExit("\nError: Can't find " + self.wdir + dname))

        ###
        # create dataframes for raw and processed data
        ###
        # find types of data and times of measurement (in hrs)
        if platereadertype == 'Tecan':
            datatypes= dfd[dfd.columns[0]].iloc[np.nonzero(dfd[dfd.columns[0]].str.startswith('Cycle Nr', na= False).values)[0]-1].values
            # if only OD data measured
            if not isinstance(datatypes[0], str): datatypes= ['OD']
            t= dfd.loc[dfd[dfd.columns[0]].str.startswith('Time [s]', na= False),
                       dfd.columns[1]:].dropna(axis= 'columns').mean().values.astype('float')/3600.0
        elif platereadertype == 'old Tecan':
            datatypes= [dfd[dfd.columns[0]].iloc[np.nonzero(dfd[dfd.columns[0]].str.startswith('Label', na= False).values)[0]].values[0].split(': ')[1]]
            t= dfd.loc[dfd[dfd.columns[0]].str.startswith('Time [s]', na= False),
                       dfd.columns[1]:].dropna(axis= 'columns').mean().values.astype('float')/3600.0
        elif platereadertype == 'Sunrise':
            datatypes= np.array(['OD'], dtype= object)
            t= gu.rmnans([float(str(ts).split('s')[0]) for ts in dfd.values[0]])/3600.0
        else:
            raise(SystemExit('\nError: ' + platereadertype + ' not recognized'))
        self.datatypes[exp]= list(datatypes)


        ###
        # extract data
        ###
        if platereadertype == 'Tecan' or platereadertype == 'old Tecan':
            df= dfd.replace('OVER', self.overflow)
            cols= df.columns
            ## metadata (not used)
            # plate reader machine
            serialnumber= df[cols[cols.str.startswith('Tecan i-control')]].iloc[0,0].split()[2]
            # date of experiment
            expdate= df[cols[1]].where(df[cols[0]] == 'Date:').dropna().iloc[0]
            # gains
            gains= df.loc[df[cols[0]].str.startswith('Gain', na= False),:].iloc[:,4].values
            ## data
            # add to dataframe
            df.index= df[cols[0]]
            rdict= []
            for x in np.arange(1,13):
                for y in 'ABCDEFGH':
                    well= y + str(x)
                    if well in df.index:
                        data= df.loc[well, cols[1]:].values.astype(np.float)
                        if data.ndim == 1: data= data[None,:]
                        if rcontents[well][0] != blank and rcontents[well][1] != blank:
                            for j in range(len(t)):
                                cons= {'experiment' : experiment, 'condition' : rcontents[well][0],
                                       'strain' : rcontents[well][1], 'time' : t[j], 'well' : well}
                                dats= {datatype : data[i,j] for i, datatype in enumerate(datatypes)}
                                cons.update(dats)
                                rdict.append(cons)

        elif platereadertype == 'Sunrise':
            ## metadata (not used)
            fcol= dfd[dfd.columns[0]]
            # plate reader machine
            serialnumber= str(fcol[fcol.str.startswith('Instrument serial', na= False)]).split(':')[1].split('\n')[0].strip()
            # date of experiment
            expdate= str(fcol[fcol.str.startswith('Date', na= False)]).split(':')[1].split('/')[0].strip()
            ## data
            # add to dataframe
            rdict= []
            for index, row in dfd.iterrows():
                if isinstance(row[-1], str) and row[-1][0] in 'ABCDEFGH':
                    well= row[-1]
                    data= row.values[:-1].astype(np.float)
                    for j in range(len(t)):
                        cons=  {'experiment' : experiment, 'condition' : rcontents[well][0],
                               'strain' : rcontents[well][1], 'time' : t[j], 'well' : well}
                        dats= {'OD': data[j]}
                        cons.update(dats)
                        rdict.append(cons)

        return pd.DataFrame(rdict)


    #####
    # Routines to display information on data and state of data processing
    #####
    def info(self):
        '''
        Displays conditions, strains, and datatypes.

        Arguments
        --
        None
        '''
        for exp in self.allexperiments:
            print('\nExperiment:', exp, '\n---')
            print('Conditions:')
            for c in self.allconditions[exp]: print('\t', c)
            print('Strains:')
            for s in self.allstrains[exp]: print('\t', s)
            print('Data types:')
            for d in self.datatypes[exp]: print('\t', d)
            if self.progress['ignoredwells']:
                print('Ignored wells:')
                if self.progress['ignoredwells'][exp]:
                    for d in self.progress['ignoredwells'][exp]: print('\t', d)
                else:
                    print('\t', 'None')


    def attributes(self):
        '''
        Displays the names of the attributes in the instance of platereader.

        Arguments
        --
        None
        '''
        ignore= ['d', 'consist', 't', 'nosamples', 'gamma', 'ODfname', 'overflow', 'nooutchannels', 'nodata', '__doc__']
        for a in self.__dict__:
            if 'corrected' not in a and 'processed' not in a and a not in ignore: print(a)



    def corrections(self, experiments= 'all', conditions= 'all', strains= 'all',
                    experimentincludes= False, experimentexcludes= False,
                    conditionincludes= False, conditionexcludes= False,
                    strainincludes= False, strainexcludes= False):
        '''
        Displays the status of corrections made for the various strains in various conditions and experiments.

        Arguments
        --
        experiments: list of experiments to include (default is 'all')
        conditions: list of conditions to include (default is 'all')
        strains: list of strains to include (default is 'all')
        experimentincludes: selects only experiments with experimentincludes in their name
        experimentexcludes: ignores experiments with experimentexcludes in their name
        conditionincludes: selects only conditions with conditionincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        strainincludes: selects only strains with strainincludes in their name
        strainexcludes: ignore strains with strainexcludes in their name
        '''
        exps, cons, strs= self.getall(experiments, experimentincludes, experimentexcludes,
                                      conditions, conditionincludes, conditionexcludes,
                                      strains, strainincludes, strainexcludes)
        df= self.sc.query('experiment == @exps and condition == @cons and strain == @strs')
        # only show corrections and not stats
        df= df[['experiment', 'strain', 'condition'] + [col for col in df.columns if 'correct' in col]]
        df= df.T
        return df



    #####
    # Routines to ignore wells that are deemed poor quality
    #####
    def makewellsdf(self):
        '''
        Internal function: makes a dataframe that stores the contents of the wells
        '''
        df= self.r[['experiment', 'condition', 'strain', 'well']].drop_duplicates()
        df= df.reset_index(drop= True)
        return df


    def contentsofwells(self, wlist):
        '''
        Displays contents of wells

        Arguments:
        --
        wlist: a well or a list of wells
        '''
        wlist= gu.makelist(wlist)
        for w in wlist:
            print('\n' + w + '\n--')
            print(self.wellsdf.query('well == @w').drop(['well'], axis= 1).to_string(index= False))


    def wellsfor(self, experiments= 'all', conditions= 'all', strains= 'all',
                 experimentincludes= False, experimentexcludes= False,
                 conditionincludes= False, conditionexcludes= False,
                 strainincludes= False, strainexcludes= False):
        '''
        Displays wells for specified conditions and strains

        Arguments
        --
        experiments: list of experiments to include (default is 'all')
        conditions: list of conditions to include (default is 'all')
        strains: list of strains to include (default is 'all')
        experimentincludes: selects only experiments with experimentincludes in their name
        experimentexcludes: ignores experiments with experimentexcludes in their name
        conditionincludes: selects only conditions with conditionincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        strainincludes: selects only strains with strainincludes in their name
        strainexcludes: ignore strains with strainexcludes in their name
        '''
        exps, cons, strs= self.getall(experiments, experimentincludes, experimentexcludes,
                                      conditions, conditionincludes, conditionexcludes,
                                      strains, strainincludes, strainexcludes, nonull= False)
        df= self.wellsdf.query('experiment == @exps and condition == @cons and strain == @strs')
        print(df.to_string(index= False))


    def ignorewells(self, exclude= [], experiments= 'all', experimentincludes= False, experimentexcludes= False,
                    clearall= False):
        '''
        Allows wells to be ignored in any future processing. If called several times, the default behaviour is for any previously ignored wells not to be re-instated.

        Arguments
        --
        exclude: list of labels (platereader locations) of wells to be excluded
        experiments: list of experiments to include (default is 'all')
        experimentincludes: selects only experiments with experimentincludes in their name
        experimentexcludes: ignores experiments with experimentexcludes in their name
        clearall: if True, all previously ignored wells are re-instated
        '''
        self.logmethod(self.logger)
        if clearall:
            # forget any previously ignoredwells
            self.r= self.origr.copy()
            self.progress['ignoredwells']= {exp: [] for exp in self.allexperiments}
            self.update_s()
            print('Warning: all corrections and analysis to raw data have been lost. No wells have been ignored.')
        else:
            if gu.islistempty(exclude):
                return
            else:
                # exclude should be a list of lists
                if isinstance(exclude, str):
                    exclude= [gu.makelist(exclude)]
                elif isinstance(exclude[0], str):
                    exclude= [exclude]
                # check consistency
                if len(self.allexperiments) == 1:
                    exps= self.allexperiments
                else:
                    exps= self.getexps(experiments, experimentincludes, experimentexcludes)
                if len(exclude) != len(exps) and not clearall:
                    print('Error: A list of wells to exclude must be given for a particular experiment or a list of experiments.')
                    return
                else:
                    # drop wells
                    for ex, exp in zip(exclude, exps):
                        # wells cannot be ignored twice
                        wex= list(set(ex) - set(self.progress['ignoredwells'][exp]))
                        # delete wells
                        df= self.r
                        df= df.set_index(['experiment', 'well'])
                        for well in wex:
                            df= df.drop((exp, well))
                        df.reset_index(level= [0, 1], inplace= True)
                        self.r= df
                        # store ignoredwells
                        self.progress['ignoredwells'][exp] += ex
                        # remove any duplicates
                        self.progress['ignoredwells'][exp]= list(set(self.progress['ignoredwells'][exp]))
                # remake summary data
                self.update_s()




    #####
    # Routines to make and update the dataframe of summary data (over wells)
    #####
    def make_s(self):
        '''
        Internal function: Calculates means and variances of all datatypes from raw data
        '''
        # find means
        df1= self.r.groupby(['experiment', 'condition', 'strain', 'time']).mean().reset_index()
        for exp in self.allexperiments:
            for dtype in self.datatypes[exp]:
                df1= df1.rename(columns= {dtype : dtype + ' mean'})
        # find variances
        df2= self.r.groupby(['experiment', 'condition', 'strain', 'time']).var().reset_index()
        for exp in self.allexperiments:
            for dtype in self.datatypes[exp]:
                df2= df2.rename(columns= {dtype : dtype + ' var'})
        return pd.merge(df1, df2)


    def update_s(self):
        '''
        Internal function: Updates means and variances of all datatypes from raw data
        '''
        self.s.update(self.make_s())



    #####
    # Routines for plotting
    #####
    def plot(self, x= 'time', y= 'OD', by= 'well', hue= None, style= None, kind= 'line', size= None,
             col= None, row= None, height= 5, aspect= 1, ymin= False,
             onefig= False, title= False, plate= False, nonull= False, messages= False,
             experiments= 'all', conditions= 'all', strains= 'all',
             experimentincludes= False, experimentexcludes= False,
             conditionincludes= False, conditionexcludes= False,
             strainincludes= False, strainexcludes= False, **kwargs):
        '''
        Plots from the underlying dataframes (chosen automatically) using Seaborn's relplot

        https://seaborn.pydata.org/generated/seaborn.relplot.html

        Arguments
        --
        x: variable for x-axis (default: 'time')
        y: variable for y-axis (default: 'OD')
        by: show data for each instance of this variable (either 'well' - default, 'experiment', 'condition', or 'strain')
        hue: (Seaborn) variable that determines the color of each line
        style: (Seaborn) variable that determines the style of each line
        kind: (Seaborn) either 'line' (default) or 'scatter'
        size: (Seaborn) variable that determines the size of each marker
        col: (Seaborn) variable that determines the columns in a multifigure plot
        row: (Seaborn) variables that determines the rows in a multifigure plot
        height: (Seaborn) height of the subplots in a multifigure plot
        aspect: (Seaborn) aspect ratio of the subplots in a multifigure plot
        ymin: if specified, sets the minimum y-value
        onefig: whether multiple or a single figure should be produced (default: False)
        title: if specified, the title of the plot (overwrites any default titles)
        plate: if True, data for each well are plotted in one figure
        nonull: if True, do not plot 'null' strains
        messsages: if True, print warnings for any data requested but not found
        experiments: list of experiments to include (default is 'all')
        conditions: list of conditions to include (default is 'all')
        strains: list of strains to include (default is 'all')
        experimentincludes: selects only experiments with experimentincludes in their name
        experimentexcludes: ignores experiments with experimentexcludes in their name
        conditionincludes: selects only conditions with conditionincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        strainincludes: selects only strains with strainincludes in their name
        strainexcludes: ignore strains with strainexcludes in their name
        '''
        if plate:
            if x == 'time':
                self.plotplate(y)
            else:
                self.plotplate(x)
        else:
            # get experiments, conditions and strains
            exps, cons, strs= self.getall(experiments, experimentincludes, experimentexcludes,
                                          conditions, conditionincludes, conditionexcludes,
                                          strains, strainincludes, strainexcludes, nonull)

            # choose the correct dataframe
            if hasattr(self, 'r') and x in self.r.columns and y in self.r.columns:
                # raw data (with wells)
                basedf= self.r
            elif x in self.s.columns and y in self.s.columns:
                # processed data (no wells)
                basedf= self.s
                if by == 'well': by= 'condition'
                nonull= True
            elif x in self.sc.columns and y in self.sc.columns:
                # summary stats
                basedf= self.sc
                df= basedf.query('experiment == @exps and condition == @cons and strain == @strs')
                if df.empty:
                    print('No data found')
                else:
                    sfig= sns.relplot(x, y, data= df, hue= hue, kind= 'scatter', style= style, size= size,
                                      col= col, row= row, aspect= aspect, height= height, **kwargs)
                    if row == None and col == None:
                        # add error bars
                        # find coordinates of points in relplot
                        xc, yc = [], []
                        for point_pair in sfig.ax.collections:
                            for xp, yp in point_pair.get_offsets():
                                xc.append(xp)
                                yc.append(yp)
                        # add error bars
                        xerr= np.sqrt(df[x + ' var'].values) if x + ' var' in df.columns else None
                        yerr= np.sqrt(df[y + ' var'].values) if y + ' var' in df.columns else None
                        if np.any(xerr): xerr= gu.rmnans(xerr)
                        if np.any(yerr): yerr= gu.rmnans(yerr)
                        sfig.ax.errorbar(xc, yc, xerr= xerr, yerr= yerr, fmt= ' ', ecolor= 'dimgray', alpha= 0.5)
                    return

            else:
                print('The variables x=', x, 'and y=', y, 'cannot be plotted against each other because they are not in the same dataframe.')
                return

            # defaults
            if col or row: onefig= True

            if onefig:
                # single plot
                df= basedf.query('experiment == @exps and condition == @cons and strain == @strs')
                if df.empty:
                    print('No data found')
                else:
                    df= self.augmentdf(df, y)
                    if by == 'well':
                        ell= df['well'].unique().size
                        if not hue: hue= 'well'
                        if not style: style == 'condition'
                        if not size: size= 'experiment'
                    elif by == 'strain':
                        ell= len(cons)
                        if not hue: hue= 'strain'
                        if not style: style= 'condition'
                        if not size: size= 'experiment'
                    elif by == 'condition':
                        ell= len(strs)
                        if not hue: hue= 'condition'
                        if not style: style= 'strain'
                        if not size: size= 'experiment'
                    elif by == 'experiment':
                        ell= len(exps)
                        if not hue: hue= 'experiment'
                        if not style: style= 'condition'
                        if not size: size= 'strain'
                    else:
                        print(by, 'not recognized')
                        return
                    if ell > 6:
                        # seaborn allows only 6 different line styles
                        print('Too many ' + by + 's for seaborn to plot with different line style')
                        return
                    else:
                        lhue= self.checkforseaborn(df, hue)
                        sfig= sns.relplot(x, y, data= df, hue= lhue, kind= kind, style= style, ci= 'sd', size= size,
                                          col= col, row= row, aspect= aspect, height= height, **kwargs)
                        if title: sfig.fig.suptitle(title)
                        if ymin is not False: plt.ylim(ymin, None)
            else:
                # multiple plots
                if by == 'well':
                    for e in exps:
                        for c in cons:
                            for s in strs:
                                df= basedf.query('experiment == @e and condition == @c and strain == @s')
                                if df.empty:
                                    if messages: print(e + '::', 'No data found for', s, 'in', c)
                                else:
                                    sfig= sns.relplot(x, y, data= df, hue= 'well', kind= kind, style= style, size= size, **kwargs)
                                    if title:
                                        sfig.fig.suptitle(title)
                                    else:
                                        sfig.fig.suptitle(e + ':: ' + s + ' in ' + c)
                                    if ymin is not False: plt.ylim(ymin, None)
                elif by == 'strain':
                    if not hue: hue= 'condition'
                    for e in exps:
                        for s in strs:
                            df= basedf.query('experiment == @e and condition == @cons and strain == @s')
                            if df.empty:
                                if messages: print(e + '::', 'No data found for', s)
                            else:
                                df= self.augmentdf(df, y)
                                lhue= self.checkforseaborn(df, hue)
                                sfig= sns.relplot(x, y, data= df, hue= lhue, kind= kind, ci= 'sd', style= style, size= size, **kwargs)
                                if title:
                                    sfig.fig.suptitle(title)
                                else:
                                    sfig.fig.suptitle(e + ':: ' + s)
                                if ymin is not False: plt.ylim(ymin, None)
                elif by == 'condition':
                    if not hue: hue= 'strain'
                    for e in exps:
                        for c in cons:
                            df= basedf.query('experiment == @e and condition == @c and strain == @strs')
                            if df.empty:
                                if messages: print(e + '::', 'No data found for', c)
                            else:
                                df= self.augmentdf(df, y)
                                lhue= self.checkforseaborn(df, hue)
                                sfig= sns.relplot(x, y, data= df, hue= lhue, kind= kind, ci= 'sd', style= style, size= size, **kwargs)
                                if title:
                                    sfig.fig.suptitle(title)
                                else:
                                    sfig.fig.suptitle(e + ':: ' + c)
                                if ymin is not False: plt.ylim(ymin, None)
                elif by == 'experiment':
                    for c in cons:
                        for s in strs:
                            df= basedf.query('condition == @c and strain == @s')
                            if df.empty:
                                if messages: print('No data found for', s, 'in', c)
                            else:
                                df= self.augmentdf(df, y)
                                sfig= sns.relplot(x, y, data= df, hue= 'experiment', kind= kind, style= style, size= size, **kwargs)
                                if title:
                                    sfig.fig.suptitle(title)
                                else:
                                    sfig.fig.suptitle(s + ' in ' + c)
                                if ymin is not False: plt.ylim(ymin, None)
                else:
                    print(by, 'not recognized')
                    return
        plt.show()



    def checkforseaborn(self, df, hue):
        # prevent seaborn crashing when called with hue of length one
        if df[hue].unique().size == 1:
            hue= None
        return hue



    def plotplate(self, dtype):
        '''
        Plots the data for each well following the layout of a 96-well plate.

        Arguments
        --
        dtype: data type to be plotted
        '''
        plt.figure()
        for pl in self.r['well'].unique():
            rowl= 'ABCDEFGH'.index(pl[0])
            coll= int(pl[1:])
            sindex= coll + 12*rowl
            plt.subplot(8, 12, sindex)
            # plot data
            wd= self.r.query('well == @pl')
            plt.plot(wd['time'].values, wd[dtype].values, '-')
            plt.tick_params(labelbottom= False, labelleft= False)
            # label well locations
            for j in range(12):
                if sindex == j+1: plt.title(j+1)
            for j, k in enumerate(np.arange(1, 96, 12)):
                if sindex == k: plt.ylabel('ABCDEFGH'[j] + ' ', rotation= 0)
        plt.suptitle(dtype)
        plt.show(block= False)




    #####
    def savefigs(self, ftype= 'pdf', onefile= True):
        '''
        Saves all current figures, either each to a separate file (default) or all to one file.

        Arguments
        --
        ftype: type of file ('pdf' is the default)
        onefile: if True (default), all figures are saved to one PDF file
        '''
        if onefile:
            gu.figs2pdf(self.wdir + ''.join(self.allexperiments) + '.pdf')
        else:
            for i in plt.get_fignums():
                plt.figure(i)
                savename= str(plt.getp(plt.gcf(), 'axes')[0].title).split("'")[1]
                savename= savename.replace(' ', '_')
                if savename == '': savename= 'Whole_plate_Figure_' + str(i)
                print('Saving', savename)
                plt.savefig(self.wdir + savename + '.' + ftype)


    #####
    def close(self):
        '''
        Close all figures.

        Arguments
        --
        None
        '''
        plt.close('all')



    #####
    # Internal functions
    #####
    def getsubset(self, type, set= 'all', includes= False, excludes= False,
                  nonull= False, nomedia= False):
        '''
        Returns a subset of all possibilities.

        Arguments
        --
        type: 'c' (conditions) or 's' (strains)
        set: list of items to include (default is 'all')
        includes: selects only items with includes in their name
        excludes: ignores items with excludes in their name
        nonull: if True, ignore null strain
        nomedia: if True, ignores 'media' condition
        '''
        if set == 'all' or includes or excludes:
            if type == 'c':
                sset= list(np.unique([con for e in self.allconditions for con in self.allconditions[e]]))
                if nomedia and 'media' in sset: sset.pop(sset.index('media'))
            elif type == 's':
                sset= list(np.unique([str for e in self.allstrains for str in self.allstrains[e]]))
                if nonull and 'null' in sset: sset.pop(sset.index('null'))
            else:
                sset= self.allexperiments
            # find those items containing keywords given in 'includes'
            if includes:
                includes= gu.makelist(includes)
                newset= []
                for s in sset:
                    gotone= 0
                    for item in includes:
                        if item in s: gotone += 1
                    if gotone == len(includes): newset.append(s)
                sset= newset
            # remove any items containing keywords given in 'excludes'
            if excludes:
                excludes= gu.makelist(excludes)
                exs= []
                for s in sset:
                    for item in excludes:
                        if item in s:
                            exs.append(s)
                            break
                for ex in exs:
                    sset.pop(sset.index(ex))
        else:
            sset= gu.makelist(set)
        if sset:
            return sorted(sset)
        else:
            if includes:
                raise(SystemExit('Nothing found for ' + ' and '.join(includes)))
            else:
                raise(SystemExit('None found'))



    def getexps(self, experiments, experimentincludes, experimentexcludes):
        '''
        Internal function: returns list of experiments
        '''
        if experimentincludes or experimentexcludes:
            exps= self.getsubset('e', includes= experimentincludes, excludes= experimentexcludes)
        elif experiments == 'all':
            exps= self.allexperiments
        else:
            exps= gu.makelist(experiments)
        return exps


    def getcons(self, conditions, conditionincludes, conditionexcludes, nomedia):
        '''
        Internal function: returns list of conditions
        '''
        if conditionincludes or conditionexcludes:
            cons= self.getsubset('c', includes= conditionincludes, excludes= conditionexcludes, nomedia= nomedia)
        elif conditions == 'all':
            cons= list(np.unique([con for e in self.allconditions for con in self.allconditions[e]]))
            if nomedia and 'media' in cons: cons.pop(cons.index('media'))
        else:
            cons= gu.makelist(conditions)
        return cons


    def getstrs(self, strains, strainincludes, strainexcludes, nonull):
        '''
        Internal function: returns list of strains
        '''
        if strainincludes or strainexcludes:
            strs= self.getsubset('s', includes= strainincludes, excludes= strainexcludes, nonull= nonull)
        elif strains == 'all':
            strs= list(np.unique([str for e in self.allstrains for str in self.allstrains[e]]))
            if nonull and 'null' in strs: strs.pop(strs.index('null'))
        else:
            strs= gu.makelist(strains)
        if nonull and 'null' in strs: strs.pop(strs.index('null'))
        return strs


    def getall(self, experiments, experimentincludes, experimentexcludes,
               conditions, conditionincludes, conditionexcludes,
               strains, strainincludes, strainexcludes,
               nonull= True, nomedia= True):
        '''
        Internal function: returns experiments, conditions, and strains
        '''
        exps= self.getexps(experiments, experimentincludes, experimentexcludes)
        cons= self.getcons(conditions, conditionincludes, conditionexcludes, nomedia)
        strs= self.getstrs(strains, strainincludes, strainexcludes, nonull)
        return exps, cons, strs


    def extractwells(self, experiment, condition, strain, datatypes):
        '''
        Internal function: extracts a list of matrices for each dtype in datatypes for the given
        experiment, condition, and strain with each column in each matrix having the data for one well
        '''
        datatypes= gu.makelist(datatypes)
        df= self.r.query('experiment == @experiment and condition == @condition and strain == @strain')
        matrices= []
        for dtype in datatypes:
            df2= df[[dtype, 'well']]
            df2well= df2.groupby('well')[dtype].apply(list)
            matrices.append(np.transpose([df2well[w] for w in df2well.index]))
        if len(datatypes) == 1:
            # return array
            return matrices[0]
        else:
            # return list of arrays
            return matrices


    def augmentdf(self, df, datatype):
        '''
        Internal function: artifically augments dataframe using 'var' (if present in the dataframe) to allow seaborn to generate errors
        '''
        if datatype + ' var' in df:
            dvar= datatype + ' var'
        elif 'mean' in datatype and datatype.split(' mean')[0] + ' var' in df:
            dvar= datatype.split(' mean')[0] + ' var'
        else:
            dvar= False
        if dvar:
            df.insert(0, 'augtype', 'mean')
            mn= df[datatype].values
            std= np.sqrt(df[dvar].values)
            # add std
            dfp= df.copy()
            dfp[datatype]= mn + std
            dfp['augtype']= '+std'
            # minus std
            dfm= df.copy()
            dfm[datatype]= mn - std
            dfm['augtype']= '-std'
            # append
            df= df.append([dfp, dfm])
        return df



    #####
    # OD correction
    #####
    def correctOD(self, figs= True, correctmedia= True, commonmedia= False, mediafigs= False,
                  experiments= 'all', experimentincludes= False, experimentexcludes= False,
                  correctmediamean= False, correctmediaiskip= False):
        '''
        Corrects OD data for a non-linear relationship between OD and cell number.
        Requires a dilution data set in the file ODfname. The default dilution curve is 'ODcorrection.txt' (measured by C Lichten in 2010).

        Arguments
        --
        figs: if True, a plot of the fit to the dilution data is produced
        correctmedia: if True (default), correct OD measurements by the OD of the media
        commonmedia: (for correctmedia) a condition containing a 'null' strain that should be used to correct other conditions
        mediafigs: if True, display figures generated by correctmedia
        experiments: list of experiments to include (default is 'all')
        experimentincludes: selects only experiments with experimentincludes in their name
        experimentexcludes: ignores experiments with experimentexcludes in their name
        correctmediamean: if True, the mean over time of the media values is used for the correction (rather than fitting with a Gaussian process)
        correctmediaiskip: use only every iskip'th data point to increase speed (must be an integer and only used for correctmedia)
        '''
        self.logmethod(self.logger)
        exps= self.getexps(experiments, experimentincludes, experimentexcludes)
        for exp in exps:
            if not self.progress['ODcorrected'][exp]:
                # correct for media
                if correctmedia:
                    self.correctmedia(datatypes= 'OD', mean= correctmediamean, figs= mediafigs, commonmedia= commonmedia,
                                      experiments= experiments, experimentincludes= experimentincludes, experimentexcludes= experimentexcludes,
                                      log= False, iskip= correctmediaiskip)
                    self.progress['ODcorrectedformedia'][exp]= True
                # fit dilution data
                if not self.progress['gc'][exp]:
                    self.findODcorrection(self.wdir + self.progress['ODfname'][exp], exp, figs)
                # correct all wells
                gc= self.progress['gc'][exp]
                gc.batchpredict(self.r.query('experiment == @exp')['OD'].values)
                self.r.loc[self.r['experiment'] == exp, 'OD']= gc.f
                self.progress['ODcorrected'][exp]= True
                # remake summary data
                self.update_s()
                # flag corrections in summary stats dataframe
                self.sc.loc[self.sc['experiment'] == exp, 'OD corrected']= True
            else:
                print('OD is already corrected')


    #####
    def findODcorrection(self, ODfname, exp, figs):
        '''
        Internal function: Uses a Gaussian process to fit serial dilution data to correct for non-linearities in the relationship between OD and cell density. The data are expected in a file ODfname.
        '''
        print('Fitting dilution data for OD correction for non-linearities')
        try:
            od, dilfac= np.loadtxt(ODfname, unpack= True)
            print('Using', ODfname)
            # process data
            dilfac= dilfac[np.argsort(od)]
            od= np.sort(od)
            # run Gaussian process
            gc= gp.sqexplinGP({0: (-4, 2), 1: (-3, 1), 2: (-6, 1), 3: (-6, 1)}, od, dilfac)
            gc.findhyperparameters(noruns= 5, exitearly= True, quiet= True)
            gc.predict(od)
            if figs:
                plt.figure()
                gc.sketch('.')
                plt.xlim([0, 1.05*np.max(od)])
                plt.ylim([0, 1.05*np.max(od)])
                plt.gca().set_aspect('equal', adjustable='box')
                plt.draw()
                plt.grid(True)
                plt.xlabel('OD')
                plt.ylabel('relative cell density')
                plt.title('Fitting ' + ODfname)
                plt.show(block= False)
            self.progress['gc'][exp]= gc
            # copy gc to experiments with the same ODfname
            for e in self.allexperiments:
                if self.progress['ODfname'][e] == self.progress['ODfname'][exp]:
                    self.progress['gc'][e]= gc
        except FileNotFoundError:
            raise(SystemExit("\nError: Can't find " + ODfname))




    #####
    # Media correction
    #####
    def correctmedia(self, datatypes= 'all', mean= False, commonmedia= False,
                     experiments= 'all', experimentincludes= False, experimentexcludes= False,
                     conditions= 'all', conditionincludes= False, conditionexcludes= False,
                     figs= True, noruns= 3, exitearly= False, bd= False, results= False,
                     log= True, iskip= False):
        '''
        Corrects OD or fluorescence for the OD or fluorescence of the media. Uses a Gaussian process to fit the time-series of measurements of all replicates of the media and subtracts this time series from the raw data.

        Arguments
        --
        conditions: list of experimental conditions to be corrected
        datatypes: data types to be corrected (default is 'OD')
        mean: if True, the mean over time of the media values is used for the correction (rather than fitting with a Gaussian process)
        commonmedia: a condition containing a 'null' strain that should be used to correct other conditions
        experiments: list of experiments to include (default is 'all')
        conditions: list of conditions to include (default is 'all')
        experimentincludes: selects only experiments with experimentincludes in their name
        experimentexcludes: ignores experiments with experimentexcludes in their name
        conditionincludes: selects only conditions with conditionincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        figs: if True, display fits of the media
        noruns: number of attempts used to fit the media
        exitearly: if True, stop at the first successful fit otherwise take the best fit from all successful fits
        bd: changes the limits on the hyperparameters for the Gaussian process, which has a squared exponential kernel, used to fit the media
        results: if True, display best-fit parameters
        log: if True, call will be added to log
        iskip: use only every iskip'th data point to increase speed (must be an integer)
        '''
        if log: self.logmethod(self.logger)
        exps= self.getexps(experiments, experimentincludes, experimentexcludes)
        cons= self.getcons(conditions, conditionincludes, conditionexcludes, nomedia= False)
        for exp in exps:
            # data types
            if datatypes == 'all':
                datatypes= self.datatypes[exp]
            else:
                datatypes= gu.makelist(datatypes)
            if self.progress['ODcorrectedformedia'][exp] and 'OD' in datatypes:
                datatypes.pop(datatypes.index('OD'))
            # correct for media
            for dtype in datatypes:
                for c in cons:
                    print(exp, ':: Correcting for media for', dtype, 'in', c)
                    cm= commonmedia if commonmedia else c
                    self.performmediacorrection(dtype, exp, c, figs, noruns, exitearly, bd, results, mean, cm, iskip)
                    self.sc.loc[(self.sc['condition'] == c) & (self.sc['experiment'] == exp),
                                dtype + ' corrected for media']= True
            if self.progress['negativevalues'][exp]:
                print('Warning: correcting media has created negative values in', exp, 'for')
                print(self.progress['negativevalues'][exp])
        # remake summary data
        self.update_s()



    #####
    def performmediacorrection(self, dtype, exp, condition, figs, noruns, exitearly, bd, results,
                               mean, commonmedia, iskip):
        '''
        Internal function: Uses a Gaussian process to fit the media over time and subtracts the best-fit values from the data.
        '''
        # find data for correction
        df= self.r.query("experiment == @exp and condition == @commonmedia and strain == 'null'")
        if df.empty:
            # no data
            print('No well annotated "null" was found for', commonmedia, 'in experiment', exp)
            print('Correcting for media abandoned!')
        else:
            # there is data
            t, data= df['time'].values, df[dtype].values
            r= self.r
            rtest= ((r['experiment'] == exp) & (r['condition'] == condition))
            # find correction
            if mean:
                corr= np.mean(data)
            else:
                if commonmedia in self.progress['mediaGP'][exp][dtype]:
                    # Gaussian process has already been fit
                    gm= self.progress['mediaGP'][exp][dtype][commonmedia]
                else:
                    # hyperparameters for Gaussian process
                    b= {0: (-6,4), 1: (-5,-2), 2: (-10,0)}
                    if bd: b= gu.mergedicts(original= b, update= bd)
                    # fit with Gaussian process
                    try:
                        if iskip:
                            t= t[::iskip]
                            data= data[::iskip]
                        f, fvar, gm= gu.smoothGP(t, data, bd= b, results= results)
                        self.progress['mediaGP'][exp][dtype][commonmedia]= gm
                        if figs:
                            plt.figure()
                            plt.plot(t, data, 'ro', np.sort(t), f, 'b-')
                            plt.xlabel('time (hours)')
                            plt.title('media correction for ' + dtype + ' in ' + condition)
                            plt.show(block= False)
                    except gp.gaussianprocessException:
                        raise(SystemExit('Fitting media failed'))
                # use GP to find correction
                gm.batchpredict(r[rtest]['time'])
                corr= gm.f
            # perform correction
            r.loc[rtest, dtype]= r[rtest][dtype] - corr
            # check for any negative values
            for s in np.unique(r[rtest]['strain'][r[rtest][dtype] < 0]):
                if s != 'null':
                    wstr= '\t' + dtype + ': ' + s + ' in ' + condition + '\n'
                    if not self.progress['negativevalues'][exp]:
                        self.progress['negativevalues'][exp]= wstr
                    else:
                        self.progress['negativevalues'][exp] += wstr


    #####
    # Statistical analysis
    #####
    def getstats(self, dtype= 'OD', bd= False, cvfn= 'matern', esterrs= False,
                 noruns= 10, exitearly= True, noinits= 100, nosamples= 100, iskip= False,
                 stats= True, plotodgr= False, figs= True,
                 experiments= 'all', experimentincludes= False, experimentexcludes= False,
                 conditions= 'all', conditionincludes= False, conditionexcludes= False,
                 strains= 'all', strainincludes= False, strainexcludes= False):
        '''
        Calls fitderiv.py to estimate the first and second time-derivatives of 'OD' (default) and corresponding summary statistics.
        Time-variables are stored in the self.s dataframe; summary statistics are stored in the self.sc dataframe

        Use, for example, p.s.columns to see what variables have been calculated.

        Arguments
        --
        dtype: type of data ('OD' - default, or, for example, 'GFP', 'c-GFPperod', or 'c-GFP')
        bd: the limits on the hyperparameters for the Gaussian process used for the fit. For example, p.getstats('1% Gal', 'GAL2', bd= {1: [-2,-2])}) fixes the flexibility to be 0.01
        cvfn: covariance function used for fit, either 'matern' (default) or 'sqexp' or 'nn' or, for example, 'sqexp : matern' to pick the covariance function with the highest maximum likelihood
        esterrs: if True, measurement errors are empirically estimated from the variance across replicates at each time point; if False, the size of the measurement error is fit from the data assuming that this size is the same at all time points
        noruns: number of attempts made for each fit (default is 5), each run is made with random initial estimates of the hyperparameters
        exitearly: if True, stop at the first successful fit otherwwise take the best fit from all successful fits
        noinits: number of attempts at finding the initial condition for the optimization
        nosamples: number of bootstrap samples to calculate errors in statistics
        iskip: use only every iskip'th data point to increase speed (must be an integer)
        stats: if False, do not calculate statistics
        plotodgr: for OD data, plots growth rate versus OD if True (default is False)
        figs: if True (default), shows the fit and inferred derivative
        experiments: list of experiments to include (default is 'all')
        conditions: list of conditions to include (default is 'all')
        strains: list of strains to include (default is 'all')
        experimentincludes: selects only experiments with experimentincludes in their name
        experimentexcludes: ignores experiments with experimentexcludes in their name
        conditionincludes: selects only conditions with conditionincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        strainincludes: selects only strains with strainincludes in their name
        strainexcludes: ignore strains with strainexcludes in their name
        '''
        self.logmethod(self.logger)
        linalgmax= 5
        exps, cons, strs= self.getall(experiments, experimentincludes, experimentexcludes,
                                      conditions, conditionincludes, conditionexcludes,
                                      strains, strainincludes, strainexcludes)
        for e in exps:
            for c in cons:
                for s in strs:
                    figtitle= e + ':: ' + s + ' in ' + c
                    print('\nFitting', dtype, 'for', figtitle)
                    if dtype in self.r.columns:
                        # raw data
                        d= self.extractwells(e, c, s, dtype)
                        if d.size == 0: break
                    elif dtype in self.s.columns:
                        # processed data
                        df= self.s.query('experiment == @e and condition == @c and strain == @s')
                        d= df[dtype].values
                        if d.size == 0: break
                    else:
                        print(dtype, 'not recognized for', figtitle)
                        return
                    t= self.s.query('experiment == @e and condition == @c and strain == @s')['time'].values
                    if dtype == 'OD':
                        snames= ['max gr', 'time of max gr', 'doubling time', 'max OD', 'lag time']
                        ylabels= ['log(OD)', 'gr']
                        logs= True
                    else:
                        esterrs= df[dtype + ' var'].values
                        snames= ['max d/dt of ' + dtype,
                                 'time of max d/dt of ' + dtype,
                                 'inverse of max d/dt of ' + dtype,
                                 'max ' + dtype, 'lag time of ' + dtype]
                        ylabels= [dtype, 'd/dt ' + dtype]
                        logs= False
                    if ' : ' in cvfn:
                        # multiple covariance functions
                        cvfns= cvfn.split(' : ')
                        fs= [fitderiv(t, d, figs= False, cvfn= cfn, logs= logs, bd= bd, esterrs= esterrs,
                                      statnames= snames, noruns= noruns, noinits= noinits, exitearly= exitearly,
                                      linalgmax= linalgmax, nosamples= nosamples, iskip= iskip)
                             for cfn in cvfns]
                        # pick the covariance function with the highest maximum likelihood
                        ci= np.argmax([f.logmaxlike for f in fs])
                        f= fs[ci]
                        print(' Using', cvfns[ci])
                        if figs:
                            plt.figure()
                            for i, cfn in enumerate(cvfns):
                                ax= plt.subplot(2, len(cvfns), i+1)
                                fs[i].plotfit('f', ylabel= ylabels[0], figtitle= figtitle)
                                ax.set_title(cfn.upper()) if i == ci else ax.set_title(cfn)
                                plt.subplot(2, len(cvfns), i+1+len(cvfns))
                                fs[i].plotfit('df', ylabel= ylabels[1])
                            plt.suptitle(figtitle, fontweight='bold')
                            plt.show(block= False)
                    else:
                        # single covariance function
                        f= fitderiv(t, d, figs= False, cvfn= cvfn, logs= logs, bd= bd, esterrs= esterrs,
                                    statnames= snames, noruns= noruns, noinits= noinits, exitearly= exitearly,
                                    linalgmax= linalgmax, nosamples= nosamples, iskip= iskip)
                        if figs:
                            plt.figure()
                            plt.subplot(2,1,1)
                            f.plotfit('f', ylabel= ylabels[0], figtitle= figtitle)
                            plt.subplot(2,1,2)
                            f.plotfit('df', ylabel= ylabels[1])
                            plt.tight_layout()
                            plt.show(block= False)
                    if dtype == 'OD':
                        outdf= pd.DataFrame({'experiment' : e, 'condition' : c, 'strain' : s, 'time' : t,
                                             'flogOD' : f.f, 'flogOD var' : f.fvar,
                                             'gr' : f.df, 'gr var' : f.dfvar,
                                             'd/dt gr' : f.ddf, 'd/dt gr var' : f.ddfvar})
                        self.progress['getstatsGP'][e][c][s]= f
                        if plotodgr:
                            gu.plotxyerr(np.exp(f.f), f.df,
                                        (np.exp(f.f + np.sqrt(f.fvar)) - np.exp(f.f - np.sqrt(f.fvar)))/2.0,
                                        np.sqrt(f.dfvar), 'OD', 'growth rate', figtitle + ' : growth rate vs OD')

                        # find summary statistics
                        fs, gs, hs= f.fitderivsample(nosamples)
                        # log2 OD ratio
                        dr= np.log2(np.exp(fs[-1,:] - fs[0,:]))
                        # find local maximum derivative
                        from scipy.signal import argrelextrema
                        da, dt= [], []
                        for gsample in np.transpose(gs):
                            tpts= argrelextrema(gsample, np.greater)[0]
                            if np.any(tpts):
                                da.append(np.max(gsample[tpts]))
                                dt.append(f.t[tpts[np.argmax(gsample[tpts])]])
                        # find area under gr vs OD
                        from scipy import integrate, interpolate
                        aa, ana= [], []
                        for fsample, gsample in zip(np.transpose(fs), np.transpose(gs)):
                            sod= np.exp(fsample)
                            integrand= lambda x: interpolate.interp1d(sod, gsample)(x)
                            iresult= integrate.quad(integrand, np.min(sod), np.max(sod),
                                                     limit= 100, full_output= 1)[0]
                            aa.append(iresult)
                            ana.append(iresult/(np.max(sod) - np.min(sod)))
                        # store results
                        statsdict= {'experiment' : e, 'condition' : c, 'strain' : s,
                                   'log2 OD ratio' : np.mean(dr), 'log2 OD ratio var' : np.var(dr),
                                   'local max gr' : np.mean(da), 'local max gr var' : np.var(da),
                                   'time of local max gr' : np.mean(dt),
                                   'time of local max gr var' : np.var(dt),
                                   'area under gr vs OD' : np.mean(aa),
                                   'area under gr vs OD var' : np.var(aa),
                                   'normalized area under gr vs OD' : np.mean(ana),
                                   'normalized area under gr vs OD var' : np.var(ana)}
                    else:
                        outdf= pd.DataFrame({'experiment' : e, 'condition' : c, 'strain' : s, 'time' : t,
                                             'f' + dtype : f.f, 'f' + dtype + ' var' : f.fvar,
                                             'd/dt ' + dtype : f.df, 'd/dt ' + dtype + ' var' : f.dfvar,
                                             'd2/dt2 ' + dtype : f.ddf, 'd2/dt2 ' + dtype + ' var' : f.ddfvar})
                        statsdict= {'experiment' : e, 'condition' : c, 'strain' : s}

                    # store results in global dataframes
                    statsdict[dtype + ' logmaxlike']= f.logmaxlike
                    statsdict[dtype + ' gp']= cvfn
                    if stats:
                        for sname in f.ds.keys(): statsdict[sname]= f.ds[sname]
                    # add growth rates, etc., to dataframe of summary data
                    if (dtype == 'OD' and 'gr' not in self.s.columns) or (dtype != 'OD' and 'f' + dtype not in self.s.columns):
                        # add new columns to dataframe
                        self.s= pd.merge(self.s, outdf, how= 'outer')
                    else:
                        # update dataframe
                        self.s= gu.absorbdf(self.s, outdf, ['experiment', 'condition', 'strain', 'time'])
                    # create or add summary stats to stats dataframe
                    statsdf= pd.DataFrame(statsdict, index= pd.RangeIndex(0, 1, 1))
                    if (dtype  + ' logmaxlike' not in self.sc.columns):
                        # add new columns to dataframe
                        self.sc= pd.merge(self.sc, statsdf, how= 'outer')
                    else:
                        # update dataframe
                        self.sc= gu.absorbdf(self.sc, statsdf, ['experiment', 'condition', 'strain'])



    #####
    def getfitnesspenalty(self, ref, com, y= 'gr', abs= False, figs= True, nosamples= 100, norm= False):
        '''
        Calculates the area between typically two growth rate versus OD curves, normalized by the length
        along the OD-axis where they overlap:

        e.g. p.getfitnesspenalty(['1.9% Raffinose 0.0µg/ml cycloheximide', '77.WT'],
                             ['1.9% Raffinose 0.5µg/ml cycloheximide', '77.WT'])

        Arguments
        --
        ref: [condition, strain] array for the reference or [experiment, condition, strain] if more than one experiment
        com: [condition, strain] array to compare to the reference or [experiment, condition, strain] if more than one experiment
        y: variable to be compared (default: 'gr')
        figs: if True, an example of the area between the curves is shown
        nosamples: determine the number bootstraps to estimate error (default: 100)
        norm: if True, returns the mean and variance of the area under the reference strain for normalisation
        '''
        self.logmethod(self.logger)
        if len(self.allexperiments) == 1:
            ref.insert(0, self.allexperiments[0])
            com.insert(0, self.allexperiments[0])
        # get and sample from Gaussian processes
        if nosamples and y == 'gr':
            # estimate errors
            try:
                # sample from Gaussian process
                f0s, g0s, h0s= self.progress['getstatsGP'][ref[0]][ref[1]][ref[2]].fitderivsample(nosamples)
                f1s, g1s, h1s= self.progress['getstatsGP'][com[0]][com[1]][com[2]].fitderivsample(nosamples)
                xsref, ysref= np.exp(f0s), g0s
                xscom, yscom= np.exp(f1s), g1s
            except KeyError:
                print("Failed: getstats('OD') needs to be run for these strains")
                return
        else:
            # no estimates of errors
            if nosamples: print("Cannot estimate errors - require y= 'gr' and a recently run getstats")
            xsref= self.s.query('experiment == @ref[0] and condition == @ref[1] and strain == @ref[2]')['OD mean'].values[:, None]
            ysref= self.s.query('experiment == @ref[0] and condition == @ref[1] and strain == @ref[2]')[y].values[:, None]
            xscom= self.s.query('experiment == @com[0] and condition == @com[1] and strain == @com[2]')['OD mean'].values[:, None]
            yscom= self.s.query('experiment == @com[0] and condition == @com[1] and strain == @com[2]')[y].values[:, None]
            if xsref.size == 0 or ysref.size == 0:
                print(ref[0] + ': Data missing for', ref[2], 'in', ref[1])
                return np.nan, np.nan
            elif xscom.size == 0 or yscom.size == 0:
                print(com[0] + ': Data missing for', com[2], 'in', com[1])
                return np.nan, np.nan
        fps= np.zeros(xsref.shape[1])
        nrm= np.zeros(xsref.shape[1])
        samples= zip(np.transpose(xsref), np.transpose(ysref), np.transpose(xscom), np.transpose(yscom))
        # process samples
        for j, (xref, yref, xcom, ycom) in enumerate(samples):
            # remove any double values in OD because of OD plateau'ing
            uxref, uiref= np.unique(xref, return_inverse= True)
            uyref= np.array([np.median(yref[np.nonzero(uiref == i)[0]]) for i in range(len(uxref))])
            uxcom, uicom= np.unique(xcom, return_inverse= True)
            uycom= np.array([np.median(ycom[np.nonzero(uicom == i)[0]]) for i in range(len(uxcom))])
            # interpolate data
            from scipy import interpolate
            iref= interpolate.interp1d(uxref, uyref, fill_value= 'extrapolate', kind= 'slinear')
            icom= interpolate.interp1d(uxcom, uycom, fill_value= 'extrapolate', kind= 'slinear')
            # find common range of x
            uxi= np.max([uxref[0], uxcom[0]])
            uxf= np.min([uxref[-1], uxcom[-1]])
            # perform integration to find normalized area between curves
            from scipy import integrate
            if abs:
                igrand= lambda x: np.abs(iref(x) - icom(x))
            else:
                igrand= lambda x: iref(x) - icom(x)
            fps[j]= integrate.quad(igrand, uxi, uxf, limit= 100, full_output= 1)[0]/(uxf - uxi)
            if norm:
                # calculate area under curve of reference strain as a normalisation
                igrand= lambda x: iref(x)
                nrm[j]= integrate.quad(igrand, uxi, uxf, limit= 100, full_output= 1)[0]/(uxf - uxi)
            # an example figure
            if figs and j == 0:
                plt.figure()
                plt.plot(uxref, uyref, 'k-', uxcom, uycom, 'b-')
                x= np.linspace(uxi, uxf, np.max([len(uxref), len(uxcom)]))
                plt.fill_between(x, iref(x), icom(x), facecolor= 'red', alpha= 0.5)
                plt.xlabel('OD')
                plt.ylabel(y)
                plt.legend([ref[0] + ':: ' + ref[2] + ' in ' + ref[1],
                            com[0] + ':: ' + com[2] + ' in ' + com[1]],
                           loc= 'upper left', bbox_to_anchor= (-0.05, 1.15))
                plt.show(block= False)
        if norm:
            return np.mean(fps), np.var(fps), np.mean(nrm), np.var(nrm)
        else:
            return np.mean(fps), np.var(fps)





    #####
    # Fluorescence corrections
    #####
    def correctauto(self, f= ['GFP', 'AutoFL'], refstrain= 'WT', figs= True, noruns= 2,
                    bd= False, no1samples= 100, results= False, minqrerr= 1.0e-6,
                    experiments= 'all', experimentincludes= False, experimentexcludes= False,
                    conditions= 'all', conditionincludes= False, conditionexcludes= False,
                    strains= 'all', strainincludes= False, strainexcludes= False):
        '''
        Corrects fluorescence data for autofluorescence in comparison with a reference strain

        Arguments
        --
        f: fluorescence measurements to be used, typically either ['mCherry'] or  ['GFP', 'AutoFL'] (default)
        refstrain: the name of the reference strain ('WT' is the default)
        figs: if True (default), display fits used in applying the correction
        noruns: number of fit attempts to be used
        bd: the limits on the hyperparameters for the Gaussian process, which has a squared exponential kernel and is used for fitting data from the reference strain
        no1samples: number of samples used to estimate errors when correcting fluorescence data measured at one wavelength
        results: if True, display best-fit parameters
        minqrerr: minimum value allowed for the estimated error in the ratio of fluorescence (AutoFL/GFP) - too small values can cause instabilities in the fitting
        experiments: list of experiments to include (default is 'all')
        conditions: list of conditions to include (default is 'all')
        strains: list of strains to include (default is 'all')
        experimentincludes: selects only experiments with experimentincludes in their name
        experimentexcludes: ignores experiments with experimentexcludes in their name
        conditionincludes: selects only conditions with conditionincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        strainincludes: selects only strains with strainincludes in their name
        strainexcludes: ignore strains with strainexcludes in their name
        '''
        self.logmethod(self.logger)
        f= gu.makelist(f)
        print('Using', refstrain, 'as the reference')
        # check have enough replicates
        if len(f) == 2:
            exps, cons, strs= self.getall(experiments, experimentincludes, experimentexcludes,
                                          conditions, conditionincludes, conditionexcludes,
                                          strains, strainincludes, strainexcludes)
            for e in exps:
                for c in cons:
                    for s in strs:
                        df= self.r.query("experiment == @e and condition == @c and strain == @s")
                        if df['well'].unique().size < 2:
                            print(e, ':: Not enough replicates to correct autofluorescence for', s, 'in', c)
                            print('Try specifying just one fluorescence measurement')
                            return
            self.correctauto2(f, refstrain, figs, noruns, bd, results, minqrerr,
                              experiments, experimentincludes, experimentexcludes,
                              conditions, conditionincludes, conditionexcludes,
                              strains, strainincludes, strainexcludes)
        elif len(f) == 1:
            self.correctauto1(f, refstrain, figs, noruns, bd, no1samples, results,
                              experiments, experimentincludes, experimentexcludes,
                              conditions, conditionincludes, conditionexcludes,
                              strains, strainincludes, strainexcludes)
        else:
            print('f must be a list of length 1 or 2')



    #####
    def correctauto1(self, f, refstrain, figs, noruns, bd, nosamples, results,
                     experiments, experimentincludes, experimentexcludes,
                     conditions, conditionincludes, conditionexcludes,
                     strains, strainincludes, strainexcludes):
        '''
        Internal function: Corrects for autofluorescence for experiments with measured emissions at one wavelength using the fluorescence of the wild-type interpolated to the OD of the tagged strain.
        '''
        # correct autofluorescence
        print('Correcting autofluorescence')
        # run through all experiments
        for e in self.getexps(experiments, experimentincludes, experimentexcludes):
            # run through all conditions
            for c in self.getcons(conditions, conditionincludes, conditionexcludes, nomedia= True):
                # process reference strain
                if f[0] not in self.progress['refGP'][e][c]:
                    self.processref1(f, refstrain, figs, noruns, bd, results, e, c)
                # GP for reference strain
                gfr, interpf, fmerr= self.progress['refGP'][e][c][f[0]]
                t= self.s.query('experiment == @e and condition == @c and strain == @refstrain')['time'].values
                for s in self.getstrs(strains, strainincludes, strainexcludes, nonull= True):
                    od, fdata= self.extractwells(e, c, s, ['OD', f[0]])
                    if od.size == 0: break
                    nr= fdata.shape[1]
                    fl, fl2, flperod, flperod2= (0,)*4
                    # run over all replicates
                    for i in range(nr):
                        # interpolate WT fluorescence errors to OD values of strain
                        try:
                            men= interpf(od[:,i])
                        except ValueError:
                            men= np.median(fmerr)*np.ones(fmerr.size)
                        # predict WT fluorescence at OD values of strain
                        gfr.predict(od[:,i], merrorsnew= men, addnoise= True)
                        # sample for estimating errors
                        fs= gfr.sample(nosamples)
                        if self.progress['ODcorrected'][e]:
                            # sample ODs corrected for non-linearities
                            self.progress['gc'][e].predict(od[:,i])
                            if False:
                                # gives high errors in flperod because sampled corrected ODs can be close to zero
                                cods= self.progress['gc'][e].sample(nosamples)
                            else:
                                cods= gu.tilec(self.progress['gc'][e].f, nosamples)
                        else:
                            cods= gu.tilec(od[:,i], nosamples)
                        # find fluorescence per cell
                        for j in range(nosamples):
                            d= fdata[:,i] - fs[:,j]
                            dperod= d/cods[:,j]
                            fl += d
                            fl2 += d**2
                            flperod += dperod
                            flperod2 += dperod**2
                    bname= 'c-' + f[0]
                    autofdict= {'experiment' : e, 'condition' : c, 'strain' : s, 'time' : t,
                                bname : fl/(nr*nosamples),
                                bname + ' var' : fl2/(nr*nosamples) - (fl/(nr*nosamples))**2,
                                bname + 'perod' : flperod/(nr*nosamples),
                                bname + 'perod var' : flperod2/(nr*nosamples) - (flperod/(nr*nosamples))**2}
                    autofdf= pd.DataFrame(autofdict)
                    if bname not in self.s.columns:
                        # extend dataframe
                        self.s= pd.merge(self.s, autofdf, how= 'outer')
                    else:
                        # update dataframe
                        self.s= gu.absorbdf(self.s, autofdf, ['experiment', 'condition', 'strain', 'time'])
                    self.sc.loc[(self.sc['experiment'] == e) & (self.sc['condition'] == c)
                                &  (self.sc['strain'] == s), f[0] + ' corrected for autofluorescence']= True



    #####
    def processref1(self, f, refstrain, figs, noruns, bd, results, experiment, condition):
        '''
        Internal function: Processes reference strain for data with one fluorescence measurement. Uses a Gaussian process to fit the fluorescence as a function of OD.

        Arguments
        --
        f: fluorescence to be corrected, such as ['mCherry']
        refstrain: the name of the reference strain, such as 'WT'
        figs: if True, display fits of the fluorescence
        noruns: number of fit attempts used
        bd: to change the limits on the hyperparameters for the Gaussian process, which has a squared exponential kernel and is used to fit the fluorescence
        results: if True, display best-fit parameters
        experiment: experiment to be corrected
        condition: condition to be corrected
        '''
        e, c= experiment, condition
        # bounds for Gaussian process for fitting reference strain
        b= {0: (3,10), 1: (-4,4), 2: (-4,2)}
        if bd: b= gu.mergedicts(original= b, update= bd)

        print(e + ':: Processing reference strain', refstrain, 'for', f[0], 'in', c)
        # fit reference strain's fluorescence as a function of OD
        df= self.r.query('experiment == @e and condition == @c and strain == @refstrain')
        if df.empty:
            print(e + ':: ' + refstrain, 'not found in', c)
            raise(SystemExit('Running correctauto failed'))
        else:
            x, y= df['OD'].values, df[f[0]].values
            # convert y into a matrix with one column per well
            ym= self.extractwells(e, c, refstrain, f[0])
            # smooth ym and set up Gaussian process
            me= gu.findsmoothvariance(ym)
            ys= y[np.argsort(x)]
            mes= np.tile(me, ym.shape[1])[np.argsort(x)]
            xs= np.sort(x)
            # fit with Gaussian process
            gfr= gp.sqexpGP(b, xs, ys, merrors= mes)
            try:
                gfr.findhyperparameters(noruns= noruns)
                if results: gfr.results()
                gfr.predict(xs)
            except gp.gaussianprocessException:
                raise(SystemExit('Fitting reference strain failed'))
            if figs:
                # plot fit
                plt.figure()
                gfr.sketch('o')
                plt.xlabel('OD')
                plt.ylabel(f[0])
                plt.title(e + ':: fitting ' + refstrain + ' for ' + c)
                plt.show(block= False)
            # interpolate WT fluorescence errors to OD values of strain
            from scipy.interpolate import interp1d
            interpf= interp1d(self.s.query('experiment == @e and condition == @c and strain == @refstrain')
                              ['OD mean'].values, me)
            # store reference strain information
            self.progress['refGP'][e][c][f[0]]= (gfr, interpf, me)





    #####
    def correctauto2(self, f, refstrain, figs, noruns, bd, results, minqrerr,
                     experiments, experimentincludes, experimentexcludes,
                     conditions, conditionincludes, conditionexcludes,
                     strains, strainincludes, strainexcludes):
        '''
        Internal function: Corrects for autofluorescence using spectral unmixing for experiments with measured emissions at two wavelengths (following Lichten et al.)
        '''
        # correct for autofluorescence
        print('Correcting autofluorescence')
        for e in self.getexps(experiments, experimentincludes, experimentexcludes):
            for c in self.getcons(conditions, conditionincludes, conditionexcludes, nomedia= True):
                # process reference strain
                if f[0] not in self.progress['refGP'][e][c]:
                    self.processref2(f, refstrain, figs, noruns, bd, results, minqrerr, e, c)
                # GP for reference strain
                gr, qrerr= self.progress['refGP'][e][c][f[0]]
                t= self.s.query('experiment == @e and condition == @c and strain == @refstrain')['time'].values
                for s in self.getstrs(strains, strainincludes, strainexcludes, nonull= True):
                    if s != refstrain:
                        f0, f1= self.extractwells(e, c, s, f)
                        if f0.size == 0 or f1.size == 0: break
                        nodata, nr= f0.shape
                        # remove autofluorescence
                        fl= self.applyautoflcorrection(gr.f, f0, f1)
                        # estimate error
                        varf= np.var(f0, 1)
                        varcf= np.var(f1, 1)
                        gr.predict(t, merrorsnew= qrerr)
                        rs= np.transpose(gr.sample(self.nosamples))
                        varfl= np.zeros(nodata)
                        # average variance over samples of ra
                        for ra in rs:
                            varfl += (varf + ra**2*varcf)/(self.gamma - ra)**2
                        varfl /= self.nosamples
                        bname= 'c-' + f[0]
                        odmean= self.s.query('experiment == @e and condition == @c and strain == @s')['OD mean'].values
                        autofdict= {'experiment' : e, 'condition' : c,'strain' : s, 'time' : t,
                                    bname : np.mean(fl,1),
                                    bname + ' var' : varfl,
                                    bname + 'perod' : np.mean(fl,1)/odmean}
                        flperodvar= np.zeros(nodata)
                        od= self.extractwells(e, c, s, 'OD')
                        for i in np.random.randint(nr, size= self.nosamples):
                            if self.progress['ODcorrected'][e]:
                                self.progress['gc'][e].predict(od[:,i])
                                cc= self.progress['gc'][e].sample(1)
                                cc= cc.flatten('F')
                                flperodvar += autofdict[bname + ' var']/cc**2
                            else:
                                flperodvar += autofdict[bname + ' var']/od[:,i]**2
                        flperodvar /= self.nosamples
                        autofdict[bname + 'perod var']= flperodvar
                        # calculate corrected levels that are above the expectation from reference strain
                        refdf= self.s.query('experiment == @e and condition == @c and strain == @refstrain')
                        reflevel= self.consist*np.sqrt(refdf[bname + ' var'].values)
                        keep= np.nonzero(refdf[bname].values > reflevel)[0]
                        autofdict['s' + bname]= np.zeros(nodata)
                        autofdict['s' + bname][keep]= autofdict[bname][keep]
                        autofdict['s' + bname + ' var']= np.zeros(nodata)
                        autofdict['s' + bname + ' var'][keep]= autofdict[bname + ' var'][keep]
                        reflevel /= self.s.query('experiment == @e and condition == @c and strain == @refstrain')['OD mean'].values
                        keep= np.nonzero(autofdict[bname+ 'perod'] > reflevel)[0]
                        autofdict['s' + bname + 'perod']= np.zeros(nodata)
                        autofdict['s' + bname + 'perod'][keep]= autofdict[bname + 'perod'][keep]
                        autofdict['s' + bname + 'perod var']= np.zeros(nodata)
                        autofdict['s' + bname + 'perod var'][keep]= autofdict[bname + 'perod var'][keep]
                        # add to dataframe
                        autofdf= pd.DataFrame(autofdict)
                        self.s= gu.absorbdf(self.s, autofdf, ['experiment', 'condition', 'strain', 'time'])
                        self.sc.loc[(self.sc['experiment'] == e) & (self.sc['condition'] == c)
                                    & (self.sc['strain'] == s), f[0] + ' corrected for autofluorescence']= True



    #####
    def processref2(self, f, refstrain, figs, noruns, bd, results, minqrerr, experiment, condition):
        '''
        Internal function: Processes reference strain data for spectral unmixing (for experiments with two fluorescence measurements). Uses a Gaussian process to fit the ratio of emitted fluorescence measurements and checks that reference strain data is itself corrected to zero.

        Arguments
        --
        f: fluorescence to be corrected, such as ['GFP', 'AutoFL']
        refstrain: the name of the reference strain, such as 'WT'
        figs: if True, display fits of the ratios and, as a check, the correction of the reference strain by itself
        noruns: number of attempts used to fit the ratio of fluorescences
        bd: the limits on the hyperparameters for the Gaussian process, which has a squared exponential kernel and is used to fit the ratio of fluorescences
        results: if True, display best-fit parameters
        minqrerr: if values for the estimated error in the fluorescence ratio fall below this value replace by this minimum value
        experiment: experiment to be corrected
        condition: condition to be corrected
        '''
        e, c= experiment, condition
        # bounds for Gaussian process for fitting reference strain
        b= {0: (-5,3), 1: (-4,-1), 2: (-4, 4)}
        if bd: b= gu.mergedicts(original= b, update= bd)

        print(e + ':: Processing reference strain', refstrain, 'for', f[0], 'in', c)
        df= self.r.query('experiment == @e and condition == @c and strain == @refstrain')
        if df.empty:
            print(e + ':: ' + refstrain, 'not found in', c)
            raise(SystemExit('Running correctauto failed'))
        f0, f1= self.extractwells(e, c, refstrain, f)
        nr= f0.shape[1]
        t= self.s.query('experiment == @e and condition == @c and strain == @refstrain')['time'].values
        # find values with zero variance
        dels= np.unique(np.nonzero(np.var(f1/f0, 1) == 0)[0])
        # remove offending values
        f0r= np.delete(f0, dels, axis= 0)
        f1r= np.delete(f1, dels, axis= 0)
        tr= np.delete(t, dels)
        # find ratio of fluorescence
        qr= (f1r/f0r).flatten('F')
        # find errors
        qrerr= np.var(f1r/f0r, 1)
        # check no errors too small
        if np.min(qrerr) < minqrerr:
            print('Warning: replacing small estimates for the error in the fluorescence ratio')
            qrerr[qrerr < minqrerr]= minqrerr
        # sort
        x= np.tile(tr, nr)
        xs= np.sort(x)
        qrs= qr[np.argsort(x)]
        qrerrs= np.repeat(qrerr, nr)[np.argsort(x)]
        # fit with Gaussian process
        gr= gp.sqexpGP(b, xs, qrs, merrors= qrerrs)
        try:
            # optimize hyperparameters of Gaussian process
            gr.findhyperparameters(noruns)
            if results: gr.results()
            # add back missing time points
            if np.any(dels): qrerr= np.interp(t, tr, qrerr)
            gr.predict(t, merrorsnew= qrerr)
        except gp.gaussianprocessException:
            raise(SystemExit('Fitting reference strain failed'))
        # store results
        self.progress['refGP'][e][c][f[0]]= gr, qrerr
        if figs:
            # plot fit
            plt.figure()
            gr.sketch('o')
            plt.xlabel('time (hours)')
            plt.ylabel(f[1] + '/' + f[0])
            plt.title('fitting ' + refstrain + ' for ' + c)
            plt.show(block= False)
        # check autofluorescence correction for reference strain
        flref= self.applyautoflcorrection(gr.f, f0, f1)
        odmean= self.s.query('experiment == @e and condition == @c and strain == @refstrain')['OD mean'].values
        bname= 'c-' + f[0]
        autofdict= {'experiment' : e, 'condition' : c, 'strain' : refstrain, 'time' : t,
                    bname : np.mean(flref, 1),
                    bname + 'perod' : np.mean(flref/odmean[:,None], 1),
                    bname + ' var' : np.var(flref, 1),
                    bname + 'perod var' : np.var(flref/odmean[:,None], 1)}
        if bname not in self.s.columns:
            self.s= pd.merge(self.s, pd.DataFrame(autofdict), how= 'outer')
        else:
            self.s= gu.absorbdf(self.s, pd.DataFrame(autofdict), ['experiment', 'condition', 'strain', 'time'])
        if figs:
            # plot correction for reference strain
            plt.figure()
            plt.plot(t, flref, '.')
            plt.plot(t, self.consist*np.sqrt(autofdict[bname + ' var']), 'r:')
            plt.plot(t, -self.consist*np.sqrt(autofdict[bname + ' var']), 'r:')
            plt.plot(t, np.zeros(np.size(t)), 'k')
            plt.ylabel('corrected ' + refstrain + ' fluorescence')
            plt.xlabel('time (hours)')
            plt.title(e + ':: ' + c + ': consistency check for reference strain ' + refstrain)
            plt.show(block= False)



    #####
    def applyautoflcorrection(self, ra, fdata, cfdata):
        '''
        Internal function: Corrects for autofluorescence returning an array of replicates.
        '''
        nr= fdata.shape[1]
        raa= np.reshape(np.tile(ra, nr), (np.size(ra), nr), order= 'F')
        return (raa*fdata - cfdata)/(raa - self.gamma*np.ones(np.shape(raa)))




    #####
    # Logging
    #####
    def log(self):
        '''
        Prints log of all methods called and their arguments.
        '''
        print(self.logstream.getvalue())


    def logmethod(self, logger):
        '''
        Internal function: logs a method and its arguments
        '''
        currframe= inspect.currentframe()
        # find frame of calling routine
        frame= inspect.getouterframes(currframe)[1].frame
        # name of calling routine
        methodname= inspect.getframeinfo(frame)[2]
        # arguments of calling routine
        args, _, _, locals= inspect.getargvalues(frame)
        # add to log
        if methodname == '__init__':
            logstring= 'p= platereader('
        else:
            logstring= 'p.' + methodname + '('
        for arg in args:
            if 'self' not in arg:
                if type(locals[arg]) is str:
                    argstr= "'" + locals[arg] + "'"
                else:
                    argstr= str(locals[arg])
                logstring += arg + '= ' + argstr + ', '
        logstring= logstring[:-2] + ')\n'
        logger.info(logstring)



    #####
    # Exporting and importing
    #####
    def exportdf(self, commonname= False):
        '''
        Exports the data as json files.
        Dataframes for the (processed) raw data, for summary data, and for summary statistics and corrections, as well as a log file, will be exported.

        Arguments
        --
        commonname: the common name for the output files
        '''
        self.logmethod(self.logger)
        if commonname:
            commonname= self.wdir + commonname
        else:
            commonname= self.wdir + ''.join(self.allexperiments)
        self.r.to_json(commonname + '_r.json', orient= 'split')
        self.s.to_json(commonname + '_s.json', orient= 'split')
        self.sc.to_json(commonname + '_sc.json', orient= 'split')
        f= open(commonname + '.log', 'w')
        f.write(self.logstream.getvalue())
        f.close()
        print('Exported to', self.wdir)


    def importdf(self, commonnames, info= True):
        '''
        Import dataframes saved as json files.

        Arguments
        --
        commonnames: a list of common names for the json files with one common name for each experiment
        '''
        self.logmethod(self.logger)
        commonnames= gu.makelist(commonnames)
        for commonname in commonnames:
            commonname= self.wdir + commonname
            for df in ['r', 's', 'sc']:
                try:
                    exec("impdf= pd.read_json(commonname + '_' + df + '.json', orient= 'split')")
                    if hasattr(self, df):
                        exec("self." + df + "= pd.merge(self." + df + ", impdf, how= 'outer')")
                    else:
                        exec("self." + df + "= impdf")
                    print('Imported', commonname + '_' + df +'.json')
                except ValueError:
                    print('No file called', commonname + '_' + df + '.json found')
        # update attributes
        self.allexperiments= list(pd.unique(self.s['experiment']))
        self.allconditions= {e : list(pd.unique(self.s.loc[self.s['experiment'] == e]['condition']))
                              for e in self.allexperiments}
        self.allstrains= {e : list(pd.unique(self.s.loc[self.s['experiment'] == e]['strain']))
                              for e in self.allexperiments}
        # find datatypes with mean in self.s
        self.datatypes= {e : list(self.s.query('experiment == @e').columns[self.s.query('experiment == @e').columns.str.contains('mean')]) for e in self.allexperiments}
        self.datatypes= {e : [dt.split(' mean')[0] for dt in self.datatypes[e]] for e in self.datatypes}
        # display info on import
        if info: self.info()


#####

if __name__ == '__main__': print(platereader.__doc__)
