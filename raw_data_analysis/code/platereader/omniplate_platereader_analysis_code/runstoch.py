#!/usr/bin/env python
import sys, subprocess
from glob import glob
import numpy as np


matlab= False
for i, opt in enumerate(sys.argv[1:]):
    if '-' in opt:
        key= opt.split('-')[1]
        if key == 'm':
            model= sys.argv[i+2]
            rootmodel= model.split('.')[0]
        elif key == 't':
            skt= sys.argv[i+2]
        elif key ==  'r':
            skr= sys.argv[i+2]
        elif key ==  'i':
            ski= sys.argv[i+2]
        elif key == 'matlab':
            matlab= True
        else:
            print(opt + ' is not recognized.')

# run sencillo
subprocess.call(['sencillo.py', '-s', model])
# run stochkit
smodel= rootmodel + '.xml'
subprocess.call(['ssa', '-m', smodel, '-t', skt, '-r', skr, '-i', ski, '--keep-trajectories', '--no-stats', '-f'])

# get species
f= open(smodel)
for line in f:
    if '<!--' in line:
        species= line.split('--')[1].split()
f.close()

# create data structure
d= {}
d['t']= []
for s in species:
    d[s]= []

# get data
ddir= rootmodel + '_output/trajectories/'
tfiles= glob(ddir + '*.txt')
for i, tfile in enumerate(tfiles):
    f= open(tfile)
    for s in species:
        d[s].append([])
    for line in f:
        sline= line.split()
        if i == 0:
            d['t'].append(float(sline[0]))
        for j, s in enumerate(species):
            d[s][i].append(int(sline[j+1]))
    f.close()
# convert to numpy arrays
d['t']= np.asarray(d['t'])
for s in species:
    d[s]= np.transpose(np.asarray(d[s]))

# save
if matlab:
    import scipy.io as sio
    sio.savemat(rootmodel + '.mat', d)
else:
    import pickle
    pickle.dump(d, open(rootmodel + '.pkl', 'wb'))

# tidy up
subprocess.call(['rm', '-f', rootmodel + '.xml'])
subprocess.call(['rm', '-rf', rootmodel + '_output'])



