#!/usr/bin/env python
from os import system
from sys import argv

if len(argv) == 4:
    datafile, oldsbmlfile, newsbmlfile= argv[1:4]
    system('java -jar tseries2sbml.jar SBMLAddTimeCourseData -csvIn ' + datafile
              + ' -sbmlIn ' + oldsbmlfile +' -sbmlOut ' + newsbmlfile)
else:
    cmd= argv[0]
    print('Usage:', cmd, 'datafile.csv originalSBMLmodel.sbml combinedSBMLmodel.sbml')
    print('Uses linear interpolation to convert a single time-varying input into SBML and then add to an existing model using the java package SBMLDataTools written by Alastair Hume and on GitHub.')
    print('datafile.csv must have two columns: time and the data.')
