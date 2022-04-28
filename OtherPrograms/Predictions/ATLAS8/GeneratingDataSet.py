#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 08:53:55 2021

@author: vla18041
"""

###########################################################################################################################
# In ATLAS we are performing a first measurement ever of Z production at 8 TeV, double differential in transverse momentum and rapidity of the the boson, without any fiducial cuts on the leptons.

# We would like to know if you (or somebody else in your group) are interested in providing us with predictions (with their uncertainties) to be compared to our measurement.
# Ideally we would be interested in NNLL(‘)+NNLO and/or N3LL(‘)+N3LO predictions, as discussed in the context of the LHC precision EW working group benchmarking.

# The timescale for publication of the measurement is in the range of 2 to 4 months.

# Here are the details of the setup for the measurement:

# sqrt(s) = 8 TeV

# Cuts:
# 80 < mll < 100 GeV
# |yll| < 3.6
# No further cuts on leptons

# Binning: 2D in ptll and yll
# pt: 
# 0, 2.5, 5.0, 8.0, 11.4, 14.9, 18.5, 22.0, 25.5, 29.0,32.6,36.4, 40.4, 44.9, 50.2, 56.4, 63.9, 73.4, 85.4, 105.0, 132.0, 173.0, 253.0, 4000
# y:
# 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6
###########################################################################################################################

import sys
sys.path.append("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/")

import DataProcessor.Point
import DataProcessor.DataSet


pt=[0., 2.5, 5.0, 8.0, 11.4, 14.9, 18.5, 22.0, 25.5, 29.0,32.6,36.4,\
    40.4, 44.9, 50.2, 56.4, 63.9, 73.4, 85.4, 105.0, 132.0, 173.0, 253.0, 4000.]

y=[0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6]

DataCurrent=DataProcessor.DataSet.DataSet('A8-predict',"DY")
DataCurrent.comment="ATLAS 8TeV 2D bining (prediction)"
DataCurrent.reference="arXiv:????"

DataCurrent.isNormalized=False
proc_current=[1,1,5]
s_current=8000.**2
Q_current=[80.,100.]
incCut=False
cutParam=[0.,0.,-3.6,3.6]
lumUncertainty=0.018
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(y)-1):
    for j in range(len(pt)-1):
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(j+i*(len(pt)-1)))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=[pt[j],pt[j+1]]
        p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
        p["Q"]=Q_current
        p["<Q>"]=91.
        p["y"]=[y[i],y[i+1]]
        p["includeCuts"]=incCut
        p["cutParams"]=cutParam
        p["xSec"]=0.
        p["uncorrErr"].append(0.00001)
        DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/OtherPrograms/Predictions/ATLAS8/"+DataCurrent.name+".csv")

