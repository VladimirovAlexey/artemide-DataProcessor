#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 16 13:54:42 2025

@author: alexey
"""

#######################################
# importing libraries
#######################################
import os
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"

SNOWFLAKE_DIR = "/data/WorkingFiles/Twist3/Snowflake/PySnowflake/"


import sys
import numpy
sys.path.append(ROOT_DIR)
sys.path.append(SNOWFLAKE_DIR)

#%%
import DataProcessor.snowInterface
import DataProcessor.DataMultiSet
import SnowFlake

#%%
path_to_INI=ROOT_DIR+"FittingPrograms/Tw3_FIT/INI/TEST.ini"
SnowFlake.initialize(path_to_INI)

NP_par=numpy.zeros(18)+0.1
SnowFlake.setNPparameters(NP_par)
SnowFlake.UpdateTables(1.0, 25.0)

#%%
### read the list of files and return the list of DataSets
def loadThisDataD2(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=ROOT_DIR+"DataLib/D2_moment/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        
    return dataCollection

#%%
##################Cut function
def cutFunc(p):
    
    return True, p

#%%
### Loading the D2 data set
theData=DataProcessor.DataMultiSet.DataMultiSet("D2set",loadThisDataD2([
    "RQCD_d2_singlet", "RQCD_d2_ud",
    "SANE18"]))

setD2=theData.CutData(cutFunc) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setD2.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setD2.sets]), 'points.') 

#%%

DataProcessor.snowInterface.PrintChi2Table(setD2,printDecomposedChi2=False)

#%%
#######################################
# Minimisation
#######################################
import time

def chi2(x):
    startT=time.time()
    #harpy.setNPparameters_uTMDFF([x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7]])
    SnowFlake.setNPparameters(x)
    SnowFlake.UpdateTables(1.0, 8.0)
    print('np set =',["{:8.3f}".format(i) for i in x])        
            
    YY=DataProcessor.snowInterface.ComputeXSec(setD2)
    ccD2,cc3=setD2.chi2(YY)    
    
    endT=time.time()
    print(':->',ccD2/setD2.numberOfPoints,"    time=",endT-startT)
    return ccD2

#%%
#### Minimize SIDIS
from iminuit import Minuit

#---- PDFbias-like row (0.083931)
initialValues=(1.0,1.0, 
                0.1,0.,0.,
                0.1,0.,0.,
                -0.1,0.,0.,
                -0.1,0.,0.,
                0.,0.,
                1.,0.)

initialErrors=(0.1,0.1, 
                0.1,0.1,0.1,
                1.,0.1,0.1,
                1.,0.1,0.1,
                1.,0.1,0.1,
                1.,0.1,1.,0.1)
searchLimits=((1.,5.),(1.,5.),
              (-50.,50.), (-50.,50.), (-50.,50.),
              (-50.,50.), (-50.,50.), (-50.,50.),
              (-50.,50.), (-50.,50.), (-50.,50.),
              (-50.,50.), (-50.,50.), (-50.,50.),
              (-50.,50.), (-50.,50.), 
              (-50.,50.), (-50.,50.))
              
# True= FIX
parametersToMinimize=(True, True,
                      False, True, True,
                      False, True, True,
                      False, True, True,
                      False, True, True,
                      True, True, 
                      False,True)

#%%

m = Minuit(chi2, initialValues)

m.errors=initialErrors
m.limits=searchLimits
m.fixed=parametersToMinimize
m.errordef=1

print(m.params)
#%%
#m.tol=0.0001*(setSIDIS.numberOfPoints+setDY.numberOfPoints)*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1
m.migrad()

print(m.params)

chi2(list(m.values))

DataProcessor.harpyInterface.PrintChi2Table(setD2,printDecomposedChi2=True)

print([round(x,1 if x >100 else 4 if x>1 else 6) for x in list(m.values)])
