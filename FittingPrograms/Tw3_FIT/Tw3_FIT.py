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
### read the list of files and return the list of DataSets
def loadThisDataG2(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=ROOT_DIR+"DataLib/G2/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        
    return dataCollection

#%%
##################Cut function
def cutFunc(p):
    if p["type"]=="G2":
        if p["<Q>"]<numpy.sqrt(2.):
            return False, p
    
    return True, p

#%%
### Loading the D2 data set
theData=DataProcessor.DataMultiSet.DataMultiSet("D2set",loadThisDataD2([
    "E143_d2","E154_d2","E155-1999_d2","E155_d2",
    "HallA-2016_d2","HERMES_d2","SANE_d2",
    "RQCD_d2_ud"
    ]))

setD2=theData.CutData(cutFunc) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setD2.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setD2.sets]), 'points.') 

#%%
### Loading the D2 data set
theData=DataProcessor.DataMultiSet.DataMultiSet("G2set",loadThisDataG2([
    #"E142.n", "E143.p", "E143.d","E143.n", 
    "E154.n",
    "E155-29.p","E155-32.p","E155-38.p",
    #"E155-29.d","E155-32.d","E155-38.d",
    #"SMC.p",
    #"HERMES",
    "HallA-2004.n","HallA-2016-4.He3","HallA-2016-5.He3",
    ]))

setG2=theData.CutData(cutFunc) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setG2.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setG2.sets]), 'points.') 

#%%
SnowFlake.setNPparameters([5.0,1.0, 
                0.05,0.,0.0,0.,-0.06,0.,0.0,0.,0.0,0.,
                0.0,0.,0.0,0.,0.0,0.])
SnowFlake.UpdateTables(1.0, 8.0)

#%%
DataProcessor.snowInterface.PrintChi2Table(setD2,printDecomposedChi2=False)
DataProcessor.snowInterface.PrintChi2Table(setG2,printDecomposedChi2=False)

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
    
    YY=DataProcessor.snowInterface.ComputeXSec(setG2)
    ccG2,cc3=setD2.chi2(YY)    
    
    chiTOTAL=(ccD2/setD2.numberOfPoints+ccG2/setG2.numberOfPoints)*(setD2.numberOfPoints+setG2.numberOfPoints)
    
    endT=time.time()
    print(':->',ccD2/setD2.numberOfPoints," ",ccG2/setG2.numberOfPoints,"    time=",endT-startT)
    return chiTOTAL

#%%
from iminuit import Minuit

#---- PDFbias-like row (0.083931)
initialValues=(1.0,1.0, 
                0.1,0.,
                0.0,0.,
                0.1,0.,
                0.0,0.,
                0.0,0.,
                0.0,0.,0.0,0.,
                0.0,0.)

initialErrors=(0.1,0.1, 
                0.1,0.1,
                0.1,0.1,
                0.1,0.1,
                0.1,0.1,
                0.1,0.1,
                0.1,0.1,0.1,0.1,0.1,0.1)
searchLimits=((1.,10.),(1.,10.),
              (-50.,50.), (-50.,50.), (-50.,50.),
              (-50.,50.), (-50.,50.), (-50.,50.),
              (-50.,50.), (-50.,50.), (-50.,50.),
              (-50.,50.), (-50.,50.), (-50.,50.),
              (-50.,50.), (-50.,50.), 
              (-50.,50.), (-50.,50.))
              
# True= FIX
parametersToMinimize=(False, True,
                      False, True,
                      True,True,
                      False, True,
                      True,True,
                      True,True,
                      True,True,True,True,
                      True,True)

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

DataProcessor.snowInterface.PrintChi2Table(setD2,printDecomposedChi2=True)
DataProcessor.snowInterface.PrintChi2Table(setG2,printDecomposedChi2=True)

print([round(x,1 if x >100 else 4 if x>1 else 6) for x in list(m.values)])
