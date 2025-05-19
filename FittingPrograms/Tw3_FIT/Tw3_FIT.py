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
dd=SnowFlake.G2List([0.1,0.2,0.3,0.4],[5.,5.,5.,5.],[100,100,100,100])
print(dd)

#%%
dd=SnowFlake.D2List([5.,5.,5.,5.],[100,100,100,100])
print(dd)