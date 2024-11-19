#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 16:40:59 2019

@author: vla18041
"""
#######################################
# importing libraries
#######################################
import os
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"

ATMDE_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide/"
HARPY_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide/harpy/"


import sys
import numpy
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)

#%%
import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet

MAINPATH=ROOT_DIR 
#%%
#######################################
#Initialize artemide
#######################################
import harpy

path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/DYpi_N4LL_resum.atmde"


harpy.initialize(path_to_constants)

inARRAY_TMDR=[1.5, 0.071624, 0.064726, 0.0]
inARRAY_PDF=[0.521462, 0.000206, 0.402948, 7.0219, 1.0, 20.4051, 1.0, 0.000123, 1.1037, 0.660734, 0.0, 0.04]


harpy.setNPparameters_TMDR(inARRAY_TMDR)
harpy.setNPparameters_uTMDPDF(inARRAY_PDF)

#%%
### read the list of files and return the list of DataSets
def loadThisDataDY(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=ROOT_DIR+"DataLib/unpolDY/"
    path_to_dataW=ROOT_DIR+"DataLib/unpolW/"
    path_to_dataPI=ROOT_DIR+"DataLib/piDY/"
    
    
    dataCollection=[]
    for name in listOfNames:
        if("W" in name):          
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataW+name+".csv")
        elif("pi" in name):
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataPI+name+".csv")
        else:
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        

    return dataCollection

#%%
##################Cut function
def cutFunc(p):

    par=0.5
    
    #  for artemide v3.    
    # p["process"]=[p["process"][0],p["process"][2],1,1]
    if(len(p["process"])==3):        
            print("UNKNOWN PROCESS IN ARTEMIDE 3"+str(p["process"]))
    
    if(p["xSec"]>0):
        err=numpy.sqrt(sum([i**2 for i in p["uncorrErr"]]))/p["xSec"]
    else:
        err=100.
    delta=p["<qT>"]/p["<Q>"]
    
    if(p["id"][0] == "E"):
        delta=p["<qT>"]/p["Q"][1] 
    
    if("run1-W" in p["id"]):
        delta=p["qT"][0]/(p["Q"][0]+5.)
    
    
    if(p["id"][0:4] == "E605"):
        if(p["Q"][0]==10.5):#UPSILON resonance-bin
            return False , p
    elif(p["id"][0:4] == "E772"):
        if(p["Q"][0]<10):#these bins seems broken
            return False , p
    elif(p["id"][0:4] == "E615"):
        if(9<p["<Q>"]<11.2):#UPSILON resonance-bin
            return False , p
    elif(p["id"][0:4] == "E228"):
        if(9<p["<Q>"]<11):#UPSILON resonance-bin
            return False , p
    else:
        if(9<p["<Q>"]<11):#UPSILON resonance-bin
            return False , p    
    
    return ((delta<0.25 and p["<qT>"]<10.) or (delta<0.25 and par/err*delta**2<1)) , p

#%%
##################Cut function
def cutFuncPI(p):
    
    delta=p["<qT>"]/p["<Q>"]
    
    if(p["id"][0] == "E"):
        delta=p["<qT>"]/p["Q"][1] 
    
    return (delta<0.1 or delta<0.3) , p

#%%
### Loading the DY data set
### This collection of dta is usual Drel-Yan process it is used for cross-check
theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY([
                          'CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                          #'A7-00y10', 'A7-10y20','A7-20y24', 
                          'A8-00y04', 'A8-04y08', 'A8-08y12', 'A8-12y16', 'A8-16y20', 'A8-20y24', 
                          'A8-46Q66', 'A8-116Q150', 
                          'A13-norm',
                          'CMS7', 'CMS8', 
                          'CMS13-00y04','CMS13-04y08','CMS13-08y12','CMS13-12y16','CMS13-16y24',
                          #'CMS13_dQ_50to76',
                          'CMS13_dQ_106to170','CMS13_dQ_170to350','CMS13_dQ_350to1000',
                          'LHCb7', 'LHCb8', 'LHCb13_dy(2021)', 
                          'PHE200', 'STAR510', 
                          'E228-200', 'E228-300', 'E228-400', 
                          'E772',
                          'E605',
                          'D0run1-W','CDFrun1-W'
                          ]))

setDY=theData.CutData(cutFunc) 

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')

#%%
### Loading the pionDY data set
### This collection of dta is usual Drel-Yan process it is used for cross-check
### Loading the data set
theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY([
                                'E615(pi)-dxF-0.0','E615(pi)-dxF-0.1','E615(pi)-dxF-0.2','E615(pi)-dxF-0.3',
                                'E615(pi)-dxF-0.4','E615(pi)-dxF-0.5','E615(pi)-dxF-0.6','E615(pi)-dxF-0.7'
                                ]))

setDY_PI=theData.CutData(cutFuncPI) 

print('Loaded ', setDY_PI.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY_PI.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setDY_PI.sets])

#%%
### this should compute the chi2 for the usual DY case
### and show some table of values.
harpy.setNPparameters([1.5, 0.071624, 0.064726, 0.0,
                       0.568631, 1.1114, 0.563709, 7.2062, 
                       4.5436, 17.3135, 0.976618, 0.006934, 
                       1.3213, 25.9727, 1.0, 1.0])

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)

#%%

### TASK 0: run this code and get the table
### TASK 1: make a plot of any cross-section (theory ad the data points)
###       some help : the points for cross-section can be taken from the liabrary here
###         e.g. setDY_PI.set[0].point[0]["xSec"]  
###             will give you the value of cross-section for the first point in given set
###      the theory prediction for the cross-section can be obtain by
###         DataProcessor.harpyInterface.ComputeXSec(setDY_PI.set[0])
###         it produces the list of predictions point-by-point for the given experiment

