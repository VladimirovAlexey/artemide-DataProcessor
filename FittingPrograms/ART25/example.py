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

path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/ART25_main.atmde"


harpy.initialize(path_to_constants)

inARRAY_TMDR=[1.5, 0.071624, 0.064726, 0.0]
inARRAY_PDF=[0.521462, 0.000206, 0.402948, 7.0219, 1.0, 20.4051, 1.0, 0.000123, 1.1037, 0.660734, 0.0, 0.04]
inARRAY_FF=[0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0,0.0,0.0,0.0,0.0]


harpy.setNPparameters_TMDR(inARRAY_TMDR)
harpy.setNPparameters_uTMDPDF(inARRAY_PDF)
harpy.setNPparameters_uTMDFF(inARRAY_FF)

#%%
### read the list of files and return the list of DataSets
def loadThisDataDY(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=ROOT_DIR+"DataLib/unpolDY/"
    path_to_dataW=ROOT_DIR+"DataLib/unpolW/"
    path_to_dataA=ROOT_DIR+"DataLib/DY_angular/"
    
    
    dataCollection=[]
    for name in listOfNames:
        if(name[-1]=="W"):
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataW+name+".csv")
        elif("_A4" in name):
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataA+name+".csv")
        elif("_Auu" in name):
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataA+name+".csv")
        else:
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        

    return dataCollection

#%%
##################Cut function
def cutFunc(p):
    
    if p["type"]=="DY":
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
        
        return p["<qT>"]<25., p
        #return ((delta<0.25 and p["<qT>"]<25.) or (delta<0.25 and par/err*delta**2<1)) , p
    
   
#%%
### Loading the DY data set
# theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY([
#                           'CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
#                           #'A7-00y10', 'A7-10y20','A7-20y24', 
#                           'A8-00y04', 'A8-04y08', 'A8-08y12', 'A8-12y16', 'A8-16y20', 'A8-20y24', 
#                           'A8-46Q66', 'A8-116Q150', 
#                           'A13-norm',
#                           'CMS7', 'CMS8', 
#                           'CMS13-00y04','CMS13-04y08','CMS13-08y12','CMS13-12y16','CMS13-16y24',
#                           #'CMS13_dQ_50to76',
#                           'CMS13_dQ_106to170','CMS13_dQ_170to350','CMS13_dQ_350to1000',
#                           'LHCb7', 'LHCb8', 'LHCb13_dy(2021)', 
#                           'PHE200', 'STAR510', 
#                           'E228-200', 'E228-300', 'E228-400', 
#                           'E772',
#                           'E605',
#                           'D0run1-W','CDFrun1-W'
#                           ]))

theDataEX=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY([                          
                          'A8-00y04'
                          ]))

#setDY=theData.CutData(cutFunc) 
setDYEX=theDataEX.CutData(cutFunc) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

#print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')

#%%
for i in range(len(setDYEX.points)):
    setDYEX.points[i]["includeCuts"]=False

#%%

####### Best fast result
harpy.setNPparameters([1.5, 0.083931, 0.030641, 0.0, 
                       0.51638, 0.002073, 0.478567, 0.373111, 
                       2.407, 22.1996, 3.7876, 0.00128, 
                       0.403343, 5e-05, 1.0, 1.0, 
                       0.69769, 0.712969, -0.133895, -0.841651, 0.846846,
                       0.774759, 1.5565, 1.1863, 0.692877, -0.569062, 
                       0.0, 0.0])
#%%

X=DataProcessor.harpyInterface.ComputeXSec(setDYEX)