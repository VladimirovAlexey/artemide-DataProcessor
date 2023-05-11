#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 15:47:01 2023

@author: alexey
"""

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

replicaFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/replicaDY.dat"
logFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/LOG.log"

import sys
import numpy
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)

### number of circles that this code will run
NumberOfReplicas=25

#%%
import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet
import DataProcessor.DataMultiSet

MAINPATH=ROOT_DIR 
#%%
#######################################
#Initialize artemide
#######################################
import harpy


path_to_constants=MAINPATH+"FittingPrograms/MVZ22/ConstantsFiles/DYonly/ForScale_N2LL"
#path_to_constants=MAINPATH+"FittingPrograms/MVZ22/ConstantsFiles/DYonly/ForJAM_N2LL"

TXT="N4LO"

harpy.initialize(path_to_constants)

initializationArray=[0.253434, 0.253434,0.253434, 0.253434,
                     0.253434, 0.253434,0.253434, 0.253434,
                     0.253434, 0.253434, 0.1,  0.1]

#################################
### Note that to vary the c3 scale I introduced the c3 as 4th NP parameter into RAD
### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#################################

harpy.setNPparameters([1.4806, 0.038969, 0.051737, 1.0, 
                       0.851645, 0.69432, 0.934676, 5.2514, 
                       0.247602, 39.6123, 0.094435, 1.8872, 
                       1.2164, 0.936465, 0.0, 0.0])

#%%
### read the list of files and return the list of DataSets
def loadThisDataDY(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=ROOT_DIR+"DataLib/unpolDY/"
    path_to_dataW=ROOT_DIR+"DataLib/unpolW/"
    
    
    dataCollection=[]
    for name in listOfNames:
        if(name[-1]=="W"):
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataW+name+".csv")
        else:
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        

    return dataCollection
#%%
##################Cut function
def cutFunc(p):
    par=0.5
    
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
    
#    return delta<0.5 and p.qT_avarage<80
    return ((delta<0.25 and p["<qT>"]<10.) or (delta<0.25 and par/err*delta**2<1)) , p

#%%
### Loading the DY data set
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
#setDYfit=theData.CutData(cutFuncFORFIT) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
#print('Loaded ', setDYfit.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDYfit.sets]), 'points.')


#%%
#(11+full data)
harpy.setNPparameters([1.4806, 0.038969, 0.051737, 1.0, 
                       0.851645, 0.69432, 0.934676, 5.2514, 
                       0.247602, 39.6123, 0.094435, 1.8872, 
                       1.2164, 0.936465, 0.0, 0.0])
DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)

#%%
import copy

#### alike A13 but with 0.5 qT-bin

DataTestA13=DataProcessor.DataSet.DataSet('A13-test',"DY")
DataTestA13.comment="test"
DataTestA13.reference="..."

refpoint=copy.deepcopy(setDY.sets[13].points[0])

for i in range(25):
    refpoint["<qT>"]=1.*(i+0.5)
    refpoint["qT"]=[1.*i,1.*(i+1)]
    refpoint["includeCuts"]=False
    DataTestA13.AddPoint(refpoint)
DataTestA13.FinalizeSet()

#### alike PHENIX but with 0.02 qT-bin

DataTestPHE=DataProcessor.DataSet.DataSet('PHE-test',"DY")
DataTestPHE.comment="test"
DataTestPHE.reference="..."

refpoint=copy.deepcopy(setDY.sets[27].points[0])

for i in range(50):
    refpoint["<qT>"]=0.04*(i+0.5)
    refpoint["qT"]=[0.04*i,0.04*(i+1)]
    DataTestPHE.AddPoint(refpoint)
DataTestPHE.FinalizeSet()

#%%
dd=DataTestA13

result={}

up=2.
down=0.75 ### because at 0.5, the mu=0.7 GeV (too small)

print("111")
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, 1.0])
harpy.varyScales(1.,1.,1.,1.)
result["111"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("+11")
harpy.varyScales(1.,2.,1.,1.)
result["+11"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("-11")
harpy.varyScales(1.,0.5,1.,1.)
result["-11"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("1+1")
harpy.varyScales(1.,1.,1.,1.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, up])
result["1+1"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("1-1")
harpy.varyScales(1.,1.,1.,1.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, down])
result["1-1"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("++1")
harpy.varyScales(1.,2.,1.,1.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, up])
result["++1"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("+-1")
harpy.varyScales(1.,2.,1.,1.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, down])
result["+-1"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("-+1")
harpy.varyScales(1.,0.5,1.,1.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, up])
result["-+1"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("--1")
harpy.varyScales(1.,0.5,1.,1.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, down])
result["--1"]=DataProcessor.harpyInterface.ComputeXSec(dd)

########################################
print("11+")
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, 1.0])
harpy.varyScales(1.,1.,1.,2.)
result["11+"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("+1+")
harpy.varyScales(1.,2.,1.,2.)
result["+1+"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("-1+")
harpy.varyScales(1.,0.5,1.,2.)
result["-1+"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("1++")
harpy.varyScales(1.,1.,1.,2.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, up])
result["1++"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("1-+")
harpy.varyScales(1.,1.,1.,2.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, down])
result["1-+"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("+++")
harpy.varyScales(1.,2.,1.,2.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, up])
result["+++"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("+-+")
harpy.varyScales(1.,2.,1.,2.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, down])
result["+-+"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("-++")
harpy.varyScales(1.,0.5,1.,2.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, up])
result["-++"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("--+")
harpy.varyScales(1.,0.5,1.,2.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, down])
result["--+"]=DataProcessor.harpyInterface.ComputeXSec(dd)

########################################
print("11-")
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, 1.0])
harpy.varyScales(1.,1.,1.,.5)
result["11-"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("+1-")
harpy.varyScales(1.,2.,1.,.5)
result["+1-"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("-1-")
harpy.varyScales(1.,0.5,1.,.5)
result["-1-"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("1+-")
harpy.varyScales(1.,1.,1.,.5)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, up])
result["1+-"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("1--")
harpy.varyScales(1.,1.,1.,.5)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, down])
result["1--"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("++-")
harpy.varyScales(1.,2.,1.,.5)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, up])
result["++-"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("+--")
harpy.varyScales(1.,2.,1.,.5)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, down])
result["+--"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("-+-")
harpy.varyScales(1.,0.5,1.,.5)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, up])
result["-+-"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("---")
harpy.varyScales(1.,0.5,1.,.5)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, down])
result["---"]=DataProcessor.harpyInterface.ComputeXSec(dd)


DUMPPATH="/data/WorkingFiles/TMD/Fit_Notes/ART23/ScaleVariations/"

from json import dump
with open(DUMPPATH+"A13-"+TXT, 'w') as fp:
    dump(result, fp)

#%%

dd=DataTestPHE

result={}

print("111")
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, 1.0])
harpy.varyScales(1.,1.,1.,1.)
result["111"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("+11")
harpy.varyScales(1.,2.,1.,1.)
result["+11"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("-11")
harpy.varyScales(1.,0.5,1.,1.)
result["-11"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("1+1")
harpy.varyScales(1.,1.,1.,1.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, up])
result["1+1"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("1-1")
harpy.varyScales(1.,1.,1.,1.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, down])
result["1-1"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("++1")
harpy.varyScales(1.,2.,1.,1.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, up])
result["++1"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("+-1")
harpy.varyScales(1.,2.,1.,1.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, down])
result["+-1"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("-+1")
harpy.varyScales(1.,0.5,1.,1.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, up])
result["-+1"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("--1")
harpy.varyScales(1.,0.5,1.,1.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, down])
result["--1"]=DataProcessor.harpyInterface.ComputeXSec(dd)

########################################
print("11+")
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, 1.0])
harpy.varyScales(1.,1.,1.,2.)
result["11+"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("+1+")
harpy.varyScales(1.,2.,1.,2.)
result["+1+"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("-1+")
harpy.varyScales(1.,0.5,1.,2.)
result["-1+"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("1++")
harpy.varyScales(1.,1.,1.,2.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, up])
result["1++"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("1-+")
harpy.varyScales(1.,1.,1.,2.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, down])
result["1-+"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("+++")
harpy.varyScales(1.,2.,1.,2.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, up])
result["+++"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("+-+")
harpy.varyScales(1.,2.,1.,2.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, down])
result["+-+"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("-++")
harpy.varyScales(1.,0.5,1.,2.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, up])
result["-++"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("--+")
harpy.varyScales(1.,0.5,1.,2.)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, down])
result["--+"]=DataProcessor.harpyInterface.ComputeXSec(dd)

########################################
print("11-")
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, 1.0])
harpy.varyScales(1.,1.,1.,.5)
result["11-"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("+1-")
harpy.varyScales(1.,2.,1.,.5)
result["+1-"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("-1-")
harpy.varyScales(1.,0.5,1.,.5)
result["-1-"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("1+-")
harpy.varyScales(1.,1.,1.,.5)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, up])
result["1+-"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("1--")
harpy.varyScales(1.,1.,1.,.5)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, down])
result["1--"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("++-")
harpy.varyScales(1.,2.,1.,.5)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, up])
result["++-"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("+--")
harpy.varyScales(1.,2.,1.,.5)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, down])
result["+--"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("-+-")
harpy.varyScales(1.,0.5,1.,.5)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, up])
result["-+-"]=DataProcessor.harpyInterface.ComputeXSec(dd)

print("---")
harpy.varyScales(1.,0.5,1.,.5)
harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, down])
result["---"]=DataProcessor.harpyInterface.ComputeXSec(dd)


DUMPPATH="/data/WorkingFiles/TMD/Fit_Notes/ART23/ScaleVariations/"

from json import dump
with open(DUMPPATH+"PHE-"+TXT, 'w') as fp:
    dump(result, fp)