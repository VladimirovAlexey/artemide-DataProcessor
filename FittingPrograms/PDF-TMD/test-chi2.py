#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 14:35:35 2021

@author: vla18041
"""

##############################
# Ploting original SV19 fit
##############################

#######################################
# importing libraries
#######################################
import os
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"

ATMDE_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide/"
#%%
import sys

import time
import numpy
sys.path.append(ROOT_DIR)
import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet

MAINPATH=ROOT_DIR
#%%
#######################################
#Initialize artemide with a replica -2
#######################################
import harpy
path_to_constants=MAINPATH+"FittingPrograms/SV19/Constants-files/"
harpy.initialize(path_to_constants+"DY_n3lo/const-NNPDF31_n3lo")
#harpy.initialize(path_to_constants+"DY_nnlo/const-HERA20_NNLO")
#harpy.initialize(path_to_constants+"DY_n3lo/const-HERA20_n3lo")
#harpy.initialize(path_to_constants+"DY_nnlo/const-MMHT14_NNLO")
#harpy.initialize(path_to_constants+"DY_nnlo/const-CT14_NNLO")
#harpy.initialize(path_to_constants+"DY_nnlo/const-PDF4LHC_NNLO")
harpy.setNPparameters([2.0340, 0.0299, 0.2512, 7.7572,334.6108, 2.4543,-4.8203, 0.1000,  0.0000])
#harpy.setNPparameters_TMDR(-2)
#harpy.setNPparameters_uTMDPDF(-2)

#%%

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(ATMDE_DIR+"Models/SV19/Replicas/"+
                                                  "DY_n3lo/DY_NNPDF31_n3lo.rep")
                                                  # "Sivers20_model9case1(noDY-n3lo).rep")

rSet.SetReplica()

#%%
#######################################
# read the list of files and return the list of DataSets
#######################################
def loadThisData(listOfNames):    
    import DataProcessor.DataSet
    
    dataCollection=[]
    for name in listOfNames:
        if( name==''): continue
        loadedData=DataProcessor.DataSet.LoadCSV(ROOT_DIR+"DataLib/unpolDY/"+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

#%%
#######################################
# Data cut function
#######################################
#### Check the point kinematics
def cutFunc(p):    
    
    par=1.0

    if(p["xSec"]>0):
        err=numpy.sqrt(sum([i**2 for i in p["uncorrErr"]]))/p["xSec"]
        #err=10000.0 ## I would like to add all points into the full set
    else:
        err=100.
    delta=p["<qT>"]/p["<Q>"]
    
    if(p["id"][0] == "E"):
        delta=p["<qT>"]/p["Q"][1]    
        
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
        
    if(p["id"][-2:]=="<u" and p["<Q>"]>10.5):
        return False,p
    if(p["id"][-2:]==">u" and p["<Q>"]<10.5):
        return False,p
    
#    return delta<0.5 and p.qT_avarage<80
    return (delta<0.1 or (delta<0.25 and par/err*delta**2<1)) , p
#%%
#######################################
# Loading the data set
#######################################

setHE=loadThisData(['CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                      'A7-00y10','A7-10y20','A7-20y24',
                      'A8-00y04','A8-04y08','A8-08y12',
                      'A8-12y16','A8-16y20','A8-20y24',
                      'A8-46Q66','A8-116Q150',
                      'CMS7', 'CMS8', 
                      'LHCb7', 'LHCb8', 'LHCb13'])
setLE=loadThisData(['PHE200', 'E228-200', 'E228-300', 'E228-400','E772','E605'])

theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",setHE+setLE)

setDY=theData.CutData(cutFunc) 

print('Loaded ', setDY.numberOfSets, 'data sets with ', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setDY.sets])

#%%

setLHCb=loadThisData(['LHCb13(2021)','LHCb13_dy(2021)'])
theDataLHCb=DataProcessor.DataMultiSet.DataMultiSet("DYset",setLHCb)

setLHCb=theDataLHCb.CutData(cutFunc) 

print('Loaded ', setLHCb.numberOfSets, 'data sets with ', sum([i.numberOfPoints for i in setLHCb.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setLHCb.sets])
                      

#%%
rSet.SetReplica()
#harpy.setNPparameters([1.93, 0.0434,0.195, 9.117, 444., 2.12, -4.89,0.,0.,0.258, 0.478, 0.484, 0.459])
#DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)
DataProcessor.harpyInterface.PrintChi2Table(setLHCb,printDecomposedChi2=True)

#%%
qq=DataProcessor.harpyInterface.ComputeXSec(setLHCb.sets[1])