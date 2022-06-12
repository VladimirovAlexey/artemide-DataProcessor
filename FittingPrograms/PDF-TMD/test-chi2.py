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
DATAPROC_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/"

#PDFinUse="SV19"
#PDFinUse="HERA20"
PDFinUse="NNPDF31"
#PDFinUse="CT18"
#PDFinUse="MSHT20"

if(PDFinUse=="SV19"):
    ATMDE_DIR = ROOT_DIR+"artemide/"
else:
    ATMDE_DIR = ROOT_DIR+"artemide-PDF/"

#%%
import sys

import numpy
sys.path.append(DATAPROC_DIR)
sys.path.remove(os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+'/artemide/harpy')
sys.path.append(ATMDE_DIR+'harpy')

import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet



#%%
#######################################
#Initialize artemide with a replica -2
#######################################
import harpy

if(PDFinUse=="SV19"):
    harpy.initialize(DATAPROC_DIR+"FittingPrograms/SV19/Constants-files/DY_n3lo/const-NNPDF31_n3lo")
    rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(ATMDE_DIR+"Models/SV19/Replicas/"+
                                                  "DY_n3lo/DY_NNPDF31_n3lo.rep")
    
else:
    harpy.initialize(ROOT_DIR
        +"artemide/Models/PDFbias22/Constants-files/const-"+PDFinUse+"_NNLO_12p")
    rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(ROOT_DIR
        +"artemide/Models/PDFbias22/REPS/SV21-"+PDFinUse+"-nnlo-EXP.rep")

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
        loadedData=DataProcessor.DataSet.LoadCSV(DATAPROC_DIR+"DataLib/unpolDY/"+name+".csv")
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

ss=setDY
print('Loaded ', ss.numberOfSets, 'data sets with ', sum([i.numberOfPoints for i in ss.sets]), 'points.')
print('Loaded experiments are', [i.name for i in ss.sets])

#%%
DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)
#%%

setLHCb=loadThisData(['LHCb13(2021)','LHCb13_dy(2021)'])
theDataLHCb=DataProcessor.DataMultiSet.DataMultiSet("DYset",setLHCb)

setLHCb=theDataLHCb.CutData(cutFunc) 

ss=setLHCb
print('Loaded ', ss.numberOfSets, 'data sets with ', sum([i.numberOfPoints for i in ss.sets]), 'points.')
print('Loaded experiments are', [i.name for i in ss.sets])
                      
#%%

setCMS13=loadThisData(['CMS13_dQ_50to76','CMS13_dQ_76to106',
                       'CMS13_dQ_106to170','CMS13_dQ_170to350',
                       'CMS13_dQ_350to1000'])
theDataCMS13=DataProcessor.DataMultiSet.DataMultiSet("DYset",setCMS13)

setCMS13=theDataCMS13.CutData(cutFunc) 

ss=setCMS13
print('Loaded ', ss.numberOfSets, 'data sets with ', sum([i.numberOfPoints for i in ss.sets]), 'points.')
print('Loaded experiments are', [i.name for i in ss.sets])

#%%
rSet.SetReplica()
DataProcessor.harpyInterface.PrintChi2Table(setLHCb,printDecomposedChi2=True)

#%%
rSet.SetReplica()
DataProcessor.harpyInterface.PrintChi2Table(setCMS13,printDecomposedChi2=True)

#%%
qq=DataProcessor.harpyInterface.ComputeXSec(setLHCb.sets[1])