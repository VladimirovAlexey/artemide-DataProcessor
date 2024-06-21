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

replicaFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/replicaCHI.dat"
replicaFileA0 =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/plotA0.dat"
replicaFileA1 =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/plotA1.dat"
replicaFileA2 =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/plotA2.dat"

import sys
import numpy
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)

### number of circles that this code will run
NumberOfReplicas=100


#%%
import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet


MAINPATH=ROOT_DIR 
#%%
#######################################
#Initialize artemide
#######################################
import harpy

path_to_constants=MAINPATH+"FittingPrograms/KPC/INI/ART23_MSHT_N4LL+BM.atmde"
#path_to_constants=MAINPATH+"FittingPrograms/ART23/ConstantsFiles/DYonly/ART23_JAM_NLL.atmde"


harpy.initialize(path_to_constants)

initializationArray=[0.253434, 0.253434,0.253434, 0.253434,
                     0.253434, 0.253434,0.253434, 0.253434,
                     0.253434, 0.253434, 0.1,  0.1]

harpy.setNPparameters_TMDR([1.584237, 0.048428,0.001,0.])

harpy.setNPparameters_uTMDPDF(initializationArray)
harpy.setNPparameters_BM_TMDPDF([0.1,0.1,0.1,0.1])

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
### read the list of files and return the list of DataSets
def loadThisDataDY_angular(listOfNames):    
    import DataProcessor.DataSet
    path_to_dataA=ROOT_DIR+"DataLib/DY_angular/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataA+name+".csv")
        dataCollection.append(loadedData)           

    return dataCollection

#%%
##################Cut function
def cutFunc_angular(p):
    if(len(p["process"])==4):
        if(p["process"][3]==3): p["process"][3]=2
    else:
        print("UNKNOWN PROCESS IN ARTEMIDE 3"+str(p["process"]))
    
    return (p["<qT>"]<10.) , p

def cutFunc_angularBM(p):
    if(len(p["process"])==4):
        if(p["process"][3]==3): p["process"][3]=2
    else:
        print("UNKNOWN PROCESS IN ARTEMIDE 3"+str(p["process"]))
    
    p["process"][3]=p["process"][3]+10
    
    return (p["<qT>"]<10.) , p

def cutFunc_angularBMPLOT(p):
    if(len(p["process"])==4):
        if(p["process"][3]==3): p["process"][3]=2
    else:
        print("UNKNOWN PROCESS IN ARTEMIDE 3"+str(p["process"]))
    
    p["process"][3]=p["process"][3]+10
    
    delta=p["<qT>"]/p["<Q>"]
    return (delta<0.25) , p

#%%
### Load the data for A0
theDataA0=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A0_0y1","A8_A0_1y2","A8_A0_2y35"]))
setA0=theDataA0.CutData(cutFunc_angular)
setA0_BM=theDataA0.CutData(cutFunc_angularBM)

### Load the data for A1
theDataA1=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A1_0y1","A8_A1_1y2"]))
setA1=theDataA1.CutData(cutFunc_angular)
setA1_BM=theDataA1.CutData(cutFunc_angularBM)
### Load the data for A2
theDataA2=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A2_0y1","A8_A2_1y2","A8_A2_2y35"]))
setA2=theDataA2.CutData(cutFunc_angular)
setA2_BM=theDataA2.CutData(cutFunc_angularBM)

#%%
setA0_BM_PLOT=theDataA0.CutData(cutFunc_angularBMPLOT)
setA1_BM_PLOT=theDataA1.CutData(cutFunc_angularBMPLOT)
setA2_BM_PLOT=theDataA2.CutData(cutFunc_angularBMPLOT)
#%%
# MAIN FIT
harpy.setNPparameters([1.5004, 0.05614, 0.03862, 0.0, 0.565, 0.0539, 0.5697, 6.64, 0.565, 20.07, 0.5697, 0.537, 1.07, 2.39, 0.0, 0.0, 0.5, 1.,1.,1.])

#%%
#### Computation of different contributions due to the unpolarized part
# XXuu=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setAuu))
# XX0=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA0))
# XX1=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA1))
# XX2=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA2))

# print("---- NORMING EXPRESSIONS---------")
# print("--------Xuu---------")
# print(XXuu)
# print("--------XX0---------")
# print(XX0)
# print("--------XX1---------")
# print(XX1)
# print("--------XX2---------")
# print(XX2)
# print("--------------------")

XXuu=numpy.array([12.78170039, 23.52787371, 25.3692582,  22.3286602,  12.59828901, 22.89768367,
 24.27085174, 21.01501173, 15.55497796, 27.35012268, 28.04748583, 23.46633026])
XX2=numpy.array([0.00133855, 0.01438535, 0.05392202, 0.11983518, 0.00139393, 0.01477209,
 0.0544654,  0.11931541, 0.00192452, 0.01966457, 0.07007986, 0.14949457])


#%%
#######################################
# Minimisation
#######################################
import time

def chi_BM(x):
    startT=time.time()
    harpy.setNPparameters_BM_TMDPDF([x[0],x[1],x[2],x[3]])
    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")        
    
    Y2=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA2_BM,method="approximate"))
    ccA2,cc3=setA2.chi2((Y2+XX2)/XXuu)
    
    chiTotal=ccA2
    cc=chiTotal/setA2.numberOfPoints
    
    endT=time.time()
    print(':->',cc,'       t=',endT-startT)
    return chiTotal

#%%
#######################################
# Generate replica of data and compute chi2
#######################################
from iminuit import Minuit

def MinForReplica():
    global setA2
        
    def repchi_2(x): 
        startT=time.time()
        harpy.setNPparameters_BM_TMDPDF([x[0],x[1],x[2],x[3]])
        print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
        
        Y2=numpy.array(DataProcessor.harpyInterface.ComputeXSec(repDataDY,method="approximate"))
        ccA2,cc3=repDataDY.chi2((Y2+XX2)/XXuu)

        cc=ccA2/repDataDY.numberOfPoints
        
        endT=time.time()
        print(':->',cc,'       t=',endT-startT)
        return ccA2
    
    repDataDY=setA2_BM.GenerateReplica()
    
    initialValues=[0.01,numpy.random.uniform(-3.,0.5),numpy.random.uniform(-0.8,1),1.]   
    
    localM = Minuit(repchi_2, initialValues)      
    
    localM.errors=(0.3,0.3,0.3,0.1)
    localM.limits=((0.01,10.),(-100,100.),(-0.95,10.),(0.3,20.))
    localM.fixed=(True, False,False,True)
    localM.errordef=1        
    localM.strategy=1

    localM.migrad()
    
    ### [chi^2, NP-parameters]
    return [localM.fval,list(localM.values)]

#%%
#######################################
# This is the main cicle. 
# It generates replica of data take random PDF and minimize it
# Save to log.
#######################################

for i in range(NumberOfReplicas):
    print('---------------------------------------------------------------')
    print('------------REPLICA ',i,' from ',NumberOfReplicas,'------------------')
    print('---------------------------------------------------------------')
    savedTime=time.time()
    
    print("Minimization started.")
    
    ## got to pseudo-data and minimization
    repRes=MinForReplica()
    print(repRes)
    print("Minimization finished.")    
    
    ## compute the chi2 for true data full
    mainChi=chi_BM(repRes[1])
    
    print("Central chi^2  computed.)")
    
    ## save to file
    f=open(replicaFile,"a+")
    print('SAVING >>  ',f.name)
    ### [total chi^2, list of NP-parameters],
    f.write(str([mainChi,repRes[1]])+"\n")
    f.close()  
    
    XX0_BM=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA0_BM_PLOT))
    XX1_BM=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA1_BM_PLOT))
    XX2_BM=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA2_BM_PLOT))
    
    f=open(replicaFileA0,"a+")
    print('SAVING PLOTS A0>>  ',f.name)
    ### [total chi^2, list of NP-parameters],
    f.write("{"+', '.join([str(x) for x in XX0_BM])+"}"+"\n")
    f.close()  
    
    f=open(replicaFileA1,"a+")
    print('SAVING PLOTS A1>>  ',f.name)
    ### [total chi^2, list of NP-parameters],
    f.write("{"+', '.join([str(x) for x in XX1_BM])+"}"+"\n")
    f.close()  
    
    f=open(replicaFileA2,"a+")
    print('SAVING PLOTS A2>>  ',f.name)
    ### [total chi^2, list of NP-parameters],
    f.write("{"+', '.join([str(x) for x in XX2_BM])+"}"+"\n")
    f.close()  



