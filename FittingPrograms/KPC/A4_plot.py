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

replicaFileAuu =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/plotAuu.dat"
replicaFileA0 =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/plotA0ff.dat"
replicaFileA1 =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/plotA1ff.dat"
replicaFileA2 =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/plotA2ff.dat"
replicaFileA3 =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/plotA3.dat"
replicaFileA4 =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/plotA4.dat"
replicaFileNU =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/plotA4.dat"

import sys
import numpy
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)

### number of circles that this code will run
NumberOfReplicas=300


#%%
import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet
import DataProcessor.DataMultiSet

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/data/WorkingFiles/TMD/Fit_Notes/ART23/REPS/ART23_run2.rep")
MAINPATH=ROOT_DIR 
#%%
#######################################
#Initialize artemide
#######################################
import harpy

path_to_constants=MAINPATH+"FittingPrograms/KPC/INI/ART23_MSHT_N4LL.atmde"


harpy.initialize(path_to_constants)

rSet.SetReplica()
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
def cutFunc_angularPLOT(p):
    if(len(p["process"])==4):
        if(p["process"][3]==3): p["process"][3]=2
        #if(p["process"][3]==202): p["process"][3]=2
    else:
        print("UNKNOWN PROCESS IN ARTEMIDE 3"+str(p["process"]))
    
    delta=p["<qT>"]/p["<Q>"]
    return (delta<0.25) , p
#%%
### Load the data for Auu
theDataAuu=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_Auu_0y1","A8_Auu_1y2","A8_Auu_2y35"]))
setAuuPLOT=theDataAuu.CutData(cutFunc_angularPLOT) 

### Load the data for A0
theDataA0=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A0_0y1","A8_A0_1y2","A8_A0_2y35"]))
setA0PLOT=theDataA0.CutData(cutFunc_angularPLOT)

### Load the data for A1
theDataA1=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A1_0y1","A8_A1_1y2"]))
setA1PLOT=theDataA1.CutData(cutFunc_angularPLOT)

### Load the data for A2
theDataA2=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A2un_0y1","A8_A2un_1y2","A8_A2un_2y35"]))
setA2PLOT=theDataA2.CutData(cutFunc_angularPLOT)

### Load the data for A3
theDataA3=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A3_0y1","A8_A3_1y2","A8_A3_2y35"]))
setA3PLOT=theDataA3.CutData(cutFunc_angularPLOT)

### Load the data for A4
theDataA4=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A4_0y1","A8_A4_1y2","A8_A4_2y35"]))
setA4PLOT=theDataA4.CutData(cutFunc_angularPLOT)


#%%
#######################################
# This is the main cicle. 
# It generates replica of data take random PDF and minimize it
# Save to log.
#######################################
import time
for i in range(NumberOfReplicas):
    
    rnd=numpy.random.randint(rSet.numberOfReplicas)
    if(rnd==2): rnd=numpy.random.randint(rSet.numberOfReplicas)
    rSet.SetReplica(rnd)
    
    
    print('---------------------------------------------------------------')
    print('------------REPLICA ',i,' from ',NumberOfReplicas,'------------')
    print('------------     R= ',rnd, '------------')
    print('---------------------------------------------------------------')
    savedTime=time.time()
    
    #XXuu=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setAuuPLOT,method="approximate"))
    #XX0=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA0PLOT,method="approximate"))
    #XX1=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA1PLOT,method="approximate"))
    #XX2=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA2PLOT,method="approximate"))
    XX3=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA3PLOT,method="approximate"))
    XX4=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA4PLOT,method="approximate"))
    
    #A0=XX0/XXuu
    #A1=XX1/XXuu[0:14]
    #A2=XX2/XXuu
    A3=XX3
    A4=XX4
    
    savedTime2=time.time()
    
    print("Computation time:",savedTime2-savedTime)
    
    # f=open(replicaFileAuu,"a+")
    # ### [total chi^2, list of NP-parameters],
    # f.write(', '.join([str(x) for x in XXuu])+"\n")
    # f.close()  
    
    # f=open(replicaFileA0,"a+")
    # ### [total chi^2, list of NP-parameters],
    # f.write(', '.join([str(x) for x in A0])+"\n")
    # f.close()  
    
    # f=open(replicaFileA1,"a+")
    # ### [total chi^2, list of NP-parameters],
    # f.write(', '.join([str(x) for x in A1])+"\n")
    # f.close()  
    
    # f=open(replicaFileA2,"a+")
    # ### [total chi^2, list of NP-parameters],
    # f.write(', '.join([str(x) for x in A2])+"\n")
    # f.close()  
    
    f=open(replicaFileA3,"a+")
    ### [total chi^2, list of NP-parameters],
    f.write(', '.join([str(x) for x in A3])+"\n")
    f.close()  
    
    f=open(replicaFileA4,"a+")
    ### [total chi^2, list of NP-parameters],
    f.write(', '.join([str(x) for x in A4])+"\n")
    f.close()  


#%%

# print("Done.  =>     Create points & append to data set ...")
# DataAN=DataProcessor.DataSet.DataSet('A8-AN',"DY")
# DataAN.comment="ATLAS8"
# DataAN.reference="1606.00689"

# DataAN.isNormalized=False
# proc_current=[1,1,1,23]
# PTcase=[0.,2.5]
# yBINS=numpy.arange(-3.5, 3.6, 0.1)
# for j in range(len(yBINS)):
#     # makeup a point
#     p=DataProcessor.Point.CreateDYPoint(DataAN.name+'.'+str(j))
#     p["process"]=proc_current
#     p["s"]=64000000.0
#     p["qT"]=PTcase
#     p["thFactor"]=1.
#     p["Q"]=[80.0, 100.0]
#     p["<Q>"]=91.2 ### to be sure that mass middle value is Z-boson mass
#     p["y"]=[yBINS[j]-0.05,yBINS[j]+0.05]
#     p["includeCuts"]=False
#     p["xSec"]=0.1
#     p["uncorrErr"].append(0.1)
#     #
#     DataAN.AddPoint(p)

# print("Done.  ")
