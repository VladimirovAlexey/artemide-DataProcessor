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

replicaFileAuu =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/plotAuu_E866.dat"
replicaFileA0 =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/plotA0ff_E866.dat"
replicaFileA2 =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/plotA2ff_E866.dat"
replicaFileA0hh =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/plotA0hh_E866.dat"
replicaFileA2hh =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/plotA2hh_E866.dat"

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
    
    return True, p

#%%
### Load the E866 data for Auu
theDataAuu=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_Qlow_nu","E866(p+p)_Qlow_dXF_nu"]))
setE1=theDataAuu.CutData(cutFunc_angularPLOT) 
theDataAuu=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_Qhigh_nu","E866(p+p)_Qhigh_dXF_nu"]))
setE2=theDataAuu.CutData(cutFunc_angularPLOT) 
theDataAuu=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_dQ_nu"]))
setEdQ=theDataAuu.CutData(cutFunc_angularPLOT) 

#%%
def ReplaceProcess(p):
    #if(p["Q"]==[4.5, 9.0]): p["Q"]=[7.5, 9.0]
    p["process"]=proc1
    p["thFactor"]=1.
    return True, p

proc1=[2,1,1,2]
setE1_uu=setE1.CutData(ReplaceProcess)
setE2_uu=setE2.CutData(ReplaceProcess)
setEdQ_uu=setEdQ.CutData(ReplaceProcess)

proc1=[2,1,1,20]
setE1_0ff=setE1.CutData(ReplaceProcess)
setE2_0ff=setE2.CutData(ReplaceProcess)
setEdQ_0ff=setEdQ.CutData(ReplaceProcess)

proc1=[2,1,1,22]
setE1_2ff=setE1.CutData(ReplaceProcess)
setE2_2ff=setE2.CutData(ReplaceProcess)
setEdQ_2ff=setEdQ.CutData(ReplaceProcess)

proc1=[2,1,1,30]
setE1_0hh=setE1.CutData(ReplaceProcess)
setE2_0hh=setE2.CutData(ReplaceProcess)
setEdQ_0hh=setEdQ.CutData(ReplaceProcess)

proc1=[2,1,1,32]
setE1_2hh=setE1.CutData(ReplaceProcess)
setE2_2hh=setE2.CutData(ReplaceProcess)
setEdQ_2hh=setEdQ.CutData(ReplaceProcess)
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
    
    xE1_uu=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1_uu,method="approximate"))
    xE2_uu=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1_uu,method="approximate"))
    xEdQ_uu=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setEdQ_uu,method="approximate"))
    
    xE1_0ff=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1_0ff,method="approximate"))
    xE2_0ff=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1_0ff,method="approximate"))
    xEdQ_0ff=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setEdQ_0ff,method="approximate"))
    
    xE1_2ff=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1_2ff,method="approximate"))
    xE2_2ff=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1_2ff,method="approximate"))
    xEdQ_2ff=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setEdQ_2ff,method="approximate"))
    
    xE_uu=xE1_uu+xE2_uu
    xE_0ff=xE1_0ff+xE2_0ff
    xE_2ff=xE1_2ff+xE2_2ff
    
    Xuu=numpy.concatenate((xE_uu,xEdQ_uu))
    X0=numpy.concatenate((xE_0ff,xEdQ_0ff))
    X2=numpy.concatenate((xE_2ff,xEdQ_2ff))
    
    savedTime2=time.time()
    
    print("Computation time:",savedTime2-savedTime)
    
    f=open(replicaFileAuu,"a+")
    ### [total chi^2, list of NP-parameters],
    f.write(', '.join([str(x) for x in Xuu])+"\n")
    f.close()  
    
    f=open(replicaFileA0,"a+")
    ### [total chi^2, list of NP-parameters],
    f.write(', '.join([str(x) for x in X0])+"\n")
    f.close()  
    
    f=open(replicaFileA2,"a+")
    ### [total chi^2, list of NP-parameters],
    f.write(', '.join([str(x) for x in X2])+"\n")
    f.close()  

#%%

lambdaCENTRAL=numpy.array([1.2,-0.35,-0.8,5.])
variations=[
    numpy.array([0,0,0,0]), ### central case
    numpy.array([0.45,0,0,0]),
    numpy.array([-0.23,0,0,0]),
    numpy.array([0,0.11,0,0]),
    numpy.array([0,-0.10,0,0]),
    numpy.array([0,0,0.07,0]),
    numpy.array([0,0,-0.07,0]),
    numpy.array([0,0,0,1.8]),
    numpy.array([0,0,0,-0.8])]

for r in variations:
    harpy.setNPparameters_BM_TMDPDF(lambdaCENTRAL+r)

    xE1_0hh=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1_0hh,method="approximate"))
    xE2_0hh=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1_0hh,method="approximate"))
    xEdQ_0hh=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setEdQ_0hh,method="approximate"))

    xE1_2hh=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1_2hh,method="approximate"))
    xE2_2hh=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1_2hh,method="approximate"))
    xEdQ_2hh=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setEdQ_2hh,method="approximate"))

    xE_0hh=xE1_0hh+xE2_0hh
    xE_2hh=xE1_2hh+xE2_2hh

    X0h=numpy.concatenate((xE_0hh,xEdQ_0hh))
    X2h=numpy.concatenate((xE_2hh,xEdQ_2hh))

    f=open(replicaFileA0hh,"a+")
    ### [total chi^2, list of NP-parameters],
    f.write(', '.join([str(x) for x in X0h])+"\n")
    f.close()  

    f=open(replicaFileA2hh,"a+")
    ### [total chi^2, list of NP-parameters],
    f.write(', '.join([str(x) for x in X2h])+"\n")
    f.close()  
    
#%%
print("Done.  =>     Create points & append to data set ...")
DataPP1=DataProcessor.DataSet.DataSet('nu p+p',"DY")

DataPP1.isNormalized=False
proc_current=[1,1,1,2]
s_current=64000000.
Q_current=[80.,100.]
y_current=[[0.05*i,0.05*(i+1)] for i in range(70)]
ptBINS=[0.,2.5]

for i in range(len(y_current)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataPP1.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=ptBINS
    p["thFactor"]=1.
    p["Q"]=Q_current
    p["y"]=y_current[i]
    p["includeCuts"]=False
    p["xSec"]=0.1
    p["uncorrErr"].append(0.1)
    #
    DataPP1.AddPoint(p)        

print("Done.  ")
    