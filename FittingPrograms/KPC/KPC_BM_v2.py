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

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/data/WorkingFiles/TMD/Fit_Notes/ART23/REPS/ART23_run2.rep")
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
harpy.setNPparameters_BM_TMDPDF([0.1,0.1,0.1,1.])

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
    
    # #### let me cut the large-qT
    # if(p["<Q>"]<40. and p["<qT>"]>2.):
    #     return False,p
    # #### let me cut the large-X
    # if(p["<Q>"]<40. and p["y"][0]>.4):
    #     return False,p
    return (p["<qT>"]<10.1) , p

def cutFunc_angularPLOT(p):
    if(len(p["process"])==4):
        if(p["process"][3]==3): p["process"][3]=2
    else:
        print("UNKNOWN PROCESS IN ARTEMIDE 3"+str(p["process"]))
    
    delta=p["<qT>"]/p["<Q>"]
    return (delta<0.25) , p

def SetProcess(p):
    p["process"]=proc1
    p["thFactor"]=1.
    return True,p

#%%
### Load the data for Auu
theDataAuu=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_Auu_0y1","A8_Auu_1y2","A8_Auu_2y35"]))
setAuu=theDataAuu.CutData(cutFunc_angular) 
setAuuPLOT=theDataAuu.CutData(cutFunc_angularPLOT)

### Load the data for A0
theDataA0=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A0_0y1","A8_A0_1y2","A8_A0_2y35"]))
setA0=theDataA0.CutData(cutFunc_angular)
setA0PLOT=theDataA0.CutData(cutFunc_angularPLOT)
proc1=[1,1,1,20]
setA0ff=setA0.CutData(SetProcess)
setA0ffPLOT=setA0PLOT.CutData(SetProcess)
proc1=[1,1,1,30]
setA0hh=setA0.CutData(SetProcess)
setA0hhPLOT=setA0PLOT.CutData(SetProcess)

### Load the data for A1
theDataA1=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A1_0y1","A8_A1_1y2"]))
setA1=theDataA1.CutData(cutFunc_angular)
setA1PLOT=theDataA1.CutData(cutFunc_angularPLOT)
proc1=[1,1,1,21]
setA1ff=setA1.CutData(SetProcess)
setA1ffPLOT=setA1PLOT.CutData(SetProcess)
proc1=[1,1,1,31]
setA1hh=setA1.CutData(SetProcess)
setA1hhPLOT=setA1PLOT.CutData(SetProcess)

### Load the data for A2
theDataA2=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A2_0y1","A8_A2_1y2","A8_A2_2y35"]))
setA2=theDataA2.CutData(cutFunc_angular)
setA2PLOT=theDataA2.CutData(cutFunc_angularPLOT)
proc1=[1,1,1,22]
setA2ff=setA2.CutData(SetProcess)
setA2ffPLOT=setA2PLOT.CutData(SetProcess)
proc1=[1,1,1,32]
setA2hh=setA2.CutData(SetProcess)
setA2hhPLOT=setA2PLOT.CutData(SetProcess)

### Load the data for A3
theDataA3=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A3_0y1","A8_A3_1y2","A8_A3_2y35"]))
setA3=theDataA3.CutData(cutFunc_angular)
setA3PLOT=theDataA3.CutData(cutFunc_angularPLOT)

### Load the data for A4
theDataA4=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A4_0y1","A8_A4_1y2","A8_A4_2y35"]))
setA4=theDataA4.CutData(cutFunc_angular)
setA4PLOT=theDataA4.CutData(cutFunc_angularPLOT)

### Load the data for A5
theDataA5=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A5_0y1","A8_A5_1y2","A8_A5_2y35"]))
setA5=theDataA5.CutData(cutFunc_angular)
setA5PLOT=theDataA5.CutData(cutFunc_angularPLOT)

### Load the data for A6
theDataA6=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A6_0y1","A8_A6_1y2"]))
setA6=theDataA6.CutData(cutFunc_angular)
setA6PLOT=theDataA6.CutData(cutFunc_angularPLOT)

#%%
### Load the E866 data for Auu
theDataAuu=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_Qlow_nu","E866(p+d)_Qlow_nu"]))
    #"E866(p+p)_Qlow_dXF_nu","E866(p+d)_Qlow_dXF_nu"]))
setE1=theDataAuu.CutData(cutFunc_angular) 
theDataAuu=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_Qhigh_nu","E866(p+d)_Qhigh_nu"]))
    #"E866(p+p)_Qhigh_dXF_nu","E866(p+d)_Qhigh_dXF_nu"]))
setE2=theDataAuu.CutData(cutFunc_angular) 
theDataAuu=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_dQ_nu","E866(p+d)_dQ_nu"]))
setEdQ=theDataAuu.CutData(cutFunc_angular) 

#%%
def ReplaceProcess(p):
    #if(p["Q"]==[4.5, 9.0]): p["Q"]=[7.5, 9.0]
    p["process"]=proc1
    p["thFactor"]=1.
    #if(p["qT"][0]>1.5): return False,p
    return True, p

proc1=[2,1,1,2]
setE1_uu=setE1.CutData(ReplaceProcess)
setE2_uu=setE2.CutData(ReplaceProcess)
setEdQ_uu=setEdQ.CutData(ReplaceProcess)

proc1=[2,1,1,20]
setE1_0ff=setE1.CutData(ReplaceProcess)
setE2_0ff=setE2.CutData(ReplaceProcess)
setEdQ_0ff=setEdQ.CutData(ReplaceProcess)

proc1=[2,1,1,21]
setE1_1ff=setE1.CutData(ReplaceProcess)
setE2_1ff=setE2.CutData(ReplaceProcess)
setEdQ_1ff=setEdQ.CutData(ReplaceProcess)

proc1=[2,1,1,22]
setE1_2ff=setE1.CutData(ReplaceProcess)
setE2_2ff=setE2.CutData(ReplaceProcess)
setEdQ_2ff=setEdQ.CutData(ReplaceProcess)

proc1=[2,1,1,30]
setE1_0hh=setE1.CutData(ReplaceProcess)
setE2_0hh=setE2.CutData(ReplaceProcess)
setEdQ_0hh=setEdQ.CutData(ReplaceProcess)

proc1=[2,1,1,31]
setE1_1hh=setE1.CutData(ReplaceProcess)
setE2_1hh=setE2.CutData(ReplaceProcess)
setEdQ_1hh=setEdQ.CutData(ReplaceProcess)

proc1=[2,1,1,32]
setE1_2hh=setE1.CutData(ReplaceProcess)
setE2_2hh=setE2.CutData(ReplaceProcess)
setEdQ_2hh=setEdQ.CutData(ReplaceProcess)


#%%
# MAIN FIT
#harpy.setNPparameters([1.5004, 0.05614, 0.03862, 0.0, 0.565, 0.0539, 0.5697, 6.64, 0.565, 20.07, 0.5697, 0.537, 1.07, 2.39, 0.0, 0.0, 0.5, 0.2,0.7,3.])
#harpy.setNPparameters([1.500, 0.057, 0.059, 0.000, 0.475, 0.369, 0.829, 5.335, -0.357, 19.847, 1.866, 0.031, 1.046, 0.066, 0.002,0.000])
harpy.setNPparameters_TMDR([1.500, 0.0565, 0.06086, 0.0])
harpy.setNPparameters_uTMDPDF([0.4223, 0.16896, 0.56410, 8.387, 5.986, 16.026, 5.646, 0.3116, 1.0575, 0.01, 0.0, 0.0])
harpy.setNPparameters_BM_TMDPDF([0.011, -0.1, 4., 2.000])

#%%
###Computation of different contributions due to the unpolarized part
XA_uu=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setAuu,method="approximate"))
XA_2ff=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA2ff,method="approximate"))


###Computation of different contributions due to the unpolarized part
XA_uuPLOT=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setAuuPLOT,method="approximate"))
XA_2ffPLOT=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA2ffPLOT,method="approximate"))

#%%
###Computation of different contributions due to the unpolarized part
XX1_uu=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1_uu,method="approximate"))
XX2_uu=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE2_uu,method="approximate"))
XXdQ_uu=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setEdQ_uu,method="approximate"))

XX1_0ff=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1_0ff,method="approximate"))
XX2_0ff=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE2_0ff,method="approximate"))
XXdQ_0ff=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setEdQ_0ff,method="approximate"))

XX1_2ff=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1_2ff,method="approximate"))
XX2_2ff=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE2_2ff,method="approximate"))
XXdQ_2ff=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setEdQ_2ff,method="approximate"))

#%%
#######################################
# Minimisation
#######################################
import time

def chi_BM(x):
    startT=time.time()
    harpy.setNPparameters_BM_TMDPDF([x[0],x[1],x[2],x[3]])
    print('np set =',["{:8.3f}".format(i) for i in x])        
    
    chiA=chi_ATLAS() 
        
    chiE=0.#chi_E866()
    
    endT=time.time()
    print(':->',chiA/setA2.numberOfPoints, '+', chiE/(setE1_2hh.numberOfPoints+setEdQ_2hh.numberOfPoints),' =',
          chiA/setA2.numberOfPoints + chiE/(setE1_2hh.numberOfPoints+setEdQ_2hh.numberOfPoints),'       t=',endT-startT)
    return 10*chiA+chiE

def chi_ATLAS():    
    Y1=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA2hh,method="approximate"))
    X1=(Y1+XA_2ff)/XA_uu    
    # print("---ATLAS---")
    # print(", ".join([str(r) for r in X1]))
    ccA2,cc3=setA2hh.chi2(X1)
    return ccA2

def chi_E866():
    Y1=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1_2hh,method="approximate"))
    Y2=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE2_2hh,method="approximate"))
    Y3=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setEdQ_2hh,method="approximate"))
    R1=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1_0hh,method="approximate"))
    R2=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE2_0hh,method="approximate"))
    R3=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setEdQ_0hh,method="approximate"))
    
    X1=2*(XX1_2ff+XX2_2ff+Y1+Y2)/(2*(XX1_uu+XX2_uu)+XX1_0ff+XX2_0ff+R1+R2)
    # print("---E866 dQT---")
    # print(", ".join([str(r) for r in X1]))
    part1,cc3=setE1_uu.chi2(X1)
    
    
    X2=2*(XXdQ_2ff+Y3)/(2*XXdQ_uu+XXdQ_0ff+R3)
    # print("---E866 dQ---")
    # print(", ".join([str(r) for r in X2]))    
    part2,cc3=setEdQ_uu.chi2(X2)
    
    #print("--->",part1/setE1_2hh.numberOfPoints,part2/setEdQ_2hh.numberOfPoints)
        
    return part1+part2

#%%
#### Minimize DY
from iminuit import Minuit

#---- PDFbias-like row
initialValues=([0.2,-0.27,9.4,5.0])

initialErrors=(0.5,0.3,0.3,0.5)
searchLimits=((0.01,100.),(-100,100.),(1.,25.),(1.,100.))
# True= FIX
parametersToMinimize=(True, False,False,False)

#%%

m = Minuit(chi_BM, initialValues)

m.errors=initialErrors
m.limits=searchLimits
m.fixed=parametersToMinimize
m.errordef=1

print(m.params)

m.strategy=1
m.migrad()

print(m.params)

valsDY=list(m.values)

chi_BM(m.values)

print([round(x,1 if x >100 else 4 if x>1 else 6) for x in list(m.values)])

sys.exit()

#%%

setA0_BM=theDataA0.CutData(cutFunc_angularBMPLOT)
setA1_BM=theDataA1.CutData(cutFunc_angularBMPLOT)
setA2_BM=theDataA2.CutData(cutFunc_angularBMPLOT)

XX0_BM=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA0_BM))
XX1_BM=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA1_BM))
XX2_BM=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA2_BM))

print("---- RESULTING EXPRESSIONS---------")
print("--------XX0_BM---------")
print(XX0_BM)
print("--------XX1_BM---------")
print(XX1_BM)
print("--------XX2_BM---------")
print(XX2_BM)
print("--------------------")

#%%
print("Done.  =>     Create points & append to data set ...")
DataPP1=DataProcessor.DataSet.DataSet('nu p+p',"DY")

DataPP1.isNormalized=False
proc_current=[1,1,1,2]
s_current=500.**2#2*800*0.938
#Q_current=[[4.7,5.5],[5.5,6.5],[6.5,7.5],[7.5,8.5],[8.5,9.0],[10.7,15.]]
#Q_current=[4.5,4.6]
Q_current=[[10.+i,11.+i] for i in range(90)]
y_current=[0.,0.5]
ptBINS=[0.5,1.0]
#ptBINS=[[i*0.25,(i+1)*0.25] for i in range(16)]

for i in range(len(Q_current)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataPP1.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=ptBINS
    p["thFactor"]=1.
    p["Q"]=Q_current[i]
    p["y"]=y_current
    p["includeCuts"]=False
    p["xSec"]=0.1
    p["uncorrErr"].append(0.1)
    #
    DataPP1.AddPoint(p)        

print("Done.  ")

#%%
print("Done.  =>     Create points & append to data set ...")
DataPP1=DataProcessor.DataSet.DataSet('nu p+p',"DY")

DataPP1.isNormalized=False
proc_current=[1,1,1,32]
s_current=[(40.+90*i)**2 for i in range(90)]
Q_current=[[10.+i,11.+i] for i in range(90)]
y_current=[0.,0.5]
ptBINS=[0.5,1.0]
#ptBINS=[[i*0.25,(i+1)*0.25] for i in range(16)]

for i in range(len(Q_current)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataPP1.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current[i]
    p["qT"]=ptBINS
    p["thFactor"]=1.
    p["Q"]=Q_current[i]
    p["y"]=y_current
    p["includeCuts"]=False
    p["xSec"]=0.1
    p["uncorrErr"].append(0.1)
    #
    DataPP1.AddPoint(p)        

print("Done.  ")
#%%
for i in range(20):
    xx=harpy.get_uTMDPDF_kT(0.7,0.1*i,1,mu=6.5)
    print(xx[2],",")
#%%
setSET=DataPP1
harpy.setNPparameters_BM_TMDPDF([0.2,-0.27,9.4,0.])
print(DataProcessor.harpyInterface.ComputeXSec(setSET,method="approximate"))
harpy.setNPparameters_BM_TMDPDF([0.2,-0.27+0.34,9.4,0.])
print(DataProcessor.harpyInterface.ComputeXSec(setSET,method="approximate"))
harpy.setNPparameters_BM_TMDPDF([0.2,-0.27-0.12,9.4,0.])
print(DataProcessor.harpyInterface.ComputeXSec(setSET,method="approximate"))
harpy.setNPparameters_BM_TMDPDF([0.2,-0.27,9.4+5.4,0.])
print(DataProcessor.harpyInterface.ComputeXSec(setSET,method="approximate"))
harpy.setNPparameters_BM_TMDPDF([0.2,-0.27,9.4-0.9,0.])
print(DataProcessor.harpyInterface.ComputeXSec(setSET,method="approximate"))

#%%

print("Done.  =>     Create points & append to data set ...")
DataPP1=DataProcessor.DataSet.DataSet('nu p+p',"DY")

DataPP1.isNormalized=False
proc_current=[1,1,1,32]
s_current=8000.**2
Q_current=[81.,101.]
y_current=[[0.,1.0],[0.,1.0],[0.,1.0],[1.,2.1],[1.,2.1],[1.,2.1]]
ptBINS=[[0.,10.],[10.,20.],[20.,35.],[0.,10.],[10.,20.],[20.,35.]]

for i in range(len(ptBINS)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataPP1.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=ptBINS[i]
    p["thFactor"]=1.
    p["Q"]=Q_current
    p["y"]=y_current[i]
    p["includeCuts"]=False
    p["xSec"]=0.1
    p["uncorrErr"].append(0.1)
    #
    DataPP1.AddPoint(p)        

print("Done.  ")

#%%

print("Done.  =>     Create points & append to data set ...")
DataPP1=DataProcessor.DataSet.DataSet('nu p+p',"DY")

DataPP1.isNormalized=False
proc_current=[1,1,1,2]
s_current=8000.**2
Q_current=[80.,100.]
y_current=[0.,2.]
ptBINS=[[0.,2.5],[2.5,5.0],[5.0,8.0],[8.0,11.4],[11.4,14.9],[14.9,18.5],[18.5,22.0],[22.0,25.5],[25.5,29.0],[29.0,32.6]]

for i in range(len(ptBINS)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataPP1.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=ptBINS[i]
    p["thFactor"]=1.
    p["Q"]=Q_current
    p["y"]=y_current
    p["includeCuts"]=False
    p["xSec"]=0.1
    p["uncorrErr"].append(0.1)
    #
    DataPP1.AddPoint(p)        

print("Done.  ")

#%%

print("Done.  =>     Create points & append to data set ...")
DataPP1=DataProcessor.DataSet.DataSet('nu p+p',"DY")

DataPP1.isNormalized=False
proc_current=[1,1,1,23]
s_current=8000.**2
Q_current=[80.,100.]
y_current=[[0.1*i-3,0.1*i-2.9] for i in range(60)]
ptBINS=[2.5,5.0]

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