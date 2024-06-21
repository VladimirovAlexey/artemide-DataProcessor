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
def cutFunc(p):
    par=0.5
    
    #  for artemide v3.    
    # p["process"]=[p["process"][0],p["process"][2],1,1]
    if(len(p["process"])==4):
        if(p["process"][3]==3): p["process"][3]=2
    if(len(p["process"])==3):
        if(p["process"][2]==1): p["process"]=[p["process"][0],1,1,1]
        elif(p["process"][2]==2): p["process"]=[p["process"][0],1,-1,1]
        elif(p["process"][2]==3): print("ERROR1")
        elif(p["process"][2]==4): print("ERROR1")
        elif(p["process"][2]==5): p["process"]=[p["process"][0],1,1,2]
        elif(p["process"][2]==6): p["process"]=[p["process"][0],1,-1,2]
        elif(p["process"][2]==7): print("ERRORW")
        elif(p["process"][2]==8): print("ERRORW")
        elif(p["process"][2]==9): print("ERRORW")
        elif(p["process"][2]==10): print("ERRORW")
        elif(p["process"][2]==11): print("ERRORW")
        elif(p["process"][2]==12): print("ERRORW")
        elif(p["process"][2]==1001): p["process"]=[p["process"][0],1,1,101]
        elif(p["process"][2]==1002): p["process"]=[p["process"][0],1,1,102]
        else:
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
    
    #ART23
    return ((delta<0.25 and p["<qT>"]<10.) or (delta<0.25 and par/err*delta**2<1)) , p
    #return ((delta<0.25 and p["<qT>"]<10.)) , p
    #return (delta<0.25) , p
    
#%%
##################Cut function
def cutFunc_angular(p):
    if(len(p["process"])==4):
        if(p["process"][3]==3): p["process"][3]=2
    else:
        print("UNKNOWN PROCESS IN ARTEMIDE 3"+str(p["process"]))
    
    return (p["<qT>"]<10.) , p

def cutFunc_angularPLOT(p):
    if(len(p["process"])==4):
        if(p["process"][3]==3): p["process"][3]=2
    else:
        print("UNKNOWN PROCESS IN ARTEMIDE 3"+str(p["process"]))
    
    delta=p["<qT>"]/p["<Q>"]
    return (delta<0.25) , p

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
                          #'CMS13_dQ_106to170','CMS13_dQ_170to350','CMS13_dQ_350to1000',
                          'LHCb7', 'LHCb8', 
                          'LHCb13_dy(2021)', 
                          'PHE200', 'STAR510', 
                          'E228-200', 'E228-300', 'E228-400', 
                          'E772','E605'
                          #'D0run1-W','CDFrun1-W'
                          ]))

setDY=theData.CutData(cutFunc) 

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')

#%%
### Load the data for Auu
theDataAuu=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_Auu_0y1","A8_Auu_1y2","A8_Auu_2y35"]))
setAuu=theDataAuu.CutData(cutFunc_angular) 
setAuuPLOT=theDataAuu.CutData(cutFunc_angularPLOT) 

### Load the data for A0
theDataA0=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A0_0y1","A8_A0_1y2","A8_A0_2y35"]))
setA0=theDataA0.CutData(cutFunc_angular)
setA0PLOT=theDataA0.CutData(cutFunc_angularPLOT)
setA0_BM=theDataA0.CutData(cutFunc_angularBM)
setA0PLOT_BM=theDataA0.CutData(cutFunc_angularBMPLOT)

### Load the data for A1
theDataA1=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A1_0y1","A8_A1_1y2"]))
setA1=theDataA1.CutData(cutFunc_angular)
setA1PLOT=theDataA1.CutData(cutFunc_angularPLOT)
setA1_BM=theDataA1.CutData(cutFunc_angularBM)
setA1PLOT_BM=theDataA1.CutData(cutFunc_angularBMPLOT)

### Load the data for A2
theDataA2=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A2un_0y1","A8_A2un_1y2","A8_A2un_2y35"]))
setA2=theDataA2.CutData(cutFunc_angular)
setA2PLOT=theDataA2.CutData(cutFunc_angularPLOT)
setA2_BM=theDataA2.CutData(cutFunc_angularBM)
setA2PLOT_BM=theDataA2.CutData(cutFunc_angularBMPLOT)

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
setA5_BM=theDataA5.CutData(cutFunc_angularBM)
setA5PLOT_BM=theDataA5.CutData(cutFunc_angularBMPLOT)

### Load the data for A6
theDataA6=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A6_0y1","A8_A6_1y2"]))
setA6_BM=theDataA6.CutData(cutFunc_angularBM)
setA6_BMPLOT=theDataA6.CutData(cutFunc_angularBMPLOT)

#%%
### Load the E866 data for Auu
theDataAuu=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_Qlow_uu","E866(p+d)_Qlow_uu",
    "E866(p+p)_Qlow_dXF_uu","E866(p+d)_Qlow_dXF_uu"]))
setE1uu=theDataAuu.CutData(cutFunc_angular) 
theDataAuu=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_Qhigh_uu","E866(p+d)_Qhigh_uu",
    "E866(p+p)_Qhigh_dXF_uu","E866(p+d)_Qhigh_dXF_uu"]))
setE2uu=theDataAuu.CutData(cutFunc_angular) 
theDataAuu=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_dQ_uu","E866(p+d)_dQ_uu"]))
setEdQuu=theDataAuu.CutData(cutFunc_angular) 

theDataA=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_Qlow_A0ff","E866(p+d)_Qlow_A0ff",
    "E866(p+p)_Qlow_dXF_A0ff","E866(p+d)_Qlow_dXF_A0ff"]))
setE1A0ff=theDataA.CutData(cutFunc_angular) 
theDataA=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_Qhigh_A0ff","E866(p+d)_Qhigh_A0ff",
    "E866(p+p)_Qhigh_dXF_A0ff","E866(p+d)_Qhigh_dXF_A0ff"]))
setE2A0ff=theDataA.CutData(cutFunc_angular) 
theDataA=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_dQ_A0ff","E866(p+d)_dQ_A0ff"]))
setEdQA0ff=theDataA.CutData(cutFunc_angular) 

theDataA=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_Qlow_A2ff","E866(p+d)_Qlow_A2ff",
    "E866(p+p)_Qlow_dXF_A2ff","E866(p+d)_Qlow_dXF_A2ff"]))
setE1A2ff=theDataA.CutData(cutFunc_angular) 
theDataA=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_Qhigh_A2ff","E866(p+d)_Qhigh_A2ff",
    "E866(p+p)_Qhigh_dXF_A2ff","E866(p+d)_Qhigh_dXF_A2ff"]))
setE2A2ff=theDataA.CutData(cutFunc_angular) 
theDataA=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_dQ_A2ff","E866(p+d)_dQ_A2ff"]))
setEdQA2ff=theDataA.CutData(cutFunc_angular)

theDataA=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_Qlow_A0hh","E866(p+d)_Qlow_A0hh",
    "E866(p+p)_Qlow_dXF_A0hh","E866(p+d)_Qlow_dXF_A0hh"]))
setE1A0hh=theDataA.CutData(cutFunc_angular) 
theDataA=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_Qhigh_A0hh","E866(p+d)_Qhigh_A0hh",
    "E866(p+p)_Qhigh_dXF_A0hh","E866(p+d)_Qhigh_dXF_A0hh"]))
setE2A0hh=theDataA.CutData(cutFunc_angular) 
theDataA=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_dQ_A0hh","E866(p+d)_dQ_A0hh"]))
setEdQA0hh=theDataA.CutData(cutFunc_angular) 

theDataA=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_Qlow_A2hh","E866(p+d)_Qlow_A2hh",
    "E866(p+p)_Qlow_dXF_A2hh","E866(p+d)_Qlow_dXF_A2hh"]))
setE1A2hh=theDataA.CutData(cutFunc_angular) 
theDataA=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_Qhigh_A2hh","E866(p+d)_Qhigh_A2hh",
    "E866(p+p)_Qhigh_dXF_A2hh","E866(p+d)_Qhigh_dXF_A2hh"]))
setE2A2hh=theDataA.CutData(cutFunc_angular) 
theDataA=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular([
    "E866(p+p)_dQ_A2hh","E866(p+d)_dQ_A2hh"]))
setEdQA2hh=theDataA.CutData(cutFunc_angular) 

#%%
# MAIN FIT
harpy.setNPparameters([1.5004, 0.05614, 0.03862, 0.0, 0.565, 0.0539, 0.5697, 6.64, 0.565, 20.07, 0.5697, 0.537, 1.07, 2.39, 0.0, 0.0, 0.01, -2.3,-0.067,1.])

#%%
### Computation of different contributions due to the unpolarized part
# XXuu=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setAuu))
# XX0=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA0))
# XX1=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA1))
# XX2=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA2))

# print("---- NORMING EXPRESSIONS FOR ATLAS---------")
# print("--------Xuu---------")
# print(XXuu)
# print("--------XX0---------")
# print(XX0)
# print("--------XX1---------")
# print(XX1)
# print("--------XX2---------")
# print(XX2)
# print("--------------------")

XXuu=numpy.array([12.81688882, 23.6235875,  25.54355697, 22.58280105, 12.62825366, 22.97866139,
 24.41624181, 21.2264618,  15.57972008, 27.41510682, 28.15863863, 23.62314219])
XX2=numpy.array([0.00134844, 0.01452126, 0.05464758, 0.12218165, 0.00140259, 0.01489148,
 0.05510013, 0.12136979, 0.00193242, 0.01977701, 0.0706681,  0.15142411])

#%% 
### PLOTS
# XXuu=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setAuuPLOT))
# XX0=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA0PLOT))
# XX1=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA1PLOT))
# XX2=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA2PLOT))
# XX3=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA3PLOT))
# XX4=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA4PLOT))

# print("---- NORMING EXPRESSIONS FOR ATLAS---------")
# print("--------Xuu---------")
# print(", ".join([str(x) for x in XXuu]))
# print("--------XX0---------")
# print(", ".join([str(x) for x in XX0]))
# print("--------XX1---------")
# print(", ".join([str(x) for x in XX1]))
# print("--------XX2---------")
# print(", ".join([str(x) for x in XX2]))
# print("--------XX3---------")
# print(", ".join([str(x) for x in XX3]))
# print("--------XX4---------")
# print(", ".join([str(x) for x in XX4]))
# print("--------------------")
#%%
###Computation of different contributions due to the unpolarized part
EEuu1=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1uu))
EEuu2=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE2uu))
EE01=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1A0ff))
EE21=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE2A0ff))
EE02=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1A2ff))
EE22=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE2A2ff))
EEuudQ=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setEdQuu))
EEA0dQ=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setEdQA0ff))
EEA2dQ=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setEdQA2ff))
EEuu=EEuu1+EEuu2
EEA0=EE01+EE02
EEA2=EE21+EE22


print("---- NORMING EXPRESSIONS FOR E866---------")
print("--------Xuu---------")
print(repr(EEuu))
print("--------XA0---------")
print(repr(EEA0))
print("--------XA2---------")
print(repr(EEA2))
print("--------------------")
print("--------Xuu_DQ---------")
print(repr(EEuudQ))
print("--------XA0_DQ---------")
print(repr(EEA0dQ))
print("--------XA2_DQ---------")
print(repr(EEA2dQ))
print("--------------------")


# EEuu=numpy.array([ 4.14671824, 10.61416931, 13.16400926, 12.31345612, 22.05868784,  4.14671824,
#  10.61416931, 13.16400926, 12.31345612, 22.05868784])
# EEA0=numpy.array([ 0.40653886, 1.24786528, 2.03160484, 2.50291939, 5.89717844, 0.40653886,
#  1.24786528, 2.03160484, 2.50291939, 5.89717844 ])
# EEA2=numpy.array([-0.00127722, -0.0036834,  -0.00569674, -0.00719406, -0.03602629, -0.00127722,
#  -0.0036834,  -0.00569674, -0.00719406, -0.03602629])

EEuu=numpy.array([1.61380481,  4.13888023,  5.16631501,  4.89574618,  9.07861257,
        1.61380481,  4.13888023,  5.16631501,  4.89574618,  9.07861257,
       11.98782748,  4.6172226 ,  3.7118367 ,  2.51886283,  2.05790751,
       11.98782748,  4.6172226 ,  3.7118367 ,  2.51886283,  2.05790751])
EEA0=numpy.array([0.08747206,  0.31107721,  0.59693076,  0.83639688,  2.23128692,
        0.08747206,  0.31107721,  0.59693076,  0.83639688,  2.23128692,
        2.97768727,  0.71148009,  0.3774653 ,  0.12608373, -0.12991693,
        2.97768727,  0.71148009,  0.3774653 ,  0.12608373, -0.12991693])
EEA2=numpy.array([-0.00091161, -0.00263399, -0.00408834, -0.00518829, -0.02637688,
       -0.00091161, -0.00263399, -0.00408834, -0.00518829, -0.02637688,
       -0.01164717, -0.00560417, -0.00589987, -0.00568882, -0.01035757,
       -0.01164717, -0.00560417, -0.00589987, -0.00568882, -0.01035757])
EEuudQ=numpy.array([15.45411749,  5.78149003,  2.32295547,  0.97416151,  0.25216414,
        0.10969687, 15.45411749,  5.78149003,  2.32295547,  0.97416151,
        0.25216414,  0.10969687])
EEA0dQ=numpy.array([ 1.61319579,  0.12261887, -0.1071475 , -0.1012104 , -0.03752075,
       -0.04519323,  1.61319579,  0.12261887, -0.1071475 , -0.1012104 ,
       -0.03752075, -0.04519323])
EEA2dQ=numpy.array([1.71565575, 0.55888501, 0.20297059, 0.07781236, 0.01870619,
       0.0059946 , 1.71565575, 0.55888501, 0.20297059, 0.07781236,
       0.01870619, 0.0059946 ])

#%%
#######################################
# Minimisation
#######################################
import time

def chi_BM(x):
    startT=time.time()
    harpy.setNPparameters_BM_TMDPDF([x[0],x[1],x[2],x[3]])
    print('np set =',["{:8.3f}".format(i) for i in x])        
    
    chiA=0.#chi_ATLAS()
    chiE=chi_E866()
    
    endT=time.time()
    print(':->',chiA/setA2_BM.numberOfPoints, '+', chiE/setE1A2hh.numberOfPoints,' =',
          chiA/setA2_BM.numberOfPoints + chiE/setE1A2hh.numberOfPoints,'       t=',endT-startT)
    return chiA+chiE

def chi_ATLAS():
    Y2=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setA2_BM,method="approximate"))
    ccA2,cc3=setA2.chi2((Y2+XX2)/XXuu)    
    return ccA2

def chi_E866():
    Y1=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1A2hh,method="approximate"))
    Y2=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE2A2hh,method="approximate"))
    R1=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE1A0hh,method="approximate"))
    R2=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setE2A0hh,method="approximate"))
    AA2=(Y1+Y2+EEA2)/EEuu
    AA0=(R1+R2+EEA0)/EEuu
           
    #print(", ".join([str(r) for r in ((2*AA2)/(2+AA0))]))
    part1,cc3=setE1A2hh.chi2((2*AA2)/(2+AA0))    
    
    Y1=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setEdQA2hh,method="approximate"))
    R1=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setEdQA0hh,method="approximate"))
    AA2=(Y1+EEA2dQ)/EEuudQ
    AA0=(R1+EEA0dQ)/EEuudQ
    
    part2,cc3=setEdQA2hh.chi2((2*AA2)/(2+AA0))    
        
    return part1+part2

#%%
#### Minimize DY
from iminuit import Minuit

#---- PDFbias-like row
initialValues=([.3, 0.22, 0.4, 2.0])

initialErrors=(0.5,0.3,0.3,0.1)
searchLimits=((0.01,100.),(-100,100.),(-0.95,10.),(0.5,20.))
# True= FIX
parametersToMinimize=(False, False,False,True)

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
proc_current=[2,1,1,30]
s_current=2*800*0.938
Q_current=[[4.7,5.5],[5.5,6.5],[6.5,7.5],[7.5,8.5],[8.5,9.0],[10.7,15.]]
y_current=[0.,0.8]
ptBINS=[0.,4.]

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
harpy.setNPparameters_BM_TMDPDF([0.5,0.118,0.739,3.])
print(DataProcessor.harpyInterface.ComputeXSec(DataPP1,method="approximate"))
harpy.setNPparameters_BM_TMDPDF([0.5,0.118+0.03,0.739,3.])
print(DataProcessor.harpyInterface.ComputeXSec(DataPP1,method="approximate"))
harpy.setNPparameters_BM_TMDPDF([0.5,0.118-0.04,0.739,3.])
print(DataProcessor.harpyInterface.ComputeXSec(DataPP1,method="approximate"))
harpy.setNPparameters_BM_TMDPDF([0.5,0.118,0.739+0.52,3.])
print(DataProcessor.harpyInterface.ComputeXSec(DataPP1,method="approximate"))
harpy.setNPparameters_BM_TMDPDF([0.5,0.118,0.739-0.32,3.])
print(DataProcessor.harpyInterface.ComputeXSec(DataPP1,method="approximate"))