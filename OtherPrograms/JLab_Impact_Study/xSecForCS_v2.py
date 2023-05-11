#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 10:38:36 2022

@author: vla18041
"""

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

import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet

#%%
#######################################
#Initialize artemide
#######################################
import harpy
path_to_constants="/data/arTeMiDe_Repository/DataProcessor/OtherPrograms/JLab_Impact_Study/PDFb_N4LL"


harpy.initialize(path_to_constants)

harpy.setNPparameters_TMDR([1.5004, 0.049098, 0.05979, 0.0])
#harpy.setNPparameters_TMDR([1.93, 0.0434, 0, 0])
harpy.setNPparameters_uTMDPDF([0.834557, 0.914247, 0.910747, 4.5973, 
                               0.004487, 38.5017, 0.001313, 1.2705, 
                               1.1989, 0.173397, 0.0, 0.0])
harpy.setNPparameters_uTMDFF([0.264,0.479,0.459,0.539]) 

originalSV19TMDR=[1.93, 0.0434, 0, 0]
originalSV19FF=[0.264,0.479,0.459,0.539]

#%%
def CreateDataSet(Qbin,xBin,zBin,ptList):
    DataN=DataProcessor.DataSet.DataSet("temp","SIDIS")
    DataN.comment="N"
    DataN.reference="Harut"
    
    proc_current=[1,1,2001]
    #s_current=2*22.*0.938+(0.938)**2
    s_current=2*42.*0.938+(0.938)**2 ### NOT JLAB!
    includeCuts=False
    cutParameters=[0.1,0.85,10.,10000.]
    
    for i in range(len(ptList)):
        # makeup a point
        p=DataProcessor.Point.CreateSIDISPoint(DataN.name+'.'+str(i))
        #print DataCurrent.name+'.'+str(i)
        p["process"]=proc_current
        p["s"]=s_current
        p["pT"]=ptList[i]        
        p["Q"]=Qbin
        p["x"]=xBin
        p["z"]=zBin
        ## cross-seciton is unknown
        p["xSec"]=0.1
        p["M_target"]=0.938
        p["M_product"]=0.139
        p["includeCuts"]=includeCuts
        p["cutParams"]=cutParameters
        #devide by (x,z,pt^2,Q^2) bin size
        #p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["Q"][1]**2-p["Q"][0]**2)/(p["z"][1]-p["z"][0])/(p["x"][1]-p["x"][0])
        p["thFactor"]=1.
        
        ### to get the subtracted cross-section
        from DataProcessor.Point import FinalizePoint
        FinalizePoint(p)
        y=p["<Q>"]**2/p["<x>"]/(p["s"]-0*p["M_target"]**2)
        gamma=0.*2.*p["<x>"]*p["M_target"]/p["<Q>"]
        eps=(1.-y-(gamma*y)**2/4.)/(1.-y+y**2/2.+(gamma*y)**2/4)
        p["thFactor"]=(1.-eps)/(y**2)/(p["pT"][1]**2-p["pT"][0]**2)/(p["Q"][1]**2-p["Q"][0]**2)/(p["z"][1]-p["z"][0])/(p["x"][1]-p["x"][0])
        
        DataN.AddPoint(p)
    
    print("Done.  ")
    
    return DataN
#%%
DUMPPATH="/data/WorkingFiles/TMD/Fit_Notes/CSdirect/aTMDeTest/"

#%%
for j in range(20):    
    Q=1.5+0.5*j
    print(Q) 

    ss=CreateDataSet([Q-0.05, Q+0.05],[0.56,0.58],[0.48,0.52],[[0.04*i,0.04*(i+1)] for i in range(150)])
    res=DataProcessor.harpyInterface.ComputeXSec(ss)
    
    dd={
        'qT':[[p["pT"][0]/p["<z>"],p["pT"][1]/p["<z>"]] for p in ss.points],
        'Q':ss.points[0]["Q"],
        'X':res
        }
    from json import dump
    with open(DUMPPATH+"case_Q="+str(Q), 'w') as fp:
        dump(dd, fp)
