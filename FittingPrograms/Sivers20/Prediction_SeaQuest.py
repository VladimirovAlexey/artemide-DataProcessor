#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 16:42:03 2023

@author: alexey
"""

#######################################
# importing libraries
#######################################
import os
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"

ATMDE_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide/"
HARPY_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide/harpy/"

#replicaFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/replicaDY.dat"
#logFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/LOG.log"

import sys
import numpy
#sys.path.remove('/data/arTeMiDe_Repository/artemide/harpy')
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)

SAVEPATH="/data/WorkingFiles/TMD/Fit_Notes/MomentsOfTMDs/data/SV19/"

#%%
import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet

MAINPATH=ROOT_DIR 

#%%
#######################################
#Initialize artemide
#######################################
import harpy

path_to_constants=MAINPATH+"FittingPrograms/Moments/ConstantsFiles/ART23_MSHT_N4LL.atmde"
#path_to_constants=MAINPATH+"FittingPrograms/Moments/ConstantsFiles/SV19_NNPDF_N3LL.atmde"

harpy.initialize(path_to_constants)

### SV19 (all=0) n3lo
#harpy.setNPparameters_TMDR([2., 0.0442327])
#harpy.setNPparameters_uTMDPDF([0.17975, 3.9081, 453.883, 2.07199, 1.60774, 0., 0.,1.,1.])

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/data/WorkingFiles/TMD/Fit_Notes/ART23/REPS/ART23_run2.rep")

rSet.SetReplica(0)
    
rSetS=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            ATMDE_DIR+"Models/BPV20/Replica-files/BPV20(n3lo).rep")
rSetS.SetReplica(-1)

rSetS.SetReplica(0)
#%%
##### This is universal point for lv1 Y=0
def MakePoint(Q,Y,qT):
    return {'type': 'DY',
     'id': 'A13-level1',  
     'process': [2, 10001, 1, 1],
     's': 15.**2,
     '<Q>': 5.47723,
     'Q': [5.47723-0.1, 5.47723+0.1],
     '<y>': Y,
     'y': [Y-0.025, Y+0.025],
     '<qT>': 1.0,
     'qT': [1.0-0.1, 1.0+0.1],
     'thFactor': 1.0,
     'includeCuts': False,
     'cutParams': [20.0, 20.0, -2.4, 2.4]}

def MakePointUn(Q,Y,qT):
    return {'type': 'DY',
     'id': 'A13-level1',  
     'process': [2, 1, 1, 1],
     's': 15.**2,
     '<Q>': 5.47723,
     'Q': [5.47723-0.1, 5.47723+0.1],
     '<y>': Y,
     'y': [Y-0.025, Y+0.025],
     '<qT>': 1.0,
     'qT': [1.0-0.1, 1.0+0.1],
     'thFactor': 1.0,
     'includeCuts': False,
     'cutParams': [20.0, 20.0, -2.4, 2.4]}

#%%
#### Compute cross-section by 
p=MakePointUn(0.1,-0.4,1.)
X=DataProcessor.harpyInterface.ComputeXSec(p)
X_lvl1=4*p["<qT>"]*p["<Q>"]*X
print(X_lvl1)

#%%
#########################################################
## Resample data with N/2
#########################################################
def Resample(dd):
    #return numpy.random.choice(dd,size=int(numpy.floor(len(dd)/2)))
    return dd[numpy.random.choice(dd.shape[0], size=int(numpy.floor(len(dd)/2)))]

#########################################################
## Determine Mean, Mode, 68%CI by resampling 
#########################################################
alpha=68
def Compute68CI(dd):    
    lowers=[]
    uppers=[]    
    for i in range(1500):
        sample=Resample(numpy.array(dd))
        lowers.append(numpy.percentile(sample,(100-alpha)/2))
        uppers.append(numpy.percentile(sample,100-(100-alpha)/2))
    
    return [numpy.mean(lowers),numpy.mean(uppers)]

#%%
yValues=[i*0.05+0.025-0.7 for i in range(16)]

rSet.SetReplica(0)

result=[]
for y in yValues:
    pS=MakePoint(15.,y,1.)
    ss=[]
    for r in range(rSetS.numberOfReplicas):
        rSetS.SetReplica(r)
        ss.append(DataProcessor.harpyInterface.ComputeXSec(pS,method="binless"))
        
    Xsivers=Compute68CI(ss)
    pU=MakePointUn(15.,y,1.)
    Xun=DataProcessor.harpyInterface.ComputeXSec(pU,method="binless")*300
    result.append([y,numpy.mean(ss)/Xun,Xsivers[0]/Xun,Xsivers[1]/Xun])
    print(y)
    
[print("{"+str(r[0])+","+str(r[1])+","+str(r[2])+","+str(r[3])+"},") for r in result]
    