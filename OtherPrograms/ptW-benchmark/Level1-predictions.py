#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 16:42:03 2023

@author: alexey
"""

#######################################
#
# Set muOPE=C0_const*c4/bT+1._dp!+2d0
#
#######################################

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
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)

SAVEPATH="/data/WorkingFiles/TMD/arTeMiDe/Benchmark/aTMDe-tables-level1-MSHT/"

#%%
import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet

MAINPATH=ROOT_DIR 

#CASE="N4LL"
#CASE="N3LLp"
#CASE="N3LL"
#CASE="NNLLp"
#CASE="NNLL"
#CASE="NLL"
CASE="NLLp"
#CASE="LL"


#%%
#######################################
#Initialize artemide
#######################################
import harpy

path_to_constants=MAINPATH+"OtherPrograms/ptW-benchmark/const-files/level1_MSHT_"+CASE+".atmde"

harpy.initialize(path_to_constants)

harpy.setNPparameters([2., 0.01, 0., 0., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

#%%

##### This is universal point for lv1 Y=0
MZ=91.15348061
def MakePoint(Q,Y,qT):
    return {'type': 'DY',
     'id': 'A13-level1',  
     'process': [1, 5, 1, 1],
     's': 169000000.0,
     '<Q>': Q,
     'Q': [Q-0.1, Q+0.1],
     '<y>': Y,
     'y': [Y-0.001, Y+0.001],
     '<qT>': qT,
     'qT': [qT-0.01, qT+0.01],
     'thFactor': 1.0,
     'includeCuts': False,
     'cutParams': [20.0, 20.0, -2.4, 2.4]}

#%%
#### Compute cross-section by 
p=MakePoint(MZ,0.,1.)
X=DataProcessor.harpyInterface.ComputeXSec(p,method="central")
X_lvl1=4*p["<qT>"]*p["<Q>"]*X
print(X_lvl1)

#%%
TYPE="A"

if(TYPE=="A"):
    harpy.setNPparameters([2., 0.01, 0., 0., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10000.0, 10000.0, 0.0, 0.0])
else:
    harpy.setNPparameters([2., 0.01, 0., 0., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

#%%
yValues=[[0.,"Y=0"],[1.2,"Y=12"],[2.4,"Y=24"],[3.6,"Y=36"]]
QValues=[[66.,"Q=66"],[MZ,"Q=MZ"],[116.,"Q=116"],[300.,"Q=300"],[1000.,"Q=1000"]]
qtValues=[float(i) for i in range(1,41)]+[float(5*i+40) for i in range(1,13)]

for yy in yValues:
    for QQ in QValues:
        print(QQ[1]+'_'+yy[1]+'_'+CASE+'_'+TYPE)
        with open(SAVEPATH+'aTMDe_'+QQ[1]+'_'+yy[1]+'_'+CASE+'_'+TYPE+'.txt', 'w') as outfile:
            outfile.write("# artemide results for level1 resummation \n")
            outfile.write("# s=13000 GeV \n")
            outfile.write("# Q=   "+str(QQ[0])+"GeV \n")
            outfile.write("# Y=   "+str(yy[0])+"\n")
            outfile.write("# order="+CASE+"\n")
            if(TYPE=="A"):
                outfile.write("# type= A (only ud flavors) \n")
            else:
                outfile.write("# type= B (5 flavors) \n")
            outfile.write("#qT   |     2qT 2Q dsigma/dqT^2 /dQ^2/dY \n")
            for qT in qtValues:
                p=MakePoint(QQ[0],yy[0],qT)
                X=DataProcessor.harpyInterface.ComputeXSec(p,method="central")
                X_lvl1=4*p["<qT>"]*p["<Q>"]*X
                outfile.write("{:8.2f}".format(qT)+'    '+"{:16.15f}".format(X_lvl1)+"\n")

#%%
TYPE="B"


if(TYPE=="A"):
    harpy.setNPparameters([2., 0.01, 0., 0., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10000.0, 10000.0, 0.0, 0.0])
else:
    harpy.setNPparameters([2., 0.01, 0., 0., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

#%%
yValues=[[0.,"Y=0"],[1.2,"Y=12"],[2.4,"Y=24"],[3.6,"Y=36"]]
QValues=[[66.,"Q=66"],[MZ,"Q=MZ"],[116.,"Q=116"],[300.,"Q=300"],[1000.,"Q=1000"]]
qtValues=[float(i) for i in range(1,41)]+[float(5*i+40) for i in range(1,13)]

for yy in yValues:
    for QQ in QValues:
        print(QQ[1]+'_'+yy[1]+'_'+CASE+'_'+TYPE)
        with open(SAVEPATH+'aTMDe_'+QQ[1]+'_'+yy[1]+'_'+CASE+'_'+TYPE+'.txt', 'w') as outfile:
            outfile.write("# artemide results for level1 resummation \n")
            outfile.write("# s=13000 GeV \n")
            outfile.write("# Q=   "+str(QQ[0])+"GeV \n")
            outfile.write("# Y=   "+str(yy[0])+"\n")
            outfile.write("# order="+CASE+"\n")
            if(TYPE=="A"):
                outfile.write("# type= A (only ud flavors) \n")
            else:
                outfile.write("# type= B (5 flavors) \n")
            outfile.write("#qT   |     2qT 2Q dsigma/dqT^2 /dQ^2/dY \n")
            for qT in qtValues:
                p=MakePoint(QQ[0],yy[0],qT)
                X=DataProcessor.harpyInterface.ComputeXSec(p,method="central")
                X_lvl1=4*p["<qT>"]*p["<Q>"]*X
                outfile.write("{:8.2f}".format(qT)+'    '+"{:16.15f}".format(X_lvl1)+"\n")
