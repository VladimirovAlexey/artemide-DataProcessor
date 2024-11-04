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
sys.path.remove('/data/arTeMiDe_Repository/artemide/harpy')
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)

SAVEPATH="/data/WorkingFiles/TMD/arTeMiDe/Benchmark/artemide-TABLES/level2_KPC/"

#%%
import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet

MAINPATH=ROOT_DIR 

CASE="N4LL"
#CASE="N3LLp"
#CASE="N3LL"
#CASE="N2LLp"
#CASE="N2LL"
#CASE="NLLp"
#CASE="NLL"
###CASE="LL"


#%%
#######################################
#Initialize artemide
#######################################
import harpy

path_to_constants=MAINPATH+"OtherPrograms/ptW-benchmark/const-files_MSHT_KPC/level2_MSHT_"+CASE+".atmde"

harpy.initialize(path_to_constants)

harpy.setNPparameters([2., 0.07, 0., 0., 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.0, 0.04])

#%%

##### This is universal point for lv1 Y=0
MZ=91.15348061
def MakePoint(Q,Y,qT):
    return {'type': 'DY',
     'id': 'A13-level1',  
     'process': [1, 1, 1, 3],
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

#sys.exit()
#%%
harpy.setNPparameters([2., 0.07, 0., 0., 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.0, 0.04])

#%%
#yValues=[[0.,"Y=0"],[1.2,"Y=12"],[2.4,"Y=24"],[3.6,"Y=36"]]
#QValues=[[MZ,"Q=MZ"],[66.,"Q=66"],[116.,"Q=116"],[300.,"Q=300"],[1000.,"Q=1000"]]
yValues=[[0.,"Y0"]]
QValues=[[MZ,"QmZ"]]
qtValues=[float(0.5*i) for i in range(1,21)]+[float(i+10) for i in range(1,31)]+[float(5*i+40) for i in range(1,13)]
#%%
Xcentral=[]
Xc2UP=[]
Xc2DOWN=[]
Xc3UP=[]
Xc3DOWN=[]
Xc4UP=[]
Xc4DOWN=[]

for yy in yValues:
    for QQ in QValues:
        for qT in qtValues:
            print("qT:",qT)
            p=MakePoint(QQ[0],yy[0],qT)
            X=DataProcessor.harpyInterface.ComputeXSec(p,method="central")
            Xcentral.append(4*p["<qT>"]*p["<Q>"]*X)
       
print("C2 UP")
harpy.varyScales(1.,2.,1.,1.)
for yy in yValues:
    for QQ in QValues:
        for qT in qtValues:
            print("qT:",qT)
            p=MakePoint(QQ[0],yy[0],qT)
            X=DataProcessor.harpyInterface.ComputeXSec(p,method="central")
            Xc2UP.append(4*p["<qT>"]*p["<Q>"]*X)
print("C3 UP")
harpy.varyScales(1.,1.,2.,1.)
for yy in yValues:
    for QQ in QValues:
        for qT in qtValues:
            print("qT:",qT)
            p=MakePoint(QQ[0],yy[0],qT)
            X=DataProcessor.harpyInterface.ComputeXSec(p,method="central")
            Xc3UP.append(4*p["<qT>"]*p["<Q>"]*X)
print("C4 UP")
harpy.varyScales(1.,1.,1.,2.)
for yy in yValues:
    for QQ in QValues:
        for qT in qtValues:
            print("qT:",qT)
            p=MakePoint(QQ[0],yy[0],qT)
            X=DataProcessor.harpyInterface.ComputeXSec(p,method="central")
            Xc4UP.append(4*p["<qT>"]*p["<Q>"]*X)
print("C2 DOWN")
harpy.varyScales(1.,0.5,1.,1.)
for yy in yValues:
    for QQ in QValues:
        for qT in qtValues:
            print("qT:",qT)
            p=MakePoint(QQ[0],yy[0],qT)
            X=DataProcessor.harpyInterface.ComputeXSec(p,method="central")
            Xc2DOWN.append(4*p["<qT>"]*p["<Q>"]*X)
print("C3 DOWN")
harpy.varyScales(1.,1.,0.5,1.)
for yy in yValues:
    for QQ in QValues:
        for qT in qtValues:
            print("qT:",qT)
            p=MakePoint(QQ[0],yy[0],qT)
            X=DataProcessor.harpyInterface.ComputeXSec(p,method="central")
            Xc3DOWN.append(4*p["<qT>"]*p["<Q>"]*X)
print("C4 DOWN")
harpy.varyScales(1.,1.,1.,0.5)
for yy in yValues:
    for QQ in QValues:
        for qT in qtValues:
            print("qT:",qT)
            p=MakePoint(QQ[0],yy[0],qT)
            X=DataProcessor.harpyInterface.ComputeXSec(p,method="central")
            Xc4DOWN.append(4*p["<qT>"]*p["<Q>"]*X)
            
#%%
DeltaResumUP=numpy.max([numpy.array(Xc2UP)-numpy.array(Xcentral),numpy.array(Xc2DOWN)-numpy.array(Xcentral),
 numpy.array(Xc3UP)-numpy.array(Xcentral),numpy.array(Xc3DOWN)-numpy.array(Xcentral)],axis=0)
DeltaResumDOWN=numpy.min([numpy.array(Xc2UP)-numpy.array(Xcentral),numpy.array(Xc2DOWN)-numpy.array(Xcentral),
 numpy.array(Xc3UP)-numpy.array(Xcentral),numpy.array(Xc3DOWN)-numpy.array(Xcentral)],axis=0)

DeltaResum=numpy.max([DeltaResumUP,numpy.abs(DeltaResumDOWN)],axis=0)

DeltaFOUP=numpy.max([numpy.array(Xc4UP)-numpy.array(Xcentral),numpy.array(Xc4DOWN)-numpy.array(Xcentral)],axis=0)
DeltaFODOWN=numpy.min([numpy.array(Xc4UP)-numpy.array(Xcentral),numpy.array(Xc4DOWN)-numpy.array(Xcentral)],axis=0)

DeltaFO=numpy.max([DeltaFOUP,numpy.abs(DeltaFODOWN)],axis=0)
        
#%%
for yy in yValues:
    for QQ in QValues:
        print(QQ[1]+'_'+yy[1]+'_'+CASE)
        with open(SAVEPATH+'aTMDe_5gen_'+QQ[1]+'_'+yy[1]+'_'+CASE+'_KPC_level2_qt.txt', 'w') as outfile:
            outfile.write("# artemide results for level 2 resummation (with KPC) \n")
            outfile.write("# s = 13000 GeV \n")
            outfile.write("# Q = 91.15348061 GeV \n")
            outfile.write("# Y = 0.0 \n")
            outfile.write("# order = "+CASE+"\n")    
            outfile.write("#  \n")
            #outfile.write("#qT   |     2qT 2Q dsigma/dqT^2 /dQ^2/dY    |  Delta_resum(down)    |  Delta_resum(up)    |  Delta_FO(down)    |  Delta_FO(up)\n")            
            # for i in range(len(qtValues)):                
            #     outfile.write("{:8.2f}".format(qtValues[i])
            #                   +'    '+"{:14.8f}".format(Xcentral[i])
            #                   +'    '+"{:14.8f}".format(DeltaResumDOWN[i])
            #                   +'    '+"{:14.8f}".format(DeltaResumUP[i])
            #                   +'    '+"{:14.8f}".format(DeltaFODOWN[i])
            #                   +'    '+"{:14.8f}".format(DeltaFOUP[i])                              
            #                   +"\n")
            outfile.write("#qT   |     2qT 2Q dsigma/dqT^2 /dQ^2/dY    |  Delta_resum    |  Delta_FO\n")
            for i in range(len(qtValues)):                
                outfile.write("{:8.2f}".format(qtValues[i])
                              +'    '+"{:14.8f}".format(Xcentral[i])
                              +'    '+"{:14.8f}".format(DeltaResum[i])                              
                              +'    '+"{:14.8f}".format(DeltaFO[i])                              
                              +"\n")