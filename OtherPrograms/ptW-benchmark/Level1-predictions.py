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

ATMDE_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide-forB_v3/"
HARPY_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide-forB_v3/harpy/"

#replicaFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/replicaDY.dat"
#logFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/LOG.log"

import sys
import numpy
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)

SAVEPATH="/data/WorkingFiles/TMD/arTeMiDe/Benchmark/aTMDe-tables-level1-MSHT/"
#SAVEPATH="/data/WorkingFiles/TMD/arTeMiDe/Benchmark/aTMDe-tables-level1-NNPDF31_luxqed/"

#%%
import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet

MAINPATH=ROOT_DIR 

#CASE="N4LL"
#CASE="N3LLp"
CASE="N3LL"
#CASE="N2LLp"
#CASE="N2LL"
#CASE="NLLp"
#CASE="NLL"
#CASE="LL"


#%%
#######################################
#Initialize artemide
#######################################
import harpy

path_to_constants=MAINPATH+"OtherPrograms/ptW-benchmark/const-files_MSHT/level1_MSHT_"+CASE+".atmde"
#path_to_constants=MAINPATH+"OtherPrograms/ptW-benchmark/const-files_NNPDF/level1_MSHT_"+CASE+".atmde"

harpy.initialize(path_to_constants)

harpy.setNPparameters([2., 0.01, 0., 0., 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.0, 0.04])

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
     'thFactor': 1.0/(1+0.5*(qT/Q)**2),
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
harpy.setNPparameters([1., 0.01, 0., 0., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0000.0, 0000.0, 0.0, 0.04])

#%%
#yValues=[[0.,"Y=0"],[1.2,"Y=12"],[2.4,"Y=24"],[3.6,"Y=36"]]
#QValues=[[MZ,"Q=MZ"],[66.,"Q=66"],[116.,"Q=116"],[300.,"Q=300"],[1000.,"Q=1000"]]
yValues=[[0.,"Y=0"]]
QValues=[[MZ,"Q=MZ"]]
qtValues=[float(0.5*i) for i in range(1,21)]+[float(i+10) for i in range(1,31)]+[float(5*i+40) for i in range(1,13)]
#%%
for yy in yValues:
    for QQ in QValues:
        print(QQ[1]+'_'+yy[1]+'_'+CASE)
        with open(SAVEPATH+'aTMDe_5gen_'+QQ[1]+'_'+yy[1]+'_'+CASE+'_level1_qt.txt', 'w') as outfile:
            outfile.write("# artemide results for level1 resummation \n")
            outfile.write("# s=13000 GeV \n")
            outfile.write("# Q=   "+str(QQ[0])+"GeV \n")
            outfile.write("# Y=   "+str(yy[0])+"\n")
            outfile.write("# order="+CASE+"\n")
            #if(TYPE=="A"):
            #    outfile.write("# type= A (only ud flavors) \n")
            #else:
            #    outfile.write("# type= B (5 flavors) \n")
            outfile.write("#qT   |     2qT 2Q dsigma/dqT^2 /dQ^2/dY \n")
            for qT in qtValues:
                p=MakePoint(QQ[0],yy[0],qT)
                X=DataProcessor.harpyInterface.ComputeXSec(p,method="central")
                X_lvl1=4*p["<qT>"]*p["<Q>"]*X
                outfile.write("{:8.2f}".format(qT)+'    '+"{:16.15f}".format(X_lvl1)+"\n")
                
#%%
##### SOME TABLES SAVED IN THE IMPROPER FORMAT --- REWRITE THEM
# SAVEPATH="/data/WorkingFiles/TMD/arTeMiDe/Benchmark/artemide-TABLES/level1/"
# CASE="LL"
# qtValues=[float(0.5*i) for i in range(1,21)]+[float(i+10) for i in range(1,31)]+[float(5*i+40) for i in range(1,13)]
# SS=[1.057945569843493 , 1.407122840082827 , 1.705268445187438 , \
# 1.916779369015682 , 2.061204385250858 , 2.153902189016873 , \
# 2.209237835340472 , 2.236546581573661 , 2.24301859744204 , \
# 2.234221839716419 , 2.214266260640891 , 2.186347461142204 , \
# 2.152495894016526 , 2.114749834147287 , 2.07416304545802 , \
# 2.031726131230393 , 1.988110839349305 , 1.944338953052895 , \
# 1.900117347597822 , 1.856225250455634 , 1.770383996553802 , \
# 1.688043493881229 , 1.609574528587811 , 1.535637818347288 , \
# 1.465743846512492 , 1.399546601598154 , 1.337462184651356 , \
# 1.278800635715146 , 1.223559863173406 , 1.171264843299709 , \
# 1.122078733115184 , 1.075849547951143 , 1.032213386886211 , \
# 0.990755865939085 , 0.951570153134404 , 0.914999113288931 , \
# 0.880237106605131 , 0.84729059526217 , 0.815944422570062 , \
# 0.786388303870785 , 0.758220004962594 , 0.731365045348567 , \
# 0.706065530243766 , 0.681492998897908 , 0.65818230134145 , \
# 0.635706093454947 , 0.614542443073051 , 0.594607907174769 , \
# 0.575351508563528 , 0.557339350709924 , 0.477146851514277 , \
# 0.411906721808112 , 0.357330410649508 , 0.312349987627678 , \
# 0.274265548200241 , 0.242207288136401 , 0.214552408710354 , \
# 0.190049472270804 , 0.169661654341702 , 0.151102604538083 , \
# 0.135259450416435 , 0.120931457344256 ]
   

# with open(SAVEPATH+'aTMDe_5gen_QmZ_Y0_'+CASE+'_level1_qt.txt', 'w') as outfile:
#     outfile.write("# artemide results for level 1 resummation \n")
#     outfile.write("# s = 13000 GeV \n")
#     outfile.write("# Q = 91.15348061 GeV \n")
#     outfile.write("# Y = 0.0 \n")
#     outfile.write("# order = "+CASE+"\n")    
#     outfile.write("#  \n")
#     outfile.write("#qT   |     2qT 2Q dsigma/dqT^2 /dQ^2/dY \n")
#     for i in range(len(qtValues)):        
#         outfile.write("{:8.2f}".format(qtValues[i])+'    '+"{:14.8f}".format(SS[i])+"\n")