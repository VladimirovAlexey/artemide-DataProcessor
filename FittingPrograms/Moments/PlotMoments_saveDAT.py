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
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)

SAVEPATH="/data/WorkingFiles/TMD/Fit_Notes/MomentsOfTMDs/data/"

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

harpy.initialize(path_to_constants)

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/data/WorkingFiles/TMD/Fit_Notes/ART23/REPS/ART23_run2.rep")
rSet.SetReplica(-1)

rSet.SetReplica(0)

#%%
DIR_PLOT=SAVEPATH+"ART23/"

#%%
# # ##################################
# # ## Save replicas of xSec
# # ##################################


rSet.SetReplica()

nn=9

nmin=nn*100+1
nmax=nn*100+1+100
if(nmin>=rSet.numberOfReplicas):
    sys.exit()
if(nmax>=rSet.numberOfReplicas): nmax=rSet.numberOfReplicas

for i in range(nmin,nmax):      
    rSet.SetReplica(i)    
    with open(SAVEPATH+'G0_ART23_mu_x.'+str(i).zfill(4), 'w') as outfile:
        for i in range(81):
            for j in range(41):
                mu=(50)**(j/40)
                x=10**(-0.05*i)
                G0=harpy.get_uTMDPDF_G0(x,mu,1)
                G0str=[str(g) for g in G0]
                outfile.write(str(mu)+",   "+str(x)+",   "+",   ".join(G0str)+"\n")

sys.exit()