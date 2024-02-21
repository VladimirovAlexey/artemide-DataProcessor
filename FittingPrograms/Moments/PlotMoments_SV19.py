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
with open(SAVEPATH+'G0_ART23_mu_x.0000', 'w') as outfile:
    for i in range(81):
        for j in range(35):
            mu=(50)**(j/40)
            x=10**(-0.05*i)
            G0=harpy.get_uTMDPDF_G0(x,mu,1)
            G0str=[str(g) for g in G0]
            outfile.write(str(mu)+",   "+str(x)+",   "+",   ".join(G0str)+"\n")
   
#%%
for i in range(81):
    mu=10.
    x=10**(-0.05*i)
    #G0=harpy.get_uTMDPDF_G0(x,mu,1)
    pdf=harpy.get_uPDF(x,mu,1)
    #print('{'+str(mu)+','+str(x*G0[7])+','+str(x*pdf[7])+'},')
    print('{'+str(x)+','+str(pdf[6])+'},')

#%%
for i in range(41):
    x=0.01
    mu=(100)**(i/40)
    #mu=10.
    #x=10**(-i/10)
    X0=harpy.get_uTMDPDF_X0(x,mu,1)
    #pdf=harpy.get_uPDF(x,mu,1)
    print('{'+str(mu)+','+str(X0[3])+','+str(X0[4])+','+str(X0[6])+','+str(X0[7])+'},')
    
#%%
for i in range(41):
    #x=0.1
    #mu=(100)**(i/40)
    mu=10.
    x=10**(-i/10)
    X0=harpy.get_uPDF(x,mu,1)
    #pdf=harpy.get_uPDF(x,mu,1)
    #print('{'+str(mu)+','+str(x*G0[7])+','+str(x*pdf[7])+'},')
    print('{'+str(x)+','+str(X0[3])+','+str(X0[4])+','+str(X0[6])+','+str(X0[7])+'},')
    
#%%
for i in range(201):
    x=0.01
    b=100*(0.00000001)**(i/200)
    X0=harpy.get_uTMDPDF(x,b,1)
    #pdf=harpy.get_uPDF(x,mu,1)
    print('{'+str(b)+','+str(X0[3])+','+str(X0[4])+','+str(X0[6])+','+str(X0[7])+'},')
    
#%%
for i in range(201):
    x=0.01
    b=100*(0.00000001)**(i/200)
    mu=1.12292/b+5
    X0=harpy.get_uPDF(x,mu,1)
    #pdf=harpy.get_uPDF(x,mu,1)
    print('{'+str(b)+','+str(X0[3])+','+str(X0[4])+','+str(X0[6])+','+str(X0[7])+'},')