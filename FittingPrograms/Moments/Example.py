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

import sys
import numpy
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)

#%%
import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet

MAINPATH=ROOT_DIR 

#%%
#######################################
#Initialize artemide
#######################################
import harpy

path_to_constants=MAINPATH+"FittingPrograms/Moments/ConstantsFiles/ART23_MSHT_N4LL+.atmde"

harpy.initialize(path_to_constants)

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/data/WorkingFiles/TMD/Fit_Notes/ART23/REPS/ART23_run2.rep")
rSet.SetReplica(-1)

rSet.SetReplica(0)

#%%
#### Plot of the zeroth-moment vs. pdf
muValues=[i*1. for i in range(1,101)]
xValues=[10**(-i/10) for i in range(31)]

x=0.01
mu=50.


#for mu in muValues:
for x in xValues:
    G0=harpy.get_uTMDPDF_G0(x,mu,1)
    pdf=harpy.get_uPDF(x,mu,1)
    #print("{", f"{x:.4} ,{G0[5+1]:.8f}, {G0[5+2]:.8f}, {G0[5-1]:.8f}, {G0[5-2]:.8f}","},")
    print("{", f"{x:.4} ,{pdf[5+1]:.8f}, {pdf[5+2]:.8f}, {pdf[5-1]:.8f}, {pdf[5-2]:.8f}","},")


#%%

mu=4.
bValues=[i*0.03 for i in range(1,101)]

rSet.SetReplica(0,part="TMDR")
## central CS-kernel
central=[harpy.get_DNP(b, mu) for b in bValues]

dist=[]
for r in range(rSet.numberOfReplicas):
    rSet.SetReplica(r,part="TMDR")
    dist.append([harpy.get_DNP(b, mu) for b in bValues])

mean=numpy.mean(dist,axis=0)
std=numpy.std(dist,axis=0)

    #%%
    import matplotlib.pyplot as plt
    
    plt.plot(bValues,central,"r",bValues,mean,"b",bValues,mean+std,"b",bValues,mean-std,"b")
    plt.ylabel('CS-kernel')
    plt.show()