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

replicaFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/replicaDY.dat"
logFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/LOG.log"

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

#PDFtoUSE="MSHT"
PDFtoUSE="NNPDF"

path_to_constants=MAINPATH+"FittingPrograms/ART23/ConstantsFiles/DYonly/ART23_"+PDFtoUSE+"_N4LL"

harpy.initialize(path_to_constants)
#%%
#%%
if(PDFtoUSE=="MSHT"):
    rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/data/WorkingFiles/TMD/Fit_Notes/ART23/REPS/ART23_run2.rep")
elif(PDFtoUSE=="NNPDF"):
    rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/data/WorkingFiles/TMD/Fit_Notes/ART23/REPS/ART23_NNPDF.rep")
else:
    raise Exception("no no no")
rSet.SetReplica(-1)

rSet.SetReplica(0)

#%%

bValues=[0.055, 0.105, 0.155, 0.205, 0.255, 0.305, 0.355, 0.405, 0.455, 0.505, \
0.555, 0.605, 0.655, 0.705, 0.755, 0.805, 0.855, 0.905, 0.955, 1.005, \
1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, \
4., 4.2, 4.4, 4.6, 4.8, 5.]

rSet.SetReplica(0, part="TMDR")

#%%
reps=[]
for i in range(rSet.numberOfReplicas+1):
    rSet.SetReplica(i, part="TMDR")
    rr=[harpy.get_DNP(b, 2.) for b in bValues]
    reps.append(rr)

reps=numpy.array(reps)
#%%
central=numpy.mean(reps,axis=0)

#%%
ss=[]
for i in range(200):
    idx=numpy.random.randint(rSet.numberOfReplicas+1, size=500)
    sample=reps[idx,:]
    up=numpy.quantile(sample, 1-0.68/2,axis=0) 
    down=numpy.quantile(sample, 0.68/2,axis=0) 
    ss.append([down,up])
    
#%%
band=numpy.mean(ss,axis=0)

#%%
print([[bValues[i],central[i],down[i],up[i]] for i in range(len(bValues))])