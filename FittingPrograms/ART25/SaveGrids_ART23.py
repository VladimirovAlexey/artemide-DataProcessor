#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 16:42:03 2023

@author: alexey
"""

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
import DataProcessor.SaveTMDGrid

MAINPATH=ROOT_DIR 

GRIDNAME1="ART23_uTMDPDF_b_Q=4"
GRIDNAME2="ART23_uTMDPDF_b_Q=10"
    
SAVEPATH1="/data/WorkingFiles/TMD/Fit_Notes/ART23/Grids/"+GRIDNAME1+"/"+GRIDNAME1
SAVEPATH2="/data/WorkingFiles/TMD/Fit_Notes/ART23/Grids/"+GRIDNAME2+"/"+GRIDNAME2

#%%
#######################################
#Initialize artemide
#######################################
import harpy

#####################################
#########   DO NOT FORGET TO RECOMPILE ARTEMIDE IN ART23 MODEL!!!!!
#####################################

path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/ART23_MSHT_N4LL_v301.atmde"
#path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/ART25_gluon.atmde"
harpy.initialize(path_to_constants)

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                        "/data/WorkingFiles/TMD/Fit_Notes/ART23/REPS/ART23_run2.rep")
    
rSet.SetReplica(0)
#%%

## I save only first 300 replicas (to save disk-space)
N_replicas_max=rSet.numberOfReplicas

#%%
#######################################
# Save Grid specification (optimal)
#######################################
with open(SAVEPATH1+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDPDF ART23 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:2302.07473"+"\n"),
    outfile.write("SetIndex: 706000"+"\n"),
    outfile.write("TMDScheme: Pavia TMDs"+"\n"),
    outfile.write("TMDType: pdf"+"\n"),
    outfile.write("CollDist: MSHT20_nnlo"+"\n"),
    outfile.write("CollDistMember:  0"+"\n"),
    outfile.write("Format: optimalTMD"+"\n"),    
    outfile.write("DataVersion: 1"+"\n"),
    outfile.write("OrderQCD: N4LL & zeta-prec."+"\n"),
    outfile.write("Regularisation:  xxx "+"\n"),
    outfile.write("NumMembers: "+str(N_replicas_max+1)+"\n"),
    outfile.write("ErrorType: Monte Carlo"+"\n"),
    outfile.write("FlavorScheme: LHAPDF style"+"\n"),
    outfile.write("Flavors: [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]"+"\n"),
    outfile.write("NumFlavors: 5"+"\n"),
    outfile.write("XMin:  0.00001"+"\n"),
    outfile.write("XMax: 1."+"\n"),
    outfile.write("QMin: 1."+"\n"),
    outfile.write("QMax: 200."+"\n"),
    outfile.write("KtoQMin: 0.0001"+"\n"),
    outfile.write("KtoQMax: 2.001"+"\n")      
    
#%%
#######################################
# Save Grid specification (optimal)
#######################################
with open(SAVEPATH2+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDPDF ART23 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:2302.07473"+"\n"),
    outfile.write("SetIndex: 706000"+"\n"),
    outfile.write("TMDScheme: Pavia TMDs"+"\n"),
    outfile.write("TMDType: pdf"+"\n"),
    outfile.write("CollDist: MSHT20_nnlo"+"\n"),
    outfile.write("CollDistMember:  0"+"\n"),
    outfile.write("Format: optimalTMD"+"\n"),    
    outfile.write("DataVersion: 1"+"\n"),
    outfile.write("OrderQCD: N4LL & zeta-prec."+"\n"),
    outfile.write("Regularisation:  xxx "+"\n"),
    outfile.write("NumMembers: "+str(N_replicas_max+1)+"\n"),
    outfile.write("ErrorType: Monte Carlo"+"\n"),
    outfile.write("FlavorScheme: LHAPDF style"+"\n"),
    outfile.write("Flavors: [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]"+"\n"),
    outfile.write("NumFlavors: 5"+"\n"),
    outfile.write("XMin:  0.00001"+"\n"),
    outfile.write("XMax: 1."+"\n"),
    outfile.write("QMin: 1."+"\n"),
    outfile.write("QMax: 200."+"\n"),
    outfile.write("KtoQMin: 0.0001"+"\n"),
    outfile.write("KtoQMax: 2.001"+"\n")  
    



#%%
rSet.SetReplica()
DataProcessor.SaveTMDGrid.SaveGrid_optimal(SAVEPATH1+"_0000.dat",Q=4.)

DataProcessor.SaveTMDGrid.SaveGrid_optimal(SAVEPATH2+"_0000.dat",Q=10.)

#%%
for i in range(N_replicas_max+1):
    print("Replica:"+ "{:04d}".format(i))
    rSet.SetReplica(i)
    DataProcessor.SaveTMDGrid.SaveGrid_optimal(SAVEPATH1+"_"+"{:04d}".format(i)+".dat",Q=4.)
    
    DataProcessor.SaveTMDGrid.SaveGrid_optimal(SAVEPATH2+"_"+"{:04d}".format(i)+".dat",Q=10.)