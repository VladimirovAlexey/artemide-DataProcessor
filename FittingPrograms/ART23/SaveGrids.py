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
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)

#%%
import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet
import DataProcessor.SaveTMDGrid

MAINPATH=ROOT_DIR 


#PDFtoUSE="MSHT"
PDFtoUSE="NNPDF"

if(PDFtoUSE=="MSHT"):
    GRIDNAME1="ART_optimal"
    GRIDNAME2="ART_kT"
else:
    GRIDNAME1="ART_"+PDFtoUSE+"_optimal"
    GRIDNAME2="ART_"+PDFtoUSE+"_kT"
    
SAVEPATH1="/data/WorkingFiles/TMD/Fit_Notes/ART23/Grids/"+GRIDNAME1+"/"+GRIDNAME1
SAVEPATH2="/data/WorkingFiles/TMD/Fit_Notes/ART23/Grids/"+GRIDNAME2+"/"+GRIDNAME2

#%%
#######################################
#Initialize artemide
#######################################
import harpy


path_to_constants=MAINPATH+"FittingPrograms/ART23/ConstantsFiles/DYonly/ART23_"+PDFtoUSE+"_N4LL"

harpy.initialize(path_to_constants)
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
#### In this model at very low Q, the zeta became negative!
#Qrange=[1.5*(200/1.5)**(n/40) for n in range(41)]
Qrange=DataProcessor.SaveTMDGrid.Qrange_default

#%%
#######################################
# Save Grid specification (optimal)
#######################################
with open(SAVEPATH1+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDPDF ART23 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:2302.?????"+"\n"),
    outfile.write("SetIndex: 706000"+"\n"),
    outfile.write("TMDScheme: Pavia TMDs"+"\n"),
    outfile.write("TMDType: pdf"+"\n"),
    outfile.write("CollDist: MSHT20_nnlo"+"\n"),
    outfile.write("CollDistMember:  0"+"\n"),
    outfile.write("Format: optimalTMD"+"\n"),    
    outfile.write("DataVersion: 1"+"\n"),
    outfile.write("OrderQCD: N4LL & zeta-prec."+"\n"),
    outfile.write("Regularisation:  xxx "+"\n"),
    outfile.write("NumMembers: "+str(rSet.numberOfReplicas+1)+"\n"),
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
# Save Grid specification (kT)
#######################################
with open(SAVEPATH2+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDPDF ART23 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:2302.?????"+"\n"),
    outfile.write("SetIndex: 706000"+"\n"),
    outfile.write("TMDScheme: Pavia TMDs"+"\n"),
    outfile.write("TMDType: pdf"+"\n"),
    outfile.write("CollDist: MSHT20_nnlo"+"\n"),
    outfile.write("CollDistMember:  0"+"\n"),    
    outfile.write("Format: TMDlib2"+"\n"),    
    outfile.write("DataVersion: 1"+"\n"),
    outfile.write("OrderQCD: N4LL & zeta-prec."+"\n"),
    outfile.write("Regularisation:  xxx "+"\n"),
    outfile.write("NumMembers: "+str(rSet.numberOfReplicas+1)+"\n"),
    outfile.write("ErrorType: Monte Carlo"+"\n"),
    outfile.write("FlavorScheme: LHAPDF style"+"\n"),
    outfile.write("Flavors: [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]"+"\n"),
    outfile.write("NumFlavors: 5"+"\n"),
    outfile.write("XMin:  0.00001"+"\n"),
    outfile.write("XMax: 1."+"\n"),
    outfile.write("QMin: "+str(Qrange[0])+"\n"),
    outfile.write("QMax: "+str(Qrange[-1])+"\n"),
    outfile.write("KtoQMin: 0.0001"+"\n"),
    outfile.write("KtoQMax: 2.001"+"\n")  
    
#%%
rSet.SetReplica()
DataProcessor.SaveTMDGrid.SaveGrid_optimal(SAVEPATH1+"_0000.dat")
DataProcessor.SaveTMDGrid.SaveGrid_kT(SAVEPATH2+"_0000.dat",Qrange=Qrange)



#%%
for i in range(rSet.numberOfReplicas+1):
    print("Replica:"+ "{:04d}".format(i))
    rSet.SetReplica(i)
    DataProcessor.SaveTMDGrid.SaveGrid_optimal(SAVEPATH1+"_"+"{:04d}".format(i)+".dat")
    DataProcessor.SaveTMDGrid.SaveGrid_kT(SAVEPATH2+"_"+"{:04d}".format(i)+".dat",Qrange=Qrange)