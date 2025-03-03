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

GRIDNAME1="ART25_uTMDPDF_optimal"
GRIDNAME2="ART25_uTMDPDF_kT"
GRIDNAME3="ART25_uTMDPDF_Q"
GRIDNAME4="ART25_uTMDPDF_optimal_kT"
GRIDNAME5="ART25_uTMDPDF_Q=2_kT"

GRIDNAME11="ART25_uTMDFF_pi_optimal"
GRIDNAME21="ART25_uTMDFF_pi_kT"
GRIDNAME31="ART25_uTMDFF_pi_Q"
GRIDNAME41="ART25_uTMDFF_pi_optimal_kT"
GRIDNAME51="ART25_uTMDFF_pi_Q=2_kT"

GRIDNAME12="ART25_uTMDFF_K_optimal"
GRIDNAME22="ART25_uTMDFF_K_kT"
GRIDNAME32="ART25_uTMDFF_K_Q"
GRIDNAME42="ART25_uTMDFF_K_optimal_kT"
GRIDNAME52="ART25_uTMDFF_K_Q=2_kT"
    
SAVEPATH1="/data/WorkingFiles/TMD/Fit_Notes/ART25/TMDgrids/"+GRIDNAME1+"/"+GRIDNAME1
SAVEPATH2="/data/WorkingFiles/TMD/Fit_Notes/ART25/TMDgrids/"+GRIDNAME2+"/"+GRIDNAME2
SAVEPATH3="/data/WorkingFiles/TMD/Fit_Notes/ART25/TMDgrids/"+GRIDNAME3+"/"+GRIDNAME3
SAVEPATH4="/data/WorkingFiles/TMD/Fit_Notes/ART25/TMDgrids/"+GRIDNAME4+"/"+GRIDNAME4
SAVEPATH5="/data/WorkingFiles/TMD/Fit_Notes/ART25/TMDgrids/"+GRIDNAME5+"/"+GRIDNAME5

SAVEPATH11="/data/WorkingFiles/TMD/Fit_Notes/ART25/TMDgrids/"+GRIDNAME11+"/"+GRIDNAME11
SAVEPATH21="/data/WorkingFiles/TMD/Fit_Notes/ART25/TMDgrids/"+GRIDNAME21+"/"+GRIDNAME21
SAVEPATH31="/data/WorkingFiles/TMD/Fit_Notes/ART25/TMDgrids/"+GRIDNAME31+"/"+GRIDNAME31
SAVEPATH41="/data/WorkingFiles/TMD/Fit_Notes/ART25/TMDgrids/"+GRIDNAME41+"/"+GRIDNAME41
SAVEPATH51="/data/WorkingFiles/TMD/Fit_Notes/ART25/TMDgrids/"+GRIDNAME51+"/"+GRIDNAME51

SAVEPATH12="/data/WorkingFiles/TMD/Fit_Notes/ART25/TMDgrids/"+GRIDNAME12+"/"+GRIDNAME12
SAVEPATH22="/data/WorkingFiles/TMD/Fit_Notes/ART25/TMDgrids/"+GRIDNAME22+"/"+GRIDNAME22
SAVEPATH32="/data/WorkingFiles/TMD/Fit_Notes/ART25/TMDgrids/"+GRIDNAME32+"/"+GRIDNAME32
SAVEPATH42="/data/WorkingFiles/TMD/Fit_Notes/ART25/TMDgrids/"+GRIDNAME42+"/"+GRIDNAME42
SAVEPATH52="/data/WorkingFiles/TMD/Fit_Notes/ART25/TMDgrids/"+GRIDNAME52+"/"+GRIDNAME52

#%%
#######################################
#Initialize artemide
#######################################
import harpy

path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/ART25_main.atmde"
harpy.initialize(path_to_constants)

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                        "/data/WorkingFiles/TMD/Fit_Notes/ART25/REPLICAS/ART25_run1.rep")
    
rSet.SetReplica(0)
#%%
#### In this model at very low Q, the zeta became negative!
#Qrange=[1.5*(200/1.5)**(n/40) for n in range(41)]
Qrange=DataProcessor.SaveTMDGrid.Qrange_default

#%%

## I save only first 300 replicas (to save disk-space)
N_replicas_max=500

#%%
#######################################
# Save Grid specification (optimal)
#######################################
with open(SAVEPATH1+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDPDF ART25 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:250?.????"+"\n"),
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
# Save Grid specification (kT)
#######################################
with open(SAVEPATH2+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDPDF ART25 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:230?.?????"+"\n"),
    outfile.write("SetIndex: 706000"+"\n"),
    outfile.write("TMDScheme: Pavia TMDs"+"\n"),
    outfile.write("TMDType: pdf"+"\n"),
    outfile.write("CollDist: MSHT20_nnlo"+"\n"),
    outfile.write("CollDistMember:  0"+"\n"),    
    outfile.write("Format: TMDlib2"+"\n"),    
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
    outfile.write("QMin: "+str(Qrange[0])+"\n"),
    outfile.write("QMax: "+str(Qrange[-1])+"\n"),
    outfile.write("KtoQMin: 0.0001"+"\n"),
    outfile.write("KtoQMax: 2.001"+"\n")  
    
#%%
#######################################
# Save Grid specification (Q)
#######################################
with open(SAVEPATH3+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDPDF ART25 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:25??.????"+"\n"),
    outfile.write("SetIndex: 706000"+"\n"),
    outfile.write("TMDScheme: Pavia TMDs"+"\n"),
    outfile.write("TMDType: pdf"+"\n"),
    outfile.write("CollDist: MSHT20_nnlo"+"\n"),
    outfile.write("CollDistMember:  0"+"\n"),    
    outfile.write("Format: TMDlib2"+"\n"),    
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
    outfile.write("QMin: "+str(Qrange[0])+"\n"),
    outfile.write("QMax: "+str(Qrange[-1])+"\n"),
    outfile.write("KtoQMin: 0.0001"+"\n"),
    outfile.write("KtoQMax: 2.001"+"\n") 
    
#%%
#######################################
# Save Grid specification (optimal) in kT
#######################################
with open(SAVEPATH4+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDPDF ART25 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:250?.????"+"\n"),
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
# Save Grid specification (optimal) in kT
#######################################
with open(SAVEPATH5+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDPDF ART25 at Q=2 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:250?.????"+"\n"),
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
with open(SAVEPATH11+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDFF for pion+ ART25 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:250?.????"+"\n"),
    outfile.write("SetIndex: 706000"+"\n"),
    outfile.write("TMDScheme: Pavia TMDs"+"\n"),
    outfile.write("TMDType: pdf"+"\n"),
    outfile.write("CollDist: MAPFF10NNLOPIp"+"\n"),
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
# Save Grid specification (kT)
#######################################
with open(SAVEPATH21+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDFF for pion+ ART25 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:230?.?????"+"\n"),
    outfile.write("SetIndex: 706000"+"\n"),
    outfile.write("TMDScheme: Pavia TMDs"+"\n"),
    outfile.write("TMDType: pdf"+"\n"),
    outfile.write("CollDist: MAPFF10NNLOPIp"+"\n"),
    outfile.write("CollDistMember:  0"+"\n"),    
    outfile.write("Format: TMDlib2"+"\n"),    
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
    outfile.write("QMin: "+str(Qrange[0])+"\n"),
    outfile.write("QMax: "+str(Qrange[-1])+"\n"),
    outfile.write("KtoQMin: 0.0001"+"\n"),
    outfile.write("KtoQMax: 2.001"+"\n")  
    
#%%
#######################################
# Save Grid specification (Q)
#######################################
with open(SAVEPATH31+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDFF for pion+ ART25 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:25??.????"+"\n"),
    outfile.write("SetIndex: 706000"+"\n"),
    outfile.write("TMDScheme: Pavia TMDs"+"\n"),
    outfile.write("TMDType: pdf"+"\n"),
    outfile.write("CollDist: MAPFF10NNLOPIp"+"\n"),
    outfile.write("CollDistMember:  0"+"\n"),    
    outfile.write("Format: TMDlib2"+"\n"),    
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
    outfile.write("QMin: "+str(Qrange[0])+"\n"),
    outfile.write("QMax: "+str(Qrange[-1])+"\n"),
    outfile.write("KtoQMin: 0.0001"+"\n"),
    outfile.write("KtoQMax: 2.001"+"\n") 
  
#%%
#######################################
# Save Grid specification (optimal) in kT
#######################################
with open(SAVEPATH41+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDFF for pion+ ART25 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:250?.????"+"\n"),
    outfile.write("SetIndex: 706000"+"\n"),
    outfile.write("TMDScheme: Pavia TMDs"+"\n"),
    outfile.write("TMDType: pdf"+"\n"),
    outfile.write("CollDist: MAPFF10NNLOPIp"+"\n"),
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
# Save Grid specification (optimal) in kT
#######################################
with open(SAVEPATH51+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDFF for pion+ ART25 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:250?.????"+"\n"),
    outfile.write("SetIndex: 706000"+"\n"),
    outfile.write("TMDScheme: Pavia TMDs"+"\n"),
    outfile.write("TMDType: pdf"+"\n"),
    outfile.write("CollDist: MAPFF10NNLOPIp"+"\n"),
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
with open(SAVEPATH12+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDFF for K+ ART25 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:250?.????"+"\n"),
    outfile.write("SetIndex: 706000"+"\n"),
    outfile.write("TMDScheme: Pavia TMDs"+"\n"),
    outfile.write("TMDType: pdf"+"\n"),
    outfile.write("CollDist: MAPFF10NNLOKAp"+"\n"),
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
# Save Grid specification (kT)
#######################################
with open(SAVEPATH22+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDFF for K+ ART25 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:230?.?????"+"\n"),
    outfile.write("SetIndex: 706000"+"\n"),
    outfile.write("TMDScheme: Pavia TMDs"+"\n"),
    outfile.write("TMDType: pdf"+"\n"),
    outfile.write("CollDist: MAPFF10NNLOKAp"+"\n"),
    outfile.write("CollDistMember:  0"+"\n"),    
    outfile.write("Format: TMDlib2"+"\n"),    
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
    outfile.write("QMin: "+str(Qrange[0])+"\n"),
    outfile.write("QMax: "+str(Qrange[-1])+"\n"),
    outfile.write("KtoQMin: 0.0001"+"\n"),
    outfile.write("KtoQMax: 2.001"+"\n")  
    
#%%
#######################################
# Save Grid specification (Q)
#######################################
with open(SAVEPATH32+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDFF for K+ ART25 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:25??.????"+"\n"),
    outfile.write("SetIndex: 706000"+"\n"),
    outfile.write("TMDScheme: Pavia TMDs"+"\n"),
    outfile.write("TMDType: pdf"+"\n"),
    outfile.write("CollDist: MAPFF10NNLOKAp"+"\n"),
    outfile.write("CollDistMember:  0"+"\n"),    
    outfile.write("Format: TMDlib2"+"\n"),    
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
    outfile.write("QMin: "+str(Qrange[0])+"\n"),
    outfile.write("QMax: "+str(Qrange[-1])+"\n"),
    outfile.write("KtoQMin: 0.0001"+"\n"),
    outfile.write("KtoQMax: 2.001"+"\n") 
    
#%%
#######################################
# Save Grid specification (optimal) in kT
#######################################
with open(SAVEPATH42+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDFF for K+ ART25 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:250?.????"+"\n"),
    outfile.write("SetIndex: 706000"+"\n"),
    outfile.write("TMDScheme: Pavia TMDs"+"\n"),
    outfile.write("TMDType: pdf"+"\n"),
    outfile.write("CollDist: MAPFF10NNLOKAp"+"\n"),
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
# Save Grid specification (optimal) in kT
#######################################
with open(SAVEPATH52+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDFF for K+ ART25 \n")
    outfile.write("Authors: V.Moos, I.Scimemi, A.Vladimirov, P.Zurita"+"\n"),
    outfile.write("Reference: arXiv:250?.????"+"\n"),
    outfile.write("SetIndex: 706000"+"\n"),
    outfile.write("TMDScheme: Pavia TMDs"+"\n"),
    outfile.write("TMDType: pdf"+"\n"),
    outfile.write("CollDist: MAPFF10NNLOKAp"+"\n"),
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
DataProcessor.SaveTMDGrid.SaveGrid_optimal(SAVEPATH1+"_0000.dat")
DataProcessor.SaveTMDGrid.SaveGrid_kT(SAVEPATH2+"_0000.dat",Qrange=Qrange)
DataProcessor.SaveTMDGrid.SaveGrid_Q(SAVEPATH3+"_0000.dat",Qrange=Qrange)
DataProcessor.SaveTMDGrid.SaveGrid_optimal_kT(SAVEPATH4+"_0000.dat")
DataProcessor.SaveTMDGrid.SaveGrid_optimal_kT(SAVEPATH5+"_0000.dat",Q=2.)

DataProcessor.SaveTMDGrid.SaveGrid_optimal(SAVEPATH11+"_0000.dat",
                                      Xrange=DataProcessor.SaveTMDGrid.XrangeFF_default,PDF="uTMDFF",h=1)
DataProcessor.SaveTMDGrid.SaveGrid_kT(SAVEPATH21+"_0000.dat",Qrange=Qrange,
                                      Xrange=DataProcessor.SaveTMDGrid.XrangeFF_default,PDF="uTMDFF",h=1)
DataProcessor.SaveTMDGrid.SaveGrid_Q(SAVEPATH31+"_0000.dat",Qrange=Qrange,
                                     Xrange=DataProcessor.SaveTMDGrid.XrangeFF_default,PDF="uTMDFF",h=1)
DataProcessor.SaveTMDGrid.SaveGrid_optimal_kT(SAVEPATH41+"_0000.dat",
                                     Xrange=DataProcessor.SaveTMDGrid.XrangeFF_default,PDF="uTMDFF",h=1)
DataProcessor.SaveTMDGrid.SaveGrid_optimal_kT(SAVEPATH51+"_0000.dat",
                                     Xrange=DataProcessor.SaveTMDGrid.XrangeFF_default,PDF="uTMDFF",h=1,Q=2.)

DataProcessor.SaveTMDGrid.SaveGrid_optimal(SAVEPATH12+"_0000.dat",
                                      Xrange=DataProcessor.SaveTMDGrid.XrangeFF_default,PDF="uTMDFF",h=2)
DataProcessor.SaveTMDGrid.SaveGrid_kT(SAVEPATH22+"_0000.dat",Qrange=Qrange,
                                      Xrange=DataProcessor.SaveTMDGrid.XrangeFF_default,PDF="uTMDFF",h=2)
DataProcessor.SaveTMDGrid.SaveGrid_Q(SAVEPATH32+"_0000.dat",Qrange=Qrange,
                                     Xrange=DataProcessor.SaveTMDGrid.XrangeFF_default,PDF="uTMDFF",h=2)
DataProcessor.SaveTMDGrid.SaveGrid_optimal_kT(SAVEPATH42+"_0000.dat",
                                     Xrange=DataProcessor.SaveTMDGrid.XrangeFF_default,PDF="uTMDFF",h=2)
DataProcessor.SaveTMDGrid.SaveGrid_optimal_kT(SAVEPATH52+"_0000.dat",
                                     Xrange=DataProcessor.SaveTMDGrid.XrangeFF_default,PDF="uTMDFF",h=2,Q=2.)

#%%
for i in range(N_replicas_max+1):
    print("Replica:"+ "{:04d}".format(i))
    rSet.SetReplica(i)
    DataProcessor.SaveTMDGrid.SaveGrid_optimal(SAVEPATH1+"_"+"{:04d}".format(i)+".dat")
    DataProcessor.SaveTMDGrid.SaveGrid_optimal(SAVEPATH11+"_"+"{:04d}".format(i)+".dat",
                                               Xrange=DataProcessor.SaveTMDGrid.XrangeFF_default,PDF="uTMDFF",h=1)
    DataProcessor.SaveTMDGrid.SaveGrid_optimal(SAVEPATH12+"_"+"{:04d}".format(i)+".dat",
                                               Xrange=DataProcessor.SaveTMDGrid.XrangeFF_default,PDF="uTMDFF",h=2)
    #DataProcessor.SaveTMDGrid.SaveGrid_kT(SAVEPATH2+"_"+"{:04d}".format(i)+".dat",Qrange=Qrange)
    #DataProcessor.SaveTMDGrid.SaveGrid_Q(SAVEPATH3+"_"+"{:04d}".format(i)+".dat",Qrange=Qrange)
    
#%%
for i in range(N_replicas_max+1):
    print("Replica:"+ "{:04d}".format(i))
    rSet.SetReplica(i)
    DataProcessor.SaveTMDGrid.SaveGrid_optimal_kT(SAVEPATH4+"_"+"{:04d}".format(i)+".dat")
    DataProcessor.SaveTMDGrid.SaveGrid_optimal_kT(SAVEPATH41+"_"+"{:04d}".format(i)+".dat",
                                               Xrange=DataProcessor.SaveTMDGrid.XrangeFF_default,PDF="uTMDFF",h=1)
    DataProcessor.SaveTMDGrid.SaveGrid_optimal_kT(SAVEPATH42+"_"+"{:04d}".format(i)+".dat",
                                               Xrange=DataProcessor.SaveTMDGrid.XrangeFF_default,PDF="uTMDFF",h=2)
    
#%%
for i in range(N_replicas_max+1):
    print("Replica:"+ "{:04d}".format(i))
    rSet.SetReplica(i)
    DataProcessor.SaveTMDGrid.SaveGrid_optimal_kT(SAVEPATH5+"_"+"{:04d}".format(i)+".dat",Q=2.)
    DataProcessor.SaveTMDGrid.SaveGrid_optimal_kT(SAVEPATH51+"_"+"{:04d}".format(i)+".dat",
                                               Xrange=DataProcessor.SaveTMDGrid.XrangeFF_default,PDF="uTMDFF",h=1,Q=2.)
    DataProcessor.SaveTMDGrid.SaveGrid_optimal_kT(SAVEPATH52+"_"+"{:04d}".format(i)+".dat",
                                               Xrange=DataProcessor.SaveTMDGrid.XrangeFF_default,PDF="uTMDFF",h=2,Q=2.)