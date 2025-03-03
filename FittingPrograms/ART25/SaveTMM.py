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

MAINPATH=ROOT_DIR 

SAVEPATH0="/data/WorkingFiles/TMD/Fit_Notes/ART25/DataForPLOTS/TMM0/"
SAVEPATH2="/data/WorkingFiles/TMD/Fit_Notes/ART25/DataForPLOTS/TMM2/"
#%%
#######################################
#Initialize artemide
#######################################
import harpy

path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/ART25_main.atmde"
#### the gluon term is needed only to compute correction, however, it requiresseparate computation because of resummation
#### os, it is fixed without resummation, and computed separately.
#path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/ART25_gluon.atmde"

harpy.initialize(path_to_constants)

#%%
import DataProcessor.ArtemideReplicaSet

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                        "/data/WorkingFiles/TMD/Fit_Notes/ART25/REPLICAS/ART25_run1.rep")
    
rSet.SetReplica(0)
#%%
## I save only first 300 replicas (to save disk-space)
N_replicas_max=501
## Like in article
Q=20. 

#%%
SAVEPATH_CS="/data/WorkingFiles/TMD/Fit_Notes/ART25/DataForPLOTS/"

bValues=[0.055, 0.105, 0.155, 0.205, 0.255, 0.305, 0.355, 0.405, 0.455, \
0.505, 0.555, 0.605, 0.655, 0.705, 0.755, 0.805, 0.855, 0.905, 0.955, \
1.005, 1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, \
3.8, 4., 4.2, 4.4, 4.6, 4.8, 5.]

with open(SAVEPATH_CS+'CS_Q=2.dat', 'w') as outfile:        
    for i in range(N_replicas_max):      
        rSet.SetReplica(i,part="TMDR")    
        print(">>"+str(i))        
        G0=[str(harpy.get_DNP(b, 2.)) for b in bValues]
        outfile.write(",   ".join(G0)+"\n")
    
#%%

################################# TMM 0 - section #############################

######## quark TMM0
for i in range(N_replicas_max):      
    rSet.SetReplica(i)    
    print(">>"+str(i))
    with open(SAVEPATH0+'G0_TMDPDF_Q=20_x.'+str(i).zfill(4), 'w') as outfile:        
        for i in range(81):                
            x=10**(-0.05*i)
            G0=harpy.get_uTMDPDF_G0(x,Q,1)
            G0str=[str(g) for g in G0]
            outfile.write(str(x)+",   "+",   ".join(G0str)+"\n")
            
for i in range(N_replicas_max):      
    rSet.SetReplica(i)    
    print(">>pi/"+str(i))
    with open(SAVEPATH0+'G0_TMDFF_pi_Q=20_x.'+str(i).zfill(4), 'w') as outfile:        
        for i in range(66):                
            x=10**(-0.02*i)
            G0=harpy.get_uTMDFF_G0(x,Q,1)
            G0str=[str(g) for g in G0]
            outfile.write(str(x)+",   "+",   ".join(G0str)+"\n")
            
for i in range(N_replicas_max):      
    rSet.SetReplica(i)    
    print(">>K/"+str(i))
    with open(SAVEPATH0+'G0_TMDFF_K_Q=20_x.'+str(i).zfill(4), 'w') as outfile:        
        for i in range(66):                
            x=10**(-0.02*i)
            G0=harpy.get_uTMDFF_G0(x,Q,2)
            G0str=[str(g) for g in G0]
            outfile.write(str(x)+",   "+",   ".join(G0str)+"\n")
            

#%%
####### only-gluon TMM0
for i in range(N_replicas_max):      
    rSet.SetReplica(i)    
    print(">>"+str(i))
    with open(SAVEPATH0+'G0_TMDPDFg_Q=20_x.'+str(i).zfill(4), 'w') as outfile:        
        for i in range(81):                
            x=10**(-0.05*i)
            G0=harpy.get_uTMDPDF_G0(x,Q,1)
            G0str=[str(G0[5])]
            outfile.write(str(x)+",   "+",   ".join(G0str)+"\n")
            
#%%
for i in range(N_replicas_max):      
    rSet.SetReplica(i)    
    print(">>pi/"+str(i))
    with open(SAVEPATH0+'G0_TMDFFg_pi_Q=20_x.'+str(i).zfill(4), 'w') as outfile:        
        for i in range(66):                
            x=10**(-0.02*i)
            G0=harpy.get_uTMDFF_G0(x,Q,1)
            G0str=[str(G0[5])]
            outfile.write(str(x)+",   "+",   ".join(G0str)+"\n")
            
for i in range(N_replicas_max):      
    rSet.SetReplica(i)    
    print(">>K/"+str(i))
    with open(SAVEPATH0+'G0_TMDFFg_K_Q=20_x.'+str(i).zfill(4), 'w') as outfile:        
        for i in range(66):                
            x=10**(-0.02*i)
            G0=harpy.get_uTMDFF_G0(x,Q,2)
            G0str=[str(G0[5])]
            outfile.write(str(x)+",   "+",   ".join(G0str)+"\n")
            
#%%
######## PDFs
for i in range(N_replicas_max):      
    rSet.SetReplica(i)    
    print(">>"+str(i))
    with open(SAVEPATH0+'PDF_Q=20_x.'+str(i).zfill(4), 'w') as outfile:        
        for i in range(81):                
            x=10**(-0.05*i)
            G0=harpy.get_uPDF(x,Q,1)
            G0str=[str(g) for g in G0]
            outfile.write(str(x)+",   "+",   ".join(G0str)+"\n")
            
for i in range(N_replicas_max):      
    rSet.SetReplica(i)    
    print(">>pi/"+str(i))
    with open(SAVEPATH0+'FF_pi_Q=20_x.'+str(i).zfill(4), 'w') as outfile:        
        for i in range(66):                
            x=10**(-0.02*i)
            G0=harpy.get_uFF(x,Q,1)
            G0str=[str(g) for g in G0]
            outfile.write(str(x)+",   "+",   ".join(G0str)+"\n")
            
for i in range(N_replicas_max):      
    rSet.SetReplica(i)    
    print(">>K/"+str(i))
    with open(SAVEPATH0+'FF_K_Q=20_x.'+str(i).zfill(4), 'w') as outfile:        
        for i in range(66):                
            x=10**(-0.02*i)
            G0=harpy.get_uFF(x,Q,2)
            G0str=[str(g) for g in G0]
            outfile.write(str(x)+",   "+",   ".join(G0str)+"\n")
            
#%%

W=[harpy.get_uTMDFF_X0(0.3,2.+mu,1)[7] for mu in range(80)]

A=[harpy.get_uTMDFF_ASX0(0.3,2.+mu,1)[7] for mu in range(80)]

#%%            
################################# TMM 2 - section #############################
######## quark TMM2
for i in range(N_replicas_max):      
    rSet.SetReplica(i)    
    print(">>"+str(i))
    with open(SAVEPATH2+'G2_TMDPDF_Q=20_x.'+str(i).zfill(4), 'w') as outfile:        
        for i in range(81):                
            x=10**(-0.05*i)
            G0=harpy.get_uTMDPDF_X0(x,Q,1)-Q**2/2*harpy.get_uTMDPDF_ASX0(x,Q,1)
            G0str=[str(g) for g in G0]
            outfile.write(str(x)+",   "+",   ".join(G0str)+"\n")
            
for i in range(N_replicas_max):      
    rSet.SetReplica(i)    
    print(">>pi/"+str(i))
    with open(SAVEPATH2+'G2_TMDFF_pi_Q=20_x.'+str(i).zfill(4), 'w') as outfile:        
        for i in range(66):                
            x=10**(-0.02*i)
            G0=harpy.get_uTMDFF_X0(x,Q,1)-Q**2/2*harpy.get_uTMDFF_ASX0(x,Q,1)
            G0str=[str(g) for g in G0]
            outfile.write(str(x)+",   "+",   ".join(G0str)+"\n")
            
for i in range(N_replicas_max):      
    rSet.SetReplica(i)    
    print(">>K/"+str(i))
    with open(SAVEPATH2+'G2_TMDFF_K_Q=20_x.'+str(i).zfill(4), 'w') as outfile:        
        for i in range(66):                
            x=10**(-0.02*i)
            G0=harpy.get_uTMDFF_X0(x,Q,2)-Q**2/2*harpy.get_uTMDFF_ASX0(x,Q,2)
            G0str=[str(g) for g in G0]
            outfile.write(str(x)+",   "+",   ".join(G0str)+"\n")