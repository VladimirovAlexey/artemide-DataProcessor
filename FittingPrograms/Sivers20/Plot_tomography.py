#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 16:40:59 2019

@author: vla18041
"""

#######################################
# importing libraries
#######################################

import sys
import numpy
#sys.path.append("/home/m/Github/artemide-DataProcessor")
sys.path.append("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor")
import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet
import DataProcessor.ArtemideReplicaSet

#MAINPATH="/home/m/Github/artemide-DataProcessor/"
MAINPATH="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"

#useOrder="nnlo"
useOrder="n3lo"

#%%
#######################################
#Initialize artemide
#######################################
import harpy
path_to_constants=MAINPATH+"FittingPrograms/Sivers20/Constants-files/"
#harpy.initialize(path_to_constants+"const-Sivers20_lo")

# harpy.initialize(path_to_constants+"const-Sivers20_nnlo_piK")

if(useOrder=="nnlo"):
    harpy.initialize(path_to_constants+"const-Sivers20_nnlo")
    
    #### All=0 Case
    harpy.setNPparameters_TMDR([2., 0.0398333])
    harpy.setNPparameters_uTMDPDF([0.185239, 6.22706, 580.946, 2.44166, -2.53161, 0.,  0.0014, 0.442, 4.14])
    harpy.setNPparameters_uTMDFF([0.279443, 0.460015, 0.435955, 0.551302])
    
    unSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(
    "/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/Replicas/DY+SIDIS/SV19_nnlo_all=0.rep")
    
    # # #### All=0 case piK
    # harpy.setNPparameters_TMDR([2., 0.0394095])
    # harpy.setNPparameters_uTMDPDF([0.180718, 4.38119, 426.208, 2.22347, -0.0646396, 0., 0.17, 0.48, 2.15])
    # harpy.setNPparameters_uTMDFF([0.293548, 0.462093, 0.442867, 0.590596, 0.427915, 0.462578, 0.304421,1.18113])
elif(useOrder=="n3lo"):
    harpy.initialize(path_to_constants+"const-Sivers20_n3lo")
    #### All=0 Case n3lo
    harpy.setNPparameters_TMDR([2., 0.0442327])
    harpy.setNPparameters_uTMDPDF([0.17975, 3.9081, 453.883, 2.07199, 1.60774, 0,  0.0014, 0.442, 4.14])
    harpy.setNPparameters_uTMDFF([0.26994, 0.456091, 0.423312, 0.615092])
    
    unSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(
    "/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/Replicas/DY+SIDIS/SV19_n3lo_all=0.rep")

############### Loading the replica distributions
if(useOrder=="nnlo"):
    
    rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/"+
                                                  "Sivers20_BPV20(nnlo).rep")### SIDIS+DY case
elif(useOrder=="n3lo"):
    rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/"+
                                                  "Sivers20_BPV20(n3lo).rep")### SIDIS+DY case

#%%
### Set unpolarized TMD according to TMD set (adds constant pion row)
def SetUnTMD(n):
    unSet.SetReplica(n,part="TMDR")    
    unSet.SetReplica(n,part="uTMDFF")
    rr=unSet.GetReplica(n,part="uTMDPDF")
    harpy.setNPparameters_uTMDPDF(rr[:-1]+[0.0014, 0.442, 4.14])

#%%
########################################
# Function that return the value of Tomographic function
# f1-px/M f1T
# mu is scale. if present return F(x,b,mu,mu^2)
# n is number of replica, is n<0, rnadom replicas selectted
########################################
def getTomography(x,px,py,mu=-1.,n=-1):
    M_proton=0.932
    p=numpy.sqrt(px**2+py**2)
    if(n<0):
        SetUnTMD(numpy.random.randint(1,unSet.numberOfReplicas+1))
        rSet.SetReplica(numpy.random.randint(1,rSet.numberOfReplicas+1))
    else:
        SetUnTMD(n)
        rSet.SetReplica(n)
    tmd=numpy.array(harpy.get_uTMDPDF_kT(x,p,1,mu=mu))-px*M_proton/p*numpy.array(harpy.get_SiversTMDPDF_kT(x,p,1,mu=mu))
    return tmd
#%%
path_to_save="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/PlotData/Tomography/"

xValues=[0.5,0.1,0.05,0.01]

for x in xValues:

    nameADD="_2GeV_rnd_x="+str(x)+".dat"
    muIn=2.
    
    dQuark=[]
    uQuark=[]
    sQuark=[]
    dbarQuark=[]
    ubarQuark=[]
    for kx in numpy.arange(-1.,1.+0.01,0.01):
        for ky in numpy.arange(-1.,1.+0.01,0.01):        
            tmd=getTomography(x, kx, ky ,n=-1,mu=muIn)
            uQuark.append([kx,ky,tmd[2+5]])
            dQuark.append([kx,ky,tmd[1+5]])
            sQuark.append([kx,ky,tmd[3+5]])
            dbarQuark.append([kx,ky,tmd[-1+5]])
            ubarQuark.append([kx,ky,tmd[-2+5]])
    
    f=open(path_to_save+"u"+nameADD,"w")
    for x in uQuark:
        f.write(str(x)+"\n")
    f.close()
    
    f=open(path_to_save+"d"+nameADD,"w")
    for x in dQuark:
        f.write(str(x)+"\n")
    f.close()
    
    f=open(path_to_save+"s"+nameADD,"w")
    for x in sQuark:
        f.write(str(x)+"\n")
    f.close()
    
    f=open(path_to_save+"dbar"+nameADD,"w")
    for x in dbarQuark:
        f.write(str(x)+"\n")
    f.close()
    
    f=open(path_to_save+"ubar"+nameADD,"w")
    for x in ubarQuark:
        f.write(str(x)+"\n")
    f.close()
