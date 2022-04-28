#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 19:51:14 2020

@author: alexey
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
    harpy.initialize(path_to_constants+"const-Sivers20_plot_nnlo")
    
elif(useOrder=="n3lo"):
    harpy.initialize(path_to_constants+"const-Sivers20_plot_n3lo")    

############### Loading the replica distributions
if(useOrder=="nnlo"):
    
    rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/"+
                                                  "Sivers20_BPV20(nnlo).rep")### SIDIS+DY case
    harpy.setNPparameters_TMDR([2., 0.0398333])
elif(useOrder=="n3lo"):
    rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/"+
                                                  #"Sivers20_model9case1(noDY-n3lo).rep")### only SIDIS case
                                                  "Sivers20_BPV20(n3lo).rep")### SIDIS+DY case
    harpy.setNPparameters_TMDR([2., 0.0442327])

#%%
#########################################################
## Resample data with N/2
#########################################################
def Resample(dd):
    #return numpy.random.choice(dd,size=int(numpy.floor(len(dd)/2)))
    return dd[numpy.random.choice(dd.shape[0], size=int(numpy.floor(len(dd)/2)))]

#########################################################
## Determine Mean, Mode, 68%CI by resampling 
#########################################################
alpha=68
def Compute68CI(dd):    
    lowers=[]
    uppers=[]    
    for i in range(1500):
        sample=Resample(numpy.array(dd))
        lowers.append(numpy.percentile(sample,(100-alpha)/2))
        uppers.append(numpy.percentile(sample,100-(100-alpha)/2))
    
    return [numpy.mean(lowers),numpy.mean(uppers)]



#%%
########################################
# Function that return the value of Sivers function distributed over replicas
## mu is scale. if present return F(x,b,mu,mu^2)
## The first element is CF value
########################################
def getSiversRow(x,b,mu=-1.):
    result=[]
    for j in range(rSet.numberOfReplicas+1):
        rSet.SetReplica(j)
        result.append(numpy.array(harpy.get_SiversTMDPDF(x,b,1,mu=mu)))
    return numpy.array(result)

#%%
########################################
# Function that return the value of Sivers function distributed over replicas
# mu is scale. if present return F(x,b,mu,mu^2)
## The first element is CF value
########################################
def getSiversRowP(x,p,mu=-1.):
    result=[]
    M2_proton=0.932**2
    for j in range(rSet.numberOfReplicas+1):
        rSet.SetReplica(j)
        result.append(numpy.array(harpy.get_SiversTMDPDF_kT(x,p,1,mu=mu))*M2_proton/p)
    return numpy.array(result)
#%%
################################################
## Computes the value of Sivers function + 68%CI band
## returns [u,d,s,sea] where u=[CF,CImin,CImax]
###############################################
def ComputeSivers(x,b,mu=-1.):
    tmd=getSiversRow(x,b,mu=mu)
    
    f=1
    CI=Compute68CI([t[f+5] for t in tmd[1:]])
    dQuark=[tmd[0][f+5],CI[0],CI[1]]
    
    f=2
    CI=Compute68CI([t[f+5] for t in tmd[1:]])
    uQuark=[tmd[0][f+5],CI[0],CI[1]]
    
    f=3
    CI=Compute68CI([t[f+5] for t in tmd[1:]])
    sQuark=[tmd[0][f+5],CI[0],CI[1]]
    
    f=-1
    CI=Compute68CI([t[f+5] for t in tmd[1:]])
    seaQuark=[tmd[0][f+5],CI[0],CI[1]]
    return [uQuark,dQuark,sQuark,seaQuark]

def ComputeSiversP(x,pT,mu=-1.):
    tmd=getSiversRowP(x,pT,mu=mu)
    
    f=1
    CI=Compute68CI([t[f+5] for t in tmd[1:]])
    dQuark=[tmd[0][f+5],CI[0],CI[1]]
    
    f=2
    CI=Compute68CI([t[f+5] for t in tmd[1:]])
    uQuark=[tmd[0][f+5],CI[0],CI[1]]
    
    f=3
    CI=Compute68CI([t[f+5] for t in tmd[1:]])
    sQuark=[tmd[0][f+5],CI[0],CI[1]]
    
    f=-1
    CI=Compute68CI([t[f+5] for t in tmd[1:]])
    seaQuark=[tmd[0][f+5],CI[0],CI[1]]
    return [uQuark,dQuark,sQuark,seaQuark]

#%%
# #########################Compute cases and save into files
# xValues=[0.001, 0.0011, 0.0013, 0.0014, 0.0016, 0.0018, 0.002, 0.0022, \
# 0.0025, 0.0028, 0.0032, 0.0035, 0.004, 0.0045, 0.005, 0.0056, 0.0063, \
# 0.0071, 0.0079, 0.0089, 0.01, 0.011, 0.013, 0.014, 0.016, 0.018, \
# 0.02, 0.022, 0.025, 0.028, 0.032, 0.035, 0.04, 0.045, 0.05, 0.056, \
# 0.063, 0.071, 0.079, 0.089, 0.1, 0.11, 0.13, 0.14, 0.16, 0.18, 0.2, \
# 0.22, 0.25, 0.28, 0.32, 0.35, 0.4, 0.45, 0.5, 0.56, 0.63, 0.71, 0.79, \
# 0.89, 1.]
# bValues=[0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, \
# 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.1, 1.2, 1.3, 1.4, 1.5, \
# 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, \
# 3., 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4., 4.5, 5., 5.5, \
# 6., 6.5, 7., 7.5, 8.]

# # bValues=[0.001, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, \
# # 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, \
# # 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, \
# # 1.85, 1.9, 1.95, 2.]



# path_to_save="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/PlotData/Sivers3D_"+useOrder+"/"
# nameADD="_b_2GeV_"+useOrder+".dat"
# muIn=2.

# ############# For EIC impact
# # bValues=[0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 1.]
# # rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(
# # "/home/vla18041/LinkData2/WorkingFiles/TMD/YR_Studies/Sivers/REPS/Sivers20_HB_opt8.rep")
# # meanReplica=rSet.GetReplica(0)
# # rSet.SetReplica()
# # path_to_save="/home/vla18041/LinkData2/WorkingFiles/TMD/YR_Studies/Sivers/3Dplots/"
# # nameADD="_b_optimal_EIC_opt8.dat"
# # muIn=-1.



# dQuark=[]
# uQuark=[]
# sQuark=[]
# seaQuark=[]
# k=0
# kmax=len(xValues)*len(bValues)
# for x in xValues:
#     for b in bValues:        
#         print(k,"/",kmax)
#         tmd=ComputeSivers(x, b,mu=muIn)
#         uQuark.append([x,b]+tmd[0])
#         dQuark.append([x,b]+tmd[1])
#         sQuark.append([x,b]+tmd[2])
#         seaQuark.append([x,b]+tmd[3])
#         k+=1

# f=open(path_to_save+"u"+nameADD,"w")
# for x in uQuark:
#     f.write(str(x)+"\n")
# f.close()

# f=open(path_to_save+"d"+nameADD,"w")
# for x in dQuark:
#     f.write(str(x)+"\n")
# f.close()

# f=open(path_to_save+"s"+nameADD,"w")
# for x in sQuark:
#     f.write(str(x)+"\n")
# f.close()

# f=open(path_to_save+"sea"+nameADD,"w")
# for x in seaQuark:
#     f.write(str(x)+"\n")
# f.close()


#%%
#########################Compute cases and save into files
xValues=[0.001, 0.0011, 0.0013, 0.0014, 0.0016, 0.0018, 0.002, 0.0022, \
0.0025, 0.0028, 0.0032, 0.0035, 0.004, 0.0045, 0.005, 0.0056, 0.0063, \
0.0071, 0.0079, 0.0089, 0.01, 0.011, 0.013, 0.014, 0.016, 0.018, \
0.02, 0.022, 0.025, 0.028, 0.032, 0.035, 0.04, 0.045, 0.05, 0.056, \
0.063, 0.071, 0.079, 0.089, 0.1, 0.11, 0.13, 0.14, 0.16, 0.18, 0.2, \
0.22, 0.25, 0.28, 0.32, 0.35, 0.4, 0.45, 0.5, 0.56, 0.63, 0.71, 0.79, \
0.89, 1.]

bValues=[0.001, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, \
0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, \
1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, \
1.85, 1.9, 1.95, 2.]

path_to_save="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/PlotData/Sivers3D_"+useOrder+"/"
nameADD="_kT_2GeV_"+useOrder+".dat"
muIn=2.

dQuark=[]
uQuark=[]
sQuark=[]
seaQuark=[]
k=0
kmax=len(xValues)*len(bValues)
for x in xValues:
    for b in bValues:        
        print(k,"/",kmax)
        tmd=ComputeSiversP(x, b,mu=muIn)
        uQuark.append([x,b]+tmd[0])
        dQuark.append([x,b]+tmd[1])
        sQuark.append([x,b]+tmd[2])
        seaQuark.append([x,b]+tmd[3])
        k+=1

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

f=open(path_to_save+"sea"+nameADD,"w")
for x in seaQuark:
    f.write(str(x)+"\n")
f.close()
