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
import time
import numpy

sys.path.append("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor")
import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet
import DataProcessor.ArtemideReplicaSet

MAINPATH="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"
SAVEPATH="/misc/data2/braun/vla18041/WorkingFiles/TMD/Fit_Notes/Predictions/JLab_grid/"

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
    unSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(
    "/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/Replicas/DY+SIDIS/SV19_nnlo_all=0.rep")
    unSet.SetReplica(0,part="TMDR")    
    unSet.SetReplica(0,part="uTMDFF")
    rr=unSet.GetReplica(0,part="uTMDPDF")
    harpy.setNPparameters_uTMDPDF(rr[:-1]+[0.0014, 0.442, 4.14])
    
    # # #### All=0 case piK
    # harpy.setNPparameters_TMDR([2., 0.0394095])
    # harpy.setNPparameters_uTMDPDF([0.180718, 4.38119, 426.208, 2.22347, -0.0646396, 0., 0.17, 0.48, 2.15])
    # harpy.setNPparameters_uTMDFF([0.293548, 0.462093, 0.442867, 0.590596, 0.427915, 0.462578, 0.304421,1.18113])
elif(useOrder=="n3lo"):
    harpy.initialize(path_to_constants+"const-Sivers20_n3lo")
    
    #### All=0 Case n3lo
    unSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(
    "/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/Replicas/DY+SIDIS/SV19_n3lo_all=0.rep")
    unSet.SetReplica(0,part="TMDR")    
    unSet.SetReplica(0,part="uTMDFF")
    rr=unSet.GetReplica(0,part="uTMDPDF")
    harpy.setNPparameters_uTMDPDF(rr[:-1]+[0.0014, 0.442, 4.14])
    
harpy.setNPparameters_SiversTMDPDF([5.2, 0.,0.,0.,0., -0.6, 15.9, 0.5, -0.2, 21.6, -0.5, -0.1, 0.4, -1.1]) 

#%%
############### Loading the replica distributions
r2Set=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/"+
                                                  "Sivers20_BPV20(n3lo).rep")

def SetUnReplica(n):
    unSet.SetReplica(n,part="TMDR")    
    unSet.SetReplica(n,part="uTMDFF")
    rr=unSet.GetReplica(n,part="uTMDPDF")
    harpy.setNPparameters_uTMDPDF(rr[:-1]+[0.0014, 0.442, 4.14])

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
#########################################
## Process the point with x,z,Q,pT
## It returns [ A , B] where A for FUU and B FUT
#########################################
def ProcessPointCentral(x,z,Q,pT):
    
        
    r2Set.SetReplica(0)
    
    ### no sea option
    # harpy.setNPparameters_SiversTMDPDF([
    #     0.5362, 5.21724, 202.749, 0.0, 0.0, -0.01661, -0.35169, -3.90804, 0.37263, -0.70153, 8.97236, 0., 2.46253, 0.])
    
    s=2*12.*0.938 + (0.938)**2.
    
    ### pion+
    FUU=harpy.SIDIS.xSecListBINLESS([[1,2,2001]],[s],[pT],[z],[x],[Q],[[0.938,0.139]])
    FUT=harpy.SIDIS.xSecListBINLESS([[1,2,12001]],[s],[pT],[z],[x],[Q],[[0.938,0.139]])
    
    ### pion-
    # FUU=harpy.SIDIS.xSecListBINLESS([[1,2,2021]],[s],[pT],[z],[x],[Q],[[0.938,0.139]])
    # FUT=harpy.SIDIS.xSecListBINLESS([[1,2,12021]],[s],[pT],[z],[x],[Q],[[0.938,0.139]])
    
    #if(FUU[0]<0): print("FUU=",FUU[0]," ",[x,z,Q,pT]," delta=",pT/z/Q)
    
    return [FUU[0],FUT[0]]

def ProcessPoint(x,z,Q,pT):
    
        
    r2Set.SetReplica(0)
    
    s=2*12.*0.938 + (0.938)**2.
    
    FUUlist=[]
    FUTlist=[]
    
    
    for i in range(unSet.numberOfReplicas):
        SetUnReplica(i+1)
        ff=harpy.SIDIS.xSecListBINLESS([[1,2,2001]],[s],[pT],[z],[x],[Q],[[0.938,0.139]])        
        FUUlist.append(ff[0])        
    
    SetUnReplica(0)
    
    
    for i in range(r2Set.numberOfReplicas):
        r2Set.SetReplica(i+1)
        ff=harpy.SIDIS.xSecListBINLESS([[1,2,12001]],[s],[pT],[z],[x],[Q],[[0.938,0.139]])
        FUTlist.append(ff[0])     
    
    fT=numpy.mean(FUTlist)
    fS=Compute68CI(FUTlist)
    
    return [numpy.mean(FUUlist),numpy.std(FUUlist),fT,fS[0]-fT,fS[1]-fT]

#%%
######################################
## grid values
######################################
############################################# Larger grid
# xValues=[]
# for i in numpy.arange(-2.,-0.02,0.05):
#     xValues.append(10.**i)

# if(xValues[-1]<0.91): xValues.append(0.91)

# zValues=[]
# for i in numpy.arange(0.1,0.9,0.025):
#     zValues.append(i)

# if(zValues[-1]<0.9): zValues.append(0.9)

# QValues=[]
# for i in numpy.arange(0,3,0.02):
#     QValues.append(10.**i)
    
# pTValues=[]
# for i in numpy.arange(0,3,0.05):
#     pTValues.append(i)

############################################# Smaller grid

# xValues=[]
# for i in numpy.arange(-2.,-0.02,0.1):
#     xValues.append(10.**i)

# if(xValues[-1]<0.91): xValues.append(0.91)

# zValues=[]
# for i in numpy.arange(0.1,0.9,0.05):
#     zValues.append(i)

# if(zValues[-1]<0.9): zValues.append(0.9)

# QValues=[]
# for i in numpy.arange(0,3,0.1):
#     QValues.append(10.**i)
    
# pTValues=[]
# for i in numpy.arange(0,3,0.2):
#     pTValues.append(i)

############################################# Harut grid

xValues=[]
for i in numpy.arange(1,91):
    xValues.append(0.01+i*0.01)

#if(xValues[-1]<0.91): xValues.append(0.91)

zValues=[]
for i in numpy.arange(2,25):
    zValues.append(0.04*i)

#if(zValues[-1]<0.9): zValues.append(0.9)

QValues=[]
for i in numpy.arange(1,40):
    QValues.append(numpy.sqrt(0.5+i*0.5))
for i in numpy.arange(1,51):
    QValues.append(numpy.sqrt(20.+i*10.))
    
pTValues=[]
for i in numpy.arange(1,61):
    pTValues.append(0.05*i)



#%%
######################################
## writing Grid
######################################
totNum=0
for x in xValues:
    for z in zValues:
        for Q in QValues:
            for pT in pTValues:
                if pT/z<=Q*0.3: totNum+=1
print("Total number of points to operate: ",totNum)

#%%
######################################
## writing Grid
######################################
currentN=0
for x in xValues:
    for z in zValues:
        line=[]
        for Q in QValues:
            for pT in pTValues:
                if pT/z>Q*0.3: continue
            
                currentN+=1
                res=ProcessPointCentral(x,z,Q,pT)
                line.append([x,z,Q,pT]+res)
        
        print("Processed: "+"{:6.3f}".format(100*float(currentN)/totNum)+" %")
        file=open(SAVEPATH+"grid_central_Harut_pim_noSEA.csv","a+")
        for l in line:
            file.write('{:g} , {:g} , {:g} , {:g} , {:g} , {:g} \n'.format(*l))
        file.close()

# currentN=0
# for x in xValues:
#     for z in zValues:
#         line=[]
#         for Q in QValues:
#             for pT in pTValues:
#                 if pT/z>Q*0.3: continue
            
#                 currentN+=1
#                 res=ProcessPoint(x,z,Q,pT)
#                 line.append([x,z,Q,pT]+res)
        
#         print("Processed: "+"{:6.3f}".format(100*float(currentN)/totNum)+" %")
#         file=open(SAVEPATH+"grid_full_Harut.csv","a+")
#         for l in line:
#             file.write('{:g} , {:g} , {:g} , {:g} , {:g} , {:g} , {:g} , {:g} , {:g} \n'.format(*l))
#         file.close()