#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  2 12:18:24 2021

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

r2Set.SetReplica(0)

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
## Process the point for DY (without Err)
#########################################
def ProcessPointDYnoERR():
    
    n=len(pTlist)    
    
    r2Set.SetReplica(0)
    SetUnReplica(0)
            
    sList=numpy.full(n,s)
    yList=[[yMin,yMax] for i in range(n)]    
    QList=[[QMin,QMax] for i in range(n)]    
    pList=[process for i in range(n)]
    iC=[False for i in range(n)]
    cP=[[0.,0.,0.,0.] for i in range(n)]
    
    
    dX=[(QList[i][1]**2-QList[i][0]**2)*(yList[i][1]-yList[i][0])*(pTlist[i][1]**2-pTlist[i][0]**2) for i in range(n)]
    
    xx=harpy.DY.xSecList(pList, sList, pTlist, QList, yList, iC, cP)
    
    return numpy.array(xx)/numpy.array(dX)

#%%
#########################################
## Process the point for DY (without Err)
#########################################
def ProcessPointSIDISnoERR():
    
    n=len(pTlist)    
    
    r2Set.SetReplica(0)
    SetUnReplica(0)
            
    sList=numpy.full(n,s)
    xList=numpy.full(n,x)
    zList=numpy.full(n,z)
    QList=numpy.full(n,Q)
    pT0list=numpy.mean(pTlist,axis=1)
    # xList=[[xMin,xMax] for i in range(n)]    
    # zList=[[zMin,zMax] for i in range(n)]    
    # QList=[[QMin,QMax] for i in range(n)]    
    
    pList=[process for i in range(n)]
    mList=[[0.938,0.139] for i in range(n)] 
    
    xx=harpy.SIDIS.xSecListBINLESS(pList, sList, pT0list, zList, xList, QList, mList)
    
    return numpy.array(xx)

#%%
################################################################
## Plot X-sec for given process in DY
################################################################

################################## CDF run1
# process=[1,2,6]
# s=(1800.)**2
# yMin=-0.1
# yMax=0.1
# QMin=91.-20.
# QMax=91.+20.
# pTlist=[[0.5*i-0.25,0.5*i+0.25] for i in range(1,120)]

################################## COMPASS piDY
process=[2,2,101]
s=357.
yMin=0.35-0.1
yMax=0.35+0.1
QMin=4.5
QMax=5.5
pTlist=[[0.05*i-0.025,0.05*i+0.025] for i in range(1,100)]

xSec=ProcessPointDYnoERR()

print([[numpy.mean(pTlist[i]),xSec[i]] for i in range(len(xSec))])

#%%
################################################################
## Plot X-sec for given process in SIDIS
################################################################

# ################################## HERMES - like
process=[1, 1, 2001]
s=52.
x=0.1
z=0.1
Q=3.
pTlist=[[0.02*i-0.01,0.02*i+0.01] for i in range(1,30)]

################################## JLAB6 - like
# process=[1, 1, 2001]
# s=23.
# x=0.4
# z=0.2
# Q=3.
# pTlist=[[0.02*i-0.01,0.02*i+0.01] for i in range(1,50)]

xSec=ProcessPointSIDISnoERR()

print([[numpy.mean(pTlist[i]),xSec[i]] for i in range(len(xSec))])