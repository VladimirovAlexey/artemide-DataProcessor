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
## Process the point with x,z,Q,pT(list)
#########################################
def ProcessPoint(x,z,Q,pTlist):
    
    n=len(pTlist)    
    
    r2Set.SetReplica(0)
    
    s=2*12.*0.938 + (0.938)**2.
    
    FUUlist=[]
    FUTlist=[]
    
    sList=numpy.full(n,s)
    xList=numpy.full(n,x)
    zList=numpy.full(n,z)
    QList=numpy.full(n,Q)
    mList=[[0.938,0.139] for i in range(n)]
    pList=[[1,2,2001] for i in range(n)]
    
    for i in range(unSet.numberOfReplicas):
        SetUnReplica(i+1)
        ff=harpy.SIDIS.xSecListBINLESS(pList,sList,pTlist,zList,xList,QList,mList)
        FUUlist.append(ff)        
        
    SetUnReplica(0)
    
    unMean=numpy.mean(FUUlist,axis=0)
    unStd=numpy.std(FUUlist,axis=0)
    
    
    pList=[[1,2,12001] for i in range(n)]
    
    for i in range(r2Set.numberOfReplicas):
        r2Set.SetReplica(i+1)
        ff=harpy.SIDIS.xSecListBINLESS(pList,sList,pTlist,zList,xList,QList,mList)
        FUTlist.append(ff)     
    
    sivMean=numpy.mean(FUTlist,axis=0)    
    siv68=[]
    for i in range(n):
        siv68.append(Compute68CI([r[i] for r in FUTlist]))
            
    res=[]
    for i in range(n):
        res.append(
            [unMean[i],unStd[i],
             sivMean[i],siv68[i][0]-sivMean[i],siv68[i][1]-sivMean[i]
             ])
    
    return res

#%%
#########################################
## Process the point with x,z,Q,pT(list) NO SEA CONTRIBUTION
#########################################
def ProcessPointNoSeaS(x,z,Q,pTlist):
    import copy
    
    n=len(pTlist)    
    
    r2Set.SetReplica(0)
    
    s=2*12.*0.938 + (0.938)**2.
    
    FUTlist=[]
    
    sList=numpy.full(n,s)
    xList=numpy.full(n,x)
    zList=numpy.full(n,z)
    QList=numpy.full(n,Q)
    mList=[[0.938,0.139] for i in range(n)]
    pList=[[1,2,12001] for i in range(n)]
    
    SetUnReplica(0)
    
    for i in range(r2Set.numberOfReplicas):
        rValue=copy.copy(r2Set.GetReplica(i+1))
        rValue[11]=0.
        rValue[13]=0.
        harpy.setNPparameters_SiversTMDPDF(rValue)
        
        ff=harpy.SIDIS.xSecListBINLESS(pList,sList,pTlist,zList,xList,QList,mList)
        FUTlist.append(ff)     
    
    sivMean=numpy.mean(FUTlist,axis=0)    
    siv68=[]
    for i in range(n):
        siv68.append(Compute68CI([r[i] for r in FUTlist]))
            
    res=[]
    for i in range(n):
        res.append(
            [sivMean[i],siv68[i][0]-sivMean[i],siv68[i][1]-sivMean[i]
             ])
    
    return res

#%%
#########################################
## Process the point with x,z,Q,pT(list) NO SEA CONTRIBUTION
#########################################
def ProcessPointNoSea(x,z,Q,pTlist):
    import copy
    n=len(pTlist)    
    
    r2Set.SetReplica(0)
    
    s=2*12.*0.938 + (0.938)**2.
    
    FUTlist=[]
    
    sList=numpy.full(n,s)
    xList=numpy.full(n,x)
    zList=numpy.full(n,z)
    QList=numpy.full(n,Q)
    mList=[[0.938,0.139] for i in range(n)]
    pList=[[1,2,12001] for i in range(n)]
    
    SetUnReplica(0)
    
    for i in range(r2Set.numberOfReplicas):
        rValue=copy.copy(r2Set.GetReplica(i+1))
        rValue[13]=0.
        harpy.setNPparameters_SiversTMDPDF(rValue)
        
        ff=harpy.SIDIS.xSecListBINLESS(pList,sList,pTlist,zList,xList,QList,mList)
        FUTlist.append(ff)     
    
    sivMean=numpy.mean(FUTlist,axis=0)    
    siv68=[]
    for i in range(n):
        siv68.append(Compute68CI([r[i] for r in FUTlist]))
            
    res=[]
    for i in range(n):
        res.append(
            [sivMean[i],siv68[i][0]-sivMean[i],siv68[i][1]-sivMean[i]
             ])
    
    return res


#%% 
######################################
## compute x-sec prefactor by comparing the structure function to the xSec.
######################################
def ComputeXsecPrefactor(x,z,Q,pT):
    
    SetUnReplica(0)
    r2Set.SetReplica(0)
    
    s1=2*12.*0.938 + (0.938)**2.
    s2=2*25.*0.938 + (0.938)**2.
    s3=4*4.*51.
    s4=4*18*275.
        
    sList=[s1,s2,s3,s4]
    xList=[x,x,x,x]
    zList=[z,z,z,z]
    QList=[Q,Q,Q,Q]
    mList=[[0.938,0.139] for i in range(4)]
    pList=[[1,2,2001] for i in range(4)]
    pXList=[[1,1,2001] for i in range(4)]
    
    f1=harpy.SIDIS.xSecListBINLESS(pList,sList,xList,zList,xList,QList,mList)
    f2=harpy.SIDIS.xSecListBINLESS(pXList,sList,xList,zList,xList,QList,mList)
    
    return list(f2/f1)

#%%
######################################
## grid values
######################################
############################################# Harut grid

xValues=[0.02, 0.05, 0.1, 0.2, 0.3, 0.4]

zValues=[0.2, 0.4, 0.6, 0.8]

QValues=[2, 4, 6,  8, 10, 20, 40]
    
pTValues=[0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6]



#%%
######################################
## writing Grid
######################################
totNum=0
for x in xValues:
    for z in zValues:
        for Q in QValues:
            totNum+=1
            # for pT in pTValues:
            #     if pT/z<=Q*0.3: totNum+=1
print("Total number of points to operate: ",totNum)

#%%
######################################
## writing Grid
######################################

import os

filename=SAVEPATH+"grid_light_Harut.csv"

if os.path.exists(filename):
  os.remove(filename)

currentN=-1
for x in xValues:
    for z in zValues:
        line=[]
        for Q in QValues:
            
            currentN+=1            
            
            pTlist=[]            
            for pT in pTValues:
                if pT/z<Q*0.3: pTlist.append(pT)
            
            if(len(pTlist)==0): continue
        
            
            res=ProcessPoint(x,z,Q,pTlist)
            
            res2=ProcessPointNoSea(x,z,Q,pTlist)
            
            res3=ProcessPointNoSeaS(x,z,Q,pTlist)
            
            xSecFactors=ComputeXsecPrefactor(x,z,Q,pTlist[0])
            
            for i in range(len(pTlist)):
                line.append([x,z,Q,pTlist[i]]+res[i]+res2[i]+res3[i]+xSecFactors)
                
            print("Processed: "+"{:6.3f}".format(100*float(currentN/totNum))+" %")
        
        if(line!=[]):
            file=open(filename,"a+")
            for l in line:
                file.write("    ".join(['{:g}'.format(k) for k in l])+" \n")
            file.close()