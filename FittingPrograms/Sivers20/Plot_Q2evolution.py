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
## Process the point computing cross-section and uncertanty bands for it
## It returns [ A , B] where A for r1Set, B for r2Set
## A~B=[CF value, lowValue, highValue]
#########################################
def ProcessPoint(p):
    import copy
        
    #### methods for DY and SIDIS different
    if p["type"]=="SIDIS":
        methodX="central"        
    elif p["type"]=="DY":
        methodX="default"        
    else:
        print("Are you crazy?")
    
    ### Compute normalization
    pNew=copy.deepcopy(p)    
    pNew["process"]=pNew["weightProcess"]
    normX=DataProcessor.harpyInterface.ComputeXSec(pNew,method=methodX)
    
    #### This is because star measures AN
    if p["id"][0:4]=="star":
        normX=-normX
        
    r2Set.SetReplica(0)
    CF2=DataProcessor.harpyInterface.ComputeXSec(p,method=methodX)/normX
    
    return [CF2]

#%%
######################################
## making a  SIDIS set
######################################


############################## COMPASS08  (Q=1.78,  0.0001 +- 0.0115)
DataCurrentSiv=DataProcessor.DataSet.DataSet('evolutionTEST',"SIDIS")
DataCurrentSiv.comment=" "
DataCurrentSiv.reference="  "
proc_current=[1,1,12011]
proc_denominator=[1,1,2011]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=False
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts
xx0=0.043
xx=[0.03, 1.0]
zz0=0.35
zz=[0.2, 1.0]
pp0=0.1547
pp=[0.1, 0.2]      

############################## HERMES  (Q=1.59,  0.0402 +- 0.0211)
# DataCurrentSiv=DataProcessor.DataSet.DataSet('evolutionTEST',"SIDIS")
# DataCurrentSiv.comment=" "
# DataCurrentSiv.reference="  "
# proc_current=[1,1,12001]
# proc_denominator=[1,1,2001]
# s_current=52.
# includeCuts=False
# cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts
# xx0=0.116
# xx=[0.098, 0.138]
# zz0=0.322
# zz=[0.28, 0.37]
# pp0=0.143
# pp=[0.0, 0.23]  

for i in range(60):    
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current    
    p1["s"]=s_current
    p1["<pT>"]=pp0
    p1["pT"]=pp
    p1["<x>"]=xx0
    p1["x"]=xx
    p1["<z>"]=zz0
    p1["z"]=zz
    if(i<18):
        p1["<Q>"]=1.5+i*0.5
    elif(i<28):
        p1["<Q>"]=10.+(i-17)
    else:
        p1["<Q>"]=1.5+(i-27)*4.0
    p1["Q"]=[p1["<Q>"]-1.,p1["<Q>"]+1.]
    p1["xSec"]=1.
    p1["M_target"]=0.938
    p1["M_product"]=0.139
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(0.1)
    p1["weightProcess"]=proc_denominator
    DataCurrentSiv.AddPoint(p1)

#%%
######################################
## making a DY set
######################################

############################## STAR  (Q=1.59,  0.0402 +- 0.0211)
DataCurrentSiv=DataProcessor.DataSet.DataSet('star',"DY")
DataCurrentSiv.comment=" "
DataCurrentSiv.reference="  "
proc_current=[1,1,10007]
proc_denominator=[1,1,7]
s_constant=80.3/numpy.sqrt(250000.) ## Q/sqrt{s}
includeCuts=False
cutParameters=[25.0, 25.0, -1.0, 1.0] #y, W^2 cuts
yy0=0.0
yy=[-1., 1.]
qT0=6.2
qq=[5.5, 7.0]

for i in range(40):    
    # makeup a point
    p1=DataProcessor.Point.CreateDYPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current        
    p1["<qT>"]=qT0
    p1["qT"]=qq
    p1["<y>"]=yy0
    p1["y"]=yy
    p1["<Q>"]=12.+i*4.
    p1["Q"]=[p1["<Q>"]-2.,p1["<Q>"]+2.]
    p1["s"]=(p1["<Q>"]/s_constant)**2
    p1["xSec"]=1.
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(0.1)
    p1["weightProcess"]=proc_denominator
    DataCurrentSiv.AddPoint(p1)

#%%
######################################
## Compute all for each point in the data set
## returns a list with [[kineamtics, xSec, deltaXSec, theory1, theory1-,theory1+,theory2, theory2-,theory2+]]
######################################
def ProcessDataSet(dataIN):
    res=[]
    for p in dataIN.points:
        asym=ProcessPoint(p)
        
        saveLine=[]
        ## saving kinematics
        if(p["type"]=="SIDIS"):
            #saveLine.append(p["<x>"])
            #saveLine.append(p["<z>"])
            saveLine.append(p["<Q>"])    
            #saveLine.append(p["<pT>"])        
        else:
            #saveLine.append(p["<y>"])
            saveLine.append(p["<Q>"])    
            #saveLine.append(p["<qT>"])    
    
        
        ## saving results of computation
        saveLine=saveLine+asym
        res.append(saveLine)
    return res

#%%
rrrr=ProcessDataSet(DataCurrentSiv)
print(rrrr)