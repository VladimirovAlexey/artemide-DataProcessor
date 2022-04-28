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

useOrder="nnlo"
#useOrder="n3lo"

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
if(useOrder=="nnlo"):
    ### only SIDIS case
    r1Set=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/"+
                                                  "Sivers20_BPV20(nnlo-noDY).rep")
    ### SIDIS+DY case
    r2Set=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/"+
                                                  "Sivers20_BPV20(nnlo).rep")
elif(useOrder=="n3lo"):
    ### only SIDIS case
    r1Set=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/"+
                                                  "Sivers20_BPV20(n3lo-noDY).rep")
    ### SIDIS+DY case
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
## Check various properties of the point
## A=bool  : point is to be included into the plot (delta<deltaMAX)
## B=bool  : was the point in the experiment)
## C=float : delta for the point
## returns (A,B,C)
#########################################
deltaMAX=0.8
def CheckPoint(setName,p):
    namesOfExperiments=['compass08.sivers.pi+.dpt', 'compass08.sivers.pi-.dpt',
                    'compass08.sivers.k+.dpt', 'compass08.sivers.k-.dpt',
                    'compass16.sivers.h+.1<z<2.dpt','compass16.sivers.h-.1<z<2.dpt',
                    'compass16.sivers.h+.z>2.dpt' ,'compass16.sivers.h-.z>2.dpt',
                    'hermes.sivers.pi+.3d','hermes.sivers.pi-.3d',
                    'hermes.sivers.k+.3d','hermes.sivers.k-.3d',
                    'jlab.sivers.pi+','jlab.sivers.pi-','jlab.sivers.k+',
                    'star.sivers.W+.dqT','star.sivers.W-.dqT',
                    'star.sivers.Z',
                    'compass.sivers.piDY.dqT']
    wasIncluded=True
    ### part of orignial cut function, that compute deltas, and check some kinematics
    if p["type"]=="DY":
        deltaTEST=0.3
        delta=p["<qT>"]/p["<Q>"]             
        if(9<p["<Q>"]<11):#UPSILON resonance-bin
            wasIncluded=False
    
    if p["type"]=="SIDIS":   
        deltaTEST=0.3
        delta=p["<pT>"]/p["<z>"]/p["<Q>"]
    
    ### Testing the inclusion
    if((setName in namesOfExperiments)
       and delta<deltaTEST and delta<deltaMAX and wasIncluded):
        wasIncluded=True
    else:
        wasIncluded=False
        
    ### testing inclusionnow
    nowInclude=(delta<deltaMAX)
    
    return nowInclude,wasIncluded,delta

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
        
    ### Compute CFvalue
    r1Set.SetReplica(0)
    CF1=DataProcessor.harpyInterface.ComputeXSec(p,method=methodX)/normX
    r2Set.SetReplica(0)
    CF2=DataProcessor.harpyInterface.ComputeXSec(p,method=methodX)/normX
    
    ### Compute xSec for each replica
    rr1=[]
    for i in range(1,r1Set.numberOfReplicas+1):
        r1Set.SetReplica(i)
        rr1.append(DataProcessor.harpyInterface.ComputeXSec(p,method=methodX)/normX)
    rr2=[]
    for i in range(1,r2Set.numberOfReplicas+1):
        r2Set.SetReplica(i)
        rr2.append(DataProcessor.harpyInterface.ComputeXSec(p,method=methodX)/normX)
        
    ### Compute 68%CI
    CI1=Compute68CI(rr1)
    CI2=Compute68CI(rr2)
    
    if(not CI1[0]<CF1<CI1[1]):
        print("Point ",p["id"]," has CF outside 68CI (noDYset) ->",[CF1,CI1[0],CI1[1]])
    if(not CI2[0]<CF2<CI2[1]):
        print("Point ",p["id"]," has CF outside 68CI (DY+SIDIS) ->",[CF2,CI2[0],CI2[1]])
    
    return [[CF1,CI1[0],CI1[1]],[CF2,CI2[0],CI2[1]]]

#%%
################################################################
## reading data base
#################################################################
import os

listToInclude=['hermes.sivers.pi+.Qint.dx.csv','hermes.sivers.pi+.Q>2.dpt.csv',
 'compass16.sivers.h-.z>2.dx.csv', 'compass16.sivers.h+.z>2.dx.csv', 'compass08.sivers.pi-.dpt.csv',
 'hermes.sivers.k+.Q<2.dpt.csv', 'hermes.sivers.pi+.Q<2.dpt.csv', 'hermes.sivers.k+.Qint.dx.csv',
 'compass.sivers.piDY.dQ.csv', 'compass08.sivers.k+.dx.csv', 'compass16.sivers.h-.z>2.dpt.csv',
 'compass16.sivers.h+.z>2.dpt.csv', 'hermes.sivers.pi-.Qint.dx.csv', 'compass16.sivers.h+.1<z<2.dpt.csv',
 'star.sivers.W-.dy.csv', 'compass16.sivers.h+.1<z<2.dz.csv', 'compass.sivers.piDY.dqT.csv',
 'compass16.sivers.h-.1<z<2.dx.csv', 'compass16.sivers.h+.z>1.dx.csv', 'hermes.sivers.k+.Qint.dpt.csv',
 'star.sivers.W-.dqT.csv', 'star.sivers.Z.csv', 'hermes.sivers.k+.Q>2.dpt.csv',
 'star.sivers.W+.dy.csv', 'hermes.sivers.k-.Qint.dx.csv', 'hermes.sivers.k+.Q>2.dz.csv',
 'compass16.sivers.h-.1<z<2.dpt.csv', 'compass08.sivers.k+.dpt.csv', 'compass08.sivers.pi-.dx.csv',
 'compass08.sivers.k+.dz.csv', 'hermes.sivers.pi+.Qint.dpt.csv', 'jlab.sivers.k+.csv',
 'hermes.sivers.pi-.Qint.dz.csv', 'compass16.sivers.h-.z>1.dx.csv', 'hermes.sivers.pi-.3d.csv',
 'hermes.sivers.k+.Q<2.dz.csv', 'compass08.sivers.k-.dz.csv', 'compass08.sivers.pi-.dz.csv',
 'compass16.sivers.h-.1<z<2.dz.csv', 'hermes.sivers.k-.3d.csv', 'jlab.sivers.pi-.csv',
 'hermes.sivers.pi-.Qint.dpt.csv', 'jlab.sivers.pi+.csv',
 'compass08.sivers.pi+.dz.csv', 'compass16.sivers.h+.z>2.dz.csv', 'compass16.sivers.h-.z>2.dz.csv',
 'compass16.sivers.h+.z>1.dz.csv', 'compass.sivers.piDY.dxF.csv', 'compass08.sivers.k-.dpt.csv',
 'compass08.sivers.k-.dx.csv', 'jlab.sivers.k-.csv', 'compass08.sivers.pi+.dx.csv',
 'hermes.sivers.pi+.3d.csv', 'hermes.sivers.k+.Qint.dz.csv', 'compass16.sivers.h+.1<z<2.dx.csv',
 'compass16.sivers.h-.z>1.dz.csv', 'star.sivers.W+.dqT.csv', 'compass08.sivers.pi+.dpt.csv',
 'hermes.sivers.pi+.Qint.dz.csv', 'hermes.sivers.k-.Qint.dz.csv', 'compass16.sivers.h+.z>1.dpt.csv',
 'hermes.sivers.pi+.Q<2.dz.csv', 'hermes.sivers.k-.Qint.dpt.csv', 'compass16.sivers.h-.z>1.dpt.csv',
 'hermes.sivers.pi+.Q>2.dz.csv', 'hermes.sivers.k+.3d.csv']

path_to_data="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/DataLib/Sivers/"
listOfFiles=[]
for f in os.listdir(path_to_data):
    if( f in listToInclude):
        listOfFiles.append(f)


#%%
######################################
## Compute all for each point in the data set
## returns a list with [[kineamtics, xSec, deltaXSec, theory1, theory1-,theory1+,theory2, theory2-,theory2+]]
######################################
def ProcessDataSet(dataIN):
    res=[]
    for p in dataIN.points:
        isIncluded,wasIncluded,delta=CheckPoint(dataIN.name,p)
        
        if(isIncluded):
            asym=ProcessPoint(p)
        else:
            continue
        
        saveLine=[]
        ## saving kinematics
        if(p["type"]=="SIDIS"):
            saveLine.append(p["<x>"])
            saveLine.append(p["<z>"])
            saveLine.append(p["<Q>"])    
            saveLine.append(p["<pT>"])        
        else:
            saveLine.append(p["<y>"])
            saveLine.append(p["<Q>"])    
            saveLine.append(p["<qT>"])    
        saveLine.append(delta)
        
        saveLine.append(wasIncluded)
        
        ## saving experimental measurement
        ## total uncorrelated uncertanty
        err=numpy.sqrt(numpy.sum([x**2 for x in p["uncorrErr"]]))
        saveLine.append(p["xSec"])
        saveLine.append(err)
        
        ## saving results of computation
        saveLine=saveLine+asym
        res.append(saveLine)
    return res

#%%
####################################################
## Process all data sets
####################################################
path_to_save="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/PlotData/"
for file in listOfFiles:
    loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+file)    
    print("Working with ",file, "   number of points =",loadedData.numberOfPoints)
    
    resultOfProcessing=ProcessDataSet(loadedData)
    
    fileName=loadedData.name+"_"+useOrder+".dat"
    print("saving to file: ",fileName)
    
    f=open(path_to_save+fileName,"w")
    if(loadedData.processType=="SIDIS"):
        f.write("## <x>, <z>, <Q>, <pT>, delta, was in fit?, experimental value, experimental uncertanty, "+
                "SIDIS fit central value, SIDIS fit low 68CI, SIDIS fit up 68CI,"+
                "SIDIS+DY fit central value, SIDIS+DY fit low 68CI, SIDIS+DY fit up 68CI"+"\n")
    else:
        
        f.write("## <y>, <Q>, <qT>, delta, was in fit?, experimental value, experimental uncertanty, "+
                "SIDIS fit central value, SIDIS fit low 68CI, SIDIS fit up 68CI,"+
                "SIDIS+DY fit central value, SIDIS+DY fit low 68CI, SIDIS+DY fit up 68CI"+"\n")
    
    for x in resultOfProcessing:
        f.write(str(x)+"\n")
    
    f.close()