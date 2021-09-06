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

#### If true fSIDIS=+fDY, (wrong)
#### if false fSIDIS=-fDY (correct)
useWrongSign=False

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
    

harpy.setNPparameters_SiversTMDPDF([5.2, 0.,0.,0.,0., -0.6, 15.9, 0.5, -0.2, 21.6, -0.5, -0.1, 0.4, -1.1]) 
#%%
### read the list of files and return the list of DataSets
def loadThisData(listOfNames):    
    import DataProcessor.DataSet
    
    #path_to_data="/home/m/Github/artemide-DataProcessor/DataLib/Sivers/"
    path_to_data="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/DataLib/Sivers/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

#%%
##################Cut function
def cutFunc(p):
    import copy
    
    if p["type"]=="DY":
        deltaTEST=0.3
        delta=p["<qT>"]/p["<Q>"]        
        
        
        if(9<p["<Q>"]<11):#UPSILON resonance-bin
            return False , p
    
    if p["type"]=="SIDIS":   
        deltaTEST=0.3
        delta=p["<pT>"]/p["<z>"]/p["<Q>"]        
    
    
    if delta<deltaTEST:
        pNew=copy.deepcopy(p)    
        pNew["process"]=pNew["weightProcess"]
        if p["type"]=="SIDIS":
            normX=DataProcessor.harpyInterface.ComputeXSec(pNew,method="central")        
        elif p["type"]=="DY":
            normX=DataProcessor.harpyInterface.ComputeXSec(pNew)        
        else:
            print("Are you crazy?")
        p["thFactor"]=p["thFactor"]/normX        
    
    #### This is because star measures AN
    if p["id"][0:4]=="star":
        p["thFactor"]=-p["thFactor"]        
    
    # ##### test sign change
    if(useWrongSign):
        if p["type"]=="DY":
            p["thFactor"]=-p["thFactor"]        
        
#    return delta<0.5 and p.qT_avarage<80
    return delta<deltaTEST, p

#%%
### Loading the data set
# theData=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisData([
#                     'compass.sivers.pi+.dpt', 'compass.sivers.pi-.dpt',
#                     'compass.sivers.k+.dpt', 'compass.sivers.k-.dpt',                  
#                     'hermes.sivers.pi+.Qint.dpt','hermes.sivers.k+.Qint.dpt',
#                     'hermes.sivers.pi-.Qint.dpt','hermes.sivers.k-.Qint.dpt',
#                     'hermes.sivers.pi+.Q<2.dpt','hermes.sivers.k+.Q<2.dpt',
#                     'hermes.sivers.pi+.Q>2.dpt','hermes.sivers.k+.Q>2.dpt',                  
#                     'hermes.sivers.pi+.3d','hermes.sivers.pi-.3d',
#                     'hermes.sivers.k+.3d','hermes.sivers.k-.3d',
#                     'compass16.sivers.h+.1<z<2.dpt',
#                     'compass16.sivers.h-.1<z<2.dpt',
#                     'compass16.sivers.h+.z>1.dpt' ,
#                     'compass16.sivers.h-.z>1.dpt' ,
#                     'compass16.sivers.h+.z>2.dpt' ,
#                     'compass16.sivers.h-.z>2.dpt' ,
#                     "jlab.sivers.pi+","jlab.sivers.pi-",
#                     "jlab.sivers.k+"]))
    

theData=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisData([
                    'compass08.sivers.pi+.dpt', 'compass08.sivers.pi-.dpt',
                    'compass08.sivers.k+.dpt', 'compass08.sivers.k-.dpt',
                    'compass16.sivers.h+.1<z<2.dpt','compass16.sivers.h-.1<z<2.dpt',
                    'compass16.sivers.h+.z>2.dpt' ,'compass16.sivers.h-.z>2.dpt',
                    'hermes.sivers.pi+.3d','hermes.sivers.pi-.3d',
                    'hermes.sivers.k+.3d','hermes.sivers.k-.3d',
                    'jlab.sivers.pi+','jlab.sivers.pi-','jlab.sivers.k+'
                    ]))

setSIDIS=theData.CutData(cutFunc) 

theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisData([
                    'star.sivers.W+.dqT','star.sivers.W-.dqT',
                    'star.sivers.W+.dy','star.sivers.W-.dy',
                    'star.sivers.Z',
                    'compass.sivers.piDY.dqT'
                    #,'compass.sivers.piDY.dQ','compass.sivers.piDY.dxF'
                    ]))

setDY=theData.CutData(cutFunc) 

print('Loaded (SIDIS)', setSIDIS.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setSIDIS.sets]), 'points.')
print('Loaded SIDIS experiments are', [i.name for i in setSIDIS.sets])

print('Loaded (DY)', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded DY experiments are', [i.name for i in setDY.sets])

print('Total number of points:',setSIDIS.numberOfPoints+setDY.numberOfPoints)

#%%
############### Just to check that the model is correct
rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/"+
                                                  "Sivers20_BPV20(n3lo).rep")
                                                  # "Sivers20_model9case1(noDY-n3lo).rep")

rSet.SetReplica()

DataProcessor.harpyInterface.PrintChi2Table(setSIDIS,method="central",printSysShift=False)

DataProcessor.harpyInterface.PrintChi2Table(setDY)

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
        
    ### Compute CFvalue
    rSet.SetReplica(0)
    CF1=DataProcessor.harpyInterface.ComputeXSec(p,method=methodX)/normX
    
    ### Compute xSec for each replica
    rr1=[]
    for i in range(1,rSet.numberOfReplicas+1):
        rSet.SetReplica(i)
        rr1.append(DataProcessor.harpyInterface.ComputeXSec(p,method=methodX)/normX)
        
    ### Compute 68%CI
    CI1=Compute68CI(rr1)
    
    if(not CI1[0]<CF1<CI1[1]):
        print("Point ",p["id"]," has CF outside 68CI (noDYset) ->",[CF1,CI1[0],CI1[1]])
    
    return [CF1,CI1[0],CI1[1]]
#%%
##################################
##writing to file
#################################
path="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Predictions/RHIC_Sivers/"

def makeFile(line,filename):
    file=open(path+filename,"w+")
    file.write("# Prediction by artemide for AN in RHIC kinematics, Based on BSV20n3lo. \n")
    file.write("# "+line+" \n")
    file.write("# yMin, yMax, pTmin, pTMax, AN(central), -delta, +delta \n")
    file.close()

def saveLine(line,filename):
    file=open(path+filename,"a+")
    file.write(",   ".join(['{:g}'.format(k) for k in line])+" \n")
    file.close()
#%%
#####################
## Original data
####################

#### output has the form [ymin, yMax, ptmin, ptMax, cental, +delta, -delta]
filename="original_binning_qT_W+.dat"
makeFile("W+ qT-differential original binning",filename)
for p in setDY.sets[0].points:
    rr=ProcessPoint(p)
    saveLine([p["y"][0],p["y"][1],p["qT"][0],p["qT"][1],rr[0],rr[1]-rr[0],rr[2]-rr[0]],filename)
 
filename="original_binning_qT_W-.dat"
makeFile("W- qT-differential original binning",filename)
for p in setDY.sets[1].points:
    rr=ProcessPoint(p)
    saveLine([p["y"][0],p["y"][1],p["qT"][0],p["qT"][1],rr[0],rr[1]-rr[0],rr[2]-rr[0]],filename)
    
filename="original_binning_y_W+.dat"
makeFile("W+ y-differential original binning",filename)
for p in setDY.sets[2].points:
    rr=ProcessPoint(p)
    saveLine([p["y"][0],p["y"][1],p["qT"][0],p["qT"][1],rr[0],rr[1]-rr[0],rr[2]-rr[0]],filename)
    
filename="original_binning_y_W-.dat"
makeFile("W- y-differential original binning",filename)
for p in setDY.sets[3].points:
    rr=ProcessPoint(p)
    saveLine([p["y"][0],p["y"][1],p["qT"][0],p["qT"][1],rr[0],rr[1]-rr[0],rr[2]-rr[0]],filename)
    
filename="original_binning_Z.dat"
makeFile("Z original binning",filename)
for p in setDY.sets[4].points:
    rr=ProcessPoint(p)
    saveLine([p["y"][0],p["y"][1],p["qT"][0],p["qT"][1],rr[0],rr[1]-rr[0],rr[2]-rr[0]],filename)

#
########################################
## Fine binning
########################################

import copy

yBins=[-1+i*0.1 for i in range(21)]
qTBins=[0.+i*0.25 for i in range(61)]

filename="fine_binning_qT_W+.dat"
makeFile("W+ qT-differential finer binning",filename)
for i in range(len(qTBins)-1):
    p=copy.deepcopy(setDY.sets[0].points[0])  
    p["qT"]=[qTBins[i],qTBins[i+1]]
    rr=ProcessPoint(p)
    saveLine([p["y"][0],p["y"][1],p["qT"][0],p["qT"][1],rr[0],rr[1]-rr[0],rr[2]-rr[0]],filename)
    
filename="fine_binning_qT_W-.dat"
makeFile("W- qT-differential finer binning",filename)
for i in range(len(qTBins)-1):
    p=copy.deepcopy(setDY.sets[1].points[0])  
    p["qT"]=[qTBins[i],qTBins[i+1]]
    rr=ProcessPoint(p)
    saveLine([p["y"][0],p["y"][1],p["qT"][0],p["qT"][1],rr[0],rr[1]-rr[0],rr[2]-rr[0]],filename)
    
filename="fine_binning_y_W+.dat"
makeFile("W+ y-differential finer binning",filename)
for i in range(len(yBins)-1):
    p=copy.deepcopy(setDY.sets[2].points[0])  
    p["y"]=[yBins[i],yBins[i+1]]
    rr=ProcessPoint(p)
    saveLine([p["y"][0],p["y"][1],p["qT"][0],p["qT"][1],rr[0],rr[1]-rr[0],rr[2]-rr[0]],filename)
    
filename="fine_binning_y_W-.dat"
makeFile("W- y-differential finer binning",filename)
for i in range(len(yBins)-1):
    p=copy.deepcopy(setDY.sets[3].points[0])  
    p["y"]=[yBins[i],yBins[i+1]]
    rr=ProcessPoint(p)
    saveLine([p["y"][0],p["y"][1],p["qT"][0],p["qT"][1],rr[0],rr[1]-rr[0],rr[2]-rr[0]],filename)

filename="fine_binning_qT_Z.dat"
makeFile("Z qT-differential finer binning",filename)
for i in range(len(qTBins)-1):
    p=copy.deepcopy(setDY.sets[4].points[0])  
    p["qT"]=[qTBins[i],qTBins[i+1]]
    rr=ProcessPoint(p)
    saveLine([p["y"][0],p["y"][1],p["qT"][0],p["qT"][1],rr[0],rr[1]-rr[0],rr[2]-rr[0]],filename)

filename="fine_binning_y_Z.dat"
makeFile("Z y-differential finer binning",filename)
for i in range(len(yBins)-1):
    p=copy.deepcopy(setDY.sets[4].points[0])  
    p["y"]=[yBins[i],yBins[i+1]]
    rr=ProcessPoint(p)
    saveLine([p["y"][0],p["y"][1],p["qT"][0],p["qT"][1],rr[0],rr[1]-rr[0],rr[2]-rr[0]],filename)