#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 16:42:03 2023

@author: alexey
"""

#######################################
# importing libraries
#######################################
import os
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"

ATMDE_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide/"
HARPY_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide/harpy/"

replicaFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/replicaDY.dat"
logFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/LOG.log"

import sys
import numpy
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)

#%%
import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet

MAINPATH=ROOT_DIR 

#%%
#######################################
#Initialize artemide
#######################################
import harpy

PDFtoUSE="MSHT"
#PDFtoUSE="NNPDF"

path_to_constants=MAINPATH+"FittingPrograms/ART23/ConstantsFiles/DYonly/ART23_"+PDFtoUSE+"_N4LL.atmde"
#path_to_constants=MAINPATH+"FittingPrograms/ART23/ConstantsFiles/DYonly/ART23_JAM_N4LL"

harpy.initialize(path_to_constants)
#%%
if(PDFtoUSE=="MSHT"):
    rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/data/WorkingFiles/TMD/Fit_Notes/ART23/REPS/ART23_run2.rep")
elif(PDFtoUSE=="NNPDF"):
    rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/data/WorkingFiles/TMD/Fit_Notes/ART23/REPS/ART23_NNPDF.rep")
else:
    raise Exception("no no no")
rSet.SetReplica(-1)

rSet.SetReplica(0)
#%%
### create the data set

import DataProcessor.Point
import DataProcessor.DataSet

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('LHC-sense',"DY")
DataCurrent.comment=" "
DataCurrent.reference=" "

DataCurrent.isNormalized=True
proc_current=[1,3,1,1]
s_current=13000.**2
Q_current=[66.,116.]
y_current=[-2.5,2.5]
incCut=False
cutParam=[27.,27.,-2.5,2.5]
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

for i in range(44,90):
    
    Q_current=[10.**(i/30.),10.**((i+1)/30.)]
    for j in range(30):
        qT_current=[j+0.,j+1.]
    
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=qT_current
        p["thFactor"]=1./(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
        p["Q"]=Q_current
        p["y"]=y_current
        p["includeCuts"]=incCut
        p["cutParams"]=cutParam
        p["xSec"]=10.
        DataCurrent.AddPoint(p)

print("Done.  ")


DataWithOutCut=DataCurrent

#%%
### create the data set

import DataProcessor.Point
import DataProcessor.DataSet

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('LHC-sense',"DY")
DataCurrent.comment=" "
DataCurrent.reference=" "

DataCurrent.isNormalized=True
proc_current=[1,3,1,1]
s_current=13000.**2
Q_current=[66.,116.]
y_current=[-2.5,2.5]
incCut=True
cutParam=[27.,27.,-2.5,2.5]
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

for i in range(44,90):
    
    Q_current=[10.**(i/30.),10.**((i+1)/30.)]
    for j in range(30):
        qT_current=[j+0.,j+1.]
    
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=qT_current
        p["thFactor"]=1./(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
        p["Q"]=Q_current
        p["y"]=y_current
        p["includeCuts"]=incCut
        p["cutParams"]=cutParam
        p["xSec"]=10.
        DataCurrent.AddPoint(p)

print("Done.  ")


DataWithCut=DataCurrent
#%%
#%%
######################
### This function compute correlaiton coefficient for the vector of replicas for the given set
### It test 6 functions (given by F1,F2,F3,F4,F5,F6)
### Returns [corr1,corr2,corr3,corr4,corr5,corr6, xSec, stdX]
### where corr is the correlation number ,stdX is (theory) uncertanty of the xSec
#####################
def CorrelationParameters(inputSet):
    
    vectorO=[]
    vectorF1=[]
    vectorF2=[]
    vectorF3=[]
    vectorF4=[]    
    vectorF5=[]
    vectorF6=[]      
    
    n1=0
    n2=0
    ## compute vectors
    #for n in range(1,rSet.numberOfReplicas+1):
    for n in range(300):
        RND=numpy.random.randint(1,high=rSet.numberOfReplicas)
        rSet.SetReplica(RND)
        NPparams=rSet.GetReplica(RND)
        ### vector of cross-section
        XX1=DataProcessor.harpyInterface.ComputeXSec(inputSet)    
        ### vector of test function
        F1=[]
        F2=[]
        F3=[]
        F4=[]
        F5=[]
        F6=[]
        for p in inputSet.points:               
            
            F1.append(NPparams[1])# CS KERNEL 1
            F2.append(NPparams[2])# CS kernel 2
            
            F3.append(NPparams[4])            
            F4.append(NPparams[5])
            
            x_current=p["<Q>"]/numpy.sqrt(p["s"])            
            TMD=harpy.get_uTMDPDF(x_current, 0.5, 1)            
            F5.append(TMD[3]*TMD[7]*4./9+TMD[4]*TMD[6]*1./9)
            
            x_current=p["<Q>"]/numpy.sqrt(p["s"])            
            TMD=harpy.get_uTMDPDF(x_current, 1., 1)            
            F6.append(TMD[3]*TMD[7]*4./9+TMD[4]*TMD[6]*1./9)
            
        vectorO.append(XX1)
        vectorF1.append(F1)
        vectorF2.append(F2)
        vectorF3.append(F3)
        vectorF4.append(F4)
        vectorF5.append(F5)
        vectorF6.append(F6)
        
        n1+=1
        n2+=1
        if(n1>(rSet.numberOfReplicas/10.)):
            n1=0
            print(n2/(rSet.numberOfReplicas+1)*100.,' %')
            
    
        
    vectorO=numpy.array(vectorO)
    vectorF1=numpy.array(vectorF1)
    vectorF2=numpy.array(vectorF2)
    vectorF3=numpy.array(vectorF3)
    vectorF4=numpy.array(vectorF4)
    vectorF5=numpy.array(vectorF5)
    vectorF6=numpy.array(vectorF6)
    
    avOF1=numpy.mean(numpy.multiply(vectorO,vectorF1),axis=0)
    avOF2=numpy.mean(numpy.multiply(vectorO,vectorF2),axis=0)
    avOF3=numpy.mean(numpy.multiply(vectorO,vectorF3),axis=0)
    avOF4=numpy.mean(numpy.multiply(vectorO,vectorF4),axis=0)
    avOF5=numpy.mean(numpy.multiply(vectorO,vectorF5),axis=0)
    avOF6=numpy.mean(numpy.multiply(vectorO,vectorF6),axis=0)
    
    avO=numpy.mean(vectorO,axis=0)
    avF1=numpy.mean(vectorF1,axis=0)
    avF2=numpy.mean(vectorF2,axis=0)
    avF3=numpy.mean(vectorF3,axis=0)
    avF4=numpy.mean(vectorF4,axis=0)
    avF5=numpy.mean(vectorF5,axis=0)
    avF6=numpy.mean(vectorF6,axis=0)
    
    
    stdO=numpy.std(vectorO,axis=0)    
    stdF1=numpy.std(vectorF1,axis=0)
    stdF2=numpy.std(vectorF2,axis=0)
    stdF3=numpy.std(vectorF3,axis=0)
    stdF4=numpy.std(vectorF4,axis=0)
    stdF5=numpy.std(vectorF5,axis=0)
    stdF6=numpy.std(vectorF6,axis=0)
    
    rho1=numpy.divide(avOF1-numpy.multiply(avO,avF1),numpy.multiply(stdO,stdF1))
    rho2=numpy.divide(avOF2-numpy.multiply(avO,avF2),numpy.multiply(stdO,stdF2))
    rho3=numpy.divide(avOF3-numpy.multiply(avO,avF3),numpy.multiply(stdO,stdF3))
    rho4=numpy.divide(avOF4-numpy.multiply(avO,avF4),numpy.multiply(stdO,stdF4))
    rho5=numpy.divide(avOF5-numpy.multiply(avO,avF5),numpy.multiply(stdO,stdF5))
    rho6=numpy.divide(avOF6-numpy.multiply(avO,avF6),numpy.multiply(stdO,stdF6))
    
    ### relative experimental error
    deltaOth=numpy.divide(stdO,avO)
    
    return numpy.transpose([rho1,rho2,rho3,rho4,rho5,rho6,deltaOth,avO])

#%%
### Compute the correlation matrix and save to file
def computeAndSaveSense(setIN,path):
    rho=CorrelationParameters(setIN)
    f=open(path,"w")        
    for i in range(setIN.numberOfPoints):    
        p=setIN.points[i]
        line=("{:f} , {:f} , {:f} , {:f} ,".format(p["Q"][0],p["Q"][1],p["qT"][0],p["qT"][1])
              +"{:f} , {:f} , {:f} , {:f} , {:f} , {:f} , {:f} , {:f}".format(
                  rho[i][0],rho[i][1],rho[i][2],rho[i][3],rho[i][4],rho[i][5],rho[i][6],rho[i][7])
              +"\n")
        f.write(line)
    
    f.close()
    
#%%
savePATH="/data/WorkingFiles/TMD/Fit_Notes/ART23/Sensetivity_for_LHC/"
computeAndSaveSense(DataWithCut,savePATH+"withCUT.csv")