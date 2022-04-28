#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: vla18041
"""

#%%
#######################################################################
# Global parameter of a run
#######################################################################

#PDFinUse="HERA20"
#PDFinUse="NNPDF31"
PDFinUse="CT18"


## automatic name generation for the run
runName="model2.2_"+PDFinUse+"_correlation2_LHCsense_"
if (runName[-1]=="_"): runName=runName[0:-1]

print(" RUN: "+runName)

#%%
#######################################
# Paths
#######################################
PathToHarpy="/home/vla18041/LinkData2/arTeMiDe_Repository/artemide-ForPDF/harpy"
PathToDataProcessor="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"
PathToLog=PathToDataProcessor+"FittingPrograms/PDF-TMD/LOGS/"
PathToConstantsFile=PathToDataProcessor+"/FittingPrograms/PDF-TMD/Constants-files/const-"+PDFinUse+"_NNLO_7p"

import sys
sys.path.remove('/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/harpy')
sys.path.append(PathToHarpy)
sys.path.append(PathToDataProcessor)

#%%
#######################################
# Output paths
#######################################
outFile=PathToLog+runName+".txt"

#%%
#######################################
# importing libraries
#######################################
import numpy

import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet

#%%
#######################################
#Initialize artemide
#######################################
import harpy

harpy.initialize(PathToConstantsFile)
initializationArray=[2.0340, 0.0299, 0.1, 2., 0.1, 2., 0.1, 1., 10000.]
harpy.setNPparameters(initializationArray)

#%%
#######################################
# Data cut function
#######################################
def cutFunc(p):    
    delta=p["<qT>"]/p["<Q>"]
    
    if(p["id"][0] == "E"):
        delta=p["<qT>"]/p["Q"][1]    
        
    if(p["id"][0:4] == "E605"):
        if(p["Q"][0]==10.5):#UPSILON resonance-bin
            return False , p
    elif(p["id"][0:4] == "E772"):
        if(p["Q"][0]<10):#these bins seems broken
            return False , p
    elif(p["id"][0:4] == "E615"):
        if(9<p["<Q>"]<11.2):#UPSILON resonance-bin
            return False , p
    elif(p["id"][0:4] == "E228"):
        if(9<p["<Q>"]<11):#UPSILON resonance-bin
            return False , p
    else:
        if(9<p["<Q>"]<11):#UPSILON resonance-bin
            return False , p
    
    return (delta<0.25) , p

#%%
#######################################
# Loading the data set
#######################################

import DataProcessor.DataSet
dataCollection=[]
loadedData=DataProcessor.DataSet.LoadCSV("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/OtherPrograms/LHC_impact_Study/LHC13TeV-Empty.csv")
dataCollection.append(loadedData)

theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",dataCollection)

setDY=theData.CutData(cutFunc) 

print('Loaded ', setDY.numberOfSets, 'data sets with ', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setDY.sets])

#%%

if(PDFinUse=="HERA20"):
    #initialValues=(2.000,  0.033, 0.230, 5.609, 0.252, 8.021, 570.223, 0.000, 0.000) #model 1.0 HERA
    #initialValues=(2.000,  0.033, 0.230, 5.609, 0.252, 8.021, 0.252, 8.021, 570.223) #model 2.0 HERA
    initialValues=(2.000, 0.034, 0.145, 10.393, 0.333, 0.000, 0.366, 10.261, 677.434) #model 2.2 HERA    
if(PDFinUse=="NNPDF31"):
    #initialValues=(2.000,  0.029, 0.345, 2.587, 0.152, 7.561, 232.544, 0. , 0.)      #model 1.0 NNPDF31
    #initialValues=(2.000,  0.029, 0.345, 2.587, 0.152, 7.561, 0.152, 7.561, 232.544) #model 2.0 NNPDF31
    initialValues=(2.000, 0.030, 0.188, 5.542, 0.200, 4.375, 0.486, 0.009, 232.793) #model 2.2 NNPDF31
if(PDFinUse=="CT18"):
    #initialValues=(2.000,  0.039, 0.161, 7.904, 0.212, 5.301, 700.642, 0. , 0.)      #model 1.0 CT18
    #initialValues=(2.000,  0.039, 0.161, 7.904, 0.212, 5.301, 0.212, 5.301, 700.642) #model 2.0 CT18
    initialValues=(2.000, 0.042, 0.094, 12.534, 0.293, 0.004, 0.003, 16.568, 819.267) #model 2.2 CT18


harpy.setNPparameters(list(initialValues))

#%%
######################
### This function compute correlaiton coefficient for the parameter pN, varied range pMin,pMax
### Returns [corr1, deltaO]
### where corr is the correlation number ,stdX is (theory) uncertanty of the xSec, deltaO is experimental unvertanty of xSec
#####################
def CorrelationParameters(pN,deltap):
    import copy
    Npoints=50
    vectorO=[]
    vectorF1=[]
    
    n1=0
    n2=0
    ## compute xSec for each of parameters 
    for n in range(Npoints):
        #param=pMin+n/(Npoints-1)*(pMin-pMax)
        
        param=numpy.abs(numpy.random.normal(initialValues[pN],deltap))
        
        NPparams=copy.copy(list(initialValues))
        NPparams[pN]=param
        harpy.setNPparameters(NPparams)
        
        ### vector of cross-section
        XX1=DataProcessor.harpyInterface.ComputeXSec(setDY)    
        ### vector of test function
        F1=[]
        for p in setDY.points:            
            F1.append(param)
            
        vectorO.append(XX1)
        vectorF1.append(F1)
        
        n1+=1
        n2+=1
        if(n1>(Npoints/10.)):
            n1=0
            print(n2/(Npoints+1)*100.,' %')
            
    
        
    vectorO=numpy.array(vectorO)
    vectorF1=numpy.array(vectorF1)
    
    ### <xSec * p>
    avOF1=numpy.mean(numpy.multiply(vectorO,vectorF1),axis=0)
    ### <xSec>
    avO=numpy.mean(vectorO,axis=0)
    ### <p>
    avF1=numpy.mean(vectorF1,axis=0)
    
    ### delta xSec
    stdO=numpy.std(vectorO,axis=0)    
    ### delta p
    stdF1=numpy.std(vectorF1,axis=0)
    
    ### correlation (<xSec * p>- <xSec><p>)/deltaXsec/deltap
    rho1=numpy.divide(avOF1-numpy.multiply(avO,avF1),numpy.multiply(stdO,stdF1))
        
    return numpy.transpose([rho1,stdO,avO])

#%%
#################################################
# run for each parameters
#################################################
p0=CorrelationParameters(1,0.02)

p1=CorrelationParameters(2,0.1)
p2=CorrelationParameters(3,2.)

p3=CorrelationParameters(4,0.1)
p4=CorrelationParameters(5,0.5)

p5=CorrelationParameters(6,0.2)
p6=CorrelationParameters(7,6.)

p7=CorrelationParameters(8,300)

#%%
#################################################
# save to file
#################################################
f=open(outFile,"w")
print('SAVING >>  ',f.name)
for i in range(setDY.numberOfPoints):    
    ppp=setDY.points[i]    
    f.write(str([ppp["id"],ppp["s"],ppp["Q"],ppp["y"],ppp["qT"],ppp["xSec"],numpy.sqrt(numpy.sum(numpy.array(ppp["uncorrErr"])**2)),
                 list(p0[i]),list(p1[i]),list(p2[i]),list(p3[i]),list(p4[i]),list(p5[i]),list(p6[i]),list(p7[i])])+"\n")
f.close()