#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 16:40:59 2019

@author: vla18041
"""

#######################################
# importing libraries
#######################################
import os
MAINPATH=os.path.join(os.path.dirname(__file__),"..","..")

import sys
import time
import numpy

sys.path.append(MAINPATH)

import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet

useOrder="n3lo"
#%%
#######################################
#Initialize artemide
#######################################
import harpy
path_to_constants=MAINPATH+"/FittingPrograms/WGT22/ConstantsFiles/"

if(useOrder=="nnlo"):
    harpy.initialize(path_to_constants+"const-WGT22_nnlo")
    
    #### All=0 Case
    harpy.setNPparameters_TMDR([2., 0.0398333])
    harpy.setNPparameters_uTMDPDF([0.185239, 6.22706, 580.946, 2.44166, -2.53161, 0.,  0.0014, 0.442, 4.14])
    harpy.setNPparameters_uTMDFF([0.279443, 0.460015, 0.435955, 0.551302])
    harpy.setNPparameters_wgtTMDPDF([0.2, 1.0])
    
elif(useOrder=="n3lo"):
    harpy.initialize(path_to_constants+"const-WGT22_n3lo")
    #### All=0 Case n3lo
    harpy.setNPparameters_TMDR([2.0, 0.04843])
    harpy.setNPparameters_uTMDPDF([0.1425, 4.8199, 580.9, 2.3889, -1.0683, 0.0,  0.0014, 0.442, 4.14])
    harpy.setNPparameters_uTMDFF([0.2797, 0.4469, 0.43215, 0.63246])
    harpy.setNPparameters_wgtTMDPDF([1.5, 1.0])
#%%
def loadThisData(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/DataLib/wgt/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

#%%
##################Cut function
def cutFunc(p):
    import copy 
    
    if p["type"]=="SIDIS":   
        deltaTEST=0.35
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
    
        
#    return delta<0.5 and p.qT_avarage<80
    return delta<deltaTEST and p["<Q>"]>1.41, p

#%%
### Loading the SIDIS data set
theData=DataProcessor.DataMultiSet.DataMultiSet("ALTset",loadThisData([
                      'hermes3D.ALT.pi+','hermes3D.ALT.pi-',
                      'hermes3D.ALT.k+','hermes3D.ALT.k-',
                      'compass16.ALT.h+.2<z.dpt','compass16.ALT.h-.2<z.dpt',
                      'compass16.ALT.h+.2<z.dz','compass16.ALT.h-.2<z.dz',
                      'compass16.ALT.h+.2<z.dx','compass16.ALT.h-.2<z.dx'#,
                      #'JLab6.ALT.pi+','JLab6.ALT.pi-'
                      ]))

setALT=theData.CutData(cutFunc) 

print('Loaded ', setALT.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setALT.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setALT.sets])

#SaveToLog('Loaded '+ str(setDY.numberOfSets) + ' data sets with '+str(sum([i.numberOfPoints for i in setDY.sets])) + ' points. \n'
#+'Loaded experiments are '+str([i.name for i in setDY.sets]))
#%%
harpy.setNPparameters_wgtTMDPDF([0.518,0.414])

DataProcessor.harpyInterface.PrintChi2Table(setALT,printDecomposedChi2=True,method="central")

#%%
def PenaltyTerm(x):
    return 0
    #if(x[0]<0.5): return 10.*(numpy.exp((0.5-x[0])**4)-1)
    #else:  return 0.
#%%
totalN=setALT.numberOfPoints

def chi_2(x):
    startT=time.time()
    harpy.setNPparameters_wgtTMDPDF(x)
    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
    
    ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(setALT,method="central")
    
    ccSIDIS2+=PenaltyTerm(x)
    
    cc=(ccSIDIS2)/totalN
    endT=time.time()
    print(':->',cc,'       t=',endT-startT)
    return ccSIDIS2


#%%

from iminuit import Minuit

initialValues=(1.5,1)

initialErrors=(0.1,  0.1)
searchLimits=((-5,5),(-10,10))

# True= FIX
parametersToMinimize=(False, False)

m = Minuit(chi_2, initialValues)

m.errors=initialErrors
m.limits=searchLimits
m.fixed=parametersToMinimize
m.errordef=1

print(m.params)

m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1
m.migrad()

print(m.params)

harpy.setNPparameters(list(m.values))

DataProcessor.harpyInterface.PrintChi2Table(setALT,printDecomposedChi2=True,method="central")

#m.minos()

#print(m.params)

#harpy.setNPparameters(list(m.values))

#DataProcessor.harpyInterface.PrintChi2Table(setALT,printDecomposedChi2=True)

#%%

import DataProcessor.harpyInterface
DataProcessor.harpyInterface.PrintPerPointContribution(setALT,method="central",output="id,<Q>,<x>,del")


#%%

def MinForReplica():
    
    
    def repchi_2(x):        
        startT=time.time()
        harpy.setNPparameters_wgtTMDPDF(x)
        #print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
            
        ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(repDataSIDIS,method="central")
        
        ccSIDIS2+=PenaltyTerm(x)
        
        cc=(ccSIDIS2)/totalNnew
        
        endT=time.time()
        #print(':->',cc,'       t=',endT-startT)
        return ccSIDIS2
    
    repDataSIDIS=setALT.GenerateReplica()
    totalNnew=repDataSIDIS.numberOfPoints    
    
    localM = Minuit(repchi_2, initialValues)
    localM.errors=initialErrors
    localM.limits=searchLimits
    localM.fixed=parametersToMinimize
    localM.errordef=1    
    localM.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
    localM.strategy=1

    localM.migrad()
    
    chi2Central=chi_2(list(localM.values))
    
    return [localM.fval,chi2Central,list(localM.values)]

#%%
#
# Generate pseudo data and minimise   100 times
#
numOfReplicas=1000
REPPATH=MAINPATH+"/FittingPrograms/WGT22/LOGS/"+"l12(d<0.35;Q2>2)-replicas.txt"
for i in range(numOfReplicas):
    print('---------------------------------------------------------------')
    print('------------REPLICA ',i,'/',numOfReplicas,'--------------------')
    print('---------------------------------------------------------------')
    repRes=MinForReplica()
    print(repRes)
    f=open(REPPATH,"a+")
    print('SAVING >>  ',f.name)
    f.write(str(repRes)+"\n")
    f.close()