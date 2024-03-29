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
usePDF="DSSV"
#usePDF="NNPDF"
#%%
#######################################
#Initialize artemide
#######################################
import harpy
path_to_constants=MAINPATH+"/FittingPrograms/WGT22/ConstantsFiles/"

if(useOrder=="nnlo"):
    harpy.initialize(path_to_constants+"const-WGT22_nnlo_"+usePDF)
    
    #### All=0 Case
    harpy.setNPparameters_TMDR([2., 0.0398333])
    harpy.setNPparameters_uTMDPDF([0.185239, 6.22706, 580.946, 2.44166, -2.53161, 0.,  0.0014, 0.442, 4.14])
    harpy.setNPparameters_uTMDFF([0.279443, 0.460015, 0.435955, 0.551302])
    harpy.setNPparameters_wgtTMDPDF([0.2, 1.0])
    
elif(useOrder=="n3lo"):
    harpy.initialize(path_to_constants+"const-WGT22_n3lo_"+usePDF)
    #### All=0 Case n3lo
    harpy.setNPparameters_TMDR([2.0, 0.04843])
    harpy.setNPparameters_uTMDPDF([0.1425, 4.8199, 580.9, 2.3889, -1.0683, 0.0,  0.0014, 0.442, 4.14])
    harpy.setNPparameters_uTMDFF([0.2797, 0.4469, 0.43215, 0.63246])
    harpy.setNPparameters_wgtTMDPDF([0.5, 1.0])
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
if(usePDF=="NNPDF"): harpy.setNPparameters_wgtTMDPDF([0.518,0.414])
elif(usePDF=="DSSV"):harpy.setNPparameters_wgtTMDPDF([0.472, 0.368])
else:raise ValueError("usePDF unknown")

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

#initialValues=(0.472, 0.368)
initialValues=(0.472, 1)

initialErrors=(1.,  1.)
searchLimits=((-5,5),(-10,10))

# True= FIX
parametersToMinimize=(False, True)

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

harpy.setNPparameters_wgtTMDPDF(list(m.values))

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

ll=[]

numOfReplicas=300
REPPATH=MAINPATH+"/FittingPrograms/WGT22/LOGS/"+"l12(d<0.35;Q2>2;NNPDF;rep)-replicas.txt"
for i in range(numOfReplicas):
    print('---------------------------------------------------------------')
    print('------------REPLICA ',i,'/',numOfReplicas,'--------------------')
    print('---------------------------------------------------------------')
    
    
    #repPDF=numpy.random.randint(1, high=100)
    #harpy.sethPDFreplica(repPDF)
    ### to speed up computation I generate  large number of data-pseudo replicas for single PDF
    ### at the next turn I will random sample it (after deleting ugly cases)
    #for j in range(20):
    repRes=MinForReplica()    
    print(repRes)
    ll.append(repRes)
        #f=open(REPPATH,"a+")
        #print('SAVING >>  ',f.name)
        #f.write(str(repRes)+", "+str(repPDF)+"\n")
        #f.close()
#%%
kk=[numpy.abs(lll[2][0]) for lll in ll]
print(numpy.mean(kk))
      
#%%
def chi_2silent(x):
    #startT=time.time()
    harpy.setNPparameters_wgtTMDPDF(x)
    #print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
    
    ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(setALT,method="central")
    
    ccSIDIS2+=PenaltyTerm(x)
    
    cc=(ccSIDIS2)/totalN
    #endT=time.time()
    #print(':->',cc)
    return ccSIDIS2
#%%
from scipy.optimize import fsolve
x0=1
for i in range(10):
    func= lambda t : (chi_2silent([i*0.2,t[0]])/totalN-0.92)
    solution=fsolve(func, x0,xtol=0.001)
    x0=solution[0]
    print("{"+"{:8.4f}, {:6.3f}".format(i*0.2,solution[0])+"},")

x0=0
for i in range(10):
    func= lambda t : (chi_2silent([1.8-i*0.2,t[0]])/totalN-.92)
    solution=fsolve(func, x0,xtol=0.001)
    x0=solution[0]
    print("{"+"{:8.4f}, {:6.3f}".format(1.8-i*0.2,solution[0])+"},")
    
#%%
for i in range(100):
    print([0.01,i*0.03],chi_2silent([0.01,i*0.03])/totalN)