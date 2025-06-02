#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 16 13:54:42 2025

@author: alexey
"""

#######################################
# importing libraries
#######################################
import os
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"

SNOWFLAKE_DIR = "/data/arTeMiDe_Repository/artemide_snowflake/harpy/"

import sys
import numpy
sys.path.remove('/data/arTeMiDe_Repository/artemide/harpy')
sys.path.append(ROOT_DIR)
sys.path.append(SNOWFLAKE_DIR)

#%%
import DataProcessor.harpyInterface
import DataProcessor.snowInterface_N2
import DataProcessor.DataMultiSet
import harpy

#%%
#######################################
#Initialize snowflake
#######################################
path_to_INI=ROOT_DIR+"FittingPrograms/Tw3_FIT/INI/TEST.ini"
harpy.initialize_snowflake(path_to_INI)

NP_par=numpy.zeros(18)+0.2
harpy.setNPparameters_tw3(NP_par)
harpy.UpdateTables(1.0, 105.0)

#%%
#######################################
#Initialize artemide
#######################################

import DataProcessor.ArtemideReplicaSet

path_to_constants=ROOT_DIR+"FittingPrograms/Tw3_FIT/INI/ART25_sivers.atmde"


harpy.initialize(path_to_constants)

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                        "/data/WorkingFiles/TMD/Fit_Notes/ART25/REPLICAS/ART25_run1.rep")
    
rSet.SetReplica(0)

#%%
### read the list of files and return the list of DataSets
def loadThisDataD2(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=ROOT_DIR+"DataLib/D2_moment/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        
    return dataCollection

#%%
### read the list of files and return the list of DataSets
def loadThisDataG2(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=ROOT_DIR+"DataLib/G2/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        
    return dataCollection

#%%
### read the list of files and return the list of DataSets
def loadThisDataSivers(listOfNames):    
    import DataProcessor.DataSet
    
    #path_to_data="/home/m/Github/artemide-DataProcessor/DataLib/Sivers/"
    path_to_data="/data/arTeMiDe_Repository/DataProcessor/DataLib/Sivers/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

#%%
##################Cut function
def cutFunc(p):
    if p["type"]=="G2":
        if p["<Q>"]<numpy.sqrt(2.):
            return False, p
    
    return True, p

#%%
#### If true fSIDIS=+fDY, (wrong)
#### if false fSIDIS=-fDY (correct)
useWrongSign=False
##################Cut function
def cutFunc_Sivers(p):
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
    #if p["id"][0:4]=="star":
    #    p["thFactor"]=-p["thFactor"]        
    #if p["type"]=="DY":
    #    p["thFactor"]=-p["thFactor"]        
    
    # ##### test sign change
    if(useWrongSign):
        if p["type"]=="DY":
            p["thFactor"]=-p["thFactor"]        
        
#    return delta<0.5 and p.qT_avarage<80
    return delta<deltaTEST, p

#%%
### Loading the D2 data set
theData=DataProcessor.DataMultiSet.DataMultiSet("D2set",loadThisDataD2([
    "E143_d2","E154_d2","E155-1999_d2","E155_d2",
    "HallA-2016_d2","HERMES_d2","SANE_d2",
    "RQCD_d2_ud"
    ]))

setD2=theData.CutData(cutFunc) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setD2.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setD2.sets]), 'points.') 

#%%
### Loading the G2 data set
theData=DataProcessor.DataMultiSet.DataMultiSet("G2set",loadThisDataG2([
    #"E142.n", "E143.p", "E143.d","E143.n", 
    "E154.n",
    "E155-29.p","E155-32.p","E155-38.p",
    #"E155-29.d","E155-32.d","E155-38.d",
    #"SMC.p",
    "HERMES",
    "HallA-2004.n","HallA-2016-4.He3","HallA-2016-5.He3",
    ]))

setG2=theData.CutData(cutFunc) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setG2.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setG2.sets]), 'points.') 

#%%
### Loading the data set for Sivers
theData=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisDataSivers([
                    'compass08.sivers.pi+.dpt', 'compass08.sivers.pi-.dpt',
                    'compass08.sivers.k+.dpt', 'compass08.sivers.k-.dpt',
                    'compass16.sivers.h+.1<z<2.dpt','compass16.sivers.h-.1<z<2.dpt',
                    'compass16.sivers.h+.z>2.dpt' ,'compass16.sivers.h-.z>2.dpt',
                    'hermes.sivers.pi+.3d','hermes.sivers.pi-.3d',
                    'hermes.sivers.k+.3d','hermes.sivers.k-.3d',
                    'jlab.sivers.pi+','jlab.sivers.pi-','jlab.sivers.k+'
                    ]))

setSivers=theData.CutData(cutFunc_Sivers) 


print('Loaded (SIDIS)', setSivers.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setSivers.sets]), 'points.')
#print('Loaded SIDIS experiments are', [i.name for i in setSivers.sets])


#%%
# harpy.setNPparameters_tw3([1.8207, 1.0, 
#                            -0.004614, 0.213014, 0.0, 0.0, 
#                            -0.001321, -0.231287, 0.0, 0.0, 
#                            0.0, 0.0, 0.0, 0.0, 
#                            0.0, 0.0, 0.0, 0.0])

# harpy.setNPparameters_SiversTMDPDF([0.5,0.0,0.0,0.0,0.0])

# harpy.UpdateTables(1.0, 105.0)
#%%
harpy.setNPparameters_tw3([2.392, 1.0, 
                           0.023151, 0.198417, -0.678277, -0.557989, 
                           -0.029809, -0.219044, 2.7968, 2.7323, 
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

harpy.UpdateTables(1.0, 105.0)

harpy.setNPparameters_SiversTMDPDF([0.003785,0.0,0.0,0.0,0.0])

# #%%
# #### PLOT d2
# Qplot=1*numpy.array(range(40))+1
# dataP=harpy.D2List(Qplot, Qplot*0+11)
# dataN=harpy.D2List(Qplot, Qplot*0+2)
# dataD=harpy.D2List(Qplot, Qplot*0+3)

#%%
DataProcessor.snowInterface_N2.PrintChi2Table(setD2,printDecomposedChi2=False)
DataProcessor.snowInterface_N2.PrintChi2Table(setG2,printDecomposedChi2=False)

DataProcessor.harpyInterface.PrintChi2Table(setSivers,method="central",printSysShift=False)
#%%
#######################################
# Minimisation
#######################################
import time

def chi2(x):
    startT=time.time()
    harpy.setNPparameters_tw3(x[0:18])
    harpy.UpdateTables(1.0, 105.0)
    harpy.setNPparameters_SiversTMDPDF([x[18],0.0,0.0,0.0,0.0])
    print('np set =',["{:8.3f}".format(i) for i in x])        
            
    YY=DataProcessor.snowInterface_N2.ComputeXSec(setD2)
    ccD2,cc3=setD2.chi2(YY)    
    
    YY=DataProcessor.snowInterface_N2.ComputeXSec(setG2)
    ccG2,cc3=setG2.chi2(YY)    
    
    YY=DataProcessor.harpyInterface.ComputeXSec(setSivers,method="central")
    ccSivers,cc3=setSivers.chi2(YY)    
    
    chiTOTAL=(ccD2/setD2.numberOfPoints+ccG2/setG2.numberOfPoints+ccSivers/setSivers.numberOfPoints)*(
        setD2.numberOfPoints+setG2.numberOfPoints+setSivers.numberOfPoints)
    
    endT=time.time()
    print(':->',ccD2/setD2.numberOfPoints," ",ccG2/setG2.numberOfPoints," ",ccSivers/setSivers.numberOfPoints,"    time=",endT-startT)
    return chiTOTAL

#%%
from iminuit import Minuit

#---- PDFbias-like row (0.083931)
initialValues=(3.0,1.0, 
                0.05,0.,0.0,0.,
                -0.06,0.,0.0,0.,
                0.0,0.,0.0,0.,
                0.0,0.,0.0,0.,
                0.5)

initialErrors=(0.1,0.1, 
                0.1,0.1,0.1,0.1,
                0.1,0.1,0.1,0.1,
                0.1,0.1,0.1,0.1,
                0.1,0.1,0.1,0.1,
                0.1)
searchLimits=((1.,5.),(1.,10.),
              (-50.,50.), (-50.,50.), (-50.,50.),
              (-50.,50.), (-50.,50.), (-50.,50.),
              (-50.,50.), (-50.,50.), (-50.,50.),
              (-50.,50.), (-50.,50.), (-50.,50.),
              (-50.,50.), (-50.,50.), (-50.,50.), (-50.,50.),
              (0.001,5.))
              
# True= FIX
parametersToMinimize=(False, True,
                      False, False,False,False,
                      True, True,
                      False,False,
                      True,True,True,True,
                      True,False,False,True,
                      True)

#%%
rSet.SetReplica(0)

m = Minuit(chi2, initialValues)

m.errors=initialErrors
m.limits=searchLimits
m.fixed=parametersToMinimize
m.errordef=1

print(m.params)
#%%
#m.tol=0.0001*(setSIDIS.numberOfPoints+setDY.numberOfPoints)*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1
m.migrad()

print(m.params)

chi2(list(m.values))

DataProcessor.snowInterface_N2.PrintChi2Table(setD2,printDecomposedChi2=True)
DataProcessor.snowInterface_N2.PrintChi2Table(setG2,printDecomposedChi2=True)
DataProcessor.harpyInterface.PrintChi2Table(setSivers,method="central",printSysShift=False)

print([round(x,1 if x >100 else 4 if x>1 else 6) for x in list(m.values)])
