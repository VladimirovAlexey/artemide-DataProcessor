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

SNOWFLAKE_DIR = "/data/arTeMiDe_Repository/artemide/harpy/"
MODEL_DIR = "/data/arTeMiDe_Repository/artemide/Models/ART25/Replica-files/"

logFile=ROOT_DIR+"FittingPrograms/Tw3_FIT/log.log"
repFile=ROOT_DIR+"FittingPrograms/Tw3_FIT/replicas_1.dat"

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
path_to_INI=ROOT_DIR+"FittingPrograms/Tw3_FIT/INI/snowflake_forRep.ini"
#path_to_INI=ROOT_DIR+"FittingPrograms/Tw3_FIT/INI/TEST.ini"
#path_to_INI=ROOT_DIR+"FittingPrograms/Tw3_FIT/INI/TEST_16x8.ini"
harpy.initialize_snowflake(path_to_INI)

NP_par=numpy.zeros(18)+0.2
harpy.setNPparameters_tw3(NP_par)
harpy.UpdateTables(1.0, 105.0)

#%%
#######################################
#Initialize artemide
#######################################

import DataProcessor.ArtemideReplicaSet

path_to_constants=ROOT_DIR+"FittingPrograms/Tw3_FIT/INI/TMD+tw3.atmde"


harpy.initialize(path_to_constants)

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(MODEL_DIR+"ART25_main.rep")
    
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


### read the list of files and return the list of DataSets
def loadThisDataG2(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=ROOT_DIR+"DataLib/G2/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        
    return dataCollection


### read the list of files and return the list of DataSets
def loadThisDataSivers(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data="/data/arTeMiDe_Repository/DataProcessor/DataLib/Sivers/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection


def loadThisDataWGT(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data="/data/arTeMiDe_Repository/DataProcessor/DataLib/wgt/"
    
    
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


##################Cut function
def cutFunc_TMD(p):
    import copy
    
    if p["type"]=="DY":
        deltaTEST=0.3
        delta=p["<qT>"]/p["<Q>"]        
        
        
        if(9<p["<Q>"]<11):#UPSILON resonance-bin
            return False , p
    
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
    
    #### This is because star measures AN
    #if p["id"][0:4]=="star":
    #    p["thFactor"]=-p["thFactor"]        
    #if p["type"]=="DY":
    #    p["thFactor"]=-p["thFactor"]        
        
#    return delta<0.5 and p.qT_avarage<80
    return delta<deltaTEST and p["<Q>"]>1.41, p


#%%
### Loading the D2 data set
theDataD2=DataProcessor.DataMultiSet.DataMultiSet("D2set",loadThisDataD2([
    #"E143_d2","E154_d2","E155-1999_d2","E155_d2",
    "HallA-2016_d2","HERMES_d2","SANE_d2",
    "RQCD_d2_ud"
    ]))

setD2=theDataD2.CutData(cutFunc) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setD2.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setD2.sets]), 'points.') 

#%%
### Loading the G2 data set
theDataG2=DataProcessor.DataMultiSet.DataMultiSet("G2set",loadThisDataG2([
    #"E142.n", "E143.p", "E143.d","E143.n", 
    "E154.n",
    "E155-29.p","E155-32.p","E155-38.p",
    #"E155-29.d","E155-32.d","E155-38.d",
    #"SMC.p",
    "HERMES",
    "HallA-2004.n","HallA-2016-4.He3","HallA-2016-5.He3",
    ]))

setG2=theDataG2.CutData(cutFunc) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setG2.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setG2.sets]), 'points.') 

#%%
### Loading the data set for Sivers
theDataS=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisDataSivers([
                    'compass08.sivers.pi+.dpt', 'compass08.sivers.pi-.dpt',
                    'compass08.sivers.k+.dpt', 'compass08.sivers.k-.dpt',
                    'compass16.sivers.h+.1<z<2.dpt','compass16.sivers.h-.1<z<2.dpt',
                    'compass16.sivers.h+.z>2.dpt' ,'compass16.sivers.h-.z>2.dpt',
                    'hermes.sivers.pi+.3d','hermes.sivers.pi-.3d',
                    'hermes.sivers.k+.3d','hermes.sivers.k-.3d',
                    'jlab.sivers.pi+','jlab.sivers.pi-','jlab.sivers.k+'
                    ]))

setSivers=theDataS.CutData(cutFunc_TMD) 


print('Loaded (SIDIS)', setSivers.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setSivers.sets]), 'points.')
#print('Loaded SIDIS experiments are', [i.name for i in setSivers.sets])

#%%
### Loading the WGT data set
theDataW=DataProcessor.DataMultiSet.DataMultiSet("ALTset",loadThisDataWGT([
                      'hermes3D.ALT.pi+','hermes3D.ALT.pi-',
                      'hermes3D.ALT.k+','hermes3D.ALT.k-',
                      'compass16.ALT.h+.2<z.dpt','compass16.ALT.h-.2<z.dpt',
                      'compass16.ALT.h+.2<z.dz','compass16.ALT.h-.2<z.dz',
                      'compass16.ALT.h+.2<z.dx','compass16.ALT.h-.2<z.dx'#,
                      #'JLab6.ALT.pi+','JLab6.ALT.pi-'
                      ]))

setALT=theDataW.CutData(cutFunc_TMD) 

print('Loaded ', setALT.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setALT.sets]), 'points.')

#%%
harpy.setNPparameters_tw3([3.112,   -0.157,   
                           0.128,    1.570,    0.127,    0.794,  
                           -0.180,   -1.594,   -5.019,   -5.872,   
                           0.000,    0.000,    0.000,    0.000,    
                           0.000,    0.000,    0.000,    0.000])
harpy.UpdateTables(1.0, 105.0)

#%%
rSet.SetReplica(330)
harpy.setNPparameters_SiversTMDPDF([0.5,0.0])
harpy.setNPparameters_wgtTMDPDF([0.5,0.0])
#%%
DataProcessor.snowInterface_N2.PrintChi2Table(setD2,printDecomposedChi2=False)
DataProcessor.snowInterface_N2.PrintChi2Table(setG2,printDecomposedChi2=False)

DataProcessor.harpyInterface.PrintChi2Table(setSivers,method="central",printSysShift=False)
DataProcessor.harpyInterface.PrintChi2Table(setALT,method="central",printSysShift=False)

#%%
sys.exit()
