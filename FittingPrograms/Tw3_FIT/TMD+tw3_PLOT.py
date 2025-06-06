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

path_to_constants=ROOT_DIR+"FittingPrograms/Tw3_FIT/INI/TMD+tw3.atmde"


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


##################Cut function for WGT data
def cutFuncWGT(p):
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
### Loading the WGT data set
theData=DataProcessor.DataMultiSet.DataMultiSet("ALTset",loadThisDataWGT([
                      'hermes3D.ALT.pi+','hermes3D.ALT.pi-',
                      'hermes3D.ALT.k+','hermes3D.ALT.k-',
                      'compass16.ALT.h+.2<z.dpt','compass16.ALT.h-.2<z.dpt',
                      'compass16.ALT.h+.2<z.dz','compass16.ALT.h-.2<z.dz',
                      'compass16.ALT.h+.2<z.dx','compass16.ALT.h-.2<z.dx'#,
                      #'JLab6.ALT.pi+','JLab6.ALT.pi-'
                      ]))

setALT=theData.CutData(cutFuncWGT) 

print('Loaded ', setALT.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setALT.sets]), 'points.')

#%%
harpy.setNPparameters_tw3([3.112, -0.157, 
                            0.128, 1.57, 0.127, 0.794, 
                            -0.180, -1.594, -5.019, -5.872, 
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

harpy.UpdateTables(1.0, 105.0)

harpy.setNPparameters_SiversTMDPDF([0.5,0.0])
harpy.setNPparameters_wgtTMDPDF([0.5,0.0])

#%%
DataProcessor.snowInterface_N2.PrintChi2Table(setD2,printDecomposedChi2=False)
DataProcessor.snowInterface_N2.PrintChi2Table(setG2,printDecomposedChi2=False)

DataProcessor.harpyInterface.PrintChi2Table(setSivers,method="central",printSysShift=False)
DataProcessor.harpyInterface.PrintChi2Table(setALT,method="central",printSysShift=False)


#%%
# #### PLOT d2
# Qplot=0.1*numpy.array(range(40))+1
# dataP=harpy.D2List(Qplot, Qplot*0+100)
# dataN=harpy.D2List(Qplot, Qplot*0+101)
# dataD=harpy.D2List(Qplot, Qplot*0+102)

#%%
def X1X2(r,phi):
    if(0<=phi<1):
        return r*(1-phi),r*phi
    elif(1<=phi<2):
        return r*(1-phi),r
    elif(2<=phi<3):
        return -r,r*(3-phi)
    elif(3<=phi<4):
        return r*(phi-4),r*(3-phi)
    elif(4<=phi<5):
        return r*(phi-4),-r
    elif(5<=phi<6):
        return r,r*(phi-6)
    else:
        raise("ERRROR")

#%%
path_to_save="/data/WorkingFiles/TMD/Fit_Notes/Twist3_25/data/Tw3PDF/"
rList=[0.01*i+0.01 for i in range(9)]+[0.05*i+0.1 for i in range(19)]
phiList=[i/8 for i in range(48)]
toSave=[]

# u-quark

for r in rList:       
    for phi in phiList:
        x1,x2=X1X2(r,phi)
        if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
        if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
    
        toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, 2., 2)]]

with open(path_to_save+"T_u", 'w') as file:
    for s in toSave:
        file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')

# d-quark
        
toSave=[]

for r in rList:       
    for phi in phiList:
        x1,x2=X1X2(r,phi)
        if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
        if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
    
        toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, 2., 1)]]

with open(path_to_save+"T_d", 'w') as file:
    for s in toSave:
        file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')

# s-quark
toSave=[]

for r in rList:       
    for phi in phiList:
        x1,x2=X1X2(r,phi)
        if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
        if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
    
        toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, 2., 3)]]

with open(path_to_save+"T_s", 'w') as file:
    for s in toSave:
        file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')

# g-quark
toSave=[]

for r in rList:       
    for phi in phiList:
        x1,x2=X1X2(r,phi)
        if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
        if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
    
        toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, 2., 0)]]

with open(path_to_save+"T_g+", 'w') as file:
    for s in toSave:
        file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
        
toSave=[]

for r in rList:       
    for phi in phiList:
        x1,x2=X1X2(r,phi)
        if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
        if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
    
        toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, 2., -10)]]

with open(path_to_save+"T_g-", 'w') as file:
    for s in toSave:
        file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
