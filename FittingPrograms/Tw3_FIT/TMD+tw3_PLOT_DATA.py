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
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/"
DATAP_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"

SNOWFLAKE_DIR = ROOT_DIR+"artemide/harpy/"
MODEL_DIR = ROOT_DIR+"artemide/Models/ART25/Replica-files/"

logFile=DATAP_DIR+"FittingPrograms/Tw3_FIT/log5.log"
repFile=DATAP_DIR+"FittingPrograms/Tw3_FIT/replicas_5.dat"

import sys
import numpy
if('/data/arTeMiDe_Repository/artemide/harpy' in sys.path):
    sys.path.remove('/data/arTeMiDe_Repository/artemide/harpy')
sys.path.append(DATAP_DIR)
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
path_to_INI=DATAP_DIR+"FittingPrograms/Tw3_FIT/INI/snowflake_forRep.ini"
#path_to_INI=DATAP_DIR+"FittingPrograms/Tw3_FIT/INI/TEST.ini"
#path_to_INI=DATAP_DIR+"FittingPrograms/Tw3_FIT/INI/TEST_16x8.ini"
harpy.initialize_snowflake(path_to_INI)

NP_par=numpy.zeros(18)+0.2
harpy.setNPparameters_tw3(NP_par)
harpy.UpdateTables(1.0, 105.0)

#%%
#######################################
#Initialize artemide
#######################################

import DataProcessor.ArtemideReplicaSet

path_to_constants=DATAP_DIR+"FittingPrograms/Tw3_FIT/INI/TMD+tw3.atmde"


harpy.initialize(path_to_constants)

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(MODEL_DIR+"ART25_main.rep")
    
rSet.SetReplica(0)

#%%
### read the list of files and return the list of DataSets
def loadThisDataD2(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=DATAP_DIR+"DataLib/D2_moment/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        
    return dataCollection


### read the list of files and return the list of DataSets
def loadThisDataG2(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=DATAP_DIR+"DataLib/G2/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        
    return dataCollection


### read the list of files and return the list of DataSets
def loadThisDataSivers(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=DATAP_DIR+"DataLib/Sivers/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection


def loadThisDataWGT(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=DATAP_DIR+"DataLib/wgt/"
    
    
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
harpy.setNPparameters_tw3([5.96232, 1.07002, -1.47505, 1.12059, 8.70425, 0.403322, 3.4056, \
-0.539385, -2.18672, -10.4606, -22.2453, -1.2878, -8.20021, 4.78427, \
1.2689, -3.53346, 1.70584, 0.])

harpy.UpdateTables(1.0, 105.0)

harpy.setNPparameters_SiversTMDPDF([0.5,0.0])
harpy.setNPparameters_wgtTMDPDF([0.5,0.0])

#%%
DataProcessor.snowInterface_N2.PrintChi2Table(setD2,printDecomposedChi2=False)
DataProcessor.snowInterface_N2.PrintChi2Table(setG2,printDecomposedChi2=False)

DataProcessor.harpyInterface.PrintChi2Table(setSivers,method="central",printSysShift=False)
DataProcessor.harpyInterface.PrintChi2Table(setALT,method="central",printSysShift=False)

#%%
repLIST=[]
#with open("/data/WorkingFiles/TMD/Fit_Notes/Twist3_25/rreps4.csv","r") as file:
with open("/data/WorkingFiles/TMD/Fit_Notes/Twist3_25/rreps5.csv","r") as file:
    lines=file.readlines()
    for line in lines:
        repLIST.append([float(p) for p in line.split(",")])       

#%%
####### Saving G2
#path_to_save="/data/WorkingFiles/TMD/Fit_Notes/Twist3_25/data/Tw3PDF_4GeV_n4/"
path_to_save="/data/WorkingFiles/TMD/Fit_Notes/Twist3_25/data_for_plots/"


for n in range(len(repLIST)):
    print("---->",n)
    harpy.setNPparameters_tw3(repLIST[n])
    harpy.UpdateTables(1.0, 8.0)
    
    # X=0.5 vs Q
    X=0.5
    f=100
    Qlist=numpy.arange(1.2,8.2,0.2)
    Xlist=[X for i in range(len(Qlist))]
    flist=[f for i in range(len(Qlist))]
        
    G2list=harpy.G2List(Xlist,Qlist,flist)

    
    with open(path_to_save+"G2_vsQ/X05."+str(n).zfill(4), 'w') as file:
        file.write(', '.join("{:.10f}".format(item) for item in G2list)+'\n')

    # X=0.3 vs Q
    X=0.3
    f=100
    Qlist=numpy.arange(1.2,8.2,0.2)
    Xlist=[X for i in range(len(Qlist))]
    flist=[f for i in range(len(Qlist))]
        
    G2list=harpy.G2List(Xlist,Qlist,flist)

    
    with open(path_to_save+"G2_vsQ/X03."+str(n).zfill(4), 'w') as file:
        file.write(', '.join("{:.10f}".format(item) for item in G2list)+'\n')
        
    # X=0.05 vs Q
    X=0.05
    f=100
    Qlist=numpy.arange(1.2,8.2,0.2)
    Xlist=[X for i in range(len(Qlist))]
    flist=[f for i in range(len(Qlist))]
        
    G2list=harpy.G2List(Xlist,Qlist,flist)

    
    with open(path_to_save+"G2_vsQ/X005."+str(n).zfill(4), 'w') as file:
        file.write(', '.join("{:.10f}".format(item) for item in G2list)+'\n')
    
    # Q=2.5 vs x
    Q=2.5
    f=100
    Xlist=numpy.arange(0.01,1,0.01 )
    Qlist=[Q for i in range(len(Xlist))]
    flist=[f for i in range(len(Qlist))]
        
    G2list=harpy.G2List(Xlist,Qlist,flist)

    
    with open(path_to_save+"G2_vsX/Q25."+str(n).zfill(4), 'w') as file:
        file.write(', '.join("{:.10f}".format(item) for item in G2list)+'\n')
        
    # Q=8. vs x
    Q=8.
    f=100
    Xlist=numpy.arange(0.01,1,0.01 )
    Qlist=[Q for i in range(len(Xlist))]
    flist=[f for i in range(len(Qlist))]
        
    G2list=harpy.G2List(Xlist,Qlist,flist)

    
    with open(path_to_save+"G2_vsX/Q80."+str(n).zfill(4), 'w') as file:
        file.write(', '.join("{:.10f}".format(item) for item in G2list)+'\n')
#%%
#### Plot G2 vs. Q
X=0.3
f=100
Qlist=numpy.arange(1.4,8.,0.2)
Xlist=[X for i in range(len(Qlist))]
flist=[f for i in range(len(Qlist))]

G2list=harpy.G2List(Xlist,Qlist,flist)

#%%
#### Plot G2 vs. Q
Q=2.5
f=100
Xlist=numpy.arange(0.01,1,0.01 )
Qlist=[Q for i in range(len(Xlist))]
flist=[f for i in range(len(Qlist))]

G2list=harpy.G2List(Xlist,Qlist,flist)
