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
path_to_INI=DATAP_DIR+"FittingPrograms/Tw3_2026/INI/snowflake_forRep.ini"
harpy.initialize_snowflake(path_to_INI)

NP_par=numpy.zeros(24)+0.2
harpy.setNPparameters_tw3(NP_par)
harpy.UpdateTables(1.0, 105.0)

#%%
#######################################
#Initialize artemide
#######################################

import DataProcessor.ArtemideReplicaSet

path_to_constants=DATAP_DIR+"FittingPrograms/Tw3_2026/INI/TMD+tw3.atmde"


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
def cutFuncD2(p):
    if p["type"]=="D2":
        if p["<Q>"]<numpy.sqrt(2.):
            return False, p
    
    if "E143" in p["id"] and p["process"]==101:
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
        if(p["<pT>"]/p["<z>"]<0.0):
            delta=0.0
            #delta=p["<pT>"]/p["<z>"]/p["<Q>"]        
        else:
            delta=p["<pT>"]/p["<z>"]/p["<Q>"]        
        
    if "compass23" in p["id"]:
        deltaTEST=0.5
    if "compass08" in p["id"]:
        deltaTEST=0.47
    if "jlab" in p["id"]:
        deltaTEST=0.45
    if "JLab" in p["id"]:
        deltaTEST=0.45
    
    if delta<deltaTEST:
        pNew=copy.deepcopy(p)    
        pNew["process"]=pNew["weightProcess"]
        if p["type"]=="SIDIS":
            normX=DataProcessor.harpyInterface.ComputeXSec(pNew,method="central")        
        elif p["type"]=="DY":
            normX=DataProcessor.harpyInterface.ComputeXSec(pNew)        
        else:
            print("Are you crazy?")
        p["thFactor"]=p["thFactor"]/normX        #### this minus is because of star
    
    #### This is because star measures AN
    if "star" in p["id"]:
        p["thFactor"]=+p["thFactor"]
    
    
        
#    return delta<0.5 and p.qT_avarage<80
    return delta<deltaTEST and p["<Q>"]>1.41, p


#%%
### Loading the D2 data set
theDataD2=DataProcessor.DataMultiSet.DataMultiSet("D2set",loadThisDataD2([
    "E143_d2",#"E154_d2",
    "E155-1999_d2","E155_d2",
    "HallA-2016_d2","HERMES_d2","SANE_d2",       
    #"RSS-2006_d2","RSS-2008_d2",    
    "RQCD_d2_ud",
    #"RQCD_d2_singlet","RQCD_d2_pn",
    #"GHMP26_d2","QCDSF_d2"
    ]))

setD2=theDataD2.CutData(cutFuncD2) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setD2.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setD2.sets]), 'points.') 

#%%
### Loading the G2 data set
theDataG2=DataProcessor.DataMultiSet.DataMultiSet("G2set",loadThisDataG2([
    #"E142.n", 
    "E143.p", "E143.d",#"E143.n", 
    "E154.n",
    "E155-29.p","E155-32.p","E155-38.p",
    "E155-29.d","E155-32.d","E155-38.d",
    #"SMC.p",
    "HERMES",
    #"HallA-2004.n",
    "HallA-2016-4.He3","HallA-2016-5.He3"
    ]))

setG2=theDataG2.CutData(cutFunc) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setG2.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setG2.sets]), 'points.') 

#%%
### Loading the data set for Sivers
theDataS=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisDataSivers([
                    'compass08.sivers.pi+.dpt', 'compass08.sivers.pi-.dpt',
                    'compass08.sivers.k+.dpt', 'compass08.sivers.k-.dpt',
                    'compass08.sivers.pi+.dx', 'compass08.sivers.pi-.dx',
                    'compass08.sivers.k+.dx', 'compass08.sivers.k-.dx',
                    'compass08.sivers.pi+.dz', 'compass08.sivers.pi-.dz',
                    'compass08.sivers.k+.dz', 'compass08.sivers.k-.dz',
                    'compass16.sivers.h+.1<z<2.dpt','compass16.sivers.h-.1<z<2.dpt',
                    'compass16.sivers.h+.z>2.dpt' ,'compass16.sivers.h-.z>2.dpt',
                    'compass16.sivers.h+.1<z<2.dz','compass16.sivers.h-.1<z<2.dz',
                    'compass16.sivers.h+.z>2.dz' ,'compass16.sivers.h-.z>2.dz',
                    'compass16.sivers.h+.1<z<2.dx','compass16.sivers.h-.1<z<2.dx',
                    'compass16.sivers.h+.z>2.dx' ,'compass16.sivers.h-.z>2.dx',
                    'compass23.sivers.h+.dpt', 'compass23.sivers.h-.dpt',
                    'compass23.sivers.h+.dx', 'compass23.sivers.h-.dx',
                    'compass23.sivers.h+.dz', 'compass23.sivers.h-.dz',
                    'hermes.sivers.pi+.3d','hermes.sivers.pi-.3d',
                    'hermes.sivers.k+.3d','hermes.sivers.k-.3d',
                    'jlab.sivers.pi+','jlab.sivers.pi-','jlab.sivers.k+','jlab.sivers.k-'
                    ]))

theDataSdy=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataSivers([
                    'star26.sivers.W-.dy', 'star26.sivers.W+.dy',
                    'star23.sivers.Z'
                    ]))

setSivers=theDataS.CutData(cutFunc_TMD) 
setSiversDY=theDataSdy.CutData(cutFunc_TMD) 


print('Loaded (SIDIS)', setSivers.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setSivers.sets]), 'points.')
print('Loaded (DY)', setSiversDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setSiversDY.sets]), 'points.')
#print('Loaded SIDIS experiments are', [i.name for i in setSivers.sets])

#%%
### Loading the WGT data set
theDataW=DataProcessor.DataMultiSet.DataMultiSet("ALTset",loadThisDataWGT([
                      'hermes3D.ALT.pi+','hermes3D.ALT.pi-',
                      'hermes3D.ALT.k+','hermes3D.ALT.k-',
                      'compass16.ALT.h+.2<z.dpt','compass16.ALT.h-.2<z.dpt',
                      'compass16.ALT.h+.2<z.dz','compass16.ALT.h-.2<z.dz',
                      'compass16.ALT.h+.2<z.dx','compass16.ALT.h-.2<z.dx',
                      'JLab6.ALT.pi+','JLab6.ALT.pi-'
                      ]))

setALT=theDataW.CutData(cutFunc_TMD) 

print('Loaded ', setALT.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setALT.sets]), 'points.')

#%%
harpy.setNPparameters_tw3([5.7019, 1.10442, 5.7019,  -1.34357,
                           1.03315, 0.0, 7.16511, 0.0,
                           -0.351023, 3.50304, 
                           -0.412292, 0.0, -1.81263, 0.0,
                           -8.88206, -16.5434, 
                           -1.19804, 0.0, -6.93412, 0.0,
                           -5.04158, -2.37763, 
                           0.832746, 0.330435])

harpy.setNPparameters_tw3([6.7942, 1.4713, 6.9114, -2.254391, 1.7656, 3.1443, 8.7584, -0.178593, \
                           -0.148088, 6.1239, -1.280427, -6.120491, -6.367843, -3.974293, \
                               -7.151409, -15.917348, 1.5086, 17.6056, 2.8707, 5.8855, -0.605963, \
                                   -2.528064, 3.9286, 8.9297, 0.432021])
    
harpy.setNPparameters_tw3([    7.5694, 2.6313, 7.8063, -1.728782, 1.9034, 3.5572, 9.8381, -0.695355, \
-1.896509, 5.7794, -1.648215, -7.742637, -8.245201, -4.971573, \
-13.808959, -11.397206, 1.7975, 21.1902, 5.924, 5.8302, -1.555601, \
-5.609505, 4.5722, 3.485, 0.699004])
    
harpy.UpdateTables(1.0, 105.0)

harpy.setNPparameters_SiversTMDPDF([0.5,0.0])
harpy.setNPparameters_wgtTMDPDF([0.5,0.0])

#%%
harpy.setNPparameters_tw3([7.81257, 2.16977, 7.89609, -1.68363, 1.82945, 4.68576, 16.0388, \
3.26094, 20.9951, 13.9676, -1.73296, -3.26, -7.93252, -1.4382, \
15.4692, -4.32768, 1.73142, 4.06223, 1.10699, 4.85654, 9.35829, \
10.0198, 4.27852, 23.0836])


harpy.UpdateTables(1.0, 105.0)
#%%
harpy.setNPparameters_SiversTMDPDF([0.82,0.0])
harpy.setNPparameters_wgtTMDPDF([0.82,0.0])


#%%
DataProcessor.snowInterface_N2.PrintChi2Table(setD2,printDecomposedChi2=False)
DataProcessor.snowInterface_N2.PrintChi2Table(setG2,printDecomposedChi2=False)

DataProcessor.harpyInterface.PrintChi2Table(setSivers,method="central",printSysShift=False)
#DataProcessor.harpyInterface.PrintChi2Table(setSiversDY,printSysShift=False)
DataProcessor.harpyInterface.PrintChi2Table(setALT,method="central",printSysShift=False)

#%%
#######################################
# Minimisation
#######################################
import time

totN=setD2.numberOfPoints+setG2.numberOfPoints+setSivers.numberOfPoints+setSiversDY.numberOfPoints+setALT.numberOfPoints

def deformation(c):
    #v=0.3
    #return numpy.exp(v*(c**(1./v)-1))
    return c

def chi2(x):
    
    startT=time.time()
    harpy.setNPparameters_tw3(x[0:24])
    harpy.UpdateTables(1.0, 105.0)
    harpy.setNPparameters_SiversTMDPDF([x[24],0.0])
    harpy.setNPparameters_wgtTMDPDF([x[24],0.0])
    print('np set =['+", ".join(["{:8.3f}".format(i) for i in x])+"]")
            
    YY=DataProcessor.snowInterface_N2.ComputeXSec(setD2)
    ccD2,cc3=setD2.chi2(YY)    
    
    YY=DataProcessor.snowInterface_N2.ComputeXSec(setG2)
    ccG2,cc3=setG2.chi2(YY)    
    
    YY=DataProcessor.harpyInterface.ComputeXSec(setSivers,method="central")
    ccSivers,cc3=setSivers.chi2(YY)
    
    #YY=DataProcessor.harpyInterface.ComputeXSec(setSiversDY)
    #ccSiversDY,cc3=setSiversDY.chi2(YY)    
    ccSiversDY=0.
    
    YY=DataProcessor.harpyInterface.ComputeXSec(setALT,method="central")
    ccWGT,cc3=setALT.chi2(YY)
    
    chiTOTAL=(deformation(ccD2/setD2.numberOfPoints)+
              deformation(ccG2/setG2.numberOfPoints)+
              deformation((ccSivers+ccSiversDY)/(setSivers.numberOfPoints+0*setSiversDY.numberOfPoints))+                  
              deformation(ccWGT/setALT.numberOfPoints))*totN
   
    endT=time.time()
    print(':->    ',"{:.4f}".format(chiTOTAL/totN),"    = ("
          "{:.2f}".format(ccD2/setD2.numberOfPoints),
          " + ","{:.2f}".format(ccG2/setG2.numberOfPoints),
          " + ","{:.2f}".format(ccSivers/setSivers.numberOfPoints),
          " + ","{:.2f}".format(ccSiversDY/setSiversDY.numberOfPoints),
          " + ","{:.2f}".format(ccWGT/setALT.numberOfPoints),
          ")    time=",endT-startT)
    return chiTOTAL

#%%
from iminuit import Minuit

#---- PDFbias-like row (0.083931)
initialValues=(5., 1.1, 5.,  -1.3,
               1.0, 0.0, 7.1, 0.0,
               -0.4, 3.5, 
               -0.4, 0.0,-1.8,0.0,
               -8.8, -14., 
               -1.2, 0.0,-6.9,0.0,
               -5.0, -2.3, 
               0.8, 0.3 ,
               0.5)

initialErrors=(0.5, 0.1, 0.5, 0.5,
                1., 1., 1., 1., 
                1., 1., 
                1., 1., 1., 1., 
                1., 1.,  
                1., 1., 1., 1., 
                1., 1.,  
                1., 1.,
                0.1)
searchLimits=((1.,10.),(0.1,10.), (1.,10.) ,(-10.,0.95),
              (-50.,50.), (-50.,50.), (-50.,50.),(-50.,50.),
              (-50.,50.), (-50.,50.), 
              (-50.,50.), (-50.,50.), (-50.,50.),(-50.,50.),
              (-50.,50.), (-50.,50.), 
              (-50.,50.), (-50.,50.), (-50.,50.),(-50.,50.),
              (-50.,50.), (-50.,50.), 
              (-50.,50.), (-50.,50.),
              (0.1,3.))
              
# True= FIX
parametersToMinimize=(False, False,False,False,
                      False, False,False,False,
                      False, False,
                      False, False,False,False,
                      False, False,
                      False, False,False,False,
                      False, False,
                      False, False,
                      False)

#Default: None. If set to None, Minuit assumes the cost function is computed in double precision. 
#If the precision of the cost function is lower (because it computes in single precision, for example) 
#set this to some multiple of the smallest relative change of a parameter that still changes the function.
precisionToMinuit=0.0001

# The convergence is detected when edm < edm_max, where edm_max is calculated as
# Migrad: edm_max = 0.002 * tol * errordef (errordef=1. by default)
# Users can set tol (default: 0.1) to a different value to either speed up convergence
### AV: This gives chi2 good up to 4 digits, if multiply it by N chi2/N is good up to 4 digits
tolToMinuit=0.1*totN

#%%
rSet.SetReplica(0)

m = Minuit(chi2, initialValues)

m.errors=initialErrors
m.limits=searchLimits
m.fixed=parametersToMinimize
m.errordef=1

#Default: None. If set to None, Minuit assumes the cost function is computed in double precision. 
#If the precision of the cost function is lower (because it computes in single precision, for example) 
#set this to some multiple of the smallest relative change of a parameter that still changes the function.
m.precision=0.0001

# The convergence is detected when edm < edm_max, where edm_max is calculated as
# Migrad: edm_max = 0.002 * tol * errordef (errordef=1. by default)
# Users can set tol (default: 0.1) to a different value to either speed up convergence
### AV: This gives chi2 good up to 4 digits, if multiply it by N chi2/N is good up to 4 digits
m.tol=0.1*totN


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
DataProcessor.harpyInterface.PrintChi2Table(setSiversDY,printSysShift=False)
DataProcessor.harpyInterface.PrintChi2Table(setALT,method="central",printSysShift=False)

print([round(x,1 if x >100 else 4 if x>1 else 6) for x in list(m.values)])