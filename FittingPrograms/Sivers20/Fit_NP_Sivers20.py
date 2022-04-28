#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 16:40:59 2019

This code fit Sivers central-data for all replicas of SV19-input

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

useOrder="nnlo"
#useOrder="n3lo"

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
    
    unSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(
    "/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/Replicas/DY+SIDIS/SV19_nnlo_all=0.rep")
    
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
    
    unSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(
    "/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/Replicas/DY+SIDIS/SV19_n3lo_all=0.rep")
    

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
        p["thFactor"]=1./normX
    
    #### This is because star measures AN
    if p["id"][0:4]=="star":
        p["thFactor"]=-p["thFactor"]
        
#    return delta<0.5 and p.qT_avarage<80
    return delta<deltaTEST, p

#%%
### Loading the data set
theData0=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisData([
                    'compass08.sivers.pi+.dpt', 'compass08.sivers.pi-.dpt',
                    'compass08.sivers.k+.dpt', 'compass08.sivers.k-.dpt',
                    'compass16.sivers.h+.1<z<2.dpt','compass16.sivers.h-.1<z<2.dpt',
                    'compass16.sivers.h+.z>2.dpt' ,'compass16.sivers.h-.z>2.dpt',
                    'hermes.sivers.pi+.3d','hermes.sivers.pi-.3d',
                    'hermes.sivers.k+.3d','hermes.sivers.k-.3d',
                    'jlab.sivers.pi+','jlab.sivers.pi-','jlab.sivers.k+'
                    ]))

setSIDIS=theData0.CutData(cutFunc) 

theData1=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisData([
                    'star.sivers.W+.dqT','star.sivers.W-.dqT',
                    'star.sivers.Z',
                    'compass.sivers.piDY.dqT'
                    ]))

setDY=theData1.CutData(cutFunc) 

print('Loaded (SIDIS)', setSIDIS.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setSIDIS.sets]), 'points.')
print('Loaded SIDIS experiments are', [i.name for i in setSIDIS.sets])

print('Loaded (DY)', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded DY experiments are', [i.name for i in setDY.sets])

print('Total number of points:',setSIDIS.numberOfPoints+setDY.numberOfPoints)

#%%
### Set unpolarized TMD according to TMD set (adds constant pion row)
def SetUnTMD(n):
    unSet.SetReplica(n,part="TMDR")    
    unSet.SetReplica(n,part="uTMDFF")
    rr=unSet.GetReplica(n,part="uTMDPDF")
    harpy.setNPparameters_uTMDPDF(rr[:-1]+[0.0014, 0.442, 4.14])

#%%
rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/"+
                                                  "Sivers20_BPV20(n3lo).rep")
rSet.SetReplica()

#harpy.setNPparameters_SiversTMDPDF([0.199479, 0.263377, 59.3553, 0., 0., -0.0263578, -0.266869, -2.83529, 0.0955288, -0.595511, 2.90686, 0.315239, 2.89282, -0.136414])

DataProcessor.harpyInterface.PrintChi2Table(setSIDIS,method="central",printSysShift=False)

DataProcessor.harpyInterface.PrintChi2Table(setDY)

#%%
#######################################
# Minimisation
#######################################
includeDY=True
if includeDY:
    totalN=setSIDIS.numberOfPoints+setDY.numberOfPoints
else :
    totalN=setSIDIS.numberOfPoints

def chi_2(x):
    startT=time.time()
    harpy.setNPparameters_SiversTMDPDF(x)
    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
    
    ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(setSIDIS,method="central")
    
    if includeDY:
        ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)
    else:
        ccDY2,cc3=0,0
    
    cc=(ccSIDIS2+ccDY2)/totalN
    endT=time.time()
    print(':->',cc,'       t=',endT-startT)
    return ccSIDIS2+ccDY2
#%%
from iminuit import Minuit

if(useOrder=="nnlo"):
    initialValues=(0.526542, 4.29217, 172.381, 0.0, 0.0, -0.0178186, -0.353326, -3.88268, 0.38586, -0.372215, 11.844, 0.996171, 2.51856, -0.539914)
elif(useOrder=="n3lo"):
    initialValues=(0.526542, 4.29217, 172.381, 0.0, 0.0, -0.0178186, -0.353326, -3.88268, 0.38586, -0.372215, 11.844, 0.996171, 2.51856, -0.539914)

initialErrors=(0.1, 10., 10. , 0.1,0.1,
               0.3, 0.5, 1., 
               1., 1., 1., 
               1., 1., 0.1)
searchLimits=((0,None),(0,None), (0,None), None, None,
              (-5.,5.),(-0.99,30), (-30.,30.),
              (-5.,5.),(-0.99,60.),(-30.,30.),
              (-5.,15.),(-0.99,30.), (-5,5))

parametersToMinimize=(False,False,False,True,True, 
                      False, False, False, 
                      False, False, False,
                      False, False, False)


m = Minuit.from_array_func(chi_2, initialValues,
      error=initialErrors, limit=searchLimits, fix=parametersToMinimize, errordef=1)

#m.get_param_states()

m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1
#%%

def MinForReplica(n):
    
    
    def repchi_2(x):        
        startT=time.time()
        harpy.setNPparameters_SiversTMDPDF(x)
        print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
        
        
        if(includeDY):
            ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(repDataDY)
        else:
            ccDY2,cc3=0,0
            
        ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(repDataSIDIS,method="central")
        cc=(ccDY2+ccSIDIS2)/totalNnew
        
        endT=time.time()
        print(':->',cc,'       t=',endT-startT)
        return ccSIDIS2+ccDY2
    
    SetUnTMD(n)
    
    repDataSIDIS=theData0.CutData(cutFunc) 
    
    if(includeDY):
        repDataDY=theData1.CutData(cutFunc) 
        totalNnew=repDataSIDIS.numberOfPoints+repDataDY.numberOfPoints
    else:
        totalNnew=repDataSIDIS.numberOfPoints    
    
    localM = Minuit.from_array_func(repchi_2, initialValues,
      error=initialErrors, limit=searchLimits, fix=parametersToMinimize, errordef=1)
    
    localM.tol=0.0001*totalNnew*10000 ### the last 0.0001 is to compensate MINUIT def
    localM.strategy=1

    localM.migrad(ncall=250,resume=False)    
    
    SetUnTMD(0)
    chi2Central=chi_2(localM.values.values())
    
    return [localM.fval,chi2Central,localM.values.values()]

#%%
#
# Generate pseudo data and minimise   100 times
#
rStart=150
rFinish=301
REPPATH=MAINPATH+"FittingPrograms/Sivers20/LOGS/"+"final(nnlo)-NPreplicas.txt"
for i in range(rStart,rFinish):
    print('---------------------------------------------------------------')
    print('------------REPLICA: ',rStart," / ",i,' / ',rFinish,'--------------------')
    print('---------------------------------------------------------------')
    repRes=MinForReplica(i)
    print(repRes)
    f=open(REPPATH,"a+")
    print('SAVING >>  ',f.name)
    f.write(str(repRes+[i])+"\n")
    f.close()