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

logFile=DATAP_DIR+"FittingPrograms/Tw3_2026/LOGS/log_1.log"
repFile=DATAP_DIR+"FittingPrograms/Tw3_2026/REPS/replicas_1.dat"

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
        if(p["<pT>"]/p["<z>"]<0.7):
            delta=0.0
            #delta=p["<pT>"]/p["<z>"]/p["<Q>"]        
        else:
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


harpy.UpdateTables(1.0, 105.0)

harpy.setNPparameters_SiversTMDPDF([0.5,0.0])
harpy.setNPparameters_wgtTMDPDF([0.5,0.0])

#%%
harpy.setNPparameters_tw3([6.02889, 1.80943, 5.94243, -1.64152, 1.38313, 1.84317, 6.55815, \
0.587521, -0.414991, 3.87324, -0.647569, -1.04194, 3.3807, -1.50536, \
-8.50641, -13.1764, -0.74414, 3.31015, -6.80776, 0.715115, -4.10758, \
-2.68135, 0.484923, 0.988009])


harpy.UpdateTables(1.0, 105.0)

harpy.setNPparameters_SiversTMDPDF([0.81,0.0])
harpy.setNPparameters_wgtTMDPDF([0.81,0.0])

#%%
DataProcessor.snowInterface_N2.PrintChi2Table(setD2,printDecomposedChi2=False)
DataProcessor.snowInterface_N2.PrintChi2Table(setG2,printDecomposedChi2=False)

DataProcessor.harpyInterface.PrintChi2Table(setSivers,method="central",printSysShift=False)
DataProcessor.harpyInterface.PrintChi2Table(setSiversDY,printSysShift=False)
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
    
    YY=DataProcessor.harpyInterface.ComputeXSec(setSiversDY)
    ccSiversDY,cc3=setSiversDY.chi2(YY)    
    
    YY=DataProcessor.harpyInterface.ComputeXSec(setALT,method="central")
    ccWGT,cc3=setALT.chi2(YY)
    
    chiTOTAL=(deformation(ccD2/setD2.numberOfPoints)+
              deformation(ccG2/setG2.numberOfPoints)+
              deformation((ccSivers+ccSiversDY)/(setSivers.numberOfPoints+setSiversDY.numberOfPoints))+                  
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
#######################################
# Generate replica of data and compute chi2
#######################################
def MinForReplica():
    global setD20,setG20,setSivers0,setSiversDY0,setALT0,\
            initialValues,initialErrors,searchLimits,parametersToMinimize,ART25replica
        
    def repchi_2(x):     
        
        if(not(isinstance(numpy.sum(x),float))):
            print("HERE WE CATCH IT!")
            sys.exit()
        
        startT=time.time()
        harpy.setNPparameters_tw3(x[0:24])
        harpy.UpdateTables(1.0, 105.0)
        #harpy.UpdateTables(1.0, 10.5)
        harpy.setNPparameters_SiversTMDPDF([x[24],0.0])
        harpy.setNPparameters_wgtTMDPDF([x[24],0.0])
        print('np set =['+", ".join(["{:8.3f}".format(i) for i in x])+"]")
                
        YY=DataProcessor.snowInterface_N2.ComputeXSec(setD2rep)
        ccD2,cc3=setD2rep.chi2(YY)    
        #ccD2=0.
        
        YY=DataProcessor.snowInterface_N2.ComputeXSec(setG2rep)
        ccG2,cc3=setG2rep.chi2(YY)    
        #ccG2=0.
        
        YY=DataProcessor.harpyInterface.ComputeXSec(setSiversrep,method="central")
        ccSivers,cc3=setSiversrep.chi2(YY)    
        #ccSivers=0.
        
        YY=DataProcessor.harpyInterface.ComputeXSec(setSiversDYrep)
        ccSiversDY,cc3=setSiversDYrep.chi2(YY)    
        #ccSiversDY=0.
        
        YY=DataProcessor.harpyInterface.ComputeXSec(setALTrep,method="central")
        ccWGT,cc3=setALTrep.chi2(YY)
        #ccWGT=0.
        
        chiTOTAL=(deformation(ccD2/setD2.numberOfPoints)+
                  deformation(ccG2/setG2.numberOfPoints)+
                  deformation((ccSivers+ccSiversDY)/(setSivers.numberOfPoints+setSiversDY.numberOfPoints))+                  
                  deformation(ccWGT/setALT.numberOfPoints))*totN

        if chiTOTAL>1000000.:
            rSet.SetReplica(ART25replica) 
        
        
        endT=time.time()
        print(':->    ',"{:.4f}".format(chiTOTAL/totN),"    = ("
              "{:.2f}".format(ccD2/setD2.numberOfPoints),
              " + ","{:.2f}".format(ccG2/setG2.numberOfPoints),
              " + ","{:.2f}".format(ccSivers/setSivers.numberOfPoints),
              " + ","{:.2f}".format(ccSiversDY/setSiversDY.numberOfPoints),
              " + ","{:.2f}".format(ccWGT/setALT.numberOfPoints),
              ")    time=",endT-startT)
        return chiTOTAL

    
    setD2rep=setD20.GenerateReplica()
    setG2rep=setG20.GenerateReplica()
    setSiversrep=setSivers0.GenerateReplica()
    setSiversDYrep=setSiversDY0.GenerateReplica()
    setALTrep=setALT0.GenerateReplica()
    
    localM = Minuit(repchi_2, initialValues)
    
    localM.errors=initialErrors
    localM.limits=searchLimits
    localM.fixed=parametersToMinimize
    localM.errordef=1    
    ### tolerance is a bit larger than for the main fit
    localM.tol=precisionToMinuit
    localM.precision=precisionToMinuit
    localM.strategy=1

    localM.migrad()
    
    ### [chi^2, NP-parameters]
    return [localM.fval,list(localM.values)]

#%%
#######################################
# LOG save function
#######################################
savedTime=time.time()
def SaveToLog(text):
    global savedTime,logFile
    newTime=time.time()
    
    import socket
    PCname=socket.gethostname()
    
    passedTime=newTime-savedTime
    hours=int(passedTime/3600)
    minutes=int((passedTime-hours*3600)/60)
    seconds=int(passedTime-hours*3600-minutes*60)
    
    with open(logFile, 'a') as file:
        file.write(PCname+ ' : ' + time.ctime()+' :  [+'+str(hours)+':'+str(minutes)+':'+str(seconds)+' ]\n')
        file.write(' --> '+text+'\n')
        file.write('\n')
    savedTime=time.time()

#%%
rSet.SetReplica(0)
setD20=theDataD2.CutData(cutFunc) 
setG20=theDataG2.CutData(cutFunc) 
setSivers0=theDataS.CutData(cutFunc_TMD)
setSiversDY0=theDataSdy.CutData(cutFunc_TMD)
setALT0=theDataW.CutData(cutFunc_TMD) 

#######################################
# This is the main cicle. 
# It generates replica of data take random PDF and minimize it
# Save to log.
#######################################

NumberOfReplicas=350

ART25replica=0
rSet.SetReplica(0)

#SaveToLog(" =============================NEW START===================================================")

for i in range(NumberOfReplicas):
    print('---------------------------------------------------------------')
    print('------------REPLICA ',i,' from ',NumberOfReplicas,'------------------')
    print('---------------------------------------------------------------')
    savedTime=time.time()
    
    ## reset PDF        
    ART25replica=numpy.random.randint(rSet.numberOfReplicas)
    print("Start computation with ART25 replica "+str(ART25replica))
    #SaveToLog("Start computation with ART25 replica "+str(ART25replica))
    
    rSet.SetReplica(ART25replica) 
    print("Minimization started.")
    
    #### this is needed to recompute the weights
    setD20=theDataD2.CutData(cutFunc) 
    setG20=theDataG2.CutData(cutFunc) 
    setSivers0=theDataS.CutData(cutFunc_TMD)
    setSiversDY0=theDataSdy.CutData(cutFunc_TMD)
    setALT0=theDataW.CutData(cutFunc_TMD) 
    
    
    ## got to pseudo-data and minimization
    repRes=MinForReplica()
    print(repRes)
    print("Minimization finished.")    
    SaveToLog(" Minization with ART25 replica "+str(ART25replica)+" completed.")
    
    harpy.UpdateTables(1.0, 105.0)
    ## compute the chi2 for true data full
    mainD2, mainD2_2 =DataProcessor.snowInterface_N2.ComputeChi2(setD2)    
    mainG2, mainG2_2 =DataProcessor.snowInterface_N2.ComputeChi2(setG2)    
    mainSiv, mainSiv_2 =DataProcessor.harpyInterface.ComputeChi2(setSivers,method="central")    
    mainSivDY, mainSivDY_2 =DataProcessor.harpyInterface.ComputeChi2(setSiversDY)    
    mainALT, mainALT_2 =DataProcessor.harpyInterface.ComputeChi2(setALT,method="central") 
    print("Central chi^2  computed.)")
    
    ## save to file
    f=open(repFile,"a+")
    print('SAVING >>  ',f.name)
    ### [total chi^2(full data DY),total chi^2(full data SIDIS), 
    ### list of chi^2 for experiments( DY),list of chi^2 for experiments( SIDIS), 
    ### PDFreplica, FFpi-replica, FFk-replica,
    ### list of NP-parameters] 
    f.write(str([mainD2,mainG2,mainSiv,mainSivDY, mainALT,
                 mainD2_2,mainG2_2,mainSiv_2,mainSivDY_2,mainALT_2,ART25replica,repRes[1]])+"\n")
    f.close()  
