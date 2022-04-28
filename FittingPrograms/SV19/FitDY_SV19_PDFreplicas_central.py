#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: vla18041
"""

#######################################
# Paths
#######################################
PathToHarpy="/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/harpy"
PathToDataProcessor="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"
PathToDataLibrary=PathToDataProcessor+"DataLib/unpolDY/"
PathToLog=PathToDataProcessor+"FittingPrograms/SV19/LOGS/"
PathToConstantsFile=PathToDataProcessor+"FittingPrograms/SV19/Constants-files/DY_nnlo/const-NNPDF31_NNLO"

import sys
#sys.path.append(PathToHarpy)
sys.path.append(PathToDataProcessor)

## if this trigger is ON the LHC data will be fit only by shape
useNormalizedLHCdata=False

#### Starting and final replica (included)
StartReplica=1
FinalReplica=3

#%%
#######################################
# importing libraries
#######################################
import numpy
import time

import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet

#%%
#######################################
# Output paths
#######################################
import socket
PCname=socket.gethostname()

replicaFile=PathToLog+"SV19_NNPDF31_PDFreplica_noA7.txt"
logFile=PathToLog+PCname+".log"

#%%
#######################################
# LOG save function
#######################################
savedTime=time.time()
def SaveToLog(text):
    global savedTime,logFile
    newTime=time.time()
    
    passedTime=newTime-savedTime
    hours=int(passedTime/3600)
    minutes=int((passedTime-hours*3600)/60)
    seconds=int(passedTime-hours*3600-minutes*60)
    
    with open(logFile, 'a') as file:
        file.write(time.ctime()+' :  [+'+str(hours)+':'+str(minutes)+':'+str(seconds)+' ]\n')
        file.write(' --> '+text+'\n')
        file.write('\n')
    savedTime=time.time()

#%%
#######################################
#Initialize artemide
#######################################
import harpy

SaveToLog('Initialization with : \n'+PathToConstantsFile)

harpy.initialize(PathToConstantsFile)
initializationArray=[2.0340, 0.0299, 0.2512, 7.7572,334.6108, 2.4543,-4.8203, 0.1000,  0.0000]
harpy.setNPparameters(initializationArray)

#%%
#######################################
# read the list of files and return the list of DataSets
#######################################
def loadThisData(listOfNames):    
    import DataProcessor.DataSet
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(PathToDataLibrary+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

#%%
#######################################
# Data cut function
#######################################
def cutFunc(p):    
    par=1.0
    if p["type"]=="DY":
        if(p["xSec"]>0):
            err=numpy.sqrt(sum([i**2 for i in p["uncorrErr"]]))/p["xSec"]
        else:
            err=100.
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
    
    if p["type"]=="SIDIS":        
        if p["<z>"]>0.8:
            return False , p
        ## bins with low z drop
        if p["<z>"]<0.2:
            return False , p
        
        par=1.0
        if p["xSec"]<0.00000001:
            err=1
            delta=1
        else:
            ##############3 I MULTIPLY THE ERROR BY 100 (so it does not affect the cuts)
            err=10000#*numpy.sqrt(p.uncorrErrorsSquare)/p.xSec    
            gamma2=(2.0*p["M_target"]*p["<x>"]/p["<Q>"])**2
            rho2=(p["M_product"]/p["<z>"]/(p["<Q>"]))**2
            qT=p["<pT>"]/p["<z>"]*numpy.sqrt((1+gamma2)/(1-gamma2*rho2))
            delta=qT/(p["<Q>"])
            
            ### compute the largest possible qT (approximate)
            gamma2WORST=(2.0*p["M_target"]*p["x"][1]/p["<Q>"])**2
            # it is definitely not a TMD point
            if gamma2WORST*rho2>1:
                return False , p
            qTWORST=p["pT"][1]/p["z"][0]*numpy.sqrt((1+gamma2WORST)/(1-gamma2WORST*rho2))
    
            ## drop if qT>Q/2
            if qTWORST>p["<Q>"]/2:
                return False , p
    
        ### drop Q<2
        if p["<Q>"]<2 :
            return False , p
    
#    return delta<0.5 and p.qT_avarage<80
    return (delta<0.1 or (delta<0.25 and par/err*delta**2<1)) , p

#%%
#######################################
# Loading the data set
#######################################
theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisData(['CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                      'A7-00y10', 'A7-10y20','A7-20y24', 
                      'A8-00y04-norm' if useNormalizedLHCdata else 'A8-00y04',
                      'A8-04y08-norm' if useNormalizedLHCdata else 'A8-04y08',
                      'A8-08y12-norm' if useNormalizedLHCdata else 'A8-08y12',
                      'A8-12y16-norm' if useNormalizedLHCdata else 'A8-12y16',
                      'A8-16y20-norm' if useNormalizedLHCdata else 'A8-16y20',
                      'A8-20y24-norm' if useNormalizedLHCdata else 'A8-20y24',
                      'A8-16y20-norm' if useNormalizedLHCdata else 'A8-16y20',
                      'A8-116Q150-norm' if useNormalizedLHCdata else 'A8-116Q150',
                      'CMS7', 'CMS8', 
                      'LHCb7', 'LHCb8', 'LHCb13', 
                      'PHE200', 'E228-200', 'E228-300', 'E228-400', 
                      'E772',
                      'E605']))

setDY=theData.CutData(cutFunc) 

print('Loaded ', setDY.numberOfSets, 'data sets with ', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setDY.sets])

SaveToLog('Loaded '+ str(setDY.numberOfSets) + ' data sets with '+str(sum([i.numberOfPoints for i in setDY.sets])) + ' points. \n'
+'Loaded experiments are '+str([i.name for i in setDY.sets]))

#%%
if useNormalizedLHCdata:
    for s in setDY.sets:
        if s.name[0:4]=='LHCb':
            s.isNormalized=True
            
        if s.isNormalized:
            s.normalizationMethod='bestChi2'
#%%
harpy.setNPparameters([2.000, 0.030, 0.254, 8.550, 373.112,  2.529, -5.910, 0.000, 0.000])

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)
#%%
#######################################
# Main chi2 formula
#######################################
totalN=setDY.numberOfPoints

def chi_2(x):
    startT=time.time()
    harpy.setNPparameters(x)
    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
    
    ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)
    
    cc=ccDY2/totalN
    endT=time.time()
    print(':->',cc,'       t=',endT-startT)
    return ccDY2

#%%
#######################################
# Setting-up Minuit
#######################################
from iminuit import Minuit

initialValues=(2.00,      0.039,      0.27,     8.9,      350.,     2.5,    -5.92, 0.,  0.)
initialErrors=(0.1,       0.01,       0.05,     0.2,      50.,      0.05,   0.05,  0.1, 0.1)
searchLimits=((1.4,4.5), (0.0001,5.0),(0.0,2.0),(0.,16.0),(0.,1000),(0.,5),  (-15,5),   None, None)

parametersToMinimize=(True,     False,    False,    False,    False,     False,  False, True, True)


#%%

m = Minuit.from_array_func(chi_2, initialValues,
      error=initialErrors, limit=searchLimits, fix=parametersToMinimize, errordef=1)

#m.get_param_states()

m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1
m.migrad()

print(m.params)

sys.exit()
#%%
#######################################
# Generate replica of data and compute chi2
#######################################
def MinForReplica():
    global totalN,setDY,initialValues,initialErrors,searchLimits,parametersToMinimize
        
    def repchi_2(x):        
        global totalN
        startT=time.time()
        harpy.setNPparameters(x)
        print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
        
        ccDY2,cc2=DataProcessor.harpyInterface.ComputeChi2(repDataDY)
        
        cc=(ccDY2)/totalN
        endT=time.time()
        print(':->',cc,'       t=',endT-startT)
        return ccDY2
    
    repDataDY=setDY
    
    localM = Minuit.from_array_func(repchi_2, initialValues,
      error=initialErrors, limit=searchLimits, fix=parametersToMinimize, errordef=1)
    
    localM.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
    localM.strategy=1

    localM.migrad()
    
    ### [chi^2, NP-parameters]
    return [localM.fval,localM.values.values()]

#%%
#######################################
# This is the main cicle. 
# It generates replica of data take random PDF and minimize it
# Save to log.
#######################################
for i in range(StartReplica,FinalReplica+1):
    print('---------------------------------------------------------------')
    print('------------REPLICA ',i,' from [',StartReplica,' , ',FinalReplica,']------------------')
    print('---------------------------------------------------------------')
    
    ## reset PDF
    harpy.setNPparameters(initializationArray)    
    harpy.setPDFreplica(i)
    SaveToLog("Start computation of replica "+str(i) +"in ["+ str(StartReplica)+','+str(FinalReplica)+"]")
    
    ## got to pseudo-data and minimization
    repRes=MinForReplica()
    print(repRes)
    SaveToLog("Minimization for replica "+str(i) +"in ["+ str(StartReplica)+','+str(FinalReplica)+"] finished.")    
    
    ## compute the chi2 for true data
    mainDY, mainDY2 =DataProcessor.harpyInterface.ComputeChi2(setDY)    
    SaveToLog("Central chi^2 for "+str(i) +"in ["+ str(StartReplica)+','+str(FinalReplica)+" computed. \n Saving to log >> "+replicaFile)
    
    ## save to file
    f=open(replicaFile,"a+")
    print('SAVING >>  ',f.name)
    ### [total chi^2(cenral), total chi^2 (pseudo data), list of chi^2 for experiments(central), number of PDF, list of NP-parameters]
    f.write(str([mainDY,repRes[0],mainDY2,i,repRes[1]])+"\n")
    f.close()
    
#%%
#######################################
# Finalizing log
#######################################
print('Computation finished')
SaveToLog("Computation finished correctly (["+ str(StartReplica)+','+str(FinalReplica)+"] range of replicas computed) +\n "
          +"----------------------------------------------------------------------------------")