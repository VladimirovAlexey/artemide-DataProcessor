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
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"

ATMDE_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide-FLAVOR/"
#%%
import sys
sys.path.remove('/data/arTeMiDe_Repository/artemide/harpy')
sys.path.append(ATMDE_DIR+"harpy")

import numpy
sys.path.append(ROOT_DIR)
import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet
import DataProcessor.DataMultiSet

MAINPATH=ROOT_DIR 
#%%
#######################################
#Initialize artemide
#######################################
import harpy
path_to_constants=MAINPATH+"FittingPrograms/MVZ22/ConstantsFiles/SV19f_N4LL"

harpy.initialize(path_to_constants)

harpy.setNPparameters_TMDR([0.0,0.0434])
harpy.setNPparameters_uTMDPDF([0.253434, .4351, 
                               0.253434, .4351,
                               0.253434, .4351,
                               0.253434, .4351,
                               0.253434, .4351, 
                               100.,1.])
harpy.setNPparameters_uTMDFF([0.264,0.479,0.459,0.539]) 

#%%
### read the list of files and return the list of DataSets
def loadThisDataDY(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=ROOT_DIR+"DataLib/unpolDY/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

def loadThisDataSIDIS(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=ROOT_DIR+"DataLib/unpolSIDIS/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

#%%
##################Cut function
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
### Loading the DY data set
theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY([
                          'CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                          #'A7-00y10', 'A7-10y20','A7-20y24', 
                          'A8-00y04', 'A8-04y08', 'A8-08y12', 'A8-12y16', 'A8-16y20', 'A8-20y24', 
                          'A8-46Q66', 'A8-116Q150', 
                          'CMS7', 'CMS8', 
                          'CMS13-00y04','CMS13-04y08','CMS13-08y12','CMS13-12y16','CMS13-16y24',
                          'CMS13_dQ_106to170','CMS13_dQ_170to350','CMS13_dQ_350to1000',
                          'LHCb7', 'LHCb8', 'LHCb13_dy(2021)', 
                          'PHE200', 'STAR510', 
                          'E228-200', 'E228-300', 'E228-400', 
                          'E772','E605'
                          ]))

setDY=theData.CutData(cutFunc) 

### Loading the SIDIS data set
theData=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisDataSIDIS([
                      'hermes.p.vmsub.zxpt.pi+','hermes.p.vmsub.zxpt.pi-',
                      'hermes.d.vmsub.zxpt.pi+','hermes.d.vmsub.zxpt.pi-',
                      'hermes.p.vmsub.zxpt.k+','hermes.p.vmsub.zxpt.k-',
                      'hermes.d.vmsub.zxpt.k+','hermes.d.vmsub.zxpt.k-',
                      'compass.d.h+','compass.d.h-']))

setSIDIS=theData.CutData(cutFunc) 

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setSIDIS.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setSIDIS.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setSIDIS.sets])

#SaveToLog('Loaded '+ str(setDY.numberOfSets) + ' data sets with '+str(sum([i.numberOfPoints for i in setDY.sets])) + ' points. \n'
#+'Loaded experiments are '+str([i.name for i in setDY.sets]))

#%%
harpy.setNPparameters_TMDR([-0.112, 0.122])
harpy.setNPparameters_uTMDPDF([0.0875,0.155,0.004327,4.401,
                               0.000133,28.14,0.0,1.098,
                               0.0469,0.0,23.165,0.0])

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)
#DataProcessor.harpyInterface.PrintChi2Table(setSIDIS,printDecomposedChi2=True)


    
#%%
#######################################
# Minimisation
#######################################
import time
totalN=setDY.numberOfPoints+setSIDIS.numberOfPoints

def chi_2DY(x):
    startT=time.time()
    harpy.setNPparameters_TMDR([x[0],x[1]])
    harpy.setNPparameters_uTMDPDF(x[2:])
    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
    
    ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)
    #ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(setSIDIS)
    
    cc=ccDY2/setDY.numberOfPoints
    endT=time.time()
    print(':->',cc,'       t=',endT-startT)
    return ccDY2

def chi_2SIDIS(x):
    startT=time.time()
    harpy.setNPparameters_uTMDFF(x)
    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
    
    #ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)
    ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(setSIDIS)
    
    cc=ccSIDIS2/setSIDIS.numberOfPoints
    endT=time.time()
    print(':->',cc,'       t=',endT-startT)
    return ccSIDIS2

def chi_2(x):
    startT=time.time()
    harpy.setNPparameters(x)
    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
    
    ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)
    ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(setSIDIS)
    
    cc=(ccDY2+ccSIDIS2)/totalN
    endT=time.time()
    print(':->',cc,'       t=',endT-startT)
    return ccDY2+ccSIDIS2

#%%
#### Minimize DY
from iminuit import Minuit

#

#initialValues=(0.0426578,0.253434, 9.04351,0.253434, 9.04351,0.253434, 9.04351,0.253434, 9.04351,0.253434, 9.04351, 346.999)
initialValues=( -0.111,   0.123, 
               0.002, 0.132, 
               0.002, 2.998, 
               0.002, 25.502, 
               0.000, 1.264, 
               0.084, 0.000, 
               25.224, 0.000)

initialErrors=(0.01,0.01, 
               0.1,  1.0, 0.1, 1.0, 
               0.1,  1.0, 0.1,  1.0, 
               0.1,  1.0, 25.0, 25.0)
searchLimits=((-0.2,2),(0.,4.), 
              (0.,10.), (0.,100.),(0.,10.), (0.,100.),
              (0.,10.), (0.,100.),(0.,10.), (0.,100.),
              (0.,10.), (0.,100.),(0.,1000),(0.,1000.))

# True= FIX
parametersToMinimize=(False,False, False, False, False,False, False, False, False, False, False, False, False ,True)

m = Minuit(chi_2DY, initialValues)

m.errors=initialErrors
m.limits=searchLimits
m.fixed=parametersToMinimize
m.errordef=1

print(m.params)

m.tol=0.0001*setDY.numberOfPoints*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1
m.migrad()

print(m.params)

valsDY=list(m.values)

harpy.setNPparameters_TMDR([valsDY[0],valsDY[1]])
harpy.setNPparameters_uTMDPDF(valsDY[2:])

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)
DataProcessor.harpyInterface.PrintChi2Table(setSIDIS,printDecomposedChi2=True)
#%%

#### Minimize SIDIS
from iminuit import Minuit

initialValues=(0.44255, 0.415158, 0.457, 0.71595)

initialErrors=(0.1,   0.1,  0.1,    0.1)
searchLimits=((0.,5.),(0.,5.),(0.,5.),None)

# True= FIX
parametersToMinimize=(False, False, False, False)

m = Minuit(chi_2SIDIS, initialValues)

m.errors=initialErrors
m.limits=searchLimits
m.fixed=parametersToMinimize
m.errordef=1

print(m.params)

m.tol=0.0001*setSIDIS.numberOfPoints*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1
m.migrad()

valsSIDIS=list(m.values)
harpy.setNPparameters_uTMDFF(valsSIDIS)

print("DY>>>>>>",valsDY)
print("SIDIS>>>",valsSIDIS)

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)
DataProcessor.harpyInterface.PrintChi2Table(setSIDIS,printDecomposedChi2=True)

#%%
####################################
# Search for minimum
###################################

# m.migrad()

# print(m.params)

# SaveToLog("MINIMIZATION FINISHED",str(m.params))
# SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))

print([harpy.get_DNP(0.01*b+0.01,2.) for b in range(400)])