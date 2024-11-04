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

ATMDE_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide/"
HARPY_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide/harpy/"


import sys
import numpy
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)


#%%
import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet

MAINPATH=ROOT_DIR 
#%%
#######################################
#Initialize artemide
#######################################
import harpy

#path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/ART25_N4LL.atmde"
path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/ART25_N4LL_DYresum.atmde"


harpy.initialize(path_to_constants)

inARRAY_TMDR=[1.5004, 0.073018, 0.038048, 0.0]
inARRAY_PDF=[0.521462, 0.000206, 0.402948, 7.0219, 1.0, 20.4051, 1.0, 0.000123, 1.1037, 0.660734, 0.0, 0.04]
inARRAY_FF=[0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0]


harpy.setNPparameters_TMDR(inARRAY_TMDR)
harpy.setNPparameters_uTMDPDF(inARRAY_PDF)
harpy.setNPparameters_uTMDFF(inARRAY_FF)

#%%
### read the list of files and return the list of DataSets
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
### Loading the SIDIS data set
theData=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisDataSIDIS([
                      'hermes.p.vmsub.zxpt.pi+','hermes.p.vmsub.zxpt.pi-',
                      'hermes.d.vmsub.zxpt.pi+','hermes.d.vmsub.zxpt.pi-',
                      'hermes.p.vmsub.zxpt.k+','hermes.p.vmsub.zxpt.k-',
                      'hermes.d.vmsub.zxpt.k+','hermes.d.vmsub.zxpt.k-',
                      'compass.d.h+','compass.d.h-']))

setSIDIS=theData.CutData(cutFunc) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setSIDIS.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setSIDIS.sets]), 'points.')

#%%
### ART23*+resum
harpy.setNPparameters_TMDR([1.5004, 0.073018, 0.038048, 0.0])
harpy.setNPparameters_uTMDPDF([0.521462, 0.000206, 0.402948, 7.0219, 1.0, 20.4051, 1.0, 0.000123, 1.1037, 0.660734, 0.0, 0.04])

### ART23*
#harpy.setNPparameters_TMDR([1.5004, 0.060236, 0.037013, 0.0])
#harpy.setNPparameters_uTMDPDF([0.602249, 0.000103, 0.445843, 7.2217, 1.0, 18.6931, 1.0, 0.000101, 1.1114, 1.3137, 0.0, 0.04])

harpy.setNPparameters_uTMDFF([0.704027, 0.780902, 0.711725, 0.372527,0.704027, 0.780902, 0.711725, 0.372527])
DataProcessor.harpyInterface.PrintChi2Table(setSIDIS,printDecomposedChi2=True)

#%%
#harpy.setNPparameters_TMDR([1.5, 0.0726487531256797, 0.03093216422113098, 1.0])
#harpy.setNPparameters_uTMDFF([0.5665109139577634, 0.3023672745403205, 0.3023672745403205, 0.3023672745403205, 0.5665109139577634, 0.2769242313483304, 0.2769242313483304, 0.2769242313483304])
#harpy.setNPparameters_uTMDPDF([0.5108596058488025, 0.3322138154268387, 0.41144675446362294, 1.5329290228129735, 0.5108596058488025, 24.772609820860936, 0.41144675446362294, 0.003779315753461415, 0.6103667356471426, 0.00789980377668142, 0.0, 0.0])

#%%
#harpy.setNPparameters_TMDR([1.5004, 0.05614, 0.03862, 1.0])
#harpy.setNPparameters_uTMDPDF([0.565, 0.0539, 0.5697, 6.64, 0.565, 20.07, 0.5697, 0.537, 1.07, 2.39, 0.0, 0.0])

#%%
#######################################
# Minimisation
#######################################
import time
totalN=setSIDIS.numberOfPoints

def chi_2SIDIS(x):
    startT=time.time()
    harpy.setNPparameters_uTMDFF([x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7]])
    print('np set =',["{:8.3f}".format(i) for i in x])        
    
    YY=DataProcessor.harpyInterface.ComputeXSec(setSIDIS)
    ccSIDIS2,cc3=setSIDIS.chi2(YY)
    
    cc=ccSIDIS2/setSIDIS.numberOfPoints    
    
    endT=time.time()
    print(':->',cc,"    time=",endT-startT)
    return ccSIDIS2

#%%
#### Minimize SIDIS
from iminuit import Minuit

#---- PDFbias-like row
initialValues=([0.704027, 0.780902, 0.711725, 0.372527,0.704027, 0.780902, 0.711725, 0.372527])

initialErrors=(0.1,0.1,0.1,0.1,
               0.1,0.1,0.1,0.1)
searchLimits=((0.0001,100.), (-100.,100.),(-100.,100.),(-100.,100.),
              (0.0001,100.), (-100.,100.),(-100.,100.),(-100.,100.))
              
# True= FIX
parametersToMinimize=(False, False, False,False,
                      False, False, False,False)

#%%

m = Minuit(chi_2SIDIS, initialValues)

m.errors=initialErrors
m.limits=searchLimits
m.fixed=parametersToMinimize
m.errordef=1

print(m.params)

m.tol=0.0001*setSIDIS.numberOfPoints*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1
m.migrad()

print(m.params)

valsSIDIS=list(m.values)

chi_2SIDIS(m.values)

DataProcessor.harpyInterface.PrintChi2Table(setSIDIS,printDecomposedChi2=True)

print([round(x,1 if x >100 else 4 if x>1 else 6) for x in list(m.values)])

#%%
sys.exit()
