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
import DataProcessor.ArtemideReplicaSet
import DataProcessor.DataMultiSet

MAINPATH=ROOT_DIR 
#%%
#######################################
#Initialize artemide
#######################################
import harpy


#path_to_constants=MAINPATH+"FittingPrograms/MVZ22/ConstantsFiles/DYonly/ART23_MSHT_N4LL"
path_to_constants=MAINPATH+"FittingPrograms/MVZ22/ConstantsFiles/DYonly/ART23_NNPDF_N4LL"

harpy.initialize(path_to_constants)

#harpy.setNPparameters([1.584237, 0.048428, 0.521983, 5.867221, 406.015479, 2.542726, -6.352752, 0.0, 0.1])

harpy.setNPparameters_TMDR([1.584237, 0.048428,0.001,0.])

harpy.setNPparameters_uTMDPDF([0.253434, 0.253434,
                                0.253434, 0.253434,
                                0.253434, 0.253434,
                                0.253434, 0.253434,
                                0.253434, 0.253434,
                                346.999,  0.1])

#%%
### read the list of files and return the list of DataSets
def loadThisDataDY(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=ROOT_DIR+"DataLib/unpolDY/"
    path_to_dataW=ROOT_DIR+"DataLib/unpolW/"
    
    
    dataCollection=[]
    for name in listOfNames:
        if(name[-1]=="W"):
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataW+name+".csv")
        else:
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        

    return dataCollection



#%%
##################Cut function
def cutFunc(p):
    par=0.5
    
    if(p["xSec"]>0):
        err=numpy.sqrt(sum([i**2 for i in p["uncorrErr"]]))/p["xSec"]
    else:
        err=100.
    delta=p["<qT>"]/p["<Q>"]
    
    if(p["id"][0] == "E"):
        delta=p["<qT>"]/p["Q"][1] 
    
    if("run1-W" in p["id"]):
        delta=p["qT"][0]/(p["Q"][0]+5.)
    
    
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
    
#    return delta<0.5 and p.qT_avarage<80
    return ((delta<0.25 and p["<qT>"]<10.) or (delta<0.25 and par/err*delta**2<1)) , p

#%%
##################Cut function
def cutFuncFORFIT(p):
    par=0.5
    
    if("CDF1" in p['id']): return False, p
    if("D01" in p['id']): return False, p
    if("D02" in p['id']): return False, p
    if("CMS7" in p['id']): return False, p
    if("CMS8" in p['id']): return False, p
    if("CDFrun1-W" in p['id']): return False, p
    
    if(p["xSec"]>0):
        err=numpy.sqrt(sum([i**2 for i in p["uncorrErr"]]))/p["xSec"]
    else:
        err=100.
    delta=p["<qT>"]/p["<Q>"]
    
    if(p["id"][0] == "E"):
        delta=p["<qT>"]/p["Q"][1] 
    
    if("run1-W" in p["id"]):
        delta=p["qT"][0]/(p["Q"][0]+5.)
    
    
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
    
#    return delta<0.5 and p.qT_avarage<80
    return ((delta<0.25 and p["<qT>"]<15.) or (delta<0.25 and par/err*delta**2<0.)) , p

#%%
### Loading the DY data set
theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY([
                          'CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                          #'A7-00y10', 'A7-10y20','A7-20y24', 
                          'A8-00y04', 'A8-04y08', 'A8-08y12', 'A8-12y16', 'A8-16y20', 'A8-20y24', 
                          'A8-46Q66', 'A8-116Q150', 
                          'A13-norm',
                          'CMS7', 'CMS8', 
                          'CMS13-00y04','CMS13-04y08','CMS13-08y12','CMS13-12y16','CMS13-16y24',
                          #'CMS13_dQ_50to76',
                          'CMS13_dQ_106to170','CMS13_dQ_170to350','CMS13_dQ_350to1000',
                          'LHCb7', 'LHCb8', 'LHCb13_dy(2021)', 
                          'PHE200', 'STAR510', 
                          'E228-200', 'E228-300', 'E228-400', 
                          'E772',
                          'E605',
                          'D0run1-W','CDFrun1-W'
                          ]))

setDY=theData.CutData(cutFunc) 
setDYfit=theData.CutData(cutFuncFORFIT) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded ', setDYfit.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDYfit.sets]), 'points.')


#%%

#(MSHT)
#harpy.setNPparameters([1.4806, 0.038969, 0.051737, 1.0, 
#                       0.851645, 0.69432, 0.934676, 5.2514, 
#                       0.247602, 39.6123, 0.094435, 1.8872, 
#                       1.2164, 0.936465, 0.0, 0.0])

#(NNPDF)
harpy.setNPparameters([1.9608, 0.051636, 0.065776, 1.0, 0.950443, 0.542395, 0.96724, 4.0706, 0.402134, 43.4538, 0.198087, 1.2893, 1.0357, 0.039715, 0.0, 0.0])
DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)




#%%
#######################################
# Minimisation
#######################################
import time
totalN=setDY.numberOfPoints

def chi_2DY(x):
    startT=time.time()
    harpy.setNPparameters_TMDR([x[0],x[1],x[2],x[3]])
    harpy.setNPparameters_uTMDPDF(x[4:])

    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
    
    ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)
    #ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(setSIDIS)
    
    cc=ccDY2/setDY.numberOfPoints
    endT=time.time()
    print(':->',cc,'       t=',endT-startT)
    return ccDY2

#%%
#### Minimize DY
from iminuit import Minuit

#---- PDFbias-like row

initialValues=([1.5004, 0.049098, 0.05979, 0.0, 
                0.834557, 0.914247, 0.910747, 4.5973, 
                0.4487, 3.85017, 0.1313, 1.2705, 
                1.1989, 0.173397, 0.0, 0.0])
initialErrors=(0.1,0.1,0.1,0.1,
                0.1,  1.0, 0.1,  1.0,
                0.1,  1.0, 0.1,  1.0,
                0.1,  1.0, 10.,  1.)
searchLimits=((0.2,2.5),(0.0,1.) ,(0.,1.), (-5.,5.),
              (0.,10.), (0.,100.),(0.,10.), (0.,100.),
              (0.,10.), (0.,100.),(0.,10.), (0.,100.),
              (0.,10.), (0.,100.),(0,5000.), (0,100.))
# True= FIX
parametersToMinimize=(False, False,False,True,
                      #True,True,True,True,
                      #True,True,True,True,
                      #True,True,True,True)
                      False, False, False,False,
                      False, False, False,False,
                      False, False, True,True)




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

chi_2DY(m.values)

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)

print([round(x,1 if x >100 else 4 if x>1 else 6) for x in list(m.values)])
