#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 16:40:59 2019

@author: vla18041
"""

#######################################
# importing libraries
#######################################

import sys
import time
import numpy
sys.path.append("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor")
import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet

MAINPATH="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"
#%%
#######################################
#Initialize artemide with a replica -2
#######################################
import harpy
path_to_constants=MAINPATH+"FittingPrograms/PDF-TMD/Constants-files/SV19/"
harpy.initialize(path_to_constants+"SV19-NNPDF31_NNLO")
#harpy.initialize(path_to_constants+"SV19-CT18_NNLO")
#harpy.initialize(path_to_constants+"SV19-MSHT20_NNLO")
harpy.setNPparameters([2.0340, 0.0299, 0.2512, 7.7572,334.6108, 2.4543,-4.8203, 0.1000,  0.0000])
#harpy.setNPparameters_TMDR(-2)
#harpy.setNPparameters_uTMDPDF(-2)

#%%
### read the list of files and return the list of DataSets
def loadThisData(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/DataLib/unpolDY/"
    
    
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
       
    
#    return delta<0.5 and p.qT_avarage<80
    return (delta<0.1 or (delta<0.25 and par/err*delta**2<1)) , p

#%%
### Loading the data set
theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisData(['CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                      'A7-00y10', 'A7-10y20','A7-20y24', 
                      'A8-00y04', 'A8-04y08', 'A8-08y12', 'A8-12y16', 'A8-16y20', 'A8-20y24', 
                      'A8-46Q66', 'A8-116Q150', 
                      'CMS7', 'CMS8',
                      'LHCb7', 'LHCb8', 'LHCb13', 
                      'PHE200', 'E228-200', 'E228-300', 'E228-400', 
                      'E772',
                      'E605']))

setDY=theData.CutData(cutFunc) 

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setDY.sets])
#%%
#######################################
# Test the main replica
#######################################
import DataProcessor.ArtemideReplicaSet
rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/Replicas/"+
                                                  "DY_nnlo/DY_NNPDF31_nnlo.rep")
                                                  # "Sivers20_model9case1(noDY-n3lo).rep")

rSet.SetReplica()

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)
    
#%%
#######################################
# Minimisation
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

from iminuit import Minuit

#initialValues=(2.29477, 0.022191,0.324112, 13.1774, 356.124, 2.05101, -10.4468, 0., 0.)  ### HERA

initialValues=(1.89278, 0.0291215,0.273023, 8.96853, 353.824, 2.50639, -5.92725, 0., 0.)  ### NNPDF

#initialValues=(1.54728,  0.04678,0.1982,   26.49689,2727.9766,    3.00668, -23.54749,   0.,0.) ##MMHT14
#initialValues=(2.34665, 0.022735,0.27706, 24.8883, 1241.34, 2.66869, -23.7589, 0., 0.) ##CT14
#initialValues=(1.92659, 0.036548,0.218079, 17.9138, 926.078, 2.5431, -15.5469, 0. ,0.) ##PDF4LHC


initialErrors=(0.1,       0.01,     0.05,      0.1,      10.,       0.05,    0.05,    0.1, 0.1)
searchLimits=((1.4,4.5), (0.001,5.0),(0.0,4.0),(8.,25.0),(0.0,400),(0.,5), None,   None, None)
parametersToMinimize=(False,     False,    False,    False,    False,     False,  False, True, True)


#%%
m = Minuit(chi_2, initialValues)

m.errors=initialErrors
m.limits=searchLimits
m.fixed=parametersToMinimize
m.errordef=1

print(m.params)

m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1
m.migrad()

## print parameters
print(m.params)

## print chi^2 table
harpy.setNPparameters(list(m.values))

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)
