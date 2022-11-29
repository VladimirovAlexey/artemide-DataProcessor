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

#path_to_constants=MAINPATH+"FittingPrograms/MVZ22/ConstantsFiles/DYonly/PDFb_DYonly_N3LL+N3LO"
path_to_constants=MAINPATH+"FittingPrograms/MVZ22/ConstantsFiles/DYonly/PDFb_DYonly_N3LL+N4LO"

#path_to_constants=MAINPATH+"FittingPrograms/MVZ22/ConstantsFiles/DYonly/SV22_DYonly_N4LL"
#path_to_constants=MAINPATH+"FittingPrograms/MVZ22/ConstantsFiles/DYonly/SV19_DYonly_N4LL"

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
    return (delta<0.10 or (delta<0.25 and par/err*delta**2<1)) , p

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
                          #'CMS13_dQ_106to170','CMS13_dQ_170to350','CMS13_dQ_350to1000',
                          'LHCb7', 'LHCb8', 'LHCb13_dy(2021)', 
                          'PHE200', 'STAR510', 
                          'E228-200', 'E228-300', 'E228-400', 
                          'E772',
                          'E605']))

setDY=theData.CutData(cutFunc) 


print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setDY.sets])

#%%

# harpy.setNPparameters([[1.584237, 0.048428, 0.521983, 5.867221, 406.015479, 2.542726, -6.352752, 0.0, 0.0]])

# #(original)
#harpy.setNPparameters([2., 0.0436918, 0.103944, 1.65603, 0.36681, 1.24648, 0.00360835, 57.5154, 0.00121168, 0.45159, 0.108631, 4.89878, 22.7981, 0])
#(1)
#harpy.setNPparameters([2.2893,2.2893, 0.044942,0., 0.054509, 1.3851, 0.12001, 1.3159, 0.0, 32.5008, 0.0, 0.044856, 0.027146, 3.1409, 18.4962,  0.0])
#(8)
#harpy.setNPparameters([1.7, 0.065702, 0.066928, 0., 0.10833, 2.7129, 0.611321, 0.333934, 0.0, 54.434, 0.0, 0.001766, 0.046555, 0.005283, 31.2679, 0.0])
#(9)
#harpy.setNPparameters([1.4264, 0.053155, 0.0654, 0.0, 0.11059, 2.9452, 0.659106, 0.991233, 0.000456, 76.373, 0., 0.298568, 0.105762, 1.1453, 19.4772, 0.0])
#(10)
#harpy.setNPparameters([1.4312, 0.045097, 0.065787, 0.0, 0.088258, 2.6835, 1.011, 0.946874, 0.000374, 64.1763, 2.4e-05, 0.05902, 0.178093, 1.2986, 11.1923, 0.0])
#(10+A13)
#harpy.setNPparameters([1.4094, 0.041582, 0.059136, 0.0, 0.117864, 1.9074, 0.964391, 1.2101, 0.039774, 39.6539, 2e-06, 0.000151, 0.346594, 0.539208, 5.7965, 0.0])
#(11)
harpy.setNPparameters([1.4107, 0.044607, 0.071958, 0.0, 0.920405, 1.2183, 1.0494, 4.095, 0.422448, 32.4762, 2e-06, 0.824624, 1.278, 2.7997, 11.1923, 0.0])


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
    # harpy.setNPparameters_TMDR([x[0],x[1],x[2],x[3]])
    # harpy.setNPparameters_uTMDPDF(x[4:])
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
initialValues=([1.4312, 0.045097, 0.065787, 0.0,
                1.0823, 1.5445, 1.0018, 0.002023,
                0.51518, 33.5879, 1e-06, 0.013953, 
                1.0528, 2.6348, 11.1923, 0.])
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
                      False, False, False,True)




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
