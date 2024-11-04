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

#path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/ART23_MSHT_N4LL_v301_resum.atmde"
path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/ART23_MSHT_N4LL_v301.atmde"


harpy.initialize(path_to_constants)

initializationArray=[0.253434, 0.253434,0.253434, 0.253434,
                     0.253434, 0.253434,0.253434, 0.253434,
                     0.253434, 0.253434, 0.1,  0.04]

harpy.setNPparameters_TMDR([1.584237, 0.048428,0.001,0.])

harpy.setNPparameters_uTMDPDF(initializationArray)

#%%
### read the list of files and return the list of DataSets
def loadThisDataDY(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=ROOT_DIR+"DataLib/unpolDY/"
    path_to_dataW=ROOT_DIR+"DataLib/unpolW/"
    path_to_dataA=ROOT_DIR+"DataLib/DY_angular/"
    
    
    dataCollection=[]
    for name in listOfNames:
        if(name[-1]=="W"):
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataW+name+".csv")
        elif("_A4" in name):
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataA+name+".csv")
        elif("_Auu" in name):
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataA+name+".csv")
        else:
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        

    return dataCollection
#%%
##################Cut function
def cutFunc(p):
    par=0.5
    
    #  for artemide v3.    
    # p["process"]=[p["process"][0],p["process"][2],1,1]
    if(len(p["process"])==3):        
            print("UNKNOWN PROCESS IN ARTEMIDE 3"+str(p["process"]))
    
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
    #return (delta<0.25) , p

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
#setDYfit=theData.CutData(cutFuncFORFIT) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
#print('Loaded ', setDYfit.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDYfit.sets]), 'points.')

#%%
### Load the data for A4
theDataA4=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY(["A8_A4_0y1","A8_A4_1y2","A8_A4_2y35"]))
setA4=theDataA4.CutData(cutFunc) 

### Load the data for Auu
theDataAuu=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY(["A8_Auu_0y1","A8_Auu_1y2","A8_Auu_2y35"]))
setAuu=theDataAuu.CutData(cutFunc) 

#%%
#(ART23+alpha=0.04)
#harpy.setNPparameters([1.5004, 0.060236, 0.037013, 0.0, 0.602249, 0.000103, 0.445843, 7.2217, 1.0, 18.6931, 1.0, 0.000101, 1.1114, 1.3137, 0.0, 0.04])
#(ART23+alpha=0.04  + resum)
harpy.setNPparameters([1.5004, 0.073018, 0.038048, 0.0, 0.521462, 0.000206, 0.402948, 7.0219, 1.0, 20.4051, 1.0, 0.000123, 1.1037, 0.660734, 0.0, 0.04])

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)


#%%
##### by Vmoos
#harpy.setNPparameters_TMDR([1.5004, 0.05614, 0.03862, 1.0])
#harpy.setNPparameters_uTMDPDF([0.565, 0.0539, 0.5697, 6.64, 0.565, 20.07, 0.5697, 0.537, 1.07, 2.39, 0.0, 0.0])
#DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)
#%%
#######################################
# Minimisation
#######################################
import time
totalN=setDY.numberOfPoints

penalty_index=[-7,-6,-5,-4,-3]

def chi_2DY(x):
    startT=time.time()
    harpy.setNPparameters_TMDR([x[0],x[1],x[2],x[3]])
    #harpy.setNPparameters_uTMDPDF(x[4:])
    harpy.setNPparameters_uTMDPDF([x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15]])
    #harpy.varyScales(1.,1.,1.,2.0)
    #harpy.varyScales(1.,1.,1.,1.0)
    # harpy.setNPparameters_TMDR([x[0],x[1],x[2],x[3]])
    # harpy.setNPparameters_uTMDPDF(x[4:])
    print('np set =',["{:8.3f}".format(i) for i in x])        
    
    YY=DataProcessor.harpyInterface.ComputeXSec(setDY)
    ccDY2,cc3=setDY.chi2(YY)
    
    penalty_array=numpy.array([max(0,abs(setDY.sets[i].DetermineAvarageSystematicShift(YY[setDY._i1[i]:setDY._i2[i]]))/setDY.sets[i].normErr[0]-1) for i in penalty_index])
    penalty_term=sum(penalty_array**6)
    
    YA4=DataProcessor.harpyInterface.ComputeXSec(setA4)
    YAuu=DataProcessor.harpyInterface.ComputeXSec(setAuu)
    ccA4,cc4=setA4.chi2(-numpy.array(YA4)/numpy.array(YAuu))
    
    #ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)
    #ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(setSIDIS)
    
    #penalty_term=0.
    #ccDY2=0.
    
    cc=ccDY2/setDY.numberOfPoints
    ccR=ccA4/setA4.numberOfPoints
    
    
    endT=time.time()
    print(':->',cc,'  ',ccR,'   +p=',penalty_term,"    time=",endT-startT)
    return ccDY2+ccA4+penalty_term

#%%
# chi2=1.163 (muOPE->2)
# [1.5, 0.049098, 0.05979, 0.0, 0.364806, 4.7486, 0.469411, 0.652846, 0.0, 33.8383, 0.0, 24.3613, 1.4598, 0.054165, 0.0, 0.0]

# chi2=1.148 (muOPE->5)
# [1.5, 0.049098, 0.05979, 0.0, 0.696635, 3.7706, 0.477978, 0.740207, 0.0, 35.6754, 0.0, 22.7604, 1.4951, 0.020935, 0.0, 0.0]
#%%
#### Minimize DY
from iminuit import Minuit

#---- PDFbias-like row
initialValues=([1.5004, 0.073018, 0.038048, 0.0, 
                0.521462, 0.000206, 0.402948, 7.0219, 
                1.0, 20.4051, 1.0, 0.000123, 
                1.1037, 0.660734, 0.0, 0.04
    ])

initialErrors=(0.1,0.1,0.1,0.1,
                0.1,  1.0, 0.1,  1.0,
                0.1,  1.0, 0.1,  1.0,
                0.1,  1.0, 10.,  1.)
searchLimits=((1.0,2.5),(0.005,0.15) ,(0.0,.2), (-5.,5.),
              (0.0001,100.), (0.0001,100.),(0.0001,100.),(0.0001,100.),
              (-100.,100.), (0.0001,100.),(-100.,100.),(0.0001,100.),
              (0.0001,100.), (0.0001,100.),(0.,100.),(0.04,100.))
              
# True= FIX
parametersToMinimize=(True, False,False,True,
                      False, False, False,False,
                      True, False, True, False,
                      False, False, True,True)

#%%

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

#%%
sys.exit()
