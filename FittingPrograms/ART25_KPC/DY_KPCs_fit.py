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

path_to_constants=MAINPATH+"FittingPrograms/ART25_KPC/ConstantsFiles/ART25_KPC.atmde"

harpy.initialize(path_to_constants)

# inARRAY_TMDR=[1.5, 0.071624, 0.064726, 0.0]
# inARRAY_PDF=[0.521462, 0.000206, 0.402948, 7.0219, 1.0, 20.4051, 1.0, 0.000123, 1.1037, 0.660734, 0.0, 0.04]

# ART 25 parameters
inARRAY_TMDR=[1.5, 0.0859, 0.0303, 0.0]
inARRAY_PDF=[0.486, 0.041, 0.569, 0.15, 5.26, 21.1, 7.7, 0.16, 0.24, 0.07, 0.0, 0.0]

harpy.setNPparameters_TMDR(inARRAY_TMDR)
harpy.setNPparameters_uTMDPDF(inARRAY_PDF)


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
    
    return ((delta<0.25 and p["<qT>"]<10.) or (delta<0.25 and par/err*delta**2<1)) , p
   

#%%
### Loading the DY data set
theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY([
                            'CDF1', 'CDF2',
                            'D01', 'D02', 'D02m', 
                           # 'A7-00y10', 'A7-10y20','A7-20y24', 
                            'A8-00y04', 'A8-04y08', 'A8-08y12', 'A8-12y16', 'A8-16y20', 'A8-20y24', 
                            'A8-46Q66', 'A8-116Q150', 'A13-norm',
                            'CMS7', 'CMS8', 
                            'CMS13-00y04','CMS13-04y08','CMS13-08y12','CMS13-12y16','CMS13-16y24',
                            # 'CMS13_dQ_50to76',
                            'CMS13_dQ_106to170',
                            # 'CMS13_dQ_170to350',
                            # 'CMS13_dQ_350to1000',
                            'LHCb7', 
                            'LHCb8', 
                            'LHCb13_dy(2021)', 
                            'PHE200', 'STAR510', 
                            'E228-200', 'E228-300', 'E228-400', 
                            'E772', 'E605',
                            # 'D0run1-W', 'CDFrun1-W'
                          ]))

setDY=theData.CutData(cutFunc) 
#setDYfit=theData.CutData(cutFuncFORFIT) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')


#%%
### Load the data for A4
theDataA4=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY(["A8_A4_0y1","A8_A4_1y2","A8_A4_2y35"]))
setA4=theDataA4.CutData(cutFunc) 

print('Loaded ', setA4.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setA4.sets]), 'points.')

### Load the data for AUU
theDataAUU=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY(["A8_Auu_0y1","A8_Auu_1y2","A8_Auu_2y35"]))
setAUU=theDataAUU.CutData(cutFunc) 

print('Loaded ', setAUU.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setAUU.sets]), 'points.')

#%%

# harpy.setNPparameters([1.5004, 0.073018, 0.038048, 0.0, 0.521462, 0.000206, 0.402948, 7.0219, 1.0, 20.4051, 1.0, 0.000123, 1.1037, 0.660734, 0.0, 0.04])

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)

#%%

# harpy.setNPparameters([1.5004, 0.073018, 0.038048, 0.0, 0.521462, 0.000206, 0.402948, 7.0219, 1.0, 20.4051, 1.0, 0.000123, 1.1037, 0.660734, 0.0, 0.04])

# DataProcessor.harpyInterface.PrintChi2Table(setA4,printDecomposedChi2=True)


#%%
#######################################
# Minimisation
#######################################
import time

# penalty_index=[-7,-6,-5,-4,-3]

def chi2DY(x):
    
    startT=time.time()
    harpy.setNPparameters_TMDR([x[0],x[1],x[2],x[3]])
    harpy.setNPparameters_uTMDPDF([x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15]])

    print('np set =',["{:8.3f}".format(i) for i in x])        
   
    # chi2 computation for the unpolDY data sets
    YY=DataProcessor.harpyInterface.ComputeXSec(setDY)
    ccDY2,cc3=setDY.chi2(YY)
   
    # chi2 computation for the A4 DY_angular data set
    YA4=DataProcessor.harpyInterface.ComputeXSec(setA4)
    YAUU=DataProcessor.harpyInterface.ComputeXSec(setAUU)
    ccA4,cc4=setA4.chi2(numpy.array(YA4)/numpy.array(YAUU))

    endT=time.time()

    print(":-> unpol_DY =", ccDY2/setDY.numberOfPoints,
          " A4 =", ccA4/setA4.numberOfPoints,
          " Total_fit =", (ccDY2+ccA4)/(setDY.numberOfPoints+setA4.numberOfPoints),
          " time =", endT-startT)
    
    return ccDY2 + ccA4


#%%
#### Minimize DY
from iminuit import Minuit

#---- PDFbias-like row

# initialValues=([1.5, 0.0859, 0.0303, 0.0, 
#                 0.486, 0.041, 0.569, 0.15, 
#                 5.26, 21.1, 7.7, 0.16,
#                 0.24, 0.07, 0.0, 0.0
#     ])

initialValues=([1.5, 0.0369, 0.0150, 0.0, 
                0.486, 0.041, 0.569, 0.15, 
                1., 2., 1., 0.16,
                0.24, 0.07, 0.0, 0.0
    ])

# initialValues=([1.5004, 0.073018, 0.038048, 0.0, 
#                 0.521462, 0.000206, 0.402948, 7.0219, 
#                 1.0, 20.4051, 1.0, 0.000123, 
#                 1.1037, 0.660734, 0.0, 0.04
#     ])

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

m = Minuit(chi2DY, initialValues)

m.errors=initialErrors
m.limits=searchLimits
m.fixed=parametersToMinimize
m.errordef=1

print(m.params)

chi2DY(m.values)

m.tol=0.0001*setDY.numberOfPoints*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1
m.migrad()

#%%

print(m.params)

valsDY=list(m.values)

chi2DY(m.values)

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)

print([round(x,1 if x >100 else 4 if x>1 else 6) for x in list(m.values)])

#%%
sys.exit()
