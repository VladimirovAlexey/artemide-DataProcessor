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

ATMDE_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide-development/"
HARPY_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide-development/harpy/"


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

#path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/ART25_main.atmde"
path_to_constants=MAINPATH+"FittingPrograms/ART25_KPC/ConstantsFiles/ART_KPC.atmde"


harpy.initialize(path_to_constants)

inARRAY_TMDR=[1.5, 0.088, 0.059, 0.0]
inARRAY_uTMDPDF=[0.85, 0.39, 0.14, 0.015, 0.01, 0.01, 0.639, 9.5, 0.0, 0.0, 0.0, 0.0]

# inARRAY_TMDR=[1.5, 0.0859, 0.0303, 0.0]
# inARRAY_uTMDPDF=[0.486, 0.041, 0.569, 0.15, 5.26, 21.1, 7.7, 0.16, 0.24, 0.07, 0.0, 0.0]

harpy.setNPparameters_TMDR(inARRAY_TMDR)

harpy.setNPparameters_uTMDPDF(inARRAY_uTMDPDF)


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
        elif("A4" in name):
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataA+name+".csv")
        elif("Auu" in name):
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
# Loading the DY data set

theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY([
                          'CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                          #'A7-00y10', 'A7-10y20','A7-20y24',
                          'A8-00y04', 'A8-04y08', 'A8-08y12', 'A8-12y16', 'A8-16y20', 'A8-20y24', 
                          'A8-46Q66', 'A8-116Q150', 
                          'A13-norm',
                          'CMS7', 'CMS8', 
                          'CMS13-00y04','CMS13-04y08','CMS13-08y12','CMS13-12y16','CMS13-16y24',
                          #'CMS13_dQ_50to76',
                          'CMS13_dQ_106to170',
                          'CMS13_dQ_170to350',
                          'CMS13_dQ_350to1000',
                          'LHCb7', 'LHCb8', 'LHCb13_dy(2021)', 
                          'PHE200', 'STAR510', 
                          'E228-200', 'E228-300', 'E228-400', 
                          'E772',
                          'E605',
                          'D0run1-W','CDFrun1-W'
                          ]))

setDY=theData.CutData(cutFunc) 

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
#print('Loaded ', setDYfit.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDYfit.sets]), 'points.')


#%%
#Table of chi2 for initial values of ART25 for unpolarized 

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)

#sys.exit()
#%%
### Load the data for A4
theDataA4=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY([
    "A8_A4_0y1","A8_A4_1y2","A8_A4_2y35",
    "CMS8-A4-00y10","CMS8-A4-10y21"]))
setA4=theDataA4.CutData(cutFunc) 

print('Loaded ', setA4.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setA4.sets]), 'points.')

### Load the data for A4UU
theDataA4UU=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY(["A8_Auu_0y1","A8_Auu_1y2","A8_Auu_2y35"]))
setA4UU=theDataA4UU.CutData(cutFunc) 

print('Loaded ', setA4UU.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setA4UU.sets]), 'points.')

#%%
#Table of chi2 for initial values of ART25 for polarized angular coefficient A4

DataProcessor.harpyInterface.PrintChi2Table(setA4,norma=setA4UU,printDecomposedChi2=True)

#Values are not correct because A4 data is normalized with unpol data so we need to consider this in our chi2 function

#%%
#######################################
# Minimisation
#######################################
import time
totalN=setDY.numberOfPoints

#penalty_index=[-7,-6,-5,-4,-3]

def chi_2DY(x):
    startT=time.time()
    
    harpy.setNPparameters_TMDR([x[0],x[1],x[2],x[3]])
    harpy.setNPparameters_uTMDPDF([x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15]])
    
    print('np set =',["{:8.3f}".format(i) for i in x])        

    # chi2 computation for the unpolDY data sets
    YY=DataProcessor.harpyInterface.ComputeXSec(setDY)
    ccDY2,cc3=setDY.chi2(YY)
    
    # chi2 computation for the A4 DY_angular data set
    YA4=DataProcessor.harpyInterface.ComputeXSec(setA4)
    YA4UU=DataProcessor.harpyInterface.ComputeXSec(setA4UU)
    A4=numpy.array(YA4)/numpy.array(YA4UU)
    ccA4,cc4=setA4.chi2(A4)
    
    cc=(ccDY2+ccA4)/(setDY.numberOfPoints+setA4.numberOfPoints)
    
    endT=time.time()
    
    print(":-> unpol_DY =", ccDY2/setDY.numberOfPoints,
          " A4 =", ccA4/setA4.numberOfPoints,
          " Total_Chi2 =", cc,
          " time =", endT-startT)
    
    return ccDY2 + ccA4
    

#%%
#### Minimize DY
from iminuit import Minuit

#---- PDFbias-like row

# initialValues=([1.5, 0.0369, 0.015, 0.0, 
#                 0.486, 0.041, 0.569, 0.15,
#                 1.0, 2., 1.0, 0.16,
#                 0.24, 0.07, 0.0, 0.0])

# initialValues=([1.5, 0.092, 0.055, 0.0, 
#                 0.393, 0.103, 0.665, 0.999,
#                 7.191, 43.949, -2.783, 0.003,
#                 0.476, 0.169, 0.0, 0.0])

initialValues=([1.5, 0.088, 0.059, 0.0, 
                0.85, 0.39, 0.14, 0.015,
                 0.01, 0.01, 0.639, 9.5,
                 0.0, 0.0, 0.0, 0.0])

# initialErrors=(0.1,0.1,0.1,0.1,
#                 0.1,  1.0, 0.1,  1.0,
#                 0.1,  1.0, 0.1,  1.0,
#                 0.1,  1.0, 1.0,  1.0)

initialErrors=(0.1,0.1,0.1,0.1,
                0.1,  1.0, 0.1,  1.0,
                0.1,  0.1, 0.1,  1.0,
                0.1,  1.0, 1.0,  1.0)

# searchLimits=((1.0,2.5),(0.005,0.15) ,(0.0,0.1), (-5.,5.),
#               (0.0001,2.), (0.0001,2.),(0.0001,2.),(0.0001,2.),
#               (0.0,10.), (0.0001,50.),(-10.,10.),(0.0001,10.),
#               (0.0001,10.), (0.0001,10.),(0.,100.),(0.,100.))

searchLimits=((1.0,2.5),(0.005,0.15) ,(0.0,0.1), (-5.,5.),
              (0.0001,2.), (0.0001,2.),(0.0001,2.),(0.0001,2.),
              (0.0001,100.), (0.0001,100.),(0.0001,10.),(0.0001,100.),
              (0.0001,10.), (0.0001,10.),(0.,100.),(0.,100.))
              
# True= FIX
# parametersToMinimize=(True, False,False,True,
#                       False, False, False,False,
#                       True, False, True, False,
#                       False, False, True,True)

parametersToMinimize=(True, False,False,True,
                      False, False, False,False,
                      False, False, False, False,
                      True, True, True,True)

#%%

m = Minuit(chi_2DY, initialValues)

print(m.params)
chi_2DY(m.values)
#%%
m = Minuit(chi_2DY, initialValues)
m.errors=initialErrors
m.limits=searchLimits
m.fixed=parametersToMinimize
m.errordef=1
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
