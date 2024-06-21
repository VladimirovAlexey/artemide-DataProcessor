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

replicaFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/replicaDY.dat"
logFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/LOG.log"

import sys
import numpy
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)

### number of circles that this code will run
NumberOfReplicas=25


#%%
import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet
import DataProcessor.DataMultiSet

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/data/WorkingFiles/TMD/Fit_Notes/ART23/REPS/ART23_run2.rep")

MAINPATH=ROOT_DIR 
#%%
#######################################
#Initialize artemide
#######################################
import harpy

path_to_constants=MAINPATH+"FittingPrograms/KPC/INI/ART23_MSHT_N4LL.atmde"
#path_to_constants=MAINPATH+"FittingPrograms/ART23/ConstantsFiles/DYonly/ART23_JAM_NLL.atmde"


harpy.initialize(path_to_constants)

initializationArray=[0.253434, 0.253434,0.253434, 0.253434,
                     0.253434, 0.253434,0.253434, 0.253434,
                     0.253434, 0.253434, 0.1,  0.1]

harpy.setNPparameters_TMDR([1.584237, 0.048428,0.001,0.])

harpy.setNPparameters_uTMDPDF(initializationArray)

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
    
    #  for artemide v3.    
    # p["process"]=[p["process"][0],p["process"][2],1,1]
    if(len(p["process"])==4):
        if(p["process"][3]==3): p["process"][3]=2
    if(len(p["process"])==3):
        if(p["process"][2]==1): p["process"]=[p["process"][0],1,1,1]
        elif(p["process"][2]==2): p["process"]=[p["process"][0],1,-1,1]
        elif(p["process"][2]==3): print("ERROR1")
        elif(p["process"][2]==4): print("ERROR1")
        elif(p["process"][2]==5): p["process"]=[p["process"][0],1,1,2]
        elif(p["process"][2]==6): p["process"]=[p["process"][0],1,-1,2]
        elif(p["process"][2]==7): print("ERRORW")
        elif(p["process"][2]==8): print("ERRORW")
        elif(p["process"][2]==9): print("ERRORW")
        elif(p["process"][2]==10): print("ERRORW")
        elif(p["process"][2]==11): print("ERRORW")
        elif(p["process"][2]==12): print("ERRORW")
        elif(p["process"][2]==1001): p["process"]=[p["process"][0],1,1,101]
        elif(p["process"][2]==1002): p["process"]=[p["process"][0],1,1,102]
        else:
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
    
    #ART23
    return ((delta<0.25 and p["<qT>"]<10.) or (delta<0.25 and par/err*delta**2<1)) , p
    #return ((delta<0.25 and p["<qT>"]<10.)) , p
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
                          'LHCb7', 'LHCb8', 
                          'LHCb13_dy(2021)', 
                          'PHE200', 'STAR510', 
                          'E228-200', 'E228-300', 'E228-400', 
                          'E772','E605'
                          #'D0run1-W','CDFrun1-W'
                          ]))

setDY=theData.CutData(cutFunc) 
#setDYfit=theData.CutData(cutFuncFORFIT) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
#print('Loaded ', setDYfit.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDYfit.sets]), 'points.')
#%%
# MAIN FIT
#harpy.setNPparameters([1.5004, 0.05614, 0.03862, 0.0, 0.565, 0.0539, 0.5697, 6.64, 0., 20.07, 0., 0.537, 1.07, 2.39, 0.0, 0.0])
harpy.setNPparameters([1.500, 0.0565, 0.06086, 0.0, 0.4223, 0.16896, 0.56410, 8.387, 5.986, 16.026, 5.646, 0.3116, 1.0575, 0.01, 0.0, 0.0])
DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True,method="approximate")

#%%
#######################################
# Minimisation
#######################################
import time
totalN=setDY.numberOfPoints

### FOR ART23
#penalty_index=[-7,-6,-5,-4,-3]

penalty_index=[-1,-2,-3,-4,-5]

def chi_2DY(x):
    startT=time.time()
    harpy.setNPparameters_TMDR([x[0],x[1],x[2],x[3]])
    #harpy.setNPparameters_uTMDPDF(x[4:])
    harpy.setNPparameters_uTMDPDF([x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15]])
    # harpy.setNPparameters_TMDR([x[0],x[1],x[2],x[3]])
    # harpy.setNPparameters_uTMDPDF(x[4:])
    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")        
    
    #YY=DataProcessor.harpyInterface.ComputeXSec(setDY)
    YY=DataProcessor.harpyInterface.ComputeXSec(setDY,method="approximate")
    #YY=RAT*numpy.array(YY0)
    ccDY2,cc3=setDY.chi2(YY)
    
    penalty_array=numpy.array([max(0,abs(setDY.sets[i].DetermineAvarageSystematicShift(YY[setDY._i1[i]:setDY._i2[i]]))/setDY.sets[i].normErr[0]-1) for i in penalty_index])
    penalty_term=sum(penalty_array**6)*10
    
    #ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)
    #ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(setSIDIS)
    
    cc=ccDY2/setDY.numberOfPoints
    endT=time.time()
    print(':->',cc,'   +p=',penalty_term/setDY.numberOfPoints,'       t=',endT-startT)
    return ccDY2+penalty_term
#%%
# chi2=1.163 (muOPE->2)
# [1.5, 0.049098, 0.05979, 0.0, 0.364806, 4.7486, 0.469411, 0.652846, 0.0, 33.8383, 0.0, 24.3613, 1.4598, 0.054165, 0.0, 0.0]

# chi2=1.148 (muOPE->5)
# [1.5, 0.049098, 0.05979, 0.0, 0.696635, 3.7706, 0.477978, 0.740207, 0.0, 35.6754, 0.0, 22.7604, 1.4951, 0.020935, 0.0, 0.0]
#%%
#### Minimize DY
from iminuit import Minuit

#---- PDFbias-like row
initialValues=([1.5004, 0.05614, 0.03862, 0.0, 
                0.565, 0.0539, 0.5697, 6.64, 
                0., 20.07, 0., 0.537, 
                1.07,2.39, 0.0, 0.0
    ])

# initialValues=([1.56142, 0.0369174, 0.0581734, 1.0,
#   0.874245, 0.913883, 0.991563, 6.05412,
#   0.353908, 46.6064, 0.115161, 1.53235,
#   1.31966, 0.434833, 0.0, 0.0])
initialErrors=(0.1,0.1,0.1,0.1,
                0.1,  1.0, 0.1,  1.0,
                0.1,  1.0, 0.1,  1.0,
                0.1,  1.0, 10.,  1.)
searchLimits=((1.0,2.5),(0.004,0.08) ,(0.0,.07), (-5.,5.),
              (0.0001,100.), (0.0001,100.),(0.0001,100.),(0.0001,100.),
              (-100.,100.), (0.0001,100.),(-100.,100.),(0.0001,100.),
              (0.0001,100.), (0.0001,100.),(0.,100.),(0.0001,100.))
              
# True= FIX
parametersToMinimize=(True, False,False,True,
                      False, False, False,False,
                      False, False, False,False,
                      False, False, False,True)

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

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True,method="approximate")

print([round(x,1 if x >100 else 4 if x>1 else 6) for x in list(m.values)])


DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)

#%%
# ### Making plots and saving to files...
# harpy.setNPparameters([1.500, 0.056, 0.061, 0.000, 0.422, 0.169, 0.564, 8.387, 5.986, 16.026, 5.646, 0.312, 1.057, 0.010, 0.000, 0.000])
# pathToPlot="/data/WorkingFiles/TMD/Fit_Notes/KPC/PlotsData/KPC_ART23+/"
# import time
# for s in setDY.sets:
#     startT=time.time()
#     XX=DataProcessor.harpyInterface.ComputeXSec(s,method="approximate")
#     endT=time.time()
#     print(':->',s.name,'       t=',endT-startT)
#     f=open(pathToPlot+s.name+".dat","w")
#     print('SAVING PLOTS>>  ',f.name)
#     ### [total chi^2, list of NP-parameters],
#     for i in range(len(XX)):
#         p=s.points[i]
#         f.write(str(p["qT"][0])+", "+str(p["qT"][1])+", "+
#                 str(p["xSec"])+", "+str(numpy.sqrt(sum(numpy.array(p["uncorrErr"])**2)))+", "+
#                 str(XX[i])+"\n")
        
#     f.close()

#%%

# def ChangeQ(p):
#     p["Q"]=[9.0,10.0]
#     return True,p


# setEE=setDY.sets[-5].CutData(ChangeQ) 