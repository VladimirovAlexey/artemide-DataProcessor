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
### read the list of files and return the list of DataSets
def loadThisDataDY_angular(listOfNames):    
    import DataProcessor.DataSet
    path_to_dataA=ROOT_DIR+"DataLib/DY_angular/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataA+name+".csv")
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
    #return ((delta<0.25 and p["<qT>"]<10.) or (delta<0.25 and par/err*delta**2<1)) , p
    #return ((delta<0.25 and p["<qT>"]<10.)) , p
    return (delta<0.25) , p

#%%
### Loading the DY data set
theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY([
                          #'CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                          #'A7-00y10', 'A7-10y20','A7-20y24', 
                          'A8-00y04', 'A8-04y08', 'A8-08y12', 'A8-12y16', 'A8-16y20', 'A8-20y24', 
                          'A8-46Q66', 'A8-116Q150', 
                          'A13-norm',
                          #'CMS7', 'CMS8', 
                          'CMS13-00y04','CMS13-04y08','CMS13-08y12','CMS13-12y16','CMS13-16y24',
                          #'CMS13_dQ_50to76',
                          #'CMS13_dQ_106to170','CMS13_dQ_170to350','CMS13_dQ_350to1000',
                          'LHCb7', 'LHCb8', 
                          'LHCb13_dy(2021)', 
                          #'PHE200', 'STAR510', 
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
theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY([
                          'A8-00y04', 'A8-04y08', 'A8-08y12', 'A8-12y16', 'A8-16y20', 'A8-20y24']))

setA8=theData.CutData(cutFunc) 

#%%
### Load the data for A4
theDataA4=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_A4_0y1","A8_A4_1y2","A8_A4_2y35"]))
setA4=theDataA4.CutData(cutFunc) 

### Load the data for Auu
theDataAuu=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY_angular(["A8_Auu_0y1","A8_Auu_1y2","A8_Auu_2y35"]))
setAuu=theDataAuu.CutData(cutFunc) 

#%%
# MAIN FIT
harpy.setNPparameters([1.5004, 0.05614, 0.03862, 0.0, 0.565, 0.0539, 0.5697, 6.64, 0.565, 20.07, 0.5697, 0.537, 1.07, 2.39, 0.0, 0.0])
#%% 
#FIT with cut facto from KPC (fast chi2=1.055)
harpy.setNPparameters([1.5, 0.053509, 0.037095, 0.0, 0.549212, 0.172144, 0.538113, 6.4105, 0.565, 19.8866, 0.5697, 0.638572, 1.0974, 10.0389, 0.004768, 0.0])

#%%
#### Computation of different contributions to the unpolarized fiducial cross-section
XXN=DataProcessor.harpyInterface.ComputeXSec(setA8)

def ToA0(p):
    p["process"][3]=20
    return True, p

setA8_A0=setA8.CutData(ToA0)

XX0=DataProcessor.harpyInterface.ComputeXSec(setA8_A0)
#%%
def ToA1(p):
    p["process"][3]=21
    return True, p

setA8_A1=setA8.CutData(ToA1)

XX1=DataProcessor.harpyInterface.ComputeXSec(setA8_A1)

def ToA2(p):
    p["process"][3]=22
    return True, p

setA8_A2=setA8.CutData(ToA2)

XX2=DataProcessor.harpyInterface.ComputeXSec(setA8_A2)

#%%
XXuu=DataProcessor.harpyInterface.ComputeXSec(setAuu)
XX4=DataProcessor.harpyInterface.ComputeXSec(setA4)

#%%
### Preparing data set for angular coefficients

print("Done.  =>     Create points & append to data set ...")
DataA4=DataProcessor.DataSet.DataSet('A8-A4',"DY")
DataA4.comment="ATLAS8"
DataA4.reference="1606.00689"

DataA4.isNormalized=False
proc_current=[1,1,1,24]
s_current=64000000.0
Q_current=[80.,100.]
y_current=[0.,.1]
lumUncertainty=0.039
DataA4.normErr.append(lumUncertainty)

ptBINS=[0.,2.5,5.0,8.0,11.4,14.9,18.5,22.0,25.5,29.0]
yBINS=[0.,1.0,2.0,3.5]

for i in range(len(ptBINS)-1):
    for j in range(len(yBINS)-1):
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataA4.name+'.'+str(i))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=[ptBINS[i],ptBINS[i+1]]
        p["thFactor"]=1.
        p["Q"]=Q_current
        p["<Q>"]=91.2 ### to be sure that mass middle value is Z-boson mass
        p["y"]=[yBINS[j],yBINS[j+1]]
        p["includeCuts"]=False
        p["xSec"]=0.1
        p["uncorrErr"].append(0.1)
        #
        DataA4.AddPoint(p)

print("Done.  ")

print("Done.  =>     Create points & append to data set ...")
DataA3=DataProcessor.DataSet.DataSet('A8-A3',"DY")
DataA3.comment="ATLAS8"
DataA3.reference="1606.00689"

DataA3.isNormalized=False
proc_current=[1,1,1,23]
s_current=64000000.0
Q_current=[80.,100.]
y_current=[0.,.1]
lumUncertainty=0.039
DataA3.normErr.append(lumUncertainty)

ptBINS=[0.,2.5,5.0,8.0,11.4,14.9,18.5,22.0,25.5,29.0]
yBINS=[0.,1.0,2.0,3.5]

for i in range(len(ptBINS)-1):
    for j in range(len(yBINS)-1):
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataA3.name+'.'+str(i))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=[ptBINS[i],ptBINS[i+1]]
        p["thFactor"]=1.
        p["Q"]=Q_current
        p["<Q>"]=91.2 ### to be sure that mass middle value is Z-boson mass
        p["y"]=[yBINS[j],yBINS[j+1]]
        p["includeCuts"]=False
        p["xSec"]=0.1
        p["uncorrErr"].append(0.1)
        #
        DataA3.AddPoint(p)

print("Done.  ")

print("Done.  =>     Create points & append to data set ...")
DataA0=DataProcessor.DataSet.DataSet('A8-A0',"DY")
DataA0.comment="ATLAS8"
DataA0.reference="1606.00689"

DataA0.isNormalized=False
proc_current=[1,1,1,20]
s_current=64000000.0
Q_current=[80.,100.]
y_current=[0.,.1]
lumUncertainty=0.039
DataA3.normErr.append(lumUncertainty)

ptBINS=[0.,2.5,5.0,8.0,11.4,14.9,18.5,22.0,25.5,29.0]
yBINS=[0.,1.0,2.0,3.5]

for i in range(len(ptBINS)-1):
    for j in range(len(yBINS)-1):
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataA0.name+'.'+str(i))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=[ptBINS[i],ptBINS[i+1]]
        p["thFactor"]=1.
        p["Q"]=Q_current
        p["<Q>"]=91.2 ### to be sure that mass middle value is Z-boson mass
        p["y"]=[yBINS[j],yBINS[j+1]]
        p["includeCuts"]=False
        p["xSec"]=0.1
        p["uncorrErr"].append(0.1)
        #
        DataA0.AddPoint(p)

print("Done.  ")

print("Done.  =>     Create points & append to data set ...")
DataA1=DataProcessor.DataSet.DataSet('A8-A3',"DY")
DataA1.comment="ATLAS8"
DataA1.reference="1606.00689"

DataA1.isNormalized=False
proc_current=[1,1,1,21]
s_current=64000000.0
Q_current=[80.,100.]
y_current=[0.,.1]
lumUncertainty=0.039
DataA1.normErr.append(lumUncertainty)

ptBINS=[0.,2.5,5.0,8.0,11.4,14.9,18.5,22.0,25.5,29.0]
yBINS=[0.,1.0,2.0,3.5]

for i in range(len(ptBINS)-1):
    for j in range(len(yBINS)-1):
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataA1.name+'.'+str(i))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=[ptBINS[i],ptBINS[i+1]]
        p["thFactor"]=1.
        p["Q"]=Q_current
        p["<Q>"]=91.2 ### to be sure that mass middle value is Z-boson mass
        p["y"]=[yBINS[j],yBINS[j+1]]
        p["includeCuts"]=False
        p["xSec"]=0.1
        p["uncorrErr"].append(0.1)
        #
        DataA1.AddPoint(p)

print("Done.  ")

print("Done.  =>     Create points & append to data set ...")
DataA2=DataProcessor.DataSet.DataSet('A8-A2',"DY")
DataA2.comment="ATLAS8"
DataA2.reference="1606.00689"

DataA2.isNormalized=False
proc_current=[1,1,1,23]
s_current=64000000.0
Q_current=[80.,100.]
y_current=[0.,.1]
lumUncertainty=0.039
DataA2.normErr.append(lumUncertainty)

ptBINS=[0.,2.5,5.0,8.0,11.4,14.9,18.5,22.0,25.5,29.0]
yBINS=[0.,1.0,2.0,3.5]

for i in range(len(ptBINS)-1):
    for j in range(len(yBINS)-1):
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataA2.name+'.'+str(i))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=[ptBINS[i],ptBINS[i+1]]
        p["thFactor"]=1.
        p["Q"]=Q_current
        p["<Q>"]=91.2 ### to be sure that mass middle value is Z-boson mass
        p["y"]=[yBINS[j],yBINS[j+1]]
        p["includeCuts"]=False
        p["xSec"]=0.1
        p["uncorrErr"].append(0.1)
        #
        DataA2.AddPoint(p)

print("Done.  ")

print("Done.  =>     Create points & append to data set ...")
DataAN=DataProcessor.DataSet.DataSet('A8-AN',"DY")
DataAN.comment="ATLAS8"
DataAN.reference="1606.00689"

DataAN.isNormalized=False
proc_current=[1,1,1,2]
DataAN.normErr.append(lumUncertainty)


for i in range(len(ptBINS)-1):
    for j in range(len(yBINS)-1):
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataAN.name+'.'+str(i))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=[ptBINS[i],ptBINS[i+1]]
        p["thFactor"]=1.
        p["Q"]=Q_current
        p["<Q>"]=91.2 ### to be sure that mass middle value is Z-boson mass
        p["y"]=[yBINS[j],yBINS[j+1]]
        p["includeCuts"]=False
        p["xSec"]=0.1
        p["uncorrErr"].append(0.1)
        #
        DataAN.AddPoint(p)

print("Done.  ")
#%%

XX4=DataProcessor.harpyInterface.ComputeXSec(DataA4)
XX3=DataProcessor.harpyInterface.ComputeXSec(DataA3)
XXN=DataProcessor.harpyInterface.ComputeXSec(DataAN)

#%%
XX0=DataProcessor.harpyInterface.ComputeXSec(DataA0)
XX1=DataProcessor.harpyInterface.ComputeXSec(DataA1)
XX2=DataProcessor.harpyInterface.ComputeXSec(DataA2)

#%%
#rSet.SetReplica(0)
harpy.setNPparameters([1.5004, 0.05614, 0.03862, 0.0, 0.565, 0.0539, 0.5697, 6.64, 0.565, 20.07, 0.5697, 0.537, 1.07, 2.39, 0.0, 0.0])
#harpy.setNPparameters([1.4806, 0.038969, 0.051737, 0.0, 0.565, 0.0539, 0.5697, 6.64, 0.565, 20.07, 0.5697, 0.537, 1.07, 2.39, 0.0, 0.0])
DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)

#%%
pathToPlot="/data/WorkingFiles/TMD/Fit_Notes/KPC/PlotsData/KPC_ART23/"
import time
for s in setDY.sets:
    startT=time.time()
    XX=DataProcessor.harpyInterface.ComputeXSec(s)
    endT=time.time()
    print(':->',s.name,'       t=',endT-startT)
    f=open(pathToPlot+s.name+".dat","w")
    print('SAVING PLOTS>>  ',f.name)
    ### [total chi^2, list of NP-parameters],
    for i in range(len(XX)):
        p=s.points[i]
        f.write(str(p["qT"][0])+", "+str(p["qT"][1])+", "+
                str(p["xSec"])+", "+str(numpy.sqrt(sum(numpy.array(p["uncorrErr"])**2)))+", "+
                str(XX[i])+"\n")
        
    f.close()

#%%
#######################################
# Minimisation
#######################################
import time
totalN=setDY.numberOfPoints

### FOR ART23
#penalty_index=[-7,-6,-5,-4,-3]

penalty_index=[-1,-2,-3,-4]

def chi_2DY(x):
    startT=time.time()
    harpy.setNPparameters_TMDR([x[0],x[1],x[2],x[3]])
    #harpy.setNPparameters_uTMDPDF(x[4:])
    harpy.setNPparameters_uTMDPDF([x[4],x[5],x[6],x[7],x[4],x[9],x[6],x[11],x[12],x[13],x[14],x[15]])
    # harpy.setNPparameters_TMDR([x[0],x[1],x[2],x[3]])
    # harpy.setNPparameters_uTMDPDF(x[4:])
    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")        
    
    #YY=DataProcessor.harpyInterface.ComputeXSec(setDY)
    YY0=DataProcessor.harpyInterface.ComputeXSec(setDY,method="approximate")
    YY=RAT*numpy.array(YY0)
    ccDY2,cc3=setDY.chi2(YY)
    
    penalty_array=numpy.array([max(0,abs(setDY.sets[i].DetermineAvarageSystematicShift(YY[setDY._i1[i]:setDY._i2[i]]))/setDY.sets[i].normErr[0]-1) for i in penalty_index])
    penalty_term=sum(penalty_array**6)
    
    #ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)
    #ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(setSIDIS)
    
    cc=ccDY2/setDY.numberOfPoints
    endT=time.time()
    print(':->',cc,'   +p=',cc+penalty_term/setDY.numberOfPoints,'       t=',endT-startT)
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
                0.565, 20.07, 0.5697, 0.537, 
                1.07, 2.39, 0.0, 0.0])

# initialValues=([1.56142, 0.0369174, 0.0581734, 1.0,
#   0.874245, 0.913883, 0.991563, 6.05412,
#   0.353908, 46.6064, 0.115161, 1.53235,
#   1.31966, 0.434833, 0.0, 0.0])
initialErrors=(0.1,0.1,0.1,0.1,
                0.1,  1.0, 0.1,  1.0,
                0.1,  1.0, 0.1,  1.0,
                0.1,  1.0, 10.,  1.)
searchLimits=((1.0,2.5),(0.001,.2) ,(0.0,.2), (-5.,5.),
              (0.,10.), (0,100.),(0.,10.), (0,100.),
              (-0.1,10.), (0.,100.),(-0.1,10.), (0.,100.),
              (0.,10.), (0.,25.),(0,5000.), (0,100.))
# True= FIX
parametersToMinimize=(True, False,False,True,
                      False, False, False,False,
                      True, False, True,False,
                      False, False, True ,True)

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

