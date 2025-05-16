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
import time
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

path_to_constants=MAINPATH+"FittingPrograms/ART25_KPC/ConstantsFiles/DY_LP_newModel.atmde"


harpy.initialize(path_to_constants)

inARRAY_TMDR=[1.5, 0.071624, 0.064726, 0.0]
inARRAY_PDF=[0.1,0.5, 1., 
             0., 0., 0., 0., 0., 0., 0., 0.,0., 0., 0., 0.,0.]


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
    
    if p["type"]=="DY":
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
                          'E605'
                          #'D0run1-W','CDFrun1-W'
                          ]))

setDY=theData.CutData(cutFunc) 
#setDYfit=theData.CutData(cutFuncFORFIT) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')

#%%

####### Best fast result
# harpy.setNPparameters([1.5, 0.083931, 0.030641, 0.0, 
#                        0.51638, 0.002073, 0.478567, 0.373111, 
#                        2.407, 22.1996, 3.7876, 0.00128, 
#                        0.403343, 5e-05, 1.0, 1.0, 
#                        0.69769, 0.712969, -0.133895, -0.841651, 0.846846,
#                        0.774759, 1.5565, 1.1863, 0.692877, -0.569062, 
#                        0.0, 0.0])


harpy.setNPparameters([1.5, 0.063075, 0.021832, 0.0, 
0.000481, 0.690827, 1.200017,
-0.5807195034513588, -0.5630230221568071, 0.11700957345993641, 
0.4584247272053746, 0.4500188050150401, -0.05771358155743821, 
-0.2466533590440217, -1.898048818628432, -5.3103787197866, 
-0.2085032276220349, 0.014927742386172421, 0.0, 0.0])

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)

#%%
setHP=setDY
penalty_index=[-5,-4,-3,-2,-1]

def chi2_GLOBAL(theory):            
    ccDY2,cc3=setHP.chi2(theory)
    
    #### This penalty term prevents low-energy DY to have extremely low normalization
    penalty_array=numpy.array([max(0,abs(setHP.sets[i].DetermineAvarageSystematicShift(theory[setHP._i1[i]:setHP._i2[i]]))/setHP.sets[i].normErr[0]-1) for i in penalty_index])
    penalty_term=sum(penalty_array**6)
    
    
    cc=ccDY2/setDY.numberOfPoints
    
    #print(':->',cc,'   +p=',penalty_term/setDY.numberOfPoints)
    return ccDY2+penalty_term

#%%
###### index enumeration of slow parameters in the common NP-row
NPslow_index=[0,1,2,3,4,5,6]
##### initial value for slow parameters (according to the indices)
NPslow_initial=[1.5, 0.063075, 0.021832, 0.0, 
                0.000481, 0.690827, 1.000]

NPslow_error=(0.1,0.01 ,0.01, 1.,
              0.1, 0.1,0.1)

NPslow_limits=((1.0,2.5),(0.005,0.2) ,(0.0,.2), (-5.,5.),
              (0.00001,100.), (0.00001,100.),(0.2,10.))

#### True=FIX
NPslow_ToMinimize=(True, False,False,True,
                   False, False, True)

###### index enumeration of fast parameters in the common NP-row
NPfast_index=[7,8,9,10,11,12,13,14,15,16,17,18,19]

##### initial value for fast parameters (according to the indices)
NPfast_initial=[-0.581, -0.563, 0.117, 
                0.458, 0.450, -0.0577, 
                -0.247, -1.898, -5.31, 
                -0.209, 0.0149, 0.0, 0.0]

NPfast_error=(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,0.2, 0.2, 0.2, 0.2,0.2)

NPfast_limits=( (-50.,50.),(-50.,50.),(-50.,50.),(-50.,50),
               (-50.,50.),(-50.,50.),(-50.,50.),(-50.,50.),
               (-50.,50.),(-50.,50.),(-50.,50.),(-50.,50.),
               (-50.,50.))
NPfast_ToMinimize=(False,False,False,False,False,False,False,False,False,False,False,True,True)

totalNPsize=len(NPslow_index)+len(NPfast_index)

NPslow_current=NPslow_initial
NPfast_current=NPfast_initial

#%%
#### sets only slow parameters
def SetBothParameters(slowP,fastP):
    
    if(len(slowP)!=len(NPslow_index)):
        raise Exception("Number of slow parameter to set, is incorrect 2")
    if(len(fastP)!=len(NPfast_index)):
        raise Exception("Number of fast parameter to set, is incorrect 2")
    
    NProw=[0 for i in range(totalNPsize)]
    ### assigning slow parameters from global variable
    for i in range(len(NPslow_index)):
        indx=NPslow_index[i]
        NProw[indx]=float(slowP[i])
    
    for i in range(len(NPfast_index)):
        indx=NPfast_index[i]
        NProw[indx]=float(fastP[i])
        
    harpy.setNPparameters(NProw)

#### sets only slow parameters
def SetSlowParameters(x):
    
    if(len(x)!=len(NPslow_index)):
        raise Exception("Number of slow parameter to set, is incorrect")
    
    NProw=[0 for i in range(totalNPsize)]
    ### assigning slow parameters from global variable
    for i in range(len(NPslow_index)):
        indx=NPslow_index[i]
        NProw[indx]=float(x[i])
    
    for i in range(len(NPfast_index)):
        indx=NPfast_index[i]
        NProw[indx]=NPfast_current[i]
        
    harpy.setNPparameters(NProw)

#### sets only fast parameters
def SetFastParameters(x):
    
    if(len(x)!=len(NPfast_index)):
        raise Exception("Number of fast parameter to set, is incorrect")
    
    NProw=[0 for i in range(totalNPsize)]
    ### assigning slow parameters from global variable
    for i in range(len(NPslow_index)):
        indx=NPslow_index[i]
        NProw[indx]=NPslow_current[i]
    
    for i in range(len(NPfast_index)):
        indx=NPfast_index[i]
        NProw[indx]=float(x[i])
        
    harpy.setNPparameters(NProw)

##### subroutines which sets up only the fast parameters i,j to 1 or -1 (rest is zero)
def setETA1(i):    
    arr=[0. for i in range(len(NPfast_index))]
    arr[i]=1.
    SetFastParameters(arr)

def setETA1m(i):
    arr=[0. for i in range(len(NPfast_index))]
    arr[i]=-1.
    SetFastParameters(arr)
    
def setETA2(i,j):
    arr=[0. for i in range(len(NPfast_index))]
    arr[i]=1.
    arr[j]=1.
    SetFastParameters(arr)
    


#%%
NUM=len(NPfast_index)

numCalls_slow=0
numCalls_fast_total=0
numCalls_XSEC=0

R00=0.
beta=[0 for i in range(NUM)]
alpha=[[0 for x in range(NUM)] for y in range(NUM)] 
#%%
def PrepareMatixALPHA():
    global R00,beta,alpha,numCalls_XSEC
    
    nullArray=[0. for i in range(NUM)]
    ### create empty array
    beta=[0 for i in range(NUM)]
    alpha=[[0 for x in range(NUM)] for y in range(NUM)] 
    
    ### restoring the matrix of the theory prediction
    
    #### constant part
    SetFastParameters(nullArray)
    R00=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setHP))
    numCalls_XSEC=numCalls_XSEC+1
    
    xSecSize=len(R00) #### this is stored for empty parameters
    
    print("Constant part is made!")
    
    #### diagonal parts
    for i in range(NUM):
        
        if(NPfast_ToMinimize[i]):
            #### if parameter is fixed it is assumed zero forever!!
            alpha[i][i]=numpy.array([0. for i in range(xSecSize)])
            beta[i]=numpy.array([0. for i in range(xSecSize)])
        else:
            print("....",i)
            setETA1(i)
            f1=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setHP))
            setETA1m(i)
            f1m=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setHP))            
            alpha[i][i]=(f1+f1m)/2-R00
            beta[i]=f1-alpha[i][i]-R00
            numCalls_XSEC=numCalls_XSEC+2
    
    print("Diagonal part is made!")
    
    ##### off-diagonal part
    for i in range(NUM-1):
        for j in range(i+1,NUM):
            if(NPfast_ToMinimize[i] or NPfast_ToMinimize[j]):
                alpha[i][j]=numpy.array([0. for i in range(xSecSize)])
            else:
                print("....",i,", ",j)
                setETA2(i,j)
                f1=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setHP))
                alpha[i][j]=f1-R00-beta[i]-beta[j]-alpha[i][i]-alpha[j][j]
                numCalls_XSEC=numCalls_XSEC+1
    
    print("Off-diagonal part is made!")

#%%
import copy
def PredictionForFastPart(x):
    global numCalls_fast_total
    ### constant part
    t=copy.copy(R00)
    ### linear part
    for i in range(NUM):
        t=t+beta[i]*x[i]
        
    ### quadratic part
    for i in range(NUM):
        for j in range(i,NUM):
            t=t+alpha[i][j]*x[i]*x[j]   
            
    numCalls_fast_total=numCalls_fast_total+1
    return t


#%%
def OptimizeFastParameters(fastP):    
    global numCalls_fast_total
    
    c1_here=numCalls_fast_total
    
    PrepareMatixALPHA()
    
    from iminuit import Minuit
    
    
    startT=time.time()
    
    def chi2_local(x):
        theory=PredictionForFastPart(x)
        return chi2_GLOBAL(theory)
    
    m_local = Minuit(chi2_local, (fastP))

    m_local.errors=NPfast_error
    m_local.limits=NPfast_limits
    m_local.fixed=NPfast_ToMinimize
    m_local.errordef=1

    m_local.migrad()
    
    endT=time.time()
    
    print("Calls for fast optimization: ",numCalls_fast_total-c1_here, "  Optimization time: ",endT-startT)
    
    return chi2_local(list(m_local.values)), list(m_local.values)

#%%
##### computes the chi2 for change of slow parameters
##### uptemizes the fast parameters, and store them globally
##### initial point for optimization is current best
def chi2_slow(slowP):
    global NPfast_current
    global NPslow_current
    global numCalls_slow
    
    startT=time.time()
    
    NPslow_current=slowP
    
    SetSlowParameters(slowP)
    ch2,listFast=OptimizeFastParameters(NPfast_current)
    
    NPfast_current=listFast
    
    numCalls_slow=numCalls_slow+1
    
    endT=time.time()
    
    print(':->',ch2/setHP.numberOfPoints,"    time=",endT-startT)
    return ch2

#%%
    
from iminuit import Minuit


m_main = Minuit(chi2_slow, NPslow_initial)

m_main.errors=NPslow_error
m_main.limits=NPslow_limits
m_main.fixed=NPslow_ToMinimize
m_main.errordef=1

m_main.strategy=1
m_main.tol=0.0001*(setHP.numberOfPoints)*10000 ### the last 0.0001 is to compensate MINUIT def

m_main.migrad()

NPslow_current=list(m_main.values)

print("---------------- SLOW PARAMETERS ------------------------")

print(m_main.params)

print([round(x,1 if x >100 else 4 if x>1 else 6) for x in NPslow_current])

print("---------------- FAST PARAMETERS ------------------------")

print(NPfast_current)

#%%
SetBothParameters(NPslow_current,NPfast_current)

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)

print("Total calls of XSEC:      ", numCalls_XSEC)
print("Total calls of slow chi2: ", numCalls_slow)
print("Total calls of fast chi2: ", numCalls_fast_total)

#%%
sys.exit()
