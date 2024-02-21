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

MAINPATH=ROOT_DIR 
#%%
#######################################
#Initialize artemide
#######################################
import harpy

path_to_constants=MAINPATH+"FittingPrograms/ART23/ConstantsFiles/DYonly/ART23_MSHT_N4LL.atmde"


harpy.initialize(path_to_constants)

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                        "/data/WorkingFiles/TMD/Fit_Notes/ART23/REPS/ART23_run2.rep")

rSet.SetReplica(0)

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
    if(True):
        if(p["process"][2]==1): p["process"]=[p["process"][0],1,1,1]
        elif(p["process"][2]==2): p["process"]=[p["process"][0],1,1,-1]
        elif(p["process"][2]==3): p["process"]=[p["process"][0],2,1,1]
        elif(p["process"][2]==4): p["process"]=[p["process"][0],2,1,-1]
        elif(p["process"][2]==5): p["process"]=[p["process"][0],3,1,1]
        elif(p["process"][2]==6): p["process"]=[p["process"][0],3,1,-1]
        elif(p["process"][2]==7): p["process"]=[p["process"][0],4,1,1]
        elif(p["process"][2]==8): p["process"]=[p["process"][0],5,1,1]
        elif(p["process"][2]==9): p["process"]=[p["process"][0],6,1,1]
        elif(p["process"][2]==10): p["process"]=[p["process"][0],4,1,-1]
        elif(p["process"][2]==11): p["process"]=[p["process"][0],5,1,-1]
        elif(p["process"][2]==12): p["process"]=[p["process"][0],6,1,-1]
        elif(p["process"][2]==1001): p["process"]=[p["process"][0],101,1,1]
        elif(p["process"][2]==1002): p["process"]=[p["process"][0],102,1,1]
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
    
#    return delta<0.5 and p.qT_avarage<80
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
                          'E605',
                          'D0run1-W','CDFrun1-W'
                          ]))

setDY=theData.CutData(cutFunc) 
#setDYfit=theData.CutData(cutFuncFORFIT) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
#print('Loaded ', setDYfit.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDYfit.sets]), 'points.')

#%%
SAVEPATH="/data/WorkingFiles/TMD/Fit_Notes/ART23/"

expSets=[setDY.sets[3],setDY.sets[13],setDY.sets[14],setDY.sets[15],setDY.sets[34],setDY.sets[35]]

listN=[]
for s in expSets:
    XX=harpy.DY.xSecList([d["process"] for d in s.points],
                            [d["s"] for d in s.points],
                            [d["qT"] for d in s.points],
                            [d["Q"] for d in s.points],
                            [d["y"] for d in s.points],
                            [d["includeCuts"] for d in s.points],
                            [d["cutParams"] for d in s.points])
    
    ff=numpy.array([(p["qT"][1]-p["qT"][0])*p["thFactor"] for p in s.points])
    listN+=[sum(XX*ff)/s._normExp]
    
with open(SAVEPATH+'norms', 'w') as outfile:
    ll=[s.name for s in expSets]
    outfile.write(",  ".join(ll)+"\n")
    
    ll=list(map(lambda n: '%.4f'%n, listN))
    outfile.write(",  ".join(ll)+"\n")
    
#%%
for k in range(300):    
    rnd=numpy.random.randint(1,1000)
    print(str(k)+"/300"+"  rep="+str(rnd))
    rSet.SetReplica(rnd)    
    listN=[]
    
    for s in expSets:
        XX=harpy.DY.xSecList([d["process"] for d in s.points],
                                [d["s"] for d in s.points],
                                [d["qT"] for d in s.points],
                                [d["Q"] for d in s.points],
                                [d["y"] for d in s.points],
                                [d["includeCuts"] for d in s.points],
                                [d["cutParams"] for d in s.points])
        
        ff=numpy.array([(p["qT"][1]-p["qT"][0])*p["thFactor"] for p in s.points])
        listN+=[sum(XX*ff)/s._normExp]
        
    with open(SAVEPATH+'norms','a') as outfile:        
        ll=list(map(lambda n: '%.4f'%n, listN))
        outfile.write(",  ".join(ll)+"\n")