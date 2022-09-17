#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 14:35:35 2021

@author: vla18041
"""

#%%
#######################################################################
# Global parameter of a run
#######################################################################

#PDFinUse="HERA20"
#PDFinUse="NNPDF31"
#PDFinUse="CT18"
PDFinUse="MSHT20"
#PDFinUse="CJ15"

#CASE="EXP"
CASE="PDF"

nn=9
if(CASE=="EXP"): nn=0

#######################################
# Paths
#######################################
import os
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"

ATMDE_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide-ForPDF/"

#PathToHarpy="/home/vla18041/LinkData2/arTeMiDe_Repository/artemide-ForPDF/harpy"
#PathToHarpy="/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/harpy"
PathToHarpy=ATMDE_DIR+"harpy"
PathToDataProcessor=ROOT_DIR
PathToDataLibrary=PathToDataProcessor+"DataLib/unpolDY/"
PathToConstantsFile=PathToDataProcessor+"/FittingPrograms/PDF-TMD/Constants-files/const-"+PDFinUse+"_NNLO_12p_EIG"

import sys
#sys.path.remove('/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/harpy')
sys.path.remove('/data/arTeMiDe_Repository/artemide/harpy')
sys.path.append(PathToHarpy)
sys.path.append(PathToDataProcessor)


#%%
#######################################
# importing libraries
#######################################
import numpy

import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet
import DataProcessor.ArtemideReplicaSet

#%%
#######################################
#Initialize artemide
#######################################
import harpy

harpy.initialize(PathToConstantsFile)
#%%
#PathToWorking='/home/vla18041/LinkData2/WorkingFiles/'
PathToWorking='/data/WorkingFiles/'
rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            PathToWorking+"TMD/Fit_Notes/TMD<->PDF/REPS/SV21-"+PDFinUse+"-nnlo-"+CASE+".rep")
rSet.SetReplica()

#%%
#######################################
# read the list of files and return the list of DataSets
#######################################
def loadThisData(listOfNames):    
    import DataProcessor.DataSet
    
    dataCollection=[]
    for name in listOfNames:
        if( name==''): continue
        loadedData=DataProcessor.DataSet.LoadCSV(PathToDataLibrary+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

#%%
##################Cut functions
includePoints=[
    'CDF2.0', 'CDF2.1', 'CDF2.2', 'CDF2.3', 'CDF2.4', 'CDF2.5', 'CDF2.6', 'CDF2.7',
    'CDF2.8', 'CDF2.9', 'CDF2.10', 'CDF2.11', 'CDF2.12', 'CDF2.13', 'CDF2.14',
    'D02m.0', 'D02m.1', 'D02m.2',
    'A8-00y04.0', 'A8-00y04.1', 'A8-00y04.2', 'A8-00y04.3', 'A8-00y04.4', 'A8-04y08.0',
    'A8-04y08.1', 'A8-04y08.2', 'A8-04y08.3', 'A8-04y08.4', 'A8-08y12.0', 'A8-08y12.1',
    'A8-08y12.2', 'A8-08y12.3', 'A8-08y12.4', 'A8-12y16.0', 'A8-12y16.1', 'A8-12y16.2',
    'A8-12y16.3', 'A8-12y16.4', 'A8-16y20.0', 'A8-16y20.1', 'A8-16y20.2', 'A8-16y20.3',
    'A8-16y20.4', 'A8-20y24.0', 'A8-20y24.1', 'A8-20y24.2', 'A8-20y24.3', 'A8-20y24.4',
    'A8-46Q66.0', 'A8-46Q66.1', 'A8-46Q66.2',
    'LHCb7.0', 'LHCb7.1', 'LHCb7.5', 'LHCb7.6', 'LHCb8.0', 'LHCb8.1',
    'LHCb8.2', 'LHCb8.3', 'LHCb8.4', 'LHCb8.5', 'LHCb8.6', 'PHE200.0', 'PHE200.1',
    'E228-200.4Q5.0', 'E228-200.4Q5.1', 'E228-200.4Q5.2', 'E228-200.4Q5.3', 'E228-200.4Q5.4',
    'E228-200.4Q5.5', 'E228-200.5Q6.0', 'E228-200.5Q6.1', 'E228-200.5Q6.2', 'E228-200.5Q6.3',
    'E228-200.5Q6.4', 'E228-200.5Q6.5', 'E228-200.5Q6.6', 'E228-200.6Q7.0', 'E228-200.6Q7.1',
    'E228-200.6Q7.2', 'E228-200.6Q7.3', 'E228-200.6Q7.4', 'E228-200.6Q7.5', 'E228-200.6Q7.6',
    'E228-200.6Q7.7', 'E228-200.6Q7.8', 'E228-200.7Q8.0', 'E228-200.7Q8.1', 'E228-200.7Q8.2',
    'E228-200.7Q8.3', 'E228-200.7Q8.4', 'E228-200.7Q8.5', 'E228-200.7Q8.6', 'E228-200.7Q8.7',
    'E228-200.7Q8.8', 'E228-200.7Q8.9', 'E228-200.8Q9.0', 'E228-200.8Q9.1', 'E228-200.8Q9.2',
    'E228-200.8Q9.3', 'E228-200.8Q9.4', 'E228-200.8Q9.7', 'E228-200.8Q9.8', 'E228-300.4Q5.0<u',
    'E228-300.4Q5.1<u', 'E228-300.4Q5.2<u', 'E228-300.4Q5.3<u', 'E228-300.4Q5.4<u', 'E228-300.4Q5.5<u',
    'E228-300.5Q6.0<u', 'E228-300.5Q6.1<u', 'E228-300.5Q6.2<u', 'E228-300.5Q6.3<u', 'E228-300.5Q6.4<u',
    'E228-300.5Q6.5<u', 'E228-300.5Q6.6<u', 'E228-300.6Q7.0<u', 'E228-300.6Q7.1<u', 'E228-300.6Q7.2<u',
    'E228-300.6Q7.3<u', 'E228-300.6Q7.4<u', 'E228-300.6Q7.5<u', 'E228-300.6Q7.6<u', 'E228-300.6Q7.7<u',
    'E228-300.6Q7.8<u', 'E228-300.7Q8.0<u', 'E228-300.7Q8.1<u', 'E228-300.7Q8.2<u', 'E228-300.7Q8.3<u',
    'E228-300.7Q8.4<u', 'E228-300.7Q8.5<u', 'E228-300.7Q8.6<u', 'E228-300.7Q8.7<u', 'E228-300.7Q8.8<u',
    'E228-300.7Q8.9<u', 'E228-300.8Q9.0<u', 'E228-300.8Q9.1<u', 'E228-300.8Q9.2<u', 'E228-300.8Q9.3<u',
    'E228-300.8Q9.4<u', 'E228-300.8Q9.5<u', 'E228-300.8Q9.6<u', 'E228-300.8Q9.7<u', 'E228-300.8Q9.8<u',
    'E228-300.8Q9.9<u', 'E228-300.8Q9.10<u', 'E228-300.11Q12.0>u', 'E228-300.11Q12.1>u',
    'E228-300.11Q12.2>u', 'E228-300.11Q12.3>u', 'E228-300.11Q12.4>u', 'E228-300.11Q12.5>u',
    'E228-300.11Q12.6>u', 'E228-300.11Q12.7>u', 'E228-300.11Q12.8>u', 'E228-300.11Q12.11>u',
    'E228-400.5Q6.0<u', 'E228-400.5Q6.1<u', 'E228-400.5Q6.2<u', 'E228-400.5Q6.3<u', 'E228-400.5Q6.4<u',
    'E228-400.5Q6.5<u', 'E228-400.5Q6.6<u', 'E228-400.6Q7.0<u', 'E228-400.6Q7.1<u', 'E228-400.6Q7.2<u',
    'E228-400.6Q7.3<u', 'E228-400.6Q7.4<u', 'E228-400.6Q7.5<u', 'E228-400.7Q8.0<u', 'E228-400.7Q8.1<u',
    'E228-400.7Q8.2<u', 'E228-400.7Q8.3<u', 'E228-400.7Q8.4<u', 'E228-400.7Q8.5<u', 'E228-400.7Q8.6<u',
    'E228-400.7Q8.7<u', 'E228-400.7Q8.8<u', 'E228-400.7Q8.9<u', 'E228-400.8Q9.0<u', 'E228-400.8Q9.1<u',
    'E228-400.8Q9.2<u', 'E228-400.8Q9.3<u', 'E228-400.8Q9.4<u', 'E228-400.8Q9.5<u', 'E228-400.8Q9.6<u',
    'E228-400.8Q9.7<u', 'E228-400.8Q9.8<u', 'E228-400.8Q9.9<u', 'E228-400.8Q9.10<u', 'E228-400.11Q12.0>u',
    'E228-400.11Q12.1>u', 'E228-400.11Q12.2>u', 'E228-400.11Q12.3>u', 'E228-400.11Q12.4>u',
    'E228-400.11Q12.5>u', 'E228-400.11Q12.6>u', 'E228-400.11Q12.7>u', 'E228-400.11Q12.8>u',
    'E228-400.11Q12.9>u', 'E228-400.11Q12.10>u', 'E228-400.11Q12.11>u', 'E228-400.11Q12.12>u',
    'E228-400.11Q12.13>u', 'E228-400.11Q12.14>u', 'E228-400.12Q13.0>u', 'E228-400.12Q13.1>u',
    'E228-400.12Q13.2>u', 'E228-400.12Q13.3>u', 'E228-400.12Q13.4>u', 'E228-400.12Q13.5>u',
    'E228-400.12Q13.6>u', 'E228-400.12Q13.7>u', 'E228-400.12Q13.8>u', 'E228-400.12Q13.9>u',
    'E228-400.12Q13.10>u', 'E228-400.12Q13.11>u', 'E228-400.12Q13.13>u', 'E228-400.12Q13.14>u',
    'E228-400.12Q13.15>u', 'E228-400.13Q14.0>u', 'E228-400.13Q14.1>u', 'E228-400.13Q14.2>u',
    'E228-400.13Q14.3>u', 'E228-400.13Q14.4>u', 'E228-400.13Q14.5>u', 'E228-400.13Q14.6>u',
    'E228-400.13Q14.7>u', 'E228-400.13Q14.8>u', 'E228-400.13Q14.9>u', 'E228-400.13Q14.11>u',
    'E228-400.13Q14.12>u', 'E772.11Q12.0', 'E772.11Q12.1', 'E772.11Q12.2', 'E772.11Q12.3',
    'E772.11Q12.4', 'E772.11Q12.5', 'E772.11Q12.6', 'E772.11Q12.7', 'E772.11Q12.8', 'E772.11Q12.9',
    'E772.11Q12.10', 'E772.12Q13.0', 'E772.12Q13.1', 'E772.12Q13.2', 'E772.12Q13.3', 'E772.12Q13.4',
    'E772.12Q13.5', 'E772.12Q13.7', 'E772.13Q14.0', 'E772.13Q14.1', 'E772.13Q14.2', 'E772.13Q14.3',
    'E772.13Q14.6', 'E772.14Q15.1', 'E605.7Q8.0<u', 'E605.7Q8.1<u', 'E605.7Q8.2<u', 'E605.7Q8.3<u',
    'E605.7Q8.4<u', 'E605.7Q8.5<u', 'E605.7Q8.6<u', 'E605.7Q8.7<u', 'E605.7Q8.8<u', 'E605.7Q8.9<u',
    'E605.8Q9.0<u', 'E605.8Q9.1<u', 'E605.8Q9.2<u', 'E605.8Q9.3<u', 'E605.8Q9.4<u', 'E605.8Q9.5<u',
    'E605.8Q9.6<u', 'E605.8Q9.7<u', 'E605.8Q9.8<u', 'E605.8Q9.9<u', 'E605.8Q9.10<u', 'E605.11Q13.0>u',
    'E605.11Q13.1>u', 'E605.11Q13.2>u', 'E605.11Q13.3>u', 'E605.11Q13.4>u', 'E605.11Q13.5>u',
    'E605.11Q13.6>u', 'E605.11Q13.7>u', 'E605.11Q13.8>u', 'E605.11Q13.9>u', 'E605.11Q13.10>u',
    'E605.11Q13.11>u', 'E605.11Q13.12>u', 'E605.11Q13.13>u', 'E605.11Q13.14>u', 'E605.11Q13.15>u',
    'E605.11Q13.16>u', 'E605.13Q18.0>u', 'E605.13Q18.1>u', 'E605.13Q18.2>u', 'E605.13Q18.3>u',
    'E605.13Q18.4>u', 'E605.13Q18.5>u', 'E605.13Q18.6>u', 'E605.13Q18.7>u', 'E605.13Q18.8>u',
    'E605.13Q18.9>u', 'E605.13Q18.10>u', 'E605.13Q18.11>u', 'E605.13Q18.12>u', 'E605.13Q18.13>u', 'E605.13Q18.14>u']


def cutFuncSV19(p):
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
            
        if(p["id"][-2:]=="<u" and p["<Q>"]>10.5):
            return False,p
        if(p["id"][-2:]==">u" and p["<Q>"]<10.5):
            return False,p
    
#    return delta<0.5 and p.qT_avarage<80
    return (delta<0.1 or (delta<0.25 and par/err*delta**2<1)) , p

### check the point against the list of dropping points.
def cutFuncSV21(p):    
    
    #### check against the presence in the reduced set.
    if(not(p["id"] in includePoints)):
        return False,p
    
    return cutFuncSV19(p)

def cutFuncPlot(p):
    if p["type"]=="DY":        
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
        
        if(p["id"][-2:]=="<u" and p["<Q>"]>10.5):
            return False,p
        if(p["id"][-2:]==">u" and p["<Q>"]<10.5):
            return False,p
    
#    return delta<0.5 and p.qT_avarage<80
    return (delta<0.35 or p["<qT>"]<10.) , p

#%%
#######################################
# Loading the data set
#######################################

setHE=loadThisData(['CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                      'A7-00y10','A7-10y20','A7-20y24', 
                      'A8-00y04','A8-04y08','A8-08y12','A8-12y16',
                      'A8-16y20','A8-20y24','A8-46Q66','A8-116Q150',
                      'CMS7', 'CMS8', 
                      'CMS13-00y04','CMS13-04y08','CMS13-08y12','CMS13-12y16','CMS13-16y24',
                      'LHCb7', 'LHCb8', 'LHCb13'])

#### I separate data above and below UPSILON, I create two copies of LE data with different names
#### the data to be split only  'E228-300', 'E228-400' and E605
setLE1=loadThisData(['PHE200', 'E228-200','E772'])
setLE2=loadThisData(['E228-300', 'E228-400','E605'])
setLE3=loadThisData(['E228-300', 'E228-400','E605'])
for s in setLE2:
    s.name+="-blwUPS"
    for p in s.points:
        p["id"]+="<u"
for s in setLE3:
    s.name+="-abvUPS"
    for p in s.points:
        p["id"]+=">u"
setLE=[setLE1[0],setLE1[1],setLE2[0],setLE3[0],setLE2[1],setLE3[1],setLE1[2],setLE2[2],setLE3[2]]


theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",setHE+setLE)

setSV21=theData.CutData(cutFuncSV21) 

setSV19=theData.CutData(cutFuncSV19) 

setPLOT=theData.CutData(cutFuncPlot)

print('SV21: Loaded ', setSV21.numberOfSets, 'data sets with ', sum([i.numberOfPoints for i in setSV21.sets]), 'points.')

print('SV19: Loaded ', setSV19.numberOfSets, 'data sets with ', sum([i.numberOfPoints for i in setSV19.sets]), 'points.')

print('PLOT: Loaded ', setPLOT.numberOfSets, 'data sets with ', sum([i.numberOfPoints for i in setPLOT.sets]), 'points.')

#%%
if(PDFinUse == "MSHT20"):
    harpy.setNPparameters([2.,0.0435856, 0.120811, 0.316864, 0.366887, 1.6959, 0.00369852, 56.4383, 0.00123501, 1.09761, 0.107047, 5.07336, 28.6827, 0])

if(PDFinUse == "HERA20"):
    harpy.setNPparameters([2.,0.0326349, 0.105422, 8.15213, 0.442909, 0.109655, 0.11605, 6.52881, 0.34749, 0.0498677, 0.491896, 5.22109, 336.538, 0])

if(PDFinUse == "NNPDF31"):
    harpy.setNPparameters([2.,0.026075, 0.284306, 2.57804, 0.398623, 1.10417, 0.0874515, 12.5858, 0.0981801, 6.05093, 0.251195, 3.00048, 180.76, 0])

if(PDFinUse == "CT18"):
    harpy.setNPparameters([2.,0.0439572, 0.0518259, 0.896614, 0.29299, 4.71699, 0.00940138, 55.5163, 0.0012351, 0.371885, 0.0122727, 9.0651, 24.2762, 0])

if(not (PDFinUse == "MSHT20" or PDFinUse == "HERA20" or PDFinUse == "NNPDF31" or PDFinUse == "CT18")):
    print("________________________________________________")
    print("______UNKNOWN SET_______________________________")
    print("________________________________________________")

#rSet.SetReplica(0)

# DataProcessor.harpyInterface.PrintChi2Table(setSV21,printDecomposedChi2=True)
DataProcessor.harpyInterface.PrintChi2Table(setSV19,printDecomposedChi2=True)

#%%

XX=[2.,0.0435856, 0.120811, 0.316864, 0.366887, 1.6959, 0.00369852, 56.4383, 0.00123501, 1.09761, 0.107047, 5.07336, 28.6827, 0]

sigma0=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setSV19))
EIG=[]

for i in range(64):
    print("----------------",i+1,"----------------------------")
    harpy.setPDFreplica(i+1)
    harpy.setNPparameters(XX)
    dS=numpy.array(DataProcessor.harpyInterface.ComputeXSec(setSV19))
    EIG.append(dS-sigma0)

#%%
EIG=numpy.array(EIG)
err0=numpy.sum(EIG*EIG,axis=0)
#%%
import copy
setWITHERR=copy.deepcopy(setSV19)

uncE=0.
corE=numpy.sqrt(1-uncE**2)

for i in range(setWITHERR.numberOfSets):
    x0=err0[setWITHERR._i1[i]:setWITHERR._i2[i]]
    ss=setWITHERR.sets[i]
    ss.numOfUncorrErr+=1
    ss.numOfCorrErr+=1
    for j in range(ss.numberOfPoints):
        p=ss.points[j]
        p["uncorrErr"].append(uncE*x0[j])
        p["corrErr"].append(corE*x0[j])
    ss.FinalizeSet()

DataProcessor.harpyInterface.PrintChi2Table(setWITHERR,printDecomposedChi2=True)

print(DataProcessor.harpyInterface.ComputeChi2(setWITHERR))
