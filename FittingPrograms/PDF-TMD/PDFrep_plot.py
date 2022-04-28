#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: vla18041
"""

#%%
#######################################################################
# Global parameter of a run
#######################################################################

constName="const-CT18_NNLO_7p"
mainName="allPDF_centralEXP_full"
replicaFileName=mainName+".dat"
pathTorepicas="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/LOGS/CT18NNLO/"

#%%
#######################################
# Paths
#######################################
PathToHarpy="/home/vla18041/LinkData2/arTeMiDe_Repository/artemide-ForPDF/harpy"
PathToDataProcessor="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"
PathToDataLibrary=PathToDataProcessor+"DataLib/unpolDY/"
PathToLog=PathToDataProcessor+"FittingPrograms/PDF-TMD/LOGS/"
PathToConstantsFile=PathToDataProcessor+"/FittingPrograms/PDF-TMD/Constants-files/"+constName

import sys
#sys.path.remove('/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/harpy')
#sys.path.append(PathToHarpy)
sys.path.append(PathToDataProcessor)

#%%
#######################################
# importing libraries and Initialize artemide
#######################################
import numpy

import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet

import harpy
harpy.initialize(PathToConstantsFile)

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
#######################################
# Data cut function
#######################################
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
    return (delta<0.1 or (delta<0.25 and par/err*delta**2<1)) , p

#%%
#######################################
# Loading the data set
#######################################
theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisData(['CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                      'A7-00y10','A7-10y20','A7-20y24',
                      'A8-00y04','A8-04y08','A8-08y12',
                      'A8-12y16','A8-16y20','A8-20y24',
                      'A8-46Q66','A8-116Q150',
                      'CMS7', 'CMS8', 
                      'LHCb7', 'LHCb8', 'LHCb13', 
                      'PHE200', 'E228-200', 'E228-300', 'E228-400', 
                      'E772',
                      'E605']))

setDY=theData.CutData(cutFunc) 

print('Loaded ', setDY.numberOfSets, 'data sets with ', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setDY.sets])

#%%
#######################################
# Read input information
#######################################
f=open(pathTorepicas+replicaFileName,"r")
PDFrep=f.readlines()
f.close()

import ast
PDFrep=[ast.literal_eval(x) for x in PDFrep]

def setReplica(n):
    line=PDFrep[n]
    harpy.setPDFreplica(line[3])
    harpy.setNPparameters(line[4])

#%%
#######################################
# Building plot for each replica
#######################################
longTable=[]
for i in range(len(PDFrep)):
    print("-------------------------------------- ", i, " / ",len(PDFrep),"---------------------------------------" )
    setReplica(i)
    dd=DataProcessor.harpyInterface.ComputeXSec(setDY)
    
    longTable.append(dd)

#%%
#######################################
# Saving result
# [set name, Xsec, XsecErr, [qT], [listOfPredictions]]
#######################################
f=open(pathTorepicas+mainName+".plt","w")
i=0
for currentSet in setDY.sets:
    for k in range(currentSet.numberOfPoints):
        currentPoint=currentSet.points[k]
        dd=[l[i] for l in longTable]
        lineToFile=[currentSet.name,
                    currentPoint["xSec"],
                    numpy.sqrt(numpy.sum(numpy.array(currentPoint["uncorrErr"])**2)),
                    currentPoint["qT"],
                    dd]
        f.write(str(lineToFile)+"\n")
        i+=1
f.close()