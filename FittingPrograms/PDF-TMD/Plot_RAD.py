#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 14:04:58 2021

@author: vla18041
"""

#%%
#######################################################################
# Global parameter of a run
#######################################################################

PDFinUse="HERA20"
#PDFinUse="NNPDF31"
#PDFinUse="CT18"
#PDFinUse="MSHT20"
#PDFinUse="CJ15"

#CASE="EXP"
CASE="PDF"

#######################################
# Paths
#######################################
PathToHarpy="/home/vla18041/LinkData2/arTeMiDe_Repository/artemide-ForPDF/harpy"
#PathToHarpy="/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/harpy"
PathToDataProcessor="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"
PathToDataLibrary=PathToDataProcessor+"DataLib/unpolDY/"
PathToConstantsFile=PathToDataProcessor+"/FittingPrograms/PDF-TMD/Constants-files/const-"+PDFinUse+"_NNLO_12p"


import sys
sys.path.remove('/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/harpy')
sys.path.append(PathToHarpy)
sys.path.append(PathToDataProcessor)

#%%
#######################################
# importing libraries
#######################################
import numpy

import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet

#%%
#######################################
#Initialize artemide
#######################################
import harpy

harpy.initialize(PathToConstantsFile)
#%%
rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/REPS/SV21-"+PDFinUse+"-nnlo-"+CASE+".rep")
rSet.SetReplica(0,part="TMDR")

#%%

#### compute the 68CI for list of arrays, axis=0
def Bootstrap68CI(dd):
    rMin=[]
    rMax=[]
    
    for i in range(2000):
        indices=numpy.random.choice(range(len(dd)),size=int(len(dd)/2))
        sample=[dd[i] for i in indices]
        rMin.append(numpy.quantile(sample, (1-.68)/2,axis=0))
        rMax.append(numpy.quantile(sample, 1-(1-.68)/2,axis=0))
        
    return (numpy.array([numpy.mean(rMin,axis=0),numpy.mean(rMax,axis=0)])).transpose()

#### compute the 68CI for sum of 2 distributions list of arrays, axis=0
def BootstrapSUM68CI(dd1,dd2):
    rMin=[]
    rMax=[]
    
    for i in range(2000):        
        indices1=numpy.random.choice(range(len(dd1)),size=numpy.min([len(dd1),len(dd2)]))
        indices2=numpy.random.choice(range(len(dd2)),size=numpy.min([len(dd1),len(dd2)]))
        sample=[dd1[i] for i in indices1]+[dd2[i] for i in indices2]
        rMin.append(numpy.quantile(sample, (1-.68)/2,axis=0))
        rMax.append(numpy.quantile(sample, 1-(1-.68)/2,axis=0))
        
    return (numpy.array([numpy.mean(rMin,axis=0),numpy.mean(rMax,axis=0)])).transpose()

#%%

#PDFinUse="HERA20"
#PDFinUse="NNPDF31"
#PDFinUse="CT18"
PDFinUse="MSHT20"

bValues=[0.03,0.055, 0.105, 0.155, 0.205, 0.255, 0.305, 0.355, 0.405, 0.455, \
0.505, 0.555, 0.605, 0.655, 0.705, 0.755, 0.805, 0.855, 0.905, 0.955, \
1.005, 1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, \
3.8, 4., 4.2, 4.4, 4.6, 4.8, 5.]

repsEXP=[]
rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/REPS/SV21-"+PDFinUse+"-nnlo-EXP.rep")
for pp in range(rSet.numberOfReplicas):   
    rSet.SetReplica(pp,part="TMDR")
    repsEXP.append([harpy.get_DNP(b, 2.) for b in bValues])
    
repsPDF=[]
rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/REPS/SV21-"+PDFinUse+"-nnlo-PDF.rep")
for pp in range(rSet.numberOfReplicas):   
    rSet.SetReplica(pp,part="TMDR")
    repsPDF.append([harpy.get_DNP(b, 2.) for b in bValues])
    
## just mean
EXPmean=numpy.mean(repsEXP,axis=0)
PDFmean=numpy.mean(repsPDF,axis=0)

## 68 CI's
EXP68=Bootstrap68CI(repsEXP)
PDF68=Bootstrap68CI(repsPDF)

## Computing weighted mean
EXPsigma2=numpy.array([((c[1]-c[0])/2)**(-2) for c in EXP68])
PDFsigma2=numpy.array([((c[1]-c[0])/2)**(-2) for c in PDF68])

EXPw=EXPsigma2/(EXPsigma2+PDFsigma2)
PDFw=PDFsigma2/(EXPsigma2+PDFsigma2)

wMean=EXPw*EXPmean+PDFw*PDFmean

## 68 CI for the sum
SUM68=BootstrapSUM68CI(repsEXP,repsPDF)

for i in range(len(bValues)):
    print("{"+str(bValues[i])+", "+'{:2.6f}'.format(wMean[i])+", "+'{:2.6f}'.format(SUM68[i][0])+", "+'{:2.6f}'.format(SUM68[i][1])+"},")