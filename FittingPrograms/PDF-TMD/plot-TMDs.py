#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 08:28:25 2020

@author: vla18041
"""

#######################################
# Paths
#######################################
PathToHarpy="/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/harpy"
PathToDataProcessor="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"

######## NNPDF
# PathToConstantsFile=PathToDataProcessor+"FittingPrograms/PDF-TMD/Constants-files/const-NNPDF31_NNLO"
# PathToReplicaFile="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/LOGS/NNPDF31/allPDF_centralEXP_full.dat"
# PathToSavePlot="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/NNPDF31/PlotData"
# numberOfLines=921
# CentralReplica=[2.000, 0.030, 0.254, 8.550, 373.112,  2.529, -5.910, 0.000, 0.000]
######## HERA
PathToConstantsFile=PathToDataProcessor+"FittingPrograms/PDF-TMD/Constants-files/const-HERA20_NNLO"
PathToReplicaFile="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/LOGS/HERA20/allPDF_centralEXP_full.dat"
PathToSavePlot="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/HERA20/PlotData"
numberOfLines=981
CentralReplica=[2., 0.0331, 0.305, 11.38, 423.8, 2.024, -8.612, 0., 0.]

import sys
#sys.path.append(PathToHarpy)
sys.path.append(PathToDataProcessor)


#%%
#######################################
# importing libraries
#######################################
import numpy
import time

import DataProcessor.harpyInterface


#%%
#######################################
#Initialize artemide
#######################################
import harpy

harpy.initialize(PathToConstantsFile)
initializationArray=[2.0340, 0.0299, 0.2512, 7.7572,334.6108, 2.4543,-4.8203, 0.1000,  0.0000]
harpy.setNPparameters(initializationArray)

#%%
############################################
## function which reads the line N from ReplicaFile
############################################
def SetReplicaEntry(n):
    import ast
    fp = open(PathToReplicaFile,"r")
    for i, line in enumerate(fp):
        if i == n-1:
            entry=ast.literal_eval(line)
            harpy.setPDFreplica(entry[3])
            harpy.setNPparameters(entry[4])  
        elif i > n:
            break
    fp.close()
    
#%%
############################################
## Save a plot of particular case N
############################################
xValues=[0.001, 0.00178, 0.00316, 0.00562, 0.01, 0.0178, 0.0316, 0.0562, 0.1, 0.126, 0.158, 0.120, 0.251, 0.316, 0.398,0.501,
         0.562, 0.631, 0.708, 0.794, 0.891,0.95]
bValues=[0.01,0.1,0.2,0.3,0.4,0.5,0.75,1.,1.5,2.0,3.0,5.0]

#fileName="NNPDF_full_central"
fileName="HERA_full_central"

def SavePlot(n):
    SetReplicaEntry(n)
    
    plot=[]
    
    for x in xValues:
        for b in bValues:
            pp=harpy.get_uTMDPDF(x,b,1)
            plot.append([x,b,list(pp)])
    
    print(PathToSavePlot+"/"+fileName+"_"+str(n))
    fp = open(PathToSavePlot+"/"+fileName+"_"+str(n)+".dat","w")
    for line in plot:
        fp.write(str(line)+"\n")
    fp.close()
    
#%%
for i in range(1,numberOfLines+1):
    SavePlot(i)
    
#%%
############################################
## Save a plot of particular case N (for central)
############################################
xValues=[0.001, 0.00178, 0.00316, 0.00562, 0.01, 0.0178, 0.0316, 0.0562, 0.1, 0.126, 0.158, 0.120, 0.251, 0.316, 0.398,0.501,
         0.562, 0.631, 0.708, 0.794, 0.891,0.95]
bValues=[0.01,0.1,0.2,0.3,0.4,0.5,0.75,1.,1.5,2.0,3.0,5.0]

#fileName="NNPDF_noFIT"
fileName="HERA_noFIT"

 
def SavePlot0(n):
    
    harpy.setPDFreplica(n)
    harpy.setNPparameters(CentralReplica) 
    
    plot=[]
    
    for x in xValues:
        for b in bValues:
            pp=harpy.get_uTMDPDF(x,b,1)
            plot.append([x,b,list(pp)])
    
    print(PathToSavePlot+"/"+fileName+"_"+str(n))
    fp = open(PathToSavePlot+"/"+fileName+"_"+str(n)+".dat","w")
    for line in plot:
        fp.write(str(line)+"\n")
    fp.close()
    
#%%
numberOfLines=1000

for i in range(numberOfLines+1):
    SavePlot0(i)