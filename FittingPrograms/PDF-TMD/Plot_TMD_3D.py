#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: vla18041
"""

#%%
#######################################################################
# Global parameter of a run
#######################################################################

#PDFinUse="HERA20"
#PDFinUse="NNPDF31"
PDFinUse="CT18"

runName="model2.2_"+PDFinUse+"_PDFrep_noA7_spltUPS"
runNameSave="model2.2_"+PDFinUse+"_heraNP"

print(" RUN: "+runName)

#%%
#######################################
# Paths
#######################################
PathToHarpy="/home/vla18041/LinkData2/arTeMiDe_Repository/artemide-ForPDF/harpy"
PathToDataProcessor="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"
PathToDataLibrary=PathToDataProcessor+"DataLib/unpolDY/"
PathToLog=PathToDataProcessor+"FittingPrograms/PDF-TMD/LOGS/"
PathToConstantsFile=PathToDataProcessor+"/FittingPrograms/PDF-TMD/Constants-files/const-"+PDFinUse+"_NNLO_7p"

import sys
sys.path.remove('/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/harpy')
sys.path.append(PathToHarpy)
sys.path.append(PathToDataProcessor)

#%%
#######################################
# Output paths
#######################################
import socket
PCname=socket.gethostname()

replicaFile=PathToLog+runName+".txt"
logFile=PathToLog+PCname+".log"

#%%
#######################################
# importing libraries
#######################################
import numpy

#%%
#######################################
#Initialize artemide
#######################################
import harpy

harpy.initialize(PathToConstantsFile)

#%%
#######################################
# loading the replicas information
#######################################
import ast

pathTodat="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/LOGS/"
f=open(pathTodat+"model2.2/"+runName+".dat","r")
replicas=f.readlines()
f.close()
replicas=[ast.literal_eval(x) for x in replicas]

#%%
#######################################
# definition of grid
#######################################
xValues=[
0.0001, 0.00011, 0.00013, 0.00014, 0.00016, 0.00018, 0.0002, 0.00022, \
0.00025, 0.00028, 0.00032, 0.00035, 0.0004, 0.00045, 0.0005, 0.00056, 0.00063, \
0.00071, 0.00079, 0.00089,\
0.001, 0.0011, 0.0013, 0.0014, 0.0016, 0.0018, 0.002, 0.0022, \
0.0025, 0.0028, 0.0032, 0.0035, 0.004, 0.0045, 0.005, 0.0056, 0.0063, \
0.0071, 0.0079, 0.0089, 0.01, 0.011, 0.013, 0.014, 0.016, 0.018, \
0.02, 0.022, 0.025, 0.028, 0.032, 0.035, 0.04, 0.045, 0.05, 0.056, \
0.063, 0.071, 0.079, 0.089, 0.1, 0.11, 0.13, 0.14, 0.16, 0.18, 0.2, \
0.22, 0.25, 0.28, 0.32, 0.35, 0.4, 0.45, 0.5, 0.56, 0.63, 0.71, 0.79, \
0.89, 0.99]
bValues=[0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, \
0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.1, 1.2, 1.3, 1.4, 1.5, \
1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, \
3., 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4., 4.5, 5., 5.5, \
6., 6.5, 7., 7.5, 8.]

    
#%%
#######################################
# compute avarage and std
# and save to file
#######################################
def SaveResults():
    print('-------------------------SAVE to FILE-----------------------------')
    meanU=numpy.mean(u,axis=0)
    meanD=numpy.mean(d,axis=0)
    meanS=numpy.mean(s,axis=0)
    meanUbar=numpy.mean(ubar,axis=0)
    meanDbar=numpy.mean(dbar,axis=0)
    meanSbar=numpy.mean(sbar,axis=0)
    
    stdU=numpy.std(u,axis=0)
    stdD=numpy.std(d,axis=0)
    stdS=numpy.std(s,axis=0)
    stdUbar=numpy.std(ubar,axis=0)
    stdDbar=numpy.std(dbar,axis=0)
    stdSbar=numpy.std(sbar,axis=0)

    if(mu<0):
        prefix="-b-optimal"
    else:
        prefix="-b-"+str(int(mu))+"GeV"
    
    fileu=open(pathTodat+"model2.2/PlotData/"+runNameSave+"-u"+prefix+".plt3d","w")
    filed=open(pathTodat+"model2.2/PlotData/"+runNameSave+"-d"+prefix+".plt3d","w")
    files=open(pathTodat+"model2.2/PlotData/"+runNameSave+"-s"+prefix+".plt3d","w")
    fileubar=open(pathTodat+"model2.2/PlotData/"+runNameSave+"-ubar"+prefix+".plt3d","w")
    filedbar=open(pathTodat+"model2.2/PlotData/"+runNameSave+"-dbar"+prefix+".plt3d","w")
    filesbar=open(pathTodat+"model2.2/PlotData/"+runNameSave+"-sbar"+prefix+".plt3d","w")
    
    i=0
    for x in xValues:
        for b in bValues:    
           fileu.write(str([x,b,meanU[i],stdU[i]])+"\n") 
           filed.write(str([x,b,meanD[i],stdD[i]])+"\n") 
           files.write(str([x,b,meanS[i],stdS[i]])+"\n") 
           fileubar.write(str([x,b,meanUbar[i],stdUbar[i]])+"\n") 
           filedbar.write(str([x,b,meanDbar[i],stdDbar[i]])+"\n") 
           filesbar.write(str([x,b,meanSbar[i],stdSbar[i]])+"\n") 
           i+=1
    
    fileu.close()
    filed.close()
    files.close()
    fileubar.close()
    filedbar.close()
    filesbar.close()
    print('-------------------------SAVED!-----------------------------')    

#%%
#######################################
# pass though each perplica and compute grid
#######################################
mu=-1.

u=[]
d=[]
s=[]
ubar=[]
dbar=[]
sbar=[]

for i in range(len(replicas)):
    print('-------------------------',i,'/',len(replicas),'-----------------------------')
    
    #harpy.setNPparameters(replicas[i][4])    
    harpy.setNPparameters([2.000, 0.036, 0.089, 8.657,  0.389, 0.000, 0.465, 6.195, 527.903])
    harpy.setPDFreplica(replicas[i][3])        
    
    u2=[]
    d2=[]
    s2=[]
    ubar2=[]
    dbar2=[]
    sbar2=[]
    
    for x in xValues:
        for b in bValues:
            tmd=harpy.get_uTMDPDF(x,b,1,mu=mu)
            u2.append(tmd[5+2])
            d2.append(tmd[5+1])
            s2.append(tmd[5+3])
            ubar2.append(tmd[5-2])
            dbar2.append(tmd[5-1])
            sbar2.append(tmd[5-3])
            
    u.append(u2)
    d.append(d2)
    s.append(s2)
    ubar.append(ubar2)
    dbar.append(dbar2)
    sbar.append(sbar2)
    
    if((i % 50)==0): SaveResults()

SaveResults()
print('-------------------------DONE!-----------------------------')

