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

import DataProcessor.ArtemideReplicaSet

path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/ART25_main.atmde"
harpy.initialize(path_to_constants)

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                        "/data/WorkingFiles/TMD/Fit_Notes/ART25/REPLICAS/ART25_run1.rep")
    
rSet.SetReplica(0)


#%%
xRange=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4]

mu=1.5
zeta=1.5**2
path_save="/data/WorkingFiles/TMD/Fit_Notes/Predictions/ART25_for_Copeland/ART25_dminusU_mu=zeta=15.csv"
path_saveBAR="/data/WorkingFiles/TMD/Fit_Notes/Predictions/ART25_for_Copeland/ART25_dBminusUB_mu=zeta=15.csv"

TMDresult=[]
TMDresultBar=[]

for x in xRange:
    
    TMDreplicas=[]
    TMDreplicasBar=[]
    for r in range(rSet.numberOfReplicas):
    #for r in range(10):
        rSet.SetReplica(r)
        
        ww=[]
        wwBar=[]
        for j in range(100):
            b=0.1*j+0.1            
            tt=harpy.get_uTMDPDF(x,b,1,mu,zeta)
            wwBar.append(tt[4]-tt[3]) ### dBar-uBar
            ww.append(tt[6]-tt[7])      ### d-u
        
        TMDreplicas.append(ww)
        TMDreplicasBar.append(wwBar)
        
    mm=numpy.mean(TMDreplicas,axis=0)
    ss=numpy.std(TMDreplicas,axis=0)
    TMDresult.append(mm)
    TMDresult.append(ss)
    
    mm=numpy.mean(TMDreplicasBar,axis=0)
    ss=numpy.std(TMDreplicasBar,axis=0)
    TMDresultBar.append(mm)
    TMDresultBar.append(ss)

with open(path_saveBAR, 'w') as file:
    file.write('dBar-uBar TMDPDF in b-space from 2503.11201 mu=sqrt[zeta]=1.5 \n')
    file.write('uncertanties are symmetric for simplicity \n')
    file.write('b(GeV), x=0.05(mean), x=0.05(std), \
               x=0.1(mean), x=0.1(std), \
               x=0.15(mean), x=0.15(std), \
                x=0.2(mean), x=0.2(std), \
                x=0.25(mean), x=0.25(std), \
                x=0.3(mean), x=0.3(std), \
                x=0.35(mean), x=0.35(std), \
                x=0.4(mean), x=0.4(std), \n')
    for j in range(100):
        b=0.1*j+0.1   
        toSave=numpy.insert(numpy.array(TMDresultBar)[:,j],0,b)
        
        file.write(', '.join("{:.8f}".format(item) for item in toSave)+'\n')

with open(path_save, 'w') as file:
    file.write('d-u TMDPDF in b-space from 2503.11201 mu=sqrt[zeta]=1.5 \n')
    file.write('uncertanties are symmetric for simplicity \n')
    file.write('b(GeV), x=0.05(mean), x=0.05(std), \
               x=0.1(mean), x=0.1(std), \
               x=0.15(mean), x=0.15(std), \
                x=0.2(mean), x=0.2(std), \
                x=0.25(mean), x=0.25(std), \
                x=0.3(mean), x=0.3(std), \
                x=0.35(mean), x=0.35(std), \
                x=0.4(mean), x=0.4(std), \n')
    for j in range(100):
        b=0.1*j+0.1   
        toSave=numpy.insert(numpy.array(TMDresult)[:,j],0,b)
        
        file.write(', '.join("{:.8f}".format(item) for item in toSave)+'\n')
        
#%%
xRange=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4]

mu=7.3
zeta=7.3**2
path_save="/data/WorkingFiles/TMD/Fit_Notes/Predictions/ART25_for_Copeland/ART25_dminusU_mu=zeta=73.csv"
path_saveBAR="/data/WorkingFiles/TMD/Fit_Notes/Predictions/ART25_for_Copeland/ART25_dBminusUB_mu=zeta=73.csv"

TMDresult=[]
TMDresultBar=[]

for x in xRange:
    
    TMDreplicas=[]
    TMDreplicasBar=[]
    for r in range(rSet.numberOfReplicas):
    #for r in range(10):
        rSet.SetReplica(r)
        
        ww=[]
        wwBar=[]
        for j in range(100):
            b=0.1*j+0.1            
            tt=harpy.get_uTMDPDF(x,b,1,mu,zeta)
            wwBar.append(tt[4]-tt[3]) ### dBar-uBar
            ww.append(tt[6]-tt[7])      ### d-u
        
        TMDreplicas.append(ww)
        TMDreplicasBar.append(wwBar)
        
    mm=numpy.mean(TMDreplicas,axis=0)
    ss=numpy.std(TMDreplicas,axis=0)
    TMDresult.append(mm)
    TMDresult.append(ss)
    
    mm=numpy.mean(TMDreplicasBar,axis=0)
    ss=numpy.std(TMDreplicasBar,axis=0)
    TMDresultBar.append(mm)
    TMDresultBar.append(ss)

with open(path_saveBAR, 'w') as file:
    file.write('dBar-uBar TMDPDF in b-space from 2503.11201 mu=sqrt[zeta]=1.5 \n')
    file.write('uncertanties are symmetric for simplicity \n')
    file.write('b(GeV), x=0.05(mean), x=0.05(std), \
               x=0.1(mean), x=0.1(std), \
               x=0.15(mean), x=0.15(std), \
                x=0.2(mean), x=0.2(std), \
                x=0.25(mean), x=0.25(std), \
                x=0.3(mean), x=0.3(std), \
                x=0.35(mean), x=0.35(std), \
                x=0.4(mean), x=0.4(std), \n')
    for j in range(100):
        b=0.1*j+0.1   
        toSave=numpy.insert(numpy.array(TMDresultBar)[:,j],0,b)
        
        file.write(', '.join("{:.8f}".format(item) for item in toSave)+'\n')

with open(path_save, 'w') as file:
    file.write('d-u TMDPDF in b-space from 2503.11201 mu=sqrt[zeta]=1.5 \n')
    file.write('uncertanties are symmetric for simplicity \n')
    file.write('b(GeV), x=0.05(mean), x=0.05(std), \
               x=0.1(mean), x=0.1(std), \
               x=0.15(mean), x=0.15(std), \
                x=0.2(mean), x=0.2(std), \
                x=0.25(mean), x=0.25(std), \
                x=0.3(mean), x=0.3(std), \
                x=0.35(mean), x=0.35(std), \
                x=0.4(mean), x=0.4(std), \n')
    for j in range(100):
        b=0.1*j+0.1   
        toSave=numpy.insert(numpy.array(TMDresult)[:,j],0,b)
        
        file.write(', '.join("{:.8f}".format(item) for item in toSave)+'\n')

#%%
xRange=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4]

mu=7.3
zeta=1.5**2
path_save="/data/WorkingFiles/TMD/Fit_Notes/Predictions/ART25_for_Copeland/ART25_dminusU_mu=73_zeta=15.csv"
path_saveBAR="/data/WorkingFiles/TMD/Fit_Notes/Predictions/ART25_for_Copeland/ART25_dBminusUB_mu=73_zeta=15.csv"

TMDresult=[]
TMDresultBar=[]

for x in xRange:
    
    TMDreplicas=[]
    TMDreplicasBar=[]
    for r in range(rSet.numberOfReplicas):
    #for r in range(10):
        rSet.SetReplica(r)
        
        ww=[]
        wwBar=[]
        for j in range(100):
            b=0.1*j+0.1            
            tt=harpy.get_uTMDPDF(x,b,1,mu,zeta)
            wwBar.append(tt[4]-tt[3]) ### dBar-uBar
            ww.append(tt[6]-tt[7])      ### d-u
        
        TMDreplicas.append(ww)
        TMDreplicasBar.append(wwBar)
        
    mm=numpy.mean(TMDreplicas,axis=0)
    ss=numpy.std(TMDreplicas,axis=0)
    TMDresult.append(mm)
    TMDresult.append(ss)
    
    mm=numpy.mean(TMDreplicasBar,axis=0)
    ss=numpy.std(TMDreplicasBar,axis=0)
    TMDresultBar.append(mm)
    TMDresultBar.append(ss)

with open(path_saveBAR, 'w') as file:
    file.write('dBar-uBar TMDPDF in b-space from 2503.11201 mu=sqrt[zeta]=1.5 \n')
    file.write('uncertanties are symmetric for simplicity \n')
    file.write('b(GeV), x=0.05(mean), x=0.05(std), \
               x=0.1(mean), x=0.1(std), \
               x=0.15(mean), x=0.15(std), \
                x=0.2(mean), x=0.2(std), \
                x=0.25(mean), x=0.25(std), \
                x=0.3(mean), x=0.3(std), \
                x=0.35(mean), x=0.35(std), \
                x=0.4(mean), x=0.4(std), \n')
    for j in range(100):
        b=0.1*j+0.1   
        toSave=numpy.insert(numpy.array(TMDresultBar)[:,j],0,b)
        
        file.write(', '.join("{:.8f}".format(item) for item in toSave)+'\n')

with open(path_save, 'w') as file:
    file.write('d-u TMDPDF in b-space from 2503.11201 mu=sqrt[zeta]=1.5 \n')
    file.write('uncertanties are symmetric for simplicity \n')
    file.write('b(GeV), x=0.05(mean), x=0.05(std), \
               x=0.1(mean), x=0.1(std), \
               x=0.15(mean), x=0.15(std), \
                x=0.2(mean), x=0.2(std), \
                x=0.25(mean), x=0.25(std), \
                x=0.3(mean), x=0.3(std), \
                x=0.35(mean), x=0.35(std), \
                x=0.4(mean), x=0.4(std), \n')
    for j in range(100):
        b=0.1*j+0.1   
        toSave=numpy.insert(numpy.array(TMDresult)[:,j],0,b)
        
        file.write(', '.join("{:.8f}".format(item) for item in toSave)+'\n')