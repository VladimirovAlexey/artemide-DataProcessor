#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 20:15:46 2019

Program that parse various DY data files to ADP-frendly formal

@author: vla18041
"""
import os
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))


import sys
sys.path.append(ROOT_DIR)

import DataProcessor.Point
import DataProcessor.DataSet

path_to_data=os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))+'/data/'
path_to_save=ROOT_DIR+"/DataLib/unpolDY/"

M_Z=91.1876### mass of Z-boson

#%%
### given in the text
proc_current=[1,1,1,3]
s_current=13000.**2
Q_current=[60.,120.]
incCut=True
cutParam=[20.,20.,2.,4.5]  # begining of sec.3
lumUncertainty=0.02 # table 1

#%%
###############################################################################
###########################CMS 13 y-differential########################################
print("Read LHcb13 yint file...")
f = open(path_to_data+"LHCb/LHCb_13_dy[2021].dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
data_here=data_from_f[7:]

for i in range(len(data_here)):
    data_here[i]=data_here[i].split(",")
    data_here[i]=[float(j) for j in data_here[i]]
#%%
print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('LHCb13_dy(2021)',"DY")
DataCurrent.comment="LHCb 13TeV 2021 update y-differential"
DataCurrent.reference="arXiv:2112.07458"

DataCurrent.isNormalized=False
y_current=[2.,4.5]
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_here)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=data_here[i][2:4]
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=data_here[i][0:2]
    p["thFactor"]=1/(p["qT"][1]-p["qT"][0])/(p["y"][1]-p["y"][0])#devide by bin size
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_here[i][4]
    ###### The uncertanty is statistic + systematic. 
    ###### Systematic is lightly correlated I ignore it
    p["uncorrErr"].append(data_here[i][5])
    p["uncorrErr"].append(data_here[i][6])
    DataCurrent.AddPoint(p)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################CMS 13 y-integrated########################################
print("Read LHcb13 yint file...")
f = open(path_to_data+"LHCb/LHCb_13[2021].dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
data_here=data_from_f[7:21]

for i in range(len(data_here)):
    data_here[i]=data_here[i].split(",")
    data_here[i]=[float(j) for j in data_here[i]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('LHCb13(2021)',"DY")
DataCurrent.comment="LHCb 13TeV 2021 update"
DataCurrent.reference="arXiv:2112.07458"

DataCurrent.isNormalized=False
y_current=[2.,4.5]
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_here)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=data_here[i][0:2]
    p["thFactor"]=1/(p["qT"][1]-p["qT"][0])#devide by bin size
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_here[i][2]
    ###### The uncertanty is statistic + systematic. 
    ###### Systematic is lightly correlated I ignore it
    p["uncorrErr"].append(data_here[i][3])
    p["uncorrErr"].append(data_here[i][4])
    DataCurrent.AddPoint(p)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")