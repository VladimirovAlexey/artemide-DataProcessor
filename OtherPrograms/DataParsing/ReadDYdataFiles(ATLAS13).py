#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 20:15:46 2019

Program that parse various DY data files to ADP-frendly formal

@author: vla18041
"""

import os
import sys
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))+"/"

sys.path.append(ROOT_DIR+"/DataProcessor/")

import DataProcessor.Point
import DataProcessor.DataSet

path_to_data=ROOT_DIR+"data/"
path_to_save=ROOT_DIR+"/DataProcessor/DataLib/unpolDY/"

M_Z=91.### mass of Z-boson

#%%
###############################################################################
###########################ATLAS 8 00y04########################################
print("Read ATLAS 13 file...")

import json
 
with open(path_to_data+"ATLAS/A13.json") as json_file:
    data_from_f = json.load(json_file)
    
#%%

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('A13-norm',"DY")
DataCurrent.comment="ATLAS 13TeV normalized to 1/sigma"
DataCurrent.reference="arXiv:1912.02844"

DataCurrent.isNormalized=True
proc_current=[1,1,5]
s_current=13000.**2
Q_current=[66.,116.]
y_current=[-2.5,2.5]
incCut=True
cutParam=[27.,27.,-2.5,2.5]
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f["values"])):
    
    p_from_d=data_from_f["values"][i]
    
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[float(p_from_d["x"][0]["low"]),
             float(p_from_d["x"][0]["high"])]
    p["thFactor"]=1./(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=float(p_from_d['y'][0]['value'])
    ## they give the uncertanty in procentage
    p["uncorrErr"].append(
        p["xSec"]*0.01*float(p_from_d['y'][0]['errors'][1]['symerror'][:-1])
        )
    p["corrErr"].append(        
            p["xSec"]*0.01*float(p_from_d['y'][0]['errors'][0]['symerror'][:-1])
            )
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
