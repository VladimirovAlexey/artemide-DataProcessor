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
path_to_save=ROOT_DIR+"/DataProcessor/DataLib/unpolW/"

M_W=80.### mass of Z-boson

#%%
###############################################################################
###########################ATLAS 8 00y04########################################
print("Read CMS W-boson-electron file...")

import json
 
with open(path_to_data+"CMS/WbosonELECTRON_2016.json") as json_file:
    data_from_f = json.load(json_file)
    
#%%

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('CMS8-W-electron',"DY")
DataCurrent.comment="CMS 8TeV normalized, W-boson production to electron"
DataCurrent.reference="arXiv:1606.05864"

DataCurrent.isNormalized=True
proc_current=[1,1,9]
s_current=8000.**2
Q_current=[20.,300.]
y_current=[-2.5,2.5]
incCut=True
cutParam=[25.,0.,-2.5,2.5]
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
    p["<Q>"]=M_W ### to be sure that mass middle value is W-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=float(p_from_d['y'][0]['value'])    
    p["uncorrErr"].append(
        float(p_from_d['y'][0]['errors'][0]['symerror'])
        )
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################ATLAS 8 00y04########################################
print("Read CMS W-boson-electron file...")

import json
 
with open(path_to_data+"CMS/WbosonMUON_2016.json") as json_file:
    data_from_f = json.load(json_file)
    
#%%

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('CMS8-W-muon',"DY")
DataCurrent.comment="CMS 8TeV normalized, W-boson production to muon"
DataCurrent.reference="arXiv:1606.05864"

DataCurrent.isNormalized=True
proc_current=[1,1,9]
s_current=8000.**2
Q_current=[20.,300.]
y_current=[-2.1,2.1]
incCut=True
cutParam=[20.,0.,-2.5,2.5]
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
    p["<Q>"]=M_W ### to be sure that mass middle value is W-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=float(p_from_d['y'][0]['value'])
    p["uncorrErr"].append(
        float(p_from_d['y'][0]['errors'][0]['symerror'])
        )
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")