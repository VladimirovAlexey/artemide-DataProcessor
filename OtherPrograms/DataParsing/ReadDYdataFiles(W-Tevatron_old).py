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

import numpy

import DataProcessor.Point
import DataProcessor.DataSet

path_to_data=ROOT_DIR+"data/"
path_to_save=ROOT_DIR+"/DataProcessor/DataLib/unpolW/"

M_W=80.### mass of Z-boson

#%%
################################################################################
print("Read D0 W-boson file...")

import json
 
with open(path_to_data+"CDF_D0/D0_W-98.json") as json_file:
    data_from_f = json.load(json_file)
    
#%%

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('D0run1-W',"DY")
DataCurrent.comment="D0 1.8TeV normalized, W-boson production to electron"
DataCurrent.reference="arXiv:hep-ph/9803003"

DataCurrent.isNormalized=True
#### not clear is it 
proc_current=[1,1,12]
s_current=1800.**2
y_current=[-1.1,1.1]
incCut=True
#### since the put restriction only on electron rapidity I do not restrict rapidity at all (it is done by y)
cutParam=[25.,25.,-20.1,20.1]
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
    #### the lower boundary is defied by M^2>E1+E2
    #### checked that upper boundrary 300 is enough
    p["Q"]=[numpy.sqrt(50.**2-p["qT"][1]**2) if p["qT"][1]<50. else 2., 300.]
    if(p["Q"][0]<10.): p["Q"][0]=10.
    p["<Q>"]=M_W ### to be sure that mass middle value is W-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=float(p_from_d['y'][0]['value'])    
    p["uncorrErr"].append(float(p_from_d['y'][0]['errors'][0]['symerror']))
    p["uncorrErr"].append(float(p_from_d['y'][0]['errors'][1]['symerror']))
    p["uncorrErr"].append(float(p_from_d['y'][0]['errors'][2]['symerror']))
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
################################################################################
print("Read CDF W-boson file...")

import json
 
with open(path_to_data+"CDF_D0/CDF_W-91.json") as json_file:
    data_from_f = json.load(json_file)
    
#%%

def qTbin(x):
    if(x<20):
        return [x-1.,x+1.]
    elif(x<50):
        return [x-2.5,x+2.5]
    elif(x<60):#### the bining for hign-qT bins is assumed from the plot
        return [50.,60.]
    elif(x<80):
        return [60.,80.]
    elif(x<130):
        return [80.,130.]
    else:
        return [130.,200.]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('CDFrun1-W',"DY")
DataCurrent.comment="CDF 1.8TeV, W-boson production to electron"
DataCurrent.reference="Phys.Rev.Lett. 66 (1991)"

DataCurrent.isNormalized=True

#### not clear is it 
proc_current=[1,1,12]
s_current=1800.**2
## no restriction on y
y_current=[-100.,100.]
incCut=False
#### since the put restriction only on electron rapidity I do not restrict rapidity at all (it is done by y)
cutParam=[20.,20.,-1.1,1.1]
lumUncertainty=0.069
##DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f["values"])):
    
    p_from_d=data_from_f["values"][i]
    
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["<qT>"]=float(p_from_d["x"][0]["value"])
    p["qT"]=qTbin(p["<qT>"])
    p["thFactor"]=1./(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
    #### the lower boundary is no defined
    p["Q"]=[40., 300.]
    p["<Q>"]=M_W ### to be sure that mass middle value is W-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=float(p_from_d['y'][0]['value'])    
    p["uncorrErr"].append(float(p_from_d['y'][0]['errors'][0]['symerror']))
    p["uncorrErr"].append(float(p_from_d['y'][0]['errors'][1]['symerror']))
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")