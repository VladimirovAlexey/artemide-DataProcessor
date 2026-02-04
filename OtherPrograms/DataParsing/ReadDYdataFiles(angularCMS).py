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

path_to_data=os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))+'/data/CMS/1504.03512/'
path_to_save=ROOT_DIR+"/DataLib/DY_angular/"
 
M_Z=91.1876### mass of Z-boson

#%%
### given in the text
#proc_current=[1,1,5]
s_current=8000.**2
Q_current=[81.,101.]
incCut=False
cutParam=[20.,20.,2.,4.5]  

#%%
def FillTable(datain):
    
    for i in range(len(datain)):
        datain[i]=datain[i].split(",")
        datain[i]=[float(j) for j in datain[i]]
    
    print("Done.  =>     Create points & append to data set ...")
       
    DataCurrent.reference="arXiv:1504.03512"
    
    DataCurrent.isNormalized=False
    
    for i in range(len(datain)):
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=[datain[i][1],datain[i][2]]
        p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
        p["Q"]=Q_current
        p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
        p["y"]=y_current
        p["includeCuts"]=incCut
        p["cutParams"]=cutParam
        p["xSec"]=datain[i][3]
        unc1=(datain[i][4]-datain[i][5])/2.
        unc2=(datain[i][6]-datain[i][7])/2.
        p["uncorrErr"].append(unc1)
        p["uncorrErr"].append(unc2)
        DataCurrent.AddPoint(p)    
    
    print("Done.  ")

#%%
###############################################################################

y_current=[0.,1.0]

print("Read CMS13 0<y<1 file...")
f = open(path_to_data+"y0to1.csv")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")
#%%

data_here=data_from_f[15:23]

proc_current=[1,1,1,3]
DataCurrent=DataProcessor.DataSet.DataSet('CMS8-Auu-00y10',"DY")
DataCurrent.comment="CMS 8TeV 0.0<|y|<1.0 Auu (just for normalization)"
FillTable(data_here)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%

data_here=data_from_f[15:23]

proc_current=[1,1,1,20]
DataCurrent=DataProcessor.DataSet.DataSet('CMS8-A0-00y10',"DY")
DataCurrent.comment="CMS 8TeV 0.0<|y|<1.0 A0"
FillTable(data_here)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%

data_here=data_from_f[30:38]

proc_current=[1,1,1,21]
DataCurrent=DataProcessor.DataSet.DataSet('CMS8-A1-00y10',"DY")
DataCurrent.comment="CMS 8TeV 0.0<|y|<1.0 A1"
FillTable(data_here)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%

data_here=data_from_f[45:53]

proc_current=[1,1,1,22]
DataCurrent=DataProcessor.DataSet.DataSet('CMS8-A2-00y10',"DY")
DataCurrent.comment="CMS 8TeV 0.0<|y|<1.0 A2"
FillTable(data_here)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%

data_here=data_from_f[60:68]

proc_current=[1,1,1,23]
DataCurrent=DataProcessor.DataSet.DataSet('CMS8-A3-00y10',"DY")
DataCurrent.comment="CMS 8TeV 0.0<|y|<1.0 A3"
FillTable(data_here)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%

data_here=data_from_f[75:83]

proc_current=[1,1,1,24]
DataCurrent=DataProcessor.DataSet.DataSet('CMS8-A4-00y10',"DY")
DataCurrent.comment="CMS 8TeV 0.0<|y|<1.0 A4"
FillTable(data_here)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################

y_current=[1.,2.1]

print("Read CMS13 0<y<1 file...")
f = open(path_to_data+"y1to21.csv")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")
#%%

data_here=data_from_f[15:23]

proc_current=[1,1,1,3]
DataCurrent=DataProcessor.DataSet.DataSet('CMS8-Auu-10y21',"DY")
DataCurrent.comment="CMS 8TeV 0.0<|y|<1.0 Auu (just for normalization)"
FillTable(data_here)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%

data_here=data_from_f[15:23]

proc_current=[1,1,1,20]
DataCurrent=DataProcessor.DataSet.DataSet('CMS8-A0-10y21',"DY")
DataCurrent.comment="CMS 8TeV 1.0<|y|<2.1 A0"
FillTable(data_here)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%

data_here=data_from_f[30:38]

proc_current=[1,1,1,21]
DataCurrent=DataProcessor.DataSet.DataSet('CMS8-A1-10y21',"DY")
DataCurrent.comment="CMS 8TeV 1.0<|y|<2.1 A1"
FillTable(data_here)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%

data_here=data_from_f[45:53]

proc_current=[1,1,1,22]
DataCurrent=DataProcessor.DataSet.DataSet('CMS8-A2-10y21',"DY")
DataCurrent.comment="CMS 8TeV 1.0<|y|<2.1 A2"
FillTable(data_here)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%

data_here=data_from_f[60:68]

proc_current=[1,1,1,23]
DataCurrent=DataProcessor.DataSet.DataSet('CMS8-A3-10y21',"DY")
DataCurrent.comment="CMS 8TeV 1.0<|y|<2.1 A3"
FillTable(data_here)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%

data_here=data_from_f[75:83]

proc_current=[1,1,1,24]
DataCurrent=DataProcessor.DataSet.DataSet('CMS8-A4-10y21',"DY")
DataCurrent.comment="CMS 8TeV 1.0<|y|<2.1 A4"
FillTable(data_here)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")