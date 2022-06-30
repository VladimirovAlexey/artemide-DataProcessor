#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 17:23:04 2022

@author: alexey
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 20:15:46 2019

Program that parse various DY data files to ADP-frendly formal

@author: vla18041
"""
import os
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/"

import sys
sys.path.append(ROOT_DIR+"DataProcessor/")

import DataProcessor.Point
import DataProcessor.DataSet
import numpy

path_to_data=ROOT_DIR+"data/"
path_to_save=ROOT_DIR+"DataProcessor/DataLib/unpolDY/"

M_Z=91.### mass of Z-boson

#%%
###############################################################################
###########################PHENIX 200##########################################
print("Read PHENIX 200 file ...")
#f = open(path_to_data+"Phenix/fig33_data_adapted.txt")
f = open(path_to_data+"Phenix/fig33.txt")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:2]
print('First line =',data_from_f[0])
#print 'last line (before)=',data_from_f[-1]
#del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]


for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split(" ")
    #del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('PHE200',"DY")
DataCurrent.comment="PHENIX 200GeV data"
DataCurrent.reference="arXiv:1805.02448"

DataCurrent.isNormalized=False
proc_current=[1,1,5]
s_current=200.**2
Q_current=[4.8,8.2]
y_current=[1.2,2.2]
lumUncertainty=0.12
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][0]-0.25,data_from_f[i][0]+0.25]
    # 0.31.. is for 1/pi from invariant cross-secX, 0.001 for nb
    # factor 2. is from the symetrization for y, is compensated by 2 bins in y.
    # 0.001 for nb->pb
    p["thFactor"]=1./(p["qT"][1]**2-p["qT"][0]**2)/(y_current[1]-y_current[0])*0.3183098861838*0.001
    print(p["thFactor"])
    p["Q"]=Q_current
    p["y"]=y_current
    p["xSec"]=data_from_f[i][1]
    p["includeCuts"]=False
    p["uncorrErr"].append(data_from_f[i][2])
    p["uncorrErr"].append((data_from_f[i][3]+data_from_f[i][4])/2.)
    p["corrErr"].append(data_from_f[i][5])
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################STAR##########################################
print("Read STAR file ...")
f = open(path_to_data+"STAR-DY/STAR_data.txt")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:2]
print('First line =',data_from_f[0])
#print 'last line (before)=',data_from_f[-1]
#del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]




for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("  ")
    #del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('STAR510',"DY")
DataCurrent.comment="STAR DY data"
DataCurrent.reference="preliminary"

bins=[0., 1.25, 2.5, 3.75, 5., 7.5, 10., 12.5, 15., 17.5, 20., 25.]

DataCurrent.isNormalized=False
proc_current=[1,1,5]
s_current=510.**2
Q_current=[73.,114.]
y_current=[-1.,1.]
cutParam=[25.,25.,-1.,1.]
lumUncertainty=0.085
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[bins[i],bins[i+1]]
    p["<qT>"]=data_from_f[i][0]
    # 0.001 for nb->pb
    p["thFactor"]=1./(p["qT"][1]-p["qT"][0])
    p["Q"]=Q_current
    p["<Q>"]=M_Z
    p["y"]=y_current
    p["xSec"]=data_from_f[i][1]
    p["includeCuts"]=True
    p["cutParams"]=cutParam
    p["uncorrErr"].append(data_from_f[i][2])
    p["uncorrErr"].append((data_from_f[i][3]+data_from_f[i][4])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")