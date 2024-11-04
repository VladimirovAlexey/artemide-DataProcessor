#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 20:15:46 2019

Program that parse various DY data files to ADP-frendly formal

@author: vla18041
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/DataProcessor/")

import numpy

import DataProcessor.Point
import DataProcessor.DataSet

path_to_data="/data/arTeMiDe_Repository/data/"
path_to_save="/data/arTeMiDe_Repository/DataProcessor/DataLib/unpolDY/"

M_Z=91.1876### mass of Z-boson

#%%
### pt bins are restored by from the covaraince matrix file
ptBins=[[0., 1.], [1., 2.], [2., 3.], [3., 4.], [4., 5.], [5., 6.], [6., 
  7.], [7., 8.], [8., 9.], [9., 10.], [10., 11.], [11., 12.], [12., 
  13.], [13., 14.], [14., 16.], [16., 18.], [18., 20.], [20., 
  22.], [22., 25.], [25., 28.], [28., 32.], [32., 37.], [37., 
  43.], [43., 52.], [52., 65.], [65., 85.], [85., 120.], [120., 
  160.], [160., 190.], [190., 220.], [220., 250.], [250., 
  300.], [300., 400.], [400., 1500.]]

def FindBin(pt_in):
    for b in ptBins:
        if(b[0]<=pt_in and pt_in<b[1]):
            return b
    raise Exception("BIN is not found")

#%%
### given in the text
proc_current=[1,1,1,3]
s_current=13000.**2
Q_current=[M_Z-15.,M_Z+15.]
incCut=True
cutParam=[25.,25.,-2.4,2.4]  # begining of sec.6
lumUncertainty=0.025 # table 2
corrSys=0.001 #### They take Background estimations as fully correlated (table 2)

#%%
###############################################################################
###########################CMS 13 y-diff########################################
print("Read CMS13 ydiff file...")
f = open(path_to_data+"CMS/CMS13-ydiff-abs-born.csv")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")
#%%
data_here=data_from_f[6:39]
for i in range(len(data_here)):
    data_here[i]=data_here[i].split(",")
    data_here[i]=[float(j) for j in data_here[i]]
 

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('CMS13-00y04',"DY")
DataCurrent.comment="CMS 13TeV 0.0<|y|<0.4 absolute"
DataCurrent.reference="arXiv:1909.04133"

DataCurrent.isNormalized=False
y_current=[0.,0.4]
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_here)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=FindBin(data_here[i][0])
    p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_here[i][1]
    ###### The uncertanty in table is total including correlated
    ###### I take it and subtract (in quadrature) correlated uncertainty. Result almost coincides with tables
    ###### The difference between result and tables is that tables for dressed, and only dimuon or di-electron
    ###### main differeence takes place at large pT
    unc=(data_here[i][2]-data_here[i][3])/2.
    p["uncorrErr"].append(numpy.sqrt(unc**2-(p["xSec"]*lumUncertainty)**2))
    p["corrErr"].append(p["xSec"]*corrSys)
    DataCurrent.AddPoint(p)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
data_here=data_from_f[41:74]
for i in range(len(data_here)):
    data_here[i]=data_here[i].split(",")
    data_here[i]=[float(j) for j in data_here[i]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('CMS13-04y08',"DY")
DataCurrent.comment="CMS 13TeV 0.4<|y|<0.8 absolute"
DataCurrent.reference="arXiv:1909.04133"

DataCurrent.isNormalized=False
y_current=[0.4,0.8]
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_here)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=FindBin(data_here[i][0])
    p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_here[i][1]
    unc=(data_here[i][2]-data_here[i][3])/2.
    p["uncorrErr"].append(numpy.sqrt(unc**2-(p["xSec"]*lumUncertainty)**2))
    p["corrErr"].append(p["xSec"]*corrSys)
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
data_here=data_from_f[76:109]
for i in range(len(data_here)):
    data_here[i]=data_here[i].split(",")
    data_here[i]=[float(j) for j in data_here[i]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('CMS13-08y12',"DY")
DataCurrent.comment="CMS 13TeV 0.8<|y|<1.2 absolute"
DataCurrent.reference="arXiv:1909.04133"

DataCurrent.isNormalized=False
y_current=[0.8,1.2]
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_here)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=FindBin(data_here[i][0])
    p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_here[i][1]
    unc=(data_here[i][2]-data_here[i][3])/2.
    p["uncorrErr"].append(numpy.sqrt(unc**2-(p["xSec"]*lumUncertainty)**2))
    p["corrErr"].append(p["xSec"]*corrSys)
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
data_here=data_from_f[111:144]
for i in range(len(data_here)):
    data_here[i]=data_here[i].split(",")
    data_here[i]=[float(j) for j in data_here[i]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('CMS13-12y16',"DY")
DataCurrent.comment="CMS 13TeV 1.2<|y|<1.6 absolute"
DataCurrent.reference="arXiv:1909.04133"

DataCurrent.isNormalized=False
y_current=[1.2,1.6]
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_here)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=FindBin(data_here[i][0])
    p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_here[i][1]
    unc=(data_here[i][2]-data_here[i][3])/2.
    p["uncorrErr"].append(numpy.sqrt(unc**2-(p["xSec"]*lumUncertainty)**2))
    p["corrErr"].append(p["xSec"]*corrSys)
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%
data_here=data_from_f[146:179]
for i in range(len(data_here)):
    data_here[i]=data_here[i].split(",")
    data_here[i]=[float(j) for j in data_here[i]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('CMS13-16y24',"DY")
DataCurrent.comment="CMS 13TeV 1.6<|y|<2.4 absolute"
DataCurrent.reference="arXiv:1909.04133"

DataCurrent.isNormalized=False
y_current=[1.6,2.4]
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_here)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=FindBin(data_here[i][0])
    p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_here[i][1]
    unc=(data_here[i][2]-data_here[i][3])/2.
    p["uncorrErr"].append(numpy.sqrt(unc**2-(p["xSec"]*lumUncertainty)**2))
    p["corrErr"].append(p["xSec"]*corrSys)
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################CMS 13 y-diff NORMALIZED########################################
print("Read CMS13 ydiff norm file...")
f = open(path_to_data+"CMS/CMS13-ydiff-norm-born.csv")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")
#%%
data_here=data_from_f[6:39]
for i in range(len(data_here)):
    data_here[i]=data_here[i].split(",")
    data_here[i]=[float(j) for j in data_here[i]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('CMS13-00y04-norm',"DY")
DataCurrent.comment="CMS 13TeV 0.0<|y|<0.4 normalized"
DataCurrent.reference="arXiv:1909.04133"

DataCurrent.isNormalized=True
y_current=[0.,0.4]
#DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_here)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=FindBin(data_here[i][0])
    p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_here[i][1]
    unc=(data_here[i][2]-data_here[i][3])/2.
    p["uncorrErr"].append(unc)
    p["corrErr"].append(p["xSec"]*corrSys)
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
data_here=data_from_f[41:74]
for i in range(len(data_here)):
    data_here[i]=data_here[i].split(",")
    data_here[i]=[float(j) for j in data_here[i]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('CMS13-04y08-norm',"DY")
DataCurrent.comment="CMS 13TeV 0.4<|y|<0.8 normalized"
DataCurrent.reference="arXiv:1909.04133"

DataCurrent.isNormalized=True
y_current=[0.4,0.8]
#DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_here)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=FindBin(data_here[i][0])
    p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_here[i][1]
    unc=(data_here[i][2]-data_here[i][3])/2.
    p["uncorrErr"].append(unc)
    p["corrErr"].append(p["xSec"]*corrSys)
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
data_here=data_from_f[76:109]
for i in range(len(data_here)):
    data_here[i]=data_here[i].split(",")
    data_here[i]=[float(j) for j in data_here[i]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('CMS13-08y12-norm',"DY")
DataCurrent.comment="CMS 13TeV 0.8<|y|<1.2 normalized"
DataCurrent.reference="arXiv:1909.04133"

DataCurrent.isNormalized=True
y_current=[0.8,1.2]
#DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_here)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=FindBin(data_here[i][0])
    p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_here[i][1]
    unc=(data_here[i][2]-data_here[i][3])/2.
    p["uncorrErr"].append(unc)
    p["corrErr"].append(p["xSec"]*corrSys)
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
data_here=data_from_f[111:144]
for i in range(len(data_here)):
    data_here[i]=data_here[i].split(",")
    data_here[i]=[float(j) for j in data_here[i]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('CMS13-12y16-norm',"DY")
DataCurrent.comment="CMS 13TeV 1.2<|y|<1.6 normalized"
DataCurrent.reference="arXiv:1909.04133"

DataCurrent.isNormalized=True
y_current=[1.2,1.6]
#DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_here)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=FindBin(data_here[i][0])
    p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_here[i][1]
    unc=(data_here[i][2]-data_here[i][3])/2.
    p["uncorrErr"].append(unc)
    p["corrErr"].append(p["xSec"]*corrSys)
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%
data_here=data_from_f[146:179]
for i in range(len(data_here)):
    data_here[i]=data_here[i].split(",")
    data_here[i]=[float(j) for j in data_here[i]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('CMS13-16y24-norm',"DY")
DataCurrent.comment="CMS 13TeV 1.6<|y|<2.4 normalized"
DataCurrent.reference="arXiv:1909.04133"

DataCurrent.isNormalized=True
y_current=[1.6,2.4]
#DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_here)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=FindBin(data_here[i][0])
    p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_here[i][1]
    unc=(data_here[i][2]-data_here[i][3])/2.
    p["uncorrErr"].append(unc)
    p["corrErr"].append(p["xSec"]*corrSys)
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")