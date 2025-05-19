#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 20:15:46 2019

Program that parse various DY data files to ADP-frendly formal

@author: vla18041
"""

import os
import sys
import numpy
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))+"/"

sys.path.append(ROOT_DIR+"/DataProcessor/")

import DataProcessor.Point
import DataProcessor.DataSet

path_to_save=ROOT_DIR+"/DataProcessor/DataLib/D2_moment/"

M_proton=0.938
#%%
print("Making D2 by Regenburg...")

DataCurrent=DataProcessor.DataSet.DataSet('RQCD_d2_ud',"D2")
DataCurrent.comment="Taken from table 4 for u and d quarks"
DataCurrent.reference="arXiv:2111.08306"

DataCurrent.isNormalized=False
Q_current=[1.95,2.05]
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(0))
p["process"]=2
p["Q"]=Q_current
p["<Q>"]=2.0
p["s"]=4.0
p["xSec"]=0.026
p["uncorrErr"]=[0.004,0.013]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)


p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(1))
p["process"]=1
p["Q"]=Q_current
p["<Q>"]=2.0
p["s"]=4.0
p["xSec"]=-0.0086
p["uncorrErr"]=[0.0026,0.0146]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%

print("Making D2 by Regenburg...")

DataCurrent=DataProcessor.DataSet.DataSet('RQCD_d2_singlet',"D2")
DataCurrent.comment="Taken from table 4 for u-d and u+d quarks"
DataCurrent.reference="arXiv:2111.08306"

DataCurrent.isNormalized=False
Q_current=[1.95,2.05]
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(0))
p["process"]=11
p["Q"]=Q_current
p["<Q>"]=2.0
p["s"]=4.0
p["xSec"]=0.034
p["uncorrErr"]=[0.004,0.011]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)


p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(1))
p["process"]=12
p["Q"]=Q_current
p["<Q>"]=2.0
p["s"]=4.0
p["xSec"]=0.018
p["uncorrErr"]=[0.005,0.022]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%

print("Making D2 by Regenburg...")

DataCurrent=DataProcessor.DataSet.DataSet('RQCD_d2_pn',"D2")
DataCurrent.comment="Taken from table 4 for p and n combination (better not to use because there are no sea-part)"
DataCurrent.reference="arXiv:2111.08306"

DataCurrent.isNormalized=False
Q_current=[1.95,2.05]
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(0))
p["process"]=100
p["Q"]=Q_current
p["<Q>"]=2.0
p["s"]=4.0
p["xSec"]=0.0105
p["uncorrErr"]=[0.019,0.065]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)


p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(1))
p["process"]=101
p["Q"]=Q_current
p["<Q>"]=2.0
p["s"]=4.0
p["xSec"]=-0.0009
p["uncorrErr"]=[0.0014,0.0069]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%

print("Making SANE...")

DataCurrent=DataProcessor.DataSet.DataSet('SANE18',"D2")
DataCurrent.comment="Taken from table 1 (they have factor -1 in definition, see eqn.(2); Q-range is guessed by us)"
DataCurrent.reference="arXiv:1805.08835"

DataCurrent.isNormalized=False
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(0))
p["process"]=100
p["Q"]=[numpy.sqrt(2.),numpy.sqrt(3.5)]
p["<Q>"]=numpy.sqrt(2.8)
p["s"]=2*M_proton*4.7+M_proton**2
p["xSec"]=-0.00414
p["uncorrErr"]=[0.00205,0.00256]
p["corrErr"]=[]
p["thFactor"]=-1.
DataCurrent.AddPoint(p)

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(1))
p["process"]=100
p["Q"]=[numpy.sqrt(3.5),numpy.sqrt(5.)]
p["<Q>"]=numpy.sqrt(4.3)
p["s"]=2*M_proton*5.9+M_proton**2
p["xSec"]=-0.00149
p["uncorrErr"]=[0.00156,0.00368]
p["corrErr"]=[]
p["thFactor"]=-1.
DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")