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
p["thFactor"]=2.
DataCurrent.AddPoint(p)


p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(1))
p["process"]=1
p["Q"]=Q_current
p["<Q>"]=2.0
p["s"]=4.0
p["xSec"]=-0.0086
p["uncorrErr"]=[0.0026,0.0146]
p["corrErr"]=[]
p["thFactor"]=2.
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
p["thFactor"]=2.
DataCurrent.AddPoint(p)


p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(1))
p["process"]=12
p["Q"]=Q_current
p["<Q>"]=2.0
p["s"]=4.0
p["xSec"]=0.018
p["uncorrErr"]=[0.005,0.022]
p["corrErr"]=[]
p["thFactor"]=2.
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
p["thFactor"]=2.
DataCurrent.AddPoint(p)


p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(1))
p["process"]=101
p["Q"]=Q_current
p["<Q>"]=2.0
p["s"]=4.0
p["xSec"]=-0.0009
p["uncorrErr"]=[0.0014,0.0069]
p["corrErr"]=[]
p["thFactor"]=2.
DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")


#%%
print("Making D2 by E143 (1995)...")

DataCurrent=DataProcessor.DataSet.DataSet('E143-1995_d2',"D2")
DataCurrent.comment="Taken from table 5a"
DataCurrent.reference="hep-ex/9511013"

DataCurrent.isNormalized=False
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

s_current=p["s"]=2*29.*M_proton+M_proton**2

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(0))
p["process"]=100
p["<Q>"]=numpy.sqrt(5.0)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=0.0054
p["uncorrErr"]=[0.005]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)


p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(1))
p["process"]=102
p["<Q>"]=numpy.sqrt(5.0)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=0.0039
p["uncorrErr"]=[0.0092]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
print("Making D2 by E143 (1998)...")

DataCurrent=DataProcessor.DataSet.DataSet('E143_d2',"D2")
DataCurrent.comment="Taken from table XXXIII"
DataCurrent.reference="hep-ph/9802357"

DataCurrent.isNormalized=False
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

s_current=p["s"]=2*29.*M_proton+M_proton**2

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(0))
p["process"]=100
p["<Q>"]=numpy.sqrt(5.0)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=0.0058
p["uncorrErr"]=[0.005]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(1))
p["process"]=102
p["<Q>"]=numpy.sqrt(5.)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=0.0051
p["uncorrErr"]=[0.0092]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(2))
p["process"]=101
p["<Q>"]=numpy.sqrt(5.0)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=0.0050
p["uncorrErr"]=[0.0210]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
print("Making D2 by E154...")

DataCurrent=DataProcessor.DataSet.DataSet('E154_d2',"D2")
DataCurrent.comment="Taken from the text after (6)"
DataCurrent.reference="hep-ex/9705017"

DataCurrent.isNormalized=False
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

s_current=p["s"]=2*48.3*M_proton+M_proton**2

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(0))
p["process"]=101
p["<Q>"]=numpy.sqrt(3.6)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=-0.004
p["uncorrErr"]=[0.038,0.005]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
print("Making D2 by E155 (1999)...")

DataCurrent=DataProcessor.DataSet.DataSet('E155-1999_d2',"D2")
DataCurrent.comment="Taken from the text in page 5"
DataCurrent.reference="hep-ex/9901006"

DataCurrent.isNormalized=False
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

s_current=p["s"]=2*38.8*M_proton+M_proton**2

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(0))
p["process"]=100
p["<Q>"]=numpy.sqrt(5.0)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=0.005
p["uncorrErr"]=[0.008]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)


p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(1))
p["process"]=102
p["<Q>"]=numpy.sqrt(5.0)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=0.008
p["uncorrErr"]=[0.005]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
print("Making D2 by E155 (2002)...")

DataCurrent=DataProcessor.DataSet.DataSet('E155_d2',"D2")
DataCurrent.comment="Taken from the text after (5)"
DataCurrent.reference="hep-ex/0204028"

DataCurrent.isNormalized=False
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

s_current=p["s"]=2*32.3*M_proton+M_proton**2

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(0))
p["process"]=100
p["<Q>"]=numpy.sqrt(5.0)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=0.0025
p["uncorrErr"]=[0.0016,0.001]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)


p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(1))
p["process"]=102
p["<Q>"]=numpy.sqrt(5.0)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=0.0054
p["uncorrErr"]=[0.0023,0.0005]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
print("Making D2 by HallA (2004)...")

DataCurrent=DataProcessor.DataSet.DataSet('HallA-2004_d2',"D2")
DataCurrent.comment="Taken from (24)"
DataCurrent.reference="hep-ex/0405006"

DataCurrent.isNormalized=False
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

s_current=p["s"]=2*5.73*M_proton+M_proton**2

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(0))
p["process"]=101
p["<Q>"]=numpy.sqrt(5.0)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=0.0062
p["uncorrErr"]=[0.0028]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
print("Making D2 by HallA (2014)...")

DataCurrent=DataProcessor.DataSet.DataSet('HallA-2014_d2',"D2")
DataCurrent.comment="Taken from table I"
DataCurrent.reference="arxiv:1404.4003"

DataCurrent.isNormalized=False
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

s_current=p["s"]=2*5.73*M_proton+M_proton**2

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(0))
p["process"]=101
p["<Q>"]=numpy.sqrt(3.21)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=-0.00421
p["uncorrErr"]=[0.00079,0.00082]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(1))
p["process"]=101
p["<Q>"]=numpy.sqrt(4.32)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=-0.00035
p["uncorrErr"]=[0.00083,0.00069]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
print("Making D2 by HallA (2016)...")

DataCurrent=DataProcessor.DataSet.DataSet('HallA-2016_d2',"D2")
DataCurrent.comment="Taken from table X"
DataCurrent.reference="arxiv:1603.03612"

DataCurrent.isNormalized=False
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

s_current=p["s"]=2*5.73*M_proton+M_proton**2

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(0))
p["process"]=101
p["<Q>"]=numpy.sqrt(3.21)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=-0.00421
p["uncorrErr"]=[0.00079,0.00082,0.00008]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(1))
p["process"]=101
p["<Q>"]=numpy.sqrt(4.32)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=-0.00035
p["uncorrErr"]=[0.00083,0.00069,0.00007]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
print("Making D2 by HERMES...")

DataCurrent=DataProcessor.DataSet.DataSet('HERMES_d2',"D2")
DataCurrent.comment="Taken from text in the last page"
DataCurrent.reference="arxiv:1112.5584"

DataCurrent.isNormalized=False
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

s_current=p["s"]=2*27.6*M_proton+M_proton**2

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(0))
p["process"]=100
p["<Q>"]=numpy.sqrt(5.0)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=0.0148
p["uncorrErr"]=[0.0096,0.0048]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")


#%%
print("Making D2 by RSS (2006)..")

DataCurrent=DataProcessor.DataSet.DataSet('RSS-2006_d2',"D2")
DataCurrent.comment="Taken from text in the last page"
DataCurrent.reference="hep-exp/0608003"

DataCurrent.isNormalized=False
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

s_current=p["s"]=2*5.75*M_proton+M_proton**2

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(0))
p["process"]=100
p["<Q>"]=numpy.sqrt(1.3)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=0.0057
p["uncorrErr"]=[0.0009,0.0007]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
print("Making D2 by RSS (2010)..")

DataCurrent=DataProcessor.DataSet.DataSet('RSS-2008_d2',"D2")
DataCurrent.comment="Taken from table II"
DataCurrent.reference="arxiv:0812:0031"

DataCurrent.isNormalized=False
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

s_current=p["s"]=2*5.75*M_proton+M_proton**2

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(0))
p["process"]=100
p["<Q>"]=numpy.sqrt(1.28)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=0.0104
p["uncorrErr"]=[0.0004,0.0013]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(1))
p["process"]=102
p["<Q>"]=numpy.sqrt(1.28)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=0.0027
p["uncorrErr"]=[0.0008,0.0017]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)

p=DataProcessor.Point.CreateD2Point(DataCurrent.name+'.'+str(2))
p["process"]=101
p["<Q>"]=numpy.sqrt(1.28)
p["Q"]=[p["<Q>"]-0.5,p["<Q>"]+0.5]
p["s"]=s_current
p["xSec"]=-0.0075
p["uncorrErr"]=[0.0021,0.000]
p["corrErr"]=[]
p["thFactor"]=1.
DataCurrent.AddPoint(p)


print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%

print("Making SANE...")

DataCurrent=DataProcessor.DataSet.DataSet('SANE_d2',"D2")
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