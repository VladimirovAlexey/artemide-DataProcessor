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
import numpy

path_to_data=os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))+'/data/ATLAS/Zangular_1606_00689/'
path_to_save=ROOT_DIR+"/DataLib/DY_angular/"

M_Z=91.1876### mass of Z-boson 

#%%
### given in the text
#proc_current=[1,1,5]
s_current=2*800*0.938
Q1_current=[4.5,9.]
#Q1_current=[5.5,9.]
Q2_current=[10.7,15.]
y_current=[0.,0.8]
incCut=False
cutParam=[20.,20.,2.,4.5]

qTBINS=[0.,0.5,1.,1.5,2.,4.]

xFbins=[0.,0.25,0.35,0.45,0.55,0.8]

qTGRAND=[0.,4.]
#qTGRAND=[0.,2.]

Q_bins=[[4.5,5.5],[5.5,6.5],[6.5,7.5],[7.5,8.5],[8.5,9.0],[10.7,15.]]


## p+p from figure2 (by Sara)
xSec1=[0.019,0.051,0.039,0.046,0.194]
un1=[0.036, 0.026, 0.031, 0.046, 0.060]

## p+d from figure2 (by Sara)
xSec2=[-0.024,0.044,0.057,0.035,0.031]
un2=[0.024, 0.017, 0.020, 0.030, 0.041]

## p+p vs XF from fig3 (by Sara)
xSec3=[0.056,0.002,0.035,0.056,0.171]
un3=[0.031, 0.029, 0.031, 0.041, 0.047]


## p+d vs XF from fig3 (by Sara)
xSec4=[0.019,0.012,0.056,0.061,-0.023]
un4=[0.021, 0.021, 0.020, 0.025, 0.029]

## p+p vs Q from fig3 (by eye)
xSec5=[0.025,0.073,0.025,0.073,0.088,-0.055]
un5=[0.032,0.03,0.027,0.037,0.062,0.08]

## p+d vs Q from fig3 (by eye)
xSec6=[0.063,0.027,0.022,0.04,0.051,0.057]
un6=[0.023,0.02,0.027,0.023,0.04,0.05]

#%%
def FillData():
    
    DataCurrent.reference="arXiv:0811.4589"
    DataCurrent.isNormalized=False
        
    # read the data
    for n in range(5):                
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(n))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=[qTBINS[n],qTBINS[n+1]]
        p["Q"]=Q_current
        p["y"]=y_current
        p["thFactor"]=1.
        p["includeCuts"]=incCut
        p["cutParams"]=cutParam
        p["xSec"]=xSec[n]
        ###### The uncertanty is statistic + systematic. 
        ###### Systematic is lightly correlated I ignore it
        p["uncorrErr"].append(un[n])
        DataCurrent.AddPoint(p)    
    
    print("Done.  ")
    
#%%
def FillDataXF():
    
    DataCurrent.reference="arXiv:0811.4589"
    DataCurrent.isNormalized=False
        
    # read the data
    for n in range(5):                
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(n))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=qTGRAND
        p["Q"]=Q_current
        p["y"]=[xFbins[n],xFbins[n+1]]
        p["thFactor"]=1.
        p["includeCuts"]=incCut
        p["cutParams"]=cutParam
        p["xSec"]=xSec[n]
        ###### The uncertanty is statistic + systematic. 
        ###### Systematic is lightly correlated I ignore it
        p["uncorrErr"].append(un[n])
        DataCurrent.AddPoint(p)    
    
    print("Done.  ")
    
#%%
def FillDataQ():
    
    DataCurrent.reference="arXiv:0811.4589"
    DataCurrent.isNormalized=False
        
    # read the data
    for n in range(6):                
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(n))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=qTGRAND
        p["Q"]=Q_bins[n]
        p["y"]=y_current
        p["thFactor"]=1.
        p["includeCuts"]=incCut
        p["cutParams"]=cutParam
        p["xSec"]=xSec[n]
        ###### The uncertanty is statistic + systematic. 
        ###### Systematic is lightly correlated I ignore it
        p["uncorrErr"].append(un[n])
        DataCurrent.AddPoint(p)    
    
    print("Done.  ")

#%%
###############################################################################
DataCurrent=DataProcessor.DataSet.DataSet('E866(p+p)_Qlow_nu',"DY")
DataCurrent.comment="Angular coefficient nu from E866 (p+p)"

proc_current=[2,1,1,212]

Q_current=Q1_current
xSec=xSec1
un=un1

FillData()
DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

###############################################################################
DataCurrent=DataProcessor.DataSet.DataSet('E866(p+p)_Qhigh_nu',"DY")
DataCurrent.comment="Angular coefficient nu from E866 (p+p)"

Q_current=Q2_current

FillData()
DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

###############################################################################
DataCurrent=DataProcessor.DataSet.DataSet('E866(p+d)_Qlow_nu',"DY")
DataCurrent.comment="Angular coefficient nu from E866 (p+d)"

Q_current=Q1_current
xSec=xSec2
un=un2

FillData()
DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

###############################################################################
DataCurrent=DataProcessor.DataSet.DataSet('E866(p+d)_Qhigh_nu',"DY")
DataCurrent.comment="Angular coefficient nu from E866 (p+d)"

Q_current=Q2_current

FillData()
DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
#####################################################################
############ XF-DIFFERENTIAL  #######################################
#####################################################################
#####################################################################
###############################################################################
DataCurrent=DataProcessor.DataSet.DataSet('E866(p+p)_Qlow_dXF_nu',"DY")
DataCurrent.comment="Angular coefficient nu from E866 (p+p)"

proc_current=[2,1,1,212]

Q_current=Q1_current
xSec=xSec3
un=un3

FillDataXF()
DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

###############################################################################
DataCurrent=DataProcessor.DataSet.DataSet('E866(p+p)_Qhigh_dXF_nu',"DY")
DataCurrent.comment="Angular coefficient nu from E866 (p+p)"

Q_current=Q2_current

FillDataXF()
DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

###############################################################################
DataCurrent=DataProcessor.DataSet.DataSet('E866(p+d)_Qlow_dXF_nu',"DY")
DataCurrent.comment="Angular coefficient nu from E866 (p+p)"

Q_current=Q1_current
xSec=xSec4
un=un4

FillDataXF()
DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

###############################################################################
DataCurrent=DataProcessor.DataSet.DataSet('E866(p+d)_Qhigh_dXF_nu',"DY")
DataCurrent.comment="Angular coefficient nu from E866 (p+p)"

Q_current=Q2_current

FillDataXF()
DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
#####################################################################
############ Q-DIFFERENTIAL  #######################################
#####################################################################
#####################################################################

DataCurrent=DataProcessor.DataSet.DataSet('E866(p+p)_dQ_nu',"DY")
DataCurrent.comment="Angular coefficient nu from E866 (p+p)"

proc_current=[2,1,1,212]

xSec=xSec5
un=un5

FillDataQ()
DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

###############################################################################
DataCurrent=DataProcessor.DataSet.DataSet('E866(p+d)_dQ_nu',"DY")
DataCurrent.comment="Angular coefficient nu from E866 (p+p)"

xSec=xSec6
un=un6

FillDataQ()
DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

