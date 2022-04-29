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

path_to_data=os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))+'/data/'
path_to_save=ROOT_DIR+"/DataLib/unpolDY/"

M_Z=91.1876### mass of Z-boson

#%%
### given in the text
proc_current=[1,1,5]
s_current=13000.**2
incCut=True
cutParam=[25.,20.,-2.4,2.4]
lumUncertainty=0.012
y_current=[-2.4,2.4]

#%%

###############################################################################
###########################CMS 13 Q1- Q2 ########################################
def ParseCMS(Q1,Q2):
    print("Read CMS13 "+str(Q1)+"to"+str(Q2)+" file...")
    f = open(path_to_data+"CMS/CMS13_"+str(Q1)+"to"+str(Q2)+".dat")
    
    data_from_f=[]
    
    for line in f:    
        data_from_f.append(line.rstrip('\n'))
    
    f.close()
    
    print("Done.  =>     Convert to numbers ...")
    
    Q_current=[float(Q1),float(Q2)]
    
    
    data_here=data_from_f[10:]
    
    for i in range(len(data_here)):
        data_here[i]=data_here[i].split(",")
        data_here[i]=[float(j) for j in data_here[i]]
    
    print("Done.  =>     Create points & append to data set ...")
    DataCurrent=DataProcessor.DataSet.DataSet('CMS13_dQ_'+str(Q1)+"to"+str(Q2),"DY")
    DataCurrent.comment="CMS 13TeV 2021 preliminary"
    DataCurrent.reference="CMS-PAS-SMP-20-003"
    
    DataCurrent.isNormalized=False
    
    DataCurrent.normErr.append(lumUncertainty)
    
    for i in range(len(data_here)):
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=data_here[i][0:2]
        p["Q"]=Q_current
        p["y"]=y_current
        p["thFactor"]=1/(p["qT"][1]-p["qT"][0])#devide by bin size
        p["includeCuts"]=incCut
        p["cutParams"]=cutParam
        p["xSec"]=data_here[i][2]
        ###### The uncertanty is statistic + systematic. 
        ###### They added lumi.uncertanty to it, I subtract (in quadrature)
        ###### Rest Systematic is lightly correlated I ignore it
        ###### Uncertnaty is given in percentage
        p["uncorrErr"].append(p["xSec"]*numpy.sqrt(data_here[i][3]**2-1.2**2)*0.01)
        DataCurrent.AddPoint(p)    
    
    print("Done.  ")
    
    DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
ParseCMS(50, 76)
ParseCMS(76, 106)
ParseCMS(106, 170)
ParseCMS(170, 350)
ParseCMS(350, 1000)