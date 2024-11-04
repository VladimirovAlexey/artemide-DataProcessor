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
proc_current=[1,1,1,3]
s_current=13000.**2
incCut=True
cutParam=[25.,20.,-2.4,2.4]
lumUncertainty=0.012
y_current=[-2.4,2.4]

#%%
#### FSR factors computed by Louis Moureaux 
#### I extract them by taking the ratio of nofsr/fsr
fsr={"50to76":
     [1.00098537, 0.97497963, 0.92974446, 0.87653573, 0.80923609, 0.73646756,
      0.66225815, 0.5801822,  0.51508677, 0.51669914, 0.59243914, 0.68497092,
      0.72740284, 0.84009018, 0.82205728, 0.9087    ],
     "76to106":
    [1.05285829, 1.05091075, 1.04498482, 1.03978998, 1.03490338, 1.02954884,
     1.02554925, 1.02210408, 1.02117097, 1.01910361, 1.01701754, 1.01836032,
     1.01750745, 1.01730813, 1.01867058, 1.02046079, 1.02091208, 1.02347862,
     1.02336165, 1.02557772, 1.0269174,  1.02751717, 1.0293184,  1.03087537,
     1.03446688, 1.0432738,  0.9932898,  1.025378,   1.02528262, 1.02589994,
     1.02191976, 1.02025379, 1.02884329, 1.02806738, 1.0233248,  1.01200851,
     1.02687311],
    "106to170":
    [1.07275747, 1.066437,   1.06056855, 1.05376852, 1.05392052, 1.04900427,
     1.04836114, 1.04849858, 1.04644792, 1.04296663, 1.04972161, 1.0434394,
     1.04831678, 1.68206705, 1.10846788, 1.04015155],
    "170to350":
    [1.08153606, 1.06157922, 1.06100099, 1.06240545, 1.039029,   1.04436067,
     1.03852346, 1.02823426, 1.02717946, 1.0260605,  1.0184391,  1.00729932,
     1.06753421],
    "350to1000":
    [1.07715219, 1.08831779, 1.03406113, 1.05885152, 1.09046837, 1.02711162,
     1.04009832, 1.03338694, 1.04638257, 1.01053939, 1.01408975, 1.00578608,
     1.04662236]
    }

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
        #devide by bin size and by fsr
        p["thFactor"]=1/(p["qT"][1]-p["qT"][0])/fsr[str(Q1)+"to"+str(Q2)][i]
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