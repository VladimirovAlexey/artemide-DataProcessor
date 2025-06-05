#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 20:15:46 2019
This program collect all the data on the SIDIS and save in  "SIDISdata_uncut.pkl"
@author: vla18041
"""
import sys
sys.path.append("/data/arTeMiDe_Repository/DataProcessor/")
import DataProcessor.Point
import DataProcessor.DataSet
import numpy


path_to_data="/data/arTeMiDe_Repository/data"
path_to_save="/data/arTeMiDe_Repository/DataProcessor/DataLib/wgt/"

totalData=[]

M_proton=0.938
m_pion=0.139
m_kaon=0.494

#%%
#### The data is taken from the table in the end of [1108.0489]
####[x, Q^2, z, pT, Api+, un1, un2, Api-, un1 ,un2 ]
table=[[0.156, 1.38, 0.50, 0.43,  0.02, 0.11, 0.03, 0.10, 0.07, 0.03],
       [0.206, 1.76, 0.52, 0.38,  0.04, 0.13, 0.06, 0.18, 0.11, 0.06],
       [0.265, 2.16, 0.54, 0.32, -0.13, 0.11, 0.05, 0.10, 0.07, 0.02],
       [0.349, 2.68, 0.58, 0.24, -0.27, 0.18, 0.13, 0.18, 0.10, 0.05]
       ]


    
#%%
############################################################################### PI+

### common variables
s_current=2*5.9*0.938+(0.938)**2
includeCuts=False
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

DataCurrent=DataProcessor.DataSet.DataSet('JLab6.ALT.pi+',"SIDIS")
DataCurrent.comment='JLab at 6 GeV data for A.LT. Pi+ off neutron'
DataCurrent.reference="1108.0489"
DataCurrent.normErr.append(0.028)#### delution factor according to polarization uncertanty


### Now we go line by line
num=0
for line in table:
    ## empty line
    
    ###### Adding point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(num))
    num+=1
    #print DataCurrent.name+'.'+str(i)
    p["process"]=[1,1,1,13003]
    p["weightProcess"]=[1,1,1,2003]
    p["s"]=s_current
    ### taken from the text
    p["Q"]=[numpy.sqrt(1.4),numpy.sqrt(2.7)]   
    p["x"]=[0.15+0.05*(num-1),0.15+0.05*num]
    p["z"]=[0.5,0.6]
    p["pT"]=[0.24,0.44] 
    
    p["<x>"]=line[0]
    p["<Q>"]=numpy.sqrt(line[1])
    p["<z>"]=line[2]
    p["<pT>"]=line[3]
    p["xSec"]=line[4]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    p["uncorrErr"].append(line[5])
    p["uncorrErr"].append(line[6])
    
    
    DataCurrent.AddPoint(p)
   
DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
############################################################################### PI-

### common variables
s_current=2*5.9*0.938+(0.938)**2
includeCuts=False
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

DataCurrent=DataProcessor.DataSet.DataSet('JLab6.ALT.pi-',"SIDIS")
DataCurrent.comment='JLab at 6 GeV data for A.LT. Pi- off neutron'
DataCurrent.reference="1108.0489"
DataCurrent.normErr.append(0.028)#### delution factor according to polarization uncertanty


### Now we go line by line
num=0
for line in table:
    ## empty line
    
    ###### Adding point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(num))
    num+=1
    #print DataCurrent.name+'.'+str(i)
    p["process"]=[1,1,-1,13003]
    p["weightProcess"]=[1,1,-1,2003]
    p["s"]=s_current
    ### taken from the text
    p["Q"]=[numpy.sqrt(1.4),numpy.sqrt(2.7)]   
    p["x"]=[0.15+0.05*(num-1),0.15+0.05*num]
    p["z"]=[0.5,0.6]
    p["pT"]=[0.24,0.44] 
    
    p["<x>"]=line[0]
    p["<Q>"]=numpy.sqrt(line[1])
    p["<z>"]=line[2]
    p["<pT>"]=line[3]
    p["xSec"]=line[7]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    p["uncorrErr"].append(line[8])
    p["uncorrErr"].append(line[9])
    
    
    DataCurrent.AddPoint(p)
   
DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")