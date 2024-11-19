#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 20:15:46 2019

Program that parse various DY data files to ADP-frendly formal

@author: vla18041
"""

import os
import DataProcessor.Point
import DataProcessor.DataSet
import numpy

ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"
path_to_data=os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))+"/data/"
path_to_savePI=ROOT_DIR+"/DataLib/piDY/"


M_Z=91.### mass of Z-boson

#%%

print("Read E537-pi dQ file ...")
f = open(path_to_data+"FNAL-537/piminus+W(dQ).dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Create points & append to data set ...")

DataCurrent=DataProcessor.DataSet.DataSet('E537(pi)-dQ',"DY")
DataCurrent.comment="E537 pi-data Q-differential"
DataCurrent.reference="Phys.Rev.D 93 (1988) 1377"

DataCurrent.isNormalized=False
proc_current=[2,2,1,103]
s_current=235.4

y_current=[-0.1,1.0]
if DataCurrent.isNormalized:
    sysError=0.
else:
    sysError=0.08
    
DataCurrent.normErr.append(sysError)

Q_name='?'
k=0
for ss in data_from_f:
    if ss[0:4] == '#: M':  
        # in this case we update the current value of Q
        endofss=ss[-8:]
        Q_name=endofss.replace('.','').replace(",","-")
        endofss=endofss.split('TO')
        Q_current=[float(i) for i in endofss]
        k=0
    elif ss[0:2] == '#:' or ss=='' or ss[0:3]=='"PT':
        #in this case we skip it
        pass
    else:
        #this is data string [pt-cetner, pt_low,pt-high, xSec , err+, err-]
        data_current=ss.split(",")
        data_current=[float(j) for j in data_current]
        #------------------------ADD POINT------------
        
        p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+Q_name+'.'+str(k))
        k+=1
        p["process"]=proc_current
        p["s"]=s_current
        # they publish pt^2
        p["qT"]=[numpy.sqrt(data_current[1]),numpy.sqrt(data_current[2])]
        #factor 0.001 is due to pb->nb
        p["thFactor"]=1/(p["qT"][1]**2-p["qT"][0]**2)/(Q_current[1]-Q_current[0])*0.001
        p["Q"]=Q_current
        p["y"]=y_current
        p["includeCuts"]=False
        p["xSec"]=data_current[3]
        p["uncorrErr"].append((data_current[4]-data_current[5])/2.)
        if(p["uncorrErr"][0]>0):
            DataCurrent.AddPoint(p)
    

DataCurrent.SaveToCSV(path_to_savePI+DataCurrent.name+".csv")

#%%

print("Read E537(pi) dxF file ...")
f = open(path_to_data+"FNAL-537/pminus+W(dxF).dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")


print("Done.  =>     Create points & append to data set ...")

DataCurrent=DataProcessor.DataSet.DataSet('E537(pi)-dxF',"DY")
DataCurrent.comment="E537(pi) data xF-differential"
DataCurrent.reference="Phys.Rev.D 93 (1988) 1377"

DataCurrent.isNormalized=False
proc_current=[2,2,1,103]
s_current=235.4

Q_current=[4.0,9.0]
if DataCurrent.isNormalized:
    sysError=0.
else:
    sysError=0.08
    
DataCurrent.normErr.append(sysError)

y_name='?'
k=0
for ss in data_from_f:
    if ss[0:5] == '#: XL':  
        # in this case we update the current value of y
        endofss=ss[-8:]
        y_name=endofss.replace('.','').replace(",","-")
        endofss=endofss.split('TO')
        y_current=[float(i) for i in endofss]
        k=0
    elif ss[0:2] == '#:' or ss=='' or ss[0:3]=='"PT':
        #in this case we skip it
        pass
    else:
        #this is data string [pt-cetner, pt_low,pt-high, xSec , err+, err-]
        data_current=ss.split(",")
        data_current=[float(j) for j in data_current]
        #------------------------ADD POINT------------
        
        p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+y_name+'.'+str(k))
        k+=1
        p["process"]=proc_current
        p["s"]=s_current
        # they publish pt^2
        p["qT"]=[numpy.sqrt(data_current[1]),numpy.sqrt(data_current[2])]
        #factor 0.001 is due to pb->nb
        p["thFactor"]=1/(p["qT"][1]**2-p["qT"][0]**2)/(y_current[1]-y_current[0])*0.001
        p["Q"]=Q_current
        p["y"]=y_current
        p["includeCuts"]=False
        p["xSec"]=data_current[3]
        p["uncorrErr"].append((data_current[4]-data_current[5])/2.)
        if(p["uncorrErr"][0]>0):
            DataCurrent.AddPoint(p)
    

DataCurrent.SaveToCSV(path_to_savePI+DataCurrent.name+".csv")

#%%

print("Read E615-pi dQ file ...")
f = open(path_to_data+"FNAL-615/FNAL-615(dQ).dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Create points & append to data set ...")

DataCurrent=DataProcessor.DataSet.DataSet('E615(pi)-dQ',"DY")
DataCurrent.comment="E615 pi-data Q-differential"
DataCurrent.reference="Phys.Rev.D 39 (1989) 92-122"

DataCurrent.isNormalized=False
proc_current=[2,2,1,103]
s_current=473.6

y_current=[0.0,1.0]
if DataCurrent.isNormalized:
    sysError=0.
else:
    sysError=0.16

DataCurrent.normErr.append(sysError)

Q_name='?'
k=0
for ss in data_from_f:
    if ss[0:4] == '#Q= ':  
        # in this case we update the current value of Q
        endofss=ss[4:]
        Q_name=endofss.replace(",","-")
        endofss=endofss.replace('[','')
        endofss=endofss.replace(']','')
        endofss=endofss.split(',')
        Q_current=[float(i) for i in endofss]
        k=0
    elif ss[0:1] == '#' or ss=='':
        #in this case we skip it
        pass
    else:
        #this is data string [pt-cetner, pt_low,pt-high, xSec , err+, err-]
        data_current=ss.split("\t")
        data_current=[float(j) for j in data_current]
        #------------------------ADD POINT------------
        
        p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+Q_name+'.'+str(k))
        k+=1
        p["process"]=proc_current
        p["s"]=s_current
        # they publish pt^2
        p["qT"]=[data_current[0]-0.125,data_current[0]+0.125]
        #factor 0.001 is due to pb->nb
        p["Q"]=Q_current
        p["y"]=y_current
        p["thFactor"]=1/(p["qT"][1]-p["qT"][0])/(Q_current[1]-Q_current[0])*0.001        
        p["xSec"]=data_current[1]
        p["includeCuts"]=False
        p["uncorrErr"].append(data_current[2])
        if(p["uncorrErr"][0]>0):
            DataCurrent.AddPoint(p)
    

DataCurrent.SaveToCSV(path_to_savePI+DataCurrent.name+".csv")

#%% 
##### We also save separate Q-bins 
for QQ in [4.05,4.50,4.95,5.40,5.85,6.75,7.65,9.00,10.35,11.70]:
    dNew=DataProcessor.DataSet.DataSet(DataCurrent.name+'-'+str(QQ),"DY")
    dNew.comment=DataCurrent.comment+' only Q='+str(QQ)
    dNew.reference=DataCurrent.reference
    dNew.normalized=DataCurrent.isNormalized
    dNew.normErr=DataCurrent.normErr
    for p in DataCurrent.points:
        if(p["Q"][0]==QQ):
            dNew.AddPoint(p)
    dNew.SaveToCSV(path_to_savePI+dNew.name+".csv")

#%%

print("Read E615-pi dxF file ...")
f = open(path_to_data+"FNAL-615/FNAL-615(dx).dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Create points & append to data set ...")

DataCurrent=DataProcessor.DataSet.DataSet('E615(pi)-dxF',"DY")
DataCurrent.comment="E615 pi-data x-differential"
DataCurrent.reference="Phys.Rev.D 39 (1989) 92-122"

DataCurrent.isNormalized=False
proc_current=[2,2,1,103]
s_current=473.6

Q_current=[4.05,8.55]
if DataCurrent.isNormalized:
    sysError=0.
else:
    sysError=0.16

DataCurrent.normErr.append(sysError)

y_name='?'
k=0
for ss in data_from_f:
    if ss[0:4] == '#X= ':  
        # in this case we update the current value of Q
        endofss=ss[4:]
        y_name=endofss.replace(",","-")
        endofss=endofss.replace('[','')
        endofss=endofss.replace(']','')
        endofss=endofss.split(',')
        y_current=[float(i) for i in endofss]
        k=0
    elif ss[0:1] == '#' or ss=='':
        #in this case we skip it
        pass
    else:
        #this is data string [pt-cetner, pt_low,pt-high, xSec , err+, err-]
        data_current=ss.split("\t")
        data_current=[float(j) for j in data_current]
        #------------------------ADD POINT------------
        
        p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+y_name+'.'+str(k))
        k+=1
        p["process"]=proc_current
        p["s"]=s_current
        # they publish pt^2
        p["qT"]=[data_current[0]-0.125,data_current[0]+0.125]
        #factor 0.001 is due to pb->nb
        p["Q"]=Q_current
        p["y"]=y_current
        p["thFactor"]=1/(p["qT"][1]-p["qT"][0])/(y_current[1]-y_current[0])*0.001        
        p["includeCuts"]=False
        p["xSec"]=data_current[1]
        p["uncorrErr"].append(data_current[2])
        if(p["uncorrErr"][0]>0):
            DataCurrent.AddPoint(p)
    

DataCurrent.SaveToCSV(path_to_savePI+DataCurrent.name+".csv")

#%% 
##### We also save separate y-bins 
for yy in [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
    dNew=DataProcessor.DataSet.DataSet(DataCurrent.name+'-'+str(yy),"DY")
    dNew.comment=DataCurrent.comment+'[y='+str(yy)+']'
    dNew.reference=DataCurrent.reference
    dNew.normalized=DataCurrent.isNormalized
    dNew.normErr=DataCurrent.normErr
    for p in DataCurrent.points:
        if(p["y"][0]==yy):
            dNew.AddPoint(p)
    dNew.SaveToCSV(path_to_savePI+dNew.name+".csv")
