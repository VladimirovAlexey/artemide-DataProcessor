#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 20:15:46 2019

Program that parse various DY data files to ADP-frendly formal

@author: vla18041
"""
import os
import sys
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"

#%%
sys.path.append(ROOT_DIR)

import DataProcessor.Point
import DataProcessor.DataSet
import numpy

path_to_data=os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))+"/data/"
path_to_save=ROOT_DIR+"/DataLib/unpolDY/"

M_Z=91.### mass of Z-boson

#%%
###############################################################################
###########################E288 200############################################
print("Read E288(200) file ...")
f = open(path_to_data+"E288/E288_200.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:9]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]


for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")    
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('E228-200',"DY")
DataCurrent.comment="E288 (200) data"
DataCurrent.reference="Phys.Rev.D 23 (1981) 604"

DataCurrent.isNormalized=False
proc_current=[2,1,1,101]
s_current=19.42**2
y_current=[0.1,0.7]
lumUncertainty=0.25
DataCurrent.normErr.append(lumUncertainty)

for j in range(7):
    Q_current=[float(4+j),float(5+j)]
    print('for Q = ',Q_current)
    for i in range(len(data_from_f)):
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(int(Q_current[0]))+'Q'+str(int(Q_current[1]))+'.'+str(i))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
        # 0.31.. is for 1/pi from invariant cross-secX, 0.001 for nb
        p["thFactor"]=1/(p["qT"][1]**2-p["qT"][0]**2)/(y_current[1]-y_current[0])*0.3183098861838*1000
        p["Q"]=Q_current
        p["y"]=y_current
        p["xSec"]=data_from_f[i][3+3*j]
        p["includeCuts"]=False
        p["uncorrErr"].append((data_from_f[i][4+3*j]-data_from_f[i][5+3*j])/2.)
        if p["xSec"] != -50:
            DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

###############################################################################
###########################E288 300############################################
print("Read E288(300) file ...")
f = open(path_to_data+"E288/E288_300.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:9]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]


for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")    
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('E228-300',"DY")
DataCurrent.comment="E288 (300) data"
DataCurrent.reference="Phys.Rev.D 23 (1981) 604"

DataCurrent.isNormalized=False
proc_current=[2,1,1,101]
s_current=23.73**2
y_current=[0.21-0.3,0.21+0.3]
lumUncertainty=0.25
DataCurrent.normErr.append(lumUncertainty)


for j in range(8):
    Q_current=[float(4+j),float(5+j)]
    print('for Q = ',Q_current)
    for i in range(len(data_from_f)):
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(int(Q_current[0]))+'Q'+str(int(Q_current[1]))+'.'+str(i))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
        # 0.31.. is for 1/pi from invariant cross-secX, 0.001 for nb
        p["thFactor"]=1/(p["qT"][1]**2-p["qT"][0]**2)/(y_current[1]-y_current[0])*0.3183098861838*1000
        p["Q"]=Q_current
        p["y"]=y_current
        p["xSec"]=data_from_f[i][3+3*j]
        p["includeCuts"]=False
        p["uncorrErr"].append((data_from_f[i][4+3*j]-data_from_f[i][5+3*j])/2.)
        if p["xSec"] != -50:
            DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")


###############################################################################
###########################E288 400############################################
print("Read E288(400) file ...")
f = open(path_to_data+"E288/E288_400.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:9]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]


for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")    
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('E228-400',"DY")
DataCurrent.comment="E288 (400) data"
DataCurrent.reference="Phys.Rev.D 23 (1981) 604"

DataCurrent.isNormalized=False
proc_current=[2,1,1,101]
s_current=27.43**2
y_current=[0.03-0.3,0.03+0.3]
lumUncertainty=0.25
DataCurrent.normErr.append(lumUncertainty)


for j in range(9):
    Q_current=[float(5+j),float(6+j)]
    print('for Q = ',Q_current)
    for i in range(len(data_from_f)):
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(int(Q_current[0]))+'Q'+str(int(Q_current[1]))+'.'+str(i))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
        # 0.31.. is for 1/pi from invariant cross-secX, 0.001 for nb
        p["thFactor"]=1/(p["qT"][1]**2-p["qT"][0]**2)/(y_current[1]-y_current[0])*0.3183098861838*1000
        p["Q"]=Q_current
        p["y"]=y_current
        p["includeCuts"]=False
        p["xSec"]=data_from_f[i][3+3*j]
        p["uncorrErr"].append((data_from_f[i][4+3*j]-data_from_f[i][5+3*j])/2.)
        if p["xSec"] != -50:
            DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")


#%%
###############################################################################
###########################E772################################################
print("Read E772  4-9 file ...")
f = open(path_to_data+"E772/E772_4to9.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:8]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]


for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")    
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('E772',"DY")
DataCurrent.comment="E772 data"
DataCurrent.reference="Phys.Rev.D 50 (1994) 3-38 + Erratum D60 (1999) 119903"

DataCurrent.isNormalized=False
proc_current=[2,1,1,102]
s_current=38.76**2
y_current=[0.1,0.3]
lumUncertainty=0.10
DataCurrent.normErr.append(lumUncertainty)

for j in range(4):
    Q_current=[float(5+j),float(6+j)]
    print('for Q = ',Q_current)
    for i in range(len(data_from_f)):
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(int(Q_current[0]))+'Q'+str(int(Q_current[1]))+'.'+str(i))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=[data_from_f[i][0]-0.125,data_from_f[i][0]+0.125]
        # 0.31.. is for 1/pi from invariant cross-secX
        p["thFactor"]=1/(p["qT"][1]**2-p["qT"][0]**2)/(y_current[1]-y_current[0])*0.3183098861838
        p["Q"]=Q_current
        p["y"]=y_current
        p["xSec"]=data_from_f[i][3+3*j]
        p["includeCuts"]=False
        p["uncorrErr"].append((data_from_f[i][4+3*j]-data_from_f[i][5+3*j])/2.)
        if p["xSec"] != -50:
            DataCurrent.AddPoint(p)

print("Read E772 11-15 file ...")
f = open(path_to_data+"E772/E772_11to15.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:8]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]


for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")    
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
for j in range(4):
    Q_current=[float(11+j),float(12+j)]
    print('for Q = ',Q_current)
    for i in range(len(data_from_f)):
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(int(Q_current[0]))+'Q'+str(int(Q_current[1]))+'.'+str(i))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=[data_from_f[i][0]-0.125,data_from_f[i][0]+0.125]
        # 0.31.. is for 1/pi from invariant cross-secX
        p["thFactor"]=1/(p["qT"][1]**2-p["qT"][0]**2)/(y_current[1]-y_current[0])*0.3183098861838
        p["Q"]=Q_current
        p["y"]=y_current
        p["xSec"]=data_from_f[i][3+3*j]
        p["includeCuts"]=False
        p["uncorrErr"].append((data_from_f[i][4+3*j]-data_from_f[i][5+3*j])/2.)        
        if p["xSec"] != -50:
            DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")


#%%
###############################################################################
###########################E605################################################
print("Read E605  7-8 file ...")
f = open(path_to_data+"E605/E605_78.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:6]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
del data_from_f[-1]
del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]


for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")    
    #del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('E605',"DY")
DataCurrent.comment="E605 data"
DataCurrent.reference="Phys.Rev.D 43 (1991) 2815"

DataCurrent.isNormalized=False
proc_current=[2,1,1,102]
s_current=38.76**2
y_current=[-0.1,0.2]
lumUncertainty=0.15
sysError=0.05
DataCurrent.normErr.append(lumUncertainty)

Q_current=[7.,8.]
print('for Q = ',Q_current)
for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.7Q8.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][0]-0.1,data_from_f[i][0]+0.1]
    # 0.31.. is for 1/pi from invariant cross-secX
    p["thFactor"]=1/(p["qT"][1]**2-p["qT"][0]**2)/(y_current[1]-y_current[0])*0.3183098861838
    p["Q"]=Q_current
    p["y"]=y_current
    p["xSec"]=data_from_f[i][3]
    p["includeCuts"]=False
    p["uncorrErr"].append((data_from_f[i][4]+data_from_f[i][5])/2.)
    p["uncorrErr"].append(sysError*p["xSec"])
    if p["xSec"] != -50:
        DataCurrent.AddPoint(p)

print("Read E605 8-9 file ...")
f = open(path_to_data+"E605/E605_89.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:6]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
del data_from_f[-1]
del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]


for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")    
    #del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")

Q_current=[8.,9.]
print('for Q = ',Q_current)
for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.8Q9.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][0]-0.1,data_from_f[i][0]+0.1]
    # 0.31.. is for 1/pi from invariant cross-secX
    p["thFactor"]=1/(p["qT"][1]**2-p["qT"][0]**2)/(y_current[1]-y_current[0])*0.3183098861838
    p["Q"]=Q_current
    p["y"]=y_current
    p["includeCuts"]=False
    p["xSec"]=data_from_f[i][3]
    p["uncorrErr"].append((data_from_f[i][4]+data_from_f[i][5])/2.)
    p["uncorrErr"].append(sysError*p["xSec"])
    if p["xSec"] != -50:
        DataCurrent.AddPoint(p)


print("Read E605 10-11 file ...")
f = open(path_to_data+"E605/E605_1011.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:6]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
del data_from_f[-1]
del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]


for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")    
    #del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")

Q_current=[10.5,11.5]
print('for Q = ',Q_current)
for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.10Q11.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][0]-0.1,data_from_f[i][0]+0.1]
    # 0.31.. is for 1/pi from invariant cross-secX
    p["thFactor"]=1/(p["qT"][1]**2-p["qT"][0]**2)/(y_current[1]-y_current[0])*0.3183098861838
    p["Q"]=Q_current
    p["y"]=y_current
    p["includeCuts"]=False
    p["xSec"]=data_from_f[i][3]
    p["uncorrErr"].append((data_from_f[i][4]+data_from_f[i][5])/2.)
    p["uncorrErr"].append(sysError*p["xSec"])
    if p["xSec"] != -50:
        DataCurrent.AddPoint(p)

print("Read E605 11-13 file ...")
f = open(path_to_data+"E605/E605_1113.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:6]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
del data_from_f[-1]
del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]


for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")    
    #del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")

Q_current=[11.5,13.5]
print('for Q = ',Q_current)
for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.11Q13.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][0]-0.1,data_from_f[i][0]+0.1]
    # 0.31.. is for 1/pi from invariant cross-secX
    p["thFactor"]=1/(p["qT"][1]**2-p["qT"][0]**2)/(y_current[1]-y_current[0])*0.3183098861838
    p["Q"]=Q_current
    p["y"]=y_current
    p["includeCuts"]=False
    p["xSec"]=data_from_f[i][3]
    p["uncorrErr"].append((data_from_f[i][4]+data_from_f[i][5])/2.)
    p["uncorrErr"].append(sysError*p["xSec"])
    if p["xSec"] != -50:
        DataCurrent.AddPoint(p)

print("Done.")

print("Read E605 13-18 file ...")
f = open(path_to_data+"E605/E605_1318.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:6]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
del data_from_f[-1]
del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]


for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")    
    #del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")

Q_current=[13.5,18.0]
print('for Q = ',Q_current)
for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.13Q18.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][0]-0.1,data_from_f[i][0]+0.1]
    # 0.31.. is for 1/pi from invariant cross-secX
    p["thFactor"]=1/(p["qT"][1]**2-p["qT"][0]**2)/(y_current[1]-y_current[0])*0.3183098861838
    p["Q"]=Q_current
    p["y"]=y_current
    p["includeCuts"]=False
    p["xSec"]=data_from_f[i][3]
    p["uncorrErr"].append((data_from_f[i][4]+data_from_f[i][5])/2.)
    p["uncorrErr"].append(sysError*p["xSec"])
    if p["xSec"] != -50:
        DataCurrent.AddPoint(p)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")


#%%

print("Read E537 dQ file ...")
f = open(path_to_data+"FNAL-537/pbar+W(dQ).dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")
del data_from_f[0:9]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]

print("Done.  =>     Create points & append to data set ...")

DataCurrent=DataProcessor.DataSet.DataSet('E537-dQ',"DY")
DataCurrent.comment="E537 data Q-differential"
DataCurrent.reference="Phys.Rev.D 93 (1988) 1377"

DataCurrent.isNormalized=False
proc_current=[2,-1,1,103]
s_current=235.4
y_current=[-0.1,1.0]
sysError=0.08
DataCurrent.normErr.append(0.08)

Q_name='?'
k=0
for ss in data_from_f:
    if ss[0:4] == '#: M':  
        # in this case we update the current value of Q
        endofss=ss[-8:]
        Q_name=endofss.replace('.','')
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
    

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%

print("Read E537 dxF file ...")
f = open(path_to_data+"FNAL-537/pbar+W(dxF).dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")
del data_from_f[0:9]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]

print("Done.  =>     Create points & append to data set ...")

DataCurrent=DataProcessor.DataSet.DataSet('E537-dxF',"DY")
DataCurrent.comment="E537 data xF-differential"
DataCurrent.reference="Phys.Rev.D 93 (1988) 1377"

DataCurrent.isNormalized=False
proc_current=[2,-1,1,103]
s_current=235.4

Q_current=[4.0,9.0]
sysError=0.08
DataCurrent.normErr.append(0.08)

y_name='?'
k=0
for ss in data_from_f:
    if ss[0:5] == '#: XL':  
        # in this case we update the current value of y
        endofss=ss[-8:]
        y_name=endofss.replace('.','')
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
    

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")