#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 20:15:46 2019

Program that parse various DY data files to ADP-frendly formal

@author: vla18041
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/DataProcessor/")

import DataProcessor.Point
import DataProcessor.DataSet

path_to_data="/data/arTeMiDe_Repository/data/"
path_to_save="/data/arTeMiDe_Repository/DataProcessor/DataLib/unpolDY/"

M_Z=91.### mass of Z-boson

#%%
###############################################################################
###########################ATLAS 8 00y04########################################
print("Read ATLAS 8 0.0<y<0.4 file...")
f = open(path_to_data+"ATLAS/ATLAS_8TeV_66to116_00y04.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:7]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")

#print data_from_f

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('A8-00y04-norm',"DY")
DataCurrent.comment="ATLAS 8TeV 0.0<|y|<0.4 normalized to 1/sigma"
DataCurrent.reference="arXiv:1512.02192"

DataCurrent.isNormalized=True
proc_current=[1,1,1,3]
s_current=8000.**2
Q_current=[66.,116.]
y_current=[0.,0.4]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.0
##DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
    p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_from_f[i][38]
    p["uncorrErr"].append((data_from_f[i][39]-data_from_f[i][40])/2.)
    p["uncorrErr"].append((data_from_f[i][41]-data_from_f[i][42])/2.)
    p["corrErr"].append((data_from_f[i][43]-data_from_f[i][44])/2.)
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################ATLAS 8 04y08########################################
print("Read ATLAS 8 0.4<y<0.8 file...")
f = open(path_to_data+"ATLAS/ATLAS_8TeV_66to116_04y08.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:7]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")

#print data_from_f

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('A8-04y08-norm',"DY")
DataCurrent.comment="ATLAS 8TeV 0.4<|y|<0.8 normalized to 1/sigma"
DataCurrent.reference="arXiv:1512.02192"

DataCurrent.isNormalized=True
proc_current=[1,1,1,3]
s_current=8000.**2
Q_current=[66.,116.]
y_current=[0.4,0.8]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.0
#DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
    p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_from_f[i][38]
    p["uncorrErr"].append((data_from_f[i][39]-data_from_f[i][40])/2.)
    p["uncorrErr"].append((data_from_f[i][41]-data_from_f[i][42])/2.)
    p["corrErr"].append((data_from_f[i][43]-data_from_f[i][44])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################ATLAS 8 08y12########################################
print("Read ATLAS 8 0.8<y<1.2 file...")
f = open(path_to_data+"ATLAS/ATLAS_8TeV_66to116_08y12.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:7]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")

#print data_from_f

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('A8-08y12-norm',"DY")
DataCurrent.comment="ATLAS 8TeV 0.8<|y|<1.2 normalized to 1/sigma"
DataCurrent.reference="arXiv:1512.02192"

DataCurrent.isNormalized=True
proc_current=[1,1,1,3]
s_current=8000.**2
Q_current=[66.,116.]
y_current=[0.8,1.2]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.0
#DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
    p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y    
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_from_f[i][38]
    p["uncorrErr"].append((data_from_f[i][39]-data_from_f[i][40])/2.)
    p["uncorrErr"].append((data_from_f[i][41]-data_from_f[i][42])/2.)
    p["corrErr"].append((data_from_f[i][43]-data_from_f[i][44])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################ATLAS 8 12y16#######################################
print("Read ATLAS 8 1.2<y<1.6 file...")
f = open(path_to_data+"ATLAS/ATLAS_8TeV_66to116_12y16.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:7]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")

#print data_from_f

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('A8-12y16-norm',"DY")
DataCurrent.comment="ATLAS 8TeV 1.2<|y|<1.6 normalized to 1/sigma"
DataCurrent.reference="arXiv:1512.02192"

DataCurrent.isNormalized=True
proc_current=[1,1,1,3]
s_current=8000.**2
Q_current=[66.,116.]
y_current=[1.2,1.6]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.0
#DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
    p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_from_f[i][38]
    p["uncorrErr"].append((data_from_f[i][39]-data_from_f[i][40])/2.)
    p["uncorrErr"].append((data_from_f[i][41]-data_from_f[i][42])/2.)
    p["corrErr"].append((data_from_f[i][43]-data_from_f[i][44])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################ATLAS 8 16y20########################################
print("Read ATLAS 8 1.6<y<2.0 file...")
f = open(path_to_data+"ATLAS/ATLAS_8TeV_66to116_16y20.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:7]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")

#print data_from_f

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('A8-16y20-norm',"DY")
DataCurrent.comment="ATLAS 8TeV 1.6<|y|<2.0 normalized to 1/sigma"
DataCurrent.reference="arXiv:1512.02192"

DataCurrent.isNormalized=True
proc_current=[1,1,1,3]
s_current=8000.**2
Q_current=[66.,116.]
y_current=[1.6,2.0]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.0
#DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
    p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_from_f[i][38]
    p["uncorrErr"].append((data_from_f[i][39]-data_from_f[i][40])/2.)
    p["uncorrErr"].append((data_from_f[i][41]-data_from_f[i][42])/2.)
    p["corrErr"].append((data_from_f[i][43]-data_from_f[i][44])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################ATLAS 8 20y24#######################################
print("Read ATLAS 8 2.0<y<2.4 file...")
f = open(path_to_data+"ATLAS/ATLAS_8TeV_66to116_20y24.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:7]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")

#print data_from_f

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('A8-20y24-norm',"DY")
DataCurrent.comment="ATLAS 8TeV 2.0<|y|<2.4 normalized to 1/sigma"
DataCurrent.reference="arXiv:1512.02192"

DataCurrent.isNormalized=True
proc_current=[1,1,1,3]
s_current=8000.**2
Q_current=[66.,116.]
y_current=[2.0,2.4]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.0
#DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
    p["thFactor"]=2/(p["qT"][1]-p["qT"][0])#devide by bin size## since we have to symmetrize y
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_from_f[i][38]
    p["uncorrErr"].append((data_from_f[i][39]-data_from_f[i][40])/2.)
    p["uncorrErr"].append((data_from_f[i][41]-data_from_f[i][42])/2.)
    p["corrErr"].append((data_from_f[i][43]-data_from_f[i][44])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")
DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################ATLAS 8 46-66#######################################
print("Read ATLAS 8 46-66 file...")
f = open(path_to_data+"ATLAS/ATLAS_8TeV_46to66.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:7]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")

#print data_from_f

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('A8-46Q66-norm',"DY")
DataCurrent.comment="ATLAS 8TeV 46<Q<66 normalized to 1/sigma"
DataCurrent.reference="arXiv:1512.02192"

DataCurrent.isNormalized=True
proc_current=[1,1,1,3]
s_current=8000.**2
Q_current=[46.,66.]
y_current=[-2.4,2.4]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.0
#DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
    p["thFactor"]=1/(p["qT"][1]-p["qT"][0])#devide by bin size
    p["Q"]=Q_current
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_from_f[i][38]
    p["uncorrErr"].append((data_from_f[i][39]-data_from_f[i][40])/2.)
    p["uncorrErr"].append((data_from_f[i][41]-data_from_f[i][42])/2.)
    p["corrErr"].append((data_from_f[i][43]-data_from_f[i][44])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################ATLAS 8 116-150#######################################
print("Read ATLAS 8 116-150 file...")
f = open(path_to_data+"ATLAS/ATLAS_8TeV_116to150.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:7]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")

#print data_from_f

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('A8-116Q150-norm',"DY")
DataCurrent.comment="ATLAS 8TeV 116<Q<150 normalized to 1/sigma"
DataCurrent.reference="arXiv:1512.02192"

DataCurrent.isNormalized=True
proc_current=[1,1,1,3]
s_current=8000.**2
Q_current=[116.,150.]
y_current=[-2.4,2.4]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.0
#DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
    p["thFactor"]=1/(p["qT"][1]-p["qT"][0])#devide by bin size
    p["Q"]=Q_current
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_from_f[i][38]
    p["uncorrErr"].append((data_from_f[i][39]-data_from_f[i][40])/2.)
    p["uncorrErr"].append((data_from_f[i][41]-data_from_f[i][42])/2.)
    p["corrErr"].append((data_from_f[i][43]-data_from_f[i][44])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

