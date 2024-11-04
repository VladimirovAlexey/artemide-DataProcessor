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
import numpy

path_to_data="/data/arTeMiDe_Repository/data/"
path_to_save="/data/arTeMiDe_Repository/DataProcessor/DataLib/unpolDY/"

M_Z=91.### mass of Z-boson

#%%
###############################################################################
###########################CDF run 1###########################################
print("Read CDF run 1 file ...")
f = open(path_to_data+"CDF_D0/CDF_run1.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:11]
del data_from_f[-1]


for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('CDF1',"DY")
DataCurrent.comment="CDF run1"
DataCurrent.reference="hep-ex/0001021"

DataCurrent.isNormalized=False
proc_current=[1,1,-1,3]
s_current=1800.**2
Q_current=[66.,116.]
y_current=[-1000.,1000.]
lumUncertainty=0.039
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
    p["thFactor"]=1/(p["qT"][1]-p["qT"][0])#devide by bin size
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=False
    p["xSec"]=data_from_f[i][3]
    p["uncorrErr"].append((data_from_f[i][4]-data_from_f[i][5])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################CDF run 2###########################################
print("Read CDF run 2 file ...")
f = open(path_to_data+"CDF_D0/CDF_run2.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:9]
#print 'First line =',data_from_f[0]
#print 'last line =',data_from_f[-1]
del data_from_f[-1]
#print 'last line =',data_from_f[-1]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")

#print data_from_f

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('CDF2',"DY")
DataCurrent.comment="CDF run2"
DataCurrent.reference="arXiv:1207.7138"

proc_current=[1,1,-1,3]
DataCurrent.isNormalized=False
s_current=1960.**2
Q_current=[66.,116.]
y_current=[-1000.,1000.]
lumUncertainty=0.058
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
    p["thFactor"]=1/(p["qT"][1]-p["qT"][0])#devide by bin size
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=False
    p["xSec"]=data_from_f[i][3]
    p["uncorrErr"].append((data_from_f[i][4]-data_from_f[i][5])/2.)
    p["corrErr"].append((data_from_f[i][6]-data_from_f[i][7])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")


###############################################################################
###########################D0 run1 ############################################
print("Read D0 run 1 file ...")
f = open(path_to_data+"CDF_D0/D0_run1.dat")

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

#print data_from_f

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('D01',"DY")
DataCurrent.comment="D0 run1"
DataCurrent.reference="hep-ex/9907009"

DataCurrent.isNormalized=False
proc_current=[1,1,-1,3]
s_current=1800.**2
Q_current=[75.,105.]
y_current=[-1000.,1000.]
lumUncertainty=0.044
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
    p["thFactor"]=1/(p["qT"][1]-p["qT"][0])#devide by bin size
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=False
    p["xSec"]=data_from_f[i][3]
    p["uncorrErr"].append((data_from_f[i][4]-data_from_f[i][5])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

###############################################################################
###########################D0 run2 ############################################
print("Read D0 run 2 file ...")
f = open(path_to_data+"CDF_D0/D0_run2.dat")

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

#print data_from_f

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('D02',"DY")
DataCurrent.comment="D0 run2"
DataCurrent.reference="arXiv:0712.0803"

DataCurrent.isNormalized=True
proc_current=[1,1,-1,3]
s_current=1960.**2
Q_current=[70.,110.]
y_current=[-1000.,1000.]
lumUncertainty=0.0
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
    p["thFactor"]=1/(p["qT"][1]-p["qT"][0])#devide by bin size
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=False
    p["xSec"]=data_from_f[i][3]
    p["uncorrErr"].append((data_from_f[i][4]-data_from_f[i][5])/2.)
    p["uncorrErr"].append((data_from_f[i][6]-data_from_f[i][7])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

###############################################################################
###########################D0 run2 (mu)#########################################
print("Read D0 run 2 (mu) file ...")
f = open(path_to_data+"CDF_D0/D0_run2m.dat")

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

#print data_from_f

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('D02m',"DY")
DataCurrent.comment="D0 run2 data for muons"
DataCurrent.reference="arXiv:1006.0618"

DataCurrent.isNormalized=True
proc_current=[1,1,-1,3]
s_current=1960.**2
Q_current=[65.,115.]
y_current=[-1000.,1000.]
incCut=True
cutParam=[15.,15.,-1.7,1.7]
lumUncertainty=0.0
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
    p["thFactor"]=1/(p["qT"][1]-p["qT"][0])#devide by bin size
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_from_f[i][3]
    p["uncorrErr"].append((data_from_f[i][4]-data_from_f[i][5])/2.)
    p["uncorrErr"].append((data_from_f[i][6]-data_from_f[i][7])/2.)
    p["corrErr"].append((data_from_f[i][8]-data_from_f[i][9])/2.)
    p["corrErr"].append((data_from_f[i][10]-data_from_f[i][11])/2.)
    p["corrErr"].append((data_from_f[i][12]-data_from_f[i][13])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")
DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%

###############################################################################
###########################ATLAS 7 00y10########################################
print("Read ATLAS 7 0.0<y<1.0 file...")
f = open(path_to_data+"ATLAS/ATLAS_7TeV_complete.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:81]
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
DataCurrent=DataProcessor.DataSet.DataSet('A7-00y10',"DY")
DataCurrent.comment="ATLAS 7TeV 0.0<|y|<1.0"
DataCurrent.reference="arXiv:1406.3660"

DataCurrent.isNormalized=True
proc_current=[1,1,1,3]
s_current=7000.**2
Q_current=[66.,116.]
y_current=[0.,1.]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.0
DataCurrent.normErr.append(lumUncertainty)

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
    p["xSec"]=data_from_f[i][3]
    p["uncorrErr"].append((data_from_f[i][3+1]-data_from_f[i][3+2])/2.)
    p["uncorrErr"].append((data_from_f[i][3+3]-data_from_f[i][3+4])/2.)
    p["corrErr"].append((data_from_f[i][3+5]-data_from_f[i][3+6])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")


###############################################################################
###########################ATLAS 7 10y20########################################
print("Read ATLAS 7 1.0<y<2.0 file...")
f = open(path_to_data+"ATLAS/ATLAS_7TeV_complete.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:81]
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
DataCurrent=DataProcessor.DataSet.DataSet('A7-10y20',"DY")
DataCurrent.comment="ATLAS 7TeV 1.0<|y|<2.0"
DataCurrent.reference="arXiv:1406.3660"

DataCurrent.isNormalized=True
proc_current=[1,1,1,3]
s_current=7000.**2
Q_current=[66.,116.]
y_current=[1.,2.]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.0
DataCurrent.normErr.append(lumUncertainty)

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
    p["xSec"]=data_from_f[i][17]
    p["uncorrErr"].append((data_from_f[i][17+1]-data_from_f[i][17+2])/2.)
    p["uncorrErr"].append((data_from_f[i][17+3]-data_from_f[i][17+4])/2.)
    p["corrErr"].append((data_from_f[i][17+5]-data_from_f[i][17+6])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")


###############################################################################
###########################ATLAS 7 20y24########################################
print("Read ATLAS 7 2.0<y<2.4 file...")
f = open(path_to_data+"ATLAS/ATLAS_7TeV_complete.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:81]
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
DataCurrent=DataProcessor.DataSet.DataSet('A7-20y24',"DY")
DataCurrent.comment="ATLAS 7TeV 2.0<|y|<2.4"
DataCurrent.reference="arXiv:1406.3660"

DataCurrent.isNormalized=True
proc_current=[1,1,1,3]
s_current=7000.**2
Q_current=[66.,116.]
y_current=[2.,2.4]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.0
DataCurrent.normErr.append(lumUncertainty)

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
    p["xSec"]=data_from_f[i][31]
    p["uncorrErr"].append((data_from_f[i][31+1]-data_from_f[i][31+2])/2.)
    p["uncorrErr"].append((data_from_f[i][31+3]-data_from_f[i][31+4])/2.)
    p["corrErr"].append((data_from_f[i][31+5]-data_from_f[i][31+6])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################ATLAS 8 00y04########################################
print("Read ATLAS 8 0.0<y<0.4 file...")
f = open(path_to_data+"ATLAS/ATLAS_8TeV_66to116_00y04_abs.dat")

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
DataCurrent=DataProcessor.DataSet.DataSet('A8-00y04',"DY")
DataCurrent.comment="ATLAS 8TeV 0.0<|y|<0.4"
DataCurrent.reference="arXiv:1512.02192"

DataCurrent.isNormalized=False
proc_current=[1,1,1,3]
s_current=8000.**2
Q_current=[66.,116.]
y_current=[0.,0.4]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.028
DataCurrent.normErr.append(lumUncertainty)

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
f = open(path_to_data+"ATLAS/ATLAS_8TeV_66to116_04y08_abs.dat")

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
DataCurrent=DataProcessor.DataSet.DataSet('A8-04y08',"DY")
DataCurrent.comment="ATLAS 8TeV 0.4<|y|<0.8"
DataCurrent.reference="arXiv:1512.02192"

DataCurrent.isNormalized=False
proc_current=[1,1,1,3]
s_current=8000.**2
Q_current=[66.,116.]
y_current=[0.4,0.8]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.028
DataCurrent.normErr.append(lumUncertainty)

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
f = open(path_to_data+"ATLAS/ATLAS_8TeV_66to116_08y12_abs.dat")

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
DataCurrent=DataProcessor.DataSet.DataSet('A8-08y12',"DY")
DataCurrent.comment="ATLAS 8TeV 0.8<|y|<1.2"
DataCurrent.reference="arXiv:1512.02192"

DataCurrent.isNormalized=False
proc_current=[1,1,1,3]
s_current=8000.**2
Q_current=[66.,116.]
y_current=[0.8,1.2]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.028
DataCurrent.normErr.append(lumUncertainty)

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
f = open(path_to_data+"ATLAS/ATLAS_8TeV_66to116_12y16_abs.dat")

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
DataCurrent=DataProcessor.DataSet.DataSet('A8-12y16',"DY")
DataCurrent.comment="ATLAS 8TeV 1.2<|y|<1.6"
DataCurrent.reference="arXiv:1512.02192"

DataCurrent.isNormalized=False
proc_current=[1,1,1,3]
s_current=8000.**2
Q_current=[66.,116.]
y_current=[1.2,1.6]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.028
DataCurrent.normErr.append(lumUncertainty)

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
f = open(path_to_data+"ATLAS/ATLAS_8TeV_66to116_16y20_abs.dat")

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
DataCurrent=DataProcessor.DataSet.DataSet('A8-16y20',"DY")
DataCurrent.comment="ATLAS 8TeV 1.6<|y|<2.0"
DataCurrent.reference="arXiv:1512.02192"

DataCurrent.isNormalized=False
proc_current=[1,1,1,3]
s_current=8000.**2
Q_current=[66.,116.]
y_current=[1.6,2.0]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.028
DataCurrent.normErr.append(lumUncertainty)

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
f = open(path_to_data+"ATLAS/ATLAS_8TeV_66to116_20y24_abs.dat")

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
DataCurrent=DataProcessor.DataSet.DataSet('A8-20y24',"DY")
DataCurrent.comment="ATLAS 8TeV 2.0<|y|<2.4"
DataCurrent.reference="arXiv:1512.02192"

DataCurrent.isNormalized=False
proc_current=[1,1,1,3]
s_current=8000.**2
Q_current=[66.,116.]
y_current=[2.0,2.4]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.028
DataCurrent.normErr.append(lumUncertainty)

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
f = open(path_to_data+"ATLAS/ATLAS_8TeV_46to66_abs.dat")

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
DataCurrent=DataProcessor.DataSet.DataSet('A8-46Q66',"DY")
DataCurrent.comment="ATLAS 8TeV 46<Q<66"
DataCurrent.reference="arXiv:1512.02192"

DataCurrent.isNormalized=False
proc_current=[1,1,1,3]
s_current=8000.**2
Q_current=[46.,66.]
y_current=[-2.4,2.4]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.028
DataCurrent.normErr.append(lumUncertainty)

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
f = open(path_to_data+"ATLAS/ATLAS_8TeV_116to150_abs.dat")

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
DataCurrent=DataProcessor.DataSet.DataSet('A8-116Q150',"DY")
DataCurrent.comment="ATLAS 8TeV 116<Q<150"
DataCurrent.reference="arXiv:1512.02192"

DataCurrent.isNormalized=False
proc_current=[1,1,1,3]
s_current=8000.**2
Q_current=[116.,150.]
y_current=[-2.4,2.4]
incCut=True
cutParam=[20.,20.,-2.4,2.4]
lumUncertainty=0.028
DataCurrent.normErr.append(lumUncertainty)

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
###########################CMS 7TeV ###########################################
print("Read CMS 7TeV file...")
f = open(path_to_data+"CMS/CMS_7.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:11]
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
DataCurrent=DataProcessor.DataSet.DataSet('CMS7',"DY")
DataCurrent.comment="CMS 7TeV"
DataCurrent.reference="arXiv:1110.4973"

DataCurrent.isNormalized=True
proc_current=[1,1,1,3]
s_current=7000.**2
Q_current=[60.,120.]
y_current=[-2.1,2.1]
incCut=True
cutParam=[20.,20.,-2.1,2.1]
lumUncertainty=0.0
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
    p["thFactor"]=1/(p["qT"][1]-p["qT"][0])#devide by bin size
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_from_f[i][9]
    p["uncorrErr"].append((data_from_f[i][10]-data_from_f[i][11])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

###############################################################################
###########################CMS 8TeV ###########################################
print("Read CMS 8TeV file...")
f = open(path_to_data+"CMS/CMS_8.dat")

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
#    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")

#print data_from_f

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('CMS8',"DY")
DataCurrent.comment="CMS 8TeV"
DataCurrent.reference="arXiv:1606.05864"

DataCurrent.isNormalized=True
proc_current=[1,1,1,3]
s_current=8000.**2
Q_current=[60.,120.]
y_current=[-2.1,2.1]
incCut=True
cutParam=[20.,20.,-2.1,2.1]
lumUncertainty=0.0
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current
    p["qT"]=[data_from_f[i][0]-data_from_f[i][1],data_from_f[i][0]+data_from_f[i][2]]
    p["thFactor"]=1/(p["qT"][1]-p["qT"][0])#devide by bin size
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_from_f[i][3]
    p["uncorrErr"].append((data_from_f[i][4]+data_from_f[i][5])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")


#%%
###############################################################################
###########################LHCb 7TeV ###########################################
print("Read LHCb 7TeV file...")
f = open(path_to_data+"LHCb/LHCb_7.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:11]
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
DataCurrent=DataProcessor.DataSet.DataSet('LHCb7',"DY")
DataCurrent.comment="LHCb 7TeV"
DataCurrent.reference="arXiv:1505.07024"

DataCurrent.isNormalized=False
proc_current=[1,1,1,3]
s_current=7000.**2
Q_current=[60.,120.]
y_current=[2.,4.5]
incCut=True
cutParam=[20.,20.,2.,4.5]
lumUncertainty=0.0
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current    
    p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
    p["thFactor"]=1.#this data is not weighted by bin size
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_from_f[i][3]
    p["uncorrErr"].append((data_from_f[i][4]-data_from_f[i][5])/2.)
    p["uncorrErr"].append((data_from_f[i][6]-data_from_f[i][7])/2.)
    p["corrErr"].append((data_from_f[i][8]-data_from_f[i][9])/2.)
    p["corrErr"].append((data_from_f[i][10]-data_from_f[i][11])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")


###############################################################################
###########################LHCb 8TeV ###########################################
print("Read LHCb 8TeV file...")
f = open(path_to_data+"LHCb/LHCb_8.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:11]
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
DataCurrent=DataProcessor.DataSet.DataSet('LHCb8',"DY")
DataCurrent.comment="LHCb 8TeV"
DataCurrent.reference="arXiv:1511.08039"

DataCurrent.isNormalized=False
proc_current=[1,1,1,3]
s_current=8000.**2
Q_current=[60.,120.]
y_current=[2.,4.5]
incCut=True
cutParam=[20.,20.,2.,4.5]
lumUncertainty=0.0
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current    
    p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
    p["thFactor"]=1.#this data is not weighted by bin size
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_from_f[i][3]
    p["uncorrErr"].append((data_from_f[i][4]-data_from_f[i][5])/2.)
    p["uncorrErr"].append((data_from_f[i][6]-data_from_f[i][7])/2.)
    p["corrErr"].append((data_from_f[i][8]-data_from_f[i][9])/2.)
    p["corrErr"].append((data_from_f[i][10]-data_from_f[i][11])/2.)
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")


###############################################################################
###########################LHCb 13TeV ###########################################
print("Read LHCb 13TeV file...")
f = open(path_to_data+"LHCb/LHCb_13.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:12]
#print 'First line =',data_from_f[0]
#print 'last line (before)=',data_from_f[-1]
#del data_from_f[-1]
#print 'last line (after) =',data_from_f[-1]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split("\t")
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")

#print data_from_f

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('LHCb13',"DY")
DataCurrent.comment="LHCb 13TeV"
DataCurrent.reference="arXiv:1607.06495"

DataCurrent.isNormalized=False
proc_current=[1,1,1,3]
s_current=13000.**2
Q_current=[60.,120.]
y_current=[2.,4.5]
incCut=True
cutParam=[20.,20.,2.,4.5]
lumUncertainty=0.0
DataCurrent.normErr.append(lumUncertainty)

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))
    p["process"]=proc_current
    p["s"]=s_current    
    p["qT"]=[data_from_f[i][1],data_from_f[i][2]]
    p["thFactor"]=1/(p["qT"][1]-p["qT"][0])#devide by bin size
    p["Q"]=Q_current
    p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
    p["y"]=y_current
    p["includeCuts"]=incCut
    p["cutParams"]=cutParam
    p["xSec"]=data_from_f[i][3]
    p["uncorrErr"].append(data_from_f[i][4])
    p["uncorrErr"].append(data_from_f[i][5])
    p["corrErr"].append(data_from_f[i][6])
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")


#%%
###############################################################################
###########################PHENIX 200##########################################
###### NOW IN A SEPARATE FILE 

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
proc_current=[2,1,-1,103]
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