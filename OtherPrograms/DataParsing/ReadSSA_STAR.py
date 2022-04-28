#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 20:15:46 2019
This program collect all the data on the SIDIS and save in  "SIDISdata_uncut.pkl"
@author: vla18041
"""
import sys
sys.path.append("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/")
import DataProcessor.Point
import DataProcessor.DataSet
import numpy

########### Alexey desktop
path_to_data="/home/vla18041/LinkData2/arTeMiDe_Repository/data"
path_to_STAR="/STAR-Sivers/"
path_to_save="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/DataLib/Sivers/"

########## Marcin laptop
#path_to_data="/home/m/Dropbox/Sivers/Data"
#path_to_COMPASS="/COMPASS08/"

M_Z=91.19
M_W=80.38

### Scale uncertanty
beamPolarizationUncertanty=0.034

#%%
###############################################################################
f = open(path_to_data+path_to_STAR+"Table4.csv")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### PW+ ######
######    d PT             ######
data_current=data_from_f[12:17]

for i in range(len(data_current)):
    #print(data_current[i])
    data_current[i]=data_current[i].split(",")    
    data_current[i]=[float(j) for j in data_current[i]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('star.sivers.W+.dqT',"DY")
DataCurrent.comment="STAR AN for W+ (differential in qT)"
DataCurrent.reference="1511.06003"

proc_current=[1,1,10007]
proc_denominator=[1,1,7]
s_current=500.**2
includeCuts=False
cutParameters=[25.,25.,-1.,1.] #pT1,pT2,yMin,yMax
DataCurrent.normErr.append(beamPolarizationUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))    
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<qT>"]=data_current[i][0]    
    p1["qT"]=[data_current[i][1],data_current[i][2]]    
    p1["<Q>"]=M_W
    p1["Q"]=[M_W-40.,M_W+40.] ### I have checked that +- 40 gives ~10^{-4} error in comparison to the total volume
    p1["y"]=[-1.,1.]
    p1["xSec"]=data_current[i][3]
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append((data_current[i][4]-data_current[i][5])/2.)
    p1["uncorrErr"].append((data_current[i][6]-data_current[i][7])/2.)
    p1["weightProcess"]=proc_denominator
    #
    DataCurrent.AddPoint(p1)    

print("Done.  ")

#DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###### PW- ######
######    d PT             ######
data_current=data_from_f[21:26]

for i in range(len(data_current)):
    #print(data_current[i])
    data_current[i]=data_current[i].split(",")    
    data_current[i]=[float(j) for j in data_current[i]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('star.sivers.W-.dqT',"DY")
DataCurrent.comment="STAR AN for W+ (differential in qT)"
DataCurrent.reference="1511.06003"

proc_current=[1,1,10008]
proc_denominator=[1,1,8]
s_current=500.**2
includeCuts=False
cutParameters=[25.,25.,-1.,1.] #pT1,pT2,yMin,yMax
DataCurrent.normErr.append(beamPolarizationUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))    
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<qT>"]=data_current[i][0]    
    p1["qT"]=[data_current[i][1],data_current[i][2]]    
    p1["<Q>"]=M_W
    p1["Q"]=[M_W-40.,M_W+40.] ### I have checked that +- 40 gives ~10^{-4} error in comparison to the total volume
    p1["y"]=[-1.,1.]
    p1["xSec"]=data_current[i][3]
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append((data_current[i][4]-data_current[i][5])/2.)
    p1["uncorrErr"].append((data_current[i][6]-data_current[i][7])/2.)
    p1["weightProcess"]=proc_denominator
    #
    DataCurrent.AddPoint(p1)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_STAR+"Table5.csv")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### PW+ ######
######    d y            ######
data_current=data_from_f[12:15]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split(",")    
    data_current[i]=[float(j) for j in data_current[i]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('star.sivers.W+.dy',"DY")
DataCurrent.comment="STAR AN for W+ (differential in y)"
DataCurrent.reference="1511.06003"

proc_current=[1,1,10007]
proc_denominator=[1,1,7]
s_current=500.**2
includeCuts=False
cutParameters=[25.,25.,-1.,1.] #pT1,pT2,yMin,yMax
DataCurrent.normErr.append(beamPolarizationUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))    
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["qT"]=[0.5,10.]    
    p1["<Q>"]=M_W
    p1["Q"]=[M_W-30.,M_W+30.] ### I have checked that +- 30 gives ~10^{-3} error in comparison to the total volume
    p1["<y>"]=data_current[i][0]
    p1["y"]=[data_current[i][1],data_current[i][2]]    
    p1["xSec"]=data_current[i][3]
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append((data_current[i][4]-data_current[i][5])/2.)
    p1["uncorrErr"].append((data_current[i][6]-data_current[i][7])/2.)
    p1["weightProcess"]=proc_denominator
    #
    DataCurrent.AddPoint(p1)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###### P-+ ######
######    d y            ######
data_current=data_from_f[19:22]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split(",")    
    data_current[i]=[float(j) for j in data_current[i]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('star.sivers.W-.dy',"DY")
DataCurrent.comment="STAR AN for W- (differential in y)"
DataCurrent.reference="1511.06003"

proc_current=[1,1,10008]
proc_denominator=[1,1,8]
s_current=500.**2
includeCuts=False
cutParameters=[25.,25.,-1.,1.] #pT1,pT2,yMin,yMax
DataCurrent.normErr.append(beamPolarizationUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))    
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["qT"]=[0.5,10.]    
    p1["<Q>"]=M_W
    p1["Q"]=[M_W-30.,M_W+30.] ### I have checked that +- 30 gives ~10^{-3} error in comparison to the total volume
    p1["<y>"]=data_current[i][0]
    p1["y"]=[data_current[i][1],data_current[i][2]]    
    p1["xSec"]=data_current[i][3]
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append((data_current[i][4]-data_current[i][5])/2.)
    p1["uncorrErr"].append((data_current[i][6]-data_current[i][7])/2.)
    p1["weightProcess"]=proc_denominator
    #
    DataCurrent.AddPoint(p1)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_STAR+"Table6.csv")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### PZ ######
######    1 -point            ######
data_current=data_from_f[13:14]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split(",")    
    data_current[i]=[float(j) for j in data_current[i]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('star.sivers.Z',"DY")
DataCurrent.comment="STAR AN for Z"
DataCurrent.reference="1511.06003"

proc_current=[1,1,10005]
proc_denominator=[1,1,5]
s_current=500.**2
includeCuts=False
cutParameters=[25.,25.,-1.,1.] #pT1,pT2,yMin,yMax
DataCurrent.normErr.append(beamPolarizationUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(i))    
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["qT"]=[0.5,10.]    
    p1["<Q>"]=M_Z
    p1["Q"]=[M_Z-30.,M_Z+30.] ### I have checked that +- 30 gives ~10^{-3} error in comparison to the total volume
    p1["<y>"]=data_current[i][0]
    p1["y"]=[data_current[i][1],data_current[i][2]]    
    p1["xSec"]=data_current[i][3]
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append((data_current[i][4]-data_current[i][5])/2.)
    p1["weightProcess"]=proc_denominator
    #
    DataCurrent.AddPoint(p1)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
