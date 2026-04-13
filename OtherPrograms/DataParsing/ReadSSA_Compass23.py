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

########### Alexey desktop
path_to_data="/data/arTeMiDe_Repository/data"
path_to_COMPASS="/COMPASS/2401.00309/"
path_to_save="/data/arTeMiDe_Repository/DataProcessor/DataLib/Sivers/"

########## Marcin laptop
#path_to_data="/home/m/Dropbox/Sivers/Data"
#path_to_COMPASS="/COMPASS08/"

totalData=[]

M_proton=0.938
m_pion=0.139
m_kaon=0.494
#%%

#### determines the limits of the Q bin
#### Qmin^2 = MAX (Q^2min, xmin y min (s-M^2), xmin/(1-xmin)*(W2min-M^2))
#### Qmax^2 = MIN (Q^2max, xmax y max (s-M^2), xmax/(1-xmax)*(W2max-M^2))
def Qbounds(xMin,xMax):
    Q2min=1.
    Q2max=10000.
    WM2min=25.-(0.938)**2
    WM2max=10000.-(0.938)**2  ## no upper limit
    yMin=0.1
    yMax=0.9
    sM2=2.*160.*0.938
    
    if xMax<1:
        return [numpy.sqrt(numpy.max([Q2min,xMin*yMin*sM2, xMin/(1-xMin)*WM2min])),
                numpy.sqrt(numpy.min([Q2max,xMax*yMax*sM2, xMax/(1-xMax)*WM2max]))]
    else:
        return [numpy.sqrt(numpy.max([Q2min,xMin*yMin*sM2, xMin/(1-xMin)*WM2min])),
                numpy.sqrt(numpy.min([Q2max,xMax*yMax*sM2, 1000*WM2max]))]

#%%
###############################################################################
###########################hermes.proton.zxpt-3D pi+###########################
print("compass.sivers file ...")
f = open(path_to_data+path_to_COMPASS+"/Sivers.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")



#%%
###### POSITIVE HADRONS H+ ######
data_current=data_from_f[6:15]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";",",").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass23.sivers.h+.dx',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers h+ (differential in x)"
DataCurrentSiv.reference="2401.00309"

proc_current=[1,1,1,12103]
proc_denominator=[1,1,1,2103]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[0.1,1.6]        
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[data_current[i][0],data_current[i][1]]
    p1["<z>"]=data_current[i][4]
    p1["z"]=[0.2,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][7])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][8]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    
    #
    DataCurrentSiv.AddPoint(p1)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
data_current=data_from_f[19:27]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";",",").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass23.sivers.h+.dz',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers h+ (differential in z)"
DataCurrentSiv.reference="2401.00309"


for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[0.1,1.6]        
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[0.003,0.7]
    p1["<z>"]=data_current[i][4]
    p1["z"]=[data_current[i][0],data_current[i][1]]
    p1["<Q>"]=numpy.sqrt(data_current[i][7])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][8]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    
    #
    DataCurrentSiv.AddPoint(p1)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS H+ ######
data_current=data_from_f[31:40]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";",",").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass23.sivers.h+.dpt',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers h+ (differential in pt)"
DataCurrentSiv.reference="2401.00309"


for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[data_current[i][0],data_current[i][1]]
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[0.003,0.7]
    p1["<z>"]=data_current[i][4]
    p1["z"]=[0.2,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][7])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][8]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    
    #
    DataCurrentSiv.AddPoint(p1)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS H- ######
data_current=data_from_f[44:53]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";",",").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass23.sivers.h-.dx',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers h- (differential in x)"
DataCurrentSiv.reference="2401.00309"

proc_current=[1,1,-1,12103]
proc_denominator=[1,1,-1,2103]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[0.1,1.6]        
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[data_current[i][0],data_current[i][1]]
    p1["<z>"]=data_current[i][4]
    p1["z"]=[0.2,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][7])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][8]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    
    #
    DataCurrentSiv.AddPoint(p1)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
data_current=data_from_f[57:65]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";",",").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass23.sivers.h-.dz',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers h- (differential in z)"
DataCurrentSiv.reference="2401.00309"


for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[0.1,1.6]        
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[0.003,0.7]
    p1["<z>"]=data_current[i][4]
    p1["z"]=[data_current[i][0],data_current[i][1]]
    p1["<Q>"]=numpy.sqrt(data_current[i][7])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][8]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    
    #
    DataCurrentSiv.AddPoint(p1)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS H- ######
data_current=data_from_f[69:78]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";",",").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass23.sivers.h-.dpt',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers h- (differential in pt)"
DataCurrentSiv.reference="2401.00309"


for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[data_current[i][0],data_current[i][1]]
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[0.003,0.7]
    p1["<z>"]=data_current[i][4]
    p1["z"]=[0.2,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][7])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][8]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    
    #
    DataCurrentSiv.AddPoint(p1)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")