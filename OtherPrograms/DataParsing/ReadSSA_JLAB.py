#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 20:15:46 2019
This program collect all the data on the SIDIS and save in  "SIDISdata_uncut.pkl"
@author: vla18041
"""
import sys
sys.path.append("/data/arTeMiDe_Repository/DataProcessor/")
#sys.path.append("/home/m/Github/artemide-DataProcessor/")
import DataProcessor.Point
import DataProcessor.DataSet
import numpy

########### Alexey desktop
path_to_data="/data/arTeMiDe_Repository/data"
path_to_JLab="/JLab-Sivers/"
path_to_save="/data/arTeMiDe_Repository/DataProcessor/DataLib/Sivers/"

########## Marcin laptop
#path_to_data="/home/m/Dropbox/Sivers/Data"
#path_to_HERMES="/HERMES09/"
#path_to_save="/home/m/Github/artemide-DataProcessor/DataLib/Sivers/"

totalData=[]

M_proton=0.938
m_pion=0.139
m_kaon=0.494

###############################################################################
#################JLab data is exeptionally badly documented####################
#########################nothing is known######################################

#%%

### bins are looked by eye from the 1404.7204 paper
### There should be taken the central VALUE!!! 
xBin=[0.1,0.4]
zBin=[0.45,0.6]
ptBin=[0.1,0.6]
QBin=[1.,1.7]
### Scale uncertanty ??
#scaleUncertanty=0.073

#%%
###############################################################################
f = open(path_to_data+path_to_JLab+"Siv_pi+")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### POSITIVE PIONS PI+ ######
data_current=data_from_f[1:5]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace(",",".").replace("\t\t\t","")
    data_current[i]=data_current[i].split("\t")    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('jlab.sivers.pi+',"SIDIS")
DataCurrent.comment="JLab SSA-Sivers pi+. The data MUST be evaluated at a point"
DataCurrent.reference="1106.0363"

proc_current=[1,1,1,12003]#neutron target
proc_denominator=[1,1,1,2003]
s_current=2*5.9*0.938+(0.938)**2
includeCuts=False
cutParameters=[0.1,0.95,2.3,10000.] #y, W^2 cuts
#DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=ptBin
    p1["<x>"]=data_current[i][1]    
    p1["x"]=xBin
    p1["<Q>"]=numpy.sqrt(data_current[i][4])
    p1["Q"]=QBin
    p1["<z>"]=data_current[i][3]
    p1["z"]=zBin
    p1["xSec"]=data_current[i][6]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][7])
    p1["uncorrErr"].append(data_current[i][8])
    p1["weightProcess"]=proc_denominator
    #
    DataCurrent.AddPoint(p1)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_JLab+"Siv_pi-")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### POSITIVE PIONS PI- ######
data_current=data_from_f[1:5]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace(",",".").replace("\t\t\t","")
    data_current[i]=data_current[i].split("\t")    
    data_current[i]=[float(j) for j in data_current[i]]
#%%
print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('jlab.sivers.pi-',"SIDIS")
DataCurrent.comment="JLab SSA-Sivers pi-. The data MUST be evaluated at a point"
DataCurrent.reference="1106.0363"

proc_current=[1,1,-1,12003]#neutraon target
proc_denominator=[1,1,-1,2003]
s_current=2*5.9*0.938+(0.938)**2
includeCuts=False
cutParameters=[0.1,0.95,2.3,10000.] #y, W^2 cuts
#DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=ptBin
    p1["<x>"]=data_current[i][1]    
    p1["x"]=xBin
    p1["<Q>"]=numpy.sqrt(data_current[i][4])
    p1["Q"]=QBin
    p1["<z>"]=data_current[i][3]
    p1["z"]=zBin
    p1["xSec"]=data_current[i][6]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][7])
    p1["uncorrErr"].append(data_current[i][8])
    p1["weightProcess"]=proc_denominator
    #
    DataCurrent.AddPoint(p1)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###### POSITIVE PIONS K+ ######

#%%
###############################################################################
f = open(path_to_data+path_to_JLab+"Siv_k+")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

data_current=data_from_f[1:5]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace(",",".").replace("\t\t\t","")
    data_current[i]=data_current[i].split("\t")    
    data_current[i]=[float(j) for j in data_current[i]]
#%%
print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('jlab.sivers.k+',"SIDIS")
DataCurrent.comment="JLab SSA-Sivers k+. The data MUST be evaluated at a point"
DataCurrent.reference="1404.7204"

proc_current=[1,1,2,12003]#neutron target
proc_denominator=[1,1,2,2003]
s_current=2*5.9*0.938+(0.938)**2
includeCuts=False
cutParameters=[0.1,0.95,2.3,10000.] #y, W^2 cuts
#DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=ptBin
    p1["<x>"]=data_current[i][1]    
    p1["x"]=xBin
    p1["<Q>"]=numpy.sqrt(data_current[i][4])
    p1["Q"]=QBin
    p1["<z>"]=data_current[i][3]
    p1["z"]=zBin
    p1["xSec"]=data_current[i][8]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][9])
    p1["uncorrErr"].append(data_current[i][11])
    p1["weightProcess"]=proc_denominator
    #
    DataCurrent.AddPoint(p1)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###### POSITIVE PIONS K- ######

###############################################################################
f = open(path_to_data+path_to_JLab+"Siv_k-")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

data_current=data_from_f[1:2]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace(",",".").replace("\t\t\t","")
    data_current[i]=data_current[i].split("\t")    
    data_current[i]=[float(j) for j in data_current[i]]
#%%
print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('jlab.sivers.k-',"SIDIS")
DataCurrent.comment="JLab SSA-Sivers k-. The data MUST be evaluated at a point"
DataCurrent.reference="1404.7204"

proc_current=[1,1,-2,12003]#neutron target
proc_denominator=[1,1,-2,2003]
s_current=2*5.9*0.938+(0.938)**2
includeCuts=False
cutParameters=[0.1,0.95,2.3,10000.] #y, W^2 cuts
#DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=ptBin
    p1["<x>"]=data_current[i][1]    
    p1["x"]=xBin
    p1["<Q>"]=numpy.sqrt(data_current[i][4])
    p1["Q"]=QBin
    p1["<z>"]=data_current[i][3]
    p1["z"]=zBin
    p1["xSec"]=data_current[i][8]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][9])
    p1["uncorrErr"].append(data_current[i][11])
    p1["weightProcess"]=proc_denominator
    #
    DataCurrent.AddPoint(p1)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")