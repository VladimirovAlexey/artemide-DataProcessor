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
path_to_COMPASS="/COMPASS/1609.07374/"
path_to_save="/data/arTeMiDe_Repository/DataProcessor/DataLib/Sivers/"

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
    WM2min=10.-(0.938)**2
    WM2max=10000.-(0.938)**2  ## no upper limit
    yMin=0.1
    yMax=0.9
    sM2=2*160*0.938
    
    if xMax<1:
        return [numpy.sqrt(numpy.max([Q2min,xMin*yMin*sM2, xMin/(1-xMin)*WM2min])),
                numpy.sqrt(numpy.min([Q2max,xMax*yMax*sM2, xMax/(1-xMax)*WM2max]))]
    else:
        return [numpy.sqrt(numpy.max([Q2min,xMin*yMin*sM2, xMin/(1-xMin)*WM2min])),
                numpy.sqrt(numpy.min([Q2max,xMax*yMax*sM2, 1000*WM2max]))]
#%%
#### determines the limits of the x bin
#### xmin = MAX { xmin, Q^2/(ymax (s-M^2)), Q^2/(Q^2+W2max-M^2) }
#### xmax = MIN { xmax, Q^2/(ymin (s-M^2)), Q^2/(Q^2+W2min-M^2) }
def xbounds(Q2min,Q2max):
    xmin=0.003
    xmax=0.9  
    WM2min=10.-(0.938)**2
    WM2max=10000.-(0.938)**2  ## no upper limit
    yMin=0.1
    yMax=0.9
    sM2=2*160*0.938
    
    return [numpy.max( [xmin,Q2min/(yMax*sM2), Q2min/(Q2min+WM2max)] ),
            numpy.min( [xmax,Q2max/(yMin*sM2), Q2max/(Q2max+WM2min)] ) ]
    
#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"Zgt01ls02_n_pt_Siv.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### NEGATIVE HADRONS h- 0.1<z<0.2, 1<Q<2 ######1
data_current=data_from_f[6:11]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")


print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass16.sivers.h-.1<z<2.dpt',"SIDIS")
DataCurrentSiv.comment="COMPASS16 SSA-Sivers h-, 0.1<z<0.2 (differential in pt)"
DataCurrentSiv.reference="1609.07374"


proc_current=[1,1,-1,12101]
proc_denominator=[1,1,-1,2101]
s_current=2*160*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

Q2_current=[1.,4.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.003,0.145]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,0.2]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- 0.1<z<0.2, 2<Q<2.5 ######2
data_current=data_from_f[15:20]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

Q2_current=[4.,6.25]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.014,0.215]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,0.2]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- 0.1<z<0.2, 2.5<Q<4 ######3
data_current=data_from_f[24:29]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[6.25,16.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.022,0.550]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,0.2]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- 0.1<z<0.2, 4<Q<9 ######4
data_current=data_from_f[33:38]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[16.,81.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.055,0.8]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,0.2]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"Zgt01ls02_n_x_Siv.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### NEGATIVE HADRONS h- 0.1<z<0.2, 1<Q<2 ######5
data_current=data_from_f[6:13]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass16.sivers.h-.1<z<2.dx',"SIDIS")
DataCurrentSiv.comment="COMPASS16 SSA-Sivers h-, 0.1<z<0.2 (differential in x)"
DataCurrentSiv.reference="1609.07374"


proc_current=[1,1,-1,12101]
proc_denominator=[1,1,-1,2101]
s_current=2*160*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

Q2_current=[1.,4.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,0.2]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- 0.1<z<0.2, 2<Q<2.5 ######6
data_current=data_from_f[17:23]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

Q2_current=[4.,6.25]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,0.2]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p) 
    
print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- 0.1<z<0.2, 2.5<Q<4 ######7
data_current=data_from_f[27:34]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[6.25,16.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]


for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,0.2]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p) 
    
print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- 0.1<z<0.2, 4<Q<9 ######8
data_current=data_from_f[38:44]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[16.,81.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]


for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,0.2]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p) 
    
print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"Zgt01ls02_n_z_Siv.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### NEGATIVE HADRONS h- 0.1<z<0.2, 1<Q<2 ######9
data_current=data_from_f[6:11]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass16.sivers.h-.1<z<2.dz',"SIDIS")
DataCurrentSiv.comment="COMPASS16 SSA-Sivers h-, 0.1<z<0.2 (differential in z)"
DataCurrentSiv.reference="1609.07374"


proc_current=[1,1,-1,12101]
proc_denominator=[1,1,-1,2101]
s_current=2*160*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

Q2_current=[1.,4.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]     
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.003,0.145]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]]   
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- 0.1<z<0.2, 2<Q<2.5 ######10
data_current=data_from_f[15:20]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

Q2_current=[4.,6.25]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.014,0.215]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]] 
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- 0.1<z<0.2, 2.5<Q<4 ######11
data_current=data_from_f[24:29]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
   
Q2_current=[6.25,16.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]      
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.022,0.550]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]]   
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- 0.1<z<0.2, 4<Q<9 ######12
data_current=data_from_f[33:38]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
   
Q2_current=[16.,81.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]      
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.055,0.8]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]]   
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"Zgt01_n_pt_Siv.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### NEGATIVE HADRONS h- z>0.1, 1<Q<2 ######13
data_current=data_from_f[6:11]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")


print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass16.sivers.h-.z>1.dpt',"SIDIS")
DataCurrentSiv.comment="COMPASS16 SSA-Sivers h-, z>0.1 (differential in pt)"
DataCurrentSiv.reference="1609.07374"


proc_current=[1,1,-1,12101]
proc_denominator=[1,1,-1,2101]
s_current=2*160*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

Q2_current=[1.,4.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.003,0.145]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- z>0.1, 2<Q<2.5 ######14
data_current=data_from_f[15:20]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

Q2_current=[4.,6.25]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.014,0.215]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- z>0.1, 2.5<Q<4 ######15
data_current=data_from_f[24:29]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[6.25,16.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.022,0.550]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- z>0.1, 4<Q<9 ######16
data_current=data_from_f[33:38]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[16.,81.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.055,0.8]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"Zgt01_n_x_Siv.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### NEGATIVE HADRONS h- z>0.1, 1<Q<2 ######17
data_current=data_from_f[6:13]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass16.sivers.h-.z>1.dx',"SIDIS")
DataCurrentSiv.comment="COMPASS16 SSA-Sivers h-, z>0.1 (differential in x)"
DataCurrentSiv.reference="1609.07374"


proc_current=[1,1,-1,12101]
proc_denominator=[1,1,-1,2101]
s_current=2*160*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

Q2_current=[1.,4.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- z>0.1, 2<Q<2.5 ######18
data_current=data_from_f[17:23]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

Q2_current=[4.,6.25]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p) 
    
print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- z>0.1, 2.5<Q<4 ######19
data_current=data_from_f[27:34]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[6.25,16.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]


for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p) 
    
print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- z>0.1, 4<Q<9 ######20
data_current=data_from_f[38:44]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[16.,81.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]


for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p) 
    
print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"Zgt01_n_z_Siv.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### NEGATIVE HADRONS h- z>0.1, 1<Q<2 ######21
data_current=data_from_f[6:11]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass16.sivers.h-.z>1.dz',"SIDIS")
DataCurrentSiv.comment="COMPASS16 SSA-Sivers h-, z>0.1 (differential in z)"
DataCurrentSiv.reference="1609.07374"


proc_current=[1,1,-1,12101]
proc_denominator=[1,1,-1,2101]
s_current=2*160*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

Q2_current=[1.,4.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]     
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.003,0.145]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]]   
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- z>0.1, 2<Q<2.5 ######22
data_current=data_from_f[15:20]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

Q2_current=[4.,6.25]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.014,0.215]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]] 
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- z>0.1, 2.5<Q<4 ######23
data_current=data_from_f[24:29]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
   
Q2_current=[6.25,16.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]      
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.022,0.550]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]]   
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- z>0.1, 4<Q<9 ######24
data_current=data_from_f[33:38]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
   
Q2_current=[16.,81.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]      
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.055,0.8]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]]   
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"Zgt02_n_pt_Siv.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### NEGATIVE HADRONS h- z>0.2, 1<Q<2 ######25
data_current=data_from_f[6:11]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")


print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass16.sivers.h-.z>2.dpt',"SIDIS")
DataCurrentSiv.comment="COMPASS16 SSA-Sivers h-, z>0.2 (differential in pt)"
DataCurrentSiv.reference="1609.07374"


proc_current=[1,1,-1,12101]
proc_denominator=[1,1,-1,2101]
s_current=2*160*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

Q2_current=[1.,4.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.003,0.145]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.2,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- z>0.2, 2<Q<2.5 ######26
data_current=data_from_f[15:20]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

Q2_current=[4.,6.25]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.014,0.215]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.2,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- z>0.2, 2.5<Q<4 ######27
data_current=data_from_f[24:29]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[6.25,16.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.022,0.550]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.2,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- z>0.2, 4<Q<9 ######28
data_current=data_from_f[33:38]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[16.,81.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.055,0.8]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.2,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"Zgt02_n_x_Siv.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### NEGATIVE HADRONS h- z>0.2, 1<Q<2 ######29
data_current=data_from_f[6:13]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass16.sivers.h-.z>2.dx',"SIDIS")
DataCurrentSiv.comment="COMPASS16 SSA-Sivers h-, z>0.2 (differential in x)"
DataCurrentSiv.reference="1609.07374"


proc_current=[1,1,-1,12101]
proc_denominator=[1,1,-1,2101]
s_current=2*160*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

Q2_current=[1.,4.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.2,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- z>0.2, 2<Q<2.5 ######30
data_current=data_from_f[17:23]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

Q2_current=[4.,6.25]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.2,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p) 
    
print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- z>0.2, 2.5<Q<4 ######31
data_current=data_from_f[27:34]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[6.25,16.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]


for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.2,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p) 
    
print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- z>0.2, 4<Q<9 ######32
data_current=data_from_f[38:44]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[16.,81.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]


for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.2,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p) 
    
print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"Zgt02_n_z_Siv.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### NEGATIVE HADRONS h- z>0.2, 1<Q<2 ######33
data_current=data_from_f[6:11]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass16.sivers.h-.z>2.dz',"SIDIS")
DataCurrentSiv.comment="COMPASS16 SSA-Sivers h-, z>0.2 (differential in z)"
DataCurrentSiv.reference="1609.07374"


proc_current=[1,1,-1,12101]
proc_denominator=[1,1,-1,2101]
s_current=2*160*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

Q2_current=[1.,4.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]     
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.003,0.145]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]]   
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- z>0.2, 2<Q<2.5 ######34
data_current=data_from_f[15:20]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

Q2_current=[4.,6.25]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.014,0.215]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]] 
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- z>0.2, 2.5<Q<4 ######35
data_current=data_from_f[24:29]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
   
Q2_current=[6.25,16.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]      
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.022,0.550]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]]   
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### NEGATIVE HADRONS h- z>0.2, 4<Q<9 ######36
data_current=data_from_f[33:38]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
   
Q2_current=[16.,81.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]      
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.055,0.8]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]]   
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"Zgt01ls02_p_pt_Siv.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### POSITIVE HADRONS h+ 0.1<z<0.2, 1<Q<2 ######1
data_current=data_from_f[6:11]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")


print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass16.sivers.h+.1<z<2.dpt',"SIDIS")
DataCurrentSiv.comment="COMPASS16 SSA-Sivers h+, 0.1<z<0.2 (differential in pt)"
DataCurrentSiv.reference="1609.07374"


proc_current=[1,1,1,12101]
proc_denominator=[1,1,1,2101]
s_current=2*160*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

Q2_current=[1.,4.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.003,0.145]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,0.2]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ 0.1<z<0.2, 2<Q<2.5 ######2
data_current=data_from_f[15:20]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

Q2_current=[4.,6.25]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.014,0.215]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,0.2]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ 0.1<z<0.2, 2.5<Q<4 ######3
data_current=data_from_f[24:29]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[6.25,16.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.022,0.550]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,0.2]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ 0.1<z<0.2, 4<Q<9 ######4
data_current=data_from_f[33:38]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[16.,81.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.055,0.8]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,0.2]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"Zgt01ls02_p_x_Siv.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### POSITIVE HADRONS h+ 0.1<z<0.2, 1<Q<2 ######5
data_current=data_from_f[6:13]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass16.sivers.h+.1<z<2.dx',"SIDIS")
DataCurrentSiv.comment="COMPASS16 SSA-Sivers h+, 0.1<z<0.2 (differential in x)"
DataCurrentSiv.reference="1609.07374"


proc_current=[1,1,1,12101]
proc_denominator=[1,1,1,2101]
s_current=2*160*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

Q2_current=[1.,4.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,0.2]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ 0.1<z<0.2, 2<Q<2.5 ######6
data_current=data_from_f[17:23]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

Q2_current=[4.,6.25]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,0.2]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p) 
    
print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ 0.1<z<0.2, 2.5<Q<4 ######7
data_current=data_from_f[27:34]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[6.25,16.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]


for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,0.2]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p) 
    
print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ 0.1<z<0.2, 4<Q<9 ######8
data_current=data_from_f[38:44]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[16.,81.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]


for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,0.2]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p) 
    
print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"Zgt01ls02_p_z_Siv.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### POSITIVE HADRONS h+ 0.1<z<0.2, 1<Q<2 ######9
data_current=data_from_f[6:11]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass16.sivers.h+.1<z<2.dz',"SIDIS")
DataCurrentSiv.comment="COMPASS16 SSA-Sivers h+, 0.1<z<0.2 (differential in z)"
DataCurrentSiv.reference="1609.07374"


proc_current=[1,1,1,12101]
proc_denominator=[1,1,1,2101]
s_current=2*160*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

Q2_current=[1.,4.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]     
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.003,0.145]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]]   
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ 0.1<z<0.2, 2<Q<2.5 ######10
data_current=data_from_f[15:20]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

Q2_current=[4.,6.25]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.014,0.215]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]] 
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ 0.1<z<0.2, 2.5<Q<4 ######11
data_current=data_from_f[24:29]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
   
Q2_current=[6.25,16.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]      
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.022,0.550]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]]   
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ 0.1<z<0.2, 4<Q<9 ######12
data_current=data_from_f[33:38]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
   
Q2_current=[16.,81.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]      
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.055,0.8]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]]   
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"Zgt01_p_pt_Siv.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### POSITIVE HADRONS h+ z>0.1, 1<Q<2 ######13
data_current=data_from_f[6:11]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")


print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass16.sivers.h+.z>1.dpt',"SIDIS")
DataCurrentSiv.comment="COMPASS16 SSA-Sivers h+, z>0.1 (differential in pt)"
DataCurrentSiv.reference="1609.07374"


proc_current=[1,1,1,12101]
proc_denominator=[1,1,1,2101]
s_current=2*160*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

Q2_current=[1.,4.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.003,0.145]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ z>0.1, 2<Q<2.5 ######14
data_current=data_from_f[15:20]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

Q2_current=[4.,6.25]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.014,0.215]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ z>0.1, 2.5<Q<4 ######15
data_current=data_from_f[24:29]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[6.25,16.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.022,0.550]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ z>0.1, 4<Q<9 ######16
data_current=data_from_f[33:38]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[16.,81.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.055,0.8]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"Zgt01_p_x_Siv.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### POSITIVE HADRONS h+ z>0.1, 1<Q<2 ######17
data_current=data_from_f[6:13]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass16.sivers.h+.z>1.dx',"SIDIS")
DataCurrentSiv.comment="COMPASS16 SSA-Sivers h+, z>0.1 (differential in x)"
DataCurrentSiv.reference="1609.07374"


proc_current=[1,1,1,12101]
proc_denominator=[1,1,1,2101]
s_current=2*160*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

Q2_current=[1.,4.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ z>0.1, 2<Q<2.5 ######18
data_current=data_from_f[17:23]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

Q2_current=[4.,6.25]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p) 
    
print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ z>0.1, 2.5<Q<4 ######19
data_current=data_from_f[27:34]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[6.25,16.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]


for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p) 
    
print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ z>0.1, 4<Q<9 ######20
data_current=data_from_f[38:44]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[16.,81.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.1,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p) 
    
print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"Zgt01_p_z_Siv.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### POSITIVE HADRONS h+ z>0.1, 1<Q<2 ######21
data_current=data_from_f[6:11]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass16.sivers.h+.z>1.dz',"SIDIS")
DataCurrentSiv.comment="COMPASS16 SSA-Sivers h+, z>0.1 (differential in z)"
DataCurrentSiv.reference="1609.07374"


proc_current=[1,1,1,12101]
proc_denominator=[1,1,1,2101]
s_current=2*160*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

Q2_current=[1.,4.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]     
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.003,0.145]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]]   
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ z>0.1, 2<Q<2.5 ######22
data_current=data_from_f[15:20]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

Q2_current=[4.,6.25]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.014,0.215]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]] 
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ z>0.1, 2.5<Q<4 ######23
data_current=data_from_f[24:29]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
   
Q2_current=[6.25,16.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]      
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.022,0.550]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]]   
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ z>0.1, 4<Q<9 ######24
data_current=data_from_f[33:38]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
   
Q2_current=[16.,81.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]      
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.055,0.8]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]]   
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"Zgt02_p_pt_Siv.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### POSITIVE HADRONS h+ z>0.2, 1<Q<2 ######25
data_current=data_from_f[6:11]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")


print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass16.sivers.h+.z>2.dpt',"SIDIS")
DataCurrentSiv.comment="COMPASS16 SSA-Sivers h+, z>0.2 (differential in pt)"
DataCurrentSiv.reference="1609.07374"


proc_current=[1,1,1,12101]
proc_denominator=[1,1,1,2101]
s_current=2*160*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

Q2_current=[1.,4.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.003,0.145]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.2,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ z>0.2, 2<Q<2.5 ######26
data_current=data_from_f[15:20]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

Q2_current=[4.,6.25]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.014,0.215]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.2,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ z>0.2, 2.5<Q<4 ######27
data_current=data_from_f[24:29]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[6.25,16.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.022,0.550]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.2,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ z>0.2, 4<Q<9 ######28
data_current=data_from_f[33:38]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[16.,81.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[data_current[i][0],data_current[i][1]]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.055,0.8]
    p["<z>"]=data_current[i][4]
    p["z"]=[0.2,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"Zgt02_p_x_Siv.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### POSITIVE HADRONS h+ z>0.2, 1<Q<2 ######29
data_current=data_from_f[6:13]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass16.sivers.h+.z>2.dx',"SIDIS")
DataCurrentSiv.comment="COMPASS16 SSA-Sivers h+, z>0.2 (differential in x)"
DataCurrentSiv.reference="1609.07374"


proc_current=[1,1,1,12101]
proc_denominator=[1,1,1,2101]
s_current=2*160*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

Q2_current=[1.,4.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.2,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ z>0.2, 2<Q<2.5 ######30
data_current=data_from_f[17:23]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

Q2_current=[4.,6.25]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.2,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p) 
    
print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ z>0.2, 2.5<Q<4 ######31
data_current=data_from_f[27:34]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[6.25,16.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]


for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.2,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p) 
    
print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ z>0.2, 4<Q<9 ######32
data_current=data_from_f[38:44]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
    
Q2_current=[16.,81.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=[data_current[i][0],data_current[i][1]]       
    p["<z>"]=data_current[i][4]
    p["z"]=[0.2,1.]
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p) 
    
print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###############################################################################
f = open(path_to_data+path_to_COMPASS+"Zgt02_p_z_Siv.dat")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### POSITIVE HADRONS h+ z>0.2, 1<Q<2 ######33
data_current=data_from_f[6:11]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass16.sivers.h+.z>2.dz',"SIDIS")
DataCurrentSiv.comment="COMPASS16 SSA-Sivers h+, z>0.2 (differential in z)"
DataCurrentSiv.reference="1609.07374"


proc_current=[1,1,1,12101]
proc_denominator=[1,1,1,2101]
s_current=2*160*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts

Q2_current=[1.,4.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]     
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.003,0.145]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]]   
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ z>0.2, 2<Q<2.5 ######34
data_current=data_from_f[15:20]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

Q2_current=[4.,6.25]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]       
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])  ##[0.014,0.215]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]] 
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ z>0.2, 2.5<Q<4 ######35
data_current=data_from_f[24:29]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
   
Q2_current=[6.25,16.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]      
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]]   
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")

#%%
###### POSITIVE HADRONS h+ z>0.2, 4<Q<9 ######36
data_current=data_from_f[33:38]


for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(";"," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")
   
Q2_current=[16.,81.]
Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]

for i in range(len(data_current)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_current[i][5]    
    p["pT"]=[0.1,10.]      
    p["<x>"]=data_current[i][2]    
    p["x"]=xbounds(Q2_current[0],Q2_current[1])   ##[0.055,0.8]
    p["<z>"]=data_current[i][4]
    p["z"]=[data_current[i][0],data_current[i][1]]   
    p["<Q>"]=numpy.sqrt(data_current[i][7])
    p["Q"]=Q_current
    p["xSec"]=data_current[i][8]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters    
    p["thFactor"]=1.         ### tobe updated
    #p["uncorrErr"].append(data_current[i][9])
    p["uncorrErr"].append(data_current[i][10])
    p["weightProcess"]=proc_denominator
    
    DataCurrentSiv.AddPoint(p)    

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
