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
path_to_COMPASS="/COMPASS/0802.2160/"
path_to_save="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/DataLib/Sivers/"

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
print("compass.proton.SSA-Sivers pi+ file ...")
f = open(path_to_data+path_to_COMPASS+"/durham_0802.2160.txt")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")



#%%
###### POSITIVE PIONS PI+ ######
del data_from_f[0:13]
data_current=data_from_f[0:9]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.pi+.dx',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers pi+ (differential in x)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.pi+.dx',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins pi+ (differential in x)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12011]
proc_denominator=[1,1,2011]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]    
    p1["pT"]=[0.1,10.]        
    p1["<x>"]=data_current[i][0]    
    p1["x"]=[data_current[i][1],data_current[i][2]]
    p1["<z>"]=data_current[i][5]
    p1["z"]=[0.2,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]
    #p["pT"]=[ptBINS[i],ptBINS[i+1]]    ### tobe updated
    p2["pT"]=[0.1,10.]        
    p2["<x>"]=data_current[i][0]    
    p2["x"]=[data_current[i][1],data_current[i][2]]
    p2["<z>"]=data_current[i][5]
    p2["z"]=[0.2,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    p2["weightProcess"]=proc_denominator
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")


#%%
data_current=data_from_f[11:20]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.pi+.dpt',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers pi+ (differential in pt)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.pi+.dpt',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins pi+ (differential in pt)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12011]
proc_denominator=[1,1,2011]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][0]       
    p1["pT"]=[data_current[i][1],data_current[i][2]]     
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][6]
    p1["z"]=[0.2,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][0]      
    p2["pT"]=[data_current[i][1],data_current[i][2]]    
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][6]
    p2["z"]=[0.2,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    p2["weightProcess"]=proc_denominator
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[22:30]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.pi+.dz',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers pi+ (differential in z)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.pi+.dz',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins pi+ (differential in z)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12011]
proc_denominator=[1,1,2011]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]  
    p1["pT"]=[0.1,10.]             
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][0]
    p1["z"]=[data_current[i][1],data_current[i][2]]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]       
    p2["pT"]=[0.1,10.]   
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][0]
    p2["z"]=[data_current[i][1],data_current[i][2]]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    p2["weightProcess"]=proc_denominator
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")


#%%
###### NEGATIVE PIONS PI- ######
data_current=data_from_f[34:43]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.pi-.dx',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers pi- (differential in x)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.pi-.dx',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins pi- (differential in x)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12031]
proc_denominator=[1,1,2031]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]    
    p1["pT"]=[0.1,10.]        
    p1["<x>"]=data_current[i][0]    
    p1["x"]=[data_current[i][1],data_current[i][2]]
    p1["<z>"]=data_current[i][5]
    p1["z"]=[0.2,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]
    p2["pT"]=[0.1,10.]    
    p2["<x>"]=data_current[i][0]    
    p2["x"]=[data_current[i][1],data_current[i][2]]
    p2["<z>"]=data_current[i][5]
    p2["z"]=[0.2,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    p2["weightProcess"]=proc_denominator
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[45:54]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.pi-.dpt',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers pi- (differential in pt)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.pi-.dpt',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins pi- (differential in pt)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12031]
proc_denominator=[1,1,2031]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][0]       
    p1["pT"]=[data_current[i][1],data_current[i][2]]     
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][6]
    p1["z"]=[0.2,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][0]      
    p2["pT"]=[data_current[i][1],data_current[i][2]] 
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][6]
    p2["z"]=[0.2,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    p2["weightProcess"]=proc_denominator
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[56:64]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.pi-.dz',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers pi- (differential in z)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.pi-.dz',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins pi- (differential in z)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12031]
proc_denominator=[1,1,2031]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]  
    p1["pT"]=[0.1,10.]             
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][0]
    p1["z"]=[data_current[i][1],data_current[i][2]]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]       
    p2["pT"]=[0.1,10.]  
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][0]
    p2["z"]=[data_current[i][1],data_current[i][2]]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    p2["weightProcess"]=proc_denominator
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

################################################
################################################
################################################
###############################################
#%%
###### POSITIVE KAONS K+ ######
data_current=data_from_f[68:77]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.k+.dx',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers k+ (differential in x)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.k+.dx',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins k+ (differential in x)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12012]
proc_denominator=[1,1,2012]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]    
    p1["pT"]=[0.1,10.]        
    p1["<x>"]=data_current[i][0]    
    p1["x"]=[data_current[i][1],data_current[i][2]]
    p1["<z>"]=data_current[i][5]
    p1["z"]=[0.2,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]
    p2["pT"]=[0.1,10.]    
    p2["<x>"]=data_current[i][0]    
    p2["x"]=[data_current[i][1],data_current[i][2]]
    p2["<z>"]=data_current[i][5]
    p2["z"]=[0.2,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    p2["weightProcess"]=proc_denominator
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[79:88]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.k+.dpt',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers k+ (differential in pt)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.k+.dpt',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins k+ (differential in pt)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12012]
proc_denominator=[1,1,2012]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][0]       
    p1["pT"]=[data_current[i][1],data_current[i][2]]     
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][6]
    p1["z"]=[0.2,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][0]      
    p2["pT"]=[data_current[i][1],data_current[i][2]]   
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][6]
    p2["z"]=[0.2,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    p2["weightProcess"]=proc_denominator
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[90:98]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.k+.dz',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers k+ (differential in z)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.k+.dz',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins k+ (differential in z)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12012]
proc_denominator=[1,1,2012]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]  
    p1["pT"]=[0.1,10.]             
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][0]
    p1["z"]=[data_current[i][1],data_current[i][2]]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]       
    p2["pT"]=[0.1,10.]   
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][0]
    p2["z"]=[data_current[i][1],data_current[i][2]]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    p2["weightProcess"]=proc_denominator
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

###################
###############
########
###
#
#%%
###### NEGATIVE KAONS K- ######
data_current=data_from_f[102:111]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.k-.dx',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers k- (differential in x)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.k-.dx',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins k- (differential in x)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12032]
proc_denominator=[1,1,2032]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]    
    p1["pT"]=[0.1,10.]        
    p1["<x>"]=data_current[i][0]    
    p1["x"]=[data_current[i][1],data_current[i][2]]
    p1["<z>"]=data_current[i][5]
    p1["z"]=[0.2,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]
    p2["pT"]=[0.1,10.]    
    p2["<x>"]=data_current[i][0]    
    p2["x"]=[data_current[i][1],data_current[i][2]]
    p2["<z>"]=data_current[i][5]
    p2["z"]=[0.2,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    p2["weightProcess"]=proc_denominator
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[113:122]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.k-.dpt',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers k- (differential in pt)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.k-.dpt',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins k- (differential in pt)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12032]
proc_denominator=[1,1,2032]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][0]       
    p1["pT"]=[data_current[i][1],data_current[i][2]]     
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][6]
    p1["z"]=[0.2,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][0]      
    p2["pT"]=[data_current[i][1],data_current[i][2]] 
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][6]
    p2["z"]=[0.2,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    p2["weightProcess"]=proc_denominator
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[124:132]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.k-.dz',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers k- (differential in z)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.k-.dz',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins k- (differential in z)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12032]
proc_denominator=[1,1,2032]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]  
    p1["pT"]=[0.1,10.]
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][0]
    p1["z"]=[data_current[i][1],data_current[i][2]]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]       
    p2["pT"]=[0.1,10.]   
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][0]
    p2["z"]=[data_current[i][1],data_current[i][2]]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    p2["weightProcess"]=proc_denominator
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
###### NEUTRAL KAONS K0 ######
data_current=data_from_f[136:141]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.k0.dx',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers k0 (differential in x)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.k0.dx',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins k0 (differential in x)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12013]
proc_denominator=[1,1,2013]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]    
    p1["pT"]=[0.1,10.]        
    p1["<x>"]=data_current[i][0]    
    p1["x"]=[data_current[i][1],data_current[i][2]]
    p1["<z>"]=data_current[i][5]
    p1["z"]=[0.2,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]
    p2["pT"]=[0.1,10.]    
    p2["<x>"]=data_current[i][0]    
    p2["x"]=[data_current[i][1],data_current[i][2]]
    p2["<z>"]=data_current[i][5]
    p2["z"]=[0.2,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    p2["weightProcess"]=proc_denominator
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[143:148]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.k0.dpt',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers k0 (differential in pt)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.k0.dpt',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins k0 (differential in pt)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12013]
proc_denominator=[1,1,2013]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][0]       
    p1["pT"]=[data_current[i][1],data_current[i][2]] 
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][6]
    p1["z"]=[0.2,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator    
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][0]      
    p2["pT"]=[data_current[i][1],data_current[i][2]]    
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][6]
    p2["z"]=[0.2,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    p2["weightProcess"]=proc_denominator
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[150:156]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.k0.dz',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers k0 (differential in z)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.k0.dz',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins k0 (differential in z)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12013]
proc_denominator=[1,1,2013]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]  
    p1["pT"]=[0.1,10.]             
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][0]
    p1["z"]=[data_current[i][1],data_current[i][2]]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["weightProcess"]=proc_denominator
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]       
    p2["pT"]=[0.1,10.]   
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][0]
    p2["z"]=[data_current[i][1],data_current[i][2]]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    p2["weightProcess"]=proc_denominator
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#
##
####
######
########
##########
############    LEADING HADRONS 
##########
########
######
####
##
#
#%%
###### POSITIVE PIONS PI+ ######
data_current=data_from_f[162:171]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.pi+leading.dx',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers pi+ leading (differential in x)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.pi+leading.dx',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins pi+ leading (differential in x)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12001]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]    
    p1["pT"]=[0.1,10.]
    p1["<x>"]=data_current[i][0]    
    p1["x"]=[data_current[i][1],data_current[i][2]]
    p1["<z>"]=data_current[i][5]
    p1["z"]=[0.25,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]
    p2["pT"]=[0.1,10.]    
    p2["<x>"]=data_current[i][0]    
    p2["x"]=[data_current[i][1],data_current[i][2]]
    p2["<z>"]=data_current[i][5]
    p2["z"]=[0.25,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[173:182]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.pi+leading.dpt',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers pi+ leading (differential in pt)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.pi+leading.dpt',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins pi+ leading (differential in pt)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12001]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][0]       
    p1["pT"]=[data_current[i][1],data_current[i][2]] 
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][6]
    p1["z"]=[0.25,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][0]      
    p2["pT"]=[data_current[i][1],data_current[i][2]]    
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][6]
    p2["z"]=[0.25,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[184:191]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.pi+leading.dz',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers pi+ leading (differential in z)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.pi+leading.dz',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins pi+ leading (differential in z)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12001]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]  
    p1["pT"]=[0.1,10.]
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][0]
    p1["z"]=[data_current[i][1],data_current[i][2]]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]       
    p2["pT"]=[0.1,10.]   
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][0]
    p2["z"]=[data_current[i][1],data_current[i][2]]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")


#%%
###### NEGATIVE PIONS PI- ######
data_current=data_from_f[195:204]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.pi-leading.dx',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers pi- leading (differential in x)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.pi-leading.dx',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins pi- leading (differential in x)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12001]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]    
    p1["pT"]=[0.1,10.]    
    p1["<x>"]=data_current[i][0]    
    p1["x"]=[data_current[i][1],data_current[i][2]]
    p1["<z>"]=data_current[i][5]
    p1["z"]=[0.25,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]
    #p["pT"]=[ptBINS[i],ptBINS[i+1]]    ### tobe updated
    p2["pT"]=[0.1,10.]    
    p2["<x>"]=data_current[i][0]    
    p2["x"]=[data_current[i][1],data_current[i][2]]
    p2["<z>"]=data_current[i][5]
    p2["z"]=[0.25,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[206:215]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.pi-leading.dpt',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers pi- leading (differential in pt)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.pi-leading.dpt',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins pi- leading (differential in pt)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12001]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][0]       
    p1["pT"]=[data_current[i][1],data_current[i][2]] 
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][6]
    p1["z"]=[0.25,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][0]      
    p2["pT"]=[data_current[i][1],data_current[i][2]] 
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][6]
    p2["z"]=[0.25,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[217:224]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.pi-leading.dz',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers pi- leading (differential in z)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.pi-leading.dz',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins pi- leading (differential in z)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12001]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]  
    p1["pT"]=[0.1,10.]     
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][0]
    p1["z"]=[data_current[i][1],data_current[i][2]]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]       
    p2["pT"]=[0.1,10.]   
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][0]
    p2["z"]=[data_current[i][1],data_current[i][2]]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

################################################
################################################
################################################
#%%
###### POSITIVE KAONS K+ ######
data_current=data_from_f[228:237]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.k+leading.dx',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers k+ leading (differential in x)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.k+leading.dx',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins k+ leading (differential in x)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12001]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]    
    p1["pT"]=[0.1,10.]    
    p1["<x>"]=data_current[i][0]    
    p1["x"]=[data_current[i][1],data_current[i][2]]
    p1["<z>"]=data_current[i][5]
    p1["z"]=[0.25,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]
    p2["pT"]=[0.1,10.]    
    p2["<x>"]=data_current[i][0]    
    p2["x"]=[data_current[i][1],data_current[i][2]]
    p2["<z>"]=data_current[i][5]
    p2["z"]=[0.25,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[239:248]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.k+leading.dpt',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers k+ leading (differential in pt)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.k+leading.dpt',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins k+ leading (differential in pt)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12001]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][0]       
    p1["pT"]=[data_current[i][1],data_current[i][2]] 
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][6]
    p1["z"]=[0.25,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][0]      
    p2["pT"]=[data_current[i][1],data_current[i][2]]    
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][6]
    p2["z"]=[0.25,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[250:257]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.k+leading.dz',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers k+ leading (differential in z)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.k+leading.dz',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins k+ leading (differential in z)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12001]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]  
    p1["pT"]=[0.1,10.]         
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][0]
    p1["z"]=[data_current[i][1],data_current[i][2]]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]       
    p2["pT"]=[0.1,10.]   
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][0]
    p2["z"]=[data_current[i][1],data_current[i][2]]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

###################
###############
########
###
#
#%%
###### NEGATIVE KAONS K- ######
data_current=data_from_f[261:270]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.k-leading.dx',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers k- leading (differential in x)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.k-leading.dx',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins k- leading (differential in x)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12001]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]    
    p1["pT"]=[0.1,10.]  
    p1["<x>"]=data_current[i][0]    
    p1["x"]=[data_current[i][1],data_current[i][2]]
    p1["<z>"]=data_current[i][5]
    p1["z"]=[0.25,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]
    p2["pT"]=[0.1,10.]    
    p2["<x>"]=data_current[i][0]    
    p2["x"]=[data_current[i][1],data_current[i][2]]
    p2["<z>"]=data_current[i][5]
    p2["z"]=[0.25,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[272:281]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.k-leading.dpt',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers k- leading (differential in pt)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.k-leading.dpt',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins k- leading (differential in pt)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12001]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][0]       
    p1["pT"]=[data_current[i][1],data_current[i][2]] 
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][6]
    p1["z"]=[0.25,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][0]      
    p2["pT"]=[data_current[i][1],data_current[i][2]]   
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][6]
    p2["z"]=[0.25,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[283:290]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.k-leading.dz',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers k- leading (differential in z)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.k-leading.dz',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins k- leading (differential in z)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12001]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]  
    p1["pT"]=[0.1,10.]         
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][0]
    p1["z"]=[data_current[i][1],data_current[i][2]]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]       
    p2["pT"]=[0.1,10.]   
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][0]
    p2["z"]=[data_current[i][1],data_current[i][2]]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
###### NEUTRAL KAONS K0 ######
data_current=data_from_f[294:299]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.k0leading.dx',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers k0 leading (differential in x)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.k0leading.dx',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins k0 leading (differential in x)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12001]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]    
    p1["pT"]=[0.1,10.] 
    p1["<x>"]=data_current[i][0]    
    p1["x"]=[data_current[i][1],data_current[i][2]]
    p1["<z>"]=data_current[i][5]
    p1["z"]=[0.25,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]
    #p["pT"]=[ptBINS[i],ptBINS[i+1]]    ### tobe updated
    p2["pT"]=[0.1,10.]    
    p2["<x>"]=data_current[i][0]    
    p2["x"]=[data_current[i][1],data_current[i][2]]
    p2["<z>"]=data_current[i][5]
    p2["z"]=[0.25,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[301:306]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.k0leading.dpt',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers k0 leading (differential in pt)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.k0leading.dpt',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins k0 leading (differential in pt)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12001]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][0]       
    p1["pT"]=[data_current[i][1],data_current[i][2]] 
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][6]
    p1["z"]=[0.25,1.]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][0]      
    p2["pT"]=[data_current[i][1],data_current[i][2]]    
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][6]
    p2["z"]=[0.25,1.]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
data_current=data_from_f[308:313]

for i in range(len(data_current)):
    data_current[i]=data_current[i].replace("[","").replace("]","").replace(","," ")
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrentSiv=DataProcessor.DataSet.DataSet('compass08.sivers.k0leading.dz',"SIDIS")
DataCurrentSiv.comment="COMPASS SSA-Sivers k0 leading (differential in z)"
DataCurrentSiv.reference="0802.2160"
DataCurrentCol=DataProcessor.DataSet.DataSet('compass08.collins.k0leading.dz',"SIDIS")
DataCurrentCol.comment="COMPASS SSA-Collins k0 leading (differential in z)"
DataCurrentCol.reference="0802.2160"


proc_current=[1,1,12001]
s_current=2.*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.] #y, W^2 cuts

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrentSiv.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][6]  
    p1["pT"]=[0.1,10.]         
    p1["<x>"]=data_current[i][5]    
    p1["x"]=[0.03,1.]
    p1["<z>"]=data_current[i][0]
    p1["z"]=[data_current[i][1],data_current[i][2]]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    #
    p2=DataProcessor.Point.CreateSIDISPoint(DataCurrentCol.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    p2["process"]=proc_current
    p2["s"]=s_current
    p2["<pT>"]=data_current[i][6]       
    p2["pT"]=[0.1,10.]   
    p2["<x>"]=data_current[i][5]    
    p2["x"]=[0.03,1.]
    p2["<z>"]=data_current[i][0]
    p2["z"]=[data_current[i][1],data_current[i][2]]
    p2["<Q>"]=numpy.sqrt(data_current[i][3])
    p2["Q"]=Qbounds(p2["x"][0],p2["x"][1])
    p2["xSec"]=data_current[i][7]
    p2["M_target"]=M_proton
    p2["M_product"]=m_pion
    p2["includeCuts"]=includeCuts
    p2["cutParams"]=cutParameters    
    p2["thFactor"]=1.         ### tobe updated
    p2["uncorrErr"].append(data_current[i][8])
    #
    DataCurrentSiv.AddPoint(p1)    
    DataCurrentCol.AddPoint(p2)  

print("Done.  ")

DataCurrentSiv.SaveToCSV(path_to_save+DataCurrentSiv.name+".csv")
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")