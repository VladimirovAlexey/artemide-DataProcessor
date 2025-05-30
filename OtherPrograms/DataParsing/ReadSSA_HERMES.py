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
path_to_HERMES="/HERMES-SSA/"
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
    yMax=0.95
    sM2=2*27.6*0.938
    
    if xMax<1:
        return [numpy.sqrt(numpy.max([Q2min,xMin*yMin*sM2, xMin/(1-xMin)*WM2min])),
                numpy.sqrt(numpy.min([Q2max,xMax*yMax*sM2, xMax/(1-xMax)*WM2max]))]
    else:
        return [numpy.sqrt(numpy.max([Q2min,xMin*yMin*sM2, xMin/(1-xMin)*WM2min])),
                numpy.sqrt(numpy.min([Q2max,xMax*yMax*sM2, 1000*WM2max]))]

#%%

### bins are presented by Gunar Schnell
### There should be taken the central VALUE!!! the bin integration is included in sys-error
xBin=[0.023,0.045,0.067,0.086,0.113,0.160,0.220,0.400]
zBin=[0.2,0.27,0.34,0.41,0.49,0.56,0.63,0.70]
ptBin=[0.00001,0.17,0.25,0.33,0.41,0.58,0.80,2.00]
### Scale uncertanty
scaleUncertanty=0.073

#%%
###############################################################################
f = open(path_to_data+path_to_HERMES+"ssa1.cgi")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
###### POSITIVE PIONS PI+ ######
######    d Z             ######
data_current=data_from_f[16:23]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.pi+.Qint.dz',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers pi+ (integrated in Q, differential in z). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"

proc_current=[1,1,1,12001]
proc_denominator=[1,1,1,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[ptBin[0],ptBin[7]]
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[xBin[0],xBin[7]]
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][4]
    p1["z"]=zBin[binN-1:binN+1]        
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

###### POSITIVE PIONS PI+ ######
######    d X             ######
data_current=data_from_f[28:35]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.pi+.Qint.dx',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers pi+ (integrated in Q, differential in x). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,1,12001]
proc_denominator=[1,1,1,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[ptBin[0],ptBin[7]]
    p1["<x>"]=data_current[i][2]    
    p1["x"]=xBin[binN-1:binN+1]        
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][4]
    p1["z"]=[zBin[0],zBin[7]]       
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

###### POSITIVE PIONS PI+ ######
######    d pT            ######
data_current=data_from_f[40:47]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.pi+.Qint.dpt',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers pi+ (integrated in Q, differential in pt). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,1,12001]
proc_denominator=[1,1,1,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=ptBin[binN-1:binN+1]        
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[xBin[0],xBin[7]]
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][4]
    p1["z"]=[zBin[0],zBin[7]]       
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
###### POSITIVE PIONS PI0 ######
######    d Z             ######
data_current=data_from_f[52:59]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.pi0.Qint.dz',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers pi0 (integrated in Q, differential in z). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,3,12001]
proc_denominator=[1,1,3,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[ptBin[0],ptBin[7]]
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[xBin[0],xBin[7]]
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][4]
    p1["z"]=zBin[binN-1:binN+1]        
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

###### POSITIVE PIONS PI0 ######
######    d X             ######
data_current=data_from_f[64:71]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.pi0.Qint.dx',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers pi0 (integrated in Q, differential in x). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,3,12001]
proc_denominator=[1,1,3,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[ptBin[0],ptBin[7]]
    p1["<x>"]=data_current[i][2]    
    p1["x"]=xBin[binN-1:binN+1]        
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][4]
    p1["z"]=[zBin[0],zBin[7]]       
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
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%

###### POSITIVE PIONS PI0 ######
######    d pT            ######
data_current=data_from_f[76:83]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.pi0.Qint.dpt',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers pi0 (integrated in Q, differential in pt). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,3,12001]
proc_denominator=[1,1,3,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=ptBin[binN-1:binN+1]        
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[xBin[0],xBin[7]]
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][4]
    p1["z"]=[zBin[0],zBin[7]]       
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
#DataCurrentCol.SaveToCSV(path_to_save+DataCurrentCol.name+".csv")

#%%
###### POSITIVE PIONS PI- ######
######    d Z             ######
data_current=data_from_f[88:95]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.pi-.Qint.dz',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers pi- (integrated in Q, differential in z). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,-1,12001]
proc_denominator=[1,1,-1,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[ptBin[0],ptBin[7]]
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[xBin[0],xBin[7]]
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][4]
    p1["z"]=zBin[binN-1:binN+1]        
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

###### POSITIVE PIONS PI- ######
######    d X             ######
data_current=data_from_f[100:107]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.pi-.Qint.dx',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers pi- (integrated in Q, differential in x). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,-1,12001]
proc_denominator=[1,1,-1,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[ptBin[0],ptBin[7]]
    p1["<x>"]=data_current[i][2]    
    p1["x"]=xBin[binN-1:binN+1]        
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][4]
    p1["z"]=[zBin[0],zBin[7]]       
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

###### POSITIVE PIONS PI- ######
######    d pT            ######
data_current=data_from_f[112:119]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.pi-.Qint.dpt',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers pi- (integrated in Q, differential in pt). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,-1,12001]
proc_denominator=[1,1,-1,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=ptBin[binN-1:binN+1]        
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[xBin[0],xBin[7]]
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][4]
    p1["z"]=[zBin[0],zBin[7]]       
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
###### POSITIVE PIONS K+  ######
######    d Z             ######
data_current=data_from_f[124:131]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.k+.Qint.dz',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers k+ (integrated in Q, differential in z). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,2,12001]
proc_denominator=[1,1,2,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[ptBin[0],ptBin[7]]
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[xBin[0],xBin[7]]
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][4]
    p1["z"]=zBin[binN-1:binN+1]        
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

###### POSITIVE PIONS K+  ######
######    d X             ######
data_current=data_from_f[136:143]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.k+.Qint.dx',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers k+ (integrated in Q, differential in x). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,2,12001]
proc_denominator=[1,1,2,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[ptBin[0],ptBin[7]]
    p1["<x>"]=data_current[i][2]    
    p1["x"]=xBin[binN-1:binN+1]        
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][4]
    p1["z"]=[zBin[0],zBin[7]]       
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
######    d pT            ######
data_current=data_from_f[148:155]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.k+.Qint.dpt',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers k+ (integrated in Q, differential in pt). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,2,12001]
proc_denominator=[1,1,2,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=ptBin[binN-1:binN+1]        
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[xBin[0],xBin[7]]
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][4]
    p1["z"]=[zBin[0],zBin[7]]       
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
###### POSITIVE PIONS K-  ######
######    d Z             ######
data_current=data_from_f[160:167]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.k-.Qint.dz',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers k- (integrated in Q, differential in z). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,-2,12001]
proc_denominator=[1,1,-2,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[ptBin[0],ptBin[7]]
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[xBin[0],xBin[7]]
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][4]
    p1["z"]=zBin[binN-1:binN+1]        
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

###### POSITIVE PIONS K-  ######
######    d X             ######
data_current=data_from_f[172:179]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.k-.Qint.dx',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers k- (integrated in Q, differential in x). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,-2,12001]
proc_denominator=[1,1,-2,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[ptBin[0],ptBin[7]]
    p1["<x>"]=data_current[i][2]    
    p1["x"]=xBin[binN-1:binN+1]        
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][4]
    p1["z"]=[zBin[0],zBin[7]]       
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

###### POSITIVE PIONS K- ######
######    d pT            ######
data_current=data_from_f[184:191]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.k-.Qint.dpt',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers k- (integrated in Q, differential in pt). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,-2,12001]
proc_denominator=[1,1,-2,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=ptBin[binN-1:binN+1]        
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[xBin[0],xBin[7]]
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][4]
    p1["z"]=[zBin[0],zBin[7]]       
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
##################################################################
########   READING FILE SSA2   ###################################
##################################################################
f = open(path_to_data+path_to_HERMES+"ssa2.cgi")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%

###### POSITIVE PIONS PI+ ######
######    d Z   (Q<2)     ######
data_current=data_from_f[16:23]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.pi+.Q<2.dz',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers pi+ (Q<2, differential in z). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,1,12001]
proc_denominator=[1,1,1,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[ptBin[0],ptBin[7]]
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[xBin[0],xBin[7]]
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=[1.,2.]    
    p1["<z>"]=data_current[i][4]
    p1["z"]=zBin[binN-1:binN+1]        
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

###### POSITIVE PIONS PI+ ######
######    d pT  (Q<2)     ######
data_current=data_from_f[28:35]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.pi+.Q<2.dpt',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers pi+ (Q<2, differential in pt). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,1,12001]
proc_denominator=[1,1,1,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=ptBin[binN-1:binN+1]        
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[xBin[0],xBin[7]]
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=[1.,2.]    
    p1["<z>"]=data_current[i][4]
    p1["z"]=[zBin[0],zBin[7]]       
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

###### POSITIVE PIONS PI+ ######
######    d Z   (Q>2)     ######
data_current=data_from_f[40:47]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.pi+.Q>2.dz',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers pi+ (Q>2, differential in z). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,1,12001]
proc_denominator=[1,1,1,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[ptBin[0],ptBin[7]]
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[xBin[0],xBin[7]]
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=[2.,2.*27.6*0.938*0.95*p1["x"][1]]    
    p1["<z>"]=data_current[i][4]
    p1["z"]=zBin[binN-1:binN+1]        
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

###### POSITIVE PIONS PI+ ######
######    d pT  (Q>2)     ######
data_current=data_from_f[52:59]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.pi+.Q>2.dpt',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers pi+ (Q>2, differential in pt). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,1,12001]
proc_denominator=[1,1,1,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=ptBin[binN-1:binN+1]        
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[xBin[0],xBin[7]]
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=[2.,2.*27.6*0.938*0.95*p1["x"][1]]    
    p1["<z>"]=data_current[i][4]
    p1["z"]=[zBin[0],zBin[7]]       
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
######    d Z   (Q<2)     ######
data_current=data_from_f[64:71]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.k+.Q<2.dz',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers k+ (Q<2, differential in z). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,2,12001]
proc_denominator=[1,1,2,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[ptBin[0],ptBin[7]]
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[xBin[0],xBin[7]]
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=[1.,2.]    
    p1["<z>"]=data_current[i][4]
    p1["z"]=zBin[binN-1:binN+1]        
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

###### POSITIVE PIONS k+ ######
######    d pT  (Q<2)     ######
data_current=data_from_f[76:83]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.k+.Q<2.dpt',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers k+ (Q<2, differential in pt). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,2,12001]
proc_denominator=[1,1,2,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=ptBin[binN-1:binN+1]        
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[xBin[0],xBin[7]]
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=[1.,2.]    
    p1["<z>"]=data_current[i][4]
    p1["z"]=[zBin[0],zBin[7]]       
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

###### POSITIVE PIONS k+ ######
######    d Z   (Q>2)     ######
data_current=data_from_f[88:95]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.k+.Q>2.dz',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers k+ (Q>2, differential in z). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,2,12001]
proc_denominator=[1,1,2,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=[ptBin[0],ptBin[7]]
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[xBin[0],xBin[7]]
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=[2.,2.*27.6*0.938*0.95*p1["x"][1]]    
    p1["<z>"]=data_current[i][4]
    p1["z"]=zBin[binN-1:binN+1]        
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

###### POSITIVE PIONS k+ ######
######    d pT  (Q>2)     ######
data_current=data_from_f[100:107]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()    
    data_current[i]=[float(j) for j in data_current[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.k+.Q>2.dpt',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers k+ (Q>2, differential in pt). The data MUST be evaluated at a point"
DataCurrent.reference="0906.3918"


proc_current=[1,1,2,12001]
proc_denominator=[1,1,2,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(data_current[0][0])))
    #print DataCurrent.name+'.'+str(i)
    binN=int(data_current[i][0])
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][5]    
    p1["pT"]=ptBin[binN-1:binN+1]        
    p1["<x>"]=data_current[i][2]    
    p1["x"]=[xBin[0],xBin[7]]
    p1["<Q>"]=numpy.sqrt(data_current[i][1])
    p1["Q"]=[2.,2.*27.6*0.938*0.95*p1["x"][1]]    
    p1["<z>"]=data_current[i][4]
    p1["z"]=[zBin[0],zBin[7]]       
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
##################################################################
########   HERMES 3D PI+  ############################################
##################################################################
f = open(path_to_data+path_to_HERMES+"SFA/hermesTMDs_data__pip_SFA_3D.txt")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
data_current=data_from_f[86:150]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()
    data_current[i][0]=[float(j) for j in data_current[i][0].split("<x<")]
    data_current[i][1]=[float(j) for j in data_current[i][1].split("<z<")]
    data_current[i][2]=[float(j) for j in data_current[i][2].split("<Pt<")]
    data_current[i][3:12]=[float(j) for j in data_current[i][3:12]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.pi+.3d',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers pi+ (3d-data). The data MUST be evaluated at a point"
DataCurrent.reference="2007.07755"


proc_current=[1,1,1,12001]
proc_denominator=[1,1,1,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(i)))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][7]    
    p1["pT"]=data_current[i][2]
    p1["<x>"]=data_current[i][4]    
    p1["x"]=data_current[i][0]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][6]
    p1["z"]=data_current[i][1]        
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["uncorrErr"].append(data_current[i][11])
    p1["weightProcess"]=proc_denominator
    #
    DataCurrent.AddPoint(p1)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")


#%%
##################################################################
########   HERMES 3D PI-  ############################################
##################################################################
f = open(path_to_data+path_to_HERMES+"SFA/hermesTMDs_data__pim_SFA_3D.txt")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
data_current=data_from_f[86:150]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()
    data_current[i][0]=[float(j) for j in data_current[i][0].split("<x<")]
    data_current[i][1]=[float(j) for j in data_current[i][1].split("<z<")]
    data_current[i][2]=[float(j) for j in data_current[i][2].split("<Pt<")]
    data_current[i][3:12]=[float(j) for j in data_current[i][3:12]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.pi-.3d',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers pi- (3d-data). The data MUST be evaluated at a point"
DataCurrent.reference="2007.07755"


proc_current=[1,1,-1,12001]
proc_denominator=[1,1,-1,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(i)))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][7]    
    p1["pT"]=data_current[i][2]
    p1["<x>"]=data_current[i][4]    
    p1["x"]=data_current[i][0]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][6]
    p1["z"]=data_current[i][1]        
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["uncorrErr"].append(data_current[i][11])
    p1["weightProcess"]=proc_denominator
    #
    DataCurrent.AddPoint(p1)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")


#%%
##################################################################
########   HERMES 3D K+  ############################################
##################################################################
f = open(path_to_data+path_to_HERMES+"SFA/hermesTMDs_data__kp_SFA_3D.txt")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
data_current=data_from_f[86:150]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()
    data_current[i][0]=[float(j) for j in data_current[i][0].split("<x<")]
    data_current[i][1]=[float(j) for j in data_current[i][1].split("<z<")]
    data_current[i][2]=[float(j) for j in data_current[i][2].split("<Pt<")]
    data_current[i][3:12]=[float(j) for j in data_current[i][3:12]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.k+.3d',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers k+ (3d-data). The data MUST be evaluated at a point"
DataCurrent.reference="2007.07755"


proc_current=[1,1,2,12001]
proc_denominator=[1,1,2,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(i)))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][7]    
    p1["pT"]=data_current[i][2]
    p1["<x>"]=data_current[i][4]    
    p1["x"]=data_current[i][0]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][6]
    p1["z"]=data_current[i][1]        
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["uncorrErr"].append(data_current[i][11])
    p1["weightProcess"]=proc_denominator
    #
    DataCurrent.AddPoint(p1)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
##################################################################
########   HERMES 3D K-  ############################################
##################################################################
f = open(path_to_data+path_to_HERMES+"SFA/hermesTMDs_data__km_SFA_3D.txt")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
data_current=data_from_f[86:150]

for i in range(len(data_current)):
    data_current[i]=data_current[i].split()
    data_current[i][0]=[float(j) for j in data_current[i][0].split("<x<")]
    data_current[i][1]=[float(j) for j in data_current[i][1].split("<z<")]
    data_current[i][2]=[float(j) for j in data_current[i][2].split("<Pt<")]
    data_current[i][3:12]=[float(j) for j in data_current[i][3:12]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.sivers.k-.3d',"SIDIS")
DataCurrent.comment="HERMESS SSA-Sivers k- (3d-data). The data MUST be evaluated at a point"
DataCurrent.reference="2007.07755"


proc_current=[1,1,-2,12001]
proc_denominator=[1,1,-2,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
DataCurrent.normErr.append(scaleUncertanty)

for i in range(len(data_current)):
    # makeup a point
    p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(i)))
    #print DataCurrent.name+'.'+str(i)
    p1["process"]=proc_current
    p1["s"]=s_current
    p1["<pT>"]=data_current[i][7]    
    p1["pT"]=data_current[i][2]
    p1["<x>"]=data_current[i][4]    
    p1["x"]=data_current[i][0]
    p1["<Q>"]=numpy.sqrt(data_current[i][3])
    p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
    p1["<z>"]=data_current[i][6]
    p1["z"]=data_current[i][1]        
    p1["xSec"]=data_current[i][9]
    p1["M_target"]=M_proton
    p1["M_product"]=m_pion
    p1["includeCuts"]=includeCuts
    p1["cutParams"]=cutParameters    
    p1["thFactor"]=1.         ### tobe updated
    p1["uncorrErr"].append(data_current[i][10])
    p1["uncorrErr"].append(data_current[i][11])
    p1["weightProcess"]=proc_denominator
    #
    DataCurrent.AddPoint(p1)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")