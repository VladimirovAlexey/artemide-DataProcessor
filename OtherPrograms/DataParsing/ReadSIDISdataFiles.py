#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 20:15:46 2019

Program that parse various SIDIS data files to ADP-frendly formal

@author: vla18041
"""

import os
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))

import sys
sys.path.append(ROOT_DIR)

import DataProcessor.Point
import DataProcessor.DataSet
import numpy

path_to_data=os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+'/DataLib/unpolSIDIS/Multiplicity_norm_byVMoos/'
path_to_save=ROOT_DIR+"/DataLib/unpolSIDIS/"


## this trigger add normalization uncertanty due to DIS normalization for multiplicities
addDISnormalizationUncertainty=True

M_proton=0.938
m_pion=0.139
m_kaon=0.494
#%%
###############################################################################
###########################hermes.proton.zxpt-3D pi+###########################
print("hermes.proton.zxpt-3D pi+ file ...")
f = open(path_to_data+"hermes.proton.zxpt-3D.no-vmsub.mults_piplus.datnew_fortran")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

del data_from_f[0:5]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split()
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")
#%%

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.p.no-vmsub.zxpt.pi+',"SIDIS")
DataCurrent.comment="HERMES proton-to-pi+ (zxpt-3D) no-vmsub. thFactor contains DIS normalization (by MSHT20nnlo)"
DataCurrent.reference="1212.5407"

proc_current=[1,1,1,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.85,10.,10000.]

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_from_f[i][12]
    p["pT"]=[data_from_f[i][13],data_from_f[i][14]]
    p["<Q>"]=numpy.sqrt(data_from_f[i][3])
    p["Q"]=[numpy.sqrt(data_from_f[i][4]),numpy.sqrt(data_from_f[i][5])]
    p["<x>"]=data_from_f[i][6]    
    p["x"]=[data_from_f[i][7],data_from_f[i][8]]
    p["<z>"]=data_from_f[i][9]
    p["z"]=[data_from_f[i][10],data_from_f[i][11]]
    p["xSec"]=data_from_f[i][0]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters
    if p["xSec"]<0.00000001:
        p["thFactor"]=1.
    else:
        p["thFactor"]=1/(p["pT"][1]-p["pT"][0])/(p["z"][1]-p["z"][0])/data_from_f[i][15]#devide by bin size multiply by DIS xSec
         
    p["uncorrErr"].append(data_from_f[i][1])
    p["uncorrErr"].append(data_from_f[i][2])
    if(addDISnormalizationUncertainty):
        if(data_from_f[i][0]==0):
            pass#p.corrErrors.append(0)
        else:
            p["corrErr"].append((data_from_f[i][16]+data_from_f[i][17])/2/data_from_f[i][15]*data_from_f[i][0])
    #
    DataCurrent.AddPoint(p)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%
###############################################################################
###########################hermes.proton.zxpt-3D pi-###########################
print("hermes.proton.zxpt-3D pi- file ..." )
f = open(path_to_data+"hermes.proton.zxpt-3D.no-vmsub.mults_piminus.datnew_fortran")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Convert to numbers ...")

del data_from_f[0:5]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split()
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")


print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.p.no-vmsub.zxpt.pi-',"SIDIS")
DataCurrent.comment="HERMES proton-to-pi- (zxpt-3D) no-vmsub. thFactor contains DIS normalization (by MSHT20nnlo)"
DataCurrent.reference="1212.5407"

proc_current=[1,1,-1,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.85,10.,10000.]

if(addDISnormalizationUncertainty):
        print('Correction for no-vmsub is not computed')

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_from_f[i][12]
    p["pT"]=[data_from_f[i][13],data_from_f[i][14]]
    p["<Q>"]=numpy.sqrt(data_from_f[i][3])
    p["Q"]=[numpy.sqrt(data_from_f[i][4]),numpy.sqrt(data_from_f[i][5])]
    p["<x>"]=data_from_f[i][6]    
    p["x"]=[data_from_f[i][7],data_from_f[i][8]]
    p["<z>"]=data_from_f[i][9]
    p["z"]=[data_from_f[i][10],data_from_f[i][11]]
    p["xSec"]=data_from_f[i][0]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters
    if p["xSec"]<0.00000001:
        p["thFactor"]=1.
    else:
        p["thFactor"]=1/(p["pT"][1]-p["pT"][0])/(p["z"][1]-p["z"][0])/data_from_f[i][15]#devide by bin size multiply by DIS xSec
         
    p["uncorrErr"].append(data_from_f[i][1])
    p["uncorrErr"].append(data_from_f[i][2])
    if(addDISnormalizationUncertainty):
        if(data_from_f[i][0]==0):
            pass#p.corrErrors.append(0)
        else:
            p["corrErr"].append((data_from_f[i][16]+data_from_f[i][17])/2/data_from_f[i][15]*data_from_f[i][0])
    #
    if p["xSec"]<0.00000001:
        print (p["z"],p["<z>"],p["xSec"], "(point dropped)")
    else:
        DataCurrent.AddPoint(p)   

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%
###############################################################################
###########################hermes.proton.zxpt-3D K+###########################
print("hermes.proton.zxpt-3D K+ file ...")
f = open(path_to_data+"hermes.proton.zxpt-3D.no-vmsub.mults_kplus.datnew_fortran")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Convert to numbers ...")

del data_from_f[0:5]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split()
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")


print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.p.no-vmsub.zxpt.k+',"SIDIS")
DataCurrent.comment="HERMES proton-to-k+ (zxpt-3D) no-vmsub. thFactor contains DIS normalization (by MSHT20nnlo)"
DataCurrent.reference="1212.5407"

proc_current=[1,1,2,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.85,10.,10000.]

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_from_f[i][12]
    p["pT"]=[data_from_f[i][13],data_from_f[i][14]]
    p["<Q>"]=numpy.sqrt(data_from_f[i][3])
    p["Q"]=[numpy.sqrt(data_from_f[i][4]),numpy.sqrt(data_from_f[i][5])]
    p["<x>"]=data_from_f[i][6]    
    p["x"]=[data_from_f[i][7],data_from_f[i][8]]
    p["<z>"]=data_from_f[i][9]    
    p["z"]=[data_from_f[i][10],data_from_f[i][11]]
    p["xSec"]=data_from_f[i][0]
    p["M_target"]=M_proton
    p["M_product"]=m_kaon
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters
    if p["xSec"]<0.00000001:
        p["thFactor"]=1.
    else:
        p["thFactor"]=1/(p["pT"][1]-p["pT"][0])/(p["z"][1]-p["z"][0])/data_from_f[i][15]#devide by bin size multiply by DIS xSec
         
    p["uncorrErr"].append(data_from_f[i][1])
    p["uncorrErr"].append(data_from_f[i][2])
    if(addDISnormalizationUncertainty):
        if(data_from_f[i][0]==0):
            pass#p.corrErrors.append(0)
        else:
            p["corrErr"].append((data_from_f[i][16]+data_from_f[i][17])/2/data_from_f[i][15]*data_from_f[i][0])
    #
    if p["xSec"]<0.00000001:
        print(p["z"],p["<z>"],p["xSec"], "(point dropped)")
    else:
        DataCurrent.AddPoint(p)   

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%
###############################################################################
###########################hermes.proton.zxpt-3D K-###########################
print("hermes.proton.zxpt-3D K- file ...")
f = open(path_to_data+"hermes.proton.zxpt-3D.no-vmsub.mults_kminus.datnew_fortran")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Convert to numbers ...")

del data_from_f[0:5]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split()
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")


print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.p.no-vmsub.zxpt.k-',"SIDIS")
DataCurrent.comment="HERMES proton-to-k- (zxpt-3D) no-vmsub. thFactor contains DIS normalization (by MSHT20nnlo)"
DataCurrent.reference="1212.5407"

proc_current=[1,1,-2,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.85,10.,10000.]

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_from_f[i][12]
    p["pT"]=[data_from_f[i][13],data_from_f[i][14]]
    p["<Q>"]=numpy.sqrt(data_from_f[i][3])
    p["Q"]=[numpy.sqrt(data_from_f[i][4]),numpy.sqrt(data_from_f[i][5])]
    p["<x>"]=data_from_f[i][6]    
    p["x"]=[data_from_f[i][7],data_from_f[i][8]]
    p["<z>"]=data_from_f[i][9]
    p["z"]=[data_from_f[i][10],data_from_f[i][11]]
    p["xSec"]=data_from_f[i][0]
    p["M_target"]=M_proton
    p["M_product"]=m_kaon
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters
    if p["xSec"]<0.00000001:
        p["thFactor"]=1.
    else:
        p["thFactor"]=1/(p["pT"][1]-p["pT"][0])/(p["z"][1]-p["z"][0])/data_from_f[i][15]#devide by bin size multiply by DIS xSec
         
    p["uncorrErr"].append(data_from_f[i][1])
    p["uncorrErr"].append(data_from_f[i][2])
    if(addDISnormalizationUncertainty):
        if(data_from_f[i][0]==0):
            pass#p.corrErrors.append(0)
        else:
            p["corrErr"].append((data_from_f[i][16]+data_from_f[i][17])/2/data_from_f[i][15]*data_from_f[i][0])
    #
    if p["xSec"]<0.00000001:
        print(p["z"],p["<z>"],p["xSec"], "(point dropped)")
    else:
        DataCurrent.AddPoint(p)     

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%
###############################################################################
###########################hermes.proton.zxpt-3D pi+###########################
print("hermes.proton.vmsub.zxpt-3D pi+ file ...")
f = open(path_to_data+"hermes.proton.zxpt-3D.vmsub.mults_piplus.datnew_fortran")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Convert to numbers ...")

del data_from_f[0:5]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split()
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.p.vmsub.zxpt.pi+',"SIDIS")
DataCurrent.comment="HERMES proton-to-pi+ (zxpt-3D) vmsub. thFactor contains DIS normalization (by MSHT20nnlo)"
DataCurrent.reference="1212.5407"

proc_current=[1,1,1,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.85,10.,10000.]

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_from_f[i][12]
    p["pT"]=[data_from_f[i][13],data_from_f[i][14]]
    p["<Q>"]=numpy.sqrt(data_from_f[i][3])
    p["Q"]=[numpy.sqrt(data_from_f[i][4]),numpy.sqrt(data_from_f[i][5])]
    p["<x>"]=data_from_f[i][6]    
    p["x"]=[data_from_f[i][7],data_from_f[i][8]]
    p["<z>"]=data_from_f[i][9]
    p["z"]=[data_from_f[i][10],data_from_f[i][11]]
    p["xSec"]=data_from_f[i][0]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters
    if p["xSec"]<0.00000001:
        p["thFactor"]=1.
    else:
        p["thFactor"]=1/(p["pT"][1]-p["pT"][0])/(p["z"][1]-p["z"][0])/data_from_f[i][15]#devide by bin size multiply by DIS xSec
         
    p["uncorrErr"].append(data_from_f[i][1])
    p["uncorrErr"].append(data_from_f[i][2])
    if(addDISnormalizationUncertainty):
        if(data_from_f[i][0]==0):
            pass#p.corrErrors.append(0)
        else:
            p["corrErr"].append((data_from_f[i][16]+data_from_f[i][17])/2/data_from_f[i][15]*data_from_f[i][0])
    #
    if p["xSec"]<0.00000001:
        print(p["z"],p["<z>"],p["xSec"], "(point dropped)")
    else:
        DataCurrent.AddPoint(p)   

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%
###############################################################################
###########################hermes.proton.zxpt-3D pi-###########################
print("hermes.proton.zxpt-3D pi- file ...")
f = open(path_to_data+"hermes.proton.zxpt-3D.vmsub.mults_piminus.datnew_fortran")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Convert to numbers ...")

del data_from_f[0:5]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split()
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")


print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.p.vmsub.zxpt.pi-',"SIDIS")
DataCurrent.comment="HERMES proton-to-pi- (zxpt-3D) vmsub. thFactor contains DIS normalization (by MSHT20nnlo)"
DataCurrent.reference="1212.5407"

proc_current=[1,1,-1,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.85,10.,10000.]

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_from_f[i][12]
    p["pT"]=[data_from_f[i][13],data_from_f[i][14]]
    p["<Q>"]=numpy.sqrt(data_from_f[i][3])
    p["Q"]=[numpy.sqrt(data_from_f[i][4]),numpy.sqrt(data_from_f[i][5])]
    p["<x>"]=data_from_f[i][6]    
    p["x"]=[data_from_f[i][7],data_from_f[i][8]]
    p["<z>"]=data_from_f[i][9]
    p["z"]=[data_from_f[i][10],data_from_f[i][11]]
    p["xSec"]=data_from_f[i][0]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters
    if p["xSec"]<0.00000001:
        p["thFactor"]=1.
    else:
        p["thFactor"]=1/(p["pT"][1]-p["pT"][0])/(p["z"][1]-p["z"][0])/data_from_f[i][15]#devide by bin size multiply by DIS xSec
         
    p["uncorrErr"].append(data_from_f[i][1])
    p["uncorrErr"].append(data_from_f[i][2])
    if(addDISnormalizationUncertainty):
        if(data_from_f[i][0]==0):
            pass#p.corrErrors.append(0)
        else:
            p["corrErr"].append((data_from_f[i][16]+data_from_f[i][17])/2/data_from_f[i][15]*data_from_f[i][0])
    #
    if p["xSec"]<0.00000001:
        print(p["z"],p["<z>"],p["xSec"], "(point dropped)")
    else:
        DataCurrent.AddPoint(p)    

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%
###############################################################################
###########################hermes.proton.zxpt-3D K+###########################
print("hermes.proton.zxpt-3D K+ file ...")
f = open(path_to_data+"hermes.proton.zxpt-3D.vmsub.mults_kplus.datnew_fortran")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Convert to numbers ...")

del data_from_f[0:5]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split()
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")


print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.p.vmsub.zxpt.k+',"SIDIS")
DataCurrent.comment="HERMES proton-to-k+ (zxpt-3D) vmsub. thFactor contains DIS normalization (by MSHT20nnlo)"
DataCurrent.reference="1212.5407"

proc_current=[1,1,2,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.85,10.,10000.]

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_from_f[i][12]
    p["pT"]=[data_from_f[i][13],data_from_f[i][14]]
    p["<Q>"]=numpy.sqrt(data_from_f[i][3])
    p["Q"]=[numpy.sqrt(data_from_f[i][4]),numpy.sqrt(data_from_f[i][5])]
    p["<x>"]=data_from_f[i][6]    
    p["x"]=[data_from_f[i][7],data_from_f[i][8]]
    p["<z>"]=data_from_f[i][9]
    p["z"]=[data_from_f[i][10],data_from_f[i][11]]
    p["xSec"]=data_from_f[i][0]
    p["M_target"]=M_proton
    p["M_product"]=m_kaon
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters
    if p["xSec"]<0.00000001:
        p["thFactor"]=1.
    else:
        p["thFactor"]=1/(p["pT"][1]-p["pT"][0])/(p["z"][1]-p["z"][0])/data_from_f[i][15]#devide by bin size multiply by DIS xSec
         
    p["uncorrErr"].append(data_from_f[i][1])
    p["uncorrErr"].append(data_from_f[i][2])
    if(addDISnormalizationUncertainty):
        if(data_from_f[i][0]==0):
            pass#p.corrErrors.append(0)
        else:
            p["corrErr"].append((data_from_f[i][16]+data_from_f[i][17])/2/data_from_f[i][15]*data_from_f[i][0])
    #
    if p["xSec"]<0.00000001:
        print(p["z"],p["<z>"],p["xSec"], "(point dropped)")
    else:
        DataCurrent.AddPoint(p)     

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%
###############################################################################
###########################hermes.proton.zxpt-3D K-###########################
print("hermes.proton.zxpt-3D K- file ...")
f = open(path_to_data+"hermes.proton.zxpt-3D.vmsub.mults_kminus.datnew_fortran")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Convert to numbers ...")

del data_from_f[0:5]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split()
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")


print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.p.vmsub.zxpt.k-',"SIDIS")
DataCurrent.comment="HERMES proton-to-k- (zxpt-3D) vmsub. thFactor contains DIS normalization (by MSHT20nnlo)"
DataCurrent.reference="1212.5407"

proc_current=[1,1,-2,2001]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.85,10.,10000.]

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_from_f[i][12]
    p["pT"]=[data_from_f[i][13],data_from_f[i][14]]
    p["<Q>"]=numpy.sqrt(data_from_f[i][3])
    p["Q"]=[numpy.sqrt(data_from_f[i][4]),numpy.sqrt(data_from_f[i][5])]
    p["<x>"]=data_from_f[i][6]    
    p["x"]=[data_from_f[i][7],data_from_f[i][8]]
    p["<z>"]=data_from_f[i][9]
    p["z"]=[data_from_f[i][10],data_from_f[i][11]]
    p["xSec"]=data_from_f[i][0]
    p["M_target"]=M_proton
    p["M_product"]=m_kaon
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters
    if p["xSec"]<0.00000001:
        p["thFactor"]=1.
    else:
        p["thFactor"]=1/(p["pT"][1]-p["pT"][0])/(p["z"][1]-p["z"][0])/data_from_f[i][15]#devide by bin size multiply by DIS xSec
         
    p["uncorrErr"].append(data_from_f[i][1])
    p["uncorrErr"].append(data_from_f[i][2])
    if(addDISnormalizationUncertainty):
        if(data_from_f[i][0]==0):
            pass#p.corrErrors.append(0)
        else:
            p["corrErr"].append((data_from_f[i][16]+data_from_f[i][17])/2/data_from_f[i][15]*data_from_f[i][0])
    
    if p["xSec"]<0.00000001:
        print(p["z"],p["<z>"],p["xSec"], "(point dropped)")
    else:
        DataCurrent.AddPoint(p)      

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

###################################################################################################
#############################HERMES deuteron########################################################
###################################################################################################
#%%
###############################################################################
###########################hermes.deuteron.zxpt-3D pi+###########################
print("hermes.deuteron.zxpt-3D pi+ file ...")
f = open(path_to_data+"hermes.deuteron.zxpt-3D.no-vmsub.mults_piplus.datnew_fortran")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Convert to numbers ...")

del data_from_f[0:5]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split()
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.d.no-vmsub.zxpt.pi+',"SIDIS")
DataCurrent.comment="HERMES deutron-to-pi+ (zxpt-3D) no-vmsub. thFactor contains DIS normalization (by MSHT20nnlo)"
DataCurrent.reference="1212.5407"

proc_current=[1,1,1,2002]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.85,10.,10000.]

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_from_f[i][12]
    p["pT"]=[data_from_f[i][13],data_from_f[i][14]]
    p["<Q>"]=numpy.sqrt(data_from_f[i][3])
    p["Q"]=[numpy.sqrt(data_from_f[i][4]),numpy.sqrt(data_from_f[i][5])]
    p["<x>"]=data_from_f[i][6]    
    p["x"]=[data_from_f[i][7],data_from_f[i][8]]
    p["<z>"]=data_from_f[i][9]
    p["z"]=[data_from_f[i][10],data_from_f[i][11]]
    p["xSec"]=data_from_f[i][0]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters
    if p["xSec"]<0.00000001:
        p["thFactor"]=1.
    else:
        p["thFactor"]=1/(p["pT"][1]-p["pT"][0])/(p["z"][1]-p["z"][0])/data_from_f[i][15]#devide by bin size multiply by DIS xSec
         
    p["uncorrErr"].append(data_from_f[i][1])
    p["uncorrErr"].append(data_from_f[i][2])
    if(addDISnormalizationUncertainty):
        if(data_from_f[i][0]==0):
            pass#p.corrErrors.append(0)
        else:
            p["corrErr"].append((data_from_f[i][16]+data_from_f[i][17])/2/data_from_f[i][15]*data_from_f[i][0])
    
    if p["xSec"]<0.00000001:
        print(p["z"],p["<z>"],p["xSec"], "(point dropped)")
    else:
        DataCurrent.AddPoint(p)   

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%
###############################################################################
###########################hermes.deuteron.zxpt-3D pi-###########################
print("hermes.deuteron.zxpt-3D pi- file ...")
f = open(path_to_data+"hermes.deuteron.zxpt-3D.no-vmsub.mults_piminus.datnew_fortran")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Convert to numbers ...")

del data_from_f[0:5]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split()
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")


print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.d.no-vmsub.zxpt.pi-',"SIDIS")
DataCurrent.comment="HERMES deutron-to-pi- (zxpt-3D) no-vmsub. thFactor contains DIS normalization (by MSHT20nnlo)"
DataCurrent.reference="1212.5407"

proc_current=[1,1,-1,2002]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.85,10.,10000.]

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_from_f[i][12]
    p["pT"]=[data_from_f[i][13],data_from_f[i][14]]
    p["<Q>"]=numpy.sqrt(data_from_f[i][3])
    p["Q"]=[numpy.sqrt(data_from_f[i][4]),numpy.sqrt(data_from_f[i][5])]
    p["<x>"]=data_from_f[i][6]    
    p["x"]=[data_from_f[i][7],data_from_f[i][8]]
    p["<z>"]=data_from_f[i][9]
    p["z"]=[data_from_f[i][10],data_from_f[i][11]]
    p["xSec"]=data_from_f[i][0]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters
    if p["xSec"]<0.00000001:
        p["thFactor"]=1.
    else:
        p["thFactor"]=1/(p["pT"][1]-p["pT"][0])/(p["z"][1]-p["z"][0])/data_from_f[i][15]#devide by bin size multiply by DIS xSec
         
    p["uncorrErr"].append(data_from_f[i][1])
    p["uncorrErr"].append(data_from_f[i][2])
    if(addDISnormalizationUncertainty):
        if(data_from_f[i][0]==0):
            pass#p.corrErrors.append(0)
        else:
            p["corrErr"].append((data_from_f[i][16]+data_from_f[i][17])/2/data_from_f[i][15]*data_from_f[i][0])
    
    if p["xSec"]<0.00000001:
        print(p["z"],p["<z>"],p["xSec"], "(point dropped)")
    else:
        DataCurrent.AddPoint(p)   

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%
###############################################################################
###########################hermes.deuteron.zxpt-3D K+###########################
print("hermes.deuteron.zxpt-3D K+ file ...")
f = open(path_to_data+"hermes.deuteron.zxpt-3D.no-vmsub.mults_kplus.datnew_fortran")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Convert to numbers ...")

del data_from_f[0:5]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split()
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")


print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.d.no-vmsub.zxpt.k+',"SIDIS")
DataCurrent.comment="HERMES deutron-to-k+ (zxpt-3D) no-vmsub. thFactor contains DIS normalization (by MSHT20nnlo)"
DataCurrent.reference="1212.5407"

proc_current=[1,1,2,2002]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.85,10.,10000.]

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_from_f[i][12]
    p["pT"]=[data_from_f[i][13],data_from_f[i][14]]
    p["<Q>"]=numpy.sqrt(data_from_f[i][3])
    p["Q"]=[numpy.sqrt(data_from_f[i][4]),numpy.sqrt(data_from_f[i][5])]
    p["<x>"]=data_from_f[i][6]    
    p["x"]=[data_from_f[i][7],data_from_f[i][8]]
    p["<z>"]=data_from_f[i][9]
    p["z"]=[data_from_f[i][10],data_from_f[i][11]]
    p["xSec"]=data_from_f[i][0]
    p["M_target"]=M_proton
    p["M_product"]=m_kaon
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters
    if p["xSec"]<0.00000001:
        p["thFactor"]=1.
    else:
        p["thFactor"]=1/(p["pT"][1]-p["pT"][0])/(p["z"][1]-p["z"][0])/data_from_f[i][15]#devide by bin size multiply by DIS xSec
         
    p["uncorrErr"].append(data_from_f[i][1])
    p["uncorrErr"].append(data_from_f[i][2])
    if(addDISnormalizationUncertainty):
        if(data_from_f[i][0]==0):
            pass#p.corrErrors.append(0)
        else:
            p["corrErr"].append((data_from_f[i][16]+data_from_f[i][17])/2/data_from_f[i][15]*data_from_f[i][0])
    #
    if p["xSec"]<0.00000001:
        print(p["z"],p["<z>"],p["xSec"], "(point dropped)")
    else:
        DataCurrent.AddPoint(p)   

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%
###############################################################################
###########################hermes.deuteron.zxpt-3D K-###########################
print("hermes.deuteron.zxpt-3D K- file ..." )
f = open(path_to_data+"hermes.deuteron.zxpt-3D.no-vmsub.mults_kminus.datnew_fortran")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Convert to numbers ...")

del data_from_f[0:5]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split()
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")


print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.d.no-vmsub.zxpt.k-',"SIDIS")
DataCurrent.comment="HERMES deutron-to-k- (zxpt-3D) no-vmsub. thFactor contains DIS normalization (by MSHT20nnlo)"
DataCurrent.reference="1212.5407"


proc_current=[1,1,-2,2002]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.85,10.,10000.]

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_from_f[i][12]
    p["pT"]=[data_from_f[i][13],data_from_f[i][14]]
    p["<Q>"]=numpy.sqrt(data_from_f[i][3])
    p["Q"]=[numpy.sqrt(data_from_f[i][4]),numpy.sqrt(data_from_f[i][5])]
    p["<x>"]=data_from_f[i][6]    
    p["x"]=[data_from_f[i][7],data_from_f[i][8]]
    p["<z>"]=data_from_f[i][9]
    p["z"]=[data_from_f[i][10],data_from_f[i][11]]
    p["xSec"]=data_from_f[i][0]
    p["M_target"]=M_proton
    p["M_product"]=m_kaon
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters
    if p["xSec"]<0.00000001:
        p["thFactor"]=1.
    else:
        p["thFactor"]=1/(p["pT"][1]-p["pT"][0])/(p["z"][1]-p["z"][0])/data_from_f[i][15]#devide by bin size multiply by DIS xSec
         
    p["uncorrErr"].append(data_from_f[i][1])
    p["uncorrErr"].append(data_from_f[i][2])
    if(addDISnormalizationUncertainty):
        if(data_from_f[i][0]==0):
            pass#p.corrErrors.append(0)
        else:
            p["corrErr"].append((data_from_f[i][16]+data_from_f[i][17])/2/data_from_f[i][15]*data_from_f[i][0])
    #
    if p["xSec"]<0.00000001:
        print(p["z"],p["<z>"],p["xSec"], "(point dropped)")
    else:
        DataCurrent.AddPoint(p)     

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%
###############################################################################
###########################hermes.deuteron.zxpt-3D pi+###########################
print("hermes.deuteron.vmsub.zxpt-3D pi+ file ...")
f = open(path_to_data+"hermes.deuteron.zxpt-3D.vmsub.mults_piplus.datnew_fortran")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Convert to numbers ...")

del data_from_f[0:5]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split()
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.d.vmsub.zxpt.pi+',"SIDIS")
DataCurrent.comment="HERMES deutron-to-pi+ (zxpt-3D) vmsub. thFactor contains DIS normalization (by MSHT20nnlo)"
DataCurrent.reference="1212.5407"

proc_current=[1,1,1,2002]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.85,10.,10000.]

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_from_f[i][12]
    p["pT"]=[data_from_f[i][13],data_from_f[i][14]]
    p["<Q>"]=numpy.sqrt(data_from_f[i][3])
    p["Q"]=[numpy.sqrt(data_from_f[i][4]),numpy.sqrt(data_from_f[i][5])]
    p["<x>"]=data_from_f[i][6]    
    p["x"]=[data_from_f[i][7],data_from_f[i][8]]
    p["<z>"]=data_from_f[i][9]
    p["z"]=[data_from_f[i][10],data_from_f[i][11]]
    p["xSec"]=data_from_f[i][0]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters
    if p["xSec"]<0.00000001:
        p["thFactor"]=1.
    else:
        p["thFactor"]=1/(p["pT"][1]-p["pT"][0])/(p["z"][1]-p["z"][0])/data_from_f[i][15]#devide by bin size multiply by DIS xSec
         
    p["uncorrErr"].append(data_from_f[i][1])
    p["uncorrErr"].append(data_from_f[i][2])
    if(addDISnormalizationUncertainty):
        if(data_from_f[i][0]==0):
            pass#p.corrErrors.append(0)
        else:
            p["corrErr"].append((data_from_f[i][16]+data_from_f[i][17])/2/data_from_f[i][15]*data_from_f[i][0])
    #
    if p["xSec"]<0.00000001:
        print(p["z"],p["<z>"],p["xSec"], "(point dropped)")
    else:
        DataCurrent.AddPoint(p)   

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%
###############################################################################
###########################hermes.deuteron.zxpt-3D pi-###########################
print("hermes.deuteron.zxpt-3D pi- file ...")
f = open(path_to_data+"hermes.deuteron.zxpt-3D.vmsub.mults_piminus.datnew_fortran")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Convert to numbers ...")

del data_from_f[0:5]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split()
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")


print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.d.vmsub.zxpt.pi-',"SIDIS")
DataCurrent.comment="HERMES deutron-to-pi- (zxpt-3D) vmsub. thFactor contains DIS normalization (by MSHT20nnlo)"
DataCurrent.reference="1212.5407"

proc_current=[1,1,-1,2002]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.85,10.,10000.]

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_from_f[i][12]
    p["pT"]=[data_from_f[i][13],data_from_f[i][14]]
    p["<Q>"]=numpy.sqrt(data_from_f[i][3])
    p["Q"]=[numpy.sqrt(data_from_f[i][4]),numpy.sqrt(data_from_f[i][5])]
    p["<x>"]=data_from_f[i][6]    
    p["x"]=[data_from_f[i][7],data_from_f[i][8]]
    p["<z>"]=data_from_f[i][9]
    p["z"]=[data_from_f[i][10],data_from_f[i][11]]
    p["xSec"]=data_from_f[i][0]
    p["M_target"]=M_proton
    p["M_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters
    if p["xSec"]<0.00000001:
        p["thFactor"]=1.
    else:
        p["thFactor"]=1/(p["pT"][1]-p["pT"][0])/(p["z"][1]-p["z"][0])/data_from_f[i][15]#devide by bin size multiply by DIS xSec
         
    p["uncorrErr"].append(data_from_f[i][1])
    p["uncorrErr"].append(data_from_f[i][2])
    if(addDISnormalizationUncertainty):
        if(data_from_f[i][0]==0):
            pass#p.corrErrors.append(0)
        else:
            p["corrErr"].append((data_from_f[i][16]+data_from_f[i][17])/2/data_from_f[i][15]*data_from_f[i][0])
    #
    if p["xSec"]<0.00000001:
        print(p["z"],p["<z>"],p["xSec"], "(point dropped)")
    else:
        DataCurrent.AddPoint(p)      

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%
###############################################################################
###########################hermes.proton.zxpt-3D K+###########################
print("hermes.deuteron.zxpt-3D K+ file ..." )
f = open(path_to_data+"hermes.deuteron.zxpt-3D.vmsub.mults_kplus.datnew_fortran")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Convert to numbers ...")

del data_from_f[0:5]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split()
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")


print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.d.vmsub.zxpt.k+',"SIDIS")
DataCurrent.comment="HERMES deutron-to-k+ (zxpt-3D) vmsub. thFactor contains DIS normalization (by MSHT20nnlo)"
DataCurrent.reference="1212.5407"

proc_current=[1,1,2,2002]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.85,10.,10000.]

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_from_f[i][12]
    p["pT"]=[data_from_f[i][13],data_from_f[i][14]]
    p["<Q>"]=numpy.sqrt(data_from_f[i][3])
    p["Q"]=[numpy.sqrt(data_from_f[i][4]),numpy.sqrt(data_from_f[i][5])]
    p["<x>"]=data_from_f[i][6]    
    p["x"]=[data_from_f[i][7],data_from_f[i][8]]
    p["<z>"]=data_from_f[i][9]
    p["z"]=[data_from_f[i][10],data_from_f[i][11]]
    p["xSec"]=data_from_f[i][0]
    p["M_target"]=M_proton
    p["M_product"]=m_kaon
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters
    if p["xSec"]<0.00000001:
        p["thFactor"]=1.
    else:
        p["thFactor"]=1/(p["pT"][1]-p["pT"][0])/(p["z"][1]-p["z"][0])/data_from_f[i][15]#devide by bin size multiply by DIS xSec
         
    p["uncorrErr"].append(data_from_f[i][1])
    p["uncorrErr"].append(data_from_f[i][2])
    if(addDISnormalizationUncertainty):
        if(data_from_f[i][0]==0):
            pass#p.corrErrors.append(0)
        else:
            p["corrErr"].append((data_from_f[i][16]+data_from_f[i][17])/2/data_from_f[i][15]*data_from_f[i][0])
    #
    if p["xSec"]<0.00000001:
        print(p["z"],p["<z>"],p["xSec"], "(point dropped)")
    else:
        DataCurrent.AddPoint(p)    
print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%
###############################################################################
###########################hermes.deuteron.zxpt-3D K-###########################
print("hermes.deuteron.zxpt-3D K- file ...")
f = open(path_to_data+"hermes.deuteron.zxpt-3D.vmsub.mults_kminus.datnew_fortran")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Convert to numbers ...")

del data_from_f[0:5]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split()
    del data_from_f[i][0]
    data_from_f[i]=[float(j) for j in data_from_f[i]]
    #k.spit("\t")


print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('hermes.d.vmsub.zxpt.k-',"SIDIS")
DataCurrent.comment="HERMES deutron-to-k- (zxpt-3D) vmsub. thFactor contains DIS normalization (by MSHT20nnlo)"
DataCurrent.reference="1212.5407"

proc_current=[1,1,-2,2002]
s_current=2*27.6*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.85,10.,10000.]

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current
    p["<pT>"]=data_from_f[i][12]
    p["pT"]=[data_from_f[i][13],data_from_f[i][14]]
    p["<Q>"]=numpy.sqrt(data_from_f[i][3])
    p["Q"]=[numpy.sqrt(data_from_f[i][4]),numpy.sqrt(data_from_f[i][5])]
    p["<x>"]=data_from_f[i][6]    
    p["x"]=[data_from_f[i][7],data_from_f[i][8]]
    p["<z>"]=data_from_f[i][9]
    p["z"]=[data_from_f[i][10],data_from_f[i][11]]
    p["xSec"]=data_from_f[i][0]
    p["M_target"]=M_proton
    p["M_product"]=m_kaon
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters
    if p["xSec"]<0.00000001:
        p["thFactor"]=1.
    else:
        p["thFactor"]=1/(p["pT"][1]-p["pT"][0])/(p["z"][1]-p["z"][0])/data_from_f[i][15]#devide by bin size multiply by DIS xSec
         
    p["uncorrErr"].append(data_from_f[i][1])
    p["uncorrErr"].append(data_from_f[i][2])
    if(addDISnormalizationUncertainty):
        if(data_from_f[i][0]==0):
            pass#p.corrErrors.append(0)
        else:
            p["corrErr"].append((data_from_f[i][16]+data_from_f[i][17])/2/data_from_f[i][15]*data_from_f[i][0])
    #
    if p["xSec"]<0.00000001:
        print(p["z"],p["<z>"],p["xSec"], "(point dropped)")
    else:
        DataCurrent.AddPoint(p)   

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%
###############################################################################
###########################COMPASS#############################################
###############################################################################

###############################################################################
###########################compass.deuteron.h+#################################
print("compass.deuteron.h+ file ...")
f = open(path_to_data+"/compass(1706.01473).d.hplus.datnew_fortran")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Convert to numbers ...")

del data_from_f[0:3]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split()    
    data_from_f[i]=[float(j) for j in data_from_f[i]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('compass.d.h+',"SIDIS")
DataCurrent.comment="COMPASS isoscalar h+. thFactor contains DIS normalization (by MSHT20nnlo)"
DataCurrent.reference="1709.07374"

proc_current=[1,1,1,2103]
s_current=2*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.]

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current    
    p["x"]=[data_from_f[i][0],data_from_f[i][1]]    
    p["Q"]=[numpy.sqrt(data_from_f[i][2]),numpy.sqrt(data_from_f[i][3])]    
    p["z"]=[data_from_f[i][4],data_from_f[i][5]]
    p["pT"]=[numpy.sqrt(data_from_f[i][6]),numpy.sqrt(data_from_f[i][7])]    
    p["xSec"]=data_from_f[i][8]
    p["M_target"]=M_proton
    p["m_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters
    if p["xSec"]<0.00000001:
        p["thFactor"]=1
    else:
        p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["z"][1]-p["z"][0])/data_from_f[i][11]#devide by bin size multiply by DIS xSec
         
    p["uncorrErr"].append(data_from_f[i][9])
    p["uncorrErr"].append(data_from_f[i][10]*data_from_f[i][8]) #### for some unknown reasons the data-file gives relative uncertanty
    if(addDISnormalizationUncertainty):
        p["corrErr"].append((data_from_f[i][12]+data_from_f[i][13])/2/data_from_f[i][11]*data_from_f[i][8])
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
#%%
###############################################################################
###########################compass.deuteron.h-#################################
print("compass.deuteron.h- file ...")
f = open(path_to_data+"/compass(1706.01473).d.hminus.datnew_fortran")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()


print("Done.  =>     Convert to numbers ...")

del data_from_f[0:3]

for i in range(len(data_from_f)):
    data_from_f[i]=data_from_f[i].split()    
    data_from_f[i]=[float(j) for j in data_from_f[i]]

print("Done.  =>     Create points & append to data set ...")
DataCurrent=DataProcessor.DataSet.DataSet('compass.d.h-',"SIDIS")
DataCurrent.comment="COMPASS isoscalar h-. thFactor contains DIS normalization (by MSHT20nnlo)"
DataCurrent.reference="1709.07374"

proc_current=[1,1,-1,2103]
s_current=2*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.]

for i in range(len(data_from_f)):
    # makeup a point
    p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(i))
    #print DataCurrent.name+'.'+str(i)
    p["process"]=proc_current
    p["s"]=s_current    
    p["x"]=[data_from_f[i][0],data_from_f[i][1]]    
    p["Q"]=[numpy.sqrt(data_from_f[i][2]),numpy.sqrt(data_from_f[i][3])]    
    p["z"]=[data_from_f[i][4],data_from_f[i][5]]
    p["pT"]=[numpy.sqrt(data_from_f[i][6]),numpy.sqrt(data_from_f[i][7])]    
    p["xSec"]=data_from_f[i][8]
    p["M_target"]=M_proton
    p["m_product"]=m_pion
    p["includeCuts"]=includeCuts
    p["cutParams"]=cutParameters
    if p["xSec"]<0.00000001:
        p["thFactor"]=1
    else:
        p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["z"][1]-p["z"][0])/data_from_f[i][11]#devide by bin size multiply by DIS xSec
         
    p["uncorrErr"].append(data_from_f[i][9])
    p["uncorrErr"].append(data_from_f[i][10]*data_from_f[i][8]) #### for some unknown reasons the data-file gives relative uncertanty
    if(addDISnormalizationUncertainty):
        p["corrErr"].append((data_from_f[i][12]+data_from_f[i][13])/2/data_from_f[i][11]*data_from_f[i][8])
    #
    DataCurrent.AddPoint(p)

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")