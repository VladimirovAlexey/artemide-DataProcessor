#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 20:15:46 2019

Program that parse various DY data files to ADP-frendly formal

@author: vla18041
"""
import os
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))


import sys
sys.path.append(ROOT_DIR)

import DataProcessor.Point
import DataProcessor.DataSet

path_to_data=os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))+'/data/ATLAS/Zangular_1606_00689/'
path_to_save=ROOT_DIR+"/DataLib/DY_angular/"

M_Z=91.1876### mass of Z-boson

#%%
### given in the text
#proc_current=[1,1,5]
s_current=8000.**2
Q_current=[80.,100.]
incCut=False
cutParam=[20.,20.,2.,4.5]  


#%%
def FillTable(fileName,reg):
    
    DataCurrent.reference="arXiv:1606.00689"
    DataCurrent.isNormalized=False
    
    f = open(path_to_data+fileName)
    
    data_from_f=[]
    
    for line in f:    
        data_from_f.append(line.rstrip('\n'))
    
    f.close()
    
    print("Done.  =>     Convert to numbers ...")
    
    
    
    #define y-bin
    Wstring=next(x for x in data_from_f if x[0:8]=="#: Y BIN")
    if(Wstring[-2:]=="<1"):
        y_current=[0.,1.]
    elif(Wstring[-2:]=="<2"):
        y_current=[1.,2.]
    elif(Wstring[-2:]==".5"):
        y_current=[2.,3.5]
    else:
        print("MISTAKE IN DETERMINATION OF Y-BIN")
        print(Wstring)
        
    
    # read the data
    n=-1
    for s in data_from_f[:-1]:
        if(s[0]=="#" or s[0]=="P" or s==''): continue    
        n+=1
        data_here=[float(j) for j in s.split(",")]
        
        # makeup a point
        p=DataProcessor.Point.CreateDYPoint(DataCurrent.name+'.'+str(n))
        p["process"]=proc_current
        p["s"]=s_current
        p["qT"]=data_here[1:3]
        p["Q"]=Q_current
        p["<Q>"]=M_Z ### to be sure that mass middle value is Z-boson mass
        p["y"]=y_current
        p["thFactor"]=1.
        p["includeCuts"]=incCut
        p["cutParams"]=cutParam
        p["xSec"]=data_here[3]
        ###### The uncertanty is statistic + systematic. 
        ###### Systematic is lightly correlated I ignore it
        p["uncorrErr"].append((data_here[4]-data_here[5])/2)
        p["uncorrErr"].append((data_here[6]-data_here[7])/2)
        if(reg):
            p["uncorrErr"].append((data_here[8]-data_here[9])/2)
        DataCurrent.AddPoint(p)    
    
    print("Done.  ")

#%%
###############################################################################
########################### A4 y-differential########################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_Auu_0y1',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent Auu (used for normalization)"

proc_current=[1,1,1,2]

print("Read ATLAS8 Auu 0y1 file...")

FillTable("Table36.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#################################################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_Auu_1y2',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A4 (used for normalization)"

print("Read ATLAS8 Auu 1y2 file...")

FillTable("Table56.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#################################################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_Auu_2y35',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A4 (used for normalization)"

print("Read ATLAS8 Auu 2y35 file...")

FillTable("Table75.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")    

#%%
###############################################################################
########################### A4 y-differential########################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A4_0y1',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A4"

proc_current=[1,1,1,24]

print("Read ATLAS8 A4 0y1 file...")

FillTable("Table36.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#################################################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A4_1y2',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A4"

print("Read ATLAS8 A4 1y2 file...")

FillTable("Table56.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#################################################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A4_2y35',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A4"

print("Read ATLAS8 A4 2y35 file...")

FillTable("Table75.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
########################### A0 y-differential########################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A0_0y1',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A0"

proc_current=[1,1,1,20]

print("Read ATLAS8 A0 0y1 file...")

FillTable("Table32.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#################################################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A0_1y2',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A0"


print("Read ATLAS8 A0 1y2 file...")

FillTable("Table52.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#################################################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A0_2y35',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A0"


print("Read ATLAS8 A0 2y35 file...")

FillTable("Table74.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
########################### A1 y-differential########################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A1_0y1',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A1"

proc_current=[1,1,1,21]

print("Read ATLAS8 A1 0y1 file...")

FillTable("Table33.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#################################################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A1_1y2',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A1"


print("Read ATLAS8 A1 1y2 file...")

FillTable("Table53.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#################################################################
##### there is no A1 for 2<y<3.5

#%%
###############################################################################
########################### A2 y-differential########################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A2_0y1',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A2"

proc_current=[1,1,1,22]

print("Read ATLAS8 A2 0y1 file...")

FillTable("Table34.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#################################################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A2_1y2',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A2"


print("Read ATLAS8 A2 1y2 file...")

FillTable("Table54.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#################################################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A2_2y35',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A2"


print("Read ATLAS8 A2 2y35 file...")

FillTable("Table66.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
########################### A3 y-differential########################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A3_0y1',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A3"

proc_current=[1,1,1,23]

print("Read ATLAS8 A3 0y1 file...")

FillTable("Table35.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#################################################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A3_1y2',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A3"


print("Read ATLAS8 A3 1y2 file...")

FillTable("Table55.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#################################################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A3_2y35',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A3"


print("Read ATLAS8 A3 2y35 file...")

FillTable("Table67.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
########################### A5 y-differential########################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A5_0y1',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A5"

proc_current=[1,1,1,25]

print("Read ATLAS8 A5 0y1 file...")

FillTable("Table37.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#################################################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A5_1y2',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A5"


print("Read ATLAS8 A5 1y2 file...")

FillTable("Table57.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#################################################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A5_2y35',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A5"


print("Read ATLAS8 A5 2y35 file...")

FillTable("Table68.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
########################### A6 y-differential########################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A6_0y1',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A6"

proc_current=[1,1,1,26]

print("Read ATLAS8 A6 0y1 file...")

FillTable("Table38.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#################################################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A6_1y2',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A6"


print("Read ATLAS8 A6 1y2 file...")

FillTable("Table58.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#################################################################
### there is no A6 measurements at 2<y<3.5

#%%
###############################################################################
########################### A7 y-differential########################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A7_0y1',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A7"

proc_current=[1,1,1,27]

print("Read ATLAS8 A7 0y1 file...")

FillTable("Table39.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#################################################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A7_1y2',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A7"


print("Read ATLAS8 A7 1y2 file...")

FillTable("Table59.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#################################################################
DataCurrent=DataProcessor.DataSet.DataSet('A8_A7_2y35',"DY")
DataCurrent.comment="ATLAS 8TeV angular coefficent A7"


print("Read ATLAS8 A7 2y35 file...")

FillTable("Table69.csv",True)

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
