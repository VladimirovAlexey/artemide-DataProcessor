#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 20 16:00:34 2025

This is a program which parses various DIS data to a standarized CVS format for easy handling and processing.

@author: guillermo
"""



import os
import sys
import numpy
import csv

ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))+"/"

##check the directory ROOT_DIR:
print("ROOT_DIR:",ROOT_DIR)

sys.path.append(ROOT_DIR+"artemide-DataProcessor/")


import DataProcessor.Point
import DataProcessor.DataSet

path_to_data="/data/arTeMiDe_Repository/data/g2Tables/"

#path_to_g2WW_data="/data/arTeMiDe_Repository/data/g2Tables/WW_Guillermo/"
path_to_g2WW_data="/data/arTeMiDe_Repository/data/g2Tables/WW_MAPPDFpol/"

path_to_save=ROOT_DIR+"/DataProcessor/DataLib/G2/"

print("The path where the file is saved is:",path_to_save)

## Relevant masses in GeV

M_proton=0.938
M_neutron=0.939

M_deuteron=(M_proton+M_neutron)/2
M_Helium3=M_neutron#2*M_proton+M_neutron

#%%
## FUNCTIONS TO GUESS BINS (RUDIMENTARY): 3 functions which take a string imput, a data table, a row and a column index and 
##compute a rough bin of the specified element assuming it's the average value of the variable in the bin. They cover 3 cases: 
    
    ## Type 1 --> Value is surrounded by data in an ordered way (usual bin)
    ## Type 2 --> Value is at the bottom of the table/ discontinuity of order begins with value
    ## Type 3 --> Value is at the top of the table/ discontinuity of order after value

## We distinguish between bins in Q² and x, because data shows <Q²> instead of <Q> we need to construct another set of functions 
## which do the same but taking the sqrt of the entries
    
"""
    Parameters:
    - var: only allows two options. "x" and "Q" to specify the variable to compute the bin for. Anything else will return an error
    - data: list of lists, each sublist represents a row..
    - row_index: index of the row containing our target element.
    - col_index: index of the column containing the data used to build the bin.

    Returns:
    - List [data_min, data_max]
"""

#--> TYPE 1: USUAL BIN.

def GuessBin_Middle(var, data, row_index, col_index):
 
 if var == "x":
        val = lambda k: float(data[k][col_index])
        
        return [val(row_index) - (val(row_index) - val(row_index - 1)) / 2,
                val(row_index) + (val(row_index + 1) - val(row_index)) / 2]
        
 elif var == "Q":
       val = lambda k: float(numpy.sqrt(data[k][col_index]))
       
       return [val(row_index) - (val(row_index) - val(row_index - 1)) / 2,
               val(row_index) + (val(row_index + 1) - val(row_index)) / 2]
 
 else:
        raise ValueError("Error: bin only admits 'x' or 'Q'")


#--> TYPE 2: VALUE AT BOTTOM/ FIRST AFTER ORDER DISCONTINUITY.

def GuessBin_Bottom(var, data, row_index, col_index):
 
 if var == "x":
        val = lambda k: float(data[k][col_index])
        
        return [val(row_index) - (val(row_index + 1) - val(row_index)) / 2,
                val(row_index) + (val(row_index + 1) - val(row_index)) / 2]
        
 elif var == "Q":
       val = lambda k: float(numpy.sqrt(data[k][col_index]))
       
       return [val(row_index) - (val(row_index + 1) - val(row_index)) / 2,
               val(row_index) + (val(row_index + 1) - val(row_index)) / 2]
 
 else:
        raise ValueError("Error: bin only admits'x' or 'Q'")


#--> TYPE 3: VALUE AT TOP/ LAST BEFORE ORDER DISCONTINUITY.

def GuessBin_Top(var, data, row_index, col_index):
 
 if var == "x":
        val = lambda k: float(data[k][col_index])
        
        return [val(row_index) - (val(row_index) - val(row_index - 1)) / 2,
                val(row_index) + (val(row_index) - val(row_index - 1)) / 2]
        
 elif var == "Q":
       val = lambda k: float(numpy.sqrt(data[k][col_index]))
       
       return [val(row_index) - (val(row_index) - val(row_index - 1)) / 2,
               val(row_index) + (val(row_index) - val(row_index - 1)) / 2]
 
 else:
        raise ValueError("Error: bin only admits'x' or 'Q'")

#%%
#############################################################
#########-------------------------------------------#########
######--------------E143 1998 proton data--------------######
#########-------------------------------------------#########
#############################################################

## NOTE: FOR THE E143 EXPERIMENTS WE NEED NOT COMPUTE THE WW TERM. THIS IS BECAUSE THEY PRESENT THE HIGHER TWIST PART FOR g2 ALONG BOTH ITS 
## UNCERTAINTIES IN THE PDF (NOT H.E.P DATA). AS A CONSEQUENCE WE JUST NEED TO IMPORT THE KINEMATIC COVERAGE, BUILD THE BIN ON Q AND 
## CONSTRUCT THE FULL TABLES FOR THE EXPERIMENT JOINING THE WHOLE DATA BUNCH.

## First we read the H.E.P data we stored, from which we can extract the kinematic coverage:
print("Reading E143 1998 proton data file ...")
f = open(path_to_data+"/E143/1998/Table31.csv")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
# Now we choose the lines in which the data is stored and pick only those:
data_current=data_from_f[14:26] 

for i in range(len(data_current)):    
    data_current[i]=data_current[i].split(",")    
    data_current[i]=[float(j) for j in data_current[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

#%%
# Now we need to guess the bin on Q :
    
def Qbounds(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current,i,3)
    
    elif 0<i<6:
        return GuessBin_Middle("Q",data_current,i,3)
    
    elif i==6:
        return GuessBin_Top("Q",data_current,i,3)
    
    elif i==7:
        return GuessBin_Bottom("Q",data_current,i,3)
    
    elif 7<i<len(data_current)-1:
        return GuessBin_Middle("Q",data_current,i,3)
    
    elif i== len(data_current)-1:
       return GuessBin_Top("Q",data_current,i,3)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

#%%
## Now we import the twist-3 data we have created from the experimental paper:
print("Reading the twist-3 information for the proton data file ...")
f = open(path_to_data+"E143/1998/E143_1998_g2bar_proton.csv")

twist3_data=[]

for line in f:    
    twist3_data.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

for i in range(len(twist3_data)):    
    twist3_data[i]=twist3_data[i].split(",")    
    twist3_data[i]=[float(j) for j in twist3_data[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

print("Done, proton's twist3 data converted to numbers")

#%%
## Finally we can build the final table by joining both the kinematic and twist3 data:
print("Making G2 by E143_1998 for the proton...")

DataCurrent=DataProcessor.DataSet.DataSet('E143.p',"G2")
DataCurrent.comment="g2bar data taken from page 86 in the paper. Bin in Q guessed by us. Additional normalization uncertainty 3.7%. The E143 experiment is the only one which presents g2-bar, therefore the uncertainties are solely the stat+syst of said g2-bar, unlike the rest of the experimental tables"
DataCurrent.reference="hep-ph/9802357"

DataCurrent.isNormalized=False   ##Donde poner esta renormalización de la incertidumbre??
lumUncertainty=0.037
s_current=(2.*29.1*M_proton+(M_proton)**2)

for i in range(len(data_current)):
    # make up a point
    p1 = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p1["process"]=100
    p1["s"]=s_current  
    p1["<Q>"]=numpy.sqrt(data_current[i][3])     
    p1["Q"]=Qbounds(i)
    p1["<x>"]=data_current[i][0]    
    p1["x"]=[data_current[i][1],data_current[i][2]]
    p1["xSec"]=twist3_data[i][0]
    ## THE UNCERTAINTY IS PRESENTED AS: STAT(G2-BAR),SYST(G2-BAR)
    p1["uncorrErr"]=[twist3_data[i][1],twist3_data[i][2]]
    p1["thFactor"]=1.        
    #
    DataCurrent.AddPoint(p1)

print("Done.")
        
DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")

#%%
###############################################################
#########---------------------------------------------#########
######--------------E143 1998 deuteron data--------------######
#########---------------------------------------------#########
###############################################################

## NOTE: FOR THE E143 EXPERIMENTS WE NEED NOT COMPUTE THE WW TERM. THIS IS BECAUSE THEY PRESENT THE HIGHER TWIST PART FOR g2 ALONG BOTH ITS 
## UNCERTAINTIES IN THE PDF (NOT H.E.P DATA). AS A CONSEQUENCE WE JUST NEED TO IMPORT THE KINEMATIC COVERAGE, BUILD THE BIN ON Q AND 
## CONSTRUCT THE FULL TABLES FOR THE EXPERIMENT JOINING THE WHOLE DATA BUNCH.

## First we read the H.E.P data we stored, from which we can extract the kinematic coverage:
print("Reading E143 1998 deuteron data file ...")
f = open(path_to_data+"/E143/1998/Table32.csv")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
# Now we choose the lines in which the data is stored and pick only those:
data_current=data_from_f[12:24] 

for i in range(len(data_current)):    
    data_current[i]=data_current[i].split(",")    
    data_current[i]=[float(j) for j in data_current[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

#%%
# Now we need to guess the bin on Q :
    
def Qbounds(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current,i,3)
    
    elif 0<i<6:
        return GuessBin_Middle("Q",data_current,i,3)
    
    elif i==6:
        return GuessBin_Top("Q",data_current,i,3)
    
    elif i==7:
        return GuessBin_Bottom("Q",data_current,i,3)
    
    elif 7<i<len(data_current)-1:
        return GuessBin_Middle("Q",data_current,i,3)
    
    elif i== len(data_current)-1:
       return GuessBin_Top("Q",data_current,i,3)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

#%%
## Now we import the twist-3 data we have created from the experimental paper:
print("Reading the twist-3 information for the deuteron data file ...")
f = open(path_to_data+"E143/1998/E143_1998_g2bar_deuteron.csv")

twist3_data=[]

for line in f:    
    twist3_data.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

for i in range(len(twist3_data)):    
    twist3_data[i]=twist3_data[i].split(",")    
    twist3_data[i]=[float(j) for j in twist3_data[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

print("Done, deuteron's twist3 data converted to numbers")

#%%
## Finally we can build the final table by joining both the kinematic and twist3 data:
print("Making G2 by E143_1998 for the deuteron...")

DataCurrent=DataProcessor.DataSet.DataSet('E143.d',"G2")
DataCurrent.comment="g2bar data taken from page 86 in the paper. Bin in Q guessed by us. Additional normalization uncertainty 4.9%. The E143 experiment is the only one which presents g2-bar, therefore the uncertainties are solely the stat+syst of said g2-bar, unlike the rest of the experimental tables"
DataCurrent.reference="hep-ph/9802357"

DataCurrent.isNormalized=False   ##Donde poner esta renormalización de la incertidumbre??
lumUncertainty=0.049
s_current=2.*29.1*M_proton+(M_proton)**2

for i in range(len(data_current)):
    # make up a point
    p1 = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p1["process"]=102
    p1["s"]=s_current  
    p1["<Q>"]=numpy.sqrt(data_current[i][3])     
    p1["Q"]=Qbounds(i)
    p1["<x>"]=data_current[i][0]    
    p1["x"]=[data_current[i][1],data_current[i][2]]
    p1["xSec"]=twist3_data[i][0]
    ## THE UNCERTAINTY IS PRESENTED AS: STAT(G2-BAR),SYST(G2-BAR)
    p1["uncorrErr"]=[twist3_data[i][1],twist3_data[i][2]]
    p1["thFactor"]=1.        
    #
    DataCurrent.AddPoint(p1)
#
print("Done.")
DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")
    
#%%
##############################################################
#########--------------------------------------------#########
######--------------E143 1998 neutron data--------------######
#########--------------------------------------------#########
##############################################################

## NOTE: FOR THE E143 EXPERIMENTS WE NEED NOT COMPUTE THE WW TERM. THIS IS BECAUSE THEY PRESENT THE HIGHER TWIST PART FOR g2 ALONG BOTH ITS 
## UNCERTAINTIES IN THE PDF (NOT H.E.P DATA). AS A CONSEQUENCE WE JUST NEED TO IMPORT THE KINEMATIC COVERAGE, BUILD THE BIN ON Q AND 
## CONSTRUCT THE FULL TABLES FOR THE EXPERIMENT JOINING THE WHOLE DATA BUNCH.

## First we read the H.E.P data we stored, from which we can extract the kinematic coverage:
print("Reading E143 1998 neutron data file ...")
f = open(path_to_data+"/E143/1998/Table33.csv")

data_from_f=[]

for line in f:    
    data_from_f.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

#%%
# Now we choose the lines in which the data is stored and pick only those:
data_current=data_from_f[14:26] 

for i in range(len(data_current)):    
    data_current[i]=data_current[i].split(",")    
    data_current[i]=[float(j) for j in data_current[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

#%%
# Now we need to guess the bin on Q :
    
def Qbounds(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current,i,3)
    
    elif 0<i<6:
        return GuessBin_Middle("Q",data_current,i,3)
    
    elif i==6:
        return GuessBin_Top("Q",data_current,i,3)
    
    elif i==7:
        return GuessBin_Bottom("Q",data_current,i,3)
    
    elif 7<i<len(data_current)-1:
        return GuessBin_Middle("Q",data_current,i,3)
    
    elif i== len(data_current)-1:
       return GuessBin_Top("Q",data_current,i,3)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

#%%
## Now we import the twist-3 data we have created from the experimental paper:
print("Reading the twist-3 information for the neutron data file ...")
f = open(path_to_data+"E143/1998/E143_1998_g2bar_neutron.csv")

twist3_data=[]

for line in f:    
    twist3_data.append(line.rstrip('\n'))

f.close()

print("Done.  =>     Convert to numbers ...")

for i in range(len(twist3_data)):    
    twist3_data[i]=twist3_data[i].split(",")    
    twist3_data[i]=[float(j) for j in twist3_data[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

print("Done, neutron's twist3 data converted to numbers")

#%%
## Finally we can build the final table by joining both the kinematic and twist3 data:
print("Making G2 by E143_1998 for the neutron...")

DataCurrent=DataProcessor.DataSet.DataSet('E143.n',"G2")
DataCurrent.comment="g2bar data taken from page 87 in the paper. Bin in Q guessed by us. Additional normalization uncertainty 4.9%. The E143 experiment is the only one which presents g2-bar, therefore the uncertainties are solely the stat+syst of said g2-bar, unlike the rest of the experimental tables"
DataCurrent.reference="hep-ph/9802357"

DataCurrent.isNormalized=False   ##Donde poner esta renormalización de la incertidumbre??
lumUncertainty=0.049
s_current=(2.*29.1*M_neutron+(M_neutron)**2)

for i in range(len(data_current)):
    # make up a point
    p1 = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p1["process"]=101
    p1["s"]=s_current  
    p1["<Q>"]=numpy.sqrt(data_current[i][3])     
    p1["Q"]=Qbounds(i)
    p1["<x>"]=data_current[i][0]    
    p1["x"]=[data_current[i][1],data_current[i][2]]
    p1["xSec"]=twist3_data[i][0]
    ## THE UNCERTAINTY IS PRESENTED AS: STAT(G2-BAR),SYST(G2-BAR)
    p1["uncorrErr"]=[twist3_data[i][1],twist3_data[i][2]]
    p1["thFactor"]=1.        
    #
    DataCurrent.AddPoint(p1)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")

#%%

##############################################################
#########--------------------------------------------#########
######--------------E143 1996 proton data --------------######
#########--------------------------------------------#########
##############################################################

print("Reading E143 1996 proton data files ...")

f1 = open(path_to_data+"/E143/1996/Table1.csv")

data_from_f1=[]

for line in f1:    
    data_from_f1.append(line.rstrip('\n'))

f1.close()

print("Done.  =>     Convert table 1 to numbers ...")


data_current1=data_from_f1[25:32]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current1)):    
    data_current1[i]=data_current1[i].split(",")    
    data_current1[i]=[float(j) for j in data_current1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 1 converted to float")


f2 = open(path_to_data+"/E143/1996/Table2.csv")

data_from_f2=[]

for line in f2:    
    data_from_f2.append(line.rstrip('\n'))

f2.close()

print("Done.  =>     Convert table 2 to numbers ...")


data_current2=data_from_f2[23:28]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current2)):    
    data_current2[i]=data_current2[i].split(",")    
    data_current2[i]=[float(j) for j in data_current2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 2 converted to float")

#%%
# NOW WE GUESS THE BINS IN EACH TABLE (WE JUST HAVE TO GUESS BOTH THE BIN ON Q)

## BINS ON Q:
    
def Qbounds1(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current1,i,3)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("Q",data_current1,i,3)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("Q",data_current1,i,3)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


def Qbounds2(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current2,i,3)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("Q",data_current2,i,3)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("Q",data_current2,i,3)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

#%%
## IMPORT THE COMPUTED WW DATA FROM THE SAVED FILE TO ADD TO THE EXPERIMENTAL TABLES
print("Reading the Wandzura-Wilczek proton data file ...")

## Table 1:
f1 = open(path_to_g2WW_data+"E143_1996_1P.csv")
WW_data1=[]
for line in f1:    
    WW_data1.append(line.rstrip('\n'))

f1.close()


for i in range(len(WW_data1)):    
    WW_data1[i]=WW_data1[i].split(",")    
    WW_data1[i]=[float(j) for j in WW_data1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

## Table 2:
f2 = open(path_to_g2WW_data+"E143_1996_2P.csv")
WW_data2=[]
for line in f2:    
    WW_data2.append(line.rstrip('\n'))

f2.close()


for i in range(len(WW_data2)):    
    WW_data2[i]=WW_data2[i].split(",")    
    WW_data2[i]=[float(j) for j in WW_data2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")  

print("Done. Both tables of WW data converted to numbers")

#%%
## CONSTRUCT THE FINAL TABLES (WE BUILD 3 SEPARATE ONES AND MERGE THEM AT THE END TO UNITE THE WHOLE E155_PROTON EXPERIMENT)
    
print("Making G2 by E143_proton 1996 ... ")

## First we create a table which contains the data from all sub-tables: bins, data_current, and WW-data:
    
Qbin1= []
for i in range(len(data_current1)):
    Qbin1.append(Qbounds1(i))
    
Qbin2= []
for i in range(len(data_current2)):
    Qbin2.append(Qbounds2(i))

data_current=data_current1+data_current2
WW_data=WW_data1+WW_data2
Qbin_tot=Qbin1+Qbin2



#%%
DataCurrent=DataProcessor.DataSet.DataSet('E143-1995.p',"G2")
DataCurrent.comment="Taken from tables 1,2, 3. Bin in Q guessed by us."
DataCurrent.reference="10.17182/hepdata.19584.v1"

DataCurrent.isNormalized=False   ##Donde poner esta renormalización de la incertidumbre??
lumUncertainty=0.0
s_current=(2.*29.1*M_proton+(M_proton)**2)

for i in range(len(data_current)):
    # make up a point
    p = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p["process"]=100
    p["s"]=s_current  
    p["<Q>"]=numpy.sqrt(data_current[i][3])     
    p["Q"]=[Qbin_tot[i][0],Qbin_tot[i][1]]
    p["<x>"]=data_current[i][0]  
    p["x"]=[data_current[i][1],data_current[i][2]]
    p["xSec"]=data_current[i][4]-WW_data[i][0]
    ## THE UNCERTAINTY IS ORDERED AS: STAT(G2),SYST(G2),STAT(WW)
    p["uncorrErr"]=[(data_current[i][5]-data_current[i][6])/2,(data_current[i][7]-data_current[i][8])/2,WW_data[i][1]]
    p["thFactor"]=1.        
    #
    DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")


#%%


################################################################
#########----------------------------------------------#########
######--------------E143 1996 deuteron data --------------######
#########----------------------------------------------#########
################################################################

print("Reading E143 1996 deuteron data files ...")

f1 = open(path_to_data+"/E143/1996/Table3.csv")

data_from_f1=[]

for line in f1:    
    data_from_f1.append(line.rstrip('\n'))

f1.close()

print("Done.  =>     Convert table 1 to numbers ...")


data_current1=data_from_f1[22:29]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current1)):    
    data_current1[i]=data_current1[i].split(",")    
    data_current1[i]=[float(j) for j in data_current1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 1 converted to float")


f2 = open(path_to_data+"/E143/1996/Table4.csv")

data_from_f2=[]

for line in f2:    
    data_from_f2.append(line.rstrip('\n'))

f2.close()

print("Done.  =>     Convert table 2 to numbers ...")


data_current2=data_from_f2[20:25]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current2)):    
    data_current2[i]=data_current2[i].split(",")    
    data_current2[i]=[float(j) for j in data_current2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 2 converted to float")

#%%
# NOW WE GUESS THE BINS IN EACH TABLE (WE JUST HAVE TO GUESS BOTH THE BIN ON Q)

## BINS ON Q:
    
def Qbounds1(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current1,i,3)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("Q",data_current1,i,3)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("Q",data_current1,i,3)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


def Qbounds2(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current2,i,3)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("Q",data_current2,i,3)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("Q",data_current2,i,3)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

#%%
## IMPORT THE COMPUTED WW DATA FROM THE SAVED FILE TO ADD TO THE EXPERIMENTAL TABLES
print("Reading the Wandzura-Wilczek proton data file ...")

## Table 1:
f1 = open(path_to_g2WW_data+"E143_1996_1D.csv")
WW_data1=[]
for line in f1:    
    WW_data1.append(line.rstrip('\n'))

f1.close()


for i in range(len(WW_data1)):    
    WW_data1[i]=WW_data1[i].split(",")    
    WW_data1[i]=[float(j) for j in WW_data1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

## Table 2:
f2 = open(path_to_g2WW_data+"E143_1996_2D.csv")
WW_data2=[]
for line in f2:    
    WW_data2.append(line.rstrip('\n'))

f2.close()


for i in range(len(WW_data2)):    
    WW_data2[i]=WW_data2[i].split(",")    
    WW_data2[i]=[float(j) for j in WW_data2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")  

print("Done. Both tables of WW data converted to numbers")

#%%
## CONSTRUCT THE FINAL TABLES (WE BUILD 3 SEPARATE ONES AND MERGE THEM AT THE END TO UNITE THE WHOLE E155_PROTON EXPERIMENT)
    
print("Making G2 by E143_deuteron 1996 ... ")

## First we create a table which contains the data from all sub-tables: bins, data_current, and WW-data:
    
Qbin1= []
for i in range(len(data_current1)):
    Qbin1.append(Qbounds1(i))
    
Qbin2= []
for i in range(len(data_current2)):
    Qbin2.append(Qbounds2(i))

data_current=data_current1+data_current2
WW_data=WW_data1+WW_data2
Qbin_tot=Qbin1+Qbin2



#%%
DataCurrent=DataProcessor.DataSet.DataSet('E143-1995.d',"G2")
DataCurrent.comment="Taken from tables 1,2, 3. Bin in Q guessed by us."
DataCurrent.reference="10.17182/hepdata.19584.v1"

DataCurrent.isNormalized=False   ##Donde poner esta renormalización de la incertidumbre??
lumUncertainty=0.0
s_current=(2.*29.1*M_deuteron+(M_deuteron)**2)

for i in range(len(data_current)):
    # make up a point
    p = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p["process"]=102
    p["s"]=s_current  
    p["<Q>"]=numpy.sqrt(data_current[i][3])     
    p["Q"]=[Qbin_tot[i][0],Qbin_tot[i][1]]
    p["<x>"]=data_current[i][0]  
    p["x"]=[data_current[i][1],data_current[i][2]]
    p["xSec"]=data_current[i][4]-WW_data[i][0]
    ## THE UNCERTAINTY IS ORDERED AS: STAT(G2),SYST(G2),STAT(WW)
    p["uncorrErr"]=[(data_current[i][5]-data_current[i][6])/2,(data_current[i][7]-data_current[i][8])/2,WW_data[i][1]]
    p["thFactor"]=1.        
    #
    DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")

#%%
##############################################################
#########--------------------------------------------#########
######--------------E142 1996 neutron data--------------######
#########--------------------------------------------#########
##############################################################

print("Reading E142 1996 Neutron data files ...")

f1 = open(path_to_data+"/E142/Table4.csv")

data_from_f1=[]

for line in f1:    
    data_from_f1.append(line.rstrip('\n'))

f1.close()

print("Done.  =>     Convert table 1 to numbers ...")


data_current1=data_from_f1[25:33]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current1)):    
    data_current1[i]=data_current1[i].split(",")    
    data_current1[i]=[float(j) for j in data_current1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 1 converted to float")

#%%
# NOW WE GUESS THE BINS IN EACH TABLE (WE JUST HAVE TO GUESS BOTH THE BIN ON Q)

## BINS ON Q:
    
def Qbounds1(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current1,i,3)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("Q",data_current1,i,3)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("Q",data_current1,i,3)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

#%%
## IMPORT THE COMPUTED WW DATA FROM THE SAVED FILE TO ADD TO THE EXPERIMENTAL TABLES
print("Reading the Wandzura-Wilczek proton data file ...")

## Table 1:
f1 = open(path_to_g2WW_data+"E142_N.csv")
WW_data1=[]
for line in f1:    
    WW_data1.append(line.rstrip('\n'))

f1.close()


for i in range(len(WW_data1)):    
    WW_data1[i]=WW_data1[i].split(",")    
    WW_data1[i]=[float(j) for j in WW_data1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")
    
print("Done. The table of WW data has been converted to numbers")

#%%
## CONSTRUCT THE FINAL TABLES (WE BUILD 3 SEPARATE ONES AND MERGE THEM AT THE END TO UNITE THE WHOLE E155_PROTON EXPERIMENT)
    
print("Making G2 by E142_Neutron ... ")

## First we create a table which contains the data from all sub-tables: bins, data_current, and WW-data:
    
Qbin1= []
for i in range(len(data_current1)):
    Qbin1.append(Qbounds1(i))


#%%
DataCurrent=DataProcessor.DataSet.DataSet('E142.n',"G2")
DataCurrent.comment="Taken from table 4. Bin in Q guessed by us. We take the three incident energies presented for the data and assume the values computed are for the average energy from those 3"
DataCurrent.reference="hep-ex/9610007"

DataCurrent.isNormalized=False   ##Donde poner esta renormalización de la incertidumbre??
lumUncertainty=0.0
s_current=(2.*22.367*M_neutron+(M_neutron)**2)

for i in range(len(data_current1)):
    # make up a point
    p = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p["process"]=101
    p["s"]=s_current  
    p["<Q>"]=numpy.sqrt(data_current1[i][3])     
    p["Q"]=[Qbin1[i][0],Qbin1[i][1]]
    p["<x>"]=data_current1[i][0]  
    p["x"]=[data_current1[i][1],data_current1[i][2]]
    p["xSec"]=data_current1[i][4]-WW_data1[i][0]
    ## THE UNCERTAINTY IS ORDERED AS: STAT(G2),SYST(G2),STAT(WW)
    p["uncorrErr"]=[(data_current1[i][5]-data_current1[i][6])/2,(data_current1[i][7]-data_current1[i][8])/2,WW_data1[i][1]]
    p["thFactor"]=1.        
    #
    DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")

#%%
#################################################################
#########-----------------------------------------------#########
######-------------- Spin Muon Collaboration --------------######
#########-----------------------------------------------#########
#################################################################

print("Reading SMC proton data files ...")

f1 = open(path_to_data+"/SMC/Table8.csv")

data_from_f1=[]

for line in f1:    
    data_from_f1.append(line.rstrip('\n'))

f1.close()

print("Done.  =>     Convert table 1 to numbers ...")


data_current1=data_from_f1[13:19]  

for i in range(len(data_current1)):    
    data_current1[i]=data_current1[i].split(",")    
    data_current1[i]=[float(j) for j in data_current1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 1 converted to float")


f2 = open(path_to_data+"/SMC/Table9.csv")

data_from_f2=[]

for line in f2:    
    data_from_f2.append(line.rstrip('\n'))

f2.close()

print("Done.  =>     Convert table 2 to numbers ...")


data_current2=data_from_f2[13:19]  

for i in range(len(data_current2)):    
    data_current2[i]=data_current2[i].split(",")    
    data_current2[i]=[float(j) for j in data_current2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 2 converted to float")

#%%
## BIN ON Q:
    
def Qbounds1(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current1,i,3)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("Q",data_current1,i,3)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("Q",data_current1,i,3)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


#%%
## CONSTRUCT THE FINAL TABLES (WE BUILD 3 SEPARATE ONES AND MERGE THEM AT THE END TO UNITE THE WHOLE E155_PROTON EXPERIMENT)
    
print("Making G2 by SMC proton...")

## First we create a table which contains the data from all sub-tables: bins, data_current, and WW-data:
    
Qbin1= []
for i in range(len(data_current1)):
    Qbin1.append(Qbounds1(i))
    
#%%
DataCurrent=DataProcessor.DataSet.DataSet('SMC.p',"G2")
DataCurrent.comment="Taken from tables 8 and 9. Bin in Q guessed by us. They neglect the systematic uncertainty on g2 and setting it to 0. Furthermore, they compute the WW-terms and add their uncertainty."
DataCurrent.reference="arXiv:hep-ex/9702005"

DataCurrent.isNormalized=False   
lumUncertainty=0.0
s_current=(2.*190*M_proton+(M_proton)**2)

for i in range(len(data_current1)):
    # make up a point
    p = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p["process"]=100
    p["s"]=s_current  
    p["<Q>"]=numpy.sqrt(data_current1[i][3])     
    p["Q"]=[Qbin1[i][0],Qbin1[i][1]]
    p["<x>"]=data_current1[i][0]  
    p["x"]=[data_current1[i][1],data_current1[i][2]]
    p["xSec"]=data_current1[i][5]-data_current2[i][5]
    ## THE UNCERTAINTY IS ORDERED AS: STAT(G2),SYST(G2),STAT(WW)
    p["uncorrErr"]=[(data_current1[i][6]-data_current1[i][7])/2,(data_current2[i][6]-data_current2[i][7])/2]
    p["thFactor"]=1.        
    #
    DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")

#%%

#######################################################
#########-------------------------------------#########
######-------------- RSSC (proton) --------------######
#########-------------------------------------#########
#######################################################

print("Reading RSCC proton data file...")

f1 = open(path_to_data+"/RSSC/2007Exp_P.csv")

data_from_f1=[]

for line in f1:    
    data_from_f1.append(line.rstrip('\n'))

f1.close()

print("Done.  =>     Convert table 1 to numbers ...")


data_current1=data_from_f1[0:26]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current1)):    
    data_current1[i]=data_current1[i].split(",")    
    data_current1[i]=[float(j) for j in data_current1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 1 converted to float")

#%%
        
## BIN ON x:
    
def xbounds1(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current1,i,0)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("x",data_current1,i,0)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("x",data_current1,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")
        
#%%
## IMPORT THE COMPUTED WW DATA FROM THE SAVED FILE TO ADD TO THE EXPERIMENTAL TABLES
print("Reading the Wandzura-Wilczek proton data file ...")

## Table 1:
f1 = open(path_to_g2WW_data+"RSSC.csv")
WW_data1=[]
for line in f1:    
    WW_data1.append(line.rstrip('\n'))

f1.close()


for i in range(len(WW_data1)):    
    WW_data1[i]=WW_data1[i].split(",")    
    WW_data1[i]=[float(j) for j in WW_data1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")
    
print("Done. The table of WW data has been converted to numbers")

#%%
## CONSTRUCT THE FINAL TABLES (WE BUILD 3 SEPARATE ONES AND MERGE THEM AT THE END TO UNITE THE WHOLE E155_PROTON EXPERIMENT)
    
print("Making G2 by RSSC (proton)... ")

## First we create a table which contains the data from all sub-tables: bins, data_current, and WW-data:
    
xbin1= []
for i in range(len(data_current1)):
    xbin1.append(xbounds1(i))

#%%
DataCurrent=DataProcessor.DataSet.DataSet('RSS.p',"G2")
DataCurrent.comment="Extracted from FIG 4 (in the paper) with the aid of ''Plot Digitalizer'' because there are no Datasets. Bin in x guessed by us. All measurements at the same energy of Q²1.36 GeV², we assume there's no bin on Q then."
DataCurrent.reference="arXiv:nucl-ex/0608003v3"

DataCurrent.isNormalized=False   ##Donde poner esta renormalización de la incertidumbre??
lumUncertainty=0.0
s_current=(2.*5.755*M_proton+(M_proton)**2)

for i in range(len(data_current1)):
    # make up a point
    p = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p["process"]=100
    p["s"]=s_current  
    p["<Q>"]=numpy.sqrt(1.3)     
    p["Q"]=[numpy.sqrt(0.8),numpy.sqrt(1.4)]
    p["<x>"]=data_current1[i][0]  
    p["x"]=[xbin1[i][0],xbin1[i][1]]
    p["xSec"]=data_current1[i][1]-WW_data1[i][0]
    ## THE UNCERTAINTY IS ORDERED AS: STAT(G2),SYST(G2),STAT(WW)
    p["uncorrErr"]=[data_current1[i][2],WW_data1[i][1]]
    p["thFactor"]=1.        
    #
    DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")



#%%
########################################################################
#########------------------------------------------------------#########
######--------------E155 2003 proton data (29.1 GeV)--------------######
#########------------------------------------------------------#########
########################################################################

print("Reading E155 2003 proton data files ...")

f1 = open(path_to_data+"/E155/2003/Table1.csv")

data_from_f1=[]

for line in f1:    
    data_from_f1.append(line.rstrip('\n'))

f1.close()

print("Done.  =>     Convert table 1 to numbers ...")


data_current1=data_from_f1[29:37]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current1)):    
    data_current1[i]=data_current1[i].split(",")    
    data_current1[i]=[float(j) for j in data_current1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 1 converted to float")


f2 = open(path_to_data+"/E155/2003/Table2.csv")

data_from_f2=[]

for line in f2:    
    data_from_f2.append(line.rstrip('\n'))

f2.close()

print("Done.  =>     Convert table 2 to numbers ...")


data_current2=data_from_f2[28:35]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current2)):    
    data_current2[i]=data_current2[i].split(",")    
    data_current2[i]=[float(j) for j in data_current2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 2 converted to float")


f3 = open(path_to_data+"/E155/2003/Table3.csv")

data_from_f3=[]

for line in f3:    
    data_from_f3.append(line.rstrip('\n'))

f3.close()

print("Done.  =>     Convert table 3 to numbers ...")


data_current3=data_from_f3[26:31]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current3)):    
    data_current3[i]=data_current3[i].split(",")    
    data_current3[i]=[float(j) for j in data_current3[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 3 converted to float")


#%%
# NOW WE GUESS THE BINS IN EACH TABLE (WE HAVE TO GUESS BOTH THE BIN ON x AND Q)

## BINS ON Q:
    
def Qbounds1(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current1,i,1)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("Q",data_current1,i,1)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("Q",data_current1,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


def Qbounds2(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current2,i,1)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("Q",data_current2,i,1)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("Q",data_current2,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")
        
        
def Qbounds3(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current3,i,1)
    
    elif 0<i<len(data_current3)-1:
        return GuessBin_Middle("Q",data_current3,i,1)
    
    elif i== len(data_current3)-1:
       return GuessBin_Top("Q",data_current3,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

## BINS ON x:
    
def xbounds1(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current1,i,0)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("x",data_current1,i,0)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("x",data_current1,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


def xbounds2(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current2,i,0)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("x",data_current2,i,0)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("x",data_current2,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")
        
        
def xbounds3(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current3,i,0)
    
    elif 0<i<len(data_current3)-1:
        return GuessBin_Middle("x",data_current3,i,0)
    
    elif i== len(data_current3)-1:
       return GuessBin_Top("x",data_current3,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

#%%
## IMPORT THE COMPUTED WW DATA FROM THE SAVED FILE TO ADD TO THE EXPERIMENTAL TABLES
print("Reading the Wandzura-Wilczek proton data file ...")

## Table 1:
    
f1 = open(path_to_g2WW_data+"E155_2003_1P(29.1GeV).csv")

WW_data1=[]

for line in f1:    
    WW_data1.append(line.rstrip('\n'))

f1.close()


for i in range(len(WW_data1)):    
    WW_data1[i]=WW_data1[i].split(",")    
    WW_data1[i]=[float(j) for j in WW_data1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

## Table 2:

f2 = open(path_to_g2WW_data+"E155_2003_2P(29.1GeV).csv")

WW_data2=[]

for line in f2:    
    WW_data2.append(line.rstrip('\n'))

f2.close()


for i in range(len(WW_data2)):    
    WW_data2[i]=WW_data2[i].split(",")    
    WW_data2[i]=[float(j) for j in WW_data2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

## Table 3:
    
f3 = open(path_to_g2WW_data+"E155_2003_3P(29.1GeV).csv")

WW_data3=[]

for line in f3:    
    WW_data3.append(line.rstrip('\n'))

f3.close()

print("Done.  =>     Convert to numbers ...")

for i in range(len(WW_data3)):    
    WW_data3[i]=WW_data3[i].split(",")    
    WW_data3[i]=[float(j) for j in WW_data3[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")    


print("Done. All 3 tables of WW data converted to numbers")

#%%
## CONSTRUCT THE FINAL TABLES (WE BUILD 3 SEPARATE ONES AND MERGE THEM AT THE END TO UNITE THE WHOLE E155_PROTON EXPERIMENT)
## Note: THE SYSTEMATIC UNCERTAINTY FOR g2 IS GIVEN BY THE FOLLOWING FORMULA (paper, page 8):
def syst_p(x):
    return 0.0016-x*0.0012
    
print("Making G2 by E155_proton 2003 (29.1 GeV)... ")

## First we create a table which contains the data from all 3 sub-tables: bins, data_current, and WW-data:

xbin1= []
for i in range(len(data_current1)):
    xbin1.append(xbounds1(i))

xbin2= []
for i in range(len(data_current2)):
    xbin2.append(xbounds2(i))

xbin3= []
for i in range(len(data_current3)):
    xbin3.append(xbounds3(i))
    
Qbin1= []
for i in range(len(data_current1)):
    Qbin1.append(Qbounds1(i))
    
Qbin2= []
for i in range(len(data_current2)):
    Qbin2.append(Qbounds2(i))
    
Qbin3= []
for i in range(len(data_current3)):
    Qbin3.append(Qbounds3(i))

data_current=data_current1+data_current2+data_current3
WW_data=WW_data1+WW_data2+WW_data3
xbin_tot=xbin1+xbin2+xbin3
Qbin_tot=Qbin1+Qbin2+Qbin3



#%%
DataCurrent=DataProcessor.DataSet.DataSet('E155-29.p',"G2")
DataCurrent.comment="Taken from tables 1,2, 3. Bin in Q and x guessed by us."
DataCurrent.reference="arXiv:hep-ex/0204028v1"

DataCurrent.isNormalized=False   ##Donde poner esta renormalización de la incertidumbre??
lumUncertainty=0.0
s_current=(2.*29.1*M_proton+(M_proton)**2)

for i in range(len(data_current)):
    # make up a point
    p = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p["process"]=100
    p["s"]=s_current  
    p["<Q>"]=numpy.sqrt(data_current[i][1])     
    p["Q"]=[Qbin_tot[i][0],Qbin_tot[i][1]]
    p["<x>"]=data_current[i][0]  
    p["x"]=[xbin_tot[i][0],xbin_tot[i][1]]  
    p["xSec"]=data_current[i][2]-data_current[i][0]*WW_data[i][0]
    ## THE UNCERTAINTY IS ORDERED AS: STAT(G2),SYST(G2),STAT(WW)
    p["uncorrErr"]=[(data_current[i][3]-data_current[i][4])/2,syst_p(data_current[i][0]),data_current[i][0]*WW_data[i][1]]
    p["thFactor"]=data_current[i][0]
    #
    DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")


#%%
##########################################################################
#########--------------------------------------------------------#########
######--------------E155 2003 deuteron data (29.1 GeV)--------------######
#########--------------------------------------------------------#########
##########################################################################

print("Reading E155 2003 deuteron data files ...")

f1 = open(path_to_data+"/E155/2003/Table1.csv")

data_from_f1=[]

for line in f1:    
    data_from_f1.append(line.rstrip('\n'))

f1.close()

print("Done.  =>     Convert table 1 to numbers ...")


data_current1=data_from_f1[56:64]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current1)):    
    data_current1[i]=data_current1[i].split(",")    
    data_current1[i]=[float(j) for j in data_current1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 1 converted to float")


f2 = open(path_to_data+"/E155/2003/Table2.csv")

data_from_f2=[]

for line in f2:    
    data_from_f2.append(line.rstrip('\n'))

f2.close()

print("Done.  =>     Convert table 2 to numbers ...")


data_current2=data_from_f2[53:60]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current2)):    
    data_current2[i]=data_current2[i].split(",")    
    data_current2[i]=[float(j) for j in data_current2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 2 converted to float")


f3 = open(path_to_data+"/E155/2003/Table3.csv")

data_from_f3=[]

for line in f3:    
    data_from_f3.append(line.rstrip('\n'))

f3.close()

print("Done.  =>     Convert table 3 to numbers ...")


data_current3=data_from_f3[47:52]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current3)):    
    data_current3[i]=data_current3[i].split(",")    
    data_current3[i]=[float(j) for j in data_current3[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 3 converted to float")


#%%
# NOW WE GUESS THE BINS IN EACH TABLE (WE HAVE TO GUESS BOTH THE BIN ON x AND Q)

## BINS ON Q:
    
def Qbounds1(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current1,i,1)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("Q",data_current1,i,1)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("Q",data_current1,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


def Qbounds2(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current2,i,1)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("Q",data_current2,i,1)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("Q",data_current2,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")
        
        
def Qbounds3(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current3,i,1)
    
    elif 0<i<len(data_current3)-1:
        return GuessBin_Middle("Q",data_current3,i,1)
    
    elif i== len(data_current3)-1:
       return GuessBin_Top("Q",data_current3,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

## BINS ON x:
    
def xbounds1(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current1,i,0)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("x",data_current1,i,0)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("x",data_current1,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


def xbounds2(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current2,i,0)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("x",data_current2,i,0)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("x",data_current2,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")
        
        
def xbounds3(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current3,i,0)
    
    elif 0<i<len(data_current3)-1:
        return GuessBin_Middle("x",data_current3,i,0)
    
    elif i== len(data_current3)-1:
       return GuessBin_Top("x",data_current3,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

#%%
## IMPORT THE COMPUTED WW DATA FROM THE SAVED FILE TO ADD TO THE EXPERIMENTAL TABLES
print("Reading the Wandzura-Wilczek proton data file ...")

## Table 1:
    
f1 = open(path_to_g2WW_data+"E155_2003_1D(29.1GeV).csv")

WW_data1=[]

for line in f1:    
    WW_data1.append(line.rstrip('\n'))

f1.close()


for i in range(len(WW_data1)):    
    WW_data1[i]=WW_data1[i].split(",")    
    WW_data1[i]=[float(j) for j in WW_data1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

## Table 2:

f2 = open(path_to_g2WW_data+"E155_2003_2D(29.1GeV).csv")

WW_data2=[]

for line in f2:    
    WW_data2.append(line.rstrip('\n'))

f2.close()


for i in range(len(WW_data2)):    
    WW_data2[i]=WW_data2[i].split(",")    
    WW_data2[i]=[float(j) for j in WW_data2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

## Table 3:
    
f3 = open(path_to_g2WW_data+"E155_2003_3D(29.1GeV).csv")

WW_data3=[]

for line in f3:    
    WW_data3.append(line.rstrip('\n'))

f3.close()

print("Done.  =>     Convert to numbers ...")

for i in range(len(WW_data3)):    
    WW_data3[i]=WW_data3[i].split(",")    
    WW_data3[i]=[float(j) for j in WW_data3[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")    


print("Done. All 3 tables of WW data converted to numbers")

#%%
## CONSTRUCT THE FINAL TABLES (WE BUILD 3 SEPARATE ONES AND MERGE THEM AT THE END TO UNITE THE WHOLE E155_PROTON EXPERIMENT)
## Note: THE SYSTEMATIC UNCERTAINTY FOR g2 IS GIVEN BY THE FOLLOWING FORMULA (paper, page 8):
def syst_d(x):
    return 0.0009-x*0.0009
    
print("Making G2 by E155_deuteron 2003 (29.1 GeV)... ")

## First we create a table which contains the data from all 3 sub-tables: bins, data_current, and WW-data:

xbin1= []
for i in range(len(data_current1)):
    xbin1.append(xbounds1(i))

xbin2= []
for i in range(len(data_current2)):
    xbin2.append(xbounds2(i))

xbin3= []
for i in range(len(data_current3)):
    xbin3.append(xbounds3(i))
    
Qbin1= []
for i in range(len(data_current1)):
    Qbin1.append(Qbounds1(i))
    
Qbin2= []
for i in range(len(data_current2)):
    Qbin2.append(Qbounds2(i))
    
Qbin3= []
for i in range(len(data_current3)):
    Qbin3.append(Qbounds3(i))

data_current=data_current1+data_current2+data_current3
WW_data=WW_data1+WW_data2+WW_data3
xbin_tot=xbin1+xbin2+xbin3
Qbin_tot=Qbin1+Qbin2+Qbin3



#%%
DataCurrent=DataProcessor.DataSet.DataSet('E155-29.d',"G2")
DataCurrent.comment="Taken from tables 1,2, 3. Bin in Q and x guessed by us."
DataCurrent.reference="arXiv:hep-ex/0204028v1"

DataCurrent.isNormalized=False   ##Donde poner esta renormalización de la incertidumbre??
lumUncertainty=0.0
s_current=(2.*29.1*M_deuteron+(M_deuteron)**2)

for i in range(len(data_current)):
    # make up a point
    p = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p["process"]=102
    p["s"]=s_current  
    p["<Q>"]=numpy.sqrt(data_current[i][1])     
    p["Q"]=[Qbin_tot[i][0],Qbin_tot[i][1]]
    p["<x>"]=data_current[i][0]  
    p["x"]=[xbin_tot[i][0],xbin_tot[i][1]]  
    p["xSec"]=data_current[i][2]-data_current[i][0]*WW_data[i][0]
    ## THE UNCERTAINTY IS ORDERED AS: STAT(G2),SYST(G2),STAT(WW)
    p["uncorrErr"]=[(data_current[i][3]-data_current[i][4])/2,syst_d(data_current[i][0]),data_current[i][0]*WW_data[i][1]]
    p["thFactor"]=data_current[i][0]        
    #
    DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")

#%%
########################################################################
#########------------------------------------------------------#########
######--------------E155 2003 proton data (32.3 GeV)--------------######
#########------------------------------------------------------#########
########################################################################

print("Reading E155 2003 proton data files ...")

f1 = open(path_to_data+"/E155/2003/Table4.csv")

data_from_f1=[]

for line in f1:    
    data_from_f1.append(line.rstrip('\n'))

f1.close()

print("Done.  =>     Convert table 1 to numbers ...")


data_current1=data_from_f1[29:37]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current1)):    
    data_current1[i]=data_current1[i].split(",")    
    data_current1[i]=[float(j) for j in data_current1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 1 converted to float")

f2 = open(path_to_data+"/E155/2003/Table5.csv")

data_from_f2=[]

for line in f2:    
    data_from_f2.append(line.rstrip('\n'))

f2.close()

print("Done.  =>     Convert table 2 to numbers ...")


data_current2=data_from_f2[28:35]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current2)):    
    data_current2[i]=data_current2[i].split(",")    
    data_current2[i]=[float(j) for j in data_current2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 2 converted to float")


f3 = open(path_to_data+"/E155/2003/Table6.csv")

data_from_f3=[]

for line in f3:    
    data_from_f3.append(line.rstrip('\n'))

f3.close()

print("Done.  =>     Convert table 3 to numbers ...")


data_current3=data_from_f3[26:31]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current3)):    
    data_current3[i]=data_current3[i].split(",")    
    data_current3[i]=[float(j) for j in data_current3[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 3 converted to float")


#%%
# NOW WE GUESS THE BINS IN EACH TABLE (WE HAVE TO GUESS BOTH THE BIN ON x AND Q)

## BINS ON Q:
    
def Qbounds1(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current1,i,1)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("Q",data_current1,i,1)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("Q",data_current1,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


def Qbounds2(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current2,i,1)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("Q",data_current2,i,1)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("Q",data_current2,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")
        
        
def Qbounds3(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current3,i,1)
    
    elif 0<i<len(data_current3)-1:
        return GuessBin_Middle("Q",data_current3,i,1)
    
    elif i== len(data_current3)-1:
       return GuessBin_Top("Q",data_current3,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

## BINS ON x:
    
def xbounds1(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current1,i,0)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("x",data_current1,i,0)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("x",data_current1,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


def xbounds2(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current2,i,0)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("x",data_current2,i,0)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("x",data_current2,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")
        
        
def xbounds3(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current3,i,0)
    
    elif 0<i<len(data_current3)-1:
        return GuessBin_Middle("x",data_current3,i,0)
    
    elif i== len(data_current3)-1:
       return GuessBin_Top("x",data_current3,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

#%%
## IMPORT THE COMPUTED WW DATA FROM THE SAVED FILE TO ADD TO THE EXPERIMENTAL TABLES
print("Reading the Wandzura-Wilczek proton data file ...")

## IMPORTANT: THE POINTS FOR x AND Q² ARE THE SAME FOR BOTH 29.1 and 32.3 GeV, THIS MEANS THAT WE ONLY NEED THE WW-TERMS FROM THE PREVIOUS
## COMPUTATION, THAT'S WHY WE LOAD THEM INSTEAD OF COMPUTING THEM AGAIN
## Table 1:
    
f1 = open(path_to_g2WW_data+"E155_2003_1P(29.1GeV).csv")

WW_data1=[]

for line in f1:    
    WW_data1.append(line.rstrip('\n'))

f1.close()


for i in range(len(WW_data1)):    
    WW_data1[i]=WW_data1[i].split(",")    
    WW_data1[i]=[float(j) for j in WW_data1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

## Table 2:

f2 = open(path_to_g2WW_data+"E155_2003_2P(29.1GeV).csv")

WW_data2=[]

for line in f2:    
    WW_data2.append(line.rstrip('\n'))

f2.close()


for i in range(len(WW_data2)):    
    WW_data2[i]=WW_data2[i].split(",")    
    WW_data2[i]=[float(j) for j in WW_data2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

## Table 3:
    
f3 = open(path_to_g2WW_data+"E155_2003_3P(29.1GeV).csv")

WW_data3=[]

for line in f3:    
    WW_data3.append(line.rstrip('\n'))

f3.close()

print("Done.  =>     Convert to numbers ...")

for i in range(len(WW_data3)):    
    WW_data3[i]=WW_data3[i].split(",")    
    WW_data3[i]=[float(j) for j in WW_data3[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")    


print("Done. All 3 tables of WW data converted to numbers")

#%%
## CONSTRUCT THE FINAL TABLES (WE BUILD 3 SEPARATE ONES AND MERGE THEM AT THE END TO UNITE THE WHOLE E155_PROTON EXPERIMENT)
## Note: THE SYSTEMATIC UNCERTAINTY FOR g2 IS GIVEN BY THE FOLLOWING FORMULA (paper, page 8):
def syst_p(x):
    return 0.0016-x*0.0012
    
print("Making G2 by E155_proton 2003 (32.3 GeV)... ")

## First we create a table which contains the data from all 3 sub-tables: bins, data_current, and WW-data:

xbin1= []
for i in range(len(data_current1)):
    xbin1.append(xbounds1(i))

xbin2= []
for i in range(len(data_current2)):
    xbin2.append(xbounds2(i))

xbin3= []
for i in range(len(data_current3)):
    xbin3.append(xbounds3(i))
    
Qbin1= []
for i in range(len(data_current1)):
    Qbin1.append(Qbounds1(i))
    
Qbin2= []
for i in range(len(data_current2)):
    Qbin2.append(Qbounds2(i))
    
Qbin3= []
for i in range(len(data_current3)):
    Qbin3.append(Qbounds3(i))

data_current=data_current1+data_current2+data_current3
WW_data=WW_data1+WW_data2+WW_data3
xbin_tot=xbin1+xbin2+xbin3
Qbin_tot=Qbin1+Qbin2+Qbin3



#%%
DataCurrent=DataProcessor.DataSet.DataSet('E155-32.p',"G2")
DataCurrent.comment="Taken from tables 4, 5, 6. Bin in Q and x guessed by us."
DataCurrent.reference="arXiv:hep-ex/0204028v1"

DataCurrent.isNormalized=False   ##Donde poner esta renormalización de la incertidumbre??
lumUncertainty=0.0
s_current=(2.*32.3*M_proton+(M_proton)**2)

for i in range(len(data_current)):
    # make up a point
    p = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p["process"]=100
    p["s"]=s_current  
    p["<Q>"]=numpy.sqrt(data_current[i][1])     
    p["Q"]=[Qbin_tot[i][0],Qbin_tot[i][1]]
    p["<x>"]=data_current[i][0]  
    p["x"]=[xbin_tot[i][0],xbin_tot[i][1]]  
    p["xSec"]=data_current[i][2]-data_current[i][0]*WW_data[i][0]
    ## THE UNCERTAINTY IS ORDERED AS: STAT(G2),SYST(G2),STAT(WW)
    p["uncorrErr"]=[(data_current[i][3]-data_current[i][4])/2,syst_p(data_current[i][0]),data_current[i][0]*WW_data[i][1]]
    p["thFactor"]=data_current[i][0]      
    #
    DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")


#%%
##########################################################################
#########--------------------------------------------------------#########
######--------------E155 2003 deuteron data (32.3 GeV)--------------######
#########--------------------------------------------------------#########
##########################################################################

print("Reading E155 2003 deuteron data files ...")

f1 = open(path_to_data+"/E155/2003/Table4.csv")

data_from_f1=[]

for line in f1:    
    data_from_f1.append(line.rstrip('\n'))

f1.close()

print("Done.  =>     Convert table 1 to numbers ...")


data_current1=data_from_f1[56:64]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current1)):    
    data_current1[i]=data_current1[i].split(",")    
    data_current1[i]=[float(j) for j in data_current1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 1 converted to float")


f2 = open(path_to_data+"/E155/2003/Table5.csv")

data_from_f2=[]

for line in f2:    
    data_from_f2.append(line.rstrip('\n'))

f2.close()

print("Done.  =>     Convert table 2 to numbers ...")


data_current2=data_from_f2[53:60]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current2)):    
    data_current2[i]=data_current2[i].split(",")    
    data_current2[i]=[float(j) for j in data_current2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 2 converted to float")


f3 = open(path_to_data+"/E155/2003/Table6.csv")

data_from_f3=[]

for line in f3:    
    data_from_f3.append(line.rstrip('\n'))

f3.close()

print("Done.  =>     Convert table 3 to numbers ...")


data_current3=data_from_f3[47:52]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current3)):    
    data_current3[i]=data_current3[i].split(",")    
    data_current3[i]=[float(j) for j in data_current3[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 3 converted to float")


#%%
# NOW WE GUESS THE BINS IN EACH TABLE (WE HAVE TO GUESS BOTH THE BIN ON x AND Q)

## BINS ON Q:
    
def Qbounds1(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current1,i,1)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("Q",data_current1,i,1)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("Q",data_current1,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


def Qbounds2(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current2,i,1)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("Q",data_current2,i,1)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("Q",data_current2,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")
        
        
def Qbounds3(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current3,i,1)
    
    elif 0<i<len(data_current3)-1:
        return GuessBin_Middle("Q",data_current3,i,1)
    
    elif i== len(data_current3)-1:
       return GuessBin_Top("Q",data_current3,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

## BINS ON x:
    
def xbounds1(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current1,i,0)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("x",data_current1,i,0)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("x",data_current1,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


def xbounds2(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current2,i,0)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("x",data_current2,i,0)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("x",data_current2,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")
        
        
def xbounds3(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current3,i,0)
    
    elif 0<i<len(data_current3)-1:
        return GuessBin_Middle("x",data_current3,i,0)
    
    elif i== len(data_current3)-1:
       return GuessBin_Top("x",data_current3,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

#%%
## IMPORT THE COMPUTED WW DATA FROM THE SAVED FILE TO ADD TO THE EXPERIMENTAL TABLES
## IMPORTANT: THE POINTS FOR x AND Q² ARE THE SAME FOR BOTH 29.1 and 32.3 GeV, THIS MEANS THAT WE ONLY NEED THE WW-TERMS FROM THE PREVIOUS
## COMPUTATION, THAT'S WHY WE LOAD THEM INSTEAD OF COMPUTING THEM AGAIN
print("Reading the Wandzura-Wilczek proton data file ...")

## Table 1:
    
f1 = open(path_to_g2WW_data+"E155_2003_1D(29.1GeV).csv")

WW_data1=[]

for line in f1:    
    WW_data1.append(line.rstrip('\n'))

f1.close()


for i in range(len(WW_data1)):    
    WW_data1[i]=WW_data1[i].split(",")    
    WW_data1[i]=[float(j) for j in WW_data1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

## Table 2:

f2 = open(path_to_g2WW_data+"E155_2003_2D(29.1GeV).csv")

WW_data2=[]

for line in f2:    
    WW_data2.append(line.rstrip('\n'))

f2.close()


for i in range(len(WW_data2)):    
    WW_data2[i]=WW_data2[i].split(",")    
    WW_data2[i]=[float(j) for j in WW_data2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

## Table 3:
    
f3 = open(path_to_g2WW_data+"E155_2003_3D(29.1GeV).csv")

WW_data3=[]

for line in f3:    
    WW_data3.append(line.rstrip('\n'))

f3.close()

print("Done.  =>     Convert to numbers ...")

for i in range(len(WW_data3)):    
    WW_data3[i]=WW_data3[i].split(",")    
    WW_data3[i]=[float(j) for j in WW_data3[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")    


print("Done. All 3 tables of WW data converted to numbers")

#%%
## CONSTRUCT THE FINAL TABLES (WE BUILD 3 SEPARATE ONES AND MERGE THEM AT THE END TO UNITE THE WHOLE E155_PROTON EXPERIMENT)
## Note: THE SYSTEMATIC UNCERTAINTY FOR g2 IS GIVEN BY THE FOLLOWING FORMULA (paper, page 8):
def syst_d(x):
    return 0.0009-x*0.0009
    
print("Making G2 by E155_deuteron 2003 (32.3 GeV)... ")

## First we create a table which contains the data from all 3 sub-tables: bins, data_current, and WW-data:

xbin1= []
for i in range(len(data_current1)):
    xbin1.append(xbounds1(i))

xbin2= []
for i in range(len(data_current2)):
    xbin2.append(xbounds2(i))

xbin3= []
for i in range(len(data_current3)):
    xbin3.append(xbounds3(i))
    
Qbin1= []
for i in range(len(data_current1)):
    Qbin1.append(Qbounds1(i))
    
Qbin2= []
for i in range(len(data_current2)):
    Qbin2.append(Qbounds2(i))
    
Qbin3= []
for i in range(len(data_current3)):
    Qbin3.append(Qbounds3(i))

data_current=data_current1+data_current2+data_current3
WW_data=WW_data1+WW_data2+WW_data3
xbin_tot=xbin1+xbin2+xbin3
Qbin_tot=Qbin1+Qbin2+Qbin3



#%%
DataCurrent=DataProcessor.DataSet.DataSet('E155-32.d',"G2")
DataCurrent.comment="Taken from tables 1,2, 3. Bin in Q and x guessed by us."
DataCurrent.reference="arXiv:hep-ex/0204028v1"

DataCurrent.isNormalized=False   ##Donde poner esta renormalización de la incertidumbre??
lumUncertainty=0.0
s_current=(2.*32.3*M_deuteron+(M_deuteron)**2)

for i in range(len(data_current)):
    # make up a point
    p = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p["process"]=102
    p["s"]=s_current  
    p["<Q>"]=numpy.sqrt(data_current[i][1])     
    p["Q"]=[Qbin_tot[i][0],Qbin_tot[i][1]]
    p["<x>"]=data_current[i][0]  
    p["x"]=[xbin_tot[i][0],xbin_tot[i][1]]  
    p["xSec"]=data_current[i][2]-data_current[i][0]*WW_data[i][0]
    ## THE UNCERTAINTY IS ORDERED AS: STAT(G2),SYST(G2),STAT(WW)
    p["uncorrErr"]=[(data_current[i][3]-data_current[i][4])/2,syst_d(data_current[i][0]),data_current[i][0]*WW_data[i][1]]
    p["thFactor"]=1.        
    #
    DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")


#%%
##############################################################
#########--------------------------------------------#########
######------------- E155 1999 proton data --------------######
#########--------------------------------------------#########
##############################################################

print("Reading E155 1999 proton data files ...")

f1 = open(path_to_data+"/E155/1999/Table1.csv")

data_from_f1=[]

for line in f1:    
    data_from_f1.append(line.rstrip('\n'))

f1.close()

print("Done.  =>     Convert table 1 to numbers ...")


data_current1=data_from_f1[30:40]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current1)):    
    data_current1[i]=data_current1[i].split(",")    
    data_current1[i]=[float(j) for j in data_current1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 1 converted to float")


f2 = open(path_to_data+"/E155/1999/Table2.csv")

data_from_f2=[]

for line in f2:    
    data_from_f2.append(line.rstrip('\n'))

f2.close()

print("Done.  =>     Convert table 2 to numbers ...")


data_current2=data_from_f2[27:34]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current2)):    
    data_current2[i]=data_current2[i].split(",")    
    data_current2[i]=[float(j) for j in data_current2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 2 converted to float")


f3 = open(path_to_data+"/E155/1999/Table3.csv")

data_from_f3=[]

for line in f3:    
    data_from_f3.append(line.rstrip('\n'))

f3.close()

print("Done.  =>     Convert table 3 to numbers ...")


data_current3=data_from_f3[25:30]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current3)):    
    data_current3[i]=data_current3[i].split(",")    
    data_current3[i]=[float(j) for j in data_current3[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 3 converted to float")


#%%
# NOW WE GUESS THE BINS IN EACH TABLE (WE HAVE TO GUESS BOTH THE BIN ON x AND Q)

## BINS ON Q:
    
def Qbounds1(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current1,i,1)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("Q",data_current1,i,1)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("Q",data_current1,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


def Qbounds2(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current2,i,1)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("Q",data_current2,i,1)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("Q",data_current2,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")
        
        
def Qbounds3(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current3,i,1)
    
    elif 0<i<len(data_current3)-1:
        return GuessBin_Middle("Q",data_current3,i,1)
    
    elif i== len(data_current3)-1:
       return GuessBin_Top("Q",data_current3,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

## BINS ON x:
    
def xbounds1(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current1,i,0)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("x",data_current1,i,0)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("x",data_current1,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


def xbounds2(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current2,i,0)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("x",data_current2,i,0)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("x",data_current2,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")
        
        
def xbounds3(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current3,i,0)
    
    elif 0<i<len(data_current3)-1:
        return GuessBin_Middle("x",data_current3,i,0)
    
    elif i== len(data_current3)-1:
       return GuessBin_Top("x",data_current3,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

#%%
## IMPORT THE COMPUTED WW DATA FROM THE SAVED FILE TO ADD TO THE EXPERIMENTAL TABLES
## IMPORTANT: THE POINTS FOR x AND Q² ARE THE SAME FOR BOTH 29.1 and 32.3 GeV, THIS MEANS THAT WE ONLY NEED THE WW-TERMS FROM THE PREVIOUS
## COMPUTATION, THAT'S WHY WE LOAD THEM INSTEAD OF COMPUTING THEM AGAIN
print("Reading the Wandzura-Wilczek proton data file ...")

## Table 1:
    
f1 = open(path_to_g2WW_data+"E155_1999_1P.csv")

WW_data1=[]

for line in f1:    
    WW_data1.append(line.rstrip('\n'))

f1.close()


for i in range(len(WW_data1)):    
    WW_data1[i]=WW_data1[i].split(",")    
    WW_data1[i]=[float(j) for j in WW_data1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

## Table 2:

f2 = open(path_to_g2WW_data+"E155_1999_2P.csv")

WW_data2=[]

for line in f2:    
    WW_data2.append(line.rstrip('\n'))

f2.close()


for i in range(len(WW_data2)):    
    WW_data2[i]=WW_data2[i].split(",")    
    WW_data2[i]=[float(j) for j in WW_data2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

## Table 3:
    
f3 = open(path_to_g2WW_data+"E155_1999_3P.csv")

WW_data3=[]

for line in f3:    
    WW_data3.append(line.rstrip('\n'))

f3.close()

print("Done.  =>     Convert to numbers ...")

for i in range(len(WW_data3)):    
    WW_data3[i]=WW_data3[i].split(",")    
    WW_data3[i]=[float(j) for j in WW_data3[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")    


print("Done. All 3 tables of WW data converted to numbers")

#%%
## CONSTRUCT THE FINAL TABLES 
    
print("Making G2 by E155_proton 1999 ... ")

## First we create a table which contains the data from all 3 sub-tables: bins, data_current, and WW-data:

xbin1= []
for i in range(len(data_current1)):
    xbin1.append(xbounds1(i))

xbin2= []
for i in range(len(data_current2)):
    xbin2.append(xbounds2(i))

xbin3= []
for i in range(len(data_current3)):
    xbin3.append(xbounds3(i))
    
Qbin1= []
for i in range(len(data_current1)):
    Qbin1.append(Qbounds1(i))
    
Qbin2= []
for i in range(len(data_current2)):
    Qbin2.append(Qbounds2(i))
    
Qbin3= []
for i in range(len(data_current3)):
    Qbin3.append(Qbounds3(i))

data_current=data_current1+data_current2+data_current3
WW_data=WW_data1+WW_data2+WW_data3
xbin_tot=xbin1+xbin2+xbin3
Qbin_tot=Qbin1+Qbin2+Qbin3



#%%
DataCurrent=DataProcessor.DataSet.DataSet('E155-38.p',"G2")
DataCurrent.comment="Taken from tables 1, 2, 3. Bin in Q and x guessed by us. No syst. uncert."
DataCurrent.reference="arXiv:hep-ex/9901006v1"

DataCurrent.isNormalized=False   ##Donde poner esta renormalización de la incertidumbre??
lumUncertainty=0.0
s_current=(2.*38.8*M_proton+(M_proton)**2)

for i in range(len(data_current)):
    # make up a point
    p = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p["process"]=100
    p["s"]=s_current  
    p["<Q>"]=numpy.sqrt(data_current[i][1])     
    p["Q"]=[Qbin_tot[i][0],Qbin_tot[i][1]]
    p["<x>"]=data_current[i][0]  
    p["x"]=[xbin_tot[i][0],xbin_tot[i][1]]  
    p["xSec"]=data_current[i][2]-data_current[i][0]*WW_data[i][0]
    ## THE UNCERTAINTY IS ORDERED AS: STAT(G2),SYST(G2),STAT(WW)
    p["uncorrErr"]=[(data_current[i][3]-data_current[i][4])/2,data_current[i][0]*WW_data[i][1]]
    p["thFactor"]=data_current[i][0]        
    #
    DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")

#%%
################################################################
#########----------------------------------------------#########
######------------- E155 1999 deuteron data --------------######
#########----------------------------------------------#########
################################################################

print("Reading E155 1999 deuteron data files ...")

f1 = open(path_to_data+"/E155/1999/Table1.csv")

data_from_f1=[]

for line in f1:    
    data_from_f1.append(line.rstrip('\n'))

f1.close()

print("Done.  =>     Convert table 1 to numbers ...")


data_current1=data_from_f1[62:72]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current1)):    
    data_current1[i]=data_current1[i].split(",")    
    data_current1[i]=[float(j) for j in data_current1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 1 converted to float")


f2 = open(path_to_data+"/E155/1999/Table2.csv")

data_from_f2=[]

for line in f2:    
    data_from_f2.append(line.rstrip('\n'))

f2.close()

print("Done.  =>     Convert table 2 to numbers ...")


data_current2=data_from_f2[53:60]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current2)):    
    data_current2[i]=data_current2[i].split(",")    
    data_current2[i]=[float(j) for j in data_current2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 2 converted to float")


f3 = open(path_to_data+"/E155/1999/Table3.csv")

data_from_f3=[]

for line in f3:    
    data_from_f3.append(line.rstrip('\n'))

f3.close()

print("Done.  =>     Convert table 3 to numbers ...")


data_current3=data_from_f3[47:52]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current3)):    
    data_current3[i]=data_current3[i].split(",")    
    data_current3[i]=[float(j) for j in data_current3[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 3 converted to float")


#%%
# NOW WE GUESS THE BINS IN EACH TABLE (WE HAVE TO GUESS BOTH THE BIN ON x AND Q)

## BINS ON Q:
    
def Qbounds1(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current1,i,1)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("Q",data_current1,i,1)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("Q",data_current1,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


def Qbounds2(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current2,i,1)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("Q",data_current2,i,1)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("Q",data_current2,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")
        
        
def Qbounds3(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current3,i,1)
    
    elif 0<i<len(data_current3)-1:
        return GuessBin_Middle("Q",data_current3,i,1)
    
    elif i== len(data_current3)-1:
       return GuessBin_Top("Q",data_current3,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

## BINS ON x:
    
def xbounds1(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current1,i,0)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("x",data_current1,i,0)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("x",data_current1,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


def xbounds2(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current2,i,0)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("x",data_current2,i,0)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("x",data_current2,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")
        
        
def xbounds3(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current3,i,0)
    
    elif 0<i<len(data_current3)-1:
        return GuessBin_Middle("x",data_current3,i,0)
    
    elif i== len(data_current3)-1:
       return GuessBin_Top("x",data_current3,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

#%%
## IMPORT THE COMPUTED WW DATA FROM THE SAVED FILE TO ADD TO THE EXPERIMENTAL TABLES
## IMPORTANT: THE POINTS FOR x AND Q² ARE THE SAME FOR BOTH 29.1 and 32.3 GeV, THIS MEANS THAT WE ONLY NEED THE WW-TERMS FROM THE PREVIOUS
## COMPUTATION, THAT'S WHY WE LOAD THEM INSTEAD OF COMPUTING THEM AGAIN
print("Reading the Wandzura-Wilczek proton data file ...")

## Table 1:
    
f1 = open(path_to_g2WW_data+"E155_1999_1D.csv")

WW_data1=[]

for line in f1:    
    WW_data1.append(line.rstrip('\n'))

f1.close()


for i in range(len(WW_data1)):    
    WW_data1[i]=WW_data1[i].split(",")    
    WW_data1[i]=[float(j) for j in WW_data1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

## Table 2:

f2 = open(path_to_g2WW_data+"E155_1999_2D.csv")

WW_data2=[]

for line in f2:    
    WW_data2.append(line.rstrip('\n'))

f2.close()


for i in range(len(WW_data2)):    
    WW_data2[i]=WW_data2[i].split(",")    
    WW_data2[i]=[float(j) for j in WW_data2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

## Table 3:
    
f3 = open(path_to_g2WW_data+"E155_1999_3D.csv")

WW_data3=[]

for line in f3:    
    WW_data3.append(line.rstrip('\n'))

f3.close()

print("Done.  =>     Convert to numbers ...")

for i in range(len(WW_data3)):    
    WW_data3[i]=WW_data3[i].split(",")    
    WW_data3[i]=[float(j) for j in WW_data3[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")    


print("Done. All 3 tables of WW data converted to numbers")

#%%
## CONSTRUCT THE FINAL TABLES 
    
print("Making G2 by E155_deuteron 1999 ... ")

## First we create a table which contains the data from all 3 sub-tables: bins, data_current, and WW-data:

xbin1= []
for i in range(len(data_current1)):
    xbin1.append(xbounds1(i))

xbin2= []
for i in range(len(data_current2)):
    xbin2.append(xbounds2(i))

xbin3= []
for i in range(len(data_current3)):
    xbin3.append(xbounds3(i))
    
Qbin1= []
for i in range(len(data_current1)):
    Qbin1.append(Qbounds1(i))
    
Qbin2= []
for i in range(len(data_current2)):
    Qbin2.append(Qbounds2(i))
    
Qbin3= []
for i in range(len(data_current3)):
    Qbin3.append(Qbounds3(i))

data_current=data_current1+data_current2+data_current3
WW_data=WW_data1+WW_data2+WW_data3
xbin_tot=xbin1+xbin2+xbin3
Qbin_tot=Qbin1+Qbin2+Qbin3



#%%
DataCurrent=DataProcessor.DataSet.DataSet('E155-38.d',"G2")
DataCurrent.comment="Taken from tables 1, 2, 3. Bin in Q and x guessed by us. No syst. uncert."
DataCurrent.reference="arXiv:hep-ex/9901006v1"

DataCurrent.isNormalized=False   ##Donde poner esta renormalización de la incertidumbre??
lumUncertainty=0.0
s_current=(2.*38.8*M_deuteron+(M_deuteron)**2)

for i in range(len(data_current)):
    # make up a point
    p = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p["process"]=102
    p["s"]=s_current  
    p["<Q>"]=numpy.sqrt(data_current[i][1])     
    p["Q"]=[Qbin_tot[i][0],Qbin_tot[i][1]]
    p["<x>"]=data_current[i][0]  
    p["x"]=[xbin_tot[i][0],xbin_tot[i][1]]  
    p["xSec"]=data_current[i][2]-data_current[i][0]*WW_data[i][0]
    ## THE UNCERTAINTY IS ORDERED AS: STAT(G2),SYST(G2),STAT(WW)
    p["uncorrErr"]=[(data_current[i][3]-data_current[i][4])/2,data_current[i][0]*WW_data[i][1]]
    p["thFactor"]=data_current[i][0]
    #
    DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")

#%%
################################################################
#########----------------------------------------------#########
######------------- E154 1997 neutron data ---------------######
#########----------------------------------------------#########
################################################################

print("Reading E154 neutron data files ...")

f1 = open(path_to_data+"/E154/Table1.csv")

data_from_f1=[]

for line in f1:    
    data_from_f1.append(line.rstrip('\n'))

f1.close()

print("Done.  =>     Convert table 1 to numbers ...")


data_current1=data_from_f1[28:38]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current1)):    
    data_current1[i]=data_current1[i].split(",")    
    data_current1[i]=[float(j) for j in data_current1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 1 converted to float")


f2 = open(path_to_data+"/E154/Table2.csv")

data_from_f2=[]

for line in f2:    
    data_from_f2.append(line.rstrip('\n'))

f2.close()

print("Done.  =>     Convert table 2 to numbers ...")


data_current2=data_from_f2[25:32]  ##Esto nos guarda los datos que queremos con el nombre 'data_current'

for i in range(len(data_current2)):    
    data_current2[i]=data_current2[i].split(",")    
    data_current2[i]=[float(j) for j in data_current2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
# data_list = [list(map(float, fila.split(','))) for fila in data_current]

print("table 2 converted to float")


#%%
# NOW WE GUESS THE BINS IN EACH TABLE (WE HAVE TO GUESS BOTH THE BIN ON x AND Q)

## BINS ON Q:
    
def Qbounds1(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current1,i,3)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("Q",data_current1,i,3)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("Q",data_current1,i,3)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


def Qbounds2(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current2,i,3)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("Q",data_current2,i,3)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("Q",data_current2,i,3)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")
        

#%%
## IMPORT THE COMPUTED WW DATA FROM THE SAVED FILE TO ADD TO THE EXPERIMENTAL TABLES

print("Reading the Wandzura-Wilczek neutron data file ...")

## Table 1:
    
f1 = open(path_to_g2WW_data+"E154_1N.csv")

WW_data1=[]

for line in f1:    
    WW_data1.append(line.rstrip('\n'))

f1.close()


for i in range(len(WW_data1)):    
    WW_data1[i]=WW_data1[i].split(",")    
    WW_data1[i]=[float(j) for j in WW_data1[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

## Table 2:

f2 = open(path_to_g2WW_data+"E154_2N.csv")

WW_data2=[]

for line in f2:    
    WW_data2.append(line.rstrip('\n'))

f2.close()


for i in range(len(WW_data2)):    
    WW_data2[i]=WW_data2[i].split(",")    
    WW_data2[i]=[float(j) for j in WW_data2[i]]  ##Esto nos quita las comas y nos convierte las entradas a float
    #k.spit("\t")

print("Done. Both tables of WW data converted to numbers")

#%%
## CONSTRUCT THE FINAL TABLES 
    
print("Making G2 by E154_neutron ... ")

## First we create a table which contains the data from all 3 sub-tables: bins, data_current, and WW-data:


Qbin1= []
for i in range(len(data_current1)):
    Qbin1.append(Qbounds1(i))
    
Qbin2= []
for i in range(len(data_current2)):
    Qbin2.append(Qbounds2(i))
    
data_current=data_current1+data_current2
WW_data=WW_data1+WW_data2
Qbin_tot=Qbin1+Qbin2



#%%
DataCurrent=DataProcessor.DataSet.DataSet('E154.n',"G2")
DataCurrent.comment="Taken from tables 1, 2, 3. Bin in Q and x guessed by us. No syst. uncert."
DataCurrent.reference="arXiv:hep-ex/9901006v1"

DataCurrent.isNormalized=False   ##Donde poner esta renormalización de la incertidumbre??
lumUncertainty=0.0
s_current=(2.*48.3*M_neutron+(M_neutron)**2)

for i in range(len(data_current)):
    # make up a point
    p = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p["process"]=101
    p["s"]=s_current  
    p["<Q>"]=numpy.sqrt(data_current[i][3])     
    p["Q"]=[Qbin_tot[i][0],Qbin_tot[i][1]]
    p["<x>"]=data_current[i][0]  
    p["x"]=[data_current[i][1],data_current[i][2]] 
    p["xSec"]=data_current[i][4]-WW_data[i][0]
    ## THE UNCERTAINTY IS ORDERED AS: STAT(G2),SYST(G2),STAT(WW)
    p["uncorrErr"]=[(data_current[i][5]-data_current[i][6])/2,(data_current[i][7]-data_current[i][8])/2,WW_data[i][1]]
    p["thFactor"]=1.        
    #
    DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")

#%%
#########################################################################################
#########-----------------------------------------------------------------------#########
######-------------- HERMES, proton (BOTH AVERAGED AND NON AVERAGED) --------------######
#########-----------------------------------------------------------------------#########
#########################################################################################

print("Reading NON AVERAGED HERMES data files...")

f1 = open(path_to_data+"/HERMES/Table1.csv")
data_from_f1=[]
for line in f1:    
    data_from_f1.append(line.rstrip('\n'))
f1.close()
print("Done.  =>     Convert table 1 to numbers ...")
data_current1=data_from_f1[12:35]  
for i in range(len(data_current1)):    
    data_current1[i]=data_current1[i].split(",")    
    data_current1[i]=[float(j) for j in data_current1[i]]  
print("table 1 converted to float")


print("Reading AVERAGED data files...")
f2 = open(path_to_data+"/HERMES/Table2.csv")
data_from_f2=[]
for line in f2:    
    data_from_f2.append(line.rstrip('\n'))
f2.close()
print("Done.  =>     Convert table 2 to numbers ...")
data_current2=data_from_f2[12:21]  
for i in range(len(data_current2)):    
    data_current2[i]=data_current2[i].split(",")    
    data_current2[i]=[float(j) for j in data_current2[i]]  
print("table 2 converted to float")

#%%
# NOW WE GUESS THE BINS FOR EACH TABLE ON Q
    
def Qbounds1(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current1,i,5)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("Q",data_current1,i,5)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("Q",data_current1,i,5)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


def Qbounds2(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current2,i,5)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("Q",data_current2,i,5)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("Q",data_current2,i,5)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")
        
def xbounds1(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current1,i,0)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("x",data_current1,i,0)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("x",data_current1,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


def xbounds2(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current2,i,0)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("x",data_current2,i,0)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("x",data_current2,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

#%%
## IMPORT THE COMPUTED WW DATA FROM THE SAVED FILE TO ADD TO THE EXPERIMENTAL TABLES
print("Reading the Wandzura-Wilczek proton data file ...")

## Table for non averaged:
f1 = open(path_to_g2WW_data+"HERMES.csv")
WW_data1=[]
for line in f1:    
    WW_data1.append(line.rstrip('\n'))
f1.close()

for i in range(len(WW_data1)):    
    WW_data1[i]=WW_data1[i].split(",")    
    WW_data1[i]=[float(j) for j in WW_data1[i]]  



## Table for averaged:
f2 = open(path_to_g2WW_data+"HERMES_AVG.csv")
WW_data2=[]
for line in f2:    
    WW_data2.append(line.rstrip('\n'))
f2.close()

for i in range(len(WW_data2)):    
    WW_data2[i]=WW_data2[i].split(",")    
    WW_data2[i]=[float(j) for j in WW_data2[i]]  

print("Done. Both tables of WW data converted to numbers")

#%%
## CONSTRUCT THE FINAL TABLES (WE BUILD 3 SEPARATE ONES AND MERGE THEM AT THE END TO UNITE THE WHOLE E155_PROTON EXPERIMENT)
    
print("Making G2 by HERMES ... ")

## First we create a table which contains the data from all sub-tables: bins, data_current, and WW-data:
    
Qbin1= []
for i in range(len(data_current1)):
    Qbin1.append(Qbounds1(i))
    
Qbin2= []
for i in range(len(data_current2)):
    Qbin2.append(Qbounds2(i))

#%%
# FIRST WE SAVE THE NON AVERAGED HERMES DATA:
DataCurrent=DataProcessor.DataSet.DataSet('HERMES',"G2")
DataCurrent.comment="Taken from table 1. Bin in Q guessed by us. Errors in g2 are correlated between bins over x. The correlation matrix is saved in: '/home/guillermo/Work/Twist_3/Exp_Data/HERMES/Correlation_Matrix_23BINS.csv' "
DataCurrent.reference="arXiv:1112.5584v3"

DataCurrent.isNormalized=False   
lumUncertainty=0.0
s_current=(2.*27.6*M_proton+(M_proton)**2)

for i in range(len(data_current1)):
    
    p = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p["process"]=100
    p["s"]=s_current  
    p["<Q>"]=numpy.sqrt(data_current1[i][5])     
    p["Q"]=[numpy.sqrt(0.18),numpy.sqrt(20.)]
    #print(p["<Q>"],p["Q"])
    p["<x>"]=data_current1[i][4]  
    p["x"]=[data_current1[i][2],data_current1[i][3]]
    p["xSec"]=data_current1[i][6]-data_current1[i][4]*WW_data1[i][0]
    ## THE UNCERTAINTY IS ORDERED AS: STAT(G2),SYST(G2),STAT(WW)
    #p["corrErr"]=[(data_current1[i][5]-data_current1[i][6])/(2*data_current1[i][4]),(data_current1[i][7]-data_current1[i][8])/(2*data_current1[i][4])]
    p["uncorrErr"]=[
        (data_current1[i][7]-data_current1[i][8])/2,(data_current1[i][9]-data_current1[i][10])/2,data_current1[i][4]*WW_data1[i][1]]
    p["thFactor"]=data_current1[i][4]
    #
    DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")

#%%

# AND FINALLY THE AVERAGED HERMES DATA:    
DataCurrent=DataProcessor.DataSet.DataSet('HERMES.av',"G2")
DataCurrent.comment="Taken from table 2. Bin in Q guessed by us. Errors in g2 are correlated between bins over x. The correlation matrix is saved in: '/home/guillermo/Work/Twist_3/Exp_Data/HERMES/Correlation_Matrix_9BINS.csv' "
DataCurrent.reference="arXiv:1112.5584v3"

DataCurrent.isNormalized=False   
lumUncertainty=0.0
s_current=(2.*27.6*M_proton+(M_proton)**2)

for i in range(len(data_current2)):
    
    p = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p["process"]=100
    p["s"]=s_current  
    p["<Q>"]=numpy.sqrt(data_current2[i][5])     
    p["Q"]=[numpy.sqrt(0.18),numpy.sqrt(20.)]
    p["<x>"]=data_current2[i][4]  
    p["x"]=[data_current2[i][2],data_current2[i][3]]
    p["xSec"]=data_current2[i][6]-data_current2[i][4]*WW_data2[i][0]
    #p["corrErr"]=[(data_current2[i][5]-data_current2[i][6])/(2*data_current2[i][4]),(data_current2[i][7]-data_current2[i][8])/(2*data_current2[i][4])]
    p["uncorrErr"]=[
        (data_current2[i][7]-data_current2[i][8])/2,(data_current2[i][9]-data_current2[i][10])/2,data_current2[i][4]*WW_data2[i][1]]
    p["thFactor"]=data_current2[i][4]   
    
    DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")

#%%


##########################################################
#########----------------------------------------#########
######-------------- Hall A 2004 data --------------######  
#########----------------------------------------#########
##########################################################

print("Reading Hall A 2004 data files ...")
# First we open the data for the He3:
f1 = open(path_to_data+"/JLab_Hall_A/2004/Table3.csv")
data_from_f1=[]
for line in f1:    
    data_from_f1.append(line.rstrip('\n'))
f1.close()
print("Done. => Converting tables to numbers...")

data_current1=data_from_f1[18:21]  
for i in range(len(data_current1)):    
    data_current1[i]=data_current1[i].split(",")    
    data_current1[i]=[float(j) for j in data_current1[i]]  

## Now the data for the neutron:
f2 = open(path_to_data+"/JLab_Hall_A/2004/Table6.csv")
data_from_f2=[]
for line in f2:    
    data_from_f2.append(line.rstrip('\n'))
f2.close()
print("Done. => Converting tables to numbers...")
#   
data_current2=data_from_f2[21:24]  
for i in range(len(data_current2)):    
    data_current2[i]=data_current2[i].split(",")    
    data_current2[i]=[float(j) for j in data_current2[i]]  
print("table converted to float")
#%%
# NOW WE GUESS THE BINS ON Q AND X.

## The one for the He3:
def Qbounds1(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current1,i,1)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("Q",data_current1,i,1)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("Q",data_current1,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")
        
def xbounds1(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current1,i,0)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("x",data_current1,i,0)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("x",data_current1,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

## The one for the neutron:
def Qbounds2(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current2,i,1)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("Q",data_current2,i,1)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("Q",data_current2,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")
    
def xbounds2(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current2,i,0)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("x",data_current2,i,0)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("x",data_current2,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

#%%
## IMPORT THE COMPUTED WW DATA FROM THE SAVED FILE TO ADD TO THE EXPERIMENTAL TABLES
print("Reading the Wandzura-Wilczek data files ...")

## Table 1:
# f1 = open(path_to_g2WW_data+"Hall_A_2004_He3.csv")
# WW_data1=[]
# for line in f1:    
#     WW_data1.append(line.rstrip('\n'))
# f1.close()
# for i in range(len(WW_data1)):    
#     WW_data1[i]=WW_data1[i].split(",")    
#     WW_data1[i]=[float(j) for j in WW_data1[i]]  

## Table 2:
f2 = open(path_to_g2WW_data+"Hall_A_2004_neutron.csv")
WW_data2=[]
for line in f2:    
    WW_data2.append(line.rstrip('\n'))
f2.close()
for i in range(len(WW_data2)):    
    WW_data2[i]=WW_data2[i].split(",")    
    WW_data2[i]=[float(j) for j in WW_data2[i]] 
    
print("Done. The WW data has been converted to numbers")

#%%
## CONSTRUCT THE FINAL TABLES 
    
print("Making G2 by Hall A ... ")

## First we create a table which contains the data from all sub-tables: bins, data_current, and WW-data:
    
Qbin1= []
for i in range(len(data_current1)):
    Qbin1.append(Qbounds1(i))
    
Qbin2= []
for i in range(len(data_current2)):
    Qbin2.append(Qbounds2(i))

xbin1= []
for i in range(len(data_current1)):
    xbin1.append(xbounds1(i))

xbin2= []
for i in range(len(data_current2)):
    xbin2.append(xbounds2(i))
    
#%%
# DataCurrent=DataProcessor.DataSet.DataSet('Hall_A_04_He3',"G2")
# DataCurrent.comment="Taken from table 3. Bin in Q guessed by us."
# DataCurrent.reference="arXiv:nucl-ex/0405006v5"

# DataCurrent.isNormalized=False   ##Donde poner esta renormalización de la incertidumbre??
# lumUncertainty=0.0
# s_current=(2.*5.7*M_Helium3+(M_Helium3)**2)

# for i in range(len(data_current1)):
#      # make up a point
#     p = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
#     print (DataCurrent.name+'.'+str(i))
#     p["process"]=103
#     p["s"]=s_current  
#     p["<Q>"]=numpy.sqrt(data_current1[i][1])     
#     p["Q"]=[Qbin1[i][0],Qbin1[i][1]]
#     p["<x>"]=data_current1[i][0]  
#     p["x"]=[xbin1[i][0],xbin1[i][1]]
#     p["xSec"]=data_current1[i][2]-WW_data1[i][0]
#     ## THE UNCERTAINTY IS ORDERED AS: STAT(G2),SYST(G2),STAT(WW)
#     p["uncorrErr"]=[(data_current1[i][3]-data_current1[i][4])/2,(data_current1[i][5]-data_current1[i][6])/2,WW_data1[i][1]]
#     p["thFactor"]=1.        
#     #
#     DataCurrent.AddPoint(p)

# print("Done.")

# DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
# print("This data is saved in:", path_to_save+DataCurrent.name+".csv")

DataCurrent=DataProcessor.DataSet.DataSet('HallA-2004.n',"G2")
DataCurrent.comment="Taken from table 6. Bin in Q guessed by us."
DataCurrent.reference="arXiv:nucl-ex/0405006v5"

DataCurrent.isNormalized=False   ##Donde poner esta renormalización de la incertidumbre??
lumUncertainty=0.0
s_current=(2.*5.7*M_neutron+(M_neutron)**2)

for i in range(len(data_current2)):
     # make up a point
    p = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p["process"]=101
    p["s"]=s_current  
    p["<Q>"]=numpy.sqrt(data_current2[i][1])     
    p["Q"]=[Qbin2[i][0],Qbin2[i][1]]
    p["<x>"]=data_current2[i][0]  
    p["x"]=[xbin2[i][0],xbin2[i][1]]
    p["xSec"]=data_current2[i][2]-WW_data2[i][0]
    ## THE UNCERTAINTY IS ORDERED AS: STAT(G2),SYST(G2),STAT(WW)
    p["uncorrErr"]=[(data_current2[i][3]-data_current2[i][4])/2,(data_current2[i][5]-data_current2[i][6])/2,WW_data2[i][1]]
    p["thFactor"]=1.        
    #
    DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")


#%%

##############################################################
#########--------------------------------------------#########
######-------------- Hall A 2016 He3 data --------------######  
#########--------------------------------------------#########
##############################################################

print("Reading Hall A 2016 data files ...")

f1 = open(path_to_data+"/JLab_Hall_A/2016/Tables_JLab_2016.csv")
data_from_f1=[]
for line in f1:    
    data_from_f1.append(line.rstrip('\n'))
f1.close()
print("Done. => Converting tables to numbers...")
## First the data with E=4.74 GeV:
data_current1=data_from_f1[2:15]  
for i in range(len(data_current1)):    
    data_current1[i]=data_current1[i].split(",")    
    data_current1[i]=[float(j) for j in data_current1[i]]  

## Now the data with E=5.89 GeV:
data_current2=data_from_f1[16:29]  
for i in range(len(data_current2)):    
    data_current2[i]=data_current2[i].split(",")    
    data_current2[i]=[float(j) for j in data_current2[i]]  
print("table converted to float")
#%%
# NOW WE GUESS THE BINS ON Q AND X.

## The one for the E=4.74 GeV:
def Qbounds1(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current1,i,1)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("Q",data_current1,i,1)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("Q",data_current1,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

## The one for the E=5.98 GeV:
def Qbounds2(i):
    if i==0: 
        return GuessBin_Bottom("Q",data_current2,i,1)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("Q",data_current2,i,1)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("Q",data_current2,i,1)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")
        
def xbounds1(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current1,i,0)
    
    elif 0<i<len(data_current1)-1:
        return GuessBin_Middle("x",data_current1,i,0)
    
    elif i== len(data_current1)-1:
       return GuessBin_Top("x",data_current1,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")


def xbounds2(i):
    if i==0: 
        return GuessBin_Bottom("x",data_current2,i,0)
    
    elif 0<i<len(data_current2)-1:
        return GuessBin_Middle("x",data_current2,i,0)
    
    elif i== len(data_current2)-1:
       return GuessBin_Top("x",data_current2,i,0)
    else :
        raise ValueError("Error: the value for the index i is outside the range of the data")

#%%
## IMPORT THE COMPUTED WW DATA FROM THE SAVED FILE TO ADD TO THE EXPERIMENTAL TABLES
print("Reading the Wandzura-Wilczek data files ...")

## Table 1:
f1 = open(path_to_g2WW_data+"Hall_A(4.74 GeV).csv")
WW_data1=[]
for line in f1:    
    WW_data1.append(line.rstrip('\n'))
f1.close()
for i in range(len(WW_data1)):    
    WW_data1[i]=WW_data1[i].split(",")    
    WW_data1[i]=[float(j) for j in WW_data1[i]]  

## Table 2:
f2 = open(path_to_g2WW_data+"Hall_A(5.89 GeV).csv")
WW_data2=[]
for line in f2:    
    WW_data2.append(line.rstrip('\n'))
f2.close()
for i in range(len(WW_data2)):    
    WW_data2[i]=WW_data2[i].split(",")    
    WW_data2[i]=[float(j) for j in WW_data2[i]] 
    
print("Done. The WW data has been converted to numbers")

#%%
## Finally we create a table which contains the data from all sub-tables: bins, data_current, and WW-data:
    
Qbin1= []
for i in range(len(data_current1)):
    Qbin1.append(Qbounds1(i))
    
Qbin2= []
for i in range(len(data_current2)):
    Qbin2.append(Qbounds2(i))

xbin1= []
for i in range(len(data_current1)):
    xbin1.append(xbounds1(i))

xbin2= []
for i in range(len(data_current2)):
    xbin2.append(xbounds2(i))
    
#%%
## TO FINISH IT OFF... CONSTRUCT THE FINAL TABLES.

## For the He3:
print("Making G2 by Hall A ... ")
DataCurrent=DataProcessor.DataSet.DataSet('HallA-2016-4.He3',"G2")
DataCurrent.comment="E=4.74 GeV. Taken from table VIII in paper. Bin in Q guessed by us."
DataCurrent.reference="arXiv:1603.03612v3"

DataCurrent.isNormalized=False   ##Donde poner esta renormalización de la incertidumbre??
lumUncertainty=0.0
s_current=(2.*4.74*M_neutron+(M_neutron)**2)

for i in range(len(data_current1)):
     # make up a point
    p = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p["process"]=102
    p["s"]=s_current  
    p["<Q>"]=numpy.sqrt(data_current1[i][1])     
    p["Q"]=[Qbin1[i][0],Qbin1[i][1]]
    p["<x>"]=data_current1[i][0]  
    p["x"]=[xbin1[i][0],xbin1[i][1]]
    p["xSec"]=data_current1[i][2]-WW_data1[i][0]
    ## THE UNCERTAINTY IS ORDERED AS: STAT(G2),SYST(G2),STAT(WW)
    p["uncorrErr"]=[data_current1[i][3],data_current1[i][4],WW_data1[i][1]]
    p["thFactor"]=1.        
    #
    DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")


## For the neutron:
DataCurrent=DataProcessor.DataSet.DataSet('HallA-2016-5.He3',"G2")
DataCurrent.comment="E=5.89 GeV. Taken from table IX in paper. Bin in Q guessed by us."
DataCurrent.reference="arXiv:1603.03612v3"

DataCurrent.isNormalized=False   ##Donde poner esta renormalización de la incertidumbre??
lumUncertainty=0.0
s_current=(2.*5.89*M_Helium3+(M_Helium3)**2)

for i in range(len(data_current2)):
     # make up a point
    p = DataProcessor.Point.CreateG2Point(DataCurrent.name + '.' + str(i))
    #print (DataCurrent.name+'.'+str(i))
    p["process"]=102
    p["s"]=s_current  
    p["<Q>"]=numpy.sqrt(data_current2[i][1])     
    p["Q"]=[Qbin2[i][0],Qbin2[i][1]]
    p["<x>"]=data_current2[i][0]  
    p["x"]=[xbin2[i][0],xbin2[i][1]]
    p["xSec"]=data_current2[i][2]-WW_data2[i][0]
    ## THE UNCERTAINTY IS ORDERED AS: STAT(G2),SYST(G2),STAT(WW)
    p["uncorrErr"]=[data_current2[i][3],data_current2[i][4],WW_data2[i][1]]
    p["thFactor"]=1.        
    #
    DataCurrent.AddPoint(p)

print("Done.")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
print("This data is saved in:", path_to_save+DataCurrent.name+".csv")






