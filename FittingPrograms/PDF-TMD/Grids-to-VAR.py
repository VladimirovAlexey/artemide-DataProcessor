#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 09:37:30 2021

Rewding of TMD-optimal grid files, and saving into mean, mean+-STD

@author: vla18041
"""

#caseName='SV21-CT18_nnlo_EXP'
#caseName='SV21-HERA20_nnlo_EXP'
#caseName='SV21-NNPDF31_nnlo_EXP'
caseName='SV21-MSHT20_nnlo_EXP'
PathToGrid='/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD-grids/Grids_optimal/'+caseName+'/'
PathToSave='/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD-grids/Grids_optimal/'+caseName+'_VAR/'

#%%
from os import listdir

#list file names from directory
listFiles=listdir(PathToGrid)

#remove .info and 0000 repica
listFiles = [el for el in listFiles if el[-3:]=='dat']
listFiles = [el for el in listFiles if el[-8:-4]!='0000']

#%%
import numpy

def ReadGridAndSetInitial(fname):
    global NUM,s15m,s14m,s13m,s12m,s11m,s11,s12,s13,s14,s15,s25m,s24m,s23m,s22m,s21m,s21,s22,s23,s24,s25,Xline,Bline
    
    
    with open(PathToGrid+fname,'r') as file:
        flist=file.readlines()
    
    Xline=flist[0]
    Bline=flist[1]
    
    ### the 3'd line has 'TMDs: ' in the begining and '/n' in the end. Remove them.
    ### Translate to dictionary
    values=eval(flist[2][6:-1])
    
    ### check that keys are as expected
    if(list(values.keys()) == [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]):
        s15m=numpy.zeros(numpy.array(values[-5]).shape)
        s25m=numpy.zeros(numpy.array(values[-5]).shape)
        s14m=numpy.zeros(numpy.array(values[-4]).shape)
        s24m=numpy.zeros(numpy.array(values[-4]).shape)
        s13m=numpy.zeros(numpy.array(values[-3]).shape)
        s23m=numpy.zeros(numpy.array(values[-3]).shape)
        s12m=numpy.zeros(numpy.array(values[-2]).shape)
        s22m=numpy.zeros(numpy.array(values[-2]).shape)
        s11m=numpy.zeros(numpy.array(values[-1]).shape)
        s21m=numpy.zeros(numpy.array(values[-1]).shape)
        s15=numpy.zeros(numpy.array(values[5]).shape)
        s25=numpy.zeros(numpy.array(values[5]).shape)
        s14=numpy.zeros(numpy.array(values[4]).shape)
        s24=numpy.zeros(numpy.array(values[4]).shape)
        s13=numpy.zeros(numpy.array(values[3]).shape)
        s23=numpy.zeros(numpy.array(values[3]).shape)
        s12=numpy.zeros(numpy.array(values[2]).shape)
        s22=numpy.zeros(numpy.array(values[2]).shape)
        s11=numpy.zeros(numpy.array(values[1]).shape)
        s21=numpy.zeros(numpy.array(values[1]).shape)
    else:
        print('File '+fname+' seems brocken. Keys are: ', list(values.keys()))
        return
    
    NUM=0

#%%

def ReadGridAndAddTerm(fname):
    global NUM,s15m,s14m,s13m,s12m,s11m,s11,s12,s13,s14,s15,s25m,s24m,s23m,s22m,s21m,s21,s22,s23,s24,s25 
    
    with open(PathToGrid+fname,'r') as file:
        flist=file.readlines()
    
    ### the 3'd line has 'TMDs: ' in the begining and '/n' in the end. Remove them.
    ### Translate to dictionary
    values=eval(flist[2][6:-1])
    
    ### check that keys are as expected
    if(list(values.keys()) == [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]):
        s15m=numpy.add(s15m,numpy.array(values[-5]))
        s14m=numpy.add(s14m,numpy.array(values[-4]))
        s13m=numpy.add(s13m,numpy.array(values[-3]))
        s12m=numpy.add(s12m,numpy.array(values[-2]))
        s11m=numpy.add(s11m,numpy.array(values[-1]))
        s15=numpy.add(s15,numpy.array(values[5]))
        s14=numpy.add(s14,numpy.array(values[4]))
        s13=numpy.add(s13,numpy.array(values[3]))
        s12=numpy.add(s12,numpy.array(values[2]))
        s11=numpy.add(s11,numpy.array(values[1]))
        
        s25m=numpy.add(s25m,numpy.square(numpy.array(values[-5])))
        s24m=numpy.add(s24m,numpy.square(numpy.array(values[-4])))
        s23m=numpy.add(s23m,numpy.square(numpy.array(values[-3])))
        s22m=numpy.add(s22m,numpy.square(numpy.array(values[-2])))
        s21m=numpy.add(s21m,numpy.square(numpy.array(values[-1])))
        s25=numpy.add(s25,numpy.square(numpy.array(values[5])))
        s24=numpy.add(s24,numpy.square(numpy.array(values[4])))
        s23=numpy.add(s23,numpy.square(numpy.array(values[3])))
        s22=numpy.add(s22,numpy.square(numpy.array(values[2])))
        s21=numpy.add(s21,numpy.square(numpy.array(values[1])))
    else:
        print('File '+fname+' seems brocken. Keys are: ', list(values.keys()))
        return
    
    NUM+=1
    

#%%
################# Long part read all and adds the expressions
ReadGridAndSetInitial(listFiles[0])
for f in listFiles:
    ReadGridAndAddTerm(f)
    
#%%
### mean values
mean5m=s15m/NUM
mean4m=s14m/NUM
mean3m=s13m/NUM
mean2m=s12m/NUM
mean1m=s11m/NUM
mean5=s15/NUM
mean4=s14/NUM
mean3=s13/NUM
mean2=s12/NUM
mean1=s11/NUM

### std's
std5m=numpy.sqrt(abs(s25m/NUM - numpy.square(mean5m)))
std4m=numpy.sqrt(abs(s24m/NUM - numpy.square(mean4m)))
std3m=numpy.sqrt(abs(s23m/NUM - numpy.square(mean3m)))
std2m=numpy.sqrt(abs(s22m/NUM - numpy.square(mean2m)))
std1m=numpy.sqrt(abs(s21m/NUM - numpy.square(mean1m)))
std5=numpy.sqrt(abs(s25/NUM - numpy.square(mean5)))
std4=numpy.sqrt(abs(s24/NUM - numpy.square(mean4)))
std3=numpy.sqrt(abs(s23/NUM - numpy.square(mean3)))
std2=numpy.sqrt(abs(s22/NUM - numpy.square(mean2)))
std1=numpy.sqrt(abs(s21/NUM - numpy.square(mean1)))

#%%

### to str
def ToStrgG(ll):
    return [['{:g}'.format(x) for x in y] for y in ll]

### save mean
valuesList = {
-5: ToStrgG(mean5m),
-4: ToStrgG(mean4m),
-3: ToStrgG(mean3m),
-2: ToStrgG(mean2m),
-1: ToStrgG(mean1m),
1: ToStrgG(mean1),
2: ToStrgG(mean2),
3: ToStrgG(mean3),
4: ToStrgG(mean4),
5: ToStrgG(mean5)
}


with open(PathToSave+caseName+'_VAR_0000.dat', 'w') as outfile:        
    outfile.write(Xline)
    outfile.write(Bline)
    outfile.write("TMDs: "+str(valuesList).replace("'","")+"\n")

print('Grid written at  : ',PathToSave+caseName+'_VAR_0000.dat')

### save mean+std
valuesList = {
-5: ToStrgG(mean5m+std5m),
-4: ToStrgG(mean4m+std4m),
-3: ToStrgG(mean3m+std3m),
-2: ToStrgG(mean2m+std2m),
-1: ToStrgG(mean1m+std1m),
1: ToStrgG(mean1+std1),
2: ToStrgG(mean2+std2),
3: ToStrgG(mean3+std3),
4: ToStrgG(mean4+std4),
5: ToStrgG(mean5+std5)
}


with open(PathToSave+caseName+'_VAR_0001.dat', 'w') as outfile:        
    outfile.write(Xline)
    outfile.write(Bline)
    outfile.write("TMDs: "+str(valuesList).replace("'","")+"\n")

print('Grid written at  : ',PathToSave+caseName+'_VAR_0001.dat')

### save mean-std
valuesList = {
-5: ToStrgG(mean5m-std5m),
-4: ToStrgG(mean4m-std4m),
-3: ToStrgG(mean3m-std3m),
-2: ToStrgG(mean2m-std2m),
-1: ToStrgG(mean1m-std1m),
1: ToStrgG(mean1-std1),
2: ToStrgG(mean2-std2),
3: ToStrgG(mean3-std3),
4: ToStrgG(mean4-std4),
5: ToStrgG(mean5-std5)
}


with open(PathToSave+caseName+'_VAR_0002.dat', 'w') as outfile:        
    outfile.write(Xline)
    outfile.write(Bline)
    outfile.write("TMDs: "+str(valuesList).replace("'","")+"\n")

print('Grid written at  : ',PathToSave+caseName+'_VAR_0002.dat')
