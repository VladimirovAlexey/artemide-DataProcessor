#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 09:37:30 2021

Redoing of TMD-optimal grid files, and saving into weighted mean, 68% CI

Replica-0) Weighted mean
Replica-1,2) 68%CI down, up for 1+1
Replica-3) Mean for PDF case
Replica-4,5) 68%CI down, up for PDF case
Replica-6) Mean for EXP case
Replica-7,8) 68%CI down, up for EXP case

@author: vla18041
"""

caseName='SV21-CT18_nnlo'
#caseName='SV21-HERA20_nnlo'
#caseName='SV21-NNPDF31_nnlo'
#caseName='SV21-MSHT20_nnlo'
PathToGrid='/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD-grids/Grids_optimal/'
PathToSave='/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD-grids/Grids_optimal/'+caseName+'_VAR/'

from os import listdir
import numpy

### to str
def ToStrgG(ll):
    return [['{:g}'.format(x) for x in y] for y in ll]
#%%

### reading PDF case

#list file names from directory
listFiles=listdir(PathToGrid+caseName+'_PDF')

#remove .info and 0000 repica
listFiles = [el for el in listFiles if el[-3:]=='dat']
listFiles = [el for el in listFiles if el[-8:-4]!='0000']

s15m_pdf=[]
s14m_pdf=[]
s13m_pdf=[]
s12m_pdf=[]
s11m_pdf=[]
s15_pdf=[]
s14_pdf=[]
s13_pdf=[]
s12_pdf=[]
s11_pdf=[]

for f in listFiles:
    with open(PathToGrid+caseName+'_PDF/'+f,'r') as file:
        flist=file.readlines()
        
    Xline=flist[0]
    Bline=flist[1]
    
    ### the 3'd line has 'TMDs: ' in the begining and '/n' in the end. Remove them.
    ### Translate to dictionary
    values=eval(flist[2][6:-1])
    
    
    
    ### check that keys are as expected
    if(list(values.keys()) == [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]):
        s15m_pdf.append(values[-5])
        s14m_pdf.append(values[-4])
        s13m_pdf.append(values[-3])
        s12m_pdf.append(values[-2])
        s11m_pdf.append(values[-1])
        s15_pdf.append(values[5])
        s14_pdf.append(values[4])
        s13_pdf.append(values[3])
        s12_pdf.append(values[2])
        s11_pdf.append(values[1])
    else:
        print('File '+f+' seems brocken. Keys are: ', list(values.keys()))
        
#%%

### reading EXP case

#list file names from directory
listFiles=listdir(PathToGrid+caseName+'_EXP')

#remove .info and 0000 repica
listFiles = [el for el in listFiles if el[-3:]=='dat']
listFiles = [el for el in listFiles if el[-8:-4]!='0000']

s15m_exp=[]
s14m_exp=[]
s13m_exp=[]
s12m_exp=[]
s11m_exp=[]
s15_exp=[]
s14_exp=[]
s13_exp=[]
s12_exp=[]
s11_exp=[]

for f in listFiles:
    with open(PathToGrid+caseName+'_EXP/'+f,'r') as file:
        flist=file.readlines()
        
    Xline=flist[0]
    Bline=flist[1]
    
    ### the 3'd line has 'TMDs: ' in the begining and '/n' in the end. Remove them.
    ### Translate to dictionary
    values=eval(flist[2][6:-1])
    
    
    
    ### check that keys are as expected
    if(list(values.keys()) == [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]):
        s15m_exp.append(values[-5])
        s14m_exp.append(values[-4])
        s13m_exp.append(values[-3])
        s12m_exp.append(values[-2])
        s11m_exp.append(values[-1])
        s15_exp.append(values[5])
        s14_exp.append(values[4])
        s13_exp.append(values[3])
        s12_exp.append(values[2])
        s11_exp.append(values[1])
    else:
        print('File '+f+' seems brocken. Keys are: ', list(values.keys()))

#%%
### computes 68CI up boundary for variable which is smapled by sample_list
def bootsrap68CIdown(sample_list,variable):
    rr=numpy.zeros(numpy.array(variable[0]).shape)
    for i in range(len(sample_list)):
        if(numpy.mod(i,100)==0): print("(up) i=",i)
        rr=numpy.add(rr,numpy.percentile([variable[j] for j in sample_list[i]],(100-68.)/2,axis=0))
    
    return rr/len(sample_list)

### computes 68CI up boundary for variable which is smapled by sample_list
def bootsrap68CIup(sample_list,variable):
    rr=numpy.zeros(numpy.array(variable[0]).shape)
    for i in range(len(sample_list)):
        if(numpy.mod(i,100)==0): print("(down) i=",i)
        rr=numpy.add(rr,numpy.percentile([variable[j] for j in sample_list[i]],100-(100-68.)/2,axis=0))
    
    return rr/len(sample_list)

#%%
#############################################################################################
###############################PDF-case######################################################
#############################################################################################
### generate sample (since it is correalted I just generate indices)
samples=[numpy.random.choice(range(len(s15m_pdf)),size=500) for i in range(2000)]

#%%

## compute ups (it takes a lot of time!)

s15m_pdf_68up=bootsrap68CIup(samples,s15m_pdf)
s14m_pdf_68up=bootsrap68CIup(samples,s14m_pdf)
s13m_pdf_68up=bootsrap68CIup(samples,s13m_pdf)
s12m_pdf_68up=bootsrap68CIup(samples,s12m_pdf)
s11m_pdf_68up=bootsrap68CIup(samples,s11m_pdf)
s15_pdf_68up=bootsrap68CIup(samples,s15_pdf)
s14_pdf_68up=bootsrap68CIup(samples,s14_pdf)
s13_pdf_68up=bootsrap68CIup(samples,s13_pdf)
s12_pdf_68up=bootsrap68CIup(samples,s12_pdf)
s11_pdf_68up=bootsrap68CIup(samples,s11_pdf)

## compute downs (it takes a lot of time!)

s15m_pdf_68down=bootsrap68CIdown(samples,s15m_pdf)
s14m_pdf_68down=bootsrap68CIdown(samples,s14m_pdf)
s13m_pdf_68down=bootsrap68CIdown(samples,s13m_pdf)
s12m_pdf_68down=bootsrap68CIdown(samples,s12m_pdf)
s11m_pdf_68down=bootsrap68CIdown(samples,s11m_pdf)
s15_pdf_68down=bootsrap68CIdown(samples,s15_pdf)
s14_pdf_68down=bootsrap68CIdown(samples,s14_pdf)
s13_pdf_68down=bootsrap68CIdown(samples,s13_pdf)
s12_pdf_68down=bootsrap68CIdown(samples,s12_pdf)
s11_pdf_68down=bootsrap68CIdown(samples,s11_pdf)

s15m_pdf_mean=numpy.mean(s15m_pdf,axis=0)
s14m_pdf_mean=numpy.mean(s14m_pdf,axis=0)
s13m_pdf_mean=numpy.mean(s13m_pdf,axis=0)
s12m_pdf_mean=numpy.mean(s12m_pdf,axis=0)
s11m_pdf_mean=numpy.mean(s11m_pdf,axis=0)
s15_pdf_mean=numpy.mean(s15_pdf,axis=0)
s14_pdf_mean=numpy.mean(s14_pdf,axis=0)
s13_pdf_mean=numpy.mean(s13_pdf,axis=0)
s12_pdf_mean=numpy.mean(s12_pdf,axis=0)
s11_pdf_mean=numpy.mean(s11_pdf,axis=0)

#%%

### save mean
valuesList = {
-5: ToStrgG(s15m_pdf_mean),
-4: ToStrgG(s14m_pdf_mean),
-3: ToStrgG(s13m_pdf_mean),
-2: ToStrgG(s12m_pdf_mean),
-1: ToStrgG(s11m_pdf_mean),
1: ToStrgG(s11_pdf_mean),
2: ToStrgG(s12_pdf_mean),
3: ToStrgG(s13_pdf_mean),
4: ToStrgG(s14_pdf_mean),
5: ToStrgG(s15_pdf_mean)
}


with open(PathToSave+caseName+'_VAR_0003.dat', 'w') as outfile:        
    outfile.write(Xline)
    outfile.write(Bline)
    outfile.write("TMDs: "+str(valuesList).replace("'","")+"\n")

print('Grid written at  : ',PathToSave+caseName+'_VAR_0003.dat')

### save 68% down
valuesList = {
-5: ToStrgG(s15m_pdf_68down),
-4: ToStrgG(s14m_pdf_68down),
-3: ToStrgG(s13m_pdf_68down),
-2: ToStrgG(s12m_pdf_68down),
-1: ToStrgG(s11m_pdf_68down),
1: ToStrgG(s11_pdf_68down),
2: ToStrgG(s12_pdf_68down),
3: ToStrgG(s13_pdf_68down),
4: ToStrgG(s14_pdf_68down),
5: ToStrgG(s15_pdf_68down)
}


with open(PathToSave+caseName+'_VAR_0004.dat', 'w') as outfile:        
    outfile.write(Xline)
    outfile.write(Bline)
    outfile.write("TMDs: "+str(valuesList).replace("'","")+"\n")

print('Grid written at  : ',PathToSave+caseName+'_VAR_0004.dat')
    
### save 68% up
valuesList = {
-5: ToStrgG(s15m_pdf_68up),
-4: ToStrgG(s14m_pdf_68up),
-3: ToStrgG(s13m_pdf_68up),
-2: ToStrgG(s12m_pdf_68up),
-1: ToStrgG(s11m_pdf_68up),
1: ToStrgG(s11_pdf_68up),
2: ToStrgG(s12_pdf_68up),
3: ToStrgG(s13_pdf_68up),
4: ToStrgG(s14_pdf_68up),
5: ToStrgG(s15_pdf_68up)
}


with open(PathToSave+caseName+'_VAR_0005.dat', 'w') as outfile:        
    outfile.write(Xline)
    outfile.write(Bline)
    outfile.write("TMDs: "+str(valuesList).replace("'","")+"\n")

print('Grid written at  : ',PathToSave+caseName+'_VAR_0005.dat')    

#%%
#############################################################################################
###############################EXP-case######################################################
#############################################################################################
### generate sample (since it is correalted I just generate indices)
samples=[numpy.random.choice(range(len(s15m_exp)),size=75) for i in range(2000)]
#%%
## compute ups (it takes a lot of time!)

s15m_exp_68up=bootsrap68CIup(samples,s15m_exp)
s14m_exp_68up=bootsrap68CIup(samples,s14m_exp)
s13m_exp_68up=bootsrap68CIup(samples,s13m_exp)
s12m_exp_68up=bootsrap68CIup(samples,s12m_exp)
s11m_exp_68up=bootsrap68CIup(samples,s11m_exp)
s15_exp_68up=bootsrap68CIup(samples,s15_exp)
s14_exp_68up=bootsrap68CIup(samples,s14_exp)
s13_exp_68up=bootsrap68CIup(samples,s13_exp)
s12_exp_68up=bootsrap68CIup(samples,s12_exp)
s11_exp_68up=bootsrap68CIup(samples,s11_exp)

## compute downs (it takes a lot of time!)

s15m_exp_68down=bootsrap68CIdown(samples,s15m_exp)
s14m_exp_68down=bootsrap68CIdown(samples,s14m_exp)
s13m_exp_68down=bootsrap68CIdown(samples,s13m_exp)
s12m_exp_68down=bootsrap68CIdown(samples,s12m_exp)
s11m_exp_68down=bootsrap68CIdown(samples,s11m_exp)
s15_exp_68down=bootsrap68CIdown(samples,s15_exp)
s14_exp_68down=bootsrap68CIdown(samples,s14_exp)
s13_exp_68down=bootsrap68CIdown(samples,s13_exp)
s12_exp_68down=bootsrap68CIdown(samples,s12_exp)
s11_exp_68down=bootsrap68CIdown(samples,s11_exp)

s15m_exp_mean=numpy.mean(s15m_exp,axis=0)
s14m_exp_mean=numpy.mean(s14m_exp,axis=0)
s13m_exp_mean=numpy.mean(s13m_exp,axis=0)
s12m_exp_mean=numpy.mean(s12m_exp,axis=0)
s11m_exp_mean=numpy.mean(s11m_exp,axis=0)
s15_exp_mean=numpy.mean(s15_exp,axis=0)
s14_exp_mean=numpy.mean(s14_exp,axis=0)
s13_exp_mean=numpy.mean(s13_exp,axis=0)
s12_exp_mean=numpy.mean(s12_exp,axis=0)
s11_exp_mean=numpy.mean(s11_exp,axis=0)

#%%

### save mean
valuesList = {
-5: ToStrgG(s15m_exp_mean),
-4: ToStrgG(s14m_exp_mean),
-3: ToStrgG(s13m_exp_mean),
-2: ToStrgG(s12m_exp_mean),
-1: ToStrgG(s11m_exp_mean),
1: ToStrgG(s11_exp_mean),
2: ToStrgG(s12_exp_mean),
3: ToStrgG(s13_exp_mean),
4: ToStrgG(s14_exp_mean),
5: ToStrgG(s15_exp_mean)
}


with open(PathToSave+caseName+'_VAR_0006.dat', 'w') as outfile:        
    outfile.write(Xline)
    outfile.write(Bline)
    outfile.write("TMDs: "+str(valuesList).replace("'","")+"\n")

print('Grid written at  : ',PathToSave+caseName+'_VAR_0006.dat')

### save 68% down
valuesList = {
-5: ToStrgG(s15m_exp_68down),
-4: ToStrgG(s14m_exp_68down),
-3: ToStrgG(s13m_exp_68down),
-2: ToStrgG(s12m_exp_68down),
-1: ToStrgG(s11m_exp_68down),
1: ToStrgG(s11_exp_68down),
2: ToStrgG(s12_exp_68down),
3: ToStrgG(s13_exp_68down),
4: ToStrgG(s14_exp_68down),
5: ToStrgG(s15_exp_68down)
}


with open(PathToSave+caseName+'_VAR_0007.dat', 'w') as outfile:        
    outfile.write(Xline)
    outfile.write(Bline)
    outfile.write("TMDs: "+str(valuesList).replace("'","")+"\n")

print('Grid written at  : ',PathToSave+caseName+'_VAR_0007.dat')
    
### save 68% up
valuesList = {
-5: ToStrgG(s15m_exp_68up),
-4: ToStrgG(s14m_exp_68up),
-3: ToStrgG(s13m_exp_68up),
-2: ToStrgG(s12m_exp_68up),
-1: ToStrgG(s11m_exp_68up),
1: ToStrgG(s11_exp_68up),
2: ToStrgG(s12_exp_68up),
3: ToStrgG(s13_exp_68up),
4: ToStrgG(s14_exp_68up),
5: ToStrgG(s15_exp_68up)
}


with open(PathToSave+caseName+'_VAR_0008.dat', 'w') as outfile:        
    outfile.write(Xline)
    outfile.write(Bline)
    outfile.write("TMDs: "+str(valuesList).replace("'","")+"\n")

print('Grid written at  : ',PathToSave+caseName+'_VAR_0008.dat')  

#%%
##############################################################################  
###########################JOIN###############################################
##############################################################################
def ComputeWeights(a1,a2,b1,b2):
    s1=(numpy.array(a2)-numpy.array(a1))/2+0.000000001
    s2=(numpy.array(b2)-numpy.array(b1))/2+0.000000001
    
    k1=1/s1**2
    k2=1/s2**2
    return k1/(k1+k2), k2/(k1+k2)

def wMean(x1,x2,w1,w2):
    xw1=[x1[i]*w1 for i in range(len(x1))]
    xw2=[x2[i]*w2 for i in range(len(x2))]
    return numpy.mean(xw1,axis=0)+numpy.mean(xw2,axis=0)

#%%
w5m_pdf, w5m_exp=ComputeWeights(s15m_pdf_68down,s15m_pdf_68up, s15m_exp_68down, s15m_exp_68up)
w4m_pdf, w4m_exp=ComputeWeights(s14m_pdf_68down,s14m_pdf_68up, s14m_exp_68down, s14m_exp_68up)
w3m_pdf, w3m_exp=ComputeWeights(s13m_pdf_68down,s13m_pdf_68up, s13m_exp_68down, s13m_exp_68up)
w2m_pdf, w2m_exp=ComputeWeights(s12m_pdf_68down,s12m_pdf_68up, s12m_exp_68down, s12m_exp_68up)
w1m_pdf, w1m_exp=ComputeWeights(s11m_pdf_68down,s11m_pdf_68up, s11m_exp_68down, s11m_exp_68up)
w5_pdf, w5_exp=ComputeWeights(s15_pdf_68down,s15_pdf_68up, s15_exp_68down, s15_exp_68up)
w4_pdf, w4_exp=ComputeWeights(s14_pdf_68down,s14_pdf_68up, s14_exp_68down, s14_exp_68up)
w3_pdf, w3_exp=ComputeWeights(s13_pdf_68down,s13_pdf_68up, s13_exp_68down, s13_exp_68up)
w2_pdf, w2_exp=ComputeWeights(s12_pdf_68down,s12_pdf_68up, s12_exp_68down, s12_exp_68up)
w1_pdf, w1_exp=ComputeWeights(s11_pdf_68down,s11_pdf_68up, s11_exp_68down, s11_exp_68up)

#%%
### save weighted center
valuesList = {
-5: ToStrgG(wMean(s15m_pdf,s15m_exp,w5m_pdf,w5m_exp)),
-4: ToStrgG(wMean(s14m_pdf,s14m_exp,w4m_pdf,w4m_exp)),
-3: ToStrgG(wMean(s13m_pdf,s13m_exp,w3m_pdf,w3m_exp)),
-2: ToStrgG(wMean(s12m_pdf,s12m_exp,w2m_pdf,w2m_exp)),
-1: ToStrgG(wMean(s11m_pdf,s11m_exp,w1m_pdf,w1m_exp)),
1: ToStrgG(wMean(s11_pdf,s11_exp,w1_pdf,w1_exp)),
2: ToStrgG(wMean(s12_pdf,s12_exp,w2_pdf,w2_exp)),
3: ToStrgG(wMean(s13_pdf,s13_exp,w3_pdf,w3_exp)),
4: ToStrgG(wMean(s14_pdf,s14_exp,w4_pdf,w4_exp)),
5: ToStrgG(wMean(s15_pdf,s15_exp,w5_pdf,w5_exp))
}


with open(PathToSave+caseName+'_VAR_0000.dat', 'w') as outfile:        
    outfile.write(Xline)
    outfile.write(Bline)
    outfile.write("TMDs: "+str(valuesList).replace("'","")+"\n")

print('Grid written at  : ',PathToSave+caseName+'_VAR_0000.dat')

#%%
def JoinedSample(a1,a2):    
    i1=numpy.random.choice(range(len(a1)),size=75)
    i2=numpy.random.choice(range(len(a2)),size=75)
    return [a1[i] for i in i1]+[a2[i] for i in i2]

### computes 68CI up boundary for variable which is smapled by sample_list
def bootsrap68CIdown_direct(a1,a2):
    NUM=2000
    rr=numpy.zeros(numpy.array(a1[0]).shape)
    for i in range(NUM):
        if(numpy.mod(i,10)==0): print("(up) i=",i , "/",NUM)
        rr=numpy.add(rr,numpy.percentile(JoinedSample(a1,a2),(100-68.)/2,axis=0))
    
    return rr/NUM

### computes 68CI up boundary for variable which is smapled by sample_list
def bootsrap68CIup_direct(a1,a2):
    NUM=2000
    rr=numpy.zeros(numpy.array(a1[0]).shape)
    for i in range(NUM):
        if(numpy.mod(i,10)==0): print("(down) i=",i , "/",NUM)
        rr=numpy.add(rr,numpy.percentile(JoinedSample(a1,a2),100-(100-68.)/2,axis=0))
    
    return rr/NUM

#%%
### save 68% down
valuesList = {
-5: ToStrgG(bootsrap68CIup_direct(s15m_pdf,s15m_exp)),
-4: ToStrgG(bootsrap68CIup_direct(s14m_pdf,s14m_exp)),
-3: ToStrgG(bootsrap68CIup_direct(s13m_pdf,s13m_exp)),
-2: ToStrgG(bootsrap68CIup_direct(s12m_pdf,s12m_exp)),
-1: ToStrgG(bootsrap68CIup_direct(s11m_pdf,s11m_exp)),
1: ToStrgG(bootsrap68CIup_direct(s11_pdf,s11_exp)),
2: ToStrgG(bootsrap68CIup_direct(s12_pdf,s12_exp)),
3: ToStrgG(bootsrap68CIup_direct(s13_pdf,s13_exp)),
4: ToStrgG(bootsrap68CIup_direct(s14_pdf,s14_exp)),
5: ToStrgG(bootsrap68CIup_direct(s15_pdf,s15_exp))
}


with open(PathToSave+caseName+'_VAR_0002.dat', 'w') as outfile:        
    outfile.write(Xline)
    outfile.write(Bline)
    outfile.write("TMDs: "+str(valuesList).replace("'","")+"\n")

print('Grid written at  : ',PathToSave+caseName+'_VAR_0002.dat')

### save 68% down
valuesList = {
-5: ToStrgG(bootsrap68CIdown_direct(s15m_pdf,s15m_exp)),
-4: ToStrgG(bootsrap68CIdown_direct(s14m_pdf,s14m_exp)),
-3: ToStrgG(bootsrap68CIdown_direct(s13m_pdf,s13m_exp)),
-2: ToStrgG(bootsrap68CIdown_direct(s12m_pdf,s12m_exp)),
-1: ToStrgG(bootsrap68CIdown_direct(s11m_pdf,s11m_exp)),
1: ToStrgG(bootsrap68CIdown_direct(s11_pdf,s11_exp)),
2: ToStrgG(bootsrap68CIdown_direct(s12_pdf,s12_exp)),
3: ToStrgG(bootsrap68CIdown_direct(s13_pdf,s13_exp)),
4: ToStrgG(bootsrap68CIdown_direct(s14_pdf,s14_exp)),
5: ToStrgG(bootsrap68CIdown_direct(s15_pdf,s15_exp))
}


with open(PathToSave+caseName+'_VAR_0001.dat', 'w') as outfile:        
    outfile.write(Xline)
    outfile.write(Bline)
    outfile.write("TMDs: "+str(valuesList).replace("'","")+"\n")

print('Grid written at  : ',PathToSave+caseName+'_VAR_0001.dat')