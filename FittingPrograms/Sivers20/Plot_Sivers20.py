#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 09:04:35 2020

@author: vla18041
"""

#######################################
# importing libraries
#######################################

import sys
import time
import numpy
import scipy.stats
#sys.path.append("/home/m/Github/artemide-DataProcessor")
sys.path.append("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor")
import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet
import DataProcessor.ArtemideReplicaSet

#MAINPATH="/home/m/Github/artemide-DataProcessor/"
MAINPATH="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"

#%%
#######################################
#Initialize artemide
#######################################
import harpy
path_to_constants=MAINPATH+"FittingPrograms/Sivers20/Constants-files/"
harpy.initialize(path_to_constants+"const-Sivers20_plot_nnlo")

#### All=0 Case
harpy.setNPparameters_TMDR([2., 0.0398333])


harpy.setNPparameters_SiversTMDPDF([5.2, 0.,0.,0.,0., -0.6, 15.9, 0.5, -0.2, 21.6, -0.5, -0.1, 0.4, -1.1]) 

#%%
##########################################
# Loading replicas
#########################################
rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(
"/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/Sivers20_model9case1.rep")
#"/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/Sivers20_model9case1(noDY)_1000.rep")
#"/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/Sivers20_model9case1(noDY)NP.rep")
#"/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/Sivers20_model9case1(noDY-n3lo).rep")
#"/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/Sivers20_model9case1(noDY-n3lo)NP.rep")
#"/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/Sivers20_model9case1(n3lo).rep")
#"/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/Sivers20_model9case1(n3lo)NP.rep")
meanReplica=rSet.GetReplica(0)
rSet.SetReplica()

#%%
########################################
# Function that return the value of Sivers function distributed over replicas
# f= flavor in LHAPDF convention (-5,....,5)
# mu is scale. if present return F(x,b,mu,mu^2)
########################################
def getSiversRow(x,b,f,mu=-1):
    result=[]
    if(mu==-1):
        for j in range(rSet.numberOfReplicas):
            rSet.SetReplica(j)
            result.append(harpy.get_SiversTMDPDF(x,b,1)[f+5])
    else:
        for j in range(rSet.numberOfReplicas):
            rSet.SetReplica(j)
            result.append(harpy.get_SiversTMDPDF(x,b,1,mu=mu)[f+5])
    return result

#%%
########################################
# Function that return the value of Sivers function distributed over replicas
# f= flavor in LHAPDF convention (-5,....,5)
# mu is scale. if present return F(x,b,mu,mu^2)
########################################
def getSiversRowP(x,p,f,mu=-1):
    result=[]
    M2_proton=0.932**2
    if(mu==-1):
        for j in range(rSet.numberOfReplicas):
            rSet.SetReplica(j)
            result.append(harpy.get_SiversTMDPDF_kT(x,p,1)[f+5]*M2_proton/p)
    else:
        for j in range(rSet.numberOfReplicas):
            rSet.SetReplica(j)
            result.append(harpy.get_SiversTMDPDF_kT(x,p,1,mu=mu)[f+5]*M2_proton/p)
    return result

#%%
#######################################
## Find a mode of the sample using HalfSampleMode method 
## https://arxiv.org/pdf/math/0505419.pdf
######################################
def HalfSampleMode(X):
    
    ### step 1
    if(len(X)==1):
        return [X[0]]
    ### step 2
    if(len(X)==2):
        return [(X[0]+X[1])/2]
    ### step 3
    if(len(X)==3):
        if(X[1]-X[0]<X[2]-X[1]):
            return [(X[0]+X[1])/2]
        elif(X[1]-X[0]>X[2]-X[1]):
            return [(X[1]+X[2])/2]
        else:
            return [X[1]]
        
    ### step 4
    n=len(X)
    nnew=int(numpy.ceil(n/2))
    j=0
    wMin=X[-1]-X[0]
    for i in range(n-nnew+1):
        w=X[i+nnew-1]-X[i]
        if(w<wMin):
            wMin=w
            j=i
        
    #print(">>",j)
    return X[j:j+nnew]

def FindMode(dd):
    
    X=numpy.sort(dd)
    
    while len(X)>1:
        X=HalfSampleMode(numpy.array(X))
    return X[0]

#%%
#########################################################
## Resample data with N/2
#########################################################
def Resample(dd):
    #return numpy.random.choice(dd,size=int(numpy.floor(len(dd)/2)))
    return dd[numpy.random.choice(dd.shape[0], size=int(numpy.floor(len(dd)/2)))]

#%%
#########################################################
## Determine Mean, Mode, 68%CI by resampling 
#########################################################
alpha=68
def ComputeParameters(dd):    
    means=[]
    modes=[]
    lowers=[]
    uppers=[]    
    for i in range(1500):
        sample=Resample(numpy.array(dd))
        means.append(numpy.mean(sample))
        modes.append(FindMode(sample))
        lowers.append(numpy.percentile(sample,(100-alpha)/2))
        uppers.append(numpy.percentile(sample,100-(100-alpha)/2))
    
    return [numpy.mean(means),numpy.mean(modes),numpy.mean(lowers),numpy.mean(uppers)]
#%%
# #########################################################
# ## Compute parameters of the NP parameters first diagonalizing the correlation matrix
# #########################################################
# ### Computing the diagonalization matrix and diagonalizing set
# #indices=[0,1,2,5,6,7,8,9,10,11,12]
# covarianceM=numpy.cov(numpy.transpose(rSet.GetParameterSet()))
# eigenM=numpy.linalg.eig(covarianceM)[1]
# eigenMinv=numpy.linalg.inv(eigenM)
# diagonalSet=[numpy.matmul(eigenMinv, x) for x in rSet.GetParameterSet()]

# ### computing prameters for diagonalized set
# modeReplica=[]
# meanReplica=[]
# downReplica=[]
# upReplica=[]
# for j in range(14):
#     example=[]
#     for i in range(len(diagonalSet)):
#         example.append(diagonalSet[i][j])
        
#     params=ComputeParameters(example)
#     modeReplica.append(params[1])
#     meanReplica.append(params[0])
#     downReplica.append(params[2])
#     upReplica.append(params[3])
    
# ### turning back
# modeReplica=numpy.matmul(eigenM,modeReplica)
# meanReplica=numpy.matmul(eigenM,meanReplica)
# downReplica=numpy.matmul(eigenM,downReplica)
# upReplica=numpy.matmul(eigenM,upReplica)

#%%
#########################################################
## Compute parameters of the NP parameters with out rotation
#########################################################
modeReplica=[]
meanReplica=[]
downReplica=[]
upReplica=[]
for j in range(14):
    example=[]
    for i in range(1,rSet.numberOfReplicas+1):
        example.append(rSet.GetReplica(i)[j])
        
    params=ComputeParameters(example)
    modeReplica.append(params[1])
    meanReplica.append(params[0])
    downReplica.append(params[2])
    upReplica.append(params[3])
    print("param=",j," >>  ",params)

#%%
for i in  range(14):
    rrr=rSet.GetReplica(0)
    print("{:2.4f},{:2.4f},{:2.4f}".format(numpy.round(rrr[i],3),
                                             numpy.round(downReplica[i]-rrr[i],3),
                                             numpy.round(upReplica[i]-rrr[i],3)))
#%%
rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/Sivers20/REPS/"+
                                                  "Sivers20_BPV20(n3lo).rep")
                                                  #"Sivers20_BPV20(nnlo).rep")
meanReplica=rSet.GetReplica(0)
rSet.SetReplica()
#%%
#########################################################
## Evaluates pp which is the list of [b,mean,low,up] for various values of x
#########################################################
pp=[]
xValues=[0.001,0.005,0.01,0.05,0.1,0.5]
f=1
for xx in xValues:
    kk=[]
    for j in range(41):
        b=0.1*j
        p0=ComputeParameters(getSiversRow(xx,b,f))
        rSet.SetReplica(0)        
        #harpy.setNPparameters_SiversTMDPDF(meanReplica)
        tmd=harpy.get_SiversTMDPDF(xx, b, 1)        
        kk.append([b,tmd[f+5],p0[0],p0[1],p0[2],p0[3]])
    pp.append(kk)
    
#%%
#########################################################
## Evaluates pp which is the list of [kT,mean,low,up] for various values of x
#########################################################
pp=[]
xValues=[0.001,0.005,0.01,0.05,0.1,0.5]
f=1
for xx in xValues:
    M2_proton=0.932**2
    kk=[]
    for j in range(21):
        kT=0.05*j
        if(j==0):
            kT=0.001
        p0=ComputeParameters(getSiversRowP(xx,kT,f,mu=2.))
        rSet.SetReplica(0)
        tmd=harpy.get_SiversTMDPDF_kT(xx, kT, 1,mu=2.)      
        kk.append([kT,M2_proton/kT*tmd[f+5],p0[1],p0[2],p0[3]])
    pp.append(kk)

#%%
print("{",end="")
for j in range(len(pp)):
    kk=pp[j]
    print("{",end="")
    for i in range(len(kk)):
        #print("{","{:2.4f},{:12.9f},{:12.9f},{:12.9f},{:12.9f}".format(kk[i][0],kk[i][1],kk[i][2],kk[i][3],kk[i][4]),"}",end="")    
        print("{","{:2.4f},{:12.9f},{:12.9f},{:12.9f}".format(kk[i][0],kk[i][2],kk[i][4],kk[i][5]),"}",end="")    
        if(i==len(kk)-1):
            print("}",end="")
        else:
            print(",")
            
    if(j==len(pp)-1):
        print("}",end="")
    else:
        print(",")
        
#%%
from matplotlib import pyplot

j=4

pyplot.plot([i[0] for i in pp[j]],[i[1] for i in pp[j]],color="black")
pyplot.plot([i[0] for i in pp[j]],[i[2] for i in pp[j]],color="red")
pyplot.plot([i[0] for i in pp[j]],[i[3] for i in pp[j]],color="green")
pyplot.plot([i[0] for i in pp[j]],[i[4] for i in pp[j]],color="orange")
pyplot.plot([i[0] for i in pp[j]],[i[5] for i in pp[j]],color="orange")
pyplot.show()


#%%
from matplotlib import pyplot

#pyplot.xlim(0,10)
pyplot.hist(example,bins=50)
pyplot.axvline(params[0],color="red")
pyplot.axvline(params[1],color="green")
pyplot.axvspan(params[2],params[3],color="green",alpha=0.2)
pyplot.show()        

#%%
pp=[]
xValues=[0.001,0.005,0.01,0.05,0.1,0.5]
f=1
for xx in xValues:
    kk=[]
    for j in range(41):
        b=0.1*j
        p0=ComputeParameters(getSiversRow(xx,b,f))
        rSet.SetReplica(0)        
        #harpy.setNPparameters_SiversTMDPDF(meanReplica)
        tmd=harpy.get_SiversTMDPDF(xx, b, 1)        
        kk.append([b,tmd[f+5],p0[0],p0[1],p0[2],p0[3]])
    pp.append(kk)
    
#%%
pp=[]
for i in range(1,rSet.numberOfReplicas+1):
    rSet.SetReplica(i)   
    tmd=harpy.get_SiversTMDPDF(0.1, 0.5, 1)
    pp.append(tmd[1+5])
    