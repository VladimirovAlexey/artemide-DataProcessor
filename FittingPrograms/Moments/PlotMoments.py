#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 16:42:03 2023

@author: alexey
"""

#######################################
# importing libraries
#######################################
import os
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"

ATMDE_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide/"
HARPY_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide/harpy/"

#replicaFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/replicaDY.dat"
#logFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/LOG.log"

import sys
import numpy
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)

SAVEPATH="/data/WorkingFiles/TMD/Fit_Notes/MomentsOfTMDs/data/"

#%%
import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet

MAINPATH=ROOT_DIR 

#%%
#######################################
#Initialize artemide
#######################################
import harpy

path_to_constants=MAINPATH+"FittingPrograms/Moments/ConstantsFiles/ART23_MSHT_N4LL.atmde"

harpy.initialize(path_to_constants)

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/data/WorkingFiles/TMD/Fit_Notes/ART23/REPS/ART23_run2.rep")
rSet.SetReplica(-1)

rSet.SetReplica(0)
#%%
rSetS=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            ATMDE_DIR+"Models/BPV20/Replica-files/BPV20(n3lo).rep")
rSetS.SetReplica(-1)

rSetS.SetReplica(0)
#%%
with open(SAVEPATH+'G0_ART23_mu_x.0000', 'w') as outfile:
    for i in range(81):
        for j in range(35):
            mu=(50)**(j/40)
            x=10**(-0.05*i)
            G0=harpy.get_uTMDPDF_G0(x,mu,1)
            G0str=[str(g) for g in G0]
            outfile.write(str(mu)+",   "+str(x)+",   "+",   ".join(G0str)+"\n")
   
#%%
for i in range(81):
    mu=10.
    x=10**(-0.05*i)
    #G0=harpy.get_uTMDPDF_G0(x,mu,1)
    pdf=harpy.get_uPDF(x,mu,1)
    #print('{'+str(mu)+','+str(x*G0[7])+','+str(x*pdf[7])+'},')
    print('{'+str(x)+','+str(pdf[6])+'},')

#%%
for i in range(1,101):
    x=0.1
    #mu=(100)**(i/40)
    mu=1.*i
    #mu=10.
    #x=10**(-i/10)
    X0=harpy.get_uTMDPDF_X0(x,mu,1)
    #pdf=harpy.get_uPDF(x,mu,1)
    print('{'+str(mu)+','+str(X0[3])+','+str(X0[4])+','+str(X0[6])+','+str(X0[7])+'},')
    
#%%
for i in range(41):
    #x=0.1
    #mu=(100)**(i/40)
    mu=10.
    x=10**(-i/10)
    X0=harpy.get_uPDF(x,mu,1)
    #pdf=harpy.get_uPDF(x,mu,1)
    #print('{'+str(mu)+','+str(x*G0[7])+','+str(x*pdf[7])+'},')
    print('{'+str(x)+','+str(X0[3])+','+str(X0[4])+','+str(X0[6])+','+str(X0[7])+'},')
    
#%%
for i in range(201):
    x=0.01
    b=100*(0.00000001)**(i/200)
    X0=harpy.get_uTMDPDF(x,b,1)
    #pdf=harpy.get_uPDF(x,mu,1)
    print('{'+str(b)+','+str(X0[3])+','+str(X0[4])+','+str(X0[6])+','+str(X0[7])+'},')
    
#%%
for i in range(201):
    x=0.01
    b=100*(0.00000001)**(i/200)
    mu=1.12292/b+5
    X0=harpy.get_uPDF(x,mu,1)
    #pdf=harpy.get_uPDF(x,mu,1)
    print('{'+str(b)+','+str(X0[3])+','+str(X0[4])+','+str(X0[6])+','+str(X0[7])+'},')
    
#%%
muValues=[i*1. for i in range(1,101)]
xValues=[10**(-i/10) for i in range(51)]

result=[]
x=0.01
mu=50.


#for mu in muValues:
for x in xValues:
    X0=harpy.get_uTMDPDF_X0(x,mu,1)
    AS=harpy.get_uTMDPDF_ASX0(x,mu,1)
    print("{",x,",",X0[5+1],',',X0[5+2],',',X0[5+3],",",X0[5-1],',',X0[5-2],',',X0[5-3],"},")
    #print("{",x,",",AS[5+1],',',AS[5+2],',',AS[5+3],",",AS[5-1],',',AS[5-2],',',AS[5-3],"},")
#%%
######################
### This computed the uncertanty for X0-bar
#####################
#########################################################
## Resample data with N/2
#########################################################
def Resample(dd):
    #return numpy.random.choice(dd,size=int(numpy.floor(len(dd)/2)))
    return dd[numpy.random.choice(dd.shape[0], size=int(numpy.floor(len(dd)/2)))]

#########################################################
## Determine Mean, Mode, 68%CI by resampling 
#########################################################
alpha=68
def Compute68CI(dd):    
    lowers=[]
    uppers=[]    
    for i in range(1500):
        sample=Resample(numpy.array(dd))
        lowers.append(numpy.percentile(sample,(100-alpha)/2,axis=0))
        uppers.append(numpy.percentile(sample,100-(100-alpha)/2,axis=0))
    
    return [numpy.mean(lowers,axis=0),numpy.mean(uppers,axis=0)]

#rSet.SetReplica(0)


muValues=[i*1. for i in range(1,21)]
xValues=[10**(-i/10) for i in range(41)]

result=[]
x=0.01
mu=10.
#a0=1.014
a0=1.021
result_u=[]
result_d=[]
result_s=[]
result_ub=[]
result_db=[]
result_sb=[]

result_u_val=[]
result_d_val=[]
result_s_val=[]

for r in range(200):
    rnd=numpy.random.randint(1,high=1000)        
    #rSet.SetReplica(rnd)
    harpy.setPDFreplica(rnd)
    
    ss_u=[]
    ss_d=[]
    ss_ub=[]
    ss_db=[]
    ss_s=[]
    ss_sb=[]
    
    ss_u_val=[]
    ss_d_val=[]
    ss_s_val=[]
    #for mu in muValues:
    for x in xValues:
        X0=harpy.get_uTMDPDF_X0(x,mu,1)
        AS=harpy.get_uTMDPDF_ASX0(x,mu,1)
        ss_sb.append(X0[5-3]-0.5*(a0*mu)**2*AS[5-3])
        ss_ub.append(X0[5-2]-0.5*(a0*mu)**2*AS[5-2])
        ss_db.append(X0[5-1]-0.5*(a0*mu)**2*AS[5-1])        
        ss_d.append(X0[5+1]-0.5*(a0*mu)**2*AS[5+1])
        ss_u.append(X0[5+2]-0.5*(a0*mu)**2*AS[5+2])
        ss_s.append(X0[5+3]-0.5*(a0*mu)**2*AS[5+3])
        
        ss_d_val.append(X0[5+1]-X0[5-1]-0.5*(a0*mu)**2*(AS[5+1]-AS[5-1]))
        ss_u_val.append(X0[5+2]-X0[5-2]-0.5*(a0*mu)**2*(AS[5+2]-AS[5-2]))
        ss_s_val.append(X0[5+3]-X0[5-3]-0.5*(a0*mu)**2*(AS[5+3]-AS[5-3]))
    
    result_u.append(ss_u)
    result_d.append(ss_d)
    result_s.append(ss_s)
    result_ub.append(ss_ub)
    result_db.append(ss_db)
    result_sb.append(ss_sb)
    
    result_u_val.append(ss_u_val)
    result_d_val.append(ss_d_val)
    result_s_val.append(ss_s_val)
    print(r)
#%%

X0_U=[numpy.mean(result_u,axis=0)]+Compute68CI(result_u)
X0_D=[numpy.mean(result_d,axis=0)]+Compute68CI(result_d)
X0_S=[numpy.mean(result_s,axis=0)]+Compute68CI(result_s)
X0_Ub=[numpy.mean(result_ub,axis=0)]+Compute68CI(result_ub)
X0_Db=[numpy.mean(result_db,axis=0)]+Compute68CI(result_db)
X0_Sb=[numpy.mean(result_sb,axis=0)]+Compute68CI(result_sb)

X0_U_val=[numpy.mean(result_u_val,axis=0)]+Compute68CI(result_u_val)
X0_D_val=[numpy.mean(result_d_val,axis=0)]+Compute68CI(result_d_val)
X0_S_val=[numpy.mean(result_s_val,axis=0)]+Compute68CI(result_s_val)

#%%
r=X0_U
[print("{"+str(xValues[i])+","+str(r[0][i])+","+str(r[1][i])+","+str(r[2][i])+"},") for i in range(len(xValues))]

r=X0_D
[print("{"+str(xValues[i])+","+str(r[0][i])+","+str(r[1][i])+","+str(r[2][i])+"},") for i in range(len(xValues))]

r=X0_Ub
[print("{"+str(xValues[i])+","+str(r[0][i])+","+str(r[1][i])+","+str(r[2][i])+"},") for i in range(len(xValues))]

r=X0_Db
[print("{"+str(xValues[i])+","+str(r[0][i])+","+str(r[1][i])+","+str(r[2][i])+"},") for i in range(len(xValues))]



#%%
#########################################################
## Resample data with N/2
#########################################################
def Resample(dd):
    #return numpy.random.choice(dd,size=int(numpy.floor(len(dd)/2)))
    return dd[numpy.random.choice(dd.shape[0], size=int(numpy.floor(len(dd)/2)))]

#########################################################
## Determine Mean, Mode, 68%CI by resampling 
#########################################################
alpha=68
def Compute68CI(dd):    
    lowers=[]
    uppers=[]    
    for i in range(1500):
        sample=Resample(numpy.array(dd))
        lowers.append(numpy.percentile(sample,(100-alpha)/2))
        uppers.append(numpy.percentile(sample,100-(100-alpha)/2))
    
    return [numpy.mean(lowers),numpy.mean(uppers)]

rSet.SetReplica(0)

xValues=[0.01+0.04*i for i in range(25)]

result=[]
f=-2
mu=10.
for x in xValues:
    ss=[]
    for r in range(rSetS.numberOfReplicas):
        rSetS.SetReplica(r)
        X0=harpy.get_SiversTMDPDF_g1(x, mu, 1)
        ss.append(X0[5+f])
        
    Xsivers=Compute68CI(ss)
    result.append([x,numpy.mean(ss),Xsivers[0],Xsivers[1]])
    print(x)
    
[print("{"+str(r[0])+","+str(r[1])+","+str(r[2])+","+str(r[3])+"},") for r in result]

#%%

import scipy.integrate as integrate

rSet.SetReplica(0)

xValues=[0.01+0.04*i for i in range(25)]

result=[]

mu=10.
ss_sea=[]
ss_u=[]
ss_d=[]

ss_tot=[]
for r in range(rSetS.numberOfReplicas):
    rSetS.SetReplica(r)    
    sea = integrate.quad(lambda x: harpy.get_SiversTMDPDF_g1(x, 10., 1)[5-1], 0.0001, 1.)
    d = integrate.quad(lambda x: harpy.get_SiversTMDPDF_g1(x, 10., 1)[5+1], 0.0001, 1.)
    u = integrate.quad(lambda x: harpy.get_SiversTMDPDF_g1(x, 10., 1)[5+2], 0.0001, 1.)
    ss_u.append(u)
    ss_d.append(d)
    ss_sea.append(sea)
    ss_tot.append(u+d+4*sea)
    

print("u: "+str(numpy.mean(ss_u))+str(Compute68CI(ss_u)))
print("d: "+str(numpy.mean(ss_d))+str(Compute68CI(ss_d)))
print("sea: "+str(numpy.mean(ss_sea))+str(Compute68CI(ss_sea)))
print("tot: "+str(numpy.mean(ss_tot))+str(Compute68CI(ss_tot)))
    