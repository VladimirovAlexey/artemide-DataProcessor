#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 16:40:59 2019

@author: vla18041
"""

#######################################
# importing libraries
#######################################
import os
MAINPATH=os.path.join(os.path.dirname(__file__),"..","..")

import sys
import time
import numpy

sys.path.append(MAINPATH)

import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet

useOrder="n3lo"
usePDF="DSSV"
#usePDF="NNPDF"
#%%
#######################################
#Initialize artemide
#######################################
import harpy
path_to_constants=MAINPATH+"/FittingPrograms/WGT22/ConstantsFiles/"

if(useOrder=="nnlo"):
    harpy.initialize(path_to_constants+"const-WGT22_nnlo_"+usePDF)
    
    #### All=0 Case
    harpy.setNPparameters_TMDR([2., 0.0398333])
    harpy.setNPparameters_uTMDPDF([0.185239, 6.22706, 580.946, 2.44166, -2.53161, 0.,  0.0014, 0.442, 4.14])
    harpy.setNPparameters_uTMDFF([0.279443, 0.460015, 0.435955, 0.551302])
    harpy.setNPparameters_wgtTMDPDF([0.2, 1.0])
    
elif(useOrder=="n3lo"):
    harpy.initialize(path_to_constants+"const-WGT22_n3lo_"+usePDF)
    #### All=0 Case n3lo
    harpy.setNPparameters_TMDR([2.0, 0.04843])
    harpy.setNPparameters_uTMDPDF([0.1425, 4.8199, 580.9, 2.3889, -1.0683, 0.0,  0.0014, 0.442, 4.14])
    harpy.setNPparameters_uTMDFF([0.2797, 0.4469, 0.43215, 0.63246])
    harpy.setNPparameters_wgtTMDPDF([1.5, 1.0])
#%%
import DataProcessor.ArtemideReplicaSet
rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/wgt22/REPS/"
                                                  +"wgt22-"+usePDF+".rep")### SIDIS+DY case

#%%
def loadThisData(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/DataLib/wgt/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

def loadThisDataSIVERS(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/DataLib/Sivers/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

#%%
##################Cut function
def cutFunc(p):
    import copy 
    
    if p["type"]=="SIDIS":   
        deltaTEST=0.35
        delta=p["<pT>"]/p["<z>"]/p["<Q>"]  
        
    if "JLab6" in p["id"]:
        return False, p
    
    
    if delta<deltaTEST:
        pNew=copy.deepcopy(p)    
        pNew["process"]=pNew["weightProcess"]
        if p["type"]=="SIDIS":
            normX=DataProcessor.harpyInterface.ComputeXSec(pNew,method="central")        
        elif p["type"]=="DY":
            normX=DataProcessor.harpyInterface.ComputeXSec(pNew)        
        else:
            print("Are you crazy?")
        p["thFactor"]=p["thFactor"]/normX        
    
        
#    return delta<0.5 and p.qT_avarage<80
    return delta<deltaTEST and p["<Q>"]>1.41, p

#%%
##################Cut function0
def cutFunc0(p):
    import copy 
    
    if p["type"]=="SIDIS":   
        deltaTEST=0.75
        delta=p["<pT>"]/p["<z>"]/p["<Q>"]        
    
    
    if delta<deltaTEST:
        pNew=copy.deepcopy(p)    
        pNew["process"]=pNew["weightProcess"]
        if p["type"]=="SIDIS":
            normX=DataProcessor.harpyInterface.ComputeXSec(pNew,method="central")        
        elif p["type"]=="DY":
            normX=DataProcessor.harpyInterface.ComputeXSec(pNew)        
        else:
            print("Are you crazy?")
        p["thFactor"]=p["thFactor"]/normX        
    
        
#    return delta<0.5 and p.qT_avarage<80
    return delta<deltaTEST, p

#%%
##################Cut function
def cutFuncDY(p):
    import copy
    
    pNew=copy.deepcopy(p)    
    pNew["process"]=pNew["weightProcess"]
    if p["type"]=="SIDIS":
        normX=DataProcessor.harpyInterface.ComputeXSec(pNew,method="central")        
    elif p["type"]=="DY":
        normX=DataProcessor.harpyInterface.ComputeXSec(pNew)        
    else:
        print("Are you crazy?")
    p["thFactor"]=p["thFactor"]/normX    
    
    if(p["process"][2]==10005):p["process"][2]=13200
    if(p["process"][2]==10007):p["process"][2]=13201
    if(p["process"][2]==10008):p["process"][2]=13202

#    return delta<0.5 and p.qT_avarage<80
    return True, p

#%%
### Loading the SIDIS data set
theData=DataProcessor.DataMultiSet.DataMultiSet("ALTset",loadThisData([
                      'hermes3D.ALT.pi+','hermes3D.ALT.pi-',
                      'hermes3D.ALT.k+','hermes3D.ALT.k-',
                      'compass16.ALT.h+.2<z.dpt','compass16.ALT.h-.2<z.dpt',
                      'compass16.ALT.h+.2<z.dz','compass16.ALT.h-.2<z.dz',
                      'compass16.ALT.h+.2<z.dx','compass16.ALT.h-.2<z.dx',
                      'JLab6.ALT.pi+','JLab6.ALT.pi-'
                      ]))

setALT=theData.CutData(cutFunc) 
setALTfull=theData.CutData(cutFunc0) 

theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataSIVERS([
                    'star.sivers.W+.dqT','star.sivers.W-.dqT',
                    'star.sivers.W+.dy','star.sivers.W-.dy',
                    'star.sivers.Z'
                    ]))

setDY=theData.CutData(cutFuncDY) 

print('Loaded ', setALT.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setALT.sets]), 'points.')
print('Loaded ', setALTfull.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setALTfull.sets]), 'points.')
print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')

#SaveToLog('Loaded '+ str(setDY.numberOfSets) + ' data sets with '+str(sum([i.numberOfPoints for i in setDY.sets])) + ' points. \n'
#+'Loaded experiments are '+str([i.name for i in setDY.sets]))
#%%
rSet.SetReplica()

DataProcessor.harpyInterface.PrintChi2Table(setALT,printDecomposedChi2=True,method="central")

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)

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

#%%
listXsec=[]
listIn=[p["id"] for p in setALT.points]

for p in setALTfull.points:
    ll=[]
    for i in range(rSet.numberOfReplicas):
        rSet.SetReplica(i+1)
        ll.append(DataProcessor.harpyInterface.ComputeXSec(p,method="central"))
    mm=numpy.mean(ll)
    ci=Compute68CI(ll)
    fitted=(p["id"] in listIn)
    listXsec.append(p["id"]+", {:9.6f}, {:9.6f}, {:9.6f}, ".format(mm,ci[0],ci[1])+str(fitted))
    
for p in setDY.points:
    ll=[]
    print(p["id"])
    for i in range(rSet.numberOfReplicas):
        rSet.SetReplica(i+1)
        ll.append(DataProcessor.harpyInterface.ComputeXSec(p))    
    mm=numpy.mean(ll)
    ci=Compute68CI(ll)
    fitted=(p["id"] in listIn)
    listXsec.append(p["id"]+", {:9.6f}, {:9.6f}, {:9.6f}, ".format(mm,ci[0],ci[1])+str(fitted))

f=open("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/wgt22/xSecs/xSec_"+usePDF+".csv","w")
for x in listXsec:
    f.write(x+"\n")
f.close()

#%%
wgtListU=[]
wgtListD=[]
wgtListS=[]
wgtListUb=[]
wgtListDb=[]
for i in range(50):    
    x=0.02*(i+1)
    print("x=",x)
    for j in range(20):
        b=(0.1*j)**2
        llU=[]
        llD=[]
        llS=[]
        llUb=[]
        llDb=[]
        for k in range(rSet.numberOfReplicas):
            rSet.SetReplica(k+1)
            rr=harpy.get_wgtTMDPDF(x,b,1)
            llU.append(rr[7])
            llD.append(rr[6])
            llS.append(rr[8])
            llUb.append(rr[3])
            llDb.append(rr[4])
        mm=numpy.mean(llU)
        ci=Compute68CI(llU)
        wgtListU.append("{:4.2f}, {:4.2f}, {:9.6f}, {:9.6f}, {:9.6f}, ".format(x,b,mm,ci[0],ci[1]))
        
        mm=numpy.mean(llD)
        ci=Compute68CI(llD)
        wgtListD.append("{:4.2f}, {:4.2f}, {:9.6f}, {:9.6f}, {:9.6f}, ".format(x,b,mm,ci[0],ci[1]))
        
        mm=numpy.mean(llS)
        ci=Compute68CI(llS)
        wgtListS.append("{:4.2f}, {:4.2f}, {:9.6f}, {:9.6f}, {:9.6f}, ".format(x,b,mm,ci[0],ci[1]))
        
        mm=numpy.mean(llUb)
        ci=Compute68CI(llUb)
        wgtListUb.append("{:4.2f}, {:4.2f}, {:9.6f}, {:9.6f}, {:9.6f}, ".format(x,b,mm,ci[0],ci[1]))
        
        mm=numpy.mean(llDb)
        ci=Compute68CI(llDb)
        wgtListDb.append("{:4.2f}, {:4.2f}, {:9.6f}, {:9.6f}, {:9.6f}, ".format(x,b,mm,ci[0],ci[1]))
        
f=open("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/wgt22/wgt-function/wgt_U_"+usePDF+".csv","w")
for x in wgtListU:
    f.write(x+"\n")
f.close()

f=open("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/wgt22/wgt-function/wgt_D_"+usePDF+".csv","w")
for x in wgtListD:
    f.write(x+"\n")
f.close()

f=open("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/wgt22/wgt-function/wgt_S_"+usePDF+".csv","w")
for x in wgtListS:
    f.write(x+"\n")
f.close()

f=open("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/wgt22/wgt-function/wgt_Ub_"+usePDF+".csv","w")
for x in wgtListUb:
    f.write(x+"\n")
f.close()

f=open("/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/wgt22/wgt-function/wgt_Db_"+usePDF+".csv","w")
for x in wgtListDb:
    f.write(x+"\n")
f.close()
