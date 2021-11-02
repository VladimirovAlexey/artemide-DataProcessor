#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 16:40:59 2019

@author: vla18041
"""

#######################################
# importing libraries
#######################################

import sys
import numpy
sys.path.append("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor")
import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet

MAINPATH="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"

tt=''
#tt='all=0'

case='h-.case4'

#%%
#######################################
#Initialize artemide
#######################################
import harpy
if(tt==''):
    path_to_constants=MAINPATH+"FittingPrograms/SV19/Constants-files/"+"DY+SIDIS_nnlo/const-DY+SIDIS_NNPDF31+DSS_nnlo"
elif(tt=='all=0'):
    path_to_constants=MAINPATH+"FittingPrograms/SV19/Constants-files/"+"DY+SIDIS_nnlo_all=0/const-DY+SIDIS_NNPDF31+DSS_nnlo_all=0"

harpy.initialize(path_to_constants)

harpy.setNPparameters_TMDR([1.93, 0.0434])
harpy.setNPparameters_uTMDPDF([0.253434, 9.04351, 346.999, 2.47992, -5.69988, 0.1, 0.])
harpy.setNPparameters_uTMDFF([0.264,0.479,0.459,0.539]) 

#%%
### read the list of files and return the list of DataSets

def loadThisDataSIDIS(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/FittingPrograms/SV19/PredictionCOMPASS/Bins/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

#%%
##################Cut function
def GetDelta(p):
    if(tt==''):
        gamma2=(2.0*p["M_target"]*p["<x>"]/p["<Q>"])**2
        rho2=(p["M_product"]/p["<z>"]/(p["<Q>"]))**2
        qT=p["<pT>"]/p["<z>"]*numpy.sqrt((1+gamma2)/(1-gamma2*rho2))
        #qT=p["<pT>"]/p["<z>"]
        return qT/(p["<Q>"])
    else:
        return p["<pT>"]/p["<z>"]/(p["<Q>"])

def cutFunc(p):   
   
    delta=GetDelta(p)
    
    if(tt==''):
        rho2=(p["M_product"]/p["<z>"]/(p["<Q>"]))**2
        
        ### compute the largest possible qT (approximate)
        gamma2WORST=(2.0*p["M_target"]*p["x"][1]/p["<Q>"])**2
        # it is definitely not a TMD point
        if gamma2WORST*rho2>1:
            return False , p
        qTWORST=p["pT"][1]/p["z"][0]*numpy.sqrt((1+gamma2WORST)/(1-gamma2WORST*rho2))
    
        ## drop if qT>Q/2
        if qTWORST>p["<Q>"]/2:
            return False , p
    
#    return delta<0.5 and p.qT_avarage<80
    return delta<0.75 , p

#%%
### Loading the SIDIS data set
theData=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisDataSIDIS([                      
                      'compass.d.'+case]))

setSIDIS=theData.CutData(cutFunc) 

print('Loaded ', setSIDIS.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setSIDIS.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setSIDIS.sets])



#%%

import DataProcessor.ArtemideReplicaSet

if(tt==''):
    rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/Replicas/DY+SIDIS/DY+SIDIS_NNPDF31+DSS_nnlo.rep")
elif(tt=='all=0'):
    rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/Replicas/DY+SIDIS/SV19_nnlo_all=0.rep")
rSet.SetReplica()


#%%
#### Compute avarage pT2 using avarage kinematics
def AvaragePT2(p):
    import copy
    
    xx=[]
    pp=[]
    
    for i in range(10):
        pN=copy.deepcopy(p)
        pN["pT"]=[(p["pT"][1]-p["pT"][0])*i/10+p["pT"][0],(p["pT"][1]-p["pT"][0])*(i+1)/10+p["pT"][0]]    
        pN["<pT>"]=(p["pT"][1]-p["pT"][0])*(i+0.5)/10+p["pT"][0]
        pN["thFactor"]=1/(pN["pT"][1]**2-pN["pT"][0]**2)/(pN["z"][1]-pN["z"][0])
        pp.append((p["pT"][1]-p["pT"][0])*(i+0.5)/10+p["pT"][0])
        xx.append(DataProcessor.harpyInterface.ComputeXSec(pN,method="binless"))
    
    pt2=0
    mm=0
    for i in range(10):
        pt2+=pp[i]**2*xx[i]
        mm+=xx[i]
    
    #print(xx)
    
    return pt2/mm
    
print(AvaragePT2(setSIDIS.points[4]))

#%%

central=DataProcessor.harpyInterface.ComputeXSec(setSIDIS)

pt2=[AvaragePT2(p) for p in setSIDIS.points]

replicas=[]
for i in range(rSet.numberOfReplicas):
    rSet.SetReplica(i)
    replicas.append(DataProcessor.harpyInterface.ComputeXSec(setSIDIS))

std=2.*numpy.std(replicas,axis=0)
    
#%%

path_to_save="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/FittingPrograms/SV19/PredictionCOMPASS/Prediction/"

if(case[-1]=='1' or case[-1]=='2' or case[-1]=='5'):

    with open(path_to_save+case+tt+".dat",'w') as file:
        file.write("#  Artemide prediction based on SV19 for "+case+". \n")
        file.write("#  xSec is the prediction for unpolarized cross-section dSigma/dz/dpT^2 in the given bin \n")
        file.write("#  Err is uncertanty of SV19 fit \n")
        file.write("#  The column delta is qT/Q-parameter of TMD factorization. Computed at the central parameters of the bins. Only qT/Q<0.25 are trustfull. \
                   Beyond this number TMD-prediction deviates from the xSec. In the table the values up to delta=0.75 are included. Most part of COMPASS data has delta>1  \n")    
        file.write("#  Q2min, Q2max, xMin, xMax, zMin, zMax, pT2min, pT2max, <pT2>, xSec, err, delta \n")
        for i in range(setSIDIS.numberOfPoints):
            p=setSIDIS.points[i]
            file.write("{:.2f}".format(p["Q"][0]**2)+', '+"{:.2f}".format(p["Q"][1]**2)+', '+\
                       str(p["x"][0])+', '+str(p["x"][1])+', '+\
                       str(p["z"][0])+', '+str(p["z"][1])+', '+\
                       '{:2f}'.format(p["pT"][0]**2)+', '+'{:2f}'.format(p["pT"][1]**2)+', '+'{:2f}'.format(pt2[i])+', '\
                       '{:g}'.format(central[i])+', '+'{:g}'.format(std[i])+', '+'{:2f}'.format(GetDelta(p))+" \n")
                
if(case[-1]=='3' or case[-1]=='4'):
    with open(path_to_save+case+tt+"-lowW.dat",'w') as file:
        file.write("#  Artemide prediction based on SV19 for "+case+" at W<12. \n")
        file.write("#  xSec is the prediction for unpolarized cross-section dSigma/dz/dpT^2 in the given bin \n")
        file.write("#  Err is uncertanty of SV19 fit \n")
        file.write("#  The column delta is qT/Q-parameter of TMD factorization. Computed at the central parameters of the bins. Only qT/Q<0.25 are trustfull. \
                   Beyond this number TMD-prediction deviates from the xSec. In the table the values up to delta=0.75 are included. Most part of COMPASS data has delta>1  \n")    
        file.write("#  Q2min, Q2max, xMin, xMax, zMin, zMax, pT2min, pT2max, <pT2>, xSec, err, delta \n")
        for i in range(setSIDIS.numberOfPoints):
            p=setSIDIS.points[i]
            if(p["cutParams"][3]<=145.):
                file.write("{:.2f}".format(p["Q"][0]**2)+', '+"{:.2f}".format(p["Q"][1]**2)+', '+\
                           str(p["x"][0])+', '+str(p["x"][1])+', '+\
                           str(p["z"][0])+', '+str(p["z"][1])+', '+\
                           '{:2f}'.format(p["pT"][0]**2)+', '+'{:2f}'.format(p["pT"][1]**2)+', '+'{:2f}'.format(pt2[i])+', '\
                           '{:g}'.format(central[i])+', '+'{:g}'.format(std[i])+', '+'{:2f}'.format(GetDelta(p))+" \n")
                    

    with open(path_to_save+case+tt+"-highW.dat",'w') as file:
        file.write("#  Artemide prediction based on SV19 for "+case+" at W<12. \n")
        file.write("#  xSec is the prediction for unpolarized cross-section dSigma/dz/dpT^2 in the given bin \n")
        file.write("#  Err is uncertanty of SV19 fit \n")
        file.write("#  The column delta is qT/Q-parameter of TMD factorization. Computed at the central parameters of the bins. Only qT/Q<0.25 are trustfull. \
                   Beyond this number TMD-prediction deviates from the xSec. In the table the values up to delta=0.75 are included. Most part of COMPASS data has delta>1  \n")    
        file.write("#  Q2min, Q2max, xMin, xMax, zMin, zMax, pT2min, pT2max, <pT2>, xSec, err, delta \n")
        for i in range(setSIDIS.numberOfPoints):
            p=setSIDIS.points[i]
            if(p["cutParams"][2]>=140.):
                file.write("{:.2f}".format(p["Q"][0]**2)+', '+"{:.2f}".format(p["Q"][1]**2)+', '+\
                           str(p["x"][0])+', '+str(p["x"][1])+', '+\
                           str(p["z"][0])+', '+str(p["z"][1])+', '+\
                           '{:2f}'.format(p["pT"][0]**2)+', '+'{:2f}'.format(p["pT"][1]**2)+', '+'{:2f}'.format(pt2[i])+', '\
                           '{:g}'.format(central[i])+', '+'{:g}'.format(std[i])+', '+'{:2f}'.format(GetDelta(p))+" \n")
