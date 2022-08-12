#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 10:38:36 2022

@author: vla18041
"""


import sys
import time
import numpy
sys.path.append("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor")
import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet

MAINPATH="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"
#%%
#######################################
#Initialize artemide
#######################################
import harpy
path_to_constants=MAINPATH+"FittingPrograms/MVZ22/ConstantsFiles/"+"const-DY+SIDIS_NNPDF31+DSS_nnlo"


harpy.initialize(path_to_constants)

harpy.setNPparameters_TMDR([1.93, 0.0434])
harpy.setNPparameters_uTMDPDF([0.253434, 9.04351, 346.999, 2.47992, -5.69988, 0.1, 0.])
harpy.setNPparameters_uTMDFF([0.264,0.479,0.459,0.539]) 

originalSV19TMDR=[1.93, 0.0434]
originalSV19FF=[0.264,0.479,0.459,0.539]

#%%
def loadThisDataSIDIS(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/JLab_Impact/PseudoData/jlab22/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

#%%
##################Cut function
def cutFunc(p):
    par=1.0   
    if p["type"]=="SIDIS":        
        if p["<z>"]>0.8:
            return False , p
        ## bins with low z drop
        if p["<z>"]<0.2:
            return False , p
        
        par=1.0
        if p["xSec"]<0.00000001:
            err=1
            delta=1
        else:
            ##############3 I MULTIPLY THE ERROR BY 100 (so it does not affect the cuts)
            err=10000#*numpy.sqrt(p.uncorrErrorsSquare)/p.xSec    
            gamma2=(2.0*p["M_target"]*p["<x>"]/p["<Q>"])**2
            rho2=(p["M_product"]/p["<z>"]/(p["<Q>"]))**2
            qT=p["<pT>"]/p["<z>"]*numpy.sqrt((1+gamma2)/(1-gamma2*rho2))
            delta=qT/(p["<Q>"])
            
            ### compute the largest possible qT (approximate)
            gamma2WORST=(2.0*p["M_target"]*p["x"][1]/p["<Q>"])**2
            # it is definitely not a TMD point
            if gamma2WORST*rho2>1:
                return False , p
            qTWORST=p["pT"][1]/p["z"][0]*numpy.sqrt((1+gamma2WORST)/(1-gamma2WORST*rho2))
    
            ## drop if qT>Q/2
            if qTWORST>p["<Q>"]/2:
                return False , p
    
        ### drop Q<2
        if p["<Q>"]<2 :
            return False , p
    
#    return delta<0.5 and p.qT_avarage<80
    return (delta<0.1 or (delta<0.25 and par/err*delta**2<1)) , p

#%%
##################Cut function
def cutFuncJLAB(p):
    if p["<z>"]>0.8:
        return False , p
    ## bins with low z drop
    if p["<z>"]<0.2:
        return False , p
    
    if p["xSec"]<0.00000001:
        delta=1
    else:
        gamma2=(2.0*p["M_target"]*p["<x>"]/p["<Q>"])**2
        rho2=(p["M_product"]/p["<z>"]/(p["<Q>"]))**2
        qT=p["<pT>"]/p["<z>"]*numpy.sqrt((1+gamma2)/(1-gamma2*rho2))
        delta=qT/(p["<Q>"])
        
        ### compute the largest possible qT (approximate)
        gamma2WORST=(2.0*p["M_target"]*p["x"][1]/p["<Q>"])**2
        # it is definitely not a TMD point
        if gamma2WORST*rho2>1:
            return False , p
        qTWORST=p["pT"][1]/p["z"][0]*numpy.sqrt((1+gamma2WORST)/(1-gamma2WORST*rho2))

        ## drop if qT>Q/2
        if qTWORST>p["<Q>"]/2:
            return False , p

    ### drop Q<2
    if p["<Q>"]<2 :
        return False , p  
    ### value of cross-section
    xSec=DataProcessor.harpyInterface.ComputeXSec(p,method="central")
    ### replace xSec by theory prediction
    p["xSec"]=xSec
     ### rescale errors by new theory
    p["uncorrErr"][0]=p["uncorrErr"][0]*numpy.abs(p["xSec"])
    p["uncorrErr"][1]=p["uncorrErr"][1]*numpy.abs(p["xSec"])
    
#    return delta<0.5 and p.qT_avarage<80
    return delta<0.25, p

#%%
### Loading the SIDIS data set
theData=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisDataSIDIS([
                      'jlab22.gen.pi+.join521',
                      'jlab22.gen.pi-.join521',
                      'jlab22.gen.k+.join521']))

setJLAB=theData.CutData(cutFuncJLAB) 

print('Loaded ', setJLAB.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setJLAB.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setJLAB.sets])

#%%
DataProcessor.harpyInterface.PrintChi2Table(setJLAB,method="central",printSysShift=False)


#%%
#######################################
# Minimisation
#######################################
totalN=setJLAB.numberOfPoints

def chi_2(x):
    startT=time.time()
    harpy.setNPparameters_uTMDFF(x)
    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
    
    ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(setJLAB,method="central")
    
    cc=ccSIDIS2/totalN
    endT=time.time()
    print(':->',cc,'       t=',endT-startT)
    return ccSIDIS2

#%%
from iminuit import Minuit


initialValues=originalSV19FF

initialErrors=(0.1,0.1,0.1,0.1)
searchLimits=((0.,None),(0.,None),(0.,None),None)

parametersToMinimize=(False,False,False,False)

m = Minuit(chi_2, initialValues)

m.errors=initialErrors
m.limits=searchLimits
m.fixed=parametersToMinimize
m.errordef=1

print(m.params)

#m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1

#%%

m.hesse()

print(m.params)
print(m.covariance)

#%%
m.minos()
print(m.params)
print(m.covariance)

#%%

# def MinForReplica():
    
    
#     def repchi_2(x):        
#         startT=time.time()
#         harpy.setNPparameters_uTMDFF(x)
#         print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
                    
#         ccSIDIS2,cc3=DataProcessor.harpyInterface.ComputeChi2(repDataSIDIS,method="central")
#         cc=ccSIDIS2/totalNnew
        
#         endT=time.time()
#         print(':->',cc,'       t=',endT-startT)
#         return ccSIDIS2
    
#     repDataSIDIS=setJLAB.GenerateReplica()
#     totalNnew=repDataSIDIS.numberOfPoints    
    
#     localM = Minuit(repchi_2, initialValues)
    
#     localM.errors=initialErrors
#     localM.limits=searchLimits
#     localM.fixed=parametersToMinimize
#     localM.errordef=1    
#     localM.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
#     localM.strategy=1
    
#     localM.migrad()
    
#     chi2Central=chi_2(list(localM.values))
    
#     return [localM.fval,chi2Central,list(localM.values)]

# #%%
# #
# # Generate pseudo data and minimise   100 times
# #
# numOfReplicas=150
# REPPATH="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/JLab_Impact/Output/unpolarized/replicas1.dat"
# for i in range(numOfReplicas):
#     print('---------------------------------------------------------------')
#     print('------------REPLICA ',i,'/',numOfReplicas,'--------------------')
#     print('---------------------------------------------------------------')
#     repRes=MinForReplica()
#     print(repRes)
#     f=open(REPPATH,"a+")
#     print('SAVING >>  ',f.name)
#     f.write(str(repRes)+"\n")
#     f.close()