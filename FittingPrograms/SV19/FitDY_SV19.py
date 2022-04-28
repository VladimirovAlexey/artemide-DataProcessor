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
import time
import numpy
sys.path.append("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor")
import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet

MAINPATH="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"
#%%
#######################################
#Initialize artemide with a replica -2
#######################################
import harpy
path_to_constants=MAINPATH+"FittingPrograms/SV19/Constants-files/"
#harpy.initialize(path_to_constants+"DY_n3lo/const-NNPDF31_n3lo")
#harpy.initialize(path_to_constants+"DY_nnlo/const-HERA20_NNLO")
#harpy.initialize(path_to_constants+"DY_n3lo/const-HERA20_n3lo")
#harpy.initialize(path_to_constants+"DY_nnlo/const-MMHT14_NNLO")
#harpy.initialize(path_to_constants+"DY_nnlo/const-CT14_NNLO")
#harpy.initialize(path_to_constants+"DY_nnlo/const-PDF4LHC_NNLO")
harpy.setNPparameters([2.0340, 0.0299, 0.2512, 7.7572,334.6108, 2.4543,-4.8203, 0.1000,  0.0000])
#harpy.setNPparameters_TMDR(-2)
#harpy.setNPparameters_uTMDPDF(-2)

#%%
### read the list of files and return the list of DataSets
def loadThisData(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/DataLib/unpolDY/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

#%%
#################### LOG save function
LOGPATH=MAINPATH+"FittingPrograms/LOGS/"+"SV19(DY)["+time.ctime()+"].log"
def SaveToLog(logTitle,text):
    with open(LOGPATH, 'a') as file:
        file.write(time.ctime())
        file.write(' --> '+logTitle+'\n')
        file.write(text)
        file.write('\n \n \n')

#%%
##################Cut function
def cutFunc(p):
    par=1.0
    if p["type"]=="DY":
        if(p["xSec"]>0):
            err=numpy.sqrt(sum([i**2 for i in p["uncorrErr"]]))/p["xSec"]
        else:
            err=100.
        delta=p["<qT>"]/p["<Q>"]
        
        if(p["id"][0] == "E"):
            delta=p["<qT>"]/p["Q"][1] 
        
        
        if(p["id"][0:4] == "E605"):
            if(p["Q"][0]==10.5):#UPSILON resonance-bin
                return False , p
        elif(p["id"][0:4] == "E772"):
            if(p["Q"][0]<10):#these bins seems broken
                return False , p
        elif(p["id"][0:4] == "E615"):
            if(9<p["<Q>"]<11.2):#UPSILON resonance-bin
                return False , p
        elif(p["id"][0:4] == "E228"):
            if(9<p["<Q>"]<11):#UPSILON resonance-bin
                return False , p
        else:
            if(9<p["<Q>"]<11):#UPSILON resonance-bin
                return False , p
    
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
### Loading the data set
theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisData(['CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                      'A7-00y10', 'A7-10y20','A7-20y24', 
                      'A8-00y04', 'A8-04y08', 'A8-08y12', 'A8-12y16', 'A8-16y20', 'A8-20y24', 
                      'A8-46Q66', 'A8-116Q150', 
                      'CMS7', 'CMS8',
                      'LHCb7', 'LHCb8', 'LHCb13', 
                      'PHE200', 'E228-200', 'E228-300', 'E228-400', 
                      'E772',
                      'E605']))

setDY=theData.CutData(cutFunc) 

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setDY.sets])

#%%
### Loading the data set
theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisData([
    'CMS13-00y04','CMS13-04y08','CMS13-08y12','CMS13-12y16','CMS13-16y24']))

setDY=theData.CutData(cutFunc) 

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setDY.sets])


#%%
#######################################
# Test the main replica
#######################################
#harpy.setNPparameters_TMDR(0)
#harpy.setNPparameters_uTMDPDF(0)
#harpy.setNPparameters([3.0666, 0.0521, 0.1567, 24.0696, 3363.3436, 2.9974, -20.6362, 0.0000, 0.0000])

#harpy.setNPparameters([2.32699, 0.0210557,0.333229, 13.1943, 352.336, 2.08224, -10.7976, 0., 0.])  ##HERA
#harpy.setNPparameters([2.5973,  0.0236,  0.3395,   14.1241, 423.8772,1.8464,  -10.4110, 0., 0.])  ##HERA N3LL

#harpy.setNPparameters([2.29, 0.0222,0.324, 13.2, 356., 2.05, -10.4, 0., 0.])  ##HERA (paper)
#harpy.setNPparameters([1.94, 0.0335,0.326, 10.1, 273., 1.70, -6.5, 0., 0.])  ##HERA N3LL (paper)

harpy.setNPparameters([1.89278, 0.0291215,0.273023, 8.96853, 353.824, 2.50639, -5.92725, 0., 0.]) ##NNPDF
#harpy.setNPparameters([1.90812, 0.0290096, 0.273555, 8.99327, 356.011, 2.51554, -5.97051, 0., 0.]) ##NNPDF N3LL

#harpy.setNPparameters([1.86, 0.0296, 0.253, 9.0, 347., 2.48, -5.7, 0., 0.]) ##NNPDF NNLL (paper)
#harpy.setNPparameters([1.62, 0.0342, 0.282, 9.7, 317., 2.42, -6.1, 0., 0.]) ##NNPDF N3LL (paper)

#harpy.setNPparameters([2.34665, 0.022735,0.27706, 24.8883, 1241.34, 2.66869, -23.7589, 0., 0.]) ##CT14
#harpy.setNPparameters([1.54728,  0.04678,0.1982,   26.49689,2727.9766,    3.00668, -23.54749,   0.,0.]) ##MMHT14
#harpy.setNPparameters([1.92659, 0.036548,0.218079, 17.9138, 926.078, 2.5431, -15.5469, 0. ,0.]) ##PDF4LHC

#harpy.setNPparameters([2.35, 0.0227,0.277, 24.9, 1241., 2.67, -23.8, 0., 0.]) ##CT14 (paper)
#harpy.setNPparameters([1.55,  0.0470,0.198,   26.4,2680.,    3.01, -23.4,   0.,0.]) ##MMHT14 (paper)
#harpy.setNPparameters([1.93, 0.0366,0.218, 17.9, 926., 2.54, -15.5, 0. ,0.]) ##PDF4LHC (paper)

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)
    
#%%
#######################################
# Minimisation
#######################################
totalN=setDY.numberOfPoints

def chi_2(x):
    startT=time.time()
    harpy.setNPparameters(x)
    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
    
    ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)
    
    cc=ccDY2/totalN
    endT=time.time()
    print(':->',cc,'       t=',endT-startT)
    return ccDY2

#%%

from iminuit import Minuit

#initialValues=(2.29477, 0.022191,0.324112, 13.1774, 356.124, 2.05101, -10.4468, 0., 0.)  ### HERA

initialValues=(1.89278, 0.0291215,0.273023, 8.96853, 353.824, 2.50639, -5.92725, 0., 0.)  ### NNPDF

#initialValues=(1.54728,  0.04678,0.1982,   26.49689,2727.9766,    3.00668, -23.54749,   0.,0.) ##MMHT14
#initialValues=(2.34665, 0.022735,0.27706, 24.8883, 1241.34, 2.66869, -23.7589, 0., 0.) ##CT14
#initialValues=(1.92659, 0.036548,0.218079, 17.9138, 926.078, 2.5431, -15.5469, 0. ,0.) ##PDF4LHC


initialErrors=(0.1,       0.01,     0.05,      0.1,      10.,       0.05,    0.05,    0.1, 0.1)
searchLimits=((1.4,4.5), (0.001,5.0),(0.0,4.0),(8.,25.0),(0.0,400),(0.,5), None,   None, None)
parametersToMinimize=(False,     False,    False,    False,    False,     False,  False, True, True)

m = Minuit.from_array_func(chi_2, initialValues,
      error=initialErrors, limit=searchLimits, fix=parametersToMinimize, errordef=1)

#m.get_param_states()

m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1

SaveToLog("MINIMIZATION STARTED",str(m.params))
#%%


m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1

m.migrad()

print(m.params)

SaveToLog("MINIMIZATION FINISHED",str(m.params))
SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))
#%%

m.hesse()

print(m.params)

SaveToLog("HESSE FINISHED",str(m.params))
SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))
#%%

m.minos()

print(m.params)

SaveToLog("MINOS FINISHED",str(m.params))
SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))
