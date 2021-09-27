#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 14:35:35 2021

@author: vla18041
"""

##############################
# Ploting SV21 various cases
##############################

#%%
#######################################################################
# Global parameter of a run
#######################################################################

#PDFinUse="HERA20"
PDFinUse="NNPDF31"
#PDFinUse="CT18"
#PDFinUse="MSHT20"
#PDFinUse="CJ15"

#######################################
# Paths
#######################################
PathToHarpy="/home/vla18041/LinkData2/arTeMiDe_Repository/artemide-ForPDF/harpy"
#PathToHarpy="/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/harpy"
PathToDataProcessor="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"
PathToDataLibrary=PathToDataProcessor+"DataLib/unpolDY/"
PathToSavings="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV21/"+PDFinUse+"/"
PathToConstantsFile=PathToDataProcessor+"/FittingPrograms/PDF-TMD/Constants-files/const-"+PDFinUse+"_NNLO_12p"

import sys
sys.path.remove('/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/harpy')
sys.path.append(PathToHarpy)
sys.path.append(PathToDataProcessor)


#%%
#######################################
# importing libraries
#######################################
import numpy
import time

import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet
import DataProcessor.ArtemideReplicaSet

#%%
#######################################
#Initialize artemide
#######################################
import harpy

harpy.initialize(PathToConstantsFile)

# rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
#                             "/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/REPS/SV21-"+PDFinUse+"-nnlo.rep")
    
rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/REPS/SV21-"+PDFinUse+"-nnlo-EXP.rep")
rSet.SetReplica()

#%%
#######################################
# read the list of files and return the list of DataSets
#######################################
def loadThisData(listOfNames):    
    import DataProcessor.DataSet
    
    dataCollection=[]
    for name in listOfNames:
        if( name==''): continue
        loadedData=DataProcessor.DataSet.LoadCSV(PathToDataLibrary+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection

#%%
##################Cut functions
includePoints=[
    'CDF2.0', 'CDF2.1', 'CDF2.2', 'CDF2.3', 'CDF2.4', 'CDF2.5', 'CDF2.6', 'CDF2.7',
    'CDF2.8', 'CDF2.9', 'CDF2.10', 'CDF2.11', 'CDF2.12', 'CDF2.13', 'CDF2.14',
    'D02m.0', 'D02m.1', 'D02m.2',
    'A8-00y04.0', 'A8-00y04.1', 'A8-00y04.2', 'A8-00y04.3', 'A8-00y04.4', 'A8-04y08.0',
    'A8-04y08.1', 'A8-04y08.2', 'A8-04y08.3', 'A8-04y08.4', 'A8-08y12.0', 'A8-08y12.1',
    'A8-08y12.2', 'A8-08y12.3', 'A8-08y12.4', 'A8-12y16.0', 'A8-12y16.1', 'A8-12y16.2',
    'A8-12y16.3', 'A8-12y16.4', 'A8-16y20.0', 'A8-16y20.1', 'A8-16y20.2', 'A8-16y20.3',
    'A8-16y20.4', 'A8-20y24.0', 'A8-20y24.1', 'A8-20y24.2', 'A8-20y24.3', 'A8-20y24.4',
    'A8-46Q66.0', 'A8-46Q66.1', 'A8-46Q66.2',
    'LHCb7.0', 'LHCb7.1', 'LHCb7.5', 'LHCb7.6', 'LHCb8.0', 'LHCb8.1',
    'LHCb8.2', 'LHCb8.3', 'LHCb8.4', 'LHCb8.5', 'LHCb8.6', 'PHE200.0', 'PHE200.1',
    'E228-200.4Q5.0', 'E228-200.4Q5.1', 'E228-200.4Q5.2', 'E228-200.4Q5.3', 'E228-200.4Q5.4',
    'E228-200.4Q5.5', 'E228-200.5Q6.0', 'E228-200.5Q6.1', 'E228-200.5Q6.2', 'E228-200.5Q6.3',
    'E228-200.5Q6.4', 'E228-200.5Q6.5', 'E228-200.5Q6.6', 'E228-200.6Q7.0', 'E228-200.6Q7.1',
    'E228-200.6Q7.2', 'E228-200.6Q7.3', 'E228-200.6Q7.4', 'E228-200.6Q7.5', 'E228-200.6Q7.6',
    'E228-200.6Q7.7', 'E228-200.6Q7.8', 'E228-200.7Q8.0', 'E228-200.7Q8.1', 'E228-200.7Q8.2',
    'E228-200.7Q8.3', 'E228-200.7Q8.4', 'E228-200.7Q8.5', 'E228-200.7Q8.6', 'E228-200.7Q8.7',
    'E228-200.7Q8.8', 'E228-200.7Q8.9', 'E228-200.8Q9.0', 'E228-200.8Q9.1', 'E228-200.8Q9.2',
    'E228-200.8Q9.3', 'E228-200.8Q9.4', 'E228-200.8Q9.7', 'E228-200.8Q9.8', 'E228-300.4Q5.0<u',
    'E228-300.4Q5.1<u', 'E228-300.4Q5.2<u', 'E228-300.4Q5.3<u', 'E228-300.4Q5.4<u', 'E228-300.4Q5.5<u',
    'E228-300.5Q6.0<u', 'E228-300.5Q6.1<u', 'E228-300.5Q6.2<u', 'E228-300.5Q6.3<u', 'E228-300.5Q6.4<u',
    'E228-300.5Q6.5<u', 'E228-300.5Q6.6<u', 'E228-300.6Q7.0<u', 'E228-300.6Q7.1<u', 'E228-300.6Q7.2<u',
    'E228-300.6Q7.3<u', 'E228-300.6Q7.4<u', 'E228-300.6Q7.5<u', 'E228-300.6Q7.6<u', 'E228-300.6Q7.7<u',
    'E228-300.6Q7.8<u', 'E228-300.7Q8.0<u', 'E228-300.7Q8.1<u', 'E228-300.7Q8.2<u', 'E228-300.7Q8.3<u',
    'E228-300.7Q8.4<u', 'E228-300.7Q8.5<u', 'E228-300.7Q8.6<u', 'E228-300.7Q8.7<u', 'E228-300.7Q8.8<u',
    'E228-300.7Q8.9<u', 'E228-300.8Q9.0<u', 'E228-300.8Q9.1<u', 'E228-300.8Q9.2<u', 'E228-300.8Q9.3<u',
    'E228-300.8Q9.4<u', 'E228-300.8Q9.5<u', 'E228-300.8Q9.6<u', 'E228-300.8Q9.7<u', 'E228-300.8Q9.8<u',
    'E228-300.8Q9.9<u', 'E228-300.8Q9.10<u', 'E228-300.11Q12.0>u', 'E228-300.11Q12.1>u',
    'E228-300.11Q12.2>u', 'E228-300.11Q12.3>u', 'E228-300.11Q12.4>u', 'E228-300.11Q12.5>u',
    'E228-300.11Q12.6>u', 'E228-300.11Q12.7>u', 'E228-300.11Q12.8>u', 'E228-300.11Q12.11>u',
    'E228-400.5Q6.0<u', 'E228-400.5Q6.1<u', 'E228-400.5Q6.2<u', 'E228-400.5Q6.3<u', 'E228-400.5Q6.4<u',
    'E228-400.5Q6.5<u', 'E228-400.5Q6.6<u', 'E228-400.6Q7.0<u', 'E228-400.6Q7.1<u', 'E228-400.6Q7.2<u',
    'E228-400.6Q7.3<u', 'E228-400.6Q7.4<u', 'E228-400.6Q7.5<u', 'E228-400.7Q8.0<u', 'E228-400.7Q8.1<u',
    'E228-400.7Q8.2<u', 'E228-400.7Q8.3<u', 'E228-400.7Q8.4<u', 'E228-400.7Q8.5<u', 'E228-400.7Q8.6<u',
    'E228-400.7Q8.7<u', 'E228-400.7Q8.8<u', 'E228-400.7Q8.9<u', 'E228-400.8Q9.0<u', 'E228-400.8Q9.1<u',
    'E228-400.8Q9.2<u', 'E228-400.8Q9.3<u', 'E228-400.8Q9.4<u', 'E228-400.8Q9.5<u', 'E228-400.8Q9.6<u',
    'E228-400.8Q9.7<u', 'E228-400.8Q9.8<u', 'E228-400.8Q9.9<u', 'E228-400.8Q9.10<u', 'E228-400.11Q12.0>u',
    'E228-400.11Q12.1>u', 'E228-400.11Q12.2>u', 'E228-400.11Q12.3>u', 'E228-400.11Q12.4>u',
    'E228-400.11Q12.5>u', 'E228-400.11Q12.6>u', 'E228-400.11Q12.7>u', 'E228-400.11Q12.8>u',
    'E228-400.11Q12.9>u', 'E228-400.11Q12.10>u', 'E228-400.11Q12.11>u', 'E228-400.11Q12.12>u',
    'E228-400.11Q12.13>u', 'E228-400.11Q12.14>u', 'E228-400.12Q13.0>u', 'E228-400.12Q13.1>u',
    'E228-400.12Q13.2>u', 'E228-400.12Q13.3>u', 'E228-400.12Q13.4>u', 'E228-400.12Q13.5>u',
    'E228-400.12Q13.6>u', 'E228-400.12Q13.7>u', 'E228-400.12Q13.8>u', 'E228-400.12Q13.9>u',
    'E228-400.12Q13.10>u', 'E228-400.12Q13.11>u', 'E228-400.12Q13.13>u', 'E228-400.12Q13.14>u',
    'E228-400.12Q13.15>u', 'E228-400.13Q14.0>u', 'E228-400.13Q14.1>u', 'E228-400.13Q14.2>u',
    'E228-400.13Q14.3>u', 'E228-400.13Q14.4>u', 'E228-400.13Q14.5>u', 'E228-400.13Q14.6>u',
    'E228-400.13Q14.7>u', 'E228-400.13Q14.8>u', 'E228-400.13Q14.9>u', 'E228-400.13Q14.11>u',
    'E228-400.13Q14.12>u', 'E772.11Q12.0', 'E772.11Q12.1', 'E772.11Q12.2', 'E772.11Q12.3',
    'E772.11Q12.4', 'E772.11Q12.5', 'E772.11Q12.6', 'E772.11Q12.7', 'E772.11Q12.8', 'E772.11Q12.9',
    'E772.11Q12.10', 'E772.12Q13.0', 'E772.12Q13.1', 'E772.12Q13.2', 'E772.12Q13.3', 'E772.12Q13.4',
    'E772.12Q13.5', 'E772.12Q13.7', 'E772.13Q14.0', 'E772.13Q14.1', 'E772.13Q14.2', 'E772.13Q14.3',
    'E772.13Q14.6', 'E772.14Q15.1', 'E605.7Q8.0<u', 'E605.7Q8.1<u', 'E605.7Q8.2<u', 'E605.7Q8.3<u',
    'E605.7Q8.4<u', 'E605.7Q8.5<u', 'E605.7Q8.6<u', 'E605.7Q8.7<u', 'E605.7Q8.8<u', 'E605.7Q8.9<u',
    'E605.8Q9.0<u', 'E605.8Q9.1<u', 'E605.8Q9.2<u', 'E605.8Q9.3<u', 'E605.8Q9.4<u', 'E605.8Q9.5<u',
    'E605.8Q9.6<u', 'E605.8Q9.7<u', 'E605.8Q9.8<u', 'E605.8Q9.9<u', 'E605.8Q9.10<u', 'E605.11Q13.0>u',
    'E605.11Q13.1>u', 'E605.11Q13.2>u', 'E605.11Q13.3>u', 'E605.11Q13.4>u', 'E605.11Q13.5>u',
    'E605.11Q13.6>u', 'E605.11Q13.7>u', 'E605.11Q13.8>u', 'E605.11Q13.9>u', 'E605.11Q13.10>u',
    'E605.11Q13.11>u', 'E605.11Q13.12>u', 'E605.11Q13.13>u', 'E605.11Q13.14>u', 'E605.11Q13.15>u',
    'E605.11Q13.16>u', 'E605.13Q18.0>u', 'E605.13Q18.1>u', 'E605.13Q18.2>u', 'E605.13Q18.3>u',
    'E605.13Q18.4>u', 'E605.13Q18.5>u', 'E605.13Q18.6>u', 'E605.13Q18.7>u', 'E605.13Q18.8>u',
    'E605.13Q18.9>u', 'E605.13Q18.10>u', 'E605.13Q18.11>u', 'E605.13Q18.12>u', 'E605.13Q18.13>u', 'E605.13Q18.14>u']


def cutFuncSV19(p):
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
            
        if(p["id"][-2:]=="<u" and p["<Q>"]>10.5):
            return False,p
        if(p["id"][-2:]==">u" and p["<Q>"]<10.5):
            return False,p
    
#    return delta<0.5 and p.qT_avarage<80
    return (delta<0.1 or (delta<0.25 and par/err*delta**2<1)) , p

### check the point against the list of dropping points.
def cutFuncSV21(p):    
    
    #### check against the presence in the reduced set.
    if(not(p["id"] in includePoints)):
        return False,p
    
    return cutFuncSV19(p)

def cutFuncPlot(p):
    if p["type"]=="DY":        
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
        
        if(p["id"][-2:]=="<u" and p["<Q>"]>10.5):
            return False,p
        if(p["id"][-2:]==">u" and p["<Q>"]<10.5):
            return False,p
    
#    return delta<0.5 and p.qT_avarage<80
    return delta<0.25 , p

#%%
#######################################
# Loading the data set
#######################################

setHE=loadThisData(['CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                      'A7-00y10','A7-10y20','A7-20y24', 
                      'A8-00y04','A8-04y08','A8-08y12','A8-12y16',
                      'A8-16y20','A8-20y24','A8-46Q66','A8-116Q150',
                      'CMS7', 'CMS8', 
                      'CMS13-00y04','CMS13-04y08','CMS13-08y12','CMS13-12y16','CMS13-16y24',
                      'LHCb7', 'LHCb8', 'LHCb13'])

#### I separate data above and below UPSILON, I create two copies of LE data with different names
#### the data to be split only  'E228-300', 'E228-400' and E605
setLE1=loadThisData(['PHE200', 'E228-200','E772'])
setLE2=loadThisData(['E228-300', 'E228-400','E605'])
setLE3=loadThisData(['E228-300', 'E228-400','E605'])
for s in setLE2:
    s.name+="-blwUPS"
    for p in s.points:
        p["id"]+="<u"
for s in setLE3:
    s.name+="-abvUPS"
    for p in s.points:
        p["id"]+=">u"
setLE=[setLE1[0],setLE1[1],setLE2[0],setLE3[0],setLE2[1],setLE3[1],setLE1[2],setLE2[2],setLE3[2]]


theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",setHE+setLE)

setSV21=theData.CutData(cutFuncSV21) 

setSV19=theData.CutData(cutFuncSV19) 

setPLOT=theData.CutData(cutFuncPlot)

print('SV21: Loaded ', setSV21.numberOfSets, 'data sets with ', sum([i.numberOfPoints for i in setSV21.sets]), 'points.')

print('SV19: Loaded ', setSV19.numberOfSets, 'data sets with ', sum([i.numberOfPoints for i in setSV19.sets]), 'points.')

print('PLOT: Loaded ', setPLOT.numberOfSets, 'data sets with ', sum([i.numberOfPoints for i in setPLOT.sets]), 'points.')

#%%
harpy.setNPparameters([2., 0.0276492, 0.271802, 2.38402, 0.352648, 1.73888, 0.0971883, \
17.6768, 0.100485, 11.6811, 0.278435, 3.35306, 149.2, 0])
#rSet.SetReplica()

DataProcessor.harpyInterface.PrintChi2Table(setSV21,printDecomposedChi2=True)
DataProcessor.harpyInterface.PrintChi2Table(setSV19,printDecomposedChi2=True)

#%%
# #######################################################
# ### Makes a file header (USE IT ONLY IF NEED TO REWRITE FILE)
# #######################################################
# #PDFinUse="HERA20"
# #PDFinUse="NNPDF31"
# #PDFinUse="CT18"
# #PDFinUse="MSHT20"

# with open('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV21/'+PDFinUse+'/chi2-'+PDFinUse+'-EXP-SV19.dat', 'w') as outfile:
#         outfile.write(str(['repNum', 'Total', [s.name for s in setSV19.sets]])+"\n")
#         outfile.write(str([0, setSV19.numberOfPoints, [s.numberOfPoints for s in setSV19.sets]])+"\n")
# with open('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV21/'+PDFinUse+'/chi2-'+PDFinUse+'-EXP-SV21.dat', 'w') as outfile:
#         outfile.write(str(['repNum', 'Total', [s.name for s in setSV21.sets]])+"\n")
#         outfile.write(str([0, setSV21.numberOfPoints, [s.numberOfPoints for s in setSV21.sets]])+"\n")

#%%
# ##################################
# ## Save replicas of NNPDFs
# ##################################

import pickle
import DataProcessor.SaveTMDGrid

rSet.SetReplica()

for i in range(rSet.numberOfReplicas):      
    rSet.SetReplica(i)
    
    #### save list of chi2 for full set
    chiList=DataProcessor.harpyInterface.ComputeChi2(setSV19)    
    with open('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV21/'+PDFinUse+'/chi2-'+PDFinUse+'-EXP-SV19.dat', 'a') as outfile:
        outfile.write(str([i]+list(chiList))+"\n")
    
    #### save list of chi2 for reduced set
    chiList=DataProcessor.harpyInterface.ComputeChi2(setSV21)    
    with open('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV21/'+PDFinUse+'/chi2-'+PDFinUse+'-EXP-SV21.dat', 'a') as outfile:
        outfile.write(str([i]+list(chiList))+"\n")
    
    #### save cross-section to plot
    listToSave=DataProcessor.harpyInterface.ComputeXSec(setPLOT)
    path='/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV21/'+PDFinUse+'/xSec-EXP/'\
        +'{:04d}'.format(i)+'.pick'
        
    with open(path, "wb") as filehandle:
        pickle.dump(listToSave,filehandle)
        
    #### save TMD-distribution        
    path='/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD-grids/Grids_optimal/SV21-'+PDFinUse+'_nnlo_EXP/SV21-'+PDFinUse+'_nnlo_EXP'\
        +'_'+'{:04d}'.format(i)+'.dat'
    DataProcessor.SaveTMDGrid.SaveGrid_optimal(path)

#%%
# # ######################################
# # ## Posteriory routine, which collects all pickle files into a single CSV using gaussian statistics
# # ######################################
# import glob
# import pickle
# import numpy

# #PDFinUse="HERA20"
# #PDFinUse="NNPDF31"
# #PDFinUse="CT18"
# PDFinUse="MSHT20"

# fileList=glob.glob('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV21/'+PDFinUse+'/xSec-EXP/*.pick')
# centralPath='/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV21/'+PDFinUse+'/xSec-EXP/0000.pick'

# with open(centralPath, "rb") as filehandle:
#     central=pickle.load(filehandle)

# reps=[]
# for pp in fileList:    
#     with open(pp, "rb") as filehandle:
#         reps.append(pickle.load(filehandle))

# central=numpy.mean(reps,axis=0)
# std=numpy.std(reps,axis=0)

# pointNames=[s["id"] for s in setPLOT.points]

# #with open('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV19-noFIT/xSec-SV19rep.dat','w') as file:
# with open('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV21/'+PDFinUse+'/xSec-'+PDFinUse+'-EXP.dat','w') as file:
#     file.write("Point id, xSec, Std \n")
#     for i in range(len(pointNames)):
#         file.write(pointNames[i]+', '+'{:g}'.format(central[i])+', '+'{:g}'.format(std[i])+" \n")

#%%
# # ######################################
# # ## Posteriory routine, which collects all pickle files into a single CSV computing 68%CI
# # ######################################
import glob
import pickle
import numpy

#########################################################
## Resample data with N/2
#########################################################
def Resample(dd):
    ss=int(numpy.floor(len(dd)/2))
    ss=max(100,ss)
    return dd[numpy.random.choice(dd.shape[0], size=ss)]

#########################################################
## Resample data for (d1+d2)/2 with N/2
#########################################################
def ResampleSUM(dd1,dd2):    
    ss=min( max(100,int(numpy.floor(len(dd1)/2))) , max(100,int(numpy.floor(len(dd2)/2))))
    return (dd1[numpy.random.choice(dd1.shape[0], size=ss)]+dd2[numpy.random.choice(dd2.shape[0], size=ss)])/2

#########################################################
## Determine Mean, 68%CI by resampling 
#########################################################
alpha=68

def ComputeParametersSUM(dd1,dd2):
    means=[]
    lowers=[]
    uppers=[]    
    for i in range(1500):
        sample=ResampleSUM(numpy.array(dd1),numpy.array(dd2))
        means.append(numpy.mean(sample))
        lowers.append(numpy.percentile(sample,(100-alpha)/2))
        uppers.append(numpy.percentile(sample,100-(100-alpha)/2))
    
    return [numpy.mean(means),numpy.mean(lowers),numpy.mean(uppers)]

def ComputeParameters(dd):    
    means=[]
    lowers=[]
    uppers=[]    
    for i in range(1500):
        sample=Resample(numpy.array(dd))
        means.append(numpy.mean(sample))
        lowers.append(numpy.percentile(sample,(100-alpha)/2))
        uppers.append(numpy.percentile(sample,100-(100-alpha)/2))
    
    return [numpy.mean(means),numpy.mean(lowers),numpy.mean(uppers)]
#%%
#PDFinUse="HERA20"
#PDFinUse="NNPDF31"
#PDFinUse="CT18"
PDFinUse="MSHT20"

fileList=glob.glob('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV21/'+PDFinUse+'/xSec/*.pick')
repsPDF=[]
for pp in fileList:    
    with open(pp, "rb") as filehandle:
        repsPDF.append(pickle.load(filehandle))
        
fileList=glob.glob('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV21/'+PDFinUse+'/xSec-EXP/*.pick')
repsEXP=[]
for pp in fileList:    
    with open(pp, "rb") as filehandle:
        repsEXP.append(pickle.load(filehandle))


pointNames=[s["id"] for s in setPLOT.points]

with open('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV21/'+PDFinUse+'/xSec-'+PDFinUse+'-ALL.dat','w') as file:
    file.write("Point id, xSec-PDF, std-PDF, 68CI-down-PDF, 68CI-up-PDF, \
               xSec-EXP, std-EXP, 68CI-down-EXP, 68CI-up-EXP, \
               xSec-ALL, std-ALL, 68CI-down-ALL, 68CI-up-ALL  \n")
    for i in range(len(pointNames)):
        print(i,'/',len(pointNames))
        d1=[r[i] for r in repsPDF]
        d2=[r[i] for r in repsEXP]
        
        centralPDF=numpy.mean(d1)
        stdPDF=numpy.std(d1)
        
        centralEXP=numpy.mean(d2)
        stdEXP=numpy.std(d2)
        
        centralALL=(centralPDF+centralEXP)/2
        stdALL=numpy.sqrt(stdEXP**2+stdPDF**2)
        
        ppPDF=ComputeParameters(d1)
        ppEXP=ComputeParameters(d2)
        ppALL=ComputeParametersSUM(d1,d2)
        
        
        file.write(pointNames[i]+', '\
                +'{:g}'.format(centralPDF)+', '+'{:g}'.format(stdPDF)+', '+'{:g}'.format(ppPDF[1])+', '+'{:g}'.format(ppPDF[2])\
                +'{:g}'.format(centralEXP)+', '+'{:g}'.format(stdEXP)+', '+'{:g}'.format(ppEXP[1])+', '+'{:g}'.format(ppEXP[2])\
                +'{:g}'.format(centralALL)+', '+'{:g}'.format(stdALL)+', '+'{:g}'.format(ppALL[1])+', '+'{:g}'.format(ppALL[2])+" \n")
