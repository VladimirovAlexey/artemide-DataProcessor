#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 14:35:35 2021

@author: vla18041
"""

##############################
# Ploting original SV19 fit
##############################

#######################################
# importing libraries
#######################################

import sys
import numpy
sys.path.append("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor")
import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet
import DataProcessor.DataMultiSet

MAINPATH="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"
#%%
#######################################
#Initialize artemide with a replica -2
#######################################
import harpy
path_to_constants=MAINPATH+"FittingPrograms/SV19/Constants-files/"
harpy.initialize(path_to_constants+"DY_n3lo/const-NNPDF31_n3lo")

harpy.setNPparameters([2.0340, 0.0299, 0.2512, 7.7572,334.6108, 2.4543,-4.8203, 0.1000,  0.0000])


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
    
#    return delta<0.5 and p.qT_avarage<80
    return (delta<0.1 or (delta<0.25 and par/err*delta**2<1)) , p

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
    
#    return delta<0.5 and p.qT_avarage<80
    return delta<0.25 , p

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
setDYplot=theData.CutData(cutFuncPlot) 

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setDY.sets])

#%%

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/Replicas/"+
                                                  "DY_n3lo/DY_NNPDF31_n3lo.rep")
                                                  # "Sivers20_model9case1(noDY-n3lo).rep")

rSet.SetReplica()


#%%
# #######################################################
# ### Makes a file header (USE IT ONLY IF NEED TO REWRITE FILE)
# #######################################################
# with open('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV19-noFIT/chi2-NNPDFrep-fit.dat', 'w') as outfile:
#         outfile.write(str(['repNum', 'Total', [s.name for s in setDY.sets]])+"\n")
#         outfile.write(str([0, setDY.numberOfPoints, [s.numberOfPoints for s in setDY.sets]])+"\n")

#%%
# # ##################################
# # ## Save replicas of NNPDFs
# # ##################################

# r=21

# import pickle

# rSet.SetReplica()

# for i in range(r*50,(r+1)*50):  
#     harpy.setPDFreplica(i)
#     rSet.SetReplica()
    
#     chiList=DataProcessor.harpyInterface.ComputeChi2(setDY)    
#     with open('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV19-noFIT/chi2-NNPDFrep.dat', 'a') as outfile:
#         outfile.write(str([i]+list(chiList))+"\n")
    
    
#     listToSave=DataProcessor.harpyInterface.ComputeXSec(setDYplot)
#     path='/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV19-noFIT/xSec-NNPDFrep/'\
#         +'{:04d}'.format(i)+'.pick'
        
#     with open(path, "wb") as filehandle:
#         pickle.dump(listToSave,filehandle)

#%%
##################################
## Save replicas of SV19
##################################
# import pickle

# rSet.SetReplica()



# for i in range(rSet.numberOfReplicas):    
#     rSet.SetReplica(i)
    
#     chiList=DataProcessor.harpyInterface.ComputeChi2(setDY)    
#     with open('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV19-noFIT/chi2-SV19rep.dat', 'a') as outfile:
#         outfile.write(str([i]+list(chiList))+"\n")
    
    
#     listToSave=DataProcessor.harpyInterface.ComputeXSec(setDYplot)
#     path='/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV19-noFIT/xSec-SV19/'\
#         +'{:04d}'.format(i)+'.pick'
        
#     with open(path, "wb") as filehandle:
#         pickle.dump(listToSave,filehandle)

#%%
##################################
## Special case: NNPDF replicas each fitted by SV19
##################################
# import pickle
# import DataProcessor.SaveTMDGrid

# r=3

# for n in range(r*100,(r+1)*100):
#     #### maximum 920
#     if(n>920):
#         break
    
#     #### read n'th line from the file
#     with open('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/LOGS/NNPDF31/allPDF_centralEXP_full.dat','r') as file:
#         for i, line in enumerate(file):
#             if i == n:
#                 ll=eval(line)    
#             elif i>n:
#                 break
#     #### last two numbers [...., number of PDF-replica,  NP-parameters ]
#     rep=ll[3]
#     param=ll[4]
    
#     #### set this parameters
#     harpy.setPDFreplica(rep)
#     harpy.setNPparameters(param)
    
#     #### 1. Save chi^2 list
#     chiList=DataProcessor.harpyInterface.ComputeChi2(setDY)    
#     with open('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV19-noFIT/chi2-NNPDFrep-fit.dat', 'a') as outfile:
#         outfile.write(str([i]+list(chiList))+"\n")
    
#     #### 2. Save values of xSec
#     listToSave=DataProcessor.harpyInterface.ComputeXSec(setDYplot)
#     path='/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV19-noFIT/xSec-NNPDFrep-fit/'\
#         +'{:04d}'.format(i)+'.pick'
        
#     with open(path, "wb") as filehandle:
#         pickle.dump(listToSave,filehandle)
    
#     #### 3. Save grid        
#     path='/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD-grids/Grids_optimal/SV19_n3lo_PDFrep-fit/SV19_n3lo_PDFrep-fit'\
#     +'_'+'{:04d}'.format(rep)+'.dat'
#     DataProcessor.SaveTMDGrid.SaveGrid_optimal(path)

#%%
# ######################################
# ## Posteriory routine, which collects all pickle files into a single CSV
# ######################################
import glob
import pickle
import numpy

# fileList=glob.glob('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV19-noFIT/xSec-SV19/*.pick')
# centralPath='/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV19-noFIT/xSec-SV19/0000.pick'

fileList=glob.glob('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV19-noFIT/xSec-NNPDFrep/*.pick')
centralPath='/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV19-noFIT/xSec-NNPDFrep/0000.pick'

with open(centralPath, "rb") as filehandle:
    central=pickle.load(filehandle)

reps=[]
for pp in fileList:    
    with open(pp, "rb") as filehandle:
        reps.append(pickle.load(filehandle))

std=numpy.std(reps,axis=0)

pointNames=[s["id"] for s in setDYplot.points]

#with open('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV19-noFIT/xSec-SV19rep.dat','w') as file:
with open('/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD<->PDF/Figures/SV19-noFIT/xSec-NNPDFrep.dat','w') as file:
    file.write("Point id, xSec, Std \n")
    for i in range(len(pointNames)):
        file.write(pointNames[i]+', '+'{:g}'.format(central[i])+', '+'{:g}'.format(std[i])+" \n")