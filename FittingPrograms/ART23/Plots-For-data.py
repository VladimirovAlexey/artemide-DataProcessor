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

replicaFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/replicaDY.dat"
logFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/LOG.log"

import sys
import numpy
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)

#%%
import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet

MAINPATH=ROOT_DIR 

#%%
#######################################
#Initialize artemide
#######################################
import harpy

PDFtoUSE="MSHT"
#PDFtoUSE="NNPDF"

path_to_constants=MAINPATH+"FittingPrograms/ART23/ConstantsFiles/DYonly/ART23_"+PDFtoUSE+"_N4LL.atmde"
#path_to_constants=MAINPATH+"FittingPrograms/ART23/ConstantsFiles/DYonly/ART23_JAM_N4LL"

harpy.initialize(path_to_constants)
#%%
if(PDFtoUSE=="MSHT"):
    rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/data/WorkingFiles/TMD/Fit_Notes/ART23/REPS/ART23_run2.rep")
elif(PDFtoUSE=="NNPDF"):
    rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                            "/data/WorkingFiles/TMD/Fit_Notes/ART23/REPS/ART23_NNPDF.rep")
else:
    raise Exception("no no no")
rSet.SetReplica(-1)

rSet.SetReplica(0)
#%%
### read the list of files and return the list of DataSets
def loadThisDataDY(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=ROOT_DIR+"DataLib/unpolDY/"
    path_to_dataW=ROOT_DIR+"DataLib/unpolW/"
    
    
    dataCollection=[]
    for name in listOfNames:
        if(name[-1]=="W"):
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataW+name+".csv")
        else:
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        

    return dataCollection
#%%
##################Cut function
def cutFunc(p):
    par=0.5
    
    #  for artemide v3.    
    # p["process"]=[p["process"][0],p["process"][2],1,1]
    if(True):
        if(p["process"][2]==1): p["process"]=[p["process"][0],1,1,1]
        elif(p["process"][2]==2): p["process"]=[p["process"][0],1,1,-1]
        elif(p["process"][2]==3): p["process"]=[p["process"][0],2,1,1]
        elif(p["process"][2]==4): p["process"]=[p["process"][0],2,1,-1]
        elif(p["process"][2]==5): p["process"]=[p["process"][0],3,1,1]
        elif(p["process"][2]==6): p["process"]=[p["process"][0],3,1,-1]
        elif(p["process"][2]==7): p["process"]=[p["process"][0],4,1,1]
        elif(p["process"][2]==8): p["process"]=[p["process"][0],5,1,1]
        elif(p["process"][2]==9): p["process"]=[p["process"][0],6,1,1]
        elif(p["process"][2]==10): p["process"]=[p["process"][0],4,1,-1]
        elif(p["process"][2]==11): p["process"]=[p["process"][0],5,1,-1]
        elif(p["process"][2]==12): p["process"]=[p["process"][0],6,1,-1]
        elif(p["process"][2]==1001): p["process"]=[p["process"][0],101,1,1]
        elif(p["process"][2]==1002): p["process"]=[p["process"][0],102,1,1]
        else:
            print("UNKNOWN PROCESS IN ARTEMIDE 3"+str(p["process"]))
        
        
    
    if(p["xSec"]>0):
        err=numpy.sqrt(sum([i**2 for i in p["uncorrErr"]]))/p["xSec"]
    else:
        err=100.
    delta=p["<qT>"]/p["<Q>"]
    
    if(p["id"][0] == "E"):
        delta=p["<qT>"]/p["Q"][1] 
    
    if("run1-W" in p["id"]):
        delta=p["qT"][0]/(p["Q"][0]+5.)
    
    
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
    return ((delta<0.25 and p["<qT>"]<10.) or (delta<0.25 and par/err*delta**2<1)) , p

def cutFuncPlot(p):
    if p["type"]=="DY":        
        #  for artemide v3.
        #p["process"]=[p["process"][0],p["process"][2],1,1]
        if(True):
            if(p["process"][2]==1): p["process"]=[p["process"][0],1,1,1]
            elif(p["process"][2]==2): p["process"]=[p["process"][0],1,1,-1]
            elif(p["process"][2]==3): p["process"]=[p["process"][0],2,1,1]
            elif(p["process"][2]==4): p["process"]=[p["process"][0],2,1,-1]
            elif(p["process"][2]==5): p["process"]=[p["process"][0],3,1,1]
            elif(p["process"][2]==6): p["process"]=[p["process"][0],3,1,-1]
            elif(p["process"][2]==7): p["process"]=[p["process"][0],4,1,1]
            elif(p["process"][2]==8): p["process"]=[p["process"][0],5,1,1]
            elif(p["process"][2]==9): p["process"]=[p["process"][0],6,1,1]
            elif(p["process"][2]==10): p["process"]=[p["process"][0],4,1,-1]
            elif(p["process"][2]==11): p["process"]=[p["process"][0],5,1,-1]
            elif(p["process"][2]==12): p["process"]=[p["process"][0],6,1,-1]
            elif(p["process"][2]==1001): p["process"]=[p["process"][0],101,1,1]
            elif(p["process"][2]==1002): p["process"]=[p["process"][0],102,1,1]
            else:
                print("UNKNOWN PROCESS IN ARTEMIDE 3"+str(p["process"]))
                
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
    return (delta<0.35 or p["<qT>"]<10.) , p

#%%
### Loading the DY data set
theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY([
                           'CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                           #'A7-00y10', 'A7-10y20','A7-20y24', 
                           'A8-00y04', 'A8-04y08', 'A8-08y12', 'A8-12y16', 'A8-16y20', 'A8-20y24', 
                           'A8-46Q66', 'A8-116Q150', 
                           'A13-norm',
                           'CMS7', 'CMS8', 
                           'CMS13-00y04','CMS13-04y08','CMS13-08y12','CMS13-12y16','CMS13-16y24',
                           #'CMS13_dQ_50to76',
                           'CMS13_dQ_106to170','CMS13_dQ_170to350','CMS13_dQ_350to1000',
                           'LHCb7', 'LHCb8', 'LHCb13_dy(2021)', 
                           'PHE200', 'STAR510', 
                           'E228-200', 'E228-300', 'E228-400', 
                           'E772',
                           'E605',
                           'D0run1-W','CDFrun1-W'
                          ]))

setDY=theData.CutData(cutFuncPlot) 

setDYFIT=theData.CutData(cutFunc) 

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
#%%

DataProcessor.harpyInterface.PrintChi2Table(setDYFIT,printDecomposedChi2=True)

#%%

A11=DataProcessor.harpyInterface.ComputeXSec(setDYFIT)
A22=[p["id"] for p in setDYFIT.points]


#%%
DIR_PLOT='/data/WorkingFiles/TMD/Fit_Notes/ART23/Data_for_plot_'+PDFtoUSE+'/'
#DIR_PLOT='/data/WorkingFiles/TMD/Fit_Notes/ART23/Data_for_plot_JLAB/'
#%%
# #######################################################
# ### Makes a file header (USE IT ONLY IF NEED TO REWRITE FILE)
# #######################################################
# with open(DIR_PLOT+'header.dat', 'w') as outfile:
#         outfile.write(str(['repNum', 'Total', [s.name for s in setDY.sets]])+"\n")
#         outfile.write(str([0, setDY.numberOfPoints, [s.numberOfPoints for s in setDY.sets]])+"\n")

#%%
# # ##################################
# # ## Save central replica
# # ##################################

# import pickle

# rSet.SetReplica()
   
    
# #### save cross-section to plot
# listToSave=DataProcessor.harpyInterface.ComputeXSec(setDY)
# path=DIR_PLOT+'/{:04d}'.format(0)+'.pick'
    
# with open(path, "wb") as filehandle:
#     pickle.dump(listToSave,filehandle)
    
#%%
# # ##################################
# # ## Save replicas of xSec
# # ##################################

import pickle

rSet.SetReplica()

nn=2

nmin=nn*100+1
nmax=nn*100+1+100
if(nmin>=rSet.numberOfReplicas):
    sys.exit()
if(nmax>=rSet.numberOfReplicas): nmax=rSet.numberOfReplicas


for i in range(nmin,nmax):      
    rSet.SetReplica(i)    
    
# for i in range(100):
#     harpy.setPDFreplica(i)
#     harpy.setNPparameters([1.56142,0.0369174,0.0581734,1.0,
#                            0.874245,0.913883,0.991563,6.05412,
#                            0.353908,46.6064,0.115161,1.53235,
#                            1.31966,0.434833,0.0,0.0])
    
    
    #### save cross-section to plot
    listToSave=DataProcessor.harpyInterface.ComputeXSec(setDY)
    path=DIR_PLOT+'{:04d}'.format(i)+'.pick'
        
    with open(path, "wb") as filehandle:
        pickle.dump(listToSave,filehandle)

sys.exit()

#%%
# ######################################
# ## Posteriory routine, which collects all pickle files into a single CSV using bootstrap 68%CI
# ## Joined result is by mean=weighted mean, 68%CI->joined
# ######################################

import glob
import pickle
import numpy


#### compute the 68CI for list of arrays, axis=0
def Bootstrap68CI(dd):
    rMin=[]
    rMax=[]
    
    for i in range(2000):
        indices=numpy.random.choice(range(len(dd)),size=int(len(dd)/2))
        sample=[dd[i] for i in indices]
        rMin.append(numpy.quantile(sample, (1-.68)/2,axis=0))
        rMax.append(numpy.quantile(sample, 1-(1-.68)/2,axis=0))
        
    return (numpy.array([numpy.mean(rMin,axis=0),numpy.mean(rMax,axis=0)])).transpose()
        


fileList=glob.glob('/data/WorkingFiles/TMD/Fit_Notes/ART23/Data_for_plot_'+PDFtoUSE+'/*.pick')
#fileList=glob.glob('/data/WorkingFiles/TMD/Fit_Notes/ART23/Data_for_plot_JLAB/*.pick')

### reading all replicas
reps=[]
for pp in fileList:    
    with open(pp, "rb") as filehandle:
        reps.append(pickle.load(filehandle))

### computing central value and the band
        
## just mean
mean=numpy.mean(reps,axis=0)

## 68 CI's
CI68=Bootstrap68CI(reps)

centralReplica=reps[0]

#%%
### determine if the point is in the fited region 
### (note it includes non-fit experiments, but under fit cuts)
point_in_FIT=[]
for pp in setDY.points:
    point_in_FIT.append(pp in setDYFIT.points)

### For each experiment, determine the chi2 computed from the mean result (total,chi2D,chiL)
### For each experiment, determine the CI68% for chi2 computed from replicas
### For each experiment, determine the shift due to the systematics
### fro included points the shift is made by DetermineShift code
### fro extra point, the resuls is continued by the last point
### (one could not use the avarage shift since extra points deviate from the fac.theorm and thus break the shift)
EXP_name=[]
EXP_points=[]
EXP_chi2=[]
EXP_chi2D=[]
EXP_chi2L=[]
EXP_chi2_mean=[]
EXP_chi2_68CI=[]
EXP_chi2_central=[]
shifts=[]
shiftsCentral=[]
avshift=[]
avshiftCentral=[]
for s in range(len(setDYFIT.sets)):
        
    s1=setDYFIT.sets[s]
    sPLOT=setDY.sets[s]
    EXP_name.append(s1.name)
    EXP_points.append(s1.numberOfPoints)
    # extracting points from the mean value (
    # in order of appearence because it may happen that the set has lacuns)
    xx=[]
    xx0=[]
    p1=s1.points
    for pp in p1:
        xx.append(mean[setDY.points.index(pp)])
        xx0.append(centralReplica[setDY.points.index(pp)])
    
    # computing chi2 for mean line
    EXP_chi2.append(s1.chi2(xx))
    # computing chi2 for centralReplica
    EXP_chi2_central.append(s1.chi2(xx0))   
    
    # computing chi2 decomposition for mean line
    chiDecomposed=s1.DecomposeChi2(xx)
    EXP_chi2D.append(chiDecomposed[0])
    EXP_chi2L.append(chiDecomposed[1])
    
    
    # for each replica compute chi2, mean and 68CI
    chi2_reps=[]
    for rr in range(len(reps)):
        xxR=[]
        for pp in p1:
            xxR.append(reps[rr][setDY.points.index(pp)])
        chi2_reps.append(s1.chi2(xxR))
    
    EXP_chi2_mean.append(numpy.mean(chi2_reps))
    EXP_chi2_68CI.append(Bootstrap68CI(chi2_reps))
    
    # computing the shift for the mean line
    shift=s1.DetermineSystematicShift(xx)
    avshift.append(s1.DetermineAvarageSystematicShift(xx))
    # computing the shift for the central replica
    shiftCentral=s1.DetermineSystematicShift(xx0)
    avshiftCentral.append(s1.DetermineAvarageSystematicShift(xx0))
    
    # now I compute the shift for the plot data
    # If point is in fit area I use the true computation value
    # if the point is beyond it, I use the % shift in the last point
    k=0
    lastshift=0.
    for pp in setDY.sets[s].points:
        index_pp=setDY.points.index(pp)
        if(pp in s1.points):
            shifts.append(shift[k])
            lastshift=shift[k]/xx[k]
            k+=1            
        else:
            shifts.append(mean[index_pp]*lastshift)
            
    # now I compute the shift for the plot data (with central replica)
    # If point is in fit area I use the true computation value
    # if the point is beyond it, I use the % shift in the last point
    k=0
    lastshift=0.
    for pp in setDY.sets[s].points:
        index_pp=setDY.points.index(pp)
        if(pp in s1.points):
            shiftsCentral.append(shift[k])
            lastshift=shiftCentral[k]/xx[k]
            k+=1            
        else:
            shiftsCentral.append(mean[index_pp]*lastshift)
    

#%%
pointNames=[s["id"] for s in setDY.points]

if(PDFtoUSE=="MSHT"): ttt="ART23"
elif(PDFtoUSE=="NNPDF"): ttt="NNPDF"
else: raise Exception("non nonon")
#ttt="JAM"

with open('/data/WorkingFiles/TMD/Fit_Notes/ART23/DataPlots/xSec-'+ttt+'.dat','w') as file:
    file.write("Point id, xSec(mean), xSec(central rep.), 68min, 68max, shift(mean), shift(central rep.), FIT \n")
    for i in range(len(pointNames)):
        file.write(pointNames[i]+', '
                    +'{:g}'.format(mean[i])+', '
                    +'{:g}'.format(centralReplica[i])+', '
                    +'{:g}'.format(CI68[i][0])+', '
                    +'{:g}'.format(CI68[i][1])+', '
                    +'{:g}'.format(shifts[i])+', '
                    +'{:g}'.format(shiftsCentral[i])+', '
                    +'{:g}'.format(point_in_FIT[i])
                    +" \n")
        
with open('/data/WorkingFiles/TMD/Fit_Notes/ART23/DataPlots/xSec-'+ttt+'_full.dat','w') as file:
    file.write("Point id, xSec, xSec(central rep.), 68min, 68max, shift, shift(central rep.), FIT, xSec(experimental), unc(experimental), <qT>, qtMIN, qtMAX \n")
    for i in range(len(pointNames)):
        file.write(pointNames[i]+', '
                    +'{:g}'.format(mean[i])+', '
                    +'{:g}'.format(centralReplica[i])+', '
                    +'{:g}'.format(CI68[i][0])+', '
                    +'{:g}'.format(CI68[i][1])+', '
                    +'{:g}'.format(shifts[i])+', '
                    +'{:g}'.format(shiftsCentral[i])+', '                   
                    +'{:g}'.format(point_in_FIT[i])+', '                    
                    +'{:g}'.format(setDY.points[i]["xSec"])+', '    
                    +'{:g}'.format(numpy.sqrt(numpy.sum(numpy.array(setDY.points[i]["uncorrErr"])**2)))+', '    
                    +'{:g}'.format(setDY.points[i]["<qT>"])+', '    
                    +'{:g}'.format(setDY.points[i]["qT"][0])+', '    
                    +'{:g}'.format(setDY.points[i]["qT"][1])
                    +" \n")


with open('/data/WorkingFiles/TMD/Fit_Notes/ART23/DataPlots/chi2-'+ttt+'.dat','w') as file:
    file.write("Set id, Npt, chi2(mean xSec), chi2(central replica), chi2(mean), chi2(68% low), chi2(68% up), <shift>, chi2D(mean xSec), chi2L(mean xSec), \n")
    for i in range(len(setDYFIT.sets)):
        file.write(setDYFIT.sets[i].name+', '
                    +'{:g}'.format(setDYFIT.sets[i].numberOfPoints)+', '
                    +'{:g}'.format(EXP_chi2[i])+', '
                    +'{:g}'.format(EXP_chi2_central[i])+', '
                    +'{:g}'.format(EXP_chi2_mean[i])+', '                    
                    +'{:g}'.format(EXP_chi2_68CI[i][0])+', '                    
                    +'{:g}'.format(EXP_chi2_68CI[i][1])+', '                    
                    +'{:g}'.format(avshift[i])+', '
                    +'{:g}'.format(EXP_chi2D[i])+', '
                    +'{:g}'.format(EXP_chi2L[i])
                    +" \n")
        
# #%%
# with open('/data/WorkingFiles/TMD/Fit_Notes/ART23/DataPlots/xSec-DATA.dat','w') as file:
#     for i in range(len(pointNames)):
#         file.write(pointNames[i]+', ')
#     file.write('\n')
#     for j in range(len(reps)):
#         for i in range(len(pointNames)):
#             file.write('{:g}'.format(reps[j][i])+', ')
#         file.write('\n')