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
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"

ATMDE_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide/"
HARPY_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide/harpy/"


import sys
import numpy
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)


#%%
import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet

MAINPATH=ROOT_DIR 
#%%
#######################################
#Initialize artemide
#######################################
import harpy

import DataProcessor.ArtemideReplicaSet


path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/ART25_main.atmde"
harpy.initialize(path_to_constants)

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(\
                        "/data/WorkingFiles/TMD/Fit_Notes/ART25/REPLICAS/ART25_run1.rep")
    
rSet.SetReplica(0)

#%%
### read the list of files and return the list of DataSets
def loadThisDataSIDIS(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=ROOT_DIR+"DataLib/unpolSIDIS/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        
    return dataCollection

#%%
### read the list of files and return the list of DataSets
def loadThisDataDY(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=ROOT_DIR+"DataLib/unpolDY/"
    path_to_dataW=ROOT_DIR+"DataLib/unpolW/"
    path_to_dataA=ROOT_DIR+"DataLib/DY_angular/"
    
    
    dataCollection=[]
    for name in listOfNames:
        if(name[-1]=="W"):
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataW+name+".csv")
        elif("_A4" in name):
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataA+name+".csv")
        elif("_Auu" in name):
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataA+name+".csv")
        else:
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        

    return dataCollection

#%%
##################Cut function
def cutFunc(p):
    
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
        
        #### in the case of MAPFF I use special sets for pi- and K-
        if('MAPFF' in path_to_constants):            
            if(p["process"][3]==2101 and p["process"][2]<0): 
                p["process"][3]=2107
                p["process"][2]=3
            if(p["process"][3]==2103 and p["process"][2]<0): 
                p["process"][3]=2108
                p["process"][2]=3
            if(p["process"][3]==2105 and p["process"][2]<0): 
                p["process"][3]=2109
                p["process"][2]=3
            if(p["process"][2]==-1): p["process"][2]=3
            if(p["process"][2]==-2): p["process"][2]=4
            
        
        return (delta<0.1 or (delta<0.25 and par/err*delta**2<1)) , p
    elif p["type"]=="DY":
        par=0.5
        
        #  for artemide v3.    
        # p["process"]=[p["process"][0],p["process"][2],1,1]
        if(len(p["process"])==3):        
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
        
        if(p["id"][0:2] == "A7"):            
                return False , p
        elif('CMS13_dQ_50to76' in p["id"]):
                return False , p
        
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
        
        return ((delta<0.25 and p["<qT>"]<10.) or (delta<0.25 and par/err*delta**2<1)) , p
    
##################Cut function
def cutFunc_PLOT(p):
    
    if p["type"]=="SIDIS":  
        
        
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
        if qTWORST>p["<Q>"]*0.55:
            return False , p
        
        # factor1=-2.*(1.-numpy.sqrt(1+gamma2*(1-delta**2)))/gamma2
        
        # if(p["x"][0]*factor1<0.002):
        #     return False , p
    
        ### drop Q<2
        if p["<Q>"]<1. :
            return False , p
        
        if p["<z>"]<0.2:
            return False , p
        
        #### in the case of MAPFF I use special sets for pi- and K-
        if('MAPFF' in path_to_constants):            
            if(p["process"][3]==2101 and p["process"][2]<0): 
                p["process"][3]=2107
                p["process"][2]=3
            if(p["process"][3]==2103 and p["process"][2]<0): 
                p["process"][3]=2108
                p["process"][2]=3
            if(p["process"][3]==2105 and p["process"][2]<0): 
                p["process"][3]=2109
                p["process"][2]=3
            if(p["process"][2]==-1): p["process"][2]=3
            if(p["process"][2]==-2): p["process"][2]=4
            
        
        return (delta<0.5 or p["<pT>"]<1.) , p
        
    
    elif p["type"]=="DY":        
        #  for artemide v3.    
        # p["process"]=[p["process"][0],p["process"][2],1,1]
        if(len(p["process"])==3):        
                print("UNKNOWN PROCESS IN ARTEMIDE 3"+str(p["process"]))
                
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
        
        return (delta<0.35 or p["<qT>"]<10.) , p
       

#%%
### Loading the SIDIS data set
theData=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisDataSIDIS([
                      'hermes.p.vmsub.zxpt.pi+','hermes.p.vmsub.zxpt.pi-',
                      'hermes.d.vmsub.zxpt.pi+','hermes.d.vmsub.zxpt.pi-',
                      'hermes.p.vmsub.zxpt.k+','hermes.p.vmsub.zxpt.k-',
                      'hermes.d.vmsub.zxpt.k+','hermes.d.vmsub.zxpt.k-',
                      'compass.d.h+','compass.d.h-']))

setSIDIS=theData.CutData(cutFunc) 
setSIDIS_PLOT=theData.CutData(cutFunc_PLOT) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setSIDIS.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setSIDIS.sets]), 'points.')
print('Loaded ', setSIDIS_PLOT.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setSIDIS_PLOT.sets]), 'points.')

#%%
### Loading the DY data set
theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY([
                          'CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                          'A7-00y10', 'A7-10y20','A7-20y24', 
                          'A8-00y04', 'A8-04y08', 'A8-08y12', 'A8-12y16', 'A8-16y20', 'A8-20y24', 
                          'A8-46Q66', 'A8-116Q150', 
                          'A13-norm',
                          'CMS7', 'CMS8', 
                          'CMS13-00y04','CMS13-04y08','CMS13-08y12','CMS13-12y16','CMS13-16y24',
                          'CMS13_dQ_50to76',
                          'CMS13_dQ_106to170','CMS13_dQ_170to350','CMS13_dQ_350to1000',
                          'LHCb7', 'LHCb8', 'LHCb13_dy(2021)', 
                          'PHE200', 'STAR510', 
                          'E228-200', 'E228-300', 'E228-400', 
                          'E772',
                          'E605',
                          'D0run1-W','CDFrun1-W'
                          ]))

setDY=theData.CutData(cutFunc) 
setDY_PLOT=theData.CutData(cutFunc_PLOT) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded ', setDY_PLOT.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY_PLOT.sets]), 'points.')

#%%
rSet.SetReplica(0)

#listToSave=DataProcessor.harpyInterface.ComputeXSec(setSIDIS_PLOT)

# DataProcessor.harpyInterface.PrintChi2Table(setSIDIS,printDecomposedChi2=True)
# DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)

#%%

###############################################################################
###   Saving the central replica value. It is better to do with the presice artemide mode
###############################################################################
# rSet.SetReplica(0)

# PATH_TO_SAVE="/data/WorkingFiles/TMD/Fit_Notes/ART25/DataForPLOTS"

# pointInFIT=[p["id"] for p in setDY.points]

# with open(PATH_TO_SAVE+'/centralX_DY.dat', "w") as filehandle:
#     filehandle.write("ID, <Q>, Qmin, Qmax, <qT>, qTmin, qTmax, inFIT, Xsec, Err, ART25(central)"+'\n')

# listToSave=DataProcessor.harpyInterface.ComputeXSec(setDY_PLOT)
   
# with open(PATH_TO_SAVE+'/centralX_DY.dat', "a") as filehandle: 
#     i=0
#     for p in setDY_PLOT.points:
        
#         line=p["id"]\
#             +', '+'{:.2f}'.format(p["<Q>"])\
#             +', '+'{:.2f}'.format(p["Q"][0])\
#             +', '+'{:.2f}'.format(p["Q"][1])\
#             +', '+'{:.2f}'.format(p["<qT>"])\
#             +', '+'{:.2f}'.format(p["qT"][0])\
#             +', '+'{:.2f}'.format(p["qT"][1])\
#             +', '+'{:.6f}'.format(p["xSec"])\
#             +', '+'{:.6f}'.format(numpy.sqrt(numpy.sum(numpy.array(p["uncorrErr"])**2)))\
#             +', '+'{0}'.format(p["id"] in pointInFIT)\
#             +', '+'{:.6f}'.format(listToSave[i])\
#             +'\n'
#         filehandle.write(line)
        
#         i+=1

#%%

# ###  
# ###   Saving the central replica value. It is better to do with the presice artemide mode
# ###
# rSet.SetReplica(0)

# PATH_TO_SAVE="/data/WorkingFiles/TMD/Fit_Notes/ART25/DataForPLOTS"

# pointInFIT=[p["id"] for p in setSIDIS.points]

# with open(PATH_TO_SAVE+'/centralX_SIDIS.dat', "w") as filehandle:
#     filehandle.write("ID, <Q>, Qmin, Qmax, <x>, xmin, xmax, <z>, zmin, zmax, <pT>, pTmin, pTmax, Xsec, Err, inFIT, ART25(central)"+'\n')

# listToSave=DataProcessor.harpyInterface.ComputeXSec(setSIDIS_PLOT)
   
# with open(PATH_TO_SAVE+'/centralX_SIDIS.dat', "a") as filehandle: 
#     i=0
#     for p in setSIDIS_PLOT.points:
        
#         line=p["id"]\
#             +', '+'{:.4f}'.format(p["<Q>"])\
#             +', '+'{:.4f}'.format(p["Q"][0])\
#             +', '+'{:.4f}'.format(p["Q"][1])\
#             +', '+'{:.4f}'.format(p["<x>"])\
#             +', '+'{:.4f}'.format(p["x"][0])\
#             +', '+'{:.4f}'.format(p["x"][1])\
#             +', '+'{:.4f}'.format(p["<z>"])\
#             +', '+'{:.4f}'.format(p["z"][0])\
#             +', '+'{:.4f}'.format(p["z"][1])\
#             +', '+'{:.4f}'.format(p["<pT>"])\
#             +', '+'{:.4f}'.format(p["pT"][0])\
#             +', '+'{:.4f}'.format(p["pT"][1])\
#             +', '+'{:.8f}'.format(p["xSec"])\
#             +', '+'{:.8f}'.format(numpy.sqrt(numpy.sum(numpy.array(p["uncorrErr"])**2)))\
#             +', '+'{0}'.format(p["id"] in pointInFIT)\
#             +', '+'{:.8f}'.format(listToSave[i])\
#             +'\n'
#         filehandle.write(line)
        
#         i+=1
    
# sys.exit()

#%%
# # ##################################
# # ## Header of the REPLICA FILES
# # ##################################

PATH_TO_SAVE="/data/WorkingFiles/TMD/Fit_Notes/ART25/DataForPLOTS"

with open(PATH_TO_SAVE+'/replicas_SIDIS.dat', "w") as filehandle:
    filehandle.write(", ".join([p["id"] for p in setSIDIS_PLOT.points])+'\n')
    
with open(PATH_TO_SAVE+'/replicas_DY.dat', "w") as filehandle:
    filehandle.write(", ".join([p["id"] for p in setDY_PLOT.points])+'\n')

#%%
# # ##################################
# # ## Save replicas of xSec
# # ## Better to compute with lower-accuracy mode.
# # ##################################
for i in range(rSet.numberOfReplicas):      
    rSet.SetReplica(i)

    listToSave=DataProcessor.harpyInterface.ComputeXSec(setDY_PLOT)
    with open(PATH_TO_SAVE+'/replicas_DY.dat', "a") as filehandle:
        filehandle.write(", ".join(['{:.6f}'.format(ss) for ss in listToSave])+"\n")

    listToSave=DataProcessor.harpyInterface.ComputeXSec(setSIDIS_PLOT)
    with open(PATH_TO_SAVE+'/replicas_SIDIS.dat', "a") as filehandle:
        filehandle.write(", ".join(['{:.6f}'.format(ss) for ss in listToSave])+"\n")
    
    print(i)
    
