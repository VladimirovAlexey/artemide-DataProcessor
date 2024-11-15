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

ATMDE_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide-playground/"
HARPY_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide-playground/harpy/"


import sys
import numpy
sys.path.remove('/data/arTeMiDe_Repository/artemide/harpy')
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

path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/ART25_N4LL_resum_MAPFF_CS.atmde"


harpy.initialize(path_to_constants)

inARRAY_TMDR=[1.5004, 0.073018, 0.038048, 0.0,1.,1.,1.]
inARRAY_PDF=[0.521462, 0.000206, 0.402948, 7.0219, 1.0, 20.4051, 1.0, 0.000123, 1.1037, 0.660734, 0.0, 0.04]
inARRAY_FF=[0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0,0.0,0.0,0.0,0.0]


harpy.setNPparameters_TMDR(inARRAY_TMDR)
harpy.setNPparameters_uTMDPDF(inARRAY_PDF)
harpy.setNPparameters_uTMDFF(inARRAY_FF)

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
    
   

#%%
### Loading the SIDIS data set
theData=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisDataSIDIS([
                      'hermes.p.vmsub.zxpt.pi+','hermes.p.vmsub.zxpt.pi-',
                      'hermes.d.vmsub.zxpt.pi+','hermes.d.vmsub.zxpt.pi-',
                      'hermes.p.vmsub.zxpt.k+','hermes.p.vmsub.zxpt.k-',
                      'hermes.d.vmsub.zxpt.k+','hermes.d.vmsub.zxpt.k-',
                      'compass.d.h+','compass.d.h-']))

setSIDIS=theData.CutData(cutFunc) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setSIDIS.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setSIDIS.sets]), 'points.')

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

setDY=theData.CutData(cutFunc) 
#setDYfit=theData.CutData(cutFuncFORFIT) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')

#%%
harpy.setNPparameters([1.5, 0.073, 0.03805, 0.0,10.0,1.,1. ,
                       0.526, 0.0001, 0.465, 1.46, 
                       1.0, 19.86, 1.0, 0.00026, 
                       0.819, 19.88, 0.0, 0.04, 
                       0.573, 0.418, 0.2450, 0.540, 
                       0.869, 1.14, -3.627, 1.4656,
                       0.0,-1.303,0.0,0.0])

#DataProcessor.harpyInterface.PrintChi2Table(setSIDIS,printDecomposedChi2=True)
DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)

#%%
chiBASE=1.2278753451919004
#tt=[0.05,0.075,0.1,0.125,0.15,0.175]
tt=[-0.05,-0.075,-0.1,-0.125,-0.15,-0.175]
for ll in tt:
    harpy.setNPparameters_TMDR([1.5, 0.073, 0.03805, -ll , 1.5,-1.,1.])
    YY=DataProcessor.harpyInterface.ComputeXSec(setDY)
    chi1,cc3=setDY.chi2(YY)
    print("{",ll,",",chi1/627-chiBASE,"},")

#%%
saveFile="/data/WorkingFiles/TMD/Fit_Notes/ART25/CS_deltaCHI.dat"
bTABLE=[0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.5,3.,4.,5.]
alphaTABLE=[0.5,0.75,1.0,1.5,2.,3.]
tt=[0.05,0.075,0.1,0.125,0.15,0.175]
for b in bTABLE:
    for alpha in alphaTABLE:
        a=0.05
        
        harpy.setNPparameters_TMDR([1.5, 0.073, 0.03805, a, b,-1.,alpha])
        YY=DataProcessor.harpyInterface.ComputeXSec(setDY)
        chi1,cc3=setDY.chi2(YY)
        RR=chi1/627-chiBASE 
        
        f=open(saveFile,"a+")
        f.write(str(b)+", "+str(alpha)+", "+str(a)+", "+str(RR)+"\n")
        f.close()  
        
        while RR<0.02:
            ss=1.
            if(b<2): ss=0.5
            if(b<1): ss=0.25
            
            if(RR<0.0001): a+=0.2*ss
            elif(RR<0.001): a+=0.1*ss
            elif(RR<0.005): a+=0.05*ss
            else: a+=0.025*ss           
            
            print(b,alpha,a)
            
            harpy.setNPparameters_TMDR([1.5, 0.073, 0.03805, a, b,-1.,alpha])
            YY=DataProcessor.harpyInterface.ComputeXSec(setDY)
            chi1,cc3=setDY.chi2(YY)
            RR=chi1/627-chiBASE 
            
            f=open(saveFile,"a+")
            f.write(str(b)+", "+str(alpha)+", "+str(a)+", "+str(RR)+"\n")
            f.close()  
            
#%%
saveFile="/data/WorkingFiles/TMD/Fit_Notes/ART25/CS_deltaCHI_negative.dat"
bTABLE=[0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.5,3.,4.,5.]
alphaTABLE=[0.5,0.75,1.0,1.5,2.,3.]
for b in bTABLE:
    for alpha in alphaTABLE:
        a=-0.01
        
        harpy.setNPparameters_TMDR([1.5, 0.073, 0.03805, a, b,-1.,alpha])
        YY=DataProcessor.harpyInterface.ComputeXSec(setDY)
        chi1,cc3=setDY.chi2(YY)
        RR=chi1/627-chiBASE 
        
        f=open(saveFile,"a+")
        f.write(str(b)+", "+str(alpha)+", "+str(a)+", "+str(RR)+"\n")
        f.close()  
        
        while RR<0.02:
            ss=0.5
            if(b<2): ss=0.25
            if(b<1): ss=0.1
            
            
            if(RR<0.0001): a-=0.2*ss
            elif(RR<0.001): a-=0.1*ss
            elif(RR<0.005): a-=0.05*ss
            else: a-=0.025*ss           
            
            print(b,alpha,a)
            
            harpy.setNPparameters_TMDR([1.5, 0.073, 0.03805, a, b,-1.,alpha])
            YY=DataProcessor.harpyInterface.ComputeXSec(setDY)
            chi1,cc3=setDY.chi2(YY)
            RR=chi1/627-chiBASE 
            
            f=open(saveFile,"a+")
            f.write(str(b)+", "+str(alpha)+", "+str(a)+", "+str(RR)+"\n")
            f.close()            
