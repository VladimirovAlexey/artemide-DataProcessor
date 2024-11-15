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

logFile=MAINPATH+"FittingPrograms/ART25/log.txt"
replicaFile ="/data/WorkingFiles/TMD/Fit_Notes/ART25/REPLICAS/attempt0_1.dat"
#%%
#######################################
#Initialize artemide
#######################################
import harpy

#path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/ART25_N4LL.atmde"
#path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/ART25_N4LL_DYresum.atmde"
path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/ART25_N4LL_resum_MAPFF.atmde"
#path_to_constants=MAINPATH+"FittingPrograms/ART25/ConstantsFiles/ART25_N4LL_NOresum_MAPFF.atmde"


harpy.initialize(path_to_constants)

inARRAY_TMDR=[1.5004, 0.073018, 0.038048, 0.0]
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
harpy.setNPparameters([1.5, 0.073, 0.038, 0.0, 
                       0.497, 0.0001, 0.4608, 1.5041, 
                       0.0, 20.24, 0.0, 0.0001, 
                       0.844, 20.20, 0.0, 0.04, 
                       0.573, 0.418, 0.2450, 0.540, 
                       0.869, 1.14, -3.627, 1.4656,
                       0.0,-1.303,0.0,0.0])

DataProcessor.harpyInterface.PrintChi2Table(setSIDIS,printDecomposedChi2=True)
DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)

#%%
#### Minimize SIDIS
from iminuit import Minuit

#---- PDFbias-like row
initialValues=([
    1.5, 0.071624, 0.064726, 0.0,
    0.568631, 1.1114, 0.563709, 7.2062, 
    4.5436, 17.3135, 0.976618, 0.006934, 
    1.3213, 25.9727, 0.0, 0.04, 
    0.674236, 0.781384, 0.441778, 0.941478, 
    0.509997, 0.262662, -3.073866, 0.198287, 
    0.0, -0.37617, 0.0, 0.0])


initialErrors=(0.1,0.1,0.1,0.1,
                0.1,  1.0, 0.1,  1.0,
                0.1,  1.0, 0.1,  1.0,
                0.1,  1.0, 10.,  1.,
                0.5,0.5,0.5,0.5,
                0.5,0.5,0.5,0.5,
                0.5,0.5,0.5,0.5)
searchLimits=((1.0,2.5),(0.005,0.15) ,(0.0,.2), (-5.,5.),
              (0.00001,100.), (0.00001,100.),(0.00001,100.),(0.00001,100.),
              (-100.,100.), (0.00001,100.),(-100.,100.),(0.00001,100.),
              (0.00001,100.), (0.00001,100.),(0.,100.),(0.04,100.),
              (0.0001,100.), (-100.,100.),(-100.,100.),(-100.,100.),
              (0.0001,100.), (-100.,100.),(-100.,100.),(-100.,100),
              (-100.,100.),(-100.,100.),(-100.,100.),(-100.,100.))
              
# True= FIX
parametersToMinimize=(True, False,False,True,
                      False, False, False,False,
                      False, False, False, False,
                      False, False, True,True,
                      False, False, False,False,
                      False, False, False,False,
                      True, False, True,True)

#%%
#######################################
# Generate replica of data and compute chi2
#######################################
import time

penalty_index=[-7,-6,-5,-4,-3]

def MinForReplica():
    global setDY,setSIDIS,initialValues,initialErrors,searchLimits,parametersToMinimize,penalty_index
        
    def repchi_2(x):     
        startT=time.time()
        #harpy.setNPparameters_uTMDFF([x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7]])
        harpy.setNPparameters(x)
        print('np set =',["{:8.3f}".format(i) for i in x])        
        
        YY=DataProcessor.harpyInterface.ComputeXSec(setSIDISrep)
        ccSIDIS2,cc3=setSIDISrep.chi2(YY)
            
        YY=DataProcessor.harpyInterface.ComputeXSec(setDYrep)
        ccDY2,cc3=setDYrep.chi2(YY)
        
        penalty_array=numpy.array([max(0,abs(setDYrep.sets[i].DetermineAvarageSystematicShift(YY[setDYrep._i1[i]:setDYrep._i2[i]]))/setDY.sets[i].normErr[0]-1) for i in penalty_index])
        penalty_term=sum(penalty_array**6)
        
        cc=[ccSIDIS2/setSIDISrep.numberOfPoints,ccDY2/setDYrep.numberOfPoints, 
            (ccSIDIS2+ccDY2)/(setSIDISrep.numberOfPoints+setDYrep.numberOfPoints)]
        
        endT=time.time()
        print(':->',cc,'   +p=',penalty_term,"    time=",endT-startT)
        return ccSIDIS2+ccDY2+penalty_term

    
    setDYrep=setDY.GenerateReplica()
    setSIDISrep=setSIDIS.GenerateReplica()
    
    localM = Minuit(repchi_2, initialValues)
    
    localM.errors=initialErrors
    localM.limits=searchLimits
    localM.fixed=parametersToMinimize
    localM.errordef=1    
    ### tolerance is a bit larger than for the main fit
    localM.tol=0.0005*(setSIDISrep.numberOfPoints+setDYrep.numberOfPoints)*10000 
    localM.strategy=1

    localM.migrad()
    
    ### [chi^2, NP-parameters]
    return [localM.fval,list(localM.values)]

#%%
#######################################
# LOG save function
#######################################
savedTime=time.time()
def SaveToLog(text):
    global savedTime,logFile
    newTime=time.time()
    
    import socket
    PCname=socket.gethostname()
    
    passedTime=newTime-savedTime
    hours=int(passedTime/3600)
    minutes=int((passedTime-hours*3600)/60)
    seconds=int(passedTime-hours*3600-minutes*60)
    
    with open(logFile, 'a') as file:
        file.write(PCname+ ' : ' + time.ctime()+' :  [+'+str(hours)+':'+str(minutes)+':'+str(seconds)+' ]\n')
        file.write(' --> '+text+'\n')
        file.write('\n')
    savedTime=time.time()

#%%
#######################################
# This is the main cicle. 
# It generates replica of data take random PDF and minimize it
# Save to log.
#######################################

NumberOfReplicas=50

SaveToLog(" ================================================================================")
SaveToLog(" ================================================================================")
SaveToLog(" ================================================================================")

for i in range(NumberOfReplicas):
    print('---------------------------------------------------------------')
    print('------------REPLICA ',i,' from ',NumberOfReplicas,'------------------')
    print('---------------------------------------------------------------')
    savedTime=time.time()
    
    ## reset PDF        
    PDFreplica=numpy.random.randint(1000)
    FFreplicaPI=numpy.random.randint(200)
    FFreplicaK=numpy.random.randint(200)
    harpy.setPDFreplica(PDFreplica)
    harpy.setFFreplica(FFreplicaPI)
    harpy.setFFreplica(FFreplicaPI)
    harpy.setFFreplica(FFreplicaK)
    harpy.setFFreplica(FFreplicaK)
    print("Start computation with PDF/FF replicas "+str(PDFreplica)+", "+str(FFreplicaPI)+", "+str(FFreplicaK))
    
    harpy.setNPparameters(initialValues)    
    print("Minimization started.")
    
    ## got to pseudo-data and minimization
    repRes=MinForReplica()
    print(repRes)
    print("Minimization finished.")    
    SaveToLog(" Minization with PDF/FF replicas "+str(PDFreplica)+", "+str(FFreplicaPI)+", "+str(FFreplicaK)+" completed.")
    
    ## compute the chi2 for true data full
    mainDY, mainDY2 =DataProcessor.harpyInterface.ComputeChi2(setDY)    
    mainSIDIS, mainSIDIS2 =DataProcessor.harpyInterface.ComputeChi2(setSIDIS)    
    print("Central chi^2  computed.)")
    
    ## save to file
    f=open(replicaFile,"a+")
    print('SAVING >>  ',f.name)
    ### [total chi^2(full data DY),total chi^2(full data SIDIS), 
    ### list of chi^2 for experiments( DY),list of chi^2 for experiments( SIDIS), 
    ### PDFreplica, FFpi-replica, FFk-replica,
    ### list of NP-parameters] 
    f.write(str([mainDY,mainSIDIS,mainDY2,mainSIDIS2,PDFreplica,FFreplicaPI,FFreplicaK, repRes[1]])+"\n")
    f.close()  
