#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: vla18041
"""

#%%
#######################################################################
# Global parameter of a run
#######################################################################

#PDFinUse="HERA20"
PDFinUse="NNPDF31"
#PDFinUse="CT18"

## if this trigger is ON the LHC data will be fit only by shape
useNormalizedLHCdata=False
## Include ATLAS 7TeV?
useA7data=False
## Split the low-energy experiment <Upsilon and >Upsilon
splitUpsilon=True
## Use the reduced set of the data
useReducedSet=True

## total number of parameters in TMDPDF (7 and 12 case are defined)
numberOfParameters=12

#### Starting and final replica (included)
StartReplica=1
FinalReplica=3

## automatic name generation for the run
runName="model4.0_"+PDFinUse+"_PDFrep_"
if (not useA7data): runName+="noA7_"
if (splitUpsilon): runName+="spltUPS_"
if (useNormalizedLHCdata): runName+="norm_"
if (useReducedSet): runName+="reducedSet_"
if (runName[-1]=="_"): runName=runName[0:-1]

print(" RUN: "+runName)

#%%
#######################################
# Paths
#######################################
PathToHarpy="/home/vla18041/LinkData2/arTeMiDe_Repository/artemide-ForPDF/harpy"
PathToDataProcessor="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/"
PathToDataLibrary=PathToDataProcessor+"DataLib/unpolDY/"
PathToLog=PathToDataProcessor+"FittingPrograms/PDF-TMD/LOGS/"
PathToSavings=PathToDataProcessor+"FittingPrograms/PDF-TMD/LOGS/"+runName+"/"
if(numberOfParameters==7):
    PathToConstantsFile=PathToDataProcessor+"/FittingPrograms/PDF-TMD/Constants-files/const-"+PDFinUse+"_NNLO_7p"
elif(numberOfParameters==12):
    PathToConstantsFile=PathToDataProcessor+"/FittingPrograms/PDF-TMD/Constants-files/const-"+PDFinUse+"_NNLO_12p"
else:
    print("UNKNOWN NUMBER OF PARAMETERS")

import sys
sys.path.remove('/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/harpy')
sys.path.append(PathToHarpy)
sys.path.append(PathToDataProcessor)

#%%
#######################################
# Output paths
#######################################
import socket
PCname=socket.gethostname()

replicaFile=PathToLog+runName+".txt"
logFile=PathToLog+PCname+".log"

#%%
#######################################
# importing libraries
#######################################
import numpy
import time

import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet

#%%
#######################################
# LOG save function
#######################################
savedTime=time.time()
def SaveToLog(text):
    global savedTime,logFile
    newTime=time.time()
    
    passedTime=newTime-savedTime
    hours=int(passedTime/3600)
    minutes=int((passedTime-hours*3600)/60)
    seconds=int(passedTime-hours*3600-minutes*60)
    
    with open(logFile, 'a') as file:
        file.write(time.ctime()+' :  [+'+str(hours)+':'+str(minutes)+':'+str(seconds)+' ]\n')
        file.write(' run '+runName+'\n')
        file.write(' --> '+text+'\n')
        file.write('\n')
    savedTime=time.time()

#%%
#######################################
#Initialize artemide
#######################################
import harpy

SaveToLog('Initialization with : \n'+PathToConstantsFile)

harpy.initialize(PathToConstantsFile)
if(numberOfParameters==7):
    initializationArray=[2.0340, 0.0299, 0.2512, 7.7572, 0.2512, 7.7572, 0.2512, 7.7572, 10000.]
elif(numberOfParameters==12):
    initializationArray=[2.0340, 0.0299, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,0.05, 0.05, 0.,0.5]
else:
    print("UNKNOWN NUMBER OF PARAMETERS")

harpy.setNPparameters(initializationArray)

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
#######################################
# Data cut function
#######################################
dropPoints=["A7-10y20.7", "A7-10y20.8", "A7-10y20.9", "A7-20y24.4", "A7-20y24.5", \
"A7-20y24.6", "A7-20y24.7", "A7-20y24.8", "A7-20y24.9", "A8-08y12.8", \
"A8-116Q150.1", "A8-116Q150.2", "A8-116Q150.3", "A8-116Q150.4", \
"A8-116Q150.5", "A8-116Q150.6", "A8-116Q150.7", "A8-116Q150.8", \
"A8-116Q150.9", "A8-12y16.8", "A8-16y20.6", "A8-16y20.7", \
"A8-16y20.8", "A8-20y24.5", "A8-20y24.6", "A8-20y24.7", "A8-20y24.8", \
"CDF1.10", "CDF1.11", "CDF1.12", "CDF1.13", "CDF1.14", "CDF1.15", \
"CDF1.16", "CDF1.17", "CDF1.18", "CDF1.19", "CDF1.20", "CDF1.21", \
"CDF1.22", "CDF1.23", "CDF1.24", "CDF1.25", "CDF1.26", "CDF1.27", \
"CDF1.28", "CDF1.29", "CDF1.3", "CDF1.30", "CDF1.31", "CDF1.32", \
"CDF1.4", "CDF1.5", "CDF1.6", "CDF1.7", "CDF1.8", "CDF1.9", \
"CDF2.15", "CDF2.16", "CDF2.17", "CDF2.18", "CDF2.19", "CDF2.20", \
"CDF2.21", "CDF2.22", "CDF2.23", "CDF2.24", "CDF2.25", "CDF2.26", \
"CDF2.27", "CDF2.28", "CDF2.29", "CDF2.30", "CDF2.31", "CDF2.32", \
"CDF2.33", "CDF2.34", "CDF2.35", "CDF2.36", "CDF2.37", "CDF2.38", \
"CDF2.39", "CDF2.40", "CDF2.41", "CDF2.42", "CDF2.43", "CDF2.44", \
"CMS7.1", "CMS7.2", "CMS7.3", "CMS7.4", "CMS7.5", "CMS7.6", "CMS7.7", \
"CMS8.0", "CMS8.1", "CMS8.2", "CMS8.3", "CMS8.4", "CMS8.5", "CMS8.6", \
"CMS8.7", "D01.10", "D01.11", "D01.12", "D01.13", "D01.14", "D01.15", \
"D01.2", "D01.3", "D01.4", "D01.5", "D01.6", "D01.7", "D01.8", \
"D01.9", "D02.1", "D02.2", "D02.3", "D02.4", "D02.5", "D02.6", \
"D02.7", "D02.8", "D02m.3", "D02m.4", "E228-200.8Q9.10", \
"E228-200.8Q9.5", "E228-200.8Q9.6", "E228-200.8Q9.9", \
"E228-300.11Q12.11", "E228-300.11Q12.3", "E228-300.11Q12.4", \
"E228-300.11Q12.5", "E228-300.11Q12.6", "E228-300.11Q12.7", \
"E228-300.11Q12.8", "E228-300.8Q9.10", "E228-400.13Q14.11", \
"E228-400.13Q14.12", "E228-400.13Q14.6", "E228-400.13Q14.7", \
"E228-400.13Q14.9", "E605.7Q8.7", "E605.7Q8.8", "E605.7Q8.9", \
"E772.12Q13.10", "E772.12Q13.6", "E772.12Q13.8", "E772.12Q13.9", \
"E772.13Q14.4", "E772.13Q14.5", "E772.14Q15.0", "E772.14Q15.3", \
"E772.14Q15.4", "E772.14Q15.5", "E772.14Q15.6", "LHCb13.1", \
"LHCb13.10", "LHCb13.2", "LHCb13.3", "LHCb13.4", "LHCb13.5", \
"LHCb13.6", "LHCb13.7", "LHCb13.8", "LHCb13.9", "LHCb7.10", \
"LHCb7.2", "LHCb7.3", "LHCb7.4", "LHCb7.7", "LHCb7.8", "LHCb7.9", \
"LHCb8.10", "LHCb8.8", "LHCb8.9", "PHE200.2",\
##### extra points by hands
"A8-116Q150.0","CDF1.0","CDF1.1","CDF1.2","D01.0", "D01.1", "D02.0","CMS7.0","LHCb13.0"]

#### Check the point kinematics
def cutFunc0(p):    
    
    par=1.0

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
def cutFunc(p):    
    
    #### check against the presence in the reduced set.
    if(useReducedSet and p["id"] in dropPoints):
        return False,p
    
    return cutFunc0(p)



#%%
#######################################
# Loading the data set
#######################################

setHE=loadThisData(['CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                      'A7-00y10' if useA7data else '',
                      'A7-10y20' if useA7data else '',
                      'A7-20y24' if useA7data else '', 
                      'A8-00y04-norm' if useNormalizedLHCdata else 'A8-00y04',
                      'A8-04y08-norm' if useNormalizedLHCdata else 'A8-04y08',
                      'A8-08y12-norm' if useNormalizedLHCdata else 'A8-08y12',
                      'A8-12y16-norm' if useNormalizedLHCdata else 'A8-12y16',
                      'A8-16y20-norm' if useNormalizedLHCdata else 'A8-16y20',
                      'A8-20y24-norm' if useNormalizedLHCdata else 'A8-20y24',
                      'A8-46Q66-norm' if useNormalizedLHCdata else 'A8-46Q66',
                      'A8-116Q150-norm' if useNormalizedLHCdata else 'A8-116Q150',
                      'CMS7', 'CMS8', 
                      'LHCb7', 'LHCb8', 'LHCb13'])

#### If I separate data above and below UPSILON, I create two copies of LE data with different names
#### the data to be split only  'E228-300', 'E228-400' and E605
if(splitUpsilon):
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
else:
    setLE=loadThisData(['PHE200', 'E228-200', 'E228-300', 'E228-400','E772','E605'])

theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",setHE+setLE)

setDY=theData.CutData(cutFunc) 

setDYfull=theData.CutData(cutFunc0) 

print('Loaded ', setDY.numberOfSets, 'data sets with ', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setDY.sets])

SaveToLog('Loaded '+ str(setDY.numberOfSets) + ' data sets with '+str(sum([i.numberOfPoints for i in setDY.sets])) + ' points. \n'
+'Loaded experiments are '+str([i.name for i in setDY.sets]))

#%%
if useNormalizedLHCdata:
    for s in setDY.sets:
        if s.name[0:4]=='LHCb':
            s.isNormalized=True
            
        if s.isNormalized:
            s.normalizationMethod='bestChi2'
#%%

if(PDFinUse=="HERA20"):
    #initialValues=(2.000,  0.033, 0.230, 5.609, 0.252, 8.021, 570.223, 0.000, 0.000) #model 1.0 HERA
    #initialValues=(2.000,  0.033, 0.230, 5.609, 0.252, 8.021, 0.252, 8.021, 570.223) #model 2.0 HERA
    #initialValues=(2.000, 0.034, 0.145, 10.393, 0.333, 0.000, 0.366, 10.261, 677.434) #model 2.2 HERA
    #initialValues=(2.000, 0.036, 0.089, 8.657,  0.389, 0.000, 0.465, 6.195, 527.903) #model 2.2 HERA reduced
    initialValues=(2.000, 0.036, 0.099, 8.578, 0.396, 0.246, 0.089, 8.604, 0.379, 0.000, 0.455, 5.355, 528.381, 0.000) #model 4.0 HERA
if(PDFinUse=="NNPDF31"):
    #initialValues=(2.000,  0.029, 0.345, 2.587, 0.152, 7.561, 232.544, 0. , 0.)      #model 1.0 NNPDF31
    #initialValues=(2.000,  0.029, 0.345, 2.587, 0.152, 7.561, 0.152, 7.561, 232.544) #model 2.0 NNPDF31
    #initialValues=(2.000, 0.030, 0.188, 5.542, 0.200, 4.375, 0.486, 0.009, 232.793) #model 2.2 NNPDF31
    #initialValues=(2.000, 0.035, 0.160, 5.139, 0.219, 3.711, 0.460, 0.038, 229.333) #model 2.2 NNPDF31 reduced
    initialValues=(2.000, 0.034, 0.330, 1.416, 0.371, 2.504, 0.128, 9.149, 0.128, 5.158, 0.195, 4.758, 192.312, 0.000) #model 4.0 NNPDF31
if(PDFinUse=="CT18"):
    #initialValues=(2.000,  0.039, 0.161, 7.904, 0.212, 5.301, 700.642, 0. , 0.)      #model 1.0 CT18
    #initialValues=(2.000,  0.039, 0.161, 7.904, 0.212, 5.301, 0.212, 5.301, 700.642) #model 2.0 CT18
    #initialValues=(2.000, 0.042, 0.094, 12.534, 0.293, 0.004, 0.003, 16.568, 819.267) #model 2.2 CT18
    #initialValues=(2.000, 0.042, 0.121, 9.626,  0.289, 0.000, 0.000, 17.820, 524.722) #model 2.2 CT18 reduced
    #initialValues=(2.000, 0.042, 0.121, 9.626,  0.289, 0.000, 0.121, 9.626,  0.289, 0.000, 0.000, 17.820, 524.722, 0.) #model 4.0 CT18
    initialValues=(2.000, 0.048, 0.063, 0.009, 0.414, 4.566, 0.049, 24.999, 0.000, 0.213, 0.002, 24.326, 0.036, 0.000) #model 4.0 CT18
harpy.setNPparameters(list(initialValues))

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)

#%%
#### translation of parameters to list in proper order

p_names=['BNP','c0','l1','l2','l3','l4','l5','l6','l7','l8','l9','l10','l11','l12']
def params_to_array(params):
    return [params[p_names[i]].value for i in range(len(params))]

#%%
#######################################
# Main chi2 formula
#######################################
totalN=setDY.numberOfPoints

def chi_2(pars):
    startT=time.time()
    
    x=params_to_array(pars)
    harpy.setNPparameters(x)
    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
    
    ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)
    
    cc=ccDY2/totalN
    endT=time.time()
    print(':->',cc,'       t=',endT-startT)
    return ccDY2


#%%
from lmfit import Parameters, fit_report, minimize

##### for test I start from CT18 case
initialValues=(2.000, 0.048, 0.063, 0.009, 0.414, 4.566, 0.049, 24.999, 0.000, 0.213, 0.002, 24.326, 0.036, 0.000)

fit_params = Parameters()
fit_params.add_many((p_names[0] , initialValues[0] , False, 1.9, 2.1, None, None),
                    (p_names[1] , initialValues[1] , True , 0.0001, 1.0, None, None),
                    (p_names[2] , initialValues[2] , True , 0.001, 3.0, None, None),
                    (p_names[3] , initialValues[3] , True , 0.001, 25.0, None, None),
                    (p_names[4] , initialValues[4] , True , 0.001, 3.0, None, None),
                    (p_names[5] , initialValues[5] , True , 0.001, 25.0, None, None),
                    (p_names[6] , initialValues[6] , True , 0.001, 3.0, None, None),
                    (p_names[7] , initialValues[7] , True , 0.001, 25.0, None, None),
                    (p_names[8] , initialValues[8] , True , 0.001, 3.0, None, None),
                    (p_names[9] , initialValues[9] , True , 0.001, 25.0, None, None),
                    (p_names[10], initialValues[10], True , 0.001, 3.0, None, None),
                    (p_names[11], initialValues[11], True , 0.001, 25.0, None, None),
                    (p_names[12],  initialValues[12], True , 0.001, 2500., None, None),
                    (p_names[13], initialValues[13], False, None, None, None, None)
                    )

#%%
#out = minimize(chi_2, fit_params, method= 'Powell',options={'ftol' : 0.0001})
#out = minimize(chi_2, fit_params, method= 'Nelder-Mead',options={'fatol' : 0.01})
out = minimize(chi_2, fit_params, method= 'COBYLA',options={'tol' : 0.0001})


DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)

print(fit_report(out))

#print('------POWELL--')
#print('------Nelder-Mead--')
print('------COBYLA--')

sys.exit()

#%%
#######################################
# Generate replica of data and compute chi2
#######################################
def MinForReplica():
    global totalN,setDY,initialValues,initialErrors,searchLimits,parametersToMinimize
        
    def repchi_2(x):        
        global totalN
        startT=time.time()
        harpy.setNPparameters(x)
        print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
        
        ccDY2,cc2=DataProcessor.harpyInterface.ComputeChi2(repDataDY)
        
        cc=(ccDY2)/totalN
        endT=time.time()
        print(':->',cc,'       t=',endT-startT)
        return ccDY2
    
    repDataDY=setDY
    
    # localM = Minuit.from_array_func(repchi_2, initialValues,
    #   error=initialErrors, limit=searchLimits, fix=parametersToMinimize, errordef=1)
    
    # localM.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
    # localM.strategy=1

    # localM.migrad()
    
    ### [chi^2, NP-parameters]
    # return [localM.fval,localM.values.values()]

#%%
#######################################
# Save plot of TMD
#######################################
xValues=[
0.0001, 0.00011, 0.00013, 0.00014, 0.00016, 0.00018, 0.0002, 0.00022, \
0.00025, 0.00028, 0.00032, 0.00035, 0.0004, 0.00045, 0.0005, 0.00056, 0.00063, \
0.00071, 0.00079, 0.00089,\
0.001, 0.0011, 0.0013, 0.0014, 0.0016, 0.0018, 0.002, 0.0022, \
0.0025, 0.0028, 0.0032, 0.0035, 0.004, 0.0045, 0.005, 0.0056, 0.0063, \
0.0071, 0.0079, 0.0089, 0.01, 0.011, 0.013, 0.014, 0.016, 0.018, \
0.02, 0.022, 0.025, 0.028, 0.032, 0.035, 0.04, 0.045, 0.05, 0.056, \
0.063, 0.071, 0.079, 0.089, 0.1, 0.11, 0.13, 0.14, 0.16, 0.18, 0.2, \
0.22, 0.25, 0.28, 0.32, 0.35, 0.4, 0.45, 0.5, 0.56, 0.63, 0.71, 0.79, \
0.89, 0.99]
bValues=[0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, \
0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.1, 1.2, 1.3, 1.4, 1.5, \
1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, \
3., 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4., 4.5, 5., 5.5, \
6., 6.5, 7., 7.5, 8.]


##### check the existance of files
from pathlib import Path
Path(PathToSavings).mkdir(parents=True, exist_ok=True)

## check file a and if absent, create it and add first lines x, and b
def TryFile(a):    
    import os.path
    if(os.path.isfile(PathToSavings+a)):
        return
    else:
        file=open(PathToSavings+a,"w")
        line=[]
        for x in xValues:
            for b in bValues:
               line.append(x)
        file.write(str(line)+"\n")
        line=[]
        for x in xValues:
            for b in bValues:
               line.append(b)
        file.write(str(line)+"\n")
        file.close()
        
TryFile("u-optimal.plt3d")
TryFile("d-optimal.plt3d")
TryFile("s-optimal.plt3d")
TryFile("uBar-optimal.plt3d")
TryFile("dBar-optimal.plt3d")
TryFile("sBar-optimal.plt3d")

TryFile("u-2GeV.plt3d")
TryFile("d-2GeV.plt3d")
TryFile("s-2GeV.plt3d")
TryFile("uBar-2GeV.plt3d")
TryFile("dBar-2GeV.plt3d")
TryFile("sBar-2GeV.plt3d")

TryFile("u-90GeV.plt3d")
TryFile("d-90GeV.plt3d")
TryFile("s-90GeV.plt3d")
TryFile("uBar-90GeV.plt3d")
TryFile("dBar-90GeV.plt3d")
TryFile("sBar-90GeV.plt3d")

def AddLine(fout,lineToAdd):
    file=open(PathToSavings+fout,"a+")
    file.write(str(lineToAdd)+"\n")
    file.close()

def SavePlot():
    
    #### optimal
    u2=[]
    d2=[]
    s2=[]
    ubar2=[]
    dbar2=[]
    sbar2=[]
    
    mu=-1.
    for x in xValues:
        for b in bValues:
            tmd=harpy.get_uTMDPDF(x,b,1,mu=mu)
            u2.append(tmd[5+2])
            d2.append(tmd[5+1])
            s2.append(tmd[5+3])
            ubar2.append(tmd[5-2])
            dbar2.append(tmd[5-1])
            sbar2.append(tmd[5-3])
    AddLine("u-optimal.plt3d",u2)
    AddLine("d-optimal.plt3d",d2)
    AddLine("s-optimal.plt3d",s2)
    AddLine("uBar-optimal.plt3d",ubar2)
    AddLine("dBar-optimal.plt3d",dbar2)
    AddLine("sBar-optimal.plt3d",sbar2)
    
    
    #### 2GeV
    u2=[]
    d2=[]
    s2=[]
    ubar2=[]
    dbar2=[]
    sbar2=[]
    
    mu=2.
    for x in xValues:
        for b in bValues:
            tmd=harpy.get_uTMDPDF(x,b,1,mu=mu)
            u2.append(tmd[5+2])
            d2.append(tmd[5+1])
            s2.append(tmd[5+3])
            ubar2.append(tmd[5-2])
            dbar2.append(tmd[5-1])
            sbar2.append(tmd[5-3])
    AddLine("u-2GeV.plt3d",u2)
    AddLine("d-2GeV.plt3d",d2)
    AddLine("s-2GeV.plt3d",s2)
    AddLine("uBar-2GeV.plt3d",ubar2)
    AddLine("dBar-2GeV.plt3d",dbar2)
    AddLine("sBar-2GeV.plt3d",sbar2)
    
    #### 90GeV
    u2=[]
    d2=[]
    s2=[]
    ubar2=[]
    dbar2=[]
    sbar2=[]
    
    mu=90.
    for x in xValues:
        for b in bValues:
            tmd=harpy.get_uTMDPDF(x,b,1,mu=mu)
            u2.append(tmd[5+2])
            d2.append(tmd[5+1])
            s2.append(tmd[5+3])
            ubar2.append(tmd[5-2])
            dbar2.append(tmd[5-1])
            sbar2.append(tmd[5-3])
    AddLine("u-90GeV.plt3d",u2)
    AddLine("d-90GeV.plt3d",d2)
    AddLine("s-90GeV.plt3d",s2)
    AddLine("uBar-90GeV.plt3d",ubar2)
    AddLine("dBar-90GeV.plt3d",dbar2)
    AddLine("sBar-90GeV.plt3d",sbar2)
    
#%%
#######################################
# Save plot of xSec
#######################################
##### check the existance of files
from pathlib import Path
Path(PathToSavings).mkdir(parents=True, exist_ok=True)
import os.path
if(os.path.isfile(PathToSavings+"xSec.plt")):
    pass
else:
    file=open(PathToSavings+"xSec.plt","w")
    line=[p["id"] for p in setDYfull.points]
    file.write(str(line)+"\n")
    line=[p["qT"] for p in setDYfull.points]
    file.write(str(line)+"\n")
    line=[p["xSec"] for p in setDYfull.points]
    file.write(str(line)+"\n")
    line=[numpy.sqrt(numpy.sum(numpy.array(p["uncorrErr"])**2)) for p in setDYfull.points]
    file.write(str(line)+"\n")
    file.close()
    
        
def SavePlotData():
    
    dd=DataProcessor.harpyInterface.ComputeXSec(setDYfull)
    AddLine("xSec.plt",dd)
#%%
#######################################
# This is the main cicle. 
# It generates replica of data take random PDF and minimize it
# Save to log.
#######################################
for i in range(StartReplica,FinalReplica+1):
    print('---------------------------------------------------------------')
    print('------------REPLICA ',i,' from [',StartReplica,' , ',FinalReplica,']------------------')
    print('---------------------------------------------------------------')
    
    ## reset PDF
    harpy.setNPparameters(initializationArray)    
    harpy.setPDFreplica(i)
    SaveToLog("Start computation of replica "+str(i) +"in ["+ str(StartReplica)+','+str(FinalReplica)+"]")
    
    ## got to pseudo-data and minimization
    repRes=MinForReplica()
    print(repRes)
    SaveToLog("Minimization for replica "+str(i) +"in ["+ str(StartReplica)+','+str(FinalReplica)+"] finished.")    
    
    ## compute the chi2 for true data full
    mainDY, mainDY2 =DataProcessor.harpyInterface.ComputeChi2(setDYfull)    
    SaveToLog("Central chi^2 for "+str(i) +"in ["+ str(StartReplica)+','+str(FinalReplica)+" computed. \n Saving to log >> "+replicaFile)
    
    ## save to file
    f=open(replicaFile,"a+")
    print('SAVING >>  ',f.name)
    ### [total chi^2(cenral), total chi^2 (pseudo data), list of chi^2 for experiments(central), number of PDF, list of NP-parameters]
    f.write(str([mainDY,repRes[0],mainDY2,i,repRes[1]])+"\n")
    f.close()
    
    print('SAVING >>  PLOTS')
    
    SavePlot()
    SavePlotData()
    
    
#%%
#######################################
# Finalizing log
#######################################
print('Computation finished')
SaveToLog("Computation finished correctly (["+ str(StartReplica)+','+str(FinalReplica)+"] range of replicas computed) +\n "
          +"----------------------------------------------------------------------------------")