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
#PDFinUse="MSHT20"
#PDFinUse="CJ15"

## if this trigger is ON the LHC data will be fit only by shape
useNormalizedLHCdata=False
## Include ATLAS 7TeV?
useA7data=False
## Split the low-energy experiment <Upsilon and >Upsilon
splitUpsilon=True
## Use the reduced set of the data
useReducedSet=True
## Use penalty term for too different flavors
usePenaltyTerm=False

## total number of parameters in TMDPDF (7 and 12 case are defined)
numberOfParameters=12

#### Starting and final replica (included)
StartReplica=1
FinalReplica=10

## automatic name generation for the run
runName="model4.0_"+PDFinUse+"_EXPrep_"
if (not useA7data): runName+="noA7_"
if (splitUpsilon): runName+="spltUPS_"
if (useNormalizedLHCdata): runName+="norm_"
if (useReducedSet): runName+="reducedSet_"
if (usePenaltyTerm): runName+="penaltyTerm_"
if (runName[-1]=="_"): runName=runName[0:-1]

print(" RUN: "+runName)

#%%
#######################################
# Paths
#######################################
PathToHarpy="/home/vla18041/LinkData2/arTeMiDe_Repository/artemide-ForPDF/harpy"
#PathToHarpy="/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/harpy"
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
#Initialize artemide
#######################################
import harpy

print('Initialization with : \n'+PathToConstantsFile)

harpy.initialize(PathToConstantsFile)
if(numberOfParameters==7):
    initializationArray=[2.0340, 0.0299, 0.2512, 7.7572, 0.2512, 7.7572, 0.2512, 7.7572, 10000.]
    #initializationArray=[2., 0.0398333,0.184739, 6.22437, 588.193, 2.44327, -2.51106, 0.,  0.] #SV19
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

#### Check the point kinematics
def cutFunc0(p):    
    
    par=1.0

    if(p["xSec"]>0):
        #err=numpy.sqrt(sum([i**2 for i in p["uncorrErr"]]))/p["xSec"]
        err=10000.0 ## I would like to add all points into the full set
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
    if(useReducedSet and not(p["id"] in includePoints)):
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


#%%
if useNormalizedLHCdata:
    for s in setDY.sets:
        if s.name[0:4]=='LHCb':
            s.isNormalized=True
            
        if s.isNormalized:
            s.normalizationMethod='bestChi2'
#%%

#initialValues=(2., 0.0398333,0.184739, 6.22437, 588.193, 2.44327, -2.51106, 0.,  0.)#SV19

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
    initialValues=(2.000, 0.048, 0.063, 0.009, 0.414, 4.566, 0.049, 24.999, 0.001, 0.213, 0.002, 24.326, 0.036, 0.000) #model 4.0 CT18
if(PDFinUse=="MSHT20"):    
    initialValues=(2.000, 0.045, 0.184, 0.097, 0.446, 0.966, 0.022, 46.341, 0.002, 2.661, 0.010, 7.659, 65.137, 0.000) #model 4.0 MHST20
if(PDFinUse=="CJ15"):    
    initialValues=(2.000, 0.037, 0.147, 0.989, 0.269, 9.893, 0.241, 70.096, 0.024, 3.355, 0.001, 9.262, 130.386, 0.000) #model 4.0 CJ15
harpy.setNPparameters(list(initialValues))

DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)
#%%
#######################################
# Main chi2 formula
#######################################
totalN=setDY.numberOfPoints


#### The penalty term accounts for possible penalty effects. 
#### This one prevents the fNP be too different for different flavors
def PenaltyTerm(x):
    r1=0.1*max(x[2]/x[4],x[4]/x[2])
    r2=0.1*max(x[3]/x[5],x[5]/x[3])
    r3=0.1*max(x[6]/x[8],x[8]/x[6])
    r4=0.1*max(x[7]/x[9],x[9]/x[7])
    
    res=0
    if(r1>1): res+=0.1*(r1-1.)
    if(r2>1): res+=0.1*(r2-1.)
    if(r3>1): res+=0.1*(r3-1.)
    if(r4>1): res+=0.1*(r4-1.)
    
    return res
    

def chi_2(x):
    startT=time.time()
    harpy.setNPparameters(x)
    print('np set =',["{:8.3f}".format(i) for i in x], end =" ")    
    
    ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)
    
    if(usePenaltyTerm): ccDY2+=totalN*PenaltyTerm(x)
    
    cc=ccDY2/totalN
    endT=time.time()
    print(':->',cc,'       t=',endT-startT)
    return ccDY2

#%%
#######################################
# Setting-up Minuit
#######################################
from iminuit import Minuit

##### model SV19
# initialErrors=(0.1,       0.002,      0.05,    0.2,  0.05,   0.2,    10.0,   0.1,   0.1)
# searchLimits=((1.4,4.5), (0.0001,5.0),(0.0,5.0),(0.,50.0),(0.0,1000.),(0.,25.0),(-30.,30.),None, None)
# parametersToMinimize=(True,     False,    False,    False,    False,     False,  False, True, True)

##### model 1.0
#initialErrors=(0.1,       0.002,      0.05,    0.2,  0.05,   0.2,    10.0,   0.1,   0.1)
#searchLimits=((1.4,4.5), (0.0001,5.0),(0.0,2.0),(0.,25.0),(0.0,2.0),(0.,25.0),(0.,1000),None, None)
#parametersToMinimize=(True,     False,    False,    False,    False,     False,  False, True, True)

##### model 2.0
# initialErrors=(0.1,       0.002,      0.05,    0.2,  0.05,    0.2,  0.05,   0.2,    10.0)
# searchLimits=((1.4,4.5), (0.0001,5.0),(0.0,10.0),(0.,50.0),(0.0,10.0),(0.,50.0),(0.0,10.0),(0.,50.0),(0.,2500))
# parametersToMinimize=(True,     False,    False,    False,    False,     False,  False, False, False)

# ##### model 4.0
initialErrors=(0.1,       0.002,      0.05,    0.2,  0.05,   0.2, 0.05,    0.2,0.05,    0.2,   0.05,    0.2, 10., 0.1)
searchLimits=((1.4,4.5), (0.0001,5.0),(0.001,2.0),(0.001,100.0),(0.001,2.0),(0.001,100.0),(0.001,2.0),
              (0.001,100.0),(0.001,2.0),(0.001,100.0),(0.001,2.0),(0.001,100.0),(0.001,2500),(0.,100))
parametersToMinimize=(True,     False,    False,    False,    False,     False,  False, False, False, False, False, False, False,True)


#%%

m = Minuit(chi_2, initialValues)

m.errors=initialErrors
m.limits=searchLimits
m.fixed=parametersToMinimize
m.errordef=1

print(m.params)

m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1
# m.migrad()

# ## print parameters
# print(m.params)

# ## print chi^2 table
# harpy.setNPparameters(list(m.values))

# print('PENALTY = ',PenaltyTerm(list(m.values)))

# DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)

# DataProcessor.harpyInterface.PrintChi2Table(setDYfull,printDecomposedChi2=True)

# sys.exit()
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
        if(usePenaltyTerm): ccDY2+=totalN*PenaltyTerm(x)
        
        cc=(ccDY2)/totalN
        endT=time.time()
        print(':->',cc,'       t=',endT-startT)
        return ccDY2
    
    repDataDY=setDY.GenerateReplica()
    
    localM = Minuit(repchi_2, initialValues)
    
    localM.errors=initialErrors
    localM.limits=searchLimits
    localM.fixed=parametersToMinimize
    localM.errordef=1    
    localM.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
    localM.strategy=1

    localM.migrad()
    
    ### [chi^2, NP-parameters]
    return [localM.fval,list(localM.values)]


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
    print("Start computation of replica "+str(i) +"in ["+ str(StartReplica)+','+str(FinalReplica)+"]")
    
    ## got to pseudo-data and minimization
    repRes=MinForReplica()
    print(repRes)
    print("Minimization for replica "+str(i) +"in ["+ str(StartReplica)+','+str(FinalReplica)+"] finished.")    
    
    ## compute the chi2 for true data full
    mainDY, mainDY2 =DataProcessor.harpyInterface.ComputeChi2(setDYfull)    
    print("Central chi^2 for "+str(i) +"in ["+ str(StartReplica)+','+str(FinalReplica)+" computed. \n Saving to log >> "+replicaFile)
    
    ## save to file
    f=open(replicaFile,"a+")
    print('SAVING >>  ',f.name)
    ### [total chi^2(full data), total chi^2 (fit data), list of chi^2 for experiments(full data), 
    ############# number of PDF, list of NP-parameters],penalty term(if used) ]
    if(usePenaltyTerm): 
        f.write(str([mainDY,repRes[0],mainDY2,i,repRes[1],PenaltyTerm(repRes[1])])+"\n")
    else:
        f.write(str([mainDY,repRes[0],mainDY2,i,repRes[1]])+"\n")
    f.close()
    