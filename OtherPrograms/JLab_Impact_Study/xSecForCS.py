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
path_to_constants=MAINPATH+"FittingPrograms/MVZ22/ConstantsFiles/"+"const-DY+SIDIS_NNPDF31+DSS_nnlo_m=0"


harpy.initialize(path_to_constants)

harpy.setNPparameters_TMDR([1.93, 0.0434])
harpy.setNPparameters_uTMDPDF([0.253434, 9.04351, 346.999, 2.47992, -5.69988, 0.1, 0.])
harpy.setNPparameters_uTMDFF([0.264,0.479,0.459,0.539]) 

originalSV19TMDR=[1.93, 0.0434]
originalSV19FF=[0.264,0.479,0.459,0.539]

#%%
def CreateDataSet(Qbin,xBin,zBin,ptList):
    DataN=DataProcessor.DataSet.DataSet("temp","SIDIS")
    DataN.comment="N"
    DataN.reference="Harut"
    
    proc_current=[1,1,2001]
    #s_current=2*22.*0.938+(0.938)**2
    s_current=2*42.*0.938+(0.938)**2 ### NOT JLAB!
    includeCuts=False
    cutParameters=[0.1,0.85,10.,10000.]
    
    for i in range(len(ptList)):
        # makeup a point
        p=DataProcessor.Point.CreateSIDISPoint(DataN.name+'.'+str(i))
        #print DataCurrent.name+'.'+str(i)
        p["process"]=proc_current
        p["s"]=s_current
        p["pT"]=ptList[i]        
        p["Q"]=Qbin
        p["x"]=xBin
        p["z"]=zBin
        ## cross-seciton is unknown
        p["xSec"]=0.1
        p["M_target"]=0.938
        p["M_product"]=0.139
        p["includeCuts"]=includeCuts
        p["cutParams"]=cutParameters
        #devide by (x,z,pt^2,Q^2) bin size
        #p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["Q"][1]**2-p["Q"][0]**2)/(p["z"][1]-p["z"][0])/(p["x"][1]-p["x"][0])
        p["thFactor"]=1.
        
        ### to get the subtracted cross-section
        from DataProcessor.Point import FinalizePoint
        FinalizePoint(p)
        y=p["<Q>"]**2/p["<x>"]/(p["s"]-0*p["M_target"]**2)
        gamma=0.*2.*p["<x>"]*p["M_target"]/p["<Q>"]
        eps=(1.-y-(gamma*y)**2/4.)/(1.-y+y**2/2.+(gamma*y)**2/4)
        p["thFactor"]=(1.-eps)/(y**2)/(p["pT"][1]**2-p["pT"][0]**2)/(p["Q"][1]**2-p["Q"][0]**2)/(p["z"][1]-p["z"][0])/(p["x"][1]-p["x"][0])
        
        DataN.AddPoint(p)
    
    print("Done.  ")
    
    return DataN
#%%
DUMPPATH="/home/vla18041/WorkingFiles/TMD/Fit_Notes/CSfromJLab/Data/ArtemideTest/"

#%%
ss=CreateDataSet([4.24264, 4.3589],[0.56,0.58],[0.48,0.52],[[0.02*i,0.02*(i+1)] for i in range(300)])
res=DataProcessor.harpyInterface.ComputeXSec(ss)

dd={
    'qT':[[p["pT"][0]/p["<z>"],p["pT"][1]/p["<z>"]] for p in ss.points],
    'Q':ss.points[0]["Q"],
    'X':res
    }
from json import dump
with open(DUMPPATH+"case11", 'w') as fp:
    dump(dd, fp)
    

ss=CreateDataSet([4., 4.12311],[0.56,0.58],[0.48,0.52],[[0.02*i,0.02*(i+1)] for i in range(300)])
res=DataProcessor.harpyInterface.ComputeXSec(ss)

dd={
    'qT':[[p["pT"][0]/p["<z>"],p["pT"][1]/p["<z>"]] for p in ss.points],
    'Q':ss.points[0]["Q"],
    'X':res
    }
from json import dump
with open(DUMPPATH+"case12", 'w') as fp:
    dump(dd, fp)
    
ss=CreateDataSet([3., 3.12311],[0.56,0.58],[0.48,0.52],[[0.02*i,0.02*(i+1)] for i in range(300)])
res=DataProcessor.harpyInterface.ComputeXSec(ss)

dd={
    'qT':[[p["pT"][0]/p["<z>"],p["pT"][1]/p["<z>"]] for p in ss.points],
    'Q':ss.points[0]["Q"],
    'X':res
    }
from json import dump
with open(DUMPPATH+"case13", 'w') as fp:
    dump(dd, fp)


ss=CreateDataSet([2.82, 3.],[0.44,0.46],[0.48,0.52],[[0.02*i,0.02*(i+1)] for i in range(300)])
res=DataProcessor.harpyInterface.ComputeXSec(ss)

dd={
    'qT':[[p["pT"][0]/p["<z>"],p["pT"][1]/p["<z>"]] for p in ss.points],
    'Q':ss.points[0]["Q"],
    'X':res
    }
from json import dump
with open(DUMPPATH+"case21", 'w') as fp:
    dump(dd, fp)


ss=CreateDataSet([2.44, 2.64],[0.44,0.46],[0.48,0.52],[[0.02*i,0.02*(i+1)] for i in range(300)])
res=DataProcessor.harpyInterface.ComputeXSec(ss)


dd={
    'qT':[[p["pT"][0]/p["<z>"],p["pT"][1]/p["<z>"]] for p in ss.points],
    'Q':ss.points[0]["Q"],
    'X':res
    }
from json import dump
with open(DUMPPATH+"case22", 'w') as fp:
    dump(dd, fp)
    
ss=CreateDataSet([1.82, 1.92],[0.44,0.46],[0.48,0.52],[[0.02*i,0.02*(i+1)] for i in range(300)])
res=DataProcessor.harpyInterface.ComputeXSec(ss)


dd={
    'qT':[[p["pT"][0]/p["<z>"],p["pT"][1]/p["<z>"]] for p in ss.points],
    'Q':ss.points[0]["Q"],
    'X':res
    }
from json import dump
with open(DUMPPATH+"case23", 'w') as fp:
    dump(dd, fp)


ss=CreateDataSet([4.24264, 4.3589],[0.2,0.3],[0.28,0.32],[[0.02*i,0.02*(i+1)] for i in range(300)])
res=DataProcessor.harpyInterface.ComputeXSec(ss)

dd={
    'qT':[[p["pT"][0]/p["<z>"],p["pT"][1]/p["<z>"]] for p in ss.points],
    'Q':ss.points[0]["Q"],
    'X':res
    }
from json import dump
with open(DUMPPATH+"case31", 'w') as fp:
    dump(dd, fp)
    

ss=CreateDataSet([4., 4.12311],[0.2,0.3],[0.28,0.32],[[0.02*i,0.02*(i+1)] for i in range(300)])
res=DataProcessor.harpyInterface.ComputeXSec(ss)

dd={
    'qT':[[p["pT"][0]/p["<z>"],p["pT"][1]/p["<z>"]] for p in ss.points],
    'Q':ss.points[0]["Q"],
    'X':res
    }
from json import dump
with open(DUMPPATH+"case32", 'w') as fp:
    dump(dd, fp)
    
ss=CreateDataSet([3., 3.12311],[0.2,0.3],[0.28,0.32],[[0.02*i,0.02*(i+1)] for i in range(300)])
res=DataProcessor.harpyInterface.ComputeXSec(ss)

dd={
    'qT':[[p["pT"][0]/p["<z>"],p["pT"][1]/p["<z>"]] for p in ss.points],
    'Q':ss.points[0]["Q"],
    'X':res
    }
from json import dump
with open(DUMPPATH+"case33", 'w') as fp:
    dump(dd, fp)
    
#%%
from json import dump
for j in range(7):
    ss=CreateDataSet([2*j*0.3+1.5, (2*j+1)*0.3+1.5],[0.44,0.46],[0.48,0.52],[[0.02*i,0.02*(i+1)] for i in range(300)])
    res=DataProcessor.harpyInterface.ComputeXSec(ss)
    
    dd={
        'qT':[[p["pT"][0]/p["<z>"],p["pT"][1]/p["<z>"]] for p in ss.points],
        'Q':ss.points[0]["Q"],
        'X':res
        }
    with open(DUMPPATH+"caseT"+str(j), 'w') as fp:
        dump(dd, fp)
        
#%%

def CreateDataSet(Qbinlist,xBin,zBin,ptbin):
    DataN=DataProcessor.DataSet.DataSet("temp","SIDIS")
    DataN.comment="N"
    DataN.reference="Harut"
    
    proc_current=[1,1,2001]
    #s_current=2*22.*0.938+(0.938)**2
    s_current=2*42.*0.938+(0.938)**2 ### NOT JLAB!
    includeCuts=False
    cutParameters=[0.1,0.85,10.,10000.]
    
    for i in range(len(Qbinlist)):
        # makeup a point
        p=DataProcessor.Point.CreateSIDISPoint(DataN.name+'.'+str(i))
        #print DataCurrent.name+'.'+str(i)
        p["process"]=proc_current
        p["s"]=s_current
        p["pT"]=ptbin
        p["Q"]=Qbinlist[i]
        p["x"]=xBin
        p["z"]=zBin
        ## cross-seciton is unknown
        p["xSec"]=0.1
        p["M_target"]=0.938
        p["M_product"]=0.139
        p["includeCuts"]=includeCuts
        p["cutParams"]=cutParameters
        #devide by (x,z,pt^2,Q^2) bin size
        #p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["Q"][1]**2-p["Q"][0]**2)/(p["z"][1]-p["z"][0])/(p["x"][1]-p["x"][0])
        p["thFactor"]=1.
        
        ### to get the subtracted cross-section
        from DataProcessor.Point import FinalizePoint
        FinalizePoint(p)
        y=p["<Q>"]**2/p["<x>"]/(p["s"]-0*p["M_target"]**2)
        gamma=0.*2.*p["<x>"]*p["M_target"]/p["<Q>"]
        eps=(1.-y-(gamma*y)**2/4.)/(1.-y+y**2/2.+(gamma*y)**2/4)
        p["thFactor"]=(1.-eps)/(y**2)/(p["pT"][1]**2-p["pT"][0]**2)/(p["Q"][1]**2-p["Q"][0]**2)/(p["z"][1]-p["z"][0])/(p["x"][1]-p["x"][0])
        
        DataN.AddPoint(p)
    
    print("Done.  ")
    
    return DataN

QQQQ=[[1.4142135623730951, 1.5811388300841898],
 [1.5811388300841898, 1.7320508075688772],
 [1.7320508075688772, 1.8708286933869707],
 [1.8708286933869707, 2.0],
 [2.0, 2.1213203435596424],
 [2.1213203435596424, 2.23606797749979],
 [2.23606797749979, 2.345207879911715],
 [2.345207879911715, 2.449489742783178],
 [2.449489742783178, 2.5495097567963922],
 [2.5495097567963922, 2.6457513110645907],
 [2.6457513110645907, 2.7386127875258306],
 [2.7386127875258306, 2.8284271247461903],
 [2.8284271247461903, 2.9154759474226504],
 [2.9154759474226504, 3.0],
 [3.0, 3.082207001484488],
 [3.082207001484488, 3.1622776601683795],
 [3.1622776601683795, 3.24037034920393],
 [3.24037034920393, 3.3166247903554],
 [3.3166247903554, 3.391164991562634],
 [3.391164991562634, 3.4641016151377544],
 [3.4641016151377544, 3.5355339059327378],
 [3.5355339059327378, 3.605551275463989],
 [3.605551275463989, 3.6742346141747673],
 [3.6742346141747673, 3.7416573867739413],
 [3.7416573867739413, 3.8078865529319543],
 [3.8078865529319543, 3.872983346207417],
 [3.872983346207417, 3.9370039370059056],
 [3.9370039370059056, 4.0],
 [4.0, 4.06201920231798],
 [4.06201920231798, 4.123105625617661],
 [4.123105625617661, 4.183300132670378],
 [4.183300132670378, 4.242640687119285],
 [4.242640687119285, 4.301162633521313],
 [4.301162633521313, 4.358898943540674],
 [4.358898943540674, 4.415880433163924],
 [4.415880433163924, 4.47213595499958],
 [4.47213595499958, 4.527692569068709],
 [4.527692569068709, 4.58257569495584],
 [4.58257569495584, 4.636809247747852],
 [4.636809247747852, 4.69041575982343]]

ss=CreateDataSet(QQQQ,[0.68, 0.72],[0.5,0.6],[0.4,0.6])
res=DataProcessor.harpyInterface.ComputeXSec(ss)