#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 20:15:46 2019

Program that parse various SIDIS data files to ADP-frendly formal

@author: vla18041
"""
import sys
sys.path.append("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/")

import DataProcessor.Point
import DataProcessor.DataSet
import numpy

path_to_save="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/FittingPrograms/SV19/PredictionCOMPASS/"

## this trigger add normalization uncertanty due to DIS normalization for multiplicities
addDISnormalizationUncertainty=True

M_proton=0.938
m_pion=0.139
m_kaon=0.494
#%%

# 1) Genral kinematic cuts

# Q^2 > 1 (GeV/c)^2
# 0.2 < y < 0.9 
# W > 5 GeV/c^2
# 0.003 < X_Bj < 0.130
# theta_gamma < 60 mrad
# 0.1 < z < 0.85
# 0.1 < P_T / (GeV/c) < 2.0

# theta_gamma is the angle between the virtual photon and the muon beam in the laboratory frame:
# the cut is not really effective, because of the W and y cut


# 2) Binning

# The P_T^2 distributions (normalised to the first bin in P_T^2 (0.02,0.06 GeV^2/c^2) ) 
# have been measured up to 3.0 GeV^2/c^2 in the following multid bins:
# 	1. 4 bins in x, 2 bins in Q^2, 4 bins in z
# 	2. 4 bins in x, 4 bins in Q^2, 4 bins in z
# 	3. 4 bins in x, 2 bins in Q^2, 4 bins in z, 2 bins in W
# 	4. 4 bins in x, 4 bins in Q^2, 4 bins in z, 2 bins in W
# 	5. 4 bins in x, 2 bins in Q^2, 7 bins in z

# For 1. and 2., the bin limits are
#    x: 0.003, 0.013, 0.020, 0.055, 0.100
#    z: 0.20,  0.30,  0.40,  0.60,  0.80
#    for 1.   Q^2: 1.0, 3.0, 16.0 (GeV/c)^2
#    for 2.   Q^2: 1.0, 1.7, 3.0, 7.0, 16.0  (GeV/c)^2

# For 3. and 4. the x, z, and Q^2 binning is the same as 1. and 2., but binning
# the data in W> and < 12 GeV/c^2

# For 5. the Q^2 binning is the same as for 1. but it is
#    x:   0.003, 0.012, 0.020, 0.038, 0.130
#    z:   0.10, 0.20, 0.25, 0.30, 0.40, 0.55, 0.70, 0.85

case1_x=[[0.003,0.013],[0.013,0.020],[0.020,0.055],[0.055,0.10]]
case2_x=case1_x
case3_x=case1_x
case4_x=case1_x
case5_x=[[0.003,0.012],[0.012,0.020],[0.020,0.038],[0.038,0.13]]

case1_z=[[0.20,  0.30],[0.30, 0.40], [0.40, 0.60], [0.60,0.80] ]
case2_z=case1_z
case3_z=case1_z
case4_z=case1_z
case5_z=[[0.10, 0.20], [0.20, 0.25], [0.25, 0.30], [0.30, 0.40], [0.40, 0.55], [0.55, 0.70], [0.70, 0.85]]

case1_Q2=[[1.0,3.0],[3.0,16.0]]
case2_Q2=[[1.0,1.7],[1.7,3.0],[3.0,7.0],[7.0,16.0]]
case5_Q2=case1_Q2

W2cut=144.

#### pt-bins (assuming it is same as earlier)
pTbins=[[0.02, 0.04], [0.04, 0.06], [0.06, 0.08], [0.08, 0.1], [0.1, 
  0.12], [0.12, 0.14], [0.14, 0.17], [0.17, 0.2], [0.2, 0.23], [0.23, 
  0.27], [0.27, 0.3], [0.3, 0.35], [0.35, 0.4], [0.4, 0.46], [0.46, 
  0.52], [0.52, 0.6], [0.6, 0.68], [0.68, 0.76], [0.76, 0.87], [0.87, 
  1.], [1., 1.12], [1.12, 1.24], [1.24, 1.38], [1.38, 1.52], [1.52, 
  1.68], [1.68, 1.85], [1.85, 2.05], [2.05, 2.35], [2.35, 
  2.65], [2.65, 3.]]

#%%
###############################################################################
###########################COMPASS    H+ CASE 1 ###############################
###############################################################################
proc_current=[1,1,2103]
s_current=2*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.]

### case 1

print("H+ case 1")
DataCurrent=DataProcessor.DataSet.DataSet('compass.d.h+.case1',"SIDIS")
DataCurrent.comment="COMPASS isoscalar h+. prediction case1"
DataCurrent.reference="SV19"

n=0

currentCaseX=case1_x
currentCaseZ=case1_z
currentCaseQ2=case1_Q2

for iQ in range(len(currentCaseQ2)):
    for ix in range(len(currentCaseX)):
        for iz in range(len(currentCaseZ)):
            for iT in range(len(pTbins)):
                # makeup a point
                p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(n))
                #print DataCurrent.name+'.'+str(i)
                p["process"]=proc_current
                p["s"]=s_current    
                p["x"]=[currentCaseX[ix][0],currentCaseX[ix][1]]    
                p["Q"]=[numpy.sqrt(currentCaseQ2[iQ][0]),numpy.sqrt(currentCaseQ2[iQ][1])]    
                p["z"]=[currentCaseZ[iz][0],currentCaseZ[iz][1]]    
                p["pT"]=[numpy.sqrt(pTbins[iT][0]),numpy.sqrt(pTbins[iT][1])]    
                p["xSec"]=1.
                p["M_target"]=M_proton
                p["m_product"]=m_pion
                p["includeCuts"]=includeCuts
                p["cutParams"]=cutParameters
                p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["z"][1]-p["z"][0])#devide by bin size multiply by DIS xSec
                     
                p["uncorrErr"].append(0.1)
                #
                DataCurrent.AddPoint(p)
                
                n+=1

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################COMPASS    H+ CASE 2 ###############################
###############################################################################
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.]

### case 2

print("H+ case 1")
DataCurrent=DataProcessor.DataSet.DataSet('compass.d.h+.case2',"SIDIS")
DataCurrent.comment="COMPASS isoscalar h+. prediction case2"
DataCurrent.reference="SV19"

n=0

currentCaseX=case2_x
currentCaseZ=case2_z
currentCaseQ2=case2_Q2

for iQ in range(len(currentCaseQ2)):
    for ix in range(len(currentCaseX)):
        for iz in range(len(currentCaseZ)):
            for iT in range(len(pTbins)):
                # makeup a point
                p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(n))
                #print DataCurrent.name+'.'+str(i)
                p["process"]=proc_current
                p["s"]=s_current    
                p["x"]=[currentCaseX[ix][0],currentCaseX[ix][1]]    
                p["Q"]=[numpy.sqrt(currentCaseQ2[iQ][0]),numpy.sqrt(currentCaseQ2[iQ][1])]    
                p["z"]=[currentCaseZ[iz][0],currentCaseZ[iz][1]]    
                p["pT"]=[numpy.sqrt(pTbins[iT][0]),numpy.sqrt(pTbins[iT][1])]    
                p["xSec"]=1.
                p["M_target"]=M_proton
                p["m_product"]=m_pion
                p["includeCuts"]=includeCuts
                p["cutParams"]=cutParameters
                p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["z"][1]-p["z"][0])#devide by bin size multiply by DIS xSec
                     
                p["uncorrErr"].append(0.1)
                #
                DataCurrent.AddPoint(p)
                
                n+=1

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################COMPASS    H+ CASE 3 ###############################
###############################################################################
includeCuts=True

### case 3

print("H+ case 1")
DataCurrent=DataProcessor.DataSet.DataSet('compass.d.h+.case3',"SIDIS")
DataCurrent.comment="COMPASS isoscalar h+. prediction case3"
DataCurrent.reference="SV19"

n=0

currentCaseX=case1_x
currentCaseZ=case1_z
currentCaseQ2=case1_Q2

cutParameters=[0.1,0.9,25.,W2cut]

for iQ in range(len(currentCaseQ2)):
    for ix in range(len(currentCaseX)):
        for iz in range(len(currentCaseZ)):
            for iT in range(len(pTbins)):
                # makeup a point
                p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(n))
                #print DataCurrent.name+'.'+str(i)
                p["process"]=proc_current
                p["s"]=s_current    
                p["x"]=[currentCaseX[ix][0],currentCaseX[ix][1]]    
                p["Q"]=[numpy.sqrt(currentCaseQ2[iQ][0]),numpy.sqrt(currentCaseQ2[iQ][1])]    
                p["z"]=[currentCaseZ[iz][0],currentCaseZ[iz][1]]    
                p["pT"]=[numpy.sqrt(pTbins[iT][0]),numpy.sqrt(pTbins[iT][1])]    
                p["xSec"]=1.
                p["M_target"]=M_proton
                p["m_product"]=m_pion
                p["includeCuts"]=includeCuts
                p["cutParams"]=cutParameters
                p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["z"][1]-p["z"][0])#devide by bin size multiply by DIS xSec
                     
                p["uncorrErr"].append(0.1)
                #
                DataCurrent.AddPoint(p)
                
                n+=1
                
cutParameters=[0.1,0.9,W2cut,10000.]

for iQ in range(len(currentCaseQ2)):
    for ix in range(len(currentCaseX)):
        for iz in range(len(currentCaseZ)):
            for iT in range(len(pTbins)):
                # makeup a point
                p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(n))
                #print DataCurrent.name+'.'+str(i)
                p["process"]=proc_current
                p["s"]=s_current    
                p["x"]=[currentCaseX[ix][0],currentCaseX[ix][1]]    
                p["Q"]=[numpy.sqrt(currentCaseQ2[iQ][0]),numpy.sqrt(currentCaseQ2[iQ][1])]    
                p["z"]=[currentCaseZ[iz][0],currentCaseZ[iz][1]]    
                p["pT"]=[numpy.sqrt(pTbins[iT][0]),numpy.sqrt(pTbins[iT][1])]    
                p["xSec"]=1.
                p["M_target"]=M_proton
                p["m_product"]=m_pion
                p["includeCuts"]=includeCuts
                p["cutParams"]=cutParameters
                p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["z"][1]-p["z"][0])#devide by bin size multiply by DIS xSec
                     
                p["uncorrErr"].append(0.1)
                #
                DataCurrent.AddPoint(p)
                
                n+=1

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################COMPASS    H+ CASE 4 ###############################
###############################################################################
includeCuts=True

### case 4

print("H+ case 1")
DataCurrent=DataProcessor.DataSet.DataSet('compass.d.h+.case4',"SIDIS")
DataCurrent.comment="COMPASS isoscalar h+. prediction case4"
DataCurrent.reference="SV19"

n=0

currentCaseX=case2_x
currentCaseZ=case2_z
currentCaseQ2=case2_Q2

cutParameters=[0.1,0.9,25.,W2cut]

for iQ in range(len(currentCaseQ2)):
    for ix in range(len(currentCaseX)):
        for iz in range(len(currentCaseZ)):
            for iT in range(len(pTbins)):
                # makeup a point
                p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(n))
                #print DataCurrent.name+'.'+str(i)
                p["process"]=proc_current
                p["s"]=s_current    
                p["x"]=[currentCaseX[ix][0],currentCaseX[ix][1]]    
                p["Q"]=[numpy.sqrt(currentCaseQ2[iQ][0]),numpy.sqrt(currentCaseQ2[iQ][1])]    
                p["z"]=[currentCaseZ[iz][0],currentCaseZ[iz][1]]    
                p["pT"]=[numpy.sqrt(pTbins[iT][0]),numpy.sqrt(pTbins[iT][1])]    
                p["xSec"]=1.
                p["M_target"]=M_proton
                p["m_product"]=m_pion
                p["includeCuts"]=includeCuts
                p["cutParams"]=cutParameters
                p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["z"][1]-p["z"][0])#devide by bin size multiply by DIS xSec
                     
                p["uncorrErr"].append(0.1)
                #
                DataCurrent.AddPoint(p)
                
                n+=1
                
cutParameters=[0.1,0.9,W2cut,10000.]

for iQ in range(len(currentCaseQ2)):
    for ix in range(len(currentCaseX)):
        for iz in range(len(currentCaseZ)):
            for iT in range(len(pTbins)):
                # makeup a point
                p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(n))
                #print DataCurrent.name+'.'+str(i)
                p["process"]=proc_current
                p["s"]=s_current    
                p["x"]=[currentCaseX[ix][0],currentCaseX[ix][1]]    
                p["Q"]=[numpy.sqrt(currentCaseQ2[iQ][0]),numpy.sqrt(currentCaseQ2[iQ][1])]    
                p["z"]=[currentCaseZ[iz][0],currentCaseZ[iz][1]]    
                p["pT"]=[numpy.sqrt(pTbins[iT][0]),numpy.sqrt(pTbins[iT][1])]    
                p["xSec"]=1.
                p["M_target"]=M_proton
                p["m_product"]=m_pion
                p["includeCuts"]=includeCuts
                p["cutParams"]=cutParameters
                p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["z"][1]-p["z"][0])#devide by bin size multiply by DIS xSec
                     
                p["uncorrErr"].append(0.1)
                #
                DataCurrent.AddPoint(p)
                
                n+=1

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################COMPASS    H+ CASE 5 ###############################
###############################################################################
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.]

### case 2

print("H+ case 5")
DataCurrent=DataProcessor.DataSet.DataSet('compass.d.h+.case5',"SIDIS")
DataCurrent.comment="COMPASS isoscalar h+. prediction case5"
DataCurrent.reference="SV19"

n=0

currentCaseX=case5_x
currentCaseZ=case5_z
currentCaseQ2=case5_Q2

for iQ in range(len(currentCaseQ2)):
    for ix in range(len(currentCaseX)):
        for iz in range(len(currentCaseZ)):
            for iT in range(len(pTbins)):
                # makeup a point
                p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(n))
                #print DataCurrent.name+'.'+str(i)
                p["process"]=proc_current
                p["s"]=s_current    
                p["x"]=[currentCaseX[ix][0],currentCaseX[ix][1]]    
                p["Q"]=[numpy.sqrt(currentCaseQ2[iQ][0]),numpy.sqrt(currentCaseQ2[iQ][1])]    
                p["z"]=[currentCaseZ[iz][0],currentCaseZ[iz][1]]    
                p["pT"]=[numpy.sqrt(pTbins[iT][0]),numpy.sqrt(pTbins[iT][1])]    
                p["xSec"]=1.
                p["M_target"]=M_proton
                p["m_product"]=m_pion
                p["includeCuts"]=includeCuts
                p["cutParams"]=cutParameters
                p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["z"][1]-p["z"][0])#devide by bin size multiply by DIS xSec
                     
                p["uncorrErr"].append(0.1)
                #
                DataCurrent.AddPoint(p)
                
                n+=1

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################COMPASS    H- CASE 1 ###############################
###############################################################################
proc_current=[1,1,2113]
s_current=2*160.*0.938+(0.938)**2
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.]

### case 1

print("H- case 1")
DataCurrent=DataProcessor.DataSet.DataSet('compass.d.h-.case1',"SIDIS")
DataCurrent.comment="COMPASS isoscalar h-. prediction case1"
DataCurrent.reference="SV19"

n=0

currentCaseX=case1_x
currentCaseZ=case1_z
currentCaseQ2=case1_Q2

for iQ in range(len(currentCaseQ2)):
    for ix in range(len(currentCaseX)):
        for iz in range(len(currentCaseZ)):
            for iT in range(len(pTbins)):
                # makeup a point
                p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(n))
                #print DataCurrent.name+'.'+str(i)
                p["process"]=proc_current
                p["s"]=s_current    
                p["x"]=[currentCaseX[ix][0],currentCaseX[ix][1]]    
                p["Q"]=[numpy.sqrt(currentCaseQ2[iQ][0]),numpy.sqrt(currentCaseQ2[iQ][1])]    
                p["z"]=[currentCaseZ[iz][0],currentCaseZ[iz][1]]    
                p["pT"]=[numpy.sqrt(pTbins[iT][0]),numpy.sqrt(pTbins[iT][1])]    
                p["xSec"]=1.
                p["M_target"]=M_proton
                p["m_product"]=m_pion
                p["includeCuts"]=includeCuts
                p["cutParams"]=cutParameters
                p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["z"][1]-p["z"][0])#devide by bin size multiply by DIS xSec
                     
                p["uncorrErr"].append(0.1)
                #
                DataCurrent.AddPoint(p)
                
                n+=1

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################COMPASS    h- CASE 2 ###############################
###############################################################################
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.]

### case 2

print("h- case 1")
DataCurrent=DataProcessor.DataSet.DataSet('compass.d.h-.case2',"SIDIS")
DataCurrent.comment="COMPASS isoscalar h-. prediction case2"
DataCurrent.reference="SV19"

n=0

currentCaseX=case2_x
currentCaseZ=case2_z
currentCaseQ2=case2_Q2

for iQ in range(len(currentCaseQ2)):
    for ix in range(len(currentCaseX)):
        for iz in range(len(currentCaseZ)):
            for iT in range(len(pTbins)):
                # makeup a point
                p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(n))
                #print DataCurrent.name+'.'+str(i)
                p["process"]=proc_current
                p["s"]=s_current    
                p["x"]=[currentCaseX[ix][0],currentCaseX[ix][1]]    
                p["Q"]=[numpy.sqrt(currentCaseQ2[iQ][0]),numpy.sqrt(currentCaseQ2[iQ][1])]    
                p["z"]=[currentCaseZ[iz][0],currentCaseZ[iz][1]]    
                p["pT"]=[numpy.sqrt(pTbins[iT][0]),numpy.sqrt(pTbins[iT][1])]    
                p["xSec"]=1.
                p["M_target"]=M_proton
                p["m_product"]=m_pion
                p["includeCuts"]=includeCuts
                p["cutParams"]=cutParameters
                p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["z"][1]-p["z"][0])#devide by bin size multiply by DIS xSec
                     
                p["uncorrErr"].append(0.1)
                #
                DataCurrent.AddPoint(p)
                
                n+=1

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################COMPASS    h- CASE 3 ###############################
###############################################################################
includeCuts=True

### case 3

print("h- case 1")
DataCurrent=DataProcessor.DataSet.DataSet('compass.d.h-.case3',"SIDIS")
DataCurrent.comment="COMPASS isoscalar h-. prediction case3"
DataCurrent.reference="SV19"

n=0

currentCaseX=case1_x
currentCaseZ=case1_z
currentCaseQ2=case1_Q2

cutParameters=[0.1,0.9,25.,W2cut]

for iQ in range(len(currentCaseQ2)):
    for ix in range(len(currentCaseX)):
        for iz in range(len(currentCaseZ)):
            for iT in range(len(pTbins)):
                # makeup a point
                p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(n))
                #print DataCurrent.name+'.'+str(i)
                p["process"]=proc_current
                p["s"]=s_current    
                p["x"]=[currentCaseX[ix][0],currentCaseX[ix][1]]    
                p["Q"]=[numpy.sqrt(currentCaseQ2[iQ][0]),numpy.sqrt(currentCaseQ2[iQ][1])]    
                p["z"]=[currentCaseZ[iz][0],currentCaseZ[iz][1]]    
                p["pT"]=[numpy.sqrt(pTbins[iT][0]),numpy.sqrt(pTbins[iT][1])]    
                p["xSec"]=1.
                p["M_target"]=M_proton
                p["m_product"]=m_pion
                p["includeCuts"]=includeCuts
                p["cutParams"]=cutParameters
                p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["z"][1]-p["z"][0])#devide by bin size multiply by DIS xSec
                     
                p["uncorrErr"].append(0.1)
                #
                DataCurrent.AddPoint(p)
                
                n+=1
                
cutParameters=[0.1,0.9,W2cut,10000.]

for iQ in range(len(currentCaseQ2)):
    for ix in range(len(currentCaseX)):
        for iz in range(len(currentCaseZ)):
            for iT in range(len(pTbins)):
                # makeup a point
                p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(n))
                #print DataCurrent.name+'.'+str(i)
                p["process"]=proc_current
                p["s"]=s_current    
                p["x"]=[currentCaseX[ix][0],currentCaseX[ix][1]]    
                p["Q"]=[numpy.sqrt(currentCaseQ2[iQ][0]),numpy.sqrt(currentCaseQ2[iQ][1])]    
                p["z"]=[currentCaseZ[iz][0],currentCaseZ[iz][1]]    
                p["pT"]=[numpy.sqrt(pTbins[iT][0]),numpy.sqrt(pTbins[iT][1])]    
                p["xSec"]=1.
                p["M_target"]=M_proton
                p["m_product"]=m_pion
                p["includeCuts"]=includeCuts
                p["cutParams"]=cutParameters
                p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["z"][1]-p["z"][0])#devide by bin size multiply by DIS xSec
                     
                p["uncorrErr"].append(0.1)
                #
                DataCurrent.AddPoint(p)
                
                n+=1

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################COMPASS    h- CASE 4 ###############################
###############################################################################
includeCuts=True

### case 4

print("h- case 1")
DataCurrent=DataProcessor.DataSet.DataSet('compass.d.h-.case4',"SIDIS")
DataCurrent.comment="COMPASS isoscalar h-. prediction case4"
DataCurrent.reference="SV19"

n=0

currentCaseX=case2_x
currentCaseZ=case2_z
currentCaseQ2=case2_Q2

cutParameters=[0.1,0.9,25.,W2cut]

for iQ in range(len(currentCaseQ2)):
    for ix in range(len(currentCaseX)):
        for iz in range(len(currentCaseZ)):
            for iT in range(len(pTbins)):
                # makeup a point
                p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(n))
                #print DataCurrent.name+'.'+str(i)
                p["process"]=proc_current
                p["s"]=s_current    
                p["x"]=[currentCaseX[ix][0],currentCaseX[ix][1]]    
                p["Q"]=[numpy.sqrt(currentCaseQ2[iQ][0]),numpy.sqrt(currentCaseQ2[iQ][1])]    
                p["z"]=[currentCaseZ[iz][0],currentCaseZ[iz][1]]    
                p["pT"]=[numpy.sqrt(pTbins[iT][0]),numpy.sqrt(pTbins[iT][1])]    
                p["xSec"]=1.
                p["M_target"]=M_proton
                p["m_product"]=m_pion
                p["includeCuts"]=includeCuts
                p["cutParams"]=cutParameters
                p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["z"][1]-p["z"][0])#devide by bin size multiply by DIS xSec
                     
                p["uncorrErr"].append(0.1)
                #
                DataCurrent.AddPoint(p)
                
                n+=1
                
cutParameters=[0.1,0.9,W2cut,10000.]

for iQ in range(len(currentCaseQ2)):
    for ix in range(len(currentCaseX)):
        for iz in range(len(currentCaseZ)):
            for iT in range(len(pTbins)):
                # makeup a point
                p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(n))
                #print DataCurrent.name+'.'+str(i)
                p["process"]=proc_current
                p["s"]=s_current    
                p["x"]=[currentCaseX[ix][0],currentCaseX[ix][1]]    
                p["Q"]=[numpy.sqrt(currentCaseQ2[iQ][0]),numpy.sqrt(currentCaseQ2[iQ][1])]    
                p["z"]=[currentCaseZ[iz][0],currentCaseZ[iz][1]]    
                p["pT"]=[numpy.sqrt(pTbins[iT][0]),numpy.sqrt(pTbins[iT][1])]    
                p["xSec"]=1.
                p["M_target"]=M_proton
                p["m_product"]=m_pion
                p["includeCuts"]=includeCuts
                p["cutParams"]=cutParameters
                p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["z"][1]-p["z"][0])#devide by bin size multiply by DIS xSec
                     
                p["uncorrErr"].append(0.1)
                #
                DataCurrent.AddPoint(p)
                
                n+=1

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")

#%%
###############################################################################
###########################COMPASS    h- CASE 5 ###############################
###############################################################################
includeCuts=True
cutParameters=[0.1,0.9,25.,10000.]

### case 2

print("h- case 5")
DataCurrent=DataProcessor.DataSet.DataSet('compass.d.h-.case5',"SIDIS")
DataCurrent.comment="COMPASS isoscalar h-. prediction case5"
DataCurrent.reference="SV19"

n=0

currentCaseX=case5_x
currentCaseZ=case5_z
currentCaseQ2=case5_Q2

for iQ in range(len(currentCaseQ2)):
    for ix in range(len(currentCaseX)):
        for iz in range(len(currentCaseZ)):
            for iT in range(len(pTbins)):
                # makeup a point
                p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(n))
                #print DataCurrent.name+'.'+str(i)
                p["process"]=proc_current
                p["s"]=s_current    
                p["x"]=[currentCaseX[ix][0],currentCaseX[ix][1]]    
                p["Q"]=[numpy.sqrt(currentCaseQ2[iQ][0]),numpy.sqrt(currentCaseQ2[iQ][1])]    
                p["z"]=[currentCaseZ[iz][0],currentCaseZ[iz][1]]    
                p["pT"]=[numpy.sqrt(pTbins[iT][0]),numpy.sqrt(pTbins[iT][1])]    
                p["xSec"]=1.
                p["M_target"]=M_proton
                p["m_product"]=m_pion
                p["includeCuts"]=includeCuts
                p["cutParams"]=cutParameters
                p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["z"][1]-p["z"][0])#devide by bin size multiply by DIS xSec
                     
                p["uncorrErr"].append(0.1)
                #
                DataCurrent.AddPoint(p)
                
                n+=1

print("Done.  ")

DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
