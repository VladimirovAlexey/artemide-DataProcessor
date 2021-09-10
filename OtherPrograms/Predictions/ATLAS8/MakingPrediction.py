#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 09:11:39 2021

@author: vla18041
"""


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
harpy.initialize(path_to_constants+"DY_nnlo/const-NNPDF31_NNLO")

harpy.setNPparameters([2.0340, 0.0299, 0.2512, 7.7572,334.6108, 2.4543,-4.8203, 0.1000,  0.0000])


#%%
### read the list of files and return the list of DataSets
import DataProcessor.DataSet
    
data=DataProcessor.DataSet.LoadCSV(
    "/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/OtherPrograms/Predictions/ATLAS8/A8-predict.csv")

############## I drop all points after ~0.25, last two bins are not trastful
def cutFunc(p):
    return p["<qT>"]<25.5 , p

setDY=data.CutData(cutFunc) 

#%%
### loading parameters of SV19

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/Replicas/"+
                                                  "DY_nnlo/DY_NNPDF31_nnlo.rep")
                                                  # "Sivers20_model9case1(noDY-n3lo).rep")

rSet.SetReplica()

#%%
## computing uncertanty due to SV19
import numpy

central=DataProcessor.harpyInterface.ComputeXSec(setDY)

dist=[]
for i in range(rSet.numberOfReplicas):
    rSet.SetReplica(i)
    dist.append(DataProcessor.harpyInterface.ComputeXSec(setDY))
    
std=numpy.std(dist,axis=0)

rSet.SetReplica()

#%%
## computing uncertanty due to scale variation

distS=[]

harpy.varyScales(1., 0.5, 1., 1.)
distS.append(DataProcessor.harpyInterface.ComputeXSec(setDY))

harpy.varyScales(1., 2., 1., 1.)
distS.append(DataProcessor.harpyInterface.ComputeXSec(setDY))

harpy.varyScales(1., 1., 1., 0.5)
distS.append(DataProcessor.harpyInterface.ComputeXSec(setDY))

harpy.varyScales(1., 1., 1., 2.)
distS.append(DataProcessor.harpyInterface.ComputeXSec(setDY))

harpy.varyScales(1., 1., 1., 1.)

#%%

scaleV=[max([abs(central[i]-distS[0][i]),abs(central[i]-distS[1][i]),abs(central[i]-distS[2][i]),abs(central[i]-distS[3][i])])
        for i in range(len(central))]

with open("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/OtherPrograms/Predictions/ATLAS8/A8-prediction.dat","w") as file:
    file.write("\
# Prediction for ATLAS 8TeV measurement of Z-boson, dSigma/dpT without cuts, by SV19 TMD extraction [1912.06532]. \
 Order = NNLL'. NP-unc. is uncertanty due to NP parameters as extracted in SV19. Scale variation is the scale variation band \
 (absolute symetric deviation from central with *2/2 variations over a box). Author: A.Vladimirov \n")
    file.write("# |ymin|, |ymax|, ptMin, ptMax, xSec, +-NP-unertainty, +-ScaleVariation")
    
    for i in range(len(central)):
        file.write(
            "{:g}".format(setDY.points[i]["y"][0])+", "\
            "{:g}".format(setDY.points[i]["y"][1])+", "\
            "{:g}".format(setDY.points[i]["qT"][0])+", "\
            "{:g}".format(setDY.points[i]["qT"][1])+", "\
            "{:g}".format(central[i])+", "\
            "{:g}".format(std[i])+", "\
            "{:g}".format(scaleV[i])+"\n ")