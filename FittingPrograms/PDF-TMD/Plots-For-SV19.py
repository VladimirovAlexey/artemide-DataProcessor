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
import os
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"

ATMDE_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide/"
#%%
import sys

import time
import numpy
sys.path.append(ROOT_DIR)
import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet

MAINPATH=ROOT_DIR
#%%
#######################################
#Initialize artemide with a replica -2
#######################################
import harpy
path_to_constants=MAINPATH+"FittingPrograms/SV19/Constants-files/"
harpy.initialize(path_to_constants+"DY_n3lo/const-NNPDF31_n3lo")
#harpy.initialize(path_to_constants+"DY_nnlo/const-HERA20_NNLO")
#harpy.initialize(path_to_constants+"DY_n3lo/const-HERA20_n3lo")
#harpy.initialize(path_to_constants+"DY_nnlo/const-MMHT14_NNLO")
#harpy.initialize(path_to_constants+"DY_nnlo/const-CT14_NNLO")
#harpy.initialize(path_to_constants+"DY_nnlo/const-PDF4LHC_NNLO")
harpy.setNPparameters([2.0340, 0.0299, 0.2512, 7.7572,334.6108, 2.4543,-4.8203, 0.1000,  0.0000])
#harpy.setNPparameters_TMDR(-2)
#harpy.setNPparameters_uTMDPDF(-2)

#%%

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(ATMDE_DIR+"Models/SV19/Replicas/"+
                                                  "DY_n3lo/DY_NNPDF31_n3lo.rep")
                                                  # "Sivers20_model9case1(noDY-n3lo).rep")

rSet.SetReplica()

#%%
# ##################################################
# # Plot of TMD vs. b with errorband
# # at given x, f
# ##################################################

# xValue=0.05
# f=2

# bValues=[0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, \
# 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.1, 1.2, 1.3, 1.4, 1.5, \
# 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, \
# 3., 3.1, 3.2, 3.3, 3.4, 3.5, 3.6]

# tmd0=[]
    
# for i in range(rSet.numberOfReplicas):
#     rSet.SetReplica(i)
    
#     tmd0.append([harpy.get_uTMDPDF(xValue,b,1)[5+f] for b in bValues])

# print(list(numpy.mean(tmd0,axis=0)))    
# print(list(numpy.std(tmd0,axis=0)))

#%%
# ##################################################
# # Plot of TMD at fnP=1
# # with variation band
# ##################################################

# xValue=0.1
# f=1

# bValues=[0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, \
# 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.1, 1.2, 1.3, 1.4, 1.5, \
# 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, \
# 3., 3.1, 3.2, 3.3, 3.4, 3.5, 3.6]

# ## setting fNP=1
# harpy.setNPparameters([2.0340, 0.000, 0.0, 0,334.6108, 2.4543,0, 0.1000,  0.0000])
# tmd0=[harpy.get_uTMDPDF(xValue,b,1)[5+f] for b in bValues]

# ## scale variation UP
# harpy.varyScales(1., 1., 1., 2.)    
# tmd1=[harpy.get_uTMDPDF(xValue,b,1)[5+f] for b in bValues]

# ## scale variation DOWN
# harpy.varyScales(1., 1., 1., 0.5)    
# tmd2=[harpy.get_uTMDPDF(xValue,b,1)[5+f] for b in bValues]



# print(list(tmd0))
# print(list(tmd1))
# print(list(tmd2))

# harpy.varyScales(1., 1., 1., 1.)


#%%
# ##################################################
# # Plot of TMD vs. x with errorband
# # at given x, f
# ##################################################

# xValues=[
# 0.0001, 0.00011, 0.00013, 0.00014, 0.00016, 0.00018, 0.0002, 0.00022, \
# 0.00025, 0.00028, 0.00032, 0.00035, 0.0004, 0.00045, 0.0005, 0.00056, 0.00063, \
# 0.00071, 0.00079, 0.00089,\
# 0.001, 0.0011, 0.0013, 0.0014, 0.0016, 0.0018, 0.002, 0.0022, \
# 0.0025, 0.0028, 0.0032, 0.0035, 0.004, 0.0045, 0.005, 0.0056, 0.0063, \
# 0.0071, 0.0079, 0.0089, 0.01, 0.011, 0.013, 0.014, 0.016, 0.018, \
# 0.02, 0.022, 0.025, 0.028, 0.032, 0.035, 0.04, 0.045, 0.05, 0.056, \
# 0.063, 0.071, 0.079, 0.089, 0.1, 0.11, 0.13, 0.14, 0.16, 0.18, 0.2, \
# 0.22, 0.25, 0.28, 0.32, 0.35, 0.4, 0.45, 0.5, 0.56, 0.63, 0.71, 0.79, \
# 0.89, 0.99]
# f=2

# bValues=0.2

# tmd0=[]
    
# for i in range(rSet.numberOfReplicas):
#     rSet.SetReplica(i)
    
#     tmd0.append([x*harpy.get_uTMDPDF(x,bValues,1)[5+f] for x in xValues])

# print(list(numpy.mean(tmd0,axis=0)))    
# print(list(numpy.std(tmd0,axis=0)))

#%%
# ##################################
# ## Save replicas of NNPDFs
# ##################################
# import DataProcessor.SaveTMDGrid

# rSet.SetReplica()

# r=9

# for i in range(r*100,(r+1)*100):
#     harpy.setPDFreplica(i)
#     rSet.SetReplica()
#     path='/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD-grids/Grids_optimal/SV19_n3lo_PDFrep/SV19_n3lo_PDFrep'\
#         +'_'+'{:04d}'.format(i)+'.dat'
#     DataProcessor.SaveTMDGrid.SaveGrid_optimal(path)
    
#%%
##################################
## Save replicas of SV19
##################################
import DataProcessor.SaveTMDGrid

rSet.SetReplica()



for i in range(rSet.numberOfReplicas):    
    rSet.SetReplica(i)
    path='/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/TMD-grids/Grids_optimal/SV19_n3lo/SV19_n3lo'\
        +'_'+'{:04d}'.format(i)+'.dat'
    DataProcessor.SaveTMDGrid.SaveGrid_optimal(path)