#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 15:47:01 2023

@author: alexey
"""

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

replicaFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/replicaDY.dat"
logFile =os.path.realpath(os.path.join(os.path.dirname(__file__))) +"/REPLICAS/LOG.log"

import sys
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)


#%%
import DataProcessor.harpyInterface

MAINPATH=ROOT_DIR 
#%%
#######################################
#Initialize artemide
#######################################
import harpy

#path_to_constants=MAINPATH+"FittingPrograms/MVZ22/ConstantsFiles/DYonly/PDFb_DYonly_N3LL+N4LO"
path_to_constants=MAINPATH+"FittingPrograms/MVZ22/ConstantsFiles/DYonly/ForScale_N4LL"

harpy.initialize(path_to_constants)

#harpy.setNPparameters_TMDR([1.48571, 0.0461374, 0.0706658,1.])

#harpy.setNPparameters_uTMDPDF(initializationArray)

#%%
harpy.setNPparameters_TMDR([1.5, 0.05, 0.07,1.])
bValues=[0.01*(n+1) for n in range(100)]
rr=[harpy.get_DNP(b, 2.) for b in bValues]
print(" ".join(["{"+str(b)+","+str(harpy.get_DNP(b, 2.))+"}," for b in bValues]))
