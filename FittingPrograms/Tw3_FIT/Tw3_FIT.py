#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 16 13:54:42 2025

@author: alexey
"""

#######################################
# importing libraries
#######################################
import os
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"

SNOWFLAKE_DIR = "/data/WorkingFiles/Twist3/Snowflake2/PySnowflake/"


import sys
import numpy
sys.path.append(ROOT_DIR)
sys.path.append(SNOWFLAKE_DIR)

#%%

import SnowFlake

#%%

path_to_INI=ROOT_DIR+"FittingPrograms/Tw3_FIT/INI/TEST.ini"
SnowFlake.initialize(path_to_INI)

#%%
NP_par=numpy.zeros(18)+0.1
SnowFlake.setNPparameters(NP_par)
SnowFlake.UpdateTables(1.0, 25.0)

#%%
dd=SnowFlake.G2List([0.1,0.2,0.3,0.4],[5.,5.,5.,5.],[100,100,100,100])
print(dd)

#%%
dd=SnowFlake.D2List([5.,5.,5.,5.],[100,100,100,100])
print(dd)