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
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/"
DATAP_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"

SNOWFLAKE_DIR = ROOT_DIR+"artemide/harpy/"
MODEL_DIR = ROOT_DIR+"artemide/Models/ART25/Replica-files/"


import sys
import numpy
if('/data/arTeMiDe_Repository/artemide/harpy' in sys.path):
    sys.path.remove('/data/arTeMiDe_Repository/artemide/harpy')
sys.path.append(DATAP_DIR)
sys.path.append(SNOWFLAKE_DIR)


#%%
import DataProcessor.harpyInterface
import DataProcessor.snowInterface_N2
import DataProcessor.DataMultiSet
import harpy

#%%
#######################################
#Initialize snowflake
#######################################
path_to_INI=DATAP_DIR+"FittingPrograms/Tw3_FIT/INI/snowflake_forPLOT.ini"
#path_to_INI=DATAP_DIR+"FittingPrograms/Tw3_FIT/INI/TEST.ini"
#path_to_INI=DATAP_DIR+"FittingPrograms/Tw3_FIT/INI/TEST_16x8.ini"
harpy.initialize_snowflake(path_to_INI)

NP_par=numpy.zeros(18)+0.2
harpy.setNPparameters_tw3(NP_par)
harpy.UpdateTables(1.0, 105.0)

#%%
#######################################
#Initialize artemide
#######################################

import DataProcessor.ArtemideReplicaSet

path_to_constants=DATAP_DIR+"FittingPrograms/Tw3_FIT/INI/TMD+tw3.atmde"


harpy.initialize(path_to_constants)

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(MODEL_DIR+"ART25_main.rep")
    
rSet.SetReplica(0)

#%%
# #### PLOT d2
# Qplot=0.1*numpy.array(range(40))+1
# dataP=harpy.D2List(Qplot, Qplot*0+100)
# dataN=harpy.D2List(Qplot, Qplot*0+101)
# dataD=harpy.D2List(Qplot, Qplot*0+102)

#%%
def X1X2(r,phi):
    if(0<=phi<1):
        return r*(1-phi),r*phi
    elif(1<=phi<2):
        return r*(1-phi),r
    elif(2<=phi<3):
        return -r,r*(3-phi)
    elif(3<=phi<4):
        return r*(phi-4),r*(3-phi)
    elif(4<=phi<5):
        return r*(phi-4),-r
    elif(5<=phi<6):
        return r,r*(phi-6)
    else:
        raise("ERRROR")
#%%
repLIST=[]
#with open("/data/WorkingFiles/TMD/Fit_Notes/Twist3_25/rreps4.csv","r") as file:
with open("/data/WorkingFiles/TMD/Fit_Notes/Twist3_25/rreps6.csv","r") as file:
    lines=file.readlines()
    for line in lines:
        repLIST.append([float(p) for p in line.split(",")])       


#%%
#path_to_save="/data/WorkingFiles/TMD/Fit_Notes/Twist3_25/data/Tw3PDF_4GeV_n4/"
path_to_save="/data/WorkingFiles/TMD/Fit_Notes/Twist3_25/data/Tw3PDF_4GeV_n6/"
rList=[0.01*i+0.01 for i in range(9)]+[0.05*i+0.1 for i in range(19)]
phiList=[i/8 for i in range(48)]
toSave=[]
Q=4.

for n in range(len(repLIST)):
    print("---->",n)
    harpy.setNPparameters_tw3(repLIST[n])
    harpy.UpdateTables(1.0, 6.0)

    # u-quark
    
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, 2)]]
    
    with open(path_to_save+"T_u."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
    # u-quark Delta
    
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, -2)]]
    
    with open(path_to_save+"deltaT_u."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
    
    # d-quark
            
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, 1)]]
    
    with open(path_to_save+"T_d."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
            
    # d-quark De;ta
            
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, -1)]]
    
    with open(path_to_save+"deltaT_d."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
    
    # s-quark
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, 3)]]
    
    with open(path_to_save+"T_s."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
            
    # s-quark Delta
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, -3)]]
    
    with open(path_to_save+"deltaT_s."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
            
    # c-quark
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, 4)]]
    
    with open(path_to_save+"T_c."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
            
    # c-quark Delta
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, -4)]]
    
    with open(path_to_save+"deltaT_c."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
    
    # g-quark
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, 0)]]
    
    with open(path_to_save+"T_g+."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
            
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, -10)]]
    
    with open(path_to_save+"T_g-."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')

#%%
#path_to_save="/data/WorkingFiles/TMD/Fit_Notes/Twist3_25/data/Tw3PDF_4GeV_n4/"
path_to_save="/data/WorkingFiles/TMD/Fit_Notes/Twist3_25/data/Tw3PDF_10GeV_n6/"
rList=[0.01*i+0.01 for i in range(9)]+[0.05*i+0.1 for i in range(19)]
phiList=[i/8 for i in range(48)]
toSave=[]
Q=10.

for n in range(len(repLIST)):
    print("---->",n)
    harpy.setNPparameters_tw3(repLIST[n])
    harpy.UpdateTables(1.0, 10.0)

    # u-quark
    
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, 2)]]
    
    with open(path_to_save+"T_u."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
    # u-quark Delta
    
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, -2)]]
    
    with open(path_to_save+"deltaT_u."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
    
    # d-quark
            
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, 1)]]
    
    with open(path_to_save+"T_d."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
            
    # d-quark De;ta
            
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, -1)]]
    
    with open(path_to_save+"deltaT_d."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
    
    # s-quark
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, 3)]]
    
    with open(path_to_save+"T_s."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
            
    # s-quark Delta
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, -3)]]
    
    with open(path_to_save+"deltaT_s."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
            
    # c-quark
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, 4)]]
    
    with open(path_to_save+"T_c."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
            
    # c-quark Delta
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, -4)]]
    
    with open(path_to_save+"deltaT_c."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
    
    # g-quark
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, 0)]]
    
    with open(path_to_save+"T_g+."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
            
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, -10)]]
    
    with open(path_to_save+"T_g-."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')

#%%
#path_to_save="/data/WorkingFiles/TMD/Fit_Notes/Twist3_25/data/Tw3PDF_4GeV_n4/"
path_to_save="/data/WorkingFiles/TMD/Fit_Notes/Twist3_25/data/Tw3PDF_100GeV_n6/"
rList=[0.01*i+0.01 for i in range(9)]+[0.05*i+0.1 for i in range(19)]
phiList=[i/8 for i in range(48)]
toSave=[]
Q=100.

for n in range(len(repLIST)):
    print("---->",n)
    harpy.setNPparameters_tw3(repLIST[n])
    harpy.UpdateTables(1.0, 100.0)

    # u-quark
    
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, 2)]]
    
    with open(path_to_save+"T_u."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
    # u-quark Delta
    
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, -2)]]
    
    with open(path_to_save+"deltaT_u."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
    
    # d-quark
            
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, 1)]]
    
    with open(path_to_save+"T_d."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
            
    # d-quark De;ta
            
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, -1)]]
    
    with open(path_to_save+"deltaT_d."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
    
    # s-quark
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, 3)]]
    
    with open(path_to_save+"T_s."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
            
    # s-quark Delta
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, -3)]]
    
    with open(path_to_save+"deltaT_s."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
            
    # c-quark
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, 4)]]
    
    with open(path_to_save+"T_c."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
            
    # c-quark Delta
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, -4)]]
    
    with open(path_to_save+"deltaT_c."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
    
    # g-quark
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, 0)]]
    
    with open(path_to_save+"T_g+."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
            
    toSave=[]
    
    for r in rList:       
        for phi in phiList:
            x1,x2=X1X2(r,phi)
            if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
            if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
            toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, -10)]]
    
    with open(path_to_save+"T_g-."+str(n).zfill(4), 'w') as file:
        for s in toSave:
            file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')

#%%
# path_to_save="/data/WorkingFiles/TMD/Fit_Notes/Twist3_25/data/Tw3PDF_4GeV/"
# rList=[0.01*i+0.01 for i in range(9)]+[0.05*i+0.1 for i in range(19)]
# phiList=[i/8 for i in range(48)]
# toSave=[]
# Q=4.

# for n in range(len(repLIST)):
#     print("---->",n)
#     harpy.setNPparameters_tw3(repLIST[n])
#     harpy.UpdateTables(1.0, 6.0)

#     # NS-quark
    
#     toSave=[]
    
#     for r in rList:       
#         for phi in phiList:
#             x1,x2=X1X2(r,phi)
#             if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
#             if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
#             toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, 2)-harpy.get_PDF_tw3(x1, x2, Q, 1)]]
    
#     with open(path_to_save+"T_NS1."+str(n).zfill(4), 'w') as file:
#         for s in toSave:
#             file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
            
#     # NS-quark
    
#     toSave=[]
    
#     for r in rList:       
#         for phi in phiList:
#             x1,x2=X1X2(r,phi)
#             if(numpy.abs(x1)<0.01 and numpy.abs(x2)<0.01): continue
#             if(numpy.abs(x1)>1. or numpy.abs(x2)>1. or numpy.abs(x1+x2)>1.): continue
        
#             toSave+=[[x1,x2,-x1-x2,harpy.get_PDF_tw3(x1, x2, Q, -2)-harpy.get_PDF_tw3(x1, x2, Q, -1)]]
    
#     with open(path_to_save+"deltaT_NS1."+str(n).zfill(4), 'w') as file:
#         for s in toSave:
#             file.write(', '.join("{:.8f}".format(item) for item in s)+'\n')
    