#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 15:11:15 2020

@author: vla18041
"""

#%%
#######################################
# importing libraries
#######################################

import sys
import time
sys.path.append("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor")
import DataProcessor.harpyInterface
import DataProcessor.ArtemideReplicaSet


#prefix="ff_K_n3lo_all=0"
prefix="n3lo_all=0"

MAINPATH="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/OtherPrograms/TMDGrids/"
SAVEPATH="/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/OtherPrograms/TMDGrids/Grids/SV19_"+prefix+"/SV19_"+prefix

#%%
#######################################
# Initializing artemide
#######################################
import harpy
path_to_constants=MAINPATH+"Constants/"
#harpy.initialize(path_to_constants+"const-TMDPDF-NNPDF")
harpy.initialize(path_to_constants+"const-TMDPDF-NNPDF-n3lo")
#harpy.initialize(path_to_constants+"const-TMDFF-DSS-n3lo")
#rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/Replicas/DY+SIDIS/SV19_nnlo_all=0.rep")
#rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/Replicas/DY+SIDIS/SV19_nnlo.rep")
#rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/Replicas/DY+SIDIS/SV19_n3lo.rep")
rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile("/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/Replicas/DY+SIDIS/SV19_n3lo_all=0.rep")
rSet.SetReplica(num=0,part="uTMDPDF")

#%%
#######################################
# Save Grid specification
#######################################
with open(SAVEPATH+'.info', 'w') as outfile:
    outfile.write("SetDesc: Unpolarized TMDFF (kaon+) SV19_"+prefix+"\n")
    outfile.write("Authors: I.Scimemi, A.Vladimirov"+"\n"),
    outfile.write("Reference: arXiv:1912.06532"+"\n"),
    outfile.write("SetIndex: 706000"+"\n"),
    outfile.write("TMDScheme: Pavia TMDs"+"\n"),
    outfile.write("TMDType: ff"+"\n"),
    outfile.write("CollDist: DSS_2017_nnlo"+"\n"),
    outfile.write("CollDistMember:  0"+"\n"),
    outfile.write("Format: TMDlib2"+"\n"),
    outfile.write("DataVersion: 1"+"\n"),
    outfile.write("OrderQCD: N3LO & zeta-prec."+"\n"),
    outfile.write("Regularisation:  xxx "+"\n"),
    outfile.write("NumMembers: "+str(rSet.numberOfReplicas+1)+"\n"),
    outfile.write("ErrorType: Monte Carlo"+"\n"),
    outfile.write("FlavorScheme: LHAPDF style"+"\n"),
    outfile.write("Flavors: [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]"+"\n"),
    outfile.write("NumFlavors: 5"+"\n"),
    outfile.write("XMin:  0.05"+"\n"),
    outfile.write("XMax: 1."+"\n"),
    outfile.write("QMin: 1."+"\n"),
    outfile.write("QMax: 200."+"\n"),
    outfile.write("KtoQMin: 0.0001"+"\n"),
    outfile.write("KtoQMax: 2.001"+"\n")
  
#%%
#######################################
# Specifying the grid parameters
#######################################
Qrange= [1., 1.11803, 1.22474, 1.4, 1.58114, 1.78885, 2., 2.23607, 2.52982, 2.82843,
         3.16228, 3.4641, 4.75, 5.09902, 6.32456, 7.1, 8., 10., 11.1803, 12.2475,
         14., 15.8114, 17.8885, 20., 22.3607, 25.2982, 28.2843, 31.6228, 34.641, 47.5,
         50.9902, 63.2456, 71, 80, 100, 111.803, 122.475, 140, 158.114, 178.885, 
         200.]
Xrange= [0.00001, 0.00002, 0.00004, 0.00006, 0.00008, 0.0001, 0.0002, 0.0004, 0.0006, 0.0008,
          0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.0045, 0.005, 0.0055,
          0.006, 0.0065, 0.007, 0.0075, 0.008, 0.0085, 0.009, 0.00925, 0.0095, 0.00975,
          0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055,
          0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.0925, 0.095, 0.0975,
          0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
          0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.925, 0.95, 0.975,
          1]
# Xrange = [0.05,0.055,0.06,0.065,0.07,0.08,0.09,0.1,
#           0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
#           0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.]
Rrange= [0.0001, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.02, 0.03, 0.04, 0.05,
         0.06, 0.07, 0.08, 0.09, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225,
         0.25, 0.275, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
         0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
         1.7, 1.8, 1.9, 2.001]

#%%
#######################################
# Make a huge dict for TMD-values (for future convinience)
#######################################
from numpy import zeros
valuesList = {
    -5: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist(),
    -4: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist(),
    -3: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist(),
    -2: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist(),
    -1: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist(),
    1: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist(),
    2: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist(),
    3: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist(),
    4: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist(),
    5: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist()
    }
#%%
#######################################
# Fill the list with values (could take long time!)
#######################################
startTime=time.time()
for i in range(len(Qrange)):
    for j in range(len(Xrange)):
        for k in range(len(Rrange)):
            Qval=float(Qrange[i])
            xval=Xrange[j]
            rval=Rrange[k]
            
            if(xval==1):
                TMDval=[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
            else:
                TMDval=harpy.get_uTMDPDF_kT(xval,rval*Qval,1,Qval,Qval**2,includeGluon=False)
                #TMDval=harpy.get_uTMDFF_kT(xval,rval*Qval,2,Qval,Qval**2,includeGluon=False)
            
            valuesList[-5][i][j][k]='{:g}'.format(xval*TMDval[0])
            valuesList[-4][i][j][k]='{:g}'.format(xval*TMDval[1])
            valuesList[-3][i][j][k]='{:g}'.format(xval*TMDval[2])
            valuesList[-2][i][j][k]='{:g}'.format(xval*TMDval[3])
            valuesList[-1][i][j][k]='{:g}'.format(xval*TMDval[4])
            #valuesList[0][i][j][k]='{:g}'.format(xval*TMDval[5])
            valuesList[1][i][j][k]='{:g}'.format(xval*TMDval[6])
            valuesList[2][i][j][k]='{:g}'.format(xval*TMDval[7])
            valuesList[3][i][j][k]='{:g}'.format(xval*TMDval[8])
            valuesList[4][i][j][k]='{:g}'.format(xval*TMDval[9])
            valuesList[5][i][j][k]='{:g}'.format(xval*TMDval[10])
            
endTime=time.time()
print('Computation time : ',endTime-startTime,' sec.')

#%%
#######################################
# Save replica 0000
#######################################
startTime=time.time()
with open(SAVEPATH+'_0000.dat', 'w') as outfile:
    outfile.write("Qg: "+str(Qrange)+"\n")
    outfile.write("xg: "+str(Xrange)+"\n")
    outfile.write("qToQg: "+str(Rrange)+"\n")
    outfile.write("TMDs: "+str(valuesList).replace("'","")+"\n")
endTime=time.time()
print('Dump time : ',endTime-startTime,' sec.')    

#%%
#######################################
# Now do it for each replica!
#######################################
for r in range(rSet.numberOfReplicas):
    #rSet.SetReplica(num=r+1,part="uTMDPDF")
    rSet.SetReplica(num=r+1,part="uTMDFF")
    print("Doing replica : ",str(r+1),"/",str(rSet.numberOfReplicas))
    startTime=time.time()
    for i in range(len(Qrange)):
        for j in range(len(Xrange)):
            for k in range(len(Rrange)):
                Qval=float(Qrange[i])
                xval=Xrange[j]
                rval=Rrange[k]
                
                if(xval==1):
                    TMDval=[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
                else:
                    #TMDval=harpy.get_uTMDPDF_kT(xval,rval*Qval,1,Qval,Qval**2,includeGluon=False)
                    TMDval=harpy.get_uTMDFF_kT(xval,rval*Qval,2,Qval,Qval**2,includeGluon=False)
                
                valuesList[-5][i][j][k]='{:g}'.format(xval*TMDval[0])
                valuesList[-4][i][j][k]='{:g}'.format(xval*TMDval[1])
                valuesList[-3][i][j][k]='{:g}'.format(xval*TMDval[2])
                valuesList[-2][i][j][k]='{:g}'.format(xval*TMDval[3])
                valuesList[-1][i][j][k]='{:g}'.format(xval*TMDval[4])
                #valuesList[0][i][j][k]='{:g}'.format(xval*TMDval[5])
                valuesList[1][i][j][k]='{:g}'.format(xval*TMDval[6])
                valuesList[2][i][j][k]='{:g}'.format(xval*TMDval[7])
                valuesList[3][i][j][k]='{:g}'.format(xval*TMDval[8])
                valuesList[4][i][j][k]='{:g}'.format(xval*TMDval[9])
                valuesList[5][i][j][k]='{:g}'.format(xval*TMDval[10])
                
    with open(SAVEPATH+'_'+'{:04d}'.format(r+1)+'.dat', 'w') as outfile:
        outfile.write("Qg: "+str(Qrange)+"\n")
        outfile.write("xg: "+str(Xrange)+"\n")
        outfile.write("qToQg: "+str(Rrange)+"\n")
        outfile.write("TMDs: "+str(valuesList).replace("'","")+"\n")
    endTime=time.time()
    print('Computation time : ',endTime-startTime,' sec.')