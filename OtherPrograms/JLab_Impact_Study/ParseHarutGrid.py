#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 08:24:58 2022

@author: vla18041
"""


import sys
sys.path.append("/home/vla18041/LinkData2/arTeMiDe_Repository/DataProcessor/")

import DataProcessor.Point
import DataProcessor.DataSet
import numpy

path_to_data="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/JLab_Impact/PseudoData/jlab22/"

path_to_save="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/JLab_Impact/PseudoData/jlab22/"

## this trigger add normalization uncertanty due to DIS normalization for multiplicities
addDISnormalizationUncertainty=False

M_proton=0.938
m_pion=0.139
m_kaon=0.494

#%%
######
### definition of bins
##c    nq2=19  q2min=1, Q2max=120, dq2=1  q2min+(i-1)*dq2<Q2_i< q2min+i*dq2    (ex bin i=6 means 6<Q^2<7)
Q2min=1.
dQ2=1
def Q2bin(i):
    return [Q2min+(i-1)*dQ2, Q2min+i*dQ2]
def Qbin(i):
    return [numpy.sqrt(Q2min+(i-1)*dQ2), numpy.sqrt(Q2min+i*dQ2)]
##c    nx=20, xmin=0.04, dx=0.04                   xmin+(j-1)*dx<x_j< xmin+j*dx (ex. bin j=10 means 0.4<x<0.44 )
xmin=0.04
dx=0.04
def xbin(i):
    return [xmin+(i-1)*dx, xmin+i*dx]
##c    nz=18, zmin=0.05, dz=0.05                   zmin+(k-1)*dz<z_k< zmin+j*dz (ex. bin k=10 means 0.5<z<0.55) 
zmin=0.05
dz=0.05
def zbin(i):
    return [zmin+(i-1)*dz, zmin+i*dz]
##c    npt=40, ptmin=0, step =0.045                ptmin+(l-1)*dpt<Pt_l< Ptmin+l*dpt (ex. bin l=10 means 0.405<z<0.45)
ptmin=0.
dpt=0.045
def ptbin(i):
    return [ptmin+(i-1)*dpt, ptmin+i*dpt]
  
#%%
##################Cut function
#### drop all points with delta>0.3, and with zero number of events
def cutFunc(p):    
    if p["uncorrErr"][0]>1.:
        return False, p
    else:
        gamma2=(2.0*p["M_target"]*p["<x>"]/p["<Q>"])**2
        rho2=(p["M_product"]/p["<z>"]/(p["<Q>"]))**2
        if((1-gamma2*rho2)<0.0001):
            return False,p
        qT=p["<pT>"]/p["<z>"]*numpy.sqrt((1+gamma2)/(1-gamma2*rho2))
        delta=qT/(p["<Q>"])
        
        ### compute the largest possible qT (approximate)
        gamma2WORST=(2.0*p["M_target"]*p["x"][1]/p["<Q>"])**2
        # it is definitely not a TMD point
        if(1-gamma2WORST*rho2<0.0001):
            return False , p
        qTWORST=p["pT"][1]/p["z"][0]*numpy.sqrt((1+gamma2WORST)/(1-gamma2WORST*rho2))

        ## drop if qT>Q/2
        if qTWORST>p["<Q>"]/2:
            return False , p

    ### drop Q<2
    if p["<Q>"]<2 :
        return False , p
    
#    return delta<0.5 and p.qT_avarage<80
    return delta<0.3 , p

#%%
#### The routine which parse the Harut's grid
#### It also save the version with cut-application (delta<0.3)
def ParseHarut(file,name,process):
    f = open(path_to_data+file)

    data_from_f=[]
    
    for line in f:    
        data_from_f.append(line.rstrip('\n'))
    
    f.close()

    print("Done.  =>     Convert to numbers ...")
    
    for i in range(len(data_from_f)):
        data_from_f[i]=data_from_f[i].split()
        #data_from_f[i]=[float(j) for j in data_from_f[i]]
        #k.spit("\t")
    
    print("Done.  =>     Create points & append to data set ...")
    DataCurrent=DataProcessor.DataSet.DataSet(name,"SIDIS")
    DataCurrent.comment="JLab22 pseudo data (only stat uncertanty) by Harut"
    DataCurrent.reference="Harut"
    
    DataCurrentCUT=DataProcessor.DataSet.DataSet(name+'.cut',"SIDIS")
    DataCurrentCUT.comment="JLab22 pseudo data (only stat uncertanty) by Harut. CUT TO TMD"
    DataCurrentCUT.reference="Harut"
    
    proc_current=process
    s_current=2*22.*0.938+(0.938)**2
    includeCuts=False
    cutParameters=[0.1,0.85,10.,10000.]
    
    for i in range(len(data_from_f)):
        # makeup a point
        p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+data_from_f[0][0])
        #print DataCurrent.name+'.'+str(i)
        p["process"]=proc_current
        p["s"]=s_current
        p["<pT>"]=float(data_from_f[i][9])
        p["pT"]=ptbin(int(data_from_f[i][5]))
        p["<Q>"]=numpy.sqrt(float(data_from_f[i][6]))
        p["Q"]=Qbin(int(data_from_f[i][2]))
        p["<x>"]=float(data_from_f[i][7])
        p["x"]=xbin(int(data_from_f[i][3]))
        p["<z>"]=float(data_from_f[i][8])
        p["z"]=zbin(int(data_from_f[i][4]))
        ## cross-seciton is unknown
        p["xSec"]=0.1
        p["M_target"]=M_proton
        p["M_product"]=m_pion
        p["includeCuts"]=includeCuts
        p["cutParams"]=cutParameters
        #devide by (x,z,pt^2,Q^2) bin size
        p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["Q"][1]**2-p["Q"][0]**2)/(p["z"][1]-p["z"][0])/(p["x"][1]-p["x"][0])
        
        #### HERE UNCERTANTY FOR 1 DAY TAKING
        ###12) Total # of events in that bin 748 below
        if(float(data_from_f[i][11])>0):
            p["uncorrErr"].append(1./numpy.sqrt(float(data_from_f[i][11])))
            ###use systematics 0.5*stat-error
            p["uncorrErr"].append(0.5/numpy.sqrt(float(data_from_f[i][11])))
        else:
            ### if bin is empty erro=1000%
            p["uncorrErr"].append(10.)
            ###use systematics 0.5*stat-error
            p["uncorrErr"].append(0.5*10)
        if(addDISnormalizationUncertainty):
            if(data_from_f[i][0]==0):
                pass#p.corrErrors.append(0)
            else:
                p["corrErr"].append(data_from_f[i][16]/data_from_f[i][15]*data_from_f[i][0])
        #
        DataCurrent.AddPoint(p)
        r1,r2=cutFunc(p)
        if(r1):
            DataCurrentCUT.AddPoint(p)
    
    print("Done.  ")
    
    DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
    DataCurrentCUT.SaveToCSV(path_to_save+DataCurrentCUT.name+".csv")

#%%
ParseHarut("pro.gen.sidis.pip.dat", "jlab22.gen.pi+", [1,1,2001])
ParseHarut("pro.gen.sidis.pim.dat", "jlab22.gen.pi-", [1,1,2021])
ParseHarut("pro.gen.sidis.kap.dat", "jlab22.gen.k+", [1,1,2002])
ParseHarut("pro.gen.sidis.kam.dat", "jlab22.gen.k-", [1,1,2022])

#%%
#### sets with joined bins
### definition of bins
##c    nx=20, xmin=0.04, dx=0.04                   xmin+(j-1)*dx<x_j< xmin+j*dx (ex. bin j=10 means 0.4<x<0.44 )
def xbinJOIN(i,j):
    return [xmin+j*(i-1)*dx, xmin+j*i*dx]
##c    nz=18, zmin=0.05, dz=0.05                   zmin+(k-1)*dz<z_k< zmin+j*dz (ex. bin k=10 means 0.5<z<0.55) 
def zbinJOIN(i,j):
    return [zmin+j*(i-1)*dz, zmin+j*i*dz]
##c    npt=40, ptmin=0, step =0.045                ptmin+(l-1)*dpt<Pt_l< Ptmin+l*dpt (ex. bin l=10 means 0.405<z<0.45)
def ptbinJOIN(i,j):
    return [ptmin+j*(i-1)*dpt, ptmin+j*i*dpt]

#%%
#### The routine which parse the Harut's grid
#### It also save the version with cut-application (delta<0.3)
def ParseHarutJOIN(file,name,process,joinX,joinZ,joinT):
    f = open(path_to_data+file)

    data_from_f=[]
    
    for line in f:    
        data_from_f.append(line.rstrip('\n'))
    
    f.close()

    print("Done.  =>     Convert to numbers ...")
    
    for i in range(len(data_from_f)):
        data_from_f[i]=data_from_f[i].split(',')
        #data_from_f[i]=[float(j) for j in data_from_f[i]]
        #k.spit("\t")
    
    print("Done.  =>     Create points & append to data set ...")
    DataCurrent=DataProcessor.DataSet.DataSet(name,"SIDIS")
    DataCurrent.comment="JLab22 pseudo data (only stat uncertanty) by Harut"
    DataCurrent.reference="Harut"
    
    DataCurrentCUT=DataProcessor.DataSet.DataSet(name+'.cut',"SIDIS")
    DataCurrentCUT.comment="JLab22 pseudo data (only stat uncertanty) by Harut. CUT TO TMD"
    DataCurrentCUT.reference="Harut"
    
    proc_current=process
    s_current=2*22.*0.938+(0.938)**2
    includeCuts=False
    cutParameters=[0.1,0.85,10.,10000.]
    
    for i in range(len(data_from_f)):
        # makeup a point
        p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+data_from_f[0][0])
        #print DataCurrent.name+'.'+str(i)
        p["process"]=proc_current
        p["s"]=s_current
        p["<pT>"]=float(data_from_f[i][9])
        p["pT"]=ptbinJOIN(int(data_from_f[i][5]),joinT)
        p["<Q>"]=numpy.sqrt(float(data_from_f[i][6]))
        p["Q"]=Qbin(int(data_from_f[i][2]))
        p["<x>"]=float(data_from_f[i][7])
        p["x"]=xbinJOIN(int(data_from_f[i][3]),joinX)
        p["<z>"]=float(data_from_f[i][8])
        p["z"]=zbinJOIN(int(data_from_f[i][4]),joinZ)
        ## cross-seciton is unknown
        p["xSec"]=0.1
        p["M_target"]=M_proton
        p["M_product"]=m_pion
        p["includeCuts"]=includeCuts
        p["cutParams"]=cutParameters
        #devide by (x,z,pt^2,Q^2) bin size
        p["thFactor"]=1/(p["pT"][1]**2-p["pT"][0]**2)/(p["Q"][1]**2-p["Q"][0]**2)/(p["z"][1]-p["z"][0])/(p["x"][1]-p["x"][0])
        
        #### HERE UNCERTANTY FOR 1 DAY TAKING
        ###12) Total # of events in that bin 748 below
        if(float(data_from_f[i][11])>0):
            p["uncorrErr"].append(1./numpy.sqrt(float(data_from_f[i][11])))
            ###use systematics 0.5*stat-error
            p["uncorrErr"].append(0.5/numpy.sqrt(float(data_from_f[i][11])))
        else:
            ### if bin is empty erro=1000%
            p["uncorrErr"].append(10.)
            ###use systematics 0.5*stat-error
            p["uncorrErr"].append(0.5*10)
        if(addDISnormalizationUncertainty):
            if(data_from_f[i][0]==0):
                pass#p.corrErrors.append(0)
            else:
                p["corrErr"].append(data_from_f[i][16]/data_from_f[i][15]*data_from_f[i][0])
        #
        DataCurrent.AddPoint(p)
        r1,r2=cutFunc(p)
        if(r1):
            DataCurrentCUT.AddPoint(p)
    
    print("Done.  ")
    
    DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
    DataCurrentCUT.SaveToCSV(path_to_save+DataCurrentCUT.name+".csv")
    
#%%
ParseHarutJOIN("pro.gen.sidis.pip.join_5_2_1.dat", "jlab22.gen.pi+.join521", [1,1,2001],5,2,1)
ParseHarutJOIN("pro.gen.sidis.pim.join_5_2_1.dat", "jlab22.gen.pi-.join521", [1,1,2021],5,2,1)
ParseHarutJOIN("pro.gen.sidis.kap.join_5_2_1.dat", "jlab22.gen.k+.join521", [1,1,2002],5,2,1)
ParseHarutJOIN("pro.gen.sidis.kam.join_5_2_1.dat", "jlab22.gen.k-.join521", [1,1,2022],5,2,1)