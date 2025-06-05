#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 20:15:46 2019
This program collect all the data on the SIDIS and save in  "SIDISdata_uncut.pkl"
@author: vla18041
"""
import sys
sys.path.append("/data/arTeMiDe_Repository/DataProcessor/")
import DataProcessor.Point
import DataProcessor.DataSet
import numpy

path_to_data="/data/arTeMiDe_Repository/data"
path_to_save="/data/arTeMiDe_Repository/DataProcessor/DataLib/wgt/"


path_to_HERMES="/HERMES-SSA/SFA/"

totalData=[]

M_proton=0.938
m_pion=0.139
m_kaon=0.494

#%%
#### determines the limits of the Q bin
#### Qmin^2 = MAX (Q^2min, xmin y min (s-M^2), xmin/(1-xmin)*(W2min-M^2))
#### Qmax^2 = MIN (Q^2max, xmax y max (s-M^2), xmax/(1-xmax)*(W2max-M^2))
def Qbounds(xMin,xMax):
    Q2min=1.
    Q2max=10000.
    WM2min=10.-(0.938)**2
    WM2max=10000.-(0.938)**2  ## no upper limit
    yMin=0.1
    yMax=0.95
    sM2=2*27.6*0.938
    
    if xMax<1:
        return [numpy.sqrt(numpy.max([Q2min,xMin*yMin*sM2, xMin/(1-xMin)*WM2min])),
                numpy.sqrt(numpy.min([Q2max,xMax*yMax*sM2, xMax/(1-xMax)*WM2max]))]
    else:
        return [numpy.sqrt(numpy.max([Q2min,xMin*yMin*sM2, xMin/(1-xMin)*WM2min])),
                numpy.sqrt(numpy.min([Q2max,xMax*yMax*sM2, 1000*WM2max]))]

#%%

### bins are presented by Gunar Schnell
### There should be taken the central VALUE!!! the bin integration is included in sys-error
xBin=[0.023,0.045,0.067,0.086,0.113,0.160,0.220,0.400]
zBin=[0.2,0.27,0.34,0.41,0.49,0.56,0.63,0.70]
ptBin=[0.00001,0.17,0.25,0.33,0.41,0.58,0.80,2.00]
### Scale uncertanty
scaleUncertanty=0.073

#%%
#### selecting only dat files
import os
listOfFiles=os.listdir(path_to_data+path_to_HERMES)
listOfFiles=[x for x in listOfFiles if '.txt' in x]


#%%
for ff in listOfFiles:
    f = open(path_to_data+path_to_HERMES+ff)
    
    data_from_f=[]
    
    for line in f:    
        data_from_f.append(line.rstrip('\n'))
    
    f.close()
    
    print("Done.  =>     Convert to numbers ...")
    
    #### fixing process
    if('pip' in ff):
        currentProcess=[1,1,1,13001]
        currentWeight=[1,1,1,2001]
        currentLabel2='pi+'
        currentD2=' pi+ from p'
    elif('pim' in ff):
        currentProcess=[1,1,-1,13001]
        currentWeight=[1,1,-1,2001]
        currentLabel2='pi-'
        currentD2=' pi- from p'
    elif('kp' in ff):
        currentProcess=[1,1,2,13001]
        currentWeight=[1,1,2,2001]
        currentLabel2='k+'
        currentD2=' k+ from p'
    elif('km' in ff):
        currentProcess=[1,1,-1,13001]
        currentWeight=[1,1,-1,2001]
        currentLabel2='k-'
        currentD2=' k- from p'
    else:
        print("CANNOT IDENTIFY HADRON! SKIP!")
        continue
    
    #### reading only lines related to the present interest
    toRead=False
    mainBulk=[]
    for ll in data_from_f:
        ### searching for interesting asymetry
        if("A_LT^cos(phi-phis)/sqrt(1-e^2) SFA DSA" in ll):
            toRead=True
            continue
        if(toRead):
            ### end the scan one another asymetry appears
            if('in 3D bins of x, z, Pt' in ll): break
            if('#' in ll): continue ## skip comented lines
            if(ll==''): continue  ## skip empty lines
            mainBulk.append(ll)
    
    for i in range(len(mainBulk)):
        mainBulk[i]=mainBulk[i].split()
        mainBulk[i][0]=[float(j) for j in mainBulk[i][0].split("<x<")]
        mainBulk[i][1]=[float(j) for j in mainBulk[i][1].split("<z<")]
        mainBulk[i][2]=[float(j) for j in mainBulk[i][2].split("<Pt<")]
        mainBulk[i][3:12]=[float(j) for j in mainBulk[i][3:12]]
     ### common variables
    s_current=2*27.6*0.938+(0.938)**2
    includeCuts=True
    cutParameters=[0.1,0.95,10.,10000.] #y, W^2 cuts
    
    
    DataCurrent=DataProcessor.DataSet.DataSet('hermes3D.ALT.'+currentLabel2,"SIDIS")
    DataCurrent.comment='HERMES 3D DSA-A.LT '+currentD2+". The data MUST be evaluated at a point"
    DataCurrent.reference="2007.07755"
    DataCurrent.normErr.append(scaleUncertanty)    


    for i in range(len(mainBulk)):
        # makeup a point
        p1=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(int(i)))
        #print DataCurrent.name+'.'+str(i)
        p1["process"]=currentProcess
        p1["s"]=s_current
        p1["<pT>"]=mainBulk[i][7]    
        p1["pT"]=mainBulk[i][2]
        p1["<x>"]=mainBulk[i][4]    
        p1["x"]=mainBulk[i][0]
        p1["<Q>"]=numpy.sqrt(mainBulk[i][3])
        p1["Q"]=Qbounds(p1["x"][0],p1["x"][1])    
        p1["<z>"]=mainBulk[i][6]
        p1["z"]=mainBulk[i][1]        
        p1["xSec"]=mainBulk[i][9]
        p1["M_target"]=M_proton
        p1["M_product"]=m_pion
        p1["includeCuts"]=includeCuts
        p1["cutParams"]=cutParameters    
        p1["thFactor"]=1.         ### tobe updated
        p1["uncorrErr"].append(mainBulk[i][10])
        p1["uncorrErr"].append(mainBulk[i][11])
        p1["weightProcess"]=currentWeight
        #
        DataCurrent.AddPoint(p1)    
    DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
    
