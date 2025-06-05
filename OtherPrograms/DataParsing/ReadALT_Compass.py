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

path_to_COMPASS="/COMPASS/1609.07374/ALT_cosHmS/"

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
    yMax=0.9
    sM2=2*160*0.938
    
    if xMax<1:
        return [numpy.sqrt(numpy.max([Q2min,xMin*yMin*sM2, xMin/(1-xMin)*WM2min])),
                numpy.sqrt(numpy.min([Q2max,xMax*yMax*sM2, xMax/(1-xMax)*WM2max]))]
    else:
        return [numpy.sqrt(numpy.max([Q2min,xMin*yMin*sM2, xMin/(1-xMin)*WM2min])),
                numpy.sqrt(numpy.min([Q2max,xMax*yMax*sM2, 1000*WM2max]))]
#%%
#### determines the limits of the x bin
#### xmin = MAX { xmin, Q^2/(ymax (s-M^2)), Q^2/(Q^2+W2max-M^2) }
#### xmax = MIN { xmax, Q^2/(ymin (s-M^2)), Q^2/(Q^2+W2min-M^2) }
def xbounds(Q2min,Q2max):
    xmin=0.003
    xmax=0.9  
    WM2min=10.-(0.938)**2
    WM2max=10000.-(0.938)**2  ## no upper limit
    yMin=0.1
    yMax=0.9
    sM2=2*160*0.938
    
    return [numpy.max( [xmin,Q2min/(yMax*sM2), Q2min/(Q2min+WM2max)] ),
            numpy.min( [xmax,Q2max/(yMin*sM2), Q2max/(Q2max+WM2min)] ) ]

#%%
#### selecting only dat files
import os
listOfFiles=os.listdir(path_to_data+path_to_COMPASS)
listOfFiles=[x for x in listOfFiles if '.dat' in x]
    
#%%
###############################################################################
for ff in listOfFiles:
    f = open(path_to_data+path_to_COMPASS+ff)
    
    data_from_f=[]
    
    for line in f:    
        data_from_f.append(line.rstrip('\n'))
    
    f.close()
    
    print("Done.  =>     Convert to numbers ...")
    
    ### parsing parameters
    lineM=data_from_f[3]
    
    #### fixing Z
    if('z > 0.2' in lineM):
        zCurrent=[0.2,1.]
        currentLabel1='2<z'
        currentD1='(0.2<z<1)'
    elif('z > 0.1' in lineM):
        zCurrent=[0.1,1.]
        currentLabel1='1<z'
        currentD1='(0.1<z<1)'
    elif('0.1 < z < 0.2' in lineM):
        zCurrent=[0.1,0.2]
        currentLabel1='1<z<2'
        currentD1='(0.1<z<0.2)'
    else:
        print("CANNOT IDENTIFY Z")
        
    #### fixing process
    if('h^+' in lineM):
        currentProcess=[1,1,1,13101]
        currentWeight=[1,1,1,2101]
        currentLabel2='h+'
        currentD2=' h+ from p'
    elif('h^-' in lineM):
        currentProcess=[1,1,-1,13101]
        currentWeight=[1,1,-1,2101]
        currentLabel2='h-'
        currentD2=' h- from p'
    else:
        print("CANNOT IDENTIFY HADRON!")
        
    #### label
    if('p_T-dependence' in lineM):
        currentLabel3='dpt'
        currentD3='(differential in pt)'
    elif('x -dependence' in lineM):
        currentLabel3='dx'
        currentD3='(differential in x)'
    elif('z -dependence' in lineM):
        currentLabel3='dz'
        currentD3='(differential in z)'
    else:
        print("CANNOT IDENTIFY CASE!")
    
    ### common variables
    s_current=2*160*0.938+(0.938)**2
    includeCuts=True
    cutParameters=[0.1,0.9,10.,10000.] #y, W^2 cuts
    
    DataCurrent=DataProcessor.DataSet.DataSet('compass16.ALT.'+currentLabel2+'.'+currentLabel1+'.'+currentLabel3,"SIDIS")
    DataCurrent.comment='COMPASS16 DSSA-A.LT '+currentD2+', '+currentD1+' '+currentD3
    DataCurrent.reference="1609.07374"
    DataCurrent.normErr.append(0.03)#### delution factor according to Bakur
    
    
    ### Now we go line by line
    mainBulk=data_from_f[5:]
    num=0
    for line in mainBulk:
        ## empty line
        if(line==''): continue 
        ## header line
        if('<x>\t<y>\t<z>\t<p_T>' in line): continue
        ## energy dependance
        if('Q^{2}/(GeV/c)^{2}' in line):
            if('1<Q^{2}/(GeV/c)^{2}<4' in line):
                Q2_current=[1.,4.]
            elif('4<Q^{2}/(GeV/c)^{2}<6.25' in line):
                Q2_current=[4.,6.25]
            elif('6.25<Q^{2}/(GeV/c)^{2}<16' in line):
                Q2_current=[6.25,16.]
            elif('16<Q^{2}/(GeV/c)^{2}<81' in line):
                Q2_current=[16.,81.]
            else:
                raise ValueError('CANNOT DETERMINE Q')
            Q_current=[numpy.sqrt(Q2_current[0]),numpy.sqrt(Q2_current[1])]
        ### normal line
        elif(line[0]=='['):
            lineParsed=line.replace("[","").replace("]","").replace(";"," ").split()  
            lineParsed=[float(j) for j in lineParsed]
            
            ###### Adding point
            p=DataProcessor.Point.CreateSIDISPoint(DataCurrent.name+'.'+str(num))
            num+=1
            #print DataCurrent.name+'.'+str(i)
            p["process"]=currentProcess
            p["weightProcess"]=currentWeight
            p["s"]=s_current
            p["Q"]=Q_current        
            
            p["x"]=xbounds(Q2_current[0],Q2_current[1])
            p["z"]=zCurrent        
            p["pT"]=[0.1,10.] 
            if(currentLabel3=='dx'):
                p["x"]=[lineParsed[0],lineParsed[1]]
            elif(currentLabel3=='dz'):
                p["z"]=[lineParsed[0],lineParsed[1]]
            elif(currentLabel3=='dpt'):
                p["pT"]=[lineParsed[0],lineParsed[1]]
            else:
                raise ValueError('ERROR1')
            p["<x>"]=lineParsed[2]
            ##[3]=<y>
            p["<z>"]=lineParsed[4]
            p["<pT>"]=lineParsed[5]
            ##[6]=<W>
            p["<Q>"]=numpy.sqrt(lineParsed[7])
            p["xSec"]=lineParsed[8]
            p["M_target"]=M_proton
            p["M_product"]=m_pion
            p["includeCuts"]=includeCuts
            p["cutParams"]=cutParameters    
            p["thFactor"]=1.         ### tobe updated
            p["uncorrErr"].append(lineParsed[10])
            
            
            DataCurrent.AddPoint(p)
            
        else:
            raise ValueError('DO NOT UNDERSTAD hte Line')
                
    DataCurrent.SaveToCSV(path_to_save+DataCurrent.name+".csv")
