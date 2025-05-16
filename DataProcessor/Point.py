#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 16:01:12 2020

Collection of methods for the operation with data points.

Each point has key "type" which define the class of a point.
"type" MUST BE preserved.

List of common keys (!=not self-filling):
    "type":     type of process (!!!)
    "id":       name of poins (!)
    "process":  process code (!)
    "xSec":     experimental cross-section (!)
    "uncorrErr":   list of uncorrelated errors [err1,err2,...] (!)
    "corrErr":     list of point-to-point correlated errors [err1,err2,...] (!)
                    WITHOUT NORMALIZATION ERRORS(they are defined in DataSet)
    "thFactor": the factor by which the theory should be multiplied 
                    to compare with the data. The theory is expected to be
                    \int d\sigma dBIN (!)
    "s":        Mandeshtam variable s (!)
    "Q":        Bin limits for the photon virtuality [Qmin,Qmax] (!)
    "<Q>":      Avarage Q (definition is not fixed)
    "includeCuts":   Inclusion of the fiducial cuts (logical) (!)
    "cutParams":    Parameters of the cut [a,b,c,d,...]  (! if include cuts=true)
    "qT":       range of qT bin [qTmin,qTmax] (! for DY)
    "<qT>":     Avarage qT (definition is not fixed)

List of DY specific keys:
    "y":        range of rapidity bin [ymin,ymax] (!)
    "<y>":      Avarage y (definition is not fixed)

List of SIDIS specific keys:
    "x":        range of x-bin [xmin,xmax] (!)
    "<x>":      Avarage x (definition is not fixed)
    "z":        range of z-bin [zmin,zmax] (!)
    "<z>":      Avarage z (definition is not fixed)
    "pT":       range of pT bin [pTmin,pTmax] (!)
    "<pT>":     Avarage pT (definition is not fixed)
    "M_target": Mass of target (if not set =proton)
    "M_product":Mass of produced hadron (if not set =pion)


@author: vla18041
"""

def CreateDYPoint(name):
    """Creates dictionary with type="DY" (preserve), and empty lists.
    
    name: identifier for the point (string)
    """
    return {
            "type":"DY",
            "id":name,
            "uncorrErr":[],
            "corrErr":[]
            }
    
def CreateSIDISPoint(name):
    """Creates dictionary with type="SIDIS" (preserve), and empty lists
    
    name: identifier for the point (string)
    """
    return {
            "type":"SIDIS",
            "id":name,
            "uncorrErr":[],
            "corrErr":[]
            }

def CreateD2Point(name):
    """Creates dictionary with type="D2" (preserve), and empty lists.
    
    name: identifier for the point (string)
    """
    return {
            "type":"D2",
            "id":name,
            "uncorrErr":[],
            "corrErr":[]
            }

def CreateG2Point(name):
    """Creates dictionary with type="G2" (preserve), and empty lists.
    
    name: identifier for the point (string)
    """
    return {
            "type":"D2",
            "id":name,
            "uncorrErr":[],
            "corrErr":[]
            }
  
def FinalizePoint(point, framework="artemide"):
    """Checks the presence and type of all must-be fields in a point-dictionary.
    
    Automatically fills not-filled fields (which are not marked by !)
    Return FALSE is there are unset fields that must be set (marked with (!))
    Return TRUE if all OK
    
    output = logical
    input  = point dictionary
    
    framework = specify the possible framework to operate with the point and
                checks the types of some variables
                default = "none"
                added = "artemide"
    """
    
    #Elementary checks
    if not isinstance(point,dict):
        print("FinalizePoint: argument must be dictionary")
        return False
    
    if not "type" in point:
        print("FinalizePoint: point-dictionary must have key 'type'")
        return False
    
    if not(point["type"] in ["DY","SIDIS","G2"]):        
        print("FinalizePoint: point-dictionary must have 'type'= ['DY','SIDIS','G2']")
        return False
    
    ######################################################################
    #Checking common fields
    #---- id ----
    if not "id" in point:
        print("FinalizePoint: 'id' is not present")
        return False
    
    if not isinstance(point.get("id"),str):
        print("FinalizePoint: 'id' is not a string")
        return False
    
    #---- "process" -----    
    if not "process" in point:
        print("FinalizePoint: 'process' is not present")
        return False
    
    ## artemide set the process by 3 integers
    if framework=="artemide":
        a=point.get("process")
        if(point["type"] in ["DY","SIDIS"]):
            if not isinstance(a,list):
                print("FinalizePoint(framework=artemide): 'process' is not a list")
                return False
            if (len(a)!=3 and len(a)!=4):
                print("FinalizePoint(framework=artemide): 'process' is not a list of length 4")
                return False
            if not all(isinstance(elem,int) for elem in a):
                print("FinalizePoint(framework=artemide): 'process' is not a list of integers")
                return False
        if(point["type"] in ["D2","G2"]):
            if not isinstance(a,int):
                print("FinalizePoint(framework=artemide): 'process' is not a list of integers")
                return False
    
    #---- xSec -----
    if not "xSec" in point:
        print("FinalizePoint: 'xSec' is not present")
        return False
    
    if not isinstance(point.get("xSec"),float):
        print("FinalizePoint: 'xSec' is not a float")
        return False
    
    #---- uncorrErr -----
    if not "uncorrErr" in point:
        print("FinalizePoint: 'uncorrErr' is not present")
        return False
    
    a=point.get("uncorrErr")
    if not isinstance(a,list):
        print("FinalizePoint: 'uncorrErr' is not a list")
        return False
    
    if len(a)>0 and not(all(isinstance(elem,float) for elem in a)):
        print("FinalizePoint: 'uncorrErr' has non-float element in a list")
        return False
    
    #---- corrErr -----
    if not "corrErr" in point:
        print("FinalizePoint: 'corrErr' is not present")
        return False
    
    a=point.get("corrErr")
    if not isinstance(a,list):
        print("FinalizePoint: 'corrErr' is not a list")
        return False
    
    if len(a)>0 and not(all(isinstance(elem,float) for elem in a)):
        print("FinalizePoint: 'corrErr' has non-float element in a list")
        return False
    
    #---- thFactor -----
    if not "thFactor" in point:
        print("FinalizePoint: 'thFactor' is not present")
        return False
    
    if not isinstance(point.get("thFactor"),float):
        print("FinalizePoint: 'thFactor' is not a float")
        return False
    
    #---- s -----
    if(point["type"] in ["DY","SIDIS"]):
        if not "s" in point:
            print("FinalizePoint: 's' is not present")
            return False
        
        if not isinstance(point.get("s"),float):
            print("FinalizePoint: 's' is not a float")
            return False
    
    #---- Q -----
    if not(_TestBin(point,"Q","<Q>")):
        return False
        
    #---- includeCuts -----
    if(point["type"] in ["DY","SIDIS"]):
        if not "includeCuts" in point:
            print("FinalizePoint: 'includeCuts' is not present")
            return False
        
        if not isinstance(point.get("includeCuts"),bool):
            print("FinalizePoint: 'includeCuts' is not a bool")
            return False
        
        #---- includeCuts -----
        if not "cutParams" in point:
            if point["includeCuts"]:
                print("FinalizePoint(includeCuts=true): 'cutParams' are not present")
                return False
            else:
                # if cuts are not needed fill with some generic stuff
                if framework=="artemide":
                    if point["type"]=="DY":
                        point["cutParams"]=[0.,0.,-1000.,+1000.]
                    elif point["type"]=="SIDIS":
                        point["cutParams"]=[0.,1.,0.,+100000.]
                else:
                    point["cutParams"]=[]
        
        if framework=="artemide":                
            a=point.get("cutParams")
            if not(isinstance(a,list)) or len(a)!=4 or not(all(isinstance(elem,float) for elem in a)):
                print("FinalizePoint(framework=artemide): 'cutParams' is not [float,float,float,float]")
                return False
            
            if point["type"]=="DY" and a[2]>a[3]:
                print("FinalizePoint(framework=artemide): 'cutParams' in DY etaMin>etaMax")
                return False
            if point["type"]=="SIDIS" and a[2]>a[3]:
                print("FinalizePoint(framework=artemide): 'cutParams' in SIDIS WMin>WMax")
                return False
            if point["type"]=="SIDIS" and a[0]>a[1]:
                print("FinalizePoint(framework=artemide): 'cutParams' in SIDIS yMin>yMax")
                return False
            
    ######################################################################
    #Checking DY fields
    if point["type"]=="DY":
    #---- qT -----
        if not(_TestBin(point,"qT","<qT>")):
            return False
    #---- y -----
        if not(_TestBin(point,"y","<y>")):
            return False
        
    ######################################################################
    #Checking SIDIS fields
    if point["type"]=="SIDIS":
    #---- x -----
        if not(_TestBin(point,"x","<x>")):
            return False
    #---- z -----
        if not(_TestBin(point,"z","<z>")):
            return False
    #---- pT -----
        if not(_TestBin(point,"pT","<pT>")):
            return False
    #---- qT -----
        if not "qT" in point:
            point["qT"]=[ point["pT"][0]/point["<z>"] , point["pT"][1]/point["<z>"] ]
        if not(_TestBin(point,"qT","<qT>")):
            return False
    #---- M_target -----
        if not "M_target" in point:
            point["M_target"]=0.938
    #---- M_product -----
        if not "M_product" in point:
            point["M_product"]=0.139
            
    ######################################################################
    #Checking D2 fields
    if point["type"]=="D2":
    #----no specific variables
        pass
    
    #Checking G2 fields
    if point["type"]=="G2":
    #----x ------------
        if not(_TestBin(point,"x","<x>")):
            return False
    
    return True

###############################################################################
#################  PRIVATE FUNCTIONS  #########################################
###############################################################################

def _TestBin(p,n1,n2):
    """Help-function that checks that in the dictionary p, the field n1 is
    present, [float,float], [a<b]. If field n2 absent fill with avarage
    """
    if not n1 in p:
        print("FinalizePoint: '",n1,"' is not present")
        return False
    
    a=p.get(n1)
    if not(isinstance(a,list)) or len(a)!=2 or not(all(isinstance(elem,float) for elem in a)):
        print("FinalizePoint: '",n1,"' is not [float,float]")
        return False
    
    if a[0]>a[1]:
        print("FinalizePoint: '",n1,"' is not ordered")
        return False
    
    #---- avarage -----
    if n2 in p:
        if not(isinstance(p[n2],float)):
            print("FinalizePoint: '",n2,"' is not a float")
            return False
    else:
        # set n2 as middle of bin
        p[n2]=(p[n1][0]+p[n1][1])/2.
    
    return True