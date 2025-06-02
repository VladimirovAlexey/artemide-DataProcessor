#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 18 09:50:56 2025

THIS IS THE SAME AS snowInterface, but for artemide+snowflake compilation

@author: alexey
"""

import harpy
import numpy
from . import DataSet
from . import DataMultiSet

def ComputeXSec(data,method="default"):
    """Computes the cross-section values for the given data and match it to data

    Parameters
    ----------
    data : DataSet or DataMultiSet or a point-like dictionary    
    method : String, optional
        Parameters for evaluation of xSec. The default is "default".
        default = usual one
        binless = evaluated with avarage values of bins only

    Returns
    -------
    YY : TYPE
        DESCRIPTION.

    """
    if isinstance(data,DataSet.DataSet) or isinstance(data,DataMultiSet.DataMultiSet):
        YY=_ComputeXSec_Data(data,method)
    elif  isinstance(data,dict):
        YY=_ComputeXSec_Point(data,method)
    else:
        raise Exception('ComputeXSec works only with DataSet or DataMultiSet or a point-like dictionary')
        
    return YY

def _ComputeXSec_Data(data,method="default"):
    """ Computes the cross-section values for the given data (matched to data)
    data can be: DataSet or DataMultiSet
    methods are:
        default = usual one
    """
    
    if method=="default":        
        if data.processType == "G2":
            XX=harpy.G2List([d["<x>"] for d in data.points],
                                [d["<Q>"] for d in data.points],
                                [d["process"] for d in data.points])
            
        elif data.processType=="D2":
            XX=harpy.D2List([d["<Q>"] for d in data.points],
                                [d["process"] for d in data.points])
    else:
        raise ValueError("Method is unknown.")
     
    YY=data.MatchWithData(XX)
    return YY

def _ComputeXSec_Point(p,method="default"):
    """ Computes the cross-section values for the given point
    methods are:
        default = usual one
    """
    
    if method=="default":        
        if p["type"] == "G2":
            XX1=harpy.G2List([p["<x>"]],[p["<Q>"]],[p["process"]])
            
        elif p["type"]=="D2":
            XX1=harpy.D2List([p["<Q>"]],[p["process"]])
        else:
            raise Exception('The dictionary is not a point.')
        
        XX=XX1[0]*p["thFactor"]
    else:
        raise ValueError("Method is unknown.")    
       
    return XX

def ComputeChi2(data,method="default"):
    """
    Computes the chi^2 values for the given data    
    data can be: DataSet or DataMultiSet

    Parameters
    ----------
    data : DataSet or DataMultiSet
        The data-set for which the chi2 computed
    method : string
        method for computation of xSec(see). The default is "default".

    Returns
    -------
    float , [float,float,...]
        The first number is chi^2 for total data-set
        The second array is the list of chi^2 for each experiment in the data set

    """
    YY=ComputeXSec(data,method)
    
    ZZ=data.chi2(YY)
    
    if isinstance(data,DataSet.DataSet):
        return ZZ, [ZZ]
    elif isinstance(data,DataMultiSet.DataMultiSet):
        return ZZ
    else:
        raise ValueError("data-argument maust be DataSet or DataMultiSet")

def PrintChi2Table(data,method="default",printSysShift=True,printDecomposedChi2=False):
    """
    Compute and print the values of chi^2/Npt for experiments

    Parameters
    ----------
    data : DataSet or DataMultiSet
        The data for which the computation is made
    method : string, optional
        The method of computation of xSection.
        The default is "default".
    printSysShift: bool, optional
        If True print the list of systematic shifts (determined by nuisance parameter)
        The default is True
    printDecomposedChi2: bool, optional
        If True the chi62 presented in the decomposed form
        The default is False

    Returns
    -------
    None.

    """
    import time
    
    startT=time.time()
    
    YY=ComputeXSec(data,method)
    
    ZZ=data.chi2(YY)
    
    if isinstance(data,DataSet.DataSet):
        chi2T, chi2Part = ZZ, [ZZ]
    else:
        chi2T, chi2Part = ZZ
    
    #chi2T, chi2Part = ComputeChi2(data,method)
    
    if printSysShift:
        shift=data.DetermineAvarageSystematicShift(YY)
    
    if printDecomposedChi2:
        decChi2=data.DecomposeChi2(YY)
    
    endT=time.time()
    
    
    
    if isinstance(data,DataSet.DataSet):
        maxLength=len(data.name)
    elif isinstance(data,DataMultiSet.DataMultiSet):
        maxLength=max([len(d.name) for d in data.sets])
    else:
        maxLength=10
    
    line="{:{width}}".format("name",width=maxLength)+' | '+'  N  '+' | '
    line2="{:-<{width}}".format("",width=maxLength)+'-|-'+'-----'+'-|-'
    if printDecomposedChi2:
        line+=' chiL^2/N '+' | '+' chiD^2/N '+' | '+' chi^2/N  '+' | '
        line2+='----------'+'-|-'+'----------'+'-|-'+'----------'+'-|-'
    else:
        line+=' chi^2/N  '+' | '
        line2+='----------'+'-|-'
    
    if printSysShift:
        line+='sys.shift%'+' | '        
        line2+='----------'+'-|-'
    print(line)
    print(line2)
   
    #Only one set in MultiSet of just Set
    if len(chi2Part)==1:
        line="{:{width}} | {:5d} |".format(data.name,data.numberOfPoints,width=maxLength)
        dataNumP=data.numberOfPoints
        if(dataNumP==0): dataNumP=1
        if printDecomposedChi2:
            line+=" {:10.3f} | {:10.3f} | {:10.3f} |".format(
                decChi2[0]/dataNumP,decChi2[1]/dataNumP,decChi2[2]/dataNumP)
        else:
            line+=" {:10.3f} |".format(chi2T/dataNumP)
        if printSysShift:
            line+=" {:10.3f} |".format(shift*100)
        print(line)
    else:
        for i in range(len(chi2Part)):
            line="{:{width}} | {:5d} |".format(data.sets[i].name,data.sets[i].numberOfPoints,width=maxLength)
            
            dataNumP=data.sets[i].numberOfPoints
            if(dataNumP==0): dataNumP=1
            
            if printDecomposedChi2:
                line+=" {:10.3f} | {:10.3f} | {:10.3f} |".format(
                    decChi2[i][0]/dataNumP,
                    decChi2[i][1]/dataNumP,
                    decChi2[i][2]/dataNumP)
            else:
                line+=" {:10.3f} |".format(chi2Part[i]/dataNumP)
            if printSysShift:
                line+=" {:10.3f} |".format(shift[i]*100)
            print(line)
            
        print(line2)
        line="{:{width}} | {:5d} |".format('Total',data.numberOfPoints, width=maxLength)
        if printDecomposedChi2:
            sumdChi=numpy.sum(decChi2, axis=0)
            line+=" {:10.3f} | {:10.3f} | {:10.3f} |".format(
                    sumdChi[0]/data.numberOfPoints,
                    sumdChi[1]/data.numberOfPoints,
                    sumdChi[2]/data.numberOfPoints)
        else:
            line+=" {:10.3f} |".format(chi2T/data.numberOfPoints)
        print(line)
        
    print("Computation time = ",endT-startT,' sec.')