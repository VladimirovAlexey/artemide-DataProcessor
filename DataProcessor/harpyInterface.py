#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 14:52:27 2020

@author: vla18041
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
        binless = evaluated with avarage values of bins only. Multiplied by area of the bin
        central = evaluated with avarage values of bins only. No further modifications.
        approximate = evaluated with approximate procedure (very few points per bin integration)
    """
    
    if method=="default":        
        if data.processType == "DY":
            XX=harpy.DY.xSecList([d["process"] for d in data.points],
                                    [d["s"] for d in data.points],
                                    [d["qT"] for d in data.points],
                                    [d["Q"] for d in data.points],
                                    [d["y"] for d in data.points],
                                    [d["includeCuts"] for d in data.points],
                                    [d["cutParams"] for d in data.points])
            
        elif data.processType=="SIDIS":
            XX=harpy.SIDIS.xSecList([d["process"] for d in data.points],
                                    [d["s"] for d in data.points],
                                    [d["pT"] for d in data.points],
                                    [d["z"] for d in data.points],
                                    [d["x"] for d in data.points],
                                    [d["Q"] for d in data.points],
                                    [d["includeCuts"] for d in data.points],
                                    [d["cutParams"] for d in data.points],
                                    [[d["M_target"],d["M_product"]] for d in data.points])
    elif method=="binless":
        if data.processType == "DY":
            XX1=harpy.DY.xSecListBINLESS([d["process"] for d in data.points],
                                    [d["s"] for d in data.points],
                                    [d["<qT>"] for d in data.points],
                                    [d["<Q>"] for d in data.points],
                                    [d["<y>"] for d in data.points],
                                    [d["includeCuts"] for d in data.points],
                                    [d["cutParams"] for d in data.points])
            XX=[]
            ## I weight the binless computation by the width of the bin
            for i in range(len(XX1)):
                XX.append(XX1[i]*
                          (data.points[i]["Q"][1]**2-data.points[i]["Q"][0]**2)*
                          (data.points[i]["qT"][1]**2-data.points[i]["qT"][0]**2)*
                          (data.points[i]["y"][1]-data.points[i]["y"][0]))
            
        elif data.processType=="SIDIS":
            XX1=harpy.SIDIS.xSecListBINLESS([d["process"] for d in data.points],
                                    [d["s"] for d in data.points],
                                    [d["<pT>"] for d in data.points],
                                    [d["<z>"] for d in data.points],
                                    [d["<x>"] for d in data.points],
                                    [d["<Q>"] for d in data.points],
                                    [[d["M_target"],d["M_product"]] for d in data.points])
            XX=[]
            for i in range(len(XX1)):
                XX.append(XX1[i]*
                          (data.points[i]["Q"][1]**2-data.points[i]["Q"][0]**2)*
                          (data.points[i]["pT"][1]**2-data.points[i]["pT"][0]**2)*
                          (data.points[i]["x"][1]-data.points[i]["x"][0])*
                          (data.points[i]["z"][1]-data.points[i]["z"][0]))
    elif method=="central":
        if data.processType == "DY":
            XX=harpy.DY.xSecListBINLESS([d["process"] for d in data.points],
                                    [d["s"] for d in data.points],
                                    [d["<qT>"] for d in data.points],
                                    [d["<Q>"] for d in data.points],
                                    [d["<y>"] for d in data.points],
                                    [d["includeCuts"] for d in data.points],
                                    [d["cutParams"] for d in data.points])
            
        elif data.processType=="SIDIS":
            XX=harpy.SIDIS.xSecListBINLESS([d["process"] for d in data.points],
                                    [d["s"] for d in data.points],
                                    [d["<pT>"] for d in data.points],
                                    [d["<z>"] for d in data.points],
                                    [d["<x>"] for d in data.points],
                                    [d["<Q>"] for d in data.points],
                                    [[d["M_target"],d["M_product"]] for d in data.points])
    elif method=="approximate":
        if data.processType == "DY":
            XX=harpy.DY.xSecListAPROX([d["process"] for d in data.points],
                                    [d["s"] for d in data.points],
                                    [d["qT"] for d in data.points],
                                    [d["Q"] for d in data.points],
                                    [d["y"] for d in data.points],
                                    [d["includeCuts"] for d in data.points],
                                    [d["cutParams"] for d in data.points])
            
        elif data.processType=="SIDIS":
            raise ValueError("Method <approximate> is not yet implemented for SIDIS.")
    else:
        raise ValueError("Method is unknown.")
     
    YY=data.MatchWithData(XX)
    return YY

def _ComputeXSec_Point(p,method="default"):
    """ Computes the cross-section values for the given point
    methods are:
        default = usual one
        binless = evaluated with avarage values of bins only. Multiplied by area of the bin
        central = evaluated with avarage values of bins only. No further modifications.
    """
    
    if method=="default":        
        if p["type"] == "DY":
            XX1=harpy.DY.xSecList([p["process"]],[p["s"]],[p["qT"]],[p["Q"]],
                                    [p["y"]],[p["includeCuts"]],[p["cutParams"]])
            
        elif p["type"]=="SIDIS":
            XX1=harpy.SIDIS.xSecList([p["process"]],[p["s"]],[p["pT"]],[p["z"]],
                                    [p["x"]],[p["Q"]],[p["includeCuts"]],
                                    [p["cutParams"]],[[p["M_target"],p["M_product"]]])
        else:
            raise Exception('The dictionary is not a point.')
        
        XX=XX1[0]*p["thFactor"]
    
    elif method=="approximate":
        if p["type"] == "DY":
            XX1=harpy.DY.xSecListAPROX([p["process"]],[p["s"]],[p["<qT>"]],[p["<Q>"]],
                                    [p["<y>"]],[p["includeCuts"]],[p["cutParams"]])
        elif p["type"]=="SIDIS":
            raise ValueError("Method <approximate> is not yet implemented for SIDIS.")
        else:
            raise Exception('The dictionary is not a point.')
        
        XX=XX1[0]*p["thFactor"]
        
    elif method=="binless":
        if p["type"] == "DY":
            XX1=harpy.DY.xSecListBINLESS([p["process"]],[p["s"]],[p["<qT>"]],[p["<Q>"]],
                                    [p["<y>"]],[p["includeCuts"]],[p["cutParams"]])
            ### Matching with data (note that the normalization is ignored)
            XX=XX1[0]*(p["Q"][1]**2-p["Q"][0]**2)*(p["qT"][1]**2-p["qT"][0]**2)* \
                          (p["y"][1]-p["y"][0])*p["thFactor"]
            
        elif p["type"]=="SIDIS":
            XX1=harpy.SIDIS.xSecListBINLESS([p["process"]],[p["s"]],[p["<pT>"]],[p["<z>"]],
                                    [p["<x>"]],[p["<Q>"]],[[p["M_target"],p["M_product"]]])            

            ### Matching with data (note that the normalization is ignored)
            XX=XX1[0]*(p["Q"][1]**2-p["Q"][0]**2)*(p["pT"][1]**2-p["pT"][0]**2)* \
                (p["x"][1]-p["x"][0])*(p["z"][1]-p["z"][0])*p["thFactor"]
        else:
            raise Exception('The dictionary is not a point.')
    elif method=="central":
        if p["type"] == "DY":
            XX1=harpy.DY.xSecListBINLESS([p["process"]],[p["s"]],[p["<qT>"]],[p["<Q>"]],
                                    [p["<y>"]],[p["includeCuts"]],[p["cutParams"]])
            ### Matching with data (note that the normalization is ignored)
            XX=XX1[0]*p["thFactor"]
            
        elif p["type"]=="SIDIS":
            XX1=harpy.SIDIS.xSecListBINLESS([p["process"]],[p["s"]],[p["<pT>"]],[p["<z>"]],
                                    [p["<x>"]],[p["<Q>"]],[[p["M_target"],p["M_product"]]])            

            ### Matching with data (note that the normalization is ignored)
            XX=XX1[0]*p["thFactor"]
        else:
            raise Exception('The dictionary is not a point.')
    
       
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
    
def PrintPerPointContribution(data,method="default",output="id",minValue=0.):
    """
    Prints a table of contribution to chi^2 per point (only diagonal unceranties are included)

    Parameters
    ----------
    data : DataSet or DataMultiSet
        DESCRIPTION.
    method : TYPE, optional
        DESCRIPTION. The default is "default".
    output : string, split by commas
        What to output, can contain, x,<x>,<Q>,del,id. The default is "id". (TO BE ADDED MORE)
        
    minValue : float
        The value of chi^2 less that which the points are nor shown

    Returns
    -------
    None.

    """
    outputlist=output.split(",")
    
    line=""
    line2=""
    if("id" in outputlist):    
        idLength=max([len(p["id"]) for p in data.points])
        if(idLength<10): idLength=10
        line+="{:{width}}".format("name",width=idLength)+' | '
        line2+="{:-<{width}}".format("",width=idLength)+'-|-'
    else:
        line+=" "
        line2+="-"
    if("x" in outputlist):
        ### determining size for x variables
        xff=min([p["x"][0] for p in data.points])
        if(xff<0.0001): xW=5
        elif(xff<0.001): xW=4
        elif(xff<0.01): xW=3
        else: xW=2
        
        ## total size of x-cell
        #[0.xW,0.xW]=2xW+7
        xLength=2*xW+7
        line+="{:{width}}".format(" x",width=xLength)+' | '
        line2+="{:-<{width}}".format("",width=xLength)+'-|-'
    if("<x>" in outputlist):        
        ### determining size for x variables
        xff=min([p["<x>"] for p in data.points])
        if(xff<0.0001): xMW=5
        elif(xff<0.001): xMW=4
        elif(xff<0.01): xMW=3
        else: xMW=2
        ## total size of x-cell
        #0.xMW=xW+2
        xMLength=xMW+2
        line+="{:{width}}".format(" <x>",width=xMLength)+' | '
        line2+="{:-<{width}}".format("",width=xMLength)+'-|-'
    if("<Q>" in outputlist):        
        ### determining size for x variables
        xff=max([p["<Q>"] for p in data.points])
        if(xff>100): QMW=3
        elif(xff<10): QMW=2
        else: QMW=1
        ## total size of x-cell
        #0.xMW=xW+2
        QMLength=QMW+2
        line+="{:{width}}".format(" <Q>",width=QMLength)+' | '
        line2+="{:-<{width}}".format("",width=QMLength)+'-|-'
    if("del" in outputlist):        
        ### determining size for del variables {:6.3f}
        delLength=6 
        line+="{:{width}}".format(" delta",width=delLength)+' | '
        line2+="{:-<{width}}".format("",width=delLength)+'-|-'
        
    
    line+=' dChi^2   '+' | '
    line2+='----------'+'-|-'

    print(line)
    print(line2)
    
    def PrintP(p,chi2):
        line=""
        if("id" in outputlist):
            line+="{:{width}} |".format(p["id"],width=idLength)
        if("x" in outputlist):
            line+=" [{:{p1}.{p2}f},{:{p1}.{p2}f}] |".format(p["x"][0],p["x"][1],p1=xW+2,p2=xW)
        if("<x>" in outputlist):            
            line+=" {:{p1}.{p2}f} |".format(p["<x>"],p1=xMW+2,p2=xMW)
        if("<Q>" in outputlist):            
            line+=" {:{p1}.{p2}f} |".format(p["<Q>"],p1=QMW+1,p2=4-QMW)
        if("del" in outputlist):
            if(p["type"]=="DY"):
                ddd=p["<qT>"]/p["<Q>"]
            elif(p["type"]=="SIDIS"):
                ddd=p["<pT>"]/p["<Q>"]/p["<z>"]
            else: raise ValueError("Unknown process")
            line+=" {:6.3f} |".format(ddd)
        line+=" {:10.3f} |".format(chi2)
        print(line)
    
    if(isinstance(data, DataSet.DataSet)):
        YY=ComputeXSec(data,method)
        ZZ=[(data.points[i]["xSec"]-YY[i])**2/(numpy.sum(numpy.array(data.points[i]["uncorrErr"])**2)) for i in range(len(YY))]
        for i in range(len(YY)):
            PrintP(data.points[i],ZZ[i])
        print(line2)
        
    elif(isinstance(data, DataMultiSet.DataMultiSet)):
        for ss in data.sets:
            YY=ComputeXSec(ss,method)
            ZZ=[(ss.points[i]["xSec"]-YY[i])**2/(numpy.sum(numpy.array(ss.points[i]["uncorrErr"])**2)) for i in range(len(YY))]
            
            print(ss.name)
            print(line2)
            
            npoints=0
            for i in range(len(YY)):
                if(ZZ[i]>minValue):
                    PrintP(ss.points[i],ZZ[i])
                    npoints+=1
            if(npoints==0):
                print("All points have chi^2 < "+str(minValue))
            print(line2)