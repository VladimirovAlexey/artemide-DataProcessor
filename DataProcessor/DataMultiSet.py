#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 17:28:13 2020

DataMultiSet is the collection of data sets.
The main purpuse of this class is to simplify the analysis of several DataSets

Fields specific for the DataMultiSet

processType     : DY or SIDIS. The sets should not mix. 
                the reason is that they are colled in parrallel, and thus
                it is better to keep them apart.
sets     : list of DataSet's included in the consideration
numberOfSets    : total number of sets in the multiSet
numberOfPoints  : total number of points in the multiSet

@author: vla18041
"""

from . import DataSet
import numpy

class DataMultiSet:
    """Collection of DataSet's
    """
    def __init__(self,name,dataSets):
        """
        

        Parameters
        ----------
        name : string
            The name for the collection of data sets
        dataSets : list of DataSet's
            These data sets will be added to the collection
        """     
        self.name=name
        self.sets=dataSets
        self.numberOfSets=len(dataSets)
        
        if all([self.sets[0].processType==ss.processType for ss in self.sets]):
            self.processType=self.sets[0].processType
        else:
            raise TypeError('processType must be same for all sets. \
                            Received value : {}'.format([ss.processType for ss in self.sets]))
        
        ## so to extract set 5, call [index1[5]:index2[5]]
        self._i1=[]#initial index of set in common list
        self._i2=[]#final index of set in common list
        i=0
        for dd in dataSets:
            self._i1.append(i)
            i+=dd.numberOfPoints
            self._i2.append(i)
            
        self.numberOfPoints=i
        
        #### I save list of points
        self.points=[]
                
        for dd in dataSets:
            for p in dd.points:
                self.points.append(p)                    
        
    def __repr__ (self):
        return "<DataMultiSet: %s with %s sets and %s points>" % (
                self.name, self.numberOfSets, self.numberOfPoints)
    
    def MatchWithData(self,theoryPrediction):
        """Reweight the data with proper factors.            
        """
        res=[]
        for i in range(self.numberOfSets):
            res.append(self.sets[i].MatchWithData(theoryPrediction[self._i1[i]:self._i2[i]]))
        return [item for sublist in res for item in sublist]
        
    def chi2(self,matchedTheory):
        """ Returns the chi2 values as total chi2, [partialchi2]
        """
        res=[]
        for i in range(self.numberOfSets):
            res.append(self.sets[i].chi2(matchedTheory[self._i1[i]:self._i2[i]]))
        return sum(res), res
    
    def DetermineSystematicShift(self,matchedTheory):
        """
        Determine the systematic shift by nuisance parameters evaluations.        

        Parameters
        ----------
        theoryPrediction : list of list of floats
            List theory predictions (matched) to be compared to the data

        Returns
        -------
        result : list of list of floats
            List of systematic shift point-per-point

        """
        res=[]
        for i in range(self.numberOfSets):
            res.append(self.sets[i].DetermineSystematicShift(matchedTheory[self._i1[i]:self._i2[i]]))
        return res
    
    def DetermineAvarageSystematicShift(self,matchedTheory):
        """
        Determine the mean systematic shift in % by nuisance parameters evaluations.

        Parameters
        ----------
        theoryPrediction : list of list of floats
             List theory predictions (matched) to be compared to the data

        Returns
        -------
        list of floats
            The avarage excess/deficit of the theory to data (due to systematic shift)

        """
        res=[]
        for i in range(self.numberOfSets):
            res.append(self.sets[i].DetermineAvarageSystematicShift(matchedTheory[self._i1[i]:self._i2[i]]))
        return res
    
    def DecomposeChi2(self,matchedTheory):
        """
        Determine the nuisance parameters and perform the decomposition of chi^2 to correlated and uncorrelated parts.

        Parameters
        ----------
        theoryPrediction : list of list of floats
             List theory predictions (matched) to be compared to the data
             
        Returns
        -------
        list of list of 3 floats
            chi2(uncorr), chi2(corr), chi2(total)

        """
        res=[]
        for i in range(self.numberOfSets):
            res.append(self.sets[i].DecomposeChi2(matchedTheory[self._i1[i]:self._i2[i]]))
        return res
    
    def FindBestNorm(self,theoryPrediction):
        """
        Evaluate the best common norm for the theory that minimizes the chi^2
        It is equl to n=(t V^{-1} xSec)/(t V^{-1} t), 
        where t is theory prediction, xSec is experimental values, and V is covariance matrix

        Parameters
        ----------
        theoryPrediction : list of floats
             List theory predictions (matched) to be compared to the data

        Returns
        -------
        Float, values of norm

        """
        res=[]
        for i in range(self.numberOfSets):
            res.append(self.sets[i].FindBestNorm(theoryPrediction[self._i1[i]:self._i2[i]]))
        return res
    
    def CutData(self,cutFunction,addName="",computeCovarianceMatrix = True):
        """
        Create an instance of MultiDataSet, which contains all sets after 
        application of cutFunction. New multisetset has name=name+addName

        Parameters
        ----------
        cutFunction : function with the interface f(Point1)=bool, Point2
            where bool states that the point2 should be included into the set.
                        
        addName : string, optional
            addendant to the name of the new sets (includind the subsets)
            Default values is ""
        
        computeCovarianceMatrix : bool, optional
            State is covariance matrix should be computed.
            Switching it False, can improve performace
            The default is True.

        Returns
        -------
        dNew : DataSet
            The cut DataSet

        """
        nSets=[s.CutData(cutFunction,addName,computeCovarianceMatrix) for s in self.sets]
        nSets= [x for x in nSets if x.numberOfPoints != 0]

        return DataMultiSet(self.name+addName,nSets)
    
    def GenerateReplica(self,includeNormInV=True):
        """Create a new multiset of pseudo-data that is replica of the given multiset.
        The reciept of generation is taken from [0808.1231] (see sec.2.4)
        Basically, it is p ->(1+r1 normErr) p(1+r2 pointErr), where 
        r1, and r2 are Gaussian randoms (r is common for the correlated error)
        
        includeNormInV: this is option that include the normalization error in V, or not.
        According to general theory it should be dropped,
        otherwice one get G. D’Agostini-bias on determination of parameters. 
        However, it is seutible only for the well-behaving data, which could affect normalization,
        In TMD case, often one have to include norm-error since it is only way to get meaningfull result.
        """
        newSets=[]
        for s in self.sets:
            newSets.append(s.GenerateReplica(includeNormInV))
            
        return DataMultiSet(self.name+'(rep)',newSets)

