#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 20:29:02 2020

DataSet is a class that contains a list of data points
It is supposed that these points have some common origin (e.g. experiment)
and thus, share correlated errors.

The points included in the list MUST have same:
type, length of error-list (both correlated and uncorrelated)

Upon finalization the object computes covariance matrix and other stuff. 
And afterwords can compute chi^2, etc.

The fields specific for the DataSet are

processType = DY or SIDIS 
points      = list of point-dictionaries
name        = the name for the set

numberOfPoints = total number of points in the set
numberOfUncorrErr = number of uncorrelated errors for each point
numberOfCorrErr = number of correlated errors for each point
numberOfNormErr = number of normalization errors for each point
normErr         = list of normalization errors 
                normalization errors are specific to the set, 
                and thus should not be placed into points.
                In relative values (0.01=1%)

isNormalized  = the set is considered normalized, 
                i.e. the theory computation will be normalized to the data
normalizationMethod = the method to normalize hte data to the theory

V             = covariance matrix
invV          = inverse covariance matrix
A             = matrix requared to compute systematic shifts
invA          = inverse A
....

@author: vla18041
"""

import numpy
from . import Point

class DataSet:
    """Collection of points with the common correlation errors.
    """
    
    def __init__(self,name,processType):
        """ name=short name
        """
        if not(processType=="DY" or processType=="SIDIS"):
            raise Exception('processType must be DY or SIDIS. Received value : {}'.format(processType))
                
        self.processType=processType
        self.name=name
        self.comment=" "
        self.reference=" "
        
        self.points=[]
        self.numberOfPoints=0
        
        ## option that set is self normalized 
        ## the theory would be normalized to experiment
        self.isNormalized=False
        ## the normalization method
        self.normalizationMethod="integral"
        ## intrinsic variable for normalization (computed during finalization)
        self._normExp=1.
        
        ## list of normalization errors
        ## normalization uncertanties are RELATIVE (0.01=1%)
        self.normErr=[]
        
        ## Numbers of errors of points (filled with finalization)
        self.numOfUncorrErr=-1
        self.numOfCorrErr=-1
        self.numOfNormErr=-1
        
        self.V=[] # covariance matrix 
        self.invV=[] # inverse covariance matrix
        self.matrixA=[] # matrix of systemtic shifts
        self.matrixAinverse=[]
        
        
    def __repr__ (self):
        return "<DataSet: %s with %s points>" % (self.name, self.numberOfPoints)
    
    def AddPoint(self,p):
        """Adds a point to the list of points.
        """
        ### Since I do not what to modify original points I must (deep)copy them
        import copy        
        pCopy = copy.deepcopy(p) 
        
        if Point.FinalizePoint(pCopy):
            
            if pCopy["type"]==self.processType:
                self.points.append(pCopy)
                self.numberOfPoints+=1
                
            else:
                print("AddPoint: Point has wrong type")
        else:
            print("AddPoint: Point is not appended")
    
    def FinalizeSet(self,computeCovarianceMatrix = True, silent = False):
        """
        Check internal consitency of the set, fill nessecary fields, 
        compute covariance matrices.
        
        Call this function only after all points are added.        

        Parameters
        ----------
        computeCovarianceMatrix : bool, optional
            State is covariance matrix should be computed.
            Switching it False, can improve performace
            The default is True.
        silent : bool, optional
            If true suppress the final message
            The default is False.

        Returns
        -------
        None.

        """        
        ### 0) Check presentce of points
        if self.numberOfPoints != len(self.points):
            print('DataSet.FinalizeSet: Mismatch in the number of points. Data set: ',self.name)
            self.numberOfPoints = len(self.points)
        if self.numberOfPoints <= 0:
            print('DataSet.Finalize: ERROR: the number of points ',self.numberOfPoints,' non-positive.')
            print('Unable to finalize data set :',self.name)
        
        ### 1) fill fields with numbers of errors + fill list of corr.error for each point
        self.numOfUncorrErr=len(self.points[0]["uncorrErr"])
        self.numOfCorrErr=len(self.points[0]["corrErr"])
        self.numOfNormErr=len(self.normErr)
        
        for p in self.points:
            if(len(p["uncorrErr"]) != self.numOfUncorrErr):
                print('DataSet.Finalize: WARNING: Mismatch in the leght of uncorrErr for point ',p.identifier,' in data set: ',self.name)
            if(len(p["corrErr"]) != self.numOfCorrErr):
                print('DataSet.Finalize:: Mismatch in the leght of corrErr for point',p.identifier,' in data set: ',self.name)
            
        
        ### 2) Calculate covariace matrices
        if computeCovarianceMatrix:
            self.CalculateV()
            self.invV=numpy.linalg.inv(self.V)
            self.CalculateA()
            self.matrixAinverse=numpy.linalg.inv(self.matrixA)
        

        ### 3) Calculate some working variables
        self._normExp=sum([(p["qT"][1]-p["qT"][0])*p["xSec"] for p in self.points])
        #self._xSecString=[p.xSec for p in self.points]     
        
        if not silent:
            print('Set ',self.name,' finalized. (',self.numberOfPoints,' points)')

    def CutData(self,cutFunction,addName):
        """ Create an instance of DataSet, which contains all points after 
            application of cutFunction. New set has name=name+addName
            
            cutFunction has interface
            f(Point1)=bool, Point2
            where bool states that the point2 should be included into the set.
        """
        import copy
        
        # create the set
        dNew=DataSet(self.name+addName,self.processType)
        
        # copy main fields
        dNew.comment=copy.copy(self.comment)
        dNew.reference=copy.copy(self.reference)
        dNew.isNormalized=self.isNormalized
        dNew.normErr=copy.deepcopy(self.normErr)
        dNew.normalizationMethod=copy.copy(self.normalizationMethod)
        
        # populate with points
        for p in self.points:
            pPass, pAdd=cutFunction(p)
            if pPass:
                dNew.AddPoint(pAdd)
        
        dNew.FinalizeSet()    
        
        return dNew       
                
        
    
    def CalculateV(self):
        """ Calculate the covariance matrix.
        
        The covariance matrix defined as, V[i,j] where
        V[i,i]=sum uncorrErr^2
        V[i,j]=sum corrErr[i] corrErr[j]
        """
        self.V=numpy.zeros(shape=(self.numberOfPoints,self.numberOfPoints))
        
        for i in range(self.numberOfPoints):
            for j in range(self.numberOfPoints):
                ## diagonal term
                if i==j:
                    a=numpy.array(self.points[i]["uncorrErr"])                    
                    self.V[i][j]+=numpy.sum(a**2)
                
                ## correlated errors 
                a=numpy.array(self.points[i]["corrErr"])
                b=numpy.array(self.points[j]["corrErr"])
                
                self.V[i][j]+=numpy.sum(a*b)
                
                ## normalization errors
                a=numpy.array(self.normErr)*self.points[i]["xSec"]
                b=numpy.array(self.normErr)*self.points[j]["xSec"]
                
                self.V[i][j]+=numpy.sum(a*b)
        
    def CalculateA(self):
        """ Evaluate the matrix A which is needed for the estimation of systematic shifts.
        
        The matrix A defined as 
        A[a,b]=delta[a,b]+sum(corrErr[a,i]corrErr[b,i]/(uncorrErr[i]^2), i in number of points)
        """
        ## determine the number of correlated errors
        sc=self.numOfCorrErr+self.numOfNormErr
        self.matrixA=numpy.identity(sc)
        
        ## building list of all correlated errors.
        cErr=[]
        for p in self.points:
            a=numpy.array(p["uncorrErr"])
            cErr.append(numpy.concatenate(
                    (numpy.array(p["corrErr"]),numpy.array(self.normErr)*p["xSec"]))
                    /numpy.sqrt(numpy.sum(a**2)))
        
         ## add contribution of a point to matrix A
        for i in range(sc):
            for j in range(sc):
                for c in cErr:          
                    self.matrixA[i][j]+=c[i]*c[j]
                    
    def chi2(self,theoryPrediction):
        """ Evaluate chi^2=x.V^{-1}.x, where x is theory prediction
        """
        diffX=[(theoryPrediction[i]-self.points[i]["xSec"]) for i in range(self.numberOfPoints)]
        return numpy.matmul(diffX,numpy.matmul(self.invV,diffX))
    
    
    def MatchWithData(self,theoryPrediction):
        """Adjust the input theory vector to the experimental data.
        
        Multiply by the proper translation factor and normalize if needed.
        """
        res=[]
        norm=1
        ## weight theory with proper factor
        for i in range(self.numberOfPoints):
            res.append(theoryPrediction[i]*self.points[i]["thFactor"])
        
        ##if nessecary normalize
        if self.isNormalized:
            if self.normalizationMethod=="integral":
                ### normalization by bin-by-bin area
                normTh=sum([(self.points[i]["qT"][1]-self.points[i]["qT"][0])*res[i] for i in range(self.numberOfPoints)])
                norm=self._normExp/normTh       
                return [v*norm for v in res]
        else:
            return res  
        
    def GenerateReplica(self,includeNormInV=True):
        """Create a new set of pseudo-data that is replica of the given set.
        The reciept of generation is taken from [0808.1231] (see sec.2.4)
        Basically, it is p ->(1+r1 normErr) p(1+r2 pointErr), where 
        r1, and r2 are Gaussian randoms (r is common for the correlated error)
        
        includeNormInV: this is option that include the normalization error in V, or not.
        According to general theory it should be dropped,
        otherwice one get G. D’Agostini-bias on determination of parameters. 
        However, it is seutible only for the well-behaving data, which could affect normalization,
        In TMD case, often one have to include norm-error since it is only way to get meaningfull result.
        """
        import copy
        
        # create the set
        dNew=DataSet(self.name+'(rep)',self.processType)
        
        # copy main fields
        dNew.comment=copy.copy(self.comment+'(replica)')
        dNew.reference=copy.copy(self.reference)
        dNew.isNormalized=self.isNormalized
        dNew.normalizationMethod=copy.copy(self.normalizationMethod)
        if includeNormInV:
            dNew.normErr=copy.deepcopy(self.normErr)
        else:
            dNew.normErr=[]

        ## RND for correlated errors
        corrRND=[numpy.random.normal() for i in range(self.numOfCorrErr)]
        
        # this is common rescaling factor due to norm uncerantity
        # it is 1+ rand * normErr, because normErr are given relative to xSec (in %)
        resFactor=1.
        for err in self.normErr:
            resFactor+=numpy.random.normal()*err
        
        # populate with points
        # the only fields to update are xSec and errors
        for p in self.points:
            # copy the point
            pNew=copy.deepcopy(p)
            # update xSec value
            # (old+uncorr*RND+corr*globalRND)*(1+norm)
            ## uncorrelated part
            for err in p["uncorrErr"]:
                pNew["xSec"]+=numpy.random.normal()*err
            # correlated part
            for i in range(self.numOfCorrErr):
                pNew["xSec"]+=corrRND[i]*p["corrErr"][i]
            # multiply by norm-error-distirbution
            pNew["xSec"]=pNew["xSec"]*resFactor
            
            ### if initial or final xSec<=0 we set the result to zero.
            # such situation happens for some very bad sets of data
            if(pNew["xSec"]<0. or p["xSec"]<=0.0):
                pNew["xSec"]=0.0
            
            ##the point-errors should be rescaled by normalization factor
            pNew["uncorrErr"]=[]
            for err in p["uncorrErr"]:
                pNew["uncorrErr"].append(err*resFactor)
            pNew["corrErr"]=[]
            for err in p["corrErr"]:
                pNew["corrErr"].append(err*resFactor)
            #
            dNew.AddPoint(pNew)
        
        dNew.FinalizeSet()   
        return dNew
    
    def SaveToCSV(self,path):
        """
        Save the current data set in the CSV-format

        Parameters
        ----------
        path : string
            Path to the CSV

        Returns
        -------
        None.

        """
        
        self.FinalizeSet(computeCovarianceMatrix = False)
        
        file=open(path,"w")
        
        ### general information
        
        file.write("Name,"+self.name+"\n")
        file.write("Comment,"+self.comment+"\n")
        file.write("Reference,"+self.reference+"\n")
        file.write("Process type,"+self.processType+"\n")
        file.write("Number of points,"+str(self.numberOfPoints)+"\n")
        file.write("Number of uncorr.errors,"+str(self.numOfUncorrErr)+"\n")
        file.write("Number of corr.errors,"+str(self.numOfCorrErr)+"\n")
        file.write("Number of norm.errors,"+str(self.numOfNormErr)+"\n")
        file.write("List of norm.errors (relative),"+",".join([str(el) for el in self.normErr])+"\n")
        file.write("Total cross-section nomalized,"+str(self.isNormalized)+"\n")
        file.write("List of points \n")
        
        ### Points table header
        
        file.write("Point id,process id,s[GeV^2],<Q>[GeV],Qmin[GeV],Qmax[GeV],")
        if self.processType=="DY":
            file.write("<y>,yMin,yMax,<qT>[GeV],qTMin[GeV],qTMax[GeV],")
        elif self.processType=="SIDIS":
            file.write("<x>,xMin,xMax,<z>,zMin,zMax,<pT>[GeV],pTMin[GeV],pTMax[GeV],")
        else:
            raise Exception("Unknown type of process")
            
        file.write("xSec,")
            
        for i in range(self.numOfUncorrErr):
            file.write("Uncorr.Err."+str(i)+",")
            
        for i in range(self.numOfCorrErr):
            file.write("Corr.Err."+str(i)+",")
            
        file.write("Th.Factor,FiducialCuts,")
        
        if self.processType=="DY":
            file.write("kCut1[GeV],kCut2[GeV],etaMin,etaMax")
        elif self.processType=="SIDIS":
            file.write("yMin,yMax,W2min[GeV^2],W2max[GeV^2],TargetMass[GeV],ProductMass[GeV]")
        else:
            raise Exception("Unknown type of process")
        
        file.write("\n")
        
        ### Points table 
        
        for p in self.points:
            file.write(p["id"])
            file.write(",")
            file.write(str(p["process"]).replace(","," -").replace("[","").replace("]",""))
            file.write(",")
            file.write(str(p["s"]))
            file.write(",")
            file.write(str(p["<Q>"]))
            file.write(",")
            file.write(str(p["Q"][0]))
            file.write(",")
            file.write(str(p["Q"][1]))
            file.write(",")
            if self.processType=="DY":
                file.write(str(p["<y>"]))
                file.write(",")
                file.write(str(p["y"][0]))
                file.write(",")
                file.write(str(p["y"][1]))
                file.write(",")
                file.write(str(p["<qT>"]))
                file.write(",")
                file.write(str(p["qT"][0]))
                file.write(",")
                file.write(str(p["qT"][1]))
            elif self.processType=="SIDIS":
                file.write(str(p["<x>"]))
                file.write(",")
                file.write(str(p["x"][0]))
                file.write(",")
                file.write(str(p["x"][1]))
                file.write(",")
                file.write(str(p["<z>"]))
                file.write(",")
                file.write(str(p["z"][0]))
                file.write(",")
                file.write(str(p["z"][1]))
                file.write(",")
                file.write(str(p["<pT>"]))
                file.write(",")
                file.write(str(p["pT"][0]))
                file.write(",")
                file.write(str(p["pT"][1]))
            else:
                raise Exception("Unknown type of process")  
                
            file.write(",")
                        
            file.write(str(p["xSec"]))
            file.write(",")
            
            for i in range(self.numOfUncorrErr):
                file.write(str(p["uncorrErr"][i]))
                file.write(",")
            
            for i in range(self.numOfCorrErr):
                file.write(str(p["corrErr"][i]))
                file.write(",")
                
            file.write(str(p["thFactor"]))
            file.write(",")
            file.write(str(p["includeCuts"]))
            file.write(",")
            file.write(str(p["cutParams"][0]))
            file.write(",")
            file.write(str(p["cutParams"][1]))
            file.write(",")
            file.write(str(p["cutParams"][2]))
            file.write(",")
            file.write(str(p["cutParams"][3]))            
            
            if self.processType=="SIDIS":
                file.write(",")
                file.write(str(p["M_target"]))
                file.write(",")
                file.write(str(p["M_product"]))
                
            file.write("\n")
        
        file.close()
        
def LoadCSV(path):
    """
    Reads the CSV file and attemt to create a DataSet from it

    Parameters
    ----------
    path : string
        Path to the file

    Returns
    -------
    The dataSet that was restored from file.

    """
    file=open(path,"r")
    
    ## read header of
    line=file.readline()    
    line=line.split(",")
    name=line[1].replace("\n","")
    
    line=file.readline()    
    line=line.split(",")
    comment=line[1].replace("\n","")
    
    
    line=file.readline()    
    line=line.split(",")
    ref=line[1].replace("\n","")
    
    line=file.readline()    
    line=line.split(",")
    processType=(line[1]).replace("\n","")
    
    line=file.readline()    
    line=line.split(",")   
    nPoints=int(line[1])
    
    line=file.readline()    
    line=line.split(",")
    nUncorrErr=int(line[1])
    
    line=file.readline()    
    line=line.split(",")
    nCorrErr=int(line[1])
    
    ### create a Data set
    dSet=DataSet(name, processType)
    
    dSet.comment=comment
    dSet.reference=ref
    
    ###    continue to read header
    line=file.readline()    
    line=line.split(",")
    nNormErr=int(line[1])
    
    line=file.readline()    
    line=line.split(",")
    for i in range(nNormErr):
        dSet.normErr.append(float(line[i+1]))
        
    line=file.readline()    
    line=line.split(",")
    dSet.isNormalized=(line[1].replace("\n","")=="True")
    
    ### skip head for points
    line=file.readline()    
    line=file.readline()
    for i in range(nPoints):
        
        k=0### number of the last line in the block
        
        line=file.readline()
        line=line.split(",")
        if processType=="DY":
            p=Point.CreateDYPoint(line[0])
        elif processType=="SIDIS":
            p=Point.CreateSIDISPoint(line[0])
        else:
            raise Exception("Unknown process")
        
        p["process"]=[int(ll) for ll in line[1].split("-")]
        p["s"]=float(line[2])
        p["<Q>"]=float(line[3])
        p["Q"]=[float(line[4]),float(line[5])]
        k=k+5
        
        if processType=="DY":
            p["<y>"]=float(line[k+1])
            p["y"]=[float(line[k+2]),float(line[k+3])]
            p["<qT>"]=float(line[k+4])
            p["qT"]=[float(line[k+5]),float(line[k+6])]
            k=k+6
            
        elif processType=="SIDIS":
            p["<x>"]=float(line[k+1])
            p["x"]=[float(line[k+2]),float(line[k+3])]
            p["<z>"]=float(line[k+4])
            p["z"]=[float(line[k+5]),float(line[k+6])]
            p["<pT>"]=float(line[k+7])
            p["pT"]=[float(line[k+8]),float(line[k+9])]
            k=k+9
            
        else:
            raise Exception("Unknown process")
        
        p["xSec"]=float(line[k+1])
        k=k+1
        
        for j in range(nUncorrErr):
            p["uncorrErr"].append(float(line[k+1]))
            k=k+1
        for j in range(nCorrErr):
            p["corrErr"].append(float(line[k+1]))
            k=k+1
        
        p["thFactor"]=float(line[k+1])
        p["includeCuts"]=((line[k+2].replace("\n",""))=="True")
        p["cutParams"]=[float(line[k+3]),float(line[k+4]),float(line[k+5]),float(line[k+6])]
        k=k+6
        
        if processType=="SIDIS":
            p["M_target"]=float(line[k+1])
            p["M_product"]=float(line[k+2])
            k=k+2
            
        dSet.AddPoint(p)
        
    file.close()
    
    dSet.FinalizeSet(computeCovarianceMatrix=False,silent=True)
    
    return dSet