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
normalizationMethod = the method to normalize the data to the theory
                currently the following methods of normalization are realised:
                    'integral'  normalize to the partial (bin-by-bin) integral of the data
                    'bestChi2'  normalize to such that the chi^2 is minimal

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
        
        ## flag that signalize the emptiness of the DataSet
        self.isEmpty=True
        
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
        self.L=[] # Cholesky decomposition of covariant matrix
        self.matrixA=[] # matrix of systemtic shifts
        self.matrixAinverse=[]
        self.isVdiagonal=False ## If V is diagonal (i.e. not correlated errors)
        
        # the list of list correlated errors for each point
        # contains point-to-point correlations and normalization correlations
        self._listOfCorrErrors=numpy.array([])
        # the list of uncorrelated errors-squared
        self._listOfVariances=numpy.array([])
        
        
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
                self.isEmpty=False
                
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
        if self.numberOfPoints == 0:
            print('DataSet.Finalize: the number of points in null.')
            print('The DataSet '+self.name+' is empty.')
            self.isEmpty=True
            return None
        
        
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
            if(self.numOfCorrErr==0 and self.numOfNormErr==0) :
                self.isVdiagonal=True
            else:
                self.isVdiagonal=False
            self._computeListOfCorrErrors()
            self._computeListOfVariances()
            self._CalculateV()
            self._CholeskyDecompositionForV()
            self._CalculateA()
            self.matrixAinverse=numpy.linalg.inv(self.matrixA)
        

        ### 3) Calculate some working variables
        self._normExp=sum([(p["qT"][1]-p["qT"][0])*p["xSec"] for p in self.points])
        #self._xSecString=[p.xSec for p in self.points]     
        
        if not silent:
            print('Set ',self.name,' finalized. (',self.numberOfPoints,' points)')

    def CutData(self,cutFunction,addName="",computeCovarianceMatrix = True):
        """
        Create an instance of DataSet, which contains all points of the original one after 
        application of cutFunction. New set has name=name+addName

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
            pPass, pAdd=cutFunction(copy.deepcopy(p))
            if pPass:
                dNew.AddPoint(pAdd)
        
        dNew.FinalizeSet(computeCovarianceMatrix)    
        
        return dNew       
                
        
    
    def _CalculateV(self):
        """ Calculate the covariance matrix.
        
        The covariance matrix defined as, V[i,j] where
        V[i,i]=sum uncorrErr^2
        V[i,j]=sum corrErr[i] corrErr[j]
        """
        self.V=numpy.zeros(shape=(self.numberOfPoints,self.numberOfPoints))
        
        for i in range(self.numberOfPoints):
            for j in range(self.numberOfPoints):
                ## diagonal term is uncorrelated error SQUARED, already precalculated
                if i==j:                    
                    self.V[i][j]+=self._listOfVariances[i]
                ## off-diagonal term is correlated errors, already precalculated
                self.V[i][j]+=numpy.dot(self._listOfCorrErrors[i],self._listOfCorrErrors[j])
       
     
    def _computeListOfCorrErrors(self):
        """
        Compute the list of correlated errors for each point. 
        It is equals to [corrErr]+[xSec*normErr]
        """
        dummy=numpy.zeros((self.numberOfPoints,self.numOfCorrErr+self.numOfNormErr))
        
        for i in range(self.numberOfPoints):
            for j in range(self.numOfCorrErr):
                dummy[i][j]=self.points[i]["corrErr"][j]
            for j in range(self.numOfNormErr):
                dummy[i][self.numOfCorrErr+j]=self.points[i]["xSec"]*self.normErr[j]
        
        self._listOfCorrErrors=dummy
            
    def _computeListOfVariances(self):
        """
        Compute the list of uncorrelated errors SQUARED for each point.
        """
        dummy=numpy.zeros(self.numberOfPoints)
        
        for i in range(self.numberOfPoints):            
            for err in self.points[i]["uncorrErr"]:
                dummy[i]+=err**2
            
        self._listOfVariances=dummy
    
    def _CholeskyDecompositionForV(self):
        """
        Compute Cholesky decomposition for covariance matrix

        Returns
        -------
        None.

        """
        #### if the V is diagonal then L=sqrt(V)
        if(self.isVdiagonal):
            self.L=numpy.sqrt(self.V)
        else:
            ### code is taken from Rosseta project
            self.L = [[0.0] * len(self.V) for _ in range(len(self.V))]
            for i, (Ai, Li) in enumerate(zip(self.V, self.L)):
                for j, Lj in enumerate(self.L[:i+1]):
                    s = sum(Li[k] * Lj[k] for k in range(j))
                    Li[j] = numpy.sqrt(Ai[i] - s) if (i == j) else \
                              (1.0 / Lj[j] * (Ai[j] - s))
                          
    def _MultiplyByLInv(self,a):
        """
        Returns the product L^({)-1)*a,
        where L is Choletsky decomposed part of the matrix, a is an array.
        Computation by 

        Parameters
        ----------
        a : list of floats            

        Returns
        -------
        L^(-1)*a

        """
        x=[]
        for i in range(self.numberOfPoints):
            x.append((a[i]-numpy.sum([self.L[i][j]*x[j] for j in range(len(x))]))/self.L[i][i])
            
        return x
        
        
    def _CalculateA(self):
        """ Evaluate the matrix A which is needed for the estimation of systematic shifts.
        
        The matrix A defined as 
        A[a,b]=delta[a,b]+sum(corrErr[a,i]corrErr[b,i]/(uncorrErr[i]^2), i in number of points)
        """
        
        if(self.isVdiagonal):
            self.matrixA=numpy.identity(0)
            return
        
        ## determine the number of correlated errors
        sc=self.numOfCorrErr+self.numOfNormErr
        
        ## building list sum_p sigma[p,j](corr)/sqrt{sigma[p](uncorr)}
        cErr=numpy.zeros((self.numberOfPoints,sc))
        for i in range(self.numberOfPoints):
            for j in range(sc):
                cErr[i][j]=self._listOfCorrErrors[i][j]/numpy.sqrt(self._listOfVariances[i])
        
        
        self.matrixA=numpy.identity(sc)
        ### add contribution of a point to matrix A
        for i in range(sc):
            for j in range(sc):
                for k in range(self.numberOfPoints):          
                    self.matrixA[i][j]+=cErr[k][i]*cErr[k][j]
    
    def MatchWithData(self,theoryPrediction):
        """
        Adjust the input theory vector to the experimental data.
        Multiply by the proper translation factor and normalize if needed.

        Parameters
        ----------
        theoryPrediction : list of floats
            The theory prediction to be compared with the data

        Returns
        -------
        list of floats
            Matched prediction

        """
        res=[]
        norm=1
        ## weight theory with proper factor
        for i in range(self.numberOfPoints):
            res.append(theoryPrediction[i]*self.points[i]["thFactor"])
        
        ##if necessary normalize
        if self.isNormalized:
            if self.normalizationMethod=="integral":
                ### normalization by bin-by-bin area
                normTh=sum([(self.points[i]["qT"][1]-self.points[i]["qT"][0])*res[i] for i in range(self.numberOfPoints)])
                norm=self._normExp/normTh       
                return [v*norm for v in res]
            elif self.normalizationMethod=="bestChi2":                
                ### normalization by best chi^2 value
                normTh=self.FindBestNorm(res)                
                return [v*normTh for v in res]
        else:
            return res                  
    
    def chi2(self,theoryPrediction):
        """
        Evaluate chi^2=x.V^{-1}.x, where x is theory prediction

        Parameters
        ----------
        theoryPrediction : list of floats
           The theory prediction (matched to the data)

        Returns
        -------
        Float
            chi^2

        """
        diffX=[(theoryPrediction[i]-self.points[i]["xSec"]) for i in range(self.numberOfPoints)]
        Lx=self._MultiplyByLInv(diffX)
        return numpy.dot(Lx,Lx)
    
    def DetermineSystematicShift(self,theoryPrediction):
        """
        Determine the systematic shift by nuisance parameters evaluations.        

        Parameters
        ----------
        theoryPrediction : list of floats
            List theory predictions (matched) to be compared to the data

        Returns
        -------
        result : list of floats
            List of systematic shift point-per-point

        """
        
        sc=self.numOfCorrErr+self.numOfNormErr
        ## empty vector for rho        
        rho=numpy.zeros(sc)
        
        for i in range(self.numberOfPoints):
            ## the vector rho[i]= (exp-theory)*corrErr[i]/uncorr^2            
            rho+=(self.points[i]["xSec"] -theoryPrediction[i])*self._listOfCorrErrors[i]/self._listOfVariances[i]
        
        ## these are nuisance parameters
        lambd=numpy.matmul(self.matrixAinverse,rho)
        ## the shifts
        dd=[]
        for ss in self._listOfCorrErrors:
            dummy=0.
            for i in range(sc):
                dummy+=lambd[i]*ss[i]
            dd.append(dummy)
        
        return dd
    
    def DetermineAvarageSystematicShift(self,theoryPrediction):
        """
        Determine the mean systematic shift in % by nuisance parameters evaluations.

        Parameters
        ----------
        theoryPrediction : list of floats
             List theory predictions (matched) to be compared to the data

        Returns
        -------
        float
            The avarage excess/deficit of the theory to data (due to systematic shift)

        """
        dd=self.DetermineSystematicShift(theoryPrediction)
        tmp=[]
        for i in range(len(dd)):
            if self.points[i]["xSec"]!=0:
                tmp.append(dd[i]/self.points[i]["xSec"])
            
        return sum(tmp)/len(dd)
    
    def DecomposeChi2(self,theoryPrediction):
        """
        Determine the nuisance parameters and perform the decomposition of chi^2 to correlated and uncorrelated parts.

        Parameters
        ----------
        theoryPrediction : list of floats
             List theory predictions (matched) to be compared to the data
             
        Returns
        -------
        list of 3 floats
            chi2(uncorr), chi2(corr), chi2(total)

        """
        sc=self.numOfCorrErr+self.numOfNormErr
        ## empty vector for rho        
        rho=numpy.zeros(sc)
        
        for i in range(self.numberOfPoints):
            ## the vector rho[i]= (exp-theory)*corrErr[i]/uncorr^2            
            rho+=(self.points[i]["xSec"] -theoryPrediction[i])*self._listOfCorrErrors[i]/self._listOfVariances[i]
        
        ## these are nuisance parameters
        lambd=numpy.matmul(self.matrixAinverse,rho)
        ## the shifts
        dd=[]
        for ss in self._listOfCorrErrors:
            dummy=0.
            for i in range(sc):
                dummy+=lambd[i]*ss[i]
            dd.append(dummy)
        
        ### uncorrelated part of chi**2
        chiD=0.
        for j in range(self.numberOfPoints):
            chiD+=(self.points[j]["xSec"]-theoryPrediction[j]-dd[j])**2/self._listOfVariances[j]
        
        chiL=0.
        for i in range(sc):
            chiL+=lambd[i]**2
        
        return [chiD,chiL,chiD+chiL]
    
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
        Ls=self._MultiplyByLInv([p["xSec"] for p in self.points])
        Lt=self._MultiplyByLInv(theoryPrediction)
        return numpy.dot(Ls,Lt)/numpy.dot(Lt,Lt)
        

        
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
        # it is prod(1+ rand * normErr), because normErr are given relative to xSec (in %)
        resFactor=1.
        for err in self.normErr:
            resFactor=resFactor*(1.+numpy.random.normal()*err)
        
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
        
        ### for weighted processes I add a column
        includeWeightProc=('weightProcess' in list(self.points[0].keys()))
        ### Points table header
        
        file.write("Point id,process id,")
        if includeWeightProc:
            file.write("weightProcess,")
        file.write("s[GeV^2],<Q>[GeV],Qmin[GeV],Qmax[GeV],")
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
            if includeWeightProc:
                file.write(str(p["weightProcess"]).replace(","," -").replace("[","").replace("]",""))            
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
    import os.path
    if not os.path.exists(path):
        raise FileNotFoundError('data-set-file '+path+' NOT FOUND')
        
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
    
    ### Check for extra columns
    line=line.split(",")    
    includeWeightProc=('weightProcess' in line)
    
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
        if includeWeightProc:
            p["weightProcess"]=[int(ll) for ll in line[2].split("-")]
            k=2
        else:
            k=1
        p["s"]=float(line[k+1])
        p["<Q>"]=float(line[k+2])
        p["Q"]=[float(line[k+3]),float(line[k+4])]
        k=k+4
        
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