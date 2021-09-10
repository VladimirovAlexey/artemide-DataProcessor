#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:06:14 2020

Class ArtemideReplicaSet incorporates the information 
about a particular replica set, encode the to parametrize TMDs in artemide,
it reads the replica files and allow to operate with it.

The object ArtemideReplicaSet has following fields:
    
numberOfReplicas:  the number of replicas


@author: vla18041
"""

def ReadRepFile(path):
    import os.path
    if not os.path.exists(path):
        raise FileNotFoundError('replica-file '+path+' NOT FOUND')
    
    
    with open(path,"r") as file:
        listFromF=file.readlines()
    
    ### I read line by line and delete lines fron the list
    ### it corresponds to the logical way I saved it in the fortran
    ### the line start with *?? indicate that the next line is the variable
    
    ### search for version
    while not listFromF[0].startswith("*V   "):
        listFromF.pop(0)
    listFromF.pop(0)
    ver=int(listFromF.pop(0))
    
    ### search for name entry
    while not listFromF[0].startswith("*A   "):
        listFromF.pop(0)
    listFromF.pop(0)
        
    name=(listFromF.pop(0)).replace("\n","")
    rSet=ArtemideReplicaSet()
    rSet.name=name
    
    ### search for indexing parameters
    while not listFromF[0].startswith("*B   "):
        listFromF.pop(0)
    ## total size
    while not listFromF[0].startswith("*0   "):
        listFromF.pop(0)
    listFromF.pop(0)
    
    rSet._totalLength=int(listFromF.pop(0))
    
    ## TMDR size
    while not listFromF[0].startswith("*3   "):
        listFromF.pop(0)
    listFromF.pop(0)
    
    line=(listFromF.pop(0)).split(",")
    
    rSet._TMDRstart=int(line[0])-2
    rSet._TMDRend=int(line[1])-1
    
    ## TMDPDF size
    while not listFromF[0].startswith("*4   "):
        listFromF.pop(0)
    listFromF.pop(0)
    
    line=(listFromF.pop(0)).split(",")
    rSet._uTMDPDFstart=int(line[0])-2
    rSet._uTMDPDFend=int(line[1])-1
    
    ## TMDFF size
    while not listFromF[0].startswith("*5   "):
        listFromF.pop(0)
    listFromF.pop(0)
    
    line=(listFromF.pop(0)).split(",")
    rSet._uTMDFFstart=int(line[0])-2
    rSet._uTMDFFend=int(line[1])-1
    
    if(ver>=3):
        ## lpTMDPDF size
        while not listFromF[0].startswith("*11  "):
            listFromF.pop(0)
        listFromF.pop(0)
        
        line=(listFromF.pop(0)).split(",")
        rSet._lpTMDPDFstart=int(line[0])-2
        rSet._lpTMDPDFend=int(line[1])-1   
    else:
        rSet._lpTMDPDFstart=-2
        rSet._lpTMDPDFend=-1   
    
    ## SiversTMDPDF appeared only in 15
    if(ver>=15):
        ## SiversTMDPDF size
        while not listFromF[0].startswith("*12  "):
            listFromF.pop(0)
        listFromF.pop(0)
        
        line=(listFromF.pop(0)).split(",")
        rSet._SiversTMDPDFstart=int(line[0])-2
        rSet._SiversTMDPDFend=int(line[1])-1
    else:
        rSet._SiversTMDPDFstart=-2
        rSet._SiversTMDPDFend=-1
    
    
    #### Some files can contain *BB section which has number of PDF replicas        
    #### it can happen only for versions >=15
    if(ver>=15):
        while not (listFromF[0].startswith("*BB  ") or listFromF[0].startswith("*C   ")):
            listFromF.pop(0)
    else:
        while (not listFromF[0].startswith("*C   ")):
            listFromF.pop(0)
    
    ### case with indicated numbers of PDFs
    if( listFromF[0].startswith("*BB  ") ):
        listFromF.pop(0)
        ## TMDPDF size
        while not listFromF[0].startswith("*4   "):
            listFromF.pop(0)
        listFromF.pop(0)
        
        line=(listFromF.pop(0)).split(",")
        rSet._c_uTMDPDFstart=int(line[0])-2
        rSet._c_uTMDPDFend=int(line[1])-1
        
        ## TMDFF size
        while not listFromF[0].startswith("*5   "):
            listFromF.pop(0)
        listFromF.pop(0)
        
        line=(listFromF.pop(0)).split(",")
        rSet._c_uTMDFFstart=int(line[0])-2
        rSet._c_uTMDFFend=int(line[1])-1
        
        ## lpTMDPDF size
        while not listFromF[0].startswith("*11  "):
            listFromF.pop(0)
        listFromF.pop(0)
        
        line=(listFromF.pop(0)).split(",")
        rSet._c_lpTMDPDFstart=int(line[0])-2
        rSet._c_lpTMDPDFend=int(line[1])-1   
        
        ## SiversTMDPDF size
        while not listFromF[0].startswith("*12  "):
            listFromF.pop(0)
        listFromF.pop(0)
        
        line=(listFromF.pop(0)).split(",")
        rSet._c_SiversTMDPDFstart=int(line[0])-2
        rSet._c_SiversTMDPDFend=int(line[1])-1
        
        ### search for number of replicas
        while not listFromF[0].startswith("*C   "):
            listFromF.pop(0)
    
    else:
        #### case with *C should be filled
        rSet._c_uTMDPDFstart=-2
        rSet._c_uTMDPDFend=-1
        rSet._c_uTMDFFstart=-2
        rSet._c_uTMDFFend=-1
        rSet._c_lpTMDPDFstart=-2
        rSet._c_lpTMDPDFend=-1 
        rSet._c_SiversTMDPDFstart=-2
        rSet._c_SiversTMDPDFend=-1
    
    #### computing the length of only NP input
    rSet._PDFlength= (rSet._c_uTMDPDFend-rSet._c_uTMDPDFstart if rSet._c_uTMDPDFstart>0 else 0)\
                    + (rSet._c_uTMDFFend-rSet._c_uTMDFFstart if rSet._c_uTMDFFstart>0 else 0)\
                    + (rSet._c_lpTMDPDFend-rSet._c_lpTMDPDFstart if rSet._c_lpTMDPDFstart>0 else 0)\
                    + (rSet._c_SiversTMDPDFend-rSet._c_SiversTMDPDFstart if rSet._c_SiversTMDPDFstart>0 else 0)
    rSet._NPLength=rSet._totalLength-rSet._PDFlength    
    
    ##########################
    ### Here we at position *C
    listFromF.pop(0)
    
    rSet.numberOfReplicas=int(listFromF.pop(0))
    
    ### search for technical replicas
    while not listFromF[0].startswith("*D   "):
        listFromF.pop(0)
    listFromF.pop(0)
    
    line=listFromF.pop(0).split(",")
    rSet.initialReplica=[float(x) for x in line[1:rSet._NPLength]]+[int(x) for x in line[rSet._NPLength:]]
        
    line=listFromF.pop(0).split(",")
    rSet.meanReplica=[float(x) for x in line[1:rSet._NPLength]]+[int(x) for x in line[rSet._NPLength:]]
    
    ### search for the start of the replica lst
    while not listFromF[0].startswith("*R   "):
        listFromF.pop(0)
    listFromF.pop(0)
    
    rSet.replicaList=[]
    for line in listFromF[:rSet.numberOfReplicas]:
        lineX=line.split(",")
        rSet.replicaList.append([float(x) for x in lineX[1:rSet._NPLength]]+[int(x) for x in lineX[rSet._NPLength:]])
    
    return rSet
        

class ArtemideReplicaSet:
    def ArtemideReplicaSet(self):
        
        self.name="NONAME"
        self.numberOfReplicas=0
                
        self._totalLength=0
        ### length of the NP part of the replica
        self._NPLength=0
        ### length of the PDF part of the replica
        self._PDFlength=0
        
        ### inidices of NP parameters
        self._TMDRstart=0
        self._uTMDPDFstart=0
        self._uTMDFFstart=0
        self._lpTMDPDFstart=0
        self._SiversTMDPDFstart=0
        
        self._TMDRend=0
        self._uTMDPDFend=0
        self._uTMDFFend=0
        self._lpTMDPDFend=0
        self._SiversTMDPDFend=0
        
        ### inidices of collinear input
        self._c_uTMDPDFstart=0
        self._c_uTMDFFstart=0
        self._c_lpTMDPDFstart=0
        self._c_SiversTMDPDFstart=0
        
        self._c_TMDRend=0
        self._c_uTMDPDFend=0
        self._c_uTMDFFend=0
        self._c_lpTMDPDFend=0
        self._c_SiversTMDPDFend=0
        
        ### flag which indicates that the replica incorporate PDF variation
        self._includesPDFvariation=False
        
        ### replica suggested for the initialization
        self.initialReplica=[]
        ### mean replica
        self.meanReplica=[]
        ### all other replicas
        self.replicaList=[]
        ########## Four integer numbers for uTMDPDF, uTMDFF, lpTMDPDF, Sivers
        ### PDF replica suggested for the initialization
        self.initialPDFReplica=[0,0,0,0]
        ### PDF mean replica
        self.meanPDFReplica=[0,0,0,0]
        ### all other PDF replicas
        self.PDFreplicaList=[0,0,0,0]
        
    def __repr__ (self):
        return "<ArtemideReplicaSet: %s with %s replicas>" % (self.name, self.numberOfReplicas)
    
    
    def SetReplica(self,num=0,part='full'):
        """
        Set the replica according to the list

        Parameters
        ----------
        num : -1 = initla replica
              0  = mean replica
              1.... = replica from list
            The default is 0.
        
        part : string, optional
            Specification which part of the replica to be set. The default is "full".
            Possible values: 'full', 'TMDR', 'uTMDPDF', 'uTMDFF', 'lpTMDPDF', 'SiversTMDPDF', etc

        Returns
        -------
        None.

        """        
        import harpy
        
        if(num==0):
            r=self.meanReplica
        elif(num==-1):
            r=self.initialReplica
        else:
            r=self.replicaList[num-1]
        
        ## resolvinf the part-condition
        if(part=="full"):
            doTMDR=True
            douTMDPDF=True
            douTMDFF=True
            dolpTMDPDF=True
            doSiversTMDPDF=True
        else:
            doTMDR=False
            douTMDPDF=False
            douTMDFF=False
            dolpTMDPDF=False
            doSiversTMDPDF=False
            
            if(part=='TMDR'):
                doTMDR=True
            if(part=='uTMDPDF'):
                douTMDPDF=True
            if(part=='uTMDFF'):
                douTMDFF=True
            if(part=='lpTMDPDF'):
                dolpTMDPDF=True
            if(part=='SiversTMDPDF'):
                doSiversTMDPDF=True
            
        
            
        
        if(doTMDR and self._TMDRend>=self._TMDRstart+1 >0):
            harpy.setNPparameters_TMDR(r[self._TMDRstart:self._TMDRend])
            
        if(douTMDPDF and self._uTMDPDFend>=self._uTMDPDFstart+1 >0):            
            if(self._c_uTMDPDFend>=self._c_uTMDPDFstart+1 >0):
                pass
                harpy.setPDFreplica(r[self._c_uTMDPDFstart])
            harpy.setNPparameters_uTMDPDF(r[self._uTMDPDFstart:self._uTMDPDFend])            
            
        if(douTMDFF and self._uTMDFFend>=self._uTMDFFstart+1 >0):
            if(self._c_uTMDFFend>=self._c_uTMDFFstart+1 >0):
                print("Modification of FF replica is not implimented")
            harpy.setNPparameters_uTMDFF(r[self._uTMDFFstart:self._uTMDFFend])
            
        if(dolpTMDPDF and self._lpTMDPDFend>=self._lpTMDPDFstart+1 >0):
            if(self._c_lpTMDPDFend>=self._c_lpTMDPDFstart+1 >0):
                print("Modification of lpPDF replica is not implimented")
            harpy.setNPparameters_lpTMDPDF(r[self._lpTMDPDFstart:self._lpTMDPDFend])
            
        if(doSiversTMDPDF and self._SiversTMDPDFend>=self._SiversTMDPDFstart+1 >0):
            if(self._SiversTMDPDFend>=self._SiversTMDPDFstart+1 >0):
                print("Modification of Sivers replica is not implimented")
            harpy.setNPparameters_SiversTMDPDF(r[self._SiversTMDPDFstart:self._SiversTMDPDFend])
            
    def GetReplica(self,num,part="full"):
        """
        Returns the values of parameters for the numbers replica

        Parameters
        ----------
        num : int
            Number of the replica (-1 = initial, 0 = mean, ...)
        part : string, optional
            Specification which part of the replica to return. The default is "full".
            Possible values: 'full', 'TMDR', 'uTMDPDF', 'uTMDFF', 'lpTMDPDF', 'SiversTMDPDF',
            'uTMDPDF-PDF', 'uTMDFF-PDF', 'lpTMDPDF-PDF', 'SiversTMDPDF-PDF'

        Returns
        -------
        array of floats

        """
        if(num==0):
            r=self.meanReplica
        elif(num==-1):
            r=self.initialReplica
        else:
            r=self.replicaList[num-1]
        
        if(part=="full"):
            return r
        elif(part=="TMDR"):
            if(self._TMDRend>=self._TMDRstart+1 >0):
                return r[self._TMDRstart:self._TMDRend]
            else:
                return []
        elif(part=="uTMDPDF"):
            if(self._uTMDPDFend>=self._uTMDPDFstart+1 >0):
                return r[self._uTMDPDFstart:self._uTMDPDFend]
            else:
                return []
        elif(part=="uTMDFF"):
            if(self._uTMDFFend>=self._uTMDFFstart+1 >0):
                return r[self._uTMDFFstart:self._uTMDFFend]
            else:
                return []
        elif(part=="lpTMDPDF"):
            if(self._lpTMDPDFend>=self._lpTMDPDFstart+1 >0):
                return r[self._lpTMDPDFstart:self._lpTMDPDFend]
            else:
                return[]
        elif(part=="SiversTMDPDF"):
            if(self._SiversTMDPDFend>=self._SiversTMDPDFstart+1 >0):
                return r[self._SiversTMDPDFstart:self._SiversTMDPDFend]
            else:
                return []
        elif(part=="uTMDPDF-PDF"):
            if(self._c_uTMDPDFend>=self._c_uTMDPDFstart+1 >0):
                return r[self._c_uTMDPDFstart:self._c_uTMDPDFend]
            else:
                return []
        elif(part=="uTMDFF-PDF"):
            if(self._c_uTMDFFend>=self._c_uTMDFFstart+1 >0):
                return r[self._c_uTMDFFstart:self._c_uTMDFFend]
            else:
                return []
        elif(part=="lpTMDPDF-PDF"):
            if(self._c_lpTMDPDFend>=self._c_lpTMDPDFstart+1 >0):
                return r[self._c_lpTMDPDFstart:self._c_lpTMDPDFend]
            else:
                return[]
        elif(part=="SiversTMDPDF-PDF"):
            if(self._c_SiversTMDPDFend>=self._c_SiversTMDPDFstart+1 >0):
                return r[self._c_SiversTMDPDFstart:self._c_SiversTMDPDFend]
            else:
                return []
        else:
            raise ValueError("part should correspond to a TMD")
        