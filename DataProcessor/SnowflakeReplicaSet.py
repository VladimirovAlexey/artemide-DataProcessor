#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:06:14 2020

Class SnowflakeReplicaSet incorporates the information 
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
    
    #print("Reading with version =", ver)
    
    ### search for name entry
    while not listFromF[0].startswith("*A   "):
        listFromF.pop(0)
    listFromF.pop(0)
        
    name=(listFromF.pop(0)).replace("\n","")
    rSet=SnowflakeReplicaSet()
    rSet.name=name
    
    ### search for indexing parameters
    while not listFromF[0].startswith("*B   "):
        listFromF.pop(0)
    ## total size
    while not listFromF[0].startswith("*0   "):
        listFromF.pop(0)
    listFromF.pop(0)
    
    rSet._totalLength=int(listFromF.pop(0))
    
    ## Main size
    while not listFromF[0].startswith("*1   "):
        listFromF.pop(0)
    listFromF.pop(0)
    
    line=(listFromF.pop(0)).split(",")
    
    rSet._MAINstart=int(line[0])-2
    rSet._MAINend=int(line[1])-1
    
    while (not listFromF[0].startswith("*C   ")):
        listFromF.pop(0)    
    ##########################
    ### Here we at position *C
    listFromF.pop(0)
    
    rSet.numberOfReplicas=int(listFromF.pop(0))
    
    ### search for technical replicas
    while not listFromF[0].startswith("*D   "):
        listFromF.pop(0)
    listFromF.pop(0)
    
    line=listFromF.pop(0).split(",")
    rSet.meanReplica=[float(x) for x in line[1:]]
    
    ### search for the start of the replica lst
    while not listFromF[0].startswith("*R   "):
        listFromF.pop(0)
    listFromF.pop(0)
    
    
    rSet.replicaList=[]
    for line in listFromF[:rSet.numberOfReplicas]:
        lineX=line.split(",")
        rSet.replicaList.append([float(x) for x in lineX[1:]])
    
    return rSet
        

class SnowflakeReplicaSet:
    def ASnowflakeReplicaSet(self):
        
        self.name="NONAME"
        self.numberOfReplicas=0
                
        self._totalLength=0
        
        ### inidices of NP parameters
        self._MAINstart=0
        
        self._MAINend=0
                
        ### mean replica
        self.meanReplica=[]
        ### all other replicas
        self.replicaList=[]
        
    def __repr__ (self):
        return "<SnowflakeReplicaSet: %s with %s replicas>" % (self.name, self.numberOfReplicas)
    
    
    def SetReplica(self,num=0,part='full'):
        """
        Set the replica according to the list

        Parameters
        ----------
        num : 0  = mean replica
              1.... = replica from list
            The default is 0.
        
        part : string, optional
            Specification which part of the replica to be set. The default is "full".

        Returns
        -------
        None.

        """        
        import harpy
        
        if(num==0):
            r=self.meanReplica
        elif(num==-1):
            r=self.meanReplica
        else:
            r=self.replicaList[num-1]
        
        harpy.setNPparameters_tw3(r[self._MAINstart:self._MAINend])
        
            
    def GetReplica(self,num,part="full"):
        """
        Returns the values of parameters for the numbers replica

        Parameters
        ----------
        num : int
            Number of the replica (-1 = initial, 0 = mean, ...)
        part : string, optional
            Specification which part of the replica to return. The default is "full".
        Returns
        -------
        array of floats

        """
        if(num==0):
            r=self.meanReplica
        elif(num==-1):
            r=self.meanReplica
        else:
            r=self.replicaList[num-1]
        