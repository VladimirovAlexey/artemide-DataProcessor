#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 06 Sep 2021

Collections of routines to write TMD grids in the TMDlib format (or similar)

@author: vla18041
"""

#############################################################
# Default values for grids grids
#############################################################

### Q[GeV]
Qrange_default= [1., 1.11803, 1.22474, 1.4, 1.58114, 1.78885, 2., 2.23607, 2.52982, 2.82843,
         3.16228, 3.4641, 4.75, 5.09902, 6.32456, 7.1, 8., 10., 11.1803, 12.2475,
         14., 15.8114, 17.8885, 20., 22.3607, 25.2982, 28.2843, 31.6228, 34.641, 47.5,
         50.9902, 63.2456, 71, 80, 100, 111.803, 122.475, 140, 158.114, 178.885, 
         200.]
### X for TMDPDFs
XrangePDF_default= [0.00001, 0.00002, 0.00004, 0.00006, 0.00008, 0.0001, 0.0002, 0.0004, 0.0006, 0.0008,
          0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.0045, 0.005, 0.0055,
          0.006, 0.0065, 0.007, 0.0075, 0.008, 0.0085, 0.009, 0.00925, 0.0095, 0.00975,
          0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055,
          0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.0925, 0.095, 0.0975,
          0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
          0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.925, 0.95, 0.975,
          1]
### X for TMDFFs
XrangeFF_default = [0.05,0.055,0.06,0.065,0.07,0.08,0.09,0.1,
          0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
          0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.]
### R=qT/Q for grids in the momentum space
Rrange_default= [0.0001, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.02, 0.03, 0.04, 0.05,
         0.06, 0.07, 0.08, 0.09, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225,
         0.25, 0.275, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
         0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
         1.7, 1.8, 1.9, 2.001]

Brange_default= [0., 0.01, 0.025, 0.05, 0.1, 
                 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
                 0.6, 0.65, 0.7, 0.75, 0.8, 0.86, 0.93, 1., 
                 1.1, 1.2, 1.3, 1.4, 1.5, 1.65, 1.8, 1.95, 2.1,
                 2.3, 2.5, 2.75, 3., 3.25, 3.5, 4., 4.5, 5., 6, 8, 10.
                 ]

KTrange_default= [0.00001, 0.01, 0.025, 0.05, 0.1, 
                 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
                 0.6, 0.65, 0.7, 0.75, 0.8, 0.86, 0.93, 1., 
                 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1,
                 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9,
                 4.1, 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 
                 10.,11.,12.,13.,14.,15.,20.,25.
                 ]

#%%
###########################################
## Save optimal TMDPDF in b-space
###########################################
def SaveGrid_optimal(path,Xrange=XrangePDF_default,Brange=Brange_default,PDF="uTMDPDF",h=1,Q=-1):
    """
    Saves the grid for optimal TMDPDF in b space, as returned from the artmide (should be setup)

    Parameters
    ----------
    path : TYPE
        DESCRIPTION.
    Xrange : TYPE, optional
        DESCRIPTION. The default is XrangePDF_default.
    Brange : TYPE, optional
        DESCRIPTION. The default is Brange_default.
    PDF : TYPE, optional
        DESCRIPTION. Which TMD to save. Default=uTMDPDF
    h : TYPE, optional
        DESCRIPTION. Number of hadron
    Q : TYPE, optional
        DESCRIPTION. The grid is saved at particular Q. If Q is negative the optimal grid is saved.

    Returns
    -------
    None.

    """

    from numpy import zeros    
    import harpy
    import time
    
    startTime=time.time()     
    
    valuesList = {
    -5: zeros([len(Xrange),len(Brange)]).tolist(),
    -4: zeros([len(Xrange),len(Brange)]).tolist(),
    -3: zeros([len(Xrange),len(Brange)]).tolist(),
    -2: zeros([len(Xrange),len(Brange)]).tolist(),
    -1: zeros([len(Xrange),len(Brange)]).tolist(),
    1: zeros([len(Xrange),len(Brange)]).tolist(),
    2: zeros([len(Xrange),len(Brange)]).tolist(),
    3: zeros([len(Xrange),len(Brange)]).tolist(),
    4: zeros([len(Xrange),len(Brange)]).tolist(),
    5: zeros([len(Xrange),len(Brange)]).tolist()
    }
        
    for i in range(len(Xrange)):
        for j in range(len(Brange)):                
                xval=float(Xrange[i])
                rval=float(Brange[j])
                
                if(xval>0.999):
                    TMDval=[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
                else:                    
                    if(PDF=="uTMDFF"):
                        TMDval=harpy.get_uTMDFF(xval,rval,h,includeGluon=False,mu=Q)
                    else:
                        TMDval=harpy.get_uTMDPDF(xval,rval,h,includeGluon=False,mu=Q)
                
                valuesList[-5][i][j]='{:g}'.format(xval*TMDval[0])
                valuesList[-4][i][j]='{:g}'.format(xval*TMDval[1])
                valuesList[-3][i][j]='{:g}'.format(xval*TMDval[2])
                valuesList[-2][i][j]='{:g}'.format(xval*TMDval[3])
                valuesList[-1][i][j]='{:g}'.format(xval*TMDval[4])
                #valuesList[0][i][j][k]='{:g}'.format(xval*TMDval[5])
                valuesList[1][i][j]='{:g}'.format(xval*TMDval[6])
                valuesList[2][i][j]='{:g}'.format(xval*TMDval[7])
                valuesList[3][i][j]='{:g}'.format(xval*TMDval[8])
                valuesList[4][i][j]='{:g}'.format(xval*TMDval[9])
                valuesList[5][i][j]='{:g}'.format(xval*TMDval[10])   
                
    startTime=time.time()
    with open(path, 'w') as outfile:        
        outfile.write("xg: "+str(Xrange)+"\n")
        outfile.write("bg: "+str(Brange)+"\n")
        outfile.write("TMDs: "+str(valuesList).replace("'","")+"\n")
    endTime=time.time()
    print('Grid written at  : ',path)
    print('Computation time  : ',endTime-startTime,' sec.')

#%%
###########################################
## Save optimal TMDPDF in kT-space
###########################################
def SaveGrid_optimal_kT(path,Xrange=XrangePDF_default,KTrange=KTrange_default,PDF="uTMDPDF",h=1,Q=-1.):
    """
    Saves the grid for optimal TMDPDF in b space, as returned from the artmide (should be setup)

    Parameters
    ----------
    path : TYPE
        DESCRIPTION.
    Xrange : TYPE, optional
        DESCRIPTION. The default is XrangePDF_default.
    Brange : TYPE, optional
        DESCRIPTION. The default is Brange_default.
    PDF : TYPE, optional
        DESCRIPTION. Which TMD to save. Default=uTMDPDF.
    h : TYPE, optional
        DESCRIPTION. Number of hadron.
    Q : TYPE, optional
        DESCRIPTION. The grid is saved at particular Q. If Q is negative the optimal grid is saved.

    Returns
    -------
    None.

    """

    from numpy import zeros    
    import harpy
    import time
    
    startTime=time.time()     
    
    valuesList = {
    -5: zeros([len(Xrange),len(KTrange)]).tolist(),
    -4: zeros([len(Xrange),len(KTrange)]).tolist(),
    -3: zeros([len(Xrange),len(KTrange)]).tolist(),
    -2: zeros([len(Xrange),len(KTrange)]).tolist(),
    -1: zeros([len(Xrange),len(KTrange)]).tolist(),
    1: zeros([len(Xrange),len(KTrange)]).tolist(),
    2: zeros([len(Xrange),len(KTrange)]).tolist(),
    3: zeros([len(Xrange),len(KTrange)]).tolist(),
    4: zeros([len(Xrange),len(KTrange)]).tolist(),
    5: zeros([len(Xrange),len(KTrange)]).tolist()
    }
        
    for i in range(len(Xrange)):
        for j in range(len(KTrange)):                
                xval=float(Xrange[i])
                rval=float(KTrange[j])
                
                if(xval>0.999):
                    TMDval=[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
                else:                    
                    if(PDF=="uTMDFF"):
                        TMDval=harpy.get_uTMDFF_kT(xval,rval,h,includeGluon=False,mu=Q)
                    else:
                        TMDval=harpy.get_uTMDPDF_kT(xval,rval,h,includeGluon=False,mu=Q)
                
                valuesList[-5][i][j]='{:g}'.format(xval*TMDval[0])
                valuesList[-4][i][j]='{:g}'.format(xval*TMDval[1])
                valuesList[-3][i][j]='{:g}'.format(xval*TMDval[2])
                valuesList[-2][i][j]='{:g}'.format(xval*TMDval[3])
                valuesList[-1][i][j]='{:g}'.format(xval*TMDval[4])
                #valuesList[0][i][j][k]='{:g}'.format(xval*TMDval[5])
                valuesList[1][i][j]='{:g}'.format(xval*TMDval[6])
                valuesList[2][i][j]='{:g}'.format(xval*TMDval[7])
                valuesList[3][i][j]='{:g}'.format(xval*TMDval[8])
                valuesList[4][i][j]='{:g}'.format(xval*TMDval[9])
                valuesList[5][i][j]='{:g}'.format(xval*TMDval[10])   
                
    startTime=time.time()
    with open(path, 'w') as outfile:        
        outfile.write("xg: "+str(Xrange)+"\n")
        outfile.write("bg: "+str(KTrange)+"\n")
        outfile.write("TMDs: "+str(valuesList).replace("'","")+"\n")
    endTime=time.time()
    print('Grid written at  : ',path)
    print('Computation time  : ',endTime-startTime,' sec.')

#%%
#######################################
# Save TMDPDF in b-space at Q, Q^2
#######################################
def SaveGrid_Q(path,Qrange=Qrange_default, Xrange=XrangePDF_default,Brange=Brange_default,PDF="uTMDPDF",h=1,includeGluon=False):
    """
    Saves the grid for TMDPDF in b space at Q, Q^2, as returned from the artmide (should be setup)

    Parameters
    ----------
    path : TYPE
        DESCRIPTION.
    Qrange : TYPE, optional
        DESCRIPTION. The default is Qrange_default.
    Xrange : TYPE, optional
        DESCRIPTION. The default is XrangePDF_default.
    Brange : TYPE, optional
        DESCRIPTION. The default is Brange_default.
    PDF : TYPE, optional
        DESCRIPTION. Which TMD to save. Default=uTMDPDF
    h : TYPE, optional
        DESCRIPTION. Number of hadron

    Returns
    -------
    None.

    """
    from numpy import zeros
    import harpy
    import time
    
    startTime=time.time()
        
    valuesList = {
    -5: zeros([len(Qrange),len(Xrange),len(Brange)]).tolist(),
    -4: zeros([len(Qrange),len(Xrange),len(Brange)]).tolist(),
    -3: zeros([len(Qrange),len(Xrange),len(Brange)]).tolist(),
    -2: zeros([len(Qrange),len(Xrange),len(Brange)]).tolist(),
    -1: zeros([len(Qrange),len(Xrange),len(Brange)]).tolist(),
    0: zeros([len(Qrange),len(Xrange),len(Brange)]).tolist(),
    1: zeros([len(Qrange),len(Xrange),len(Brange)]).tolist(),
    2: zeros([len(Qrange),len(Xrange),len(Brange)]).tolist(),
    3: zeros([len(Qrange),len(Xrange),len(Brange)]).tolist(),
    4: zeros([len(Qrange),len(Xrange),len(Brange)]).tolist(),
    5: zeros([len(Qrange),len(Xrange),len(Brange)]).tolist()
    }
    for i in range(len(Qrange)):
        for j in range(len(Xrange)):
            for k in range(len(Brange)):
                Qval=float(Qrange[i])
                xval=float(Xrange[j])
                rval=float(Brange[k])
                
                if(xval==0.999):
                    TMDval=[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
                else:
                    if(PDF=="uTMDFF"):
                        TMDval=harpy.get_uTMDFF(xval,rval,h,Qval,Qval**2,includeGluon=includeGluon)
                    else:
                        TMDval=harpy.get_uTMDPDF(xval,rval,h,Qval,Qval**2,includeGluon=includeGluon)
                
                valuesList[-5][i][j][k]='{:g}'.format(xval*TMDval[0])
                valuesList[-4][i][j][k]='{:g}'.format(xval*TMDval[1])
                valuesList[-3][i][j][k]='{:g}'.format(xval*TMDval[2])
                valuesList[-2][i][j][k]='{:g}'.format(xval*TMDval[3])
                valuesList[-1][i][j][k]='{:g}'.format(xval*TMDval[4])
                if(includeGluon): valuesList[0][i][j][k]='{:g}'.format(xval*TMDval[5])
                valuesList[1][i][j][k]='{:g}'.format(xval*TMDval[6])
                valuesList[2][i][j][k]='{:g}'.format(xval*TMDval[7])
                valuesList[3][i][j][k]='{:g}'.format(xval*TMDval[8])
                valuesList[4][i][j][k]='{:g}'.format(xval*TMDval[9])
                valuesList[5][i][j][k]='{:g}'.format(xval*TMDval[10])                

    with open(path, 'w') as outfile:
        outfile.write("Qg: "+str(Qrange)+"\n")
        outfile.write("xg: "+str(Xrange)+"\n")
        outfile.write("bg: "+str(Brange)+"\n")
        outfile.write("TMDs: "+str(valuesList).replace("'","")+"\n")
    endTime=time.time()
    print('Grid written at  : ',path)
    print('Computation time  : ',endTime-startTime,' sec.')

#%%
#######################################
# Save TMDPDF in kT-space at Q, Q^2
#######################################
def SaveGrid_kT(path,Qrange=Qrange_default, Xrange=XrangePDF_default,Rrange=Rrange_default,PDF="uTMDPDF",h=1):
    """
    Saves the grid for optimal TMDPDF in kT space at Q, Q^2, as returned from the artmide (should be setup)

    Parameters
    ----------
    path : TYPE
        DESCRIPTION.
    Qrange : TYPE, optional
        DESCRIPTION. The default is Qrange_default.
    Xrange : TYPE, optional
        DESCRIPTION. The default is XrangePDF_default.
    Rrange : TYPE, optional
        DESCRIPTION. The default is Rrange_default.
    PDF : TYPE, optional
        DESCRIPTION. Which TMD to save. Default=uTMDPDF.
    h : TYPE, optional
        DESCRIPTION. Number of hadron

    Returns
    -------
    None.

    """
    from numpy import zeros
    import harpy
    import time
    
    startTime=time.time()
    
    valuesList = {
    -5: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist(),
    -4: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist(),
    -3: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist(),
    -2: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist(),
    -1: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist(),
    1: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist(),
    2: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist(),
    3: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist(),
    4: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist(),
    5: zeros([len(Qrange),len(Xrange),len(Rrange)]).tolist()
    }
        
    for i in range(len(Qrange)):
        for j in range(len(Xrange)):
            for k in range(len(Rrange)):
                Qval=float(Qrange[i])
                xval=float(Xrange[j])
                rval=float(Rrange[k])
                
                if(xval>0.999):
                    TMDval=[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
                else:
                    if(PDF=="uTMDFF"):
                        TMDval=harpy.get_uTMDFF_kT(xval,rval*Qval,h,Qval,Qval**2,includeGluon=False)
                    else:
                        TMDval=harpy.get_uTMDPDF_kT(xval,rval*Qval,h,Qval,Qval**2,includeGluon=False)
                    
                
                valuesList[-5][i][j][k]='{:g}'.format(xval*TMDval[0])
                valuesList[-4][i][j][k]='{:g}'.format(xval*TMDval[1])
                valuesList[-3][i][j][k]='{:g}'.format(xval*TMDval[2])
                valuesList[-2][i][j][k]='{:g}'.format(xval*TMDval[3])
                valuesList[-1][i][j][k]='{:g}'.format(xval*TMDval[4])
                #valuesList[0][i][j][k]='{:g}'.format(xval*TMDval[5])
                valuesList[1][i][j][k]='{:g}'.format(xval*TMDval[6])
                valuesList[2][i][j][k]='{:g}'.format(xval*TMDval[7])
                valuesList[3][i][j][k]='{:g}'.format(xval*TMDval[8])
                valuesList[4][i][j][k]='{:g}'.format(xval*TMDval[9])
                valuesList[5][i][j][k]='{:g}'.format(xval*TMDval[10])
          
    with open(path, 'w') as outfile:
        outfile.write("Qg: "+str(Qrange)+"\n")
        outfile.write("xg: "+str(Xrange)+"\n")
        outfile.write("qToQg: "+str(Rrange)+"\n")
        outfile.write("TMDs: "+str(valuesList).replace("'","")+"\n")
    endTime=time.time()
    print('Grid written at  : ',path)
    print('Computation time  : ',endTime-startTime,' sec.')
