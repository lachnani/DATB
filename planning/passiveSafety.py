# -*- coding: utf-8 -*-
"""
Created on Sun Apr 20 19:07:30 2025

@author: Hakim Lachnani
"""

import numpy as np
from guidance import linearizedModels as lm

def korim(L, oth, kovLim):
    """
    Assesses passive safety using the Keep Out Region Intersection Method (KORIM) 
    of CLROE intersections.

    Parameters
    ----------
    L : 6x1 double
        CLROE state
    oth : double
        Operational time horizon (look-ahead), given in chief orbit periods
    KovLim : nx1 double
        Keep out volume limits.

    Returns
    -------
    psa : boolean
        Passive safety assessment

    """
    
    A0 = L[0]
    alpha = L[1]
    xOff = L[2]
    yOff = L[3]
    B0 = L[4]
    beta = L[5]
    
    KR = kovLim[0]
    KI = kovLim[1]
    KC = kovLim[2]
    
    
    ##########################
    # Initial Checks
    ##########################
    
    # Check xOff to A0 ratio
    xOffMag = np.abs(xOff)
    if xOffMag >= A0 + KR:
        # Radial outside of KR for all time
        return True
    
    # Check secular drift
    secularStart = yOff
    secularEnd = yOff - 3*oth*np.pi*xOff
    if np.abs(secularStart) >= KI + 2*A0 and np.abs(secularEnd) >= KI + 2*A0 and np.sign(secularStart) == np.sign(secularEnd):
        return True
    
    ##########################
    # Radial/Cross-Track
    ##########################
    
    # Check radial
    ntRFull, num_ntR, ntR = korInterval(A0, alpha, xOff, KR)
        
    # Check cross-track
    ntCFull, num_ntC, ntC = korInterval(B0, beta, 0, KC)
        
    # Determine ntRC
    ntRCFull, num_ntRC, ntRC = rcIntersection(ntRFull, num_ntR, ntR, ntCFull, num_ntC, ntC, beta)             
    if num_ntRC == 0:
        return True
    
    ##########################
    # In-Track
    ########################## 
    
    # If ntRC is full, we will necessarily breach the KOV if in-track crosses 
    # zero
    if ntRCFull and np.sign(secularStart) != np.sign(secularEnd):
        return False
    
    if np.any(ntRC < 0):
        # We need to pay attention to the start and end intervals carefully...
        # We need to additionally check the interval start
        IofZero = -2*A0*np.sin(alpha) + yOff
        if np.abs(IofZero) < KI:
           return False
    
    # Check whether in-track is monotonic
    if A0 == 0:
        iMonotonic = True
    elif xOffMag/A0 >= 4/3:
        iMonotonic = True
    else:
        iMonotonic = False
        
    # Compute extrema if not monotonic
    ntRC = ntRC[0:num_ntRC,:]
    IofNtRC = -2*A0*np.sin(ntRC + alpha) + yOff - 1.5*ntRC*xOff 
    ntExtrema = np.zeros((2))
    extremaBounds = np.zeros((2,2))
    extremaContained = np.array([False,False])
    if not iMonotonic:
        ntExtrema[0] = np.arccos(-0.75*xOff/A0)
        ntExtrema[1] = 2*np.pi - ntExtrema[0]
        ntExtrema = wrapIntTo2Pi(ntExtrema - alpha)
        for i in range(2):
            for j in range(num_ntRC):
                extremaContained[i] = np.bitwise_and(ntExtrema[i] > ntRC[j,0], ntExtrema[i] < ntRC[j,1])
                if extremaContained[i]:
                    extremaBounds[i,0] = ntRC[j,0]
                    extremaBounds[i,1] = ntRC[j,1]
        IofNtExtrema = -2*A0*np.sin(ntExtrema + alpha) + yOff - 1.5*ntExtrema*xOff 
        IofNtExtremaBounds = -2*A0*np.sin(extremaBounds + alpha) + yOff - 1.5*extremaBounds*xOff 
    
    for i in range(int(oth)):
           
        # Check the normal cases
        if not iMonotonic and np.any(extremaContained):
            if np.any(np.bitwise_and(np.abs(IofNtExtrema) < KI , extremaContained)):
               return False
            for j in range(2):
                if extremaContained[j]:
                    if np.sign(IofNtExtremaBounds[j,0]) != np.sign(IofNtExtrema[j]) or \
                        np.sign(IofNtExtremaBounds[j,1]) != np.sign(IofNtExtrema[j]):
                       return False
            IofNtExtrema += - 3*np.pi*xOff
            IofNtExtremaBounds += - 3*np.pi*xOff
        if containedOrSignFlip(IofNtRC[:,0],IofNtRC[:,1],KI):
           return False
        IofNtRC += - 3*np.pi*xOff
    return True

def korInterval(A,phi,b,K):
    """
    Returns intervals in which the oscillator A*cos(x+phi) is between -K and +K.
    Since there are two possible intervals, xK is a 2x2 array.
    Also returns the number of intervals and the fullness flag

    """
    
    xK = np.zeros((2,2))
    
    bMag = np.abs(b)
    if bMag + A <= K:
        # Oscillator inside the KOR for all time
        xKFull = True
        num_xK = 1
        xK[0,0] = 0 
        xK[0,1] = 2*np.pi 
    elif np.abs(bMag - A) < K:
        # Oscillator inside the KOR for one interval
        xKFull = False
        num_xK = 1
        if b > 0:
            # Starting outside the KOR
            xK[0,0] = np.arccos((K-b)/A)
            xK[0,1] = 2*np.pi - xK[0,0]  
        else:
            # Starting inside the KOR
            xK[0,0] = -np.arccos((-K-b)/A)
            xK[0,1] = -xK[0,0]
        xK[0,:] = wrapIntTo2Pi(xK[0,:] - phi)
    else:
        # Oscillator inside the KOR for two intervals
        xKFull = False
        num_xK = 2 
        xK[0,0] = np.arccos((K-b)/A)
        xK[0,1] = np.arccos((-K-b)/A)
        # Since we have two intervals, the whole interval must be crossed before pi
        xK[1,:] = mirrorToPi(xK[0,:]) # Will always be in [pi,2pi]
        xK = xK - phi
        if xK[0,0] < 0 and xK[0,1] < 0:
            # With the phase shift, the first interval has gone negative. Swap
            # and shift by 2*pi
            xK[[0,1]] = xK[[1,0]]
            xK[1,:] += 2*np.pi   
            
    return xKFull, num_xK, xK
    

def containedOrSignFlip(l,r,K):
    """
    Determines whether the left and right hand sides are contained in an interval +-K

    """    
    if np.min(np.abs(l)) < K:
       return True
    elif np.min(np.abs(r)) < K:
       return True
    elif not np.all(np.sign(l) == np.sign(r)):
       return True
    return False
    
    
def wrapIntTo2Pi(interval):
    interval = np.sort(interval)
    if interval[0] < 0 and interval[1] < 0:
        interval += 2*np.pi
    return interval

def mirrorToPi(interval):
    """
    Mirrors angular interval relative to pi

    """
    
    dPi = np.pi - interval
    return np.sort(np.pi + dPi)

def rcIntersection(ntRFull, num_ntR, ntR, ntCFull, num_ntC, ntC, beta):
    """ 
    Compute the intersection between radial and cross-track oscillatory intervals.
    Leverages the fact that the cross-track oscillations will always be symmetric
    about -beta. Therefore the entire problem is shifted, and a duplicate of the lower
    ntR interval is made if it crosses the zero. We therefore end up with up to three
    radial intervals to compare to up to two cross-track intervals.
    
    Assumptions:
        1. ntR has been ordered such no interval is entirely negative or
           entirely greater than 2pi
        2. beta is in [0, 2pi]
    """
    
    ntRC = np.zeros((4,2))
    ntRCFull = ntRFull and ntCFull
    if ntRCFull:
        num_ntRC = 1
        ntRC[0,0] = 0
        ntRC[0,1] = 2*np.pi
    elif ntRFull:
        ntRC[:2,:] = ntC
        num_ntRC = num_ntC
    elif ntCFull:
        ntRC[:2,:] = ntR
        num_ntRC = num_ntR
    else:    
        # Shift both ntR and ntC by beta
        ntR += beta
        ntC += beta
        
        # Compute the interval intersections for all four permutations
        ind = 0
        if num_ntR == 1:
            # Single radial interval; compare to C1 and C2
            for i in range(2):
                if intervalIntersectionCheck(ntR[0,:], ntC[i,:]):
                    ntRC[i,:] = intervalIntersectionCompute(ntR[0,:], ntC[i,:]) - beta
                    ind += 1
            
            # If the first interval of ntR crosses zero, duplicate it right for the final interval
            if ntR[0,0] < 0:
                if intervalIntersectionCheck(ntR[0,:] + 2*np.pi, ntC[1,:]):
                    ntRC[ind,:] = intervalIntersectionCompute(ntR[0,:] + 2*np.pi, ntC[1,:]) - beta
                    ind += 1
                    
        elif num_ntR == 2:
            # R1 to C1
            if intervalIntersectionCheck(ntR[0,:], ntC[0,:]):
                ntRC[ind,:] = intervalIntersectionCompute(ntR[0,:], ntC[0,:]) - beta
                ind += 1
                
            # Depending on where R1 ends, either increment to R2 or test R1 to C2.
            # This ensures we only ever have four intervals
            if ntR[0,1] <= np.pi:
                # Test R2 to C1
                if intervalIntersectionCheck(ntR[1,:], ntC[0,:]):
                    ntRC[ind,:] = intervalIntersectionCompute(ntR[1,:], ntC[0,:]) - beta
                    ind += 1
            else:
                # Test R1 to C2
                if intervalIntersectionCheck(ntR[0,:], ntC[1,:]):
                    ntRC[ind,:] = intervalIntersectionCompute(ntR[0,:], ntC[1,:]) - beta
                    ind += 1
                    
            # R2 to C2
            if intervalIntersectionCheck(ntR[1,:], ntC[1,:]):
                ntRC[ind,:] = intervalIntersectionCompute(ntR[1,:], ntC[1,:]) - beta
                ind += 1
                
            # If the first interval of ntR crosses zero, duplicate it right for the final interval
            if ntR[0,0] < 0:
                if intervalIntersectionCheck(ntR[0,:] + 2*np.pi, ntC[1,:]):
                    ntRC[ind,:] = intervalIntersectionCompute(ntR[0,:] + 2*np.pi, ntC[1,:]) - beta
                    ind += 1
        num_ntRC = ind
        
        # Fix the shift
        ntR -= beta
        ntC -= beta
    
    return ntRCFull, num_ntRC, ntRC
            
    
        
def intervalIntersectionCompute(a,b):
    """
    Computes the intersection of two intervals 
    """
    return np.array([max(a[0],b[0]),min([a[1],b[1]])])
        
    
def intervalIntersectionCheck(a,b):
    """
    Checks for the intersection of two intervals 
    """
    
    if b[0] > a[1] or a[0] > b[1]:
        return False
    else:
        return True


def clroeProp(L, n, oth, tStep, kovLim, kovType = 'PRISM'):
    """
    Assesses passive safety using CLROE propagation

    Parameters
    ----------
    L : 6x1 double
        CLROE state
    n : double
        Chief orbit mean motion
    oth : double
        Operational time horizon (sec)
    tStep: double
        time step
    KovLim : nx1 double
        Keep out volume limits.
    KovType : TYPE, optional
        Keep out volume type. The default is 'PRISM'.

    Returns
    -------
    psa : boolean
        Passive safety assessment

    """
    
    t = 0 
    
    while t < oth:
        t = t + tStep       
        rho = clroeEval(L, n*t)
        kovBreach = checkKovBreach(rho, kovLim, kovType)
        if kovBreach:
            return False
        
    return True
                    
        
def clroeEval(L, nt):
    """
    Evaluate CLROE state

    Parameters
    ----------
    L : 6x1 double
        CLROE state
    nt : double

    Returns
    -------
    rho

    """
    
    rho = np.zeros((3,))
    
    #A0 = L[0]
    #alpha = L[1]
    #xOff = L[2]
    #yOff = L[3]
    #B0 = L[4]
    #beta = L[5]
    
    rho[0] = L[0]*np.cos(nt + L[1]) + L[2]
    rho[1] = -2*L[0]*np.sin(nt + L[1]) + L[3] - 1.5*nt*L[2]
    rho[2] = L[4]*np.cos(nt + L[5])
    
    return rho    

def hcwProp(r, v, n, oth, tStep, kovLim, kovType = 'PRISM'):
    """
    Assesses passive safety using HCW propagation

    Parameters
    ----------
    r : 3x1 double
        Initial relative position in RIC
    v : 3x1 double
        Initial relative velocity in RIC
    n : double
        Chief orbit mean motion
    oth : double
        Operational time horizon (sec)
    tStep: double
        time step
    KovLim : nx1 double
        Keep out volume limits.
    KovType : TYPE, optional
        Keep out volume type. The default is 'PRISM'.

    Returns
    -------
    psa : boolean
        Passive safety assessment

    """
    
    t = 0 
    
    stm = lm.hcw(n).Phi(0, tStep)
    
    while t < oth:
        t = t + tStep    
        r, v = lm.multPhi(stm, r, v)
        kovBreach = checkKovBreach(r, kovLim, kovType)
        if kovBreach:
            return False
        
    return True
    
    
def checkKovBreach(rho, kovLim, kovType = 'PRISM'):
    """
    Checks whether relative position is in KOV

    Parameters
    ----------
    rho : 3x1 double
        relative RIC position.
    kovLim : nx1 double
        Keep out volume.
    kovType : string, optional
        KOV type. The default is 'PRISM'.

    Returns
    -------
    boolean.

    """
    if kovType == "PRISM":
        if all(np.greater(kovLim,np.abs(rho))) == 1:
            return True
    elif kovType == "SPHERE":
        if np.linalg.norm(rho) < kovLim:
            return True
    elif kovType == "ELIPSOID":
        # TODO
        return True
    
    return False
    