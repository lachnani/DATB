# -*- coding: utf-8 -*-
"""
Created on Fri Jan 16 17:22:01 2026

@author: Hakim Lachnani
"""

from guidance import linearizedModels as lm
import numpy as np
from numpy import linalg as la

def twoImpulseBurn(stm, r0, v0, rf, vf):
    """
    Computes RIC frame burn vectors for two-impluse rendezvous. 
    
    Separates the in-plane and cross-track burn components to avoid 
    singularities.

    Parameters
    ----------
    stm : State transition matrix (HCW, sGA, GA, or YA).
    r0 : initial relative position in RIC.
    v0 : initial relative velocity in RIC.
    rf : final relative position in RIC.
    vf : final relative velocity in RIC.

    Returns
    -------
    dv0 : initial delta-v in RIC.
    dvf : final delta-v in RIC.

    """
    
    # Initialize Outputs
    dv0 = np.zeros((3,))
    dvf = np.zeros((3,))
    
    # Split the STM into in-plane and cross track components
    Phi_IP, Phi_OP = splitStm(stm)
    
    # In-plane 
    Phi_IP_rr = Phi_IP[0:2,0:2]
    Phi_IP_rv = Phi_IP[0:2,2:4]
    Phi_IP_vr = Phi_IP[2:4,0:2]
    Phi_IP_vv = Phi_IP[2:4,2:4]
    
    if la.det(Phi_IP_rv) > np.finfo(np.float64).eps:
        dv0_IP = np.matmul(la.inv(Phi_IP_rv),(rf[0:2] - np.matmul(Phi_IP_rr,r0[0:2]))) - v0[0:2]
        dvf_IP = vf[0:2] - np.matmul(Phi_IP_vr,r0[0:2]) - np.matmul(Phi_IP_vv,(v0[0:2]+dv0_IP))
    else:
        dv0_IP = np.zeros((2,))
        dvf_IP = np.zeros((2,))
        
    # Out of plane
    if np.abs(Phi_OP[0,1]) > np.finfo(np.float64).eps:
        dv0_OP = (rf[2] - Phi_OP[0,0]*r0[2])/Phi_OP[0,1] - v0[2]
        dvf_OP = vf[2] - Phi_OP[1,0]*r0[2] - Phi_OP[1,1]*(v0[2]+dv0_OP)
    else:
        dv0_OP = 0.0;
        dvf_OP = 0.0;
        
    # Populate Outputs
    dv0[0:2] = dv0_IP
    dv0[2] = dv0_OP
    dvf[0:2] = dvf_IP
    dvf[2] = dvf_OP
    
    return dv0, dvf

def wptBurn(stm, r0, v0, rf):
    """
    Computes RIC frame burn vector for a position targeting waypoint burn. 
    
    Separates the in-plane and cross-track burn components to avoid 
    singularities.

    Parameters
    ----------
    stm : State transition matrix (HCW, sGA, GA, or YA).
    r0 : initial relative position in RIC.
    v0 : initial relative velocity in RIC.
    rf : final relative position in RIC.

    Returns
    -------
    dv0 : initial delta-v in RIC.

    """
    
    # Initialize Outputs
    dv0 = np.zeros((3,))
    
    # Split the STM into in-plane and cross track components
    Phi_IP, Phi_OP = splitStm(stm)
    
    # In-plane 
    Phi_IP_rr = Phi_IP[0:2,0:2]
    Phi_IP_rv = Phi_IP[0:2,2:4]
    
    if la.det(Phi_IP_rv) > np.finfo(np.float64).eps:
        dv0_IP = np.matmul(la.inv(Phi_IP_rv),(rf[0:2] - np.matmul(Phi_IP_rr,r0[0:2]))) - v0[0:2]
    else:
        dv0_IP = np.zeros((2,))
        
    # Out of plane
    if np.abs(Phi_OP[0,1]) > np.finfo(np.float64).eps:
        dv0_OP = (rf[2] - Phi_OP[0,0]*r0[2])/Phi_OP[0,1] - v0[2]
    else:
        dv0_OP = 0.0;
        
    # Populate Outputs
    dv0[0:2] = dv0_IP
    dv0[2] = dv0_OP
    
    return dv0

def wptBurnTable(stmClass, wptTbl, t0, r0, v0, lastAction = 'IDLE'):
    """
    Computes series of burns for a waypoint table

    Parameters
    ----------
    stmClass : State transition matrix class (HCW, sGA, GA, or YA).
    wptTbl : Waypoint table class.
    t0 : initial time
    r0 : initial relative position in RIC.
    v0 : initial relative velocity in RIC.

    Returns
    -------
    dvTbl : nominal waypoint burn table.

    """
    
    # Initialize Output
    # Number of burns is equal to the number of waypoints plus one. The first 
    # burn is at the initial time targeting the first waypoint, etc
    burnTable = np.zeros((3,wptTbl.numWpts+1))
    
    # First Burn
    stm = stmClass.Phi(t0, wptTbl.t[0])
    burnTable[:,0] = wptBurn(stm, r0, v0, wptTbl.relPosRic[:,0])
    
    # Propagated velocity at second waypoint
    rProp, vProp = lm.multPhi(stm, r0, v0 + burnTable[:,0])
    
    # Sanity Check that rProp matches the next waypoint position
    if ~np.all(np.isclose(rProp, wptTbl.relPosRic[:,0])):
        raise ValueError("wptGuidance.py: initial waypoint burn calculation error")
    
    # Second Burn to Penultimate Burn
    for burn in range(1, wptTbl.numWpts):
        stm =  stmClass.Phi(wptTbl.t[burn-1], wptTbl.t[burn])
        burnTable[:,burn] = wptBurn(
            stm, 
            wptTbl.relPosRic[:,burn-1], 
            vProp, 
            wptTbl.relPosRic[:,burn])
        
        # Propagated velocity at next waypoint
        rProp, vProp = lm.multPhi(stm, wptTbl.relPosRic[:,burn-1], vProp + burnTable[:,burn])
        
        # Sanity Check that rProp matches the next waypoint position
        if ~np.all(np.isclose(rProp, wptTbl.relPosRic[:,burn])):
            raise ValueError("wptGuidance.py: incremental waypoint burn calculation error")
            
    # Last burn
    if lastAction == 'IDLE':
        # Coast
        burnTable[:,-1] = np.zeros((3,))
    elif lastAction == 'STOP':
        burnTable[:,-1] = -vProp
        
    return np.insert(wptTbl.t, 0, t0, axis=0), burnTable
    
    
    
def splitStm(stm):
    """
    Splits state transition matrix into in-plane and cross-track components

    Parameters
    ----------
    stm : State transition matrix (6x6).

    Returns
    -------
    Phi_IP : In-Plane State transition matrix (4x4).
    Phi_OP : Out-of-Plane State transition matrix (2x2).

    """
    
    Phi_IP = np.array(((stm[0,0],stm[0,1],stm[0,3],stm[0,4]),
                       (stm[1,0],stm[1,1],stm[1,3],stm[1,4]),
                       (stm[3,0],stm[3,1],stm[3,3],stm[3,4]),
                       (stm[4,0],stm[4,1],stm[4,3],stm[4,4])))
    
    Phi_OP = np.array(((stm[2,2],stm[2,5]),
                       (stm[5,2],stm[5,5])))
    
    return Phi_IP, Phi_OP