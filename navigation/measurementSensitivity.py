# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 12:21:16 2026

@author: Hakim Lachnani
"""

import numpy as np
from numpy import linalg as la

"""
Measurement sensitivity on state. Measurement types include:
    - angles: Azimuth and Elevation angles (as measured in the RIC frame)
    - relRange: Range
    - anglesRange: Azimuth, Elevation, and Range
    - relRangeRate: Range and Range-Rate
    - anglesRangeRate: Azimuth, Elevation, Range, and Range-Rate
    - cv3dof: Relative Position (in RIC)
    - cv6dof: Relative Position and Velocity (in RIC)
    - pnt: Inertial deputy PNT

State types include:
    - dualInertial: 12 DOF inertial (rc, vc, rd, vd)
    - rectRic: 6 DOF relative rectilinear RIC (relPosRectRic, relVelRectRic)
    - curvRic: 6 DOF relative curvilinear RIC (relPosCurvRic, relVelCurvRic)
    - deputyInertial: 6 DOF deputy inertial (rd, vd)
    
Each function returns the measurement sensitivity matrix H. The input to the 
function is the estimated state. H is (numMeas) x (numStates).

Note: Only pnt maps to the deputyInertial state

"""

def angles_dualInertial(x):
    return 0.0

def angles_rectRic(x):
    H = np.zeros((2,6))
    ip2 = x[0]**2 + x[1]**2
    if ip2 > np.eps:
        H[0,0] = -x[1]/ip2 # daz/dx
        H[0,1] = x[0]/ip2  # daz/dy
        ip = np.sqrt(ip2)
    rng2 = ip2 + x[2]**2 
    H[1,0] = -x[0]*x[2]/(rng2*ip) # del/dx
    H[1,1] = -x[1]*x[2]/(rng2*ip) # del/dy
    H[1,2] = ip/rng2 # del/dz
    return H
        
def relRange_rectRic(x):
    H = np.zeros((1,6))
    H[0,0:3] = np.transpose(x[0:3]/la.norm(x[0:3]))
    return H

def anglesRange_rectRic(x):
    H = np.zeros((3,6))
    H[0:2,:] = angles_rectRic(x)
    H[2,:] = relRange_rectRic(x)
    return H
        
def relRangeRate_rectRic(x):
    H = np.zeros((2,6))
    H[0,:] = relRange_rectRic(x)
    rng2 = x[0]**2 + x[1]**2 + x[2]**2
    rng = np.sqrt(rng2)
    rngRate = np.dot(x[0:3],x[3:6])/rng
    for ii in range(3):
        H[1,ii] = x[ii+3] - x[ii]*rngRate/rng2
        H[1,ii+3] = x[ii]/rng
    return H
        
def anglesRangeRate_rectRic(x):
    H = np.zeros((4,6))
    H[0:2,:] = angles_rectRic(x)
    H[2:4,:] = relRangeRate_rectRic(x)
    return H
        
def cv3dof_rectRic(x = None):
    return np.eye(3,6)

def cv6dof_rectRic(x = None):
    return np.eye(6)

def pnt_dualInertial(x = None):
    return np.eye(3,12,6)

def pnt_deputyInertial(x = None):
    return np.eye(3,6)
              
    