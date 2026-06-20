# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 12:21:16 2026

@author: Hakim Lachnani
"""

import numpy as np

"""
Measurement sensitivity on state. Measurement types include:
    - angles: Azimuth and Elevation angles (as measured in the RIC frame)
    - relRange: Range
    - anglesRange: Azimuth, Elevation, and Range
    - relRangeRate: Range and Range-Rate
    - anglesRangeRate: Azimuth, Elevation, Range, and Range-Rate
    - cv3dof: Relative Position (in RIC)
    - cv6dof: Relative Position and Velocity (in RIC)

State types include:
    - dualInertial: 12 DOF inertial (rc, vc, rd, vd)
    - rectRic: 6 DOF relative rectilinear RIC (relPosRectRic, relVelRectRic)
    - curvRic: 6 DOF relative curvilinear RIC (relPosCurvRic, relVelCurvRic)
    
Each function returns the measurement sensitivity matrix H. The input to the 
function is the estimated state. H is (numMeas) x (numStates)

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
        
def relRange(x):
    return 0.0

def anglesRange(x):
    return 0.0 
        
def relRangeRate(x):
    return 0.0 
        
def anglesRangeRate(x):
    return 0.0
        
def cv3dof_rectRic(x = None):
    return np.eye(3,6)

def cv6dof_rectRic(x = None):
    return np.eye(6)
              
    