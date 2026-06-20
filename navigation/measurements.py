# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 12:21:16 2026

@author: Hakim Lachnani
"""

import numpy as np
from numpy import random as rand

"""
Measurements derived from truth states and given covariances. Measurement types
include:
    - angles: Azimuth and Elevation angles (as measured in the RIC frame)
    - relRange: Range
    - anglesRange: Azimuth, Elevation, and Range
    - relRangeRate: Range and Range-Rate
    - anglesRangeRate: Azimuth, Elevation, Range, and Range-Rate
    - cv3dof: Relative Position (in RIC)
    - cv6dof: Relative Position and Velocity (in RIC)
    
For each function, the inputs are the formation struct frm, and the measurement 
covariance R.

"""

def angles(frm, R):
    return np.array([frm.az,frm.el]) + \
        rand.multivariate_normal(np.zeros(2,),R) 
        
def relRange(frm, var):
    return frm.rng + rand.normal(0, np.sqrt(var))

def anglesRange(frm, R):
    return np.array([frm.az,frm.el,frm.rng]) + \
        rand.multivariate_normal(np.zeros(3,),R) 
        
def relRangeRate(frm, R):
    return np.array([frm.rng,frm.rngRate]) + \
        rand.multivariate_normal(np.zeros(2,),R) 
        
def anglesRangeRate(frm, R):
    return np.array([frm.az,frm.el,frm.rng,frm.rngRate]) + \
        rand.multivariate_normal(np.zeros(4,),R) 
        
def cv3dof(frm, R):
    return frm.relPosRectRic + rand.multivariate_normal(np.zeros(3,),R) 

def cv6dof(frm, R):
    return np.array([frm.relPosRectRic,frm.relVelRectRic]) + \
        rand.multivariate_normal(np.zeros(2,),R)  
              
    