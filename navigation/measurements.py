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
    - angles: Azimuth and Elevation angles (as measured in the LOS frame)
    - relRange: Range
    - anglesRange: Azimuth, Elevation, and Range
    - relRangeRate: Range and Range-Rate
    - anglesRangeRate: Azimuth, Elevation, Range, and Range-Rate
    - cv3dof: Relative Position (in RIC)
    - cv6dof: Relative Position and Velocity (in RIC)
    - pnt: Inertial deputy PNT
    
For each function, the inputs are the simulation struct sim, and the 
measurement covariance R.

"""

def get(sim, measType, R):
    if measType == "angles":
        return angles(sim, R)
    elif measType == "relRange":
        return relRange(sim, R)
    elif measType == "anglesRange":
        return anglesRange(sim, R)
    elif measType == "relRangeRate":
        return relRangeRate(sim, R)
    elif measType == "anglesRangeRate":
        return anglesRangeRate(sim, R)
    elif measType == "cv3dof":
        return cv3dof(sim, R)
    elif measType == "cv6dof":
        return cv6dof(sim, R) 
    elif measType == "pnt":
        return pnt(sim, R)
    return

def angles(sim, R):
    return np.array([sim.az,sim.el]) + \
        rand.multivariate_normal(np.zeros(2,),R) 
        
def relRange(sim, var):
    return sim.frm.rng + rand.normal(0, np.sqrt(var))

def anglesRange(sim, R):
    return np.array([sim.az,sim.el,sim.frm.rng]) + \
        rand.multivariate_normal(np.zeros(3,),R) 
        
def relRangeRate(sim, R):
    return np.array([sim.frm.rng,sim.frm.rngRate]) + \
        rand.multivariate_normal(np.zeros(2,),R) 
        
def anglesRangeRate(sim, R):
    return np.array([sim.az,sim.el,sim.frm.rng,sim.frm.rngRate]) + \
        rand.multivariate_normal(np.zeros(4,),R) 
        
def cv3dof(sim, R):
    return sim.frm.relPosRectRic + rand.multivariate_normal(np.zeros(3,),R) 

def cv6dof(sim, R):
    return np.array([sim.frm.relPosRectRic,sim.frm.relVelRectRic]) + \
        rand.multivariate_normal(np.zeros(2,),R)  
        
def pnt(sim, R):
    return sim.frm.deputy.r + rand.multivariate_normal(np.zeros(3,),R) 
              
def calcAzEl(rc, rd, dcmInr2Los):
    losInr = rc - rd
    losUnitInr = losInr / np.linalg.norm(losInr)
    losUnitLos = np.matmul(dcmInr2Los, losUnitInr)
    az = np.arctan2(losUnitLos[1],losUnitLos[0])
    el = np.arcsin(losUnitLos[2])
    return az, el
    