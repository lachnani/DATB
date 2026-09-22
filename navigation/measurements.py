# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 12:21:16 2026

@author: Hakim Lachnani
"""

import numpy as np
from numpy import random as rand

# Define measurement types and their indexing
measType = {
    "none":             np.array([0,0,0,0], dtype='bool'),
    "angles":           np.array([1,1,0,0], dtype='bool'),
    "range":            np.array([0,0,1,0], dtype='bool'),
    "rangeRate":        np.array([0,0,0,1], dtype='bool'),
    "anglesRange":      np.array([1,1,1,0], dtype='bool'),
    "rangeRR":          np.array([0,0,1,1], dtype='bool'),
    "anglesRangeRR":    np.array([1,1,1,1], dtype='bool')
    }

"""
Measurements derived from truth states and given covariances. Measurement types
include azimuth and elevation angles (as measured in the LOS frame), range, and
range-rate
"""
        
def get(frm, R):
    return np.array([frm.az,frm.el,frm.rng,frm.rngRate]) + \
        rand.multivariate_normal(np.zeros(4,),R) 

def inertialMeasurementSensitivity(fltr):
    """
    Calculates measurement sensitivity in the inertial frame. Returns a 4x6 
    matrix with the positive convention.
    
    Ref: Woffinden, David Charles, "Angles-Only Navigation for Autonomous 
    Orbital Rendezvous" (2008). All Graduate Theses and Dissertations. 12.
    https://digitalcommons.usu.edu/etd/12

    Parameters
    ----------
    fltr : class
        filter class; see estimator.py

    Returns
    -------
    H : 4x6 double
        measurement sensitivity matrix

    """
    ### Azimuth and Elevation
    # Partials of the LOS unit vector wrt az and el in LOS frame ([1] Eq 6.86)
    P_az_los = np.array([-np.cos(fltr.el)*np.sin(fltr.az),np.cos(fltr.el)*np.cos(fltr.az),0.0])
    P_el_los = np.array([-np.sin(fltr.el)*np.cos(fltr.az),-np.sin(fltr.el)*np.sin(fltr.az),0.0])
    # Partials of the LOS unit vector wrt az and el in Inr frame 
    P_az_inr = np.matmul(np.transpose(fltr.dcmInr2Los),P_az_los)
    P_el_inr = np.matmul(np.transpose(fltr.dcmInr2Los),P_el_los)
    # Partials with respect to position ([1] Eq 6.104)
    H_az_r = np.transpose(P_az_inr)/(fltr.rng*np.cos(fltr.el)**2)
    H_el_r = np.transpose(P_el_inr)/fltr.rng
    ### Range 
    # Unit LOS vector in the inertial frame
    losUnitInr = (fltr.chiefPosInr - fltr.deputyPosInr)/fltr.rng
    # Partial with respect to position ([1] Eq 6.104)
    H_rng_r = np.transpose(losUnitInr)
    ### Range-Rate
    # Partial with respect to position 
    H_rngRate_r = np.transpose(fltr.chiefVelInr - fltr.deputyVelInr - fltr.rngRate*losUnitInr)/fltr.rng
    # Partial with respect to velocity
    H_rngRate_v = np.transpose(losUnitInr)
    # Construct full H matrix
    return np.block([
                    [H_az_r,      np.zeros((1,3))],
                    [H_el_r,      np.zeros((1,3))],
                    [H_rng_r,     np.zeros((1,3))],
                    [H_rngRate_r, H_rngRate_v    ]])
              
def calcAzEl(rc, rd, dcmInr2Los):
    losInr = rc - rd
    losUnitInr = losInr / np.linalg.norm(losInr)
    losUnitLos = np.matmul(dcmInr2Los, losUnitInr)
    az = np.arctan2(losUnitLos[1],losUnitLos[0])
    el = np.arcsin(losUnitLos[2])
    return az, el
    