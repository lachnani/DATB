# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 08:43:04 2026

@author: Hakim Lachnani
"""

import numpy as np

"""
Feedback control drives the relative state to the guidance trajectory. The 
guidance trajectory may be defined as any set of relative states (pos/vel or 
elements). The output of all functions is a 3x1 delta-v vector. The frame of
delta-v vector depends on the function. In most cases, the function is frame
agnostic, so the output delta-v is in the same frame as the inputs. Controllers
include:
    - pd: proportional-derivative (in RIC rect/curv or ECI)

"""

def pd(r, v, rDes, vDes, Kr, Kv):
    """
    Proportional-Derivative Feedback Control

    Parameters
    ----------
    r : 3x1 double
        current position
    v : 3x1 double
        current velocity
    rDes : 3x1 double
        desired position
    vDes : 3x1 double
        desired velocity
    Kr : 3x3 double
        Proportional gain
    Kv : 3x3 double
        Derivative gain

    Returns
    -------
    dv: 3x1 double
        control delta-v

    """
    return np.matmul(Kr, rDes-r) + np.matmul(Kv, vDes-v)
    

