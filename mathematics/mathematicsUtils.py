# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 10:50:54 2026

@author: Hakim Lachnani
"""

import numpy as np
from numpy import linalg as la

def abaT(a, b):
    return np.matmul(a,np.matmul(b,np.transpose(a))) 

def rotateVec(dcm, v):
    return dcm @ v

def rotateMat(dcm, M):
    return abaT(dcm, M)

def vecMag(v):
    return la.norm(v)

def unitVec(v):
    vMag = vecMag(v)
    return v/vMag, vMag

def zeroVec(n = 3):
    return np.zeros(n,)