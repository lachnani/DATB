# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 10:50:54 2026

@author: Hakim Lachnani
"""

import numpy as np

def abaT(a, b):
    return np.matmul(a,np.matmul(b,np.transpose(a))) 

def rotateVec(dcm, v):
    return dcm @ v

def rotateMat(dcm, M):
    return abaT(dcm, M)