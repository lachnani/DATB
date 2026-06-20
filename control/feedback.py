# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 08:43:04 2026

@author: Hakim Lachnani
"""

import numpy as np

def pd(r, v, rDes, vDes, Kr, Kv):
    return np.matmul(Kr, rDes-r) + np.matmul(Kv, vDes-v)
    

