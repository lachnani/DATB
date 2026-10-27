# -*- coding: utf-8 -*-
"""
Created on Sun Dec 14 16:35:34 2025

@author: Hakim Lachnani
"""

import numpy as np

class deltaVTable():
    """ 
    Base Delta-V table class
    
    """
    
    def __init__(self):
        self.numDv = 0
        self.t = np.empty( shape=(0, 0) )
        self.dv = np.empty( shape=(3, 0) )
        self.frame = []
        
    def increment(self):
        if self.numDv > 0:
            self.numDv += -1
            self.t = self.t[1:]
            self.dv = self.dv[:,1:]
            self.frame = self.frame[1:]
        
    def appendDvs(self, dvTimes, dvVectors, dvFrames):
        self.t = np.append(self.t, dvTimes)
        self.dv = np.append(self.dv, dvVectors, axis=1)
        self.frame.append(dvFrames)
        self.numDv += dvTimes.size
        
    def truncate(self, numDv):
        self.t = self.t[0:numDv]
        self.dv = self.dv[:,0:numDv]
        self.frame = self.frame[0:numDv]
        self.numDv = numDv
        
    def asInr(self, dvIdx = 0, dcmInrToRic = np.zeros((3,3))):
        """
        Returns inertial delta-v vector
        """
        if self.frame[dvIdx].upper() == 'INR':
            # Return DV as is
            return self.dv[:,dvIdx]
        elif self.frame[dvIdx].upper() == 'RIC':
            # Rotate from RIC frame to Inertial frame
            return dcmInrToRic.T @ self.dv[:,dvIdx]
        else:
            return np.zeros(3,)
        
    def asRic(self, dvIdx = 0, dcmInrToRic = np.zeros((3,3))):
        """
        Returns RIC delta-v vector
        """
        if self.frame[dvIdx].upper() == 'RIC':
            # Return DV as is
            return self.dv[:,dvIdx]
        elif self.frame[dvIdx].upper() == 'INR':
            # Rotate from Inertial frame to RIC frame
            return dcmInrToRic @ self.dv[:,dvIdx]
        else:
            return np.zeros(3,)