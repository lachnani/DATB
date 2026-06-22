# -*- coding: utf-8 -*-
"""
Created on Sun Dec 14 16:35:34 2025

@author: Hakim Lachnani
"""

import numpy as np
from kinematics import kinematicsUtils as uKin

class waypointTable():
    """ 
    Base waypoint table class
    
    """
    
    def __init__(self):
        self.numWpts = 0
        self.t = np.empty( shape=(0, 0) )
        self.relPosRic = np.empty( shape=(3, 0) )
        
    def increment(self):
        if self.numWpts > 0:
            self.numWpts += -1
            self.t = self.t[1:]
            self.relPosRic = self.relPosRic[:,1:]
        
    def appendWpts(self, wptTimes, wptPos):
        self.t = np.append(self.t, wptTimes)
        self.relPosRic = np.append(self.relPosRic, wptPos, axis=1)
        self.numWpts += wptTimes.size
        
    def truncate(self, numWpts):
        self.t = self.t[0:numWpts]
        self.relPosRic = self.relPosRic[:,0:numWpts]
        self.numWpts = numWpts
        
    def genClroeWpts(self, L, n, t0, tf, numWpts):
        """
        Generates waypoints according to CLROEs
        """
        wptTimes = np.linspace(t0, tf, numWpts)
        wptPos = np.zeros((3,numWpts))
        wptVel = np.zeros((3,numWpts))
        relPosRectRic = np.zeros((3,))
        relVelRectRic = np.zeros((3,))
        for wpt in range(numWpts):
            uKin.clroe2ric(L, n, wptTimes[wpt], relPosRectRic, relVelRectRic) 
            wptPos[:,wpt] = relPosRectRic
            wptVel[:,wpt] = relVelRectRic
        self.t = np.append(self.t, wptTimes)
        self.relPosRic = np.append(self.relPosRic, wptPos, axis=1)
        self.numWpts += numWpts
        
    def genRptWpts(self, pos, t0, dt, numWpts):
        """
        Generates repeated waypoints
        """
        for wpt in range(numWpts):
            self.t = np.append(self.t, t0 + wpt*dt)
            pos = np.reshape(pos, (-1, 1))
            if self.numWpts == 0:
                self.relPosRic = pos
            else:
                self.relPosRic = np.append(self.relPosRic, pos, axis=1)
            self.numWpts += 1
            
    def genFmcWpts(self, L, n, t0, speedSf, numWpts, numRevs=1, constantRange=False):
        """
        Generates FMC Waypoints
        """
        # Force xOff = 0
        L[2] = 0.0 
        if constantRange:
            L[3] = 0.0
        if numWpts > 2:
            # Need at least two waypoints for Cfmc
            revTime = 2*np.pi/n/speedSf
            dt = revTime/(numWpts - 1)
            pos = np.zeros((3,))
            vel = np.zeros((3,))
            for rev in range(numRevs):
                for wpt in range(numWpts):
                    t = t0 + wpt*dt + rev*revTime
                    self.t = np.append(self.t, t)
                    uKin.clroe2ric(L, n*speedSf, t, pos, vel)    
                    if constantRange:
                        pos = L[0]*pos/np.linalg.norm(pos)
                    self.relPosRic = np.append(self.relPosRic, np.reshape(pos, (-1, 1)), axis=1)
                    self.numWpts += 1
