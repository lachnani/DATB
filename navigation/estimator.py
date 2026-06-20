# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 23:18:57 2026

@author: Hakim Lachnani
"""

import numpy as np

class ExtendedKalmanFilter:
    """ Base EKF class
    
    Attributes
    ----------
    t : float
        time 
    x: nx1 float
        state
    P: nxn float
        covariance
    Q: nxn float
        process noise
    f: function
        state update function f(dt, x, u)
    Phi: function
        state transition matrix function Phi(dt)
    d: function
        state correction function d(dt, x, u)
    D: function
        state correction matrix function D(dt, x)
    S: nxn float
        state correction covariance matrix
        
    Variables are based on the following:
        1. Woffinden, David Charles, "Angles-Only Navigation for Autonomous 
        Orbital Rendezvous" (2008). All Graduate Theses and Dissertations. 12.
        https://digitalcommons.usu.edu/etd/12
    
    """
    def __init__(
            self, 
            t, x, P, Q, f, Phi, d, D, S
            ):
        
        # Initialize the filter state
        self.t = t
        self.x = x
        self.P = P
        self.Q = Q
        self.f = f
        self.Phi = Phi
        self.d = d
        self.D = D
        self.S = S
        self.n = np.size(x,0)
        
    def propagate(self, dt, u):
        """
        Propagate state and covariance using dyamic equations

        Parameters
        ----------
        dt : float
            propagation delta-time
        u : 3x1 float
            control acceleration

        """
        self.t = self.t + dt
        self.x = self.f(dt, self.x, u)
        Phi = self.Phi(dt)
        self.P = np.matmul(Phi,np.matmul(self.P,np.transpose(Phi))) + self.Q
        
    def correct(self, dt, u):
        """
        Correct state and covariance for impulsive control u

        Parameters
        ----------
        dt : float
            propagation delta-time
        u : 3x1 float
            impulsive control force

        """
        self.x = self.x + self.d(dt, self.x, u)
        Id = np.eye(self.n) + self.D(dt, self.x)
        self.P = self.P + \
            np.matmul(Id,np.matmul(self.P,np.transpose(Id))) + self.Q
            
    def update(self, z, zHat, H, R):
        """
        Update state and covariance with measurements z

        Parameters
        ----------
        z : mx1 float
            measurement
        zHat : mx1 float
            expected measurement
        H : mxn float
            measurement sensitivity matrix
        R : mxm float
            measurement covariance

        """
        K = kalmanGain(self.P, H, R)
        self.x = self.x + np.matmul(K, z-zHat)
        self.P = np.matmul(np.eye(self.n)-np.matmul(K,H),self.P)
        
def kalmanGain(P, H, R):
    return np.matmul(P,np.matmul(np.transpose(H),\
                     np.inv(np.matmul(H,np.matmul(P,np.transpose(H))) + R)))
