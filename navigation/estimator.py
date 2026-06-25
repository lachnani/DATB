# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 23:18:57 2026

@author: Hakim Lachnani
"""

import numpy as np
from numpy import linalg as la

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
    F: function
        state transition matrix function F(dt, x)
    S: nxn float
        state correction covariance matrix
        
    Variables are based on the following:
        1. Yaakov Bar-Shalom, X.-Rong Li, Thiagalingam Kirubarajan, "Estimation 
        with Applications to Tracking and Navigation: Theory, Algorithms and 
        Software" (2002). 
        https://onlinelibrary.wiley.com/doi/book/10.1002/0471221279
        2. Woffinden, David Charles, "Angles-Only Navigation for Autonomous 
        Orbital Rendezvous" (2008). All Graduate Theses and Dissertations. 12.
        https://digitalcommons.usu.edu/etd/12
    
    """
    def __init__(
            self, 
            t, x, P, Q, f, F, S
            ):
        
        # Initialize the filter state
        self.t = t
        self.x = x
        self.P = P
        self.Q = Q
        self.f = f
        self.F = F
        self.S = S
        self.n = np.size(x,0)
        
    def propagate(self, dt, u):
        """
        Propagate state and covariance using dyamic equations. Follows Figure 
        10.3.3-1 of [1]. If the control force is non-zero, additional velocity 
        covariance is added per [2].

        Parameters
        ----------
        dt : float
            propagation delta-time
        u : 3x1 float
            control acceleration

        """
        self.t = self.t + dt
        F = self.F(dt, self.x)
        self.x = self.f(dt, self.x, u)
        if la.norm(u) > 0.0:
            self.P = np.matmul(F,np.matmul(self.P,np.transpose(F))) + self.Q + self.S
        else:
            self.P = np.matmul(F,np.matmul(self.P,np.transpose(F))) + self.Q
            
    def update(self, nu, H, R):
        """
        Update state and covariance with measurements z. Follows Figure 
        10.3.3-1 of [1].

        Parameters
        ----------
        nu : mx1 float
            measurement residual (z - zHat)
        H : mxn float
            measurement sensitivity matrix
        R : mxm float
            measurement covariance

        """
        S = residualCov(self.P, H, R)
        W = kalmanGain(self.P, H, S)
        self.x = self.x + np.matmul(W, nu)
        self.P = self.P - np.matmul(W,np.matmul(S,np.transpose(W)))
        
def residualCov(P, H, R):
    return np.matmul(H,np.matmul(P,np.transpose(H))) + R
        
def kalmanGain(P, H, S):
    return np.matmul(P,np.matmul(np.transpose(H),la.inv(S)))
