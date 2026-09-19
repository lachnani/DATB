# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 19:52:45 2026

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
    Q: function
        process noise funtion Q(dt)
    f: function
        state update function f(dt, x, u, param)
    F: function
        state transition matrix function F(dt, x)
    S: nxn float
        state correction covariance matrix
        
    Variables are based on the following:
        1. Yaakov Bar-Shalom, X.-Rong Li, Thiagalingam Kirubarajan, "Estimation 
        with Applications to Tracking and Navigation: Theory, Algorithms and 
        Software" (2002). 
        https://onlinelibrary.wiley.com/doi/book/10.1002/0471221279
        2. Crassidis and Junkins, "Optimal Estimation of Dynamic Systems" 
        (2012).
        https://www.routledge.com/Optimal-Estimation-of-Dynamic-Systems/Crassidis-Junkins/p/book/9781032917610
        
    
    """
    def __init__(
            self, 
            t, x, P, Q, f, Phi, Qu, param = None
            ):
        
        # Initialize the filter state
        self.t = t          # time
        self.x = x          # state
        self.P = P          # covariance
        self.Q = Q          # process noise function
        self.f = f          # state equation function
        self.Phi = Phi      # state transition matrix function
        self.Qu = Qu        # maneuver process noise function
        self.param = param  # assitional parameters
        self.n = np.size(x,0) # state size
        
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
        self.x = self.f(dt, self.x, u, self.param)
        if la.norm(u) > 0.0:
            self.P = propagateCov(self.P, self.Phi(dt, self.x), self.Q(dt) + self.Qu(dt,u))
        else:
            self.P = propagateCov(self.P, self.Phi(dt, self.x), self.Q(dt))
            
    def update(self, y, H, R):
        """
        Update state and covariance with measurements z. Follows Figure 
        10.3.3-1 of [1].

        Parameters
        ----------
        y : mx1 double
            measurement residual (z - zHat)
        H : mxn double
            measurement sensitivity matrix
        R : mxm double
            measurement covariance

        """
        self.x, self.P = measurementUpdate(self.x, self.P, self.n, y, H, R)
        
def propagateCov(P, Phi, Q):
    """
    Propagates covariance matrix and adds process noise using first order model

    Parameters
    ----------
    P : nxn double
        state covariance
    Phi : nxn double
        state transition matrix
    Q : nxn double
        process noise matrix

    Returns
    -------
    PProp : nxn double
        propagated state covariance

    """
    return symmetrize(np.matmul(Phi,np.matmul(P,np.transpose(Phi))) + Q)

def measurementUpdate(x, P, n, y, H, R):
    """
    Update state and covariance with measurements z. Follows Figure 
    10.3.3-1 of [1].

    Parameters
    ----------
    x : nx1 double
        state vector
    P : nxn double
        state covariance
    n : int
        state dimension
    y : mx1 double
        measurement residual (z - zHat)
    H : mxn double
        measurement sensitivity matrix
    R : mxm double
        measurement covariance

    Returns
    -------
    xUpd : nx1 double
        updated state vector
    PUpd : nxn double
        updated state covariance

    """
    S = residualCov(P, H, R)
    K = kalmanGain(P, H, S)
    xUpd = x + np.matmul(K, y)
    PUpd = np.matmul(np.eye(n) - np.matmul(K,H),P)
    PUpd = symmetrize(PUpd)
    return xUpd, PUpd
        
def residualCov(P, H, R):
    return np.matmul(H,np.matmul(P,np.transpose(H))) + R
        
def kalmanGain(P, H, S):
    return np.matmul(P,np.matmul(np.transpose(H),la.inv(S)))

def symmetrize(P):
    return 0.5*(P + np.transpose(P))

def normEstErrSqr(dx, P):
    return np.matmul(np.transpose(dx),np.matmul(np.inv(P),dx))

def normInvnSqr(nu, S):
    return np.matmul(np.transpose(nu),np.matmul(np.inv(S),nu))