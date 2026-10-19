# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 19:52:45 2026

@author: Hakim Lachnani
"""

import sys
import numpy as np
from numpy import linalg as la

COND_LIM = 1.0 / sys.float_info.epsilon

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
        self.x, self.P, self.nis = measurementUpdate(self.x, self.P, self.n, y, H, R)
        
class Buffer:
    """
    Kalman Filter Buffer Class. Stores data for consistency testing.
    
    Attributes
    ----------
    n : integer
        buffer depth 
    n : integer
        buffer head 
    full : boolean
        buffer full status
    t : n double
        time history
    est: nState x n double
        state estimate history
    res: nMeas x n double
        measurement residual history
    nis: n double
        normalized innovations squared history
        
    """
    def __init__(self, nBuffer, nMeas, nState):
        """
        Initialize Buffer.

        Parameters
        ----------
        nBuffer : integer
            buffer depth.
        nMeas : integer
            measurement dimension.
        nState : integer
            state dimension.

        """
        
        self.n = nBuffer
        self.head = 0
        self.full = False
        self.t = np.zeros((nBuffer,))
        self.est = np.zeros((nState, nBuffer))
        self.res = np.zeros((nMeas, nBuffer))
        self.nis = np.zeros((nBuffer,))
        
    def update(self, t, est, res, nis):
        """
        Update Buffer.

        Parameters
        ----------
        t : double 
            time.
        est : nState x 1 double
            state estimate.
        res : nMeas x 1 double
            measurement residual.
        nis : double 
            normalized innovations squared.

        """
        
        # Roll buffer if full
        if self.full == True:
            self.t   = np.roll(self.t,   -1)
            self.est = np.roll(self.est, -1)
            self.res = np.roll(self.res, -1)
            self.nis = np.roll(self.nis, -1)
        
        # Update buffer in head index
        self.t[self.head] = t
        self.est[:,self.head] = est
        self.res[:,self.head] = res
        self.nis[self.head] = nis
        
        # Increment head if not full
        if self.full == False:
            self.head = self.head + 1 
            # Check to see if buffer is now full
            if self.head > self.n - 1:
                self.head = self.n - 1 
                self.full = True
    
        
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
    return symmetrize((abaT(Phi, P)) + Q)

def measurementUpdate(x, P, n, y, H, R, beta = 0):
    """
    Update state and covariance with residuals y. Follows Figure 
    10.3.3-1 of [1]. Uses the stable Joseph form update.
    
    Uderweighting performed using Lear's method. Ref: Zanetti, DeMars, and
    Bishop, "Underweighting Nonlinear Measurements" 
    https://arc.aiaa.org/doi/10.2514/1.50596

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
    beta : double
        underweighting factor

    Returns
    -------
    xUpd : nx1 double
        updated state vector
    PUpd : nxn double
        updated state covariance
    nisUpd : double
        measurement normalized innovarions squared

    """
    # Residual Covariance
    S = residualCov(P, H, R)
    # Condition check to ensure invertibility
    if la.cond(S) < COND_LIM:
        # NIS
        nis = normInvnSqr(y, S)
        # Underweighting factor
        HPHT = abaT(H, P)
        betaHPHT = beta*HPHT
        if beta > 0:
            S = S + betaHPHT
        # Kalman gain
        K = kalmanGain(P, H, S)
        # State Update
        xUpd = x + np.matmul(K, y)
        # Covariance Update
        I_KH = np.eye(n) - np.matmul(K,H)
        PUpd = abaT(I_KH, P) + abaT(K, R + betaHPHT)
        PUpd = symmetrize(PUpd)
        return xUpd, PUpd, nis
    else:
        return x, P, 0.0
    
def abaT(a, b):
    return np.matmul(a,np.matmul(b,np.transpose(a)))   
        
def residualCov(P, H, R):
    return abaT(H, P) + R
        
def kalmanGain(P, H, S):
    return np.matmul(P,np.matmul(np.transpose(H),la.inv(S)))

def symmetrize(P):
    return 0.5*(P + np.transpose(P))

def normEstErrSqr(dx, P):
    return np.matmul(np.transpose(dx),np.matmul(la.inv(P),dx))

def normInvnSqr(nu, S):
    return np.matmul(np.transpose(nu),np.matmul(la.inv(S),nu))

def checkPositiveDefinite(P, tol = 1e-10):
    # Step 1: Ensure it is a square matrix
    if P.shape[0] != P.shape[1]:
        return False
        
    # Step 2: Check for symmetry
    if not np.allclose(P, P.T):
        return False
        
    # Step 3: Calculate eigenvalues and verify they are strictly positive
    return np.all(la.eigvals(P) > tol)

def uduComposition(U, D):
    return U.dot(np.diag(D)).dot(U.T)

def uduDecomposition(P):
    """
    UDU^T Decomposition per the Bierman-Thornton filter.
    
    Ref: https://github.com/pykalman/pykalman/blob/main/pykalman/sqrt/bierman.py

    Parameters
    ----------
    P : nxn double
        Covariance.

    Returns
    -------
    U : nxn double
        Upper triangular matrix U.
    D : n double
        Vector of diagonal entries D.
    """
    n = P.shape[0]

    # Make a working copy and extract upper triangle
    P = np.triu(P).copy()
    U = np.eye(n)
    D = np.zeros(n)
  
    for j in reversed(range(2, n + 1)):
      D[j - 1] = P[j - 1, j - 1]
      if D[j - 1] > 0:
        alpha = 1.0 / D[j - 1]
      else:
        alpha = 0.0
  
      for k in range(1, j):
        beta = P[k - 1, j - 1]
        U[k - 1, j - 1] = alpha * beta
        P[0:k, k - 1] = P[0:k, k - 1] - beta * U[0:k, j - 1]
  
    D[0] = P[0, 0]
    return U, D

def uduScalarMeasurementUpdate(U, D, h, r, beta = 0):
    """
    Measurement update of the Bierman-Thornton filter for scalar measurements.
    
    Ref: https://github.com/pykalman/pykalman/blob/main/pykalman/sqrt/bierman.py
    
    TODO: Fix underweighting in covariance update!

    Parameters
    ----------
    U : nxn double
        Upper triangular matrix U from covariance UDU decomposition.
    D : n double
        Vector of diagonal entries D from covariance UDU decomposition.
    h : n double
        Row of sensitivity matrix.
    r : double
        measurement variance.
    beta : double
        underweighting factor

    Returns
    -------
    UBar : nxn double
        Updated upper triangular matrix U from covariance UDU decomposition.
    DBar : n double
        Updated vector of diagonal entries D from covariance UDU decomposition.
    k : double
        Kalman gain for the measurement.
    """
    
    n = len(h)

    f = h.dot(U)  # len(f) == n
    g = np.diag(D).dot(f)  # len(g) == n_dim-state
    alpha0 = f.dot(g)
    alpha = (1 + beta)*f.dot(g) + r 

    gamma = np.zeros(n)
    UBar = np.zeros((n, n))
    DBar = np.zeros(n)
    k = np.zeros(n)

    gamma[0] = r + g[0] * f[0]
    DBar[0] = D[0] * r / gamma[0]
    k[0] = g[0]

    UBar[0, 0] = 1
    for j in range(1, n):
        gamma[j] = gamma[j - 1] + g[j] * f[j]
        DBar[j] = D[j] * gamma[j - 1] / gamma[j]
        UBar[:, j] = U[:, j] - (f[j] / gamma[j - 1]) * k
        k = k + g[j] * U[:, j]
        
    k = k/alpha

    return UBar, DBar, k

def uduMeasurementUpdate(x, U, D, y, H, R, beta = 0):
    """
    Measurement update of the Bierman-Thornton filter for vector measurements.
    
    Ref: https://github.com/pykalman/pykalman/blob/main/pykalman/sqrt/bierman.py

    Parameters
    ----------
    x : nx1 double
        state vector
    U : nxn double
        Upper triangular matrix U from covariance UDU decomposition.
    D : n double
        Vector of diagonal entries D from covariance UDU decomposition.
    n : int
        state dimension
    y : mx1 double
        measurement residual (z - zHat)
    H : mxn double
        measurement sensitivity matrix
    R : mxm double
        measurement covariance
    beta : double
        underweighting factor

    Returns
    -------
    xUpd : nx1 double
        updated state vector
    UBar : nxn double
        Updated upper triangular matrix U from covariance UDU decomposition.
    DBar : n double
        Updated vector of diagonal entries D from covariance UDU decomposition.

    """
    
    # extract size of state space
    m = len(y)

    # calculate corrected state mean, covariance
    xUpd = x.copy()
    UUpd = U.copy()
    DUpd = D.copy() 
    for i in range(m):
        # extract components for updating i-th coordinate's covariance
        h = H[i]
        r = R[i, i]

        # calculate new UDU' decomposition for corrected state covariance
        # and a new column of the kalman gain
        UUpd, DUpd, k = uduScalarMeasurementUpdate(UUpd, DUpd, h, r, beta)

        # update corrected state mean
        xUpd = xUpd + k.dot(y[i])
        
    return xUpd, UUpd, DUpd