# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 23:18:57 2026

@author: Hakim Lachnani
"""

import numpy as np
from numpy import linalg as la
from dynamics import dynamicsUtils as uDyn
from dynamics import orbit as orb

class InertialEKF:
    """
    Inertial EKF is based on:
    Woffinden, David Charles, "Angles-Only Navigation for Autonomous 
    Orbital Rendezvous" (2008). All Graduate Theses and Dissertations. 12.
    https://digitalcommons.usu.edu/etd/12
    """
    
    def __init__(
            self,
            t, rc, vc, Pc, rd, vd, Pd, procVar, dvVar, sensor = None
            ):
        
        # Initialize the nav states
        self.t = t 
        self.rc = rc
        self.vc = vc
        self.Pc = Pc
        self.rd = rd
        self.vd = vd
        self.Pd = Pd
        self.procVar = procVar
        self.dvVar = dvVar
        self.sensor = sensor
        
        # Initialize the filter
        x = np.array([self.rc, self,vc, self.rd, self.vd])
        P = np.block([
                [Pc,               np.zeros((6, 6))],
                [np.zeros((6, 6)), Pd              ]])
        S = dvVar*np.eye(3)
        
        
    def stateUpdate(self, dt, x, u):
        """
        State update function. RK4 orbit propagator for chief and deputy 
        inertial states

        """
        rc = x[0:3]
        vc = x[3:6]
        rd = x[6:9]
        vd = x[9:12]
        uDyn.Orbit_rk4(self.pert["solarGrav"], self.pert["lunarGrav"], self.pert["drag"], self.pert["jnum"], \
                       self.moon.r, self.sun.r, self.Cd, self.normA, \
                       np.zeros((3,)), dt, rc, vc)
        uDyn.Orbit_rk4(self.pert["solarGrav"], self.pert["lunarGrav"], self.pert["drag"], self.pert["jnum"], \
                       self.moon.r, self.sun.r, self.Cd, self.normA, \
                       u, dt, rd, vd)
        return np.array([rc, vc, rd, vd])
    
    def stateTransitionMatrix(self, dt, x):
        """
        State transition function. Assumes point mass gravity

        """
        # State vectors
        rc = x[0:3]
        rd = x[6:9]
        rcMag = la.norm(rc)
        rdMag = la.norm(rd)
        rcHat = rc/rcMag
        rdHat = rd/rdMag 
        # State transition matrix
        FC1 = np.eye(3)
        FC2 = -(orb.MU_EARTH/(rcMag**3))*(np.eye(3)-3*np.matmul(rcHat,np.transpose(rcHat)))
        FD1 = np.eye(3)
        FD2 = -(orb.MU_EARTH/(rdMag**3))*(np.eye(3)-3*np.matmul(rdHat,np.transpose(rdHat)))
        F = np.block([
                [np.zeros((3,3)),FC1,np.zeros((3,3)),np.zeros((3,3))],
                [FC2,np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3))],
                [np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3)),FD1],
                [np.zeros((3,3)),np.zeros((3,3)),FD2,np.zeros((3,3))]])
        return np.eye(12) + F*dt
    
    def processNoise(self, dt):
        ncvQ = ncvProcessNoise(dt)
        Q = self.procVar*np.block([
                [ncvQ,             np.zeros((6, 6))],
                [np.zeros((6, 6)), ncvQ            ]])
        return Q
        

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
            self.P = np.matmul(F,np.matmul(self.P,np.transpose(F))) + self.Q(dt) + self.S
        else:
            self.P = np.matmul(F,np.matmul(self.P,np.transpose(F))) + self.Q(dt)
            
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
        self.P = np.matmul(np.eye(self.n) - np.matmul(W,H),self.P)
        
def residualCov(P, H, R):
    return np.matmul(H,np.matmul(P,np.transpose(H))) + R
        
def kalmanGain(P, H, S):
    return np.matmul(P,np.matmul(np.transpose(H),la.inv(S)))

def normEstErrSqr(dx, P):
    return np.matmul(np.transpose(dx),np.matmul(np.inv(P),dx))

def normInvnSqr(nu, S):
    return np.matmul(np.transpose(nu),np.matmul(np.inv(S),nu))

def ncvProcessNoise(dt):
    return np.block([
                    [np.eye(3)*(dt**3)/3,np.eye(3)*(dt**2)/2],
                    [np.eye(3)*(dt**2)/2,np.eye(3)*dt]])
