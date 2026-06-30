# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 23:18:57 2026

@author: Hakim Lachnani
"""

import numpy as np
from numpy import linalg as la
from dynamics import dynamicsUtils as uDyn
from kinematics import kinematicsUtils as uKin
from dynamics import orbit as orb
from dynamics import ephemerides as eph
import measurements as meas

class DualInertialEKF:
    """
    Dual Inertial EKF is based on:
    Woffinden, David Charles, "Angles-Only Navigation for Autonomous 
    Orbital Rendezvous" (2008). All Graduate Theses and Dissertations. 12.
    https://digitalcommons.usu.edu/etd/12
    """
    
    def __init__(
            self,
            tJ2000, rc, vc, Pc, rd, vd, Pd, procVar, dvVar, pert = None
            ):
        
        # Initialize the inertial nav states nav states
        self.tJ2000 = tJ2000
        self.rsoPosInr = rc
        self.rsoVelInr = vc
        self.rsoCovInr = Pc
        self.svPosInr = rd
        self.svVelInr = vd
        self.svCovInr = Pd
        
        # Initialize DCMs
        self.dcmInr2Ric = np.zeros((3,3))
        self.dcmRic2Los = np.zeros((3,3))
        self.dcmInr2Los = np.zeros((3,3))
        uKin.dcmInr2Ric(self.rsoPosInr, self.rsoVelInr, self.dcmInr2Ric)
        
        # Initialize sun and moon ephemeris
        self.sun = eph.SunEphemeris(self.tJ2000)
        self.moon = eph.MoonEphemeris(self.tJ2000)    
        
        # Initialize the filter
        x = np.concatenate([self.rsoPosInr, self.rsoVelInr, self.svPosInr, self.svVelInr])
        P = np.block([
                [Pc,               np.zeros((6, 6))],
                [np.zeros((6, 6)), Pd              ]])
        S = np.block([
            [np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3))],
            [np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3))],
            [np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3))],
            [np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3)),dvVar*np.eye(3)]])
        def processNoise(dt):
            ncvQ = ncvProcessNoise(dt)
            Q = procVar*np.block([
                    [ncvQ,             np.zeros((6, 6))],
                    [np.zeros((6, 6)), ncvQ            ]])
            return Q
        self.ekf = ExtendedKalmanFilter(
            self.tJ2000, 
            x, 
            P, 
            processNoise, 
            stateUpdateDualInertial, 
            stmDualInertial, 
            S)
        
        # Derived relative states
        self.relPosRectRic = np.zeros((3,))
        self.relVelRectRic = np.zeros((3,))
        uKin.rv2ric(self.rsoPosInr, self.rsoVelInr, self.svPosInr, self.svVelInr, self.relPosRectRic, self.relVelRectRic)
        self.relCovRectRic = covInrToRic(P,self.dcmInr2Ric)
        uKin.dcmRic2Los(self.relPosRectRic, self.dcmRic2Los)
        self.dcmInr2Los = np.matmul(self.dcmRic2Los,self.dcmInr2Ric)
        
        # Compute measurement parameters
        self.az, self.el = meas.calcAzEl(self.rsoPosInr, self.svPosInr, self.dcmInr2Los)
        self.rng = la.norm(self.relPosRectRic)
        self.rngRate = np.dot(self.relPosRectRic, self.relVelRectRic) / self.rng
        
    def propagate(self, dt, u):
        self.ekf.propagate(dt, u)
        
    def update(self, meas, measType):
        # Determine expected measurement
        self.az, self.el = meas.calcAzEl(self.rsoPosInr, self.svPosInr, self.dcmInr2Los)
        self.rng = la.norm(self.relPosRectRic)
        self.rngRate = np.dot(self.relPosRectRic, self.relVelRectRic) / self.rng
        self.measExpected = np.array([self.az, self.el, self.rng, self.rngRate])
        # Residual
        self.meas = meas
        self.measType = measType
        self.measResidual = self.meas - self.measExpected
        # Sensitivity matrix
        self.measSensititivityMat = meas.sensitivityDualInertial(self)
        # Index based on measurement type
        self.measIndx = meas.measType[self.measType]
        # HL TODO: For now, hardcoding R
        R = (30e-6)**2*np.eye(2) # 30 urad^2
        # Call base EKF
        self.ekf.update(
            self.measResidual[self.measIndx], 
            self.measSensititivityMat[self.measIndx,:], 
            R)
        
    def sync(self):
        # Time
        self.tJ2000 = self.ekf.t
        # Absolute states from ekf
        self.rsoPosInr = self.ekf.x[0:3]
        self.rsoVelInr = self.ekf.x[3:6]
        self.rsoCovInr = self.ekf.P[0:6,0:6]
        self.svPosInr = self.ekf.x[6:9]
        self.svVelInr = self.ekf.x[9:12]
        self.svCovInr = self.ekf.P[6:12,6:12]
        # Inertial to RIC DCM
        uKin.dcmInr2Ric(self.rsoPosInr, self.rsoVelInr, self.dcmInr2Ric)
        # Relative states from absolute states
        uKin.rv2ric(self.rsoPosInr, self.rsoVelInr, self.svPosInr, self.svVelInr, self.relPosRectRic, self.relVelRectRic)
        self.relCovRectRic = covInrToRic(self.ekf.P,self.dcmInr2Ric)
        uKin.dcmRic2Los(self.relPosRectRic, self.dcmRic2Los)
        self.dcmInr2Los = np.matmul(self.dcmRic2Los,self.dcmInr2Ric)
        # Compute measurement parameters
        self.az, self.el = meas.calcAzEl(self.rsoPosInr, self.svPosInr, self.dcmInr2Los)
        self.rng = la.norm(self.relPosRectRic)
        self.rngRate = np.dot(self.relPosRectRic, self.relVelRectRic) / self.rng
        
        

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
        2. Woffinden, David Charles, "Angles-Only Navigation for Autonomous 
        Orbital Rendezvous" (2008). All Graduate Theses and Dissertations. 12.
        https://digitalcommons.usu.edu/etd/12
    
    """
    def __init__(
            self, 
            t, x, P, Q, f, F, S, param = None
            ):
        
        # Initialize the filter state
        self.t = t
        self.x = x
        self.P = P
        self.Q = Q
        self.f = f
        self.F = F
        self.S = S
        self.param = param
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
        self.x = self.f(dt, self.x, u, self.param)
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

def stateUpdateInertial(dt, x, u, param = None):
    r = x[0:3]
    v = x[3:6]
    if param == None:
        pert = {
            "jnum": 0,
            "solarGrav": False,
            "lunarGrav": False,
            "SRP": False,
            "drag": False,
            "Cd": 0.0,
            "normalizedArea": 0.0
            }
        
    uDyn.Orbit_rk4(pert["solarGrav"], pert["lunarGrav"], pert["drag"], pert["jnum"], \
                   np.zeros((3,)), np.zeros((3,)), pert["Cd"], pert["normalizedArea"], \
                   u, dt, r, v)
    return np.concatenate([r, v])

def stateUpdateDualInertial(dt, x, u, param = None):
    xc = x[0:6]
    xd = x[6:12]   
    return np.concatenate([stateUpdateInertial(dt, xc, np.zeros((3,)), param), 
                           stateUpdateInertial(dt, xd, u, param)])

def stmInertial(dt, x):
    # State parameters
    r = x[0:3]
    rMag = la.norm(r)
    rHat = r/rMag
    # State transition matrix
    F1 = np.eye(3)
    F2 = -(orb.MU_EARTH/(rMag**3))*(np.eye(3)-3*np.matmul(rHat,np.transpose(rHat)))
    F = np.block([
            [np.zeros((3,3)),F1],
            [F2,np.zeros((3,3))]])
    return np.eye(6) + F*dt

def stmDualInertial(dt, x):
    xc = x[0:6]
    xd = x[6:12]
    return np.block([
            [stmInertial(dt,xc), np.zeros((6, 6))  ],
            [np.zeros((6, 6))  , stmInertial(dt,xd)]])

def covInrToRic(Pi,RN):
    """
    Converts dual inertial convariance to relative RIC covariance.
    Per Eq. 9.11, 9.12 of Woffinden.

    Parameters
    ----------
    Pi : 12x12 double
        Dual inertial covariance.
    RN : 3x3 double
        Inertial to RIC DCM.

    Returns
    -------
    Pric: 6x6 double
        Relative Covariance in RIC.

    """
    Hr = np.block([-np.eye(3),np.zeros((3,3)),np.eye(3),np.zeros((3,3))])
    Hv = np.block([np.zeros((3,3)),-np.eye(3),np.zeros((3,3)),np.eye(3)])
    PrRic = np.matmul(RN,
                      np.matmul(Hr,
                                np.matmul(Pi,
                                          np.matmul(np.transpose(Hr),
                                                    np.transpose(RN)))))
    PvRic = np.matmul(RN,
                      np.matmul(Hv,
                                np.matmul(Pi,
                                          np.matmul(np.transpose(Hv),
                                                    np.transpose(RN)))))
    return np.block([
                    [np.zeros((3,3)),PrRic],
                    [PvRic,np.zeros((3,3))]])

