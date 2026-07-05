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
import measurements

class RelativeEKF:
    """
    Relative EKF is a two parrallel EKFs with an absolute state and a relative state
    """
    
    def __init__(
            self,
            tJ2000, rc, vc, Pc, rd, vd, Pd, 
            procVarInr, procVarRel, dvVar, measCov, pert = None,
            anchor = "DEPUTY", frame = "RECTRIC", stm = "HCW"
            ):
        
        # Save settings
        self.anchor = anchor
        self.frame = frame
        self.stm = stm
        
        # Initialize the inertial nav states nav states
        self.tJ2000 = tJ2000
        self.chiefPosInr = rc
        self.chiefVelInr = vc
        self.chiefCovInr = Pc
        self.deputyPosInr = rd
        self.deputyVelInr = vd
        self.deputyCovInr = Pd
        
        # Initialize DCMs
        self.dcmInr2Ric = np.zeros((3,3))
        self.dcmRic2Los = np.zeros((3,3))
        self.dcmInr2Los = np.zeros((3,3))
        uKin.dcmInr2Ric(self.chiefPosInr, self.chiefVelInr, self.dcmInr2Ric)
        
        # Initialize sun and moon ephemeris
        self.sun = eph.SunEphemeris(self.tJ2000)
        self.moon = eph.MoonEphemeris(self.tJ2000)  
        
        # Common filter paramters
        S = np.block([
            [np.zeros((3,3)),np.zeros((3,3))],
            [np.zeros((3,3)),dvVar*np.eye(3)]])
        self.R = measCov
        
        # Initialize the inertial filter
        if self.anchor == "CHIEF":
            xInr = np.concatenate([self.chiefPosInr, self.chiefVelInr])
            PInr = Pc
        elif self.anchor == "DEPUTY":
            xInr = np.concatenate([self.deputyPosInr, self.deputyVelInr])
            PInr = Pd
        def processNoiseInr(dt):
            ncvQ = ncvProcessNoise(dt)
            Q = procVarInr*ncvQ
            return Q       
        self.ekfInr = ExtendedKalmanFilter(
            self.tJ2000, 
            xInr, 
            PInr, 
            processNoiseInr, 
            stateUpdateInertial, 
            stmInertial, 
            S)
        
        # Derived relative states
        self.relPosRectRic = np.zeros((3,))
        self.relVelRectRic = np.zeros((3,))
        uKin.rv2ric(self.chiefPosInr, self.chiefVelInr, self.deputyPosInr, self.deputyVelInr, self.relPosRectRic, self.relVelRectRic)
        self.relCovRectRic = covInrToRic(np.block([
                [Pc,               np.zeros((6, 6))],
                [np.zeros((6, 6)), Pd              ]]),self.dcmInr2Ric)
        uKin.dcmRic2Los(self.relPosRectRic, self.dcmRic2Los)
        self.dcmInr2Los = np.matmul(self.dcmRic2Los,self.dcmInr2Ric)
        
        # Initialize the relative filter
        xRel = np.concatenate([self.relPosRectRic, self.relVelRectRic])
        def processNoiseRel(dt):
            ncvQ = ncvProcessNoise(dt)
            Q = procVarRel*ncvQ
            return Q    
        # HL TODO: Replace inertial state update and stm with relative versions
        self.ekfRel = ExtendedKalmanFilter(
            self.tJ2000, 
            xRel, 
            self.relCovRectRic, 
            processNoiseRel, 
            stateUpdateInertial, 
            stmInertial, 
            S)
        
        # Compute measurement parameters
        self.az, self.el = measurements.calcAzEl(self.chiefPosInr, self.deputyPosInr, self.dcmInr2Los)
        self.rng = la.norm(self.relPosRectRic)
        self.rngRate = np.dot(self.relPosRectRic, self.relVelRectRic) / self.rng
        
    def propagate(self, dt, aCtrlInEci):
        if self.anchor == "CHIEF":
            self.ekfInr.propagate(dt, np.zeros((3,)))
        elif self.anchor == "DEPUTY":
            self.ekfInr.propagate(dt, aCtrlInEci)
        self.ekfRel.propagate(dt, np.matmul(self.dcmInr2Ric,aCtrlInEci))
        
    def update(self, meas, measType):
        # Determine expected measurement
        self.az, self.el = measurements.calcAzEl(self.chiefPosInr, self.deputyPosInr, self.dcmInr2Los)
        self.rng = la.norm(self.relPosRectRic)
        self.rngRate = np.dot(self.relPosRectRic, self.relVelRectRic) / self.rng
        self.measExpected = np.array([self.az, self.el, self.rng, self.rngRate])
        # Residual
        self.meas = meas
        self.measType = measType
        self.measResidual = self.meas - self.measExpected
        # Sensitivity matrix
        self.measSensititivityMat = measurements.sensitivityRelative(self)
        # Index based on measurement type
        self.measIndx = measurements.measType[self.measType]
        # Call relative EKF
        self.ekfRel.update(
            self.measResidual[self.measIndx], 
            self.measSensititivityMat[self.measIndx,:], 
            self.R[self.measIndx,self.measIndx])
        
    def sync(self):
        # Time
        self.tJ2000 = self.ekf.t
        # Relative states from relative EKF
        self.relPosRectRic = self.ekfRel.x[0:3]
        self.relVelRectRic = self.ekfRel.x[3:6]
        self.relCovRectRic = self.ekfRel.P
        # Absolute states dependent on anchor choice
        if self.anchor == "CHIEF":
            # Chief states from inertial EKF
            self.chiefPosInr = self.ekfInr.x[0:3]
            self.chiefVelInr = self.ekfInr.x[3:6]
            self.chiefCovInr = self.ekfInr.P
            # Deputy states as derived from chief and relative states
            uKin.ric2rv(self.chiefPosInr, self.chiefVelInr, self.relPosRectRic, self.relVelRectRic, self.deputyPosInr, self.deputyVelInr)
        elif self.anchor == "DEPUTY":
            self.deputyPosInr = self.ekfInr.x[0:3]
            self.deputyVelInr = self.ekfInr.x[3:6]
            self.deputyCovInr = self.ekfInr.P
            # Chief states as derived from deputy and relative states
            uKin.ric2rv(self.deputyPosInr, self.deputyVelInr, -1*self.relPosRectRic, -1*self.relVelRectRic, self.chiefPosInr, self.chiefVelInr)
        # Inertial to RIC DCM
        uKin.dcmInr2Ric(self.chiefPosInr, self.chiefVelInr, self.dcmInr2Ric)
        # Covariance of the remaining state
        if self.anchor == "CHIEF":
            self.deputyCovInr = covRicToInr(self.chiefCovInr, self.relCovRectRic, self.dcmInr2Ric)
        elif self.anchor == "DEPUTY":
            self.chiefCovInr = covRicToInr(self.deputyCovInr, self.relCovRectRic, self.dcmInr2Ric)
        # LOS DCMs
        uKin.dcmRic2Los(self.relPosRectRic, self.dcmRic2Los)
        self.dcmInr2Los = np.matmul(self.dcmRic2Los,self.dcmInr2Ric)
        # Compute measurement parameters
        self.az, self.el = measurements.calcAzEl(self.chiefPosInr, self.deputyPosInr, self.dcmInr2Los)
        self.rng = la.norm(self.relPosRectRic)
        self.rngRate = np.dot(self.relPosRectRic, self.relVelRectRic) / self.rng
        

class DualInertialEKF:
    """
    Dual Inertial EKF is based on:
    Woffinden, David Charles, "Angles-Only Navigation for Autonomous 
    Orbital Rendezvous" (2008). All Graduate Theses and Dissertations. 12.
    https://digitalcommons.usu.edu/etd/12
    """
    
    def __init__(
            self,
            tJ2000, rc, vc, Pc, rd, vd, Pd, procVar, dvVar, measCov, pert = None
            ):
        
        # Initialize the inertial nav states nav states
        self.tJ2000 = tJ2000
        self.chiefPosInr = rc
        self.chiefVelInr = vc
        self.chiefCovInr = Pc
        self.deputyPosInr = rd
        self.deputyVelInr = vd
        self.deputyCovInr = Pd
        
        # Initialize DCMs
        self.dcmInr2Ric = np.zeros((3,3))
        self.dcmRic2Los = np.zeros((3,3))
        self.dcmInr2Los = np.zeros((3,3))
        uKin.dcmInr2Ric(self.chiefPosInr, self.chiefVelInr, self.dcmInr2Ric)
        
        # Initialize sun and moon ephemeris
        self.sun = eph.SunEphemeris(self.tJ2000)
        self.moon = eph.MoonEphemeris(self.tJ2000)    
        
        # Initialize the filter
        x = np.concatenate([self.chiefPosInr, self.chiefVelInr, self.deputyPosInr, self.deputyVelInr])
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
        self.R = measCov
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
        uKin.rv2ric(self.chiefPosInr, self.chiefVelInr, self.deputyPosInr, self.deputyVelInr, self.relPosRectRic, self.relVelRectRic)
        self.relCovRectRic = covInrToRic(P,self.dcmInr2Ric)
        uKin.dcmRic2Los(self.relPosRectRic, self.dcmRic2Los)
        self.dcmInr2Los = np.matmul(self.dcmRic2Los,self.dcmInr2Ric)
        
        # Compute measurement parameters
        self.az, self.el = measurements.calcAzEl(self.chiefPosInr, self.deputyPosInr, self.dcmInr2Los)
        self.rng = la.norm(self.relPosRectRic)
        self.rngRate = np.dot(self.relPosRectRic, self.relVelRectRic) / self.rng
        
    def propagate(self, dt, aCtrlInEci):
        self.ekf.propagate(dt, aCtrlInEci)
        
    def update(self, meas, measType):
        # Determine expected measurement
        self.az, self.el = measurements.calcAzEl(self.chiefPosInr, self.deputyPosInr, self.dcmInr2Los)
        self.rng = la.norm(self.relPosRectRic)
        self.rngRate = np.dot(self.relPosRectRic, self.relVelRectRic) / self.rng
        self.measExpected = np.array([self.az, self.el, self.rng, self.rngRate])
        # Residual
        self.meas = meas
        self.measType = measType
        self.measResidual = self.meas - self.measExpected
        # Sensitivity matrix
        self.measSensititivityMat = measurements.sensitivityDualInertial(self)
        # Index based on measurement type
        self.measIndx = measurements.measType[self.measType]
        # Call base EKF
        self.ekf.update(
            self.measResidual[self.measIndx], 
            self.measSensititivityMat[self.measIndx,:], 
            self.R[self.measIndx,self.measIndx])
        
    def sync(self):
        # Time
        self.tJ2000 = self.ekf.t
        # Absolute states from ekf
        self.chiefPosInr = self.ekf.x[0:3]
        self.chiefVelInr = self.ekf.x[3:6]
        self.chiefCovInr = self.ekf.P[0:6,0:6]
        self.deputyPosInr = self.ekf.x[6:9]
        self.deputyVelInr = self.ekf.x[9:12]
        self.deputyCovInr = self.ekf.P[6:12,6:12]
        # Inertial to RIC DCM
        uKin.dcmInr2Ric(self.chiefPosInr, self.chiefVelInr, self.dcmInr2Ric)
        # Relative states from absolute states
        uKin.rv2ric(self.chiefPosInr, self.chiefVelInr, self.deputyPosInr, self.deputyVelInr, self.relPosRectRic, self.relVelRectRic)
        self.relCovRectRic = covInrToRic(self.ekf.P,self.dcmInr2Ric)
        uKin.dcmRic2Los(self.relPosRectRic, self.dcmRic2Los)
        self.dcmInr2Los = np.matmul(self.dcmRic2Los,self.dcmInr2Ric)
        # Compute measurement parameters
        self.az, self.el = measurements.calcAzEl(self.chiefPosInr, self.deputyPosInr, self.dcmInr2Los)
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
                    [PrRic,np.zeros((3,3))],
                    [np.zeros((3,3)),PvRic]])

def covRicToInr(Pi,Pric,RN):
    """
    Converts relative RIC covariance to dual inertial convariance.

    Parameters
    ----------
    Pi : 6x6 double
        Inertial covariance for either chief or deputy in Inertial frame.
    Pric : 6x6
        Relative covariance in the RIC frame.
    RN : 3x3 double
        Inertial to RIC DCM.

    Returns
    -------
    Pi2: 6x6 double
        Inertial covariance for either chief or deputy in Inertial frame.

    """
    # Second inertial covariance is relative cov minus first inertial cov
    Pi2 = np.matmul(np.transpose(RN),np.matmul(Pric,RN)) - Pi
    # HL TODO: Force positive definite
    # Pi2[Pi2 < 0] = 0
    return Pi2