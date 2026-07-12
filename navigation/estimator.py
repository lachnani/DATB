# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 23:18:57 2026

@author: Hakim Lachnani
"""

import numpy as np
from numpy import linalg as la
from dynamics import dynamicsUtils as uDyn
from kinematics import kinematicsUtils as uKin
from dynamics import ephemerides as eph
import measurements

class RelativeEKF:
    """
    Relative EKF is based on:
    Hablani, "Autonomous Inertial Relative Navigation with 
    Sight-Line-Stabilized Sensors for Spacecraft Rendezvous," 2009
    """
    
    def __init__(
            self,
            tJ2000, rc, vc, Pc, rd, vd, Pd, 
            procVar, dvVar, measCov, pert = None,
            anchor = "DEPUTY"
            ):
        
        # Save settings
        self.anchor = anchor
        
        # Initialize the inertial nav states nav states
        self.tJ2000 = tJ2000
        self.chiefPosInr = rc
        self.chiefVelInr = vc
        self.chiefCovInr = Pc
        self.deputyPosInr = rd
        self.deputyVelInr = vd
        self.deputyCovInr = Pd
        self.crossCovInr  = np.zeros((6,6))
        
        # Initialize DCMs
        self.dcmInr2Ric = np.zeros((3,3))
        self.dcmRic2Los = np.zeros((3,3))
        self.dcmInr2Los = np.zeros((3,3))
        uKin.dcmInr2Ric(self.chiefPosInr, self.chiefVelInr, self.dcmInr2Ric)
        self.omegaRicWrtInrInInr = np.cross(self.chiefPosInr, self.chiefVelInr) / np.dot(self.chiefPosInr,self.chiefPosInr)
        
        # Initialize sun and moon ephemeris
        self.sun = eph.SunEphemeris(self.tJ2000)
        self.moon = eph.MoonEphemeris(self.tJ2000)    
        
        # Relative inertial states
        self.relPosInr = self.deputyPosInr - self.chiefPosInr
        self.relVelInr = self.deputyVelInr - self.chiefVelInr
        self.relCovInr = absCovToRelCov(self.chiefCovInr, self.deputyCovInr, self.crossCovInr)
        
        # Relative RIC states
        self.relPosRectRic = np.zeros((3,))
        self.relVelRectRic = np.zeros((3,))
        uKin.rv2ric(self.chiefPosInr, self.chiefVelInr, self.deputyPosInr, self.deputyVelInr, self.relPosRectRic, self.relVelRectRic)
        self.relCovRectRic = rotateCov(self.relCovInr, self.dcmInr2Ric, self.omegaRicWrtInrInInr)
        uKin.dcmRic2Los(self.relPosRectRic, self.dcmRic2Los)
        self.dcmInr2Los = np.matmul(self.dcmRic2Los,self.dcmInr2Ric)
        
        # Compute measurement parameters
        self.az, self.el = measurements.calcAzEl(self.chiefPosInr, self.deputyPosInr, self.dcmInr2Los)
        self.rng = la.norm(self.relPosRectRic)
        self.rngRate = np.dot(self.relPosRectRic, self.relVelRectRic) / self.rng
        
        # Initialize the filter
        if self.anchor == "CHIEF":
            x = np.concatenate([self.chiefPosInr, self.chiefVelInr, self.relPosInr, self.relVelInr])
            P = np.block([
                    [self.chiefCovInr, np.zeros((6, 6))],
                    [np.zeros((6, 6)), self.relCovInr]])
            S = np.block([
                [np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3))],
                [np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3))],
                [np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3))],
                [np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3)),dvVar*np.eye(3)]])
            stateUpdate = stateUpdateChiefAnchor
        elif self.anchor == "DEPUTY":
            x = np.concatenate([self.deputyPosInr, self.deputyVelInr, self.relPosInr, self.relVelInr])
            P = np.block([
                    [self.deputyCovInr, np.zeros((6, 6))],
                    [np.zeros((6, 6)), self.relCovInr]])
            S = np.block([
                [np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3))],
                [np.zeros((3,3)),dvVar*np.eye(3),np.zeros((3,3)),np.zeros((3,3))],
                [np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3))],
                [np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3)),dvVar*np.eye(3)]])
            stateUpdate = stateUpdateDeputyAnchor
        # HL TODO: Is the relative process noise the same as the absolute?
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
            stateUpdate, 
            stmRelative, 
            S)
        
        
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
        self.measSensititivityMat = measurements.sensitivityRelative(self)
        # Index based on measurement type
        self.measIndx = measurements.measType[self.measType]
        # Call relative EKF
        self.ekf.update(
            self.measResidual[self.measIndx], 
            self.measSensititivityMat[self.measIndx,:], 
            self.R[self.measIndx,self.measIndx])
        
    def sync(self):
        # Time
        self.tJ2000 = self.ekf.t
        # Ephemeris
        self.sun.update(self.tJ2000)
        self.moon.update(self.tJ2000)
        # Relative inertial states from relative EKF
        self.relPosInr = self.ekf.x[6:9]
        self.relVelInr = self.ekf.x[9:12]
        self.relCovInr = self.ekf.P[6:12,6:12]
        # Re-assert decoupling
        self.crossCovInr  = np.zeros((6,6))
        # Absolute states dependent on anchor choice
        if self.anchor == "CHIEF":
            # Chief states from inertial EKF
            self.chiefPosInr = self.ekf.x[0:3]
            self.chiefVelInr = self.ekf.x[3:6]
            self.chiefCovInr = self.ekf.P[0:6,0:6]
            # Deputy states as derived from chief and relative states
            self.deputyPosInr = self.chiefPosInr + self.relPosInr
            self.deputyVelInr = self.chiefVelInr + self.relVelInr
            self.deputyCovInr = self.chiefCovInr + self.relCovInr # Assuming no cross-correlation
        elif self.anchor == "DEPUTY":
            self.deputyPosInr = self.ekf.x[0:3]
            self.deputyVelInr = self.ekf.x[3:6]
            self.deputyCovInr = self.ekf.P[0:6,0:6]
            # Chief states as derived from deputy and relative states
            self.chiefPosInr = self.deputyPosInr - self.relPosInr
            self.chiefVelInr = self.deputyVelInr - self.relVelInr
            self.chiefCovInr = self.deputyCovInr + self.relCovInr # Assuming no cross-correlation
        # Inertial to RIC DCM
        uKin.dcmInr2Ric(self.chiefPosInr, self.chiefVelInr, self.dcmInr2Ric)
        self.omegaRicWrtInrInInr = np.cross(self.chiefPosInr, self.chiefVelInr) / np.dot(self.chiefPosInr,self.chiefPosInr)
        # Relative RIC states
        uKin.rv2ric(self.chiefPosInr, self.chiefVelInr, self.deputyPosInr, self.deputyVelInr, self.relPosRectRic, self.relVelRectRic)
        self.relCovRectRic = rotateCov(self.relCovInr, self.dcmInr2Ric, self.omegaRicWrtInrInInr)
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
        self.crossCovInr  = np.zeros((6,6))
        
        # Initialize DCMs
        self.dcmInr2Ric = np.zeros((3,3))
        self.dcmRic2Los = np.zeros((3,3))
        self.dcmInr2Los = np.zeros((3,3))
        uKin.dcmInr2Ric(self.chiefPosInr, self.chiefVelInr, self.dcmInr2Ric)
        self.omegaRicWrtInrInInr = np.cross(self.chiefPosInr, self.chiefVelInr) / np.dot(self.chiefPosInr,self.chiefPosInr)
        
        # Initialize sun and moon ephemeris
        self.sun = eph.SunEphemeris(self.tJ2000)
        self.moon = eph.MoonEphemeris(self.tJ2000)    
        
        # Relative inertial states
        self.relPosInr = self.deputyPosInr - self.chiefPosInr
        self.relVelInr = self.deputyVelInr - self.chiefVelInr
        self.relCovInr = absCovToRelCov(self.chiefCovInr, self.deputyCovInr, self.crossCovInr)
        
        # Relative RIC states
        self.relPosRectRic = np.zeros((3,))
        self.relVelRectRic = np.zeros((3,))
        uKin.rv2ric(self.chiefPosInr, self.chiefVelInr, self.deputyPosInr, self.deputyVelInr, self.relPosRectRic, self.relVelRectRic)
        self.relCovRectRic = rotateCov(self.relCovInr, self.dcmInr2Ric, self.omegaRicWrtInrInInr)
        uKin.dcmRic2Los(self.relPosRectRic, self.dcmRic2Los)
        self.dcmInr2Los = np.matmul(self.dcmRic2Los,self.dcmInr2Ric)
        
        # Compute measurement parameters
        self.az, self.el = measurements.calcAzEl(self.chiefPosInr, self.deputyPosInr, self.dcmInr2Los)
        self.rng = la.norm(self.relPosRectRic)
        self.rngRate = np.dot(self.relPosRectRic, self.relVelRectRic) / self.rng
        
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
        # Ephemeris
        self.sun.update(self.tJ2000)
        self.moon.update(self.tJ2000)
        # Absolute states from ekf
        self.chiefPosInr = self.ekf.x[0:3]
        self.chiefVelInr = self.ekf.x[3:6]
        self.chiefCovInr = self.ekf.P[0:6,0:6]
        self.deputyPosInr = self.ekf.x[6:9]
        self.deputyVelInr = self.ekf.x[9:12]
        self.deputyCovInr = self.ekf.P[6:12,6:12]
        self.crossCovInr  = self.ekf.P[0:6,6:12]
        # Inertial to RIC DCM
        uKin.dcmInr2Ric(self.chiefPosInr, self.chiefVelInr, self.dcmInr2Ric)
        self.omegaRicWrtInrInInr = np.cross(self.chiefPosInr, self.chiefVelInr) / np.dot(self.chiefPosInr,self.chiefPosInr)
        # Relative inertial states
        self.relPosInr = self.deputyPosInr - self.chiefPosInr
        self.relVelInr = self.deputyVelInr - self.chiefVelInr
        self.relCovInr = absCovToRelCov(self.chiefCovInr, self.deputyCovInr, self.crossCovInr)
        # Relative RIC states
        uKin.rv2ric(self.chiefPosInr, self.chiefVelInr, self.deputyPosInr, self.deputyVelInr, self.relPosRectRic, self.relVelRectRic)
        self.relCovRectRic = rotateCov(self.relCovInr, self.dcmInr2Ric, self.omegaRicWrtInrInInr)
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
            "jnum": 2,
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

def stateUpdateChiefAnchor(dt, x, u, param = None):
    # Extract chief and deputy state
    xc0 = x[0:6]
    xd0 = x[6:12] + xc0
    # Propagate 
    xc = stateUpdateInertial(dt, xc0, np.zeros((3,)), param)
    xd = stateUpdateInertial(dt, xd0, u, param)
    # Reform relative state
    return np.concatenate([xc,xd-xc])

def stateUpdateDeputyAnchor(dt, x, u, param = None):
    # Extract chief and deputy state
    xd0 = x[0:6]
    xc0 = xd0 - x[6:12]
    # Propagate 
    xc = stateUpdateInertial(dt, xc0, np.zeros((3,)), param)
    xd = stateUpdateInertial(dt, xd0, u, param)
    # Reform relative state
    return np.concatenate([xd,xd-xc])

def stmInertial(dt, x):
    """
    Inertial State Transition Matrix. Assumes Earth gravity plus J2. Uses first
    order Taylor series expansion, so it is only valid for small time steps.
    
    Refs: 
        1. Markley, "Approximate Cartesian State Transition Matrix"
        2. Hablani, "Autonomous Inertial Relative Navigation with 
        Sight-Line-Stabilized Sensors for Spacecraft Rendezvous"

    Parameters
    ----------
    dt : double
        propagation delta-time.
    x : 6x1 double
        Inertial position/velocity state.

    Returns
    -------
    F : 
        state transition matrix.

    """
    # State transition matrix
    F2 = np.zeros((3,3))
    uDyn.gravPartial(x[0:3],F2)
    F = np.block([
            [np.zeros((3,3)),np.eye(3)],
            [F2,np.zeros((3,3))]])
    return np.eye(6) + F*dt

def stmDualInertial(dt, x):
    xc = x[0:6]
    xd = x[6:12]
    return np.block([
            [stmInertial(dt,xc), np.zeros((6, 6))  ],
            [np.zeros((6, 6))  , stmInertial(dt,xd)]])

def stmRelative(dt, x):
    # Linearize the inertial motion about the anchor
    stmI = stmInertial(dt,x[0:6])
    return np.block([
            [stmI,               np.zeros((6, 6))  ],
            [np.zeros((6, 6))  , stmI              ]])

def absCovToRelCov(Pc,Pd,Px):
    """
    Converts absolute covariances to relative covariances. Does not rotate 
    frames.
    
    HL TODO: Confirm this math is correct!

    Parameters
    ----------
    Pc : 6x6 double
        Chief inertial covariance.
    Pd : 6x6 double
        Deputy inertial covariance.
    Px : 6x6 double
        Chief/deputy cross covariance.

    Returns
    -------
    Prel : 6x6 double
        Chief to deputy inertial relative covariance.

    """
    # Initialize
    Prel = np.zeros((6,6))
    
    # Construct block covariance
    Prel[0:3,0:3] = Pc[0:3,0:3] + Pd[0:3,0:3] - np.transpose(Px[0:3,0:3]) - Px[0:3,0:3]
    Prel[3:6,3:6] = Pc[3:6,3:6] + Pd[3:6,3:6] - np.transpose(Px[3:6,3:6]) - Px[3:6,3:6]
    Prel[0:3,3:6] = Pc[0:3,3:6] + Pd[0:3,3:6] - np.transpose(Px[3:6,0:3]) - Px[0:3,3:6]
    Prel[3:6,0:3] = np.transpose(Prel[0:3,3:6])
    
    return Prel

def rotateCov(Pa,BA,omegaBwrtAinA):
    """
    Rotates 6-DOF position velocity covariance from frame A to B.
    Ref: Drotar, "Transformation of Covariance Matrices Between Inertial and 
    Earth-Fixed Coordinates"

    Parameters
    ----------
    Pa : 6x6 double
        Covariance in frame A.
    BA : 3x3 double
        Rotation matrix from A to B.
    omegaBwrtAinA : 3x1 double
        angular velocity of frame B with respect to A, expressed in A.

    Returns
    -------
    Pa : 6x6 double
        Covariance in frame B.

    """
    
    # Skew-symmentric operator from omaga
    wx = np.cross(np.eye(3),omegaBwrtAinA)
    
    # Construct 6x6 Jacobian matrix
    J = np.block([[BA,np.zeros((3,3))],
                  [-np.matmul(BA,wx),BA]])
    
    # Rotate covariance
    return np.matmul(J, np.matmul(Pa, np.transpose(J)))

def initCovFromRic(varRic,dcmInr2Ric,omegaRicWrtInrInInr):
    """
    Initializes inertial covariance from RIC uncertainties. Assumes a radial-
    in-track correlation coefficient derived from HCW dynamics.
    
    Ref: Woffinden, David Charles, "Angles-Only Navigation for Autonomous 
    Orbital Rendezvous" (2008). All Graduate Theses and Dissertations. 12.
    https://digitalcommons.usu.edu/etd/12 

    Parameters
    ----------
    varRic : 6x1 double 
        RIC frame position and velocity variances.
    dcmInr2Ric : 3x3 double
        Inertial to RIC DCM.
    omegaRicWrtInrInInr : 3x1 double
        Angular velocity of RIC frame w.r.t Inertial frame.

    Returns
    -------
    PInr : 6x6 double
        Initial covariance in the inertial frame.

    """
    # Correlation coefficient
    f = -np.sqrt(3)/2
    
    # Construct RIC STM (Eq 6.25)
    PRic = varRic*np.eye(6)
    PRic[0,4] = f*np.sqrt(varRic[0])*np.sqrt(varRic[4])
    PRic[4,0] = PRic[0,4]
    PRic[1,3] = f*np.sqrt(varRic[1])*np.sqrt(varRic[3])
    PRic[3,1] = PRic[1,3]
    
    # Rotate to Inertial frame
    dcmRic2Inr = np.transpose(dcmInr2Ric)
    omegaInrWrtRicInRic = -np.matmul(dcmInr2Ric,omegaRicWrtInrInInr)
    return rotateCov(PRic, dcmRic2Inr, omegaInrWrtRicInRic)