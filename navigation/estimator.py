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
import kalmanFilter as kf


class DualInertialFilter:
    """
    Dual Inertial Filter class for translational spacecraft rendezvous. Primary 
    filter states are:
        deputyPosInr
        deputyVelInr
        chiefPosInr
        deputyVelInr
    Class also populates remaining navigation states of interest
    
    Ref: Hakim Lachnani and Kevin Schroeder, "Comparative Analysis of 
    Navigation Filter Formulations for Spacecraft Rendezvous"
    
    Parameters
    ----------
    tJ2000 : double
        time since J2000 epoch.
    rc : 3x1 double
        Chief inertial position.
    vc : 3x1 double
        Chief inertial velocity.
    Pc : 6x6 double
        Chief inertial covariance.
    rd : 3x1 double
        Deputy inertial position.
    vd : 3x1 double
        Deputy inertial velocity.
    Pd : 6x6 double
        Deputy inertial covariance.
    Qc : 3x3 double
        Chief process noise power spectral density in RIC.
    Qd : 3x3 double
        Deputy process noise power spectral density in RIC.
    dvVar : 3x1 double
        Delta-V process noise array (scale factor, quantization, pointing).
    measCov : 4x4 double
        Measurement covariance matrix (az, el, rng, rngRate).
    pertc: dictionary
        Chief perturbation dictionary.
    pertd: dictionary
        Deputy perturbation dictionary.
    coupling: boolean
        Filter coupling.
        
    """
    
    def __init__(
            self,
            tJ2000, 
            rc, vc, Pc, 
            rd, vd, Pd, 
            Qc, Qd, dvVar, measCov, 
            pertc = None, pertd = None,
            coupling = True
            ):
        
        # Filter constants
        self.numStates = 12
        self.coupling = coupling
        
        # Initialize the inertial nav states nav states
        self.tJ2000 = tJ2000
        self.chiefPosInr = rc
        self.chiefVelInr = vc
        self.chiefCovInr = Pc
        self.deputyPosInr = rd
        self.deputyVelInr = vd
        self.deputyCovInr = Pd
        self.deputyChiefCrossCovInr  = np.zeros((6,6))
        
        # Initialize DCMs
        self.dcmInr2Ric = np.zeros((3,3))
        self.dcmInr2DepRic = np.zeros((3,3))
        self.dcmRic2Los = np.zeros((3,3))
        self.dcmInr2Los = np.zeros((3,3))
        uKin.dcmInr2Ric(self.chiefPosInr, self.chiefVelInr, self.dcmInr2Ric)
        uKin.dcmInr2Ric(self.deputyPosInr, self.deputyVelInr, self.dcmInr2DepRic)
        self.omegaRicWrtInrInInr = np.cross(self.chiefPosInr, self.chiefVelInr) / np.dot(self.chiefPosInr,self.chiefPosInr)
        
        # Initialize sun and moon ephemeris
        self.sun = eph.SunEphemeris(self.tJ2000)
        self.moon = eph.MoonEphemeris(self.tJ2000)    
        
        # Relative inertial states
        self.relPosInr = self.deputyPosInr - self.chiefPosInr
        self.relVelInr = self.deputyVelInr - self.chiefVelInr
        self.relCovInr = self.deputyCovInr + self.chiefCovInr - self.deputyChiefCrossCovInr - np.transpose(self.deputyChiefCrossCovInr)
        
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
        self.measCov = measCov
        
        # Save process noise matrices
        self.deputyProcNoiseInRic = Qd
        self.chiefProcNoiseInRic = Qc
        self.dvProcNoise = dvVar
        
        # Save perturbation libraries
        self.chiefPerturbations = pertc
        self.deputyPerturbations = pertd
        
        # Initialize the filter states
        self.x = np.concatenate([self.deputyPosInr, self.deputyVelInr, self.chiefPosInr, self.chiefVelInr])
        self.P = np.block([
                [Pd,               np.zeros((6, 6))],
                [np.zeros((6, 6)), Pc              ]])

        
    def propagate(self, dt, aCtrlInEci):        
        # Compute process nosie
        deputyProcNoiseInr = ncvProcessNoise(dt, np.matmul(np.transpose(self.dcmInr2DepRic),self.deputyProcNoiseInRic))
        chiefProcNoiseInr = ncvProcessNoise(dt, np.matmul(np.transpose(self.dcmInr2Ric),self.chiefProcNoiseInRic))
        stateProcNoiseInr = np.block([
                                     [deputyProcNoiseInr, np.zeros((6, 6))],
                                     [np.zeros((6, 6))  , chiefProcNoiseInr]])
        if la.norm(aCtrlInEci) > 0.0:
            # Add maneuver noise
            maneuverProcNoise = dvProcessNoise(dt*aCtrlInEci,self.dvVar[0],self.dvVar[1],self.dvVar[2])
            maneuverProcNoiseInr = np.block([
                                            [maneuverProcNoise, np.zeros((6, 6))],
                                            [np.zeros((6, 6)) , np.zeros((6, 6))]])
            stateProcNoiseInr = stateProcNoiseInr + maneuverProcNoiseInr
            
        # Compute state transition matrix
        stm = np.block([
                       [stmInertial(dt,self.x[0:6]), np.zeros((6, 6))            ],
                       [np.zeros((6, 6))           , stmInertial(dt,self.x[6:12])]])
        
        # Update timestep
        self.tJ2000 = self.tJ2000 + dt
        
        # Update ephemerides
        self.sun.update(self.tJ2000)
        self.moon.update(self.tJ2000)
        
        # Propagate States
        self.x[0:6] = stateUpdateInertial(dt, self.x[0:6], aCtrlInEci, self.deputyPerturbations)
        self.x[6:12] = stateUpdateInertial(dt, self.x[6:12], np.zeros((3,)), self.chiefPerturbations)
        
        # Propagate Covariance
        self.P = kf.propagateCov(self.P, stm, stateProcNoiseInr)
        
        # Sync filter to update all intermediate states
        self.sync()
        
        
    def update(self, meas, measType):
        # Determine expected measurement
        self.measExpected = np.array([self.az, self.el, self.rng, self.rngRate])
        
        # Compute residual
        self.meas = meas
        self.measType = measType
        self.measResidual = self.meas - self.measExpected
        
        # Compute sensitivity matrix
        H = measurements.inertialMeasurementSensitivity(self)
        if self.coupling:
            # Coupled filter: measurements affect both states
            self.measSensititivityMat = np.block([-H,H])
        else:
            # Decoupled filter: measurements affect chief state only
            self.measSensititivityMat = np.block([np.zeros((4,6)),H])
        
        # Index based on measurement type
        self.measIndx = measurements.measType[self.measType]
        
        # Perform state update
        self.x, self.P = kf.measurementUpdate(self.x, 
                                              self.P, 
                                              self.numStates, 
                                              self.measResidual[self.measIndx], 
                                              self.measSensititivityMat[self.measIndx,:],
                                              self.measCov[self.measIndx,self.measIndx])
        
        # Sync filter to update all intermediate states
        self.sync()

        
    def sync(self):
        # Deputy and Chief inertial states
        self.deputyPosInr = self.x[0:3]
        self.deputyVelInr = self.x[3:6]
        self.deputyfCovInr = self.P[0:6,0:6]
        self.chiefPosInr = self.x[6:9]
        self.chiefVelInr = self.x[9:12]
        self.chiefCovInr = self.P[6:12,6:12]
        self.deputyChiefCrossCovInr  = self.P[0:6,6:12]
        # Inertial to RIC DCMs
        uKin.dcmInr2Ric(self.chiefPosInr, self.chiefVelInr, self.dcmInr2Ric)
        uKin.dcmInr2Ric(self.deputyPosInr, self.deputyVelInr, self.dcmInr2DepRic)
        self.omegaRicWrtInrInInr = np.cross(self.chiefPosInr, self.chiefVelInr) / np.dot(self.chiefPosInr,self.chiefPosInr)
        # Relative inertial states
        self.relPosInr = self.deputyPosInr - self.chiefPosInr
        self.relVelInr = self.deputyVelInr - self.chiefVelInr
        self.relCovInr = self.deputyCovInr + self.chiefCovInr - self.deputyChiefCrossCovInr - np.transpose(self.deputyChiefCrossCovInr)
        # Relative RIC states
        uKin.rv2ric(self.chiefPosInr, self.chiefVelInr, self.deputyPosInr, self.deputyVelInr, self.relPosRectRic, self.relVelRectRic)
        self.relCovRectRic = rotateCov(self.relCovInr, self.dcmInr2Ric, self.omegaRicWrtInrInInr)
        uKin.dcmRic2Los(self.relPosRectRic, self.dcmRic2Los)
        self.dcmInr2Los = np.matmul(self.dcmRic2Los,self.dcmInr2Ric)
        # Compute measurement parameters
        self.az, self.el = measurements.calcAzEl(self.chiefPosInr, self.deputyPosInr, self.dcmInr2Los)
        self.rng = la.norm(self.relPosRectRic)
        self.rngRate = np.dot(self.relPosRectRic, self.relVelRectRic) / self.rng   
        
class InertialRelativeFilter:
    """
    Inertial Relative Filter class for translational spacecraft rendezvous. 
    Primary filter states are:
        deputyPosInr
        deputyVelInr
        relPosInr
        relVelInr
    Class also populates remaining navigation states of interest
    
    Ref: Hakim Lachnani and Kevin Schroeder, "Comparative Analysis of 
    Navigation Filter Formulations for Spacecraft Rendezvous"
    
    Parameters
    ----------
    tJ2000 : double
        time since J2000 epoch.
    rc : 3x1 double
        Chief inertial position.
    vc : 3x1 double
        Chief inertial velocity.
    Pc : 6x6 double
        Chief inertial covariance.
    rd : 3x1 double
        Deputy inertial position.
    vd : 3x1 double
        Deputy inertial velocity.
    Pd : 6x6 double
        Deputy inertial covariance.
    Qc : 3x3 double
        Chief process noise power spectral density in RIC.
    Qd : 3x3 double
        Deputy process noise power spectral density in RIC.
    dvVar : 3x1 double
        Delta-V process noise array (scale factor, quantization, pointing).
    measCov : 4x4 double
        Measurement covariance matrix (az, el, rng, rngRate).
    pertc: dictionary
        Chief perturbation dictionary.
    pertd: dictionary
        Deputy perturbation dictionary.
    coupling: boolean
        Filter coupling.
        
    """
    
    def __init__(
            self,
            tJ2000, 
            rc, vc, Pc, 
            rd, vd, Pd, 
            Qc, Qd, dvVar, measCov, 
            pertc = None, pertd = None,
            coupling = True
            ):
        
        # Filter constants
        self.numStates = 12
        self.coupling = coupling
        
        # Initialize the inertial nav states nav states
        self.tJ2000 = tJ2000
        self.chiefPosInr = rc
        self.chiefVelInr = vc
        self.chiefCovInr = Pc
        self.deputyPosInr = rd
        self.deputyVelInr = vd
        self.deputyCovInr = Pd
        self.deputyChiefCrossCovInr  = np.zeros((6,6))
        
        # Initialize DCMs
        self.dcmInr2Ric = np.zeros((3,3))
        self.dcmInr2DepRic = np.zeros((3,3))
        self.dcmRic2Los = np.zeros((3,3))
        self.dcmInr2Los = np.zeros((3,3))
        uKin.dcmInr2Ric(self.chiefPosInr, self.chiefVelInr, self.dcmInr2Ric)
        uKin.dcmInr2Ric(self.deputyPosInr, self.deputyVelInr, self.dcmInr2DepRic)
        self.omegaRicWrtInrInInr = np.cross(self.chiefPosInr, self.chiefVelInr) / np.dot(self.chiefPosInr,self.chiefPosInr)
        
        # Initialize sun and moon ephemeris
        self.sun = eph.SunEphemeris(self.tJ2000)
        self.moon = eph.MoonEphemeris(self.tJ2000)    
        
        # Relative inertial states
        self.relPosInr = self.deputyPosInr - self.chiefPosInr
        self.relVelInr = self.deputyVelInr - self.chiefVelInr
        self.relCovInr = self.deputyCovInr + self.chiefCovInr - self.deputyChiefCrossCovInr - np.transpose(self.deputyChiefCrossCovInr)
        
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
        self.measCov = measCov
        
        # Save process noise matrices
        self.deputyProcNoiseInRic = Qd
        self.chiefProcNoiseInRic = Qc
        self.dvProcNoise = dvVar
        
        # Save perturbation libraries
        self.chiefPerturbations = pertc
        self.deputyPerturbations = pertd
        
        # Initialize the filter states
        self.x = np.concatenate([self.deputyPosInr, self.deputyVelInr, self.relPosInr, self.relVelInr])
        if self.coupling:
            self.P = np.block([
                              [Pd, Pd     ],
                              [Pd, Pd + Pc]])
        else:
            self.P = np.block([
                              [Pd,               np.zeros((6,6))],
                              [np.zeros((6,6)) , Pd + Pc        ]])

        
    def propagate(self, dt, aCtrlInEci):        
        # Compute process nosie
        deputyProcNoiseInr = ncvProcessNoise(dt, np.matmul(np.transpose(self.dcmInr2DepRic),self.deputyProcNoiseInRic))
        chiefProcNoiseInr = ncvProcessNoise(dt, np.matmul(np.transpose(self.dcmInr2Ric),self.chiefProcNoiseInRic))
        stateProcNoiseInr = np.block([
                                     [deputyProcNoiseInr, deputyProcNoiseInr                    ],
                                     [deputyProcNoiseInr, deputyProcNoiseInr + chiefProcNoiseInr]])
        if la.norm(aCtrlInEci) > 0.0:
            # Add maneuver noise
            maneuverProcNoise = dvProcessNoise(dt*aCtrlInEci,self.dvVar[0],self.dvVar[1],self.dvVar[2])
            maneuverProcNoiseInr = np.block([
                                            [maneuverProcNoise, np.zeros((6, 6))],
                                            [np.zeros((6, 6)) , np.zeros((6, 6))]])
            stateProcNoiseInr = stateProcNoiseInr + maneuverProcNoiseInr
            
        # Compute state transition matrix
        deputyStmInr = stmInertial(dt,self.x[0:6])
        chiefStmInr = stmInertial(dt,np.concatenate((self.chiefPosInr, self.chiefVelInr), axis=0))
        stm = np.block([
                       [deputyStmInr,               np.zeros((6, 6))],
                       [deputyStmInr - chiefStmInr, chiefStmInr     ]])
        
        # Update timestep
        self.tJ2000 = self.tJ2000 + dt
        
        # Update ephemerides
        self.sun.update(self.tJ2000)
        self.moon.update(self.tJ2000)
        
        # Propagate States
        self.x[0:6] = stateUpdateInertial(dt, self.x[0:6], aCtrlInEci, self.deputyPerturbations)
        self.x[6:12] = self.x[0:6] - stateUpdateInertial(
            dt, np.concatenate((self.chiefPosInr, self.chiefVelInr)), np.zeros((3,)), self.chiefPerturbations)
        
        # Propagate Covariance
        self.P = kf.propagateCov(self.P, stm, stateProcNoiseInr)
        
        # Sync filter to update all intermediate states
        self.sync()
        
        
    def update(self, meas, measType):
        # Determine expected measurement
        self.measExpected = np.array([self.az, self.el, self.rng, self.rngRate])
        
        # Compute residual
        self.meas = meas
        self.measType = measType
        self.measResidual = self.meas - self.measExpected
        
        # Compute sensitivity matrix
        H = measurements.inertialMeasurementSensitivity(self)
        self.measSensititivityMat = np.block([np.zeros((4,6)),-H])
        
        # Index based on measurement type
        self.measIndx = measurements.measType[self.measType]
        
        # Perform state update
        self.x, self.P = kf.measurementUpdate(self.x, 
                                              self.P, 
                                              self.numStates, 
                                              self.measResidual[self.measIndx], 
                                              self.measSensititivityMat[self.measIndx,:],
                                              self.measCov[self.measIndx,self.measIndx])
        
        # Sync filter to update all intermediate states
        self.sync()

        
    def sync(self):
        # Inertial Deputy State
        self.deputyPosInr = self.x[0:3]
        self.deputyVelInr = self.x[3:6]
        self.deputyCovInr = self.P[0:6,0:6]
        # Inertial Relative States
        self.relPosInr = self.x[6:9]
        self.relVelInr = self.x[9:12]
        self.relCovInr = self.P[6:12,6:12]
        # Chief states as derived from deputy and relative states
        self.chiefPosInr = self.deputyPosInr - self.relPosInr
        self.chiefVelInr = self.deputyVelInr - self.relVelInr
        self.chiefCovInr = self.deputyCovInr + self.relCovInr - self.P[0:6,6:12] - np.transpose(self.P[0:6,6:12])
        # Deputy and Chief cross covariance
        self.deputyChiefCrossCovInr  = self.deputyCovInr - self.P[0:6,6:12]
        # Inertial to RIC DCMs
        uKin.dcmInr2Ric(self.chiefPosInr, self.chiefVelInr, self.dcmInr2Ric)
        uKin.dcmInr2Ric(self.deputyPosInr, self.deputyVelInr, self.dcmInr2DepRic)
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
    
    
def ncvProcessNoise(dt, Q = np.eye(3)):
    return np.block([
                    [Q*(dt**3)/3, Q*(dt**2)/2],
                    [Q*(dt**2)/2, Q*dt       ]])

def stateUpdateInertial(dt, x, u, pert = None):
    r = x[0:3]
    v = x[3:6]
    if pert == None:
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

def stmRelative(dt, x):
    # Linearize the inertial motion about the anchor
    stmI = stmInertial(dt,x[0:6])
    return np.block([
            [stmI,               np.zeros((6, 6))  ],
            [np.zeros((6, 6))  , stmI              ]])

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

def dvProcessNoise(dv,varSf,varQ,varP):
    """
    Process noise due to delta-v maneuvers.
    
    Ref: Vaughn, Andrew Thomas, "A Monte-Carlo Performance Analysis of Kalman 
    Filter and Targeting Algorithms for Autonomous Orbital Rendezvous" (2004).
    https://dspace.mit.edu/entities/publication/a094c648-1ad3-451b-b4ae-839dfda52095

    Parameters
    ----------
    dv : 3x1 double
        Delta-V vector.
    varSf : double
        Scale factor variance.
    varQ : double
        Quantization variance.
    varP : double
        Pointing variance.

    Returns
    -------
    Qv : 6x6 double
        Delta-V process noise.

    """
    dvdvT = np.matmul(dv,np.transpose(dv))
    QSf = varSf*dvdvT
    QQ = varQ*np.eye(3)
    QP = varP*(np.matrix.trace(dvdvT)*np.eye(3) - dvdvT)
    return np.block([
            [np.zeros((3, 3)), np.zeros((3, 3))],
            [np.zeros((3, 3)), QSf + QQ + QP   ]])