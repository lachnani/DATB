# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 11:01:08 2026

@author: Hakim Lachnani
"""

import numpy as np
from numpy import linalg as la
from dynamics import dynamicsUtils as uDyn
from dynamics import orbit as orb
from kinematics import kinematicsUtils as uKin
from navigation import estimator
from navigation import kalmanFilter

class Main:
    
    def __init__(
            self,
            tJ2000, dt, # Common parameters
            bufferDepth # Navigation parameters
            ):
        
        # Initialize Time
        self.tJ2000 = tJ2000
        self.dt = dt
        
        # Initialize Navigation Struct
        self.nav = Navigation(self.tJ2000, bufferDepth)
        
        # Initialize Guidance Struct
        self.guid = 0 
        
        # Initialize Control Struct
        self.ctrl = Control(self.tJ2000)        
        
    def cycle(self, 
              measAvailable = False, meas = np.zeros((4,)), measType = "anglesRange"):
        
        # Propagate time
        self.tJ2000 = self.tJ2000 + self.dt
        
        # Navigation
        if self.nav.fltrInit == True:
            self.nav.propagate(self.dt, self.ctrl.aCtrlInEci)
            if measAvailable == True:
                self.nav.update(meas, measType)
                # self.nav.runChecks()
            self.nav.sync()
        
        # Guidance
        
        # Control
        
    
        
    
class Navigation:   
    
    def __init__(
            self,
            tJ2000, bufferDepth):
        
        # Initialize Time
        self.tJ2000 = tJ2000
        
        # Set default constants
        self.mu = orb.MU_EARTH
        
        # Set filter status
        self.fltrInit = False
        self.fltrConverged = False
        self.fltrDiverged = False
        self.fltrCorrupted = False
        self.fltrConsistent = True
        
        # Initialize deputy states
        self.deputyPosInr = np.zeros((3,))
        self.deputyVelInr = np.zeros((3,))
        self.deputyEqEl = np.zeros((6,))
        self.deputyOrbEl = np.zeros((6,))
        self.deputyMeanMotion = 0.0
        self.deputySmaVar = 0.0 
        self.deputyInEclipse = False
        
        # Initialize chief states
        self.chiefPosInr = np.zeros((3,))
        self.chiefVelInr = np.zeros((3,))
        self.chiefEqEl = np.zeros((6,))
        self.chiefOrbEl = np.zeros((6,))
        self.chiefMeanMotion = 0.0
        self.chiefSmaVar = 0.0 
        self.chiefInEclipse = False
        
        # Initialize RIC frame
        self.dcmInr2Ric = np.zeros((3,3))
        self.dcmRic2Los = np.zeros((3,3))
        self.dcmInr2Los = np.zeros((3,3))
        
        # Initialize Relative States
        self.relPosRectRic = np.zeros((3,))
        self.relVelRectRic = np.zeros((3,))
        self.relPosCurvRic = np.zeros((3,))
        self.relvelCurvRic = np.zeros((3,))
        self.diffOrbEl = np.zeros((6,))
        self.diffEqEl = np.zeros((6,))
        self.relOrbEl = np.zeros((6,))
        self.rectClroe = np.zeros((6,))
        self.curvClroe = np.zeros((6,))
        
        # Initialize environmental parameters
        self.losEarthAng = 0.0 
        self.losMoonAng = 0.0 
        self.losSunAng = 0.0 
        
        # Initialize measurement parameters
        self.rng = 0.0 
        self.rngRate = 0.0 
        
        # Initialize Buffer
        self.bffr = kalmanFilter.Buffer(bufferDepth, 4, 12)
        
        
    def configureFilter(self, 
                        deputyProcNoiseRic, chiefProcNoiseRic, relProcNoiseRic,
                        dvVar, measCov):
        
        self.fltrDeputyProcNoiseRic = deputyProcNoiseRic
        self.fltrChiefProcNoiseRic = chiefProcNoiseRic
        self.fltrRelProcNoiseRic = relProcNoiseRic
        self.fltrDvVar = dvVar
        self.fltrMeasCov = measCov
        

    def initFilter(self, filterType,
                   deputyPosInr, deputyVelInr, deputyCovInr, 
                   chiefPosInr, chiefVelInr, chiefCovInr
                   ):
        
        # Assign filter type
        self.fltrType = filterType
        
        # Initialize filter
        if self.filterType == "DualInertial":
            self.fltr = estimator.DualInertialFilter(
                self.tJ2000, 
                chiefPosInr, chiefVelInr, chiefCovInr, 
                deputyPosInr, deputyVelInr, deputyCovInr, 
                self.fltrChiefProcNoiseRic, self.fltrDeputyProcNoiseRic, 
                self.fltrDvVar, self.fltrMeasCov)
        elif self.filterType == "InertialRelative":
            self.fltr =  estimator.InertialRelativeFilter(
                self.tJ2000, 
                chiefPosInr, chiefVelInr, chiefCovInr, 
                deputyPosInr, deputyVelInr, deputyCovInr, 
                self.fltrDeputyProcNoiseRic, self.fltrRelProcNoiseRic,
                self.fltrDvVar, self.fltrMeasCov)
        elif self.filterType == "ChiefDecoupledDualInertial":
            self.fltr = estimator.ChiefDecoupledDualInertialFilter(
                self.tJ2000, 
                chiefPosInr, chiefVelInr, chiefCovInr, 
                deputyPosInr, deputyVelInr, deputyCovInr, 
                self.fltrChiefProcNoiseRic, self.fltrDeputyProcNoiseRic, 
                self.fltrDvVar, self.fltrMeasCov)
        elif self.filterType == "RelativeDecoupledInertialRelative":
            self.fltr =  estimator.RelativeDecoupledInertialRelativeFilter(
                self.tJ2000, 
                chiefPosInr, chiefVelInr, chiefCovInr, 
                deputyPosInr, deputyVelInr, deputyCovInr, 
                self.fltrDeputyProcNoiseRic, self.fltrRelProcNoiseRic,
                self.fltrDvVar, self.fltrMeasCov)
        
        # Mark the filter as initialized
        self.fltrInit = True
        
        # Sync filter states
        self.fltr.sync()
        
        
    def propagate(self, dt, aCtrlInEci):
        
        # Update time
        self.tJ2000 = self.tJ2000 + dt
        
        # Propagate Nav Filter
        self.fltr.propagate(dt, aCtrlInEci)
        
        # Update output states
        self.sync()
        
    def update(self, meas, measType):
        
        # Update Nav Filter
        self.fltr.update(meas, measType)
        
        # Update buffer
        self.bffr.update(self.fltr.tJ2000, self.fltr.x, self.fltr.measResidual, self.fltr.nis)
        
        # Update output states
        self.sync()
        
    def checkConvergence(self):
        
        # TODO: Add a check
        self.fltrConverged = False
        
    def checkDivergence(self):
        
        # TODO: Add a check
        self.fltrDiverged = False
        
    def checkCorruption(self):
        
        # TODO: Add a check
        self.fltrCorrupted = False
        
    def checkConsistency(self):
        
        # TODO: Add a check
        self.fltrConsistent = True
        
    def sync(self):
        
        # Deputy and Chief inertial states
        self.deputyPosInr = self.fltr.deputyPosInr
        self.deputyVelInr = self.fltr.deputyVelInr
        self.chiefPosInr = self.fltr.chiefPosInr
        self.chiefVelInr = self.fltr.chiefVelInr 
        
        # Deputy and Chief orbit elements
        uKin.rv2ee(self.mu, self.deputyPosInr, self.deputyVelInr, self.deputyEqEl)
        uKin.rv2oe(self.mu, self.deputyPosInr, self.deputyVelInr, self.deputyOrbEl)
        uKin.rv2ee(self.mu, self.chiefPosInr, self.chiefVelInr, self.chiefEqEl)
        uKin.rv2oe(self.mu, self.chiefPosInr, self.chiefVelInr, self.chiefOrbEl)
        
        # Mean Motion 
        self.deputyMeanMotion = uKin.meanMotion(self.mu, self.deputyOrbEl[0])
        self.chiefMeanMotion = uKin.meanMotion(self.mu, self.chiefOrbEl[0]) 
            
        # Semimajor Axis variance
        self.deputySmaVar = uKin.smaVariance(self.mu, self.deputyPosInr, self.deputyVelInr, self.deputyOrbEl[0], self.fltr.deputyCovInr)
        self.chiefSmaVar = uKin.smaVariance(self.mu, self.chiefPosInr, self.chiefVelInr, self.chiefOrbEl[0], self.fltr.chiefCovInr)
    
        # Eclipse status
        self.deputyInEclipse = uDyn.eclipse(self.deputyPosInr, self.fltr.sun.rUnit)
        self.chiefInEclipse = uDyn.eclipse(self.chiefPosInr, self.fltr.sun.rUnit)
        
        # Convert relative states
        uKin.rv2ric(self.chiefPosInr, self.chiefVelInr, self.deputyPosInr, self.deputyVelInr, self.relPosRectRic, self.relVelRectRic)
        uKin.rectRic2curvRic(self.chiefPosInr, self.chiefVelInr, self.relPosRectRic, self.relVelRectRic, self.relPosCurvRic, self.relVelCurvRic)
        # TODO: Fix to wrap angles appropriately...
        self.diffOrbEl = self.deputyOrbEl - self.chiefOrbEl
        self.diffEqEl = self.deputyEqEl - self.chiefEqEl
        uKin.oe2roe(self.chiefOrbEl, self.deputyOrbEl, self.relOrbEl)
        uKin.ric2clroe(self.relPosRectRic, self.relVelRectRic, self.chiefMeanMotion, 0, self.rectClroe)
        uKin.ric2clroe(self.relPosCurvRic, self.relVelCurvRic, self.chiefMeanMotion, 0, self.curvClroe)
        
        # Frames
        uKin.dcmInr2Ric(self.chiefPosInr, self.chiefVelInr, self.dcmInr2Ric)
        uKin.dcmRic2Los(self.relPosRectRic, self.dcmRic2Los)
        self.dcmInr2Los = np.matmul(self.dcmRic2Los,self.dcmInr2Ric)
            
        # Compute Environment Parameters
        self.losEarthAng, self.losMoonAng, self.losSunAng = \
            uKin.envAngles(self.chiefPosInr, self.deputyPosInr, self.fltr.moon.r, self.fltr.sun.r)
        self.visibilityAng = uKin.visibilityAngle(self.chiefPosInr, self.deputyPosInr, self.fltr.sun.r)
          
        # Compute measurement parameters
        self.rng = la.norm(self.relPosRectRic)
        self.rngRate = np.dot(self.relPosRectRic, self.relVelRectRic) / self.rng
        
    
class Control:
    
    def __init__(
            self,
            tJ2000):
        
        self.tJ2000 = tJ2000
        
        self.aCtrlInEci = np.zeros((3,))