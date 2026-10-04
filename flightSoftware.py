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

class Navigation:
    
    
    def __init__(
            self,
            tJ2000):
        
        # Initialize Time
        self.tJ2000 = tJ2000
        
        # Set default constants
        self.mu = orb.MU_EARTH
        
        # Set filter status
        self.fltrInit = False
        self.fltrConverged = False
        self.fltrDiverged = False
        self.fltrCorrupted = False
        
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
        
        # Update output states
        self.sync()
        
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
        
        # Mean Motion (TODO: Make uKin util)
        if self.deputyOrbEl[0] > 0:
            self.deputyMeanMotion = np.sqrt(self.mu/self.deputyOrbEl[0]**3)
        else:
            self.deputyMeanMotion = 0 
            
        if self.chiefOrbEl[0] > 0:
            self.chiefMeanMotion = np.sqrt(self.mu/self.chiefOrbEl[0]**3)
        else:
            self.chiefMeanMotion = 0 
            
        # Semimajor Axis variance
        self.deputySmaVar = self.smaVariance(self.mu, self.deputyPosInr, self.deputyVelInr, self.deputyOrbEl[0], self.fltr.deputyCovInr)
        self.chiefSmaVar = self.smaVariance(self.mu, self.chiefPosInr, self.chiefVelInr, self.chiefOrbEl[0], self.fltr.chiefCovInr)
    
        # Eclipse status
        self.deputyInEclipse = uDyn.eclipse(self.deputyPosInr, self.fltr.sun.rUnit)
        self.chiefInEclipse = uDyn.eclipse(self.chiefPosInr, self.fltr.sun.rUnit)
        
        # RIC Frame
        uKin.dcmInr2Ric(self.chiefPosInr, self.chiefVelInr, self.dcmInr2Ric)
        
        # Convert relative states
        uKin.rv2ric(self.chiefPosInr, self.chiefVelInr, self.deputyPosInr, self.deputyVelInr, self.relPosRectRic, self.relVelRectRic)
        uKin.rectRic2curvRic(self.chiefPosInr, self.chiefVelInr, self.relPosRectRic, self.relVelRectRic, self.relPosCurvRic, self.relVelCurvRic)
        # TODO: Fix to wrap angles appropriately...
        self.diffOrbEl = self.deputyOrbEl - self.chiefOrbEl
        self.diffEqEl = self.deputyEqEl - self.chiefEqEl
        uKin.oe2roe(self.chiefOrbEl, self.deputyOrbEl, self.relOrbEl)
        uKin.ric2clroe(self.relPosRectRic, self.relVelRectRic, self.chiefMeanMotion, 0, self.rectClroe)
        uKin.ric2clroe(self.relPosCurvRic, self.relVelCurvRic, self.chiefMeanMotion, 0, self.curvClroe)
            
        # Compute Environment Parameters
        self.losEarthAng, self.losMoonAng, self.losSunAng = \
            uKin.envAngles(self.chiefPosInr, self.deputyPosInr, self.fltr.moon.r, self.fltr.sun.r)
        self.visibilityAng = uKin.visibilityAngle(self.chiefPosInr, self.deputyPosInr, self.fltr.sun.r)
          
        # Compute measurement parameters
        self.rng = la.norm(self.relPosRectRic)
        self.rngRate = np.dot(self.relPosRectRic, self.relVelRectRic) / self.rng
        
        
    def smaVariance(mu, r, v, a, P):
        """
        Computes semimajor axis variance from state estimates
        
        Ref: Carpenter and D'Souza, "Navigation Filter Best Practices"
        
        TODO: Move to kinematicsUtils.c

        Parameters
        ----------
        mu : double
            gravitational parameter.
        r : 3x1 double
            inertial position.
        v : 3x1 double
            inertial velocity.
        a : double
            semimajor axis estimate.
        P : 6x6 double
            inertial covariance.

        Returns
        -------
        smaVar: double
            semi-major axis variance.

        """
        
        Fa = 2*a**2*np.block([np.transpose(r)/la.norm(r)**3,np.transpose(v)/mu])
        return np.matmul(Fa, np.matmul(P, np.transpose(Fa)))