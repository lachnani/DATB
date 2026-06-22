# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:47:48 2024

@author: Hakim Lachnani
"""

import numpy as np
from dynamics import orbit as orb
from kinematics import kinematicsUtils as uKin


class Formation():
    """ Base Formation Class
    
    Attributes
    ----------
    t : float
        time 
    chief: orbit
        chief orbit
    deputy: orbit
        deputy orbit
    relativeState: 6x1 float
        relative state
    """
    
    def __init__(
            self, 
            chief, deputy, relState, frmType = "FORMATION_DUAL_ABS",
            relStateType = "RELSTATE_RECTRIC", pert = None, settings = None
            ):
        
        # Declare relative states
        self.relPosRectRic = np.zeros((3,))
        self.relVelRectRic = np.zeros((3,))
        self.relPosCurvRic = np.zeros((3,))
        self.relVelCurvRic = np.zeros((3,))
        self.doe = np.zeros((6,))
        self.dee = np.zeros((6,))
        self.roe = np.zeros((6,))
        self.rectClroe = np.zeros((6,))
        self.curvClroe = np.zeros((6,))
        
        # Declare DCM
        self.dcmInr2Ric = np.zeros((3,3))
        self.dcmRic2Los = np.zeros((3,3))
        
        # Define default perturbations
        if (pert == None):
            pert = {
                "jnum": 0,
                "solarGrav": False,
                "lunarGrav": False,
                "SRP": False,
                "drag": False,
                "Cd": 0.0,
                "normalizedArea": 0.0
                }
            
        # Define default settings
        if (settings == None):
            orbitSettings = {
                "environments": True,
                "elements": True
                }
            settings = {
                "relStates": True,
                "environments": False,
                "measurements": False,
                "orbit": orbitSettings
                }
        
        # Initialize formation
        if frmType == "FORMATION_DUAL_ABS":
            # Dual inertial defined by two orbits
            self.t = chief.t
            self.chief = chief
            self.deputy = deputy
            self.deputy.t = self.t # Sync the epochs to the chief
            uKin.rv2ric(self.chief.r, self.chief.v, self.deputy.r, self.deputy.v, self.relPosRectRic, self.relVelRectRic)
            uKin.rectRic2curvRic(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic, self.relPosCurvRic, self.relVelCurvRic)
            relStateType = "RELSTATE_RECTRIC"
            
        elif frmType == "FORMATION_CHIEF_ANCHOR":
            # Chief anchor defined with chief and relative state
            self.t = chief.t
            self.chief = chief
            if relStateType == "RELSTATE_RECTRIC":
                self.relPosRectRic = relState[0:3]
                self.relVelRectRic = relState[3:6]
                uKin.rectRic2curvRic(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic, self.relPosCurvRic, self.relVelCurvRic)
                self.deputy = orb.Orbit(self.chief.tJ2000, np.concatenate(ric2rv(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic)), pert = pert)
            elif relStateType == "RELSTATE_CURVRIC":
                # Convert the curvilinear ric state to rectilinear
                self.relPosCurvRic = relState[0:3]
                self.relVelCurvRic = relState[3:6]
                uKin.curvRic2rectRic(self.chief.r, self.chief.v, self.relPosCurvRic, self.relVelCurvRic, self.relPosRectRic, self.relVelRectRic)
                self.deputy = orb.Orbit(self.chief.tJ2000, np.concatenate(ric2rv(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic)), pert = pert)
            elif relStateType == "RELSTATE_DOE":
                self.doe = relState
                self.deputy = orb.Orbit(self.chief.tJ2000, self.chief.oe + self.doe, stateType = "STATE_KEPEL", pert = pert)
            elif relStateType == "RELSTATE_DEE":
                self.dee = relState
                self.deputy = orb.Orbit(self.chief.tJ2000, self.chief.ee + self.dee, stateType = "STATE_EQEL", pert = pert)
            elif relStateType == "RELSTATE_RECT_CLROE":
                self.rectClroe = relState
                uKin.clroe2ric(self.rectClroe, self.chief.meanMotion, 0, self.relPosRectRic, self.relVelRectRic)
                uKin.rectRic2curvRic(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic, self.relPosCurvRic, self.relVelCurvRic)
                self.deputy = orb.Orbit(self.chief.tJ2000, np.concatenate(ric2rv(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic)), pert = pert)
            elif relStateType == "RELSTATE_CURV_CLROE":
                self.curvClroe = relState
                uKin.clroe2ric(self.curvClroe, self.chief.meanMotion, 0, self.relPosCurvRic, self.relVelCurvRic)
                uKin.curvRic2rectRic(self.chief.r, self.chief.v, self.relPosCurvRic, self.relVelCurvRic, self.relPosRectRic, self.relVelRectRic)
                self.deputy = orb.Orbit(self.chief.tJ2000, np.concatenate(ric2rv(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic)), pert = pert)
            
        elif frmType == "FORMATION_DEPUTY_ANCHOR":
            # Deputy anchor defined with deputy and relative state
            self.t = deputy.t
            self.deputy = deputy
            if relStateType == "RELSTATE_RECTRIC":
                self.relPosRectRic = relState[0:3]
                self.relVelRectRic = relState[3:6]
                uKin.rectRic2curvRic(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic, self.relPosCurvRic, self.relVelCurvRic)
                self.chief = orb.Orbit(self.deputy.tJ2000, np.concatenate(ric2rv(self.deputy.r, self.deputy.v, -1*self.relPosRectRic, -1*self.relVelRectRic)), pert = pert)
            elif relStateType == "RELSTATE_CURVRIC":
                # Convert the curvilinear ric state to rectilinear
                self.relPosCurvRic = relState[0:3]
                self.relVelCurvRic = relState[3:6]
                uKin.curvRic2rectRic(self.chief.r, self.chief.v, self.relPosCurvRic, self.relVelCurvRic, self.relPosRectRic, self.relVelRectRic)
                self.chief = orb.Orbit(self.deputy.tJ2000, np.concatenate(ric2rv(self.deputy.r, self.deputy.v, -1*self.relPosRectRic, -1*self.relVelRectRic)), pert = pert)
            elif relStateType == "RELSTATE_DOE":
                self.doe = relState
                self.deputy = orb.Orbit(self.deputy.tJ2000, self.deputy.oe - self.doe, stateType = "STATE_KEPEL", pert = pert)
            elif relStateType == "RELSTATE_DEE":
                self.dee = relState
                self.deputy = orb.Orbit(self.deputy.tJ2000, self.deputy.ee - self.dee, stateType = "STATE_EQEL", pert = pert)
            elif relStateType == "RELSTATE_RECT_CLROE":
                self.rectClroe = relState
                uKin.clroe2ric(self.rectClroe, self.chief.meanMotion, 0, self.relPosRectRic, self.relVelRectRic)
                uKin.rectRic2curvRic(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic, self.relPosCurvRic, self.relVelCurvRic)
                self.chief = orb.Orbit(self.deputy.tJ2000, np.concatenate(ric2rv(self.deputy.r, self.deputy.v, -1*self.relPosRectRic, -1*self.relVelRectRic)), pert = pert)
            elif relStateType == "RELSTATE_CURV_CLROE":
                self.curvClroe = relState
                uKin.clroe2ric(self.curvClroe, self.chief.meanMotion, 0, self.relPosCurvRic, self.relVelCurvRic)
                uKin.curvRic2rectRic(self.chief.r, self.chief.v, self.relPosCurvRic, self.relVelCurvRic, self.relPosRectRic, self.relVelRectRic)
                self.chief = orb.Orbit(self.deputy.tJ2000, np.concatenate(ric2rv(self.deputy.r, self.deputy.v, -1*self.relPosRectRic, -1*self.relVelRectRic)), pert = pert)
            
        # Convert relative representations
        if relStateType == "RELSTATE_DOE" or relStateType == "RELSTATE_DEE":
            uKin.rv2ric(self.chief.r, self.chief.v, self.deputy.r, self.deputy.v, self.relPosRectRic, self.relVelRectRic)
            uKin.rectRic2curvRic(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic, self.relPosCurvRic, self.relVelCurvRic)
        
        if relStateType != "RELSTATE_DEE":
            self.dee = self.deputy.ee - self.chief.ee
        
        if relStateType != "RELSTATE_DOE":
            self.doe = self.deputy.oe - self.chief.oe
            
        if relStateType != "RELSTATE_RECT_CLROE":
            uKin.ric2clroe(self.relPosRectRic, self.relVelRectRic, self.chief.meanMotion, 0, self.rectClroe)
            
        if relStateType != "RELSTATE_CURV_CLROE":
            uKin.ric2clroe(self.relPosCurvRic, self.relVelCurvRic, self.chief.meanMotion, 0, self.curvClroe)
            
        if relStateType != "RELSTATE_ROE":
            uKin.oe2roe(self.chief.oe, self.deputy.oe, self.roe)
            
        # Set perturbations
        self.chief.pert = pert
        self.deputy.pert = pert
        
        # Synch epochs
        self.deputy.tJ2000 = self.chief.tJ2000
        
        # Set settings
        self.settings = settings
        self.chief.settings = self.settings["orbit"]
        self.deputy.settings = self.settings["orbit"]
        
        # Define RIC frame 
        uKin.dcmInr2Ric(self.chief.r, self.chief.v, self.dcmInr2Ric)
        
        # Define the LOS frame
        uKin.dcmRic2Los(self.relPosRectRic, self.dcmRic2Los)
        
        # Compute environment angles
        self.losEarthAng, self.losMoonAng, self.losSunAng = \
            uKin.envAngles(self.chief.r, self.deputy.r, self.deputy.moon.r, self.deputy.sun.r)
        
        # Compute measurement parameters
        self.rng, self.rngRate, self.az, self.el = \
            uKin.measParams(self.relPosRectRic, self.relVelRectRic)
        

    def propagate(self, dt):
        """
        Propagates formation

        Parameters
        ----------
        dt : float
            delta time for propagation

        Returns
        -------
        None.

        """
        
        # Update Time
        self.t = self.t + dt
        
        # Increment Delta-V counter if accelerating
        self.chief.dvTot += np.linalg.norm(self.chief.aCtrlInEci) * dt
        self.deputy.dvTot += np.linalg.norm(self.deputy.aCtrlInEci) * dt
        
        # Update State
        self.chief.propagate(dt)
        self.deputy.propagate(dt)
        
        # Update RIC frame 
        uKin.dcmInr2Ric(self.chief.r, self.chief.v, self.dcmInr2Ric)
        
        # Convert relative states
        uKin.rv2ric(self.chief.r, self.chief.v, self.deputy.r, self.deputy.v, self.relPosRectRic, self.relVelRectRic)
        uKin.rectRic2curvRic(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic, self.relPosCurvRic, self.relVelCurvRic)
        if (self.settings["relStates"] == True):
            # TODO: Fix to wrap angles appropriately...
            self.doe = self.deputy.oe - self.chief.oe
            self.dee = self.deputy.ee - self.chief.ee
            uKin.oe2roe(self.chief.oe, self.deputy.oe, self.roe)
            uKin.ric2clroe(self.relPosRectRic, self.relVelRectRic, self.chief.meanMotion, 0, self.rectClroe)
            uKin.ric2clroe(self.relPosCurvRic, self.relVelCurvRic, self.chief.meanMotion, 0, self.curvClroe)
            
        # Compute Environment Parameters
        if (self.settings["environments"] == True):
            self.losEarthAng, self.losMoonAng, self.losSunAng = \
                uKin.envAngles(self.chief.r, self.deputy.r, self.deputy.moon.r, self.deputy.sun.r)
            
        # Compute Measurement Parameters
        if (self.settings["measurements"] == True):
            self.rng, self.rngRate, self.az, self.el = \
                uKin.measParams(self.relPosRectRic, self.relVelRectRic)
            uKin.dcmRic2Los(self.relPosRectRic, self.dcmRic2Los)
        

def ric2rv(r, v, relPosRectRic, relVelRectRic):
    """
    Map the deputy position and velocity vector relative to the chief RIC frame to inertial frame.

    :param r: chief inertial position
    :param v: chief inertial velocity
    :param relPosRectRic: RIC frame relative position
    :param relVelRectRic: RIC frame relative velocity
    :return: r_d: deputy inertial position
             v_d: deputy inertial velocity
    """

    r_d = np.zeros((3,))
    v_d = np.zeros((3,))
    uKin.ric2rv(r, v, relPosRectRic, relVelRectRic, r_d, v_d)
    return r_d, v_d