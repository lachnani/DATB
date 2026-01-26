# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:47:48 2024

@author: Hakim Lachnani
"""

import numpy as np
from dynamics import orbit as orb
from dynamics import dynamicsUtils as uDyn


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
        self.dAmico = np.zeros((6,))
        self.rectClroe = np.zeros((6,))
        self.curvClroe = np.zeros((6,))
        
        # Declare DCM
        self.dcmInr2Ric = np.zeros((3,3))
        
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
            uDyn.rv2ric(self.chief.r, self.chief.v, self.deputy.r, self.deputy.v, self.relPosRectRic, self.relVelRectRic)
            uDyn.rectRic2curvRic(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic, self.relPosCurvRic, self.relVelCurvRic)
            relStateType = "RELSTATE_RECTRIC"
            
        elif frmType == "FORMATION_CHIEF_ANCHOR":
            # Chief anchor defined with chief and relative state
            self.t = chief.t
            self.chief = chief
            if relStateType == "RELSTATE_RECTRIC":
                self.relPosRectRic = relState[0:3]
                self.relVelRectRic = relState[3:6]
                uDyn.rectRic2curvRic(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic, self.relPosCurvRic, self.relVelCurvRic)
                self.deputy = orb.Orbit(self.chief.tJ2000, np.concatenate(ric2rv(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic)), pert = pert)
            elif relStateType == "RELSTATE_CURVRIC":
                # Convert the curvilinear ric state to rectilinear
                self.relPosCurvRic = relState[0:3]
                self.relVelCurvRic = relState[3:6]
                uDyn.curvRic2rectRic(self.chief.r, self.chief.v, self.relPosCurvRic, self.relVelCurvRic, self.relPosRectRic, self.relVelRectRic)
                self.deputy = orb.Orbit(self.chief.tJ2000, np.concatenate(ric2rv(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic)), pert = pert)
            elif relStateType == "RELSTATE_DOE":
                self.doe = relState
                self.deputy = orb.Orbit(self.chief.tJ2000, self.chief.oe + self.doe, stateType = "STATE_KEPEL", pert = pert)
            elif relStateType == "RELSTATE_DEE":
                self.dee = relState
                self.deputy = orb.Orbit(self.chief.tJ2000, self.chief.ee + self.dee, stateType = "STATE_EQEL", pert = pert)
            elif relStateType == "RELSTATE_RECT_CLROE":
                self.rectClroe = relState
                uDyn.clroe2ric(self.rectClroe, self.chief.meanMotion, 0, self.relPosRectRic, self.relVelRectRic)
                uDyn.rectRic2curvRic(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic, self.relPosCurvRic, self.relVelCurvRic)
                self.deputy = orb.Orbit(self.chief.tJ2000, np.concatenate(ric2rv(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic)), pert = pert)
            elif relStateType == "RELSTATE_CURV_CLROE":
                self.curvClroe = relState
                uDyn.clroe2ric(self.curvClroe, self.chief.meanMotion, 0, self.relPosCurvRic, self.relVelCurvRic)
                uDyn.curvRic2rectRic(self.chief.r, self.chief.v, self.relPosCurvRic, self.relVelCurvRic, self.relPosRectRic, self.relVelRectRic)
                self.deputy = orb.Orbit(self.chief.tJ2000, np.concatenate(ric2rv(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic)), pert = pert)
            
        elif frmType == "FORMATION_DEPUTY_ANCHOR":
            # Deputy anchor defined with deputy and relative state
            self.t = deputy.t
            self.deputy = deputy
            if relStateType == "RELSTATE_RECTRIC":
                self.relPosRectRic = relState[0:3]
                self.relVelRectRic = relState[3:6]
                uDyn.rectRic2curvRic(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic, self.relPosCurvRic, self.relVelCurvRic)
                self.chief = orb.Orbit(self.deputy.tJ2000, np.concatenate(ric2rv(self.deputy.r, self.deputy.v, -1*self.relPosRectRic, -1*self.relVelRectRic)), pert = pert)
            elif relStateType == "RELSTATE_CURVRIC":
                # Convert the curvilinear ric state to rectilinear
                self.relPosCurvRic = relState[0:3]
                self.relVelCurvRic = relState[3:6]
                uDyn.curvRic2rectRic(self.chief.r, self.chief.v, self.relPosCurvRic, self.relVelCurvRic, self.relPosRectRic, self.relVelRectRic)
                self.chief = orb.Orbit(self.deputy.tJ2000, np.concatenate(ric2rv(self.deputy.r, self.deputy.v, -1*self.relPosRectRic, -1*self.relVelRectRic)), pert = pert)
            elif relStateType == "RELSTATE_DOE":
                self.doe = relState
                self.deputy = orb.Orbit(self.deputy.tJ2000, self.deputy.oe - self.doe, stateType = "STATE_KEPEL", pert = pert)
            elif relStateType == "RELSTATE_DEE":
                self.dee = relState
                self.deputy = orb.Orbit(self.deputy.tJ2000, self.deputy.ee - self.dee, stateType = "STATE_EQEL", pert = pert)
            elif relStateType == "RELSTATE_RECT_CLROE":
                self.rectClroe = relState
                uDyn.clroe2ric(self.rectClroe, self.chief.meanMotion, 0, self.relPosRectRic, self.relVelRectRic)
                uDyn.rectRic2curvRic(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic, self.relPosCurvRic, self.relVelCurvRic)
                self.chief = orb.Orbit(self.deputy.tJ2000, np.concatenate(ric2rv(self.deputy.r, self.deputy.v, -1*self.relPosRectRic, -1*self.relVelRectRic)), pert = pert)
            elif relStateType == "RELSTATE_CURV_CLROE":
                self.curvClroe = relState
                uDyn.clroe2ric(self.curvClroe, self.chief.meanMotion, 0, self.relPosCurvRic, self.relVelCurvRic)
                uDyn.curvRic2rectRic(self.chief.r, self.chief.v, self.relPosCurvRic, self.relVelCurvRic, self.relPosRectRic, self.relVelRectRic)
                self.chief = orb.Orbit(self.deputy.tJ2000, np.concatenate(ric2rv(self.deputy.r, self.deputy.v, -1*self.relPosRectRic, -1*self.relVelRectRic)), pert = pert)
            
        # Convert relative representations
        if relStateType == "RELSTATE_DOE" or relStateType == "RELSTATE_DEE":
            uDyn.rv2ric(self.chief.r, self.chief.v, self.deputy.r, self.deputy.v, self.relPosRectRic, self.relVelRectRic)
            uDyn.rectRic2curvRic(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic, self.relPosCurvRic, self.relVelCurvRic)
        
        if relStateType != "RELSTATE_DEE":
            self.dee = self.deputy.ee - self.chief.ee
        
        if relStateType != "RELSTATE_DOE":
            self.doe = self.deputy.oe - self.chief.oe
            
        if relStateType != "RELSTATE_RECT_CLROE":
            uDyn.ric2clroe(self.relPosRectRic, self.relVelRectRic, self.chief.meanMotion, 0, self.rectClroe)
            
        if relStateType != "RELSTATE_CURV_CLROE":
            uDyn.ric2clroe(self.relPosCurvRic, self.relVelCurvRic, self.chief.meanMotion, 0, self.curvClroe)
            
        if relStateType != "RELSTATE_DAMICO":
            uDyn.oe2dAmico(self.chief.oe, self.deputy.oe, self.dAmico)
            
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
        uDyn.dcmInr2Ric(self.chief.r, self.chief.v, self.dcmInr2Ric)
        
        # Compute environment angles
        self.losEarthAng, self.losMoonAng, self.losSunAng = \
            uDyn.envAngles(self.chief.r, self.deputy.r, self.deputy.moon.r, self.deputy.sun.r)
        
        # Compute measurement parameters
        self.rng, self.rngRate, self.az, self.el = \
            uDyn.measParams(self.relPosRectRic, self.relVelRectRic)
        

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
        
        # Update State
        if np.linalg.norm(self.chief.aCtrlInRic) > np.finfo(float).eps:
            self.chief.aCtrlInEci = np.matmul(np.transpose(self.dcmInr2Ric), self.chief.aCtrlInRic)
            self.chief.dvTot += np.linalg.norm(self.chief.aCtrlInRic) * dt
        else:
            self.chief.aCtrlInEci = np.zeros((3,))
        if np.linalg.norm(self.deputy.aCtrlInRic) > np.finfo(float).eps:
            self.deputy.aCtrlInEci = np.matmul(np.transpose(self.dcmInr2Ric), self.deputy.aCtrlInRic)
            self.deputy.dvTot += np.linalg.norm(self.deputy.aCtrlInRic) * dt
        else:
            self.deputy.aCtrlInEci = np.zeros((3,))
        self.chief.propagate(dt)
        self.deputy.propagate(dt)
        
        # Update RIC frame 
        uDyn.dcmInr2Ric(self.chief.r, self.chief.v, self.dcmInr2Ric)
        
        # Convert relative states
        uDyn.rv2ric(self.chief.r, self.chief.v, self.deputy.r, self.deputy.v, self.relPosRectRic, self.relVelRectRic)
        uDyn.rectRic2curvRic(self.chief.r, self.chief.v, self.relPosRectRic, self.relVelRectRic, self.relPosCurvRic, self.relVelCurvRic)
        if (self.settings["relStates"] == True):
            # TODO: Fix to wrap angles appropriately...
            self.doe = self.deputy.oe - self.chief.oe
            self.dee = self.deputy.ee - self.chief.ee
            uDyn.oe2dAmico(self.chief.oe, self.deputy.oe, self.dAmico)
            uDyn.ric2clroe(self.relPosRectRic, self.relVelRectRic, self.chief.meanMotion, 0, self.rectClroe)
            uDyn.ric2clroe(self.relPosCurvRic, self.relVelCurvRic, self.chief.meanMotion, 0, self.curvClroe)
            
        # Compute Environment Parameters
        if (self.settings["environments"] == True):
            self.losEarthAng, self.losMoonAng, self.losSunAng = \
                uDyn.envAngles(self.chief.r, self.deputy.r, self.deputy.moon.r, self.deputy.sun.r)
            
        # Compute Measurement Parameters
        if (self.settings["measurements"] == True):
            self.rng, self.rngRate, self.az, self.el = \
                uDyn.measParams(self.relPosRectRic, self.relVelRectRic)
        

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
    uDyn.ric2rv(r, v, relPosRectRic, relVelRectRic, r_d, v_d)
    return r_d, v_d