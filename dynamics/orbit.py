"""
Created on 02/08/2025

@author: Hakim Lachnani

Orbit.py contains the orbit class
    
"""

# import required modules:
import math
import numpy as np
from dynamics import ephemerides as eph
from dynamics import dynamicsUtils as uDyn
from kinematics import kinematicsUtils as uKin

MU_EARTH = 398600.436 #km^3/s^2
D2R = np.pi/180.0

class Orbit:
    def __new__(cls, *args, **kwargs):
        # Create a new orbit
        return super().__new__(cls)

    def __init__(self, 
              epoch, state, stateType = "STATE_POSVEL", pert = None, settings = None):
        # Initialize a new orbit from position and velocity vectors, classic orbital elements, or equinoctial elements
        self.r = np.zeros((3,))
        self.v = np.zeros((3,))
        self.oe = np.zeros((6,))
        self.ee = np.zeros((6,))
        
        # Define default settings
        if (settings == None):
            settings = {
                "environments": True,
                "elements": True
                }
        self.settings = settings
        
        # Define the central body
        self.mu = MU_EARTH

        # Initialize the state and convert to other representations
        if stateType == "STATE_POSVEL":
            self.r = state[0:3]
            self.v = state[3:6]
            uKin.rv2ee(self.mu, self.r, self.v, self.ee)
            uKin.rv2oe(self.mu, self.r, self.v, self.oe)
        elif stateType == "STATE_KEPEL":
            self.oe = state
            uKin.oe2ee(self.oe, self.ee)
            uKin.oe2rv(self.mu, self.oe, self.r, self.v)
        elif stateType == "STATE_EQEL":
            self.ee = state
            uKin.ee2rv(self.mu, self.ee, self.r, self.v)
            uKin.rv2oe(self.mu, self.r, self.v, self.oe)      
            
        self.meanMotion = np.sqrt(self.mu/self.oe[0]**3)

        #Define the epoch
        self.t = 0
        self.tJ2000 = epoch
        
        #Define whether perturbations are enabled
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
        self.pert = pert
        
        # Define Drag settings
        if self.pert["drag"]:
            self.Cd = self.pert["Cd"]
            self.normA = self.pert["normalizedArea"]
        else:
            self.Cd = 0.0
            self.normA = 0.0
        
        #Initialize sun and moon ephemeris
        self.sun = eph.SunEphemeris(self.tJ2000)
        self.moon = eph.MoonEphemeris(self.tJ2000)    
        
        if (self.settings["environments"] == True):
            self.inEclipse = uDyn.eclipse(self.r, self.sun.rUnit)
        else:
            self.inEclipse = False
        
        #Set the control acceleration to 0
        self.aCtrlInEci = np.zeros(3)
        
        #Set the accumulated Delta-V to 0
        self.dvTot = 0.0
        
    def propagate(self, dt):
        """
        Propagates orbit using RK4

        Parameters
        ----------
        dt : float
            delta time for propagation

        Returns
        -------
        None.

        """
        
        uDyn.Orbit_rk4(self.pert["solarGrav"], self.pert["lunarGrav"], self.pert["drag"], self.pert["jnum"], \
                       self.moon.r, self.sun.r, self.Cd, self.normA, \
                       self.aCtrlInEci, dt, self.r, self.v)
        self.t = self.t + dt
        self.tJ2000 = self.tJ2000 + dt
        
        if (self.settings["elements"] == True):
            uKin.rv2ee(self.mu, self.r, self.v, self.ee)
            uKin.rv2oe(self.mu, self.r, self.v, self.oe)
            if self.oe[0] > 0:
                self.meanMotion = np.sqrt(self.mu/self.oe[0]**3)
            else:
                self.meanMotion = 0 
        
        self.sun.update(self.tJ2000)
        self.moon.update(self.tJ2000)
        if (self.settings["environments"] == True):
            self.inEclipse = uDyn.eclipse(self.r, self.sun.rUnit)
