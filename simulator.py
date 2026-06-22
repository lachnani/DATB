# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 15:38:52 2025

@author: Hakim Lachnani
"""

from scipy.spatial.transform import Rotation as R
import numpy as np
import time

from dynamics import formation as frm
from dynamics import orbit as orb
from planning import passiveSafety as ps
import logger


class Simulator():
    """ Base Simulation Class
    
    Attributes
    ----------
    t : float
        time 
    frm: Formation
        Formation with chief and deputy
    fsw: FlightSoftware
        deputy flight software
    """
    
    def __init__(
            self, settings,
            formation, flightsoftware, quiet = False
            ):
        
        # Initialize simulation
        self.t = 0 
        self.cycle = 1
        self.quiet = quiet
        
        # Save start time
        self.startTime = time.time()
        
        self.settings = settings
        self.kovType = self.settings["dynamics"]["kov"]["type"]
        if self.kovType == "PRISM" or self.kovType == "PRISM_CURV":
            self.kovLim = np.array([self.settings["dynamics"]["kov"]["R"], \
                                    self.settings["dynamics"]["kov"]["I"], \
                                    self.settings["dynamics"]["kov"]["C"]])
        elif self.kovType == "SPHERE":
            self.kovLim = self.settings["dynamics"]["kov"]["range"]
        elif self.kovType == "ELIPSOID":
            # TODO
            self.kovLim = 0
        
        self.frm = formation
        self.frm.settings = self.settings["formation"]
        self.fsw = flightsoftware 
        
        # Burn perturbations
        self.burnPertAngDeg = 0.0 
        self.burnPertMagPerc = 0.0
        
        self.initStatus()
        if self.status and self.settings["log"]["status"]:
            # Initialize log
            self.log = logger.Log(np.int64(np.ceil(self.settings["simDuration"]/self.settings["log"]["dt"]))+1, self.settings)
            self.log.kovType = self.kovType
            self.log.kovLim = self.kovLim
            self.updateLog()
            
        if self.status == True and not self.quiet:
            print("SIM:", settings["name"], "initialized")
        
    def initStatus(self):
        """
        Check initial sim status

        Returns
        -------
        None.

        """
        
        self.status = True
        
        # Confirm dynamics are on
        if self.settings["dynamics"]["status"] == False:
            self.status = False
        
        # Confirm dynamics dt is the smallest
        if ((self.settings["dynamics"]["dt"] > self.settings["log"]["dt"]) or
            (self.settings["dynamics"]["dt"] > self.settings["fsw"]["dt"])):
            self.status = False
            
        # Confirm all dts are a multiple of dynamics dt
        if ((self.settings["log"]["dt"] % self.settings["dynamics"]["dt"] != 0) or
            (self.settings["fsw"]["dt"] % self.settings["dynamics"]["dt"] != 0)):
            self.status = False
            
    
    def checkStatus(self):
        """
        Check sim status

        Returns
        -------
        None.

        """
        
        # Confirm all the times agree
        if ((self.t != self.frm.t) or
            (self.t != self.frm.chief.t) or
            (self.t != self.frm.deputy.t)):
            self.status = False
            
        # End the sim if beyond the pre-defined sim duration
        if self.t > self.settings["simDuration"]:
            self.terminate()
        
    
    def run(self, t):
        """
        Runs simulation to given time

        Parameters
        ----------
        t : float
            time to run until

        Returns
        -------
        None.

        """
        
        while self.t < t and self.status == True:
            
            self.cycle += 1
            
            # Propagate dynamics
            if self.settings["dynamics"]["status"] == True:
                self.frm.propagate(self.settings["dynamics"]["dt"])
                self.t = self.t + self.settings["dynamics"]["dt"]
                # Check whether the KOV is breached
                if not self.log.kovBreach:
                    if self.kovType == "PRISM_CURV":
                        if ps.checkKovBreach(self.frm.relPosCurvRic, self.kovLim, "PRISM"):
                            self.log.kovBreach = True
                            self.log.tKovBreach = self.t
                    else:
                        if ps.checkKovBreach(self.frm.relPosRectRic, self.kovLim, self.kovType):
                            self.log.kovBreach = True
                            self.log.tKovBreach = self.t
                
            # Propagate FSW
            if self.settings["fsw"]["status"] == True and (self.t % self.settings["fsw"]["dt"] == 0):
                self.fsw.propagate(self.fsw.dt)
                
            # Log
            if self.settings["log"]["status"] == True and (self.t % self.settings["log"]["dt"] == 0):
                self.log.i += 1
                self.updateLog()
            
            self.checkStatus()
                
    def updateLog(self):
        
        # Times
        self.log.time[self.log.i] = np.array((self.t, self.frm.chief.tJ2000))
        
        # Ephemerides
        if (self.settings["formation"]["orbit"]["environments"] == True):
            self.log.sunPosEci[self.log.i] = self.frm.deputy.sun.r
            self.log.moonPosEci[self.log.i] = self.frm.deputy.moon.r
        
        # Vehicle Parameters
        self.log.posVehEci[self.log.i]      = self.frm.deputy.r
        self.log.velVehEci[self.log.i]      = self.frm.deputy.v
        self.log.accThrVehEci[self.log.i]   = self.frm.deputy.aCtrlInEci
        self.log.accThrVehRic[self.log.i]   = np.matmul(self.frm.dcmInr2Ric, self.frm.deputy.aCtrlInEci)
        self.log.dvTotVeh[self.log.i]       = self.frm.deputy.dvTot
        if (self.settings["formation"]["orbit"]["elements"] == True):
            self.log.oeVeh[self.log.i]          = self.frm.deputy.oe
            self.log.eeVeh[self.log.i]          = self.frm.deputy.ee
        if (self.settings["formation"]["orbit"]["environments"] == True):
            self.log.eclipseVeh[self.log.i]     = self.frm.deputy.inEclipse
        
        # RSO Paramaters
        self.log.posRsoEci[self.log.i]      = self.frm.chief.r
        self.log.velRsoEci[self.log.i]      = self.frm.chief.v
        self.log.accThrRsoEci[self.log.i]   = self.frm.chief.aCtrlInEci
        self.log.accThrRsoRic[self.log.i]   = np.matmul(self.frm.dcmInr2Ric, self.frm.chief.aCtrlInEci)
        self.log.dvTotRso[self.log.i]       = self.frm.chief.dvTot
        if (self.settings["formation"]["orbit"]["elements"] == True):
            self.log.oeRso[self.log.i]          = self.frm.chief.oe 
            self.log.eeRso[self.log.i]          = self.frm.chief.ee
        if (self.settings["formation"]["orbit"]["environments"] == True):
            self.log.eclipseRso[self.log.i]     = self.frm.chief.inEclipse
        
        # Relative Parameters
        self.log.relPosRectRic[self.log.i]  = self.frm.relPosRectRic
        self.log.relVelRectRic[self.log.i]  = self.frm.relVelRectRic
        self.log.relPosCurvRic[self.log.i]  = self.frm.relPosCurvRic
        self.log.relVelCurvRic[self.log.i]  = self.frm.relVelCurvRic
        self.log.qInrToRic[self.log.i]      = R.from_matrix(self.frm.dcmInr2Ric).as_quat()
        if (self.settings["formation"]["relStates"] == True):
            self.log.doe[self.log.i]            = self.frm.doe
            self.log.dee[self.log.i]            = self.frm.dee
            self.log.roe[self.log.i]            = self.frm.roe
            self.log.rectClroe[self.log.i]      = self.frm.rectClroe
            self.log.curvClroe[self.log.i]      = self.frm.curvClroe
        if (self.settings["formation"]["environments"] == True):
            self.log.losAngles[self.log.i] = np.array((self.frm.losEarthAng, self.frm.losSunAng, self.frm.losMoonAng))
        if (self.settings["formation"]["measurements"] == True):
            self.log.measParams[self.log.i] = np.array((self.frm.rng, self.frm.rngRate, self.frm.az, self.frm.el))
            self.log.qRicToLos[self.log.i]  = R.from_matrix(self.frm.dcmRic2Los).as_quat()
        
        if self.settings["fsw"]["status"] == True:
            # Log FSW
            temp = 0
    
    def terminate(self):
        
        if self.status == True:
            self.status = False
            if not self.quiet:
                print("SIM:", self.settings["name"], "terminating successfully")
        else:
            print("SIM:", self.settings["name"], "failed!")
            print("SIM:", self.settings["name"], "terminating")
            
        if not self.quiet:
            endTime = time.time()
            elapsed_time = endTime - self.startTime
            print("SIM: Simulation completed in %.1f seconds" % elapsed_time)
            print("SIM: Realtime factor: %.2f" % (self.t/elapsed_time))
            print("-----------------------------------------------------")
            
    def impulsiveBurn(self, dvEci, sv = "Deputy"):
        """
        Executes impulsive burn
        """
        # Perturb Delta-V
        dvEciPert = pertVec(dvEci, self.burnPertAngDeg, self.burnPertMagPerc)
        
        # Apply delta-v by moving the inertial velocity and incrementing the
        # accumulated delta-v
        if sv == "Deputy":
            self.frm.deputy.v       += dvEciPert
            self.frm.deputy.dvTot   += np.linalg.norm(dvEciPert)
        else:
            self.frm.chief.v        += dvEciPert
            self.frm.chief.dvTot    += np.linalg.norm(dvEciPert)        
            
    
def pertVec(v, angDeg, magPerc):
    """
    Randomly perturbs vector. Assumes Gaussian perturbation distributions, 
    and 3-sigma input magnitudes.
    """
    if magPerc > 0.0:
        magPert = np.random.normal(0.0, 0.01*magPerc/3) + 1 # Magnitude perturbation as ratio
    else:
        magPert = 1.0;
    
    if angDeg > 0.0:
        angPert = np.random.normal(0.0, angDeg/3) # Angle perturbation in degrees
        randVec = np.random.rand(3,)
        rotVec = np.cross(v, randVec) # Forces an orthogonal rotation vector
        rotVec = rotVec / np.linalg.norm(rotVec) # Normalize
        r = R.from_rotvec(angPert * rotVec, degrees=True)   
        return magPert*r.apply(v)
    else:
        return magPert*v
    
    
    
    
    