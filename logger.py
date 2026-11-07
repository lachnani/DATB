# -*- coding: utf-8 -*-
"""
Created on Sat Mar 22 20:30:54 2025

@author: Hakim Lachnani
"""

import numpy as np

class Log():
    """ Base Log Class
    
    Attributes
    ----------
    i : integer
        index
    df: DataFrame
    """
    
    def __init__(self, size, settings):
        """
        Initialize the log

        """ 
        
        self.i = 0
        
        # Times
        self.simTime = np.zeros((size,))
        self.tJ2000  = np.zeros((size,))
        
        # Keep-out volume
        self.kovType = "PRISM"
        self.kovLim = 0
        self.kovBreach = False
        self.tKovBreach = 0
        
        # Ephemerides
        if (settings["formation"]["orbit"]["environments"] == True):
            self.sunPosEci = np.zeros((3,size))
            self.moonPosEci = np.zeros((3,size))
      
        # Vehicle Parameters
        self.posVehEci      = np.zeros((3,size))
        self.velVehEci      = np.zeros((3,size))
        self.accThrVehEci   = np.zeros((3,size))
        self.accThrVehRic   = np.zeros((3,size))
        self.dvTotVeh       = np.zeros((size,))
        if (settings["formation"]["orbit"]["elements"] == True):
            self.oeVeh          = np.zeros((6,size))
            self.eeVeh          = np.zeros((6,size))
        if (settings["formation"]["orbit"]["environments"] == True):
            self.eclipseVeh     = np.zeros((size,))
        
        # RSO Paramaters
        self.posRsoEci      = np.zeros((3,size))
        self.velRsoEci      = np.zeros((3,size))
        self.accThrRsoEci   = np.zeros((3,size))
        self.accThrRsoRic   = np.zeros((3,size))
        self.dvTotRso       = np.zeros((size,))
        if (settings["formation"]["orbit"]["elements"] == True):
            self.oeRso          = np.zeros((6,size))
            self.eeRso          = np.zeros((6,size))
        if (settings["formation"]["orbit"]["environments"] == True):
            self.eclipseRso     = np.zeros((size,))
        
        # Relative Parameters
        self.relPosRectRic  = np.zeros((3,size))
        self.relVelRectRic  = np.zeros((3,size))
        self.relPosCurvRic  = np.zeros((3,size))
        self.relVelCurvRic  = np.zeros((3,size))
        self.qInrToRic      = np.zeros((4,size))
        if (settings["formation"]["relStates"] == True):
            self.doe            = np.zeros((6,size))
            self.dee            = np.zeros((6,size))
            self.roe            = np.zeros((6,size))
            self.rectClroe      = np.zeros((6,size))
            self.curvClroe      = np.zeros((6,size))
        if (settings["formation"]["environments"] == True):
            self.losAngEarth = np.zeros((size,))
            self.losAngSun   = np.zeros((size,))
            self.losAngMoon  = np.zeros((size,))
            self.visibility  = np.zeros((size,))
        if (settings["formation"]["measurements"] == True):
            self.measParams = np.zeros((4,size))
            self.qRicToLos  = np.zeros((4,size))
            self.qInrToLos  = np.zeros((4,size))
        
        # Navigation Flight Software
        if (settings["fsw"]["status"] == True):
            # Estimated States
            self.fswNavPosVehEci      = np.zeros((3,size))
            self.fswNavVelVehEci      = np.zeros((3,size))
            self.fswNavOeVeh          = np.zeros((6,size))
            self.fswNavPosRsoEci      = np.zeros((3,size))
            self.fswNavVelRsoEci      = np.zeros((3,size))
            self.fswNavOeRso          = np.zeros((6,size))
            self.fswNavRelPosRectRic  = np.zeros((3,size))
            self.fswNavRelVelRectRic  = np.zeros((3,size))
            self.fswNavRelPosCurvRic  = np.zeros((3,size))
            self.fswNavRelVelCurvRic  = np.zeros((3,size))
            self.fswNavDoe            = np.zeros((6,size))
            self.fswNavDee            = np.zeros((6,size))
            self.fswNavRoe            = np.zeros((6,size))
            self.fswNavRectClroe      = np.zeros((6,size))
            self.fswNavCurvClroe      = np.zeros((6,size))
            # Filter Statuses
            self.fswNavFltrInit       = np.zeros((size,))
            self.fswNavFltrConverged  = np.zeros((size,))
            self.fswNavFltrDiverged   = np.zeros((size,))
            self.fswNavFltrCorrupted  = np.zeros((size,))
            self.fswNavFltrConsistent = np.zeros((size,))
            self.fswNavFltrUnderweight= np.zeros((size,))
            # Filter Covariances
            self.fswNavCovVehEciDiag  = np.zeros((6,size))
            self.fswNavCovVehEciDet   = np.zeros((2,size))
            self.fswNavCovVehEciTrace = np.zeros((2,size))
            self.fswNavCovVehEciFro   = np.zeros((2,size))
            self.fswNavCovVehEciEig   = np.zeros((6,size))
            self.fswNavCovRsoEciDiag  = np.zeros((6,size))
            self.fswNavCovRsoEciDet   = np.zeros((2,size))
            self.fswNavCovRsoEciTrace = np.zeros((2,size))
            self.fswNavCovRsoEciFro   = np.zeros((2,size))
            self.fswNavCovRsoEciEig   = np.zeros((6,size))
            self.fswNavCovRelRicDiag  = np.zeros((6,size))
            self.fswNavCovRelRicDet   = np.zeros((2,size))
            self.fswNavCovRelRicTrace = np.zeros((2,size))
            self.fswNavCovRelRicFro   = np.zeros((2,size))
            self.fswNavCovRelRicEig   = np.zeros((6,size))
            self.fswNavVehSmaVar      = np.zeros((size,))
            self.fswNavRsoSmaVar      = np.zeros((size,))
            # Measurement Residuals
            self.fswNavMeasType       = ['']*size
            self.fswNavMeasRes        = np.zeros((4,size))
            # Errors
            self.fswNavPosVehEciErr      = np.zeros((3,size))
            self.fswNavVelVehEciErr      = np.zeros((3,size))
            self.fswNavPosRsoEciErr      = np.zeros((3,size))
            self.fswNavVelRsoEciErr      = np.zeros((3,size))
            self.fswNavRelPosRectRicErr  = np.zeros((3,size))
            self.fswNavRelVelRectRicErr  = np.zeros((3,size))
            self.fswNavRelPosCurvRicErr  = np.zeros((3,size))
            self.fswNavRelVelCurvRicErr  = np.zeros((3,size))
            self.fswNavDoeErr            = np.zeros((6,size))
            self.fswNavDeeErr            = np.zeros((6,size))
            self.fswNavRoeErr            = np.zeros((6,size))
            self.fswNavRectClroeErr      = np.zeros((6,size))
            self.fswNavCurvClroeErr      = np.zeros((6,size))
            # Filter Consistency Metrics
            self.fswNavMeasNrmInvnSqr   = np.zeros((size,))
            self.fswNavVehNrmEstErrSqr  = np.zeros((size,))
            self.fswNavRsoNrmEstErrSqr  = np.zeros((size,))
            self.fswNavRelNrmEstErrSqr  = np.zeros((size,))
            
        
