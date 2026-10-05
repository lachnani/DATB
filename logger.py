# -*- coding: utf-8 -*-
"""
Created on Sat Mar 22 20:30:54 2025

@author: Hakim Lachnani
"""

import numpy as np
import pandas as pd

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
        self.time = pd.DataFrame(np.zeros((size,2)), columns=['simTime', 'J2000']).T
        
        # Keep-out volume
        self.kovType = "PRISM"
        self.kovLim = 0
        self.kovBreach = False
        self.tKovBreach = 0
        
        # Ephemerides
        if (settings["formation"]["orbit"]["environments"] == True):
            self.sunPosEci = pd.DataFrame(np.zeros((size,3)), columns=['X', 'Y', 'Z']).T
            self.moonPosEci = pd.DataFrame(np.zeros((size,3)), columns=['X', 'Y', 'Z']).T
      
        # Vehicle Parameters
        self.posVehEci      = pd.DataFrame(np.zeros((size,3)), columns=['X', 'Y', 'Z']).T
        self.velVehEci      = pd.DataFrame(np.zeros((size,3)), columns=['X', 'Y', 'Z']).T
        self.accThrVehEci   = pd.DataFrame(np.zeros((size,3)), columns=['X', 'Y', 'Z']).T
        self.accThrVehRic   = pd.DataFrame(np.zeros((size,3)), columns=['R', 'I', 'C']).T
        self.dvTotVeh       = pd.DataFrame(np.zeros((size,1)), columns=['DV']).T
        if (settings["formation"]["orbit"]["elements"] == True):
            self.oeVeh          = pd.DataFrame(np.zeros((size,6)), columns=['a', 'e', 'i', 'RAAN', 'argP', 'M']).T
            self.eeVeh          = pd.DataFrame(np.zeros((size,6)), columns=['a', 'l', 'P1', 'P2', 'Q1', 'Q2']).T
        if (settings["formation"]["orbit"]["environments"] == True):
            self.eclipseVeh     = pd.DataFrame(np.zeros((size,1)), columns=['Eclipse']).T
        
        # RSO Paramaters
        self.posRsoEci      = pd.DataFrame(np.zeros((size,3)), columns=['X', 'Y', 'Z']).T
        self.velRsoEci      = pd.DataFrame(np.zeros((size,3)), columns=['X', 'Y', 'Z']).T
        self.accThrRsoEci   = pd.DataFrame(np.zeros((size,3)), columns=['X', 'Y', 'Z']).T
        self.accThrRsoRic   = pd.DataFrame(np.zeros((size,3)), columns=['R', 'I', 'C']).T
        self.dvTotRso       = pd.DataFrame(np.zeros((size,1)), columns=['DV']).T
        if (settings["formation"]["orbit"]["elements"] == True):
            self.oeRso          = pd.DataFrame(np.zeros((size,6)), columns=['a', 'e', 'i', 'RAAN', 'argP', 'M']).T
            self.eeRso          = pd.DataFrame(np.zeros((size,6)), columns=['a', 'l', 'P1', 'P2', 'Q1', 'Q2']).T
        if (settings["formation"]["orbit"]["environments"] == True):
            self.eclipseRso     = pd.DataFrame(np.zeros((size,1)), columns=['Eclipse']).T
        
        # Relative Parameters
        self.relPosRectRic  = pd.DataFrame(np.zeros((size,3)), columns=['R', 'I', 'C']).T
        self.relVelRectRic  = pd.DataFrame(np.zeros((size,3)), columns=['R', 'I', 'C']).T
        self.relPosCurvRic  = pd.DataFrame(np.zeros((size,3)), columns=['R', 'I', 'C']).T
        self.relVelCurvRic  = pd.DataFrame(np.zeros((size,3)), columns=['R', 'I', 'C']).T
        self.qInrToRic      = pd.DataFrame(np.zeros((size,4)), columns=['qs', 'qx', 'qy', 'qz']).T
        if (settings["formation"]["relStates"] == True):
            self.doe            = pd.DataFrame(np.zeros((size,6)), columns=['da', 'de', 'di', 'dRAAN', 'dargP', 'dM']).T
            self.dee            = pd.DataFrame(np.zeros((size,6)), columns=['da', 'dl', 'dP1', 'dP2', 'dQ1', 'dQ2']).T
            self.roe            = pd.DataFrame(np.zeros((size,6)), columns=['da', 'dlambda', 'dex', 'dey', 'dix', 'diy']).T
            self.rectClroe      = pd.DataFrame(np.zeros((size,6)), columns=['A0', 'alpha', 'xOff', 'yOff', 'B0', 'beta']).T
            self.curvClroe      = pd.DataFrame(np.zeros((size,6)), columns=['A0', 'alpha', 'xOff', 'yOff', 'B0', 'beta']).T
        if (settings["formation"]["environments"] == True):
            self.losAngles  = pd.DataFrame(np.zeros((size,3)), columns=['Earth', 'Sun', 'Moon']).T
            self.visibility = pd.DataFrame(np.zeros((size,1)), columns=['Angle']).T
        if (settings["formation"]["measurements"] == True):
            self.measParams = pd.DataFrame(np.zeros((size,4)), columns=['Rng', 'Rng-Rate', 'Az', 'El']).T
            self.qRicToLos  = pd.DataFrame(np.zeros((size,4)), columns=['qs', 'qx', 'qy', 'qz']).T
            self.qInrToLos  = pd.DataFrame(np.zeros((size,4)), columns=['qs', 'qx', 'qy', 'qz']).T
        
        # Navigation Flight Software
        if (settings["fsw"]["status"] == True):
            # Estimated States
            self.fswNavPosVehEci      = pd.DataFrame(np.zeros((size,3)), columns=['X', 'Y', 'Z']).T
            self.fswNavVelVehEci      = pd.DataFrame(np.zeros((size,3)), columns=['X', 'Y', 'Z']).T
            self.fswNavPosRsoEci      = pd.DataFrame(np.zeros((size,3)), columns=['X', 'Y', 'Z']).T
            self.fswNavVelRsoEci      = pd.DataFrame(np.zeros((size,3)), columns=['X', 'Y', 'Z']).T
            self.fswNavRelPosRectRic  = pd.DataFrame(np.zeros((size,3)), columns=['R', 'I', 'C']).T
            self.fswNavRelVelRectRic  = pd.DataFrame(np.zeros((size,3)), columns=['R', 'I', 'C']).T
            self.fswNavRelPosCurvRic  = pd.DataFrame(np.zeros((size,3)), columns=['R', 'I', 'C']).T
            self.fswNavRelVelCurvRic  = pd.DataFrame(np.zeros((size,3)), columns=['R', 'I', 'C']).T
            self.fswNavDoe            = pd.DataFrame(np.zeros((size,6)), columns=['da', 'de', 'di', 'dRAAN', 'dargP', 'dM']).T
            self.fswNavDee            = pd.DataFrame(np.zeros((size,6)), columns=['da', 'dl', 'dP1', 'dP2', 'dQ1', 'dQ2']).T
            self.fswNavRoe            = pd.DataFrame(np.zeros((size,6)), columns=['da', 'dlambda', 'dex', 'dey', 'dix', 'diy']).T
            self.fswNavRectClroe      = pd.DataFrame(np.zeros((size,6)), columns=['A0', 'alpha', 'xOff', 'yOff', 'B0', 'beta']).T
            self.fswNavCurvClroe      = pd.DataFrame(np.zeros((size,6)), columns=['A0', 'alpha', 'xOff', 'yOff', 'B0', 'beta']).T
            # Filter Statuses
            self.fswNavFltrInit       = pd.DataFrame(np.zeros((size,1)), columns=['Status']).T
            self.fswNavFltrConverged  = pd.DataFrame(np.zeros((size,1)), columns=['Status']).T
            self.fswNavFltrDiverged   = pd.DataFrame(np.zeros((size,1)), columns=['Status']).T
            self.fswNavFltrCorrupted  = pd.DataFrame(np.zeros((size,1)), columns=['Status']).T
            self.fswNavFltrConsistent = pd.DataFrame(np.zeros((size,1)), columns=['Status']).T
            # Filter Covariances
            self.fswNavCovVehEciDiag  = pd.DataFrame(np.zeros((size,6)), columns=['X', 'Y', 'Z', 'Xv', 'Yv', 'Zv']).T
            self.fswNavCovRsoEciDiag  = pd.DataFrame(np.zeros((size,6)), columns=['X', 'Y', 'Z', 'Xv', 'Yv', 'Zv']).T
            self.fswNavCovRelRicDiag  = pd.DataFrame(np.zeros((size,6)), columns=['R', 'I', 'C', 'Rv', 'Iv', 'Cv']).T
            self.fswNavVehSmaVar      = pd.DataFrame(np.zeros((size,1)), columns=['var']).T
            self.fswNavRsoSmaVar      = pd.DataFrame(np.zeros((size,1)), columns=['var']).T
            # Measurement Residuals
            self.fswNavMeasType       = pd.DataFrame(np.zeros((size,1)), columns=['type']).T
            self.fswNavMeasRes        = pd.DataFrame(np.zeros((size,4)), columns=['az', 'el', 'rng', 'rngRate']).T
            # Errors
            self.fswNavPosVehEciErr      = pd.DataFrame(np.zeros((size,3)), columns=['X', 'Y', 'Z']).T
            self.fswNavVelVehEciErr      = pd.DataFrame(np.zeros((size,3)), columns=['X', 'Y', 'Z']).T
            self.fswNavPosRsoEciErr      = pd.DataFrame(np.zeros((size,3)), columns=['X', 'Y', 'Z']).T
            self.fswNavVelRsoEciErr      = pd.DataFrame(np.zeros((size,3)), columns=['X', 'Y', 'Z']).T
            self.fswNavRelPosRectRicErr  = pd.DataFrame(np.zeros((size,3)), columns=['R', 'I', 'C']).T
            self.fswNavRelVelRectRicErr  = pd.DataFrame(np.zeros((size,3)), columns=['R', 'I', 'C']).T
            self.fswNavRelPosCurvRicErr  = pd.DataFrame(np.zeros((size,3)), columns=['R', 'I', 'C']).T
            self.fswNavRelVelCurvRicErr  = pd.DataFrame(np.zeros((size,3)), columns=['R', 'I', 'C']).T
            self.fswNavDoeErr            = pd.DataFrame(np.zeros((size,6)), columns=['da', 'de', 'di', 'dRAAN', 'dargP', 'dM']).T
            self.fswNavDeeErr            = pd.DataFrame(np.zeros((size,6)), columns=['da', 'dl', 'dP1', 'dP2', 'dQ1', 'dQ2']).T
            self.fswNavRoeErr            = pd.DataFrame(np.zeros((size,6)), columns=['da', 'dlambda', 'dex', 'dey', 'dix', 'diy']).T
            self.fswNavRectClroeErr      = pd.DataFrame(np.zeros((size,6)), columns=['A0', 'alpha', 'xOff', 'yOff', 'B0', 'beta']).T
            self.fswNavCurvClroeErr      = pd.DataFrame(np.zeros((size,6)), columns=['A0', 'alpha', 'xOff', 'yOff', 'B0', 'beta']).T
            
        
