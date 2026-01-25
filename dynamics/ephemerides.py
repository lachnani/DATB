"""
Created on 02/08/2025

@author: Hakim Lachnani

ephemerides.py contains:
    - The SunEphemeris Class
    - The MoonEphemeris Class
    
"""

# import required modules:
import numpy as np
from dynamics import dynamicsUtils as uDyn


class SunEphemeris:
    def __new__(cls, *args, **kwargs):
        # Create a new orbit
        return super().__new__(cls)

    def __init__(self, 
              tJ2000):
        # Initialize an ephemeris state
        self.tJ2000 = tJ2000
        
        # Compute position relative to the Earth
        self.rUnit = np.zeros((3,))
        self.r = np.zeros((3,))
        uDyn.sunEph(self.tJ2000, self.rUnit, self.r)
        
    def update(self,tJ2000):
        self.tJ2000 = tJ2000
        
        # Compute position relative to the Earth
        uDyn.sunEph(self.tJ2000, self.rUnit, self.r)
        
        
class MoonEphemeris:
    def __new__(cls, *args, **kwargs):
        # Create a new orbit
        return super().__new__(cls)

    def __init__(self, 
              tJ2000):
        # Initialize an ephemeris state
        self.tJ2000 = tJ2000
        
        # Compute position relative to the Earth
        self.rUnit = np.zeros((3,))
        self.r = np.zeros((3,))
        uDyn.moonEph(self.tJ2000, self.rUnit, self.r)
        
    def update(self,tJ2000):
        self.tJ2000 = tJ2000
        
        # Compute position relative to the Earth
        uDyn.moonEph(self.tJ2000, self.rUnit, self.r)
        