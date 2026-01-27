"""
Created on 02/08/2025

@author: Hakim Lachnani

Orbit.py contains:
    - Parameters for various central bodies
    - The Orbit class
        - oeMeanOscMap: J2 Osculating to/from mean kelplerian elements
    
"""

# import required modules:
import math
import numpy as np
from dynamics import ephemerides as eph
from dynamics import dynamicsUtils as uDyn

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
            uDyn.rv2ee(self.mu, self.r, self.v, self.ee)
            uDyn.rv2oe(self.mu, self.r, self.v, self.oe)
        elif stateType == "STATE_KEPEL":
            self.oe = state
            uDyn.oe2ee(self.oe, self.ee)
            uDyn.oe2rv(self.mu, self.oe, self.r, self.v)
        elif stateType == "STATE_EQEL":
            self.ee = state
            uDyn.ee2rv(self.mu, self.ee, self.r, self.v)
            uDyn.rv2oe(self.mu, self.r, self.v, self.oe)      
            
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
            uDyn.rv2ee(self.mu, self.r, self.v, self.ee)
            uDyn.rv2oe(self.mu, self.r, self.v, self.oe)
            if self.oe[0] > 0:
                self.meanMotion = np.sqrt(self.mu/self.oe[0]**3)
            else:
                self.meanMotion = 0 
        
        self.sun.update(self.tJ2000)
        self.moon.update(self.tJ2000)
        if (self.settings["environments"] == True):
            self.inEclipse = uDyn.eclipse(self.r, self.sun.rUnit)
                               

def oeMeanOscMap(req, J2, oe, oep, sign):
    """
    First-order J2 Mapping Between Mean and Osculating Orbital Elements

    Analytical Mechanics of Space Systems, Hanspeter Schaub, John L. Junkins, 4th edition.
    [m] or [km] should be the same both for req and a

    :param req: equatorial radius
    :param J2:
    :param oe: classical orbit element set
    :param oep:
    :param sign: sgn=1:mean to osc, sgn=-1:osc to mean

    """
    a = oe[0]
    e = oe[1]
    i = oe[2]
    RAAN = oe[3]
    argP = oe[4]
    M = oe[5]
    E = uDyn.M2E(M, e)
    f = uDyn.E2f(E, e)
    gamma2  = sign*J2/2*((req/oe.a)**2)
    eta     = math.sqrt(1-oe.e*oe.e)
    gamma2p = gamma2/(eta**4)
    a_r     = (1+oe.e*math.cos(oe.f))/(eta**2)
    # calculate oep.a
    ap = oe.a + oe.a*gamma2*((3*(math.cos(oe.i))**2-1)*(a_r**3-1/(eta**3)) \
       +3*(1-(math.cos(oe.i))**2)*(a_r**3)*math.cos(2*oe.argP+2*oe.f))  # (F.7)

    de1 = gamma2p/8*e*(eta**2)*(1-11*((math.cos(i))**2)-40*((math.cos(i)) **4) \
        /(1-5*((math.cos(i))**2)))*math.cos(2*argP)  # (F.8)

    de = de1 + (eta ** 2) / 2 * (gamma2 *((3 * ((math.cos(i)) ** 2) - 1) / (eta ** 6) \
       *(e * eta + e / (1 + eta) + 3 * math.cos(f) + 3 * e * ((math.cos(f)) ** 2) + (e ** 2) \
       *((math.cos(f)) ** 3)) + 3 * (1 - ((math.cos(i)) ** 2)) / (eta ** 6) \
       *(e + 3 * math.cos(f) + 3 * e * ((math.cos(f)) ** 2) + (e ** 2) * ((math.cos(f)) ** 3)) * math.cos(2 * argP + 2 * f)) \
       - gamma2p * (1 - ((math.cos(i)) ** 2)) *(3 * math.cos(2 * argP + f) + math.cos(2 * argP + 3 * f)))  # (F.9)

    di = -e*de1/(eta**2)/math.tan(i) + gamma2p/2*math.cos(i)*math.sqrt(1-((math.cos(i))**2)) \
       *(3*math.cos(2*argP+2*f) + 3*e*math.cos(2*argP+f)+e*math.cos(2*argP+3*f))  # (F.10)

    MpopOp = M + argP + RAAN + gamma2p / 8 * (eta ** 3) * (1 - 11 * ((math.cos(i)) ** 2) \
           - 40 *((math.cos(i)) ** 4) / (1 - 5 * ((math.cos(i)) ** 2))) * math.sin(2 * argP) \
           - gamma2p / 16 * (2 + (e ** 2) - 11 * (2 + 3 * (e ** 2)) * ((math.cos(i)) ** 2) - 40 * (2 + 5 * (e ** 2)) \
           *((math.cos(i)) ** 4) / (1 - 5 * ((math.cos(i)) ** 2)) - 400 * (e ** 2) * ((math.cos(i)) ** 6) \
           /((1 - 5 * ((math.cos(i)) ** 2)) ** 2)) * math.sin(2 * argP) \
           + gamma2p / 4 * (-6 *(1 - 5 * ((math.cos(i)) ** 2)) * (f - M + e * math.sin(f)) + (3 - 5 * ((math.cos(i)) ** 2)) \
           *(3 * math.sin(2 * argP + 2 * f) + 3 * e * math.sin(2 * argP + f) + e * math.sin(2 * argP + 3 * f))) \
           - gamma2p / 8 * (e ** 2) * math.cos(i) * (11 + 80 *((math.cos(i)) ** 2) / (1 - 5 * ((math.cos(i)) ** 2)) \
           + 200 * ((math.cos(i)) ** 4) / ((1 - 5 * ((math.cos(i)) ** 2)) ** 2)) * math.sin(2 * argP) \
           - gamma2p / 2 * math.cos(i) * (6 *(f - M + e * math.sin(f)) - 3 * math.sin(2 * argP + 2 * f) \
           - 3 * e * math.sin(2 * argP + f) - e * math.sin(2 * argP + 3 * f))  # (F.11)

    edM = gamma2p / 8 * e * (eta ** 3) * (1 - 11 * ((math.cos(i)) ** 2) - 40 * ((math.cos(i)) ** 4) \
        /(1 - 5 * ((math.cos(i)) ** 2))) * math.sin(2 * argP) - gamma2p / 4 * (eta ** 3) * (2 * (3 * ((math.cos(i)) ** 2) - 1) \
        * ((a_r * eta) ** 2 + a_r + 1) * math.sin(f) + 3 *(1 - ((math.cos(i)) ** 2)) *((-(a_r * eta) ** 2 - a_r + 1) \
        * math.sin(2 * argP + f) + ((a_r * eta) ** 2 + a_r + 1 / 3) * math.sin(2 * argP + 3 * f)))  # (F.12)

    dRAAN = -gamma2p / 8 * (e ** 2) * math.cos(i) * (11 + 80 * ((math.cos(i)) ** 2) /(1 - 5 * ((math.cos(i)) ** 2)) \
           + 200 * ((math.cos(i)) ** 4) /((1 - 5 * ((math.cos(i)) ** 2)) ** 2)) * math.sin(2 * argP) \
           - gamma2p / 2 * math.cos(i) * (6 *(f - M + e * math.sin(f)) - 3 * math.sin(2 * argP + 2 * f) \
           - 3 * e * math.sin(2 * argP + f) - e* math.sin(2 * argP + 3 * f))  # (F.13)

    d1 = (e+de)*math.sin(M) + edM*math.cos(M)  # (F.14)
    d2 = (e+de)*math.cos(M) - edM*math.sin(M)  # (F.15)

    Mp = math.atan2(d1, d2)  # (F.16)
    ep = math.sqrt(d1**2+d2**2)  # (F.17)

    d3 = (math.sin(i/2)+math.cos(i/2)*di/2)*math.sin(RAAN) + math.sin(i/2)*dRAAN*math.cos(RAAN)  # (F.18)
    d4 = (math.sin(i/2)+math.cos(i/2)*di/2)*math.cos(RAAN) - math.sin(i/2)*dRAAN*math.sin(RAAN)  # (F.19)

    RAANp = math.atan2(d3, d4)  # (F.20)
    ip = 2*math.asin(math.sqrt(d3**2+d4**2))  # (F.21)
    argPp = MpopOp - Mp - RAANp  # (F.22)

    Ep = uDyn.M2E(Mp, ep)
    fp = uDyn.E2f(Ep, ep)

    oep.a = ap
    oep.e = ep
    oep.i = ip
    oep.RAAN = RAANp
    oep.argP = argPp
    oep.f = fp
    return
