"""
Created on 02/08/2025

@author: Hakim Lachnani

Orbit.py contains:
    - Parameters for various central bodies
    - The Orbit class
    - Orbit utilities:
        - E2f: Eccentric to true anomaly
        - E2M: Eccentric to mean anomaly
        - f2E: True to eccentric anomaly
        - M2E: Mean to eccentric anomaly
        - oe2rv: Keplerian elements to inertial states
        - rv2oe: Inertial states to keplerian elements
        - oe2ee: Keplerian elements to equinoctial elements
        - rv2ee: Inertial states to equinoctial elements
        - ee2rv: Equinoctial elements to inertial states
        - grav: Gravitational acceleration
        - atmosphericDensity
        - atmosphericDrag
        - jPerturb
        - solarRad: solar radiation acceleration
        - oeMeanOscMap: J2 Osculating to/from mean kelplerian elements
    
"""

# import required modules:
import math
import numpy as np
from numpy import linalg as la
from dynamics import ephemerides as eph
from dynamics import dynamicsUtils as uDyn

DB0_EPS = 1e-30
eps = 1e-13
maxIteration = 200
tolerance = 1e-15

AU = 149597870.693  # astronomical unit in units of kilometers #
D2R = (np.pi / 180.)

# Gravitational Constants mu = G*m, where m is the body of the attracting body.  All units are km^3/s^2.
# Values are obtained from SPICE kernels in http://naif.jpl.nasa.gov/pub/naif/generic_kernels/
MU_SUN = 132712440023.310
MU_EARTH = 398600.436
MU_MOON = 4902.799

# body information for major solar system bodies. Units are in km.
# data taken from http://nssdc.gsfc.nasa.gov/planetary/planets.html
# Sun #
REQ_SUN = 695000.  # km #
SUN_ANGLE_FROM_EARTH = 0.00436332312998582

# Earth #
REQ_EARTH = 6378.1366  # km, from SPICE #
RP_EARTH = 6356.7519  # km, from SPICE #
J2_EARTH = 1082.616e-6
J3_EARTH = -2.53881e-6
J4_EARTH = -1.65597e-6
J5_EARTH = -0.15e-6
J6_EARTH = 0.57e-6
SMA_EARTH = 1.00000011 * AU
I_EARTH = 0.00005 * D2R
E_EARTH = 0.01671022
OMEGA_EARTH = 0.00007292115  # Earth's planetary rotation rate, rad/sec #

# Moon #
REQ_MOON = 1737.4
J2_MOON = 202.7e-6
SMA_MOON = 0.3844e6
I_MOON = 5.154 * D2R 
E_MOON = 0.0549


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
        self.aCtrlInRic = np.zeros(3)
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
                               


""" 
Define all the utility functions:
    Accelerations:
        grav
        
    
"""

def atmosphericDensity(alt):
    """
    This program computes the atmospheric density based on altitude
    supplied by user.  This function uses a curve fit based on
    atmospheric data from the Standard Atmosphere 1976 Data. This
    function is valid for altitudes ranging from 100km to 1000km.

    .. note::

        This code can only be applied to spacecraft orbiting the Earth

    :param alt: altitude in km
    :return:  density at the given altitude in kg/m^3
    """
    # Smooth exponential drop-off after 1000 km #
    if alt > 1000.:
        logdensity = (-7E-05) * alt - 14.464
        density = math.pow(10., logdensity)
        return density

    # Calculating the density based on a scaled 6th order polynomial fit to the log of density #
    val = (alt - 526.8000) / 292.8563
    logdensity = 0.34047 * math.pow(val, 6) - 0.5889 * math.pow(val, 5) - 0.5269 * math.pow(val, 4) \
                 + 1.0036 * math.pow(val, 3) + 0.60713 * math.pow(val, 2) - 2.3024 * val - 12.575

    # Calculating density by raising 10 to the log of density #
    density = math.pow(10., logdensity)

    return density


def atmosphericDrag(Cd, A, m, r, v):
    """
     This program computes the atmospheric drag acceleration
     vector acting on a spacecraft.
     Note the acceleration vector output is inertial, and is
     only valid for altitudes up to 1000 km.
     Afterwards the drag force is zero. Only valid for Earth.

     :param Cd:  drag coefficient of the spacecraft
     :param A: cross-sectional area of the spacecraft in m^2
     :param m: mass of the spacecraft in kg
     :param r: Inertial position vector of the spacecraft in km  [x;y;z]
     :param v: Inertial velocity vector of the spacecraft in km/s [vx;vy;vz]
     :return: The inertial acceleration vector due to atmospheric drag in km/sec^2
    """
    # find the altitude and velocity #
    rMag = la.norm(r)
    vMag = la.norm(v)
    alt = rMag - REQ_EARTH
    advec = np.zeros(3)

    # Checking if user supplied a orbital position is insede the earth #
    if alt <= 0.:
        print("ERROR: atmosphericDrag() received r = [{} {} {}].". \
            format(str(r[1]), str(r[2]), str(r[3])))
        print('The value of r should produce a positive altitude for the Earth.')
        advec.fill(np.NaN)
        return

    # get the Atmospheric density at the given altitude in kg/m^3 #
    density = atmosphericDensity(alt)

    # compute the magnitude of the drag acceleration #
    ad = ((-0.5) * density * (Cd * A / m) * (math.pow(vMag * 1000., 2))) / 1000.

    # computing the vector for drag acceleration #
    advec = (ad / vMag) * v

    return advec


def solarRad(A, m, sunvec):
    """
    Computes the inertial solar radiation force vectors
    based on cross-sectional Area and mass of the spacecraft
    and the position vector of the body to the sun.

    .. note::

        It is assumed that the solar radiation pressure decreases quadratically with distance from sun (in AU)

    Solar Radiation Equations obtained from
    Earth Space and Planets Journal Vol. 51, 1999 pp. 979-986

    :param A: Cross-sectional area of the spacecraft that is facing the sun in m^2.
    :param m: The mass of the spacecraft in kg.
    :param sunvec: Position vector to the Sun in units of AU. Earth has a distance of 1 AU.
    :return:   arvec, The inertial acceleration vector due to the effects of Solar Radiation pressure in km/sec^2.  The vector
               components of the output are the same as the vector
               components of the sunvec input vector.
    """
    # Solar Radiation Flux #
    flux = 1372.5398

    # Speed of light #
    c = 299792458.

    # Radiation pressure coefficient #
    Cr = 1.3

    # Magnitude of position vector #
    sundist = la.norm(sunvec)

    # Computing the acceleration vector #
    arvec = ((-Cr * A * flux) / (m * c * math.pow(sundist, 3)) / 1000.) * sunvec

    return arvec


def v3Normalize(v):
    result = np.zeros(3)
    norm = la.norm(v)
    if norm > DB0_EPS:
        result = (1. / norm) * v
    return result


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
