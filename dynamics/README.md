# Dynamics

The DATB dynamics module contains models and utilities for orbital mechanics and translational aspects of relative motion.

There are four dynamics classes in DATB:

* Moon Ephemeris
* Sun Ephemeris
* Orbit
* Formation

## Ephemerides

There are classes for moon and sun ephemerides contained in ephemerides.py. Each classes contains three attributes:

* tJ2000: epoch time since J2000 in seconds
* r: position in ECI
* rUnit: position direction in ECI

The update(tJ2000) method updates the time and position of the ephemeris.

## Orbit

The orbit class is contained in orbit.py. The orbit class contains the following attributes:

* t: time since simulation start
* tJ2000: epoch time since J2000 in seconds
* r: position in ECI
* v: velocity in ECI
* oe: Keplerian orbit elements
* ee: Equinoctial elements
* mu: Earth gravitational acceleration (nominally constant)
* meanMotion: orbit mean motion
* moon: moon ephemeris class
* sun: sun ephemeris class
* aCtrlInEci: Control acceleration in ECI
* dvTot: accumulated delta-v over scenario
* inEclipse: eclipse status
* Cd: coefficient of drag
* normA: normalized area (area/mass)
* settings: simulation settings
* pert: pertubation settings

The propagate(dt) method uses and RK4 routine to propagate the orbit and ephemerides.

## Formation

The formation class is contained in formation.py. The formation class contains the following attributes:

* t: time since simulation start
* tJ2000: epoch time since J2000 in seconds
* chief: chief orbit class
* deputy: deputy orbit class
* relPosRectRic: relative position in the rectilinear RIC frame
* relVelRectRic: relative velocity in the rectilinear RIC frame
* relPosCurvRic: relative position in the curvilinear RIC frame
* relVelCurvRic: relative velocity in the curvilinear RIC frame
* doe: differential Keplerian orbit elements
* dee: differential Equinoctial elements
* dAmico: d'Amico relative orbital elements
* rectClroe: rectilinear CLROEs
* curvClroe: curvilinear CLROEs
* dcmInr2Ric: ECI to RIC frame DCM
* dcmRic2Los: RIC to Line of Sight frame DCM
* settings: simulation settings
* losEarthAng: angle between LOS vector and Earth center
* losMoonAng: angle between LOS vector and Moon center
* losSunAng: angle between LOS vector and Sun center
* rng: range
* rngRate: range rate
* az: azimuth in RIC frame
* el: elevation in RIC frame

The propagate(dt) method propagates both the chief and deputy orbits and recalculates relative parameters.



