// file: dynamicsUtils.c
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <stdbool.h>
#include <stddef.h>

#include "dynamicsConstants.h"
#include "dynamicsUtils.h"

// Helper vector math functions
static double clamp(double x, double min, double max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}
static void v3_zero(double v[3]) {
    v[0] = v[1] = v[2] = 0.0;
}
static void v3_set(const double x, const double y, const double z, double out[3]) {
    out[0] = x; out[1] = y; out[2] = z;
}
static void v3_copy(const double src[3], double dest[3]) {
    dest[0] = src[0]; dest[1] = src[1]; dest[2] = src[2];
}
static double v3_dot(const double a[3], const double b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
static void v3_cross(const double a[3], const double b[3], double out[3]) {
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}
static double v3_norm(const double v[3]) {
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}
static void v3_scale(const double v[3], double s, double out[3]) {
    out[0] = v[0] * s; out[1] = v[1] * s; out[2] = v[2] * s;
}
static void v3_sub(const double a[3], const double b[3], double out[3]) {
    out[0] = a[0] - b[0]; out[1] = a[1] - b[1]; out[2] = a[2] - b[2];
}
static void v3_add(const double a[3], const double b[3], double out[3]) {
    out[0] = a[0] + b[0]; out[1] = a[1] + b[1]; out[2] = a[2] + b[2];
}
static void v3_add_inplace(double a[3], const double b[3]) {
    a[0] += b[0]; a[1] += b[1]; a[2] += b[2];
}
static void v3_normalize(const double v[3], double out[3]) {
    double n = v3_norm(v);
    if (n > 1e-30) {
        out[0] = v[0] / n; out[1] = v[1] / n; out[2] = v[2] / n;
    }
    else {
        out[0] = out[1] = out[2] = 0.0;
    }
}
static double v3_angle(const double a[3], const double b[3]) {
    // Returns the angle in radians between vectors a and b
    double dot = v3_dot(a, b);
    double norm_a = v3_norm(a);
    double norm_b = v3_norm(b);
    if (norm_a < 1e-30 || norm_b < 1e-30) {
        // Undefined angle if either vector is zero
        return 0.0;
    }
    double cos_theta = dot / (norm_a * norm_b);
    // Clamp to [-1, 1] to avoid NaN due to floating point errors
    cos_theta = clamp(cos_theta, -1.0, 1.0);
    return acos(cos_theta);
}
static double wrapTo2Pi(double x) {
    x = fmod(x, 2.0 * M_PI);
    if (x < 0) x += 2.0 * M_PI;
    return x;
}
static void mat3_v3_mul(const double mat[3][3], const double vec[3], double out[3]) {
    for (int i = 0; i < 3; ++i) {
        out[i] = mat[i][0] * vec[0] + mat[i][1] * vec[1] + mat[i][2] * vec[2];
    }
}
static void mat3_transpose(const double M[3][3], double MT[3][3]) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            MT[i][j] = M[j][i];
}

static double atmosphericDensity(double alt)
{
    /*
	Computes atmospheric density at a given altitude using a curve fit to the Standard Atmosphere 1976 data.

	Parameters
    ----------
    alt : double
		Altitude in kilometers (valid range: 100 km to 1000 km)

    Returns
    -------
	density : double
	    Atmospheric density at the given altitude in kg/m^3
    */

    double logdensity;
    double density;
    double val;

    /* Smooth exponential drop-off after 1000 km */
    if (alt > 1000.) {
        logdensity = (-7e-05) * alt - 14.464;
        density = pow(10., logdensity);
        return density;
    }

    /* Calculating the density based on a scaled 6th order polynomial fit to the log of density */
    val = (alt - 526.8000) / 292.8563;
    logdensity = 0.34047 * pow(val, 6) - 0.5889 * pow(val, 5) - 0.5269 * pow(val, 4)
        + 1.0036 * pow(val, 3) + 0.60713 * pow(val, 2) - 2.3024 * val - 12.575;

    /* Calculating density by raising 10 to the log of density */
    density = pow(10., logdensity);

    return density;
}

/*
* Orbit Functions
*/

void acceleration(
    const int solarGrav,
    const int lunarGrav,
    const int drag,
	const int jnum,
    const double moon_r[3], 
    const double sun_r[3], 
	const double Cd,
	const double normA,
    const double r[3], 
	const double v[3],
	const double aCtrlInEci[3],
    double a[3]) {
    /*
    Orbit acceleration due to gravity, perturbations, and thrust

    Parameters
    ----------
    solarGrav: int
        Flag to include solar gravity perturbation (1 to include, 0 to exclude)
    lunarGrav: int
        Flag to include lunar gravity perturbation (1 to include, 0 to exclude)
    drag: int
        Flag to include atmospheric drag perturbation (1 to include, 0 to exclude)
    jnum: int
        Number of zonal harmonics to include in Earth's gravity model (0 for point mass, 2 for J2, etc.)
    moon_r: Array
        Position vector of the Moon in kilometers [x;y;z]. Required if lunarGrav is 1.
    sun_r: Array
        Position vector of the Sun in kilometers [x;y;z]. Required if solarGrav is 1.
    Cd: double
        Drag coefficient for atmospheric drag calculation. Required if drag is 1.
    normA: double
        Spacecraft cross-sectional area divided by mass (A/m) for atmospheric drag calculation. Required if drag is 1.
    r: Array
		Position vector in kilometers [x;y;z] at the current time step.
    v: Array
		Velocity vector in kilometers per second [vx;vy;vz] at the current time step.
    aCtrlInEci: Array
        Control acceleration vector in ECI frame in km/s^2. Can be zero if no control acceleration.

    Returns
    -------
	a: Array
		Total acceleration vector in kilometers per second squared [ax;ay;az] at the current time step, 
        including gravity, perturbations, and control acceleration.
    */

    double aEarth[3], aMoon[3], aSun[3], aDrag[3], aCtrl[3], tmp[3];
	double earth_r[3] = { 0, 0, 0 }; // Earth at origin

    grav(r, CELESTIAL_EARTH, earth_r, aEarth);
    if (jnum > 1) {
        jPerturb(r, jnum, CELESTIAL_EARTH, tmp);
        v3_add_inplace(aEarth, tmp);
    }

    if (lunarGrav) {
        double rMoon[3];
        v3_sub(r, moon_r, rMoon);
        grav(rMoon, CELESTIAL_MOON, moon_r, aMoon);
        jPerturb(rMoon, 2, CELESTIAL_MOON, tmp);
        v3_add_inplace(aMoon, tmp);
    }
    else {
        v3_zero(aMoon);
    }

    if (solarGrav) {
        double rSun[3];
        v3_sub(r, sun_r, rSun);
        grav(rSun, CELESTIAL_SUN, sun_r, aSun);
    }
    else {
        v3_zero(aSun);
    }

    if (drag) {
        dragPerturb(Cd, normA, r, v, aDrag);
    }
    else {
        v3_zero(aDrag);
	}

    v3_copy(aCtrlInEci, aCtrl);

    // a = aEarth + aMoon + aSun + aDrag + aCtrl
    a[0] = aEarth[0] + (aMoon[0] + aSun[0] + aDrag[0] + aCtrl[0]);
    a[1] = aEarth[1] + (aMoon[1] + aSun[1] + aDrag[1] + aCtrl[1]);
    a[2] = aEarth[2] + (aMoon[2] + aSun[2] + aDrag[2] + aCtrl[2]);
}

void stateUpdate(
    const int solarGrav,
    const int lunarGrav,
    const int drag,
    const int jnum,
    const double moon_r[3],
    const double sun_r[3],
    const double Cd,
    const double normA,
    const double aCtrlInEci[3], 
    const double y[6], 
    double yDot[6]) {
    /*
	State update equations for RK4 integration. 
    Computes the time derivative of the state vector y = [r(3), v(3)] given the current state and acceleration.

    Parameters
    ----------
    solarGrav: int
        Flag to include solar gravity perturbation (1 to include, 0 to exclude)
    lunarGrav: int
        Flag to include lunar gravity perturbation (1 to include, 0 to exclude)
    drag: int
        Flag to include atmospheric drag perturbation (1 to include, 0 to exclude)
    jnum: int
        Number of zonal harmonics to include in Earth's gravity model (0 for point mass, 2 for J2, etc.)
    moon_r: Array
        Position vector of the Moon in kilometers [x;y;z]. Required if lunarGrav is 1.
    sun_r: Array
        Position vector of the Sun in kilometers [x;y;z]. Required if solarGrav is 1.
    Cd: double
        Drag coefficient for atmospheric drag calculation. Required if drag is 1.
    normA: double
        Spacecraft cross-sectional area divided by mass (A/m) for atmospheric drag calculation. Required if drag is 1.
    aCtrlInEci: Array
        Control acceleration vector in ECI frame in km/s^2. Can be zero if no control acceleration.
    y: Array
		State vector [r(3), v(3)] at the current time step.

    Returns
    -------
	yDot: Array
		State derivative vector [rDot(3), vDot(3)] where rDot = v and vDot = acceleration(r) at the current time step.
    */

    /* State update: y = [r(3), v(3)], yDot = [v, a] */

    double r[3], v[3], vDot[3];
    v3_copy(y, r);
    v3_copy(y + 3, v);

    // rDot = v
    v3_copy(v, yDot);

    // vDot = acceleration(r)
	acceleration(solarGrav, lunarGrav, drag, jnum, moon_r, sun_r, Cd, normA, r, v, aCtrlInEci, vDot);
    v3_copy(vDot, yDot + 3);
}

void Orbit_rk4(
    const int solarGrav,
    const int lunarGrav,
    const int drag,
    const int jnum,
    const double moon_r[3],
    const double sun_r[3],
    const double Cd,
    const double normA,
    const double aCtrlInEci[3],
    double dt,
    double r[3],
    double v[3]) {
    /*
	Propagates orbit using 4th order Runge-Kutta numerical integration method with the state update function defined above.

    Parameters
    ----------
    solarGrav: int
		Flag to include solar gravity perturbation (1 to include, 0 to exclude)
	lunarGrav: int
		Flag to include lunar gravity perturbation (1 to include, 0 to exclude)
	drag: int
	    Flag to include atmospheric drag perturbation (1 to include, 0 to exclude)
	jnum: int
	    Number of zonal harmonics to include in Earth's gravity model (0 for point mass, 2 for J2, etc.)
	moon_r: Array
	    Position vector of the Moon in kilometers [x;y;z]. Required if lunarGrav is 1.
	sun_r: Array
	    Position vector of the Sun in kilometers [x;y;z]. Required if solarGrav is 1.
	Cd: double
	    Drag coefficient for atmospheric drag calculation. Required if drag is 1.
	normA: double
	    Spacecraft cross-sectional area divided by mass (A/m) for atmospheric drag calculation. Required if drag is 1.
	aCtrlInEci: Array
	    Control acceleration vector in ECI frame in km/s^2. Can be zero if no control acceleration.
    dt: double
		Time step for propagation in seconds.

    Returns
    -------
    r: Array
		Updated position vector in kilometers [x;y;z] after propagation.
	v: Array
		Updated velocity vector in kilometers per second [vx;vy;vz] after propagation.
    */

    double y0[6], k1[6], k2[6], k3[6], k4[6], y[6], tmp[6];
    // y0 = [r, v]
    v3_copy(r, y0);
    v3_copy(v, y0 + 3);

    // k1 = f(y0)
    stateUpdate(solarGrav, lunarGrav, drag, jnum, moon_r, sun_r, Cd, normA, aCtrlInEci, y0, k1);

    // k2 = f(y0 + dt*k1/2)
    for (int i = 0; i < 6; ++i) tmp[i] = y0[i] + dt * k1[i] / 2.0;
	stateUpdate(solarGrav, lunarGrav, drag, jnum, moon_r, sun_r, Cd, normA, aCtrlInEci, tmp, k2);

    // k3 = f(y0 + dt*k2/2)
    for (int i = 0; i < 6; ++i) tmp[i] = y0[i] + dt * k2[i] / 2.0;
    stateUpdate(solarGrav, lunarGrav, drag, jnum, moon_r, sun_r, Cd, normA, aCtrlInEci, tmp, k3);

    // k4 = f(y0 + dt*k3)
    for (int i = 0; i < 6; ++i) tmp[i] = y0[i] + dt * k3[i];
    stateUpdate(solarGrav, lunarGrav, drag, jnum, moon_r, sun_r, Cd, normA, aCtrlInEci, tmp, k4);

    // y = y0 + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
    for (int i = 0; i < 6; ++i)
        y[i] = y0[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);

    // Update orbit state
    v3_copy(y, r);
    v3_copy(y + 3, v);
    
}


int eclipse(const double rSc[3], const double sun_rUnit[3]) {
    /*
	Determines whether Earth-orbiting spacecraft is in eclipse by the Sun using geometric calculations 
    based on the apparent size of the Earth and Sun from the spacecraft's perspective.

    Parameters
    ----------
    rSc: Array
		Spacecraft position vector in kilometers [x;y;z].
    sun_rUnit: Array
		Unit vector from Earth to Sun in the inertial frame [x;y;z]

    Returns
    -------
    eclipseFlag: int
		Flag for spacecraft eclipse status (1 if in eclipse, 0 if not in eclipse)
    */

    double rScMag = v3_norm(rSc);

    // rScUnit = rSc / rScMag
    double rScUnit[3] = { 0, 0, 0 };
    if (rScMag > 0) {
        rScUnit[0] = rSc[0] / rScMag;
        rScUnit[1] = rSc[1] / rScMag;
        rScUnit[2] = rSc[2] / rScMag;

        // Angle between Sun and Earth from SC
        double sunEarthAng = acos(
            clamp(-1.0 * v3_dot(sun_rUnit, rScUnit), -1.0, 1.0)
        );

        // Geometric calculations
        double cosg = sqrt(rScMag * rScMag - REQ_EARTH * REQ_EARTH) / rScMag;
        double acosg = acos(clamp(cosg, -1.0, 1.0));
        double apparentEarthRadius = sin(acosg);
        double cosRatio = 0.0;
        double apparentSunRadius = 0.0;
        double apparentEarthSunDist = 0.0;

        cosRatio = cosg / COS_SUN_ANGLE_FROM_EARTH;
        apparentSunRadius = cosRatio * SIN_SUN_ANGLE_FROM_EARTH;
        apparentEarthSunDist = sqrt(cosRatio * cosRatio + 1.0 - 2.0 * cosRatio * cos(sunEarthAng));

        if (apparentEarthSunDist > (apparentEarthRadius + apparentSunRadius)) {
            return 0;
        }
        else {
            return 1;
        }
    }
    else {
        return 1;
    }
}

void grav(const double r[3], const int body, const double rBody[3], double agrav[3]) {
    /*
    Computes the gravitational acceleration

    Parameters
    ----------
	r: Array
	    Cartesian Position vector in kilometers [x;y;z].
	body: int
	    Celestial body for which to calculate gravity (CELESTIAL_EARTH, CELESTIAL_MOON, CELESTIAL_SUN)
	rBody: Array
	    Position vector of the body exerting gravity in kilometers [x;y;z]

    Returns
    -------
	agrav: Array
	    The total acceleration vector due to the gravity in km/sec^2 [accelx;accely;accelz]
    */
    double mu = 0.0;
    switch (body)
    {
    case CELESTIAL_EARTH:
        mu = MU_EARTH;
		break;
    case CELESTIAL_MOON:
		mu = MU_MOON;
		break;
	case CELESTIAL_SUN:
		mu = MU_SUN;
		break;
    default:
        mu = MU_EARTH;
        break;
    }

    double rMag = v3_norm(r);
	double aGravBody[3] = { 0, 0, 0 };
    if (rMag > 0.0) {
		aGravBody[0] = -mu * r[0] / (rMag * rMag * rMag);
		aGravBody[1] = -mu * r[1] / (rMag * rMag * rMag);
		aGravBody[2] = -mu * r[2] / (rMag * rMag * rMag);
    }

	double rBodyMag = v3_norm(rBody);
	double aGravEarthToBody[3] = { 0, 0, 0 };
	if ((rBodyMag > 0.0) && (body != CELESTIAL_EARTH)) {
		aGravEarthToBody[0] = -mu * rBody[0] / (rBodyMag * rBodyMag * rBodyMag);
		aGravEarthToBody[1] = -mu * rBody[1] / (rBodyMag * rBodyMag * rBodyMag);
		aGravEarthToBody[2] = -mu * rBody[2] / (rBodyMag * rBodyMag * rBodyMag);
	}

    v3_add(aGravBody, aGravEarthToBody, agrav);
}

void jPerturb(const double r[3], int num, const int body, double ajtot[3])
{
    /*
	J-perturbation acceleration for Earth and Moon. Sun is not oblate, so no J perturbations are applied.

    Parameters
    ----------
    r: Array
		Cartesian Position vector in kilometers [x;y;z].
    num: int
		Number of zonal harmonics to include (2-6)
	body: int
	    Celestial body for which to calculate J perturbations (CELESTIAL_EARTH, CELESTIAL_MOON, CELESTIAL_SUN)

    Returns
    -------
    ajtot: Array
		Total acceleration vector due to J perturbations in km/sec^2 [accelx;accely;accelz]
    */

    double mu;
    double req;
    double J2, J3, J4, J5, J6;
    double x;
    double y;
    double z;
    double rMag;
    double temp[3];
    double temp2[3];

    v3_zero(ajtot);

    switch (body)
    {
    case CELESTIAL_EARTH:
        mu = MU_EARTH;
        req = REQ_EARTH;
        J2 = J2_EARTH;
        J3 = J3_EARTH;
        J4 = J4_EARTH;
        J5 = J5_EARTH;
        J6 = J6_EARTH;
        break;
    case CELESTIAL_MOON:
        mu = MU_MOON;
        req = REQ_MOON;
        J2 = J2_MOON;
        J3 = 0.0;
        J4 = 0.0;
        J5 = 0.0;
        J6 = 0.0;
        break;
    case CELESTIAL_SUN:
        return;
    default:
        mu = MU_EARTH;
        req = REQ_EARTH;
        J2 = J2_EARTH;
        J3 = J3_EARTH;
        J4 = J4_EARTH;
        J5 = J5_EARTH;
        J6 = J6_EARTH;
        break;
    }

    /* Calculate the J perturbations */
    x = r[0];
    y = r[1];
    z = r[2];
    rMag = v3_norm(r);

    /* Calculating the total acceleration based on user input */
    if (num >= 2) {
        v3_set((1.0 - 5.0 * pow(z / rMag, 2.0)) * (x / rMag),
            (1.0 - 5.0 * pow(z / rMag, 2.0)) * (y / rMag),
            (3.0 - 5.0 * pow(z / rMag, 2.0)) * (z / rMag), ajtot);
        v3_scale(ajtot, -3.0 / 2.0 * J2 * (mu / pow(rMag, 2.0)) * pow(req / rMag, 2.0), ajtot);
    }
    if (num >= 3) {
        v3_set(5.0 * (7.0 * pow(z / rMag, 3.0) - 3.0 * (z / rMag)) * (x / rMag),
            5.0 * (7.0 * pow(z / rMag, 3.0) - 3.0 * (z / rMag)) * (y / rMag),
            -3.0 * (10.0 * pow(z / rMag, 2.0) - (35.0 / 3.0) * pow(z / rMag, 4.0) - 1.0), temp);
        v3_scale(temp, 1.0 / 2.0 * J3 * (mu / pow(rMag, 2.0)) * pow(req / rMag, 3.0), temp2);
        v3_add(ajtot, temp2, ajtot);
    }
    if (num >= 4) {
        v3_set((3.0 - 42.0 * pow(z / rMag, 2.0) + 63.0 * pow(z / rMag, 4.0)) * (x / rMag),
            (3.0 - 42.0 * pow(z / rMag, 2.0) + 63.0 * pow(z / rMag, 4.0)) * (y / rMag),
            (15.0 - 70.0 * pow(z / rMag, 2) + 63.0 * pow(z / rMag, 4.0)) * (z / rMag), temp);
        v3_scale(temp, 5.0 / 8.0 * J4 * (mu / pow(rMag, 2.0)) * pow(req / rMag, 4.0), temp2);
        v3_add(ajtot, temp2, ajtot);
    }
    if (num >= 5) {
        v3_set(3.0 * (35.0 * (z / rMag) - 210.0 * pow(z / rMag, 3.0) + 231.0 * pow(z / rMag, 5.0)) * (x / rMag),
            3.0 * (35.0 * (z / rMag) - 210.0 * pow(z / rMag, 3.0) + 231.0 * pow(z / rMag, 5.0)) * (y / rMag),
            -(15.0 - 315.0 * pow(z / rMag, 2.0) + 945.0 * pow(z / rMag, 4.0) - 693.0 * pow(z / rMag, 6.0)), temp);
        v3_scale(temp, 1.0 / 8.0 * J5 * (mu / pow(rMag, 2.0)) * pow(req / rMag, 5.0), temp2);
        v3_add(ajtot, temp2, ajtot);
    }
    if (num >= 6) {
        v3_set((35.0 - 945.0 * pow(z / rMag, 2) + 3465.0 * pow(z / rMag, 4.0) - 3003.0 * pow(z / rMag, 6.0)) * (x / rMag),
            (35.0 - 945.0 * pow(z / rMag, 2.0) + 3465.0 * pow(z / rMag, 4.0) - 3003.0 * pow(z / rMag, 6.0)) * (y / rMag),
            -(3003.0 * pow(z / rMag, 6.0) - 4851.0 * pow(z / rMag, 4.0) + 2205.0 * pow(z / rMag, 2.0) - 245.0) * (z / rMag), temp);
        v3_scale(temp, -1.0 / 16.0 * J6 * (mu / pow(rMag, 2.0)) * pow(req / rMag, 6.0), temp2);
        v3_add(ajtot, temp2, ajtot);
    }
}

void dragPerturb(const double Cd, const double normA, const double r[3], const double v[3], double ad[3]) {

    /*
    Drag perturbation acceleration

    Parameters
    ----------
    Cd : double
        drag coefficient
    normA : double
        cross-sectional area of the spacecraft in m^2 divided by the mass of the spacecraft in kg
    r : Array
        Inertial position vector of the spacecraft in km  [x;y;z]
    v : Array
        Inertial velocity vector of the spacecraft in km/s [vx;vy;vz]

    Returns
    -------
    ad : Array
        Inertial acceleration vector due to atmospheric drag in km/sec^2
    */

    double rMag;
    double vMag;
    double alt;
    double density;
	double adMag;

    /* find the altitude and velocity */
    rMag = v3_norm(r);
    vMag = v3_norm(r);
    alt = rMag - REQ_EARTH;

    /* Checking if user supplied a orbital position is inside the earth or if the altitude is above 1000 km */
    if ((alt <= 0.) || (alt > 1000)) {
        v3_set(0.0, 0.0, 0.0, ad);
        return;
    }

    /* get the Atmospheric density at the given altitude in kg/m^3 */
    density = atmosphericDensity(alt);

    /* compute the magnitude of the drag acceleration */
    adMag = ((-0.5) * density * (Cd * normA) * (pow(vMag * 1000., 2))) / 1000.;

    /* computing the vector for drag acceleration */
    v3_scale(v, adMag / vMag, ad);
}

/*
* Ephemeris Functions
*/

void moonEph(double tJ2000, double rUnit[3], double r[3])
{
    /*
    Moon Ephemeris relative to Earth

    Ref: David Simpson, "An Alternative Lunar Ephemeris Model for On-board Flight Software Use"

    Parameters
	----------
	tJ2000 : double
	    Time in seconds since J2000 epoch (January 1, 2000, 12:00 TT)

	Returns
	-------
    rUnit : Array
		Unit vector from Earth to Moon in ECI frame [x; y; z]
	r : Array
	    Vector from Earth to Moon in ECI frame in km [x; y; z]
    */

    double T = 0.0;
    double eclipticLongitude = 0.0;
    double eclipticLongitudeJ2K = 0.0;
    double eclipticLattitude = 0.0;
    double eclipticLattitudeJ2K = 0.0;
    double horizontalParallax = 0.0;
    double rMagEarthRad = 0.0;
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;

    const double moon_a[3][7] = { {
        383.0, 31.5, 10.6, 6.2, 3.2, 2.3, 0.8},{
        351.0, 28.9, 13.7, 9.7, 5.7, 2.9, 2.1},{
        153.2, 31.5, 12.5, 4.2, 2.5, 3.0, 1.8} };

    const double moon_w[3][7] = { {
        8399.685, 70.990, 16728.377, 1185.622, 7143.070, 15613.745, 8467.263},{
        8399.687, 70.997, 8433.466, 16728.380, 1185.667, 7143.058, 15613.755},{
        8399.672, 8433.464, 70.996, 16728.364, 1185.645, 104.881, 9399.116} };

    const double moon_delta[3][7] = { {
        5.381, 6.169, 1.453, 0.481, 5.017, 0.857, 1.010},{
        3.811, 4.596, 4.766, 6.165, 5.164, 0.300, 5.565},{
        3.807, 1.629, 4.595, 6.162, 5.167, 2.555, 6.248} };
    
    // Centuries since J2000 
    T = SEC2CENTURY * tJ2000;
    
    // Ecliptic Longitude in radians
    eclipticLongitude = MOON_EL0 + (MOON_EL1 * T);
    eclipticLongitude += MOON_EL_S_0 * sin(MOON_EL_S0_0 + MOON_EL_S1_0 * T);
    eclipticLongitude += MOON_EL_S_1 * sin(MOON_EL_S0_1 + MOON_EL_S1_1 * T);
    eclipticLongitude += MOON_EL_S_2 * sin(MOON_EL_S0_2 + MOON_EL_S1_2 * T);
    eclipticLongitude += MOON_EL_S_3 * sin(MOON_EL_S0_3 + MOON_EL_S1_3 * T);
    eclipticLongitude += MOON_EL_S_4 * sin(MOON_EL_S0_4 + MOON_EL_S1_4 * T);
    eclipticLongitude += MOON_EL_S_5 * sin(MOON_EL_S0_5 + MOON_EL_S1_5 * T);
    eclipticLongitude = fmod(eclipticLongitude,(2 * M_PI));
    
    // Ecliptic Lattitude in radians
    eclipticLattitude = 0;
    eclipticLattitude += MOON_ELA_S_0 * sin(MOON_ELA_S0_0 + MOON_ELA_S1_0 * T);
    eclipticLattitude += MOON_ELA_S_1 * sin(MOON_ELA_S0_1 + MOON_ELA_S1_1 * T);
    eclipticLattitude += MOON_ELA_S_2 * sin(MOON_ELA_S0_2 + MOON_ELA_S1_2 * T);
    eclipticLattitude += MOON_ELA_S_3 * sin(MOON_ELA_S0_3 + MOON_ELA_S1_3 * T);
    eclipticLattitude = fmod(eclipticLattitude,(2 * M_PI));
        
    // Horizontal parralax in radians
    horizontalParallax = MOON_HP0;
    horizontalParallax += MOON_HP_C_0 * cos(MOON_HP_C0_0 + MOON_HP_C1_0 * T);
    horizontalParallax += MOON_HP_C_1 * cos(MOON_HP_C0_1 + MOON_HP_C1_1 * T);
    horizontalParallax += MOON_HP_C_2 * cos(MOON_HP_C0_2 + MOON_HP_C1_2 * T);
    horizontalParallax += MOON_HP_C_3 * cos(MOON_HP_C0_3 + MOON_HP_C1_3 * T);
    horizontalParallax = fmod(horizontalParallax,(2 * M_PI));
    
    // Semidiameter in radians
    //SD = 0.2724 * moonPi

    // Precession Constants
    double T2 = pow(T, 2);
    a = MOON_PRCSN_P1_1 * T + MOON_PRCSN_P2_1 * T2;
    b = MOON_PRCSN_P1_2 * T + MOON_PRCSN_P2_2 * T2;
    c = MOON_PRCSN_P0_3 + MOON_PRCSN_P1_3 * T + MOON_PRCSN_P2_3 * T2;

	// Precession of Ecliptic Longitude and Lattitude
	eclipticLongitudeJ2K = eclipticLongitude - a + b * cos(eclipticLongitude + c) * tan(eclipticLattitude);
	eclipticLattitudeJ2K = eclipticLattitude - b * sin(eclipticLongitude + c);
    
    // Distance in Earth Radii
    rMagEarthRad = 1/sin(horizontalParallax);
    
    // Conversion to ECI
    rUnit[0] = cos(eclipticLattitudeJ2K) * cos(eclipticLongitudeJ2K);
    rUnit[1] = COS_OBLIQUEITY_ANGLE * cos(eclipticLattitudeJ2K) * sin(eclipticLongitudeJ2K) - SIN_OBLIQUEITY_ANGLE * sin(eclipticLattitudeJ2K);
    rUnit[2] = SIN_OBLIQUEITY_ANGLE * cos(eclipticLattitudeJ2K) * sin(eclipticLongitudeJ2K) + COS_OBLIQUEITY_ANGLE * sin(eclipticLattitudeJ2K);
    
    r[0] = rUnit[0] * rMagEarthRad * REQ_EARTH;
    r[1] = rUnit[1] * rMagEarthRad * REQ_EARTH;
	r[2] = rUnit[2] * rMagEarthRad * REQ_EARTH;

}

void sunEph(double tJ2000, double rUnit[3], double r[3])
{
    /*
    Sun Ephemeris relative to Earth

    Ref: 2017 Astronomical Almanac Section C5    

	Parameters
	----------
	tJ2000 : double
	    Time in seconds since J2000 epoch (January 1, 2000, 12:00 TT)

	Returns
	-------
    rUnit : Array
		Unit vector from Earth to Sun in ECI frame [x; y; z]
	r : Array
	    Vector from Earth to Sun in ECI frame in km [x; y; z]
    */
    
    double n;
    double L;
    double g;
    double R;
    double epsilon;
	double sunLambda;
    double sunLambdaRad;
    double epsilonRad;
    
    // Days since J2000
    n = tJ2000 * SEC2DAY;
    
    // Mean longitude of the sun in deg
    L = fmod((280.460 + 0.9856474 * n),360);
    
    // Mean anomaly in degrees
    g = fmod((357.528 + 0.9856003 * n),360);
    
    // Ecliptic Longitude in degrees
    sunLambda = fmod((L + 1.915 * sin(D2R * g) + 0.0020 * sin(D2R * 2 * g)),360);
    
    // Ecliptic Lattitude
    // beta = 0
    
    // Distance in AU
    R = 1.00014 - 0.01671 * cos(D2R * g) - 0.00014 * cos(D2R * 2 * g);
    
    // Obliquity of the ecliptic in degrees
    epsilon = fmod((23.439 - 0.0000004 * n),360);
    
    // Conversion to ECI
    sunLambdaRad = D2R * sunLambda;
    epsilonRad = D2R * epsilon;
    rUnit[0] = cos(sunLambdaRad);
    rUnit[1] = cos(epsilonRad)  *  sin(sunLambdaRad);
    rUnit[2] = sin(epsilonRad)  *  sin(sunLambdaRad);
    
    r[0] = rUnit[0] * R * AU;
	r[1] = rUnit[1] * R * AU;
	r[2] = rUnit[2] * R * AU;
}


