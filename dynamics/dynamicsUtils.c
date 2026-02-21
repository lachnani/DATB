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

// Helper vector math functions (implement or use your own library)
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
     *   This function uses a curve fit based on
     *   atmospheric data from the Standard Atmosphere 1976 Data. This
     *   function is valid for altitudes ranging from 100km to 1000km.
     *
     *   Note: This code can only be applied to spacecraft orbiting the Earth
     * Inputs:
     *   alt = altitude in km
     * Outputs:
     *   density = density at the given altitude in kg/m^3
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
** Formation Functions
*/

void dcmInr2Ric(const double r[3], const double v[3], double RN[3][3]) {
    /*
    Compute the RIC frame DCM RN
    :param r: inertial position
    :param v: inertial velocity
    :return: RN: DCM that maps from the inertial frame N to the RIC (i.e. orbit) frame
    */
    double ir[3], h[3], ih[3], itheta[3];
    double norm_r = v3_norm(r);
    for (int i = 0; i < 3; ++i) ir[i] = r[i] / norm_r;
    v3_cross(r, v, h);
    double norm_h = v3_norm(h);
    for (int i = 0; i < 3; ++i) ih[i] = h[i] / norm_h;
    v3_cross(ih, ir, itheta);
    for (int i = 0; i < 3; ++i) {
        RN[0][i] = ir[i];
        RN[1][i] = itheta[i];
        RN[2][i] = ih[i];
    }
}

void dcmRic2Los(const double relPosRectRic[3], double LR[3][3]) {
    /*
    Compute the RIC to LOS DCM LR
    Assumes the LOS frame is defined with the z-axis along the line of sight vector,
	the x-axis in the chief orbital plane, and the y-axis completing the right-handed system.

	Parameters
	----------
    relPosRectRic : Array
		Rectilinear RIC frame relative position

	Returns
	-------
	LR : 3x3 double
	    DCM that maps from the RIC frame to the LOS frame
    */
	double rho = v3_norm(relPosRectRic);
    double ir[3];
    if (rho > 0) {
        for (int i = 0; i < 3; ++i)
            ir[i] = relPosRectRic[i] / rho;
    }
    else {
        v3_zero(ir);
    }
    double itheta[3] = { 0, 0, 1 }; // LOS z-axis
    double ih[3];
    v3_cross(itheta, ir, ih);
    double norm_ih = v3_norm(ih);
    if (norm_ih > 0) {
        for (int i = 0; i < 3; ++i)
            ih[i] /= norm_ih;
    }
    else {
        v3_zero(ih);
    }
    v3_cross(ir, ih, itheta);
    for (int i = 0; i < 3; ++i) {
        LR[0][i] = ih[i];
        LR[1][i] = itheta[i];
        LR[2][i] = ir[i];
	}
}

void rv2ric(const double r[3], const double v[3], const double r_d[3], const double v_d[3],
    double relPosRectRic[3], double relVelRectRic[3]) {
    /*
    Express the deputy position and velocity vector as chief by the chief RIC frame.

    :param r: chief inertial position
    :param v: chief inertial velocity
    :param r_d: deputy inertial position
    :param v_d: deputy inertial velocity
    :return: relRvRic: RIC frame relative position and velocity vectors
    */
    double RN[3][3];
    dcmInr2Ric(r, v, RN);
    double norm_r = v3_norm(r);
    double cross_rv[3];
    v3_cross(r, v, cross_rv);
    double fDot = v3_norm(cross_rv) / (norm_r * norm_r);
    double argP_RN_R[3] = { 0, 0, fDot };
    double dr[3], dv[3], tmp[3];
    for (int i = 0; i < 3; ++i) {
        dr[i] = r_d[i] - r[i];
        dv[i] = v_d[i] - v[i];
    }
    mat3_v3_mul(RN, dr, relPosRectRic);
    mat3_v3_mul(RN, dv, tmp);
    double cross_argP[3];
    v3_cross(argP_RN_R, relPosRectRic, cross_argP);
    for (int i = 0; i < 3; ++i)
        relVelRectRic[i] = tmp[i] - cross_argP[i];
}

void ric2rv(const double r[3], const double v[3], const double relPosRectRic[3], const double relVelRectRic[3],
    double r_d[3], double v_d[3]) {
    /*
    Map the deputy position and velocity vector relative to the chief RIC frame to inertial frame.

    :param r: chief inertial position
    :param v: chief inertial velocity
    :param relPosRectRic: RIC frame relative position
    :param relVelRectRic: RIC frame relative velocity
    :return: r_d: deputy inertial position
             v_d: deputy inertial velocity
    */
    double RN[3][3], NR[3][3];
    dcmInr2Ric(r, v, RN);
    mat3_transpose(RN, NR);
    double norm_r = v3_norm(r);
    double cross_rv[3];
    v3_cross(r, v, cross_rv);
    double fDot = v3_norm(cross_rv) / (norm_r * norm_r);
    double argP_HN_H[3] = { 0, 0, fDot };
    double tmp1[3], tmp2[3], cross_argP[3];
    mat3_v3_mul(NR, relPosRectRic, tmp1);
    for (int i = 0; i < 3; ++i)
        r_d[i] = r[i] + tmp1[i];
    v3_cross(argP_HN_H, relPosRectRic, cross_argP);
    for (int i = 0; i < 3; ++i)
        tmp2[i] = relVelRectRic[i] + cross_argP[i];
    mat3_v3_mul(NR, tmp2, tmp1);
    for (int i = 0; i < 3; ++i)
        v_d[i] = v[i] + tmp1[i];
}

void curvRic2rectRic(const double r[3], const double v[3],
    const double relPosCurvRic[3], const double relVelCurvRic[3],
    double relPosRectRic[3], double relVelRectRic[3]) {
    /*
    Converts the elements in the curvilinear RIC frame to the rectilinear RIC frame

    Parameters
    ----------
    r : Array
        Chief position
    v : Array
        Chief velocity
    relPosCurvRic : Array
        Curvilinear RIC frame relative position
    relVelCurvRic : Array
        Curvilinear RIC frame relative vel

    Returns
    -------
    relPosRectRic : Array
        Rectilinear RIC frame relative position
    relVelRectRic : Array
        Rectilinear RIC frame relative velocity
    */

	/* Chief Parameters */ 
    double rMag = v3_norm(r);
    double rMagDot = v3_dot(r, v) / rMag;

    /* Curvilinear Paramters */
    double dr = relPosCurvRic[0];
    double dTheta = relPosCurvRic[1] / rMag;
    double dPhi = relPosCurvRic[2] / rMag;
    double drDot = relVelCurvRic[0];
    double dThetaDot = (relVelCurvRic[1] - rMagDot * dTheta) / rMag;
    double dPhiDot = (relVelCurvRic[2] - rMagDot * dPhi) / rMag;

    /* Position Transformation */
    double rP = rMag + dr;
    double rPDot = rMagDot + drDot;
    double cdTheta = cos(dTheta), sdTheta = sin(dTheta);
    double cdPhi = cos(dPhi), sdPhi = sin(dPhi);

    relPosRectRic[0] = rP * cdTheta * cdPhi - rMag;
    relPosRectRic[1] = rP * sdTheta * cdPhi;
    relPosRectRic[2] = rP * sdPhi;
    relVelRectRic[0] = rPDot * cdTheta * cdPhi - rP * (dThetaDot * sdTheta * cdPhi + dPhiDot * cdTheta * sdPhi) - rMagDot;
    relVelRectRic[1] = rPDot * sdTheta * cdPhi + rP * (dThetaDot * cdTheta * cdPhi - dPhiDot * sdTheta * sdPhi);
    relVelRectRic[2] = rPDot * sdPhi + rP * dPhiDot * cdPhi;
}

void rectRic2curvRic(const double r[3], const double v[3],
    const double relPosRectRic[3], const double relVelRectRic[3],
    double relPosCurvRic[3], double relVelCurvRic[3]) {
    /*
    Parameters
    ----------
    r : Array
        Chief position
    v : Array
        Chief velocity
    relPosRectRic : Array
        Rectilinear RIC frame relative position
    relVelRectRic : Array
        Rectilinear RIC frame relative velocity

    Returns
    -------
    relPosCurvRic : Array
        Curvilinear RIC frame relative position
    R : Array
        Curvilinear RIC frame relative vel
    */

    /* Chief Parameters */
    double rMag = v3_norm(r);
    double rMagDot = v3_dot(r, v) / rMag;

    /* Position Transformation */
	double rMagX = rMag + relPosRectRic[0];
	double rMagX2 = rMagX * rMagX;
	double Y2 = relPosRectRic[1] * relPosRectRic[1];
	double dr = sqrt(rMagX2 + Y2 + relPosRectRic[2] * relPosRectRic[2]) - rMag;
	double dTheta = atan2(relPosRectRic[1], rMagX);
	double dPhi = atan2(relPosRectRic[2], sqrt(rMagX2 + Y2));

	/* Velocity Transformation */
    double drDot = (rMagX * relVelRectRic[0] + relPosRectRic[1] * relVelRectRic[1] + relPosRectRic[2] * relVelRectRic[2]) / (rMag + dr);
	double dThetaDot = (rMagX * relVelRectRic[1] - relPosRectRic[1] * relVelRectRic[0]) / (rMagX2 + Y2);
	double dPhiDot = (relVelRectRic[2] * (rMagX2 + Y2) - relPosRectRic[2] * (rMagX * relVelRectRic[0] + relPosRectRic[1] * relVelRectRic[1])) / ((rMag + dr) * (rMag + dr) * sqrt(rMagX2 + Y2));

    /* Output */
    relPosCurvRic[0] = dr;
    relPosCurvRic[1] = rMag * dTheta;
    relPosCurvRic[2] = rMag * dPhi;
    relVelCurvRic[0] = drDot;
    relVelCurvRic[1] = rMagDot * dTheta + rMag * dThetaDot;
    relVelCurvRic[2] = rMagDot * dPhi + rMag * dPhiDot;
}

void clroe2ric(const double L[6], double n, double t, double relPosRic[3], double relVelRic[3]) {
    /*
    Parameters
    ----------
    L : 6x1 double
        Circular Linearized Relative Orbit Element state.
        A0 = L[0]
        alpha = L[1]
        xOff = L[2]
        yOff = L[3]
        B0 = L[4]
        beta = L[5]
    n : double
        Orbit mean motion.
    t : double
        Time since CLROE time of validity

    Returns
    -------
    relPosRic : Array
        RIC frame relative position
    relVelRic : Array
        RIC frame relative velocity
    */
    double nt = n * t;
    relPosRic[0] = L[0] * cos(nt + L[1]) + L[2];
    relPosRic[1] = -2 * L[0] * sin(nt + L[1]) + L[3] - 1.5 * nt * L[2];
    relPosRic[2] = L[4] * cos(nt + L[5]);
    relVelRic[0] = -n * L[0] * sin(nt + L[1]);
    relVelRic[1] = -2 * n * L[0] * cos(nt + L[1]) - 1.5 * n * L[2];
    relVelRic[2] = -n * L[4] * sin(nt + L[5]);
}

void ric2clroe(const double relPosRic[3], const double relVelRic[3], double n, double t, double L[6]) {
    /*
    Parameters
    ----------
    relPosRic : Array
        RIC frame relative position
    relVelRic : Array
        RIC frame relative velocity
    n : double
        Orbit mean motion.
    t : double
        Time since CLROE time of validity

    Returns
    -------
    L : 6x1 double
        Circular Linearized Relative Orbit Element state.
        A0 = L[0]
        alpha = L[1]
        xOff = L[2]
        yOff = L[3]
        B0 = L[4]
        beta = L[5]
    */
    memset(L, 0, 6 * sizeof(double));
    if (n > 0) {
        double n2 = n * n;
        double nt = n * t;
        L[0] = sqrt(9 * n2 * relPosRic[0] * relPosRic[0]
            + relVelRic[0] * relVelRic[0]
            + 12 * n * relPosRic[0] * relVelRic[1]
            + 4 * relVelRic[1] * relVelRic[1]) / n;
        if (fabs(L[0]) > 0) {
            L[1] = atan2(-relVelRic[0], -3 * n * relPosRic[0] - 2 * relVelRic[1]) - nt;
        }
        L[2] = 4 * relPosRic[0] + 2 * relVelRic[1] / n;
        L[3] = -2 * relVelRic[0] / n + relPosRic[1] + (6 * n * relPosRic[0] + 3 * relVelRic[1]) * t;
        L[4] = sqrt(n2 * relPosRic[2] * relPosRic[2] + relVelRic[2] * relVelRic[2]) / n;
        if (fabs(L[4]) > 0) {
            L[5] = atan2(-relVelRic[2], n * relPosRic[2]) - nt;
        }
    }
}

void oe2dAmico(const double oeChief[6], const double oeDeputy[6], double dAmico[6]) {
    /*
	Compute the d'Amico relative orbital elements from chief and deputy osculating orbital elements.

	Parameters
	----------
	oeChief : Array
	    Chief osculating orbital elements [a, e, i, RAAN, argPer, trueAnom]
	oeDeputy : Array
	    Deputy osculating orbital elements [a, e, i, RAAN, argPer, trueAnom]

	Returns
	-------
	dAmico : Array
	    d'Amico relative orbital elements [da, dlambda, de_x, de_y, di_x, di_y]
    */
	double aC = oeChief[0];
	double eC = oeChief[1];
	double iC = oeChief[2];
	double raanC = oeChief[3];
	double argPerC = oeChief[4];
	double trueAnomC = oeChief[5];
	double aD = oeDeputy[0];
	double eD = oeDeputy[1];
	double iD = oeDeputy[2];
	double raanD = oeDeputy[3];
	double argPerD = oeDeputy[4];
	double trueAnomD = oeDeputy[5];
	dAmico[0] = (aD - aC) / aC; // da
	double uC = trueAnomC + argPerC; // chief argument of latitude
	double uD = trueAnomD + argPerD; // deputy argument of latitude
	double lambdaC = uC + raanC * cos(iC); // chief mean longitude
	double lambdaD = uD + raanD * cos(iC); // deputy mean longitude (note: cos(iC), not cos(iD))
	dAmico[1] = lambdaD - lambdaC; // dlambda
	dAmico[2] = eD * cos(argPerD) - eC * cos(argPerC); // de_x
	dAmico[3] = eD * sin(argPerD) - eC * sin(argPerC); // de_y
	dAmico[4] = iD - iC; // di_x
    dAmico[5] = (raanD - raanC) * sin(iC); // di_y
}

void measParams(const double relPosRic[3], const double relVelRic[3], double *rng, double *rngRate, double *az, double *el) {
    /*
    Parameters
    ----------
    relPosRic : Array
        RIC frame relative position
    relVelRic : Array
        RIC frame relative velocity

    Returns
    -------
    rng : double
		Range
	rngRate : double
	    Range rate
	az : double
	    Azimuth angle
	el : double
	    Elevation angle
    */
	*rng = v3_norm(relPosRic);
    if (rng > 0) {
        *rngRate = v3_dot(relPosRic, relVelRic) / *rng;
    }
    else {
        *rngRate = 0.0;
	}
    *az = atan2(relPosRic[1], relPosRic[0]);
	*el = atan2(relPosRic[2], sqrt(pow(relPosRic[0],2) + pow(relPosRic[1],2)));
}

/*
* Orbit Functions
*/

// Orbit acceleration (gravity, perturbations, thrust)
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
    /*Orbit Acceleration (Gravity, Perturbations, and Thrust)*/
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

// Propagate orbit using RK4
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

	/* Propagate orbit using RK4 method */

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

    :param r: Cartesian Position vector in kilometers [x;y;z].
    :param body: body variable
	:param rBody: position vector of the body exerting gravity
    :return: agrav, The total acceleration vector due to the gravity in 
                km/sec^2 [accelx;accely;accelz]
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

void    oe2rv(double mu, const double oe[6], double r[3], double v[3])
{

    /*
    Translates the orbit elements :

    =======================================
    a   semi - major axis         km
    e   eccentricity
    i   inclination               rad
    AN  ascending node            rad
    AP  argument of periapses     rad
    M   mean anomaly angle        rad
    =======================================

    to the inertial Cartesian position and velocity vectors.
    The attracting body is specified through the supplied
    gravitational constant mu(units of km ^ 3 / s ^ 2).

     */


    double e;                   /* eccentricty */
    double a;                   /* semi-major axis */
	double M;                   /* mean anomaly */
    double f;                   /* true anomaly */
    double rMag;                /* orbit radius */
    double i;                   /* orbit inclination angle */
    double p;                   /* the parameter or the semi-latus rectum */
    double argP;                /* argument of perigee */
    double RAAN;                /* argument of the ascending node */
    double theta;               /* true latitude theta = argP + f */
    double h;                   /* orbit angular momentum magnitude */
    double eps;                 /* small numerical value parameter */

    /* define what is a small numerical value */
    eps = 1e-12;

    /* map classical elements structure into local variables */
    a = oe[0];
    e = oe[1];
    i = oe[2];
    RAAN = oe[3];
    argP = oe[4];
    M = oe[5];

	f = E2f(M2E(M, e), e); /* true anomaly from mean anomaly */

    if ((fabs(e - 1.0) < eps) && (fabs(a) > eps)) { /* 1D rectilinear elliptic/hyperbolic orbit case */

        /* oe2rv does not support rectilinear orbits*/
    }
    else { /* general 2D orbit case */
        /* evaluate semi-latus rectum */
        if (fabs(a) > eps) {
            p = a * (1 - e * e);        /* elliptic or hyperbolic */
        }
        else {
            /* oe2rv does not support parabolic orbits*/
        }

        rMag = p / (1 + e * cos(f));   /* orbit radius */
        theta = argP + f;             /* true latitude angle */
        h = sqrt(mu * p);           /* orbit ang. momentum mag. */

        r[0] = rMag * (cos(RAAN) * cos(theta) - sin(RAAN) * sin(theta) * cos(i));
        r[1] = rMag * (sin(RAAN) * cos(theta) + cos(RAAN) * sin(theta) * cos(i));
        r[2] = rMag * (sin(theta) * sin(i));

        v[0] = -mu / h * (cos(RAAN) * (sin(theta) + e * sin(argP)) + sin(RAAN) * (cos(
            theta) + e * cos(argP)) * cos(i));
        v[1] = -mu / h * (sin(RAAN) * (sin(theta) + e * sin(argP)) - cos(RAAN) * (cos(
            theta) + e * cos(argP)) * cos(i));
        v[2] = -mu / h * (-(cos(theta) + e * cos(argP)) * sin(i));
    }
}

void rv2oe(double mu, const double r[3], const double v[3], double oe[6])
{
    /*
    Translates the orbit elements inertial Cartesian position
    vector r and velocity vector v into the corresponding
    classical orbit elements where

    === ========================= =======
    a   semi-major axis           km
    e   eccentricity
    i   inclination               rad
    AN  ascending node            rad
    AP  argument of periapses     rad
    f   true anomaly angle        rad
    === ========================= =======

    If the orbit is rectilinear, then this will be the eccentric or hyperbolic anomaly

    The attracting body is specified through the supplied
    gravitational constant mu (units of km^3/s^2).
    */

    double hVec[3];             /* orbit angular momentum vector */
    double ihHat[3];            /* normalized orbit angular momentum vector */
    double h;                   /* orbit angular momentum magnitude */
	double e;                   /* orbit eccentricity */
    double a;                   /* semi-major axis */
    double M;                   /* mean anomaly */
    double f;                   /* true anomaly */
    double i;                   /* orbit inclination angle */
    double p;                   /* the parameter or the semi-latus rectum */
    double argP;                /* argument of perigee */
    double RAAN;                /* argument of the ascending node */
	double rPeriap;             /* perigee radius */
	double rApoap;              /* apogee radius */
	double alpha;               /* semi-major axis parameter */
    double v3[3];               /* temp vector */
    double n1Hat[3];            /* 1st inertial frame base vector */
    double n3Hat[3];            /* 3rd inertial frame base vector */
    double nVec[3];             /* line of nodes vector */
    double inHat[3];            /* normalized line of nodes vector */
    double irHat[3];            /* normalized position vector */
    double rMag;                /* current orbit radius */
    double vMag;                   /* orbit velocity magnitude */
    double eVec[3];             /* eccentricity vector */
    double ieHat[3];            /* normalized eccentricity vector */
    double eps;                 /* small numerical value parameter */

    /* define what is a small numerical value */
    eps = 1e-12;

    /* define inertial frame axes */
    v3_set(1.0, 0.0, 0.0, n1Hat);
    v3_set(0.0, 0.0, 1.0, n3Hat);

    /* norms of position and velocity vectors */
    rMag = v3_norm(r);
    vMag = v3_norm(v);
    v3_normalize(r, irHat);

    /* Calculate the specific angular momentum and its magnitude */
    v3_cross(r, v, hVec);
    h = v3_norm(hVec);
    v3_normalize(hVec, ihHat);
    p = h * h / mu;

    /* Calculate the line of nodes */
    v3_cross(n3Hat, hVec, nVec);
    if (v3_norm(nVec) < eps) {
        /* near equatorial orbits */
        v3_copy(n1Hat, inHat);
    }
    else {
        v3_normalize(nVec, inHat);
    }

    /* Orbit eccentricity vector */
    v3_scale(r, vMag * vMag / mu - 1.0 / rMag, eVec);
    v3_scale(v, v3_dot(r, v) / mu, v3);
    v3_sub(eVec, v3, eVec);
    e = v3_norm(eVec);
    rPeriap = p / (1.0 + e);

    /* Orbit eccentricity unit vector */
    if (e > eps) {
        v3_normalize(eVec, ieHat);
    }
    else {
        /* near circular orbits, make i_e_hat equal to line of nodes */
        v3_copy(inHat, ieHat);
    }

    /* compute semi-major axis */
    alpha = 2.0 / rMag - vMag * vMag / mu;
    if (fabs(alpha) > eps) {
        /* elliptic or hyperbolic case */
        a = 1.0 / alpha;
        rApoap = p / (1.0 - e);
    }
    else {
        /* parabolic case */
        a = 0.0;
        rApoap = 0.0;
    }

    /* Calculate the inclination */
    i = acos(hVec[2] / h);

    /* Calculate Ascending Node RAAN */
    v3_cross(n1Hat, inHat, v3);
    RAAN = atan2(v3[2], inHat[0]);
    if (RAAN < 0.0) {
        RAAN += 2 * M_PI;
    }

    /* Calculate Argument of Periapses argP */
    v3_cross(inHat, ieHat, v3);
    argP = atan2(v3_dot(ihHat, v3), v3_dot(inHat, ieHat));
    if (argP < 0.0) {
        argP += 2 * M_PI;
    }

    /* Calculate true anomaly angle f */
    v3_cross(ieHat, irHat, v3);
    f = atan2(v3_dot(ihHat, v3), v3_dot(ieHat, irHat));
    if (f < 0.0) {
        f += 2 * M_PI;
    }

	M = E2M(f2E(f, e), e); /* mean anomaly from true anomaly */

	oe[0] = a;          /* semi-major axis */
	oe[1] = e;          /* eccentricity */
	oe[2] = i;          /* inclination angle */
	oe[3] = RAAN;       /* right ascension of ascending node */
	oe[4] = argP;       /* argument of periapses */
	oe[5] = M;          /* mean anomaly angle */

    return;
}

void oe2ee(const double oe[6], double ee[6])
{
    /*
    Translates the orbit elements:

    === ========================= =======
    a   semi-major axis           km
    e   eccentricity
    i   inclination               rad
    AN  ascending node            rad
    AP  argument of periapses     rad
    M   mean anomaly angle        rad
    === ========================= =======
    
    into the corresponding
    equinoctial elements where

    === =================================== =======
    a   semi-major axis                     km
    l   mean longitude (RAAN + argP + M)    rad
    P1  e * sin(RAAN + argP)             
    P2  e * cos(RAAN + argP)           
    Q1  tan(i/2) * sin(RAAN)   
    Q2  tan(i/2) * cos(RAAN)   
    === =================================== =======

    If the orbit is rectilinear, then this will be the eccentric or hyperbolic anomaly

    The attracting body is specified through the supplied
    gravitational constant mu (units of km^3/s^2).
    */

    double a = oe[0];
    double e = oe[1];
    double i = oe[2];
    double RAAN = oe[3];
    double argP = oe[4];
    double M = oe[5];

    double omegaBar = RAAN + argP;
    double l = omegaBar + M;
    double P1 = e * sin(omegaBar);
    double P2 = e * cos(omegaBar);
    double Q1 = tan(i / 2.0) * sin(RAAN);
    double Q2 = tan(i / 2.0) * cos(RAAN);

    ee[0] = a;
    ee[1] = l;
    ee[2] = P1;
    ee[3] = P2;
    ee[4] = Q1;
    ee[5] = Q2;
}


void rv2ee(double mu, const double r[3], const double v[3], double ee[6])
{
    /*
    Translates the orbit elements inertial Cartesian position
    vector r and velocity vector v into the corresponding
    equinoctial elements where

    === =================================== =======
    a   semi-major axis                     km
    l   mean longitude (RAAN + argP + M)    rad
    P1  e * sin(RAAN + argP)             
    P2  e * cos(RAAN + argP)           
    Q1  e * tan(i/2) * sin(RAAN)   
    Q2  e * tan(i/2) * cos(RAAN)   
    === =================================== =======

    The attracting body is specified through the supplied
    gravitational constant mu (units of km^3/s^2).
    
    Ref: "Satellite Orbits," Montenbruck & Gill
    */

    const double eps = 1e-13;
    double rMag = v3_norm(r);
    double vMag = v3_norm(v);

    // Check for NaN input
    if (isnan(r[0] + r[1] + r[2]) || isnan(v[0] + v[1] + v[2])) {
        for (int i = 0; i < 6; ++i) ee[i] = NAN;
        return;
    }

    // Specific angular momentum
    double hVec[3];
    v3_cross(r, v, hVec);
    double h = v3_norm(hVec);

    // W = hVec / h
    double W[3];
    v3_scale(hVec, 1.0 / h, W);

    // Inclination
    double i = atan2(sqrt(W[0] * W[0] + W[1] * W[1]), W[2]);

    double a, l, P1, P2, Q1, Q2;
    if (fabs(i - M_PI) > eps) {
        // Semimajor axis
        a = 1.0 / (2.0 / rMag - vMag * vMag / mu);

        // Q1, Q2
        Q1 = W[0] / (1.0 + W[2]);
        Q2 = -W[1] / (1.0 + W[2]);

        // f and g vectors
        double denom = 1.0 + Q1 * Q1 + Q2 * Q2;
        double fvec[3] = {
            (1.0 - Q1 * Q1 + Q2 * Q2) / denom,
            (2.0 * Q1 * Q2) / denom,
            (-2.0 * Q1) / denom
        };
        double gvec[3] = {
            (2.0 * Q1 * Q2) / denom,
            (1.0 + Q1 * Q1 - Q2 * Q2) / denom,
            (2.0 * Q2) / denom
        };
        v3_normalize(fvec, fvec);
        v3_normalize(gvec, gvec);

        // A = cross(v, hVec) - mu*r/rMag
        double vh[3], mur[3], A[3];
        v3_cross(v, hVec, vh);
        v3_scale(r, mu / rMag, mur);
        v3_sub(vh, mur, A);

        // P1, P2
        P1 = v3_dot(A, gvec) / mu;
        P2 = v3_dot(A, fvec) / mu;

        // In-plane coordinates
        double X1 = v3_dot(r, fvec);
        double Y1 = v3_dot(r, gvec);
        double rCalc[3] = { 0, 0, 0 };
        for (int j = 0; j < 3; ++j) rCalc[j] = X1 * fvec[j] + Y1 * gvec[j];
        if (v3_norm(rCalc) < 1e-30 || fabs(v3_norm(rCalc) - rMag) / rMag > eps) {
            for (int i = 0; i < 6; ++i) ee[i] = NAN;
            return;
        }

        // Eccentric longitude
        double e2 = P1 * P1 + P2 * P2;
        double Beta = 1.0 / (1.0 + sqrt(1.0 - e2));
        double sqrt1me2 = sqrt(1.0 - e2);
        double cosF = P2 + ((1.0 - P2 * P2 * Beta) * X1 - P1 * P2 * Beta * Y1) / (a * sqrt1me2);
        double sinF = P1 + ((1.0 - P1 * P1 * Beta) * Y1 - P1 * P2 * Beta * X1) / (a * sqrt1me2);
        double F = atan2(sinF, cosF);

        // Mean longitude
        l = F - P2 * sinF + P1 * cosF;
        // Wrap l to [0, 2pi)
        l = fmod(l, 2.0 * M_PI);
        if (l < 0) l += 2.0 * M_PI;
    }
    else {
        a = l = P1 = P2 = Q1 = Q2 = NAN;
    }

    ee[0] = a;
    ee[1] = l;
    ee[2] = P1;
    ee[3] = P2;
    ee[4] = Q1;
    ee[5] = Q2;
}

void ee2rv(double mu, const double ee[6], double r[3], double v[3])
{
    /*
    Translates the equinoctial elements:

    === =================================== =======
    a   semi-major axis                     km
    l   mean longitude (RAAN + argP + M)    rad
    P1  e * sin(RAAN + argP)             
    P2  e * cos(RAAN + argP)           
    Q1  tan(i/2) * sin(RAAN)   
    Q2  tan(i/2) * cos(RAAN)   
    === =================================== =======

    to the inertial Cartesian position and velocity vectors.
    The attracting body is specified through the supplied
    gravitational constant mu (units of km^3/s^2).
    
    Ref: Battin page 494
    */

    double a = ee[0];
    double l = ee[1];
    double P1 = ee[2];
    double P2 = ee[3];
    double Q1 = ee[4];
    double Q2 = ee[5];

    l = wrapTo2Pi(l);

    double omegaBar = atan2(P1, P2);

    // Eccentricity
    double e2 = P1 * P1 + P2 * P2;
    double e = sqrt(e2);

    // Anomalies
    double M = wrapTo2Pi(l - omegaBar);
    double E = M2E(M, e);
    double f = E2f(E, e);

    // True longitude
    double L = omegaBar + f;

    // Orbit angular momentum and parameter
    double h = sqrt(mu * a * (1.0 - e2));
    double p = h * h / mu;

    // Position magnitude
    double rMag = a * (1.0 - e2) / (1.0 + e * cos(f));

    // Position and velocity in the equinoctial frame
    double rEq[3] = { rMag * cos(L), rMag * sin(L), 0.0 };
    double vEq[3] = { (h / p) * (-P1 - sin(L)), (h / p) * (P2 + cos(L)), 0.0 };

    // Rotation matrix from equinoctial to inertial frame
    double Q1_2 = Q1 * Q1, Q2_2 = Q2 * Q2;
    double denom = 1.0 + Q1_2 + Q2_2;
    double Eq2Inr[3][3] = {
        { (1.0 - Q1_2 + Q2_2) / denom, (2.0 * Q1 * Q2) / denom, (2.0 * Q1) / denom },
        { (2.0 * Q1 * Q2) / denom, (1.0 + Q1_2 - Q2_2) / denom, (-2.0 * Q2) / denom },
        { (-2.0 * Q1) / denom, (2.0 * Q2) / denom, (1.0 - Q1_2 - Q2_2) / denom }
    };

    // Transform to inertial frame
    mat3_v3_mul(Eq2Inr, rEq, r);
    mat3_v3_mul(Eq2Inr, vEq, v);
}

double E2f(double Ecc, double e)
{
    double f;

    if((e >= 0) && (e < 1)) {
        f = 2 * atan2(sqrt(1 + e) * sin(Ecc / 2), sqrt(1 - e) * cos(Ecc / 2));
    } else {
        f = NAN;
    }

    return f;
}

double E2M(double Ecc, double e)
{
    double M;

    if((e >= 0) && (e < 1)) {
        M = Ecc - e * sin(Ecc);
    } else {
        M = NAN;
    }

    return M;
}

double f2E(double f, double e)
{
    double Ecc;

    if((e >= 0) && (e < 1)) {
        Ecc = 2 * atan2(sqrt(1 - e) * sin(f / 2), sqrt(1 + e) * cos(f / 2));
    } else {
        Ecc = NAN;
    }

    return Ecc;
}

double M2E(double M, double e)
{
    double eps = 1e-13;
    double dE = 10 * eps;
    double EPrev = M;
    double E;
    int    count = 0;
    int    maxIteration = 200;

    if((e >= 0) && (e < 1)) {
        while(fabs(dE) > eps) {
            E = M + e * sin(EPrev);
            dE = E - EPrev;
            EPrev = E;
            count += 1;
            if(++count > maxIteration) {
                dE = 0.;
            }
        }
    } else {
        E = NAN;
    }

    return E;
}

void dragPerturb(const double Cd, const double normA, const double r[3], const double v[3], double ad[3]) {

    /*
     * Inputs:
     *   Cd = drag coefficient of the spacecraft
	 *   normA = cross-sectional area of the spacecraft in m^2 divided by the mass of the spacecraft in kg
     *   r = Inertial position vector of the spacecraft in km  [x;y;z]
     *   v = Inertial velocity vector of the spacecraft in km/s [vx;vy;vz]
     * Outputs:
     *   as = The inertial acceleration vector due to atmospheric
     *             drag in km/sec^2
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
    Ref: "An Alternative Lunar Ephemeris Model for On-board Flight Software Use" by David Simpson
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

void envAngles(const double rChief[3], const double rDeputy[3], const double rMoon[3], const double rSun[3],
    double *losEarthAng, double *losMoonAng, double *losSunAng) {
    /*
    Parameters
    ----------
    rChief : Array
		Chief position in ECI frame
	rDeputy : Array
		Deputy position in ECI frame
	rMoon : Array
	    Moon position in ECI frame
	rSun : Array
	    Sun position in ECI frame

    Returns
    -------
    losEarthAng : double
		Line-of-sight angle between deputy-chief vector and deputy-Earth vector in radians
	losMoonAng : double
		Line-of-sight angle between deputy-chief vector and deputy-Moon vector in radians
	losSunAng : double
		Line-of-sight angle between deputy-chief vector and deputy-Sun vector in radians
    */
	double losEci[3];
    double earthVecEci[3];
	double moonVecEci[3];
	double sunVecEci[3];
    
	v3_sub(rChief, rDeputy, losEci);
	v3_set(-rDeputy[0], -rDeputy[1], -rDeputy[2], earthVecEci) ; // Earth is at origin
	v3_sub(rMoon, rDeputy, moonVecEci);
	v3_sub(rSun, rDeputy, sunVecEci);

	*losEarthAng = v3_angle(losEci, earthVecEci);
	*losMoonAng = v3_angle(losEci, moonVecEci);
	*losSunAng = v3_angle(losEci, sunVecEci);
}


