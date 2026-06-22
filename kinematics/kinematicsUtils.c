// file: kinematicsUtils.c
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <stdbool.h>
#include <stddef.h>

#include "kinematicsUtils.h"

/* common conversions */
#ifndef M_PI
#define M_PI            3.141592653589793
#endif
#ifndef M_PI_2
#define M_PI_2          1.5707963267948966
#endif
#ifndef M_E
#define M_E             2.718281828459045
#endif
#ifndef D2R
#define D2R             (M_PI/180.)
#endif
#ifndef R2D
#define R2D             180./M_PI
#endif
#ifndef NAN
#define NAN             sqrt((float)-1)
#endif
#ifndef RPM
#define RPM             (2.*M_PI/60.)
#endif
#ifndef DAY2SEC
#define DAY2SEC         (24.*3600.)
#endif
#ifndef SEC2DAY
#define SEC2DAY         (1.0 / 24. / 3600.)
#endif
#ifndef CENTURY2DAY
#define CENTURY2DAY     36525.
#endif
#ifndef DAY2CENTURY
#define DAY2CENTURY     1.0 / CENTURY2DAY
#endif
#ifndef SEC2CENTURY
#define SEC2CENTURY     3.1688087814028947e-10
#endif

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

/*
** Formation Functions
*/

void dcmInr2Ric(const double r[3], const double v[3], double RN[3][3]) {
    /*
    Compute the RIC frame DCM RN

    Parameters
    ----------
    r : Array
		Inertial position vector
    v : Array
		Inertial velocity vector

    Returns
    -------
    RN : 3x3 double
		DCM that maps from the inertial frame N to the RIC (i.e. orbit) frame
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
    Compute the RIC to Line-of-Sight frame (LOS) DCM LR
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

    Parameters
    ----------
    r : Array
		Chief inertial position
	v : Array
	    Chief inertial velocity
    r_d : Array
        Deputy inertial position
    v_d : Array
        Deputy inertial velocity

    Returns
    -------
	relPosRectRic : Array
	    Relative position in rectilinear RIC frame
	relVelRectRic : Array
	    Relative velocity in rectilinear RIC frame
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

    Parameters
    ----------
    r : Array
		Chief inertial position
	v : Array
	    Chief inertial velocity
	relPosRectRic : Array
	    Relative position in rectilinear RIC frame
	relVelRectRic : Array
	    Relative velocity in rectilinear RIC frame

    Returns
    -------
	r_d : Array
	    Deputy inertial position
	v_d : Array
	    Deputy inertial velocity
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

void oe2roe(const double oeChief[6], const double oeDeputy[6], double roe[6]) {
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
	roe : Array
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
	roe[0] = (aD - aC) / aC; // da
	double uC = trueAnomC + argPerC; // chief argument of latitude
	double uD = trueAnomD + argPerD; // deputy argument of latitude
	double lambdaC = uC + raanC * cos(iC); // chief mean longitude
	double lambdaD = uD + raanD * cos(iC); // deputy mean longitude (note: cos(iC), not cos(iD))
	roe[1] = lambdaD - lambdaC; // dlambda
	roe[2] = eD * cos(argPerD) - eC * cos(argPerC); // de_x
	roe[3] = eD * sin(argPerD) - eC * sin(argPerC); // de_y
	roe[4] = iD - iC; // di_x
    roe[5] = (raanD - raanC) * sin(iC); // di_y
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

void envAngles(const double rChief[3], const double rDeputy[3], const double rMoon[3], const double rSun[3],
    double *losEarthAng, double *losMoonAng, double *losSunAng) {
    /*
	Environmental Line-of-Sight Angles
   
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

/*
* Orbit Functions
*/

void oe2rv(double mu, const double oe[6], double r[3], double v[3])
{
    /*
    Translates the Keplerian orbit elements :

    =======================================
    a   semi - major axis         km
    e   eccentricity
    i   inclination               rad
    AN  ascending node            rad
    AP  argument of periapses     rad
    M   mean anomaly angle        rad
    =======================================

    to the inertial Cartesian position and velocity vectors.
    
    Parameters
    ----------
    mu : double
        Gravitational parameter in  km^3/s^2
    oe : Array
        Keplerian orbit elements (as defined above)

    Returns
    -------
    r : Array
        Inertial position vector of the spacecraft in km  [x;y;z]
    v : Array
        Inertial velocity vector of the spacecraft in km/s [vx;vy;vz]
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

    /* map Keplerian elements structure into local variables */
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
    Keplerian orbit elements where

    === ========================= =======
    a   semi-major axis           km
    e   eccentricity
    i   inclination               rad
    AN  ascending node            rad
    AP  argument of periapses     rad
    f   true anomaly angle        rad
    === ========================= =======

    Parameters
    ----------
    mu : double
        Gravitational parameter in  km^3/s^2
    r : Array
        Inertial position vector of the spacecraft in km  [x;y;z]
    v : Array
        Inertial velocity vector of the spacecraft in km/s [vx;vy;vz]

    Returns
    -------
    oe : Array
        Keplerian orbit elements (as defined above)
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
    eps = 1e-16;

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
    }
    else {
        /* parabolic case */
        a = 0.0;
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

    Parameters
    ----------
    oe : Array
        Keplerian orbit elements (as defined above)

    Returns
    -------
    ee : Array
        Equinoctial elements (as defined above)
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

    Parameters
    ----------
    mu : double
        Gravitational parameter
    r : Array
        Inertial position vector of the spacecraft in km  [x;y;z]
    v : Array
        Inertial velocity vector of the spacecraft in km/s [vx;vy;vz]

    Returns
    -------
    ee : Array
        Equinoctial elements (as defined above)
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

    Parameters
    ----------
    mu : double
        Gravitational parameter
    ee : Array
        Equinoctial elements (as defined above)

    Returns
    -------
    r : Array
        Inertial position vector of the spacecraft in km  [x;y;z]
    v : Array
        Inertial velocity vector of the spacecraft in km/s [vx;vy;vz]
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
    /*
    Eccentric to True Anomaly

    Parameters
    ----------
    E : double
        Eccentric anomaly
    e : double
        eccentricity

    Returns
    -------
    f : double
        True anomaly
    */

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
    /*
    Eccentric to Mean Anomaly 

    Parameters
    ----------
    E : double
        Eccentric anomaly
    e : double
        eccentricity

    Returns
    -------
    M : double
        Mean anomaly
    */

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
    /*
    True to Eccentric Anomaly

    Parameters
    ----------
    f : double
        True anomaly
    e : double
        eccentricity

    Returns
    -------
    E : double
        Eccentric anomaly
    */

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
    /*
    Mean to Eccentric Anomaly (Solves Kepler's equation using Newton's method)

    Parameters
    ----------
    M : double
        Mean anomaly
    e : double
        eccentricity

    Returns
    -------
    E : double
        Eccentric anomaly
    */

    double eps = 1e-16;
    double dE;
    double EPrev = M;
    double E;
    int    count = 0;
    int    maxIteration = 200;

    if((e >= 0) && (e < 1)) {
        for (int j = 0; j < maxIteration; ++j) {
            E = M + e * sin(EPrev);
            dE = E - EPrev;
            EPrev = E;
            count += 1;
            if(fabs(dE) < eps) {
                return E;
            }
        }
    } else {
        E = NAN;
    }

    return E;
}
