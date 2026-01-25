/*
 *  dynamicsConstants.h
 */

#ifndef _DYNAMICS_CONSTANTS_H_
#define _DYNAMICS_CONSTANTS_H_

#include <math.h>

#ifndef G_UNIVERSIAL
#define G_UNIVERSIAL    6.67259e-20 /* universial Gravitational Constant, units are in km^3/s^2/kg */
#endif
#ifndef AU
#define AU              149597870.693 /* astronomical unit in units of kilometers */
#endif
#ifndef AU2KM
#define AU2KM           149597870.693 /* convert astronomical unit to kilometers */
#endif

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

#ifndef EARTH_GRAV
#define EARTH_GRAV 9.80665
#endif

/* Gravitational Constants mu = G*m, where m is the planet of the attracting planet.  All units are km^3/s^2.
 * Values are obtained from SPICE kernels in http://naif.jpl.nasa.gov/pub/naif/generic_kernels/
 */
#define MU_SUN    132712440023.310
#define MU_MERCURY       22032.080
#define MU_VENUS        324858.599
#define MU_EARTH        398600.436
#define MU_MOON           4902.799
#define MU_MARS          42828.314
#define MU_JUPITER   126712767.881
#define MU_SATURN     37940626.068
#define MU_URANUS      5794559.128
#define MU_NEPTUNE     6836534.065
#define MU_PLUTO           983.055

/* planet information for major solar system bodies. Units are in km.
 * data taken from http://nssdc.gsfc.nasa.gov/planetary/planets.html
 */

/* Sun */
#define REQ_SUN      695000.                /* km */

/* Earth */
#define REQ_EARTH      6378.1366            /* km, from SPICE */
#define RP_EARTH       6356.7519            /* km, from SPICE */
#define J2_EARTH        1082.616e-6
#define J3_EARTH          -2.53881e-6
#define J4_EARTH          -1.65597e-6
#define J5_EARTH          -0.15e-6
#define J6_EARTH           0.57e-6
#define SMA_EARTH         1.00000011*AU
#define I_EARTH           0.00005*D2R
#define E_EARTH           0.01671022
#define OMEGA_EARTH       0.00007292115     /* Earth's planetary rotation rate, rad/sec */

/* Moon */
#define REQ_MOON       1737.4
#define J2_MOON         202.7e-6
#define SMA_MOON          0.3844e6
#define E_MOON            0.0549
#define I_MOON          5.154 * D2R 

// Sun ephemeris parameters from 2017 Astronomical Almanac
#define SUN_ANGLE_FROM_EARTH 0.00436332312998582
#define COS_SUN_ANGLE_FROM_EARTH (cos(SUN_ANGLE_FROM_EARTH))
#define SIN_SUN_ANGLE_FROM_EARTH (sin(SUN_ANGLE_FROM_EARTH))

// Moon ephemeris parameters from David Simpson
#define MOON_EL0    D2R*218.32
#define MOON_EL1    D2R*481267.883
#define MOON_EL_S_0    D2R*+6.29
#define MOON_EL_S_1    D2R*-1.27
#define MOON_EL_S_2    D2R*+0.66
#define MOON_EL_S_3    D2R*+0.21
#define MOON_EL_S_4    D2R*-0.19
#define MOON_EL_S_5    D2R*-0.11
#define MOON_EL_S0_0    D2R*135.9
#define MOON_EL_S0_1    D2R*259.2
#define MOON_EL_S0_2    D2R*235.7
#define MOON_EL_S0_3    D2R*269.9
#define MOON_EL_S0_4    D2R*357.5
#define MOON_EL_S0_5    D2R*186.5
#define MOON_EL_S1_0    D2R*+477198.85
#define MOON_EL_S1_1    D2R*-413335.38
#define MOON_EL_S1_2    D2R*+890534.23
#define MOON_EL_S1_3    D2R*+954397.70
#define MOON_EL_S1_4    D2R*+35999.05
#define MOON_EL_S1_5    D2R*+966404.05
#define MOON_ELA_S_0    D2R*+5.13
#define MOON_ELA_S_1    D2R*+0.28
#define MOON_ELA_S_2    D2R*-0.28
#define MOON_ELA_S_3    D2R*-0.17
#define MOON_ELA_S0_0    D2R*93.3
#define MOON_ELA_S0_1    D2R*228.2
#define MOON_ELA_S0_2    D2R*318.3
#define MOON_ELA_S0_3    D2R*217.6
#define MOON_ELA_S1_0    D2R*+483202.03
#define MOON_ELA_S1_1    D2R*+960400.87
#define MOON_ELA_S1_2    D2R*+6003.18
#define MOON_ELA_S1_3    D2R*-407332.20
#define MOON_HP0    D2R*0.9508
#define MOON_HP_C_0    D2R*+0.0518
#define MOON_HP_C_1    D2R*+0.0095
#define MOON_HP_C_2    D2R*+0.0078
#define MOON_HP_C_3    D2R*+0.0028
#define MOON_HP_C0_0    D2R*134.9
#define MOON_HP_C0_1    D2R*259.2
#define MOON_HP_C0_2    D2R*235.7
#define MOON_HP_C0_3    D2R*269.9
#define MOON_HP_C1_0    D2R*+477198.85
#define MOON_HP_C1_1    D2R*-413335.38
#define MOON_HP_C1_2    D2R*+890534.23
#define MOON_HP_C1_3    D2R*+954397.70
#define MOON_PRCSN_P0_1	0.
#define MOON_PRCSN_P0_2	0.
#define MOON_PRCSN_P0_3	D2R*+5.12362
#define MOON_PRCSN_P1_1	D2R*+1.396971
#define MOON_PRCSN_P1_2	D2R*+0.013056
#define MOON_PRCSN_P1_3	D2R*-1.155358
#define MOON_PRCSN_P2_1	D2R*+0.0003086
#define MOON_PRCSN_P2_2	D2R*-0.0000092
#define MOON_PRCSN_P2_3	D2R*-0.0001964
#define OBLIQUEITY_ANGLE 0.409092804222328 /* radians, from SPICE kernels, 23.439291 degrees */
#define COS_OBLIQUEITY_ANGLE 0.917482137072805
#define SIN_OBLIQUEITY_ANGLE 0.397777155931913

#endif
