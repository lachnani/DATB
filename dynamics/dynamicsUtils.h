// file: dynamicsUtils.h

#define CELESTIAL_EARTH 0
#define CELESTIAL_MOON  1
#define CELESTIAL_SUN   2

//Formation
void    dcmInr2Ric(const double r[3], const double v[3], double RN[3][3]);
void    rv2ric(const double r[3], const double v[3], const double r_d[3], const double v_d[3],
        double relPosRectRic[3], double relVelRectRic[3]);
void    ric2rv(const double r[3], const double v[3], const double relPosRectRic[3], const double relVelRectRic[3],
        double r_d[3], double v_d[3]);
void    curvRic2rectRic(const double r[3], const double v[3],
        const double relPosCurvRic[3], const double relVelCurvRic[3],
        double relPosRectRic[3], double relVelRectRic[3]);
void    rectRic2curvRic(const double r[3], const double v[3],
        const double relPosRectRic[3], const double relVelRectRic[3],
        double relPosCurvRic[3], double relVelCurvRic[3]);
void    clroe2ric(const double L[6], double n, double t, double relPosRic[3], double relVelRic[3]);
void    ric2clroe(const double relPosRic[3], const double relVelRic[3], double n, double t, double L[6]);
void    oe2dAmico(const double oeChief[6], const double oeDeputy[6], double dAmico[6]);
void    measParams(const double relPosRic[3], const double relVelRic[3], double *rng, double *rngRate, double *az, double *el);

//Orbit 
void    acceleration(
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
    double a[3]);

void     stateUpdate(
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
    double yDot[6]);

void    Orbit_rk4(
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
    double v[3]);

int		eclipse(const double rSc[3], const double sun_rUnit[3]);
void	grav(const double r[3], const int body, const double rBody[3], double agrav[3]);
void    jPerturb(const double r[3], int num, const int body, double ajtot[3]);
void    dragPerturb(const double Cd, const double normA, const double r[3], const double v[3], double ad[3]);

void    oe2rv(double mu, const double oe[6], double r[3], double v[3]);
void    rv2oe(double mu, const double r[3], const double v[3], double oe[6]);
void	oe2ee(const double oe[6], double ee[6]);
void    ee2rv(double mu, const double ee[6], double r[3], double v[3]);
void    rv2ee(double mu, const double r[3], const double v[3], double ee[6]);

double  E2f(double E, double e);
double  E2M(double E, double e);
double  f2E(double f, double e);
double  M2E(double M, double e);

//Ephemerides
void    moonEph(double tJ2000, double rUnit[3], double r[3]);
void    sunEph(double tJ2000, double rUnit[3], double r[3]);
void    envAngles(const double rChief[3], const double rDeputy[3], const double rMoon[3], const double rSun[3], 
        double *losEarthAng, double *losMoonAng, double *losSunAng);