// file: dynamicsUtils.h

#define CELESTIAL_EARTH 0
#define CELESTIAL_MOON  1
#define CELESTIAL_SUN   2

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

//Ephemerides
void    moonEph(double tJ2000, double rUnit[3], double r[3]);
void    sunEph(double tJ2000, double rUnit[3], double r[3]);