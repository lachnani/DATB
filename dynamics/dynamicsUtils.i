/* file: dynamicsUtils.i */
%module dynamicsUtils
%{
/* Everything in this block will be copied in the wrapper file. We include the C header file nececounterary to compile the interface
*/
#define SWIG_FILE_WITH_INIT
#include "dynamicsConstants.h"
#include "dynamicsUtils.h"
%}

%include "numpy.i"
%init %{
import_array();
%}

/* typemaps.i allows input and output pointer arguments to be specified using the names INPUT, OUTPUT, or INOUT */
%include "typemaps.i"

/* list functions to be interfaced: */

void    dcmInr2Ric(const double IN_ARRAY1[3], const double IN_ARRAY1[3], double INPLACE_ARRAY2[3][3]);
void    rv2ric(const double IN_ARRAY1[3], const double IN_ARRAY1[3], const double IN_ARRAY1[3], const double IN_ARRAY1[3],
        double INPLACE_ARRAY1[3], double INPLACE_ARRAY1[3]);
void    ric2rv(const double IN_ARRAY1[3], const double IN_ARRAY1[3], const double IN_ARRAY1[3], const double IN_ARRAY1[3],
        double INPLACE_ARRAY1[3], double INPLACE_ARRAY1[3]);
void    curvRic2rectRic(const double IN_ARRAY1[3], const double IN_ARRAY1[3],
        const double IN_ARRAY1[3], const double IN_ARRAY1[3],
        double INPLACE_ARRAY1[3], double INPLACE_ARRAY1[3]);
void    rectRic2curvRic(const double IN_ARRAY1[3], const double IN_ARRAY1[3],
        const double IN_ARRAY1[3], const double IN_ARRAY1[3],
        double INPLACE_ARRAY1[3], double INPLACE_ARRAY1[3]);
void    clroe2ric(const double IN_ARRAY1[6], double n, double t, double INPLACE_ARRAY1[3], double INPLACE_ARRAY1[3]);
void    ric2clroe(const double IN_ARRAY1[3], const double IN_ARRAY1[3], double n, double t, double INPLACE_ARRAY1[6]);
void    oe2dAmico(const double IN_ARRAY1[6], const double IN_ARRAY1[6], double INPLACE_ARRAY1[6]);
void    measParams(const double IN_ARRAY1[3], const double IN_ARRAY1[3], double *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT);

void    Orbit_rk4(
    const int solarGrav,
    const int lunarGrav,
    const int drag,
    const int jnum,
    const double IN_ARRAY1[3],
    const double IN_ARRAY1[3],
    const double Cd,
    const double normA,
    const double IN_ARRAY1[3],
    double dt,
    double INPLACE_ARRAY1[3],
    double INPLACE_ARRAY1[3]);

int		eclipse(const double IN_ARRAY1[3], const double IN_ARRAY1[3]);
void	grav(const double IN_ARRAY1[3], const int body, const double IN_ARRAY1[3], double INPLACE_ARRAY1[3]);
void    jPerturb(const double IN_ARRAY1[3], int num, const int body, double INPLACE_ARRAY1[3]);

void    oe2rv(double mu, const double IN_ARRAY1[6], double INPLACE_ARRAY1[3], double INPLACE_ARRAY1[3]);
void    rv2oe(double mu, const double IN_ARRAY1[3], const double IN_ARRAY1[3], double INPLACE_ARRAY1[6]);
void	oe2ee(const double IN_ARRAY1[6], double INPLACE_ARRAY1[6]);
void    ee2rv(double mu, const double IN_ARRAY1[6], double INPLACE_ARRAY1[3], double INPLACE_ARRAY1[3]);
void    rv2ee(double mu, const double IN_ARRAY1[3], const double IN_ARRAY1[3], double INPLACE_ARRAY1[6]);

double  E2f(double E, double e);
double  E2M(double E, double e);
double  f2E(double f, double e);
double  M2E(double M, double e);

void moonEph(double tJ2000, double INPLACE_ARRAY1[3], double INPLACE_ARRAY1[3]);
void sunEph(double tJ2000, double INPLACE_ARRAY1[3], double INPLACE_ARRAY1[3]);
void    envAngles(const double IN_ARRAY1[3], const double IN_ARRAY1[3], const double IN_ARRAY1[3], const double IN_ARRAY1[3], 
        double *OUTPUT, double *OUTPUT, double *OUTPUT);