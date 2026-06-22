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

void moonEph(double tJ2000, double INPLACE_ARRAY1[3], double INPLACE_ARRAY1[3]);
void sunEph(double tJ2000, double INPLACE_ARRAY1[3], double INPLACE_ARRAY1[3]);