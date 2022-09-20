
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief declarations of usefull constants
*/
/*----------------------------------------------------------------------------*/


#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <math.h> // for M_PI
#include <stddef.h> //for size_t


/**
\todo 2012-01-12 Damien Allain : clean up redundant uses of the following constants
\note lib_proj.h does define these constants, but to a lower precision ...
*/
const double r2d = 180. / M_PI;///<multiply by this to convert radians to degrees
const double d2r = M_PI / 180.;///<multiply by this to convert degrees to radians
const double M_PIx2 = 2. * M_PI;

/// Von Karman's constant.
const double kappa = 0.4;

/// Minimum horizontal speed authorized.
const double u0_min = 0.01;

/// Classical value of the drag coefficient.
const double z0_def = 2e-3;

/// Square of von Karman's constant.
const double kappa2 = kappa * kappa;

/// Constant that has the same role as von Karman's constant.
const double K = 0.1;

/// Square of K.
const double K2 = K * K;

/// K pow 4.
const double K4 = K2 * K2;

/// Rotation rate of Earth.
const double Omega = 7.292e-5;

/// Value used in computations.
const double sixteenOmega2 = 16. * Omega * Omega;

/// Radius of a sphere having the same volume as the Earth.
const double MeanEarthRadius=6371e3;

/// square of above
const double MeanEarthRadius2=MeanEarthRadius*MeanEarthRadius;

/// degree to meter convertion
const double d2m=d2r*MeanEarthRadius;

/// meter to degree convertion
const double m2d=r2d/MeanEarthRadius;

/// \todo The following two constants MUST become parameters.

/// Asked precision (for instance in Newton's method).
const double EPSILON = 1e-6;
/// Maximum number of iteration authorized in Newton's method.
const size_t MAXIT = 2;

#ifdef MAIN_SOURCE
#define GLOBAL
#define INITIALVALUE(x) =x
#else
#define GLOBAL extern
#define INITIALVALUE(x)
#endif

//double love_h2 INITIALVALUE(0.6), love_k2 INITIALVALUE(0.3);
GLOBAL double love_h2 INITIALVALUE(0.6032), love_k2 INITIALVALUE(0.2980); // D. C. Agnew, University of California San Diego, San Diego, CA, USA
GLOBAL double linear_h INITIALVALUE(.1), explicit_loading INITIALVALUE(1.0);
GLOBAL double CheckInterval INITIALVALUE(3600.0);
GLOBAL double R INITIALVALUE(6.3675e6);//only used in (s)tugo, thankfully

GLOBAL double P_g INITIALVALUE(9.80665);

GLOBAL double rho_water INITIALVALUE(1026.);

GLOBAL double tau_r INITIALVALUE(600.0);
GLOBAL int ContinueCount INITIALVALUE(0);
GLOBAL int ContinueCount3D INITIALVALUE(0);

/* NOTE: DO NOT ENABLE CODE BELOW AS IT WILL BREAK POCViP !!!
#undef GLOBAL
#undef INITIALVALUE
*/

#endif // CONSTANTS_H
