/**
 * \file friction-additional.cpp
 * \author Yoann LE BARS
 * \version 1.0
 * \date December 15th, 2006
 * \date December 20th, 2006
 * \date August 11th, 2008
 *
 * Additional function used for friction with Karman condition and using a
 * zonal z0.
 */

// see 2007.09.26 version

/**************************************************************************

  Nonlinear finite element time stepping model

***************************************************************************/

#include <assert.h>

#include "tugo-prototypes.h"
#include "constants.h"
#include "drag.hpp"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void karman_Cd(parameter_t  data, double *H_u, int nnodes, float *minCd)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------
     coefficient de CHEZY variable en fonction de la profondeur 
     d'apres ##soulsby 93## Coastal engineering vol 21 1993 p41--69 
      formule (12) page 52 
     Cd=[ 0.40 / (ln(h/z0))-1 ] 
      thierry le 16May2006

     Yoann, Decembre 15th, 2006 : total rewriting, see internal report for
                                  details.

    Decembre 20th, 2006 : computing He using Newton method. Nowadays, EPSILON
			  and MAXIT are constants, they have to become
			  parameters given by the user.
------------------------------------------------------------------------*/

  /**
   * \todo EPSILON and MAXIT are constants, they have to become parameters
   * given by users.
   */

  for (int n = 0 ; n < nnodes ; n ++) {
/**----------------------------------------------------------------------------
    rugosity length */
    double z0 = data.z0[n];
#ifdef DRY
    const double H = max(H_u[n], z0 + EPSILON) ;
#else
    const double H = H_u[n];      /// HERE
//    const double H = H_u[n]/2 ; /// HERE
    if (z0 >= H) {
      data.Cd[n] = 0.;
      continue;
      }
#endif

    if (z0 == 0.) {  
/**----------------------------------------------------------------------------
      Viscous case */
      data.FrictionRatio[n] = 0.0;
      z0=gz0_default;
      }
    else {
      data.FrictionRatio[n] = 1.0;
      }
/**----------------------------------------------------------------------------
    Constant for the computation, see intenal report. */
    const double l = log(1e4) - 1.;
    
/**----------------------------------------------------------------------------
    Temporal value for computing. */
    const double temp = (kappa * (H - z0)) / (H * log(H / z0) + z0 - H);
     
/**----------------------------------------------------------------------------
    Value of the friction condition.*/
    const double Cd = MAX(minCd[n], temp * temp);
    data.Cd[n] = Cd;
//     data.Cd[n] = 0.0;   /// HERE !!!
    }
}
