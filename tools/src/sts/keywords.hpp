/*****************************************************************************
*                        C O N S T A N T E S . H P P                         *
* -------------------------------------------------------------------------- *
* Header file which contains all keywords either used by tugo and tugogui.   *
* -------------------------------------------------------------------------- *
* Created     : december 5th, 2006                                           *
* Last update : december 6th, 2006                                           *
* Version     : 1.5                                                          *
* Contribuers : Yoann LE BARS                                                *
*****************************************************************************/

#ifndef KEYWORDS_HPP
#define KEYWORDS_HPP

/**
 * \file keywords.hpp
 * \brief Header file which contains all keywords either used by tugo and
 * tugogui.
 * \author Yoann LE BARS
 * \version 1.5
 * \date December 5th, 2006
 * \date December 6th, 2006
 * \date August 13th, 2008
 */

#include <string>

/**
 * \namespace Keywords
 * \brief Namespace for every keywords.
 */
namespace Keywords {
  using std::string;

  /* Constantes for the configuration file. */

  // bool constants

  const string KEY_TRUE  = "TRUE" ;
  const string KEY_FALSE = "FALSE" ;

  // Type of elements

  const string KEY_INTRINSIC    = "INTRINSIC" ;
  const string KEY_AUTO         = "AUTO" ;
  const string KEY_LGP0         = "LGP0" ;
  const string KEY_LGP1         = "LGP1" ;
  const string KEY_DGP1         = "DGP1" ;
  const string KEY_LGP2         = "LGP2" ;
  const string KEY_NCP1         = "NCP1" ;

  const string KEY_LGP1_CG      = "LGP1-CG" ;
  const string KEY_LGP1_DG      = "LGP1-DG" ;
  const string KEY_LGP2_CG      = "LGP2-CG" ;
  const string KEY_NCP1_CG      = "NCP1-CG" ;
  const string KEY_NCP1_DG      = "NCP1-DG" ;
  
  const string KEY_CQP0         = "CQP0" ;
  const string KEY_CQP1         = "CQP1" ;
  const string KEY_CQN1         = "CQN1" ;


  const string KEY_P1P1         = "P1-P1" ;
  const string KEY_P1NCP1       = "P1-NCP1" ;
  const string KEY_P1NCP1FAKE   = "P1-NCP1-fake" ;

  const string KEY_GWE_CLX      = "GWE-CLX" ;
  const string KEY_GWE_XXX      = "GWE-XXX" ;
  const string KEY_WE_EXPLICIT  = "WE-EXPLICIT" ;
  const string KEY_WE_IMPLICIT  = "WE-IMPLICIT" ;
  const string KEY_FV_IMPLICIT  = "FV-IMPLICIT" ;

  const string KEY_NQUAD        = "QUADRATURE" ;
  const string KEY_INTGL        = "INTEGRALE" ;

  const string KEY_FORWARD      = "FORWARD" ;
  const string KEY_LEAPFROG     = "LEAPFROG" ;
  const string KEY_UNSPECIFIED  = "UNSPECIFIED" ;

  // Type of rugosity

  const string KEY_LINEAR    = "LINEAR" ;
  const string KEY_QUADRATIC = "QUADRATIC" ;
  const string KEY_KARMAN    = "KARMAN" ;

  const int FRICTION_LINEAR    = 0 ;
  const int FRICTION_QUADRATIC = 1 ;
  const int FRICTION_KARMAN    = 2 ;

  // Utilitary keyword.

  const string KEY_DEFAULT = "<default> " ;
  const string KEY_NONE = "NONE" ;

  const string KEY_MKS = "MKS" ;
  const string KEY_CGS = "CGS" ;
}

#endif
