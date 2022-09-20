

/**************************************************************************

  POC-SOLVERS interface, 2006-2012

  Part of the Unstructured Ocean Grid initiative and SYMPHONIE suite

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      PhD, LEGOS, Toulouse, France

***************************************************************************/

#include <string.h>
#include <omp.h>
#include <complex.h>
#include <vector>

#include "solvers-functions.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void ptr_delete_void(void* & p, void (*F)(void *))

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  delete a void pointer and set it to NULL
------------------------------------------------------------------------------*/
{
  if(!p) return;
  
  F(p);
  
  p=0;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void deletep_void(void **p, void (*F)(void *))

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  delete a void pointer and set it to NULL
------------------------------------------------------------------------------*/
{
  if(!*p) return;
  
  F(*p);
  
  *p=NULL;
}





