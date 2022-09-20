
/**************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Cyril Nguyen       LA/CNRS,    Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief FE quadrangle integration
*/
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <cmath>

//#define MAIN_FE_INTEGRATION_SOURCE

#include "fe.h"
#include "geo.h"
#include "gauss.h"

#ifndef M_PI
#   define M_PI 3.14159265358979323844
#endif

#define OPTIMAL



// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   template <typename T> T fe_integrale_QP1xQP1_2D_template(int i, int j)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   T sum;
//   sum = 0;
//
//   sum = fe_QP1xQP1_2D[i][j];
//
//   return (sum);
// }
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   double fe_integrale_QP1xQP1_2D(int i, int j)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   double sum;
//   sum = fe_integrale_QP1xQP1_2D_template( i,  j);
//
//   return (sum);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> T fe_integrale_QP1xQP1_2D_template(T* p, T* q)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T sum;
  sum = 0;

  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      sum += p[i] * q[j] * fe_QP1xQP1_2D[i][j];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_integrale_QP1xQP1_2D(double *p, double *q)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
  sum = fe_integrale_QP1xQP1_2D_template( p,  q);
  
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> T fe_integrale_QP1xQP1xQP1_2D_template(int i, int j, T *r)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T sum;
  sum = 0;

  for(int k = 0; k < 4; k++)
    sum += r[k] * fe_QP1xQP1xQP1_2D[i][j][k];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_integrale_QP1xQP1xQP1_2D(int i, int j, double *r)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
  sum = fe_integrale_QP1xQP1xQP1_2D_template( i,  j, r);
  
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> T fe_integrale_QP1xQP1xQP1_2D_template(int i, T* q, T *r)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T sum;
  sum = 0;

  for(int j = 0; j < 4; j++)
    for(int k = 0; k < 4; k++)
      sum += q[j] * r[k] * fe_QP1xQP1xQP1_2D[i][j][k];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_integrale_QP1xQP1xQP1_2D(int i, double *q, double *r)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
  sum = fe_integrale_QP1xQP1xQP1_2D_template( i,  q, r);
  
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> T fe_integrale_QP1xQP1xQP1_2D_template(T* p, T* q, T *r)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T sum;
  sum = 0;

  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      for(int k = 0; k < 4; k++)
        sum += p[i] * q[j] * r[k] * fe_QP1xQP1xQP1_2D[i][j][k];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_integrale_QP1xQP1xQP1_2D( double *p, double *q, double *r)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
  sum = fe_integrale_QP1xQP1xQP1_2D_template( p,  q, r);
  
  return (sum);
}


