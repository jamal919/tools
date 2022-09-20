
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

\brief FE gauss integration
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



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_quadrature_LGP2xIPG7_2D(int i, int k, gauss_t gauss)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
  sum=gauss.base[LGP2][k][i]*gauss.w[k];
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_quadrature_dxP2xIPG_2D(int i, int k, gauss_t gauss)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
  sum=gauss.base_x[LGP2][k][i]*gauss.w[k];
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_quadrature_dyP2xIPG_2D(int i, int k, gauss_t gauss)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
  sum=gauss.base_y[LGP2][k][i]*gauss.w[k];
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_quadrature_dxP2xLGP2xIPG7_2D(double *p, double *q, int k, gauss_t gauss)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, l, status;
  double gauss_p[16],gauss_q[16];
  int    gauss_n=16;
  double sum, **base;
  
  sum = 0;

  for (l=0;l<gauss.nnodes;l++) {
    gauss_p[l]=0;
    gauss_q[l]=0;
    for (i=0;i<6;i++) {
      gauss_p[l]+=p[i]*gauss.base_x[LGP2][l][i];
      gauss_q[l]+=q[i]*gauss.base_x[LGP2][l][i];
      }
    }

  for(j = 0; j < 6; j++) {
    for(i = 0; i < 6; i++) {
      sum += p[i] * q[j] * fe_dxP2xLGP2xLGP1_2D[i][j][k];
      }
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> T fe_QIntegraleXXXXxIPG7_2D_template(int i, T *q, gauss_t gauss, int discretisation)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int l, status;
  T sum;
  double **base=gauss.base[discretisation];
  sum = 0;

  for (l=0;l<gauss.nnodes;l++) {
    sum+=gauss.w[l]*q[l]*base[l][i];
    }

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> T fe_QIntegraleLGP2xIPG7_2D_template(int i, T *q, gauss_t gauss)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int l, status;
  T sum;
  double **base=gauss.base[LGP2];
  sum = 0;

  for (l=0;l<gauss.nnodes;l++) {
    sum+=gauss.w[l]*q[l]*base[l][i];
    }

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float fe_QIntegraleLGP2xIPG7_2D(int i, float *q, gauss_t gauss)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  return (fe_QIntegraleLGP2xIPG7_2D_template( i, q, gauss));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_QIntegraleLGP2xIPG7_2D(int i, double *q, gauss_t gauss)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  return (fe_QIntegraleLGP2xIPG7_2D_template( i, q, gauss));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> fe_QIntegraleLGP2xIPG7_2D(int i, complex<double> *q, gauss_t gauss)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  return (fe_QIntegraleLGP2xIPG7_2D_template( i, q, gauss));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gauss_init_Q(int gauss_n, double *gauss_x, double *gauss_y,double *gauss_w)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*----------------------------------------------------------------------
  gauss quadrature on square*/

  int i, j, k, l;
  double a, b, c, d;
  double c0, c1, c2, c3, c4, c5, c6;
  double ri[4], sj[4], ai[4], bj[4];
  double x[10],w[10],chk;
  int npg;
  int n;
  
  switch (gauss_n) {

    case 1:
      gauss_w[0] = 1.0;
      gauss_x[0] = 0.5;
      gauss_y[0] = 0.5;
      break;
      
    case 16:
      x[0]= -0.861136311594053;
      x[1]= -0.339981043584856;
      x[2]= +0.339981043584856;
      x[3]= +0.861136311594053;
      w[0]=  0.347854845137454;
      w[1]=  0.652145154862546;
      w[2]=  0.652145154862546;
      w[3]=  0.347854845137454;
      for(k=0;k<4;k++) {
        x[k]=0.5*(1.+x[k]);
        w[k]=0.5*w[k];
        }
      n=0;
      chk=0.;
      for(l=0;l<4;l++) {
        for(k=0;k<4;k++) {
          gauss_x[n]=x[k];
          gauss_y[n]=x[l];
          gauss_w[n]=w[k]*w[l];
          chk+=gauss_w[n];
          n++;
          }
        }
      break;
    default:
      return(-1);
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gauss_init(int gauss_n, double *gauss_x, double *gauss_y,double *gauss_w)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*----------------------------------------------------------------------
  gauss quadrature on triangle*/
  int i, j, k;
  double a, b, c, d;
  double c0, c1, c2, c3, c4, c5, c6;
  double ri[4], sj[4], ai[4], bj[4];
  int npg;
  int n;

/*
  !****************************************************************
  ! Cette routine cree les tables de coordonnees et de poids des
  ! points de Gauss pour une integration numerique dans des
  ! elements triangles avec :
  !  en entree :
  !     -IPG            : nombre de points d'integration
  !     -IPGKED=MILIEU  : si on a 3 points aux mileux des cotes
  !     -      =INTE    : si on a 3 points a l'interieur
  !  en sortie :
  !     -Gauss_XY           : coordonnees des points d'integration
  !     -Gauss_W           : poids des points d'integration
  !     -IPG            : nombre total de points d'integration
  !****************************************************************
*/

  c3 = 1. / 3.;
  c5 = 1. / 5.;
  c6 = 1. / 6.;

  /*---------------------- 1 point au centre de gravite -------------------*/

  switch (gauss_n) {

    case (1):
      gauss_w[0] = 0.5;
      gauss_x[0] = c3;
      gauss_y[0] = c3;
      break;

  /*------------------------------- 3 points ------------------------------*/

    case (3):
      gauss_w[0] = c6;
      gauss_w[1] = c6;
      gauss_w[2] = c6;
   /*-------------------- 3 points au tiers des medianes -------------------*/

      gauss_x[0] = c6;
      gauss_y[0] = c6;
      gauss_x[1] = 2. * c3;
      gauss_y[1] = c6;
      gauss_x[2] = c6;
      gauss_y[2] = 2. * c3;
      break;

  /*--------------------- 3 points au milieu des cotes --------------------*/

      /*
         IF[ipgked == 'MILIEU'] THEN
         gauss_x[1]=0.5
         gauss_y[1]=0.5
         gauss_x[2]=0.
         gauss_y[2]=0.5
         gauss_x[3]=0.5
         gauss_y[3]=0.
       */
  /*------- 4 points [au tiers des medianes et a leur intersection] -------*/

    case (4):
      c0 = -0.28125E+0;
      c1 = 0.260417E+0;
      gauss_w[0] = c0;
      gauss_w[1] = c1;
      gauss_w[2] = c1;
      gauss_w[3] = c1;
      gauss_x[0] = c3;
      gauss_y[0] = c3;
      gauss_x[1] = c5;
      gauss_y[1] = c5;
      gauss_x[2] = 3. * c5;
      gauss_y[2] = c5;
      gauss_x[3] = c5;
      gauss_y[3] = 3. * c5;
      break;

    case (6):
      c0 = 0.111691E+0;
      c1 = 0.054876E+0;
      c2 = 0.445948E+0;
      c4 = 0.091576E+0;
      gauss_w[0] = c0;
      gauss_w[1] = c0;
      gauss_w[2] = c0;
      gauss_w[3] = c1;
      gauss_w[4] = c1;
      gauss_w[5] = c1;
      gauss_x[0] = c2;
      gauss_y[0] = c2;
      gauss_x[1] = 1. - 2. * c2;
      gauss_y[1] = c2;
      gauss_x[2] = c2;
      gauss_y[2] = 1. - 2. * c2;
      gauss_x[3] = c4;
      gauss_y[3] = c4;
      gauss_x[4] = 1. - 2. * c4;
      gauss_y[4] = c4;
      gauss_x[5] = c4;
      gauss_y[5] = 1. - 2. * c4;
      break;

 /*------------------------------- 7 points ------------------------------*/

    case (7):
      c0 = (155. + sqrt(15.)) / 2400.;
      c1 = 31. / 240. - c0;
      c2 = (6. + sqrt(15.)) / 21.;
      c4 = 4. / 7. - c2;
      gauss_w[0] = 9. / 80.;
      gauss_w[1] = c0;
      gauss_w[2] = c0;
      gauss_w[3] = c0;
      gauss_w[4] = c1;
      gauss_w[5] = c1;
      gauss_w[6] = c1;
      gauss_x[0] = c3;
      gauss_y[0] = c3;
      gauss_x[1] = c2;
      gauss_y[1] = c2;
      gauss_x[2] = 1. - 2. * c2;
      gauss_y[2] = c2;
      gauss_x[3] = c2;
      gauss_y[3] = 1. - 2. * c2;
      gauss_x[4] = c4;
      gauss_y[4] = c4;
      gauss_x[5] = 1. - 2. * c4;
      gauss_y[5] = c4;
      gauss_x[6] = c4;
      gauss_y[6] = 1. - 2. * c4;
      break;

    /*------------------------------ 12 points ------------------------------*/
    case (12):
      c0 = 0.025422453185103;
      c1 = 0.058393137863189;
      c2 = 0.041425537809187;
      a = 0.063089014491502;
      b = 0.249286745170910;
      c = 0.310352451033785;
      d = 0.053145049844816;
      gauss_w[0] = c0;
      gauss_w[1] = c0;
      gauss_w[2] = c0;
      gauss_w[3] = c1;
      gauss_w[4] = c1;
      gauss_w[5] = c1;
      gauss_w[6] = c2;
      gauss_w[7] = c2;
      gauss_w[8] = c2;
      gauss_w[9] = c2;
      gauss_w[10] = c2;
      gauss_w[11] = c2;
      gauss_x[0] = a;
      gauss_y[0] = a;
      gauss_x[1] = 1 - 2 * a;
      gauss_y[1] = a;
      gauss_x[2] = a;
      gauss_y[2] = 1 - 2 * a;
      gauss_x[3] = b;
      gauss_y[3] = b;
      gauss_x[4] = 1 - 2 * b;
      gauss_y[4] = b;
      gauss_x[5] = b;
      gauss_y[5] = 1 - 2 * b;
      gauss_x[6] = c;
      gauss_y[6] = d;
      gauss_x[7] = d;
      gauss_y[7] = c;
      gauss_x[8] = 1 - (c + d);
      gauss_y[8] = c;
      gauss_x[9] = 1 - (c + d);
      gauss_y[9] = d;
      gauss_x[10] = c;
      gauss_y[10] = 1 - (c + d);
      gauss_x[11] = d;
      gauss_y[11] = 1 - (c + d);
      break;

/*------------------------------ 16 points ------------------------------*/

    case (16):
      ri[0] = 0.0694318422;
      ri[1] = 0.3300094782;
      ri[2] = 0.6699905218;
      ri[3] = 0.9305681558;
      sj[0] = 0.0571041961;
      sj[1] = 0.2768430136;
      sj[2] = 0.5835904324;
      sj[3] = 0.8602401357;
      ai[0] = 0.1739274226;
      ai[1] = 0.3260725774;
      ai[2] = 0.3260725774;
      ai[3] = 0.1739274226;
      bj[0] = 0.1355069134;
      bj[1] = 0.2034645680;
      bj[2] = 0.1298475476;
      bj[3] = 0.0311809709;
      for(i = 0; i < 4; i++)
        for(j = 0; j < 4; j++) {
          gauss_w[j + i * 4] = ai[i] * bj[j];
          gauss_x[j + i * 4] = sj[j];
          gauss_y[j + i * 4] = ri[i] * (1.0 - sj[j]);
        }
      break;

    default:
      return(-1);
    }
  return(0);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   void gauss_init_R4(int gauss_n, float *gauss_x, float *gauss_y, float *gauss_w)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
// /*----------------------------------------------------------------------*/
//   int i, j, k;
//   float a, b, c, d;
//   float c0, c1, c2, c3, c4, c5, c6;
//   float ri[4], sj[4], ai[4], bj[4];
//   int npg;
//   int n;
// /*----------------------------------------------------------------------*/
//
// /*
//   !****************************************************************
//   ! Cette routine cree les tables de coordonnees et de poids des
//   ! points de Gauss pour une integration numerique dans des
//   ! elements triangles avec :
//   !  en entree :
//   !     -IPG            : nombre de points d'integration
//   !     -IPGKED=MILIEU  : si on a 3 points aux mileux des cotes
//   !     -      =INTE    : si on a 3 points a l'interieur
//   !  en sortie :
//   !     -Gauss_XY           : coordonnees des points d'integration
//   !     -Gauss_W           : poids des points d'integration
//   !     -IPG            : nombre total de points d'integration
//   !****************************************************************
// */
//
//   c3 = 1. / 3.;
//   c5 = 1. / 5.;
//   c6 = 1. / 6.;
//
//   /*---------------------- 1 point au centre de gravite -------------------*/
//
//   switch (gauss_n) {
//
//     case (1):
//       gauss_w[0] = 0.5;
//       gauss_x[0] = c3;
//       gauss_y[0] = c3;
//       break;
//
//   /*------------------------------- 3 points ------------------------------*/
//
//     case (3):
//       gauss_w[0] = c6;
//       gauss_w[1] = c6;
//       gauss_w[2] = c6;
//    /*-------------------- 3 points au tiers des medianes -------------------*/
//
//       gauss_x[0] = c6;
//       gauss_y[0] = c6;
//       gauss_x[1] = 2. * c3;
//       gauss_y[1] = c6;
//       gauss_x[2] = c6;
//       gauss_y[2] = 2. * c3;
//       break;
//
//   /*--------------------- 3 points au milieu des cotes --------------------*/
//
//       /*
//          IF[ipgked == 'MILIEU'] THEN
//          gauss_x[1]=0.5
//          gauss_y[1]=0.5
//          gauss_x[2]=0.
//          gauss_y[2]=0.5
//          gauss_x[3]=0.5
//          gauss_y[3]=0.
//        */
//   /*------- 4 points [au tiers des medianes et a leur intersection] -------*/
//
//     case (4):
//       c0 = -0.28125E+0;
//       c1 = 0.260417E+0;
//       gauss_w[0] = c0;
//       gauss_w[1] = c1;
//       gauss_w[2] = c1;
//       gauss_w[3] = c1;
//       gauss_x[0] = c3;
//       gauss_y[0] = c3;
//       gauss_x[1] = c5;
//       gauss_y[1] = c5;
//       gauss_x[2] = 3. * c5;
//       gauss_y[2] = c5;
//       gauss_x[3] = c5;
//       gauss_y[3] = 3. * c5;
//       break;
//
//     case (6):
//       c0 = 0.111691E+0;
//       c1 = 0.054876E+0;
//       c2 = 0.445948E+0;
//       c4 = 0.091576E+0;
//       gauss_w[0] = c0;
//       gauss_w[1] = c0;
//       gauss_w[2] = c0;
//       gauss_w[3] = c1;
//       gauss_w[4] = c1;
//       gauss_w[5] = c1;
//       gauss_x[0] = c2;
//       gauss_y[0] = c2;
//       gauss_x[1] = 1. - 2. * c2;
//       gauss_y[1] = c2;
//       gauss_x[2] = c2;
//       gauss_y[2] = 1. - 2. * c2;
//       gauss_x[3] = c4;
//       gauss_y[3] = c4;
//       gauss_x[4] = 1. - 2. * c4;
//       gauss_y[4] = c4;
//       gauss_x[5] = c4;
//       gauss_y[5] = 1. - 2. * c4;
//       break;
//
//  /*------------------------------- 7 points ------------------------------*/
//
//     case (7):
//       c0 = (155. + sqrt(15.)) / 2400.;
//       c1 = 31. / 240. - c0;
//       c2 = (6. + sqrt(15.)) / 21.;
//       c4 = 4. / 7. - c2;
//       gauss_w[0] = 9. / 80.;
//       gauss_w[1] = c0;
//       gauss_w[2] = c0;
//       gauss_w[3] = c0;
//       gauss_w[4] = c1;
//       gauss_w[5] = c1;
//       gauss_w[6] = c1;
//       gauss_x[0] = c3;
//       gauss_y[0] = c3;
//       gauss_x[1] = c2;
//       gauss_y[1] = c2;
//       gauss_x[2] = 1. - 2. * c2;
//       gauss_y[2] = c2;
//       gauss_x[3] = c2;
//       gauss_y[3] = 1. - 2. * c2;
//       gauss_x[4] = c4;
//       gauss_y[4] = c4;
//       gauss_x[5] = 1. - 2. * c4;
//       gauss_y[5] = c4;
//       gauss_x[6] = c4;
//       gauss_y[6] = 1. - 2. * c4;
//       break;
//
//     /*------------------------------ 12 points ------------------------------*/
//     case (12):
//       c0 = 0.025422453185103;
//       c1 = 0.058393137863189;
//       c2 = 0.041425537809187;
//       a = 0.063089014491502;
//       b = 0.249286745170910;
//       c = 0.310352451033785;
//       d = 0.053145049844816;
//       gauss_w[0] = c0;
//       gauss_w[1] = c0;
//       gauss_w[2] = c0;
//       gauss_w[3] = c1;
//       gauss_w[4] = c1;
//       gauss_w[5] = c1;
//       gauss_w[6] = c2;
//       gauss_w[7] = c2;
//       gauss_w[8] = c2;
//       gauss_w[9] = c2;
//       gauss_w[10] = c2;
//       gauss_w[11] = c2;
//       gauss_x[0] = a;
//       gauss_y[0] = a;
//       gauss_x[1] = 1 - 2 * a;
//       gauss_y[1] = a;
//       gauss_x[2] = a;
//       gauss_y[2] = 1 - 2 * a;
//       gauss_x[3] = b;
//       gauss_y[3] = b;
//       gauss_x[4] = 1 - 2 * b;
//       gauss_y[4] = b;
//       gauss_x[5] = b;
//       gauss_y[5] = 1 - 2 * b;
//       gauss_x[6] = c;
//       gauss_y[6] = d;
//       gauss_x[7] = d;
//       gauss_y[7] = c;
//       gauss_x[8] = 1 - (c + d);
//       gauss_y[8] = c;
//       gauss_x[9] = 1 - (c + d);
//       gauss_y[9] = d;
//       gauss_x[10] = c;
//       gauss_y[10] = 1 - (c + d);
//       gauss_x[11] = d;
//       gauss_y[11] = 1 - (c + d);
//       break;
//
// /*------------------------------ 16 points ------------------------------*/
//
//     case (16):
//       ri[0] = 0.0694318422;
//       ri[1] = 0.3300094782;
//       ri[2] = 0.6699905218;
//       ri[3] = 0.9305681558;
//       sj[0] = 0.0571041961;
//       sj[1] = 0.2768430136;
//       sj[2] = 0.5835904324;
//       sj[3] = 0.8602401357;
//       ai[0] = 0.1739274226;
//       ai[1] = 0.3260725774;
//       ai[2] = 0.3260725774;
//       ai[3] = 0.1739274226;
//       bj[0] = 0.1355069134;
//       bj[1] = 0.2034645680;
//       bj[2] = 0.1298475476;
//       bj[3] = 0.0311809709;
//       for(i = 0; i < 4; i++)
//         for(j = 0; j < 4; j++) {
//           gauss_w[j + i * 4] = ai[i] * bj[j];
//           gauss_x[j + i * 4] = sj[j];
//           gauss_y[j + i * 4] = ri[i] * (1.0 - sj[j]);
//         }
//       break;
//   }
// }
