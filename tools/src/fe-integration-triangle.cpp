
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
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

\brief FE triangle integration
*/
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>

#define MAIN_FE_INTEGRATION_SOURCE

#include "fe.h"
#include "geo.h"
 
int fe_chk_integrales(FILE *);

double (* fe_sproduct_ptrLGP1xLGP1_2D_a)  (double *, double *);
double (* fe_sproduct_ptrLGP1xLGP1_2D_b)  (double *, int);

double (* fe_sproduct_ptrLGP1xLGP1xLGP1_2D_a) (double *, double *, double *);
double (* fe_sproduct_ptrLGP1xLGP1xLGP1_2D_b) (double *, double *, int);

double (* fe_sproduct_ptrLGP1xLGP1xLGP1xLGP1_2D) (double *, double *, double *, int);

double (* fe_sproduct_ptrLGP1xNCP1_2D)      (double *, double *);
 
double (* fe_sproduct_ptrLGP1xLGP1xNCP1_2D_a)  (double *, double *, double *);
double (* fe_sproduct_ptrLGP1xLGP1xNCP1_2D_b)  (double *, int     , double *);

//double (* fe_sproductLGP1xLGP1_2D)  (double *, double *);
//double (* fe_sproductLGP1xLGP1xLGP1_2D) (double *, double *, double *);
//double (* fe_sproductLGP1xLGP1xLGP1_2D) (double *, double *, int);

//template <typename T> double (* fe_sproductLGP1xLGP1xLGP1_2D) (double *, double *, T );

//double (* fe_sproductLGP1xLGP1xLGP1_2D) (double *, double *, int);

//double (* fe_sproductLGP1xNCP1_2D)      (double *, double *);
//double (* fe_sproductLGP1xLGP1xNCP1_2D) (double *, double *, double *);

/**-----------------------------------------------------------------------
2D nodal coordinates*/
double fe_LGP1_x[3]={0.0, 1.0, 0.0};
double fe_LGP1_y[3]={0.0, 0.0, 1.0};

double fe_NCP1_x[3]={0.5, 0.0, 0.5};
double fe_NCP1_y[3]={0.5, 1.0, 0.0};

double fe_LGP2_x[6]={0.0, 0.5, 1.0, 0.5, 0.0, 0.0};
double fe_LGP2_y[6]={0.0, 0.0, 0.0, 0.5, 1.0, 0.5};

/*------------------------------------------------------------------------------
3D integration coefficients*/
// double fe_LGP1_3D [6];
// double fe_LGP1xLGP1_3D [6][6];
// double fe_LGP1xLGP1xLGP1_3D [6][6][6];

const double fe_LGP1xLGP1xNCP1_2D_c1= -0.016666666666667;
const double fe_LGP1xLGP1xNCP1_2D_c2=  0.008333333333333;
const double fe_LGP1xLGP1xNCP1_2D_c3=  0.050000000000000;
const double fe_LGP1xLGP1xNCP1_2D_c4=  0.025000000000000;

polynomial_t baseLGP2_x_2D[6];
polynomial_t baseLGP2_y_2D[6];

//#define NOMINAL
#define OPTIMAL
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorial(int val)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(val > 1)
    return factorial(val - 1) * val;
  return 1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int polyinit(monomial2D_t * canonic, int order, polynomial_t * p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l;

  p->order = order;
  switch (p->order) {
    case 0:
      p->n = 1;
      break;

    case 1:
      p->n = 3;
      break;

    case 2:
      p->n = 6;
      break;

    case 3:
      p->n = 10;
      break;

    case 4:
      p->n = 15;
      break;

    case 5:
      p->n = 21;
      break;

    case 6:
      p->n = 28;
      break;

    default:
      check_error(-1, "illicit order (not implemented)", __LINE__, __FILE__, 1);
      break;
    }

  p->m = new monomial2D_t *[p->n];
  p->c = new double *[p->order + 1];
  for(l = 0; l < p->order + 1; l++)
    p->c[l] = new double[p->order + 1];

  for(k = 0; k < p->order + 1; k++)
    for(l = 0; l < p->order + 1; l++) {
      p->c[k][l] = 0;
      }
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int polyinit3D(int order, polynomial3D_t * p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j, k, l;

  p->order = order;
  switch (p->order) {
    case 0:
      p->n = 1;
      break;

    case 1:
      p->n = 4;
      break;

    case 2:
      p->n = 10;
      break;

    case 3:
      p->n = 20;
      break;

    case 4:
      p->n = 35;
      break;

    case 5:
      p->n = 56;
      break;

    case 6:
      p->n = 84;
      break;

    default:
      check_error(-1, "illicit order (not implemented)", __LINE__, __FILE__, 1);
      break;
    }

//  p->m = new monomial3D_t *[p->n];
  p->m = 0;
  p->c = new double **[p->order + 1];
  for(l = 0; l < p->order + 1; l++) {
    p->c[l] = new double *[p->order + 1];
    for(k = 0; k < p->order + 1; k++) {
      p->c[l][k] = new double [p->order + 1];
      }
    }

  for(j = 0; j < p->order + 1; j++)
    for(k = 0; k < p->order + 1; k++)
      for(l = 0; l < p->order + 1; l++) {
        p->c[j][k][l] = 0;
        }
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int polyfree(polynomial_t *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int l;
  for(l = 0; l < p->order + 1; l++)
     delete[] p->c[l];
  delete[] p->c;
  delete[] p->m;

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int polyfree(polynomial3D_t *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l;
  
  for(l = 0; l < p->order + 1; l++) {
    for(k = 0; k < p->order + 1; k++) {
      delete[] p->c[l][k];
      }
    delete[] p->c[l];
    }
  delete[] p->c;
  
  if(p->m!=0)
    delete[] p->m;

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int polymul(monomial2D_t * canonic, polynomial_t p1, polynomial_t p2,
                  polynomial_t p3)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l, m, n;

  for(k = 0; k < p3.order + 1; k++)
    for(l = 0; l < p3.order + 1; l++) {
      p3.c[k][l] = 0;
      }

  for(k = 0; k < p1.n; k++) {
    for(l = 0; l < p2.n; l++) {
      m = canonic[k].x + canonic[l].x;
      if(m>p3.order) {
        printf("polynomial order overflow: %d %\n",p3.order+1,m);
        }
      n = canonic[k].y + canonic[l].y;
      if(n>p3.order) {
        printf("polynomial order overflow: %d %\n",p3.order+1,n);
        }
      p3.c[m][n] += p1.c[canonic[k].x][canonic[k].y] * p2.c[canonic[l].x][canonic[l].y];
      }
    }

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int polymul3D(monomial3D_t * canonic, polynomial3D_t p1, polynomial3D_t p2,
                  polynomial3D_t p3)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l, m, n, o;

/**----------------------------------------------------------------------
  initialise resulting polynomial*/
  for(k = 0; k < p3.order + 1; k++) {
    for(l = 0; l < p3.order + 1; l++) {
      for(m = 0; m < p3.order + 1; m++) {
        p3.c[k][l][m] = 0;
        }
      }
    }

  for(k = 0; k < p1.n; k++) {
    for(l = 0; l < p2.n; l++) {
/**----------------------------------------------------------------------
      */
      m = canonic[k].x + canonic[l].x;
      n = canonic[k].y + canonic[l].y;
      o = canonic[k].z + canonic[l].z;
      p3.c[m][n][o] += p1.c[canonic[k].x][canonic[k].y][canonic[k].z] * p2.c[canonic[l].x][canonic[l].y][canonic[l].z];
      }
    }

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double polyint1D(monomial2D_t * canonic, polynomial_t p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------
  x-integration along a [0;1] segment*/
  int l, m, n;
  double sum = 0.;

  for(l = 0; l < p.n; l++) {
    m = canonic[l].x;
    n = canonic[l].y;
    sum += p.c[m][n] / double (m + 1);
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double polyint2D(monomial2D_t * canonic, polynomial_t p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------
  (x,y)-integration over a reference triangle*/
  int l, m, n;
  double sum = 0.;

  for(l = 0; l < p.n; l++) {
    m = canonic[l].x;
    n = canonic[l].y;
    sum += p.c[m][n] * factorial(m) * factorial(n) / factorial(m + n + 2);
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double polyint_Q2D(monomial2D_t * canonic, polynomial_t p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------
  (x,y)-integration over a reference quadrangle:
  
       | (1-y) (t2-t1) + y (t3-t4)    (1-y) (p2-p1) + y (p3-p4) |
  J =  |                                                        |
       | (1-x) (t4-t1) + x (t3-t2)    (1-x) (p4-p1) + x (p3-p2) |

  jacobian = det(J)
  
  (t2-t1) (p4-p1) - (t4-t1) (p2-p1)
  
----------------------------------------------------------------------*/
  int l, m, n;
  double sum = 0.;

  for(l = 0; l < p.n; l++) {
    m = canonic[l].x;
    n = canonic[l].y;
    sum += p.c[m][n] / (m + 1) / (n + 1);
    }
  return (sum);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double polyint3D(monomial3D_t * canonic, polynomial3D_t p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------
  (x,y,z)-integration over a reference prism*/
  int /*k,*/ l, m, n, o;
  double value, sum = 0.;
  double test[3];

  for(l = 0; l < 3; l++) test[l]=0.0;

  for(l = 0; l < p.n; l++) {
    m = canonic[l].x;
    n = canonic[l].y;
    o = canonic[l].z;
    if(p.c[m][n][o]!=0) {
      value = p.c[m][n][o] * factorial(m) * factorial(n) / factorial(m + n + 2)/ double (o + 1);
      sum += value;
      test[o]+=p.c[m][n][o] * factorial(m) * factorial(n) / factorial(m + n + 2);
//      printf("%d %d %d %d %lf %lf %lf\n",l,m,n,o,p.c[m][n][o],value,sum);
      }
    }
//   for(k = 0; k < 3; k++) {
//     printf("%d\n",k);
//     sum=0.0;
//     for(l = 0; l < 3; l++) test[l]=0.0;
//     for(l = 0; l < p.n; l++) {
//       m = canonic[l].x;
//       n = canonic[l].y;
//       o = canonic[l].z;
//       if((p.c[m][n][o]!=0) && (o==k)){
//         value = p.c[m][n][o] * factorial(m) * factorial(n) / factorial(m + n + 2)/ double (o + 1);
//         sum += value;
//         test[o]+=p.c[m][n][o] * factorial(m) * factorial(n) / factorial(m + n + 2);
//         printf("%3d %d %d %d %lf %lf %lf\n",l,m,n,o,p.c[m][n][o],value,sum);
//         }
//       }
//     }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_integrale_init(void)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n, status;
//   double sum2x2[3][3];

  FILE *out=NULL;
/*------------------------------------------------------------------------------
  see integration tabulations*/

  polynomial_t p1, p2, p3, p4, p5;

  polynomial_t baseLGP1_1D[3], baseNCP1_1D[3];
  polynomial_t baseLGP1_2D[3], baseNCP1_2D[3];

  polynomial_t baseLGP2_1D[3];
  polynomial_t baseLGP2_2D[6];
  
  polynomial_t baseQP1_2D[4];
  
  polynomial3D_t p3D1, p3D2, p3D3, p3D4, p3D5;

  polynomial3D_t baseLGP1_3D[6];

  monomial2D_t canonic2D[100];
  monomial3D_t canonic3D[100];
/*------------------------------------------------------------------------------
  define the canonical polynomials by prescribing the x and y exponant*/

/**----------------------------------------------------------------------
  degree 0 */
/*------------------------------------------------------------------------------
  p(x,y)=1*/
  canonic2D[0].x = 0;
  canonic2D[0].y = 0;

/**----------------------------------------------------------------------
  degree 1 */
/*------------------------------------------------------------------------------
  p(x,y)=x*/
  canonic2D[1].x = 1;
  canonic2D[1].y = 0;

/*------------------------------------------------------------------------------
  p(x,y)=y*/
  canonic2D[2].x = 0;
  canonic2D[2].y = 1;

/**----------------------------------------------------------------------
  degree 2 */
/*------------------------------------------------------------------------------
  p(x,y)=xy*/
  canonic2D[3].x = 1;
  canonic2D[3].y = 1;

/*------------------------------------------------------------------------------
  p(x,y)=x²*/
  canonic2D[4].x = 2;
  canonic2D[4].y = 0;

/*------------------------------------------------------------------------------
  p(x,y)=y²*/
  canonic2D[5].x = 0;
  canonic2D[5].y = 2;

/**----------------------------------------------------------------------
  degree 3 */
/*------------------------------------------------------------------------------
  p(x,y)=xy²*/
  canonic2D[6].x = 1;
  canonic2D[6].y = 2;

/*------------------------------------------------------------------------------
  p(x,y)=x²y*/
  canonic2D[7].x = 2;
  canonic2D[7].y = 1;

/*------------------------------------------------------------------------------
  p(x,y)=x³*/
  canonic2D[8].x = 3;
  canonic2D[8].y = 0;

/*------------------------------------------------------------------------------
  p(x,y)=y³*/
  canonic2D[9].x = 0;
  canonic2D[9].y = 3;

/**----------------------------------------------------------------------
  degree 4 */
/*------------------------------------------------------------------------------
  p(x,y)=x³y*/
  canonic2D[10].x = 3;
  canonic2D[10].y = 1;

/*------------------------------------------------------------------------------
  p(x,y)=x²y²*/
  canonic2D[11].x = 2;
  canonic2D[11].y = 2;

/*------------------------------------------------------------------------------
  p(x,y)=xy³*/
  canonic2D[12].x = 1;
  canonic2D[12].y = 3;

/*------------------------------------------------------------------------------
  p(x,y)=x⁴*/
  canonic2D[13].x = 4;
  canonic2D[13].y = 0;

/*------------------------------------------------------------------------------
  p(x,y)=y⁴*/
  canonic2D[14].x = 0;
  canonic2D[14].y = 4;

/**----------------------------------------------------------------------
  degree 5 */
/*------------------------------------------------------------------------------
  p(x,y)=x⁴y*/
  canonic2D[15].x = 4;
  canonic2D[15].y = 1;

/*------------------------------------------------------------------------------
  p(x,y)=y⁵*/
  canonic2D[16].x = 3;
  canonic2D[16].y = 2;

/*------------------------------------------------------------------------------
  p(x,y)=y⁵*/
  canonic2D[17].x = 2;
  canonic2D[17].y = 3;

/*------------------------------------------------------------------------------
  p(x,y)=y⁵*/
  canonic2D[18].x = 1;
  canonic2D[18].y = 4;

/*------------------------------------------------------------------------------
  p(x,y)=x⁵*/
  canonic2D[19].x = 5;
  canonic2D[19].y = 0;

/*------------------------------------------------------------------------------
  p(x,y)=y⁵*/
  canonic2D[20].x = 0;
  canonic2D[20].y = 5;

/**----------------------------------------------------------------------
  degree 6 */

  n=21;

/*------------------------------------------------------------------------------
  p(x,y)=x⁵y*/
  canonic2D[n].x = 5;
  canonic2D[n].y = 1;
  n++;

/*------------------------------------------------------------------------------
  p(x,y)=x⁴y²*/
  canonic2D[n].x = 4;
  canonic2D[n].y = 2;
  n++;

/*------------------------------------------------------------------------------
  p(x,y)=x³y³*/
  canonic2D[n].x = 3;
  canonic2D[n].y = 3;
  n++;

/*------------------------------------------------------------------------------
  p(x,y)=x²y⁴*/
  canonic2D[n].x = 2;
  canonic2D[n].y = 4;
  n++;

/*------------------------------------------------------------------------------
  p(x,y)=xy⁵*/
  canonic2D[n].x = 1;
  canonic2D[n].y = 5;
  n++;

/*------------------------------------------------------------------------------
  p(x,y)=x⁶*/
  canonic2D[n].x = 6;
  canonic2D[n].y = 0;
  n++;

/*------------------------------------------------------------------------------
  p(x,y)=y⁶*/
  canonic2D[n].x = 0;
  canonic2D[n].y = 6;
  n++;

/**----------------------------------------------------------------------
  degree 0 */
/*------------------------------------------------------------------------------
  p(x,y,z)=1*/
  canonic3D[0].x = 0;
  canonic3D[0].y = 0;
  canonic3D[0].z = 0;

/**----------------------------------------------------------------------
  degree 1 */
/*------------------------------------------------------------------------------
  p(x,y,z)=x*/
  canonic3D[1].x = 1;
  canonic3D[1].y = 0;
  canonic3D[1].z = 0;

/*------------------------------------------------------------------------------
  p(x,y,z)=y*/
  canonic3D[2].x = 0;
  canonic3D[2].y = 1;
  canonic3D[2].z = 0;

/*------------------------------------------------------------------------------
  p(x,y,z)=z*/
  canonic3D[3].x = 0;
  canonic3D[3].y = 0;
  canonic3D[3].z = 1;

/**----------------------------------------------------------------------
  degree 2 */
/*------------------------------------------------------------------------------
  p(x,y,z)=x² */
  canonic3D[4].x = 2;
  canonic3D[4].y = 0;
  canonic3D[4].z = 0;

/*------------------------------------------------------------------------------
  p(x,y,z)=y² */
  canonic3D[5].x = 0;
  canonic3D[5].y = 2;
  canonic3D[5].z = 0;

/*------------------------------------------------------------------------------
  p(x,y,z)=z² */
  canonic3D[6].x = 0;
  canonic3D[6].y = 0;
  canonic3D[6].z = 2;

/*------------------------------------------------------------------------------
  p(x,y,z)=xy */
  canonic3D[7].x = 1;
  canonic3D[7].y = 1;
  canonic3D[7].z = 0;

/*------------------------------------------------------------------------------
  p(x,y,z)=xz */
  canonic3D[8].x = 1;
  canonic3D[8].y = 0;
  canonic3D[8].z = 1;

/*------------------------------------------------------------------------------
  p(x,y,z)=yz */
  canonic3D[9].x = 0;
  canonic3D[9].y = 1;
  canonic3D[9].z = 1;

/**----------------------------------------------------------------------
  degree 3 */
/*------------------------------------------------------------------------------
  p(x,y,z)=x³ */
  canonic3D[10].x = 3;
  canonic3D[10].y = 0;
  canonic3D[10].z = 0;

/*------------------------------------------------------------------------------
  p(x,y,z)=y³ */
  canonic3D[11].x = 0;
  canonic3D[11].y = 3;
  canonic3D[11].z = 0;

/*------------------------------------------------------------------------------
  p(x,y,z)=z³ */
  canonic3D[12].x = 0;
  canonic3D[12].y = 0;
  canonic3D[12].z = 3;

/*------------------------------------------------------------------------------
  p(x,y,z)=x²y */
  canonic3D[13].x = 2;
  canonic3D[13].y = 1;
  canonic3D[13].z = 0;

/*------------------------------------------------------------------------------
  p(x,y,z)=x²z */
  canonic3D[14].x = 2;
  canonic3D[14].y = 0;
  canonic3D[14].z = 1;

/*------------------------------------------------------------------------------
  p(x,y,z)=xy² */
  canonic3D[15].x = 1;
  canonic3D[15].y = 2;
  canonic3D[15].z = 0;

/*------------------------------------------------------------------------------
  p(x,y,z)= y²z */
  canonic3D[16].x = 0;
  canonic3D[16].y = 2;
  canonic3D[16].z = 1;

/*------------------------------------------------------------------------------
  p(x,y,z)= xz² */
  canonic3D[17].x = 1;
  canonic3D[17].y = 0;
  canonic3D[17].z = 2;

/*------------------------------------------------------------------------------
  p(x,y,z)= yz² */
  canonic3D[18].x = 0;
  canonic3D[18].y = 1;
  canonic3D[18].z = 2;

/*------------------------------------------------------------------------------
  p(x,y,z)= xyz */
  canonic3D[19].x = 1;
  canonic3D[19].y = 1;
  canonic3D[19].z = 1;

/**----------------------------------------------------------------------
  degree 4 */
/*------------------------------------------------------------------------------
  p(x,y,z)=x⁴ */
  canonic3D[20].x = 4;
  canonic3D[20].y = 0;
  canonic3D[20].z = 0;

/*------------------------------------------------------------------------------
  p(x,y,z)=y⁴ */
  canonic3D[21].x = 0;
  canonic3D[21].y = 4;
  canonic3D[21].z = 0;

/*------------------------------------------------------------------------------
  p(x,y,z)=z⁴ */
  canonic3D[22].x = 0;
  canonic3D[22].y = 0;
  canonic3D[22].z = 4;

/*------------------------------------------------------------------------------
  p(x,y,z)=x³y */
  canonic3D[23].x = 3;
  canonic3D[23].y = 1;
  canonic3D[23].z = 0;

/*------------------------------------------------------------------------------
  p(x,y,z)=x³z */
  canonic3D[24].x = 3;
  canonic3D[24].y = 0;
  canonic3D[24].z = 1;

/*------------------------------------------------------------------------------
  p(x,y,z)=xy³ */
  canonic3D[25].x = 1;
  canonic3D[25].y = 3;
  canonic3D[25].z = 0;

/*------------------------------------------------------------------------------
  p(x,y,z)= y³z */
  canonic3D[26].x = 0;
  canonic3D[26].y = 3;
  canonic3D[26].z = 1;

/*------------------------------------------------------------------------------
  p(x,y,z)= xz³ */
  canonic3D[27].x = 1;
  canonic3D[27].y = 0;
  canonic3D[27].z = 3;

/*------------------------------------------------------------------------------
  p(x,y,z)= yz³ */
  canonic3D[28].x = 0;
  canonic3D[28].y = 1;
  canonic3D[28].z = 3;

/*------------------------------------------------------------------------------
  p(x,y,z)= x²y² */
  canonic3D[29].x = 2;
  canonic3D[29].y = 2;
  canonic3D[29].z = 0;

/*------------------------------------------------------------------------------
  p(x,y,z)= x²z² */
  canonic3D[30].x = 2;
  canonic3D[30].y = 0;
  canonic3D[30].z = 2;

/*------------------------------------------------------------------------------
  p(x,y,z)= y²z² */
  canonic3D[31].x = 0;
  canonic3D[31].y = 2;
  canonic3D[31].z = 2;

  n=19;

  m=4;
  for(i=0;i<=m;i++) {
    for(j=0;j<=m-i;j++) {
      k=m-i-j;
      n++;
      canonic3D[n].x = i;
      canonic3D[n].y = j;
      canonic3D[n].z = k;
      }
   }
  m=5;
  for(i=0;i<=m;i++) {
    for(j=0;j<=m-i;j++) {
      k=m-i-j;
      n++;
      canonic3D[n].x = i;
      canonic3D[n].y = j;
      canonic3D[n].z = k;
      }
   }
  m=6;
  for(i=0;i<=m;i++) {
    for(j=0;j<=m-i;j++) {
      k=m-i-j;
      n++;
      canonic3D[n].x = i;
      canonic3D[n].y = j;
      canonic3D[n].z = k;
      }
   }

/*------------------------------------------------------------------------------
  allocate 1D Lagrange-P1 base functions*/
  for(k = 0; k < 2; k++) {
    status = polyinit(canonic2D, 1, &(baseLGP1_1D[k]));
    }

/*------------------------------------------------------------------------------
  allocate 1D non conforming-P1 base functions*/
  for(k = 0; k < 1; k++) {
    status = polyinit(canonic2D, 1, &(baseNCP1_1D[k]));
    }

/*------------------------------------------------------------------------------
  allocate 2D Lagrange-P1 base functions*/
  for(k = 0; k < 3; k++) {
    status = polyinit(canonic2D, 1, &(baseLGP1_2D[k]));
    status = polyinit(canonic2D, 1, &(baseNCP1_2D[k]));
    }

/*------------------------------------------------------------------------------
  allocate 1D Lagrange-P2 base functions*/
  for(k = 0; k < 3; k++) {
    status = polyinit(canonic2D, 2, &(baseLGP2_1D[k]));
    }

/*------------------------------------------------------------------------------
  allocate 2D Lagrange-P2 base functions*/
  for(k = 0; k < 6; k++) {
    status = polyinit(canonic2D, 2, &(baseLGP2_2D[k]));
    status = polyinit(canonic2D, 1, &(baseLGP2_x_2D[k]));
    status = polyinit(canonic2D, 1, &(baseLGP2_y_2D[k]));
    }

  for(k = 0; k < 4; k++) {
    status = polyinit(canonic2D, 2, &(baseQP1_2D[k]));
    }
    
  for(k = 0; k < 6; k++) {
    status = polyinit3D(2, &(baseLGP1_3D[k]));
    }
/*------------------------------------------------------------------------------
  define 1D Lagrange-P1 base functions, coefficients of cannonics */
  baseLGP1_1D[0].c[0][0] =  1; /*  1 */
  baseLGP1_1D[0].c[1][0] = -1; /* -x */

  baseLGP1_1D[1].c[0][0] = 0; /*  0 */
  baseLGP1_1D[1].c[1][0] = 1; /*  x */

/*------------------------------------------------------------------------------
  define 1D non conforming-P1 base functions, coefficients of cannonics */
  baseNCP1_1D[0].c[0][0] = -1; /* -1 */
  baseNCP1_1D[0].c[1][0] =  2; /* 2x */

/*------------------------------------------------------------------------------
  define 2D Lagrange-P1 base functions, coefficients of cannonics */
  baseLGP1_2D[0].c[0][0] =  1; /*  1 */
  baseLGP1_2D[0].c[1][0] = -1; /*  x */
  baseLGP1_2D[0].c[0][1] = -1; /*  y */

  baseLGP1_2D[1].c[0][0] = 0; /*  1 */
  baseLGP1_2D[1].c[1][0] = 1; /*  x */
  baseLGP1_2D[1].c[0][1] = 0; /*  y */

  baseLGP1_2D[2].c[0][0] = 0; /*  1 */
  baseLGP1_2D[2].c[1][0] = 0; /*  x */
  baseLGP1_2D[2].c[0][1] = 1; /*  y */

/*------------------------------------------------------------------------------
  define 2D non conforming-P1 base functions, coefficients of cannonics */
  baseNCP1_2D[0].c[0][0] = -1;
  baseNCP1_2D[0].c[1][0] =  2;
  baseNCP1_2D[0].c[0][1] =  2;

  baseNCP1_2D[1].c[0][0] =  1;
  baseNCP1_2D[1].c[1][0] = -2;
  baseNCP1_2D[1].c[0][1] =  0;

  baseNCP1_2D[2].c[0][0] =  1;
  baseNCP1_2D[2].c[1][0] =  0;
  baseNCP1_2D[2].c[0][1] = -2;
/*------------------------------------------------------------------------------
  define 1D Lagrange-P2 base functions, coefficients of cannonics */
/*
  beta[0] = 2. * (x - 0.5) * (x - 1.)=2(x²+0.5+-x-0.5x)=2x²+1+-3x
  beta[1] = -4. * x * (x - 1.)=-4x²+4x
  beta[2] = 2. * x * (x - 0.5)=2x²-1
*/
  baseLGP2_1D[0].c[0][0] =  1; /*  1 */
  baseLGP2_1D[0].c[1][0] = -3; /*  x */
  baseLGP2_1D[0].c[0][1] =  0; /*  y */
  baseLGP2_1D[0].c[1][1] =  0; /* xy */
  baseLGP2_1D[0].c[2][0] =  2; /* x² */
  baseLGP2_1D[0].c[0][2] =  0; /* y² */

  baseLGP2_1D[1].c[0][0] =  0; /*  1 */
  baseLGP2_1D[1].c[1][0] =  4; /*  x */
  baseLGP2_1D[1].c[0][1] =  0; /*  y */
  baseLGP2_1D[1].c[1][1] =  0; /* xy */
  baseLGP2_1D[1].c[2][0] = -4; /* x² */
  baseLGP2_1D[1].c[0][2] =  0; /* y² */

  baseLGP2_1D[2].c[0][0] = -1; /*  1 */
  baseLGP2_1D[2].c[1][0] =  0; /*  x */
  baseLGP2_1D[2].c[0][1] =  0; /*  y */
  baseLGP2_1D[2].c[1][1] =  0; /* xy */
  baseLGP2_1D[2].c[2][0] =  2; /* x² */
  baseLGP2_1D[2].c[0][2] =  0; /* y² */


/*------------------------------------------------------------------------------
  define 2D Lagrange-P2 base functions, coefficients of cannonics */
//   beta[0] =  2. * (x + y - 0.5) * (x + y - 1.);   /*  2x²+2y²+4xy-3x-3y+1 */
//   beta[1] = -4. * x * (x + y - 1.);               /* -4x²    -4xy+4x      */
//   beta[2] =  2. * x * (x - 0.5);                  /*  2x²        - x      */
//   beta[3] =  4. * x * y;                          /*         +4xy         */
//   beta[4] =  2. * y * (y - 0.5);                  /*     +2y²       - y   */
//   beta[5] = -4. * y * (x + y - 1.);               /*     -4y²-4xy   +4y   */

  baseLGP2_2D[0].c[0][0] =  1; /*    1 */
  baseLGP2_2D[0].c[1][0] = -3; /*  -3x */
  baseLGP2_2D[0].c[0][1] = -3; /*  -3y */
  baseLGP2_2D[0].c[1][1] =  4; /*  4xy */
  baseLGP2_2D[0].c[2][0] =  2; /*  2x² */
  baseLGP2_2D[0].c[0][2] =  2; /*  2y² */

  baseLGP2_2D[1].c[0][0] =  0; /*  1 */
  baseLGP2_2D[1].c[1][0] =  4; /*  x */
  baseLGP2_2D[1].c[0][1] =  0; /*  y */
  baseLGP2_2D[1].c[1][1] = -4; /* xy */
  baseLGP2_2D[1].c[2][0] = -4; /* x² */
  baseLGP2_2D[1].c[0][2] =  0; /* y² */

  baseLGP2_2D[2].c[0][0] =  0; /*  1 */
  baseLGP2_2D[2].c[1][0] = -1; /*  x */
  baseLGP2_2D[2].c[0][1] =  0; /*  y */
  baseLGP2_2D[2].c[1][1] =  0; /* xy */
  baseLGP2_2D[2].c[2][0] =  2; /* x² */
  baseLGP2_2D[2].c[0][2] =  0; /* y² */

  baseLGP2_2D[3].c[0][0] =  0; /*  1 */
  baseLGP2_2D[3].c[1][0] =  0; /*  x */
  baseLGP2_2D[3].c[0][1] =  0; /*  y */
  baseLGP2_2D[3].c[1][1] =  4; /* xy */
  baseLGP2_2D[3].c[2][0] =  0; /* x² */
  baseLGP2_2D[3].c[0][2] =  0; /* y² */

  baseLGP2_2D[4].c[0][0] =  0; /*  1 */
  baseLGP2_2D[4].c[1][0] =  0; /*  x */
  baseLGP2_2D[4].c[0][1] = -1; /*  y */
  baseLGP2_2D[4].c[1][1] =  0; /* xy */
  baseLGP2_2D[4].c[2][0] =  0; /* x² */
  baseLGP2_2D[4].c[0][2] =  2; /* y² */

  baseLGP2_2D[5].c[0][0] =  0; /*  1 */
  baseLGP2_2D[5].c[1][0] =  0; /*  x */
  baseLGP2_2D[5].c[0][1] =  4; /*  y */
  baseLGP2_2D[5].c[1][1] = -4; /* xy */
  baseLGP2_2D[5].c[2][0] =  0; /* x² */
  baseLGP2_2D[5].c[0][2] = -4; /* y² */

//   beta_x[0] =   4*x+4*y-3;   /*  4x+4y-3 = -3*beta_0 +   beta_1 +  beta_2*/
//   beta_x[1] =  -8*x-4*y+4;   /* -8x-4y+4 =  4*beta_0 - 4*beta_1          */
//   beta_x[2] =   4*x    -1;   /*  4x   -1 =   -beta_0 + 3*beta_1 -  beta_2*/
//   beta_x[3] =      +4*y  ;   /*    +4y   =                       4*beta_2*/
//   beta_x[4] =           0;   /*        0 */
//   beta_x[5] =      -4*y  ;   /*    -4y   =                      -4*beta_2*/
//
//   beta_y[0] =   4*x+4*y-3;   /*  4x+4y-3 */
//   beta_y[1] =  -4*x      ;   /* -4x      */
//   beta_y[2] =           0;   /*        0 */
//   beta_y[3] =   4*x      ;   /*  4x      */
//   beta_y[4] =      +4*y-1;   /*    +4y-1 */
//   beta_y[5] =  -4*x-8*y+4;   /* -4x-8y+4 */

  baseLGP2_x_2D[0].c[0][0] = -3; /*  1 */
  baseLGP2_x_2D[0].c[1][0] =  4; /*  x */
  baseLGP2_x_2D[0].c[0][1] =  4; /*  y */

  baseLGP2_x_2D[1].c[0][0] =  4; /*  1 */
  baseLGP2_x_2D[1].c[1][0] = -8; /*  x */
  baseLGP2_x_2D[1].c[0][1] = -4; /*  y */

  baseLGP2_x_2D[2].c[0][0] = -1; /*  1 */
  baseLGP2_x_2D[2].c[1][0] =  4; /*  x */
  baseLGP2_x_2D[2].c[0][1] =  0; /*  y */

  baseLGP2_x_2D[3].c[0][0] =  0; /*  1 */
  baseLGP2_x_2D[3].c[1][0] =  0; /*  x */
  baseLGP2_x_2D[3].c[0][1] =  4; /*  y */

  baseLGP2_x_2D[4].c[0][0] =  0; /*  1 */
  baseLGP2_x_2D[4].c[1][0] =  0; /*  x */
  baseLGP2_x_2D[4].c[0][1] =  0; /*  y */

  baseLGP2_x_2D[5].c[0][0] =  0; /*  1 */
  baseLGP2_x_2D[5].c[1][0] =  0; /*  x */
  baseLGP2_x_2D[5].c[0][1] = -4; /*  y */

  baseLGP2_y_2D[0].c[0][0] = -3; /*  1 */
  baseLGP2_y_2D[0].c[1][0] =  4; /*  x */
  baseLGP2_y_2D[0].c[0][1] =  4; /*  y */

  baseLGP2_y_2D[1].c[0][0] =  0; /*  1 */
  baseLGP2_y_2D[1].c[1][0] = -4; /*  x */
  baseLGP2_y_2D[1].c[0][1] =  0; /*  y */

  baseLGP2_y_2D[2].c[0][0] =  0; /*  1 */
  baseLGP2_y_2D[2].c[1][0] =  0; /*  x */
  baseLGP2_y_2D[2].c[0][1] =  0; /*  y */

  baseLGP2_y_2D[3].c[0][0] =  0; /*  1 */
  baseLGP2_y_2D[3].c[1][0] =  4; /*  x */
  baseLGP2_y_2D[3].c[0][1] =  0; /*  y */

  baseLGP2_y_2D[4].c[0][0] = -1; /*  1 */
  baseLGP2_y_2D[4].c[1][0] =  0; /*  x */
  baseLGP2_y_2D[4].c[0][1] =  4; /*  y */

  baseLGP2_y_2D[5].c[0][0] =  4; /*  1 */
  baseLGP2_y_2D[5].c[1][0] = -4; /*  x */
  baseLGP2_y_2D[5].c[0][1] = -8; /*  y */

/*------------------------------------------------------------------------------
  define 3D Lagrange-P1 base functions, coefficients of cannonics */
  baseLGP1_3D[0].c[0][0][1] = -1; /* -1 * z */
  baseLGP1_3D[0].c[1][0][1] =  1; /*  x * z */
  baseLGP1_3D[0].c[0][1][1] =  1; /*  y * z */

  baseLGP1_3D[1].c[0][0][1] =  0; /*  0 * z */
  baseLGP1_3D[1].c[1][0][1] = -1; /* -x * z */
  baseLGP1_3D[1].c[0][1][1] =  0; /*  0 * z */

  baseLGP1_3D[2].c[0][0][1] =  0; /*  0 * z */
  baseLGP1_3D[2].c[1][0][1] =  0; /*  0 * z */
  baseLGP1_3D[2].c[0][1][1] = -1; /* -y * z */

  baseLGP1_3D[0].c[0][0][0] =  1; /*  1 */
  baseLGP1_3D[0].c[1][0][0] = -1; /* -x */
  baseLGP1_3D[0].c[0][1][0] = -1; /* -y */

  baseLGP1_3D[1].c[0][0][0] =  0; /*  0 */
  baseLGP1_3D[1].c[1][0][0] =  1; /*  x */
  baseLGP1_3D[1].c[0][1][0] =  0; /*  0 */

  baseLGP1_3D[2].c[0][0][0] =  0; /*  0 */
  baseLGP1_3D[2].c[1][0][0] =  0; /*  0 */
  baseLGP1_3D[2].c[0][1][0] =  1; /*  y */

  baseLGP1_3D[3].c[0][0][1] =  1; /*  1 * z */
  baseLGP1_3D[3].c[1][0][1] = -1; /* -x * z */
  baseLGP1_3D[3].c[0][1][1] = -1; /* -y * z */

  baseLGP1_3D[4].c[0][0][1] =  0; /*  0 * z */
  baseLGP1_3D[4].c[1][0][1] =  1; /*  x * z */
  baseLGP1_3D[4].c[0][1][1] =  0; /*  0 * z */

  baseLGP1_3D[5].c[0][0][1] =  0; /*  0 * z */
  baseLGP1_3D[5].c[1][0][1] =  0; /*  0 * z */
  baseLGP1_3D[5].c[0][1][1] =  1; /*  y * z */

/*------------------------------------------------------------------------------
  define 2D QP1 base functions, coefficients of cannonics */
//   beta[0] =  (1-x) (1-y) ;      /*  1 -x -y + xy  */
//   beta[1] =     x  (1-y) ;      /*     x    - xy  */
//   beta[2] =     x y      ;      /*            xy  */
//   beta[3] =   (1-x) y    ;      /*        y - xy  */

  baseQP1_2D[0].c[0][0] =  1; /*  1 */
  baseQP1_2D[0].c[1][0] = -1; /* -x */
  baseQP1_2D[0].c[0][1] = -1; /* -y */
  baseQP1_2D[0].c[1][1] =  1; /* xy */
  baseQP1_2D[0].c[2][0] =  0; /* x² */
  baseQP1_2D[0].c[0][2] =  0; /* y² */

  baseQP1_2D[1].c[0][0] =  0; /*  1 */
  baseQP1_2D[1].c[1][0] =  1; /*  x */
  baseQP1_2D[1].c[0][1] =  0; /*  y */
  baseQP1_2D[1].c[1][1] = -1; /* xy */
  baseQP1_2D[1].c[2][0] =  0; /* x² */
  baseQP1_2D[1].c[0][2] =  0; /* y² */

  baseQP1_2D[2].c[0][0] =  0; /*  1 */
  baseQP1_2D[2].c[1][0] =  0; /*  x */
  baseQP1_2D[2].c[0][1] =  0; /*  y */
  baseQP1_2D[2].c[1][1] =  1; /* xy */
  baseQP1_2D[2].c[2][0] =  0; /* x² */
  baseQP1_2D[2].c[0][2] =  0; /* y² */

  baseQP1_2D[3].c[0][0] =  0; /*  1 */
  baseQP1_2D[3].c[1][0] =  0; /*  x */
  baseQP1_2D[3].c[0][1] =  1; /*  y */
  baseQP1_2D[3].c[1][1] = -1; /* xy */
  baseQP1_2D[3].c[2][0] =  0; /* x² */
  baseQP1_2D[3].c[0][2] =  0; /* y² */

/*------------------------------------------------------------------------------
  allocate buffer polynomials for product operation*/
  status = polyinit(canonic2D, 2, &p1);
  status = polyinit(canonic2D, 3, &p2);
  status = polyinit(canonic2D, 4, &p3);
  status = polyinit(canonic2D, 5, &p4);
  status = polyinit(canonic2D, 6, &p5);

/*------------------------------------------------------------------------------
  allocate buffer polynomials for product operation*/
  status = polyinit3D( 2, &p3D1);
  status = polyinit3D( 3, &p3D2);
  status = polyinit3D( 4, &p3D3);
  status = polyinit3D( 5, &p3D4);
  status = polyinit3D( 6, &p3D5);

  out=fopen("integration.out","w");
  if(out==0) {
      TRAP_ERR_EXIT(-1,"Probleme a l ouverture du fichier : integration.out \n");
    }
/**----------------------------------------------------------------------------
  compute 1D LGP1 only related integrales*/
  for(i = 0; i < 2; i++) {
    for(j = 0; j < 2; j++) {
      status = polymul(canonic2D, baseLGP1_1D[i], baseLGP1_1D[j], p1);
      fe_LGP1xLGP1_1D[i][j] = polyint1D(canonic2D, p1);
      for(k = 0; k < 2; k++) {
        status = polymul(canonic2D, baseLGP1_1D[k], p1, p2);
        fe_LGP1xLGP1xLGP1_1D[i][j][k] = polyint1D(canonic2D, p2);
        for(l = 0; l < 2; l++) {
          status = polymul(canonic2D, baseLGP1_1D[l], p2, p3);
          fe_LGP1xLGP1xLGP1xLGP1_1D[i][j][k][l] = polyint1D(canonic2D, p3);
          for(m = 0; m < 2; m++) {
            status = polymul(canonic2D, baseLGP1_1D[m], p3, p4);
            fe_LGP1xLGP1xLGP1xLGP1xLGP1_1D[i][j][k][l][m] = polyint1D(canonic2D, p4);
            }
          }
        }
      for(k = 0; k < 1; k++) {
        status = polymul(canonic2D, baseNCP1_1D[k], p1, p2);
        fe_LGP1xLGP1xNCP1_1D[i][j][k] = polyint1D(canonic2D, p2);
        }
      }
    }
  for(i = 0; i < 2; i++) {
    for(j = 0; j < 2; j++) {
//      fe_LGP1xLGP1_1D[i][j] = polyint1D(canonic2D, p1);
      fprintf(out,"1D LGP1xLGP1: %d %d %f\n", i, j, fe_LGP1xLGP1_1D[i][j]);
      }
    }
  for(i = 0; i < 2; i++) {
    for(j = 0; j < 2; j++) {
      for(k = 0; k < 2; k++) {
        fprintf(out,"1D LGP1xLGP1xLGP1: %d %d %d %f\n", i, j, k, fe_LGP1xLGP1xLGP1_1D[i][j][k]);
        }
      }
    }
  for(i = 0; i < 2; i++) {
    for(j = 0; j < 2; j++) {
      for(k = 0; k < 1; k++) {
        fprintf(out,"1D LGP1xLGP1xNCP1: %d %d %d %f\n", i, j, k, fe_LGP1xLGP1xNCP1_1D[i][j][k]);
        }
      }
    }

/**----------------------------------------------------------------------
  compute 2D LGP1 only related integrales*/
  for(i = 0; i < 3; i++) {
    fe_LGP1_2D[i] = polyint2D(canonic2D, baseLGP1_2D[i]);
    for(j = 0; j < 3; j++) {
      status = polymul(canonic2D, baseLGP1_2D[i], baseLGP1_2D[j], p1);
      fe_LGP1xLGP1_2D[i][j] = polyint2D(canonic2D, p1);
      for(k = 0; k < 3; k++) {
        status = polymul(canonic2D, baseLGP1_2D[k], p1, p2);
        fe_LGP1xLGP1xLGP1_2D[i][j][k] = polyint2D(canonic2D, p2);
        for(l = 0; l < 3; l++) {
          status = polymul(canonic2D, baseNCP1_2D[l], p2, p3);
          fe_LGP1xLGP1xLGP1xNCP1_2D[i][j][k][l] = polyint2D(canonic2D, p3);
          }
        for(l = 0; l < 3; l++) {
          status = polymul(canonic2D, baseLGP1_2D[l], p2, p3);
          fe_LGP1xLGP1xLGP1xLGP1_2D[i][j][k][l] = polyint2D(canonic2D, p3);
          }
        }
      for(k = 0; k < 3; k++) {
        status = polymul(canonic2D, baseNCP1_2D[k], p1, p2);
        fe_LGP1xLGP1xNCP1_2D[i][j][k] = polyint2D(canonic2D, p2);
        for(l = 0; l < 3; l++) {
          status = polymul(canonic2D, baseNCP1_2D[l], p2, p3);
          fe_LGP1xLGP1xNCP1xNCP1_2D[i][j][k][l] = polyint2D(canonic2D, p3);
          }
        }
      for(k = 0; k < 3; k++) {
        status = polymul(canonic2D, p1, baseNCP1_2D[k], p2);
/**----------------------------------------------------------------------
        test*/
        test_fe_LGP1xLGP1xNCP1_2D[i][j][k] = polyint2D(canonic2D, p2);
        }
      }
    }
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      fprintf(out,"2D LGP1xLGP1: %d %d %f\n", i, j,fe_LGP1xLGP1_2D[i][j]);
      }
    }
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      for(k = 0; k < 3; k++) {
        fprintf(out,"2D LGP1xLGP1xLGP1: %d %d %d %lf\n", i, j, k, fe_LGP1xLGP1xLGP1_2D[i][j][k]);
        }
      }
    }
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      for(k = 0; k < 3; k++) {
        fprintf(out,"2D LGP1xLGP1xNCP1: %d %d %d %lf\n", i, j, k, fe_LGP1xLGP1xNCP1_2D[i][j][k]);
        }
      }
    }

/**----------------------------------------------------------------------
  compute 2D LGP2 only related integrales*/
  for(i = 0; i < 6; i++) {
    fe_LGP2_2D[i] = polyint2D(canonic2D, baseLGP2_2D[i]);
    for(j = 0; j < 6; j++) {
      status = polymul(canonic2D, baseLGP2_2D[i], baseLGP2_2D[j], p3);
      fe_LGP2xLGP2_2D[i][j] = polyint2D(canonic2D, p3);
      for(k = 0; k < 6; k++) {
        status = polymul(canonic2D, baseLGP2_2D[k], p3, p5);
        fe_LGP2xLGP2xLGP2_2D[i][j][k] = polyint2D(canonic2D, p5);
        }
      for(k = 0; k < 3; k++) {
        status = polymul(canonic2D, baseLGP1_2D[k], p3, p4);
        fe_LGP2xLGP2xLGP1_2D[i][j][k] = polyint2D(canonic2D, p4);
        }
      for(k = 0; k < 3; k++) {
        status = polymul(canonic2D, baseNCP1_2D[k], p3, p4);
        fe_LGP2xLGP2xNCP1_2D[i][j][k] = polyint2D(canonic2D, p4);
        }
      }
    }

/**----------------------------------------------------------------------
  compute 2D LGP2 x NCP1 related integrales*/
  for(i = 0; i < 6; i++) {
//    fe_LGP2_2D[i] = polyint2D(canonic2D, baseLGP2_2D[i]);
    for(j = 0; j < 3; j++) {
      status = polymul(canonic2D, baseLGP2_2D[i], baseNCP1_2D[j], p2);
      fe_LGP2xNCP1_2D[i][j] = polyint2D(canonic2D, p2);
      for(k = 0; k < 3; k++) {
        status = polymul(canonic2D, baseNCP1_2D[k], p2, p3);
        fe_LGP2xNCP1xNCP1_2D[i][j][k] = polyint2D(canonic2D, p3);
        }
      }
    for(j = 0; j < 3; j++) {
      status = polymul(canonic2D, baseLGP2_2D[i], baseLGP1_2D[j], p2);
      fe_LGP2xLGP1_2D[i][j] = polyint2D(canonic2D, p2);
      for(k = 0; k < 3; k++) {
        status = polymul(canonic2D, baseLGP1_2D[k], p2, p3);
        fe_LGP2xLGP1xLGP1_2D[i][j][k] = polyint2D(canonic2D, p3);
        }
      }
    }

/**----------------------------------------------------------------------
  compute 2D LGP1xLGP2 derivative related integrales*/
  for(i = 0; i < 6; i++) {
    for(j = 0; j < 3; j++) {
      status = polymul(canonic2D, baseLGP2_x_2D[i], baseLGP1_2D[j], p1);
      fe_dxP2xLGP1_2D[i][j] = polyint2D(canonic2D, p1);
      for(k = 0; k < 3; k++) {
        status = polymul(canonic2D, p1, baseLGP1_2D[k], p2);
        fe_dxP2xLGP1xLGP1_2D[i][j][k] = polyint2D(canonic2D, p2);
        }
      }
    }

  for(i = 0; i < 6; i++) {
    for(j = 0; j < 3; j++) {
      status = polymul(canonic2D, baseLGP2_y_2D[i], baseLGP1_2D[j], p1);
      fe_dyP2xLGP1_2D[i][j] = polyint2D(canonic2D, p1);
      for(k = 0; k < 3; k++) {
        status = polymul(canonic2D, p1, baseLGP1_2D[k], p2);
        fe_dyP2xLGP1xLGP1_2D[i][j][k] = polyint2D(canonic2D, p2);
        }
      }
    }

  for(i = 0; i < 6; i++) {
    for(j = 0; j < 6; j++) {
      status = polymul(canonic2D, baseLGP2_x_2D[i], baseLGP2_2D[j], p2);
      for(k = 0; k < 3; k++) {
        status = polymul(canonic2D, p2, baseLGP1_2D[k], p3);
        fe_dxP2xLGP2xLGP1_2D[i][j][k] = polyint2D(canonic2D, p3);
        }
      for(k = 0; k < 6; k++) {
        status = polymul(canonic2D, p2, baseLGP2_2D[k], p4);
        fe_dxP2xLGP2xLGP2_2D[i][j][k] = polyint2D(canonic2D, p4);
        }
      }
    }

  for(i = 0; i < 6; i++) {
    for(j = 0; j < 6; j++) {
      status = polymul(canonic2D, baseLGP2_y_2D[i], baseLGP2_2D[j], p2);
      for(k = 0; k < 3; k++) {
        status = polymul(canonic2D, p2, baseLGP1_2D[k], p3);
        fe_dyP2xLGP2xLGP1_2D[i][j][k] = polyint2D(canonic2D, p3);
        }
      for(k = 0; k < 6; k++) {
        status = polymul(canonic2D, p2, baseLGP2_2D[k], p4);
        fe_dyP2xLGP2xLGP2_2D[i][j][k] = polyint2D(canonic2D, p4);
        }
      }
    }

/**----------------------------------------------------------------------
  compute 2D NCP1xLGP2 derivative related integrales*/
  for(i = 0; i < 6; i++) {
    for(j = 0; j < 3; j++) {
      status = polymul(canonic2D, baseLGP2_x_2D[i], baseNCP1_2D[j], p1);
      fe_dxP2xNCP1_2D[i][j] = polyint2D(canonic2D, p1);
      }
    }

  for(i = 0; i < 6; i++) {
    for(j = 0; j < 3; j++) {
      status = polymul(canonic2D, baseLGP2_y_2D[i], baseNCP1_2D[j], p1);
      fe_dyP2xNCP1_2D[i][j] = polyint2D(canonic2D, p1);
      }
    }

  for(i = 0; i < 6; i++) {
    for(j = 0; j < 6; j++) {
      status = polymul(canonic2D, baseLGP2_x_2D[i], baseLGP2_2D[j], p2);
      for(k = 0; k < 3; k++) {
        status = polymul(canonic2D, p2, baseNCP1_2D[k], p3);
        fe_dxP2xLGP2xNCP1_2D[i][j][k] = polyint2D(canonic2D, p3);
        }
      }
    }
  for(i = 0; i < 6; i++) {
    for(j = 0; j < 6; j++) {
      status = polymul(canonic2D, baseLGP2_y_2D[i], baseLGP2_2D[j], p2);
      for(k = 0; k < 3; k++) {
        status = polymul(canonic2D, p2, baseNCP1_2D[k], p3);
        fe_dyP2xLGP2xNCP1_2D[i][j][k] = polyint2D(canonic2D, p3);
        }
      }
    }

/**----------------------------------------------------------------------
  compute 2D NCP1 only related integrales*/
  for(i = 0; i < 3; i++) {
    fe_NCP1_2D[i] = polyint2D(canonic2D, baseNCP1_2D[i]);
    for(j = 0; j < 3; j++) {
      status = polymul(canonic2D, baseNCP1_2D[i], baseNCP1_2D[j], p1);
      fe_NCP1xNCP1_2D[i][j] = polyint2D(canonic2D, p1);
      for(k = 0; k < 3; k++) {
        status = polymul(canonic2D, baseNCP1_2D[k], p1, p2);
        fe_NCP1xNCP1xNCP1_2D[i][j][k] = polyint2D(canonic2D, p2);
        }
      }
    }
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      fprintf(out,"2D NCP1xNCP1: %d %d %f\n", i, j, fe_NCP1xNCP1_2D[i][j]);
      }
    }
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      for(k = 0; k < 3; k++) {
        fprintf(out,"2D NCP1xNCP1xNCP1: %d %d %d %lf\n", i, j, k, fe_NCP1xNCP1xNCP1_2D[i][j][k]);
        }
      }
    }
/*------------------------------------------------------------------------------
  compute 2D NCP1 and LGP1 cross integrales*/
//   for(i = 0; i < 3; i++) {
//     for(j = 0; j < 3; j++) {
//       status = polymul(canonic2D, baseLGP1_2D[i], baseNCP1_2D[j], p1);
//       sum2x2[i][j] = polyint2D(canonic2D, p1);
//       fe_LGP1xNCP1[i][j] = polyint2D(canonic2D, p1);
//       fprintf(out,"2D LGP1xNCP1: %d %d %f\n", i, j, sum2x2[i][j]);
//       }
//     }
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      status = polymul(canonic2D, baseLGP1_2D[i], baseNCP1_2D[j], p1);
      fe_LGP1xNCP1_2D[i][j] = polyint2D(canonic2D, p1);
      for(k = 0; k < 3; k++) {
        status = polymul(canonic2D, baseLGP1_2D[k], p1, p2);
        fe_LGP1xLGP1xNCP1_2D[i][k][j] = polyint2D(canonic2D, p2);
        }
      }
    }
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      fprintf(out,"2D LGP1xNCP1: %d %d %f\n", i, j, fe_LGP1xNCP1_2D[i][j]);
      }
    }
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      for(k = 0; k < 3; k++) {
        fprintf(out,"2D LGP1xLGP1xNCP1: %d %d %d %lf\n", i, j, k, fe_LGP1xLGP1xNCP1_2D[i][j][k]);
        }
      }
    }
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      status = polymul(canonic2D, baseLGP1_2D[i], baseNCP1_2D[j], p1);
      for(k = 0; k < 3; k++) {
        status = polymul(canonic2D,p1, baseNCP1_2D[k],p2);
        fe_LGP1xNCP1xNCP1_2D[i][j][k] = polyint2D(canonic2D, p2);
        }
      }
    }
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      for(k = 0; k < 3; k++) {
        fprintf(out,"2D LGP1xNCP1xNCP1: %d %d %d %lf\n", i, j, k, fe_LGP1xNCP1xNCP1_2D[i][j][k]);
        }
      }
    }
/*------------------------------------------------------------------------------
  compute 3D LGP1  integrales*/
  for(i = 0; i < 6; i++) {
    fe_LGP1_3D[i] = polyint3D(canonic3D, baseLGP1_3D[i]);
    for(j = 0; j < 6; j++) {
      status = polymul3D(canonic3D, baseLGP1_3D[i], baseLGP1_3D[j], p3D3);
      fe_LGP1xLGP1_3D[i][j] = polyint3D(canonic3D, p3D3);
//      fprintf(out,"3D LGP1xLGP1: %d %d %f\n", i, j, fe_LGP1xLGP1_3D[i][j]);
      for(k = 0; k < 6; k++) {
        status = polymul3D(canonic3D,p3D3, baseLGP1_3D[k],p3D5);
        fe_LGP1xLGP1xLGP1_3D[i][j][k] = polyint3D(canonic3D, p3D5);
//        fprintf(out,"3D LGP1xLGP1xLGP1: %d %d %d %f\n", i, j, k, fe_LGP1xLGP1xLGP1_3D[i][j][k]);
        }
      }
    }

  for(i = 0; i < 6; i++) {
    fprintf(out,"3D LGP1: %d %f\n", i, fe_LGP1_3D[i]);
    }
  for(i = 0; i < 6; i++) {
    for(j = 0; j < 6; j++) {
      fprintf(out,"3D LGP1xLGP1: %d %d %f\n", i, j, fe_LGP1xLGP1_3D[i][j]);
      }
    }
  for(i = 0; i < 6; i++) {
    for(j = 0; j < 6; j++) {
      for(k = 0; k < 6; k++) {
        fprintf(out,"3D LGP1xLGP1xLGP1: %d %d %d %f\n", i, j, k, fe_LGP1xLGP1xLGP1_3D[i][j][k]);
        }
      }
    }

/**----------------------------------------------------------------------
  compute 2D QP1 related integrales*/
  for(i = 0; i < 4; i++) {
    fe_QP1_2D[i] = polyint_Q2D(canonic2D, baseQP1_2D[i]);
    for(j = 0; j < 4; j++) {
      status = polymul(canonic2D, baseQP1_2D[i], baseQP1_2D[j], p3);
      fe_QP1xQP1_2D[i][j] = polyint_Q2D(canonic2D, p3);
      for(k = 0; k < 4; k++) {
        status = polymul(canonic2D, baseQP1_2D[k], p3, p5);
        fe_QP1xQP1xQP1_2D[i][j][k] = polyint_Q2D(canonic2D, p5);
        }
      }
    }

  status=fe_chk_integrales(out);
  fclose(out);

/*------------------------------------------------------------------------------
  free 1D Lagrange-P1 base functions*/
  for(k = 0; k < 2; k++) {
    status = polyfree(&(baseLGP1_1D[k]));
    }

/*------------------------------------------------------------------------------
  free 1D non conforming-P1 base functions*/
  for(k = 0; k < 1; k++) {
    status = polyfree(&(baseNCP1_1D[k]));
    }

/*------------------------------------------------------------------------------
  free 2D Lagrange-P1 base functions*/
  for(k = 0; k < 3; k++) {
    status = polyfree(&(baseLGP1_2D[k]));
    status = polyfree(&(baseNCP1_2D[k]));
    }

/*------------------------------------------------------------------------------
  free 1D Lagrange-P2 base functions*/
  for(k = 0; k < 3; k++) {
    status = polyfree(&(baseLGP2_1D[k]));
    }

/*------------------------------------------------------------------------------
  free 2D Lagrange-P2 base functions*/
  for(k = 0; k < 6; k++) {
    status = polyfree(&(baseLGP2_2D[k]));
    status = polyfree(&(baseLGP2_x_2D[k]));
    status = polyfree(&(baseLGP2_y_2D[k]));
    }

  status = polyfree(&p1);
  status = polyfree(&p2);
  status = polyfree(&p3);
  status = polyfree(&p4);
  status = polyfree(&p5);

  status = polyfree(&p3D1);
  status = polyfree(&p3D2);
  status = polyfree(&p3D3);
  status = polyfree(&p3D4);
  status = polyfree(&p3D5);
  for(k = 0; k < 6; k++) {
    status = polyfree(&(baseLGP1_3D[k]));
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T1,typename T2,typename T3> T1 fe_integraleLGP1xLGP1xLGP1_2D_template(T1 *p, T2 *q, T3 *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
  see integration tabulations*/
#ifdef NOMINAL
  int i, j, k;
  T1 sum;
  sum = 0;

  for(k = 0; k < 3; k++) {
    for(j = 0; j < 3; j++) {
      for(i = 0; i < 3; i++) {
        sum += p[i] * q[j] * r[k] * fe_LGP1xLGP1xLGP1_2D[i][j][k];
        }
      }
    }
  return (sum);
#endif

#ifdef OPTIMAL
  T1 tmp1,tmp2,tmp3,sum2;
  tmp1  = p[0]*q[0]*r[0]+p[1]*q[1]*r[1]+p[2]*q[2]*r[2];
  tmp2  = p[0]*q[1]*r[1]+p[0]*q[2]*r[2]+p[1]*q[0]*r[0]+p[1]*q[2]*r[2]+p[2]*q[0]*r[0]+p[2]*q[1]*r[1];
  tmp2 += p[1]*q[0]*r[1]+p[2]*q[0]*r[2]+p[0]*q[1]*r[0]+p[2]*q[1]*r[2]+p[0]*q[2]*r[0]+p[1]*q[2]*r[1];
  tmp2 += p[1]*q[1]*r[0]+p[2]*q[2]*r[0]+p[0]*q[0]*r[1]+p[2]*q[2]*r[1]+p[0]*q[0]*r[2]+p[1]*q[1]*r[2];
  tmp3  = p[0]*q[1]*r[2]+p[0]*q[2]*r[1]+p[1]*q[0]*r[2]+p[1]*q[2]*r[0]+p[2]*q[0]*r[1]+p[2]*q[1]*r[0];

  sum2=tmp1*fe_LGP1xLGP1xLGP1_2D[0][0][0]+tmp2*fe_LGP1xLGP1xLGP1_2D[0][1][0]+tmp3*fe_LGP1xLGP1xLGP1_2D[0][1][2];

  return (sum2);
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_integraleLGP1xLGP1xLGP1_2D(double *p, double *q, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=fe_integraleLGP1xLGP1xLGP1_2D_template(p,q,r);
  
  return result;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> fe_integraleLGP1xLGP1xLGP1_2D(double *p, complex<double> *q, complex<double> *r)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n;
  complex<double> sum,tmp1,tmp2,tmp3;
/*----------------------------------------------------------------------
  see integration tabulations*/
  tmp1  = p[0]*q[0]*r[0]+p[1]*q[1]*r[1]+p[2]*q[2]*r[2];
  tmp2  = p[0]*q[1]*r[1]+p[0]*q[2]*r[2]+p[1]*q[0]*r[0]+p[1]*q[2]*r[2]+p[2]*q[0]*r[0]+p[2]*q[1]*r[1];
  tmp2 += p[1]*q[0]*r[1]+p[2]*q[0]*r[2]+p[0]*q[1]*r[0]+p[2]*q[1]*r[2]+p[0]*q[2]*r[0]+p[1]*q[2]*r[1];
  tmp2 += p[1]*q[1]*r[0]+p[2]*q[2]*r[0]+p[0]*q[0]*r[1]+p[2]*q[2]*r[1]+p[0]*q[0]*r[2]+p[1]*q[1]*r[2];
  tmp3  = p[0]*q[1]*r[2]+p[0]*q[2]*r[1]+p[1]*q[0]*r[2]+p[1]*q[2]*r[0]+p[2]*q[0]*r[1]+p[2]*q[1]*r[0];

  sum=tmp1*fe_LGP1xLGP1xLGP1_2D[0][0][0]+tmp2*fe_LGP1xLGP1xLGP1_2D[0][1][0]+tmp3*fe_LGP1xLGP1xLGP1_2D[0][1][2];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP1xLGP1_2D(complex<double> *p, double *q, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xLGP1xLGP1_2D_template(p,q,r);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP1xLGP1_2D(complex<double> *p, complex<double> *q, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xLGP1xLGP1_2D_template(p,q,r);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP1xLGP1_2D(complex<double> *p, complex<double> *q, complex<double> *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xLGP1xLGP1_2D_template(p,q,r);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP1xLGP0_2D(double *p, complex<double> *q, complex<double> *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  complex<double> sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(j = 0; j < 3; j++)
    for(i = 0; i < 3; i++)
      sum += p[i] * q[j] * fe_LGP1xLGP1_2D[i][j];

  sum*=r[0];
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T1,typename T2>
            T1 fe_integraleLGP1xLGP1xLGP1_2D_template(T1 *p, T2 *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
  see integration tabulations*/
#ifdef NOMINAL
  int i, j;
  T1 sum;
  sum = 0;

  for(j = 0; j < 3; j++)
    for(i = 0; i < 3; i++)
      sum += p[i] * q[j] * fe_LGP1xLGP1xLGP1_2D[i][j][k];

  return (sum);
#endif

#ifdef OPTIMAL
  T1 sum2;
  sum2 = 0;

  switch(k) {
    case 0:

      sum2+=fe_LGP1xLGP1xLGP1_2D[0][0][0]*(p[0]*q[0]); // 0.050000
      sum2+=fe_LGP1xLGP1xLGP1_2D[0][1][0]*(p[0]*q[1]+p[0]*q[2]+p[1]*q[0]+p[1]*q[1]+p[2]*q[0]+p[2]*q[2]); // 0.016667
      sum2+=fe_LGP1xLGP1xLGP1_2D[1][2][0]*(p[1]*q[2]+p[2]*q[1]); // 0.008333
      break;

    case 1:

      sum2+=fe_LGP1xLGP1xLGP1_2D[0][0][1]*(p[0]*q[0]+p[0]*q[1]+p[1]*q[0]+p[1]*q[2]+p[2]*q[1]+p[2]*q[2]); // 0.016667
      sum2+=fe_LGP1xLGP1xLGP1_2D[0][2][1]*(p[0]*q[2]+p[2]*q[0]); // 0.008333
      sum2+=fe_LGP1xLGP1xLGP1_2D[1][1][1]*(p[1]*q[1]); // 0.050000
      break;

    case 2:

      sum2+=fe_LGP1xLGP1xLGP1_2D[0][0][2]*(p[0]*q[0]+p[0]*q[2]+p[1]*q[1]+p[1]*q[2]+p[2]*q[0]+p[2]*q[1]); // 0.016667
      sum2+=fe_LGP1xLGP1xLGP1_2D[0][1][2]*(p[0]*q[1]+p[1]*q[0]); // 0.008333
      sum2+=fe_LGP1xLGP1xLGP1_2D[2][2][2]*(p[2]*q[2]); // 0.050000
      break;
    }
  return (sum2);
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1xLGP1_2D(double *p, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=fe_integraleLGP1xLGP1xLGP1_2D_template(p,q,k);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP1xLGP1_2D(double *p, complex<double> *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xLGP1xLGP1_2D_template(q,p,k);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP1xLGP1_2D(complex<double> *p, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xLGP1xLGP1_2D_template(p,q,k);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP1xLGP1_2D(complex<double> *p, complex<double> *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xLGP1xLGP1_2D_template(p,q,k);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_integraleLGP1xLGP1xLGP1xLGP1_2D(double *p, double *q, double *r, int l)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum2;
/*------------------------------------------------------------------------------
  see integration tabulations*/
#ifdef NOMINAL
  int i, j, k;
  double sum,tmp;
  sum = 0;

  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      tmp=p[i] * q[j];
      for(k = 0; k < 3; k++) {
//        sum += p[i] * q[j] * r[k] * fe_LGP1xLGP1xLGP1xLGP1_2D[i][j][k][l];
        sum += tmp * r[k] * fe_LGP1xLGP1xLGP1xLGP1_2D[i][j][k][l];
        }
      }
    }

  return (sum);
#endif

#ifdef OPTIMAL
  sum2 = 0;

  switch (l) {
    case 0:

      sum2+=fe_LGP1xLGP1xLGP1xLGP1_2D[0][0][0][0]*(p[0]*q[0]*r[0]); // 0.033333
      sum2+=fe_LGP1xLGP1xLGP1xLGP1_2D[0][0][1][0]*(p[0]*q[0]*r[1]+p[0]*q[0]*r[2]+p[0]*q[1]*r[0]+p[0]*q[2]*r[0]
                                                  +p[1]*q[0]*r[0]+p[1]*q[1]*r[1]+p[2]*q[0]*r[0]
                                                  +p[2]*q[2]*r[2]); // 0.008333
      sum2+=fe_LGP1xLGP1xLGP1xLGP1_2D[0][1][1][0]*(p[0]*q[1]*r[1]+p[0]*q[2]*r[2]+p[1]*q[0]*r[1]+p[1]*q[1]*r[0]
                                                  +p[2]*q[0]*r[2]+p[2]*q[2]*r[0]); // 0.005556
      sum2+=fe_LGP1xLGP1xLGP1xLGP1_2D[0][1][2][0]*(p[0]*q[1]*r[2]+p[0]*q[2]*r[1]+p[1]*q[0]*r[2]+p[1]*q[1]*r[2]
                                                  +p[1]*q[2]*r[0]+p[1]*q[2]*r[1]+p[1]*q[2]*r[2]
                                                  +p[2]*q[0]*r[1]+p[2]*q[1]*r[0]+p[2]*q[1]*r[1]
                                                  +p[2]*q[1]*r[2]+p[2]*q[2]*r[1]); // 0.002778
      break;

    case 1:

      sum2+=fe_LGP1xLGP1xLGP1xLGP1_2D[0][0][0][1]*(p[0]*q[0]*r[0]+p[0]*q[1]*r[1]+p[1]*q[0]*r[1]+p[1]*q[1]*r[0]
                                                  +p[1]*q[1]*r[2]+p[1]*q[2]*r[1]+p[2]*q[1]*r[1]
                                                  +p[2]*q[2]*r[2]); // 0.008333
      sum2+=fe_LGP1xLGP1xLGP1xLGP1_2D[0][0][1][1]*(p[0]*q[0]*r[1]+p[0]*q[1]*r[0]+p[1]*q[0]*r[0]+p[1]*q[2]*r[2]
                                                  +p[2]*q[1]*r[2]+p[2]*q[2]*r[1]); // 0.005556
      sum2+=fe_LGP1xLGP1xLGP1xLGP1_2D[0][0][2][1]*(p[0]*q[0]*r[2]+p[0]*q[1]*r[2]+p[0]*q[2]*r[0]+p[0]*q[2]*r[1]
                                                  +p[0]*q[2]*r[2]+p[1]*q[0]*r[2]+p[1]*q[2]*r[0]
                                                  +p[2]*q[0]*r[0]+p[2]*q[0]*r[1]+p[2]*q[0]*r[2]
                                                  +p[2]*q[1]*r[0]+p[2]*q[2]*r[0]); // 0.002778
      sum2+=fe_LGP1xLGP1xLGP1xLGP1_2D[1][1][1][1]*(p[1]*q[1]*r[1]); // 0.033333
      break;

    case 2:

      sum2+=fe_LGP1xLGP1xLGP1xLGP1_2D[0][0][0][2]*(p[0]*q[0]*r[0]+p[0]*q[2]*r[2]+p[1]*q[1]*r[1]+p[1]*q[2]*r[2]
                                                  +p[2]*q[0]*r[2]+p[2]*q[1]*r[2]+p[2]*q[2]*r[0]
                                                  +p[2]*q[2]*r[1]); // 0.008333
      sum2+=fe_LGP1xLGP1xLGP1xLGP1_2D[0][0][1][2]*(p[0]*q[0]*r[1]+p[0]*q[1]*r[0]+p[0]*q[1]*r[1]+p[0]*q[1]*r[2]
                                                  +p[0]*q[2]*r[1]+p[1]*q[0]*r[0]+p[1]*q[0]*r[1]
                                                  +p[1]*q[0]*r[2]+p[1]*q[1]*r[0]+p[1]*q[2]*r[0]
                                                  +p[2]*q[0]*r[1]+p[2]*q[1]*r[0]); // 0.002778
      sum2+=fe_LGP1xLGP1xLGP1xLGP1_2D[0][0][2][2]*(p[0]*q[0]*r[2]+p[0]*q[2]*r[0]+p[1]*q[1]*r[2]+p[1]*q[2]*r[1]
                                                  +p[2]*q[0]*r[0]+p[2]*q[1]*r[1]); // 0.005556
      sum2+=fe_LGP1xLGP1xLGP1xLGP1_2D[2][2][2][2]*(p[2]*q[2]*r[2]); // 0.033333
      break;
    }
  return (sum2);
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T1,typename T2>  T1 fe_integraleLGP1xLGP0_2D_template(T1 *p, T2 *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  T1 sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(i = 0; i < 3; i++)
    sum += p[i] * fe_LGP1_2D[i];

  sum*=q[0];
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP0_2D(double *p, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=fe_integraleLGP1xLGP0_2D_template(p,q);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP0_2D(complex<double> *p, complex<double> *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xLGP0_2D_template(p,q);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T1,typename T2>  T1 fe_integraleLGP1xNCP1_2D_template(T1 *p, T2 *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T1 sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/

#ifdef NOMINAL
  int i, j;
  sum = 0;
  for(j = 0; j < 3; j++)
    for(i = 0; i < 3; i++)
      if(i != j)
        sum += p[i] * q[j] / (double) 12.;

#endif

#ifdef OPTIMAL
   sum = (p[0]*(q[1]+q[2]) + p[1]*(q[0]+q[2]) + p[2]*(q[0]+q[1])) / (double) 12.;
#endif

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xNCP1_2D(double *p, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=fe_integraleLGP1xNCP1_2D_template(p,q);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xNCP1_2D(complex<double> *p, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xNCP1_2D_template(p,q);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xNCP1_2D(complex<double> *p, complex<double> *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xNCP1_2D_template(p,q);
  
  return result;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T1,typename T2,typename T3>  void fe_integraleLGP1xNCP1_2D_template(T1 *p, T2 *q, T3 & sum)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
  see integration tabulations*/

#ifdef NOMINAL
  int i, j;
  sum = 0;
  for(j = 0; j < 3; j++)
    for(i = 0; i < 3; i++)
      if(i != j)
        sum += p[i] * q[j] / (double) 12.;

#endif

#ifdef OPTIMAL
   sum = (p[0]*(q[1]+q[2]) + p[1]*(q[0]+q[2]) + p[2]*(q[0]+q[1])) / (double) 12.;
#endif

}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xNCP1_2D(double *p, complex<double> *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  fe_integraleLGP1xNCP1_2D_template(p,q,result);
  
  return result;
}


#if 0
//On one hand, these functions give the same results ...
#define fe_integraleNCP1_2D fe_integraleLGP1_2D
//... but on the other, this might confuse people when expanding to products involving NCP1 discretisations.
#else
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> T fe_integraleNCP1_2D_template(T *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  T sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/

  sum = 0;

  for(i = 0; i < 3; i++)
     sum += p[i] / (T) 6.;
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float fe_integraleNCP1_2D(float *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float result;
  
  result=fe_integraleNCP1_2D_template(p);
  
  return result;
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_integraleNCP1_2D(double *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=fe_integraleNCP1_2D_template(p);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<float> fe_integraleNCP1_2D(complex<float> *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<float> result;
  
  result=fe_integraleNCP1_2D_template(p);
  
  return result;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> fe_integraleNCP1_2D(complex<double> *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleNCP1_2D_template(p);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_integraleNCP1_2D(int i)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum =(double) 1./(double) 6.;

  return (sum);
}
#endif

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T1> T1 fe_integraleLGP1xNCP1_2D_template(T1 *p, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T1 sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/

//   int i;
//   sum = 0;
//
//   for(i = 0; i < 3; i++)
//     if(i != j)
//       sum += p[i] / (double) 12.;

  switch(j) {
    case 0:
      sum=(p[1]+p[2])/(double) 12.0;
      break;
    case 1:
      sum=(p[0]+p[2])/(double) 12.0;
      break;
    case 2:
      sum=(p[0]+p[1])/(double) 12.0;
      break;
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xNCP1_2D(double *p, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=fe_integraleLGP1xNCP1_2D_template(p,j);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xNCP1_2D(complex<double> *p, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xNCP1_2D_template(p,j);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T1> T1 fe_integraleLGP1xNCP1_2D_template(int i, T1 *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j;
  T1 sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;
  for(j = 0; j < 3; j++)
    if(i != j)
        sum +=  q[j] / (double) 12.;

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xNCP1_2D(int i, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=fe_integraleLGP1xNCP1_2D_template(i,q);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xNCP1_2D(int i, complex<double> *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xNCP1_2D_template(i,q);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> T fe_integraleLGP1_2D_template(T *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  T sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(i = 0; i < 3; i++)
    sum += p[i]/(double) 6.;

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1_2D(double *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=fe_integraleLGP1_2D_template(p);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1_2D(complex<double> *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1_2D_template(p);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_integraleLGP1_2D(int i)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum =(double) 1./(double) 6.;

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T1,typename T2>
                  T1 fe_integraleLGP1xLGP1_2D_template(T1 *p, T2 *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  T1 sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0.;

  for(j = 0; j < 3; j++)
    for(i = 0; i < 3; i++)
      sum += p[i] * q[j] * fe_LGP1xLGP1_2D[i][j];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1_2D(double *p, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=fe_integraleLGP1xLGP1_2D_template(p,q);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP1_2D(complex<double> *p, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xLGP1_2D_template(p,q);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP1_2D(complex<double> *p,complex<double> *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xLGP1_2D_template(p,q);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double integraleP1xP1(double *p, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Notes

  22/04/2008:

    same as above function, but more optimal numerically (is it?)

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  sum = 0;
  for(j = 0; j < 3; j++)
    sum += p[j] * q[j] / (double) 12.;
  for(j = 0; j < 3; j++)
    for(i = 0; i < 3; i++)
      if(i != j)
        sum += p[i] * q[j] / (double) 24.;
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> T fe_integraleLGP1xLGP1_2D_template(T *p, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  T sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(i = 0; i < 3; i++)
    sum += p[i] * fe_LGP1xLGP1_2D[i][j];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1_2D(double *p, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=fe_integraleLGP1xLGP1_2D_template(p,j);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP1_2D(complex<double> *p, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xLGP1_2D_template(p,j);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1_2D(int i, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum =fe_LGP1xLGP1_2D[i][j];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP2_2D(double *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  double sum;
  sum = 0;

  for(i = 0; i < 6; i++)
    sum += p[i] * fe_LGP2_2D[i];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP2xNCP1_2D(double *p, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  double sum;
  sum = 0;

  for(i = 0; i < 6; i++)
    for(j = 0; j < 3; j++)
      sum += p[i] * q[j] * fe_LGP2xNCP1_2D[i][j];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP2xNCP1_2D(double *p, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  double sum;
  sum = 0;

  for(i = 0; i < 6; i++)
      sum += p[i] * fe_LGP2xNCP1_2D[i][j];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP2xNCP1_2D(int i, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j;
  double sum;
  sum = 0;

  for(j = 0; j < 3; j++)
    sum += q[j] * fe_LGP2xNCP1_2D[i][j];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP2xNCP1_2D(complex<double> *p, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  complex<double> sum;
  sum = 0;

  for(i = 0; i < 6; i++)
    sum += p[i] * fe_LGP2xNCP1_2D[i][j];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP2xNCP1_2D(int j,complex<double> *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  complex<double> sum;
  sum = 0;

  for(i = 0; i < 6; i++)
    sum += p[j] * fe_LGP2xNCP1_2D[i][j];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP2xLGP1_2D(double *p, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  double sum;
  sum = 0;

  for(i = 0; i < 6; i++)
    for(j = 0; j < 3; j++)
      sum += p[i] * q[j] * fe_LGP2xLGP1_2D[i][j];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP2xLGP1_2D(int i, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j;
  double sum;
  sum = 0;

  for(j = 0; j < 3; j++)
    sum += q[j] * fe_LGP2xLGP1_2D[i][j];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP2xLGP1_2D(double *p, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  double sum;
  sum = 0;

  for(i = 0; i < 6; i++)
    sum += p[i] * fe_LGP2xLGP1_2D[i][j];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP2xLGP1_2D(complex<double> *p, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  complex<double> sum;
  sum = 0;

  for(i = 0; i < 6; i++)
    sum += p[i] * fe_LGP2xLGP1_2D[i][j];

  return (sum);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   complex<double> fe_integraleLGP2xLGP1_2D(complex<double> *p, int j)
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int i;
//   complex<double> sum;
//   sum = 0;
//
//   for(i = 0; i < 6; i++)
//     sum += p[i] * fe_LGP2xLGP1_2D[i][j];
//
//   return (sum);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP2xLGP1_2D(int j,complex<double> *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  complex<double> sum;
  sum = 0;

  for(i = 0; i < 6; i++)
    sum += p[j] * fe_LGP2xLGP1_2D[i][j];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP2xLGP1xLGP1_2D(double *p, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(j = 0; j < 3; j++) {
    for(i = 0; i < 6; i++) {
      sum += p[i] * q[j] * fe_LGP2xLGP1xLGP1_2D[i][j][k];
      }
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP2xLGP1xLGP1_2D(double *p, double *q, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j,k;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(k = 0; k < 3; k++) {
    for(j = 0; j < 3; j++) {
      for(i = 0; i < 6; i++) {
        sum += p[i] * q[j] * r[k] * fe_LGP2xLGP1xLGP1_2D[i][j][k];
        }
      }
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP2xLGP1xLGP1_2D(double *p, complex<double> *q, complex<double> *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j,k;
  complex<double> sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(k = 0; k < 3; k++) {
    for(j = 0; j < 3; j++) {
      for(i = 0; i < 6; i++) {
        sum += p[i] * q[j] * r[k] * fe_LGP2xLGP1xLGP1_2D[i][j][k];
        }
      }
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double>  fe_integraleLGP2xLGP2xLGP1_2D(double *p, complex<double>  *q, complex<double>  *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j,k;
  complex<double>  sum;

  sum = 0;

  for(k = 0; k < 3; k++) {
    for(j = 0; j < 6; j++) {
      for(i = 0; i < 6; i++) {
        sum += p[i] * q[j] * r[k] * fe_LGP2xLGP2xLGP1_2D[i][j][k];
        }
      }
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double>  fe_integraleLGP2xLGP2xNCP1_2D(double *p, complex<double>  *q, complex<double>  *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j,k;
  complex<double>  sum;

  sum = 0;

  for(k = 0; k < 3; k++) {
    for(j = 0; j < 6; j++) {
      for(i = 0; i < 6; i++) {
        sum += p[i] * q[j] * r[k] * fe_LGP2xLGP2xNCP1_2D[i][j][k];
        }
      }
    }
  return (sum);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> T fe_integraleLGP2xLGP2_2D_template(T *p, T *q)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n;
  T sum,sum2;
  double tmp1,tmp2,tmp3;
/*----------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(j = 0; j < 6; j++) {
    for(i = 0; i < 6; i++) {
      sum += p[i] * q[j] * fe_LGP2xLGP2_2D[i][j];
      }
    }
  return (sum);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_integraleLGP2xLGP2_2D(double *p, double *q)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
  sum=fe_integraleLGP2xLGP2_2D_template(p, q);
  return (sum);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> fe_integraleLGP2xLGP2_2D(complex<double>  *p, complex<double>  *q)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> sum;
  sum=fe_integraleLGP2xLGP2_2D_template(p, q);
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T1, typename T2> T1 fe_integraleLGP2xLGP2xLGP2_2D_template(T2 *p, T1 *q, T1 *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  int i, j, k;
  T1 sum;
  sum = 0;

  for(k = 0; k < 6; k++) {
    for(j = 0; j < 6; j++) {
      for(i = 0; i < 6; i++) {
        sum += p[i] * q[j] * r[k] * fe_LGP2xLGP2xLGP2_2D[i][j][k];
        }
      }
    }
  return (sum);
// #endif

// #ifdef OPTIMAL
//   T1 tmp1,tmp2,tmp3,sum2;
//   tmp1  = p[0]*q[0]*r[0]+p[1]*q[1]*r[1]+p[2]*q[2]*r[2];
//   tmp2  = p[0]*q[1]*r[1]+p[0]*q[2]*r[2]+p[1]*q[0]*r[0]+p[1]*q[2]*r[2]+p[2]*q[0]*r[0]+p[2]*q[1]*r[1];
//   tmp2 += p[1]*q[0]*r[1]+p[2]*q[0]*r[2]+p[0]*q[1]*r[0]+p[2]*q[1]*r[2]+p[0]*q[2]*r[0]+p[1]*q[2]*r[1];
//   tmp2 += p[1]*q[1]*r[0]+p[2]*q[2]*r[0]+p[0]*q[0]*r[1]+p[2]*q[2]*r[1]+p[0]*q[0]*r[2]+p[1]*q[1]*r[2];
//   tmp3  = p[0]*q[1]*r[2]+p[0]*q[2]*r[1]+p[1]*q[0]*r[2]+p[1]*q[2]*r[0]+p[2]*q[0]*r[1]+p[2]*q[1]*r[0];
//
//   sum2=tmp1*fe_LGP1xLGP1xLGP1_2D[0][0][0]+tmp2*fe_LGP1xLGP1xLGP1_2D[0][1][0]+tmp3*fe_LGP1xLGP1xLGP1_2D[0][1][2];
//
//   return (sum2);
// #endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP2xLGP2xLGP2_2D(double *p, double *q, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum=fe_integraleLGP2xLGP2xLGP2_2D_template(p, q, r);
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP2xLGP2xLGP2_2D(complex<double> *p, complex<double> *q, complex<double> *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> sum=fe_integraleLGP2xLGP2xLGP2_2D_template(p, q, r);
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP2xLGP2xLGP2_2D(double *p, complex<double> *q, complex<double> *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> sum=fe_integraleLGP2xLGP2xLGP2_2D_template(p, q, r);
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integrale_dxP2xLGP1_2D(double *p, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(i = 0; i < 6; i++) {
    sum += p[i] * fe_dxP2xLGP1_2D[i][j];
    }

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integrale_dxP2xLGP1xLGP1_2D(double *p, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(j = 0; j < 3; j++) {
    for(i = 0; i < 6; i++) {
      sum += p[i] * q[j] * fe_dxP2xLGP1xLGP1_2D[i][j][k];
      }
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_integrale_dxP2xLGP2xLGP1_2D(double *p, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(j = 0; j < 6; j++) {
    for(i = 0; i < 6; i++) {
      sum += p[i] * q[j] * fe_dxP2xLGP2xLGP1_2D[i][j][k];
      }
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex <double> fe_integrale_dxP2xLGP2xLGP1_2D(complex <double> *p, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  complex <double> sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(j = 0; j < 6; j++) {
    for(i = 0; i < 6; i++) {
      sum += p[i] * q[j] * fe_dxP2xLGP2xLGP1_2D[i][j][k];
      }
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex <double> fe_integrale_dxP2xLGP2xLGP1_2D(complex <double> *p, complex <double> *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  complex <double> sum,sum2;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(j = 0; j < 6; j++) {
    for(i = 0; i < 6; i++) {
      sum += p[i] * q[j] * fe_dxP2xLGP2xLGP1_2D[i][j][k];
      }
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex <double> fe_integrale_dxP2xLGP2xLGP1_2D(int i, complex <double> *q, complex <double> *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j,k;
  complex <double> sum,sum2;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(k = 0; k < 3; j++) {
    for(j = 0; j < 6; j++) {
      sum += q[j] * r[k] * fe_dxP2xLGP2xLGP1_2D[i][j][k];
      }
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integrale_dxP2xLGP2xLGP1_2D(int i, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(j = 0; j < 6; j++) {
    sum += q[j] * fe_dxP2xLGP2xLGP1_2D[i][j][k];
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integrale_dyP2xLGP2xLGP1_2D(double *p, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(j = 0; j < 6; j++) {
    for(i = 0; i < 6; i++) {
      sum += p[i] * q[j] * fe_dyP2xLGP2xLGP1_2D[i][j][k];
      }
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex <double> fe_integrale_dyP2xLGP2xLGP1_2D(complex <double> *p, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  complex <double> sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(j = 0; j < 6; j++) {
    for(i = 0; i < 6; i++) {
      sum += p[i] * q[j] * fe_dyP2xLGP2xLGP1_2D[i][j][k];
      }
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex <double> fe_integrale_dyP2xLGP2xLGP1_2D(complex <double> *p, complex <double> *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  complex <double> sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(j = 0; j < 6; j++) {
    for(i = 0; i < 6; i++) {
      sum += p[i] * q[j] * fe_dyP2xLGP2xLGP1_2D[i][j][k];
      }
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex <double> fe_integrale_dyP2xLGP2xLGP1_2D(int i, complex <double> *q, complex <double> *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j,k;
  complex <double> sum,sum2;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(k = 0; k < 3; j++) {
    for(j = 0; j < 6; j++) {
      sum += q[j] * r[k] * fe_dyP2xLGP2xLGP1_2D[i][j][k];
      }
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integrale_dyP2xLGP2xLGP1_2D(int i, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(j = 0; j < 6; j++) {
    sum += q[j] * fe_dyP2xLGP2xLGP1_2D[i][j][k];
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_integraleLGP2xNCP1xNCP1_2D(double *p, double *q, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(k = 0; k < 3; k++) {
    for(j = 0; j < 3; j++) {
      for(i = 0; i < 6; i++) {
        sum += p[i] * q[j] * r[k] * fe_LGP2xNCP1xNCP1_2D[i][j][k];
        }
      }
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP2xNCP1xNCP1_2D(double *p, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;
  for(j = 0; j < 3; j++) {
    for(i = 0; i < 6; i++) {
      sum += p[i] * q[j] * fe_LGP2xNCP1xNCP1_2D[i][j][k];
      }
    }
  return (sum);
    }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integrale_dxP2xNCP1_2D(double *p, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(i = 0; i < 6; i++) {
    sum += p[i] * fe_dxP2xNCP1_2D[i][j];
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integrale_dyP2xNCP1_2D(double *p, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(i = 0; i < 6; i++) {
    sum += p[i] * fe_dyP2xNCP1_2D[i][j];
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> T  fe_integrale_dxP2xLGP2xNCP1_2D_template(T *p, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
  see integration tabulations*/
#ifdef NOMINAL
  int i, j;
  T sum;
//#if 1
  sum = 0;

  for(j = 0; j < 6; j++) {
    for(i = 0; i < 6; i++) {
      sum += p[i] * q[j] * fe_dxP2xLGP2xNCP1_2D[i][j][k];
      }
    }
  return (sum);
#else
  T sum2=0;
  switch(k) {
    case 0:

      sum2+=3.3333333333332e-02*(p[0]*q[0]+p[0]*q[1]+p[0]*q[5]+p[5]*q[0]); // 0.033333
      sum2+=5.5555555555558e-03*(p[0]*q[2]+p[0]*q[4]); // 0.005556
      sum2+=5.5555555555556e-02*(p[0]*q[3]+p[2]*q[1]); // 0.055556
      sum2+=-3.3333333333331e-02*(p[1]*q[0]+p[3]*q[0]); // -0.033333
      sum2+=-8.8888888888889e-02*(p[1]*q[1]+p[5]*q[5]); // -0.088889
      sum2+=-5.5555555555556e-02*(p[1]*q[2]); // -0.055556
      sum2+=-1.3333333333333e-01*(p[1]*q[3]); // -0.133333
      sum2+=2.2222222222222e-02*(p[1]*q[4]); // 0.022222
      sum2+=-4.4444444444445e-02*(p[1]*q[5]+p[5]*q[1]); // -0.044444
      sum2+=-1.1102230246252e-16*(p[2]*q[0]+p[4]*q[0]+p[4]*q[1]+p[4]*q[2]+p[4]*q[3]+p[4]*q[4]+p[4]*q[5]); // -0.000000
      sum2+=5.0000000000000e-02*(p[2]*q[2]); // 0.050000
      sum2+=7.7777777777778e-02*(p[2]*q[3]); // 0.077778
      sum2+=-2.7777777777778e-02*(p[2]*q[4]); // -0.027778
      sum2+=1.1111111111111e-02*(p[2]*q[5]+p[5]*q[2]); // 0.011111
      sum2+=4.4444444444445e-02*(p[3]*q[1]); // 0.044444
      sum2+=-1.1111111111111e-02*(p[3]*q[2]); // -0.011111
      sum2+=1.7777777777778e-01*(p[3]*q[3]); // 0.177778
      sum2+=6.6666666666667e-02*(p[3]*q[4]); // 0.066667
      sum2+=8.8888888888889e-02*(p[3]*q[5]); // 0.088889
      sum2+=-1.7777777777778e-01*(p[5]*q[3]); // -0.177778
      sum2+=-6.6666666666667e-02*(p[5]*q[4]); // -0.066667
      break;

    case 1:

      sum2+=-5.0000000000000e-02*(p[0]*q[0]); // -0.050000
      sum2+=-5.5555555555555e-02*(p[0]*q[1]+p[2]*q[5]); // -0.055556
      sum2+=1.1102230246252e-16*(p[0]*q[2]+p[4]*q[0]+p[4]*q[1]+p[4]*q[2]+p[4]*q[3]+p[4]*q[4]+p[4]*q[5]); // 0.000000
      sum2+=-1.1111111111111e-02*(p[0]*q[3]+p[3]*q[0]); // -0.011111
      sum2+=2.7777777777778e-02*(p[0]*q[4]); // 0.027778
      sum2+=-7.7777777777778e-02*(p[0]*q[5]); // -0.077778
      sum2+=5.5555555555555e-02*(p[1]*q[0]); // 0.055556
      sum2+=8.8888888888889e-02*(p[1]*q[1]+p[3]*q[3]); // 0.088889
      sum2+=3.3333333333333e-02*(p[1]*q[2]+p[5]*q[2]); // 0.033333
      sum2+=4.4444444444444e-02*(p[1]*q[3]+p[3]*q[1]); // 0.044444
      sum2+=-2.2222222222222e-02*(p[1]*q[4]); // -0.022222
      sum2+=1.3333333333333e-01*(p[1]*q[5]); // 0.133333
      sum2+=-5.5555555555555e-03*(p[2]*q[0]+p[2]*q[4]); // -0.005556
      sum2+=-3.3333333333333e-02*(p[2]*q[1]+p[2]*q[2]+p[2]*q[3]+p[3]*q[2]); // -0.033333
      sum2+=6.6666666666667e-02*(p[3]*q[4]); // 0.066667
      sum2+=1.7777777777778e-01*(p[3]*q[5]); // 0.177778
      sum2+=1.1111111111111e-02*(p[5]*q[0]); // 0.011111
      sum2+=-4.4444444444444e-02*(p[5]*q[1]); // -0.044444
      sum2+=-8.8888888888889e-02*(p[5]*q[3]); // -0.088889
      sum2+=-6.6666666666667e-02*(p[5]*q[4]); // -0.066667
      sum2+=-1.7777777777778e-01*(p[5]*q[5]); // -0.177778
      break;

    case 2:

      sum2+=-5.0000000000000e-02*(p[0]*q[0]); // -0.050000
      sum2+=-7.7777777777778e-02*(p[0]*q[1]+p[1]*q[2]); // -0.077778
      sum2+=2.7777777777778e-02*(p[0]*q[2]); // 0.027778
      sum2+=-1.1111111111111e-02*(p[0]*q[3]+p[5]*q[0]+p[5]*q[2]); // -0.011111
      sum2+=1.1102230246252e-16*(p[0]*q[4]+p[1]*q[1]+p[1]*q[4]+p[2]*q[4]+p[3]*q[3]+p[3]*q[5]+p[4]*q[0]+p[4]*q[1]+p[4]*q[2]+p[4]*q[3]+p[4]*q[4]+p[4]*q[5]+p[5]*q[3]+p[5]*q[5]); // 0.000000
      sum2+=-5.5555555555556e-02*(p[0]*q[5]); // -0.055556
      sum2+=7.7777777777777e-02*(p[1]*q[0]+p[2]*q[1]); // 0.077778
      sum2+=-4.4444444444445e-02*(p[1]*q[3]+p[5]*q[1]); // -0.044444
      sum2+=4.4444444444445e-02*(p[1]*q[5]+p[3]*q[1]); // 0.044444
      sum2+=-2.7777777777778e-02*(p[2]*q[0]); // -0.027778
      sum2+=5.0000000000000e-02*(p[2]*q[2]); // 0.050000
      sum2+=5.5555555555556e-02*(p[2]*q[3]); // 0.055556
      sum2+=1.1111111111111e-02*(p[2]*q[5]+p[3]*q[0]+p[3]*q[2]); // 0.011111
      sum2+=-6.6666666666667e-02*(p[3]*q[4]); // -0.066667
      sum2+=6.6666666666667e-02*(p[5]*q[4]); // 0.066667
      break;
    }
  return (sum2);
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integrale_dxP2xLGP2xNCP1_2D(double *p, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
  sum=fe_integrale_dxP2xLGP2xNCP1_2D_template(p, q, k);
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integrale_dxP2xLGP2xNCP1_2D(complex<double> *p, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> sum;
  sum=fe_integrale_dxP2xLGP2xNCP1_2D_template(p, q, k);
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integrale_dxP2xLGP2xNCP1_2D(int i, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(j = 0; j < 6; j++) {
    sum += q[j] * fe_dxP2xLGP2xNCP1_2D[i][j][k];
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex <double> fe_integrale_dxP2xLGP2xNCP1_2D(int i, complex <double> *q, complex <double> *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j,k;
  complex <double> sum,sum2;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(k = 0; k < 3; k++) {
    for(j = 0; j < 6; j++) {
      sum += q[j] * r[k] * fe_dxP2xLGP2xNCP1_2D[i][j][k];
      }
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> T  fe_integrale_dyP2xLGP2xNCP1_2D_template(T *p, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
  see integration tabulations*/
#ifdef NOMINAL
  T sum;
  int i, j;
  sum = 0;

  for(j = 0; j < 6; j++) {
    for(i = 0; i < 6; i++) {
      sum += p[i] * q[j] * fe_dyP2xLGP2xNCP1_2D[i][j][k];
      }
    }
  return (sum);
#else
  T sum2=0;
  switch(k) {
    case 0:

      sum2+=3.3333333333332e-02*(p[0]*q[0]+p[0]*q[1]+p[0]*q[5]+p[1]*q[0]); // 0.033333
      sum2+=5.5555555555558e-03*(p[0]*q[2]+p[0]*q[4]); // 0.005556
      sum2+=5.5555555555556e-02*(p[0]*q[3]+p[4]*q[5]); // 0.055556
      sum2+=-8.8888888888889e-02*(p[1]*q[1]+p[5]*q[5]); // -0.088889
      sum2+=-6.6666666666667e-02*(p[1]*q[2]); // -0.066667
      sum2+=-1.7777777777778e-01*(p[1]*q[3]); // -0.177778
      sum2+=1.1111111111111e-02*(p[1]*q[4]+p[4]*q[1]); // 0.011111
      sum2+=-4.4444444444445e-02*(p[1]*q[5]+p[5]*q[1]); // -0.044444
      sum2+=0.0000000000000e+00*(p[2]*q[0]+p[2]*q[1]+p[2]*q[2]+p[2]*q[3]+p[2]*q[4]+p[2]*q[5]+p[4]*q[0]); // 0.000000
      sum2+=-3.3333333333333e-02*(p[3]*q[0]+p[5]*q[0]); // -0.033333
      sum2+=8.8888888888889e-02*(p[3]*q[1]); // 0.088889
      sum2+=6.6666666666667e-02*(p[3]*q[2]); // 0.066667
      sum2+=1.7777777777778e-01*(p[3]*q[3]); // 0.177778
      sum2+=-1.1111111111111e-02*(p[3]*q[4]); // -0.011111
      sum2+=4.4444444444445e-02*(p[3]*q[5]); // 0.044444
      sum2+=-2.7777777777778e-02*(p[4]*q[2]); // -0.027778
      sum2+=7.7777777777778e-02*(p[4]*q[3]); // 0.077778
      sum2+=5.0000000000000e-02*(p[4]*q[4]); // 0.050000
      sum2+=2.2222222222222e-02*(p[5]*q[2]); // 0.022222
      sum2+=-1.3333333333333e-01*(p[5]*q[3]); // -0.133333
      sum2+=-5.5555555555556e-02*(p[5]*q[4]); // -0.055556
      break;

    case 1:

      sum2+=-5.0000000000000e-02*(p[0]*q[0]); // -0.050000
      sum2+=-5.5555555555555e-02*(p[0]*q[1]); // -0.055556
      sum2+=1.1102230246252e-16*(p[0]*q[2]+p[1]*q[1]+p[1]*q[3]+p[2]*q[0]+p[2]*q[1]+p[2]*q[2]+p[2]*q[3]+p[2]*q[4]+p[2]*q[5]+p[3]*q[1]+p[3]*q[3]+p[4]*q[2]+p[5]*q[2]+p[5]*q[5]); // 0.000000
      sum2+=-1.1111111111111e-02*(p[0]*q[3]+p[1]*q[0]+p[1]*q[4]); // -0.011111
      sum2+=2.7777777777778e-02*(p[0]*q[4]); // 0.027778
      sum2+=-7.7777777777778e-02*(p[0]*q[5]+p[5]*q[4]); // -0.077778
      sum2+=6.6666666666667e-02*(p[1]*q[2]); // 0.066667
      sum2+=-4.4444444444444e-02*(p[1]*q[5]+p[5]*q[3]); // -0.044444
      sum2+=1.1111111111111e-02*(p[3]*q[0]+p[3]*q[4]+p[4]*q[1]); // 0.011111
      sum2+=-6.6666666666667e-02*(p[3]*q[2]); // -0.066667
      sum2+=4.4444444444444e-02*(p[3]*q[5]+p[5]*q[1]); // 0.044444
      sum2+=-2.7777777777778e-02*(p[4]*q[0]); // -0.027778
      sum2+=5.5555555555556e-02*(p[4]*q[3]); // 0.055556
      sum2+=5.0000000000000e-02*(p[4]*q[4]); // 0.050000
      sum2+=7.7777777777778e-02*(p[4]*q[5]+p[5]*q[0]); // 0.077778
      break;

    case 2:

      sum2+=-5.0000000000000e-02*(p[0]*q[0]); // -0.050000
      sum2+=-7.7777777777778e-02*(p[0]*q[1]); // -0.077778
      sum2+=2.7777777777778e-02*(p[0]*q[2]); // 0.027778
      sum2+=-1.1111111111111e-02*(p[0]*q[3]+p[3]*q[0]); // -0.011111
      sum2+=1.1102230246252e-16*(p[0]*q[4]+p[2]*q[0]+p[2]*q[1]+p[2]*q[2]+p[2]*q[3]+p[2]*q[4]+p[2]*q[5]); // 0.000000
      sum2+=-5.5555555555556e-02*(p[0]*q[5]+p[4]*q[1]); // -0.055556
      sum2+=1.1111111111111e-02*(p[1]*q[0]); // 0.011111
      sum2+=-1.7777777777778e-01*(p[1]*q[1]); // -0.177778
      sum2+=-6.6666666666667e-02*(p[1]*q[2]); // -0.066667
      sum2+=-8.8888888888889e-02*(p[1]*q[3]); // -0.088889
      sum2+=3.3333333333333e-02*(p[1]*q[4]+p[5]*q[4]); // 0.033333
      sum2+=-4.4444444444444e-02*(p[1]*q[5]); // -0.044444
      sum2+=1.7777777777778e-01*(p[3]*q[1]); // 0.177778
      sum2+=6.6666666666667e-02*(p[3]*q[2]); // 0.066667
      sum2+=8.8888888888889e-02*(p[3]*q[3]+p[5]*q[5]); // 0.088889
      sum2+=-3.3333333333333e-02*(p[3]*q[4]+p[4]*q[3]+p[4]*q[4]+p[4]*q[5]); // -0.033333
      sum2+=4.4444444444444e-02*(p[3]*q[5]+p[5]*q[3]); // 0.044444
      sum2+=-5.5555555555556e-03*(p[4]*q[0]+p[4]*q[2]); // -0.005556
      sum2+=5.5555555555556e-02*(p[5]*q[0]); // 0.055556
      sum2+=1.3333333333333e-01*(p[5]*q[1]); // 0.133333
      sum2+=-2.2222222222222e-02*(p[5]*q[2]); // -0.022222
      break;
    }
  return (sum2);
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integrale_dyP2xLGP2xNCP1_2D(double *p, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
  sum=fe_integrale_dyP2xLGP2xNCP1_2D_template(p, q, k);
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integrale_dyP2xLGP2xNCP1_2D(complex<double> *p, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> sum;
  sum=fe_integrale_dyP2xLGP2xNCP1_2D_template(p, q, k);
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integrale_dyP2xLGP2xNCP1_2D(int i, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(j = 0; j < 6; j++) {
    sum += q[j] * fe_dyP2xLGP2xNCP1_2D[i][j][k];
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex <double> fe_integrale_dyP2xLGP2xNCP1_2D(int i, complex <double> *q, complex <double> *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j,k;
  complex <double> sum,sum2;
/*------------------------------------------------------------------------------
  see integration tabulations*/
/*#ifdef NOMINAL*/
  sum = 0;

  for(k = 0; k < 3; k++) {
    for(j = 0; j < 6; j++) {
      sum += q[j] * r[k] * fe_dyP2xLGP2xNCP1_2D[i][j][k];
      }
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T,typename T1,typename T2,typename T3> T fe_integraleLGP1xLGP1xNCP1_2D_template(T1 *p, T2 *q, T3 *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
  see integration tabulations*/

#ifdef NOMINAL
  int i, j, k;
  T sum;
  
  sum = 0;

  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++) {
  T1 tmp;
      tmp=p[i] * q[j];
      for(k = 0; k < 3; k++)
        sum += tmp * r[k] * fe_LGP1xLGP1xNCP1_2D[i][j][k];
      }
  return (sum);
#endif

#ifdef OPTIMAL
  T sum2;
  
  sum2 = 0;

  sum2+=fe_LGP1xLGP1xNCP1_2D[0][0][0]*(p[0]*q[0]*r[0] +p[1]*q[1]*r[1] +p[2]*q[2]*r[2]); // -0.016667
  sum2+=fe_LGP1xLGP1xNCP1_2D[0][0][1]*(p[0]*q[0]*r[1] +p[0]*q[0]*r[2] +p[1]*q[1]*r[0] +p[1]*q[1]*r[2] +p[2]*q[2]*r[0] +p[2]*q[2]*r[1]); // 0.050000
  sum2+=fe_LGP1xLGP1xNCP1_2D[0][1][0]*(p[0]*q[1]*r[0] +p[0]*q[1]*r[1] +p[0]*q[2]*r[0] +p[0]*q[2]*r[2] +p[1]*q[0]*r[0] +p[1]*q[0]*r[1]
                                      +p[1]*q[2]*r[1] +p[1]*q[2]*r[2] +p[2]*q[0]*r[0] +p[2]*q[0]*r[2] +p[2]*q[1]*r[1] +p[2]*q[1]*r[2]); // 0.008333
  sum2+=fe_LGP1xLGP1xNCP1_2D[0][1][2]*(p[0]*q[1]*r[2] +p[0]*q[2]*r[1] +p[1]*q[0]*r[2] +p[1]*q[2]*r[0] +p[2]*q[0]*r[1] +p[2]*q[1]*r[0]); // 0.025000

  return (sum2);
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1xNCP1_2D(double *p, double *q, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=fe_integraleLGP1xLGP1xNCP1_2D_template<double,double,double,double>(p,q,r);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP1xNCP1_2D(double *p, double *q, complex<double> *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xLGP1xNCP1_2D_template<dcomplex,double,double,dcomplex>(p,q,r);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP1xNCP1_2D(complex<double> *p, double *q, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xLGP1xNCP1_2D_template<dcomplex,dcomplex,double,double>(p,q,r);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP1xNCP1_2D(complex<double> *p, double *q, complex<double> *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xLGP1xNCP1_2D_template<dcomplex,dcomplex,double,dcomplex>(p,q,r);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> T fe_integraleLGP1xLGP1xNCP1_2D_template(double *p, int j, T *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, k;
  T sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(i = 0; i < 3; i++) {
    for(k = 0; k < 3; k++) {
      sum += p[i] * r[k] * fe_LGP1xLGP1xNCP1_2D[i][j][k];
      }
    }

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1xNCP1_2D(double *p, int j, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=fe_integraleLGP1xLGP1xNCP1_2D_template(p,j,r);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP1xNCP1_2D(double *p, int j, complex<double> *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xLGP1xNCP1_2D_template(p,j,r);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_integraleLGP1xLGP1xNCP1_2D(double *p, double *q,int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum2;
/*------------------------------------------------------------------------------
  equivalent to the previous one when r=0 except for k where r=1*/
#ifdef NOMINAL
  int i, j;
//   double sum;
  sum = 0;

  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++) {
      sum += p[i] * q[j] * fe_LGP1xLGP1xNCP1_2D[i][j][k];
      }
//  return (sum);
#endif

#ifdef OPTIMAL
  sum2 = 0;

//   switch(k) {
//     case 0:
//
//       sum2+=fe_LGP1xLGP1xNCP1_2D[0][0][0]*(p[0]*q[0]); // -0.016667
//       sum2+=fe_LGP1xLGP1xNCP1_2D[0][1][0]*(p[0]*q[1]+p[0]*q[2]+p[1]*q[0]+p[2]*q[0]); // 0.008333
//       sum2+=fe_LGP1xLGP1xNCP1_2D[1][1][0]*(p[1]*q[1]+p[2]*q[2]); // 0.050000
//       sum2+=fe_LGP1xLGP1xNCP1_2D[1][2][0]*(p[1]*q[2]+p[2]*q[1]); // 0.025000
//       break;
//
//     case 1:
//
//       sum2+=fe_LGP1xLGP1xNCP1_2D[0][0][1]*(p[0]*q[0]+p[2]*q[2]); // 0.050000
//       sum2+=fe_LGP1xLGP1xNCP1_2D[0][1][1]*(p[0]*q[1]+p[1]*q[0]+p[1]*q[2]+p[2]*q[1]); // 0.008333
//       sum2+=fe_LGP1xLGP1xNCP1_2D[0][2][1]*(p[0]*q[2]+p[2]*q[0]); // 0.025000
//       sum2+=fe_LGP1xLGP1xNCP1_2D[1][1][1]*(p[1]*q[1]); // -0.016667
//       break;
//
//     case 2:
//
//       sum2+=fe_LGP1xLGP1xNCP1_2D[0][0][2]*(p[0]*q[0]+p[1]*q[1]); // 0.050000
//       sum2+=fe_LGP1xLGP1xNCP1_2D[0][1][2]*(p[0]*q[1]+p[1]*q[0]); // 0.025000
//       sum2+=fe_LGP1xLGP1xNCP1_2D[0][2][2]*(p[0]*q[2]+p[1]*q[2]+p[2]*q[0]+p[2]*q[1]); // 0.008333
//       sum2+=fe_LGP1xLGP1xNCP1_2D[2][2][2]*(p[2]*q[2]); // -0.016667
//       break;
//     }
//   return (sum2);

  switch(k) {
    case 0:

      sum2+=fe_LGP1xLGP1xNCP1_2D_c1*(p[0]*q[0]);                               // -0.016667
      sum2+=fe_LGP1xLGP1xNCP1_2D_c2*(p[0]*q[1]+p[0]*q[2]+p[1]*q[0]+p[2]*q[0]); // 0.008333
      sum2+=fe_LGP1xLGP1xNCP1_2D_c3*(p[1]*q[1]+p[2]*q[2]);                     // 0.050000
      sum2+=fe_LGP1xLGP1xNCP1_2D_c4*(p[1]*q[2]+p[2]*q[1]);                     // 0.025000
      break;

    case 1:

      sum2+=fe_LGP1xLGP1xNCP1_2D_c3*(p[0]*q[0]+p[2]*q[2]);                     // 0.050000
      sum2+=fe_LGP1xLGP1xNCP1_2D_c2*(p[0]*q[1]+p[1]*q[0]+p[1]*q[2]+p[2]*q[1]); // 0.008333
      sum2+=fe_LGP1xLGP1xNCP1_2D_c4*(p[0]*q[2]+p[2]*q[0]);                     // 0.025000
      sum2+=fe_LGP1xLGP1xNCP1_2D_c1*(p[1]*q[1]);                               // -0.016667
      break;

    case 2:

      sum2+=fe_LGP1xLGP1xNCP1_2D_c3*(p[0]*q[0]+p[1]*q[1]);                     // 0.050000
      sum2+=fe_LGP1xLGP1xNCP1_2D_c4*(p[0]*q[1]+p[1]*q[0]);                     // 0.025000
      sum2+=fe_LGP1xLGP1xNCP1_2D_c2*(p[0]*q[2]+p[1]*q[2]+p[2]*q[0]+p[2]*q[1]); // 0.008333
      sum2+=fe_LGP1xLGP1xNCP1_2D_c1*(p[2]*q[2]);                               // -0.016667
      break;
    }
  return (sum2);
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void fe_multiple_integraleLGP1xLGP1xNCP1_2D(double *p, double *q, double *sum)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  sum[0] =fe_LGP1xLGP1xNCP1_2D_c1*(p[0]*q[0]);                               // -0.016667
  sum[0]+=fe_LGP1xLGP1xNCP1_2D_c2*(p[0]*q[1]+p[0]*q[2]+p[1]*q[0]+p[2]*q[0]); // 0.008333
  sum[0]+=fe_LGP1xLGP1xNCP1_2D_c3*(p[1]*q[1]+p[2]*q[2]);                     // 0.050000
  sum[0]+=fe_LGP1xLGP1xNCP1_2D_c4*(p[1]*q[2]+p[2]*q[1]);                     // 0.025000

  sum[1] =fe_LGP1xLGP1xNCP1_2D_c3*(p[0]*q[0]+p[2]*q[2]);                     // 0.050000
  sum[1]+=fe_LGP1xLGP1xNCP1_2D_c2*(p[0]*q[1]+p[1]*q[0]+p[1]*q[2]+p[2]*q[1]); // 0.008333
  sum[1]+=fe_LGP1xLGP1xNCP1_2D_c4*(p[0]*q[2]+p[2]*q[0]);                     // 0.025000
  sum[1]+=fe_LGP1xLGP1xNCP1_2D_c1*(p[1]*q[1]);                               // -0.016667

  sum[2] =fe_LGP1xLGP1xNCP1_2D_c3*(p[0]*q[0]+p[1]*q[1]);                     // 0.050000
  sum[2]+=fe_LGP1xLGP1xNCP1_2D_c4*(p[0]*q[1]+p[1]*q[0]);                     // 0.025000
  sum[2]+=fe_LGP1xLGP1xNCP1_2D_c2*(p[0]*q[2]+p[1]*q[2]+p[2]*q[0]+p[2]*q[1]); // 0.008333
  sum[2]+=fe_LGP1xLGP1xNCP1_2D_c1*(p[2]*q[2]);                               // -0.016667

//   sum[0] =fe_LGP1xLGP1xNCP1_2D[0][0][0]*(p[0]*q[0]); // -0.016667
//   sum[0]+=fe_LGP1xLGP1xNCP1_2D[0][1][0]*(p[0]*q[1]+p[0]*q[2]+p[1]*q[0]+p[2]*q[0]); // 0.008333
//   sum[0]+=fe_LGP1xLGP1xNCP1_2D[1][1][0]*(p[1]*q[1]+p[2]*q[2]); // 0.050000
//   sum[0]+=fe_LGP1xLGP1xNCP1_2D[1][2][0]*(p[1]*q[2]+p[2]*q[1]); // 0.025000
//
//   sum[1] =fe_LGP1xLGP1xNCP1_2D[0][0][1]*(p[0]*q[0]+p[2]*q[2]); // 0.050000
//   sum[1]+=fe_LGP1xLGP1xNCP1_2D[0][1][1]*(p[0]*q[1]+p[1]*q[0]+p[1]*q[2]+p[2]*q[1]); // 0.008333
//   sum[1]+=fe_LGP1xLGP1xNCP1_2D[0][2][1]*(p[0]*q[2]+p[2]*q[0]); // 0.025000
//   sum[1]+=fe_LGP1xLGP1xNCP1_2D[1][1][1]*(p[1]*q[1]); // -0.016667
//
//   sum[2] =fe_LGP1xLGP1xNCP1_2D[0][0][2]*(p[0]*q[0]+p[1]*q[1]); // 0.050000
//   sum[2]+=fe_LGP1xLGP1xNCP1_2D[0][1][2]*(p[0]*q[1]+p[1]*q[0]); // 0.025000
//   sum[2]+=fe_LGP1xLGP1xNCP1_2D[0][2][2]*(p[0]*q[2]+p[1]*q[2]+p[2]*q[0]+p[2]*q[1]); // 0.008333
//   sum[2]+=fe_LGP1xLGP1xNCP1_2D[2][2][2]*(p[2]*q[2]); // -0.016667
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1xNCP1xNCP1_2D(double *p, double *q, double *r, double *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
#ifdef NOMINAL
  int i, j, k, l;
  double tmp,tmp2;
  sum = 0;

  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      tmp=p[i] * q[j];
      for(k = 0; k < 3; k++) {
        tmp2=tmp * r[k];
        for(l = 0; l < 3; l++)
          sum += tmp2 * s[l] * fe_LGP1xLGP1xNCP1xNCP1_2D[i][j][k][l];
        }
      }
    }
  return (sum);
#endif

#ifdef OPTIMAL
  sum = 0;

  sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][0][0][0]*(p[0]*q[0]*r[0]*s[0]+p[1]*q[1]*r[1]*s[1]+p[2]*q[2]*r[2]*s[2]); // 0.016667
  sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][0][0][1]*(p[0]*q[0]*r[0]*s[1]+p[0]*q[0]*r[0]*s[2]+p[0]*q[0]*r[1]*s[0]+p[0]*q[0]*r[2]*s[0]
                                             +p[1]*q[1]*r[0]*s[1]+p[1]*q[1]*r[1]*s[0]+p[1]*q[1]*r[1]*s[2]
                                             +p[1]*q[1]*r[2]*s[1]+p[2]*q[2]*r[0]*s[2]+p[2]*q[2]*r[1]*s[2]
                                             +p[2]*q[2]*r[2]*s[0]+p[2]*q[2]*r[2]*s[1]); // -0.016667
  sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][0][1][1]*(p[0]*q[0]*r[1]*s[1]+p[0]*q[0]*r[2]*s[2]+p[1]*q[1]*r[0]*s[0]+p[1]*q[1]*r[2]*s[2]
                                             +p[2]*q[2]*r[0]*s[0]+p[2]*q[2]*r[1]*s[1]); // 0.038889
  sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][0][1][2]*(p[0]*q[0]*r[1]*s[2]+p[0]*q[0]*r[2]*s[1]+p[1]*q[1]*r[0]*s[2]+p[1]*q[1]*r[2]*s[0]
                                             +p[2]*q[2]*r[0]*s[1]+p[2]*q[2]*r[1]*s[0]); // 0.027778
  sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][1][0][0]*(p[0]*q[1]*r[0]*s[0]+p[0]*q[1]*r[1]*s[1]+p[0]*q[2]*r[0]*s[0]+p[0]*q[2]*r[2]*s[2]
                                             +p[1]*q[0]*r[0]*s[0]+p[1]*q[0]*r[1]*s[1]+p[1]*q[2]*r[1]*s[1]
                                             +p[1]*q[2]*r[2]*s[2]+p[2]*q[0]*r[0]*s[0]+p[2]*q[0]*r[2]*s[2]
                                             +p[2]*q[1]*r[1]*s[1]+p[2]*q[1]*r[2]*s[2]); // 0.008333
  sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][1][0][1]*(p[0]*q[1]*r[0]*s[1]+p[0]*q[1]*r[1]*s[0]+p[0]*q[2]*r[0]*s[2]+p[0]*q[2]*r[2]*s[0]
                                             +p[1]*q[0]*r[0]*s[1]+p[1]*q[0]*r[1]*s[0]+p[1]*q[2]*r[1]*s[2]
                                             +p[1]*q[2]*r[2]*s[1]+p[2]*q[0]*r[0]*s[2]+p[2]*q[0]*r[2]*s[0]
                                             +p[2]*q[1]*r[1]*s[2]+p[2]*q[1]*r[2]*s[1]); // -0.002778
  sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][1][0][2]*(p[0]*q[1]*r[0]*s[2]+p[0]*q[1]*r[1]*s[2]+p[0]*q[1]*r[2]*s[0]+p[0]*q[1]*r[2]*s[1]
                                            +p[0]*q[2]*r[0]*s[1]+p[0]*q[2]*r[1]*s[0]+p[0]*q[2]*r[1]*s[2]
                                            +p[0]*q[2]*r[2]*s[1]+p[1]*q[0]*r[0]*s[2]+p[1]*q[0]*r[1]*s[2]
                                            +p[1]*q[0]*r[2]*s[0]+p[1]*q[0]*r[2]*s[1]+p[1]*q[2]*r[0]*s[1]
                                            +p[1]*q[2]*r[0]*s[2]+p[1]*q[2]*r[1]*s[0]+p[1]*q[2]*r[2]*s[0]
                                            +p[2]*q[0]*r[0]*s[1]+p[2]*q[0]*r[1]*s[0]+p[2]*q[0]*r[1]*s[2]
                                            +p[2]*q[0]*r[2]*s[1]+p[2]*q[1]*r[0]*s[1]+p[2]*q[1]*r[0]*s[2]
                                            +p[2]*q[1]*r[1]*s[0]+p[2]*q[1]*r[2]*s[0]); // 0.002778
  sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][1][2][2]*(p[0]*q[1]*r[2]*s[2]+p[0]*q[2]*r[1]*s[1]+p[1]*q[0]*r[2]*s[2]+p[1]*q[2]*r[0]*s[0]
                                            +p[2]*q[0]*r[1]*s[1]+p[2]*q[1]*r[0]*s[0]); // 0.019444
  return (sum);
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1xNCP1xNCP1_2D(double *p, double *q, double *r, int l)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;
#ifdef NOMINAL
  int i, j, k;
  double tmp;
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++) {
      tmp=p[i] * q[j];
      for(k = 0; k < 3; k++) {
        sum += tmp * r[k] * fe_LGP1xLGP1xNCP1xNCP1_2D[i][j][k][l];
        }
      }
#endif

#ifdef OPTIMAL
  switch(l) {
    case 0:

      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][0][0][0]*(p[0]*q[0]*r[0]); // 0.016667
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][0][1][0]*(p[0]*q[0]*r[1]+p[0]*q[0]*r[2]+p[1]*q[1]*r[1]+p[2]*q[2]*r[2]); // -0.016667
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][1][0][0]*(p[0]*q[1]*r[0]+p[0]*q[2]*r[0]+p[1]*q[0]*r[0]+p[2]*q[0]*r[0]); // 0.008333
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][1][1][0]*(p[0]*q[1]*r[1]+p[0]*q[2]*r[2]+p[1]*q[0]*r[1]+p[2]*q[0]*r[2]); // -0.002778
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][1][2][0]*(p[0]*q[1]*r[2]+p[0]*q[2]*r[1]+p[1]*q[0]*r[2]+p[1]*q[2]*r[1]
                                                 +p[1]*q[2]*r[2]+p[2]*q[0]*r[1]+p[2]*q[1]*r[1]+p[2]*q[1]*r[2]); // 0.002778
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[1][1][0][0]*(p[1]*q[1]*r[0]+p[2]*q[2]*r[0]); // 0.038889
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[1][1][2][0]*(p[1]*q[1]*r[2]+p[2]*q[2]*r[1]); // 0.027778
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[1][2][0][0]*(p[1]*q[2]*r[0]+p[2]*q[1]*r[0]); // 0.019444
      break;

    case 1:

      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][0][0][1]*(p[0]*q[0]*r[0]+p[1]*q[1]*r[0]+p[1]*q[1]*r[2]+p[2]*q[2]*r[2]); // -0.016667
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][0][1][1]*(p[0]*q[0]*r[1]+p[2]*q[2]*r[1]); // 0.038889
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][0][2][1]*(p[0]*q[0]*r[2]+p[2]*q[2]*r[0]); // 0.027778
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][1][0][1]*(p[0]*q[1]*r[0]+p[1]*q[0]*r[0]+p[1]*q[2]*r[2]+p[2]*q[1]*r[2] ); // -0.002778
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][1][1][1]*(p[0]*q[1]*r[1]+p[1]*q[0]*r[1]+p[1]*q[2]*r[1]+p[2]*q[1]*r[1]); // 0.008333
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][1][2][1]*(p[0]*q[1]*r[2]+p[0]*q[2]*r[0]+p[0]*q[2]*r[2]+p[1]*q[0]*r[2]
                                                 +p[1]*q[2]*r[0]+p[2]*q[0]*r[0]+p[2]*q[0]*r[2]+p[2]*q[1]*r[0]); // 0.002778
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][2][1][1]*(p[0]*q[2]*r[1]+p[2]*q[0]*r[1]); // 0.019444
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[1][1][1][1]*(p[1]*q[1]*r[1]); // 0.016667
      break;

    case 2:

      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][0][0][2]*(p[0]*q[0]*r[0]+p[1]*q[1]*r[1]+p[2]*q[2]*r[0]+p[2]*q[2]*r[1]); // -0.016667
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][0][1][2]*(p[0]*q[0]*r[1]+p[1]*q[1]*r[0]); // 0.027778
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][0][2][2]*(p[0]*q[0]*r[2]+p[1]*q[1]*r[2]); // 0.038889
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][1][0][2]*(p[0]*q[1]*r[0]+p[0]*q[1]*r[1]+p[0]*q[2]*r[1]+p[1]*q[0]*r[0]
                                                 +p[1]*q[0]*r[1]+p[1]*q[2]*r[0]+p[2]*q[0]*r[1]+p[2]*q[1]*r[0]); // 0.002778
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][1][2][2]*(p[0]*q[1]*r[2]+p[1]*q[0]*r[2]); // 0.019444
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][2][0][2]*(p[0]*q[2]*r[0]+p[1]*q[2]*r[1]+p[2]*q[0]*r[0]+p[2]*q[1]*r[1]); // -0.002778
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[0][2][2][2]*(p[0]*q[2]*r[2]+p[1]*q[2]*r[2]+p[2]*q[0]*r[2]+p[2]*q[1]*r[2]); // 0.008333
      sum+=fe_LGP1xLGP1xNCP1xNCP1_2D[2][2][2][2]*(p[2]*q[2]*r[2]); // 0.016667
      break;
    }
#endif

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1xLGP1xNCP1_2D(double *p, double *q, int k, double *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j/*, k*/, l;
  double tmp,sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++) {
      tmp=p[i] * q[j];
//      for(k = 0; k < 3; k++)
        for(l = 0; l < 3; l++)
          sum += tmp * s[l] * fe_LGP1xLGP1xLGP1xNCP1_2D[i][j][k][l];
      }

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1xLGP1xNCP1_2D(double *p, double *q, double *r, double *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l;
  double tmp,sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++) {
      tmp=p[i] * q[j];
      for(k = 0; k < 3; k++)
        for(l = 0; l < 3; l++)
          sum += tmp * r[k]* s[l] * fe_LGP1xLGP1xLGP1xNCP1_2D[i][j][k][l];
      }

  return (sum);
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// double fe_integraleLGP1xNCP1xNCP1_2D(double *p, double *q, double *r)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   double sum;
// /*------------------------------------------------------------------------------
//   see integration tabulations*/
//   sum = 0;
// #ifdef NOMINAL
//   int i, j, k;
//   for(i = 0; i < 3; i++)
//     for(j = 0; j < 3; j++)
//       for(k = 0; k < 3; k++)
//         sum += p[i] * q[j] * r[k] * fe_LGP1xNCP1xNCP1_2D[i][j][k];
// 
// #endif
// #ifdef OPTIMAL
// /**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// 
//   Development notes:
// 
//   Check : MANDATORY !!!
// 
//   Note:
// 
//   19/10/2009
// 
//     serious bug found : wrong optimal formula
// 
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
// //   sum+=fe_LGP1xNCP1xNCP1_2D[0][0][0]*(p[0]*q[0]*r[0] +p[1]*q[1]*r[1] +p[2]*q[2]*r[2]); // -0.016667
// //   sum+=fe_LGP1xNCP1xNCP1_2D[0][0][1]*(p[0]*q[0]*r[1] +p[0]*q[0]*r[2] +p[1]*q[1]*r[0] +p[1]*q[1]*r[2] +p[2]*q[2]*r[0] +p[2]*q[2]*r[1]); // 0.050000
// //   sum+=fe_LGP1xNCP1xNCP1_2D[0][1][0]*(p[0]*q[1]*r[0] +p[0]*q[1]*r[1] +p[0]*q[2]*r[0] +p[0]*q[2]*r[2] +p[1]*q[0]*r[0] +p[1]*q[0]*r[1]
// //                                      +p[1]*q[2]*r[1] +p[1]*q[2]*r[2] +p[2]*q[0]*r[0] +p[2]*q[0]*r[2] +p[2]*q[1]*r[1] +p[2]*q[1]*r[2]); // 0.008333
// //   sum+=fe_LGP1xNCP1xNCP1_2D[0][1][2]*(p[0]*q[1]*r[2] +p[0]*q[2]*r[1] +p[1]*q[0]*r[2] +p[1]*q[2]*r[0] +p[2]*q[0]*r[1] +p[2]*q[1]*r[0]); // 0.025000
//  sum+=fe_LGP1xNCP1xNCP1_2D[0][0][0]*(p[0]*q[0]*r[0]+p[0]*q[1]*r[2]+p[0]*q[2]*r[1]+p[1]*q[0]*r[2]+p[1]*q[1]*r[1]+p[1]*q[2]*r[0]
//                                     +p[2]*q[0]*r[1]+p[2]*q[1]*r[0]+p[2]*q[2]*r[2]); // 0.033333
//  sum+=fe_LGP1xNCP1xNCP1_2D[0][0][1]*(p[0]*q[0]*r[1]+p[0]*q[0]*r[2]+p[0]*q[1]*r[0]+p[0]*q[2]*r[0]+p[1]*q[0]*r[1]+p[1]*q[1]*r[0]
//                                     +p[1]*q[1]*r[2]+p[1]*q[2]*r[1]+p[2]*q[0]*r[2]+p[2]*q[1]*r[2]+p[2]*q[2]*r[0]+p[2]*q[2]*r[1]); // -0.016667
//  sum+=fe_LGP1xNCP1xNCP1_2D[0][1][1]*(p[0]*q[1]*r[1]+p[0]*q[2]*r[2]+p[1]*q[0]*r[0]+p[1]*q[2]*r[2]+p[2]*q[0]*r[0]+p[2]*q[1]*r[1]); // 0.066667
// #endif
// 
//   return (sum);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// double fe_integraleLGP1xNCP1xNCP1_2D(double *p, double *q, int k)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   double sum;
// /*------------------------------------------------------------------------------
//   see integration tabulations*/
//   sum = 0;
// #ifdef NOMINAL
//   int i, j;
//   for(i = 0; i < 3; i++)
//     for(j = 0; j < 3; j++)
//       sum += p[i] * q[j] * fe_LGP1xNCP1xNCP1_2D[i][j][k];
// #endif
// 
// #ifdef OPTIMAL
//   switch(k){
// /**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// 
//   Development notes:
// 
//   Check : MANDATORY !!!
// 
//   Note:
// 
//   19/10/2009
// 
//     serious bug found : wrong optimal formula
// 
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
// //     case 0:
// //
// //       sum+=fe_LGP1xNCP1xNCP1_2D[0][0][0]*(p[0]*q[0]); // -0.016667
// //       sum+=fe_LGP1xNCP1xNCP1_2D[0][1][0]*(p[0]*q[1]+p[0]*q[2]+p[1]*q[0]+p[2]*q[0]); // 0.008333
// //       sum+=fe_LGP1xNCP1xNCP1_2D[1][1][0]*(p[1]*q[1]+p[2]*q[2]); // 0.050000
// //       sum+=fe_LGP1xNCP1xNCP1_2D[1][2][0]*(p[1]*q[2]+p[2]*q[1]); // 0.025000
// //       break;
// //
// //     case 1:
// //
// //       sum+=fe_LGP1xNCP1xNCP1_2D[0][0][1]*(p[0]*q[0]+p[2]*q[2]); // 0.050000
// //       sum+=fe_LGP1xNCP1xNCP1_2D[0][1][1]*(p[0]*q[1]+p[1]*q[0]+p[1]*q[2]+p[2]*q[1]); // 0.008333
// //       sum+=fe_LGP1xNCP1xNCP1_2D[0][2][1]*(p[0]*q[2]+p[2]*q[0]); // 0.025000
// //       sum+=fe_LGP1xNCP1xNCP1_2D[1][1][1]*(p[1]*q[1]); // -0.016667
// //       break;
// //
// //     case 2:
// //
// //       sum+=fe_LGP1xNCP1xNCP1_2D[0][0][2]*(p[0]*q[0]+p[1]*q[1]); // 0.050000
// //       sum+=fe_LGP1xNCP1xNCP1_2D[0][1][2]*(p[0]*q[1]+p[1]*q[0]); // 0.025000
// //       sum+=fe_LGP1xNCP1xNCP1_2D[0][2][2]*(p[0]*q[2]+p[1]*q[2]+p[2]*q[0]+p[2]*q[1]); // 0.008333
// //       sum+=fe_LGP1xNCP1xNCP1_2D[2][2][2]*(p[2]*q[2]); // -0.016667
// //       break;
// 
// /**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// 
//   Development notes:
// 
//   Check : MANDATORY !!!
// 
//   Note:
// 
//   27/03/2010
// 
//     optimization
// 
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
// //     case 0:
// //
// //       sum+=fe_LGP1xNCP1xNCP1_2D[0][0][0]*(p[0]*q[0]+p[1]*q[2]+p[2]*q[1]); // 0.033333
// //       sum+=fe_LGP1xNCP1xNCP1_2D[0][1][0]*(p[0]*q[1]+p[0]*q[2]+p[1]*q[1]+p[2]*q[2]); // -0.016667
// //       sum+=fe_LGP1xNCP1xNCP1_2D[1][0][0]*(p[1]*q[0]+p[2]*q[0]); // 0.066667
// //       break;
// //
// //     case 1:
// //
// //       sum+=fe_LGP1xNCP1xNCP1_2D[0][0][1]*(p[0]*q[0]+p[1]*q[0]+p[1]*q[2]+p[2]*q[2]); // -0.016667
// //       sum+=fe_LGP1xNCP1xNCP1_2D[0][1][1]*(p[0]*q[1]+p[2]*q[1]); // 0.066667
// //       sum+=fe_LGP1xNCP1xNCP1_2D[0][2][1]*(p[0]*q[2]+p[1]*q[1]+p[2]*q[0]); // 0.033333
// //       break;
// //
// //     case 2:
// //
// //       sum+=fe_LGP1xNCP1xNCP1_2D[0][0][2]*(p[0]*q[0]+p[1]*q[1]+p[2]*q[0]+p[2]*q[1]); // -0.016667
// //       sum+=fe_LGP1xNCP1xNCP1_2D[0][1][2]*(p[0]*q[1]+p[1]*q[0]+p[2]*q[2]); // 0.033333
// //       sum+=fe_LGP1xNCP1xNCP1_2D[0][2][2]*(p[0]*q[2]+p[1]*q[2]); // 0.066667
// //       break;
// 
//     case 0:
// 
//       sum+=3.3333333333333e-02*(p[0]*q[0]+p[1]*q[2]+p[2]*q[1]); // 0.033333
//       sum+=-1.6666666666667e-02*(p[0]*q[1]+p[0]*q[2]+p[1]*q[1]+p[2]*q[2]); // -0.016667
//       sum+=6.6666666666667e-02*(p[1]*q[0]+p[2]*q[0]); // 0.066667
//       break;
// 
//     case 1:
// 
//       sum+=-1.6666666666667e-02*(p[0]*q[0]+p[1]*q[0]+p[1]*q[2]+p[2]*q[2]); // -0.016667
//       sum+=6.6666666666667e-02*(p[0]*q[1]+p[2]*q[1]); // 0.066667
//       sum+=3.3333333333333e-02*(p[0]*q[2]+p[1]*q[1]+p[2]*q[0]); // 0.033333
//       break;
// 
//     case 2:
// 
//       sum+=-1.6666666666667e-02*(p[0]*q[0]+p[1]*q[1]+p[2]*q[0]+p[2]*q[1]); // -0.016667
//       sum+=3.3333333333333e-02*(p[0]*q[1]+p[1]*q[0]+p[2]*q[2]); // 0.033333
//       sum+=6.6666666666667e-02*(p[0]*q[2]+p[1]*q[2]); // 0.066667
//       break;
//     }
// #endif
//   return (sum);
// }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> T fe_integraleLGP1xNCP1xNCP1_2D_template(double *p, T *q, T *r)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n;
  T sum;
/*----------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;
#ifdef NOMINAL
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      for(k = 0; k < 3; k++)
        sum += p[i] * q[j] * r[k] * fe_LGP1xNCP1xNCP1_2D[i][j][k];

#endif
#ifdef OPTIMAL
 sum+=fe_LGP1xNCP1xNCP1_2D[0][0][0]*(p[0]*q[0]*r[0]+p[0]*q[1]*r[2]+p[0]*q[2]*r[1]+p[1]*q[0]*r[2]+p[1]*q[1]*r[1]+p[1]*q[2]*r[0]
                                    +p[2]*q[0]*r[1]+p[2]*q[1]*r[0]+p[2]*q[2]*r[2]); // 0.033333
 sum+=fe_LGP1xNCP1xNCP1_2D[0][0][1]*(p[0]*q[0]*r[1]+p[0]*q[0]*r[2]+p[0]*q[1]*r[0]+p[0]*q[2]*r[0]+p[1]*q[0]*r[1]+p[1]*q[1]*r[0]
                                    +p[1]*q[1]*r[2]+p[1]*q[2]*r[1]+p[2]*q[0]*r[2]+p[2]*q[1]*r[2]+p[2]*q[2]*r[0]+p[2]*q[2]*r[1]); // -0.016667
 sum+=fe_LGP1xNCP1xNCP1_2D[0][1][1]*(p[0]*q[1]*r[1]+p[0]*q[2]*r[2]+p[1]*q[0]*r[0]+p[1]*q[2]*r[2]+p[2]*q[0]*r[0]+p[2]*q[1]*r[1]); // 0.066667
#endif

  return (sum);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_integraleLGP1xNCP1xNCP1_2D(double *p, double *q, double *r)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
  sum=fe_integraleLGP1xNCP1xNCP1_2D_template(p, q, r);
  return (sum);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> fe_integraleLGP1xNCP1xNCP1_2D(double *p, complex<double> *q, complex<double> *r)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> sum;
  sum=fe_integraleLGP1xNCP1xNCP1_2D_template(p, q, r);
  return (sum);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_integraleLGP1xNCP1xNCP1_2D(double *p, double *q, int k)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, l, m, n;
  double sum;
/*----------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;
#ifdef NOMINAL
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      sum += p[i] * q[j] * fe_LGP1xNCP1xNCP1_2D[i][j][k];
#endif

#ifdef OPTIMAL
  switch(k){
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check : MANDATORY !!!

  Note:

  27/03/2010

    optimization 

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
//     case 0: 
// 
//       sum+=fe_LGP1xNCP1xNCP1_2D[0][0][0]*(p[0]*q[0]+p[1]*q[2]+p[2]*q[1]); // 0.033333
//       sum+=fe_LGP1xNCP1xNCP1_2D[0][1][0]*(p[0]*q[1]+p[0]*q[2]+p[1]*q[1]+p[2]*q[2]); // -0.016667
//       sum+=fe_LGP1xNCP1xNCP1_2D[1][0][0]*(p[1]*q[0]+p[2]*q[0]); // 0.066667
//       break;
// 
//     case 1: 
// 
//       sum+=fe_LGP1xNCP1xNCP1_2D[0][0][1]*(p[0]*q[0]+p[1]*q[0]+p[1]*q[2]+p[2]*q[2]); // -0.016667
//       sum+=fe_LGP1xNCP1xNCP1_2D[0][1][1]*(p[0]*q[1]+p[2]*q[1]); // 0.066667
//       sum+=fe_LGP1xNCP1xNCP1_2D[0][2][1]*(p[0]*q[2]+p[1]*q[1]+p[2]*q[0]); // 0.033333
//       break;
// 
//     case 2: 
// 
//       sum+=fe_LGP1xNCP1xNCP1_2D[0][0][2]*(p[0]*q[0]+p[1]*q[1]+p[2]*q[0]+p[2]*q[1]); // -0.016667
//       sum+=fe_LGP1xNCP1xNCP1_2D[0][1][2]*(p[0]*q[1]+p[1]*q[0]+p[2]*q[2]); // 0.033333
//       sum+=fe_LGP1xNCP1xNCP1_2D[0][2][2]*(p[0]*q[2]+p[1]*q[2]); // 0.066667
//       break;

    case 0:

      sum+=3.3333333333333e-02*(p[0]*q[0]+p[1]*q[2]+p[2]*q[1]); // 0.033333
      sum+=-1.6666666666667e-02*(p[0]*q[1]+p[0]*q[2]+p[1]*q[1]+p[2]*q[2]); // -0.016667
      sum+=6.6666666666667e-02*(p[1]*q[0]+p[2]*q[0]); // 0.066667
      break;

    case 1:

      sum+=-1.6666666666667e-02*(p[0]*q[0]+p[1]*q[0]+p[1]*q[2]+p[2]*q[2]); // -0.016667
      sum+=6.6666666666667e-02*(p[0]*q[1]+p[2]*q[1]); // 0.066667
      sum+=3.3333333333333e-02*(p[0]*q[2]+p[1]*q[1]+p[2]*q[0]); // 0.033333
      break;

    case 2:

      sum+=-1.6666666666667e-02*(p[0]*q[0]+p[1]*q[1]+p[2]*q[0]+p[2]*q[1]); // -0.016667
      sum+=3.3333333333333e-02*(p[0]*q[1]+p[1]*q[0]+p[2]*q[2]); // 0.033333
      sum+=6.6666666666667e-02*(p[0]*q[2]+p[1]*q[2]); // 0.066667
      break;
    }
#endif
  return (sum);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_integraleNCP1xNCP1_2D(double *p, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/

  sum = 0;

  for(j = 0; j < 3; j++)
     sum += p[j] * q[j] * fe_NCP1xNCP1_2D[j][j];

  return (sum);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T1, typename T2> T1 fe_integraleNCP1xNCP1xNCP1_2D_template(T2 *p, T1 *q, T1 *r)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n;
  T1 sum,sum2;
  T1 tmp1,tmp2,tmp3;

  sum = 0;

  for(k = 0; k < 3; k++) {
    for(j = 0; j < 3; j++) {
      for(i = 0; i < 3; i++) {
        sum += p[i] * q[j] * r[k] * fe_NCP1xNCP1xNCP1_2D[i][j][k];
        }
      }
    }
  return (sum);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_integraleNCP1xNCP1xNCP1_2D(double *p, double *q, double *r)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum=fe_integraleNCP1xNCP1xNCP1_2D_template(p, q, r);
  return (sum);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> fe_integraleNCP1xNCP1xNCP1_2D(complex<double> *p, complex<double> *q, complex<double> *r)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> sum=fe_integraleNCP1xNCP1xNCP1_2D_template(p, q, r);
  return (sum);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> fe_integraleNCP1xNCP1xNCP1_2D(double *p, complex<double> *q, complex<double> *r)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> sum=fe_integraleNCP1xNCP1xNCP1_2D_template(p, q, r);
  return (sum);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_integraleLGP1xLGP1xLGP1_1D_nonoptimal(double *p, double *q, double *r) /// JUST SLOWER THAN fe_integraleLGP1xLGP1xLGP1_1D

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(i = 0; i < 2; i++)
    for(j = 0; j < 2; j++)
      for(k = 0; k < 2; k++)
        sum += p[i] * q[j] * r[k] * fe_LGP1xLGP1xLGP1_1D[i][j][k];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T1,typename T2,typename T3> T1 fe_integraleLGP1xLGP1xLGP1_1D_template(T1 *p, T2 *q, T3 *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T1 sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/

  sum  =p[0]*(q[0]*r[1]+q[1]*(r[0]+r[1]));
  sum +=p[1]*(q[0]*(r[0]+r[1])+q[1]*r[0]);
  sum /=12.0;
  sum +=(p[0]*q[0]*r[0]+p[1]*q[1]*r[1])/4.0;

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1xLGP1_1D(double *p, double *q, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=fe_integraleLGP1xLGP1xLGP1_1D_template(p,q,r);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP1xLGP1_1D(complex<double> *p, double *q, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xLGP1xLGP1_1D_template(p,q,r);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> fe_integraleLGP1xLGP1xLGP1_1D(double *p, double *q, complex<double> *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/

  sum  =p[0]*(q[0]*r[1]+q[1]*(r[0]+r[1]));
  sum +=p[1]*(q[0]*(r[0]+r[1])+q[1]*r[0]);
  sum /=12.0;
  sum +=(p[0]*q[0]*r[0]+p[1]*q[1]*r[1])/4.0;


  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP1xLGP1_1D(complex<double> *p, complex<double> *q, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xLGP1xLGP1_1D_template(p,q,r);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1xLGP1_1D(double *p, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(i = 0; i < 2; i++)
    for(j = 0; j < 2; j++)
      sum += p[i] * q[j] * fe_LGP1xLGP1xLGP1_1D[i][j][k];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1xLGP1xLGP1_1D_nonoptimal(double *p, double *q, double *r, double *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(l = 0; l < 2; l++)
    for(k = 0; k < 2; k++)
      for(j = 0; j < 2; j++)
        for(i = 0; i < 2; i++)
          sum += p[i] * q[j] * r[k] * s[l] * fe_LGP1xLGP1xLGP1xLGP1_1D[i][j][k][l];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1xLGP1xLGP1_1D(double *p, double *q, double *r, double *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(i = 0; i < 2; i++) {
    for(j = 0; j < 2; j++) {
      for(k = 0; k < 2; k++) {
        for(l = 0; l < 2; l++) {
          sum += p[i] * q[j] * r[k] * s[l] * fe_LGP1xLGP1xLGP1xLGP1_1D[i][j][k][l];
          }
        }
      }
    }

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1xLGP1xLGP1_1D(double *p, double *q, double *r, int l)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(i = 0; i < 2; i++) {
    for(j = 0; j < 2; j++) {
      for(k = 0; k < 2; k++) {
        sum += p[i] * q[j] * r[k] *  fe_LGP1xLGP1xLGP1xLGP1_1D[i][j][k][l];
        }
      }
    }

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1xLGP1xLGP1xLGP1_1D(double *p, double *q, double *r, double *s, double *t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(i = 0; i < 2; i++) {
    for(j = 0; j < 2; j++) {
      for(k = 0; k < 2; k++) {
        for(l = 0; l < 2; l++) {
          for(m = 0; m < 2; m++) {
            sum += p[i] * q[j] * r[k] * s[l] * t[m] * fe_LGP1xLGP1xLGP1xLGP1xLGP1_1D[i][j][k][l][m];
            }
          }
        }
      }
    }

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1xLGP1xLGP1xLGP1_1D(double *p, double *q, double *r, double *s, int m)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(i = 0; i < 2; i++) {
    for(j = 0; j < 2; j++) {
      for(k = 0; k < 2; k++) {
        for(l = 0; l < 2; l++) {
          sum += p[i] * q[j] * r[k] * s[l] * fe_LGP1xLGP1xLGP1xLGP1xLGP1_1D[i][j][k][l][m];
          }
        }
      }
    }

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T1,typename T2> T1 fe_integraleLGP1xLGP1_1D_template(T1 *p, T2 *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  T1 sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(j = 0; j < 2; j++)
    for(i = 0; i < 2; i++)
      sum += p[i] * q[j] * fe_LGP1xLGP1_1D[i][j];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1_1D(double *p, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=fe_integraleLGP1xLGP1_1D_template(p,q);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP1_1D(complex<double> *p, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xLGP1_1D_template(p,q);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1xLGP1_1D(complex<double> *p, complex<double> *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1xLGP1_1D_template(p,q);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> T fe_integraleLGP1_1D_template(T *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  T sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(i = 0; i < 2; i++)
    sum += p[i] * 0.5;

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1_1D(double *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=fe_integraleLGP1_1D_template(p);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

complex<double> fe_integraleLGP1_1D(complex<double> *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=fe_integraleLGP1_1D_template(p);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1_1D(double *p, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(i = 0; i < 2; i++)
    sum += p[i] * fe_LGP1xLGP1_1D[i][j];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1xNCP1_1D(double *p, double *q, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(k = 0; k < 1; k++)
    for(j = 0; j < 2; j++)
      for(i = 0; i < 2; i++)
        sum += p[i] * q[j] * r[k] * fe_LGP1xLGP1xNCP1_1D[i][j][k];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double integraleL1xL1(double *p, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/

  sum = 0;
  for(j = 0; j < 3; j++)
    sum += p[j] * q[j] / (double) 12.;
  for(j = 0; j < 3; j++)
    for(i = 0; i < 3; i++)
      if(i != j)
        sum += p[i] * q[j] / (double) 24.;
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleP0_1D(double *x, double *y, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/

  sum = 0;
  for(i = 0; i < n - 1; i++) {
    sum += (x[i + 1] - x[i]) * y[i];
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleP1_1D(double *x, double *y, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;
  for(i = 0; i < n - 1; i++) {
    sum += (x[i + 1] - x[i]) * (y[i] + y[i+1]) / (double) 2.;
    }
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_euclidianLGP1xLGP1_2D(double *p, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/

  sum=(p[0] * q[0] + p[1] * q[1] + p[2] * q[2])/6.0;
/*
  int j;
  sum = 0;
  for(j = 0; j < 3; j++)
    sum += p[j] * q[j] / (double) 6.;
*/
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_euclidianLGP1xLGP1_2D(double *p, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/

  sum = p[j] / (double) 6.;
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_euclidianLGP1xLGP1xLGP1_2D(double *p, double *q, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/

  sum = 0;
  for(j = 0; j < 3; j++)
    sum += p[j] * q[j] * r[j] / (double) 6.;
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_euclidianLGP1xLGP1xLGP1_2D(double *p, double *q, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/

  sum = p[j] * q[j] / (double) 6.;
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_euclidianLGP1xLGP1xLGP1xLGP1_2D(double *p, double *q, double *r, double *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/

  sum = 0;
  for(j = 0; j < 3; j++)
    sum += p[j] * q[j] * r[j] * s[j] / (double) 6.;
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_euclidianLGP1xLGP1xLGP1xLGP1_2D(double *p, double *q, double *r, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/

  sum = 0;
  for(j = 0; j < 3; j++)
    sum += p[j] * q[j] * r[j] / (double) 6.;
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

inline  double fe_euclidianLGP1xNCP1_2D(double *p, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum,tmp[3];
/*------------------------------------------------------------------------------
  see integration tabulations*/

  tmp[0]=q[1]+q[2]-q[0];
  tmp[1]=q[0]+q[2]-q[1];
  tmp[2]=q[0]+q[1]-q[2];

  sum=(p[0] * tmp[0] + p[1] * tmp[1] + p[2] * tmp[2])/6.0;
/*
  sum = 0;
  for(j = 0; j < 3; j++)
    sum += p[j] * tmp[j] / (double) 6.;
*/
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_euclidianLGP1xLGP1xNCP1_2D(double *p, double *q, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j;
  double sum,tmp[3];
/*------------------------------------------------------------------------------
  see integration tabulations*/

  tmp[0]=r[1]+r[2]-r[0];
  tmp[1]=r[0]+r[2]-r[1];
  tmp[2]=r[0]+r[1]-r[2];

  sum = 0;
  for(j = 0; j < 3; j++)
    sum += p[j] * q[j] * tmp[j] / (double) 6.;
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_euclidianLGP1xLGP1xNCP1_2D(double *p, int j, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum,tmp[3];
/*------------------------------------------------------------------------------
  see integration tabulations*/

  tmp[0]=r[1]+r[2]-r[0];
  tmp[1]=r[0]+r[2]-r[1];
  tmp[2]=r[0]+r[1]-r[2];

  sum = p[j] * tmp[j] / (double) 6.;
  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1_3D(double *p, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(j = 0; j < 3; j++)
    for(i = 0; i < 3; i++)
      sum += p[i] * q[j] * fe_LGP1xLGP1_3D[i][j];

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1xLGP1_3D(double *p, double *q, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(k = 0; k < 6; k++)
    for(j = 0; j < 6; j++)
      for(i = 0; i < 6; i++)
        sum += p[i] * q[j] * r[k] * fe_LGP1xLGP1xLGP1_3D[i][j][k];

  return (sum);
}


#if 0 /* thankfully unused */
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_integraleLGP1xLGP1xLGP1_3D(double *p, int j, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//   int i;
  double sum;
/*------------------------------------------------------------------------------
  see integration tabulations*/
  sum = 0;

  for(j = 0; j < 6; j++)
    sum += p[i] * fe_LGP1xLGP1xLGP1_3D[i][j][k];

  return (sum);
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_sproductLGP1xLGP1_2D(double *p, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
/*------------------------------------------------------------------------------
  */
  sum = fe_sproduct_ptrLGP1xLGP1_2D_a(p,q);

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_sproductLGP1xLGP1_2D(double *p, int j)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
/*------------------------------------------------------------------------------
  */
  sum = fe_sproduct_ptrLGP1xLGP1_2D_b(p,j);

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_sproductLGP1xLGP1xLGP1_2D(double *p, double *q, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
/*------------------------------------------------------------------------------
  */
  sum = fe_sproduct_ptrLGP1xLGP1xLGP1_2D_a(p,q,r);

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_sproductLGP1xLGP1xLGP1_2D(double *p, double *q, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
/*------------------------------------------------------------------------------
  */
  sum = fe_sproduct_ptrLGP1xLGP1xLGP1_2D_b(p,q,k);

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_sproductLGP1xLGP1xLGP1xLGP1_2D(double *p, double *q, double *r, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
/*------------------------------------------------------------------------------
  */
  sum = fe_sproduct_ptrLGP1xLGP1xLGP1xLGP1_2D(p,q,r,k);

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_sproductLGP1xNCP1_2D(double *p, double *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
/*------------------------------------------------------------------------------
  */
  sum = fe_sproduct_ptrLGP1xNCP1_2D(p,q);

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_sproductLGP1xLGP1xNCP1_2D(double *p, double *q, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
/*------------------------------------------------------------------------------
  */
  sum = fe_sproduct_ptrLGP1xLGP1xNCP1_2D_a(p,q,r);

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_sproductLGP1xLGP1xNCP1_2D(double *p, int j, double *r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sum;
/*------------------------------------------------------------------------------
  */
  sum = fe_sproduct_ptrLGP1xLGP1xNCP1_2D_b(p,j,r);

  return (sum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_sproduct_init(int sproduct)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  switch (sproduct) {
    case NQUAD:
      fe_sproduct_ptrLGP1xLGP1_2D_a=fe_euclidianLGP1xLGP1_2D;
      fe_sproduct_ptrLGP1xLGP1_2D_b=fe_euclidianLGP1xLGP1_2D;

      fe_sproduct_ptrLGP1xLGP1xLGP1_2D_a=fe_euclidianLGP1xLGP1xLGP1_2D;
      fe_sproduct_ptrLGP1xLGP1xLGP1_2D_b=fe_euclidianLGP1xLGP1xLGP1_2D;

      fe_sproduct_ptrLGP1xNCP1_2D=fe_euclidianLGP1xNCP1_2D;

      fe_sproduct_ptrLGP1xLGP1xNCP1_2D_a=fe_euclidianLGP1xLGP1xNCP1_2D;
      fe_sproduct_ptrLGP1xLGP1xNCP1_2D_b=fe_euclidianLGP1xLGP1xNCP1_2D;

      fe_sproduct_ptrLGP1xLGP1xLGP1xLGP1_2D=fe_euclidianLGP1xLGP1xLGP1xLGP1_2D;
      break;

    case INTGL:
      fe_sproduct_ptrLGP1xLGP1_2D_a=fe_integraleLGP1xLGP1_2D;
      fe_sproduct_ptrLGP1xLGP1_2D_b=fe_integraleLGP1xLGP1_2D;

      fe_sproduct_ptrLGP1xLGP1xLGP1_2D_a=fe_integraleLGP1xLGP1xLGP1_2D;
      fe_sproduct_ptrLGP1xLGP1xLGP1_2D_b=fe_integraleLGP1xLGP1xLGP1_2D;

      fe_sproduct_ptrLGP1xNCP1_2D=fe_integraleLGP1xNCP1_2D;

      fe_sproduct_ptrLGP1xLGP1xNCP1_2D_a=fe_integraleLGP1xLGP1xNCP1_2D;
      fe_sproduct_ptrLGP1xLGP1xNCP1_2D_b=fe_integraleLGP1xLGP1xNCP1_2D;

      fe_sproduct_ptrLGP1xLGP1xLGP1xLGP1_2D=fe_integraleLGP1xLGP1xLGP1xLGP1_2D;

      break;

     default:
       break;
     }
}

