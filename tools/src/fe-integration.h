
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief finite element integration functions declarations
*/
/*----------------------------------------------------------------------------*/

#ifndef FE_INTEGRATION_H
#define FE_INTEGRATION_H

#include "fe-classes.h" /* for mesh_t */

extern int fe_sproduct_init(int);

extern double fe_sproductLGP1xLGP1_2D      (double *, double *);
extern double fe_sproductLGP1xLGP1_2D      (double *, int);

extern double fe_sproductLGP1xLGP1xLGP1_2D (double *, double *, double *);
extern double fe_sproductLGP1xLGP1xLGP1_2D (double *, double *, int);

extern double fe_sproductLGP1xLGP1xLGP1xLGP1_2D (double *, double *,double *, int);

extern double fe_sproductLGP1xNCP1_2D      (double *, double *);
extern double fe_sproductLGP1xNCP1_2D      (double *, double *);

extern double fe_sproductLGP1xLGP1xNCP1_2D (double *, double *, double *);
extern double fe_sproductLGP1xLGP1xNCP1_2D (double *, double *, int);
extern double fe_sproductLGP1xLGP1xNCP1_2D (double *, int     , double *);


#undef PREFIX
#ifdef MAIN_FE_INTEGRATION_SOURCE
#define PREFIX
#else
#define PREFIX extern
#endif

/**-----------------------------------------------------------------------
1D integration coefficients*/

PREFIX double fe_LGP1xLGP1_1D [2][2];
PREFIX double fe_LGP1xLGP1xLGP1_1D [2][2][2];
PREFIX double fe_LGP1xLGP1xNCP1_1D [2][2][2];
PREFIX double fe_LGP1xLGP1xLGP1xLGP1_1D [2][2][2][2];
PREFIX double fe_LGP1xLGP1xLGP1xLGP1xLGP1_1D [2][2][2][2][2];

/**-----------------------------------------------------------------------
2D integration coefficients*/

PREFIX double fe_LGP1_2D [3];
PREFIX double fe_LGP1xLGP1_2D [3][3];
PREFIX double fe_LGP1xLGP1xLGP1_2D [3][3][3];
PREFIX double fe_LGP1xLGP1xLGP1xLGP1_2D [3][3][3][3];

PREFIX double fe_NCP1_2D[3];
PREFIX double fe_NCP1xNCP1_2D[3][3];
PREFIX double fe_NCP1xNCP1xNCP1_2D[3][3][3];

PREFIX double fe_LGP1xNCP1_2D[3][3];
PREFIX double fe_LGP1xLGP1xNCP1_2D[3][3][3];
PREFIX double fe_LGP1xLGP1xNCP1xNCP1_2D[3][3][3][3];
PREFIX double fe_LGP1xLGP1xLGP1xNCP1_2D[3][3][3][3];
PREFIX double test_fe_LGP1xLGP1xNCP1_2D[3][3][3];
PREFIX double fe_LGP1xNCP1xNCP1_2D[3][3][3];

//double fe_LGP1xLGP1xNCP1[3][3][3];

PREFIX double fe_LGP2_2D [6];
PREFIX double fe_LGP2xLGP2_2D [6][6];
PREFIX double fe_LGP2xLGP2xLGP2_2D [6][6][6];

PREFIX double fe_LGP2xNCP1_2D[6][3];
PREFIX double fe_LGP2xNCP1xNCP1_2D[6][3][3];

PREFIX double fe_LGP2xLGP1_2D[6][3];
PREFIX double fe_LGP2xLGP1xLGP1_2D[6][3][3];

PREFIX double fe_LGP2xLGP2xLGP1_2D[6][6][3];
PREFIX double fe_LGP2xLGP2xNCP1_2D[6][6][3];

PREFIX double fe_dxP2xLGP2xLGP1_2D[6][6][3];
PREFIX double fe_dyP2xLGP2xLGP1_2D[6][6][3];

PREFIX double fe_dxP2xLGP2xLGP2_2D[6][6][6];
PREFIX double fe_dyP2xLGP2xLGP2_2D[6][6][6];

PREFIX double fe_dxP2xLGP1_2D[6][3];
PREFIX double fe_dyP2xLGP1_2D[6][3];

PREFIX double fe_dxP2xLGP1xLGP1_2D[6][3][3];
PREFIX double fe_dyP2xLGP1xLGP1_2D[6][3][3];

PREFIX double fe_dxP2xLGP2xNCP1_2D[6][6][3];
PREFIX double fe_dyP2xLGP2xNCP1_2D[6][6][3];

PREFIX double fe_dxP2xLGP2xNCP1xNCP1_2D[6][6][3][3];
PREFIX double fe_dyP2xLGP2xNCP1xNCP1_2D[6][6][3][3];

PREFIX double fe_dxP2xNCP1_2D[6][3];
PREFIX double fe_dyP2xNCP1_2D[6][3];

PREFIX double fe_QP1_2D [4];
PREFIX double fe_QP1xQP1_2D [4][4];
PREFIX double fe_QP1xQP1xQP1_2D [4][4][4];

//commented out definition
//extern double integraleP1xP1_NQ(double *p, double *q);

/**-----------------------------------------------------------------------
3D integration coefficients*/

PREFIX double fe_LGP1_3D [6];
PREFIX double fe_LGP1xLGP1_3D [6][6];
PREFIX double fe_LGP1xLGP1xLGP1_3D [6][6][6];

#undef PREFIX

/**-----------------------------------------------------------------------------
Integrale 2D (triangle)*/
extern double fe_integraleLGP1_2D(double *p);
extern complex<double> fe_integraleLGP1_2D(complex<double> *p);
extern double fe_integraleLGP1_2D  (int i=-1);
#if 0
//On one hand, these functions give the same results ...
#define fe_integraleNCP1_2D fe_integraleLGP1_2D
//... but on the other, this might confuse people when expanding to products involving NCP1 discretisations.
#else
extern float fe_integraleNCP1_2D(float *p);
extern double fe_integraleNCP1_2D(double *p);
extern complex<float> fe_integraleNCP1_2D(complex<float> *p);
extern complex<double> fe_integraleNCP1_2D(complex<double> *p);
extern double fe_integraleNCP1_2D  (int i=-1);
#endif

extern double fe_integraleNCP1xNCP1_2D  (double *, double *);
extern complex<double> fe_integraleNCP1xNCP1_2D  (complex<double> *, complex<double> *);

extern double fe_integraleNCP1xNCP1_2D  (double *, int);
extern complex<double> fe_integraleNCP1xNCP1_2D  (complex<double> *, int);

extern double          fe_integraleNCP1xNCP1xNCP1_2D (double *, double *, double *);
extern complex<double> fe_integraleNCP1xNCP1xNCP1_2D (complex<double> *, complex<double> *, complex<double> *);
extern complex<double> fe_integraleNCP1xNCP1xNCP1_2D (double *, complex<double> *, complex<double> *);

extern double fe_integraleLGP1xLGP0_2D  (double *, double *);

extern double fe_integraleLGP1xNCP1_2D(double *, double *);
extern double fe_integraleLGP1xNCP1_2D(double *p, int j);
extern double fe_integraleLGP1xNCP1_2D(int i, double *q);

extern complex<double> fe_integraleLGP1xNCP1_2D(double *p, complex<double> *q);
extern complex<double> fe_integraleLGP1xNCP1_2D(complex<double> *p, double *q);
extern complex<double> fe_integraleLGP1xNCP1_2D(complex<double> *p, complex<double> *q);
extern complex<double> fe_integraleLGP1xNCP1_2D(complex<double> *p, int j);
extern complex<double> fe_integraleLGP1xNCP1_2D(int i, complex<double> *q);

extern double fe_integraleLGP1xLGP1_2D  (double *, double *);
extern double fe_integraleLGP1xLGP1_2D  (double *, int);
extern double fe_integraleLGP1xLGP1_2D  (int, int);

extern complex<double> fe_integraleLGP1xLGP1_2D(complex<double> *, complex<double> *);
extern complex<double> fe_integraleLGP1xLGP1_2D(double *, complex<double> *);
extern complex<double> fe_integraleLGP1xLGP1_2D(complex<double> *,double *);
extern complex<double> fe_integraleLGP1xLGP1_2D(complex<double> *, int);

extern double integraleP1xP1    (double *, double *);

extern double fe_integraleLGP1xLGP1xLGP1_2D (double *, double *, double *);
extern double fe_integraleLGP1xLGP1xLGP1_2D (double *, double *, int);

extern complex<double> fe_integraleLGP1xLGP1xLGP1_2D(double *p, complex<double> *q, complex<double> *r);
extern complex<double> fe_integraleLGP1xLGP1xLGP1_2D(double *p, complex<double> *q, int);

extern complex<double> fe_integraleLGP1xLGP1xLGP1_2D(complex<double> *p, double *q, double *r);
extern complex<double> fe_integraleLGP1xLGP1xLGP1_2D(complex<double> *p, complex<double> *q, double *r);
extern complex<double> fe_integraleLGP1xLGP1xLGP1_2D(complex<double> *p, complex<double> *q, complex<double> *r);
extern complex<double> fe_integraleLGP1xLGP1xLGP1_2D(complex<double> *p, complex<double> *q, int);

extern double fe_integraleLGP1xLGP1xLGP1xLGP1_2D (double *, double *,double *, int);

extern complex<double> fe_integraleLGP1xLGP1xLGP0_2D(double *p, complex<double> *q, complex<double> *r);

extern double fe_integraleLGP1xLGP1xNCP1_2D (double *, double *, double *);
extern double fe_integraleLGP1xLGP1xNCP1_2D (double *, double *, int);
extern double fe_integraleLGP1xLGP1xNCP1_2D (double *, int, double *);
extern double fe_integraleLGP1xLGP1xNCP1_2D (double *, int, int);


extern complex<double> fe_integraleLGP1xLGP1xNCP1_2D(double *p, double *q, complex<double> *r);
extern complex<double> fe_integraleLGP1xLGP1xNCP1_2D(complex<double> *p, double *q, double *r);
extern complex<double> fe_integraleLGP1xLGP1xNCP1_2D(complex<double> *p, double *q, complex<double> *r);
extern complex<double> fe_integraleLGP1xLGP1xNCP1_2D(double *p, int j, complex<double> *r);

extern void fe_multiple_integraleLGP1xLGP1xNCP1_2D (double *, double *, double*);

extern double fe_integraleLGP1xLGP1xNCP1xNCP1_2D (double *, double *, double *, double *);
extern double fe_integraleLGP1xLGP1xNCP1xNCP1_2D (double *, double *, double *, int);

extern double fe_integraleLGP1xNCP1xNCP1_2D (double *, double *, double *);
extern double fe_integraleLGP1xNCP1xNCP1_2D (double *, double *, int);

extern complex<double> fe_integraleLGP1xNCP1xNCP1_2D(double *p, complex<double> *q, complex<double> *r);

extern double fe_integraleLGP1xLGP1xLGP1xNCP1_2D(double *, double *, double *, double *);
extern double fe_integraleLGP1xLGP1xLGP1xNCP1_2D(double *, double *, int , double *);


extern  double fe_integraleLGP2_2D(double *);

extern double fe_integraleLGP2xLGP2_2D (double *, double *);
extern complex<double> fe_integraleLGP2xLGP2_2D (complex<double> *, complex<double> *);

extern double          fe_integraleLGP2xLGP2xLGP2_2D (double *, double *, double *);
extern complex<double> fe_integraleLGP2xLGP2xLGP2_2D (complex<double> *, complex<double> *, complex<double> *);
extern complex<double> fe_integraleLGP2xLGP2xLGP2_2D (double *, complex<double> *, complex<double> *);

extern complex<double>  fe_integraleLGP2xLGP2xLGP1_2D(double *p, complex<double>  *q, complex<double>  *r);
extern complex<double>  fe_integraleLGP2xLGP2xNCP1_2D(double *p, complex<double>  *q, complex<double>  *r);

extern double fe_integraleLGP2xLGP1_2D(double *, double *);
extern double fe_integraleLGP2xLGP1_2D(int,      double *);
extern double fe_integraleLGP2xLGP1_2D(double *, int);

extern complex<double> fe_integraleLGP2xLGP1_2D(complex<double> *p, int j);
extern complex<double> fe_integraleLGP2xLGP1_2D(int j, complex<double> *p);

extern double fe_integraleLGP2xLGP1xLGP1_2D(double *, double *, double*);
extern double fe_integraleLGP2xLGP1xLGP1_2D(double *, double *,int);
extern complex<double> fe_integraleLGP2xLGP1xLGP1_2D(double *p, complex<double> *q, complex<double> *r);

extern double fe_integrale_dxP2xLGP1_2D(double *, int);
extern double fe_integrale_dyP2xLGP1_2D(double *, int);

extern double fe_integrale_dxP2xLGP1xLGP1_2D(double *, double *, int);
extern double fe_integrale_dyP2xLGP1xLGP1_2D(double *, double *, int);

extern double fe_integrale_dxP2xLGP2xLGP1_2D(double *, double *, int);
extern double fe_integrale_dyP2xLGP2xLGP1_2D(double *, double *, int);

extern complex<double> fe_integrale_dxP2xLGP2xLGP1_2D(complex<double> *, complex<double> *, int);
extern complex<double> fe_integrale_dxP2xLGP2xLGP1_2D(complex<double> *, double *, int);
extern complex<double> fe_integrale_dxP2xLGP2xLGP1_2D(int , complex<double> *, complex<double> *);

extern complex<double> fe_integrale_dyP2xLGP2xLGP1_2D(complex<double> *, complex<double> *, int);
extern complex<double> fe_integrale_dyP2xLGP2xLGP1_2D(complex<double> *, double *, int);
extern complex<double> fe_integrale_dyP2xLGP2xLGP1_2D(int , complex<double> *, complex<double> *);

extern double fe_integrale_dxP2xLGP2xLGP1_2D(int, double *, int);
extern double fe_integrale_dyP2xLGP2xLGP1_2D(int, double *, int);

extern complex<double> fe_integrale_dxP2xLGP2xLGP2_2D(complex<double> *, double *, int);

extern complex<double> fe_integrale_dyP2xLGP2xLGP2_2D(complex<double> *, double *, int);

extern double fe_integraleLGP2xNCP1_2D(double *, double *);
extern double fe_integraleLGP2xNCP1_2D(double *, int);
extern double fe_integraleLGP2xNCP1_2D(int, double *);

extern complex<double> fe_integraleLGP2xNCP1_2D(complex<double> *p, int j);
extern complex<double> fe_integraleLGP2xNCP1_2D(int j, complex<double> *p);

extern double fe_integraleLGP2xNCP1xNCP1_2D(double *, double *, double *);
extern double fe_integraleLGP2xNCP1xNCP1_2D(double *, double *, int);

extern double fe_integrale_dxP2xNCP1_2D(double *,int );
extern double fe_integrale_dyP2xNCP1_2D(double *,int );

extern double fe_integrale_dxP2xLGP2xNCP1_2D(double *, double *, int);
extern double fe_integrale_dyP2xLGP2xNCP1_2D(double *, double *, int);
extern complex<double> fe_integrale_dxP2xLGP2xNCP1_2D(complex<double> *, double *, int);
extern complex<double> fe_integrale_dyP2xLGP2xNCP1_2D(complex<double> *, double *, int);

extern double fe_integrale_dxP2xLGP2xNCP1_2D(int, double *, int);
extern double fe_integrale_dyP2xLGP2xNCP1_2D(int, double *, int);

extern complex<double> fe_integrale_dxP2xLGP2xNCP1_2D(int i, complex <double> *q, complex <double> *r);
extern complex<double> fe_integrale_dyP2xLGP2xNCP1_2D(int i, complex <double> *q, complex <double> *r);

extern double fe_integrale_QP1xQP1_2D(double *p, double *q);
extern double fe_integrale_QP1xQP1xQP1_2D(int i, int j, double *r);

extern complex<double> fe_integrale_dxP2xLGP2xNCP1xNCP1_2D(complex<double> *, double *, double *, int);
extern complex<double> fe_integrale_dyP2xLGP2xNCP1xNCP1_2D(complex<double> *, double *, double *, int);

extern double fe_integrale_dxP2xLGP2xNCP1xNCP1_2D(int i, double *q, double *r, int l);
extern double fe_integrale_dyP2xLGP2xNCP1xNCP1_2D(int i, double *q, double *r, int l);

/**-----------------------------------------------------------------------------
Integrale 1D (segment)*/
extern double          fe_integraleLGP1xLGP1xLGP1_1D(double   *p, double   *q, double *r);
extern complex<double>        fe_integraleLGP1xLGP1xLGP1_1D(complex<double> *p, double   *q, double *r);
extern complex<double>        fe_integraleLGP1xLGP1xLGP1_1D(complex<double> *p, complex<double> *q, double *r);
extern double fe_integraleLGP1xLGP1xLGP1_1D (double *, double *, int);

extern complex<double> fe_integraleLGP1xLGP1xLGP1_1D(double *p, double *q, complex<double> *r);

extern double fe_integraleLGP1xLGP1xLGP1xLGP1_1D(double *p, double *q, double *r, double *s);
extern double fe_integraleLGP1xLGP1xLGP1xLGP1_1D(double *p, double *q, double *r, int);

extern double fe_integraleLGP1xLGP1xLGP1xLGP1xLGP1_1D(double *p, double *q, double *r, double *s, double *t);
extern double fe_integraleLGP1xLGP1xLGP1xLGP1xLGP1_1D(double *p, double *q, double *r, double *s, int);

extern double fe_integraleLGP1xLGP1xNCP1_1D (double *, double *, double *);
extern double fe_integraleLGP1xLGP1_1D(double *p, double *q);
extern complex<double> fe_integraleLGP1xLGP1_1D(complex<double> *p, double *q);
extern complex<double> fe_integraleLGP1xLGP1_1D(complex<double> *p, complex<double> *q);
extern double fe_integraleLGP1xLGP1_1D      (double *, int);

extern double fe_integraleLGP1_1D(double *p);
extern complex<double> fe_integraleLGP1_1D(complex<double> *p);

/**-----------------------------------------------------------------------------
Integrale 3D (prism) */
extern double fe_integraleLGP1xLGP1xLGP1_3D (double *, double *, double *);
extern double fe_integraleLGP1xLGP1xLGP1_3D (double *,int, int);

/**-----------------------------------------------------------------------------
Quadrature 2D */
//extern inline double fe_euclidianLGP1xLGP1_2D (double *, double *)  __attribute__((always_inline));
extern double fe_euclidianLGP1xLGP1_2D (double *, double *);
extern double fe_euclidianLGP1xLGP1_2D (double *, int);

extern double fe_euclidianLGP1xLGP1xLGP1_2D (double *, double *, double *);
extern double fe_euclidianLGP1xLGP1xLGP1_2D (double *, double *, int);

extern double fe_euclidianLGP1xLGP1xLGP1xLGP1_2D (double *, double *, double *, double *);
extern double fe_euclidianLGP1xLGP1xLGP1xLGP1_2D (double *, double *, double *, int);

extern double fe_euclidianLGP1xNCP1_2D (double *, double *);
//extern inline double fe_euclidianLGP1xNCP1_2D (double *, double *)  __attribute__((always_inline));

extern double fe_euclidianLGP1xLGP1xNCP1_2D (double *, double *, double *);
extern double fe_euclidianLGP1xLGP1xNCP1_2D (double *, int, double *);

extern void   integrale(mesh_t mesh, float *buffer, double *sum);
extern void fe_integrale_init(void);

extern double fe_integraleP0_1D(double *, double *, int);
extern double fe_integraleP1_1D(double *, double *, int);

/**-----------------------------------------------------------------------------
Line integrales 2D */
extern  double fe_LineIntegraleLGP1xLGP1_2D(double *p, int j, int k);
extern  double fe_LineIntegraleLGP1xLGP1_2D(double p[2][3], int j, int k, int l, int m);
extern  double fe_LineIntegraleLGP1xLGP1_2D(int j, double p[2][3], int k, int l, int m);
extern  double fe_LineAverageLGP1_2D(double *p, int i, double *q, int j);

extern  double fe_LineIntegraleNCP1xNCP1_2D(double p[2][3], int j, int k, int l, int m);
extern  double fe_LineIntegraleNCP1xNCP1_2D(int j, double M[2][3], int k, int l, int m);
extern  double fe_LineIntegraleNCP1xLGP1_2D(double p[2][3], int j, int k, int l, int m);
extern  double fe_LineIntegraleNCP1xLGP0_2D(double p[2][3], int j, int k, int l, int m);



/**-----------------------------------------------------------------------------
Gaussian quadrature 2D */

extern int gauss_init(int gauss_n, double *gauss_x, double *gauss_y,double *gauss_w);
extern int gauss_init_Q(int gauss_n, double *gauss_x, double *gauss_y,double *gauss_w);
//extern void gauss_init_R4(int gauss_n, float *gauss_x, float *gauss_y, float *gauss_w);


#endif /* FE_INTEGRATION_H */
