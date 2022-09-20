
#include "filter.h" /* for filter_t */

#if LANCZOS_H == 0
#define LANCZOS_H 1


extern void lanczos1D_init(double f_c, int l_f, double *weight);

extern void lanczos1D_init(float step, int nval, float scale, float **weight, int *nweight);
extern void lanczos1D_init(double step, int nval, double scale, double **weight, int *nweight);

extern void lanczos1D_BF(int n, double *h, double mask, double *weight, int m, double *hf);
extern int  lanczos1D(int n, double dx, double l_c, double *h, double mask,double *lf);
extern int  lanczos1D(int n, double dx, double *l_c, double *h, double mask,double *lf);

extern int Lanczos1D_BF(double  L, int a, double *x, int nvalues, complex<float> *buffer, complex<float> mask, complex<float>* & BF);
extern int Lanczos1D_BF(double *L, int a, double *x, int nvalues, complex<float> *buffer, complex<float> mask, complex<float>* & BF);

extern int Lanczos2D_weights(double L, double dx, double dy, int a, filter_t<float> & filter);

#endif
