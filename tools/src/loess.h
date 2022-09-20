

#if LOESS_H == 0
#define LOESS_H 1

#include "tools-structures.h"
#include "map.h"
extern void Loess1D_init(double dx, double l_c, int m, double *weight, int *);
extern void Loess1D_BF(int n,const double *h, double mask, int m,const double *weight, double *hf);

extern void Loess1D_BF(float *, float, int, float *, int, float *);

extern void Loess1D_BF_nomask(const float *buf, int nval, const float *weight, int nweight, float *out);
extern void Loess1D_BF_nomask(const double *buf, int nval,const double *weight, int nweight, double *out);

extern void Loess1D(int n, double dx, double l_c,const double *h, double mask,double *lf,int *rstatus);
extern void Loess1D_init(float step, int nval, float scale, float **weight, int *);

extern int loess1d_irregular(int n, double l_c, double *h, double *x, double mask,double *lf);
extern int loess1d_irregular(int n, float l_c, float *h, float *x, float mask,float *lf);

extern void loess1d(float *buf, float mask, float step, int nval, float scale, float *out);
extern int  loess1d(int n, double l_c, complex<float> *h, double *x, complex<float> mask,complex<float> *lf);

extern void  Loess2D(grid_t grid, float *buf, float mask, float scale, float *out);
extern float *Loess2D_BF(grid_t grid, float *buf, float mask, float scale, float *weight, float *out);
extern void  Loess2D_init(grid_t grid, float scale, float **weight);

extern int Loess2D_init_cartesian(double scale_x, double scale_y, double azimuth, double dx, double dy, int nxmax, int nymax, filter_t<float> & filter);

#endif
