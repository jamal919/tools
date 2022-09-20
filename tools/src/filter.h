
#if FILTER_H == 0
#define FILTER_H 1

#include "tools-structures.h"

#define LOESS     0
#define LANCZOS   1
#define MEDIAN    2

template <typename T> class filter_t
{
private:
public:
  int nx, ny;
  T   *weights;
  double scale_x, scale_y, azimuth;
  
  filter_t() {
    nx=0;
    ny=0;
    weights=0;
    scale_x=scale_y=azimuth=0;
    }
    
  void clear() {
    nx=0;
    ny=0;
    if(weights!=0) delete[] weights;
    weights=0;
    scale_x=scale_y=azimuth=0;
    }
};

// functions de fourier.cpp
#include "fourier.h"

// functions de loess.cpp
#include "loess.h"

// functions de lanczos.cpp
#include "lanczos.h"

extern int filter1D_BF(int mode, double  L, double *x, int nvalues, complex<float> *buffer, complex<float> mask, complex<float>* & BF);
extern int filter1D_BF(int mode, double *L, double *x, int nvalues, complex<float> *buffer, complex<float> mask, complex<float>* & BF);

extern int filter2D_InitWeigths(int mode, double L, double dx, double dy, filter_t<float> & filter);
extern int filter2D_InitWeigths(int mode, double scale_x, double scale_y, double dx, double dy, int nxmax, int nymax, filter_t<float> & filter);

extern statistic_t filtering(const char *save,const double *h,const double *t, double mask, int n, double dt);
extern void season(char *save, double *h, double *t, double mask, int n, double dt);
extern void fourier(const char *save, double *h, double mask, int n, double dt, int option);

extern double *tide_filters(const char *filtername, double *h, double mask, int nvalues );
extern double *tide_filters(const char *filtername, double, double *h, double mask, int nvalues );

extern int HarmonicFilter_init(double *t, const vector<double> & f, double* &H, double* &LF, size_t ndata);
#define HarmonicFilterSpectrumType 4
extern int HarmonicFilter(double *t, double *buffer, double *filtered, double mask, size_t ndata);
extern int HarmonicFilter_experiment();

#endif
