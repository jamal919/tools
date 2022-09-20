
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Yoann Le Bars      LEGOS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief hconstant_t and harmonic_t clases
*/
/*----------------------------------------------------------------------------*/

#if HARMONIC_H == 0
#define HARMONIC_H 1

#include "tools-structures.h"
#include "statistic.h"

class hconstant_t {
private :
public :
  float *a,*G; /* amplitude and phase */
  complex<float> *z;
  size_t size; /* number of waves */
  spectrum_t *s;

  hconstant_t() {
    a=G=NULL;
    z=NULL;
    size=0;
    s=NULL;
    }
  
  hconstant_t(mgr_t mgr) {
    a=G=NULL;
    z=NULL;
    size=0;
    s=NULL;
    init_polar(mgr.nwave);
    for(int k=0;k<size;k++) {
      a[k]=mgr.data[k].amp;
      G[k]=mgr.data[k].phi;
      }
    }

  hconstant_t(int nwave) {
    a=G=NULL;
    z=NULL;
    size=0;
    s=NULL;
    init_polar(nwave);
    for(int k=0;k<size;k++) {
      a[k]=0.0;
      G[k]=0.0;
      }
    }

  void init_polar(size_t newSize) {
    deletep(&a);
    deletep(&G);
    size=newSize;
    a=new float[size];
    G=new float[size];
    }

  void init_complex(size_t newSize) {
    deletep(&z);
    size=newSize;
    z=new complex<float>[size];
    }

  void allocate(size_t newSize) {
    init_polar(newSize);
    init_complex(newSize);
    }

  void set_complex(bool cleanUp=false) {
    if(z!=NULL) delete[] z;
    z=new complex<float>[size];
    for(size_t i=0;i<size;i++) {
      z[i]=polar<float>(a[i],-G[i]*d2r);
      };
    
    if(cleanUp==false)
      return;
    deletep(&a);
    deletep(&G);
    }

  void set_polar(int degrees=0) {
    if(a!=NULL) delete[] a;
    if(G!=NULL) delete[] G;
    a=new float[size];
    G=new float[size];
    for(size_t i=0;i<size;i++) {
      a[i]=abs(z[i]);
      if(degrees)
        G[i]=-arg(z[i])*r2d;
      else
        G[i]=-arg(z[i]);
      };
    }

  void destroy() {
    deletep(&a);
    deletep(&G);
    deletep(&z);
    size=0;
    }
    
  bool hasMaskValue(float mask){
    for(size_t i = 0; i < size; ++i){
      if( (is_equal(a[i], mask, float(1.e-3)))
        || (is_equal(G[i], mask, float(1.e-3))) ){
        return(true);
      }
    }
    return(false);
  }
  
  bool isEmpty(){
    return(size == 0);
  };
    
};


#include <admittance-mts.h>

class harmonic_t {
private :
public :
  double *A, *M;
  double ***rhs;        ///< right-hand side (RHS) vector : [nrhs][nndes[]][neq]
  double **cs, **sn;
  int *pivot;
  int nframes;          ///< number of frames
  int neq;
  int nrhs;             ///< number of variables
  size_t *nndes;           ///< number of points in the grid : [nrhs]
  spectrum_t spectrum;
  bool **mask;          ///< mask : [nrhs][nndes[]] false is out of domain
  
  harmonic_t() {
    A=M=0;
    cs=sn=0;
    pivot=0;
//    n=0;
    neq=0;
    rhs=0;
    nndes=0;
    mask=0;
    }
  
  int allocate_coefficients() {
    int n,nw;
    
    nw=spectrum.n;
    if(nw<=0 || nframes<=0) return -1;
    
    neq=2*nw;
    
    cs=new double*[nw];
    sn=new double*[nw];
    for(n = 0; n < nw; n++) {
      cs[n]=new double[nframes];
      sn[n]=new double[nframes];
      }
    return(0);
    }
  
  int allocate_matrix() {
    int nw,neq2;
    
    nw=spectrum.n;
    if(nw<=0) return -1;
    
    neq=2*nw;
    
    neq2=neq*neq;
    
    A = new double[neq2];
    for(int k = 0; k < neq2; k++) {
      A[k] = 0;
      }
    return(0);
    }
  
  int unity(int row) {
    int nw;
    
    nw=spectrum.n;
    if(nw<=0 || A==0) return -1;
    
    neq=2*nw;
    
    for(int k = 0; k < neq; k++) {
      A[k*neq+row] = 0;
      }
    
    A[row*neq+row] = 1.0;
    
    return(0);
    }
  
  int substitute(int row, double *c) {
    int nw;
    
    nw=spectrum.n;
    if(nw<=0 || A==0) return -1;
    
    neq=2*nw;
    
    for(int k = 0; k < neq; k++) {
      A[k*neq+row] = c[k];
      }
    
    return(0);
    }
  
  int factorize(int destructive);
  
  int solve(double *b, int nrhs, int destructive);
  
  void set_rhs(double *serie, double* & b, int nrhs, int destructive) {
    deletep(&b);
    b=new double[neq];
    for(size_t k = 0; k < neq; k++) b[k]=0.0;
    for(size_t k = 0; k < spectrum.n; k++) {
      for(size_t m = 0; m < nframes; m++) {
        b[2 * k]     += cs[k][m] * serie[m];
        b[2 * k + 1] += sn[k][m] * serie[m];
        }
      }
    }
  
  void solvable(int *keep) {
    double *tmp=new double[nframes], mask=nan("?");
    statistic_t s,s1,s2;
    
    for(size_t i = 0; i < spectrum.n; i++){
      for(size_t m = 0; m < nframes; m++){
        tmp[m]=cs[i][m];
        }
      s=get_statistics(tmp,mask,nframes,0);
      if(s.std<0.1) {
        printf("component %d %s : std = %lf\n",i, spectrum.waves[i].name,s.std);
        keep[i]=0;
        continue;
        }
      vector<double> sub1,sub2;
      for(size_t m = 0; m < nframes; m++){
        if(tmp[m]>s.mean) sub1.push_back(tmp[m]);
        else sub2.push_back(tmp[m]);
        }
      s1=get_statistics(sub1, mask, 0);
      s2=get_statistics(sub2, mask, 0);
      if((s1.std<0.1) && (s2.std<0.1)) {
        printf("component %d %s : std = %lf %lf\n",i, spectrum.waves[i].name,s1.std,s2.std);
        keep[i]=0;
        continue;
        }
      }
    }
  
  int error(double *residuals, int n, double* & error);

  void set_MatrixAdmittance(tidal_wave & wave, int k, admittance_t & admittance) {
    double c1[neq], c2[neq];
    for(size_t j = 0; j < neq; j++) {
      c1[j] = 0.;
      c2[j] = 0.;
      }
    c1[2 * k]     = 1. / wave.Ap;
    c2[2 * k + 1] = 1. / wave.Ap;

    size_t species = wave.nT;
    double coef[3] = {1., 1., 1.};
    switch(admittance.method[species]){
          
      case ADMITTANCE_UNAVAILABLE:
//         return(0);
        break;

      case ADMITTANCE_SPLINE:
        admittance_mts_sweightP(admittance, wave, coef);
        for(size_t l = 0; l < 3; l++) {
          size_t j = admittance.windex_s[species][l];
          if(j==-1)TRAP_ERR_EXIT(ENOEXEC,"admittance waves index for nT=%d not initialized\n",k);
          c1[2 * j]     = -coef[l];
          c2[2 * j + 1] = -coef[l];
          }
        break;

      case ADMITTANCE_LINEAR:
        admittance_mts_lweightP(admittance, wave, coef);
        for(size_t l = 0; l < 2; l++) {
          size_t j = admittance.windex_l[species][l];
          if(j==-1)TRAP_ERR_EXIT(ENOEXEC,"admittance waves index for nT=%d not initialized\n",k);
          c1[2 * j]     = -coef[l];
          c2[2 * j + 1] = -coef[l];
          }
        break;

      default:
        check_error(-1, "admittance method is not known", __LINE__, __FILE__, 1);
        break;
      }
    this->substitute(2*k,  c1);
    this->substitute(2*k+1,c2);
    }
  
  void set_RhsAdmittance(double *b, double *zr, double *zi, tidal_wave & wave, int k, admittance_t & admittance) {
//     double c1[neq], c2[neq];
//     for(size_t j = 0; j < neq; j++) {
//       c1[j] = 0.;
//       c2[j] = 0.;
//       }
//     c1[2 * k]     = zr[k] / wave.Ap;
//     c2[2 * k + 1] = zi[k] / wave.Ap;
//
//     size_t species = wave.nT;
//     double coef[3] = {1., 1., 1.};
//     switch(admittance.method[species]){
//
//       case ADMITTANCE_UNAVAILABLE:
// //         return(0);
//         break;
//
//       case ADMITTANCE_SPLINE:
//         admittance_mts_sweightP(admittance, wave, coef);
//         for(size_t l = 0; l < 3; l++) {
//           size_t j = admittance.windex_s[species][l];
//           if(j==-1)TRAP_ERR_EXIT(ENOEXEC,"admittance waves index for nT=%d not initialized\n",k);
//           c1[2 * j]     = -zr[j]*coef[l];
//           c2[2 * j + 1] = -zi[j]*coef[l];
//           }
//         break;
//
//       case ADMITTANCE_LINEAR:
//         admittance_mts_lweightP(admittance, wave, coef);
//         for(size_t l = 0; l < 2; l++) {
//           size_t j = admittance.windex_l[species][l];
//           if(j==-1)TRAP_ERR_EXIT(ENOEXEC,"admittance waves index for nT=%d not initialized\n",k);
//           c1[2 * j]     = -zr[j]*coef[l];
//           c2[2 * j + 1] = -zi[j]*coef[l];
//           }
//         break;
//
//       default:
//         check_error(-1, "admittance method is not known", __LINE__, __FILE__, 1);
//         break;
//       }
    b[2*k]  =0.0;
    b[2*k+1]=0.0;
    for(size_t j = 0; j < neq; j++) {
      b[2*k]  -=A[j*neq +2*k]  *zr[j/2];
      b[2*k+1]-=A[j*neq +2*k+1]*zi[j/2];
      }
    }
  
  void destroy() {
    const int n=neq/2;
    
    if(M!=A) delete[] M;
    deletep(&A);
    M=0;
    
    deletep(&pivot);
    
    deletep2D(&cs,n);
    deletep2D(&sn,n);
    
    neq=0;
    }
  
  };

#endif
