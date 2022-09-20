#include "tools-structures.h"
#include "tides.h"
#include "archive.h"

#ifndef TUGO_CLASSES_H
#define TUGO_CLASSES_H


typedef struct {
  double *zmin,*zmax;            /* min elevation, max elevation                                     */
  double *MKE,*MPE,*MBF;         /* mean kinetic energy, mean potential energy, mean bottom friction */
  double count;
} climatology2D_t;


class state_t {
private :
public :
  double **u, **v, **w;           /* velocities                           */
  double **U, **V, **W;           /* transport                            */
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check : MANDATORY !!!

  Note:

  26/10/2009

    U was massic velocities, it is now transport. It is now Us (for
    U star).

    Many routines used the former definition, it will be corrected 
    progressively

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  double **us, **vs, **ws;        /* massic velocities  (u star= u x rho) */
  double **Us, **Vs, **Ws;        /* massic transport   (U star= U x rho) */
  double **omega, **omicro;       /* vertical sigma velocities            */
  double **q2, **q2l;             /* tTKE                                 */
  double **T,  **S, **rho;        /* tracers, density                     */
  double **p, **q;                /* total and baroclinic pressure        */
  double **Ex,  **Ey;             /* Horizontal diffusion                 */
  double **Evx, **Evy;            /* Vertical diffusion                   */
  double **Ax, **Ay;              /* Horizontal advection                 */
  double **Kh3D;                  /* Horizontal diffusion coefficient     */
  double **px, **py;              /* total pressure gradient              */
  double **qx, **qy;              /* baroclinic pressure gradient         */
  double **z,**dsdt;              /* Immersion                            */
  double **z__;                   /*                                      */
  double **dsdx,**dsdy;           /* Upper/lower faces slope              */
  double **weight3D;              /* edge 3D integration weight           */
  double **volume;                /* prism volume                         */
  double *umean,*vmean;           /* depth average 3D velocities          */
  double *Exmean,*Eymean;         /* depth average 3D horizontal diffusion*/
  double *Axmean,*Aymean;         /* depth average 3D horizontal advection*/
  double **advection;             /* depth average 3D horizontal advection*/
  double *uexternal,*vexternal;   /* mean 2D velocities over a 3D step    */
  double *elevation,*r;           /* sea surface elevation, linear friction coefficient                */
  int    umode,wmode,tmode;
  int    udim;

  template <typename T> void u_erase(T ***v) {
    deletep2D(v,udim);
    }

  state_t() {
    u=0; v=0; w=0;
    U=0; V=0; W=0;
    omega=0; omicro=0;
    us=0; vs=0; ws=0;
    Us=0; Vs=0; Ws=0;
    q2=0; q2l=0;
    T=0; S=0; rho=0;
    p=0; q=0;
    Ex=0; Ey=0;
    Evx=0; Evy=0;
    Ax=0; Ay=0;
    Kh3D=0;
    px=0; py=0;
    qx=0; qy=0;
    z=0; dsdt=0;
    z__=0;
    dsdx=0; dsdy=0;
    weight3D=0;
    }
 
  void destroy() {
    u_erase(&u);
    u_erase(&v);
//     erase(&w);
//     erase(&us);erase(&vs);erase(&ws);
//     erase(&U); erase(&V); erase(&W);
//     erase(&Us); erase(&Vs); erase(&Ws);
//     erase(&omega); erase(&omicro);
    u_erase(&q2);
    u_erase(&q2l);
//     erase(&T); erase(&S); erase(&rho);
//     erase(&p);erase(&q);
//     erase(&Ex); erase(&Ey);
//     erase(&Evx);erase(& Evy);
//     erase(&Ax); erase(&Ay);
//     erase(&Kh3D);
//     erase(&px);erase(& py);
//     erase(&qx);erase(&qy);
//     erase(&z); erase(&dsdt);
//     erase(&z__);
//     erase(&dsdx); erase(&dsdy);
//     erase(&weight3D);
    }
//     double * &operator [].u (size_t i) {
//       return u[i];
//       }
};

#include "specification.h"

class wave_t {
  private:
  public:
    double *Sxx, *Sxy, *Syy; /* Radiation stress tensor components             */
    double *RSdvg_x, *RSdvg_y;    /* Radiation stress divergence components (wave forcing vector)              */

    wave_t() {
      Sxx = 0;
      Sxy = 0;
      Syy = 0;
      RSdvg_x = 0;
      RSdvg_y = 0;
    }

    int allocate(int dim) {
      Sxx  = new double[dim];
      if(Sxx==0) return(-1);
      Sxy  = new double[dim];
      if(Sxy==0) return(-1);
      Syy  = new double[dim];
      if(Syy==0) return(-1);
      RSdvg_x  = new double[dim];
      if(RSdvg_x==0) return(-1);
      RSdvg_y  = new double[dim];
      if(RSdvg_y==0) return(-1);
      return(0);
    }

    void destroy() {
      deletep(&Sxx);
      deletep(&Sxy);
      deletep(&Syy);
      deletep(&RSdvg_x);
      deletep(&RSdvg_y);
    }
};

#include "polygones.h"

#include "poc-netcdf.hpp"


class hrhs_t {
private :
public :
  double **z,**u,**v;       /*               */
  int zdim,udim,neq;

  hrhs_t() {
     z=u=v=0;
     }

  void initialise() {
     int k,n;
     for(n=0;n<zdim;n++) {
       for(k=0;k<neq;k++) {
         z[n][k]=0.0;
         }
       }
     for(n=0;n<udim;n++) {
       for(k=0;k<neq;k++) {
         u[n][k]=0.0;
         v[n][k]=0.0;
         }
       }
     }

  void destroy() {
     int n;
     if(z!=0) {
       for(n=0;n<zdim;n++) {
         delete[] z[n];
         }
       }
     delete[] z;
     z=0;
     if(u!=0) {
       for(n=0;n<udim;n++) {
         delete[] u[n];
         }
       }
     delete[] u;
     u=0;
     if(v!=0) {
       for(n=0;n<udim;n++) {
         delete[] v[n];
         }
       }
     delete[] v;
     v=0;
     }
};

#include "init_config.hpp"

#endif /* TUGO_CLASSES_H */
