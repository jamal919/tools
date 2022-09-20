
/******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief map operation (mainly gradients and diffusions) definitions
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>

#include "tools-structures.h"

#include "maths.h"
#include "geo.h"
#include "constants.h"
#include "map.def"
#include "map.h"
#include "filter.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int map_gradient_00(const grid_t &grid,int n, T *dum, T mask, int ref, T *dumx, T *dumy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  returns gradient in input units per metres
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
#define RADIUS 6300.e+03
#define d2r M_PI/180.0
  int  status=0;
  size_t lmin = 0;
  size_t lmax = grid.ny;
  size_t kmin = 0;
  size_t kmax = grid.nx;

  int nRequestedProcs=-1;
  int nprocs=initialize_OPENMP(nRequestedProcs);
  
#pragma omp parallel for if(nprocs>1)
  for(size_t l = 0; l < grid.ny; l++) {
    for(size_t k = 0; k < grid.nx; k++) {
      size_t m = l * grid.nx + k;
      dumx[m] = mask;
      dumy[m] = mask;
      }
    }

#pragma omp parallel for if(nprocs>1)
  for(int l = lmin; l < lmax; l++) {
//     size_t lm1=max(0,l-1);
//     size_t lp1=min(grid.ny-1,l+1);
    int lm1=max(0,l-1);
    int lp1=min(grid.ny-1,l+1);
    double dy = grid.dy*(lp1-lm1);
    for(int k = kmin; k < kmax; k++) {
//       size_t km1=max(0,k-1);
//       size_t kp1=min(grid.nx-1,k+1);
      int km1=max(0,k-1);
      int kp1=min(grid.nx-1,k+1);
      double dx = grid.dx*(kp1-km1);
/* *----------------------------------------------------------------------------
      gradient along x*/
//       size_t m = l * grid.nx + k;
//       size_t m1 = l * grid.nx + km1;
//       size_t m2 = l * grid.nx + kp1;
      int m = l * grid.nx + k;
      int m1 = l * grid.nx + km1;
      int m2 = l * grid.nx + kp1;
      if( (dum[m1]!=mask)  && (dum[m2]!=mask) ) {
        T dz = dum[m2] - dum[m1];
        dumx[m] = dz / (T) dx;
        }
/* *----------------------------------------------------------------------------
      gradient along y*/
      m1 = lm1 * grid.nx + k;
      m2 = lp1 * grid.nx + k;
      if( (dum[m1]!=mask)  && (dum[m2]!=mask) ) {
       T dz = dum[m2] - dum[m1];
        dumy[m] = dz / (T) dy;
        }
      }
    }

  if(ref == GEOCENTRIC) {
    const double f=d2r * RADIUS;
#pragma omp parallel for if(nprocs>1)
    for(size_t l = lmin; l < lmax; l++) {
      double yy=map_grid_y (grid, 0, l);
      double cosine=cos(yy * d2r);
      for(size_t k = kmin; k < kmax; k++) {
        size_t m = l * grid.nx + k;
        if(dumx[m]!=mask) dumx[m] /= f * cosine;
        if(dumy[m]!=mask) dumy[m] /= f;
        }
      }
    }

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_gradient_02(const grid_t &grid,int n, T *dum, T mask, int ref, T *dumx, T *dumy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  returns gradient in input units per metres
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
#define RADIUS 6300.e+03
#define d2r M_PI/180.0
  int status=0;
  int lmin, lmax, kmin, kmax;

  lmin = 0;
  lmax = grid.ny;
  kmin = 0;
  kmax = grid.nx;

#pragma omp parallel for
  for(int l = 0; l < grid.ny; l++) {
    for(int k = 0; k < grid.nx; k++) {
      int64_t m = grid.Hindex(k,l);
      dumx[m] = mask;
      dumy[m] = mask;
      }
    }

#pragma omp parallel for
  for(int l = lmin; l < lmax; l++) {
    int km1,kp1,lm1,lp1;
    int64_t m, m1, m2;
    double xx, yy;
    double xx1, yy1;
    double xx2, yy2;
    double dx, dy;
    T dz;
    lm1=max(0,l-1);
    lp1=min(grid.ny-1,l+1);
    for(int k = kmin; k < kmax; k++) {
      km1=max(0,k-1);
      kp1=min(grid.nx-1,k+1);
/* *----------------------------------------------------------------------------
      gradient along x*/
      m = grid.Hindex(k,l);
      if(dum[m]==mask) continue;
      xx=map_grid_x (grid, k, l);
      yy=map_grid_y (grid, k, l);
      m1 = grid.Hindex(km1,l);
      m2 = grid.Hindex(kp1,l);
      if(dum[m1]==mask) continue;
      if(dum[m2]==mask) continue;
      xx1 = map_grid_x (grid, km1, l);
      yy1 = map_grid_y (grid, km1, l);
      xx2 = map_grid_x (grid, kp1, l);
      yy2 = map_grid_y (grid, kp1, l);
      dx = xx2 - xx1;
      dy = yy2 - yy1;
      dz = dum[m2] - dum[m1];
      dumx[m] = dz / (T) dx;
/* *----------------------------------------------------------------------------
      gradient along y*/
      if(lmin==lmax-1) {
        dumy[m] = 0.0;
        }
      else {
        m1 = grid.Hindex(k,lm1);
        m2 = grid.Hindex(k,lp1);
        if(dum[m1]==mask) continue;
        if(dum[m2]==mask) continue;
        xx1 = map_grid_x (grid, k, lm1);
        yy1 = map_grid_y (grid, k, lm1);
        xx2 = map_grid_x (grid, k, lp1);
        yy2 = map_grid_y (grid, k, lp1);
        dx = xx2 - xx1;
        dy = yy2 - yy1;
        dz = dum[m2] - dum[m1];
        dumy[m] = dz / (T) dy;
        }
      if(ref == GEOCENTRIC) {
        dumx[m] /= d2r * RADIUS * cos(yy * d2r);
        dumy[m] /= d2r * RADIUS;
        }
      }
    }

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int map_gradient_template(const grid_t &grid,int n, T *dum, T mask, int ref, T *dumx, T *dumy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  switch (grid.modeH) {
    case 0:
      status = map_gradient_00(grid, n, dum, mask, ref, dumx, dumy);
      break;

    case 1:
      status = map_gradient_02(grid, n, dum, mask, ref, dumx, dumy);
      break;

    case 2:
      status = map_gradient_02(grid, n, dum, mask, ref, dumx, dumy);
      break;
    }
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_gradient(const grid_t &grid,int n, float *dum, float mask, int ref, float *dumx, float *dumy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status = map_gradient_template(grid, n, dum, mask, ref, dumx, dumy);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_gradient(const grid_t &grid,int n, complex<float> *dum, complex<float> mask, int ref, complex<float> *dumx, complex<float> *dumy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status = map_gradient_template(grid, n, dum, mask, ref, dumx, dumy);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_divergence_02(const grid_t &grid, T *dumx, T *dumy, T mask, int ref, T *divergence)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

  returns gradient in input units per metres
  
-----------------------------------------------------------------------------*/
{
#define RADIUS 6300.e+03
#define d2r M_PI/180.0
  int k, l, m, status=0;
  int m1, m2;
  int lmin, lmax, kmin, kmax;
  int km1,kp1,lm1,lp1;
  double xx, yy;
  double xx1, yy1;
  double xx2, yy2;
  double dx, dy;
  T dz;

  lmin = 0;
  lmax = grid.ny;
  kmin = 0;
  kmax = grid.nx;

  for(l = 0; l < grid.ny; l++) {
    for(k = 0; k < grid.nx; k++) {
      m = l * grid.nx + k;
      divergence[m] = mask;
      }
    }

  for(l = lmin; l < lmax; l++) {
    lm1=max(0,l-1);
    lp1=min(grid.ny-1,l+1);
    for(k = kmin; k < kmax; k++) {
      double dvx,dvy;
      km1=max(0,k-1);
      kp1=min(grid.nx-1,k+1);
/* *----------------------------------------------------------------------------
      gradient along x*/
      m = l * grid.nx + k;
      if(dumx[m]==mask) continue;
      if(dumy[m]==mask) continue;
      xx=map_grid_x (grid, k, l);
      yy=map_grid_y (grid, k, l);
      m1 = l * grid.nx + km1;
      m2 = l * grid.nx + kp1;
      if(dumx[m1]==mask) continue;
      if(dumx[m2]==mask) continue;
      xx1 = map_grid_x (grid, km1, l);
      yy1 = map_grid_y (grid, km1, l);
      xx2 = map_grid_x (grid, kp1, l);
      yy2 = map_grid_y (grid, kp1, l);
      dx = xx2 - xx1;
      dy = yy2 - yy1;
      dz = dumx[m2] - dumx[m1];
      dvx = dz / (T) dx;
/* *----------------------------------------------------------------------------
      gradient along y*/
      if(lmin==lmax-1) {
        dvy = 0.0;
        }
      else {
        m1 = lm1 * grid.nx + k;
        m2 = lp1 * grid.nx + k;
        if(dumy[m1]==mask) continue;
        if(dumy[m2]==mask) continue;
        xx1 = map_grid_x (grid, k, lm1);
        yy1 = map_grid_y (grid, k, lm1);
        xx2 = map_grid_x (grid, k, lp1);
        yy2 = map_grid_y (grid, k, lp1);
        dx = xx2 - xx1;
        dy = yy2 - yy1;
        dz = dumy[m2] - dumy[m1];
        dvy = dz / (T) dy;
        }
      if(ref == GEOCENTRIC) {
        dvx /= d2r * RADIUS * cos(yy * d2r);
        dvy /= d2r * RADIUS;
        }
      divergence[m]=dvx+dvy;
      }
    }

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int map_divergence_template(const grid_t &grid, T *dumx, T *dumy, T mask, int ref, T *divergence)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  switch (grid.modeH) {
    case 0:
//       status = map_gradient_00(grid, n, dum, mask, ref, dumx, dumy);
      break;

    case 1:
      status = map_divergence_02(grid, dumx, dumy, mask, ref, divergence);
      break;
    case 2:
      status = map_divergence_02(grid, dumx, dumy, mask, ref, divergence);
      break;
  }
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_divergence(const grid_t &grid, double *dumx, double *dumy, double mask, int ref, double *divergence)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=map_divergence_template(grid, dumx, dumy, mask, ref, divergence);
  
  return (status);
}




/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_curve00(const grid_t &grid, T *dum, T mask, int ref, T *dumx, T *dumy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*returns second derivative in input units per square metres */
{
  /* not finished properly, very approximative */
  int k,l,m,status=0;
  int m1,m2;
  int lmin,lmax,kmin,kmax;
  double yy;
  double dx,dy;
  T dz;
  double c;

  lmin=1;
  lmax=grid.ny-1;
  kmin=1;
  kmax=grid.nx-1;

  for( l=0; l<grid.ny; l++) {
    for( k=0; k <grid.nx; k++) {
      m=l*grid.nx+k;
      dumx[m]=mask;
      dumy[m]=mask;
      }
    }

  for( l=lmin; l<lmax; l++) {
    for( k=kmin; k <kmax; k++) {
      m=l*grid.nx+k;
      yy=grid.y[m];
      m1=l*grid.nx+k-1;
      m2=l*grid.nx+k+1;
      dx=grid.dx;
      dz=(T) (dum[m2]-2.*dum[m]+dum[m1])/dx;
      dumx[m]=dz/dx;
      m1=(l-1)*grid.nx+k;
      m2=(l+1)*grid.nx+k;
      dy=grid.dy;
      dz=(T) (dum[m2]-2.*dum[m]+dum[m1])/dy;
      dumy[m]=dz/dy;
      switch (ref) {
        case GEOCENTRIC:
          c=cos(yy*d2r);
          dumx[m]/=c*c;
          break;

//         case GEOMETRIC:
        case CARTESIAN:
          c=cos(yy*d2r);
          dumx[m]/=MeanEarthRadius2*c*c;
          dumy[m]/=MeanEarthRadius2;
          break;
        }
      }
    }


  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_curve02(const grid_t &grid, T *dum, T mask, int ref, T *dumx, T *dumy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*returns second derivative in input units per square metres */
/* not finished properly, very approximative */
  int k,l,m,status=0;
  int m1,m2;
  int lmin,lmax,kmin,kmax;
  double xx,yy;
  double xx1,yy1;
  double xx2,yy2;
  double dx,dy;
  T dz;
  double c;

  lmin=1;
  lmax=grid.ny-1;
  kmin=1;
  kmax=grid.nx-1;

  for( l=0; l<grid.ny; l++) {
    for( k=0; k <grid.nx; k++) {
      m=l*grid.nx+k;
      dumx[m]=mask;
      dumy[m]=mask;
      }
    }

  for( l=lmin; l<lmax; l++) {
    for( k=kmin; k <kmax; k++) {
      m=l*grid.nx+k;
      if(dum[m]==mask) continue;
      xx=grid.x[m];
      yy=grid.y[m];
      m1=l*grid.nx+k-1;
      if(dum[m1]==mask) continue;
      m2=l*grid.nx+k+1;
      if(dum[m2]==mask) continue;
      xx1=grid.x[m1];
      xx2=grid.x[m2];
      dx=0.5*(xx2-xx1);
      dz=(T) (dum[m2]-dum[m])/(xx2-xx)-(dum[m]-dum[m1])/(xx-xx1);
      dumx[m]=dz/dx;
      m1=(l-1)*grid.nx+k;
      if(dum[m1]==mask) continue;
      m2=(l+1)*grid.nx+k;
      if(dum[m2]==mask) continue;
      yy1=grid.y[m1];
      yy2=grid.y[m2];
      dy=0.5*(yy2-yy1);
      dz=(T) (dum[m2]-dum[m])/(yy2-yy)-(dum[m]-dum[m1])/(yy-yy1);
      dumy[m]=dz/dy;
//       if(ref==GEOCENTRIC) {
//         c=cos(yy*d2r);
//         dumx[m]/=MeanEarthRadius2*c*c;
//         dumy[m]/=MeanEarthRadius2;
//         }
      }
    }
 
  size_t count;
  
  count=occurence(mask, dumx, grid.Hsize());
  count=occurence(mask, dumy, grid.Hsize());
  
  double factor=180.*180./M_PI/M_PI/MeanEarthRadius2;
  
  if(ref==GEOCENTRIC) {
    for( l=lmin; l<lmax; l++) {
      for( k=kmin; k <kmax; k++) {
        m=l*grid.nx+k;
        yy=grid.y[m];
        c=cos(yy*d2r);
        if(dumx[m]!=mask) dumx[m]*=factor/c/c;
        if(dumy[m]!=mask) dumy[m]*=factor;
        }
      }
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_curve_template(const grid_t &grid, T *dum, T mask, int ref, T *dumx, T *dumy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  switch(grid.modeH)
    {
    case 0:
      status=map_curve00(grid, dum, mask, ref, dumx, dumy);
      break;

    case 2:
      status=map_curve02(grid, dum, mask, ref, dumx, dumy);
      break;
    }

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_curve(const grid_t &grid, float *dum, float mask, int ref, float *dumx, float *dumy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=map_curve_template(grid, dum, mask, ref, dumx, dumy);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_curve(const grid_t &grid, double *dum, double mask, int ref, double *dumx, double *dumy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=map_curve_template(grid, dum, mask, ref, dumx, dumy);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_curve(const grid_t &grid, complex<double> *dum, complex<double> mask, int ref, complex<double> *dumx, complex<double> *dumy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=map_curve_template(grid, dum, mask, ref, dumx, dumy);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_laplacian_template(const grid_t &grid, T *dum, T mask, int ref, T* & laplacian)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  T *dumx, *dumy;
  
  dumx=new T[grid.Hsize()];
  dumy=new T[grid.Hsize()];
  
  if(laplacian==0) laplacian=new T[grid.Hsize()];
  
  status=map_curve(grid, dum, mask, ref, dumx, dumy);
  
  for(size_t m=0; m<grid.Hsize();m++) {
    if(dumx[m]!=mask and dumy[m]!=mask) laplacian[m]=dumx[m]+dumy[m];
//     if(dumx[m]!=mask and dumy[m]!=mask) laplacian[m]=dumx[m];
    else laplacian[m]=mask;
    }
    
  delete[] dumx;
  delete[] dumy;
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_laplacian(const grid_t &grid, complex<double> *dum, complex<double> mask, int ref, complex<double> * & laplacian)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=map_laplacian_template(grid, dum, mask, ref, laplacian);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename TYPE> int map_smooth_local_template(const grid_t &grid, TYPE *buf, TYPE mask, float scale, int k, int l, TYPE *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j;
  int imin,imax,jmin,jmax;
  float w,r,s,xx,yy;
  float lx,ly,q,dum;
  int nx,ny,offset_j,offset_l;
  float *weight,sum;
  TYPE value,z,tmp;

  nx=int(NINT(scale/grid.dx+1)+1);
  ny=int(NINT(scale/grid.dy+1)+1);

  nx=min(nx,grid.nx);
  ny=min(ny,grid.ny);

  weight=new float[nx*ny];

  for(j=0;j<ny;j++)
    for (i=0;i<nx;i++)
      weight[j*nx+i]=0.;

  for(j=0;j<ny;j++)
    for(i=0;i<nx;i++) {
      lx=i*grid.dx/scale;
      ly=j*grid.dy/scale;
      q=lx*lx+ly*ly;
      if (q <= 1.) {
        dum=1-q*q*q;
        weight[j*nx+i]=dum*dum*dum;
        }
      }

  tmp=0.;
  sum=0.;

  imin=max(0,k-nx+1);
  imax=min(grid.nx,k+nx);
  jmin=max(0,l-ny+1);
  jmax=min(grid.ny,l+ny);

  for(j=jmin;j<jmax;j++) {
    offset_j=j*grid.nx;
    for(i=imin;i<imax;i++) {
      z=buf[offset_j+i];
      if(z!=mask) {
        w=weight[abs(j-l)*nx+abs(i-k)];
        sum+=w;
        tmp+=z*w;
        }
      }
    }

  if(sum!=0.) {
    value=tmp/sum;
    }
  else {
    value=mask;
    }

  *out=value;

  delete[] weight;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename TYPE> int map_smooth_local_template(const grid_t &grid, const TYPE *buf, const TYPE mask, const float *weight, const int nx, const int ny, const int k, const int l, TYPE *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j;
  int imin,imax,jmin,jmax;
  float w;
  float sumW;
  TYPE sumZ;
  TYPE value,z;

  sumW=0.;
  sumZ=0.;

  imin=max(0,k-nx+1);
  imax=min(grid.nx,k+nx);
  jmin=max(0,l-ny+1);
  jmax=min(grid.ny,l+ny);

  for(j=jmin;j<jmax;j++) {
    const TYPE *pointer=&buf[j*grid.nx];
    const int offset=abs(j-l)*nx;
    for(i=imin;i<imax;i++) {
      z=pointer[i];
      if(z!=mask) {
        w=weight[offset+abs(i-k)];
        sumW+=w;
        sumZ+=z*w;
        }
      }
    }

  if(sumW!=0.) {
    value=sumZ/sumW;
    }
  else {
    value=mask;
    }

  *out=value;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_smooth_weigths(const grid_t & grid, const float scale, const int k, const int l, float * &weight, int & nx, int & ny)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j;
  float lx,ly,q,dum;
  float c,scale_x,scale_y;
  bool spherical=(grid.projection==0);
  
  c=1.0;
  if(spherical) c=cos(map_grid_y(grid,k,l)*d2r);

  scale_x=scale*c;
  scale_y=scale;
  
  nx=int(NINT(scale_x/grid.dx+1)+1);
  ny=int(NINT(scale_y/grid.dy+1)+1);

  nx=min(nx,grid.nx);
  ny=min(ny,grid.ny);

  weight=new float[nx*ny];

  for(j=0;j<ny;j++)
    for (i=0;i<nx;i++)
      weight[j*nx+i]=0.;

  for(j=0;j<ny;j++)
    for(i=0;i<nx;i++) {
      lx=i*grid.dx/scale_x;
      ly=j*grid.dy/scale_y;
      q=sqrt(lx*lx+ly*ly);
      if (q <= 1.) {
        dum=1-q*q*q;
        weight[j*nx+i]=dum*dum*dum;
        }
      }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_smooth_weigths(int mode, const grid_t & grid, const float scale, const int k, const int l, filter_t<float> & filter)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int nxmax,nymax;
  double dx,  dy;
  double c,scale_x,scale_y;
  bool spherical=(grid.projection==0);
  
  c=1.0;
  if(spherical) {
    c=cos(map_grid_y(grid,k,l)*d2r);
    status=map_resolution(grid, k, l, dx, dy, 0);
    }
  else {
    dx=grid.dx;
    dy=grid.dy;
    }
  
  scale_x=scale*c;
  scale_y=scale;
  
  nxmax=grid.nx/2;
  nymax=grid.ny/2;
  
  status=filter2D_InitWeigths(mode, scale_x, scale_y, dx, dy, nxmax, nymax, filter);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_smooth_local(const grid_t &grid, float *buf, float mask, float scale, int k, int l, float *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int nx,ny;
  float *weight=0;
  status=map_smooth_weigths(grid, scale, k, l, weight, nx, ny);
  status=map_smooth_local_template(grid, buf, mask, weight, nx, ny, k, l, out);
  delete[] weight;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_smooth_local(const grid_t &grid, complex<float> *buf, complex<float> mask, float scale, int k, int l, complex<float> *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int nx,ny;
  float *weight=0;
  status=map_smooth_weigths(grid, scale, k, l, weight, nx, ny);
  status=map_smooth_local_template(grid, buf, mask, weight, nx, ny, k, l, out);
  delete[] weight;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int map_smooth_template(int mode, const grid_t & grid, const T *buf, const T mask, const float scale, T *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int j;
  
  int nRequestedProcs=-1;
  int nprocs=initialize_OPENMP(nRequestedProcs, 0);
  
#pragma omp parallel for private(status) if(nprocs>1)
  for(j=0;j<grid.ny;j++) {
//     float *weight=0;
    int i=0/*, nx, ny*/;
    filter_t<float> filter;
    
//     status=map_smooth_weigths(grid, scale, i, j, weight, nx, ny);
    
    status=map_smooth_weigths(mode, grid, scale, i, j, filter);
    
    for(i=0;i<grid.nx;i++) {
      const int m=j*grid.nx+i;
      if(buf[m]!=mask) {

//         status=map_smooth_local_template(grid, buf, mask, (const float *) weight, (const int) nx, (const int) ny, (const int) i, (const int) j, &out[m]);

        status=map_smooth_local_template(grid, buf, mask, filter.weights, filter.nx, filter.ny, (const int) i, (const int) j, &out[m]);

        }
      else {
        out[m]=mask;
        }
      }
      
//     if(weight!=0) delete[] weight;
    
    filter.clear();
    
    }
    
  return 0;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_smooth(int mode, const grid_t & grid, const float *buf, const float mask, const float scale, float *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=map_smooth_template(mode, grid, buf, mask, scale, out);
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_smooth(int mode, const grid_t & grid, const double *buf, const double mask, const float scale, double *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=map_smooth_template(mode, grid, buf, mask, scale, out);
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_smooth(int mode, const grid_t & grid, const complex<float> *buf, const complex<float> mask, const float scale, complex<float> *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=map_smooth_template(mode, grid, buf, mask, scale, out);
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int map_smooth_template(const grid_t & grid, const T *buf, const T mask, const float *scale, T *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int j;
  
  int nRequestedProcs=-1;
  int nprocs=initialize_OPENMP(nRequestedProcs, 0);
  
#pragma omp parallel for private(status) if(nprocs>1)
  for(j=0;j<grid.ny;j++) {
    float *weight=0;
    for(int i=0;i<grid.nx;i++) {
      const int m=j*grid.nx+i;
      int nx, ny;
      if(buf[m]==mask) {
        out[m]=mask;
        continue;
        }
      if( (i>1) && (scale[m]==scale[m-1]) ) {
        }
      else {
        delete[] weight;
        status=map_smooth_weigths(grid, scale[m], i, j, weight, nx, ny);
        }
      status=map_smooth_local_template(grid, buf, mask, (const float *) weight, (const int) nx, (const int) ny, (const int) i, (const int) j, &out[m]);
      }
    delete[] weight;
    }
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_smooth(const grid_t & grid, const float *buf, const float mask, const float *scale, float *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=map_smooth_template(grid, buf, mask, scale, out);
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_smooth(const grid_t & grid, const double *buf, const double mask, const float *scale, double *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=map_smooth_template(grid, buf, mask, scale, out);
  return(status);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int map_smooth_latThenLong_template(const grid_t & grid, const T *buf, const T mask, const float scale, T *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j;
  struct timeval stv;   //start timeval, for progression
  
  int nRequestedProcs=1;
  int nprocs=initialize_OPENMP(nRequestedProcs, 0);
  
  gettimeofday(&stv);
  
  T *b=new T[grid.nx*grid.ny];
  
#pragma omp parallel for private(j) schedule(dynamic,1) if(nprocs>1)
  for(j=0;j<grid.ny;j++){
    int i,ii,i3,m,m2;
    double w;                        //weight
    T s;               //sum
//     #pragma omp critical
//     {STDERR_BASE_LINE_FUNC("1/2,%d/%d:%04.3gs; \r",j+1,grid.ny,difftime(stv));}

    double dx=max(MeanEarthRadius*grid.dx*cos(map_grid_y(grid,0,j)*d2r)*d2r,1000.);
    int n=int(NINT(scale/dx)+1);
    n=min(n,grid.nx);
    double *ws=new double[n];         //weights

/**----------------------------------------------------------------------------
    filter in longitudinal direction */
    for(i=0;i<n;i++){
      w=i*dx/scale;
      if(w<=1.){
        w=1-w*w*w;
        w=w*w*w;
        ws[i]=w;
        }
      else{
        ws[i]=0.;
        }
      //ws[i]=i?0:1;//deactivation
      }
    
    for(i=0;i<grid.nx;i++){
      m=j*grid.nx+i;
      if(buf[m]!=mask) {
        s=0.;
        w=0.;
        for(ii=1-n;ii<n;ii++){
          i3=i+ii;
/**----------------------------------------------------------------------------
          periodicity */
          if (grid.circular==1) {
            while (i3<0)  i3+=grid.nx-1;
            }
          else {
            if (i3<0) continue;
            }
            
          if (grid.circular==1) {
            while (i3>=grid.nx) i3-=grid.nx-1;
            }
          else {
            if (i3>=grid.nx) continue;
            }
          m2=j*grid.nx+i3;
          if(buf[m2]!=mask){
            i3=abs(ii);
            s+=buf[m2]*(float)ws[i3];
            w+=ws[i3];
            }
          }
        b[m]=s/(float)w;
        if(w==0.0) {
          out[m]=mask;
          }
        }
      else {
        b[m]=mask;
        }
      //b[m]=buf[m];//deactivation
      //b[m]=m==(grid.ny/2)*grid.nx?1:0;//debug
      }
    
    delete[] ws;
    }
  
/**----------------------------------------------------------------------------
  filter in latitudinal direction, weigth is not dependent upon latitude */
  double dy=max(MeanEarthRadius*grid.dy*d2r,1000.);
  int n=int(NINT(scale/dy+1)+1);
  n=min(n,grid.ny);
  double *ws=new double[n];//weights
  double w;//weight
  int i;
  for(i=0;i<n;i++){
    w=i*dy/scale;
    if(w<=1.){
      w=1-w*w*w;
      w=w*w*w;
      ws[i]=w;
      }
    else{
      ws[i]=0.;
      }
    //ws[i]=i>2?0:1;//debug
    }
  
#pragma omp parallel for private(j,i,w) schedule(dynamic,1) if(nprocs>1)
  for(i=0;i<grid.nx;i++){
//     #pragma omp critical
//     {STDERR_BASE_LINE_FUNC("2/2,%d/%d:%04.3gs; \r",j+1,grid.nx,difftime(stv));}
    int j,jj,j3,m,m2;
    T s;//sum
    
    for(j=0;j<grid.ny;j++){
      m=j*grid.nx+i;
      if(b[m]!=mask) {
        s=0.;
        w=0.;
        for(jj=1-n;jj<n;jj++){
          j3=j+jj;
          if( (j3<0) || (j3>=grid.ny) ) continue;
          m2=j3*grid.nx+i;
          if(b[m2]!=mask){
            j3=abs(jj);
            s+=b[m2]*(float)ws[j3];
            w+=ws[j3];
            }
          }
        out[m]=s/(float)w;
        if(w==0.0) {
          out[m]=mask;
          }
        }
      else {
        out[m]=mask;
        }
      //out[m]=b[m];//deactivation
      }
    }
  delete[] ws;
  delete[] b;
  
  fprintf(stderr,"\n");
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_smooth_latThenLong(const grid_t & grid, const complex<float> *buf, const complex<float> mask, const float scale, complex<float> *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=map_smooth_latThenLong_template(grid, buf, mask, scale, out);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_smooth_latThenLong(const grid_t & grid, const float *buf, const float mask, const float scale, float *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=map_smooth_latThenLong_template(grid, buf, mask, scale, out);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_diffusion_leapfrog(const grid_t & grid, const complex<float> *buf, const complex<float> mask, const float scale, complex<float> *out, int niteration)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------

  Smoothing by diffusion operator:
  --------------------------------
  
    dC/dt = nu (d²C/dx²+d²C/dy²)
      
-----------------------------------------------------------------------------*/
{
  int status=0;
  int j,k;
  complex<double> zz;
  const float R=6.3e+06,dx=R*grid.dx*d2r;
  const float nu=scale;
  complex<double> *state[3],*swap,dmask=mask;
  const float L=2*M_PI*sqrt(nu);
  const float CFL=1./(8.*nu/dx/dx);
//  const float tau=180.;
//  const float tau=240;
  const double tau=CFL/3.13;
  const double coef=nu*tau/(dx*dx);
  
  for(k=0;k<3;k++) state[k]=new complex<double> [grid.nx*grid.ny];
  
  printf("typical lengths=%f m smoothing lengths=%f km\n\n",L,L*sqrt(tau*niteration)/1000.);
  printf("CFL=%f time step=%f\n\n",CFL,tau);

  int nRequestedProcs=-1;
  int nprocs=initialize_OPENMP(nRequestedProcs);
  
#pragma omp parallel for if(nprocs>1)
  for(j=0;j<grid.ny;j++) {
    int i;
    for(i=0;i<grid.nx;i++) {
      const int m=j*grid.nx+i;
      state[0][m]=buf[m];
      state[1][m]=buf[m];
      state[2][m]=buf[m];
      out[m]=buf[m];
      }
    }
  for(int iteration=0;iteration<niteration;iteration++) {
    printf("iteration %6d\n",iteration);
#pragma omp parallel for private(zz) if(nprocs>1)
    for(j=1;j<grid.ny-1;j++) {
      int i;
      for(i=1;i<grid.nx-1;i++) {
        const int m=j*grid.nx+i;
        if(out[m]!=mask) {
          const int m1=m-1;
          const int m2=m+1;
          const int n1=m-grid.nx;
          const int n2=m+grid.nx;
#if 1
          const complex<double> z1= 0.5*(state[0][m1]+state[1][m1]);
          if(z1==dmask) continue;
          const complex<double> z2= 0.5*(state[0][m2]+state[1][m2]);
          if(z2==dmask) continue;
          const complex<double> z3= 0.5*(state[0][n1]+state[1][n1]);
          if(z3==dmask) continue;
          const complex<double> z4= 0.5*(state[0][n2]+state[1][n2]);
          if(z4==dmask) continue;
          zz=z1+z2+z3+z4-2.*(state[0][m]+state[1][m]);
#else
          const int frame=0;
          const complex<double> z1= state[frame][m1];
          if(z1==dmask) continue;
          const complex<double> z2= state[frame][m2];
          if(z2==dmask) continue;
          const complex<double> z3= state[frame][n1];
          if(z3==dmask) continue;
          const complex<double> z4= state[frame][n2];
          if(z4==dmask) continue;
          zz=z1+z2+z3+z4-4.* state[frame][m];
#endif
          state[2][m]=state[0][m]+2.*coef*zz;
          }
        }
      }
#pragma omp parallel for if(nprocs>1)
    for(j=1;j<grid.ny-1;j++) {
      int i;
      for(i=1;i<grid.nx-1;i++) {
        const int m=j*grid.nx+i;
        if(out[m]!=mask) {
          state[1][m]= 0.9*state[1][m]+ 0.05*(state[0][m]+state[2][m]);
          }
        }
      }
    swap=state[0];
    state[0]=state[1];
    state[1]=state[2];
    state[2]=swap;
    }
    
#pragma omp parallel for if(nprocs>1)
  for(j=0;j<grid.ny;j++) {
    int i;
    for(i=0;i<grid.nx;i++) {
      const int m=j*grid.nx+i;
      out[m]=state[2][m];
      }
    }
    
  for(k=0;k<3;k++) delete[] state[k];
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_diffusion_forward01_template(const grid_t & grid, const T *buf, const T mask, const float scale, T *out, int niteration)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  Smoothing by diffusion operator:
  --------------------------------
  
    dC/dt = nu (d²C/dx²+d²C/dy²)
  
    C+ = C- + tau * nu * (d²C/dx²+d²C/dy²)
  
  CFL condition:
  --------------
  
    dt < 1/(2*nu/dx²)
    
    Leapfrog : 2*dt*nu/dx²<1/4
    
  Length scale:
  -------------
  
    C=c(t) exp(ikx) -->  dc/dt -nu k² c=0 --> c=c(0,k) exp(-nu k² t)
    
    k²nu=1 --> (2pi/L)²nu=1
    
    L=2pi sqrt(nu)
    
------------------------------------------------------------------------------*/
{
  int status=0;
  int j,k;
  T zz;
  const double R=6.3e+06,dx=R*grid.dx*d2r;
  const double nu=scale;
  T *state[3],*swap,dmask=mask;
  const double L=2*M_PI*sqrt(nu);
  const double CFL=1./(4.*nu/dx/dx);
//  const float tau=180.;
  const double tau=CFL/1.13;
  const double coef=nu*tau/(dx*dx);
  
  for(k=0;k<2;k++) state[k]=new T [grid.nx*grid.ny];  
    
  printf("\ntypical lengths=%f m smoothing lengths=%f km\n",L,L*sqrt(tau*niteration)/1000.0);
  printf("CFL=%f s time step=%f s\n\n",CFL,tau);

  int nRequestedProcs=-1;
  int nprocs=initialize_OPENMP(nRequestedProcs);
  
#pragma omp parallel for if(nprocs>1)
  for(j=0;j<grid.ny;j++) {
    int i;
    for(i=0;i<grid.nx;i++) {
      const int m=j*grid.nx+i;
      state[0][m]=(T) buf[m];
      state[1][m]=(T) buf[m];
      out[m]=buf[m];
//       if(buf[m]!=mask) {
//         state[0][m]=buf[m];
//         state[1][m]=buf[m];
//         out[m]=buf[m];
//         }
//       else {
//         state[0][m]=0;
//         state[1][m]=0;
//         out[m]=0;
//         }
      }
    }
  for(int iteration=0;iteration<niteration;iteration++) {
    printf("iteration %6d\n",iteration);
#pragma omp parallel for private(zz) if(nprocs>1)
    for(j=1;j<grid.ny-1;j++) {
      int i;
      for(i=1;i<grid.nx-1;i++) {
        const int m=j*grid.nx+i;
        int done=1;
        if(out[m]!=mask) {
          done=0;
          const int m1=m-1;
          const int m2=m+1;
          const int n1=m-grid.nx;
          const int n2=m+grid.nx;
          const int frame=0;
          const T z1= state[frame][m1];
          if(z1==dmask) continue;
          const T z2= state[frame][m2];
          if(z2==dmask) continue;
          const T z3= state[frame][n1];
          if(z3==dmask) continue;
          const T z4= state[frame][n2];
          if(z4==dmask) continue;
          zz=z1+z2+z3+z4-(T)4.* state[frame][m];
          state[1][m]=state[0][m]+(T)coef*zz;
          done=1;
//           const int mm1=m-grid.nx-1;
//           const int mm2=m+grid.nx+1;
//           const int nn1=m-grid.nx+1;
//           const int nn2=m+grid.nx-1;
//           const T zz1= state[frame][mm1];
//           if(zz1==dmask) continue;
//           const T zz2= state[frame][mm2];
//           if(zz2==dmask) continue;
//           const T zz3= state[frame][nn1];
//           if(zz3==dmask) continue;
//           const T zz4= state[frame][nn2];
//           if(zz4==dmask) continue;
//           zz=zz1+zz2+zz3+zz4-4.* state[frame][m];
//           state[1][m]+=coef*zz/2.0;
          }
        if(done==0) {
          double w=0.0;
          zz=0;
          const int m1=m-1;
          const int m2=m+1;
          const int n1=m-grid.nx;
          const int n2=m+grid.nx;
          const int frame=0;
          const T z1= state[frame][m1];
          if(z1!=dmask) {
            zz+=z1;
            w++;
            }
          const T z2= state[frame][m2];
          if(z2!=dmask) {
            zz+=z2;
            w++;
            }
          const T z3= state[frame][n1];
          if(z3!=dmask) {
            zz+=z3;
            w++;
            }
          const T z4= state[frame][n2];
          if(z4!=dmask) {
            zz+=z4;
            w++;
            }
          const int mm1=m-grid.nx-1;
          const int mm2=m+grid.nx+1;
          const int nn1=m-grid.nx+1;
          const int nn2=m+grid.nx-1;
          const T zz1= state[frame][mm1];
          if(zz1!=dmask) {
            zz+=zz1;
            w++;
            }
          const T zz2= state[frame][mm2];
          if(zz2!=dmask) {
            zz+=zz2;
            w++;
            }
          const T zz3= state[frame][nn1];
          if(zz3!=dmask) {
            zz+=zz3;
            w++;
            }
          const T zz4= state[frame][nn2];
          if(zz4!=dmask) {
            zz+=zz4;
            w++;
            }
          zz+=state[frame][m];
          w++;
          state[1][m]=zz/(T) w;
          }
        }
      }
    swap=state[0];
    state[0]=state[1];
    state[1]=swap;
    }
    
#pragma omp parallel for if(nprocs>1)
  for(j=0;j<grid.ny;j++) {
    int i;
    for(i=0;i<grid.nx;i++) {
      const int m=j*grid.nx+i;
      out[m]=state[0][m];
      }
    }
    
  for(k=0;k<2;k++) delete[] state[k];
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_diffusion_forward02_template(const grid_t & grid, const T *buf, const T mask, const float scale, T *out, int niteration)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  Smoothing by diffusion operator:
  --------------------------------
  
    dC/dt = nu (d²C/dx²+d²C/dy²)
  
    C+ = C- + tau * nu * (d²C/dx²+d²C/dy²)
  
  CFL condition:
  --------------
  
    dt < 1/(2*nu/dx²)
    
    Leapfrog : 2*dt*nu/dx²<1/4
    
  Length scale:
  -------------
  
    C=c(t) exp(ikx) -->  dc/dt -nu k² c=0 --> c=c(0,k) exp(-nu k² t)
    
    k²nu=1 --> (2pi/L)²nu=1
    
    L=2pi sqrt(nu)
    
------------------------------------------------------------------------------*/
{
  int status=0;
  int j,k;
  T zz;
  const double nu=scale;
  T *state[3],*swap,dmask=mask;
  const double L=2*M_PI*sqrt(nu);
  double *dx, *dy;
  
  status=map_spherical_resolution(grid, &dx, &dy);
  
  range_t<double> rx;
  range_t<double> ry;
  for(size_t m=0;m<grid.Hsize();m++) {
    if(buf[m]==mask) continue;
    rx<<dx[m];
    ry<<dy[m];
    }
    
  const double dxmin=min(rx.min,ry.min);
  const double CFL=1./(4.*nu/dxmin/dxmin);
  const double tau=CFL/1.13;
  const double coef=nu*tau;
  
  for(k=0;k<2;k++) state[k]=new T [grid.nx*grid.ny];  
    
  printf("\ntypical lengths=%f m smoothing lengths=%f km\n",L,L*sqrt(tau*niteration)/1000.0);
  printf("CFL=%f s time step=%f s\n\n",CFL,tau);

  int nRequestedProcs=-1;
  int nprocs=initialize_OPENMP(nRequestedProcs);
  int milestone=max(1,niteration/100);
  
#pragma omp parallel for if(nprocs>1)
  for(j=0;j<grid.ny;j++) {
    int i;
    for(i=0;i<grid.nx;i++) {
      const int m=j*grid.nx+i;
      state[0][m]=(T) buf[m];
      state[1][m]=(T) buf[m];
      out[m]=buf[m];
      }
    }
  for(int iteration=0;iteration<niteration;iteration++) {
    const int frame=0;
    if(iteration % milestone ==0) printf("iteration %6d\n",iteration);
#pragma omp parallel for private(zz) if(nprocs>1)
    for(j=0;j<grid.ny;j++) {
      int i;
      for(i=0;i<grid.nx;i++) {
        const int m=j*grid.nx+i;
        const T z=state[frame][m];
        if(z!=mask) {
          T dx2=(T) (dx[m]*dx[m]);
          T dy2=(T) (dy[m]*dy[m]);
          zz=0;
          const int m1=m-1;
          const int m2=m+1;
          const int n1=m-grid.nx;
          const int n2=m+grid.nx;
//           const T z1= state[frame][m1];
          T z1 (i!=0 ?  z1= state[frame][m1] :  z1= dmask);
          if(z1!=dmask) {
            zz+=(z1-z)/dx2;
            }
//           const T z2= state[frame][m2];
          T z2 (i!=grid.nx-1 ?  z2= state[frame][m2] :  z2= dmask);
          if(z2!=dmask) {
            zz+=(z2-z)/dx2;
            }
//           const T z3= state[frame][n1];
          T z3 (j!=0 ?  z3= state[frame][n1] :  z3= dmask);
          if(z3!=dmask) {
            zz+=(z3-z)/dy2;
            }
//           const T z4= state[frame][n2];
          T z4 (j!=grid.ny-1 ?  z4= state[frame][n2] :  z4= dmask);
          if(z4!=dmask) {
            zz+=(z4-z)/dy2;
            }
          state[1][m]=state[0][m]+(T)coef*zz;
          }
        }
      }
    swap=state[0];
    state[0]=state[1];
    state[1]=swap;
    }
    
#pragma omp parallel for if(nprocs>1)
  for(j=0;j<grid.ny;j++) {
    int i;
    for(i=0;i<grid.nx;i++) {
      const int m=j*grid.nx+i;
      out[m]=state[0][m];
      }
    }
    
  for(k=0;k<2;k++) delete[] state[k];
  delete[] dx;
  delete[] dy;
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_diffusion_forward(const grid_t & grid, const complex<float> *buf, const complex<float> mask, const float scale, complex<float> *out, int niteration)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  switch(grid.modeH)  {
    case 0:
      status=map_diffusion_forward01_template(grid, buf, mask, scale, out, niteration);
      break;
    default:
      status=map_diffusion_forward02_template(grid, buf, mask, scale, out, niteration);
      break;
    }
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_diffusion_forward(const grid_t & grid, const float *buf, const float mask, const float scale, float *out, int niteration)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  switch(grid.modeH)  {
    case 0:
      status=map_diffusion_forward01_template(grid, buf, mask, scale, out, niteration);
      break;
    default:
      status=map_diffusion_forward02_template(grid, buf, mask, scale, out, niteration);
      break;
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_diffusion(const grid_t grid, const complex<float> *buf, const complex<float> mask, const float scale, complex<float> *out, int niteration)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------

  Smoothing by diffusion operator:
  --------------------------------
  
    dC/dt = nu (d²C/dx²+d²C/dy²)
  
    C+ = C- + tau * nu * (d²C/dx²+d²C/dy²)
  
  CFL condition:
  --------------
  
    dt < 1/(2*nu/dx²)
    
    Leapfrog : 2*dt*nu/dx²<1/4
    
  Length scale:
  -------------
  
    C=c(t) exp(ikx) -->  dc/dt -nu k² c=0 --> c=c(0,k) exp(-nu k² t)
    
    k²nu=1 --> (2pi/L)²nu=1
    
    L=2pi sqrt(nu)
    
    
    
-----------------------------------------------------------------------------*/
{
  int status;
  status=map_diffusion_forward(grid,  buf, mask, scale, out, niteration);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int map_average_template (const grid_t &grid, T *buffer, T mask, grid_t target, int m0, int n0, T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *-----------------------------------------------------------------------------

  return average value of buffer over a target cell

------------------------------------------------------------------------------*/
{
  int k,l,m,n;
  int i0,j0,mm;
  int status;
  double count;
  double x,y,xx,yy;
  T tmp;
  int i_range, j_range;

  count=0;
  tmp=0;
  
/* *-----------------------------------------------------------------------------
  get position of target(m0,n0)*/
  x=map_grid_x(target,m0,n0);
  y=map_grid_y(target,m0,n0);
/* *-----------------------------------------------------------------------------
  get grid indices of (x,y)*/
  status=map_index(grid, x, y, &i0, &j0);
  
/* *-----------------------------------------------------------------------------
  fix range, should be computed from grid and target resolutions*/
  i_range=max(target.dx/grid.dx,1.);
  j_range=max(target.dy/grid.dy,1.);
  
  for(l=max(0,j0-j_range);l<min(grid.ny-1,j0+j_range);l++) {
    for(k=max(0,i0-i_range);k<min(grid.nx-1,i0+i_range);k++) {
      xx=map_grid_x(grid,k,l);
      yy=map_grid_y(grid,k,l);
/* *-----------------------------------------------------------------------------
      get target indices of (xx,yy)*/
      status=map_index(target,xx,yy,&m,&n);
      if(status==0) {
        if((m!=m0)||(n!=n0)) continue;
        mm=l*grid.nx+k;
        if(buffer[mm]!=mask) {
          tmp+=buffer[mm];
          count++;
          }
        }
      }
    }

  if(count!=0) {
    tmp/=count;
    *z=tmp;
    status=0;
    }
  else {
    *z=mask;
    status=-1;
    }
    
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_average (const grid_t &grid, float *buffer, float mask, grid_t target, int m0, int n0, float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=map_average_template (grid, buffer, mask, target, m0, n0, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_average (const grid_t &grid, double *buffer, double mask, grid_t target, int m0, int n0, double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=map_average_template (grid, buffer, mask, target, m0, n0, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_average (const grid_t &grid, complex<float> *buffer, complex<float> mask, grid_t target, int m0, int n0, complex<float> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=map_average_template (grid, buffer, mask, target, m0, n0, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_average (const grid_t &grid, complex<double> *buffer, complex<double> mask, grid_t target, int m0, int n0, complex<double> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=map_average_template (grid, buffer, mask, target, m0, n0, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> void distance_to_nearest_unmasked_template(const grid_t &grid,const T *buffer,T mask,double *distance, T *value=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// distance in pixels to nearest unmasked node
/** \note This assumes the grid is orthonormal!!!
*/
/*----------------------------------------------------------------------------*/
{
  const int n=grid.Hsize();
  int k,i0,j0,ii,ji,ie,je,i,j,m,id,jd,md;
  
  for(m=0;m<n;m++){
    double *dm=&distance[m];
    const T * bm=&buffer[m];
    
    if(*bm!=mask and not isnan(*bm)){
      *dm=0;
      if(value!=0)
        value[m]=*bm;
      }
    else{
      *dm=+INFINITY;
      if(value!=0)
        value[m]=mask;
      }
    
    }
  
  double dd,d00,d01,d10,d11,X,Y,b,c,D;
  bool updated;
  
  for(k=0;k<4;k++){
    
    if(k&1){
      i0=0;
      ii=1;
      ie=grid.nx;
      }
    else{
      i0=grid.nx-1;
      ii=-1;
      ie=-1;
      }
    
    if(k&2){
      j0=0;
      ji=1;
      je=grid.ny;
      }
    else{
      j0=grid.ny-1;
      ji=-1;
      je=-1;
      }
    
    for(j=j0;(je-j)*ji>0;j+=ji)
      for(i=i0;(ie-i)*ii>0;i+=ii){
        
        m=i+j*grid.nx;
        double *dm=&distance[m];
  /*------------------------------------------------------------------------------
        find masked */
        if(*dm<=0.)
          continue;
        
  /*------------------------------------------------------------------------------
        compute distance */
        d00=+INFINITY;
        
        for(jd=j-ji;(j+ji-jd)*ji>0;jd+=ji){
          
          if(jd<0 or jd>=grid.ny)
            continue;
          
          for(id=i-ii;(i+ii-id)*ii>0;id+=ii){
          
            if(id<0 or id>=grid.nx)
              continue;
            
            if(id==i and jd==j)
              continue;
            
            md=id+jd*grid.nx;
            const double *dmd=&distance[md];
            
            if(id!=i and jd!=j)
              d00=*dmd;
            else if(id==i)
              d01=*dmd;
            else if(jd==j)
              d10=*dmd;
            
            if(*dm<*dmd)
              continue;
            
            if(id==i or jd==j)
              dd=1;
            else
              dd=M_SQRT2;
            
            updated=updatemin(dm,*dmd+dd);
            
            if(updated and value!=0)
              value[m]=value[md];
            }
          }
        
        if(not isfinite(d00))
          continue;
  /*------------------------------------------------------------------------------
        make sure distance is better:
        
        d00 d01
        d10 d11
        
        The gradient is such as:
          ((d10+d11)*.5-(d00+d01)*.5)^2+((d01+d11)*.5-(d00+d10)*.5)^2=1
          <==> (d10+d11-d00-d01)^2+(d01+d11-d00-d10)^2=4
        with X=d10-d00-d01 and Y=d01-d00-d10 :
          <==> (X+d11)^2+(Y+d11)^2=4=X^2+Y^2+2*d11*(X+Y)+2*d11^2
          <==> 0=X^2+Y^2-4+2*(X+Y)*d11+2*d11^2
        with b=2*(X+Y) and c=X^2+Y^2-4 and D=b^2-8*c :
          d11=(-b + sqrt(D) )/4
        */
        X=d10-d00-d01;
        Y=d01-d00-d10;
        b=2*(X+Y);
        c=X*X+Y*Y-4;
        D=b*b-8*c;
        d11=(-b + sqrt(D) )/4;
        
        if(false and d00==0 and d10>0 and d01>0)
          STDERR_BASE_LINE_FUNC("[%d,%d]:%g %g, %g:%g\n",i0,j0,d00,d01,d10,d11);
        updatemin(dm,d11);
        }
    
    //TRAP_ERR_RETURN(,1,"testing\n");
    }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void distance_to_nearest_unmasked(const grid_t &grid,const double *buffer,double mask,double *distance,double *value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  distance_to_nearest_unmasked_template(grid,buffer,mask,distance,value);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void distance_to_nearest_unmasked(const grid_t &grid,const complex<double> *buffer,complex<double> mask,double *distance,complex<double> *value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  distance_to_nearest_unmasked_template(grid,buffer,mask,distance,value);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void distance_to_nearest_unmasked(const grid_t &grid,const float *buffer,float mask,double *distance,float *value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  distance_to_nearest_unmasked_template(grid,buffer,mask,distance,value);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void distance_to_nearest_unmasked(const grid_t &grid,const char *buffer,char mask,double *distance,char *value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  distance_to_nearest_unmasked_template(grid,buffer,mask,distance,value);
}
