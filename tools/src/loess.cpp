
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Detides comodo-compliant NetCDF outputs and produces tidal atlases.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "filter.h"

#include "solverlib.h"

#include "tools-define.h"

#include "loess.def"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void Loess1D_init(double dx, double l_c, int m, double *weight, int *nweight)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
Schlax, M. G., et D. B. Chelton (1992), Frequency Domain Diagnostics for Linear Smoothers,
Journal of the American Statistical Association, 87(420), 1070-1081.
-----------------------------------------------------------------------*/
{
  int i;
  double p,q,r;

  for(i=0;i<2*m+1;i++)
    weight[i]=0.;

  for(i=0;i<2*m+1;i++) {
    r=(double)(i-m)*dx/l_c;
/* *----------------------------------------------------------------------------
    Bug, correction 09/12/2010; see reference above */
//    q=r*r;
    q=fabs(r);
    if (q <= 1.) {
      p=1.-q*q*q;
      weight[i]=p*p*p;
      }
    else
      weight[i]=0.;
    }
  *nweight=2*m+1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void Loess1D_BF(int n,const double *h, double mask, int m,const double *weight, double *hf)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int i,j;
  double sum_w,sum_wh,w,z;

  for(i=0;i<n;i++) {
    sum_w=0.;
    sum_wh=0.;
    for(j=max(i-m,0);j<min(i+m,n);j++) {
      z=h[j];
      if(z != mask) {
        w=weight[m+i-j];
        sum_w =sum_w +w;
        sum_wh=sum_wh+w*z;
        }
      }
    if(sum_w != 0.0)
      hf[i]=sum_wh/sum_w;
    else
      hf[i]=mask;
    }
/* *-----------------------------------------------------------------------
  force lf to mask value if initial array is masked*/
  for(j=0;j<n;j++) if(h[j] == mask) hf[j]=mask;
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void Loess1D(int n, double dx, double l_c,const double *h, double mask,double *lf,int *rstatus)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *-----------------------------------------------------------------------
  return low-passed filtered serie*/
{
  int  m;
  double *weight=NULL;
  int nweight;

  m=min((int)(l_c/dx)+1,n);

  exitIfNull(
    weight=(double *) malloc((2*m+1)*sizeof(double))
    );

  Loess1D_init(dx,l_c,m,weight,&nweight);

  Loess1D_BF(n,h,mask,m,weight,lf);

  free(weight);
  *rstatus=0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void loess1d(float *buf, float mask, float step, int nval, float scale, float *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
Schlax, M. G., et D. B. Chelton (1992), Frequency Domain Diagnostics for Linear Smoothers,
Journal of the American Statistical Association, 87(420), 1070-1081.
-----------------------------------------------------------------------------*/
{
  int i,k;
  int imin,imax;
  float w;
  float lx,q,dum,sum,z;
  int n;
  float *weight=NULL,tmp;

  n=int(NINT(scale/step));
  n=min(nval,n);

  exitIfNull(
    weight=(float *) malloc(n*sizeof(float))
    );

  for (i=0;i<n;i++) weight[i]=0.;

  for(i=0;i<n;i++) {
      lx=i*step/scale;
/* *----------------------------------------------------------------------------
      Bug, correction 09/12/2010; see reference above */
//      q=lx*lx;
      q=fabs(lx);
      if (q <= 1.) {
        dum=1-q*q*q;
        weight[i]=dum*dum*dum;
        }
      }

  for(k=0;k<nval;k++) {
    tmp=0.;
    sum=0.;
    imin=max(0,k-n+1);
    imax=min(nval,k+n);
    for(i=imin;i<imax;i++) {
      z=buf[i];
      if(z!=mask) {
        w=weight[abs(i-k)];
        sum+=w;
        tmp+=z*w;
        }
      }
    if(sum!=0.) out[k]=tmp/sum;
    else out[k]=mask;
    }

  free(weight);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void Loess1D_BF(float *buf, float mask, int nval, float *weight, int nweight, float *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k;
  int imin,imax;
  float w;
  float sum,z;
  int n;
  float tmp;

/*-----------------------------------------------------------------------------

buf    : serie
mask   : flag de validit�nval   : nombre de valeur max de la s�ie
weight : le tableau des poids

en sortie out BF de buf

-----------------------------------------------------------------------------*/

  n=nweight;

  for(k=0;k<nval;k++) {
    tmp=0.;
    sum=0.;
    imin=max(0,k-n+1);
    imax=min(nval,k+n);
    for(i=imin;i<imax;i++) {
      z=buf[i];
      if(z!=mask) {
        w=weight[abs(i-k)];
        sum+=w;
        tmp+=z*w;
        }
      }
    if(sum!=0.) out[k]=tmp/sum;
    else out[k]=mask;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void Loess1D_BF_nomask(const float *buf, int nval, const float *weight, int nweight, float *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

buf    : serie
mask   : flag de validit�nval   : nombre de valeur max de la s�ie
weight : le tableau des poids

en sortie out BF de buf

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

{
  int i,k;
  int m,n;
  float sum, *ww=NULL;
  const float *z=NULL;

  n=nweight;
  m=2*n-1;

//  ww=new float[m];
  exitIfNull(
    ww=(float *) malloc(m*sizeof(float))
    );
  for(i=0;i<n;i++) {
    ww[i+n-1]=weight[i];
    ww[n-1-i]=weight[i];
    }

  sum=0;
  for(i=0;i<m;i++) sum+=ww[i];
  for(i=0;i<m;i++) ww[i]/=sum;;

  for(k=n-1;k<nval-n;k++) {
    z=&(buf[k-n+1]);
#ifdef BLASF_
    sum = sdot(m, z, 1, ww, 1);
#elif BLASF
    sum = sdot(m, z, 1, ww, 1);
#elif BLASC
    sum = sdot(m, z, 1, ww, 1);
#elif CBLAS
    sum = cblas_sdot(m, z, 1, ww, 1);
#else
    sum=0;
    for(i=0;i<m;i++) {
      sum+=ww[i]*z[i];
      }
#endif
    out[k]=sum;
    }
//  delete [] ww;
  free(ww);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void Loess1D_BF_nomask(const double *buf, int nval,const double *weight, int nweight, double *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

buf    : serie
mask   : flag de validit�nval   : nombre de valeur max de la s�ie
weight : le tableau des poids

en sortie out BF de buf

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

{
  int i,k;
  int m,n;
  double sum, *ww=NULL;
  const double *z=NULL;

  n=nweight;
  m=2*n-1;

//  ww=new float[m];
  exitIfNull(
    ww=(double *) malloc(m*sizeof(double))
    );
  for(i=0;i<n;i++) {
    ww[i+n-1]=weight[i];
    ww[n-1-i]=weight[i];
    }

  sum=0;
  for(i=0;i<m;i++) sum+=ww[i];
  for(i=0;i<m;i++) ww[i]/=sum;;

  for(k=n-1;k<nval-n;k++) {
    z=&(buf[k-n+1]);
#ifdef BLASF_
    sum = ddot(m, z, 1, ww, 1);
#elif BLASF
    sum = ddot(m, z, 1, ww, 1);
#elif BLASC
    sum = ddot(m, z, 1, ww, 1);
#elif CBLAS
    sum = cblas_ddot(m, z, 1, ww, 1);
#else
    sum=0;
    for(i=0;i<m;i++) {
      sum+=ww[i]*z[i];
      }
#endif
    out[k]=sum;
    }
//  delete [] ww;
  free(ww);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void Loess1D_init(float step, int nval, float scale, float **weight, int *nweight)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  int n;
  float lx,q,dum;

/*-----------------------------------------------------------------------------

step  :  pas d'�hantillonnage
nval  :  nombre de valeur max de la s�ie
scale :  demi �helle de coupure

en sortie: le tableau des poids

*-----------------------------------------------------------------------------*/
  n=int(NINT(scale/step));
  n=min(nval,n);

  *weight=aset(n,0.f);
  if( *weight == NULL ) n=-1;

  for(i=0;i<n;i++) {
    if(scale!=0) lx=i*step/scale;
    q=fabs(lx);
    if (q <= 1.) {
      dum=1-q*q*q;
      (*weight)[i]=dum*dum*dum;
      }
    }
  *nweight=n;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void Loess2D(grid_t grid, float *buf, float mask, float scale, float *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l;
  int imin,imax,jmin,jmax;
  float w;
  float lx,ly,q,dum,sum,z;
  int nx,ny,offset_j,offset_l;
  float *weight=NULL,tmp;

  nx=int(NINT(scale/grid.dx+1)+1);
  ny=int(NINT(scale/grid.dy+1)+1);

  nx=min(nx,grid.nx);
  ny=min(ny,grid.ny);

  exitIfNull(
    weight=(float *) malloc(nx*ny*sizeof(float))
    );

  for(j=0;j<ny;j++)
    for (i=0;i<nx;i++)
      weight[j*nx+i]=0.;

  for(j=0;j<ny;j++)
    for(i=0;i<nx;i++) {
      lx=i*grid.dx/scale;
      ly=j*grid.dy/scale;
/* *----------------------------------------------------------------------------
      Bug, correction 09/12/2010; see reference above */
//      q=lx*lx+ly*ly;
      q=sqrt(lx*lx+ly*ly);
      if (q <= 1.) {
        dum=1-q*q*q;
        weight[j*nx+i]=dum*dum*dum;
        }
      }

  for(l=0;l<grid.ny;l++) {
    offset_l=l*grid.nx;
    jmin=max(0,l-ny+1);
    jmax=min(grid.ny,l+ny);
    for(k=0;k<grid.nx;k++) {
      tmp=0.;
      sum=0.;
      imin=max(0,k-nx+1);
      imax=min(grid.nx,k+nx);
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
      if(sum!=0.) out[offset_l+k]=tmp/sum;
      else out[offset_l+k]=mask;
      }
    }

  free(weight);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *Loess2D_BF(grid_t grid, float *buf, float mask, float scale, float *weight, float *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l;
  int imin,imax,jmin,jmax;
  float w;
  float sum,z;
  int nx,ny,offset_j,offset_l,offset_w;
  float tmp;

  if (out==0) out=new float[grid.ny*grid.nx];

  nx=int(NINT(scale/grid.dx+1)+1);
  ny=int(NINT(scale/grid.dy+1)+1);

  nx=min(nx,grid.nx);
  ny=min(ny,grid.ny);

  for(l=0;l<grid.ny;l++) {
    offset_l=l*grid.nx;
    jmin=max(0,l-ny+1);
    jmax=min(grid.ny,l+ny);
    for(k=0;k<grid.nx;k++) {
      if(buf[offset_l+k]==mask) {
        out[offset_l+k]=mask;
        continue;
        }
      tmp=0.;
      sum=0.;
      imin=max(0,k-nx+1);
      imax=min(grid.nx,k+nx);
      for(j=jmin;j<jmax;j++) {
        offset_w=abs(j-l)*nx;
        offset_j=j*grid.nx;
        for(i=imin;i<imax;i++) {
          z=buf[offset_j+i];
          if(z!=mask) {
            w=weight[offset_w+abs(i-k)];
            sum+=w;
            tmp+=z*w;
            }
          }
        }
      if(sum!=0.) out[offset_l+k]=tmp/sum;
      else out[offset_l+k]=mask;
      }
    }
  return(out);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int Loess2D_init_cartesian(double scale_x, double scale_y, double azimuth, double dx, double dy, int nxmax, int nymax, float* & weight, int & nx, int & ny)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

  assume rectangular,cartesian grid with regular spacing (unfiform dx dy)

  scales must be given in grid coordinates' units

------------------------------------------------------------------------------*/
{
  int i,j;
  double lx,ly,q,dum;

  nx=int(NINT(scale_x/dx+1)+1);
  ny=int(NINT(scale_y/dy+1)+1);

  nx=min(nx,nxmax);
  ny=min(ny,nymax);

  weight=new float[nx*ny];

  exitIfNull(weight);
  
  for(j=0;j<ny;j++)
    for (i=0;i<nx;i++)
      weight[j*nx+i]=0.;

  for(j=0;j<ny;j++) {
    for(i=0;i<nx;i++) {
      lx=i*dx/scale_x;
      ly=j*dy/scale_y;
      q=sqrt(lx*lx+ly*ly);
      if (q <= 1.) {
        dum=1-q*q*q;
        weight[j*nx+i]=dum*dum*dum;
        }
      }
    }
  return(0);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int Loess2D_init_cartesian(double scale_x, double scale_y, double azimuth, double dx, double dy, int nxmax, int nymax, filter_t<float> & filter)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=Loess2D_init_cartesian(scale_x, scale_y, azimuth, dx, dy, nxmax, nymax, filter.weights, filter.nx, filter.ny);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void Loess2D_init_spherical(grid_t grid, float scale_x, float scale_y, float azimuth, float **weight)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

  assume rectangular,spherical grid with regular spacing (unfiform dlon dlat)

  scales must be given in meters, grid coordinates in degrees

------------------------------------------------------------------------------*/
{
  int i,j;
  float lx,ly,q,dum;
  int nx,ny;

  nx=int(NINT(scale_x/grid.dx+1)+1);
  ny=int(NINT(scale_y/grid.dy+1)+1);

  nx=min(nx,grid.nx);
  ny=min(ny,grid.ny);

  exitIfNull(
    *weight=(float *) malloc(nx*ny*sizeof(float))
    );

  for(j=0;j<ny;j++)
    for (i=0;i<nx;i++)
      (*weight)[j*nx+i]=0.;

  for(j=0;j<ny;j++) {
    for(i=0;i<nx;i++) {
      lx=i*grid.dx/scale_x;
      ly=j*grid.dy/scale_y;
      q=sqrt(lx*lx+ly*ly);
      if (q <= 1.) {
        dum=1-q*q*q;
        (*weight)[j*nx+i]=dum*dum*dum;
        }
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void Loess2D_init(grid_t grid, float scale, float **weight)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j;
  float lx,ly,q,dum;
  int nx,ny;

  nx=int(NINT(scale/grid.dx+1)+1);
  ny=int(NINT(scale/grid.dy+1)+1);

  nx=min(nx,grid.nx);
  ny=min(ny,grid.ny);

  exitIfNull(
    *weight=(float *) malloc(nx*ny*sizeof(float))
    );

  for(j=0;j<ny;j++)
    for (i=0;i<nx;i++)
      (*weight)[j*nx+i]=0.;

  for(j=0;j<ny;j++)
    for(i=0;i<nx;i++) {
      lx=i*grid.dx/scale;
      ly=j*grid.dy/scale;
      q=sqrt(lx*lx+ly*ly);
      if (q <= 1.) {
        dum=1-q*q*q;
        (*weight)[j*nx+i]=dum*dum*dum;
        }
      }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int loess1d_template(int n, double l_c, T *h, double *x, T mask,T *lf)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j;
  int holes=0;
  double s,q,*w=NULL;

  if(n==1) lf[0]=h[0];/* cut is NaN there! */
  else {
    exitIfNull(
      w=(double *) malloc(n*sizeof(double))
      );
    l_c/=2;
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) {
        q=(x[j]-x[i])/l_c;
/* *----------------------------------------------------------------------------
        Bug, correction 09/12/2010; see reference above */
//        q*=q;
        q=fabs(q);
        s=1-q*q*q;
        if(q>1) w[j]=0;
        else w[j]=s*s*s;
        }
      s=0;
      lf[i]=0;
      for(j=0;j<n;j++)
        if(h[j]!=mask) {
          s+=w[j];
          lf[i]+=((float) w[j])*h[j];
          }
      if(s!=0) {
        lf[i]/=s;
        }
      else {
        lf[i]==mask;
        holes++;
        }
      }
    free(w);
    }
  return(holes);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int loess1d_irregular(int n, double l_c, double *h, double *x, double mask,double *lf)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, first=1,j0;
  
  int nRequestedProcs=-1;
  int nprocs __attribute__((unused))=initialize_OPENMP(nRequestedProcs,0);

  if(n==1) lf[0]=h[0];/* cut is NaN there! */
  else {
    l_c/=2;
    first=1;
    j0=0;
#pragma omp parallel for firstprivate(first,j0)
    for(i=0; i<n ; i++) {
      lf[i]=0;
      if(first==1) {
        j0=0;
        first=0;
        }
      double w,wtot=0,s,q;
      for(int j=j0;j<n;j++) {
        if(h[j]==mask) continue;
        q=(x[j]-x[i])/l_c;
        if(q<-1.) continue;
        if(q> 1.) break;
        if(wtot==0.0) j0=j;
        q=fabs(q);
        s=1-q*q*q;
        w=s*s*s;
        lf[i]+=w*h[j];
        wtot+=w;
        }
//       if(wtot!=0) {
      if(wtot>0.5) {
        lf[i]/=wtot;
        }
      else {
        lf[i]=mask;
        }
      }
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int loess1d_irregular(int n, float l_c, float *h, float *x, float mask,float *lf)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j;
  float s,q,*w=NULL;

/* *----------------------------------------------------------------------------
  n   : number of value in vector h (to be filtered)
  l_c : cut distance (100 km -> l_c=50 km)
  h   : input vector
  x   : position vector
  mask:
  lf  : low-pass fliter output
  rstatus: flag
------------------------------------------------------------------------------*/

  if(n==1) lf[0]=h[0];/* cut is NaN there! */
  else {
    exitIfNull(
      w=(float *) malloc(n*sizeof(double))
      );
    l_c/=2;
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) {
        q=(x[j]-x[i])/l_c;
/* *----------------------------------------------------------------------------
        Bug, correction 09/12/2010; see reference above */
//        q*=q;
        q=fabs(q);
        s=1-q*q*q;
        if(q>1) w[j]=0;
        else w[j]=s*s*s;
        }
      s=0;
      for(j=0;j<n;j++)  s+=w[j];
      lf[i]=0;
      for(j=0;j<n;j++) lf[i]+=w[j]*h[j];
      lf[i]/=s;
      }
    free(w);
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int loess1d(int n, double l_c, complex<float> *h, double *x, complex<float> mask, complex<float> *lf)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int holes;
  holes=loess1d_template(n, l_c, h,x, mask,lf);
  return(holes);
}


