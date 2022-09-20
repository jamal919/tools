
/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "constants.h"

#include "map.h"
#include "geo.h"
#include "statistic.h"
#include "map.def"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double map_grid_x (const grid_t & grid, int k, int l)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m;
  double x;

  switch(grid.modeH) {
    case 0:
      x=grid.xmin+k*grid.dx;
      break;

    case 1:
      x=grid.x[k];
      break;

    case 2:
      m=l*grid.nx+k;
      x=grid.x[m];
      break;

   default:
      TRAP_ERR_EXIT(-1,"exiting\n");
      break;
    }

  return (x);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double map_grid_y (const grid_t & grid, int k, int l)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m;
  double y;

  switch(grid.modeH) {
    case 0:
      y=grid.ymin+l*grid.dy;
      break;

    case 1:
      y=grid.y[l];
      break;

    case 2:
      m=l*grid.nx+k;
      y=grid.y[m];
      break;

   default:
      TRAP_ERR_EXIT(-1,"exiting\n");
      break;
    }

  return (y);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_persistence_template(const grid_t & grid, T *buf, T mask, double radius)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j;
  
  T *tmp=new T[grid.Hsize()];

#pragma omp parallel for
  for(j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      int m=j*grid.nx+i;
      tmp[m]=buf[m];
      }
    }
    
#pragma omp parallel for
  for(j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      int m=j*grid.nx+i;
      if(buf[m]==mask) {
        T count=0.0;
        T r=0;
        for(int jj=max(0,j-1);jj<min(grid.ny,j+2);jj++) {
          for(int ii=max(0,i-1);ii<min(grid.nx,i+2);ii++) {
            int n=jj*grid.nx+ii;
            if(buf[n]!=mask) {
              r+=buf[n];
              count+=1.;
              }
            }
//           if(tmp[m]!=mask) break;
          }
        if(count!=(T) 0) tmp[m]=r/count;
        }
      else {
        tmp[m]=buf[m];
        }
      }
    }
    
#pragma omp parallel for
  for(j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      int m=j*grid.nx+i;
      buf[m]=tmp[m];
      }
    }
    
  delete[] tmp;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_persistence(const grid_t & grid, signed char *buf, signed char mask, double radius)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=map_persistence_template(grid, buf, mask, radius);
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_persistence(const grid_t & grid, float *buf, float mask, double radius)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=map_persistence_template(grid, buf, mask, radius);
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_persistence(const grid_t & grid, double *buf, double mask, double radius)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=map_persistence_template(grid, buf, mask, radius);
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_persistence(const grid_t & grid, complex<float> *buf, complex<float> mask, double radius)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=map_persistence_template(grid, buf, mask, radius);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_extendgrid(grid_t *grid, int *modified)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------
  check if extending the grid is necessary, i.e. if grid covers the sphere
  but first column is not repeated (uneasy for interpolation)

  Ideally this should be imposed by the calling sequence
-----------------------------------------------------------------------*/
  int k, l, m, n, status;
  double *tmpx,*tmpy;
  double dx,x1,x2;

  *modified=0;

  if(grid->modeH!=0) {
    dx=grid->x[1]-grid->x[0];
    }
  else {
    dx=grid->dx;
    }
  x1=grid->x[0];
  x2=grid->x[grid->nx-1]+dx;

  if(fmod(x2-x1, (double) 360.0)!=0.0) {
    *modified=0;
    return(0);
    }

  *modified=1;

  grid->nx++;

  tmpx = new double[grid->nx*grid->ny];
  tmpy = new double[grid->nx*grid->ny];

/*------------------------------------------------------------------------------
  translate buffer*/
  for(l=grid->ny-1; l>-1; l--) {
    for(k=grid->nx-2; k>-1; k--) {
      m=l*(grid->nx-1)+k; /*original position*/
      n=l*grid->nx+k;     /*extended position*/
      tmpx[n] = grid->x[m];
      tmpy[n] = grid->y[m];
//       double x,y;
//       grid->xy(k,l,x,y);
//       tmpx[n] = x;
//       tmpy[n] = y;
      }
    }
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check : MANDATORY !!!

  Notes:

  20/04/2010:

    bug found:

  for(l=0; l<grid->ny; l++) {
    m=l*grid->nx;
    n=l*grid->nx+grid->nx-1;
    tmpx[n] = grid->x[m]+360.0;
    tmpy[n] = grid->y[m];
    }

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*------------------------------------------------------------------------------
  repeat first column in the last one*/
  for(l=0; l<grid->ny; l++) {
    m=l*(grid->nx-1);         /*first column*/
    n=l*grid->nx+grid->nx-1;  /*last column*/
    tmpx[n] = grid->x[m]+360.0;
    tmpy[n] = grid->y[m];
//     double x,y;
//     grid->xy(0,l,x,y);
//     tmpx[n] = x+360.0;
//     tmpy[n] = y;
    }

  delete[] grid->x;
  delete[] grid->y;

  grid->x=tmpx;
  grid->y=tmpy;

  grid->xmax+=grid->dx;
  
  status = 0;
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int map_extendbuffer_template(const grid_t & grid, T **buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l, m, n, status;
  T  *tmp;
    
  tmp =new T [grid.nx*grid.ny];

/*------------------------------------------------------------------------------
  translate buffer*/
  for(l=0; l<grid.ny; l++) {
    for(k=0; k<grid.nx-1; k++) {
      m=l*(grid.nx-1)+k; /*original position*/
      n=l*grid.nx+k;     /*extended position*/
      tmp[n] = (*buffer)[m];
      }
    }

/*------------------------------------------------------------------------------
  repeat first column in the last one*/
  for(l=0; l<grid.ny; l++) {
    m=l*(grid.nx-1);         /*first column*/
    n=l*grid.nx+grid.nx-1;  /*last column*/
    tmp[n] = (*buffer)[m];
    }

  delete[] (*buffer);

  *buffer=tmp;

  status = 0;
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int map_extendfield_template(grid_t *grid, T **buffer, int *modified)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  *modified=0;
  
  status=map_extendgrid(grid, modified);
  if(*modified==0) {
    status = 0;
    return (status);
    }
    
  status=map_extendbuffer_template(*grid, buffer);

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_extendfield(grid_t *grid, complex<float> **buffer, int *modified)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=map_extendfield_template(grid, buffer, modified);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int map_extendfield_template(grid_t *grid, T **buffer, int nbuffers, int *modified)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------
  
-----------------------------------------------------------------------*/
  int k, status;

  *modified=0;
  
  status=map_extendgrid(grid, modified);
  if(*modified==0) {
    status = 0;
    return (status);
    }
  for(k=0;k<nbuffers;k++) {
    status=map_extendbuffer_template(*grid, &(buffer[k]));
    }
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_extendfield(grid_t *grid, complex<float> **buffer, int nbuffers, int *modified)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=map_extendfield_template(grid, buffer, nbuffers, modified);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_mirror (const grid_t & grid, float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  float tmp;
  int i,j,m1,m2;

  for (j=0;j<grid.ny/2;j++)
    for (i=0;i<grid.nx;i++) {
      m1=grid.nx*j+i;
      m2=grid.nx*(grid.ny-j-1)+i;
      tmp=buffer[m1];
      buffer[m1]=buffer[m2];
      buffer[m2]=tmp;
      }
  return (0);
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_integraleR(const grid_t & grid, float *buf, float mask, double *sum,double *area)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,m;
  double w=0,a=0;
  double p;
  double s=0;

  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      p=map_grid_y(grid,i,j);
      if(buf[m]!=mask) {
        w=cos(p*d2r);
        s+=w* buf[m];
        a+=w;
        }
      }
    }

  *sum=s;
  *area=a;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_integraleC(const grid_t & grid, fcomplex *buf, fcomplex mask, dcomplex *sum,dcomplex *area)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,m;
  double w=0;
  dcomplex a=dcomplex(0.,0.);
  double p;
  dcomplex s=0,z;
  

  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      p=map_grid_y(grid,i,j);
      if(buf[m]!=mask) {
        w=cos(p*d2r);
        z=(dcomplex) buf[m];
        s+=z*dcomplex(w,0.);
        a+=w;
        }
      }
    }

  *sum=s;
  *area=a;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_export_template(const grid_t & g_source, T *b_source, T m_source, grid_t g_target, T *b_target, T m_target, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j,status;
  T z;
  int nprocs;
  
  if(mode<0){
    mode=0;
    }
  
  nprocs=initialize_OPENMP(-1);

  switch(mode) {
    case 0:
#pragma omp parallel for private(z,status) if(nprocs>1)
      for(j=0;j<g_target.ny;j++) {
        double t,p;
        for(int i=0;i<g_target.nx;i++) {
          int m=g_target.nx*j+i;
          g_target.xy(i,j,t,p);
          t=map_recale(g_source,t);
          status=map_interpolation(g_source,b_source,m_source,t,p,&z);
          if(status ==0) {
            b_target[m]=z;
            }
          else {
            b_target[m]=m_target;
            }
          }
        }
      break;
    case 1:
#pragma omp parallel for private(z,status) if(nprocs>1)
      for(j=0;j<g_target.ny;j++) {
        for(int i=0;i<g_target.nx;i++) {
          int m=g_target.nx*j+i;
//           double t,p;
//           g_target.xy(i,j,t,p);
          status=map_average(g_source,b_source,m_source,g_target,i,j,&z);
          if(status ==0) {
            b_target[m]=z;
            }
          else {
            b_target[m]=m_target;
            }
          }
        }
      for(j=0;j<g_target.ny;j++) {
        int m1=g_target.nx*j;
        int m2=g_target.nx*j+g_target.nx-1;
        b_target[m2]=b_target[m1];
        }
      break;
    default:
      break;
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_export(const grid_t & g_source, float *b_source, float m_source, grid_t g_target, float *b_target, float m_target, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=map_export_template(g_source, b_source, m_source, g_target, b_target, m_target, mode);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_export(const grid_t & g_source, double *b_source, double m_source,grid_t g_target, double *b_target, double m_target, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=map_export_template(g_source, b_source, m_source, g_target, b_target, m_target, mode);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_export(const grid_t & g_source, complex<float> *b_source, complex<float> m_source,grid_t g_target, complex<float> *b_target, complex<float> m_target, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=map_export_template(g_source, b_source, m_source, g_target, b_target, m_target, mode);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_export(const grid_t & g_source, complex<double> *b_source, complex<double> m_source,grid_t g_target, complex<double> *b_target, complex<double> m_target, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=map_export_template(g_source, b_source, m_source, g_target, b_target, m_target, mode);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_subgrid(grid_t & g_source, grid_t *g_target, int incx, int incy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  grid sub-sampling
  
  to be further verified and fixed
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  int i,j,m,n;
  status=0;
  
  *g_target=g_source;
  
  switch(g_source.modeH) {
    case 0:
      g_target->nx=g_source.nx/incx;
      g_target->ny=g_source.ny/incy;
      
      g_target->dx=g_source.dx*incx;
      g_target->dy=g_source.dy*incy;
      break;
    case 1:
      g_target->nx=g_source.nx/incx+1;
      g_target->ny=g_source.ny/incy+1;
      
      g_target->x=new double[g_target->nx];
      g_target->y=new double[g_target->ny];
      
      for(j=0;j<g_target->ny;j++) {
        n=j*incy;
        g_target->y[j]=g_source.y[n];
        }
      for(i=0;i<g_target->nx;i++) {
        n=i*incx;
        g_target->x[i]=g_source.x[n];
        }
      g_target->dx=g_source.dx*incx;
      g_target->dy=g_source.dy*incy;
      break;
    case 2:
      g_target->nx=g_source.nx/incx+1;
      g_target->ny=g_source.ny/incy+1;
      
      g_target->x=new double[g_target->nx*g_target->ny];
      g_target->y=new double[g_target->nx*g_target->ny];
      
      for(j=0;j<g_target->ny;j++) {
        for(i=0;i<g_target->nx;i++) {
          m=g_target->nx*j+i;
          n=g_source.nx*j*incy+i*incx;
          g_target->x[m]=g_source.x[n];
          g_target->y[m]=g_source.y[n];
          }
        }
      g_target->dx=g_source.dx*incx;
      g_target->dy=g_source.dy*incy;
      break;
    }
  
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_supgrid(grid_t & g_source, grid_t *g_target, int incx, int incy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  grid super-sampling
  
  to be further verified and fixed
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  int i,j,m,n;
  status=0;
  
  *g_target=g_source;
     
  g_target->nx=(g_source.nx-1)*incx+1;
  g_target->ny=(g_source.ny-1)*incy+1;
  
  g_target->dx=g_source.dx/(double) incx;
  g_target->dy=g_source.dy/(double) incy;
 
  switch(g_source.modeH) {
    case 0:
      break;
      
    case 1:
     
      g_target->x=new double[g_target->nx];
      g_target->y=new double[g_target->ny];
      
      for(j=0;j<g_target->ny;j++) {
        n=j*incy;
        g_target->y[j]=g_source.y[n];
        }
      for(i=0;i<g_target->nx;i++) {
        n=i*incx;
        g_target->x[i]=g_source.x[n];
        }
      break;
      
    case 2:
      
      g_target->x=new double[g_target->nx*g_target->ny];
      g_target->y=new double[g_target->nx*g_target->ny];
      
      for(j=0;j<g_target->ny;j++) {
        for(i=0;i<g_target->nx;i++) {
          m=g_target->nx*j+i;
          n=g_source.nx*j*incy+i*incx;
          g_target->x[m]=g_source.x[n];
          g_target->y[m]=g_source.y[n];
          }
        }
      break;
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_remap_template (grid_t & g_source, T *b_source, T m_source, grid_t *g_target, T **b_target, T *m_target, int incr, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int incx=incr, incy=incr;
  
  status=map_subgrid(g_source, g_target, incx, incy);
  
  *m_target=m_source;
  *b_target=new T[g_target->nx*g_target->ny];
  status=map_export_template( g_source, b_source, m_source, *g_target, *b_target, *m_target,mode);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_remap (grid_t & g_source, complex<float> *b_source, complex<float> m_source, grid_t *g_target, complex<float> **b_target, complex<float> *m_target, int incr,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=map_remap_template (g_source, b_source, m_source, g_target, b_target, m_target,incr, mode);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_remap (grid_t & g_source, float *b_source, float m_source, grid_t *g_target, float **b_target, float *m_target, int incr,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=map_remap_template (g_source, b_source, m_source, g_target, b_target, m_target,incr, mode);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  statistic_c_t cget_geostatistics(const grid_t & grid,fcomplex *h, fcomplex mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k,count,status;
  double variance;
  complex<double> mean,z;
  statistic_c_t s;
  dcomplex csum,carea;
  double   dsum,darea;
  float *buffer=NULL;
  int n;

  count=0;
  variance=0;

  status=map_integraleC(grid, h, mask, &csum,&carea);

  mean=csum/carea;

  exitIfNull(
    buffer=(float *) malloc(grid.nx*grid.ny*sizeof(float))
    );

  n=grid.nx*grid.ny;

  for(k=0;k<n;k++) {
    if(h[k] != mask) {
      buffer[k]=norm((complex<double>) h[k]-mean)/2.;
      count++;
      }
    else
      buffer[k]=1.e+10;
    }

  status=map_integraleR(grid, buffer, 1.e+10, &dsum,&darea);
  variance=dsum/darea;

  print_stats(n, count, mean, variance, 'r');

  free(buffer);
  return(s);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_projection_forward(grid_t & grid, const char *proj4)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  range_t<double> r;
  
  status=geo_to_projection (proj4, grid.x, grid.y, grid.nx*grid.ny, &grid.projection);
  
  r.init(grid.x, grid.nx*grid.ny);
  grid.xmin=r.min;
  grid.xmax=r.max;
  
  r.init(grid.y, grid.nx*grid.ny);
  grid.ymin=r.min;
  grid.ymax=r.max;
  
  grid.dx=(grid.xmax-grid.xmin)/((double)grid.nx-1.);
  grid.dy=(grid.ymax-grid.ymin)/((double)grid.ny-1.);
  
  if(proj4!=0) grid.proj4options=strdup(proj4);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_projection_backward(grid_t & grid, const char *proj4)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  range_t<double> r;
  double x;
  
/*------------------------------------------------------------------------------
  any mode 0 or 1 grid in cartesian will need mode 2 in spherical coords */
  status=map_completegridaxis(&grid, 2, 0);
  
  status=projection_to_geo (proj4, grid.x, grid.y, grid.nx*grid.ny);

  x=grid.x[0]+180.0;
  for(int j=0;j<grid.ny;j++) {
    size_t m=grid.nx*j;
    x=degree_recale(grid.x[m],x);
    for(int i=0;i<grid.nx;i++) {
      grid.x[m+i]=degree_recale(grid.x[m+i], x);
      }
    }

//   switch(grid.modeH) {
//     case 0:
//       status=projection_to_geo (proj4, grid.x, grid.y, grid.nx*grid.ny);
//       break;
//     case 1:
//       status=projection_to_geo (proj4, grid.x, grid.y, grid.nx*grid.ny);
//       break;
//     case 2:
//       status=projection_to_geo (proj4, grid.x, grid.y, grid.nx*grid.ny);
//       break;
//     }
//   
//   switch(grid.modeH) {
//     case 0:
//       break;
//     case 1:
//       x=grid.x[0]+180.0;
//       for(int i=0;i<grid.nx;i++) {
//         grid.x[i]=degree_recale(grid.x[i], x);
//         }
//       break;
//     case 2:
//       x=grid.x[0]+180.0;
//       for(int j=0;j<grid.ny;j++) {
//         size_t m=grid.nx*j;
//         x=degree_recale(grid.x[m],x);
//         for(int i=0;i<grid.nx;i++) {
//           grid.x[m+i]=degree_recale(grid.x[m+i], x);
//           }
//         }
//       break;
//     }

  r.init(grid.x, grid.nx*grid.ny);
  grid.xmin=r.min;
  grid.xmax=r.max;
  
  r.init(grid.y, grid.nx*grid.ny);
  grid.ymin=r.min;
  grid.ymax=r.max;
  
  grid.dx=(grid.xmax-grid.xmin)/((double)grid.nx-1.);
  grid.dy=(grid.ymax-grid.ymin)/((double)grid.ny-1.);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t map_get_spherical(projPJ projection,grid_t cgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    j;
  int    status;
  grid_t sgrid;

  sgrid.dx=cgrid.dx;
  sgrid.dy=cgrid.dy;
  sgrid.nx=cgrid.nx;
  sgrid.ny=cgrid.ny;
  sgrid.nz=cgrid.nz;
  sgrid.modeH=2;

  sgrid.x=new double[sgrid.nx*sgrid.ny];
  sgrid.y=new double[sgrid.nx*sgrid.ny];

  sgrid.ymin= +1.e+35;
  sgrid.ymax= -1.e+35;
  sgrid.xmin= +1.e+35;
  sgrid.xmax= -1.e+35;

  const size_t n=sgrid.Hsize();
#pragma omp parallel for
  for(size_t m=0;m<n;m++) {
    sgrid.x[m]=cgrid.x[m];
    sgrid.y[m]=cgrid.y[m];
    }

  status=projection_to_geo (cgrid.proj4options, sgrid.x, sgrid.y, n);

  for(j=0;j<sgrid.ny;j++) {
    for(int i=0;i<sgrid.nx;i++) {
      int k=sgrid.nx*j+i;
      updatemin(&sgrid.ymin,sgrid.y[k]);
      updatemax(&sgrid.ymax,sgrid.y[k]);
      updatemin(&sgrid.xmin,sgrid.x[k]);
      updatemax(&sgrid.xmax,sgrid.x[k]);
      }
    }

  return(sgrid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

grid_t map_get_spherical(geo_t projection,grid_t cgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k;
  int    status;
  grid_t sgrid;

  sgrid.dx=cgrid.dx;
  sgrid.dy=cgrid.dy;
  sgrid.nx=cgrid.nx;
  sgrid.ny=cgrid.ny;
  sgrid.nz=cgrid.nz;
  sgrid.modeH=2;

  exitIfNull(
    sgrid.x=(double *) malloc(sgrid.nx*sgrid.ny*sizeof(double))
    );
  exitIfNull(
    sgrid.y=(double *) malloc(sgrid.nx*sgrid.ny*sizeof(double))
    );

  sgrid.ymin= +1.e+35;
  sgrid.ymax= -1.e+35;
  sgrid.xmin= +1.e+35;
  sgrid.xmax= -1.e+35;

  for(j=0;j<sgrid.ny;j++) {
    for(i=0;i<sgrid.nx;i++) {
      k=sgrid.nx*j+i;
      status=geo_mercator_inverse(projection,&(sgrid.x[k]),&(sgrid.y[k]),cgrid.x[k],cgrid.y[k]);
      updatemin(&sgrid.ymin,sgrid.y[k]);
      updatemax(&sgrid.ymax,sgrid.y[k]);
      updatemin(&sgrid.xmin,sgrid.x[k]);
      updatemax(&sgrid.xmax,sgrid.x[k]);
      }
    }
  sgrid.circular=cgrid.circular;
  sgrid.overlapped=cgrid.overlapped;
  sgrid.connex=cgrid.connex;

  return(sgrid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

grid_t map_get_cartesian(projPJ projection,grid_t sgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k/*,status*/;
  grid_t cgrid;

  cgrid.dx=sgrid.dx;
  cgrid.dy=sgrid.dy;
  cgrid.nx=sgrid.nx;
  cgrid.ny=sgrid.ny;
  cgrid.nz=sgrid.nz;
  cgrid.modeH=2;

  cgrid.x=new double[cgrid.nx*cgrid.ny];
  cgrid.y=new double[cgrid.nx*cgrid.ny];

  cgrid.ymin= +1.e+35;
  cgrid.ymax= -1.e+35;
  cgrid.xmin= +1.e+35;
  cgrid.xmax= -1.e+35;
  
  
//   int n=sgrid.Hsize();
// #pragma omp parallel for
//   for(size_t m=0;m<n;m++) {
//     double x,y;
//     sgrid.xy(i,j,x,y);
//     sgrid.x[m]=cgrid.x[m];
//     sgrid.y[m]=cgrid.y[m];
//     }
//   status=geo_to_projection (projection, cgrid.x, cgrid.y, n);

  for(j=0;j<cgrid.ny;j++) {
    for(i=0;i<cgrid.nx;i++) {
      k=cgrid.nx*j+i;
      double x,y;
      sgrid.xy(i,j,x,y);
      geo_to_projection(projection,y,x,&(cgrid.x[k]),&(cgrid.y[k]));
      updatemin(&cgrid.ymin,cgrid.y[k]);
      updatemax(&cgrid.ymax,cgrid.y[k]);
      updatemin(&cgrid.xmin,cgrid.x[k]);
      updatemax(&cgrid.xmax,cgrid.x[k]);
      }
    }

  return(cgrid);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

grid_t map_get_cartesian(geo_t projection,grid_t sgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k;
  int    status;
  grid_t cgrid;

  cgrid.dx=sgrid.dx;
  cgrid.dy=sgrid.dy;
  cgrid.nx=sgrid.nx;
  cgrid.ny=sgrid.ny;
  cgrid.nz=sgrid.nz;
  cgrid.modeH=2;

  cgrid.x=new double[cgrid.nx*cgrid.ny];
  cgrid.y=new double[cgrid.nx*cgrid.ny];

  cgrid.ymin= +1.e+35;
  cgrid.ymax= -1.e+35;
  cgrid.xmin= +1.e+35;
  cgrid.xmax= -1.e+35;

  for(j=0;j<cgrid.ny;j++) {
    for(i=0;i<cgrid.nx;i++) {
      k=cgrid.nx*j+i;
      status=geo_mercator_directe(projection,sgrid.x[k],sgrid.y[k],&(cgrid.x[k]),&(cgrid.y[k]));
      updatemin(&cgrid.ymin,cgrid.y[k]);
      updatemax(&cgrid.ymax,cgrid.y[k]);
      updatemin(&cgrid.xmin,cgrid.x[k]);
      updatemax(&cgrid.xmax,cgrid.x[k]);
      }
    }

  return(cgrid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t map_Rgrid(frame_t frame, double dx, double dy, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k;
  grid_t grid;

  grid.dx=dx;
  grid.dy=dy;

  grid.nx=(int) floor((frame.xmax-frame.xmin)/dx+0.5)+1;
  grid.ny=(int) floor((frame.ymax-frame.ymin)/dy+0.5)+1;

  grid.xmin=frame.xmin;
  grid.xmax=frame.xmax;
  grid.ymin=frame.ymin;
  grid.ymax=frame.ymax;

  grid.modeH=mode;

  switch(grid.modeH) {
    case 2:
      grid.x=new double[grid.nx*grid.ny];
      grid.y=new double[grid.nx*grid.ny];
      for(j=0;j<grid.ny;j++) {
        for(i=0;i<grid.nx;i++) {
          k=grid.nx*j+i;
          grid.y[k]=grid.ymin+(double) j*dy;
          grid.x[k]=grid.xmin+(double) i*dx;
          }
        }
    }
  return(grid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t map_Rgrid(frame_t frame, size_t nx, size_t ny, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  grid_t grid;

  grid.dx=(frame.xmax-frame.xmin)/(nx-1);
  grid.dy=(frame.ymax-frame.ymin)/(ny-1);

  grid.nx=nx;
  grid.ny=ny;

  grid.xmin=frame.xmin;
  grid.xmax=frame.xmax;
  grid.ymin=frame.ymin;
  grid.ymax=frame.ymax;

  grid.modeH=0;

  return(grid);
}
