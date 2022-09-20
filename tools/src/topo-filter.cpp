

/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

// #define MAIN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "rutin.h"
#include "geo.h"
#include "polygones.h"
#include "netcdf-proto.h"
#include "functions.h"
#include "grd.h"
#include "filter.h"
#include "map.h"
#include "fe-proto.h"

#include "bmg.h"
#include "topo.h"

#define XYZ 0
#define YXZ 1
#include "map.def"



// c ***********************************************************	
// 	subroutine medfil (z, m, n, irad)
// c
// c two-dimensional median filter 
// c running window is disk-shaped if mesh is square
// c mirror images are assumed beyond boundaries and corners
// c
// c input/output: z(m,n) -- column size= maxx (parameter)
// c input: irad, disk radius in gridpoints
// c
// c de mey '88
// c
// 	parameter (maxx=301, mz1= -maxx+2, mz2= 2*maxx-1)
// 	parameter (maxs= 9*maxx*maxx-12*maxx+4)
// 	parameter (eps= 1e-6)
// c
// 	real z(m,n), zt(mz1:mz2,mz1:mz2), zs(maxs)
// c
// 	if (irad.gt.max(m,n)) return
// 	radius= irad-1
// 	if (m.gt.maxx.or.n.gt.maxx) return
// c
// c make test matrix
// c
// 	m1= -m+2
// 	m2= 2*m-1
// 	n1= -n+2
// 	n2= 2*n-1
// 	do 10 i= 1, m
// 	do 10 j= 1, n
// 10	  zt(i,j)= z(i,j)
// 	do 20 i= 0, m1, -1
// 	do 20 j= 1, n
// 20	  zt(i,j)= z(-i+2,j)
// 	do 30 i= m2, m+1, -1
// 	do 30 j= 1, n
// 30	  zt(i,j)= z(m2-i+1,j)
// 	do 40 i= m1, m2
// 	do 40 j= 0, n1, -1
// 40	  zt(i,j)= zt(i,-j+2)
// 	do 50 i= m1, m2
// 	do 50 j= n2, n+1, -1
// 50	  zt(i,j)= zt(i,n2-j+1)
// c
// c main loop
// c the filtered field is written back onto z
// c
// 	do 70 i= 1, m
// 	do 70 j= 1, n
// 	  ns= 0
// 	  do 60 ii= i-irad+1, i+irad-1
// 	  do 60 jj= j-irad+1, j+irad-1
// 	    dist= sqrt(float(ii-i)**2+float(jj-j)**2)+eps
// 	    if (dist.gt.radius) go to 60
// 	    ns= ns+1
// 	    zs(ns)= zt(ii,jj)
// 60	    continue
// 	  if (ns.eq.0) go to 70
// 	  call shell(zs,ns)
// 	  median= ns/2
// 	  z(i,j)= zs(median)
// 70	  continue
// c
// 	return
// 	end



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int  median_filter_spherical(grid_t & grid, T *buffer, T mask, double radius, int i, int j, T & out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int width;
  double radius2=radius*radius;
  int m=j*grid.nx+i;
  
//   width=radius/grid.dx/110.+1;
  width=radius/grid.dx+1;
  
  int n=j*grid.nx+i;
  double t,p;
  grid.xy(i,j,t,p);
  
  float *values=new float[(2*width+1)*(2*width+1)];
  int nvalues=0;
  
  range_t<int> i_range,j_range;
  i_range.init(i,width,0,grid.nx-1);
  j_range.init(j,width,0,grid.ny-1);
  
  for(int jj=j_range.min;jj<j_range.max;jj++) {
    for(int ii=i_range.min;ii<i_range.max;ii++) {
      int nn=jj*grid.nx+ii;
      if(buffer[nn]==mask) continue;
      double tt,pp;
      grid.xy(ii,jj,tt,pp);
//       double d=geo_haversin_km(t,p,tt,pp);
      double d=(t-tt)*(t-tt)+(p-pp)*(p-pp);
//      if(d<radius*110.) {
      if(d<radius2) {
        values[nvalues]=buffer[nn];
        nvalues++;
        }
      }
    }

  if(nvalues==0) {
    out=buffer[m];
    }
  else if(nvalues==1) {
    out=values[0];
    }
  else {
    median_shell(values,nvalues);
    int median= nvalues/2;
    out=values[median];
    }
    
  delete[] values;


  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int loess_filter_init(const grid_t & grid, float scale, float* & weight, int & nx, int & ny)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j;

  float lx,ly,q,dum;

  nx=int(NINT(scale/grid.dx+1)+1);
  ny=int(NINT(scale/grid.dy+1)+1);

  nx=MIN(nx,grid.nx);
  ny=MIN(ny,grid.ny);

  weight=new float[nx*ny];

  for(j=0;j<ny;j++)
    for (i=0;i<nx;i++)
      weight[j*nx+i]=0.;

  for(j=0;j<ny;j++) {
    for(i=0;i<nx;i++) {
      lx=i*grid.dx/scale;
      ly=j*grid.dy/scale;
      q=lx*lx+ly*ly;
      if (q <= 1.) {
        dum=1-q*q*q;
        weight[j*nx+i]=dum*dum*dum;
        }
      }
    }

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename TYPE> int loess_filter_template(const grid_t & grid, TYPE *buf, TYPE mask, float *weight, int nx, int ny, int k, int l, TYPE & out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j, offset_j;
  int imin,imax,jmin,jmax;
  float w,r;
  float sum,z;
  float tmp,value;

  tmp=0.;
  sum=0.;

  imin=MAX(0,k-nx+1);
  imax=MIN(grid.nx,k+nx);
  jmin=MAX(0,l-ny+1);
  jmax=MIN(grid.ny,l+ny);

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

  out=value;

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int loess_filter(const grid_t & grid, float *buf, float mask, float *weight, int nx, int ny, int k, int l, float & out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=loess_filter_template(grid, buf, mask, weight, nx, ny, k, l, out);
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int loess_filter(const grid_t & grid, float *buf, float mask, float scale, int k, int l, float & out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int nx,ny;
  float *weight;
  
  status=loess_filter_init(grid, scale, weight, nx,  ny);
  
  status=loess_filter_template(grid, buf, mask, weight, nx, ny, k, l, out);
  
  delete[] weight;
  return(status);
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// template <typename TYPE> int window_filter(grid_t grid, TYPE *buf, TYPE mask, float scale, int k, int l, TYPE & out)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int i,j;
//   int imin,imax,jmin,jmax;
//   float w,r,s,xx,yy;
//   float lx,ly,q,dum,sum,z;
//   int nx,ny,offset_j,offset_l;
//   float *weight,tmp,value;
// 
//   nx=int(NINT(scale/grid.dx+1)+1);
//   ny=int(NINT(scale/grid.dy+1)+1);
// 
//   nx=MIN(nx,grid.nx);
//   ny=MIN(ny,grid.ny);
// 
//   weight=new float[nx*ny];
// 
//   for(j=0;j<ny;j++)
//     for (i=0;i<nx;i++)
//       weight[j*nx+i]=0.;
// 
//   for(j=0;j<ny;j++)
//     for(i=0;i<nx;i++) {
//       lx=i*grid.dx/scale;
//       ly=j*grid.dy/scale;
//       q=lx*lx+ly*ly;
//       if (q <= 1.) {
//         dum=1-q*q*q;
//         weight[j*nx+i]=dum*dum*dum;
//         }
//       }
// 
//   tmp=0.;
//   sum=0.;
// 
//   imin=MAX(0,k-nx+1);
//   imax=MIN(grid.nx,k+nx);
//   jmin=MAX(0,l-ny+1);
//   jmax=MIN(grid.ny,l+ny);
// 
//   for(j=jmin;j<jmax;j++) {
//     offset_j=j*grid.nx;
//     for(i=imin;i<imax;i++) {
//       z=buf[offset_j+i];
//       if(z!=mask) {
//         w=weight[abs(j-l)*nx+abs(i-k)];
//         sum+=w;
//         tmp+=z*w;
//         }
//       }
//     }
// 
//   if(sum!=0.) {
//     value=tmp/sum;
//     }
//   else {
//     value=mask;
//     }
// 
//   out=value;
// 
//   delete[] weight;
//   return(0);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int window_filter(grid_t grid, float *buf, float mask, float scale, int k, int l, float & out)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   
//   status=loess_filter(grid, buf, mask, scale, k, l, out);
//   
//   return(status);
// }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename TYPE> int topo_smooth_template(grid_t grid, TYPE *topo, TYPE topomask, float maxscale, range_t<float> range, const char *poly, const char *trusted)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,l,m,count,status;
  int   imin,imax,jmin,jmax;
  float *buffer=0,*radius=0,mask,z,d;
  signed char   *selected=0;
  double x,y;
  grid_t chk_grid;
  float *chk_topo,chk_mask;
  bool debug=false;
  
  status=topo_loadfield("/home/data/topography/regional/argentina/meshing/mapxyz-unstructured.grd", &chk_grid, &chk_topo, &chk_mask, debug);

  buffer=new TYPE[grid.nx*grid.ny];
  radius=new float[grid.nx*grid.ny];

  selected=new signed char[grid.nx*grid.ny];
  for(int j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      selected[m]=1;
      buffer[m]=topo[m];
      radius[m]=0.0;
      }
    }

  if(poly!=0) {
    count=topo_PolygonsSelection(grid, poly, selected, true, 0);
    if(count<0) return(-1);
    }
  if (range.width()!=0.) {
    status=topo_RangeSelection(grid, topo, topomask, range, selected, 0);
    }

/* *------------------------------------------------------------------------------
  */
  for(int j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      if(selected[m]==0) continue;
      double t,p;
      grid.xy(i,j,t,p);
      float zz;
      t=map_recale(chk_grid,t);
      status=map_interpolation(chk_grid, chk_topo, chk_mask,t,p,&zz);
      if(zz==chk_mask) continue;
      if(zz<-50.) {
        float ratio=fabs(topo[m]/zz);
        if((ratio-0.8)*(ratio-1.2)>0) topo[m]=zz;
        }
      }
    }

/* *------------------------------------------------------------------------------
  to smooth selected points with prescribed lengthscale*/
  for(int j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      if (topo[m]<range.min) selected[m]=0;
      if(selected[m]==0) continue;
      if (topo[m]!=topomask) {
        radius[m]=maxscale;
        }
      }
    }

  int nprocs=initialize_OPENMP(-1);

/* *------------------------------------------------------------------------------
  to smooth surrounding points with decreasing lengthscale*/
  if(poly!=0) 
  for(int j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      if (radius[m]==maxscale) {
        imin=MAX(i-maxscale/grid.dx,0);
        imax=MIN(i+maxscale/grid.dx+1,grid.nx);
        jmin=MAX(j-maxscale/grid.dy,0);
        jmax=MIN(j+maxscale/grid.dy+1,grid.ny);
#pragma omp parallel for if(nprocs>1)  
        for(int jj=jmin;jj<jmax;jj++) {
          for(int ii=imin;ii<imax;ii++) {
            int n=jj*grid.nx+ii;
            double dx=(ii-i)*grid.dx;
            double dy=(jj-j)*grid.dy;
            double d=MIN(1.,1.-sqrt(dx*dx+dy*dy)/maxscale);
            d=MAX(d,0.0);
            radius[n]=MAX(radius[n],maxscale*d);
            }
          }
        }
      }
    }

//#pragma omp parallel for private(status) if(nprocs>1)  
  for(int j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      size_t m=j*grid.nx+i;
      if(selected[m]==0) continue;
      if (radius[m]!=0) {
//        status=loess_filter(grid, topo, topomask, radius[m], i, j, buffer[m]);
        status=median_filter_spherical(grid, topo, topomask, radius[m], i, j, buffer[m]);
        }
      }
    }

  for(int j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      if(selected[m]==0) continue;
      topo[m]=buffer[m];
      }
    }

  delete[] buffer;
  delete[] radius;
  delete[] selected;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_smooth_spherical(grid_t grid, float *topo, float topomask, float maxscale, range_t<float> range, const char *poly, const char *trusted)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=topo_smooth_template(grid, topo, topomask, maxscale, range, poly, trusted);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename TYPE> int topo_smooth_template(grid_t grid, TYPE *topo, TYPE topomask, float maxscale, signed char*selected, const char *trusted, int filter)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,l,m,count,status;
  int   imin,imax,jmin,jmax;
  float *buffer=0,*radius=0,mask,z,d;
  double x,y;
  grid_t chk_grid;
  float *chk_topo,chk_mask;
  bool debug=false;
  
  buffer=new TYPE[grid.nx*grid.ny];
  radius=new float[grid.nx*grid.ny];

/*------------------------------------------------------------------------------
  replace outliers with trusted bathymetry*/
//   status=topo_loadfield("/home/data/topography/regional/argentina/meshing/mapxyz-unstructured.grd", &chk_grid, &chk_topo, &chk_mask, debug);
//   for(int j=0;j<grid.ny;j++) {
//     for(int i=0;i<grid.nx;i++) {
//       m=j*grid.nx+i;
//       if(selected[m]==0) continue;
//       double t,p;
//       grid.xy(i,j,t,p);
//       float zz;
//       t=map_recale(chk_grid,t);
//       status=map_interpolation(chk_grid, chk_topo, chk_mask,t,p,&zz);
//       if(zz==chk_mask) continue;
//       if(zz<-50.) {
//         float ratio=fabs(topo[m]/zz);
//         if((ratio-0.8)*(ratio-1.2)>0) topo[m]=zz;
//         }
//       }
//     }

/*------------------------------------------------------------------------------
  to smooth selected points with prescribed lengthscale*/
  for(int j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      buffer[m]=topo[m];
      radius[m]=0.0;
      if(selected[m]==0) continue;
      if (topo[m]!=topomask) {
        radius[m]=maxscale;
        }
      }
    }

  int nprocs=initialize_OPENMP(-1);

/*------------------------------------------------------------------------------
  to smooth surrounding points with decreasing lengthscale*/
//   if(poly!=0) 
//   for(int j=0;j<grid.ny;j++) {
//     for(int i=0;i<grid.nx;i++) {
//       m=j*grid.nx+i;
//       if (radius[m]==maxscale) {
//         imin=MAX(i-maxscale/grid.dx,0);
//         imax=MIN(i+maxscale/grid.dx+1,grid.nx);
//         jmin=MAX(j-maxscale/grid.dy,0);
//         jmax=MIN(j+maxscale/grid.dy+1,grid.ny);
// #pragma omp parallel for if(nprocs>1)  
//         for(int jj=jmin;jj<jmax;jj++) {
//           for(int ii=imin;ii<imax;ii++) {
//             int n=jj*grid.nx+ii;
//             double dx=(ii-i)*grid.dx;
//             double dy=(jj-j)*grid.dy;
//             double d=MIN(1.,1.-sqrt(dx*dx+dy*dy)/maxscale);
//             d=MAX(d,0.0);
//             radius[n]=MAX(radius[n],maxscale*d);
//             }
//           }
//         }
//       }
//     }

  if(filter==LOESS) {
    printf("%s : use LOESS filter\n", __func__);
    }
  else {
    printf("%s : use MEDIAN filter\n", __func__);
    }

// #pragma omp parallel for private(status) if(nprocs>1)  
  for(int j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      size_t m=j*grid.nx+i;
      if(selected[m]==0) continue;
      if (radius[m]!=0) {
        if(filter==LOESS) {
          status=loess_filter(grid, topo, topomask, radius[m], i, j, buffer[m]);
          }
        else {
          status=median_filter_spherical(grid, topo, topomask, radius[m], i, j, buffer[m]);
          }
        }
      }
    }

  for(int j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      if(selected[m]==0) continue;
      topo[m]=buffer[m];
      }
    }

  delete[] buffer;
  delete[] radius;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_smooth_cartesian(const grid_t & grid, float *topo, float topomask, float maxscale, range_t<float> range, bool keep_masked, vector<plg_t> & polygons, const char *trusted, int filter)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  size_t count;
  grid_t cgrid;
  double dx;
  signed char*selected=0;
  float *ctopo;
  bool check_polygons=true;
  vector<plg_t> limits;
    
  if(polygons.size()==0) {
    frame_t frame(grid.xmin,grid.xmax,grid.ymin,grid.ymax);
    plg_t p(frame, PLG_SPHERICAL);
    polygons.push_back(p);
    check_polygons=false;
    }
    
  dx=grid.dx*110000.;
  status=fe_defgrid(&cgrid, polygons, limits, dx);
  cgrid.modeH=0;
  
  ctopo=new float[cgrid.Hsize()];

  int nprocs=initialize_OPENMP(-1);

#pragma omp parallel for private(status) if(nprocs>1)  
  for(int m=0;m<cgrid.Hsize();m++) {
    double x, y, t, p;
    cgrid.xy(m,x,y);
    projection_to_geo(cgrid.projection, &p, &t, x, y);
    status=map_interpolation(grid, topo, topomask, t, p, &ctopo[m]);
    }
    
  selected=new signed char[cgrid.Hsize()];
  for(int j=0;j<cgrid.ny;j++) {
    for(int i=0;i<cgrid.nx;i++) {
      size_t m=j*cgrid.nx+i;
      selected[m]=1;
      if(keep_masked and ctopo[m]==topomask) selected[m]=0;
      }
    }

  if(check_polygons) {
    count=topo_PolygonsSelection(cgrid, polygons, selected, false, 0);
    if(count<0) return(-1);
    }
    
  if (range.min!=NAN and range.max!=NAN and range.min!=range.max) {
    status=topo_RangeSelection(cgrid, ctopo, topomask, range, selected, 0);
    }
  
  count=occurence(topomask, ctopo, cgrid.Hsize());
  status=topo_smooth_template(cgrid, ctopo, topomask, maxscale*110000, selected, trusted, filter);
  
#pragma omp parallel for private(status) if(nprocs>1)  
  for(int m=0;m<grid.Hsize();m++) {
    if(keep_masked and topo[m]==topomask) continue;
//    if(selected[m]==0) continue;
    double x, y, t, p;
    float z;
    grid.xy(m,t,p);
    geo_to_projection(cgrid.projection, p, t, &x, &y);
    status=map_interpolation(cgrid, ctopo, topomask, x, y, &z);
    if(status==0) topo[m]=z;
    }

  delete[] selected;

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_smooth_cartesian(const grid_t & grid, float *topo, float topomask, float maxscale, range_t<float> range, bool keep_masked, const char *polygonfile, const char *trusted, int filter)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  vector<plg_t> polygons;
  
  if(polygonfile!=0) {
    status=plg_load(polygonfile, PLG_FORMAT_UNKNOWN, polygons);
    if(status!=0) return(-1);
    }
  
  status=topo_smooth_cartesian(grid, topo, topomask, maxscale, range, keep_masked, polygons, trusted, filter);
  
  return(status);
}

