
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
#include <string.h>
#include <stdarg.h>

#include "tools-structures.h"

#include "map.h"
#include "grd.h"
#include "ascii.h"
#include "netcdf-proto.h"
#include "poc-time.h"
#include "map.def"
// #include "sts/tugo-commons.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**
\brief Purpose : create vorticity grid from z grid

Note:

    zgrid to be TRUE zgrid (no funny additional row nor column)
    
    inside :v-points are barycenters of z-points
    
    limits : extrapolation needed
    
    
      v           v           v
    
            z           z
    
      v           v           v
             
            z           z
    
      v           v           v


@param zgrid 	input z_grid	
@return grid_t vorticity grid
*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t map_vgrid_2(const grid_t & zgrid, float *ztopo, float topomask, const char *output)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  construct vorticity grid from tracer grid (or vertex grd from cell center grid)
  
  needs mode 2 coordinates
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int    i,j,k,l,m,n;
  int    status;
  double x,y;
  float *vtopo;
  pocgrd_t ncgrid;
  grid_t grid;
//  const char *output="vorticity-grid.nc";
  double mask;
  cdfvar_t variable;

  grid.nx=zgrid.nx+1;
  grid.ny=zgrid.ny+1;
  grid.nz=zgrid.nz+1;

  grid.x=new double[grid.nx*grid.ny];
  grid.y=new double[grid.nx*grid.ny];
  if(zgrid.z!=0) grid.z=new double[grid.nx*grid.ny*grid.nz];

  grid.zmask=zgrid.zmask;
  mask=grid.zmask;

  if(zgrid.z!=0)
  for(k=0;k<grid.nz;k++) {
    for(j=0;j<grid.ny;j++) {
      for(i=0;i<grid.nx;i++) {
        n=grid.nx*grid.ny*k+grid.nx*j+i;
        grid.z[n]=grid.zmask;
        }
      }
    }

/* *-----------------------------------------------------------------------------
  interpolate mid-points */
  for(l=1;l<grid.ny-1;l++) {
    for(k=1;k<grid.nx-1;k++) {
      x= 0;
      y= 0;
      i=k-1;
      j=l-1;
      m=zgrid.nx*j+i;
      x+= zgrid.x[m];
      y+= zgrid.y[m];
      i=k;
      j=l-1;
      m=zgrid.nx*j+i;
      x+= zgrid.x[m];
      y+= zgrid.y[m];
      i=k;
      j=l;
      m=zgrid.nx*j+i;
      x+= zgrid.x[m];
      y+= zgrid.y[m];
      i=k-1;
      j=l;
      m=zgrid.nx*j+i;
      x+= zgrid.x[m];
      y+= zgrid.y[m];
      n=grid.nx*l+k;
      grid.x[n]=x/4.0;
      grid.y[n]=y/4.0;
      }
    }

/* *-----------------------------------------------------------------------------
  linear extrapolation */
  for(l=1;l<grid.ny-1;l++) {
    k=0;
    n=grid.nx*l+k;
    grid.x[n]=2.*grid.x[grid.nx*l+k+1]-grid.x[grid.nx*l+k+2];
    grid.y[n]=2.*grid.y[grid.nx*l+k+1]-grid.y[grid.nx*l+k+2];
    k=grid.nx-1;
    n=grid.nx*l+k;
    grid.x[n]=2.*grid.x[grid.nx*l+k-1]-grid.x[grid.nx*l+k-2];
    grid.y[n]=2.*grid.y[grid.nx*l+k-1]-grid.y[grid.nx*l+k-2];
    }

  for(k=1;k<grid.nx-1;k++) {
    l=0;
    n=grid.nx*l+k;
    grid.x[n]=2.*grid.x[grid.nx*(l+1)+k]-grid.x[grid.nx*(l+2)+k];
    grid.y[n]=2.*grid.y[grid.nx*(l+1)+k]-grid.y[grid.nx*(l+2)+k];
    l=grid.ny-1;
    n=grid.nx*l+k;
    grid.x[n]=2.*grid.x[grid.nx*(l-1)+k]-grid.x[grid.nx*(l-2)+k];
    grid.y[n]=2.*grid.y[grid.nx*(l-1)+k]-grid.y[grid.nx*(l-2)+k];
    }

  k=0;
  l=0;
  n=grid.nx*l+k;
  grid.x[n]=2.*grid.x[grid.nx*(l+1)+k+1]-grid.x[grid.nx*(l+2)+k+2];
  grid.y[n]=2.*grid.y[grid.nx*(l+1)+k+1]-grid.y[grid.nx*(l+2)+k+2];
  
  k=grid.nx-1;
  l=0;
  n=grid.nx*l+k;
  grid.x[n]=2.*grid.x[grid.nx*(l+1)+k-1]-grid.x[grid.nx*(l+2)+k-2];
  grid.y[n]=2.*grid.y[grid.nx*(l+1)+k-1]-grid.y[grid.nx*(l+2)+k-2];
    
  k=grid.nx-1;
  l=grid.ny-1;
  n=grid.nx*l+k;
  grid.x[n]=2.*grid.x[grid.nx*(l-1)+k-1]-grid.x[grid.nx*(l-2)+k-2];
  grid.y[n]=2.*grid.y[grid.nx*(l-1)+k-1]-grid.y[grid.nx*(l-2)+k-2];
    
  k=0;
  l=grid.ny-1;
  n=grid.nx*l+k;
  grid.x[n]=2.*grid.x[grid.nx*(l-1)+k+1]-grid.x[grid.nx*(l-2)+k+2];
  grid.y[n]=2.*grid.y[grid.nx*(l-1)+k+1]-grid.y[grid.nx*(l-2)+k+2];
    
  grid.modeH=2;
  
  grid.xmax=grid.x[0];
  grid.xmin=grid.x[0];
  grid.ymax=grid.y[0];
  grid.ymin=grid.y[0];
  for(int i=0; i<(grid.nx*grid.ny);i++){
    if (grid.x[i] > grid.xmax) grid.xmax=grid.x[i];
    if (grid.x[i] < grid.xmin) grid.xmin=grid.x[i];
    if (grid.y[i] > grid.ymax) grid.ymax=grid.y[i];
    if (grid.y[i] < grid.ymin) grid.ymin=grid.y[i];
    }
    
  grid.projection=zgrid.projection;
  if(zgrid.proj4options!=0)
    grid.proj4options=poc_strdup(zgrid.proj4options);
   
  if(ztopo!=0 && output!=0) {
    vtopo=map_extrapolate_z2v(zgrid, grid, ztopo, topomask);
    status= poc_createfile(output);
    status=poc_sphericalgrid_xy(output,"",grid,&ncgrid);
    mask=1.e+35;
    poc_standardvariable_xy(& variable,"dummy",mask,"m",1., 0.,"dummy","dummy","dummy",ncgrid);
    status=create_ncvariable(output, &variable);
    status = poc_write_xy(output, grid, variable.id, vtopo);
    variable.destroy();
    }
      
  return(grid);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t map_vgrid_1(const grid_t & zgrid, float *ztopo, float topomask, const char *output)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  construct vorticity grid from tracer grid (or vertex grd from cell center grid)
  
  needs mode 1 coordinates
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int    i,j,k,l,m,n;
  int    status;
  grid_t grid;
  
  grid.nx=zgrid.nx+1;
  grid.ny=zgrid.ny+1;
  grid.nz=zgrid.nz+1;

  grid.x=new double[grid.nx];
  grid.y=new double[grid.ny];
//   if(zgrid.z!=0) grid.z=new double[grid.nx*grid.ny*grid.nz];

//   grid.zmask=zgrid.zmask;
//   mask=grid.zmask;

//   if(zgrid.z!=0)
//   for(k=0;k<grid.nz;k++) {
//     for(j=0;j<grid.ny;j++) {
//       for(i=0;i<grid.nx;i++) {
//         n=grid.nx*grid.ny*k+grid.nx*j+i;
//         grid.z[n]=grid.zmask;
//         }
//       }
//     }
  
  for(k=1;k<grid.nx-1;k++) {
    grid.x[k]=0.5*(zgrid.x[k-1]+zgrid.x[k]);
    }
  k=0;
  grid.x[k]=grid.x[k+1]-(grid.x[k+2]-grid.x[k+1]);
  k=grid.nx-1;
  grid.x[k]=grid.x[k-1]+(grid.x[k-1]-grid.x[k-2]);
  
  for(k=1;k<grid.ny-1;k++) {
    grid.y[k]=0.5*(zgrid.y[k-1]+zgrid.y[k]);
    }
  k=0;
  grid.y[k]=grid.y[k+1]-(grid.y[k+2]-grid.y[k+1]);
  k=grid.ny-1;
  grid.y[k]=grid.y[k-1]+(grid.y[k-1]-grid.y[k-2]);
  
  grid.modeH=1;
  
  return(grid);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t map_vgrid(const grid_t & zgrid, float *ztopo, float topomask, const char *output)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  construct vorticity grid from tracer grid
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  grid_t grid;
  switch (zgrid.modeH) {
    case 0:
      TRAP_ERR_EXIT(-1,"not yet implemented\n");
      break;
    case 1:
      grid=map_vgrid_1(zgrid, ztopo, topomask, output);
      break;
    case 2:
      grid=map_vgrid_2(zgrid, ztopo, topomask, output);
      break;
    }
  return(grid);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t map_z2ugrid(const grid_t & zgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k,k1,k2;
  double dx;
  grid_t ugrid;

  ugrid.dx=zgrid.dx;
  ugrid.dy=zgrid.dy;
  ugrid.nx=zgrid.nx;
  ugrid.ny=zgrid.ny;
  ugrid.nz=zgrid.nz;
  ugrid.modeH=0;

  ugrid.x=new double[ugrid.nx*ugrid.ny];
  ugrid.y=new double[ugrid.nx*ugrid.ny];

  ugrid.ymin= +1.e+35;
  ugrid.ymax= -1.e+35;
  ugrid.xmin= +1.e+35;
  ugrid.xmax= -1.e+35;

  for(j=0;j<ugrid.ny;j++) {
    for(i=0;i<ugrid.nx;i++) {
      k=ugrid.nx*j+i;
      if(i==0) {
        k1=zgrid.nx*j+i;
        k2=zgrid.nx*j+i+1;
        }
      else {
        k1=zgrid.nx*j+i-1;
        k2=zgrid.nx*j+i;
        }
      dx=zgrid.x[k2]-zgrid.x[k1];
      ugrid.x[k]=zgrid.x[k]-0.5*dx;
      ugrid.y[k]=zgrid.y[k];
      updatemin(&ugrid.ymin,ugrid.y[k]);
      updatemax(&ugrid.ymax,ugrid.y[k]);
      updatemin(&ugrid.xmin,ugrid.x[k]);
      updatemax(&ugrid.xmax,ugrid.x[k]);
      }
    }

  return(ugrid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  resize_t map_CheckZgrid(grid_t & zgrid, char *zlandmask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,m,n;
  bool masked;
  bool remove_first_row=false,remove_first_col=false;
  bool remove_last_row=false,remove_last_col=false;
  grid_t grid=zgrid;
  resize_t resize;
  char OceanFlag=1;
  
  resize.former.push_back(grid.nx);
  resize.former.push_back(grid.ny);

/*------------------------------------------------------------------------------
  check if firts row has only mask values */  
  i=0;
  masked=true;
  for(j=0;j<zgrid.ny;j++) {
    n=zgrid.nx*j+i;
    if(zlandmask[n]==OceanFlag) {
      masked=false;
      break;
      }
    }
  if(masked) remove_first_row=true;
  
/*------------------------------------------------------------------------------
  check if last row has only mask values */  
  i=zgrid.nx-1;
  masked=true;
  for(j=0;j<zgrid.ny;j++) {
    n=zgrid.nx*j+i;
    if(zlandmask[n]==OceanFlag) {
      masked=false;
      break;
      }
    }
  if(masked) remove_last_row=true;
  
/*------------------------------------------------------------------------------
  check if first column has only mask values */  
  j=0;
  masked=true;
  for(i=0;i<zgrid.nx;i++) {
    n=zgrid.nx*j+i;
    if(zlandmask[n]==OceanFlag) {
      masked=false;
      break;
      }
    }
  if(masked) remove_first_col=true;
  
/*------------------------------------------------------------------------------
  check if last column has only mask values */  
  j=zgrid.ny-1;
  masked=true;
  for(i=0;i<zgrid.nx;i++) {
    n=zgrid.nx*j+i;
    if(zlandmask[n]==OceanFlag) {
      masked=false;
      break;
      }
    }
  if(masked) remove_last_col=true;
  
  
  if(!remove_first_col && !remove_first_row && !remove_last_col && !remove_last_row) return(resize);
  
  printf("map_CheckZgrid: remove_last_col =%d remove_last_row =%d\n",remove_last_col, remove_last_row);
  printf("map_CheckZgrid: remove_first_col=%d remove_first_row=%d\n",remove_first_col,remove_first_row);
  
  if(remove_last_col) grid.nx--;
  if(remove_last_row) grid.ny--;
  
  if(remove_first_col) {
    grid.nx--;
    resize.start.push_back(1);
    }
  else {
    resize.start.push_back(0);
    }

  if(remove_first_row) {
    grid.ny--;
    resize.start.push_back(1);
    }
  else {
    resize.start.push_back(0);
    }
  
  grid.x=new double[grid.Hsize()];
  grid.y=new double[grid.Hsize()];
    
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      n=zgrid.nx*(j+resize.start[1])+i+resize.start[0];
      m=grid.nx*j+i;
      grid.x[m]=zgrid.x[n];
      grid.y[m]=zgrid.y[n];
      }
    }
    
  grid.proj4options=strdup(zgrid.proj4options);
  grid.projection=zgrid.projection;
  
  zgrid.free();
  
  zgrid=grid;
  
//   resize.start.push_back(0);
//   resize.start.push_back(0);

  resize.end.push_back(grid.nx);
  resize.end.push_back(grid.ny);
  
  resize.neutral=false;

  return(resize);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_resize_grid(grid_t & zgrid, resize_t resize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int i,j,m,n;
  grid_t grid=zgrid;
  range_t<double> range;
  
  if(resize.neutral) return(0);
    
  grid.nx=resize.end[0]-resize.start[0]+1;
  grid.ny=resize.end[1]-resize.start[1]+1;
  
  grid.x=new double[grid.Hsize()];
  grid.y=new double[grid.Hsize()];
 
  for(j=0;j<grid.ny;j++) {
    int jj=j+resize.start[1];
    for(i=0;i<grid.nx;i++) {
      int ii=i+resize.start[0];
      n=zgrid.nx*jj+ii;
      m=grid.nx*j+i;
      grid.x[m]=zgrid.x[n];
      grid.y[m]=zgrid.y[n];
      }
    }
 
  status=map_minmax(&grid);
  
  zgrid.free();
  
  zgrid=grid;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_resize_template(grid_t & zgrid, const resize_t & resize, T * & buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  int i,j,m,n;
  grid_t grid;
  T *tmp;
  
  if(resize.neutral) return(0);
  
  grid.nx=resize.former[0];
  grid.ny=resize.former[1];
  
  tmp=new T[zgrid.Hsize()];
  
  for(j=0;j<zgrid.ny;j++) {
    int jj=j+resize.start[1];
    for(i=0;i<zgrid.nx;i++) {
      int ii=i+resize.start[0];
      n=grid.nx*jj+ii;
      m=zgrid.nx*j+i;
      tmp[m]=buffer[n];
      }
    }
  
  delete[] buffer;
  
  buffer=tmp;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_resize(grid_t & zgrid, const resize_t & resize, char * & buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=map_resize_template(zgrid, resize, buffer);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_resize(grid_t & zgrid, const resize_t & resize, float * & buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=map_resize_template(zgrid, resize, buffer);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_resize(grid_t & zgrid, const resize_t & resize, double * & buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=map_resize_template(zgrid, resize, buffer);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t map_z2vgrid(const grid_t & zgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k,k1,k2;
  double dy;
  grid_t vgrid;

  vgrid.dx=zgrid.dx;
  vgrid.dy=zgrid.dy;
  vgrid.nx=zgrid.nx;
  vgrid.ny=zgrid.ny;
  vgrid.nz=zgrid.nz;
  vgrid.modeH=0;

  vgrid.x=new double[vgrid.nx*vgrid.ny];
  vgrid.y=new double[vgrid.nx*vgrid.ny];

  vgrid.ymin= +1.e+35;
  vgrid.ymax= -1.e+35;
  vgrid.xmin= +1.e+35;
  vgrid.xmax= -1.e+35;

  for(j=0;j<vgrid.ny;j++) {
    for(i=0;i<vgrid.nx;i++) {
      k=vgrid.nx*j+i;
      if(j==0) {
        k1=zgrid.nx*j+i;
        k2=zgrid.nx*(j+1)+i;
        }
      else {
        k1=zgrid.nx*(j-1)+i;
        k2=zgrid.nx*j+i;
        }
      vgrid.x[k]=zgrid.x[k];
      dy=zgrid.y[k2]-zgrid.y[k1];
      vgrid.y[k]=zgrid.y[k]-0.5*dy;
      updatemin(&vgrid.ymin,vgrid.y[k]);
      updatemax(&vgrid.ymax,vgrid.y[k]);
      updatemin(&vgrid.xmin,vgrid.x[k]);
      updatemax(&vgrid.xmax,vgrid.x[k]);
      }
    }

  return(vgrid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t map_f2zgrid(const grid_t & fgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k;
  grid_t zgrid;

  zgrid.nx=fgrid.nx-1;
  zgrid.ny=fgrid.ny-1;
  zgrid.nz=fgrid.nz;
  zgrid.modeH=fgrid.modeH;

  zgrid.x=new double[zgrid.Hsize()];
  zgrid.y=new double[zgrid.Hsize()];

  zgrid.ymin= +1.e+35;
  zgrid.ymax= -1.e+35;
  zgrid.xmin= +1.e+35;
  zgrid.xmax= -1.e+35;

  for(j=0;j<zgrid.ny;j++) {
    for(i=0;i<zgrid.nx;i++) {
      k=zgrid.nx*j+i;
      int m1=fgrid.Hindex(i,j);
      int m2=fgrid.Hindex(i+1,j);
      int m3=fgrid.Hindex(i+1,j+1);
      int m4=fgrid.Hindex(i,j+1);
      zgrid.x[k]=0.25*(fgrid.x[m1]+fgrid.x[m2]+fgrid.x[m3]+fgrid.x[m4]);
      zgrid.y[k]=0.25*(fgrid.y[m1]+fgrid.y[m2]+fgrid.y[m3]+fgrid.y[m4]);
      updatemin(&zgrid.ymin,zgrid.y[k]);
      updatemax(&zgrid.ymax,zgrid.y[k]);
      updatemin(&zgrid.xmin,zgrid.x[k]);
      updatemax(&zgrid.xmax,zgrid.x[k]);
      }
    }
  
  if(fgrid.proj4options!=0) {
    zgrid.proj4options=strdup(fgrid.proj4options);
    zgrid.projection=fgrid.projection;
    }

  return(zgrid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t map_f2vgrid(const grid_t & fgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k;
  grid_t grid;

  grid.nx=fgrid.nx-1;
  grid.ny=fgrid.ny;
  grid.nz=fgrid.nz;
  grid.modeH=fgrid.modeH;

  grid.x=new double[grid.Hsize()];
  grid.y=new double[grid.Hsize()];

  grid.ymin= +1.e+35;
  grid.ymax= -1.e+35;
  grid.xmin= +1.e+35;
  grid.xmax= -1.e+35;

  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      k=grid.nx*j+i;
      int m1=fgrid.Hindex(i,j);
      int m2=fgrid.Hindex(i+1,j);
      grid.x[k]=0.5*(fgrid.x[m1]+fgrid.x[m2]);
      grid.y[k]=0.5*(fgrid.y[m1]+fgrid.y[m2]);
      updatemin(&grid.ymin,grid.y[k]);
      updatemax(&grid.ymax,grid.y[k]);
      updatemin(&grid.xmin,grid.x[k]);
      updatemax(&grid.xmax,grid.x[k]);
      }
    }

  return(grid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t map_f2ugrid(const grid_t & fgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k;
  grid_t grid;

  grid.nx=fgrid.nx;
  grid.ny=fgrid.ny-1;
  grid.nz=fgrid.nz;
  grid.modeH=fgrid.modeH;

  grid.x=new double[grid.Hsize()];
  grid.y=new double[grid.Hsize()];

  grid.ymin= +1.e+35;
  grid.ymax= -1.e+35;
  grid.xmin= +1.e+35;
  grid.xmax= -1.e+35;

  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      k=grid.nx*j+i;
      int m1=fgrid.Hindex(i,j);
      int m2=fgrid.Hindex(i,j+1);
      grid.x[k]=0.5*(fgrid.x[m1]+fgrid.x[m2]);
      grid.y[k]=0.5*(fgrid.y[m1]+fgrid.y[m2]);
      updatemin(&grid.ymin,grid.y[k]);
      updatemax(&grid.ymax,grid.y[k]);
      updatemin(&grid.xmin,grid.x[k]);
      updatemax(&grid.xmax,grid.x[k]);
      }
    }

  return(grid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float *map_interpolate_v2z(grid_t zgrid, grid_t vgrid, float *v_field, float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes

  Purpose : create z field from vorticity field

  Check :

  Note:

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int k,l,n,status;
  double x,y;
  size_t count;
  
  float *z_field=new float[zgrid.Hsize()];
  
  count=occurence(mask, v_field, vgrid.Hsize());
  
  for(l=0;l<zgrid.ny;l++) {
    for(k=0;k<zgrid.nx;k++) {
      n=zgrid.nx*l+k;
      z_field[n]=mask;
      vgrid.xy(k,l,x,y);
      status=map_interpolation(vgrid, v_field, mask, x, y, &(z_field[n]));
      }
    }
    
  count=occurence(mask, z_field, zgrid.Hsize());
  
  return(z_field);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *map_extrapolate_z2v(grid_t zgrid, grid_t vgrid, float *z_field, float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes

  Purpose : create vorticity field from z field

  Check :

  Note:

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int    i,j,k,l,m,n;
  float *count;
  float *v_field;
  pocgrd_t ncgrid;

  cdfvar_t variable;
  int shift=0;

  v_field=new float[vgrid.Hsize()];
  count  =new float[vgrid.Hsize()];
  
  if(vgrid.nx==zgrid.nx) shift=1;

  for(l=0;l<vgrid.ny;l++) {
    for(k=0;k<vgrid.nx;k++) {
      n=vgrid.nx*l+k;
      v_field[n]=0.;
      count[n]=0.;
      }
    }
    
  for(l=shift;l<zgrid.ny;l++) {
    for(k=shift;k<zgrid.nx;k++) {
      n=zgrid.nx*l+k;
      if(z_field[n]>10000.) {
//        printf("map_extrapolate_z2v: %d %f\n",n,z_field[n]);
        continue;
        }
      i=k-shift;
      j=l-shift;
      m=vgrid.nx*j+i;
      v_field[m]+=z_field[n];
      count[m]++;
      i=k-shift;
      j=l+1-shift;
      m=vgrid.nx*j+i;
      v_field[m]+=z_field[n];
      count[m]++;
      i=k+1-shift;
      j=l+1-shift;
      m=vgrid.nx*j+i;
      v_field[m]+=z_field[n];
      count[m]++;
      i=k+1-shift;
      j=l-shift;
      m=vgrid.nx*j+i;
      v_field[m]+=z_field[n];
      count[m]++;
      }
    }

  for(l=0;l<vgrid.ny;l++) {
    for(k=0;k<vgrid.nx;k++) {
      n=vgrid.nx*l+k;
      v_field[n]/=count[n];
      if(count[n]==0) {
//        printf("map_extrapolate_z2v failure: %d %d\n",k,l);
        continue;
        }
      }
    }
    
  delete [] count;
  return(v_field);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int parse_GridOptions(const string  & options,  metagrid_t & meta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<string> tokens,keys,values;
  int k;
  string delimiter=" ";
    
  delimiter=" ";
  tokens=string_split(options, delimiter);
  
  delimiter="=";
  for(k=0;k<tokens.size();k++) {
    vector<string> tmp=string_split(tokens[k], delimiter);
    keys.push_back(tmp[0]);
    values.push_back(tmp[1]);
    }
    
  for(k=0;k<keys.size();k++) {
    const string *keyk=&keys[k];
    char *valuek=poc_strdup(values[k].c_str());
    
    if(*keyk=="gridfile") {
      meta.gridfile=valuek;
      continue;
      }
      
    if(*keyk=="maskfile") {
      meta.maskfile=valuek;
      continue;
      }
    
    if(*keyk=="tag") {
      meta.tag=(float) atof(valuek);
      continue;
      }
    
    if(*keyk=="target") {
      meta.target=valuek;
      continue;
      }
    
    const size_t keyl=keyk->length();
    if(keyl<2) continue;
    
    const string
      keystart=keyk->substr(0,keyl-1),
      keyend=keyk->substr(keyl-1);
    gnames_t *gname=0;
    
    /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    Keys can be any combination of :
    - vlon_, vlat_, vmask_ or vtopo_
    - followed by t, u, v or f
    Examples : vlon_t=... vlat_t=... vlon_f=... vlat_f=...
    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
    
    if(keyend=="t")
      gname=&meta.z_gnames;
    if(keyend=="u")
      gname=&meta.u_gnames;
    if(keyend=="v")
      gname=&meta.v_gnames;
    if(keyend=="f")
      gname=&meta.f_gnames;
    if(gname==0) continue;
    
    if(keystart=="vlon_") {
      gname->vlon=valuek;
      continue;
      }
    if(keystart=="vlat_") {
      gname->vlat=valuek;  /* set but never used */
      continue;
      }
    if(keystart=="vmask_") {
      gname->vmask=valuek;
      continue;
      }
    if(keystart=="vtopo_") {
      gname->vtopo=valuek;
      continue;
      }
    
//     if(keystart=="ha_") {
//       gname->ha=valuek;
//       continue;
//       }
//     if(keystart=="hg_") {
//       gname->hg=valuek;
//       continue;
//       }
//     if(keystart=="ua_") {
//       gname->ua=valuek;
//       continue;
//       }
//     if(keystart=="ug_") {
//       gname->ug=valuek;
//       continue;
//       }
//     if(keystart=="va_") {
//       gname->va=valuek;
//       continue;
//       }
//     if(keystart=="vg_") {
//       gname->vg=valuek;
//       continue;
//       }
//     
//     if(keystart=="LSAa_") {
//       gname->LSAa=valuek;
//       continue;
//       }
//     if(keystart=="LSAg_") {
//       gname->LSAg=valuek;
//       continue;
//       }
    
    STDERR_BASE_LINE_FUNC("key \""+*keyk+"\" not recognised\n");
    delete[]valuek;
    }
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t map_duplicategrid(grid_t zgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes

  Purpose : extend Symphonie elevation grid (MECO+1,NECO+1) to cover
            the (u,v) grid

  Check :

  Note:

  zgrid should be (MECO,NECO), but due to Symphonie netcdf outputs,
  it has an additional (last) row and column
  extension adds an additional (first) row and column

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int    i,j,k,l,m,n;
  int    status;
  double x,y;
  double *dx,*dy/*,*dz*/;
  pocgrd_t ncgrid;
  grid_t grid;
  double mask;

  cdfvar_t variable;

  grid.nx=zgrid.nx+1;
  grid.ny=zgrid.ny+1;
  grid.nz=zgrid.nz+1;

  grid.x=new double[grid.nx*grid.ny];
  grid.y=new double[grid.nx*grid.ny];
  grid.z=new double[grid.nx*grid.ny*grid.nz];

  grid.zmask=zgrid.zmask;
  mask=grid.zmask;

  status=map_cartesian_resolution(zgrid, &dx, &dy/*, &dz*/);

  for(l=0;l<grid.ny;l++) {
    for(k=0;k<grid.nx;k++) {
      i=max(k-1,0);
      j=max(l-1,0);
      m=zgrid.nx*j+i;
      x= zgrid.x[m];
      y= zgrid.y[m];
      n=grid.nx*l+k;
      grid.x[n]=x-(i-k+1)*dx[m];
      grid.y[n]=y-(j-l+1)*dy[m];
      }
    }

  grid.modeH=2;
  m=0;
  grid.xmin=grid.x[m];
  grid.ymin=grid.y[m];

  m=grid.nx*grid.ny-1;
  grid.xmax=grid.x[m];
  grid.ymax=grid.y[m];

  return(grid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float *map_duplicate3D(grid_t zgrid, float *buffer, grid_t extended)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes

  Purpose : extend Symphonie elevation grid (MECO+1,NECO+1) to cover
            the (u,v) grid

  Check :

  Note:

  zgrid should be (MECO,NECO), but due to Symphonie netcdf outputs,
  it has an additional (last) row and column
  extension adds an additional (first) row and column

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int    i,j,k,l,m,n,r,s;
  float *z;

  z=new float[extended.nx*extended.ny*extended.nz];

  for(l=0;l<extended.ny;l++) {
    for(k=0;k<extended.nx;k++) {
      i=max(k-1,0);
      j=max(l-1,0);
      for(r=0;r<zgrid.nz;r++) {
        s=r;
        m=zgrid.nx*zgrid.ny*r+zgrid.nx*j+i;
        n=extended.nx*extended.ny*s+extended.nx*l+k;
        z[n]=buffer[m];
        }
      }
    }

  return(z);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *map_duplicate2D(grid_t zgrid, float *buffer, grid_t extended)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes

  Purpose : extend Symphonie elevation grid (MECO+1,NECO+1) to cover
            the (u,v) grid

  Check :

  Note:

  zgrid should be (MECO,NECO), but due to Symphonie netcdf outputs,
  it has an additional (last) row and column
  extension adds an additional (first) row and column

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int    i,j,k,l,m,n;
  float *z;

  z=new float[extended.nx*extended.ny];

  for(l=0;l<extended.ny;l++) {
    for(k=0;k<extended.nx;k++) {
      i=max(k-1,0);
      j=max(l-1,0);
      m=zgrid.nx*j+i;
      n=extended.nx*l+k;
      z[n]=buffer[m];
      }
    }

  return(z);
}
