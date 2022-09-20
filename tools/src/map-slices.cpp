
/*******************************************************************************

  T-UGO tools, 2006-2015

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief map slices definitions
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "tools-structures.h"

#include "map.h"
#include "map.def"
#include "geo.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mapc_index(const grid_t & grid, double x, double y, int *ktrue, int *ltrue)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#warning mapc_index is not used but map_index02 is used. Both are badly undocumented.
{
  int k0,l0;
  double ux,uy,un;
  int status,m0;
  int attempt;
  int x_out=0,y_out=0;
  

  if(x < grid.xmin) {
    x_out=1;
    x=grid.xmin;
    }
    
  if (x > grid.xmax) {
    x_out=1;
    x=grid.xmax;
    }
 
  if(y < grid.ymin) {
    y_out=1;
    y=grid.ymin;
    }
 
  if(y > grid.ymax) {
    y_out=1;
    y=grid.ymax;
    }

  k0=grid.nx/2;
  l0=grid.ny/2;

  for(attempt=1;attempt<6;attempt++) {
    m0=grid.nx*l0+k0;
    ux= grid.x[m0+1]-grid.x[m0];
    uy= grid.y[m0+1]-grid.y[m0];
    un=sqrt(ux*ux+uy*uy);
    ux=ux/un;
    uy=uy/un;
    k0= k0+int( ((x-grid.x[m0])*ux+(y-grid.y[m0])*uy)/un );
    if(k0<0){
      k0=0;
      }
    if(k0>grid.nx-1){
      k0=grid.nx-1;
    }

    ux= grid.x[grid.nx+m0]-grid.x[m0];
    uy= grid.y[grid.nx+m0]-grid.y[m0];
    un=sqrt(ux*ux+uy*uy);
    ux=ux/un;
    uy=uy/un;

    l0= l0+int( ((x-grid.x[m0])*ux+(y-grid.y[m0])*uy)/un  );
    if(l0<0){
      l0=0;
      }
    if(l0>grid.ny-1){
      l0=grid.ny-1;
      }
    }
  
  m0=grid.nx*l0+k0;

  ux= grid.x[m0+1]-grid.x[m0];
  uy= grid.y[m0+1]-grid.y[m0];
  un=sqrt(ux*ux+uy*uy);
  ux=ux/un;
  uy=uy/un;

  status=0;

  k0= k0+int( ((x-grid.x[m0])*ux+(y-grid.y[m0])*uy)/un );
  if(k0<0) {
    k0=-1;
    status=-1;
    }
  if(k0>grid.nx-1) {
    k0=grid.nx;
    status=-1;
    }

  ux= grid.x[grid.nx+m0]-grid.x[m0];
  uy= grid.y[grid.nx+m0]-grid.y[m0];
  un=sqrt(ux*ux+uy*uy);
  ux=ux/un;
  uy=uy/un;

  l0= l0+int( ((x-grid.x[m0])*ux+(y-grid.y[m0])*uy)/un  );
  if(l0<0) {
    l0=-1;
    status=-1;
    }
  if(l0>grid.ny-1) {
    l0=grid.ny;
    status=-1;
    }

  *ktrue=k0;
  *ltrue=l0;

 return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int check_vertical_direction(const double *z, double zmask, size_t HSize, int nz, bool *increasing, double *factor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// check vertical direction of depths
/**
\param *z array[HSize*nz] of depths
\param zmask mask value of \c *z
\param HSize
\param nz number of levels
\param *increasing set whether the depths are increasing with the level indexes
\param *factor set to the sign of the depths
\returns 0 on success or -1 if all z[] are masked.
*/
/*----------------------------------------------------------------------------*/
{
  int k,l,m,status=-1;

  double *depths=NULL,h,mean;
  
  exitIfNull(
    depths=new double[nz]
    );

  for(m=0;m<HSize;m++){
    mean=0;
    l=0;
    for(k=0;k<nz;k++){
      h=z[m+k*HSize];
      if(h != zmask){
        depths[l]=h;
        mean+=h;
        l++;
        }
      }
    if(l>1){
//       *factor=copysign(1.,depths[0]);
      *factor=copysign(1.,mean);
      //if(depths[0]==depths[1]) /* FOR SOME REASON THIS FAILS WHEN DEPTH IS STORED AS float IN THE FILE */
      if(fabs(depths[0]-depths[1])<0.001)
        continue;
//         TRAP_ERR_EXIT(-1,"invalid vertical grid: 2 identical consecutive levels at %g\n",depths[0]);
      *increasing= (*factor)*depths[0] < (*factor)*depths[1] ;
      status=0;
      goto cleanUp;
      }
    }
  
cleanUp:
  delete[]depths;
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int check_vertical_direction(const grid_t & grid, bool *increasing, double *factor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  size_t Hsize;
  
  switch(grid.modeV){
  case 1:
    Hsize=1;
    break;
  case 3:
    Hsize=grid.Hsize();
    break;
  default:
    TRAP_ERR_EXIT(ENOEXEC,"%s(grid,) with grid.modeV=%d\n",__func__,grid.modeV);
    }
  
  status=check_vertical_direction(grid.z,grid.zmask,Hsize,grid.nz,increasing,factor);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int map_sliceH_template(const grid_t & grid, double level, T *buffer, T mask, T *slice)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  bool increasing;
  double factor;

  status=check_vertical_direction(grid,&increasing,&factor);
  if(status)TRAP_ERR_RETURN(status,1,"check_vertical_direction() error\n");
  
/*------------------------------------------------------------------------------
  compute the horizontal slice array */

  #pragma omp parallel
  {
  int j,i;
  T *vector;
  double *depths,h;
  int k,l,m;
  
  vector=new T[grid.nz];
  depths=new double[grid.nz];
  
  #pragma omp for
  for(j=0;j<grid.ny;j++){
    
    for(i=0;i<grid.nx;i++){
      m=j*grid.nx+i;
      l=0;
      
      if(increasing){
        for(k=0;k<grid.nz;k++){
          h=grid.GetZ(m,k);
          if(h != grid.zmask){
            depths[l] =factor*h;
            vector[l] =buffer[m+k*grid.nx*grid.ny];
            l++;
            }
          }
        }
      else{
        for(k=grid.nz-1;k>=0;k--){
          h=grid.GetZ(m,k);
          if(h != grid.zmask){
            depths[l] =factor*h;
            vector[l] =buffer[m+k*grid.nx*grid.ny];
            l++;
            }
          }
        }
      
      if(l> 0){
        status=map_interpolate1D(vector, depths, mask,l, level, &slice[m]);
        }
      else{
        slice[m]=mask;
        }
      
      }
    }
  
  delete[]vector;
  delete[]depths;
  }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_sliceH(const grid_t & grid, double level, float *buffer, float mask, float *slice)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=map_sliceH_template(grid,level,buffer,mask,slice);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_sliceH(const grid_t & grid, double level, double *buffer, double mask, double *slice)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=map_sliceH_template(grid,level,buffer,mask,slice);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_sliceH(const grid_t & grid, double level, double *buffer, double mask, complex<double> *slice)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  int n=grid.Hsize();
  double *d=new double[n];
  
  result=map_sliceH_template(grid,level,buffer,mask,d);
  
  valcpy(slice,d,n);
  
  delete[]d;
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_sliceH(const grid_t & grid, double level, complex<double> *buffer, complex<double> mask, complex<double> *slice)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=map_sliceH_template(grid,level,buffer,mask,slice);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int map_sliceV_template(const grid_t & grid, double *x, double *y, int nloc, T *buffer, T mask,grid_t *slice_grid, T **slice, int verbose=1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int i,k,m,status;
  bool increasing=false;

  T *bufferm=NULL,value;
  T *local=NULL;
  double *depths=NULL,h,factor;
  
  range_t<double> xR,yR,zR;
  xR.init(x,nloc);
  zR.init(y,nloc);

  status=check_vertical_direction(grid,&increasing,&factor);
  if(status)TRAP_ERR_RETURN(status,verbose,"check_vertical_direction() error\n");
  
  if(verbose){
    printf("increasing=%d,factor=%g\n",(int)increasing,factor);
    }

/*---------------------------------------------------------------------*//**<h1>
  define the section grid </h1>*/

  slice_grid->nx=nloc;
  slice_grid->ny=grid.nz;
  /** \note For easier interpolation later :
  - the grid_t::x, grid_t::x and grid_t::z parameters of \c slice_grid
  - and **slice
  are transposed, therefore the grid_t::modeH is set to -2. */
  slice_grid->modeH=-2;
  slice_grid->nz=1;
  slice_grid->modeV=-1;
  
  const int nxy=slice_grid->Hsize();

  exitIfNull(
    slice_grid->x=new double [nxy]
    );
  exitIfNull(
    slice_grid->y=new double [nxy]
    );
  exitIfNull(
    slice_grid->z=new double [nxy]
    );

  exitIfNull(
    local=new T[nxy]
    );

  int64_t accel;
  double lon,lat;
  
/*---------------------------------------------------------------------*//**<h1>
  get depths and values along the section </h1>*/
  for(k=0;k<grid.nz;k++){
    m=grid.nx*grid.ny*k;
    switch(grid.modeV){
    case 1:
      depths=&grid.z[k];
      break;
    case 3:
      depths=&grid.z[m];
      break;
    default:
      TRAP_ERR_EXIT(ENOEXEC,"%s(grid,) with grid.modeV=%d\n",__func__,grid.modeV);
      }
    bufferm=&buffer[m];
    
    if(grid.list){
      accel=-1;
      }
    
    for(i=0;i<slice_grid->nx;i++){
      if(increasing)
        m=i*slice_grid->ny+k;
      else
        m=i*slice_grid->ny+grid.nz-k-1;
      
      lon=x[i];
      lat=y[i];
      
      slice_grid->x[m]=lon;
      slice_grid->z[m]=lat;
      
      switch(grid.modeV){
      case 1:
        h=*depths;
        break;
      case 3:
        if(grid.list)
          index_interpolation(grid,lon,lat,&accel,depths,grid.zmask,&h,verbose);
        else
          status=map_interpolation(grid,depths,grid.zmask,lon,lat,&h);
        break;
        }
      if(grid.modeH!=2 or grid.list!=0)
        index_interpolation(grid,lon,lat,&accel,bufferm,mask,&value,verbose);
      else
        status=map_interpolation(grid,bufferm,mask,lon,lat,&value);
      if(verbose>=2)STDERR_BASE_LINE("[%d]: %g %g %g %g\n",i,lon,lat,h,value);
      if(h != grid.zmask){
        h*=factor;
        yR<<h;
        }
      else{
        h=NAN;
        value=mask;
        }
      slice_grid->y[m]=h;
      local[m]=value;
      }
    }
  
/*---------------------------------------------------------------------*//**<h1>
  fill holes in the depths </h1>*/
  for(i=0;i<slice_grid->nx;i++){
    m=i*slice_grid->ny;
    
    for(k=0;k<grid.nz && isnan(slice_grid->y[m+k]);k++);
    if(k==grid.nz){
      //STDERR_BASE_LINE("fool-proofing for %d\n",i);
      slice_grid->y[m]=INFINITY;
      }
    h=slice_grid->y[m];
    
    for(k=0;k<grid.nz;k++,m++){
      if(isnan(slice_grid->y[m])){
        slice_grid->y[m]=h;
        }
      else{
        h=slice_grid->y[m];
        }
      }
    
    }
  
  slice_grid->xmax=xR.max;
  slice_grid->xmin=xR.min;
  slice_grid->ymax=yR.max;
  slice_grid->ymin=yR.min;
  slice_grid->zmax=zR.max;
  slice_grid->zmin=zR.min;
  
  slice_grid->dx=(slice_grid->xmax-slice_grid->xmin)/(slice_grid->nx-1);
  slice_grid->dy=(slice_grid->ymax-slice_grid->ymin)/(slice_grid->ny-1);
  slice_grid->dz=(slice_grid->zmax-slice_grid->zmin)/(slice_grid->nz-1);
  
  if(verbose){
    slice_grid->print();
    }

  status=0;

  *slice=local;

  return(status);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_sliceV(const grid_t & grid, double *x, double *y, int nloc, float *buffer, float mask,grid_t *slice_grid, float **slice, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=map_sliceV_template(grid,x,y,nloc,buffer,mask,slice_grid,slice,verbose);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_sliceV(const grid_t & grid, double *x, double *y, int nloc, double *buffer, double mask,grid_t *slice_grid, double **slice, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=map_sliceV_template(grid,x,y,nloc,buffer,mask,slice_grid,slice,verbose);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_sliceV(const grid_t & grid, double *x, double *y, int nloc, double *buffer, double mask,grid_t *slice_grid, complex<double> **slice, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  double *d;
  
  status=map_sliceV_template(grid,x,y,nloc,buffer,mask,slice_grid,&d,verbose);
  
  int n=slice_grid->Hsize();
  
  *slice=new complex<double>[n];
  
  valcpy(*slice,d,n);
  
  delete[]d;
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_sliceV(const grid_t & grid, double *x, double *y, int nloc, complex<double> *buffer, complex<double> mask,grid_t *slice_grid, complex<double> **slice, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=map_sliceV_template(grid,x,y,nloc,buffer,mask,slice_grid,slice,verbose);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_profile01(const grid_t & grid, double x, double y, float *buffer, float mask,
                 grid_t *slice_grid, float **profile_x, float **profile_y, int *nloc)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int i,j,k,l,m,status;
  int negative=0,increasing=0,z_levels;

  float *vector2=NULL,factor;
  float *physic=NULL,value;
  double *depths=NULL,h,href,*vector1=NULL,*levels=NULL;

  exitIfNull(
    vector1=(double *) malloc(grid.nz*sizeof(double))
    );
  exitIfNull(
    vector2=(float *)  malloc(grid.nz*sizeof(float))
    );

/*------------------------------------------------------------------------------
  check the vertical arrangement */

  exitIfNull(
    depths=(double *) malloc(grid.nz*sizeof(double))
    );
  exitIfNull(
    levels=(double *) malloc(grid.nz*sizeof(double))
    );

/*------------------------------------------------------------------------------
  horizontal levels? */
  z_levels=1;
  for(k=0;k<grid.nz;k++) {
    for(m=0;m<grid.ny*grid.nx;m++) {
      h=grid.z[m+k*grid.nx*grid.ny];
      if (h!=grid.zmask) {
        href=h;
        levels[k]=h;
        break;
        }
      }
    for(m=m;m<grid.ny*grid.nx;m++) {
      h=grid.z[m+k*grid.nx*grid.ny];
      if ((h!=grid.zmask) && (h!=href)) {
        z_levels=0;
        k=grid.nz-1;
        break;
        }
      }
    }

/*------------------------------------------------------------------------------
  increasing/descending positive/negative levels? */
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      l=0;
      for(k=0;k<grid.nz;k++) {
        h=grid.z[m+k*grid.nx*grid.ny];
        if(h != grid.zmask) {
          depths[l]=grid.z[m+k*grid.nx*grid.ny];
          l++;
          }
        }
      if(l>1) {
        if (depths[0] < depths[1]) increasing=1;
        if(increasing==1) {
          if (depths[0] < 0.)      negative=1;
          }
        else
          {
          if (depths[1] < 0.)      negative=1;
          }
        goto next;
        }
      }
    }

/*------------------------------------------------------------------------------
  no valid depth available */

  status=-1;
  goto finished;

 next:

  free(depths);
/*------------------------------------------------------------------------------
  depth taken as negative for interpolation*/

  factor=1.;

  if(negative==0) {
    factor=-1.;
    increasing=1-increasing;
    }

  printf("negative= %d,increasing= %d,z_levels= %d\n",negative,increasing,z_levels);

/*------------------------------------------------------------------------------
  define the section grid*/


  slice_grid->nx=400;
  slice_grid->ny=400;

  slice_grid->xmin =1.e+10;
  slice_grid->xmax=-1.e+10;
  slice_grid->ymin =1.e+10;
  slice_grid->ymax=-1.e+10;

  slice_grid->x= NULL;
  slice_grid->y= NULL;
  slice_grid->z= NULL;

/*------------------------------------------------------------------------------
  get depth min/max along the section */
  for(k=0;k<grid.nz;k++) {
    m=grid.nx*grid.ny*k;
    depths=&(grid.z[m]);
    physic=&(buffer[m]);
    status= map_interpolation(grid,depths,grid.zmask,x,y, &h);
    status= map_interpolation(grid,physic,mask,x,y, &value);
    if((h != grid.zmask) && (value != mask)) {
      updatemin(&slice_grid->ymin,factor*h);
      updatemax(&slice_grid->ymax,factor*h);
      updatemin(&slice_grid->xmin,value);
      updatemax(&slice_grid->xmax,value);
      }
    }

  slice_grid->dx=(slice_grid->xmax-slice_grid->xmin)/(slice_grid->nx-1);
  slice_grid->dy=(slice_grid->ymax-slice_grid->ymin)/(slice_grid->ny-1);

/*------------------------------------------------------------------------------
  first create a vertical profile of buffer at subx[i],suby[i]*/
  l=0;
  if(increasing==1) {
    for(k=0;k<grid.nz;k++) {
      m=grid.nx*grid.ny*k;
      depths=&(grid.z[m]);
      physic=&(buffer[m]);
      status= map_interpolation(grid,depths,grid.zmask,x,y, &h);
      if(h != grid.zmask) {
        vector1[l]=factor*h;
        status= map_interpolation(grid,physic,mask,x,y, &vector2[l]);
        l++;
        }
      }
    }
  else
    {
    for(k=grid.nz-1;k>=0;k--) {
      m=grid.nx*grid.ny*k;
      depths=&(grid.z[m]);
      physic=&(buffer[m]);
      status= map_interpolation(grid,depths,grid.zmask,x,y, &h);
      if(h != grid.zmask) {
        vector1[l]=factor*h;
        status= map_interpolation(grid,physic,mask,x,y, &vector2[l]);
        l++;
        }	
      }	
    }

  *nloc=l;

  status=-1;
  if(*nloc==0) goto finished;

  exitIfNull(
    *profile_x=(float *) malloc(l*sizeof(float))
    );
  exitIfNull(
    *profile_y=(float *) malloc(l*sizeof(float))
    );
  for(l=0;l<*nloc;l++) {
    (*profile_x)[l]=vector2[l];
    (*profile_y)[l]=vector1[l];
    }

  slice_grid->zmin=0;
  slice_grid->zmax=0;
  slice_grid->dz=0;
  slice_grid->nz=1;

  slice_grid->modeH= 0;
  slice_grid->modeV=-1;

  slice_grid->circular=0;
  slice_grid->overlapped=0;
  slice_grid->connex=1;

  mapc_printgrid3d(*slice_grid);

  status=0;

 finished:

  if (levels   != NULL) free(levels);
  if (vector1  != NULL) free(vector1);
  if (vector2  != NULL) free(vector2);

  return(status);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_chklevels(const grid_t & grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m;
  int z_levels;

  double h,href,*levels=NULL;

  exitIfNull(
    levels=(double *) malloc(grid.nz*sizeof(double))
    );

/*------------------------------------------------------------------------------
  horizontal levels? */
  z_levels=1;
  for(k=0;k<grid.nz;k++) {
    for(m=0;m<grid.ny*grid.nx;m++) {
      h=grid.z[m+k*grid.nx*grid.ny];
      if (h!=grid.zmask) {
        href=h;
        break;
        }
      }
    for(m=m;m<grid.ny*grid.nx;m++) {
      h=grid.z[m+k*grid.nx*grid.ny];
      if ((h!=grid.zmask) && (h!=href)) {
        z_levels=0;
        k=grid.nz-1;
        break;
        }
      }
    levels[k]=h;
    }

  return(z_levels);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_diffusion01(const grid_t & grid, float *buffer, float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int i,j,k,m,n,status;
  float z;
  double x,y;
  int *count[2];
  float *tmp=NULL;

  exitIfNull(
    count[0]=(int *) malloc(grid.nz*sizeof(int))
    );
  exitIfNull(
    count[1]=(int *) malloc(grid.nz*sizeof(int))
    );

  exitIfNull(
    tmp=(float *) malloc(grid.nx*grid.ny*grid.nz*sizeof(float))
    );

  for (k=0;k<grid.nz;k++) {
    count[0][k]=0;
    for(j=0;j<grid.ny;j++) {
      for(i=0;i<grid.nx;i++) {
        m=k*grid.nx*grid.ny+j*grid.nx+i;
        tmp[m]=buffer[m];
        if(buffer[m]==mask) count[0][k]++;
        }
      }
    }

  for (k=0;k<grid.nz;k++) {
    for(j=0;j<grid.ny;j++) {
      for(i=0;i<grid.nx;i++) {
        n=j*grid.nx+i;
        m=k*grid.nx*grid.ny+j*grid.nx+i;
        if(buffer[m]==mask) {
          x=grid.x[n]+0.05;
          y=grid.y[n];
          status= map_interpolation(grid,&(buffer[k*grid.nx*grid.ny]),mask,x,y,&z);
          if(z!=mask) {
            tmp[m]=z;
            }
          }
        else continue;
        if(buffer[m]==mask) {
          x=grid.x[n]-0.05;
          y=grid.y[n];
          status= map_interpolation(grid,&(buffer[k*grid.nx*grid.ny]),mask,x,y,&z);
          if(z!=mask) {
            tmp[m]=z;
            }
          }
        if(buffer[m]==mask) {
          x=grid.x[n];
          y=grid.y[n]+0.05;
          status= map_interpolation(grid,&(buffer[k*grid.nx*grid.ny]),mask,x,y,&z);
          if(z!=mask) {
            tmp[m]=z;
            }
          }
        if(buffer[m]==mask) {
          x=grid.x[n];
          y=grid.y[n]-0.05;
          status= map_interpolation(grid,&(buffer[k*grid.nx*grid.ny]),mask,x,y,&z);
          if(z!=mask) {
            tmp[m]=z;
            }
          }
        }
      }
    }

  for (k=0;k<grid.nz;k++) {
    count[1][k]=0;
    for(j=0;j<grid.ny;j++) {
      for(i=0;i<grid.nx;i++) {
        m=k*grid.nx*grid.ny+j*grid.nx+i;
        if(buffer[m]==mask) {
          buffer[m]=tmp[m];
          }
        if(buffer[m]==mask) count[1][k]++;
        }
      }
//    printf("layer %d: %d %d mask values\n",k,count[0][k],count[1][k]);
    }

  free(count[0]);
  free(count[1]);
  free(tmp);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_profile02(const grid_t & grid, double x, double y, float *buffer, float mask,
                float **profile_x, float **profile_y, int *nloc)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,status;
  int negative=0,increasing=0/*,z_levels*/;

  float *vector2=NULL,factor;
  float *physic=NULL;
  double *depths=NULL,h,*vector1=NULL,*levels=NULL;

  exitIfNull(
    vector1=new double[grid.nz]
    );
  exitIfNull(
    vector2=new float[grid.nz]
    );

/*------------------------------------------------------------------------------
  check the vertical arrangement */

  exitIfNull(
    depths=new double[grid.nz]
    );
  exitIfNull(
    levels=new double[grid.nz]
    );

/*------------------------------------------------------------------------------
  increasing/descending positive/negative levels? */
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
/*------------------------------------------------------------------------------
      focus check on valid cells */
      if(buffer[m]==mask) continue;
      l=0;
      for(k=0;k<grid.nz;k++) {
        h=grid.z[m+k*grid.nx*grid.ny];
        if(h != grid.zmask) {
          depths[l]=grid.z[m+k*grid.nx*grid.ny];
          l++;
          }
        }
      if(l>1) {
        if (depths[0] < depths[1]) increasing=1;
        if(increasing==1) {
          if (depths[0] < 0.)      negative=1;
          }
        else {
          if (depths[1] < 0.)      negative=1;
          }
        goto next;
        }
      }
    }

/*------------------------------------------------------------------------------
  no valid depth available */
  status=-1;
  goto error;

 next:

  free(depths);
/*------------------------------------------------------------------------------
  depth taken as negative for interpolation*/

  factor=1.;

  if(negative==0) {
    factor=-1.;
    increasing=1-increasing;
    }

  /* printf("negative= %d,increasing= %d,z_levels= %d\n",negative,increasing,z_levels); */

/*------------------------------------------------------------------------------
  first create a vertical profile of buffer at subx[i],suby[i]*/
  l=0;
  if(increasing==1) {
    for(k=0;k<grid.nz;k++) {
      m=grid.nx*grid.ny*k;
      depths=&(grid.z[m]);
      physic=&(buffer[m]);
      status= map_interpolation(grid,depths,grid.zmask,x,y, &h);
      if(status!=0) goto error;
      if(h != grid.zmask) {
        vector1[l]=factor*h;
        status= map_interpolation(grid,physic,mask,x,y, &vector2[l]);
/*         if(status!=0) goto error; */
        l++;
        }
      }
    }
  else {
    for(k=grid.nz-1;k>=0;k--) {
      m=grid.nx*grid.ny*k;
      depths=&(grid.z[m]);
      physic=&(buffer[m]);
      status= map_interpolation(grid,depths,grid.zmask,x,y, &h);
      if(status!=0) goto error;
      if(h != grid.zmask) {
        vector1[l]=factor*h;
        status= map_interpolation(grid,physic,mask,x,y, &vector2[l]);
/*         if(status!=0) goto error; */
        l++;
        }	
      }	
    }

  *nloc=l;

  status=-1;
  if(*nloc==0) goto finished;


  if(*profile_x!=NULL) free(*profile_x);
  if(*profile_y!=NULL) free(*profile_y);

  exitIfNull(
    *profile_x=(float *) malloc(l*sizeof(float))
    );
  exitIfNull(
    *profile_y=(float *) malloc(l*sizeof(float))
    );

  for(l=0;l<*nloc;l++) {
    (*profile_x)[l]=vector2[l];
    (*profile_y)[l]=vector1[l];
    if(vector2[l]!=mask) status=0;
    }

 finished:

  if (levels   != NULL) free(levels);
  if (vector1  != NULL) free(vector1);
  if (vector2  != NULL) free(vector2);
  return(status);

 error:

  status=-1;
  if (levels   != NULL) free(levels);
  if (vector1  != NULL) free(vector1);
  if (vector2  != NULL) free(vector2);
  return(status);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_sliceV3D(const grid_t & grid, double *x, double *y, int nloc, float *buffer, float mask, grid_t *slice_grid, float **slice, float **bottom,int xoption,int yoption)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int i,j,k,l,m,status;
  int negative=0,increasing=0,z_levels;
  int n0,n1;
  int check=0;
  int compute_hmin,compute_hmax;

  double t,p;
  double tx,ty,tn;
  double *subx=NULL,*suby=NULL;

  float z,*vector2=NULL,*local=NULL,factor,level;
  float *physic=NULL,value;
  double *depths=NULL,h,href,*vector1=NULL,*levels=NULL;

  exitIfNull(
    vector1=(double *) malloc(grid.nz*sizeof(double))
    );
  exitIfNull(
    vector2=(float *)  malloc(grid.nz*sizeof(float))
    );

/*------------------------------------------------------------------------------
  initialize the section grid to empty grid*/
  slice_grid->nx=0;
  slice_grid->ny=0;

  slice_grid->xmin =1.e+10;
  slice_grid->xmax=-1.e+10;

/*------------------------------------------------------------------------------
  removed */

/*
  slice_grid->ymin =1.e+10;
  slice_grid->ymax=-1.e+10;
*/

  slice_grid->x= NULL;
  slice_grid->y= NULL;
  slice_grid->z= NULL;

/*------------------------------------------------------------------------------
  check the vertical arrangement */

  exitIfNull(
    depths=(double *) malloc(grid.nz*sizeof(double))
    );
  exitIfNull(
    levels=(double *) malloc(grid.nz*sizeof(double))
    );

/*------------------------------------------------------------------------------
  horizontal levels? */
  z_levels=1;
  for(k=0;k<grid.nz;k++) {
    for(m=0;m<grid.ny*grid.nx;m++) {
      h=grid.z[m+k*grid.nx*grid.ny];
      if (h!=grid.zmask) {
        href=h;
        levels[k]=h;
        break;
        }
      }
    for(m=m;m<grid.ny*grid.nx;m++) {
      h=grid.z[m+k*grid.nx*grid.ny];
      if ((h!=grid.zmask) && (h!=href)) {
        z_levels=0;
        k=grid.nz-1;
        break;
        }
      }
    }

/*------------------------------------------------------------------------------
  increasing/descending positive/negative levels? */
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      t=grid.x[m];
      p=grid.y[m];
      l=0;
      for(k=0;k<grid.nz;k++) {
        h=grid.z[m+k*grid.nx*grid.ny];
        if(h != grid.zmask) {
          depths[l]=grid.z[m+k*grid.nx*grid.ny];
          l++;
          }
        }
      if(l>1) {
        if (depths[0] < depths[1]) increasing=1;
        if(increasing==1) {
          if (depths[0] < 0.)      negative=1;
          }
        else {
          if (depths[1] < 0.)      negative=1;
          }
        goto next;
        }
      }
    }

/*------------------------------------------------------------------------------
  no valid depth available */

  status=-1;
  goto finished;

 next:

  free(depths);
/*------------------------------------------------------------------------------
  depth taken as negative for interpolation*/

  factor=1.;

  if(negative==0) {
    factor=-1.;
    increasing=1-increasing;
    }

  if(check==1) printf("negative= %d,increasing= %d,z_levels= %d\n",negative,increasing,z_levels);

/*------------------------------------------------------------------------------
  define the section grid*/
  slice_grid->ny=400;

  if(nloc==2) slice_grid->nx=400;
  else slice_grid->nx=nloc;

  exitIfNull(
    subx=(double *) malloc(slice_grid->nx*sizeof(double))
    );
  exitIfNull(
    suby=(double *) malloc(slice_grid->nx*sizeof(double))
    );

  tx=(x[nloc-1]-x[0])/(slice_grid->nx-1);
  ty=(y[nloc-1]-y[0])/(slice_grid->nx-1);
  tn=sqrt(tx*tx+ty*ty);

/*------------------------------------------------------------------------------
  select first and last following abscisse option*/
  switch (xoption) {
    case 0:
    if(x[nloc-1]<x[0]) {
      n0=nloc-1;
      n1=0;
      }
    else {
      n0=0;
      n1=nloc-1;
      }
    break;

    case 1:
    if(y[nloc-1]<y[0]) {
      n0=nloc-1;
      n1=0;
      }
    else {
      n0=0;
      n1=nloc-1;
      }
    break;

    default:
    n0=0;
    n1=nloc-1;
    break;
    }

  tx=(x[n1]-x[n0])/(slice_grid->nx-1);
  ty=(y[n1]-y[n0])/(slice_grid->nx-1);
  tn=sqrt(tx*tx+ty*ty);

  if(nloc==2) {
/*------------------------------------------------------------------------------
    sample a 2-position entry*/
    for(i=0;i<slice_grid->nx;i++) {
      t=x[n0]+i*tx;
      t=map_recale(grid,t);
      p=y[n0]+i*ty;
      subx[i]=t;
      suby[i]=p;
      }
    }
  else {
/*------------------------------------------------------------------------------
    use a multiple-position entry*/
    for(i=0;i<slice_grid->nx;i++) {
      subx[i]=x[i];
      suby[i]=y[i];
      }
    }

  exitIfNull(
    slice_grid->x=(double *) malloc(slice_grid->nx*slice_grid->ny*sizeof(double))
    );
  exitIfNull(
    slice_grid->y=(double *) malloc(slice_grid->nx*slice_grid->ny*sizeof(double))
    );
  slice_grid->z= NULL;

  exitIfNull(
    local=(float *) malloc(slice_grid->nx*slice_grid->ny*sizeof(float))
    );
  exitIfNull(
    *bottom=(float *) malloc(slice_grid->nx*sizeof(float))
    );

/*------------------------------------------------------------------------------
  commented, must be done by calling routine
  slice_grid->ymin =1.e+10;
  slice_grid->ymax=-1.e+10;
*/

  if (slice_grid->ymin== 1.e+10)
    compute_hmin=1;
  else
    compute_hmin=0;

  if (slice_grid->ymax==-1.e+10)
    compute_hmax=1;
  else
    compute_hmax=0;

/*------------------------------------------------------------------------------
  get depth overall min/max and local  minimum along the section */
  //for all positions
  for(i=0;i<slice_grid->nx;i++) {
    (*bottom)[i]=1.e+5;
    //for all levels
    for(k=0;k<grid.nz;k++) {
      m=grid.nx*grid.ny*k;
      depths=&(grid.z[m]);
      physic=&(buffer[m]);
      status= map_interpolation(grid,depths,grid.zmask,subx[i],suby[i], &h);
      status= map_interpolation(grid,physic,mask,subx[i],suby[i], &value);
      if((h != grid.zmask) && (value != mask)) {
        if(compute_hmin==1) updatemin(&slice_grid->ymin,factor*h);
        if(compute_hmax==1) updatemax(&slice_grid->ymax,factor*h);
        updatemin(&(*bottom)[i],factor*h);
        }
      }
    }

/*------------------------------------------------------------------------------
  regular sampling for y*/
  slice_grid->dy=(slice_grid->ymax-slice_grid->ymin)/(slice_grid->ny-1);

/*------------------------------------------------------------------------------
  Build section horizontal and vertical axis */
  for(j=0;j<slice_grid->ny;j++) {
    tn=0;
    for(i=0;i<slice_grid->nx;i++) {
      m=i+slice_grid->nx*j;
/*------------------------------------------------------------------------------
      horizontal axis*/
      switch (xoption) {
/*------------------------------------------------------------------------------
        plot along longitudes*/
        case 0:
        slice_grid->x[m]=subx[i];
        break;

/*------------------------------------------------------------------------------
        plot along latitudes*/
        case 1:
        slice_grid->x[m]=suby[i];
        break;

/*------------------------------------------------------------------------------
        plot along distance in degrees*/
        case 2:
          if(i==0)
            tn=0;
          else {
            tx=subx[i]-subx[0];
            ty=suby[i]-suby[0];
            tn+=sqrt(tx*tx+ty*ty);
            }
          slice_grid->x[m]=tn;
/*         slice_grid->x[m]=i*tn;  */
        break;

/*------------------------------------------------------------------------------
        plot along distance in km*/
        case 3:
          if(i==0)
            tn=0;
          else
            tn+=geo_km(subx[i],suby[i],subx[i-1],suby[i-1]);
          slice_grid->x[m]=tn;
        break;
        }
/*------------------------------------------------------------------------------
      vertical axis*/
      if(yoption==0) slice_grid->y[m]=slice_grid->ymin+j*slice_grid->dy;
      else           slice_grid->y[m]=j*slice_grid->dy;
/*
      if(z_levels==0) slice_grid->y[m]=slice_grid->ymin+j*slice_grid->dy;
      else            slice_grid->y[m]=levels[i];
*/
      }
    }

/*------------------------------------------------------------------------------
  Interpolation*/
  for(i=0;i<slice_grid->nx;i++) {
/*------------------------------------------------------------------------------
    first create a vertical profile of buffer at subx[i],suby[i]*/
    l=0;
    if(increasing==1) {
      for(k=0;k<grid.nz;k++) {
        m=grid.nx*grid.ny*k;
        depths=&(grid.z[m]);
        physic=&(buffer[m]);
        status= map_interpolation(grid,depths,grid.zmask,subx[i],suby[i], &h);
        if(h != grid.zmask) {
          vector1[l]=factor*h;
          status= map_interpolation(grid,physic,mask,subx[i],suby[i], &vector2[l]);
          l++;
          }
/*
        else
          {
          if(k==0) printf("#map_sliceV3D: %f %f\n",subx[i],suby[i]);
          }
*/
        }
      }
    else {
      for(k=grid.nz-1;k>=0;k--) {
        m=grid.nx*grid.ny*k;
        depths=&(grid.z[m]);
        physic=&(buffer[m]);
        status= map_interpolation(grid,depths,grid.zmask,subx[i],suby[i], &h);
        if(h != grid.zmask) {
          vector1[l]=factor*h;
          status= map_interpolation(grid,physic,mask,subx[i],suby[i], &vector2[l]);
          l++;
          }	
        }	
      }
/*------------------------------------------------------------------------------
    then interpolate on the slice profile*/
    if(l !=0) {
      for(j=0;j<slice_grid->ny;j++) {
        m=slice_grid->nx*j+i;
        if(yoption==0) level=slice_grid->y[m];
        else           level=slice_grid->y[m]+(*bottom)[i];
        status=map_interpolate1D(vector2, vector1, mask,l, level, &z);
        local[m]=z;
        }
      }
    else
      for(j=0;j<slice_grid->ny;j++) {
        m=slice_grid->nx*j+i;
        local[m]=mask;
        }
    }

  if(yoption!=0) {
    for(i=0;i<slice_grid->nx;i++) {
      (*bottom)[i]*=-1.;
      }
    }

  slice_grid->xmin=slice_grid->x[0];
  slice_grid->xmax=slice_grid->x[slice_grid->nx-1];
  slice_grid->dx=(slice_grid->xmax-slice_grid->xmin)/(slice_grid->nx-1);

  slice_grid->ymin=slice_grid->y[0];
  slice_grid->ymax=slice_grid->y[slice_grid->nx*slice_grid->ny-1];
  slice_grid->dy=(slice_grid->ymax-slice_grid->ymin)/(slice_grid->ny-1);

  slice_grid->zmin=0;
  slice_grid->zmax=0;
  slice_grid->dz=0;
  slice_grid->nz=1;

  slice_grid->modeH= 2;
  slice_grid->modeV=-1;

  slice_grid->circular=0;
  slice_grid->overlapped=0;
  slice_grid->connex=1;

  if(check==1) mapc_printgrid3d(*slice_grid);

  status=0;

 finished:

  if (levels   != NULL) free(levels);
  if (vector1  != NULL) free(vector1);
  if (vector2  != NULL) free(vector2);
  if (subx     != NULL) free(subx);
  if (suby     != NULL) free(suby);

  *slice=local;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_sliceV3Dvector(const grid_t & grid, double *x, double *y, int nloc, float *bufferx, float *buffery, float mask,grid_t *slice_grid, float **slice, float **bottom,int xoption, int yoption, int voption)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int i,j,k,l,m,status;
  int negative=0,increasing=0,z_levels;
  int n0,n1;
  int compute_hmin,compute_hmax;

  double t,p;
  double tx,ty,tn,nx,ny;
  double *subx=NULL,*suby=NULL;

  float z,u,v,*vector2=NULL,*local=NULL,factor,level;
  float *physic=NULL;
  double *depths=NULL,h,href,*vector1=NULL,*levels=NULL;

  exitIfNull(
    vector1=(double *) malloc(grid.nz*sizeof(double))
    );
  exitIfNull(
    vector2=(float *)  malloc(grid.nz*sizeof(float))
    );

/*------------------------------------------------------------------------------
  initialize the section grid to empty grid*/

  slice_grid->nx=0;
  slice_grid->ny=0;

  slice_grid->xmin =1.e+10;
  slice_grid->xmax=-1.e+10;
/*------------------------------------------------------------------------------
  removed */

/*
  slice_grid->ymin =1.e+10;
  slice_grid->ymax=-1.e+10;
*/

  slice_grid->x= NULL;
  slice_grid->y= NULL;
  slice_grid->z= NULL;

/*------------------------------------------------------------------------------
  check the vertical arrangement */

  exitIfNull(
    depths=(double *) malloc(grid.nz*sizeof(double))
    );
  exitIfNull(
    levels=(double *) malloc(grid.nz*sizeof(double))
    );

/*------------------------------------------------------------------------------
  horizontal levels? */
  z_levels=1;
  for(k=0;k<grid.nz;k++) {
    for(m=0;m<grid.ny*grid.nx;m++) {
      h=grid.z[m+k*grid.nx*grid.ny];
      if (h!=grid.zmask) {
        href=h;
        levels[k]=h;
        break;
        }
      }
    for(m=m;m<grid.ny*grid.nx;m++) {
      h=grid.z[m+k*grid.nx*grid.ny];
      if ((h!=grid.zmask) && (h!=href)) {
        z_levels=0;
        k=grid.nz-1;
        break;
        }
      }
    }

/*------------------------------------------------------------------------------
  increasing/descending positive/negative levels? */
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      t=grid.x[m];
      p=grid.y[m];
      l=0;
      for(k=0;k<grid.nz;k++) {
        h=grid.z[m+k*grid.nx*grid.ny];
        if(h != grid.zmask) {
          depths[l]=grid.z[m+k*grid.nx*grid.ny];
          l++;
          }
        }
      if(l>1) {
        if (depths[0] < depths[1]) increasing=1;
        if(increasing==1) {
          if (depths[0] < 0.)      negative=1;
          }
        else
          {
          if (depths[1] < 0.)      negative=1;
          }
       goto next;
        }
      }
    }

/*------------------------------------------------------------------------------
  no valid depth available */

  status=-1;
  goto finished;

 next:

  free(depths);
/*------------------------------------------------------------------------------
  depth taken as negative for interpolation*/

  factor=1.;

  if(negative==0) {
    factor=-1.;
    increasing=1-increasing;
    }

  printf("negative= %d,increasing= %d,z_levels= %d\n",negative,increasing,z_levels);
/*------------------------------------------------------------------------------
  define the section grid*/
  slice_grid->ny=400;

  if(nloc==2) slice_grid->nx=400;
  else slice_grid->nx=nloc;

  exitIfNull(
    subx=(double *) malloc(slice_grid->nx*sizeof(double))
    );
  exitIfNull(
    suby=(double *) malloc(slice_grid->nx*sizeof(double))
    );

  tx=(x[nloc-1]-x[0])/(slice_grid->nx-1);
  ty=(y[nloc-1]-y[0])/(slice_grid->nx-1);
  tn=sqrt(tx*tx+ty*ty);

  switch (xoption)
    {
    case 0:
    if(x[nloc-1]<x[0]) {
      n0=nloc-1;
      n1=0;
      }
    else
      {
      n0=0;
      n1=nloc-1;
      }
    break;

    case 1:
    if(y[nloc-1]<y[0]) {
      n0=nloc-1;
      n1=0;
      }
    else
      {
      n0=0;
      n1=nloc-1;
      }
    break;

    default:
    n0=0;
    n1=nloc-1;
    break;
    }

  tx=(x[n1]-x[n0])/(slice_grid->nx-1);
  ty=(y[n1]-y[n0])/(slice_grid->nx-1);
  tn=sqrt(tx*tx+ty*ty);

  if(nloc==2) {
/*------------------------------------------------------------------------------
    sample a 2-position entry*/
    for(i=0;i<slice_grid->nx;i++) {
      t=x[n0]+i*tx;
      t=map_recale(grid,t);
      p=y[n0]+i*ty;
      subx[i]=t;
      suby[i]=p;
      }
    }
  else
    {
/*------------------------------------------------------------------------------
    use a multiple-position entry*/
    for(i=0;i<slice_grid->nx;i++) {
      subx[i]=x[i];
      suby[i]=y[i];
      }
    }

  exitIfNull(
    slice_grid->x=(double *) malloc(slice_grid->nx*slice_grid->ny*sizeof(double))
    );
  exitIfNull(
    slice_grid->y=(double *) malloc(slice_grid->nx*slice_grid->ny*sizeof(double))
    );
  slice_grid->z= NULL;

  exitIfNull(
    local=(float *) malloc(slice_grid->nx*slice_grid->ny*sizeof(float))
    );
  exitIfNull(
    *bottom=(float *) malloc(slice_grid->nx*sizeof(float))
    );
/*------------------------------------------------------------------------------
  commented, must be done by calling routine
  slice_grid->ymin =1.e+10;
  slice_grid->ymax=-1.e+10;
*/

  if (slice_grid->ymin== 1.e+10)
    compute_hmin=1;
  else
    compute_hmin=0;

  if (slice_grid->ymax==-1.e+10)
    compute_hmax=1;
  else
    compute_hmax=0;

/*------------------------------------------------------------------------------
  get depth min/max along the section */
  for(i=0;i<slice_grid->nx;i++) {
    (*bottom)[i]=1.e+5;
    for(k=0;k<grid.nz;k++) {
      m=grid.nx*grid.ny*k;
      depths=&(grid.z[m]);
      status= map_interpolation(grid,depths,grid.zmask,subx[i],suby[i], &h);
      physic=&(bufferx[m]);
      status= map_interpolation(grid,physic,mask,subx[i],suby[i], &u);
      physic=&(buffery[m]);
      status= map_interpolation(grid,physic,mask,subx[i],suby[i], &v);
      if((h != grid.zmask) && (u != mask) && (v != mask)) {
        if(compute_hmin==1) updatemin(&slice_grid->ymin,factor*h);
        if(compute_hmax==1) updatemax(&slice_grid->ymax,factor*h);
        updatemin(&(*bottom)[i],factor*h);
        }
      else
        {
          printf("i=%d k=%d h=%f u=%f v=%f mask=%f\n",i,k,h,u,v,mask);
        }
      }
    }

  slice_grid->dy=(slice_grid->ymax-slice_grid->ymin)/(slice_grid->ny-1);

/*------------------------------------------------------------------------------
  Build section axis */

  for(j=0;j<slice_grid->ny;j++)
    for(i=0;i<slice_grid->nx;i++) {
      m=i+slice_grid->nx*j;
      switch (xoption)
        {
/*------------------------------------------------------------------------------
        plot along longitudes*/
        case 0:
        slice_grid->x[m]=subx[i];
        break;

/*------------------------------------------------------------------------------
        plot along latitudes*/
        case 1:
        slice_grid->x[m]=suby[i];
        break;

/*------------------------------------------------------------------------------
        plot along distance in degrees*/
        case 2:
          if(i==0)
            tn=0;
          else
            {
            tx=subx[i]-subx[0];
            ty=suby[i]-suby[0];
            tn+=sqrt(tx*tx+ty*ty);
            }
          slice_grid->x[m]=tn;
/*         slice_grid->x[m]=i*tn;  */
        break;

/*------------------------------------------------------------------------------
        plot along distance in km*/
        case 3:
          if(i==0)
            tn=0;
          else
            tn+=geo_km(subx[i],suby[i],subx[i-1],suby[i-1]);
          slice_grid->x[m]=tn;
        break;
        }
      if(yoption==0) slice_grid->y[m]=slice_grid->ymin+j*slice_grid->dy;
      else           slice_grid->y[m]=j*slice_grid->dy;
/*
      if(z_levels==0) slice_grid->y[m]=slice_grid->ymin+j*slice_grid->dy;
      else            slice_grid->y[m]=levels[i];
*/
      }


/*------------------------------------------------------------------------------
  Interpolation*/
  for(i=0;i<slice_grid->nx;i++) {
    tx=subx[max(i+1,slice_grid->nx-1)]-subx[min(i-1,0)];
    ty=suby[max(i+1,slice_grid->nx-1)]-suby[min(i-1,0)];
    tn=sqrt(tx*tx+ty*ty);
    nx=-ty/tn;
    ny= tx/tn;
/*------------------------------------------------------------------------------
    first create a vertical profile of buffer at subx[i],suby[i]*/
    l=0;
    if(increasing==1) {
      for(k=0;k<grid.nz;k++) {
        m=grid.nx*grid.ny*k;
        depths=&(grid.z[m]);
        status= map_interpolation(grid,depths,grid.zmask,subx[i],suby[i], &h);
        if(h != grid.zmask) {
          vector1[l]=factor*h;
          physic=&(bufferx[m]);
          status= map_interpolation(grid,physic,mask,subx[i],suby[i], &u);
          physic=&(buffery[m]);
          status= map_interpolation(grid,physic,mask,subx[i],suby[i], &v);
          if((u != mask) && (v != mask)) {
            switch (voption)
              {
              case 0:
                vector2[l]=sqrt(u*u+v*v);
                break;
              case 1:
                vector2[l]=nx*u+ny*v;
                break;
              case 2:
                vector2[l]=ny*u-nx*v;
                break;
              }
            }
          else
            vector2[l]=mask;
          l++;
          }
        else
          {
          if(k==0) printf("%f %f\n",subx[i],suby[i]);
          }
        }
      }
    else
      {
      for(k=grid.nz-1;k>=0;k--) {
        m=grid.nx*grid.ny*k;
        depths=&(grid.z[m]);
        status= map_interpolation(grid,depths,grid.zmask,subx[i],suby[i], &h);
        if(h != grid.zmask) {
          vector1[l]=factor*h;
          physic=&(bufferx[m]);
          status= map_interpolation(grid,physic,mask,subx[i],suby[i], &u);
          physic=&(buffery[m]);
          status= map_interpolation(grid,physic,mask,subx[i],suby[i], &v);
          if((u != mask) && (v != mask)) {
            switch (voption)
              {
              case 0:
                vector2[l]=sqrt(u*u+v*v);
                break;
              case 1:
                vector2[l]=nx*u+ny*v;
                break;
              case 2:
                vector2[l]=ny*u-nx*v;
                break;
              }
            }
          else
            vector2[l]=mask;
          l++;
          }	
        }	
      }
/*------------------------------------------------------------------------------
    then interpolate on the slice profile*/
    if(l !=0) {
      for(j=0;j<slice_grid->ny;j++) {
        m=slice_grid->nx*j+i;
        if(yoption==0) level=slice_grid->y[m];
        else          level=slice_grid->y[m]+(*bottom)[i];
        status=map_interpolate1D(vector2, vector1, mask,l, level, &z);
        local[m]=z;
        }
      }
    else
      for(j=0;j<slice_grid->ny;j++) {
        m=slice_grid->nx*j+i;
        local[m]=mask;
        }
    }

  if(yoption!=0) {
    for(i=0;i<slice_grid->nx;i++) {
      (*bottom)[i]*=-1.;
      }
    }

  slice_grid->xmin=slice_grid->x[0];
  slice_grid->xmax=slice_grid->x[slice_grid->nx-1];
  slice_grid->dx=(slice_grid->xmax-slice_grid->xmin)/(slice_grid->nx-1);

  slice_grid->ymin=slice_grid->y[0];
  slice_grid->ymax=slice_grid->y[slice_grid->nx*slice_grid->ny-1];
  slice_grid->dy=(slice_grid->ymax-slice_grid->ymin)/(slice_grid->ny-1);

  slice_grid->zmin=0;
  slice_grid->zmax=0;
  slice_grid->dz=0;

  slice_grid->nz=1;

  slice_grid->modeH= 2;
  slice_grid->modeV=-1;

  slice_grid->circular=0;
  slice_grid->overlapped=0;
  slice_grid->connex=1;

  mapc_printgrid3d(*slice_grid);

  status=0;

 finished:

  if (levels   != NULL) free(levels);
  if (vector1  != NULL) free(vector1);
  if (vector2  != NULL) free(vector2);
  if (subx     != NULL) free(subx);
  if (suby     != NULL) free(suby);

  *slice=local;

  return(status);
}
