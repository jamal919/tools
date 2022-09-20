
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

/******************************************************************************/
/** \file

\brief various map utilities functions
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


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_get_format(const char *format)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int id;
  toponyms_t map_formats;
  
  if(format==0) return(-1);

//   map_formats["unknown"]=0;

  map_formats[MAP_FORMAT_NAME_GRD]        =MAP_FORMAT_CODE_GRD;
  map_formats[MAP_FORMAT_NAME_ASCII]      =MAP_FORMAT_CODE_ASCII;
  map_formats[MAP_FORMAT_NAME_ASCII_GIS]  =MAP_FORMAT_CODE_ASCII_GIS;
  map_formats[MAP_FORMAT_NAME_ASCII_SLIM] =MAP_FORMAT_CODE_ASCII_SLIM;

  int bkp=map_formats.size();
  
  id=map_formats[format];
  
  if(map_formats.size()!=bkp){
    id=-1;
    printf("available formats: %s %d\n",get_key_list(map_formats).c_str(), map_formats.size());
    }
//   printf("%s\n",get_key_list(map_formats).c_str());

  return(id);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  frame_t limits(const grid_t & grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  frame_t frame;
  
  frame=frame_t(grid.xmin,grid.xmax,grid.ymin,grid.ymax);
  
  return frame;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  frame_t limits(const grid_t *grids,int ngrids)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k;
  frame_t frame;
  
  frame=limits(grids[0]);
  
  for(k=1;k<ngrids;k++){
    double x;
    x=degree_recale(grids[k].xmin, grids[0].xmin);
    updatemin(&frame.xmin,x);
    x=degree_recale(grids[k].xmax, grids[0].xmax);
    updatemax(&frame.xmax,x);
    updatemin(&frame.ymin,grids[k].ymin);
    updatemax(&frame.ymax,grids[k].ymax);
    }
  
  return(frame);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> T map_recale_template(const grid_t &grid,T x)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T xmid,shifted;

//  xmid=grid.xmin+(T)180.;
  xmid=0.5*(grid.xmin+grid.xmax);
  shifted=degree_recale(x,xmid);

  return(shifted);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float map_recale(const grid_t &grid,float x)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float shifted;

  shifted=map_recale_template(grid,x);

  return(shifted);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double map_recale(const grid_t &grid,double x)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double shifted;

  shifted=map_recale_template(grid,x);

  return(shifted);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void map_recale_all_x(grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// recale all x coordinates with previous one
/**
\param[in,out] *grid
*/
/*----------------------------------------------------------------------------*/
{
  int i,j,m,n=-1;
  
  switch(grid->modeH){
  case 0:
    return;
  case 1:
    n=1;
    break;
  case 2:
    n=grid->ny;
    break;
  case MODEH_REDUCED_GG:
    return;
  default:
    TRAP_ERR_EXIT(ENOEXEC,"not coded yet for modeH=%d\n",grid->modeH);
    }
  
  for(j=0;j<n;j++){
    double *xm;
    
    m=j*grid->nx;
    
    for(i=0;i<grid->nx;i++,m++){
      
      xm=&grid->x[m];
      
      if(i==0){
        
        if(j>1){
          degree_recale(xm,grid->x[m-grid->nx]);
          }
        
        continue;
        }
      
      degree_recale(xm,grid->x[m-1]);
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_resolution00(const grid_t &grid, double *dx, double *dy, double *dz)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,m,m0,m1;

  for(j=1;j<grid.ny;j++) {
    m=j;
    m0=j-1;
    m1=j+1;
    dy[m]=0.5*(grid.y[m1]-grid.y[m0]);
    }

  for(i=1;i<grid.nx;i++) {
    m=i;
    m0=i-1;
    m1=i+i;
    dx[m]=0.5*(grid.x[m1]-grid.x[m0]);
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double map_cartesian_distance(const grid_t &grid, int m0, int m1, double step)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double d, dx, dy;
  double x0,x1,y0,y1;
  
  dx=grid.x[m1]-grid.x[m0];
  dy=grid.y[m1]-grid.y[m0];
  
  d=hypot(dx,dy);
  d/=step;
  
  return(d);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double map_spherical_distance(const grid_t &grid, int m0, int m1, double step)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double t1,p1,t2,p2;
  double d;
  
  t2=grid.x[m1];
  t1=grid.x[m0];
  
  p2=grid.y[m1];
  p1=grid.y[m0];
  
  d=geo_haversin(t1,p1,t2,p2);
  d/=step;
  
  return(d);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_resolution02(const grid_t &grid, double *dx, double *dy, double *dz, double (*map_distance)(const grid_t &, int, int, double))

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,l,m,m0,m1;
  
//  double (*map_distance)(const grid_t &, int, int, double);
  
//   map_distance=map_cartesian_distance;

  for(j=1;j<grid.ny-1;j++) {
    for(i=1;i<grid.nx-1;i++) {
      m=grid.nx*j+i;
      m0=grid.nx*j+i-1;
      m1=grid.nx*j+i+1;
      dx[m]=map_distance(grid, m0, m1, 2.0);
      m0=grid.nx*(j-1)+i;
      m1=grid.nx*(j+1)+i;
      dy[m]=map_distance(grid, m0, m1, 2.0);
      }
    }

  i=0;
  for(j=1;j<grid.ny-1;j++) {
    m=grid.nx*j+i;
    m0=grid.nx*j+i;
    m1=grid.nx*j+i+1;
    dx[m]=map_distance(grid, m0, m1, 1.0);
    m0=grid.nx*(j-1)+i;
    m1=grid.nx*(j+1)+i;
    dy[m]=map_distance(grid, m0, m1, 2.0);
    }

  i=grid.nx-1;
  for(j=1;j<grid.ny-1;j++) {
    m=grid.nx*j+i;
    m0=grid.nx*j+i-1;
    m1=grid.nx*j+i;
    dx[m]=map_distance(grid, m0, m1, 1.0);
    m0=grid.nx*(j-1)+i;
    m1=grid.nx*(j+1)+i;
    dy[m]=map_distance(grid, m0, m1, 2.0);
    }

  j=0;
  for(i=1;i<grid.nx-1;i++) {
    m=grid.nx*j+i;
    m0=grid.nx*j+i-1;
    m1=grid.nx*j+i+1;
    dx[m]=map_distance(grid, m0, m1, 2.0);
    m0=grid.nx*j+i;
    m1=grid.nx*(j+1)+i;
    dy[m]=map_distance(grid, m0, m1, 1.0);
    }

  j=grid.ny-1;
  for(i=1;i<grid.nx-1;i++) {
    m=grid.nx*j+i;
    m0=grid.nx*j+i-1;
    m1=grid.nx*j+i+1;
    dx[m]=map_distance(grid, m0, m1, 2.0);
    m0=grid.nx*(j-1)+i;
    m1=grid.nx*j+i;
    dy[m]=map_distance(grid, m0, m1, 1.0);
    }

  i=0;
  j=0;
  m=grid.nx*j+i;
  m0=grid.nx*j+i;
  m1=grid.nx*j+i+1;
  dx[m]=map_distance(grid, m0, m1, 1.0);
  m0=grid.nx*j+i;
  m1=grid.nx*(j+1)+i;
  dy[m]=map_distance(grid, m0, m1, 1.0);

  i=grid.nx-1;
  j=0;
  m=grid.nx*j+i;
  m0=grid.nx*j+i-1;
  m1=grid.nx*j+i;
  dx[m]=map_distance(grid, m0, m1, 1.0);
  m0=grid.nx*j+i;
  m1=grid.nx*(j+1)+i;
  dy[m]=map_distance(grid, m0, m1, 1.0);

  i=0;
  j=grid.ny-1;
  m=grid.nx*j+i;
  m0=grid.nx*j+i;
  m1=grid.nx*j+i+1;
  dx[m]=map_distance(grid, m0, m1, 1.0);
  m0=grid.nx*(j-1)+i;
  m1=grid.nx*j+i;
  dy[m]=map_distance(grid, m0, m1, 1.0);

  i=grid.nx-1;
  j=grid.ny-1;
  m=grid.nx*j+i;
  m0=grid.nx*j+i-1;
  m1=grid.nx*j+i;
  dx[m]=map_distance(grid, m0, m1, 1.0);
  m0=grid.nx*(j-1)+i;
  m1=grid.nx*j+i;
  dy[m]=map_distance(grid, m0, m1, 1.0);

  if(grid.z==0) return(0);
  if(dz==0) return(0);

  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      for(l=1;l<grid.nz-1;l++) {
        m=grid.ny*grid.nx*l+grid.nx*j+i;
        m0=grid.ny*grid.nx*(l-1)+grid.nx*j+i;
        m1=grid.ny*grid.nx*(l+1)+grid.nx*j+i;
        dz[m]=0.5*(grid.z[m1]-grid.z[m0]);
//        dz[m]=0.0;
        }
      l=0;
      m=grid.ny*grid.nx*l+grid.nx*j+i;
      m0=grid.ny*grid.nx*l+grid.nx*j+i;
      m1=grid.ny*grid.nx*(l+1)+grid.nx*j+i;
      dz[m]=grid.z[m1]-grid.z[m0];
//      dz[m]=0.0;
      l=grid.nz-1;
      m=grid.ny*grid.nx*l+grid.nx*j+i;
      m0=grid.ny*grid.nx*(l-1)+grid.nx*j+i;
      m1=grid.ny*grid.nx*l+grid.nx*j+i;
      dz[m]=grid.z[m1]-grid.z[m0];
//      dz[m]=0.0;
      }
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_cartesian_resolution(const grid_t &grid, double **dx, double **dy, double **dz)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  const size_t size=grid.Hsize();
  double *dz0;
  
  if(dz!=0){
    *dz=new double[size*grid.ny];
    dz0=*dz;
    }
  else
    dz0=NULL;

  switch (grid.modeH) {
    case 0:
      status=-1;
      break;

    case 1:
      *dx=new double[grid.nx];
      *dy=new double[grid.ny];
      status=map_resolution00(grid, *dx, *dy, dz0);
      break;

    case 2:
      *dx=new double[size];
      *dy=new double[size];
      status=map_resolution02(grid, *dx, *dy, dz0, map_cartesian_distance);
      break;

    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_spherical_resolution(const grid_t &grid, double **dx, double **dy, double **dz)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  const size_t size=grid.Hsize();
  double *dz0;
  
  if(dz!=0){
    *dz=new double[size*grid.ny];
    dz0=*dz;
    }
  else
    dz0=NULL;

  switch (grid.modeH) {
    case 0:
      status=-1;
      break;

    case 1:
      status=-1;
      break;

    case 2:
      *dx=new double[size];
      *dy=new double[size];
      status=map_resolution02(grid, *dx, *dy, dz0, map_spherical_distance);
      break;

    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double map_resolution(const grid_t &grid, double t, double p, double *d, double *dy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,status __attribute__((unused));
  double rmin/*,r[4],c*/;
  double x1,x2,y1,y2;
  int map_index(const grid_t &grid, double x, double y, int *ktrue, int *ltrue);

  status=map_index(grid,t,p,&k,&l);

//   c=cos(p*d2r);
  
  if(k==0) k=1;
  if(k==grid.nx-1) k=grid.nx-2;

  if(l==0) l=1;
  if(l==grid.ny-1) l=grid.ny-2;

  x1=map_grid_x(grid,k-1,l);
  x2=map_grid_x(grid,k+1,l);

  y1=map_grid_y(grid,k,l-1);
  y2=map_grid_y(grid,k,l+1);

//   r[0]=geo_distance(x1,y1,x2,y1,'m');
//   r[1]=geo_distance(x1,y2,x2,y2,'m');
// 
//   r[2]=geo_distance(x1,y1,x1,y2,'m');
//   r[3]=geo_distance(x2,y1,x2,y2,'m');

  rmin=1000.0*geo_distance(x1,y1,x2,y2)/sqrt(2.0);

  switch (grid.modeH) {
    case 0:
      break;

    case 1:
      break;

    case 2:
      break;

    }
  if(rmin==0.) {
    printf("%s : trouble\n",__func__);
    }

  return(rmin);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double map_resolution(const grid_t & grid, int k, int l, double & dx, double & dy, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  give resolution of a grid
   
  mode: if 0, use geo_distance() with m. if 1, use hypot().
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double rmin;
  double x1,x2,y1,y2;
  
  if(k==0) k=1;
  if(k==grid.nx-1) k=grid.nx-2;

  if(l==0) l=1;
  if(l==grid.ny-1) l=grid.ny-2;
  
  grid.xy(k-1,l,x1,y1);
  grid.xy(k+1,l,x2,y2);
  
  switch (mode) {
    case 0:
      dx=geo_distance(x1,y1,x2,y2,'m');
      break;
    
    case 1:
      dx=hypot(x2-x1,y2-y1);
      break;
    
    }
  
  grid.xy(k,l-1,x1,y1);
  grid.xy(k,l+1,x2,y2);
  
  switch (mode) {
    case 0:
      dy=geo_distance(x1,y1,x2,y2,'m');
      break;
    
    case 1:
      dy=hypot(x2-x1,y2-y1);
      break;
    
    }
  
  rmin=sqrt(dx*dx+dy*dy)/sqrt(2.0);

  if(rmin==0.) {
    printf("trouble\n");
    }

  return(rmin);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mapc_checkconnexity02(const grid_t &grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//   double local=0.;
  int connex=1;

  const double dmax=500000.;

#pragma omp parallel for
  for(int l=0;l<grid.ny;l++) {
    double ux,uy,uz;
    double a1,b1;
    int m;
    m=grid.nx*l;
    a1=grid.x[m]*d2r;
    b1=grid.y[m]*d2r;
    ux=cos(a1)*cos(b1);
    uy=sin(a1)*cos(b1);
    uz=sin(b1);
    for(int k=0;k<grid.nx-1;k++) {
      double vx,vy,vz;
      double a2,b2,d;
      m=grid.nx*l+k+1;
      a2=grid.x[m]*d2r;
      b2=grid.y[m]*d2r;
      vx=cos(a2)*cos(b2);
      vy=sin(a2)*cos(b2);
      vz=sin(b2);
      d=acos(ux*vx+uy*vy+uz*vz)*MeanEarthRadius;
      ux=vx;
      uy=vy;
      uz=vz;
//       updatemax(&local,d);
      if(d>dmax) {
#pragma omp critical(not_connex)
          {
          connex=0;
          }
        }
      }
    }
  
  return(connex);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mapc_checkconnexity01(const grid_t &grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    l;
  double a1,a2=0.,b1,b2,d,local=0.;
  double ux,uy,uz,vx=NAN,vy=NAN,vz=NAN;

  const double dmax=500000.;

  a1=grid.x[0]*d2r;
  
  for(l=0;l<grid.ny;l++) {
    b1=grid.y[l]*d2r;
    ux=vx;
    uy=vy;
    uz=vz;
    vx=cos(a1)*cos(b1);
    vy=sin(a1)*cos(b1);
    vz=sin(b1);
    if(l==0)continue;
    d=acos(ux*vx+uy*vy+uz*vz)*MeanEarthRadius;
    local=max(local,d);
    if(d>dmax) {
      return(0);
      }
    b2=min(fabs(b1),b2);
    }
  
  for(l=0;l<grid.nx;l++) {
    a2=grid.x[l]*d2r;
    ux=vx;
    uy=vy;
    uz=vz;
    vx=cos(a2)*cos(b2);
    vy=sin(a2)*cos(b2);
    vz=sin(b2);
    if(l==0)continue;
    d=acos(ux*vx+uy*vy+uz*vz)*MeanEarthRadius;
    local=max(local,d);
    if(d>dmax) {
      return(0);
      }
    }
  
  return(1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mapc_checkconnexity(const grid_t &grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  switch(grid.modeH){
  case 0:
    status=1;
    break;
  case 1:
    status=mapc_checkconnexity01(grid);
    break;
  case 2:
    status=mapc_checkconnexity02(grid);
    break;
  default:
    TRAP_ERR_EXIT(ENOEXEC,"modeH=%d\n",grid.modeH);
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void next_m (const grid_t & grid,int *m,int i0,int im,int j0,int jm,bool checkNan)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// next valid border grid index, matrix-clockwise
/*----------------------------------------------------------------------------*/
{
  int i,j,border=0;
#define BORDER_I0 1
#define BORDER_IM 2
#define BORDER_J0 4
#define BORDER_JM 8
  
  grid.ij(*m,&i,&j);
  
  if(i==i0)
    border|=BORDER_I0;
  
  if(i==im)
    border|=BORDER_IM;
  
  if(j==j0)
    border|=BORDER_J0;
  
  if(j==jm)
    border|=BORDER_JM;
  
  if(border==0 and not checkNan)
    TRAP_ERR_EXIT(ENOEXEC,"%s:programming error\n",__func__);
  
  while(checkNan){
    if(isnan(grid.x[*m]) and isnan(grid.y[*m]))
      TRAP_ERR_EXIT(ENOEXEC,"%s:programming error [%d]\n",__func__,*m);
    
    int m1;
    
    m1=*m-1;
    if(i>i0 and isnan(grid.x[m1]) and isnan(grid.y[m1]))
      border|=BORDER_I0;
    
    m1=*m+1;
    if(i<im and isnan(grid.x[m1]) and isnan(grid.y[m1]))
      border|=BORDER_IM;
    
    m1=*m-grid.nx;
    if(j>j0 and isnan(grid.x[m1]) and isnan(grid.y[m1]))
      border|=BORDER_J0;
    
    m1=*m+grid.nx;
    if(j<jm and isnan(grid.x[m1]) and isnan(grid.y[m1]))
      border|=BORDER_JM;
    
    if(border!=0)
      break;
    
    /* inner corners : fake equivalent border */
    
    m1=*m-grid.nx-1;
    if(i>i0 and j>j0 and isnan(grid.x[m1]) and isnan(grid.y[m1]))
      border|=BORDER_I0;
    
    m1=*m-grid.nx+1;
    if(i<im and j>j0 and isnan(grid.x[m1]) and isnan(grid.y[m1]))
      border|=BORDER_J0;
    
    m1=*m+grid.nx+1;
    if(i<im and j<jm and isnan(grid.x[m1]) and isnan(grid.y[m1]))
      border|=BORDER_IM;
    
    m1=*m+grid.nx-1;
    if(i>i0 and j<jm and isnan(grid.x[m1]) and isnan(grid.y[m1]))
      border|=BORDER_JM;
    
    if(border!=BORDER_I0 and border!=BORDER_IM and
       border!=BORDER_J0 and border!=BORDER_JM)
      TRAP_ERR_EXIT(ENOEXEC,"%s:programming error\n",__func__);
    
    break;
    }
  
  if(border & BORDER_J0 and border & BORDER_JM)
    TRAP_ERR_EXIT(ENOEXEC,"%s:programming error: %d\n",__func__,border);
  if(border & BORDER_I0 and border & BORDER_IM)
    TRAP_ERR_EXIT(ENOEXEC,"%s:programming error: %d\n",__func__,border);
  
  if(border & BORDER_J0 and not (border & BORDER_IM) ){
    *m+=1;
    return;
    }
  
  if(border & BORDER_IM and not (border & BORDER_JM) ){
    *m+=grid.nx;
    return;
    }
  
  if(border & BORDER_JM and not (border & BORDER_I0) ){
    *m-=1;
    return;
    }
  
  if(border & BORDER_I0 and not (border & BORDER_J0) ){
    *m-=grid.nx;
    return;
    }
  
  TRAP_ERR_EXIT(ENOEXEC,"%s:programming error: [%d]%d\n",__func__,*m,border);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *map_get_perimeter_indexes(const grid_t & grid,int *n,const double *buffer,double mask,int detectMaskedLine,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// get grid perimeter indexes
/**
\param grid
\param *n set to the number of indexes
\param *buffer
\param mask
\param detectMaskedLine detect bogus, and sadly most current, masked line on one or several border of the field
\return array[*n] of indexes
*/
/*----------------------------------------------------------------------------*/
{
  int n_,m;
  if(n==0)
    n=&n_;
  
  bool i0m=false,inm=false,j0m=false,jnm=false; /*< whether any of the lines are masked */
  
  if(detectMaskedLine!=0){
    if(verbose) STDERR_BASE_LINE_FUNC("detecting masked lines...");
    
    int i;
    
    /* first column */
    m=0;
    for(i=0;i<grid.ny;i++,m+=grid.nx)
      if(buffer[m]!=mask) break;
    i0m= i>=grid.ny;
    
    /* last column */
    m=grid.nx-1;
    for(i=0;i<grid.ny;i++,m+=grid.nx)
      if(buffer[m]!=mask) break;
    inm= i>=grid.ny;
    
    /* first row */
    m=0;
    for(i=0;i<grid.nx;i++,m++)
      if(buffer[m]!=mask) break;
    j0m= i>=grid.nx;
    
    /* last row */
    m=(grid.ny-1)*grid.nx;
    for(i=0;i<grid.nx;i++,m++)
      if(buffer[m]!=mask) break;
    jnm= i>=grid.nx;
    
    if(verbose) fprintf(stderr," i:%d%d j:%d%d\n",i0m,inm,j0m,jnm);
    }
  
  const int
    i0=i0m,
    imax=grid.nx-inm-1,
    ni=imax-i0,
    j0=j0m,
    jmax=grid.ny-jnm-1,
    nj=jmax-j0;
  int
    m0=grid.nx*j0+i0;
  
  vector<int> vi;
  
  *n=2*(ni+nj)-4;
  
  m=m0;
  if(grid.modeH==2){
    do{
      if(isfinite(grid.x[m]) and isfinite(grid.y[m]))
        break;
      
      next_m(grid,&m,i0,imax,j0,jmax,false);
      
      }while(m!=m0);
    m0=m;
    }
  
  do{
    vi.push_back(m);
    
    if(m==3846)
    STDERR_BASE_LINE_FUNC("%d\n",m);
    next_m(grid,&m,i0,imax,j0,jmax,true);
    
    }while(m!=m0 and vi.size()<2**n);
  
  int *indexes=0;
  
  *n=vi.size();
  indexes=copy(vi);
  
  return indexes;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

grid_t map_getgrid2d(const grid_t &grid3d)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  return(grid3d);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

grid_t map_getgrid3d(const grid_t &grid2d)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  grid_t grid;

  grid.xmin=grid2d.xmin;
  grid.xmax=grid2d.xmax;

  grid.ymin=grid2d.ymin;
  grid.ymax=grid2d.ymax;

  grid.zmin=0.0;
  grid.zmax=0.0;

  grid.nx=grid2d.nx;
  grid.ny=grid2d.ny;
  grid.nz=1;

  grid.dx=grid2d.dx;
  grid.dy=grid2d.dy;
  grid.dz=0.0;

  grid.x=grid2d.x;
  grid.y=grid2d.y;
  grid.z=NULL;

  grid.mask=0;

  grid.modeH=grid2d.modeH;
  grid.modeV=-1;
/*
  grid.connex=grid2d.connex;
  grid.overlapped=grid2d.overlapped;
  grid.circular=grid2d.circular;
*/
  return(grid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_set_dxdy(grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Initialises grid_t::dx and grid_t::dy
/**
\returns 0
*/
/*----------------------------------------------------------------------------*/
{
  int status=0;
  
  /** \sa map_set2Dgrid() */
  grid->dx=(grid->xmax-grid->xmin)/(grid->nx-1);
  grid->dy=(grid->ymax-grid->ymin)/(grid->ny-1);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void map_allocate_x_y(grid_t *grid,int modeH,int nx,int ny)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///allocates x and y
/*----------------------------------------------------------------------------*/
{
  if(modeH>=0)
    grid->modeH=modeH;
  
  if(nx>0)
    grid->nx=nx;
  if(ny>0)
    grid->ny=ny;
  
  switch (grid->modeH){
  
  case 0:
    break;
  
  case 1:
    grid->x=new double[grid->nx];
    grid->y=new double[grid->ny];
    break;
  
  case 2:
    int n=grid->Hsize();
    grid->x=new double[n];
    grid->y=new double[n];
    break;
  
  }
  
}


/*----------------------------------------------------------------------------*/
///Initialises x and y as vectors
/**
\param *grid grid to complete
\returns 0
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_completegridaxis_1(grid_t *grid, int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,status=0;

  switch (grid->modeH) {
/*----------------------------------------------------------------------------*/
    case 0:
      map_allocate_x_y(grid,1);
      for(j=0;j<grid->ny;j++)
        grid->y[j]=grid->ymin+j*grid->dy;
      for(i=0;i<grid->nx;i++)
        grid->x[i]=grid->xmin+i*grid->dx;
      break;

/*----------------------------------------------------------------------------*/
    case 1:
      break;

/*----------------------------------------------------------------------------*/
    case 2:{
      bool compressible=true;
      int c,m,m0,mj;
      
/*------------------------------------------------------------------------------
      test compressibility */
      for(c=0;c<2 and compressible;c++){
        const double *data=0;
        switch(c){
        case 0: data=grid->x; break;
        case 1: data=grid->y; break;
          }
        
        m=grid->nx;
        for(j=1;j<grid->ny and compressible;j++){
          mj=m;
          
          for(i=0;i<grid->nx;i++,m++){
            
            if(c==0)
              m0=i;
            else
              m0=mj;
            
            if(data[m0]==data[m]) continue;
            
            compressible=false;
            break;
            }
          }
        
        }
      
/*------------------------------------------------------------------------------
      compress */
      if(not compressible)
        return -1;
      
      if(verbose>0) STDERR_BASE_LINE_FUNC("COMODO-COMPRESSING %dx%d grid\n",grid->nx,grid->ny);
      double
        *x=grid->x,
        *y=grid->y;
      map_allocate_x_y(grid,1);
      for(j=0,m=0;j<grid->ny;j++,m+=grid->nx)
        grid->y[j]=y[m];
      for(i=0;i<grid->nx;i++)
        grid->x[i]=x[i];
      delete [] x;
      delete [] y;
      
      }break;

    }

  grid->modeH=1;
  return(status);
}


/*----------------------------------------------------------------------------*/
///Initialises x and y as vectors or matrices
/**
\date reviewed 25 Jul 2011
\author Damien Allain

\param *grid grid to complete
\param mode If 1, initialises as vectors. If 2(default), initialises as matrices.
\returns 0 or -1 if nothing needed to be done (which is not a great deal)
Calls either map_completegridaxis_1() or map_completegridaxis_2().
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_completegridaxis(grid_t *grid, int target_mode, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;

  switch (target_mode) {
    case 0:
      if(grid->modeH!=0) return(status);
      status=map_completegridaxis_1(grid);
      break;

    case 1:
      status=map_completegridaxis_1(grid,verbose);
      break;

    case 2:
      status=map_completegridaxis_2(grid);
      break;

    }

  return(status);
}

/*----------------------------------------------------------------------------*/
///Sets xmin,ymin,xmax,ymax,dx and dy.
/**
\date created 25 Jul 2011
\author Damien Allain

\param *grid grid to set
\param xmin
\param ymin
\param xmax
\param ymax
\param dx Default:0.
\param dy Default:dx
\returns 0 if success, 1 if error (e.g. if dx or dy are 0. or negative)
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_set2Dgrid(grid_t *grid, double xmin, double ymin, double xmax, double ymax, double dx, double dy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  
  if(dy==0.){
    dy=dx;
    }
  grid->xmin=xmin;
  grid->ymin=ymin;
  grid->xmax=xmax;
  grid->ymax=ymax;
  grid->dx=dx;
  grid->dy=dy;

  if(grid->dx<=0. || grid->dy<=0.)
  ///If dx or dy are 0. or negative, neither nx nor ny will be calculated.
    return 1;
    
  ///nx and ny are otherwise calculated according to the other parameters.
  grid->nx  = (int)( (grid->xmax-grid->xmin)/grid->dx )+1;
  grid->ny  = (int)( (grid->ymax-grid->ymin)/grid->dy )+1;
  
  grid->modeH=0;

  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_minmax(grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///set xmin,ymin,xmax,ymax
/**
\return 0 if success, 1 if error (i.e. never, because it #TRAP_ERR_EXIT's instead)
*/
/*----------------------------------------------------------------------------*/
{
  switch (grid->modeH){
  
  case 1:
    poc_minmax(grid->x,grid->nx,NAN,&grid->xmin,&grid->xmax);
    poc_minmax(grid->y,grid->ny,NAN,&grid->ymin,&grid->ymax);
    break;
  
  case 2:{
    const int n=grid->Hsize();
    poc_minmax(grid->x,n,NAN,&grid->xmin,&grid->xmax);
    poc_minmax(grid->y,n,NAN,&grid->ymin,&grid->ymax);
    }break;
  
  default:
    TRAP_ERR_EXIT(ENOEXEC,"%s should not have been called with grid_t::modeH=%d\n",__func__,grid->modeH);
    }
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_set2Dgrid(grid_t *grid, frame_t frame, double dx, double dy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=map_set2Dgrid(grid, frame.xmin, frame.ymin, frame.xmax, frame.ymax, dx, dy);

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_set2Dgrid(grid_t *grid, frame_t frame, int nx, int ny)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  double dx, dy;
  
  dx=(frame.xmax-frame.xmin)/(double) nx;
  dy=(frame.ymax-frame.ymin)/(double) ny;

  status=map_set2Dgrid(grid, frame.xmin, frame.ymin, frame.xmax, frame.ymax, dx, dy);

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_set2Dgrid(grid_t *grid, frame_t frame, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  double dx, dy;
  
  dx=(frame.xmax-frame.xmin)/(double) n;
  dy=(frame.ymax-frame.ymin)/(double) n;
  
  if(dx>dy) dx=dy;
  else dy=dx;

  status=map_set2Dgrid(grid, frame.xmin, frame.ymin, frame.xmax, frame.ymax, dx, dy);

  return status;
}


/*----------------------------------------------------------------------------*/
///Initialises x and y as matrices
/**
\param *grid grid to complete
\returns 0
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_completegridaxis_2(grid_t *grid, int transposed)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,m,status=0;

  switch (grid->modeH) {
    case 0:
      map_allocate_x_y(grid,2);
      for(j=0;j<grid->ny;j++) {
        for(i=0;i<grid->nx;i++) {
          if(transposed)
            m=j+i*grid->ny;
          else
            m=j*grid->nx+i;
          grid->x[m]=grid->xmin+i*grid->dx;
          grid->y[m]=grid->ymin+j*grid->dy;
          }
        }
      break;

    case 1:{
      double
        *x=grid->x,
        *y=grid->y;
      map_allocate_x_y(grid,2);
      for(j=0;j<grid->ny;j++) {
        for(i=0;i<grid->nx;i++) {
          if(transposed)
            m=j+i*grid->ny;
          else
            m=j*grid->nx+i;
          grid->x[m]=x[i];
          grid->y[m]=y[j];
          }
        }
      delete [] x;
      delete [] y;
      }break;

    case 2:
      if(grid->x==0 && grid->y==0) {
        status=-1;
        }
      break;

    }
  
  if(transposed){
    swapValues(&grid->nx,&grid->ny);
    grid->modeH=-2;
    }
  else
    grid->modeH=2;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map3d_completegridaxis(grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;

  switch (grid->modeH) {
    case 0:
      status=map_completegridaxis_2(grid);
      break;

    case 1:
    case 2:
      break;

    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void map_printgrid_head (const grid_t &grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double dx=NAN,dy=NAN;
  bool approximated=false;
  
  switch (grid.modeH) {
    case 0:
      dx=grid.dx;
      dy=grid.dy;
      break;

    case 1:
    case 2:
      if(grid.dx==0 or isnan(grid.dx)) {
        dx=(grid.xmax-grid.xmin)/(grid.nx-1);
        approximated=true;
        }
      else dx=grid.dx;
      if(grid.dy==0 or isnan(grid.dy)) {
        dy=(grid.ymax-grid.ymin)/(grid.ny-1);
        approximated=true;
        }
      else dy=grid.dy;
      break;

    }
  printf("xmin = %f, xmax= %f\n",grid.xmin,grid.xmax);
  printf("ymin = %f, ymax= %f\n",grid.ymin,grid.ymax);
  printf("dx = %lf, dy= %lf",dx,dy);
  if(approximated) {
    printf(" (approximated)\n");
    }
  else {
    printf("\n");
    }
  printf("nx = %d, ny= %d",grid.nx,grid.ny);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void map_printgrid (const grid_t &grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  map_printgrid_head(grid);
  printf(", mode= %d\n",grid.modeH);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void mapc_printgrid3d (const grid_t &grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  map_printgrid_head(grid);
  printf("\n");
  printf("modeH= %d, modeV= %d\n",grid.modeH,grid.modeV);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int wheader(FILE *out,char **comment, grid_t grid, float spec, int n, float *average)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   M;
  size_t nitems __attribute__((unused)),size;

  size=256;
  nitems=fwrite(comment[0],size,1,out);
  nitems=fwrite(comment[1],size,1,out);
  nitems=fwrite(comment[2],size,1,out);
  nitems=fwrite(comment[3],size,1,out);

  size=sizeof(grid);
  nitems=fwrite(&grid,size,1,out);
  size=sizeof(float);
  nitems=fwrite(&spec,size,1,out);
  
  M=grid.nx*grid.ny;
  size=sizeof(int);
  nitems=fwrite(&M,size,1,out);
  size=sizeof(int);
  nitems=fwrite(&n,size,1,out);
  size=n*sizeof(float);
  nitems=fwrite(average,size,1,out);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int rheader( FILE *out, grid_t *grid, float *spec, int *M, int *n, float *average)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  recordsize;
  size_t nitems __attribute__((unused)),size;
  char comment[256];

  recordsize=256;
  nitems=fread(comment,recordsize,1,out);
  printf("%s\n",comment);
  nitems=fread(comment,recordsize,1,out);
  printf("%s\n",comment);
  nitems=fread(comment,recordsize,1,out);
  printf("%s\n",comment);
  nitems=fread(comment,recordsize,1,out);
  printf("%s\n",comment);

  recordsize=sizeof(grid);
  nitems=fread(grid,recordsize,1,out);
  recordsize=sizeof(float);
  nitems=fread(spec,recordsize,1,out);
  
  recordsize=sizeof(int);
  nitems=fread(M,recordsize,1,out);
  recordsize=sizeof(int);
  nitems=fread(n,recordsize,1,out);
  size=(*n)*sizeof(float);
  exitIfNull(
    average=(float *)malloc(size)
    );
  nitems=fread(average,size,1,out);
  return(0);
  
}


#if 0
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool map_IsRectangular(grid_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,status;

  switch (grid.modeH) {
    case 0:
      return(true);
      break;

    case 1:
      for(j=1;j<grid.ny;j++) {
        if(grid.y[j]!=grid.y[0]) {
          return(false);
          }
        }
      for(i=1;i<grid.nx;i++) {
        if(grid.x[i]!=grid.x[0]) {
          return(false);
          }
        }
      return(true);
      break;

    case 2:
      double x,y;
      y=grid.y[0];
      for(j=1;j<grid.ny;j++) {
        size_t m=j*grid.nx;
        if(grid.y[m]!=grid.y[0]) {
          return(false);
          }
        }
      y=grid.y[grid.nx-1];
      for(j=1;j<grid.ny;j++) {
        size_t m=j*grid.nx+grid.nx-1;
        if(grid.y[m]!=y) {
          return(false);
          }
        }
      x=grid.x[0];
      for(i=1;i<grid.nx;i++) {
        size_t m=i;
        if(grid.x[i]!=x) {
          return(false);
          }
        }
      x=grid.x[(grid.ny-1)*grid.nx];
      for(i=1;i<grid.nx;i++) {
        size_t m=(grid.ny-1)*grid.nx+i;
#error HERE!!!
        if(grid.x[i]!=x) {
          return(false);
          }
        }
      return(true);
      break;

    }
  
  
  return(false);
}
#endif
