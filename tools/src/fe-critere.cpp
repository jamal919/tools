
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\date 2D mesh generator functions
*/
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <unistd.h> //unlink

#include "version-macros.def" //for VERSION and REVISION

#include "tools-structures.h"

#include "functions.h"
#include "constants.h"
#include "geo.h"
#include "fe.h"
#include "archive.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "poc-netcdf-data.hpp"
#include "list.h"
#include "maths.h"
#include "map.def"
#include "topo.h"

#include "zapper.h"     /*  rutin.h contains common utility routines  */

#define MSH_MAXSIZE_OPEN     1
#define MSH_MAXSIZE_SHELF    2
#define MSH_TIDAL_WAVELENGTH 3
#define MSH_TOPO_SLOPE       4
#define MSH_SURFWAVE_LENGTH  5
#define MSH_SURFWAVE_DZ      6
#define MSH_SURFWAVE_CFL     7


#include "statistic.h"


const double km2m=1e3;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_defgrid(grid_t *grid, vector<plg_t> & polygons, vector<plg_t> & limits, double resolution)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  return a cartesian grid at the required resolution (in meters)
  
  assessing proper cartesian min/max setting is a recurrent issue in polar regions
  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  frame_t frame;
  projPJ projection;
  double ref_lat,ref_lon,x,y;
  char *pj_parameters=new char[512];
  bool debug=false, subdivide;;

  frame=plg_recale(polygons);
  
  bool equalize=false;
  double radius=1.0;
  
  subdivide=true;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  tackling polar regions issues (orthodromic versus loxodromic)
  
  flag test to be revised (class method?)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(subdivide) {
    for(int s=0;s<polygons.size();s++) {
      plg_t *tmp=new plg_t;
      if(polygons[s].flag==0) {
        tmp->duplicate(polygons[s]);
        limits.push_back(*tmp);
        }
      else if(polygons[s].flag[0]=='T') {
        tmp->duplicate(polygons[s]);
        limits.push_back(*tmp);
        }
      else {
        *tmp=plg_subdivide(polygons[s], radius, equalize, PLG_LOXODROMIC, debug);
        limits.push_back(*tmp);
        }
      delete tmp;
      }
    }
  else plg_duplicate(polygons, limits);
    
  ref_lat=(frame.ymax+frame.ymin)/2;
  ref_lon=(frame.xmax+frame.xmin)/2;
  
  printf("%s : reference lat=%lf lon=%lf\n", __func__, ref_lat,ref_lon);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  define stereographic oblique projection, set cartesian grid frame and compute
  polygons cartesian coordinates
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  projection=assign_StereoOblique(ref_lat,ref_lon,pj_parameters);
  geo_to_projection(projection, ref_lat, ref_lon, &x, &y);
    
  printf("\ndefine background working grid: resolution=%lf m, projection=%s (x=%lf y=%lf)\n", resolution, pj_parameters,x,y);
    
  status=plg_cartesian(projection,polygons);
  frame=plg_cartesian_minmax(polygons);
  
  status=plg_cartesian(projection,limits);

  frame=plg_cartesian_minmax(limits);
  
//   status=plg_destroy_entries(limits);
  
//   y=(frame.ymax+frame.ymin)/2;
//   x=(frame.xmax+frame.xmin)/2;
//     
//   projection_to_geo(projection, &ref_lat, &ref_lon, x, y);
  
  frame.dilatation(0.05);

  grid->xmin=frame.xmin;
  grid->xmax=frame.xmax;

  grid->ymin=frame.ymin;
  grid->ymax=frame.ymax;

  grid->dx=resolution;
  grid->dy=resolution;

  grid->xmin=NINT(grid->xmin/grid->dx)*grid->dx-grid->dx;
  grid->xmax=NINT(grid->xmax/grid->dx)*grid->dx+grid->dx;

  grid->ymin=NINT(grid->ymin/grid->dy)*grid->dy-grid->dy;
  grid->ymax=NINT(grid->ymax/grid->dy)*grid->dy+grid->dy;

  grid->nx=(int) NINT((grid->xmax-grid->xmin)/grid->dx)+1;
  grid->ny=(int) NINT((grid->ymax-grid->ymin)/grid->dy)+1;
  
  grid->modeH=0;

  status=map_completegridaxis(grid,2);
  if(status!=0) return(status);

  grid->projection=projection;
  grid->proj4options=poc_strdup(pj_parameters);
  
  deletep(&pj_parameters);
  
  map_printgrid(*grid);
  printf("\n");

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_defgrid(grid_t *grid, plg_t *polygons, int npolygons, vector<plg_t> & limits, double step)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  vector<plg_t> p;
  
  p=plg_array2vector(polygons, npolygons);
  
  status=fe_defgrid(grid, p, limits, step);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int defproj(grid_t *grid, plg_t *polygons, int npolygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  frame_t frame;
  projPJ projection;
  double ref_lat[2], ref_lon[2], x, y;
  char *pj_parameters=new char[512];

//   frame=plg_spherical_minmax(polygons, npolygons);
  frame= plg_recale(polygons, npolygons);

  ref_lat[0]=(frame.ymax+frame.ymin)/2;
  ref_lon[0]=(frame.xmax+frame.xmin)/2;
  
  printf("%s : reference lat=%lf lon=%lf\n", __func__, ref_lat[0],ref_lon[0]);

/*------------------------------------------------------------------------------
  define Oblique Stereographic projection */
  projection=assign_StereoOblique(ref_lat[0],ref_lon[0],pj_parameters);

  status=plg_cartesian(projection,polygons, npolygons);

  frame=plg_cartesian_minmax(polygons, npolygons);

  grid->projection=projection;
  grid->proj4options=poc_strdup(pj_parameters);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_ReloadMeshsize(const char *filename, grid_t *grid, float **density, float *mask, plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  cdfgbl_t global;
  int id;
  cdfvar_t varinfo;
  double dmask;
  frame_t frame;
  int ptype=0;
  projPJ projection;
  double ref_lat[2],ref_lon[2];


  status= cdf_globalinfo(filename,&global,0);
  if(status !=0) goto error;

/*------------------------------------------------------------------------------
  load netcdf grid */
  id=cdf_identify(global,"density");
  status=cdf_loadvargrid_2d (filename,id, grid);
  if(status !=0) goto error;

  id=cdf_identify(global,"x");
  status=cdf_varinfo(filename,"x",&varinfo,1);
  status= poc_getvar2d (filename, id,(int) 0, grid->x, &dmask,varinfo);
  if(status !=0) goto error;

  id=cdf_identify(global,"y");
  status=cdf_varinfo(filename,"y",&varinfo,1);
  status= poc_getvar2d (filename, id,(int)  0, grid->y, &dmask,varinfo);
  if(status !=0) goto error;

  {
  const int nxy=grid->Hsize();
  poc_minmax(grid->x,nxy,dmask,&grid->xmin,&grid->xmax);
  poc_minmax(grid->y,nxy,dmask,&grid->ymin,&grid->ymax);

  grid->dx=(grid->xmax-grid->xmin) / ((double) grid->nx-1.);
  grid->dy=(grid->ymax-grid->ymin) / ((double) grid->ny-1.);

  *density=new float[nxy];
    
  }
  
/*------------------------------------------------------------------------------
  load netcdf variable */
  id=cdf_identify(global,"density");
  status=cdf_varinfo(filename,"density",&varinfo,1);
  status= poc_getvar2d (filename, id, 0, *density, mask,varinfo);
  if(status !=0) goto error;
 
//   grid->modeH=0;
//   
//   delete[] grid->x;grid->x=0;
//   delete[] grid->y;grid->y=0;
//   
//   status=map_completegridaxis(grid,2);
//   if(status!=0) return(status);
  
  map_printgrid(*grid);
  
  printf("\n");

  return(0);

error:
  return(-1);

}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int set_density00(criteria_t criteria, grid_t grid, float *density, float*topo, float mask, float *buffer, float *constrain)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
 maxsize criterion
  
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int    i,j,k;
  float  H,l;
  
  const float
    maxsize=criteria.maxsize*km2m,
    shelf_maxsize=criteria.shelf_maxsize*km2m;
  
  printf("mesh density, step #1, deep sea  maxsize=%lf m\n", criteria.maxsize*km2m);
  printf("mesh density, step #1, shelf sea maxsize=%lf m\n", criteria.shelf_maxsize*km2m);

  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      k=grid.nx*j+i;
      if(density[k]!=mask) {
        l=maxsize;
        H=-topo[k];
        if(constrain!=0) constrain[k]=MSH_MAXSIZE_OPEN;
        if(H<criteria.open_shelf_limit){
          updatemin(&l,shelf_maxsize);
          if(constrain!=0) constrain[k]=MSH_MAXSIZE_SHELF;
          }
        if(buffer!=0) buffer[k]=l;
        density[k]=l;
        }
      }
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int set_density01(criteria_t criteria,grid_t grid, float *density, float*topo, float mask, float *buffer, float *constrain)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
  (tidal) wavelength criterion
  
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int    i,j,k;
  float  H,l;
  
  printf("mesh density, step #2, factor=%lf (dimensionless)\n", criteria.factor1);

  if(buffer) aset(buffer,grid.Hsize(),mask);
  
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      k=grid.nx*j+i;
      if(density[k]!=mask) {
        H=-topo[k];
        if(H<0) {
          H=2.0;
          }
//        density[k]=min(1000.*criteria.factor1*sqrt(g*H),density[k]);
        l=1000.*criteria.factor1*sqrt(H);
/* *----------------------------------------------------------------------------
        warning, temporary change in criteria computation */
//         if(H<250) updatemin(&l,35000.);
//         if(H<450) updatemin(&l,35000.);
        if(buffer!=0) buffer[k]=l;
        density[k]=min(l,density[k]);
        if(density[k]==l and constrain !=0) constrain[k]=MSH_TIDAL_WAVELENGTH;
        }
      }
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int set_density02(criteria_t criteria,grid_t grid, float *density, float*topo, float mask, float *buffer, float *constrain)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
  bottom topography slope criterion
  
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int    i,j,k,status;
  float  H;
  float *dHdx,*dHdy,slope,ratio,l;
  const float factor=criteria.factor2;

  printf("mesh density, step #3, factor=%lf (dimensionless)\n", criteria.factor2);

  if(buffer!=0) for(size_t m=0;m<grid.Hsize();m++) buffer[m]=mask;

  dHdx=new float[grid.nx*grid.ny];
  dHdy=new float[grid.nx*grid.ny];
  status= map_gradient(grid,grid.nx, topo, mask,CARTESIAN, dHdx, dHdy);

  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      k=grid.nx*j+i;
      if(density[k]!=mask) {
        H=-topo[k];
        if(H<0) {
          H=2.0;
          }
        if(dHdx[k]==mask) continue;
        if(dHdy[k]==mask) continue;
/*------------------------------------------------------------------------------
        do not prescribe slope criterion in very shallow water*/
        if(H<criteria.hmin2) continue;
        slope=hypot(dHdx[k],dHdy[k]);
        if(slope!=0) {
          ratio=MAX(1000.*criteria.minratio,H/slope);
          l=factor*ratio;
          if(buffer!=0) buffer[k]=l;
          density[k]=min(l,density[k]);
          if(density[k]==l and constrain!=0) constrain[k]=MSH_TOPO_SLOPE;
          }
        }
      }
    }
  delete[] dHdx;
  delete[] dHdy;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int set_density03(criteria_t criteria, grid_t grid, float *density, float mask,  plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
    boundary criterion
    
    This routine intend to have the mesh density criterion consistent
    with the sampling of the boundaries. It works so that a) within a band
    of 1.5 x length of subsgement (line between 2 consecutive boundary
    points) the criterion is equal the 2xlength of subsgement b) within
    a band of 1.5 x length of subsgement the criterion is equal the
    4xlength of subsgement :

    Outside of Domain
                     L
                 <------->
                 1       2
    -------------X-------X-------------X-------------
           |                   |  |
           |                   |  |
           |  Area modified    |  | 1.5x L
           |                   |  |
           |                   | \|/
           ---------------------  V

    Inside of Domain

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int   s;
  
  double limits_consistency_extent=criteria.limits_consistency_extent;
  double limits_consistency_factor=criteria.limits_consistency_factor;

  float *local=new float[grid.Hsize()];
  
  for (int m=0; m< grid.Hsize(); m++) local[m]=1.e+10;
  
  for(s=0;s<npolygones;s++) {
    if(polygones[s].npt==0) continue;
#pragma omp parallel for
    for(size_t k=0;k<polygones[s].npt-1;k++) {
      plg_t plg;
      plg.init(5);
      if(polygones[s].flag!=0) {
        if(criteria.resample_openlimits  && polygones[s].flag[k]=='M') continue;
        if(criteria.resample_rigidlimits && polygones[s].flag[k]=='T') continue;
        }
      char *interior;
      double d,tx,ty,nx,ny;
      tx=polygones[s].x[k+1]-polygones[s].x[k];
      ty=polygones[s].y[k+1]-polygones[s].y[k];
      nx=-ty;
      ny= tx;
      plg.x[0]=polygones[s].x[k]-0.5*tx-limits_consistency_extent*nx;
      plg.y[0]=polygones[s].y[k]-0.5*ty-limits_consistency_extent*ny;
      plg.x[1]=plg.x[0]+2.*tx;
      plg.y[1]=plg.y[0]+2.*ty;
      plg.x[2]=plg.x[1]+2.*limits_consistency_extent*nx;
      plg.y[2]=plg.y[1]+2.*limits_consistency_extent*ny;
      plg.x[3]=plg.x[2]-2.*tx;
      plg.y[3]=plg.y[2]-2.*ty;
      plg.x[4]=plg.x[0];
      plg.y[4]=plg.y[0];
      
      d=sqrt(tx*tx+ty*ty);
      d*=limits_consistency_factor;
      
      interior=plg_test_grid(grid,&plg,(int) 1);
      for (size_t m=0;m<grid.ny*grid.nx;m++) {
        if(density[m]!=mask) {
          if(interior[m]==1) {
/*------------------------------------------------------------------------------
            warning, temporary change in criteria computation */
#pragma omp critical
              {
              local[m]=min((float) d,local[m]);
//               density[m]=d;
              }
            }
          }
        }
      delete[] interior;
      plg.destroy();
      }
    }
    
  for (int m=0; m< grid.Hsize(); m++) {
    if(local[m]!=1.e+10) {
      if(density[m]!=mask) density[m]=local[m];
      }
    }
  
  delete[] local;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename TYPE> int surface_wavelength(double T, float depth, TYPE & l)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double k,f,eps,g,sig;
  double X,Y,FD,H,error;

  f=1.0/T;

  eps=0.0001;
  g=9.81;
  sig=2.*M_PI*f;

/*------------------------------------------------------------------------------

  sig²=g*k* tanh (kH)
  
  X = kH
  
  Y = Hsig²/g

  Y = X tanh(X)

------------------------------------------------------------------------------*/

  Y=depth*sig*sig/g;    /** a is the squared adimensional frequency*/
  X=sqrt(Y);

  error=1.;
  while (abs(error) > eps) {
    H=tanh(X);
    error=Y-X*H;
    FD=-H-X/cosh(X)/cosh(X);
    X=X-error/FD;
    }

  k=X/depth;
  l=2.*M_PI/k;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int set_density04_a(criteria_t criteria,grid_t grid, float *density, float*topo, float mask, float *buffer, float *constrain)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  Pour le critère, c'est tanh (kH), ou H est la profondeur et k le
  nombre d'onde des vagues défini par :
 
       sig²=g*k* tanh (kH),avec g=9.81
       sig =2*pi/T avec T la période des vagues
 
  k dépend de la profondeur et de la fréquence des vagues (f=1/T, ou T est
  la période des vagues).
  Pour obtenir k, pour notre application on peut fixer T=10s donc f=0.1Hz.
  Il faut donc ensuite calculer k pour chaque profondeur associé aux
  triangles. k doit être calculé en inversion la relation de dispersion
  (sig²=g*k* tanh (kH), comme dans l'exemple matlab fourni ci-dessous.
  Pour gagner du temps il est possible de précalculer tout ça dans un
  table (k vs H).
  
  En fixant la taille du triangle à la cote Dx0 et celle au large Dx1, on
  peut avoir la taille du triangle locale Dx par la formule :
  Dx=Dx0+(Dx1-Dx0)*tanh(kH)
 
  function dispNewtonTH=dispNewtonTH(f,dep)
  % inverts the linear dispersion relation sig^2=g*k*tanh(k*dep) to get
  % k from T=2*pi/sig
  eps=0.0001;
  g=9.81;
     sig=2.*pi.*f;
     Y=dep.*sig.^2./g ;   %a is the squared adimensional frequency
      X=sqrt(Y);
      I=1;
     F=1.;
     while abs(max(F)) > eps
          H=tanh(X);
          F=Y-X.*H;
          FD=-H-X./cosh(X).^2;
          X=X-F./FD;
     end
     dispNewtonTH=X./dep;
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  

  surface wavelength deformation criterion
  
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int    i,j,n,status,count=0;
  float  H;
  float *dHdx,*dHdy,slope;
  double T;
  float k,kk,lamda,l,d;

  printf("mesh density, step #4, surface wave wavelength wave_period=%lf  factor=%lf  (dimensionless)\n", criteria.surface_wave_period, criteria.factor4);

  if(buffer!=0) for(size_t m=0;m<grid.Hsize();m++) buffer[m]=mask;
  
  if(criteria.factor4<=0.0) {
    printf("surface wave factor not consistent (actually set to zero), abort\n");
     TRAP_ERR_EXIT(-1," you might check the input criteria file to solve the issue\n");
    }

  dHdx=new float[grid.nx*grid.ny];
  dHdy=new float[grid.nx*grid.ny];
  status= map_gradient(grid,grid.nx, topo, mask, CARTESIAN, dHdx, dHdy);
  
/*------------------------------------------------------------------------------
  Surface wave period (sec) */
  T=criteria.surface_wave_period;
  
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      n=grid.nx*j+i;
      if(density[n]!=mask) {
        H=-topo[n];
        if(H<0) {
          H=2.0;
          }
        if(dHdx[n]==mask) continue;
        if(dHdy[n]==mask) continue;
        slope=sqrt(dHdx[n]*dHdx[n]+dHdy[n]*dHdy[n]);
        if(slope!=0) {
          status=surface_wavelength(T,H,lamda);
/*------------------------------------------------------------------------------
          Surface wave wavenumber */
          k=2*M_PI/lamda;
          kk=slope*k/(0.5*sinh(2*k*H)+k*H);
/*------------------------------------------------------------------------------
          deformation length */
          l=2*M_PI/kk;
/*------------------------------------------------------------------------------
          adequate resolution */
          l/=30.0;
/*------------------------------------------------------------------------------
          resolution constraint */
          d=criteria.factor4*l;
          if(buffer!=0) buffer[n]=d;
          density[n]=min(d,density[n]);
          if(density[n]==d and constrain!=0) constrain[n]=MSH_SURFWAVE_LENGTH;
          if(H<20) {
            count++;
            }
          }
        }
      }
    }
  
  delete[] dHdx;
  delete[] dHdy;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int set_density04_b(criteria_t criteria,grid_t grid, float *density, float*topo, float mask, float *buffer, float *constrain)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  surface wavelength criterion, dZmax type (Aaron)
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int    j,status;
  float *dHdx,*dHdy;
  float T;
  float dZmax=4.;

  printf("mesh density, step #4, surface wave dZmax wave_period=%lf  factor=%lf  (dimensionless)\n", criteria.surface_wave_period, criteria.factor4);

  if(buffer!=0) for(size_t m=0;m<grid.Hsize();m++) buffer[m]=mask;

  if(criteria.factor4<=0.0) {
    printf("surface wave factor not consistent (actually set to zero), abort\n");
     TRAP_ERR_EXIT(-1," you might check the input criteria file to solve the issue\n");
    }

  dHdx=new float[grid.nx*grid.ny];
  dHdy=new float[grid.nx*grid.ny];
  status= map_gradient(grid,grid.nx, topo, mask, CARTESIAN, dHdx, dHdy);

/*------------------------------------------------------------------------------
  Surface wave period (sec) */
  T=criteria.surface_wave_period;
  
#pragma omp parallel for
  for(j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      size_t n=grid.nx*j+i;
      if(density[n]!=mask) {
        float H=-topo[n];
        if(H<0) {
          H=2.0;
          }
        if(dHdx[n]==mask) continue;
        if(dHdy[n]==mask) continue;
        float slope=sqrt(dHdx[n]*dHdx[n]+dHdy[n]*dHdy[n]);
        if(slope!=0) {
          float lamda;
          status=surface_wavelength(T,H,lamda);
/*------------------------------------------------------------------------------
          Surface wave wavenumber */
          float k=2*M_PI/lamda;
          float l = dZmax / (slope * ( (4.*k*H) / ( 2.*k*H + sinh(2.*k*H) ) ) );
          if(buffer!=0) buffer[n]=l;
          density[n]=min(l,density[n]);
          if(density[n]==l and constrain !=0) constrain[n]=MSH_SURFWAVE_DZ;
          }
        }
      }
    }

  delete[] dHdx;
  delete[] dHdy;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int set_density04_c(criteria_t criteria,grid_t grid, float *density, float*topo, float mask, float *buffer, float *constrain)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  surface wavelength CFL criterion (Aaron)
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int   i,j,n, status;
  float H,g=9.81;
  float Cg,G,c,dt;
  double T;
  float k,lamda,l,d;

  printf("mesh density, step #4, surface wave CFL desired_max_timestep=%lf  \n", criteria.surface_wave_dt);
/*------------------------------------------------------------------------------
  Surface wave period (sec) */
  T=criteria.surface_wave_period;
  dt=criteria.surface_wave_dt;

  if(buffer!=0) for(size_t m=0;m<grid.Hsize();m++) buffer[m]=mask;
  
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      n=grid.nx*j+i;
      if(density[n]!=mask) {
        H=-topo[n];
        if(H<0) {
          H=2.0;
          }
        status=surface_wavelength(T,H,lamda);
/*------------------------------------------------------------------------------
        Surface wave wavenumber */
        k=2*M_PI/lamda;
        G=2.*H*k/sinh(2*k*H);
        c=sqrt(g*tanh(k*H)/k);
        Cg=0.5*c*(1+G);
        l=Cg*dt;
        d=l;
        density[n]=max(d,density[n]);
        if(buffer!=0) buffer[n]=l;
        if(density[n]==d and constrain !=0) constrain[n]=MSH_SURFWAVE_CFL;
        }
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int set_density05(criteria_t criteria,grid_t grid, float *density, float*topo, float mask, float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------
  (tidal) wavelength criterion from prior solutions*/
{
  int    i,j,k,status;
  float  l,*L=0,Lmask,z;
  grid_t Lgrid;
  double x,y,lon,lat;

  if(buffer!=0) for(size_t m=0;m<grid.Hsize();m++) buffer[m]=mask;
  
  status=map_loadfield_cdf("/home/softs/data/tides/FES2004/netcdf/lamda.nc", "lamda4_15", Lgrid, L, Lmask);
  
  if(status!=0) return(-1);

  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      k=grid.nx*j+i;
      if(density[k]!=mask) {
        x=grid.x[k];
        y=grid.y[k];
        projection_to_geo(grid.projection,&lat,&lon,x,y);
        lon=map_recale(Lgrid,lon);
        status=map_interpolation(Lgrid, L,Lmask,lon,lat,&z);
        if (z!=Lmask) {
          l=z*km2m;
          if(buffer!=0) buffer[k]=l;
          density[k]=min(l,density[k]);
          }
        }
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int set_density06(criteria_t criteria,grid_t grid, float *density, float *topo, float mask, float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------
  baroclinic wavelength criterion */
{
  int    i,j,k,status;
  grid_t Lgrid;
  double x,y;
  double N=2.e-03, cmin;
  int nlevels=40;
  double omega=2*M_PI/12./3600;
  
  float *dHdx=new float[grid.nx*grid.ny];
  float *dHdy=new float[grid.nx*grid.ny];
  status= map_gradient(grid, grid.nx, topo, mask,CARTESIAN, dHdx, dHdy);

  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      k=grid.nx*j+i;
      if(density[k]!=mask) {
        x=grid.x[k];
        y=grid.y[k];
        double H=-topo[k];
        cmin=N*H/(nlevels-1)/M_PI;
        double K=omega/cmin;
        double L=2.*M_PI/K;
        float slope=sqrt(dHdx[k]*dHdx[k]+dHdy[k]*dHdy[k]);
        if(slope>0.01)
          density[k]=min((float) L,density[k]);
        }
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int smooth_density_clx(const grid_t & grid, float *density, float mask, float slopemax)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k,m,n,status,count;
  int    ii,jj;
  float *dHdx,*dHdy,slope;
  int    proxies[9][2]={1,0,1,1,0,1,-1,1,-1,0,-1,-1,0,-1,1,-1,0,0};
  double dx[9],x1,y1,x2,y2;

  dHdx=new float[grid.nx*grid.ny];
  dHdy=new float[grid.nx*grid.ny];

  status= map_gradient(grid,grid.nx, density, mask,CARTESIAN, dHdx, dHdy);

  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      k=grid.nx*j+i;
      if(density[k]!=mask) {
        slope=sqrt(dHdx[k]*dHdx[k]+dHdy[k]*dHdy[k]);
//        density[k]=slope;
        }
      }
    }

  i=1;
  j=1;
  m=grid.nx*j+i;
  x1=grid.xmin+grid.dx;
  y1=grid.ymin+grid.dy;
  for(k=0;k<8;k++) {
    m=grid.nx*j+i;
    ii=i+proxies[k][0];
    jj=j+proxies[k][1];
    n=grid.nx*jj+ii;
    x2=grid.xmin+ii*grid.dx;
    y2=grid.ymin+jj*grid.dy;
    dx[k]=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    }

/*------------------------------------------------------------------------------
  count is the number of intervention at each loop */
//  count=0;
//   while (count<niterations) {
//   count++;
//  count=0;
  count=1000;
  while (count>10) {
    count=0;
    for(j=1;j<grid.ny-1;j++) {
      for(i=1;i<grid.nx-1;i++) {
        m=grid.nx*j+i;
        if(density[m]!=mask) {
          for(k=0;k<8;k++) {
            ii=i+proxies[k][0];
            jj=j+proxies[k][1];
            n=grid.nx*jj+ii;
            if(density[n]!=mask) {
              slope=fabs(density[m]-density[n])/dx[k];
              if(isnan(slope)) {
                slope=0;
                }
              if(1+fabs(slope)>slopemax) { /// HERE !!!
                slope=0.99*slopemax-1; /// HERE !!!
                count++;
                if(density[m]>density[n]) {
                  density[m]=density[n]+slope*dx[k];
                  }
                else {
                  density[n]=density[m]+slope*dx[k];
                  }
                }
              }
            }
          }
        }
      }
    }


  delete[] dHdx;
  delete[] dHdy;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int smooth_density(const grid_t & grid, float *density, float mask, float slopemax)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k,m,n,/*status,*/count;
  int    ii,jj;
  float /**dHdx,*dHdy,*/slope;
  int    proxies[9][2]={1,0,1,1,0,1,-1,1,-1,0,-1,-1,0,-1,1,-1,0,0};
  double dx[9],x1,y1,x2,y2,factor;

//   dHdx=new float[grid.nx*grid.ny];
//   dHdy=new float[grid.nx*grid.ny];
//  status= map_gradient(grid,grid.nx, density, mask,CARTESIAN, dHdx, dHdy);


/*------------------------------------------------------------------------------
  grid is given in cartesian coordinates, resolution is fixed */
  i=1;
  j=1;
  m=grid.nx*j+i;
  x1=grid.xmin+grid.dx;
  y1=grid.ymin+grid.dy;
  for(k=0;k<8;k++) {
    m=grid.nx*j+i;
    ii=i+proxies[k][0];
    jj=j+proxies[k][1];
    n=grid.nx*jj+ii;
    x2=grid.xmin+ii*grid.dx;
    y2=grid.ymin+jj*grid.dy;
    dx[k]=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    }

    
/*------------------------------------------------------------------------------
  count is the number of intervention at each loop */
  count=1000;
  while (count>10) {
    count=0;
    for(j=1;j<grid.ny-1;j++) {
      for(i=1;i<grid.nx-1;i++) {
        m=grid.nx*j+i;
        if(density[m]!=mask) {
          for(k=0;k<8;k++) {
            ii=i+proxies[k][0];
            jj=j+proxies[k][1];
            n=grid.nx*jj+ii;
            if(density[n]!=mask) {
              slope=fabs(density[m]-density[n])/dx[k];
              factor=density[m]/dx[k];
              if(isnan(slope)) {
                slope=0;
                }
              if(factor*fabs(slope)>slopemax) {
                slope=0.99*slopemax/factor;
                count++;
                if(density[m] > density[n]) {
                  updatemin(&density[m],density[n]+slope*dx[k]);
                  }
                else {
                  updatemin(&density[n],density[m]+slope*dx[k]);
                  }
                }
              }
            }
          }
        }
      }
    }

//    #pragma omp critical(threadUnsafeNetcdf)

//   delete[] dHdx;
//   delete[] dHdy;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int smooth_density_new(const grid_t & grid, float *density, float mask, float slopemax)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count;
  int    proxies[8][2]={1,0,1,1,0,1,-1,1,-1,0,-1,-1,0,-1,1,-1};
  double d[8];

  float *buffer=0;
  
  if(slopemax>1) {
    printf("\n\WARNING !!! reworked smoothing algorithm, request 0.15x(ancient slope limiter=%lf?) setting\n\n",slopemax);
    }
  poc_copy(buffer,density,grid.Hsize());

/*------------------------------------------------------------------------------
  grid is given in cartesian coordinates, resolution is fixed */
  for(int k=0;k<8;k++) {
    const int
    ii=proxies[k][0],
    jj=proxies[k][1];
    const double
    dx=ii*grid.dx,
    dy=jj*grid.dy;
    d[k]=hypot(dx,dy);
    }
  
#define DEBUG_smooth_density_new 0
#if DEBUG_smooth_density_new
  const string path="smooth_density_new.nc";
  poc_var_t var;
  var.init("smoothing",NC_FLOAT);
  var
    <<poc_att_t(_FillValue,mask)
    <<poc_dim_t("time",0)
    <<poc_dim_t("ny",grid.ny)
    <<poc_dim_t("nx",grid.nx);
  poc_create(path.c_str(),poc_global_t()<<var);
  int frame=0;
#endif
  
/*------------------------------------------------------------------------------
  count is the number of intervention at each loop */
  do {
#pragma omp parallel for
    for(int j=1;j<grid.ny-1;j++) {
      size_t m,n;
      int i,ii,jj;
      float bm1;
      
      for(i=1;i<grid.nx-1;i++) {
        m=grid.nx*j+i;        
        bm1=density[m];
        if(bm1==mask)
          continue;
        
        if(isnan(bm1))
          TRAP_ERR_EXIT(ENOEXEC,"NAN found\n");
        
        for(int k=0;k<8;k++) {
          ii=i+proxies[k][0];
          jj=j+proxies[k][1];
          n=grid.nx*jj+ii;
          if(density[n]==mask)
            continue;
          updatemin(&bm1,density[n]+slopemax*d[k]);
          }
        buffer[m]=bm1;
        }
      }
    
    count=0;
    for(size_t m=0;m<grid.Hsize();m++) {
      if(buffer[m]==density[m])
        continue;
      count++;
      density[m]=buffer[m];
      }
    
#if DEBUG_smooth_density_new
    poc_put_vara(path,var,frame++,buffer);
#endif
    
    } while (count>10);

  delete[] buffer;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_selectedges(mesh_t mesh, grid_t grid, float *density, float mask, int *selected)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,m,n,n1,n2,status;
  int   count=0;
  int   nelts,nndes,nedges;
  edge_t *edges;
  triangle_t *elt;
  double t,p,d[3],t1,t2,p1,p2;
  double dmax;
  float z[3];

  elt=mesh.triangles;

  nelts =mesh.ntriangles;
  nndes =mesh.nvtxs;
  nedges=mesh.nedges;
  edges =mesh.edges;

  printf ("check size criterion...\n");

/*-----------------------------------------------------------------------------
  distance criterion */
  for (m=0;m<mesh.ntriangles;m++) {
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].edges[k];
      n1=edges[n].extremity[0];
      n2=edges[n].extremity[1];
      t1=mesh.vertices[n1].lon;
      t2=mesh.vertices[n2].lon;
      p1=mesh.vertices[n1].lat;
      p2=mesh.vertices[n2].lat;
      d[k]=sqrt((t2-t1)*(t2-t1)+(p2-p1)*(p2-p1));
      }
    status=map_interpolation(grid, density,mask,t,p,&z[k]);
    if(d[k]>dmax) {
      selected[n]=1;
      count++;
      }
    else {
      selected[n]=0;
      }
    }
  printf ("number of edges selected %d (/%d)\n",count,nedges);

  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_DecimateInterior(criteria_t & criteria, mesh_t & mesh, grid_t & grid, float *density, float mask, plg_t *polygons, int npolygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  int n1,n2;
  char *keep,*passed;
  float *Lsize;
  double x,y;
  
  double Lmin=min(criteria.minsize,criteria.shelf_minsize)*km2m;
  double Lmax=max(criteria.minsize,criteria.shelf_minsize)*km2m;

  if(Lmin<=0.0) {
    printf("prescribed element minimum size=%lf, abort\n", Lmin);
     TRAP_ERR_EXIT(-1," you might check the input criteria file to solve the issue\n");
    }

  passed =new char [mesh.nvtxs];
  keep   =new char [mesh.nvtxs];
  Lsize  =new float[mesh.nvtxs];

  for(n=0;n<mesh.nvtxs;n++) {
    keep[n]=0;
    passed[n]=0;
    mesh.vertices[n].nedges=0;
    }

  for(n=0;n<mesh.nedges;n++) {
    n1=mesh.edges[n].extremity[0];
    n2=mesh.edges[n].extremity[1];
    mesh.vertices[n1].nedges++;
    mesh.vertices[n2].nedges++;
    }

  for(n=0;n<mesh.nvtxs;n++) {
    if(mesh.vertices[n].edges!=0) delete[] mesh.vertices[n].edges;
    mesh.vertices[n].edges=new int[mesh.vertices[n].nedges];
    mesh.vertices[n].nedges=0;
    x=mesh.vertices[n].lon;
    y=mesh.vertices[n].lat;
    status=map_interpolation(grid, density,mask,x,y,&(Lsize[n]));
    if(Lsize[n]==mask) Lsize[n]=criteria.maxsize*km2m;
    }

  for(n=0;n<mesh.nedges;n++) {
    n1=mesh.edges[n].extremity[0];
    n2=mesh.edges[n].extremity[1];
    mesh.vertices[n1].edges[mesh.vertices[n1].nedges]=n;
    mesh.vertices[n2].edges[mesh.vertices[n2].nedges]=n;
    mesh.vertices[n1].nedges++;
    mesh.vertices[n2].nedges++;
    }

  for(n=0;n<mesh.nvtxs;n++) {
    if(passed[n]==1) continue;
/*------------------------------------------------------------------------------
    keep all boundary nodes */
    if(mesh.vertices[n].code!=0) {
      keep[n]=1;
      passed[n]=1;
      continue;
      }
/*------------------------------------------------------------------------------
    ??? */
//     if(Lsize[n]<Lmin) {
//       keep[n]=1;
//       passed[n]=1;
//       continue;
//       }
    for(k=0;k<mesh.vertices[n].nedges;k++) {
      m=mesh.vertices[n].edges[k];
      n1=mesh.edges[m].extremity[0];
      n2=mesh.edges[m].extremity[1];
      
      plg_t p=plg_t(mesh.vertices[n1].lon,mesh.vertices[n1].lat,mesh.vertices[n2].lon,mesh.vertices[n2].lat);
      double radius=criteria.cellsize*km2m;
      int equalize=1;
      plg_t q=plg_resample(p, radius, equalize);
      float L=1.e+10;
      for(size_t l=0;l<q.npt;l++) {
        float LL;
        status=map_interpolation(grid,density,mask,q.x[l],q.y[l],&LL);
//         LL=5000.0;
        if(LL==mask) {
          LL=criteria.maxsize*km2m;
          }
        updatemin(&L,LL);
        }
//      L=min(Lsize[n1],Lsize[n2]);
      q.destroy();
      p.destroy();
      if(mesh.edges[m].L>0.75*L) {
        keep[n]=1;
        break;
        }
      }
//    if(keep[n]==1) {
/*------------------------------------------------------------------------------
      do not touch neighbours */
      for(k=0;k<mesh.vertices[n].nngh;k++) {
        m=mesh.vertices[n].ngh[k];
        keep[m]=1;
        passed[m]=1;
        }
//      continue;
//      }
    passed[n]=1;
    }

  for(n=0;n<mesh.nvtxs;n++) {
    if(keep[n]==0) {
      mesh.vertices[n].code=-1;
      }
    }
    
#ifdef EXTERN_TRIANGLE
  status=fe_savenodes("tmp.nod", NODE_FILE_FORMAT_TRIGRID, mesh);
#endif
  
  delete[] passed;
  delete[] keep;
  delete[] Lsize;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t fe_nodit(criteria_t criteria, grid_t sgrid, float *density, float mask, plg_t *polygones, int npolygones, bool debug, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,n,step,status;
  int ii,jj,kk;
  int count=0;
  double x,y;
  char *interior;
  grid_t *grid;
  point_t *points;
  mesh_t *mesh,*spherical,*previous;
  double Lmin=min(criteria.minsize,criteria.shelf_minsize)*km2m;
  int previous_nvtxs=0;
  
  grid=&sgrid;
  
/*------------------------------------------------------------------------------
  flag interior cells*/
  grid->projection=sgrid.projection;
  interior=plg_test_grid(*grid,polygones,npolygones);

  count=0;
  for(j=0;j<grid->ny;j++) {
    for(i=0;i<grid->nx;i++) {
      k=grid->nx*j+i;
      if(interior[k]==1) count++;
      }
    }

  printf("initial number of inner nodes: %d\n", count);
    
/*------------------------------------------------------------------------------
  eliminate some unnecessary points*/
  count=0;
  for(j=2;j<grid->ny-2;j+=2) {
    for(i=2;i<grid->nx-2;i+=2) {
      k=grid->nx*j+i;
      if(density[k]==mask) continue;
      if(interior[k]!=1)   continue;
      if(density[k]>3.*grid->dx) {
        for(jj=j-2;jj<j+3;jj++) {
          for(ii=i-2;ii<i+3;ii++) {
            kk=grid->nx*jj+ii;
//            if(interior[kk]!=1)   goto next;
            if(density[kk]==mask) goto next;
            if(density[kk]<3.*grid->dx) {
              goto next;
              }
            }
          }
        for(jj=j-1;jj<j+2;jj++) {
          for(ii=i-1;ii<i+2;ii++) {
            kk=grid->nx*jj+ii;
            if(kk!=k) interior[kk]=-1;
            }
          }
next:
        count++;
        }
      }
    }

  count=0;
  for(j=0;j<grid->ny;j++) {
    for(i=0;i<grid->nx;i++) {
      k=grid->nx*j+i;
      if(interior[k]==1) count++;
      }
    }
    
  printf("first filter applied, number of inner nodes: %d\n", count);

  points=new point_t[count];

  count=0;
  for(j=0;j<grid->ny;j++) {
    for(i=0;i<grid->nx;i++) {
      k=grid->nx*j+i;
      if(interior[k]==1) {
        points[count].t=grid->x[k];
        points[count].p=grid->y[k];
        count++;
        }
      }
    }

  delete[] interior;
   
/*------------------------------------------------------------------------------
  create node set by merging boundary nodes and interior nodes sets */
  mesh=fe_createnodes(polygones, npolygones, points, count);
  delete[] points;
  
  if(debug) {
    for (n=0;n<mesh->nvtxs;n++) {
      x=mesh->vertices[n].lon;
      y=mesh->vertices[n].lat;
      projection_to_geo(sgrid.projection,&(mesh->vertices[n].lat),&(mesh->vertices[n].lon),x,y);
      }
    status=fe_savenodes("tmp.nod", NODE_FILE_FORMAT_TRIGRID, *mesh);
    for (n=0;n<mesh->nvtxs;n++) {
      x=mesh->vertices[n].lon;
      y=mesh->vertices[n].lat;
      geo_to_projection(sgrid.projection,y,x,&(mesh->vertices[n].lon),&(mesh->vertices[n].lat));
      }
    }

#ifdef EXTERN_TRIANGLE
  status=fe_savenodes("tmp.nod", NODE_FILE_FORMAT_TRIGRID, *mesh);
  status=system("node2poly.exe tmp.nod tmp.poly");
  if(status!=0) return(NULL);
  status=system("triangle.exe -p tmp.poly");
  if(status!=0) return(NULL);
#endif

  spherical=fe_triangulate(*mesh, sgrid.projection, 0, debug);
  mesh->destroy();

//   for (n=0;n<spherical->nvtxs;n++) {
//     x=spherical->vertices[n].lon;
//     y=spherical->vertices[n].lat;
//     projection_to_geo(sgrid.projection,&(spherical->vertices[n].lat),&(spherical->vertices[n].lon),x,y);
//     }
//   status=fe_savemesh("mesh-critere-raw.nei",MESH_FILE_FORMAT_TRIGRID, *spherical);


  if(verbose) printf(
    "#################################################################\n"
    "decimate node set according to mesh resolution\n\n");
  
  for (step=0;step<criteria.niterations;step++) {
/*------------------------------------------------------------------------------
    eliminate nodes*/
#ifdef EXTERN_TRIANGLE
    spherical=new mesh_t;
    status=fe_init(spherical);
    status=fe_readmesh_TGL("tmp.1.ele","tmp.1.node",spherical);
#endif
    if(verbose) printf("step %3d: nvtxs=%6d\n", step,spherical->nvtxs);

    spherical->type=1;

    status=fe_geometry(spherical);

    status=fe_edgetable(spherical,0,0);
    if(status!=0) {
      for (n=0;n<spherical->nvtxs;n++) {
        x=spherical->vertices[n].lon;
        y=spherical->vertices[n].lat;
        projection_to_geo(sgrid.projection,&(spherical->vertices[n].lat),&(spherical->vertices[n].lon),x,y);
        }
      status=fe_savemesh("mesh-abort.nei",MESH_FILE_FORMAT_TRIGRID, *spherical);
      STDOUT_BASE_LINE("abandon mesh generation, current mesh in %s\n", "mesh-abort.nei");
      exit(-1);
      }
    status=fe_vertex_crosstables02(spherical);
    status=fe_codetable2(spherical,0,1,1);
    if(status!=0) {
      for (n=0;n<spherical->nvtxs;n++) {
        x=spherical->vertices[n].lon;
        y=spherical->vertices[n].lat;
        projection_to_geo(sgrid.projection,&(spherical->vertices[n].lat),&(spherical->vertices[n].lon),x,y);
        }
      status=fe_savemesh("mesh-abort.nei",MESH_FILE_FORMAT_TRIGRID, *spherical);
      STDOUT_BASE_LINE("abandon mesh generation, current mesh in %s\n", "mesh-abort.nei");
      exit(-1);
      }

    status=fe_reshapeall(*spherical,3);

    if(debug) {
      for (n=0;n<spherical->nvtxs;n++) {
        x=spherical->vertices[n].lon;
        y=spherical->vertices[n].lat;
        projection_to_geo(sgrid.projection,&(spherical->vertices[n].lat),&(spherical->vertices[n].lon),x,y);
        }
      char filename[1024];
      sprintf(filename,"tmp.%2.2d.nei",step);
      status=fe_savemesh(filename,MESH_FILE_FORMAT_TRIGRID, *spherical);
      for (n=0;n<spherical->nvtxs;n++) {
        x=spherical->vertices[n].lon;
        y=spherical->vertices[n].lat;
        geo_to_projection(sgrid.projection,y,x,&(spherical->vertices[n].lon),&(spherical->vertices[n].lat));
        }
      }

/*------------------------------------------------------------------------------
    to be optimized */
    status=fe_DecimateInterior(criteria, *spherical, sgrid, density, mask, polygones, npolygones);
    
    status=fe_savenodes("tmp.nod", NODE_FILE_FORMAT_TRIGRID, *spherical);
    spherical->destroy();
    
    status=fe_readnodes("tmp.nod", NODE_FILE_FORMAT_TRIGRID, spherical);

#ifdef EXTERN_TRIANGLE
    status=system("node2poly.exe tmp.nod tmp.poly");
    if(status!=0) return(NULL);
    status=system("triangle.exe -p tmp.poly");
    if(status!=0) return(NULL);
#endif

    status=system("rm -f tmp.nod");
      
    previous=spherical;
    spherical=fe_triangulate(*previous, 0, 0, debug);
    previous->destroy();
    delete previous;
    
    if(spherical->nvtxs==previous_nvtxs) {
      break;
      }
    previous_nvtxs=spherical->nvtxs;
    }

#ifdef EXTERN_TRIANGLE
  spherical=new mesh_t;
  status=fe_init(spherical);
  status=fe_readmesh_TGL ("tmp.1.ele","tmp.1.node",spherical);
#endif
  
  status=fe_edgetable(spherical,0,0);
  status=fe_vertex_crosstables02(spherical);
  status=fe_codetable2(spherical,0,1,1);
  
/*------------------------------------------------------------------------------
  finalize mesh*/
  for (n=0;n<spherical->nvtxs;n++) {
    x=spherical->vertices[n].lon;
    y=spherical->vertices[n].lat;
    projection_to_geo(sgrid.projection,&(spherical->vertices[n].lat),&(spherical->vertices[n].lon),x,y);
    }
  
  spherical->type=0;
  
  return(*spherical);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int initialise_meshgrid(const char *output,const criteria_t & criteria, grid_t *grid, float **density, float *mask, plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  create and initialise a cartesian grid at the required resolution

 *output UNUSED
 criteria
 *grid
 **density allocated to what will be \c grid->Hsize() and set to 1, 0.5 or what will be \c *mask
 *mask
 *polygones used with fe_defgrid()
 npolygones

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int count,status;
  int i,j,k;
  int ii,jj,kk;
  char *interior;
  vector<plg_t> limits;
  
  if(polygones==0) {
    printf("mesh limits not set\n");
    return(-1);
    }
    
  if(criteria.cellsize<=0) {
    printf("density grid resolution not set\n");
    return(-1);
    }
    
  printf("#################################################################\n");
  printf("detect interior cells\n");
  
/*------------------------------------------------------------------------------
  criteria size are given in km, need to be converted in m*/
  status=fe_defgrid(grid, polygones, npolygones, limits, criteria.cellsize*km2m);
  if(status!=0) return(-1);
  
  grid->nz=1;
  grid->z=NULL;

  interior=plg_test_grid(*grid,limits);

  *mask=-99999.0;
  
/*------------------------------------------------------------------------------
  initialise density with 1 for interior (polygon-wise) grid points, *mask otherwise */
  *density=new float[grid->Hsize()];
  for(j=0;j<grid->ny;j++) {
    for(i=0;i<grid->nx;i++) {
      k=grid->nx*j+i;
      if(interior[k]==1) (*density)[k]=interior[k];
      else (*density)[k]=*mask;
      }
    }

/*------------------------------------------------------------------------------
  initialise density with 0.5 for nearly interior (polygon-wise) grid points*/
  count=0;
  for(j=0;j<grid->ny;j++) {
    for(i=0;i<grid->nx;i++) {
      k=grid->nx*j+i;
      if((*density)[k]!=1) {
        for(jj=max(0,j-1);jj<min(j+2,grid->ny);jj++) {
          for(ii=max(0,i-1);ii<min(i+2,grid->nx);ii++) {
            kk=grid->nx*jj+ii;
            if((*density)[kk]==1) {
              (*density)[k]=0.5;
              goto next;
              }
            }
          }
next:   count++;
        }
      }
    }
  
  delete[] interior;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void rescale_density(const grid_t & grid, float *density, float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// rescale density to correct for the projection effects away from the center
/*----------------------------------------------------------------------------*/
{
  const int n=grid.Hsize();
  const double a=2*MeanEarthRadius;/* adjacent */
  
#pragma omp parallel
  {
  int m;
  double x,y;
  double b,h,c;
  float *densitym;
  
#pragma omp for
  for(m=0;m<n;m++){
    densitym=&density[m];
    
    if(*densitym==mask)
      continue;
    
    grid.xy(m,x,y);
    b=hypot(x,y);
    h=hypot(a,b);
    c=h/a;
    *densitym *=square(c);
    }
  
  }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_ComputeMeshsize(const char *output, const criteria_t & criteria, grid_t *grid, float **density, float *mask, plg_t *polygones, int npolygones, int initialize, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  grid_t topogrid;
  float *topo=0,*topobase=0,topomask;
  float *buffer=0, *constrain=0;
  poc_var_t var;
  pocgrd_t ncgrid;
  bool failed=false;
  bool write=(output!=0 and strlen(output)!=0);
  
  if(initialize==1) {
/*------------------------------------------------------------------------------
    initialise the working cartesian grid at the required resolution */
    status=initialise_meshgrid(output, criteria, grid, density, mask, polygones, npolygones);
    if(status!=0) {
      return(-1);
      }
    }
  
  if(write) {
    status=unlink(output);
    status=poc_sphericalgrid_xy(output,"",*grid,&ncgrid);
    status=poc_def_att(output,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
    status=poc_inq_var(output,ncgrid.lat->name,&var);
    var.attributes.clear();
    var<<poc_att_t("coordinates",(string)ncgrid.lat->name+" "+ncgrid.lon->name);
    }

/*------------------------------------------------------------------------------
  load topography */
  printf("#################################################################\n");
  printf("load topography: %s\n", criteria.regulardepth);
  
  status=topo_loadfield(criteria.regulardepth, &topogrid, &topobase, &topomask, debug);
  if(status!=0){
    printf("WARNING not a workable topography !!! Setting to default %g m.\n",criteria.hmax);
    
    topogrid=get_zonegrid("30.");
    topogrid.print(cout,"topo:  ");
    
    topobase=aset(topogrid.Hsize(),(float)criteria.hmax);
    topomask=-1.;
    }


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  interpolate depth on node density grid
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  topo= new float[grid->nx*grid->ny];
  
  double *lon=new double[grid->Hsize()];
  double *lat=new double[grid->Hsize()];
  status=projection_to_geo(grid->proj4options, grid->x, grid->y, lon, lat, grid->Hsize());
  
#pragma omp parallel for private(status)
  for(int j=0;j<grid->ny;j++) {
    for(int i=0;i<grid->nx;i++) {
      int k=j*grid->nx+i;
      float z;
      status=map_interpolation(topogrid, topobase,topomask,lon[k],lat[k],&z);
      if (z!=topomask)
        topo[k]=z;
      else {
        topo[k]=*mask;
        }
      }
    }
  delete[] lon;
  delete[] lat;
  
  for(int j=0;j<grid->ny;j++) {
    for(int i=0;i<grid->nx;i++) {
      int k=j*grid->nx+i;
      if (topo[k]==*mask){
        if((*density)[k]!=*mask) {
          failed=true;
          }
        }
      }
    }
  
  topogrid.free();
  if(topobase!=0) delete[] topobase;
  
  if(failed) {
    printf("warning: %s does not cover all the interior cells in grid\n", criteria.regulardepth);
    }

  if(write) status=poc_put_var(output,var.init("topography",NC_FLOAT,"topography","m",topomask),topo);
   var.attributes.erase(_FillValue);


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  apply criteria and store for checks
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("compute mesh density criteria\n");
  if(write) printf("save density setting progression to %s\n",output);

  printf("\nopen waters / shelf limit set to %lf (m) \n",criteria.open_shelf_limit);
  
  buffer=new float [grid->Hsize()];
  constrain=aset(grid->Hsize(),*mask);

  printf("\nmesh density, step #1 (maxsize limitation)\n");
  status= set_density00(criteria, *grid, *density, topo, *mask, buffer, constrain);
  
  if(write) status=poc_put_var(output,var.init("density-0",NC_FLOAT,"","m",*mask,"maximum cell size"),*density);

/*------------------------------------------------------------------------------
  barotropic gravity wave, elevation wavelength constraint (Kelvin wave)*/
  if(criteria.use(0)) {
    printf("\nmesh density, step #2 (gravity wavelength)\n");
    status= set_density01(criteria, *grid, *density, topo, *mask, buffer, constrain);
    }
  if(write) status=poc_put_var(output,var.init("density-1",NC_FLOAT,"","m",*mask,"gravity wavelength"),*density);

/*------------------------------------------------------------------------------
  barotropic gravity wave, velocity constraint (Kelvin wave) */
  if(criteria.use(1)) {
    printf("\nmesh density, step #3 (depth gradient)\n");
    status= set_density02(criteria, *grid, *density, topo, *mask, buffer, constrain);
    }
  if(write) status=poc_put_var(output,var.init("density-2",NC_FLOAT,"","m",*mask,"depth gradient"),*density);

/*------------------------------------------------------------------------------
  ocean surface wave, LEGOS typical deformation length constraint */
  if(criteria.use(3)) {
    printf("\nmesh density, step #4 (ocean waves wavelength deformation)\n");
    status= set_density04_a(criteria, *grid, *density, topo, *mask, buffer, constrain);
    if(buffer!=0) {
      status=poc_put_var(output,var.init("density-4a-alone",NC_FLOAT,"","m",*mask,"ocean waves wavelength"),buffer);
      }
    }
  if(write) status=poc_put_var(output,var.init("density-4a",NC_FLOAT,"","m",*mask,"ocean waves wavelength"),*density);
  
/*------------------------------------------------------------------------------
  ocean surface wave, typical deformation length constraint */
  if(criteria.use(6)) {
    printf("\nmesh density, step #4 (ocean waves dZ max)\n");
    status= set_density04_b(criteria, *grid, *density, topo, *mask, buffer, constrain);
    if(buffer!=0) {
      if(write) status=poc_put_var(output,var.init("density-4b-alone",NC_FLOAT,"","m",*mask,"ocean waves dZ max"),buffer);
      }
    }
  if(write) status=poc_put_var(output,var.init("density-4b",NC_FLOAT,"","m",*mask,"ocean waves dZ max"),*density);

/*------------------------------------------------------------------------------
  ocean surface wave, CFL WW3 constraint */
//   if(criteria.surface_wave_cfl) {
  if(criteria.use(7)) {
    printf("\nmesh density, step #4 (ocean waves CFL)\n");
    status= set_density04_c(criteria, *grid, *density, topo, *mask, buffer, constrain);
    if(buffer!=0) {
      if(write) status=poc_put_var(output,var.init("density-4c-alone",NC_FLOAT,"","m",*mask,"ocean waves CFL"),buffer);
      }
    }
  if(write) status=poc_put_var(output,var.init("density-4c",NC_FLOAT,"","m",*mask,"ocean waves CFL"),*density);

/*------------------------------------------------------------------------------
  barotropic gravity wave, actual typical length constraint */
  if(criteria.use(4)) {
    printf("\nmesh density, step #5 (true tidal wavelength)\n");
    status= set_density05(criteria, *grid, *density, topo, *mask, buffer);
    }
  if(write) status=poc_put_var(output,var.init("density-5",NC_FLOAT,"","m",*mask,"true tidal wavelength"),*density);

  
/*------------------------------------------------------------------------------
  baroclinic gravity wave */
  if(criteria.use(5)) {
    printf("\nmesh density, step #6 (baroclinic wavelength)\n");
    status= set_density06(criteria, *grid, *density, topo, *mask, buffer);
    }
  if(write) status=poc_put_var(output,var.init("density-6",NC_FLOAT,"","m",*mask,"baroclinic wavelength"),*density);

  
  if(constrain!=0) {
    poc_var_t crit_num=var;
    
    crit_num.init("crit_num",NC_FLOAT,"","",*mask);
    
#define poc_att_t_MACRO(x) poc_att_t(#x "_value",x)
    crit_num
      <<poc_att_t_MACRO(MSH_MAXSIZE_OPEN)
      <<poc_att_t_MACRO(MSH_MAXSIZE_SHELF)
      <<poc_att_t_MACRO(MSH_TIDAL_WAVELENGTH)
      <<poc_att_t_MACRO(MSH_TOPO_SLOPE)
      <<poc_att_t_MACRO(MSH_SURFWAVE_LENGTH)
      <<poc_att_t_MACRO(MSH_SURFWAVE_DZ)
      <<poc_att_t_MACRO(MSH_SURFWAVE_CFL)
      ;
#undef poc_att_t_MACRO
    
    if(write) status=poc_put_var(output,crit_num,constrain);
    }
  
/*------------------------------------------------------------------------------
  projection rescaling */
  rescale_density(*grid, *density, *mask);
  
/*------------------------------------------------------------------------------
  external constraint (for mesh splitting/assembly) */
  if(criteria.resample_obsolete==0) {
    printf("\nmesh density, step #7 (model limits resolution consistency)\n");
    status= set_density03(criteria, *grid, *density, *mask, polygones, npolygones);
    if(write) status=poc_put_var(output,var.init("density-3",NC_FLOAT,"","m",*mask,"model limits resolution consistency"),*density);
    }
    
/*------------------------------------------------------------------------------
  smooth density field*/
  if(criteria.smooth==1) {
    printf("\nmesh density, step #8 (smoothing)\n");
//     status= smooth_density(*grid, *density, *mask, criteria.maxrate, criteria.niterations);
    status= smooth_density_new(*grid, *density, *mask, criteria.maxrate);
    }

  if(write) status=poc_put_var(output,var.init("density",NC_FLOAT,"","m",*mask,"final density"),*density);

  printf("\nmesh density sucessfully completed\n");
  
  deletep(&topo);
  deletep(&buffer);
  deletep(&constrain);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_ModifyMeshsize(const char *output, criteria_t criteria, grid_t grid, float *prior, float mask, plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,status;
  int i,j,k;
  grid_t topogrid;
  char *interior;
  cdfvar_t variable;
  pocgrd_t ncgrid;
  float *density;
  string merged;
  
  if(output!=0 and strlen(output)!=0) merged="merged.nc";
  else merged="";

  status=plg_cartesian((projPJ )grid.projection,polygones, npolygones);

  interior=plg_test_grid(grid,polygones,npolygones);

/*------------------------------------------------------------------------------
  initialise density with 1 for interior (polygon-wise) grid points, -1 otherwise */
  density=new float[grid.nx*grid.ny];
  
  count=0;
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      if(prior[k]!=mask) {
        density[k]=interior[k];
        if(density[k]!=mask) count++;
        }
      else density[k]=mask;
      }
    }
    
  criteria.resample_obsolete=1;
  status=fe_ComputeMeshsize(output, criteria, &grid, &density, &mask, polygones, npolygones, 0, false);
  
  if(merged!="") {
    status=poc_createfile(merged.c_str(), NC_64BIT_OFFSET);
    status=poc_sphericalgrid_xy(merged.c_str(),"",grid,&ncgrid);

    poc_standardvariable_xy(&variable, "prior",mask,"dimensionless",1., 0.,"density","density","density", ncgrid);
    status=create_ncvariable(merged.c_str(), &variable);
    status=poc_write_xy(merged.c_str(), grid, variable.id, prior);
    variable.destroy();
    }
  
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      if( (density[k]!=mask) && (prior[k]!=mask) ) {
        prior[k]=min(density[k],prior[k]);
        }
      }
    }
  
  delete[]density;
  
  if(merged!="") {
    poc_standardvariable_xy(&variable, "merged",mask,"dimensionless",1., 0.,"density","density","density", ncgrid);
    status=create_ncvariable(merged.c_str(), &variable);
    status=poc_write_xy(merged.c_str(), grid, variable.id, prior);
    variable.destroy();
    }
  
/*------------------------------------------------------------------------------
  smooth density field*/
  status= smooth_density(grid, prior, mask, criteria.maxrate);

  if(merged!="") {
    poc_standardvariable_xy(&variable, "density",mask,"dimensionless",1., 0.,"density","density","density", ncgrid);
    status=create_ncvariable(merged.c_str(), &variable);
    status=poc_write_xy(merged.c_str(), grid, variable.id, prior);
    variable.destroy();
    }

  printf("%s sucessfully completed\n",__FUNCTION__);

  return(0);
}


