/**************************************************************************

  Nonlinear finite element time stepping model

  David Greenberg    Bedford Institute of Oceanography   October 1991

***************************************************************************/

#define MAIN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "tools-structures.h"
#include "constants.h"

#include "geo.h"
#include "fe.h"
#include "archive.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "netcdf-proto.h"
#include "functions.h"
#include "mass.h"
#include "matrix.h"

#include "map.def"

#include "zapper.h"     /*  rutin.h contains common utility routines  */

extern  void fe_LGP0_to_LGP1(double *in, double *out, mesh_t mesh);

double earth_radius=6.3e+06;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float correlation_distance_01(float slope,float lmax)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *----------------------------------------------------------------------
  slope of 1.e-01 is already a large slope value (100m / km)
  associated to a correlation length of 5km ?
-----------------------------------------------------------------------*/
  float tmp;

  tmp=MIN(lmax,500./slope);

  return(tmp);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int correlation_spherical(const char *output, const char *depthfile, mesh_t mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,pos,status;
  int i,j,k,m,n;
  double x,y,lon,lat;
  grid_t topogrid,spherical;
  float *topo,*topobase,topomask,h;
  float ux,uy,l;
  cdfvar_t variable;
  pocgrd_t ncgrid;
  double km2m=1000.0;
  float *dHdx,*dHdy;
  float dx,dy;
  float *direction_x,*direction_y,*slope,*length,*correlation;
  float *azimuth,*Lx,*Ly;
  float d1,d2,corr;
  float a,lmax=25000.0;

/* *----------------------------------------------------------------------
  load bathymetry grid*/
  status=grd_loadgrid(depthfile,&topogrid);
  if(status !=0) {
    printf("cannot load bathymetry file=%s\n",depthfile);
    goto error;
    }

/* *----------------------------------------------------------------------
  load bathymetry data*/
  topobase= new float[topogrid.nx*topogrid.ny];
  topogrid.modeH=0;
  status=grd_loadr1(depthfile,topogrid,topobase,&topomask);

  status=map_completegridaxis(&topogrid,2);

  dHdx=new float[topogrid.nx*topogrid.ny];
  dHdy=new float[topogrid.nx*topogrid.ny];
/* *----------------------------------------------------------------------
  returns gradient in field units (actually meters) per meters*/
  status= map_gradient(topogrid,topogrid.nx, topobase, topomask,GEOCENTRIC, dHdx, dHdy);

  status= poc_createfile(output);

  status=poc_sphericalgrid_xy(output,"",topogrid,&ncgrid);

  poc_standardvariable_xy(&variable, "bathymetry",topomask,"m",1., 0.,"bathymetry","bathymetry","bathymetry",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  topogrid, variable.id, topobase);
  variable.destroy();

/* *----------------------------------------------------------------------
  */
  poc_standardvariable_xy(&variable, "gradient_x",topomask,"dimensionless",1., 0.,"gradient_x","gradient_x","gradient_x",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  topogrid, variable.id, dHdx);
  variable.destroy();

  poc_standardvariable_xy(&variable, "gradient_y",topomask,"dimensionless",1., 0.,"gradient_y","gradient_y","gradient_y",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  topogrid, variable.id, dHdy);
  variable.destroy();

/* *----------------------------------------------------------------------
  memory allocation for buffers*/
  direction_x=new float[topogrid.nx*topogrid.ny];
  direction_y=new float[topogrid.nx*topogrid.ny];
  slope      =new float[topogrid.nx*topogrid.ny];
  length     =new float[topogrid.nx*topogrid.ny];
  azimuth    =new float[topogrid.nx*topogrid.ny];
  Lx         =new float[topogrid.nx*topogrid.ny];
  Ly         =new float[topogrid.nx*topogrid.ny];

  for(j=0;j<topogrid.ny;j++) {
    for(i=0;i<topogrid.nx;i++) {
      k=j*topogrid.nx+i;
      slope[k]=sqrt(dHdx[k]*dHdx[k]+dHdy[k]*dHdy[k]);
      if (slope[k] != 0.0) {
        direction_x[k]=dHdx[k]/slope[k];
        direction_y[k]=dHdy[k]/slope[k];
/* *----------------------------------------------------------------------
        correlation length as function of slope*/
        length[k]=correlation_distance_01(slope[k],lmax);
        azimuth[k]=atan2(direction_y[k],direction_x[k])*r2d;
        }
      else {
        direction_x[k]=0;
        direction_y[k]=0;
        length[k]=lmax;
        azimuth[k]=0;
        }
      if(topobase[k]>5.0) {
        slope[k]=topomask;
        direction_x[k]=topomask;
        direction_y[k]=topomask;
        length[k]=topomask;
        azimuth[k]=topomask;
        }
      }
    }

  for(j=0;j<topogrid.ny;j++) {
    for(i=0;i<topogrid.nx;i++) {
      k=j*topogrid.nx+i;
      a=fabs(dHdx[k]);
      if (a != 0.0) {
/* *----------------------------------------------------------------------
        correlation length as function of slope*/
        Lx[k]=correlation_distance_01(a,lmax);
        }
      else {
        Lx[k]=lmax;
        }
      a=fabs(dHdy[k]);
      if (a != 0.0) {
/* *----------------------------------------------------------------------
        correlation length as function of slope*/
        Ly[k]=correlation_distance_01(a,lmax);
        }
      else {
        Ly[k]=lmax;
        }
      if(topobase[k]>5.0) {
        Lx[k]=topomask;
        Ly[k]=topomask;
        }
      }
    }

  poc_standardvariable_xy(&variable, "direction_x",topomask,"dimensionless",1., 0.,"direction_x","direction_x","direction_x",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  topogrid, variable.id, direction_x);
  variable.destroy();

  poc_standardvariable_xy(&variable, "direction_y",topomask,"dimensionless",1., 0.,"direction_y","direction_y","direction_y",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  topogrid, variable.id, direction_y);
  variable.destroy();

  poc_standardvariable_xy(&variable, "slope",topomask,"dimensionless",1., 0.,"slope","slope","slope",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  topogrid, variable.id, slope);
  variable.destroy();

  poc_standardvariable_xy(&variable, "length",topomask,"m",1., 0.,"length","length","length",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  topogrid, variable.id, length);
  variable.destroy();

  poc_standardvariable_xy(&variable, "azimuth",topomask,"degrees",1., 0.,"azimuth","azimuth","azimuth",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  topogrid, variable.id, azimuth);
  variable.destroy();

  poc_standardvariable_xy(&variable, "L_x",topomask,"degrees",1., 0.,"L_x","L_x","L_x",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  topogrid, variable.id, Lx);
  variable.destroy();

  poc_standardvariable_xy(&variable, "L_y",topomask,"degrees",1., 0.,"L_y","L_y","L_y",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  topogrid, variable.id, Ly);
  variable.destroy();

/* *----------------------------------------------------------------------
  example*/
  correlation=new float[topogrid.nx*topogrid.ny];

/* *----------------------------------------------------------------------
  point de +forte pente*/
  pos=maxpos(slope,topogrid.nx*topogrid.ny);
  for(j=0;j<topogrid.ny;j++) {
    for(i=0;i<topogrid.nx;i++) {
      k=j*topogrid.nx+i;
/* *----------------------------------------------------------------------
      dx et dy en kilomètre*/
      dx=(topogrid.x[k]-topogrid.x[pos])*cos(topogrid.y[pos]*d2r)*earth_radius*d2r;
      dy=(topogrid.y[k]-topogrid.y[pos])*earth_radius*d2r;
      ux=direction_x[pos];
      uy=direction_y[pos];
/* *----------------------------------------------------------------------
      normalised absciss along gradient direction*/
      d1=(ux*dx+uy*dy)/length[pos];
/* *----------------------------------------------------------------------
      normalised absciss across gradient direction*/
      d2=(-uy*dx+ux*dy)/lmax;
      correlation[k]=exp(-(d1*d1+d2*d2)/2.0);  /** cf. http://fr.wikipedia.org/wiki/Loi_normale */
      }
    }

  poc_standardvariable_xy(&variable, "correlation",topomask,"dimensionless",1., 0.,"correlation","correlation","correlation",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  topogrid, variable.id, correlation);
  variable.destroy();

/* *----------------------------------------------------------------------
  compute on mesh nodes*/
  for(n=0;n<mesh.nvtxs;n++) {
    lon=mesh.vertices[n].lon;
    lat=mesh.vertices[n].lat;
    status=map_interpolation(topogrid, direction_x,topomask,lon,lat,&ux);
    status=map_interpolation(topogrid, direction_y,topomask,lon,lat,&uy);
    status=map_interpolation(topogrid, length,topomask,lon,lat,&l);
    for(k=0;k<mesh.vertices[n].nngh;k++) {
      m=mesh.vertices[n].ngh[k];
      dx=(mesh.vertices[m].lon-lon)*cos(lat*d2r)*earth_radius*d2r;
      dy=(mesh.vertices[m].lat-lat)*earth_radius*d2r;
/* *----------------------------------------------------------------------
      normalised absciss along gradient direction*/
      d1=(ux*dx+uy*dy)/l;
/* *----------------------------------------------------------------------
      normalised absciss across gradient direction*/
      d2=(-uy*dx+ux*dy)/lmax;
      corr=exp(-(d1*d1+d2*d2));
      }
    }

  printf("mesh correlation sucessfully completed\n");
  return(0);

error:
  return(-1);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int defgrid(grid_t *grid, grid_t topogrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  frame_t frame;
  int ptype=0;
  projPJ projection;
  double ref_lat[2],ref_lon[2];

  ref_lat[0]=(topogrid.ymax+topogrid.ymin)/2;
  ref_lat[0]*=M_PI/180.0;

  projection=assign_projection(ptype,ref_lat,ref_lon);

  map_completegridaxis(grid,2);

  grid->projection=projection;

  return(0);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int correlation_cartesian(char *output, char *depthfile, mesh_t mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,pos,status;
  int i,j,k,m,n;
  double x,y,lon,lat;
  grid_t topogrid,spherical;
  float *topo,*topobase,topomask,h;
  float ux,uy,l;
  cdfvar_t variable;
  pocgrd_t ncgrid;
  double km2m=1000.0;
  float *dHdx,*dHdy;
  float dx,dy;
  float *direction_x,*direction_y,*slope,*length,*correlation;
  float d1,d2,corr;
  float lmax=25000.0;

/* *----------------------------------------------------------------------
  load bathymetry grid*/
  status=grd_loadgrid(depthfile,&topogrid);
  if(status !=0) {
    printf("cannot load bathymetry file=%s\n",depthfile);
    goto error;
    }

/* *----------------------------------------------------------------------
  load bathymetry data*/
  topobase= new float[topogrid.nx*topogrid.ny];
  topogrid.modeH=0;
  status=  grd_loadr1(depthfile,topogrid,topobase,&topomask);

  printf("mesh correlation sucessfully completed\n");
  return(0);

error:
  return(-1);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nndes,status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int option,channels=0;
  FILE *file;
  FILE *out;
  char *keyword,*zone=NULL,*gridfile=NULL;
  char *meshfile=NULL,*depthfile=NULL,*output=NULL,*poly=NULL;
  mesh_t mesh;
  double dmax=0,cmax=0;
  int *selected,*targeted;
  grid_t topogrid;
  float *topo,topomask,h;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'b' :
          depthfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'd' :
          sscanf(argv[n+1],"%lf",&dmax);
          n++;
          n++;
          break;

        case 'e' :
          sscanf(argv[n+1],"%lf",&cmax);
          n++;
          n++;
          break;

        case 'c' :
          channels=1;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          poly= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        n++;
        break;
      }
    free(keyword);
    }

/* *----------------------------------------------------------------------
  load and initialize mesh structure*/
//   status=fe_init(&mesh);
//
//   if(meshfile != NULL) {
//     status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
//     if(status!=0) {
//       printf("unable to read the original mesh in %s\n",meshfile);
//       goto error;
//       }
//     status=fe_list(&mesh);
//     if(status!=0) {
//       printf("unable to build the element list from the original mesh\n");
//       goto error;
//       }
//     }
//  else {
//    printf("no mesh file specified; abort...\n");
//    goto error;
//    }
//
//   status= fe_edgetable(&mesh,0);
//   if(status!=0) {
//     printf("unable to build the edge list from the original mesh\n");
//     goto error;
//     }

/* *----------------------------------------------------------------------
  */
  status=correlation_spherical(output,depthfile, mesh);

end: __OUT_BASE_LINE__("end of mesh-topo ... \n");
  exit(0);

error:
  __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}
