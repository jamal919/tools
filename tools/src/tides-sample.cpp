
/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include "tools-define.h"
#include "tools-structures.h"

#include "poc-time.h"
#include "list.h"
#include "map.h"
#include "mgr.h"
#include "grd.h"
#include "filter.h"
#include "functions.h"
#include "tides.h"
#include "netcdf-proto.h"

extern int mgr_create(int nmgr, mgr_t ***mgr, const spectrum_t& s);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int decode_atlasname(char *atlas_directory,char *atlas_convention,char *wave, char **filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=tide_decode_atlasname(atlas_directory,atlas_convention,wave, 0,filename,1);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int interpolate(grid_t grid, int n, complex<float> *dum, fcomplex mask,
                            double x, double y, fcomplex *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
 Loading and self-attraction from diagnostic fields
 Read the amplitude (meter) and phase lag, convert into
 real and imaginary part
----------------------------------------------------------------------*/
{
  int i, j, k, l, m, status;
  int kk,ll,mm;
  int nextk[24]={-1,+1, 0, 0,-1,+1,-1,+1,-2,+2, 0, 0,-2,+2,-2,+2,-1,-1,+1,+1,-2,-2,+2,+2};
  int nextl[24]={ 0, 0,-1,+1,-1,+1,+1,-1, 0, 0,-2,+2,-1,-1,+1,+1,-2,+2,-2,+2,-2,+2,-2,+2};
  float zr, zi;
  float *amp, *pha, *zx, *zy, factor;
  float *buf[6],spec[6];
  fcomplex zz,cmask;

  status = map_interpolation(grid, dum, cmask, x, y,&zz);
  if(status != 0) {
    status=map_index( grid,  x,  y, &k, &l);
    if(status == 0) {
      mm=0;
      do {
        kk=k+nextk[mm];
        ll=l+nextl[mm];
        mm++;
        status=-1;
        if(kk<0) continue;
        if(ll<0) continue;
        if(kk>grid.nx-1) continue;
        if(ll>grid.ny-1) continue;
        x=map_grid_x(grid,kk,ll);
        y=map_grid_y(grid,kk,ll);
        status = map_interpolation(grid, dum, cmask, x, y,&zz);
        } while((status!=0)&&(mm<24));
      mm=0;
      while((status!=0)&&(mm<24)) {
        kk=k+2*nextk[mm];
        ll=l+2*nextl[mm];
        mm++;
        status=-1;
        if(kk<0) continue;
        if(ll<0) continue;
        if(kk>grid.nx-1) continue;
        if(ll>grid.ny-1) continue;
        x=map_grid_x(grid,kk,ll);
        y=map_grid_y(grid,kk,ll);
        status = map_interpolation(grid, dum, cmask, x, y,&zz);
        }
      }
    }
  *z=zz;
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int read_atlas_HUVd(char *filename, vector<mgr_t> mgr, int nmgr,double *a,double *G, double mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
 Loading and self-attraction from diagnostic fields
 Read the amplitude (meter) and phase lag, convert into
 real and imaginary part
\sa read_atlas()
----------------------------------------------------------------------*/
#warning similar to read_atlas()
{
  int i, j, k, l, m, n, items,node, status;
  int kk,ll,mm;
  int nextk[24]={-1,+1, 0, 0,-1,+1,-1,+1,-2,+2, 0, 0,-2,+2,-2,+2,-1,-1,+1,+1,-2,-2,+2,+2};
  int nextl[24]={ 0, 0,-1,+1,-1,+1,+1,-1, 0, 0,-2,+2,-1,-1,+1,+1,-2,+2,-2,+2,-2,+2,-2,+2};
  float zr, zi;
  float *buf[6],spec[6],*depth,h,rmask;
  complex<float> z,u,v,cmask;
  fcomplex *tide[3];
  FILE *in;
  double x, y, dum,A,B,I;
  char units[256];
  char *sunit;
  size_t count[3], start[3], var;
  ellipse_t e;

  grid_t grid;

  cdfgbl_t global;
  int id;
  variable_t varinfo;

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = 1;
  count[2] = 1;

  status= cdf_globalinfo(filename,&global,0);
  if(status !=0) goto error;

  id=cdf_identify(global,"Ha");
  status=cdf_loadvargrid_2d (filename,id, &grid);
  if(status !=0) goto error;

  for(k=0;k<6;k++) {
    buf[k] = new float[grid.nx*grid.ny];
    }
  depth = new float[grid.nx*grid.ny];

  for(k=0;k<3;k++) {
    tide[k] =new complex<float> [grid.nx*grid.ny];
    }

/*-----------------------------------------------------------------------
  load netcdf variable */
  id=cdf_identify(global,"Ha");
  status= cdf_loadvar_r1_2d (filename, id, 0, 0, grid, grid.nx, buf[0], &spec[0] ,&varinfo);
  id=cdf_identify(global,"Hg");
  status= cdf_loadvar_r1_2d (filename, id, 0, 0, grid, grid.nx, buf[1], &spec[1] ,&varinfo);
  id=cdf_identify(global,"Ua");
  status= cdf_loadvar_r1_2d (filename, id, 0, 0, grid, grid.nx, buf[2], &spec[2] ,&varinfo);
  id=cdf_identify(global,"Ug");
  status= cdf_loadvar_r1_2d (filename, id, 0, 0, grid, grid.nx, buf[3], &spec[3] ,&varinfo);
  id=cdf_identify(global,"Va");
  status= cdf_loadvar_r1_2d (filename, id, 0, 0, grid, grid.nx, buf[4], &spec[4] ,&varinfo);
  id=cdf_identify(global,"Vg");
  status= cdf_loadvar_r1_2d (filename, id, 0, 0, grid, grid.nx, buf[5], &spec[5] ,&varinfo);
  id=cdf_identify(global,"depth");
  status= cdf_loadvar_r1_2d (filename, id, 0, 0, grid, grid.nx, depth,  &spec[5] ,&varinfo);

  cmask=fcomplex(9999.,9999.);

  for (j=0;j<grid.ny;j++)
    for (i=0;i<grid.nx;i++) {
      n=i+grid.nx*j;
      if((buf[0][n]!=spec[0]) && (buf[1][n]!=spec[1])) {
        tide[0][n]=fcomplex(buf[0][n]*cos(buf[1][n]*d2r),-buf[0][n]*sin(buf[1][n]*d2r));
        tide[1][n]=fcomplex(buf[2][n]*cos(buf[3][n]*d2r),-buf[2][n]*sin(buf[3][n]*d2r));
        tide[2][n]=fcomplex(buf[4][n]*cos(buf[5][n]*d2r),-buf[4][n]*sin(buf[5][n]*d2r));
        }
      else {
        tide[0][n]=cmask;
        }
      }
  for(k=0;k<6;k++) {
    delete[] buf[k];
    }

  printf("node lon lat h : amp pha a b i G p\n");
  for(n = 0; n < nmgr; n++) {
    if(strcmp(mgr[n].loc.units,"degrees")==0) {
      x = mgr[n].loc.lon;
      y = mgr[n].loc.lat;
      }
    else {
      x = mgr[n].loc.lon * r2d;
      y = mgr[n].loc.lat * r2d;
      }
    if(x < grid.xmin - grid.dx / 2.)
      x = x + 360.0;
    if(x > grid.xmax + grid.dx / 2.)
      x = x - 360.0;
    status = interpolate(grid, grid.nx, tide[0], cmask, x, y,&z);
    status = interpolate(grid, grid.nx, tide[1], cmask, x, y,&u);
    status = interpolate(grid, grid.nx, tide[2], cmask, x, y,&v);
    status = map_interpolation(grid, depth,   rmask, x, y,&h);
    if(status != 0) {
      printf("interpolation error: node %d lon=%lf lat=%lf \n", n, x, y);
      a[n] =  mask;
      G[n] =  mask;
      }
    else {
      a[n] =  abs(z);
/*----------------------------------------------------------------------
      convert into phase lag*/
      G[n] = -arg(z)* r2d;
      if(G[n]<0.) G[n]+=360.;
      e=ellipse_parameter(u,v);
      }
    A=mgr[n].loc.depth*e.a;
    B=mgr[n].loc.depth*e.b;
    I=e.inclination*r2d;
    if(I<0.0) I+=180.0;
    printf("node %d lon=%lf lat=%lf h=%lf : %8.3lf %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf\n",
            n, x, y,mgr[n].loc.depth,a[n],G[n],A,B,I,e.phase*r2d,e.polarisation);
    A=h*e.a;
    B=h*e.b;
    I=e.inclination*r2d;
    if(I<0.0) I+=180.0;
    printf("node %d lon=%lf lat=%lf h=%lf : %8.3lf %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf\n",
           n, x, y,h,a[n],G[n],A,B,I,e.phase*r2d,e.polarisation);
    }

  free(grid.x);
  free(grid.y);
  for(k=0;k<3;k++) {
    delete[] tide[k];
    }
  return (status);

error:
  map_printgrid(grid);
  printf("interpolation error: node %d lon=%lf lat=%lf \n", n, x, y);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double mask=999.9;
  double *residual,*buffer;
  double  **a,**G;
  int i,j,k,l,m,n,fmt,status,tg;
  FILE *in,*out;
  char *model=NULL,*data, *input,*output=NULL, *keyword, *path=NULL,*mgrfile=NULL;
  fcomplex meanc;
  spectrum_t spectrum,reference_spectrum;
  astro_angles_t astro_angles;
  char *wave[100];
  int nwave=0;
  date_t date;
  vector<mgr_t> mgr;
  int nmgr;
  char *atlas_directory=NULL,*atlas_convention=NULL,*posfile=NULL;
  double tau;

  char *bathymetry=NULL;
  grid_t topogrid;
  float *topo,topomask,h,hlim=15000.;
  double x,y;
  float factor=10.0;

  fct_echo( argc, argv);

//  fprintf(stderr,"%s  -starting computation *********\n",argv[0]);
  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

/* *----------------------------------------------------------------------
        naming convention for tidal atlases*/
        case 'm' :
          atlas_convention= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        path for tidal atlases*/
        case 'p' :
          path= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        output file name*/
        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        positions file name*/
        case 'l' :
          posfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        harmonic data input file*/
        case 'c' :
          mgrfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
          break;
        }
        break;

      default:
/* *----------------------------------------------------------------------
          tidal wave list*/
          wave[nwave]= strdup(argv[n]);
          spectrum.n=nwave+1;
          nwave++;
          n++;
        break;
      }
      free(keyword);
    }

  if(path!=NULL) atlas_directory=strdup(path);

  date.year=1950;
  date.month=1;
  date.day=1;
  date.second=0.0;
  reference_spectrum = initialize_tide(&astro_angles,date);

  spectrum.nmax=nwave+1;
  spectrum.n=0;
  spectrum.waves=new tidal_wave[spectrum.nmax];

  for(j = 0; j < nwave; j++) {
    for(k = 0; k < reference_spectrum.n; k++) {
      if(strcmp(reference_spectrum.waves[k].name, wave[j]) == 0) {
        break;
        }
      }
    if(k != reference_spectrum.n) {
      addwave(&spectrum, reference_spectrum.waves[k]);
      }
    }

  for (i=0; i<spectrum.n; i++) {
    spectrum.waves[i].init();
//    printf ("wave: %10s, pulsation: %12.6f degrees/h \n", spectrum.waves[i].name,spectrum.waves[i].omega);
    }

  status=mgr_create(1,mgr,spectrum);
  for (i=0; i<spectrum.n; i++) {
    for (j=i+1; j<spectrum.n; j++) {
      tau=deltaOmega2separation(spectrum.waves[j].omega-spectrum.waves[i].omega);
      printf ("wave: %10s %10s, separation: %9.3f days \n",
        spectrum.waves[i].name,spectrum.waves[j].name,tau);
      }
    }

/* *-----------------------------------------------------------------------------
  read tide gauge database */
  in=fopen(posfile,"r");
  fscanf(in,"%d",&nmgr);
  status=mgr_create(nmgr,mgr,spectrum);
  for(n=0;n<nmgr;n++) {
    fscanf(in,"%lf %lf %lf",&(mgr[n].loc.lon),&(mgr[n].loc.lat),&(mgr[n].loc.depth));
    mgr[n].loc.units=strdup("degrees");
    }
  fclose(in);
  if(status==-1) {
    printf("error\n");
    goto error;
    }

  a=(double **) malloc(nwave*sizeof(double));
  G=(double **) malloc(nwave*sizeof(double));

  for (k=0;k<nwave;k++) {
    a[k]=(double *) malloc(nmgr*sizeof(double));
    G[k]=(double *) malloc(nmgr*sizeof(double));
    }

/* *-----------------------------------------------------------------------------
  read tidal atlas database */
  for (k=0;k<nwave;k++) {
    status=decode_atlasname(atlas_directory,atlas_convention,wave[k],&model);
    printf("validate: treating %s wave from %s \n",wave[k],model);
    status=read_atlas_HUVd(model,mgr,nmgr,a[k],G[k],mask);
    }

  for (k=0;k<nwave;k++) {
    free(a[k]);
    free(G[k]);
    }
  free(a);
  free(G);

  __ERR_BASE_LINE__("exiting\n");exit(0);
error:
  __OUT_BASE_LINE__("validate: error detected, quit ... \n");
  exit(-1);
}
