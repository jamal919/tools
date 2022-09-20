
/*******************************************************************************
T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative
  
E-mail: florent.lyard@legos.obs-mip.fr
*******************************************************************************/
/**
\file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada
\author  Yves Soufflet      LEGOS, Toulouse, France

VERSION :
\date  28/06/2011 : Yves Soufflet: Add input help function

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Validates tidal atlases.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>

#include <errno.h>

#include "config.h"

#include "tools-define.h"
#include "tools-structures.h"
#include "constants.h"

#include "fe.h"
#include "poc-time.h"
#include "list.h"
#include "map.h"
#include "map.def"
#include "grd.h"
#include "geo.h"
#include "filter.h"
#include "functions.h"
#include "tides.h"
#include "netcdf-proto.h"

extern int crozet(const char *, const char *, int, int);
extern int levitus(void);
extern int WOA2005(void);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : argv[0]
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    " NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    " USE\n"
    "   %s -g file.mgr -a WAVE-file.nc wave1 [ wave2 ...]\n",prog_name);
  printf("\n"
    " DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "   Validates tidal atlases.\n"
    "   More description eventually.\n"
    "\n"
    " OPTIONS :\n"
   " -b        bathymetry database\n"
   " -d        print data/models differences \n"
   " -l        output format in latex \n"
   " -m,-a     naming convention for tidal atlases\n"
   " -p        path for tidal atlases\n"
   " -o        output file name\n"
   " -c,-g     harmonic data input file\n"
   " -v        followed by 2 variables names. Default: Ha Hg\n"
   " -x        output order: FILE, ALPHA, LAT, LON, IMMERSION\n"
   " -h        Print this help\n"
   " -unstructured LGP1 or LGP2 if the netcdf file is unstructured: please provide the discretisation type"
    "\n"
    "NOTE : -m and -c are provided for backward compatibility only and are deprecated.\n"
    ); /** \endcode */

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int newAinit(T* & buffer, T spec, int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
///\returns always 0, whether something was done or not!
{
  if(buffer!=0) return(0);
  
  buffer=aset(nvalues,spec);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LengthScale(const char **varnames, char **wave, char *output, char *atlas_directory, char *atlas_convention,
                  char *bathymetry, int structured, char *discretisation, int iteration)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count;
  int i,j,k,l,m,n,fmt,status,tg;
  FILE *in=NULL,*out=NULL;

  char c;
  spectrum_t spectrum;
  date_t date;
  char *model=NULL;

  int *list=NULL;
  
  mesh_t mesh;
  grid_t grid, topogrid;
  complex<float> *tides, *dumx, *dumy, mask, z;
  float *topo=0,topomask,h,hlim=15000.,rmask;
  float kx, ky, *lamda1=0, *lamda2=0, *lamda3=0, *lamda4=0, length;
  double x,y;
  float factor=10.0,spec=1e+37;

  date.year=1950;
  date.month=1;
  date.day=1;
  date.second=0.0;
  astro_angles_t astro_angles;
  spectrum.init(initialize_tide(&astro_angles,date),wave);

  for (i=0; i<spectrum.n; i++) {
    spectrum.waves[i].init();
    printf ("wave: %10s, pulsation: %12.6f degrees/h \n", spectrum.waves[i].name,spectrum.waves[i].omega);
    }

/* *-----------------------------------------------------------------------------
  read topo database */
  if(bathymetry!=NULL) {
//    bathymetry=strdup("/home/softs/genesis/data/topography/data/legos.2.grd");
    status=grd_loadgrid(bathymetry,&topogrid);
    if(status !=0) {
      __OUT_BASE_LINE__("cannot load bathymetry file=%s\n",bathymetry);
      exit(-1);
      }
    exitIfNull(
      topo=(float *) malloc(topogrid.nx*topogrid.ny*sizeof(float))
      );
    topogrid.modeH=0;
    status=  grd_loadr1(bathymetry,topogrid,topo,&topomask);
    }
  else {
    topo=0;
    }


/* *-----------------------------------------------------------------------------
  read tidal atlas database and interpolate*/
  if(structured==1) {
    for (k=0;k<spectrum.n;k++) {
      status=tide_decode_atlasname(atlas_directory,atlas_convention,wave[k], 0, &model);
      printf("validate: treating %s wave from %s \n",wave[k],model);
      status=read_SGatlas(model,varnames,grid,tides,&mask);
      dumx=new complex<float>[grid.ny*grid.nx];
      dumy=new complex<float>[grid.ny*grid.nx];
//       if(lamda1==0) {
//         lamda1=new float[grid.ny*grid.nx];
//         for(j=0;j<grid.ny;j++) {
//           for(i=0;i<grid.nx;i++) {
//             m=j*grid.nx+i;
//             lamda1[m]=1.e+10;
//             }
//           }
//         }
//       if(lamda2==0) {
//         lamda2=new float[grid.ny*grid.nx];
//         for(j=0;j<grid.ny;j++) {
//           for(i=0;i<grid.nx;i++) {
//             m=j*grid.nx+i;
//             lamda2[m]=1.e+10;
//             }
//           }
//         }
//       if(lamda3==0) {
//         lamda3=new float[grid.ny*grid.nx];
//         for(j=0;j<grid.ny;j++) {
//           for(i=0;i<grid.nx;i++) {
//             m=j*grid.nx+i;
//             lamda3[m]=1.e+10;
//             }
//           }
//         }
//       if(lamda4==0) {
//         lamda4=new float[grid.ny*grid.nx];
//         for(j=0;j<grid.ny;j++) {
//           for(i=0;i<grid.nx;i++) {
//             m=j*grid.nx+i;
//             lamda4[m]=1.e+10;
//             }
//           }
//         }
      status=newAinit(lamda1,spec,grid.ny*grid.nx);
      status=newAinit(lamda2,spec,grid.ny*grid.nx);
      status=newAinit(lamda3,spec,grid.ny*grid.nx);
      status=newAinit(lamda4,spec,grid.ny*grid.nx);
      status=map_gradient(grid, grid.nx, tides, mask,  GEOCENTRIC, dumx, dumy);
      rmask=mask.real();
      for(j=0;j<grid.ny;j++) {
        for(i=0;i<grid.nx;i++) {
          m=j*grid.nx+i;
          if((dumx[m]!=mask) && (dumy[m]!=mask)) {
            x=map_grid_x (grid, i, j);
            y=map_grid_y (grid, i, j);
            if(topo!=0) {
              x=map_recale(topogrid,x);
              status=map_interpolation(topogrid, topo,topomask,x,y,&h);
              h=fabs(h);
              }
            else {
              h=15000;
              }
            z=dumx[m]/tides[m];
            kx=z.real();
            z=dumy[m]/tides[m];
            ky=z.real();
            length=2.*M_PI/sqrt(kx*kx+ky*ky);
            length/=(1000.*15);
            lamda1[m]=min(length,lamda1[m]);
            z=dumx[m]/tides[m];
            kx=z.imag();
//            kx=abs(z);
            z=dumy[m]/tides[m];
            ky=z.imag();
//            ky=abs(z);
            length=2.*M_PI/sqrt(kx*kx+ky*ky);
            length/=(1000.*15);
            lamda2[m]=min(length,lamda2[m]);
            length=3600.*sqrt(9.81*h)*360.0/spectrum.waves[k].omega;
            length/=(1000.*15);
            lamda3[m]=min(length,lamda3[m]);
            lamda4[m]=min(lamda1[m],lamda4[m]);
            lamda4[m]=min(lamda2[m],lamda4[m]);
            lamda4[m]=min(lamda3[m],lamda4[m]);
            }
          else {
            lamda1[m]=rmask;
            lamda2[m]=rmask;
            lamda3[m]=rmask;
            lamda4[m]=rmask;
            }
          }
        }
      delete[] dumx;
      delete[] dumy;
      }
    status=save_SG("lamda.nc", grid, lamda1, rmask, "lamda1_15", "kilometers", "tidal_length_scale", 1);
    status=save_SG("lamda.nc", grid, lamda2, rmask, "lamda2_15", "kilometers", "tidal_length_scale", 0);
    status=save_SG("lamda.nc", grid, lamda3, rmask, "lamda3_15", "kilometers", "tidal_length_scale", 0);
    status=save_SG("lamda.nc", grid, lamda4, rmask, "lamda4_15", "kilometers", "tidal_length_scale", 0);
    delete[] lamda1;
    delete[] lamda2;
    delete[] lamda3;
    delete[] lamda4;
    }
  else {
    discretisation_t z_descriptor;
    for (k=0;k<spectrum.n;k++) {
      status=tide_decode_atlasname(atlas_directory,atlas_convention,wave[k], 0, &model);
      printf("validate: treating %s wave from %s \n",wave[k],model);
      status=read_UG2Datlas(model,mesh,tides,mask,discretisation, iteration, 0);
      z_descriptor=mesh.LGP2descriptor;
      int nvalues=z_descriptor.nnodes;
      complex<float> *gx, *gy;
      status=newAinit(lamda1,spec,nvalues);
      status=newAinit(lamda2,spec,nvalues);
      status=newAinit(lamda3,spec,nvalues);
      status=newAinit(lamda4,spec,nvalues);
      gx=new complex<float>[z_descriptor.nnpe];
      gy=new complex<float>[z_descriptor.nnpe];
      for(m=0;m<mesh.ntriangles;m++) {
        fe_gradient(mesh, z_descriptor, tides, m, gx, gy);
        for(k=0;k<z_descriptor.nnpe;k++) {
          n=z_descriptor.NIbE[m][k];
          z=gx[k]/tides[n];
          kx=z.real();
          z=gy[k]/tides[n];
          ky=z.real();
          length=2.*M_PI/sqrt(kx*kx+ky*ky);
          length/=(1000.*15);
          lamda1[n]=min(length,lamda1[n]);
          z=gx[k]/tides[n];
          kx=z.imag();
          z=gy[k]/tides[n];
          ky=z.imag();
          length=2.*M_PI/sqrt(kx*kx+ky*ky);
          length/=(1000.*15);
          lamda2[n]=min(length,lamda2[n]);
//           length=3600.*sqrt(9.81*h)*360.0/spectrum.waves[k].omega;
//           length/=(1000.*15);
//           lamda3[n]=min(length,lamda3[n]);
          lamda4[n]=min(lamda1[n],lamda4[n]);
          lamda4[n]=min(lamda2[n],lamda4[n]);
//           lamda4[n]=min(lamda3[n],lamda4[n]);
          }
        }
      }
    rmask=spec;
    grid=get_zonegrid("global-HR");
    status=map_completegridaxis_2(&grid);
    int *z_elts=fe_scan_elements(mesh,grid,0);
    if(z_elts==NULL) return(-1);
    float *SGbuf=new float[grid.nx*grid.ny];
    status=fe_map(mesh, lamda1, z_descriptor.type, grid, z_elts, SGbuf, rmask);
    status=save_SG("lamda.nc", grid, SGbuf, rmask, "lamda1_15", "kilometers", "tidal_length_scale", 1);
    status=fe_map(mesh, lamda2, z_descriptor.type, grid, z_elts, SGbuf, rmask);
    status=save_SG("lamda.nc", grid, SGbuf, rmask, "lamda2_15", "kilometers", "tidal_length_scale", 0);
//     status=fe_map(mesh, lamda1, z_descriptor.type, grid, z_elts, SGbuf, rmask);
//     status=save_SG("lamda.nc", grid, SGbuf, rmask, "lamda3_15", "kilometers", "tidal_length_scale", 0);
    status=fe_map(mesh, lamda4, z_descriptor.type, grid, z_elts, SGbuf, rmask);
    status=save_SG("lamda.nc", grid, SGbuf, rmask, "lamda4_15", "kilometers", "tidal_length_scale", 0);
    delete[] lamda1;
    delete[] lamda2;
    delete[] lamda3;
    delete[] lamda4;
    }

  free(model);
  free(spectrum.waves);
    
  return(0);
error:
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gravitational_constant_01(const grid_t & grid, float *g)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Computes latitude-dependant gravity of Earth
/**
\param grid structured grid
\param *g latitude-dependant gravity of Earth
\returns always 0
*/
/*----------------------------------------------------------------------------*/
{
  int n, status;
//  printf("use latitude-dependant gravitational constant\n");
  for(n=0;n<grid.ny*grid.nx;n++) {
    g[n]=gravitational_constant_01deg(grid.y[n]);
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int Loading(const char **varnames, char **wave, char *root_output, char *atlas_directory,char *atlas_convention,
              char *loading_directory,char *loading_convention,  char *bathymetry, int structured, char *discretisation, int iteration)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count;
  int i,j,k,l,m,n,fmt,status,tg;
  FILE *in=NULL,*out=NULL;

  char c;
  spectrum_t spectrum;
  date_t date;
  char *model=NULL;
  string output;

  int *list=NULL;
  
  mesh_t mesh;
  grid_t grid, loading_grid;
  complex<float> *tides, *loading, *ratio, *smoothed, mask, loading_mask, z, zz, mean;
  float rmask;
  float kx, ky, *lamda1=0, *lamda2=0, *lamda3=0, *lamda4=0, *g=0, *gprime=0, length;
  double x,y;
  float factor=10.0;
  grid_t topogrid;
  float *topo=0,topomask,h;

  date.year=1950;
  date.month=1;
  date.day=1;
  date.second=0.0;
  astro_angles_t astro_angles;
  spectrum.init(initialize_tide(&astro_angles,date),wave);

  for (i=0; i<spectrum.n; i++) {
    spectrum.waves[i].init();
    printf ("wave: %10s, pulsation: %12.6f degrees/h \n", spectrum.waves[i].name,spectrum.waves[i].omega);
    }

/* *-----------------------------------------------------------------------------
  read topo database */
  if(bathymetry!=NULL) {
//    bathymetry=strdup("/home/softs/genesis/data/topography/data/legos.2.grd");
    status=grd_loadgrid(bathymetry,&topogrid);
    if(status !=0) {
      __OUT_BASE_LINE__("cannot load bathymetry file=%s\n",bathymetry);
      exit(-1);
      }
    exitIfNull(
      topo=(float *) malloc(topogrid.nx*topogrid.ny*sizeof(float))
      );
    topogrid.modeH=0;
    status=  grd_loadr1(bathymetry,topogrid,topo,&topomask);
    }
  else {
    topo=0;
    }

/* *-----------------------------------------------------------------------------
  read tidal atlas database and loading and interpolate*/
  if(structured==1) {
    for (k=0;k<spectrum.n;k++) {
      count=0;
      mean=0;
      output=(string)wave[k]+(string)"."+(string)root_output;
      status=tide_decode_atlasname(atlas_directory,atlas_convention,wave[k], 0,&model);
      printf("validate: treating %s wave from %s \n",wave[k],model);
      status=read_SGatlas(model,varnames,grid,tides,&mask);
      status=tide_decode_atlasname(loading_directory,loading_convention,wave[k], 0,&model);
      printf("validate: treating %s wave from %s \n",wave[k],model);
      status=read_SGatlas(model,varnames,loading_grid,loading,&loading_mask);
      ratio=new complex<float>[grid.ny*grid.nx];
      smoothed=new complex<float>[grid.ny*grid.nx];
      lamda1=new float[grid.ny*grid.nx];
      for(j=0;j<grid.ny;j++) {
        for(i=0;i<grid.nx;i++) {
          m=j*grid.nx+i;
          ratio[m]=mask;
          lamda1[m]=1.e+10;
          }
        }
      rmask=mask.real();
      for(j=0;j<grid.ny;j++) {
        for(i=0;i<grid.nx;i++) {
          m=j*grid.nx+i;
          if(tides[m]!=mask) {
            x=map_grid_x (grid, i, j);
            y=map_grid_y (grid, i, j);
            x=map_recale(loading_grid,x);
            status=map_interpolation(loading_grid, loading,loading_mask,x,y,&zz);
            if(zz==loading_mask) continue;
            if(topo!=0) {
              x=map_recale(topogrid,x);
              status=map_interpolation(topogrid, topo,topomask,x,y,&h);
              h=fabs(h);
              }
            else {
              h=15000;
              }
            z=zz/tides[m];
            ratio[m]=z;
            mean+=ratio[m];
            lamda1[m]=z.real();
            lamda1[m]=max(lamda1[m],(float) 0.01);
            lamda1[m]=min(lamda1[m],(float) 0.15);
            count++;
            }
          else {
            lamda1[m]=rmask;
            }
          }
        }
      
//      status=map_diffusion(grid, ratio, mask, 1.e+05, smoothed,30000);
      grid_t subgrid;
      complex< float > *subratio,*subsmoothed,submask;
      status=map_remap(grid, ratio, mask, &subgrid, &subratio, &submask, 8, 1);
#ifdef CHK
      {
      const char *varnames[2]  ={"ratio_a","ratio_G"};
      const char *standards[2] ={"ratio_a","ratio_G"};
      const char *units[2]={"dimensionless","dimensionless"};
      status=save_SG("subratio.nc", subgrid, subratio, submask, varnames, units, standards, 1);
      }
#endif
//      status=map_smooth_latThenLong(grid, ratio, mask, 1.e+05, smoothed);
      subsmoothed=new complex<float>[subgrid.ny*subgrid.nx];
      float scale=7.5e+05;
      status=map_smooth_latThenLong(subgrid, subratio,    submask, scale, subsmoothed);
      status=map_smooth_latThenLong(subgrid, subsmoothed, submask, scale, subsmoothed);
      status=map_smooth_latThenLong(subgrid, subsmoothed, submask, scale, subsmoothed);
      status=map_smooth_latThenLong(subgrid, subsmoothed, submask, scale, subsmoothed);
#ifdef CHK
      {
      const char *varnames[2]  ={"ratio_smoothed_a","ratio_smoothed_G"};
      const char *standards[2] ={"ratio_smoothed_a","ratio_smoothed_G"};
      const char *units[2]={"dimensionless","dimensionless"};
      status=save_SG("subratio.nc", subgrid, subsmoothed, submask, varnames, units, standards, 0);
      }
#endif
      status=map_export(subgrid, subsmoothed, submask, grid, smoothed, mask, 0);
#ifdef CHK
      {
      const char *varnames[2]  ={"ratio_a_raw","ratio_G_raw"};
      const char *standards[2] ={"ratio_a_raw","ratio_G_raw"};
      const char *units[2]={"dimensionless","dimensionless"};
      status=save_SG(output.c_str(), grid, ratio, mask, varnames, units, standards, 1);
      }
#endif
      {
      const char *varnames[2]  ={"ratio_a","ratio_G"};
      const char *standards[2] ={"ratio_a","ratio_G"};
      const char *units[2]={"dimensionless","dimensionless"};
      status=save_SG(output.c_str(), grid, smoothed, mask, varnames, units, standards, 1);
      }
      status=save_SG(output.c_str(), grid, lamda1, rmask, "factor_raw", "dimensionless", "loading_tide_ratio_raw", 0);
      for(j=0;j<grid.ny;j++) {
        for(i=0;i<grid.nx;i++) {
          m=j*grid.nx+i;
          if(smoothed[m]!=mask) {
            lamda1[m]=smoothed[m].real();
            }
          else {
            lamda1[m]=rmask;
            }
          }
        }
      status=save_SG(output.c_str(), grid, lamda1, rmask, "factor", "dimensionless", "loading_tide_ratio", 0);
      complex< float > *distant=new complex<float>[grid.ny*grid.nx];
      complex< float > *local=new complex<float>[grid.ny*grid.nx];
      for(j=0;j<grid.ny;j++) {
        for(i=0;i<grid.nx;i++) {
          m=j*grid.nx+i;
          if(tides[m]!=mask) {
            x=map_grid_x (grid, i, j);
            y=map_grid_y (grid, i, j);
            x=map_recale(loading_grid,x);
            status=map_interpolation(loading_grid, loading,loading_mask,x,y,&zz);
            if(zz==loading_mask) continue;
            distant[m]=zz-lamda1[m]*tides[m];
            local[m]=lamda1[m]*tides[m];
            }
          else {
            distant[m]=mask;
            local[m]=mask;
            }
          }
        }
      {
      const char *varnames[2]  ={"distant_a","distant_G"};
      const char *standards[2] ={"distant_a","distant_G"};
      const char *units[2]={"m","m"};
      status=save_SG(output.c_str(), grid, distant, mask, varnames, units, standards, 0);
      }
      {
      const char *varnames[2]  ={"local_a","local_G"};
      const char *standards[2] ={"local_a","local_G"};
      const char *units[2]={"m","m"};
      status=save_SG(output.c_str(), grid, local, mask, varnames, units, standards, 0);
      }
      g=new float[grid.ny*grid.nx];
      status=gravitational_constant_01(grid,g);
      gprime=new float[grid.ny*grid.nx];
      for(j=0;j<grid.ny;j++) {
        for(i=0;i<grid.nx;i++) {
          m=j*grid.nx+i;
          if(lamda1[m]!=rmask) {
            gprime[m]=(1.-lamda1[m])*g[m];
            }
          else {
            gprime[m]=rmask;
            }
          }
        }
      status=save_SG(output.c_str(), grid, g,      rmask, "g",      "m/s²", "g", 0);
      status=save_SG(output.c_str(), grid, gprime, rmask, "gprime", "m/s²", "gprime", 0);
      delete[] g;
      delete[] gprime;
      delete[] ratio;
      delete[] tides;
      delete[] loading;
      delete[] distant;
      delete[] local;
      delete[] subratio;
      delete[] subsmoothed;
      mean/=count;
      }
    delete[] lamda1;
    }
  else {
    for (k=0;k<spectrum.n;k++) {
      status=tide_decode_atlasname(atlas_directory,atlas_convention,wave[k], 0,&model);
      printf("validate: treating %s wave from %s \n",wave[k],model);
      status=read_UG2Datlas(model,mesh,tides,mask,discretisation, iteration, 0);
      }
    }



  free(model);
  free(spectrum.waves);
    
  return(0);
error:
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int smoother(const char **varnames, char **wave, char *root_output, char *atlas_directory,char *atlas_convention,
              int structured, char *discretisation, int iteration)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count;
  int i,j,k,l,m,n,fmt,status,tg;
  FILE *in=NULL,*out=NULL;

  char c;
  spectrum_t spectrum;
  date_t date;
  char *model=NULL;
  string output;

  int *list=NULL;
  
  mesh_t mesh;
  grid_t grid, loading_grid;
  complex<float> *tides, *loading, *ratio, *smoothed, mask, loading_mask, z, zz, mean;
  float rmask;
  float kx, ky, *lamda1=0, *lamda2=0, *lamda3=0, *lamda4=0, length;
  double x,y;
  float factor=10.0;
  grid_t topogrid;
  float *topo=0,topomask,h;

  date.year=1950;
  date.month=1;
  date.day=1;
  date.second=0.0;
  astro_angles_t astro_angles;
  spectrum.init(initialize_tide(&astro_angles,date),wave);

  for (i=0; i<spectrum.n; i++) {
    spectrum.waves[i].init();
    printf ("wave: %10s, pulsation: %12.6f degrees/h \n", spectrum.waves[i].name,spectrum.waves[i].omega);
    }


/* *-----------------------------------------------------------------------------
  read tidal atlas database and loading and interpolate*/
  if(structured==1) {
    for (k=0;k<spectrum.n;k++) {
      count=0;
      mean=0;
      output=(string)wave[k]+(string)"."+(string)root_output;
      status=tide_decode_atlasname(atlas_directory,atlas_convention,wave[k], 0,&model);
      printf("validate: treating %s wave from %s \n",wave[k],model);
      status=read_SGatlas(model,varnames,grid,tides,&mask);
      smoothed=new complex<float>[grid.ny*grid.nx];
      lamda1=new float[grid.ny*grid.nx];

      grid_t subgrid;
      complex< float > *subsmoothed,*subtides,submask;
      status=map_remap(grid, tides, mask, &subgrid, &subtides, &submask, 8, 1);

      subsmoothed=new complex<float>[subgrid.ny*subgrid.nx];
      float scale=4.e+05;
      status=map_smooth_latThenLong(subgrid, subtides,    submask, scale, subsmoothed);
      status=map_smooth_latThenLong(subgrid, subsmoothed, submask, scale, subsmoothed);
      status=map_smooth_latThenLong(subgrid, subsmoothed, submask, scale, subsmoothed);
      status=map_smooth_latThenLong(subgrid, subsmoothed, submask, scale, subsmoothed);
      status=map_export(subgrid, subsmoothed, submask, grid, smoothed, mask, 0);
      {
      const char *varnames[2]  ={"raw_a","raw_G"};
      const char *standards[2] ={"raw_a","raw_G"};
      const char *units[2]={"dimensionless","dimensionless"};
      status=save_SG(output.c_str(), grid, tides, mask, varnames, units, standards, 1);
      }
      {
      const char *varnames[2]  ={"smoothed_a","smoothed_G"};
      const char *standards[2] ={"smoothed_a","smoothed_G"};
      const char *units[2]={"dimensionless","dimensionless"};
      status=save_SG(output.c_str(), grid, smoothed, mask, varnames, units, standards,0);
      }
      delete[] tides;
      delete[] subtides;
      delete[] smoothed;
      delete[] subsmoothed;
      mean/=count;
      }
    delete[] lamda1;
    }
  else {
    for (k=0;k<spectrum.n;k++) {
      status=tide_decode_atlasname(atlas_directory,atlas_convention,wave[k], 0,&model);
      printf("validate: treating %s wave from %s \n",wave[k],model);
      status=read_UG2Datlas(model,mesh,tides,mask,discretisation, iteration, 0);
      }
    }

  free(model);
  free(spectrum.waves);
    
  return(0);
error:
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
  function that will be ran when the executable is started
  See the print_help <a href=#func-members>function</a> for its use.
----------------------------------------------------------------------*/
{
  int i,j,k,l,m,n,fmt,status,tg;
  char *output=NULL,*mgroutput=NULL, *keyword=NULL, *mgrfile=NULL;
  char *polygone=0;
  char *ordering=NULL,*regions=NULL;
  char *wave[101];
  int nwave=0;//wave index
  date_t date;
  char *atlas_directory=NULL,*atlas_convention=NULL;
  char show_diff=0,latex=0;

  char *bathymetry=NULL;
  char *discretization=NULL;
  int structured=1;

  const int nvars=2;
  char *varnames[nvars]={"Ha","Hg"},*discretisation;
  int defaultVars=1,iteration=-1;

  fct_echo( argc, argv);

//  fprintf(stderr,"%s  -starting computation *********\n",argv[0]);
  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        if(strcmp(keyword,"-unstructured")==0) {
          discretisation= strdup(argv[n+1]);
          structured=0;
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"--regions")==0) {
          regions= strdup(argv[n+1]);
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"--polygons")==0) {
          polygone= strdup(argv[n+1]);
          n++;
          n++;
          break;
          }
        switch (keyword[1]) {

/* *----------------------------------------------------------------------
        bathymetry database*/
        case 'b' :
          bathymetry= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        print data/models differences*/
        case 'd' :
          show_diff=1;
          n++;
          break;

/* *----------------------------------------------------------------------
        iteration index*/
        case 'i' :
          sscanf(argv[n+1],"%d",&iteration);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        output format*/
        case 'l' :
          latex=1;
          n++;
          break;

/* *----------------------------------------------------------------------
        naming convention for tidal atlases*/
        case 'm' : fprintf(stderr,"Warning: deprecated %s option. Use %s instead\n",keyword,"-a");
        case 'a' :
          atlas_convention= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        path for tidal atlases*/
        case 'p' :
          atlas_directory= strdup(argv[n+1]);
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
        output order: FILE, ALPHA, LAT, LON, IMMERSION*/
        case 'x' :
          ordering= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'v' :
          if(varnames[0]){
            if(defaultVars)
              defaultVars=0;
            else{
              __FILE_LINE__(stdout,"multiple use of -v : the last one overrides the previous.");
              for(k=0;k<nvars;k++){
                free(varnames[k]);
                }
              }
            }
          for(k=0;k<nvars;k++){
            varnames[k]=strdup(argv[n+1+k]);
            }
          n++;
          n+=nvars;
          break;

        case 'h' :
          print_help(argv[0]);
          exit(0);
        default:
          printf("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
          break;
        }
        break;

      default:
/* *----------------------------------------------------------------------
          tidal wave list*/
          wave[nwave]= strdup(argv[n]);
//          printf("input wave=%s\n",wave[nwave]);
          //spectrum.n=nwave+1;
          nwave++;
          n++;
        break;
      }
      free(keyword);
    }

  if(atlas_directory==NULL) atlas_directory=strdup(".");

/**----------------------------------------------------------------------------
  tidal vmax computation */
  const char *filename="/home/models/crozet/zoom-v2/ctbto/ORCA-4/10LayersTKEKv06-again/M2.spectral3D.nc";
  status=crozet(filename, (const char *) polygone,19,5);
//   status=crozet(filename, (const char *) polygone,19,5);
//   const char *filename="/home/models/crozet/mini-zoom-v1/ctbto/homogeneous/15LayersTKEKv06-barotropic/M2.spectral3D.nc";
//  const char *polygone="/home/models/crozet/documents/area-north.plg";
//  status=crozet(filename, (const char *) polygone,5,5);


/**----------------------------------------------------------------------------
  tidal wavelength computation */
  wave[nwave]=NULL;
  discretisation=strdup("LGP2");
  structured=0;
  iteration=-1;
//  status=LengthScale((const char**)varnames,  wave, output, atlas_directory,atlas_convention, bathymetry, structured, discretisation, iteration);

/**----------------------------------------------------------------------------
  tidal loading linear approximation */
  char *loading_convention=strdup("WAVE.load.nc");
  char *loading_directory=strdup("/home/softs/data/loading/FES2004");
//  status=Loading((const char**)varnames,  wave, "ratio.nc", atlas_directory,atlas_convention,loading_directory,loading_convention, bathymetry, structured, discretisation, iteration);

/**----------------------------------------------------------------------------
  smoother for GLORYS analyses */
  wave[0]=strdup("Ssa");
  wave[1]=strdup("Sa");
  wave[2]=0;
//  char *varnames2[2]={"sossheig_a", "sossheig_G"};
  atlas_convention=strdup("WAVE-sossheig-atlas-reformatted.nc");
  atlas_directory=strdup("/home/softs/data/climatology/GLORYS-v2");
//  status=smoother((const char**)varnames,  wave, "GLORYS.nc", atlas_directory,atlas_convention, structured, discretisation, iteration);
  
/**----------------------------------------------------------------------------
  vertical modes */
//  status=levitus();
//  status=WOA2005();
  
  __ERR_BASE_LINE__("%s -computation complete ^^^^^^^^^\n",argv[0]);
  exit(0);
error:
  __OUT_BASE_LINE__("validate: error detected, quit ... \n");
  exit(-1);
}
