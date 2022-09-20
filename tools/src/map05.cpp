
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

#include <config.h>

#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include "tools-structures.h"
#include "tools-define.h"

#include "map.h"
#include "grd.h"
#include "ascii.h"
#include "netcdf-proto.h"
#include "poc-time.h"
#include "map.def"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_LengthScale(grid_t grid, complex<float> *tides, complex<float> mask, float* & lamda1, float* & lamda2, float* & lamda3)

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
  grid_t topogrid;
  complex<float> *dumx, *dumy, z;
  float *topo=0,topomask,h,hlim=15000.,rmask;
  float kx, ky, length;
  double x,y;
  float factor=10.0;

  dumx=new complex<float>[grid.ny*grid.nx];
  dumy=new complex<float>[grid.ny*grid.nx];
  if(lamda1==0) {
    lamda1=new float[grid.ny*grid.nx];
    for(j=0;j<grid.ny;j++) {
      for(i=0;i<grid.nx;i++) {
        m=j*grid.nx+i;
        lamda1[m]=1.e+10;
        }
      }
    }
  if(lamda2==0) {
    lamda2=new float[grid.ny*grid.nx];
    for(j=0;j<grid.ny;j++) {
      for(i=0;i<grid.nx;i++) {
        m=j*grid.nx+i;
        lamda2[m]=1.e+10;
        }
      }
    }

  if(lamda3==0) {
    lamda3=new float[grid.ny*grid.nx];
    for(j=0;j<grid.ny;j++) {
      for(i=0;i<grid.nx;i++) {
        m=j*grid.nx+i;
        lamda3[m]=1.e+10;
        }
      }
    }

  status=map_gradient(grid, grid.nx, tides, mask,  GEOCENTRIC, dumx, dumy);
  rmask=mask.real();
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      if((dumx[m]!=mask) && (dumy[m]!=mask)) {
        x=map_grid_x (grid, i, j);
        y=map_grid_y (grid, i, j);
        z=dumx[m]/tides[m];
        kx=z.real();
        z=dumy[m]/tides[m];
        ky=z.real();
        length=2.*M_PI/sqrt(kx*kx+ky*ky);
        length/=(1000.*15);
        lamda1[m]=MIN(length,lamda1[m]);
        z=dumx[m]/tides[m];
        kx=z.imag();
//        kx=abs(z);
        z=dumy[m]/tides[m];
        ky=z.imag();
//        ky=abs(z);
        length=2.*M_PI/sqrt(kx*kx+ky*ky);
        length/=(1000.*15);
        lamda2[m]=MIN(length,lamda2[m]);
        lamda3[m]=min(lamda1[m],lamda3[m]);
        lamda3[m]=min(lamda2[m],lamda3[m]);
        }
      else {
        lamda1[m]=rmask;
        lamda2[m]=rmask;
        lamda3[m]=rmask;
        }
      }
    }
  delete[] dumx;
  delete[] dumy;
    
  return(0);
error:
  return(-1);
}

