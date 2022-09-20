
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

/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  rectify aims to correct bathymetric database; here antartctic region.

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

#define MAIN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "tools-structures.h"

#include "rutin.h"
#include "geo.h"
#include "polygones.h"
#include "grd.h"
#include "map.h"
#include "functions.h"

#define nstat 12

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  z;
  geo_t projection;
  float  mask=1.e+10;
  double  x,y,t,p;
  int i,j,k,l,m,n,status;
  char *keyword,*zone;
  char *rootname=NULL,*output=NULL,*input=NULL,*pathname=NULL,*bathymetry=NULL;
  char *bathy=NULL,*notebook=NULL,*meshbook=NULL,*poly=NULL;
  grid_t topogrid;
  grid_t grid;
  short *topo,smask=256*127+255;
  float *reference,*tmp;
  plg_t *polygones=NULL;
  int npolygones=0;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'b' :
          bathymetry= strdup(argv[n+1]);
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
        if(input==NULL) {
          input= strdup(argv[n]);
          n++;
          }
        else {
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
          }
        break;
      }
      free(keyword);
    }

/*-----------------------------------------------------------------------------
  read topo modifiable database */
  status=grd_loadgrid(input,&topogrid);
  if(status !=0) {
    __OUT_BASE_LINE__("cannot load bathymetry file=%f\n",input);
    exit(-1);
    }

  topo=(short *) malloc(topogrid.nx*topogrid.ny*sizeof(short));
  status= grd_loads1(input,topogrid,topogrid.nx,topo,&smask);
  topogrid.modeH=0;
/*-----------------------------------------------------------------------------
  alreday done in grd_loads1
  status= grd_mirror_r(topogrid,topogrid.nx,topo,mask); */

/*-----------------------------------------------------------------------------
  read topo reference database */
  status=grd_loadgrid(bathymetry,&grid);
  if(status !=0) {
    __OUT_BASE_LINE__("cannot load bathymetry file=%s\n",bathymetry);
    exit(-1);
    }

  reference=(float *)malloc(grid.nx*grid.ny*sizeof(float));
  status= grd_extract(bathymetry,grid,grid.nx,reference,&mask);
  grid.modeH=0;
  status= grd_mirror_r(grid,grid.nx,reference,mask);

/*------------------------------------------------------------------------------
  rectify antarctic regions (-85<phi<-60) and store topo */
  mask=-99999.;
  for(j=0;j<topogrid.ny;j++) {
    p=topogrid.ymin+j*topogrid.dy;
    if(p<-85.5) continue; /*bottom of Ross sea*/
    if(p<-60.0) {
      for(i=0;i<topogrid.nx;i++) {
        m=topogrid.nx*(topogrid.ny-j-1)+i;
/*
      m=topogrid.nx*j+i;
      t=topogrid.x[m];
      p=topogrid.y[m];
*/
        if (topo[m]==1) {
/*        if (topo[m]>0) { */
          t=topogrid.xmin+i*topogrid.dx;
          status=map_interpolation(grid,reference,mask,t,p,&z);
          if(status !=0) {
             __OUT_BASE_LINE__("interpolation error at %f %f\n",t,p);
             exit(-1);
             }
          topo[m]=(short)(z);
          }
        }
      }
    else
      break;
    }

/*------------------------------------------------------------------------------
  save rectified topo */
  status=grd_save(input,topogrid,topogrid.nx, topo, smask);

  free(topo);

  __OUT_BASE_LINE__("bathymetry sucessfully completed\n");

  exit(0);

 error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
  
}
