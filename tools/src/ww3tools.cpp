
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

#include "tools-structures.h"

#include "rutin.h"
#include "geo.h"
#include "polygones.h"
#include "ascii.h"
#include "netcdf-proto.h"
#include "grd.h"
#include "map.h"
#include "sym-io.h"

extern void ww3_parse_notebook_grid (const char *filename);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ww3_topo(grid_t grid, const char *bathymetry, const char *poly)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,status;
  float  z;
  geo_t projection;
  float  topo_mask=1.e+10,mask=1.e+10;
  double  x,y;
  size_t size;
  FILE *in,*out;
  char *keyword,*zone,*s;
  char *rootname=NULL,*output=NULL,*input=NULL,*pathname=NULL;
  char *bathy=NULL,*notebook=NULL,*gridfile=NULL;
  grid_t sgrid;
  grid_t topogrid;
  float  *landmask,*topo,*topobase;
  plg_t *polygones=NULL;
  int npolygones=0;
  float *resolution;
  
  if(rootname==NULL) rootname= strdup("wave");

/*-----------------------------------------------------------------------------
  set working grid */

  if((zone == NULL) && (gridfile == NULL) && (notebook == NULL)) goto error;
  if((zone != NULL) && (gridfile != NULL)) goto error;
  if((zone != NULL) && (notebook != NULL)) goto error;
  if((gridfile != NULL) && (notebook != NULL)) goto error;

  if(zone != NULL) {
/*-----------------------------------------------------------------------------
    get grid from pre-compiled zone*/
    sgrid=get_zonegrid(zone);
    sgrid.modeH=0;
    status=map_completegridaxis_2(&sgrid);
    }
  else {
/*-----------------------------------------------------------------------------
    get grid from pre-built grid file*/
    zone=strdup("grid");
    status=cdf_loadvargrid (gridfile,0,&sgrid);
    }

/*-----------------------------------------------------------------------------
  get landmask from polygons */

  sprintf(input,"%s.plg",poly);
  if(poly!=NULL) {
    sprintf(input,"%s.plg",poly);
    status=plg_load_scan(input, &polygones, &npolygones);
    if(status !=0) return(-1);
    }
  else __ERR_BASE_LINE__("exiting\n");exit(-1);
  landmask=set_landmask02(sgrid,polygones,npolygones);
  printf("land mask sucessfully completed\n");

/*-----------------------------------------------------------------------------
  read topo database */
  status=grd_loadgrid(bathymetry,&topogrid);
  if(status !=0) return(-1);

  topobase= new float[topogrid.nx*topogrid.ny];
  topogrid.modeH=0;
  status= grd_loadr1(bathymetry,topogrid,topobase,&topo_mask);

/*------------------------------------------------------------------------------
  Interpolate and store topo on symphonie grid*/
//  status=copy_grid_to_grid1d(topogrid,&topogrid1d);

  mask=-99999.;
  topo=(float *)malloc(sgrid.nx*sgrid.ny*sizeof(float));
  for(j=0;j<sgrid.ny;j++)
    for(i=0;i<sgrid.nx;i++) {
      k=j*sgrid.nx+i;
      x=sgrid.x[k];
      y=sgrid.y[k];
      status=map_interpolation(topogrid, topobase,topo_mask,x,y,&z);
      if (z!=topo_mask)
        topo[k]=z;
      else
        topo[k]=mask;
      }
  printf("bathymetry sucessfully completed\n");

  sprintf(output,"%s.spherical.nc",rootname);
  status=symio_createfile(output,(size_t) sgrid.nx,(size_t) sgrid.ny, sgrid);
  status=symio_writefile (output, landmask, topo, 0, 0);

  for(j=0;j<sgrid.ny;j++)
    for(i=0;i<sgrid.nx;i++) {
      k=j*sgrid.nx+i;
      if (topo[k]!=mask) {
/*------------------------------------------------------------------------------
        set ocean positive depth*/
        topo[k]=-topo[k];
/*------------------------------------------------------------------------------
        force 0 elevation on land*/
        if (landmask[k]==-1.) topo[k]=0.0;
/*------------------------------------------------------------------------------
        fix land point with ocean depth*/
        if ((landmask[k]==-1.) && (topo[k]>0.0))  topo[k]=0.0;
/*------------------------------------------------------------------------------
        fix ocean point with land elevation*/
        if ((landmask[k]==1.) && (topo[k]<2.5))  topo[k]=2.5;
        }
      }

  sprintf(output,"%s.shame_on_me.asc",rootname);
  status= ascii_saver1 (output, sgrid, topo, mask,"%5.0f ",sgrid.nx);

  free(topo);
  free(landmask);

  return(0);

error:
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  z;
  geo_t projection;

  float  topo_mask=1.e+10,mask=1.e+10;
  double  x,y;

  int i,j,k,l,m,n,status;

  size_t size;
  FILE *in,*out;
  char *keyword,*zone,*s;

  char *rootname=NULL,*output=NULL,*input=NULL,*pathname=NULL,*bathymetry=NULL;
  char *bathy=NULL,*notebook=NULL,*poly=NULL,*gridfile=NULL;

  grid_t sgrid,topogrid1d;
  grid_t topogrid;

  float  *landmask,*topo,*topobase,*tmp;
  plg_t *polygones=NULL;
  int npolygones=0;

  fct_echo(argc,argv);
 
  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1])
        {
        case 'r' :
          rootname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'd' :
          pathname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'b' :
          bathymetry= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'g' :
          gridfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'n' :
          notebook= strdup(argv[n+1]);
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
        __OUT_BASE_LINE__("unknown option %s\n",keyword);
        exit(-1);
        break;
      }
      free(keyword);
    }

  input=(char *)malloc(1024);
  output=(char *)malloc(1024);

  if(rootname==NULL) rootname= strdup("wave");

  ww3_parse_notebook_grid ((const char*) notebook);


/*-----------------------------------------------------------------------------
  set working grid */

  if((zone == NULL) && (gridfile == NULL) && (notebook == NULL)) goto error;
  if((zone != NULL) && (gridfile != NULL)) goto error;
  if((zone != NULL) && (notebook != NULL)) goto error;
  if((gridfile != NULL) && (notebook != NULL)) goto error;

  if(zone != NULL) {
/*-----------------------------------------------------------------------------
    get grid from pre-compiled zone*/
    sgrid=get_zonegrid(zone);
    sgrid.modeH=0;
    status=map_completegridaxis_2(&sgrid);
    }
  else {
/*-----------------------------------------------------------------------------
    get grid from pre-built grid file*/
    zone=strdup("grid");
    status=cdf_loadvargrid (gridfile,0,&sgrid);
    }

  status= ww3_topo(sgrid, (const char *) bathymetry, (const char *) poly);

  sprintf(output,"%s.shame_on_me.asc",rootname);
  status= ascii_saver1 (output, sgrid, topo, mask,"%5.0f ",sgrid.nx);

  free(topo);
  free(landmask);

  __ERR_BASE_LINE__("exiting\n");exit(0);

error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
}
