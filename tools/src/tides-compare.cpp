
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

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"
#include "tides.h"
#include "map.h"
#include "netcdf-proto.h"
#include "bmg.h"
#include "grd.h"
#include "filter.h"
#include "statistic.h"
#include "functions.h"
#include "mgr.h"
#include "polygones.h"

using namespace std;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int select(grid_t grid, char *poly, signed char **selected)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,i,j,m,n,status;
  int   count=0;
  int   inside;
  edge_t *edges;
  double t,p;
  plg_t *polygones=NULL;
  int npolygones=0;
  frame_t frame;

  *selected=new signed char[grid.ny*grid.nx];

  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      (*selected)[k]=0;
      }
    }
/*-----------------------------------------------------------------------------
  area criterion */
  if(poly!=NULL) plg_load_scan(poly, &polygones, &npolygones);
  else {
    for(j=0;j<grid.ny;j++) {
      for(i=0;i<grid.nx;i++) {
        k=j*grid.nx+i;
        (*selected)[k]=1;
        }
      }
    goto end;
    }

 frame=plg_cartesian_minmax(polygones,npolygones);

  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i++) {
      n=grid.nx*j+i;
      t=grid.x[n];
      p=grid.y[n];
      if ((t<frame.xmin) || (t>frame.xmax)) continue;
      if ((p<frame.ymin) || (p>frame.ymax)) continue;
      inside=plg_TestInterior(t,p,polygones,npolygones);
      if (inside==PLG_POINT_INTERIOR) {
        count++;
        (*selected)[n]=1;
        }
      }
    }

end:
  printf ("number of point selected %d\n",count);
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double x,y;
  float  rmask,*rbufe,*rbufa,*rbufG, *buf[2],spec[2],time;
  int status;
  int i,j,k,l,m,n;
  FILE *out;
  char *output=NULL,*keyword,wave[10],*bathymetry=NULL,*mgrfile=NULL,*zone=NULL,*poly=NULL;
  char *input[10]={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  grid_t grid[2];
  grid_t grid3d,topogrid;
  cdfgbl_t global[10];
  variable_t varinfo;
  fcomplex *tide[10]={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL},cmask[10];
  fcomplex *delta[3]={NULL,NULL,NULL};
  fcomplex *cbuf=NULL,z[2],e;
  int fmt=0;
  int variable;
  int linear=0,count=0;
  float a,G;
  float *topo,topomask,h;
  statistic_c_t s4;
  statistic_t s1,s2,s3;
  vector<mgr_t> mgr;
  int nmgr=0;
  signed char *selected;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'b' :
          bathymetry= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'w' :
          sscanf(argv[n+1],"%s",wave);
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

        case 'p' :
          poly= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        input[count]= strdup(argv[n]);
        count++;
        n++;
        break;
      }
      free(keyword);
    }

  if(input[1] == NULL) goto error;

/* *-----------------------------------------------------------------------------
  read tide gauge database */
  if(mgrfile!=0) {
    nmgr=mgr_load(mgrfile, mgr);
    if(nmgr==-1) {
      printf("error in reading %s\n",mgrfile);
      goto error;
      }
    }

/* *-----------------------------------------------------------------------
  load tidal solutions*/
  for(k=0;k<count;k++) {
    status= cdf_globalinfo(input[k],&global[k],0);
    if(status!=0) fmt=1;
    switch (fmt) {
      case 0:
        status= cdf_globalinfo(input[k],&global[k],0);
        if(status !=0) goto error;
        variable=cdf_identify(global[k],"Ha");
        status=cdf_loadvargrid_2d (input[k],variable, &grid3d);
        if(status !=0) goto error;
        grid[k]= map_getgrid2d(grid3d);
        buf[0]  = new float[grid[k].nx*grid[k].ny];
        buf[1]  = new float[grid[k].nx*grid[k].ny];
        tide[k] = new fcomplex[grid[k].nx*grid[k].ny];
/*-----------------------------------------------------------------------
        load netcdf variable */
        variable=cdf_identify(global[k],"Ha");
        status= cdf_loadvar_r1_2d (input[k], variable, 0, 0, grid3d, grid3d.nx, buf[0], &spec[0] ,&varinfo);
        variable=cdf_identify(global[k],"Hg");
        status= cdf_loadvar_r1_2d (input[k], variable, 0, 0, grid3d, grid3d.nx, buf[1], &spec[1] ,&varinfo);
        cmask[k]=fcomplex(9999.,9999.);
        for (j=0;j<grid[k].ny;j++)
          for (i=0;i<grid[k].nx;i++) {
            n=i+grid[k].nx*j;
            if((buf[0][n]!=spec[0]) && (buf[1][n]!=spec[1])) {
              tide[k][n]=fcomplex(buf[0][n]*cos(buf[1][n]*d2r),-buf[0][n]*sin(buf[1][n]*d2r));
              }
            else {
              tide[k][n]=cmask[0];
              }
            }
        free(buf[0]);
        free(buf[1]);
        break;

      case 1:
        status= bmg_loadc1(input[k],1,1,1,&grid[k],&(tide[k]),&(cmask[k]),&time);
        break;
      }
    }

/*-----------------------------------------------------------------------------
  read topo database */
  if(bathymetry==NULL) bathymetry=strdup("/home/softs/genesis/data/topography/data/gridone.grd");
  status=grd_loadgrid(bathymetry,&topogrid);
  if(status !=0) {
    __OUT_BASE_LINE__("cannot load bathymetry file=%f\n",bathymetry);
    exit(-1);
    }
  topo= new float[topogrid.nx*topogrid.ny];
  topogrid.modeH=0;
  status=  grd_loadr1(bathymetry,topogrid,topo,&topomask);

  if(fabs(topogrid.xmax-180.0)<1.e-03) topogrid.xmax=(double) 180.0;
  if(fabs(topogrid.ymax-90.0)<1.e-03)  topogrid.ymax=(double) 90.0;

  rbufe  =new float[grid[1].nx*grid[1].ny];
  rbufa  =new float[grid[1].nx*grid[1].ny];
  rbufG  =new float[grid[1].nx*grid[1].ny];

  tide[2] =new fcomplex[grid[1].nx*grid[1].ny];
  delta[0]=new fcomplex[grid[1].nx*grid[1].ny];

  status= select(grid[1], poly, &selected);

/*-----------------------------------------------------------------------------
  compute differences*/
  for (j=0;j<grid[1].ny;j++)
  for (i=0;i<grid[1].nx;i++) {
    n=i+grid[1].nx*j;
    if(grid[1].modeH==0) m=j;
    else m=n;
    y=grid[1].y[m];
    if(grid[1].modeH==0) m=i;
    else m=n;
    x=map_recale(grid[0],grid[1].x[m]);
    z[1]=tide[1][n];
    status=map_interpolation(grid[0], tide[0], cmask[0], x,y,&(z[0]));
    x=map_recale(topogrid,x);
    status=map_interpolation(topogrid, topo,topomask,x,y,&h);
    if((z[0]!=cmask[0]) && (z[1]!=cmask[1])&& (h!=topomask)&& (selected[n]==1)) {
      tide[2][n]=z[1]-z[0];
      a=abs(tide[2][n]);
      G=arg(tide[2][n]);
      rbufe[n]=a;
      rbufa[n]=abs(z[1])-abs(z[0]);
      rbufG[n]=(arg(z[1])-arg(z[0]))*r2d;
      if(rbufG[n] >  180.0) rbufG[n]-=180.;
      if(rbufG[n] < -180.0) rbufG[n]+=180.;
      }
    else {
      tide[2][n]=cmask[1];
      rbufe[n]=9999.;
      rbufa[n]=9999.;
      rbufG[n]=9999.;
      }
    }
  for (n=0;n<grid[1].ny*grid[1].nx;n++) {
    delta[0][n]=tide[2][n];
    }

  if(output==0) {
    output=new char[1024];
    if(wave!=0) {
      sprintf(output,"%s.vector-difference.nc",wave);
      }
    else {
      sprintf(output,"vector-difference.nc");
      }
    }

//  status= bmg_savec1(output,1,1,1,grid[1],tide[2],time, cmask[1]);
  
  status=tides_savec1(output, "a", "G", "vector difference", "m", grid[1],tide[2], cmask[1]);
  
  printf("# comparisons : norm, amplitude, phase, complex error\n");

  printf("# global comparisons------------------------------\n");
  printf("total difference,     ");
  s1=get_statistics(rbufe, 9999.,grid[1].nx*grid[1].ny,1);
  printf("amplitude difference, ");
  s2=get_statistics(rbufa, 9999.,grid[1].nx*grid[1].ny,1);
  printf("phase lag difference, ");
  s3=get_statistics(rbufG, 9999.,grid[1].nx*grid[1].ny,1);

  s4=cget_statistics(tide[2], cmask[1],grid[1].nx*grid[1].ny);
  s4=cget_geostatistics(grid[1],tide[2], cmask[1]);

/*-----------------------------------------------------------------------------
  compute differences*/
  for (j=0;j<grid[1].ny;j++)
  for (i=0;i<grid[1].nx;i++) {
    n=i+grid[1].nx*j;
    if(grid[1].modeH==0) m=j;
    else m=n;
    y=grid[1].y[m];
    if(grid[1].modeH==0) m=i;
    else m=n;
    x=map_recale(topogrid,grid[1].x[m]);
    status=map_interpolation(topogrid, topo,topomask,x,y,&h);
    if(abs(h)<500.) {
      tide[2][n]=cmask[1];
      rbufe[n]=9999.;
      rbufa[n]=9999.;
      rbufG[n]=9999.;
      }
    }

  printf("# h>500m comparisons ------------------------------\n");
  s1=get_statistics(rbufe, 9999.,grid[1].nx*grid[1].ny,1);
  s2=get_statistics(rbufa, 9999.,grid[1].nx*grid[1].ny,1);
  s3=get_statistics(rbufG, 9999.,grid[1].nx*grid[1].ny,1);
  s4=cget_statistics(tide[2], cmask[1],grid[1].nx*grid[1].ny);
  s4=cget_geostatistics(grid[1],tide[2], cmask[1]);

/*-----------------------------------------------------------------------------
  compute differences*/
  for (j=0;j<grid[1].ny;j++)
  for (i=0;i<grid[1].nx;i++) {
    n=i+grid[1].nx*j;
    if(grid[1].modeH==0) m=j;
    else m=n;
    y=grid[1].y[m];
    if(grid[1].modeH==0) m=i;
    else m=n;
    x=map_recale(topogrid,grid[1].x[m]);
    status=map_interpolation(topogrid, topo,topomask,x,y,&h);
    if(abs(h)<1000.) {
      tide[2][n]=cmask[1];
      rbufe[n]=9999.;
      rbufa[n]=9999.;
      rbufG[n]=9999.;
      }
    }

  printf("# h>1000m comparisons------------------------------\n");
  s1=get_statistics(rbufe, 9999.,grid[1].nx*grid[1].ny,1);
  s2=get_statistics(rbufa, 9999.,grid[1].nx*grid[1].ny,1);
  s3=get_statistics(rbufG, 9999.,grid[1].nx*grid[1].ny,1);
  s4=cget_statistics(tide[2], cmask[1],grid[1].nx*grid[1].ny);
  s4=cget_geostatistics(grid[1],tide[2], cmask[1]);

/*-----------------------------------------------------------------------------
  compute differences*/
  for (n=0;n<grid[1].ny*grid[1].nx;n++) {
    tide[2][n]=delta[0][n];
    }

  for (j=0;j<grid[1].ny;j++)
  for (i=0;i<grid[1].nx;i++) {
    n=i+grid[1].nx*j;
    if(grid[1].modeH==0) m=j;
    else m=n;
    y=grid[1].y[m];
    if(grid[1].modeH==0) m=i;
    else m=n;
    x=map_recale(topogrid,grid[1].x[m]);
    if(tide[2][n]==cmask[1]) continue;
    rbufe[n]=9999.;
    rbufa[n]=9999.;
    rbufG[n]=9999.;
    status=map_interpolation(topogrid, topo,topomask,x,y,&h);
    if(status !=0) {
      printf("cannot interpolate bathymetry\n");
      goto error;
      }
    if(abs(h)>500.) {
      tide[2][n]=cmask[1];
      rbufe[n]=9999.;
      rbufa[n]=9999.;
      rbufG[n]=9999.;
      }
    else {
      a=abs(tide[2][n]);
      G=arg(tide[2][n]);
      rbufe[n]=a;
      rbufa[n]=9999.;
      rbufG[n]=9999.;
      }
    }

  printf("# h<500 comparisons------------------------------\n");
  s1=get_statistics(rbufe, 9999.,grid[1].nx*grid[1].ny,1);
  s2=get_statistics(rbufa, 9999.,grid[1].nx*grid[1].ny,1);
  s3=get_statistics(rbufG, 9999.,grid[1].nx*grid[1].ny,1);
  s4=cget_statistics(tide[2], cmask[1],grid[1].nx*grid[1].ny);
  s4=cget_geostatistics(grid[1],tide[2], cmask[1]);

  for (n=0;n<nmgr;n++) {
    x=mgr[n].loc.lon;
    y=mgr[n].loc.lat;
    x=map_recale(grid[0],x);
    status=map_interpolation(grid[0], tide[0], cmask[0], x,y,&(z[0]));
    x=map_recale(grid[1],x);
    status=map_interpolation(grid[1], tide[1], cmask[1], x,y,&(z[1]));
    if((z[0]!=cmask[0]) && (z[1]!=cmask[1])) {
      tide[2][n]=z[1]-z[0];
      a=abs(tide[2][n]);
      G=arg(tide[2][n]);
      rbufe[n]=a;
      rbufa[n]=abs(z[1])-abs(z[0]);
      rbufG[n]=(arg(z[1])-arg(z[0]))*r2d;
      if(rbufG[n] >  180.0) rbufG[n]-=180.;
      if(rbufG[n] < -180.0) rbufG[n]+=180.;
      }
    else {
      tide[2][n]=cmask[1];
      rbufe[n]=9999.;
      rbufa[n]=9999.;
      rbufG[n]=9999.;
      }
    }

  if(nmgr!=0) {
    printf("# TG positions comparisons------------------------------\n");
    s1=get_statistics(rbufe, 9999.,nmgr,1);
    s2=get_statistics(rbufa, 9999.,nmgr,1);
    s3=get_statistics(rbufG, 9999.,nmgr,1);
    s4=cget_statistics(tide[2], cmask[1],nmgr);
    }

end: __OUT_BASE_LINE__("end of compare ... \n");
  exit(0);
error:
  __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}
