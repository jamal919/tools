#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "matrix.h"

#include "geo.h"
#include "map.h"
#include "grd.h"
#include "bmg.h"
#include "polygones.h"
#include "tides.h"
#include "tides.def"
#include "filter.h"
#include "statistic.h"
#include "rutin.h"

#define nstat 12

//theses functions are defined is mss-io.cpp CPP file
extern int mss_createfile(char *filename,size_t x_len, size_t y_len, grid_t grid);
extern int writefile(char *filename, float *mss, float *cls, float *ino);
extern int savenodes(char *filename, serie_t metadata);
extern serie_t load_metadata_ref(char *input, int point);
extern int decale_metadata(data_t *data,int *keep,int nmes);
extern serie_t load_metadata_raw(char *input,  plg_t *polygones, int npolygones, char *zone);
extern void write_list_header(FILE *out,int n);
extern void save_meantracks_CLS(serie_t TPdata,char *track,float *mss,plg_t *polygones,int npolygones);
extern void save_meantracks_CTOH(serie_t TPdata,char *sat, char *track,float *mss,plg_t *polygones,int npolygones,char *root);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int compute_annual(float* h, double *t, float mask, int count)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

{
  int    i,j,k,n,status;
  double *time,*serie,*residuals,mean;
  spectrum_t s,solved;
  statistic_t sdum;

/*-----------------------------------------------------------------------------
  remove annual signal*/

  status=spectrum_init(&s);
  status=addwave(&s, wSa);
  spectrum_terminate(s);

  serie    =(double *)malloc(count*sizeof(double));
  residuals=(double *)malloc(count*sizeof(double));
  time=(double *)malloc(count*sizeof(double));

  n=0;
  for(k=0;k<count;k++) {
    if(h[k]!=mask) {
      serie[n]=h[k];
      time[n]=t[k];
      n++;
      }
    }

  mask=1.e+35;
  harmonic_analysis_with_parsing(serie,mask,time,residuals,&mean,n, s,solved, 0);
 
  sdum=get_statistics(serie,     (double) 9999.9, n);
  sdum=get_statistics(residuals, (double) 9999.9, n);

  free(serie);
  free(residuals);
  free(time);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t get_mssgrid(char *zone)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  grid_t grid;
  int zone_initialised=0;

  grid.xmin=    0.0 ;
  grid.ymin=  -80.0;
  grid.xmax= +360.0 ;
  grid.ymax=  +82.0;
  grid.dx  =    1./30.;
  grid.dy  =    1./30.;
  grid.nx  = 10801;
  grid.ny  =  4861;
  zone_initialised=1;

  return(grid);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t get_zonegrid(char *zone)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  grid_t grid;
  int zone_initialised=0;

/* NEA grid */
  if(strcmp(zone,"NEA")==0) {
  map_set2Dgrid(&grid,-20.0,+29.5,+15.0,+65.0,1./10.,1./10.);
  zone_initialised=1;
  }

/* NEA_hr grid haute resolution */
  if(strcmp(zone,"NEA_hr")==0) {
  map_set2Dgrid(&grid,-20.0,+29.5,+15.0,+65.0,1./50.,1./50.);
  zone_initialised=1;
  }

/* Shom grille de REMY BARAILLE */
  if(strcmp(zone,"shom_remy")==0) {
  map_set2Dgrid(&grid,-16.0,+42.0,+4.0,+52.0,1./60.,1./60.);
  zone_initialised=1;
  }

/* Manche_hr grid haute resolution */
  if(strcmp(zone,"Manche_hr")==0) {
  map_set2Dgrid(&grid,-6.0,+6.0,+48.0,+53.0,1./100.,1./100.);
  zone_initialised=1;
  }

 /* NEA_shom grid distribution shom  haute resolution */
  if(strcmp(zone,"NEA_shom")==0) {
  map_set2Dgrid(&grid,-20.0,+40.0,+15.0,+55.0,1./30.,1./30.);
  zone_initialised=1;
  }
      
/* Iroise grid */
  if(strcmp(zone,"iroise")==0) {
  map_set2Dgrid(&grid,-6.0,+47.25,-3.5,+49.0,1./60.,1./60.);
  zone_initialised=1;
  }
    
/* global grid */
  if(strcmp(zone,"global")==0) {
  map_set2Dgrid(&grid,0.0,-80.0,+360.0,+80.0,1./4.,1./4.);
  zone_initialised=1;
  }
  
/* global grid */
  if(strcmp(zone,"global-loren")==0) {
  map_set2Dgrid(&grid,0.0,-90.0,+360.0,+90.0,1.,1.);
  zone_initialised=1;
  }
  
    
/* mesdsea grid, normal resolution */
  if(strcmp(zone,"medsea")==0) {
  map_set2Dgrid(&grid,-10.0,27.5,40.0,47.5,1./10.,1./10.);
  zone_initialised=1;
  }

/* alobran grid, normal resolution */
  if(strcmp(zone,"alboran")==0) {
  map_set2Dgrid(&grid,-9.0,32,-2+1./300,38+1./300,1./60.,1./60.);
  zone_initialised=1;
  }
  
/* alobran grid, normal resolution */
  if(strcmp(zone,"strait-of-sicily")==0) {
  map_set2Dgrid(&grid,-9.0,32,-2+1./300,38+1./300,1./60.,1./60.);
  zone_initialised=1;
  }
 
/* Peru/chile grid, normal resolution */
  if(strcmp(zone,"amsud")==0) {
  map_set2Dgrid(&grid,-115.0,-50.0,-70.0,0,1./30.,1./30.);
  zone_initialised=1;
  }

/* caspian grid, normal resolution, LR, add 04/07/05 */
 if(strcmp(zone,"caspian")==0) {
    /* Frame global */
    map_set2Dgrid(&grid,46.0,36.0,55.0,47.5,0.05,0.05);
    grid.modeH=0;
    zone_initialised=1;
    }

  if(zone_initialised!=1) {
    __OUT_BASE_LINE__("no valid region specified (%s); abort...\n",zone);
    exit(-1);
    }
  /* LR, modif 04/07/05 */

  return(grid);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t get_trackgrid(double *lon,double *lat)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  grid_t grid;
  grid_t cgrid;
  double dt,dp,dn;
  double tx,ty,nx,ny;
  int i,j,k,dx,dy;
  double t0,p0;


  int zone_initialised=0;

/* cartesian, rotated grid */

  cgrid.ymin=   +0.0;
  cgrid.ymax= +200.0 ;
  cgrid.xmin=   -2.0 ;
  cgrid.xmax=   +2.0 ;
  cgrid.dx  =    1.0;
  cgrid.dy  =    2.0;

  dt=lon[1]-lon[0];
  dp=lat[1]-lat[0];
  dn=sqrt(dt*dt+dp*dp);
  dt/=dn;
  dp/=dn;

  tx=dt;
  ty=dp;

  nx=-dp;
  ny= dt;

  grid.nx  =   9;
  dx=4;
  grid.dx  =  0.5/100.;

  grid.nx  =   5;
  dx=2;
  grid.dx  =  1.25/100.;

  grid.ny  = 200;
  dy=20;
  grid.dy  =  2.0/100.;

  grid.modeH  = 2;

  t0=lon[0]+dx*nx*grid.dx-dy*tx*grid.dy;
  p0=lat[0]+dx*ny*grid.dx-dy*ty*grid.dy;

  grid.x=(double *)malloc(grid.nx*grid.ny*sizeof(double));
  grid.y=(double *)malloc(grid.nx*grid.ny*sizeof(double));

  grid.ymin= +1.e+35;
  grid.ymax= -1.e+35;
  grid.xmin= +1.e+35;
  grid.xmax= -1.e+35;

  for (j=0;j<grid.ny;j++)
    for (i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      grid.x[k]=t0-i*nx*grid.dx+j*tx*grid.dy;
      grid.y[k]=p0-i*ny*grid.dx+j*ty*grid.dy;
      grid.ymin=MIN(grid.ymin,grid.y[k]);
      grid.ymax=MAX(grid.ymax,grid.y[k]);
      grid.xmin=MIN(grid.xmin,grid.x[k]);
      grid.xmax=MAX(grid.xmax,grid.x[k]);
      }

  zone_initialised=1;
   
  if(zone_initialised!=1) {__ERR_BASE_LINE__("exiting\n");exit(-1);}

  return(grid);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void get_extremities(char *sat,char *xover,serie_t TPdata,double *lon,double *lat)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,count;
  double x,y;

  sscanf(xover,"%d",&k);

  if (strcasecmp(sat,"GFO")!=0)
    /*-----------------------------------------------------------------------------
      all satellites except GFO */
    {
      if(k%2==0) /*descending*/
        {
          lon[0]=+1.e+35;
          lat[0]=-1.e+35;
          lon[1]=-1.e+35;
          lat[1]=+1.e+35;
          for(k=0;k<TPdata.count;k++) {
              lon[0]=MIN(TPdata.data[k].lon,lon[0]);
              lon[1]=MAX(TPdata.data[k].lon,lon[1]);
              lat[0]=MAX(TPdata.data[k].lat,lat[0]);
              lat[1]=MIN(TPdata.data[k].lat,lat[1]);
            }
        }
      else
        {
          lon[0]=+1.e+35;
          lat[0]=+1.e+35;
          lon[1]=-1.e+35;
          lat[1]=-1.e+35;
          for(k=0;k<TPdata.count;k++) {
              lon[0]=MIN(TPdata.data[k].lon,lon[0]);
              lon[1]=MAX(TPdata.data[k].lon,lon[1]);
              lat[0]=MIN(TPdata.data[k].lat,lat[0]);
              lat[1]=MAX(TPdata.data[k].lat,lat[1]);
            }
        }

    }
  else
    /*-----------------------------------------------------------------------------
      GFO */
    {
      if(k%2==0) /*descending*/
        {
          lon[0]=+1.e+35;
          lat[0]=+1.e+35;
          lon[1]=-1.e+35;
          lat[1]=-1.e+35;
          for(k=0;k<TPdata.count;k++) {
              lon[0]=MIN(TPdata.data[k].lon,lon[0]);
              lon[1]=MAX(TPdata.data[k].lon,lon[1]);
              lat[0]=MIN(TPdata.data[k].lat,lat[0]);
              lat[1]=MAX(TPdata.data[k].lat,lat[1]);
            }
        }
      else
        {
          lon[0]=+1.e+35;
          lat[0]=-1.e+35;
          lon[1]=-1.e+35;
          lat[1]=+1.e+35;
          for(k=0;k<TPdata.count;k++) {
              lon[0]=MIN(TPdata.data[k].lon,lon[0]);
              lon[1]=MAX(TPdata.data[k].lon,lon[1]);
              lat[0]=MAX(TPdata.data[k].lat,lat[0]);
              lat[1]=MIN(TPdata.data[k].lat,lat[1]);
            }
        }

    }

/*   mssgrid=get_trackgrid(lon,lat); */
  printf("extrimities :  %f %f %f %f\n",lon[0],lat[0],lon[1],lat[1]);


  /*------------------------------------------------------------------------------
    smooth extremities computation */


    count=0;
    x=0;
    y=0;
    for(k=0;k<TPdata.count;k++) {
      if((fabs(lon[0]-TPdata.data[k].lon) < 0.05) || (fabs(TPdata.data[k].lat-lat[0]) < 0.05)) {
        x+=TPdata.data[k].lon;
        y+=TPdata.data[k].lat;
        count++;
        }
      }
    lon[0]=x/count;
    lat[0]=y/count;

    count=0;
    x=0;
    y=0;
    for(k=0;k<TPdata.count;k++) {
      if((fabs(lon[1]-TPdata.data[k].lon) < 0.05) || (fabs(TPdata.data[k].lat-lat[1]) < 0.05)) {
        x+=TPdata.data[k].lon;
        y+=TPdata.data[k].lat;
        count++;
        }
      }
    lon[1]=x/count;
    lat[1]=y/count;

/*
  projection=geo_mercator_init(lon[1],lat[1], (double) 0.0, (double) 90.0);
  status=geo_mercator_directe(projection,lon[0],lat[0],&x,&y);
  status=geo_mercator_directe(projection,lon[1],lat[1],&x,&y);
  x=geo_distance(lon[0],lat[0],lon[1],lat[1]);
*/

  printf("extrimities :  %f %f %f %f\n",lon[0],lat[0],lon[1],lat[1]);


}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t get_trackgridsioux(double *lon,double *lat, serie_t data)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  grid_t grid;
  grid_t cgrid;
  double dt,dp,dn;
  double tx,ty,nx,ny;
  double x,y,*d,*alpha;
  double sum1,sum2,a,c,kappa;
  grid_t backbone;
  int i,j,k,dx,dy;
  double t0,p0;


  int zone_initialised=0;

/* cartesian, rotated grid */

  cgrid.ymin=   +0.0;
  cgrid.ymax= +200.0 ;
  cgrid.xmin=   -2.0 ;
  cgrid.xmax=   +2.0 ;
  cgrid.dx  =    1.0;
  cgrid.dy  =    2.0;

  dt=lon[1]-lon[0];
  dp=lat[1]-lat[0];
  dn=sqrt(dt*dt+dp*dp);
  dt/=dn;
  dp/=dn;

  tx=dt;
  ty=dp;

  nx=-dp;
  ny= dt;

  grid.dx  =  1.0/100.;
  grid.nx  =   5;
  dx=2;

  grid.dx  =  0.25/100.;
  grid.nx  =   9;
  dx=4;

  grid.dy  =  0.5/100.;
  dy=(int)( 0.1/grid.dy );
  grid.ny  = (int)( 2*dy+dn/grid.dy );

  grid.dx  =  1.5/100.;
  grid.nx  =   3;
  dx=(int)( NINT(grid.dx-1)/2 );
  dx=1;

  grid.dx  =   1.0/100.;
  grid.nx  =   5;
  dx=2;

  grid.dy  =  1.0/100.;
  dy=(int)( 0.1/grid.dy );
  grid.ny  = (int)( 2*dy+dn/grid.dy );

/* test !!!!!!  */
/*   grid.dx  =  3.0/100.; */
/*   grid.nx  =   3; */
/*   dx=1; */
/*   grid.dy  =  0.5/100.; */
/*   dy=0.1/grid.dy; */
/*   grid.ny  = 2*dy+dn/grid.dy; */

  
  grid.modeH  = 2;

  t0=lon[0]+dx*nx*grid.dx-dy*tx*grid.dy;
  p0=lat[0]+dx*ny*grid.dx-dy*ty*grid.dy;

  grid.x=(double *)malloc(grid.nx*grid.ny*sizeof(double));
  grid.y=(double *)malloc(grid.nx*grid.ny*sizeof(double));

  grid.ymin= +1.e+35;
  grid.ymax= -1.e+35;
  grid.xmin= +1.e+35;
  grid.xmax= -1.e+35;

  for (j=0;j<grid.ny;j++)
    for (i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      grid.x[k]=t0-i*nx*grid.dx+j*tx*grid.dy;
      grid.y[k]=p0-i*ny*grid.dx+j*ty*grid.dy;
      grid.ymin=MIN(grid.ymin,grid.y[k]);
      grid.ymax=MAX(grid.ymax,grid.y[k]);
      grid.xmin=MIN(grid.xmin,grid.x[k]);
      grid.xmax=MAX(grid.xmax,grid.x[k]);
      }

  t0=lon[0]-dy*tx*grid.dy;
  p0=lat[0]-dy*ty*grid.dy;

  d=(double *)malloc(data.count*sizeof(double));
  alpha=(double *)malloc(data.count*sizeof(double));
  for (k=0;k<data.count;k++) {
    x=data.data[k].lon;
    y=data.data[k].lat;
    d[k]=(x-t0)*nx+(y-p0)*ny;
    alpha[k]=((x-t0)*tx+(y-p0)*ty)/dn;
    }

  sum1=0;
  sum2=0;

  for (k=0;k<data.count;k++) {
    c=(1.-alpha[k])*alpha[k];
    sum2+=d[k]*c;
    sum1+=c*c;
    }

  kappa=sum2/sum1;

  backbone.ny=grid.ny;
  backbone.x=(double *)malloc(backbone.ny*sizeof(double));
  backbone.y=(double *)malloc(backbone.ny*sizeof(double));
  for (j=0;j<backbone.ny;j++) {
    a=(double) j/(double)backbone.ny;
    c=(1.-a)*a;
    backbone.x[j]=t0+j*tx*grid.dy+c*kappa*nx;
    backbone.y[j]=p0+j*ty*grid.dy+c*kappa*ny;
    }

  grid.ymin= +1.e+35;
  grid.ymax= -1.e+35;
  grid.xmin= +1.e+35;
  grid.xmax= -1.e+35;

  for (j=0;j<grid.ny;j++)
    for (i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      grid.x[k]=backbone.x[j]-(i-dx)*nx*grid.dx;
      grid.y[k]=backbone.y[j]-(i-dx)*ny*grid.dx;
      grid.ymin=MIN(grid.ymin,grid.y[k]);
      grid.ymax=MAX(grid.ymax,grid.y[k]);
      grid.xmin=MIN(grid.xmin,grid.x[k]);
      grid.xmax=MAX(grid.xmax,grid.x[k]);
      }

  zone_initialised=1;
   
  if(zone_initialised!=1) {__ERR_BASE_LINE__("exiting\n");exit(-1);}

  return(grid);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int set_topographic(double *Cm, int dimM, grid_t mssgrid, float *mss_cls)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*-----------------------------------------------------------------------------
  simple test : set gaussian covariances*/
  int i,j,k,l,m;
  int status;
  float  *P2,*D,dx2;
  double x,y;
  float *vector,dum;
  float z,mask;
  float *topo,*topox,*topoy,*tmp;
  char input[1024],output[1024];
  grid_t topogrid3d;

  double norm;

/*------------------------------------------------------------------------------
  load bottom topography */

  sprintf(input,"/data/ocean/produits/topography/data/albicocca.grd");

  status=grd_loadgrid(input,&topogrid3d);
  //status=map_completegridaxis_2(&topogrid3d);
  status=map3d_completegridaxis(&topogrid3d);

  tmp=(float *)malloc(topogrid3d.nx*topogrid3d.ny*sizeof(float));
  status= grd_loadr1(input,topogrid3d,topo,&mask);

/*-----------------------------------------------------------------------------
  surface curvature constraint

  minimize a weighted (s")

  s C s = |s"| = (sP") D P"s

where:

  D  : diagonale weighting matrix
  P" : second derivative operator

  P" is a highly banded matrix

-----------------------------------------------------------------------------*/

  topox=(float *)malloc(topogrid3d.nx*topogrid3d.ny*sizeof(float));
  topoy=(float *)malloc(topogrid3d.nx*topogrid3d.ny*sizeof(float));

//  status=copy_grid_to_grid1d(topogrid3d,&topogrid3d);

  status=map_curve(topogrid3d, topo, mask,(int) 0, topox,topoy);

  sprintf(output,"../data/topo.crv");
  status=bmg_saver2(output, 1, 1, 1, topogrid3d, topox, topoy, (float) 0., mask);


/*-----------------------------------------------------------------------------
  surface seconde derivative coefficients:
  ----------------------------------------

  1st order approximation:

  s"=[s(m+1)+s(m-1)-2*s(m)]/dx

  D is a modified second derivative
  equivalent to a correlation, or an inverse L/range

  range=0.1 metres, L=0.1 D=10.


-----------------------------------------------------------------------------*/

  P2=(float *)malloc(dimM*dimM*sizeof(float));
  D=(float *)malloc(dimM*sizeof(float));

  for(k=0;k<dimM*dimM;k++) P2[k]=0;
  for(k=0;k<dimM;k++) D[k]=0;

  for(j=0;j<mssgrid.ny;j++)
    for(i=0;i<mssgrid.nx;i++) {
      k=j*mssgrid.nx+i;
      x=mssgrid.x[k];
      y=mssgrid.y[k];
      status=map_interpolation(topogrid3d, topox,mask,x,y,&z);
      if (z!=mask) D[k]=z*z;
      status=map_interpolation(topogrid3d, topoy,mask,x,y,&z);
      if (z!=mask) D[k]+=z*z;
/*-----------------------------------------------------------------------------
      appropriate if unknown is mss itself, not mss error*/
      D[k]=10*sqrt(D[k])*mss_cls[k]/fabs (topo[k]);
      D[k]=1000.;
      }

/*-----------------------------------------------------------------------------
  compute D x P2 along x direction directly in P2*/

  dx2=mssgrid.dx*mssgrid.dx;

  for(j=0;j<mssgrid.ny;j++) {
    for(i=0;i<mssgrid.nx;i++) {
      k=j*mssgrid.nx+i;
      P2[k*dimM+k]=-2./dx2/D[k];
      if(i< mssgrid.nx-1) P2[k*dimM+k+1]=1./dx2/D[k];
      if(i>0)             P2[k*dimM+k-1]=1./dx2/D[k];
      }
    }

/*-----------------------------------------------------------------------------
  compute P2 x D x P2  along x direction

   ( k, l) =   [D x P2] (k,m)  x [D x P2]  (m,l)
   ( k, l) =   [D x P2]  (m,k)  x [D x P2]  (m,l)

  it is k col of [D x P2] x l col of [D x P2] product
-----------------------------------------------------------------------------*/

  printf("upper bw of P2 %d\n",gmatrix_ubw(P2,dimM,dimM));
  printf("lower bw of P2 %d\n",gmatrix_lbw(P2,dimM,dimM));

  vector =(float *)malloc(dimM*sizeof(float));
  for(l=0;l<dimM;l++) {
    for(m=0;m<dimM;m++) vector[m]=P2[l*dimM+m];

#if ATLAS == 1
norm=cblas_dsdot(dimM, vector, 1, vector, 1);
#elif CBLAS == 1
norm=cblas_dsdot(dimM, vector, 1, vector, 1);
#else
#error "stuff missing"
#endif

    if(norm==0.) {
      continue;
      }
    for(k=MAX(0,l-4);k<=MIN(dimM-1,l+4);k++) {
#if ATLAS == 1
      dum=cblas_dsdot(dimM, &(P2[l*dimM]), 1, &(P2[k*dimM]), 1);
#elif CBLAS == 1
      dum=cblas_dsdot(dimM, &(P2[l*dimM]), 1, &(P2[k*dimM]), 1);
#else
#error "stuff missing"
#endif
      Cm[l*dimM+k]+=dum;
      }
    }
/*-----------------------------------------------------------------------------
  internal consistency chek, Cm must be symmetric */
  status=dgmatrix_symmetry(Cm, dimM);


/*-----------------------------------------------------------------------------
  compute D x P2 along y direction directly in P2*/

  for(k=0;k<dimM*dimM;k++) P2[k]=0;

  dx2=mssgrid.dy*mssgrid.dy;

  for(j=0;j<mssgrid.ny;j++) {
    for(i=0;i<mssgrid.nx;i++) {
      k=j*mssgrid.nx+i;
      P2[k*dimM+k]=-2./dx2/D[k];
      if(j< mssgrid.ny-1) P2[k*dimM+k+dimM]=1./dx2/D[k];
      if(j>0)             P2[k*dimM+k-dimM]=1./dx2/D[k];
      }
    }

/*-----------------------------------------------------------------------------
  compute P2 x D x P2  along y direction

   ( k, l) =   [D x P2] (k,m)  x [D x P2]  (m,l)
   ( k, l) =   [D x P2]  (m,k)  x [D x P2]  (m,l)

  it is k col of [D x P2] x l col of [D x P2] product
-----------------------------------------------------------------------------*/

  for(l=0;l<dimM;l++) {
    for(m=0;m<dimM;m++) vector[m]=P2[l*dimM+m];
#if ATLAS == 1
norm=cblas_dsdot(dimM, vector, 1, vector, 1);
#elif CBLAS == 1
norm=cblas_dsdot(dimM, vector, 1, vector, 1);
#else
#error "stuff missing"
#endif

    if(norm==0.) {
      continue;
      }
    for(k=MAX(0,l-4);k<=MIN(dimM-1,l+4);k++) {
#if ATLAS == 1
      dum=cblas_dsdot(dimM, &(P2[l*dimM]), 1, &(P2[k*dimM]), 1);
#elif CBLAS == 1
      dum=cblas_dsdot(dimM, &(P2[l*dimM]), 1, &(P2[k*dimM]), 1);
#else
#error "stuff missing"
#endif
      Cm[l*dimM+k]+=dum;
      }
    }

  free(vector);
  free(tmp);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int set_gaussian(double *Cm, int dimM, grid_t mssgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*-----------------------------------------------------------------------------
  simple test : set gaussian covariances*/

  double *A,*B,*C1,*C2;
  int     neq,job,*pivot,status,nrhs;
  int     i,j,k,l,m,n,kk,ll;
  double *window,lengthx,lengthy;
  double  dx,dy;
  int     hlwx,hlwy,lwx,lwy;
  char c1=CblasRowMajor;
  char c2=CblasUpper;


  lwx=11;
  lwy=21;

  hlwx=(lwx-1)/2;
  hlwy=(lwy-1)/2;

  window=(double *)malloc(lwx*lwy*sizeof(double));

/*-----------------------------------------------------------------------------
  simple test : set gaussian covariances*/
  lengthx=0.025;
  lengthy=0.025;

/*   lengthx=0.02; */
/*   lengthy=0.1; */

  for(j=0;j<lwy;j++) {
    for(i=0;i<lwx;i++) {
      k=j*lwx+i;
      dx=(i-hlwx)*mssgrid.dx/lengthx;
      dy=(j-hlwy)*mssgrid.dy/lengthy;
      window[k]=0.1*exp(-0.5*(dx*dx+dy*dy));
/*       window[k]=0.01*1./(1+(dx*dx+dy*dy)); */
      }
    }

  A=(double *)malloc(dimM*dimM*sizeof(double));
  B=(double *)malloc(dimM*dimM*sizeof(double));
  C2=(double *)malloc(dimM*dimM*sizeof(double));
  for(k=0;k<dimM*dimM;k++) A[k]=0;
  for(k=0;k<dimM*dimM;k++) C2[k]=0;

  for(j=0;j<mssgrid.ny;j++) {
    for(i=0;i<mssgrid.nx;i++) {
      k=j*mssgrid.nx+i;
      for(ll=MAX(j-hlwy,0);ll<=MIN(mssgrid.ny-1,j+hlwy);ll++) {
        for(kk=MAX(i-hlwx,0);kk<=MIN(mssgrid.nx-1,i+hlwx);kk++) {
          l=ll*mssgrid.nx+kk;
          A[l*dimM+k]=window[(ll-j+hlwy)*lwx+(kk-i+hlwx)];
          }
        }
      }
    }


/*-----------------------------------------------------------------------------
  internal consistency check, A must be symmetric */
  status=dgmatrix_symmetry(A, dimM);

  status= dgmatrix_transpose(A, B, dimM, dimM);
  C1=dmatrix_product(A, B, dimM, dimM, dimM);

  neq=dimM;

  goto skip_test;

/*-----------------------------------------------------------------------------
  internal consistency check, C1 must be symmetric and positive (+definite ?) */


#if LAPACKF_ == 1
//    status=dpotrf_(&c1, &neq,C1,&neq);
#elif ATLAS == 1
    status=clapack_dpotrf(CblasRowMajor, CblasUpper, neq,C1,neq);
#else
#error "stuff missing"
#endif
  if(status==0) printf("C1 found to be symmetric and positive; good job!!!\n");
  else
    {
    __OUT_BASE_LINE__("C1 found to be non-symmetric or non-positive (status=%d); bad boy!!!\n",status);
    exit(-1);
    }

/*-----------------------------------------------------------------------------
  C1 has been factored to perform tests, we must rebuild it*/
  free(C1);
  C1=dmatrix_product(A, B, dimM, dimM, dimM);

 skip_test:

/*-----------------------------------------------------------------------------
  internal consistency chek, C1 must be symmetric */
  status=dgmatrix_symmetry(C1, dimM);

/*-----------------------------------------------------------------------------
  inverse the covariances*/
  pivot =(int *)malloc(dimM*sizeof(int));
  job=0;

  status=poc_getrf(neq, C1, pivot);

  for(l=0;l<dimM;l++) {
     /*     for(k=0;k<dimM;k++) C2[l*dimM+k]=0.; */
    C2[l*dimM+l]=1.;
    status=poc_getrs(neq,nrhs,C1,pivot,&C2[l*dimM]);
    }

/*-----------------------------------------------------------------------------
  internal consistency chek, Cm must be symmetric */
  status=dgmatrix_symmetry(C2, dimM);

/*-----------------------------------------------------------------------------
  copy C2 in Cm, in such a way that Cm is built to be symmetric */
  for (l=0;l<dimM;l++)
    for (k=l;k<dimM;k++) {
      m=dimM*k+l;
      n=dimM*l+k;
      Cm[n]=0.5*(C2[m]+C2[n]);
      Cm[m]=Cm[n];
      }

  free(C1);
  free(A);
  free(B);
  free(C2);
  free(pivot);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int set_byderivative(double *Cm, int dimM, grid_t mssgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*-----------------------------------------------------------------------------
  s*/

  double *A,*B,*C1,*C2;
  int     neq,job,*pivot,status;
  int     i,j,k,l,m,n,kk,ll;
  double *window,lengthx,lengthy;
  double  dx,dy,x,y,variance;
  double  dx2,dy2;

/*-----------------------------------------------------------------------------
  correlation distance*/
  lengthx=0.1;
  lengthy=0.05;


/*-----------------------------------------------------------------------------
  normalize*/
  lengthx/=mssgrid.dx;
  lengthy/=mssgrid.dy;

  dx2=lengthx*lengthx;
  dy2=lengthy*lengthy;

  variance=0.1;

  A=(double *)malloc(dimM*dimM*sizeof(double));
  B=(double *)malloc(dimM*dimM*sizeof(double));
  
  for(k=0;k<dimM*dimM;k++) A[k]=0;

  for(j=0;j<mssgrid.ny;j++) {
    for(i=0;i<mssgrid.nx;i++) {
      k=j*mssgrid.nx+i;
      x=mssgrid.x[k];
      y=mssgrid.x[k];

      if((i!=0) && (i!=mssgrid.nx-1)) {
        l=j*mssgrid.nx+i-1;
        A[l*dimM+k]=dx2;
        l=j*mssgrid.nx+i+1;
        A[l*dimM+k]=dx2;
        A[k*dimM+k]-=2.*dx2;
        }

/*-----------------------------------------------------------------------------
      shifted condition*/
/*
      if(i==0) {
        l=j*mssgrid.nx+i;
        A[l*dimM+k]=dx2*100;
        l=j*mssgrid.nx+i+2;
        A[l*dimM+k]=dx2*100;
        l=j*mssgrid.nx+i+1;
        A[k*dimM+k]-=2*dx2*100;
        }
      if(i==mssgrid.nx-1) {
        l=j*mssgrid.nx+i-2;
        A[l*dimM+k]=dx2*100;
        l=j*mssgrid.nx+i;
        A[l*dimM+k]=dx2*100;
        l=j*mssgrid.nx+i-1;
        A[k*dimM+k]-=2*dx2*100;
        }

*/
/*-----------------------------------------------------------------------------
      continuty condition*/
/*

      if(i==0) {
        l=j*mssgrid.nx+i+1;
        A[l*dimM+k]=dx2;
        A[k*dimM+k]-=2*dx2;
        }
      if(i==mssgrid.nx-1) {
        l=j*mssgrid.nx+i-1;
        A[l*dimM+k]=dx2;
        A[k*dimM+k]-=2*dx2;
        }
*/
/*-----------------------------------------------------------------------------
  slope condition*/

      if(i==0) {
        l=j*mssgrid.nx+i;
        A[l*dimM+k]=-1./mssgrid.dx;
        l=j*mssgrid.nx+i+1;
        A[l*dimM+k]=+1./mssgrid.dx;
        }

      if(i==mssgrid.nx-1) {
        l=j*mssgrid.nx+i-1;
        A[l*dimM+k]=+1./mssgrid.dx;
        l=j*mssgrid.nx+i;
        A[l*dimM+k]=-1./mssgrid.dx;
        }


      }
    }

  status= dgmatrix_transpose(A, B, dimM, dimM);
  C1=dmatrix_product02(A, B, dimM, dimM, dimM);

  for(k=0;k<dimM*dimM;k++) A[k]=0;

  for(j=0;j<mssgrid.ny;j++) {
    for(i=0;i<mssgrid.nx;i++) {
      k=j*mssgrid.nx+i;
      x=mssgrid.x[k];
      y=mssgrid.x[k];

      if((j!=0) && (j!=mssgrid.ny-1)) {
        l=(j-1)*mssgrid.nx+i;
        A[l*dimM+k]=dy2;
        l=(j+1)*mssgrid.nx+i;
        A[l*dimM+k]=dy2;
        A[k*dimM+k]-=2*dy2;
        }

      }
    }

  status= dgmatrix_transpose(A, B, dimM, dimM);
  C2=dmatrix_product02(A, B, dimM, dimM, dimM);

/*-----------------------------------------------------------------------------
  internal consistency check, A must be symmetric */
  status=dgmatrix_symmetry(C1, dimM);
  status=dgmatrix_symmetry(C2, dimM);

/*-----------------------------------------------------------------------------
  copy A in Cm, in such a way that Cm is built to be symmetric */

  for(k=0;k<dimM*dimM;k++) Cm[k]=0;

  for (l=0;l<dimM;l++) {
    m=dimM*l+l;
    Cm[m]=1./variance;
    }

  for (l=0;l<dimM;l++)
    for (k=l;k<dimM;k++) {
    m=dimM*l+k;
    Cm[m]+=C1[m]+C2[m];
    }

  for (l=0;l<dimM;l++)
    for (k=l;k<dimM;k++) {
    m=dimM*k+l;
    n=dimM*l+k;
    Cm[m]=Cm[n];
    }

  free(A);
  free(B);
  free(C1);
  free(C2);

  return(0);
}

static bool expired=expire(20130922,20131122);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#define rsize 1000

  float  *buf,*error_cls,dum,z,r,hh;
  geo_t projection;
  float  meanh,meanp;
  double time,*tobs,*tdum;
  float ftime;
  double weights[4],latmin,latmax,lonmean,latmean;
  float  georef[100],ignref[100];
  int    nodes[4];

  float  mask=1.e+10;
  double  x,y;
  double *lon,*lat;
  double avg[2][nstat],var[2][nstat];
  double  *A,*b;

  int i,j,k,l,L,m,n;
  int nvalue,nframe,nobs,nsample,count;
  int neq,job,*pivot,status;
  int dimM,dimD;
  int ldaM,ldaD;
  int hlw,lw;
  int diagM,diagD;
  int TPpoint;

  size_t size;
  FILE *in,*out;
  char *keyword,*zone,*s;

  char *sat=NULL;
  char *rootname=NULL,*output=NULL,*outroot=NULL,*input=NULL,*mssfile=NULL;
  char *xover=NULL,*gauge=NULL,*poly=NULL;
  char sample_path[1024]="./";

  covariance_t c;
  regression_t r1,r2;
  statistic_t ss;

  grid_t grid,mssgrid;

  serie_t TPdata,TGdata;

  int    trk[rsize],type[rsize];
  float  *mss_cls,*model,*analysis,*innovation,*observation;
  double *error;
  float  *tmp;
  double *Cm;
  float  *Cd,*G,*tG;
  float  *tmp1,*tmp2,*tmp3;
  float  *svector,snorm;
  double *dvector,dnorm;
  float  z1,z2,z3,z4;
  float  *oceanic1,*oceanic2,*oceanic3,*oceanic4;
  float  *hobs,*mog,*tidesTG,*tides;
  float  *delta1,*delta2;
  float  *mss1,*mss2,*mss3,*mss4,*hf;
  float  *rawTG,*raw;
  float  *variance1,*variance2;
  int nindex,*index;
  
  fct_echo(argc,argv);
 
  plg_t *polygones=NULL;
  int npolygones=0;

  n=1;
  while (n < argc)
    {
    keyword=strdup(argv[n]);
    switch (keyword[0])
      {
      case '-':
        switch (keyword[1])
        {
        case 'r' :
          rootname= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*         case 's' : */
/*           sprintf(sample_path,"%s",argv[n+1]); */
/*           n++; */
/*           n++; */
/*           break; */

        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          sat= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'x' :
          xover= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'g' :
          gauge= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          poly= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          outroot= strdup(argv[n+1]);
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
        else
          {
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
          }
        break;
      }
      free(keyword);
    }

  input=(char *)malloc(1024);
  output=(char *)malloc(1024);

  if(outroot==NULL) outroot=strdup("../data");

  sprintf(input,"../data/polygones/%s.plg",poly);
  if(poly!=NULL) plg_load_scan(input, &polygones, &npolygones);

/*-----------------------------------------------------------------------------
  Read altimeter data */

  sprintf(input,"../data/track-raw/track-raw.%s.%s.dat",sat,xover);
/*   sprintf(input,"/calcul/maewo/roblou/processingTP2/track-raw.%s.%s.dat",sat,xover);  */
  TPdata=load_metadata_raw(input, polygones, npolygones,poly);

  if(TPdata.count==0) {__ERR_BASE_LINE__("exiting\n");exit(-1);}

/*------------------------------------------------------------------------------
  Read tide gauge data */

  sprintf(input,"../data/gauges/%s.TG.dat",gauge);
  TGdata=load_metadata_ref(input,1);

  if(TGdata.count==0) {__ERR_BASE_LINE__("exiting\n");exit(-1);}

/*-----------------------------------------------------------------------------
  Select ascending or descending, depends on satellite:
  a GFO ascending track is different than a T/P ascending one!! */
  lon=(double *) calloc(2,sizeof(double));
  lat=(double *) calloc(2,sizeof(double));

  get_extremities(sat,xover,TPdata,lon,lat);

/*------------------------------------------------------------------------------
  Create mss grid */
  mssgrid=get_trackgridsioux(lon,lat,TPdata);

  sprintf(output,"%s/track.%s.%s.%s.grid.nc",outroot,sat,xover,poly);
  status=mss_createfile(output,(size_t) mssgrid.nx,(size_t) mssgrid.ny, mssgrid);

/*------------------------------------------------------------------------------
  Load CLS mss */
  mssfile=(char *)malloc(1024);
  sprintf(mssfile,"../data/%s%s%s","mss.",zone,".bmg");

  status=bmg_loadgrid (mssfile,&grid);

  status=map_completegridaxis_2(&grid);

  buf=(float *)malloc(grid.nx*grid.ny*sizeof(float));
  if (buf == NULL) {__ERR_BASE_LINE__("exiting\n");exit(8);}
  error_cls=(float *)malloc(grid.nx*grid.ny*sizeof(float));
  if (error_cls == NULL)exit(11) ;

  status= bmg_loadr1 (mssfile,1,1,1,grid,buf,&mask,&ftime);
  if (status != 0) exit(9) ;

  status= bmg_loadr1 (mssfile,2,1,1,grid,error_cls,&mask,&ftime);
  if (status != 0) exit(10) ;


/*------------------------------------------------------------------------------
  Interpolate and store mss on track grid*/

  mss_cls=(float *)malloc( mssgrid.nx*mssgrid.ny*sizeof(float));
  for(j=0;j<mssgrid.ny;j++)
    for(i=0;i<mssgrid.nx;i++) {
      x=mssgrid.x[j*mssgrid.nx+i];
      y=mssgrid.y[j*mssgrid.nx+i];
      status=map_interpolation(grid, buf,mask,x,y,&z);
      if (z!=mask) mss_cls[j*mssgrid.nx+i]=z/1000.;
      else mss_cls[j*mssgrid.nx+i]=1.0e+35;
      }

  sprintf(output,"%s/track.%s.%s.%s.mss.nc",outroot,sat,xover,poly);
  status=mss_createfile(output,(size_t) mssgrid.nx,(size_t) mssgrid.ny, mssgrid);


/*-----------------------------------------------------------------------------

Inverse problem solution
------------------------

For details and foramlism, see Tarantola (Elsevier), pp197

Penalty function J:

J=(d-Gm)Cd(d-Gm)+(m-m0)Cm(m-m0)

Minimizing J (analysis step):

m(optimal)=m'=m0+[GxCdxG+Cm]xGxCd(d-Gm)

Posterior covariance of inverse analysis:

Cm'=[GxCdxG+Cm]

where:

 denotes inverse
 denotes transpose

Cd  = prior error covariance on data
Cm  = prior error covariance on model
Cm' = posterior error covariance on analysis
G   = interpolation/observation operator

Cd=inverse of Cd
Cm=inverse of Cm
G = transpose of G

Cd is dimDxdimD square matrix (huge in this application)
Cm is dimMxdimM square matrix (reasonable in this application)

GxCd is dimMxdimG rectangular matrix

I. Quick analysis problem
xxxxxxxxxxxxxxxxxxxxxxxx

3 steps technique:

1) form GxCd and [GxCdxG+Cm]

2) solve [GxCdxG+Cm]X=GxCd(d-Gm)

3) m'=m0+X

posterior covariances not available !!!
---------------------------------------

II. Full analysis problem
xxxxxxxxxxxxxxxxxxxxxxxx

3 steps technique:

1) form GxCd and [GxCdxG+Cm]

2) explicitely compute Cm'=[GxCdxG+Cm]

3) m'=m0+Cm'x GxCd(d-Gm)

posterior covariances available !!!
-----------------------------------

III. Consistency check
xxxxxxxxxxxxxxxxxxxxxx

Jminimum=(d-Gm')Cd(d-Gm')+(m'-m0)Cm(m'-m0) ~ Nobs

-----------------------------------------------------------------------------*/
/*   TPdata.count=1; */
/*-----------------------------------------------------------------------------
  matrix leading dimensions */
  dimM=mssgrid.nx*mssgrid.ny;
  dimD=TPdata.count;

/*-----------------------------------------------------------------------------
  square inverse covariance matrix initialisation.
  Cm and Cd are actually the inverse matrices */
  Cm=(double *)malloc(dimM*dimM*sizeof(double));
  for(k=0;k<dimM*dimM;k++) Cm[k]=0.;

/*-----------------------------------------------------------------------------
  Cd is (horizontal) band matrix; diagonal at row 1+ldaD/2 */
  ldaD=1;
  diagD=ldaD/2;
  Cd=(float *)malloc(dimD*ldaD*sizeof(float));
  for(k=0;k<dimD*ldaD;k++) Cd[k]=0.;

/*-----------------------------------------------------------------------------
  diagonal (inverse variance) constraint */
  for(k=0;k<dimM;k++) Cm[k*dimM+k]=1./0.01; /* 100 cm =0.01 m */
  for(k=0;k<dimD;k++) Cd[k*ldaD+diagD]=1./0.002; /* 20 cm =0.002 m */

/*-----------------------------------------------------------------------------
  observation matrix initialisation */
  if( (G =(float *)malloc(dimM*dimD*sizeof(float))) == NULL) gmerror("can not allocate G matrix");
  if( (tG=(float *)malloc(dimM*dimD*sizeof(float))) == NULL) gmerror("can not allocate tG matrix");
  for(k=0;k<dimM*dimD;k++) G[k]=0;

/*-----------------------------------------------------------------------------
  get data point interpolation coeffcients (observation operator)*/
  for(k=0;k<dimD;k++) {
    x=TPdata.data[k].lon;
    y=TPdata.data[k].lat;
    status=map_coeffcients02(mssgrid,x,y,weights,nodes,&count);
    if(status==-1) {
      printf("out of grid: %5d %5.3f %6.3f\n",k,x,y);
      }
    for(m=0;m<count;m++) G[k+nodes[m]*dimD]=weights[m];
    }

  for(k=0;k<dimD;k++) {
    dum=0.;
    for(m=0;m<dimM;m++) dum+=G[m*dimD+k]*mss_cls[m];
    if(fabs(dum-TPdata.data[k].values[24]) > 0.01) {
      printf("%6d %lf %6.3f %6.3f %6.3f %6.3f %3d \n",k,TPdata.data[k].time,
                                                            TPdata.data[k].lon,TPdata.data[k].lat,
                                                            TPdata.data[k].values[24],dum,
                                                            TPdata.data[k].cycle);
      }
    }

/*-----------------------------------------------------------------------------
  compute G transpose, i.e. G */

  status=gmatrix_transpose(G, tG, dimD, dimM);


/*-----------------------------------------------------------------------------
  simple test : set diagonal covariance matrix Cm
  goto next;
*/


/*-----------------------------------------------------------------------------
  simple test : set diagonal covariance matrix Cm*/

  status=set_byderivative(Cm,  dimM,  mssgrid);
  goto next;

  status=set_gaussian(Cm,  dimM,  mssgrid);
  goto next;

  status=set_topographic(Cm,  dimM,  mssgrid, mss_cls);
  goto next;

 next:

  printf("model covariance matrix done...\n");

/*-----------------------------------------------------------------------------
  internal consistency chek, Cm must be symmetric */
  status=dgmatrix_symmetry(Cm, dimM);

  printf("rank of Cm     = %d\n",dimM);
  printf("upper bw of Cm = %d\n",dgmatrix_ubw(Cm,dimM,dimM));
  printf("lower bw of Cm = %d\n",dgmatrix_lbw(Cm,dimM,dimM));

  variance1 =(float *)malloc(dimM*sizeof(float));
  for(k=0;k<dimM;k++) variance1[k]=Cm[k*dimM+k];

/*-----------------------------------------------------------------------------
  compute G x Cd */

  tmp1 =(float *)malloc(dimM*dimD*sizeof(float));
  for(k=0;k<dimM*dimD;k++) tmp1[k]=0;

  for(l=0;l<dimD;l++) {
    for(k=0;k<dimM;k++) {
      for(m=l-diagD;m<=l+diagD;m++) {
        n=l-m+diagD;                    /* m,l */
/*      tmp1[l*dimM+k]+=tG[m*dimM+k]*Cd[l*dimD+m]; obsolete */
        tmp1[l*dimM+k]+=tG[m*dimM+k]*Cd[l*ldaD+n];
        }
      }
    }
  status=gmatrix_ubw(tG,dimM,dimD);
  status=gmatrix_ubw(tmp1,dimM,dimD);

/*-----------------------------------------------------------------------------
  compute G x Cd x G */

  tmp2 =(float *)malloc(dimM*dimM*sizeof(float));
  for(k=0;k<dimM*dimM;k++) tmp2[k]=0;

  svector =(float *)malloc(dimD*sizeof(float));
  index  =(int *)malloc(dimD*sizeof(int));

  for(k=0;k<dimM;k++) {
/*-----------------------------------------------------------------------------
    get row k of tmp1 */
    nindex=0;
    for(m=0;m<dimD;m++) {
      svector[m]=tmp1[m*dimM+k];
      if(svector[m] !=0. ) {
        nindex++;
        index[nindex-1]=m;
        }
      }
    if(nindex==0) {
      continue;
      }

    if(nindex < 100) {
      for(l=0;l<dimM;l++) {
/*-----------------------------------------------------------------------------
      scalar product with col l of G */
        for (n=0;n<nindex;n++) {
          m=index[n];
          tmp2[l*dimM+k]+= svector[m]*G[l*dimD+m];
          }
        }
      }
    else {
      for(l=0;l<dimM;l++)
#if ATLAS == 1
   tmp2[l*dimM+k]=cblas_dsdot(dimD, svector, 1,&(G[l*dimD]) , 1);
#elif CBLAS == 1
   tmp2[l*dimM+k]=cblas_dsdot(dimD, svector, 1,&(G[l*dimD]) , 1);
#else
  #error "stuff missing"
#endif
      }
    }

/*-----------------------------------------------------------------------------
  old method
  for(k=0;k<dimM;k++) {
    for(m=0;m<dimD;m++) vector[m]=tmp1[m*dimM+k];
    norm= sdot( dimD, vector, 1, vector, 1);
    if(norm==0.) {
      continue;
      }
    for(l=0;l<dimM;l++) {
      tmp2[l*dimM+k]= sdot( dimD, vector, 1, &(G[l*dimD]), 1);
      }
    }
*/
  free(svector);
  free(index);

/*-----------------------------------------------------------------------------
  internal consistency chek, tmp2 must be symmetric */
  status=sgmatrix_symmetry(tmp2, dimM);

  printf("upper bw of  G x Cd x G %d\n",gmatrix_ubw(tmp2,dimM,dimM));
  printf("lower bw of  G x Cd x G %d\n",gmatrix_lbw(tmp2,dimM,dimM));

/*-----------------------------------------------------------------------------
  compute G x Cd x G + Cm */

  A =(double *)malloc(dimM*dimM*sizeof(double));
  for(k=0;k<dimM*dimM;k++) A[k]=tmp2[k]+Cm[k];
/*   for(k=0;k<dimM*dimM;k++) A[k]=tmp2[k]; */
/*   for(k=0;k<dimM*dimM;k++) A[k]=Cm[k]; */

/*-----------------------------------------------------------------------------
  internal consistency chek, tmp2 must be symmetric */
  status=dgmatrix_symmetry(A, dimM);

/*-----------------------------------------------------------------------------
  compute d-Gm  */

  model=mss_cls;

  observation=(float *)malloc(dimD*sizeof(float));
  for(k=0;k<dimD;k++) {
/*     observation[k]=TPdata.data[k].values[0]-TPdata.data[k].values[1]-TPdata.data[k].values[12]; LR, 07/02/05, zone AMSUD!!! */
    observation[k]=TPdata.data[k].values[0]-TPdata.data[k].values[1]-TPdata.data[k].values[18];
    observation[k]=TPdata.data[k].values[29];
    }
 
  error =(double *)malloc(dimD*sizeof(double));
  for(k=0;k<dimD;k++) {
    error[k]=observation[k];
    dum=0.;
/*-----------------------------------------------------------------------------
    compute row k, col m  of G x model  */
    for(m=0;m<dimM;m++) dum+=G[m*dimD+k]*model[m];
    error[k]-=dum;
    }

#if ATLAS == 1
   dnorm= ( cblas_ddot(dimD, error, 1,error, 1) / dimD );
#elif CBLAS == 1
   dnorm= ( cblas_ddot(dimD, error, 1,error, 1) / dimD );
#else
#error "stuff missing"
#endif

  printf("prior error euclidian norm (rms) : %f \n",sqrt(dnorm));

/*-----------------------------------------------------------------------------
  compute G x Cd x (d-Gm)  */

  b =(double *)malloc(dimM*sizeof(double));
  for(k=0;k<dimM;k++) {
    b[k]=0.;
/*-----------------------------------------------------------------------------
    compute row k, col m  of (G x Cd) x row m of error */
    for(m=0;m<dimD;m++) b[k]+=tmp1[m*dimM+k]*error[m];
    }

/*-----------------------------------------------------------------------------
  compute prior J */

  dvector =(double *)malloc(dimD*sizeof(double));
  for(k=0;k<dimD;k++) {
    dvector[k]=0;
    for(m=k-diagD;m<=k+diagD;m++) {
      n=k-m+diagD;
      dvector[k]+=Cd[k*ldaD+n]*error[m];
      }
    }

#if ATLAS == 1
 dnorm=cblas_ddot(dimD, dvector, 1,error, 1);
#elif CBLAS == 1
 dnorm=cblas_ddot(dimD, dvector, 1,error, 1);
#else
#error "stuff missing"
#endif
  free(dvector);

  printf("prior error penalty function : %f %f \n",dnorm,dnorm/dimD);

/*-----------------------------------------------------------------------------
  convert linear problem matrix into band matrix shape*/

/*-----------------------------------------------------------------------------
  solve the linear problem using ??? later should be Cholesky*/

  job=0;

/*   printf("upper bw of A %d\n",gmatrix_ubw(A,dimM,dimM)); */
/*   printf("lower bw of A %d\n",gmatrix_lbw(A,dimM,dimM)); */

  neq=dimM;
  pivot =(int *)malloc(dimM*sizeof(int));

  int nrhs=1;
  status=poc_getrs(neq,nrhs,A,pivot,b);
  free(A);
  free(pivot);


/*-----------------------------------------------------------------------------
  compute analysis  */

  analysis   =(float *)malloc(dimM*sizeof(float));
  innovation =(float *)malloc(dimM*sizeof(float));
  for(k=0;k<dimM;k++) {
    innovation[k]=b[k];
    analysis[k]=model[k]+b[k];
    }

/*-----------------------------------------------------------------------------
  compute posterior d-Gm  */

  model=analysis;

  for(k=0;k<dimD;k++) {
    error[k]=observation[k];
    dum=0.;
/*-----------------------------------------------------------------------------
    compute row k, col m  of G x model  */
    for(m=0;m<dimM;m++) dum+=G[m*dimD+k]*model[m];
    error[k]-=dum;
    }

#if ATLAS == 1
 dnorm=cblas_ddot(dimD, error, 1,error, 1) / dimD;
#elif CBLAS == 1
 dnorm=cblas_ddot(dimD, error, 1,error, 1) / dimD;
#else
#error "stuff missing"
#endif
  printf("posterior error euclidian norm (rms) : %f \n",sqrt(dnorm));

/*-----------------------------------------------------------------------------
  compute posterior J */

/*-----------------------------------------------------------------------------
  data departure = (d-Gm') Cd (d-Gm') */
  dvector =(double *)malloc(dimD*sizeof(double));
  for(k=0;k<dimD;k++) {
    dvector[k]=0;
    for(m=k-diagD;m<=k+diagD;m++) {
      n=k-m+diagD;
      dvector[k]+=Cd[k*ldaD+n]*error[m];
      }
    }

#if ATLAS == 1
 dnorm=cblas_ddot(dimD,dvector , 1,error, 1);
#elif CBLAS == 1
 dnorm=cblas_ddot(dimD,dvector , 1,error, 1);
#else
#error "stuff missing"
#endif
  free(dvector);
  printf("posterior error penalty function : %f %f \n",dnorm,dnorm/dimD);

/*-----------------------------------------------------------------------------
  model departure = (m0-m') Cm (m0-m') */
  dvector =(double *)malloc(dimM*sizeof(double));
  for(k=0;k<dimM;k++) {
    dvector[k]=0;
    for(m=0;m<dimM;m++) {
      dvector[k]+=Cm[m*dimM+k]*b[m];
      }
    }

#if ATLAS == 1
 dnorm+=cblas_ddot(dimD,dvector , 1,b, 1);
#elif CBLAS == 1
 dnorm+=cblas_ddot(dimD,dvector , 1,b, 1);
#else
#error "stuff missing"
#endif
  free(dvector);

  printf("posterior error penalty function : %f %f \n",dnorm,dnorm/dimD);

/*-----------------------------------------------------------------------------
  save optimal mss*/

  for(k=0;k<dimM;k++) {
    innovation[k]*=100; /* innovtion in cm */
    }

  sprintf(output,"%s/track.%s.%s.%s.mss.nc",outroot,sat,xover,poly);
  status=mss_createfile(output,(long unsigned int) mssgrid.nx,(long unsigned int) mssgrid.ny, mssgrid);
  status=writefile(output, analysis,mss_cls,innovation);

  oceanic1=(float *)malloc(dimD*sizeof(float));
  oceanic2=(float *)malloc(dimD*sizeof(float));
  oceanic3=(float *)malloc(dimD*sizeof(float));
  oceanic4=(float *)malloc(dimD*sizeof(float));

  raw=(float *)malloc(dimD*sizeof(float));

  delta1=(float *)malloc(dimD*sizeof(float));
  delta2=(float *)malloc(dimD*sizeof(float));

  mss1=(float *)malloc(dimD*sizeof(float));
  mss2=(float *)malloc(dimD*sizeof(float));
  mss3=(float *)malloc(dimD*sizeof(float));
  mss4=(float *)malloc(dimD*sizeof(float));
  hf=(float *)malloc(dimD*sizeof(float));
  tdum=(double *)malloc(dimD*sizeof(double));
  tides=(float *)malloc(dimD*sizeof(float));

  hobs=(float *)malloc(TGdata.count*sizeof(float));
  mog=(float *)malloc(TGdata.count*sizeof(float));
  tobs=(double *)malloc(TGdata.count*sizeof(double));
  tidesTG=(float *)malloc(TGdata.count*sizeof(float));
  rawTG=(float *)malloc(TGdata.count*sizeof(float));

  for(k=0;k<TGdata.count;k++) {
    hobs[k]=TGdata.data[k].values[29];
    mog[k]=TGdata.data[k].values[1];
    tidesTG[k]=TGdata.data[k].values[21];
    rawTG[k]=TGdata.data[k].values[28];
    tobs[k]=TGdata.data[k].time;
    }

  mask=9999.9;
  for(k=0;k<dimD;k++) {
    tdum[k]=TPdata.data[k].time;
    dum=0.;
    for(m=0;m<dimM;m++) dum+=G[m*dimD+k]*analysis[m];
    mss1[k]=dum;
    dum=0.;
    for(m=0;m<dimM;m++) dum+=G[m*dimD+k]*mss_cls[m];
    mss2[k]=dum;
    mss3[k]=TPdata.data[k].values[24];
    x=TPdata.data[k].lon;
    y=TPdata.data[k].lat;
    status=map_interpolation(grid, buf,mask,x,y,&dum);
    mss4[k]=dum/1000.0;
    hf[k]=TPdata.data[k].values[1]+TPdata.data[k].values[10]+TPdata.data[k].values[11]+TPdata.data[k].values[12];
    status=map_interpolate1D(mog, tobs, TGdata.mask, TGdata.count, TPdata.data[k].time, &(hf[k]));
    status=map_interpolate1D(tidesTG, tobs, TGdata.mask, TGdata.count, TPdata.data[k].time, &(tides[k]));
    oceanic1[k]=TPdata.data[k].values[29]-mss1[k];
    oceanic2[k]=TPdata.data[k].values[29]-mss2[k];
    oceanic3[k]=TPdata.data[k].values[29]-mss3[k];
    status=map_interpolate1D(hobs, tobs, TGdata.mask, TGdata.count, TPdata.data[k].time, &(oceanic4[k]));
    status=map_interpolate1D(rawTG, tobs, TGdata.mask, TGdata.count, TPdata.data[k].time, &(raw[k]));
/*
    if(status<-1) {
        printf("%d %lf %f %d\n",k,TPdata.data[k].time,oceanic4[k],status);
      }
*/
    status=map_interpolation(grid, error_cls,mask,x,y,&dum);
    error[k]=dum;
    if(fabs(mss4[k]-mss3[k]) > 0.01) {
      printf("%6d %lf %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %4d \n",
                                                            k,TPdata.data[k].time,TPdata.data[k].lon,TPdata.data[k].lat,
                                                            mss1[k],mss2[k],mss3[k],mss4[k],
                                                            TPdata.data[k].cycle);
      }
    }
  free(tobs);


  /*-----------------------------------------------------------------------------
    create time series on mean tracks */
/*   save_meantracks_CLS(TPdata,xover,mss1,polygones,npolygones);  */
  save_meantracks_CTOH(TPdata,sat,xover,mss1,polygones,npolygones,outroot);


  printf("TP (mss=LEGOS) \n");
  ss=get_statistics(oceanic1, mask, dimD, 1);
  printf("TP (mss=CLS) \n");
  ss=get_statistics(oceanic2, mask, dimD, 1);
  printf("TP (mss=CLS) \n");
  ss=get_statistics(oceanic3, mask, dimD, 1);
  printf("TG \n");
  ss=get_statistics(oceanic4, TGdata.mask, dimD, 1);
 
/*-----------------------------------------------------------------------------

  1,2  : time, latitude
  3..6 : signal altim�rique dealias�- mss LEGOS, -mss CLS01 , -mss GDR, observation dealias�  7    : corrections
  8..10: mss CLS-01 - LEGOS, CLS01 -GDR, CLS01 -GDR, CLS01 (grille LEGOS - interpolation directe)
  9    : erreur associ� a CLS01
  10   : cycle
-----------------------------------------------------------------------------*/

  sprintf(output,"%s/track.%s.%s.%s.mss.dat",outroot,sat,xover,poly);
  out=fopen(output,"w");
  for(k=0;k<dimD;k++) {
    fprintf(out,"%lf %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %d \n",
                                                            TPdata.data[k].time,TPdata.data[k].lat,
                                                            oceanic1[k],oceanic2[k],oceanic3[k],oceanic4[k],
                                                            TPdata.data[k].values[0]-TPdata.data[k].values[29],
                                                            mss2[k]-mss1[k],mss2[k]-mss3[k],mss4[k]-mss3[k],
                                                            error[k],TPdata.data[k].cycle);
    if(oceanic1[k]< -1.0) {
        continue;
      }
    }
  fclose(out);

/*-----------------------------------------------------------------------------
  
  1,2  : time, latitude
  3..5 : signal altim�rique dealias�- mss LEGOS, -mss CLS01 , observation dealias�  6,7  : signal altim�rique dealias�-observation - mss LEGOS, -mss CLS01
  8    : ???
  9,10 : MOG global, MOG global + LSA + solid + FES
  11,12: SLA, tides
  13   : correction
  14   : cycle
 -----------------------------------------------------------------------------*/

#define MOG 12

  sprintf(output,"%s/track.TP.%s.%s.mss.dat",outroot,xover,poly);
  sprintf(output,"%s/track.TP.%s.%s.%s.cal.dat",outroot,xover,poly,gauge);
  out=fopen(output,"w");
  for(k=0;k<dimD;k++) {
    if(oceanic4[k]!=TGdata.mask)
    fprintf(out,"%lf %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %d \n",
                                                            TPdata.data[k].time,TPdata.data[k].lat,
                                                            oceanic1[k],oceanic2[k],oceanic4[k],
                                                            oceanic1[k]-oceanic4[k],oceanic2[k]-oceanic4[k],
                                                            oceanic1[k]-oceanic4[k]+TPdata.data[k].values[MOG]-tides[k],
                                                            TPdata.data[k].values[1],hf[k],
                                                            TPdata.data[k].values[12],tides[k],
                                                            TPdata.data[k].values[0]-TPdata.data[k].values[29],
                                                            TPdata.data[k].cycle);
    }
  fclose(out);

  for(k=0;k<dimD;k++) {
    if(oceanic4[k]!=TGdata.mask)  oceanic3[k]=TPdata.data[k].values[28]-mss1[k];
    else oceanic3[k]=TGdata.mask;
    }

  status=compute_annual(oceanic3,tdum,TGdata.mask,dimD);

  for(k=0;k<dimD;k++) {
    if(oceanic4[k]!=TGdata.mask)  oceanic3[k]=TPdata.data[k].values[28]-mss2[k];
    else oceanic3[k]=TGdata.mask;
    }

  status=compute_annual(oceanic3,tdum,TGdata.mask,dimD);

  for(k=0;k<dimD;k++) {
    if(oceanic4[k]!=TGdata.mask)  oceanic3[k]=raw[k];
    else oceanic3[k]=TGdata.mask;
    }

  status=compute_annual(oceanic3,tdum,TGdata.mask,dimD);

  for(latmin=mssgrid.ymin-0.05;latmin<mssgrid.ymax;latmin+=.05) {
    latmax=latmin+.1;
    lonmean=0.;
    latmean=0.;
    count=0;
    for(k=0;k<dimD;k++) {
      if(oceanic4[k]!=TGdata.mask)  oceanic3[k]=oceanic1[k]-oceanic4[k];
      else oceanic3[k]=mask;
      }
    for(k=0;k<dimD;k++) {
      if((TPdata.data[k].lat>latmax) || (TPdata.data[k].lat<latmin)) {
        oceanic3[k]=mask;
        }
      }
    for(k=0;k<dimD;k++) {
      if(oceanic3[k]!=mask) {
        lonmean+=TPdata.data[k].lon;
        latmean+=TPdata.data[k].lat;
        count++;
        }
      }
    if(count==0) continue;
    lonmean/=count;
    latmean/=count;
    printf("%9.3f %9.3f %9.3f\n",latmin,latmax,geo_km(TGdata.data[0].lon,TGdata.data[0].lat,lonmean,latmean));
    for(k=0;k<dimD;k++) delta1[k]=oceanic3[k];
    printf("TP - TG (mss=LEGOS) \n");
    ss=sget_filteredstatistics(oceanic3, mask, dimD, 2.);
/*
    for(k=0;k<dimD;k++) {
      if(oceanic4[k]!=TGdata.mask)  oceanic3[k]=TPdata.data[k].values[28]-mss1[k]-raw[k];
      else oceanic3[k]=mask;
      }
    for(k=0;k<dimD;k++) {
      if((TPdata.data[k].lat>latmax) || (TPdata.data[k].lat<latmin)) {
        oceanic3[k]=mask;
        }
      }
    ss=sget_filteredstatistics(oceanic3, mask, dimD, 2.);
*/
    for(k=0;k<dimD;k++) {
      if(oceanic4[k]!=TGdata.mask)  oceanic3[k]=oceanic2[k]-oceanic4[k];
      else oceanic3[k]=mask;
      }
    for(k=0;k<dimD;k++) {
      if((TPdata.data[k].lat>latmax) || (TPdata.data[k].lat<latmin)) {
        oceanic3[k]=mask;
        }
      }
/*
    ss=sget_filteredstatistics(oceanic3, mask, dimD, 2.);
*/

    for(k=0;k<dimD;k++) {
      if(oceanic4[k]!=TGdata.mask)  oceanic3[k]=oceanic1[k];
      else oceanic3[k]=mask;
      }
    for(k=0;k<dimD;k++) {
      if((TPdata.data[k].lat>latmax) || (TPdata.data[k].lat<latmin)) {
        oceanic3[k]=mask;
        }
      }
    for(k=0;k<dimD;k++) delta2[k]=oceanic3[k];
    printf("TP (mss=LEGOS) \n");
    ss=sget_filteredstatistics(oceanic3, mask, dimD, 2.);
    printf("TP -TG (mss=LEGOS,CLS) \n");
    ss=sget_filteredstatistics2(delta1,delta2, mask, dimD, 3.);
    }

  free(analysis);
  free(mss_cls);
  free(innovation);
  free(model);
  free(tdum);

  free(tmp1);
  free(tmp2);

  __ERR_BASE_LINE__("exiting\n");exit(0);
  
}
