
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

/**-************************************************************************

  altimetry tidal data error

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
#include "filter.h"
#include "tides.h"
#include "netcdf-proto.h"
#include "geo.h"
#include "legend.h"
#include "polygones.h"
#include "statistic.h"
#include "filter.h"
#include "functions.h"

#include "legend.def"

#include "grd.h"

#include "matrix.h"

#include "altimetry.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int reorder(double *x, float *z, int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,m,n;
  double tmp,dum;

  for(k=0;k<nvalues;k++) {
    l=minpos(&(x[k]),nvalues-k)+k;
    if(k!=l) {
       tmp=x[k];
       x[k]=x[l];
       x[l]=tmp;
       dum=z[k];
       z[k]=z[l];
       z[l]=dum;
       }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int map_OI(grid_t grid, T *buffer, double *xdata, double *ydata, T *zdata, T mask,int ndata,grid_t topogrid,float *topo,float topomask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,status=0;
  double d,scale=100.,x,y;
  float h;
  float w,c;
  T s;
  double dx,dy,rx,ry,ddx,ddy,scale2,cs;
  int   imin,imax,jmin,jmax;
  int **list,*card,ii,jj,mm,nlist,partition=20;
  frame_t frame;

/* *----------------------------------------------------------------------
  use partioning to optimize mapping*/

  partition=MIN(20,(int) ((grid.xmax-grid.xmin)/4.0));
  partition=MIN(20,(int) ((grid.ymax-grid.ymin)/4.0));

  rx=4;
  ry=4;

  dx=(grid.xmax-grid.xmin)/(double) partition;
  dy=(grid.ymax-grid.ymin)/(double) partition;

  frame.xmin=grid.xmin-0.1*dx;
  frame.xmax=grid.xmax+0.1*dx;
  frame.ymin=grid.ymin-0.1*dy;
  frame.ymax=grid.ymax+0.1*dy;

  dx=(frame.xmax-frame.xmin)/(double) partition;
  dy=(frame.ymax-frame.ymin)/(double) partition;

//  boxes=partition_frame(frame,partition);

  nlist=partition*partition;

  card=new int[nlist];
  for(l=0;l<nlist;l++) card[l]=0;

  for(n=0;n<ndata;n++) {
    x=xdata[n];
    y=ydata[n];
    i=(int) floor((x-frame.xmin)/dx);
    j=(int) floor((y-frame.ymin)/dy);
    imin=(int) floor((x-rx-frame.xmin)/dx-0.1);
    imax=(int) floor((x+rx-frame.xmin)/dx+0.1);
    jmin=(int) floor((y-ry-frame.ymin)/dy-0.1);
    jmax=(int) floor((y+ry-frame.ymin)/dy+0.1);
    for (jj=MAX(0,jmin);jj<=MIN(jmax,partition-1);jj++) {
      for (ii=MAX(0,imin);ii<=MIN(imax,partition-1);ii++) {
        l=jj*partition+ii;
        card[l]++;
        }
      }
    }

  list=new int*[nlist];
  for(l=0;l<nlist;l++) list[l]=new int[card[l]];

  for(l=0;l<nlist;l++) card[l]=0;
  for(n=0;n<ndata;n++) {
    x=xdata[n];
    y=ydata[n];
    i=(int) floor((x-frame.xmin)/dx);
    j=(int) floor((y-frame.ymin)/dy);
    imin=(int) floor((x-rx-frame.xmin)/dx-0.1);
    imax=(int) floor((x+rx-frame.xmin)/dx+0.1);
    jmin=(int) floor((y-ry-frame.ymin)/dy-0.1);
    jmax=(int) floor((y+ry-frame.ymin)/dy+0.1);
    for (jj=MAX(0,jmin);jj<=MIN(jmax,partition-1);jj++) {
      for (ii=MAX(0,imin);ii<=MIN(imax,partition-1);ii++) {
        l=jj*partition+ii;
        list[l][card[l]]=n;
        card[l]++;
        }
      }
    }

  for(n=0;n<grid.nx*grid.ny;n++) buffer[n]=mask;

  scale2=scale*scale/110./110.;

  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx; i++) {
      x=map_grid_x(grid,i,j);
      y=map_grid_y(grid,i,j);
      cs=cos(y*d2r);
      status=map_interpolation(topogrid,topo,topomask,x,y,&h);
      if(h>0) continue;
      s=0;
      w=0;
      ii=(int) floor((x-frame.xmin)/dx);
      jj=(int) floor((y-frame.ymin)/dy);
      l=jj*partition+ii;
      for(mm=0;mm<card[l];mm++) {
        k=list[l][mm];
//      for(k=0;k<ndata;k++) {
        if(zdata[k]!=mask) {
          ddx=fabs(x-xdata[k])*cs;
          if(ddx>4.0) continue;
          ddy=fabs(y-ydata[k]);
          if(ddy>4.0) continue;
/* *----------------------------------------------------------------------
         sacrifice some precision to optimize mapping*/
          d=(ddx*ddx+ddy*ddy)/scale2;
          c=exp(-d);
//           d=geo_km(x,y,xdata[k],ydata[k])/scale;
//           c=exp(-d*d);
          s+=zdata[k]*c;
          w+=c;
          }
        }
      n=j*grid.nx+i;
      if(w!=0) {
        buffer[n]=s/w;
        }
      else {
        buffer[n]=mask;
        }
      }
    }

  delete[] card;
  for(l=0;l<nlist;l++) delete[] list[l];
  delete [] list;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int coastal_distance(grid_t grid, T *buffer, T mask,float *topo,float topomask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,status=0;
  int imin,imax,jmin,iimin,iimax,jmax,delta_i,delta_j;
  double c,d,scale=100.,x,y;
  char *flag;
  float h;
  int nextk[24]={-1,+1, 0, 0,-1,+1,-1,+1,-2,+2, 0, 0,-2,+2,-2,+2,-1,-1,+1,+1,-2,-2,+2,+2};
  int nextl[24]={ 0, 0,-1,+1,-1,+1,+1,-1, 0, 0,-2,+2,-1,-1,+1,+1,-2,+2,-2,+2,-2,+2,-2,+2};
  T s;


  flag=new char[grid.ny*grid.nx];

  for(n=0;n<grid.nx*grid.ny;n++) buffer[n]=1.e+10;

  for(n=0;n<grid.nx*grid.ny;n++) flag[n]=0;

  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx; i++) {
      n=j*grid.nx+i;
      if(topo[n]>0.) {
        m=0;
        status=-1;
        do {
          k=i+nextk[m];
          l=j+nextl[m];
          m++;
          if(k<0) continue;
          if(l<0) continue;
          if(k>grid.nx-1) continue;
          if(l>grid.ny-1) continue;
          if(topo[l*grid.nx+k]<0.) status=0;
          } while((status!=0)&&(m<9));
        if(status==0) {
          flag[n]=1;
          }
        }
      }
    }

  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx; i++) {
      x=map_grid_x(grid,i,j);
      y=map_grid_y(grid,i,j);
      n=j*grid.nx+i;
      if(topo[n]<0.) {
        imin=i-1;
        while (topo[j*grid.nx+imin]<0.0){
          if(imin==0) break;
          imin--;
          }
        imax=i+1;
        while (topo[j*grid.nx+imax]<0.0){
          if(imax==grid.nx-1) break;
          imax++;
          }
        jmin=j-1;
        while (topo[jmin*grid.nx+i]<0.0){
          if(jmin==0) break;
          jmin--;
          }
        jmax=j+1;
        while (topo[jmax*grid.nx+i]<0.0){
          if(jmax==grid.ny-1) break;
          jmax++;
          }
        delta_i=MIN(imax-i, i-imin);
        delta_j=MIN(jmax-j, j-jmin);
//         for(l=j-delta_j;l<j+delta_j+1;l++) {
//           for(k=i-delta_i;k<i+delta_i+1;k++) {
        iimin=imin;
        iimax=imax;
        for(l=j;l>=jmin;l--) {
          for(k=i;k>=iimin;k--) {
            if(flag[l*grid.nx+k]==1) {
              d=geo_km(x,y,map_grid_x(grid,k,l),map_grid_y(grid,k,l));
              buffer[n]=MIN(d,buffer[n]);
              iimin=k;
//              jmin=l;
              }
            }
          for(k=i;k<iimax+1;k++) {
            if(flag[l*grid.nx+k]==1) {
              d=geo_km(x,y,map_grid_x(grid,k,l),map_grid_y(grid,k,l));
              buffer[n]=MIN(d,buffer[n]);
              iimax=k;
//              jmin=l;
              }
            }
          }
        iimin=imin;
        iimax=imax;
        for(l=j;l<jmax+1;l++) {
          for(k=i;k>=iimin;k--) {
            if(flag[l*grid.nx+k]==1) {
              d=geo_km(x,y,map_grid_x(grid,k,l),map_grid_y(grid,k,l));
              buffer[n]=MIN(d,buffer[n]);
              iimin=k;
//              jmax=l;
              }
            }
          for(k=i;k<iimax+1;k++) {
            if(flag[l*grid.nx+k]==1) {
              d=geo_km(x,y,map_grid_x(grid,k,l),map_grid_y(grid,k,l));
              buffer[n]=MIN(d,buffer[n]);
              iimax=k;
//              jmax=l;
              }
            }
          }
        }
      else {
        buffer[n]=mask;
        }
//      d=geo_km(x,y,xdata[k],ydata[k])/scale;
      }
    }

  delete[] flag;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int map_OI(grid_t grid, T *buffer, double *xdata, double *ydata, T *zdata, double *buble,T mask,int ndata,grid_t topogrid,float *topo,float topomask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,status=0;
  double c,d,scale=100.,x,y;
  double cs,dx,dy,scale2;
  float h;
  T s;

  for(n=0;n<grid.nx*grid.ny;n++) buffer[n]=mask;

  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx; i++) {
      x=map_grid_x(grid,i,j);
      y=map_grid_y(grid,i,j);
      cs=cos(y*d2r);
      n=j*grid.nx+i;
      status=map_interpolation(topogrid,topo,topomask,x,y,&h);
      if(h>0) continue;
      for(k=0;k<ndata;k++) {
        if(zdata[k]!=mask) {
//          scale=buble[k];
          dx=fabs(x-xdata[k])*cs;
          if(dx>4.0) continue;
          dy=fabs(y-ydata[k]);
          if(dy>4.0) continue;
          scale2=buble[k]*buble[k]/110./110.;
//           if(fabs(x-xdata[k])>4.0) continue;
//           if(fabs(y-ydata[k])>4.0) continue;
/* *----------------------------------------------------------------------
         sacrifice some precision to optimize mapping*/
          d=(dx*dx+dy*dy)/scale2;
          c=exp(-d);
//           d=geo_km(x,y,xdata[k],ydata[k])/scale;
//           c=exp(-d*d);
          s=zdata[k]*c;
          if(buffer[n]==mask) {
            buffer[n]=s;
            }
          else {
            buffer[n]=MAX(s,buffer[n]);
            }
          }
        }
      }
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int compute_alias_error(char *wave, int kk, vector<mgr_t> mgr,int nmgr,xover_t *xover,int nxovers, float *alias_error, int **list, int *card)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int *track,*subtrack,ntracks,ntics;
  int npt,i,j,m,n,nn,k,l,ll,target,w,status;
  int *start,*ends;
  int n1,n2,n3;
  int ncycle,ncyclemax;
  float *xover_error,*error;
  double *x,*y,*h,*t,*d,*lon, *lat,tt;
  float *noise, error_min, rmask=-1;
  std::complex<float> *c;

  error   = new float [nmgr];

  subtrack = new int   [nmgr];
  track    = new int   [nmgr];
  x     = new double[nmgr];
  y     = new double[nmgr];

/* *----------------------------------------------------------------------------
  count number of different tracks and register them*/
  ntracks=0;
  for(n = 0; n < nmgr; n++) {
    for (k=0;k<ntracks;k++) {
      if(mgr[n].track==subtrack[k]) break;
      }
    if(k==ntracks) {
      subtrack[ntracks]=mgr[n].track;
      printf("track %d with index %d\n",ntracks,subtrack[ntracks]);
      ntracks++;
      }
    }

/* *----------------------------------------------------------------------------
  initialisation*/
  t=new double[nxovers];
  xover_error=new float [nxovers];

  for(n = 0; n < nmgr; n++) {
    x[n]=geo_recale(mgr[n].loc.lon,0.,180.);
    y[n]=mgr[n].loc.lat;
    alias_error[n]=0.;
    }

/* *----------------------------------------------------------------------------
  examin xover errors*/
  for (k=0;k<ntracks;k++) {
    ntics=0;
    n1=-1;
    n2=-1;
    n3=-1;
    for(l = 0; l < nxovers; l++) {
      if (xover[l].node[2]==-1) continue;
      if (xover[l].node[3]==-1) continue;
/* *----------------------------------------------------------------------------
      if error too large, find an alternative xover*/
      if (abs(xover[l].error[kk])>.05) {
        ncyclemax=0;
        target=-1;
        for(ll = 0; ll < nxovers; ll++) {
          if(xover[ll].id==xover[l].id) {
            ncycle=(mgr[xover[ll].node[0]].duree+mgr[xover[ll].node[1]].duree)/2;
            printf("xover %d error=%f ncycle=%d\n",ll,abs(xover[ll].error[kk]),ncycle);
            if(ll==l) continue;
            if(ncycle>ncyclemax) {
              target=ll;
              ncyclemax=ncycle;
              }
            }
          }
        if(target!=-1) {
          printf("substitue %d (%f) with %d (%f)\n",l,abs(xover[l].error[kk]),target,abs(xover[target].error[kk]));
          xover[l].error[kk]=xover[target].error[kk];
          }
        }
/* *----------------------------------------------------------------------------
      detect true track number when at least 2 xovers available*/
      if(mgr[xover[l].node[0]].track==subtrack[k]) {
        if(abs(xover[l].error[kk])==1.e+10) continue;
        if(abs(xover[l].error[kk])>0.02) continue;
        xover_error[ntics]=abs(xover[l].error[kk]);
        t[ntics]=geo_km(x[list[k][0]],y[list[k][0]],xover[l].t,xover[l].p);
        ntics++;
        if(n1==-1) {
          n1=xover[l].track[0];
          n2=xover[l].track[1];
          }
        else {
          if(n1==xover[l].track[0]) n3=n1;
          if(n2==xover[l].track[1]) n3=n2;
          }
        }
      if(mgr[xover[l].node[1]].track==subtrack[k]) {
        if(abs(xover[l].error[kk])==1.e+10) continue;
        if(abs(xover[l].error[kk])>0.02) continue;
        xover_error[ntics]=abs(xover[l].error[kk]);
        t[ntics]=geo_km(x[list[k][0]],y[list[k][0]],xover[l].t,xover[l].p);
        ntics++;
        if(n1==-1) {
          n1=xover[l].track[0];
          n2=xover[l].track[1];
          }
        else {
          if(n1==xover[l].track[0]) n3=n1;
          if(n2==xover[l].track[1]) n3=n2;
          }
        }
      }
/* *----------------------------------------------------------------------------
    linear interpolation and persistence*/
    if (ntics !=0) {
/* *----------------------------------------------------------------------------
      set position in increasing order*/
      status= reorder(t, xover_error, ntics);
      for(nn = 0; nn < card[k]; nn++) {
        n=list[k][nn];
        track[n]=n3;
        tt=geo_km(x[list[k][0]],y[list[k][0]],x[n],y[n]);
        status=map_interpolate1D(xover_error, t, rmask, ntics, tt, &alias_error[n],1);
//        printf("%d : %lf %lf %f \n",n,x[n],y[n],abs(error[n]));
        }
      printf("#--------------------------------- track portion %d done (part of track %d) ...\n",k,n3);
      }
    else{
      printf("#--------------------------------- track portion %d skipped (part of track %d) ...\n",k,n3);
      }
    }

  delete[] t;
  delete[] xover_error;

  delete[] x;
  delete[] y;
  delete[] error;
  delete[] subtrack;
  delete[] track;


  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int data_error(char *poly,char *wave, int kk, char *zone, vector<mgr_t> mgr,int nmgr,xover_t *xover,int nxovers,grid_t topogrid,float *topo,float topomask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int *track,*subtrack,ntracks,ntics;
  int **linked,*nlinked,**list,*card,count,npt,i,j,m,mm,n,nn,k,l,ll,w,status;
  int n1,n2,n3;
  float *alias_error,*coastal_error,*analysis_error,*error;
  double *x,*y,*h,*t,*d,*lon, *lat,*buble,tt;
  float error_min, *rbuffer[10], rmask=data_rmask,a,p;
  std::complex<float> *cbuffer,cmask=-1;
  std::complex<float> *c;
  grid_t grid;
  float **covariance;
  char obsname[1024];

  cdfvar_t variable;
  pocgrd_t ncgrid;
  char output[1024];

  xerror_t xerror;

  buble   = new double [nmgr];
  error   = new float [nmgr];
  alias_error   = new float [nmgr];
  coastal_error = new float [nmgr];
  analysis_error = new float [nmgr];
  subtrack = new int   [nmgr];
  track    = new int   [nmgr];
  x     = new double[nmgr];
  y     = new double[nmgr];

/* *----------------------------------------------------------------------------
  count number of different tracks and register them*/
  ntracks=0;
  for(n = 0; n < nmgr; n++) {
    for (k=0;k<ntracks;k++) {
      if(mgr[n].track==subtrack[k]) break;
      }
    if(k==ntracks) {
      subtrack[ntracks]=mgr[n].track;
      printf("track %d with index %d\n",ntracks,subtrack[ntracks]);
      ntracks++;
      }
    }

/* *----------------------------------------------------------------------------
  inventory of altimeter gauges by tracks*/
  list=new int*[ntracks];
  card=new int[ntracks];

  for (k=0;k<ntracks;k++) {
    count=0;
    for(n = 0; n < nmgr; n++) {
      if(mgr[n].track==subtrack[k]) {
        count++;
        }
      }
    list[k]=new int[count];
    card[k]=count;
    count=0;
    for(n = 0; n < nmgr; n++) {
      if(mgr[n].track==subtrack[k]) {
        list[k][count]=n;
        count++;
        }
      }
    }

/* *----------------------------------------------------------------------------
  initialisation*/
  for(n = 0; n < nmgr; n++) {
    x[n]=geo_recale(mgr[n].loc.lon,0.,180.);
    y[n]=mgr[n].loc.lat;
    analysis_error[n]=rmask;
    alias_error[n]=rmask;
    coastal_error[n]=0.;
    error[n]=rmask;
    buble[n]=20.;
    }

  xerror.n=nmgr;
  xerror.alias=alias_error;
  xerror.alongtrack=coastal_error;
  xerror.harmonic=analysis_error;

  xerror.cut=new float[nmgr];

  xerror.zraw=new complex<float>[nmgr];
  xerror.zlwf=new complex<float>[nmgr];
  xerror.zhgf=new complex<float>[nmgr];
  xerror.smoothed=new complex<float>[nmgr];

  for(n = 0; n < nmgr; n++) {
    xerror.zraw[n]=data_cmask;
    xerror.zlwf[n]=data_cmask;
    xerror.zhgf[n]=data_cmask;
    }

/* *----------------------------------------------------------------------------
  error estimate: harmonic analysis error bar*/
  status= compute_analysis_error(wave, mgr, nmgr, analysis_error);
  for(n = 0; n < nmgr; n++) {
    if(analysis_error[n]!=rmask) {
      error[n]=analysis_error[n];
      }
    }

/* *----------------------------------------------------------------------------
  error estimate: xover misfit*/
  status= compute_alias_error  (wave, kk, mgr, nmgr, xover, nxovers, alias_error,list,card);
  for(n = 0; n < nmgr; n++) {
    if(alias_error[n]!=rmask) {
      error[n]=MAX(error[n],alias_error[n]);
      }
    }

/* *----------------------------------------------------------------------------
  error estimate: along-track HF variability*/
  status= compute_coastal_error(wave, kk, mgr, nmgr, xerror ,buble,list,card);
  for(n = 0; n < nmgr; n++) {
    if(coastal_error[n]!=rmask) {
      error[n]=MAX(error[n],coastal_error[n]);
      }
    }

/* *----------------------------------------------------------------------------
  Coastal distance*/
  rbuffer[0]=new float[topogrid.nx*topogrid.ny];

/* *----------------------------------------------------------------------------
  Optimal interpolation mapping*/
//   status=map_completegridaxis(&topogrid,2);
//   status=coastal_distance(topogrid, rbuffer[0],rmask,topo,topomask);
//
//   sprintf(output,"distance.nc",wave);
//   status= poc_createfile(output);
//   status=poc_sphericalgrid_xy(output,"",topogrid,&ncgrid);
//
//   poc_standardvariable_xy(&variable,"coastal_distance",rmask,"km",1., 0.,"coastal distance","coastal distance","coastal distance",ncgrid);
//   status=create_ncvariable(output, &variable);
//   status=poc_write_xy(output, topogrid, variable.id, rbuffer[0]);
//   status=free_ncvariable(&variable);

  delete[] rbuffer[0];

//  goto covariances;

/* *----------------------------------------------------------------------------
  Map errors*/
//  zone=strdup("persian");
  grid=get_zonegrid(zone);
  status=map_completegridaxis(&grid,2);

  rbuffer[0]=new float[grid.nx*grid.ny];
  rbuffer[1]=new float[grid.nx*grid.ny];
  rbuffer[2]=new float[grid.nx*grid.ny];
  rbuffer[3]=new float[grid.nx*grid.ny];
  cbuffer=new std::complex<float> [grid.nx*grid.ny];

/* *----------------------------------------------------------------------------
  Optimal interpolation mapping*/
  status= map_OI(grid, rbuffer[0],x,y,alias_error,rmask,nmgr,topogrid,topo,topomask);

  sprintf(output,"%s.error.%s.nc",wave,zone);
  status= poc_createfile(output);
  status=poc_sphericalgrid_xy(output,"",grid,&ncgrid);

  poc_standardvariable_xy(&variable,"error_alias",rmask,"m",1., 0.,"error","error","error",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output, grid, variable.id, rbuffer[0]);
  variable.destroy();

covariances:

/* *----------------------------------------------------------------------------
  prepare track correlation*/
  linked=new int*[nmgr];
  covariance=new float*[nmgr];
  nlinked=new int[nmgr];
  for(n = 0; n < nmgr; n++) {
    nlinked[n]=0;
    if(alias_error[n]=rmask) {
      status=map_interpolation(grid,rbuffer[0],rmask,x[n],y[n],&alias_error[n]);
      }
    }
  for(n = 0; n < nmgr; n++) {
    if(alias_error[n]!=rmask) {
/* *----------------------------------------------------------------------------
      initialize with alias error */
      error[n]=alias_error[n];
      }
    else {
/* *----------------------------------------------------------------------------
      if no alias error, initialize with arbitrary error */
      error[n]=0.05;
      }
    if(coastal_error[n]!=rmask) error[n]=MAX(error[n],coastal_error[n]);
//    if(coastal_error[n]!=rmask)  error[n]=coastal_error[n];
//    if(analysis_error[n]!=rmask) error[n]=analysis_error[n];
    w=mgr[n].wave_index(wave);
    if(w!=-1) {
      mgr[n].data[w].error=100.*error[n]/mgr[n].data[w].amp;
      }
    }
/* *----------------------------------------------------------------------------
  filter data with less than 5% error tolerance*/
  status=mgr_save("clean.mgr" ,mgr, 0.05);

/* *----------------------------------------------------------------------------
  variance and covariance computation*/
  for (k=0;k<ntracks;k++) {
    for(nn = 0; nn < card[k]; nn++) {
      n=list[k][nn];
      linked[n]=new int[card[k]];
      covariance[n]=new float[card[k]];
      linked[n][0]=n;
      covariance[n][0]=error[n]*error[n];
      nlinked[n]=1;
      for(mm = 0; mm < card[k]; mm++) {
        m=list[k][mm];
        if(m!=n) {
          linked[n][nlinked[n]]=m;
//          tt=geo_km(x[m],y[m],x[n],y[n])/100.;
          tt=geo_km(x[m],y[m],x[n],y[n])/sqrt(buble[n]*buble[m]);
          covariance[n][nlinked[n]]=error[m]*error[n]*exp(-tt);
          nlinked[n]++;
          }
        }
      }
    }

  sprintf(obsname,"%s.obs",wave);
/* *----------------------------------------------------------------------------
  COMAPI : no covariance !!! */
//  status=mgr_save_obs4assim(obsname,poly,wave, nmgr ,mgr, covariance, linked, nlinked);
  status=mgr_save_obs4assim(obsname,poly,wave, nmgr ,mgr, covariance, 0, 0);

  sprintf(obsname,"%s.lgd",wave);
  {
  float *all[4];
  all[0]=analysis_error;
  all[1]=alias_error;
  all[2]=coastal_error;
  all[3]=error;
  status=mgr_Constants2Legend_00(obsname,wave, nmgr ,mgr, all, covariance, track, linked, nlinked);
  }
/* *----------------------------------------------------------------------------
  Optimal interpolation mapping of harmonic analysis errors*/
  status= map_OI(grid, rbuffer[2],x,y,analysis_error,buble,rmask,nmgr,topogrid,topo,topomask);

  poc_standardvariable_xy(&variable,"error_analysis",rmask,"m",1., 0.,"error","error","error",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output, grid, variable.id, rbuffer[2]);
  variable.destroy();

/* *----------------------------------------------------------------------------
  Optimal interpolation mapping of along-track errors*/
//  status= map_OI(grid, rbuffer[1],x,y,coastal_error,buble,(float) 0.0,nmgr,topogrid,topo,topomask);
  status= map_OI(grid, rbuffer[1],x,y,coastal_error,buble,rmask,nmgr,topogrid,topo,topomask);

  poc_standardvariable_xy(&variable,"error_alongtrack",rmask,"m",1., 0.,"error","error","error",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output, grid, variable.id, rbuffer[1]);
  variable.destroy();

  for (n=0;n<grid.nx*grid.ny;n++) {
    if((rbuffer[0][n]!=rmask) && (rbuffer[1][n]!=rmask)) rbuffer[3][n]=MAX(rbuffer[0][n],rbuffer[1][n]);
    else rbuffer[3][n]=rmask;
    }
  poc_standardvariable_xy(&variable,"error_integral",rmask,"m",1., 0.,"error","error","error",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output, grid, variable.id, rbuffer[3]);
  variable.destroy();

  for(n = 0; n < nmgr; n++) {
    w=mgr[n].wave_index(wave);
    if(w!=-1) {
      a=abs(xerror.zlwf[n]);
      p=arg(xerror.zlwf[n])*r2d;
      mgr[n].data[w].amp=a;
      mgr[n].data[w].phi=p;
      mgr[n].data[w].error=abs(xerror.zhgf[n])/sqrt(xerror.cut[n]/6.0);
      }
    }

  delete[] x;
  delete[] y;
  delete[] buble;
  delete[] error;
  delete[] coastal_error;
  delete[] alias_error;
  delete[] analysis_error;
  delete[] subtrack;
  delete[] track;

  for(k=0;k<4;k++) delete[] rbuffer[k];

  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double mask=999.9;
  double tau;
  double *residual,*buffer;
  double  zi,zr,mean[10],rms[10],count;
  double  **a=NULL,**G=NULL,a1,p1,a2,p2,d,da,dG;
  int i,j,k,kk,l,m,n,fmt,status,tg;
  FILE *in=NULL,*out=NULL;
  char *report,hfile[1204],tmp[256], datafile[256], lgdfile[1024];
  char *data, *input,*output=NULL, *keyword, *poly=NULL,*mgrfile=NULL,*zone=NULL,*xover_base=NULL;
  fcomplex meanc;
  char c;
  spectrum_t spectrum,reference_spectrum;
  char *wave[100];
  int found,nwave=0;
  date_t date;
  vector<mgr_t> mgr;
  int i1,i2,m1,m2,n1,n2,nmgr;
  char *atlas_directory=NULL,*atlas_convention=NULL;
  fcomplex cmisfit;
  int *list,indice[4];
  double t1,t2,t_xover,p_xover;
//  double p1,p2;
  double dmin;
  int neighbour,nlegends;
  legend_t   *legends;
  legend02_t *legends02;

  char *bathymetry=NULL;
  grid_t topogrid;
  float *topo=NULL,topomask,h;
  double x,y,x1[2],y1[2],x2[2],y2[2],z,radius=10.0;
  std::complex<float> tide[2],c1,c2,cmask;
  int maxcycles=1;
  xover_t *xover;
  xoverbase_t full;
  int nxover;
  char *comments;

  fprintf(stderr,"%s  -starting computation *********\n",argv[0]);

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        if(strcmp(keyword,"-maxcycles")==0) {
          sscanf(argv[n+1],"%d",&maxcycles);
          n++;
          n++;
          continue;
          }
        switch (keyword[1]) {
        case 'm' :
          atlas_convention= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'i' :
          input= strdup(argv[n+1]);
          n++;
          n++;
          break;

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

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'c' :
          mgrfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'r' :
          sscanf(argv[n+1],"%lf",&radius);
          n++;
          n++;
          break;

        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'x' :
          xover_base= strdup(argv[n+1]);
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
        wave[nwave]= strdup(argv[n]);
        printf("input wave=%s\n",wave[nwave]);
        spectrum.n=nwave+1;
        nwave++;
        n++;
        break;
      }
      free(keyword);
    }

  full=read_xover(xover_base);

/* *----------------------------------------------------------------------------
  Initialise tidal library */
  date.year=1950;
  date.month=1;
  date.day=1;
  date.second=0.0;
  astro_angles_t astro_angles;
  reference_spectrum = initialize_tide(&astro_angles,date);

  spectrum.nmax=nwave+1;
  spectrum.n=0;
  exitIfNull(spectrum.waves=new tidal_wave[spectrum.nmax]);
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
    printf ("wave: %10s, pulsation: %12.6f degrees/h \n", spectrum.waves[i].name,spectrum.waves[i].omega);
    }

  for (i=0; i<spectrum.n; i++) {
    for (j=i+1; j<spectrum.n; j++) {
      tau=deltaOmega2separation(spectrum.waves[j].omega-spectrum.waves[i].omega);
      printf ("wave: %10s %10s, separation: %9.3f days \n",
        spectrum.waves[i].name,spectrum.waves[j].name,tau);
      }
    }

/* *----------------------------------------------------------------------------
  Load tide gauge database */
  if(mgrfile==NULL) goto error;
  strcpy(datafile,mgrfile);
  nmgr=mgr_load(datafile, mgr);

/* *----------------------------------------------------------------------------
  Warning : code will failed if track number not in station name*/
  status=mgr_save("filtered.mgr", mgr, (float) 20.0);

  nmgr=mgr_load("filtered.mgr", mgr);

  list= mgr_order(mgr, nmgr);

  exitIfNull(
    a=(double **) malloc(nwave*sizeof(double))
    );
  
  exitIfNull(
    G=(double **) malloc(nwave*sizeof(double))
    );

  for (k=0;k<nwave;k++) {
    exitIfNull(
      a[k]=(double *) malloc(nmgr*sizeof(double))
      );
    exitIfNull(
      G[k]=(double *) malloc(nmgr*sizeof(double))
      );
    }

  if(output==NULL) {
    report=strdup("report.out");
    }
  else {
    report=strdup(output);
    }

/* *----------------------------------------------------------------------------
  Load bottom topography */
  if(bathymetry==NULL) bathymetry=strdup("/home/softs/genesis/data/topography/data/legos.2.grd");
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

/* *----------------------------------------------------------------------------
  initialise legend for display outputs */
  legends=new legend_t[nwave];

  l=0;
  legends[l].ID=l;
  legends[l].Type=LGD_GRAPH;
  legends[l].T=0;
  legends[l].P=0;
  legends[l].X=0;
  legends[l].Y=0;

  legends02=new legend02_t;
  legends02->np=nmgr;
  legends02->nz=5;
  legends02->pen0=0;
  legends02->pen1=0;
  legends02->points=new point_t[nmgr];

  legends[l].ptr=(char *) legends02;

  nlegends=1;

  comments=new char[1024];
  comments[0]=0;

  strcat(comments,"#  0: xover index \n");
  strcat(comments,"#  1: longitude of xover \n");
  strcat(comments,"#  2: latitude  of xover \n");
  strcat(comments,"#  3->0: available cycles \n");
  strcat(comments,"#  4->1: amplitude \n");
  strcat(comments,"#  5->2: phase lag \n");
  strcat(comments,"#  5->3: analysis error \n");
  strcat(comments,"#  6->4: available cycles (%)\n");

  for (m=0;m<nmgr;m++) {
    legends02->points[m].z=new float[legends02->nz];
    }

/* *----------------------------------------------------------------------------
  Interpolate depth at tide gauge*/
  for (m=0;m<nmgr;m++) {
    x=map_recale(topogrid,mgr[m].loc.lon);
    y=mgr[m].loc.lat;
    status=map_interpolation(topogrid,topo,topomask,x,y,&h);
    mgr[m].loc.depth=h;
    }

/* *----------------------------------------------------------------------------
  create legend with all valid space gauges*/
  for (k=0;k<nwave;k++) {
    l=0;
    for (m=0;m<nmgr;m++) {
      i=mgr[m].wave_index(wave[k]);
      if(i==-1) continue;
      for(kk = 0; kk < reference_spectrum.n; kk++) {
        if(strcmp(reference_spectrum.waves[kk].name, wave[k]) == 0) {
          break;
          }
        }
      if(kk != reference_spectrum.n) {
        mgr[m].data[i].constituent=reference_spectrum.waves[kk];
        }
      j=0;
      legends02->points[l].t=mgr[m].loc.lon;
      legends02->points[l].p=mgr[m].loc.lat;
      a1=mgr[m].data[i].amp;
      p1=mgr[m].data[i].phi;
      legends02->points[l].z[j++]=mgr[m].duree;
      legends02->points[l].N=mgr[m].track;
      legends02->points[l].z[j++]=a1;
      legends02->points[l].z[j++]=p1;
      legends02->points[l].z[j++]=mgr[m].data[i].error;
      legends02->points[l].z[j++]=100.*mgr[m].duree/(float) maxcycles;
      l++;
      }
    sprintf(lgdfile,"along-track-%s.lgd",wave[k]);
    legends02->np=l;
    status=lgd_save(lgdfile, legends, nlegends,NULL,NULL);
    }

/* *----------------------------------------------------------------------------
  find xover primary nodes */
  xover=new xover_t[nmgr];
  nxover=0;
  for (m=0;m<nmgr;m++) {
    t1=mgr[m].loc.lon;
    p1=mgr[m].loc.lat;
//     i1=mgr[m].wave_index(wave[k]);
//     if(i1==-1) continue;
    dmin=1.e+10;
    neighbour=-1;
    for (n=m+1;n<nmgr;n++) {
      t2=mgr[n].loc.lon;
      t2=geo_recale(t2,t1,(double) 180.0);
//       i2=mgr[n].wave_index(wave[k]);
//       if(i2==-1) continue;
/* *----------------------------------------------------------------------------
     to be revised for very high latitudes*/
      if(fabs(t2-t1) > 1.0)  continue;
      p2=mgr[n].loc.lat;
      if(fabs(p2-p1) > 1.0)  continue;
      d=geo_km(t1,p1,t2,p2);
      if((d<dmin) && (n!=m) && (n!=m+1) && (n!=m-1)) {
        neighbour=n;
        dmin=d;
        t_xover=0.5*(t1+t2);
        p_xover=0.5*(p1+p2);
        }
      }
    if((dmin<radius) && (neighbour!=m) && (neighbour!=m+1) && (neighbour!=m-1)) {
      nxover++;
      }
    else neighbour=-1;
    if(neighbour==-1) continue;

    n=neighbour;
    l=nxover-1;
    xover[l].node[0]=m;
    xover[l].node[1]=n;
    xover[l].node[2]=-1;
    xover[l].node[3]=-1;
    xover[l].t=t_xover;
    xover[l].p=p_xover;
    xover[l].d=dmin;
    xover[l].error=new std::complex<float> [nwave];
    xover[l].c[0]=new complex<float>[nwave];
    xover[l].c[1]=new complex<float>[nwave];
    status= xover_id(full,&(xover[l]));
/* *----------------------------------------------------------------------------
    desperate attempt to get the track number */
    if(status==-1) {
      k=strlen(mgr[m].name)-7;
      sscanf(&(mgr[m].name[k]),"%3d",&(xover[l].track[0]));
      k=strlen(mgr[n].name)-7;
      sscanf(&(mgr[n].name[k]),"%3d",&(xover[l].track[1]));
      if(xover[l].track[0]<xover[l].track[1]) {
        xover[l].id=252*xover[l].track[1]+xover[l].track[0];
        }
      else {
        xover[l].id=252*xover[0].track[1]+xover[l].track[1];
        }
      }
    printf("basic xover (%d, CTOH id=%d) found for tracks %d %d (%d %d) \n",l,xover[l].id,xover[l].track[0],xover[l].track[1],mgr[m].track,mgr[n].track);
    }
  printf("%d basic xover found\n",nxover);

/* *----------------------------------------------------------------------------
  create legend with all valid primary xovers */
  legends02->nz=13;
  legends02->pen0=0;
  legends02->pen1=0;

  for (m=0;m<nmgr;m++) {
    delete [] legends02->points[m].z;
    legends02->points[m].z=new float[legends02->nz];
    }

  if ((out = fopen(report, "w")) == NULL) {
      __ERR_BASE_LINE__("file opening issue : %s \n",report);
      exit(-1);
    }

  comments=new char[1024];
  comments[0]=0;

  strcat(comments,"#  0: xover index \n");
  strcat(comments,"#  1: longitude of xover \n");
  strcat(comments,"#  2: latitude  of xover \n");
  strcat(comments,"#  3: amplitude 1st track \n");
  strcat(comments,"#  4: phase lag 1st track \n");
  strcat(comments,"#  5: amplitude 2nd track \n");
  strcat(comments,"#  6: phase lag 2nd track \n");
  strcat(comments,"#  7: absolute misfit (m) \n");
  strcat(comments,"#  8: xover depth \n");
  strcat(comments,"#  9: available cycles node 1 \n");
  strcat(comments,"# 10: available cycles node 2 \n");
  strcat(comments,"# 11: available cycles node 3 \n");
  strcat(comments,"# 12: available cycles node 4 \n");
  strcat(comments,"# 13: relative misfit (%) \n");
  strcat(comments,"# 14: mean amplitude\n");

  for (k=0;k<nwave;k++) {
    legends02->np=0;
    for (l=0;l<nxover;l++) {
      i1=-1;
      i2=-1;
      n1=xover[l].node[0];
      i1=mgr[n1].wave_index(wave[k]);
      n2=xover[l].node[1];
      i2=mgr[n2].wave_index(wave[k]);
      if((i1==-1) || (i2==-1)) {
/* *----------------------------------------------------------------------------
        mask error if wave not available */
        xover[l].error[k]=1.e+10;
//        xover[l].error[k]=-1;
        continue;
        }
      legends02->np++;
      n=legends02->np-1;
      legends02->points[n].t=xover[l].t;
      legends02->points[n].p=xover[l].p;
      a1=mgr[n1].data[i1].amp;
      p1=    mgr[n1].data[i1].phi;
      if(p1<  0) p1+=360.0;
      if(p1>360) p1-=360.0;
      a2=mgr[n2].data[i2].amp;
      p2=    mgr[n2].data[i2].phi;
      if(p2<  0) p2+=360.0;
      if(p2>360) p2-=360.0;
      zr=a2*cos(-p2*d2r)-a1*cos(-p1*d2r);
      zi=a2*sin(-p2*d2r)-a1*sin(-p1*d2r);
      d=sqrt(zr*zr+zi*zi);
      x=map_recale(topogrid,xover[l].t);
      y=xover[l].p;
      status=map_interpolation(topogrid,topo,topomask,x,y,&h);
#ifdef ECHO
      printf("%5d %5d lon= %6.2lf lat= %6.2lf depth= %7.1lf d=%6.2lf km: a1= %7.1f G1= %7.1f a2= %7.1f G2= %7.1f - e= %7.1lf\n",
             n1,n2,x,y,h,xover[l].d,a1,p1,a2,p2,d);
#endif
      fprintf(out,"%5d %5d lon= %6.2lf lat= %6.2lf depth= %7.1lf d=%6.2lf km: a1= %7.1f G1= %7.1f a2= %7.1f G2= %7.1f - e= %7.1lf\n",
             n1,n2,x,y,h,xover[l].d,a1,p1,a2,p2,d);
      xover[l].c[0][k]=fct_polar2complex(a1,p1);
      xover[l].c[1][k]=fct_polar2complex(a2,p2);
      xover[l].error[k]=xover[l].c[1][k]-xover[l].c[0][k];
      j=0;
      legends02->points[n].N=l;
      legends02->points[n].z[j++]=a1;
      legends02->points[n].z[j++]=p1;
      legends02->points[n].z[j++]=a2;
      legends02->points[n].z[j++]=p2;
      legends02->points[n].z[j++]=d;
      legends02->points[n].z[j++]=xover[l].d;
      legends02->points[n].z[j++]=mgr[xover[l].node[0]].duree;
      legends02->points[n].z[j++]=mgr[xover[l].node[1]].duree;
      legends02->points[n].z[j++]=mgr[xover[l].node[0]].duree;
      legends02->points[n].z[j++]=mgr[xover[l].node[1]].duree;
      legends02->points[n].z[j++]=100.*d/fabs(a1+a2);
      legends02->points[n].z[j++]=0.5*fabs(a1+a2);
      legends02->points[n].z[j++]=fabs(arg(xover[l].error[k])*r2d);
      }
    sprintf(lgdfile,"xover-%s.lgd",wave[k]);
    status=lgd_save(lgdfile, legends, nlegends,NULL,comments);
    }
  fclose(out);


/* *----------------------------------------------------------------------------
  find secondary nodes */
  for (l=0;l<nxover;l++) {
    m1=xover[l].node[0];
    m2=xover[l].node[1];
    xover[l].node[2]=-1;
    xover[l].node[3]=-1;
    x1[0]=mgr[m1].loc.lon;
    y1[0]=mgr[m1].loc.lat;
    x2[0]=mgr[m2].loc.lon;
    x2[0]=geo_recale(x2[0],x1[0],(double) 180.0);
    y2[0]=mgr[m2].loc.lat;
    d=geo_km(x1[0],y1[0],x2[0],y2[0]);
    if(d<0.1) {
      x=0.5*(x1[0]+x1[1]);
      y=0.5*(y1[0]+y1[1]);
      xover[l].node[2]=m1+1;
      xover[l].node[3]=m2+1;
      xover[l].t=x;
      xover[l].p=y;
//       printf("%5d %5d %5d lon= %8.4lf lon= %8.4lf lon= %8.4lf lat= %8.4lf lat= %8.4lf lat= %8.4lf\n",l,m1,n1,x1[0],x1[1],x,y1[0],y1[1],y);
//       printf("%5d %5d %5d lon= %8.4lf lon= %8.4lf lon= %8.4lf lat= %8.4lf lat= %8.4lf lat= %8.4lf\n",l,m2,n2,x2[0],x2[1],x,y2[0],y2[1],y);
      }
    for (j=-1;j<3;j+=2) {
      n2=m2+j;
      if( (n2<0) || (n2==nmgr) ) continue;
      x2[1]=mgr[n2].loc.lon;
      x2[1]=geo_recale(x2[1],x2[0],(double) 180.0);
      y2[1]=mgr[n2].loc.lat;
      d=geo_km(x2[0],y2[0],x2[1],y2[1]);
      if(d>10.) continue;
      for (i=-1;i<3;i+=2) {
        n1=m1+i;
        if( (n1<0) || (n1==nmgr) ) continue;
        x1[1]=mgr[n1].loc.lon;
        x1[1]=geo_recale(x1[1],x1[0],(double) 180.0);
        y1[1]=mgr[n1].loc.lat;
        d=geo_km(x1[0],y1[0],x1[1],y1[1]);
        if(d>10.) continue;
        if(plg_secantpoint(x1,y1,x2,y2,&x,&y)==PLG_LINES_SECANT) {
          xover[l].node[2]=n1;
          xover[l].node[3]=n2;
          xover[l].t=x;
          xover[l].p=y;
#ifdef ECHO
          printf("%5d %5d %5d lon= %8.4lf lon= %8.4lf lon= %8.4lf lat= %8.4lf lat= %8.4lf lat= %8.4lf\n",l,m1,n1,x1[0],x1[1],x,y1[0],y1[1],y);
          printf("%5d %5d %5d lon= %8.4lf lon= %8.4lf lon= %8.4lf lat= %8.4lf lat= %8.4lf lat= %8.4lf\n",l,m2,n2,x2[0],x2[1],x,y2[0],y2[1],y);
#endif
          printf("optimal xover (%d, CTOH id=%d) found for tracks %d %d (%d %d %d %d) \n",l,xover[l].id,xover[l].track[0],xover[l].track[1],m1,n1,m2,n2);
          }
        }
      }
    }

/* *----------------------------------------------------------------------------
  remove redundancy*/
  for (l=0;l<nxover;l++) {
    m1=xover[l].node[0];
    m2=xover[l].node[1];
    n1=xover[l].node[2];
    n2=xover[l].node[3];
    if(n1==-1) continue;
    if(n2==-1) continue;
    for (k=l+1;k<nxover;k++) {
      if(((xover[k].node[0]==n1) && (xover[k].node[2]==m1)) ||
         ((xover[k].node[1]==n2) && (xover[k].node[3]==m2))) {
        xover[k].node[2]=-1;
        xover[k].node[3]=-1;
        break;
        }
      }
    }

/* *----------------------------------------------------------------------------
  create legend with all valid secondary xovers */

  if ((out = fopen(report, "w")) == NULL) {
      __ERR_BASE_LINE__("file opening issue : %s \n",report);
      exit(-1);
    }

  comments=new char[1024];
  comments[0]=0;

  strcat(comments,"#  0: xover index \n");
  strcat(comments,"#  1: longitude of xover \n");
  strcat(comments,"#  2: latitude  of xover \n");
  strcat(comments,"#  3: amplitude 1st track \n");
  strcat(comments,"#  4: phase lag 1st track \n");
  strcat(comments,"#  5: amplitude 2nd track \n");
  strcat(comments,"#  6: phase lag 2nd track \n");
  strcat(comments,"#  7: absolute misfit (m) \n");
  strcat(comments,"#  8: xover depth \n");
  strcat(comments,"#  9: available cycles node 1 \n");
  strcat(comments,"# 10: available cycles node 2 \n");
  strcat(comments,"# 11: available cycles node 3 \n");
  strcat(comments,"# 12: available cycles node 4 \n");
  strcat(comments,"# 13: relative misfit (%) \n");
  strcat(comments,"# 14: mean amplitude\n");

  for (k=0;k<nwave;k++) {
    legends02->ID=k;
    legends02->np=0;
    for (l=0;l<nxover;l++) {
      if(xover[l].node[2]==-1) continue;
      for(j=0;j<4;j++) {
        n=xover[l].node[j];
        indice[j]=mgr[n].wave_index(wave[k]);
        }
      if((indice[0]==-1) || (indice[1]==-1)|| (indice[2]==-1)|| (indice[3]==-1)) continue;
      legends02->np++;
      n=legends02->np-1;
      legends02->points[n].t=xover[l].t;
      legends02->points[n].p=xover[l].p;
      n1=xover[l].node[0];
/* *----------------------------------------------------------------------------
      interpolate tide at xover exact location*/
      tide[0]=fct_polar2complex(mgr[n1].data[indice[0]].amp,mgr[n1].data[indice[0]].phi);
      n2=xover[l].node[2];
      tide[1]=fct_polar2complex(mgr[n2].data[indice[2]].amp,mgr[n2].data[indice[2]].phi);
      x1[0]=mgr[n1].loc.lon;
      x1[1]=mgr[n2].loc.lon;
      x1[1]=geo_recale(x1[1],x1[0],(double) 180.0);
      x=xover[l].t;
      x=geo_recale(x,x1[0],(double) 180.0);
      status=map_interpolate1D(tide, x1, cmask, 2, x, &c1);
      n1=xover[l].node[1];
      tide[0]=fct_polar2complex(mgr[n1].data[indice[1]].amp,mgr[n1].data[indice[1]].phi);
      n2=xover[l].node[3];
      tide[1]=fct_polar2complex(mgr[n2].data[indice[3]].amp,mgr[n2].data[indice[3]].phi);
      x1[0]=mgr[n1].loc.lon;
      x1[1]=mgr[n2].loc.lon;
      x1[1]=geo_recale(x1[1],x1[0],(double) 180.0);
/* *----------------------------------------------------------------------------
      interpolate tide at xover exact location*/
      x=geo_recale(x,x1[0],(double) 180.0);
      status=map_interpolate1D(tide, x1, cmask, 2, x, &c2);
      d=abs(c2-c1);
      x=map_recale(topogrid,xover[l].t);
      y=xover[l].p;
      status=map_interpolation(topogrid,topo,topomask,x,y,&h);
#ifdef ECHO
      printf("%5d %5d lon= %6.2lf lat= %6.2lf depth= %7.1lf d=%6.2lf km: a1= %7.1f G1= %7.1f a2= %7.1f G2= %7.1f - e= %7.1lf\n",
             n1,n2,x,y,h,xover[l].d,a1,p1,a2,p2,d);
      fprintf(out,"%5d %5d lon= %6.2lf lat= %6.2lf depth= %7.1lf d=%6.2lf km: a1= %7.1f G1= %7.1f a2= %7.1f G2= %7.1f - e= %7.1lf\n",
             n1,n2,x,y,h,xover[l].d,a1,p1,a2,p2,d);
#endif
      xover[l].c[0][k]=c1;
      xover[l].c[1][k]=c2;
/* *----------------------------------------------------------------------------
      update legend*/
      j=0;
      legends02->points[n].N=xover[l].id;
      a1=abs(c1);
      legends02->points[n].z[j++]=a1;
      p1=atan2(c1.imag(),c1.real())*r2d;
      if(p1<0.0) p1+=360.0;
      legends02->points[n].z[j++]=p1;
      a2=abs(c2);
      legends02->points[n].z[j++]=a2;
      p2=atan2(c2.imag(),c2.real())*r2d;
      if(p2<0.0) p2+=360.0;
      legends02->points[n].z[j++]=p2;
      legends02->points[n].z[j++]=d;
      legends02->points[n].z[j++]=h;
      legends02->points[n].z[j++]=mgr[xover[l].node[0]].duree;
      legends02->points[n].z[j++]=mgr[xover[l].node[2]].duree;
      legends02->points[n].z[j++]=mgr[xover[l].node[1]].duree;
      legends02->points[n].z[j++]=mgr[xover[l].node[3]].duree;
      legends02->points[n].z[j++]=100.*abs(c2-c1)/fabs(a1+a2);
      legends02->points[n].z[j++]=0.5*fabs(a1+a2);
      legends02->points[n].z[j++]=fabs(arg(c2-c1)*r2d);
      xover[l].error[k]=c2-c1;
      }
    sprintf(lgdfile,"xover-improved-%s.lgd",wave[k]);
    status=lgd_save(lgdfile, legends, nlegends,NULL,comments);
    status=data_error(poly,wave[k],k,zone,mgr,nmgr,xover,nxover,topogrid,topo,topomask);
    }
  fclose(out);

  for (k=0;k<nwave;k++) {
    free(a[k]);
    free(G[k]);
    }
  free(a);
  free(G);
  free(topo);

  status=mgr_save("smoothed.mgr", mgr, (float) 0.10);

  __ERR_BASE_LINE__("%s -computation complete ^^^^^^^^^\n",argv[0]);
  exit(0);
error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
}
