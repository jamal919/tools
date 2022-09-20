
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
//#include <complex.h>
 
#include "tools-define.h"
#include "tools-structures.h"

#include "poc-time.h"
#include "list.h"
#include "map.h"
#include "mgr.h"
#include "filter.h"
#include "tides.h"
#include "geo.h"
#include "legend.h"
#include "polygones.h"
#include "functions.h"

#include "legend.def"

#include "grd.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double mask=999.9;
  double tau;
  double *residual,*buffer;
  double  zi,zr,mean[10],rms[10],count;
  double  **a,**G,a1,p1,a2,p2,d,da,dG;
  int i,j,k,l,m,n,fmt,status,tg;
  FILE *in,*out;
  char *report,hfile[1204],tmp[256], datafile[256], lgdfile[1024];
  char *data, *input,*output=NULL, *keyword, *path=NULL,*mgrfile=NULL,*xover_base=NULL;
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
  float *topo,topomask,h;
  double x,y,x1[2],y1[2],x2[2],y2[2],z,radius=10.0;
  std::complex<float> tide[2],c1,c2,cmask;

  xover_t *xover;
  xoverbase_t full;
  int nxover;
  int ncycle,ncyclemax=0,target,ll;
  char *comments;
  int  *redundant,nredundants;

  fprintf(stderr,"%s  -starting computation *********\n",argv[0]);

  fct_echo( argc, argv);

//  full=read_xover("/home/softs/data/gauges/all_xing.dat");

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
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
          path= strdup(argv[n+1]);
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

  if(xover_base==NULL) {
    __OUT_BASE_LINE__("cannot read cross-over location database ( missing -x [database path/name] ?)...\n");
    exit(-1);
    }
  full=read_xover(xover_base);

  if(path==NULL) path=getenv("PWD");
  atlas_directory=strdup(path);

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
  strcpy(datafile,mgrfile);
  nmgr=mgr_load(datafile, mgr);

  list= mgr_order(mgr, nmgr);

  a=(double **) malloc(nwave*sizeof(double));
  G=(double **) malloc(nwave*sizeof(double));

  for (k=0;k<nwave;k++) {
    a[k]=(double *) malloc(nmgr*sizeof(double));
    G[k]=(double *) malloc(nmgr*sizeof(double));
    }

  if(output==NULL) {
    report=strdup("xover.out");
    }
  else {
    report=strdup(output);
    }

/* *----------------------------------------------------------------------------
  Load bottom topography */
  if(bathymetry==NULL) bathymetry=strdup("/home/softs/genesis/data/topography/data/legos.2.grd");
  status=grd_loadgrid(bathymetry,&topogrid);
  if(status !=0) {
    __OUT_BASE_LINE__("cannot load bathymetry file=%f\n",bathymetry);
    exit(-1);
    }
  topo= (float *) malloc(topogrid.nx*topogrid.ny*sizeof(float));
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
  legends02->nz=3;
  legends02->pen0=0;
  legends02->pen1=0;
  legends02->points=new point_t[nmgr];

  legends[l].ptr=(char *) legends02;

  nlegends=1;

  for (m=0;m<nmgr;m++) {
    legends02->points[m].z=new float[2*nwave];
    }

/* *----------------------------------------------------------------------------
  create legend with all valid space gauge*/
  comments=new char[1024];
  comments[0]=0;

  strcat(comments,"#  0: track index \n");
  strcat(comments,"#  1: longitude of space gauge \n");
  strcat(comments,"#  2: latitude  of space gauge \n");
  strcat(comments,"#  3: amplitude\n");
  strcat(comments,"#  4: phase lag\n");
  strcat(comments,"#  5: number of valid measurements\n");

  for (k=0;k<nwave;k++) {
    l=0;
    for (m=0;m<nmgr;m++) {
      i=mgr[m].wave_index(wave[k]);
      if(i==-1) {
        continue;
        }
      else {
        j=0;
        legends02->points[l].t=mgr[m].loc.lon;
        legends02->points[l].p=mgr[m].loc.lat;
        a1=mgr[m].data[i].amp;
        p1=mgr[m].data[i].phi;
        legends02->points[l].N=mgr[m].track;
        legends02->points[l].z[j++]=a1;
        legends02->points[l].z[j++]=p1;
        legends02->points[l].z[j++]=mgr[m].duree;
        l++;
        }
      }
    sprintf(lgdfile,"along-track-%s.lgd",wave[k]);
    legends02->np=l;
    status=lgd_save(lgdfile, legends, nlegends,NULL,comments);
    }
  delete[] comments;

/* *----------------------------------------------------------------------------
  find xover primary nodes */
  xover=new xover_t[nmgr];
  nxover=0;
  int nRequestedProcs=-1;
  int nprocs=initialize_OPENMP(nRequestedProcs);
  for (m=0;m<nmgr;m++) {
//   for (m=0;m<1000;m++) {
    t1=mgr[m].loc.lon;
    p1=mgr[m].loc.lat;
    dmin=1.e+10;
    neighbour=-1;
#pragma omp parallel for
    for (n=m+1;n<nmgr;n++) {
      double t2=mgr[n].loc.lon;
      t2=geo_recale(t2,t1,(double) 180.0);
/* *----------------------------------------------------------------------------
     to be revised for very high latitudes*/
      if(fabs(t2-t1) > 1.0)  continue;
      double p2=mgr[n].loc.lat;
      if(fabs(p2-p1) > 1.0)  continue;
      double d=geo_haversin_km(t1,p1,t2,p2);
      if((d<dmin) && (n!=m) && (n!=m+1) && (n!=m-1)) {
#pragma omp critical(mininum_distance)
        {
        neighbour=n;
        dmin=d;
        t_xover=0.5*(t1+t2);
        p_xover=0.5*(p1+p2);
        }
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
//    printf("tracks %d %d \n",mgr[m].track,mgr[n].track);
    }

/* *----------------------------------------------------------------------------
  create legend with all valid primary xovers */
  legends02->nz=12;
  legends02->pen0=0;
  legends02->pen1=0;

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

  for (m=0;m<nmgr;m++) {
    delete [] legends02->points[m].z;
    legends02->points[m].z=new float[legends02->nz];
    }

  out=fopen(report,"w");

  for (k=0;k<nwave;k++) {
    legends02->np=0;
    for (l=0;l<nxover;l++) {
      n1=xover[l].node[0];
      i1=mgr[n1].wave_index(wave[k]);
      n2=xover[l].node[1];
      i2=mgr[n2].wave_index(wave[k]);
      if((i1==-1) || (i2==-1)) continue;
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
      printf("%5d %5d lon= %6.2lf lat= %6.2lf depth= %7.1lf d=%6.2lf km: a1= %7.1f G1= %7.1f a2= %7.1f G2= %7.1f - e= %7.1lf\n",
             n1,n2,x,y,h,xover[l].d,a1,p1,a2,p2,d);
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
      }
    sprintf(lgdfile,"xover-%s.lgd",wave[k]);
    status=lgd_save(lgdfile, legends, nlegends,NULL,comments);
    }
  fclose(out);
  delete[] comments;


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
/* *----------------------------------------------------------------------------
    miserable patch for colocated xovers */
    if(d<0.1) {
      x=0.5*(x1[0]+x1[1]);
      y=0.5*(y1[0]+y1[1]);
      xover[l].node[2]=m1+1;
      xover[l].node[3]=m2+1;
      xover[l].t=x;
      xover[l].p=y;
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
        if((plg_secantpoint(x1,y1,x2,y2,&x,&y)==PLG_LINES_SECANT)) {
          xover[l].node[2]=n1;
          xover[l].node[3]=n2;
          xover[l].t=x;
          xover[l].p=y;
          printf("%5d %5d %5d lon= %8.4lf lon= %8.4lf lon= %8.4lf lat= %8.4lf lat= %8.4lf lat= %8.4lf\n",l,m1,n1,x1[0],x1[1],x,y1[0],y1[1],y);
          printf("%5d %5d %5d lon= %8.4lf lon= %8.4lf lon= %8.4lf lat= %8.4lf lat= %8.4lf lat= %8.4lf\n",l,m2,n2,x2[0],x2[1],x,y2[0],y2[1],y);
          }
        }
      }
    }

/* *----------------------------------------------------------------------------
  remove redundancy*/
  redundant=new int[nxover];
  for (l=0;l<nxover;l++) redundant[l]=0;
  nredundants=0;

  for (l=0;l<nxover;l++) {
    m1=xover[l].node[0];
    m2=xover[l].node[1];
    n1=xover[l].node[2];
    n2=xover[l].node[3];
    if(redundant[l]==1) continue;
    if(n1==-1) continue;
    if(n2==-1) continue;
    for (k=l+1;k<nxover;k++) {
      if((xover[k].node[0]==n1) && (xover[k].node[2]==m1)) {
        if(((xover[k].node[1]==m2) && (xover[k].node[3]==n2)) || ((xover[k].node[1]==n2) && (xover[k].node[3]==m2))) {
          redundant[k]=1;
          nredundants++;
          }
        }
      if((xover[k].node[0]==m1) && (xover[k].node[2]==n1)) {
        if(((xover[k].node[1]==m2) && (xover[k].node[3]==n2)) || ((xover[k].node[1]==n2) && (xover[k].node[3]==m2))) {
          redundant[k]=1;
          nredundants++;
          }
        }
      if(xover[k].id==-1) continue;
      if(xover[k].id==xover[l].id) {
        redundant[k]=1;
	nredundants++;
        }
      }
    }

/* *----------------------------------------------------------------------------
  create legend with all valid secondary xovers */
  out=fopen(report,"w");

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
/* *----------------------------------------------------------------------------
      increment legend cardinal*/
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
      a1=abs(c1);
      p1=atan2(c1.imag(),c1.real())*r2d;
      if(p1<0.0) p1+=360.0;
      a2=abs(c2);
      p2=atan2(c2.imag(),c2.real())*r2d;
      if(p2<0.0) p2+=360.0;
      printf("%5d %5d lon= %6.2lf lat= %6.2lf depth= %7.1lf d=%6.2lf km: a1= %7.2f G1= %7.1f a2= %7.2f G2= %7.1f - e= %7.1lf\n",
             n1,n2,x,y,h,xover[l].d,a1,p1,a2,p2,d);
      fprintf(out,"%5d %5d lon= %6.2lf lat= %6.2lf depth= %7.1lf d=%6.2lf km: a1= %7.2f G1= %7.1f a2= %7.2f G2= %7.1f - e= %7.1lf\n",
             n1,n2,x,y,h,xover[l].d,a1,p1,a2,p2,d);
      xover[l].c[0][k]=c1;
      xover[l].c[1][k]=c2;
      xover[l].error[k]=c2-c1;
      }
    }


  for (k=0;k<nwave;k++) {
    legends02->ID=k;
    legends02->np=0;
int dum;
    for (dum=0;dum<nxover;dum++) {
      l=dum;
      if(xover[l].node[2]==-1) continue;
      if(redundant[l]==1) continue;
      for(j=0;j<4;j++) {
        n=xover[l].node[j];
        indice[j]=mgr[n].wave_index(wave[k]);
        }
      if((indice[0]==-1) || (indice[1]==-1)|| (indice[2]==-1)|| (indice[3]==-1)) continue;
/* *----------------------------------------------------------------------------
      if error too large, find an alternative xover*/
      if (abs(xover[l].error[k])>.05) {
        ncyclemax=0;
        target=-1;
        for(ll = 0; ll < nxover; ll++) {
          if(xover[ll].id==xover[l].id) {
            ncycle=NINT( 0.5*(mgr[xover[ll].node[0]].duree+mgr[xover[ll].node[1]].duree));
            printf("xover %d error=%f ncycle=%d\n",ll,abs(xover[ll].error[k]),ncycle);
            if(ll==l) continue;
            if(ncycle>ncyclemax) {
              target=ll;
              ncyclemax=ncycle;
              }
            }
          }
        if(target!=-1) {
          printf("substitue %d (%f) with %d (%f)\n",l,abs(xover[l].error[k]),target,abs(xover[target].error[k]));
          l=target;
          if(xover[l].node[2]==-1) xover[l].node[2]=xover[l].node[0];
          if(xover[l].node[3]==-1) xover[l].node[3]=xover[l].node[1];
//          xover[l].error[k]=xover[target].error[k];
          }
        }
/* *----------------------------------------------------------------------------
      increment legend cardinal*/
      legends02->np++;
      n=legends02->np-1;
      legends02->points[n].t=xover[l].t;
      legends02->points[n].p=xover[l].p;
      c1=xover[l].c[0][k];
      c2=xover[l].c[1][k];
      d=abs(c2-c1);
      x=map_recale(topogrid,xover[l].t);
      y=xover[l].p;
      status=map_interpolation(topogrid,topo,topomask,x,y,&h);
/* *----------------------------------------------------------------------------
      update legend*/
      j=0;
      legends02->points[n].N=xover[l].id;
/* *----------------------------------------------------------------------------
      4 first records: amplitude,phase for ascending and descending tracks*/
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
/* *----------------------------------------------------------------------------
      5th record: vector difference*/
      legends02->points[n].z[j++]=d;
/* *----------------------------------------------------------------------------
      6th record: mean depth*/
      legends02->points[n].z[j++]=h;
/* *----------------------------------------------------------------------------
      7 to 10th record: number of valid observations for surrounding track points*/
      legends02->points[n].z[j++]=mgr[xover[l].node[0]].duree;
      legends02->points[n].z[j++]=mgr[xover[l].node[2]].duree;
      legends02->points[n].z[j++]=mgr[xover[l].node[1]].duree;
      legends02->points[n].z[j++]=mgr[xover[l].node[3]].duree;
/* *----------------------------------------------------------------------------
      11th record: vector difference in percent*/
      legends02->points[n].z[j++]=100.*abs(c2-c1)/fabs(a1+a2);
/* *----------------------------------------------------------------------------
      12th record: mean amplitude*/
      legends02->points[n].z[j++]=0.5*fabs(a1+a2);
      }
    sprintf(lgdfile,"xover-improved-%s.lgd",wave[k]);
    status=lgd_save(lgdfile, legends, nlegends,NULL,comments);
    }

  fclose(out);

  delete[] comments;

  for (k=0;k<nwave;k++) {
    free(a[k]);
    free(G[k]);
    }
  free(a);
  free(G);

  __ERR_BASE_LINE__("%s -computation complete ^^^^^^^^^\n",argv[0]);
  exit(0);
error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
}
