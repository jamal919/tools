
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Detides comodo-compliant NetCDF outputs and produces tidal atlases.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#include "tools-structures.h"

#include "tools-define.h"

#include "statistic.h"
#include "geo.h"

#include "altimetry.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int loess1d_variable(int n, double l_c, double *tuning, T *h, double *x, T mask,T *lf, float *cut)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j;
  int holes=0;
  double s,q,*w=NULL;
  double width;

  if(n==1) lf[0]=h[0];/* cut is NaN there! */
  else {
    exitIfNull(
      w=(double *) malloc(n*sizeof(double))
      );
//    l_c/=2;
    for(i=0;i<n;i++) {
/* *----------------------------------------------------------------------------
      L=T*sqrt(gh)~12.5*3600*3.3/1000 sqrt(h)~140 sqrt(h) w=*/
//      width=MAX(l_c/sqrt(2000./fabs(tuning[i])),30.0);
//      width=MAX(sqrt(fabs(tuning[i]))*l_c/50.,30.0);

      width=MAX(sqrt(fabs(tuning[i]))*5.0*l_c,50.0);
//      width=MAX(sqrt(fabs(tuning[i]))*5.0*l_c,25.0);

      width=MIN(width,120.0);
      cut[i]=width;
      for(j=0;j<n;j++) {
        q=(x[j]-x[i])/width;
//        q*=q;
        q=abs(q);
        s=1-q*q*q;
        if(q>1) w[j]=0;
        else w[j]=s*s*s;
        }
      s=0;
      lf[i]=0;
      for(j=0;j<n;j++) {
        if(h[j]!=mask) {
          s+=w[j];
          lf[i]+=((float) w[j])*h[j];
          }
        }
      if(s!=0) {
        lf[i]/=s;
        }
      else {
        lf[i]==mask;
        holes++;
        }
      }
    free(w);
    }
  return(holes);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int print_diag(complex<float> *h, double *t, double *duration, double *depth, double *lon, double *lat,
               complex<float> mask,float **noise, complex<float> *lwf, complex<float> *hgf,
                float *cut, int n, int track,int segment, char *wave)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k,count,status;
  int i,j,holes;
  FILE *out=NULL;
  char filename[1024];

  sprintf(filename,"%s.track-%3.3d-%2.2d.gnu",wave,track,segment);

  if ((out = fopen(filename, "w")) == NULL) {
      __ERR_BASE_LINE__("file opening issue : %s \n",filename);
      exit(-1);
    }

  for (i=0; i< n; i++) {
    if(hgf[i] != mask) {
/*  1->  5*/      fprintf(out,"%lf %lf %lf %lf %f ",      t[i],lon[i],lat[i],depth[i],cut[i]);
/*  6->  9*/      fprintf(out,"%lf %lf %lf %lf ",         abs(h[i]),abs(hgf[i]),abs(lwf[i]),abs(h[i])-abs(lwf[i]));
/* 10-> 12*/      fprintf(out,"%lf %lf %lf ",             arg(h[i]),arg(hgf[i]),arg(lwf[i]));
/* 13-> 18*/      fprintf(out,"%lf %lf %lf %lf %lf %lf ", h[i].real(), h[i].imag(),lwf[i].real(),lwf[i].imag(), hgf[i].real(),hgf[i].imag());
/* 19-> 21*/      fprintf(out,"%lf %lf %lf \n",           noise[0][i], noise[1][i], noise[2][i]);
      }
    }
  fclose(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int *check_outlayers(complex<float> *h, double *t, double *duration, double *depth, double *lon, double *lat,
                     complex<float> mask,float *noise, float **misfit, complex<float> **lwf, complex<float> **hgf,
                     float **cut, int n, int track,int segment, char *wave, double pulsation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k,count,status,noutlayers;
  int *outlayers;
  int i,j,holes;
  double mean,variance,z;
  complex<float> *hf,*lf, *hh;
  float delta;
  double time,limit;
  bool add_line=True;
  FILE *out=NULL;
  statistic_c_t shf;
  statistic_t   s;
  double std,L0;
  char filename[1024];

  hh =new complex<float>[n];
  hf =new complex<float>[n];
  lf =new complex<float>[n];

  *cut=new float[n];
  s=get_statistics(duration,(double) -1.,n);

/*-----------------------------------------------------------------------------
  copy original vector*/
  for (i=0; i< n; i++) {
/* *-----------------------------------------------------------------------------
    mask short records*/
    limit=MIN(s.mean-2*s.std,150.0);
    limit=MAX(limit,100.0);
    if(duration[i] > limit) {
      hh[i]=h[i];
      }
    else {
      hh[i]=mask;
      }
    }

/* *-----------------------------------------------------------------------------
  make a first selection from full vector std*/
//  cut=60;  // km
//  holes=loess1d (n,cut,hh,t,mask,lf);
//  L0=120.;
  L0=30.0/pulsation;
  holes=loess1d_variable (n,L0,depth,hh,t,mask,lf,*cut);

  for (i=0; i< n; i++) {
    if(h[i] != mask) {
      hf[i]=h[i]-lf[i];
      }
    else {
      hf[i]=mask;
      }
    }

  shf=cget_statistics(hf,mask,n);

  std=MAX(shf.std,0.02);

  for (i=0; i< n; i++) {
    if(hf[i] != mask) {
      delta=abs(hf[i]);
      if(delta>3.*std) {
        hh[i]=mask;
        }
      }
    }

/* *-----------------------------------------------------------------------------
  compute HF std deviation from selected (valid) values*/
//  cut=60;  // km
//  holes=loess1d_variable (n,cut,hh,t,mask,lf);
//  L0=120.;
  L0=30.0/pulsation;
  holes=loess1d_variable (n,L0,depth,hh,t,mask,lf,*cut);

  for (i=0; i< n; i++) {
    if(hh[i] != mask) {
      hf[i]=h[i]-lf[i];
      }
    else {
      hf[i]=mask;
      }
    }

  shf=cget_statistics(hf,mask,n);
  std=MAX(shf.std,0.005);

  for (i=0; i< n; i++) {
    if(h[i] != mask) {
      hf[i]=h[i]-lf[i];
      }
    else {
      hf[i]=mask;
      }
    }

/* *-----------------------------------------------------------------------------
  make the final selection*/
  noutlayers=0;
  for (i=0; i< n; i++) {
    if(hf[i] != mask) {
      delta=abs(hf[i]);
      if(delta>3.*std) {
//        printf("outlayer\n");
        noutlayers++;
        }
      }
    else {
      noutlayers++;
      }
    }

  outlayers=new int[noutlayers+1];
  *misfit=new float[noutlayers+1];

  noutlayers=0;
  for (i=0; i< n; i++) {
    if(hf[i] != mask) {
      delta=abs(hf[i]);
      noise[i]=delta;
      if(delta>3.*std) {
//        printf("outlayer\n");
        outlayers[noutlayers]=i;
        (*misfit)[noutlayers]=delta;
        noutlayers++;
        }
      }
    else {
      outlayers[noutlayers]=i;
      (*misfit)[noutlayers]=1.0;
      noutlayers++;
      }
    }

  outlayers[noutlayers]=-1;

//   sprintf(filename,"%s.track-%3.3d-%2.2d.dat",wave,track,segment);
//   out=fopen(filename,"w");
//   for (i=0; i< n; i++) {
//     if(hf[i] != mask) {
//       fprintf(out,"%lf %lf %lf %lf %f ",      t[i],lon[i],lat[i],depth[i],(*cut)[i]);
//       fprintf(out,"%lf %lf %lf %lf ",         abs(h[i]),abs(hf[i]),abs(lf[i]),abs(h[i])-abs(lf[i]));
//       fprintf(out,"%lf %lf %lf ",             arg(h[i]),arg(hf[i]),arg(lf[i]));
//       fprintf(out,"%lf %lf %lf %lf %lf %lf ", h[i].real(), h[i].imag(),lf[i].real(),lf[i].imag(), hf[i].real(),hf[i].imag());
//       fprintf(out,"%lf \n", noise[i]);
//       }
//     }
//   fclose(out);

  *lwf=lf;
  *hgf=hf;

  delete[] hh;
//  delete[] hf;
//  delete[] lf;

  return(outlayers);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int compute_coastal_error(char *wave, int kk, vector<mgr_t> mgr,int nmgr, xerror_t xerror, double *buble, int **list, int *card)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int *track,*subtrack,ntracks,ntics;
  int count,npt,i,j,m,mm,n,nn,k,l,ll,w,status;
  int *start,*ends;
  int n1,n2,n3;
  int *outlayers;
  float *misfit;
  double *x,*y,*h,*d, *t,*lon, *lat,tt,pulsation;
  float *noise[10], *cut, error_min, rmask=data_rmask;
  std::complex<float> cmask=-1;
  std::complex<float> *c,*lf,*hf;

  subtrack = new int   [nmgr];
  track    = new int   [nmgr];
  x     = new double[nmgr];
  y     = new double[nmgr];

/* *----------------------------------------------------------------------------
  count number of different tracks and register them*/
  ntracks=0;
  for(n = 0; n < nmgr; n++) {
    track[n]=mgr[n].track;
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
  for(n = 0; n < nmgr; n++) {
    x[n]=geo_recale(mgr[n].loc.lon,0.,180.);
    y[n]=mgr[n].loc.lat;
    xerror.alongtrack[n]=0.;
    }

/* *----------------------------------------------------------------------------
  examin along-track errors*/
  for (k=0;k<ntracks;k++) {
/* *----------------------------------------------------------------------------
    detect holes in track sequence*/
/* *----------------------------------------------------------------------------
    first count*/
    count=1;
    for(nn = 1; nn < card[k]; nn++) {
      m=list[k][nn-1];
      n=list[k][nn];
      tt=geo_km(x[m],y[m],x[n],y[n]);
//      printf("%d %d -> %d : %lf %lf %lf \n",k,m,n,x[n],y[n],tt);
      if(tt>10.) count++;
      }
/* *----------------------------------------------------------------------------
    then build arrays*/
    start=new int[count];
    ends=new int[count];
    start[0]=0;
    ends[0]=card[k]-1;
    count=1;
    for(nn = 1; nn < card[k]; nn++) {
      m=list[k][nn-1];
      n=list[k][nn];
      tt=geo_km(x[m],y[m],x[n],y[n]);
      if(tt>12.) {
        ends[count]=ends[count-1];
        ends[count-1]=nn-1;
        start[count] =nn;
        count++;
        }
      }
/* *----------------------------------------------------------------------------
    check data consistency by track's sub-segments*/
    for(l=0;l<count;l++) {
      npt=ends[l]-start[l]+1;
      d=new double[npt];
      h=new double[npt];
      t=new double[npt];
      lon=new double[npt];
      lat=new double[npt];
      noise[0]=new float[npt];
      noise[1]=new float[npt];
      noise[2]=new float[npt];
      c=new complex<float>[npt];
      m=0;
      for(nn = start[l]; nn <= ends[l]; nn++) {
        n=list[k][nn];
        w=mgr[n].wave_index(wave);
        h[m]=mgr[n].loc.depth;
        lon[m]=mgr[n].loc.lon;
        lat[m]=mgr[n].loc.lat;
        d[m]=mgr[n].duree;
        t[m]=geo_km(x[list[k][0]],y[list[k][0]],x[n],y[n]);
        if(w==-1) {
          c[m]=cmask;
          }
        else {
          c[m]=fct_polar2complex(mgr[n].data[w].amp,mgr[n].data[w].phi);
          pulsation=mgr[n].data[w].constituent.omega;
          }
        m++;
        }
      outlayers= check_outlayers(c,t,d,h,lon,lat,cmask,noise[0],&misfit,&lf,&hf,&cut,m,track[n],l,wave,pulsation);
/* *----------------------------------------------------------------------------
      background noise (internal tide signature, etc..)*/
      m=0;
      for(nn = start[l]; nn <= ends[l]; nn++) {
        n=list[k][nn];
        xerror.alongtrack[n]=noise[0][m];
        noise[1][m]=xerror.harmonic[n];
        noise[2][m]=xerror.alias[n];
        xerror.cut[n]=cut[m];
        xerror.zraw[n]=c[m];
        xerror.zlwf[n]=lf[m];
        xerror.zhgf[n]=hf[m];
/* *----------------------------------------------------------------------------
        covariance radius (km)*/
        buble[n]=25.0;
        m++;
        }
      status= print_diag(c,t,d,h,lon,lat,cmask,noise,lf,hf,cut,m,track[n],l,wave);
      nn=0;
/* *----------------------------------------------------------------------------
      special for outlayers*/
      while(outlayers[nn]!=-1) {
        n=list[k][outlayers[nn]+start[l]];
//        w=mgr[n].wave_index(wave);
        w=mgr[n].wave_index(wave);
        if(w==-1) {
          printf("track %d (%d), outlayer %3.3d gauge=%3.3d : %lf %lf masked \n",k,track[n],outlayers[nn],n,x[n],y[n]);
          xerror.alongtrack[n]=data_rmask;
          }
        else {
          printf("track %d (%d), outlayer %3.3d gauge=%3.3d : %9.3lf %9.3lf %9.3f %9.3f %9.3f \n",
                 k,track[n],outlayers[nn],n,x[n],y[n],mgr[n].data[w].amp,mgr[n].data[w].phi,misfit[nn]);
          xerror.alongtrack[n]=2.*misfit[nn];
          }
/* *----------------------------------------------------------------------------
        covariance radius (km)*/
        buble[n]=5.0;
        nn++;
        }
      delete[] outlayers;
      delete[] misfit;
      delete[] cut;
      delete[] lf;
      delete[] hf;
      delete[] h;
      delete[] t;
      delete[] c;
      delete[] d;
      delete[] lon;
      delete[] lat;
      delete[] noise[0];
      delete[] noise[1];
      delete[] noise[2];
      }
    delete[] start;
    delete[] ends;
    }

  delete[] x;
  delete[] y;
  delete[] subtrack;
  delete[] track;

  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int compute_analysis_error(char *wave, vector<mgr_t> mgr,int nmgr,float *analysis_error)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,m,n,k,l,w,status;

/* *----------------------------------------------------------------------------
  initialisation*/
  for(n = 0; n < nmgr; n++) {
    analysis_error[n]=0.;
    w=mgr[n].wave_index(wave);
    if(w!=-1) {
      analysis_error[n]=mgr[n].data[w].error;
      }
    }


  return(0);
  }
