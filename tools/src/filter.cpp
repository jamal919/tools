
/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "tools-define.h"
#include "poc-time.h"
#include "filter.h"
#include "statistic.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int filter1D_BF_template(int mode, double L, double *x, int nvalues,T *buffer, T mask, T* & BF)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int a=4;
  
  switch(mode) {
    case LOESS:
//      status=Lanczos2D_weights(L, dx, dy, a, filter);
      break;
    case LANCZOS:
      status=Lanczos1D_BF(L, a, x, nvalues, buffer, mask, BF);
      break;
    }
    
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int filter1D_BF(int mode, double L, double *x, int nvalues, complex<float> *buffer, complex<float> mask, complex<float>* & BF)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=filter1D_BF_template(mode, L, x, nvalues, buffer, mask, BF);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int filter1D_BF_template(int mode, double *L, double *x, int nvalues,T *buffer, T mask, T* & BF)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int a=4;
  
  switch(mode) {
    case LOESS:
//      status=Lanczos2D_weights(L, dx, dy, a, filter);
      break;
    case LANCZOS:
      status=Lanczos1D_BF(L, a, x, nvalues, buffer, mask, BF);
      break;
    }
    
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int filter1D_BF(int mode, double *L, double *x, int nvalues, complex<float> *buffer, complex<float> mask, complex<float>* & BF)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=filter1D_BF_template(mode, L, x, nvalues, buffer, mask, BF);
  
  return(status);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int filter2D_InitWeigths(int mode, double scale, double dx, double dy, filter_t<float> & filter)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int a=4;
  double L;
  
  switch(mode) {
    case LOESS:
      status=Lanczos2D_weights(L, dx, dy, a, filter);
      break;
    case LANCZOS:
      status=Lanczos2D_weights(L, dx, dy, a, filter);
      break;
    }
    
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int filter2D_InitWeigths(int mode, double scale_x, double scale_y, double dx, double dy, int nxmax, int nymax, filter_t<float> & filter)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int a=4;
  double L;
  double azimuth;
  
  switch(mode) {
    case LOESS:
      status=Loess2D_init_cartesian(scale_x, scale_y, azimuth, dx, dy, nxmax, nymax, filter);
      break;
    case LANCZOS:
      status=Lanczos2D_weights(L, dx, dy, a, filter);
      break;
    }
    
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

statistic_t filtering(const char *save,const double *h,const double *t, double mask, int n, double dt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int i;
  double z;
  double *hf=NULL,*mhf=NULL,*mlf=NULL,*lf=NULL,*vlf=NULL,*vhf=NULL,*var=NULL,*uhf=NULL,*chf=NULL;
  double *topex=NULL,*ers=NULL,*env=NULL,*chk=NULL;
  double cut,time;
  bool add_line=True;
  FILE *out=NULL;
  statistic_t s,sraw,s_topex;

  exitIfNull(
    hf=new double[n]
    );
  exitIfNull(
    lf=new double[n]
    );
  exitIfNull(
    mhf=new double[n]
    );
  exitIfNull(
    mlf=new double[n]
    );
  exitIfNull(
    vhf=new double[n]
    );
  exitIfNull(
    uhf=new double[n]
    );
  exitIfNull(
    chf=new double[n]
    );
  exitIfNull(
    vlf=new double[n]
    );

  exitIfNull(
    var=new double[n]
    );

  exitIfNull(
    env=new double[n]
    );
  exitIfNull(
    topex=new double[n]
    );
  exitIfNull(
    ers=new double[n]
    );

  exitIfNull(
    chk=new double[n]
    );

  cut=30;
  Loess1D(n,dt,cut,h,mask,vlf,&status);

  cut=15;
  Loess1D(n,dt,cut,h,mask,lf,&status);

  cut=5;
  Loess1D(n,dt,cut,h,mask,mlf,&status);

  cut=2.5;
  Loess1D(n,dt,cut,h,mask,mhf,&status);

/* *----------------------------------------------------------------------------
  limit of significant period (as inferred from time sampling in meteo forcing)*/
  cut=0.25;
  Loess1D(n,dt,cut,h,mask,hf,&status);

/* *----------------------------------------------------------------------------
  limit of significant period (as inferred from time sampling in meteo forcing)*/
  cut=0.125;
//  cut=0.25;
  Loess1D(n,dt,cut,h,mask,uhf,&status);
  
  const double chf0=0.25;
  Loess1D(n,dt,chf0,h,mask,chf,&status);

  cut=10;
  Loess1D(n,dt,cut,h,mask,topex,&status);

  for (i=0; i< n; i++) {
    if(h[i] != mask) {
      vhf[i]=h[i]-hf[i];
      uhf[i]=h[i]-uhf[i];
      chf[i]=h[i]-chf[i];
      topex[i]=h[i]-topex[i]; /* 0 to 20 days */
      hf[i]=hf[i]-mhf[i];
      mhf[i]=mhf[i]-mlf[i];
      mlf[i]=mlf[i]-lf[i];
      lf[i]=lf[i]-vlf[i];
      }
    else {
      vhf[i]=mask;
      uhf[i]=mask;
      chf[i]=mask;
      hf[i]=mask;
      mhf[i]=mask;
      mlf[i]=mask;
      lf[i]=mask;
      topex[i]=mask;
      }
    }

 for (i=0; i< n; i++) {
    if(h[i] != mask) {
      var[i]=topex[i]*topex[i];
      }
    else {
      var[i]=mask;
      }
    }

  cut=5;
  Loess1D(n,dt,cut,var,mask,env,&status);

 for (i=0; i< n; i++) {
    if(env[i] != mask) {
      env[i]=sqrt(env[i]);
      }
    else {
      env[i]=mask;
      }
    }

/*printf("raw signal--------------\n");*/
  printf("  0<T<inf    -- ");
  sraw=get_statistics(h,mask,n);

/*printf("very low frequency------\n");*/
  printf(" 60<T<inf    -- ");
  s=get_statistics(vlf,mask,n);

/*printf("low frequency-----------\n");*/
  printf(" 30<T< 60    -- ");
  s=get_statistics(lf,mask,n);

/*printf("medium/low frequency----\n");*/
  printf(" 10<T< 30    -- ");
  s=get_statistics(mlf,mask,n);

/*printf("medium/high frequency---\n");*/
  printf("  5<T< 10    -- ");
  s=get_statistics(mhf,mask,n);

/*printf("high frequency----------\n");*/
  printf("0.5<T<  5    -- ");
  s=get_statistics(hf,mask,n);

/*printf("very high frequency-----\n");*/
  printf("  0<T<0.5    -- ");
  s=get_statistics(vhf,mask,n);

/*printf("ultra high frequency-----\n");*/
  printf("  0<T<0.25   -- ");
  s=get_statistics(uhf,mask,n);

/*printf("high frequency energy---\n");*/
  printf("  0<T<  5    -- ");
  s=get_statistics(env,mask,n);

  for (i=0; i< n; i++) {
    if(env[i] != mask) {
      if(abs(env[i]-s.mean) < 3*s.std) chk[i]=env[i];
      else chk[i]=mask;
      }
    else {
      chk[i]=mask;
      }
    }
/*printf("high frequency energy---\n");*/
  printf("  0<T<  5    -- ");
  s=get_statistics(chk,mask,n);
 
  for (i=0; i< n; i++) {
    if(env[i] != mask) {
      if(abs(chk[i]-s.mean) > 3*s.std) chk[i]=mask;
      }
    else {
      chk[i]=mask;
      }
    }
/*printf("high frequency enrrgy---\n");*/
  printf("  0<T<  5    -- ");
  s=get_statistics(chk,mask,n);

  if(strlen(save) != 0) {
    out=fopen(save,"w");
    if (out == NULL) TRAP_ERR_EXIT(errno,"fopen(\"%s\",\"w\") error (%d %s)\n",save,errno,strerror(errno));
    int nskipped=0;
    printf("saving %d samples to %s: ",n,save);
    for (i=0; i< n; i++) {
      if(h[i] != mask) {
        time=t[i];
        z=h[i]-vhf[i]-sraw.mean;
        fprintf(out,"%12.6f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
        time,h[i],vlf[i],lf[i],mlf[i],mhf[i],hf[i],vhf[i],mhf[i]+hf[i],z,topex[i],env[i],topex[i]-vhf[i]);
  /*    1    2    3      4      5      6     7     8      9           10  11      12 */
        add_line=True;
        }
      else {
/*         if(add_line) fprintf(out,"\n");  */
/*         time=t[i]; */
/*         fprintf(out,"%12.6f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n", */
/*         time,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.); */
        add_line=False;
        nskipped++;
        }
      }
    printf("%d skipped.\n",nskipped);
    fclose(out);
    }

  cut=10;
  Loess1D(n,dt,cut,h,mask,mlf,&status);

  cut=5;
  Loess1D(n,dt,cut,h,mask,mhf,&status);

  for (i=0; i< n; i++) {
    if(h[i] != mask) {
      mhf[i]=h[i]-mhf[i]-vhf[i];
      mlf[i]=h[i]-mlf[i]-vhf[i];
      }
    else {
      mhf[i]=mask;
      mlf[i]=mask;
      }
    }

  printf("0.5<T< 20    -- ");
  s_topex=get_statistics(mlf,mask,n);

  printf("0.5<T< 10    -- ");
  s=get_statistics(mhf,mask,n);

  cut=10;
  Loess1D(n,dt,cut,h,mask,mlf,&status);

#define REMOVE_CUSTOM_HF 0
  
  for (i=0; i< n; i++) {
    if(h[i] != mask) {
#if REMOVE_CUSTOM_HF
      mlf[i]=h[i]-mlf[i]-chf[i];
#else
      mlf[i]=h[i]-mlf[i]-uhf[i];
#endif
      }
    else {
      mlf[i]=mask;
      }
    }

#if REMOVE_CUSTOM_HF
  printf("%g<T<20  -- ",2*chf0);
#else
  printf("0.25<T< 20  -- ");
#endif
  s_topex=get_statistics(mlf,mask,n);


  delete[]uhf;
  delete[]vhf;
  delete[]chf;
  delete[]hf;
  delete[]mhf;
  delete[]vlf;
  delete[]lf;
  delete[]mlf;
  delete[]var;
  delete[]env;
  delete[]topex;
  delete[]ers;
  delete[]chk;

  return(s_topex);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void season(char *save, double *h, double *t, double mask, int n, double dt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,m;
  double *buffer=NULL;
  double time1,time2;
  date_t date;


  exitIfNull(
    buffer=new double[n]
    );

  date.day=1;
  date.month=5;
  date.year=1998;
  date.second=0;
  time1=cnes_time(date,'d');
  date.day=1;
  date.month=3;
  date.year=1999;
  date.second=0;
  time2=cnes_time(date,'d');

  m=-1;
  for (i=0;i<n;i++) {
    if((t[i] >= time1) && (t[i] <= time2)) {
      m++;
      buffer[m]=h[i];
      }
    }
  printf("reduced analysis: %lf to %lf, %d values \n",time1,time2,m);
  filtering(save,buffer,t,mask,m+1,dt);
  free(buffer);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fourier(const char *save, double *h, double mask, int n, double dt, int option)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,nf;
  int i;
  double *power=NULL,*phase=NULL,*frqcy=NULL,*vhf=NULL;
  double cut,df;
  FILE *out=NULL;
  statistic_t s,sraw;
//   char *axex=(char *)("d/h"),*axey=(char *)("cm"),*title=(char *)("NO TITLE"),*psfile=NULL;

  if(n==0) return;

  exitIfNull(
    power=new double[n]
    );
  exitIfNull(
    phase=new double[n]
    );
  exitIfNull(
    frqcy=new double[n]
    );

/* FHL 10/09/2002:
    psd_r1 (&dt,h,&mask,&n,power,phase,frqcy,&nf,&option,&status); */
  psd_r1_c (dt,h,mask,n,power,phase,frqcy,&nf,option&3,&status);
  exitIfNull(
    vhf=new double[n]
    );
  df=1./(n*dt);
  cut=1./60.;
  Loess1D(nf,df,cut,power,mask,vhf,&status);

  if(strlen(save) != 0) {
    double sum=0.;
    out=fopen(save,"w");
    if(out==0) TRAP_ERR_EXIT(-1,"Probleme a l ouverture du fichier : %s \n",save);
    if(option&4)
      fprintf(out,"#frequency(deg/sampling),power,phase(deg),filtered power,cumulated power\n");
    for (i=1; i< nf; i++) {
      sum+=power[i];
      fprintf(out,"%f %f %f %f %f\n",frqcy[i]*360.,power[i],phase[i],vhf[i],sum);
      }
    fclose(out);
    }
/*
  plot_regularpsd_1d (frqcy,power,&nf,save,psfile,axex,axey,title);
*/
  free(power);
  free(phase);
  free(vhf);
  free(frqcy);
}

