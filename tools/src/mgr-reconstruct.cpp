
/*******************************************************************************
 
  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
/**
\file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada
\author  Clement MAYET      LEGOS, Toulouse, France (PhD)
\author  Yves Soufflet      LEGOS, Toulouse, France

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Convert between different sealevel time serie formats.
*/
/*----------------------------------------------------------------------------*/

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <map>

#include "tools-structures.h"

#include "tides.h"
#include "fe.h"
#include "archive.h"
#include "poc-netcdf.hpp"
#include "poc-time.h"
#include "mgr-converter.h"
#include "mgr.h"
#include "functions.h"
#include "spectrum.h"
#include "filter.h"

using namespace std;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LRC_analytic(double *u, double *t, int n, double sampling, double R, double C, double L)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  u=exp(-lamda t) cos(omega t - G)
  
  lamda=R/(2L)
  
  omega0²=1/LC
  
  omega=sqrt(omega0²-lamda²)
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double time;
  size_t size=n;
  double *LRC=new double[n];
  double lamda, omega, t0;
  int max_positive, max_negative;
  FILE *out;
  double Tauopt,Topt,t0opt;
  int posopt;
  
  double T  =20.*60;
  double Tau=62.*60;
  
  double maxcor=-1.e+10;
  
  lamda=1./Tau;
  omega=2*M_PI/T;
  
  for (T=5;T<25;T+=1)
  for (Tau=10;Tau<30;Tau+=2)
  for(t0=0; t0<300.0; t0+=10.0) {
    lamda=1./Tau/60.;
    omega=2*M_PI/T/60.;
    for(int k=0; k<n; k++) {
      time=(t[k]-t[0])*24*3600.-t0;
      if(time<0) LRC[k]=0.0;
      else LRC[k]=exp(-lamda*time)*cos(omega*time);
      }
    
    double *correlations=get_maxcorrelation(LRC, u, 600, 300, 1000, max_positive, max_negative, 0);
  
    int pos=maxpos(correlations, 450);
    
    if(correlations[pos]> maxcor){
      printf("t0=%lf T=%lf lamda=%lf maxcor=%lf\n", t0, T, Tau, correlations[pos]);
      Tauopt=Tau;
      Topt=T;
      t0opt=t0;
      maxcor=correlations[pos];
      posopt=pos;
      }
    
    delete[] correlations;
    }
    
  lamda=1./Tauopt/60.;
  omega=2*M_PI/Topt/60.;
  for(int k=0; k<n; k++) LRC[k]=0.0;
  for(int k=posopt; k<n; k++) {
    time=(t[k-posopt]-t[0])*24*3600.-t0opt;
    if(time<0) LRC[k]=0.0;
    else LRC[k]=exp(-lamda*time)*cos(omega*time);
    }
    
  out=fopen("LRC.gnu","w");
  for(int k=0;k<300;k++) {
    fprintf(out,"%lf %lf %lf\n", t[k], LRC[k], 300000.*u[k]);
    }
 
  fclose(out);
    
  
  return(0);
  
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void LRC_integrator(double *u, double *v, double *E, double sampling, double R, double C, double L)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  Ldi/dt+Ri=Ur
  
  i=dQ/dt, Uc=CQ
  
 
  LC d²u/dt² + RC du/dt + u = f
  
  LC/dt²+ RC/2/dt
  
  http://www.sciences.univ-nantes.fr/sites/jacques_charrier/tp/rcrlrlc/index.html
  
  https://fr.wikipedia.org/wiki/Circuit_RC
  
  https://fr.wikiversity.org/wiki/%C3%89tude_des_syst%C3%A8mes_%C3%A9lectriques/Condensateur_et_circuit_RC
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double *f, *s, value, T;
  double dt;
  size_t size;
  double count;
  FILE* out;
  
  size=10000000;
  
//   R=1;
//   C=5.0;
  
  T=5.0;
  
  dt=0.01;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  integrate signal
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(u==0) u=new double[size];
  
  u[0]=0.0;
  u[1]=0.0;
  
  for(int k=0;k<size;k++) {
    double time=k*dt;
    int m=time/sampling;
    double a= L*C/dt/dt+R*C/2.0/dt;
    double b=-L*C/dt/dt+R*C/2.0/dt;
    double c=2*L*C/dt/dt;
    u[k+1]=u[k-1]*b/a+(E[m]-u[k]+c*u[k])/a;
    u[k]=0.9*u[k]+0.05*(u[k-1]+u[k+1]);
    }
    
  out=fopen("LRC-integrator.gnu","w");
  for(int k=0;k<size;k++) {
    double time=k*dt;
    if(fmod(time,sampling)==0.0) {
      int m=time/sampling;
      fprintf(out,"%lf %lf %lf %lf\n", time/3600., u[k], v[m], E[m]);
      }
    }
 
  fclose(out);
  
  printf("end of integration\n");
  
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void RC_integrator(double *u, double *v, double *E, double sampling, double R, double C)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  RC du/dt + u = f
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double *f, *s, value, T;
  double dt;
  size_t size;
  double count;
  FILE* out;
  
  size=10000000;
  
//   R=1;
//   C=5.0;
  
  T=5.0;
  
  dt=0.1;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  integrate signal
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(u==0) u=new double[size];
  
  u[0]=0.0;
  u[1]=0.0;
  
  for(int k=0;k<size;k++) {
    double time=k*dt;
    int m=time/sampling;
    u[k+1]=u[k-1]+2.0*dt*(E[m]-u[k])/R/C;
    u[k]=0.9*u[k]+0.05*(u[k-1]+u[k+1]);
    }
    
  out=fopen("RC-integrator.gnu","w");
  for(int k=0;k<size;k++) {
    double time=k*dt;
    if(fmod(time,sampling)==0.0) {
      int m=time/sampling;
      fprintf(out,"%lf %lf %lf %lf\n", time/3600., u[k], v[m], E[m]);
      }
    }
 
  fclose(out);
  
  printf("end of integration\n");
  
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int seek_bestRC(double *u, double *E, double dt, int n, double & RC)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  RC du/dt + u = f
  
  J=sum[(f-(RC du/dt + u))²]
  
  dJ/dRC=-2 sum[(f-(RC du/dt + u))*du/dt]
  
  dJ/dRC=0 -> RC sum[(du/dt)²]=sum[(f-u)*du/dt]
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double sum,a;
  double *d=new double[n];
  double R, C;
  double *predicted=0;
  
  for(int k=1;k<n-1;k++) {
    d[k]=(u[k+1]-u[k-1])/2./dt;
    }
    
  a=0;
  sum=0;
  
  for(int k=1;k<n-1;k++) {
    a+=d[k]*d[k];
    sum+=(E[k]-u[k])*d[k];
    }
  
  RC=sum/a;
  
  R=1;
  C=1000*RC;
  
  RC_integrator(predicted, u, E, dt, R, C);
  
  C=1000*RC;
  LRC_integrator(predicted, u, E, dt, R, C, 5000.0);

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int seek_bestLRC(double *u, double *E, double dt, int n, double & RC, double & LC)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  LC d²u/dt² + RC du/dt + u = f
  
  J=sum[(f-(LC d²u/dt² + RC du/dt + u))²]
  
  
  dJ/dRC=-2 sum[(f-(LC d²u/dt² + RC du/dt + u))*du/dt]
  
  dJ/dRC=0 -> LC sum[d²u/dt² du/dt]+RC sum[(du/dt)²]=sum[(f-u)*du/dt]
  
 
  dJ/dLC=-2 sum[(f-(LC d²u/dt² + RC du/dt + u))*d²u/dt²]
  
  dJ/dLC=0 -> LC sum[(d²u/dt²)²]+RC sum[d²u/dt² du/dt]=sum[(f-u)*d²u/dt²]
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  double sum1, sum2,a,b,c,d;
  double *d1=new double[n], *d2=new double[n];
  double *matrix=new double[4], *inverse;
  double R, C;
  double *predicted=0;
  
  for(int k=1;k<n-1;k++) {
    d1[k]=(u[k+1]-u[k-1])/2./dt;
    d2[k]=(u[k+1]-2.*u[k-1]+u[k-1])/dt/dt;
    }
    
  a=0;
  b=0;
  c=0;
  d=0;
  sum1=0;
  sum2=0;
  
  for(int k=1;k<n-1;k++) {
    a+=d1[k]*d1[k];
    b+=d1[k]*d2[k];
    c+=d1[k]*d2[k];
    d+=d2[k]*d2[k];
    sum1+=(E[k]-u[k])*d1[k];
    sum2+=(E[k]-u[k])*d2[k];
    }
  
  matrix[0]=a;
  matrix[1]=b;
  matrix[2]=c;
  matrix[3]=d;
  
//   status=matrix_inverse(matrix, 2, inverse, false, 1);
    
  R=1;
  C=1000*RC;
  
  RC_integrator(predicted, u, E, dt, R, C);
  
  C=1000*RC;
  LRC_integrator(predicted, u, E, dt, R, C, 5000.0);

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void RC_simulator(void)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  RCdu/dt+u=f
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double *u, *f, *s, *E, dt, value, T;
  double R, C, L;
  size_t size;
  double count;
  FILE* out;
  
  size=10000;
  
  R=1;
  C=1.0;
  
  dt=0.01;
  T=5.0;
  
  u=new double[size+1];
  f=new double[size];
  s=new double[size];
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  generate step-wise forcing
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  count=0;
  value=0.0;
  for(int k=0;k<size;k++) {
    f[k]=value;
    s[k]=value;
    count+=dt;
    if(count>T) {
      value=1.0-value;
      count=0;
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  smooth step-wise forcing
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  int half_window=10;
  
  for(int k=half_window;k<size-half_window;k++) {
    double w=0.0;
    value=0.0;
    for(int l=-half_window;l<half_window;l++) {
      value+=f[k+l];
      w+=1.0;
      }
    s[k]=value/w;
    f[k]=value/w;
    }
  
  for(int k=half_window;k<size-half_window;k++) {
    double w=0.0;
    value=0.0;
    for(int l=-half_window;l<half_window;l++) {
      value+=f[k+l];
      w+=1.0;
      }
    s[k]=value/w;
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  generate step-wise forcing
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  count=0;
  value=1.0;
  for(int k=0;k<size;k++) {
    double r=(T-count)/T;
    double v=value*r*r-0.5;
    f[k]=v;
    s[k]=v;
    count+=dt;
    if(count>T) {
      value=1.0;
      count=0;
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  integrate signal
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  E=s;
  
  u[0]=0.0;
  u[1]=0.0;
  
  for(int k=0;k<size;k++) {
    u[k+1]=u[k-1]+2.0*dt*(E[k]-u[k])/R/C;
    u[k]=0.9*u[k]+0.05*(u[k-1]+u[k+1]);
    }
    
  out=fopen("RC.gnu","w");
  for(int k=0;k<size;k++) {
    fprintf(out,"%lf %lf %lf\n", k*dt, u[k], E[k]);
    }
 
  fclose(out);
  
  printf("end of integration\n");
}
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void RCQ_simulator(void)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  RCdu/dt+u²=f
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double *u, *f, *s, *E, dt, value, offset;
  double T1,T2,T;
  double amplitude;
  double R, C, L;
  size_t size;
  double count;
  FILE* out;
  
  size=100000;
  
  R=4.0;
  C=1.0;
  
  amplitude=1.25;
  
  dt=0.001;
  T=5.0;
  T1=2.85,T2=T-T1;
  offset=0.0;
  
  u=new double[size+1];
  f=new double[size];
  s=new double[size];
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  generate step-wise forcing
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  count=0;
  value=0.0;
  int flag=0;
  for(int k=0;k<size;k++) {
    f[k]=value-offset;
    s[k]=value-offset;
    count+=dt;
    switch (flag) {
      case 0:
        if(count>T1) {
          value=amplitude-value;
          count=0;
          flag=1-flag;
          }
        break;
      case 1:
        if(count>T2) {
          value=amplitude-value;
          count=0;
          flag=1-flag;
          }
        break;
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  smooth step-wise forcing
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  int half_window=500;
  
  for(int k=half_window;k<size-half_window;k++) {
    double w=0.0,sum=0.0;
    value=0.0;
    if(f[k]>f[k-half_window]) continue;
//     for(int l=-half_window;l<half_window+1;l++) {
    for(int l=0;l<half_window+1;l++) {
      double r,q=1.- abs(double (l))/(double) half_window;
      q=abs(q);
      r=(1-q*q*q);
      w=r*r*r;
      value+=f[k+l]*w;
      sum+=w;
      }
    if(s[k]==amplitude) s[k]=value/sum;
//     f[k]=value/w;
    }
//   
//   for(int k=half_window;k<size-half_window;k++) {
//     double w=0.0;
//     value=0.0;
//     for(int l=-half_window;l<half_window;l++) {
//       value+=f[k+l];
//       w+=1.0;
//       }
//     s[k]=value/w;
//     }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  generate step-wise forcing
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   count=0;
//   value=1.0;
//   for(int k=0;k<size;k++) {
//     double r=(T-count)/T;
//     double v=value*r*r-0.5;
//     f[k]=v;
//     s[k]=v;
//     count+=dt;
//     if(count>T) {
//       value=1.0;
//       count=0;
//       }
//     }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  integrate signal
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  E=s;
  
  u[0]=0.0;
  u[1]=0.0;
  
  for(int k=0;k<size;k++) {
    u[k+1]=u[k-1]+2.0*dt*(E[k]-u[k])/R/C/abs(u[k]+.05);
    u[k]=0.9*u[k]+0.05*(u[k-1]+u[k+1]);
    }
    
  out=fopen("RCQ.gnu","w");
  for(int k=0;k<size;k++) {
    fprintf(out,"%lf %lf %lf %lf\n", k*dt, u[k], E[k], f[k]);
    }
 
  fclose(out);
  
  printf("end of integration\n");
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void RC_old(void)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  RCdu/dt+u=f
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double *u, *f, *s, *E, dt, value, T;
  double R, C, L;
  size_t size;
  double count;
  FILE* out;
  
  tseries_t *serieA, *serieB;
  string filename;
  int ncolumns=1, status;
  const char *header_format=0;
  mooring_t mooring;
  char unit='m';
  const char *format=FORMAT_NAME_GNU;
  
  filename="HE_Rouen_20160601-20170731_20160531_20170731.gnu";
  mgr_load_timeserie(filename.c_str(), &serieA, unit, format, ncolumns, header_format, &mooring, &status);
          
  filename="HE_Honfleur_20160601-20170731_20160531_20170731.gnu";
  mgr_load_timeserie(filename.c_str(), &serieB, unit, format, ncolumns, header_format, &mooring, &status);
  
  int max_positive, max_negative;
  double *correlations=get_maxcorrelation(serieA[0].x[0],serieB[0].x[0], 100000, 1000, 1000, max_positive, max_negative);
  
  double RC;
  status=seek_bestRC(serieA[0].x[0], &(serieB[0].x[0][max_positive]), 300.0, 100000, RC);
  
  double *i1=new double[serieA[0].n];
  double *pred=new double[serieA[0].n];
  double *i1HF=new double[serieA[0].n];
  double *i1BF=new double[serieA[0].n];
  double *d1=new double[serieA[0].n];
  double *d2=new double[serieA[0].n];
  
  dt=serieA[0].t[1]-serieA[0].t[0];
  dt*=3600.*24.;
  
  C=5.e+03;
  L=100.0;
  R=2.0;
  
  statistic_t stat;
  stat=get_statistics(serieA[0].x[0],serieA[0].n);
  
  i1[0]=0;
  for(int k=1;k<serieA[0].n-1;k++) {
    d1[k]=-(serieA[0].x[0][k+1]-serieA[0].x[0][k-1])*R*C/2./dt+serieA[0].x[0][k];
    d2[k]=d1[k]+(serieA[0].x[0][k+1]-2.*serieA[0].x[0][k-1]+serieA[0].x[0][k-1])*L*C/dt/dt;
    
    i1[k]=i1[k-1]+(serieA[0].x[0][k]-stat.mean)*dt/(365.*24.*60.);
    d1[k]=(serieA[0].x[0][k+1]-serieA[0].x[0][k-1])/2./dt;
    d2[k]=(serieA[0].x[0][k+1]-2.*serieA[0].x[0][k]+serieA[0].x[0][k-1])/dt/dt;
    }
  double sampling;
//   status=LRC_analytic(d2, serieA[0].t, serieA[0].n, sampling, R, C, L);
  
  status=LRC_analytic(&(d2[10000]), &(serieA[0].t[10000]), serieA[0].n, sampling, R, C, L);

  
  double mean, window=1.5;
  double *residuals=new double[serieA[0].n];
  int nodal=1;
  bool Z0=false;
/*------------------------------------------------------------------------------
  long period needed because of non-linearities (Mf=M2-N2, MSf=M2-S2, etc...) */
  spectrum_t solved, spectrum = spectrum_init_ref("ESTUARINE", Z0);
  spectrum.remove(wSa);
  spectrum.remove(wSsa);
  
  status=loess1d_irregular(serieA[0].n, window, i1, serieA[0].t, serieA[0].mask, i1BF);
    
  for(int k=1;k<serieA[0].n-1;k++) {
    i1HF[k]=i1[k]-i1BF[k];
    }
       
  mgr_data_t *tmp= harmonic_analysis_with_parsing(i1HF,serieA[0].mask,serieA[0].t,residuals, &mean, serieA[0].n, spectrum, solved, nodal, stdout);
 
  for(int k=1;k<serieA[0].n-1;k++) {
    pred[k]=10000.*((i1HF[k+1]-residuals[k+1])-(i1HF[k-1]-residuals[k-1]))/2./dt;
    }
    
  out=fopen("RC-derived-A.gnu","w");
  for(int k=1;k<serieA[0].n-1;k++) {
    fprintf(out,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", serieA[0].t[k], serieA[0].x[0][k], 1000.*d1[k], 1000.*300.*d2[k], i1[k], i1HF[k], i1BF[k], residuals[k], i1HF[k]-residuals[k], pred[k]);
    }
    
  fclose(out);
  
  for(int k=1;k<serieB[0].n-1;k++) {
    d1[k]=-(serieB[0].x[0][k+1]-serieB[0].x[0][k-1])*R*C/2./dt+serieB[0].x[0][k];
    d2[k]=d1[k]+(serieB[0].x[0][k+1]-2.*serieB[0].x[0][k-1]+serieB[0].x[0][k-1])*L*C/dt/dt;
    d1[k]=(serieB[0].x[0][k+1]-serieB[0].x[0][k-1])/2./dt;
    d2[k]=(serieB[0].x[0][k+1]-2.*serieB[0].x[0][k]+serieB[0].x[0][k-1])/dt/dt;
    }
    
  out=fopen("RC-derived-B.gnu","w");
  for(int k=1;k<serieB[0].n-1;k++) {
    fprintf(out,"%lf %lf %lf %lf\n", serieB[0].t[k], serieB[0].x[0][k], 1000.*d1[k], 1000.*300.*d2[k]);
    }
 
  fclose(out);
  
  printf("end of integration\n");
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *derivative_1(tseries_t & serie, double sampling, int target, bool normalize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  statistic_t s;
  
  double *derivative=new double[serie.n];
  for(int k=0;k<serie.n;k++)  derivative[k]=serie.mask;
  
  for(int k=1;k<serie.n-1;k++) {
    if(serie.x[target][k+1]==serie.mask or serie.x[target][k-1]==serie.mask) continue;
    derivative[k]=(serie.x[target][k+1]-serie.x[target][k-1])/2./sampling;
    }
    
//   if(normalize) {
//     s=get_statistics(derivative, serie.mask, serie.n, 0);
//     }
  return(derivative);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *derivative_2(tseries_t & serie, double sampling, int target, bool normalize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  statistic_t s;
  
  double *derivative=new double[serie.n];
  for(int k=0;k<serie.n;k++)  derivative[k]=serie.mask;
  
  for(int k=1;k<serie.n-1;k++) {
    if(serie.x[target][k+1]==serie.mask or serie.x[target][k]==serie.mask or serie.x[target][k-1]==serie.mask) continue;
    derivative[k]=(serie.x[target][k+1]-2*serie.x[target][k]+serie.x[target][k-1])/sampling/sampling;
    }
    
//   if(normalize) {
//     s=get_statistics(derivative, serie.mask, serie.n, 0);
//     for(int k=0;k<serie.n;k++) {
//       if(serie.x[target][k]==serie.mask) continue;
//       serie.x[target][k]=(serie.x[target][k]-s.mean)/s.std;
//       }
//     }
    
  return(derivative);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *tidal_bore(double omega, double length, double dx, float h)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  examin propagation distorsion in shallow depths
  
  found to be quite limited

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  statistic_t s;
  double u,t=0,g=9.81,a=5.0,v=1.0,dt=300, T;
  FILE *out;
  
  T=12.5*3600;
  omega=2*M_PI/(12.5*3600);
  
  omega=2.0*omega;
  
  length=2.e+05;
  dx=100.;
  h=10;
  
  size_t size=length/dx+1;
  size_t nframes=T/dt+1;
  
  double *bore=new double[size];
  double *x=new double[size];
  double *linear=new double[size];
  double *k=new double[size];
  
  double *mgr=new double[nframes];
  double *mgr0=new double[nframes];
  
  double clinear=sqrt(g*h);
  double klinear=omega/clinear;
  
  for(int j=0;j<nframes;j++) {
    double t=j*dt;
    u=0;
    bore[0]=a*cos(omega*t);
    linear[0]=a*cos(omega*t);
    for(int n=1;n<size;n++) {
      x[n]=n*dx;
      double c=sqrt(g*(h+bore[n-1]));
      k[n]=omega/c;
      u+=k[n]*dx;
      bore[n]=a*cos(omega*t-u);
      linear[n]=a*cos(omega*t-klinear*x[n]);
      }
    mgr[j]=bore[size-1];
    mgr0[j]=linear[size-1];
    } 
//   out=fopen("bore.gnu","w");
//   for(int i=0;i<size;i++) {
//     fprintf(out,"%lf %lf %lf %lf %lf\n",x[i],bore[i],linear[i],bore[i]-linear[i],k[i]);
//     }
    
  out=fopen("bore.gnu","w");
  for(int i=0;i<nframes;i++) {
    fprintf(out,"%lf %lf %lf %lf\n",i*dt,mgr[i],mgr0[i],mgr[i]-mgr0[i]);
    }
    
  fclose(out);
        
  return(bore); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int getstep(tseries_t & serie, double sampling, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  statistic_t s;
  size_t pos, size, absolute, shift;
  double *p, value;
  double diff;
  FILE *out;
  
  vector<size_t> positions;
  vector<double> values;
  
  size=13.5/24.0/sampling;
  
  shift=size/10.0;
  for(int k=0;k<serie.n;k++) {
    if(serie.x[target][k]==serie.mask) serie.x[target][k]=NAN;
    }
    
  p=serie.x[target];
  absolute=0;
  while(true) {
    pos=maxpos(p, size, &value);
    if(pos!=-1) {
      positions.push_back(absolute+pos);
      values.push_back(value);
      if(positions.size()==1) diff=0.;
      else diff=serie.t[positions[positions.size()-1]]-serie.t[positions[positions.size()-2]];
      printf("max: %d %lf %lf %lf\n", absolute+pos, serie.t[absolute+pos], diff, value);
      }
    p=p+pos+shift;
    absolute+=pos+shift;
    if(positions[positions.size()-1]+size>serie.n) break;
    }
  
  for(int k=0;k<serie.n;k++) {
    if(isnan(serie.x[target][k])) serie.x[target][k]=serie.mask;
    }

  double *t=new double[positions.size()];
  double *x=new double[positions.size()];
  double *y=new double[positions.size()];
  
  for(int k=0;k<positions.size();k++) {
    x[k]=values[k];
    if(isnan(x[k])) x[k]=serie.mask;
    t[k]=serie.t[positions[k]];
    }
  for(int k=1;k<positions.size();k++) {
    y[k]=serie.t[positions[k]]-serie.t[positions[k-1]];
    if(y[k]<720/(24.*60.) or y[k]>790/(24.*60.)) {
      x[k]=serie.mask;
      y[k]=serie.mask;
      }
    }
  y[0]=y[1];
    
  double *residuals=new double[positions.size()];
  int nodal=1;
  double  mean;
  bool Z0=false;

  spectrum_t solved, spectrum = spectrum_init_ref("ESTUARINE", Z0);

  mgr_data_t *tmp= harmonic_analysis_with_parsing(x,serie.mask,t,residuals, &mean, positions.size(), spectrum, solved, nodal, stdout);
  for(int k=0;k<positions.size();k++) {
    if(residuals[k]==serie.mask) residuals[k]=NAN;
    if(x[k]==serie.mask) x[k]=NAN;
    if(y[k]==serie.mask) y[k]=NAN;
    }
  
  out=fopen("max.gnu","w");
  for(int k=1;k<positions.size();k++) {
    fprintf(out,"%lf %lf %lf %lf %lf\n",serie.t[positions[k]], y[k]*24.*60., x[k], residuals[k], x[k]-residuals[k]);
    }
    
  fclose(out);
  
  return(0); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_try(tseries_t *serie,const mooring_t & mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
 
  re-create missing values in time series
  
  Principles:
  
    low frequency contribution reconstructed from time series filtering
  
    high frequency contribution reconstructed from tidal prediction
    
  known to be failing in extreme conditions (estuarine flooding or draft)
  
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
  date_t start,final;
  extern int remove_outlayers(double *t, double *h, double mask, double sampling, int n);
  extern double fill_gaps(double *t, double *h, double mask, int n, double w, bool fit);

  std::string line;
  cout << endl << "-------- starting reconstruction --------" << endl << endl;
  
  start=poctime_getdatecnes(serie->t[0], 'd');
  final=poctime_getdatecnes(serie->t[serie->n-1], 'd');
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  mask flat signal
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  int nmax=6;
  for (int i=0; i< serie->n-nmax; i++) {
    bool  keep=false;
    for (int k=1; k<nmax;k++) {
      if(serie->x[0][i]!=serie->x[0][i+k]) {
        keep=true;
        break;
        }
      }
    if(not keep) {
      double value=serie->x[0][i];
      int k=0;
      while(serie->x[0][i+k]==value) {
        serie->x[0][i+k]=serie->mask;
        k++;
        }
      i+=k-1;
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  perform harmonic analysis to extract sinusoidal tidal contribution to HF sifnal
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  double *residuals=new double[serie->n];
  double *LF=new double[serie->n];
  int nodal=1;
  double  mean;
  bool Z0=false;
  double window=3.0;

//   spectrum_t solved, spectrum = spectrum_init_ref("SHELF", Z0);
//   spectrum_t solved, spectrum = spectrum_init_ref("COASTAL", Z0);
//   spectrum_t solved, spectrum = spectrum_init_ref("ESTUARINE-HF", Z0);
/*------------------------------------------------------------------------------
  long period needed because of non-linearities (Mf=M2-N2, MSf=M2-S2, etc...) */
  spectrum_t solved, spectrum = spectrum_init_ref("ESTUARINE", Z0);
//   spectrum.remove(wSa);
//   spectrum.remove(wSsa);
  
  status=loess1d_irregular(serie->n, window, serie->x[0], serie->t, serie->mask, LF);
  for (int i=0; i< serie->n; i++) {
    if(serie->x[0][i]!=serie->mask) serie->x[0][i]-=LF[i];
    }

  mgr_data_t *tmp= harmonic_analysis_with_parsing(serie->x[0],serie->mask,serie->t,residuals, &mean, serie->n, spectrum, solved, nodal, stdout);
 
  if(tmp==0) return(-1);
  
  tseries_t *extended=new tseries_t;

  extended->nparam=11;
 
  extended->x=new double*[extended->nparam];
  extended->x[0]=serie->x[0];
  extended->x[1]=residuals;
  
  for(int k=2;k<extended->nparam;k++) {
    extended->x[k]=new double[serie->n];
    }

  extended->t=serie->t;
  extended->mask=serie->mask;
  extended->n=serie->n;
  extended->mooring=serie->mooring;
 
  *serie=*extended;
  
  for (int i=0; i< serie->n; i++) {
    if(serie->x[0][i] == serie->mask) {
      for(int k=1;k<serie->nparam;k++) {
        serie->x[k][i]=serie->mask;
        }
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  reconstruct time serie : resample time series at fixed interval
  
  interpolate missing values when gaps do not exceed a given duration
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  identify alleged time sampling */
  double sampling=get_timesampling(serie->t, serie->mask, serie->n);
  sampling=floor(sampling*24.*3600.+0.5)/24/3600.;
  serie->sampling=sampling;
  printf("sampling appears to be %lf mn (%d data over %d tics,  %d missing\n",sampling*24.*60., serie->n, serie->maxsize(), serie->maxsize()-serie->n);
  
  printf("number of valid values before resampling: %d (over %d)\n",serie->n-serie->size(0,serie->mask),serie->n);
  tseries_t *reconstructed=new tseries_t;
  
/*------------------------------------------------------------------------------
  do not re-interpolate if 2 successive time larger than max_gap (assign mask)*/
  double max_gap=sampling+5.e-05;
  *reconstructed=serie->resample(sampling, max_gap);
  
  *serie=*reconstructed;
  printf("number of valid values after  resampling: %d (over %d),  %d missing, max gap=%lf hours \n",serie->n-serie->size(0,serie->mask),serie->n, serie->size(0,serie->mask), max_gap*24.);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  compute 1st derivatives
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  serie->x[2]=derivative_1(*serie, sampling, 0, false);
  serie->x[3]=derivative_1(*serie, sampling, 1, false);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  compute value and occurence time of max 1st derivatives for total and residual
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=getstep(*serie, sampling, 2);
  status=getstep(*serie, sampling, 3);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  compute 2nd derivatives
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  serie->x[4]=derivative_2(*serie, sampling, 0, false);
  serie->x[5]=derivative_2(*serie, sampling, 1, false);
  
  status=getstep(*serie, sampling, 4);
  status=getstep(*serie, sampling, 5);
  
  bool normalize=true;
  if(normalize) {
    for(int k=2;k<6;k++) serie->normalize(k);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  filtering of tidal residuals
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   double window;
  window=(25./8)/24./2.;
  window=1./8.;
  window=1./10.;
//   window=1./12.;
//   window=1./24.;
  status=loess1d_irregular(serie->n, window, serie->x[1], serie->t, serie->mask, serie->x[6]);
    
  for (int i=0; i< serie->n; i++) {
    if(serie->x[1][i] != serie->mask and serie->x[6][i] != serie->mask) {
      serie->x[7][i]=serie->x[1][i]-serie->x[6][i];
      }
    }
    
  status=harmonic_predictionTS(serie->x[8], serie->t, serie->n, tmp, solved.n, nodal);

  return(1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void RC(void)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 

  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double *u, *f, *s, *E, dt, value, T, *bore;
  double R, C, L;
  size_t size;
  double count;
  FILE* out;
  
  
  double dx, omega, length;
  float h;
  bore=tidal_bore(omega, length, dx, h);
  
  RCQ_simulator();
  
  tseries_t *serieA, *serieB;
  string filename;
  int ncolumns=1, status;
  const char *header_format=0;
  mooring_t mooring;
  char unit='m';
  const char *format=FORMAT_NAME_GNU;
  
  filename="HE_Rouen_20160601-20170731_20160531_20170731.gnu";
  mgr_load_timeserie(filename.c_str(), &serieA, unit, format, ncolumns, header_format, &mooring, &status);
        
  status=mgr_try(serieA, mooring);
	
  status=mgr_save_timeserie("try.gnu", mooring, serieA[0], 'm', format, true);
  
  printf("end of integration\n");
}
