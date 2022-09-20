
/*******************************************************************************

  T-UGO tools, 2006-2018

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief miscallenous tidal definitions
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>

#include <string.h>
#include <iostream>

#include "tools-define.h"
#include "tools-structures.h"

#include "poc-time.h"
#include "statistic.h"

#include "tides.def"

#include "spectrum.h"
#include "zapper.h"
#include "maths.h"
#include "statistic.h"

#include <gsl/gsl_statistics_double.h>

#ifndef VERBOSE
#define VERBOSE -1
#endif

/*------------------------------------------------------------------------------
 
http://en.wikipedia.org/wiki/Singular_value

http://en.wikipedia.org/wiki/Condition_number

------------------------------------------------------------------------------*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_analysis_old(double *serie,double *time,double *residuals, double *mean,int n, spectrum_t s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

time in cnes days

------------------------------------------------------------------------------*/
{
  int k,l,m,neq,status,step,nw,stepmax=100;
  double *A=NULL,*b=NULL;
  int *pivot=NULL;
  double **cs=NULL,**sn=NULL;
  tidal_wave wave;
  double V,V0,omega;
  float  *a=NULL,*G=NULL;
  double *zr=NULL,*zi=NULL,hmean,std;
  date_t start;
  extern date_t cnesdate(double t, char fmt);
  FILE *out=NULL;
  astro_angles_t astro_angles;

  nw=s.n;
  hmean=0.;
  for(m=0; m<n; m++) {
    residuals[m]=serie[m];
    hmean+=residuals[m];
    }

  hmean/=n;
  for(m=0; m<n; m++) {
    residuals[m]-=hmean;
    }

  start=poctime_getdatecnes(time[0], 'd');
  init_argument(&astro_angles,start);

  neq=2*s.n;
  exitIfNull(
    A=(double *) malloc(neq*neq*sizeof(double))
    );
  for(k=0; k<neq; k++)
    for(l=0; l<neq; l++) A[k*neq+l]=0;

  exitIfNull(
    b=(double *) malloc(neq*sizeof(double))
    );
  exitIfNull(
    zr=(double *) malloc(nw*sizeof(double))
    );
  exitIfNull(
    zi=(double *) malloc(nw*sizeof(double))
    );
  for(k=0; k<nw; k++) zr[k]=0;
  for(k=0; k<nw; k++) zi[k]=0;

  exitIfNull(
    a=(float *) malloc(nw*sizeof(float))
    );
  for(k=0; k<nw; k++) a[k]=0;
  exitIfNull(
    G=(float *) malloc(nw*sizeof(float))
    );
  for(k=0; k<nw; k++) G[k]=0;

  exitIfNull(
    cs=(double **) malloc(s.n*sizeof(double))
    );
  exitIfNull(
    sn=(double **) malloc(s.n*sizeof(double))
    );

  for(k=0; k < s.n; k++) {
      exitIfNull(
        cs[k]=(double *) malloc(n*sizeof(double))
        );
      exitIfNull(
        sn[k]=(double *) malloc(n*sizeof(double))
        );
    }

  exitIfNull(
    pivot=(int *) malloc(neq*sizeof(int))
    );

  for(k=0; k<s.n; k++) b[k]=0;

  for(k=0; k < s.n; k++) {
    wave=s.waves[k];
    omega=wave.omega*dph2rpd;
    V0=greenwhich_argument(astro_angles,wave);
    V0=0;
    for(m=0; m<n; m++) {
      V=omega*(time[m]-time[0])+V0;
      cs[k][m]=cos(V);
      sn[k][m]=sin(V);
      }
    }

  for(k=0; k < s.n; k++) {
    for(l=0; l<s.n; l++) {
      for(m=0; m<n; m++) {
        A[2*l*neq    +2*k]  +=cs[k][m]*cs[l][m];
        A[2*l*neq    +2*k+1]+=sn[k][m]*cs[l][m];
        A[(2*l+1)*neq+2*k]  +=cs[k][m]*sn[l][m];
        A[(2*l+1)*neq+2*k+1]+=sn[k][m]*sn[l][m];
        }
      }
    }

//   for(k=0; k < s.n; k++) {
//     A[2*k*neq    +2*k]   +=1.e-04*n*n;
//     A[(2*k+1)*neq+2*k+1] +=1.e-04*n*n;
//     }

  status=poc_getrf(neq, A, pivot);

  for(step=0;step<stepmax;step++) {
    for(k=0; k<neq; k++) b[k]=0;
    for(k=0; k < s.n; k++) {
      for(m=0; m<n; m++) {
        b[2*k]  +=cs[k][m]*residuals[m];
        b[2*k+1]+=sn[k][m]*residuals[m];
        }
      }

  status=poc_getrs(neq,1,A,pivot,b);


    /*    printf("step %d -----------------------\n",step);*/
    for(k=0; k < s.n; k++) {
      zr[k]+=b[2*k];
      zi[k]+=b[2*k+1];
      a[k]=sqrt(zr[k]*zr[k]+zi[k]*zi[k]);
      G[k]=360.0*atan2(zi[k],zr[k])/pi2;
      if(G[k] < 0.) G[k]+=360.;
/*       printf("%s wave: %f %f \n",s.waves[k].name,100*a[k],G[k]); */
      }

    for(k=0; k < s.n; k++) {
      for(m=0; m<n; m++) {
        residuals[m]-=cs[k][m]*b[2*k]+sn[k][m]*b[2*k+1];
        }
      }
    std=0;
    for(m=0; m<n; m++) {
      std+=residuals[m]*residuals[m];
      }
    std=sqrt(std/n);
    printf("step %d: std=%lf \n",step,std);
    }

  for(k=0; k < s.n; k++) {
    serie[2*k]=zr[k];
    serie[2*k+1]=zi[k];
    a[k]=sqrt(zr[k]*zr[k]+zi[k]*zi[k]);
    G[k]=360.0*atan2(zi[k],zr[k])/pi2;
    if(G[k] < 0.) G[k]+=360.;
//    printf("%s wave: %5.1f %5.1f %5.1f \n",s.waves[k].name,100*a[k],G[k],360./(s.waves[k].omega*24.));
    }

  out=fopen("spectrum.txt","w");
  if(out==0) {
    TRAP_ERR_EXIT(-1,"file opening issue : spectrum.txt \n");
    }

  for(k=0; k < s.n; k++) {
    serie[2*k]=zr[k];
    serie[2*k+1]=zi[k];
    a[k]=sqrt(zr[k]*zr[k]+zi[k]*zi[k]);
    G[k]=360.0*atan2(zi[k],zr[k])/pi2;
    if(G[k] < 0.) G[k]+=360.;
    fprintf(out,"%s %5.1f %5.1f %5.1f \n",s.waves[k].name,100*a[k],G[k],360./(s.waves[k].omega*24.));
    }

  fclose(out);

  out=fopen("residuals.txt","w");
  if(out==0) {
    TRAP_ERR_EXIT(-1,"file opening issue : residuals.txt \n");
    }

  for(k=0; k < n; k++) {
    fprintf(out,"%1f %lf\n",time[k],residuals[k]);
    }

  fclose(out);

  free((double *)A);
  free((double *)b);
  free(pivot);

  for(k=0; k < s.n; k++) free((double *)cs[k]);
  for(k=0; k < s.n; k++) free((double *)sn[k]);

  free(cs);
  free(sn);

  free(a);
  free(G);

  free(zr);
  free(zi);

  *mean=hmean;

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_analysis_old(float *serie,double *time,float *residuals, float *mean,int n, spectrum_t s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

time in cnes days

----------------------------------------------------------------------*/
{
  int k,l,m,neq,status,step,nw,stepmax=2,nrhs=1;
  double *A=NULL,*b=NULL;
  int *pivot=NULL;
  double **cs=NULL,**sn=NULL;
  tidal_wave wave;
  double V,V0,omega;
  float  *a=NULL,*G=NULL;
  double *zr=NULL,*zi=NULL,hmean;
  date_t start;
  extern date_t cnesdate(double t, char fmt);
  astro_angles_t astro_angles;

  nw=s.n;
  hmean=0.;
  for(m=0; m<n; m++) {
    residuals[m]=serie[m];
    hmean+=residuals[m];
    }
  hmean/=n;
  for(m=0; m<n; m++) {
    residuals[m]-=hmean;
    }

  start=poctime_getdatecnes(time[0], 'd');
  init_argument(&astro_angles,start);

  neq=2*s.n;
  if( (A=(double *) malloc(neq*neq*sizeof(double))) ==NULL) return(-1);
  for(k=0; k<neq; k++)
    for(l=0; l<neq; l++) A[k*neq+l]=0;

  if( (b=(double *) malloc(neq*sizeof(double))) ==NULL) return(-1);
  if( (zr=(double *) malloc(nw*sizeof(double))) ==NULL) return(-1);
  if( (zi=(double *) malloc(nw*sizeof(double))) ==NULL) return(-1);
  for(k=0; k<nw; k++) zr[k]=0;
  for(k=0; k<nw; k++) zi[k]=0;

  if( (a=(float *) malloc(nw*sizeof(float))) ==NULL) return(-1);
  for(k=0; k<nw; k++) a[k]=0;
  if( (G=(float *) malloc(nw*sizeof(float))) ==NULL) return(-1);
  for(k=0; k<nw; k++) G[k]=0;

  if( (cs=(double **) malloc(s.n*sizeof(double))) ==NULL) return(-1);
  if( (sn=(double **) malloc(s.n*sizeof(double))) ==NULL) return(-1);
  for(k=0; k < s.n; k++) if( (cs[k]=(double *) malloc(n*sizeof(double))) ==NULL) return(-1);
  for(k=0; k < s.n; k++) if( (sn[k]=(double *) malloc(n*sizeof(double))) ==NULL) return(-1);

  if( (pivot=(int *) malloc(neq*sizeof(int))) ==NULL) return(-1);

  for(k=0; k<s.n; k++) b[k]=0;

  for(k=0; k < s.n; k++) {
    wave=s.waves[k];
    omega=wave.omega*dph2rpd;
    V0=greenwhich_argument(astro_angles,wave);
    for(m=0; m<n; m++) {
      V=omega*(time[m]-time[0])+V0;
      cs[k][m]=cos(V);
      sn[k][m]=sin(V);
      }
    }

  for(k=0; k < s.n; k++) {
    for(l=0; l<s.n; l++) {
      for(m=0; m<n; m++) {
        A[2*l*neq    +2*k]  +=cs[k][m]*cs[l][m];
        A[2*l*neq    +2*k+1]+=sn[k][m]*cs[l][m];
        A[(2*l+1)*neq+2*k]  +=cs[k][m]*sn[l][m];
        A[(2*l+1)*neq+2*k+1]+=sn[k][m]*sn[l][m];
        }
      }
    }

  status=poc_getrf(neq, A, pivot);

  for(step=0;step<stepmax;step++) {
    for(k=0; k<neq; k++) b[k]=0;
    for(k=0; k < s.n; k++) {
      for(m=0; m<n; m++) {
         b[2*k]  +=cs[k][m]*residuals[m];
         b[2*k+1]+=sn[k][m]*residuals[m];
        }
      }

    status=poc_getrs(neq, nrhs, A, pivot, b);


    /*    printf("step %d -----------------------\n",step);*/
    for(k=0; k < s.n; k++) {
      zr[k]+=b[2*k];
      zi[k]+=b[2*k+1];
      a[k]=sqrt(zr[k]*zr[k]+zi[k]*zi[k]);
      G[k]=360.0*atan2(zi[k],zr[k])/pi2;
      if(G[k] < 0.) G[k]+=360.;
/*       printf("%s wave: %f %f \n",s.waves[k].name,100*a[k],G[k]); */
      }

    for(k=0; k < s.n; k++) {
      for(m=0; m<n; m++) {
        residuals[m]-=cs[k][m]*b[2*k]+sn[k][m]*b[2*k+1];
        }
      }
    }

  for(k=0; k < s.n; k++) {
    serie[2*k]=zr[k];
    serie[2*k+1]=zi[k];
    a[k]=sqrt(zr[k]*zr[k]+zi[k]*zi[k]);
    G[k]=360.0*atan2(zi[k],zr[k])/pi2;
    if(G[k] < 0.) G[k]+=360.;
//     printf("%s wave: %5.1f %5.1f \n",s.waves[k].name,a[k],G[k]);
    }

  free((double *)A);
  free((double *)b);
  free(pivot);

  for(k=0; k < s.n; k++) free((double *)cs[k]);
  for(k=0; k < s.n; k++) free((double *)sn[k]);

  free(cs);
  free(sn);

  free(a);
  free(G);

  free(zr);
  free(zi);

  *mean=hmean;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double tide_alias(tidal_wave w, double Trepet)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double omega = tide_omega(w,"dph"); // true pulsation in deg/hour
  double alias = omega;

  if(360. / omega  < 2 * Trepet * 24.){                  // comparing periods in hour
    double nRev = omega / 360. * Trepet * 24.;
    double Ta = Trepet/ fabs(nRev  - floor(nRev +0.5));  // aliased period in days
    alias = 360. / (Ta * 24.);                           // aliased pulsation in deg/hour
#ifdef DEBUG
    printf("%10s period %12.6f (alias)\n",w.name, Ta);
#endif
    }
  else{
#ifdef DEBUG
    printf("%10s period %12.6f (true)\n",w.name, tide_period(w) / 24.);
#endif
    }

  return(alias);// aliased/true pulsation in deg/hour
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void spectrum_define_alias(spectrum_t *s, double Trepet, FILE *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  for(int i = 0; i < s->n; i++){
    s->waves[i].aliased = tide_alias(s->waves[i],Trepet);
    }
  spectrum_print(*s,out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void spectrum_redef(spectrum_t *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------------
  check redundant (aliased) frequencies  - need XXX.aliased field already set */
/**----------------------------------------------------------------------------
  should not be used anymore */

  int *keep = new int[s->n];
  for(int i = 0; i < s->n; i++){
    keep[i] = 1;
    }

  spectrum_isRedundant(*s, keep);

  // clean spectrum
  spectrum_t tmp;
  tmp.n = 0;
  tmp.nmax = s->n;
  size_t size = sizeof(tidal_wave);
  for(int i = 0; i < s->n; i++){
    if(keep[i] != 0) {
      tmp.n++;
      }
    }
  tmp.waves=new tidal_wave[tmp.n];

  for(int i = 0, n = 0; i < s->n; i++){
    if(keep[i] != 0) {
      memmove(&(tmp.waves[n++]), &(s->waves[i]), size);
      }
    else {
#if 1 || VERBOSE > 0
      printf("wave: %10s discarded\n", s->waves[i].name);
#endif
      }
    }

  *s = tmp;

  tmp.destroy();
  delete [] keep;
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void spectrum_isSolvable(spectrum_t s, double duration, int *keep, int *resolve, FILE *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double tau;

  for(int i = 0; i < s.n; i++) {
    if(keep[i] == 0) {
      continue;
      }
    if(s.waves[i].aliased == 0.0){
      tau=1.e+10;
      }
    if(strcmp(s.waves[i].name,"Z0") == 0.0){
      continue; /*Z0 special case */
      }
    tau = (1. / s.waves[i].aliased) * 360. / 24.; // true/aliased period in days
    if(fabs(tau) > duration) {
      resolve[i] = 0;
      keep[i]    = 0;
      if(out!=0) {
        printf("wave: %10s not directly solvable, period=%9.3f days duration=%9.3f days\n", s.waves[i].name,tau, duration);
        }
      }
  }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void spectrum_isRedundant(spectrum_t s, int *keep)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
check redundant (aliased) frequencies  - poorly usable, redundant with separatiob check */
{
double tau;

  for(int i = 0; i < s.n; i++) {
    for(int j = i + 1; j < s.n; j++) {
      tau = fabs(s.waves[j].aliased - s.waves[i].aliased);
      if(tau == 0.) {
        keep[j] = 0;
        }
      }
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void spectrum_checkRayleigh(spectrum_t s, double duration, int *keep, int *resolve, FILE *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double tau,delta;
  double dph2d = 360. / 24.;
  double aliasT_i = 0, aliasT_j = 0;
  bool excluded;
  
  if(out!=0) {
    fprintf(out, "separation issues (duration=%9.3f years)\n",duration/365.25);
    }
  for(int i = 0; i < s.n; i++) {
    aliasT_i = dph2d / s.waves[i].aliased ; // aliased period in days
    if(keep!=0) {
      if(keep[i] == 0)
        continue;
      }
    for(int j = i + 1; j < s.n; j++) {
      aliasT_j = dph2d / s.waves[j].aliased ;
      delta = fabs(s.waves[i].aliased - s.waves[j].aliased);
      excluded=false;
      if(delta==0) { // non-separable constituents
        if(keep!=0) {
          keep[j] = 0;
          }
        excluded=true;
        }
      if(excluded) {
        if(out!=0) {
          fprintf(out, "waves: %10s %10s, separation: %9.3f \n",s.waves[i].name, s.waves[j].name, 1./0.);
          continue;
          }
        }
      delta = fabs(aliasT_i - aliasT_j);
      tau = (aliasT_i * aliasT_j) / delta;    // separation period in days
      excluded=false;
      if(tau > duration) { // non-separable constituents
        if(keep!=0) {
          keep[j] = 0;
          }
        excluded=true;
        }
      if(excluded) {
        if(out!=0) {
         if(tau>1000.*365.25) {
            fprintf(out, "waves: %10s %10s, separation: ~%9.3f \n",s.waves[i].name, s.waves[j].name, 1./0.);
            }
         else if(tau>100.*365.25) {
            fprintf(out, "waves: %10s %10s, separation:  %9.3f (centuries) \n",s.waves[i].name, s.waves[j].name, tau/365.25/100.);
            }
          else if(tau>365.25) {
            fprintf(out, "waves: %10s %10s, separation:  %9.3f (years) \n",s.waves[i].name, s.waves[j].name, tau/365.25);
            }
          else {
            fprintf(out, "waves: %10s %10s, separation:  %9.3f (days) \n",s.waves[i].name, s.waves[j].name, tau);
            }
          }
        }
      }
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void spectrum_isAdmittanceSolvable(spectrum_t s, double duration, int *keep, int *resolve, int *deduce, FILE *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
   */
{
  admittance_t admittance;
  
  int *method = admittance_mts_check(&admittance, s, resolve);

  for(int i = 0; i < s.n; i++) {
    if(keep[i] == 0) {
      if(s.waves[i].Ap != 0.0) {
/*------------------------------------------------------------------------------
        admittance can do ? */
        if(admittance_mts_verify(s.waves[i],method) != 0) {
          if(out!=0) {
            fprintf(out, "wave: %10s lost for admittance mode\n",s.waves[i].name);
            }
          continue;
          }
/*------------------------------------------------------------------------------
        yes it can */
        keep[i]   = 1;
        deduce[i] = 1;
        if(out!=0) {
          fprintf(out, "wave: %10s pushed back in admmittance mode\n",s.waves[i].name);
          }
        }
      else {
        if(out!=0) {
          fprintf(out, "wave: %10s lost for admittance mode\n",s.waves[i].name);
          }
        }
      }
    }
  
  delete[] method;
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void spectrum_reset(spectrum_t& s, const string spectrumType, const double Trepet, FILE *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  s.destroy();
  
  bool Z0=false;
  s = spectrum_init_ref(spectrumType,Z0);
  
  spectrum_define_alias(&s, Trepet, out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mgr_data_t *harmonic_analysis_with_parsing(const double *serie, double mask,const double *time, double *residuals, double *mean, int n, spectrum_t prior, spectrum_t & s, int nodal, FILE *log)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**-----------------------------------------------------------------------------
 
  Extremely unsafe !!! as the proposed spectrum has been already parsed to allow
  an optimal analysis, distinguishing unsolvable, directly solvable and
  admittance solvable constituents.
  
  Hence re-processing spectrum parsing can lead to highly unstable and poorly
  manageable code.
  
  time in cnes days
  
-----------------------------------------------------------------------------*/

/* *-----------------------------------------------------------------------------
  remove mean*/
  statistic_t stat;
  int verbose=(log!=0);
  stat=get_statistics(serie, mask, n, verbose);
  
  *mean=stat.mean;

  for(size_t m = 0; m < n; m++) {
    if(serie[m]!=mask) {
      residuals[m] = serie[m] - *mean;
      }
    else {
      residuals[m] = mask;
      }
    }
//   int nmasked=occurence<double>(mask, residuals, n);

/* *----------------------------------------------------------------------
  check spectrum */
  double duration = time[n-1] - time[0];
  s=prior;
  
  int *keep    = aset(s.n,1);
  int *deduce  = aset(s.n,0);
  int *direct  = new int[s.n];
  int *resolve = aset(s.n,1);
  
//  FILE *debuglog  =stdout;
  FILE *debuglog  =0;

/**----------------------------------------------------------------------------
  spectrum parsing */
/*----------------------------------------------------------------------------
  check aliased period against serie duration */
  spectrum_isSolvable(s, duration, keep, resolve, debuglog);
/*----------------------------------------------------------------------------
  check aliased period separation */
  spectrum_checkRayleigh(s, duration, keep, direct, debuglog);
  valcpy(direct,keep,s.n);
/*----------------------------------------------------------------------------
  check admittance salvation*/
  spectrum_isAdmittanceSolvable(s, duration, keep, resolve, deduce, debuglog);
/*----------------------------------------------------------------------------
  remove unsolved constituents */
  spectrum_reduce(s, keep, deduce, debuglog);
  
  if(s.n <= 0){
    return(0);
    }
  admittance_t admittance;
  int *method=admittance_mts_check(&admittance, s, direct);

  int nsolved = 0;
  int nvalid  = 0;
  for(size_t k = 0; k < s.n; k++){
    size_t m = s.waves[k].nT;
    if( (keep[k] == 1) && (deduce[k] == 0) ){
      nsolved++;
      }
    if(m > 2){
      continue;
      }
    if(deduce[k] == 1){
      if(admittance_mts_verify(s.waves[k],method) != 0){
        keep[k] = 0;
        deduce[k] = 0;
//#if VERBOSE > 2
        printf("remove %s from analysis list\n", s.waves[k].name);
//#endif
        }
      }
    }

/* *----------------------------------------------------------------------
  start analysis process */

//  if( (method[0] == 0) && (method[1] == 0) && (method[2] == 0) ) return(0);
  delete[] method;

  int nw = s.n;
  int neq = 2 * s.n;

  double *A = 0, *b = 0;
  if( (A = new double[neq * neq]) == 0) return(0);
  aset(A, neq * neq, 0.);

  if( (b = new double[neq]) == 0) return(0);
  aset(b, neq, 0.);

  double **cs = 0, **sn = 0;
  if( (cs = new double*[nw]) == 0) return(0);
  if( (sn = new double*[nw]) == 0) return(0);
  for(size_t k = 0; k < nw; k++){
    if( (cs[k] = new double[n]) == 0) return(0);
    if( (sn[k] = new double[n]) == 0) return(0);
    }

/* *-----------------------------------------------------------------------------
  compute sinus and cosinus at sampling time, with or without nodal corrections*/
  

  if(nodal == 1){
#pragma omp parallel for
    for(size_t m = 0; m < n; m++) {
      astro_angles_t astro_angles;
      tidal_wave wave;
      double f, u, V;
      init_argument(&astro_angles,poctime_getdatecnes(time[m], 'd'));
      for(size_t k = 0; k < s.n; k++) {
        wave = s.waves[k];
        V = greenwhich_argument(astro_angles,wave);
        u = nodal_phase(astro_angles,wave);
        f = nodal_factor(astro_angles,wave.formula);
        cs[k][m] = f * cos(V + u);
        sn[k][m] = f * sin(V + u);
        }
      }
    }
  else{
    for(size_t k = 0; k < s.n; k++) {
      astro_angles_t astro_angles;
      tidal_wave wave;
      double V;
      wave = s.waves[k];
      for(size_t m = 0; m < n; m++) {
        init_argument(&astro_angles,poctime_getdatecnes(time[m], 'd'));
        V = greenwhich_argument(astro_angles,wave);
        cs[k][m] = cos(V);
        sn[k][m] = sin(V);
        }
      }
    }

/* *-----------------------------------------------------------------------------
  fill harmonic analysis matrix*/
#pragma omp parallel for
  for(size_t k = 0; k < s.n; k++) {
    for(size_t l = 0; l < s.n; l++) {
      for(size_t m = 0; m < n; m++) {
        if(serie[m]==mask) continue;
        A[2 * l * neq       + 2 * k]     += cs[k][m] * cs[l][m];
        A[2 * l * neq       + 2 * k + 1] += sn[k][m] * cs[l][m];
        A[(2 * l + 1) * neq + 2 * k]     += cs[k][m] * sn[l][m];
        A[(2 * l + 1) * neq + 2 * k + 1] += sn[k][m] * sn[l][m];
        }
      }
    }

/* *--------------------------------------------------------------------
  adjust matrix for unresolved/admitted waves*/

  double coef[3] = {1., 1., 1.};

  for(size_t i = 0; i < s.n; i++) {
    if(keep[i] == 0) {
      for(size_t j = 0; j < s.n; j++)
        A[2 * j * neq + 2 * i] = 0.;
      for(size_t j = 0; j < s.n; j++)
        A[2 * j * neq + 2 * i + 1] = 0.;
      for(size_t j = 0; j < s.n; j++)
        A[(2 * j + 1) * neq + 2 * i] = 0.;
      for(size_t j = 0; j < s.n; j++)
        A[(2 * j + 1) * neq + 2 * i + 1] = 0.;
      A[2 * i * neq + 2 * i] = 1.0;
      A[(2 * i + 1) * neq + 2 * i + 1] = 1.0;
      if(log!=0) {
        fprintf(log,"wave: %10s discarded\n", s.waves[i].name);
        }
      }
    else {
      if(log!=0) {
        printf("wave: %10s analysed", s.waves[i].name);
        }

/**--------------------------------------------------------------------

      h=()/potential;

      h-a cos(f-f0) -b sin (f-f0) -c = 0

      a= ri[0][0] h0 + ri[1][0] h1 + ri[2][0] h2
      b= ri[0][1] h0 + ri[1][1] h1 + ri[2][1] h2
      c= ri[0][2] h0 + ri[1][2] h1 + ri[2][2] h2

      h-(ri[0][0] cos(f-f0) + ri[0][1] sin (f-f0) ri[0][2]) h0
       -(ri[1][0] cos(f-f0) + ri[1][1] sin (f-f0) ri[1][2]) h1
       -(ri[2][0] cos(f-f0) + ri[2][1] sin (f-f0) ri[2][2]) h2
       =h- (c0 h0 +c1 h1 + c2 h2)
       =0

 --------------------------------------------------------------------*/

      if(deduce[i] == 1) {
        if(log!=0) {
          fprintf(log," in admittance mode\n");
          }
        for(size_t j = 0; j < neq; j++)
          A[j * neq + 2 * i] = 0.;
        for(size_t j = 0; j < neq; j++)
          A[j * neq + 2 * i + 1] = 0.;
        A[2 * i * neq + 2 * i] = 1. / s.waves[i].Ap;
        A[(2 * i + 1) * neq + 2 * i + 1] = 1. / s.waves[i].Ap;

        size_t k = s.waves[i].nT;
        switch(admittance.method[k]){
          case ADMITTANCE_UNAVAILABLE: {
            return(0);
            }
          case ADMITTANCE_SPLINE: {
            admittance_mts_sweightP(admittance, s.waves[i], coef);
            for(size_t l = 0; l < 3; l++) {
              size_t j = admittance.windex_s[k][l];
              A[2 * j * neq + 2 * i] = -coef[l];
              A[(2 * j + 1) * neq + 2 * i + 1] = -coef[l];
              }
            break;
            }
          case ADMITTANCE_LINEAR: {
            admittance_mts_lweightP(admittance, s.waves[i], coef);
            for(size_t l = 0; l < 2; l++) {
              size_t j = admittance.windex_l[k][l];
              A[2 * j * neq + 2 * i] = -coef[l];
              A[(2 * j + 1) * neq + 2 * i + 1] = -coef[l];
              }
            break;
            }
          default: return(0);
          }
        }
      else {
        if(log!=0) {
          printf("\n");
          }
        }
      }
    }

//  admittance_mts_terminate(admittance);

  int *pivot = 0;
  if( (pivot = new int[neq]) == 0) return(0);

  int status = neq+2;
  status=poc_getrf(neq, A, pivot);
  
  if(status!=0) return(0);

  double *zr = 0, *zi = 0;
  if( (zr = new double[nw]) == 0) return(0);
  aset(zr, nw, 0.);
  if( (zi = new double[nw]) == 0) return(0);
  aset(zi, nw, 0.);

  float  *a = 0, *G = 0;
  if( (a = new float[nw]) == 0) return(0);
  aset(a, nw, 0.f);
  if( (G = new float[nw]) == 0) return(0);
  aset(G, nw, 0.f);

  int stepmax = 2;
  int nrhs = 1;
  for(size_t step = 0; step < stepmax; step++) {
    aset(b, neq, 0.);
    for(size_t k = 0; k < s.n; k++) {
      for(size_t m = 0; m < n; m++) {
        if(residuals[m]!=mask) {
          b[2 * k]     += cs[k][m] * residuals[m];
          b[2 * k + 1] += sn[k][m] * residuals[m];
          }
        }
      }

/* *--------------------------------------------------------------------
    adjust right-hand side members for unresolved/admitted waves*/
    for(size_t i = 0; i < s.n; i++) {
/*------------------------------------------------------------------------------
      solution forced to zero*/
      if(keep[i] == 0) {
        b[2 * i]     = 0.;
        b[2 * i + 1] = 0.;
        }
/*------------------------------------------------------------------------------
      admittance rhs is zero*/
      if(deduce[i] == 1) {
        b[2 * i]     = 0.;
        b[2 * i + 1] = 0.;
        }
      }
    status=poc_getrs(neq, nrhs, A, pivot, b);

    for(size_t k = 0; k < nw; k++) {
      zr[k] += b[2 * k];
      zi[k] += b[2 * k + 1];
      a[k] = sqrt(zr[k] * zr[k] + zi[k] * zi[k]);
      G[k] = 360.0 * atan2(zi[k],zr[k]) / (2 * M_PI);
      if(G[k] < 0.) G[k] += 360.;
      }
    for(size_t k = 0; k < nw; k++) {
      for(size_t m = 0; m < n; m++) {
        if(residuals[m]!=mask) {
          residuals[m] -= cs[k][m] * b[2 * k] + sn[k][m] * b[2 * k + 1];
          }
        }
      }
    } // end iterative analysis

/* *----------------------------------------------------------------------------
  compute amplitudes and phase lags*/
  for(size_t k = 0; k < nw; k++){
    if(keep[k] == 0) {
      a[k] = 0.;
      G[k] = 0.;
      continue;
      }
    a[k] = sqrt(zr[k] * zr[k] + zi[k] * zi[k]);
    G[k] = 360.0 * atan2(zi[k],zr[k]) / (2 * M_PI);
    if(G[k] < 0.) G[k] += 360.;
    }

/* *----------------------------------------------------------------------------
  compute error on estimates, based on:

  Data Analysis Methods in Physical Oceanography
  by William J.; Thomson Emery
  Pub date: 1998

  -----------------------------------------------------------------------------*/

#if defined(ATLAS) && ATLAS == 1
  double s2 = cblas_ddot(n, residuals, 1, residuals, 1);
#else
  double s2 = 0;
  for(size_t i = 0; i < n; i++){
    if(residuals[i]!=mask) {
      s2 += residuals[i] * residuals[i];
      }
    }
#endif
/* *----------------------------------------------------------------------------
  unbiased estimator of the residual variance */
  s2 /= (n - 2 * nsolved - 1);

  vector<double> f_error(0, neq);
  double c = 0.;
  double *mu    = new double[neq];
  double *bck   = new double[neq];
  for(size_t i = 0; i < neq; i++){
    aset(mu, neq, 0.);
/* *----------------------------------------------------------------------------
    commentaires ? */
    mu[i] = 1.f;
    memcpy(bck, mu, neq * sizeof(double));
    
    for(size_t i = 0; i < s.n; i++) {
/*------------------------------------------------------------------------------
      solution forced to zero*/
      if(keep[i] == 0) {
        mu[2 * i]     = 0.;
        mu[2 * i + 1] = 0.;
        }
/*------------------------------------------------------------------------------
      admittance rhs is zero*/
      if(deduce[i] == 1) {
        mu[2 * i]     = 0.;
        mu[2 * i + 1] = 0.;
        }
      }

/* *----------------------------------------------------------------------------
    modified normal equations */
    status=poc_getrs(neq, nrhs, A, pivot, mu);
    
    for(size_t j = 0; j < neq; j++){
      c += bck[j] * mu[j];
      }
/* *----------------------------------------------------------------------------
    error on i-th element estimates */
    f_error.push_back(sqrt(c * s2));
    }
  delete[] mu;
  delete[] bck;

  
/* *----------------------------------------------------------------------------
  overwrite error for waves in admittance list  */
  status = admittance_mts_error(admittance, s, deduce, f_error);
  if(status!=0){
    return(0);
    }
  
  admittance_mts_terminate(admittance);
  
/* *----------------------------------------------------------------------------
  store elements */
  mgr_data_t *data;
  for(size_t k = 0; k < nw; k++) {
    if(keep[k]==1) nvalid++;
    }
  data = new mgr_data_t[nw];
  for(size_t k = 0; k < nw; k++) {
    data[k].amp = a[k];
    data[k].phi = G[k];
    data[k].constituent=s.waves[k];
    data[k].error = sqrt( f_error[2 * k] * f_error[2 * k]
        + f_error[2 * k + 1] * f_error[2 * k + 1] ); // set err = sqrt( s²(i) + s²(i+1) )
    }
  f_error.clear();
  stat=get_statistics(residuals, mask, n, verbose);

  delete[] A;
  delete[] b;
  delete[] pivot;

  for(size_t k = 0; k < nw; k++) delete[] cs[k];
  for(size_t k = 0; k < nw; k++) delete[] sn[k];

  delete[] cs;
  delete[] sn;

  delete[] a;
  delete[] G;

  delete[] zr;
  delete[] zi;

  delete[] keep;
  delete[] deduce;
  
//  stat=get_statistics(residuals, mask, n);

  return(data);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_parsing(const double *A, const spectrum_t & s, int *keep, int *deduce, double maxratio, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/// http://www.math.ualberta.ca/~wenchen/reference/akram_ir_siam_rev.pdf
{
  
  maxratio=0.3;
  
  if(verbose){
/*------------------------------------------------------------------------------
    print matrix */
    printf("%9s ","wave ");
    for(size_t k = 0; k < s.n; k++) {
      printf("%6s ",s.waves[k].name);
      }
    printf("\n");
  
    for(size_t k = 0; k < s.n; k++) {
      printf("%6s %3d",s.waves[k].name,k);
      double rk,ik;
      matrix_auto_correlation(A,s,k,&rk,&ik);
      
      for(size_t l = 0; l < s.n; l++) {
        double m;
        m=matrix_correlation(A,s,k,l,rk,ik);
        printf(" %6.3f",m);
        }
      printf("\n");
      }
    }
  
  for(size_t k = 0; k < s.n; k++) {
/*------------------------------------------------------------------------------
    auto-correlation */
    double rk,ik;
    matrix_auto_correlation(A,s,k,&rk,&ik);
    
    for(size_t l = k+1; l < s.n; l++) {
/*------------------------------------------------------------------------------
      cross-correlation */
      double m;
      m=matrix_correlation(A,s,k,l,rk,ik);
      if(m>1.1){
        string path="matrix_gives_error.csv";
        save_full_matrix(path,0,-1,s,A,keep, deduce);
        TRAP_ERR_EXIT(ENOEXEC,"programming error : %g with %s and %s\n Matrix saved to "+path+"\n",m, s.waves[k].name, s.waves[l].name);
        }
      
      double tau=deltaOmega2separation(s.waves[k].omega-s.waves[l].omega);
      double tau_aliased=deltaOmega2separation(s.waves[k].aliased-s.waves[l].aliased);
      if(m>0.35){
        if(verbose==1) printf("%6s %6s %d %d %d %d  %6.3f separation:  true=%6.1f  aliased=%8.1f\n", s.waves[k].name, s.waves[l].name, keep[k], keep[l], deduce[k], deduce[l],m, tau, tau_aliased);
        }
      if(keep[k]==0) continue;
      
/*------------------------------------------------------------------------------
      ratio */
      if(m>0.3){
        if(keep[k]==1) {
//        deduce[l]=1;
          keep[l]=0;
          printf("%6s has pushed away %6s : %d %d %d %d normalized correlation=%6.3f separation:  true=%6.1f  aliased=%8.1f\n",
          s.waves[k].name, s.waves[l].name, keep[k], keep[l], deduce[k], deduce[l], m, tau, tau_aliased);
          }
        }
      }
    }
  
  const int k=s.wave_index("K1");
  for(size_t i=0;k>=0 && i<2;i++){
    const char *name;
    switch(i){
    case 0:
//     if(mission=="ERS")
      name="Sa";
//     else if(mission=="T/P")
//       name="Ssa";
//     else
//      return 0;
      break;
    case 1:
      name="S2";
      break;
    default:
      return 0;
      }
    
    int l=s.wave_index(name);
    if(l<0) return 0;
    
    double m;
    double rkrk,ikik;
    matrix_member(A,s,k,k,&rkrk,&ikik);
    const double
      rk=sqrt(rkrk),
      ik=sqrt(ikik);
    m=matrix_correlation(A,s,k,l,rk,ik);
    
    printf("%s %.6f correlated with K1\n",name,m);
    /* TODO : statistics of the output of the above line to choose the right factor below.
     COMMAND TIP: sed -re 's/[^ ]+ +([^ ]+) correlated with K1/\1/g;t;d' */
    if(m>0.3) return 1;
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double matrix_condition(double *A, spectrum_t s, int *keep, int *deduce, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int neq = 2 * s.n;
  int status, dim=neq;
  
  int LWORK  =1 + 6*dim + 2*dim*dim;
  int LIWORK =3 + 5*dim;
  int ITYPE=1;
  char JOB='V', UPLO='U';

//#define STATICS

#ifdef STATICS
   double WORK[LWORK],W[dim];
   int IWORK[LIWORK];
   double A[dim*dim],B[dim*dim];
#else
  double *WORK,*W;
  int    *IWORK;
//   double *A;
  double *M,*B;
  
  WORK =     new double [LWORK];
  IWORK =    new int    [LIWORK];
  W =        new double [dim];
  M =        new double [dim*dim];
  B =        new double [dim*dim];
#endif
  for(int m=0;m<dim*dim;m++) B[m]=0;
  for(int m=0;m<dim*dim;m++) M[m]=A[m];
  for(int m=0;m<dim;m++)    B[m*dim+m]=1.0;
  dsygvd_(&ITYPE, &JOB, &UPLO, &dim, M, &dim, B, &dim, W, WORK, &LWORK, IWORK, &LIWORK, &status);
  
  double condition=W[0]/W[dim-1];
  
//   for(size_t k = 0; k < s.n; k++) {
//     double a=M[2 * k * neq       + 2 * k];
//     double b=M[2 * k * neq       + 2 * k + 1];
//     double c=sqrt(a*a+b*b);
//     printf("%d ", k);
//     for(size_t l = 0; l < s.n; l++) {
//       double aa=M[2 * k * neq       + 2 * l];
//       double bb=M[2 * k * neq       + 2 * l + 1];
//       double cc=sqrt(aa*aa+bb*bb)/c;
//       if(cc>1.e-02){
//         double tau=deltaOmega2separation(s.waves[k].omega-s.waves[l].omega);
//         double tau_aliased=deltaOmega2separation(s.waves[k].aliased-s.waves[l].aliased);
//         printf("%6s ", s.waves[l].name);
//         }
//       }
//     printf("\n");
//     }
  
  delete[] M;
  delete[] B;
  delete[] W;
  delete[] WORK;
  delete[] IWORK;
 
  return(condition);
 
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cumulate(const spectrum_t & AnalysisList, mgr_data_t *data, const hconstant_t & prior)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{ /** extremely unsafe: size of 'data' array must be length of list 'AnalysisList' */

  size_t j = -1;
  int anomaly=0;
  
  complex<float> zPrior = 0, zAnalysis = 0;
  
  for(size_t w = 0; w < AnalysisList.n; w++) {
    if( (j = index_from_name(AnalysisList,AnalysisList.waves[w].name)) == w ){
      zAnalysis = polar(data[j].amp,-data[j].phi*d2r);
      }
    if( (j = index_from_name(*(prior.s),AnalysisList.waves[w].name)) != -1 ){
      zPrior = polar((double) prior.a[j],-prior.G[j]*d2r);
      zAnalysis += zPrior;
      }
    
    data[w].amp =  abs(zAnalysis);
    data[w].phi = -arg(zAnalysis) * r2d;
    if(data[w].phi < 0){
      data[w].phi += 360;
      }
    if(data[w].amp>10.0) {
      printf("warning %s : restored=%f (analysed=%f prior=%f)\n",AnalysisList.waves[w].name,data[w].amp,abs(zAnalysis),abs(zPrior));
      anomaly=1;
      }
    }
  if(anomaly==0) {
    return(0);
    }
  else {
    return(-1);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int initialise(const spectrum_t & spectrum, const hconstant_t & prior, double *zr, double *zi)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{ /** extremely unsafe: size of 'data' array must be length of list 'AnalysisList' */

  size_t j = -1;
  
  complex<double> zPrior = 0;
  
  for(size_t w = 0; w < spectrum.n; w++) {
    zPrior=0.0;
    if(prior.s!=0) {
      if( (j = index_from_name(*(prior.s),spectrum.waves[w].name)) != -1 ){
        zPrior = polar((double) prior.a[j], prior.G[j]*d2r);
        }
      }
    
    zr[w]=real(zPrior);
    zi[w]=imag(zPrior);
    }
  return(0);
}


int matrix_csv_i=0;
const char *matrixPath="matrix.csv";


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mgr_data_t *harmonic_analysis(const double *serie, double *residuals,const double *time, int n, const spectrum_t & s, double averaging_time, int* & keep, int* & deduce, const hconstant_t & prior, int nodal, double maxratio, FILE *log, int tag)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  double *b=0;
  mgr_data_t *analysed;
  int *resolve;
  harmonic_t x;
  admittance_t admittance;
  double condition[2];
//   double averaging_time=0.0;
  
  keep=new int[s.n];
  deduce=new int[s.n];
  resolve=new int[s.n];
  
  for(size_t k = 0; k < s.n; k++) {
    keep[k]=1;
    resolve[k]=1;
    deduce[k]=0;
    }
  
  double duration = time[n - 1] - time[0];
// /*----------------------------------------------------------------------------
//   check aliased period against serie duration */
//   spectrum_isSolvable(s, duration, keep, resolve, log);
// /*----------------------------------------------------------------------------
//   check aliased period separation */
//   spectrum_checkRayleigh(s, duration, keep,resolve, log);
    
  status=harmonic_init_new(x, time, n, s, nodal, averaging_time);
  if(status!=0) TRAP_ERR_EXIT(-1, "failure\n");
  
  condition[0]=matrix_condition(x.A, x.spectrum, keep, deduce,0);
  
/*----------------------------------------------------------------------------
  check solvable frequencies*/
  x.solvable(keep);
  
/*----------------------------------------------------------------------------
  check component interaction*/
  status=matrix_parsing(x.A, x.spectrum, keep, deduce, maxratio, 0);
  if(0/*||tag==11683*/){
    matrix_csv_i=0;
    save_full_matrix(matrixPath, &matrix_csv_i, tag, x.spectrum, x.A, keep, deduce);
    TRAP_ERR_EXIT(ENOEXEC,"testing %d\n",status);
    }
  if(status!=0){
    printf("cross correlation too high: abort analysis\n");
    save_full_matrix(matrixPath, &matrix_csv_i, tag, x.spectrum, x.A, keep, deduce);
    return(0);
    }
  
/*----------------------------------------------------------------------------
  check admittance availability*/
  for(size_t k = 0; k < s.n; k++) {
    resolve[k]=(keep[k]==1) and (deduce[k]==0);
    }
  int *method=admittance_mts_check(&admittance, s, resolve);
  
//   for(int k=0; k<s.n;k++) {
//     if(s.waves[k].nT==0) {
//       }
//     }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  not necessary to stop process, should be mangaed in spectrum_isAdmittanceSolvable
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(method[0]==ADMITTANCE_UNAVAILABLE) {
    printf("long period admittance mode not available (%d data) (warning)\n",n);
//     printf("long period admittance mode not available (%d data), abort analysis\n",n);
//     save_full_matrix(matrixPath, &matrix_csv_i, tag, x.spectrum, x.A, keep, deduce);
//     return(0);
    }
    
  if(method[1]==ADMITTANCE_UNAVAILABLE) {
    printf("diurnal period admittance mode not available (%d data) (warning)\n",n);
//     printf("diurnal admittance mode not available (%d data), abort analysis\n",n);
//     save_full_matrix(matrixPath, &matrix_csv_i, tag, x.spectrum, x.A, keep, deduce);
//     return(0);
    }
    
  if(method[2]==ADMITTANCE_UNAVAILABLE) {
    printf("semi-diurnal period admittance mode not available (%d data) (warning)\n",n);
//     printf("semi-diurnal admittance mode not available (%d data), abort analysis\n",n);
//     save_full_matrix(matrixPath, &matrix_csv_i, tag, x.spectrum, x.A, keep, deduce);
//     return(0);
    }
  
/*----------------------------------------------------------------------------
  check admittance salvation*/
  spectrum_isAdmittanceSolvable(s, duration, keep, resolve, deduce, log);
/*----------------------------------------------------------------------------
  remove unsolved constituents */
//   spectrum_reduce(s, keep, deduce, stdout);
//
//   x.destroy();
//   status=harmonic_init_new(x, time, n, s, nodal);
  
  for(size_t k = 0; k < s.n; k++) {
    if(keep[k]==0) {
      x.unity(2*k);
      x.unity(2*k+1);
      }
    }
  
/*----------------------------------------------------------------------------
  adjust matrix for admittance*/
  for(size_t k = 0; k < s.n; k++) {
    if(deduce[k] == 1) {
      if(log!=0) {
        fprintf(log,"%s in admittance mode\n",s.waves[k].name);
        }
      x.set_MatrixAdmittance(s.waves[k], k, admittance);
      }
    }
  
  condition[1]=matrix_condition(x.A, x.spectrum, keep, deduce,0);
  
  printf("conditioning : %le, %lf\n", condition[0], condition[1]);
  
  statistic_t stat=get_statistics(serie,n);
  
  for(size_t m = 0; m < n; m++) {
    residuals[m] = serie[m]-stat.mean;
    }
//   for(size_t m = 0; m < n; m++) {
//     residuals[m] = serie[m];
//     }
    
    
  status=x.factorize(0);
  if(status!=0) TRAP_ERR_EXIT(-1, "failure\n");
  
  double *zr = 0, *zi = 0;
  if( (zr = new double[s.n]) == 0) return(0);
  if( (zi = new double[s.n]) == 0) return(0);
  
  status=initialise(s, prior, zr, zi);
  
  int stepmax = 1;
  int nrhs = 1;
  for(size_t step = 0; step < stepmax; step++) {
    x.set_rhs(residuals,b,nrhs,0);
    for(size_t k = 0; k < s.n; k++) {
      if(deduce[k] == 1) {
/*----------------------------------------------------------------------------
        adjust rhs for admittance*/
        x.set_RhsAdmittance(b, zr, zi, s.waves[k], k, admittance);
//         b[2*k]=b[2*k+1]=0;
        }
      if(keep[k] == 0) {
        b[2*k]=b[2*k+1]=0;
        }
      }
    status=x.solve(b, nrhs, 0);
    if(status!=0) TRAP_ERR_EXIT(-1, "failure\n");
    for(size_t k = 0; k < s.n; k++) {
      zr[k] += b[2 * k];
      zi[k] += b[2 * k + 1];
      }
    for(size_t k = 0; k < s.n; k++) {
      for(size_t m = 0; m < n; m++) {
        residuals[m] -= x.cs[k][m] * b[2 * k] + x.sn[k][m] * b[2 * k + 1];
        }
      }
    }
  
/*----------------------------------------------------------------------------
  compute analysis formal error*/
  double *FormalError=0;
  status=x.error(residuals, n, FormalError);
  
  if(status!=0) {
    printf("harmonic_t::error() %d.\n",status);
    save_full_matrix(matrixPath, &matrix_csv_i, tag, x.spectrum, x.A, keep, deduce, b, zr,zi, FormalError);
    return 0;
    }
  
  analysed = new mgr_data_t[s.n];
  
  for(size_t k = 0; k < s.n; k++) {
    mgr_data_t *analysedk=&analysed[k];
    
    fct_xy2polar(zr[k], zi[k], analysedk->amp, analysedk->phi);
    
    analysedk->error = hypot(FormalError[2 * k], FormalError[2 * k + 1]);
    
    if(!isfinite(analysedk->error)){
      printf("wrong error %g for at least wave %s\n",analysedk->error,s.waves[k].name);
      save_full_matrix(matrixPath, &matrix_csv_i, tag, x.spectrum, x.A, keep, deduce, b, zr,zi, FormalError);
      deletep(&analysed);
      break;
      }
    
    analysedk->constituent=s.waves[k];
    }
  
  delete[] FormalError;
  
/*----------------------------------------------------------------------------
  clean-ups */
  x.destroy();
  
  delete[] zr;
  delete[] zi;
  delete[] resolve;
  
  admittance_mts_terminate(admittance);
  
  return(analysed);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mgr_data_t *harmonic_analysis(const double *serie, double *residuals,const double *time, int n,const spectrum_t & prior,double averaging_time, double repetition, spectrum_t & s, int nodal, double maxratio, FILE *log)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  wrapper for the core harmonic_analysis
------------------------------------------------------------------------------*/
{
  spectrum_t work=prior;
  
  int* keep;
  int* deduce;
  hconstant_t atlas;
//   double averaging_time=0.0;
  
//  double repetition=20.86460;
  
  double duration=time[n-1]-time[0];
  
  if(repetition==0) repetition=get_timesampling(time, 1.e+10, n);
  
  if(repetition!=0) {
    spectrum_define_alias(&work, repetition, log);
    spectrum_checkRayleigh(work, duration, 0, 0, log);
    }
  
  mgr_data_t *mgr=harmonic_analysis(serie, residuals, time, n, work, averaging_time, keep, deduce, atlas, nodal, maxratio, log);
  
  delete[] keep;
  delete[] deduce;
  
  s=work;
  
  if(mgr==0) s.n=0;
  
  return(mgr);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mgr_data_t *check_parsing(double *time, int n, const spectrum_t & s, int nodal, FILE *log)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// check parsing with any sampling
/**
\param *time times of the samples
\param n number of samples
\param s prediction spectrum
*/
/*----------------------------------------------------------------------------*/
{
  int status;
  double *serie, *residuals, *tides;
  mgr_data_t *data, *analysed;
  int *keep, *deduce;
  hconstant_t prior(s.n);
  double dtr=M_PI/180.;
  double averaging_time=0.0;

     
  serie=new double[n];
  tides=new double[n];
  residuals=new double[n];
  double maxratio=0.3;
  
  data=new mgr_data_t[s.n];

  for(size_t k = 0; k < s.n; k++) {
    data[k].amp = s.waves[k].Ap;
    data[k].phi = 0.0;
    data[k].constituent=s.waves[k];
    }
  status=harmonic_predictionTS(serie, time, n, data, s.n, nodal);

  prior.s=new spectrum_t;
  (*(prior.s)).duplicate(s);
  
  for(size_t k = 0; k < s.n; k++) {
    float random=sin(k*1.0)*sin(k*1.0);
    if(s.waves[k].Ap==0) prior.a[k] =  random*0.1*0.05;
    else prior.a[k] =  random*0.1*s.waves[k].Ap;
    prior.G[k] = 15.0;
    }
  status=harmonic_predictionTS(tides, time, n, prior, nodal);
  
  for(size_t m = 0; m < n; m++) {
    serie[m] -= tides[m];
    }
  
  analysed=harmonic_analysis(serie, residuals, time, n, s, averaging_time, keep, deduce, prior, nodal, maxratio, log);
  for(size_t k = 0; k < s.n; k++) {
    float random=sin(k*1.0)*sin(k*1.0);
    if(s.waves[k].Ap==0) prior.a[k] =  random*0.1*0.05;
    else prior.a[k] =  random*0.1*s.waves[k].Ap;
    prior.G[k] = 15.0;
    }

  delete[] serie;
  delete[] tides;
  delete[] residuals;
  
  if(analysed==0) goto finish;
  
  for(size_t k = 0; k < s.n; k++) {
    if(keep[k]==0) continue;
    complex< double > z1=polar(data[k].amp, dtr*data[k].phi)/ (double) s.waves[k].Ap;
    complex< double > z2=polar(analysed[k].amp, dtr*analysed[k].phi)/ (double) s.waves[k].Ap;
    double error=100.*abs(z2-z1);
    double a1=data[k].amp/s.waves[k].Ap;
    double G1=data[k].phi;
    double a2=analysed[k].amp/s.waves[k].Ap;
    double G2=analysed[k].phi;
    double T=360./s.waves[k].aliased/24.;
    printf("%6s %3d %6.3lf %6.3lf %6.2lf %6.2lf  error=%6.1lf (%) aliased = %6.1lf days\n",
            s.waves[k].name, deduce[k], a1, G1, a2, G2, error, T);
    }
  delete[] analysed;
  
finish:
  delete[] keep;
  delete[] deduce;
  
  delete[] data;
  
  prior.destroy();
  
  return(analysed);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mgr_data_t *harmonic_analysis_without_parsing(double *serie, double *time, double *residuals,
                                      double *mean, int n, spectrum_t s, int nodal, int *keep, int *deduce, FILE *log)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *-----------------------------------------------------------------------------
 
  spectrum has been already parsed to allow an optimal analysis, distinguishing
  unsolvable, directly solvable and admittance solvable constituents.
  
  time in cnes days
   
-----------------------------------------------------------------------------*/
{
  int status = -1;
  admittance_t admittance;
  int anomaly=0;
  
/* *-----------------------------------------------------------------------------
  remove mean*/
  for(size_t m = 0; m < n; m++) {
    if(abs(serie[m])> 1.e+03) {
      anomaly=1;
      }
    }
  *mean = gsl_stats_mean(serie, 1, n);
  for(size_t m = 0; m < n; m++) {
    residuals[m] = serie[m] - *mean;
    }

  delete[] admittance_mts_check(&admittance, s);
  
  int nsolved = 0;
  for(size_t k = 0; k < s.n; k++){
    if( (keep[k] == 1) || (deduce[k] == 1) ){
      nsolved++;
      }
    }

/* *----------------------------------------------------------------------
  start analysis process */
  int nw = s.n;
  int neq = 2 * s.n;

  double *A = 0, *b = 0;
  if( (A = new double[neq * neq]) == 0) return(0);
  aset(A, neq * neq, 0.);

  if( (b = new double[neq]) == 0) return(0);
  aset(b, neq, 0.);

  double **cs = 0, **sn = 0;
  if( (cs = new double*[nw]) == 0) return(0);
  if( (sn = new double*[nw]) == 0) return(0);
  for(size_t k = 0; k < nw; k++){
    if( (cs[k] = new double[n]) == 0) return(0);
    if( (sn[k] = new double[n]) == 0) return(0);
    }

/* *-----------------------------------------------------------------------------
  compute sinus and cosinus at sampling time, with or without nodal corrections*/
  double f, u, V;
  tidal_wave wave;
  astro_angles_t astro_angles;

  if(nodal == 1){
    for(size_t k = 0; k < s.n; k++) {
      wave = s.waves[k];
      for(size_t m = 0; m < n; m++) {
        init_argument(&astro_angles,poctime_getdatecnes(time[m], 'd'));
        V = greenwhich_argument(astro_angles,wave);
        u = nodal_phase(astro_angles,wave);
        f = nodal_factor(astro_angles,wave.formula);
        cs[k][m] = f * cos(V + u);
        sn[k][m] = f * sin(V + u);
        }
      }
    }
  else{
    for(size_t k = 0; k < s.n; k++) {
      wave = s.waves[k];
      for(size_t m = 0; m < n; m++) {
        init_argument(&astro_angles,poctime_getdatecnes(time[m], 'd'));
        V = greenwhich_argument(astro_angles,wave);
        cs[k][m] = cos(V);
        sn[k][m] = sin(V);
        }
      }
    }

/* *-----------------------------------------------------------------------------
  fill harmonic analysis matrix*/
  for(size_t k = 0; k < s.n; k++) {
    for(size_t l = 0; l < s.n; l++) {
      for(size_t m = 0; m < n; m++) {
        A[2 * l * neq       + 2 * k]     += cs[k][m] * cs[l][m];
        A[2 * l * neq       + 2 * k + 1] += sn[k][m] * cs[l][m];
        A[(2 * l + 1) * neq + 2 * k]     += cs[k][m] * sn[l][m];
        A[(2 * l + 1) * neq + 2 * k + 1] += sn[k][m] * sn[l][m];
        }
      }
    }

//   status=matrix_parsing(A, s, keep, deduce);
  
/* *--------------------------------------------------------------------
  adjust matrix for unresolved/admitted waves*/

  double coef[3] = {1., 1., 1.};

  for(size_t i = 0; i < s.n; i++) {
    if(keep[i] == 0) {
      for(size_t j = 0; j < s.n; j++)
        A[2 * j * neq + 2 * i] = 0.;
      for(size_t j = 0; j < s.n; j++)
        A[2 * j * neq + 2 * i + 1] = 0.;
      for(size_t j = 0; j < s.n; j++)
        A[(2 * j + 1) * neq + 2 * i] = 0.;
      for(size_t j = 0; j < s.n; j++)
        A[(2 * j + 1) * neq + 2 * i + 1] = 0.;
      A[2 * i * neq + 2 * i] = 1.0;
      A[(2 * i + 1) * neq + 2 * i + 1] = 1.0;
      if(log!=0) {
        fprintf(log,"wave: %10s discarded\n", s.waves[i].name);
        }
      }
    else {
      if(log!=0) {
        printf("wave: %10s analysed", s.waves[i].name);
        }

/**--------------------------------------------------------------------

      h=()/potential;

      h-a cos(f-f0) -b sin (f-f0) -c = 0

      a= ri[0][0] h0 + ri[1][0] h1 + ri[2][0] h2
      b= ri[0][1] h0 + ri[1][1] h1 + ri[2][1] h2
      c= ri[0][2] h0 + ri[1][2] h1 + ri[2][2] h2

      h-(ri[0][0] cos(f-f0) + ri[0][1] sin (f-f0) ri[0][2]) h0
       -(ri[1][0] cos(f-f0) + ri[1][1] sin (f-f0) ri[1][2]) h1
       -(ri[2][0] cos(f-f0) + ri[2][1] sin (f-f0) ri[2][2]) h2
       =h- (c0 h0 +c1 h1 + c2 h2)
       =0

 --------------------------------------------------------------------*/

      if(deduce[i] == 1) {
        if(log!=0) {
          fprintf(log," in admittance mode\n");
          }
        for(size_t j = 0; j < neq; j++)
          A[j * neq + 2 * i] = 0.;
        for(size_t j = 0; j < neq; j++)
          A[j * neq + 2 * i + 1] = 0.;
        A[2 * i * neq + 2 * i] = 1. / s.waves[i].Ap;
        A[(2 * i + 1) * neq + 2 * i + 1] = 1. / s.waves[i].Ap;

        size_t k = s.waves[i].nT;
        switch(admittance.method[k]){
          
          case ADMITTANCE_UNAVAILABLE:
            return(0);
            break;

          case ADMITTANCE_SPLINE:
            admittance_mts_sweightP(admittance, s.waves[i], coef);
            for(size_t l = 0; l < 3; l++) {
              size_t j = admittance.windex_s[k][l];
              if(j==-1) TRAP_ERR_EXIT(ENOEXEC,"admittance waves index for nT=%d not initialized\n",k);
              A[2 * j * neq + 2 * i] = -coef[l];
              A[(2 * j + 1) * neq + 2 * i + 1] = -coef[l];
              }
            break;

          case ADMITTANCE_LINEAR:
            admittance_mts_lweightP(admittance, s.waves[i], coef);
            for(size_t l = 0; l < 2; l++) {
              size_t j = admittance.windex_l[k][l];
              if(j==-1) TRAP_ERR_EXIT(ENOEXEC,"admittance waves index for nT=%d not initialized\n",k);
              A[2 * j * neq + 2 * i] = -coef[l];
              A[(2 * j + 1) * neq + 2 * i + 1] = -coef[l];
              }
            break;

          default:
            TRAP_ERR_EXIT(-1, "admittance method is not known\n");
            break;
          }
        }
      else {
        if(log!=0) {
          printf("\n");
          }
        }
      }
    }

  admittance_mts_terminate(admittance);

  int *pivot = 0;
  if( (pivot = new int[neq]) == 0) return(0);

  status=poc_getrf(neq, A, pivot);

  double *zr = 0, *zi = 0;
  if( (zr = new double[nw]) == 0) return(0);
  aset(zr, nw, 0.);
  if( (zi = new double[nw]) == 0) return(0);
  aset(zi, nw, 0.);

  float  *a = 0, *G = 0;
  if( (a = new float[nw]) == 0) return(0);
  aset(a, nw, 0.f);
  if( (G = new float[nw]) == 0) return(0);
  aset(G, nw, 0.f);

  int stepmax = 2;
  int nrhs = 1;
  for(size_t step = 0; step < stepmax; step++) {
    aset(b, neq, 0.);
    for(size_t k = 0; k < s.n; k++) {
      for(size_t m = 0; m < n; m++) {
        b[2 * k]     += cs[k][m] * residuals[m];
        b[2 * k + 1] += sn[k][m] * residuals[m];
        }
      }

/* *--------------------------------------------------------------------
    adjust right-hand side members for unresolved/admitted waves*/
    for(size_t i = 0; i < s.n; i++) {
/*------------------------------------------------------------------------------
      solution forced to zero*/
      if(keep[i] == 0) {
        b[2 * i]     = 0.;
        b[2 * i + 1] = 0.;
        }
/*------------------------------------------------------------------------------
      admittance rhs is zero*/
      if(deduce[i] == 1) {
        b[2 * i]     = 0.;
        b[2 * i + 1] = 0.;
        }
      }
  status=poc_getrs(neq, nrhs, A, pivot, b);

    for(size_t k = 0; k < nw; k++) {
      zr[k] += b[2 * k];
      zi[k] += b[2 * k + 1];
      a[k] = sqrt(zr[k] * zr[k] + zi[k] * zi[k]);
      G[k] = 360.0 * atan2(zi[k],zr[k]) / (2 * M_PI);
      if(G[k] < 0.) G[k] += 360.;
      }
    for(size_t k = 0; k < nw; k++) {
      for(size_t m = 0; m < n; m++) {
        residuals[m] -= cs[k][m] * b[2 * k] + sn[k][m] * b[2 * k + 1];
        }
      }
    } // end iterative analysis

/* *----------------------------------------------------------------------------
  compute amplitudes and phase lags*/
  for(size_t k = 0; k < nw; k++){
    if(keep[k] == 0) {
      a[k] = 0.;
      G[k] = 0.;
      continue;
      }
    a[k] = sqrt(zr[k] * zr[k] + zi[k] * zi[k]);
    G[k] = 360.0 * atan2(zi[k],zr[k]) / (2 * M_PI);
    while(G[k] <   0.) G[k] += 360.;
    while(G[k] > 360.) G[k] -= 360.;
    if(a[k]>10.0) {
      anomaly=1;
      }
    }

/* *----------------------------------------------------------------------------
  compute error on estimates, based on:

  Data Analysis Methods in Physical Oceanography
  by William J.; Thomson Emery
  Pub date: 1998

  -----------------------------------------------------------------------------*/

#if defined(ATLAS) && ATLAS == 1
  double s2 = cblas_ddot(n, residuals, 1, residuals, 1);
#else
  double s2 = 0;
  for(size_t i = 0; i < n; i++){
    s2 += residuals[i] * residuals[i];
    }
#endif
/* *----------------------------------------------------------------------------
  unbiased estimator of the residual variance */
  s2 /= (n - 2 * nsolved - 1);

  vector<double> f_error(0, neq);
  double c = 0.;
  double *mu    = new double[neq];
  double *bck   = new double[neq];
  for(size_t i = 0; i < neq; i++){
    aset(mu, neq, 0.);
/* *----------------------------------------------------------------------------
    commentaires ? */
    mu[i] = 1.f;
    memcpy(bck, mu, neq * sizeof(double));
/* *----------------------------------------------------------------------------
    modified normal equations */
    status=poc_getrs(neq, nrhs, A, pivot, mu);
    for(size_t j = 0; j < neq; j++){
      c += bck[j] * mu[j];
      }
/* *----------------------------------------------------------------------------
    error on i-th element estimates */
    f_error.push_back(sqrt(c * s2));
    }
  delete[] mu;
  delete[] bck;

  
/* *----------------------------------------------------------------------------
  overwrite error for waves in admittance list  */
  status = admittance_mts_error(admittance, s, deduce, f_error);
  if( status != 0){
    return(0);
    }
  
  
/* *----------------------------------------------------------------------------
  store elements */
  mgr_data_t *data;
  data = new mgr_data_t[nw];
  for(size_t k = 0; k < nw; k++) {
    data[k].amp = a[k];
    data[k].phi = G[k];
    data[k].constituent=s.waves[k];
    data[k].error = sqrt( f_error[2 * k] * f_error[2 * k]
                        + f_error[2 * k + 1] * f_error[2 * k + 1] ); // set err = sqrt( s²(i) + s²(i+1) )
    }
  f_error.clear();

  delete[] A;
  delete[] b;
  delete[] pivot;

  for(size_t k = 0; k < nw; k++) delete[] cs[k];
  for(size_t k = 0; k < nw; k++) delete[] sn[k];

  delete[] cs;
  delete[] sn;

  delete[] a;
  delete[] G;

  delete[] zr;
  delete[] zi;

  return(data);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_compute(double *serie, double *residual, int n, harmonic_t x, spectrum_t s, int maxstep, double *a, double *G)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------

  Perform a harmonic analysis with matrix given in arguments. Usefull for
  muliple analysis over a given period

-----------------------------------------------------------------------*/
{
  int k, m, status, step, nw, nrhs=1;
  double *h, *b, bb1, bb2;
  double *zr, *zi, mean, dum;
  double *wcs,*wsn;

  nw = x.neq / 2;
  h = new double[n];
  mean = 0.;
  for(m = 0; m < n; m++) {
    h[m] = serie[m];
    mean += h[m];
    }
  mean /= n;
  for(m = 0; m < n; m++) {
    h[m] -= mean;
    }

//   rms=0.0;
//   for(m = 0; m < n; m++) {
//     rms += h[m]*h[m];
//     }
//   rms=sqrt(rms/n);
//   printf("prior mean=%lf,rms=%lf\n",mean,rms);

  b  = new double[x.neq];
  zr = new double[nw];
  zi = new double[nw];
  for(k = 0; k < nw; k++) {
    zr[k] = 0;
    zi[k] = 0;
    a[k]  = 0;
    G[k]  = 0;
    }
  
  for(step = 0; step < maxstep; step++) {
    for(k = 0; k < x.neq; k++)
      b[k] = 0;
    for(k = 0; k < nw; k++) {
      wcs=x.cs[k];
      wsn=x.sn[k];
      for(m = 0; m < n; m++) {
//         b[2 * k]     += x.cs[k][m] * h[m];
//         b[2 * k + 1] += x.sn[k][m] * h[m];
        b[2 * k]     += wcs[m] * h[m];
        b[2 * k + 1] += wsn[m] * h[m];
        }
      }

    status=poc_getrs(x.neq, nrhs, x.A, x.pivot, b);

//    printf("step %d -----------------------\n",step);
//    if((maxstep>1) || (residual!=0)) {
    for(k = 0; k < nw; k++) {
      zr[k] += b[2 * k];
      zi[k] += b[2 * k + 1];
      }
    
    if((maxstep>1) || (residual!=0)) {
    for(k = 0; k < nw; k++) {
      bb1=b[2 * k];
      bb2=b[2 * k + 1];
      wcs=x.cs[k];
      wsn=x.sn[k];
#ifdef USE_BLAS
      daxpy_ (&n, &bb1, wcs, &unit, h, &unit);
      daxpy_ (&n, &bb2, wsn, &unit, h, &unit);
#else
      for(m = 0; m < n; m++) {
//        h[m] -= x.cs[k][m] * bb1 + x.sn[k][m] * bb2;
        h[m] -= wcs[m] * bb1 + wsn[m] * bb2;
        }
#endif
      }
    dum = 0.;
    for(m = 0; m < n; m++) {
      dum += h[m];
      }
    dum /= n;
    for(m = 0; m < n; m++) {
      h[m] -= dum;
      }
    mean += dum;
    }
    }
  
  for(k = 0; k < nw; k++) {
    a[k] = sqrt(zr[k] * zr[k] + zi[k] * zi[k]);
    G[k] = atan2(zi[k], zr[k]);
//     G[k] = 360.0 * atan2(zi[k], zr[k]) / (2 * pi);
//     if(G[k] < 0.)
//       G[k] += 360.;
//     printf("%s wave: %f %f \n",s.waves[k].name,a[k],G[k]);
    }


//   rms=0.0;
//   for(m = 0; m < n; m++) {
//     rms += h[m]*h[m];
//     }
//   rms=sqrt(rms/n);
//
//   printf("posterior mean=%lf,rms=%lf\n",mean,rms);


  if(residual!=0) {
    for(m = 0; m < n; m++) {
      residual[m] = h[m] + mean;
      }
    }
  zaparr(b);
  zaparr(h);
  zaparr(zr);
  zaparr(zi);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void harmonic_compute(double *serie, double mask,  double *residual, int n, harmonic_t x, spectrum_t s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,status,step,nw,stepmax=2;
  double *h=NULL,*b=NULL;
  int count=0;
  float  *a=NULL,*G=NULL;
  double *zr=NULL,*zi=NULL,mean,dum;

  nw=x.neq/2;

  exitIfNull(
    h=(double *) malloc(n*sizeof(double))
    );
  for(m=0; m<n; m++) h[m]=serie[m];

  for(m=0; m<n; m++) {
    if(h[m]!=mask) count++;
    }

  mean=0.;
  for(m=0; m<n; m++) {
    if(h[m]!=mask) mean+=h[m];
    }
  mean/=count;
  for(m=0; m<n; m++) {
    if(h[m]!=mask) h[m]-=mean;
    }

  exitIfNull(
    b=(double *) malloc(x.neq*sizeof(double))
    );
  exitIfNull(
    zr=(double *) malloc(nw*sizeof(double))
    );
  exitIfNull(
    zi=(double *) malloc(nw*sizeof(double))
    );
  for(k=0; k<nw; k++) zr[k]=0;
  for(k=0; k<nw; k++) zi[k]=0;

  exitIfNull(
    a=(float *) malloc(nw*sizeof(float))
    );
  for(k=0; k<nw; k++) a[k]=0;
  exitIfNull(
    G=(float *) malloc(nw*sizeof(float))
    );
  for(k=0; k<nw; k++) G[k]=0;

  for(k=0; k<nw; k++) b[k]=0;

  for(step=0;step<stepmax;step++) {
    for(k=0; k<x.neq; k++) b[k]=0;
    for(k=0; k < nw; k++) {
      for(m=0; m<n; m++) {
        if(h[m]!=mask) {
          b[2*k]  +=x.cs[k][m]*h[m];
          b[2*k+1]+=x.sn[k][m]*h[m];
          }
        }
      }

    x.nrhs=1;
    status=poc_getrs(x.neq, x.nrhs, x.A, x.pivot, b);

    /*    printf("step %d -----------------------\n",step);*/
    for(k=0; k < nw; k++) {
      zr[k]+=b[2*k];
      zi[k]+=b[2*k+1];
      a[k]=sqrt(zr[k]*zr[k]+zi[k]*zi[k]);
      G[k]=360.0*atan2(zi[k],zr[k])/pi2;
      if(G[k] < 0.) G[k]+=360.;
/*      printf("%s wave: %f %f \n",s.waves[k].name,a[k],G[k]);*/
      }

    for(k=0; k < nw; k++) {
      for(m=0; m<n; m++) {
        if(h[m]!=mask)h[m]-=x.cs[k][m]*b[2*k]+x.sn[k][m]*b[2*k+1];
        }
      }
    dum=0.;
    for(m=0; m<n; m++) {
      if(h[m]!=mask) dum+=h[m];
      }
    dum/=count;
    for(m=0; m<n; m++) {
      if(h[m]!=mask) h[m]-=dum;
      }
    mean+=dum;
    }

  for(k=0; k < nw; k++) {
    serie[2*k]=zr[k];
    serie[2*k+1]=zi[k];
    a[k]=sqrt(zr[k]*zr[k]+zi[k]*zi[k]);
    G[k]=360.0*atan2(zi[k],zr[k])/pi2;
    if(G[k] < 0.) G[k]+=360.;
/*    fprintf(stderr,"%s wave: %5.1f %5.1f \n",s.waves[k].name,a[k],G[k]);*/
    }

  for(m=0; m<n; m++) {
    if(h[m]!=mask) residual[m]=h[m]+mean;
    else           residual[m]=mask;
    }

  free((double *)b);
  free(h);
  free(a);
  free(G);
  free(zr);
  free(zi);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void harmonic_init(double *h, double mask, int n, date_t start, double dt, spectrum_t s, harmonic_t *ptr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,m,status,nw,step=0;
  tidal_wave wave;
  double V,V0,omega,norm;
  double *buffer=NULL;
  harmonic_t x;
  astro_angles_t astro_angles;

  nw=s.n;

  init_argument(&astro_angles,start);

  x.neq=2*s.n;
  x.nframes=n;

  exitIfNull(
    x.cs=(double **) malloc(s.n*sizeof(double))
    );
  exitIfNull(
    x.sn=(double **) malloc(s.n*sizeof(double))
    );

  for(k=0; k < s.n; k++) {
      exitIfNull(
        x.cs[k]=(double *) malloc(n*sizeof(double))
        );
      exitIfNull(
        x.sn[k]=(double *) malloc(n*sizeof(double))
        );
    }

  exitIfNull(
    x.pivot=(int *) malloc(x.neq*sizeof(int))
    );

  for(step=0;step<=0;step++) {
    exitIfNull(
      x.A=(double *) malloc(x.neq*x.neq*sizeof(double))
      );
    for(k=0; k<x.neq; k++)
      for(l=0; l<x.neq; l++) x.A[k*x.neq+l]=0;

    for(k=0; k<nw; k++) {
      wave=s.waves[k];
      omega=wave.omega*dph2rps;
      V0=greenwhich_argument(astro_angles,wave);
      for(m=0; m<n; m++) {
        if(h[m]!=mask) {
          V=(double) m*omega*dt+V0;
          x.cs[k][m]=cos(V);
          x.sn[k][m]=sin(V);
          }
        }
      }

    for(k=0; k<nw; k++) {
      for(l=0; l<nw; l++) {
        for(m=0; m<n; m++) {
          if(h[m]!=mask) {
/*---------- line 2k and 2k+1, column 2l and 2l+1 -----------*/
            x.A[2*l*x.neq    +2*k]  +=x.cs[k][m]*x.cs[l][m];
            x.A[2*l*x.neq    +2*k+1]+=x.sn[k][m]*x.cs[l][m];
            x.A[(2*l+1)*x.neq+2*k]  +=x.cs[k][m]*x.sn[l][m];
            x.A[(2*l+1)*x.neq+2*k+1]+=x.sn[k][m]*x.sn[l][m];
            }
          }
        }
      }

    for(k=0; k<nw; k++) {
      norm=max(s.waves[k].Ap,5.f);
/*---------- line 2k column 2k and line 2k+1 column 2k+1 -----------*/
      x.A[2*k*x.neq    +2*k]  +=(step*5.)*(step*5.)/(norm*norm);
      x.A[(2*k+1)*x.neq+2*k+1]+=(step*5.)*(step*5.)/(norm*norm);
      }

    exitIfNull(
      buffer=(double *) malloc(x.neq*sizeof(double))
      );


    status=poc_getrf(x.neq, x.A, x.pivot);
    }
  *ptr=x;
  free(buffer);

}
