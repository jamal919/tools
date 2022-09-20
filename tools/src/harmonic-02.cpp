
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief definition of mainly obsolete spectrum_t and harmonic operations
Other harmonic operations are defined in harmonic-01.cpp
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "tides.def"
#include "tides.h"


int nsave=0;

double *Aharm=0,**bharm[3];
// spectrum_t AnalysisList;

#include "admittance.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_optimize(spectrum_t spectrum, double duration, double sampling, double *minimal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double tau;
  int i,j,status=0;

/*-----------------------------------------------------------------------------
  check for aliased frequencies */
//   for(i = 0; i < spectrum.n; i++) {
//     if(spectrum.waves[i].omega == 0.0)
//       continue; /*Z0 special case */
//     tau = (1. / spectrum.waves[i].omega) * 360. / 24.;
//     if(fabs(tau) > duration / 3600. / 24.) {
//       printf("wave: %10s, period: %9.3f days \n", spectrum.waves[i].name,tau);
//       status=-1;
//       }
//     }

  *minimal=0;
  
/*-----------------------------------------------------------------------------
  check for unresolved frequencies */
  for(i = 0; i < spectrum.n; i++) {
    if(spectrum.waves[i].omega == 0.0)
      continue; /*Z0 special case */
    tau = (1. / spectrum.waves[i].omega) * 360. / 24.;
    if(fabs(tau) > duration / 3600. / 24.) {
      printf("wave: %10s, period: %9.3f days \n", spectrum.waves[i].name,fabs(tau));
      }
    updatemax(minimal,fabs(tau)*d2s);
    }

/*-----------------------------------------------------------------------------
  check for redundant frequencies and non-separable constituents*/
  for(i = 0; i < spectrum.n; i++) {
    if(spectrum.waves[i].omega == 0.0)
      continue; /*Z0 special case */
    for(j = i + 1; j < spectrum.n; j++) {
      if(spectrum.waves[j].omega == 0.0)
        continue; /*Z0 special case */
      tau = spectrum.waves[j].omega - spectrum.waves[i].omega;
      if(tau == 0.) {
        printf("wave: %10s %10s, separation IMPOSSIBLE!!!\n",
          spectrum.waves[i].name, spectrum.waves[j].name);
        }
      tau = (1. / tau) * 360. / 24.;
      if(fabs(tau) > 2 * duration / 3600. / 24.) {
        printf("wave: %10s %10s, separation: %9.3f days \n",
          spectrum.waves[i].name, spectrum.waves[j].name, fabs(tau));
        }
      updatemax(minimal,fabs(tau)*d2s);
      }
    }
  
  *minimal *=1.1;
  updatemax(minimal,5.0*d2s);
  printf("optimal analysis duration : %9.3lf days\n",*minimal / 3600. / 24.);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_check(spectrum_t spectrum, double duration, double sampling)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double tau;
  int i,j,status=0;

/*-----------------------------------------------------------------------------
  check for aliased frequencies */
//   for(i = 0; i < spectrum.n; i++) {
//     if(spectrum.waves[i].omega == 0.0)
//       continue; /*Z0 special case */
//     tau = (1. / spectrum.waves[i].omega) * 360. / 24.;
//     if(fabs(tau) > duration / 3600. / 24.) {
//       printf("wave: %10s, period: %9.3f days \n", spectrum.waves[i].name,tau);
//       status=-1;
//       }
//     }
  
/*-----------------------------------------------------------------------------
  check for unresolved frequencies */
  for(i = 0; i < spectrum.n; i++) {
    if(spectrum.waves[i].omega == 0.0)
      continue; /*Z0 special case */
    tau = (1. / spectrum.waves[i].omega) * 360. / 24.;
    if(fabs(tau) > duration / 3600. / 24.) {
      printf("wave: %10s, period: %9.3f days \n", spectrum.waves[i].name,fabs(tau));
      status=-1;
      }
    }

/*-----------------------------------------------------------------------------
  check for redundant frequencies and non-separable constituents*/
  for(i = 0; i < spectrum.n; i++) {
    if(spectrum.waves[i].omega == 0.0)
      continue; /*Z0 special case */
    for(j = i + 1; j < spectrum.n; j++) {
      if(spectrum.waves[j].omega == 0.0)
        continue; /*Z0 special case */
      tau = spectrum.waves[j].omega - spectrum.waves[i].omega;
      if(tau == 0.) {
        printf("wave: %10s %10s, separation IMPOSSIBLE!!!\n",
          spectrum.waves[i].name, spectrum.waves[j].name);
        }
      tau = (1. / tau) * 360. / 24.;
      if(fabs(tau) > 2 * duration / 3600. / 24.) {
        printf("wave: %10s %10s, separation: %9.3f days \n", spectrum.waves[i].name, spectrum.waves[j].name, fabs(tau));
        status=-1;
        }
      }
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int addwave(spectrum_t *list, tidal_wave wave)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,nadd, k, n, nwave, nwmax;
  size_t size;

  size=sizeof(tidal_wave);

  n=list->n+1;
  nwmax=list->nmax;

  if(n > nwmax) {
    nwmax++;
    exitIfNull(
      list->waves=(tidal_wave *) realloc(list->waves,nwmax*size)
      );
    }
  memcpy( &(list->waves[n-1]), &wave, size);

  list->nmax=nwmax;
  list->n=n;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int removewave(spectrum_t *list, int i)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nAfter;
  list->n--;
  nAfter=list->n-i;
  memcpy( &(list->waves[i]), &(list->waves[i+1]), nAfter*sizeof(tidal_wave));

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int spectrum_init(spectrum_t *list)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,nadd, k, n, nwave, nwmax;
  char *in=NULL,*out=NULL;
  size_t size;

  list->n=0;
  list->nmax=10;

  size=sizeof(tidal_wave);
  exitIfNull(list->waves=new tidal_wave[list->nmax]);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void spectrum_terminate(harmonic_t harmonic, spectrum_t list)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  allocate memory and set to zero the matrix A and RHS vector b for
  harmonic analysis purposes
----------------------------------------------------------------------*/
{
  int   i,k,n,neq;

  neq=2*list.n;

  for (i=0; i<list.n; i++) {
    if(list.waves[i].omega == 0.) {
      list.waves[i].init();
      }
    printf ("wave: %10s, pulsation: %12.6f degrees/h \n", list.waves[i].name,list.waves[i].omega);
    }

  harmonic.A =new double[neq*neq];
  for(n=0;n<neq*neq;n++) harmonic.A[n]=0.0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void spectrum_terminate(spectrum_t list)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  allocate memory and set to zero the matrix A and RHS vector b for
  harmonic analysis purposes
----------------------------------------------------------------------*/
{
  int   i,k,n,neq;

  neq=2*list.n;

  for (i=0; i<list.n; i++) {
    if(list.waves[i].omega == 0.) {
      list.waves[i].init();
      }
    printf ("wave: %10s, pulsation: %12.6f degrees/h \n", list.waves[i].name,list.waves[i].omega);
    }

  Aharm =new double[neq*neq];
  for(n=0;n<neq*neq;n++) Aharm[n]=0.0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void harmonic_init(spectrum_t WaveList, spectrum_t  & AnalysisList, int nndes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  allocate memory and set to zero the matrix A and RHS vector b for
  harmonic analysis purposes
----------------------------------------------------------------------*/
{
  int   i,k,n,neq;
  int   nadd;

  AnalysisList.n=WaveList.n;

  exitIfNull(AnalysisList.waves=new tidal_wave[AnalysisList.n]);
  for(k=0; k < WaveList.n; k++) {
    AnalysisList.waves[k]=WaveList.waves[k];
    }

  neq=2*AnalysisList.n;

  for (i=0; i<AnalysisList.n; i++) {
    if(AnalysisList.waves[i].omega == 0.) {
      AnalysisList.waves[i].init();
      }
    printf ("wave: %10s, pulsation: %12.6f degrees/h \n", AnalysisList.waves[i].name,AnalysisList.waves[i].omega);
    }

  Aharm =new double[neq*neq];

  for(i=0;i<3;i++) {
    exitIfNull(
      bharm[i]=(double **) malloc(nndes*sizeof(double *))
      );
    for(n=0; n<nndes; n++) {
      exitIfNull(
        bharm[i][n]=(double *) malloc(neq*sizeof(double))
        );
      for(k=0; k < neq; k++) {
        bharm[i][n][k]=0;
        }
      }
    }

  for(k=0; k < neq*neq; k++) {
    Aharm[k]  =0;
    }

  printf("#harmonic_start OK ...\n");
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void harmonic_save(spectrum_t  & AnalysisList, double **b[3], int nndes, date_t start,date_t final)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  save the harmonic constants analysed at mesh nodes

----------------------------------------------------------------------*/
{
  int   k,n;
  float a,G,a1,G1,a2,G2;
  tidal_wave wave;
  double zi,zr;
  char filename[1024];
  char s_start[1024],s_final[1024];
  FILE *out=NULL;
  char *sdate=NULL;
  date_t actual;

  sprintf(s_start,"%4.4d-%2.2d-%2.2d",start.year,start.month,start.day);
  sprintf(s_final,"%4.4d-%2.2d-%2.2d",final.year,final.month,final.day);


  k=0;

  wave=AnalysisList.waves[k];
  sprintf(filename,"%s.ele.%s-%s.s2r",wave.name,s_start,s_final);
  out=fopen(filename,"w");
  if(out==0) {
      TRAP_ERR_EXIT(-1,"file opening issue : %s \n",filename);
    }

  fprintf(out,"analysed from mOG2D archives [%s to %s]\n",s_start,s_final);
  fprintf(out,"mean elevation (m)\n");
  fprintf(out,"%d nodes\n", nndes);
  for(n=0; n<nndes; n++) {
      zr=b[0][n][2*k];
/*-------------------- elevation in meters ------------------------*/
      fprintf(out,"%d %f\n",n,zr);
      }
  fclose(out);

/*-------------------- Tidal currents ------------------------*/
  sprintf(filename,"%s.uv.%s-%s.v2r",wave.name,s_start,s_final);
  out=fopen(filename,"w");
  if(out==0) {
      TRAP_ERR_EXIT(-1,"file opening issue : %s \n",filename);
    }
  fprintf(out,"analysed from mOG2D archives [%s to %s]\n",s_start,s_final);
  fprintf(out,"mean current (m/s)\n");
  fprintf(out,"%d nodes\n", nndes);
  for(n=0; n<nndes; n++) {
      zr=b[1][n][2*k];
/*-------------------- Amplitude in meters/s-----------------------*/
      a1=zr;
      zr=b[2][n][2*k];
      a2=zr;
      fprintf(out,"%d %f %f \n",n,a1,a2);
      }
  fclose(out);


  for(k=0; k < AnalysisList.n; k++) {
    wave=AnalysisList.waves[k];
/*-------------------- Tidal elevations ------------------------*/
    sprintf(filename,"%s.ele.%s-%s.s2c",wave.name,s_start,s_final);
    out=fopen(filename,"w");
    if(out==0) {
        TRAP_ERR_EXIT(-1,"file opening issue : %s \n",filename);
      }
    fprintf(out,"analysed from mOG2D archives [%s to %s]\n",s_start,s_final);
    fprintf(out,"harmonic constants (m)\n");
    fprintf(out,"%d nodes\n", nndes);
    for(n=0; n<nndes; n++) {
      zr=b[0][n][2*k];
      zi=b[0][n][2*k+1];
/*-------------------- Amplitude in meters ------------------------*/
      a=sqrt(zr*zr+zi*zi);
      G=atan2(zi,zr)*r2d;
      fprintf(out,"%d %f %f\n",n,a,G);
      }
    fclose(out);
/*-------------------- Tidal currents ------------------------*/
    sprintf(filename,"%s.uv.%s-%s.v2c",wave.name,s_start,s_final);
    out=fopen(filename,"w");
    if(out==0) {
        TRAP_ERR_EXIT(-1,"file opening issue : %s \n",filename);
      }
    fprintf(out,"analysed from mOG2D archives [%s to %s]\n",s_start,s_final);
    fprintf(out,"harmonic constants (m/s)\n");
    fprintf(out,"%d nodes\n", nndes);
    for(n=0; n<nndes; n++) {
      zr=b[1][n][2*k];
      zi=b[1][n][2*k+1];
/*-------------------- Amplitude in meters/s-----------------------*/
      a1=sqrt(zr*zr+zi*zi);
      G1=atan2(zi,zr)*r2d;
      zr=b[2][n][2*k];
      zi=b[2][n][2*k+1];
      a2=sqrt(zr*zr+zi*zi);
      G2=atan2(zi,zr)*r2d;
      fprintf(out,"%d %f %f %f %f\n",n,a1,G1,a2,G2);
      }
    fclose(out);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void harmonic_end(spectrum_t  & AnalysisList, date_t start,date_t final,int nndes, int count)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  perform the harmonic analysis at mesh nodes

----------------------------------------------------------------------*/
{
  int  k,l,n,neq=2*AnalysisList.n;
  int  status,step,nrhs=1;
  int *pivot=NULL;
  double *A=NULL,**b[3],*tmp=NULL,*x=NULL,*y=NULL;

  printf("#harmonic analysis, nnodes=%d, neq=%d ...\n",nndes,neq);
  exitIfNull(
    pivot=(int *)  malloc(neq*sizeof(int))
    );
  exitIfNull(
    A=(double *)  malloc(neq*neq*sizeof(double))
    );
  exitIfNull(
    b[0]=(double **) malloc(nndes*sizeof(double *))
    );
  exitIfNull(
    b[1]=(double **) malloc(nndes*sizeof(double *))
    );
  exitIfNull(
    b[2]=(double **) malloc(nndes*sizeof(double *))
    );

  for(n=0; n<nndes; n++) {
    exitIfNull(
      b[0][n]=(double *) malloc(neq*sizeof(double))
      );
    exitIfNull(
      b[1][n]=(double *) malloc(neq*sizeof(double))
      );
    exitIfNull(
      b[2][n]=(double *) malloc(neq*sizeof(double))
      );
    for(k=0; k < neq; k++) {
      b[0][n][k]=bharm[0][n][k];
      b[1][n][k]=bharm[1][n][k];
      b[2][n][k]=bharm[2][n][k];
      }
    }

  for(k=0; k < neq*neq; k++) {
    A[k]=Aharm[k];
    }

  for(k=0; k < neq; k++) {
    A[k*neq+1]=0;
    }

  A[neq+1]=1.0;

  status=poc_getrf(neq, A, pivot);

  for(n=0; n<nndes; n++) {
    for(step=0;step<1;step++) {
      for(k=0;k<3;k++) status=poc_getrs(neq, nrhs, A, pivot, b[k][n]);
      }
    }

  printf("#resolution OK ...\n");

  harmonic_save(AnalysisList, b,nndes, start, final);

  free(pivot);
  free(A);
  for(n=0; n<nndes; n++){
    free(b[0][n]);
    free(b[1][n]);
    free(b[2][n]);
    }
  free(b[0]);
  free(b[1]);
  free(b[2]);

  free(Aharm);
  for(n=0; n<nndes; n++) {
    free(bharm[0][n]);
    free(bharm[1][n]);
    free(bharm[2][n]);
    }
  free(bharm[0]);
  free(bharm[1]);
  free(bharm[2]);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void harmonic_storage_obsolete(spectrum_t  & AnalysisList, double t,int nndes,float *hmean,float *buffer[3],int *count, const astro_angles_t &astro_angles)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

  update the harmonic matrix and righthand side vector
  \sa harmonic_storage(double,harmonic_t,int,float**)

-----------------------------------------------------------------------------*/
{
  int   k,l,n,neq=2*AnalysisList.n;
  float a,G;
  tidal_wave wave;
  double V,V0,omega;
  double h,*cs=NULL,*sn=NULL,u,v;

  if(AnalysisList.n==0) return;

  for(n=0; n<nndes; n++){
    h=buffer[0][n];
    hmean[n]+=h;
    }

  (*count)++;

  exitIfNull(
    cs=(double *) malloc(AnalysisList.n*sizeof(double))
    );
  exitIfNull(
    sn=(double *) malloc(AnalysisList.n*sizeof(double))
    );

  harmonic_coefficients(t, AnalysisList, cs, sn, 0, astro_angles);

  for(k=0; k < AnalysisList.n; k++) {
    for(n=0; n<nndes; n++) {
      h=buffer[0][n];
      u=buffer[1][n];
      v=buffer[2][n];
      bharm[0][n][2*k]  +=h*cs[k];
      bharm[0][n][2*k+1]+=h*sn[k];
      bharm[1][n][2*k]  +=u*cs[k];
      bharm[1][n][2*k+1]+=u*sn[k];
      bharm[2][n][2*k]  +=v*cs[k];
      bharm[2][n][2*k+1]+=v*sn[k];
      }

/*-----------------------------------------------------------------------------
 line 2k is derivation with respect to real part
 line 2k+1 is derivation with respect to imaginary part
 column 2l is coefficient with respect to real part
 column 2l+1 is derivation with respect to imaginary part
------------------------------------------------------------------------------*/
    for(l=0; l < AnalysisList.n; l++) {
      Aharm[2*l*neq    +2*k]  +=cs[k]*cs[l];
      Aharm[2*l*neq    +2*k+1]+=sn[k]*cs[l];
      Aharm[(2*l+1)*neq+2*k]  +=cs[k]*sn[l];
      Aharm[(2*l+1)*neq+2*k+1]+=sn[k]*sn[l];
      }
    }
  free(cs);
  free(sn);

}

