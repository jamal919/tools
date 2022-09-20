
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief definition of most harmonic operations
Other harmonic operations are defined in harmonic-02.cpp
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "functions.h"
#include "zapper.h"

#include "tides.def"
#include "constants.h"
#include "tides.h"
#include "poc-time.h"

#define nspecies 3
#include "admittance-mts.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int harmonic_coefficients_template(double t, spectrum_t spectrum, T *cs, T *sn, int nodal,const astro_angles_t &astro_angles)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** Compute tidal astronomical phase (sinus and cosinus)
at the given time.

\param t time in seconds \bug The reference must be the same for each related calls of this function and for the prior call to init_argument().
\param spectrum waves to compute the phase for
\param cs <b> &lt;float or double&gt;[spectrum.n] </b> array of cosines
\param sn <b> &lt;float or double&gt;[spectrum.n] </b> array of sines
\param nodal whether nodal corrections are carried out
\returns 0

\sa harmonic_correction(), harmonic_storage(double,harmonic_t,int,float**) and  harmonic_storage(double,int,float*,float*[3],int*)

It fills \b cs and \b sn as the real and imaginary parts of \f$ \nu_n e^{j\left(\omega_n t + \phi_n\right)} \f$,
with \f$ \nu_n \f$ the complex nodal correction,
\f$ \omega_n \f$ the wave pulsation and
\f$ \phi_n \f$ the Geenwich argument.
*/
/*----------------------------------------------------------------------------*/
{
  int j;
  double f;
  double V, V0, u, omega;
  tidal_wave wave;

  for(j = 0; j < spectrum.n; j++) {
    wave = spectrum.waves[j];
    V0 = greenwhich_argument(astro_angles,wave);
    omega = wave.omega * dph2rps;
    V = omega * t + V0;
    if(nodal == 1) {
      u = nodal_phase(astro_angles,wave);
      f = nodal_factor(astro_angles,wave.formula);
/*       printf(" %6s %9.3f %9.3f \n",wave.name,u,f); */
      cs[j] = f * cos(V + u);
      sn[j] = f * sin(V + u);
    } else {
      cs[j] = cos(V);
      sn[j] = sin(V);
    }
  }
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_coefficients(double t, spectrum_t spectrum, float *cs, float *sn, int nodal, const astro_angles_t &astro_angles)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** Compute tidal astronomical phase (sinus and cosinus)
See harmonic_coefficients_template().
*/
/*----------------------------------------------------------------------------*/
{return harmonic_coefficients_template(t,spectrum,cs,sn,nodal,astro_angles);}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_coefficients(double t, spectrum_t spectrum, double *cs, double *sn, int nodal, const astro_angles_t &astro_angles)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** Compute tidal astronomical phase (sinus and cosinus)
See harmonic_coefficients_template().
*/
/*----------------------------------------------------------------------------*/
{return harmonic_coefficients_template(t,spectrum,cs,sn,nodal,astro_angles);}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double deltaOmega2separation(double deltaOmega)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/**
\param deltaOmega difference of pulsation (deg/h). MUST BE NON-ZERO!
\returns separation in days :
\code 1. / fabs(deltaOmega) * 360. / 24. \endcode
*/
/*----------------------------------------------------------------------------*/
{
  return 1. / fabs(deltaOmega) * 360. / 24.;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_coefficients(const double *time, int n,const spectrum_t & spectrum, harmonic_t & x, int nodal, double averaging_time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
 
  Compute tidal astronomical phase (sinus and cosinus)

------------------------------------------------------------------------------*/
{
  astro_angles_t *astro_angles;
  bool averaged=(averaging_time!=0.0);
  double omega;
  
  astro_angles =new astro_angles_t[n];
  
#pragma omp parallel for
  for(int m=0; m<n; m++) {
    date_t start=poctime_getdatecnes(time[m], 'd');
    init_argument(&astro_angles[m],start);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  compute time-dependent coefficients:
  
  instantaneous observations:
  ---------------------------
  
    c = f * cos(V + u)
    s = f * sin(V + u)

  time-averaged observations:
  ---------------------------
  
    c =  f * [sin(V + u) (t+dT) - sin(V + u) (t)] / (omega*dT)
    
    s = -f * [cos(V + u) (t+dT) - cos(V + u) (t)] / (omega*dT)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
    
  for(int k=0; k < spectrum.n; k++) {
    tidal_wave wave;
    wave=spectrum.waves[k];
    omega=wave.omega*dph2rps;
#pragma omp parallel for
    for(int m=0; m<n; m++) {
      double f,u;
      double V0=greenwhich_argument(astro_angles[m],wave);
      if(nodal==1) {
        u=nodal_phase(astro_angles[m],wave);
        f=nodal_factor(astro_angles[m],wave.formula);
        }
      else {
        u=0.;
        f=1.;
        }
      if(averaged) {
        double V=V0 + u;
        x.cs[k][m] =  f * sin(V);
        x.sn[k][m] = -f * cos(V);
        V=V0+ u - omega*averaging_time;
        x.cs[k][m] -=  f * sin(V);
        x.sn[k][m] -= -f * cos(V);
        x.cs[k][m] /= omega*averaging_time;
        x.sn[k][m] /= omega*averaging_time;
        }
      else {
        double V=V0+u;
        x.cs[k][m] = f * cos(V);
        x.sn[k][m] = f * sin(V);
        }
      }
    }
  
  delete[] astro_angles;

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double printAndGetSpectrumDetails(spectrum_t AnalysisList)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints the critical latitude, pulsation and separations of tidal waves

\param AnalysisList a pre-initialised spectrum_t
\return minimum number of days for separation
*/
/*----------------------------------------------------------------------------*/
{
  int i, j;//wave indexes
  const double two_Omega_deg_h=two_Omega/dph2rps;//in deg/h
  double omegai,
    critical,//< critical latitude in degrees
    tau;//<separation in deg/h or in days
  double maxNeeded=-1.;

/*------------------------------------------------------------------------------
  ensure all AnalysisList.waves[].omega are initialised */
  AnalysisList.waves_init();
  
/*------------------------------------------------------------------------------
  wave list */
  printf("Wave list for analysis: %d \n", AnalysisList.n);
  for(i = 0; i < AnalysisList.n; i++) {
    omegai=AnalysisList.waves[i].omega;
    printf("wave: %10s, pulsation: %12.6f degrees/h",AnalysisList.waves[i].name, omegai);
    if(omegai <= two_Omega_deg_h){
      critical = asin(omegai / two_Omega_deg_h)*r2d;
      printf(", critical latitude: %6.1f degrees",critical);
      }
    printf("\n");
    }
  
  printf("\n");
  
/*------------------------------------------------------------------------------
  separation */
  printf("Check separation for analysis:\n");
  for(i = 0; i < AnalysisList.n; i++) {
    for(j = i + 1; j < AnalysisList.n; j++) {
      tau = AnalysisList.waves[j].omega - AnalysisList.waves[i].omega;
      if(tau == 0.) {
        printf("wave: %10s %10s, separation IMPOSSIBLE!!!\n",
                AnalysisList.waves[i].name, AnalysisList.waves[j].name);
        maxNeeded=INFINITY;
        }
      else {
        tau = deltaOmega2separation(tau);
        updatemax(&maxNeeded,tau);
        if(tau > 14.)
          printf("wave: %10s %10s, separation: %9.3f days \n",
                  AnalysisList.waves[i].name, AnalysisList.waves[j].name, tau);
        }
      }
    }
  
  return maxNeeded;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_start(harmonic_t *harmonic,const spectrum_t AnalysisList)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief allocates and zeroes the harmonic matrix

\param *harmonic the structure to modify
\param AnalysisList a pre-initialised spectrum_t
\returns 0 if success, -1 if AnalysisList is not initialised

\sa harmonic_start(spectrum_t,int,int,int)
*/
/*----------------------------------------------------------------------------*/
{
  int k,neq,neq2;

  neq = 2 * AnalysisList.n;
  if(neq<=0){
    harmonic->A=NULL;
    return -1;
    }
  neq2=neq*neq;
  harmonic->spectrum=AnalysisList;
  harmonic->neq   =neq;
  harmonic->A = new double[neq2];
  for(k = 0; k < neq2; k++) {
    harmonic->A[k] = 0;
    }
  
  ///It also sets to NULL cs, sn, pivot and mask and to 0 nrhs, nndes and nframes
  harmonic->cs=NULL;
  harmonic->sn=NULL;
  harmonic->pivot=NULL;
  harmonic->mask=NULL;
  harmonic->nrhs=0;
  harmonic->nndes=0;
  harmonic->nframes=0;
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_start(harmonic_t *harmonic,const size_t *nndes, int nrhs)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief allocates and zeroes the rhs

\param *harmonic the structure to modify, with an already initialised neq
\param nndes number of points in the grid
\param nrhs number of variables. See the note about harmonic_t.
\returns 0 if success, -1 if neq is not initialised or if nndes<=0 or if nrhs<=0

\sa harmonic_start(spectrum_t,int,int,int)
*/
/*----------------------------------------------------------------------------*/
{
  int i,n,k;
  int neq=harmonic->neq;
  if(neq<=0 || nndes<=0 || nrhs<=0)
    return -1;

  harmonic->nrhs=nrhs;
  exitIfNull(harmonic->rhs = new double **[nrhs]);
  for(i = 0; i < nrhs; i++) {
    exitIfNull(harmonic->rhs[i] = new double *[nndes[i]]);
    for(n = 0; n < nndes[i]; n++) {
      exitIfNull(harmonic->rhs[i][n] = new double[neq]);
      for(k = 0; k < harmonic->neq; k++) {
        harmonic->rhs[i][n][k] = 0;
        }
      }
    }
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_start(harmonic_t *harmonic,int nframe)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief allocates the sines and cosines

\param *harmonic the structure to modify, with an already initialised spectrum
\param nframes number of frames. Default:0. If specified, allocates sines and cosines.
\returns 0 if success, -1 if the spectrum is not initialised or if nframe<=0

\sa harmonic_start(spectrum_t,int,int,int)
*/
/*----------------------------------------------------------------------------*/
{
  int n,nw;
  nw=harmonic->spectrum.n;
  if(nw<=0 || nframe<=0)
    return -1;
  
  harmonic->nframes=nframe;
  exitIfNull(harmonic->cs=new double*[nw]);
  exitIfNull(harmonic->sn=new double*[nw]);
  for(n = 0; n < nw; n++) {
    exitIfNull(harmonic->cs[n]=new double[nframe]);
    exitIfNull(harmonic->sn[n]=new double[nframe]);
    }
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  harmonic_t harmonic_start(spectrum_t AnalysisList, const size_t *nndes, int nrhs, int nframe)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief returns a partially allocated and zeroed harmonic_t structure

\param AnalysisList a pre-initialised spectrum_t
\param *nndes array[nrhs] number of points in the grid
\param nrhs number of variables. See the note about harmonic_t.
\param nframes number of frames. Default:0. If specified, allocates sines and cosines.
\returns the partially allocated and zeroed harmonic_t structure

\sa harmonic_storage()
*/
/*----------------------------------------------------------------------------*/
{
  int status;
  harmonic_t harmonic;
  
  /// It allocates and zeroes the harmonic matrix with harmonic_start(harmonic_t*,spectrum_t)
  status=harmonic_start(&harmonic,AnalysisList);
  if(status!=0){
    STDERR_BASE_LINE("error with harmonic_start() because AnalysisList.n=%d<=0\n",AnalysisList.n);
    goto end;
    }
  
  harmonic.nndes=copy(nndes,nrhs);
  
  /// It allocates and zeroes the rhs with harmonic_start(harmonic_t*,int,int)
  harmonic_start(&harmonic,nndes,nrhs);
  
  if(nframe){
    ///If required to do so, it allocates the sines and cosines with harmonic_start(harmonic_t*,int)
    harmonic_start(&harmonic,nframe);
    }
  
end:
  return(harmonic);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 extern "C" void harmonic_start_(harmonic_t *harmonic, spectrum_t *spectrum, size_t *nndes, int *nframe)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Fortran wrapper for harmonic_start(spectrum_t AnalysisList,const int *nndes, int nrhs, int nframe)
/*----------------------------------------------------------------------------*/
{
  *harmonic=harmonic_start(*spectrum,nndes,1,*nframe);
  if(harmonic->A==NULL) TRAP_ERR_EXIT(-1,"error with harmonic_start\n");
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_init(harmonic_t *x,double *time, int nodal,const astro_angles_t &astro_angles, double averaging_time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief initialises the harmonic matrix

\date 2011-09-16 Damien ALLAIN : creation from harmonic_init01()

\param *x the structure to modify, with an already initialised spectrum. See harmonic_start().
\param *times CNES times of frames in seconds
\param nodal whether nodal corrections are carried out
\param astro_angles astronomic angles
\returns -1 if \c x->A is NULL, 0 otherwise
*/
/*----------------------------------------------------------------------------*/
{
  date_t t_schureman(1900,1,1,0.0);               //time reference for astronomical abaques
  date_t t_CNES(1950,1,1,0.0);                    //time reference for CNES times
  double t0=(astro_angles.t_julian*jc2d+julian_day(t_schureman)-julian_day(t_CNES))*24.*3600.;//time in seconds since 1950-01-01 00:00
  int k, l, m;
  tidal_wave wave;
  double f,u,V, V0, omega;
  // double mask;
  bool averaged=(averaging_time!=0.0);
  
//   printf("%s : initialise harmonic matrix (%d frames)\n",__func__,x->nframes);
  
  int nframe=x->nframes;
  int nw=x->spectrum.n;//number of waves
  
  /// \bug 2011-09-16 Damien ALLAIN : this does not check for proper initialisation
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  compute time-dependent coefficients:
  
  instantaneous observations:
  ---------------------------
  
    c = f * cos(V + u)
    s = f * sin(V + u)

  time-averaged observations:
  ---------------------------
  
    c =  f * [sin(V + u) (t+dT) - sin(V + u) (t)] / (omega*dT)
    
    s = -f * [cos(V + u) (t+dT) - cos(V + u) (t)] / (omega*dT)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(k=0; k<nw; k++) {
    wave=x->spectrum.waves[k];
    omega=wave.omega*dph2rps;
    V0=greenwhich_argument(astro_angles,wave);
    if(nodal == 1) {
      u = nodal_phase (astro_angles,wave);
      f = nodal_factor(astro_angles,wave.formula);
      }
    else {
      u = 0;
      f = 1.;
      }
    if(averaged and omega!=0) {
      for(m=0; m<nframe; m++) {
//         if(time[m]!=mask) continue;
        V=(double) omega*(time[m]-t0)+V0;
        x->cs[k][m] =  f * sin(V + u);
        x->sn[k][m] = -f * cos(V + u);
        V=(double) omega*(time[m]-averaging_time-t0)+V0;
        x->cs[k][m] -=  f * sin(V + u);
        x->sn[k][m] -= -f * cos(V + u);
        x->cs[k][m] /= omega*averaging_time;
        x->sn[k][m] /= omega*averaging_time;
        }
      }
    else {
      for(m=0; m<nframe; m++) {
//         if(time[m]!=mask) continue;
        V=(double) omega*(time[m]-t0)+V0;
        x->cs[k][m] = f * cos(V + u);
        x->sn[k][m] = f * sin(V + u);
        }
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  assembly matrix

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(x->A==0)
    return -1;
  
  for(k = 0; k < nw; k++) {
    for(l = 0; l < nw; l++) {
      for(m = 0; m < nframe; m++) {
/*-----------------------------------------------------------------------------
        row 2k   is derivation with respect to real part
        row 2k+1 is derivation with respect to imaginary part
        col 2l   is coefficient with respect to real part
        col 2l+1 is derivation with respect to imaginary part
------------------------------------------------------------------------------*/
        x->A[2 * l * x->neq + 2 * k]           += x->cs[k][m] * x->cs[l][m];
        x->A[2 * l * x->neq + 2 * k + 1]       += x->sn[k][m] * x->cs[l][m];
        x->A[(2 * l + 1) * x->neq + 2 * k]     += x->cs[k][m] * x->sn[l][m];
        x->A[(2 * l + 1) * x->neq + 2 * k + 1] += x->sn[k][m] * x->sn[l][m];
        }
      }
/*-----------------------------------------------------------------------------
    permanent tide Z0, imaginary parts needs a special treatment*/
    if(x->spectrum.waves[k].omega == 0.0)
      x->A[(2 * k + 1) * x->neq + 2 * k + 1] = 1.0;
    }
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 extern "C" void harmonic_init_(harmonic_t *harmonic, double *time,int *nodal,astro_angles_t *astro_angles, double averaging_time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  Fortran wrapper for harmonic_init
------------------------------------------------------------------------------*/
{
  harmonic_init(harmonic,time,*nodal,*astro_angles, averaging_time);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_init_new(harmonic_t & x,const double *time, int nframes,const spectrum_t & spectrum, int nodal, double averaging_time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int k;
  
  x.spectrum=spectrum;
  x.nframes=nframes;
  
  int nw=x.spectrum.n;
  
  status=x.allocate_coefficients();
  if(status!=0) TRAP_ERR_EXIT(-1, "failure\n");
  
  status=harmonic_coefficients(time, nframes, x.spectrum, x, nodal, averaging_time);
  if(status!=0) TRAP_ERR_EXIT(-1, "failure\n");
  
  x.allocate_matrix();
  
  for(k = 0; k < nw; k++) {
#pragma omp parallel for
    for(int l = 0; l < nw; l++) {
      for(int m = 0; m < nframes; m++) {
/*-----------------------------------------------------------------------------
        line 2k   is derivation with respect to real part
        line 2k+1 is derivation with respect to imaginary part
        column 2l   is coefficient with respect to real part
        column 2l+1 is coefficient with respect to imaginary part
------------------------------------------------------------------------------*/
        x.A[2 * l * x.neq + 2 * k]           += x.cs[k][m] * x.cs[l][m];
        x.A[2 * l * x.neq + 2 * k + 1]       += x.sn[k][m] * x.cs[l][m];
        x.A[(2 * l + 1) * x.neq + 2 * k]     += x.cs[k][m] * x.sn[l][m];
        x.A[(2 * l + 1) * x.neq + 2 * k + 1] += x.sn[k][m] * x.sn[l][m];
        }
      }
/*-----------------------------------------------------------------------------
    permanent tide Z0, imaginary parts needs a special treatment*/
    if(x.spectrum.waves[k].omega == 0.0)
      x.A[(2 * k + 1) * x.neq + 2 * k + 1] = 1.0;
    }
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_init01(int nframe, double *time, spectrum_t s, harmonic_t *x, int nodal,const astro_angles_t &astro_angles)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  double averaging_time=0;

  harmonic_start(x,s);
  harmonic_start(x,nframe);

  harmonic_init(x, time, nodal, astro_angles, averaging_time);

  x->pivot = new int[x->neq];
  status=poc_getrf(x->neq, x->A, x->pivot);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int harmonic_init01(int nframe, date_t start, double *time, spectrum_t s, harmonic_t *x, int nodal,const astro_angles_t &astro_angles)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=harmonic_init01(nframe,time,s,x,nodal,astro_angles);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void harmonic_free(harmonic_t harmonic)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  zaparr(harmonic.A);
  harmonic_free_rhs(harmonic);
  if(harmonic.cs!=NULL)
    zaparr(harmonic.cs);
  if(harmonic.sn!=NULL)
    zaparr(harmonic.sn);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void harmonic_free_rhs(harmonic_t harmonic)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,n;
  for(i = 0; i < harmonic.nrhs; i++) {
    for(n = 0; n < harmonic.nndes[i]; n++) {
      zaparr(harmonic.rhs[i][n]);
      }
    zaparr(harmonic.rhs[i]);
    }
}


#if 0
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void harmonic_save(double t, double **b[3], int *keep, int idx)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  save the harmonic constants analysed at mesh nodes
----------------------------------------------------------------------*/
{
  int i, k, m, n, missing = 0,status;
  float a, G, a1, G1, a2, G2, x, y;
  tidal_wave wave;
  double zi, zr;
  char filename[1024];
  FILE *out;
  char *sdate,*sstart;

}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

hconstant_t **harmonic_analysis_core(harmonic_t harmonic, double duration, int option)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/**  perform the harmonic analysis at grid nodes
\param harmonic the harmonic_t structure to solve, previously built with harmonic_storage().
\param duration acquisition duration in seconds
\param option if \c option & 2, print debug informations.
  If \c option & 1, use the admittance method for poorly separated waves.
\returns a hconstant_t[harmonic.nrhs][harmonic.nndes[]] array. See #harmonic_analysis_constants_USE_polar about whether it is initialised in polar or in complex

See the documentation of harmonic_storage(double,harmonic_t,int,float**,astro_angles_t) for the maths of this.
\todo Document here the maths done here about wave separation.
*/
/*----------------------------------------------------------------------------*/
{
  int i, j, k, l, m, n;//<2 wave indexes, A value index,, variable index, node index
  int neq = harmonic.neq;
  int status, step;
  int *pivot, *keep, *deduce;
  double *A, ***b, tau, factor;
  double coef[3];
  double zi, zr;
  int destructive=1;
  hconstant_t **constants;
  ///\note The separation criteria is twice the analysis period, which is OK for model outputs.
  const double days=duration/d2s,separation=days*2.;
  
//  extern basewave_s;
  if(option & 2){/// <h1>if \c option & 2</h1>
    ///print debug informations
    STDERR_BASE_LINE("A[%d]:",neq);
    for(k=0;k<neq;k++){
      fprintf(stderr," %g",harmonic.A[k]);
      }
    fprintf(stderr,"\n");
    STDERR_BASE_LINE("rhs[0<%d][%d][%d]:",harmonic.nrhs,harmonic.nndes,neq);
    if(harmonic.nndes[0]<neq){
      for(n=0;n<harmonic.nndes[0];n++){
        for(k=0;k<neq;k++){
          fprintf(stderr," %g",harmonic.rhs[0][n][k]);
          }
        }
      }
    else{
      fprintf(stderr,"...");
      }
    fprintf(stderr,"\n");
    STDERR_BASE_LINE("%g %d\n",duration,option);
    }
  
  const spectrum_t &AnalysisList=harmonic.spectrum;
  admittance_t admittance;
/*----------------------------------------------------------------------------
  allow/disallow admittances */
  if(option & 1) {
    delete[] admittance_mts_check(&admittance,AnalysisList);
    }
  factor = 24. / 360.;

/*---------------------------------------------------------------------*//**<h1>
  sort out what can and can not be analysed </h1>*/
  STDOUT_BASE_LINE("#harmonic analysis (%d,period : %gd,separation limit : %gd) ...\n",option & 1,days,separation);

  keep = new int[AnalysisList.n];
  deduce = new int[AnalysisList.n];
  for(i = 0; i < AnalysisList.n; i++)
    keep[i] = 1;
  for(i = 0; i < AnalysisList.n; i++)
    deduce[i] = 0;

/*------------------------------------------------------------------------------
  check for unresolved frequencies */
  for(i = 0; i < AnalysisList.n; i++) {
    deduce[i] = 0;
    if(keep[i] == 0)
      continue;
    if(AnalysisList.waves[i].omega == 0.0)
      continue; /*Z0 special case */
    tau = (1. / AnalysisList.waves[i].omega) * 360. / 24.;
    if(fabs(tau) > days) {
      keep[i] = 0;
      printf("wave: %10s, period: %9.3f days \n", AnalysisList.waves[i].name,tau);
      }
    }

/*------------------------------------------------------------------------------
  check for redondant frequencies and non-separable constituents*/
  for(i = 0; i < AnalysisList.n; i++) {
    if(keep[i] == 0)
      continue;
    for(j = i + 1; j < AnalysisList.n; j++) {
      tau = AnalysisList.waves[j].omega - AnalysisList.waves[i].omega;
      if(tau == 0.) {
        keep[j] = 0;
        printf("wave: %10s %10s, separation IMPOSSIBLE!!!\n",
          AnalysisList.waves[i].name, AnalysisList.waves[j].name);
        }
      if(option & 1) {
        tau = deltaOmega2separation(tau);
        if(tau > separation) {
          printf("wave: %10s %10s, separation: %9.3f days.\n",
            AnalysisList.waves[i].name, AnalysisList.waves[j].name, tau);
          if( AnalysisList.waves[i].Ap > AnalysisList.waves[j].Ap ) {
            keep[j] = 0;
            }
          else {
            keep[i] = 0;
            }
          }
        }
      }
    }

/*------------------------------------------------------------------------------
  check for non-separable constituents but admittance solvable*/
  for(i = 0; i < AnalysisList.n; i++) {
    if(keep[i] == 0) {
      tau = deltaOmega2separation(AnalysisList.waves[i].omega);
      if((AnalysisList.waves[i].Ap != 0.0)
         && (tau < separation)) {
/*------------------------------------------------------------------------------
        admittance can do ? */
        if(admittance_mts_verify(AnalysisList.waves[i],admittance.method) != 0)
          continue;
/*------------------------------------------------------------------------------
        yes it can */
        keep[i] = 1;
        deduce[i] = 1;
        printf("wave: %10s, deduced.\n", AnalysisList.waves[i].name);
        }
      }
    }

  pivot = new int[neq];

  if(destructive==1) {
    A = harmonic.A;
    b = harmonic.rhs;
    }
  else {
    pivot = new int[neq];
    A = copy(harmonic.A, neq*neq);
    b = new double**[harmonic.nrhs];
    for(m = 0; m < harmonic.nrhs; m++) {
      b[m] = new double *[harmonic.nndes[m]];
      for(n = 0; n < harmonic.nndes[m]; n++) {
        b[m][n] = copy(harmonic.rhs[m][n],neq);
        }
      }
    }
/*------------------------------------------------------------------------------
  adjust right-hand side members for unresolved/admitted waves*/

  for(m = 0; m < harmonic.nrhs; m++) {
    for(i = 0; i < AnalysisList.n; i++) {
/*------------------------------------------------------------------------------
      solution forced to zero*/
      /* or */
/*------------------------------------------------------------------------------
      admittance rhs is zero*/
      if(keep[i] == 0 or deduce[i] == 1) {
        for(n = 0; n < harmonic.nndes[m]; n++) {
          b[m][n][2 * i]     = 0.;
          b[m][n][2 * i + 1] = 0.;
          }
        }
      }
    }

/*------------------------------------------------------------------------------
  adjust matrix for unresolved/admitted waves*/

  for(i = 0; i < AnalysisList.n; i++) {
    if(keep[i] == 0) {
      /**If an unresolved wave can not be analysed with the admittance method, its column in the harmonix matrix harmonic_t::A is replaced with 0s but for its row, replaced with a 1, which takes this wave out of the equation. */
      for(j = 0; j < AnalysisList.n; j++)
        A[2 * j * neq + 2 * i] = 0.;
      for(j = 0; j < AnalysisList.n; j++)
        A[2 * j * neq + 2 * i + 1] = 0.;
      for(j = 0; j < AnalysisList.n; j++)
        A[(2 * j + 1) * neq + 2 * i] = 0.;
      for(j = 0; j < AnalysisList.n; j++)
        A[(2 * j + 1) * neq + 2 * i + 1] = 0.;
      A[2 * i * neq + 2 * i] = 1.0;
      A[(2 * i + 1) * neq + 2 * i + 1] = 1.0;
      }
    else
      printf("wave: %10s analysed \n", AnalysisList.waves[i].name);
/*------------------------------------------------------------------------------

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
      printf("wave: %10s, analysed in admittance mode \n", AnalysisList.waves[i].name);
      ///If an unresolved wave can be analysed with the admittance method, its column in the harmonix matrix harmonic_t::A is replaced with weights calculated from the astronomical potential of the waves, i.e. :
      ///-0s for non-astronomical waves
      for(j = 0; j < neq; j++)
        A[j * neq + 2 * i] = 0.;
      for(j = 0; j < neq; j++)
        A[j * neq + 2 * i + 1] = 0.;
      ///- weights calculated from the inverse of the astronomical potential wave_t::Ap of this wave
      A[2 * i * neq + 2 * i] = 1. / AnalysisList.waves[i].Ap;
      A[(2 * i + 1) * neq + 2 * i + 1] = 1. / AnalysisList.waves[i].Ap;
      ///- weights calculated by admittance_sweightP()
      admittance_mts_sweightP(admittance,AnalysisList.waves[i], coef);
      k = AnalysisList.waves[i].nT;
      for(l = 0; l < 3; l++) {
        j = admittance.windex_s[k][l];
        if(j==-1)TRAP_ERR_EXIT(ENOEXEC,"admittance waves index for nT=%d not initialized\n",k);
        /*        printf(" %d %d %s \n",j,windex_s[k][l],AnalysisList.waves[j].name); */
        A[2 * j * neq + 2 * i] = -coef[l];
        A[(2 * j + 1) * neq + 2 * i + 1] = -coef[l];
        }
      }
    }

  status = -1;

  fprintf(stderr,"\r\e[K");//clear line
/*---------------------------------------------------------------------*//**<h1>
  solve the equation </h1>*/
  STDERR_BASE_LINE("#harmonic analysis: factorisation");

  status=poc_getrf(neq, A, pivot);

  if(status != 0) {
    check_error(-1, "singular harmonic matrix", __LINE__, __FILE__);
    }

  fprintf(stderr,"\r\e[K");//clear line
  STDERR_BASE_LINE("#harmonic analysis: solve");

  for(m = 0; m < harmonic.nrhs; m++) {
    for(n = 0; n < harmonic.nndes[m]; n++) {
      for(step = 0; step < 1; step++) {
        int nrhs=1;
        status=poc_getrs(neq, nrhs, A, pivot, b[m][n]);
        }
      }
    }

  fprintf(stderr,"\r\e[K");//clear line
/*---------------------------------------------------------------------*//**<h1>
  fill \c constants </h1>*/
  STDERR_BASE_LINE("#harmonic analysis: constants alloc");

  constants=new hconstant_t *[harmonic.nrhs];
  for(m = 0; m < harmonic.nrhs; m++) {
    constants[m] = new hconstant_t [harmonic.nndes[m]];
    for(n = 0; n < harmonic.nndes[m]; n++) {
      ///depending on #harmonic_analysis_constants_USE_polar
      #if harmonic_analysis_constants_USE_polar
      constants[m][n].init_polar(AnalysisList.n);
      #else
      constants[m][n].init_complex(AnalysisList.n);
      #endif
      for(k = 0; k < AnalysisList.n; k++) {
        zr = b[m][n][2 * k];
        zi = b[m][n][2 * k + 1];
        #if harmonic_analysis_constants_USE_polar
        a = sqrt(zr * zr + zi * zi);
        G = atan2(zi, zr);
        constants[m][n].a[k] = a;
        constants[m][n].G[k] = G*r2d;
        #else
        constants[m][n].z[k]=complex<float>(zr,-zi);
        #endif
        }
      }
    }

  fprintf(stderr,"\r\e[K");//clear line
  STDERR_BASE_LINE("#harmonic analysis: clean up");
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check : MANDATORY !!!

  Notes


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  if(pivot!=0) delete[] pivot;
  if(destructive==0) {
    delete[]A;
    for(m = 0; m < harmonic.nrhs; m++) {
      for(n = 0; n < harmonic.nndes[m]; n++) {
        delete[] b[m][n];
        }
      delete[] b[m];
      }
    delete[]b;
    }

  delete[] keep;

  fprintf(stderr,"\r\e[K");//clear line
  STDERR_BASE_LINE("#harmonic analysis: done\n");

  return(constants);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 extern "C" void harmonic_analysis_(harmonic_t *harmonic, double *duration, int *option, float *rs, float *is)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Fortran wrapper for harmonic_analysis_core(harmonic_t harmonic, double duration, int option)
/**
\param rs real parts
\param is imaginary parts
*/
/*----------------------------------------------------------------------------*/
{
  hconstant_t **constants,*constantsm;
  int n,k;//<node and waves indexes
  int c;//<constant index
  int nw=harmonic->spectrum.n,nn=harmonic->nndes[0];
  
  constants=harmonic_analysis_core(*harmonic,*duration,*option);
  constantsm=constants[0];//NOTE:nrhs=1
  
  for(n=0;n<nn;n++){
    for(k=0;k<nw;k++){
      #if 0
      ///\note give (wave,node) as Fortran indexes
      c=n*nw+k;
      #else
      ///\note give (node,wave) as Fortran indexes
      c=k*nn+n;
      #endif
      #if harmonic_analysis_constants_USE_polar
      TRAP_ERR_EXIT(ENOEXEC,"not coded yet\n");
      #else
      rs[c]=real(constantsm[n].z[k]);
      is[c]=imag(constantsm[n].z[k]);
      #endif
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void harmonic_correction(double time, spectrum_t AnalysisList, hconstant_t **constants,float **buffer,const size_t *nndes, int nrhs, int nodal_corrections, const astro_angles_t &astro_angles)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  perform the harmonic correction at grid nodes
----------------------------------------------------------------------*/
{
  int k, m, n, status;
  float *cs, *sn;
  float tide;
  date_t start;

  if(AnalysisList.n == 0)
    return;

  cs = new float[AnalysisList.n];
  sn = new float[AnalysisList.n];

/* *-----------------------------------------------------------------------------
  to be secured later for nodal corrections precision*/
//   start=cnesdate(time, 'd');
//   init_argument(start);
//   time=0.0;

/* *-----------------------------------------------------------------------------
  compute tidal coefficients*/
  status = harmonic_coefficients(time, AnalysisList, cs, sn, nodal_corrections, astro_angles);

  for(n=0;n<nrhs;n++) {
    for(m=0; m<nndes[n]; m++) {
      tide=0;
      for(k=0; k < AnalysisList.n; k++) {
        tide+=constants[n][m].z[k].real()*cs[k]-constants[n][m].z[k].imag()*sn[k];
        }
      buffer[n][m]-=tide;
      }
    }

  delete[] cs;
  delete[] sn;
}


/*----------------------------------------------------------------------------*/
/// Updates the harmonic matrix harmonic_t::A and right-hand side vector harmonic_t::rhs
/** of a harmonic_t structure
at the given time and with the given values.

\param t time in seconds as \b double. The reference must be the same for each related calls of this function : see harmonic_coefficients().
\param harmonic harmonic_t structure to update
\param nodal_corrections whether nodal corrections are carried out.
\param buffer <tt>float[harmonic.nrhs][harmonic.nndes]</tt> buffer of values.

\bug This should be used way more often. See bug note of harmonic_coefficients()

We have :
\f{eqnarray*}{
\left[ \nu_n e^{j\left(\omega_n t + \phi_n\right)} \right] \left[ x_n \right] &=& \left[ h_t \right] \\
\left[ \nu_n e^{j\left(\omega_n t + \phi_n\right)} \right]^* \left[ \nu_n e^{j\left(\omega_n t + \phi_n\right)} \right] \left[ x_n \right] &=& \left[ \nu_n e^{j\left(\omega_n t + \phi_n\right)} \right]^* \left[ h_t \right]
\f}
with :
- \f$ \nu_n \f$ a complex number giving the nodal correction (in amplitude and in phase)
- \f$ w_n \f$ the pulsation of the wave
- \f$ t \f$ the time since the reference
- \f$ \phi_n \f$ the astronomic angle of the wave
- the modulus and argument of \f$ x_n \f$ are respectively the amplitudes and phases of the analysed waves
- \f$ A \equiv \left[ \nu_n e^{j\left(\omega_n t + \phi_n\right)} \right]^* \left[ \nu_n e^{j\left(\omega_n t + \phi_n\right)} \right] \f$ the harmonic matrix harmonic_t::A
- \f$ b \equiv \left[ \nu_n e^{j\left(\omega_n t + \phi_n\right)} \right]^* \left[ h_t \right] \f$ the right-hand side vector harmonic_t::rhs
- \f$ x \equiv \left[ x_n \right] \f$ the harmonic coefficients

The solution of the equation is calculated by harmonic_analysis_core().
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> void harmonic_storage_template(double t, harmonic_t harmonic, int nodal_corrections, T **buffer, const astro_angles_t &astro_angles)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l, m, n, neq = harmonic.neq, status;
  float *cs, *sn;

  struct timeval before;
  gettimeofday(&before);
  spectrum_t AnalysisList=harmonic.spectrum;

  if(AnalysisList.n == 0)
    return;

  cs = new float[AnalysisList.n];
  sn = new float[AnalysisList.n];

  ///It calculates \f$ \nu_n e^{j\left(\omega_n t + \phi_n\right)} \f$ with harmonic_coefficients().
  status = harmonic_coefficients(t, AnalysisList, cs, sn, nodal_corrections, astro_angles);

  for(m = 0; m < harmonic.nrhs; m++) {
    T *bufferm=buffer[m];
    double **rhsm=harmonic.rhs[m];
    #pragma omp parallel for private(k,n)
    for(n = 0; n < harmonic.nndes[m]; n++) {
      T h = bufferm[n];
      double *rhsmn=rhsm[n];//pointer : used for speed
      for(k = 0; k < AnalysisList.n; k++) {
        rhsmn[2 * k]     += h * cs[k];
        rhsmn[2 * k + 1] += h * sn[k];
        }
      }
    }
/*-----------------------------------------------------------------------------
 line 2k is derivation with respect to real part
 line 2k+1 is derivation with respect to imaginary part
 column 2l is coefficient with respect to real part
 column 2l+1 is derivation with respect to imaginary part
------------------------------------------------------------------------------*/
  for(k = 0; k < AnalysisList.n; k++) {
    for(l = 0; l < AnalysisList.n; l++) {
      harmonic.A[2 * l * neq + 2 * k]           += cs[k] * cs[l];
      harmonic.A[2 * l * neq + 2 * k + 1]       += sn[k] * cs[l];
      harmonic.A[(2 * l + 1) * neq + 2 * k]     += cs[k] * sn[l];
      harmonic.A[(2 * l + 1) * neq + 2 * k + 1] += sn[k] * sn[l];
      }
/*-----------------------------------------------------------------------------
    permanent tide Z0, imaginary parts needs a special treatment*/
    if(AnalysisList.waves[k].omega == 0.0) {
      harmonic.A[(2 * k + 1) * neq + 2 * k + 1] = 1.0;
      }
    }
  delete[] cs;
  delete[] sn;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void harmonic_storage(double t, harmonic_t harmonic, int nodal_corrections, double **buffer, const astro_angles_t &astro_angles)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{harmonic_storage_template(t,harmonic,nodal_corrections,buffer,astro_angles);}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void harmonic_storage(double t, harmonic_t harmonic, int nodal_corrections, float **buffer, const astro_angles_t &astro_angles)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{harmonic_storage_template(t,harmonic,nodal_corrections,buffer,astro_angles);}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

extern "C" void harmonic_storage4_t_(harmonic_t *harmonic, double *t, float *buffer, int *nodal_corrections, astro_angles_t *astro_angles)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Fortran wrapper for harmonic_storage(double t, harmonic_t harmonic, int nodal_corrections, float **buffer, const astro_angles_t &astro_angles)
/*----------------------------------------------------------------------------*/
{
  harmonic_storage(*t,*harmonic,*nodal_corrections,&buffer,*astro_angles);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

extern "C" void harmonic_storage8_t_(harmonic_t *harmonic, double *t, double *buffer, int *nodal_corrections, astro_angles_t *astro_angles)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Fortran wrapper for harmonic_storage(double t, harmonic_t harmonic, int nodal_corrections, float **buffer, const astro_angles_t &astro_angles)
/*----------------------------------------------------------------------------*/
{
  harmonic_storage(*t,*harmonic,*nodal_corrections,&buffer,*astro_angles);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

extern "C" void harmonic_get_a_and_rhs_(harmonic_t *harmonic, double *A, double *rhs)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Fortran wrapper for reading harmonic_t::A and harmonic_t::rhs
/**
\param *A array[harmonic_t::neq*harmonic_t::neq] of double
\param *rhs array[harmonic_t::nndes*harmonic_t::neq] of double

\sa harmonic_set_a_and_rhs_()
*/
/*----------------------------------------------------------------------------*/
{
  int n,k,i;//<node or wave index, wave index, array index
  
  int nn=harmonic->nndes[0],neq=harmonic->neq,neq2=neq*neq;
  
  double *hA=harmonic->A,**hrhs=harmonic->rhs[0],*rhsn;//NOTE:nrhs=1
  
  for(i=0;i<neq2;i++){
    A[i]=hA[i];
    }
  for(n=0,i=0;n<nn;n++,i+=neq){
    rhsn=hrhs[n];
    for(k=0;k<neq;k++){
      rhs[i+k]=rhsn[k];
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 extern "C" void harmonic_set_a_and_rhs_(harmonic_t *harmonic, double *A, double *rhs)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Fortran wrapper for writing to harmonic_t::A and harmonic_t::rhs
/**
\param *A array[harmonic_t::neq*harmonic_t::neq] of double
\param *rhs array[harmonic_t::nndes*harmonic_t::neq] of double

\sa harmonic_get_a_and_rhs_()
*/
/*----------------------------------------------------------------------------*/
{
  int n,k,i;//<node or wave index, wave index, array index
  
  int nn=harmonic->nndes[0],neq=harmonic->neq,neq2=neq*neq;
  
  double *hA=harmonic->A,**hrhs=harmonic->rhs[0],*rhsn;//NOTE:nrhs=1
  
  for(i=0;i<neq2;i++){
    hA[i]=A[i];
    }
  for(n=0,i=0;n<nn;n++,i+=neq){
    rhsn=hrhs[n];
    for(k=0;k<neq;k++){
      rhsn[k]=rhs[i+k];
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void harmonic_correction(const harmonic_t harmonic, int frame, hconstant_t **constants, float **buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// perform the harmonic correction of the given values at the given time
/**
\param harmonic for the sines and cosines
\param frame index of the frame
\param **constants <tt>hconstant_t[harmonic.nrhs][harmonic.nndes]</tt> array from harmonic_analysis_core()
\param **buffer <tt>float[harmonic.nrhs][harmonic.nndes]</tt> buffer of values.
*/
/*----------------------------------------------------------------------------*/
{
  int k, m, n;
  float tide;//tidal prediction
  float *cs, *sn;//sine and cosine buffer : used for speed
  hconstant_t *constn;//pointer : used for speed
  bool *maskn=0;//pointer : used for speed
  
  const int nw=harmonic.spectrum.n;
  
  if(nw == 0)
    return;
  cs = new float[nw];
  sn = new float[nw];
  for(k = 0; k < nw; k++) {
    cs[k]=harmonic.cs[k][frame];
    sn[k]=harmonic.sn[k][frame];
    }
  
  if(constants[0][0].z==0)
    TRAP_ERR_EXIT(ENOEXEC,"programming error : %s() must be called with a constants with z allocated\n",__func__);
  
  for(n=0;n<harmonic.nrhs;n++) {
    constn=constants[n];
    float *buffern=buffer[n];//pointer : used for speed
    if(harmonic.mask)
      maskn=harmonic.mask[n];
    
    #pragma omp parallel for private(m,tide,k)
    for(m=0; m<harmonic.nndes[n]; m++) {
      if(harmonic.mask and not maskn[m])
        continue;
      
      complex<float> *constnmz=constn[m].z;//pointer : used for speed
      tide=0.;
      for(k=0; k < nw; k++) {
        tide+=constnmz[k].real()*cs[k]-constnmz[k].imag()*sn[k];
        }
      buffern[m]-=tide;
      }
    }
  
  delete[] cs;
  delete[] sn;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> void harmonic_storage_template(const harmonic_t &harmonic, int frame, T **buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Updates the harmonic matrix harmonic_t::A and right-hand side vector harmonic_t::rhs
/** of a harmonic_t structure
at the given time and with the given values.

\param harmonic harmonic_t structure to update. See harmonic_init() for its initialisation.
\param frame index of the frame
\param buffer <tt>float[harmonic.nrhs][harmonic.nndes]</tt> buffer of values.

\bug This should be used way more often.

See the documentation of harmonic_storage(double,harmonic_t,int,float**,astro_angles_t) for the maths of this.

The solution of the equation is calculated by harmonic_analysis_core().
*/
/*----------------------------------------------------------------------------*/
{
  int k, m, n;
  double *cs, *sn;//sine and cosine buffer : used for speed
  double **rhsm;//pointer : used for speed
  T *bufferm;//pointer : used for speed
  bool *maskm=0;//pointer : used for speed
  
  if(harmonic.spectrum.n == 0)
    return;
  cs = new double[harmonic.spectrum.n];
  sn = new double[harmonic.spectrum.n];
  for(k = 0; k < harmonic.spectrum.n; k++) {
    cs[k]=harmonic.cs[k][frame];
    sn[k]=harmonic.sn[k][frame];
    }
  /// \bug 2011-09-16 Damien ALLAIN : this does not check for proper initialisation
  
  for(m = 0; m < harmonic.nrhs; m++) {
    rhsm=harmonic.rhs[m];
    bufferm=buffer[m];
    if(harmonic.mask)
      maskm=harmonic.mask[m];
    
    #pragma omp parallel for private(k,n)
    for(n = 0; n < harmonic.nndes[m]; n++) {
      T h;
      if(harmonic.mask and not maskm[n])
        continue;
      
      h = bufferm[n];
      double *rhsmn=rhsm[n];//pointer : used for speed
      for(k = 0; k < harmonic.spectrum.n; k++) {
        rhsmn[2 * k]     += h * cs[k];
        rhsmn[2 * k + 1] += h * sn[k];
        }
      }
    }
  
  delete[] cs;
  delete[] sn;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void harmonic_storage(const harmonic_t &harmonic, int frame, double **buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{harmonic_storage_template(harmonic,frame,buffer);}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void harmonic_storage(const harmonic_t &harmonic, int frame, float **buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{harmonic_storage_template(harmonic,frame,buffer);}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 extern "C" void harmonic_storage4_i_(harmonic_t *harmonic, int *frame, float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Fortran wrapper for harmonic_storage(const harmonic_t &harmonic, int frame, float **buffer)
/*----------------------------------------------------------------------------*/
{
  harmonic_storage(*harmonic,*frame-1,&buffer);//NOTE:nrhs=1
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 extern "C" void harmonic_storage8_i_(harmonic_t *harmonic, int *frame, double *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Fortran wrapper for harmonic_storage(const harmonic_t &harmonic, int frame, double **buffer)
/*----------------------------------------------------------------------------*/
{
  harmonic_storage(*harmonic,*frame-1,&buffer);//NOTE:nrhs=1
}
