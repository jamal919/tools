
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
#include <string.h>

#include "tools-structures.h"
#include "fe.def"
#include "fe.h"
#include "map.h"
#include "geo.h"
#include "tides.h"
#include "tides.def"

#include "admittance.h"

#ifndef VERBOSE
#define VERBOSE 1
#endif

static double response_s[nspecies][3][3];
static tidal_wave basewave_s[nspecies][3];
int windex_s[nspecies][3]={{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

static double response_l[nspecies][2][2];
static tidal_wave basewave_l[nspecies][2];
int windex_l[nspecies][2]={{-1,-1},{-1,-1},{-1,-1}};

static int method[nspecies];


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void admittance_test(spectrum_t WaveList,spectrum_t AdmittanceList, int nby)
//unused
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j, k, l;
  double delta[3] = {0, 0, 0}, cs[3] = {0, 0, 0}, sn[3] = {0, 0, 0};
  double hx = 0, hy = 0;
  double c1[3] = {0, 0, 0}, c2[3] = {0, 0, 0}, matrix[9], ri[3][3];
  tidal_wave wave[3], w;
  fcomplex elevation[3];
  float amp, pha;

  k = nby - 1;
  for(j = 0; j < WaveList.n; j++) {
    printf("given %s %f %f \n", WaveList.waves[j].name, amp, pha * r2d);
    if(strcmp(WaveList.waves[j].name, "M2") == 0) {
      elevation[0]= polar<float>(amp,pha);
      }
    if(strcmp(WaveList.waves[j].name, "N2") == 0) {
      elevation[1]= polar<float>(amp,pha);
      }
    if(strcmp(WaveList.waves[j].name, "K2") == 0) {
      elevation[2]= polar<float>(amp,pha);
      }
    }

  size_t size = sizeof(tidal_wave);

  ///\date 2011-10-10 Damien Allain : corrected memmove mistake on argument ordering
  memmove(&wave[0], &wM2, size);
  memmove(&wave[1], &wN2, size);
  memmove(&wave[2], &wK2, size);

  for(k = 0; k < 3; k++)
    wave[k].init();
  
  for(k = 0; k < 3; k++)
    delta[k] = (wave[k].omega - wave[0].omega) * dph2cpd;
  

  for(k = 0; k < 3; k++) {
    cs[k] = cos(4. * M_PI * delta[k]);
    sn[k] = sin(4. * M_PI * delta[k]);
    }

  for(k = 0; k < 3; k++) {
    matrix[k] = cs[k];
    matrix[k + 3] = sn[k];
    matrix[k + 6] = 1.0;
    }
  
  int neq = 3;

#ifdef LINPACKF_
  dgefa_(matrix, &neq, &neq, pivot, &status);
#elif LINPACKF
  dgefa(matrix, &neq, &neq, pivot, &status);
#elif LINPACKC
  int rcond;
  dgeco(matrix, neq, neq, pivot, &rcond);
/*   dgefa(matrix, neq, neq, pivot,&status); */
#endif

  for(k = 0; k < 3; k++){
    c1[k] = real(elevation[k]) / wave[k].Ap;
  }
  
  int pivot[3], status = 0, job = 0;
#ifdef LINPACKF_
  for(k = 0; k < 3; k++)
    dgesl_(matrix, &neq, &neq, pivot, c1, &job);
#elif LINPACKF
  for(k = 0; k < 3; k++)
    dgesl(matrix, &neq, &neq, pivot, c1, &job);
#elif LINPACKC
  for(k = 0; k < 3; k++)
    dgesl(matrix, neq, neq, pivot, c1, job);
#endif

  for(k = 0; k < 3; k++)
    c2[k] = imag(elevation[k]) / wave[k].Ap;
#ifdef LINPACKF_
  for(k = 0; k < 3; k++)
    dgesl_(matrix, &neq, &neq, pivot, c2, &job);
#elif LINPACKF
  for(k = 0; k < 3; k++)
    dgesl(matrix, &neq, &neq, pivot, c2, &job);
#elif LINPACKC
  for(k = 0; k < 3; k++)
    dgesl(matrix, neq, neq, pivot, c2, job);
#endif
  
  double gg;

  for(j = 0; j < AdmittanceList.n; j++) {
    bcopy((char *) &w, (char *) &AdmittanceList.waves[j], sizeof(tidal_wave));
    gg = (w.omega - wave[0].omega) * dph2cpd * 4. * M_PI;
    hx = w.Ap * (c1[0] * cos(gg) + c1[1] * sin(gg) + c1[2]);
    hy = w.Ap * (c2[0] * cos(gg) + c2[1] * sin(gg) + c2[2]);
    amp = sqrt(hx * hx + hy * hy);
    pha = atan2(hy, hx);
    printf("deduced %s %f %f \n", AdmittanceList.waves[j].name, amp, pha * r2d);
    }

  for(k = 0; k < 3; k++)
    for(l = 0; l < 3; l++)
      ri[k][l] = 0.0;

  for(k = 0; k < 3; k++)
    ri[k][k] = 1.0;

  job = 0;
#ifdef LINPACKF_
  for(k = 0; k < 3; k++)
    dgesl_(matrix, &neq, &neq, pivot, ri[k], &job);
#elif LINPACKF
  for(k = 0; k < 3; k++)
    dgesl(matrix, &neq, &neq, pivot, ri[k], &job);
#elif LINPACKC
  for(k = 0; k < 3; k++)
    dgesl(matrix, neq, neq, pivot, ri[k], job);
#endif

  for(j = 0; j < AdmittanceList.n; j++) {
    bcopy((char *) &w, (char *) &AdmittanceList.waves[j], sizeof(tidal_wave));
    for(k = 0; k < 3; k++){
      c1[k] = 0.;
      c2[k] = 0.;
      }
    for(k = 0; k < 3; k++){
      for(l = 0; l < 3; l++){
        c1[l] += ri[k][l] * real(elevation[k]) / wave[k].Ap;
        c2[l] += ri[k][l] * imag(elevation[k]) / wave[k].Ap;
        }
      }
    gg = (w.omega - wave[0].omega) * dph2cpd * 4. * M_PI;
    hx = w.Ap * (c1[0] * cos(gg) + c1[1] * sin(gg) + c1[2]);
    hy = w.Ap * (c2[0] * cos(gg) + c2[1] * sin(gg) + c2[2]);
    amp = sqrt(hx * hx + hy * hy);
    pha = atan2(hy, hx);
    printf("deduced %s %f %f \n", AdmittanceList.waves[j].name, amp, pha * r2d);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int admittance_init(int option)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
   size_t size = sizeof(tidal_wave);

/* *-------------------------------------------------------------------------
   spline interpolation*/
   bcopy((char *) &(wM2),(char *) &(basewave_s[2][0]),  size);
   bcopy((char *) &(wN2),(char *) &(basewave_s[2][1]),  size);
   switch(option){
     case 0:
       bcopy((char *) &(wK2),(char *) &(basewave_s[2][2]),  size);
       break;
     case 1:
       bcopy((char *) &(wK2),(char *) &(basewave_s[2][2]),  size);
       break;
     }
   admittance_scoefficient(basewave_s[2], response_s[2]);

   bcopy((char *) &(wK1),(char *) &(basewave_s[1][0]),  size);
   bcopy((char *) &(wO1),(char *) &(basewave_s[1][1]),  size);
   bcopy((char *) &(wQ1),(char *) &(basewave_s[1][2]),  size);
   admittance_scoefficient(basewave_s[1], response_s[1]);

   bcopy((char *) &(wMf), (char *) &(basewave_s[0][0]),  size);
   bcopy((char *) &(wMm), (char *) &(basewave_s[0][1]),  size);
   bcopy((char *) &(wMtm),(char *) &(basewave_s[0][2]),  size);
   admittance_scoefficient(basewave_s[0], response_s[0]);

/* *-------------------------------------------------------------------------
   linear interpolation*/
   bcopy((char *) &(wM2),(char *) &(basewave_l[2][0]),  size);
   bcopy((char *) &(wK2),(char *) &(basewave_l[2][1]),  size);
   admittance_lcoefficient(basewave_l[2], response_l[2]);

   bcopy((char *) &(wK1),(char *) &(basewave_l[1][0]),  size);
   bcopy((char *) &(wO1),(char *) &(basewave_l[1][1]),  size);
   admittance_lcoefficient(basewave_l[1], response_l[1]);

   bcopy((char *) &(wMf),(char *) &(basewave_l[0][0]),  size);
   bcopy((char *) &(wMm),(char *) &(basewave_l[0][1]),  size);
   admittance_lcoefficient(basewave_l[0], response_l[0]);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int admittance_init(const spectrum_t &basisList, const tidal_wave &w)

 /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  size_t size = sizeof(tidal_wave);

  switch (basisList.n){
  case 3:{

/* *-------------------------------------------------------------------------
  spline interpolation*/
  switch(w.nT){
    case 2:{
      bcopy((char *) &(basisList.waves[0]),(char *) &(basewave_s[2][0]),  size);
      bcopy((char *) &(basisList.waves[1]),(char *) &(basewave_s[2][1]),  size);
      bcopy((char *) &(basisList.waves[2]),(char *) &(basewave_s[2][2]),  size);
      admittance_scoefficient(basewave_s[2], response_s[2]);
      return(0);
      }
    case 1:{
      bcopy((char *) &(basisList.waves[0]),(char *) &(basewave_s[1][0]),  size);
      bcopy((char *) &(basisList.waves[1]),(char *) &(basewave_s[1][1]),  size);
      bcopy((char *) &(basisList.waves[2]),(char *) &(basewave_s[1][2]),  size);
      admittance_scoefficient(basewave_s[1], response_s[1]);
      return(0);
      }
    case 0:{
      bcopy((char *) &(basisList.waves[0]),(char *) &(basewave_s[0][0]),  size);
      bcopy((char *) &(basisList.waves[1]),(char *) &(basewave_s[0][1]),  size);
      bcopy((char *) &(basisList.waves[2]),(char *) &(basewave_s[0][2]),  size);
      admittance_scoefficient(basewave_s[0], response_s[0]);
      return(0);
      }
    default: return(1);
    }

  }
  case 2:{
/* *-------------------------------------------------------------------------
  linear interpolation*/
  switch(w.nT){
    case 2:{
      bcopy((char *) &(basisList.waves[0]),(char *) &(basewave_l[2][0]),  size);
      bcopy((char *) &(basisList.waves[1]),(char *) &(basewave_l[2][1]),  size);
      admittance_lcoefficient(basewave_l[2], response_l[2]);
      return(0);
    }

    case 1:{
      bcopy((char *) &(basisList.waves[0]),(char *) &(basewave_l[1][0]),  size);
      bcopy((char *) &(basisList.waves[1]),(char *) &(basewave_l[1][1]),  size);
      admittance_lcoefficient(basewave_l[1], response_l[1]);
      return(0);
    }

    case 0:{
      bcopy((char *) &(basisList.waves[0]),(char *) &(basewave_l[0][0]),  size);
      bcopy((char *) &(basisList.waves[1]),(char *) &(basewave_l[0][1]),  size);
      admittance_lcoefficient(basewave_l[0], response_l[0]);
      return(0);
    }
    default: return(1);
  }

  }// end case 2 (linear)

  default: return(1);
  }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void admittance_terminate()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  for(int k = 0; k< nspecies; k++){ // species
    // spline admittance functions case
    for(int l = 0; l < 3; l++){
      for(int i=0; i < 3; i++) response_s[k][l][i] = 0;
      }

    // linear admittance functions case
    for(int l = 0; l < 2; l++){
      for(int i = 0; i < 2; i++) response_l[k][l][i] = 0;
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *admittance_check(const spectrum_t &WaveList, int option)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///check admittance availability, assuming WaveList spectrum is fully solvable
/**
\param WaveList list of all waves
\param option given to admittance_init(int), where it has no effect
\returns an array[WaveList.n] : #ADMITTANCE_UNAVAILABLE, #ADMITTANCE_LINEAR, #ADMITTANCE_SPLINE
*/
/*----------------------------------------------------------------------------*/
{
  int count_s[3] = {0, 0, 0};
  int count_l[3] = {0, 0, 0};

  for(size_t k = 0; k < 3; k++) {
    for(size_t l = 0; l < 3; l++) windex_s[k][l] = -1;
    for(size_t l = 0; l < 2; l++) windex_l[k][l] = -1;
    }

  admittance_init(option);

  for(size_t j = 0; j < WaveList.n; j++) {
    for(size_t k = 0; k < 3; k++) {
      for(size_t l = 0; l < 3; l++) {
        if(strcmp(WaveList.waves[j].name, basewave_s[k][l].name) == 0) {
          windex_s[k][l] = j;
          count_s[k]++;
          }
        }
      for(size_t l = 0; l < 2; l++) {
        if(strcmp(WaveList.waves[j].name, basewave_l[k][l].name) == 0) {
          windex_l[k][l] = j;
          count_l[k]++;
          }
        }
      }
    }

  for(size_t k = 0; k < 3; k++) {
    if(count_s[k] != 3)
      method[k] = ADMITTANCE_LINEAR;
    else
      method[k] = ADMITTANCE_SPLINE;
    }

  for(size_t k = 0; k < 3; k++) {
    if((count_l[k] != 2) && (method[k] == ADMITTANCE_LINEAR))
      method[k] = ADMITTANCE_UNAVAILABLE;
    }

#if VERBOSE > 1
  for(size_t k = 0; k < 3; k++) {
    switch (method[k]) {
      case ADMITTANCE_UNAVAILABLE:
        printf("specie %d : admittance not available \n", k);
        break;
      case ADMITTANCE_LINEAR:
        printf("specie %d : admittance linear \n", k);
        break;
      case ADMITTANCE_SPLINE:
        printf("specie %d : admittance spline \n", k);
        break;
      default:
        TRAP_ERR_EXIT(-1, "admittance method is not known\n");
        break;
      }
    }
#endif

  return (method);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *admittance_check(const spectrum_t &WaveList, const tidal_wave &w)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int count_s[3] = {0, 0, 0};
  int count_l[3] = {0, 0, 0};

  for(size_t k = 0; k < 3; k++) { // LR, init
    for(size_t l = 0; l < 3; l++) windex_s[k][l] = -1;
    for(size_t l = 0; l < 2; l++) windex_l[k][l] = -1;
    }

  admittance_init(WaveList, w);

  for(size_t j = 0; j < WaveList.n; j++) {
    for(size_t k = 0; k < 3; k++) {
      for(size_t l = 0; l < 3; l++) {
        if(strcmp(WaveList.waves[j].name, basewave_s[k][l].name) == 0) {
          windex_s[k][l] = j;
          count_s[k]++;
          }
        }
      for(size_t l = 0; l < 2; l++) {
        if(strcmp(WaveList.waves[j].name, basewave_l[k][l].name) == 0) {
          windex_l[k][l] = j;
          count_l[k]++;
          }
        }
      }
    }

  for(size_t k = 0; k < 3; k++) {
    if(count_s[k] != 3)
      method[k] = ADMITTANCE_LINEAR;
    else
      method[k] = ADMITTANCE_SPLINE;
    }

  for(size_t k = 0; k < 3; k++) {
    if((count_l[k] != 2) && (method[k] == ADMITTANCE_LINEAR))
      method[k] = ADMITTANCE_UNAVAILABLE;
    }

#if VERBOSE > 1
  for(size_t k = 0; k < 3; k++) {
    switch (method[k]) {
      case ADMITTANCE_UNAVAILABLE:
        printf("specie %d : admittance not available \n", k);
        break;
      case ADMITTANCE_LINEAR:
        printf("specie %d : admittance linear \n", k);
        break;
      case ADMITTANCE_SPLINE:
        printf("specie %d : admittance spline \n", k);
        break;
      default:
        TRAP_ERR_EXIT(-1, "admittance method is not known\n");
        break;
     }
   }
#endif

  return (method);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void admittance_sweightP(tidal_wave w, double coef[3])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k = w.nT;
  admittance_sweight(w, basewave_s[k], response_s[k], coef);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int admittance_verify(tidal_wave w)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k = w.nT;
  
  if(w.Ap == 0.0) return (-1);

  if(method[k] == ADMITTANCE_UNAVAILABLE)
    return (-1);
  else
    return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void admittance_lweightP(tidal_wave w, double coef[2])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k = w.nT;
  admittance_lweight(w, basewave_l[k], response_l[k], coef);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void admittance_compute(spectrum_t WaveList, hconstant_t ztide, int *available)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/**
\param ztide list of polar tidal constants, with phases in degrees
*/
/*----------------------------------------------------------------------------*/
{
  int j, k, l;
//   tidal_wave w;
  fcomplex elevation_s[3][3];
  fcomplex elevation_l[3][2];
  float amp, pha;

  for(j = 0; j < WaveList.n; j++) {
     amp = ztide.a[j];
     pha = ztide.G[j] * d2r; // LR, add
//      printf("given %s %f %f \n", WaveList.waves[j].name, amp, pha * r2d);
/*------------------------------------------------------------------------------
    loop on the 3 tidal species*/
    for(k = 0; k < 3; k++) {
/*------------------------------------------------------------------------------
      loop on the 3 base wave*/
      for(l = 0; l < 3; l++) {
        if(strcmp(WaveList.waves[j].name, basewave_s[k][l].name) == 0) {
          elevation_s[k][l] = polar<float>(amp,-pha);
          }
        }
/*------------------------------------------------------------------------------
      loop on the 2 base wave*/
      for(l = 0; l < 2; l++) {
        if(strcmp(WaveList.waves[j].name, basewave_l[k][l].name) == 0) {
          elevation_l[k][l] = polar<float>(amp,-pha);
          }
        }
      }
    }

  for(j =0; j < WaveList.n; j++) {
    if(available[j]) continue;
    //bcopy((char *) &w, (char *) &WaveList.waves[j], sizeof(tidal_wave));
    //bcopy((void *) &w, (void *) &(WaveList.waves[j]), sizeof(tidal_wave)); // LR, test
    //k = w.nT;
    k = WaveList.waves[j].nT; // LR, test
    switch (method[k]) {
      case ADMITTANCE_SPLINE:
//        admittance_scompute(w, basewave_s[k], elevation_s[k], response_s[k],&amp, &pha);
        admittance_scompute(WaveList.waves[j], basewave_s[k], elevation_s[k], response_s[k],&amp, &pha);// LR, test
        ztide.a[j] = (float) amp;
        ztide.G[j] = (float) pha;
        break;

      case ADMITTANCE_LINEAR:
//        admittance_lcompute(w, basewave_l[k], elevation_l[k], response_l[k], &amp, &pha);
        admittance_lcompute(WaveList.waves[j], basewave_l[k], elevation_l[k], response_l[k], &amp, &pha);// LR, test
        ztide.a[j] = (float) amp;
        ztide.G[j] = (float) pha;
        break;

      case ADMITTANCE_UNAVAILABLE:
        ztide.a[j] = 0.0;
        ztide.G[j] = 0.0;
        break;

      default:
        TRAP_ERR_EXIT(-1, "admittance method is not known\n");
        break;
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int admittance_error_unsafe(const spectrum_t& s, const int *deduce, vector<double>& error)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  // init base waves and responses coefficients
  // and assign admittance methods
  admittance_t admittance;
  int *method = admittance_mts_check(&admittance, s);
  delete[]method;

  // compute linear combination of basis waves' error
  double coef[3] = {0.,0.,0.};

  tidal_wave wave;
  for(size_t k = 0; k < s.n; k++) {
    
    if(deduce[k] == 1){
      wave = s.waves[k];
      error[k] = 0;
      
      size_t species = s.waves[k].nT;
      switch(admittance.method[species]){

        case ADMITTANCE_UNAVAILABLE:
          return(1);
          break;

        case ADMITTANCE_SPLINE:
          admittance_mts_sweightP(admittance, wave, coef);
          for(size_t w = 0; w < 3; ++w){
            int j = admittance.windex_s[species][w];
            error[k] += coef[w] * error[j];
            }
          error[k] *= s.waves[k].Ap;
          break;
        
        case ADMITTANCE_LINEAR:
          admittance_mts_lweightP(admittance, wave, coef);
          for(size_t w = 0; w < 2; ++w){
            int j = admittance.windex_l[species][w];
            error[k] += coef[w] / s.waves[j].Ap * error[j];
            }
          break;
          
        default:
          TRAP_ERR_EXIT(-1, "admittance method is not known\n");
          break;
        }
      }
    }

  admittance_mts_terminate(admittance);
 
  return(0);
}

