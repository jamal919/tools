
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

/**----------------------------------------------------------------------------

  MULTI-THREAD SAFE VERSION OF ADMITTANCE

-----------------------------------------------------------------------------*/

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

#include "admittance-mts.h"

#ifndef VERBOSE
#define VERBOSE 2
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int admittance_mts_init(admittance_t *admittance, int option)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *-------------------------------------------------------------------------
   spline interpolation*/
  admittance->basewave_s[2][0]=wM2;
  admittance->basewave_s[2][1]=wN2;
  switch(option){
    case 0:
      admittance->basewave_s[2][2]=wK2;
      break;
    case 1:
      admittance->basewave_s[2][2]=wK2;
      break;
    }
  admittance_scoefficient(admittance->basewave_s[2], admittance->response_s[2]);

  admittance->basewave_s[1][0]=wK1;
  admittance->basewave_s[1][1]=wO1;
  admittance->basewave_s[1][2]=wQ1;
  admittance_scoefficient(admittance->basewave_s[1], admittance->response_s[1]);

  admittance->basewave_s[0][0]=wMf;
  admittance->basewave_s[0][1]=wMm;
  admittance->basewave_s[0][2]=wMtm;
  admittance_scoefficient(admittance->basewave_s[0], admittance->response_s[0]);

/* *-------------------------------------------------------------------------
  linear interpolation*/
#define ALTERNATIVE_LINEAR_PIVOTS 0
  
  admittance->basewave_l[2][0]=wM2;
  admittance->basewave_l[2][1]=wK2;
  admittance_lcoefficient(admittance->basewave_l[2], admittance->response_l[2]);

#if ALTERNATIVE_LINEAR_PIVOTS
#warning see below
STDERR_BASE_LINE("*** USING Q1 INSTEAD OF K1 AS LINEAR DIURNAL PIVOT ***\n");
  admittance->basewave_l[1][0]=wQ1;
#else
  admittance->basewave_l[1][0]=wK1;
#endif
  admittance->basewave_l[1][1]=wO1;
  admittance_lcoefficient(admittance->basewave_l[1], admittance->response_l[1]);

  admittance->basewave_l[0][0]=wMf;
  admittance->basewave_l[0][1]=wMm;
  admittance_lcoefficient(admittance->basewave_l[0], admittance->response_l[0]);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void admittance_mts_init(admittance_t *admittance,const spectrum_t &WaveList,const int *available)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///wrong initialisation
/**
\param *admittance to be initialised
\param WaveList list of all waves
\param *available array[WaveList.n] : 0 when not available, 1 otherwise (default).
*/
/*----------------------------------------------------------------------------*/
{
  int i,j;
  
  expire(20140207,20150207);
  
  for(i=0;i<WaveList.n;i++){
    if(available!=NULL && !available[i])continue;
    const tidal_wave *w=&WaveList.waves[i];
    if(!w->Ap || w->nT>3)continue;
    range_t<double> ApRange;
    
/*--------------------------------------------------------------------*//**<h1>
  append the available waves to the spline list </h1>*/
    for(j=0;j<3 &&
      admittance->basewave_s[w->nT][j].name[0] &&
      strcmp(admittance->basewave_s[w->nT][j].name,w->name)
      ;j++)ApRange << admittance->basewave_s[w->nT][j].Ap;
    if(j<3){
      }
    ///over-writing the weakest when the list is full
    else if(w->Ap>ApRange.min){
      for(j=0;j<3 &&
        admittance->basewave_s[w->nT][j].Ap>ApRange.min;j++);
      }
    else continue;
    admittance->basewave_s[w->nT][j]=*w;
    admittance->windex_s[w->nT][j]=i;
    }
  
/*--------------------------------------------------------------------*//**<h1>
  calculate coefficients </h1>*/
  for(i=0;i<3;i++){
    for(j=0;j<3 &&
      admittance->basewave_s[i][j].name[0]
      ;j++);
    
    switch(j){
    case 3:
      admittance->method[i] = ADMITTANCE_SPLINE;
      admittance_scoefficient(admittance->basewave_s[i],admittance->response_s[i]);
      break;
    
    case 2:
      admittance->method[i] = ADMITTANCE_LINEAR;
      ///setting the linear list if necessary
      for(j=0;j<2;j++){
        admittance->basewave_l[i][j]=admittance->basewave_s[i][j];
        admittance->windex_l[i][j]=admittance->windex_s[i][j];
        }
      
      admittance_lcoefficient(admittance->basewave_l[i],admittance->response_l[i]);
      break;
    
    default:
      admittance->method[i] = ADMITTANCE_UNAVAILABLE;
      continue;
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *admittance_mts_check(admittance_t *admittance,const spectrum_t &WaveList,int *available)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int *method=new int [nspecies];

  for(size_t k = 0; k < 3; k++) {
    for(size_t l = 0; l < 3; l++) admittance->windex_s[k][l] = -1;
    admittance->count_s[k]=0;
    for(size_t l = 0; l < 2; l++) admittance->windex_l[k][l] = -1;
    admittance->count_l[k]=0;
    }

  admittance_mts_init(admittance, 0);

  for(size_t k = 0; k < 3; k++) {
/*------------------------------------------------------------------------------
    spline admittance */
    for(size_t l = 0; l < 3; l++) {
      int j=WaveList.wave_index(admittance->basewave_s[k][l].name);
      if(j==-1)
        continue;
      if(!available or available[j]==1) {
        admittance->windex_s[k][l] = j;
        admittance->count_s[k]++;
        }
      }
/*------------------------------------------------------------------------------
    linear admittance */
    for(size_t l = 0; l < 2; l++) {
      int j=WaveList.wave_index(admittance->basewave_l[k][l].name);
      if(j==-1)
        continue;
      if(!available or available[j]==1) {
        admittance->windex_l[k][l] = j;
        admittance->count_l[k]++;
        }
      }
    }

  for(size_t k = 0; k < 3; k++) {
    if(admittance->count_s[k] != 3) {
      admittance->method[k] = ADMITTANCE_LINEAR;
      }
    else
      admittance->method[k] = ADMITTANCE_SPLINE;
    }

  for(size_t k = 0; k < 3; k++) {
    if((admittance->count_l[k] != 2) && (admittance->method[k] == ADMITTANCE_LINEAR)) {
      admittance->method[k] = ADMITTANCE_UNAVAILABLE;
      }
    method[k]=admittance->method[k];
    }

#if VERBOSE > 2
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
      }
    }
#endif

  return (method);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void admittance_mts_compute(const admittance_t &admittance,spectrum_t WaveList,hconstant_t ztide,int *available, float mask, bool force, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j, k, l;
  complex<float> elevation_s[3][3];
  complex<float> elevation_l[3][2];
  float amp, pha;
//   bool *used=new bool[WaveList.n]:

  for(j = 0; j < WaveList.n; j++) {
    amp = ztide.a[j];
    pha = ztide.G[j] * d2r;
//     used[j]=false;
/*------------------------------------------------------------------------------
    loop on the 3 tidal species*/
    for(k = 0; k < 3; k++) {
/*------------------------------------------------------------------------------
      loop on the 3 base wave*/
      for(l = 0; l < 3; l++) {
        if(strcmp(WaveList.waves[j].name, admittance.basewave_s[k][l].name) == 0) {
          elevation_s[k][l] = polar<float>(amp,-pha);
//           used[j]=true;
          }
        }
/*------------------------------------------------------------------------------
      loop on the 2 base wave*/
      for(l = 0; l < 2; l++) {
        if(strcmp(WaveList.waves[j].name, admittance.basewave_l[k][l].name) == 0) {
          elevation_l[k][l] = polar<float>(amp,-pha);
//           used[j]=true;
          }
        }
      }
    }

  for(j =0; j < WaveList.n; j++) {
    if(available[j]==1 and !force) continue;
    k = WaveList.waves[j].nT;
    switch (admittance.method[k]) {
      case ADMITTANCE_SPLINE:
        admittance_scompute(WaveList.waves[j], admittance.basewave_s[k], elevation_s[k], admittance.response_s[k],&amp, &pha);
        ztide.a[j] = (float) amp;
        ztide.G[j] = (float) pha;
        if (verbose==1) printf("admittance_mts_compute, %s deduced from spline admittance\n",WaveList.waves[j].name);
        break;

      case ADMITTANCE_LINEAR:
        admittance_lcompute(WaveList.waves[j], admittance.basewave_l[k], elevation_l[k], admittance.response_l[k], &amp, &pha);
        ztide.a[j] = (float) amp;
        ztide.G[j] = (float) pha;
        if (verbose==1) printf("admittance_mts_compute, %s deduced from linear admittance\n",WaveList.waves[j].name);
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

  void admittance_mts_terminate(admittance_t admittance)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  for(int k = 0; k< nspecies; k++){ // species
    // spline admittance functions case
    for(int l = 0; l < 3; l++){
      for(int i=0; i < 3; i++) admittance.response_s[k][l][i] = 0;
      }

    // linear admittance functions case
    for(int l = 0; l < 2; l++){
      for(int i = 0; i < 2; i++) admittance.response_l[k][l][i] = 0;
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void admittance_mts_sweightP(const admittance_t &admittance,const tidal_wave &w, double coef[3])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k = w.nT;
  admittance_sweight(w, admittance.basewave_s[k], admittance.response_s[k], coef);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void admittance_mts_lweightP(const admittance_t &admittance,const tidal_wave &w, double coef[2])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k = w.nT;
  admittance_lweight(w, admittance.basewave_l[k], admittance.response_l[k], coef);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int admittance_mts_verify(const tidal_wave &w, const int *method)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k = w.nT;
  if(k > 2) return(-1);

  float Ap = w.Ap;
  if(method[k] == ADMITTANCE_UNAVAILABLE)
    return (-1);
  else
    if(Ap < 1.e-4){
      return (-1);
      }
    else{
      return (0);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int admittance_mts_verify(const admittance_t &admittance,const spectrum_t & WaveList,int *available,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j,status;
  int notOk=0;
  
  for(j =0; j < WaveList.n; j++) {
    if(available[j]==1) continue;
    tidal_wave *wave=&WaveList.waves[j];
    status=admittance_mts_verify(*wave,admittance.method);
    if(status!=0){
      notOk++;
      if(verbose) STDERR_BASE_LINE("Can not resolve safely wave %s\n",wave->name);
      }
    }
  
  return notOk;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int admittance_mts_error(const admittance_t &admittance,const spectrum_t& s,const int *deduce, vector<double>& error)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*-----------------------------------------------------------------------------
  compute linear combination of basis waves' error */

  double coef[3] = {0., 0., 0.};

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
 
  return(0);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void admittance_scoefficient(tidal_wave wave[3], double ri[3][3])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l;
  double delta[3] = {0, 0, 0}, cs[3] = {0, 0, 0}, sn[3] = {0, 0, 0};
  double matrix[9];

  for(k = 0; k < 3; k++){
    wave[k].init();
    }
  
  for(k = 0; k < 3; k++) {
    delta[k] = (wave[k].omega - wave[0].omega) * dph2cpd * 4. * M_PI;
    }

  for(k = 0; k < 3; k++) {
    cs[k] = cos(delta[k]);
    sn[k] = sin(delta[k]);
    }

  for(k = 0; k < 3; k++) {
    matrix[k] = cs[k];
    matrix[k + 3] = sn[k];
    matrix[k + 6] = 1.0;
    }

/* matrice factorisation */
  int neq = 3;
  int pivot[3], status = 0;
  status=poc_getrf(neq, matrix, pivot);

/*------------------------------------------------------------------------------
  impulses response (or matrix inverse) :
  
  
  h=a cos(f-f0) +b sin (f-f0) +c
  
  cos (f0-f0) sin(f0-f0) 1  a     h0
  cos (f1-f0) sin(f1-f0) 1  b  =  h1
  cos (f2-f0) sin(f2-f0) 1  c     h2

  cos (f0-f0) sin(f0-f0) 1  r[0][0]     1
  cos (f1-f0) sin(f1-f0) 1  r[0][1]  =  0
  cos (f2-f0) sin(f2-f0) 1  r[0][2]     0
  
  cos (f0-f0) sin(f0-f0) 1  r[1][0]     0
  cos (f1-f0) sin(f1-f0) 1  r[1][1]  =  1
  cos (f2-f0) sin(f2-f0) 1  r[1][2]     0

  cos (f0-f0) sin(f0-f0) 1  r[2][0]     0
  cos (f1-f0) sin(f1-f0) 1  r[2][1]  =  0
  cos (f2-f0) sin(f2-f0) 1  r[2][2]     1
  
  
  a= ri[0][0] h0 + ri[1][0] h1 + ri[2][0] h2
  b= ri[0][1] h0 + ri[1][1] h1 + ri[2][1] h2
  c= ri[0][2] h0 + ri[1][2] h1 + ri[2][2] h2
  
  
--------------------------------------------------------------------*/
  for(k = 0; k < 3; k++){
    for(l = 0; l < 3; l++){
      ri[k][l] = 0.0;
      }
    }

  for(k = 0; k < 3; k++){
    ri[k][k] = 1.0;
    }

  int nrhs = 1;
  for(k = 0; k < 3; k++)
    status=poc_getrs(neq, nrhs, matrix, pivot, ri[k]);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void admittance_lcoefficient(tidal_wave wave[2], double ri[2][2])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l;
  double delta[2] = {0, 0};
  double matrix[4] = {0, 0, 0, 0};

  for(k = 0; k < 2; k++){
    wave[k].init(); // omega [deg/h]
    }
  
  for(k = 0; k < 2; k++) {
    delta[k] = (wave[k].omega - wave[0].omega) * dph2cpd * 4 * M_PI;
    }

  for(k = 0; k < 2; k++) {
    matrix[k] = delta[k];
    matrix[k + 2] = 1.0;
//     STDERR_BASE_LINE("%s %g\n",wave[k].name,wave[k].omega);
//     STDERR_BASE_LINE("%g %g\n",matrix[k],matrix[k+2]);
    }

  int neq = 2;
  int pivot[2] = {0, 0}, status = -1;
  status=poc_getrf(neq, matrix, pivot);

  for(k = 0; k < 2; k++){
    for(l = 0; l < 2; l++){
      ri[k][l] = 0.0;
      }
    }

  for(k = 0; k < 2; k++){
    ri[k][k] = 1.0;
    }

  int nrhs = 1;
  for(k = 0; k < 2; k++)
    status=poc_getrs(neq, nrhs, matrix, pivot, ri[k]);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void admittance_sweight(const tidal_wave &w,const tidal_wave wave[3],const double ri[3][3], double coef[3])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**--------------------------------------------------------------------
 
 admittance coefficients, normalized by the potential's amplitude
 of the pivot constituents
 

 h = a cos(f-f0) +b sin (f-f0) +c

 a= ri[0][0] h0 + ri[1][0] h1 + ri[2][0] h2
 b= ri[0][1] h0 + ri[1][1] h1 + ri[2][1] h2
 c= ri[0][2] h0 + ri[1][2] h1 + ri[2][2] h2

  h=(ri[0][0] cos(f-f0) + ri[0][1] sin (f-f0) ri[0][2]) h0
   +(ri[1][0] cos(f-f0) + ri[1][1] sin (f-f0) ri[1][2]) h1
   +(ri[2][0] cos(f-f0) + ri[2][1] sin (f-f0) ri[2][2]) h2
   =c0 h0 +c1 h1 + c2 h2

--------------------------------------------------------------------*/

  double gg = (w.omega - wave[0].omega) * dph2cpd * 4. * M_PI;

  for(int k = 0; k < 3; k++) {
    if((ri[k][0]==0.) && (ri[k][1]==0.) && (ri[k][1]==0.)) {
      TRAP_ERR_EXIT(-1, "admittance coefficients not initialized\n");
      }
    coef[k] = ri[k][0] * cos(gg) + ri[k][1] * sin(gg) + ri[k][2];
    coef[k] /= wave[k].Ap;
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void admittance_lweight(const tidal_wave &w,const tidal_wave wave[2],const double ri[2][2], double coef[2])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------

  h = a (f-f0) +b

  a= ri[0][0] h0 + ri[1][0] h1
  b= ri[0][1] h0 + ri[1][1] h1

   h=(ri[0][0] (f-f0) + ri[0][1] ) h0
    +(ri[1][0] (f-f0) + ri[1][1] ) h1
    =c0 h0 +c1 h1

--------------------------------------------------------------------*/

  double gg = (w.omega - wave[0].omega) * dph2cpd  * 4. * M_PI;

  for(int k = 0; k < 2; k++) {
    coef[k] = ri[k][0] * gg + ri[k][1];
    coef[k] /= wave[k].Ap;
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void admittance_s_or_l_weight(const tidal_wave &w,const int n,const tidal_wave wave[3],const double ri[3][3],double c[3])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(n!=3){TRAP_ERR_EXIT(ENOEXEC,"should never reach here!");}
  
  admittance_sweight(w,wave,ri,c);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void admittance_s_or_l_weight(const tidal_wave &w,const int n,const tidal_wave wave[2],const double ri[2][2],double c[2])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(n!=2){TRAP_ERR_EXIT(ENOEXEC,"should never reach here!");}
  
  admittance_lweight(w,wave,ri,c);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<int n> void admittance_template_compute(const tidal_wave &w,const tidal_wave wave[n],const complex<float> elevation[n],const double ri[n][n],float *amp,float *pha)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  complex<float> h=0;
  double c[n];
  
  admittance_s_or_l_weight(w,n,wave,ri,c);
  
  for(i=0;i<n;i++){
    h+=(float)c[i]*elevation[i];
    }
  h*=w.Ap;
  
  *amp = abs(h);
  *pha = -arg(h) * r2d;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void admittance_lcompute(const tidal_wave &w,const tidal_wave wave[2],const complex<float> elevation[2],const double ri[2][2],float *amp,float *pha)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  admittance_template_compute(w,wave,elevation,ri,amp,pha);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void admittance_scompute(const tidal_wave &w,const tidal_wave wave[3],const complex<float> elevation[3],const double ri[3][3],float *amp,float *pha)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  admittance_template_compute(w,wave,elevation,ri,amp,pha);
}
