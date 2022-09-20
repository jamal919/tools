
/*******************************************************************************

  T-UGO tools, 2006-2013

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Sara Fleury        LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Detide altimetric data series. USED BY CTOH.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include "config.h"

#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include "tools-structures.h"

#include "fe.h"
#include "archive.h"
#include "poc-time.h"
#include "mgr.h"
#include "functions.h"
#include "map.h"
#include "map.def"
#include "xtrack-io.h"
#include "spectrum.h"
#include "filter.h"
#include "statistic.h"

#ifndef VERBOSE
#define VERBOSE 1
#endif

//#define OPEN_MP


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

#define __REPORT_THIS__(s,msg) do {cout << msg; (s << msg).flush();} while(0)
#define __REPORT_THIS_AND_EXIT__(s,msg) do {__REPORT_THIS__(s,msg); TRAP_ERR_EXIT(1,"exiting\n");} while(0)
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

class satellite_t{
private :
public :
  double repetitivity;
  double duration;
  double inclination;
  int    ncycles;
  date_t launch;              /** launch date                        */
  date_t start,termination;   /** operational start, operational end */
  date_t current;             /** last operational available records */

  satellite_t() {
    repetitivity =-1;
    duration     =-1;
    ncycles      =-1;
    }

  satellite_t(double dT, double T) {
    repetitivity =dT;
    duration     = T;
    ncycles      =(int) floor(1.+T/dT);
    }

  void destroy() {
    }
};
/*------------------------------------------------------------------------------
 T/P  : 343 cycles    Jason 1 : 239 cycles    Jason 2 : 127 cycles     709 cycles = 19.2 years
 T/P  : 124 cycles    Jason 1 : 100 cycles                             224 cycles =  6.1 years
 ERS1 :  26 cycles    ERS2    :  77 cycles   EnviSAT  :  84 cycles     187 cycles = 17.9 years
 GFO  : 186 cycles                                                     186 cycles =  8.7 years
 */
const satellite_t topex ( 9.91564, 365.25*19.2);
const satellite_t topex2( 9.91564, 365.25* 6.1);
const satellite_t ers   (35.00000, 365.25*17.9);
const satellite_t gfo   (17.05064, 365.25* 8.7);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  satellite_t initilize_satellite(const string & satellite)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  satellite_t sat;
/*    int satellite_id;*/
  
  if( satellite == "topex") {
/*    satellite_id = TOPEX;*/
    sat=topex;
    }
  else if( satellite == "topex2") {
/*    satellite_id = TOPEX;*/
    sat=topex2;
    }
  else if( satellite == "ers") {
/*    satellite_id = TOPEX;*/
    sat=ers;
    }
  else if( satellite == "gfo") {
/*    satellite_id = TOPEX;*/
    sat=gfo;
    }
  else {
    ;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void missions_info(string spectrum)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  double Trepet;
  bool Z0=false;
  spectrum_t AnalysisList = spectrum_init_ref(spectrum.c_str(), Z0);
  
  double duration=18*365;
  Trepet=35.;
  spectrum_define_alias(&AnalysisList, Trepet, stdout);
  spectrum_checkRayleigh(AnalysisList, duration, 0, 0, stdout);
  
/*------------------------------------------------------------------------------
  check ERS */
  Trepet=35.;
  duration=18*365;
  printf("\n");
  printf("#################################################################\n");
  printf("ERS/EnviSat: T=%lf duration=%lf\n",Trepet,duration);
  spectrum_reset(AnalysisList, "COASTAL", Trepet, 0);
  spectrum_define_alias(&AnalysisList, Trepet, stdout);
  spectrum_checkRayleigh(AnalysisList, duration, 0, 0, stdout);
/*------------------------------------------------------------------------------
  check ERS Xovers */
  Trepet=35./2;
  duration=18*365;
  printf("\n");
  printf("#################################################################\n");
  printf("ERS/EnviSat Xovers: T=%lf duration=%lf\n",Trepet,duration);
  spectrum_reset(AnalysisList, "COASTAL", Trepet, 0);
  spectrum_define_alias(&AnalysisList, Trepet, stdout);
  spectrum_checkRayleigh(AnalysisList, duration, 0, 0, stdout);
  
/*------------------------------------------------------------------------------
  check GFO */
  Trepet=17.05064;
  duration=8*365;
  printf("\n");
  printf("#################################################################\n");
  printf("GFO: T=%lf duration=%lf\n",Trepet,duration);
  spectrum_reset(AnalysisList, "COASTAL", Trepet, 0);
  spectrum_define_alias(&AnalysisList, Trepet, stdout);
  spectrum_checkRayleigh(AnalysisList, duration, 0, 0, stdout);
/*------------------------------------------------------------------------------
  check GFO Xovers */
  Trepet=17.05064/2.0;
  duration=8.7*365;
  printf("\n");
  printf("#################################################################\n");
  printf("GFO Xovers: T=%lf duration=%lf\n",Trepet,duration);
  spectrum_reset(AnalysisList, "COASTAL", Trepet, 0);
  spectrum_define_alias(&AnalysisList, Trepet, stdout);
  spectrum_checkRayleigh(AnalysisList, duration, 0, 0, stdout);
  
/*------------------------------------------------------------------------------
  check TOPEX */
  Trepet=9.91564;
  duration=19*365;
  printf("\n");
  printf("#################################################################\n");
  printf("TOPEX: T=%lf duration=%lf\n",Trepet,duration);
  spectrum_reset(AnalysisList, "COASTAL", Trepet, 0);
  spectrum_define_alias(&AnalysisList, Trepet, stdout);
  spectrum_checkRayleigh(AnalysisList, duration, 0, 0, stdout);
  
/*------------------------------------------------------------------------------
  check TOPEX interleaved */
  Trepet=9.91564;
  duration=6*365;
  printf("\n");
  printf("#################################################################\n");
  printf("TOPEX interleaved: T=%lf duration=%lf\n",Trepet,duration);
  spectrum_reset(AnalysisList, "COASTAL", Trepet, 0);
  spectrum_define_alias(&AnalysisList, Trepet, stdout);
  spectrum_checkRayleigh(AnalysisList, duration, 0, 0, stdout);
/*------------------------------------------------------------------------------
  check TOPEX interleaved */
  Trepet=9.91564/2.;
  duration=6*365;
  printf("\n");
  printf("#################################################################\n");
  printf("TOPEX interleaved Xovers: T=%lf duration=%lf\n",Trepet,duration);
  spectrum_reset(AnalysisList, "COASTAL", Trepet, 0);
  spectrum_define_alias(&AnalysisList, Trepet, stdout);
  spectrum_checkRayleigh(AnalysisList, duration, 0, 0, stdout);
  TRAP_ERR_EXIT(-1,"exiting\n");
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  spectrum_t atlas_define(string atlas_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  spectrum_t atlas;
  int i = 0;

  if(atlas_name=="GLORYS-v2"){
    atlas.n = atlas.nmax = 2;
    atlas.waves = new tidal_wave[atlas.n];

    atlas.waves[i++] = wSa;
    atlas.waves[i++] = wSsa;

    initialize_omega(&atlas);
    return(atlas);
    }
  
  if(atlas_name=="FES2004"){
    atlas.n = atlas.nmax = 14;
    atlas.waves = new tidal_wave[atlas.n];

    atlas.waves[i++] = wMm;  // FES2004 atlas
    atlas.waves[i++] = wMf;  // FES2004 atlas
    atlas.waves[i++] = wMtm; // FES2004 atlas
    atlas.waves[i++] = wMSqm;// FES2004 atlas

    atlas.waves[i++] = wQ1;  // FES2004 atlas
    atlas.waves[i++] = wO1;  // FES2004 atlas
    atlas.waves[i++] = wP1;  // FES2004 atlas
    atlas.waves[i++] = wK1;  // FES2004 atlas
    atlas.waves[i++] = wS1;  // FES2004 atlas

    atlas.waves[i++] = w2N2; // FES2004 atlas
    atlas.waves[i++] = wN2;  // FES2004 atlas
    atlas.waves[i++] = wM2;  // FES2004 atlas
    atlas.waves[i++] = wS2;  // FES2004 atlas
    atlas.waves[i++] = wK2;  // FES2004 atlas

    initialize_omega(&atlas);
    return(atlas);
    }

  if(atlas_name=="FES2012"){
    atlas.n = atlas.nmax = 13;
    atlas.waves = new tidal_wave[atlas.n];

    atlas.waves[i++] = wMm;  // FES2012 atlas
    atlas.waves[i++] = wMf;  // FES2012 atlas
    atlas.waves[i++] = wMtm; // FES2012 atlas
//    atlas.waves[i++] = wMSqm;// FES2012 atlas

    atlas.waves[i++] = wQ1;  // FES2012 atlas
    atlas.waves[i++] = wO1;  // FES2012 atlas
    atlas.waves[i++] = wP1;  // FES2012 atlas
    atlas.waves[i++] = wK1;  // FES2012 atlas
    atlas.waves[i++] = wS1;  // FES2012 atlas

    atlas.waves[i++] = w2N2; // FES2012 atlas
    atlas.waves[i++] = wN2;  // FES2012 atlas
    atlas.waves[i++] = wM2;  // FES2012 atlas
    atlas.waves[i++] = wS2;  // FES2012 atlas
    atlas.waves[i++] = wK2;  // FES2012 atlas

    initialize_omega(&atlas);
    return(atlas);
    }

  if(atlas_name=="GOT4.7"){
    atlas.n = atlas.nmax = 9;
    atlas.waves =  new tidal_wave[atlas.n];

    atlas.waves[i++] = wQ1;  // GOT4.7 atlas
    atlas.waves[i++] = wO1;  // GOT4.7 atlas
    atlas.waves[i++] = wP1;  // GOT4.7 atlas
    atlas.waves[i++] = wK1;  // GOT4.7 atlas

    atlas.waves[i++] = wN2;  // GOT4.7 atlas
    atlas.waves[i++] = wM2;  // GOT4.7 atlas
    atlas.waves[i++] = wS2;  // GOT4.7 atlas
    atlas.waves[i++] = wK2;  // GOT4.7 atlas

    atlas.waves[i++] = wM4;  // GOT4.7 atlas

    initialize_omega(&atlas);
    return(atlas);
    }

  if(atlas_name=="MEDSEA"){
    initialize_omega(&atlas);
    return(atlas);
    }

  if(atlas_name=="NEA"){
    atlas.n = atlas.nmax = 19;
    atlas.waves = new tidal_wave[atlas.n];

    atlas.waves[i++] = wMm;  // T-UGOm 2D atlas
    atlas.waves[i++] = wMf;  // T-UGOm 2D atlas
    atlas.waves[i++] = wMtm; // T-UGOm 2D atlas
    atlas.waves[i++] = wMSqm;// T-UGOm 2D atlas

    atlas.waves[i++] = wQ1;  // T-UGOm 2D atlas
    atlas.waves[i++] = wO1;  // T-UGOm 2D atlas
    atlas.waves[i++] = wP1;  // T-UGOm 2D atlas
    atlas.waves[i++] = wK1;  // T-UGOm 2D atlas
    atlas.waves[i++] = wS1;  // FES98 atlas

    atlas.waves[i++] = w2N2;// T-UGOm 2D atlas
    atlas.waves[i++] = wN2; // T-UGOm 2D atlas
    atlas.waves[i++] = wM2; // T-UGOm 2D atlas
    atlas.waves[i++] = wS2; // T-UGOm 2D atlas
    atlas.waves[i++] = wK2; // T-UGOm 2D atlas

    atlas.waves[i++] = wN4;  // T-UGOm 2D atlas
    atlas.waves[i++] = wMN4; // T-UGOm 2D atlas
    atlas.waves[i++] = wM4;  // T-UGOm 2D atlas
    atlas.waves[i++] = wMS4; // T-UGOm 2D atlas
    atlas.waves[i++] = wS4;  // T-UGOm 2D atlas

    initialize_omega(&atlas);
    return(atlas);
    }

  if(atlas_name=="PERSIAN"){
    atlas.n = atlas.nmax = 6;
    atlas.waves =  new tidal_wave[atlas.n];

    //     atlas.waves[i++] = wMm;  // T-UGOm 2D atlas
    //     atlas.waves[i++] = wMf;  // T-UGOm 2D atlas
    //     atlas.waves[i++] = wMtm; // T-UGOm 2D atlas
    //     atlas.waves[i++] = wMSqm;// T-UGOm 2D atlas

    //     atlas.waves[i++] = wQ1;  // T-UGOm 2D atlas
    atlas.waves[i++] = wO1;  // T-UGOm 2D atlas
    //     atlas.waves[i++] = wP1;  // T-UGOm 2D atlas
    atlas.waves[i++] = wK1;  // T-UGOm 2D atlas
    //     atlas.waves[i++] = wS1;  // FES98 atlas

    //     atlas.waves[i++] = w2N2;// T-UGOm 2D atlas
    atlas.waves[i++] = wN2;// T-UGOm 2D atlas
    atlas.waves[i++] = wM2;// T-UGOm 2D atlas
    atlas.waves[i++] = wS2;// T-UGOm 2D atlas
    //     atlas.waves[i++] = wK2;// T-UGOm 2D atlas

    //     atlas.waves[i++] = wN4;// T-UGOm 2D atlas
    //     atlas.waves[i++] = wMN4;// T-UGOm 2D atlas
    atlas.waves[i++] = wM4;// T-UGOm 2D atlas
    //     atlas.waves[i++] = wMS4;// T-UGOm 2D atlas
    //     atlas.waves[i++] = wS4;// T-UGOm 2D atlas

    initialize_omega(&atlas);
    return(atlas);
    }

  // default
  atlas.n = atlas.nmax = 0;
  atlas.waves =  0;
  return(atlas);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int atlas_expand(spectrum_t *atlas, string atlas_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  if(atlas_name=="GLORYS-v2"){

    return(0);
    }
  
  else if(atlas_name=="FES2004"){

    tide_addwave(atlas, wMSm);  // spline admittance
    tide_addwave(atlas, wMSf);  // spline admittance
    tide_addwave(atlas, wMStm); // spline admittance
    tide_addwave(atlas, wMqm);  // spline admittance
    tide_addwave(atlas, wSa);   // equilibrium
    tide_addwave(atlas, wSsa);  // equilibrium
    tide_addwave(atlas, w2Q1);  // spline admittance
    tide_addwave(atlas, wSig1); // spline admittance
    tide_addwave(atlas, wRo1);  // spline admittance
    tide_addwave(atlas, wKi1);  // spline admittance
    tide_addwave(atlas, wPi1);  // spline admittance
    tide_addwave(atlas, wPsi1); // spline admittance
    tide_addwave(atlas, wPhi1); // spline admittance
    tide_addwave(atlas, wTta1); // spline admittance
    tide_addwave(atlas, wJ1);   // spline admittance
    tide_addwave(atlas, wOO1);  // spline admittance
    tide_addwave(atlas, wR2);   // spline admittance
    tide_addwave(atlas, wLa2);  // spline admittance
    tide_addwave(atlas, wL2);   // spline admittance
    tide_addwave(atlas, wT2);   // spline admittance
    tide_addwave(atlas, wNu2);  // spline admittance
    tide_addwave(atlas, wMu2);  // spline admittance
    tide_addwave(atlas, wKJ2);  // spline admittance

    initialize_omega(atlas);
    return(0);
    }

  if(atlas_name=="GOT4.7"){

    tide_addwave(atlas, w2Q1); // spline admittance
    tide_addwave(atlas, wSig1); // spline admittance
    tide_addwave(atlas, wRo1); // spline admittance
    tide_addwave(atlas, wKi1); // spline admittance
    tide_addwave(atlas, wPi1); // spline admittance
    tide_addwave(atlas, wPsi1); // spline admittance
    tide_addwave(atlas, wPhi1); // spline admittance
    tide_addwave(atlas, wTta1); // spline admittance
    tide_addwave(atlas, wJ1);  // spline admittance
    tide_addwave(atlas, wOO1); // spline admittance
    tide_addwave(atlas, wR2);  // spline admittance
    tide_addwave(atlas, wLa2); // spline admittance
    tide_addwave(atlas, wL2);  // spline admittance
    tide_addwave(atlas, wT2);  // spline admittance
    tide_addwave(atlas, wNu2); // spline admittance
    tide_addwave(atlas, wMu2); // spline admittance
    tide_addwave(atlas, wKJ2); // spline admittance

    initialize_omega(atlas);
    return(0);
    }

  if(atlas_name=="MEDSEA"){
    initialize_omega(atlas);
    return(1);
    }

  if(atlas_name=="NEA"){

    tide_addwave(atlas, wMSm); // spline admittance
    tide_addwave(atlas, wMSf); // spline admittance
    tide_addwave(atlas, wMStm);// spline admittance
    tide_addwave(atlas, wMqm); // spline admittance
    tide_addwave(atlas, wSsa); // equilibrium
    tide_addwave(atlas, w2Q1); // spline admittance
    tide_addwave(atlas, wSig1); // spline admittance
    tide_addwave(atlas, wRo1); // spline admittance
    tide_addwave(atlas, wKi1); // spline admittance
    tide_addwave(atlas, wPi1); // spline admittance
    tide_addwave(atlas, wPsi1); // spline admittance
    tide_addwave(atlas, wPhi1); // spline admittance
    tide_addwave(atlas, wTta1); // spline admittance
    tide_addwave(atlas, wJ1);  // spline admittance
    tide_addwave(atlas, wOO1); // spline admittance
    tide_addwave(atlas, wR2);  // spline admittance
    tide_addwave(atlas, wLa2); // spline admittance
    tide_addwave(atlas, wL2);  // spline admittance
    tide_addwave(atlas, wT2);  // spline admittance
    tide_addwave(atlas, wNu2); // spline admittance
    tide_addwave(atlas, wMu2); // spline admittance
    tide_addwave(atlas, wKJ2); // spline admittance

    initialize_omega(atlas);
    return(1);
    }

  if(atlas_name=="PERSIAN"){

    initialize_omega(atlas);
    return(1);
  }// PERSIAN

  return(1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

hconstant_t extended_atlas_constants(string atlasName, spectrum_t base, spectrum_t *extended, double x, double y, hconstant_t elevation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int anomaly;

/* *----------------------------------------------------------------------
  must be consistent with hconstant_t elevation (potentially bugged)*/ // call to same function!
//  *s = atlas_define(atlasName); // LR: should not be necessary

/*------------------------------------------------------------------------------
  extend prediction spectrum*/
  extended->duplicate(base);
  status=atlas_expand(extended, atlasName);
  
  if(extended->n==0) {
    anomaly=1;
    }

  hconstant_t ele;
  ele.a = new float[extended->n];
  ele.G = new float[extended->n];
  ele.size=extended->n;
  
  int *available = new int[extended->n];
  for(int k = 0; k < base.n; k++){
    ele.a[k] = elevation.a[k];
    ele.G[k] = elevation.G[k];
    available[k] = 1;
    }
  delete [] elevation.a;
  delete [] elevation.G;
  for(int k = base.n; k < extended->n; k++){
    available[k] = 0;
    }

/*------------------------------------------------------------------------------
  infer missing constituents using admmittance*/
//#ifdef OPEN_MP
  admittance_t admittance;
  status=admittance_mts_init(&admittance, 0);
  int *method = admittance_mts_check(&admittance, *extended);
  admittance_mts_compute(admittance, *extended, ele, available,0.0,false,0);
  admittance_mts_terminate(admittance);
// #else
//   status=admittance_init(0);
//   int *method = admittance_check(*extended, 0);
//   admittance_compute(*extended, ele, available);
//   admittance_terminate();
// #endif
  delete[] method;

/*------------------------------------------------------------------------------
  infer missing constituents using equilibrium */
  for (int k = base.n; k < extended->n; k++) {
    if( (available[k] == 0) && (strcmp(extended->waves[k].name,"Ssa") == 0) ) {
      tidal_equilibrium(extended->waves[k], x, y, &(ele.a[k]), &(ele.G[k]));
      }
    if( (available[k] == 0) && (strcmp(extended->waves[k].name,"Sa") == 0) ) {
      tidal_equilibrium(extended->waves[k], x, y, &(ele.a[k]), &(ele.G[k]));
      }
    }

  delete [] available;

  ele.s=extended;
  
  return(ele);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int restore(const spectrum_t& AnalysisList, mgr_data_t *data, const spectrum_t& PredictionList, const hconstant_t& prior)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{ /** extremely unsafe: size of 'data' array must be length of list 'AnalysisList' */

  size_t j = -1;
  int anomaly=0;
  
  complex<float> zPrior = 0, zAnalysis = 0;
  
  for(size_t w = 0; w < AnalysisList.n; w++) {
    if( (j = index_from_name(AnalysisList,AnalysisList.waves[w].name)) == w ){ // internal test
      zAnalysis = polar(data[j].amp,-data[j].phi*d2r); // analysis() returns phase lags in 0 < deg < 360
      }
    if( (j = index_from_name(PredictionList,AnalysisList.waves[w].name)) != -1 ){
      zPrior = polar(static_cast<double> (prior.a[j]),static_cast<double> (-prior.G[j])*d2r); // atlas_constants() returns phase lags in 0 < deg < 360
      zAnalysis += zPrior;
      } // else use analysed values
    
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

  int check_timeseries(double Trepet, int n, double *time, double *z, double *residuals)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *debugFile = fopen("debug.time_series.dat", "w");
  for(size_t i = 0; i < n; i++){
    fprintf(debugFile, "%12.6f %9.6f %9.6f\n", time[i], z[i], residuals[i]);
    }
  fclose(debugFile);
  for(size_t i = 0; i < n; i++){
    time[i] *= 24;
    }
  fourier("debug.fft_z.dat", z, 99.9999, n, Trepet,0);
  fourier("debug.fft_res.dat", residuals, 99.9999, n, Trepet,0);
  for(size_t i = 0; i < n; i++){
    time[i] /= 24;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int register_analysis(int index, vector<mgr_t> & mgr, mgr_data_t *data,const char *input, const serie_t & currentRecord, double Trepet, const spectrum_t & analysed, int nValid, double *time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
    store harmonic constants
    
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
    mgr_t gauge;
    
    gauge.number = index;
    
//     strcpy(gauge.name, input);
    sprintf(gauge.name, "%s-%4.4d",input,index);
    
    if( (Trepet == 9.9156) || (Trepet == 17.048860) || (Trepet == 35.) ){
      strcpy(gauge.origine, "X-TRACK SLA");
      strcpy(gauge.validation, "Statistical QC@CTOH");
      }
    else {
      strcpy(gauge.origine, "Unknown");
      strcpy(gauge.validation, "Unknown");
      }
    gauge.data = data;
    gauge.nwave = analysed.n;
    gauge.loc.units = strdup("degrees");
    gauge.loc.lon = currentRecord.data[0].lon;
    gauge.loc.lat = currentRecord.data[0].lat;
    gauge.loc.depth = 0.;
    gauge.mindex = 1;
    gauge.duree = nValid;
    date_t start = poctime_getdatecnes(time[0], 'd');
    date_t final = poctime_getdatecnes(time[nValid - 1], 'd');
    sprintf(gauge.debut,"%02d/%02d/%4d",start.day,start.month,start.year);
    sprintf(gauge.fin,  "%02d/%02d/%4d",final.day,final.month,final.year);
    mgr.push_back(gauge);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int detide(vector<serie_t> track, int record, hconstant_t *OceanTide_base, float atlas_fmask, int nValid, double *z, double *time, int account_for_nodal_args)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//       if(OceanTide_base[record].hasMaskValue(atlas_fmask) == true){
//         printf(" skipped\n");
//         for(int i = 0, j = 0; i < track.count; i++){
//           track.data[i].values[X_TRACK_SERIE_T_TIDE_ANALYSIS] = track.mask;
//           track.data[i].values[X_TRACK_SERIE_T_RESIDUAL]      = track.mask;
//           }
//         return(-1);
//         }
//
//       double x = track.data[0].lon;
//       double y = track.data[0].lat;
//       double sla;
//       int anomaly=0;
// /*------------------------------------------------------------------------------
//       used to be multi-thread unsafe - fixed */
//       OceanTide_extended = extended_atlas_constants(atlasName, OceanPredictionBaseList, &OceanPredictionExtendedList, x, y, OceanTide_base[record]);
//
//       double *tides = new double[nValid];
// //      (void) harmonic_predictionTS(tides, time, nValid, OceanPredictionExtendedList, OceanTide_extended.a, OceanTide_extended.G, account_for_nodal_args);
//       (void) harmonic_predictionTS(tides, time, nValid, OceanTide_extended, account_for_nodal_args);
//
//       for(size_t j = 0 ; j < nValid ; j++) {
//         sla=z[j];
//         z[j] -= tides[j];
//         if(abs(z[j])> 1.e+03) {
//           anomaly=1;
//           printf("alert %d: t=%.6f sla=%.6f detided sla=%.6f tides=%.6f %.6f  (FES2004) %.6f (GOT4.7)\n",record,time[j], sla, z[j], tides[j],
//                track.data[j].values[X_TRACK_SERIE_T_TIDE_FES2004],
//                track.data[j].values[X_TRACK_SERIE_T_TIDE_GOT4_7]);
//           }
//         }
//       delete [] tides;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : argv[0]
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    "NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s OPTIONS input\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Detide altimetric data series from X-TRACK as input. USED BY CTOH.\n"
    "  May also produce a residual .mgr.stat file.\n"
    "\n"
    "OPTIONS :\n"
    "  --help,-h  Show this help and exit.\n"
    "  --info  show Rayleigh for several misssions and exit\n"
    "  -gnu=yes  put ascii data (time, height, 11 DAC parameters, 12 tidal atlas parameters, residual [see X_TRACK_SERIE_T_* in xtrack-io.h]) in a .gnu file for use with, e.g., gnuplot. Disabled by -np.\n"
    "  -r  followed by the repetitivity of measurements in days\n"
    "  -satellite  followed by the name of the satellite (UNUSED!)\n"
    "  -spectrum  name of spectrum\n"
    "  -latmin  minimum of absolute value of latitude: exclude lower latitudes\n"
    "    Defaults are /home/softs/data/climatology/GLORYS-v2/WAVE.GLORYS.nc smoothed_a smoothed_G\n"
    "  -np  followed by the number of cores. Default: 1. Anything above 1 will enable OpenMP parallelisation and will disable the output of the residual. Giving -1 will take all the available cores on your machine. See also the ENVIRONMENT section below.\n"
    "  --dac=no  do not apply GDR DAC correction\n"
    "  --tide=no  do not apply GDR tide correction\n"
    "  -o : followed by 3 parameters for ocean tide: convention, amplitude variable name, phaselag variable name\n"
    "  -l : followed by 3 parameters for load tide: convention, amplitude variable name, phaselag variable name. ONLY FOR GEOCENTRIC SERIES!\n"
    "  -ogcm_parameters : followed by 3 parameters for GLORYS corrections: convention, amplitude variable name, phaselag variable name.\n"
    "  -cm : use if atlases are in cm\n"
    "  -percent,-m  following by the minimum number of elements in the time series. Default: 160.\n"
    "  --auto-parse  use matrix instead of Rayleigh\n"
    "  -verbose : followed by a number, that, if >=2, will enable printing some stuff\n"
    "  --debug  more verbose for debugging\n"
    "  -s  followed by the start date in dd/mm/yyyy format\n"
    "  -e  followed by the end date in dd/mm/yyyy format\n"
    "  -a  followed by one of these atlas names for spectrum expansion: GLORYS-v2 FES2004 FES2012 GOT4.7 MEDSEA NEA PERSIAN\n"
    "  -f  followed by one of these tide gauge formats: ASCII (Default) NETCDF\n"
    "  -b  followed by one of these error budgets: psd. Anything else, like formal, the default, does nothing in particular.\n"
    );
  print_OPENMP_help(prog_name); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  ofstream log("altimetry-detidor.log");
  string logged;

  fct_echo(argc,argv);

  // init program args default values
  const char *input = 0;
  string atlasName = "";
  const char *oceanTide_convention = 0, *oceanTide_ampName=0, *oceanTide_phaName=0;
  const char *loading_convention = 0, *loading_ampName=0, *loading_phaName=0;
  double scaleAtlas=1.;
  const char *ogcm_convention = 0, *ogcm_amplitude = 0, *ogcm_phaselag = 0;
  string satellite="undefined", spectrum="COASTAL";
  bool atlas_detiding = false, gdr_detiding = true;
  bool deloading  = false;
  bool dealiasing = false, auto_parse = false;
  int use_ib   = 1;
  double Trepet = 0.;
  string mgrFormat = "ASCII";
  date_t reduced_start = NADate, reduced_final = NADate;
  string error_budget = "formal";
  int nmesMin = 160;
  int nRequestedProcs=1;
  bool save_gnu = false;
  statistic_t stats_input, stats_residuals;
  double maxratio=0.3;
  double latmin=0.0;
  double averaging_time=0.0;
  bool recompose=true;

  int n = 1;
  const char *keyword = 0;
  FILE *debuglog=0, *foutname;
  int verbose=0;
  
  
  while (n < argc) {
    keyword = argv[n];
    switch (keyword[0]) {
    case '-':
      if(strcmp(keyword,"--help")==0 || strcmp(keyword,"-h")==0) {
        print_help(argv[0]);
        exit(0);
        }
      if(strcmp(keyword,"--info")==0) {
        missions_info(spectrum);
        }
      if(strcmp(keyword,"--auto-parse")==0) {
        auto_parse=true;
        n++;
        break;
        }
      if(strcmp(keyword,"--recompose=no")==0) {
        recompose=false;
        n++;
        break;
        }
      if(strcmp(keyword,"-gnu=yes")==0) {
        save_gnu = true;
        n++;
        break;
        }
      if(strcmp(keyword,"-satellite")==0) {
        satellite = argv[n+1];
        n++;
        n++;
        break;
        }
      if(strcmp(keyword,"-spectrum")==0) {
        spectrum = argv[n+1];
        n++;
        n++;
        break;
        }
      if(strcmp(keyword,"-verbose")==0) {
        sscanf(argv[n+1],"%d",&verbose);
        n++;
        n++;
        break;
        }
     if(strcmp(keyword,"-latmin")==0) {
        sscanf(argv[n+1],"%lf",&latmin);
        n++;
        n++;
        break;
        }
      if(strcmp(keyword,"-averaging")==0) {
        int nitems __attribute__((unused))= sscanf(argv[n+1],"%lf",&averaging_time);
        n++;
        n++;
        break;
        }
      if(strcmp(keyword,"--debug")==0) {
        debuglog=stdout;
        n++;
        break;
        }
      if(strcmp(keyword,"-ogcm_parameters")==0) {
        ogcm_convention = argv[n+1];
        ogcm_amplitude  = argv[n+2];
        ogcm_phaselag   = argv[n+3];
        n++;
        n+=3;
        dealiasing = true;
        break;
        }
      if(strcmp(keyword,"-cm")==0) {
        scaleAtlas=0.01;
        n++;
        break;
        }
      switch (keyword[1]) {

/*------------------------------------------------------------------------------
      starting date for reduced analysis*/
      case 's' :{
        const char *s = argv[n+1];
        n++;
        n++;
        sscanf(s, "%d/%d/%d", &reduced_start.day, &reduced_start.month, &reduced_start.year);
        reduced_start.second = 0.;
        break;
        }

/*------------------------------------------------------------------------------
      final date for reduced analysis*/
      case 'e' :{
        const char *s = argv[n+1];
        n++;
        n++;
        sscanf(s,"%d/%d/%d", &reduced_final.day, &reduced_final.month, &reduced_final.year);
        reduced_final.second = 0.;
        break;
        }

/*------------------------------------------------------------------------------
      generic name for tidal atlas*/
      case 'a' :
        atlasName = argv[n+1];
        n++;
        n++;
        break;

/*------------------------------------------------------------------------------
      ocean tide atlases*/
      case 'o' :
        oceanTide_convention = argv[n+1];
        oceanTide_ampName    = argv[n+2];
        oceanTide_phaName    = argv[n+3];
        n++;
        n+=3;
        atlas_detiding=true;
        break;

/*------------------------------------------------------------------------------
      load tide atlases*/
      case 'l' :
        loading_convention = argv[n+1];
        loading_ampName =  argv[n+2];
        loading_phaName =  argv[n+3];
        n++;
        n+=3;
        deloading=true;
        break;

/*------------------------------------------------------------------------------
      format for mgr output*/
      case 'f' :
        mgrFormat = argv[n+1];
        n++;
        n++;
        break;

/*------------------------------------------------------------------------------
      error budget option*/
      case 'b' :
        error_budget = argv[n+1];
        n++;
        n++;
        break;
        
/*------------------------------------------------------------------------------
      repeat period in days*/
      case 'r' :{
        Trepet = atof(argv[n+1]);
        n++;
        n++;
        break;
        }

      default:
        if(strcmp(keyword,"-np")==0) {
          nRequestedProcs = atoi(argv[n+1]);;
          n++;
          n++;
          break;
          }
        else if (strcmp(keyword,"--dac=no")==0) {
          use_ib = 0;
          n++;
          break;
          }
        else if (strcmp(keyword,"--tide=no")==0) {
          gdr_detiding = false;
          n++;
          break;
          }
/*------------------------------------------------------------------------------
        minimum size of time series */
        else if(strcmp(keyword,"-percent")==0 || strcmp(keyword,"-m")==0) {
          nmesMin = atoi(argv[n+1]);
          n++;
          n++;
          break;
          }
        else {
          asprintf(logged="","ERROR : redundant option %s\n",keyword);
          __REPORT_THIS_AND_EXIT__(log, logged);
          }
        }
      break;

/*------------------------------------------------------------------------------
    input file (X-TRACK SLA ref) */
    default:
      if(input == NULL) {
        input = argv[n];
        n++;
        }
      else {
        asprintf(logged="","ERROR : redundant option %s\n",keyword);
        __REPORT_THIS_AND_EXIT__(log, logged);
        }
      break;
      }
    
    }

  if(input == NULL) {
    printf("*** no input provided ***\n");
    print_help(argv[0]);
    exit(-1);
    }

  printf(" ^^^^^^^^^^^^^ start computation ^^^^^^^^^^^^^\n");
  
  string analysisMethod("HARMONIC");
  cout << "use "<< analysisMethod << " analysis" << endl;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  default configuration initialization 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  size_t pos = -1;
  string rootname = input;
  
  if( (pos = rootname.find_last_of(".")) != -1){
    n = rootname.length() - pos;
    rootname.erase(pos,n);
    }
  string residualsfile = rootname + ".res.dat";

  if( mgrFormat=="") {
    mgrFormat = "ASCII";
    }

  std::vector<int> indexRecord = XTRACK_ref_get_indexRecords(input);
  int nRecords = indexRecord.size();
  printf("file=%s, nrecords=%d\n",input,nRecords);
  
  if(nRecords == 0) __REPORT_THIS_AND_EXIT__(log, "ERROR : XTRACK_ref_get_indexRecords() failed ");

  if(Trepet < 0.) __REPORT_THIS_AND_EXIT__(log, "ERROR : repeat period is not valid or non initialised ");
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  initialize optional de-tiding, de-loading and de-aliasing
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  spectrum_t OceanPredictionBaseList;
  spectrum_t LoadingPredictionBaseList;
  spectrum_t ogcmPredictionBaseList;
  hconstant_t *OceanTide_base   = 0;
  hconstant_t *LoadingTide_base = 0;
  hconstant_t *ogcmTide_base = 0;
  float atlas_fmask = 9999.f;

  if(( atlas_detiding == true ) or ( deloading == true ) or ( dealiasing == true ) or true) {
    double *x = new double[nRecords];
    double *y = new double[nRecords];
    int *nrecords = new int[nRecords];
    int status = -1;
    XTRACK_ref_get_AllHeaders(input, indexRecord, x, y, nrecords, status);
    if( atlas_detiding == true ) {
      cout << "use "<< atlasName <<" atlas for detiding" << endl;
      OceanPredictionBaseList=atlas_define(atlasName);
      OceanTide_base = tide_atlas2positions(oceanTide_convention, OceanPredictionBaseList, oceanTide_ampName, oceanTide_phaName, x, y, nRecords, atlas_fmask, verbose);
      if(OceanTide_base == 0) {
        __REPORT_THIS_AND_EXIT__(log, "ERROR : tide_atlas2positions() failed for ocean tide");
        }
      scale_constants(OceanTide_base,nRecords,OceanPredictionBaseList.n,scaleAtlas);
      }
    else printf("no prior de-tiding from external atlas done....\n");
    if( deloading == true ) {
      cout << "use "<< atlasName <<" atlas for deloading" << endl;
      LoadingPredictionBaseList=atlas_define(atlasName);
      LoadingTide_base = tide_atlas2positions(loading_convention, LoadingPredictionBaseList, loading_ampName, loading_phaName, x, y, nRecords, atlas_fmask, verbose);
      if(LoadingTide_base == 0) {
        __REPORT_THIS_AND_EXIT__(log, "ERROR : tide_atlas2positions() failed for load tide");
        }
      scale_constants(LoadingTide_base,nRecords,LoadingPredictionBaseList.n,scaleAtlas);
      }
    else printf("no prior de-loading from external atlas done....\n");
    if( dealiasing == true ) {
      cout << "use "<< "GLORYS-v2" <<" atlas for deloading" << endl;
//      if(ogcm_convention==0) ogcm_convention=strdup("WAVE-sossheig-atlas.nc");
      if(ogcm_convention==0) ogcm_convention=strdup("/home/softs/data/climatology/GLORYS-v2/WAVE.GLORYS.nc");
      if(ogcm_amplitude==0)  ogcm_amplitude=strdup("smoothed_a");
      if(ogcm_phaselag==0)   ogcm_phaselag =strdup("smoothed_G");
      ogcmPredictionBaseList=atlas_define(atlasName);
      ogcmTide_base = tide_atlas2positions(ogcm_convention, ogcmPredictionBaseList, ogcm_amplitude, ogcm_phaselag, x, y, nRecords, atlas_fmask, verbose);
      if(ogcmTide_base == 0) {
        __REPORT_THIS_AND_EXIT__(log, "ERROR : tide_atlas2positions() failed ");
        }
      }
    else printf("no prior de-aliasing from external atlas done....\n");
    delete[] x;
    delete[] y;
    delete[] nrecords;
    cout << endl;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  MAIN LOOP ON TIME SERIES
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  string outname;
  outname=rootname+".mgr.stat";
  foutname=fopen(outname.c_str(),"w");

  vector<mgr_t> mgr;

  FILE *residualsFile = 0;
  int header_ok = 0;
  int account_for_nodal_args = 1;
  
  int status = -1;
  vector<serie_t> tracks = ctoh_load_metadata(input, status);
  
  if(tracks.size()!=nRecords) TRAP_ERR_EXIT(ENOEXEC,"programming error %u!=%d\n",tracks.size(),nRecords);
  
  int nprocs=initialize_OPENMP(nRequestedProcs);
#ifdef OPEN_MP
  #pragma omp parallel for private(status) if(nprocs>1)
#endif
  
  for(size_t record = 0; record < nRecords; record++) {
    int proc_id=omp_get_thread_num();
    ostringstream oss;
    spectrum_t OceanPredictionExtendedList;
    spectrum_t LoadingPredictionExtendedList;
    serie_t &track=tracks[record];
    
    status=XTRACK_CheckRecords(track, indexRecord[record],rootname.c_str());
    
//    int nprocs=omp_get_num_threads();

/*------------------------------------------------------------------------------
    load time serie + corrections */
    if(track.count<=0){
      //STDOUT_BASE_LINE("tracks[%u].count=%d <= 0 : skipping\n",record,track.count);
      continue;
      }
    if(fabs(track.data[0].lat)<latmin) {
      for(int i = 0; i < track.count; i++){
        track.data[i].values[X_TRACK_SERIE_T_TIDE_ANALYSIS] = track.mask;
        track.data[i].values[X_TRACK_SERIE_T_RESIDUAL]      = track.mask;
        }
      continue;
      }
    
/*------------------------------------------------------------------------------
    reduce time frame */
    if(isad(reduced_start) or isad(reduced_final)) {
      track = aktarus_reduce(track, reduced_start, reduced_final);
      }
    
    if(track.count == 0) {
      __REPORT_THIS__(log, "WARNING: analysisof aborted... no data\n ");
      continue;
      }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    build corrected SLAs -> X_TRACK_SERIE_T_RESIDUAL
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    int use_gdr_tide;
    if(atlas_detiding) {
/*------------------------------------------------------------------------------
      de-tiding performed dynamically, no GDR ocean tide correction applied */
      use_gdr_tide = 0;
      }
    else if(gdr_detiding) {
/*------------------------------------------------------------------------------
      apply GDR ocean tide correction */
      use_gdr_tide = 1;
      }
    else {
/*------------------------------------------------------------------------------
      de-tiding performed dynamically, no GDR ocean tide correction applied */
      use_gdr_tide = 0;
      }
    
/*------------------------------------------------------------------------------
    process SLA*/
    XTRACK_ref_build_SLA(&track, use_gdr_tide, use_ib);
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    clean up time series and check number of valida data
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    double *z    = new double[track.count];
    double *time = new double[track.count];
    std::vector<int> maskIndex;
    size_t nValid = 0;

    for(size_t i = 0; i < track.count; i++){
      if( !is_equal(track.data[i].values[X_TRACK_SERIE_T_RESIDUAL], track.mask,float(1e-3)) ){
        maskIndex.push_back(i);
        z[nValid]    = track.data[i].values[X_TRACK_SERIE_T_RESIDUAL];
        time[nValid] = track.data[i].time;
        nValid++;
        }
      }

    if(nValid < nmesMin) {
#if VERBOSE > 2
      oss << "WARNING: analysis aborted... not enough data (" << nValid << " mes.)" << endl;
      __REPORT_THIS__(log, oss.str());
#endif
      if(verbose>=2) printf("%d/%d: record #%d %4d valid measures (min=%4d), skipped\n",proc_id,nprocs,indexRecord[record],nValid, nmesMin);
      for(int i = 0; i < track.count; i++){
        track.data[i].values[X_TRACK_SERIE_T_TIDE_ANALYSIS] = track.mask;
        track.data[i].values[X_TRACK_SERIE_T_RESIDUAL]      = track.mask;
        }
      continue;
      }

/*------------------------------------------------------------------------------
    set optimal spectrum */
    spectrum_t analysed;
    spectrum_reset(analysed, spectrum, Trepet, debuglog);

    double duration = time[nValid - 1] - time[0];
    int *keep    = new int[analysed.n];
    int *deduce  = new int[analysed.n];
    int *resolve = new int[analysed.n];
    for(size_t w = 0; w < analysed.n; ++w){
      keep[w]    = 1;
      resolve[w] = 1;
      deduce[w]  = 0;
      }

    if(auto_parse == false) {
/*------------------------------------------------------------------------------
      spectrum parsing */
/*----------------------------------------------------------------------------
      check aliased period against serie duration */
      spectrum_isSolvable(analysed, duration, keep, resolve, debuglog);
/*----------------------------------------------------------------------------
      check aliased period separation */
      spectrum_checkRayleigh(analysed, duration, keep,resolve, debuglog);
/*----------------------------------------------------------------------------
      check admittance salvation*/
      spectrum_isAdmittanceSolvable(analysed, duration, keep, resolve, deduce, debuglog);
/*----------------------------------------------------------------------------
      remove unsolved constituents */
      spectrum_reduce(analysed, keep, deduce, debuglog);
      }
    
    int nWaves = 0;
    for(size_t i = 0; i < analysed.n ; i++){
      if( (keep[i] == 1) && (deduce[i] == 0) ){
        nWaves++;
        }
      }

/* *---------------------------------------------------------------------------
    limit spectrum size as a function of the number of observations*/
    if (2 * nWaves > nValid) {
#if VERBOSE > 0
      oss << "WARNING: analysis aborted... too short or too many waves ("
          << track.count << " mes. and " << nWaves << " waves)" << endl << endl;
      __REPORT_THIS__(log, oss.str());
#endif
      for(int i = 0; i < track.count; i++){
        track.data[i].values[X_TRACK_SERIE_T_TIDE_ANALYSIS] = track.mask;
        track.data[i].values[X_TRACK_SERIE_T_RESIDUAL]      = track.mask;
        }
      printf(" skipped\n");
      continue;
      }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    preliminary de-tiding (optional) from user prescribed atlas
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    hconstant_t OceanTide_extended;
    if( atlas_detiding == true ) {
      if(OceanTide_base[record].hasMaskValue(atlas_fmask) == true){
#if VERBOSE > 0
        __REPORT_THIS__(log, "dediting failed: no prior model for some constituents\n\n");
#endif
        printf(" skipped\n");
        for(int i = 0; i < track.count; i++){
          track.data[i].values[X_TRACK_SERIE_T_TIDE_ANALYSIS] = track.mask;
          track.data[i].values[X_TRACK_SERIE_T_RESIDUAL]      = track.mask;
          }
        continue;
        }
      
      double x = track.data[0].lon;
      double y = track.data[0].lat;
      double sla;
      int anomaly=0;
/*------------------------------------------------------------------------------
      used to be multi-thread unsafe - fixed */
      OceanTide_extended = extended_atlas_constants(atlasName, OceanPredictionBaseList, &OceanPredictionExtendedList, x, y, OceanTide_base[record]);

      double *tides = new double[nValid];
      harmonic_predictionTS(tides, time, nValid, OceanTide_extended, account_for_nodal_args);
      
      size_t j;
      
      for(j = 0 ; j < nValid ; j++) {
        sla=z[j];
        z[j] -= tides[j];
        if(abs(z[j])> 1.e+03 or not isfinite(z[j])) {
          anomaly=1;
          printf("alert %d: t=%.6f, sla=%.6f, detided sla=%.6f, tides=%.6f %.6f  (FES2004) %.6f (GOT4.7)\n",record,time[j], sla, z[j], tides[j],
               track.data[j].values[X_TRACK_SERIE_T_TIDE_FES2004],
               track.data[j].values[X_TRACK_SERIE_T_TIDE_GOT4_7]);
          fflush(stdout);
          if(not isfinite(tides[j]) ){
            STDERR_BASE_LINE("tracks[%u] tide[%u]=%g prediction failed\n",record,j,tides[j]);
            break;
            }
          }
        }
      
      delete [] tides;
      
      if(j < nValid)
        continue;
      }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    preliminary deloading (optional) from user prescribed atlas
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    hconstant_t LoadingTide_extended;
    if( deloading == true ) {
      if(LoadingTide_base[record].hasMaskValue(atlas_fmask) == true){
#if VERBOSE > 0
        __REPORT_THIS__(log, "deloading failed: no value for some constituents\n\n");
#endif
        printf(" skipped\n");
        for(int i = 0; i < track.count; i++){
          track.data[i].values[X_TRACK_SERIE_T_TIDE_ANALYSIS] = track.mask;
          track.data[i].values[X_TRACK_SERIE_T_RESIDUAL]      = track.mask;
          }
        continue;
        }
      
      double x = track.data[0].lon;
      double y = track.data[0].lat;
/*------------------------------------------------------------------------------
      used to be multi-thread unsafe - fixed */
      LoadingTide_extended =  extended_atlas_constants(atlasName, LoadingPredictionBaseList, &LoadingPredictionExtendedList, x, y, LoadingTide_base[record]);

      double *loading = new double[nValid];
      harmonic_predictionTS(loading, time, nValid, LoadingPredictionExtendedList, LoadingTide_extended.a, LoadingTide_extended.G, account_for_nodal_args);

      for(size_t j = 0 ; j < nValid ; j++) {
        z[j] -= loading[j];
        }
      delete [] loading;
      }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    preliminary dealiasing (optional) from user prescribed atlas (non-tidal)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    if( dealiasing == true ) {
      if(ogcmTide_base[record].hasMaskValue(atlas_fmask) == true){
#if VERBOSE > 0
        __REPORT_THIS__(log, "de-aliasing failed: no value for some constituents\n\n");
#endif
        printf(" skipped\n");
        for(int i = 0; i < track.count; i++){
          track.data[i].values[X_TRACK_SERIE_T_TIDE_ANALYSIS] = track.mask;
          track.data[i].values[X_TRACK_SERIE_T_RESIDUAL]      = track.mask;
          }
        continue;
        }
      
//       double x = track.data[0].lon;
//       double y = track.data[0].lat;

      double *ogcm = new double[nValid];
      harmonic_predictionTS(ogcm, time, nValid, ogcmPredictionBaseList, ogcmTide_base[record].a, ogcmTide_base[record].G, account_for_nodal_args);

      for(size_t j = 0 ; j < nValid ; j++) {
        z[j] -= ogcm[j];
        }
      delete [] ogcm;
      }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    tidal harmonic analysis
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    mgr_data_t *data;
    double *residuals = new double[nValid];

    printf("%d/%d: record #%d %4d valid measures (min=%4d), analysed\n",proc_id,nprocs,indexRecord[record],nValid, nmesMin);
    date_t start = poctime_getdatecnes(time[0], 'd');
    date_t final = poctime_getdatecnes(time[nValid - 1], 'd');
    printf("%02d/%02d/%4d %02d/%02d/%4d\n",start.day,start.month,start.year, final.day,final.month,final.year);

    if( analysisMethod=="HARMONIC" ) {
      double mean=0;
//       data=check_parsing(time, nValid, analysed, account_for_nodal_args, debuglog);
      if(auto_parse == false) {
        data = harmonic_analysis_without_parsing(z, time, residuals, &mean, nValid, analysed, account_for_nodal_args,keep,deduce,debuglog);
        }
      else {
        data=harmonic_analysis(z, residuals, time, nValid, analysed, averaging_time, keep, deduce, OceanTide_extended, account_for_nodal_args, maxratio, debuglog, indexRecord[record]);
        }
      if( data == NULL) {
        __REPORT_THIS__(log, analysisMethod + " analysis failed\n\n");
        for(int i = 0; i < track.count; i++){
          track.data[i].values[X_TRACK_SERIE_T_TIDE_ANALYSIS] = track.mask;
          track.data[i].values[X_TRACK_SERIE_T_RESIDUAL]      = track.mask;
          }
        continue;
        }
#if defined(DEBUG)
/*------------------------------------------------------------------------------
      multi-thread unsafe, also should be moved inside extern routine */
//      status= check_timeseries(Trepet, time, z, residuals);
#endif
      }
    else {
      __REPORT_THIS_AND_EXIT__(log, "ERROR : analysis method not available");
      }
    
    stats_input     = get_statistics(z,          track.mask, nValid);
    stats_residuals = get_statistics(residuals,  track.mask, nValid);
    
    fprintf(foutname,"%d\t %lf\t %lf\t %d\t %8.4lf\t  %8.4lf\t  %8.4lf\t %8.4lf\t  %8.4lf\t  %8.4lf\t  %8.4lf\n",
      record, track.data[0].lon, track.data[0].lat, nValid,
      stats_input.mean, square(stats_input.std), stats_input.std,
      stats_residuals.mean, square(stats_residuals.std), stats_residuals.std,
      square(stats_input.std) - square(stats_residuals.std));
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    set error budget estimate
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    if(error_budget=="psd"){
      vector<double> bg_error = spectrum_psd_average(analysed, residuals, time, nValid, Trepet, false);
      for(size_t w = 0; w < analysed.n; w++){
        data[w].bg_contamination_error = bg_error[w];
        if(mgrFormat=="ASCII"){
          data[w].error = bg_error[w]; // overwrite formal inversion error in mgr ASCII output
          }
        }
      }
    
    delete[] keep;
    delete[] deduce;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    (possibly) restore prior constants
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    if(atlas_detiding == true){
      if(recompose) status=restore(analysed, data, OceanPredictionExtendedList, OceanTide_extended);
      OceanTide_extended.destroy();
      }
      
    if(deloading == true){
      LoadingTide_extended.destroy();
      LoadingPredictionExtendedList.destroy();
      }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    store constants in mgr array
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    int& initialIndex = indexRecord[record];
    status=register_analysis(initialIndex, mgr, data, input, track, Trepet, analysed, nValid, time);
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    write processed ASCII data, multi-thread unsafe
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//#pragma omp critical(threadUnsafe)
    if(nprocs==1) {
      if(!header_ok) {
        residualsFile = fopen((char *)residualsfile.c_str(),"w");
        if( residualsFile == 0 ){
          __REPORT_THIS_AND_EXIT__(log, "ERROR : can not create harmonic prediction file \n");
          }
        XTRACK_ref_write_header(residualsFile,nRecords);
        header_ok = 1;
        }
      int DAC = X_TRACK_SERIE_T_DAC_GLOBAL, TIDE = X_TRACK_SERIE_T_TIDE_FES2004;

      if(use_ib==0){
        DAC = XTRACK_get_meteo(&track);
        }
        
      if(use_gdr_tide){
        TIDE = XTRACK_get_tide(&track);
        }

      for(int i = 0, j = 0; i < maskIndex.size(); i++){
        j = maskIndex.at(i);
        track.data[j].values[X_TRACK_SERIE_T_TIDE_ANALYSIS] = track.data[j].values[X_TRACK_SERIE_T_RESIDUAL] - residuals[i];
        track.data[j].values[X_TRACK_SERIE_T_RESIDUAL]      = residuals[i];
        if(use_ib==0){
/*------------------------------------------------------------------------------
          remove GDR DAC correction to residual as it was not applied to SLA*/
          track.data[j].values[X_TRACK_SERIE_T_RESIDUAL] -= track.data[j].values[DAC];
          }
        if(use_gdr_tide) {
/*------------------------------------------------------------------------------
          add GDR tide correction to tide as it was applied to SLA */
          track.data[j].values[X_TRACK_SERIE_T_TIDE_ANALYSIS] += track.data[j].values[TIDE];
          }
        if(atlas_detiding) {
/*------------------------------------------------------------------------------
          TODO : add ATLAS tide correction to tide as it was applied to SLA */
//           track.data[j].values[X_TRACK_SERIE_T_TIDE_ANALYSIS] += track.data[j].values[TIDE];
          }
        }
/*------------------------------------------------------------------------------
      write residual time series */
      if(save_gnu == true) {
        XTRACK_WriteAsciiRecords(residualsFile, &track, &(indexRecord[record]), 1, rootname.c_str());
        }
      else{
        XTRACK_WriteAsciiRecords(residualsFile, &track, &(indexRecord[record]), 1, 0);
        }
      }

/*------------------------------------------------------------------------------
    clean loop variables and reset analysis list*/
    maskIndex.clear();
    delete [] residuals;
    residuals = 0;
    delete [] z;
    z = 0;
    delete [] time;
    time = 0;
    }

/*------------------------------------------------------------------------------
  complete harmonic predictions */
  if(residualsFile != 0){
    fclose(residualsFile);
    }

  fclose (foutname);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    save residuals
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  string ncout=input;
  pos=ncout.find(".ref");
  if(pos==string::npos) pos=ncout.find("_sla");
  if(pos!=string::npos) {
    ncout.insert(pos,".empiric");
    status=ctoh_save_metadata_NetCDF(ncout, input, tracks);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    save constants
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  string mgrFileName = rootname + ".mgr";
  if(mgrFormat=="NETCDF"){
    mgrFileName = mgrFileName + ".nc";
    }
  if(mgr.size()<1){
    printf("NOTHING TO SAVE TO "+mgrFileName+"\n");
    goto end;
    }
  printf("save constants to "+mgrFileName+"\n");
  
  status=mgr_save(mgrFileName.c_str() ,mgr, mgrFormat);
  if(status != 0) {
    __REPORT_THIS_AND_EXIT__(log, "ERROR : can not save constants \n");
    }

  indexRecord.clear();
  for(size_t m = 0; m < mgr.size(); m++) {
    mgr[m].clear();
    }

end:
  cout << " ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^" << endl;
  return(0);
}
