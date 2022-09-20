
/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Thierry Letellier  LEGOS, Toulouse, France (PhD)
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

#include <config.h>

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <list>
#include <vector>
#include <algorithm>

#include "tools-structures.h"

#include "spectrum.h"
#include "rutin.h"

#include <gsl/gsl_matrix.h>

#ifndef VERBOSE
#define VERBOSE 1
#endif


const char* spectrum_list_default[] = SPECTRUM_LIST_DEFAULT;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void spectrum_print(spectrum_t s, FILE *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double dph2dpd=24.;
  double T;
  
  if(out!=0) {
    fprintf(out, "Spectrum, aliased frequencies:\n");
    }
  
  for(int i = 0; i < s.n; i++){
    if(out!=0) {
      T=360./s.waves[i].aliased/24.;
      if(isinf(T)==FP_INFINITE) {
        fprintf(out, "%10s : %12.2lf (degree per days)  %8.1lf\n", s.waves[i].name,s.waves[i].aliased*dph2dpd, T);
        }
      else if(T>1000.*365.25) {
        fprintf(out, "%10s : %12.2lf (degree per days) ~%8.1lf\n", s.waves[i].name,s.waves[i].aliased*dph2dpd, 1./0.);
        }
      else if(T>100.*365.25) {
        fprintf(out, "%10s : %12.2lf (degree per days)  %8.1lf (centuries)\n", s.waves[i].name,s.waves[i].aliased*dph2dpd, T/365.25/100.);
        }
      else if(T>365.25) {
        fprintf(out, "%10s : %12.2lf (degree per days)  %8.1lf (years)\n", s.waves[i].name,s.waves[i].aliased*dph2dpd, T/365.25);
        }
      else {
        fprintf(out, "%10s : %12.2lf (degree per days)  %8.1lf (days)\n", s.waves[i].name,s.waves[i].aliased*dph2dpd, T);
        }
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

spectrum_t spectrum_init_ref_deep(bool Z0, bool LF)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  spectrum_t spectrum;
  spectrum.waves = new tidal_wave[75];

  size_t i = 0;
/* *----------------------------------------------------------------------------
  major semi-diurnal*/
  spectrum.waves[i++] = wM2;
  spectrum.waves[i++] = wN2;
  spectrum.waves[i++] = wK2;
  spectrum.waves[i++] = wS2;
  spectrum.waves[i++] = w2N2;

/* *----------------------------------------------------------------------------
  major diurnal*/
  spectrum.waves[i++] = wO1;
  spectrum.waves[i++] = wK1;
  spectrum.waves[i++] = wQ1;
  spectrum.waves[i++] = wP1;

/* *----------------------------------------------------------------------------
  major long period*/
  spectrum.waves[i++] = wMf;
  spectrum.waves[i++] = wMm;
  spectrum.waves[i++] = wMtm;

/* *----------------------------------------------------------------------------
  major non-linear*/
  spectrum.waves[i++] = wM4;
  spectrum.waves[i++] = wMS4;

/* *----------------------------------------------------------------------------
  second rank semi-diurnal*/
  spectrum.waves[i++] = wNu2;
  spectrum.waves[i++] = wMu2;
  spectrum.waves[i++] = wL2;
  spectrum.waves[i++] = wT2;
  spectrum.waves[i++] = wR2;

/* *----------------------------------------------------------------------------
  second rank long period*/
  spectrum.waves[i++] = wMSqm;
  spectrum.waves[i++] = wMSm;
//  spectrum.waves[i++] = wSsa;
  spectrum.waves[i++] = wSa;
  
  spectrum.n = spectrum.nmax = i;
  for (size_t k = 0; k < spectrum.n; k++) {
    if(spectrum.waves[k].omega == 0.) {
      spectrum.waves[k].init();
      spectrum.waves[k].aliased = spectrum.waves[k].omega;
      }
    }

  if(not LF)
  for(int k=0; k<spectrum.n;k++) {
    if(spectrum.waves[k].omega<10.0) {
      spectrum.remove(k);
      k--;
      }
    }

  return(spectrum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

spectrum_t spectrum_init_ref_reduced(void)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  spectrum_t spectrum;
  spectrum.waves = new tidal_wave[75];

  size_t i = 0;
/* *----------------------------------------------------------------------------
  major semi-diurnal*/
  spectrum.waves[i++] = wM2;
  spectrum.waves[i++] = wN2;
  spectrum.waves[i++] = wK2;
  spectrum.waves[i++] = wS2;
  spectrum.waves[i++] = w2N2;

/* *----------------------------------------------------------------------------
  major diurnal*/
  spectrum.waves[i++] = wO1;
  spectrum.waves[i++] = wK1;
  spectrum.waves[i++] = wQ1;

/* *----------------------------------------------------------------------------
  major long period*/
  spectrum.waves[i++] = wMf;

/* *----------------------------------------------------------------------------
  major non-linear*/
  spectrum.waves[i++] = wM4;
  
  spectrum.n = spectrum.nmax = i;
  for (size_t k = 0; k < spectrum.n; k++) {
    if(spectrum.waves[k].omega == 0.) {
      spectrum.waves[k].init();
      spectrum.waves[k].aliased = spectrum.waves[k].omega;
      }
    }

  return(spectrum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

spectrum_t spectrum_init_ref_shelf(bool Z0, bool LF)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  spectrum_t spectrum;
  spectrum.waves = new tidal_wave[75];

  size_t i = 0;
/* *----------------------------------------------------------------------------
  major semi-diurnal*/
  spectrum.waves[i++] = wM2;
  spectrum.waves[i++] = wN2;
  spectrum.waves[i++] = wK2;
  spectrum.waves[i++] = wS2;
  spectrum.waves[i++] = w2N2;

/* *----------------------------------------------------------------------------
  major diurnal*/
  spectrum.waves[i++] = wO1;
  spectrum.waves[i++] = wK1;
  spectrum.waves[i++] = wQ1;
  spectrum.waves[i++] = wP1;

/* *----------------------------------------------------------------------------
  major long period*/
  spectrum.waves[i++] = wMf;
  spectrum.waves[i++] = wMm;
  spectrum.waves[i++] = wMtm;

/* *----------------------------------------------------------------------------
  major non-linear*/
  spectrum.waves[i++] = wM4;
  spectrum.waves[i++] = wMS4;
  spectrum.waves[i++] = wM6;

/* *----------------------------------------------------------------------------
  second rank semi-diurnal*/
  spectrum.waves[i++] = wNu2;
  spectrum.waves[i++] = wMu2;
  spectrum.waves[i++] = wL2;
  spectrum.waves[i++] = wT2;

/* *----------------------------------------------------------------------------
  second rank long period*/
  spectrum.waves[i++] = wMSqm;
//  spectrum.waves[i++] = wSsa;
//  spectrum.waves[i++] = wSa;

  spectrum.waves[i++] = wMN4;
  spectrum.waves[i++] = wLa2;
  spectrum.waves[i++] = w2MS6;
  spectrum.waves[i++] = w2SM2;
  spectrum.waves[i++] = wMK4;
//  spectrum.waves[i++] = wMNS2;
//   spectrum.waves[i++] = w2MN6;
  spectrum.waves[i++] = wM_SK_2;
//   spectrum.waves[i++] = w3MS8;
//   spectrum.waves[i++] = w3MS4;
  spectrum.waves[i++] = wMSK2;
  spectrum.waves[i++] = wN4;
  spectrum.waves[i++] = wM_KS_2;
//   spectrum.waves[i++] = wS1;
  spectrum.waves[i++] = w2MK6;
  spectrum.waves[i++] = w2SM6;
  spectrum.waves[i++] = wE2;
//   spectrum.waves[i++] = wOQ2;
  spectrum.waves[i++] = wMKS2;
//  spectrum.waves[i++] = wSN4;
//   spectrum.waves[i++] = wMSN6;
  spectrum.waves[i++] = wS4;
//   spectrum.waves[i++] = wS3;
  spectrum.waves[i++] = wM3;
  spectrum.waves[i++] = wR2;
//   spectrum.waves[i++] = w2MK3;
//   spectrum.waves[i++] = wSK4;
  spectrum.waves[i++] = wMSf;
//   spectrum.waves[i++] = wMSK6;
//   spectrum.waves[i++] = wM1;
  spectrum.waves[i++] = wPsi1;
//   spectrum.waves[i++] = wSO3;
//   spectrum.waves[i++] = wRo1;
//   spectrum.waves[i++] = wMSN2;
//   spectrum.waves[i++] = wMK3;
//   spectrum.waves[i++] = wSK3;
//   spectrum.waves[i++] = wMSm;
//   spectrum.waves[i++] = wKJ2;
  spectrum.waves[i++] = wMStm;
//   spectrum.waves[i++] = wJ1;
  spectrum.waves[i++] = wSig1;
//   spectrum.waves[i++] = w2Q1;
//   spectrum.waves[i++] = wSO1;
//   spectrum.waves[i++] = wMP1;
//   spectrum.waves[i++] = wOO1;
//   spectrum.waves[i++] = wMqm;
  spectrum.waves[i++] = wKi1;
/*  spectrum.waves[i++] = w2MK2;*/
  spectrum.waves[i++] = wPhi1;
//   spectrum.waves[i++] = wMO3;
  spectrum.waves[i++] = wTta1;
//   spectrum.waves[i++] = wKQ1;
  spectrum.waves[i++] = wPi1;

  spectrum.n = spectrum.nmax = i;
  for (size_t k = 0; k < spectrum.n; k++) {
    if(spectrum.waves[k].omega == 0.) {
      spectrum.waves[k].init();
      spectrum.waves[k].aliased = spectrum.waves[k].omega;
      }
    }
  if(not LF)
  for(int k=0; k<spectrum.n;k++) {
    if(spectrum.waves[k].omega<10.0) {
      spectrum.remove(k);
      k--;
      }
    }
  

  return(spectrum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

spectrum_t spectrum_init_coastal(bool Z0, bool LF)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
 // reference spectrum based on Calais harmonic analysis
  spectrum_t spectrum;
  spectrum.nmax = 100;

  spectrum.waves = new tidal_wave[spectrum.nmax];

  size_t i = 0;
/* *----------------------------------------------------------------------------
  major semi-diurnal*/
  spectrum.waves[i++] = wM2;
  spectrum.waves[i++] = wN2;
  spectrum.waves[i++] = wK2;
  spectrum.waves[i++] = wS2;
  spectrum.waves[i++] = w2N2;

/* *----------------------------------------------------------------------------
  major diurnal*/
  spectrum.waves[i++] = wO1;
  spectrum.waves[i++] = wK1;
  spectrum.waves[i++] = wQ1;
  spectrum.waves[i++] = wP1;

/* *----------------------------------------------------------------------------
  major long period*/
  spectrum.waves[i++] = wMf;
  spectrum.waves[i++] = wMm;
  spectrum.waves[i++] = wMtm;

/* *----------------------------------------------------------------------------
  major non-linear*/
  spectrum.waves[i++] = wM4;
  spectrum.waves[i++] = wMS4;
  spectrum.waves[i++] = wM6;

/* *----------------------------------------------------------------------------
  second rank semi-diurnal*/
  spectrum.waves[i++] = wNu2;
  spectrum.waves[i++] = wMu2;
  spectrum.waves[i++] = wL2;
  spectrum.waves[i++] = wT2;

/* *----------------------------------------------------------------------------
  second rank long period*/
  spectrum.waves[i++] = wMSqm;
  spectrum.waves[i++] = wSsa;
  spectrum.waves[i++] = wSa;

  spectrum.waves[i++] = wMN4;
  spectrum.waves[i++] = wLa2;
  spectrum.waves[i++] = w2MS6;
  spectrum.waves[i++] = w2SM2;
  spectrum.waves[i++] = wMK4;
  spectrum.waves[i++] = wMNS2;
  spectrum.waves[i++] = w2MN6;
  spectrum.waves[i++] = wM_SK_2;
  spectrum.waves[i++] = w3MS8;
  spectrum.waves[i++] = w3MS4;
  spectrum.waves[i++] = wMSK2;
  spectrum.waves[i++] = wN4;
  spectrum.waves[i++] = wM_KS_2;
  spectrum.waves[i++] = wS1;
  spectrum.waves[i++] = w2MK6;
  spectrum.waves[i++] = w2SM6;
  spectrum.waves[i++] = wE2;
  spectrum.waves[i++] = wOQ2;
  spectrum.waves[i++] = wMKS2;
  spectrum.waves[i++] = wSN4;
  spectrum.waves[i++] = wMSN6;
  spectrum.waves[i++] = wS4;
  spectrum.waves[i++] = wS3;
  spectrum.waves[i++] = wM3;
  spectrum.waves[i++] = wR2;
  spectrum.waves[i++] = w2MK3;
  spectrum.waves[i++] = wSK4;
  spectrum.waves[i++] = wMSf;
  spectrum.waves[i++] = wMSK6;
  spectrum.waves[i++] = wM1;
  spectrum.waves[i++] = wPsi1;
  spectrum.waves[i++] = wSO3;
  spectrum.waves[i++] = wRo1;
  spectrum.waves[i++] = wMSN2;
  spectrum.waves[i++] = wMK3;
  spectrum.waves[i++] = wSK3;
  spectrum.waves[i++] = wMSm;
  spectrum.waves[i++] = wKJ2;
  spectrum.waves[i++] = wMStm;
  spectrum.waves[i++] = wJ1;
  spectrum.waves[i++] = wSig1;
  spectrum.waves[i++] = w2Q1;
  spectrum.waves[i++] = wSO1;
  spectrum.waves[i++] = wMP1;
  spectrum.waves[i++] = wOO1;
  spectrum.waves[i++] = wMqm;
  spectrum.waves[i++] = wKi1;
  spectrum.waves[i++] = w2MK2;
  spectrum.waves[i++] = wPhi1;
  spectrum.waves[i++] = wMO3;
  spectrum.waves[i++] = wTta1;
  spectrum.waves[i++] = wKQ1;
  spectrum.waves[i++] = wPi1;

  spectrum.n = i;
  for (size_t k = 0; k < spectrum.n; k++) {
    if(spectrum.waves[k].omega == 0.) {
      spectrum.waves[k].init();
      spectrum.waves[k].aliased = spectrum.waves[k].omega;
      }
    }
    
  if(not LF)
  for(int k=0; k<spectrum.n;k++) {
    if(spectrum.waves[k].omega<10.0) {
      spectrum.remove(k);
      k--;
      }
    }

  return(spectrum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

spectrum_t spectrum_init_estuarine(bool Z0, bool LF)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  spectrum_t spectrum;
  compound_tools_t compound;
//   int status;
  
  spectrum=initialize_tide(Z0);
  
  if(not LF)
  for(int k=0; k<spectrum.n;k++) {
    if(spectrum.waves[k].omega<10.0) {
      spectrum.remove(k);
      k--;
      }
    }
  
//   status=spectrum.add(wN8);
//   status=spectrum.add(wM8);
//   status=spectrum.add(wS8);
  
//   status=spectrum.add(wMSf);

//  status=spectrum.add(compound.init("NKM2").add(wN2,1).add(wK2,1).add(wM2,-1).realize());

  return(spectrum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

spectrum_t spectrum_init_ref(string refClass, bool Z0, bool exitOnError)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  if(refClass=="ESTUARINE"){
    bool LF=true;
    return(spectrum_init_estuarine(Z0, LF));
    }

  if(refClass=="ESTUARINE-HF"){
    bool LF=false;
    return(spectrum_init_estuarine(Z0, LF));
    }

  if(refClass=="COMODO"){
    return(initialize_tide());
    }

  if(refClass=="COASTAL"){
    bool LF=true;
    return spectrum_init_coastal(Z0, LF);
    }

  if(refClass=="COASTAL-HF"){
    bool LF=false;
    return spectrum_init_coastal(Z0, LF);
    }

  if(refClass=="SHELF"){
    bool LF=true;
    return(spectrum_init_ref_shelf(Z0, LF));
    }

  if(refClass=="SHELF-HF"){
    bool LF=false;
    return(spectrum_init_ref_shelf(Z0, LF));
    }

  if(refClass=="DEEP"){
    bool LF=true;
    return(spectrum_init_ref_deep(Z0, LF));
    }
    
  if(refClass=="DEEP-HF"){
    bool LF=false;
    return(spectrum_init_ref_deep(Z0, LF));
    }
    
  if(refClass=="REDUCED"){
    return(spectrum_init_ref_reduced());
    }
    
  return(spectrum_init_from_file(refClass,exitOnError));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  spectrum_t spectrum_init_from_file(const string & filename,bool exitOnError)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,code,i,nwaves,addnul=0;
  char name[64];
  double omega;
  FILE  *in=NULL;
  spectrum_t basic=initialize_tide();
  spectrum_t spectrum;
  
  in=fopen(filename.c_str(),"r");
  if(in==0) {
    if(exitOnError)
      TRAP_ERR_EXIT(-1,"spectrum file not found : "+filename+"\n");
    TRAP_ERR_RETURN(spectrum,1,"spectrum file not found : "+filename+"\n");
    }

  status=fscanf(in,"%d ",&nwaves);
  spectrum.init(nwaves);

  for(i=0;i<nwaves;i++) {
    status=fscanf(in,"%d %s %*c %lf",&code,name,&omega);
    if(code==999) {
      spectrum.waves[i]=wNUL;
      spectrum.waves[i].nT=addnul;
      addnul++;
      }
    else {
      status=spectrum.add (basic,name,0);
      if(status==-1) {
        TRAP_ERR_EXIT(-1,"cannot register wave : %s\n",name);
        }
      if(spectrum.n-1!=i) {
        STDERR_BASE_LINE("i=%d, status=%d, wave %s duplicated\n",i,status,name);
        }
      }
    // sara - what's for ? Need at least the feof test that was not there !
    // do {a=fgetc(in);} while ( a!= '\n'&& !feof(in));
    }
  fclose(in);

  for(i=0;i<spectrum.n;i++) {
    if(spectrum.waves[i].omega == 0.) {
      spectrum.waves[i].init();
      spectrum.waves[i].aliased = spectrum.waves[i].omega;
      }
#if VERBOSE > 1
      printf ("wave: %10s, period: %15.11f days \n",
              spectrum.waves[i].name,  tide_period(spectrum.waves[i]) / 24.);
#endif
    }
  return spectrum;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  spectrum_t spectrum_match(spectrum_t s, spectrum_t reference)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  spectrum_t matching;
  matching.waves = new tidal_wave[reference.n];
  matching.n = 0;
  for(size_t i = 0; i != s.n; i++){
    for(size_t j = 0; j != reference.n ; j++){
      if(strcmp(s.waves[i].name, reference.waves[j].name) == 0){
        memmove(&(matching.waves[matching.n++]), &(reference.waves[j]), sizeof(tidal_wave));
        break;
        }
      }
    }

 matching.nmax = matching.n;
 return(matching);

}


inline bool wcompare(tidal_wave first, tidal_wave second){
  return(first.omega < second.omega);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 void spectrum_order(spectrum_t *s, string orderCrit)

 /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

    list<tidal_wave> spectrum;
    for(size_t k = 0; k < s->n ; k++){
      spectrum.push_back(s->waves[k]);
      }

    if(orderCrit=="INVERSE"){
      spectrum.reverse();
      }
    if(orderCrit=="FREQUENCY"){
      spectrum.sort(wcompare);
      }

    size_t size = sizeof(tidal_wave);
    size_t k = 0;
    for(list<tidal_wave>::iterator it = spectrum.begin(); it != spectrum.end(); it++, k++){
      memmove(&(s->waves[k]), &(it), size);
      }

}

 /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 void spectrum_order(spectrum_t *s, string orderCrit, std::vector<size_t> &perm)

 /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

    list<tidal_wave> spectrum;
    for(size_t k = 0; k < s->n ; k++){
      spectrum.push_back(s->waves[k]);
      }
    list<tidal_wave> bck = spectrum;

    if(orderCrit=="REVERSE"){
      spectrum.reverse();
      }
    if(orderCrit=="FREQUENCY"){
      spectrum.sort(wcompare);
      }

    size_t size = sizeof(tidal_wave);
    size_t k = 0;
    for(list<tidal_wave>::iterator it = spectrum.begin(); it != spectrum.end(); it++){
      memmove(&(s->waves[k++]), &(*it), size);
      }

    for(list<tidal_wave>::iterator it1 = spectrum.begin(); it1 != spectrum.end(); it1++){
      k = 0;
      for(list<tidal_wave>::iterator it2 = bck.begin(); it2 != bck.end(); it2++, k++){
        if(strcmp((*it1).name, (*it2).name) == 0){
          perm.push_back(k);
          }
        }
      }

#if defined(DEBUG)
    for(size_t r = 0; r != perm.size(); r++){
      printf("%10s: %2d <- %02d\n", s->waves[r].name, r, perm[r]);
      }
#endif

    spectrum.clear();
    bck.clear();

  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void spectrum_reduce(spectrum_t &s, int *keep, int *deduce, FILE *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  spectrum_t tmp;
  tmp.waves = new tidal_wave[s.n];
  int *LocalKeep    = new int[s.n];
  int *LocalDeduce  = new int[s.n];
  
  for(size_t w = 0; w < s.n; w++){
    LocalKeep[w]  =keep[w];
    LocalDeduce[w]=deduce[w];
    }
 
  for(size_t w = 0; w < s.n; w++){
    if(LocalKeep[w] == 1){
      tmp.waves[tmp.n] = s.waves[w];
      keep[tmp.n]=LocalKeep[w];
      deduce[tmp.n]=LocalDeduce[w];
      tmp.n++;
      }
    }
  
  s = tmp;
  spectrum_print(s, out);
  
  delete[] LocalKeep;
  delete[] LocalDeduce;
}


 /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 gsl_matrix *adjacencyMatrix_build(spectrum_t s, double duration)

 /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

  build adjacency matrix:

  matrix(i,j] = 1 if wave i and j not separable
  matrix(i,j] = 0 if wave i and j separable

  separation : at least 1 cycle difference over duration

  tau such that :

    abs(fi -fj) tau = 2 PI = abs(2PI/Ti -2PI/Tj) tau

    abs(Ti-Tj)/Ti*Tj = 1/tau

-----------------------------------------------------------------------------*/
 {
  double tau = 0;
  double dph2d = static_cast<double> (360. / 24.);
  double aliasT_i = 0, aliasT_j = 0;

  gsl_matrix *adjacencyMatrix = gsl_matrix_alloc(s.n,s.n);

  gsl_matrix_set_all(adjacencyMatrix, 0.);

  for(size_t i = 0; i < s.n; i++) {
    for (size_t j = 0 ; j < s.n; j++){
      if(j == i) continue;
      aliasT_i = dph2d / s.waves[i].aliased ; // aliased period in days
      aliasT_j = dph2d / s.waves[j].aliased ;
      tau =  fabs((aliasT_i * aliasT_j) / (aliasT_i - aliasT_j)); // separation time for aliased periods

      if(duration < tau){
        gsl_matrix_set(adjacencyMatrix,i,j,1);
#if VERBOSE > 2
        cout << s.waves[i].name <<" and "<< s.waves[j].name << endl;
#endif
        }
      else {
        gsl_matrix_set(adjacencyMatrix,i,j,0);
        }
      }
    }

 return(adjacencyMatrix);
 }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 std::vector<EnvirWave> adjacencyList_build(spectrum_t s, double duration)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
 {
  double *vertexWeight1 = 0, *vertexWeight2 = 0 ;
  gsl_matrix *edgeWeight = 0;

/* *----------------------------------------------------------------------------
  separation */
  edge_setWeight_byRaileigh(s,duration,&edgeWeight);
/* *----------------------------------------------------------------------------
  base wave for adimmtance, with potential, without potential etc... */
  vertex_setWeight_asIndividual(s,vertexWeight1);
/* *----------------------------------------------------------------------------
  Magnitude criteria based on spectrum order? */
  vertex_setWeight_byMagnitude(s,vertexWeight2);

  std::vector<EnvirWave> separations (s.n);

  for(size_t i = 0; i < s.n; i++) {
    const double w =  vertexWeight2[i];
    const AdjacencyRule tmpPair(s.waves[i], 1.);
/* *----------------------------------------------------------------------------
    initialize ith element of separation list with wave i */
    separations[i].first.push_back(tmpPair);
    separations[i].second = w;

    for (size_t j = 0 ; j < s.n; j++){
      if(j == i) continue;
      const double dph2d = static_cast<double> (360. / 24.);
      const double aliasT_i = dph2d / s.waves[i].aliased ; // aliased period in days
      const double aliasT_j = dph2d / s.waves[j].aliased ;
      const double tau =  fabs((aliasT_i * aliasT_j) / (aliasT_i - aliasT_j)); // separation time for aliased periods
      if(duration < tau){
        const double wij = gsl_matrix_get(edgeWeight,i,j) * vertexWeight1[j] * vertexWeight2[j];
/* *----------------------------------------------------------------------------
        add wave j to ith speparation list*/
        const AdjacencyRule tmpPair(s.waves[j], wij);
        separations[i].first.push_back(tmpPair);
        }
      }// wave j
    }// wave i

  gsl_matrix_free(edgeWeight);
  edgeWeight = 0;
  delete [] vertexWeight1;
  vertexWeight1 = 0;
  delete [] vertexWeight2;
  vertexWeight2 = 0;

  return(separations);
  }

inline bool operator < (const AdjacencyRule &w1, const AdjacencyRule &w2) {
  return(w1.second < w2.second);
  }
inline bool operator > (const AdjacencyRule &w1, const AdjacencyRule &w2) {
  return(w1.second > w2.second);
  }
inline bool operator == (const AdjacencyRule &w1, const AdjacencyRule &w2) {
  return(w1.first.code == w2.first.code);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 std::vector<tidal_wave> adjacencyList_reduce(std::vector<EnvirWave> &source, int nmes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------
   remove adjacent wave *j of k-th constituent
-----------------------------------------------------------------------------*/
{
  std::vector<tidal_wave> admittance;

  for(size_t k = 0; k < source.size(); k++){
     while(source[k].first.size() > 1){
       AdjacencyRuleList::iterator j = source[k].first.begin();
       ++j; // skip 1st sample, aka itself

/* *----------------------------------------------------------------------------
       remove in wave j from waves adjacency list (starting from wave k+1)*/
       for(size_t w = k + 1; w < source.size(); w++){
         AdjacencyRuleList::iterator start = source[w].first.begin();
         AdjacencyRuleList::iterator p = find(++start,source[w].first.end(),(*j));
         if(p != source[w].first.end()){
           source[k].first.erase(p);
           }
         }
/* *----------------------------------------------------------------------------
       move *j from source spectrum to admittance list*/
       for(size_t i = k + 1; i < source.size(); i++){
         const int code = j->first.code;
         if ((source[i].first.begin())->first.code == code){
/* *----------------------------------------------------------------------------
           move to admittance */
           admittance.push_back((source[i].first.begin())->first);
#if VERBOSE > 2
           printf("%s pushed to admittance to accommodate %s %d %d\n",(source[i].first.begin())->first.name,(source[k].first.begin())->first.name,i,k);
#endif	
/* *----------------------------------------------------------------------------
           remove from source */
           source.erase(source.begin() + i);
           break;
           }
         }
/* *----------------------------------------------------------------------------
       remove in wave j from wave k's adjacency list*/
       source[k].first.erase(j);
       }
     }

/* *----------------------------------------------------------------------------
   reduce list again to respect wave/observation cardinal ratio (2)*/
//   while(2 * source.size() > nmes){
   while(6 * source.size() > nmes){
     admittance.push_back(((source[source.size()-1]).first).begin()->first);
     source.erase(source.end());
     }

   return(admittance);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 void adjacencyMatrix_plot(gsl_matrix *adjacencyMatrix)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
 {
 FILE *testFile = NULL;
 testFile = fopen("adjacencyMatrix.dat","w");
 if(testFile==0) {
     TRAP_ERR_EXIT(-1,"Probleme a l ouverture du fichier : adjacencyMatrix.dat \n");
   }

 for(size_t i = 0; i < adjacencyMatrix->size1; i++) {
   for (size_t j = 0 ; j < adjacencyMatrix->size2; j++){
     if(gsl_matrix_get(adjacencyMatrix,i,j) == 1.){
       fprintf(testFile, "* ");
       }
     if(gsl_matrix_get(adjacencyMatrix,i,j) == 0.){
       fprintf(testFile, "  ");
       }
     }
   fprintf(testFile, "\n");
   }
 fclose(testFile);

 }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 void adjacencyMatrix_weight(gsl_matrix *adjacencyMatrix, spectrum_t s, double duration)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
 {
   double *vertexWeight1 = 0, *vertexWeight2 = 0 ;
   gsl_matrix *edgeWeight = 0;
   gsl_matrix *modified = gsl_matrix_alloc(s.n,s.n);
   gsl_matrix_set_zero(modified);

   edge_setWeight_byRaileigh(s,duration,&edgeWeight);
   vertex_setWeight_asIndividual(s,vertexWeight1);
   vertex_setWeight_byMagnitude(s,vertexWeight2);

   for(size_t i = 0; i < adjacencyMatrix->size1; i++) {
     double weightedContrib = 0.;
     size_t j;
     for (j = 0 ; j < adjacencyMatrix->size2; j++){
       weightedContrib += (gsl_matrix_get(adjacencyMatrix,i,j) * vertexWeight1[j] * vertexWeight2[j]);
       }
     gsl_matrix_set(modified,i,j,weightedContrib);
     }

   gsl_matrix_memcpy(adjacencyMatrix,modified);
   gsl_matrix_free(modified);

 } // end adjacencyMatrix_weight()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 void edge_setWeight_byRaileigh(spectrum_t s, double duration, gsl_matrix **weight)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

  build adjacency matrix:

  matrix(i,j] = 1 if wave i and j not separable
  matrix(i,j] = 0 if wave i and j separable

  separation : at least 1 cycle difference over duration

  tau such that :

    abs(fi -fj) tau = 2 PI = abs(2PI/Ti -2PI/Tj) tau

    abs(Ti-Tj)/Ti*Tj = 1/tau

-----------------------------------------------------------------------------*/
 {
  double tau = 0;
  double dph2d = static_cast<double> (360. / 24.);
  double aliasT_i = 0, aliasT_j = 0;

  *weight = gsl_matrix_alloc(s.n,s.n);
  gsl_matrix_set_identity(*weight);

  for(size_t i = 0; i < s.n; i++) {
    aliasT_i = dph2d / s.waves[i].aliased ; // aliased period in days
    for (size_t j = 0 ; j < s.n; j++){
      if(j == i) {
        if(s.waves[i].aliased==0) {
          gsl_matrix_set(*weight,i,i,1.);
          }
        else {
          continue;
          }
        }
      aliasT_j = dph2d / s.waves[j].aliased ;
      tau =  fabs((aliasT_i * aliasT_j) / (aliasT_i - aliasT_j)); // separation time for aliased periods
      if(duration < tau){
        gsl_matrix_set(*weight,i,j,1.);
        }
      else{
        gsl_matrix_set(*weight,i,j,0.);
        }
      }
    }

 }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 void vertex_setWeight_byMagnitude(spectrum_t s, double* &weight)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(weight == 0){
    weight = new double[s.n];
    }
  for(size_t i = 0; i < s.n; i++) {
    weight[i] = 1./(0.5 * (i + 1.));
    }

}


 /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 void vertex_setWeight_gotPotential(spectrum_t s, double* &weight)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  string wName = "";
  size_t m;

  wName = "2Q1";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }
  wName = "Rho1";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }
  wName = "Sig1";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }
  wName = "J1";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }
  wName = "Ki1";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }
  wName = "Tta1";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }
  wName = "Pi1";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }
  wName = "P1";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }
  wName = "Psi1";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }
  wName = "Phi1";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }
  wName = "OO1";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }
  wName = "MP1";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }
  wName = "KQ1";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }
  wName = "M1";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.001;
    }
  wName = "E2";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }
  wName = "2N2";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.01;
    }
  wName = "Mu2";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.01;
    }
  wName = "Nu2";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.01;
    }
  wName = "L2";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.01;
    }
  wName = "La2";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.01;
    }
  wName = "T2";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }
  wName = "S2";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.01;
    }
  wName = "KJ2";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }
  wName = "R2";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }

 } // end vertex_setWeight_gotPotential()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 void vertex_setWeight_gotNoPotential(spectrum_t s, double* &weight)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
 {
  string wName = "";
  size_t m;

    // S1
  wName = "S1";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }

 // compounds tides (no potential amplitude)
 wName = "SO1";
 if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
   weight[m] = 0.0001;
   }
 wName = "OQ2";
 if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
   weight[m] = 0.001;
   }
 wName = "2MK2";
 if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
   weight[m] = 0.0001;
   }
 wName = "MSK2";
 if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
   weight[m] = 0.01;
   }
 wName = "MSN2";
 if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
   weight[m] = 0.0001;
   }
 wName = "2SM2";
 if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
   weight[m] = 0.01;
   }
 wName = "M(SK)2";
 if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
   weight[m] = 0.01;
   }
 wName = "M(KS)2";
 if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
   weight[m] = 0.01;
   }
 wName = "MKS2";
 if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
   weight[m] = 0.01;
   }

 } // end vertex_setWeight_gotNoPotential()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 void vertex_setWeight_isBasewave(spectrum_t s, double* &weight)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
 {
   string wName = "";
   size_t m;

  wName = "M2";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 1.;
    }
  wName = "N2";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 1.;
    }
  wName = "K2";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.95;
    }

  wName = "K1";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.75;
    }
  wName = "O1";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.75;
    }
  wName = "Q1";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.7;
    }

  wName = "Mf";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.25;
    }
  wName = "Mm";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.25;
    }
  wName = "Mtm";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.25;
    }

 }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 void vertex_setWeight_asIndividual(spectrum_t s, double* &weight)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  weight = new double[s.n];

  for(size_t i = 0; i < s.n; i++) { // init weights
    weight[i] = 1.;
    }

   // set weights by (independent) categories
  vertex_setWeight_isBasewave(s,weight);
  vertex_setWeight_gotPotential(s,weight);
  vertex_setWeight_gotNoPotential(s,weight);

   // additional constituents
  string wName = "";
  size_t m;

   // long period constituents
  wName = "Sa";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }
  wName = "Ssa";
  if( (m = index_from_name(s,(char *) wName.c_str())) != -1) {
    weight[m] = 0.0001;
    }

 }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 spectrum_t convert2spectrum_t(std::vector<tidal_wave> list)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  spectrum_t s;

  s.n = s.nmax = list.size();
  s.waves = new tidal_wave[s.n];

  for(size_t k = 0; k< list.size(); k++){
    s.waves[k] = list[k];
    }

  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 spectrum_t convert2spectrum_t(std::vector<EnvirWave> list)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  spectrum_t s;

  s.n = s.nmax = list.size();
  s.waves = new tidal_wave[s.n];

  for(size_t k = 0; k< list.size(); k++){
    s.waves[k] = (list[k].first.begin())->first;
    }

  return(s);
}

#if 0
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 void spectrum_check_bugged(spectrum_t *optimal, double duration, int nmes, int **keep, int **deduce, FILE *out)

 /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  if(out!=0) {
    fprintf(out,"PRIOR spectrum: \n");
    for(size_t i = 0; i < optimal->n; i++){
      fprintf(out,"%s \n",optimal->waves[i].name);
      }
    }

/* *----------------------------------------------------------------------------
  build adjacency matrix for plotting */
  gsl_matrix *adjacencyMatrix = adjacencyMatrix_build(*optimal, duration);
  adjacencyMatrix_plot(adjacencyMatrix);
  gsl_matrix_free(adjacencyMatrix);

/* *----------------------------------------------------------------------------
  build list of adjacencies for each waves */
  std::vector<EnvirWave> directList = adjacencyList_build(*optimal,duration);

/* *----------------------------------------------------------------------------
  reduce adjacencies */
  std::vector<tidal_wave> admittanceList = adjacencyList_reduce(directList,nmes);

/* *----------------------------------------------------------------------------
  set admittance extension list */
  int *keepLocal    = new int[admittanceList.size()];
  int *deduceLocal  = new int[admittanceList.size()];
  int *resolveLocal = new int[admittanceList.size()];
  for(size_t i = 0; i < admittanceList.size(); i++){
    deduceLocal[i]  = 0;
    keepLocal[i]    = 1;
    resolveLocal[i] = 1;
    }

  spectrum_t directSpectrum     = convert2spectrum_t(directList);
  spectrum_t admittanceSpectrum = convert2spectrum_t(admittanceList);
  if(out!=0) {
    fprintf(out,"Direct spectrum: \n");
    for(size_t i = 0; i < directSpectrum.n; i++){
      fprintf(out,"%s \n",directSpectrum.waves[i].name);
      }
    }
  if(out!=0) {
    fprintf(out,"Admittance spectrum: \n");
    for(size_t i = 0; i < admittanceSpectrum.n; i++){
      fprintf(out,"%s \n",admittanceSpectrum.waves[i].name);
      }
    }

/* *----------------------------------------------------------------------------
  init admittance */
  int *method = admittance_check(directSpectrum,0);
  for(size_t k = 0; k < admittanceSpectrum.n; k++){
    if(admittance_mts_verify(admittanceSpectrum.waves[k],method) == 0){
      deduceLocal[k] = 1;
      keepLocal[k]   = 1;
      }
    }
  admittance_terminate();


/* *----------------------------------------------------------------------------
  check admittance extension list */
  spectrum_isAdmittanceSolvable(admittanceSpectrum,duration,keepLocal,resolveLocal,deduceLocal);

  // set optimal spectrum := direct + admittance spectra
  // set keep and deduce arrays (compatibility)
  optimal->n = optimal->nmax = (directSpectrum.n + admittanceSpectrum.n);
  int *ukeep = new int[optimal->n];
  int *udeduce = new int[optimal->n];

  for(size_t i = 0; i < directList.size(); i++){ // direct mode
    memmove(&(optimal->waves[i]), &(directSpectrum.waves[i]), sizeof(tidal_wave));
    ukeep[i] = 1;
    udeduce[i] = 0;
#if VERBOSE > 2
    cout << directSpectrum.waves[i].name << " solved in direct mode" << endl;
#endif
    }
  size_t n = directSpectrum.n;
  for(size_t i = 0; i < admittanceSpectrum.n; i++){ // admittance list (try)
    if(keepLocal[i]){
      ukeep[n] = 1;
      udeduce[n] = 1;
      if(admittance_mts_verify(admittanceSpectrum.waves[i],method) == 0){
        memmove(&(optimal->waves[n++]), &(admittanceSpectrum.waves[i]), sizeof(tidal_wave));
#if VERBOSE > 2
        cout << admittanceSpectrum.waves[i].name << " added in admittance list" << endl;
#endif
        }
#if VERBOSE > 2
      else{
        cout << "can not add " << admittanceSpectrum.waves[i].name << " in admittance list" << endl;
        ukeep[n] = 0;
        udeduce[n] = 0;
        }
#endif
      }
    else{ // trash
      ukeep[n] = 0;
      udeduce[n] = 0;
#if VERBOSE > 2
      cout << admittanceSpectrum.waves[i].name << " discarded" << endl;
#endif
      }
    }
 
  delete [] keepLocal;
  delete [] deduceLocal;

/**----------------------------------------------------------------------------
  reconstruct optimal spectrum */
  optimal->n = optimal->nmax = n;

  std::vector<size_t> perm;
  //spectrum_order(optimal, "FREQUENCY", perm); /// PERMUTATION
  *keep = new int[optimal->n];
  *deduce = new int[optimal->n];
  for(size_t i = 0; i < optimal->n; i++){
    perm.push_back(i); /// NO PERMUTATION
    (*keep)[i]   = ukeep[perm[i]];
    (*deduce)[i] = udeduce[perm[i]];
    }

  perm.clear();
  delete[] ukeep;
  delete[] udeduce;

  if(out!=0) {
    fprintf(out,"OPTIMAL spectrum: \n");
    for(size_t i = 0; i < optimal->n; i++){
      fprintf(out,"%s keep=%d deduced=%d\n",optimal->waves[i].name, (*keep)[i], (*deduce)[i]);
      }
    }

  delete [] admittanceSpectrum.waves;
  delete [] directSpectrum.waves;
}

#endif /* 0 */



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void init_spectrum(spectrum_t *spectrum, int *n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for spectrum_t::init(int)
/*----------------------------------------------------------------------------*/
{
  /* As spectrum_t::NULLPointers() is private, we MUST call it from
  a friend and, therefore, C++ function */
  spectrum->NULLPointers();
  
  spectrum->init(*n);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

extern "C" void init_spectrum_(spectrum_t *spectrum, int *n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Fortran wrapper for init_spectrum()
/*----------------------------------------------------------------------------*/
{
  init_spectrum(spectrum,n);
}


/*----------------------------------------------------------------------------*/
/// global reference spectrum, for Fortran wrappers
/**
This variable will be initialised with init_fullspectrum()
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

spectrum_t *fullspectrum=NULL;
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void init_fullspectrum()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// NOT FOR FORTRAN: called by Fortran wrappers to initialise ::fullspectrum
/*----------------------------------------------------------------------------*/
{
  if(fullspectrum==NULL){/// <h1>If needed</h1>
    ///it allocates fullspectrum
    fullspectrum=new spectrum_t;
    ///then it initialises it with initialize_tide()
    *fullspectrum=initialize_tide();
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void cleanfullspectrum()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Fortran wrapper for spectrum_t::clear()
/**
To be called when the ::fullspectrum is redundant, eventually,
to clean up the allocation carried-out by init_fullspectrum()
*/
/*----------------------------------------------------------------------------*/
{
  if(fullspectrum!=NULL){
    fullspectrum->destroy();
    delete fullspectrum;
    fullspectrum=NULL;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 extern "C" void add_wave_to_spectrum_(spectrum_t *spectrum,const char *wavename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Fortran wrapper for spectrum_t::add(spectrum_t,const char *)
/**
\param wavename SPACE-TERMINATED STRING BECAUSE OF C TO FORTRAN MESS !
\note actually Fortran strings are space padded anyway :)
*/
/*----------------------------------------------------------------------------*/
{
  ///sort out Fortran string \c wavename with poc_fortran_strdup()
  char *wavename_=poc_fortran_strdup(wavename);
  int n;
  
  ///It calls init_fullspectrum()
  init_fullspectrum();
  
  n=spectrum->add(*fullspectrum,wavename_,0);
  if(n<0)STDERR_BASE_LINE("WARNING : did not find wave: %s\n",wavename_);
  
  delete[]wavename_;
}
