
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

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

#define __TIDE_DEF_MAIN
#include "tides.def"

#include "spectrum.h"
#include "zapper.h"
#include "maths.h"

#ifndef VERBOSE
#define VERBOSE -1
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tide_addwave(spectrum_t *list, tidal_wave wave)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n = list->n, nwmax = list->nmax;
  size_t size = sizeof(tidal_wave);

  if(n >= nwmax) {
    nwmax += 10;
    tidal_wave *tmpWave = new tidal_wave[nwmax];
    for (size_t j = 0; j < list->n; j++) {
      tmpWave[j] = list->waves[j];
      }
    delete[] list->waves;
    list->waves = tmpWave;
    list->nmax = nwmax;
    }

  memmove(&(list->waves[n]), &wave, size);
  list->n = n+1;

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int harmonic_coefficients(int nframe, date_t start, double *time, spectrum_t s, harmonic_t *x, int nodal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, m, status, nw;
  tidal_wave wave;
  double f,u,V, V0, omega, mask=-1;

  astro_angles_t astro_angles;
  init_argument(&astro_angles,start);
  //  init_argument(start);

  nw = s.n;

  x->neq = 2 * s.n;
  x->nframes = nframe;

  x->cs = new double *[s.n];
  x->sn = new double *[s.n];

  for(k = 0; k < s.n; k++)
    x->cs[k] = new double[nframe];
  for(k = 0; k < s.n; k++)
    x->sn[k] = new double[nframe];

  x->A = new double[x->neq * x->neq];

  for(k = 0; k < x->neq*x->neq; k++)
    x->A[k] = 0;

  for(k=0; k<nw; k++) {
    wave=s.waves[k];
    omega=wave.omega*dph2rps;
    V0=greenwhich_argument(astro_angles,wave);
    if(nodal == 1) {
      u = nodal_phase(astro_angles,wave);
      f = nodal_factor(astro_angles,wave.formula);
      }
    else {
      u = 0;
      f = 1.;
      }
    for(m=0; m<nframe; m++) {
    if(time[m]!=mask) {
      V=(double) omega*time[m]+V0;
      x->cs[k][m] = f * cos(V + u);
      x->sn[k][m] = f * sin(V + u);
      }
    }
  }

  return(status);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void tides_compound()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  compound_tools_t ST1,ST2,ST3,SNK2;
  
  ST1.init("ST1");
  ST1.add(wN2,2);
  ST1.add(wK2,1);
  ST1.add(wS2,-2);
  
  ST2.init("ST2");
  ST2.add(wM2, 1);
  ST2.add(wN2, 1);
  ST2.add(wK2, 1);
  ST2.add(wS2,-2);
  
  ST3.init("ST3");
  ST3.add(wM2, 2);
  ST3.add(wS2, 1);
  ST3.add(wK2,-2);
  
  SNK2.init("SNK2");
  SNK2.add(wS2, 1);
  SNK2.add(wN2, 1);
  SNK2.add(wK2,-1);
  
  wST1=ST1.realize();
  wST2=ST2.realize();
  wST3=ST3.realize();
  wSNK2=SNK2.realize();
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  spectrum_t initialize_tide(astro_angles_t *astro_angles, date_t reference)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///wrapper for init_argument() and initialize_tide()
/*----------------------------------------------------------------------------*/
{
  spectrum_t s;
  
  init_argument(astro_angles,reference);
  
  s=initialize_tide();
  
  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  spectrum_t initialize_tide(bool Z0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Gives a spectrum containing all waves registered in TUGO
/*----------------------------------------------------------------------------*/
{
  int i;
  spectrum_t spectrum;
  compound_tools_t compound;

  spectrum.init(200/*+122*/);

  i=0;
  if(Z0){
    spectrum.waves[i++]= wZ0;
    spectrum.waves[i++]= wM0;
    }
  
/*------------------------------------------------------------------------------
  declare first admittance supporting waves to ease harmonic analysis parsing */
  spectrum.waves[i++]= wN2;
  spectrum.waves[i++]= wM2;
  spectrum.waves[i++]= wK2;

  spectrum.waves[i++]= wQ1;
  spectrum.waves[i++]= wO1;
  spectrum.waves[i++]= wK1;

  spectrum.waves[i++]= wMm;
  spectrum.waves[i++]= wMm_1;
  spectrum.waves[i++]= wMm_2;
  spectrum.waves[i++]= wMf;
  spectrum.waves[i++]= wMf_1;
  spectrum.waves[i++]= wMf_2;
  spectrum.waves[i++]= wMtm;
  
  
  spectrum.waves[i++]= wSa;
  spectrum.waves[i++]= wSa_1;
  spectrum.waves[i++]= wSsa;
  spectrum.waves[i++]= wSsa_1;
  spectrum.waves[i++]= wSta;
  spectrum.waves[i++]= wMSm;
//   spectrum.waves[i++]= wMm;
  spectrum.waves[i++]= wMSf;
//   spectrum.waves[i++]= wMf;
  spectrum.waves[i++]= wMStm;
  spectrum.waves[i++]= wMStm_1;
//   spectrum.waves[i++]= wMtm;
  spectrum.waves[i++]= wMqm;
  spectrum.waves[i++]= wMSqm;

//   spectrum.waves[i++]= wRH5;

  spectrum.waves[i++]= w2Q1;
  spectrum.waves[i++]= wSig1;
//   spectrum.waves[i++]= wQ1;
  spectrum.waves[i++]= wRo1;
//   spectrum.waves[i++]= wO1;
  spectrum.waves[i++]= wM1;
#if USE_M1_12
  spectrum.waves[i++]= wM1_1;
  spectrum.waves[i++]= wM1_2;
#endif
  spectrum.waves[i++]= wM1_3;
  spectrum.waves[i++]= wKi1;
  spectrum.waves[i++]= wPi1;
  spectrum.waves[i++]= wP1;
//   spectrum.waves[i++]= wK1;
  spectrum.waves[i++]= wPsi1;
  spectrum.waves[i++]= wPhi1;
  spectrum.waves[i++]= wS1;
  spectrum.waves[i++]= wTta1;
  spectrum.waves[i++]= wJ1;
  spectrum.waves[i++]= wOO1;
  spectrum.waves[i++]= wMP1;
  spectrum.waves[i++]= wSO1;
  spectrum.waves[i++]= wKQ1;

  spectrum.waves[i++]= wE2;
  spectrum.waves[i++]= w2N2;
  spectrum.waves[i++]= wMu2;
//   spectrum.waves[i++]= wN2;
  spectrum.waves[i++]= wNu2;
//   spectrum.waves[i++]= wM2;
  spectrum.waves[i++]= wLa2;
  spectrum.waves[i++]= wL2;
  spectrum.waves[i++]= wT2;
  spectrum.waves[i++]= wS2;
  spectrum.waves[i++]= wR2;
//   spectrum.waves[i++]= wK2;
  spectrum.waves[i++]= wKJ2;
  spectrum.waves[i++]= wOQ2;
  spectrum.waves[i++]= w2MK2;
  spectrum.waves[i++]= wMSK2;
  ///\note : MSK2 and OP2 have the same pulsation but a different nodal correction and are both in the list
  spectrum.waves[i++]= wMSN2;
  spectrum.waves[i++]= w2SM2;
  spectrum.waves[i++]= wM_SK_2;
  spectrum.waves[i++]= wM_KS_2;
  spectrum.waves[i++]= wMKS2;
  spectrum.waves[i++]= wOP2;
  spectrum.waves[i++]= wMNS2;
  spectrum.waves[i++]= w2SN2;
  spectrum.waves[i++]= w2NS2;

  spectrum.waves[i++]= wM3;
  spectrum.waves[i++]= wS3;
  spectrum.waves[i++]= wMK3;
  spectrum.waves[i++]= w2MK3;
  ///\note : 2MK3 and MO3 have the same pulsation but a different nodal correction and are both in the list
  spectrum.waves[i++]= wSO3;
  spectrum.waves[i++]= wSK3;
  spectrum.waves[i++]= wMO3;

  spectrum.waves[i++]= wN4;
  spectrum.waves[i++]= wMN4;
  spectrum.waves[i++]= wM4;
  spectrum.waves[i++]= wMS4;
  spectrum.waves[i++]= wMK4;
  spectrum.waves[i++]= wS4;
  spectrum.waves[i++]= wSN4;
  spectrum.waves[i++]= w3MS4;
  spectrum.waves[i++]= wSK4;
  //NOTE : wR4 and wSK4 seam to have the same pulsation

  spectrum.waves[i++]= w2MN6;
  spectrum.waves[i++]= wM6;
  spectrum.waves[i++]= w2MS6;
  spectrum.waves[i++]= w2MK6;
  spectrum.waves[i++]= wMSN6;
  spectrum.waves[i++]= w2SM6;
  spectrum.waves[i++]= wMSK6;

  spectrum.waves[i++]= w3MS8;

  tides_compound();

  spectrum.waves[i++]= wST1;
  spectrum.waves[i++]= wST2;
  spectrum.waves[i++]= wST3;
  spectrum.waves[i++]= wSNK2;

  //in SHOM list and missing from tides.def
  ///\note K1+O1=M2,K1+P1=S2,M1+O1=N2
  spectrum.waves[i++]= compound.init("3OK1").add(wO1,3).add(wK1,-2).realize();
//  spectrum.waves[i++]= compound.init("3OK1").add(wO1,3).add(wK1,-2).realize(); ???
  spectrum.waves[i++]= compound.init("MS1").add(wM2,1).add(wS1,-1).realize();
  //"A19"
  spectrum.waves[i++]= compound.init("2MN2S2").add(wM2,2).add(wN2,1).add(wS2,-2).realize();
  spectrum.waves[i++]= compound.init("3M2S2").add(wM2,3).add(wS2,-2).realize();
  spectrum.waves[i++]= compound.init("MNuS2").add(wM2,1).add(wNu2,1).add(wS2,-1).realize();
  //spectrum.waves[i++]= compound.init("2MK2").add(wM2,2).add(wK2,-1).realize();
  spectrum.waves[i++]= compound.init("NKM2").add(wN2,1).add(wK2,1).add(wM2,-1).realize();
  spectrum.waves[i++]= compound.init("SKM2").add(wS2,1).add(wK2,1).add(wM2,-1).realize();

//  spectrum.waves[i++]= compound.init("2S2N2").add(wS2,2).add(w2N2,-1).realize();
  spectrum.waves[i++]= compound.init("2SMu2").add(wS2,2).add(wMu2,-1).realize();
  
  spectrum.waves[i++]= compound.init("MNu4").add(wM2,1).add(wNu2,1).realize();
  spectrum.waves[i++]= compound.init("ML4").add(wM2,1).add(wL2,1).realize();
  
  spectrum.waves[i++]= compound.init("S5").add(wS1,1).add(wS2,2).realize();
  
  spectrum.waves[i++]= compound.init("2ML6").add(wM2,2).add(wL2,1).realize();
  spectrum.waves[i++]= compound.init("MSL6").add(wM2,1).add(wS2,1).add(wL2,1).realize();
  spectrum.waves[i++]= compound.init("MNK6").add(wM2,1).add(wN2,1).add(wK2,1).realize();

  spectrum.waves[i++]= compound.init("3MSN6").add(wM2,3).add(wS2,1).add(wN2,-1).realize();
  spectrum.waves[i++]= compound.init("3MNS6").add(wM2,3).add(wS2,-1).add(wN2,1).realize();
  
  spectrum.waves[i++]= compound.init("3MSK6").add(wM2,3).add(wS2,1).add(wK2,-1).realize();
  spectrum.waves[i++]= compound.init("3MKS6").add(wM2,3).add(wS2,-1).add(wK2,1).realize();

  spectrum.waves[i++]= compound.init("3MNK6").add(wM2,3).add(wN2,1).add(wK2,-1).realize();
  spectrum.waves[i++]= compound.init("3MKN6").add(wM2,3).add(wN2,-1).add(wK2,1).realize();

  spectrum.waves[i++]= compound.init("3MNL6").add(wM2,3).add(wN2,1).add(wL2,-1).realize();
  spectrum.waves[i++]= compound.init("3MLN6").add(wM2,3).add(wN2,-1).add(wL2,1).realize();

  spectrum.waves[i++]= compound.init("2MNu6").add(wM2,2).add(wNu2,1).realize();
  spectrum.waves[i++]= compound.init("MSNu6").add(wM2,1).add(wS2,1).add(wNu2,1).realize();

  spectrum.waves[i++]= compound.init("S6").add(wS2,3).realize();
  
  spectrum.waves[i++]= compound.init("S7").add(wS1,1).add(wS2,3).realize();
  
  spectrum.waves[i++]= compound.init("M8").add(wM2,4).realize();
  spectrum.waves[i++]= compound.init("S8").add(wS2,4).realize();
  spectrum.waves[i++]= compound.init("N8").add(wN2,4).realize();
  
  spectrum.waves[i++]= compound.init("S9").add(wS1,1).add(wS2,4).realize();
  
  spectrum.waves[i++]= compound.init("3MS8").add(wM2,3).add(wS2,1).realize();
  spectrum.waves[i++]= compound.init("3MN8").add(wM2,3).add(wN2,1).realize();
  spectrum.waves[i++]= compound.init("4MNS8").add(wM2,4).add(wN2,1).add(wS2,-1).realize();
  spectrum.waves[i++]= compound.init("3MNu8").add(wM2,3).add(wNu2,1).realize();
  spectrum.waves[i++]= compound.init("3MK8").add(wM2,3).add(wK2,1).realize();
  spectrum.waves[i++]= compound.init("3ML8").add(wM2,3).add(wL2,1).realize();
  
  spectrum.waves[i++]= compound.init("2M2S8").add(wM2,2).add(wS2,2).realize();
  
  spectrum.waves[i++]= compound.init("2MSN8").add(wM2,2).add(wS2,1).add(wN2,1).realize();
  spectrum.waves[i++]= compound.init("2MSK8").add(wM2,2).add(wS2,1).add(wK2,1).realize();

  spectrum.waves[i++]= compound.init("2M2S8").add(wM2,2).add(wS2,2).realize();
  spectrum.waves[i++]= compound.init("2M2N8").add(wM2,2).add(wN2,2).realize();
  spectrum.waves[i++]= compound.init("2M2K8").add(wM2,2).add(wK2,2).realize();
  
  spectrum.waves[i++]= compound.init("MSNK8").add(wM2,1).add(wS2,1).add(wN2,1).add(wK2,1).realize();
  spectrum.waves[i++]= compound.init("2SKN8").add(wS2,2).add(wK2,1).add(wN2,1).realize();

//  spectrum.waves[i++]= compound.init("2MSN8").add(wM2,2).add(wS2,1).add(wN2,1).realize();

  spectrum.waves[i++]= compound.init("S10").add(wS2,5).realize();
  
  spectrum.waves[i++]= compound.init("S11").add(wS1,1).add(wS2,5).realize();
  
  spectrum.waves[i++]= compound.init("MQ3").add(wM2,1).add(wQ1,1).realize();
  spectrum.waves[i++]= compound.init("MS3").add(wM2,1).add(wS1,1).realize();
  //"A87"
  spectrum.waves[i++]= compound.init("SP3").add(wS2,1).add(wP1,1).realize();
  spectrum.waves[i++]= compound.init("K3").add(wK1,3).realize();
  
  spectrum.waves[i++]= compound.init("MA4").add(wM4,1).add(wSa,-1).realize();
  spectrum.waves[i++]= compound.init("2MKS4").add(wM2,2).add(wK2,1).add(wS2,-1).realize();
  
  spectrum.waves[i++]= compound.init("2MQ5").add(wM2,2).add(wQ1,1).realize();
  spectrum.waves[i++]= compound.init("2MO5").add(wM2,2).add(wO1,1).realize();
  spectrum.waves[i++]= compound.init("2NK5").add(wN2,2).add(wK1,1).realize();
  spectrum.waves[i++]= compound.init("3MS5").add(wM2,3).add(wS1,-1).realize();
  spectrum.waves[i++]= compound.init("3MP5").add(wM2,3).add(wP1,-1).realize();
  spectrum.waves[i++]= compound.init("M5").add(wM1,5).realize();
  spectrum.waves[i++]= compound.init("2MP5").add(wM2,2).add(wP1,1).realize();
  spectrum.waves[i++]= compound.init("2MS5").add(wM2,2).add(wS1,1).realize();
  spectrum.waves[i++]= compound.init("2MK5").add(wM2,2).add(wK1,1).realize();
  
  spectrum.waves[i++]= compound.init("NSK5").add(wN2,1).add(wS2,1).add(wK1,1).realize();
  spectrum.waves[i++]= compound.init("3MQ5").add(wM2,3).add(wQ1,-1).realize();
  spectrum.waves[i++]= compound.init("MSP5").add(wM2,1).add(wS2,1).add(wP1,1).realize();
  spectrum.waves[i++]= compound.init("MSK5").add(wM2,1).add(wS2,1).add(wK1,1).realize();

  spectrum.waves[i++]= compound.init("2MNS4").add(wM2,2).add(wN2,1).add(wS2,-1).realize();
  spectrum.waves[i++]= compound.init("NK4").add(wN2,1).add(wK2,1).realize();
  spectrum.waves[i++]= compound.init("2MSN4").add(wM2,2).add(wS2,1).add(wN2,-1).realize();

  spectrum.waves[i++]= compound.init("4MK6").add(wM2,4).add(wK2,-1).realize();
  spectrum.waves[i++]= compound.init("4MS6").add(wM2,4).add(wS2,-1).realize();
  //spectrum.waves[i++]= compound.init("").add(w,1).add(w,1).realize();

  if(spectrum.nmax<i) TRAP_ERR_EXIT(ENOEXEC,"wrong prior number of waves spectrum.init(%d<%d)\n",spectrum.nmax-1,i-1);
  spectrum.n=i;

  for (i=0; i<spectrum.n; i++) {
    if(spectrum.waves[i].omega == 0.) {
      spectrum.waves[i].init();
/* *----------------------------------------------------------------------------
      */
      spectrum.waves[i].aliased = spectrum.waves[i].omega;
      }
/* printf ("wave: %10s, pulsation: %12.6f degrees/h \n",
                                     spectrum.waves[i].name,spectrum.waves[i].omega);*/
    }

  return(spectrum);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void tide_wavearray_initC(const char* filename, spectrum_t *spectrum_list)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int code,i,tmp,addnul=0;
  char a;
  FILE  *spectrefile=NULL;

  if ((spectrefile=fopen(filename,"r")) == NULL) {
    STDERR_BASE_LINE("wave file not found : %s \n",filename);
    exit(-1);
    }

  fscanf(spectrefile,"%d ",&tmp);
  spectrum_list->n=tmp;

  for(i=0;i<spectrum_list->n;i++) {
    fscanf(spectrefile,"%d ",&code);
    if(code==999) {
      spectrum_list->waves[i]=wNUL;
      spectrum_list->waves[i].nT=addnul;
      addnul++;
      fscanf(spectrefile,"%s %c %lf ",&spectrum_list->waves[i].name,&a,&spectrum_list->waves[i].omega);
      }
    do {a=fgetc(spectrefile);} while ( a!= '\n');

/*-------------------------- Longues periodes ---------------------------*/

   if(code == 47) (spectrum_list->waves[i]=wSa);
   else if(code == 46) (spectrum_list->waves[i]=wSsa);
   else if(code == 41) (spectrum_list->waves[i]=wMSm);
   else if(code == 38) (spectrum_list->waves[i]=wMm);
   else if(code == 39) (spectrum_list->waves[i]=wMSf);
   else if(code == 40) (spectrum_list->waves[i]=wMf);
   else if(code == 44) (spectrum_list->waves[i]=wMStm);
   else if(code == 42) (spectrum_list->waves[i]=wMtm);
   else if(code == 43) (spectrum_list->waves[i]=wMqm);
   else if(code == 45) (spectrum_list->waves[i]=wMSqm);

/*------------------------------- Diurnal -------------------------------*/

   else if(code == 67) (spectrum_list->waves[i]=w2Q1);
   else if(code == 68) (spectrum_list->waves[i]=wSig1);
   else if(code == 27) (spectrum_list->waves[i]=wQ1);
   else if(code == 69) (spectrum_list->waves[i]=wRo1);
   else if(code == 1) (spectrum_list->waves[i]=wO1);

 /*  else if(code == 121) (spectrum_list->waves[i]=wTau1);*/

   else if(code == 74) (spectrum_list->waves[i]=wM1);
   else if(code == 70) (spectrum_list->waves[i]=wKi1);
   else if(code == 71) (spectrum_list->waves[i]=wPi1);
   else if(code == 2) (spectrum_list->waves[i]=wP1);
   else if(code == 3) (spectrum_list->waves[i]=wK1);
   else if(code == 0) (spectrum_list->waves[i]=wPsi1);
   else if(code == 72) (spectrum_list->waves[i]=wPhi1);
   else if(code == 26) (spectrum_list->waves[i]=wS1);
   else if(code == 73) (spectrum_list->waves[i]=wTta1);
   else if(code == 29) (spectrum_list->waves[i]=wJ1);
   else if(code == 28) (spectrum_list->waves[i]=wOO1);

  /* else if(code == 170) (spectrum_list->waves[i]=wNu1);*/

   /* else if(code == 75) (spectrum_list->waves[i]=wM12);*/

   else if(code == 49) (spectrum_list->waves[i]=wMP1);
   else if(code == 48) (spectrum_list->waves[i]=wSO1);
   else if(code == 76) (spectrum_list->waves[i]=wKQ1);

/*---------------------------- Semi-diurnal -----------------------------*/

   /*else if(code == 200) (spectrum_list->waves[i]=wEps2);*/

   else if(code == 4) (spectrum_list->waves[i]=wE2);
   else if(code == 5) (spectrum_list->waves[i]=w2N2);
   else if(code == 6) (spectrum_list->waves[i]=wMu2);
   else if(code == 7) (spectrum_list->waves[i]=wN2);
   else if(code == 8) (spectrum_list->waves[i]=wNu2);

   /*else if(code == 230) (spectrum_list->waves[i]=wGamma2);*/

   /*else if(code == 231) (spectrum_list->waves[i]=wA2);*/

   else if(code == 9) (spectrum_list->waves[i]=wM2);

   /*else if(code == 233) (spectrum_list->waves[i]=wDelta2);*/

   else if(code == 10) (spectrum_list->waves[i]=wLa2);
   else if(code == 11) (spectrum_list->waves[i]=wL2);
   else if(code == 12) (spectrum_list->waves[i]=wT2);
   else if(code == 13) (spectrum_list->waves[i]=wS2);
   else if(code == 50) (spectrum_list->waves[i]=wR2);
   else if(code == 14) (spectrum_list->waves[i]=wK2);
   else if(code == 77) (spectrum_list->waves[i]=wKJ2);

   /*else if(code == 260) (spectrum_list->waves[i]=wEta2);*/

   else if(code == 51) (spectrum_list->waves[i]=wOQ2);
   else if(code == 65) (spectrum_list->waves[i]=w2MK2);
   else if(code == 31) (spectrum_list->waves[i]=wMSK2);
   else if(code == 15) (spectrum_list->waves[i]=wMSN2);
   else if(code == 16) (spectrum_list->waves[i]=w2SM2);
   else if(code == 37) (spectrum_list->waves[i]=wM_SK_2);
   else if(code == 36) (spectrum_list->waves[i]=wM_KS_2);
   else if(code == 30) (spectrum_list->waves[i]=wMKS2);
   else if(code == 100) (spectrum_list->waves[i]=wOP2);
   else if(code == 101) (spectrum_list->waves[i]=wMNS2);

/*---------------------------- 3rd-diurnal -----------------------------*/

   else if(code == 34) (spectrum_list->waves[i]=wM3);
   else if(code == 35) (spectrum_list->waves[i]=wS3);
   else if(code == 24) (spectrum_list->waves[i]=wMK3);
   else if(code == 25) (spectrum_list->waves[i]=w2MK3);
   else if(code == 53) (spectrum_list->waves[i]=wSO3);
   else if(code == 54) (spectrum_list->waves[i]=wSK3);
   else if(code == 102) (spectrum_list->waves[i]=wMO3);

/*---------------------------- 4th-diurnal -----------------------------*/

   else if(code == 33) (spectrum_list->waves[i]=wN4);
   else if(code == 17) (spectrum_list->waves[i]=wMN4);
   else if(code == 18) (spectrum_list->waves[i]=wM4);
   else if(code == 19) (spectrum_list->waves[i]=wMS4);
   else if(code == 20) (spectrum_list->waves[i]=wMK4);
   else if(code == 56) (spectrum_list->waves[i]=wS4);
   else if(code == 55) (spectrum_list->waves[i]=wSN4);
   else if(code == 58) (spectrum_list->waves[i]=w3MS4);
   else if(code == 103) (spectrum_list->waves[i]=wSK4);

   /*else if(code == 500) (spectrum_list->waves[i]=wM5);*/

/*---------------------------- 6th-diurnal -----------------------------*/

   else if(code == 21) (spectrum_list->waves[i]=w2MN6);
   else if(code == 22) (spectrum_list->waves[i]=wM6);
   else if(code == 59) (spectrum_list->waves[i]=w2MS6);
   else if(code == 60) (spectrum_list->waves[i]=w2MK6);
   else if(code == 23) (spectrum_list->waves[i]=wMSN6);
   else if(code == 104) (spectrum_list->waves[i]=w2SM6);
   else if(code == 105) (spectrum_list->waves[i]=wMSK6);

   /*---------------------------- 8th-diurnal -----------------------------*/

   else if(code == 61) (spectrum_list->waves[i]=w3MS8);
   else if(code == 187) (spectrum_list->waves[i]=wmean);

   else {
     printf ("wave: %10s not available \n", spectrum_list->waves[i].name);
     }

   }
  fclose(spectrefile);

  for(i=0;i<spectrum_list->n;i++) {
    if(spectrum_list->waves[i].omega == 0.) {
      spectrum_list->waves[i].init();
      spectrum_list->waves[i].aliased = spectrum_list->waves[i].omega;
      }
#if VERBOSE > 1
      printf ("wave: %10s, period: %15.11f days \n",
              spectrum_list->waves[i].name,  tide_period(spectrum_list->waves[i]) / 24.);
#endif
  }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void tide_wavearray_initC(spectrum_t *spectrum_list, int nwave)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  char filename[256] ;
  sprintf(filename,"ondes%d.dut",nwave);
  tide_wavearray_initC(filename, spectrum_list);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int index_from_name(const spectrum_t & s, const char *name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  for(int k = 0; k < s.n; k++) {
    if(strcmp(s.waves[k].name, name) == 0) {
      return(k);
    }
  }
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  tidal_wave wave_from_name(const char *name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  tidal_wave wave = wNUL;
  spectrum_t reference_spectrum = initialize_tide();

  int k;
  
  k=index_from_name(reference_spectrum,name);
  
  if(k>=0)
    wave = reference_spectrum.waves[k];
  
  reference_spectrum.destroy();
  
  return(wave);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

tidal_wave wave_from_code(int code)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  tidal_wave wave;

/*-------------------------- Longues periodes ---------------------------*/

  switch(code) {

   case 47:
     wave=wSa;
     break;

   case 46:
     wave=wSsa;
     break;

   case 41:
     wave=wMSm;
     break;

   case 38:
     wave=wMm;
     break;

   case 39:
     wave=wMSf;
     break;

   case 40:
     wave=wMf;
     break;

   case 44:
     wave=wMStm;
     break;

   case 42:
     wave=wMtm;
     break;

   case 43:
     wave=wMqm;
     break;

   case 45:
     wave=wMSqm;
     break;

/*------------------------------- Diurnal -------------------------------*/

   case 67:
     wave=w2Q1;
     break;

   case 68:
     wave=wSig1;
     break;

   case 27:
     wave=wQ1;
     break;

   case 69:
     wave=wRo1;
     break;

   case 1:
     wave=wO1;
     break;

 /*  case 121:
     wave=wTau1;
     break;*/

   case 74:
     wave=wM1;
     break;

   case 70:
     wave=wKi1;
     break;

   case 71:
     wave=wPi1;
     break;

   case 2:
     wave=wP1;
     break;

   case 3:
     wave=wK1;
     break;

   case 0:
     wave=wPsi1;
     break;

   case 72:
     wave=wPhi1;
     break;

   case 26:
     wave=wS1;
     break;

   case 73:
     wave=wTta1;
     break;

   case 29:
     wave=wJ1;
     break;

   case 28:
     wave=wOO1;
     break;

  /* case 170:
     wave=wNu1;
     break;*/

   /* case 75:
     wave=wM12;
     break;*/

   case 49:
     wave=wMP1;
     break;

   case 48:
     wave=wSO1;
     break;

   case 76:
     wave=wKQ1;
     break;

/*---------------------------- Semi-diurnal -----------------------------*/

   /*case 200:
     wave=wEps2;
     break;*/

   case 4:
     wave=wE2;
     break;

   case 5:
     wave=w2N2;
     break;

   case 6:
     wave=wMu2;
     break;

   case 7:
     wave=wN2;
     break;

   case 8:
     wave=wNu2;
     break;

   /*case 230:
     wave=wGamma2;
     break;*/

   /*case 231:
     wave=wA2;
     break;*/

   case 9:
     wave=wM2;
     break;

   /*case 233:
     wave=wDelta2;
     break;*/

   case 10:
     wave=wLa2;
     break;

   case 11:
     wave=wL2;
     break;

   case 12:
     wave=wT2;
     break;

   case 13:
     wave=wS2;
     break;

   case 50:
     wave=wR2;
     break;

   case 14:
     wave=wK2;
     break;

   case 77:
     wave=wKJ2;
     break;

   /*case 260:
     wave=wEta2;
     break;*/

   case 51:
     wave=wOQ2;
     break;

   case 65:
     wave=w2MK2;
     break;

   case 31:
     wave=wMSK2;
     break;

   case 15:
     wave=wMSN2;
     break;

   case 16:
     wave=w2SM2;
     break;

   case 37:
     wave=wM_SK_2;
     break;

   case 36:
     wave=wM_KS_2;
     break;

   case 30:
     wave=wMKS2;
     break;

   case 100:
     wave=wOP2;
     break;

   case 101:
     wave=wMNS2;
     break;

/*---------------------------- 3rd-diurnal -----------------------------*/

   case 34:
     wave=wM3;
     break;

   case 35:
     wave=wS3;
     break;

   case 24:
     wave=wMK3;
     break;

   case 25:
     wave=w2MK3;
     break;

   case 53:
     wave=wSO3;
     break;

   case 54:
     wave=wSK3;
     break;

   case 102:
     wave=wMO3;
     break;

/*---------------------------- 4th-diurnal -----------------------------*/

   case 33:
     wave=wN4;
     break;

   case 17:
     wave=wMN4;
     break;

   case 18:
     wave=wM4;
     break;

   case 19:
     wave=wMS4;
     break;

   case 20:
     wave=wMK4;
     break;

   case 56:
     wave=wS4;
     break;

   case 55:
     wave=wSN4;
     break;

   case 58:
     wave=w3MS4;
     break;

   case 103:
     wave=wSK4;
     break;

   /*case 500:
     wave=wM5;
     break;*/

/*---------------------------- 6th-diurnal -----------------------------*/

   case 21:
     wave=w2MN6;
     break;

   case 22:
     wave=wM6;
     break;

   case 59:
     wave=w2MS6;
     break;

   case 60:
     wave=w2MK6;
     break;

   case 23:
     wave=wMSN6;
     break;

   case 104:
     wave=w2SM6;
     break;

   case 105:
     wave=wMSK6;
     break;

   /*---------------------------- 8th-diurnal -----------------------------*/

   case 61:
     wave=w3MS8;
     break;


   case 187:
     wave=wmean;
     break;

   }

  return(wave);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int initialize_omega(spectrum_t *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  spectrum_t spectrum;

  for (i=0; i<s->n; i++) {
    if(s->waves[i].omega == 0.) {
      s->waves[i].init();
      s->waves[i].aliased = s->waves[i].omega;
      }
/*printf ("wave: %10s, pulsation: %12.6f degrees/h \n",
                                     spectrum.waves[i].name,spectrum.waves[i].omega);*/
    }

  return(0);
}
