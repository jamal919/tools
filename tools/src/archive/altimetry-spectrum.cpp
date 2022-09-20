
/**************************************************************************

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

***************************************************************************/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>

#include "tools-structures.h"

#include "spectrum.h"
#include "rutin.h"

#ifndef VERBOSE
#define VERBOSE 1
#endif

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

spectrum_t spectrum_init_ERS(void)

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

  return(spectrum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

spectrum_t spectrum_init_GFO(void)

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

  return(spectrum);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

spectrum_t spectrum_init_TP_interleaved(void)

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

  return(spectrum);
}


