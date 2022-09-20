
/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2012

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France

***************************************************************************/

#include "config.h"
#include "tugo-prototypes.h"
#include "solverlib.h"

int harmonic_nsave = 0;

extern const tidal_wave wS1;

extern const tidal_wave wM3;

extern const tidal_wave wM4;
extern const tidal_wave wN4;
extern const tidal_wave wM4;
extern const tidal_wave wS4;
extern const tidal_wave wMN4;

extern const tidal_wave wMS4;
extern const tidal_wave wMK4;
extern const tidal_wave wSN4;
extern const tidal_wave w3MS;
extern const tidal_wave wSK4;

extern const tidal_wave wM6;

extern const double omega_T;
extern const double omega_s; 
extern const double omega_h;
extern const double omega_p;
extern const double omega_n;
extern const double omega_p1;

extern int tide_complex(double, spectrum_t, float *, float *, int);
extern int *admittance_check();
extern double tide_omega(tidal_wave,const char *);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int addwave(spectrum_t *list, tidal_wave wave, int prescribed)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, nadd, k, n, nwave, nwmax;
  char *in, *out;
  size_t size;

/* to be implemanted later
  for (k=0;k<list->n;k++) {
    if(strcmp(list->wave[k].name,wave.name)==0) return(1);
    }
*/
  size = sizeof(tidal_wave);

  n = list->n;
  nwmax = list->nmax;

  /* /!\ Very bad choice for the data structure: it would be more efficient,
     easyer to implement and induce fewer bugs to use a chained list from the
     standard library. Or, if it is really a rezisable array that is wanted,
     use a vector from the standard library. /!\ */
  if(n >= nwmax) {
    nwmax+=10;
    /* Never use realloc with a pointer that has been defined by new: it is
       impossible to predict what will happen. */
    /*list->wave       = (tidal_wave *) realloc(list->wave, nwmax * size);
    list->prescribed = (int *) realloc(list->prescribed, nwmax * sizeof(int));*/
    tidal_wave *tmpWave = new tidal_wave[nwmax];
    int *tmpPrescribed  = new int[nwmax];
    for (size_t j = 0; j < list->n; j++) {
      tmpWave[j] = list->waves[j];
      tmpPrescribed[j] = list->prescribed[j];
      }
    delete[] list->waves;
    list->waves = tmpWave;
    delete[] list->prescribed;
    list->prescribed = tmpPrescribed;
    }

//  memmove(&(list->waves[n]), &wave, size);
  list->waves[n]=wave;
  list->prescribed[n]=prescribed;

  list->nmax = nwmax;
  list->n = n+1;

  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LoadAnalysisList()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, idum, j, k, n, nitems;
  float fdum;
  char name[1024];
  FILE *in;
  tidal_wave *dum;
  size_t size;

  size = sizeof(tidal_wave);

  if((in = fopen(WaveFile, "r")) == NULL) {
    printf("#LoadAnalysisList: cannot open %s \n", WaveFile);
    return (-1);
    }

  nitems = fscanf(in, "%d", &n);
  if(nitems != 1) {
    printf("#LoadAnalysisList: cannot read properly %s \n", WaveFile);
    return (-1);
    }

  fprintf(stdout, " %d tidal analysed waves.\n", n);

  AnalysisList.waves = new tidal_wave[n];
  AnalysisList.n = n;

  j = -1;
  while(!feof(in) && j < n - 1) {
    nitems = fscanf(in, "%s", name);
    if(nitems != 1) {
      printf("#LoadAnalysisList: cannot read properly %s \n", WaveFile);
      return (-1);
      }
    for(k = 0; k < reference_spectrum.n; k++) {
      if(strcmp(reference_spectrum.waves[k].name, name) == 0) {
        break;
        }
      }
    j++;
    printf("wave %3d on %d :  %s \n", j + 1, n, name);
    memmove(&(AnalysisList.waves[j]), &(reference_spectrum.waves[k]),size);
    printf("wave %3d:  %s %s \n", j, AnalysisList.waves[j].name, reference_spectrum.waves[k].name);
    }

  fclose(in);
  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int load_wave_list(char *filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, idum, j, k, n, nitems;
  float fdum;
  char name[1024], text[1024];
  FILE *in;
  tidal_wave *dum;
  size_t size;
  int *method;

  size = sizeof(tidal_wave);

  WaveList.n    = 0;
  WaveList.nmax = 0;

  if((in = fopen(filename, "r")) == NULL) {
    printf("load_wave_list: cannot open %s \n", filename);
    return (-1);
    }

  while(!feof(in)) {
    nitems = fscanf(in, "%s", text);
    if(strcmp(text, "#SPECTRUM") == 0)
      goto next;
    }
  printf("load_wave_list: missing <#SPECTRUM> keyword in %s \n",filename);
  return (-1);

next:
  nitems = fscanf(in, "%d", &n);
  if(nitems != 1) {
    printf("load_wave_list: cannot read properly %s \n", filename);
    return (-1);
    }

  WaveList.n    = 0;
  WaveList.nmax = n+100;
  WaveList.waves = new tidal_wave[WaveList.nmax];
  WaveList.prescribed = new int[WaveList.nmax];

  for(j = 0; j < n; j++) {
    nitems = fscanf(in, "%s", name);
    if(nitems != 1) {
      printf("load_wave_list: cannot read properly %s \n", WaveFile);
      return (-1);
      }
    for(k = 0; k < reference_spectrum.n; k++) {
      if(strcmp(reference_spectrum.waves[k].name, name) == 0) {
        break;
        }
      }
    if(k != reference_spectrum.n) {
      addwave(&WaveList, reference_spectrum.waves[k],-1);
      }
    else {
      printf("load_wave_list: cannot include %s \n", name);
      }
    }
  fclose(in);
  return (0);
}
