
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

#include  "tugo-prototypes.h"
#include  "poc-netcdf.hpp"
#include "map.h"

#include "legend.def"
#include "legend.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int diag_BC(tidaldata_t *data, int ndata, int nwaves)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int neighbour,nlegends;
  int k,j,l,m,n,status;
  legend_t   *legends;
  legend02_t *legends02;
  char *comments;
  char lgdfile[1024];
  edge_t *edgePtr;
  int orientation;

  legends=new legend_t[1];

  l=0;
  legends[l].ID=l;
  legends[l].Type=LGD_GRAPH;
  legends[l].T=0;
  legends[l].P=0;
  legends[l].X=0;
  legends[l].Y=0;

  comments=new char[1024];
  comments[0]=0;

  strcat(comments,"#  0: node index \n");
  strcat(comments,"#  1: longitude of data \n");
  strcat(comments,"#  2: latitude  of data \n");
  strcat(comments,"#  3->0: normal-x\n");
  strcat(comments,"#  4->1: normal-y\n");

  legends02=new legend02_t[1];

  n=0;
  legends02[n].ID=n;
  legends02[n].np=ndata;
  legends02[n].nz=2*nwaves;
  legends02[n].pen0=0;
  legends02[n].pen1=0;
  legends02[n].points=new point_t[legends02[n].np];
  legends[n].ptr=(char *) &(legends02[n]);
  for (m=0;m<legends02[n].np;m++) {
    legends02[n].points[m].z=new float[legends02[n].nz];
    }
  for (l=0;l<ndata;l++) {
    legends02[n].points[l].N=l;
    legends02[n].points[l].t=TidalData[l].lon*180./M_PI;
    legends02[n].points[l].p=TidalData[l].lat*180./M_PI;
    j=0;
/**----------------------------------------------------------------------------
    right-hande side normal*/
    for(k=0;k<nwaves;k++) {
      legends02[n].points[l].z[j++]=TidalData[l].Ha[k];
      legends02[n].points[l].z[j++]=TidalData[l].Hg[k];
      }
    }

  sprintf(lgdfile,"diagnostic.lgd");

  nlegends=1;

  status=lgd_save(lgdfile, legends, nlegends,NULL,comments);

  for (m=0;m<legends02[n].np;m++) {
    delete[] legends02[n].points[m].z;
    }
  delete[] legends02[n].points;

  delete[] legends02;
  delete[] legends;

  delete[] comments;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int GetBdyObs(const char *filename, double scale, int verbose, int debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, idum, j, k, nitems,status=0;
  float fdum;
  char name[1024],units[256],msg[1024],char_dum[1024];
  char *s;
  double factor;
  FILE *in1;
  tidal_wave *dum;
  
  verbose=( verbose && (gCPU_ID==gCPU_MASTER));
  
  if(verbose) {
    printf(SEPARATOR_1);
    printf("load tidal open boundary conditions\n");
    }

  in1 = fopen(filename, "r");
  
  if(in1 == NULL) {
    sprintf(msg,"Cannot open tidal boundary file: %s\n", filename);
    status = -1;
    check_error(status, msg, __LINE__, __FILE__, 1);
    }

  s=fgets(char_dum, 1024, in1);
  nitems=sscanf(char_dum,"%d %d %s", &NTidalData, &NTidalWave,units);

  if(nitems==2) {
    if(verbose) printf("\nno units given in %s, assumes cm\n",filename);
    factor=0.01;
    }
  else {
    if((strcmp(units,"cm")==0) || (strcmp(units,"CM")==0)) {
      factor=0.01;
      }
    if((strcmp(units,"m")==0)||(strcmp(units,"M")==0)) {
      factor=1.0;
      }
    if(verbose) printf("\nunits given in %s: %s \n",filename,units);
    }

  if(verbose) fprintf(stdout, "%d tidal data points. %d tidal input waves\n", NTidalData, NTidalWave);

  if(NTidalData != 0) {
    TidalData = new tidaldata_t [NTidalData];
    if(TidalData == NULL)
      check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
    for(i = 0; i < NTidalData; i++) {
      TidalData[i].Ha = new float[NTidalWave] ();
      TidalData[i].Hg = new float[NTidalWave] ();
      TidalData[i].Ua = new float[NTidalWave] ();
      TidalData[i].Ug = new float[NTidalWave] ();
      TidalData[i].Va = new float[NTidalWave] ();
      TidalData[i].Vg = new float[NTidalWave] ();
      }
    }

  WaveList.n     = NTidalWave;
  WaveList.nmax  = 200;
  
  WaveList.waves      = new tidal_wave[WaveList.nmax];
  WaveList.prescribed = new int[WaveList.nmax];
  
  for(i = 0; i < WaveList.nmax; i++) {
    WaveList.prescribed[i]=-1;
    }

  j = -1;
  
  if(verbose) printf("Tidal boundary conditions found for waves : ");
  
  while(!feof(in1) && j < NTidalWave - 1) {
    s=fgets(char_dum, 1024, in1);
    if(s==0) {
      check_error(-1, "missing wave in OBC tidal file", __LINE__, __FILE__, 1);
      }
    sscanf(char_dum,"%s", name);
/**----------------------------------------------------------------------------
    Patch for name mismatch in FES atlas */
    if(strcmp("Msqm", name)==0) strcpy(name,"MSqm"); 
/**----------------------------------------------------------------------------
    indentify tidal wave from spectrum list */
    k=reference_spectrum.wave_index(name);
    j++;
/*------------------------------------------------------------------------------
    designed for non tidal frequency simulation*/
    if(k<0){
      nitems=sscanf(name,"%lg",&WaveList.waves[j].omega);
      if(nitems==1){
        //this is safer than a strcpy
        snprintf(WaveList.waves[j].name,TIDAL_WAVE_NAME_LEN,"%gdph",WaveList.waves[j].omega);
        k=reference_spectrum.n;
	WaveList.waves[j].Ap=-1;
        }
      }
    if(k<0) {
      fprintf(stdout, "->>>%s<<<<- not identified %d %d/%d\n", name,nitems,k,reference_spectrum.n);
      fflush(stdout);
      status=-1;
      check_error(status, "unknown tidal wave in OBC's file", __LINE__, __FILE__, 1);
      }

    if(k<reference_spectrum.n)
      memmove(&(WaveList.waves[j]), &(reference_spectrum.waves[k]),sizeof(tidal_wave));
    
    if(verbose) printf("%s ", WaveList.waves[j].name);

    gBoundary_FullTideAvailable=1;
    for(i = 0; i < NTidalData; i++) {
      s=fgets(char_dum, 1024, in1);
      nitems=sscanf(char_dum,"%lf %lf %f %f %f %f %f %f", &TidalData[i].lon, &TidalData[i].lat,
                                                     &(TidalData[i].Ha[j]), &(TidalData[i].Hg[j]),
                                                     &(TidalData[i].Ua[j]), &(TidalData[i].Ug[j]),
                                                     &(TidalData[i].Va[j]), &(TidalData[i].Vg[j]));
      TidalData[i].lon *= d2r;
      TidalData[i].lat *= d2r;
      TidalData[i].Ha[j] *= factor;
      TidalData[i].Hg[j] *= d2r;
      switch (nitems) {
        case 4:
          TidalData[i].Ua[j] = 999.9;
          TidalData[i].Ug[j] = 999.9;
          TidalData[i].Va[j] = 999.9;
          TidalData[i].Vg[j] = 999.9;
          gBoundary_FullTideAvailable=0;
          break;
        case 8:
          TidalData[i].Ua[j] *= factor;
          TidalData[i].Ug[j] *= d2r;
          TidalData[i].Va[j] *= factor;
          TidalData[i].Vg[j] *= d2r;
          gBoundary_FullTideAvailable=1;
          break;
        default:
          check_error(-1, "illegal format", __LINE__, __FILE__, 1);
          break;
        }
      }
    }

  if(verbose) printf("\n");
  fclose(in1);

  if(verbose) status=diag_BC(TidalData, NTidalData,NTidalWave);

  if(verbose)  {
    printf(SEPARATOR_2);
    }
    
  return(status);

}
