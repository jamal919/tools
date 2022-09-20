
/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "tides.h"
#include "fe.h"
#include "archive.h"
#include "poc-time.h"
#include "mgr.h"
#include "functions.h"

extern int sadcp_info(char *filename, double **time);
extern int sadcp_read(char *filename, float **u, float **v, int level);
extern int sadcp_readall(char *filename, float **u, float **v);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double dT;
  double tau;

  int i,j,k,l,nitems,w;
  int n,status;
  FILE *in;
  char *input=0,*keyword,*pathname=NULL,*s;
  char *atlas_directory=NULL,*atlas_convention=NULL;
  char *onde[10];
  int nonde;
  date_t start,final,reference,start_date;
  double *time, *z, *residuals, mean,mask=1.e+10;
  float  *u, *v;
  statistic_t stat;
  spectrum_t AnalysisList,PredictionList,solved;
  int detiding=0;
  hconstant_t prior;
  int nodal=1;
  double x,y,*tides;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 's' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d/%d",&start.day,&start.month,&start.year);
          break;

       case 'f' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d/%d",&final.day,&final.month,&final.year);
          break;

/* *----------------------------------------------------------------------
        naming convention for tidal atlases*/
        case 'c' :
          atlas_convention= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'w' :
          i=1;
          while(n+i< argc)  {
            onde[i-1]= strdup(argv[n+i]);
            i++;
            }
          n=n+i;
          nonde=i-1;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        if(input==NULL) {
          input= strdup(argv[n]);
          n++;
          }
        else {
          __OUT_BASE_LINE__("redundant option %s\n",keyword);
          exit(-1);
          }
        break;
      }
      free(keyword);
    }

string rootname;

  rootname=input;

size_t pos=rootname.find(".nc");

  rootname.erase(pos,4);

  if(pathname==NULL) pathname=strdup(".");

/* *-------------------------------------------------------------------------------------
  load time serie*/
//   n= sadcp_info(input, &time);
//   n= sadcp_read(input, &u, &v, 0);
//   n= sadcp_readall(input, &u, &v);


  AnalysisList.waves=new tidal_wave[100];
  AnalysisList.nmax=100;

  tide_wavearray_initC(&AnalysisList,20);

/*-------------------------------------------------------------------------------------
  HARMONIC ANALYSIS */
  residuals=new double[n];

  mgr_data_t *data;

  if(detiding==1) {
    prior=atlas_constants(atlas_convention, PredictionList,  x,  y);
    tides=new double[n];
    status=harmonic_predictionTS(tides, time, n, PredictionList, prior.a, prior.G, nodal);
    for (k=0;k<n;k++) {
      z[k]-=tides[k];
      }
    }

  nodal=1;
  double repetition=0.0;
  double maxratio=0.3;
  double averaging_time=0.0;
  data=harmonic_analysis(z, residuals, time, n, AnalysisList, averaging_time, repetition, solved, nodal, maxratio, stdout);

/*------------------------------------------------------------------------------
  recomposition*/
  if(detiding==1) {
    for(w=0;w<PredictionList.n;w++) {
//      k=index_from_name(spectrum_t s,char *name)
      }
    }

  vector<mgr_t> mgr;
  int nmgr=1;

  start=poctime_getdatecnes(time[0],   'd');
  final=poctime_getdatecnes(time[n-1], 'd');

  for(i=0;i<nmgr;i++) {
    mgr_t gauge;
    gauge.number=i;
    strcpy(gauge.name,"no_name");
    strcpy(gauge.origine,"undocumented");
    gauge.data=data;
    gauge.nwave=AnalysisList.n;
    gauge.loc.units=strdup("degrees");    sprintf(gauge.debut,"%2d/%2d/%4d",start.day,start.month,start.year);

    gauge.loc.lon=0.0;
    gauge.loc.lat=0.0;
    gauge.mindex=1;
    sprintf(gauge.debut,"%2d/%2d/%4d",start.day,start.month,start.year);
    sprintf(gauge.fin,"%2d/%2d/%4d",final.day,final.month,final.year);
    mgr.push_back(gauge);
    }


  string mgrfile;
  mgrfile=rootname+".mgr";
  status=mgr_save_ascii((char *) mgrfile.c_str() ,mgr);

  string residualsfile;
  residualsfile=rootname+".res";
  status=save_timeserie((char *)residualsfile.c_str(), residuals, mask, time, n);

  goto end;
error:
end:
  __OUT_BASE_LINE__(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(0);
}
