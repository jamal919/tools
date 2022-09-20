

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

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "poc-time.h"
#include "filter.h"
#include "tides.h"
#include "functions.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double mask=999.9;
  double dT[4],t0,origine,t1=0,t2=1.e+10;
  double *h[4],*t[4],*buffer;

  float  count, hours,w,dum;
  double start,finish,time,epsilon;

  int i,j,k,l,m,n,shift,status;
  int nr[4],nst,test=0;
  FILE *in,*out,*statfile;
  char name[256],tmp[256];
  char *model=NULL,*data=NULL, *input, *keyword,*translation=NULL;

  pressure_station *sample;
  statistic_t stat[3];
  int  year,month,day,hour;
  int start_year=1900,last_year=2100,start_month=1,last_month=12;

  fct_echo( argc, argv);

  fprintf(stderr,"%s  -starting computation *********\n",argv[0]);

  n=1;
  while (n < argc) {
    if(strcmp(argv[n],"-start") ==0) {
      status=sscanf(argv[n+1],"%d/%d",&start_month,&start_year);
      n++;
      n++;
      continue;
      }
    if(strcmp(argv[n],"-end") == 0) {
      status=sscanf(argv[n+1],"%d/%d",&last_month,&last_year);
      n++;
      n++;
      continue;
      }
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
/* *-----------------------------------------------------------------------------
        observation time serie (file)*/
        case 'd' :
          data= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *-----------------------------------------------------------------------------
        model time serie (file list)*/
        case 'm' :
          model= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
          break;
        }
        break;

      default:
        __OUT_BASE_LINE__("unknown option %s\n",keyword);
        exit(-1);
        break;
      }
      free(keyword);
    }

  statfile=fopen("statistics","a");

/* *-----------------------------------------------------------------------------
  read the station list*/
  in=fopen(model,"r");
  fscanf(in, "%d ", &nst);
  sample=(pressure_station *)malloc(nst*sizeof(pressure_station));
  for(k=0; k<nst;k++) {
    fscanf(in,"%s",sample[k].name);
    printf("add %s \n",sample[k].name);
    }
  fclose(in);

  origine=julian_day(1,1,1950);
  t1=julian_day(start_month,1,start_year)-origine;
  if(last_month==12)
    t2=julian_day(1,1,last_year+1)-origine;
  else
    t2=julian_day(last_month+1,1,last_year)-origine;

/* *-----------------------------------------------------------------------------
  read the observation file  */
  strcpy(name,data);
  printf("treating observed sea level from %s\n",name);
  status=matroos_load(name, t1, t2, &h[0], &t[0], &mask, &nr[0],&dT[0]);
  if(status!=0) {
    __OUT_BASE_LINE__("aborted (%s)\n",name);
    exit(-1);
    }

  for (k=0;k<nst;k++) {
/* *-----------------------------------------------------------------------------
    read the model elevation file  */
    strcpy(name,sample[k].name);
    status=matroos_load(name, t1, t2, &h[1], &t[1], &mask, &nr[1],&dT[1]);
    if(status!=0) {
      printf("aborted (%s)\n",name);
      continue;
      }
    buffer=(double *)malloc(nr[1]*sizeof(double));
    for (i=0; i< nr[1]; i++) {
/* *----------------------------------------------------------------------------
      convert in centimetres  */
      if(h[1][i] != mask) buffer[i]=h[1][i]*100.0;
      else buffer[i]=mask;
      }
    strcpy(tmp,name);
    strcat(tmp,".model.fft");
    fourier(tmp,buffer,mask,nr[1],dT[1],0);
    free(buffer);

/* *-----------------------------------------------------------------------------
    start analysis  */
    start =-1.e+10;
    finish=+1.e+10;

    for (i=0;i<2;i++) start =MAX(start,t[i][0]);
    for (i=0;i<2;i++) finish=MIN(finish,t[i][nr[i]-1]);

    if(start > finish) continue;

    dT[2]=MAX(dT[0],dT[1]);
    nr[2]=(int)( NINT((finish-start)/dT[2])+1 );
    h[2]=(double *)malloc(nr[2]*sizeof(double));
    t[2]=(double *)malloc(nr[2]*sizeof(double));

    epsilon=dT[2]/10.;

/* *-----------------------------------------------------------------------------
    build the common time serie: observations */
    n=-1;
    for (time=start; time <= finish+epsilon; time+=dT[2]) {
      n++;
      i=(int)( NINT((time-t[1][0])/dT[1]) );
      m=(int)( NINT((time-t[0][0])/dT[0]) );
      if(h[0][m] != mask) {
/* *----------------------------------------------------------------------------
        convert in centimetres */
        h[2][n]=h[0][m]*100.0;
        }
      else h[2][n]=mask;
      t[2][n]=t[1][i];
      }
    printf("sea level observation at %s: no correction -------------\n",sample[k].name);
    strcpy(name,sample[k].name);
    strcat(name,".obs");
    stat[0]=filtering(name,h[2],t[2],mask,nr[2],dT[2]);
    printf("#topex 0 %10s samples: %6d; mean: %8.4lf; variance: %8.4lf\n",
          sample[k].name,n,stat[0].mean,(stat[0].std)*(stat[0].std));

    strcpy(tmp,name);
    strcat(tmp,".fft");
    fourier(tmp,h[2],mask,nr[2],dT[2],0);

/* *-----------------------------------------------------------------------------
    build the common time serie: observations minus model */
    n=-1;
    for (time=start; time <= finish+epsilon; time+=dT[2]) {
      n++;
      i=(int)( NINT((time-t[1][0])/dT[1]) );
      m=(int)( NINT((time-t[0][0])/dT[2]) );
      if((h[0][m] != mask)&&(h[1][i] != mask)) {
        h[2][n]=(h[0][m]-h[1][i])*100.0;
        }
      else {
        h[2][n]=mask;
        }
      }
    printf("sea level observation at %s: model output corrected-------\n",sample[k].name);
    strcpy(name,sample[k].name);
    strcat(name,".dif");
    stat[1]=filtering(name,h[2],t[2],mask,nr[2],dT[2]);
    printf("#topex 1 %10s samples: %6d; mean: %8.4lf; variance: %8.4lf\n",
          sample[k].name,n,stat[1].mean,(stat[1].std)*(stat[1].std));

    strcpy(tmp,name);
    strcat(tmp,".fft");
    fourier(tmp,h[2],mask,nr[2],dT[2],0);

    fprintf(statfile,"%10s %6d %8.4lf %8.4lf %8.4lf %8.4lf\n",
          sample[k].name,n, stat[0].std,stat[1].std,
          (stat[0].std)*(stat[0].std), (stat[1].std)*(stat[1].std));
    for (i=1;i<2;i++)  free (h[i]);
    for (i=1;i<2;i++)  free (t[i]);
    }

  fclose(statfile);
  __ERR_BASE_LINE__("%s -computation complete ^^^^^^^^^\n",argv[0]);
  exit(0);
}
