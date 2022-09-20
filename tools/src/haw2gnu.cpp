

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "poc-time.h"
#include "filter.h"
#include "tides.h"
#include "mgr.h"
#include "functions.h"


/*-------------------------------------------------------------------------------*/

void save_hawaii(FILE *hawaii,double lat,double lon,double *residual,double *time,int numstation,char *stationname,date_t firstdate,int nmes,double t1)

/*-------------------------------------------------------------------------------*/
{
  int ndays,*gap,i,j,k,n,delay,year,firsthour;
  double *res,mask,*t,t0,data[12],skip;
  char filename[30],strcode[20],*newstationname;
  date_t dum;
  double epsilon=1.e-08,chk,delta;

  gap=(int *)malloc(nmes*sizeof(int));
  newstationname=(char *)malloc(20*sizeof(char));

  sprintf(strcode,"%03d",numstation);
  for(i=0;i<strlen(stationname);i++) newstationname[i]=stationname[i];
  for(i=strlen(stationname);i<8;i++) newstationname[i]=' ';
  newstationname[8]='\0';
  strcat(strcode,newstationname);
  mask=9999;

  free(newstationname);

  /*----------------------------------------------*/
  /* replacing blanks in time data by mask value  */
  /*----------------------------------------------*/

//  getcnesdate(time[0],&dum);    /* 1st date in timezone = GMT!*/
  dum=poctime_getdatecnes(time[0],'h');  /* 1st date in timezone = GMT!*/

  if (firstdate.year>dum.year) firsthour=0;
  else firsthour=floor(dum.second/3600.0+0.5);

/* *----------------------------------------------------------------------------
  detect holes in time serie*/
  n=nmes;
  for (i=1;i<nmes;i++) {
    delta=time[i]-time[i-1];
    chk=(delta-1.-epsilon)*(delta-1.+epsilon);
    if(chk>0) gap[i]=0;
#warning The only one usefull line in this loop follows
    gap[i]=1;
//    gap[i]=(int)(time[i]-time[i-1]);
    if (gap[i]!=1) {
      n+=gap[i];
      }
    }

  if(n>nmes) n=n+(12-n%12)+12+24;
  else n+=24;

  res=(double *)malloc(n*sizeof(double));
  t=  (double *)malloc(n*sizeof(double));

  for(i=0;i<n;i++) res[i]=mask;/*!!!!*/
  for (i=0;i<nmes;i++) if (residual[i]==99999) residual[i]=mask; /* for the conversion sonel -> hawaii */

  delay=0;
  i=1;
/*   res[0]=residual[0]; */

  i=1;
  for(j=1;j<nmes;j++) {
    if(gap[i]==1) res[j+delay+firsthour]=residual[i];
    else {
      for(k=j;k<(j+gap[i]-1);k++) res[k+delay+firsthour]=mask;
      res[k+delay+firsthour]=residual[i];
      delay+=gap[i]-1;
      }
    i++;
    }

  skip =julian_day(firstdate.month,firstdate.day,firstdate.year)-julian_day(1,1,1950);/* 1st date in timezone = GMT+...*/
  ndays=julian_day(firstdate.month,firstdate.day,firstdate.year)-julian_day(1,1,firstdate.year);/*in days since 1/1/year */

  if (firstdate.year==dum.year) {
    delay=0;   /* no shift with GMT */
    skip=julian_day(1,1,dum.year)-julian_day(1,1,1950);/* skip between cnes ref and 1/1/year (if serie start after 1/1) */
    ndays=julian_day(dum.month,dum.day,dum.year)-julian_day(1,1,firstdate.year);/*in days since 1/1/year */
    res[firsthour]=residual[0];
    }
  else
    if (firstdate.year>dum.year) delay=(int)(dum.second-firstdate.second)/3600.0-24;  /*GMT+delay*/
    else {
      delay=(int) 24-(firstdate.second-dum.second)/3600.0;  /*GMT+delay, delay<0*/
      for(i=0;i<delay;i++) res[i]=mask;
      res[delay]=residual[0];
      }

  for(i=0;i<n;i++) t[i]=(skip+ndays)*24.0+i;

  /*---------------------------------------------*/
  /* filling of the beginning of the first year  */
  /*---------------------------------------------*/

  if (ndays!=0) hawaii_write_header(strcode,firstdate.year,lat,lon,hawaii);
  for(i=0;i<12;i++) data[i]=mask;
  for(i=0;i<2*ndays;i++) {
    t0=skip*24+i*12;
    year=hawaii_write_data(strcode,t0,data,mask,hawaii,lat,lon,firstdate.year);
    }

  /*--------------------*/
  /*    writing data    */
  /*--------------------*/

  t0=(skip-ndays)*24+i*12+firstdate.second/3600;
  i=0;
  while(t0<(t1-12)&&(i+12<=n)) { /* modif:'i+12<n' before, LR, 30/04/04 */
    t0=t[i];
    for(j=0,k=0;j<12;j++,k++) data[k]=res[i-delay+j]; /* +delay*/
    year=hawaii_write_data(strcode,t0,data,mask,hawaii,lat,lon,year);
    i+=12;
    }

  free(res);
  free(gap);
  free(t);
}


//extern int load_hawaii(char *filename, double **elevation, double **time, double *mask, int *n, char *outfile);
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double mask=999.9;
  double t0,origine,t1=0,t2=1.e+10;
  double *h,*t;

  double start,finish,time,epsilon;

  int i,n,status;
  int filter,test=0;
  char *data=NULL, *input=NULL, *keyword, *output=NULL, *translation=NULL, *format=NULL;
  char duplicated[1024],*s;
  bool add_line=False;
  pressure_station *sample;
  statistic_t stat[3],sdum;
  int  year,month,day,hour;
  int nvalues;
  int start_year=1900,last_year=2100,start_month=1,last_month=12;
  mooring_t mooring;
  int decimating=0,duplicate=0,truncate=0;

  fct_echo( argc, argv);

  fprintf(stderr,"%s  -starting computation *********\n",argv[0]);

  n=1;
  while (n < argc) {
    if(strcmp(argv[n],"-start") ==0) {
      status=sscanf(argv[n+1],"%d/%d",&start_month,&start_year);
      n++;
      n++;
      decimating=1;
      continue;
      }
    if(strcmp(argv[n],"-end") == 0) {
      status=sscanf(argv[n+1],"%d/%d",&last_month,&last_year);
      n++;
      n++;
      decimating=2;
      continue;
      }
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 'd' :
          duplicate=1;
          n++;
          break;

        case 't' :
          truncate=1;
          n++;
          break;

        case 'c' :
          test=1;
          n++;
          break;

        case 'f' :
          format= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
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
        if(input==0) {
          input= strdup(argv[n]);
          n++;
          }
        else {
          __OUT_BASE_LINE__("input file already given: %s\n",input);
          exit(-1);
          }
        break;
      }
      free(keyword);
    }

  if(format==0) format=strdup("obs");

/* *-----------------------------------------------------------------------------
  read the input file*/
  if(strcmp(format,"obs")==0) {
    if(output==NULL) {
      output=strdup(input);
      output[strlen(output)-4]=0;
      strcat(output,".obs");
      }

    date_t firstdate,start,finish,test;
    double t1;

/* *-----------------------------------------------------------------------------
    some test on time functions*/
//     test.year=1916;
//     test.month=12;
//     test.day=31;
//     test.second=3600;
//
//     test.year=1917;
//     test.month=1;
//     test.day=1;
//     test.second=3600.;
//     t1=cnes_time(test,'d');
//     getcnesdate( t1,&firstdate,'d');
//
//     test.year=1917;
//     test.month=1;
//     test.day=1;
//     test.second=0.;
//     t1=cnes_time(test,'d');
//     getcnesdate(t1,&firstdate,'d');

    nvalues= hawaii_getlength(input,&firstdate,&finish,&mooring);
    t1=cnes_time(finish,'h')+11.0;
    t=new double[nvalues];
    h=new double[nvalues];
    FILE *hawaii=fopen(input,"r");
    nvalues=read_hawaii(hawaii,&(mooring.lat),&(mooring.lon),t,h, firstdate,nvalues);
    fclose(hawaii);

    int numstation;
    //char *stationname="DUM";//warning and unused
    if(duplicate==1) {
      strcpy(duplicated,input);
      s=strstr(duplicated,".dat");
      *s=0;
      strcat(duplicated,".duplicated.dat\0");

      hawaii=fopen(duplicated,"w");
      save_hawaii(hawaii,mooring.lat,mooring.lon,h,t,mooring.code,mooring.name, firstdate, nvalues, t1);
      fclose(hawaii);
      }
    if(decimating==1) {
      start.year=start_year;
      }
    else {
      start.year=1800;
      }
    if(decimating==2) {
      start.year =start_year;
      finish.year=last_year;
      }
    else {
      start.year=1800;
      finish.year=2100;
      }
    start.month=1;
    start.day=1;
    start.second=0;
    status=mgr_loadGLOSS_local(input, start, finish, &mooring, &h, &t, &mask, &nvalues, output);
    firstdate=poctime_getdatecnes(t[0],'d');

    if(truncate==1) {
      strcpy(duplicated,input);
      s=strstr(duplicated,".dat");
      *s=0;
      strcat(duplicated,".truncated.dat\0");

      hawaii=fopen(duplicated,"w");
      poctime_convert(t,nvalues,'d','h');
      for (i=0; i< nvalues; i++) {
        if(h[i]!=mask) {
          h[i]=h[i]*1000.0;
          }
        else h[i]=mask;
        }
      save_hawaii(hawaii,mooring.lat,mooring.lon,h,t,mooring.code,mooring.name, firstdate, nvalues, t1);
      fclose(hawaii);
      }
    }
  else if (strcmp(format,"matroos")==0) {
    if(output==NULL) {
      output=strdup(input);
      output[strlen(output)-4]=0;
      strcat(output,".slv");
      }
    status=hawaii_load(input, &mooring, &h, &t, &mask, &nvalues, (char *) 0);
    mooring.type=0;
    status=matroos_save(output, mooring, h, t, mask, nvalues);
    }

  __ERR_BASE_LINE__("%s -computation complete ^^^^^^^^^\n",argv[0]);
  exit(0);
}
