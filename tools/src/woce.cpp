
/*******************************************************************************

  T-UGO tools, 2006-2016

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Thierry Letellier  LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief input/output for some ascii timeseries formats (MATROOS, HAWAII, senetosa, aktarus and maybe others).
*/
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "tools-define.h"
#include "poc-time.h"
#include "tides.h"

#ifndef VERBOSE
#define VERBOSE -1
#endif

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int matroos_load(const char *filename, double t1, double t2, double **h, double **t, double *mask, int *n, double *dt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  Lit la serie temporelle a etudier
 -----------------------------------------------------------------------
  input: filename

-----------------------------------------------------------------------
   output:
      H:   mesures
      T:   timing des mesures, in hours, from 1/1/1950 0 hours
      N:   nombre total de mesures
-----------------------------------------------------------------------*/
{
  FILE *in=NULL;
  int  k,nitems,nskip;
  int N;
  int year,month,day,hour,minute;
  double scale=1.,local=9.999;
  double dT,time,start,finish=-1.,previous,dum;
  double d1,arrondi;
  char a,line[256];

  dT=1.e+10;
  N=0;

  if ( !(in=fopen(filename,"r"))) return(-1);

  nskip=-1;
  do {
    fgets(line,256,in);
    nskip++;
    }  while ((line[0] == '#') && !feof(in));
  rewind(in);

  for(k=0;k<nskip;k++) fgets(line,256,in);

  time=-1.e+10;
  while(time<t1) {
//    fgets(line,256,in);
    nitems=fscanf(in,"%4d%2d%2d%2d%2d %lf",&year,&month,&day,&hour,&minute,&d1);
    time=julian_day(month,day,year)+hour/24.+minute/24./60.-CNES0jd;
    if (nitems != 6) goto error;
    do { a=fgetc(in);   }  while ((a != '\n') && !feof(in));
    }

  start=time;
  if (start>=t2) goto error;

  previous=start;
  while(!feof(in) && (time<t2)) {
    nitems=fscanf(in,"%4d%2d%2d%2d%2d %lf",&year,&month,&day,&hour,&minute,&d1);
    time=julian_day(month,day,year)+hour/24.+minute/24./60.-CNES0jd;
    do {a=fgetc(in);}  while ((a != '\n') && !feof(in));
    if(time > previous) {
/*------------------------------------------------------------------------------
      minute precison */
      dum=floor((time-previous)*24.*60.0+0.5)/24.0/60.0;
      if ((dT > dum) && (dum !=0)) {
        dT=dum;
        }
      }
    previous=time;
    finish=max(time,finish);
    }

  arrondi=5.; /*arrondi en minute*/
  dT=floor(dT*24.*60.0/arrondi+0.5)*arrondi/24.0/60.0;

  N=(int)( floor((finish-start)/dT+1+0.5) );

  printf("#observations: %d, sampling: %g s,  %s to  %s \n",N,dT*24*3600,sgetcnesdate(start*24.),sgetcnesdate(time*24.));

  exitIfNull(
    *h=(double *) malloc(N*sizeof(double))
    );
  exitIfNull(
    *t=(double *) malloc(N*sizeof(double))
    );

  for(k=0;k<N;k++) (*t)[k]=start+(double) k*dT;
  for(k=0;k<N;k++) (*h)[k]=*mask;

  rewind(in);

  for(k=0;k<nskip;k++) fgets(line,256,in);

  while(!feof(in)) {
    nitems=fscanf(in,"%4d%2d%2d%2d%2d %lf",&year,&month,&day,&hour,&minute,&d1);
    time=julian_day(month,day,year)+hour/24.+minute/24./60.-CNES0jd;
    if(time<t1) continue;
    if(time>t2) break;
    k=(int)( floor((time-start)/dT+0.5) );
    if(d1 != local) {
      (*h)[k]=d1/scale;
      }
    else {
      (*h)[k]=*mask;
      }
    }
  fclose(in);

  *n=N;
  *dt=dT;
  return(0);

error:
  fclose(in);
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

mooring_t hawaii_decode(char *line)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//001POHNPEI 1975  LAT=06 59.2N  LONG=158 14.6E  TIMEZONE=GMT
//002BETIO   1992  LAT=01 21.7N  LONG=172 55.8E  TIMEZONE=GMT
  int   n;
  int   lon, lat, year;
  float lonmin,latmin;
  char  dum[256],name[256],north,east;
  mooring_t local;

  sscanf(line,"%3d",&n);

  sscanf(&(line[3]),"%s",name);
  local.name=strdup(name);
  sscanf(&(line[10]),"%d",&year);
  sscanf(&(line[17]),"%4s",dum);
  sscanf(&(line[21]),"%d",&lat);
  sscanf(&(line[24]),"%f",&latmin);
  sscanf(&(line[28]),"%c",&north);
  switch (north) {
    case 'N':
      local.lat=lat+latmin/60.0;
      break;
    case 'S':
      local.lat=-(lat+latmin/60.0);
      break;
    }

  sscanf(&(line[31]),"%5s",dum);
  sscanf(&(line[36]),"%d",&lon);
  sscanf(&(line[40]),"%f",&lonmin);
  sscanf(&(line[44]),"%c",&east);
  switch (east) {
    case 'E':
      local.lon=lon+lonmin/60.0;
      break;
    case 'W':
      local.lon=-(lon+lonmin/60.0);
      break;
    }

  return(local);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int hawaii_load(const char *filename, mooring_t *mooring, double **elevation, double **time, double *mask, int *n, char *outfile)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *-----------------------------------------------------------------------
  Lit la serie temporelle a etudier

  Le temps est en heure compte a partir de l'origine des temps
  t0= 1er janvier 1950 a 0 heure

  Les mesures sont lues en MM "format HAWAII"  et converties en METRES
-----------------------------------------------------------------------
  input: filename
 
-----------------------------------------------------------------------
   output:

      H:   mesures
      T:   timing des mesures, in hours, from 1/1/1950 0 hours
      N:   nombre total de mesures
      Code , Longitude et latitude de la station
      Timezone
-----------------------------------------------------------------------*/
{
  FILE *in=NULL,*out=NULL;
  char line[256],*info=0;
  int  year_p,year,month,day,ind;
  char name[11];
  float z[12],w;
  int i,j,count,N,fdepth,nitems;
  date_t start;
  double *buf=NULL,*lf=NULL,*mf=NULL,average;
  double t0;
  double *h=NULL, *t=NULL;

  in=fopen(filename,"r");
  if(in==0) TRAP_ERR_EXIT(errno,"fopen(\"%s\",\"r\") error (%d %s)\n",filename,errno,strerror(errno));

  year_p=-1;
  year=0;
  N=0;
  while(!feof(in)) {
//    fscanf(in,"%11c%4d %s %s %s %s %s\n",name,&year,s1,s2,s3,s4,s5);
    fgets(line,256,in);
    if(info==0) {
      info=strdup(line);
      }
    nitems=sscanf(line,"%11c%4d",name,&year);
    if(year != year_p) {
      if(year<1950) continue;
      year_p=year;
      printf("year: %d \n",year);
      }
    do {
      if (feof(in)) break;
      nitems=fscanf(in,"%11c%4d",name,&year);
      if(year != year_p) {
        year_p=year;
        printf("year: %d \n",year);
        continue;
        }
      fgets(line,3,in);
      month=atoi(line);
      fgets(line,3,in);
      day=atoi(line);
      fgets(line,2,in);
      ind=atoi(line);
      /*      printf("%s %4d %2d %2d %1d \n",name,year,month,day,ind);*/
      if(start.year == 0) {
        start.year=year;
        start.month=month;
        start.day=day;
        start.second=(ind-1)*12*3600.0;
        }
      for (i=0; i<12;i++) fscanf(in,"%f",&z[i]);
      fscanf(in,"\n");
      N=N+12;
      } while (month != 12 || day != 31 || ind != 2);
    }
  rewind(in);

  printf("Number of observations: %d (~ %f years)\n",N,N/24./30.5/12.);

  buf=new double[N];

  N=0;
  while(!feof(in)) {
//    fscanf(in,"%11c%4d %s %s %s %s %s\n",name,&year,s1,s2,s3,s4,s5);
    fgets(line,256,in);
    sscanf(line,"%11c%4d",name,&year);
    if(year != year_p) {
      if(year<1950) continue;
      year_p=year;
      printf("year: %d \n",year);
      }
    do {
      if (feof(in)) break;
      fscanf(in,"%11c%4d",name,&year);
      if(year != year_p) {
        year_p=year;
        printf("year: %d \n",year);
        continue;
        }
      fgets(line,3,in);
      month=atoi(line);
      fgets(line,3,in);
      day=atoi(line);
      fgets(line,2,in);
      ind=atoi(line);
      /*      printf("%s %4d %2d %2d %1d \n",name,year,month,day,ind);*/
      for (i=0; i<12;i++) fscanf(in,"%lf",&buf[i+N]);
      fscanf(in,"\n");
      N=N+12;
      } while (month != 12 || day != 31 || ind != 2);
    }
  fclose(in);

/*
  t0=(julian_day(start)-CNES0jd)*24*3600.+start.second;
*/

  h=new double[N];
  *mask=9999.0;

/* *----------------------------------------------------------------------------
  convert to meters*/
  for (i=0; i< N; i++) {
    if(buf[i]!=*mask) {
      h[i]=buf[i]/1000.0;
      }
    else h[i]=*mask;
    }

  average=0;
  count=0;
  for (i=0; i< N; i++) {
    if(h[i]!=*mask) {
      average+=h[i];
      count++;
      }
    }
  average /= count;

  fdepth=48;
  lf=new double[N];
  for (i=0; i< N; i++) {
    lf[i]=0;
    w=0;
    for (j=max(0,i-fdepth); j< min(N,i+fdepth); j++) {
      if(h[j]!=*mask) {
        lf[i]+= h[j];
        w=w+1;
        }
      }
    if (w != 0) lf[i]=lf[i]/w-average;
    else lf[i]=0.;
    }

  fdepth=6;
  mf=new double[N];
  for (i=0; i< N; i++) {
    mf[i]=0;
    w=0;
    for (j=max(0,i-fdepth); j< min(N,i+fdepth); j++) {
      if(h[j]!=*mask) {
        mf[i]+= h[j];
        w=w+1;
        }
      }
    if (w != 0) mf[i]=mf[i]/w-average;
    else mf[i]=0.;
    }

  t=new double[N];
  t0=(julian_day(start)-CNES0jd)+start.second/24./3600.;
  for (i=0; i< N; i++) {
    t[i]=t0+(double)i/24.;
    }

  *mooring=hawaii_decode(info);


  if(outfile!=0) {
    out=fopen(outfile,"w");
    if(out==0) TRAP_ERR_EXIT(-1,"unable to open file: %s \n",outfile);
    for (i=0; i< N; i++) {
      if(h[i]!=*mask) {
//        fprintf(out,"%6d %10.3f %10.3f %10.3f %10.3f %10.3f\n",
//              i,t[i],h[i]-average,lf[i],h[i]-average-lf[i],mf[i]-lf[i]);
        fprintf(out,"%10.3f %10.3f %10.3f %10.3f %10.3f\n",
              t[i],h[i]-average,lf[i],h[i]-average-lf[i],mf[i]-lf[i]);
        }
      }
    fclose(out);
    }

  delete[] mf;
  delete[] lf;
  delete[] buf;

  *n=N;
  *elevation=h;
  *time=t;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int senetosa_load(const char *filename, double t1, double t2, double **h, double **t, double mask, int *n, double *dt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  Lit la serie temporelle a etudier
 -----------------------------------------------------------------------
  input: filename

-----------------------------------------------------------------------
   output:
        
      H:   mesures
      T:   timing des mesures, in hours, from 1/1/1950 0 hours
      N:   nombre total de mesures
      Code , Longitude et latitude de la station
      Timezone
-----------------------------------------------------------------------*/
{
  FILE *in=NULL;
  int  k,nitems;
  int N;
  const double scale=1.,local=9.999;
  double dT,time,start,finish=-1.,previous,dum;
  double d1,arrondi;
  char a,line[256];

  dT=1.e+10;
  N=0;

  in=fopen(filename,"r");
  
  if (in==0) return(errno);

  time=-1.e+10;
  while(time<t1) {
    nitems=fscanf(in,"%lf %f",&time,&d1);
    if (nitems != 2) {
      goto error;
      }
    do { a=fgetc(in);   }  while ((a != '\n') && !feof(in));
    }

  start=time;
  if (start>=t2) goto error;

  previous=start;
  while(!feof(in) && (time<t2)) {
    fscanf(in,"%lf %f",&time,&d1);
    do {a=fgetc(in);}  while ((a != '\n') && !feof(in));
    if(time > previous) {
/*------------------------------------------------------------------------------
      minute precison */
      dum=floor((time-previous)*24.*60.0+0.5)/24.0/60.0;
      if ((dT > dum) && (dum !=0)) {
        dT=dum;
        }
      }
    previous=time;
    finish=max(time,finish);
    }

  arrondi=5.; /*arrondi en minute*/
  dT=floor(dT*24.*60.0/arrondi+0.5)*arrondi/24.0/60.0;

/*------------------------------------------------------------------------------
  DANGER, PROVISOIRE */
//  dT=1./24.; //LR, MODIFICATION: 12/12/2006

  N=(int)( floor((finish-start)/dT+1+0.5) );

  printf("#observations: %d, sampling: %g s,  %s to  %s \n",N,dT*24*3600,sgetcnesdate(start*24.),sgetcnesdate(time*24.));

  exitIfNull(
    *h=(double *) malloc(N*sizeof(double))
    );
  exitIfNull(
    *t=(double *) malloc(N*sizeof(double))
    );
  for(k=0;k<N;k++) (*t)[k]=start+(double) k*dT;
  for(k=0;k<N;k++) (*h)[k]=mask;

  rewind(in);

  while(!feof(in)) {
    fgets(line,sizeof(line),in);
    d1=local;
    nitems=sscanf(line,"%lf %lf",&time,&d1);
/*    if (nitems != 2)  */
/*       goto error;    */
    if(time<t1) continue;
    if(time>t2) break;
    k=(int)( floor((time-start)/dT+0.5) );
    if(d1 != local) {
      (*h)[k]=d1/scale;
      }
    else {
      (*h)[k]=mask;
      }
/*     (*t)[k]=time; */
    }
  fclose(in);

  *n=N;
  *dt=dT;
  return(0);

  error:
  fclose(in);
  return(-1);
}

/*-------------------------------------------------------------------------------*/

int skipHeader(FILE *fic_data)

/*-------------------------------------------------------------------------------*/
{

  // copy of header() from woce03.cpp
  // woce03.cpp has some aktarus/tools incompatibilities
  
  int nXpoints;
  char line[300];

  do  fgets(line,sizeof(line),fic_data);
  while (strncmp(line,"#-- HEADER END ---------------------------------------",20) != 0 );

  fgets(line,sizeof(line),fic_data);
  fgets(line,sizeof(line),fic_data);
  sscanf(line, "%d", &nXpoints);
  return(nXpoints);
}

/*--------------------------------------------------------------------------------*/

int aktarus_load_serie(double *time, double *sealevel, double *lon, double *lat,FILE *fic_data,int heure)

/*--------------------------------------------------------------------------------*/
{
  char a;
  int i, Mes = -1;
  char line[300];

  fgets(line,sizeof(line),fic_data); /*#--------------------- */
  fgets(line,sizeof(line),fic_data); /* Pt  :         */

  do {
    a=fgetc(fic_data);
    } while ( a!= ':');
  fscanf(fic_data, "%lf",lon);
  if (*lon<0) *lon+=360;
  do {a=fgetc(fic_data);} while ( a!= '\n');

  do {a=fgetc(fic_data);} while ( a!= ':');
  fscanf(fic_data, "%lf",lat);
  do {a=fgetc(fic_data);} while ( a!= '\n');

  do {a=fgetc(fic_data);} while ( a!= ':');
  fscanf(fic_data, "%d",&Mes);
  do {a=fgetc(fic_data);} while ( a!= '\n');

  for (i=0; i<Mes;i++) {
    fgets(line,sizeof(line),fic_data);
    sscanf(line, " %lf %lf", &time[i],&sealevel[i]);
    if (heure==0)       time[i]/=3600;   /*temps en secondes -> temps en heures*/
    else if (heure==2)  time[i]*=24;     /*temps jours --> temps en heures*/
    }

  return(Mes);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  serie_t aktarus_reduce(serie_t & source,date_t start, date_t end)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  serie_t reduce;
  data_t *tmp = NULL;

  double t1;
  double t2;
  
  if(isad(start))
    t1 = julian_day(start)-CNES0jd;
  else
    t1 = -INFINITY;
  
  if(isad(end))
    t2 = julian_day(end)-CNES0jd;
  else
    t2 = +INFINITY;
  
  tmp = new data_t[source.count]; //over alloc
  reduce.count=0;
  for(size_t i = 0; i < source.count; i++){
    if(source.data[i].values[0] != source.mask)
      if((source.data[i].time >= t1) && (source.data[i].time <= t2)){
        tmp[reduce.count++] = source.data[i];
        }
    }

  if(source.data!=0) {
    delete[] source.data;
    source.data = 0;
    }
  if( reduce.count!=0) {
    reduce.data = new data_t[reduce.count];
    }
  else {
    reduce.data = 0;
    }
  
  for(size_t i = 0; i < reduce.count; i++){
    reduce.data[i] = tmp[i];
    }
  reduce.mask = source.mask;
  
  delete[] tmp;

#if VERBOSE > 0
  printf("reduced time series starts %s, finishes %s \n",sgetcnesdate(reduce.data[0].time*24.),sgetcnesdate(reduce.data[reduce.count-1].time*24.));
#endif
  return(reduce);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int load_serie(const char *filename, double **time, double **z, char unit)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, n,ndata;
  char line[300];
  FILE *in=NULL;

  in=fopen(filename,"r");
  if(in==0) TRAP_ERR_EXIT(errno,"fopen(\"%s\",\"r\") error (%d %s)\n",filename,errno,strerror(errno));

  n=0;
  ndata=0;
  while(!feof(in)) {
    fgets(line,sizeof(line),in);
    if(line[0]!='#') ndata++;
    n++;
    }

  n--;
  ndata--;

  rewind(in);

  *time=new double[ndata];
  *z=   new double[ndata];
  
  ndata=0;
  for (i=0; i<n;i++) {
    fgets(line,sizeof(line),in);
    if(line[0]!='#') {
      sscanf(line, " %lf %lf", &(*time)[ndata],&(*z)[ndata]);
      ndata++;
      }
    }
  fclose(in);
  return(ndata);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int load_serie1D(const char *filename, double **time, double **z, char unit, mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  char *s;
  int i, n,nitems,ndata;
  char line[300];
  FILE *in=NULL;

  in=fopen(filename,"r");
  if(in==0) TRAP_ERR_EXIT(errno,"fopen(\"%s\",\"r\") error (%d %s)\n",filename,errno,strerror(errno));

  n=0;
  ndata=0;
  while(!feof(in)) {
    s=fgets(line,sizeof(line),in);
    if(s==0) break;
    if(line[0]!='#') ndata++;
    n++;
    }

  n--;
  ndata--;

  rewind(in);

  *time=new double[ndata];
  *z=   new double[ndata];
  
  ndata=0;
  for (i=0; i<n;i++) {
    fgets(line,sizeof(line),in);
    if(line[0]!='#') {
      nitems=sscanf(line, " %lf %lf", &(*time)[ndata],&(*z)[ndata]);
      if(nitems!=2) {
        printf("erroneous line (%d), shoulbd be: time value (2 items)\n",i);
        fclose(in);
        return(-1);
        }
      ndata++;
      }
    else {
      char key[32];
      mooring->name=new char[256];
      nitems=sscanf(line, "%s %s %lf %lf %f", key, mooring->name, &(mooring->lon),&(mooring->lat),&(mooring->depth));
      if(nitems!=5) {
        printf("erroneous header, shoulbd be: #X...X Y...Y longitude latitude depth (5 items)\n");
        fclose(in);
        return(-1);
        }
      }
    }
  fclose(in);
  return(ndata);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int load_serie2D(const char *filename, double **time, double ***z, char unit, mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, n,ndata;
  char line[300];
  FILE *in=NULL;
  char key[32];

  in=fopen(filename,"r");
  if(in==0) TRAP_ERR_EXIT(errno,"fopen(\"%s\",\"r\") error (%d %s)\n",filename,errno,strerror(errno));

  n=0;
  ndata=0;
  while(!feof(in)) {
    fgets(line,sizeof(line),in);
    if(line[0]!='#') ndata++;
    n++;
    }

  n--;
  ndata--;

  rewind(in);

  *time=   new double[ndata];
  (*z)[0]=   new double[ndata];
  (*z)[1]=   new double[ndata];
  
  ndata=0;
  for (i=0; i<n;i++) {
    fgets(line,sizeof(line),in);
    if(line[0]!='#') {
      sscanf(line, " %lf %lf %lf", &(*time)[ndata],&(((*z)[0])[ndata]),&(((*z)[1])[ndata]));
      ndata++;
      }
    else {
      mooring->name=new char[256];
      sscanf(line, "%s %s %lf %lf %lf", key, mooring->name, &(mooring->lon),&(mooring->lat),&(mooring->depth));
      }
    }
  fclose(in);
  return(ndata);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int sample_load(const char *filename, double **h, double **u, double **v, double **p,double **t, double *mask, int *n, double *dt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  Lit la serie temporelle a etudier
 -----------------------------------------------------------------------
  input: filename
-----------------------------------------------------------------------
  output:
      H:   mesures
      T:   timing des mesures, in hours, from 1/1/1950 0 hours
      N:   nombre total de mesures
      Code , Longitude et latitude de la station
      Timezone
-----------------------------------------------------------------------*/
{
  FILE *in=NULL;
  int  k;
  int i,N=0,shift=0,nitems;
  double d1,d2,d3,d4,arrondi;
  double dT,time,start,finish=0.,previous,dum;
  char a,line[256];
  
  *mask=999.9;

  dT=1.e+10;

  if ( !(in=fopen(filename,"r"))) return(0);

#define _COLDSTART
#ifdef  _COLDSTART
  for (i=0;i<shift*24; i++) do { a=fgetc(in);   }  while ((a != '\n') && !feof(in));
#endif

  nitems=fscanf(in,"%lf %lf",&start,&d1);
  if(nitems!=2) {
    goto error;
    }

  do { a=fgetc(in);   }  while ((a != '\n') && !feof(in));

  previous=start;
  while(!feof(in)) {
    nitems=fscanf(in,"%lf %f",&time,&d1);
    if(nitems!=2) {
      break;
      }
    do { a=fgetc(in);   }  while ((a != '\n') && !feof(in));
    if(time > previous) {
/*------------------------------------------------------------------------------
      minute precison */
      dum=floor((time-previous)*24.*60.0+0.5)/24.0/60.0;
      if ((dT > dum) && (dum !=0)) {
        dT=dum;
        }
      }
    previous=time;
    finish=max(time,finish);
    }

  arrondi=5.; /*arrondi en minute*/
  dT=floor(dT*24.*60.0/arrondi+0.5)*arrondi/24.0/60.0;

  N=(int)( floor((finish-start)/dT+1+0.5) );
  printf("Number of observations: %d\n",N);

  exitIfNull(
    *h=new double[N]
    );
  exitIfNull(
    *u=new double[N]
    );
  exitIfNull(
    *v=new double[N]
    );
  exitIfNull(
    *p=new double[N]
    );
  exitIfNull(
    *t=new double[N]
    );

  for(k=0;k<N;k++) (*t)[k]=start+(double) k*dT;
  for(k=0;k<N;k++) (*h)[k]=*mask;
  for(k=0;k<N;k++) (*u)[k]=*mask;
  for(k=0;k<N;k++) (*v)[k]=*mask;
  for(k=0;k<N;k++) (*p)[k]=*mask;

  rewind(in);
#ifdef _COLDSTART
  for (i=0;i<shift*24; i++) do { a=fgetc(in);   }  while ((a != '\n') && !feof(in));
#endif

  while(!feof(in)) {
    fgets(line,sizeof(line),in);
    nitems=sscanf(line,"%lf %lf %lf %lf %lf",&time,&d1,&d2,&d3,&d4);
    if(nitems!=5) {
      break;
      }
    k=(int)( floor((time-start)/dT+0.5) );
    (*h)[k]=d1;
    (*u)[k]=d2;
    (*v)[k]=d3;
    (*p)[k]=d4;
    (*t)[k]=time;
    }
  fclose(in);

  *n=N;
  *dt=dT;
  return(1);

error:
  *n=N;
  *dt=-1.;
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int sample_load(const char *filename, tseries_t & serie, double *dt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=sample_load(filename, &serie.x[0], &serie.x[1], &serie.x[2], &serie.x[3], &serie.t, &serie.mask, &serie.n, dt);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int sample_save(const char *filename, double *h, double *u, double *v, double *p, double mask, double *t, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  double time;
  bool add_line=True;
  FILE *out=NULL;

  printf("create %s\n",filename);

  out=fopen(filename,"w");
  if(out==0) TRAP_ERR_EXIT(errno,"fopen(\"%s\",\"w\") error (%d %s)\n",filename,errno,strerror(errno));
  for (i=0; i< n; i++) {
    if(h[i] != mask) {
      time=t[i];
      fprintf(out,"%12.6lf %6.3lf \n",time,h[i],u[i],v[i],p[i]);
      add_line=True;
      }
    else {
      if(add_line) fprintf(out,"\n");
      add_line=False;
      }
    }
  fclose(out);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_timeserie(const char *filename, double *h, double mask, double *t, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  double time;
  bool add_line=True;
  FILE *out=NULL;

//  printf("create %s\n",filename);

  out=fopen(filename,"w");
  if(out==0) TRAP_ERR_EXIT(errno,"fopen(\"%s\",\"w\") error (%d %s)\n",filename,errno,strerror(errno));

  for (i=0; i< n; i++) {
    if(h[i] != mask) {
      time=t[i];
      fprintf(out,"%12.6lf %6.3lf \n",time,h[i]);
      add_line=True;
      }
    else {
      if(add_line) fprintf(out,"\n");
      add_line=False;
      }
    }
  fclose(out);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matroos_save(const char *filename, mooring_t mooring, double *h, double *t, double mask, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
-----------------------------------------------------------------------*/
{
  FILE *out=NULL;
  int  k,nitems;
  int hour,minute;
  date_t date;
  const char *nature=NULL;

  if ( !(out=fopen(filename,"w"))) return(-1);

  switch (mooring.type) {
    case 0:
      nature="observed";
      break;
    case 1:
      nature="model";
      break;
    case -1:
    default:
      nature="unknown";
      break;
   }

  fprintf(out,"#------------------------------------------------------\n"//
              "# Timeseries retrieved from the MATROOS series database\n"//
              "# Created at Mon Dec 21 12:28:35 CET 2009\n"//
              "#------------------------------------------------------\n"//
              "# Location    : %s\n"//
              "# Position    : (%f,%f)\n"//
              "# Source      : %s\n"//
              "# Unit        : waterlevel_surge\n"//
              "# Analyse time: 197001010000\n"//
              "# Timezone    : GMT\n"//
              "#------------------------------------------------------\n",
               mooring.name,mooring.lon,mooring.lat,nature);

  for(k=0;k<n;k++) {
    date=poctime_getdatecnes(t[k],'d');
    hour=date.second/3600.;
    minute=(date.second-3600.*hour)/60.;
    nitems=fprintf(out,"%4d%2.2d%2.2d%2.2d%2.2d %9.4lf\n",date.year,date.month,date.day,hour,minute,h[k]);
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int timeserie_save(const char *filename, mooring_t mooring, double *h, double mask, double *t, int n, int fmt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  switch (fmt) {
    case TIMESRIES_FORMAT_MATROOS:
      status= matroos_save(filename, mooring, h, t, mask, n);
      break;

    case TIMESRIES_FORMAT_CLASSIC:
      status= save_timeserie(filename, h, mask, t, n);
      break;
    }
  return(status);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
// int timeserie_save(char *filename, mooring_t mooring, serie_t serie, char *fmt)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   switch (fmt) {
//     case TIMESRIES_FORMAT_MATROOS:
//       status= matroos_save(filename, mooring, h, t, mask, n);
//       break;
//
//     case TIMESRIES_FORMAT_CLASSIC:
//       status= save_timeserie(filename, h, mask, t, n);
//       break;
//     }
//   return(status);
// }

