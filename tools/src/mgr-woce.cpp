
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file
tide gauges stuff : woce ?
*/
/*----------------------------------------------------------------------------*/

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "poc-netcdf-data.hpp"
#include "filter.h"
#include "tides.h"
#include "mgr.h"
#include "functions.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int hawaii_decode_position_p(char *line, mooring_t *local, int *year)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//001POHNPEI 1975  LAT=06 59.2N  LONG=158 14.6E  TIMEZONE=GMT
//002BETIO   1992  LAT=01 21.7N  LONG=172 55.8E  TIMEZONE=GMT
  int   n;
  int   lon, lat, year_;
  if(year==0)
    year=&year_;
  mooring_t local_;
  if(local==0)
    local=&local_;
  float lonmin,latmin;
  const int nameLength=256;
  char  name[nameLength],north,east;
  
  aset(name,nameLength,'\0');
  
#if 0
  char  dum[256];
  
  sscanf(line,"%3d",&n);

  local->code=n;

  sscanf(&(line[3]),"%8c",name);

  name[8]=0;
  local->name=strdup(name);

  sscanf(&(line[10]),"%d",year);
  sscanf(&(line[17]),"%4s",dum);
  sscanf(&(line[21]),"%d",&lat);
  sscanf(&(line[24]),"%f",&latmin);
  sscanf(&(line[28]),"%c",&north);
  switch (north) {
    case 'N':
      local->lat=lat+latmin/60.0;
      break;
    case 'S':
      local->lat=-(lat+latmin/60.0);
      break;
    default:
      TRAP_ERR_EXIT(-1, "format error\n");
    }

  sscanf(&(line[31]),"%5s",dum);
  sscanf(&(line[36]),"%d",&lon);
  sscanf(&(line[40]),"%f",&lonmin);
  sscanf(&(line[44]),"%c",&east);
  switch (east) {
    case 'E':
      local->lon=lon+lonmin/60.0;
      break;
    case 'W':
      local->lon=-(lon+lonmin/60.0);
      break;
    default:
      TRAP_ERR_EXIT(-1, "format error\n");
    }
#else
  char lonS[nameLength],latS[nameLength];
  
  aset(lonS,nameLength,'\0');
  aset(latS,nameLength,'\0');
  
  n=sscanf(line,"%03d%8c%4d %*4c%d %4c%c %*5c%d %4c%c",
    &local->code,name,year,
    &lat,latS,&north,
    &lon,lonS,&east);
  
  if(n>=3 and local==&local_)
    return 0;
  
  if(n==9){
    sscanf(latS,"%f",&latmin);
    sscanf(lonS,"%f",&lonmin);
    }
  else{
    /* JASL format */
    //001C Pohnpei-C          Fd. St. Micronesia  2014 06587N 158118E 0000 1 00000R MM
    //002D Tarawa-D,Betio     Rep. of Kiribati    2014 01218N 172558E 0000 1 00000R MM
    
    aset(lonS,nameLength,'\0');
    aset(latS,nameLength,'\0');
    
    n=sscanf(line,"%03d%*2c%19c%*20c%d %2c%3f%c%*c%3c%3f%c",
      &local->code,name,year,
      latS,&latmin,&north,
      lonS,&lonmin,&east);
    
    if(n!=9)
      TRAP_ERR_EXIT(-1, "format error: %d!=9 items scanned in:%s",n,line);
    
    sscanf(latS,"%d",&lat);
    sscanf(lonS,"%d",&lon);
    
    latmin*=0.1;
    lonmin*=0.1;
    }
  
  local->name=strdup(name);
  
  bool isError=false;
  
  switch (north) {
    case 'N':
      local->lat=lat+latmin/60.0;
      break;
    case 'S':
      local->lat=-(lat+latmin/60.0);
      break;
    default:
      isError=true;
    }
  
  switch (east) {
    case 'E':
      local->lon=lon+lonmin/60.0;
      break;
    case 'W':
      local->lon=-(lon+lonmin/60.0);
      break;
    default:
      isError=true;
    }
  
  if(isError)
    TRAP_ERR_EXIT(-1, "format error: north=`%c` and south=`%c`:%s\n",north,east,line);
#endif

  return(0);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

date_t hawaii_decode_data(FILE *in, char *line)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//001POHNPEI 1975  LAT=06 59.2N  LONG=158 14.6E  TIMEZONE=GMT
//002BETIO   1992  LAT=01 21.7N  LONG=172 55.8E  TIMEZONE=GMT
  int   k, nitems;
  char  name[256];
  mooring_t local;
  int  year,month,day,ind;
  char  s_year[8],s_month[8],s_day[8],s_code[8],s_ind[8];
  date_t start;

  for(k=0;k<8;k++) {
    s_year[k]=0;
    s_month[k]=0;
    s_day[k]=0;
    s_code[k]=0;
    s_ind[k]=0;
    }

  nitems=sscanf(line,"%3c%8c%4c%2c%2c%1c",s_code,name,s_year,s_month,s_day,s_ind);
  year=atoi(s_year);
  month=atoi(s_month);
  day=atoi(s_day);
  ind=atoi(s_ind);

  start.year=year;
  start.month=month;
  start.day=day;
  start.second=(ind-1)*12*3600.0;

  return(start);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_loadGLOSS_local(const char *filename, date_t first, date_t last, mooring_t *mooring, double **elevation, double **time, double *mask, int *n, const char *outfile)

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
  float w;
  int i,j,count,N,pN,fdepth;
  date_t start(0,-1,-1,0.);
  double *buf=NULL,*lf=NULL,*mf=NULL,average;
  double t0;
  double *h=NULL, *t=NULL;

  in=fopen(filename,"r");
  if (in == NULL) {
    TRAP_ERR_EXIT(errno,"error opening %s (%d %s)\n",filename,errno,strerror(errno));
    }
  fgets(line,256,in);
  info=strdup(line);
  rewind(in);

  hawaii_decode_position_p(info,mooring);
  if(elevation==0) return(0);
  
/* -----------------------------------------------------------------------------
  read number of data */
  N=0;
  while(!feof(in)) {
    fgets(line,256,in);
    N++;
    }
  rewind(in);

  N *= 12;
  printf("mgr_loadGLOSS_local: Number of lines: %d (~ %f years)\n",N,N/24./30.5/12.);

  buf=new double[N];

  N=0;
  year_p = -1;
  pN=0;
  while(!feof(in)) {
    fgets(line,256,in);
    hawaii_decode_position_p(line,0,&year);
    if(year != year_p) {
      if(year<first.year) continue;
      if(year>last.year)  break;
      printf("year: %d (%d obs)\n",year_p, N-pN);
      pN = N;
      year_p=year;
      }
    do {
      if (feof(in)) {
        printf("year: %d  (%d obs)\n",year, N-pN);
        break;
        }
      fscanf(in,"%11c%4d",name,&year);
      if(year != year_p) {
        printf("year: %d  (%d obs)\n",year_p, N-pN);
        pN = N;
        year_p=year;
        continue;
        }
      fgets(line,3,in);
      month=atoi(line);
      fgets(line,3,in);
      day=atoi(line);
      fgets(line,2,in);
      ind=atoi(line);
      if(start.year == 0) {
        start.year=year;
        start.month=month;
        start.day=day;
        start.second=(ind-1)*12*3600.0;
        }
      for (i=0; i<12;i++) fscanf(in,"%lf",&buf[i+N]);
      fscanf(in,"\n");
      N=N+12;
      } while (month != 12 || day != 31 || ind != 2);
    }
  fclose(in);

/* *----------------------------------------------------------------------------
  set cnes time*/
  h=new double[N];
  t=new double[N];
  *mask=9999.0;

  t0=(julian_day(start)-CNES0jd)+start.second/24./3600.;
  
  for (i=0; i< N; i++) {
    t[i]=t0+(double)i/24.;
    }

/* *----------------------------------------------------------------------------
  convert to meters*/
  for (i=0; i< N; i++) {
    if(buf[i]!=*mask) {
      h[i]=buf[i]/1000.0;
      }
    else h[i]=*mask;
    }

  *n=N;
  *elevation=h;
  *time=t;
    
  if(outfile==0) return(0);
  
/* *----------------------------------------------------------------------------
  some optional checks*/
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
  
  double cut=10*24.;
  double dt=1.0;
  int status;
  Loess1D(N,dt,cut,h,*mask,lf,&status);
  for (i=0; i< N; i++) {
    if(lf[i]!=*mask) lf[i]-=average;
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
  
  double *detided = tide_filters("demerliac",h, *mask, N);
  double *godin   = tide_filters("godin"    ,h, *mask, N);
  double *munk    = tide_filters("munk",     h, *mask, N);
  
  out=fopen(outfile,"w");
  if (out == NULL) {
    TRAP_ERR_EXIT(errno,"error opening %s (%d %s)\n",outfile,errno,strerror(errno));
    }
  for (i=0; i< N; i++) {
    if(h[i]!=*mask) {
      if(detided[i]==*mask) detided[i]=average;
      if(detided[i]!=*mask) detided[i]-=average;
      if(godin[i]==*mask) godin[i]=average;
      if(godin[i]!=*mask) godin[i]-=average;
      if(munk[i]==*mask) munk[i]=average;
      if(munk[i]!=*mask) munk[i]-=average;
      fprintf(out,"%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
            t[i],h[i]-average,lf[i],h[i]-average-lf[i],mf[i]-lf[i],detided[i],godin[i],munk[i]);
      }
    }
  fclose(out);

  delete[] detided;
  delete[] mf;
  delete[] lf;
  delete[] buf;
  free(info);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_loadGLOSS(const char *filename, tseries_t *serie, mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nvalues,status;

  const char *outfile=0;
  date_t first(1850,1,1,0.0);
  date_t last (2050,1,1,0.0);

  double *elevation;
  double *time;
  double mask;
  
  if(serie!=0) {
    status=mgr_loadGLOSS_local(filename, first, last, mooring, &elevation, &time, &mask, &nvalues, outfile);
    }
  else {
    status=mgr_loadGLOSS_local(filename, first, last, mooring, 0, 0, &mask, &nvalues, outfile);
    return(0);
    }
  
  serie->t=time;
  serie->x[0]=elevation;
  serie->mask=mask;
  serie->n=nvalues;
  serie->lat=mooring->lat;
  serie->lon=mooring->lon;
  
  return(nvalues);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_buggy_mgr_var(const char *filename, const poc_var_t & var, double *coord)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  size_t npoints;
  
  status=poc_get_var_length(var,&npoints);
  
  if(npoints!=1){
    double *coords=new double[npoints];
    
    status=poc_get_var(filename,var,coords);
    
    range_t<double> coordR(coords,npoints);
    
    delete[]coords;
    
    if(coordR.min!=coordR.max)
      TRAP_ERR_EXIT(ENOEXEC,"not coded yet for buggy files : %g<="+var.name+"<=%g...\n",coordR.min,coordR.max);
    
    *coord=coordR.min;
    }
  else{
    status=poc_get_var(filename,var,coord);
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_loadHAWAII_NC(const char *filename, tseries_t *serie, mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,m,status;
  
  int verbose=1;
  
  poc_data_t<double> elevation;
  poc_global_t glob;
  
  status=poc_inq(filename,&glob,verbose);
  const poc_var_t
    *lonVar=glob.variables.findP("longitude","LONGITUDE"),
    *latVar=glob.variables.findP("latitude","LATITUDE"),
    *sshVar=glob.variables.findP(
      "sea_surface_height_above_reference_level",
      "SLEV");
  
  elevation.info=*sshVar;
  
  size_t npoints=-1;
  int nframes=-1;
  status=poc_get_var_length(elevation.info,&npoints,0,&nframes);
  
/*------------------------------------------------------------------------------
  coordinates */
  if(npoints==1){
    status=poc_get_buggy_mgr_var(filename,*latVar,&mooring->lat);
    status=poc_get_buggy_mgr_var(filename,*lonVar,&mooring->lon);
    }
  else{
    mooring->lat=NAN;
    mooring->lon=NAN;
    }
  
  if(serie==0)
    return status;
  
  serie->n=nframes;
  serie->nparam=npoints;
  serie->lat=mooring->lat;
  serie->lon=mooring->lon;
  
/*------------------------------------------------------------------------------
  time */
  
  status=poc_gettime(filename,&serie->origin,&serie->t,&nframes);
  
  /* convert to days */
  for(j=0;j<nframes;j++){
    serie->t[j]/=86400.;
    }
  
/*------------------------------------------------------------------------------
  elevation */
  
  elevation.init("");
  
  /* read */
  status=poc_get_var(filename,elevation.info,elevation.data,verbose);
  
  /* scale */
  poc_att_t *att;
  att=elevation.info.attributes.findP("units");
  if(att->as_string()=="millimeters"){
    elevation.scale*=1e-3;
    }
  poc_scale_data(elevation.data,elevation.length,elevation.scale,elevation.offset,&elevation.mask,verbose);
  
  /* dispatch */
  serie->x=new double*[npoints];
  for(i=0;i<npoints;i++)
    serie->x[i]=new double[nframes];
  
  for(m=0,j=0;j<nframes;j++){
    for(i=0;i<npoints;i++,m++)
      serie->x[i][j]=elevation.data[m];
    }
  
  /* clean-up */
  if(npoints==1){
    i=0;
    double *seriexi=serie->x[i];
    
    for(m=0,j=0;j<nframes;j++){
      
      double seriexij=seriexi[j];
      
      if(seriexij==elevation.mask)
        continue;
      
      seriexi[m]=seriexij;
      serie->t[m]=serie->t[j];
      
      m++;
      }
    
    serie->n=m;
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int hawaii_getlength(const char *filename,date_t *start,date_t *finish,mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *-----------------------------------------------------------------------
  Estimate the number of data to be read
-----------------------------------------------------------------------*/
{
  FILE *in=NULL;
  char line[1024],*s=NULL;
  int n;
//  mooring_t mooring;
  int fisrt_header=1,first_data=1;

  in=fopen(filename,"r");
  if(in==0) {
      TRAP_ERR_EXIT(errno,"error opening %s (%d %s)\n",filename,errno,strerror(errno));
    }

  n=0;
  while(!feof(in)) {
    fgets(line,1024,in);
    if(feof(in)) break;
    s=strstr(line,"LAT");
    if(s==0) {
/* *----------------------------------------------------------------------------
      data */
      n+=12;
      if(first_data!=0) {
        *start=hawaii_decode_data(in,line);
        first_data=0;
        }
     *finish=hawaii_decode_data(in,line);
     }
    else {
/* *----------------------------------------------------------------------------
      header */
      if(fisrt_header!=0) {
        hawaii_decode_position_p(line,mooring);
        fisrt_header=0;
        }
      }
    }
  fclose(in);

  printf("hawaii_getlength: Number of observations: %d (~ %f years)\n",n,n/24./30.5/12.);

  return(n);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int read_hawaii_positions(const char *input, pressure_station **sample)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//001POHNPEI 1975  LAT=06 59.2N  LONG=158 14.6E  TIMEZONE=GMT
//002BETIO   1992  LAT=01 21.7N  LONG=172 55.8E  TIMEZONE=GMT
  int k, l, n, nst,length;
  int lon, lat, year;
  float lonmin,latmin;
  char line[300], dum[256],north,east,*s=NULL;
  FILE *in=NULL;
  pressure_station *local=NULL;
//  pressure_station sample;

  in=fopen(input,"r");
  if(in==NULL) {
    STDOUT_BASE_LINE("cannot open file %s, abort... \n",input);
    exit(-1);
    }
  fscanf(in, "%d ", &nst);
  local=new pressure_station[nst];
  for(k=0; k<nst;k++) {
//    fscanf(in," %lf %lf %s",&sample[k]->t,&sample[k]->p,sample[k]->name);
    fgets(line,sizeof(line)-1,in);
    sscanf(line,"%3d",&n);
    s=strstr(line,"LAT");
//    sscanf(&(line[3]),"%s",&(local[k].name));
    length=(int) (s-line) -5-3-1;
    strncpy(local[k].name,&(line[3]),length);
    sscanf(&(line[10]),"%d",&year);
    sscanf(&(line[17]),"%4s",dum);
    sscanf(&(line[21]),"%d",&lat);
    sscanf(&(line[24]),"%f",&latmin);
    sscanf(&(line[28]),"%c",&north);
    switch (north) {
      case 'N':
        local[k].p=lat+latmin/60.0;
        break;
      case 'S':
        local[k].p=-(lat+latmin/60.0);
        break;
      }
    sscanf(&(line[31]),"%5s",dum);
    sscanf(&(line[36]),"%d",&lon);
    sscanf(&(line[40]),"%f",&lonmin);
    sscanf(&(line[44]),"%c",&east);
    switch (east) {
      case 'E':
        local[k].t=lon+lonmin/60.0;
        break;
      case 'W':
        local[k].t=-(lon+lonmin/60.0);
        break;
      }
    sscanf(&(line[45]),"%s",&(local[k].time_zone));
    sscanf(&(line[60]),"%s",&(local[k].filename));
    l=strlen(local[k].name);
    while (local[k].name[l-1]==' '){
      local[k].name[l-1]=0;
      l--;
      }
    for(l=0;l<strlen(local[k].name);l++) {
      if(local[k].name[l]==' ') {
        local[k].name[l]='-';
        }
      }
    s=strstr(local[k].filename,".");
    *s=0;
    }
  fclose(in);
  *sample=local;
  return(nst);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int read_hawaii(FILE *hawaii,double *lat,double *lon,double *time,double *sealevel,date_t firstdate,int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  year,i,l,m,nitems,cpt/*,iskip*/;
  int  cnt,numstation;
  int tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12;
//  int nitems;
  char *s=NULL;
  float firsthalf,skip ;

  char stationname[60],strlat[60],strlon[60],timezone[60];
  char strlat_min[60],strlon_min[60],*direction=NULL,min[60];
  char line[300],*header=NULL,strref[20]="TIMEZONE";
  char strmonth[30],strday[10],strhalf;

  exitIfNull(
    header=(char *) malloc(200*sizeof(char))
    );

  m=0;
  i=0;
  cnt=0;
  cpt=0;
  while(!(feof(hawaii))) {
    s=fgets(line,300,hawaii);
    if(s==0) break;
    header=strstr(line,strref);
    
    if(header != NULL) {

/*-----------------------------------------------------------------------------
     header */
      nitems=sscanf(line,"%3d%8c%4d LAT=%2c %5c LONG=%3c %5c TIMEZONE=%s",&numstation,stationname,&year,strlat,strlat_min,strlon,strlon_min,timezone);

      if (nitems==8) {
        cnt++;
        cpt++;

        *lat=atof(strlat);
        strncpy(min,strlat_min,strlen(strlat_min)-1);
        *lat+=atof(min)/60.0;
        direction=strchr(strlat_min,'S');
        if (direction != NULL ) *lat=0.0-(*lat);

        *lon=atof(strlon);
        strncpy(min,strlon_min,strlen(strlon_min)-1);
        *lon+=atof(min)/60.0;
        direction=strchr(strlon_min,'W');
        if (direction != NULL ) *lon=360.0-(*lon);
        }
      else TRAP_ERR_EXIT(-1,"bad header reading");
      }
    else {

/*-----------------------------------------------------------------------------
   data */
      while(cnt==i) {  /*prevent from not full year (last year!)*/
        nitems=sscanf(line,"%3d%8c%4d%2c%2c%1c",&numstation,stationname,&year,strmonth,strday,&strhalf);
        if(nitems!=6) TRAP_ERR_EXIT(-1,"bad data reading");
        if(cnt==1)
          if(strhalf=='2') firsthalf=0.5;/*prevent from series starting after noon*/
          else firsthalf=0.0;

        nitems=sscanf(line+20,"%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d",&tmp1,&tmp2,&tmp3,&tmp4,&tmp5,&tmp6,&tmp7,&tmp8,&tmp9,&tmp10,&tmp11,&tmp12);

        if(nitems==12) {
          cnt++; /* line counting */
          l=(i-cpt)*12+0; sealevel[l]=(tmp1);
          l=(i-cpt)*12+1; sealevel[l]=(tmp2);
          l=(i-cpt)*12+2; sealevel[l]=(tmp3);
          l=(i-cpt)*12+3; sealevel[l]=(tmp4);
          l=(i-cpt)*12+4; sealevel[l]=(tmp5);
          l=(i-cpt)*12+5; sealevel[l]=(tmp6);
          l=(i-cpt)*12+6; sealevel[l]=(tmp7);
          l=(i-cpt)*12+7; sealevel[l]=(tmp8);
          l=(i-cpt)*12+8; sealevel[l]=(tmp9);
          l=(i-cpt)*12+9; sealevel[l]=(tmp10);
          l=(i-cpt)*12+10;sealevel[l]=(tmp11);
          l=(i-cpt)*12+11;sealevel[l]=(tmp12);
      /* shift in level array: i=12 -> cnt+12 */
          }
        else TRAP_ERR_EXIT(-1,"bad data reading");
        }
      }
    i++;
    if(feof(hawaii)) break;
    if(s==0) break;
    }; /* end of reading file*/

  cnt--;      /* last incrementation before exiting "while" */
  cnt-=cpt;   /* discount cpt headers */
  cnt*=12;    /* cnt -> nmes */

/* *----------------------------------------------------------------------------
  convert time in CNES days */
  //useless lines commented out
  //skip=julian_day(12,31,1949)-julian_day(1,1,1950);
  //iskip=julian_day(firstdate)-julian_day(1,1,1950);
  skip=julian_day(firstdate)+firsthalf-julian_day(1,1,1950);
  //date_t dum;
  //dum=poctime_getdatecnes(0.0,'d'); /* 1st date in timezone = GMT!*/
  //dum=poctime_getdatecnes(skip,'d'); /* 1st date in timezone = GMT!*/

/* *----------------------------------------------------------------------------
  convert time in CNES hours after having corrected for time zone shift */
  if(strcmp(timezone,"GMT+4")==0) for(i=0;i<cnt;i++) time[i]=skip*24.0+i-4;
  if(strcmp(timezone,"GMT+3")==0) for(i=0;i<cnt;i++) time[i]=skip*24.0+i-3;
  if(strcmp(timezone,"GMT+2")==0) for(i=0;i<cnt;i++) time[i]=skip*24.0+i-2;
  if(strcmp(timezone,"GMT+1")==0) for(i=0;i<cnt;i++) time[i]=skip*24.0+i-1;

  if(strcmp(timezone,"GMT")==0)   for(i=0;i<cnt;i++) time[i]=skip*24.0+i;

  if(strcmp(timezone,"GMT-1")==0) for(i=0;i<cnt;i++) time[i]=skip*24.0+i+1;
  if(strcmp(timezone,"GMT-2")==0) for(i=0;i<cnt;i++) time[i]=skip*24.0+i+2;
  if(strcmp(timezone,"GMT-3")==0) for(i=0;i<cnt;i++) time[i]=skip*24.0+i+3;
  if(strcmp(timezone,"GMT-4")==0) for(i=0;i<cnt;i++) time[i]=skip*24.0+i+4;

  if (header!=NULL) free(header);

  return(cnt);

}/*end read_hawaii*/


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void hawaii_write_header(const char *strcode, int year,double lat,double lon, FILE *hawaii)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double dpart,ipart;
  char strlon[20],strlat[20];

  dpart=modf(lat,&ipart);
  
  if(lat>=0)
    sprintf(strlat,"%02.0lf %04.1fN",fabs(ipart),fabs(dpart*60));
  else
    sprintf(strlat,"%02.0lf %04.1fS",fabs(ipart),fabs(dpart*60));

  if(lon<0.0) lon+=360.0;

  if(lon<=180) {
    dpart=modf(lon,&ipart);
    sprintf(strlon,"%03.0lf %04.1fE",fabs(ipart),fabs(dpart*60));
    }
  else {
    lon=360-lon;
    dpart=modf(lon,&ipart);
    sprintf(strlon,"%03.0lf %04.1fW",fabs(ipart),fabs(dpart*60));
    }

  fprintf(hawaii,"%8s%4d  LAT=%s  LONG=%s  TIMEZONE=GMT\n",strcode,year,strlat,strlon);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int hawaii_write_data(const char *strcode,double time, double *data,double mask , FILE *hawaii,double lat,double lon,int refyear)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int half;
  date_t linedate;

  linedate=poctime_getdatecnes(time,'h');

  if(linedate.second==0)
    half=1;
  else
    half=2;

  if(linedate.year==refyear)
    fprintf(hawaii,"%8s%4d%2.2d%2.2d%1d%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f\n",
      strcode,linedate.year,linedate.month,linedate.day,half,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9],data[10],data[11]);
  else {
    hawaii_write_header(strcode,linedate.year,lat,lon,hawaii);
    fprintf(hawaii,"%8s%4d%2.2d%2.2d%1d%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f\n",
      strcode,linedate.year,linedate.month,linedate.day,half,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9],data[10],data[11]);
    }

  return(linedate.year);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void mgr_savehawaii(FILE *hawaii,double lat,double lon,double *residual,double *time,int numstation,const char *stationname,date_t firstdate,int nmes,double t1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int ndays,*gap=NULL,i,j,k,n,delay,year,firsthour;
  double *res=NULL,mask,*t=NULL,t0,data[12],skip;
  char strcode[20],*newstationname=NULL;
  date_t dum;
  double epsilon=1.e-08,chk,delta;

  exitIfNull(
    gap=(int *)malloc(nmes*sizeof(int))
    );
  exitIfNull(
    newstationname=(char *)malloc(20*sizeof(char))
    );

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
/* *----------------------------------------------------------------------------
  mon patch */
  firstdate=poctime_getdatecnes(time[0], 'h'); /* 1st date in timezone = GMT!*/
  firstdate.month=1;
  firstdate.day=1;
  firstdate.second=0;
/* *--------------------------------------------------------------------------*/

  if (firstdate.year>dum.year) firsthour=0;
  else firsthour=floor(dum.second/3600.0+0.5);

/* *----------------------------------------------------------------------------
  detect holes in time serie*/
  n=nmes;
  for (i=1;i<nmes;i++) {
    delta=time[i]-time[i-1];
    chk=(delta-1.-epsilon)*(delta-1.+epsilon);
    if(chk>0) gap[i]=0;
//    gap[i]=1;
    gap[i]=(int)  floor(time[i]-time[i-1]+0.5);
    if (gap[i]!=1) {
      n+=gap[i]-1;
      }
    }

  if(n>nmes) n=n+(12-n%12)+12+24;
  else n+=24;

  exitIfNull(
    res=(double *)malloc(n*sizeof(double))
    );
  exitIfNull(
    t=(double *)malloc(n*sizeof(double))
    );

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

  skip =julian_day(firstdate)-julian_day(1,1,1950);/* 1st date in timezone = GMT+...*/
  ndays=julian_day(firstdate)-julian_day(1,1,firstdate.year);/*in days since 1/1/year */

  if (firstdate.year==dum.year) {
    delay=0;   /* no shift with GMT */
    skip=julian_day(1,1,dum.year)-julian_day(1,1,1950);/* skip between cnes ref and 1/1/year (if serie start after 1/1) */
    ndays=julian_day(dum)-julian_day(1,1,firstdate.year);/*in days since 1/1/year */
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


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_savehawaii(const char *filename,double lat,double lon,double *h,double *t,int code,char *name,date_t firstdate,int nvalues,double t1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  double mask=9.999;
  FILE *hawaii=NULL;

  hawaii=fopen(filename,"w");
  if(hawaii==0) {
    TRAP_ERR_EXIT(errno,"error opening %s (%d %s)\n",filename,errno,strerror(errno));
    }

  poctime_convert(t,nvalues,'d','h');
  for (i=0; i< nvalues; i++) {
    if(h[i]!=mask) {
      h[i]=h[i]*1000.0;
      }
    else h[i]=mask;
    }
  mgr_savehawaii(hawaii,lat,lon,h,t,code,name, firstdate, nvalues, t1*24.);
  fclose(hawaii);
  return(0);
}



