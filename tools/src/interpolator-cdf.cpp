

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

#include "tools-structures.h"

#include "poc-netcdf.def"

#include "fe.h"
#include "map.h"
#include "geo.h"
#include "sym-io.h"
#include "polygones.h"
#include "grd.h"
#include "netcdf-proto.h"
#include "poc-time.h"

typedef struct
  {
  char    name[60];
  double  t,p;
  double  time;
  date_t  date;
  float   h,depth,quality;
  } sample_t;
// 2007 298.628657 5.516000 -51.094100 108.939217 0.298290

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

static double dday_in_year(date_t date)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *-----------------------------------------------------------------------

  compute time elapsed in year in decimal days

!-----------------------------------------------------------------------*/

{
  int jd;
  double days;
  jd=julian_day(date.month,date.day,date.year)- julian_day(1, 1,date.year);

  days=(double) jd + date.second/3600./24.;
  return(days);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  date_t getdate_fromstring(char *text)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int option, status, ncid, k;
  size_t lengthhp, index[1], *dimlgth;
  int vtime, d;
  int day, month, year, jday;
  int hour, minute, seconds;
  nc_type xtypep;
  char name[1024];
  char smonth[1024];
  char *pointer, *tmp;
  double factor;
  date_t date;

  tmp=strndup(text,256);

  sscanf(strtok(tmp, "-"), "%d", &year);
  if((pointer = strtok(NULL, "-")) == NULL)
    goto error;
  sscanf(pointer, "%3s", &smonth);
  if((pointer = strtok(NULL, " ")) == NULL)
    goto error;
  sscanf(pointer, "%d", &day);

  hour    = 0;
  minute  = 0;
  seconds = 0;

//  tmp=strndup(text,256);

  pointer = strtok(NULL, ":");
  if(pointer != NULL)
    sscanf(pointer, "%d", &hour);
  pointer = strtok(NULL, ":");
  if(pointer != NULL)
    sscanf(pointer, "%d", &minute);
  pointer = strtok(NULL, "\0");
  if(pointer != NULL)
    sscanf(pointer, "%d", &seconds);

  for(k = 0; k < 12; k++)
    if(strcmp(smonth, uc_names[k]) == 0)
      break;
  month = k + 1;
/*------------------------------------------------------------------------
  patch: month given as number and not text*/
  if(k == 12)
    sscanf(smonth, "%d", &month);

  date.year=year;
  date.month=month;
  date.day=day;
  date.second=hour*3600+minute*60+seconds;

  return (date);

error:
  return (date);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float z;
  double x,y;
  float mask=0;

  float  *buffer[2], **serie;
  double sampling,*time,days;
  int fmt,i,j,k,l,n;
  int year,month,frame,status,nst;

  FILE *in,*out;
  char text[256];
  char outname[256]="\0",*rootname=NULL;
  sample_t *sample;
  char *option,*keyword,*s,*varname=NULL,*input,*datafile;
  cdfvar_t info;
  cdfgbl_t data_info,grid_info;
  date_t origine,reference;
  grid_t grid;

  fct_echo(argc,argv);
 
  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'i' :
          input= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'r' :
          rootname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          sscanf(argv[n+1],"%lf",&sampling);
          n++;
          n++;
          break;

        case 'v' :
          varname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        datafile= strdup(argv[n]);
        n++;
        break;
      }
      free(keyword);
    }

  l=strlen(datafile);
  if(strrchr(datafile,(int) '/') != NULL) {
    l=strlen(strrchr(datafile,(int) '/'));
    rootname=new char[1024];
    strncpy(rootname,datafile,l-l+1);
    }
  else  rootname= strdup("./");

 // printf("rootname: %s \n",rootname);

//  status= poc_gettime   (datafile,  info,  grid_info, &origine, &time, &nframes);

  in=fopen(input,"r");
  if(in==NULL) {__ERR_BASE_LINE__("exiting\n");exit(-1);}
  fscanf(in, "%d ", &nst);
  sample=new sample_t[nst];

//   for(k=0; k<nst;k++) {
//     fscanf(in," %lf %lf",&sample[k].t,&sample[k].p);
//     fgets(text,256,in);
//     sample[k].date=getdate_fromstring(text);
//     sample[k].time=cnes_time(sample[k].date,'s');
//     printf(" %f %f %s \n",sample[k].t,sample[k].p,sample[k].name);
//     }
  for(k=0; k<nst;k++) {
    fscanf(in,"%d %lf %lf %lf %f %f",&year,&days,&sample[k].p,&sample[k].t,&sample[k].depth,&sample[k].quality);
    fgets(text,256,in);
    reference.year=year;
    reference.month=1;
    reference.day=1;
    reference.second=0.;
    sample[k].date=poctime_getdate(days,reference,'d', &(sample[k].time));
    sample[k].time=cnes_time(sample[k].date,'s');
//    printf(" %f %f %lf \n",sample[k].t,sample[k].p,sample[k].time);
    }
  fclose(in);

  month=sample[0].date.month;

  status= cdf_globalinfo(datafile,&grid_info,0);
  status= cdf_varinfo  ((const char*) datafile,"elevation",&info,0);

  status= poc_getgrid2d (datafile, grid_info, info, &grid);

  buffer[0]=new float[grid.nx*grid.ny];
//  buffer[1]=new float[grid.nx*grid.ny*grid.nz];

  serie=new float*[nst];
  for(k=0; k<nst;k++) {
    serie[k]=new float[grid.nt];
    }

  for(frame=0;frame<grid.nt;frame++) {
    status= poc_getvar2d (datafile, info.id, frame,(float *) buffer[0], &mask ,info);
    for (n=0;n<nst;n++) {
      x=sample[n].t;
      y=sample[n].p;
      status= map_interpolation(grid, buffer[0], mask,x,y,&z);
//      status= map_interpolation(grid, buffer[1], mask,x,y,&z);
//      sample[n].h=z;
      serie[n][frame]=z;
      }
    }

  for (k=0; k<nst; k++) {
    status= map_interpolate1D(serie[k], grid.time, mask, grid.nt, sample[k].time, &z);
//    printf("%lf %lf %lf %f\n",sample[k].time,sample[k].t,sample[k].p,z);
    printf("%d %lf %lf %lf %f %f %f\n",sample[k].date.year,dday_in_year(sample[k].date),sample[k].p,sample[k].t,sample[k].depth,sample[k].quality,z);
    }

  __ERR_BASE_LINE__("%s -computation complete ^^^^^^^^^\n",argv[0]);
  exit(0);
}




