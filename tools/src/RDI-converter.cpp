

/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Cyril Nguen        LA, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/
/* *-----------------------------------------------------------------------

  It reads RDI ADCP ascii files and extract bathymetry

-----------------------------------------------------------------------*/


#include <stdio.h>
#include <string.h>

#include "functions.h"
#include "poc-assertions.h"
#include "poc-time.h"

#include "tools-define.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int RDI_read(const char *input, const char *output)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in,*out;
  double dyear,dmonth,dday,dhour,dminute,dsec,temp;
  double cnestime,lon,lat;
  int year,month,hour,day,minute,second,jday,nlayers;
  float dum,h1,h2,h3,h4,h;
  float depth,velocity,cap,u,v;
  float beamstatus[4];
  int n,nitems;
  char line[128], *s;
  int count,GPS_loss;

  in =fopen(input,"r");
  out=fopen(output,"w");

  printf("converting %s\n",input);

  s=fgets(line,128,in);
  if(s==0) goto finished;
  s=fgets(line,128,in);
  if(s==0) goto finished;
  s=fgets(line,128,in);
  if(s==0) goto finished;

  count=0;
  while (!feof(in)) {
    GPS_loss=0;
/* *----------------------------------------------------------------------------
    Decode header*/
    s=fgets(line,128,in);
    if(s==0) goto finished;
    sscanf(line,"%d %d %d %d %d %d",&year,&month,&day,&hour,&minute,&second);
    year+=2000;
    s=fgets(line,128,in);
    if(s==0) goto finished;
    sscanf(line,"%f %f %f %f %f %f %f %f %f %f %f %f",&dum,&dum,&dum,&dum,&dum,&dum,&dum,&dum,&h1,&h2,&h3,&h4);
    s=fgets(line,128,in);
    if(s==0) goto finished;
    s=fgets(line,128,in);
    if(s==0) goto finished;
    sscanf(line,"%lf %lf",&lat,&lon);
    if(lat==30000.) {
/* *----------------------------------------------------------------------------
      LGPS position unavailble*/
//      printf("GPS loss at %d\n",count);
      GPS_loss=1;
      }
    s=fgets(line,128,in);
    if(s==0) goto finished;
    s=fgets(line,128,in);
    if(s==0) goto finished;
/* *----------------------------------------------------------------------------
    Loop on ADCP vertical cells*/
    sscanf(line,"%d",&nlayers);
    for (n=0;n<nlayers;n++) {
      s=fgets(line,128,in);
      if(s==0) goto finished;
      nitems=sscanf(line,"%f %f %f %f %f %f %f %f %f %f %f %f %f",&depth,&velocity,&cap,&u,&v,&dum,&dum,
                                                                  &beamstatus[0],&beamstatus[1],&beamstatus[2],&beamstatus[3],&dum,&dum);
      }
    h=h1;
    h=MIN(h,h2);
    h=MIN(h,h3);
    h=MIN(h,h4);
    if((h!=0.0) && (GPS_loss==0)) {
      nitems=fprintf(out," %lf %lf %f\n",lon,lat,-h);
      }
    count++;
    }

finished:
  fclose(in);
  fclose(out);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in,*out,*out2;
  int n,status;
  char line[128], *s;

  fct_echo(argc,argv);
 
  status=RDI_read(argv[1],argv[2]);

  __ERR_BASE_LINE__("exiting\n");wexit(0);
}
