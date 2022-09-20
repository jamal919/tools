/* ----------------------------------------------------------------------
 * Copyright (c) 2001-2009 LEGOS
 * All rights reserved.
 *
 * Redistribution and use  in source  and binary  forms,  with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *   1. Redistributions of  source  code must retain the  above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice,  this list of  conditions and the following disclaimer in
 *      the  documentation  and/or  other   materials provided  with  the
 *      distribution.
 *
 * THIS  SOFTWARE IS PROVIDED BY  THE  COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND  ANY  EXPRESS OR IMPLIED  WARRANTIES,  INCLUDING,  BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES  OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR  PURPOSE ARE DISCLAIMED. IN  NO EVENT SHALL THE COPYRIGHT
 * HOLDERS OR      CONTRIBUTORS  BE LIABLE FOR   ANY    DIRECT, INDIRECT,
 * INCIDENTAL,  SPECIAL,  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN  CONTRACT, STRICT LIABILITY, OR
 * TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE   OF THIS SOFTWARE, EVEN   IF ADVISED OF   THE POSSIBILITY OF SUCH
 * DAMAGE.
 *
 * DATE OF CREATION: 2001
 * AUTHOR: Thierry LETELLIER, Laurent ROBLOU
 * CONTACT: thierry.letellier@free.fr,laurent.roblou@legos.obs-mip.fr
 *
 * MESSAGE FROM THE AUTHORS:
 * In consideration of this free of charge access, users shall communicate
 * the results obtained directly  from the software. They shall inform the
 * authors of such publications.
 * Thanks!
 *
 ---------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "analyse_alti.h"
#include "aktarus.h"

/* Functions list: */
/* void ERROR: generates PT_#.error file containing some informations about point #*/
/* void sortie_MGR: generates a file in mgr format for Xscan */
/* void sortieIOS_std: generates a .ios.arm file containing harmonic constants in the standard ios format*/
/* void sortie_BAR: generates a .bar.arm file containing harmonic constants and the corresponding error bars*/
/* float module: returns the modulus of a complex */
/* float phase: returns the phase of a complex */
/* int diffmodel: makes statistics on harmonic constants/GOT2000 harmonic constants difference */
/* void save_hawaii: saves a 1D, hourly sampled serie*/
/* void  hawaii_write_header: writes the GLOSS/Hawaii header*/
/* void  hawaii_write_data: writes the 12 data of a line */


/*-------------------------------------------------------------------------------*/
 void ERROR(int nopt,int nmes,int nbond,double lat,double lon,FILE *error)
/*-------------------------------------------------------------------------------*/
{

  char xfilenamex[256];
  
   
  fprintf(error,"nombre de mesure=%d inferieur à 2*nombre d'ondes à analyser=%d",nmes,nbond);
  fprintf(error,"situation lat %lf  lon %lf",lat,lon);
  fclose(error);
}


/*-------------------------------------------------------------------------------*/
void sortie_MGR(double lat, double lon, wavelist work_s, float *amplitude, float *retard,FILE *mgrfile, char form, double Trepet, int numstation, char *stationname,int startyear,int endyear,double trecord,char *version)
/*-------------------------------------------------------------------------------*/
{
  int i;

 for(i=0;i<80;i++) fprintf(mgrfile,"_");
 fprintf(mgrfile,"\n \n");

/*ecriture selon format d'entree et type de donnees*/
/*###################################################*/
 if ((form=='l')&&(Trepet==35))
   fprintf(mgrfile," STATION No%5.4d  : %s_%05d\n ORIGINE          : ERS (%s)\n",numstation,stationname,numstation,version);
 else
   if ((form=='l')&&(Trepet==9.9156))
     fprintf(mgrfile," STATION No%5.4d  : %s_%05d\n ORIGINE          : T/P (%s)\n",numstation,stationname,numstation,version);
 else
   if ((form=='l')&&(Trepet==0.041666))
     fprintf(mgrfile," STATION No%5.4d  : %s_%05d\n ORIGINE          : TG (%s)\n",numstation,stationname,numstation,version);
 else
   if ((form=='l'))
     fprintf(mgrfile," STATION No%5.4d  : %s_%05d\n ORIGINE          : unknown (%s)\n",numstation,stationname,numstation,version);

 else
   if (form=='h')
     fprintf(mgrfile," STATION No%5.4d  : %8s\n ORIGINE          : hawai (%s)\n",numstation,stationname,version);
/*###################################################*/


 fprintf(mgrfile," LOCALISATION     :%10.3fN %8.3fE    0.00m Triangle :     0\n",lat,lon);

 fprintf(mgrfile," ENREGISTREMENT   : Debut: %d        Fin: %d         Duree: %5.0f     \n",startyear,endyear,trecord);
 fprintf(mgrfile," VALIDATION       : no   \n");
 if (form=='h') for (i=0;i<work_s.n;i++)  fprintf(mgrfile,"%3d  %-6s        %9.6f   %10.6f\n",work_s.wave[i].code,work_s.wave[i].name,amplitude[i]/1000.0,retard[i]);
 else  for (i=0;i<work_s.n;i++)  fprintf(mgrfile,"%3d  %-6s        %9.6f   %10.6f\n",work_s.wave[i].code,work_s.wave[i].name,amplitude[i],retard[i]);
}



/*-------------------------------------------------------------------------------*/
void sortieIOS_std(char form, int nopt,int nmes,double lat,double lon,wavelist work_s,float *amplitude,float *retard,FILE *IOS,char *stationname,int numstation,double t0,double t1,double t2,double rms1,double rms2,double mean)
/*-------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------*/
/* this subroutine tries to respect the standard ios.arm format                   */
/* ios.arm generated readable by xscan anyway                                    */
/*-------------------------------------------------------------------------------*/
{
  int i;
  double dpart,ipart;
  char strlon[20],strlat[20];
  date_t date0,date1,date2;

  fprintf(IOS,"%s\n",stationname);

  dpart=modf(lat,&ipart);
  if(dpart>=0) sprintf(strlat,"%2.0lf %4.1fN",fabs(ipart),fabs(dpart*60));
  else sprintf(strlat,"%2.0lf %4.1fS",fabs(ipart),fabs(dpart*60));
  fprintf(IOS,"%s\n",strlat);

  if(lon<=180) {
      dpart=modf(lon,&ipart);
      sprintf(strlon,"%3.0lf %4.1fE",fabs(ipart),fabs(dpart*60));
    }
  else
    {
      lon=360-lon;
      dpart=modf(lon,&ipart);
      sprintf(strlon,"%3.0lf %4.1fW",fabs(ipart),fabs(dpart*60));
    }
  fprintf(IOS,"%s\n",strlon);


  getcnesdate(t0,&date0);/*time of the first record*/
  date0.year-=(date0.year/100)*100;
  getcnesdate(t1,&date1);/*time of the last record*/
  date1.year-=(date1.year/100)*100;
  getcnesdate(t2,&date2);/*time of the mid record (approximatively)*/
  date2.year-=(date2.year/100)*100;

  fprintf(IOS,"1analysis of hourly tidal heights  stn  %3d    %2.0fh %02d/%02d/%02d to  %2.0fh %02d/%02d/%02d\n",numstation,floor(date0.second/3600),date0.day,date0.month,date0.year,floor(date1.second/3600),date1.day,date1.month,date1.year);

  fprintf(IOS,"no.obs.=%6d  no.pts.anal.=%6d midpt=%1.0fh %02d/%02d/%02d separation=0.00\n",nmes,nmes,floor(date2.second/3600),date2.day,date2.month,date2.year);

  if(form =='h') fprintf(IOS,"RMS before Analysis : %-9.4g cm RMS after Analysis : % -9.4g cm\n",rms1/10.0,rms2/10.0);
  else fprintf(IOS,"RMS before Analysis : %-9.4g U. RMS after Analysis : % -9.4g U.\n",rms1,rms2);

  fprintf(IOS,"  no name    frequency   stn  m-y/ m-y  a  g\n");
  if(form =='h') fprintf(IOS,"  0 S0     .0000000000    %3d        %12.4f\n",numstation,mean/10);
  else fprintf(IOS,"  0 S0     .0000000000    %3d        %12.4f\n",numstation,mean*100);/* .list in m!*/
  
 for (i=0;i<work_s.n;i++) {
    if(form == 'h')  fprintf(IOS,"%3d %-6s%11.10f    %3d            %8.4f %6.2f\n",work_s.wave[i].code,work_s.wave[i].name,work_s.wave[i].omega/360.0,numstation,amplitude[i]/10,retard[i]);
      else fprintf(IOS,"%3d %-6s%11.10f    %3d            %8.4f %6.2f\n",work_s.wave[i].code,work_s.wave[i].name,work_s.wave[i].omega/360.0,numstation,amplitude[i]*100,retard[i]);
    }
  fclose(IOS);
}

/*-------------------------------------------------------------------------------*/
void sortie_BAR(char form, int nopt,int nmes,double lat,double lon,wavelist work_s,float *amplitude,float *retard,FILE *BAR,char *stationname,int numstation,double t0,double t1,double t2,double rms1,double rms2,double mean,double *noise)
/*-------------------------------------------------------------------------------*/
{
  int i;
  double dpart,ipart;
  char strlon[20],strlat[20];
  date_t date0,date1,date2;

  fprintf(BAR,"%s\n",stationname); /* line 1*/

  dpart=modf(lat,&ipart);
  if(dpart>=0) sprintf(strlat,"%2.0lf %4.1fN",fabs(ipart),fabs(dpart*60));
  else sprintf(strlat,"%2.0lf %4.1fS",fabs(ipart),fabs(dpart*60));
  fprintf(BAR,"%s\n",strlat);  /* line 2 */

  if(lon<=180) {
      dpart=modf(lon,&ipart);
      sprintf(strlon,"%3.0lf %4.1fE",fabs(ipart),fabs(dpart*60));
    }
  else
    {
      lon=360-lon;
      dpart=modf(lon,&ipart);
      sprintf(strlon,"%3.0lf %4.1fW",fabs(ipart),fabs(dpart*60));
    }
  fprintf(BAR,"%s\n",strlon); /* line 3 */


  getcnesdate(t0,&date0);/*time of the first record*/
  date0.year-=(date0.year/100)*100;
  getcnesdate(t1,&date1);/*time of the last record*/
  date1.year-=(date1.year/100)*100;
  getcnesdate(t2,&date2);/*time of the mid record (approximatively)*/
  date2.year-=(date2.year/100)*100;

  fprintf(BAR,"1analysis of hourly tidal heights  stn  %3d    %2.0fh %02d/%02d/%02d to  %2.0fh %02d/%02d/%02d\n",numstation,floor(date0.second/3600),date0.day,date0.month,date0.year,floor(date1.second/3600),date1.day,date1.month,date1.year);  /* line 4 */

  fprintf(BAR,"no.obs.=%6d  no.pts.anal.=%6d midpt=%1.0fh %02d/%02d/%02d separation=0.00\n",nmes,nmes,floor(date2.second/3600),date2.day,date2.month,date2.year); /* line 5 */

  if(form =='h') fprintf(BAR,"RMS before Analysis : %-9.4g cm RMS after Analysis : % -9.4g cm\n",rms1/10.0,rms2/10.0); /* line 6 */
  else fprintf(BAR,"RMS before Analysis : %-9.4g U. RMS after Analysis : % -9.4g U.\n",rms1,rms2);

  fprintf(BAR,"# Name | freq (1/day)  |  A (cm)  |  G (deg)  | Noise level (cm) | ratio\n");
  fprintf(BAR,"#-----------------------------------------------------------------------\n");
  for(i=0;i<work_s.n;i++)
    if(form == 'h')  fprintf(BAR,"%6s   %11.10f    %9.4f    %6.2f    %9.4f       %9.2f\n",work_s.wave[i].name,work_s.wave[i].omega/360.0, amplitude[i]/10.0,retard[i],noise[i]/10,amplitude[i]/noise[i]);
    else fprintf(BAR,"%6s   %11.10f    %9.4f    %6.2f    %9.4f       %9.2f\n",work_s.wave[i].name,work_s.wave[i].omega/360.0, amplitude[i]*100,retard[i],noise[i]*100,amplitude[i]/noise[i]);
 
  fclose(BAR);
}

/*-------------------------------------------------------------------------------*/
int diffmodel(double lat, double lon, fcomplex *tabZ_deti, FILE **ecart,int nopt, int nmes, double *Delta,int i ,wavelist work_s,int *ultro,char *gpath)
/*-------------------------------------------------------------------------------*/
{

char input[256],zz[56],test[256];
int k,verbose,found,ultra;
fcomplex Z,tmp,mask;

 double temp;

 found=0;
 ultra=999;

 for(k=0;k<256;k++) test[k]=gpath[k];
if(work_s.wave[i].code==9) {sprintf(zz,"M2");found=1;ultra=0;}
else if(work_s.wave[i].code==13){sprintf(zz,"S2");found=1;ultra=1;}
else if(work_s.wave[i].code==7) {sprintf(zz,"N2");found=1;ultra=2;}
else if(work_s.wave[i].code==3) {sprintf(zz,"K1");found=1;ultra=3;}
else if(work_s.wave[i].code==1) {sprintf(zz,"O1");found=1;ultra=4;}
else if(work_s.wave[i].code==2) {sprintf(zz,"P1");found=1;ultra=5;}
else if(work_s.wave[i].code==14){sprintf(zz,"K2");found=1;ultra=6;}
else if(work_s.wave[i].code==27){sprintf(zz,"Q1");found=1;ultra=7;}
  
 if (found) {
     sprintf(input,"%s/%s.den",test,zz);
     loadingEC(input,lat,lon,&Z,&mask);
     
     if((Z.reel!=mask.reel)&&(Z.imag!=mask.imag)) {
         
         
         Delta[0]=module(Z)-module(tabZ_deti[i]);
         
         
         Delta[1]=-phase(Z)-phase(tabZ_deti[i]);
         while(Delta[1]<-180)Delta[1] +=360.;
         while(Delta[1]>180)Delta[1] -=360.;
         
     
      
         tmp.reel=Z.reel-tabZ_deti[i].reel;
         tmp.imag=-Z.imag-tabZ_deti[i].imag;
         Delta[2]=module(tmp);
     
         fprintf(ecart[ultra]," %8.3f %8.3f %9.4f %9.4f %9.4f  %d  %d\n",lat,lon,Delta[0] , Delta[1],Delta[2] ,nopt,nmes);
         verbose=0;
       }
     else
       {
         for(k=0;k<8;k++) Delta[k]=0;
         verbose=1;
       }
     
   }/*if found*/
 *ultro=ultra;
     return(verbose);
}/*end.*/


/*-------------------------------------------------------------------------------*/

void save_hawaii(FILE *hawaii,double lat,double lon,double *residual,double *time,int numstation,char *stationname,date_t firstdate,int nmes,double t1)

/*-------------------------------------------------------------------------------*/
{
  int ndays,*gap,i,j,k,n,delay,year,firsthour;
  double *res,mask,*t,t0,data[12],skip;
  char filename[30],strcode[20],*newstationname;
  date_t dum;

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

  getcnesdate(time[0],&dum); /* 1st date in timezone = GMT!*/
  if (firstdate.year>dum.year) firsthour=0;
  else firsthour=floor(dum.second/3600.0+0.5);

  n=nmes;
  for (i=1;i<nmes;i++) {
      gap[i]=(int)(time[i]-time[i-1]);
      if (gap[i]!=1) n+=gap[i];
    }

  if(n>nmes) n=n+(12-n%12)+12+24;
  else n+=24;

  res=(double *)malloc(n*sizeof(double));
  t=(double *)malloc(n*sizeof(double));

  for(i=0;i<n;i++) res[i]=mask;/*!!!!*/
  for (i=0;i<nmes;i++) if (residual[i]==99999) residual[i]=mask; /* for the conversion sonel -> hawaii */

  delay=0;
  i=1;
/*   res[0]=residual[0]; */

  i=1;
  for(j=1;j<nmes;j++) {
    if(gap[i]==1) res[j+delay+firsthour]=residual[i];
    else
      {
      for(k=j;k<(j+gap[i]-1);k++) res[k+delay+firsthour]=mask;
      res[k+delay+firsthour]=residual[i];
      delay+=gap[i]-1;
      }
    i++;
    }
 
  skip=julian_day(firstdate.month,firstdate.day,firstdate.year)-julian_day(1,1,1950);/* 1st date in timezone = GMT+...*/
  ndays=julian_day(firstdate.month,firstdate.day,firstdate.year)-julian_day(1,1,firstdate.year);/*in days since 1/1/year */

  if (firstdate.year==dum.year) {
      delay=0;   /* no shift with GMT */
      skip=julian_day(1,1,dum.year)-julian_day(1,1,1950);/* skip between cnes ref and 1/1/year (if serie start after 1/1) */
      ndays=julian_day(dum.month,dum.day,dum.year)-julian_day(1,1,firstdate.year);/*in days since 1/1/year */
      res[firsthour]=residual[0];
    }
  else
    if (firstdate.year>dum.year) delay=(int)(dum.second-firstdate.second)/3600.0-24;  /*GMT+delay*/
    else
      {
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
  while(t0<(t1-12)&&(i+12<=n)) /* modif:'i+12<n' before, LR, 30/04/04 */
    {
      t0=t[i];
      for(j=0,k=0;j<12;j++,k++) data[k]=res[i-delay+j]; /* +delay*/
      year=hawaii_write_data(strcode,t0,data,mask,hawaii,lat,lon,year);
      i+=12;
    }

  free(res);
  free(gap);
  free(t);

}/*end save_hawaii*/


/*-------------------------------------------------------------------------------*/

  void hawaii_write_header_obsolete(char *strcode, int year,double lat,double lon, FILE *hawaii)

/*-------------------------------------------------------------------------------*/
{
#warning buggy duplicate
  double dpart,ipart;
  char strlon[20],strlat[20];

  dpart=modf(lat,&ipart);
  if(dpart>=0) sprintf(strlat,"%02.0lf %04.1fN",fabs(ipart),fabs(dpart*60));
  else sprintf(strlat,"%02.0lf %04.1fS",fabs(ipart),fabs(dpart*60));

  if(lon<=180) {
    dpart=modf(lon,&ipart);
    sprintf(strlon,"%03.0lf %04.1fE",fabs(ipart),fabs(dpart*60));
    }
  else  {
    lon=360-lon;
    dpart=modf(lon,&ipart);
    sprintf(strlon,"%03.0lf %04.1fW",fabs(ipart),fabs(dpart*60));
    }

  fprintf(hawaii,"%8s%4d  LAT=%s  LONG=%s  TIMEZONE=GMT\n",strcode,year,strlat,strlon);
}


/*-------------------------------------------------------------------------------*/

int hawaii_write_data(char *strcode,double time, double *data,double mask , FILE *hawaii,double lat,double lon,int refyear)

/*-------------------------------------------------------------------------------*/
{
  int half;
  date_t linedate;

  getcnesdate(time,&linedate);
  if(linedate.second==0) half=1;
  else half=2;

  if(linedate.year==refyear)
  fprintf(hawaii,"%8s%4d%2d%2d%1d%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f\n",strcode,linedate.year,linedate.month,linedate.day,half,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9],data[10],data[11]);
else
  {
  hawaii_write_header(strcode,linedate.year,lat,lon,hawaii);
  fprintf(hawaii,"%8s%4d%2d%2d%1d%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f\n",strcode,linedate.year,linedate.month,linedate.day,half,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9],data[10],data[11]);
  }

  return(linedate.year);

}


/*-----------------------------------------------------------------------------------------------*/

void list_save_serie(int pos, double lon, double lat, int nopt, FILE *out,int n,double *time, double *sealevel)

/*-----------------------------------------------------------------------------------------------*/
{
  int i;
  date_t startdate,enddate;
  
  
  fprintf(out,"#-- HEADER -------------------------------------------\n");
  fprintf(out,"# Column 1 : date in days referred to\n");
  fprintf(out,"#            CNES date (01-JAN-1950 00:00:00.0)\n");
  fprintf(out,"# Column 2 : sea surface height (in meters)\n");
  fprintf(out,"#-- HEADER END ---------------------------------------\n");
  fprintf(out,"Number of crossover points :\n");
  fprintf(out,"%d\n",1);
  fprintf(out,"#-----------------------------------------------------\n");
  fprintf(out,"Pt  : %d\n",nopt);
  fprintf(out,"lon : %lf\n",lon);
  fprintf(out,"lat : %lf\n",lat);
  fprintf(out,"Mes : %d\n",n);

  for(i=pos;i<pos+n;i++)  fprintf(out,"%13.6f  %9.4f\n",time[i]/24,sealevel[i]);

}


