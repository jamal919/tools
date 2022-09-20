
/*******************************************************************************

  T-UGO tools, 2006-2016

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Analyses filtered/detided TUGOm sample files.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main <a href=#func-members>function</a> for how this works
and the print_help <a href=#func-members>function</a> for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "poc-time.h"
#include "filter.h"
#include "tides.h"
#include "mgr.h"
#include "mgr-converter.h"
#include "legend.h"
#include "functions.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : as[0]
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    "NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Analyses filtered/detided TUGOm sample files.\n"
    "\n"
    "OPTIONS\n"
    "  -h,--help  Show this help and exit.\n"
    "  -start  followed by start date in mm/yyyy format\n"
    "  -end  followed by end date in mm/yyyy format\n"
    "  -cm  scale to centimeters\n"
    "  -check  produce filtered time series\n"
    "  --no-fft  disable costly production of FFTs\n"
    "  -d  followed by observations directory. Default: ./\n"
    "  -e  followed by observations extension. Default: .res.65.mlf.obs\n"
    "  -m  followed by samples directory. Default: ./\n"
    "  -t  followed by translation file\n"
    "  -i  followed by stations list\n"
    );
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double mask=999.9,mask2;
  double dT[4],origine,t1=0,t2=1.e+10;
  double *h[4],*t[4],*buffer;/* [p,mod,obs,dif] */

  double start,finish,time,epsilon,factor=100.0;

  int i,k,l,m,n,status;
  int nr[4],nst,test=0;
  FILE *in,*out,*statfile;
  char *name,tmp[256];
  char *model=NULL,*data=NULL, *input, *translation=NULL;
  char *extension=NULL;
  bool add_line=false,doFourier=true;
  pressure_station *sample;
  statistic_t stat[3],sdum;
  int start_year=1900,last_year=2100,start_month=1,last_month=12;

  char  location[1000][256],sample_name[1000][256],station[1000][256];
  float georef[1000],ignref[1000];
  int   TGid,nitems,nsample;
  float TGlon[1000],TGlat[1000];
  
 // vector<statistic_t> gstat[3];
  vector<point_t> points;

  fct_echo( argc, argv);

  fprintf(stderr,"%s  -starting computation *********\n",argv[0]);

  n=1;
  while (n < argc) {
    const char *keyword=argv[n];
    
    if(strcmp(keyword,"--help")==0 || strcmp(keyword,"-h")==0) {
      print_help(argv[0]);
      exit(0);
      }
    
    if(strcmp(keyword,"-start") ==0) {
      status=sscanf(argv[n+1],"%d/%d",&start_month,&start_year);
      n++;
      n++;
      continue;
      }
    
    if(strcmp(keyword,"-end") == 0) {
      status=sscanf(argv[n+1],"%d/%d",&last_month,&last_year);
      n++;
      n++;
      continue;
      }
    
    if(strcmp(keyword,"-cm") == 0) {
      factor=1.0;
      n++;
      continue;
      }
    
    if(strcmp(keyword,"-check") == 0) {
      test=1;
      n++;
      continue;
      }
    
    if(strcmp(keyword,"--no-fft") == 0) {
      doFourier=false;
      n++;
      continue;
      }
    
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'd' :
          data= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'e' :
          extension= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          model= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 't' :
          translation= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'i' :
          input= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
          break;
        }
        break;

      default:
        STDOUT_BASE_LINE("unknown option %s\n",keyword);
        exit(-1);
        break;
      }
    
    }

  if(extension==NULL) {
    extension=strdup(".res.65.mlf.obs");
    printf("use <%s> for default observation extension\n",extension);
    }
 
  if(data==NULL) {
    data=strdup(".");
    printf("use <%s/> for default input directory\n",data);
    }
  data=(char *)realloc(data,strlen(data)+2);
  strcat(data,"/");

  if( model != NULL) {
    model=(char *)realloc(model,strlen(model)+2);
    strcat(model,"/");
    }
  else model=strdup("./");

  origine=julian_day(1,1,1950);
  t1=julian_day(start_month,1,start_year)-origine;
  if(last_month==12)
    t2=julian_day(1,1,last_year+1)-origine;
  else
    t2=julian_day(last_month+1,1,last_year)-origine;

/* *-----------------------------------------------------------------------------
  read the input file*/
//   in=fopen(input,"r");
// 
//   fscanf(in, "%d ", &nst);
// 
//   sample=(pressure_station *)malloc(nst*sizeof(pressure_station));
// 
// /* *-----------------------------------------------------------------------------
//   read the station list, and retains stations in latitude limits*/
//   n=0;
//   for(k=0; k<nst;k++) {
//     fscanf(in," %lf %lf %s",&sample[n].t,&sample[n].p,sample[n].name);
//     printf(" %f %f %s \n",sample[n].t,sample[n].p,sample[n].name);
//     if((abs(sample[n].p) < latmin) ) {
//       continue;
//       }
//     else if((abs(sample[n].p) > latmax) ) {
//       continue;
//       }
//     else {
//       n++;
//       }
//     }
//   fclose(in);
// 
//   nst=n;

  nst=mgr_read_stations(input,&sample);

/* *-----------------------------------------------------------------------------
  read the translation table for names equivalence*/
  if(translation != NULL) {
    in=fopen(translation,"r");
    if(in==0)
      TRAP_ERR_EXIT(errno,"fopen(\"%s\",\"r\") error (%d %s)\n",translation,errno,strerror(errno));
    nitems=fscanf(in,"%d",&nsample);
    k=0;
    while (!feof(in)) {
      nitems=fscanf(in,"%s %s %f %f %f %f %s",
                  &location[k],&sample_name[k],&TGlon[k],&TGlat[k],&georef[k],&ignref[k],&station[k]);
//       printf("%s %s %f %f %f %f %s\n",
//                   location[k],sample_name[k],TGlon[k],TGlat[k],georef[k],ignref[k],station[k]);
      if(nitems!=7) break;
      k++;
      }
    fclose(in);
    nsample=k;
    }
  else {
    for(k=0; k<nst;k++) {
      strcpy(location[k],sample[k].name);
      strcpy(station[k],sample[k].name);
      }
    nsample=nst;
    }

  statfile=fopen("statistics","a");
  
  name=new char[1024];

  for (k=0;k<nst;k++) {
    string inputPathsAssertion;
/* *-----------------------------------------------------------------------------
    read the model IB file */
    for (l=0;l<nsample;l++) {
      if(strcmp(sample[k].name,location[l])==0) break;
      }
    if(l==nsample) {
      printf("could not find %s in translation file, skip station\n",sample[k].name);
      continue;
      }
    TGid=l;
    strcpy(name,model);
//    strcpy(name,"sample.");
    strcat(name,sample[k].name);
    strcat(name,".p");
    printf("%d/%d: lon=%g lat=%g\n",k+1,nst,sample[k].t,sample[k].p);
    fflush(stdout);
    printf("sea level at %15s: IBC -- ",sample[k].name);
    asprintf(inputPathsAssertion,"(from %s",name);
    status=senetosa_load(name, t1, t2, &h[0], &t[0], mask, &nr[0],&dT[0]);
    if(status!=0) {
      STDOUT_BASE_LINE("aborted (%s)\n",name);
      continue;
      }
    buffer=(double *)malloc(nr[0]*sizeof(double));
    for (i=0; i< nr[0]; i++) {
/* *----------------------------------------------------------------------------
      convert in centimetres if needed */
      if(h[0][i] != mask) buffer[i]=h[0][i]*factor;
      else buffer[i]=mask;
      }
    if(test) {
      sprintf(tmp,"%s.1",name);
      sdum=filtering(tmp,buffer,t[0],mask,nr[0],dT[0]);
      }
    if(doFourier){
      sprintf(tmp,"%s.fft",name);
      fourier(tmp,buffer,mask,nr[0],dT[0],0);
      }
    free(buffer);

/* *-----------------------------------------------------------------------------
    read the model elevation file if needed */
    sprintf(name,"%s%s.h",model,sample[k].name);
    fflush(stdout);
    printf("sea level at %15s: MOD -- ",sample[k].name);
    asprintf(inputPathsAssertion,", %s",name);
    status=senetosa_load(name, t1, t2, &h[1], &t[1], mask, &nr[1],&dT[1]);
    if(status!=0) {
      STDOUT_BASE_LINE("aborted (%s)\n",name);
      continue;
      }
    buffer=(double *)malloc(nr[1]*sizeof(double));
    for (i=0; i< nr[1]; i++) {
/* *----------------------------------------------------------------------------
      convert in centimetres if needed */
      if(h[1][i] != mask) buffer[i]=h[1][i]*factor;
      else buffer[i]=mask;
      }
    if(test) {
      sprintf(tmp,"%s.1",name);
      sdum=filtering(tmp,buffer,t[1],mask,nr[1],dT[1]);
      }
    if(doFourier){
      sprintf(tmp,"%s.fft",name);
      fourier(tmp,buffer,mask,nr[1],dT[1],0);
      }
    free(buffer);
 
/* *-----------------------------------------------------------------------------
    read the observation file  */
    strcpy(name,data);
/*     strcat(name,sample[k].name); */
/*     strcat(name,".obs"); */
    strcat(name,station[TGid]);
    fflush(stdout);
    printf("sea level at %15s: DAT -- ",sample[k].name);
    asprintf(inputPathsAssertion," and %s)",name);
    if(extension[0]!='\0'){
      strcat(name,extension);
      status=senetosa_load(name, t1, t2, &h[2], &t[2], mask, &nr[2],&dT[2]);
      if(status!=0) {
        STDOUT_BASE_LINE("aborted (%s)\n",name);
        continue;
        }
      mask2=mask;
      }
    else{
      tseries_t *serie,reduced,resampled;
      mooring_t mooring;
      int ncol=2;
      
      /* read */
      mgr_load_timeserie(name, &serie, 'm', NULL, ncol, NULL, &mooring, &status);
      if(status!=0){
        STDOUT_BASE_LINE("mgr_load_timeserie(\"%s\",...) error (%d %s)\n",name,status,strerror(status));
        continue;
        }
      
      /* select column */
      for(i=0;i<ncol;i++){
        if(i==1)
          continue;
        delete[]serie->x[i];
        }
      serie->x[0]=serie->x[1];
      serie->nparam=1;
      
      /* select dates */
      const date_t
        start_date=poctime_getdatecnes(t1,'d'),
        end_date  =poctime_getdatecnes(t2,'d');
      
      reduced=serie->reduce(start_date,end_date);
      
      if(reduced.n>0){
        double sampling=get_timesampling(reduced.t, reduced.mask, reduced.n);
        sampling=round(sampling*86400.)/86400.;
        resampled=reduced.resample(sampling, sampling*2);
        
        mask2=resampled.mask;
        h[2]=resampled.x[0];
        t[2]=resampled.t;
        nr[2]=resampled.n;
        dT[2]=sampling;
        printf("%d obs %gs-sampled in [%s ; %s]\n",resampled.n,dT[2]*86400.,poctime_sdate_cnes(resampled.t[0],'d'),poctime_sdate_cnes(resampled.t[resampled.n-1],'d'));
        }
      else{
        nr[2]=0;
        printf("*** 0 obs ***\n");
        }
      
      /* clean-up */
      serie->destroy();
      delete[]serie;
      reduced.destroy();
      delete[]resampled.x;
      }
    printf("%s\n",inputPathsAssertion.c_str());
    
    if(nr[2]<=0)
      continue;
    
    buffer=(double *)malloc(nr[2]*sizeof(double));
    for (i=0; i< nr[2]; i++) {
      if(h[2][i] != mask2) buffer[i]=h[2][i]*factor;
      else buffer[i]=mask;
      }

    if(test) {
/*
      strcpy(tmp,name);
      strcat(tmp,".0");
      sdum=filtering(tmp,buffer,t[2],mask,nr[2],dT[2]);
*/
      }
/*
    strcpy(tmp,name);
    strcat(tmp,".0.fft");
    if(doFourier)
      fourier(tmp,buffer,t[2],mask,nr[2],dT[2],0);
*/
    free(buffer);

/* *-----------------------------------------------------------------------------
    start analysis  */
    start =-1.e+10;
    finish=+1.e+10;

    for (i=0;i<3;i++) {
      updatemax(&start,t[i][0]);
      updatemin(&finish,t[i][nr[i]-1]);
      }

    if(start > finish){
      printf("%s > %s : no stat\n",sgetcnesdate(start*24.),sgetcnesdate(finish*24.));
      continue;
      }

    dT[3]=dT[0];
    nr[3]=(int)( floor( (finish-start)/dT[3] )+1 );
    h[3]=(double *)malloc(nr[3]*sizeof(double));
    t[3]=(double *)malloc(nr[3]*sizeof(double));

    out=fopen("out","w");
    for (time=start; time <= finish; time+=dT[3]) {
      i=(int)( NINT((time-t[0][0])/dT[0]) );
      m=(int)( NINT((time-t[2][0])/dT[2]) );
      if(0<=m and m<nr[2] and h[2][m] != mask2) {
        fprintf(out,"%12.6f %12.6f %6.3f %6.3f %6.3f %s\n",
//           t[0][i],t[2][m],h[2][m]-h[0][i],h[2][m]-h[1][i],h[1][i]-h[0][i],
          t[0][i],t[2][m],h[0][i],h[1][i],h[2][m],
          poctime_sdate_cnes(t[0][i],'d'));
        add_line=True;
        }
      else {
        if(add_line) fprintf(out,"\n");
        add_line=False;
        }
      }
    fclose(out);

    epsilon=dT[3]/100.;

/* *-----------------------------------------------------------------------------
    build the common time serie: observations */
    n=-1;
    for (time=start; time <= finish+epsilon and n<nr[3]-1; time+=dT[3]) {
      n++;
      i=(int)( NINT((time-t[0][0])/dT[0]) );
      m=(int)( NINT((time-t[2][0])/dT[2]) );
      if(0<=m and m<nr[2] and h[2][m] != mask2) {
/*-----------------------------------PROVISOIRE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
          /* h[2][m]=h[2][m]+h[0][i]; */ /* Capraia ENEA = BPR */
/* *----------------------------------------------------------------------------
        convert in centimetres */
        h[3][n]=h[2][m]*factor;
        }
      else{
        h[3][n]=mask;
        }
      t[3][n]=t[0][i];
      }
    
    printf("sea level observation at %s: no correction -------------\n",sample[k].name);
    fflush(stdout);
    strcpy(name,model);
    strcat(name,sample[k].name);
    strcat(name,".obs.1");
    stat[0]=filtering(name,h[3],t[3],mask,nr[3],dT[3]);
    printf("#topex 1 %10s samples: %6d; mean: %8.4lf; variance: %8.4lf\n", sample[k].name,n,stat[0].mean,stat[0].variance());

    if(doFourier){
      strcpy(name,model);
      strcat(name,sample[k].name);
      strcat(name,".obs");
      sprintf(tmp,"%s.fft",name);
      fourier(tmp,h[3],mask,nr[3],dT[3],0);
      }
/*     season("",h[3],t[3],mask,nr[3],dT[3]); */

/* *-----------------------------------------------------------------------------
    build the common time serie: observations minus IB */
    n=-1;
    for (time=start; time <= finish+epsilon and n<nr[3]-1; time+=dT[3]) {
      n++;
      i=(int)( NINT((time-t[0][0])/dT[0]) );
      m=(int)( NINT((time-t[2][0])/dT[2]) );
      if(0<=m and m<nr[2] and h[2][m] != mask2 and h[0][i] != mask)
        h[3][n]=(h[2][m]-h[0][i])*factor;
      else
        h[3][n]=mask;
      }
    printf("sea level observation at %s: inverse barometer corrected--\n",sample[k].name);
    fflush(stdout);
    strcpy(name,model);
    strcat(name,sample[k].name);
    strcat(name,".p.dif");
    stat[1]=filtering(name,h[3],t[3],mask,nr[3],dT[3]);
    printf("#topex 2 %10s samples: %6d; mean: %8.4lf; variance: %8.4lf\n", sample[k].name,n,stat[1].mean,stat[1].variance());
    if(doFourier){
      sprintf(tmp,"%s.fft",name);
      fourier(tmp,h[3],mask,nr[3],dT[3],0);
      }
/*     season("senetosa.dif",h[3],t[3],mask,nr[3],dT[3]); */

/* *-----------------------------------------------------------------------------
    build the common time serie: observations minus model */
    n=-1;
    for (time=start; time <= finish+epsilon and n<nr[3]-1; time+=dT[3]) {
      n++;
      i=(int)( NINT((time-t[0][0])/dT[0]) );
      m=(int)( NINT((time-t[2][0])/dT[2]) );
      if(0<=m and m<nr[2] and h[2][m] != mask2 and h[1][i] != mask) {
        h[3][n]=(h[2][m]-h[1][i])*factor;
        }
      else {
//         if(h[1][i]==mask)
//           STDERR_BASE_LINE("t[<dif>][%d]=%.2f;h[<mod>][%d]=%g\n",i,t[3][i],i,h[1][i]);
        h[3][n]=mask;
        }
      }
    printf("sea level observation at %s: model output corrected-------\n",sample[k].name);
    fflush(stdout);
    strcpy(name,model);
    strcat(name,sample[k].name);
    strcat(name,".h.dif");
    stat[2]=filtering(name,h[3],t[3],mask,nr[3],dT[3]);
    printf("#topex 3 %10s samples: %6d; mean: %8.4lf; variance: %8.4lf\n", sample[k].name, n, stat[2].mean, stat[2].variance());
    fflush(stdout);
    if(doFourier){
      sprintf(tmp,"%s.fft",name);
      fourier(tmp,h[3],mask,nr[3],dT[3],0);
      }
/*     season("senetosa.dif",h[3],t[3],mask,nr[3],dT[3]); */

/* *-----------------------------------------------------------------------------
    save statistics */
    fprintf(statfile,"%10s %8.3lf %8.3lf %6d %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf\n", sample[k].name, sample[k].t, sample[k].p, n,
            stat[0].std,stat[1].std,stat[2].std, stat[0].variance(),stat[1].variance(),stat[2].variance());
    point_t point(k, sample[k].t, sample[k].p, 3);
    for(l=0;l<3;l++) point.z[l]=stat[l].variance();
    points.push_back(point);
    
//    for (k=0;k<3;k++) gstat[k].push_back(stat[k]);
    
/* *-----------------------------------------------------------------------------
    build the common time serie: model divided by IB */
//     out=fopen("ibtest","w");
//     for (i=0; i< nr[0]; i++) {
//       if((h[1][i]==mask)||(h[0][i] ==mask)) {
// /*         printf("%f %f \n",h[0][i],h[1][i]); */
//         }
//       if((fabs(h[1][i]) > 0.01)&&(fabs(h[0][i]) > 0.01)) {
//         h[0][i]=h[1][i]/h[0][i];
//         fprintf(out,"%f\n",h[0][i]);
//         }
//       else {
//         h[0][i]=mask;
//         }
//       }
//     fclose(out);
//     printf("model ouput:difference with IB correction----------------\n");
// /*    filtering("senetosa.dif",h[0],t[0],mask,nr[0],dT[0]);*/
//     strcpy(name,model);
//     strcat(name,sample[k].name);
//     strcat(name,".adm");
//     sdum=filtering(name,h[0],t[0],mask,nr[0],dT[0]);
    for (i=0;i<4;i++){
      delete[] h[i];
      delete[] t[i];
      }
    }
  

  legend_t *legends=lgd_import(points);
  status=lgd_save("statistics.lgd", legends, 1, NULL, NULL);
  

  fclose(statfile);
  TRAP_ERR_EXIT(0,"%s -computation complete ^^^^^^^^^\n",argv[0]);
}
