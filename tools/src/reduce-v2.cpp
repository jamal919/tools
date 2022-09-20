
/*******************************************************************************

  T-UGO tools, 2006-2013

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Filter Mog2D archives. USED BY CTOH.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "tools-structures.h"

#include "tides.h"
#include "poc-netcdf-data.hpp"
#include "tides.def"
#include "functions.h"
#include "fe.h"
#include "poc-time.h"
#include "archive.h"
#include "filter.h"

#define NONE       -1
#define HARMONIC    0
#define CLIMATOLOGY 1


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int defineCurrentDate(int index, int month,int year, int *currentMonth, int *currentYear)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  if( (month < 1) || ( month > 12) ) return(1);

  switch(index){
    case 0: // define previous month
      if(month>1) {
        *currentMonth=month-1;
        *currentYear=year;
        return(0);
        }
      else {
        *currentMonth=12;
        *currentYear=year-1;
        return(0);
        }
      break;

    case 1:// define curent month
      *currentMonth=month;
      *currentYear=year;
      return(0);
      break;

    case 2:// define next month
      if(month<12) {
        *currentMonth=month+1;
        *currentYear=year;
        return(0);
        }
      else {
        *currentMonth=1;
        *currentYear=year+1;
        return(0);
        }
      break;

    default:
      return(1);
      break;

  }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : argv[0]
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    "NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Filter Mog2D archives. USED BY CTOH.\n"
    "Compute and write clx-formatted filtered Mog2D archives.\n"
    "\n"
    "OPTIONS\n"
    "If a default is not specified, the option is mandatory.\n"
    "  -h,--help : show this help and exit\n"
    "  -start : followed by starting month for filtering request.\n"
    "  -end : followed by ending month for filtering request.\n"
    "  --pressure-only : skip elevation\n"
    "  -h : followed by root name for elevation archives. Default : analysis\n"
    "  -p : followed by root name for pressure archives. Default : forcing\n"
    "  -r : followed by ratio factor for archive reduction. Default : 1\n"
    "  -f : followed by input archive format: clx or cdf\n"
    "  -oF : followed by output archive format: clx or cdf\n"
    "  -i : followed by path name for archive. Default : ./\n"
    "  -o : followed by path name for filtered archive. Default : ./\n"
    "  -no_default : disable default values for -h and -p\n"
    "  -pressure_factor : followed by pressure factor. Default: 1.0\n"
    "  -d : OBSOLETE! Use -m instead\n"
    "  -m : followed by deting method: analysis, climatology or none. Default: analysis\n"
    "    The climatology method currently requires S2-elevation-atlas.nc or S2.ele.s2c\n"
    "  -w : followed by filtering window, in days. Default: 20.0\n"
    "  -hf : also remove HF\n"
    );
  print_OPENMP_help(prog_name);
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,reI,k;
  int n,status;
  int year,month;
  int currentMonth,currentYear;
  int nndes;
  int start_year,last_year,start_month,last_month;
  int skip=0,nframes=0,nbFrames[3];
  int nWeightBF,nWeightHF;
  int resampling=1;

  double *htime=NULL;
  float **h=NULL,**buffer=NULL;
  float sampling=0;
  float *weightBF=NULL,*weightHF=NULL;
  bool removeHF=false;
//   float mask=1.e+10;
  float half;
  float delta;
  float tStart=0,tEnd=0;

  double time;
/* *----------------------------------------------------------------------------
  filter window, default is Topex*/
  double window=20.0;
  char *h_root=NULL,*p_root=NULL,*pathname=NULL,*format=NULL;
  int inputFormat,outputFormat=TUGO_UG_BINARY_FORMAT;
  char *output=NULL,*default_dir=NULL;
  char *hfile=NULL,*pfile=NULL,*method=NULL;
  bool skipElevation=false;
  const char *keyword;

//   int detide_s2=0,detide_k2=0,detide_s1=0,detide_sa=0;
  int detiding=HARMONIC;

  int default_root_allowed=1;
  fcomplex **constants=NULL,*cbuffer=NULL;
  char  tidefile[1024]="\0";
  float factor=1.0;

  date_t date;
  meta_archive_t info,filtered;
  poc_var_t var;
  spectrum_t spectrum;

  fct_echo( argc, argv);

/*------------------------------------------------------------------------------
  quick how-to of the program  */
  if(argc==1){
    print_help(argv[0]);
    return(0);
    }

  /*---------------------------------------------------------------------
    test input arguments */

  n=1;
  while (n < argc) {
    keyword=argv[n];
    
    if( ( n==argc-1 and
          strcmp(keyword,"-h")==0
          ) or
        strcmp(keyword,"--help")==0 ){
      print_help(argv[0]);
      return(0);
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

    if(strcmp(keyword,"-no_default") == 0) {
      default_root_allowed=0;
      n++;
      continue;
      }

    if(strcmp(keyword,"--pressure-only") == 0) {
      skipElevation=true;
      n++;
      continue;
      }

    if(strcmp(keyword,"-oF") == 0) {
      if(strcmp(argv[n+1],"cdf")==0)
        outputFormat=TUGO_UG_NETCDF_FORMAT;
      n++;
      n++;
      continue;
      }

    if(strcmp(keyword,"-pressure_factor") == 0) {
      status=sscanf(argv[n+1],"%f",&factor);
      n++;
      n++;
      continue;
      }

    if(strcmp(keyword,"-hf") == 0) {
      removeHF=true;
      n++;
      continue;
      }

    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
/*------------------------------------------------------------------------------
          detiding method */
          case 'd' :
            method= strdup(argv[n+1]);
            n++;
            n++;
            STDOUT_BASE_LINE("obsolete key, use -m\n");
            exit(-1);
            break;

/*------------------------------------------------------------------------------
          detiding method */
          case 'm' :
            method= strdup(argv[n+1]);
            n++;
            n++;
            break;

/*------------------------------------------------------------------------------
          root name for default ocean archive */
          case 'h' :
            h_root= strdup(argv[n+1]);
            n++;
            n++;
            break;

/*------------------------------------------------------------------------------
          root name for default forcing archive */
          case 'p' :
            p_root= strdup(argv[n+1]);
            n++;
            n++;
            break;

/*------------------------------------------------------------------------------
          pathname for default output directory */
          case 'o' :
            default_dir= strdup(argv[n+1]);
            n++;
            n++;
            break;

/*------------------------------------------------------------------------------
          pathname for default archive directory */
          case 'i' :
            pathname= strdup(argv[n+1]);
            n++;
            n++;
            break;

/*------------------------------------------------------------------------------
          format option  */
          case 'f' :
            format= strdup(argv[n+1]);
            n++;
            n++;
            break;

/*------------------------------------------------------------------------------
          resampling option  */
          case 'r' :
            sscanf(argv[n+1],"%d",&resampling);
            n++;
            n++;
            break;

/*------------------------------------------------------------------------------
          filter window */
          case 'w' :
            sscanf(argv[n+1],"%lf",&window);
            n++;
            n++;
            break;

          default:
            STDOUT_BASE_LINE("unknown option %s\n",keyword);
            exit(-1);
          }
        break;

      default:
        n++;
        break;
      }
    
    }

  fprintf(stderr,"\n %s -starting filtering T-UGO archives ^^^^^^^^^\n",argv[0]);

/*------------------------------------------------------------------------------
  default initialization  */
  if(method==NULL) {
    method=strdup("analysis");
    printf("use harmonic analysis default for detiding method\n");
    }

  if(default_dir==NULL) {
    default_dir=strdup("./");
    printf("use <./> as root name for default output directorys\n");
    }
  
  if(format==0){
    printf("*** format (option -f) not given ***\n");
    print_help(argv[0]);
    exit(-1);
    }
  else if(strcmp(format,"cdf") == 0) {
    inputFormat=TUGO_UG_NETCDF_FORMAT;
    }
  else if(strcmp(format,"clx") == 0) {
    inputFormat=TUGO_UG_BINARY_FORMAT;
    }
  else{
    printf("*** unknown format (option -f) %s ***\n",format);
    print_help(argv[0]);
    exit(-1);
    }
  
  if(default_root_allowed) {
    if(h_root==NULL) {
      h_root= strdup("analysis");
      }
    printf("use %s as root name for sea state archive files\n",h_root);
    if(p_root==NULL) {
      p_root= strdup("forcing");
      }
    printf("use %s as root name for forcing archive files\n",p_root);
    }

  if(pathname==NULL) {
    pathname=(char *)  strdup("./");
    printf("use <./> as default path for input directory\n");
    }
  else {
    printf("use %s as path for input directory\n",pathname );
    }
  

  if( (h_root==NULL) && (p_root==NULL) ) {
    printf("*** -pressure_factor given with neiter -h nor -p  ***\n",resampling);
    print_help(argv[0]);
    exit(-1);
    }

/*------------------------------------------------------------------------------
  Basic checks */
  if(resampling<0) {
    printf("*** ratio factor (option -r) %d<0  ***\n",resampling);
    print_help(argv[0]);
    exit(-1);
    }

  date.year=start_year;
  date.month=start_month;
  date.day=1;
  date.second=0;
  tStart=(float)cnes_time(date,'d');

  date.year=last_year;
  date.month=last_month;
  date.day=1;
  date.second=0;
  tEnd=(float) cnes_time(date,'d');

  if(tStart>tEnd)
    TRAP_ERR_EXIT(-1,"start date %s after end date %s\n",
      poctime_sdate_cnes(tStart,'d'),
      poctime_sdate_cnes(tEnd  ,'d'));

  if(strcmp(method,"climatology")==0) {
    detiding=CLIMATOLOGY;
    }
  else if(strcmp(method,"analysis")==0) {
    detiding=HARMONIC;
    }
  else if(strcmp(method,"none")==0) {
    detiding=NONE;
    }
  else {
    printf("*** deting method (option -m) must be none, analysis or climatology, not %s ***\n",method);
    print_help(argv[0]);
    exit(-1);
    }
  printf("deting method: %s\n",method);

//   switch (detiding) {
//     case CLIMATOLOGY:
//       if(detide_s2==1) {
//         status=tide_addwave(&spectrum, wS2);
//        }
//       if(detide_k2==1) {
//         status=tide_addwave(&spectrum, wK2);
//         }
//       if(detide_s1==1) {
//         status=tide_addwave(&spectrum, wS1);
//         }
//       if(detide_sa==1) {
//         status=tide_addwave(&spectrum, wSa);
//         }
//       break;
//     }


/*------------------------------------------------------------------------------
    Purposes:

    remove S2 signal
    filter signal by keeping 12.5 hours < period < 20 days for h,u and v
    filter signal by keeping 20 days  < period for ib
    leave wind stress unchanged

    then resample archive at N hours sampling rate (N is given in input)

---------------------------------------------------------------------*/

  exitIfNull(
    hfile=(char*)   malloc(1024)
    );
  exitIfNull(
    pfile=(char*)   malloc(1024)
    );
  exitIfNull(
    output=(char *) malloc(1024)
    );
  
  char separator;
  if(inputFormat==TUGO_UG_NETCDF_FORMAT) {
    if(h_root!=NULL) {
      if(strcmp(h_root,"analysis")==0)
        separator='-';
      else
        separator='.';
      sprintf(hfile,"%s/%s%c%4d.%02d.nc",pathname,h_root,separator,start_year,start_month);
      strcpy(pfile,hfile);
      }
    }
  else if(inputFormat==TUGO_UG_BINARY_FORMAT) {
    if(h_root!=NULL) sprintf(hfile,"%s/%s-%4d.%02d"   ,pathname,h_root,start_year,start_month);
    if(p_root!=NULL) sprintf(pfile,"%s/%s-%4d.%02d"   ,pathname,p_root,start_year,start_month);
    }
  else TRAP_ERR_EXIT(ENOEXEC,"unknown format code %d\n",inputFormat);

  if(h_root!=NULL) {
    status=archive_info(hfile,inputFormat,&info);
    if(status!=0) TRAP_ERR_EXIT(status,"archive_info(\"%s\",%d,) error (status=%d)\n",hfile,inputFormat,status);
    }
  else if(p_root!=NULL){
    status=archive_info(pfile,inputFormat,&info);
    if(status!=0) TRAP_ERR_EXIT(status,"archive_info(\"%s\",%d,) error (status=%d)\n",pfile,inputFormat,status);
    }

  nndes=info.mesh.nvtxs;
  nframes=info.nframe+4*24;// ceil nframes/month for allocations
  sampling=info.sampling;
  archive_freeinfo(&info);

  exitIfNull(
    buffer=(float **) malloc(3*sizeof(float *))
    );
  for(i=0;i<3;i++)
    exitIfNull(
      buffer[i]=(float *)malloc(nndes*sizeof(float))
      );

/* *---------------------------------------------------------------------
  initialize window x days loess filter*/

  half=window/2.0*d2s;
  lanczos1D_init(sampling,3*nframes,half,&weightBF,&nWeightBF);
  if(nWeightBF==-1) TRAP_ERR_EXIT(nWeightBF,"lanczos1D_init error (nWeightBF=%d)\n",nWeightBF);

/*------------------------------------------------------------------------------
  init 12.5 hours loess filter*/
  if(removeHF){
    half=6.25*3600.;
    lanczos1D_init(sampling,3*nframes,half,&weightHF,&nWeightHF);
    if(nWeightBF==-1) TRAP_ERR_EXIT(nWeightBF,"lanczos1D_init error (nWeightBF=%d)\n",nWeightBF);
    }


  exitIfNull(
    h=(float **) malloc(nndes*sizeof(float *))
    );
  for(n=0;n<nndes;n++)
    exitIfNull(
      h[n]=(float *) malloc(3*nframes*sizeof(float))
      );

  exitIfNull(
    htime=(double *) malloc(3*nframes*sizeof(double))
    );

/* *---------------------------------------------------------------------
  hardcoded, init tidal spectrum*/
  spectrum.n=1;
  spectrum.waves=new tidal_wave[spectrum.n];
  i=0;
  spectrum.waves[i++] = wS2;

  for (i=0; i<spectrum.n; i++) {
    spectrum.waves[i].init();
    if(h_root!=NULL) printf ("remove wave: %10s, pulsation: %12.6f degrees/h \n", spectrum.waves[i].name,spectrum.waves[i].omega);
    }

/* *------------------------------------------------------------------------
  !!! not secure if consecutive archive files are from different FE mesh */
  if(detiding==CLIMATOLOGY) {
    if((constants==NULL) && (spectrum.n !=0)) {
      status=initialize_omega(&spectrum);
      constants=new fcomplex*[nndes];
      for(n=0;n<nndes;n++) {
        constants[n]=new fcomplex[spectrum.n];
        }
      cbuffer=new fcomplex[nndes];
      for(k=0;k<spectrum.n;k++) {
        sprintf(tidefile,"%s-elevation-atlas.nc",spectrum.waves[k].name);
        printf("Trying climatology from %s\n",tidefile);fflush(stdout);
        status=poc_get_cvara(tidefile,"elevation_a","elevation_G",0,cbuffer);
        if(status!=0){
          NC_CHKERR_BASE_LINE(status,"poc_get_cvara error");
          sprintf(tidefile,"%s.ele.s2c",spectrum.waves[k].name);
          printf("Trying climatology from %s\n",tidefile);fflush(stdout);
          status=quoddy_loadc1(tidefile, nndes,cbuffer);
          }
        if(status!=0) TRAP_ERR_EXIT(status,"quoddy_loadc1(\"%s\",,) error (%d %s)\n",tidefile,status,strerror(status));
        for(n=0;n<nndes;n++) {
          constants[n][k]=cbuffer[n];
          }
        }
      }
    }
/** !!! not secure if consecutive archive files are from different FE mesh
------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
  start temporal loop*/
  for (year=start_year;year<=last_year;year++) {
    for (month=1;month<=12;month++) {
      if((year==start_year) && (month < start_month)) continue;
      if((year==last_year)  && (month > last_month))  break;
      if(skipElevation) goto pressure_only;
/*------------------------------------------------------------------------------
      running window centered on current month*/
      skip=0;
      for (k=0;k<3;k++) {
        status=defineCurrentDate(k,month,year,&currentMonth,&currentYear);
        if(status!=0) TRAP_ERR_EXIT(status,"defineCurrentDate error (status=%d)\n",status);
        if(inputFormat==TUGO_UG_NETCDF_FORMAT) {
          sprintf(hfile,"%s/%s%c%4d.%02d.nc",pathname,h_root,separator,currentYear,currentMonth);
          strcpy(pfile,hfile);
          }
        else if(inputFormat==TUGO_UG_BINARY_FORMAT) {
          sprintf(hfile,"%s/%s-%4d.%02d",pathname,h_root,currentYear,currentMonth);
          sprintf(pfile,"%s/%s-%4d.%02d",pathname,p_root,currentYear,currentMonth);
          }
        else TRAP_ERR_EXIT(ENOEXEC,"unknown format code %d\n",inputFormat);

/*------------------------------------------------------------------------------
        treat elevation*/
        status=archive_info(hfile,inputFormat,&info);
        if(status!=0) TRAP_ERR_EXIT(status,"archive_info(\"%s\",%d,) error (status=%d)\n",hfile,inputFormat,status);
        if(k==1){
          status=archive_info(hfile,inputFormat,&filtered);
          if(status!=0) TRAP_ERR_EXIT(status,"archive_info(\"%s\",%d,) error (status=%d)\n",hfile,inputFormat,status);
          }
        nbFrames[k]=info.nframe;
        printf("read %d frames elevation archives from %s\n",info.nframe,hfile);
        delta=(float) cnes_time(info.reference,'d');

        info.flag=FLAG_STATE;
        for (i=0;i<nbFrames[k];i++) {
          status=archive_read(hfile,inputFormat,info,i,buffer,0,&time);
          if(status!=0) TRAP_ERR_EXIT(status,"archive_read error (status=%d)\n",status);
          if(i%24==0){
            printf("t= %6.1f hours status= %d %s\n",time/3600.0,status, sgetnewdate(info.reference,time));
            fflush(stdout);
            }
          for(n=0;n<nndes;n++) h[n][i+skip]=buffer[0][n];
          htime[i+skip]=(delta+time/d2s);// t in CNES days for harmonic analysis
          }

        skip+=nbFrames[k];
        archive_freeinfo(&info);
        }

      for (k=0,nframes=0;k<3;k++) nframes+=nbFrames[k]; // total nb of frames during 3 months

      skip=nbFrames[0];
      
#define DO_FILTER 1
      fprintf(stderr,"filter elevation archives\n");
#pragma omp parallel
      {
      int n,i;
      float *residuals=new float[nframes];
      float *tide=0,mean;
      
      if(detiding!=NONE)
        tide=new float[nframes];
      
#if DO_FILTER
      float *lf=NULL,*hf=NULL;
      
      lf=new float[nframes];
      if(removeHF)
        hf=new float[nframes];
#endif
      
#pragma omp for
      for(n=0;n<nndes;n++) {
        
        if(detiding==HARMONIC) {
/* *---------------------------------------------------------------------
          remove S2 by harmonic analysis */
          status=harmonic_analysis_old(h[n],htime,residuals,&mean,nframes,spectrum);
          if(status!=0) TRAP_ERR_EXIT(status,"harmonic_analysis_old error (status=%d)\n",status);
          for(i=0;i<nframes;i++) residuals[i]+=mean;
          }
        else if(detiding==CLIMATOLOGY) {
/* *---------------------------------------------------------------------
          remove S2 by climatoloy prediction*/
          status=harmonic_predictionTS(tide, htime, nframes, spectrum, constants[n], 0);
          if(status!=0) TRAP_ERR_EXIT(status,"harmonic_predictionTS error (status=%d)\n",status);
          for(i=0;i<nframes;i++) residuals[i]=h[n][i]-tide[i];
          }
        else if(detiding==NONE) {
#if 1
          for(i=0;i<nframes;i++) residuals[i]=h[n][i];
#endif
          }
        
#if DO_FILTER
        
        if(removeHF){
/* *---------------------------------------------------------------------

  additional filtering, to be reworked and tested (Demerliac?)

---------------------------------------------------------------------**/
/* *---------------------------------------------------------------------
          first eliminate elevation time series high frequency (~12h)*/
          //Loess1D_BF(residuals,mask,nframes,weightHF,nWeightHF,lf);
          Loess1D_BF_nomask(residuals,nframes,weightHF,nWeightHF,lf);
/* *---------------------------------------------------------------------
          second, eliminate elevation time series low frequency (20 days)*/
          //Loess1D_BF(lf,mask,nframes,weightBF,nWeightBF,hf);
          Loess1D_BF_nomask(lf,nframes,weightBF,nWeightBF,hf);
          for(i=0;i<nbFrames[1];i++) h[n][skip+i]=lf[skip+i]-hf[skip+i];
          }
        else{
/* *---------------------------------------------------------------------

  usual treatment

---------------------------------------------------------------------**/
/* *---------------------------------------------------------------------
          eliminate elevation time series low frequency (20 days)*/
          //Loess1D_BF(residuals,mask,nframes,weightBF,nWeightBF,lf);
          Loess1D_BF_nomask(residuals,nframes,weightBF,nWeightBF,lf);
          for(i=0;i<nbFrames[1];i++) h[n][skip+i]=residuals[skip+i]-lf[skip+i];
          }
        
#else
        for(i=0;i<nbFrames[1];i++) h[n][skip+i]=residuals[skip+i];
#endif
        }
      
      if(detiding!=NONE)
        delete[]tide;
#if DO_FILTER
      deletep(&hf);
      deletep(&lf);
#endif
      delete[]residuals;
      }

/*------------------------------------------------------------------------------
      write filtered current month archive*/
      if(outputFormat==TUGO_UG_NETCDF_FORMAT) {
        if (resampling>1) sprintf(output,"%s/%s-filt-%d-%4d.%02d.nc",default_dir,h_root,resampling,year,month);
        else sprintf(output,"%s/%s-filt-%4d.%02d.nc",default_dir,h_root,year,month);
        
        var.init("elevation",NC_FLOAT,"filtered elevation","m");
        }
      else if(outputFormat==TUGO_UG_BINARY_FORMAT) {
        if (resampling>1) sprintf(output,"%s/%s-filt-%d-%4d.%02d",default_dir,h_root,resampling,year,month);
        else sprintf(output,"%s/%s-filt-%4d.%02d",default_dir,h_root,year,month);
        }
      else TRAP_ERR_EXIT(ENOEXEC,"unknown format code %d\n",inputFormat);
      printf("write %d frames elevation archives to %s\n",nbFrames[1],output);

      filtered.sampling*=resampling;
      filtered.nframe/=resampling;
      
      for(i=0,reI=0;i<nbFrames[1];i+=resampling,reI++) {
/* *---------------------------------------------------------------------
        caution: buffer[1] and buffer[2] are not filtered !!! */
        for(n=0;n<nndes;n++) buffer[0][n]=h[n][skip+i];
        time=(htime[skip+i]-delta)*d2s;
        
        if(i%24==0){
          printf("t= %6.1f hours status= %d %s\n",time/3600.0,status, sgetnewdate(filtered.reference,time));
          fflush(stdout);
          }
        
        if(outputFormat==TUGO_UG_NETCDF_FORMAT){
          status=poc_put_mesh_vara(output,filtered.mesh,LGP1,&var,reI,buffer[0],0);
          status=poc_put_vara(output,"time",reI,&time);
          continue;
          }
        
        if(i==0) {
          status=clx_archivewriteheader(output,filtered,buffer);
          if(status!=0) TRAP_ERR_EXIT(status,"clx_archivewriteheader error (status=%d)\n",status);
          }
        else {
          status=archive_update2(output, filtered, buffer, time);
          if(status!=0) TRAP_ERR_EXIT(status,"archive_update2 error (status=%d)\n",status);
          }
        }
      archive_freeinfo(&filtered);

      if(p_root==NULL) continue;

/* *---------------------------------------------------------------------
      treat pressure time series */
pressure_only:

      skip=0;
      for (k=0;k<3;k++) {
        status=defineCurrentDate(k,month,year,&currentMonth,&currentYear);
        if(status!=0) TRAP_ERR_EXIT(status,"defineCurrentDate error (status=%d)\n",status);
        
        if(inputFormat==TUGO_UG_NETCDF_FORMAT)
          sprintf(pfile,"%s/%s%c%4d.%02d.nc",pathname,h_root,separator,currentYear,currentMonth);
        else if(inputFormat==TUGO_UG_BINARY_FORMAT)
          sprintf(pfile,"%s/%s-%4d.%02d"   ,pathname,p_root,currentYear,currentMonth);
        else TRAP_ERR_EXIT(ENOEXEC,"unknown format code %d\n",inputFormat);
        
        status=archive_info(pfile,inputFormat,&info);
        if(status!=0) TRAP_ERR_EXIT(status,"archive_info(\"%s\",%d,) error (status=%d)\n",pfile,inputFormat,status);
        if(k==1){
          status=archive_info(pfile,inputFormat,&filtered);
          if(status!=0) TRAP_ERR_EXIT(status,"archive_info(\"%s\",%d,) error (status=%d)\n",pfile,inputFormat,status);
          }
        nbFrames[k]=info.nframe;
        printf("read %d frames pressure archives from %s\n",info.nframe,pfile);
        delta=(float) cnes_time(info.reference,'d');

        info.flag=FLAG_FORCING;
        for (i=0;i<nbFrames[k];i++) {
          status=archive_read(pfile,inputFormat,info,i,buffer,0,&time);
          if(status!=0) TRAP_ERR_EXIT(status,"archive_read error (status=%d)\n",status);
          if(i%24==0){
            printf("t= %6.1f hours status= %d %s\n",time/3600.0,status, sgetnewdate(info.reference,time));
            fflush(stdout);
            }
          for(n=0;n<nndes;n++) h[n][i+skip]=buffer[0][n];
          htime[i+skip]=(delta+time/d2s);// t in CNES days
          }
        skip+=nbFrames[k];
        archive_freeinfo(&info);
        }

      for (k=0,nframes=0;k<3;k++) nframes+=nbFrames[k];// total nb of frames during 3 months

      skip=nbFrames[0];
      
/* *---------------------------------------------------------------------
      eliminate pressure time series low frequency*/
#if DO_FILTER
      fprintf(stderr,"filter pressure archives\n");
#pragma omp parallel
      {
      float *lf=NULL;
      int n,i;
      
      lf=new float[nframes];
#pragma omp for
      for(n=0;n<nndes;n++) {
//      Loess1D_BF(h[n],mask,nframes,weightBF,nWeightBF,lf);
        Loess1D_BF_nomask(h[n],nframes,weightBF,nWeightBF,lf);
        for(i=0;i<nbFrames[1];i++) h[n][skip+i]=lf[skip+i];
        }
      
      delete[]lf;
      }
#endif

/*------------------------------------------------------------------------------
      write filtered archive*/
      if(outputFormat==TUGO_UG_NETCDF_FORMAT) {
        if(resampling>1) sprintf(output,"%s/%s-filt-%d-%4d.%02d.nc",default_dir,h_root,resampling,year,month);
        else sprintf(output,"%s/%s-filt-%4d.%02d.nc",default_dir,h_root,year,month);
        
        var.init("pressure",NC_FLOAT,"filtered pressure","m");
        }
      else if(outputFormat==TUGO_UG_BINARY_FORMAT) {
        if(resampling>1) sprintf(output,"%s/%s-filt-%d-%4d.%02d",default_dir,p_root,resampling,year,month);
        else sprintf(output,"%s/%s-filt-%4d.%02d",default_dir,p_root,year,month);
        }
      else TRAP_ERR_EXIT(ENOEXEC,"unknown format code %d\n",outputFormat);
      printf("write %d frames pressure archives to %s\n",nbFrames[1],output);

      filtered.sampling*=resampling;
      filtered.nframe/=resampling;

      for(i=0,reI=0;i<nbFrames[1];i+=resampling,reI++) {
        for(n=0;n<nndes;n++) buffer[0][n]=factor*h[n][skip+i];
        time=(htime[skip+i]-delta)*d2s;
        
        if(i%24==0){
          printf("t= %6.1f hours status= %d %s\n",time/3600.0,status, sgetnewdate(filtered.reference,time));
          fflush(stdout);
          }
        
        if(outputFormat==TUGO_UG_NETCDF_FORMAT){
          status=poc_put_mesh_vara(output,filtered.mesh,LGP1,&var,reI,buffer[0],0);
          continue;
          }
        
        if(i==0) {
          status=clx_archivewriteheader(output,filtered,buffer);
          if(status!=0) TRAP_ERR_EXIT(status,"clx_archivewriteheader error (status=%d)\n",status);
          }
        else {
          status=archive_update2(output, filtered, buffer, time);
          if(status!=0) TRAP_ERR_EXIT(status,"archive_update2 error (status=%d)\n",status);
          }
        }
      archive_freeinfo(&filtered);
      } // end loop on month
    }// end loop on year

  /*---------------------------------------------------------------------
    free memory*/
  for(n=0;n<nndes;n++) free(h[n]);
  free(h);
  for(n=0;n<3;n++) free(buffer[n]);
  free(buffer);
  free(hfile);
  free(pfile);
  free(output);

  free(h_root);
  free(p_root);
  free(format);

  free(weightBF);
  free(weightHF);
  free(htime);

  fprintf(stderr,"%s -computation complete ^^^^^^^^^\n",argv[0]);
  return(0);
}
