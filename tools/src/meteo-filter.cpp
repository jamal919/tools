
/*******************************************************************************

  T-UGO tools, 2006-2013

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Filter meteo forcing.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>

#include "config.h"

#include "functions.h" /* for print_OPENMP_help() */
#include "poc-assertions.h" /* for wexit() */
#include "poc-time.h" /* for date_t */

#include "lanczos.h"

#include "mgr-converter.h" /* for mgr_save_timeserie() */

#include "datastream.h"
#include "poc-netcdf-data.hpp" /* for poc_get_grid() */
#include "tides.h" /* for d2s */

#include "poc-grib.h"


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
    "  Filter meteo forcing.\n"
#ifndef HAVE_LIBGRIB_API
#warning WILL NOT WORK : see below
    "  THIS WILL NOT WORK BECAUSE YOU HAVE NOT COMPILED THE TOOLS WITH GRIB-API!\n"
#endif
    "\n"
    "OPTIONS\n"
    "If a default is not specified, the option is mandatory.\n"
    "  -h,--help : show this help and exit\n"
    "  -p : followed by archive directory. Default: ./\n"
    "  -c : followed by archive file name convention. See CONVENTION below. Default: EA_YYYYMM.EC\n"
    "  -s : followed by the start date. See DATE FORMATS below\n"
    "  -f : followed by the end date. See DATE FORMATS below\n"
    );
  print_poctime_scan_date_help(1);
  print_stream_DecodeName_help();
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
  parse arguments */
  int i,k,n,status;
  const char *keyword,*path=".",*convention="EA_YYYYMM.EC";
  date_t startDate=NADate,endDate=NADate;
  char **varnames=0;
  int nvars=0;
  
  if(argc<=1) {
    print_help(argv[0]);
    wexit(-1);
    }
  
  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=argv[n];
    if( strcmp(keyword,"--help")==0 or
        strcmp(keyword,"-h")==0 ) {
      print_help(argv[0]);
      wexit(0);
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 'p' :
          path= argv[n+1];
          n++;
          n++;
          break;

        case 'c' :
          convention= argv[n+1];
          n++;
          n++;
          break;

        case 's' :
          status=poctime_scan_date(argv[n+1],&startDate,&endDate);
          n++;
          n++;
          break;

        case 'f' :
          status=poctime_scan_date(argv[n+1],0,&endDate);
          n++;
          n++;
          break;

        case 'v' :
          for(k=1;n+k<argc && argv[n+k][0]!='-';k++);
          nvars=k-1;
          if(varnames){
            __FILE_LINE__(stdout,"multiple use of -v : the last one overrides the previous.");
            for(k=0;varnames[k]!=NULL;k++){
              free(varnames[k]);
              }
            delete[] varnames;
            }
          varnames=new char*[nvars+1];
          for(k=0;k<nvars;k++){
            varnames[k]=strdup(argv[n+1+k]);
            }
          varnames[nvars]=NULL;
          n++;
          n+=nvars;
          break;

        default:
          printf("unknown option %s\n",keyword);
          print_help(argv[0]);
          wexit(-1);
        }
        break;

      default:
        printf("unknown option %s\n",keyword);
        print_help(argv[0]);
        wexit(-1);
      }
    
    }
  
  bool do_print_help=false;
  printf("#number of variables : %d\n",nvars);
  if(nvars==0) {
    printf("*** Please give the names of the variables to filter with -v ***\n");
    do_print_help=true;
    }
  
  char *startS,*endS;
  startS=sgetdate(startDate);
  endS  =sgetdate(  endDate);
  printf("About to process [%s %s]:\n",startS,endS);
  delete[]startS;delete[]endS;
  if(isnad(startDate)) {
    printf("*** Please give a valid start date with -s ***\n");
    do_print_help=true;
    }
  else if(isnad(endDate)) {
    printf("*** Please give a valid end date with -f ***\n");
    do_print_help=true;
    }
  
  if(do_print_help){
    print_help(argv[0]);
    wexit(-1);
    }
  
/*------------------------------------------------------------------------------
  initialise sampling */
  double cut,sampling;
#if HarmonicFilterSpectrumType == 3
  cut=3./24.;
#else
  cut=6./24.;
#endif
  /* This is :
    + right for 1h ERA5
    + wrong for all others, that are >=3h and therefore can at best hardly be filtered anyway */
  sampling=1./24;
  
/*------------------------------------------------------------------------------
  initialise filter */
  double *filter,*LF,*HF,*times,*HFC;
  const int n_2=(3+cut)/sampling;
  
  n=n_2*2+1;
  
  times=new double[n];
  LF=aset(n,0.);
  HF=new double[n];
  HFC=new double[n];
  
  /* impulsion */
  filter=aset(n,0.);
  filter[n/2]=1.;
  
  /* low frequency part */
  double *weightLF,half;
  int nWeightLF;
  
  half=cut/2.0;
  lanczos1D_init(sampling,n,half,&weightLF,&nWeightLF);
  
  Loess1D_BF_nomask(filter,n,weightLF,nWeightLF,LF);
  
  /* LF => HF */
  for(i=0;i<n;i++)
    HF[i]=filter[i]-LF[i];
  
  /* remove harmonics from HF */
  for(i=0;i<n;i++)
    times[i]=sampling*i;
  
  status=HarmonicFilter(times, HF, HFC, NAN, n);
  
  /* finalise filter */
  for(i=0;i<n;i++){
    filter[i]=LF[i]+HF[i]-HFC[i];
    }
  
#if 0
  tseries_t saved(4);
  saved.n=n;
  saved.t=times;
  saved.x[0]=filter;
  saved.x[1]=LF;
  saved.x[2]=HF;
  saved.x[3]=HFC;
  status=mgr_init_formats();
  status=mgr_save_timeserie("filter" TOSTRING(HarmonicFilterSpectrumType), saved.mooring, saved, 'm', FORMAT_NAME_GNU, true);
  TRAP_ERR_EXIT(ENOEXEC,"not finished\n");
#endif
  
  delete[]times;
  delete[]LF;
  delete[]HF;
  delete[]HFC;
  
/*------------------------------------------------------------------------------
  initialise sequences and grib messages */
#ifdef HAVE_LIBGRIB_API
  const double startD=cnes_time(startDate,'s');
  const double   endD=cnes_time(  endDate,'s');
  const double dt=d2s*sampling;
  const double filter_dt=n_2*dt;
  double time,time0;
  
  sequence_t< SGfield_t<float> > *sequences;
  poc_dim_t dim;
  size_t s,smax=0;
  
  poc_global_t *global=0;
  FILE *F;
  grib_handle **Hs;
  
  sequences=new sequence_t< SGfield_t<float> >[nvars];
  Hs=new grib_handle*[nvars];
  
  time0=startD-filter_dt;
  
  for(i=0;i<nvars;i++) {
    sequence_t< SGfield_t<float> > *sequencei=&sequences[i];
    const char *varnamei=varnames[i];
    
    printf("#variable %d: %s\n",i,varnamei);
    
    /* initialise file stream */
    filestream_t<float> *filestream=new filestream_t<float>(path, convention, varnamei,0);
    status=filestream->finalize(time0);
    if(status != 0) NC_TRAP_ERROR(wexit,status,1,"filestream->finalize(%g) error",time0);
    
    /* initialise grid of frames */
    grid_t *grid=new grid_t;
    status=poc_get_grid(filestream->filename,varnamei,grid);
    if(status != 0) NC_TRAP_ERROR(wexit,status,1,"poc_get_grid(\""+filestream->filename+"\",%s) error",varnamei);
    
    /* initialise sequence */
    sequencei->allocate(n);
    
    for(k=0;k<sequencei->nframes;k++) {
      SGfield_t<float> *field=new SGfield_t<float>(filestream,grid,stream_SGxNETCDF);
      sequencei->frames[k]=*field;
      }
    
    /* initialise grib message */
    if(global==0)global=&filestream->glob;
    poc_var_t *var=global->variables.findP(varnamei);
    poc_att_t *att=var->attributes.findP(POC_GRIB_FILE_POS_VATT_NAME);
    s=(*att)[0];
    
    F=fopen(filestream->filename.c_str(),"r");
    if(F==0)TRAP_ERR_EXIT(errno,"fopen(\""+filestream->filename+"\",[r]) error (%d %s)\n",errno,strerror(errno));
    fseek(F,s,SEEK_SET);
    
    Hs[i]=grib_handle_new_from_file(0, F, &status);
    
    fclose(F);
    F=0;
    
    dim=var->dimensions.back();
    s=dim.len;
    smax=max(smax,s);
    }
  
/*------------------------------------------------------------------------------
  filter */
  double *buffer;
  struct timeval walltime;
  string outputPath;
  date_t frameDate;
  char *frameS,*tmp;
  long int param;
  const void *Hbuf;
  size_t Hbufsize;
  
  buffer=new double[smax];
  
  for(time=startD;time<endD;time+=dt){
    
    stream_DecodeName(time, "", "filtered-EA_YYYYMM.EC", &tmp);
    
    if(outputPath!=tmp or F==0){
      if(F!=0)fclose(F);
      outputPath=tmp;
      F=fopen(outputPath.c_str(),"w");
      }
    
    delete[]tmp;
    
    frameDate=poctime_getdatecnes(time,'s');
    frameS=sgetdate(frameDate);
    printf("Frame %s in "+outputPath+"\n",frameS);fflush(stdout);
    delete[]frameS;
    
    for(i=0;i<nvars;i++) {
      sequence_t< SGfield_t<float> > *sequencei=&sequences[i];
      const char *varnamei=varnames[i];
      
      poc_var_t *var=global->variables.findP(varnamei);
      dim=var->dimensions.back();
      s=dim.len;
      
      printf("#variable %d: %s (%u values). ",i,varnamei,s);fflush(stdout);
      gettimeofday(&walltime);
      
      time0=time-filter_dt;
      
      /* update sequence */
      status=sequencei->init(time0);
      if(status != 0) TRAP_ERR_EXIT(status,"sequence_t<>::init(%g) error\n",time0);
      
      printf("done in %g s. ",difftime(&walltime));
      
      /* filter */
      printf("Filtering... ");fflush(stdout);
      
      aset(buffer,s,0.);
      
// #warning comment this line and the one below
//       if(false)
      for(k=0;k<n;k++){
        const float *xk=sequencei->frames[k].x;
        const double w=filter[k];
        
        //#pragma omp parallel for /* slower with openMP */
        for(int m=0;m<s;m++){
          buffer[m]+=w*xk[m];
          }
        
        }
      
      printf("done in %g s. ",difftime(&walltime));
      
      /* output */
      printf("Writing "+outputPath+" ... ");fflush(stdout);
      
      const SGfield_t<float> *framen=&sequencei->frames[n_2];
//       filestream_t<float> *filestream=(filestream_t<float> *)framen->stream; /* same for all frames */
      
      grib_handle *Hi=Hs[i];
      
      if(framen->time!=time) TRAP_ERR_EXIT(ENOEXEC,"programming error : %g!=%g\n",framen->time,time);
      
      param=(frameDate.year*100+frameDate.month)*100+frameDate.day;
      status=grib_set_long(Hi,"dataDate",param);
      param=frameDate.second/3600;
      status=grib_set_long(Hi,"hour",param);
      status=grib_set_double_array(Hi,"values",buffer,s);
      
      status=grib_get_message(Hi,&Hbuf,&Hbufsize);
      status=fwrite(Hbuf,1,Hbufsize,F);
      if(status!=Hbufsize) TRAP_ERR_EXIT(errno,"fwrite(,1,,(\""+outputPath+"\")) error (%d %s)\n",errno,strerror(errno));
      
      printf("done in %g s.\n",difftime(&walltime));
      
      /* rotate frames */
      sequencei->rotate_frames();
      }
    
    }
  
  if(F!=0)fclose(F);
#else
  TRAP_ERR_EXIT(ENOEXEC,"Not compiled with GRIB-API\n");
#endif
  
  return 0;
}
