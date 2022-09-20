
/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\date 2011-10-25 Damien ALLAIN : creation

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Corrects soundings with given atlases

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/


#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "tools-structures.h"

#include "tides.h"
#include "functions.h" //safely includes omp.h
#include "matrix.h"
#include "poc-time.h"
#include "filter.h"
#include "map.h"
//#include "tides.def"

#include "zapper.h"


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
    "USE\n"
    "  %s -p lon_lat_list -a atlas_convention [-s start -f end] -w wave1 [wave2 ... ] [OPTIONS] \n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Corrects soundings with given atlases.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  show this help\n"
    "  -p  followed by the path of the list of soundings: it is an ascii file with the year month day hour minute second latitude longitude sounding.\n"
    "  -a  followed by the atlas file name convention. See below.\n"
    "  -v  followed by the variables names for the amplitude and the phase. Default: Ha Hg\n"
    "  -w  followed by the list of waves to correct for\n"
    "  -o  followed by the path of the ascii output. Default: detided-soundings.dat\n"
    "\n"
    "TIP\n"
    "  To get all available atlases :\n"
    "f=(*.spectral.nc);%s -p soundings.dat -a WAVE.spectral.nc -v a_eta_LGP2 G_eta_LGP2 -w ${f[@]/.spectral.nc}\n",prog_name);
  print_tide_decode_atlasname_help();
    /** \endcode */
}

typedef struct {
  double time,lat,lon,inD,outD;//<time, lat and lon, input and output depth
  } epoch_t;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// function that will be ran when the executable is started
/** See the print_help <a href=#func-members>function</a> for its use.
*/
/*----------------------------------------------------------------------------*/
{
  int i,j,status;//<epoch and argument index, wave and variable index, argument index, status

  char *inP,*outP;//<paths of the input and the output files
  FILE *inF,*outF;
  
  vector<epoch_t> epochs;
  range_t<double> dateRange;
  
  double a,G;//<amplitudes and phases in the atlas at the point
  hconstant_t *constants;//<constants at the points
  
  double *buffer;
  double GA;//<Greenwich argument

  char *aPC=NULL,*aP;//<atlas path convention,atlas path
  char *waveNames[100];//<list of wave names
  spectrum_t WaveList;//<list of waves
  const int nodal_corrections=1;
  astro_angles_t astro_angles;
  
  char *keyword,*s;//<option and the argument that follows
  int verbose=1;
  
  const int nvars=2;
  char *varnames[nvars]={"Ha","Hg"};
  char *aN,*gN;//<decoded names of amplitudes and phases
  int defaultVars=1;
  
  string cmd;
  cmd=fct_echo( argc, argv);
  
  j=0;
  waveNames[j]=NULL;

  i=1;
  while (i < argc) {
    keyword=strdup(argv[i]);
    switch (keyword[0]) {
    case '-':
      if( strcmp(keyword,"-h")==0 or
          strcmp(keyword,"--help")==0 ){
        print_help(argv[0]);
        return 0;
        }
      switch (keyword[1]) {

      case 'p' :
        inP=strdup(argv[i+1]);
        i++;
        i++;
        break;

      case 'a' :
        aPC=strdup(argv[i+1]);
        i++;
        i++;
        break;

      case 'v' :
        if(varnames[0]){
          if(defaultVars)
            defaultVars=0;
          else{
            __FILE_LINE__(stdout,"multiple use of -v : the last one overrides the previous.");
            for(j=0;j<nvars;j++){
              free(varnames[j]);
              }
            }
          }
        for(j=0;j<nvars;j++){
          varnames[j]=strdup(argv[i+1+j]);
          }
        i++;
        i+=nvars;
        break;

      case 'w' :
        i++;
        j=pos((char*)NULL,waveNames,100);
        for(;i<argc;j++,i++) {
          waveNames[j]=strdup(argv[i]);
          }
        waveNames[j]=NULL;
        break;

      case 'o' :
        outP=strdup(argv[i+1]);
        i++;
        i++;
        break;

      default:
        printf("unknown option %s\i",keyword);
        print_help(argv[0]);
        exit(-1);
      }
      break;

      }
    free(keyword);
    }
  
  if(aPC==NULL){
    fprintf(stderr,"*** Please provide a convention for the name of the atlases. ***\n");
    print_help(argv[0]);
    exit(EINVAL);
    }

/*-----------------------------------------------------------------------------
  initialise : */
  
  //- list of waves
  if(waveNames[0]==NULL){
    fprintf(stderr,"*** Please provide a list of waves. ***\n");
    print_help(argv[0]);
    exit(EINVAL);
    }
  WaveList.init(initialize_tide(),waveNames);
  printf("# number of waves : %d \n",WaveList.n);

  //- input
  if(inP==NULL){
    fprintf(stderr,"*** Please provide an input file. ***\n");
    print_help(argv[0]);
    exit(EINVAL);
    }
  inF=fopen(inP,"r");
  if(inF==NULL){
    __OUT_BASE_LINE__("Can not open %s : %s",inP,strerror(errno));
    exit(errno);
    }
  
  //- output
  if(outP==NULL){
    outP=strdup("detided-soundings.dat");
    printf("Default output file : %s\n",outP);
    }
  outF=fopen(outP,"w");
  if(outF==NULL){
    __OUT_BASE_LINE__("Can not open %s : %s",outP,strerror(errno));
    exit(errno);
    }
  
/*-----------------------------------------------------------------------------
  read input : */
  
  for(i=0;true;i++){
    int read;
    date_t time;
    int hour,minute;
    epoch_t e;
    
    __ERR_BASE_LINE__("reading %d : %d;\r",i,read);
    
    //- time
    read=fscanf(inF,"%*[^0-9]%d",&time.year);if(read<1)break;
    read=fscanf(inF,"%*[^0-9]%d",&time.month);if(read<1)break;
    read=fscanf(inF,"%*[^0-9]%d",&time.day);if(read<1)break;
    read=fscanf(inF,"%*[^0-9]%d",&hour);if(read<1)break;
    read=fscanf(inF,"%*[^0-9]%d",&minute);if(read<1)break;
    read=fscanf(inF,"%*[^0-9]%g",&time.second);if(read<1)break;
    time.second+=((hour*60.)+minute)*60.;
    
    //- position
    e.lat=e.lon=NAN;
    read=fscanf(inF,"%*[^-0-9]%lg",&e.lat);if(read<1)break;
    read=fscanf(inF,"%*[^-0-9]%lg",&e.lon);if(read<1)break;
    
    //- sounding
    read=fscanf(inF,"%*[^-0-9]%lg",&e.inD);if(read<1)break;
    
    e.time=cnes_time(time,'s');
    
    epochs.push_back(e);
    dateRange<<e.time;
    }
  fclose(inF);
  const int n=i;//<number of epochs
  
  printf("# number of points found in %s : %d. Dates from %s to %s.\n",inP,n,
    poctime_sdate_cnes(dateRange.min,'s'),
    poctime_sdate_cnes(dateRange.max,'s')
    );
  
  {double *lat=new double[n];
  double *lon=new double[n];
  for(i=0;i<n;i++){
    lat[i]=epochs[i].lat;
    lon[i]=epochs[i].lon;
    }
  
/*-----------------------------------------------------------------------------
  read atlases */
  
  constants=tide_atlas2positions(aPC,WaveList,varnames[0],varnames[1],lon,lat,n,NAN,verbose);
  delete[]lat;
  delete[]lon;}
  
/*-----------------------------------------------------------------------------
  predict, correct and write output */
  if(verbose)__OUT_BASE_LINE__("Calling init_argument():\n");
  const double t0=epochs[0].time;
  init_argument(&astro_angles,t0,verbose);
  
  if(verbose)__OUT_BASE_LINE__("Writing %s...\n",outP);
  free(outP);
  
  fprintf(outF,"#file produced with : "+cmd+"\n");
  fprintf(outF,"#time latitude longitude sounding prediction depth");
  for(j=0;j<WaveList.n;j++){
    const tidal_wave *wave=&WaveList.waves[j];
    fprintf(outF," %s_amplitude %s_phase",wave->name,wave->name);
    }
  fprintf(outF,"\n");
  
  for(i=0;i<n;i++){
    epoch_t *e=&epochs[i];
    char *sdate;
    double prediction;
    double t=e->time-t0;
    
    sdate=poctime_sdate_cnes(e->time,'s');
    fprintf(outF,"%s",sdate);
    delete[]sdate;
    
    fprintf(outF," %-9.9g %-9.9g %-9.9g",e->lat,e->lon,e->inD);
    
    harmonic_prediction(astro_angles,&prediction,t,1,WaveList,&constants[i],nodal_corrections);
    e->outD=e->inD-prediction;
    
    fprintf(outF," %-9.9g %-9.9g",prediction,e->outD);
    
    for(j=0;j<WaveList.n;j++){
      fprintf(outF," %-9.9g %-9.9g",constants[i].a[j],constants[i].G[j]);
      }
    fprintf(outF,"\n");
    }
  
  fclose(outF);
  
  if(verbose)__OUT_BASE_LINE__("done.\n");
  
  return 0;
}
