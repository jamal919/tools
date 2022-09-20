
/*******************************************************************************
T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative
  
E-mail: florent.lyard@legos.obs-mip.fr
*******************************************************************************/
/**
\file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada
\author  Clement MAYET      LEGOS, Toulouse, France (PhD)
\author  Yves Soufflet      LEGOS, Toulouse, France

VERSION :
\date  11/04/2011 : Clement MAYET: add SONEL_HR input format, add a print_help function and some Doxygen comments.
\date  24/04/2011 : Clement MAYET: Bugfix in arguments reading
\date  23/06/2011 : Yves Soufflet: Add input format for profilers

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief This program is a converter between different harmonic analyses formats

\todo 12/04/2011 : Clement MAYET: re-organize the functions (group the read functions, the write functions, the others... maybe in different files ) ?

*******************************************************************************/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <map>

#include "tools-structures.h"

#include "fe.h"
#include "archive.h"
#include "netcdf-proto.h"
#include "poc-time.h"
#include "mgr.h"
#include "functions.h"
#include "spectrum.h"

using namespace std;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double adjust_phase(double phase)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double G;
  
  G=fmod(phase,360.0);
  
  if (G<0.0) G+=360.;
  
  return(G);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int timezone_correction(double shift, vector<mgr_t> & mgr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
 
  shift in minutes

------------------------------------------------------------------------------*/
{
  int n;

  bool Z0=true;
  spectrum_t WaveList=spectrum_init_ref("ESTUARINE",Z0);
  
  printf("%s : apply %lf mn time zone correction\n",__func__,shift);
  
  for(n=0;n<mgr.size();n++) {
    for(int w=0;w<mgr[n].nwave;w++) {
      string name=mgr[n].data[w].constituent.name;
      tidal_wave wave=WaveList.wave(name.c_str());
      wave.init();
/* *----------------------------------------------------------------------------
      compute phase lag difference*/
      double dG=shift*wave.omega/60.;
      mgr[n].data[w].phi+=dG;
      mgr[n].data[w].phi=adjust_phase(mgr[n].data[w].phi);
      }
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(char *prog_name)

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
    "  %s FILE [OPTIONS] \n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  This program is a converter between different harmonic analyses formats\n"
    "\n"
    "OPTIONS"
    "  -s  time shift (in minutes) to apply to original data phase-lag\n"
    "  -o  specify the output filename.\n"
    "      Default is the input file name with the extension changed to match that of the output format.\n"
    "  --header  contains information about the mgr station that can be printed in the output file.\n"
    "            Keywords are lon, lat, code, name and depth.\n"
    +isFormat_help("  ","",". Default: " MGR_FORMAT_NAME_LEGOS_ASCII)+
    "\n"
    "FORMATS\n"
    );
  print_mgrh_formats();
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  convert tide gauges harmonic data file from one format to another

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int n,status;
  char *input=0;
  char *output=0, *header=0, *wave=0, *poly=0;
  string input_format="",output_format=MGR_FORMAT_NAME_LEGOS_ASCII;
  string rootname="", keyword;
  double shift=0.0;
  bool exclude=true;

  fct_echo( argc, argv);
  
  status= init_mgrh_formats();

  n=1;
  while (n < argc) {
    const char *argvn=argv[n];
    
    if(isFormat(argv,&n,&input_format,&output_format)){
      continue;
      }
    keyword=argv[n];
    if(keyword=="--include"){
      poly= strdup(argv[n+1]);
      exclude=false;
      n++;
      n++;
      continue;
      }
    
    if(keyword=="--exclude"){
      poly= strdup(argv[n+1]);
      exclude=true;
      n++;
      n++;
      continue;
      }
    
    switch (argvn[0]) {
      case '-':
        switch (argvn[1]) {
        case 'p' :
          poly= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          sscanf(argv[n+1],"%lf",&shift);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;
          
        case 'I' :
          input_format=argv[n+1];
          n++;
          n++;
          break;
          
        case 'O' :
          output_format=argv[n+1];
          n++;
          n++;
          break;
          
        case 'w' :
          wave= strdup(argv[n+1]);
          n++;
          n++;
          break;
          
        case 'h' :
          print_help(argv[0]);
          wexit(0);

        case '-' :
          if(strstr(argvn,"--header=")!=0){
            header=strdup(argvn);
            n++;
            break;
            }
          break;
          if(strcmp(argvn,"--help")==0){
            print_help(argv[0]);
            wexit(0);
            }
          break;

        default:
          printf("unknown option %s\n",argvn);
          print_help(argv[0]);
          wexit(-1);
        }
        break;

      default:
        if(input==NULL) {
          input= strdup(argvn);
          n++;
          }
        else {
          printf("redundant option %s\n",argvn);
          print_help(argv[0]);
          wexit(-1);
          }
        break;
      }
    }


  printf("#################################################################\n");
  printf("load tide gauge file : %s\n",input);
  vector<mgr_t> mgr;
  int nmgr;

  int informat=MgrHarmonic_formats[input_format];
  nmgr=mgr_load(input, mgr, informat);
  printf("number of station found : %d\n",nmgr);
  
  if(poly!=0) {
    if(exclude) {
      status=mgr_exclude(mgr, poly, PLG_POINT_INTERIOR);
      }
    else {   
      status=mgr_exclude(mgr, poly, PLG_POINT_EXTERIOR);
      }
    }

  if(output==0) {
    rootname=input;
    size_t l=rootname.length();
    
    if(strrncmp(rootname,".mgr")==0)
      rootname.erase(l-4);
    else if(strrncmp(rootname,".nc")==0)
      rootname.erase(l-3);
    }

  if(shift!=0.0)
    status=timezone_correction(shift, mgr);
  
  if(output==0) {
    if(output_format=="ASCII" or output_format=="LEGOS-ASCII"){
      rootname+=".mgr";
      }
    else if(output_format=="NETCDF" or output_format=="LEGOS-NETCDF"){
      rootname+=".nc";
      }
    else if(output_format=="LEGOS-OBS"){
      rootname+=".obs";
      }
    output=strdup(rootname.c_str());
    }

  printf("#################################################################\n");
  printf("save tide gauge  file : %s\n",output);
  status=mgr_save(output, mgr, output_format, wave);

  STDOUT_BASE_LINE(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  wexit(0);
}
