
/*******************************************************************************
T-UGO tools, 2006-2012

  Unstructured Ocean Grid initiative
  
E-mail: florent.lyard@legos.obs-mip.fr
*******************************************************************************/
/**
\file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada
\author  Clement MAYET      LEGOS, Toulouse, France (PhD)
\author  Yves Soufflet      LEGOS, Toulouse, France

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Validate time series ?

*******************************************************************************/


#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <map>

#include "tools-define.h"
#include "tools-structures.h"

#include "fe.h"
#include "archive.h"
#include "poc-netcdf.hpp"
#include "poc-time.h"
#include "mgr-converter.h"
#include "mgr.h"
#include "spectrum.h"
#include "filter.h"
#include "functions.h"
#include "statistic.h"

using namespace std;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int parse_InputOptions(const string options, char* & filename, char* & format, double & mask, char & units, int & column)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<string> tokens,keys,values;
  int k;
  string delimiter=" ";
  
  filename=0;
  
  delimiter=" ";
  tokens=string_split(options, delimiter);
  
  delimiter="=";
  for(k=0;k<tokens.size();k++) {
    vector<string> tmp=string_split(tokens[k], delimiter);
    keys.push_back(tmp[0]);
    values.push_back(tmp[1]);
    }
    
  for(k=0;k<keys.size();k++) {
    if(keys[k]=="file") {
      filename=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="format") {
      format=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="mask") {
      sscanf(values[k].c_str(),"%lf", &mask);
      continue;
      }
    if(keys[k]=="units") {
      sscanf(values[k].c_str(),"%c", &units);
      continue;
      }
    if(keys[k]=="column") {
      sscanf(values[k].c_str(),"%d", &column);
      continue;
      }
    }
  return(0);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int validate(tseries_t serie, tseries_t & data, tseries_t & simulation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  statistic_t s;
  double limit=0.1;
  
  s=OutliersFrequency(serie.x[2], serie.mask, serie.n, limit);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int validate(tseries_t & data, tseries_t & simulation, double offset, mooring_t mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  date_t start,final;
  statistic_t s;
  
  std::string line;
    
  tseries_t serie;

  serie.nparam=3;
  serie.n=data.n;
  serie.mask=data.mask;
  serie.t=data.t;
 
  serie.x=new double*[serie.nparam];
  
  for(int k=0;k<serie.nparam;k++) {
    serie.x[k]=new double[serie.n];
    for (int i=0; i< serie.n; i++) {
      serie.x[k][i]=serie.mask;
      }
    }

  for (int i=0; i< serie.n; i++) {
    if(data.x[0][i]!=data.mask) serie.x[0][i]=data.x[0][i]-offset;
    }
    
  int nRS __attribute__((unused)) =serie.size(1,serie.mask);
  int nLF __attribute__((unused)) =serie.size(2,serie.mask);
    
/*------------------------------------------------------------------------------
  identify alleged time sampling */
  double sampling=get_timesampling(serie.t, serie.mask, serie.n);
  sampling=floor(sampling*24.*3600.+0.5)/24/3600.;
  printf("sampling appears to be %lf mn\n",sampling*24.*60.);
  
/*------------------------------------------------------------------------------
  resample simulation time serie*/
  printf("number of valid simulation values before reconstruction: %d (over %d)\n",simulation.n-simulation.size(0,simulation.mask),simulation.n);
  tseries_t tmp;
  
/*------------------------------------------------------------------------------
  do not re-interpolate if 2 successive time larger than max_gap (assign mask)*/
  double max_gap=1./24/1.+0.01;
  tmp=simulation.resample(serie.t, serie.n, max_gap);
  printf("number of valid simulation values after reconstruction: %d (over %d)\n",tmp.n-tmp.size(0,tmp.mask),tmp.n);
  
  for (int i=0; i< serie.n; i++) {
    serie.x[1][i]=tmp.x[0][i];
    serie.x[2][i]=serie.mask;
    }
    
  double limit=0.1;
  
  for (int i=0; i< serie.n; i++) {
    if( (serie.x[0][i]!=serie.mask) && (serie.x[1][i]!=serie.mask) ) {
      serie.x[2][i]=serie.x[1][i]-serie.x[0][i];
      }
    }
  int nValid __attribute__((unused)) =serie.size(2,serie.mask);
  
  s=get_statistics(serie.x[2], serie.mask, serie.n);
  
  s=OutliersFrequency(serie.x[2], serie.mask, serie.n, limit);
  
  return(1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int proccess(const char *DataInput, const char *ModelInput, date_t start, date_t final, double offset)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  tseries_t *data, *simulation;
  char *filename=0, *header_format=0, *format=0, units;
  int nb_col_GNU=10, nseries, column;
  mooring_t mooring;
  string options;
  double mask;
  
  status= mgr_init_formats();
  
  string DataString="file=/home/data/sealevel/seine/11-CAUDEBEC_Caudebec_20061231_20110101.gnu format=GNU units=m";
  string ModelString="file=/home/models/seine/fleuve_trace_alti/simulation-tides+rivers+DAC/Cd-1.5-NEA-smag0.05/sample.8 format=GNU units=m";
  
  DataString="file=/home/data/sealevel/seine/14-ST_LEONARD_Saint_Leonard_20061231_20110101.gnu format=GNU units=m";
  ModelString="file=/home/models/seine/fleuve_trace_alti/simulation-tides+rivers+DAC/Cd-1.5-NEA-smag0.05/sample.5 format=GNU units=m";

  DataString="file=/home/data/sealevel/seine/16-TANCARVILLE_Tancarville_20061231_20110101.gnu format=GNU units=m";
  ModelString="file=/home/models/seine/fleuve_trace_alti/simulation-tides+rivers+DAC/Cd-1.5-NEA-smag0.05+Hcorr-v1/sample.4 format=GNU units=m";

/*------------------------------------------------------------------------------
  load observation time series */

  options=DataString;
  status=parse_InputOptions(options, filename, format, mask, units, column);
  
  nseries = mgr_load_timeserie(filename, &data, units, format, nb_col_GNU, header_format, &mooring);
  
/*------------------------------------------------------------------------------
  load simulation time series */

//  options="file=/home/models/seine/fleuve_trace_alti/simulation-tides+rivers+DAC/Cd-1.5-NEA-smag0.05+Hcorr/sample.8 format=GNU units=m";
  options=ModelString;
  status=parse_InputOptions(options, filename, format, mask, units, column);
  
  nseries = mgr_load_timeserie(filename, &simulation, units, format, nb_col_GNU, header_format, &mooring);
  
  tseries_t data_reduced = data[0].reduce(start, final);
  
  status=validate(data_reduced, simulation[0], offset, mooring);
  
  return(1);
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
         "NAME AND VERSION:\n"
         "  %s version " VERSION " " REVISION "\n", prog_name);

  printf("USAGE:\n"
         "  %s <file> [options]\n"
         "\n", prog_name);

   printf("DESCRIPTION:\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
          "  Validate time series ?\n"
          "\n"
          "NONE OF THESE OPTIONS HAVE ANY EFFECT :\n"
          +isFormat_help("   ",". Default: GNU",". Default: GNU")+
          "   -n <nb>                  : nb of columns for input GNU format (default=1)\n"
          "   -s <dd/mm/yyyy>          : start date\n"
          "   -f <dd/mm/yyyy>          : final date\n"
          "   --step=<n_years>y        : split time series per nb years or nb month\n"
          "   --mooring_info=<info>    : input info on mooring, eg: \"lon=<lon> lat=<lat> code=<code> name=<name depth=<depth>\"\n"
          "   --header_format=<format> : format to decode mooring info on file\n"
          "   -o <name>                : to save output in files with prefix <name>\n"
          "   --gnu_nice               : nice gnu (allways set to True)\n"
          );

   mgr_print_formats();

  TRAP_ERR_EXIT(-1,"exiting\n"); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int ncol=1;
  date_t step(0,0,0,0.);
  int n,status;
  char *input=NULL,*s;
  char *output=NULL,*mooring_info=NULL,*header_format=NULL;
  string input_format;
  string output_format="GNU";
  string rootname;
  bool gnu_nice=1;
  date_t start(0,0,0,0.),final(0,0,0,0.);
  int reconstruct=0;
  double time_correction=0;

  string DataInput, ModelInput;
  
  fct_echo( argc, argv);

  if (argc==1)
    print_help(argv[0]);

  n=1;
  while (n < argc) {
    const char *keyword=argv[n];
    
    if(isFormat(argv,&n,&input_format,&output_format)){
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 's' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          start = poctime_scan_start_date(s);
          if (start.year == 0) {
            TRAP_ERR_EXIT(-1,"mgr_converter: poctime_scan_start_date failed\n");
            }
          break;

       case 'f' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          final = poctime_scan_final_date(s);
          if (final.year == 0) {
            TRAP_ERR_EXIT(-1,"mgr_converter: poctime_scan_final_date failed\n");
            }
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'n' :
          sscanf(argv[n+1],"%d",&ncol);
          n++;
          n++;
          break;
          
        case 'h' :
          print_help(argv[0]);

        case '-' :
          if(strncmp(keyword,"--nice=")==7 || strcmp(keyword,"--nice=yes")==0 || strcmp(keyword,"--gnu_nice")==0) {
            gnu_nice=1;
            printf ("gnu_nice is now always on\n");
            n++;
            break;
            }
          if(strcmp(keyword,"--reconstruct")==0 ) {
            reconstruct=1;
            printf ("series reconstruction on\n");
            n++;
            break;
            }
          if(strcmp(keyword,"--time-correction")==0 ) {
            sscanf(argv[n+1],"%lf",&time_correction);
            printf ("time correction (additive)=%lf seconds\n",time_correction);
            n++;
            n++;
            break;
            }
          if(strncmp(argv[n],"--header=",9)==0){
            mooring_info=&(argv[n][9]);
            n++;
            break;
            }
          if(strncmp(argv[n],"--mooring_info=",15)==0){
            mooring_info=&(argv[n][15]);
            n++;
            break;
            }
          if(strncmp(argv[n],"--header_format=",16)==0){
            header_format = &(argv[n][16]);
            n++;
            break;
            }
          if(strncmp(argv[n],"--step=",7)==0){
            int nb_step;
            char unit;
            sscanf(&(argv[n][7]),"%d%c",&nb_step,&unit);
            if (unit=='y') step.year = nb_step;
            else if (unit=='m') step.month = nb_step;
            else if (unit=='d') step.day = nb_step;
            else {
              print_help(argv[0]);
              STDOUT_BASE_LINE("bad syntax %s for option --step\n",argv[n]);
              exit(-1);
              }
            n++;
            break;
            }
          printf("unknown option %s\n",argv[n]);
          print_help(argv[0]);
          exit(-1);
          break;

        default:
          printf("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
        if(input==NULL) {
          input= strdup(argv[n]);
          n++;
          }
        else {
          printf("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
          }
        break;
      }
    
    }

  double offset=5.0;

  status= proccess(DataInput.c_str(), ModelInput.c_str(), start, final, offset);
  
  STDOUT_BASE_LINE(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(0);
}
