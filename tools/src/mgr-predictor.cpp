
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\date 2011-10-25 Damien ALLAIN : creation

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Predicts tides from given mgr (tidal constituents) file

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
#include "mgr.h"
#include "mgr-converter.h"
//#include "tides.def"
#include "poc-time.h"
#include "poc-netcdf-data.hpp"
#include "fe.h"
#include "archive.h"

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
    "  %s -m input_mgr_file -s start_date -f end_date -i time_increment UNITS -w wave1 [wave2 ... ] [OPTIONS] \n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Predicts tides from given mgr (tidal constituents) file .\n"
    "  If start and end dates are provided, predicts tides (and spring/neap times with --extrema).\n"
    "  Input file is ascii_mgr or netcdf format, and contains tidal waves amplitude/phase lag for one or multiple stations.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  show this help\n"
    "  -m  followed by the name of input mgr file (tidal waves amplitude and phase from which prediction is computed)\n"
    "  -s  followed by the start date. See DATE FORMATS below\n"
    "  -f  followed by the end date. See DATE FORMATS below\n"
    "  -i  followed by the date increment concatenated with the unit (WITHOUT SPACE! : -i 600s): s (default), m, h or d (day). Default increment: 3600s=60m=1h.\n"
    "  -w  followed by the list of waves to predict for\n"
    "  --nodal=no do not do nodal corrections\n"
    "  --extrema NOT WORKING  activate extrema and spring neap tide output\n"
    "  -o  NOT USED ... followed by the path of the ascii output. Default: predictions.dat\n"
    );
  print_poctime_scan_date_help(0);
    /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int extract(const char *output, const char *name_template, const char *varname, const char *tname, double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  poc_global_t global;
  poc_data_t<double> var,time;
  char *input;
  mesh_t mesh;
  int discretisation=LGP1;
  discretisation_t descriptor;
  char *comments[2];
  size_t i;
  
  date_t date=poctime_getdatecnes(t, 'd');
  
  status=decode_name(date, 0, name_template, &input);
  status=fe_readgeometry(input, &mesh);
 
  discretisation_init(&mesh,discretisation,0);
  descriptor=get_descriptor(mesh,discretisation);
  
  status=poc_inq(input,&global);

//   if(status!=0) break;

  status=var.init(global,varname);
      
  status=time.init(global,tname);
  status=time.read_data(input,0,0,1);
  
  double delta=fabs(t-time.data[0]/d2s), nearest;
  nearest=time.data[0]/d2s;
  int frame=0;
  for (i=0;i<time.nframes;i++) {
    status=time.read_data(input,i,0,1);
    time.data[0]/=d2s;
    if(fabs(t-time.data[0])<delta) {
      delta=fabs(t-time.data[0]);
      nearest=time.data[0];
      frame=i;
      }
    }

  printf("targeted/nearest cnes time %s %s\n", poctime_sdate_cnes(t,'d',' '), poctime_sdate_cnes(nearest,'d',' '));
  
  status=var.read_data(input,frame,0,1);
 
 
  comments[0]=new char[256];
  comments[1]=new char[256];
  
  sprintf(comments[0],"%c POC/Noveltis, extracted from %s",169,input);
  sprintf(comments[1],"variable %s",varname);
  status=quoddy_saver1(output, mesh.nvtxs, var.data,comments);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int local_extrema(tseries_t & serie, vector<int> & imin, vector<int> & imax)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int p=0, status=0;
  size_t n;
  double *s=new double[serie.n], *laplacian=new double[serie.n];
  int state=0;
//   vector<double> time;
  bool debug=false;

/*------------------------------------------------------------------------------
  first step : find all local min and max (extrema) */
  for(n=1;n<serie.n-1;n++) {
    s[n]=(serie.x[p][n+1]-serie.x[p][n])*(serie.x[p][n]-serie.x[p][n-1]);
    laplacian[n]=serie.x[p][n+1]+serie.x[p][n-1]-2.*serie.x[p][n];
    }
    
  for(n=1;n<serie.n-1;n++) {
    if(s[n]<0) {
      if(laplacian[n]<0) {
        if(debug) printf("maximum at %lf, %lf\n",serie.t[n], serie.x[p][n]);
        imax.push_back(n);
        if(state==1) {
          printf("anomaly at %lf, %lf\n",serie.t[n], serie.x[p][n]);
          }
        state=1;
        }
      if(laplacian[n]>0) {
        if(debug) printf("minimum at %lf, %lf\n",serie.t[n], serie.x[p][n]);
//         time.push_back(serie.t[n]);
        imin.push_back(n);
        if(state==-1) {
          printf("anomaly at %lf, %lf\n",serie.t[n], serie.x[p][n]);
          }
//         if(serie.x[p][n]>0.) {
//           printf("minimum at %lf, %lf\n",serie.t[n], serie.x[p][n]);
//           }
        state=-1;
        }
      }
    }
  delete[] laplacian;
  delete[] s;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int global_extrema_index(tseries_t & serie, vector<int> & pos, int minmax, double lap=-1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int p=0, status;
  size_t k, m, n;
  vector<double> delta;
  bool debug=false;

  for(k=0;k<pos.size()-1;k++) {
    m=pos[k];
    n=pos[k+1];
    delta.push_back(serie.t[n]-serie.t[m]);
    }
  status=median_shell (delta);
  
  double median=delta[delta.size()/2];
  printf("median = %lf\n",median);
  
  updatemax(&median,lap);
  
/*------------------------------------------------------------------------------
  third step : select main extrema */
  for(k=0;k<pos.size()-1;k++) {
    m=pos[k];
    n=pos[k+1];
    double d=serie.t[n]-serie.t[m];
    if(d<0.60*median) {
      if(debug) printf("local duplicated extrema at %lf %lf, %lf %lf\n",serie.t[m],serie.t[n], serie.x[p][m], serie.x[p][n]);
      switch (minmax) {
        case -1:
          if(serie.x[p][m] > serie.x[p][n]) {
            pos.erase(pos.begin()+k);
            }
          else {
            pos.erase(pos.begin()+k+1);
            }
          break;
        case 1:
          if(serie.x[p][m] < serie.x[p][n]) {
            pos.erase(pos.begin()+k);
            }
          else {
            pos.erase(pos.begin()+k+1);
            }
          break;
        }
      k--;
      }
    }
    
  for(k=0;k<pos.size();k++) {
    n=pos[k];
    if(debug) printf("extrema at %lf, %lf\n",serie.t[n], serie.x[p][n]);
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  tseries_t global_extrema(tseries_t & serie, vector<int> & pos, int minmax, double lap=-1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int k,n,p=0;
  status=global_extrema_index(serie, pos, minmax, lap);
  
  tseries_t out(pos.size(),1);
  out.mask=1.e+10;

  for(k=0;k<pos.size();k++) {
    n=pos[k];
    out.t[k]=serie.t[n];
    out.x[0][k]=serie.x[p][n];
    }
 
  return(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int extrema(tseries_t & serie)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int p=0, status;
  size_t k, n;
  vector<int> imin, imax;
  bool debug=false;
  tseries_t hightide, lowtide;
  tseries_t neap, spring;

/*------------------------------------------------------------------------------
  first step : find all local min and max (extrema) forlow/high tides */
  status=local_extrema(serie, imin, imax);

/*------------------------------------------------------------------------------
  second step : eliminate secundary extrema */
  lowtide  = global_extrema(serie, imin, -1);
  hightide = global_extrema(serie, imax, +1);
    
  for(k=0;k<imin.size();k++) {
    n=imin[k];
    if(debug) printf("minimum at %lf, %lf\n",serie.t[n], serie.x[p][n]);
    }
   
  imin.clear();
  imax.clear();
 
/*------------------------------------------------------------------------------
  third step : locate spring/neap tides */
  
/*------------------------------------------------------------------------------
  spring/neap tides from high tides*/
  status = local_extrema(hightide, imin, imax);
  neap   = global_extrema(hightide, imin, -1, 14.0);
  spring = global_extrema(hightide, imax, +1, 14.0);
 
  mooring_t mooring;
  mgr_save_timeserie("hightide",  mooring, hightide, 'm', (string) FORMAT_NAME_GNU, false);
  mgr_save_timeserie("spring-ht", mooring, spring,   'm', (string) FORMAT_NAME_GNU, false);
  mgr_save_timeserie("neap-ht",   mooring, neap,     'm', (string) FORMAT_NAME_GNU, false);
  
  spring.destroy();
  neap.destroy();
 
  imin.clear();
  imax.clear();
 
/*------------------------------------------------------------------------------
  spring/neap tides from low tides*/
  status = local_extrema(lowtide, imin, imax);
  spring = global_extrema(lowtide, imin, -1, 14.0);
  neap   = global_extrema(lowtide, imax, +1, 14.0);
 
  mgr_save_timeserie("lowtide",   mooring, lowtide, 'm', (string) FORMAT_NAME_GNU, false);
  mgr_save_timeserie("spring-lt", mooring, spring,  'm', (string) FORMAT_NAME_GNU, false);
  mgr_save_timeserie("neap-lt",   mooring, neap,    'm', (string) FORMAT_NAME_GNU, false);
  
  spring.destroy();
  neap.destroy();
   
  imin.clear();
  imax.clear();
   
/*------------------------------------------------------------------------------
  spring/neap tides from range */
  tseries_t range(min(lowtide.n, hightide.n), 1);
  for(n=0;n<range.n;n++) {
    range.x[0][n]=hightide.x[0][n]-lowtide.x[0][n];
    range.t[n]=0.5*(hightide.t[n]+lowtide.t[n]);
    }
    
  status = local_extrema(range, imin, imax);
  neap   = global_extrema(range, imin, -1, 14.0);
  spring = global_extrema(range, imax, +1, 14.0);
 
  mgr_save_timeserie("range",     mooring, range,  'm', (string) FORMAT_NAME_GNU, false);
  mgr_save_timeserie("spring-tr", mooring, spring, 'm', (string) FORMAT_NAME_GNU, false);
  mgr_save_timeserie("neap-tr",   mooring, neap,   'm', (string) FORMAT_NAME_GNU, false);
  
    
  for(k=0;k<imin.size();k++) {
    printf("neap tide   : %10.5lf julian day/cnes time %s, range = %9.2lf m\n",neap.t[k], poctime_sdate_cnes(neap.t[k],'d',' '), neap.x[0][k]);
    }
 
  for(k=0;k<imax.size();k++) {
    printf("spring tide : %10.5lf julian day/cnes time %s, range = %9.2lf m\n",spring.t[k], poctime_sdate_cnes(spring.t[k],'d',' '),spring.x[0][k]);
    }
    
  size_t *springlist = sort(spring.x[0],spring.n);
  size_t *neaplist   = sort(neap.x[0],  neap.n);
  
  k=springlist[0];
  printf("lowest spring tide  : %10.5lf julian day/cnes time %s, range = %9.2lf m\n",spring.t[k], poctime_sdate_cnes(spring.t[k],'d',' '), spring.x[0][k]);
  n=imax[k];
  printf("high tide : %10.5lf julian day/cnes time %s, range = %9.2lf m\n",hightide.t[n],   poctime_sdate_cnes(hightide.t[n],'d',' '), hightide.x[0][n]);
  printf("low  tide : %10.5lf julian day/cnes time %s, range = %9.2lf m\n\n",lowtide.t[n],  poctime_sdate_cnes(lowtide.t[n], 'd',' '), lowtide.x[0][n]);
  
  k=springlist[spring.n-1];
  printf("highest spring tide : %10.5lf julian day/cnes time %s, range = %9.2lf m\n",spring.t[k], poctime_sdate_cnes(spring.t[k],'d',' '), spring.x[0][k]);
  n=imax[k];
  printf("high tide : %10.5lf julian day/cnes time %s, range = %9.2lf m\n",hightide.t[n],   poctime_sdate_cnes(hightide.t[n],'d',' '), hightide.x[0][n]);
  printf("low  tide : %10.5lf julian day/cnes time %s, range = %9.2lf m\n\n",lowtide.t[n],  poctime_sdate_cnes(lowtide.t[n], 'd',' '), lowtide.x[0][n]);
  
  status=extract("spring-hightide.s2r", "/home/models/seine/fleuve_trace_alti-v2/simulation-tides+rivers+DAC/out_debit_reel_alti_GIP/tugo.state-YYYY.MM.nc", "elevation", "time", hightide.t[n]);
  status=extract("spring-lowtide.s2r", "/home/models/seine/fleuve_trace_alti-v2/simulation-tides+rivers+DAC/out_debit_reel_alti_GIP/tugo.state-YYYY.MM.nc", "elevation", "time", lowtide.t[n]);
  
  k=neaplist[0];
  printf("lowest neap tide    : %10.5lf julian day/cnes time %s, range = %9.2lf m\n",neap.t[k],   poctime_sdate_cnes(neap.t[k],  'd',' '), neap.x[0][k]);
  n=imin[k];
  printf("high tide : %10.5lf julian day/cnes time %s, range = %9.2lf m\n",hightide.t[n],   poctime_sdate_cnes(hightide.t[n],'d',' '), hightide.x[0][n]);
  printf("low  tide : %10.5lf julian day/cnes time %s, range = %9.2lf m\n\n",lowtide.t[n],  poctime_sdate_cnes(lowtide.t[n], 'd',' '), lowtide.x[0][n]);
  
  status=extract("neap-hightide.s2r", "/home/models/seine/fleuve_trace_alti-v2/simulation-tides+rivers+DAC/out_debit_reel_alti_GIP/tugo.state-YYYY.MM.nc", "elevation", "time", hightide.t[n]);
  status=extract("neap-lowtide.s2r", "/home/models/seine/fleuve_trace_alti-v2/simulation-tides+rivers+DAC/out_debit_reel_alti_GIP/tugo.state-YYYY.MM.nc", "elevation", "time", lowtide.t[n]);
  
  k=neaplist[neap.n-1];
  printf("highest neap tide   : %10.5lf julian day/cnes time %s, range = %9.2lf m\n",neap.t[k],   poctime_sdate_cnes(neap.t[k],  'd',' '), neap.x[0][k]);
  n=imin[k];
  printf("high tide : %10.5lf julian day/cnes time %s, range = %9.2lf m\n",hightide.t[n],   poctime_sdate_cnes(hightide.t[n],'d',' '), hightide.x[0][n]);
  printf("low  tide : %10.5lf julian day/cnes time %s, range = %9.2lf m\n\n",lowtide.t[n],  poctime_sdate_cnes(lowtide.t[n], 'd',' '), lowtide.x[0][n]);

  spring.destroy();
  neap.destroy();
 
  imin.clear();
  imax.clear();
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void SpringNeap_cycle(astro_angles_t astro_angles, spectrum_t & s, hconstant_t & constants, int nodal, date_t start, date_t final)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  tidal_wave wave,wM2,wS2;
  double V0,omega;
  double f,Vu;
  double VM2,VS2;
  double time=0;
  
/**--------------------------------------------------------------------------
  compute spring/neap tide calendar */

  int M2id=s.wave_index("M2");
  int S2id=s.wave_index("S2");
  bool abortThis=false;
  if(M2id<0){
    STDERR_BASE_LINE("Can not find M2 in the list of waves:\n");
    abortThis=true;
    }
  if(S2id<0){
    STDERR_BASE_LINE("Can not find S2 in the list of waves:\n");
    abortThis=true;
    }
  if(abortThis){
    STDERR_BASE_LINE("aborting %s\n",__func__);
    return;
    }
  
  wM2=s.wave("M2");
//  wM2.init();
  wS2=s.wave("S2");
  
  double d=wS2.omega-wM2.omega;
//  double r=2*(wS2.omega-wM2.omega)/(wS2.omega+wM2.omega);
  
  double duration=360./24./d;
  printf("spring/neap tides cycle: %lf days\n",duration);
  
  wave=wM2;
  omega=wave.omega*dph2rpd;
  V0=greenwhich_argument(astro_angles,wave);
  if(nodal==1) {
    Vu=nodal_phase(astro_angles,wave);
    f=nodal_factor(astro_angles,wave.formula);
    }
  else {
    Vu=0.;
    f=1.;
    }
  VM2=omega*time/3600.+V0+Vu;
  
  wave=wS2;
  omega=wave.omega*dph2rpd;
  V0=greenwhich_argument(astro_angles,wave);
  if(nodal==1) {
    Vu=nodal_phase(astro_angles,wave);
    f=nodal_factor(astro_angles,wave.formula);
    }
  else {
    Vu=0.;
    f=1.;
    }
  VS2=omega*time/3600.+V0+Vu;
  
/**--------------------------------------------------------------------------
  M2 & S2 phases at start (time=0) */
  VM2-=constants.G[M2id]*d2r;
  VS2-=constants.G[S2id]*d2r;
  
  double dV=fmod(VM2-VS2,2*M_PI);
  dV-=2*M_PI;
  
/**--------------------------------------------------------------------------
  convert seconds into julian days */
  double cnes_start=cnes_time(start,'s');
  cnes_start/=d2s;
  double cnes_final=cnes_time(final,'s');
  cnes_final/=d2s;
  
  double t=cnes_start;
  while(t<cnes_final) {
/**--------------------------------------------------------------------------
    in-phase M2 S2 */
    double spring=dV/(d*dph2rpd);
    printf("spring tide : %10.5lf %10.5lf julian day/cnes time %s\n",spring,spring+cnes_start, poctime_sdate_cnes(spring+cnes_start,'d',' '));
    dV+=M_PI;
/**--------------------------------------------------------------------------
    in-opposition M2 S2 */
    double neap=dV/(d*dph2rpd);
    printf("neap tide   : %10.5lf %10.5lf julian day/cnes time %s\n",neap,neap+cnes_start, poctime_sdate_cnes(neap+cnes_start,'d',' '));
    dV+=M_PI;
    t=neap+cnes_start;
    }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// function that will be ran when the executable is started
/** See the print_help <a href=#func-members>function</a> for its use.
*/
/*----------------------------------------------------------------------------*/
{
  int wI,n,status,k;    //wave index, argument index, status, variable index

  date_t start=NADate,final=NADate;
  double increment=NAN;
  char *unit;                          //unit of time increment
  char *waveNames[100];                //list of wave names
  int nodal_corrections=1;
  astro_angles_t astro_angles;
  char *output=NULL;                   //path of the ascii output
  const char *keyword,*s;                    //option and the argument that follows

  vector<string> stations;
  char *mgrfile=NULL;
  bool doSpringNeap=false;
  char *time_template=NULL;

  string cmd;
  cmd=fct_echo( argc, argv);

  wI=0;
  waveNames[wI]=NULL;

  n=1;
  while (n < argc) {
    keyword=argv[n];
    if(strcmp(keyword,"--nodal=no")==0) {
      nodal_corrections=0;
      n++;
      }
    else if(strcmp(keyword,"--extrema")==0) {
      doSpringNeap=true;
      n++;
      }
    else {
      switch (keyword[0]) {
      case '-':
        if( strcmp(keyword,"-h")==0 or
            strcmp(keyword,"--help")==0 ){
          print_help(argv[0]);
          return 0;
          }
        if(strcmp(keyword,"--time")==0) {
          time_template=strdup(argv[n+1]);
          n++;
          n++;
          continue;
          }
        switch (keyword[1]) {

        case 'm' :
          mgrfile=strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          s=argv[n+1];
          n++;
          n++;
          status=poctime_scan_date(s,&start,0);
          break;

        case 'f' :
          s=argv[n+1];
          n++;
          n++;
          status=poctime_scan_date(s,0,&final);
          break;

        case 'i' :
          s=argv[n+1];
          n++;
          n++;
          increment=strtod(s,&unit);
          switch(tolower(*unit)){
                case 'd':
                  increment*=24.;
                case 'h':
                  increment*=60.;
                case 'm':
                  increment*=60.;
                case 0:
                case 's':
                  printf("Time increment of %g s.\n",increment);
                  break;
                default:
                  fprintf(stderr,"unit %s not recognised.\n");
                  print_help(argv[0]);
                  exit(1);
                }
          break;

        case 'w' :
          n++;
          wI=pos((char*)NULL,waveNames,100);
          for(;n<argc;wI++,n++) {
                waveNames[wI]=strdup(argv[n]);
                }
          waveNames[wI]=NULL;
          break;

        case 'o' :
          output=strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          printf("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
          break;
        }
      break;
    
      default:
        stations.push_back(argv[n]);
        n++;
        break;

        }
      }
    
    }

  if(isad(start)!=isad(final)){
    fprintf(stderr,"*** Provide start AND end dates for predictions ***\n");
    print_help(argv[0]);
    exit(-1);
    }

/*---------------------------------------------------------------------*//**<h1>
  get amplitude and phase () </h1>*/
  vector<hconstant_t> vconstants;
  vector<spectrum_t> vspectrum;
  vector<mooring_t> vmoorings;
  status=mgr_load((const char *) mgrfile, vmoorings, vspectrum, vconstants, MGR_FORMAT_CODE_LEGOS_ASCII);
  
/*---------------------------------------------------------------------*//**<h1>
  if there are start and end dates, predict with harmonic_prediction()</h1>*/
  if(isnad(start) || isnad(final)){
    TRAP_ERR_EXIT(0,"*** Provide start and end dates for predictions ***\n");
    }
  if(isnan(increment)){
    increment=3600.;
    printf("Default time increment of %g s.\n",increment);
    }
    
  init_argument(&astro_angles,start,1);
  status=mgr_init_formats();

  printf("Selected dates from %s to %s.\n",sgetdate(start),sgetdate(final));
    
  for(k=0;k<vmoorings.size();k++) {
    bool go_ahead;
    if(stations.size()!=0) {
      go_ahead=false;
      for(int l=0; l<stations.size(); l++) {
        if((string) vmoorings[k].name==stations[l]) {
          go_ahead=true;
          }
        }
      }
    else go_ahead=true;
    if(!go_ahead) continue;
    
    tseries_t serie;
    vconstants[k].s=&vspectrum[k];
    
    harmonic_predictionTS(serie, start, final, increment, vconstants[k], nodal_corrections);
    
    string s=(string) vmoorings[k].name+"-prediction";
    if(time_template==0) time_template=strdup("CNES");
    mgr_save_timeserie(s, vmoorings[k], serie, 'm', (string) FORMAT_NAME_GNU, false, time_template);
    
/*---------------------------------------------------------------------*//**<h1>
    SpringNeap_cycle()</h1>*/
    if(doSpringNeap){
      status=extrema(serie);
      SpringNeap_cycle(astro_angles, vspectrum[k], vconstants[k], nodal_corrections, start, final);
      }
    serie.destroy();
    };

    
  
  STDOUT_BASE_LINE(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(0);
}
