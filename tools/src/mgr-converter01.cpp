
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada
\author  Sara Fleury        LEGOS/CTOH, Toulouse, France

\brief sealevel time serie input/output
*/
/*----------------------------------------------------------------------------*/

#include "config.h"

#include <unistd.h> // for access
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <string>
#include <map>

using namespace std;

#include "tools-structures.h"

#include "poc-netcdf-data.hpp"
#include "poc-time.h"
#include "mgr-converter.h"
#include "statistic.h"
#include "functions.h"
#include "spectrum.h"
#include "filter.h"
#include "geo.h"

using namespace std;

toponyms_t timeseries_formats;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void mgr_print_formats()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints list of available mgr time series format.

\param prog_name name of the program to be printed by this help function : argv[0]
*/
/*----------------------------------------------------------------------------*/
{
  printf(
   "\n FILE FORMAT : \n"
   "      -SONEL_HR   :    'YYYY-MM-DD HH:MM:SS  slv_float' (ex: 1995-08-10 11:50:00  3.050 ). No header is read (for now ....)\n"
   "                        Extension file: .slv This format is from the 10mn SONEL data (website)\n"
   "      -YMD   'YYYY:MM:DD HH:MM    float_value_of_slv'   ( ex : 1995:08:10 11:50   3.050 ).\n"
   "                           No header is read (for now ....), \n"
   "      -LIST           'cnes_date    float_value_of_slv'  with cnes_date: julian date from 1950/01/01\n"
   "                           ? HEADER ?   NOT AVAILABLE YET  \n"
   "      -RMN            ? DESCRIBE RMN FORMAT ?\n"
   "      -PAB            ? DESCRIBE PAB FORMAT ?\n"
   "      -PROFILERS       Reads a BINARY file with format: 'nbin ntime YYYY MM DD HH MM SS  bin1 bin2 etc.'\n"
   "                         (All values are Float) No header is read (for now ....), \n"
   "      GNU      : ASCII with at least 2 columns: 'time  slv1'. For more columns use option -n<ncolumns>.\n"
   "      GLOSS or HAWAII\n"
   "      SONEL_HR\n"
   "      SEINE\n"
   "      RADAR\n"
   "      FOREMAN\n"
   "      MINH\n"
   "      LIST2\n"
   "      JMA\n"
   "      LEGEND\n"
   "      BODC\n"
   "      DART\n"
  );
/**
* \todo 11/04/2011 : Clement : complete the help function of mgr-converter.cpp, describe RMN and PAB formats
*/
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_remove_masked_data( tseries_t **series, int nseries)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int s, ndata, ind, ncol, c;
  double mask;
  int num_s_out = 0;
  int new_ind=0;

/*------------------------------------------------------------------------------
  loop over series */
  for (s=0; s<nseries; s++) {

    ndata = series[s]->n;
    mask  = series[s]->mask;
    ncol  = series[s]->nparam;

/*------------------------------------------------------------------------------
    loop over data */
    new_ind=0;
    for (ind=0; ind<ndata; ind++) {
/*------------------------------------------------------------------------------
      check mask for each param */
      for (c=0; c<ncol;c++)
        if (series[s]->x[c][ind] == mask)
          break;
      if (c<ncol) continue;
/*------------------------------------------------------------------------------
      record data to its new index */
      if (new_ind!=ind) {
        for (c=0; c<ncol;c++)
          series[num_s_out]->x[c][new_ind] = series[s]->x[c][ind];
        series[num_s_out]->t[new_ind] = series[s]->t[ind];
        }
      new_ind++;
      }

    // if (new_ind!=ind)
    STDOUT_BASE_LINE("remove %d/%d masked data\n", ind-new_ind,ind);

/*------------------------------------------------------------------------------
    record new number of data */
    series[num_s_out]->n = new_ind;

    /* empty serie : jump it : NO SO SIMPLE, NEED TO GET ALL PARAMS*/
    // if (new_ind<=1) {
    //  STDOUT_BASE_LINE("remove empty serie %d\n", num_s);
    //  continue;
    //}
    num_s_out++;
    }

  return num_s_out;
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_converter(const char *file_in, const char *format_in, string delimiter, string time_format, int ncolumns_GNU,
  const char *header_format, const char *mooring_info,
  date_t date_start, date_t date_final, date_t step, int reconstruct,double shift,
  const char *file_out_prefix,const char *format_out, tseries_t ***out_series,
  double resampling, bool gnu_nice)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief main function to read and convert time serie files

\param file_in       : name of the input ascii file
\param format_in     : format of the file: RMN PAB GNU GLOSS SONEL_HR RADAR HAWAII FOREMAN PROFILERS MINH LIST2 JMA LEGEND BODC DART  (default GNU)
\param ncolumns_GNU  : 0 or, nb of columns for GNU format (default=1)
\param mooring_info  : NULL or, input info on mooring "lon=<lon> lat=<lat> code=<code> name=<name depth=<depth>"
\param header_format : NULL or, format to read first line mooring info in file, eg: "# %s <lon> <lat>  %s <code>  %d <name> <depth>"
\param date_start    : {0,0,0} or, starting date
\param date_final    : {0,0,0} or, final date
\param step          : duration step for analysis
\param file_out_prefix      : NULL, file name to save data  (without extension!)
\param format_out    : NULL or if file_out_prefix, wanted format for file out (default GNU)

*/
/*----------------------------------------------------------------------------*/
{
  tseries_t *serie=NULL;
  int nseries;
  mooring_t mooring;
  date_t real_start, real_final;
  date_t start, final, next;
  int status;
  int j,k,n;
  int *nreduced;
  char output_dated[1024];
  char output[1024], rootname[1024];
  const char *cpos;
  char *pos;
  char mooring_name[128];
  static tseries_t **reduced;
//   bool gnu_nice=1;
  int nb_series = 0;

  string format_str="";
  if (format_out)
    format_str = format_out;

  if (file_in == NULL) {
    STDERR_BASE_LINE("input file name required \n");
    return 0;
    }

/*------------------------------------------------------------------------------
  init */
  mooring.lon=0.0;
  mooring.lat=0.0;
  
  if(mooring_info) mooring=mgr_decode_mooring_info(mooring_info);
  
  status= mgr_init_formats();
  
  if (file_out_prefix) {
    sprintf(rootname, "%s", file_out_prefix);
    }
  else {
    cpos = strrchr(file_in, '/');
    if (!cpos) cpos = file_in;
    else cpos++;
    sprintf(rootname, "%s", cpos);
    if ((pos = strrchr(rootname,'.')))
      *pos = '\0';
    }
  

/*-----------------------------------------------------------------------------
  load time series */
  nseries = mgr_load_timeserie(file_in, &serie, 'm', format_in, ncolumns_GNU, header_format, &mooring);
  if(nseries<=0)
    return nseries;
  
  for(k=0; k<nseries; k++) {
    tseries_t *seriek=&serie[k];
    
    const double
      offset=
        shift/24./3600.
        +cnes_time(seriek->origin,'d');
    
    if(offset==0.)
      continue;
    
    for(int l=0;l<seriek->n;l++) {
      seriek->t[l]+=offset;
      }
    
    }
  
//  status=mgr_save_timeserie_NetCDF("test.nc", mooring, serie, nseries);

// /* *----------------------------------------------------------------------------
//   add mooring name to file out */
//   if (mooring.name) {
//     sprintf(mooring_name, "%s", mooring.name);
//     char *p=mooring_name;
//     while ((p=strchr(p,int(' '))))
//       p[0] = '_';
//     sprintf(rootname, "%s_%s", rootname, mooring_name);
//     }

/*------------------------------------------------------------------------------
  split series */
  reduced  = new tseries_t*[nseries];
  nreduced = new int[nseries];

/*------------------------------------------------------------------------------
  build k output series */
  int nb_tot_series=0;
  
  for(k=0; k<nseries; k++) {

/*------------------------------------------------------------------------------
    start build output file names */
    if (file_out_prefix && (nseries != 1) )
      sprintf(output, "%s.%4.4d", rootname, k);
    else
      sprintf(output, "%s", rootname);
    
/* *----------------------------------------------------------------------------
    add mooring name to file out */
    if (serie[k].mooring.name) {
      sprintf(mooring_name, "%s", serie[k].mooring.name);
      char *p=mooring_name;
      while ((p=strchr(p,int(' '))))
        p[0] = '_';
      sprintf(output, "%s_%s", output, mooring_name);
      }

/*------------------------------------------------------------------------------
    split or clip time series : get start / end dates of the serie */
    real_start = poctime_getdatecnes(serie[k].t[0], 'd');
    real_final = poctime_getdatecnes(serie[k].t[serie[k].n-1],'d');
    
/*------------------------------------------------------------------------------
    select starting date */
    if (date_start.year) {
      if (poctime_is_after(date_start, real_final)) {
        STDERR_BASE_LINE("analysis period starts (%s) after time series finishes (%s), no data for %s \n", sgetdate(date_start), sgetdate(real_final), file_in);
        return 0;
        }
      real_start = poctime_latest(real_start, date_start);
      }
/*------------------------------------------------------------------------------
    if no starting date but steps: then start at begining of the year ... */
    if (step.year) {
      real_start.month = 1;
      real_start.day = 1;
      /* always same years spliting */
      if (step.year!=1) {
        if (date_start.year) {
          int year_modulo;
          year_modulo = floor((real_start.year - date_start.year) / step.year);
          real_start.year = date_start.year + year_modulo *  step.year;
          }
        }
      }
/*------------------------------------------------------------------------------
     ... or of the month */
    else if (step.month) {
      real_start.day = 1;
      }

/*------------------------------------------------------------------------------
    select ending date */
    if (date_final.year) {
      if (poctime_is_before(date_final, real_start)) {
        STDERR_BASE_LINE("analysis period finishes (%s) after time series starts (%s), no data for %s \n", sgetdate(date_final), sgetdate(real_start), file_in);
        return 0;
        }
      real_final = poctime_soonest(real_final, date_final);
      }
/*------------------------------------------------------------------------------
     no ending date: go until end of year ... */
    if (step.year) {
      real_final.day = 31;
      real_final.month = 12;
      }
/*------------------------------------------------------------------------------
    ... or end of month */
    else if (step.month) {
      real_final.day = 31;
      }

/*------------------------------------------------------------------------------
    count number of split series */
    if (step.year || step.month)
      nreduced[k] = poctime_count_steps(real_start, real_final, step);
    else {
      nreduced[k] = 1;
      step.year = 100;
      }

/*------------------------------------------------------------------------------
    alloc time series */
    reduced[k] = new tseries_t[nreduced[k]];

/*------------------------------------------------------------------------------
    loop over time spliting */
    start = real_start;
    
    for(j = 0; j < nreduced[k]; j++) {
      tseries_t *reducedkj= &reduced[k][j];
      
/*------------------------------------------------------------------------------
      get interval */
      next = poctime_add_step(start, step, real_final, &final);
      printf("next %04d/%02d/%02d start %04d/%02d/%02d %fh final %04d/%02d/%02d %fh\n",
              next.year, next.month, next.day,
              start.year, start.month, start.day, start.second/3600.,
              final.year, final.month, final.day, final.second/3600.);

/*------------------------------------------------------------------------------
      extract */
      reducedkj->n = 0;
      *reducedkj = serie[k].reduce(start, final);

/*------------------------------------------------------------------------------
      jump if serie empty */
      if (reducedkj->n <= 1) {
        start = next;
        nreduced[k]--;
        STDOUT_BASE_LINE("Remove empty serie %s \n", file_in);
        continue;
        }

/*------------------------------------------------------------------------------
      resample */
      if(resampling>0){
        int count=reducedkj->resample_count(resampling);
        tseries_t resampled;
        struct timeval before;
        STDERR_BASE_LINE("resampling at %gs, giving %d time steps\n",resampling*86400.,count);
        gettimeofday(&before);
        resampled=reducedkj->resample(resampling,.5);
        STDERR_BASE_LINE("took %gs computation time\n",difftime(before));
        TRAP_ERR_EXIT(ENOEXEC,"only testing\n");
        resampled.destroy();
        }

/*------------------------------------------------------------------------------
      save data */
      if (file_out_prefix) {
        sprintf(output_dated, "%s_%04d%02d%02d_%04d%02d%02d", output,
                start.year, start.month, start.day,
                final.year, final.month, final.day);
        if(reconstruct==1)
        status=mgr_reconstruct(reducedkj, &mooring);
        status=mgr_save_timeserie(string(output_dated), reducedkj->mooring, *reducedkj, 'm', format_str, gnu_nice);
        }
      
/*------------------------------------------------------------------------------
      next */
      start = next;

      nb_tot_series += nreduced[k];
      }

    serie[k].destroy();
    } /* loop over k series */

/*------------------------------------------------------------------------------
  return series */
  if (out_series) {
    *out_series = new tseries_t*[nb_tot_series];
    n=0;
    for(k=0; k<nseries; k++) {
      for (j=0; j<nreduced[k]; j++) {
        (*out_series)[n] = &(reduced[k][j]);
        n++;
        }
      }
    nb_series = n;
    }

  return nb_series;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

string mgr_build_out_prefix(string f_in, string out_prefix, string *path, string *basename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///\brief build-up a prefix path
/**
\param f_in
\param out_prefix
\param *path set to the directory name. Without effect if \c NULL.
\param *basename set to the base name. Without effect if \c NULL.
\returns built-up prefix path
The built-up prefix path is made of the path and the base name.
The directory name is that of \c out_prefix, unless it is empty, in which case it is that of \c f_in.
The base name is that of \c out_prefix, unless it does not have one, in which case it is that of \c f_in.
*/
/*----------------------------------------------------------------------------*/
{
  string out;
  string base_f_in;
  size_t pos;

  /*  remove postfix from file_in  */
  base_f_in = f_in;
  pos = base_f_in.rfind(".");
  if (pos != string::npos) {
    base_f_in.erase(pos);
    }

  /* no out_prefix */
  if (out_prefix.empty()) {
    mgr_get_path_basename(base_f_in, path, basename);
    return base_f_in;
    }

  /* out_prefix is only a path */
  if ((*out_prefix.rbegin()) == '/') {
    /* remove path from base_f_in */
    mgr_get_path_basename(base_f_in, NULL, basename);
    if (path) *path = out_prefix;
    out=out_prefix;
    if(basename) out += *basename;
    return out;
    }

  /* out_prefix provides every things */
  mgr_get_path_basename(out_prefix, path, basename);
  out = out_prefix;
  return out;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void mgr_get_path_basename(string in, string *path, string *base)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///\brief get path and basename
/**
\param in directory or file path
\param *path set to the directory name. Without effect if \c NULL.
\param *base set to the file name. Without effect if \c NULL.
*/
/*----------------------------------------------------------------------------*/
{
  size_t pos;

  pos = in.rfind("/");
  if (pos == string::npos) {
    if (path) *path = "./";
    if (base) *base = in;
    return;
    }
  if (path){
    *path = in;
    path->erase(pos+1);
    }
  if (base) *base = in.substr(pos+1,string::npos);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_init_formats(void)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int chk;
  
  chk=timeseries_formats.size();
  if(chk>9)
    return 0;
  
  timeseries_formats[FORMAT_NAME_RMN]       =FORMAT_CODE_RMN;
  timeseries_formats[FORMAT_NAME_RMN]       =FORMAT_CODE_RMN;
  timeseries_formats[FORMAT_NAME_PAB]       =FORMAT_CODE_PAB;
  timeseries_formats[FORMAT_NAME_GNU]       =FORMAT_CODE_GNU;
  timeseries_formats[FORMAT_NAME_GLOSS]     =FORMAT_CODE_GLOSS;
  timeseries_formats[FORMAT_NAME_SONEL_HR]  =FORMAT_CODE_SONEL_HR;
  timeseries_formats[FORMAT_NAME_RADAR]     =FORMAT_CODE_RADAR;
  timeseries_formats[FORMAT_NAME_HAWAII]    =FORMAT_CODE_HAWAII;
  timeseries_formats[FORMAT_NAME_HAWAII_NC] =FORMAT_CODE_HAWAII_NC;
  timeseries_formats[FORMAT_NAME_FOREMAN]   =FORMAT_CODE_FOREMAN;
  timeseries_formats[FORMAT_NAME_PROFILERS] =FORMAT_CODE_PROFILERS;
  timeseries_formats[FORMAT_NAME_MINH]      =FORMAT_CODE_MINH;
  timeseries_formats[FORMAT_NAME_LIST2]     =FORMAT_CODE_LIST2;
  timeseries_formats[FORMAT_NAME_DMY]       =FORMAT_CODE_DMY;
  timeseries_formats[FORMAT_NAME_YMD]       =FORMAT_CODE_YMD;
  timeseries_formats[FORMAT_NAME_LEGEND]    =FORMAT_CODE_LEGEND;
  timeseries_formats[FORMAT_NAME_BODC]      =FORMAT_CODE_BODC;
  timeseries_formats[FORMAT_NAME_BODC2]     =FORMAT_CODE_BODC2;
  timeseries_formats[FORMAT_NAME_PUERTOS]   =FORMAT_CODE_PUERTOS;
  timeseries_formats[FORMAT_NAME_DART]      =FORMAT_CODE_DART;
  timeseries_formats[FORMAT_NAME_SEINE]     =FORMAT_CODE_SEINE;
  timeseries_formats[FORMAT_NAME_REFMAR]    =FORMAT_CODE_REFMAR;
  timeseries_formats[FORMAT_NAME_RADAR_RAW] =FORMAT_CODE_RADAR_RAW;
  timeseries_formats[FORMAT_NAME_NETCDF]    =FORMAT_CODE_NETCDF;
  timeseries_formats[FORMAT_NAME_CANADA]    =FORMAT_CODE_CANADA;
  timeseries_formats[FORMAT_NAME_USER]      =FORMAT_CODE_USER;
  timeseries_formats[FORMAT_NAME_LIST]      =FORMAT_CODE_LIST;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mooring_t mgr_decode_mooring_info(string header)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nitems;
  string s,dum;
  size_t pos;
  mooring_t mooring;

  s=header;

//rootname.erase(pos,4);

  if(s.c_str()==0) return(mooring);

  pos=s.find("--header=");
  if(pos!=string::npos) s.erase(pos,9);

  pos=s.find("lat=");
  if(pos!=string::npos) {
    dum=s.substr(pos+4,string::npos);
    nitems=sscanf(dum.c_str(),"%lf",&mooring.lat);
    }

  pos=s.find("lon=");
  if(pos!=string::npos) {
    dum=s.substr(pos+4,string::npos);
    nitems=sscanf(dum.c_str(),"%lf",&mooring.lon);
    }

  pos=s.find("depth=");
  if(pos!=string::npos) {
    dum=s.substr(pos+6,string::npos);
    nitems=sscanf(dum.c_str(),"%lf",&mooring.depth);
    }

  pos=s.find("code=");
  if(pos!=string::npos) {
    dum=s.substr(pos+5,string::npos);
    nitems=sscanf(dum.c_str(),"%d",&mooring.code);
    }

  pos=s.find("name=");
  if(pos!=string::npos) {
    dum=s.substr(pos+5,string::npos);
    mooring.name=new char[256];
    nitems=sscanf(dum.c_str(),"%s",mooring.name);
    }

  mooring.initialized=true;
  
  return(mooring);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_parse_mooring_header(const char *header_format,const char *line, mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  // "lon=lon lat=lat code=code name=name depth=depth"
  // "lon lat code name depth tmp tmp"
{
  char * pos_l;
  char * pos_h;
  char name[256];
  int res;

  if(!header_format) return 0;

  char *hf=strdup(header_format);
  char *l=strdup(line);
  pos_l = strtok (l," ");
  pos_h = strtok (hf," ");

  while (pos_l && pos_h) {
    
    if (strncmp("<lat>", pos_h, 5)==0)
      res = sscanf(pos_l, "%lf", &mooring->lat);
    else if (strncmp("<lon>", pos_h, 5)==0)
      res = sscanf(pos_l, "%lf", &mooring->lon);
    else if (strncmp("<depth>", pos_h, 7)==0)
      res = sscanf(pos_l, "%lf", &mooring->depth);
    else if (strncmp("<code>", pos_h, 6)==0)
      res = sscanf(pos_l, "%d", &mooring->code);
    else if (strncmp("<name>", pos_h, 6)==0) {
      res = sscanf(pos_l, "%s", &name);
      mooring->name = strdup(name);
      }
    pos_l = strtok (NULL," ");
    pos_h = strtok (NULL," ");
    }
  
  free(l);
  free(hf);
  
  return 1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_load_serie_nD(const char *filename, int & ncol, double* & time, double** & z, char unit,const char *header_format, mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, n,ndata,nitems;/* <index, number of lines, number of time steps */
  char line[300];
  FILE *in=NULL;
  char key[32];
  char * pos;
  
  in=fopen(filename,"r");
  if (in == NULL) TRAP_ERR_EXIT(-1,"unable to open file: %s \n",filename);

/*------------------------------------------------------------------------------
  get max number of data in file*/
  n=0;
  ndata=0;
  while(!feof(in)) {
    fgets(line,sizeof(line),in);
    n++;
    if(line[0]=='#') continue;
    ndata++;
    
    i=count_token(line, " \t");
    updatemax(&ncol,i-1);
    }
  n--;
  ndata--;
  
  rewind(in);

/*------------------------------------------------------------------------------
  allocate memory*/
  time=new double[ndata];
  z=new double*[ncol];
  for  (i=0; i<ncol; i++) {
    z[i]=   new double[ndata];
    }
  
/*------------------------------------------------------------------------------
  read data */
  ndata=0;
  for (i=0; i<n;i++) {
    fgets(line,sizeof(line),in);
/*------------------------------------------------------------------------------
    first line may be header */
    if (i==0 and (mooring==NULL or mooring->name==NULL or mooring->code == -1)) {
      if (header_format)
        mgr_parse_mooring_header(header_format, line, mooring);
      else if(line[0]=='#') {
        if (mooring->name==NULL) mooring->name=new char[256];
        nitems=sscanf(line, "%s %s %d %lf %lf %lf", key, mooring->name, &(mooring->code), &(mooring->lon),&(mooring->lat),&(mooring->depth));
	if(nitems==6) mooring->initialized=true;
        }
      }
/*------------------------------------------------------------------------------
    standard line */
    if(line[0]!='#') {
/*------------------------------------------------------------------------------
      read time */
      if(line[0]=='\n') continue;
      for(int k=0;k<300;k++) {
        if(line[k]=='\t') line[k]=' ';
        if(line[k]=='\n') break;
        }
      pos = strtok (line," ");
      int nitems=sscanf(pos, " %lf", &time[ndata]);
/*------------------------------------------------------------------------------
      loop over data */
      int j = 0;
      do {
        pos = strtok (NULL," ");
        if (pos==0) pos = strtok (NULL,"\t");
        if (pos==0) {
          TRAP_ERR_EXIT(-1,"** ERROR: bad format file: %s \n",filename);
          }
        nitems=sscanf(pos, " %lf", &z[j][ndata]);
        j++;
        } while (j<ncol);
      ndata++;
      }

    /* comment line, may be header ? */
//     else if (mooring==NULL || mooring->name==NULL || mooring->code == -1) {
//       if (header_format)
//         mgr_parse_mooring_header(header_format, line, mooring);
//       else {
//         int nb;
//         if (mooring->name==NULL)
//           mooring->name=new char[256];
//         nb=sscanf(line, "# %s %s %lf %lf %lf", key, mooring->name, &(mooring->lon),&(mooring->lat),&(mooring->depth));
//         if (nb != 5) {
//           nb=sscanf(line, "# %s %s %d %lf %lf %lf", key, mooring->name, &(mooring->code), &(mooring->lon),&(mooring->lat),&(mooring->depth));
//           if (nb != 6)
//             // sara: NOT POSSIBLE !?! we are in the case : line[0] == '#'
//             nb=sscanf(line, "%s %s %d %lf %lf %lf", key, mooring->name, &(mooring->code), &(mooring->lon),&(mooring->lat),&(mooring->depth));
//           if (nb != 6)
//             TRAP_ERR_EXIT(-1,"** ERROR: unknown format file: %s, line %s\n",filename,line);
//           }
//         }
//       }
    }

  fclose(in);
  return(ndata);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_line_count(string filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  std::ifstream file(filename.c_str());
  int lines;
  
  if ( file ){
    lines = std::count(std::istreambuf_iterator<char>( file ), std::istreambuf_iterator<char>(),'\n' );
    return(lines);
    }
  else return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void mgr_list_save_serie(const string & output,const tseries_t & serie,const mooring_t & mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *file = fopen(output.c_str(), "w");

  fprintf(file,"#-- HEADER -------------------------------------------\n");
  fprintf(file,"# Column 1 : date in days referred to\n");
  fprintf(file,"#            CNES date (01-JAN-1950 00:00:00.0)\n");
  fprintf(file,"# Column 2 : sea surface height (in meters)\n");
  fprintf(file,"#-- HEADER END ---------------------------------------\n");
  fprintf(file,"Number of crossover points :\n");
  fprintf(file,"%d\n",1);
  fprintf(file,"#-----------------------------------------------------\n");
  fprintf(file,"Pt  : 1\n");
  fprintf(file,"lon : %lf\n",mooring.lon);
  fprintf(file,"lat : %lf\n",mooring.lat);
  fprintf(file,"Mes : %d\n",serie.n);

  for(int i = 0; i < serie.n; i++){
    fprintf(file,"%13.6f %9.4f\n", serie.t[i], serie.x[1][i]);
    }

  fclose(file);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void mgr_list2_save_serie(const string & output,const tseries_t & serie)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *file = fopen(output.c_str(), "w");

  fprintf(file,"#-- HEADER -------------------------------------------\n");
  fprintf(file,"# Column 1 : date in days referred to\n");
  fprintf(file,"#            CNES date (01-JAN-1950 00:00:00.0)\n");
  fprintf(file,"# Column 2 : sea surface height (in meters)\n");
  fprintf(file,"#-- HEADER END ---------------------------------------\n");
  fprintf(file,"#lon : %lf\n",serie.lon);
  fprintf(file,"#lat : %lf\n",serie.lat);
  fprintf(file,"#Mes : %d\n",serie.n);

  for(int i = 0; i < serie.n; i++){
    fprintf(file,"%13.6f %9.4f\n", serie.t[i], serie.x[0][i]);
    }

  fclose(file);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

bool mgr_RMNinfo(string filename, double &lon, double &lat)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int deg;
  float min;
  char str_lat[20],str_lon[20];
  bool gotit = false;

  string shortname(filename);
  size_t pos = string::npos;
  if( (pos = shortname.find_last_of("/")) != string::npos){
    shortname = shortname.substr(pos + 1);
    }

  if(shortname=="ancona.csv") { sprintf(str_lon,"13 30.13"); sprintf(str_lat,"43 37.46"); gotit = true; }
  if(shortname=="bari.csv") { sprintf(str_lon,"16 52.00"); sprintf(str_lat,"41 08.24"); gotit = true; }
  if(shortname=="cagliari.csv") { sprintf(str_lon,"09 06.53"); sprintf(str_lat,"39 12.35"); gotit = true; }
  if(shortname=="carloforte.csv") { sprintf(str_lon,"08 18.29"); sprintf(str_lat,"39 08.37"); gotit = true; }
  if(shortname=="catania.csv") { sprintf(str_lon,"15 05.37"); sprintf(str_lat,"37 29.50"); gotit = true; }
  if(shortname=="civitavecchia.csv") { sprintf(str_lon,"11 47.00"); sprintf(str_lat,"42 05.22"); gotit = true; }
  if(shortname=="crotone.csv") { sprintf(str_lon,"17 08.09"); sprintf(str_lat,"39 04.41"); gotit = true; }
  if(shortname=="genova.csv") { sprintf(str_lon,"08 55.33"); sprintf(str_lat,"44 24.31"); gotit = true; }
  if(shortname=="imperia.csv") { sprintf(str_lon,"08 01.09"); sprintf(str_lat,"43 52.37"); gotit = true; }
  if(shortname=="lampedusa.csv") { sprintf(str_lon,"12 37.00"); sprintf(str_lat,"35 29.00"); gotit = true; }
  if(shortname=="livorno.csv") { sprintf(str_lon,"10 17.60"); sprintf(str_lat,"43 32.41"); gotit = true; }
  if(shortname=="messina.csv") { sprintf(str_lon,"15 33.54"); sprintf(str_lat,"38 11.21"); gotit = true; }
  if(shortname=="napoli.csv") { sprintf(str_lon,"14 16.09"); sprintf(str_lat,"40 50.23"); gotit = true; }
  if(shortname=="ortona.csv") { sprintf(str_lon,"14 24.56"); sprintf(str_lat,"42 21.20"); gotit = true; }
  if(shortname=="otranto.csv") { sprintf(str_lon,"18 29.49"); sprintf(str_lat,"40 08.47"); gotit = true; }
  if(shortname=="palermo.csv") { sprintf(str_lon,"13 22.14"); sprintf(str_lat,"38 07.14"); gotit = true; }
  if(shortname=="palinuro.csv") { sprintf(str_lon,"15 16.29"); sprintf(str_lat,"40 01.50"); gotit = true; }
  if(shortname=="pescara.csv") { sprintf(str_lon,"14 13.39"); sprintf(str_lat,"42 28.21"); gotit = true; }
  if(shortname=="portoempedocle.csv") { sprintf(str_lon,"13 31.28"); sprintf(str_lat,"37 17.25"); gotit = true; }
  if(shortname=="portotorres.csv") { sprintf(str_lon,"08 24.16"); sprintf(str_lat,"40 50.27"); gotit = true; }
  if(shortname=="ravenna.csv") { sprintf(str_lon,"12 16.48"); sprintf(str_lat,"44 29.48"); gotit = true; }
  if(shortname=="reggiocalabria.csv") { sprintf(str_lon,"15 38.57"); sprintf(str_lat,"38 07.14"); gotit = true; }
  if(shortname=="salerno.csv") { sprintf(str_lon,"14 44.54"); sprintf(str_lat,"40 40.53"); gotit = true; }
  if(shortname=="taranto.csv") { sprintf(str_lon,"17 13.29"); sprintf(str_lat,"40 28.31"); gotit = true; }
  if(shortname=="tremiti.csv") { sprintf(str_lon,"15 29.00"); sprintf(str_lat,"42 07.00"); gotit = true; }
  if(shortname=="trieste.csv") { sprintf(str_lon,"13 45.22"); sprintf(str_lat,"45 39.16"); gotit = true; }
  if(shortname=="venezia.csv") { sprintf(str_lon,"12 25.25"); sprintf(str_lat,"45 25.22"); gotit = true; }
  if(shortname=="vieste.csv") { sprintf(str_lon,"16 10.29"); sprintf(str_lat,"41 53.34"); gotit = true; }

  sscanf(str_lat,"%d %f",&deg,&min);
  lat=deg+min/60;

  sscanf(str_lon,"%d %f",&deg,&min);
  lon=deg+min/60;

  return(gotit);
}


/* Generated with :
sed -re 's/^ *(int mgr_load.+$)/extern \1;/;t;d' mgr-loaders.cpp
*/
extern int mgr_loadBODC(const char *filename, tseries_t **serie, mooring_t *mooring);
extern int mgr_loadPuertos(const char *filename, tseries_t **serie, mooring_t *mooring);
extern int mgr_loadDART(const char *filename, tseries_t **serie, mooring_t *mooring);
extern int mgr_loadRMN(char *filename);
extern int mgr_loadRadar(const char *filename, tseries_t **serie);
extern int mgr_loadRadarRaw(const char *filename, tseries_t **serie);
extern int mgr_loadList2(const char *filename, tseries_t **serie, mooring_t *mooring);
extern int mgr_loadProfilers(const char *filename, tseries_t **serie);
extern int mgr_loadPAB(const char *filename, tseries_t *serie);
extern int mgr_loadMINH(const char *filename, tseries_t *serie);
extern int mgr_loadSONEL_HR(const char *filename, tseries_t *serie, mooring_t *mooring);
extern int mgr_loadREFMAR(const char *filename, tseries_t *serie, mooring_t *mooring);
extern int mgr_loadBODC2(const char *filename, tseries_t *serie, mooring_t *mooring);
extern int mgr_loadOceanFisheriesCA(const char *filename, tseries_t *serie, mooring_t *mooring);
extern int mgr_loadUser(const char *filename, tseries_t *serie, mooring_t *mooring);
extern int mgr_loadDMY(const char *filename, tseries_t *serie);
extern int mgr_loadYMD(const char *filename, tseries_t *serie);
extern int mgr_loadSeine(const char *filename, tseries_t *serie, mooring_t *mooring);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_reconstruct(tseries_t *serie,const mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
 
  re-create missing values in time series
  
  Principles:
  
    low frequency contribution reconstructed from time series filtering
  
    high frequency contribution reconstructed from tidal prediction
    
  known to be failing in extreme conditions (estuarine flooding or draft)
  
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
  double *LF, window;
  date_t start,final;
  extern int remove_outlayers(double *t, double *h, double mask, double sampling, int n);
  extern double fill_gaps(double *t, double *h, double mask, int n, double w, bool fit);

  std::string line;
  cout << endl << "-------- starting reconstruction --------" << endl << endl;
  
  start=poctime_getdatecnes(serie->t[0], 'd');
  final=poctime_getdatecnes(serie->t[serie->n-1], 'd');
  
    
/*------------------------------------------------------------------------------
  sub-sample original series, to be done somewhere else */
//   tseries_t *reduced=new tseries_t;
//   int increment=1;
//   *reduced=serie->subsample(increment);
//   reduced->first = poctime_getdatecnes(reduced->t[0],'d');
//   reduced->last  = poctime_getdatecnes(reduced->t[reduced->n-1],'d');
//   *serie=*reduced;

/*------------------------------------------------------------------------------
  reduce original series, to be done somewhere else */
#if 0
  tseries_t *reduced=new tseries_t;
  *reduced=serie->reduce(poctime_getdatecnes(serie->t[0],'d'),poctime_getdatecnes(serie->t[0]+366,'d'));
  reduced->first = poctime_getdatecnes(reduced->t[0],'d');
  reduced->last  = poctime_getdatecnes(reduced->t[reduced->n-1],'d');
  *serie=*reduced;
#endif
  
  double *residuals=new double[serie->n];
  int nodal=1;
  double  mean;
  bool Z0=false;
//   spectrum_t solved, spectrum = spectrum_init_ref("COASTAL");
//   spectrum_t solved, spectrum = spectrum_init_ref("ESTUARINE-HF", Z0);
/*------------------------------------------------------------------------------
  long period needed because of non-linearities (Mf=M2-N2, MSf=M2-S2, etc...) */
  spectrum_t solved, spectrum = spectrum_init_ref("ESTUARINE", Z0);
  spectrum.remove(wSa);
  spectrum.remove(wSsa);
  mgr_data_t *tmp= harmonic_analysis_with_parsing(serie->x[0],serie->mask,serie->t,residuals, &mean, serie->n, spectrum, solved, nodal, stdout);

  mgr_t gauge;
  vector<mgr_t> mgr;
  
  gauge.number=0;
  if(mooring->name!=0)
    strcpy(gauge.name,mooring->name);
  else
    strcpy(gauge.name,"no_name");
  gauge.number = mooring->code;
  strcpy(gauge.origine,"undocumented");
  gauge.data=tmp;
  gauge.nwave=solved.n;
  gauge.loc.units=strdup("degrees");
  gauge.loc.lon=mooring->lon;
  gauge.loc.lat=mooring->lat;
  gauge.loc.depth=mooring->depth;
  gauge.mindex=1;
  gauge.duree=serie->t[serie->n-1]-serie->t[0];
  sprintf(gauge.debut,"%2.2d/%2.2d/%4d",start.day,start.month,start.year);
  sprintf(gauge.fin,  "%2.2d/%2.2d/%4d",final.day,final.month,final.year);
  mgr.push_back(gauge);
  status=mgr_save_ascii("test.mgr", mgr);
  
  if(tmp==0) return(-1);
  
  tseries_t *extended=new tseries_t;

  extended->nparam=11;
 
  extended->x=new double*[extended->nparam];
  extended->x[0]=serie->x[0];
  extended->x[1]=residuals;
  
  for(int k=2;k<extended->nparam;k++) {
    extended->x[k]=new double[serie->n];
    }

  extended->t=serie->t;
  extended->mask=serie->mask;
  extended->n=serie->n;
  extended->mooring=serie->mooring;
 
  *serie=*extended;
    
// /*------------------------------------------------------------------------------
//   smooth residuals */
//   status=loess1d_irregular(serie->n, 1.0, serie->x[1], serie->t, serie->mask,serie->x[2]);

#if 0
  for (int i=0; i< serie->n; i++) {
    int kmin=max(i-50,0);
    int kmax=min(i+50,serie->n);
//    statistic_t s=get_statistics(&(residuals[kmin]),serie->mask,kmax-kmin-1,0);
    statistic_t s=get_statistics(&(serie->x[2][kmin]),serie->mask,kmax-kmin-1,0);
    if(serie->x[0][i] != serie->mask) {
//      serie->x[2][i]=serie->x[0][i]-residuals[i]-mean;
      if(abs(residuals[i]-s.mean) > 6*s.std) {
        serie->x[0][i]=serie->mask;
        }
      }
    }
#endif

  int nRS=serie->size(1,serie->mask);
  int nLF=serie->size(2,serie->mask);
  
  for (int i=0; i< serie->n; i++) {
    if(serie->x[0][i] == serie->mask) {
      for(int k=1;k<serie->nparam;k++) {
        serie->x[k][i]=serie->mask;
        }
      }
    }
  
/*------------------------------------------------------------------------------
  identify alleged time sampling */
  double sampling=get_timesampling(serie->t, serie->mask, serie->n);
  sampling=floor(sampling*24.*3600.+0.5)/24/3600.;
  serie->sampling=sampling;
  printf("sampling appears to be %lf mn (%d data over %d tics,  %d missing\n",sampling*24.*60., serie->n, serie->maxsize(), serie->maxsize()-serie->n);
  
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  reconstruct time serie : resample time series at fixed interval
  
  interpolate missing values when gaps do not exceed a given duration
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("number of valid values before resampling: %d (over %d)\n",serie->n-serie->size(0,serie->mask),serie->n);
  tseries_t *reconstructed=new tseries_t;
  
/*------------------------------------------------------------------------------
  do not re-interpolate if 2 successive time larger than max_gap (assign mask)*/
  double max_gap=1./24/1.+0.01;
  *reconstructed=serie->resample(sampling, max_gap);
  
  *serie=*reconstructed;
  printf("number of valid values after  resampling: %d (over %d),  %d missing, max gap=%lf hours \n",serie->n-serie->size(0,serie->mask),serie->n, serie->size(0,serie->mask), max_gap*24.);
  
// /*------------------------------------------------------------------------------
//   try Demerliac filter to isolate residual tidal signal */
//   serie->x[5] = tide_filters("demerliac", sampling, serie->x[1], serie->mask, serie->n);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  try to remove observation outlayers 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   int nmasked_prior=serie->size(1,serie->mask);

  status=remove_outlayers(serie->t, serie->x[1], serie->mask, sampling, serie->n);
  
//   int nmasked_posterior=serie->size(1,serie->mask);
  
//   tmp= harmonic_analysis_with_parsing(serie->x[1],serie->mask,serie->t,residuals, &mean, serie->n, spectrum, solved, nodal, stdout);
//   
//   gauge.number=1;
//   if(mooring->name!=0)
//     strcpy(gauge.name,mooring->name);
//   else
//     strcpy(gauge.name,"no_name");
//   gauge.number = mooring->code;
//   strcpy(gauge.origine,"undocumented");
//   gauge.data=tmp;
//   gauge.nwave=solved.n;
//   gauge.loc.units=strdup("degrees");
//   gauge.loc.lon=mooring->lon;
//   gauge.loc.lat=mooring->lat;
//   gauge.loc.depth=mooring->depth;
//   gauge.mindex=1;
//   gauge.duree=serie->t[serie->n-1]-serie->t[0];
//   sprintf(gauge.debut,"%2.2d/%2.2d/%4d",start.day,start.month,start.year);
//   sprintf(gauge.fin,  "%2.2d/%2.2d/%4d",final.day,final.month,final.year);
//   mgr.push_back(gauge);
//   status=mgr_save_ascii("test.mgr", mgr);
//   
//   if(tmp==0) return(-1);

/*------------------------------------------------------------------------------
  smooth residuals */
  window=0.125/2;
  status=loess1d_irregular(serie->n, window, serie->x[1], serie->t, serie->mask, serie->x[2]);
  for (int i=0; i< serie->n; i++) {
    serie->x[8][i]=serie->x[1][i];
    if(serie->x[1][i] != serie->mask) serie->x[2][i]=serie->x[1][i];
    }

  nLF=serie->size(2,serie->mask);
  
  LF=new double[serie->n];
  
  window=0.;
  while(nLF!=0) {
    window+=0.125/2.0;
    status=loess1d_irregular(serie->n, window, serie->x[2], serie->t, serie->mask, LF);
    for (int i=0; i< serie->n; i++) {
      if(serie->x[2][i] == serie->mask) serie->x[2][i]=LF[i];
      }
    nLF=serie->size(2,serie->mask);
    }

/*------------------------------------------------------------------------------
  try Demerliac filter to isolate residual tidal signal */
  serie->x[5] = tide_filters("demerliac", sampling, serie->x[2], serie->mask, serie->n);
  
  status=harmonic_predictionTS(serie->x[3], serie->t, serie->n, tmp, solved.n, nodal);

//   serie->x[5] = tide_filters("demerliac", sampling, serie->x[1], serie->mask, serie->n);
//   for (int i=0; i< serie->n; i++) {
//     if(serie->x[1][i] != serie->mask and serie->x[5][i] != serie->mask) {
//       double y=serie->x[3][i]+serie->x[1][i]-serie->x[5][i];
//       serie->x[8][i]=y*y;
//       }
//     else {
//       serie->x[8][i]=serie->mask;
//       }
//     }
//   window=10.0;
//   status=loess1d_irregular(serie->n, window, serie->x[8], serie->t, serie->mask, serie->x[9]);
//   for (int i=0; i< serie->n; i++) {
//     if(serie->x[3][i] != serie->mask) {
//       double y=serie->x[3][i];
//       serie->x[8][i]=y*y;
//       }
//     else {
//       serie->x[8][i]=serie->mask;
//       }
//     }
//   window=10.0;
//   status=loess1d_irregular(serie->n, window, serie->x[8], serie->t, serie->mask, serie->x[10]);
//   for (int i=0; i< serie->n; i++) {
//     if(serie->x[9][i] != serie->mask and serie->x[10][i] != serie->mask) {
//       double y=serie->x[10][i]/serie->x[9][i];
//       serie->x[8][i]=sqrt(y);
//       }
//     else {
//       serie->x[8][i]=0;
//       }
//     if(serie->x[9][i]  == serie->mask) serie->x[9][i] =0.0;
//     if(serie->x[10][i] == serie->mask) serie->x[10][i]=0.0;
//     }
  
  nRS=serie->size(1,serie->mask);
  nLF=serie->size(2,serie->mask);
  
  status=fill_gaps(serie->t, serie->x[8], serie->mask, serie->n, 1.0, false);
  
  for (int i=0; i< serie->n; i++) {
    if(serie->x[0][i] == serie->mask) serie->x[0][i]=mean;
    if(serie->x[1][i] == serie->mask) serie->x[1][i]=0;
    if(serie->x[0][i] == serie->mask) {
      serie->x[0][i]=serie->x[3][i]+mean;
      if(serie->x[2][i]!=serie->mask) serie->x[0][i]+=serie->x[2][i];
      serie->x[1][i]=0;
      if(serie->x[2][i]!=serie->mask) serie->x[1][i]+=serie->x[2][i];
      }
    if(serie->x[1][i] == serie->mask) {
      serie->x[1][i]=0;
      if(serie->x[2][i]!=serie->mask) serie->x[1][i]+=serie->x[2][i];
      }
    if(serie->x[2][i] == serie->mask) {
      serie->x[2][i]=0;
      }
/*------------------------------------------------------------------------------
    reconstruction */
    serie->x[4][i]=serie->x[3][i]+serie->x[2][i]+mean;
/*------------------------------------------------------------------------------
    reconstruction */
    serie->x[8][i]+=serie->x[3][i]+mean;
    if(serie->x[5][i] == serie->mask) {
      serie->x[5][i]=0.;
      serie->x[6][i]=0.;
      serie->x[7][i]=0.;
      }
    else {
/*------------------------------------------------------------------------------
      non-harmonic tidal signal; mean added to ease plots */
      serie->x[6][i]=serie->x[2][i]-serie->x[5][i]+mean;
/*------------------------------------------------------------------------------
      total tidal signal */
      serie->x[7][i]=serie->x[3][i]+serie->x[2][i]-serie->x[5][i];
      }
    }
    
  
  printf("number of valid values after reconstruction: %d (over %d)\n",serie->n-serie->size(3,serie->mask),serie->n);
  printf("start %04d/%02d/%02d %fh final %04d/%02d/%02d %fh\n", start.year, start.month, start.day, start.second/3600.,
                                                                final.year, final.month, final.day, final.second/3600.);

  printf("0 : original; 1 : harmonic tidal residuals; 2 : low frequency; 3 : harmonic tides; 4 : recontructed;\
          5 : residuals; 6 : non-harmonic tide; 7 : cumulated tide; 8 : cumulated tide\n");
  printf("mean=%lf window=%lf\n",mean,window);

  return(1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_decode_position(const char *filename, const char *format, const char *header_format, mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int fmt, nseries=0;
  int status;

  status= mgr_init_formats();
  
/*------------------------------------------------------------------------------
  get format code from format name */
  fmt=timeseries_formats[format];
  
  switch(fmt) {
    case FORMAT_CODE_PAB:
      status=-1;
      break;

    case FORMAT_CODE_GNU:
      status=-1;
      break;
      
    case FORMAT_CODE_GLOSS:
    case FORMAT_CODE_HAWAII:
      nseries=mgr_loadGLOSS( filename, 0, mooring);
      break;

    case FORMAT_CODE_HAWAII_NC:
      status=mgr_loadHAWAII_NC( filename, 0, mooring);
      break;

    case FORMAT_CODE_SONEL_HR:
      status=-1;
      break;

    case FORMAT_CODE_RADAR:
      status=-1;
      break;
      
    case FORMAT_CODE_PROFILERS:
      status=-1;
      break;

    case FORMAT_CODE_MINH:
      status=-1;
      break;
      
    case FORMAT_CODE_DMY:
      status=-1;
      break;

    case FORMAT_CODE_YMD:
      status=-1;
      break;

    case FORMAT_CODE_BODC:
      status=-1;
      break;

    case FORMAT_CODE_DART:
      status=-1;
      break;

    case FORMAT_CODE_SEINE:
      status=-1;
      break;

    case FORMAT_CODE_REFMAR:
      nseries=mgr_loadREFMAR(filename, 0, mooring);
      break;

    case FORMAT_CODE_CANADA:
      nseries=mgr_loadOceanFisheriesCA(filename, 0, mooring);
      break;

    default:
      TRAP_ERR_EXIT(-1, "input format %s unknown by %s, but %s is poorly updated although easy to update!\n", format, __func__, __func__);
    }
    
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int mgr_ParseFormat(const char *filename, const char *format_in)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int fmt;
  const char *default_format="GNU";
  const char *format=NULL;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  parse input format setting
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (format_in!=NULL) {
/*------------------------------------------------------------------------------
    format provided */
    format = format_in;
    }    
  else {
/*------------------------------------------------------------------------------
    try to guess format from file suffix  */
    if (strrncmp(filename,  ".gnu") == 0)
      format = FORMAT_NAME_GNU;
    else if (strrncmp(filename,  ".list") == 0)
      format = FORMAT_NAME_LIST2;
    else if (strrncmp(filename,  ".dfo") == 0)
      format = FORMAT_NAME_FOREMAN;
    else if (strrncmp(filename,  ".dat") == 0)
      format = FORMAT_NAME_GLOSS;
    else if (strrncmp(filename,  ".nc") == 0)
      format = FORMAT_NAME_HAWAII_NC;
    }
  
  if (!format) {
    STDOUT_BASE_LINE("Unknown file format, try default GNU format. To change format use option -iF");
    format = default_format;
    }
  
/*------------------------------------------------------------------------------
  get format code from format name */
  mgr_init_formats();
  fmt=timeseries_formats[format];

  return(fmt);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int mgr_load_timeserie(const char *filename, tseries_t **serie, char unit, const char *format_in, int & ncolumns, const char *header_format, mooring_t *mooring, int *status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int fmt, n, i, nseries=0;
  const char *default_format="GNU";
  const char *format=NULL;
  int status_;
  
  if(status==0)
    status=&status_;
  else{
    if(access(filename,R_OK)!=0){
      *status=errno;
      return 0;
      }
    *status=0;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  parse input format setting
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  fmt=mgr_ParseFormat(filename, format_in);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  process file
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  switch(fmt) {
    case FORMAT_CODE_PAB:
      *serie=new tseries_t[1];
      (*serie)[0].x=new double*[1];
      (*serie)[0].nparam=1;
      nseries=1;
      n=mgr_loadPAB(filename, &(*serie)[0]);
      (*serie)[0].mooring = *mooring;
      break;

    case FORMAT_CODE_LIST2:
      nseries=mgr_loadList2(filename, serie, mooring);
      /* copy same mooring in each series... */
//       for (i=0;i<nseries;i++)
//         (*serie)[i].mooring = *mooring;
//       break;
      break;
      
    case FORMAT_CODE_GNU:
      *serie=new tseries_t[1];
//       (*serie)[0].x=new double*[ncolumns];
//       (*serie)[0].nparam=ncolumns;
      serie[0]->n=mgr_load_serie_nD(filename, ncolumns, (*serie)[0].t, (*serie)[0].x, 'm', header_format, mooring);
      (*serie)[0].mooring = *mooring;
      (*serie)[0].nparam=ncolumns;
      nseries=1;
      break;
      
    case FORMAT_CODE_GLOSS:
    case FORMAT_CODE_HAWAII:
      *serie=new tseries_t[1];
      (*serie)[0].x=new double*[1];
      (*serie)[0].nparam=1;
      (*serie)[0].mask=9999.0;
      (*serie)[0].n=mgr_loadGLOSS(filename, &(*serie)[0], mooring);
//      serie[0]->n=mgr_loadGLOSS_local(filename, date_t first, date_t last, mooring_t *mooring, double **elevation, double **time, double *mask, int *n, char *outfile)
//      serie[0]->n=hawaii_load(filename, mooring_t *mooring, double **elevation, double **time, double *mask, int *n, char *outfile)
      (*serie)[0].mooring = *mooring;
      nseries=1;
      break;

    case FORMAT_CODE_HAWAII_NC:
      *serie=new tseries_t[1];
      *status=mgr_loadHAWAII_NC(filename, &(*serie)[0], mooring);
      (*serie)[0].mooring = *mooring;
      nseries=1;
      break;

    case FORMAT_CODE_SONEL_HR:
      *serie=new tseries_t[1];
      (*serie)[0].x=new double*[1];
      (*serie)[0].nparam=1;
      nseries=1;
      n=mgr_loadSONEL_HR(filename, &(*serie)[0], mooring);
      (*serie)[0].mooring = *mooring;
      break;

    case FORMAT_CODE_RADAR:
      nseries=mgr_loadRadar(filename, serie);
      /* copy same mooring in each series... */
      for (i=0;i<nseries;i++)
        (*serie)[i].mooring = *mooring;
      break;
      
    case FORMAT_CODE_RADAR_RAW:
      nseries=mgr_loadRadarRaw(filename, serie);
      /* copy same mooring in each series... */
      for (i=0;i<nseries;i++)
        (*serie)[i].mooring = *mooring;
      break;
      
    case FORMAT_CODE_PROFILERS:
      nseries=mgr_loadProfilers(filename, serie);
      /* copy same mooring in each series... */
      for (i=0;i<nseries;i++)
        (*serie)[i].mooring = *mooring;
      break;

    case FORMAT_CODE_MINH:
      *serie=new tseries_t[1];
//      (*serie)[0].x=new double*[1];
      nseries=1;
      n=mgr_loadMINH(filename, &(*serie)[0]);
      (*serie)[0].mooring = *mooring;
      break;
      
    case FORMAT_CODE_DMY:
      *serie=new tseries_t[1];
      (*serie)[0].x=new double*[1];
      (*serie)[0].nparam=1;
      nseries=1;
      n=mgr_loadDMY(filename, &(*serie)[0]);
//       mooring->lon=0.0;
//       mooring->lat=0.0;
      (*serie)[0].mooring = *mooring;
      break;

    case FORMAT_CODE_YMD:
      *serie=new tseries_t[1];
      (*serie)[0].x=new double*[1];
      (*serie)[0].nparam=1;
      nseries=1;
//       mooring->lon=0.0;
//       mooring->lat=0.0;
      n=mgr_loadYMD(filename, &(*serie)[0]);
      if(n==0) nseries=0;
      (*serie)[0].mooring = *mooring;
      break;

    case FORMAT_CODE_BODC:
      *serie=new tseries_t[1];
      (*serie)[0].x=new double*[1];
      (*serie)[0].nparam=1;
      nseries=mgr_loadBODC(filename, serie, mooring);
      (*serie)[0].mooring = *mooring;
      break;

    case FORMAT_CODE_BODC2:
      *serie=new tseries_t[1];
      (*serie)[0].x=new double*[1];
      (*serie)[0].nparam=1;
      nseries=1;
      n=mgr_loadBODC2(filename, &(*serie)[0], mooring);
      (*serie)[0].mooring = *mooring;
      break;

    case FORMAT_CODE_PUERTOS:
      *serie=new tseries_t[1];
      (*serie)[0].x=new double*[1];
      (*serie)[0].nparam=1;
      nseries=mgr_loadPuertos(filename, serie, mooring);
      (*serie)[0].mooring = *mooring;
      break;

    case FORMAT_CODE_DART:
      *serie=new tseries_t[1];
      (*serie)[0].x=new double*[1];
      (*serie)[0].nparam=1;
      nseries=1;
      n=mgr_loadDART(filename, serie, mooring);
      (*serie)[0].mooring = *mooring;
      break;

    case FORMAT_CODE_SEINE:
      *serie=new tseries_t[1];
      (*serie)[0].x=new double*[1];
      (*serie)[0].nparam=1;
      nseries=mgr_loadSeine(filename, &(*serie)[0], mooring);
      (*serie)[0].mooring = *mooring;
      break;

    case FORMAT_CODE_REFMAR:
      *serie=new tseries_t[1];
      (*serie)[0].x=new double*[1];
      (*serie)[0].nparam=1;
      nseries=1;
      n=mgr_loadREFMAR(filename, &(*serie)[0], mooring);
      (*serie)[0].mooring = *mooring;
      break;

    case FORMAT_CODE_CANADA:
      *serie=new tseries_t[1];
      (*serie)[0].x=new double*[1];
      (*serie)[0].nparam=1;
      nseries=1;
      n=mgr_loadOceanFisheriesCA(filename, &(*serie)[0], mooring);
      (*serie)[0].mooring = *mooring;
      break;

    case FORMAT_CODE_USER:
      *serie=new tseries_t[1];
      (*serie)[0].x=new double*[1];
      (*serie)[0].nparam=1;
      nseries=1;
      n=mgr_loadUser(filename, &(*serie)[0], mooring);
      (*serie)[0].mooring = *mooring;
      break;

    default:
      TRAP_ERR_EXIT(-1, "input format unknown");
    }

/*------------------------------------------------------------------------------
  replace ' ' by '_' in mooring name */
  char *p=(*serie)[0].mooring.name;
  if(p!=0) {
    while ((p=strchr(p,int(' ')))) {
      p[0] = '_';
      }
    }
  
  return(nseries);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_saveforeman(const char *filename,const tseries_t & serie, char unit)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

(i) One card for the variables IOUT1,RAYOPT,ZOFF,ICHK,OBSFAC,INDPR,NSTRP in the format
    (I2,2X,F4.2,2X,F10.0,I2,3X,F10.7,215).
IOUT1  = 6 if the only output desired is a line printer listing of results,
       = 2 if both analysis output and listing are desired;
RAYOPT = Rayleigh criterion constant value if different from 1.0;
ZOFF   = constant to be subtracted from all the hourly heights;
ICHK   = 0 if the hourly height input data is to be checked for format errors,
       = otherwise if this checking to be waived;
OBSFAC = scaling factor, if different from 0.01, which will multiply the hourly
         observations, in order to produce the desired units for the final constituent
         amplitudes. (e.g. if the hourly observations are in mm/s and the final units
         are to be ft/sec, then this variable would be set to 0.0032808.);
INDPR = 1 if hourly height predictions based on the analysis results are to be
        calculated and written onto device number 10. If there is inference, this
        parameter value will also give the rms residual error after inference
        adjustments have been made,
      = 0 if no such predictions are desired;
NSTRP = number of successive moving average filters that have been applied to the
        original data.

(ii) One card for each possible inference pair. The format is (2(4X,A5,E16.10),2F10.3) and
     the respective values read are:
KONAN & SIGAN = name and frequency of the analysed constituent to be used for the
                inference;
KONIN & SIGIN = name and frequency of the inferred constituent;
                R = amplitude ratio of KONIN to KONAN;
ZETA          = Greenwich phase lag of the inferred constituent subtracted from the
                Greenwich phase lag of the analysed constituent.
These are terminated by one blank card.

(iii) One card for each shallow water constituent, other than those in the standard 69 constituent
data package, to be considered for inclusion in the analysis. The Rayleigh comparison
constituent is also required and the additional shallow water constituent must be found in
data type (i) of logical unit 8, but have a blank data field where the Rayleigh comparison
constituent is expected. The format is (6X,A5,4X,A5) and a blank card is required at the
end.

(iv) One card in the format (I1,1X,10I2) specifying the following information on the period
of the analysis:
INDY = 8 indicates an analysis is desired for the upcoming period;
= 0 indicates no further analyses are required;
IHH1,IDD1,IMM1,IYY1,ICC1 = hour, day, month, year and century of the beginning of
the analysis (measured in time ITZONE of input data (v));
IHHL,IDDL,IMML,IYYL,ICCL = hour, day, month, year and century of the end of the
analysis.

(v) One card in the format (I1,4X,A5,3A6,A4,A3,1X,2I2,I3,I2,5X,A5) containing the fol-
lowing information on the tidal station:
INDIC = 1 if J card output is desired (no longer used),
      = otherwise if not;
KSTN = tidal station number;
(NA(J),J=1,4) = tidal station name (22 characters maximum length);
ITZONE = time zone of the hourly observations;
LAD,LAM = station latitude in degrees and minutes;
LOD,LOM = station longitude in degrees and minutes;
IREF = reference station number.

(vi) The hourly height data cards contain the following information in the format
(I1,1X,I5,4X,I2,1X,3I2,12A4).

------------------------------------------------------------------------------*/
{
  int i,k,l,ll,nvalues,ndays;
  int card,century,id=100;
  int lad,lam,lod,lom;
  int mask=9999;
  double t0,t;
  double hmin=+1.e+10;
  int h[24];
  date_t first,last,instant;
  FILE *hawaii=fopen(filename,"w");

  poctime_convert(serie.t,serie.n,'d','h');

  for (i=0; i< serie.n; i++) {
    if(serie.x[0][i]!=serie.mask) {
//      serie.x[i]=serie.x[i]*100.0;
      updatemin(&hmin,serie.x[0][i]);
      }
    }

  for (i=0; i< serie.n; i++) {
    if(serie.x[0][i]!=serie.mask) {
      serie.x[0][i]-=hmin;
      }
    }

  first=poctime_getdatecnes(serie.t[0], 'h');
  first.second=0;

  t0=cnes_time(first,'h');

  last=poctime_getdatecnes(serie.t[serie.n-1], 'h');
  last.second=24*3600.;

/*------------------------------------------------------------------------------
  line #1*/
  int    IOUT1=2,ICHK=0,INDPR=1,NSTRP=0;
  float  RAYOPT=1.0,ZOFF=0.0,OBSFAC=1.00;
  fprintf(hawaii,"%2d  %4.2f  %10.0f%2d   %10.7f%5d%5d\n", IOUT1,RAYOPT,ZOFF,ICHK,OBSFAC,INDPR,NSTRP);
//    (I2,2X,F4.2,2X,F10.0,I2,3X,F10.7,215).

/*------------------------------------------------------------------------------
  line #2*/
//  fprintf(hawaii,"\n");

/*------------------------------------------------------------------------------
  line #3*/
  char *line;
  line=strdup("    K1       0.0417807462    P1       0.0415525871   0.33093   -7.07");
  fprintf(hawaii,"%s\n",line);
  line=strdup("    S2       0.0833333333    K2       0.0835614924   0.27215   -22.40");
  fprintf(hawaii,"%s\n",line);
  fprintf(hawaii,"\n");

/*------------------------------------------------------------------------------
  line #4*/
  line=strdup("      M10      M8");
  fprintf(hawaii,"%s\n",line);
  fprintf(hawaii,"\n");

/*------------------------------------------------------------------------------
  line #5*/
  card=8;
  fprintf(hawaii,"%1d ", card);
  century=first.year/100;
  fprintf(hawaii,"%2d%2d%2d%2d%2d", (int) NINT(first.second/3600.),first.day,first.month,first.year-100*century,century);
  fprintf(hawaii,"%2d%2d%2d%2d%2d", (int) NINT(last.second/3600.),last.day,last.month,last.year-100*century,century);
  fprintf(hawaii,"\n");

/*------------------------------------------------------------------------------
  line #6*/
  lad=lam=lod=lom=0;
  card=1;
  fprintf(hawaii,"%1d    %5d%18s%4s%3s %2d%2d%3d%2d     %5s\n", card,id,"name","----"," +0",lad,lam,lod,lom,"ref.");

/*------------------------------------------------------------------------------
  line #(data)*/
  nvalues=(int) floor(ellapsed_time(first,last,'h')+0.5);
  ndays=nvalues/24;

  i=0;
  for(k=0;k<ndays;k++) {
    for(l=0;l<24;l++) h[l]=mask;
    t=t0+k*24.+1; /**start at 1 am or pm*/
    instant=poctime_getdatecnes(t, 'h');
    while((serie.t[i]-t+1.)<1.e-03) {
      i++;
      }
    for(l=0;l<24;l++) {
      if(i+l==serie.n) break;
      ll=floor(serie.t[i+l]-t+0.5);
      if(ll>=24) break;
      h[ll]=NINT(serie.x[0][i+l]);
      }
    century=instant.year/100;
    card=1;
    fprintf(hawaii,"%1d %5d    %2d %2d%2d%2d", card,id,century,instant.day,instant.month,instant.year-100*century);
    for(l=0;l<12;l++) fprintf(hawaii,"%4d",h[l]);
    fprintf(hawaii,"\n");
    card=2;
    fprintf(hawaii,"%1d %5d    %2d %2d%2d%2d", card,id,century,instant.day,instant.month,instant.year-100*century);
    for(l=12;l<24;l++) fprintf(hawaii,"%4d",h[l]);
    fprintf(hawaii,"\n");
    }

  card=0;
  fprintf(hawaii,"%1d \n", card);
  fclose(hawaii);
  poctime_convert(serie.t,serie.n,'h','d');
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_save_timeserie_gnu(const char *filename,const tseries_t & serie,const mooring_t & mooring, char unit, bool gnu_safe, string time_template, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k=0;
  int i;
  double time,startd;
  bool add_line=gnu_safe;
  FILE *out;
  vector<string> token;
  char *startc;

/*
  mooring->lon=serie.lon;
  mooring->lat=serie.lat;
  mooring->depth=0.;
*/
  if(time_template!="") {
    token=string_split((const  string) time_template, " ");
    }
    
  double sampling=get_timesampling(serie.t, serie.mask, serie.n);
  
  sampling=floor(sampling*24.*3600.+0.5)/24/3600.;

  startd=cnes_time(serie.origin, 'd');

  out=fopen(filename,"w");
  
  if(mooring.initialized) {
    fprintf(out,"# %s %d %9.3lf %9.3lf %9.3f %d %s\n",
            mooring.name, mooring.code,
            mooring.lon, mooring.lat,
            mooring.depth, serie.n, filename);
    }
  else {
    fprintf(out,"# <NAME> <CODE> <LON> <LAT> <DEPTH> <NDATA> <FILENAME>\n");
    }
    
  for (i=0; i< serie.n; i++) {
    bool masked=true;
    for(k=0;k<serie.nparam;k++) {
       if(serie.x[k][i]!=serie.mask) {
         masked=false;
         break;
         }
      }
    time=serie.t[i];
    if(not masked) {
      if(i>0) {
        if(fabs(time-serie.t[i-1]-sampling)*d2s > 10.0) {
          if(add_line) {
            fprintf(out,"\n");
            if(verbose==1) fprintf(stderr,"time break at: %6d %lf %s\n",i, time, sgetcnesdate(time*24.));
            }
          add_line=0;
          }
        }
      if(token.size()==0) {
        fprintf(out,"%12.6lf",time);
        }
      else {
        for(k=0;k<token.size();k++) {
          if(token[k]=="CNES") {
//             fprintf(out,"%.6f ",(time+startd)/d2s);
            fprintf(out,"%.6f ",time+startd);
            }
          if(token[k]=="CALENDAR") {
            startc=poctime_sdate_cnes(time+startd,serie.time_unit,' ');
            fprintf(out,"%s ",startc);
            free(startc);
            }
          }
        }
      for(k=0;k<serie.nparam;k++) {
        fprintf(out," %7.4lf",serie.x[k][i]);
        }
      fprintf(out,"\n");
      if(gnu_safe) add_line=1;
      }
    else {
      if(add_line) {
        time=serie.t[i];
        fprintf(out,"\n");
        fprintf(stderr,"masked data break at: %6d %lf %s\n",i, time, sgetcnesdate(time*24.));
        }
      add_line=0;
      continue;
      }
    }
  
  fclose(out);
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_save_timeserie_NetCDF(const string & filename,const tseries_t *serie, int nseries)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status = NC_NOERR + 1;
 
  poc_dim_t nframes, npoints;
  
  poc_list_t<poc_dim_t> dimX, dimT, dimXT;

  poc_global_t info;
  
  status=poc_create(filename.c_str(), info, 1);
  if( status != NC_NOERR) NC_TRAP_ERROR(wexit,status,1,"poc_create(\""+filename+"\",) failed");
  
  nframes=poc_dim_t("T",serie[0].n);
  npoints=poc_dim_t("X",nseries);
  
  dimX.push_back(npoints);
  dimT.push_back(nframes);

  dimXT.push_back(npoints);
  dimXT.push_back(nframes);
 
/**----------------------------------------------------------------------------
  put longitudes */
  poc_data_t<double> lon;
  lon.info.init("lon",NC_DOUBLE,"","degrees",serie[0].mask).dimensions=dimX;
  lon.init("");
  
  lon.data=new double[nseries];
  for(int m=0;m<nseries;m++) lon.data[m]=serie[m].lon;
  
  lon.write_data(filename);

/**----------------------------------------------------------------------------
  put latitudes */
  poc_data_t<double> lat;
  lat.info.init("lat",NC_DOUBLE,"","degrees",serie[0].mask).dimensions=dimX;
  lat.init("");
  
  lat.data=new double[nseries];
  for(int m=0;m<nseries;m++) lat.data[m]=serie[m].lat;
  
  lat.write_data(filename);

/**----------------------------------------------------------------------------
  put time */
  poc_data_t<double> time;
  time.info.init("time",NC_DOUBLE,"","days since 1950/01/01 00:00:00",serie[0].mask).dimensions=dimT;
  time.init("");
  
  deletep(&time.data);
  time.data=new double[serie[0].n];
  for(int m=0;m<serie[0].n;m++) time.data[m]=serie[0].t[m];
  
  time.write_data(filename,-2);

/**----------------------------------------------------------------------------
  put data */
  poc_data_t<double> u;
  u.info.init("",NC_FLOAT,"","m",serie[0].mask).dimensions=dimXT;
  if(serie->nparam==1)
    u.info.name="elevation";
  else
    u.info.name="u";
  u.init("");
  
  u.data=new double[nseries*serie[0].n];
  for(int m=0;m<nseries;m++) {
    for(int n=0;n<serie[m].n;n++) {
      u.data[m*serie[m].n+n]=serie[m].x[0][n];
      }
    }
  u.write_data(filename);

/**----------------------------------------------------------------------------
  put data */
  if(serie->nparam>=2){
  poc_data_t<double> v;
  v.info.init("v",NC_FLOAT,"","m",serie[0].mask).dimensions=dimXT;
  v.init("");
  
  v.data=new double[nseries*serie[0].n];
  for(int m=0;m<nseries;m++) {
    for(int n=0;n<serie[m].n;n++) {
      v.data[m*serie[m].n+n]=serie[m].x[1][n];
      }
    }
  v.write_data(filename);
  }

/**----------------------------------------------------------------------------
  put data */
  poc_data_t<double> valids;
  valids.info.init("nvalids",NC_DOUBLE,"","",serie[0].mask).dimensions=dimX;
  valids.init("");
  
  valids.data=new double[nseries];
  for(int m=0;m<nseries;m++) {
    int nmasked=serie[m].size(0,serie->mask);
    valids.data[m]=serie[m].n-nmasked;
    }
  valids.write_data(filename);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_save_timeserie(const string & rootname,const mooring_t & mooring, const tseries_t & serie, char unit,const string & format, bool gnu_nice, string time_format, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int fmt, status;
  string filename="";
  
  fmt=timeseries_formats[format.c_str()];

  switch(fmt) {

    case FORMAT_CODE_GNU:
      filename = rootname+ ".gnu";
      status=mgr_save_timeserie_gnu(filename.c_str(), serie, mooring, unit,gnu_nice, time_format, verbose);
      break;

    case FORMAT_CODE_LIST:
      filename = rootname+ ".list";
      mgr_list_save_serie(filename, serie, mooring);
      break;

    case FORMAT_CODE_LIST2:
      filename = rootname+ ".list";
      mgr_list2_save_serie(filename, serie);
      break;

    case FORMAT_CODE_FOREMAN:
      filename = rootname+ ".dfo";
      status=mgr_saveforeman(filename.c_str(), serie,unit);
      break;

    case FORMAT_CODE_GLOSS:
    case FORMAT_CODE_HAWAII:
      filename = rootname+ ".dat";
      status=mgr_savehawaii(filename.c_str(),mooring.lat,mooring.lon,serie.x[0],serie.t,mooring.code,mooring.name, serie.first, serie.n, serie.t[serie.n-1]+0.5);
      break;

    case FORMAT_CODE_LEGEND:
      filename = rootname+ ".dat";
//      status=mgr_savehawaii(filename.c_str(),mooring.lat,mooring.lon,serie.x[0],serie.t,mooring.code,mooring.name, serie.first, serie.n, serie.t[serie.n-1]+0.5);
        status=mgr_Timeserie2Legend(filename.c_str(), mooring.lat,mooring.lon, serie.t, serie.n, serie.x, serie.nparam);
      break;
      
    case FORMAT_CODE_NETCDF:
      filename = rootname+ ".nc";
      status=mgr_save_timeserie_NetCDF(filename.c_str(),&serie,1);
      break;
    default:
      TRAP_ERR_EXIT(-1, "output format unknown");
    }

  printf ("written %s\n", filename.c_str());
  return(status);
}
