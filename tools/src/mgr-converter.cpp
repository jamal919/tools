
/*******************************************************************************
 
  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

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
\brief Convert between different sealevel time serie formats.
*/
/*----------------------------------------------------------------------------*/

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

#include "tools-structures.h"

#include "tides.h"
#include "fe.h"
#include "archive.h"
#include "poc-netcdf.hpp"
#include "poc-time.h"
#include "mgr-converter.h"
#include "mgr.h"
#include "functions.h"
#include "spectrum.h"
#include "filter.h"
#include "topo.h"

using namespace std;

extern void RC(void);

class metafield_t {
private:
public:
  string variable;
  string format;
  string gridfile, maskfile, datafile;
  string type;                           // structured/unstructured           
  string frame;
  float tag;
  double factor;
  
  metafield_t() {
    factor=1.0;
  }
    
  void destroy() {
    variable.clear();
    format.clear();
    gridfile.clear();
    maskfile.clear();
    datafile.clear();
    type.clear();
    frame.clear();
   }
};

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int parse_FieldOptions(const string & options, metafield_t * meta, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
  parse field input setting
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  vector<string> tokens,keys,values;
  int k;
  string delimiter=" ";
      
/*------------------------------------------------------------------------------
  first separate space-separated tokens */
  delimiter=" ";
  tokens=string_split(options, delimiter);
  
  if(debug) printf("tokens founds: %d\n",tokens.size());
  
  if(tokens.size()==1) {
    meta->datafile=options;
    return(0);
    }
  
/*------------------------------------------------------------------------------
  then identify key=value paires */
  delimiter="=";
  for(k=0;k<tokens.size();k++) {
    vector<string> tmp=string_split(tokens[k], delimiter);
    if(debug) printf("token %d, paire sizes: %d\n",k,tmp.size());
    keys.push_back(tmp[0]);
    if(tmp.size()==1) {
/*------------------------------------------------------------------------------
      old-fashioned input, no more supported */
      printf("wrong input setting %s; please set <file=%s format=NETCDF variable=xxx etc...\n",options.c_str(),options.c_str());
      TRAP_ERR_EXIT(-1, "parsing of %s failed\n", options.c_str());
      }
    values.push_back(tmp[1]);
    }
  
  for(k=0;k<keys.size();k++) {
    const string *keyk=&keys[k];
    char *valuek=poc_strdup(values[k].c_str());
    
    if(*keyk=="file") {
      meta->datafile=valuek;
      delete[]valuek;
      continue;
      }
    if(*keyk=="grid") {
      meta->gridfile=valuek;
      delete[]valuek;
      continue;
      }
    if(*keyk=="mask") {
      meta->maskfile=valuek;
      delete[]valuek;
      continue;
      }
    if(*keyk=="format") {
      meta->format=valuek;
      delete[]valuek;
      continue;
      }
    if(*keyk=="variable") {
      meta->variable=valuek;
      delete[]valuek;
      continue;
      }
    if(*keyk=="type") {
      meta->type=valuek;
      delete[]valuek;
      continue;
      }
    if(*keyk=="frame") {
      meta->frame=valuek;
      delete[]valuek;
      continue;
      }
    if(*keyk=="factor") {
      meta->factor=atof(valuek);
      delete[]valuek;
      continue;
      }
    
    STDERR_BASE_LINE_FUNC("key \""+*keyk+"\" not recognised\n");
    delete[]valuek;
    }
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<string> string_split2(const string & line, const  string  & separator)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nitems;
  size_t pointer, next;
  vector<string> array;
  string substring;
  string tmp;
  
  substring=line;
  
/*------------------------------------------------------------------------------
  remove initial occurence of delimiter */
  pointer = substring.find(separator);
  while(pointer==0) {
    tmp=substring.substr(pointer+separator.length());
    substring=tmp;
    pointer=substring.find(separator);
    }
    
  while(true) {
    pointer = substring.find(separator);
/*------------------------------------------------------------------------------
    remove repetitive occurence of delimiter */
    while(pointer==0) {
      tmp=substring.substr(pointer+separator.length());
      substring=tmp;
      pointer=substring.find(separator);
      }
    if(substring=="") break;
/*------------------------------------------------------------------------------
    at this point, substring first character is valid (i.e. NOT a delimiter) */
    next=substring.find(separator);
    if(next == string::npos) {
      string *token=new string;
      *token=substring.substr(0);
      array.push_back(*token);
      }
    else {
      string *token=new string;
      *token=substring.substr(0, next);
      array.push_back(*token);
      }
    if(next == string::npos)  break;
    tmp=substring.substr(next);
    substring=tmp;
    }

  return (array);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_SlpitLine(const string & line, const string & time_format, string & delimiter, bool disciplined, string & time_string, string & data_string)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  scan the most generic and standard ascii line of tide-gauge time series
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  int verbose=1;
  size_t pos;
  string tmp;
  
  if(disciplined) {
/*-----------------------------------------------------------------------------
    assume all time and data column separated with a unique delimiter */
    tmp=line;
    pos=tmp.find(delimiter);
    if(pos == string::npos) return(-1);
    }
  else {
    if(time_format=="CNES") {
      pos=line.find_first_not_of(" ");
      tmp=line.substr(pos);
      pos=tmp.find_first_of(" ");
      if(pos == string::npos) return(-1);
      }
    else {
/*-----------------------------------------------------------------------------
      ok if fix time format (i.e. calendar) */
      pos=line.find_first_not_of(" ");
      if(pos == string::npos) return(-1);
      tmp=line.substr(pos);
      pos=time_format.length();
      if(pos>tmp.length()) return(-1);
      } 
    }
   
  time_string=tmp.substr(0,pos);
  data_string=tmp.substr(pos);
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int scan_t(const string & time_string, const string & time_format, double & t, int recordI=-1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  scan the most generic and standard ascii line of tide-gauge time series
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int cnt, status;
  int year, month, day, hour, minute;
  double second, value;
  bool echo=true;
  
  string c_format;
  
  if(time_format=="YYYY-MM-DD HH-MM-SS") {
    c_format="%4d%*[-/ ]%2d%*[-/ ]%2d %2d%*[: ]%2d%*[: ]%lf";
    cnt=sscanf(time_string.c_str(), c_format.c_str(), &year, &month, &day, &hour, &minute, &second);
    if(cnt!=6) return(-1);
    t= julian_day(month, day, year)  - CNES0jd + hour /24. + minute / (24. * 60.) + second / (24.*3600.);
    }
  else if(time_format=="YYYY-MM-DD HH-MM") {
    c_format="%4d%*[-/ ]%2d%*[-/ ]%2d %2d%*[: ]%2d";
    cnt=sscanf(time_string.c_str(), c_format.c_str(), &year, &month, &day, &hour, &minute);
    second=0;
    if(cnt!=5) return(-1);
    t= julian_day(month, day, year)  - CNES0jd + hour /24. + minute / (24. * 60.) + second / (24.*3600.);
    }
  else if(time_format=="DD-MM-YYYY HH-MM-SS") {
    c_format="%2d%*[-/ ]%2d%*[-/ ]%4d %2d%*[: ]%2d%*[: ]%lf";
    cnt=sscanf(time_string.c_str(), c_format.c_str(), &day, &month, &year, &hour, &minute, &second);
    if(cnt!=6) return(-1);
    t= julian_day(month, day, year)  - CNES0jd + hour /24. + minute / (24. * 60.) + second / (24.*3600.);
    }
  else if(time_format=="DDMMYYYYHHMM") {
    c_format="%2d%2d%4d%2d%2d";
    cnt=sscanf(time_string.c_str(), c_format.c_str(), &day, &month, &year, &hour, &minute);
    if(cnt!=5) return(-1);
    second=0.0;
    t= julian_day(month, day, year)  - CNES0jd + hour /24. + minute / (24. * 60.) + second / (24.*3600.);
    }
  else if(time_format=="YYYYMMDDHHMM") {
    c_format="%4d%2d%2d%2d%2d";
    cnt=sscanf(time_string.c_str(), c_format.c_str(), &year, &month, &day, &hour, &minute);
    if(cnt!=5) return(-1);
    second=0.0;
    t= julian_day(month, day, year)  - CNES0jd + hour /24. + minute / (24. * 60.) + second / (24.*3600.);
    }
  else if(time_format=="CNES") {
    c_format="%lf";
    cnt=sscanf(time_string.c_str(), c_format.c_str(), &t);
    if(cnt!=1) return(-1);
    }

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int scan_z(const string & data, const string & delimiter, bool disciplined, vector<double> & z, double mask, int recordI=-1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  scan the most generic and standard ascii line of tide-gauge time series
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int cnt,status;
  int year, month, day, hour, minute;
  double second, value;
  bool echo=false;
    
  vector<string> tokens;
  
  z.clear();
  
  if(delimiter!="" and disciplined) {
    tokens=string_split(data, delimiter);
    for(int k=0; k<tokens.size()-1; k++) {
      cnt=sscanf(tokens[k+1].c_str(), "%lf", &value);
      if(cnt!=1) value=mask;
      z.push_back(value);
      if(echo) printf(" %lf\n",value);
      }
//     printf("%s : found data columns with delimiter <%s> in <%s>, discplined=%d \n", __func__, delimiter.c_str(), data.c_str(), tokens.size()-1);
    }
  else {
    tokens=string_split2(data, " ");
    for(int k=0; k<tokens.size(); k++) {
      cnt=sscanf(tokens[k].c_str(), "%lf", &value);
      if(cnt!=1) value=mask;
      z.push_back(value);
      if(echo) printf(" %lf\n",value);
      }
//     printf("%s : found data columns with delimiter <%s> in <%s>, undiscplined=%d \n", __func__, delimiter.c_str(), data.c_str(), tokens.size());
    }
        
  if(echo) printf(" \n");
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_ParseLineFormat(string & data, string & delimiter, bool disciplined, bool guess, int & ncolumns)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  scan the most generic and standard ascii line of tide-gauge time series
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int cnt;
  int year, month, day, hour, minute;
  double second, value;
  vector<string> tokens, tokens2;
  int verbose=1;
  string ascii_match;
  size_t pos;
  
  if(delimiter!="" and disciplined) {
    tokens=string_split(data, delimiter);
    ncolumns=tokens.size()-1;
    printf("%s : found data columns with delimiter <%s> in <%s>, discplined=%d \n", __func__, delimiter.c_str(), data.c_str(), tokens.size()-1);
    }
  else {
    tokens=string_split2(data, " ");
    ncolumns=tokens.size();
    printf("%s : found data columns with delimiter <%s> in <%s>, undiscplined=%d \n", __func__, delimiter.c_str(), data.c_str(), tokens.size());
    }

#if 0
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  <space> delimiter special case :
  
  - could be undisciplined ascci
  - could be true separator (hence handling missing data)
  
  following lines disabled, user MUST specify if disciplined or indisciplined
  format when delimiter is <space>
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(delimiter!="") tokens=string_split(data, delimiter);
  
  if(delimiter==" ") {
    string test;
    tokens2=string_split2(data, delimiter);
    int l1=tokens.size();
    int l2=tokens2.size();
    int missing=0;
    for(int k=0;k<tokens.size(); k++) {
      if(tokens[k].length()==0) missing++;
      }
    if(l1>l2) {
      printf("%s : found data columns with delimiter <%s> in <%s>, discplined=%d (%d missing) undisciplined=%d \n", __func__, delimiter.c_str(), data.c_str(), l1, missing, l2);
      ascii_match=test;
      }
    if(l1-missing==l2) {
      test="  ";
      tokens=string_split(data, test);
      }
    l1=tokens.size();
    missing=0;
    for(int k=0;k<tokens.size(); k++) {
      if(tokens[k].length()==0) {
        missing++;
        }
      }
    if(l1>l2) {
      printf("%s : found data columns with delimiter <%s> in <%s>, discplined=%d (%d missing) undisciplined=%d \n", __func__, test.c_str(), data.c_str(), l1, missing, l2);
      ascii_match=test;
      }
    if(l1-missing==l2) {
      test="   ";
      tokens=string_split(data, test);
      }
    l1=tokens.size();
    missing=0;
    for(int k=0;k<tokens.size(); k++) {
      if(tokens[k].length()==0) {
        missing++;
        }
      }
    if(l1>l2) {
      printf("%s : found data columns with delimiter <%s> in <%s>, discplined=%d (%d missing) undisciplined=%d \n", __func__, test.c_str(), data.c_str(), l1, missing, l2);
      ascii_match=test;
      }
    tokens=string_split(data, ascii_match);
    l1=tokens.size();
    missing=0;
    for(int k=0;k<tokens.size(); k++) {
      if(tokens[k].length()==0) {
        missing++;
        }
      }
    printf("%s : found data columns with delimiter <%s> in <%s>, discplined=%d (%d missing) undisciplined=%d \n", __func__, ascii_match.c_str(), data.c_str(), l1, missing, l2);
    }
  
  if(tokens.size()<2) {
    ncolumns=tokens.size()-1;
    if(verbose) printf("%s : found %d data columns with delimiter <%s> in %s, try other delimiters\n", __func__, ncolumns, delimiter.c_str(), line.c_str());
    delimiter=="";
    }

  if(delimiter=="" or tokens.size()<2) {
#endif

  if(guess and disciplined and tokens.size()<2) {
    while(true) {
/*------------------------------------------------------------------------------
      typical csv delimiter */
      delimiter=";";
      tokens=string_split(data, delimiter);
      if(tokens.size() !=0) break;
/*------------------------------------------------------------------------------
      typical csv delimiter */
      delimiter="\t";
      tokens=string_split(data, delimiter);
      if(tokens.size() !=0) break;
/*------------------------------------------------------------------------------
      typical ascii delimiter */
      delimiter="  ";
      tokens=string_split(data, delimiter);
      if(tokens.size() !=0) break;
/*------------------------------------------------------------------------------
      typical ascii delimiter */
      delimiter=" ";
      tokens=string_split(data, delimiter);
      if(tokens.size() !=0) break;
      break;
      }
    ncolumns=tokens.size()-1;
    }
  
  
  if(verbose) printf("%s : found %d data columns with delimiter <%s> in <%s> \n", __func__, ncolumns, delimiter.c_str(), data.c_str());
    
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_ParseOpenFormat(const char *filename, const string & time_format, string & delimiter, bool disciplined, bool guess, int & nrecords, int & ncolumns)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  2010/01/01 00:00:00  1.68
  2010/01/01 00:05:00  1.74

  or

  2006-01-01 00:00:00   1.020
  2006-01-01 00:05:00   0.970
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i,status=-1;
  vector<double> values;
  date_t date;
  double mask=-9999.9, t;
  int count;
  
  std::string line="",time_string, data_string;

  std::ifstream input(filename);
  if(input.fail())
    TRAP_ERR_RETURN(0,1,"%s: ifstream(\"%s\") error (%d %s)\n",__func__,filename,errno,strerror(errno));
  
  nrecords=ncolumns=0;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  count data lines (timeserie length)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nrecords=mgr_line_count((string) filename);

// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//   
//   parse time format
//   
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// 
//   string c_format;
//   
//   if(time_format=="YYYY-MM-DD HH-MM-SS") {
//     c_format="%4d%*[-/ ]%2d%*[-/ ]%2d %2d%*[: ]%2d%*[: ]%lf";
//     }
//   else if(time_format=="DD-MM-YYYY HH-MM-SS") {
//     c_format="%2d%*[-/ ]%2d%*[-/ ]%4d %2d%*[: ]%2d%*[: ]%lf";
//     }
//   else if(time_format=="CNES") {
//     c_format="%lf";
//     }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  count data columns
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  count=-1;
  while (count<nrecords) {
    std::getline( input, line );
    count++;
    if(line.c_str()[0]=='#') continue;
    status=mgr_SlpitLine(line, time_format, delimiter, disciplined, time_string, data_string);
    status=scan_t(time_string, time_format, t);
    if(status!=0) {
      printf("rejected line : <%s> \n", line.c_str());
      continue;
      }
    status=mgr_ParseLineFormat(data_string, delimiter, disciplined, guess, ncolumns);
    if(status!=0) continue;
    break;
    }
  
  input.close();
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_loadOpenFormat(const char *filename, tseries_t *serie, const string & time_format, string & delimiter, bool disciplined, bool guess, bool no_redundancy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  2010/01/01 00:00:00  1.68
  2010/01/01 00:05:00  1.74

  or

  2006-01-01 00:00:00   1.020
  2006-01-01 00:05:00   0.970
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  vector<double> values;
  double t;
  date_t date;
  double mask=-9999.9;
  
  std::string line="",time_string, data_string;
  cout << endl << "-------- starting conversion --------" << endl << endl;

  std::ifstream input(filename);
  if(input.fail())
    TRAP_ERR_RETURN(0,1,"%s: ifstream(\"%s\") error (%d %s)\n",__func__,filename,errno,strerror(errno));

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  parse datafile content
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  int nrecords, ncolumns;
  status=mgr_ParseOpenFormat(filename, time_format, delimiter, disciplined, guess, nrecords, ncolumns);

  serie->nparam=ncolumns;
  serie->x=new double*[ncolumns];
  for(int k=0;k<ncolumns;k++) serie->x[k]=new double[nrecords];
  serie->t=new double[nrecords];
  serie->n=nrecords;
  serie->mask=99999;

  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  read datafile content
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  int count = 0;
  
  for(int i = 0; i < nrecords; i++){
    std::getline( input, line );
    if(line.c_str()[0]=='#') continue;
    status=mgr_SlpitLine(line, time_format, delimiter, disciplined, time_string, data_string);
    status=scan_t(time_string, time_format, t, i);
    if(status!=0) {
      printf("rejected line : <%s> \n", line.c_str());
      i--;
      continue;
      }
    status=scan_z(data_string, delimiter, disciplined, values, serie->mask, i);
    if(status!=0) {
      printf("rejected line : <%s> \n", line.c_str());
      i--;
      continue;
      }
    if(no_redundancy and count>0) {
      if(t<=serie->t[count-1]) continue;
      }
    for(int k=0;k<values.size();k++) serie->x[k][count]=values[k];
    serie->t[count] = t;
    count++;
    }

  serie->n=count;

  serie->first=poctime_getdatecnes(serie->t[0],'d');

  return(nrecords);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_converter_local(const char *file_in, const string format_in, string delimiter, string time_format, int ncolumns_GNU, bool split,
    const char *header_format, const char *mooring_info,
    date_t date_start, date_t date_final, date_t step, int reconstruct, double time_shift, double units_correction, double local_datum_correction, string datum_correction,
    const char *file_out_prefix, string suffix, const char *format_out, tseries_t ***out_series,
    double resampling, bool gnu_nice)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  tseries_t *serie=NULL;
  int nseries, fmt;
  mooring_t mooring;
  date_t real_start, real_final;
  date_t start, final, next;
  int status;
  int j,k,n;
  int *nreduced;
  char output_dated[1024];
  char *output=new char[1024], rootname[1024];
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

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  initialisation
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  mooring.lon=NAN;
  mooring.lat=NAN;
    
  status= mgr_init_formats();
  
  if (file_out_prefix!=0) {
    sprintf(rootname, "%s", file_out_prefix);
    }
  else {
    cpos = strrchr(file_in, '/');
    if (!cpos) cpos = file_in;
    else cpos++;
    sprintf(rootname, "%s", cpos);
    if ((pos = strrchr(rootname,'.'))) *pos = '\0';
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  load time series 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  fmt=mgr_ParseFormat(file_in, format_in.c_str());

  switch(fmt) {
      
    case FORMAT_CODE_DMY:
    case FORMAT_CODE_YMD:
    case FORMAT_CODE_GNU:
      {
      bool discplined=(delimiter!="");
      bool no_redundancy=true;
      tseries_t *tmp=new tseries_t;
      bool guess=true;
      status=mgr_loadOpenFormat(file_in, tmp, time_format, delimiter, discplined, guess, no_redundancy);
      if(split) {
        tmp->split(serie);
        nseries=tmp->nparam;
        tmp->destroy();
        }
      else {
        nseries=1;
        serie=tmp;
        }
      }
      break;

    default:
      nseries = mgr_load_timeserie(file_in, &serie, 'm', format_in.c_str(), ncolumns_GNU, header_format, &mooring);
    }
  
  if(nseries<=0)
    return nseries;
  
  if(mooring_info) mooring=mgr_decode_mooring_info(mooring_info);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  apply time shift; origin setting must be controlled
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(k=0; k<nseries; k++) {
    double offset;
    offset=time_shift/24./3600.+cnes_time(serie[k].origin,'d');
    if(offset==0.) continue;
    for(int l=0;l<serie[k].n;l++) {
      serie[k].t[l]+=offset;
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  apply units factor
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(k=0; k<nseries; k++) {
    if(units_correction==1.) continue;
    for(int j=0;j<serie[k].nparam;j++) {
      serie[k].stat(j);
      for(int l=0;l<serie[k].n;l++) {
        if(serie[k].x[j][l]!=serie[k].mask) serie[k].x[j][l]*=units_correction;
        }
      serie[k].stat(j);
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  apply datum shift
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(k=0; k<nseries; k++) {
    if(local_datum_correction==0.) continue;
    for(int j=0;j<serie[k].nparam;j++) {
      serie[k].stat(j);
      for(int l=0;l<serie[k].n;l++) {
        if(serie[k].x[j][l]!=serie[k].mask) serie[k].x[j][l]+=local_datum_correction;
        }
      serie[k].stat(j);
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  apply geodetic shift
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  string topofile;
  grid_t topogrid;
  float *topo, mask, value;
  bool debug;
  
  if(datum_correction!="") {
    metafield_t meta;
    status=parse_FieldOptions(datum_correction, &meta, debug);

    status=topo_loadfield(meta.datafile.c_str(), meta.variable.c_str(), &topogrid, &topo, &mask, debug);
  
    for(k=0; k<nseries; k++) {
      double x,y;
      x=map_recale(topogrid,mooring.lon);
      y=mooring.lat;
      value=NAN;
      status=map_interpolation(topogrid, topo, mask, x,y,&value);
      if(status!=0) continue;
      value*=meta.factor;
      printf("serie %d : correction=%lf\n",value);
      for(int j=0;j<serie[k].nparam;j++) {
        serie[k].stat(j);
        for(int l=0;l<serie[k].n;l++) {
          if(serie[k].x[j][l]!=serie[k].mask) serie[k].x[j][l]+=value;
          }
        serie[k].stat(j);
        }
      }
    }

/*------------------------------------------------------------------------------
  split series */
  reduced  = new tseries_t*[nseries];
  nreduced = new int[nseries];

/*------------------------------------------------------------------------------
  build k output series */
  int nb_tot_series=0;
  
  for(k=0; k<nseries; k++) {

    if(mooring_info) serie[k].mooring=mooring;
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
    construct output filename
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
    start build output file names */
    if (file_out_prefix!=0 and (nseries != 1) )
      sprintf(output, "%s.%4.4d", rootname, k);
    else
      sprintf(output, "%s", rootname);
    
/*------------------------------------------------------------------------------
    add mooring name to file out */
    if (serie[k].mooring.initialized and serie[k].mooring.name!=0) {
      if(strstr(output, serie[k].mooring.name)==0) {
        sprintf(mooring_name, "%s", serie[k].mooring.name);
        char *p=mooring_name;
        while ((p=strchr(p,int(' '))))
          p[0] = '_';
        sprintf(output, "%s_%s", output, mooring_name);
        }
      }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
    set splitting limits
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
    split or clip time series : get start / end dates of the serie */
    real_start = poctime_getdatecnes(serie[k].t[0], 'd');
    real_final = poctime_getdatecnes(serie[k].t[serie[k].n-1],'d');
    
/*------------------------------------------------------------------------------
    select starting date */
    if (date_start.year!=0) {
      if (poctime_is_after(date_start, real_final)) {
        STDERR_BASE_LINE("analysis period starts (%s) after time series finishes (%s), no data for %s \n", sgetdate(date_start), sgetdate(real_final), file_in);
        return 0;
        }
      real_start = poctime_latest(real_start, date_start);
      }
/*------------------------------------------------------------------------------
    if no starting date but steps: then start at begining of the year ... */
    if (step.year!=0) {
      real_start.month = 1;
      real_start.day = 1;
      /* always same years spliting */
      if (step.year!=1) {
        if (date_start.year!=0) {
          int year_modulo;
          year_modulo = floor((real_start.year - date_start.year) / step.year);
          real_start.year = date_start.year + year_modulo *  step.year;
          }
        }
      }
/*------------------------------------------------------------------------------
    ... or of the month */
    else if (step.month!=0) {
      real_start.day = 1;
      }

/*------------------------------------------------------------------------------
    select ending date */
    if (date_final.year!=0) {
      if (poctime_is_before(date_final, real_start)) {
        STDERR_BASE_LINE("analysis period finishes (%s) after time series starts (%s), no data for %s \n", sgetdate(date_final), sgetdate(real_start), file_in);
        return 0;
        }
      real_final = poctime_soonest(real_final, date_final);
      }
/*------------------------------------------------------------------------------
    no ending date: go until end of year ... */
    if (step.year!=0) {
      real_final.day = 31;
      real_final.month = 12;
      }
/*------------------------------------------------------------------------------
    ... or end of month */
    else if (step.month!=0) {
      real_final.day = 31;
      }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
    construct splitted consecutive timeseries sequence
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
    count number of split series */
    if (step.year or step.month)
      nreduced[k] = poctime_count_steps(real_start, real_final, step);
    else {
      nreduced[k] = 1;
//       step.year = 100;
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
      if(step.year!=0 or step.month!=0) next = poctime_add_step(start, step, real_final, &final);
      else final=real_final;
      printf("process timeseries : start %04d/%02d/%02d %fh final %04d/%02d/%02d %fh\n",
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
        string filename=output_dated;
        if(suffix!="") filename+="-"+suffix;
        if(reconstruct==1)
          status=mgr_reconstruct(reducedkj, &mooring);
        status=mgr_save_timeserie(filename, reducedkj->mooring, *reducedkj, 'm', format_str, gnu_nice);
        }
      
/*------------------------------------------------------------------------------
      next sequence start */
      start = next;

      nb_tot_series += nreduced[k];
      }

    serie[k].destroy();
    } /* loop over k series */

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  return series
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

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
    "NAME AND VERSION:\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s <file> [options]\n", prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Convert between different sealevel time serie formats.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  show this help and exit\n"
    +isFormat_help("   ",". Default: GNU",". Default: GNU")+
    "   -n <nb>                  : nb of columns for input GNU format (default=1)\n"
    "  -s  followed by the start date. See DATE FORMATS below\n"
    "  -f  followed by the end date. See DATE FORMATS below\n"
    "   --step=<n_years>y        : split time series per nb years or nb month\n"
    "   --reconstruct            : enable reconstruction\n"
    "   --resampling             : followed by resampling period in seconds. If not given, no resampling will be done.\n"
    "   --time-correction        : followed by the offset to correct by, in seconds\n"
    "   --mooring_info=<info>    : input info on mooring, eg: \"lon=<lon> lat=<lat> code=<code> name=<name> depth=<depth>\"\n"
    "   --header_format=<format> : format to decode mooring info on file\n"
    "   -o <name>                : to save output in files with prefix <name>\n"
    "   --gnu_nice               : nice gnu (allways set to True)\n"
    "\n"
    "TIP\n"
    "  In `ipython --pylab=...', look at the reconstruction with:\n"
    "d=loadtxt(...);plot_date(datestr2num('1950-01-01')+d[:,0],d[:,1:],'-');legend(('?','input','filtered','prediction','reconstructed'))\n"
    );
  print_poctime_scan_date_help(0);
  mgr_print_formats();
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  convert (usually tide gauges) timeseries file from one format to another

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int ncol=1;
  int nprocs=omp_get_num_procs(); 
  date_t step(0,0,0,0.);
  int n,status;
  const char *keyword,*s;
  const char *input=0, *output=0, *mooring_info=0, *header_format=0;
  string rootname;
  bool gnu_nice=true;
  date_t start(0,0,0,0.), final(0,0,0,0.);
  int reconstruct=0;
  double time_correction=0, local_datum_correction=0.0, units_correction=1.0, resampling=NAN;
  string input_format;
  string output_format="GNU";
  string suffix;
  string time_format="YYYY-MM-DD HH-MM-SS";
//   string delimiter=" ";
  string delimiter="";
  string datum_correction;
  bool split=true;

  fct_echo( argc, argv);
//   int nprocs=omp_get_num_procs();
  printf("nmprocs=%d\n",nprocs);
  
//   RC();

  if (argc==1){
    print_help(argv[0]);
    exit(0);
    }

  n=1;
  while (n < argc) {
    keyword=argv[n];
    
    if(isFormat(argv,&n,&input_format,&output_format)){
      continue;
      }
    
    switch (keyword[0]) {
      case '-':
        if( strcmp(keyword,"-h")==0 or
            strcmp(keyword,"--help")==0 ){
          print_help(argv[0]);
          exit(0);
          }

        switch (keyword[1]) {

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

        case 'o' :
          output=argv[n+1];
          n++;
          n++;
          break;

        case 'n' :
          sscanf(argv[n+1],"%d",&ncol);
          n++;
          n++;
          break;

        case '-' :
          if(strcmp(keyword,"--suffix")==0 ) {
            suffix=argv[n+1];
            n++;
            n++;
            break;
            }
          if(strcmp(keyword,"--nice=yes")==0 || strcmp(keyword,"--gnu_nice")==0) {
            gnu_nice=1;
            n++;
            break;
            }
          if(strcmp(keyword,"--nice=no")==0) {
            gnu_nice=0;
            n++;
            break;
            }
          if(strcmp(keyword,"--reconstruct")==0 ) {
            reconstruct=1;
            printf ("series reconstruction on\n");
            n++;
            break;
            }
          if(strcmp(keyword,"--resampling")==0 ) {
            sscanf(argv[n+1],"%lf",&resampling);
            resampling/=86400.;
            n++;
            n++;
            break;
            }
          if(strcmp(keyword,"--split=yes")==0 ) {
            split=true;
            n++;
            break;
            }
          if(strcmp(keyword,"--split=no")==0 ) {
            split=false;
            n++;
            break;
            }
          if(strcmp(keyword,"--delimiter")==0 ) {
            delimiter=argv[n+1];
            n++;
            n++;
            break;
            }
          if(strcmp(keyword,"--time-format")==0 ) {
            time_format=argv[n+1];
            n++;
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
          if(strcmp(keyword,"--local-datum-correction")==0 ) {
            sscanf(argv[n+1],"%lf",&local_datum_correction);
            printf ("datum correction (additive, applied after units correction if any)=%lf m\n",local_datum_correction);
            n++;
            n++;
            break;
            }
          if(strcmp(keyword,"--datum-correction")==0 ) {
            datum_correction=argv[n+1];
//             printf ("datum correction (additive, applied after units correction if any)=%lf m\n",datum_correction);
            n++;
            n++;
            break;
            }
          if(strcmp(keyword,"--units-correction")==0 ) {
            sscanf(argv[n+1],"%lf",&units_correction);
            printf ("units correction (multiplicative)=%lf m\n",units_correction);
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
              STDOUT_BASE_LINE("bad syntax %s for option --step\n",argv[n]);
              print_help(argv[0]);
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
          input=argv[n];
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

/*-----------------------------------------------------------------------------
  force output for interactive mode */
  if(output==0) output="";
  rootname = mgr_build_out_prefix(input, output, NULL, NULL);

/*-----------------------------------------------------------------------------
  let's convert ! */
  int nb_series;
  nb_series = mgr_converter_local(input, input_format, delimiter, time_format, ncol, split, header_format, mooring_info,
                            start,  final, step, reconstruct, time_correction, units_correction, local_datum_correction, datum_correction,
                            &(rootname[0]), suffix, output_format.c_str(), NULL, resampling, gnu_nice);

//  filename=rootname+".dfo";
//  status=mgr_saveforeman((const char *)filename.c_str(), serie,'m');

// Verdon   45 32   1 02
//   mooring.lat=45.53333;
//   mooring.lon=-1.03333;

// Pauillac 45 12   0 45
//   mooring.lat=45.20;
//   mooring.lon=-0.75;

//  filename=rootname+".dat";
//  mgr_savehawaii(filename.c_str(),mooring.lat,mooring.lon,serie.x,serie.t,mooring.code,mooring.name, serie.first, serie.n, serie.t[serie.n-1]+0.5);

  STDOUT_BASE_LINE(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(0);
}

// example:
// r /gsa/fles/gloss/hourly-20120625/obs/origin/h418.dat -s 2010 -f 2010 -o /tmp/ -iF GLOSS --step=1y
