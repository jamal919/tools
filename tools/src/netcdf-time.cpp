
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief poc-netcdf decoding definitions
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <string.h>

#include "tools-structures.h"
#include "poc-netcdf.def"
#include "netcdf-proto.h"
#include "poc-time.h"
#include "fe.h"
#include "zapper.h"
#include "poc-grib.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_writetime(const char *filename, int frame, int id, double time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
   int  ncid;			/* netCDF id */
   size_t index[3];
   int status=0;

   index[0]=frame;

   status=nc_open(filename,NC_WRITE,&ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status=nc_put_var1_double(ncid,id,index,&time);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status = nc_close(ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
 error:
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_writetime(const char *filename, int frame, const char *name, double time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
   int  ncid,id;			/* netCDF id */
   size_t index[3];
   int status=0;

   index[0]=frame;

   status=nc_open(filename,NC_WRITE,&ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   
        status=nc_inq_varid(ncid,name,&id);
        if(status!=0) id=-1;
        
   status=nc_put_var1_double(ncid,id,index,&time);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status = nc_close(ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
 error:
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_settime(const char *filename, int id, double *time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Fills the time buffer of a NetCDF file
/** that has been initialised to have a time variable,
e.g. with poc_createfile() and create_ncvariable().

\date 2011-09-26 Damien Allain : review and documentation

\param *filename NetCDF file path
\param id time variable id
\param *time time values
*/
/*----------------------------------------------------------------------------*/
{
   int  ncid;
   int status;

   status=nc_open(filename,NC_WRITE,&ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   
   status=nc_put_var_double(ncid,id,time);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   
   status = nc_close(ncid);
   nc_check_error(status,__LINE__,__FILE__);
   
 error:
  return(status);
}


/*----------------------------------------------------------------------------*/
/** Get times of frames of a given NetCDF file
with the origin from the "time_origin" attribute.
It gets the units from the "units" or "time_origin" attribute.
\param ncid
\param tvid time variable id
\param origine pointer to the time origin. If there is neither a "time_origin" nor a "units" attribute, default to CNES origin of 1950/01/01 00:00. If neither the "time_origin" nor the "units" attribute give the date in a clear manner, default to NADate.
\param time pointer to array of times of frames in seconds. Can be NULL.
\param nframes pointer to the number of frames

\todo 2011-09-29 Damien Allain : use info:/libc/General%20Time%20String%20Parsing getdate_r
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_gettime(int ncid, int tvid, date_t *origine, double **time, size_t *nframes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, k;
  size_t lengthhp;
  int tdid;//time var id and time dim id
  int day, month, year;
  int hour, minute, seconds;
  nc_type xtypep;
  char smonth[1024];
  char *text, *pointer, *time_units;
  double factor;
  year=month=day=0;
  hour=minute=seconds=0;//date may be specified, but not time, so default to 0

/*------------------------------------------------------------------------------
  time origin */
  status = nc_inq_att(ncid, tvid, "time_origin", &xtypep, &lengthhp);
  if(status != 0)
/*------------------------------------------------------------------------------
    no time origin attribute found, let's try units*/
    goto units;

  text = new char[lengthhp + 1];
  status = nc_get_att_text(ncid, tvid, "time_origin", text);
  if(status != NC_NOERR)
    goto error;

  text[lengthhp] = '\0';
  sscanf(strtok(text, "-"), "%d", &year);
  pointer = strtok(NULL, "-");
  if(pointer==0)
    goto error;
  sscanf(pointer, "%3s", &smonth);
  pointer = strtok(NULL, " ");
  if(pointer==0)
    goto error;
  sscanf(pointer, "%d", &day);

  pointer = strtok(NULL, ":");
  if(pointer != NULL)
    sscanf(pointer, "%d", &hour);
  pointer = strtok(NULL, ":");
  if(pointer != NULL)
    sscanf(pointer, "%d", &minute);
  pointer = strtok(NULL, "\0");
  if(pointer != NULL)
    sscanf(pointer, "%d", &seconds);

  for(k = 0; k < 12; k++)
    if(strcasecmp(smonth, uc_names[k]) == 0)
      break;
  month = k + 1;
/*------------------------------------------------------------------------------
  patch: month given as number and not text*/
  if(k == 12)
    sscanf(smonth, "%d", &month);
  if(month==13){
    status=NC_ECHAR;/* Attempt to convert between text & numbers */
    nc_check_error(status,__LINE__,__FILE__,"Can not interpret month name %s",smonth);
    return status;
    }

units:
  factor = 1.;
  status = nc_inq_att(ncid, tvid, "units", &xtypep, &lengthhp);
  if(status != NC_NOERR) {
/*------------------------------------------------------------------------------
    assume undocumented time to be Cnes time */
    year=1950;
    month=1;
    day=1;
    hour=minute=seconds=0;
    factor = 24 * 3600.;
    goto load;
    }
  
  text = new char[lengthhp + 1];
  status = nc_get_att_text(ncid, tvid, "units", text);
  if(status != NC_NOERR)
    goto error;

  text[lengthhp] = '\0';
  time_units = strdup(strtok(text, " "));
  
  factor=1.;
  switch(tolower(time_units[0])){
    case 'd':
      factor*=24.;
    case 'h':
      factor*=60.;
    case 'm':
      factor*=60.;
    case 's':
      break;
    default:
      status=ENOSYS;/* Function not implemented */
      NC_CHKERR_LINE_FILE(status,"Unknown time unit \"%s\"",time_units);
      return status;
    }

  status = nc_get_att_text(ncid, tvid, "units", text);
  if(status != NC_NOERR)
    goto error;
  text[lengthhp] = '\0';
  pointer = strstr(text, "since");
  if(pointer!=0) {
    pointer += 6;
    sscanf(pointer, "%d%*c%d%*c%d%*c%d%*c%d%*c%d", &year, &month, &day, &hour, &minute, &seconds);
    }
  zaparr(text);

/*------------------------------------------------------------------------------
  assume time is 1D buffer */
/*
  jday = julian_day(month, day, year);
*/
load:

  if(year==0){//year 0 does not exit : year 1 follows year -1
    *origine=NADate;
    }
  else{
    origine->year=year;
    origine->month=month;
    origine->day=day;
    origine->second=hour*3600+minute*60+seconds;
    }

  //get NUMBER OF dimensions ...
  status=nc_inq_varndims(ncid,tvid,&tdid);
  if(status!=NC_NOERR){nc_check_error(status,__LINE__,__FILE__);return status;}
  if(tdid!=1)//time should only have 1 dimension
    return NC_EVARSIZE;
  //.. before getting the ID of the dimension ...
  status=nc_inq_vardimid(ncid,tvid,&tdid);
  if(status!=NC_NOERR){nc_check_error(status,__LINE__,__FILE__);return status;}
  //so that you can (at last) have the number of frames
  status=nc_inq_dimlen(ncid,tdid,nframes);
  if(status!=NC_NOERR){nc_check_error(status,__LINE__,__FILE__);return status;}
  if(time!=NULL){
    *time = new double[*nframes];
    status = nc_get_var_double(ncid, tvid, *time);
    if(status!=NC_NOERR){nc_check_error(status,__LINE__,__FILE__);return status;}
  /*       if(time[t-1]*factor/3600./24.>1.e+10) time[t-1]=0.; */
  /*       info->time=time[t-1]*factor/3600./24.+jday-julian_day(1,1,1950); */
  //  zaparr(time);
  /*------------------------------------------------------------------------
    return time in seconds */
    for(size_t k=0;k<*nframes;k++) {
      (*time)[k]*=factor;
      }
    }

  return (0);

error:
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_gettime(int ncid, date_t *origine, double **time, size_t *nframes, cdfvar_t *timevar)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Get times of frames of a given NetCDF file
/** assuming the time variable name is one of ::timeNames
See poc_gettime(int,int,date_t*,double**,size_t*) for more help, including the one about the other parameters.

\return NC_ENOTVAR if not time variable has been found
*/
/*----------------------------------------------------------------------------*/
{
  int status;
  int tvid=NC_ENOTVAR;//time var id
  int i,nvars;//variable index
  char name[NC_MAX_NAME+1];

  ///It check all the names of ::timeNames
  status=nc_inq_nvars(ncid,&nvars);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"nc_inq_nvars() error");
  
  for(i=0;i<nvars;i++){
    status=nc_inq_varname(ncid,i,name);
    if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"nc_inq_varname(,%d,) error",i);
    if(isT(name)){
      tvid=i;
      break;
      }
    }
  
  if(tvid<0){
    if(time!=0)
      *time=0;
    *nframes=1;
    NC_TRAP_ERROR(return,tvid,1,"no time variable found");
    }

  if(timevar!=NULL){
    cdf_varinfo(ncid,tvid,timevar,0);
    if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"error on variable %s",timeNames[i]);
    }

  return poc_gettime(ncid,tvid,origine,time,nframes);
}


/*----------------------------------------------------------------------------*/
/** Get times of frames of a given NetCDF file
\param filename path of the NetCDF file
See poc_gettime(int,date_t*,double**,size_t*,cdfvar_t*) for more help, including the one about the other parameters.
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_gettime(const char *filename, date_t *origine, double **time, size_t *nframes, cdfvar_t *timevar)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int ncid,status;

  if(is_grib(filename)){
#ifdef HAVE_LIBGRIB_API
    status=poc_grib_gettime(filename,origine,time,nframes);
    return status;
#else
    TRAP_ERR_RETURN(-1,1,"Compile with grip_api\n");
#endif
    }
  
  status = nc_open(filename, 0, &ncid);
  if(status != NC_NOERR)
    return status;

  status=poc_gettime(ncid, origine, time, nframes, timevar);
  NC_CHKERR_BASE_LINE(status,"poc_gettime((\"%s\"),...) error",filename);
  nc_close(ncid);

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_gettime_template(T file, date_t *origine, double **time, int *nframes, cdfvar_t *timevar=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** Get times of frames of a given NetCDF file
See poc_gettime(const char*,date_t*,double**,size_t*,cdfvar_t*) for more help.
Declared because of size_t<->int conflict
*/
/*----------------------------------------------------------------------------*/
{
  size_t snf;
  int status;
  
  status=poc_gettime(file, origine, time, &snf, timevar);
  *nframes=snf;
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_gettime(int file, date_t *origine, double **time, int *nframes, cdfvar_t *timevar)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_gettime_template(file,origine,time,nframes,timevar);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_gettime(char *file, date_t *origine, double **time, int *nframes, cdfvar_t *timevar)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_gettime_template(file,origine,time,nframes,timevar);
  return status;
}


/*----------------------------------------------------------------------------*/
/** Get times of frames of a given NetCDF file
\param ncid
\param variable variable on which poc_decode_associates() will be called, so that decoded_t::vt will be used as the time variable id
\param global also necessary for poc_decode_associates()
See poc_gettime(int,int,date_t*,double**,size_t*) for more help, including the one about the other parameters
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_gettime(int ncid, cdfvar_t variable, cdfgbl_t global, date_t *origine, double **time, size_t *nframes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int tvid;//time var id
  decoded_t decoded;
  bool debug=false;

  status=poc_decode_associates( variable, global, &decoded, 1, debug);
  if(status!=NC_NOERR)return status;

  tvid=decoded.vt;
  if(tvid == -1)
    return NC_ENOTVAR;
  
  return poc_gettime(ncid,tvid,origine,time,nframes);
}


/*----------------------------------------------------------------------------*/
/** Get times of frames of a given NetCDF file
\param filename path of the NetCDF file
See poc_gettime(int,cdfvar_t,cdfgbl_t,date_t*,double**,size_t*) for more help, including the one about the other parameters
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_gettime(const char *filename, cdfvar_t variable, cdfgbl_t global, date_t *origine, double **time, size_t *nframes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int ncid,status;
  
  status = nc_open(filename, 0, &ncid);
  if(status != NC_NOERR)
    return status;
  status=poc_gettime(ncid, variable, global, origine, time, nframes);
  nc_close(ncid);

  //if(status != NC_NOERR)
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_timeorigin(const char *filename, int vtime, date_t * reference)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int status, ncid, k;
  size_t lengthhp;
  int day, month, year;
  int hour, minute, seconds;
  nc_type xtypep;
  char smonth[1024];
  char *text=NULL, *pointer=NULL;

  status = nc_open(filename, 0, &ncid);
  if(status != NC_NOERR)
    goto error;

/*------------------------------------------------------------------------------
  time origin */
  status = nc_inq_att(ncid, vtime, "time_origin", &xtypep, &lengthhp);
  if(status != 0)
    goto error;

  text = new char[lengthhp + 1];
  status = nc_get_att_text(ncid, vtime, "time_origin", text);
  if(status != NC_NOERR)
    goto error;

  text[lengthhp] = '\0';
  sscanf(strtok(text, "-"), "%d", &year);
  pointer = strtok(NULL, "-");
  if(pointer==0)
    goto error;
  sscanf(pointer, "%3s", &smonth);
  pointer = strtok(NULL, "-");
  if(pointer==0)
    goto error;
  sscanf(pointer, "%d", &day);

  hour = 0;
  minute = 0;
  seconds = 0;
  pointer = strtok(NULL, ":");
  if(pointer != NULL)
    sscanf(pointer, "%d", &hour);
  pointer = strtok(NULL, ":");
  if(pointer != NULL)
    sscanf(pointer, "%d", &minute);
  pointer = strtok(NULL, "\0");
  if(pointer != NULL)
    sscanf(pointer, "%d", &seconds);

  for(k = 0; k < 12; k++)
    if(strcmp(smonth, uc_names[k]) == 0)
      break;
  month = k + 1;
/*------------------------------------------------------------------------------
  patch: month given as number and not text*/
  if(k == 12)
    sscanf(smonth, "%d", &month);

  reference->day = day;
  reference->month = month;
  reference->year = year;
  reference->second = seconds + 60. * minute + 3600. * hour;

  status = nc_close(ncid);
  if(status != NC_NOERR)
    goto error;

  return (0);

error:
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_timefilterlist(vector<string> *pathlist,double *start,double *final,double **ts,double origd)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Keep only the files between a start and a finish date in a given filelist
/** and updates the given start and finish date to the dates of the first and last frames within these times.
\param pathlist pointer to original list of file paths, modified to exclude files that have only frames outside of the given time range
\param start pointer to required start date, inclusive, in CNES seconds, modified to reflect the time of the first frame at or after the given start time
\param final pointer to required end date, exclusive, in CNES seconds, modified to reflect the time of the last frame before the given end time
\param **ts=NULL a pointer to the (sorted if the files also are) array of the times of the frames in CNES seconds
\param origd=NAN if specified, default time origin in CNES seconds if not in the files
\returns the number of frames kept
*/
/*----------------------------------------------------------------------------*/
{
  int i,pn=pathlist->size();   //file index, number of paths, time index
  int j,n;                     //time index,number of frames,
  int status;                  //NetCDF status
  double **fts,**ftsi,*ftsij;  /* times of frames of files */
  size_t *ns,tn=0,tfi=-1;         //numbers of frames and total thereof,total time index
  date_t orig;                 //time origin
  double t0;                   //time origin
  fts=new double*[pn];
  ns=new size_t[pn];
  struct timeval stv;          //start timeval (for progression)
  
  gettimeofday(&stv);
  
/*------------------------------------------------------------------------------
  For all files: */
  for(i=0;i<pn;i++){
    const string &path=(*pathlist)[i];
    if(i>0) printf("%s%s%s%s",cr,el,cuu1,el);
    printf("(%d/%d;%04.3g)frames of "+path+":\n",i+1,pn,difftime(stv));
    ///- it first gets the times of all the frames with poc_gettime()
    status=poc_gettime(path.c_str(),&orig,NULL,&ns[i]);
    if(status!=0){
      if(ts!=0)
        *ts=0;
      NC_CHKERR_BASE_LINE(status,"poc_gettime(\""+path+"\",,0,) error: assuming number of frames is the number of files, i.e. %d",pn);
      return pn;
      }
    printf("%d",ns[i]);fflush(stdout);
    
    ftsi=&fts[i];
    *ftsi=NULL;
    status=poc_gettime(path.c_str(),&orig,ftsi,&ns[i]);
    t0=cnes_time(orig,'s');
    if(isnan(t0))t0=origd;
    if(status!=NC_NOERR || isnan(t0)){
      fprintf(stderr,"\n");
      if(*ftsi!=NULL) deletep(ftsi);
      if(status!=NC_NOERR)
        NC_CHKERR_BASE_LINE(status,"error getting times from "+path);
      else if(isnan(t0))
        NC_CHKERR_BASE_LINE(NC_ENOTATT,"error getting time origin from "+path);
      continue;
      }
    printf(" in[%s ; %s], t0=%.10g",poctime_sdate_cnes((*ftsi)[0]+t0,'s'),poctime_sdate_cnes((*ftsi)[ns[i]-1]+t0,'s'),t0);fflush(stdout);
    n=0;
    for(j=0;j<ns[i];j++){
      ftsij=&(*ftsi)[j];
      (*ftsij)+=t0;
      bool tooSoon,tooLate;
      /* allow for *start or *final to be NAN */
      tooSoon=*start>*ftsij;
      tooLate=*ftsij>=*final;
      if(not tooSoon and not tooLate)
/*------------------------------------------------------------------------------
        it counts the frames that are within the time boundaries */
        n++;
      else
/*------------------------------------------------------------------------------
        and discards the frames that are outside */
        *ftsij=NAN;
      }
/*------------------------------------------------------------------------------
    and, if none are kept, marks the file as discarded */
    if(n==0){
      deletep(ftsi);
      continue;
      }
    tn+=n;
    }
  printf("\n");
  if(ts){
    *ts=new double[tn];
    tfi=tn-1;
    }
  
/*------------------------------------------------------------------------------
  It sets start and final to the earliest and latest kept time frames, respectively */
  range_t<double> tR;
  for(i=pn-1;i>=0;i--){//Loop in reverse order...
    if(fts[i]==NULL){
      //...so that deletion does not mess up the values you will process later.
      (*pathlist).erase((*pathlist).begin()+i);
      continue;
      }
    for(j=ns[i]-1;j>=0;j--){//Loop again in reverse order...
      if(ts && !isnan(fts[i][j])){
        (*ts)[tfi]=fts[i][j];//...so that times are sorted if the files are.
        tfi--;
        }
      tR << fts[i][j];
      }
    deletep(&fts[i]);
    }
  deletep(&fts);
  deletep(&ns);
  
  *start=tR.min;
  *final=tR.max;
  
  return tn;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int getSkippedFrames(size_t nt,double *ts,int handleRepeatedFrames,bool **skipFrame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  when times of frames are duplicated, compute which frame to skip
  
Parameters:
  
  + nt : number of frames
  + ts : times of the frames
  + **skipFrame
  + handleRepeatedFrames :
      
      - if 0 :
          
          * if some of the time frames are repeated, return -1
          * otherwise return 0
          
      - if +1: take the last one of duplicates and return 0
      - if -1: take the first one of duplicates and return 0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  if(ts==0)
    return 0;
  
  if(dupos(ts,nt)<0)
    return 0;
  
  int i,j;
  
  if(handleRepeatedFrames)
    *skipFrame=aset(nt,false);
  
  STDOUT_BASE_LINE("Simultaneous time frames :\n");
  for(i=0;i<nt-1;i++) {
    for(j=i+1;j<nt;j++) {
      if(ts[i]==ts[j]) {
        if(*skipFrame!=0){
          if(handleRepeatedFrames<0)/* take first */
            (*skipFrame)[j]=true;
          if(handleRepeatedFrames>0)/* take last */
            (*skipFrame)[i]=true;
          }
        else
          printf("[%d==%d]%s = %.10g\n",i,j,poctime_sdate_cnes(ts[j],'s'),ts[i]);
        break;
        }
      }
    }
  
  if(*skipFrame!=0){
    
    for(i=0;i<nt-1;i++)
      if((*skipFrame)[i])
        printf("will skip [%d]%s = %.10g\n",i,poctime_sdate_cnes(ts[i],'s'),ts[i]);
    
    fflush(stdout);
    return 0;
    }
  
  fflush(stdout);
  return -1;
}
