
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

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Detides radar data. POORLY CODED! DO NOT USE!

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
#include <unistd.h>

#include <stdint.h> //uint64_t
#include <sys/time.h> //gettimeofday

#include "tools-structures.h"

#include "functions.h" //safely includes omp.h
#include "matrix.h"
#include "poc-netcdf.def"
#include "poc-netcdf-data.hpp"
#include "poc-time.h"
#include "filter.h"
#include "map.h"
#include "spectrum.h"
#include "statistic.h"

struct timeval stv;///<start timeval, for progression
///These global variables are used to count the amount of time this programme took for I/O tasks
struct timeval before;
double rt=0.,rt2=0.,rrt=0.,wt=0.;// read, 2nd read, redundant read and write times
double hst=0.,hat=0.,hct=0.;// harmonic storage, analysis and correction times.

#define WAVE_AS_SEPARATE_FILE 0
#define WAVE_AS_DIMENSION 1

#define GLOBAL      0
#define PER_YEAR    1
#define PER_MONTH   2
#define PER_WEEK    3
#define PER_DAY     4
#define PER_HOUR    5
#define PER_MINUTE  6
#define PER_SECOND  7
#define CUSTOM      8

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void increment_month(date_t & actual)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(actual.month<12) 
    actual.month++;
  else {
    actual.month=1;
    actual.year++;
    }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void increment_day(date_t & actual)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(actual.day<poctime_dpm(actual.month,actual.year))
    actual.day++;
  else {
    actual.day=1;
    increment_month(actual);
    }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void increment_seconds(date_t & actual, double step)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(actual.second<24.36000-step)
    actual.second+=step;
  else {
    actual.second+=step-24.*3600.;
    increment_day(actual);
    }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int date_incremental(date_t & actual, int packing)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  status;
  
  switch(packing) {
    case PER_YEAR:
      actual.year++;
      break;
      
    case PER_MONTH:
      increment_month(actual);
      break;
      
    case PER_WEEK:
      for(int k=0;k<7;k++) increment_day(actual);
      break;
      
    case PER_DAY:
      increment_day(actual);
      break;
      
    case PER_HOUR:
      increment_seconds(actual, 3600.);
      break;
      
    case PER_MINUTE:
      increment_seconds(actual, 60.);
      break;
      
    case PER_SECOND:
      increment_seconds(actual, 1.);
      break;
      
    }
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int decode_name(date_t actual, const char *name_template, const char *varname, char **filename, int *packing)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  status;
  char dummy[256], *pointer;
  date_t cdf_reference;
  FILE *out;

/* *----------------------------------------------------------------------------
  Reconstruct the input file name*/
  (*filename) = new char[strlen(name_template) + 256];
  sprintf((*filename), "%s", name_template);

  out = fopen(*filename, "r");
  status = (out == NULL);

  switch (status) {
    case 0:
/* *----------------------------------------------------------------------------
      file exists, do nothing more*/
      fclose(out);
      *packing=GLOBAL;
      break;

    default:
/*-----------------------------------------------------------------------------
      use format convention information*/
/* *----------------------------------------------------------------------------
      substitute YYYY with current year*/
      pointer = strstr((*filename), "YYYY");
      if(pointer != NULL) {
        sprintf(dummy, "%4d", actual.year);
        strncpy(pointer, dummy, 4);
        *packing=PER_YEAR;
        }
/* *----------------------------------------------------------------------------
      substitute MM with current month*/
      pointer = strstr((*filename), "MM");
      if(pointer != NULL) {
        sprintf(dummy, "%2.2d", actual.month);
        strncpy(pointer, dummy, 2);
        *packing=PER_MONTH;
        }
/* *----------------------------------------------------------------------------
      substitute DD with current day*/
      pointer = strstr((*filename), "DD");
      if(pointer != NULL) {
        sprintf(dummy, "%2.2d", actual.day);
        strncpy(pointer, dummy, 2);
        *packing=PER_DAY;
        }
/* *----------------------------------------------------------------------------
      substitute DD with current day*/
      pointer = strstr((*filename), "NNN");
      if(pointer != NULL) {
        int day=day_in_year(actual);
        sprintf(dummy, "%3.3d", day);
        strncpy(pointer, dummy, 3);
        *packing=PER_DAY;
        }
/* *----------------------------------------------------------------------------
      substitute HH with current hour*/
      pointer = strstr((*filename), "HH");
      if(pointer != NULL) {
        sprintf(dummy, "%2.2d", (int) floor(actual.second/3600.));
        strncpy(pointer, dummy, 2);
        *packing=PER_HOUR;
        }
/* *----------------------------------------------------------------------------
      substitute MN with current hour*/
      pointer = strstr((*filename), "MN");
      if(pointer != NULL) {
        double hour=floor(actual.second/3600.);
        sprintf(dummy, "%2.2d", (int) floor((actual.second-3600.*hour)/60.));
        strncpy(pointer, dummy, 2);
        *packing=PER_HOUR;
        }
// /* *----------------------------------------------------------------------------
//       lookup file fitting ???? with minute and seconds index*/
//       pointer = strstr((*filename), "????");
//       if(pointer != NULL) {
//         int i,j;
//         for(j=0;j<60;j++) {
//           for(i=0;i<60;i++) {
//             sprintf(dummy, "%2.2d%2.2d", j,i);
//             strncpy(pointer, dummy, 4);
//             status= get_file_size(*filename, 0);
// /* *----------------------------------------------------------------------------
//             status=0 if file found*/
//             if(status==0) break;
//             }
//           if(status==0) break;
//           }
//         }
/* *----------------------------------------------------------------------------
      substitute VARNAME*/
      pointer = strstr((*filename), "VARNAME");
      if(pointer != NULL) {
        sprintf(dummy, "%s", varname);
        strncpy(pointer, dummy, strlen(varname));
        }
      break;

      break;
    }

  return (status);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int save_SGXYT(const char *output, grid_t grid, T *buffer, T mask, const char *varname,const char *unit, const char *standard, int create)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,status;
  fcomplex z;
  cdfvar_t variable;
  pocgrd_t ncgrid;
  date_t origin;
  
  printf("#################################################################\n");
  printf("write output file : %s\n",output);
  if(create==1) {
    status= poc_createfile(output);
    grid.nz=1;
    grid.z=NULL;
    status=poc_sphericalgrid_xyzt(output,"","",grid,&ncgrid);
    if(grid.time!=0) {
      variable= poc_standardtime("time","nt","seconds", origin);
      status=create_ncvariable(output, &variable);
      for(count=0; count<grid.nt; count++) status=poc_writetime(output, count, variable.id,  grid.time[count]);
      variable.destroy();
      }
    }
  else {
//    status=poc_sphericalgrid_xyzt(output,&ncgrid);
    }

  status=poc_standardvariable_xyt(varname,mask,unit,(T)1.,(T)0.,standard,standard,standard,ncgrid,variable);
  status=create_ncvariable(output, &variable);
  status=poc_put_var(output, variable.id, buffer);
  variable.destroy();

}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int save_SGXYT(const char *output, grid_t grid, T *buffer, T mask, const char *varname,const char *unit, const char *standard, 
				     int createfile, int creategrid, int createvariable, int frame, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,status;
  fcomplex z;
  cdfvar_t variable;
  pocgrd_t ncgrid;
  
  poc_data_t<double> lon,lat,time;
  poc_dim_t ncycles, nx, ny;
  poc_var_t vdum,vresiduals;
    
  date_t origin;
  
  if(verbose==1) {
    printf("#################################################################\n");
    printf("write %s (frame %d) in output file : %s\n",varname, frame, output);
    }
  if(createfile==1) {
    status= poc_createfile(output);
    }
    
  if(creategrid==1) {
    grid.nz=1;
    grid.z=NULL;
    status=poc_sphericalgrid_xyt(output,"",grid,&ncgrid);
    if(grid.time!=0) {
      variable= poc_standardtime("time","t","seconds", origin);
      status=create_ncvariable(output, &variable);
      for(count=0; count<grid.nt; count++) status=poc_writetime(output, count, variable.id,  grid.time[count]);
      variable.destroy();
      }
//     lon.write_data(output);
//     lat.write_data(output);
//     for(size_t m=0; m<ncycles.len;m++) {
//       time.write_data(output,m);
//       }
    }
  else {
    if(createvariable==1) status=poc_inq_grid(output,"","",&ncgrid);
    }

  if(createvariable==1) {
//     vdum.clear();
//     vdum.init(varname,NC_FLOAT,standard,"units",mask);
//     vdum.dimensions<<vresiduals.dimensions[1];
//     vdum.dimensions<<vresiduals.dimensions[2];
    status=poc_standardvariable_xyt(varname,mask,unit,(T)1.,(T)0.,standard,standard,standard,ncgrid,variable);
    status=create_ncvariable(output, &variable);
    }
  else {
    status=cdf_varinfo(output,varname,&variable,0);
//     status=poc_inq_var(output,varname,&vdum,0);
    }

  status=poc_write_xyt(output,  grid, frame, variable.id, buffer);
  variable.destroy();
//     status=poc_put_var(output,vdum,buffer,frame,1);
//     vdum.clear();
    
    return(status);
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// template <typename T> int save_SGXYT(const char *output, grid_t grid, T *buffer, T mask, const char *varname,const char *unit, const char *standard, 
// 				     int createfile, int creategrid, int createvariable, int frame, int verbose)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int count,status;
//   fcomplex z;
//   cdfvar_t variable;
//   pocgrd_t ncgrid;
//   
//   poc_data_t<double> lon,lat,time;
//   poc_dim_t ncycles, nx, ny;
//   poc_var_t vdum,vresiduals;
//     
//   date_t origin;
//   
//   if(verbose==1) {
//     printf("#################################################################\n");
//     printf("write %s (frame %d) in output file : %s\n",varname, frame, output);
//     }
//   if(createfile==1) {
//     status= poc_createfile(output);
//     }
//     
//   if(creategrid==1) {
//     lon.write_data(output);
//     lat.write_data(output);
//     for(size_t m=0; m<ncycles.len;m++) {
//       time.write_data(output,m);
//       }
//     }
//   else {
//     if(createvariable==1) status=poc_inq_grid(output,"","",&ncgrid);
//     }
// 
//   if(createvariable==1) {
//     vdum.clear();
//     vdum.init(varname,NC_FLOAT,standard,"units",mask);
//     vdum.dimensions<<vresiduals.dimensions[1];
//     vdum.dimensions<<vresiduals.dimensions[2];
//     }
//   else {
//     status=poc_inq_var(output,varname,&vdum,0);
//     }
// 
//     status=poc_put_var(output,vdum,buffer,frame,1);
//     vdum.clear();
//     
//     return(status);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int detide_radar(const string & rootname,const string & source, spectrum_t s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status = NC_NOERR + 1;
  int ncid;
 
  string filename;
  poc_dim_t ncycles, nx, ny;
  
  vector<poc_dim_t> dim1D, dim2D;

  grid_t grid;
  float *amp,*pha;
  
  vector<serie_t> track;
  
  poc_data_t<double> lon,lat,azimuth,valids;
  poc_var_t vtime;
  poc_var_t var_u, var_v;
  poc_global_t info;
  
  status=poc_inq(source,&info);
  if(status!=0) NC_TRAP_ERROR(wexit,status,1,"poc_inq(\""+source+"\",) error");
 
  status=lon.init(info,"lon");
  status=lat.init(info,"lat");
  status=azimuth.init(info,"azimuth");
  status=valids.init(info,"count");
  status=info.variables.find("time",  &vtime);
  
  poc_print(vtime);
      
  ncycles=vtime.dimensions[0];
  
  ny=lon.info.dimensions[0];
  nx=lon.info.dimensions[1];

//   dim1D.push_back(npoints);
//   
//   dim2D.push_back(ncycles);
//   dim2D.push_back(npoints);
 
/**----------------------------------------------------------------------------
  get longitudes */
  lon.read_data(source);
  
/**----------------------------------------------------------------------------
  get latitudes */
  lat.read_data(source);
    
/**----------------------------------------------------------------------------
  get latitudes */
  azimuth.read_data(source);
  
//return 0;
    
/**----------------------------------------------------------------------------
  get nvalids */
  valids.read_data(source);
  
  printf("\n\n");
  poc_print(valids.info);
  
  int verbose=0;
  status=nc_open(source.c_str(),NC_NOWRITE,&ncid);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"nc_open(\""+filename+"\",NC_NOWRITE,) error");
    

/**----------------------------------------------------------------------------
  get time */
  poc_data_t<double> time;
  
  time.data=new double[ncycles.len];
  size_t start[1]={0};
  size_t count[1]={ncycles.len};
  status=poc_get_vara(ncid,vtime.id,start,count,time.data);

/**----------------------------------------------------------------------------
  get u */
  poc_data_t<double> u;
  status=poc_inq_var(source,"u",  &var_u);

  u.init(var_u);
  
  u.data=new double[ncycles.len];
  
  grid.nx=nx.len;
  grid.ny=ny.len;
  grid.nt=ncycles.len;
  
  grid.x=lon.data;
  grid.y=lat.data;
  grid.time=time.data;
  
  grid.modeH=2;
  float fspec=u.mask;

  double *residuals=new double[grid.Hsize()*ncycles.len];
  
  for(int k=0; k<ncycles.len*grid.Hsize(); k++) residuals[k]=u.mask;
 
  amp=new float[s.n*grid.Hsize()];
  pha=new float[s.n*grid.Hsize()];
  
  for(int k=0; k<s.n*grid.Hsize(); k++) amp[k]=fspec;
  for(int k=0; k<s.n*grid.Hsize(); k++) pha[k]=fspec;

  float *reduction =new float[grid.Hsize()];
  
  float *std_raw   = new float[grid.Hsize()];
  float *std_res   = new float[grid.Hsize()];
  float *mean_raw  = new float[grid.Hsize()];
  float *mean_res  = new float[grid.Hsize()];
  
  
  int analysed=0;
  
  int nprocs __attribute__((unused)) =initialize_OPENMP(-1);
  for(int j=0;j<ny.len;j++) {
// #pragma omp parallel for reduction(+:analysed) if(nprocs>1)
    for(int i=0;i<nx.len;i++) {
      int n=j*nx.len+i;
      
      amp[n]=fspec;
      pha[n]=fspec;
      reduction[n]=fspec;
      std_raw[n]=fspec;
      std_res[n]=fspec;
      mean_raw[n]=fspec;
      mean_res[n]=fspec;
     
      double mean;
      spectrum_t solved;
      int nodal=1;
      size_t start[3]={0,j,i};
      size_t count[3]={ncycles.len,1,1};
      if(valids.data[n]>2500) {
        statistic_t s1,s2;
        printf("i=%3d j=%3d n=%6d / %6d analysed=%6d\n",i,j,n,grid.Hsize(),analysed);
        double *serie=new double[ncycles.len];
//       #if NETCDF_CAN_PARALLEL_IO == 0
//       #pragma omp critical(threadUnsafeNetcdf)
//       #endif
        {
        status=poc_get_vara(ncid,var_u.id,start,count,serie);
        }
        s1=get_statistics(serie,u.mask, ncycles.len,0);
        double *buffer=new double[ncycles.len];
        mgr_data_t *tmp=harmonic_analysis_with_parsing(serie,u.mask,time.data, buffer, &mean, ncycles.len, s, solved, nodal, 0);
        for(int k=0; k<ncycles.len; k++) if(buffer[k]!=u.mask) residuals[k*grid.Hsize()+n]=buffer[k]+mean;
        s2=get_statistics(buffer,u.mask, ncycles.len,0);
        for(int k=0;k<solved.n;k++) {
          int kk=s.wave_index(solved.waves[k].name);
          amp[kk*grid.Hsize()+n]=tmp[k].amp;
          pha[kk*grid.Hsize()+n]=tmp[k].phi;
          }
        reduction[n]=100*(s1.std-s2.std)/s1.std;
        std_raw[n]=s1.std;
        std_res[n]=s2.std;
        mean_raw[n]=s1.mean;
        mean_res[n]=s2.mean+mean;
        delete[] buffer;
        delete[] serie;
//         amp[n]=n;
//         pha[n]=n;
        analysed++;
        delete[] tmp;
        }
      }
    }
    
  filename=rootname+"-analysis.nc";
  status=save_SGXYT(filename.c_str(), grid, residuals, u.mask, "residuals", "m/s", "residuals", 1);
  
  status=save_SGXY(filename.c_str(), grid, valids.data, u.mask, "count", "none", "count", 0);
  status=save_SGXY(filename.c_str(), grid, reduction, fspec, "reduction", "none", "reduction", 0);
  status=save_SGXY(filename.c_str(), grid, mean_raw,  fspec, "mean_raw", "m/s", "mean_raw", 0);
  status=save_SGXY(filename.c_str(), grid, std_raw,   fspec, "std_raw", "m/s", "std_raw", 0);
  status=save_SGXY(filename.c_str(), grid, mean_res,  fspec, "mean_res", "m/s", "mean_res", 0);
  status=save_SGXY(filename.c_str(), grid, std_res,   fspec, "std_res", "m/s", "std_res", 0);
  
  azimuth.write_data(filename.c_str());
  
  for(int k=0;k<s.n;k++) {
    bool go_ahead=false;
    for(int n=0; n<grid.Hsize(); n++) {
      if(amp[k*grid.Hsize()+n]!=fspec) {
        go_ahead=true;
        break;
        }
      }
    if(!go_ahead) continue;
    filename=rootname+"-" +s.waves[k].name +".nc";
    status=save_SGXY(filename.c_str(), grid, &amp[k*grid.Hsize()], fspec, "Ua", "m/s",     "Ua", 1);
    status=save_SGXY(filename.c_str(), grid, &pha[k*grid.Hsize()], fspec, "Ug", "degrees", "Ug", 0);
    azimuth.write_data(filename.c_str());
    }
    
  nc_close(ncid);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int decode_time(const poc_var_t & vtime, date_t & origine, char & units, double & factor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, k;
  int day, month, year;
  int hour, minute, seconds;
  char smonth[1024];
  char *text, *pointer, *time_units;
  
  year=month=day=0;
  hour=minute=seconds=0;

/*------------------------------------------------------------------------
  time origin */
  const poc_att_t *att;
  att=vtime.attributes.findP("time_origin");
  if(att!=0){
    text=strdup(att->as_charp());
    }
  else
/*------------------------------------------------------------------------
    no time origin attribute found, let's try units*/
    goto units;
  
/*------------------------------------------------------------------------
  parse time origin attribute */
  sscanf(strtok(text, "-"), "%d", &year);
  if((pointer = strtok(NULL, "-")) == NULL)
    goto error;
  sscanf(pointer, "%3s", &smonth);
  if((pointer = strtok(NULL, " ")) == NULL)
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
/*------------------------------------------------------------------------
  patch: month given as number and not text*/
  if(k == 12)
    sscanf(smonth, "%d", &month);
  if(month==13){
    status=NC_ECHAR;/* Attempt to convert between text & numbers */
    nc_check_error(status,__LINE__,__FILE__,"Can not interpret month name %s",smonth);
    return -1;
    }

  free(text);
    
units:
  factor = 1.;
  att=vtime.attributes.findP("units");
  if(att!=0) {
    text=strdup(att->as_charp());
    }
  else {
/*------------------------------------------------------------------------
    assume undocumented time to be Cnes time */
    year=1950;
    month=1;
    day=1;
    hour=minute=seconds=0;
    factor = 24 * 3600.;
    goto load;
    }
    
  time_units = strdup(strtok(text, " "));
  
//   units=tolower(time_units[0]);
  
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
      return -1;
    }
    
  text=strdup(att->as_charp());

  if( (pointer = strstr(text, "since")) != NULL) {
    pointer += 6;
    pointer = strtok(pointer, "-");
    if(pointer != NULL)
      sscanf(pointer, "%d", &year);
    pointer = strtok(NULL, "-");
    if(pointer != NULL)
      sscanf(pointer, "%d", &month);
    pointer = strtok(NULL, " ");
    if(pointer != NULL)
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
    }
    
  free(text);

/*------------------------------------------------------------------------
  assume time is 1D buffer */
load:

  if(year==0){//year 0 does not exit : year 1 follows year -1
    origine=NADate;
    }
  else{
    origine.year=year;
    origine.month=month;
    origine.day=day;
    origine.second=hour*3600+minute*60+seconds;
    }

  return (0);

error:
  return (-1);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int detide_radar_new(const string & rootname,const string & source, spectrum_t s, int nvalidmin, const char *vflags, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status = NC_NOERR + 1;
  int ncid;
  bool use_flags;
  string filename;
  poc_dim_t ncycles, nx, ny;
  
  vector<poc_dim_t> dim1D, dim2D;

  grid_t grid;
  float *amp=0,*pha=0;
  bool  *gflags=0;
  
  vector<serie_t> track;
  
  poc_data_t<double> lon,lat,azimuth,flags;
  poc_var_t vtime;
  poc_var_t var_u, var_v;
  poc_global_t info;
  
  status=poc_inq(source,&info);
  if(status!=0) NC_TRAP_ERROR(wexit,status,1,"poc_inq(\""+source+"\",) error");
 
  status=lon.init(info,"longitude");
  status=lat.init(info,"latitude");
  status=azimuth.init(info,"hcdt");
  
  if(vflags!=0) {
    status=flags.init(info,vflags);
    if(status==0) use_flags=true;
    else TRAP_ERR_EXIT(status,"fatal error: could not init variable %s from %s\n\n",vflags,source.c_str());
    }
  else {
    use_flags=false;
    }
   
  status=info.variables.find("time",  &vtime);
  
  if(verbose) poc_print(vtime);
      
  ncycles=vtime.dimensions[0];
  
  if(lon.info.dimensions.size()==2) {
    ny=lon.info.dimensions[0];
    nx=lon.info.dimensions[1];
    }
  else {
    ny=lat.info.dimensions[0];
    nx=lon.info.dimensions[0];
    }

/**----------------------------------------------------------------------------
  get longitudes */
  lon.read_data(source);
  
/**----------------------------------------------------------------------------
  get latitudes */
  lat.read_data(source);
    
/**----------------------------------------------------------------------------
  get azimuth */
  azimuth.read_data(source);
      
/**----------------------------------------------------------------------------
  get flags at all frames*/
  if(use_flags) {
    gflags=new bool[nx.len*ny.len*ncycles.len];
    for(size_t k=0; k<ncycles.len;k++) {
      flags.read_data(source,k);
      for(size_t n=0; n<nx.len*ny.len; n++) gflags[k*nx.len*ny.len+n]=(flags.data[n]==1);
      }
    
    if(verbose) {
      printf("\n\n");
      poc_print(flags.info);
      }
    }
  
  status=nc_open(source.c_str(),NC_NOWRITE,&ncid);
  if(status!=NC_NOERR) {
    NC_TRAP_ERROR(return,status,verbose,"nc_open(\""+filename+"\",NC_NOWRITE,) error");
    }
    

/**----------------------------------------------------------------------------
  get time */
  poc_data_t<double> time;
  
  time.info=vtime;
  
  time.data=new double[ncycles.len];
  size_t start[1]={0};
  size_t count[1]={ncycles.len};
  status=poc_get_vara(ncid,vtime.id,start,count,time.data);
  
  date_t origin;
  double factor;
  char units='\0';
  status=decode_time(vtime, origin, units, factor);
  
/**----------------------------------------------------------------------------
  harmonic analysis expect cnes time in days*/
  double t0=cnes_time(origin,'d');
  double *cnestime=new double[ncycles.len];
  for(size_t m=0;m<ncycles.len;m++) cnestime[m]=time.data[m];
  poctime_convert(cnestime, ncycles.len, units, 'd');
  for(size_t m=0;m<ncycles.len;m++) cnestime[m]+=t0;
  
/**----------------------------------------------------------------------------
  get u */
  poc_data_t<float> u;
  status=poc_inq_var(source,"hcsp",  &var_u);

  u.init(var_u);
  
  u.data=new float[ncycles.len];
  
  grid.nx=nx.len;
  grid.ny=ny.len;
  grid.nt=ncycles.len;
  
  grid.x=lon.data;
  grid.y=lat.data;
//   grid.time=time.data;
//  grid.time=cnestime;
  grid.time=0;
  
  if(lon.info.dimensions.size()==2) {
    grid.modeH=2;
    }
  else {
    grid.modeH=1;
    }
    
  float  fspec=u.mask;
  double mask =u.mask;

  float *residuals=new float[grid.Hsize()*ncycles.len];
  
  for(int k=0; k<ncycles.len*grid.Hsize(); k++) residuals[k]=fspec;
 
  amp=new float[s.n*grid.Hsize()];
  pha=new float[s.n*grid.Hsize()];
  
  for(int k=0; k<s.n*grid.Hsize(); k++) amp[k]=fspec;
  for(int k=0; k<s.n*grid.Hsize(); k++) pha[k]=fspec;

  float *gAzimuth = new float[grid.Hsize()];
  for(int k=0; k<grid.Hsize(); k++) gAzimuth[k]=fspec;
  
  float *reduction = new float[grid.Hsize()];
  for(int k=0; k<grid.Hsize(); k++) reduction[k]=fspec;
  
  float *nvalids   = new float[grid.Hsize()];
  for(int k=0; k<grid.Hsize(); k++) nvalids[k]=fspec;
 
  float *std_raw   = new float[grid.Hsize()];
  float *std_res   = new float[grid.Hsize()];
  float *mean_raw  = new float[grid.Hsize()];
  float *mean_res  = new float[grid.Hsize()];
  
  for(int k=0; k<grid.Hsize(); k++) std_raw[k] =fspec;
  for(int k=0; k<grid.Hsize(); k++) std_res[k] =fspec;
  for(int k=0; k<grid.Hsize(); k++) mean_raw[k]=fspec;
  for(int k=0; k<grid.Hsize(); k++) mean_res[k]=fspec;
  
  int analysed=0;
  int *valid=new int[ncycles.len];
  
  int nprocs __attribute__((unused)) =initialize_OPENMP(-1);
  
  for(int j=0;j<ny.len;j++) {
//   for(int j=40;j<50;j++) {
// #pragma omp parallel for reduction(+:analysed) if(nprocs>1)
    for(int i=0;i<nx.len;i++) {
      int n=j*nx.len+i;
      
      amp[n]=fspec;
      pha[n]=fspec;
      reduction[n]=fspec;
      std_raw[n]=fspec;
      std_res[n]=fspec;
      mean_raw[n]=fspec;
      mean_res[n]=fspec;
     
      double mean;
      spectrum_t solved;
      statistic_t s1,s2;
      int nodal=1;
      size_t start[3]={0,j,i};
      size_t count[3]={ncycles.len,1,1};
      nvalids[n]=0;
      if(use_flags) {
/*------------------------------------------------------------------------------
        first control on #valid observations */
//       status=poc_get_vara(ncid,flags.info.id,start,count,valid);
//       for(size_t k=0; k<ncycles.len;k++) if (flags[k*nx.len*ny.len+n]==1) nvalids++;
        for(size_t k=0; k<ncycles.len;k++) valid[k]=(gflags[k*nx.len*ny.len+n]==1);
        for(size_t k=0; k<ncycles.len;k++) if (valid[k]==1) nvalids[n]++;
        if(nvalids[n]<nvalidmin) continue;
        }
      double *serie=new double[ncycles.len];
      double *a=    new double[ncycles.len];
//       #if NETCDF_CAN_PARALLEL_IO == 0
//       #pragma omp critical(threadUnsafeNetcdf)
//       #endif
      {
      status=poc_get_vara(ncid,var_u.id,start,count,serie);
      status=poc_get_vara(ncid,azimuth.info.id,start,count,a);
      }
      if(use_flags)
        for(size_t k=0; k<ncycles.len;k++) {
          if (valid[k]==0 and serie[k]!=mask) {
            serie[k]=mask;
            }
          }
      if(!use_flags)
/*------------------------------------------------------------------------------
        second control on #valid observations */
        for(size_t k=0; k<ncycles.len;k++) {
          if (serie[k]!=mask) nvalids[n]++;
          }
      if(nvalids[n]<nvalidmin) continue;
      printf("i=%3d j=%3d n=%6d / %6d analysed=%6d\r",i,j,n,grid.Hsize(),analysed);
      fflush(stdout);
      
/*------------------------------------------------------------------------------
      aggregate valid azimuth  */
      for(size_t k=0; k<ncycles.len;k++) if (a[k]!=mask) gAzimuth[n]=a[k];
      
/*------------------------------------------------------------------------------
      full signal stats  */
      s1=get_statistics(serie, mask, ncycles.len,0);
      
/*------------------------------------------------------------------------------
      process harmonic analysis */
      double *buffer=new double[ncycles.len];
      mgr_data_t *tmp=harmonic_analysis_with_parsing(serie, mask, cnestime, buffer, &mean, ncycles.len, s, solved, nodal, 0);
      
      for(int k=0; k<ncycles.len; k++) if(buffer[k]!=u.mask) residuals[k*grid.Hsize()+n]=buffer[k]+mean;
/*------------------------------------------------------------------------------
      detided signal stats  */
      s2=get_statistics(buffer,u.mask, ncycles.len,0);
      for(int k=0;k<solved.n;k++) {
        int kk=s.wave_index(s.waves[k].name);
        amp[kk*grid.Hsize()+n]=tmp[k].amp;
        pha[kk*grid.Hsize()+n]=tmp[k].phi;
        }
/*------------------------------------------------------------------------------
      signal reduction (percent) */
      reduction[n]=100*(s1.std-s2.std)/s1.std;
/*------------------------------------------------------------------------------
      signal statistics */
      std_raw[n]=s1.std;
      std_res[n]=s2.std;
      mean_raw[n]=s1.mean;
      mean_res[n]=s2.mean+mean;
      delete[] buffer;
      delete[] serie;
      delete[] a;
//         amp[n]=n;
//         pha[n]=n;
      analysed++;
      delete[] tmp;
      }
    }
  delete[] valid;
  delete[] gflags;
    
  filename=rootname+"-analysis.nc";
  
  poc_var_t vresiduals=var_u;
  poc_var_t vdum;
    
  status=poc_createfile(filename.c_str());
  
  lon.write_data(filename.c_str());
  lat.write_data(filename.c_str());

  for(size_t m=0; m<ncycles.len;m++) {
    time.read_data(ncid,m);
    time.write_data(filename.c_str(),m);
    }
    
  for(size_t m=0; m<ncycles.len;m++) {
    azimuth.read_data(ncid,m);
    azimuth.write_data(filename.c_str(),m);
    }
 
  vresiduals.init("residuals",NC_FLOAT,"residuals","m/s",u.mask);
  status=poc_put_var(filename.c_str(),vresiduals,residuals,1);
  
  vdum.clear();
  vdum.init("count",NC_FLOAT,"count","none",fspec);
  vdum.dimensions<<vresiduals.dimensions[1];
  vdum.dimensions<<vresiduals.dimensions[2];
  status=poc_put_var(filename.c_str(),vdum,nvalids,1);
  
  vdum.clear();
  vdum.init("mean_raw",NC_FLOAT,"mean_raw","m/s",fspec);
  vdum.dimensions<<vresiduals.dimensions[1];
  vdum.dimensions<<vresiduals.dimensions[2];
  status=poc_put_var(filename.c_str(),vdum,mean_raw,1);
  
  vdum.clear();
  vdum.init("std_raw",NC_FLOAT,"std_raw","m/s",fspec);
  vdum.dimensions<<vresiduals.dimensions[1];
  vdum.dimensions<<vresiduals.dimensions[2];
  status=poc_put_var(filename.c_str(),vdum,std_raw,1);
  
  vdum.clear();
  vdum.init("mean_res",NC_FLOAT,"mean_res","m/s",fspec);
  vdum.dimensions<<vresiduals.dimensions[1];
  vdum.dimensions<<vresiduals.dimensions[2];
  status=poc_put_var(filename.c_str(),vdum,mean_res,1);
  
  vdum.clear();
  vdum.init("std_res",NC_FLOAT,"std_res","m/s",fspec);
  vdum.dimensions<<vresiduals.dimensions[1];
  vdum.dimensions<<vresiduals.dimensions[2];
  status=poc_put_var(filename.c_str(),vdum,std_res,1);
  
  vdum.clear();
  vdum.init("azimuth",NC_FLOAT,"azimuth","degree",fspec);
  vdum.dimensions<<vresiduals.dimensions[1];
  vdum.dimensions<<vresiduals.dimensions[2];
  status=poc_put_var(filename.c_str(),vdum,gAzimuth,1);
  
  vdum.clear();
  vdum.init("reduction",NC_FLOAT,"reduction","percent",fspec);
  vdum.dimensions<<vresiduals.dimensions[1];
  vdum.dimensions<<vresiduals.dimensions[2];
  status=poc_put_var(filename.c_str(),vdum,reduction,1);
  
  for(int k=0;k<s.n;k++) {
    bool go_ahead=false;
    for(int n=0; n<grid.Hsize(); n++) {
      if(amp[k*grid.Hsize()+n]!=fspec) {
        go_ahead=true;
        break;
        }
      }
    if(!go_ahead) continue;
    filename=rootname+"-" +s.waves[k].name +".nc";
    status=save_SGXY(filename.c_str(), grid, &amp[k*grid.Hsize()], fspec, "Ua", "m/s",     "Ua", 1);
    status=save_SGXY(filename.c_str(), grid, &pha[k*grid.Hsize()], fspec, "Ug", "degrees", "Ug", 0);
    status=save_SGXY(filename.c_str(), grid, gAzimuth, fspec, "azimuth", "degrees", "azimuth", 0);
    }
    
  nc_close(ncid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int combine_tide(const string & a, const string & b, const char *filename, const char *gridfile, const char *varname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, status;
  grid_t RadarGrid[2], grid;
  complex<float> *Unative[2], mask;
  complex<float> *U, *V;
  const string source[2]={a,b};
  double *Anative[2], dmask;
//   poc_var_t vlon, vlat;
  poc_global_t info;
  poc_data_t<double> azimuth[2];
  
  mask=complex<float> (9999,9999);
  dmask=9999;
  
  status=poc_inq(a,&info);
  if(status!=0) NC_TRAP_ERROR(wexit,status,1,"poc_inq(\""+a+"\",) error");
 
//   status=poc_inq_var(a,"lon",   &vlon);
//   status=poc_inq_var(a,"lat",   &vlat);
  
  status=poc_get_grid(a,"Ua",&RadarGrid[0],1);
  status=poc_get_grid(b,"Ua",&RadarGrid[1],1);
  
  if(gridfile==0) {
    grid.xmin=min(RadarGrid[0].xmin,RadarGrid[1].xmin);
    grid.xmax=min(RadarGrid[0].xmax,RadarGrid[1].xmax);

    grid.ymin=min(RadarGrid[0].ymin,RadarGrid[1].ymin);
    grid.ymax=min(RadarGrid[0].ymax,RadarGrid[1].ymax);
  
    map_set2Dgrid(&grid, grid.xmin, grid.ymin, grid.xmax, grid.ymax,1./120.,1./120.);
    }
  else {
    status=poc_get_grid(gridfile,varname,&grid,1);
    }
    
  for(int k=0;k<2;k++) {
    poc_data_t<float> Ua, Ug;
    status=poc_inq_var(source[k],"Ua",  &Ua.info);
    Ua.init();
    Ua.read_data(source[k]);
    status=poc_inq_var(source[k],"Ug",  &Ug.info);
    Ug.init();
    Ug.read_data(source[k]);
    status=poc_inq_var(source[k],"azimuth",  &azimuth[k].info);
    azimuth[k].init();
    azimuth[k].read_data(source[k]);
    Unative[k]=new complex<float> [RadarGrid[k].Hsize()];
    Anative[k]=new double [RadarGrid[k].Hsize()];
    for(n=0;n<RadarGrid[k].Hsize();n++) {
      if(Ua.data[n]!=Ua.mask) {
        Unative[k][n]=Ua.data[n]*complex<float>(cos(-Ug.data[n]*M_PI/180.),sin(-Ug.data[n]*M_PI/180.));
        Anative[k][n]=azimuth[k].data[n]*M_PI/180.0;
        }
      else {
        Unative[k][n]=mask;
//        Anative[k][n]=mask;
        Anative[k][n]=azimuth[k].data[n]*M_PI/180.0;
        }
      }
    Ua.destroy_data();
    Ug.destroy_data();
    }


  U=new complex<float> [grid.Hsize()];
  V=new complex<float> [grid.Hsize()];

  for(int j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      complex<float> Uradial[2];
      double Aradial[2];
      n=grid.Hindex(i,j);
      double x,y;
      grid.xy(i,j,x,y);
/*------------------------------------------------------------------------------
      interpolation on regular grid from radar sampling grid*/
      status=map_interpolation(RadarGrid[0], Unative[0], mask,  x, y, &Uradial[0]);
      status=map_interpolation(RadarGrid[0], Anative[0], dmask, x, y, &Aradial[0]);
      status=map_interpolation(RadarGrid[1], Unative[1], mask,  x, y, &Uradial[1]);
      status=map_interpolation(RadarGrid[1], Anative[1], dmask, x, y, &Aradial[1]);
      if( (Uradial[0]!=mask) && (Uradial[1]!=mask)) {
/*------------------------------------------------------------------------------
        reconstruction of the 2D vector*/
/*------------------------------------------------------------------------------
        eastward-northward direction = pi/2 - azimuth (north/east=90°/0° azimuth=0°/90° */
        vector2D_t u(M_PI/2.-Aradial[0]), v(M_PI/2.-Aradial[1]);
/*------------------------------------------------------------------------------
        apparently velocity was given toward the radar, now corrected */
//         cvector2D_t w=math_vector_coordinates03((complex<double>) -Uradial[0], u, (complex<double>) -Uradial[1], v);
        cvector2D_t w=math_vector_coordinates03((complex<double>) Uradial[0], u, (complex<double>) Uradial[1], v);
        U[n]=w.x;
        V[n]=w.y;
        }
      else {
        U[n]=mask;
        V[n]=mask;
        }
      }
    }
    
  pocgrd_t ncgrid;
    
  status= poc_createfile(filename);
  
  status=map_completegridaxis(&grid,2);
  
  status=poc_sphericalgrid_xy(filename,"",grid,&ncgrid);
  
  status=tides_savec1(filename, "Ua", "Ug", "currents", "m/s", grid, U, mask, ncgrid);
  status=tides_savec1(filename, "Va", "Vg", "currents", "m/s", grid, V, mask, ncgrid);

  return 0;
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int loadRadarRaw_ascii(const char *filename, const char *name, double ref_lon, double ref_lat, double azimuth_shift, date_t start, date_t end, int nx)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

  Load RADAR-HF ascii data and rewrite them under netcdf format
  
-----------------------------------------------------------------------------*/
{
  int k, l, i, n, m, nitems, keep, status;
  int ny,nt;
  double *u, *mask, *count, spec=9999.;
  double *time, *lon, *lat,t0, *azimuth, *radius;
  grid_t grid;
  FILE *in;
  
  cout << endl << "-------- starting conversion --------" << endl << endl;
  
//  ny=mgr_line_count(filename);
  
  ny=71;
  
  double delta=ellapsed_time(start,end,'s');
  nt=delta/1200+1;
  
//  nt=70*24*3;

  lon=new double[nx*ny];
  lat=new double[nx*ny];
  time=new double[nt];
    
  u=new double[nx*ny*nt];

  for(n=0;n<nt*nx*ny;n++) {
    u[n]=spec;
    }
  
  mask=new double[nx*ny];
  
  count=new double[nx*ny];
  
  azimuth=new double[ny];
  radius=new double[nx];
 
  double x0,y0;
  
//  end  =date_t(2007,10,27,0.);
  
  t0=cnes_time(start,'s');
  
  string name_template="YYYYNNNHHMN_"+(string)name+".coc";
  
  for(i=0;i<nt;i++) {
    time[i]=t0+i*1200.;
    date_t present = poctime_getdatecnes(time[i], 's');
    char *frame, keyword[1024];
    int packing;
    status=decode_name(present, name_template.c_str(), 0, &frame, &packing);
    time[i]/=24.*3600.;
    in=fopen(frame,"r");
    if(in==0) continue;
//    printf("reading %s\n",frame);
    fscanf(in,"%s",keyword);
    for(k=0;k<ny;k++) {
      fscanf(in,"%lf",&azimuth[k]);
      }
   
    for(l=0;l<nx;l++) {
      nitems=fscanf(in,"%lf",&radius[l]);
      if(nitems!=1) break;
      for(k=0;k<ny;k++) {
        n=k*nx+l;
        m=i*nx*ny+n;
        fscanf(in,"%s",keyword);
        nitems=sscanf(keyword,"%lf",&u[m]);
        if(nitems!=1) {
          u[m]=spec;
          }
        else if(isnan(u[m])!=0) {
          u[m]=spec;
          }
        else
/*------------------------------------------------------------------------------
          minus sign apply to reverse direction (initially radard-ward) */
          u[m]=-u[m];
        }
      }
    fclose(in);
    delete[] frame;
    }
  
  projPJ projection;
  projection=assign_StereoOblique(ref_lat,ref_lon);
  
  geo_to_projection(projection, ref_lat,ref_lon,&x0,&y0);
  
  for(k=0;k<ny;k++) {
    for(l=0;l<nx;l++) {
/*------------------------------------------------------------------------------
      +180.0 added to reverse direction (initially radard-ward) */
      double alpha=(azimuth[k]+azimuth_shift+180.)*M_PI/180.0;
      double x=x0+radius[l]*cos(alpha)*1000.;
      double y=y0+radius[l]*sin(alpha)*1000.;
      n=k*nx+l;
      double t,p;
      projection_to_geo(projection, &p, &t,x,y);
      lon[n]=t;
      lat[n]=p;
      }
    }
  grid.nx=nx;
  grid.ny=ny;
  grid.nt=nt;
  grid.modeH=2;
  
  grid.x=lon;
  grid.y=lat;
  grid.time=time;
  
  status=save_SGXYT(filename, grid, u, spec, "u", "m/s", "u", 1);
    
  for(n=0;n<nx*ny;n++) {
    keep=0;
    for(size_t k=0;k<nt;k++) {
      if(u[k*nx*ny+n]!=spec) {
        keep++;
        }
      }
    count[n]=keep;
    }
  status=save_SGXY(filename, grid, count, spec, "count", "none", "count", 0);

  for(k=0;k<ny;k++) {
    for(l=0;l<nx;l++) {
      n=k*nx+l;
/*------------------------------------------------------------------------------
      direction relative to the radar orientation (in math. direct rotation) */
      mask[n]=azimuth[k];
      }
    }
  status=save_SGXY(filename, grid, mask, spec, "relative_direction", "degree", "azimuth", 0);

  for(k=0;k<ny;k++) {
    for(l=0;l<nx;l++) {
      n=k*nx+l;
/*------------------------------------------------------------------------------
      +180.0 added to reverse direction (initially radard-ward) */
      double alpha=(azimuth[k]+azimuth_shift+180.);
/*------------------------------------------------------------------------------
      azimuth (relative to north, anti-direct rotation) */      
      alpha=90.-alpha+360;
      mask[n]=alpha;
      }
    }
  status=save_SGXY(filename, grid, mask, spec, "azimuth", "degree", "azimuth", 0);

  for(k=0;k<ny;k++) {
    for(l=0;l<nx;l++) {
      n=k*nx+l;
      mask[n]=radius[l];
      }
    }
  status=save_SGXY(filename, grid, mask, spec, "radius", "km", "radius", 0);

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int loadRadarRaw_netcdf(const char *filename, const char *name, double ref_lon, double ref_lat, double azimuth_shift, date_t start, date_t end, int nx)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

  Load RADAR-HF ascii data and rewrite them under netcdf format
  
-----------------------------------------------------------------------------*/
{
  int k, l, i, n, m, nitems, keep, status;
  int ny,nt;
  double *u, *mask, *count, spec=9999.;
  double *time, *lon, *lat,t0, *azimuth, *radius;
  grid_t grid;
  FILE *in;
  
  cout << endl << "-------- starting conversion --------" << endl << endl;
  
//  ny=mgr_line_count(filename);
  
  ny=71;
  
  double delta=ellapsed_time(start,end,'s');
  nt=delta/1200+1;
  
//  nt=70*24*3;

  lon=new double[nx*ny];
  lat=new double[nx*ny];
  time=new double[nt];
    
  u=new double[nx*ny*nt];

  for(n=0;n<nt*nx*ny;n++) {
    u[n]=spec;
    }
  
  mask=new double[nx*ny];
  
  count=new double[nx*ny];
  
  azimuth=new double[ny];
  radius=new double[nx];
 
  double x0,y0;
  
//  end  =date_t(2007,10,27,0.);
  
  t0=cnes_time(start,'s');
  
  string name_template="YYYYNNNHHMNSS_MMDD_VIG_TR_C2D_CRAD_FILT_"+(string)name+".nc";
  
  for(i=0;i<nt;i++) {
    time[i]=t0+i*1200.;
    date_t present = poctime_getdatecnes(time[i], 's');
    char *source, keyword[1024];
    int packing;
    status=decode_name(present, name_template.c_str(), 0, &source, &packing);
    time[i]/=24.*3600.;
    
    poc_data_t<double> lon,lat,azimuth,valids;
    poc_var_t vtime;
    poc_var_t var_u, var_v;
    poc_global_t info;
    poc_dim_t nx,ny,ncycles;
  
    status=poc_inq(source,&info);
    if(status!=0) NC_TRAP_ERROR(wexit,status,1,"poc_inq(\""+(string) source+"\",) error");
 
    status=lon.init(info,"lon");
    status=lat.init(info,"lat");
    status=azimuth.init(info,"azimuth");
    status=valids.init(info,"count");
    status=info.variables.find("time",  &vtime);
  
    poc_print(vtime);
      
    ncycles=vtime.dimensions[0];
  
    ny=lon.info.dimensions[0];
    nx=lon.info.dimensions[1];
 
/**----------------------------------------------------------------------------
    get longitudes */
    lon.read_data(source);
  
/**----------------------------------------------------------------------------
    get latitudes */
    lat.read_data(source);
    
/**----------------------------------------------------------------------------
    get latitudes */
    azimuth.read_data(source);
      
/**----------------------------------------------------------------------------
    get nvalids */
    valids.read_data(source);
  
    printf("\n\n");
    poc_print(valids.info);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int loadRadarRaw(const char *filename, date_t s, date_t e)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

/*------------------------------------------------------------------------------
//  Brezellec: -4.6667  48.069  -21.00 (Master) x=100.4556; y= 63.2617 km (Ref.: lon0=-6; lat0=47.5)
//  Garchine:  -4.7757  48.5030  20.00 (Slave)  x= 92.2432; y=111.5141 en km
//      cartes axe 0° de bas en haut  (n'est plus valable; 90 deg ajoutes)
//         x0=0; y0=0; theta0=0; lon0=-6; lat0=47.5;  % pt de reference pour la grille en km
//         if i_sta==1;
//           theta0=-21+90;
//           [x0 y0]=lonlat2km(lon0,lat0, -4.6667,  48.069);
//           teta=fliplr(teta')';disp('RETOURNEMENT DES AZIMUTHS POUR BREZELLEC (19/10/09)');
//         else
//           theta0=20+90;
//           [x0 y0]=lonlat2km(lon0,lat0, -4.7757,  48.503);
//         end
//         alphal=pi*(theta0+teta+90)/180.;
//         x=x0+dist*cos(alphal');
//         y=y0+dist*sin(alphal');
------------------------------------------------------------------------------*/
{
  int status;
  
  int nx[2]={75,77};
  
  date_t start[2], end;
  
//   double azimuth_shift[2]={-21.+90,+20.+90};
  double azimuth_shift[2]={-21.,+20.};
  double ref_lat[2],ref_lon[2];
  
  ref_lon[0]= -4.6667;
  ref_lat[0]= 48.069;
  
  ref_lon[1]= -4.7757;
  ref_lat[1]= 48.5030;

//   start[0]=date_t(2007, 8,23,600.);
//   start[1]=date_t(2007, 8,23,0.);
  start[0]=s;
  start[0].second=600.0;
  start[1]=s;
  
//  end  =date_t(2007,10,27,0.);
  end  =e;
  
#if 1
  const char *name[2]={"bre","gar"};
  status=loadRadarRaw_ascii("brezellec.nc",name[0], ref_lon[0], ref_lat[0], azimuth_shift[0], start[0], end, nx[0]);
  status=loadRadarRaw_ascii("garchine.nc", name[1], ref_lon[1], ref_lat[1], azimuth_shift[1], start[1], end, nx[1]);
#else
  const char *name[2]={"BRE","GAR"};
  status=loadRadarRaw_netcdf("brezellec.nc",name[0], ref_lon[0], ref_lat[0], azimuth_shift[0], start[0], end, nx[0]);
  status=loadRadarRaw_netcdf("garchine.nc", name[1], ref_lon[1], ref_lat[1], azimuth_shift[1], start[1], end, nx[1]);
#endif
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int set_mask(const char *a, const char *b, const char *c, grid_t & grid, float * & mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, nitems, status;
  FILE *in;
  
  grid.nx=83;
  grid.ny=72;

  grid.x=new double[grid.nx];
  grid.y=new double[grid.ny];
  
  grid.modeH=1;
  
  in=fopen(a,"r");
  for(n=0;n<grid.nx;n++) fscanf(in,"%lf", &grid.x[n]);
  fclose(in);
  
  in=fopen(b,"r");
  for(n=grid.ny-1; n>=0;n--) fscanf(in,"%lf", &grid.y[n]);
  fclose(in);
  
  grid.xmin=grid.x[0];
  grid.ymin=grid.y[0];

  grid.xmax=grid.x[grid.nx-1];
  grid.ymax=grid.y[grid.ny-1];

  mask=new float[grid.Hsize()];

  in=fopen(c,"r");
//   for(n=0;n<grid.Hsize();n++) {
//     nitems=fscanf(in,"%f", &(mask[n]));
//     }
  for(int j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      n=(grid.ny-j-1)*grid.nx+i;
      nitems=fscanf(in,"%f", &(mask[n]));
      }
    }
  fclose(in);
  status=save_SGXY("mask.nc", grid, mask, (float) -1.0, "mask", "none", "mask", 1);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int combine_tide(const char *a, const char *b, const char *c, const char *d, const char *filename, const char *gridfile, const char *varname, float *rmin)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, status;
  grid_t RadarGrid[2], grid;
  float *mask_value, flag;
  complex<float> *Unative[2], mask;
  complex<float> *U, *V;
  const char *source[2]={a,b};
  const char *series[2]={c,d};
  double *Anative[2], dmask;
  poc_global_t info;
  poc_data_t<double> azimuth[2];
  poc_data_t<float> reduction[2], radius[2];
  grid_t mask_grid;
  
  printf("reconstruct full 2D vectors from %s %s\n",a,b);
  
  mask=complex<float> (9999,9999);
  dmask=9999;

  for(int k=0;k<2;k++) {
    string name=(string) series[k] +"-analysis.nc";
    status=reduction[k].init(name,"reduction");
    if(status!=0) return(-1);
    status=reduction[k].read_data(name,0);
    if(status!=0) return(-1);
    }
  
  for(int k=0;k<2;k++) {
    string name=(string) series[k] +".nc";
    status=radius[k].init(name, "radius");
    if(status!=0) return(-1);
    status=radius[k].read_data(name,0);
    if(status!=0) return(-1);
    }
  
  status=set_mask("../gri_lon_deg.dat", "../gri_lat_deg.dat", "../mask11.dat", mask_grid, mask_value);
  
  status=poc_inq(a,&info);
  if(status!=0) return(-1);
 
//   status=poc_inq_var(a,"lon",   &vlon);
//   status=poc_inq_var(a,"lat",   &vlat);
  
  status=poc_get_grid(a,"Ua",&RadarGrid[0],1);
  status=poc_get_grid(b,"Ua",&RadarGrid[1],1);
  
  if(gridfile==0) {
    grid.xmin=min(RadarGrid[0].xmin,RadarGrid[1].xmin);
    grid.xmax=min(RadarGrid[0].xmax,RadarGrid[1].xmax);

    grid.ymin=min(RadarGrid[0].ymin,RadarGrid[1].ymin);
    grid.ymax=min(RadarGrid[0].ymax,RadarGrid[1].ymax);
  
    map_set2Dgrid(&grid, grid.xmin, grid.ymin, grid.xmax, grid.ymax,1./120.,1./120.);
    }
  else {
    status=poc_get_grid(gridfile,varname,&grid,1);
    }
    
  for(int k=0;k<2;k++) {
    poc_data_t<float> Ua, Ug;
    status=poc_inq_var(source[k],"Ua",  &Ua.info);
    Ua.init();
    Ua.read_data(source[k]);
    status=poc_inq_var(source[k],"Ug",  &Ug.info);
    Ug.init();
    Ug.read_data(source[k]);
    status=poc_inq_var(source[k],"azimuth",  &azimuth[k].info);
    azimuth[k].init();
    azimuth[k].read_data(source[k]);
    Unative[k]=new complex<float> [RadarGrid[k].Hsize()];
    Anative[k]=new double [RadarGrid[k].Hsize()];
    for(n=0;n<RadarGrid[k].Hsize();n++) {
//       if (reduction[k].data[n]<rmin[k]) {
//         Ua.data[n]=Ua.mask;
//         }
//       if (radius[k].data[n]>75.) {
//         Ua.data[n]=Ua.mask;
//         }
      status=map_interpolation(mask_grid, mask_value, (float) -1, RadarGrid[k].x[n], RadarGrid[k].y[n], &flag);
      if (NINT(flag)!=1.0) {
        Ua.data[n]=Ua.mask;
        }
      if(Ua.data[n]!=Ua.mask) {
        Unative[k][n]=Ua.data[n]*complex<float>(cos(-Ug.data[n]*M_PI/180.),sin(-Ug.data[n]*M_PI/180.));
        Anative[k][n]=azimuth[k].data[n]*M_PI/180.0;
        }
      else {
        Unative[k][n]=mask;
//        Anative[k][n]=mask;
        Anative[k][n]=azimuth[k].data[n]*M_PI/180.0;
        }
      }
    Ua.destroy_data();
    Ug.destroy_data();
    }


  U=new complex<float> [grid.Hsize()];
  V=new complex<float> [grid.Hsize()];

  for(int j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      complex<float> Uradial[2];
      double Aradial[2];
      n=grid.Hindex(i,j);
      double x,y;
      grid.xy(i,j,x,y);
      status=map_interpolation(RadarGrid[0], Unative[0], mask,  x, y, &Uradial[0]);
      status=map_interpolation(RadarGrid[0], Anative[0], dmask, x, y, &Aradial[0]);
      status=map_interpolation(RadarGrid[1], Unative[1], mask,  x, y, &Uradial[1]);
      status=map_interpolation(RadarGrid[1], Anative[1], dmask, x, y, &Aradial[1]);
      if( (Uradial[0]!=mask) && (Uradial[1]!=mask)) {
/*------------------------------------------------------------------------------
        reconstruction of the 2D vector*/
/*------------------------------------------------------------------------------
        eastward-northward direction = pi/2 - azimuth (north/east=90°/0° azimuth=0°/90° */
        vector2D_t u(M_PI/2.-Aradial[0]), v(M_PI/2.-Aradial[1]);
/*------------------------------------------------------------------------------
        apparently velocity was given toward the radar, now corrected */
//         cvector2D_t w=math_vector_coordinates03((complex<double>) -Uradial[0], u, (complex<double>) -Uradial[1], v);
        cvector2D_t w=math_vector_coordinates03((complex<double>) Uradial[0], u, (complex<double>) Uradial[1], v);
        U[n]=w.x;
        V[n]=w.y;
        }
      else {
        U[n]=mask;
        V[n]=mask;
        }
      }
    }
    
  pocgrd_t ncgrid;
    
  status= poc_createfile(filename);
  
  status=map_completegridaxis(&grid,2);
  
  status=poc_sphericalgrid_xy(filename,"",grid,&ncgrid);
  
  status=tides_savec1(filename, "Ua", "Ug", "currents", "m/s", grid, U, mask, ncgrid);
  status=tides_savec1(filename, "Va", "Vg", "currents", "m/s", grid, V, mask, ncgrid);

  return 0;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int combine_tide_new(const char *a, const char *b, const char *c, const char *d, 
                       const char *filename, const char *gridfile, const char *varname, float *rmin, float *max_distance, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, status;
  grid_t RadarGrid[2], grid;
  float *mask_value, flag;
  complex<float> *Unative[2], mask;
  complex<float> *U, *V;
  const char *source[2]={a,b};
  const char *series[2]={c,d};
  double *Anative[2], dmask;
  poc_global_t info;
  poc_data_t<double> azimuth[2];
  poc_data_t<float> reduction[2], radius[2];
  grid_t mask_grid;
  
  printf("reconstruct full 2D vectors from %s %s\n",a,b);
  
  mask=complex<float> (9999,9999);
  dmask=9999;
  
  for(int k=0;k<2;k++) {
    string name=(string) series[k];
    status=reduction[k].init(name,"reduction");
    if(status!=0) return(-1);
    status=reduction[k].read_data(name,0);
    if(status!=0) return(-1);
    }
  
  for(int k=0;k<2;k++) {
    string name=(string) series[k] +".nc";
    status=radius[k].init(name, "radius");
//     if(status!=0) return(-1);
    status=radius[k].read_data(name,0);
//     if(status!=0) return(-1);
    }
  
//   status=set_mask("../gri_lon_deg.dat", "../gri_lat_deg.dat", "../mask11.dat", mask_grid, mask_value);
  
  status=poc_inq(a,&info);
  if(status!=0) return(-1);
 
//   status=poc_inq_var(a,"lon",   &vlon);
//   status=poc_inq_var(a,"lat",   &vlat);
  
  status=poc_get_grid(a,"Ua",&RadarGrid[0],0);
  status=poc_get_grid(b,"Ua",&RadarGrid[1],0);
  
  if(gridfile==0) {
    grid.xmin=min(RadarGrid[0].xmin,RadarGrid[1].xmin);
    grid.xmax=min(RadarGrid[0].xmax,RadarGrid[1].xmax);

    grid.ymin=min(RadarGrid[0].ymin,RadarGrid[1].ymin);
    grid.ymax=min(RadarGrid[0].ymax,RadarGrid[1].ymax);
  
    map_set2Dgrid(&grid, grid.xmin, grid.ymin, grid.xmax, grid.ymax,1./120.,1./120.);
    }
  else {
    status=poc_get_grid(gridfile,varname,&grid,verbose);
    }
    
  for(int k=0;k<2;k++) {
    poc_data_t<float> Ua, Ug;
    status=poc_inq_var(source[k],"Ua",  &Ua.info);
    if(status!=0) TRAP_ERR_EXIT(status,"fatal error: could not init variable Ua from %s\n\n",source[k]);
    Ua.init();
    Ua.read_data(source[k]);
    status=poc_inq_var(source[k],"Ug",  &Ug.info);
    if(status!=0) TRAP_ERR_EXIT(status,"fatal error: could not init variable Ug from %s\n\n",source[k]);
    Ug.init();
    Ug.read_data(source[k]);
    status=poc_inq_var(source[k],"azimuth",  &azimuth[k].info);
    if(status!=0) TRAP_ERR_EXIT(status,"fatal error: could not init variable azimuth from %s\n\n",source[k]);
    azimuth[k].init();
    azimuth[k].read_data(source[k]);
    Unative[k]=new complex<float> [RadarGrid[k].Hsize()];
    Anative[k]=new double [RadarGrid[k].Hsize()];
    for(n=0;n<RadarGrid[k].Hsize();n++) {
      if (reduction[k].data[n]<rmin[k]) {
        Ua.data[n]=Ua.mask;
        }
      if (radius[k].data!=0)
	if (radius[k].data[n]>max_distance[k]) {
          Ua.data[n]=Ua.mask;
          }
//       status=map_interpolation(mask_grid, mask_value, (float) -1, RadarGrid[k].x[n], RadarGrid[k].y[n], &flag);
//       if (NINT(flag)!=1.0) {
//         Ua.data[n]=Ua.mask;
//         }
      if(Ua.data[n]!=Ua.mask) {
        Unative[k][n]=Ua.data[n]*complex<float>(cos(-Ug.data[n]*M_PI/180.),sin(-Ug.data[n]*M_PI/180.));
        Anative[k][n]=azimuth[k].data[n]*M_PI/180.0;
        }
      else {
        Unative[k][n]=mask;
//        Anative[k][n]=mask;
        Anative[k][n]=azimuth[k].data[n]*M_PI/180.0;
        }
      }
    Ua.destroy_data();
    Ug.destroy_data();
    }


  U=new complex<float> [grid.Hsize()];
  V=new complex<float> [grid.Hsize()];

  for(int j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      complex<float> Uradial[2];
      double Aradial[2];
      n=grid.Hindex(i,j);
      double x,y;
      grid.xy(i,j,x,y);
      status=map_interpolation(RadarGrid[0], Unative[0], mask,  x, y, &Uradial[0]);
      status=map_interpolation(RadarGrid[0], Anative[0], dmask, x, y, &Aradial[0]);
      status=map_interpolation(RadarGrid[1], Unative[1], mask,  x, y, &Uradial[1]);
      status=map_interpolation(RadarGrid[1], Anative[1], dmask, x, y, &Aradial[1]);
      if( (Uradial[0]!=mask) && (Uradial[1]!=mask)) {
/*------------------------------------------------------------------------------
        reconstruction of the 2D vector*/
/*------------------------------------------------------------------------------
        eastward-northward direction = pi/2 - azimuth (north/east=90°/0° azimuth=0°/90° */
        vector2D_t u(M_PI/2.-Aradial[0]), v(M_PI/2.-Aradial[1]);
        cvector2D_t w=math_vector_coordinates03((complex<double>) Uradial[0], u, (complex<double>) Uradial[1], v);
        U[n]=w.x;
        V[n]=w.y;
        }
      else {
        U[n]=mask;
        V[n]=mask;
        }
      }
    }
    
  pocgrd_t ncgrid;
    
  status= poc_createfile(filename);
  
  if(grid.modeH==0)  status=map_completegridaxis(&grid,1);
  
  status=poc_sphericalgrid_xy(filename,"",grid,&ncgrid);
  
  status=tides_savec1(filename, "Ua", "Ug", "currents", "m/s", grid, U, mask, ncgrid);
  status=tides_savec1(filename, "Va", "Vg", "currents", "m/s", grid, V, mask, ncgrid);

  return 0;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int combine_mean(const char *a, const char *b, const char *filename, const char *gridfile, const char *varname, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, status;
  grid_t RadarGrid[2], grid;
  float *Unative[2], mask;
  float *U, *V;
  const char *source[2]={a,b};
  double *Anative[2], dmask;
  poc_global_t info;
  poc_data_t<float> Ua[2];
  poc_data_t<double> azimuth[2];
  poc_data_t<double> lon,lat;
  
  mask=9999;
  dmask=NC_FILL_DOUBLE;
  
  status=poc_inq(a,&info);
  if(status!=0) return(-1);
 
  status=poc_get_grid(a,"mean_res",&RadarGrid[0],verbose);
  status=poc_get_grid(b,"mean_res",&RadarGrid[1],verbose);
  
  if(gridfile==0) {
    grid.xmin=min(RadarGrid[0].xmin,RadarGrid[1].xmin);
    grid.xmax=min(RadarGrid[0].xmax,RadarGrid[1].xmax);

    grid.ymin=min(RadarGrid[0].ymin,RadarGrid[1].ymin);
    grid.ymax=min(RadarGrid[0].ymax,RadarGrid[1].ymax);
  
    map_set2Dgrid(&grid, grid.xmin, grid.ymin, grid.xmax, grid.ymax,1./120.,1./120.);
    }
  else {
    status=poc_get_grid(gridfile,varname,&grid,verbose);
    }
    
  for(int k=0;k<2;k++) {
//     poc_data_t<float> Ua;
    status=poc_inq_var(source[k],"mean_res",  &Ua[k].info);
    if(status!=0) TRAP_ERR_EXIT(status,"fatal error: could not init variable mean_res from %s\n\n",source[k]);
    Ua[k].init();
    Ua[k].read_data(source[k]);
  mask=Ua[0].mask;
    status=poc_inq_var(source[k],"azimuth",  &azimuth[k].info);    
    if(status!=0) TRAP_ERR_EXIT(status,"fatal error: could not init variable azimuth from %s\n\n",source[k]);
    azimuth[k].init();
    azimuth[k].read_data(source[k]);
    Unative[k]=new float  [RadarGrid[k].Hsize()];
    Anative[k]=new double [RadarGrid[k].Hsize()];
    for(n=0;n<RadarGrid[k].Hsize();n++) {
      if(Ua[k].data[n]!=Ua[k].mask) {
        Unative[k][n]=Ua[k].data[n];
        Anative[k][n]=azimuth[k].data[n]*M_PI/180.0;
        }
      else {
        Unative[k][n]=mask;
//        Anative[k][n]=mask;
        Anative[k][n]=azimuth[k].data[n]*M_PI/180.0;
        }
      }
    Ua[k].destroy_data();
    }

  U=new float [grid.Hsize()];
  V=new float [grid.Hsize()];

  
  for(int j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      float  Uradial[2];
      double Aradial[2];
      n=grid.Hindex(i,j);
      double x,y;
      grid.xy(i,j,x,y);
      status=map_interpolation(RadarGrid[0], Unative[0], mask,  x, y, &Uradial[0]);
      status=map_interpolation(RadarGrid[0], Anative[0], dmask, x, y, &Aradial[0]);
      status=map_interpolation(RadarGrid[1], Unative[1], mask,  x, y, &Uradial[1]);
      status=map_interpolation(RadarGrid[1], Anative[1], dmask, x, y, &Aradial[1]);
      if( (Uradial[0]!=mask) && (Uradial[1]!=mask)) {
/*------------------------------------------------------------------------------
        reconstruction of the 2D vector*/
/*------------------------------------------------------------------------------
        eastward-northward direction = pi/2 - azimuth (north/east=90°/0° azimuth=0°/90° */
        vector2D_t u(M_PI/2.-Aradial[0]), v(M_PI/2.-Aradial[1]);
        vector2D_t w=math_vector_coordinates03((double) Uradial[0], u, (double) Uradial[1], v);
        U[n]=w.x;
        V[n]=w.y;
        }
      else {
        U[n]=mask;
        V[n]=mask;
        }
      }
    }
    
#if 1
  poc_data_t<float> vU=Ua[0], vV=Ua[0];
  vU.info.name="u_residual_average";
  vV.info.name="v_residual_average";
  vU.data=U;
  vV.data=V;
  status=vU.write_data(filename);
  status=vV.write_data(filename);
#else          
  status=save_SGXY(filename, grid, U, mask, "u_residual_average", "m/s", "u", 0);
  status=save_SGXY(filename, grid, V, mask, "v_residual_average", "m/s", "v", 0);
#endif    

  return 0;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int combine_residuals(const char *a, const char *b, const char *filename, const char *gridfile, const char *varname, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, status;
  grid_t RadarGrid[2], grid;
  float *Unative[2], mask;
  float *U, *V;
  const char *source[2]={a,b};
  double *Anative[2], dmask;
  poc_global_t info[2];
  poc_data_t<double> azimuth[2];
  poc_data_t<float>  Ua[2];
  poc_data_t<double> time[2];
  size_t nframes;
  int createfile, creategrid, createvariable, range[2];

  poc_data_t<double> lon,lat;
  
  mask=9999;
  dmask=NC_FILL_DOUBLE;
  
  status=poc_inq(a,&info[0]);
  if(status!=0) return(-1);
  status=poc_inq(b,&info[1]);
  if(status!=0) return(-1);
 
  status=poc_get_grid(a,"residuals",&RadarGrid[0],verbose);
  status=poc_get_grid(b,"residuals",&RadarGrid[1],verbose);
  
//   status=info[0].variables.find("time",  &vtime[0]);
//   status=info[1].variables.find("time",  &vtime[1]);
 
  if(gridfile==0) {
    grid.xmin=min(RadarGrid[0].xmin,RadarGrid[1].xmin);
    grid.xmax=min(RadarGrid[0].xmax,RadarGrid[1].xmax);

    grid.ymin=min(RadarGrid[0].ymin,RadarGrid[1].ymin);
    grid.ymax=min(RadarGrid[0].ymax,RadarGrid[1].ymax);
  
    map_set2Dgrid(&grid, grid.xmin, grid.ymin, grid.xmax, grid.ymax,1./120.,1./120.);
    }
  else {
    status=poc_get_grid(gridfile,varname,&grid,verbose);
    }
    
  for(int k=0;k<2;k++) {
    status=poc_inq_var(source[k],"residuals",  &Ua[k].info);
    Ua[k].init();
    status=poc_inq_var(source[k],"azimuth",  &azimuth[k].info);
    azimuth[k].init();
    status=poc_inq_var(source[k],"time",  &time[k].info);
    time[k].init();
    Unative[k]=new float  [RadarGrid[k].Hsize()];
    Anative[k]=new double [RadarGrid[k].Hsize()];
    azimuth[k].read_data(source[k]);
    nframes=Ua[k].nframes;
    time[k].read_data(source[k]);
    }
    
  U=new float [grid.Hsize()];
  V=new float [grid.Hsize()];
    
  mask=Ua[0].mask;
    
#if 1
  status=poc_createfile(filename);
  
  status=lon.init(info[0],"longitude");
  status=lat.init(info[0],"latitude");

/**----------------------------------------------------------------------------
  get longitudes */
  lon.read_data(source[0]);
  lon.write_data(filename);
  
/**----------------------------------------------------------------------------
  get latitudes */
  lat.read_data(source[0]);
  lat.write_data(filename);

  poc_data_t<float> vU=Ua[0], vV=Ua[0];
  vU.info.name="u_residual";
  vV.info.name="v_residual";
  vU.data=U;
  vV.data=V;
#endif    
    
  for(size_t m=0;m<nframes;m++) {
    printf("frame=%6d / %6d\r",m,nframes);
    for(int k=0;k<2;k++) {
      Ua[k].read_data(source[k],m);
      azimuth[k].read_data(source[k],m);
      time[k].read_data(source[k],m);
      for(n=0;n<RadarGrid[k].Hsize();n++) {
        if(Ua[k].data[n]!=Ua[k].mask) {
          Unative[k][n]=Ua[k].data[n];
          Anative[k][n]=azimuth[k].data[n]*M_PI/180.0;
          }
        else {
          Unative[k][n]=mask;
          Anative[k][n]=azimuth[k].data[n]*M_PI/180.0;
          }
        }
//       Ua[k].destroy_data();
//       azimuth[k].destroy_data();
      }

/*------------------------------------------------------------------------------
    reconstruction of the 2D vector*/
    for(int j=0;j<grid.ny;j++) {
      for(int i=0;i<grid.nx;i++) {
        float Uradial[2];
        double Aradial[2];
        n=grid.Hindex(i,j);
        double x,y;
        grid.xy(i,j,x,y);
        status=map_interpolation(RadarGrid[0], Unative[0], mask,  x, y, &Uradial[0]);
        status=map_interpolation(RadarGrid[0], Anative[0], dmask, x, y, &Aradial[0]);
        status=map_interpolation(RadarGrid[1], Unative[1], mask,  x, y, &Uradial[1]);
        status=map_interpolation(RadarGrid[1], Anative[1], dmask, x, y, &Aradial[1]);
        if( (Uradial[0]!=mask) && (Uradial[1]!=mask)) {
/*------------------------------------------------------------------------------
          eastward-northward direction = pi/2 - azimuth (north/east=90°/0° azimuth=0°/90° */
          vector2D_t u(M_PI/2.-Aradial[0]), v(M_PI/2.-Aradial[1]);
          vector2D_t w=math_vector_coordinates03((double) Uradial[0], u, (double) Uradial[1], v);
          U[n]=w.x;
          V[n]=w.y;
          }
        else {
          U[n]=mask;
          V[n]=mask;
          }
        }
      }
#if 1
    status=time[0].write_data(filename,m,1);
    status=vU.write_data(filename,m,1);
    status=vV.write_data(filename,m,1);
#else
    createfile=(m==0);
    creategrid=(m==0);
    createvariable=(m==0);
    status=save_SGXYT(filename, grid, U, mask, "u_residual", "m/s", "u", createfile, creategrid, createvariable, m, 0);
    createfile=0;
    creategrid=0;
    createvariable=(m==0);
    status=save_SGXYT(filename, grid, V, mask, "v_residual", "m/s", "v", createfile, creategrid, createvariable, m, 0);
#endif
    }
          

  return 0;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int concatenator(vector<string> & list, const char *outPath)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,status;
  string cmd;
  
  for(int k=0;k<list.size(); k++) cmd+=list[k]+" ";
  
  bool do_help=false;
    
/*------------------------------------------------------------------------------
  check arguments */
  
  if(list.size()<1){
    return 0;
    }
  
  if(access(outPath,F_OK)==0){
    fprintf(stderr,"*** %s already exists ***\n",outPath);
    do_help=true;
    }
  
  for(i=1;i<list.size();i++){
    if(access(list[i].c_str(),R_OK)==0) continue;
    fprintf(stderr,"*** can not read "+list[i]+" : %s ***\n",strerror(errno));
    do_help=true;
    }
  
  if(do_help){
    return -1;
    }
  
/*------------------------------------------------------------------------------
  header */
  
  poc_global_t gin,gout("constructed around " __LINE_FILE_PACKAGE_REVISION);
  bool withTimeDim;
  int maxLen;
  
  for(i=1;i<list.size();i++){
    status=poc_inq(list[i],&gin);
    NC_CHKERR_BASE_LINE(status,"poc_inq(\""+list[i]+"\",) error");
    
    for(j=0;j<gin.variables.size();j++){
      poc_var_t *var=&gin.variables[j];
      
      size_t len;
      status=poc_get_var_length(*var,&len);
      updatemax(&maxLen,len);
      
      withTimeDim=unlimitTimeDim(var);
      gout.variables<<*var;
      }
    
    for(j=0;j<gin.attributes.size();j++){
      poc_att_t *att=&gin.attributes[j];
      gout.attributes<<*att;
      }
    
    }
  
  gout<<poc_att_t("history",cmd);
  printf("creating %s ...",outPath);fflush(stdout);
  status=poc_create(outPath,gout,1);
  if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,1,"poc_create(\"%s\",,) error");
  printf(" done (max length=%d).\n",maxLen);
  
/*------------------------------------------------------------------------------
  data */
  
  /* frame indexes */
  int ofi=0,ifi;
  double *data=new double[maxLen];
  
  for(i=1;i<list.size();i++){
    
    printf("processing "+list[i]);fflush(stdout);
    
    status=poc_inq(list[i],&gin);
    NC_CHKERR_BASE_LINE(status,"poc_inq(\""+list[i]+"\",) error");
    
    /* number of input frames */
    int nif;
    poc_dim_t *dim;
    dim=findTimeDim(gin);
    if(dim!=0){
      nif=dim->len;
      printf(" (%d frames)",nif);fflush(stdout);
      }
    else
      nif=1;
    
    for(j=0;j<gin.variables.size();j++){
      poc_var_t *var=&gin.variables[j];
      printf("%s"+var->name,j==0?": ":", ");fflush(stdout);
      dim=findTimeDim(*var);
      for(ifi=0;ifi<max(1,nif);ifi++){
        poc_get_vara(list[i],*var,ifi    ,data);
        poc_put_vara(outPath,*var,ifi+ofi,data);
        if(dim==0)break;
        }
      }
    
    if(dim!=0)
      ofi+=nif;
    
    printf(", done.\n");
    }
  
  return 0;
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
    "USE\n"
    "  %s [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Detides radar data. POORLY CODED! DO NOT USE!\n"
    "  It takes garchine*.nc and brezellec*.nc as input.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  Show this help and exit.\n"
    "  --analyse=no  do not take garchine.nc and brezellec.nc as input\n"
    "  --reformat=yes  read YYYYNNNHHMN_.coc to produce garchine.nc and brezellec.nc\n"
    "  --spectrum  followed by spectrum, that may be a file\n"
    "  -g   followed by the path of the grid file. This is only necessary when the coordinates are not available in the files to analyse and you want to produce atlases or use control points.\n"
    "  -s   followed by the start date in dd/mm/yyyy format\n"
    "  -f   followed by the end date in dd/mm/yyyy format\n"
    );
  print_OPENMP_help(prog_name); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// function that will be ran when the executable is started
/** See the print_help <a href=#func-members>function</a> for its use. */
/*----------------------------------------------------------------------------*/
{
  int n,status;//argument index or number of files, NetCDF status
  char *gridfile=NULL, *keyword, *spectrum=NULL, *dum, *vflag=0;;
  date_t start=NADate,final=NADate;
  int analyse=1, reformat=0, long_periods=0;
  int verbose=0;
  vector<string> convention,input,residuals,rootname;
  vector<string> list[2],filelist;
  int ninput=0;
  date_t test(1950,1,1);
  float min_reduction[2];
  float max_distance[2];
  bool use_flags=false;
  
  min_reduction[0]=min_reduction[1]=0;
  max_distance[0]=max_distance[1]=1.e+10;
  
  struct timeval mainbefore;
  
   #if _OPENMP < 200505  // Version 2.5 May 2005
  #error Check whether your version _OPENMP (below) of openmp can do nesting properly. If so, update the above line accordingly.
  #pragma message "_OPENMP=" TOSTRING(_OPENMP)
  #endif

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"--help")==0) {
      print_help(argv[0]);
      wexit(0);
      }
    else if(strcmp(keyword,"-flag")==0) {
      use_flags=true;
      vflag=strdup(argv[n+1]);
      n++;
      n++;
      continue;
      }
    else if(strcmp(keyword,"--analysis=no")==0) {
      analyse=0;
      n++;
      continue;
      }
    else if(strcmp(keyword,"--analysis=yes")==0) {
      analyse=1;
      n++;
      continue;
      }
    else if(strcmp(keyword,"--analysis=no")==0) {
      analyse=0;
      n++;
      continue;
      }
    else if(strcmp(keyword,"--verbose")==0) {
      verbose=1;
      n++;
      continue;
      }
    else if(strcmp(keyword,"--reformat=yes")==0) {
      reformat=1;
      n++;
      continue;
      }
    else if(strcmp(keyword,"--reformat=no")==0) {
      reformat=0;
      n++;
      continue;
      }
    else if(strcmp(keyword,"--spectrum")==0) {
      spectrum=strdup(argv[n+1]);
      n++;
      n++;
      continue;
      }
    else if(strcmp(keyword,"--min_reduction")==0) {
      {
      string tmp = argv[n+1];
      vector<double> values=get_tokens(tmp," ");
      tmp.clear();
      min_reduction[0]=values[0];
      min_reduction[1]=values[1];
      }
      n++;
      n++;
      continue;
      }
    else if(strcmp(keyword,"--max_distance")==0) {
      {
      string tmp = argv[n+1];
      vector<double> values=get_tokens(tmp," ");
      tmp.clear();
      max_distance[0]=values[0];
      max_distance[1]=values[1];
      }
      n++;
      n++;
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 'c' :
          {
          string tmp = argv[n+1];
          convention=string_split(tmp," ");
	  tmp.clear();
          }
          ninput++;
          n++;
          n++;
          break;

        case 'i' :
          {
          string tmp = argv[n+1];
          input=string_split(tmp," ");
	  tmp.clear();
          }
          ninput++;
          n++;
          n++;
          break;

        case 'r' :
          {
          string tmp = argv[n+1];
          rootname=string_split(tmp," ");
	  tmp.clear();
          }
          ninput++;
          n++;
          n++;
          break;

        case 'g' :
          gridfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          dum= strdup(argv[n+1]);
          n++;
          n++;
          status=sscanf(dum,"%d/%d/%d",&start.day,&start.month,&start.year);
          if(status==3)
            start.second=0.;
          free(dum);
          break;

       case 'f' :
          dum= strdup(argv[n+1]);
          n++;
          n++;
          status=sscanf(dum,"%d/%d/%d",&final.day,&final.month,&final.year);
          if(status==3)
            final.second=0.;
          free(dum);
          break;

        case 'h' :
          print_help(argv[0]);
          wexit(0);
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

  gettimeofday(&mainbefore);
  
  if(rootname.size() ==0) {
    rootname.push_back("radar#1");
    rootname.push_back("radar#2");
    }

  if(spectrum==0) spectrum=strdup("DEEP");
  bool Z0=false;
  spectrum_t s=spectrum_init_ref(spectrum,Z0);
  
  if(!long_periods) {
    for(int k=0;k<s.n;k++) {
      tidal_wave w=s.waves[k];
      if(w.nT==0) {
        s.remove(w);
        k--;
        }
      }
    }
      
  if(analyse) {
#if 0
    if(reformat) {
/*------------------------------------------------------------------------------
      original data were given in ascii format, unefficient, convert them in 
      netcdf format; on output, still distance/azimuth data*/
      status=loadRadarRaw("",start, final);
      }
/*------------------------------------------------------------------------------
    radar data (distance/azimuth radial currents) harmonic analysis; on output,
    distance/azimuth harmonic constants */
    status=detide_radar("garchine", "garchine.nc", s);
    status=detide_radar("brezellec","brezellec.nc",s);
#else
    for(int k=0; k<convention.size();k++) {
      if(reformat) {
        int packing;
        char *filename, *varname=0;
        status=decode_name(test, convention[k].c_str(), varname, &filename, &packing);
        test=start;
        while(test<final) {
          status=decode_name(test, convention[k].c_str(), varname, &filename, &packing);
          FILE *in;
          in=fopen(filename,"r");
          if(in!=0) {
            list[k].push_back(filename);
            fclose(in);
            }
          delete[] filename;
          date_incremental(test,packing);
          }
        char *format="YYYY.MM.DD";
        char *s1, *s2;
        status=decode_name(start, format, varname, &s1, &packing);
        status=decode_name(final, format, varname, &s2, &packing);
        string tmp=rootname[k]+"-"+s1+"-"+s2+".nc";
        input.push_back(tmp);
        concatenator(list[k],tmp.c_str());
        list[k].clear();
        delete[] s1;
        delete[] s2;
        }
      int nvalidmin=15*24*6;
      status=detide_radar_new(rootname[k].c_str(), input[k].c_str(), s, nvalidmin, vflag, verbose);
      }
#endif
    }
    
    
//   if(filelist!=0) {
//     list=load_filelist(filelist);
//     }
  if(input.size()==0) {
    if(convention.size()==2) {
      for(int k=0; k<convention.size();k++) {
        char *format="YYYY.MM.DD";
        char *s1, *s2;
        int packing;
        status=decode_name(start, format, 0, &s1, &packing);
        status=decode_name(final, format, 0, &s2, &packing);
        string tmp=rootname[k]+"-"+s1+"-"+s2+".nc";
        input.push_back(tmp);
        delete[] s1;
        delete[] s2;
        }
      }
    }
    
    
//   if(analyse) {
// #if 0
// /*------------------------------------------------------------------------------
//     radar data (distance/azimuth radial currents) harmonic analysis; on output,
//     distance/azimuth harmonic constants */
//     status=detide_radar("garchine", "garchine.nc", s);
//     status=detide_radar("brezellec","brezellec.nc",s);
// #else
//     status=detide_radar_new("garchine", "garchine.nc", s, verbose);
//     status=detide_radar_new("brezellec","brezellec.nc",s, verbose);
// #endif
//     }

    residuals.push_back(rootname[0]+"-analysis.nc");
    residuals.push_back(rootname[1]+"-analysis.nc");

/*------------------------------------------------------------------------------
  create 2D velocity maps */
  for(int k=0;k<s.n;k++) {
//     string a="brezellec-"+(string)s.waves[k].name+(string)".nc";
//     string b="garchine-"+(string)s.waves[k].name+(string)".nc";
    string a=rootname[0]+"-"+(string)s.waves[k].name+(string)".nc";
    string b=rootname[1]+"-"+(string)s.waves[k].name+(string)".nc";
#if 0
//    string c="radar-"+(string)s.waves[k].name+(string)".nc";
//    status=combine_tide(a.c_str(),b.c_str(),c.c_str(),gridfile,"Ua");
    string e="radar-masked-"+(string)s.waves[k].name+(string)".nc";
    status=combine_tide(a.c_str(),b.c_str(), residuals[0].c_str(), residuals[1].c_str(), e.c_str(),gridfile,"Ua",rmin);
#else
    string e="radar-"+(string)s.waves[k].name+(string)".nc";
    status=combine_tide_new(a.c_str(),b.c_str(), residuals[0].c_str(), residuals[1].c_str(), e.c_str(),a.c_str(),"Ua",min_reduction,max_distance,verbose);
#endif
    if(status!=0) {
      printf("combine_tide error, wave=%s\n",s.waves[k].name);
      }
    }
    
//   string a="brezellec-analysis.nc";
//   string b="garchine-analysis.nc";
    
  string a=rootname[0]+"-analysis.nc";
  string b=rootname[1]+"-analysis.nc";
  
//   gridfile=strdup("brezellec.nc");
  
  gridfile=strdup(input[0].c_str());
  
  string c="radar-detided.nc";
  status=combine_residuals(a.c_str(),b.c_str(),c.c_str(),gridfile,"hcdt",verbose);
  
  c="radar-detided.nc";
  status=combine_mean(a.c_str(),b.c_str(),c.c_str(),gridfile,"hcdt",verbose);
  
  STDOUT_BASE_LINE(" ^^^^^^^^^^^^^ end of computation, that took %#g s with rt=%#g,rt2=%#g,rrt=%g,wt=%#g,hst=%#g,hat=%#g,detide=%#g ^^^^^^^^^^^^^\n",difftime(mainbefore),rt,rt2,rrt,wt,hst,hat,hct);
  exit(0);
}
