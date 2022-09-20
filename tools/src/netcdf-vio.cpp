
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief netcdf data input/output
*/
/*----------------------------------------------------------------------------*/

#include "poc-netcdf.hpp"


/*
#########
# INPUT #
#########
*/


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_vara(int ncid, int id, const size_t start[], const size_t count[],char *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_get_vara_text(ncid, id, start, count, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_vara(int ncid, int id, const size_t start[], const size_t count[],signed char *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_get_vara_schar(ncid, id, start, count, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_vara(int ncid, int id, const size_t start[], const size_t count[],short *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_get_vara_short(ncid, id, start, count, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_vara(int ncid, int id, const size_t start[], const size_t count[],int *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_get_vara_int(ncid, id, start, count, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_vara(int ncid, int id, const size_t start[], const size_t count[],float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_get_vara_float(ncid, id, start, count, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_vara(int ncid, int id, const size_t start[], const size_t count[],double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_get_vara_double(ncid, id, start, count, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int poc_get_vara_template(int ncid,const cdfvar_t &info,int frame, T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_get_vara(int,int,const size_t[],const size_t[],T*)
/*----------------------------------------------------------------------------*/
{
  int status,iDI;//,index of info.dim
  const cdfdim_t *dim;
  size_t *start,*count;
  
  start=new size_t[info.ndim];
  count=new size_t[info.ndim];
  
  for(iDI=0;iDI<info.ndim;iDI++){
    dim=&info.dim[iDI];
    
    if(isT(*dim)){
      while(frame<0){
        frame+=dim->length;
        }
      start[iDI]=frame;
      count[iDI]=1;
      }
    else{
      start[iDI]=0;
      count[iDI]=dim->length;
      }
    }
  
  status=poc_get_vara(ncid,info.id,start,count,z);
  
  delete[]start;
  delete[]count;
  
  return status;
}
#if 0
extern int poc_get_vara(int ncid,const cdfvar_t &info,int frame, T *z);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(int ncid,const cdfvar_t &info,int frame, T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(ncid,info,frame,z);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int poc_get_vara_template(const string &filename, const cdfvar_t & info,int frame, T *z, int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_get_vara(int,V,size_t,T*)
/*----------------------------------------------------------------------------*/
{
  int ncid,status;
  
  status=nc_open(filename.c_str(),NC_NOWRITE,&ncid);
  if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,verbose,"nc_open(\""+filename+"\",NC_NOWRITE,) error");
  
  status=poc_get_vara(ncid, info, frame, z);
  if(status!=NC_NOERR && verbose) NC_CHKERR_LINE_FILE(status,"poc_get_vara() error in "+filename);
  
  nc_close(ncid);
  
  return(status);
}
#if 0
extern int poc_get_vara(const string &filename, const cdfvar_t &info,int frame, T *z, int verbose=0);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(const string &filename, const cdfvar_t &info,int frame, T *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(filename,info,frame,z,verbose);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_var(int ncid, int id ,char *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_get_var_text(ncid, id, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_var(int ncid, int id ,int8_t *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_get_var_schar(ncid, id, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_var(int ncid, int id ,short *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=nc_get_var_short(ncid, id, z);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_var(int ncid, int id ,float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_get_var_float(ncid, id, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_var(int ncid, int id ,double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_get_var_double(ncid, id, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int poc_get_var_template(int ncid,const string & var, T *z, int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_get_var(int,int,T*)
/*----------------------------------------------------------------------------*/
{
  int status,varid;
  
  status=nc_inq_varid(ncid,var.c_str(),&varid);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"nc_inq_varid(,%s,) error",var.c_str());
  
  status=poc_get_var(ncid,varid,z);
  
  return status;
}
#if 0
extern int poc_get_var(int ncid,const string & var, T *z, int verbose=0);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_var(int ncid,const string & var, T *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_var_template(ncid,var,z,verbose);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_var(int ncid,const string & var, double *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_var_template(ncid,var,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_var(int ncid,const string & var, float *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_var_template(ncid,var,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_var(int ncid,const string & var, short *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_var_template(ncid,var,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int poc_get_var_template(const string & filename,const string & var, T *z, int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_get_var(int,string,T*,int)
/*----------------------------------------------------------------------------*/
{
  int ncid,status;
  
  status=nc_open(filename.c_str(),NC_NOWRITE,&ncid);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"nc_open(\""+filename+"\",NC_NOWRITE,) error");
  
  status=poc_get_var(ncid, var, z);
  if(status!=NC_NOERR && verbose)NC_CHKERR_LINE_FILE(status,"poc_get_var() error in "+filename);
  
  nc_close(ncid);
  
  return(status);
}
#if 0
extern int poc_get_var(const string & filename,const string & var, T *z, int verbose=0);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_var(const string & filename,const string & var, T *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_var_template(filename,var,z,verbose);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_var(const string & filename,const string & var, double *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_var_template(filename,var,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_var(const string & filename,const string & var, float *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_var_template(filename,var,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_var(const string & filename,const string & var, short *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_var_template(filename,var,z,verbose);
  return status;
}


/*
##########
# OUTPUT #
##########
*/


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_put_vara(int ncid, int varid, const size_t *start, const size_t *count,const char *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_put_vara_text(ncid,varid,start,count,out);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_put_vara(int ncid, int varid, const size_t *start, const size_t *count,const signed char *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_put_vara_schar(ncid,varid,start,count,out);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_put_vara(int ncid, int varid, const size_t *start, const size_t *count,const short *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_put_vara_short(ncid,varid,start,count,out);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_put_vara(int ncid, int varid, const size_t *start, const size_t *count,const int *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_put_vara_int(ncid,varid,start,count,out);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_put_vara(int ncid, int varid, const size_t *start, const size_t *count,const double *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_put_vara_double(ncid,varid,start,count,out);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_put_vara(int ncid, int varid, const size_t *start, const size_t *count,const float *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_put_vara_float(ncid,varid,start,count,out);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_put_var(int ncid, int varid, const char *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_put_var_text(ncid,varid,out);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_put_var(int ncid, int varid, const signed char *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_put_var_schar(ncid,varid,out);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_put_var(int ncid, int varid, const short *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_put_var_short(ncid,varid,out);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_put_var(int ncid, int varid, const int *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_put_var_int(ncid,varid,out);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_put_var(int ncid, int varid, const double *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_put_var_double(ncid,varid,out);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_put_var(int ncid, int varid, const float *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_put_var_float(ncid,varid,out);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_put_var(const char *filename, int varid, const float *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,result=0;
  int ncid;
  
  status=nc_open(filename,NC_WRITE,&ncid);
  nc_check_error(status,__LINE__,__FILE__);
  if(status !=0) return -1;
 
  status=nc_put_var_float(ncid,varid,out);
  nc_check_error(status,__LINE__,__FILE__,"poc_put_var error");
  if(status !=0) return -1;
 
  status = nc_close(ncid);
  nc_check_error(status,__LINE__,__FILE__);
  if(status !=0) return -1;
  
 return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_put_var(const char *filename, int varid, const short *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,result=0;
  int ncid;
  
  status=nc_open(filename,NC_WRITE,&ncid);
  nc_check_error(status,__LINE__,__FILE__);
  if(status !=0) return -1;
 
  status=nc_put_var_short(ncid,varid,out);
  nc_check_error(status,__LINE__,__FILE__,"poc_put_var error");
  if(status !=0) return -1;
 
  status = nc_close(ncid);
  nc_check_error(status,__LINE__,__FILE__);
  if(status !=0) return -1;
  
 return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_put_var(const char *filename, int varid, const double *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,result=0;
  int ncid;
  
  status=nc_open(filename,NC_WRITE,&ncid);
  nc_check_error(status,__LINE__,__FILE__);
  if(status !=0) return -1;
 
  status=nc_put_var_double(ncid,varid,out);
  nc_check_error(status,__LINE__,__FILE__,"poc_put_var error");
  if(status !=0) return -1;
 
  status = nc_close(ncid);
  nc_check_error(status,__LINE__,__FILE__);
  if(status !=0) return -1;
  
 return result;
}
