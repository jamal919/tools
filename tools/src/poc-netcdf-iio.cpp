
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

\brief info input/output for the new (as of 2012-01-24) poc-netcdf classes
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#if TUGO

#include "tugo-prototypes.h"

#endif

#include "poc-netcdf-data.hpp" /* for poc_get_axes() */

#include "poc-grib.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_inq_dim(int ncid,int dimid, poc_dim_t *dim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Get dimension
/**
\param ncid NetCDF file ID
\param dimid NetCDF dimension ID
\param *att pointer to the initialised dimension
\returns NC_NOERR if success or the NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  int status,unlimdimid;
  char name[NC_MAX_NAME+1];
  size_t len;
  
  status=nc_inq_dim(ncid,dimid,name,&len);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"nc_inq_dim(,%d,) error",dimid);
  
  status=nc_inq_unlimdim(ncid,&unlimdimid);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"nc_inq_unlimdim() error");
  
  dim->init(name,len,dimid==unlimdimid);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_inq_dim(const string &path,int dimid, poc_dim_t *dim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
{
  int status,close_status,ncid;
  int verbose=0;
  
  status=nc_open(path.c_str(),NC_NOWRITE,&ncid);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"nc_open(\""+path+"\",NC_NOWRITE,) error");
  
  status=poc_inq_dim(ncid,dimid,dim);
  if(status!=NC_NOERR){
    if(verbose)NC_CHKERR_BASE_LINE(status,"poc_inq_dim error");
    }
  
  close_status=nc_close(ncid);
  NC_CHKERR_BASE_LINE(close_status,"nc_close() error with "+path);
  
  if(status==NC_NOERR)
    status=close_status;

  return status;
}


//from #include <sys/types.h>
inline int nc_get_att_(int ncid, int varid, const char name[], char in[])
{return nc_get_att_text(ncid,varid,name,in);}
inline int nc_get_att_(int ncid, int varid, const char name[], unsigned char in[])
{return nc_get_att_uchar(ncid,varid,name,in);}
inline int nc_get_att_(int ncid, int varid, const char name[], signed char in[])
{return nc_get_att_schar(ncid,varid,name,in);}
inline int nc_get_att_(int ncid, int varid, const char name[], short in[])
{return nc_get_att_short(ncid,varid,name,in);}
inline int nc_get_att_(int ncid, int varid, const char name[], int in[])
{return nc_get_att_int(ncid,varid,name,in);}
inline int nc_get_att_(int ncid, int varid, const char name[], int64_t in[])
{return nc_get_att_long(ncid,varid,name,(long*)in);}
inline int nc_get_att_(int ncid, int varid, const char name[], float in[])
{return nc_get_att_float(ncid,varid,name,in);}
inline int nc_get_att_(int ncid, int varid, const char name[], double in[])
{return nc_get_att_double(ncid,varid,name,in);}
#ifdef NC_USHORT
inline int nc_get_att_(int ncid, int varid, const char name[], unsigned short in[])
{return nc_get_att_ushort(ncid,varid,name,in);}
/* Contrary to what man:/netcdf says, nc_get_att_int64 and nc_get_att_uint64 do not exist !!! */
inline int nc_get_att_(int ncid, int varid, const char name[], unsigned int in[])
{return nc_get_att_uint(ncid,varid,name,in);}
//On 64 bits machines, int is 32 bits and long and long long are both 64 bits.
inline int nc_get_att_(int ncid, int varid, const char name[], uint64_t in[])
{return nc_get_att_ulonglong(ncid,varid,name,(unsigned long long*)in);}
inline int nc_get_att_(int ncid, int varid, const char name[], char *in[])
{return nc_get_att_string(ncid,varid,name,in);}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_inq_att(int ncid,int varid,const char*attname,poc_att_t *att,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Get variable or global attribute
/**
\param ncid NetCDF file ID
\param varid NetCDF variable ID. If NC_GLOBAL, get global attribute
\param attname attribute name
\param *att pointer to the initialised attribute
\returns NC_NOERR if success or the NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  int status;
  size_t len;
  nc_type type;
  
  status=nc_inq_att(ncid,varid,attname,&type,&len);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"nc_inq_att(,,\"%s\",,) error",attname);
  
  status=nc_inq_attid(ncid,varid,attname,&att->id);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"nc_inq_attid(,,\"%s\",) error",attname);
  
  switch(type){
    case NC_CHAR:            /* ISO/ASCII character */
      {
      char *data=new char[len+1];
      status=nc_get_att_(ncid,varid,attname,data);
      if(status!=NC_NOERR){NC_CHKERR_BASE_LINE(status,"nc_get_att_(,,\"%s\",) error",attname);delete[]data;return status;}
      data[len]=0;
      att->init(attname,data,len);
      }
      break;
  #define unionDataType(t) \
    {t *data=new t[len]; \
    status=nc_get_att_(ncid,varid,attname,data); \
    if(status!=NC_NOERR){NC_CHKERR_BASE_LINE(status,"nc_get_att_(,%d,\"%s\",) error",varid,attname);delete[]data;return status;} \
    att->init(attname,data,len);}break;
    case NC_BYTE:
      unionDataType(int8_t)   /* signed 1 byte integer */
    case NC_SHORT:
      unionDataType(int16_t)  /* signed 2 byte integer */
    case NC_INT:
      unionDataType(int32_t)  /* signed 4 byte integer */
    case NC_FLOAT:
      unionDataType(float)    /* single precision floating point number */
    case NC_DOUBLE:
      unionDataType(double)   /* double precision floating point number */
    #ifdef NC_UBYTE
    case NC_UBYTE:
      unionDataType(uint8_t)  /* unsigned 1 byte int */
    case NC_USHORT:
      unionDataType(uint16_t) /* unsigned 2-byte int */
    case NC_UINT:
      unionDataType(uint32_t) /* unsigned 4-byte int */
    case NC_INT64:
      unionDataType(int64_t)  /* signed 8-byte int */
    case NC_UINT64:
      unionDataType(uint64_t) /* unsigned 8-byte int */
    case NC_STRING:
      {char **data=new char*[len];
      status=nc_get_att_(ncid,varid,attname,data);
      if(status!=NC_NOERR){NC_CHKERR_BASE_LINE(status,"nc_get_att_(,%d,\"%s\",) error",varid,attname);deletep2D(&data,len,free);return status;}
      att->init(attname,data,len);}break;
    default:
      NC_CHKERR_BASE_LINE(NC_EBADTYPE,"unknown type %d",type);wexit(NC_EBADTYPE);
    #else
    default:
      NC_CHKERR_LINE_FILE(NC_EBADTYPE,"unknown type %d. You may need to UPGRADE YOUR NetCDF!",type);wexit(NC_EBADTYPE);
    #endif
    }
  #undef unionDataType
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_inq_att(int ncid,const char*attname,poc_att_t *att,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Get global attribute
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  return poc_inq_att(ncid,NC_GLOBAL,attname,att,verbose); /** \endcode  So see poc_inq_att(int,int,const char*,poc_att_t*)  */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_inq_att(const string &path,const char*attname,poc_att_t *att,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,close_status,ncid;
  
  status=nc_open(path.c_str(),NC_NOWRITE,&ncid);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"nc_open(\""+path+"\",NC_NOWRITE,) error");
  
  status=poc_inq_att(ncid,attname,att,verbose);
  if(status!=NC_NOERR){
    if(verbose>=0)NC_CHKERR_BASE_LINE(status,"poc_inq_att() error");
    }
  
  close_status=nc_close(ncid);
  NC_CHKERR_BASE_LINE(close_status,"nc_close() error with "+path);
  
  if(status==NC_NOERR)
    status=close_status;

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_inq_var(int ncid,int varid,poc_var_t *var)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Get variable informations
/**
\param ncid NetCDF file ID
\param varid NetCDF variable ID
\param *var if NULL, only check existence of variable
\returns NC_NOERR if success or the NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  int status,i;
  char name[NC_MAX_NAME+1];
  char attname[NC_MAX_NAME+1];
  nc_type type;
  int ndims,*dimids,natts;
  poc_dim_t dim;
  poc_att_t att;
  
  status=nc_inq_var(ncid,varid,name,&type,&ndims,NULL,&natts);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"nc_inq_var() error");
  if(var==0)return status;
  *var=poc_var_t(name,type,varid);
  
  ///-dimensions with poc_inq_dim()
  dimids=new int[ndims];
  status=nc_inq_vardimid(ncid,varid,dimids);
  if(status!=NC_NOERR){NC_CHKERR_BASE_LINE(status,"nc_inq_vardimid() error with %s",name);goto error;}
  for(i=0;i<ndims;i++){
    status=poc_inq_dim(ncid,dimids[i],&dim);
    if(status!=NC_NOERR){NC_CHKERR_BASE_LINE(status,"poc_inq_dim(,%d:%d/%d,) error with %s",dimids[i],i+1,ndims,name);goto error;}
    *var<<dim;
    }
  
  poc_get_axes("",var);
  
  ///-attributes with poc_inq_att()
  for(i=0;i<natts;i++){
    status=nc_inq_attname(ncid,varid,i,attname);
    if(status!=NC_NOERR) {
      NC_CHKERR_BASE_LINE(status,"nc_inq_attname(,,%d:<%d,) error with %s",i,natts,name);goto error;
      }
    status=poc_inq_att(ncid,varid,attname,&att);
    if(status!=NC_NOERR) {
      NC_CHKERR_BASE_LINE(status,"poc_inq_att(,,\"%s\":%d/%d,) error with %s",attname,i+1,natts,name);goto error;
      }
    *var<<att;
    }
  
  error:
  delete[]dimids;
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_inq_var(int ncid,const string &varname,poc_var_t *var,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Get variable informations
/**
\param ncid NetCDF file ID
\param varname variable name
\param *var
\param verbose
\returns NC_NOERR if success or the NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  int status,varid;
  
/*------------------------------------------------------------------------------
  get the variable ID */
  status=nc_inq_varid(ncid,varname.c_str(),&varid);
  NC_TRAP_ERROR(return,status,verbose,"nc_inq_varid(,\""+varname+"\",) error");
  
  ///-calls poc_inq_var(int,int,poc_var_t*,int)
  status=poc_inq_var(ncid,varid,var);
  if(status!=NC_NOERR)
    NC_CHKERR_LINE_FILE(status,"poc_inq_var() error: "+varname+" ncid=%d varid=%d", ncid,varid);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_inq_var(const string &path,const string &varname,poc_var_t *var,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Get variable informations
/**
\param path NetCDF file path
\param varname variable name
\param *var
\param verbose
\returns NC_NOERR if success or the NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  int status,close_status,ncid;
  
/*------------------------------------------------------------------------------
  open */
  status=nc_open(path.c_str(),NC_NOWRITE,&ncid);
  if(status!=NC_NOERR){
    
    if(status==NC_ENOTNC) do{
/*------------------------------------------------------------------------------
      fail safe */
      poc_global_t global;
      int i;
      
      status=poc_inq(path,&global,verbose);
      if(status!=NC_NOERR) break;
      
      i=global.variables.find(varname,var);
      if(i>=0)
        return status;
      
/*----------------------------------------------------------------------------*/
      }while(0);
    
    NC_TRAP_ERROR(return,status,verbose,"nc_open(\""+path+"\",NC_NOWRITE,) error");
    }
  
  ///-calls poc_inq_var(int,string,poc_var_t*,int)
  status=poc_inq_var(ncid,varname,var,verbose);
  if(status!=NC_NOERR){
    if(verbose>=0)NC_CHKERR_BASE_LINE(status,"poc_inq_var(,\""+path+"\",var=\""+varname+"\",,) error");
    }
  
  ///-closes the file
  close_status=nc_close(ncid);
  if(verbose>=0)NC_CHKERR_BASE_LINE(close_status,"nc_close() error with "+path);
  
  status=(status==NC_NOERR?close_status:status);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_inq(int ncid,poc_global_t *global, int *format)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// get header informations about a file
/**
\param ncid NetCDF file ID
\returns NC_NOERR if success or the NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  int status,i;
  char attname[NC_MAX_NAME+1];
  int ndims,nvars,natts;
  poc_dim_t dim;
  poc_var_t var;
  poc_att_t att;
  
  status=nc_inq(ncid,&ndims,&nvars,&natts,&i);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"nc_inq() error");
  
  global->dimensions.clear(ndims);
  global->attributes.clear(natts);
  global->variables.clear(nvars);
  
/*------------------------------------------------------------------------------
  dimensions */
  for(i=0;i<ndims;i++){
    status=poc_inq_dim(ncid,i,&dim);
    if(status!=NC_NOERR){
      /* THIS CAN HAPPEN AND NOT BE A BIG DEAL ! */
      NC_CHKERR_BASE_LINE(status,"poc_inq_dim(,%d:<%d,) error",i,ndims);
      continue;
      }
    *global<<dim;
    }
  
/*------------------------------------------------------------------------------
  attributes */
  for(i=0;i<natts;i++){
    status=nc_inq_attname(ncid,NC_GLOBAL,i,attname);
    if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"nc_inq_attname(,NC_GLOBAL,%d:<%d,) error",i,natts);
    status=poc_inq_att(ncid,attname,&att);
    if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"poc_inq_att(,\"%s\":%d/%d,) error",attname,i+1,natts);
    *global<<att;
    }
  
/*------------------------------------------------------------------------------
  variables */
  for(i=0;i<nvars;i++){
    status=poc_inq_var(ncid,i,&var);
    if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"poc_inq_var(,%d:<%d,) error",i,nvars);
    *global<<var;
    }
  
/*------------------------------------------------------------------------------
  format */
  if(format!=0){
    int status1;
    status1=nc_inq_format(ncid,format);
    NC_CHKERR_BASE_LINE(status1,"nc_inq_format() error");
    }
  
  return status;
}


extern int poc_qoddy_inq(const string &path,poc_global_t *global,int verbose);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_inq(const string &path,poc_global_t *global,int verbose, int *format)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// get header informations about a file
/**
\param path
\param global
\param verbose
\returns NC_NOERR if success or the NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  int status,close_status,ncid;
  
/*------------------------------------------------------------------------------
  open the file */
  status=nc_open(path.c_str(),NC_NOWRITE,&ncid);
  if(status!=NC_NOERR){
    
/*----------------------------------------------------------------------------*/
    if(is_grib(path)){
#ifdef HAVE_LIBGRIB_API
      status=poc_grib_inq(path,global,verbose);
      
      if(status!=0 and verbose>=0)
        STDERR_BASE_LINE("poc_grib_inq(\""+path+"\",) error (%d %s)\n",status,strerror(status));
      
      return status;
#else
      TRAP_ERR_RETURN(ENOEXEC,1,"Please compile with GRIB_API\n");
#endif
      }
    
/*----------------------------------------------------------------------------*/
    const string extension=get_extension(path);
    
    if( extension=="s2r" or
        extension=="s2c" or
        extension=="v2r" or
        extension=="v2c" ){
      status=poc_qoddy_inq(path,global,verbose);
      
      if(status!=0 and verbose>=0)
        STDERR_BASE_LINE("poc_qoddy_inq(\""+path+"\",) error (%d %s)\n",status,strerror(status));
      
      return status;
      }
    
/*----------------------------------------------------------------------------*/
    NC_TRAP_ERROR(return,status,verbose,"nc_open(\""+path+"\",NC_NOWRITE,) error");
    }

  ///-calls poc_inq(int,poc_global_t*)
  status=poc_inq(ncid,global,format);
  if(verbose>=0) NC_CHKERR_BASE_LINE(status,"poc_inq() error with "+path);
  
  ///-closes the file
  close_status=nc_close(ncid);
  if(verbose>=0) NC_CHKERR_BASE_LINE(close_status,"nc_close() error with "+path);
  
  if(status==NC_NOERR)
    status=close_status;

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_open_and_may_create(const string &path, int *ncid, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
{
  int status,i;
  
  for(i=0;i<2;i++){
    status=nc_open(path.c_str(),NC_WRITE,ncid);
    if(status!=ENOENT or i>0)
      break;
    poc_global_t global;
    global<<poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION);
    status=poc_create(path,global,verbose);
    NC_TRAP_ERROR(return,status,verbose,"poc_create(\""+path+"\",,,) error");
    }
  
  NC_TRAP_ERROR(return,status,verbose,"nc_open(\""+path+"\",NC_WRITE,) error");
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_def_dim(int ncid,const poc_dim_t & dim,int *dimid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Define a dimension
/**
\param ncid NetCDF file ID
\param dim dimension
\param *dimid optional pointer to the NetCDF ID of the new dimension
\returns NC_NOERR if success or the NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  int status;
  size_t len=dim.isunlimited?NC_UNLIMITED:dim.len;
  
  ///The file must be in define mode
  //nc_redef(ncid);
  
  status=nc_inq_dimid(ncid,dim.name.c_str(),dimid);
  switch(status){
  case NC_EBADDIM:
    status=nc_def_dim(ncid,dim.name.c_str(),len,dimid);//man:/netcdf ``If dimid is not a NULL pointer then upon successful completion dimid will contain the dimension ID of the newly created dimension.''
    if(status!=NC_NOERR)NC_CHKERR_BASE_LINE(status,"nc_def_dim() error");
    break;
  case NC_NOERR:{
    poc_dim_t dim0;
    status=poc_inq_dim(ncid,*dimid,&dim0);
    if(status!=NC_NOERR)NC_CHKERR_BASE_LINE(status,"poc_inq_dim() error");
    if(dim0!=dim){
      status=NC_ENAMEINUSE; /* String match to name in use */
      NC_CHKERR_BASE_LINE(status,"dimension \""+dim.name+"\" has been defined with %d!=%d",dim0.len,dim.len);
      }
    }break;
  default:
    NC_CHKERR_BASE_LINE(status,"nc_inq_dimid(,\""+dim.name+"\",) error");
    }
  
  return status;
}


inline int nc_put_att_(int ncid, int varid, const char name[], nc_type xtype, size_t len, const char out[])
{return nc_put_att_text(ncid,varid,name,len,out);}
inline int nc_put_att_(int ncid, int varid, const char name[], nc_type xtype, size_t len, const unsigned char out[])
{return nc_put_att_uchar(ncid,varid,name,xtype,len,out);}
inline int nc_put_att_(int ncid, int varid, const char name[], nc_type xtype, size_t len, const signed char out[])
{return nc_put_att_schar(ncid,varid,name,xtype,len,out);}
inline int nc_put_att_(int ncid, int varid, const char name[], nc_type xtype, size_t len, const short out[])
{return nc_put_att_short(ncid,varid,name,xtype,len,out);}
inline int nc_put_att_(int ncid, int varid, const char name[], nc_type xtype, size_t len, const int out[])
{return nc_put_att_int(ncid,varid,name,xtype,len,out);}
inline int nc_put_att_(int ncid, int varid, const char name[], nc_type xtype, size_t len, const int64_t out[])
{return nc_put_att_long(ncid,varid,name,xtype,len,(long*)out);}
inline int nc_put_att_(int ncid, int varid, const char name[], nc_type xtype, size_t len, const float out[])
{return nc_put_att_float(ncid,varid,name,xtype,len,out);}
inline int nc_put_att_(int ncid, int varid, const char name[], nc_type xtype, size_t len, const double out[])
{return nc_put_att_double(ncid,varid,name,xtype,len,out);}
#ifdef NC_USHORT
inline int nc_put_att_(int ncid, int varid, const char name[], nc_type xtype, size_t len, const unsigned short out[])
{return nc_put_att_ushort(ncid,varid,name,xtype,len,out);}
inline int nc_put_att_(int ncid, int varid, const char name[], nc_type xtype, size_t len, const unsigned int out[])
{return nc_put_att_uint(ncid,varid,name,xtype,len,out);}
/* Contrary to what man:/netcdf says, nc_put_att_int64 and nc_put_att_uint64 do not exist !!! */
inline int nc_put_att_(int ncid, int varid, const char name[], nc_type xtype, size_t len, const uint64_t out[])
{return nc_put_att_ulonglong(ncid,varid,name,xtype,len,(unsigned long long*)out);}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_def_att(int ncid,int varid,const poc_att_t & att)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Define a variable or a global attribute
/**
\param ncid NetCDF file ID
\param varid NetCDF variable ID. If NC_GLOBAL, define global attribute
\param att attribute
\returns NC_NOERR if success or the NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  int def_status,status;
  
  /* BECAUSE OF A COMPILATION BUG WITH poc_var_t::init(,,,,) WHEN USING -Ofast */
  static int verbose=1;
  if(att.name==_FillValue and att.type==NC_DOUBLE and isnan(*att.data_as_double)){
    const int oldVerbose=verbose;
    verbose=-1;
    TRAP_ERR_RETURN(0,oldVerbose,"FAIL SAFE: poc_att_t(\""+att.name+"\",%g)\n",*att.data_as_double);
    }
  
/*-----------------------------------------------------------------------------
  puts the file in define mode */
  def_status=nc_redef(ncid);
  if(def_status!=NC_NOERR && def_status!=NC_EINDEFINE)NC_TRAP_ERROR(return,def_status,1,"nc_redef() error");
  
  switch(att.type){
    #define unionDataType(t) status=nc_put_att_(ncid,varid,att.name.c_str(),att.type,att.len,att.data_as_ ## t);break
    case NC_CHAR:
      unionDataType(char);     /* signed 1 byte integer */
    case NC_BYTE:
      unionDataType(int8_t);   /* signed 1 byte integer */
    case NC_SHORT:
      unionDataType(int16_t);  /* signed 2 byte integer */
    case NC_INT:
      unionDataType(int32_t);  /* signed 4 byte integer */
    case NC_FLOAT:
      unionDataType(float);    /* single precision floating point number */
    case NC_DOUBLE:
      unionDataType(double);   /* double precision floating point
number */
    #ifdef NC_UBYTE
    case NC_UBYTE:
      unionDataType(uint8_t);  /* unsigned 1 byte int */
    case NC_USHORT:
      unionDataType(uint16_t); /* unsigned 2-byte int */
    case NC_UINT:
      unionDataType(uint32_t); /* unsigned 4-byte int */
    case NC_INT64:
      unionDataType(int64_t);  /* signed 8-byte int */
    case NC_UINT64:
      unionDataType(uint64_t); /* unsigned 8-byte int */
    default:
      NC_CHKERR_BASE_LINE(NC_EBADTYPE,"unknown type %d",att.type);exit(-2);
    #else
    default:
      NC_CHKERR_BASE_LINE(NC_EBADTYPE,"unknown type %d. You may need to UPGRADE YOUR NetCDF!",att.type);exit(-2);
    #warning UPGRADE YOUR NetCDF : NC_UBYTE NC_USHORT NC_UINT NC_INT64 NC_UINT64 NC_STRING are undefined
    #endif
    #undef unionDataType
    }
  if(status!=NC_NOERR)
    NC_CHKERR_BASE_LINE(status,"nc_put_att_(,,(\""+att.name+"\"),%d,%d,) error",att.type,att.len);
  if(status==NC_EBADTYPE and att.type>=NC_UBYTE)
    STDERR_BASE_LINE("Maybe this is because you are trying to put a NetCDF-4-typed attribute or variable in a NetCDF 3 file...\n");
  
/*-----------------------------------------------------------------------------
  takes the file out of define mode if it was not so already */
  if(def_status==NC_NOERR){
    status=nc_enddef(ncid);
    NC_CHKERR_BASE_LINE(status,"nc_enddef() error");
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_def_att(int ncid, const poc_att_t & att)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Put global attribute
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  return poc_def_att(ncid,NC_GLOBAL,att); /** \endcode  So see poc_def_att(int,int,const poc_att_t&)  */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_def_att(const string &path, const poc_att_t &att,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Put global attribute
/*----------------------------------------------------------------------------*/
{
  int status,end_status,close_status,ncid;
  
  ///open the file in define mode
  status=nc_open(path.c_str(),NC_WRITE,&ncid);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"nc_open(\""+path+"\",NC_WRITE,) error");
  
  status=nc_redef(ncid);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"nc_redef() error");
  
  ///call poc_def_att(int,const poc_att_t&)
  status=poc_def_att(ncid,att);
  if(verbose)NC_CHKERR_BASE_LINE(status,"poc_def() error with "+path);
  
  ///-takes the file out of define mode and closes it
  end_status=nc_enddef(ncid);
  if(verbose)NC_CHKERR_BASE_LINE(end_status,"nc_enddef() error with "+path);
  
  close_status=nc_close(ncid);
  if(verbose)NC_CHKERR_BASE_LINE(close_status,"nc_close() error with "+path);
  
  return status==NC_NOERR?(end_status==NC_NOERR?close_status:end_status):status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_update_history(const string &path, const string &cmd, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Put global attribute
/*----------------------------------------------------------------------------*/
{
  int status;
  
  const time_t
    now=time(0);
  struct tm gmt;
  const int dateStrLen=99+cmd.length();
  char dateStr[dateStrLen];
  string historyStr;
  
  gmtime_r(&now,&gmt);
  strftime(dateStr,dateStrLen,"%F %T : ",&gmt);
  historyStr=dateStr+cmd;
  
  poc_att_t history;
  status=poc_inq_att(path,"history",&history,-1);
  
  if(status!=0){
    status=poc_def_att(path,poc_att_t("history",historyStr),verbose);
    }
  else{
    status=poc_def_att(path,poc_att_t("history",historyStr+'\n'+history.as_string()),verbose);
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_def_var(int ncid,const poc_var_t & var,int *varid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Define a variable
/**
\param ncid NetCDF file ID
\param var variable
\param *varid optional pointer to the NetCDF ID of the new variable
\returns NC_NOERR if success or the NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  int status,def_status,varid0,
    ndims=var.dimensions.size(),
    *dimids=new int[ndims],
    i,format;
  if(varid==NULL){
    varid=&varid0;
    }
  
  if(var.name==""){
    STDERR_BASE_LINE("FAIL-SAFE : skipping empty-named variable...\n");
    return 0;
    }
  
/*-----------------------------------------------------------------------------
  puts the file in define mode */
  def_status=nc_redef(ncid);
  if(def_status!=NC_NOERR && def_status!=NC_EINDEFINE)NC_TRAP_ERROR(return,def_status,1,"nc_redef() error");
  
/*-----------------------------------------------------------------------------
  defines the dimensions with poc_def_dim() */
  for(i=0;i<ndims;i++){
    const poc_dim_t *dim=&var.dimensions[i];
    status=poc_def_dim(ncid,*dim,&dimids[i]);
    if(status!=NC_NOERR)
      NC_TRAP_ERROR(return,status,1,"poc_def_dim(,(\""+dim->name+"\",%d),) error",dim->len);
    }
  
/*-----------------------------------------------------------------------------
  defines the variable */
  status=nc_def_var(ncid,var.name.c_str(),var.type,var.dimensions.size(),dimids,varid);
  if(status==NC_ENAMEINUSE){
    poc_var_t var0;
    poc_inq_var(ncid,var.name,&var0);
    if(var0.dimensions==var.dimensions){
      status=0;
      *varid=var0.id;
      goto defatt;
      }
    NC_TRAP_ERROR(return,status,1,"nc_def_var(,\""+var.name+"\",%d,%d dimensions,) error with different %d dimensions",var.type,var.dimensions.size(),var0.dimensions.size());
    }
  if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,1,"nc_def_var(,\""+var.name+"\",%d,%d,) error",var.type,var.dimensions.size());
  
/*-----------------------------------------------------------------------------
  if possible, turn on compression, level 1, because anything higher seams to be much slower with little extra space gain */
  status=nc_inq_format(ncid,&format);
  if(status!=NC_NOERR){
    NC_CHKERR_BASE_LINE(status,"nc_inq_format() error");
    format=0;
    }
  
  switch(format){
    case NC_FORMAT_NETCDF4:
    case NC_FORMAT_NETCDF4_CLASSIC:
      #ifdef NC_NETCDF4
      status=nc_def_var_deflate(ncid,*varid,1,1,1);
      if(status!=NC_NOERR)NC_CHKERR_BASE_LINE(status,"nc_def_var_deflate() error with "+var.name);
      #else
      #warning UPGRADE YOUR NetCDF : nc_def_var_deflate() not available
      NC_CHKERR_BASE_LINE(ENOEXEC,"nc_def_var_deflate() not available on this machine !!! UPGRADE YOUR NetCDF !!!");
      #endif
    }
  
/*-----------------------------------------------------------------------------
  defines the attributes with poc_def_att() */
defatt:
  for(i=0;i<var.attributes.size();i++){
    const poc_att_t *att=&var.attributes[i];
    if(att->name==POC_GRIB_FILE_POS_VATT_NAME)
      continue;
    if(att->name==_FillValue and att->type!=var.type){
      STDERR_BASE_LINE_FUNC("*** \""+att->name+"\" HAS TYPE %d WHEN VAR \""+var.name+"\" HAS TYPE %d: SKIPPING ***\n",att->type,var.type);
      continue;
      }
    status=poc_def_att(ncid,*varid,*att);
    NC_TRAP_ERROR(return,status,1,"poc_def_att(,,(\""+att->name+"\")) error");
    }
  
/*-----------------------------------------------------------------------------
  takes the file out of define mode if it was not so already */
  if(def_status==NC_NOERR){
    status=nc_enddef(ncid);
    NC_CHKERR_BASE_LINE(status,"nc_enddef() error");
    if(status==NC_EBADTYPE)
      poc_print(var);
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_def_var(const string & path, const poc_var_t & var, int *varid, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Define a variable
/**
\param path
\param var variable
\param verbose Default:0
\returns NC_NOERR if success or the NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  int status,end_status,close_status,ncid;
  
  ///-opens the file in define mode
  status=nc_open(path.c_str(),NC_WRITE,&ncid);
  if(status!=NC_NOERR) {
    NC_TRAP_ERROR(return,status,verbose,"nc_open(\""+path+"\",NC_WRITE,) error");
    }
    
  status=nc_redef(ncid);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"nc_redef() error");
  
  ///-calls poc_def_var(int,const poc_var_t&)
  status=poc_def_var(ncid,var,varid);
  if(verbose)NC_CHKERR_BASE_LINE(status,"poc_def_var(,(\""+var.name+"\"),) error with "+path);
  
  ///-takes the file out of define mode and closes it
  end_status=nc_enddef(ncid);
  if(verbose)NC_CHKERR_BASE_LINE(end_status,"nc_enddef() error with "+path);
  close_status=nc_close(ncid);
  if(verbose)NC_CHKERR_BASE_LINE(close_status,"nc_close() error with "+path);
  
  if(status==NC_NOERR){
    if(end_status==NC_NOERR)
      return close_status;
    return end_status;
    }
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_def(int ncid,const poc_global_t & global)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Define everything
/**
\param ncid NetCDF file ID
\param global
\returns NC_NOERR if success or the NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  int status,def_status,i;
  
  ///-puts the file in define mode
  def_status=nc_redef(ncid);
  if(def_status!=NC_NOERR && def_status!=NC_EINDEFINE)NC_TRAP_ERROR(return,def_status,1,"nc_redef() error");
  
  ///-defines the dimensions with poc_def_dim()
  for(i=0;i<global.dimensions.size();i++){
    status=poc_def_dim(ncid,global.dimensions[i]);
    if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"poc_def_dim() error");
    }
  
  ///-defines the variables with poc_def_var()
  for(i=0;i<global.variables.size();i++){
    status=poc_def_var(ncid,global.variables[i]);
    if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"poc_def_var() error");
    }
  
  ///-defines the attributes with poc_def_att()
  for(i=0;i<global.attributes.size();i++){
    status=poc_def_att(ncid,global.attributes[i]);
    if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"poc_def_att() error");
    }
  
  ///-takes the file out of define mode if it was not so already
  if(def_status==NC_NOERR){
    status=nc_enddef(ncid);
    NC_CHKERR_BASE_LINE(status,"nc_enddef() error");
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_create(const string &path,const poc_global_t & global,int verbose,int cmode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Create a new netcdf file
/**
\param path
\param global
\param verbose Default:0
\param cmode Default:0. You can set this to NC_NETCDF4|NC_CLASSIC_MODEL
\returns NC_NOERR if success or the NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  int status,close_status,ncid;
  poc_global_t fast=global;
  int i;
  
  #if 0
  if(compression){
    #ifdef NC_NETCDF4
    compression=NC_NETCDF4|NC_CLASSIC_MODEL;
    #else
    #warning UPGRADE YOUR NetCDF : NC_NETCDF4 NC_CLASSIC_MODEL are undefined
    NC_CHKERR_BASE_LINE(ENOEXEC,"compression not available on this machine !!! UPGRADE YOUR NetCDF !!!");
    #endif
    }
  #endif
  
  ///-creates the files with nc_create()
  /* When a netCDF dataset is created, is is opened NC_WRITE */
  status=nc_create(path.c_str(),cmode,&ncid);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"nc_create(\""+path+"\",%d,) error",cmode);
  
  fast.variables=poc_list_t<poc_var_t>();
  
  ///-calls poc_def(int,const poc_dim_t&)
  status=poc_def(ncid,fast);
  if(verbose)NC_CHKERR_BASE_LINE(status,"poc_def() error with "+path);
  
  ///-closes the file
  close_status=nc_close(ncid);
  if(verbose)NC_CHKERR_BASE_LINE(close_status,"nc_close() error with "+path);
  
  ///-defines the variables with poc_def_var(const string&,const poc_var_t&)
  for(i=0;i<global.variables.size();i++){
    const poc_var_t *var=&global.variables[i];
    status=poc_def_var(path,*var);
    if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose>=0,"poc_def_var(\""+path+"\",) error with \""+var->name+"\"");
    }
  
  if(status==NC_NOERR)
    status=close_status;

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_create(const string &path,const string & production,int verbose,int cmode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper
/**
\param path
\param production
\param verbose Default:0
\returns NC_NOERR if success or the NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  int status;
  poc_global_t global;
  
  global<<poc_att_t("production",production);
  
  status=poc_create(path,global,verbose,cmode);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void poc_print(const poc_dim_t & dim,ostream &out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///prints the properties of a poc_dim_t
/*----------------------------------------------------------------------------*/
{
  out<<dim.id<<":"<<dim.name<<"="<<dim.len;
  if(dim.isunlimited)out<<" currently";
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<class T> inline void poc_print_type_template(const T &arg,ostream &out=cout)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///prints the type of a poc_att_t or a poc_var_t
/*----------------------------------------------------------------------------*/
{
  switch(arg.type){
    #define caseType(x) case x:out<<#x;break
    caseType(NC_NAT);
    caseType(NC_BYTE);
    caseType(NC_CHAR);
    caseType(NC_SHORT);
    caseType(NC_INT);
    caseType(NC_FLOAT);
    caseType(NC_DOUBLE);
    #ifdef NC_UBYTE
    caseType(NC_UBYTE);
    caseType(NC_USHORT);
    caseType(NC_UINT);
    caseType(NC_INT64);
    caseType(NC_UINT64);
    caseType(NC_STRING);
    default:
      NC_CHKERR_BASE_LINE(NC_EBADTYPE,"unknown type %d",arg.type);exit(-2);
    #else
    default:
      NC_CHKERR_BASE_LINE(NC_EBADTYPE,"unknown type %d. You may need to UPGRADE YOUR NetCDF!",arg.type);exit(-2);
    #endif
    #undef caseType
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void poc_print_type(const poc_att_t &arg,ostream &out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_print_type_template(arg,out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void poc_print_type(const poc_var_t &arg,ostream &out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_print_type_template(arg,out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void poc_print(const poc_att_t & att,ostream &out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///prints the properties of a poc_att_t
/*----------------------------------------------------------------------------*/
{
  int i;
  
  out<<att.id<<":"<<att.name<<"=("<<att.len<<"x";
  poc_print_type(att,out);
  out<<")";
  switch(att.type){
    case NC_CHAR:
    case NC_STRING:
      out<<'"'<<att.as_charp()<<'"';
      break;
    case NC_BYTE:
    case NC_SHORT:
    case NC_INT:
    case NC_FLOAT:
    case NC_DOUBLE:
      for(i=0;i<att.len;i++){
        sprintf(out,"%.11Lg ",att[i]);
        }
      break;
    #ifdef NC_UBYTE
    case NC_UBYTE:
    case NC_USHORT:
    case NC_UINT:
    case NC_INT64:
    case NC_UINT64:
      for(i=0;i<att.len;i++){
        sprintf(out,"%.11Lg ",att[i]);
        }
      break;
    default:
      NC_CHKERR_BASE_LINE(NC_EBADTYPE,"unknown type %d",att.type);exit(-2);
    #else
    default:
      NC_CHKERR_BASE_LINE(NC_EBADTYPE,"unknown type %d. You may need to UPGRADE YOUR NetCDF!",att.type);exit(-2);
    #endif
    }
  out<<";";
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void poc_print(const poc_var_t & var,ostream &out,int verbose)

/*----------------------------------------------------------------------------*/
///prints the properties of a poc_var_t
/*----------------------------------------------------------------------------*/
{
  poc_print_type(var,out);
  out<<" "<<var.id<<":"<<var.name<<"(";
  poc_print(",",var.dimensions,out);
  out<<")";
  if(var.attributes.size()>0 and verbose>0){
    out<<"\n  ";
    poc_print("\n  ",var.attributes,out);
    }
  out<<"\n";
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void poc_print(const poc_global_t &global,ostream &out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///prints the properties of a poc_global_t
/*----------------------------------------------------------------------------*/
{
  poc_print(" ",global.dimensions,out);
  out<<"\n";
  poc_print("",global.variables,out);
  poc_print("",global.attributes,out);
}

