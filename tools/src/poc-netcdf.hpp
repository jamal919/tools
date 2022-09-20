
/*******************************************************************************

  T-UGO tools, 2006-2013

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief new (as of 2012-01-24) poc-netcdf classes

Trying to set up another set of poc-netcdf classes
\dot
digraph TUGOm {
  rankdir=LR;
  name -> variable
  name -> dimension
  name -> attribute
  dimensions -> variable
  attributes -> variable
  dimensions -> global
  attributes -> global
  }
\enddot
*/
/*----------------------------------------------------------------------------*/
/** \page tips Coding tips

\section var add attributes to a variable
\code
  poc_var_t var;
  ...
  var << poc_att_t("offset",-2713.15) << poc_att_t("comment","example");
\endcode
*/
/*----------------------------------------------------------------------------*/

#if POC_NETCDF_HPP == 0
#define POC_NETCDF_HPP 1

#include <stdio.h> //asprintf
#include <stdlib.h> //exit
#include <stdint.h> //uint64_t
#include <string.h>
#include "netcdf-proto.h"

#if TUGO
#include "poc-netcdf.h"
#include <deque>
#endif

#include <iostream> //cout

#include "functions.h"
#include "poc-list.hpp"

#include "tools-define.h"
#include "version-macros.def"
#ifndef __LINE_FILE_PACKAGE_REVISION
#define __LINE_FILE_PACKAGE_REVISION  "line " TOSTRING(__LINE__) " of " __FILE__
#endif

#define NC_FILL_AMPLITUDE -1.


class poc_dim_t : public poc_name_id_t{
public:
  size_t len;      ///<dimension length
  bool   isunlimited;

private:
  
  void init(){
    id=NC_EBADDIM;/* Invalid dimension id or name */
    len=0;
    isunlimited=false;
    }
  
  void init(size_t len0,bool isunlimited0){
    len=len0;
    isunlimited=isunlimited0;
    }
  
  void init(const poc_dim_t &src){
    init_name_and_id(src);
    len=src.len;
    isunlimited=src.isunlimited;
    }
  
public:
  
  void init(size_t len0){
    init(len0,len0==NC_UNLIMITED);
    }
  
  poc_dim_t(){
    init();
    }
  
  poc_dim_t(const poc_dim_t &src){
    init(src);
    }
  
  poc_dim_t &operator = (const poc_dim_t &src){
    init(src);
    return *this;
    }
  
  /*
  void init(const string name0,int id0,int nisunlimited){
    init_name_and_id(name0,id0);
    isunlimited=nisunlimited;
    }
  */
  
  poc_dim_t(const string name0,size_t len0){
    name=name0;
    init(len0);
    }
  
  poc_dim_t(const string name0,size_t len0,bool isunlimited0){
    init(name0,len0,isunlimited0);
    }
  
  poc_dim_t(size_t len0){
    asprintf(name,"n%u",len0);
    init(len0);
    }
  
  void init(const string name0,size_t len0,bool isunlimited0){
    name=name0;
    init(len0,isunlimited0);
    }
  
  void init(const string name0,size_t len0){
    name=name0;
    init(len0);
    }
  
  bool operator==(const poc_dim_t &src) const{
    return isunlimited==src.isunlimited && (isunlimited || len==src.len);
    }
  
  bool operator!=(const poc_dim_t &src) const{
    return !operator==(src);
    }
  };

class poc_att_t : public poc_name_id_t{
public:
  nc_type type;
  size_t len; ///<number of elements
  
private:
  
  union{
    #define unionDataType(t) t *data_as_ ## t
    unionDataType(char);
    unionDataType(int8_t);
    unionDataType(uint8_t);
    unionDataType(int16_t);
    unionDataType(uint16_t);
    unionDataType(int32_t);
    unionDataType(uint32_t);
    unionDataType(int64_t);
    unionDataType(uint64_t);
    unionDataType(float);
    unionDataType(double);
    #undef unionDataType
    char **data_as_string;
    void *data;
    };
  
  void init(){
    type=NC_NAT; /* NAT = 'Not A Type' (c.f. NaN) */
    id=NC_ENOTATT;
    len=0;
    data=NULL;
    }
  
  void init(const poc_att_t &src){
    init_name_and_id(src);
    
    type=src.type;
    
    len=src.len;
    
    switch(type){
    case NC_NAT: /* SADLY THIS CAN HAPPEN ... */
      data=0; /* fail safe */
      break;
    case NC_CHAR:                            /* ISO/ASCII character */
      data_as_char=poc_strdup(src.data_as_char);break;
      #define copy_data_as__t(t) data_as_ ## t=copy(src.data_as_ ## t,len);break
    case NC_BYTE:   copy_data_as__t(int8_t);   /* signed 1 byte integer */
    case NC_SHORT:  copy_data_as__t(int16_t);  /* signed 2 byte integer */
    case NC_INT:    copy_data_as__t(int32_t);  /* signed 4 byte integer */
    case NC_FLOAT:  copy_data_as__t(float);    /* single precision floating point number */
    case NC_DOUBLE: copy_data_as__t(double);   /* double precision floating point
number */
    #ifdef NC_UBYTE
    case NC_UBYTE:  copy_data_as__t(uint8_t);  /* unsigned 1 byte int */
    case NC_USHORT: copy_data_as__t(uint16_t); /* unsigned 2-byte int */
    case NC_UINT:   copy_data_as__t(uint32_t); /* unsigned 4-byte int */
    case NC_INT64:  copy_data_as__t(int64_t);  /* signed 8-byte int */
    case NC_UINT64: copy_data_as__t(uint64_t); /* unsigned 8-byte int */
    case NC_STRING:
      data_as_string=0;
      copy_data(src.data_as_string,len);
      break;
    default: NC_CHKERR_BASE_LINE(NC_EBADTYPE,"unknown type %d",type);wexit(NC_EBADTYPE);
    #else
    #warning UPGRADE YOUR NetCDF : NC_UBYTE NC_USHORT NC_UINT NC_INT64 NC_UINT64 NC_STRING are undefined
    default: NC_CHKERR_BASE_LINE(NC_EBADTYPE,"unknown type %d. You may need to UPGRADE YOUR NetCDF!",type);wexit(NC_EBADTYPE);
    #endif
      #undef copy_data_as__t
      }
    }
  

public:
  
  poc_att_t(){
    init();
    }
  
  poc_att_t(const poc_att_t &src){
    init(src);
    }
  
  poc_att_t &operator = (const poc_att_t &src){
    destroy();
    init(src);
    return *this;
    }
  
  void init_data(char* data0, const int & nn=-1){
    destroy();
    if(nn<0){
      len=strlen(data0);
      }
    else{
      len=nn;
      }
    data_as_char=data0;
    type=NC_CHAR;
    }
  
  void copy_data(const char*data0,const int &nn=-1){
    destroy();
    if(nn<0){
      len=strlen(data0);
      }
    else{
      len=nn;
      }
    data_as_char=copy(data0,len+1);
    type=NC_CHAR;
    }
  
  void copy_data(const string &data0){
    copy_data(data0.c_str(),data0.length());
    }
  
  #define NC_UnionDataType(NC_T,t) \
    void init_data(      t *data0,int nlen){destroy();len=nlen;type=NC_T;data_as_ ## t=data0;} \
    void copy_data(const t *data0,int nlen){destroy();len=nlen;type=NC_T;data_as_ ## t=copy(data0,len);} \
    void copy_data(const t &data0         ){destroy();len=1   ;type=NC_T;data_as_ ## t=copy(&data0,1);}
  NC_UnionDataType(NC_BYTE,  int8_t)   /* signed 1 byte integer */
  NC_UnionDataType(NC_SHORT, int16_t)  /* signed 2 byte integer */
  NC_UnionDataType(NC_INT,   int32_t)  /* signed 4 byte integer */
  NC_UnionDataType(NC_FLOAT, float)    /* single precision floating point number */
  NC_UnionDataType(NC_DOUBLE,double)   /* double precision floating point number */
  #ifdef NC_UBYTE
  NC_UnionDataType(NC_UBYTE, uint8_t)  /* unsigned 1 byte int */
  NC_UnionDataType(NC_USHORT,uint16_t) /* unsigned 2-byte int */
  NC_UnionDataType(NC_UINT,  uint32_t) /* unsigned 4-byte int */
  NC_UnionDataType(NC_INT64, int64_t)  /* signed 8-byte int */
  NC_UnionDataType(NC_UINT64,uint64_t) /* unsigned 8-byte int */
  void init_data(char **data0,int nlen){destroy();len=nlen;type=NC_STRING;data_as_string=data0;}
  void copy_data(char **data0,int nlen){
    destroy();len=nlen;type=NC_STRING;
    data_as_string=new char *[len];
    for(int i=0;i<len;i++)
      data_as_string[i]=strdup(data0[i]);
    }
  #endif
  #undef NC_UnionDataType
  
template<typename T> void init_copy(const string &name0,const T &data0){
    name=name0;
    copy_data(data0);
    }
  
template<typename T> poc_att_t(const string &name0,const T &data0){
    data=NULL;
    init_copy(name0,data0);
    }
  
  poc_att_t(const string &name0,const char *data0){
    data=NULL;
    init_copy(name0,data0);
    }
  
template<typename T> void init(const string &name0,T *data0,const int &nn){
    name=name0;
    init_data(data0,nn);
    }
  
template<typename T> void init_copy(const string &name0,const T *data0,const int &nn){
    name=name0;
    copy_data(data0,nn);
    }
  
  long double operator[](const unsigned int i) const{
    if(i>=len) TRAP_ERR_EXIT(ENOEXEC,"programming error : i=%u>=%u\n",i,len);
    switch(type){
      #define return_data_as__t(t) return data_as_ ## t[i];
    case NC_CHAR:   return_data_as__t(char)     /* ISO/ASCII character */
    case NC_BYTE:   return_data_as__t(int8_t)   /* signed 1 byte integer */
    case NC_SHORT:  return_data_as__t(int16_t)  /* signed 2 byte integer */
    case NC_INT:    return_data_as__t(int32_t)  /* signed 4 byte integer */
    case NC_FLOAT:  return_data_as__t(float)    /* single precision floating point number */
    case NC_DOUBLE: return_data_as__t(double)   /* double precision floating point number */
    #ifdef NC_UBYTE
    case NC_UBYTE:  return_data_as__t(uint8_t)  /* unsigned 1 byte int */
    case NC_USHORT: return_data_as__t(uint16_t) /* unsigned 2-byte int */
    case NC_UINT:   return_data_as__t(uint32_t) /* unsigned 4-byte int */
    case NC_INT64:  return_data_as__t(int64_t)  /* signed 8-byte int */
    case NC_UINT64: return_data_as__t(uint64_t) /* unsigned 8-byte int */
    case NC_STRING: return data_as_string[0][i];
    default: NC_CHKERR_BASE_LINE(NC_EBADTYPE,"unknown type %d",type);wexit(NC_EBADTYPE);
    #else
    default: NC_CHKERR_BASE_LINE(NC_EBADTYPE,"unknown type %d. You may need to UPGRADE YOUR NetCDF!",type);wexit(NC_EBADTYPE);
    #endif
      #undef return_data_as__t
      }
    }
  
  const string as_string() const{
    if(type==NC_CHAR){
      return string(data_as_char);
      }
    if(type==NC_STRING){
      if(len==1)
        return string(data_as_string[0]);
      }
    else if(len==1){
      ostringstream oss;
      
      oss << this->operator[](0u);
      
      return string(oss.str());
      }
    return string();
    }
  
  const char *as_charp() const{
    if(type==NC_CHAR){
      return data_as_char;
      }
    if(type==NC_STRING){
      if(len==1)
        return data_as_string[0];
      }
    return NULL;
    }
  
  void destroy(){
    if(data!=NULL){
      if(type==NC_STRING){
        deletep2D(&data_as_string,len,free);
        }
      else{
        //delete[]data;//warning: deleting ‘void*’ is undefined [enabled by default]
        /* IF IT CRASHES BELOW, MAKE SURE YOU USE valgrind */
        delete[]data_as_char;
        }
      data=NULL;
      }
    }
  
  ~poc_att_t(){
    destroy();
    }
  
  friend int poc_def_att(int ncid,int varid,const poc_att_t &att);
  };


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<class T> void poc_print(const string & s,const poc_list_t<T> &list,ostream &out=cout)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  for(int i=0;i<list.size();i++){
    if(i>0)
      out << s;
    poc_print(list[i],out);
    }
}


#define CF_CONVENTIONS_URL "http://cfconventions.org/cf-conventions/cf-conventions.html"
#define standard_names_URL CF_CONVENTIONS_URL "#standard-name"
extern string comodo_standard_name(const string & varName="",const string waveName="");
extern string long_name_from_varname(const string & varName="");


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

class poc_var_t : public poc_name_id_t

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
public:
  nc_type type;
  
  poc_list_t<poc_att_t> attributes;
  poc_list_t<poc_dim_t> dimensions;
  
  string axes;
  
private:
  
  void init(){
    type=NC_NAT; /* NAT = 'Not A Type' (c.f. NaN) */
    id=NC_ENOTVAR;
    }
  
  void init(const poc_var_t &src){
    init_name_and_id(src);
    
    type=src.type;
    
    attributes=src.attributes;
    dimensions=src.dimensions;
    
    axes=src.axes;
    }

public:
  
  poc_var_t(){
    init();
    }
  
  poc_var_t(const poc_var_t &src){
    init(src);
    }
  
  poc_var_t &operator = (const poc_var_t &src){
    init(src);
    return *this;
    }
  
  poc_var_t(const string &name0,nc_type type0,int id0=NC_ENOTVAR){
    type=type0;
    init_name_and_id(name0,id0);
    }
  
  poc_var_t & init(const string &name0,nc_type type0,const string &standard_name="",const string &units="",double mask=NAN,const string &long_name=""){
    double scale,offset,mask0=NAN;
    int status;
    
    extern int poc_decode_mask(const poc_var_t &var, double *scale, double *offset, double    *spec, double    mask_value=NC_FILL_DOUBLE);
    
    /* if there are no _FillValue attribute, the following line will not change mask */
    status=poc_decode_mask(*this,&scale,&offset,&mask0,mask);
    if(isnan(mask))
      mask=mask0;
    
    name=name0;
    type=type0;
    /// see #standard_names_URL
    if(standard_name!=""){
      if(standard_name[0]=='_'){
        string computed_standard_name=comodo_standard_name(name,standard_name.substr(1));
        attributes << poc_att_t("standard_name",computed_standard_name);
        }
      else
        attributes << poc_att_t("standard_name",standard_name);
      }
    else
      attributes.erase("standard_name");
    
    if(units!="")
      attributes << poc_att_t("units",units);
    else
      attributes.erase("units");
    
    /* clean-up useless attributes */
    
    if(scale==1.)
      attributes.erase("scale_factor");
    
    if(offset==0.)
      attributes.erase("add_offset");
    
    /* keep or set _FillValue attribute only if necessary */
    
    attributes.erase(_FillValue);
    attributes.erase("missing_value");
    
    if(not isnan(mask))
      switch(type){
        case NC_BYTE:  if(mask!=NC_FILL_BYTE  )attributes << poc_att_t(_FillValue,(int8_t)mask);break;
        case NC_CHAR:  if(mask!=NC_FILL_CHAR  )attributes << poc_att_t(_FillValue,(char)mask);break;
        case NC_SHORT: if(mask!=NC_FILL_SHORT )attributes << poc_att_t(_FillValue,(int16_t)mask);break;
        case NC_INT:   if(mask!=NC_FILL_INT   )attributes << poc_att_t(_FillValue,(int32_t)mask);break;
        case NC_FLOAT: if(mask!=NC_FILL_FLOAT )attributes << poc_att_t(_FillValue,(float)mask);break;
        case NC_DOUBLE:if(mask!=NC_FILL_DOUBLE)attributes << poc_att_t(_FillValue,mask);break;
        #ifdef NC_UBYTE
        case NC_UBYTE: if(mask!=NC_FILL_UBYTE )attributes << poc_att_t(_FillValue,(uint8_t)mask);break;
        case NC_USHORT:if(mask!=NC_FILL_USHORT)attributes << poc_att_t(_FillValue,(uint16_t)mask);break;
        case NC_UINT:  if(mask!=NC_FILL_UINT  )attributes << poc_att_t(_FillValue,(uint32_t)mask);break;
        case NC_INT64: if(mask!=NC_FILL_INT64 )attributes << poc_att_t(_FillValue,(int64_t)mask);break;
        case NC_UINT64:if(mask!=NC_FILL_UINT64)attributes << poc_att_t(_FillValue,(uint64_t)mask);break;
        case NC_STRING: TRAP_ERR_EXIT(ENOEXEC,"not coded for NC_STRING\n");
        default: NC_CHKERR_BASE_LINE(NC_EBADTYPE,"unknown type %d",type);wexit(NC_EBADTYPE);
        #else
        default: NC_CHKERR_BASE_LINE(NC_EBADTYPE,"unknown type %d. You may need to UPGRADE YOUR NetCDF!",type);wexit(NC_EBADTYPE);
        #endif
        }
    
    /* for Ferret */
    if(long_name!=""){
      if(long_name=="_"){
        string computed_long_name=long_name_from_varname(name);
        attributes << poc_att_t("long_name",computed_long_name);
        }
      else
        attributes << poc_att_t("long_name",long_name);
      }
    else
      attributes.erase("long_name");
    
    return *this;
    }
  
  poc_var_t & operator << (const poc_att_t &src) {
    attributes.set(src);
    return *this;
    }
  
  poc_var_t & operator << (const poc_dim_t &src) {
    dimensions.set(src);
    return *this;
    }
  
  void init(const poc_dim_t &dim,nc_type type0) {
    ///dimension variable
    init(dim.name,type0);
    dimensions.clear();
    attributes.clear();
    *this << dim;
    }
  
  void clear(){
    name.clear();
    dimensions.clear();
    attributes.clear();
    }
  };


class poc_global_t{
public:
  
  poc_list_t<poc_att_t> attributes;
  poc_list_t<poc_dim_t> dimensions;
  
  poc_list_t<poc_var_t> variables;
  
private:
  
  void init(const poc_global_t &src){
    attributes=src.attributes;
    dimensions=src.dimensions;
    
    variables=src.variables;
    }

public:
  
  poc_global_t(const poc_global_t &src){
    init(src);
    }
  
  void init(const string prod="initialised around " __LINE_FILE_PACKAGE_REVISION){
    *this << poc_att_t("production",prod);
    
    const int dateL=32;
    char date[dateL];
    ///-get the date and time with get_date()
    get_date(date,dateL);
    *this << poc_att_t("creation_date",date);
    
    *this << poc_att_t("Conventions","CF-1.5 COMODO-0.1");
    }
  
  poc_global_t(const string prod="constructed around " __LINE_FILE_PACKAGE_REVISION){
    init(prod);
    }
  
  poc_global_t & operator = (const poc_global_t &src){
    init(src);
    return *this;
    }
  
  poc_global_t & operator << (const poc_att_t &src){
    attributes.set(src);
    return *this;
    }
  
  poc_global_t & operator << (const poc_dim_t &src){
    dimensions.set(src);
    return *this;
    }
  
  poc_global_t & operator << (const poc_var_t &src){
    variables.set(src);
    return *this;
    }
  };
  
class poc_grid_t{
  private:
  public:
  int      id;
  char     *name;
  poc_var_t lon,lat,time;

  poc_grid_t(){
    id = -1;
    name = 0;
    }

/*    ~poc_grid_t(){
      id = -1;
      if(name != 0) delete[] name; name = 0;
      if(filename != 0) delete[] filename; filename = 0;
    }*/
};


//from poc-netcdf-iio.cpp
extern int poc_inq_dim(int ncid,int dimid,poc_dim_t *dim);
extern int poc_inq_dim(const string &path,int dimid, poc_dim_t *dim);

extern int poc_inq_att(int ncid,int varid,const char*attname,poc_att_t *att,int verbose=0);
extern int poc_inq_att(int ncid,const char*attname,poc_att_t *att,int verbose=0);
extern int poc_inq_att(const string &path,const char*attname,poc_att_t *att,int verbose=0);
extern int poc_inq_var(int ncid,int varid,poc_var_t *var);
extern int poc_inq_var(int ncid,const string &varname,poc_var_t *var,int verbose=0);
extern int poc_inq_var(const string &path,const string &varname,poc_var_t *var,int verbose=0);

extern int poc_inq(int ncid,poc_global_t *global, int *format=0);
extern int poc_inq(const string &path,poc_global_t *global,int verbose=0, int *format=0);

extern int poc_def_dim(int ncid,const poc_dim_t &dim,int *dimid=NULL);
extern int poc_def_att(int ncid, int varid, const poc_att_t& att);
extern int poc_def_att(int ncid,const poc_att_t &att);
extern int poc_def_att(const string &path,const poc_att_t &att,int verbose=0);

extern int poc_update_history(const string &path, const string &cmd, int verbose=0);

extern int poc_def_var(int ncid,const poc_var_t &var,int *varid=NULL);
extern int poc_def_var(const string &path,const poc_var_t &var, int *varid=NULL, int verbose=0);
extern int poc_def(int ncid,const poc_global_t &global);
extern int poc_create(const string &path,const poc_global_t &global,int verbose=0,int cmode=0);
extern int poc_create(const string &path,const string & production,int verbose=0,int cmode=0);

extern void poc_print_type(const poc_att_t &arg,ostream &out=cout);
extern void poc_print_type(const poc_var_t &arg,ostream &out=cout);

extern void poc_print(const poc_dim_t &dim,ostream &out=cout);
extern void poc_print(const poc_att_t& att,ostream &out=cout);
extern void poc_print(const poc_var_t& var,ostream &out=cout,int verbose=1);
extern void poc_print(const poc_global_t& global,ostream &out=cout);


//from poc-netcdf-comodo-iio.cpp
extern int poc_get_frame_names(const string &filename,poc_list_t<poc_name_id_t> *frame_names,int verbose=0,string *framevarname=NULL);
extern int poc_get_frame_from_name(const string &filename,const string &name,int *frame,int verbose=0);
extern bool isT(const poc_dim_t &dim);
extern int poc_find_timevarid(int ncid,int *tvid);
extern int poc_find_timevarid(const string & filename,int *tvid);
extern int poc_gettime(int ncid, date_t *origine, double **time, size_t *nframes, poc_var_t *timevar);
extern int poc_gettime(const string & filename, date_t *origine, double **time, size_t *nframes, poc_var_t *timevar);
#define Dimensionless_Vertical_Coordinates_URL CF_CONVENTIONS_URL "#dimensionless-v-coord"
extern bool is_Dimensionless_Vertical_Coordinates(const string &sN);
extern int findStandardName(const poc_global_t &global,const string *sNs,const int sNC,const poc_list_t<poc_dim_t> &dimensions,int verbose);
extern bool nameS_start_or_end_with(const poc_var_t &var,const string &prefix);
extern int find1DVariableThasIs(const poc_global_t&global,const char*prefix,int verbose);
extern int find2DVariableThasIs(const poc_global_t&global,int dimlen,int verbose);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_gettime(T file, date_t *origine, double **time, int *nframes, poc_var_t *timevar=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** Get times of frames of a given NetCDF file
See poc_gettime(const char*,date_t*,double**,size_t*,poc_var_t*) for more help.
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

extern int poc_get_var_length(const poc_var_t &info, size_t *length, int *nspacedim=NULL, int *nframes=NULL);
extern int poc_get_var_length(const string & path,const string & varname,size_t *length, int *nspacedim=NULL, int *nframes=NULL, int verbose=0);

extern const poc_dim_t * findTimeDim(const poc_var_t & info,int *index=NULL);
extern poc_dim_t * findTimeDim(poc_var_t & info,int *index=NULL);
extern poc_dim_t * findTimeDim(poc_global_t & glob,int *index=0);
extern bool unlimitTimeDim(poc_var_t *info);

extern int poc_decode_mask(const poc_var_t &var, double *scale, double *offset, char      *spec, char      mask_value=NC_FILL_CHAR);
extern int poc_decode_mask(const poc_var_t &var, double *scale, double *offset, int8_t    *spec, int8_t    mask_value=NC_FILL_BYTE);
extern int poc_decode_mask(const poc_var_t &var, double *scale, double *offset, short int *spec, short int mask_value=NC_FILL_SHORT);
extern int poc_decode_mask(const poc_var_t &var, double *scale, double *offset, int       *spec, int       mask_value=NC_FILL_INT);
extern int poc_decode_mask(const poc_var_t &var, double *scale, double *offset, float     *spec, float     mask_value=NC_FILL_FLOAT);
extern int poc_decode_mask(const poc_var_t &var, double *scale, double *offset, double    *spec, double    mask_value=NC_FILL_DOUBLE);
/* avoid going into <complex> while debugging */
static complex<double> NC_FILL_COMPLEX=NC_FILL_DOUBLE;
extern int poc_decode_mask(const poc_var_t &var, double *scale, double *offset, complex<double> *spec, const complex<double> & mask_value=NC_FILL_COMPLEX);

extern int poc_get_var_length(const poc_global_t & global, const string & varname, poc_var_t *info, size_t *length, int *nspacedim=0);
extern int poc_get_var_length_and_mask(const poc_global_t & global, const string & varname, poc_var_t *info, size_t *length, int *nspacedim, float *spec);

#endif
