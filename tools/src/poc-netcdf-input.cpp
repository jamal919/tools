
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief poc-netcdf variable data input/output definitions
*/
/*----------------------------------------------------------------------------*/

#include "constants.h"
#include "poc-netcdf-io.h"
#include "poc-netcdf-data.hpp"

#include "poc-grib.h" /* for poc_grib_get_vara() */


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename Td,typename Ts> bool poc_scale_data_template(Td *data, int length, Ts scale, Ts offset, Td *mask, Td specific, int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// scales the data
/**
\param[in,out] data
\param[in,out] mask set to \c specific if scaling occurs
*/
/*----------------------------------------------------------------------------*/
{
  int i,parallel=length>10000000;
  timeval before;
  
  if(scale==1. && offset==0.)
    return false;
  
  if(verbose>0){
    STDERR_BASE_LINE("%s:length=%d => parallel=%d\n",__func__,length,parallel);
    gettimeofday(&before);
    }
  
#pragma omp parallel for if(parallel)
  for(i=0;i<length;i++){
    Td *datai=&data[i];
    if(*datai==*mask){
      *datai=specific;
      }
    else{
      *datai=*datai*scale+offset;
      }
    }
  
  if(verbose>0)STDERR_BASE_LINE("%s:%gs\n",__func__,difftime(before));
  
  *mask=specific;
  
  return true;
}
#if 0


extern bool poc_scale_data(Td *data, int length, Ts scale, Ts offset, Td *mask, int verbose=0);
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool poc_scale_data(Td *data, int length, Ts scale, Ts offset, Td *mask, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  Td specific=NC_FILL_;
  bool result;
  result=poc_scale_data_template(data, length, scale, offset, mask, specific, verbose);
  return result;
}


extern bool poc_scale_data(Td *data, int length, Ts scale, Ts offset, Td *mask, Td specific, int verbose=0);
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool poc_scale_data(Td *data, int length, Ts scale, Ts offset, Td *mask, Td specific, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  bool result;
  result=poc_scale_data_template(data, length, scale, offset, mask, specific, verbose);
  return result;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool poc_scale_data(double *data, int length, double scale, double offset, double *mask, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double specific=NC_FILL_DOUBLE;
  bool result;
  result=poc_scale_data_template(data, length, scale, offset, mask, specific, verbose);
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool poc_scale_data(float *data, int length, double scale, double offset, float *mask, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float specific=NC_FILL_FLOAT;
  bool result;
  result=poc_scale_data_template(data, length, scale, offset, mask, specific, verbose);
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool poc_scale_data(float *data, int length, double scale, double offset, float *mask, float specific, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  bool result;
  result=poc_scale_data_template(data, length, scale, offset, mask, specific, verbose);
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool poc_scale_data(short int *data, int length, double scale, double offset, short int *mask, short int specific, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  bool result;
  result=poc_scale_data_template(data, length, scale, offset, mask, specific, verbose);
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool poc_scale_data(signed char *data, int length, double scale, double offset, signed char *mask, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  signed char specific=NC_FILL_CHAR;
  bool result;
  result=poc_scale_data_template(data, length, scale, offset, mask, specific, verbose);
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int poc_scale_data_template(T *data, int length, const poc_var_t & info, T *mask=0, int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_scale_data(Td*,in,Ts,Ts,Td*,int)
/*----------------------------------------------------------------------------*/
{
  int status;
  double scale,offset;
  T _mask;
  if(mask==0)
    mask=&_mask;
  
  status=poc_decode_mask(info,&scale,&offset,mask);
  
  poc_scale_data(data,length,scale,offset,mask,verbose);
  
  return status;
}
#if 0
extern int poc_scale_data(T *data, int length, const poc_var_t & info, T *mask=0, int verbose=0);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_scale_data(T *data, int length, const poc_var_t & info, T *mask, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_scale_data_template(data,length,info,mask,verbose);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_scale_data(double *data, int length, const poc_var_t & info, double *mask, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_scale_data_template(data,length,info,mask,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_get_vara_template(int ncid,const poc_var_t &info,int frame, T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_get_vara(int,int,const size_t[],const size_t[],T*)
/*----------------------------------------------------------------------------*/
{
  int status,i;
  const int dimC=info.dimensions.size();
  const poc_dim_t *dim;
  size_t *start,*count;
  
  start=new size_t[dimC];
  count=new size_t[dimC];
  
  for(i=0;i<dimC;i++){
    dim=&info.dimensions[i];
    
    if(isT(*dim)){
      while(frame<0){
        frame+=dim->len;
        }
      start[i]=frame;
      count[i]=1;
      }
    else{
      start[i]=0;
      count[i]=dim->len;
      }
    }
  
  status=poc_get_vara(ncid,info.id,start,count,z);
  
  delete[]start;
  delete[]count;
  
  return status;
}
#if 0
extern int poc_get_vara(int ncid,const poc_var_t &info,int frame, T *z);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(int ncid,const poc_var_t &info,int frame, T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(ncid,info,frame,z);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(int ncid,const poc_var_t &info,int frame, double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(ncid,info,frame,z);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(int ncid,const poc_var_t &info,int frame, float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(ncid,info,frame,z);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(int ncid,const poc_var_t &info,int frame, short int *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(ncid,info,frame,z);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(int ncid,const poc_var_t &info,int frame, char *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(ncid,info,frame,z);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(int ncid,const poc_var_t &info,int frame, signed char *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(ncid,info,frame,z);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int poc_get_vara_template(const string &filename, const poc_var_t & info,int frame, T *z, int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_get_vara(int,V,size_t,T*)
/*----------------------------------------------------------------------------*/
{
  int ncid,status;
  
  status=nc_open(filename.c_str(),NC_NOWRITE,&ncid);
  
  if(status!=NC_NOERR){
    
/*----------------------------------------------------------------------------*/
    if(is_grib(filename)){
#ifdef HAVE_LIBGRIB_API
      status=poc_grib_get_vara(filename,info,frame,z,verbose);
      
      if(status!=0 and verbose>=0)
        STDERR_BASE_LINE("poc_grib_get_vara(\""+filename+"\",(\""+info.name+"\"),) error (%d %s)\n",status,grib_get_error_message(status));
      
      return status;
#else
      TRAP_ERR_RETURN(ENOEXEC,1,"Please compile with GRIB_API\n");
#endif
      }

/*----------------------------------------------------------------------------*/
    if(strrncmp(filename,".s2r")==0){
      status=quoddy_loadr1(filename.c_str(),NC_MAX_INT,z);
      return status;
      }
    
/*----------------------------------------------------------------------------*/
    
    NC_TRAP_ERROR(return,status,verbose,"nc_open(\""+filename+"\",NC_NOWRITE,) error");
    }
  
  status=poc_get_vara(ncid, info, frame, z);
  if(status!=NC_NOERR && verbose) NC_CHKERR_LINE_FILE(status,"poc_get_vara() error in "+filename);
  
  nc_close(ncid);
  
  return(status);
}
#if 0
extern int poc_get_vara(const string &filename, const poc_var_t &info,int frame, T *z, int verbose=0);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(const string &filename, const poc_var_t &info,int frame, T *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(filename,info,frame,z,verbose);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(const string &filename, const poc_var_t &info,int frame, double *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(filename,info,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(const string &filename, const poc_var_t &info,int frame, float *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(filename,info,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(const string &filename, const poc_var_t &info,int frame, short int *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(filename,info,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(const string &filename, const poc_var_t &info,int frame, char *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(filename,info,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(const string &filename, const poc_var_t &info,int frame, signed char *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(filename,info,frame,z,verbose);
  return status;
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_get_vara_template(int ncid,const poc_var_t &info,T *z,const string & axes="",const indexes_t & indexes= indexes_t() )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_get_vara(int,int,const size_t[],const size_t[],T*)
/*----------------------------------------------------------------------------*/
{
  const int iC=indexes.size();
  if(axes!="" and axes.size()!=iC) TRAP_ERR_EXIT(ENOEXEC,"programming error: %s called with axes and indexes of different sizes when axes is not empty: %d and %d\n",__func__,axes.size(),iC);
  
  int status,i;
  const int dimC=info.dimensions.size();
  const poc_dim_t *dim;
  size_t *start,*count;
  
  start=new size_t[dimC];
  count=new size_t[dimC];
  
  for(i=0;i<dimC;i++){
    dim=&info.dimensions[i];
    char iai;
    iai=info.axes[i];
    const size_t ai=axes.find(iai);
    
    if(ai==string::npos){
      start[i]=0;
      count[i]=dim->len;
      }
    else{
      int ii;
      ii=indexes[ai];
      while(ii<0)
        ii+=dim->len;
      start[i]=ii;
      count[i]=1;
      }
    }
  
  status=poc_get_vara(ncid,info.id,start,count,z);
  
  delete[]start;
  delete[]count;
  
  return status;
}
#if 0
extern int poc_get_vara(int ncid,const poc_var_t &info,T *z,const string & axes="",const indexes_t & indexes= indexes_t() );


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(int ncid,const poc_var_t &info,T *z,const string & axes,const indexes_t & indexes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(ncid,info,z,axes,indexe);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(int ncid,const poc_var_t &info,double *z,const string & axes,const indexes_t & indexes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(ncid,info,z,axes,indexes);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(int ncid,const poc_var_t &info,float *z,const string & axes,const indexes_t & indexes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(ncid,info,z,axes,indexes);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(int ncid,const poc_var_t &info,int8_t *z,const string & axes,const indexes_t & indexes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(ncid,info,z,axes,indexes);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(int ncid,const poc_var_t &info,char *z,const string & axes,const indexes_t & indexes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(ncid,info,z,axes,indexes);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_get_vara_template(const string &filename,const poc_var_t &info,T *z,const string & axes="",const indexes_t & indexes= indexes_t() , int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_get_vara(int,poc_var_t,T*,string,indexes_t)
/*----------------------------------------------------------------------------*/
{
  int ncid,status;
  
  status=nc_open(filename.c_str(),NC_NOWRITE,&ncid);
  if(status!=NC_NOERR) {
    
/*----------------------------------------------------------------------------*/
    if(is_grib(filename) and axes=="T"){
#ifdef HAVE_LIBGRIB_API
      status=poc_grib_get_vara(filename,info,indexes[0],z,verbose);
      
      if(status!=0 and verbose>=0)
        STDERR_BASE_LINE("poc_grib_get_vara(\""+filename+"\",(\""+info.name+"\"),) error (%d %s)\n",status,strerror(status));
      
      return status;
#else
      TRAP_ERR_RETURN(ENOEXEC,1,"Please compile with GRIB_API\n");
#endif
      }
/*----------------------------------------------------------------------------*/
    
    NC_TRAP_ERROR(return,status,verbose,"nc_open(\""+filename+"\",NC_NOWRITE,) error");
    }
  status=poc_get_vara(ncid, info, z, axes, indexes);
  if(status!=NC_NOERR && verbose) NC_CHKERR_LINE_FILE(status,"poc_get_vara() error in "+filename);
  
  nc_close(ncid);
  
  return(status);
}
#if 0
extern int poc_get_vara(const string &filename,const poc_var_t &info,T *z,const string & axes="",const indexes_t & indexes= indexes_t() , int verbose=0);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(const string &filename,const poc_var_t &info,T *z,const string & axes,const indexes_t & indexes, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(filename,info,z,axes,indexes,verbose);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(const string &filename,const poc_var_t &info,double *z,const string & axes,const indexes_t & indexes, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(filename,info,z,axes,indexes,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(const string &filename,const poc_var_t &info,float *z,const string & axes,const indexes_t & indexes, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(filename,info,z,axes,indexes,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(const string &filename,const poc_var_t &info,int8_t *z,const string & axes,const indexes_t & indexes, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(filename,info,z,axes,indexes,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_vara(const string &filename,const poc_var_t &info,char *z,const string & axes,const indexes_t & indexes, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_vara_template(filename,info,z,axes,indexes,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_get_cvara(int ncid,const poc_var_t &avar,const poc_var_t &gvar,int frame,complex<T> *z,int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// read a complex
/*----------------------------------------------------------------------------*/
{
  if(z==0) TRAP_ERR_EXIT(ENOEXEC,"NULL output buffer\n");
  
  int status,i;
  poc_data_t<T> a,g;//<amplitude and phase variables
  double angleFactor=r2d;
  const poc_att_t *gunit;
  
  a.init(avar);
  g.init(gvar);
  
  const string
    &gname=g.info.name,
    ncatted_help=
      "===> To correct this, run the following command (you can replace "+gname+" by a regexp):\n"
      "  ncatted -a \"units,"+gname+",o,char,deg\" -O [<backupfile>] <inputfile>\n";
  
  gunit=g.info.attributes.findP("units");
  if(gunit==NULL or
      ( gunit->type!=NC_CHAR and gunit->type!=NC_STRING )
      ){
    status=NC_ENOTATT;
    if(gunit==0)NC_CHKERR_LINE_FILE(status,"`units' attribute not found for phase variable "+gname);
    else if(gunit->type!=NC_CHAR)NC_CHKERR_LINE_FILE(status,"`units' attribute not a string for phase variable "+gname);
    fprintf(stderr,""+ncatted_help);
    return status;
    }
  
  switch(tolower((*gunit)[0])){
  case 'r':
    angleFactor=1.;break;
  //case '°'://warning: multi-character character constant [-Wmultichar]
  case 'd':
    angleFactor=d2r;break;
  case 'g':
    angleFactor=M_PI/200.;break;//gradians, if ever anyone uses this
  default:
    NC_TRAP_ERROR(return,NC_ENOTATT,1,"phase unit not understood: "+gunit->as_string()+".\n"+ncatted_help);
    }
  g.scale*=angleFactor;
  
  status=a.read_data(ncid,frame,verbose-1,1);
  if(status) NC_TRAP_ERROR(return,status,verbose,"poc_data_t<>::read_data() error");
  status=g.read_data(ncid,frame,verbose-1,1);
  if(status) NC_TRAP_ERROR(return,status,verbose,"poc_data_t<>::read_data() error");
  
  for(i=0;i<a.length;i++){
    if(a.data[i]==a.mask || g.data[i]==g.mask){
      z[i]=NC_FILL_COMPLEX;
      continue;
      }
    z[i]=polar<T>(a.data[i],-g.data[i]);
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_get_cvara_template(int ncid,const string &aname,const string &gname,int frame,complex<T> *z,int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// read a complex
/*----------------------------------------------------------------------------*/
{
  int status;
  poc_var_t avar,gvar;//<amplitude and phase variables
  
  status=poc_inq_var(ncid,aname,&avar);
  if(status) NC_TRAP_ERROR(return,status,verbose,"poc_inq_var(,\""+aname+"\",) error");
  status=poc_inq_var(ncid,gname,&gvar);
  if(status) NC_TRAP_ERROR(return,status,verbose,"poc_inq_var(,\""+gname+"\",) error");
  
  status=poc_get_cvara(ncid,avar,gvar,frame,z,verbose);
  if(status!=NC_NOERR){
    if(verbose)NC_CHKERR_LINE_FILE(status,"poc_get_cvara(,\""+aname+"\",\""+gname+"\",%d) error",frame);
    }
  
  return status;
}
#if 0
extern int poc_get_cvara(int ncid,const string &aname,const string &gname,int frame,complex<T> *z,int verbose=0);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_cvara(int ncid,const string &aname,const string &gname,int frame,complex<T> *z,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_cvara_template(ncid,aname,gname,frame,z,verbose);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_cvara(int ncid,const string &aname,const string &gname,int frame,complex<double> *z,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_cvara_template(ncid,aname,gname,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_cvara(int ncid,const string &aname,const string &gname,int frame,complex<float> *z,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_cvara_template(ncid,aname,gname,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<class V,typename T> int poc_get_cvara_template(const string &path,const V &aname,const V &gname,int frame,complex<T> *z,int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// read a complex
/*----------------------------------------------------------------------------*/
{
  int status;
  
  if(strrncmp(path,".s2c")==0){
    status=quoddy_loadc1(path.c_str(),NC_MAX_INT,z);
    return status;
    }
  
/*----------------------------------------------------------------------------*/
  int close_status,ncid;
  
  status=nc_open(path.c_str(),NC_NOWRITE,&ncid);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"nc_open(\""+path+"\",NC_NOWRITE,) error");
  
  status=poc_get_cvara(ncid,aname,gname,frame,z,verbose);
  if(status!=NC_NOERR){
    if(verbose)NC_CHKERR_LINE_FILE(status,"poc_get_cvara(\""+path+"\",,,%d) error",frame);
    }
  
  close_status=nc_close(ncid);
  NC_CHKERR_LINE_FILE(close_status,"nc_close() error with "+path);
  
  if(status==NC_NOERR)
    status=close_status;

  return status;
}
#if 0
extern int poc_get_cvara(const string &path,const V &aname,const V &gname,int frame,complex<T> *z,int verbose=0);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_cvara(const string &path,const V &aname,const V &gname,int frame,complex<T> *z,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_cvara_template(path,aname,gname,frame,z,verbose);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_cvara(const string &path,const string &aname,const string &gname,int frame,complex<float> *z,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_cvara_template(path,aname,gname,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_cvara(const string &path,const string &aname,const string &gname,int frame,complex<double> *z,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_cvara_template(path,aname,gname,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_cvara(const string &path,const poc_var_t &aname,const poc_var_t &gname,int frame,complex<double> *z,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_cvara_template(path,aname,gname,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_get_var_template(int ncid,const poc_var_t &info,T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_get_var(int,int,T*)
/*----------------------------------------------------------------------------*/
{
  int status;
  
  status=poc_get_var(ncid,info.id,z);
  
  return status;
}
#if 0
extern int poc_get_var(int ncid,const poc_var_t &info,T *z);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_var(int ncid,const poc_var_t &info,T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_var_template(ncid,info,z);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_var(int ncid,const poc_var_t &info,double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_var_template(ncid,info,z);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_get_var_template(const string & path,const poc_var_t &info,T *z,int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_get_var(int,poc_var_t,T*)
/*----------------------------------------------------------------------------*/
{
  int status,close_status,ncid;
  
  status=nc_open(path.c_str(),NC_NOWRITE,&ncid);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"nc_open(\""+path+"\",NC_NOWRITE,) error");
  
  status=poc_get_var(ncid,info,z);
  NC_CHKERR_BASE_LINE(status,"poc_get_var((\""+path+"\")) error");
  
  close_status=nc_close(ncid);
  NC_CHKERR_BASE_LINE(close_status,"nc_close((\""+path+"\")) error");
  
  if(status==NC_NOERR)
    status=close_status;
  
  return status;
}
#if 0
extern int poc_get_var(const string & path,const poc_var_t &info,T *z,int verbose=0);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_var(const string & path,const poc_var_t &info,T *z,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_var_template(path,info,z,verbose);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_var(const string & path,const poc_var_t &info,double *z,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_get_var_template(path,info,z,verbose);
  return status;
}
