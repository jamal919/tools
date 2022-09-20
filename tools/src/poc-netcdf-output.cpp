
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


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename Td,typename Ts> Td *poc_unscale_data_template(const Td *data, int length, Ts scale, Ts offset, Td mask=NC_FILL_DOUBLE)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// un-scales the data for writing
/*----------------------------------------------------------------------------*/
{
  int i;
  Td *u=new Td[length];
  
  scale=1./scale;
  
  for(i=0;i<length;i++){
    if(data[i]==mask)
      u[i]=data[i];
    else
      u[i]=(data[i]-offset)*scale;
    }
  
  return u;
}
#if 0
extern Td *poc_unscale_data(const Td *data, int length, Ts scale, Ts offset, Td mask=NC_FILL_DOUBLE);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  Td *poc_unscale_data(const Td *data, int length, Ts scale, Ts offset, Td mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  Td *result;
  result=poc_unscale_data_template(data,length,scale,offset,mask);
  return result;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *poc_unscale_data(const double *data, int length, double scale, double offset, double mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double *result;
  result=poc_unscale_data_template(data,length,scale,offset,mask);
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float *poc_unscale_data(const float *data, int length, double scale, double offset, float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float *result;
  result=poc_unscale_data_template(data,length,scale,offset,mask);
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *poc_unscale_data(const int *data, int length, double scale, double offset, int mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int *result;
  result=poc_unscale_data_template(data,length,scale,offset,mask);
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *poc_unscale_data(const char *data, int length, double scale, double offset, char mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  char *result;
  result=poc_unscale_data_template(data,length,scale,offset,mask);
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  signed char *poc_unscale_data(const signed char *data, int length, double scale, double offset, signed char mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  signed char *result;
  result=poc_unscale_data_template(data,length,scale,offset,mask);
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int poc_put_vara_template(int ncid,const poc_var_t &info, int frame,const T *z, int unscale=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_put_vara(int,int,const size_t[],const size_t[],T*)
/**
\param unscale Set to 1 if you want to unscale \c z before writing.
*/
/*----------------------------------------------------------------------------*/
{
  int status,varid,i,len=1;
  const int dimC=info.dimensions.size();
  poc_var_t savedInfo;
  double scale,offset;
  T mask;
  
  status=nc_inq_varid(ncid,info.name.c_str(),&varid);
  if(status==NC_ENOTVAR){
    poc_def_var(ncid,info,&varid);
    savedInfo=info;
    }
  else if(unscale){
    poc_inq_var(ncid,varid,&savedInfo);
    }
  if(unscale){/// <h1>if \c unscale is set to 1</h1>
    ///get \c scale and \c offset values with poc_decode_mask()
    poc_decode_mask(savedInfo,&scale,&offset,&mask);
    }

  size_t *start,*count;
  
  start=new size_t[dimC];
  count=new size_t[dimC];
  
  for(i=0;i<dimC;i++){
    const poc_dim_t *dim=&info.dimensions[i];
    if( isT(*dim) and
        frame!=-2 ){
      start[i]=frame;
      count[i]=1;
      }
    else{
      start[i]=0;
      count[i]=dim->len;
      len*=count[i];
      }
    }
  
  if(unscale && (scale!=1. || offset!=0.)){/// <h1>if \c unscale is set to 1 and if \c scale and \c offset values make it necessary to do so</h1>
    ///safely unscale \c z
    T *u=poc_unscale_data(z,len,scale,offset,mask);
    status=poc_put_vara(ncid,varid,start,count,u);
    delete[]u;
    }
  else{
    status=poc_put_vara(ncid,varid,start,count,z);
    }
  delete[]start;
  delete[]count;
  return status;
}
#if 0
extern int poc_put_vara(int ncid,const poc_var_t &info, int frame,const T *z, int unscale=0);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_vara(int ncid,const poc_var_t &info, int frame,const T *z, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_vara_template(ncid,info,frame,z,unscale);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_vara(int ncid,const poc_var_t &info, int frame,const double *z, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_vara_template(ncid,info,frame,z,unscale);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_vara(int ncid,const poc_var_t &info, int frame,const float *z, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_vara_template(ncid,info,frame,z,unscale);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_vara(int ncid,const poc_var_t &info, int frame,const int *z, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_vara_template(ncid,info,frame,z,unscale);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_vara(int ncid,const poc_var_t &info, int frame,const char *z, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_vara_template(ncid,info,frame,z,unscale);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_vara(int ncid,const poc_var_t &info, int frame,const signed char *z, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_vara_template(ncid,info,frame,z,unscale);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int poc_put_vara_template(const string &path,const poc_var_t &info, int frame,const T *z, int verbose=0, int unscale=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_put_vara(int,poc_var_t,size_t,T*)
/*----------------------------------------------------------------------------*/
{
  int status,close_status,ncid;
  
/*------------------------------------------------------------------------------
  opens the file */
  status=poc_open_and_may_create(path,&ncid,verbose);
  NC_TRAP_ERROR(return,status,verbose,"poc_open_and_may_create(\""+path+"\",,) error");
  
/*------------------------------------------------------------------------------
  calls poc_put_vara(int,poc_var_t,size_t,T*) */
  status=poc_put_vara(ncid, info, frame, z,unscale);
  if(verbose) NC_CHKERR_LINE_FILE(status,"poc_put_vara() error with "+info.name+" in "+path);
  
/*------------------------------------------------------------------------------
  closes the file */
  close_status=nc_close(ncid);
  if(verbose) NC_CHKERR_LINE_FILE(close_status,"nc_close() error with "+path);
  
  if(status==NC_NOERR)
    status=close_status;
  
  return status;
}
#if 0
extern int poc_put_vara(const string &path,const poc_var_t &info, int frame,const T *z, int verbose=0, int unscale=0);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_vara(const string &path,const poc_var_t &info, int frame,const T *z, int verbose, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_vara_template(path,info,frame,z,verbose,unscale);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_vara(const string &path,const poc_var_t &info, int frame,const double *z, int verbose, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_vara_template(path,info,frame,z,verbose,unscale);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_vara(const string &path,const poc_var_t &info, int frame,const float *z, int verbose, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_vara_template(path,info,frame,z,verbose,unscale);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_vara(const string &path,const poc_var_t &info, int frame,const int *z, int verbose, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_vara_template(path,info,frame,z,verbose,unscale);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_vara(const string &path,const poc_var_t &info, int frame,const char *z, int verbose, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_vara_template(path,info,frame,z,verbose,unscale);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_vara(const string &path,const poc_var_t &info, int frame,const signed char *z, int verbose, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_vara_template(path,info,frame,z,verbose,unscale);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename F,typename T> int poc_put_vara_template(const F &file,const string &var, int frame, const T *z, int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// saves a complex
/*----------------------------------------------------------------------------*/
{
  poc_var_t info;
  int status;
  
  status=poc_inq_var(file,var,&info);
  if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,verbose,"poc_inq_var(,"+var+") error");
  
  status=poc_put_vara(file,info,frame,z);
  
  return status;
}
#if 0
extern int poc_put_vara(F file,const string &var, int frame, const T *z, int verbose=0);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_vara(F file,const string &var, int frame, const T *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_vara_template(file,var,frame,z,verbose);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_vara(const string &file,const string &var, int frame, const double *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_vara_template(file,var,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_vara(const string &file,const string &var, int frame, const float *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_vara_template(file,var,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_vara(const string &file,const string &var, int frame, const int *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_vara_template(file,var,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int poc_put_cvara_template(int ncid,const poc_var_t &ainfo,const poc_var_t &ginfo, size_t frame, const complex<T> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  
  saves a complex buffer
  
------------------------------------------------------------------------------*/
{
  int status,avarid,gvarid,i;
  size_t *start,*count,length=1,nlimited=0;
  T *a,*g;
  
  status=nc_inq_varid(ncid,ainfo.name.c_str(),&avarid);
  if(status==NC_ENOTVAR){
    poc_def_var(ncid,ainfo,&avarid);
    }
  status=nc_inq_varid(ncid,ginfo.name.c_str(),&gvarid);
  if(status==NC_ENOTVAR){
    poc_def_var(ncid,ginfo,&gvarid);
    }
  
  status=nc_inq_attid(ncid,gvarid,"units",0);
  if(status==NC_ENOTATT){
    status=poc_def_att(ncid,gvarid,poc_att_t("units","degrees"));
    if(status!=0) NC_CHKERR_BASE_LINE(status,"poc_def_att(,,poc_att_t(\"units\",)) error");
    }
  
  T amask,gmask;
  status=poc_decode_mask(ainfo,NULL,NULL,&amask);
  status=poc_decode_mask(ginfo,NULL,NULL,&gmask);

  if(ainfo.dimensions.size()!=ginfo.dimensions.size()){
    return NC_EDIMSIZE; /* Invalid dimension size */
    }
  
  start=new size_t[ainfo.dimensions.size()];
  count=new size_t[ainfo.dimensions.size()];
  
  for(i=0;i<ainfo.dimensions.size();i++){
    if(isT(ainfo.dimensions[i])){
      start[i]=frame;
      count[i]=1;
      }
    else{
      if(ainfo.dimensions[i].len!=ginfo.dimensions[i].len){
        status=NC_EDIMSIZE; /* Invalid dimension size */
        goto deleteStartAndCount;
        }
      start[i]=0;
      count[i]=ainfo.dimensions[i].len;
      length*=count[i];
      nlimited++;
      }
    }
  if(nlimited==0){
    length=0;
    }
  
  a=new T[length];
  g=new T[length];
  for(i=0;i<length;i++){
    if(real(z[i])==amask){
      a[i]=amask;
      g[i]=gmask;
      continue;
      }
    a[i]=abs(z[i]);
    g[i]=-arg(z[i])*r2d;
    }
  
  status=poc_put_vara(ncid,avarid,start,count,a);
  status=poc_put_vara(ncid,gvarid,start,count,g);
  
  delete[]a;
  delete[]g;
  
deleteStartAndCount:

  delete[]start;
  delete[]count;
  return status;
}
#if 0
extern int poc_put_cvara(int ncid,const poc_var_t &ainfo,const poc_var_t &ginfo, size_t frame, const complex<T> *z);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_cvara(int ncid,const poc_var_t &ainfo,const poc_var_t &ginfo, size_t frame, const complex<T> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_cvara_template(ncid,ainfo,ginfo,frame,z);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_cvara(int ncid,const poc_var_t &ainfo,const poc_var_t &ginfo, size_t frame, const complex<double> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_cvara_template(ncid,ainfo,ginfo,frame,z);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_cvara(int ncid,const poc_var_t &ainfo,const poc_var_t &ginfo, size_t frame, const complex<float> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_cvara_template(ncid,ainfo,ginfo,frame,z);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int poc_put_cvara_template(const string &path, const poc_var_t &ainfo,const poc_var_t &ginfo, int frame, const complex<T> *z, int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// saves a complex
/*----------------------------------------------------------------------------*/
{
  int status,close_status,ncid;
  
  ///-opens the file
  status=poc_open_and_may_create(path,&ncid,verbose);
  NC_TRAP_ERROR(return,status,verbose,"poc_open_and_may_create(\""+path+"\",,) error");
  
  ///-calls poc_put_vara(int,poc_var_t,poc_var_t,size_t,T*)
  status=poc_put_cvara(ncid, ainfo, ginfo, frame, z);
  if(verbose) NC_CHKERR_LINE_FILE(status,"poc_put_vara() error with "+ainfo.name+" and "+ginfo.name+" in "+path);
  
  ///-closes the file
  close_status=nc_close(ncid);
  if(verbose) NC_CHKERR_LINE_FILE(close_status,"nc_close() error with "+path);
  
  if(status==NC_NOERR)
    status=close_status;

  return status;
}
#if 0
extern int poc_put_cvara(const string &path, const poc_var_t &ainfo,const poc_var_t &ginfo, int frame, const complex<T> *z, int verbose=0);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_cvara(const string &path, const poc_var_t &ainfo,const poc_var_t &ginfo, int frame, const complex<T> *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_cvara_template(path,ainfo,ginfo,frame,z,verbose);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_cvara(const string & path, const poc_var_t &ainfo,const poc_var_t &ginfo, int frame, const complex<double> *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_cvara_template(path,ainfo,ginfo,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_cvara(const string & path, const poc_var_t &ainfo,const poc_var_t &ginfo, int frame, const complex<float> *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_cvara_template(path,ainfo,ginfo,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename F,typename T> int poc_put_cvara_template(const F &file,const string &aname,const string &gname, int frame, const complex<T> *z, int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// saves a complex
/*----------------------------------------------------------------------------*/
{
  poc_var_t ainfo,ginfo;
  int status;
  
  status=poc_inq_var(file,aname,&ainfo);
  if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,verbose,"poc_inq_var(,\""+aname+"\") error");
  status=poc_inq_var(file,gname,&ginfo);
  if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,verbose,"poc_inq_var(,\""+gname+"\") error");
  
  status=poc_put_cvara(file,ainfo,ginfo,frame,z);
  
  return status;
}
#if 0
extern int poc_put_cvara(const F &file,const string &aname,const string &gname, int frame, const complex<T> *z, int verbose=0);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_cvara(const F &file,const string &aname,const string &gname, int frame, const complex<T> *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_cvara_template(file,aname,gname,frame,z,verbose);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_cvara(const string &file,const string &aname,const string &gname, int frame, const complex<double> *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_cvara_template(file,aname,gname,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_cvara(const string &file,const string &aname,const string &gname, int frame, const complex<float> *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_cvara_template(file,aname,gname,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int poc_put_cvara_template(const string & path, const poc_var_t & infoTemplate, int frame, const complex<T> *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  poc_var_t ainfo=infoTemplate,ginfo=infoTemplate;
  
  string amp_standard_name,pha_standard_name,units;
  const poc_att_t *att;
  
  att=infoTemplate.attributes.findP("standard_name");
  if(att!=0){
    size_t pos;
    amp_standard_name=att->as_string();
    pha_standard_name=amp_standard_name;
    pos=amp_standard_name.find("_due_to_");
    if(pos==string::npos)
      pos=amp_standard_name.length();
    amp_standard_name.insert(pos,"_amplitude");
    pha_standard_name.insert(pos,"_phase_lag");
    }
  
  att=infoTemplate.attributes.findP("units");
  if(att!=0)
    units=att->as_string();
  
  ainfo.init(infoTemplate.name+"_a",infoTemplate.type,amp_standard_name,units);
  ginfo.init(infoTemplate.name+"_G",infoTemplate.type,pha_standard_name,"degrees");
  
  status=poc_put_cvara_template(path,ainfo,ginfo,frame,z,verbose);
  
  return status;
}
#if 0
extern int poc_put_cvara(const string & path, const poc_var_t & infoTemplate, int frame, const complex<T> *z, int verbose);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_cvara(const string & path, const poc_var_t & infoTemplate, int frame, const complex<T> *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_cvara_template(path,infoTemplate,frame,z,verbose);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_cvara(const string & path, const poc_var_t & infoTemplate, int frame, const complex<double> *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_cvara_template(path,infoTemplate,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_cvara(const string & path, const poc_var_t & infoTemplate, int frame, const complex<float> *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_cvara_template(path,infoTemplate,frame,z,verbose);
  return status;
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int poc_put_var_template(int ncid,const poc_var_t &info, const T *z, int unscale=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_put_var(int,int,T*)
/**
\param unscale Set to 1 if you want to unscale \c z before writing.
*/
/*----------------------------------------------------------------------------*/
{
  int status,varid;
  poc_var_t savedInfo;
  double scale,offset;
  T mask;
  
  status=nc_inq_varid(ncid,info.name.c_str(),&varid);
  if(status==NC_ENOTVAR){
    poc_def_var(ncid,info,&varid);
    savedInfo=info;
    }
  else if(unscale){
    poc_inq_var(ncid,varid,&savedInfo);
    }
  if(unscale){/// <h1>if \c unscale is set to 1</h1>
    ///get \c scale and \c offset values with poc_decode_mask()
    poc_decode_mask(savedInfo,&scale,&offset,&mask);
    }

  if(unscale && (scale!=1. || offset!=0.)){/// <h1>if \c unscale is set to 1 and if \c scale and \c offset values make it necessary to do so</h1>
    T *u;
    int i,len=1;
    vector<poc_dim_t> &dimensions=savedInfo.dimensions;
    
    for(i=0;i<dimensions.size();i++){
      len*=dimensions[i].len;
      }
    
    ///safely unscale \c z
    u=poc_unscale_data(z,len,scale,offset,mask);
    status=poc_put_var(ncid,varid,u);
    delete[]u;
    }
  else{
    status=poc_put_var(ncid,varid,z);
    }
  
  return status;
}
#if 0
extern int poc_put_var(int ncid,const poc_var_t &info, const T *z, int unscale=0);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_var(int ncid,const poc_var_t &info, const T *z, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_var_template(ncid,info,z,unscale);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_var(int ncid,const poc_var_t &info, const double *z, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_var_template(ncid,info,z,unscale);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_var(int ncid,const poc_var_t &info, const float *z, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_var_template(ncid,info,z,unscale);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_var(int ncid,const poc_var_t &info, const int *z, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_var_template(ncid,info,z,unscale);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int poc_put_var_template(const string &path, const poc_var_t &info, const T *z, int verbose=0, int unscale=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_put_var(int,poc_var_t,size_t,T*)
/*----------------------------------------------------------------------------*/
{
  int status,close_status,ncid;
  
  status=poc_open_and_may_create(path,&ncid,verbose);
  NC_TRAP_ERROR(return,status,verbose,"poc_open_and_may_create(\""+path+"\",,) error");
  
  ///-calls poc_put_var(int,poc_var_t,size_t,T*)
  status=poc_put_var(ncid, info, z,unscale);
  if(verbose)NC_CHKERR_LINE_FILE(status,"poc_put_var() error with "+info.name+" in "+path);
  
  ///-closes the file
  close_status=nc_close(ncid);
  if(verbose)NC_CHKERR_LINE_FILE(close_status,"nc_close() error with "+path);
  
  if(status==NC_NOERR)
    status=close_status;

  return status;
}
#if 0
extern int poc_put_var(const string &path, const poc_var_t &info, const T *z, int verbose=0, int unscale=0);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_var(const string &path, const poc_var_t &info, const T *z, int verbose, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_var_template(path,info,z,verbose,unscale);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_var(const string &path, const poc_var_t &info, const double *z, int verbose, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_var_template(path,info,z,verbose,unscale);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_var(const string &path, const poc_var_t &info, const float *z, int verbose, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_var_template(path,info,z,verbose,unscale);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_var(const string &path, const poc_var_t &info, const int *z, int verbose, int unscale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_var_template(path,info,z,verbose,unscale);
  return status;
}

