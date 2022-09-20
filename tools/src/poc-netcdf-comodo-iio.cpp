
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief functions definitions for COMODO (COmmunauté de MODélisation Océanique) info input/output
*/
/*----------------------------------------------------------------------------*/

#if TUGO

#include "tugo-prototypes.h"

#else

#include "version-macros.def" //for VERSION and REVISION

#include "fe-proto.h"

#endif

#include "poc-netcdf-data.hpp"
#include "poc-grib.h"            /* for poc_grib_get_time() */


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_frame_names(const string &filename,poc_list_t<poc_name_id_t> *frame_names,int verbose,string *framevarname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// gets frame names, mainly for atlases, where the frame name is the name of the wave
/**
\param filename
\param *framevarname set to the name of the variable containing the frame names
\returns 0 on full success, -1 on failure to find any frame name variable, \c -n on finding \c n>1 frame name variable
*/
/*----------------------------------------------------------------------------*/
{
  poc_global_t global;
  int i,j,n,l,status;
  const poc_var_t *var;
  poc_name_id_t name_and_id;
  char *name;
  
/*-----------------------------------------------------------------------------
  parse the file */
  status=poc_inq(filename,&global,verbose);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"poc_inq(\""+filename+"\") error");
  
  n=global.variables.size();
  j=-1;
  
  for(i=0;i<n;i++){
    var=&global.variables[i];
    
    if(var->type!=NC_CHAR ||
       var->dimensions.size()!=2 ||
       !isT(var->dimensions[0]))
      continue;
    
    if(j>=0){
      j=-2;
      }
    else if(j<-1){
      j--;
      }
    else{
      j=i;
      }
    }
  
  if(j<0)TRAP_ERR_RETURN(j,j<-1 && verbose,"*** found %d frame name variables ***\n",-j);
  
  var=&global.variables[j];
  
  if(framevarname!=NULL)*framevarname=var->name;
  
/*-----------------------------------------------------------------------------
  read the names */
  n=var->dimensions[0].len;
  l=var->dimensions[1].len;
  name=new char[l+1];
  name[l]='\0';
  
  for(i=0;i<n;i++){
    status=poc_get_vara(filename,*var,i,name,1);
    if(status!=NC_NOERR)NC_TRAP_ERROR(wexit,status,1,"poc_get_vara(\""+filename+"\",\""+var->name+"\",%d,) error",i);
    
    name_and_id.id=i;
    name_and_id.name=name;
    
    frame_names->push_back(name_and_id);
    }
  
  delete[]name;

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_frame_from_name(const string &filename,const string &name,int *frame,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// sets frame from its name
/**
\param filename
\param name name of the frame
\param *frame Set to the frame number if found. Unchanged otherwise.
\returns 0 if poc_get_frame_names() fails to find a frame variable name or if the wave is found, -1 if the wave is not found, \c -n on finding \c n>1 frame name variable
*/
/*----------------------------------------------------------------------------*/
{
  int status;
  poc_list_t<poc_name_id_t> framelist;
  const poc_name_id_t *frame_name_and_id;
  string framevarname;
  
  status=poc_get_frame_names(filename,&framelist,verbose,&framevarname);
  switch(status){
    case 0:
      status=-1;
      break;
    case -1:
      return 0;
    default:
      return status;
    }
  
  frame_name_and_id=framelist.findP(name);
  
  if(frame_name_and_id){
    *frame=frame_name_and_id->id;
    status=0;
    }
  
  if(verbose && status)STDERR_BASE_LINE("could not find frame "+name+" in "+filename+"/"+framevarname+"\n");
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool isT(const poc_dim_t &dim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  Checks if a dimension is a time dimension

  dim : reference to the dimension
  returns see isT(const char*)
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  bool unlimitedDim,timeDim;
  bool unlimited;
  
  unlimitedDim=dim.isunlimited;
  timeDim=isT(dim.name.c_str());
  
  unlimited=(unlimitedDim or timeDim);
  
  return(unlimited);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_find_timevarid(int ncid,int *tvid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,i,n;
  
  n=sizeof(timeNames)/sizeof(char *);
  
  /* check all the names of ::timeNames */
  for(i=0;i<n;i++){
    status=nc_inq_varid(ncid,timeNames[i],tvid);
    if(status==NC_NOERR)
      return status;
    if(status!=NC_ENOTVAR)NC_TRAP_ERROR(return,status,1,"nc_inq_varid() error on variable %s",timeNames[i]);
    }

  poc_global_t global;
  status=poc_inq(ncid,&global);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"poc_inq() error");
  
  n=global.variables.size();
  
  for(i=0;i<n;i++){
    if(isT(global.variables[i].name.c_str())){
      *tvid=i;
      return 0;
      }
    }
  
#warning should also check unlimited dimension variable
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_find_timevarid(const string & filename,int *tvid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int ncid,status;

  status = nc_open(filename.c_str(), 0, &ncid);
  if(status != NC_NOERR)
    return status;

  status=poc_find_timevarid(ncid,tvid);
  nc_close(ncid);

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_gettime(int ncid, date_t *origine, double **time, size_t *nframes, poc_var_t *timevar)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Get times of frames of a given NetCDF file
/** assuming the time variable name is one of ::timeNames.

Parameters
----------

:timevar: can be null

See poc_gettime(int,int,date_t*,double**,size_t*) for more help, including the one about the other parameters.
*/
/*----------------------------------------------------------------------------*/
{
  int status,tvid;
  
  status=poc_find_timevarid(ncid,&tvid);

  if(timevar!=NULL){
    status=poc_inq_var(ncid,tvid,timevar);
    if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"poc_inq_var(,%d,) error on variable %s",tvid);
    }

  ///then calls poc_gettime(int,int,date_t*,double**,size_t*)
  return poc_gettime(ncid,tvid,origine,time,nframes);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_gettime(const string & filename, date_t *origine, double **time, size_t *nframes, poc_var_t *timevar)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** Get times of frames of a given NetCDF file
\param filename path of the NetCDF file
See poc_gettime(int,date_t*,double**,size_t*,poc_var_t*) for more help, including the one about the other parameters.
*/
/*----------------------------------------------------------------------------*/
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

  status = nc_open(filename.c_str(), 0, &ncid);
  if(status != NC_NOERR)
    return status;

  status=poc_gettime(ncid, origine, time, nframes, timevar);
  nc_close(ncid);

  //if(status != NC_NOERR)
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_var_length(const poc_var_t &info, size_t *length, int *nspacedim, int *nframes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Gives length of time frame and number of space dimensions

 info
  *length    : pointer to length of time frame. If NULL, without effect.
  *nspacedim : pointer to number of space dimensions. If NULL, without effect. Default: NULL. Can be set to 0.
  *nframes   : pointer to number of frames. If NULL, without effect. Default: NULL.

  returns NC_NOERR

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int iDI;//index of info.dim
  int nspacedim_=0,nframes_=1,status=NC_NOERR;
  size_t length_=1;
  
  for(iDI=0;iDI<info.dimensions.size();iDI++){
    const poc_dim_t *dim=&info.dimensions[iDI];
    
    if(isT(*dim)){
      nframes_=dim->len;
      continue;
      }
    
    /* PATCH for 1 level pseudo-3D data */
    if(dim->len==1 && dim->name=="lev")
      continue;
    
    length_*=dim->len;
    nspacedim_++;
    }
  
  if(length!=NULL)
    *length=length_;
  
  if(nspacedim!=NULL)
    *nspacedim=nspacedim_;
  
  if(nframes!=NULL)
    *nframes=nframes_;
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_var_length(const string & path,const string & varname,size_t *length, int *nspacedim, int *nframes, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
   wrapper for poc_get_var_length(const poc_var_t &,int*,int*,int*)
   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  poc_var_t var;
  
  status=poc_inq_var(path,varname,&var);
  if(status)NC_TRAP_ERROR(return,status,verbose,"poc_inq_var(\""+path+"\",\""+varname+"\",) error");
  
  status=poc_get_var_length(var,length,nspacedim,nframes);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<class D,class V>  D findTimeDim_template(V & info,int *index=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// find time dimension
/**
\param *info variable
\param *index set to the index of the dimension
\returns a pointer to the dimension or NULL if none found
*/
/*----------------------------------------------------------------------------*/
{
  int i;
  const int n=info.dimensions.size();
  D dim;
  
  for(i=0;i<n;i++){
    dim=&info.dimensions[i];
    
    if(!isT(*dim))continue;
    
    if(index!=0)
      *index=i;
    
    return dim;
    }
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  const poc_dim_t * findTimeDim(const poc_var_t & info,int *index)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const poc_dim_t *dim;
  
  dim=findTimeDim_template<const poc_dim_t*>(info,index);
  
  return dim;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_dim_t * findTimeDim(poc_var_t & info,int *index)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_dim_t *dim;
  
  dim=findTimeDim_template<poc_dim_t*>(info,index);
  
  return dim;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_dim_t * findTimeDim(poc_global_t & glob,int *index)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_dim_t *dim;
  
  dim=findTimeDim_template<poc_dim_t*>(glob,index);
  
  return dim;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool unlimitTimeDim(poc_var_t *info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// set time dimension, if any, to unlimited
/**
\param *info variable
*/
/*----------------------------------------------------------------------------*/
{
  poc_dim_t *dim;
  dim=findTimeDim(*info);
  
  if(dim==0) return false;
  
  dim->init(0);
  
  return true;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename Tr,typename Td> inline int poc_decode_mask_template(const poc_var_t & var, Tr *scale, Tr *offset, Td *spec, Td mask_value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/** \brief gets scale, offset and mask value of a variable

\param var variable
\param *scale pointer to scale value. Defaults to 1 if not found in the attributes. Without effect if NULL
\param *offset pointer to offset value. Defaults to 0 if not found in the attributes. Without effect if NULL
\param *spec pointer to mask value. Defaults to \c mask_value if not found in the attributes. Without effect if NULL
\param mask_value mask value. Defaults to NC_FILL_CHAR NC_FILL_SHORT NC_FILL_INT NC_FILL_FLOAT or NC_FILL_DOUBLE
\returns NC_ENOTATT if one of the requested attributes is not found (which is not a big deal) or NC_NOERR otherwise.
*/
{
  int i,status=NC_NOERR;
  Tr scale_,offset_;
  Td spec_;
  int has_scale,has_offset,has_spec,has_min,has_max;
  
  has_scale=has_offset=has_spec=has_min=has_max=0;
  
  scale_=1.;
  offset_=0.;
  ///\date 2012-02-29 Damien Allain : added custom setting for default value
  spec_=mask_value;
  
  const poc_att_t *att;
  
  att=var.attributes.findP("scale_factor");
  if(att!=0){
    scale_=(*att)[0];
    has_scale=1;
    }
  
  att=var.attributes.findP("add_offset");
  if(att!=0){
    offset_=(*att)[0];
    has_offset=1;
    }
  
  /* find last attribute for mask value */
  for(i=var.attributes.size()-1;i>=0 && !(
    var.attributes[i].name=="missing_value" ||
    //var.attributes[i].name!="mask" ||
    var.attributes[i].name=="_FillValue" )
    ;i--);
  if(i>=0){
    att=&var.attributes[i];
    if(att->type!=var.type){
      STDERR_BASE_LINE_FUNC("*** \""+att->name+"\" HAS TYPE %d WHEN VAR \""+var.name+"\" HAS TYPE %d: IGNORING ***\n",att->type,var.type);
      }
    else{
      spec_=(*att)[0];
      has_spec=1;
      }
    }
  
  if(scale!=NULL){
    if(!has_scale)  status=NC_ENOTATT;/* Attribute not found */
    *scale=scale_;
    }
  if(offset!=NULL){
    if(!has_offset) status=NC_ENOTATT;/* Attribute not found */
    *offset=offset_;
    }
  if(spec!=NULL){
    if(!has_spec)   status=NC_ENOTATT;/* Attribute not found */
    *spec=spec_;
    }

  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_decode_mask(const poc_var_t &var, double *scale, double *offset, char *spec, char mask_value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_decode_mask_template(var,scale,offset,spec,mask_value);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_decode_mask(const poc_var_t &var, double *scale, double *offset, signed char *spec, signed char mask_value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_decode_mask_template(var,scale,offset,spec,mask_value);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_decode_mask(const poc_var_t &var, double *scale, double *offset, short int *spec, short int mask_value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_decode_mask_template(var,scale,offset,spec,mask_value);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_decode_mask(const poc_var_t &var, double *scale, double *offset, int *spec, int mask_value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_decode_mask_template(var,scale,offset,spec,mask_value);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_decode_mask(const poc_var_t &var, double *scale, double *offset, float *spec, float mask_value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_decode_mask_template(var,scale,offset,spec,mask_value);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_decode_mask(const poc_var_t &var, double *scale, double *offset, double *spec, double mask_value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_decode_mask_template(var,scale,offset,spec,mask_value);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_decode_mask(const poc_var_t &var, double *scale, double *offset, complex<double> *spec, const complex<double> & mask_value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  double rmask=real(mask_value);
  
  status=poc_decode_mask(var,scale,offset,&rmask,rmask);
  *spec=rmask;
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_var_length(const poc_global_t & global, const string & varname, poc_var_t *info, size_t *length, int *nspacedim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
{
  int status;
  
  status=global.variables.find(varname,info);
  if(status<0)NC_TRAP_ERROR(return,NC_ENOTVAR,1,"poc_global_t::variables.find(\""+varname+"\",) error\n");
  
  status=poc_get_var_length(*info,length,nspacedim);
  if(status!=NC_NOERR)
    NC_TRAP_ERROR(return,status,1,"poc_get_var_length error with "+varname+"\n");
  STDERR_BASE_LINE(""+varname+" : %d values",*length);
  if(nspacedim!=0)
    fprintf(stderr," in %d dimensions",*nspacedim);
  fprintf(stderr,"\n");
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_get_var_length_and_mask_template(const poc_global_t & global, const string & varname, poc_var_t *info, size_t *length, int *nspacedim, T *spec)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_global_t::variables::find(), poc_get_var_length() and poc_decode_mask()
/**
\returns the same as poc_decode_mask()
*/
/*----------------------------------------------------------------------------*/
{
  int status;
  
  status=poc_get_var_length(global,varname,info,length,nspacedim);
  
  status=poc_decode_mask(*info,NULL,NULL,spec);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_var_length_and_mask(const poc_global_t & global, const string & varname, poc_var_t *info, size_t *length, int *nspacedim, float *spec)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=poc_get_var_length_and_mask_template(global,varname,info,length,nspacedim,spec);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_inq_dimvars(const string &filename,const poc_list_t<poc_dim_t> &dims,poc_data_t<double> *dimvars,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// reads the dimension variables
/**
\param filename
\param dims list of dimensions to read the variables for
\param *dimvars array of dimension variables. Must be at least \c dims.size() allocated
\param verbose
\returns the number of variables found or a NetCDF error code
*/
/*----------------------------------------------------------------------------*/
{
  int i,j,n=0,status;
  poc_global_t global;
  
  /// <h2>Call poc_inq(string,poc_global_t*,int)</h2>
  status=poc_inq(filename,&global,verbose);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"poc_inq(\""+filename+"\",) error");
  
  for(i=0;i<dims.size();i++){
    j=global.variables.find(dims[i].name);
    if(j<0)
      continue;
    dimvars[n].init(global.variables[j]);
    dimvars[n].read_data(filename.c_str());
    n++;
    }

  return n;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void modifyFoundIndex(int *found,int i,const poc_list_t<poc_var_t> &variables,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///modify the found value
/**
\param[in,out] *found value to modify
\param i
\param global
\param global
\param verbose
\returns The index of the variable if only one complies, -1 if none comply, \c -n if \c n comply
*/
/*----------------------------------------------------------------------------*/
{
  if(*found==-1)
    *found=i;
  else if(*found>=0){
    if(verbose>=0){
      STDERR_BASE_LINE("\""+variables[*found].name+"\" complies but \n");
      STDERR_BASE_LINE("\""+variables[i].name+"\" also complies\n");
      }
    *found=-2;
    }
  else{
    (*found)--;
    if(verbose>0)STDERR_BASE_LINE("\""+variables[i].name+"\" also complies\n");
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

bool is_Dimensionless_Vertical_Coordinates(const string &sN)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///is one of #Dimensionless_Vertical_Coordinates_URL
/**
\param sN value of the \c "standard_name" attribute
*/
/*----------------------------------------------------------------------------*/
{
  bool result;
  
  result=strncmp("ocean_",sN)==0 or strncmp("atmosphere_",sN)==0;
  
  if(result)
    result=sN.find("_coordinate")!=string::npos;
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int findStandardName(const poc_global_t &global,const string *sNs,const int sNC,const poc_list_t<poc_dim_t> &dimensions,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Find the variable that has the standard name
/**
\param global
\param *sNs array[sNC] values of the \c "standard_name" attribute
\param sNC
\param verbose
\returns The index of the variable if only one complies, -1 if none comply, \c -n if \c n comply
*/
/*----------------------------------------------------------------------------*/
{
  int i,j,found=-1,status;//variable and name indexex,index of found
  const poc_att_t *att;
  
  for(i=0;i<global.variables.size();i++){/// <h1>For all variables</h1>
    
    const poc_var_t &var=global.variables[i];
    
    if(var.dimensions.size()!=1)
      continue;
    for(j=0;j<dimensions.size();j++){
      if(dimensions[j]==var.dimensions[0])
        break;
      }
    if(j>=dimensions.size())
      continue;
    
    att=var.attributes.findP("standard_name");
    
    if(!att){
      continue;
      }
    
    const char *standard_name=att->as_charp();
    status=is_Dimensionless_Vertical_Coordinates(standard_name);
    if(status==false){
      continue;
      }
    
    /// <h2>if OK, call modifyFoundIndex()</h2>
    modifyFoundIndex(&found,i,global.variables,verbose);
    }
  
  return found;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool nameS_start_or_end_with(const poc_var_t &var,const string &prefix)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Check whether any of the nameS of a variable start with a prefix
/** discarding the other names that contain \c "grid_index"
\param var variable
\param prefix prefix
*/
/*----------------------------------------------------------------------------*/
{
  const poc_att_t *att;
  int i;//<attribute index
  /** \note start and end are checked only for the variable name */
  bool endsWith=strrncasecmp(var.name,"_"+prefix)==0;
  bool startsWith=strncasecmp(prefix,var.name)==0;
  bool result=endsWith || startsWith;
  
  /** \note for the attributes only the start of the content is checked */
  if(!result){//name does not start with prefix
    for(i=0;i<var.attributes.size();i++){
      att=&var.attributes[i];
      if(strrncasecmp(att->name,"name"))//att is not a name
        continue;
      const string data(att->as_string());
      if(data.find("grid_index")!=string::npos)
        continue;
      if(!strncasecmp(prefix,data))
        break;//other name start with prefix
      }
    result=i<var.attributes.size();
    }
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int findVariableThasIs(const poc_global_t &global,const poc_list_t<poc_dim_t> &list,int nDimMax,const char *prefix,int verbose,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Find the variable that has all its dimensions in a list and that has its nameS starting with a prefix
/**
\param global
\param list of dimensions
\param nDimMax maximum number of dimensions
\param prefix
\param verbose
\param mode 0(default): by name. 1: by size.
\returns The index of the variable if only one complies, -1 if none comply, \c -n if \c n comply
*/
/*----------------------------------------------------------------------------*/
{
  int i,j,k,found=-1;//general indexes,index of found
  const poc_var_t *var;
  const int ndim=list.size();
  
  for(i=0;i<global.variables.size();i++){/// <h1>For all variables</h1>
    var=&global.variables[i];
    const int nvdim=var->dimensions.size();
    
/*---------------------------------------------------------------------*//**<h2>
    check if all its dimensions are in \c list </h2> */
    if(nvdim==0 || nvdim>nDimMax)//it has no dimension or too many
      continue;//so it is not interesting
    
    if(nvdim>ndim)//it has more dimensions than in list
      continue;//so one of them is not in list
    
    for(j=0;j<nvdim;j++){
      const poc_dim_t *dim=&var->dimensions[j];
      
      switch(mode){
      case 0:
        k=list.find(dim->name);
        break;
      case 1:
        for(k=0;k<ndim && list[k]!=*dim;k++);
        break;
      default:
        TRAP_ERR_EXIT(ENOEXEC,"%s should not have been called with mode=%d!\n",__func__,mode);
        }
      
      if(k<0 || k>=ndim)//dimension is not in list
        break;
      
      if(list[k]!=*dim)//dimension is of different size
        break;
      }
    
    if(j<nvdim)
      continue;//one of the dimensions is not in list
    
/*---------------------------------------------------------------------*//**<h2>
    call nameS_start_or_end_with()</h2> */
    if(!nameS_start_or_end_with(*var,prefix))//no name start with prefix
      continue;
    
/*---------------------------------------------------------------------*//**<h2>
    if OK, call modifyFoundIndex()</h2> */
    modifyFoundIndex(&found,i,global.variables,verbose);
    }
  
  return found;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int findVariableThasIs(const poc_global_t &global,const poc_list_t<poc_dim_t> &list,int nDimMax,const char (*prefixes)[16],int verbose,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  
  i=findVariableThasIs(global, list, nDimMax, prefixes[0], verbose, mode);
  if(i<-1){
    if(verbose>0)STDERR_BASE_LINE_FUNC("trying tighter search criterion (%s)...\n",prefixes[1]);
    i=findVariableThasIs(global, list, nDimMax, prefixes[1], verbose, mode);
    if(verbose>0 && i>=0)STDERR_BASE_LINE("\""+global.variables[i].name+"\" is the only one that complies with the tighter search criterion.\n");
    }
  else if(i==-1){
    if(verbose>0)STDERR_BASE_LINE_FUNC("trying human search criterion (%s)...\n",prefixes[2]);
    i=findVariableThasIs(global, list, nDimMax, prefixes[2], verbose, mode);
    if(verbose>0 && i>=0)STDERR_BASE_LINE("\""+global.variables[i].name+"\" is the only one that complies with the human search criterion.\n");
    }
  else if(verbose>0)STDERR_BASE_LINE_FUNC("found \""+global.variables[i].name+"\" with 1st search criterion (%s).\n",prefixes[0]);
  
  return i;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int findAnyVariableThasIs(const poc_global_t &global,int ndim,const char *suffix,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,status;
  
  for(i=0;i<global.variables.size();i++){
    const poc_var_t *var=&global.variables[i];
    const int nvdim=var->dimensions.size();
    
    if(nvdim!=ndim)
      continue;
    
    const poc_att_t *att=var->attributes.findP("standard_name");
    
    if(att==0)
      continue;
    
    const string standard_name=att->as_string();
    
    status=strrncasecmp(standard_name,suffix);
    if(status!=0)
      continue;
    
    return i;
    }
  
  return -1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int find1DVariableThasIs(const poc_global_t &global,const char *prefix,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Find the 1 dimensional variable that has its nameS starting with a prefix
/**
\param global
\param prefix
\param verbose
\returns The index of the variable if only one complies, -1 if none comply, \c -n if \c n comply
*/
/*----------------------------------------------------------------------------*/
{
  int i,found=-1;//general indexes,index of found
  const poc_var_t *var;
  
  for(i=0;i<global.variables.size();i++){/// <h1>For all variables</h1>
    var=&global.variables[i];
    
    /// <h2>check whether 1 dimension </h2>
    if(var->dimensions.size()!=1)//not 1 dimensional
      continue;
    
    /// <h2>call nameS_start_or_end_with()</h2>
    if(!nameS_start_or_end_with(*var,prefix))//no name start with prefix
      continue;
    
    /// <h2>if OK, call modifyFoundIndex()</h2>
    modifyFoundIndex(&found,i,global.variables,verbose);
    }
  
  return found;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int find2DVariableThasIs(const poc_global_t &global,int dimlen,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Find the 2 dimensional variable that has its second dimension equal to \c dimlen
/**
\param global
\param dimlen
\param verbose
\returns The index of the variable if only one complies, -1 if none comply, \c -n if \c n comply
*/
/*----------------------------------------------------------------------------*/
{
  int i,found=-1;//general indexes,index of found
  const poc_list_t<poc_dim_t> *dimensions;
  
  for(i=0;i<global.variables.size();i++){/// <h1>For all variables</h1>
    dimensions=&global.variables[i].dimensions;
    
    /// <h2>check whether 2 dimensions </h2>
    if(dimensions->size()!=2)//not 2 dimensional
      continue;
    
    /// <h2>check second dimension</h2>
    if((*dimensions)[1].len!=dimlen)//2nd dimension not equal to dimlen
      continue;
    
    /// <h2>if OK, call modifyFoundIndex()</h2>
    modifyFoundIndex(&found,i,global.variables,verbose);
    }
  
  return found;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

const poc_dim_t *findAxisDim(const poc_var_t &cvar,const poc_global_t &global,char axis,const poc_var_t *var,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Find the axis dimension of a variable from the description of a comodo-compliant NetCDF file
/**
\param cvar coordinate variable
\param global file description
\param axis letter code : X, Y, Z, T, F or whatever...
\param *var pointer to variable, used for the patch. Default: pointer to \c cvar
\param verbose
\return pointer to the axis dimension or NULL
*/
/*----------------------------------------------------------------------------*/
{
  int i;
  const poc_var_t *dv;
  const poc_att_t *att;
  const poc_dim_t *dim;
  const char *axisAtt,*c;
  char axisString[2]={0,0};
  
  axis=toupper(axis);
  
/*---------------------------------------------------------------------*//**<h1>
  For all dimensions </h1>*/
  
  for(i=0;i<cvar.dimensions.size();i++){
    dim=&cvar.dimensions[i];
    dv=global.variables.findP(dim->name);
    if(dv==NULL)continue;
    /** if the dimension has a variable */
    att=dv->attributes.findP("axis");
    if(att==NULL)continue;
    /** that has an \c "axis" attribute */
    axisAtt=att->as_charp();
    if(axisAtt==NULL || strncasecmp(axisAtt,&axis,1))
      continue;
    /** that is the same as \c axis,
    return this dimension */
    return dim;
    }
  
/*---------------------------------------------------------------------*//**<h1>
  Patched for old NetCDF files </h1>*/
  
  /// <h2>It scans the \c "axis" or \c "content" attribute, if any</h2>
  if(var==NULL)var=&cvar;
  att=var->attributes.findP("axis","content");
  if(att){
    axisAtt=att->as_charp();
    c=strchr(axisAtt, axis);
    if(c==NULL)
      c=strchr(axisAtt, tolower(axis));
    if(c!=NULL){
      if(verbose>0) STDERR_BASE_LINE("WARNING: using patch with "+att->name+" attribute\n");
      i=c-axisAtt;
      return &var->dimensions[i];
      }
    }

  /// <h2>recognise dimension name</h2>
  *axisString=axis;
  dim=global.dimensions.findP(axisString);
  if(dim==NULL){
    *axisString=tolower(axis);
    dim=global.dimensions.findP(axisString);
    }
  for(i=0;dim==NULL && i<var->dimensions.size();i++){
    dim=&var->dimensions[i];
    if(dim->name[0]==axis)
      break;
    switch(axis){
    case 'X':
      if(strncasecmp("xi_",dim->name))
        dim=NULL;
      break;
    case 'Y':
      if(strncmp("eta_",dim->name))
        dim=NULL;
      break;
    default:
      dim=NULL;
      }
    }
  if(dim){
    if(verbose>0) STDERR_BASE_LINE("WARNING: using patch with \""+dim->name+"\" dimension recognised.\n");
    return dim;
    }

  if(verbose>0) STDERR_BASE_LINE("WARNING: PATCH FAILED to find %c axis dimension for "+cvar.name+" !!!\n",axis);
  return NULL;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 void fakeDimVar(poc_data_t<double> *dv,const poc_dim_t &dim,const char *axis,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///stuff an axis variable with NC_FILL_DOUBLE
/**
\param dim dimension
\param axis letter code : X, Y, Z, T or F (or whatever, actually)
\param verbose
\returns nothing
*/
/*----------------------------------------------------------------------------*/
{
  if(verbose)STDERR_BASE_LINE("*** FAKING "+dim.name+" WITH NC_FILL_DOUBLE ***\n");
  dv->info=poc_var_t(dim.name,NC_FLOAT);
  dv->info << dim;
  dv->info << poc_att_t("axis",axis);
  dv->info << poc_att_t("production","stuffed around " __LINE_FILE_PACKAGE_REVISION);
  dv->init("");
  aset(dv->data,dim.len,NC_FILL_DOUBLE);
}
