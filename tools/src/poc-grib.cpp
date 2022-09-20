
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief grib input function definitions

See:
https://software.ecmwf.int/wiki/display/GRIB/keys_iterator.c
https://software.ecmwf.int/wiki/display/GRIB/grib_api.h+File+Reference
*/
/*----------------------------------------------------------------------------*/

#include <stdio.h>

#include "functions.h"          /* for cr, el       */
#include "poc-netcdf-data.hpp"  /* for poc_global_t */
#include "poc-grib.h"           /* for #POC_GRIB_   */
#include "poc-time.h"           /* for cnes_time()  */


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool is_grib(const string & path)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  return
    strrncasecmp(path,".grib2")==0 or
    strrncasecmp(path,".grib")==0 or
    strrncasecmp(path,".grb2")==0 or
    strrncasecmp(path,".grb")==0 or
    strrncasecmp(path,".ec")==0;
}

#ifdef HAVE_LIBGRIB_API

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_grib_get_string(grib_handle *H,const string & name,string *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// parse a grib file
/**
\param path
\param global
\param verbose
\returns NC_NOERR if success or the NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  int status;
  char *c;
  const char *cname=name.c_str();
  
  size_t size;
  status=grib_get_length(H,cname,&size);
//   updatemax(&size,2048); /* fails to prevent an error when name=="pl" */
  c=aset(size+1,'\0');
  
  status=grib_get_string(H,cname,c,&size);
  *s=c;
  delete[]c;
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_grib_inq(const string &path,poc_global_t *global,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// parse a grib file
/**
\param path
\param global
\param verbose
\returns NC_NOERR if success or the NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  int status,i,j;
  static string lastPath;
  
/*------------------------------------------------------------------------------
  make sure file has not been parsed already because this actually takes ages */
  
  poc_att_t *att;
  att=global->attributes.findP(POC_GRIB_FILE_PATH_ATT_NAME);
  if(att!=0 and att->as_string()==path){
    if(lastPath==path)
      verbose=-1;
    else
      lastPath=path;
    TRAP_ERR_RETURN(0,verbose,"skipping already parsed "+path+"\n");
    }
  
  global->dimensions.clear();
  global->attributes.clear(1);
  global->variables.clear();
  
/*------------------------------------------------------------------------------
  open file */
  FILE *F;
  
  F=fopen(path.c_str(),"r");
  if(F==0)
    TRAP_ERR_RETURN(errno,verbose,"fopen(\""+path+"\",\"r\") error (%d %s)\n",errno,strerror(errno));
  
/*------------------------------------------------------------------------------
  parse messages */
  
  struct timeval before;
  gettimeofday(&before);
  int nmsg,m;
  
  size_t *seeks=0;
  string *names=0;
  
  nmsg=999;/* reaseonable prior */
  
  seeks=new size_t[nmsg];
  names=new string[nmsg];
  
  nmsg=-1;
  
  for(m=0;m<=nmsg or nmsg==-1;m++){
    
    size_t *seekm=&seeks[m];
    *seekm=ftell(F);
    
    grib_handle *H = NULL;
    H = grib_handle_new_from_file(0, F, &status);
    if(H==0 or status!=0){
      if(status!=0){
        fclose(F);
        if(H!=0)
          grib_handle_delete(H);
        TRAP_ERR_RETURN(status,verbose,"grib_handle_new_from_file(0,(\""+path+"\"[%d]),) error (%d %s)\n",m,status,grib_get_error_message(status));
        }
      /* H==0 and status!=0 mean EOF */
      break;
      }
    
    if(m==0 and verbose>1){
/*------------------------------------------------------------------------------
      parse keys of message */
      grib_keys_iterator *K;
      const unsigned long
        key_iterator_filter_flags=
    //       GRIB_KEYS_ITERATOR_ALL_KEYS;/* skip nothing, not so bad because some keys are hidden if anything is skipped */
    //       0x7fL;/* skip everything ! */
      /*
      From ``for f in 01.out 02.out 04.out 08.out 10.out 20.out 40.out;do vd $f 00.out;done``:
        + 0x01L GRIB_KEYS_ITERATOR_SKIP_READ_ONLY skip
          - numberOfDataPoints and *NumberOf*Values that is acceptable because of grib_get_size(,"values",)
          - year but neither yearOfCentury nor dataDate (not easy to accept)
        + 0x02L GRIB_KEYS_ITERATOR_SKIP_OPTIONAL skip nothing (funny)
        + 0x04L GRIB_KEYS_ITERATOR_SKIP_EDITION_SPECIFIC skip :
          - yearOfCentury and julianDay, but not year, month, day, hour, minute, second
          - scaleValuesBy and offsetValuesBy (maybe acceptable)
          - getNumberOfValues but none of the other *NumberOf*Values
          - indicatorOfTypeOfLevel (maybe acceptable)
        + 0x08L GRIB_KEYS_ITERATOR_SKIP_CODED skip
          - month, day, hour, minute but neither year, dataDate, dataTime nor second (not easy to accept)
        + 0x10L GRIB_KEYS_ITERATOR_SKIP_COMPUTED skip *unit (BAD) and missingValue (EVEN WORSE)
        + 0x20L GRIB_KEYS_ITERATOR_SKIP_DUPLICATES skip nothing (funny)
        + 0x40L GRIB_KEYS_ITERATOR_SKIP_FUNCTION skip
          - numberOfDataPoints, not numberOfValues but all other *NumberOf*Values
          - scaleValuesBy and offsetValuesBy (maybe acceptable)
          - maximum, minimum, average, standardDeviation
      */
          GRIB_KEYS_ITERATOR_SKIP_EDITION_SPECIFIC | GRIB_KEYS_ITERATOR_SKIP_FUNCTION; /* keep 144/199 keys ... */
      unsigned long flags;
      
      for(unsigned long i=0x1L;i<0x7fL or (i&0x7fL);i<<=1){
        if(verbose<9)
          flags=key_iterator_filter_flags;
        else
          flags=i&0x7fL;
        string fileName;
        asprintf(fileName,"%02x.out",flags);
        FILE *out=fopen(fileName.c_str(),"w");
        STDERR_BASE_LINE(""+fileName+"\n");
        
        K=grib_keys_iterator_new(H,flags,0);
        if(K==0)
          TRAP_ERR_EXIT(ENOEXEC,"grib_keys_iterator_new((\""+path+"\"[%d]),,0) error\n",m);
        
        string value;
        
        while(grib_keys_iterator_next(K)){
          const char* name = grib_keys_iterator_get_name(K);
          poc_grib_get_string(H,name,&value);
          if(out==0)
            printf("%s = %s\n",name,value.c_str());
          else
            fprintf(out,"%s = %s\n",name,value.c_str());
          }
        
        grib_keys_iterator_delete(K);
        
        if(out==0)
          break;
        
        fclose(out);
        
        if(flags==key_iterator_filter_flags)
          break;
        }
      
      if(flags!=key_iterator_filter_flags)
        TRAP_ERR_EXIT(ENOEXEC,"testing different flags\n");
/*----------------------------------------------------------------------------*/
      }
    
    size_t size;
    grib_get_size(H,"values",&size);
    string *name=&names[m];
    poc_grib_get_string(H,"shortName",name);
    
    if(verbose>0){
      STDERR_BASE_LINE("at byte %lu, message %d/%d: "+*name+"[%u]",*seekm,m,nmsg,size);
      if(verbose>1)
        fprintf(stderr,"\n");
      else
        fprintf(stderr,"%s%s",el,cr);
      }
    
    poc_var_t *var;
    var=global->variables.findP(*name);
    if(var==0){
/*------------------------------------------------------------------------------
      1st time step : initialise variables */
      poc_var_t var0(*name,NC_DOUBLE);
      global->variables<<var0;
      var=global->variables.findP(*name);
      
      double mask;
      grib_get_double(H,"missingValue",&mask);
      *var<<poc_att_t(_FillValue,mask);
      
      string attString;
      
      poc_grib_get_string(H,"units",&attString);
      *var<<poc_att_t("units",attString);
      
      poc_grib_get_string(H,"name",&attString);
      *var<<poc_att_t("standard_name",attString);
      
      *var<<poc_dim_t(POC_GRIB_TIME_NAME,1);
      
      poc_grib_get_string(H,"gridType",&attString);
      *var<<poc_dim_t(attString,size);
      
      var->axes="T3";
      }
    else{
/*------------------------------------------------------------------------------
      next time steps : initialise file seeks */
      poc_dim_t *dim;
      dim=var->dimensions.findP(POC_GRIB_TIME_NAME);
      dim->len++;
      
      const size_t
        nvar=global->variables.size();
      
      if(m>2*nvar){
        size_t timeStepSize,fileSize;
        int timeStepCount,remainder;
        status=get_file_size(path,&fileSize);
        
        int m2;
        
        if(verbose>0) STDERR_BASE_LINE("checking file seeks\n");
        range_t<size_t> timeStepSizeR(-1uL,0uL);
        
        for(i=0;i<nvar;i++){
          m2=m-m%nvar+i;
          if(m2>m)m2-=nvar;
          if(verbose>1) STDERR_BASE_LINE("i=%d ; m2=%d\n",i,m2);
          timeStepSize=seeks[m2]-seeks[m2-nvar];
          timeStepSizeR<<timeStepSize;
          }
        
        if(timeStepSizeR.min!=timeStepSizeR.max){
/*----------------------------------------------------------------------------*/
          if(verbose>0) STDERR_BASE_LINE("time step sizes between %lu and %lu so far\n",timeStepSizeR.min,timeStepSizeR.max);
          
          vector<size_t> newSeeks;
          copy(&newSeeks,seeks,m+1);
          delete[]seeks;
          m2=m;
          size_t seekm2;
          if(verbose>0) STDERR_BASE_LINE("%d=%u messages so far\n",m2,newSeeks.size());
          
          while(true){
            status=grib_handle_delete(H);
            if(status!=0){
              fclose(F);
              TRAP_ERR_RETURN(status,verbose,"grib_handle_delete((\""+path+"\"[%d])) error (%d %s)\n",m2,status,grib_get_error_message(status));
              }
            
            m2++;
            seekm2=ftell(F);
            
            H = grib_handle_new_from_file(0, F, &status);
            if(H==0 or status!=0){
              if(status!=0){
                fclose(F);
                if(H!=0)
                  grib_handle_delete(H);
                TRAP_ERR_RETURN(status,verbose,"grib_handle_new_from_file(0,(\""+path+"\"[%d]),) error (%d %s)\n",m2,status,grib_get_error_message(status));
                }
              break;
              }
            
            newSeeks.push_back(seekm2);
            }
          
          seeks=copy(newSeeks);
          
          nmsg=newSeeks.size();
          timeStepCount=nmsg/nvar;
          if(verbose>0) STDERR_BASE_LINE("%d messages and %d variables so %d time steps\n",nmsg,nvar,timeStepCount);
          }
        else{
/*----------------------------------------------------------------------------*/
          if(verbose>0) STDERR_BASE_LINE("extrapolating file seeks from file size of %lu and time step size of %lu\n",fileSize,timeStepSize);
          timeStepCount=round((double)fileSize/timeStepSize);
          remainder=fileSize-timeStepSize*timeStepCount;
          nmsg=timeStepCount*nvar;
          if(remainder!=0 and verbose>0)
            STDERR_BASE_LINE("file size of %lu, time step size of %lu so %d time steps (and %d messages) with a remainder of %d\n",
              fileSize,timeStepSize,timeStepCount,nmsg,remainder);
          
          remainder=seeks[nvar]-timeStepSize;
          if(remainder!=0 and verbose>0)
            STDERR_BASE_LINE("file position of %lu and time step size of %lu so remainder of %d\n",seeks[nvar],timeStepSize,remainder);
          
          size_t *oldSeeks=seeks;
          oldSeeks[0]=remainder;
          seeks=new size_t[nmsg];
          
          for(m2=0;m2<nmsg;m2++){
            i=m2%nvar;
            j=m2/nvar;
            seeks[m2]=oldSeeks[i]+(size_t)j*timeStepSize;
            }
          
          delete[]oldSeeks;
          seeks[0]=0;
/*----------------------------------------------------------------------------*/
          }
        
        delete[]names;
        names=new string[nmsg];
        
        for(m2=0;m2<nmsg;m2++){
          i=m2%nvar;
          if(verbose>1 and m2!=0){
            j=m2/nvar;
            STDERR_BASE_LINE("at byte %lu, expecting to find message %d/%d, var[%d], time[%d]\n",seeks[m2],m2,nmsg,i,j);
            }
          names[m2]=global->variables[i].name;
          }
        
        for(i=0;i<nvar;i++){
          var=&global->variables[i];
          dim=var->dimensions.findP(POC_GRIB_TIME_NAME);
          dim->len=timeStepCount;
          }
        
        status=grib_handle_delete(H);
        if(status!=0){
          fclose(F);
          TRAP_ERR_RETURN(status,verbose,"grib_handle_delete((\""+path+"\"[%d])) error (%d %s)\n",m,status,grib_get_error_message(status));
          }
        
        break;
        }
/*----------------------------------------------------------------------------*/
      }
    
    status=grib_handle_delete(H);
    if(status!=0){
      fclose(F);
      TRAP_ERR_RETURN(status,verbose,"grib_handle_delete((\""+path+"\"[%d])) error (%d %s)\n",m,status,grib_get_error_message(status));
      }
    
    }
  
  if(m==0){/* empty file */
    status=feof(F);
    fclose(F);
    if(status!=0)
      TRAP_ERR_RETURN(status,verbose,"\""+path+"\" appears to be empty\n");
    TRAP_ERR_EXIT(ENOEXEC,"error that is not understandable!\n");
    }
  
/*------------------------------------------------------------------------------
  save file positions for speed */
  for(i=0;i<global->variables.size();i++){
    poc_var_t *var=&global->variables[i];
    const poc_dim_t
      *dim=var->dimensions.findP(POC_GRIB_TIME_NAME);
    
    if(i==0) *global<<*dim;
    
    size_t
      *varSeeks=new size_t[dim->len];
    
    j=0;
    for(m=0;m<nmsg;m++){
      const string *name=&names[m];
      if(*name!=var->name)
        continue;
      if(j<dim->len)
        varSeeks[j]=seeks[m];
      else if(verbose>=0)
        STDERR_BASE_LINE("WARNING: AT BYTE %d MESSAGE %d IS %d/%d FOR "+*name+" : FAIL SAFING...\n",seeks[m],m,j,dim->len);
      j++;
      }
    
    *var<<poc_att_t(POC_GRIB_FILE_POS_VATT_NAME,0u);
    att=var->attributes.findP(POC_GRIB_FILE_POS_VATT_NAME);
    att->init_data(varSeeks,dim->len);
    }
  
  delete[]seeks;
  
  global->dimensions[0].isunlimited=true;
  global->attributes<<poc_att_t(POC_GRIB_FILE_PATH_ATT_NAME,path);
  
/*------------------------------------------------------------------------------
  close file */
  
  status=fclose(F);
  if(status!=0)
    TRAP_ERR_RETURN(status,verbose,"fclose((\""+path+"\")) error (%d %s)\n",status,strerror(status));
  
  return status;
}


int poc_grib_get_array (grib_handle* h, const char* key, char* vals, size_t *length){
  return grib_get_string(h,key,vals,length);}
int poc_grib_get_array (grib_handle* h, const char* key, unsigned char* vals, size_t *length){
  return grib_get_bytes(h,key,vals,length);}
int poc_grib_get_array (grib_handle* h, const char* key, long int* vals, size_t *length){
  return grib_get_long_array(h,key,vals,length);}
int poc_grib_get_array (grib_handle* h, const char* key, double* vals, size_t *length){
  return grib_get_double_array(h,key,vals,length);}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_grib_get_grid_data(const string &path,const poc_var_t &var,poc_grid_data_t *gdata, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Load grid data for a variable from a GRIB file
/**
\param var variable
\param *gdata grid data
\param path GRIB file path
\param verbose
\returns 0 if success or the GRIB error code if error.

The \c *gdata will have:
  - for structured grids, \c gdata->xd and \c gdata->yd :
    + always with one dimension
    + optionally with a name
    + with data that may be fake
  - for reduced Gaussian grids :
    + \c gdata->xd set to \c "pl" values
    + both \c gdata->xd and \c gdata->yd with the same \c "pl" dimension
*/
/*----------------------------------------------------------------------------*/
{
  int status,i;
  
  poc_dim_t dim=var.dimensions.back();
  
  for(i=0;i<gdata->vc/2;i++){
    if(i==2 or i==5)
      continue;
    poc_data_t<double> *vdata=&(*gdata)[i];
    vdata->info.dimensions<<dim;
    if(i>1)
      continue;
    vdata->init();
    }
  
/*------------------------------------------------------------------------------
  open file */
  FILE *F;
  
  F=fopen(path.c_str(),"r");
  if(F==0)
    TRAP_ERR_RETURN(errno,verbose,"fopen(\""+path+"\",\"r\") error (%d %s)\n",errno,strerror(errno));
  
/*------------------------------------------------------------------------------
  seek message */
  {
  const poc_att_t *att;
  poc_global_t global;
  att=var.attributes.findP(POC_GRIB_FILE_POS_VATT_NAME);
  if(att==0){
#if 0
    fclose(F);
    status=NC_ENOTATT;
    NC_TRAP_ERROR(return,status,verbose,"file seek attribute " POC_GRIB_FILE_POS_VATT_NAME " not found");
#else
    status=poc_grib_inq(path,&global,verbose);
    const poc_var_t *var2;
    const poc_dim_t *dim2;
    for(i=0;i<global.variables.size();i++){
      var2=&global.variables[i];
      dim2=&var2->dimensions.back();
      if(dim2->len==dim.len)
        break;
      }
    att=var2->attributes.findP(POC_GRIB_FILE_POS_VATT_NAME);
#endif
    }
  const int
    position=(*att)[0];
  status=fseek(F,position,SEEK_SET);
  if(status!=0){
    fclose(F);
    TRAP_ERR_RETURN(errno,verbose,"fseek((\""+path+"\"),%d,SEEK_SET) error (%d %s)\n",position,errno,strerror(errno));
    }
  }
/*------------------------------------------------------------------------------
  open message */
  
  grib_handle *H = NULL;
  H = grib_handle_new_from_file(0, F, &status);
  if(H==0 or status!=0){
    fclose(F);
    if(H!=0)
      grib_handle_delete(H);
    TRAP_ERR_RETURN(status,verbose,"grib_handle_new_from_file(0,(\""+path+"\",\""+var.name+"\"),) error (%d %s)\n",status,grib_get_error_message(status));
    }

/*------------------------------------------------------------------------------
  CORE: call to GRIB API interface funtion */
  
  double *dummy=new double[dim.len];
  size_t size=dim.len;
  status=grib_get_data(H,gdata->yv.data,gdata->xv.data,dummy,&size);
  delete[]dummy;
  if(status!=0){
    fclose(F);
    grib_handle_delete(H);
    TRAP_ERR_RETURN(status,verbose,"grib_get_data(0,(\""+path+"\",\""+var.name+"\"),) error (%d %s)\n",status,grib_get_error_message(status));
    }
  
  /* reduced row length : rrl */
  const char *rrlName="pl";
  
  status=grib_get_size(H,rrlName,&size);
  if(status==GRIB_NOT_FOUND){
    /* regular mode 1 */
    const int
      dimCount=2,
      vcount=gdata->vc;
    
    const char
      *dimNames[dimCount]={"Ni","Nj"},
      *axes[dimCount]={"X","Y"};
    long j,sizei;
    
    for(j=0;j<dimCount;j++){
      (*gdata)[j].info.dimensions.pop_back();
      }
    
    for(i=0;i<dimCount;i++){
      const char *dimName=dimNames[i];
      
      status=grib_get_long(H,dimName,&sizei);
      if(status!=0){
        fclose(F);
        grib_handle_delete(H);
        TRAP_ERR_RETURN(status,verbose,"grib_get_long((\""+path+"\",\""+var.name+"\"),\"%s\",) error (%d %s)\n",dimName,status,grib_get_error_message(status));
        }
      
      dim.init(dimName,sizei);
      
      poc_data_t<double> *gi=&(*gdata)[i+vcount/2];
      
      gi->info.dimensions.clear();
      gi->info<<dim;
      
      fakeDimVar(gi,dim,axes[i]);
      
      for(j=0;j<dimCount;j++){
        (*gdata)[j].info.dimensions<<dim;
        }
      }
    
    }
  else if(status!=0){
    /* GRIB error */
    fclose(F);
    grib_handle_delete(H);
    TRAP_ERR_RETURN(status,verbose,"grib_get_size((\""+path+"\",\""+var.name+"\"),\"%s\",) error (%d %s)\n",rrlName,status,grib_get_error_message(status));
    }
  else{
    /* reduced Gaussian */
    const poc_dim_t rrlDim(rrlName,size);
    
    gdata->xd.info.dimensions.clear();/* fail safe */
    gdata->xd.info.dimensions<<rrlDim;
    gdata->xd.init();
    status=poc_grib_get_array(H,rrlName,gdata->xd.data,&size);
    if(status!=0){
      fclose(F);
      grib_handle_delete(H);
      TRAP_ERR_RETURN(status,verbose,"poc_grib_get_array((\""+path+"\",\""+var.name+"\"),,\"%s\",,) error (%d %s)\n",rrlName,status,grib_get_error_message(status));
      }
    
    gdata->yd.info.dimensions.clear();/* fail safe */
    gdata->yd.info.dimensions<<rrlDim;
    }
  
/*------------------------------------------------------------------------------
  close message and file */
  grib_handle_delete(H);
  
  status=fclose(F);
  if(status!=0)
    TRAP_ERR_RETURN(status,verbose,"fclose((\""+path+"\")) error (%d %s)\n",status,strerror(status));
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_grib_get_vara_template(const string &path, const poc_var_t &var,int frame, T *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
/*------------------------------------------------------------------------------
  open file */
  FILE *F;
  
  F=fopen(path.c_str(),"r");
  if(F==0)
    TRAP_ERR_RETURN(errno,verbose,"fopen(\""+path+"\",\"r\") error (%d %s)\n",errno,strerror(errno));
  
/*------------------------------------------------------------------------------
  seek message */
  
  const poc_att_t *att;
  att=var.attributes.findP(POC_GRIB_FILE_POS_VATT_NAME);
  const int64_t
    position=(*att)[frame];
  status=fseek(F,position,SEEK_SET);
  if(status!=0){
    fclose(F);
    TRAP_ERR_RETURN(errno,verbose,"fseek((\""+path+"\"),%d,SEEK_SET) error (%d %s)\n",position,errno,strerror(errno));
    }
  
/*------------------------------------------------------------------------------
  open message */
  
  grib_handle *H = NULL;
  H = grib_handle_new_from_file(0, F, &status);
  if(H==0 or status!=0){
    fclose(F);
    if(H!=0)
      grib_handle_delete(H);
    TRAP_ERR_RETURN(status,verbose,"grib_handle_new_from_file(0,(\""+path+"\",\""+var.name+"\"),) error (%d %s)\n",status,grib_get_error_message(status));
    }

/*------------------------------------------------------------------------------
  CORE */
  size_t size;
  const char *name="values";
  status=grib_get_size(H,name,&size);
  if(status!=0){
    fclose(F);
    grib_handle_delete(H);
    TRAP_ERR_RETURN(status,verbose,"grib_get_size((\""+path+"\",\"%s\"),) error (%d %s)\n",name,status,grib_get_error_message(status));
    }
  status=poc_grib_get_array(H,name,z,&size);
  if(status!=0){
    fclose(F);
    grib_handle_delete(H);
    TRAP_ERR_RETURN(status,verbose,"poc_grib_get_array((\""+path+"\",\"%s\"),) error (%d %s)\n",name,status,grib_get_error_message(status));
    }
  
/*------------------------------------------------------------------------------
  close message and file */
  grib_handle_delete(H);
  
  status=fclose(F);
  if(status!=0)
    TRAP_ERR_RETURN(status,verbose,"fclose((\""+path+"\")) error (%d %s)\n",status,strerror(status));
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T,typename T0> int poc_grib_get_vara_template2(const string &path, const poc_var_t &var,int frame, T *z, T0 dummy, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  size_t length;
  
  status=poc_get_var_length(var,&length,0,0);
  
  T0 *z0=new T0[length];
  
  status=poc_grib_get_vara_template(path,var,frame,z0,verbose);
  
  valcpy(z,z0,length);
  
  delete[]z0;
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_grib_get_vara(const string &path, const poc_var_t &var,int frame, char *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_grib_get_vara_template(path,var,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_grib_get_vara(const string &path, const poc_var_t &var,int frame, unsigned char *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_grib_get_vara_template(path,var,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_grib_get_vara(const string &path, const poc_var_t &var,int frame, long int *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_grib_get_vara_template(path,var,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_grib_get_vara(const string &path, const poc_var_t &var,int frame, double *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_grib_get_vara_template(path,var,frame,z,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_grib_get_vara(const string &path, const poc_var_t &var,int frame, signed char *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_grib_get_vara_template2(path,var,frame,z,0L,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_grib_get_vara(const string &path, const poc_var_t &var,int frame, float *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_grib_get_vara_template2(path,var,frame,z,0.,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_grib_get_vara(const string &path, const poc_var_t &var,int frame, short *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_grib_get_vara_template2(path,var,frame,z,0L,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_grib_gettime(const string &path, date_t *origine, double **time, size_t *nframes, poc_global_t *global,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  poc_global_t global_;
  if(global==0)
    global=&global_;
  
  origine->init();
  
  status=poc_grib_inq(path,global,verbose);
  
  const poc_var_t *var=&global->variables[0];
  const poc_att_t *att=var->attributes.findP(POC_GRIB_FILE_POS_VATT_NAME);
  
  *nframes=att->len;
  if(time!=0)
    *time=new double[*nframes];
  
/*------------------------------------------------------------------------------
  open file */
  FILE *F;
  
  F=fopen(path.c_str(),"r");
  if(F==0)
    TRAP_ERR_RETURN(errno,verbose,"fopen(\""+path+"\",\"r\") error (%d %s)\n",errno,strerror(errno));
  
  size_t frame;
  
  for(frame=0;frame<min(2uL,*nframes);frame++){
/*------------------------------------------------------------------------------
    seek message */
    const int64_t
      position=(*att)[frame];
    status=fseek(F,position,SEEK_SET);
    if(status!=0){
      fclose(F);
      TRAP_ERR_RETURN(errno,verbose,"fseek((\""+path+"\"),%d,SEEK_SET) error (%d %s)\n",position,errno,strerror(errno));
      }
  
/*------------------------------------------------------------------------------
    open message */
    grib_handle *H = NULL;
    H = grib_handle_new_from_file(0, F, &status);
    if(H==0 or status!=0){
      fclose(F);
      if(H!=0)
        grib_handle_delete(H);
      TRAP_ERR_RETURN(status,verbose,"grib_handle_new_from_file(0,(\""+path+"\",\""+var->name+"\"),) error (%d %s)\n",status,grib_get_error_message(status));
      }

/*------------------------------------------------------------------------------
    get date of frame */
    size_t i;
    const size_t partCount=6u+1u;
    long dateParts[partCount],
      &year=dateParts[0],&month=dateParts[1],&day=dateParts[2],
      &hour=dateParts[3],&minute=dateParts[4],&second=dateParts[5],
      &step=dateParts[6];
    const char
      *partNames[partCount]={
        "year","month","day","hour","minute","second",
        "startStep"};
    
    for(i=0;i<partCount;i++){
      const char *partName=partNames[i];
      long *datePart=&dateParts[i];
      
      status=grib_get_long(H,partName,datePart);
      if(status!=0){
        fclose(F);
        grib_handle_delete(H);
        TRAP_ERR_RETURN(status,verbose,"grib_get_long((\""+path+"\",\""+var->name+"[%u]\"),\"%s\",) error (%d %s)\n",frame,partName,status,grib_get_error_message(status));
        }
      }
    
    const date_t
      date(
        year,month,day,
        ((hour+step)*60.+minute)*60.+second );
    
    if(time!=0){
      double *frametime=&(*time)[frame];
      *frametime=cnes_time(date,'s');
      }
    
/*------------------------------------------------------------------------------
    close message */
    grib_handle_delete(H);
    }
  
  if(*nframes>2uL and time!=0){
/*------------------------------------------------------------------------------
    extrapolate */
    const double
      t0=(*time)[0],
      dt=(*time)[1]-t0;
    for(;frame<*nframes;frame++){
      double *frametime=&(*time)[frame];
      *frametime=t0+frame*dt;
      }
    
    }
  
/*------------------------------------------------------------------------------
  close file */
  status=fclose(F);
  if(status!=0)
    TRAP_ERR_RETURN(status,verbose,"fclose((\""+path+"\")) error (%d %s)\n",status,strerror(status));
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_grib_gettime(const string &path, date_t *origine, double **time, int *nframes, poc_global_t *global,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  size_t nframes0;
  status=poc_grib_gettime(path, origine, time, &nframes0, global, verbose);
  *nframes=(int)nframes0;
  return status;
}

#endif
