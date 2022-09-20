
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Makes statistics on NetCDF outputs.

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

#include "tides.h"
#include "functions.h" //safely includes omp.h
#include "matrix.h"
#include "poc-netcdf-data.hpp"
#include "poc-time.h"
#include "map.h"

struct timeval stv;///<start timeval, for progression
///These global variables are used to count the amount of time this programme took for I/O tasks
struct timeval before;
double rt=0.,rt2=0.,rrt=0.,wt=0.;// read, 2nd read, redundant read and write times
double ft=0.;// filtering time

string cmd;///<command line
size_t nt;///<number of samples

const char *argv0=NULL;
void print_help(const char *prog_name);

poc_att_t period_att;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int atlas_(int nvars, double **average, double **stddev, double **minval, double **maxval, const bool *mask, poc_grid_data_t *coordinates, const poc_var_t *info, const char *outputDir=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,i,nsd,status;
  size_t length;
  poc_global_t global("constructed around " __LINE_FILE_PACKAGE_REVISION);
  poc_var_t aveV[nvars],stdV[nvars],minV[nvars],maxV[nvars];
  poc_att_t *att;
  string coordinateNames,s;
  double defaultSpec,spec;//fill value
  
  for(k=0;k<nvars;k++){
///----------------------------------------------------------------------------
/// It initialises the file with the coordinates variables
    status=poc_get_var_length(info[k],&length,&nsd);
    for(i=0;i<coordinates[k].vc;i++){
      if(coordinates[k][i].info.name.length()==0){continue;}
      global.variables << coordinates[k][i].info;
      }
    coordinateNames="";
    for(i=poc_grid_data_t::vc/2-1;i>=0;i--){
      if(coordinateNames!="")
        coordinateNames+=" ";
      coordinateNames+=coordinates[k][i].info.name;
      }
    
///----------------------------------------------------------------------------
/// and with the average variable
    aveV[k].init(info[k].name+"_mean",NC_DOUBLE);
    stdV[k].init(info[k].name+"_stddev",NC_DOUBLE);
    minV[k].init(info[k].name+"_min",info[k].type);
    maxV[k].init(info[k].name+"_max",info[k].type);
    
    /* dimensions */
    for(i=0;i<info[k].dimensions.size();i++){
      if(isT(info[k].dimensions[i]))continue;
      aveV[k] << info[k].dimensions[i];
      }
    stdV[k].dimensions=aveV[k].dimensions;
    minV[k].dimensions=aveV[k].dimensions;
    maxV[k].dimensions=aveV[k].dimensions;
    
    /* attributes */
    if(coordinateNames!="")
      aveV[k] << poc_att_t("coordinates",coordinateNames);
    
    minV[k].attributes=aveV[k].attributes;
    maxV[k].attributes=aveV[k].attributes;
    
    switch(aveV[k].type){
    #define unionDataType(t) case NC_ ## t:defaultSpec=NC_FILL_ ## t;break
      unionDataType(BYTE);
      unionDataType(SHORT);
      unionDataType(INT);
      unionDataType(FLOAT);
      unionDataType(DOUBLE);
    #ifdef NC_UBYTE
      unionDataType(UBYTE);
      unionDataType(USHORT);
      unionDataType(UINT);
      unionDataType(INT64);
      unionDataType(UINT64);
      default:
        status=NC_EBADTYPID;
    #else
      default:
        status=NC_EAXISTYPE;
      #warning UPGRADE YOUR NetCDF : NC_UBYTE NC_USHORT NC_UINT NC_INT64 NC_UINT64 are undefined
    #endif
      NC_CHKERR_LINE_FILE(status,"type error with type %d\n",info[k].type);exit(status);
    #undef unionDataType
      }
    status=poc_decode_mask(info[k],NULL,NULL,&spec,defaultSpec);
    switch(aveV[k].type){
    #define unionDataType(t) aveV[k] << poc_att_t("_FillValue",(t)spec);break
      case NC_BYTE:    unionDataType(int8_t);
      case NC_SHORT:   unionDataType(int16_t);
      case NC_INT:     unionDataType(int32_t);
      case NC_FLOAT:   unionDataType(float);
      case NC_DOUBLE:  aveV[k] << poc_att_t("_FillValue",spec); break;
    #ifdef NC_UBYTE
      case NC_UBYTE:   unionDataType(uint8_t);
      case NC_USHORT:  unionDataType(uint16_t);
      case NC_UINT:    unionDataType(uint32_t);
      case NC_INT64:   unionDataType(int64_t);
      case NC_UINT64:  unionDataType(uint64_t);
    #else
      #warning UPGRADE YOUR NetCDF : NC_UBYTE NC_USHORT NC_UINT NC_INT64 NC_UINT64 are undefined
    #endif
    #undef unionDataType
      }
    
    stdV[k].attributes=aveV[k].attributes;
    
    for(i=0;i<info[k].attributes.size();i++){
      att=(poc_att_t*)&info[k].attributes[i];
      
      if(att->name=="long_name"){
        aveV[k] << poc_att_t(att->name,"average of "+att->as_string());
        stdV[k] << poc_att_t(att->name,"standard deviation of "+att->as_string());
        minV[k] << poc_att_t(att->name,"minimum of "+att->as_string());
        maxV[k] << poc_att_t(att->name,"maximum of "+att->as_string());
        continue;
        }
      
      if(att->name=="standard_name"){
        aveV[k] << poc_att_t(att->name,"standard_deviation_of_"+att->as_string());
        stdV[k] << poc_att_t(att->name,"mean_of_"+att->as_string());
        minV[k] << poc_att_t(att->name,"minimum_of_"+att->as_string());
        maxV[k] << poc_att_t(att->name,"maximum_of_"+att->as_string());
        continue;
        }
      
      if(att->name=="coordinates" || //already specified above
        att->name=="associate" ||
        att->name=="missing_value"){
        continue;
        }
        
      if(att->name=="content") {
        string s=att->as_string();
        printf("%s %s\n",att->name.c_str(),s.c_str());
        if(s[0]=='T') {
          s.erase(0, 1);
          printf("%s %s\n",att->name.c_str(),s.c_str());
          char *value=strdup(s.c_str());
          int len=s.size();
          att->init_data(value, len);
          }
        }
      
      minV[k] << *att;
      maxV[k] << *att;
      
      if(att->name==_FillValue) //already specified above for aveV and stdV
        continue;
      
      aveV[k] << *att;
      
      if(att->name=="add_offset")
        continue;
      
      stdV[k] << *att;
      }
    
    global.variables << aveV[k];
    global.variables << stdV[k];
    global.variables << minV[k];
    global.variables << maxV[k];
    }
  
  global << poc_att_t("history",cmd);
  if(period_att.len>0)
    global << period_att;

  ///it creates the file
  if(outputDir!=NULL)s=string(outputDir)+"/";
  s+=info[0].name+"-stats.nc";
  STDERR_BASE_LINE("(over)writing %s ...",s.c_str());
  unlink(s.c_str());
  status=poc_create(s,global,1,NC_NETCDF4);
  if(status!=NC_NOERR)NC_TRAP_ERROR(wexit,status,1,"poc_create(\""+s+"\",) error\n");
  
  for(k=0;k<nvars;k++){
    ///it writes the coordinates
    for(i=0;i<coordinates[k].vc;i++){
      if(coordinates[k][i].info.name.length()==0){continue;}
      fprintf(stderr," "+coordinates[k][i].info.name+",");
      status=coordinates[k][i].write_data(s.c_str(),0,0,1);//descaling
      NC_CHKERR_LINE_FILE(status,"poc_data_t::write_data(\"%s\") error with %s on %s",s.c_str(),coordinates[k][i].info.name.c_str());
      }
    
    ///and it writes the data
    fprintf(stderr," "+aveV[k].name+",");
    status=poc_put_vara(s.c_str(),aveV[k],0,average[k]);
    NC_CHKERR_LINE_FILE(status,"poc_put_vara error with "+aveV[k].name+" on "+s);
    fprintf(stderr," "+stdV[k].name+",");
    status=poc_put_vara(s.c_str(),stdV[k],0,stddev[k]);
    NC_CHKERR_LINE_FILE(status,"poc_put_vara error with "+stdV[k].name+" on "+s);
    fprintf(stderr," "+minV[k].name+",");
    status=poc_put_vara(s.c_str(),minV[k],0,minval[k]);
    NC_CHKERR_LINE_FILE(status,"poc_put_vara error with "+minV[k].name+" on "+s);
    fprintf(stderr," "+maxV[k].name+",");
    status=poc_put_vara(s.c_str(),maxV[k],0,maxval[k]);
    NC_CHKERR_LINE_FILE(status,"poc_put_vara error with "+maxV[k].name+" on "+s);
    fprintf(stderr," done.\n");
    }
  #undef T
  #undef NC_T

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int scalar_statistics(const char **varnames, int nvars, vector<string> filelist, char *gridfile, size_t nt,double *ts,const double origd=NAN,int handleRepeatedFrames=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// RMS of the filtered data
/**
\param **varnames list of variable names
\param nvars number of variables
\param filelist list of files
\param *gridfile path to grid file
\param *controlfile path to list of control points
\param nt number of frames
\param ts times of the frames in second
\param atlas whether an atlas is done
\param origd if specified, default time origin if not in the files
\param handleRepeatedFrames if 0: return NULL if some of the time frames are repeated. If 1:take the last one. If -1:take the first one.
\returns the RMS of the filtered data, or NULL if error

\date 2011-09-19 Damien Allain : added frame times check
*/
/*----------------------------------------------------------------------------*/
{
  double **inBuf,**sum,**sum2,**minval,**maxval;
  int j,k,i;//variable and file indexes
  int status;
  int ifi,dfi;//input frame index,duplicate frame index
  poc_global_t global;
  poc_var_t *info;
  double *spec;
  bool *mask=NULL;
  size_t *nvalues,ntotal;
  int *n_spacedims;
  
  if(nvars<=0)return NC_EINVAL;
//   if(nvars!=1){
//     CHKERR_LINE_FILE(ENOEXEC,"no more than one variable accepted yet because of historical bug");
//     TRAP_ERR_EXIT(ENOEXEC,"exiting\n");/* Exec format error */
//     }
  nvalues=new size_t[nvars];
  n_spacedims=new int[nvars];
  spec=new double[nvars];
  
/*------------------------------------------------------------------------------
  handle repeated frames */
  bool *skipFrame=NULL;
  status=getSkippedFrames(nt,ts,handleRepeatedFrames,&skipFrame);
  if(status!=0){
    printf("\n");
    #warning clean this up with all the tfi,dfi
    STDOUT_BASE_LINE(" *** simultaneous time frames: give either --take-first or --take-last ***\n");
    print_help(argv0);
    printf("\n");
    return NC_EINVAL;
    }

  if(ts!=0){
    double startd,finald;
    minpos(ts,nt,&startd);
    maxpos(ts,nt,&finald);
    period_att.init_copy("analysis_period",(string)"from "+poctime_sdate_cnes(startd,'s')+" to "+poctime_sdate_cnes(finald,'s'));
    }
  
/*---------------------------------------------------------------------*//**<h1>
  gets variable dimensions </h1>*/
  info=new poc_var_t[nvars];
  gettimeofday(&before);
  status=poc_inq(filelist[0],&global,0);
  rrt+=difftime(before);
  ntotal=0;
  for(k=0;k<nvars;k++){
    status=global.variables.find(varnames[k],&info[k]);
    if(status<0){
      status=NC_ENOTVAR;NC_CHKERR_LINE_FILE(status,"poc_global_t::variables.find(\"%s\",) error in %s",varnames[k],filelist[0].c_str());return status;
      }
    status=poc_get_var_length(info[k],&nvalues[k],&n_spacedims[k]);
    if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,1,"poc_get_var_length error with %s in %s",varnames[k],filelist[0].c_str());
    STDERR_BASE_LINE("%s : %d values in %d dimensions\n",varnames[k],nvalues[k],n_spacedims[k]);
    ntotal+=nvalues[k];
    status=poc_decode_mask(info[k],NULL,NULL,&spec[k]);
    if(status==NC_ENOTATT) NC_CHKERR_LINE_FILE(status,"Use default mask value of %g for %s in %s",spec[k],varnames[k],filelist[0].c_str());//warning only
    }
  
  ///to allocate buffers
  
  inBuf=new double*[nvars];
  
  sum   = new double*[nvars];
  sum2  = new double*[nvars];
  minval= new double*[nvars];
  maxval= new double*[nvars];
  
  for(k=0;k<nvars;k++){
    exitIfNull(inBuf[k]  = new double[nvalues[k]]);
    exitIfNull(sum[k]    = new double[nvalues[k]]);
    exitIfNull(sum2[k]   = new double[nvalues[k]]);
    exitIfNull(minval[k] = new double[nvalues[k]]);
    exitIfNull(maxval[k] = new double[nvalues[k]]);
    for(j=0;j<nvalues[k];j++){
      sum[k][j]    = 0.;
      sum2[k][j]   = 0.;
      minval[k][j] = +INFINITY;
      maxval[k][j] = -INFINITY;
      }
    }
  
/*---------------------------------------------------------------------*//**<h1>
  compute the sums </h1>*/
  /// reading the input files one by one
  struct timeval b4;
  
  gettimeofday(&b4);
  
  for(i=0;i<filelist.size();i++){//for all files
    int nframes;//number of input frames
    date_t origine;
    double *time=0,od;
    
    gettimeofday(&before);
    if(ts!=0)
      status= poc_gettime(filelist[i].c_str(), &origine, &time, &nframes);
    else
      nframes=1;
    rrt+=difftime(before);
    //status=poc_decode_axis(info[0], global, &(decoded[0]));
    od=cnes_time(origine,'s');
    if(isnan(od))od=origd;
    printf("\r(%d/%d;%04.3g)reading %s ...",i,filelist.size(),difftime(b4),filelist[i].c_str());fflush(stdout);
    
    /// and frame by frame
    for(ifi=0;ifi<nframes;ifi++) {//for all frames
      int tfi;
      if(time!=0){
        tfi=pos(time[ifi]+od,ts,nt);
        if(tfi<0)
          continue;
        if(skipFrame!=0 and skipFrame[tfi])
          continue;
        }
      else
        tfi=ifi;
      ///\bug 2011-09-02 Damien Allain : assuming variables have the same sizes !!!
      printf("\r(%d/%d;%04.3g)reading %s (frame %d/%d,%g,%g) ",i,filelist.size(),difftime(b4),filelist[i].c_str(),tfi+1,nt,rt,ft);fflush(stdout);
      for(k=0;k<nvars;k++){
        gettimeofday(&before);
        status=poc_get_vara(filelist[i],info[k],ifi,inBuf[k]);
        rt+=difftime(before);
        if(status){NC_CHKERR_LINE_FILE(status,"poc_get_vara(\""+filelist[i]+"\",,,) error with "+info[0].name);exit(1);}
        }
      if(mask==NULL){
        int nmasked=0;
        /// if necessary allocating the mask to the number of points and setting its values from the values of the first variable
        mask=new bool[nvalues[0]];
        for(k=0;k<nvalues[0];k++) {
          mask[k]=(inBuf[0][k]!=spec[0]);
          nmasked+=!mask[k];
          }
        fprintf(stderr,"\n");
        STDOUT_BASE_LINE("%d / %d masked : %d\n",nmasked,nvalues[0],mask[nvalues[0]/2]);
        }
      /// and doing the sums
      gettimeofday(&before);
      for(k=0;k<nvars;k++){
        double *inBufk=inBuf[k],*sumk=sum[k],*sum2k=sum2[k],*minvalk=minval[k],*maxvalk=maxval[k];
#define debug_scalar_statistics 0
#if !debug_scalar_statistics
        #pragma omp parallel for private(j)
#endif
        for(j=0;j<nvalues[k];j++){
          //note : j may well be optimised out here
          const double value=inBufk[j];
          if(value==spec[k]){
            sumk[j]=value;
            continue;
            }
          sumk[j]+=value;
          sum2k[j]+=square(value);
          updatemin(&minvalk[j],value);
          updatemax(&maxvalk[j],value);
          }
        }
      ft+=difftime(before);
      }//EO for all frames
    deletep(&time);
    }//EO for all files
  printf("\n");

/*---------------------------------------------------------------------*//**<h1>
  compute masks, mean and stddev </h1>*/
  for(k=0;k<nvars;k++){
    double *sumk=sum[k],*sum2k=sum2[k],*minvalk=minval[k],*maxvalk=maxval[k];
#if !debug_scalar_statistics
    #pragma omp parallel for private(j)
#endif
    for(j=0;j<nvalues[k];j++){
      if(sumk[j]==spec[k]){
        sum2k[j]=spec[k];
        minvalk[j]=spec[k];
        maxvalk[j]=spec[k];
        continue;
        }
      sumk[j]/=nt;
      sum2k[j]=sqrt(sum2k[j]/nt-sumk[j]*sumk[j]);
      /**
      We have the mean \f$ \mu \f$ and standard deviation \f$ \sigma \f$ :
      \f{eqnarray*}{
      \mu &=& \frac{\sum x_i}{n} \\
      \sigma &=& \sqrt{\frac{\sum (x_i-\mu)^2}{n}} \\
        &=& \sqrt{\frac{\sum (x_i^2- 2 x_i \mu + \mu^2)}{n}}  \\
        &=& \sqrt{\frac{\sum x_i^2}{n}-2\frac{\sum x_i}{n}\mu+\mu^2}  \\
        &=& \sqrt{\frac{\sum x_i^2}{n}-\mu^2}
      \f}
      */
      }
    }
  
/*---------------------------------------------------------------------*//**<h1>
  read the coordinates </h1>*/
  poc_grid_data_t *coordinates=NULL;
  STDERR_BASE_LINE("#atlases:\n");
  
  coordinates=new poc_grid_data_t[nvars];
  
  for(k=0;k<nvars;k++){
    status=ENOENT;/* No such file or directory */
    if(gridfile!=NULL)
      status=poc_get_grid_data(gridfile,info[k],&coordinates[k]);
    if(status!=NC_NOERR)
      status=poc_get_grid_data(filelist[0],info[k],&coordinates[k]);
    if(status!=NC_NOERR){
      NC_CHKERR_LINE_FILE(status,"poc_get_grid_data() error (see failsafe info below) on %s with both %s and "+filelist[0],varnames[0],gridfile);
      STDERR_BASE_LINE("THE GRID VARIABLES, IF ANY, WILL NOT BE SAVED.\n");
      }
    }

/*---------------------------------------------------------------------*//**<h1>
  make atlases </h1>*/
  /// by calling atlas_()
  status=atlas_(nvars, sum,sum2,minval,maxval, mask, coordinates, info);
  STDERR_BASE_LINE("#atlases: done\n");
  
  //clean-up
  for(k=0;k<nvars;k++){
    delete[]inBuf[k];
    delete[]sum[k];
    delete[]sum2[k];
    }
  delete[]inBuf;
  delete[]sum;
  delete[]sum2;
  
  delete[] info;
  delete[] nvalues;
  delete[] n_spacedims;
  return NC_NOERR;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int scalar_statisticsUG(const char **varnames, int nvars, vector<string> filelist, char *gridfile, size_t nt,double *ts,const double origd=NAN)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// RMS of the filtered data
/**
\param **varnames list of variable names
\param nvars number of variables
\param filelist list of files
\param *gridfile path to grid file
\param *controlfile path to list of control points
\param nt number of frames
\param ts times of the frames in second
\param atlas whether an atlas is done
\param origd=NAN if specified, default time origin if not in the files
\returns the RMS of the filtered data, or NULL if error

\date 2011-09-19 Damien Allain : added frame times check
*/
/*----------------------------------------------------------------------------*/
{
  double **inBuf,**sum,**sum2,**minval,**maxval;
  int j,k,i;//variable and file indexes
  int status;
  int ifi,dfi;//input frame index,duplicate frame index
  poc_global_t global;
  poc_var_t *info;
  double *spec;
  bool *mask=NULL;
  size_t *nvalues,ntotal;
  int *n_spacedims;
  
  if(nvars<=0)return NC_EINVAL;

  nvalues=new size_t[nvars];
  n_spacedims=new int[nvars];
  spec=new double[nvars];
  
/*---------------------------------------------------------------------*//**<h1>
  if there are duplicates in the frame times given </h1>*/
  if(dupos(ts,nt)>=0){
    STDERR_BASE_LINE("Duplicate time frames :\n");
    for(ifi=0;ifi<nt-1;ifi++) {
      for(dfi=ifi+1;dfi<nt;dfi++) {
        if(ts[ifi]==ts[dfi]) {
          printf("[%d==%d]%s = %.10g\n",ifi,dfi,poctime_sdate_cnes(ts[dfi],'s'),ts[ifi]);
          break;
          }
        }
      }
    return NC_EINVAL;
    }
  
/*---------------------------------------------------------------------*//**<h1>
  gets variable dimensions </h1>*/
  info=new poc_var_t[nvars];
  gettimeofday(&before);
  status=poc_inq(filelist[0],&global,0);
  rrt+=difftime(before);
  ntotal=0;
  for(k=0;k<nvars;k++){
    status=global.variables.find(varnames[k],&info[k]);
    if(status<0){
      status=NC_ENOTVAR;NC_CHKERR_LINE_FILE(status,"poc_global_t::variables.find(\"%s\",) error in %s",varnames[k],filelist[0].c_str());return status;
      }
    status=poc_get_var_length(info[k],&nvalues[k],&n_spacedims[k]);
    if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,1,"poc_get_var_length error with %s in %s",varnames[k],filelist[0].c_str());
    STDERR_BASE_LINE("%s : %d values in %d dimensions\n",varnames[k],nvalues[k],n_spacedims[k]);
    ntotal+=nvalues[k];
    status=poc_decode_mask(info[k],NULL,NULL,&spec[k]);
    if(status==NC_ENOTATT) NC_CHKERR_LINE_FILE(status,"Use default mask value of %g for %s in %s",spec[k],varnames[k],filelist[0].c_str());//warning only
    }
  
  inBuf=new double*[nvars];
  
  sum   = new double*[nvars];
  sum2  = new double*[nvars];
  minval= new double*[nvars];
  maxval= new double*[nvars];
  
  for(k=0;k<nvars;k++){
    exitIfNull(inBuf[k]  = new double[nvalues[k]]);
    exitIfNull(sum[k]    = new double[nvalues[k]]);
    exitIfNull(sum2[k]   = new double[nvalues[k]]);
    exitIfNull(minval[k] = new double[nvalues[k]]);
    exitIfNull(maxval[k] = new double[nvalues[k]]);
    for(j=0;j<nvalues[k];j++){
      sum[k][j]    = 0.;
      sum2[k][j]   = 0.;
      minval[k][j] = +INFINITY;
      maxval[k][j] = -INFINITY;
      }
    }
  
/*---------------------------------------------------------------------*//**<h1>
  compute the sums </h1>*/
  /// reading the input files one by one
  struct timeval b4;
  
  gettimeofday(&b4);
  
  for(i=0;i<filelist.size();i++){//for all files
    int nframes;//number of input frames
    date_t origine;
    double *time,od;
    
    gettimeofday(&before);
    status= poc_gettime(filelist[i].c_str(), &origine, &time, &nframes);
    rrt+=difftime(before);
    od=cnes_time(origine,'s');
    if(isnan(od))od=origd;
    printf("\r(%d/%d;%04.3g)reading %s ...",i,filelist.size(),difftime(b4),filelist[i].c_str());fflush(stdout);
    
    /// and frame by frame
    for(ifi=0;ifi<nframes;ifi++) {//for all frames
      int tfi=pos(time[ifi]+od,ts,nt);
      if(tfi<0)
        continue;
      ///\bug 2011-09-02 Damien Allain : assuming variables have the same sizes !!!
      printf("\r(%d/%d;%04.3g)reading %s (frame %d/%d,%g,%g) ",i,filelist.size(),difftime(b4),filelist[i].c_str(),tfi+1,nt,rt,ft);fflush(stdout);
      for(k=0;k<nvars;k++){
        gettimeofday(&before);
        status=poc_get_vara(filelist[i],info[k],ifi,inBuf[k]);
        rt+=difftime(before);
        if(status){NC_CHKERR_LINE_FILE(status,"poc_get_vara(\""+filelist[i]+"\",,,) error with "+info[0].name);exit(1);}
        }
      if(mask==NULL){
        int nmasked=0;
        /// if necessary allocating the mask to the number of points and setting its values from the values of the first variable
        mask=new bool[nvalues[0]];
        for(k=0;k<nvalues[0];k++) {
          mask[k]=(inBuf[0][k]!=spec[0]);
          nmasked+=!mask[k];
          }
        fprintf(stderr,"\n");
        STDOUT_BASE_LINE("%d / %d masked : %d\n",nmasked,nvalues[0],mask[nvalues[0]/2]);
        }
      /// and doing the sums
      gettimeofday(&before);
      for(k=0;k<nvars;k++){
        double *inBufk=inBuf[k],*sumk=sum[k],*sum2k=sum2[k],*minvalk=minval[k],*maxvalk=maxval[k];
        #pragma omp parallel for private(j)
        for(j=0;j<nvalues[k];j++){
          //note : j may well be optimised out here
          const double value=inBufk[j];
          if(value==spec[k]){
            sumk[j]=value;
            continue;
            }
          sumk[j]+=value;
          sum2k[j]+=square(value);
          updatemin(&minvalk[j],value);
          updatemax(&maxvalk[j],value);
          }
        }
      ft+=difftime(before);
      }//EO for all frames
    delete[] time;
    }//EO for all files
  printf("\n");

/*---------------------------------------------------------------------*//**<h1>
  compute masks, mean and stddev </h1>*/
  for(k=0;k<nvars;k++){
    double *sumk=sum[k],*sum2k=sum2[k],*minvalk=minval[k],*maxvalk=maxval[k];
    #pragma omp parallel for private(j)
    for(j=0;j<nvalues[k];j++){
      if(sumk[j]==spec[k]){
        sum2k[j]=spec[k];
        minvalk[j]=spec[k];
        maxvalk[j]=spec[k];
        continue;
        }
      sumk[j]/=nt;
      sum2k[j]=sqrt(sum2k[j]/nt-sumk[j]*sumk[j]);
      /**
      We have the mean \f$ \mu \f$ and standard deviation \f$ \sigma \f$ :
      \f{eqnarray*}{
      \mu &=& \frac{\sum x_i}{n} \\
      \sigma &=& \sqrt{\frac{\sum (x_i-\mu)^2}{n}} \\
        &=& \sqrt{\frac{\sum (x_i^2- 2 x_i \mu + \mu^2)}{n}}  \\
        &=& \sqrt{\frac{\sum x_i^2}{n}-2\frac{\sum x_i}{n}\mu+\mu^2}  \\
        &=& \sqrt{\frac{\sum x_i^2}{n}-\mu^2}
      \f}
      */
      }
    }
  
// /*---------------------------------------------------------------------*//**<h1>
//   read the coordinates </h1>*/
//   poc_grid_data_t *coordinates=NULL;
//   STDERR_BASE_LINE("#atlases:\n");
//   
//   coordinates=new poc_grid_data_t[nvars];
//   
//   for(k=0;k<nvars;k++){
//     status=ENOENT;/* No such file or directory */
//     if(gridfile!=NULL)
//       status=poc_get_grid_data(gridfile,info[k],&coordinates[k]);
//     if(status!=NC_NOERR)
//       status=poc_get_grid_data(filelist[0],info[k],&coordinates[k]);
//     if(status!=NC_NOERR){
//       NC_CHKERR_LINE_FILE(status,"poc_get_grid_data() error (see failsafe info below) on %s with both %s and "+filelist[0],varnames[0],gridfile);
//       STDERR_BASE_LINE("THE GRID VARIABLES, IF ANY, WILL NOT BE SAVED.\n");
//       }
//     }

// /*---------------------------------------------------------------------*//**<h1>
//   make atlases </h1>*/
//   /// by calling atlas_()
//   status=atlas_(nvars, sum,sum2,minval,maxval, mask, coordinates, info);
//   STDERR_BASE_LINE("#atlases: done\n");
  
  //clean-up
  for(k=0;k<nvars;k++){
    delete[]inBuf[k];
    delete[]sum[k];
    delete[]sum2[k];
    }
  delete[]inBuf;
  delete[]sum;
  delete[]sum2;
  
  delete[] info;
  delete[] nvalues;
  delete[] n_spacedims;
  return NC_NOERR;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int vector_statistics(const char **varnames, int nvars, vector<string> filelist, char *gridfile, size_t nt,double *ts,const double origd=NAN)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// RMS of the filtered data
/**
\param **varnames list of variable names
\param nvars number of variables
\param filelist list of files
\param *gridfile path to grid file
\param *controlfile path to list of control points
\param nt number of frames
\param ts times of the frames in second
\param atlas whether an atlas is done
\param origd=NAN if specified, default time origin if not in the files
\returns the RMS of the filtered data, or NULL if error

\date 2011-09-19 Damien Allain : added frame times check
*/
/*----------------------------------------------------------------------------*/
{
  double **inBuf,**sum,**sum2,**minval,**maxval;
  int j,k,i;//variable and file indexes
  int status;
  int ifi,dfi;//input frame index,duplicate frame index
  poc_global_t global;
  poc_var_t *info;
  double *spec;
  bool *mask=NULL;
  size_t *nvalues,ntotal;
  int *n_spacedims;
  
  if(nvars<=0)return NC_EINVAL;
  
//   if(nvars!=1){
//     CHKERR_LINE_FILE(ENOEXEC,"no more than one variable accepted yet because of historical bug");
//     TRAP_ERR_EXIT(ENOEXEC,"exiting\n");/* Exec format error */
//     }
    
  nvalues     =new size_t[nvars];
  n_spacedims =new int[nvars];
  spec        =new double[nvars];
  
/*---------------------------------------------------------------------*//**<h1>
  if there are duplicates in the frame times given </h1>*/
  if(dupos(ts,nt)>=0){
    STDERR_BASE_LINE("Duplicate time frames :\n");
    for(ifi=0;ifi<nt-1;ifi++) {
      for(dfi=ifi+1;dfi<nt;dfi++) {
        if(ts[ifi]==ts[dfi]) {
          printf("[%d==%d]%s = %.10g\n",ifi,dfi,poctime_sdate_cnes(ts[dfi],'s'),ts[ifi]);
          break;
          }
        }
      }
    ///return NC_EINVAL
    return NC_EINVAL;
    }
  
/*---------------------------------------------------------------------*//**<h1>
  gets variable dimensions </h1>*/
  info=new poc_var_t[nvars];
  gettimeofday(&before);
  status=poc_inq(filelist[0],&global,0);
  rrt+=difftime(before);
  
  ntotal=0;
  for(k=0;k<nvars;k++){
    status=global.variables.find(varnames[k],&info[k]);
    if(status<0){status=NC_ENOTVAR;NC_CHKERR_LINE_FILE(status,"poc_global_t::variables.find(\"%s\",) error in %s",varnames[k],filelist[0].c_str());return status;}
    status=poc_get_var_length(info[k],&nvalues[k],&n_spacedims[k]);
    if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"poc_get_var_length error with %s in %s",varnames[k],filelist[0].c_str());
    STDERR_BASE_LINE("%s : %d values in %d dimensions\n",varnames[k],nvalues[k],n_spacedims[k]);
    ntotal+=nvalues[k];
    status=poc_decode_mask(info[k],NULL,NULL,&spec[k]);
    if(status==NC_ENOTATT)NC_CHKERR_LINE_FILE(status,"Use default mask value of %g for %s in %s",spec[k],varnames[k],filelist[0].c_str());//warning only
    }
  nvalues[nvars]=1;//see normalisation
  
  ///to allocate buffers
  inBuf  =new double*[nvars];
  sum    =new double*[nvars];
  sum2   =new double*[nvars];
  minval =new double*[nvars];
  maxval =new double*[nvars];
  for(k=0;k<=nvars;k++){
    exitIfNull( inBuf[k]  =new double[nvalues[k]] );
    exitIfNull( sum[k]    =new double[nvalues[k]] );
    exitIfNull( sum2[k]   =new double[nvalues[k]] );
    exitIfNull( minval[k] =new double[nvalues[k]] );
    exitIfNull( maxval[k] =new double[nvalues[k]] );
    for(j=0;j<nvalues[k];j++){
      sum[k][j] =0.;
      sum2[k][j]=0.;
      minval[k][j]=+INFINITY;
      maxval[k][j]=-INFINITY;
      }
    }
  
/*---------------------------------------------------------------------*//**<h1>
  compute the sums </h1>*/
  /// reading the input files one by one
  struct timeval b4;
  
  gettimeofday(&b4);
  
  for(i=0;i<filelist.size();i++){//for all files
    int nframes;//number of input frames
    date_t origine;
    double *time,od;
    
    gettimeofday(&before);
    status= poc_gettime(filelist[i].c_str(), &origine, &time, &nframes);
    rrt+=difftime(before);
    //status=poc_decode_axis(info[0], global, &(decoded[0]));
    od=cnes_time(origine,'s');
    if(isnan(od))od=origd;
    printf("\r(%d/%d;%04.3g)reading %s ...",i,filelist.size(),difftime(b4),filelist[i].c_str());fflush(stdout);
    
    /// and frame by frame
    for(ifi=0;ifi<nframes;ifi++) {//for all frames
      int tfi=pos(time[ifi]+od,ts,nt);
      if(tfi<0)
        continue;
      ///\bug 2011-09-02 Damien Allain : assuming variables have the same sizes !!!
      printf("\r(%d/%d;%04.3g)reading %s (frame %d/%d,%g,%g) ",i,filelist.size(),difftime(b4),filelist[i].c_str(),tfi+1,nt,rt,ft);fflush(stdout);
      for(k=0;k<nvars;k++){
        gettimeofday(&before);
        status=poc_get_vara(filelist[i],info[0],ifi,inBuf[k]);
        rt+=difftime(before);
        if(status){NC_CHKERR_LINE_FILE(status,"poc_get_vara(\""+filelist[i]+"\",,,) error with "+info[0].name);exit(1);}
        }
      inBuf[nvars][0]=(tfi==0)?1.:0.;//normalisation input
      if(mask==NULL){
        int nmasked=0;
        /// if necessary allocating the mask to the number of points and setting its values from the values of the first variable
        mask=new bool[nvalues[0]];
        for(k=0;k<nvalues[0];k++) {
          mask[k]=(inBuf[0][k]!=spec[0]);
          nmasked+=!mask[k];
          }
        fprintf(stderr,"\n");
        STDOUT_BASE_LINE("%d / %d masked : %d\n",nmasked,nvalues[0],mask[nvalues[0]/2]);
        }
      /// and doing the sums
      gettimeofday(&before);
      for(k=0;k<nvars;k+=2){
        double *inBufkx=inBuf[k];
        double *sumk=sum[k],  *sum2k=sum2[k],  *minvalk=minval[k],  *maxvalk=maxval[k];
        double *inBufky=inBuf[k+1];
        
//#pragma omp parallel for private(j)
        for(j=0;j<nvalues[k];j++){
          const double value_x=inBufkx[j];
          const double value_y=inBufky[j];
          if ( (value_x==spec[k]) || (value_y==spec[k+1]) ) {
            sumk[j]=value_x;
            continue;
            }
          const double value=sqrt(inBufkx[j]*inBufkx[j]+inBufky[j]*inBufky[j]);
          sumk[j]  +=value;
          sum2k[j] +=value*value;
          if(minvalk[j]>value){
            minvalk[j]=value;
            }
          if(maxvalk[j]<value){
            maxvalk[j]=value;
            }
          }
        }
      ft+=difftime(before);
      }
    delete[] time;
    }
  printf("\n");

/*---------------------------------------------------------------------*//**<h1>
  compute masks, mean and stddev </h1>*/
  for(k=0;k<nvars;k+=2){
    double *sumk=sum[k],*sum2k=sum2[k],*minvalk=minval[k],*maxvalk=maxval[k];
    #pragma omp parallel for private(j)
    for(j=0;j<nvalues[k];j++){
      if(sumk[j]==spec[k]){
        sum2k[j]   =spec[k];
        minvalk[j] =spec[k];
        maxvalk[j] =spec[k];
        continue;
        }
      sumk[j]/=nt;
      sum2k[j]=sqrt(sum2k[j]/nt-sumk[j]*sumk[j]);
      /**
      We have the mean \f$ \mu \f$ and standard deviation \f$ \sigma \f$ :
      \f{eqnarray*}{
      \mu &=& \frac{\sum x_i}{n} \\
      \sigma &=& \sqrt{\frac{\sum (x_i-\mu)^2}{n}} \\
        &=& \sqrt{\frac{\sum (x_i^2- 2 x_i \mu + \mu^2)}{n}}  \\
        &=& \sqrt{\frac{\sum x_i^2}{n}-2\frac{\sum x_i}{n}\mu+\mu^2}  \\
        &=& \sqrt{\frac{\sum x_i^2}{n}-\mu^2}
      \f}
      */
      }
    }
  
/*---------------------------------------------------------------------*//**<h1>
  read the coordinates </h1>*/
  poc_grid_data_t *coordinates=NULL;
  STDERR_BASE_LINE("#atlases:\n");
  
  coordinates=new poc_grid_data_t[nvars];
  
  for(k=0;k<nvars;k++){
    status=ENOENT;/* No such file or directory */
    if(gridfile!=NULL)
      status=poc_get_grid_data(gridfile,info[k],&coordinates[k]);
    if(status!=NC_NOERR)
      status=poc_get_grid_data(filelist[0],info[k],&coordinates[k]);
    if(status!=NC_NOERR){
      NC_CHKERR_LINE_FILE(status,"poc_get_grid_data() error (see failsafe info below) on %s with both %s and "+filelist[0],varnames[0],gridfile);
      STDERR_BASE_LINE("THE GRID VARIABLES, IF ANY, WILL NOT BE SAVED.\n");
      }
    }

  nvars--;
  status= atlas_(nvars, sum, sum2, minval, maxval, mask, coordinates, info);
  STDERR_BASE_LINE("#atlases: done\n");
  
  for(k=0;k<nvars;k++){
    delete[]inBuf[k];
    delete[]sum[k];
    delete[]sum2[k];
    }
  delete[]inBuf;
  delete[]sum;
  delete[]sum2;
  
  delete[] info;
  delete[] nvalues;
  delete[] n_spacedims;
  return NC_NOERR;
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
    "  %s file1 [ file2 ... ] [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Makes statistics on NetCDF outputs.\n"
    "  It takes a file or a list of files as input and produces the statistics atlases.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  Show this help and exit.\n"
    "  --take-first  if some time frames are simultaneous, take the first one\n"
    "  --take-last  if some time frames are simultaneous, take the last one\n"
    "  -c  followed by convention CURRENTLY WITHOUT EFFECT!\n"
    "  -l  followed by the path of the list of files to analyse. This list will override the list given as arguments.\n"
    "  -g  followed by the path of the grid\n"
    "  -s  followed by the start date. See DATE FORMATS below\n"
    "  -f  followed by the end date. See DATE FORMATS below\n"
    "  -o  followed by the default date origin. See DATE FORMATS below\n"
    "  -v  followed by the list of variables to make stats upon\n");
  print_poctime_scan_date_help(0);
  print_OPENMP_help(prog_name); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// function that will be ran when the executable is started
/** See the print_help <a href=#func-members>function</a> for its use.
*/
/*----------------------------------------------------------------------------*/
{
  int k;//variable index
  int n,status;//argument index or number of files, NetCDF status
  char *gridfile=NULL,*keyword;
  char *convention=NULL,*list=0;
  date_t start=NADate,final=NADate,orig=NADate;
  double startd,finald,origd,*ts;
  harmonic_t harmonic;
  int nvars=0;
  int handleRepeatedFrames=0;
  char **varnames=NULL;
  vector<string> filelist;
  
  struct timeval mainbefore;
  
  #if _OPENMP < 200505  // Version 2.5 May 2005
  #error Check whether your version _OPENMP (below) of openmp can do nesting properly. If so, update the above line accordingly.
  #pragma message "_OPENMP=" TOSTRING(_OPENMP)
  #endif

  cmd=fct_echo( argc, argv);
  argv0=argv[0];

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"--help")==0) {
      print_help(argv[0]);
      exit(0);
      }
    else if(strcmp(keyword,"--take-first")==0) {
      if(handleRepeatedFrames>0){
        STDOUT_BASE_LINE(" *** either --take-first or --take-last allowed ***\n");
        print_help(argv[0]);
        wexit(-1);
        }
      handleRepeatedFrames=-1;
      n++;
      continue;
      }
    else if(strcmp(keyword,"--take-last")==0) {
      if(handleRepeatedFrames<0){
        STDOUT_BASE_LINE(" *** either --take-first or --take-last allowed ***\n");
        print_help(argv[0]);
        wexit(-1);
        }
      handleRepeatedFrames=1;
      n++;
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 'c' : /// \bug 2011-09-02 Damien Allain : -c option currently without effect
          convention= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'l' :
          list= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'g' :
          gridfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          status=poctime_scan_date(argv[n+1],&start,0);
          n++;
          n++;
          break;

        case 'f' :
          status=poctime_scan_date(argv[n+1],0,&final);
          n++;
          n++;
          break;

        case 'o' :
          status=poctime_scan_date(argv[n+1],&orig,0);
          n++;
          n++;
          break;

        case 'v' :
          for(k=1;n+k<argc && argv[n+k][0]!='-';k++);
          nvars=k-1;
          if(varnames){
            __FILE_LINE__(stdout,"multiple use of -v : the last one overrides the previous.");
            for(k=0;varnames[k]!=NULL;k++){
              free(varnames[k]);
              }
            delete[] varnames;
            }
          varnames=new char*[nvars+1];
          for(k=0;k<nvars;k++){
            varnames[k]=strdup(argv[n+1+k]);
            }
          varnames[nvars]=NULL;
          n++;
          n+=nvars;
          break;

        case 'h' :
          print_help(argv[0]);
          exit(0);
          break;

        default:
          printf("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
        filelist.push_back((string) argv[n]);
        n++;
        break;
      }
      free(keyword);
    }

  printf("#number of variables : %d\n",nvars);
  if(nvars==0) {
    printf("\n%s:%d: *** Please give the names of the variables to filter ***\n",__FILE__,__LINE__);
    print_help(argv[0]);
    exit(-1);
    }
  for(k=0;k<nvars;k++) {
    printf("#variable %d: %s\n",k,varnames[k]);
    }

  if(list!=0) {
    filelist=load_filelist(list);
    }

  if(filelist.size()==0) {
    printf("\n%s:%d: *** Please give a NetCDF file to analyse ***\n",__FILE__,__LINE__);
    print_help(argv[0]);
    exit(-1);
    }

  if(gridfile!=NULL && access(gridfile,R_OK)) {//NOTE: access rights may change during the run
    //This is to avoid crashing without output after doing the spectral analysis
    check_error(errno,__LINE__,__FILE__,"No read access to %s",gridfile);
    exit(1);
    }
  
  gettimeofday(&mainbefore);
  
/*------------------------------------------------------------------------------
  It checks the list of files against the start and end dates with poc_timefilterlist() */
  startd=cnes_time(start,'s');
  finald=cnes_time(final,'s');
  origd=cnes_time(orig,'s');
  printf("Selected data from %s to %s.\n",poctime_sdate_cnes(startd,'s'),poctime_sdate_cnes(finald,'s'));
  
  n=filelist.size();
  nt=poc_timefilterlist(&filelist,&startd,&finald,&ts,origd);
  
  printf("\n");
  printf("Kept %d/%d files with %d frames from %s to %s.\n\n",filelist.size(),n,nt, poctime_sdate_cnes(startd,'s'),poctime_sdate_cnes(finald,'s'));
  n=filelist.size();
  if(nt==0)
    TRAP_ERR_EXIT(1,"no frames within time boundaries\n");
  //start=poctime_getdatecnes(startd,'s');//not usefull

/*------------------------------------------------------------------------------
  It finally calls variables_statistics() */
  status=scalar_statistics((const char**)varnames,nvars,filelist,gridfile,nt,ts,origd,handleRepeatedFrames);
//  status=vector_statistics((const char**)varnames,nvars,filelist,gridfile,nt,ts,origd);
  if(status!=NC_NOERR){
    STDOUT_BASE_LINE("statistics failed\n");
    exit(-1);
    }

  STDOUT_BASE_LINE(" ^^^^^^^^^^^^^ end of computation, that took %#g s with rt=%#g,rt2=%#g,rrt=%g,wt=%#g,filtering=%#g ^^^^^^^^^^^^^\n",difftime(mainbefore),rt,rt2,rrt,wt,ft);
  exit(0);
}
