
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Detides comodo-compliant NetCDF outputs and produces tidal atlases.

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
#include <unistd.h> //unlink
#include <sys/stat.h> //mkdir

#include <stdint.h> //uint64_t
#include <sys/time.h> //gettimeofday

#include "tools-structures.h"

#include "tides.h"
#include "functions.h" //safely includes omp.h
#include "matrix.h"
#include "poc-netcdf-data.hpp"
#include "poc-time.h"
#include "filter.h"
#include "map.h"

struct timeval stv;///<start timeval, for progression
///These global variables are used to count the amount of time this programme took for I/O tasks
struct timeval before;
double rt=0.,rt2=0.,rrt=0.,wt=0.;// read, 2nd read, redundant read and write times
double hst=0.,hat=0.,hct=0.;// harmonic storage, analysis and correction times.

string cmd;/*< command line */

const char *argv0=NULL;
void print_help(const char *prog_name);

poc_att_t period_att;

#define WAVE_AS_SEPARATE_FILE 0
#define WAVE_AS_DIMENSION 1


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  string atlas_file_name(const char *outputDir,const int format,const char * waveName,const string & varName)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// write atlas
/**
\param *constants array[] of hconstant_t
\param *mask
\param spectrum
\param coordinates will be descaled at write time if necessary.
\param info
\param *outputDir
\param format
*/
/*----------------------------------------------------------------------------*/
{
  string s;
  
  if(outputDir!=NULL)
    s=string(outputDir)+"/";
  
  if(format==WAVE_AS_SEPARATE_FILE){
    s=s+waveName+"-";
    }
  s+=varName+"-atlas.nc";
  
  return s;
}


#if 0
#define T double
#define NC_T NC_DOUBLE
#define NC_FILL_T NC_FILL_DOUBLE
#else
#define T float
#define NC_T NC_FLOAT
#define NC_FILL_T NC_FILL_FLOAT
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int atlas_(const hconstant_t *constants, const bool *mask, const spectrum_t &spectrum,poc_grid_data_t &coordinates,const poc_var_t &info, const char *outputDir,const int format)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// write atlas
/**
\param *constants array[] of hconstant_t
\param *mask
\param spectrum
\param coordinates will be descaled at write time if necessary.
\param info
\param *outputDir
\param format
*/
/*----------------------------------------------------------------------------*/
{
  int i,w,status;
  size_t length;
  poc_global_t global("constructed around " __LINE_FILE_PACKAGE_REVISION);
  poc_var_t aV,gV,wV,wNV("waveNames",NC_CHAR);//<variables for amplitude, phase, and waves frequencies and names
  poc_dim_t wDim("W",spectrum.n);
  poc_att_t spectrumA("note","Darwin spectrum. See " WAVELISTREF);
  poc_att_t *att;
  string s="";
  double scale,offset;
  const T spec=-NC_FILL_T;
  T *aData,*gData;

/*---------------------------------------------------------------------*//**<h1>
  initialise the file </h1>*/
  global << period_att;
  global << poc_att_t("history",cmd);
  
  /// with the coordinates variables
  string coordinates_vars;
  for(i=0;i<coordinates.vc;i++){
    poc_var_t *var=&coordinates[i].info;
    if(var->name=="")
      continue;
    global.variables << *var;
    if(i<coordinates.vc/2)
      /* re-do coordinates atttribute to make sure:
        1) it is there
        2) it contains no time variable */
      coordinates_vars.insert(0,var->name+
        (coordinates_vars==""?"":" ") );
    }
  status=poc_get_var_length(info,&length);
  /// with the amplitude and phase variables.
  s=string(info.name)+AMPSUF;
  aV.init(s.c_str(),NC_T);
  s=string(info.name)+PHASUF;
  gV.init(s.c_str(),NC_T);
  
  switch(format){
  case WAVE_AS_SEPARATE_FILE:
    break;
  
  case WAVE_AS_DIMENSION:
    aV.dimensions << wDim;
    gV.dimensions << wDim;
    
    wV.init(wDim,NC_DOUBLE);
    wV << poc_att_t("units","deg/h");
    wV.attributes << spectrumA;
    
    wNV.dimensions << wDim << poc_dim_t("maxWaveNameLen",MAXIMUM_WAVE_NAME_LENGTH);
    wNV << poc_att_t("standard_name","tidal_frequency_name");
    wNV.attributes << spectrumA;
    
    global.variables << wV << wNV;
    break;
  
    //NOTE: waves could also be on different variable names
  
  default:
    TRAP_ERR_EXIT(ENOEXEC,"programming error");
    }
  
  for(i=0;i<info.dimensions.size();i++){
    if(isT(info.dimensions[i]))continue;
    aV << info.dimensions[i];
    gV << info.dimensions[i];
    }
  
  status=poc_decode_mask(info,&scale,&offset,(double*)0);
  offset/=scale;
  aV << poc_att_t("_FillValue",spec);
  gV << poc_att_t("_FillValue",spec);
  aV << poc_att_t("coordinates",coordinates_vars);
  gV << poc_att_t("coordinates",coordinates_vars);
  for(i=0;i<info.attributes.size();i++){
    att=(poc_att_t*)&info.attributes[i];
    if(att->name=="long_name"){
      aV << poc_att_t(att->name,"tidal amplitude of "+att->as_string());
      gV << poc_att_t(att->name,"tidal phase of "+att->as_string());
      continue;
      }
    if(att->name=="standard_name"){
      aV << poc_att_t(att->name,att->as_string()+"_amplitude_due_to_non_equilibrium_ocean_tide");
      gV << poc_att_t(att->name,att->as_string()+"_phase_due_to_non_equilibrium_ocean_tide");
      continue;
      }
    if(att->name=="units" ||
       att->name=="scale_factor"){
      aV << *att;
      continue;
      }
    if(att->name=="associate" ||
       att->name=="coordinates" || /* already specified by default */
       att->name=="content" ||
       att->name=="valid_min" ||
       att->name=="valid_max" ||
       att->name=="valid_range" ||
       att->name=="add_offset" || //offset will be numerically removed for Z0
       att->name=="_FillValue" || //already specified by default
       att->name=="missing_value"){
      continue;
      }
    aV << *att;
    gV << *att;
    }
  gV << poc_att_t("units","degrees");
  global.variables << aV << gV;
  //poc_print(global);

  aData=new T[length];
  gData=new T[length];
  for(w=0;w<spectrum.n;w++){/// <h1>For each wave</h1>
    
    const tidal_wave *wave=&spectrum.waves[w];
    int frame=0;
    
    ///create the file
    if(format==WAVE_AS_SEPARATE_FILE || w==0){
      s=atlas_file_name(outputDir,format,wave->name,info.name);
      STDERR_BASE_LINE("(over)writing %s ...",s.c_str());
      unlink(s.c_str());
      status=poc_create(s,global,1,NC_NETCDF4);
      if(status!=NC_NOERR){NC_CHKERR_BASE_LINE(status,"poc_create(\""+s+"\",) error");poc_print(global);printf("\n");continue;}
      ///write the grid variables
      for(i=0;i<coordinates.vc;i++){
        const string *varName=&coordinates[i].info.name;
        if(*varName==""){
          continue;}
        fprintf(stderr," "+*varName+",");
        status=coordinates[i].write_data(s.c_str(),0,0,1);//descaling
        NC_CHKERR_BASE_LINE(status,"poc_data_t(\""+*varName+"\")::write_data(\""+s+"\") error");
        }
      
      }
    
    if(format==WAVE_AS_DIMENSION){
      frame=w;
      char waveName[MAXIMUM_WAVE_NAME_LENGTH+1];
      
      fprintf(stderr,"\n%s: "+wV.name+",",wave->name);
      poc_put_vara(s.c_str(),wV,frame,&wave->omega);
      
      memset(waveName,0,MAXIMUM_WAVE_NAME_LENGTH+1);
      strncpy(waveName,wave->name,MAXIMUM_WAVE_NAME_LENGTH+1);
      fprintf(stderr," "+wNV.name+",");
      poc_put_vara(s.c_str(),wNV,frame,waveName);
      }
    
    ///and write the data with poc_put_vara(string,poc_var_t,size_t,T*,int,int)
    for(i=0;i<length;i++){
      if(mask[i]){
        aData[i]=constants[i].a[w];
        ///removing the offset for wave Z0
        if(wave->omega==0.){
          aData[i]+=offset;
          }
        gData[i]=constants[i].G[w];
        }
      else{
        aData[i]=spec;
        gData[i]=spec;
        }
      }
    
    fprintf(stderr," "+aV.name+",");
    status=poc_put_vara(s,aV,frame,aData,0,0);//no descaling
    NC_CHKERR_BASE_LINE(status,"poc_put_vara(\""+s+"\",\""+aV.name+"\",,) error");
    
    fprintf(stderr," "+gV.name+",");
    status=poc_put_vara(s,gV,frame,gData,0,0);//no descaling
    NC_CHKERR_BASE_LINE(status,"poc_put_vara(\""+s+"\",\""+gV.name+"\",,) error");
    
    if(format==WAVE_AS_SEPARATE_FILE){
      fprintf(stderr," done.\n");
      }
    }
  
  if(format==WAVE_AS_DIMENSION){
    fprintf(stderr,"\ndone.\n");
    }
  
  delete[]aData;
  delete[]gData;
  #undef T
  #undef NC_T

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void atlas_(int nvars,hconstant_t **constants, bool **mask, const spectrum_t &spectrum,poc_grid_data_t *coordinates,const poc_var_t *info, const char *outputDir,int format)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// write atlases
/**
\param nvars
\param[in] **constants array[nvars][] of hconstant_t
\param *mask
\param spectrum
\param[in,out] *coordinates array[nvars] of coordinates. Will be descaled at write time if necessary.
\param *info array[nvars]
\param *outputDir
\param format
*/
/*----------------------------------------------------------------------------*/
{
  int k;
  
  for(k=0;k<nvars;k++)
    atlas_(constants[k],mask[k],spectrum,coordinates[k],info[k],outputDir,format);
}


#define DO_STATISTICS 0


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int detide(vector<string> filelist, poc_var_t *info, size_t nt,double *ts, const harmonic_t &harmonic, hconstant_t **constants, int *controlIs,double origd=NAN, char *outputDir=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Detiding.
/**
\param filelist list of files
\param *info list of time and variable informations
\param nt number of frames
\param *ts times of the frames
\param harmonic contains list of waves, number of values and variables, mask, etc...
\param **constants harmonic constants
\param *controlIs control indexes
\param origd=NAN if specified, default time origin if not in the files
\param *outputDir=NULL if specified, path to the directory where the output files will be written. IT IS STRONGLY ADVISED TO OUTPUT ON A DIFFERENT HARD DRIVE WHEN DOING DETIDING.
\returns NC_NOERR if success or the NetCDF error code if error.

\date 2011-09-19 Damien Allain : added frame times check, started documentation
\date 2011-09-27 Damien Allain : parallelisation
*/
/*----------------------------------------------------------------------------*/
{
  const int nvars=harmonic.nrhs;
  const size_t *nvalues=harmonic.nndes;
  
  float ***buffers,**buffer;
  int nproc;//maximum number of threads
  int k,n,i;//indexes of variables or samples, values or samples, input files or process
  int ifi;//input frame index
  int status;
  cdfatt_t aB;//attribute buffer
  int nsamples;
  double **residuals=NULL,**signal=NULL,*stime;//2 arrays[nsamples][nt] for the control samples and control time axis
#if DO_STATISTICS
  double **mean=NULL,**rms=NULL;
#endif
  FILE *out;//control output files
  char filename[1024];//control output file name
  int tfi;//total frame index
  double dt,startd,finald;//sampling(h), start and end (s since 1950/01/01)
  int parallelfiles;
  if(outputDir==NULL)outputDir=strdup(".");
  
  ///If there are duplicates in the frame times given, return NC_EEXIST
  if(dupos(ts,nt)>=0)return NC_EEXIST;
  
  minpos(ts,nt,&startd);
  maxpos(ts,nt,&finald);
  dt=(finald-startd)/(nt-1)/3600.;
  
  //check sample indexes
  nsamples=0;
  if(controlIs!=NULL){
    for(;-1<controlIs[nsamples] && controlIs[nsamples]<nvalues[0];nsamples++);//initialises nsamples
    }
  if(nsamples>0){
    //allocate samples
    residuals=new double*[nsamples];
    signal   =new double*[nsamples];
    for(n=0;n<nsamples;n++) {
      residuals[n]=new double[nt];
      signal[n]   =new double[nt];
      }
    }
  
  //allocate memory for vectors
  #ifdef OMP_H
  nproc=omp_get_max_threads();//1 set of buffers per process
  #else
  nproc=1;
  #endif
  buffers=new float**[nproc];
  for(i=0;i<nproc;i++){
    buffers[i]=new float*[nvars];
    for(k=0;k<nvars;k++){
      exitIfNull(buffers[i][k]=new float[nvalues[k]]);
      }
    }
  #ifndef OMP_H
  buffer=buffers[0];
  #endif
  
  stime    =new double[nt];
  
#if DO_STATISTICS
  mean  =new double*[nvars];
  rms   =new double*[nvars];
  for(k=0;k<nvars;k++) {
    mean[k] =new double[nvalues];
    rms [k] =new double[nvalues];
    for(n=0;n<nvalues;n++) {
      mean[k][n] =0.;
      rms [k][n] =0.;
      }
    }
#endif
  
/*------------------------------------------------------------------------------
  for each input files */
  gettimeofday(&stv);
  // \note 2011-09-22 Damien Allain : had a look at http://bisqwit.iki.fi/story/howto/openmp/#LoopNesting
  parallelfiles=(nt/filelist.size()<nproc*9);
  #pragma omp parallel for private(i,status,ifi,k, tfi,buffer,n) schedule(dynamic,1) if(nproc>1 && parallelfiles)
  for(i=0; i<filelist.size(); i++){
    poc_global_t global;
    int nif,nof;//number of input and output frames
    date_t origine;
    double *itime,od,*otime;
    poc_var_t detided[nvars],detided_vtime,vtime;
    double scale,offset;
    float fmask[nvars];
    string output;
    int ofi;//output frame index
    #ifdef OMP_H
    int fti=omp_get_thread_num();//file thread index
    #endif
    #if NETCDF_CAN_PARALLEL_IO == 0
    #pragma omp critical(threadUnsafeNetcdf)
    #endif
    {
    //status= cdf_globalinfo (filelist[i].c_str(),&global,0);
    status = poc_gettime(filelist[i].c_str(), &origine, &itime, &nif, &vtime);
    NC_CHKERR_BASE_LINE(status,"poc_gettime(\""+filelist[i]+"\") error");
    }//END OF #pragma omp critical(threadUnsafeNetcdf)
    if(status!=NC_NOERR)continue;
    od=cnes_time(origine,'s');
    if(isnan(od))od=origd;
    otime=new double[nif];//don't bother about over-allocating for a few 100s of values when we are allocating for 10000s
    for(ofi=0,ifi=0;ifi<nif;ifi++){
      tfi=pos(itime[ifi]+od,ts,nt);
      if(tfi<0)continue;
      otime[ofi]=itime[ifi];
      ofi++;
      }
    nof=ofi;
    
    if(vtime.name==""){/* happens with grib files */
      vtime.init("time",NC_DOUBLE);
      
      }
    detided_vtime=vtime;
    
    /* for extracts that are only partly covering a file */
    if(detided_vtime.dimensions.size()>0)
      detided_vtime.dimensions[0].len=nof;
    else
      detided_vtime << poc_dim_t("time",nof, true);
    
    char *sdate=sgetdate(origine);
    detided_vtime << poc_att_t("units","seconds since "+(string)sdate);
    delete[]sdate;
    
    ///It OVERWRITES OR creates a corresponding output file
    output=(string)outputDir+"/detided-"+strrchr0(filelist[i].c_str(),'/');
    unlink(output.c_str());
    ///with the detided variables
    for(k=0;k<nvars;k++) {
      detided[k]=info[k];
      detided[k].name=detided[k].name+"_detided";
      //setting the time dimension to nof
      detided[k] << detided_vtime.dimensions[0];
      poc_decode_mask(info[k],&scale,&offset,&fmask[k],-NC_FILL_FLOAT);
      global.variables << detided[k];
      }
    /// and the time variable.
    global.variables << detided_vtime;
    #if NETCDF_CAN_PARALLEL_IO == 0
    #pragma omp critical(threadUnsafeNetcdf)
    #endif
    {
    printf("\r(%d/%d;%04.3g)writing "+output+" ...%s",i+1,filelist.size(),difftime(stv),el);fflush(stdout);
    gettimeofday(&before);
    status=poc_create(output,global,1,NC_NETCDF4);
    NC_CHKERR_BASE_LINE(status,"poc_create(\""+output+"\") error");
    wt+=difftime(before);
    }//END OF #pragma omp critical(threadUnsafeNetcdf)
    if(status!=NC_NOERR)continue;
    /// <h3>For each frame within this file</h3>
    #pragma omp parallel for private(ifi,ofi,tfi,status,buffer,k,n) schedule(dynamic,1) if(nproc>1 && !parallelfiles)
    for(ifi=0;ifi<nif;ifi++) {
      int bi;//buffer index
      tfi=pos(itime[ifi]+od,ts,nt);//total frame index
      if(tfi<0)continue;
      ofi=pos(itime[ifi],otime,nof);//output frame index
      if(ofi<0)continue;
      #ifdef OMP_H
      if(parallelfiles)
        bi=fti;
      else
        bi=omp_get_thread_num();
      buffer=buffers[bi];
      #endif
      #if NETCDF_CAN_PARALLEL_IO == 0
      #pragma omp critical(threadUnsafeNetcdf)
      #endif
      {
      printf("\r(%d/%d;%04.3g)writing "+output+" (frame %d/%d)...%s",i+1,filelist.size(),difftime(stv),tfi+1,nt,el);fflush(stdout);
      ///read the frame (for a second time) in the input file,
      gettimeofday(&before);
      for(k=0;k<nvars;k++) {
        status=poc_get_vara(filelist[i],info[k],ifi,buffer[k]);
        if(status)NC_TRAP_ERROR(wexit,status,1,"error while reading "+info[k].name+" on "+filelist[i]);
        }
      rt2+=difftime(before);
      }//END OF #pragma omp critical(threadUnsafeNetcdf)
      for(n=0;n<nsamples;n++) {
        signal[n][tfi]=buffer[0][controlIs[n]];
        }
      stime[tfi]=itime[ifi]+od-startd;
      
      /* detide */
      harmonic_correction(harmonic,tfi,constants,buffer);
      
      for(n=0;n<nsamples;n++) {
        residuals[n][tfi]=buffer[0][controlIs[n]];
        }
      #if NETCDF_CAN_PARALLEL_IO == 0
      #pragma omp critical(threadUnsafeNetcdf)
      #endif
      {
      if(harmonic.mask==0) TRAP_ERR_EXIT(ENOEXEC,"programming error : harmonic.mask=%p\n",harmonic.mask);
      for(k=0;k<nvars;k++) {
        if(harmonic.mask[k]==0) TRAP_ERR_EXIT(ENOEXEC,"programming error : harmonic.mask[%d]=%p\n",k,harmonic.mask[k]);
        for(n=0;n<nvalues[k];n++) {
          if(not harmonic.mask[k][n]) {
            buffer[k][n]=fmask[k];
            }
#if DO_STATISTICS
          else {
            ///if enabled (it is not) it will also do some statistics
            /// \bug 2011-09-22 Damien Allain : there should be one set of buffers per thread (see \c buffer), that should be merged later
            /// \bug 2011-09-22 Damien Allain : the statistics are not given in any output anyway
            mean[k][n] +=buffer[k][n];
            rms [k][n] +=buffer[k][n]*buffer[k][n];
            }
#endif
          }
        /// and writes it to the corresponding output file.
        printf("\r(%d/%d;%04.3g)writing "+output+" (frame %d/%d) th.%d "+info[k].name+"...%s",i+1,filelist.size(),difftime(stv),tfi+1,nt,bi,el);fflush(stdout);
        gettimeofday(&before);
        status=poc_put_vara(output,detided_vtime,ofi,&otime[ofi]);
        NC_CHKERR_BASE_LINE(status,"poc_put_vara(\""+output+"\",\""+detided_vtime.name+"\",,) error");
        status=poc_put_vara(output,detided[k],ofi,buffer[k]);
        if(status!=NC_ERANGE)NC_CHKERR_BASE_LINE(status,"poc_put_vara(\""+output+"\",\""+detided[k].name+"\",,) error");
        wt+=difftime(before);
        //STDERR_BASE_LINE("status=%d;\n",status);
        }
      }//END OF #pragma omp critical(threadUnsafeNetcdf)
      }
    delete[]itime;
    }
  printf("\n");

#if DO_STATISTICS
  for(k=0;k<nvars;k++) {
    double nvalid=0;
    for(n=0;n<nvalues;n++) {
      if(harmonic.mask[k][n]) {
        nvalid+=1;;
        }
      }
    for(n=0;n<nvalues;n++) {
      if(harmonic.mask[k][n]) {
        mean[k][n] /=nvalid;
        rms [k][n] =sqrt(rms[k][n]/nvalid-mean[k][n]*mean[k][n]);
        }
      }
    }
#endif
//   output=string("stat-");
//   output+=*p;
//   status= poc_createfile(output.c_str(),global);

  for(i=0;i<nproc;i++){
    for(k=0;k<nvars;k++){
      delete[]buffers[i][k];
      }
    delete[]buffers[i];
    }
  delete[]buffers;

  /// <h3>It also saves the values at control points :</h3>
  STDOUT_BASE_LINE("%d;\n",residuals!=NULL && signal!=NULL);
  if(residuals!=NULL && signal!=NULL){
    STDOUT_BASE_LINE("saving %d sample files\n",nsamples);
    for(n=0;n<nsamples;n++) {
      STDOUT_BASE_LINE("saving *-%s-%d.txt\n",info[0].name.c_str(),controlIs[n]);
      hconstant_t &c=constants[0][controlIs[n]];
      c.set_polar();
      /** - harmonic constants to \code /**/ // COMPILED CODE BELOW !!!
      sprintf(filename,"constants-%s-%d.txt",info[0].name.c_str(),controlIs[n]); /** \endcode */
      out=fopen(filename,"w");
      fprintf(out,"#wave amplitude phase(deg) real imag\n");
      for(i=0;i<harmonic.spectrum.n;i++){
        fprintf(out,"%10s %lf %lf %lf %lf\n", harmonic.spectrum.waves[i].name, c.a[i],c.G[i], real(c.z[i]),imag(c.z[i]));
        }
      fclose(out);
      ///\bug 2011-10-14 Damien Allain : saving constants should be in analyse() ?
      /** - time series to \code /**/ // COMPILED CODE BELOW !!!
      sprintf(filename,"series-%s-%d.txt",info[0].name.c_str(),controlIs[n]); /** \endcode */
      out=fopen(filename,"w");
      fprintf(out,"#time(days) signal residual\n");
      for(tfi=0;tfi<nt;tfi++) fprintf(out,"%lf %lf %lf\n", stime[tfi]/24./3600., signal[n][tfi], residuals[n][tfi]);
      fclose(out);
      /** - signal spectrum to \code /**/ // COMPILED CODE BELOW !!!
      sprintf(filename,"signal-fft-%s-%d.txt",info[0].name.c_str(),controlIs[n]); /** \endcode */
      fourier(filename, signal[n], -9999.9, nt, dt, 4);
      /** - residual spectrum to \code /**/ // COMPILED CODE BELOW !!!
      sprintf(filename,"residuals-fft-%s-%d.txt",info[0].name.c_str(),controlIs[n]); /** \endcode */
      fourier(filename, residuals[n], -9999.9, nt, dt, 4);
      /** You can easily plot these files in gnuplot using the code below : \code
load 'control.gplt' #generated by detide(). Edit this file if you want to plot but the last point.
set term wxt 0 size 900,500 title "series"
unset logscale y
plot 'series-'.v.'-'.p.'.txt' u 1:2 w l t 'signal' , '' u 1:3 w l t 'residuals'
set term wxt 1 size 900,500 title "FFTs"
set logscale y
plot 'signal-fft-'.v.'-'.p.'.txt' u ($1):(sqrt($2)) w l t 'signal' , 'residuals-fft-'.v.'-'.p.'.txt' u ($1):(sqrt($2)) w l t 'residuals'
      \endcode
      */
      delete[] residuals[n];
      delete[] signal[n];
      }
    delete[] residuals;
    delete[] signal;
    
    STDOUT_BASE_LINE("#########################################################\n"
      "  You can easily plot the spectra by running these commands in `ipython --pylab=...':\n"
      "%run control.gplt\n"
      "f,s=loadtxt('signal-fft-'+v+'-'+p+'.txt',usecols=(0,1),unpack=1)\n"
      "r=loadtxt('residuals-fft-'+v+'-'+p+'.txt',usecols=[1])\n"
      "plot(f,c_[s,r]**.5)\n"
      "xlabel('frequency (deg/sample)');legend(('signal','residuals'))\n"
      );
    STDOUT_BASE_LINE("#########################################################\n");
    }

#if DO_STATISTICS
  delete[] stime;
  for(k=0;k<nvars;k++) {
    delete[]mean[k];
    delete[]rms [k];
    }
  delete[]mean;
  delete[]rms;
#endif
#undef DO_STATISTICS

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  hconstant_t **analyse(const char **varnames,int nvars,vector<string> filelist,char *gridfile,char *controlfile,
                        size_t nt,double *ts,spectrum_t spectrum,int nodal_corrections, double averaging_time,
                        int variable_mask,const astro_angles_t &astro_angles,int atlas,double origd=NAN,
                        char *outputDir=NULL,const int format=WAVE_AS_SEPARATE_FILE,int handleRepeatedFrames=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Harmonic analysis, eventually atlas, and detiding.
/**
\param **varnames list of variable names
\param nvars number of variables
\param filelist list of files
\param *gridfile path to grid file
\param *controlfile path to list of control points
\param nt number of frames
\param ts times of the frames
\param spectrum list of waves
\param nodal_corrections whether nodal corrections are carried out
\param variable_mask whether mask values are variable
\param atlas whether an atlas is done
\param origd if specified, default time origin if not in the files
\param *outputDir if specified, path to the directory where the output files will be written. IT IS STRONGLY ADVISED TO OUTPUT ON A DIFFERENT HARD DRIVE WHEN DOING DETIDING.
\param handleRepeatedFrames see getSkippedFrames()
\returns the harmonic analysis, or NULL if error

\date 2011-09-19 Damien Allain : added frame times check
*/
/*----------------------------------------------------------------------------*/
{
  float ***buffers,**buffer;
  double ****rhS;//right-hand sides
  int nproc;//maximum number of threads
  double duration,startd,finald;
  int m,k,i;/*< node, variable and file indexes */
  int status;
  int ifi,sfi;//input frame index,simultaneous frame index
  poc_global_t global;
  poc_var_t *info;
  float *spec;
  bool **mask=NULL;
  size_t *nvalues,ntotal;
  int *n_spacedims;
  harmonic_t harmonic;
  hconstant_t **constants;
  int parallelfiles;
  
  if(nvars<=0)return NULL;
  nvalues=new size_t[nvars];
  n_spacedims=new int[nvars];
  spec=new float[nvars];
  buffer=new float*[nvars];
  
  /* handle repeated frames */
  
  bool *skipFrame=NULL;
  status=getSkippedFrames(nt,ts,handleRepeatedFrames,&skipFrame);
  if(status!=0){
    printf("\n");
    STDOUT_BASE_LINE(" *** simultaneous time frames: give either --take-first or --take-last ***\n");
    print_help(argv0);
    printf("\n");
    return NULL;
    }
  

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  calculate the duration of the analysis from the times of the frames
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
  minpos(ts,nt,&startd);
  maxpos(ts,nt,&finald);
  duration=finald-startd;
  period_att.init_copy("analysis_period",(string)"from "+poctime_sdate_cnes(startd,'s')+" to "+poctime_sdate_cnes(finald,'s'));
  
  /// <h3>get variable dimensions</h3>
  info=new poc_var_t[nvars];
  gettimeofday(&before);
  status=poc_inq(filelist[0],&global,0);
  rrt+=difftime(before);
  ntotal=0;
  for(k=0;k<nvars;k++){
    status=poc_get_var_length_and_mask(global,varnames[k],&info[k],&nvalues[k],&n_spacedims[k],&spec[k]);
    if(status==NC_ENOTATT)NC_CHKERR_BASE_LINE(status,"Use default mask value of %g for %s in "+filelist[0]+"\n",spec[k],varnames[k]);/* warning only */
    else if(status!=NC_NOERR){NC_CHKERR_BASE_LINE(status,"poc_get_var_length_and_mask(,\"%s\",...) error in "+filelist[0]+"\n",varnames[k]);return NULL;}
    
    ntotal+=nvalues[k];
    }
  
#ifdef OMP_H
  nproc=omp_get_max_threads();//1 set of buffers per process
#else
  nproc=1;
#endif
  
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  harmonic analysis and detiding
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
  if(spectrum.n<20)
    printAndGetSpectrumDetails(spectrum);
  omp_set_nested(1);
  STDERR_BASE_LINE("OMP:dynamic=%d,nested=%d,version=%d;\n",omp_get_dynamic(),omp_get_nested(),_OPENMP);
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  man:/sysconf says :
    long sysconf(int name)
  with name as :
  _SC_PAGESIZE : Size of a page in bytes.
  _SC_PHYS_PAGES : The number of pages of physical memory.
  _SC_AVPHYS_PAGES : The number of currently available pages of physical memory.
  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
#ifdef _SC_AVPHYS_PAGES
  STDERR_BASE_LINE("remains %ld of %ld *%ld %g%% of RAM\n",sysconf(_SC_AVPHYS_PAGES),sysconf(_SC_PHYS_PAGES),sysconf(_SC_PAGE_SIZE),100.*sysconf(_SC_AVPHYS_PAGES)/sysconf(_SC_PHYS_PAGES));
#endif
  #include "tides.def"
  spectrum.add(wZ0);
  
  /* allocate */
  harmonic=harmonic_start(spectrum, nvalues, nvars, nt);
  
#ifdef OMP_H
  do{///being carefull to limit the number of threads if necessary to avoid overfilling the RAM.
    uint64_t mem_proc=(harmonic.neq*sizeof(double)+sizeof(void*)+sizeof(float))*ntotal,
    ram_size=sysconf(_SC_PHYS_PAGES)*sysconf(_SC_PAGE_SIZE);
    int nprocmax=floor(.75*ram_size/mem_proc);
    STDERR_BASE_LINE("Needing %lub ~ (waves*2)*sizeof(double)*ntotal=%d*%d*%db of memory per process.\n",mem_proc,harmonic.neq,sizeof(double),ntotal);
    STDERR_BASE_LINE("There's %lub of RAM.\n",ram_size);
    if(nproc<nprocmax)break;
    STDERR_BASE_LINE("May only afford %d processes.\n",nprocmax);
    if(nprocmax>2){
      /** \note 2012-02 Damien Allain : I tested 2 shallow processes :
      - for many values, one for reading and one for writing : it is faster
      - for smaller amount of values : it is slightly slower */
      nprocmax=2;
      }
    STDERR_BASE_LINE("Using only %d of %d processors at bottom parallelisation level to avoid overfilling the RAM.\n",nprocmax,nproc);
    nproc=nprocmax;
    if(nproc>1){
      if(nproc!=2){//because omp_set_dynamic() disables one level of parallelisation
        omp_set_dynamic(1);
        }
      omp_set_num_threads(nproc);
      }
    } while(0);
#endif
  
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  allocate matrix and arrays
  
*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  buffers=new float**[nproc];
  rhS=new double***[nproc];//right-hand sides
  for(i=0;i<nproc;i++){
    buffers[i]=new float*[nvars];
    for(k=0;k<nvars;k++){
      exitIfNull(buffers[i][k]=new float[nvalues[k]]);
      }
    if(i>0)
      /// allocate and zeroe one rhs per extra CPU with harmonic_start(harmonic_t*,int,int)
      harmonic_start(&harmonic, nvalues, nvars);
    rhS[i]=harmonic.rhs;
#ifdef _SC_AVPHYS_PAGES
    STDERR_BASE_LINE("%d of %d : remains %ld of %ld *%ld %g%% of RAM\n",i+1,nproc,sysconf(_SC_AVPHYS_PAGES),sysconf(_SC_PHYS_PAGES),sysconf(_SC_PAGE_SIZE),100.*sysconf(_SC_AVPHYS_PAGES)/sysconf(_SC_PHYS_PAGES));
    if(0 && sysconf(_SC_AVPHYS_PAGES)<sysconf(_SC_PHYS_PAGES)/5) TRAP_ERR_EXIT(ENOMEM,"more than 80%% of RAM taken : aborting.\n");
#endif
    }
  
#ifndef OMP_H
  buffer=buffers[0];
  harmonic.rhs=rhS[0];
#endif
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  initialise the harmonic matrix
  
*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  STDERR_BASE_LINE("Initialising harmonic matrix for %d waves...",spectrum.n);
  harmonic_init(&harmonic,ts,nodal_corrections,astro_angles,averaging_time);
  fprintf(stderr," done.\n");
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  construct the harmonic matrix and RHS vectors
  
*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  parallelfiles=(nt/filelist.size()<nproc*9);
  STDERR_BASE_LINE("parallelfiles=(%d/%d<%d*%d)=%d;\n",nt,filelist.size(),nproc,9,parallelfiles);
  
  struct timeval b4;
  gettimeofday(&b4);
  
#pragma omp parallel for private(i) schedule(dynamic,1) if(nproc>1 && parallelfiles)
  for(i=0;i<filelist.size();i++){//for all files
    int nif;//number of input frames
    date_t origine;
    double *time,od;
#ifdef OMP_H
    int fti=omp_get_thread_num();//file thread index
#endif
#if NETCDF_CAN_PARALLEL_IO == 0
#pragma omp critical(threadUnsafeNetcdf)
#endif
    {
    gettimeofday(&before);
    status= poc_gettime(filelist[i].c_str(), &origine, &time, &nif);
    rrt+=difftime(before);
    }
    od=cnes_time(origine,'s');
    if(isnan(od))od=origd;
    printf("\r(%d/%d;%04.3g)reading %s ...",i+1,filelist.size(),difftime(b4),filelist[i].c_str());fflush(stdout);
    
    /// and frame by frame
#pragma omp parallel for private(ifi,buffer,k) firstprivate(harmonic) schedule(dynamic,1) if(nproc>1 && !parallelfiles)
    for(ifi=0;ifi<nif;ifi++) {//for all frames
      
      int tfi=pos(time[ifi]+od,ts,nt);
      if(tfi<0) /* frame is not in the selected times */
        continue;
      if(skipFrame!=0 and skipFrame[tfi])
        continue;
      
      int bi;//buffer index
#ifdef OMP_H
      if(parallelfiles)
        bi=fti;
      else
        bi=omp_get_thread_num();
      harmonic.rhs=rhS[bi];
      buffer=buffers[bi];
#else
#error not finished
#endif
      
#if NETCDF_CAN_PARALLEL_IO == 0
#pragma omp critical(threadUnsafeNetcdf)
#endif
      {
      for(k=0;k<nvars;k++) {
        const poc_var_t *infok=&info[k];
        printf("\r(%d/%d;%04.3g)reading %s (frame %d/%d) th.%d "+infok->name+"...%s",i+1,filelist.size(),difftime(b4),filelist[i].c_str(),tfi+1,nt,bi,el);fflush(stdout);
        gettimeofday(&before);
        //NOTE: poc_get_vara(int,const V,size_t frame,T*) IS NOT FASTER THAN poc_get_vara(string,const V,size_t frame,T*) !!!
        status=poc_get_vara(filelist[i],*infok,ifi,buffer[k]);
        rt+=difftime(before);
        if(status){NC_CHKERR_BASE_LINE(status,"poc_get_vara(\""+filelist[i]+"\",(\""+infok->name+"\"),,) error");exit(1);}
        }
      
      bool assertedMask=false;
      
      /* setting mask from the values of the first frame */
      if(mask==NULL){
        mask=new bool*[nvars];
        
        for(k=0;k<nvars;k++) {
          int nmasked=0;
          const size_t *nvaluesk=&nvalues[k];
          float
            *bufferk=buffer[k],
            speck=spec[k];
          
          bool *maskk;
          mask[k]=new bool[*nvaluesk];
          maskk=mask[k];
          
#pragma omp parallel for reduction(+:nmasked)
          for(m=0;m<*nvaluesk;m++) {
            bool *maskkm=&maskk[m];
            *maskkm=(bufferk[m]!=speck);
            nmasked+= not *maskkm;
            }
          
          if(nmasked){
            printf("%s"+info[k].name+" %d/%d",assertedMask?", ":"\nmasked: ",nmasked,*nvaluesk);
            assertedMask=true;
            }
          
          }
        
        }
      else if(variable_mask){
        /* Because some nodes may be masked on some but not all the frames */
        
        for(k=0;k<nvars;k++) {
          int moreMasked=0;
          const size_t *nvaluesk=&nvalues[k];
          float
            *bufferk=buffer[k],
            speck=spec[k];
          
          bool *maskk=mask[k];
          
#pragma omp parallel for reduction(+:moreMasked)
          for(m=0;m<*nvaluesk;m++) {
            bool *maskkm=&maskk[m];
            
            if(*maskkm && bufferk[m]==speck){
              *maskkm=0;
              moreMasked++;
              }
            
            }
          
          if(moreMasked){
            printf("%s"+info[k].name+" %d",assertedMask?", ":"more masked: ",moreMasked);
            assertedMask=true;
            }
          
          }
        
        }
      
      if(assertedMask)
        printf("\n");
      
      }//END OF #pragma omp critical(threadUnsafeNetcdf)
      /// calling harmonic_storage(harmonic_t,int,float**) at each frame.
      harmonic.mask=mask;//for harmonic_storage()
      if(nproc==2){
#pragma omp critical(EITHER_read_OR_harmonic_storage)
        {
        omp_set_num_threads(omp_get_num_procs());     //faster
        //omp_set_num_threads(omp_get_num_procs()-1); //than this
        harmonic_storage(harmonic, tfi, buffer);
        }
        }
      else{
        harmonic_storage(harmonic, tfi, buffer);
        }
      }//EO for all frames
    delete[] time;
    }//EO for all files
    
  omp_set_num_threads(omp_get_num_procs());
  
  harmonic.mask=mask;//because it has been done inside a parallelised loop without lastprivate
  printf("\n");
  hst+=difftime(b4);

  //buffers clean-up
  for(i=0;i<nproc;i++){
    deletep2D(&buffers[i],nvars);
    }
  delete[]buffers;
  deletep(&skipFrame);
  
  //rhs merge
  for(i=1;i<nproc;i++){
    int m,n;
    double **rhS0m,**rhSim;//pointers : used for speed
    for(m = 0; m < harmonic.nrhs; m++) {
      rhS0m=rhS[0][m];
      rhSim=rhS[i][m];
      //#pragma omp parallel for private(n,k)
      for(n = 0; n < harmonic.nndes[m]; n++) {
        for(k = 0; k < harmonic.spectrum.n; k++) {
          rhS0m[n][2 * k]    +=rhSim[n][2 * k]    ;
          rhS0m[n][2 * k + 1]+=rhSim[n][2 * k + 1];
          }
        }
      }
    harmonic.rhs=rhS[i];
    harmonic_free_rhs(harmonic);
    }
  harmonic.rhs=rhS[0];

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  compute solution of the harmonic analysis equation
  
*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  gettimeofday(&before);
  constants=harmonic_analysis_core(harmonic, duration,1);
  hat+=difftime(before);
  harmonic_free_rhs(harmonic);
  
  harmonic.spectrum.n=harmonic.spectrum.nmax-1;
  spectrum.n=spectrum.nmax-1;

#if harmonic_analysis_constants_USE_polar
  for(k=0;k<nvars;k++) {
    for(i=0;i<nvalues;i++) {
      constants[k][i].set_complex();
      }
    }
#endif

/*------------------------------------------------------------------------------
  check the control file */
  FILE *in=NULL,*out=NULL;//control input and output files
  if(controlfile!=0) {
    in=fopen(controlfile,"r");
    out=fopen("control.gplt","w");
    if(in==NULL || out==NULL) {
      if(in==NULL)printf("cannot open file %s\n",controlfile);
      if(out==NULL)printf("cannot write to file control.gplt\n");
      if(in!=NULL){fclose(in);in=NULL;}
      if(out!=NULL){fclose(out);out=NULL;}
      }
    }

/*------------------------------------------------------------------------------
  If there is a control file or if atlases are requested */
  poc_grid_data_t *coordinates=0;
  
  if(in!=NULL || atlas>0){
    int coordI;//coordinate index
    
    coordinates=new poc_grid_data_t[nvars];
    
    for(k=0;k<nvars;k++){
      poc_grid_data_t &coordinatesk=coordinates[k];
      
      STDERR_BASE_LINE("#Coordinates of "+info[k].name+" : ");
      
      status=ENOENT;/* No such file or directory */
      
      if(gridfile!=NULL)
        status=poc_get_grid_data(gridfile,info[k],&coordinatesk,0);
      
      if(status!=NC_NOERR)
        status=poc_get_grid_data(filelist[0],info[k],&coordinatesk,0);
      
      if(status!=NC_NOERR){
        NC_CHKERR_BASE_LINE(status,"poc_get_grid_data(,(\"%s\"),) error (see failsafe info below) with both %s and "+filelist[0],varnames[0],gridfile);
        STDERR_BASE_LINE("THE GRID VARIABLES, IF ANY, WILL NOT BE SAVED.\n");
        coordI=0;
        }
      else{
        for(coordI=0;coordI<coordinatesk.vc/2 && coordinatesk[coordI].data;coordI++);
        fprintf(stderr,"%d spacial dimensions.\n",coordI);
        }
      
      n_spacedims[k]=coordI;
      
      if(atlas<=0)
        break;
      }
    
    }
  
/*-------------------------------------------------------------------------------------*/
  if(atlas>0) {/// <h3>If requested to do so, make atlases</h3>
    
#if !harmonic_analysis_constants_USE_polar
    for(k=0;k<nvars;k++) {
      for(i=0;i<nvalues[k];i++) {
        constants[k][i].set_polar(1);
        }
      }
#endif
    
    atlas_(nvars, constants, mask, spectrum, coordinates, info, outputDir,format);
    STDERR_BASE_LINE("#atlases: done\n");
    if(atlas>1){
      return constants;
      }
    }

/*-------------------------------------------------------------------------------------*/
  /// <h3> If there is a control file, initialise the control samples</h3>
  int *controlIs=NULL;
  int nsamples=0;
  if(in!=NULL){
    ///\bug 2011-12-15 Damien Allain : there are only samples for the first variable !!!
    const poc_grid_data_t &coordinates0=coordinates[0];
    
    int coordI;//coordinate index
    double controlCoords[n_spacedims[0]];
    double minDist2,dist2;//SQUARES of distances
    
    fscanf(in,"%d\n",&nsamples);
    controlIs=new int[nsamples+1];
    controlIs[nsamples]=-1;//show end of the array
    for(k=0;k<nsamples;k++) {//for each sample
      printf("point %d/%d:",k,nsamples);
      for(coordI=0;coordI<n_spacedims[0];coordI++){
        fscanf(in,"%lf",&controlCoords[coordI]);
        printf(" %s~%g",coordinates0[coordI].info.name.c_str(),controlCoords[coordI]);
        }
      ///and interpolating to the nearest neighbour
      minDist2=INFINITY;
      controlIs[k]=-1;
      for(i=0;i<nvalues[0];i++){
        dist2=0.;
        if(mask[i]==0)continue;
        for(coordI=0;coordI<n_spacedims[0];coordI++){
          ///\bug 2011-12-15 Damien Allain : considering the coordinates0 are rectangular!
          dist2+=pow(controlCoords[coordI]-coordinates0[coordI].data[i],2);
          }
        if(minDist2>dist2){
          minDist2=dist2;
          controlIs[k]=i;
          }
        }
      printf(" => index=%d(%g)",controlIs[k],sqrt(minDist2));
      fprintf(out,"p='%d' # point %d:",controlIs[k],k);
      if(controlIs[k]<0 || nvalues[0]<=controlIs[k]){printf(" ***ERROR : %d out[0;%d[***\n",controlIs[k],nvalues[0]);fprintf(out,"\n");break;}
      for(coordI=0;coordI<n_spacedims[0];coordI++){
        printf(" %s=%g",coordinates0[coordI].info.name.c_str(),coordinates0[coordI].data[controlIs[k]]);
        fprintf(out," %s=%g",coordinates0[coordI].info.name.c_str(),coordinates0[coordI].data[controlIs[k]]);
        }
      printf("\n");fprintf(out,"\n");
      }//EO for each sample
    fprintf(out,"v='%s' # variable name\n",info[0].name.c_str());
    fclose(in);
    fclose(out);
    }
  
/*-------------------------------------------------------------------------------------*/
  /// <h3>finally call detide().</h3>
  gettimeofday(&b4);
  status=detide(filelist, info, nt,ts, harmonic, constants,controlIs,origd,outputDir);
  hct+=difftime(b4);

  delete[] info;
  delete[] nvalues;
  delete[] n_spacedims;
  return constants;
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
    "  %s file1 [ file2 ... ] [OPTIONS] -d wave1 [ wave2 ... ]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Detides comodo-compliant NetCDF outputs and produces tidal atlases.\n"
    "  It takes a file or a list of files as input, carries out a spectral analysis on the given list of wave and produces the atlases and the detided output.\n"
    "  The time variable must comply with the CF conventions : " CF_CONVENTIONS_URL "#time-coordinate\n"
    "\n"
    "OPTIONS\n"
    "  --help,-h : show this help and exit\n"
    "  --nodal=no : do not do nodal corrections\n"
    "  -c,--control : followed by the path of the list of control points: it is an ascii file with the number of control points followed by their coordinates (longitude latitude [layer]).\n"
    "  --only-atlases : produce only atlases : no detiding. It of course implies -a\n"
    "  --load-atlases : only load atlases : no analysis. Currently, this also de-activates --control option.\n"
    "  --drying : take into account that the mask MAY vary A LITTLE at each frame. This takes a tiny bit of CPU, without taking any extra time when the speed is limited by the hard drive.\n"
    "  --flagged : take into account that the mask WILL vary A LOT at each frame. Currently equivalent to --drying.\n"
    "  --take-first : if some time frames are simultaneous, take the first one\n"
    "  --take-last : if some time frames are simultaneous, take the last one\n"
    "  --band-cut : followed by 2 frequencies in deg/h\n"
    "  --averaging : NOT DOCUMENTED. Default : 0\n"
    "  -a : produce atlas\n"
    "  -1 : put all atlases in one file. YOU WILL NOT BE ABLE TO USE comodo-admittance ON THIS FILE.\n"
    "  -t : show separation tables of all harmonics and list of files within time boundaries\n"
    "  -l : followed by the path of the list of files to analyse. This list will override the list given as arguments.\n"
    "  -g : followed by the path of the grid file. This is only necessary when the coordinates are not available in the files to analyse and you want to produce atlases or use control points.\n"
    "  -p : followed by an output folder path. IT IS STRONGLY ADVISED TO OUTPUT ON A DIFFERENT HARD DRIVE WHEN DOING DETIDING.\n"
    "  -s : followed by the start date. See DATE FORMATS below\n"
    "  -f : followed by the end date. See DATE FORMATS below\n"
    "  -o : followed by the default date origin. See DATE FORMATS below\n"
    "  -v : followed by the names of the variables to detide\n"
    "  -d : followed by the list of waves to analyse. A good start is Q1 O1 P1 K1 N2 M2 S2 K2 M4 MS4\n"
    );
  print_poctime_scan_date_help(0);
  print_OPENMP_help(prog_name); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// function that will be ran when the executable is started
/** See the print_help <a href=#func-members>function</a> for its use.
For example :
\code comodo-detidor -s 15/12/2008 -f 28/02/2009 -l files.list-test -control control.dat -g grille_domconcat.nc -v ssh tem sal -a \
  -d M2 S2 N2 K2 2N2 La2 O1 K1 P1 Q1 M4 S4 N4 MS4 MN4 M3 M6 S1 MSK2 2SM2 SN4 MNS2 MSN2 MK3 2SN2 MK4 2NS2 ST2 ST3 SNK2 \endcode
\code comodo-detidor -s 02/01/2000  -control control.dat -v XE -a -d K1 O1 Q1 M2 S2 N2 K2 L2 2N2 2NS2 2SM2 2SN2 MK3 M4 SN4 MS4 3MS4 MN4 S4 M6 2MS6 \endcode

See also \code ncks -s "%.12g\n" -v time,latitude,longitude,XE -d ni,217 -d nj,181 champs_ionc4.nc \endcode
See also \code ncks -mv time detided-XE-champs_ionc4.nc \endcode
*/
/*----------------------------------------------------------------------------*/
{
  double tau;
  int i,j,k;//2 wave indexes, variable index
  int n,status;//argument index or number of files, NetCDF status
  char *gridfile=NULL,*controlfile=NULL,*keyword,*outputDir=NULL;
  int format=WAVE_AS_SEPARATE_FILE;
  char *list=0;
  char *onde[100];
  int nonde=0;
  date_t start=NADate,final=NADate,orig=NADate;
  double startd,finald,origd,*ts;
  size_t tn;
  spectrum_t WaveList;
  int nodal_corrections=1,variable_mask=0,atlas=0,separationTables=0,
    handleRepeatedFrames=0;
  astro_angles_t astro_angles;
  harmonic_t harmonic;
  hconstant_t **constants;
  int nvars=0;
  char **varnames=NULL;
  vector<string> filelist;
  double bcf0=NAN,bcf1=NAN; /*< band-cut frequencies */
  double averaging_time=0.0;
  
  struct timeval mainbefore;
  
  #if 0
  #define nTasks 32
  double cts[nTasks];
  int threadIs[nTasks];
  int ncpu=omp_get_max_threads();
  int bla=-1;
  struct timespec b4;
  gettimeofday(&before);
  printf("%d threads max, nesting : %d\n",ncpu,omp_get_nested());
  #if 1
  //clock_gettime(CLOCK_REALTIME,&b4);
  #pragma omp parallel for private(i)
  for(i=0;i<nTasks;i++){
    cts[i]=difftime(before);
    threadIs[i]=omp_get_thread_num();
    }
  for(i=0;i<nTasks;i++){
    printf("%d %g\n",threadIs[i],cts[i]);
    }
  printf("\n");
  #endif
  #ifdef omp_get_ancestor_thread_num
  bla=1;
  #pragma omp parallel if(0) {
  #pragma omp parallel for if(bla) private(i,j)
  for(i=0;i<ncpu;i++){
    printf("%d:%d\n",i,omp_get_thread_num());
    #pragma omp parallel for if(!bla) private(j)
    for(j=0;j<ncpu;j++){
      printf("%d,%d:%d,%d => %d %d\n",i,j,omp_get_ancestor_thread_num(omp_get_level()-1),omp_get_thread_num(),omp_get_active_level(),omp_get_level());
      }
    }
  }
  printf("\n");
  #endif
  STDOUT_BASE_LINE("%d %d %d\n",sizeof(void*),sizeof(double),sizeof(float));
  STDOUT_BASE_LINE("%d\n",_OPENMP);
  TRAP_ERR_EXIT(ENOEXEC,"testing\n");
  #endif
  #if _OPENMP < 200505  // Version 2.5 May 2005
  #error Check whether your version _OPENMP (below) of openmp can do nesting properly. If so, update the above line accordingly.
  #pragma message "_OPENMP=" TOSTRING(_OPENMP)
  #endif

  fct_echo( argc, argv);
  argv0=argv[0];

  i=0;
  onde[i]=NULL;

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if( strcmp(keyword,"--help")==0 or
        strcmp(keyword,"-h")==0 ) {
      print_help(argv[0]);
      exit(0);
      }
    else if(strcmp(keyword,"--nodal=no")==0) {
      nodal_corrections=0;
      n++;
      continue;
      }
    else if(strcmp(keyword,"--drying")==0) {
      variable_mask=1;
      n++;
      continue;
      }
    else if(strcmp(keyword,"--flagged")==0) {
      variable_mask=2;
      n++;
      continue;
      }
    else if(strcmp(keyword,"--averaging")==0) {
      int nitems;
      nitems=sscanf(argv[n+1],"%lf",&averaging_time);
      n++;
      n++;
      continue;
      }
    else if(strcmp(keyword,"--control")==0) {
      if(controlfile!=NULL)
        free(controlfile);
      controlfile=strdup(argv[n+1]);
      n++;
      n++;
      continue;
      }
    else if(strcmp(keyword,"--only-atlases")==0) {
      atlas|=2;
      n++;
      continue;
      }
    else if(strcmp(keyword,"--load-atlases")==0) {
      atlas|=4;
      n++;
      continue;
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
    else if(strcmp(keyword,"--band-cut")==0) {
      bcf0=atof(argv[n+1]);
      bcf1=atof(argv[n+2]);
      n++;
      n++;
      n++;
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 'c' :
          if(controlfile!=NULL)
            free(controlfile);
          controlfile=strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'a' :
          atlas|=1;
          n++;
          break;

        case 't' :
          separationTables=1;
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

        case 'p' :
          outputDir= strdup(argv[n+1]);  /* directory */
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

        case 'd' :
          n++;
          i=pos((char*)NULL,onde,100);
          for(;n<argc;i++,n++) {
            onde[i]= strdup(argv[n]);
            }
          onde[i]=NULL;
          nonde=i;
          break;

        case '1' :
          n++;
          format=WAVE_AS_DIMENSION;
          break;

        default:
          printf("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
        filelist.push_back(argv[n]);
        n++;
        break;
      }
      free(keyword);
    }
  
  if(outputDir==NULL) outputDir=strdup(".");
  status=mkdir(outputDir,0777);
  if(status!=0 && errno!=EEXIST) TRAP_ERR_EXIT(errno,"mkdir(\"%s\",) error (%d %s)\n",outputDir,errno,strerror(errno));
  
  status=0;
  
  printf("# number of waves to analyse: %d\n",nonde);
  if(nonde==0) {
    STDOUT_BASE_LINE(" *** Please give at least one wave ***\n");
    status=-1;
    }
  
  printf("#number of variables : %d\n",nvars);
  if(nvars==0) {
    STDOUT_BASE_LINE(" *** Please give the names of the variables to filter ***\n");
    status=-1;
    }
  for(k=0;k<nvars;k++) {
    printf("#variable %d: %s\n",k,varnames[k]);
    }

  if(list) {
    filelist=load_filelist(list);
    }

  if(filelist.size()==0) {
    if(list)
      STDOUT_BASE_LINE(" *** the file list %s is empty ***\n",list);
    STDOUT_BASE_LINE(" *** Please give a NetCDF file to analyse ***\n");
    status=-1;
    }
  
  if( (atlas & 4) and atlas!=4) {
    STDOUT_BASE_LINE(" *** --load-atlases is incompatible with both -a and --only-atlases ***\n");
    status=-1;
    }
  
  if(status!=0) {
    print_help(argv[0]);
    wexit(-1);
    }

  if(gridfile!=NULL && access(gridfile,R_OK)) {//NOTE: access rights may change during the run
    //This is to avoid crashing without output after doing the spectral analysis
    TRAP_ERR_EXIT(errno,"no read access to %s (%d %s)\n",gridfile,errno,strerror(errno));
    }
  
  gettimeofday(&mainbefore);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  check list of files against start and end dates
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  startd=cnes_time(start,'s');
  finald=cnes_time(final,'s');
  origd=cnes_time(orig,'s');
  printf("Selected data from %s to %s.\n",
    poctime_sdate_cnes(startd,'s'),poctime_sdate_cnes(finald,'s'));
  n=filelist.size();
  tn=poc_timefilterlist(&filelist,&startd,&finald,&ts,origd);
  if(separationTables){
    printf("Kept files :");
    for(i=0;i<filelist.size();i++){
      printf(" %s",filelist[i].c_str());
      }
    }
  printf("\n");
  printf("Kept %d/%d files with %d frames from %s to %s.\n\n",filelist.size(),n,tn,
    poctime_sdate_cnes(startd,'s'),poctime_sdate_cnes(finald,'s'));
  n=filelist.size();
  if(tn==0)
    TRAP_ERR_EXIT(1,"no frames within time boundaries\n");

  start=poctime_getdatecnes(startd,'s');
  init_argument(&astro_angles,start,0);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  init data for harmonic analysis
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
  const double duration=(finald-startd)/3600.;
  double df=NAN;
  
  if( isfinite(bcf0) and isfinite(bcf1) ){
    if(bcf1<bcf0)
      swapValues(&bcf1,&bcf0);
    df=360./duration;
    WaveList.nmax=(bcf1-bcf0)/df+1;
    }

  WaveList.init(initialize_tide(),onde,-WaveList.nmax);
  printf("# number of waves to analyse: %d \n",WaveList.n);
  
  if( isfinite(df) ){
    tidal_wave w;
    for(w.omega=bcf0;w.omega<bcf1;w.omega+=df){
      sprintf(w.name,"%.9gdph",w.omega);
      status=WaveList.add(w,0);
      if(status<0) TRAP_ERR_EXIT(ENOEXEC,"spectrum_t::add() error\n");
      }
    }
  
  if(WaveList.n==0) {
    STDOUT_BASE_LINE(" *** Please give at least one VALID wave ***\n");
    print_help(argv[0]);
    exit(-1);
    }

  if(separationTables){
    printf("Pulsations and frequencies of waves:\n");
    for (i=0; i<WaveList.n; i++) {
      if(WaveList.waves[i].omega == 0.) {
        WaveList.waves[i].init();
        }
      printf ("%10s %12.6f deg/h %12.9f h-1 %g\n", WaveList.waves[i].name,WaveList.waves[i].omega,WaveList.waves[i].omega/360.,greenwhich_argument(astro_angles,WaveList.waves[i]));
      }
  
    printf("\nWaves separations in days :\n\n");
    for (i=0; i<WaveList.n-1; i++) {
      printf ("%10s %9s %10s %9s\n",
        "wave",WaveList.waves[i].name,
        "wave",WaveList.waves[i].name);
      k=0;
      for (j=i+1; j<WaveList.n; j++) {
        tau=deltaOmega2separation(WaveList.waves[j].omega-WaveList.waves[i].omega);
        printf ("%10s %9.3f ",WaveList.waves[j].name,tau);
        if(k){
          printf("\n");
          k=0;
          }
        else
          k=1;
        }
      if(k){
        printf("\n");
        }
      printf ("\n");
      }
    TRAP_ERR_EXIT(-1,"exiting\n");
    }

//   for(k=0;k<nvars;k++) {
//     constants=analyse((const char**) &(varnames[k]),1,filelist,gridfile,controlfile, start,final,WaveList,nodal_corrections,atlas);
//     }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  harmonic analysis and detiding
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(atlas==4){ /* NO analysis : detiding only */
    
    /* load atlases and variable headers */
    poc_global_t global;
    poc_var_t *info=new poc_var_t[nvars];
    status=poc_inq(filelist[0],&global,0);
    
    constants=new hconstant_t *[nvars];
    
    size_t *nvalues=new size_t[nvars];
    harmonic.mask=new bool *[nvars];
    
    for(k=0;k<nvars;k++){
      const string varName=varnames[k];
      const string fileConvention=atlas_file_name(outputDir,format,"WAVE",varName);
      const string aName=varName+AMPSUF;
      const string gName=varName+PHASUF;
      
      constants[k]=load_atlas(fileConvention,aName,gName,WaveList);
      
      status=poc_get_var_length(global,varName,&info[k],&nvalues[k]);
      NC_TRAP_ERROR(wexit,status,1,"poc_get_var_length(,\"%s\",...) error in "+filelist[0]+"\n",varnames[k]);
      
      hconstant_t *constantsk=constants[k];
      
      const size_t nvaluesk=nvalues[k];
      bool *maskk=new bool[nvaluesk];
      
      for(i=0;i<nvaluesk;i++){
        hconstant_t *constantski=&constantsk[i];
        constantski->set_complex(true);
        maskk[i]= not isnan(constantski->z[0]);
        }
      
      harmonic.mask[k]=maskk;
      }
    
    /* initialise harmonic */
    harmonic.spectrum=WaveList;
    harmonic.neq   =2*WaveList.n;
    
    harmonic.nndes=copy(nvalues,nvars);
    
    harmonic.nrhs=nvars;
    
    harmonic_start(&harmonic,tn);
    
    STDERR_BASE_LINE("Initialising harmonic matrix for %d waves...",WaveList.n);
    harmonic_init(&harmonic,ts,nodal_corrections,astro_angles,averaging_time);
    fprintf(stderr," done.\n");
    
    /* detide */
    status=detide(filelist, info, tn,ts, harmonic, constants,0,origd,outputDir);
    
    }
  else{
    constants=analyse((const char**) varnames,nvars,filelist,gridfile,controlfile, tn,ts,WaveList,nodal_corrections, averaging_time,variable_mask,astro_angles,atlas,origd,outputDir,format,handleRepeatedFrames);
    
    if(constants==NULL) TRAP_ERR_EXIT(-1,"analyse failed\n");
    }
  
  STDOUT_BASE_LINE(" ^^^^^^^^^^^^^ end of computation, that took %#g s with rt=%#g,rt2=%#g,rrt=%g,wt=%#g,hst=%#g,hat=%#g,detide=%#g ^^^^^^^^^^^^^\n",difftime(mainbefore),rt,rt2,rrt,wt,hst,hat,hct);
  return 0;
}
