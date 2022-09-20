
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
\brief Filters comodo-compliant NetCDF outputs and produces energy atlases.

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
#include "tides.def"
#include "functions.h" //safely includes omp.h
#include "matrix.h"
#include "poc-netcdf-data.hpp"
#include "poc-time.h"
#include "filter.h"
#include "map.h"
#include "maths.h"

struct timeval stv;///<start timeval, for progression
///These global variables are used to count the amount of time this programme took for I/O tasks
struct timeval before;
double rt=0.,rt2=0.,rrt=0.,wt=0.;// read, 2nd read, redundant read and write times
double ft=0.;// filtering time

string cmd;///<command line
size_t nt;///<number of samples

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int atlas_(int nvars, double **RMS, const bool *mask, poc_grid_data_t *coordinates, const poc_var_t *info, const char *outputPrefix)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,i,li,nsd,status;
  size_t length;
  poc_global_t global("constructed around " __LINE_FILE_PACKAGE_REVISION);
  poc_var_t a;
  poc_dim_t wdim("waves",NC_UNLIMITED);
  poc_att_t *att;
  string s="";
  poc_deque_t<string> outVarNames;
  char *coordinateNames;
  #if 0
  #define T double
  #define NC_T NC_DOUBLE
  #define NC_FILL_T NC_FILL_DOUBLE
  #else
  #define T float
  #define NC_T NC_FLOAT
  #define NC_FILL_T NC_FILL_FLOAT
  #endif
  double scale;//scale factor
  T spec;//fill value

/*------------------------------------------------------------------------------
  compose header */
  outVarNames.clear(nvars);
  
  for(k=0;k<nvars;k++){
    const poc_var_t *var=&info[k];
    
    ///It initialises the file with the coordinates variables
    status=poc_get_var_length(*var,&length,&nsd);
    STDERR_BASE_LINE("***WARNING : nsd=%d IS THIS THE RIGHT VALUE ???***\n",nsd);
    for(i=0;i<nsd;i++){
      if(coordinates[k][i].info.name!="")
        global.variables << coordinates[k][i].info;
      s+=(string)" "+coordinates[k][i].info.name;
      }
    coordinateNames=poc_strdup(s.substr(1).c_str());
    ///and with the RMS variable
    s=var->name+"_RMS";
    outVarNames.push_back(s);
    a.init(s,NC_T);
    li=0;//set to 1 if ou want to have wave as a dimension
    switch(li){
    case 1:
      a.dimensions << wdim;
    case 0:
      break;
    default:
      CHKERR_LINE_FILE(ENOEXEC,"programming error");exit(-2);
      }
    for(i=0;i<var->dimensions.size();i++){
      if(isT(var->dimensions[i]))continue;
      a << var->dimensions[i];
      }
    status=poc_decode_mask(*var,&scale,NULL,&spec,-NC_FILL_T);
    a << poc_att_t("_FillValue",spec);
    ///setting the scale factor so that the user can toggle it
    a << poc_att_t("scale_factor",scale);
    a << poc_att_t("coordinates",coordinateNames);
    delete[]coordinateNames;
    for(i=0;i<var->attributes.size();i++){
      att=(poc_att_t*)&var->attributes[i];
      if(att->name=="long_name"){
        a << poc_att_t(att->name,"RMS of filtered "+att->as_string());
        continue;
        }
      if(att->name=="standard_name"){
        a << poc_att_t(att->name,"root_mean_square_of_filtered_"+att->as_string());
        continue;
        }
      if(att->name=="units"){
        a << poc_att_t(*att);
        continue;
        }
      if(att->name=="associate" ||
        att->name=="coordinates" || //already specified above
        att->name=="content" ||
        att->name=="valid_min" ||
        att->name=="valid_max" ||
        att->name=="valid_range" ||
        ///\bug 2011-12-14 Damien Allain : this does not put add_offset attribute EVEN FOR LOW-PASS-FILTERED DATA (a feature not available yet)
        att->name=="add_offset" || //do not put add_offset for RMS of non-low-pass-filtered data
        att->name=="_FillValue" || //already specified by default
        att->name=="scale_factor" || //already specified by default
        att->name=="missing_value"){
        continue;
        }
      a << *att;
      }
    global.variables << a;
    //cdf_print_varinfo(g);
    }
  global << poc_att_t("number_of_samples",(int32_t)nt);
  global << poc_att_t("RMS_normalisation_scale_factor",RMS[nvars][0]);
  global << poc_att_t("history",cmd);

  
/*------------------------------------------------------------------------------
  create file and write header */
  if(outputPrefix!=0)s=outputPrefix;
  else s="";
  s+=info[0].name+"-filtered.nc";
  STDERR_BASE_LINE("(over)writing "+s+" ...");
  unlink(s.c_str());
  status=poc_create(s,global,1);
  if(status!=NC_NOERR)NC_TRAP_ERROR(wexit,status,1,"poc_create(\""+s+"\",) error");
  
/*------------------------------------------------------------------------------
  write data */
  for(k=0;k<nvars;k++){
    ///it writes the coordinates
    for(i=0;i<nsd;i++){
      if(coordinates[k][i].info.name=="")
        continue;
      fprintf(stderr," "+coordinates[k][i].info.name+",");
      status=coordinates[k][i].write_data(s);
      NC_CHKERR_LINE_FILE(status,"poc_data_t::write_data(\""+s+"\") error with "+coordinates[k][i].info.name);
      }
    ///and it writes the data
    fprintf(stderr," "+outVarNames[k]+",");
    status=poc_put_vara(s,outVarNames[k],0,RMS[k]);
    NC_CHKERR_LINE_FILE(status,"poc_put_vara() error with "+outVarNames[k]+" on "+s);
    }
  fprintf(stderr," done.\n");
  #undef T
  #undef NC_T

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double **filter(const char **varnames, int nvars, vector<string> filelist, char *gridfile, char *controlfile, size_t nt,double *ts, double frequency, double q, int atlas,double origd=NAN, char *outputPrefix=NULL)

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
\param frequency pulsation of the central frequency in degrees per hour. If NAN, take inertial frequency.
\param q Q factor
\param atlas whether an atlas is done
\param origd=NAN if specified, default time origin if not in the files
\param *outputPrefix=NULL if specified, path to the directory where the output files will be written. IT IS STRONGLY ADVISED TO OUTPUT ON A DIFFERENT HARD DRIVE WHEN DOING FILTERING.
\returns the RMS of the filtered data, or NULL if error

\date 2011-09-19 Damien Allain : added frame times check
*/
/*----------------------------------------------------------------------------*/
{
  float **inBuf,*outBuf,
    **B0,**B1,**B2,**B20;
  double **RMS;
  double duration,startd,finald;
  int j,k,i;//variable and file indexes
  int status;
  int ifi,dfi;//input frame index,duplicate frame index
  poc_global_t global;
  poc_var_t *info;
  float *spec;
  bool *mask=NULL;
  size_t *nvalues,nMax,ntotal;
  int *n_spacedims;
  
  double sampling,w0,b0,b1,b2,a1,a2;//filter parameters
  
  if(nvars<=0)return NULL;
  nvalues=new size_t[nvars+1];//see normalisation
  n_spacedims=new int[nvars];
  spec=new float[nvars];
  
  ///If there are duplicates in the frame times given, it returns NULL.
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
    return NULL;
    }
  
  ///It calculates the duration of the analysis from the times of the frames
  minpos(ts,nt,&startd);
  maxpos(ts,nt,&finald);
  duration=finald-startd;
  
  /// <h3>It gets variable dimensions</h3>
  info=new poc_var_t[nvars];
  gettimeofday(&before);
  status=poc_inq(filelist[0],&global,0);
  
  poc_var_t timeVar;
  timeVar.init("time",NC_DOUBLE,"time","seconds since 1950-01-01 00:00");
  for(k=0;k<global.dimensions.size();k++){
    poc_dim_t *dim=&global.dimensions[k];
    if(not isT(*dim))
      continue;
    timeVar << *dim;
    break;
    }
  
  rrt+=difftime(before);
  ntotal=0;
  nMax=1;
  for(k=0;k<nvars;k++){
    status=poc_get_var_length_and_mask(global,varnames[k],&info[k],&nvalues[k],&n_spacedims[k],&spec[k]);
    if(status==NC_ENOTATT)NC_CHKERR_LINE_FILE(status,"Use default mask value of %g for %s in "+filelist[0],spec[k],varnames[k]);/* warning only */
    else if(status!=NC_NOERR){NC_CHKERR_LINE_FILE(status,"poc_get_var_length_and_mask(,\"%s\",...) error in "+filelist[0]+"\n",varnames[k]);return NULL;}
    
    updatemax(&nMax,nvalues[k]);
    ntotal+=nvalues[k];
    
    unlimitTimeDim(&info[k]);
    }
  nvalues[nvars]=1;//see normalisation
  
/*-------------------------------------------------------------------------------------*/
  if(isnan(frequency)){
    const tidal_wave K2=wK2;
    astro_angles_t angles;
    init_argument(&angles,0.);
    TRAP_ERR_EXIT(ENOEXEC,"not finished\n");
    }
  
/*-------------------------------------------------------------------------------------*/
  /* initialize filter */
  sampling=duration/(nt-1);
  w0=frequency*dph2rps*sampling;//central frequency in rad/sample
  STDOUT_BASE_LINE("sampling period : %gs;central frequency : %gdeg/h =%g= %g sampling periods;q=%g\n",sampling,frequency,pi2/w0,1/(frequency*dph2rps/pi2*sampling),q);
  b0=w0/(2*q);b1=0;b2=-b0;a1=2*(1-b0)*cos(w0);a2=2*b0-1; // Band pass
  STDOUT_BASE_LINE("b0=%g,b1=%g,b2=%g,a1=%g,a2=%g\n",b0,b1,b2,a1,a2);
  
  /* allocate buffers */
  inBuf=new float*[nvars+1];
  B0=new float*[nvars+1];
  B1=new float*[nvars+1];
  B2=new float*[nvars+1];
  RMS=new double*[nvars+1];
  for(k=0;k<=nvars;k++){
    exitIfNull(inBuf[k]=new float[nvalues[k]]);
    exitIfNull(B0[k]=aset(nvalues[k],0.f));
    exitIfNull(B1[k]=aset(nvalues[k],0.f));
    exitIfNull(B2[k]=aset(nvalues[k],0.f));
    exitIfNull(RMS[k]=aset(nvalues[k],0.));
    }
  exitIfNull(outBuf=new float[nMax]);
  
  /* initialise output file */
  string outPath;
  if(outputPrefix!=0)
    outPath=outputPrefix;
  outPath+="filtered.nc";
  if(atlas<2){
    if(access(outPath.c_str(),W_OK)!=0){
      printf("creating "+outPath+"\n");
      poc_global_t global("constructed around " __LINE_FILE_PACKAGE_REVISION);
      global << poc_att_t("history",cmd);
      global << timeVar;
      status=poc_create(outPath,global,1);
      NC_CHKERR_LINE_FILE(status,"poc_create(\""+outPath+"\",,) error");
      }
    else
      printf(outPath+" already created\n");
    }
  
  FILE *f=fopen("filter.dat","w");
  
/*-------------------------------------------------------------------------------------*/
  /* filter */
  struct timeval b4;
  gettimeofday(&b4);
  
  for(i=0;i<filelist.size();i++){//for all files
    int nif;//number of input frames
    date_t origine;
    double *time,od;
    gettimeofday(&before);
    status= poc_gettime(filelist[i].c_str(), &origine, &time, &nif);
    rrt+=difftime(before);
    od=cnes_time(origine,'s');
    if(isnan(od))od=origd;
    //printf("\r(%d/%d;%04.3g)reading "+filelist[i]+" ...",i,filelist.size(),difftime(b4));fflush(stdout);
    
    for(ifi=0;ifi<nif;ifi++) {//for all frames
      int tfi=pos(time[ifi]+od,ts,nt);
      if(tfi<0)
        continue;
      ///\bug 2011-09-02 Damien Allain : assuming variables have the same sizes !!!
      if(timeIsOld()){printf("\r(%d/%d;%04.3g)reading "+filelist[i]+" (frame %d/%d,%g,%g,%g) ",i,filelist.size(),difftime(b4),tfi+1,nt,rt,ft,RMS[nvars][0]);fflush(stdout);}
      
      for(k=0;k<nvars;k++){
        gettimeofday(&before);
        status=poc_get_vara(filelist[i],info[k],ifi,inBuf[k]);
        rt+=difftime(before);
        if(status){NC_CHKERR_LINE_FILE(status,"poc_get_vara(\""+filelist[i]+"\",,%d,) error with "+info[k].name,ifi);exit(1);}
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
        STDOUT_BASE_LINE("%d / %d masked : %d",nmasked,nvalues[0],mask[nvalues[0]/2]);
        }
      
      /// delaying registers
      B20=B2;B2=B1;B1=B0;B0=B20;
      
      /// and doing the weighted sums
      gettimeofday(&before);
      for(k=0;k<=nvars;k++){
        float *inBufk=inBuf[k],*B0k=B0[k],*B1k=B1[k],*B2k=B2[k];
        double *RMSk=RMS[k];
        
        #pragma omp parallel for private(j)
        for(j=0;j<nvalues[k];j++){
          /* skip mask values */
          if(inBufk[j]==spec[k]){
            outBuf[j]=spec[k];
            continue;
            }
          
          ///\f$ B0=B1*a1+B2*a2+inBuf \f$
          //note : j may well be optimised out here
          B0k[j]=B1k[j]*a1+B2k[j]*a2+inBufk[j];
          ///\f$ O=B0*b0+B1*b1+B2*b2 \f$
          outBuf[j]=B0k[j]*b0+B1k[j]*b1+B2k[j]*b2;
          ///\f$ RMS+=O**2 \f$
          RMSk[j]+=square((double)outBuf[j]);
          }
        
        if(k==nvars){
          fprintf(f,"%g %g %g %g %g %g\n",inBufk[0],B0k[0],B1k[0],B2k[0],outBuf[0],RMSk[0]);
          }
        else if(atlas<2){
          if(timeIsOld()){printf("\r(%d/%d;%04.3g)writing "+outPath+" (frame %d/%d) ",i,filelist.size(),difftime(b4),tfi+1,nt);fflush(stdout);}
          status=poc_put_vara(outPath,info[k],tfi,outBuf);
          NC_CHKERR_LINE_FILE(status,"poc_put_vara(\""+outPath+"\",(\""+info[k].name+"\"),%d,) error",tfi);
          }
        
        }
      
      if(atlas<2){
        status=poc_put_vara(outPath,timeVar,tfi,&ts[tfi]);
        NC_CHKERR_LINE_FILE(status,"poc_put_vara(\""+outPath+"\",,%d,) error with "+timeVar.name,tfi);
        }
      
      ft+=difftime(before);
      }//EO for all frames
    
    delete[] time;
    }//EO for all files
  
  printf("\n");
  fclose(f);

  //buffers clean-up
  for(k=0;k<nvars;k++){
    float *inBufk=inBuf[k];
    double *RMSk=RMS[k];
    
    for(j=0;j<nvalues[k];j++){
      
      if(inBufk[j]==spec[k]){
        RMSk[j]=spec[k];
        continue;
        }
      
      RMSk[j]=sqrt(RMSk[j]/nt);
      }
    delete[]inBufk;
    delete[]B0[k];
    delete[]B1[k];
    delete[]B2[k];
    }
  STDOUT_BASE_LINE("energy of response to impulse : %g\n",RMS[nvars][0]);
  RMS[nvars][0]=pow(RMS[nvars][0]*nt,-.5)*2;
  STDOUT_BASE_LINE("RMS normalisation factor : %g\n",RMS[nvars][0]);
  delete[]inBuf;
  delete[]B0;
  delete[]B1;
  delete[]B2;
  
/*-------------------------------------------------------------------------------------*/
  /// <h3> It checks the control file</h3>
  FILE *in=NULL,*out;//control input and output files
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

  /// <h3> If there is a control file or if atlases are requested, it reads the coordinates</h3>
  poc_grid_data_t *coordinates=NULL;
  if(in!=NULL || atlas>0){
    STDERR_BASE_LINE("#atlases:\n");
    coordinates=new poc_grid_data_t[nvars];
    for(k=0;k<nvars;k++){
      status=ENOENT;/* No such file or directory */
      if(gridfile!=NULL)
        status=poc_get_grid_data(gridfile,info[k],&coordinates[k]);
      if(status!=NC_NOERR)
        status=poc_get_grid_data(filelist[0],info[k],&coordinates[k]);
      if(status!=NC_NOERR){NC_CHKERR_LINE_FILE(status,"poc_get_grid_data() error on %s with both %s and "+filelist[0],varnames[k],gridfile);exit(-1);}
      }
    }

  /// <h3>If requested to do so, it makes atlases</h3>
  if(atlas>0) {
    ///by calling atlas_()
    status= atlas_(nvars, RMS, mask, coordinates, info, outputPrefix);
    STDERR_BASE_LINE("#atlases: done\n");
    if(atlas>1){
      delete[]coordinates;
      return RMS;
      }
    }

  /// <h3> If there is a control file, it initialises the control samples</h3>
  ///\bug 2011-12-15 Damien Allain : there are only samples for the first variable !!!
  int *controlIs=NULL;
  int nsamples=0;
  if(in!=NULL){
    ///by reading the sample coordinates
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
        printf(" %s~%g",coordinates[0][coordI].info.name.c_str(),controlCoords[coordI]);
        }
      ///and interpolating to the nearest neighbour
      minDist2=INFINITY;
      controlIs[k]=-1;
      for(i=0;i<nvalues[0];i++){
        dist2=0;
        if(mask[i]==0)continue;
        for(coordI=0;coordI<n_spacedims[0];coordI++){
          ///\bug 2011-12-15 Damien Allain : considering the coordinates are rectangular!
          dist2+=pow(controlCoords[coordI]-coordinates[0][coordI].data[i],2);
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
        printf(" %s=%g",coordinates[0][coordI].info.name.c_str(),coordinates[0][coordI].data[controlIs[k]]);
        fprintf(out," %s=%g",coordinates[0][coordI].info.name.c_str(),coordinates[0][coordI].data[controlIs[k]]);
        }
      printf("\n");fprintf(out,"\n");
      }//EO for each sample
    fprintf(out,"v='%s' # variable name\n",info[0].name.c_str());
    fclose(in);
    fclose(out);
    }
  if(coordinates!=NULL){
    delete[]coordinates;
    }
  
  delete[] info;
  delete[] nvalues;
  delete[] n_spacedims;
  return RMS;
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
    "  Filters comodo-compliant NetCDF outputs and produces energy atlases.\n"
    "  It takes a file or a list of files as input, carries out a DF-II IIR filtering and produces the energy atlases and the filtered output.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  Show this help and exit.\n"
    "  --control followed by the path of the list of control points: it is an ascii file with the number of control points followed by their coordinates (longitude latitude [layer]).\n"
    "  --only-atlases produce only atlases : do not output the filtered data\n"
    "  -a  produce atlases\n"
    "  -c  followed by convention CURRENTLY WITHOUT EFFECT!\n"
    "  -l  followed by the path of the list of files to analyse. This list will override the list given as arguments.\n"
    "  -g  followed by the path of the grid\n"
    "  -p  followed by an output prefix. It can be a folder path followed by /. IT IS STRONGLY ADVISED TO OUTPUT ON A DIFFERENT HARD DRIVE WHEN DOING FILTERING.\n"
    "  -s  followed by the start date in dd/mm/yyyy format\n"
    "  -f  followed by the end date in dd/mm/yyyy format\n"
    "  -o  followed by the default date origin in dd/mm/yyyy format\n"
    "  -v  followed by the list of variables to filter\n"
    "  -w  followed by the frequency at which to filter\n"
    "  -q  followed by the Q factor\n"
    "\n"
    "BUG\n"
    "  When giving more than one variable, atlases are not produced.\n");
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
  char *gridfile=NULL,*controlfile=NULL,*keyword,*outputPrefix=NULL,*s;
  char *convention=NULL,*list=0;
  date_t start=NADate,final=NADate,orig=NADate;
  double startd,finald,origd,*ts;
  spectrum_t WaveList;
  double omega=NAN,q;//Q factor
  int atlas=0;
  harmonic_t harmonic;
  double **constants;
  int nvars=0;
  char **varnames=NULL;
  vector<string> filelist;
  
  struct timeval mainbefore;
  
  #if _OPENMP < 200505  // Version 2.5 May 2005
  #error Check whether your version _OPENMP (below) of openmp can do nesting properly. If so, update the above line accordingly.
  #pragma message "_OPENMP=" TOSTRING(_OPENMP)
  #endif

  cmd=fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"--help")==0) {
      print_help(argv[0]);
      exit(0);
      }
    else if(strcmp(keyword,"--control")==0) {
      controlfile= strdup(argv[n+1]);
      n++;
      n++;
      continue;
      }
    else if(strcmp(keyword,"--only-atlases")==0) {
      atlas=2;
      n++;
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 'a' :
          atlas=1;
          n++;
          break;

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

        case 'p' :
          outputPrefix= strdup(argv[n+1]);  /* directory */
          n++;
          n++;
          break;

        case 's' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          status=sscanf(s,"%d/%d/%d",&start.day,&start.month,&start.year);
          if(status==3)
            start.second=0.;
          free(s);
          break;

       case 'f' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          status=sscanf(s,"%d/%d/%d",&final.day,&final.month,&final.year);
          if(status==3)
            final.second=0.;
          free(s);
          break;

       case 'o' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          status=sscanf(s,"%d/%d/%d",&orig.day,&orig.month,&orig.year);
          if(status==3)
            orig.second=0.;
          free(s);
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

        case 'w' :
          WaveList.init(initialize_tide(),argv[n+1]);
          n++;
          n++;
          break;

        case 'q' :
          q=atof(argv[n+1]);
          n++;
          n++;
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
  
  ///It checks the list of files against the start and end dates with poc_timefilterlist()
  startd=cnes_time(start,'s');
  finald=cnes_time(final,'s');
  origd=cnes_time(orig,'s');
  printf("Selected data from %s to %s.\n",
    poctime_sdate_cnes(startd,'s'),poctime_sdate_cnes(finald,'s'));
  n=filelist.size();
  nt=poc_timefilterlist(&filelist,&startd,&finald,&ts,origd);
  printf("\n");
  printf("Kept %d/%d files with %d frames from %s to %s.\n\n",filelist.size(),n,nt,
    poctime_sdate_cnes(startd,'s'),poctime_sdate_cnes(finald,'s'));
  n=filelist.size();
  if(nt==0)
    TRAP_ERR_EXIT(1,"no frames within time boundaries\n");
  //start=poctime_getdatecnes(startd,'s');//not usefull
  
  if(WaveList.n>0){
    omega=WaveList.waves[0].omega;
    }

  /// It finally calls filter().
  constants=filter((const char**) varnames,nvars,filelist,gridfile,controlfile, nt,ts,omega,q,atlas,origd,outputPrefix);
  if(constants==NULL){
    STDOUT_BASE_LINE("filtering failed\n");
    exit(-1);
    }

  STDOUT_BASE_LINE(" ^^^^^^^^^^^^^ end of computation, that took %#g s with rt=%#g,rt2=%#g,rrt=%g,wt=%#g,filtering=%#g ^^^^^^^^^^^^^\n",difftime(mainbefore),rt,rt2,rrt,wt,ft);
  exit(0);
}
