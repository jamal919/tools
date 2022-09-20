
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
\brief Extracts control-point values.

\date 2012-07-18 Damien Allain : creation based on comodo-detidor

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
#include "functions.h"
#include "matrix.h"
#include "poc-netcdf-data.hpp"
#include "poc-time.h"
#include "filter.h"
#include "map.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int extractControlPoints(const char **varnames, int nvars, vector<string> filelist,const char *gridfile,const char *controlfile, size_t nt,double *ts,double origd=NAN,const string output="control.nc")

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Harmonic analysis, eventually atlas, and detiding.
/**
\param **varnames list of variable names
\param nvars number of variables
\param filelist list of files
\param *controlfile path to list of control points
\param nt number of frames
\param ts times of the frames in CNES seconds
\param origd=NAN if specified, default time origin if not in the files
\returns NC_NOERR on success or the NetCDF error status on failure
*/
/*----------------------------------------------------------------------------*/
{
  double startd;
  double **signal=NULL; //< array[nsamples][nt] for the control samples
  int k,i,n;//variable, file and sample indexes
  int coordI;//coordinate index
  int status;
  int ifi,dfi;//input frame index,duplicate frame index
  int tfi;//total frame index
  poc_global_t global;
  int format=0;
  poc_var_t *info,detided_vtime;
  float *spec;
  size_t *nvalues,ntotal;
  int *n_spacedims;
  
/*-------------------------------------------------------------------------------------*/
  /// <h3> It checks the control file</h3>
  FILE *in;
  
  if(controlfile==0)
    in=stdin;
  else{
    in=fopen(controlfile,"r");//control input file
    if(in==NULL) TRAP_ERR_EXIT(errno,"Can not open %s (%s)",controlfile,strerror(errno));
    }

  if(nvars<=0)return NC_NOERR;
  if(nvars!=1){
    CHKERR_LINE_FILE(ENOEXEC,"no more than one variable accepted yet because of harmonic_t bug");
    exit(-2);/* Exec format error */
    }
  
  ///If there are duplicates in the frame times given, it returns NC_EEXIST.
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
    return NC_EEXIST;
    }
  minpos(ts,nt,&startd);
  
  /// <h3>It gets variable dimensions</h3>
  nvalues=new size_t[nvars];
  spec=new float[nvars];
  n_spacedims=new int[nvars];
  info=new poc_var_t[nvars];
  status=poc_inq(filelist[0],&global,0,&format);
  ntotal=0;
  for(k=0;k<nvars;k++){
    status=poc_get_var_length_and_mask(global,varnames[k],&info[k],&nvalues[k],&n_spacedims[k],&spec[k]);
    if(status==NC_ENOTATT)NC_CHKERR_BASE_LINE(status,"Use default mask value of %g for %s in "+filelist[0],spec[k],varnames[k]);/* warning only */
    else if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"poc_get_var_length_and_mask(,\"%s\",...) error in "+filelist[0]+"\n",varnames[k]);
    
    ntotal+=nvalues[k];
    }
  
  /// <h3> It reads the coordinates</h3>
  poc_grid_data_t coordinates;
  status=ENOENT;/* No such file or directory */
  if(gridfile!=NULL)
    status=poc_get_grid_data(gridfile,info[0],&coordinates);
  if(status!=NC_NOERR)
    status=poc_get_grid_data(filelist[0],info[0],&coordinates);
  if(status!=NC_NOERR){NC_CHKERR_BASE_LINE(status,"poc_get_grid_data() error on %s with both %s and "+filelist[0],varnames[0],gridfile);exit(-1);}
  
  int modeH=-9;
  if(nvalues[0]==coordinates.xv.length){
    for(k=0;k<nvars;k++)
      n_spacedims[k]=2;
    }
  else if(nvalues[0]==coordinates.xv.length*coordinates.yv.length){
    modeH=1;
    }
  else{
    TRAP_ERR_EXIT(ENOEXEC,"Not coded for lengths mismatch %s[%d] and lon[%d] x lat[%d] = [%d]\n",
      varnames[0],nvalues[0],
      coordinates.xv.length,coordinates.yv.length,
      coordinates.xv.length*coordinates.yv.length);
    }
  
  /// <h3> It initialises the control samples</h3>
  ///\bug 2011-12-15 Damien Allain : there are only samples for the first variable !!!
  int *controlIs=NULL;
  int nsamples=0;
  double **controlCoords;
  ///by reading the sample coordinates
  double minDist2,dist2;//SQUARES of distances
  status=fscanf(in,"%d\n",&nsamples);
  if(status==EOF) TRAP_ERR_EXIT(errno,"fscanf(,\"%%d\\n\",) error with %s : end of file.\n",controlfile);
  if(status!=1) TRAP_ERR_EXIT(errno,"fscanf(,\"%%d\\n\",) error %d with %s (%s)\n",errno,controlfile,strerror(errno));
  if(nsamples<=0){STDOUT_BASE_LINE("0 samples!\n");exit(-1);}
  controlIs=new int[nsamples];
  controlCoords=new double*[n_spacedims[0]];
  for(coordI=0;coordI<n_spacedims[0];coordI++){
    controlCoords[coordI]=new double[nsamples];
    }
  for(k=0;k<nsamples;k++) {//for each sample
    STDOUT_BASE_LINE("point %d/%d:",k,nsamples);fflush(stdout);
    for(coordI=0;coordI<n_spacedims[0];coordI++){
      fscanf(in,"%lf",&controlCoords[coordI][k]);
      printf(" %s~%g",coordinates[coordI].info.name.c_str(),controlCoords[coordI][k]);fflush(stdout);
      }
    ///and interpolating to the nearest neighbour
    controlIs[k]=-1;
    
    if(modeH==1){
      int Is[2];
      for(coordI=0;coordI<2;coordI++){
        minDist2=INFINITY;
        for(i=0;i<coordinates[coordI].length;i++){
          ///\bug 2011-12-15 Damien Allain : considering the coordinates are rectangular!
          dist2=pow(controlCoords[coordI][k]-coordinates[coordI].data[i],2);
          if(minDist2>dist2){
            minDist2=dist2;
            Is[coordI]=i;
            }
          }
        controlCoords[coordI][k]=coordinates[coordI].data[Is[coordI]];
        printf(" "+coordinates[coordI].info.name+"=%g",controlCoords[coordI][k]);
        }
      controlIs[k]=Is[0]+Is[1]*coordinates.xv.length;
      }
    else{
      minDist2=INFINITY;
      for(i=0;i<coordinates.xv.length;i++){
        dist2=0.;
        ///\note 2012-07-18 Damien Allain : even if it is masked
        //if(mask[i]==0)continue;
        for(coordI=0;coordI<2;coordI++){
          ///\bug 2011-12-15 Damien Allain : considering the coordinates are rectangular!
          dist2+=pow(controlCoords[coordI][k]-coordinates[coordI].data[i],2);
          }
        if(minDist2>dist2){
          minDist2=dist2;
          controlIs[k]=i;
          }
        }
      }
    
    int i2D=controlIs[k];
    if(n_spacedims[0]>=3){
      controlIs[k]+=controlCoords[2][k]*coordinates.xv.length;
      printf(" => index=%d/%d(%g)",i2D,controlIs[k],sqrt(minDist2));fflush(stdout);
      }
    else
      printf(" => index=%d(%g)",controlIs[k],sqrt(minDist2));fflush(stdout);
    if(controlIs[k]<0 || nvalues[0]<=controlIs[k]){printf(" ***ERROR : %d out[0;%d[***\n",controlIs[k],nvalues[0]);break;}
    for(coordI=0;coordI<n_spacedims[0];coordI++){
      if(modeH==1)break;
      if(coordI==2)
        i2D=controlIs[k];
      controlCoords[coordI][k]=coordinates[coordI].data[i2D];
      printf(" "+coordinates[coordI].info.name+"=%g",controlCoords[coordI][k]);
      }
    printf("\n");
    }//EO for each sample
  fclose(in);
  for(coordI=0;coordI<coordinates.vc;coordI++){
    coordinates[coordI].destroy_data();
    }
  
  /// <h3>It allocates memory for the samples and buffers</h3>
  if(nsamples>0){
    signal=new double*[nsamples];
    //continuous memory buffer for easier writing
    signal[0]=new double[nsamples*nt];
    for(n=1;n<nsamples;n++) {
      signal[n]=&signal[n-1][nt];
      }
    }
  
/* As per testing of 12-13 Oct 2009, reading frames is slower
  BUT FOR CUMBERSOMELY WRITTEN NetCDF 4 FILES !!! */
#define read_frame 0
  
#if read_frame==1
  float **buffer;
  buffer=new float*[nvars];
  for(k=0;k<nvars;k++){
    exitIfNull(buffer[k]=new float[nvalues[k]]);
    }
#endif

/*-------------------------------------------------------------------------------------*/
  /// <h2>For each <i>input</i> files</h2>
  struct timeval b4;
  gettimeofday(&b4);
  double od0;           //< date origine of first file
  double *stime;        //< times since od0
  stime=new double[nt];
  
  for(i=0; i<filelist.size(); i++){
    poc_global_t global;
    int nif;//number of input frames
    date_t origine;
    double *time,od;
    poc_var_t vtime;
    #if NETCDF_CAN_PARALLEL_IO == 0
    #pragma omp critical(threadUnsafeNetcdf)
    #endif
    {
    status = poc_gettime(filelist[i].c_str(), &origine, &time, &nif, &vtime);
    }//END OF #pragma omp critical(threadUnsafeNetcdf)
    od=cnes_time(origine,'s');
    if(isnan(od))od=origd;
    if(i==0){
      detided_vtime=vtime;
      od0=od;
      }
    
#if read_frame==0
    tfi=pos(time[0]+od,ts,nt);//total frame index
    if(tfi<0)
      tfi=0;
/*------------------------------------------------------------------------------
    read slice */
    #if NETCDF_CAN_PARALLEL_IO == 0
    #pragma omp critical(threadUnsafeNetcdf)
    #endif
    {
    int ncid;
    status=nc_open(filelist[i].c_str(),NC_NOWRITE,&ncid);
    for(k=0;k<nvars;k++)
    for(n=0;n<nsamples;n++) {
      poc_var_t *infok=&info[k];
      printf("\r(%d/%d;%04.3g)reading, in "+filelist[i]+", "+infok->name+"[",i+1,filelist.size(),difftime(b4));fflush(stdout);
      status=nc_inq_varid(ncid,varnames[k],&infok->id);
      
      const int dimC=infok->dimensions.size();
      const poc_dim_t *dim;
      size_t *start,*count;
      int
        nvaluesk=nvalues[k],
        controlIsn=controlIs[n];
      
      start=new size_t[dimC];
      count=new size_t[dimC];
      
      for(int i=0;i<dimC;i++){
        dim=&infok->dimensions[i];
        
        if(isT(*dim)){
          start[i]=pos(ts[0]-od,time,nif);
          count[i]=min(nif-start[i],nt);
          }
        else{
          nvaluesk/=dim->len;
          
          start[i]=controlIsn/nvaluesk;
          
          count[i]=1;
          
          controlIsn%=nvaluesk;
          }
        printf("[%d:%d]",start[i],count[i]);fflush(stdout);
        }
      
      printf("] (sample %d/%d)...",n+1,nsamples);fflush(stdout);
      status=poc_get_vara(ncid,infok->id,start,count,&signal[n][tfi]);
      
      delete[]start;
      delete[]count;
      }
    
    status=nc_close(ncid);
    }//END OF #pragma omp critical(threadUnsafeNetcdf)
#endif
    
    for(ifi=0;ifi<nif;ifi++) {
/*------------------------------------------------------------------------------
      read frame */
      tfi=pos(time[ifi]+od,ts,nt);//total frame index
      if(tfi<0)
        continue;
#if read_frame==1
      #if NETCDF_CAN_PARALLEL_IO == 0
      #pragma omp critical(threadUnsafeNetcdf)
      #endif
      {
      printf("\r(%d/%d;%04.3g)reading "+filelist[i]+" (frame %d/%d)...",i+1,filelist.size(),difftime(b4),tfi+1,nt);fflush(stdout);
      ///It reads the frame in the input file
      for(k=0;k<nvars;k++) {
        status=poc_get_vara(filelist[i],info[k],ifi,buffer[k]);
        if(status){NC_CHKERR_BASE_LINE(status,"poc_get_vara error with "+info[k].name+" in "+filelist[i]);exit(1);}
        }
      }//END OF #pragma omp critical(threadUnsafeNetcdf)
      for(n=0;n<nsamples;n++) {
        signal[n][tfi]=buffer[0][controlIs[n]];
        }
#endif
      stime[tfi]=ts[tfi]-od0;
      }
    
    delete[]time;
    }
  printf("\r(%d/%d;%04.3g)read complete.%s\n",i,filelist.size(),difftime(b4),el);
/*-------------------------------------------------------------------------------------*/

  //buffer clean-up
  delete[]controlIs;
#if read_frame==1
  for(k=0;k<nvars;k++){
    delete[]buffer[k];
    }
  delete[]buffer;
#endif
  
/*------------------------------------------------------------------------------
  compose output header */
  poc_dim_t samples_dim("samples",nsamples);
  global=poc_global_t("initialised around " __LINE_FILE_PACKAGE_REVISION);
  string coordVarNames;
  /* time variable */
  detided_vtime.dimensions[0].init(nt);
  detided_vtime<<poc_att_t("units","seconds since "+ (string)poctime_sdate_cnes(od0) );
  coordVarNames=detided_vtime.name;
  global.variables<<detided_vtime;
  detided_vtime.id=global.variables.back().id;
  /* coordinates */
  for(coordI=0;coordI<n_spacedims[0];coordI++){
    coordinates[coordI].info.dimensions=poc_list_t<poc_dim_t>(samples_dim);
    coordVarNames+=" "+coordinates[coordI].info.name;
    global.variables<<coordinates[coordI].info;
    coordinates[coordI].info.id=global.variables.back().id;
    }
  /* samples */
  info[0].dimensions=poc_list_t<poc_dim_t>(samples_dim);
  info[0].dimensions<<detided_vtime.dimensions[0];
  global<<poc_att_t("coordinate_variable_names",coordVarNames);
  global.variables<<info[0];
  info[0].id=global.variables.back().id;
  
/*------------------------------------------------------------------------------
  write output */
  STDOUT_BASE_LINE("Writing samples to "+output+"... ");fflush(stdout);
  gettimeofday(&b4);
  
  if(strrncmp(output,".nc")!=0){
    FILE *f;
    char *sdate;
    
    f=fopen(output.c_str(),"w");
    
    for(tfi=0;tfi<nt;tfi++){
      
      sdate=poctime_sdate_cnes(ts[tfi],'s','_');
      fprintf(f,"%s",sdate);
      delete[]sdate;
      
      for(k=0;k<nsamples;k++){
        fprintf(f," %g",signal[k][tfi]);
        }
      
      fprintf(f,"\n");
      }
    
    fclose(f);
    goto cleanUp;
    }
  
  switch(format){
  case NC_FORMAT_CLASSIC:
    format=NC_CLASSIC_MODEL;
    break;
  case NC_FORMAT_64BIT:
    format=NC_64BIT_OFFSET;
    break;
  case NC_FORMAT_NETCDF4:
    format=NC_NETCDF4;
    break;
  case NC_FORMAT_NETCDF4_CLASSIC:
    format=NC_NETCDF4|NC_CLASSIC_MODEL;
    break;
  default:
    TRAP_ERR_EXIT(ENOEXEC,"not coded for format=%d\n",format);
    }
  
  status=poc_create(output,global,0,format);
  NC_CHKERR_BASE_LINE(status,"poc_create(\""+output+"\",,,%d) error ",format);
  
  status=poc_put_var(output,detided_vtime,stime);
  NC_CHKERR_BASE_LINE(status,"poc_put_var(\""+output+"\",\""+detided_vtime.name+"\",) error");
  
  status=poc_put_var(output,info[0],signal[0]);
  NC_CHKERR_BASE_LINE(status,"poc_put_var(\""+output+"\",\""+info[0].name+"\",) error");
  
  for(coordI=0;coordI<n_spacedims[0];coordI++){
    status=poc_put_var(output,coordinates[coordI].info,controlCoords[coordI]);
    NC_CHKERR_BASE_LINE(status,"poc_put_var(\""+output+"\",\""+coordinates[coordI].info.name+"\",) error");
    }
  
cleanUp:
  deletep2D(&controlCoords,n_spacedims[0]);
  deletep2D(&signal,1);
  
  printf("took %gs\n",difftime(b4));
  
  delete[]info;
  delete[]nvalues;
  delete[]n_spacedims;
  
  return status;
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
    "  %s file1 [ file2 ... ] [OPTIONS] ...\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Extracts control-point values.\n"
    "\n"
    "OPTIONS :\n"
    "  -v   followed by the list of variables to control\n"
    "  -c   followed by the path of the list of control points: it is an ascii file with the number of control points followed by their coordinates (longitude latitude [layer]).\n"
    "  -l   followed by the path of the list of files to analyse. This list will override the list given as arguments.\n"
    "  -g   followed by the path of the grid\n"
    "  -s : followed by the start date. See DATE FORMATS below\n"
    "  -f : followed by the end date. See DATE FORMATS below\n"
    "  -o : followed by the default date origin. See DATE FORMATS below\n"
    "  -O : followed by the path of the output file.\n"
    "      If the extension is ``.nc``, the format will be NetCDF; otherwise, ASCII.\n"
    "      Default: control.nc\n"
    "  --test : followed by axes, then indexes. MUST BE LAST OPTION IF USED.\n"
    );
  print_poctime_scan_date_help(0);
  /** \endcode */
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
  int k;//file index, variable index
  int n,status;//argument index or number of files, NetCDF status
  const char *gridfile=NULL,*controlfile=NULL,*keyword;
  string outputFile="control.nc";
  const char *list=0;
  date_t start=NADate,final=NADate,orig=NADate;
  double startd,finald,origd,*ts;
  size_t tn;
  int nvars=0;
  char **varnames=NULL;
  vector<string> filelist;
  bool doHelp=false;

  string cmd;
  cmd=fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=argv[n];
    
    if(strcmp(keyword,"--help")==0) {
      print_help(argv[0]);
      exit(0);
      }
    
    if(strcmp(keyword,"--test")==0) {
      
      int status;
      poc_data_t<double> extract;
      
      status=poc_inq_var(filelist[0],varnames[0],&extract.info);
      if(gridfile==0)
        gridfile=filelist[0].c_str();
      status=poc_get_axes(gridfile,&extract.info,0xf);
      STDOUT_BASE_LINE("variable axes:"+extract.info.axes+".\n");
      
      n++;
      if(n>=argc) TRAP_ERR_EXIT(-1,"only showing axes of variable\n");
      
      const string axes=argv[n];
      STDERR_BASE_LINE("indexes axes:"+axes+".\n");
      
      poc_deque_t<int> indexes;
      int i;
      n++;
      for(i=0;n<argc && i<axes.length();n++,i++)
        indexes<<atoi(argv[n]);
      
      FILE *ascii=0;
      int ncid=0;
      STDOUT_BASE_LINE("Ouput to "+outputFile+" :\n");
      if(strrncmp(outputFile,".nc")==0){
        status = nc_open(outputFile.c_str(), NC_WRITE, &ncid);
        if(status==ENOENT){
          status = nc_create(outputFile.c_str(), NC_NOCLOBBER, &ncid);
          NC_TRAP_ERROR(wexit,status,1,"nc_create(\""+outputFile+"\",NC_NOCLOBBER,) error");
          status = nc_enddef(ncid);
          NC_TRAP_ERROR(wexit,status,1,"nc_enddef((\""+outputFile+"\")) error");
          }
        else
          NC_TRAP_ERROR(wexit,status,1,"nc_open(\""+outputFile+"\",NC_WRITE,) error");
        }
      else{
        ascii=fopen(outputFile.c_str(),"a");
        if(ascii==0) TRAP_ERR_EXIT(errno,"fopen(\""+outputFile+"\",\"a\") error (%d %s)\n",errno,strerror(errno));
        write_header_if_empty(ascii,"from numpy import array");
        fprintf(ascii,"# %s\n",cmd.c_str());
        fprintf(ascii,"# ");BASE_LINE(ascii,"\n");
        }
      
      for(int j=0;j<nvars;j++){
        
        string extractAxes;
        poc_var_t outputVar;
        int k=0,m=0;
        int col=0;
        
        for(i=0;i<filelist.size();i++){
          status=poc_inq_var(filelist[i],varnames[j],&extract.info);
          
          if(i==0){
            status=poc_get_axes(gridfile,&extract.info,0);
            extractAxes=extract.info.axes;
            if(ascii==0){
              outputVar=extract.info;
              for(k=outputVar.axes.size()-1;k>=0;k--){
                char c=outputVar.axes[k];
                if(strchr(axes.c_str(),c)==0)
                  continue;
                outputVar.axes.erase(k);
                outputVar.dimensions.erase(k);
                }
              outputVar<<poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION);
              status=poc_def_var(ncid,outputVar);
              }
            }
          else
            extract.info.axes=extractAxes;
          
          status=extract.init(axes);
          
          struct timeval before;
          STDOUT_BASE_LINE("reading %d from "+filelist[i]+" : ",extract.length);fflush(stdout);
          gettimeofday(&before);
          status=extract.read_data(filelist[i],indexes);
          printf("done in %g s.\n",difftime(before));
          
          if(ascii==0)
           for(k=0;k<extract.length;k++,m+=1){
            status=poc_put_vara(ncid,outputVar,m,&extract.data[k]);
            NC_TRAP_ERROR(wexit,status,1,"poc_put_vara((\""+outputFile+"\"),(\""+outputVar.name+"\"),%d,) error",m);
            }
          else
           for(k=0;k<extract.length;k++,m+=1){
            
            if(m==0)
              col+=fprintf(ascii,string(varnames[j])+"=array([");
            else
              col+=fprintf(ascii,",");
            
            if(col>=110){
              fprintf(ascii,"\n");
              col=0;
              }
            
            /* hourly time in h since 1950-01-01 needs 8 digits */
            col+=fprintf(ascii,"%.9g",extract.data[k]);
            }
          
          extract.destroy_data();
          }
        
        if(ascii!=0)
          fprintf(ascii,"])\n");
        }
      
      if(ascii==0)
        status=nc_close(ncid);
      else
        fclose(ascii);
      STDOUT_BASE_LINE(" "+outputFile+" completed.\n");
      
      TRAP_ERR_EXIT(ENOEXEC,"testing\n");
      }
    
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 'c' :
          controlfile= argv[n+1];
          n++;
          n++;
          break;

        case 'l' :
          list= argv[n+1];
          n++;
          n++;
          break;

        case 'g' :
          gridfile= argv[n+1];
          n++;
          n++;
          break;

        case 's' :
          status=poctime_scan_date(argv[n+1],&start,0,-1);
          if(status!=0){
            STDOUT_BASE_LINE(" *** Could not scan start date \"%s\" ***\n",argv[n+1]);
            doHelp=true;
            }
          n++;
          n++;
          break;

        case 'f' :
          status=poctime_scan_date(argv[n+1],0,&final,-1);
          if(status!=0){
            STDOUT_BASE_LINE(" *** Could not scan final date \"%s\" ***\n",argv[n+1]);
            doHelp=true;
            }
          n++;
          n++;
          break;

        case 'o' :
          status=poctime_scan_date(argv[n+1],&orig,0,-1);
          if(status!=0){
            STDOUT_BASE_LINE(" *** Could not scan date origin \"%s\" ***\n",argv[n+1]);
            doHelp=true;
            }
          n++;
          n++;
          break;

        case 'O' :
          outputFile= argv[n+1];
          n++;
          n++;
          break;

        case 'v' :
          for(k=1;n+k<argc && argv[n+k][0]!='-';k++);
          nvars=k-1;
          if(varnames){
            STDOUT_BASE_LINE("multiple use of -v : the last one overrides the previous.");
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
    
    }
  
  printf("#number of variables : %d\n",nvars);
  if(nvars==0) {
    STDOUT_BASE_LINE(" *** Please give the names of the variables to control ***\n");
    doHelp=true;
    }
  for(k=0;k<nvars;k++) {
    printf("#variable %d: %s\n",k,varnames[k]);
    }
  
  if(list!=0) {
    filelist=load_filelist(list);
    }
  
  if(filelist.size()==0) {
    STDOUT_BASE_LINE(" *** Please give a NetCDF file to analyse ***\n");
    }
  
  if(doHelp==true){
    print_help(argv[0]);
    exit(-1);
    }
  
/*------------------------------------------------------------------------------
  It checks the list of files against the start and end dates with poc_timefilterlist() */
  startd=cnes_time(start,'s');
  finald=cnes_time(final,'s');
  origd=cnes_time(orig,'s');
  printf("Selected data from %s to %s.\n",  poctime_sdate_cnes(startd,'s'),poctime_sdate_cnes(finald,'s'));
  
  n=filelist.size();
  tn=poc_timefilterlist(&filelist,&startd,&finald,&ts,origd);
  printf("\n");
  printf("Kept %d/%d files with %d frames from %s to %s.\n\n",filelist.size(),n,tn, poctime_sdate_cnes(startd,'s'),poctime_sdate_cnes(finald,'s'));
  
  n=filelist.size();
  if(tn==0)
    TRAP_ERR_EXIT(1,"no frames within time boundaries");

/*------------------------------------------------------------------------------
  It finally calls extractControlPoints() */
  status=extractControlPoints((const char**) varnames,nvars,filelist,gridfile,controlfile, tn,ts,origd,outputFile);
  if(status!=NC_NOERR)NC_CHKERR_BASE_LINE(status,"extractControlPoints() error");

  return status;
}
