
/*******************************************************************************

  T-UGO tools, 2006-2018

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Interpolates TUGOm outputs at given points and times for altimetry.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <unistd.h> /* for access */
#include <stdio.h>
#include <glob.h>
#include <time.h> /* for time(), gmtime_r(), strftime() */

#include "poc-netcdf-data.hpp"
#include "functions.h" //safely includes omp.h
#include "poc-time.h"
#include "netcdf-proto.h" /* for poc_get_lon_lat_time() */

#include "fe-proto.h" /* for discretisation_init(), get_descriptor_address() */
#include "datastream.h"


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
    "  %s [-p lon_lat_list] -a atlas_convention [-s start -f end] -w wave1 [wave2 ... ] [OPTIONS] \n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Interpolates TUGOm outputs at given points and times for altimetry.\n"
    "  If start and end dates are provided, predicts tides.\n"
    "  It activates OpenMP parallelisation only if they are many points per time step.\n"
    "  If the output format is NetCDF, the 'units' attribute of the output variable is set to \"m\".\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help : show this help and exit\n"
    "  -p : followed by the pattern of NetCDF files with latitude, longitude and time variables\n"
    "  -c : followed by the convention of the time series.\n"
    "  -iV,--input-vars : followed by the input variables names for the elevation and the presure. Default: elevation pressure\n"
    "  -o : followed by the path of the output. Default: %%s-extraction.nc, with %%s the base name of the patterned input file.\n"
    "        If empty, add prediction to the patterned input file.\n"
    "  -oV,--output-vars : followed by the output variables names for the elevation and the presure. Default: DAC IB\n"
    );
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
  int n,status,vI;
  
  const char *inputConvention=0;
  const char *position_file=0;  //path of the lat and lon of the control points
  
  string output="%s-extraction.nc";     //path of the output
  const char *keyword;      //option
  
  const int nvars=2;
  string inputNames[]={"elevation","pressure"};
  string outputNames[]={"DAC","IB"};
  bool inputDefaults=1,outputDefaults=1;
  
  string cmd;
  cmd=fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=argv[n];
    
    if( strcmp(keyword,"-h")==0 or
        strcmp(keyword,"--help")==0 ){
      print_help(argv[0]);
      wexit(0);
      }
    
    if( strcmp(keyword,"-iV")==0 or
        strcmp(keyword,"--input-vars")==0 ) {
      if(inputDefaults)
        inputDefaults=false;
      else __FILE_LINE__(stdout,"multiple use of -iV or --input-vars : the last one overrides the previous.");
      
      for(vI=0;vI<nvars;vI++){
        inputNames[vI]=argv[n+1+vI];
        }
      
      n++;
      n+=nvars;
      continue;
      }

    if( strcmp(keyword,"-oV")==0 or
        strcmp(keyword,"--output-vars")==0 ) {
      if(outputDefaults)
        outputDefaults=false;
      else __FILE_LINE__(stdout,"multiple use of -oV or --output-vars : the last one overrides the previous.");
      
      for(vI=0;vI<nvars;vI++){
        outputNames[vI]=argv[n+1+vI];
        }
      
      n++;
      n+=nvars;
      continue;
      }
    
    switch (keyword[0]) {
      case '-':
        
        switch (keyword[1]) {

        case 'p' :
          position_file=argv[n+1];
          n++;
          n++;
          break;

        case 'c' :
          inputConvention=argv[n+1];
          n++;
          n++;
          break;

        case 'o' :
          output=argv[n+1];
          n++;
          n++;
          break;

        default:
          printf("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
        printf("unknown option %s\n",keyword);
        print_help(argv[0]);
        exit(-1);
      }
    
    }
  
  
/*------------------------------------------------------------------------------
  output variables */
  const float mask=1e4;
  poc_var_t bv[nvars];
  
  for(vI=0;vI<nvars;vI++){
    poc_var_t *bvi=&bv[vI];
    /* making sure _FillValue will appear as an attribute
        because of python netCDF4 bug... */
    bvi->init(outputNames[vI],NC_FLOAT,"","m",mask);
    *bvi << poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION);
    }
  
/*------------------------------------------------------------------------------
  input elevations */
  mesh_t mesh;
  int discretisation;
  discretisation_t *descriptor;
  grid_t grid;
  list_t list;
  
  filestream_t<float> *filestream=0;
  filestream=new filestream_t<float>("",inputConvention,inputNames[0].c_str(),0);
  
  sequence_t< UGfield_t<float> > streams[nvars];
  
/*------------------------------------------------------------------------------
  input tracks */
  glob_t globbed;
  globbed.gl_pathv=0;
  int fI;
  
  status=glob(position_file,0,0,&globbed);
  
  if(status!=0 or globbed.gl_pathv==0){
    const char *msg;
    /* see man:/glob */
    switch(status){
    case GLOB_NOSPACE: msg="out of memory";break;
    case GLOB_ABORTED: msg="read error";break;
    case GLOB_NOMATCH: msg="no match";break;
    default: msg="unknown error code";
      }
    TRAP_ERR_EXIT(status,"glob(\"%s\",0,0,) error (%d %s)\n",position_file,status,msg);
    }
  
  int l=0,dl=0;
  
  const int ncoord=3;
  double *coords[ncoord]={0,0,0},
    *&lon  =coords[0],
    *&lat  =coords[1],
    *&times=coords[2];
  float *buffers[nvars]={0,0};
  int pn,lastN=-1;
  double timei,loni,lati,r;
  int e,last;
  const int nframes=2;
  
  string outputI;
  
  for(fI=0;;fI++){
    position_file=globbed.gl_pathv[fI];
    if(position_file==0)
      break;
    
    do_cuu1(&l);
    printf("[%d/%u]",fI+1,globbed.gl_pathc);
    pn=poc_get_lon_lat_time(position_file); l+=2;
    
    /* (re-)allocate buffers */
    if(times==0 or lastN<pn){
      for(vI=0;vI<ncoord;vI++)
        deletep(&coords[vI]);
      
      for(vI=0;vI<nvars;vI++)
        deletep(&buffers[vI]);
      
      lastN=1;
      while(lastN<pn)
        lastN*=2;
      
      for(vI=0;vI<ncoord;vI++)
        coords[vI]=new double[lastN];
      
      for(vI=0;vI<nvars;vI++)
        buffers[vI]=new float[lastN];
      }
    
    do_cuu1(&l);
    printf("[%d/%u]",fI+1,globbed.gl_pathc);
    pn=poc_get_lon_lat_time(position_file,lon,lat,times,&bv[0].dimensions); l+=2;
    bv[1].dimensions=bv[0].dimensions;
    
    /* check time */
    timei=times[0];
    status=filestream->finalize(timei,-1);
    if(status!=0){
      timei=times[pn-1];
      status=filestream->finalize(timei,-1);
      }
    if(status!=0)
      continue;
    
    if(mesh.vertices==0){
      poc_var_t var;
      
      status=poc_inq_var(filestream->filename,filestream->varname,&var);
      status=poc_get_mesh(filestream->filename,var,&mesh,-1,0,&discretisation);
      
      fe_Allocate_and_CreateList(mesh,&grid,&list);
      
      discretisation_init(&mesh,discretisation,0);
      descriptor=get_descriptor_address(mesh,discretisation);
      
      streams[0]=initialise_UG_sequence(timei,filestream,nframes,&mesh,descriptor,stream_UGxNETCDF);
      dl=2;do_cuu1(&dl);
      streams[1]=initialise_UG_sequence(timei,"",inputConvention,inputNames[1],nframes,&mesh,descriptor,stream_UGxNETCDF,0);
      l+=2;
      }
    
    int pI,i;
    float z[nframes];
    
    dl=0;
    
    for(pI=0;pI<pn;pI++){
      
      /* update streams */
      timei=times[pI];
      
      if(timei > streams[0].frames[1].time){
        printf("%s",el);fflush(stdout);
        dl=2;
        }
      
      status=streams[0].check(timei,-1);
      if(status!=0){
        dl=0;
        for(vI=0;vI<nvars;vI++){
          buffers[vI][pI]=mask;
          }
        continue;
        }
      
      r=streams[1].get_r(timei);
      
      /* interpolate */
      loni =lon[pI];
      lati =lat[pI];
      
      e=fe_detect(mesh,mesh.triangles,list,grid,loni,lati,last);
      
      for(vI=0;vI<nvars;vI++){
        float *bufferii=&buffers[vI][pI];
        sequence_t< UGfield_t<float> > *streami=&streams[vI];
        
        for(i=0;i<nframes;i++)
          status=fe_interpolate2D(mesh,discretisation,streami->frames[i].x,mask,loni,lati,e,&z[i]);
        
        *bufferii=streami->interpolate_order(r,z,mask);
        
        /* inverse barometer is the additive inverse of the pressure */
        if(vI==1)
          (*bufferii)*=-1.;
        }
      
      }
    
    l+=dl;
    
    /* output */
    if(output!=""){
      string inputName;
      inputName=strrchr0(globbed.gl_pathv[fI],'/');
      inputName.erase(inputName.size()-3);
      
      asprintf(outputI="",output,inputName.c_str());
      }
    else{
      outputI=globbed.gl_pathv[fI];
      }
    
    status=access(outputI.c_str(),W_OK);
    
    if(status!=0){
      poc_global_t global;
      
      for(vI=0;vI<nvars;vI++)
        global << bv[vI];
      
      printf("[%d/%u:%d]Creating "+outputI+"\n",fI+1,globbed.gl_pathc,pn);fflush(stdout);l++;
      status=poc_create(outputI,global);
      }
    else{
      printf("[%d/%u:%d]Writing "+outputI+"\n",fI+1,globbed.gl_pathc,pn);fflush(stdout);l++;
      }
    
    status=poc_update_history(outputI,cmd);
    
    for(vI=0;vI<nvars;vI++)
      status=poc_put_var(outputI,bv[vI],buffers[vI]);
    
    }
  
}
