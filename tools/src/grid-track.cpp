
/*******************************************************************************

  T-UGO tools, 2006-2018

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Grid statistics about track parameter(s).

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
    "  %s [OPTIONS] \n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Grid statistics about track parameter(s).\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help : show this help and exit\n"
    "  -p : followed by the pattern of NetCDF input files with latitude, longitude and variables to make stats upon\n"
    "  -f : followed by (a suffix of) the first path to process. Default : whatever the pattern gives\n"
    "  -l : followed by (a suffix of) the last path to process. Default : whatever the pattern gives\n"
    "  -v : followed by the names of variables to make stats upon\n"
    "  -r : followed by the resoution in degrees. Default : 2\n"
    "  -o : followed by the path of the output\n"
    "  --debug-lon-lat : followed by lon and lat of a debug point. Default : none\n"
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
  
  const char *inputPath=0;  //path of the lat and lon of the control points
  const char *firstPath=0,*lastPath=0;
  const char *output=0;     //path of the output
  const char *keyword;      //option
  
  poc_deque_t<string> varNames;
  double resolution=2.;
  
  double debugLon=NAN,debugLat=NAN;
  
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
    
    if( strcmp(keyword,"--debug-lon-lat")==0 ){
      debugLon=atof(argv[n+1]);
      debugLat=atof(argv[n+2]);
      n++;
      n+=2;
      continue;
      }
    
    switch (keyword[0]) {
      case '-':
        
        switch (keyword[1]) {

        case 'p' :
          inputPath=argv[n+1];
          n++;
          n++;
          break;

        case 'f' :
          firstPath=argv[n+1];
          n++;
          n++;
          break;

        case 'l' :
          lastPath=argv[n+1];
          n++;
          n++;
          break;

        case 'v' :
          for(n++;n<argc && *argv[n]!='-';n++){
            varNames.push_back(argv[n]);
            }
          break;

        case 'r' :
          resolution=atof(argv[n+1]);
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
  
  const int vN=varNames.size();
  
  bool doHelp=false;
  
  if(inputPath==0) {
    STDOUT_BASE_LINE(" *** Please give the pattern of NetCDF input files with -v ***\n");
    doHelp=true;
    }
  
  if(vN==0) {
    STDOUT_BASE_LINE(" *** Please give the names of the variables to make stats upon with -v ***\n");
    doHelp=true;
    }
  
  if(output==0) {
    STDOUT_BASE_LINE(" *** Please give the path of the output with -o ***\n");
    doHelp=true;
    }
  
  if(doHelp==true){
    print_help(argv[0]);
    exit(-1);
    }
  
/*------------------------------------------------------------------------------
  pattern */
  glob_t globbed;
  globbed.gl_pathv=0;
  int fI;
  
  status=glob(inputPath,0,0,&globbed);
  
  if(status!=0 or globbed.gl_pathv==0){
    const char *msg;
    /* see man:/glob */
    switch(status){
    case GLOB_NOSPACE: msg="out of memory";break;
    case GLOB_ABORTED: msg="read error";break;
    case GLOB_NOMATCH: msg="no match";break;
    default: msg="unknown error code";
      }
    TRAP_ERR_EXIT(status,"glob(\"%s\",0,0,) error (%d %s)\n",inputPath,status,msg);
    }
  
/*------------------------------------------------------------------------------
  grid */
  grid_t grid;
  
  map_set2Dgrid(&grid,0.,-90.,+360.,+90.,resolution);
  n=grid.Hsize();
  
  int64_t debugM=-1;
  string debugPath;
  FILE *debugF=0;
  
  if(isfinite(debugLon) and isfinite(debugLat)){
    
    index_interpolation(grid,debugLon,debugLat,&debugM,0,NAN,0,0);
    
    if(0<=debugM and debugM<n){
      grid.xy((size_t)debugM,debugLon,debugLat);
      debugPath=output;
      replace(&debugPath,".nc","");
      asprintf(debugPath,"%+gE%+gN.dat",debugLon,debugLat);
      debugF=fopen(debugPath.c_str(),"w");
      }
    }
  
  double **sums,**sum2s,**mins,**maxs;
  int **counts;
  
  sums  =new double*[vN];
  sum2s =new double*[vN];
  mins  =new double*[vN];
  maxs  =new double*[vN];
  counts=new int   *[vN];
  
  for(vI=0;vI<vN;vI++){
    sums  [vI]=aset(n,0.);
    sum2s [vI]=aset(n,0.);
    mins  [vI]=aset(n,(double)+INFINITY);
    maxs  [vI]=aset(n,(double)-INFINITY);
    counts[vI]=aset(n,0);
    }
  
/*------------------------------------------------------------------------------
  process */
  
  int pn,pI;
  double *lon=0,*lat=0,*time=0;
  poc_data_t<double> data;
  
  for(fI=0;;fI++){
    inputPath=globbed.gl_pathv[fI];
    if(inputPath==0)
      break;
    
    if(firstPath!=0 and strrncmp(inputPath,firstPath)<0)
      continue;
    if(lastPath!=0 and strrncmp(lastPath,inputPath)<0)
      continue;
    
    printf("[%d/%u]",fI+1,globbed.gl_pathc);
    
    pn=poc_get_lon_lat_time(inputPath);
    
    /* (re-)allocate buffers */
    if(data.data==0 or data.length<pn){
      deletep(&lon);
      deletep(&lat);
      deletep(&data.data);
      
      data.length=1;
      while(data.length<pn)
        data.length*=2;
      
      lon=new double[data.length];
      lat=new double[data.length];
      data.data=new double[data.length];
      
      if(debugF!=0){
        deletep(&time);
        time=new double[data.length];
        }
      }
    
    printf("[%d/%u]",fI+1,globbed.gl_pathc);
    pn=poc_get_lon_lat_time(inputPath,lon,lat,time,0);
    
    for(vI=0;vI<vN;vI++){
      status=poc_inq_var(inputPath,varNames[vI],&data.info);
      NC_TRAP_ERROR(wexit,status,1,"poc_inq_var() error");
      
      status=poc_get_var(inputPath,data.info,data.data);
      NC_TRAP_ERROR(wexit,status,1,"poc_get_var() error");
      
      status=data.decode_mask();
      status=data.scale_data();
      
      int64_t m=-1;
      
      double *sum =sums  [vI];
      double *sum2=sum2s [vI];
      double *min =mins  [vI];
      double *max =maxs  [vI];
      int   *count=counts[vI];
      
      for(pI=0;pI<pn;pI++){
        const double datai=data.data[pI];
        if(datai==data.mask)
          continue;
        
        index_interpolation(grid,lon[pI],lat[pI],&m,0,NAN,0,0);
        if(m<0 or n<=m)
          continue;
        if(m==debugM)
          fprintf(debugF,"%.9g %.9g\n",time[pI],datai);
        
        count[m]++;
        sum  [m]+=datai;
        sum2 [m]+=square(datai);
        
        bool updated;
        updated=updatemin(&min[m],datai);
        updated=updatemax(&max[m],datai);
        }
      
      }
    
    }
  
/*------------------------------------------------------------------------------
  output */
  if(debugF!=0){
    STDOUT_BASE_LINE("closing "+debugPath+"\n");
    fclose(debugF);
    }
  
  map_completegridaxis(&grid,1);
  
  poc_var_t var;
  var=poc_save_grid(output,grid,__FILE__,__LINE__);
  status=poc_def_att(output,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
  status=poc_def_att(output,poc_att_t("history",cmd));
  
  for(vI=0;vI<vN;vI++){
    double *sum =sums  [vI];
    double *sum2=sum2s [vI];
    double *min =mins  [vI];
    double *max =maxs  [vI];
    int   *count=counts[vI];
    
    /* fix wrap around */
    int64_t m0=0,m1=grid.nx-1;
    for(pI=0;pI<grid.ny;pI++,m0+=grid.nx,m1+=grid.nx){
      int *countm0 =&count[m0];
      int *countm1 =&count[m1];
      *countm0+=*countm1;
      *countm1 =*countm0;
      
      double *summ0 =&sum [m0];
      double *summ1 =&sum [m1];
      *summ0+=*summ1;
      *summ1 =*summ0;
      
      double *sum2m0=&sum2[m0];
      double *sum2m1=&sum2[m1];
      *sum2m0+=*sum2m1;
      *sum2m1 =*sum2m0;
      
      double *minm0 =&min [m0];
      double *minm1 =&min [m1];
      updatemin(minm0,*minm1);
      *minm1 =*minm0;
      
      double *maxm0 =&max [m0];
      double *maxm1 =&max [m1];
      updatemax(maxm0,*maxm1);
      *maxm1 =*maxm0;
      }
    
    /* compute average and std */
    for(pI=0;pI<n;pI++){
      const int counti=count[pI];
      double *sumi =&sum [pI];
      double *sum2i=&sum2[pI];
      
      if(counti==0){
        *sumi =NC_FILL_DOUBLE;
        *sum2i=NC_FILL_DOUBLE;
        min[pI]=NC_FILL_DOUBLE;
        max[pI]=NC_FILL_DOUBLE;
        continue;
        }
      
      *sumi/=counti;
      *sum2i=sqrt(*sum2i/counti-square(*sumi));
      }
    
    var.init(varNames[vI]+"_count",NC_DOUBLE);
    status=poc_put_var(output,var,count);
    
    var.init(varNames[vI]+"_average",NC_DOUBLE);
    status=poc_put_var(output,var,sum);
    
    var.init(varNames[vI]+"_std",NC_DOUBLE);
    status=poc_put_var(output,var,sum2);
    
    var.init(varNames[vI]+"_min",NC_DOUBLE);
    status=poc_put_var(output,var,min);
    
    var.init(varNames[vI]+"_max",NC_DOUBLE);
    status=poc_put_var(output,var,max);
    }
}
