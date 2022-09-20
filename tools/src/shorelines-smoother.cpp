
/**************************************************************************

  T-UGO tools, 2006-2013

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Convert between different coastline formats.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main <a href=#func-members>function</a> for how this works
and the print_help <a href=#func-members>function</a> for how to use this.
*/
/*----------------------------------------------------------------------------*/


#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <string.h>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <map>

#include "tools-structures.h"

#include "fe.h"
#include "fe-proto.h"
#include "geo.h"
#include "archive.h"
#include "netcdf-proto.h"
#include "poc-time.h"
#include "functions.h"
#include "datastream.h"
#include "filter.h"
#include "map.h"

#include "polygons.hpp"
#include "exceptions.hpp"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  SGfield_t<float> distance_to_shorelines(string filename, vector<plg_t> & polygons, double dx, int mode, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  float mask=-2;
  float mean_bkp,mean;
  plg_array_t shorelines;
  SGfield_t<float> field;
  float *r=0;
  int    proxies[9][2]={1,0,1,1,0,1,-1,1,-1,0,-1,-1,0,-1,1,-1,0,0};
  size_t count;
  vector<plg_t> limits;
    
  field.grid=new grid_t;
  grid_t & grid=*field.grid;
  
  debug=true;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  define cartesian working grid and compute landmask
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  build working grid (cartesian, structured), apply projection to shorelines */
  status=fe_defgrid(&grid, polygons, limits, dx);
  
/*------------------------------------------------------------------------------
  set interior/exterior flag */
  printf("#------------------------------------------------------\n");
  printf("step 1: set landmask flag\n");
  float *landmask=set_landmask01(grid, polygons);
  
/*------------------------------------------------------------------------------
  in case where polygons are limits, not shorelines */
  if(mode!=0) {
    printf("boundary mode: use inverted landmask\n");
    for (size_t m=0; m<grid.Hsize();m++) landmask[m]=-landmask[m];
    }
  else {
    printf("shoreline mode: use native landmask\n");
    }
  
  if(filename!="") remove(filename.c_str());

  if(filename!="") {
    status=save_SGXY(filename.c_str(), grid, landmask, mask, "landmask-0", "NONE", "landmask-0", 1, 1);
    }
    
  float *flag=0;

  float *distance = new float[grid.Hsize()];
  float *delta    = new float[grid.Hsize()];

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  initialize "0" distance (i.e. nodes around shorelines)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(int j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      size_t m=j*grid.nx+i;
      distance[m]=mask;
      double dmin=+INFINITY;
      for(int k=0; k<9; k++) {
        int ii=proxies[k][0];
        int jj=proxies[k][1];
        size_t mm=jj*grid.nx+ii;
        if(landmask[m]==landmask[mm]) continue;
        double dx=grid.dx*(i-ii);
        double dy=grid.dy*(j-jj);
        double d=dx*dx+dy*dy;
        if(d<dmin) dmin=d/2;
        }
      if(dmin!=+INFINITY) distance[m]=dmin;
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  compute square distance 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  do {
#pragma omp parallel for
    for(int j=0;j<grid.ny;j++) {
      for(int i=0;i<grid.nx;i++) {
        size_t m=j*grid.nx+i;
        if(distance[m]!=mask) continue;
        double dmin=+INFINITY;
        for(int k=0; k<9; k++) {
          int ii=proxies[k][0];
          int jj=proxies[k][1];
          size_t mm=jj*grid.nx+ii;
          if(distance[mm]==mask) continue;
          double dx=grid.dx*(i-ii);
          double dy=grid.dy*(j-jj);
          double d=distance[mm]+dx*dx+dy*dy;
          if(d<dmin) dmin=d;
          }
        if(dmin!=+INFINITY) distance[m]=dmin;
        }
      }
      count=occurence(mask, distance, grid.Hsize());
    } while (count!=0);
  
  for(int m=0;m<grid.Hsize();m++) {
    if(distance[m]!=mask) distance[m]=landmask[m]*sqrt(distance[m]);
    }
   
  if(filename!="") {
    status=save_SGXY(filename.c_str(), grid, distance, mask, "distance", "NONE", "distance", 0, 0);
    }
  
  delete[] landmask;
  if(flag!=0) delete[] flag;
  if(r!=0)    delete[] r;
  delete[] delta;

  field.x=distance;
  field.mask=mask;
  
  delete[] shorelines.p;
  
  return(field);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(char *prog_name)

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
    "  Convert between different coastline formats.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  Show this help and exit.\n"
    "  -i  followed by the path of the input file.\n"
    "  -o  followed by the path of the output file. Default from the input file and output format if both are given.\n"
    "  -I  followed by the input format. Default from the extension of the input file.\n"
    "  -O  followed by the output format. Default from the extension of the output file.\n"
    "  -d  followed by the decimation factor. Default : 1 (no decimation).\n"
    "  --proj=...  projection parameters. Default : no projection.\n");
  printf("\n"
    "FILE FORMATS :\n");
  plg_print_formats();
  /** \endcode */
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int format(string filename, string formatname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int fmt;
  
  if(formatname!="") {
    fmt=plg_format_from_name(formatname);
    }
  else if(filename!="") {
    fmt=plg_find_format(filename);
    
    }
  else fmt=PLG_FORMAT_SCAN;
  
  if(fmt==PLG_FORMAT_UNKNOWN){
    plg_print_formats();
    exit(-1);
    }

  return(fmt);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l;
  int n,status;
  string keyword,input="",PolygonFile="", FrameString="", rootname="", output="", zone="";
  statistic_t stat;
  string inFormatName,outFormatName;
  string ShorelinesFile="";
  int inFormat,outFormat,decimation=1;
  plg_t *polygones;
  int npolygones;
  char *proj4=NULL;
  char *s;
  int extract=0;
  double ocean_t,ocean_p;
  vector<plg_t> extraction, limits, shorelines, smoothed;
  frame_t frame;
  point2D_t point;
  bool debug=false;
  double isocontour=-0.4;
  bool tuned=true;
  bool enforcement=true, show_map=false, subdivide=false;
  double lengthscale=2000;
  double dx=-1;
  int nitems;
  double factor;
  bool show_mesh;
  plg_filter_t plg_filter;
  
  
  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=argv[n];
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 'l' :
          nitems=sscanf(argv[n+1],"%lf",&lengthscale);
          if(nitems!=1) {
            printf("invalid option argument"+keyword+" "+argv[n+1]+"\n");
            print_help(argv[0]);
            exit(-1);
            }
          n++;
          n++;
          break;
	  
        case 'r' :
          nitems=sscanf(argv[n+1],"%lf",&dx);
          if(nitems!=1) {
            printf("invalid option argument"+keyword+" "+argv[n+1]+"\n");
            print_help(argv[0]);
            exit(-1);
            }
          n++;
          n++;
          break;
	  
        case 'p' :
          PolygonFile=argv[n+1];
          n++;
          n++;
          break;
        
        case 'o' :
          output=argv[n+1];
          n++;
          n++;
          break;
        
        case 's' :
          ShorelinesFile=argv[n+1];
          n++;
          n++;
          break;
        
        case 'z' :
          zone=argv[n+1];
          n++;
          n++;
          break;
        
        case 'I' :
          inFormatName= argv[n+1];
          n++;
          n++;
          break;
        
        case 'O' :
          outFormatName= argv[n+1];
          n++;
          n++;
          break;
        
        case 'd' :
          decimation=strtol(argv[n+1],&s,10);
          if(s==argv[n+1]){
            printf("unknown decimation factor %s\n",argv[n+1]);
            print_help(argv[0]);
            exit(-1);
            }
          n++;
          n++;
          break;
        
        case 'h' :
          print_help(argv[0]);
          exit(0);

        case '-' :
          if(keyword=="--help"){
            print_help(argv[0]);
            exit(0);
            }
          else if(strncmp("--proj=",keyword)==0){
            proj4=strdup(argv[n]+7);
//             proj4=strdup(keyword.c_str());
            n++;
            }
          else if(strncmp("--extract",keyword)==0){
            extract=1;
            n++;
            }
          else if(strncmp("--frame",keyword)==0){
            FrameString= argv[n+1];
            n++;
            n++;
            }
          else if("--isocontour"==keyword){
            sscanf(argv[n+1],"%f",&isocontour);
            n++;
            n++;
            }
          else if("--no-enforcement"==keyword){
            enforcement=false;
            n++;
            }
          else if(strncmp("--debug",keyword)==0){
            debug=true;
            n++;
            }
          else {
            printf("unknown option "+keyword+"\n");
            print_help(argv[0]);
            exit(-1);
	    }
          break;

        default:
          printf("unknown option "+keyword+"\n");
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
        printf("unknown option "+keyword+"\n");
        print_help(argv[0]);
        exit(-1);
      }
    }
    
// string ShorelinesFile="/home/data/shorelines/data/coasts-2011-resaved.cst";
  if(ShorelinesFile=="") ShorelinesFile="/home/data/shorelines/gshhs-2.0/gshhs_f.cst";
  
  if(dx==-1) dx=lengthscale/4.0;
  
//  zone="korea";  // debugging

  if(zone!="") {
    status=plg_SetZone(zone.c_str(), frame, rootname, point);
    limits.push_back(plg_t(frame, PLG_INIT_SEPARATE));
    status=plg_cartesian((projPJ) 0, limits);
    extract=1;
    }
//  extract=0;  // debugging

  if(rootname=="") rootname="anonymous";
//   FrameString="[-10.0:-1.0;42.5:50]";
//   FrameString="[350.0:360.0;42.5:50]";
  
  if(FrameString!="") {
    status=plg_DecodeFrame(FrameString, frame);
    limits.push_back(plg_t(frame, PLG_INIT_SEPARATE));
    status=plg_cartesian((projPJ) 0, limits);
    extract=1;
    }
    
  if(PolygonFile!="") {
    inFormat=plg_find_format(PolygonFile);
    status=plg_load(PolygonFile, inFormat, limits);
    if(status!=NC_NOERR) __NC_TRAP_ERROR__(wexit,status,1,"plg_load(\""+input+"\",%d,,) error",inFormat);
    extract=1;
    }

  if(output=="") output=rootname+"-shorelines-smoothed.plg";
  if(input=="") input=rootname+"-extraction.plg";
    
  frame_t tmp=frame;
  tmp.dilatation(-0.025);

  vector<float>  isovalues;
  SGfield_t<float> resolution;

  if(extract==1) {
/*------------------------------------------------------------------------------
    extract shorelines dataset*/
    printf("#################################################################\n");
    printf("extract shorelines\n");

    vector<plg_t> dilated=plg_dilatation_cartesian(limits, (double) 1.05, debug);
    extraction=plg_extract(ShorelinesFile, dilated, rootname, PLG_SPHERICAL, 0, debug);
    if(extraction.size()==0) goto error;
    printf("#################################################################\n");
    printf("store extracted shorelines in %s\n",input.c_str());
    status=plg_save(input.c_str(), PLG_FORMAT_SCAN, extraction);
    inFormatName="SCAN";
    }
  else {
    input=ShorelinesFile;
    }
    
  printf("#################################################################\n");
  printf("smooth shorelines\n");

/*------------------------------------------------------------------------------
  load new shorelines data set*/
  inFormat=format(input, inFormatName);
  status=plg_load(input, inFormat, shorelines);
  if(status!=0) goto error;
  
/*------------------------------------------------------------------------------
  define isovalues for smoothed contours*/
  isovalues.clear();
  for(int k=0;k<7;k++) isovalues.push_back(0.2-k*0.1);
  show_map=debug;
  
  factor=2.0;
  show_mesh=false;
  
  plg_filter.init(lengthscale, resolution, factor, dx, subdivide, tuned, enforcement, show_map, show_mesh);
  
  smoothed=plg_smooth(rootname, shorelines, plg_filter, isovalues, isocontour, point, 0, debug);
    
  if(shorelines.size()==0) goto error;
  
  for(int k=0;k<2;k++) {
    plg_smoothlines(smoothed);
    }
/*------------------------------------------------------------------------------
  restore original extraction frame*/
  plg_destroy_entries(extraction);
  extraction=plg_extract(smoothed, limits, rootname, PLG_SPHERICAL, 0, debug);

#if 0    
    status=plg_cartesian(projection, contours);
    vector<plg_t> sampled;
    for(int s=0;s<contours.size();s++) {
      contours[s].flag=new char[contours[s].npt-1];
      for(int kk=0;kk<contours[s].npt-1;kk++) contours[s].flag[kk]='T';
      plg_t test=plg_resample_rigid_obsolete(contours[s], projection, (double) 5.0, 1, 0);
      sampled.push_back(test);
      }
    output=rootname+"-resampled-" + tmp.str() +".plg";
    status=plg_save(output.c_str(), PLG_FORMAT_SCAN, sampled);
    sampled.clear();
#endif

  printf("#################################################################\n");
  printf("store smooth shorelines in %s\n",output.c_str());
  outFormat=format(output, outFormatName);
  status=plg_save(output.c_str(), outFormat, extraction);    

  
  __OUT_BASE_LINE__(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(0);
  
error:
  __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}
