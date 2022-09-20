
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
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

#include "polygons.hpp"
#include "exceptions.hpp"

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int plg_QuickFrameDeletion(vector<plg_t> & polygons, vector<plg_t> & outside, vector<plg_t> & temporary, frame_t frame, int mode)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//  
// 
//  
// cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
// {
//   int intersection, status;
//   int i,j,pivot;
//   double x,y;
//   bool ok;
//   
//   bool *keep=new bool[polygons.size()];
//   double center=frame.x_center();
//   
//   for(i=0;i<polygons.size();i++) {
//     intersection=0;
//     for(j=0;j<polygons[i].npt;j++) {
//       x=polygons[i].x[j];
//       y=polygons[i].y[j];
//       if(mode==PLG_SPHERICAL) x=geo_recale(x, center, 180.0);
//       ok=frame.inside(x, y);
//       if(ok) {
//         intersection++;
//         pivot=j;
//         break;
//         }
//       }
//     if(intersection!=0) {
//       temporary.push_back(polygons[i]);
//       keep[i]=false;
//       }
//     else {
//       keep[i]=true;
//       }
//     }
//     
//   for(i=0;i<polygons.size();i++) {
//     if(keep[i]) outside.push_back(polygons[i]);
//     }
//    
//   status=0;
//   return(status);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int plg_SplitPolygons(vector<plg_t> & polygons, const vector<plg_t> & selection, vector<plg_t> & outside, vector<plg_t> & inside, vector<plg_t> & intersecting, bool strict, int coordinates)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//  
//   split polygon set into inside, outside and intersecting sets
//  
// cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
// {
//   int i, j, flag, intersection, status;
//   int check=0;
//   int is_inside;
//   frame_t frame;
//   double x,y;
//   vector<plg_t>  temporary;
// 
// /*------------------------------------------------------------------------------
//   do first quick selection based on min/max frame                             */
//   frame=plg_spherical_minmax(selection);
//   status=plg_QuickFrameDeletion(polygons, outside, temporary, frame, coordinates);
//   
// /*------------------------------------------------------------------------------
//   then remove shorelines with no point inside selection polygons              */
//   double center=frame.x_center();
//   for(i=0;i<temporary.size();i++) {
//     intersection=0;
//     for(j=0;j<temporary[i].npt;j++) {
//       x=temporary[i].x[j];
//       y=temporary[i].y[j];
//       if(coordinates==PLG_SPHERICAL) x=geo_recale(x, center, 180.0);
//       is_inside=0;
//       flag=plg_single(selection[0], x, y, &is_inside, coordinates, check);
//       switch (strict) {
//         case true:
//           if(flag==PLG_POINT_INTERIOR) {
//             intersection++;
//             }
//           break;
//         case false:
//           if(flag!=PLG_POINT_EXTERIOR) {
//             intersection++;
//             }
//           break;
//         }
//       }
//     if(intersection==temporary[i].npt) {
//       inside.push_back(temporary[i]);
//       }
//     else if(intersection==0) {
//       outside.push_back(temporary[i]);
//       }
//     else {
//       intersecting.push_back(temporary[i]);
//       }
//     }
//   
//   return(status);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int plg_DeleteSelection(vector<plg_t> & polygons, const vector<plg_t> & selection, bool strict, int coordinates)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//  
//   wrapper
//  
// cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
// {
//   int i, j, flag, intersection, status;
//   int check=0;
//   int is_inside;
//   frame_t frame;
//   double x,y;
//   vector<plg_t>  temporary;
//   vector<plg_t> inside, intersecting,outside;
//   
//   status=plg_SplitPolygons(polygons, selection, outside, inside, intersecting, strict, coordinates);
// 
//   inside.clear();
//   intersecting.clear();
//   polygons.clear();
//   
//   polygons=outside;
//  
//   return(status);
// }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int defgrid(grid_t & grid, vector<plg_t> & polygons, double step)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  return a cartesian grid at the required resolution
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
  frame_t frame;
  projPJ projection;
  double ref_lat,ref_lon,x,y;
  char *pj_parameters=new char[512];

  frame= plg_recale(polygons);

  ref_lat=(frame.ymax+frame.ymin)/2;
  ref_lon=(frame.xmax+frame.xmin)/2;
  
/*------------------------------------------------------------------------------
  define stereographic oblique projection */
  projection=assign_StereoOblique(ref_lat,ref_lon,pj_parameters);
  
  geo_to_projection(projection, ref_lat, ref_lon, &x, &y);
    
  status=plg_cartesian(projection,polygons);
    
  frame= plg_cartesian_minmax(polygons);
    
  y=(frame.ymax+frame.ymin)/2;
  x=(frame.xmax+frame.xmin)/2;
    
  projection_to_geo(projection, &ref_lat, &ref_lon, x, y);

  printf("\ndefine background working grid: resolution=%lf m, projection=%s\n", step, pj_parameters);
  
  frame.dilatation(0.05);

  grid.xmin=frame.xmin;
  grid.xmax=frame.xmax;

  grid.ymin=frame.ymin;
  grid.ymax=frame.ymax;

  grid.dx=step;
  grid.dy=step;

  grid.xmin=NINT(grid.xmin/grid.dx)*grid.dx-grid.dx;
  grid.xmax=NINT(grid.xmax/grid.dx)*grid.dx+grid.dx;

  grid.ymin=NINT(grid.ymin/grid.dy)*grid.dy-grid.dy;
  grid.ymax=NINT(grid.ymax/grid.dy)*grid.dy+grid.dy;

  grid.nx=(int) NINT((grid.xmax-grid.xmin)/grid.dx)+1;
  grid.ny=(int) NINT((grid.ymax-grid.ymin)/grid.dy)+1;
  
  grid.modeH=0;

  status=map_completegridaxis(&grid,1);
  if(status!=0) return(status);

  grid.projection=projection;
  grid.proj4options=poc_strdup(pj_parameters);
  
  map_printgrid(grid);
  printf("\n");

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int FootPrint(vector<plg_t> & polygons, double step, string rootname, grid_t & grid, vector<size_t>* & pixels, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  pixelize
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
//   grid_t grid;
//   double step;
//   vector<size_t> *pixels;
  vector<size_t> active;
  int *passed;
  
  frame_t frame=plg_spherical_minmax(polygons);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  define the footprint master grid
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  status=defgrid(grid, polygons, step);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  detect polygons pixel
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  pixels=new vector<size_t> [polygons.size()];
  
  passed=new int[grid.Hsize()];
  
  for(size_t m=0; m<grid.Hsize(); m++) passed[m]=0;
    
  for(size_t p=0;p<polygons.size();p++) {
    for(size_t k=0;k<polygons[p].npt;k++) {
      double x=polygons[p].x[k];
      double y=polygons[p].y[k];
      int i,j;
      size_t m, pos;
      status=map_index(grid, x, y, &i, &j);
      m=j*grid.nx+i;
      passed[m]=1;
      pos=vpos(m,pixels[p]);
      if (pos==-1) {
        pixels[p].push_back(m);
        }
//       pos=vpos(m,active);
//       if (pos==-1) {
//         active.push_back(m);
//         }
      }
    if(debug) printf("%s : polygon %d #pixels=%d\n", __func__, p, pixels[p].size());
    }
  
  
  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int vcommon(const vector<T> & p, const vector<T> & q, vector<T> & r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  for(size_t k=0;k<p.size(); k++) {
    size_t pos;
    pos=vpos(p[k],q);
    if (pos!=-1) {
      r.push_back(p[k]);
      }
    }
  
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_CheckOverlapped(vector<plg_t> & polygons, bool repair, string rootname, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  detect and fix overlapping polygons
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
  int repaired=0;
  vector<size_t> *pixels;
  size_t count;
  double step;
  grid_t grid;
    
  step=5.e+03;
  status=FootPrint(polygons, step, rootname, grid, pixels, debug);

  count=0;
  for(size_t p=0;p<polygons.size()-1;p++) {
    frame_t frame=plg_spherical_minmax(&polygons[p], 1);
    for(size_t q=p+1;q<polygons.size();q++) {
      vector<size_t> common;
      status=vcommon(pixels[p],pixels[q],common);
      if(common.size()==0) continue;
      if(debug) printf("%s : polygons %d %d #pixels=%d\n", __func__, p, q, common.size());
      count++;
      for(int k=0;k<polygons[p].npt;k++) {
        size_t pos,posX,posY;
        int i,j;
        size_t m;
        status=map_index(grid, polygons[p].x[k], polygons[p].y[k], &i, &j);
        m=j*grid.nx+i;
        pos=vpos(m, common);
        if(pos==-1) continue;
        posX=vpos(polygons[p].x[k], polygons[q].x, polygons[q].npt);
        if(posX==-1) continue;
        posY=vpos(polygons[p].y[k], polygons[q].y, polygons[q].npt);
        if(posY!=posX) continue;
        printf("%s : polygons %d %d k %d pos %d %d\n", __func__, p, q, k, posX, posY);
        }
      }
    }
    
  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t> plg_extract_out_0(vector<plg_t> & shorelines_base, const vector<plg_t> & selection, string rootname, int mode, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
  reduce shorelines database by removing lines inside given selection polygons 
   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  vector<plg_t> extraction;
  size_t size;
  bool repair=true;
  
  vector<plg_t> outside, inside, intersecting, limits;
 
  size=shorelines_base.size();

  status=plg_CheckClosure(shorelines_base, true, 1, debug);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  make a first quick deletion of base polygons fully in exclusion zone

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
    
  bool strict=true;
  status=plg_SplitPolygons(shorelines_base, selection, outside, inside, intersecting, strict, mode);

  status=plg_CheckDuplicated(intersecting, true, rootname, false);
  
  size=intersecting.size();

  status=plg_CheckClosure(intersecting, true, 1, false);
  
  if(status!=0) return (extraction);
  
  extraction=outside;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  cut base polygons fully in exclusion zone

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
    
  printf("#cut polygons at frame limits\n");
 
  plg_t tmp;
  tmp.duplicate(selection[0]);
  
  limits.push_back(tmp);
 
  for(int k=0;k<intersecting.size();k++){
    vector<plg_t> pp;
    plg_t tmpp;
    tmpp.duplicate(intersecting[k]);
    pp.push_back(tmpp);
    status=plg_CutAtIntersections(pp, limits, PLG_SPHERICAL);    
/*------------------------------------------------------------------------------
    polygon was cut, removes outside sections*/
    status=plg_DeleteSelection(pp, selection, strict, mode);
    if(status!=0) return(extraction);
    for(int kk=0;kk<pp.size();kk++) {
      if(pp[kk].npt>1) extraction.push_back(pp[kk]);
      }
    pp.clear();
    }
    
  return (extraction);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t> plg_extract_out(const string & ShorelinesFile, const vector<plg_t> & selection, string rootname, int mode, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   wrapper
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  vector<plg_t> shorelines_base;
  vector<plg_t> extraction;
  size_t size;
  bool repair=true;

  int inFormat=plg_find_format(ShorelinesFile.c_str());
  
  status=plg_load(ShorelinesFile.c_str(), inFormat, shorelines_base);
  if(status!=0) return(extraction);
  
  extraction=plg_extract_out_0(shorelines_base, selection, rootname, mode, debug);
  
  return (extraction);
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
    "  merge a polygon set with an existing shoreline inside a polygon-defined extent \n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  Show this help and exit.\n"
    "  --debug  \n"
    "  --frame  followed by frame. WILL BE OVERWRITTEN TO [-10.0:0.0;42.5:50] !\n"
    "  -r  followed by the root name\n"
    "  -o  followed by the path of the output file. Default from the root name.\n"
    "  -O  followed by the output format\n"
    "  -p  followed by a path to a limiting polygone\n"
    "  -s  followed by the path to the shorelines. Default is /home/data/shorelines/gshhs-2.0/gshhs_f.cst\n"
    "  -z  followed by some zone name\n"
    );
  printf("\n"
    "FILE FORMATS :\n");
  plg_print_formats();
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  string keyword,PolygonFile="", FrameString="", rootname="", zone="", output;
  string outFormatName="";
  string ShorelinesFile[2];
  int inFormat,outFormat;
  vector<plg_t> extraction[2], merge;
  frame_t frame;
  vector<plg_t> limits;
  point2D_t point;
  bool debug=false, repair;
  int size[2];
  vector<int> list[2];
  double rmax=500.0;
  
  fct_echo( argc, argv);
  
  ShorelinesFile[0]="";
  ShorelinesFile[1]="";
  
  n=1;
  while (n < argc) {
    keyword=argv[n];
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

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
          ShorelinesFile[0]=argv[n+1];
          n++;
          n++;
          break;
        
        case 'i' :
          ShorelinesFile[1]=argv[n+1];
          n++;
          n++;
          break;
        
        case 'r' :
          rootname=argv[n+1];
          n++;
          n++;
          break;
        
        case 'z' :
          zone=argv[n+1];
          n++;
          n++;
          break;
        
        case 'O' :
          outFormatName= argv[n+1];
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
          else if(strncmp("--oF",keyword)==0){
            outFormatName= argv[n+1];
            n++;
            n++;
            }
          else if(strncmp("--frame",keyword)==0){
            FrameString= argv[n+1];
            n++;
            n++;
            }
          else if(strncmp("--debug",keyword)==0){
            debug=true;
            n++;
            }
          else if(strncmp("--rmax",keyword)==0){
            sscanf(argv[n+1],"%lf",&rmax);
            n++;
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
  
  if(outFormatName!="") {
    outFormat=plg_format_from_name(outFormatName);
    if(outFormat==PLG_FORMAT_UNKNOWN){
      plg_print_formats();
      exit(-1);
      }
    }
  else {
//     outFormat=PLG_FORMAT_SCAN;
    printf("\nWARNING : output format no specified, use %s format (can be specified with --oF XXX option)\n\n",ShorelinesFile[0].c_str());
    outFormat=plg_find_format(ShorelinesFile[0]);
    if(outFormat==PLG_FORMAT_UNKNOWN){
      plg_print_formats();
      exit(-1);
      }
    }
    
  if(rootname=="") {
    rootname="anonymous";
    printf("\nWARNING : rootname for auxilliary outputs files no specified, use %s (can be specified with -r XXX option)\n\n",rootname.c_str());
    }
  
  if(ShorelinesFile[0]=="") {
    ShorelinesFile[0]="/home/data/shorelines/gshhs-2.2/gshhs_f.cst";
    printf("\nWARNING : use %s default 1st shoreline file (as none was definded with -s option)\n\n",ShorelinesFile[0].c_str());
    }
  
  if(ShorelinesFile[1]=="") {
    ShorelinesFile[1]="/home/data/shorelines/TCHistolitt-2008/TCHistolitt-V1-0-RGF93-polygone.cst";
    printf("\nWARNING : use %s default 2nd shoreline file (as none was definded with -i option)\n\n",ShorelinesFile[1].c_str());
    }
  
  if(output=="") {
    output=rootname+"-merge.plg";
    printf("\nWARNING : rootname for main output files no specified, use %s (can be specified with -o XXX option)\n\n",output.c_str());
    }
    
  if(zone!="") {
    status=plg_SetZone(zone.c_str(), frame, rootname, point);
    limits.push_back(plg_t(frame, PLG_INIT_SEPARATE));
    status=plg_cartesian((projPJ) 0, limits);
    }

//   FrameString= "[-10.0:0.0;42.5:50]";
//   __ERR_BASE_LINE__(" *** OVERWROTE FRAME TO "+FrameString+" ***\n");
  if(FrameString!="") {
    status=plg_DecodeFrame(FrameString, frame);
    limits.push_back(plg_t(frame, PLG_INIT_SEPARATE));
    status=plg_cartesian((projPJ) 0, limits);
    }

  if(PolygonFile!="") {
    printf("#################################################################\n");
    printf("load selection polygons : %s\n", PolygonFile.c_str());
    inFormat=plg_find_format(PolygonFile);
    status=plg_load(PolygonFile, inFormat, limits);
    if(status!=0) TRAP_ERR_EXIT(status,"plg_load(\""+PolygonFile+"\",%d,,) error",inFormat);
    }
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  extract shorelines dataset in limits
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(limits.size()!=0) {
    printf("#################################################################\n");
    printf("extract shorelines\n");

    printf("-----------------------> treating %s (keep shorelines polygons #1 exterior to selection polygons)\n",ShorelinesFile[0].c_str());
    extraction[0]=plg_extract_out(ShorelinesFile[0], limits, rootname, PLG_SPHERICAL, debug);
    if(extraction[0].size()==0) {
      printf("no primary polygons selected, abort\n");
      goto error;
      }
  
    printf("-----------------------> treating %s (keep shorelines polygons #2 interior to selection polygons)\n",ShorelinesFile[1].c_str());
    extraction[1]=plg_extract(ShorelinesFile[1], limits, rootname, PLG_SPHERICAL, 2, debug);
    if(extraction[1].size()==0) {
      printf("no secondary polygons selected, abort\n");
      goto error;
      }
    status=plg_CheckOverlapped(extraction[1], true, rootname, debug);
    }
  else {
    status=plg_load(ShorelinesFile[0], inFormat, extraction[0]);
    status=plg_load(ShorelinesFile[1], inFormat, extraction[1]);
    }

  size[0]=extraction[0].size();
  size[1]=extraction[1].size();
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  detect polylines that will need to be assemblied
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(int k=0;k<extraction[0].size();k++) {
    if(plg_isclosed(extraction[0][k])==0) {
      if(debug) {
        printf("initial shorelines polygon %d not closed\n",k);
        printf("t=%lf p=%lf\n",extraction[0][k].t[0],extraction[0][k].p[0]);
        }
      list[0].push_back(k);
      }
    }
    
  for(int k=0;k<extraction[1].size();k++) {
    if(plg_isclosed(extraction[1][k])==0) {
      if(debug) {
        printf("initial shorelines polygon %d not closed\n",k);
        printf("t=%lf p=%lf\n",extraction[0][k].t[0],extraction[1][k].p[0]);
        }
      list[1].push_back(k);
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  concat extracted dataset
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=plg_merge(extraction[0],extraction[1]);
 
  if(debug) {
    string filename=rootname+"-extraction.plg";
    status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, extraction[0]);
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  join extremities of both set
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(int l=0;l<list[0].size();l++) {
    int s=list[0][l];
    plg_point_t point;
    paire_t paire;
    double d;
    point=plg_point_t(extraction[0][s],0);
    paire=plg_find_polygon(extraction[1], list[1], point, d);
    d*=1000.;
    if(d<500.0) {
      int r=paire.value[0]+size[0];
      status=plg_join(extraction[0][s],extraction[0][r],1.01*d);
      if(status==0) {
        printf("polygon %d and %d joined\n",s,r);
        }
      }
    point=plg_point_t(extraction[0][s],extraction[0][s].npt-1);
    paire=plg_find_polygon(extraction[1], list[1], point, d);
    d*=1000.;
    if(d<rmax) {
      int r=paire.value[0]+size[0];
      status=plg_join(extraction[0][s],extraction[0][r],1.01*d);
      if(status==0) {
/*------------------------------------------------------------------------------
        to be implemented in a light mode */
// 	status=plg_checkAutoSecant(&(extraction[0][s]),1,0);
        printf("polygon %d and %d joined\n",s,r);
        }
      }
    }

    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  check closure, perform concatenation when needed (repair=true)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("-----------------------> checking %s\n",output.c_str());
  repair=true;
  status=plg_CheckClosure(extraction[0], repair, 0, debug);
  
  for(int k=0;k<extraction[0].size();k++) {
//     printf("polygon %d : size=%d ",k, extraction[0][k].npt);
    if(plg_isclosed(extraction[0][k])==0) {
      printf("polygon %d : size=%d ",k, extraction[0][k].npt);
      printf("not closed -> t=%lf p=%lf",extraction[0][k].t[0],extraction[0][k].p[0]);
      printf("\n");
      }
//     printf("\n");
    }
    
  
  status=plg_save(output.c_str(), outFormat, extraction[0]);
  
  __OUT_BASE_LINE__(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(0);
error:
  __OUT_BASE_LINE__("\nerror detected, quit ... \n");
  exit(-1);
}
