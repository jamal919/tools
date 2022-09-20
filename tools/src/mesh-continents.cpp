
/*******************************************************************************

  T-UGO tools, 2006-2015

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
\brief Extends a 2D mesh on land.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main <a href=#func-members>function</a> for how this works
and the print_help <a href=#func-members>function</a> for how to use this.
*/
/*----------------------------------------------------------------------------*/


#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unistd.h> //access

#include "tools-structures.h"

#include "functions.h"
#include "geo.h"
#include "fe.h"
#include "archive.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "list.h"
#include "maths.h"
#include "map.def"
#include "mass.h"
#include "netcdf-proto.h"
#include "constants.h" //d2m

#include "zapper.h"     /*  rutin.h contains common utility routines  */

// #define MSH_MAXSIZE_OPEN     1
// #define MSH_MAXSIZE_SHELF    2
// #define MSH_TIDAL_WAVELENGTH 3
// #define MSH_TOPO_SLOPE       4
// #define MSH_SURFWAVE_LENGTH  5
// #define MSH_SURFWAVE_DZ      6
// #define MSH_SURFWAVE_CFL     7

#include "statistic.h"


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
    "  %s [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Extends a 2D mesh on land.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  Show this help and exit.\n"
    "  --rootname  root name of the MANY intermediate output files. Default: <BASE PATH OF THE MESH TO EXTEND>-continents\n"
    "  --debug  produce extra output files for debugging\n"
    "  --recover  use files produced with --debug to recover after a crash. Implies --debug\n"
    "  -m  followed by path to the mesh to extend\n"
    "  -o  followed by path to the output mesh. Default: <ROOTNAME>.nei\n"
    );
  /** \endcode */
}


#define STAR_LINE "********************************************************************************\n"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,j,status;
  bool debug=false,recover=false;
  char *keyword;
  string seaPath,landPath,meshPath,rootname,path;
  mesh_t sea,land,mesh;
  frame_t frame;
  grid_t grid;
  float *density,mask;
  vector<plg_t> islands;
  criteria_t criteria;
  range_t<double> edgeL;

  if(argc<=1){
    print_help(argv[0]);
    exit(0);
    }
  
  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=argv[n];
    
    switch (keyword[0]) {
      case '-':
        if(strcmp(keyword,"--debug")==0) {
          debug=true;
          n++;
          break;
          }
        if(strcmp(keyword,"--recover")==0) {
          recover=true;
          n++;
          break;
          }
        if(strcmp("--rootname",keyword)==0){
          rootname= argv[n+1];
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"--help")==0){
          print_help(argv[0]);
          exit(0);
          }
        switch (keyword[1]) {
        case 'h' :
          print_help(argv[0]);
          exit(0);

        case 'm' :
          seaPath= argv[n+1];
          n++;
          n++;
          break;

        case 'o' :
          meshPath= argv[n+1];
          n++;
          n++;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      }
    
    }
  
  if(recover) debug=true;
  
  printf("#################################################################\n");
  printf("reading "+seaPath+", ");fflush(stdout);
  
  status=fe_readmesh(seaPath.c_str(),MESH_FILE_FORMAT_TRIGRID,&sea);
  
  printf("compiling:\n");
  status=fe_list(&sea);//elements
  status=fe_edgetable(&sea,0,0);
  status=fe_codetable2(&sea,0,1,0,1);
  
  for(j=0;j<sea.nedges;j++){
    const edge_t *edge=&sea.edges[j];
    const vertex_t
      *v0=&sea.vertices[edge->extremity[0]],
      *v1=&sea.vertices[edge->extremity[1]];
    const double
      L=geo_km(v0->lon,v0->lat,v1->lon,v1->lat);
    
    edgeL << L;
    }
  
  printf("lengths in [%g;%g]km\n",edgeL.min,edgeL.max);
  
  printf("#################################################################\n");
  printf("computing islands\n");
  
  fe_limits2poly(sea,islands,0,false);
  n=islands.size();
  
  plg_degree_recale(&islands);
  frame=plg_spherical_minmax(islands);
  const double
    xDelta=frame.x_size();
  printf("longitude from %g to %g, delta=%g\n",frame.xmin,frame.xmax,xDelta);
  
  if(rootname==""){
    rootname=seaPath.substr(0,seaPath.length()-4)+"-continents";
    }

  if(debug) {
    path=rootname+"-islands.plg";
    STDOUT_BASE_LINE("saving to "+path+"\n");
    status=plg_save (path.c_str(), PLG_FORMAT_UNKNOWN, islands);
    }
  
  printf("#################################################################\n");
  printf("initialise default resolution criteria\n");
  
  const double m2km=0.001;
  criteria.shelf_minsize=edgeL.min;
  /* fail safe */
  updatemax(&criteria.shelf_minsize,xDelta*d2m*m2km/2e3);
  updatemax(&criteria.shelf_minsize,frame.y_size()*d2m*m2km/1e3);
  
  status=fe_InitCriteria("",criteria);
  
  criteria.cellsize=criteria.shelf_minsize;
  criteria.hmax=-1e4;
  criteria.maxrate=0.2;
  
  criteria.maxsize=edgeL.max*m2km;
  updatemax(&criteria.maxsize,100.);
  
  /* see set_density01() */
  criteria.factor1=criteria.maxsize/sqrt(-criteria.hmax);

  if(debug) {
    path=rootname+".crt";
    STDOUT_BASE_LINE("saving to "+path+"\n");
    status= fe_save_criteria(path.c_str(), criteria);
    }
  
  printf("#################################################################\n");
  printf("computing meshes of %d islands\n",n);
  
  if(meshPath=="") {
    meshPath=rootname+".nei";
    printf("no filename specified for output (-o [filename]); using default name ("+meshPath+")...\n");
    }
  
  mesh=sea;
  
/*------------------------------------------------------------------------------
  loop */
  const int nWidth=floor(log10(n))+1;
  
#undef TEST_MESH_CONTINENTS
//#define TEST_MESH_CONTINENTS 122
#ifdef TEST_MESH_CONTINENTS
  j=TEST_MESH_CONTINENTS;
  if(recover){
    STDERR_BASE_LINE("*** DISABLING --recover FOR TESTING PURPOSES ***\n");
    recover=false;
    }
  if(not debug){
    STDERR_BASE_LINE("*** IMPLYING --debug FOR TESTING PURPOSES ***\n");
    debug=true;
    }
#else
  j=0;
#endif
  
  for(;j<n;j++){
    
    if(j==0){
      bool skipped=xDelta<270;
      
      if(skipped) printf("Skipping");
      else printf("Including");
      printf(" polygone number %d.\n",j);
      
      if(skipped) continue;
      }
    
    bool recoverJ=false;
    char jStr[10];
    sprintf(jStr,"%0*d",nWidth,j);
    
    if(debug or recover){
      landPath=rootname+"-"+jStr+".nei";
      }
    
    if(recover){
      status=access(landPath.c_str(),R_OK);
      
      if(status==0)
        recoverJ=true;
      else
        STDOUT_BASE_LINE("access(\""+landPath+"\",R_OK) error (%d %s)\n",errno,strerror(errno));
      }
    
    if(recoverJ){
/*----------------------------------------------------------------------------*/
      printf(STAR_LINE "RECOVERY of "+landPath+"\n" STAR_LINE);fflush(stdout);
      
      status=fe_readmesh(landPath.c_str(),MESH_FILE_FORMAT_TRIGRID,&land);
      }
    else{
/*----------------------------------------------------------------------------*/
      printf(STAR_LINE "%s/%d:GENERATION\n" STAR_LINE,jStr,n);fflush(stdout);
      
//       if(j==12)
//         scanf("%d",&status);
      
      path=rootname+"-"+jStr+"-criteria.nc";
      status=fe_ComputeMeshsize(path.c_str(), criteria, &grid, &density, &mask, &islands[j], 1, 1, false);
      
      land=fe_nodit(criteria, grid, density, mask, &islands[j], 1, false, 0);
      
      if(debug) {
        STDOUT_BASE_LINE("saving to "+landPath+" (%d)\n",status);
        status= fe_savemesh(landPath.c_str(),MESH_FILE_FORMAT_TRIGRID, land);
        }
      }
    
    status=fe_geometry(&land);
#if 0
    path=rootname+"-"+jStr+".log";
    FILE *f=fopen(path.c_str(),"w");
    int i,k;
    fprintf(f,"%d vertices\n",land.nvtxs);
    for(i=0;i<land.nvtxs;i++){
      const vertex_t *vertex=&land.vertices[i];
      fprintf(f,"[%d] (%g;%g)%d",i,vertex->lon,vertex->lat,vertex->code);
      for(k=0;k<vertex->nngh;k++)fprintf(f," %d",vertex->ngh[k]);
      fprintf(f,"\n");
      }
    fprintf(f,"%d edges\n",land.nedges);
    for(i=0;i<land.nedges;i++){
      const edge_t *edge=&land.edges[i];
      fprintf(f,"[%d] %d<->%d %d|%d\n",i,edge->extremity[0],edge->extremity[1],edge->shared[0],edge->shared[1]);
      }
    fprintf(f,"%d limits\n",land.nlimits);
    for(i=0;i<land.nlimits;i++){
      const limit_t *limit=&land.limits[i];
      fprintf(f,"[%d] ",i);
      for(k=0;k<limit->nvertex;k++)fprintf(f," %d",limit->vertex[k]);
      fprintf(f,"\n");
      }
    fprintf(f,"%d triangles\n",land.ntriangles);
    for(i=0;i<land.ntriangles;i++){
      const triangle_t *triangle=&land.triangles[i];
      fprintf(f,"[%d] %d<->%d<->%d %d|%d|%d\n",i,
        triangle->vertex[0],triangle->vertex[1],triangle->vertex[2],
        triangle->edges[0],triangle->edges[1],triangle->edges[2]);
      }
    fclose(f);
#endif
    
#ifdef TEST_MESH_CONTINENTS
    TRAP_ERR_EXIT(ENOEXEC,"testing\n");
#endif
    
/*----------------------------------------------------------------------------*/
    printf(STAR_LINE "%s/%d:MERGE\n" STAR_LINE,jStr,n);fflush(stdout);
    
    mesh_t previous=mesh;
    
    mesh=fe_merge(land,mesh,edgeL.min/2);
    
    if(debug){
      landPath=rootname+"-"+jStr+"-merged.nei";
      STDOUT_BASE_LINE("saving to "+landPath+" (%d)\n",status);
      status=fe_savemesh(landPath.c_str(),MESH_FILE_FORMAT_TRIGRID, mesh);
      }
    
    fflush(stdout);
    previous.destroy();
    land.destroy();
    grid.free();
    delete[]density;
    }
  
  status=fe_savemesh(meshPath.c_str(),MESH_FILE_FORMAT_TRIGRID, mesh);
  
  return status;
}
