
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
\brief Set landmask and interpolate topography for a given grid.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main <a href=#func-members>function</a> for how this works
and the print_help <a href=#func-members>function</a> for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "config.h"

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <unistd.h> // for getpid
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "functions.h"
#include "constants.h"
#include "map.h"
#include "geo.h"
#include "polygones.h"
#include "grd.h"
#include "fe-proto.h"
#include "sym-io.h"
#include "ascii.h"
#include "netcdf-proto.h"
#include "poc-netcdf-data.hpp"
#include "poc-time.h"
#include "datastream.h"


extern int map_loadfield_grd(const char *input,grid_t *grid, float **buffer, float *mask);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_loadfield_cdf(const char *input,grid_t *grid, float **buffer, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    v,status;
  int    verbose,frame=0;
  cdfgbl_t data_info,grid_info;
  cdfvar_t h_info;
  
  *buffer=0;

  verbose=0;
  
  status=cdf_globalinfo(input,&data_info,verbose);
  if(status!=0) return(-1);
  
  if(verbose) {
    for (v=0;v<data_info.nvarsp;v++) {
      printf("variable %3d: name %s, type %d,ndim %d \n",v,(data_info.variable[v]).name,data_info.variable[v].type,data_info.variable[v].ndim);
      }
    }
  status=cdf_varinfo(input,"bathymetry",&h_info,verbose);
  if(status!=0) return(-1);

  status=cdf_globalinfo(input,&grid_info,verbose);
  if(status!=0) return(-1);
  status= poc_getgrid2d (input, grid_info, h_info, grid);
  if(status!=0) return(-1);

  *buffer=new float[grid->nx*grid->ny];
  status= poc_getvar2d (input, h_info.id, frame, *buffer, mask, h_info);
  if(status!=0) return(-1);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_loadgrid(const char *input,const char *variable, grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    v,status;
  int    verbose;
  cdfgbl_t data_info,grid_info;
  cdfvar_t h_info;

  verbose=0;
  status=cdf_globalinfo(input,&data_info,verbose);
  for (v=0;v<data_info.nvarsp;v++) {
    printf("variable %3d: name %s, type %d,ndim %d \n",v,(data_info.variable[v]).name,data_info.variable[v].type,data_info.variable[v].ndim);
    }
  status=cdf_varinfo(input,variable,&h_info,1);
  if(status!=0) return(-1);

  status=cdf_globalinfo(input,&grid_info,verbose);
  status= poc_getgrid2d (input, grid_info, h_info, grid);
  if(status!=0) return(-1);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_loadfield(const char *filename,grid_t *grid, float **buffer, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;

  if(!strrncmp(filename,".grd")) {
    status=map_loadfield_grd(filename,grid, buffer, mask);
    }
  else if(!strrncmp(filename,".nc")) {
    status=map_loadfield_cdf(filename,grid, buffer, mask);
    if(status!=0) return(-1);
//     for(m=0;m<grid->ny*grid->nx;m++) if((*buffer)[m]!=*mask) (*buffer)[m]*=-1.;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_interpolate(const vector<string> & databases, grid_t grid, float * & topo, float mask, int average)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,n, status;
  grid_t topo_grid;
  float  *topo_base,topo_mask, z;
  double x,y;
  
/* *-----------------------------------------------------------------------------
  read and interpolate topo database */
  topo=new float[grid.Hsize()];
  
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      topo[k]=mask;
      }
    }

  for(n=0;n<databases.size();n++) {
/* *------------------------------------------------------------------------------
    Load bottom topography*/
//     printf("#################################################################\n");
    printf("#----------------------------------------------------------------\n");
    printf("load bathymetry file: %s\n",databases[n].c_str());
    status=map_loadfield(databases[n].c_str(),&topo_grid,&topo_base,&topo_mask);
    if(status !=0) goto error;
/* *------------------------------------------------------------------------------
    Interpolate and store topo on symphonie grid*/
    for(j=0;j<grid.ny;j++) {
      for(i=0;i<grid.nx;i++) {
        k=j*grid.nx+i;
        x=grid.x[k];
        y=grid.y[k];
        switch (average) {
          case 0:
            status=map_interpolation(topo_grid, topo_base,topo_mask,x,y,&z);
            break;
          case 1:
            status=map_average(topo_grid, topo_base,topo_mask,grid,i,j,&z);
            break;
          }
        if (z!=topo_mask) {
          topo[k]=z;
          }
        }
      }
    topo_grid.free();
    delete[] topo_base;
    }
  return(0);
  
error:
  printf("file loading error, return with bad status\n");
  return(-1);

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

  int inFormat=plg_find_format(ShorelinesFile.c_str());
  
  status=plg_load(ShorelinesFile.c_str(), inFormat, shorelines_base);
  if(status!=0) return(extraction);
  
  extraction=plg_extract_out(shorelines_base, selection, rootname, mode, debug);
  
  return (extraction);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int get_polygon(const string & rootname, const string & PolygonsFile, vector<plg_t> & polygons, const plg_t & GridLimits,
                  grid_t & sgrid, SGfield_t<float> & fresolution, double wsize, double factor, double isocontour, bool filter, point2D_t point,
                  bool reduce, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  load a dedicated polygon file, possibly reduce and filer processed
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;

  printf("#################################################################\n");
  printf("load landmask polygons file: %s\n",PolygonsFile.c_str());fflush(stdout);
  status=plg_load(PolygonsFile, polygons);
  if(status) TRAP_ERR_EXIT(status,"plg_load(\"%s\",,) error\n",PolygonsFile.c_str());

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  reduce polygon set to grid extent to minimize computational costs at next use
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(reduce) {
    frame_t GridFrame=plg_spherical_minmax(&GridLimits,1);
    double dilatation=0.05;
    dilatation=1./(1.* min(sgrid.nx,sgrid.ny) );
    
    GridFrame.dilatation(dilatation);

    plg_t ExtendedLimits=plg_t(GridFrame, PLG_INIT_SEPARATE);
    
    int target=0;

    vector<plg_t> dilated;
/*------------------------------------------------------------------------------
  grid limits issue when covered area is large*/
//     dilated.push_back(ExtendedLimits);
    dilated.push_back(GridLimits);
    dilated=plg_dilatation_cartesian(dilated,1.05,debug);
    
    vector<plg_t> extraction=plg_extract(polygons, dilated, rootname, PLG_SPHERICAL, target, debug);
    plg_destroy_entries(polygons);
    polygons=extraction;
    
    string TmpPoly=(string) rootname +"-boundaries-reduced.plg";
    printf("#################################################################\n");
    printf("save reduced boundary polygons[%u] : %s\n",polygons.size(),TmpPoly.c_str());fflush(stdout);
    status=plg_save(TmpPoly.c_str(), PLG_FORMAT_SCAN, polygons);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  filter polygons (removing unwanted details and merging micro blocks)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(filter) {
    double lengthscale=0.0;

    vector<float> isovalues;
    for(int k=0;k<7;k++) isovalues.push_back(0.2-k*0.1);
        
//     SGfield_t<float> fresolution;
//     fresolution.x=sresolution;
//     fresolution.grid=&sgrid;
//     fresolution.mask=-1;
    
    int target=2;
    bool enforcement=false, tuned=true, show_map=false, show_mesh=false, subdivide=false;
    if(debug) show_map=true;

    plg_filter_t plg_filter;
    plg_filter.init(lengthscale, fresolution, factor, wsize, subdivide, tuned, enforcement, show_map, show_mesh);
    
    vector<plg_t> smoothed=plg_smooth(rootname, polygons, plg_filter, isovalues, isocontour, point, target, debug);
    
    plg_destroy_entries(polygons);
    polygons=smoothed;
    
    printf("#################################################################\n");
    STDOUT_BASE_LINE("remove swimming pools\n");
    status=plg_RemoveInlandLimits(polygons, point.x, point.y, (projPJ) 0, PLG_POINT_EXTERIOR, debug);
    string TmpPoly=(string) rootname +"-boundaries-filtered.plg";
    printf("#################################################################\n");
    printf("save filtered boundary polygons : %s\n",TmpPoly.c_str());
    status=plg_save(TmpPoly.c_str(), PLG_FORMAT_SCAN, polygons);
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int create_polygon(const string & rootname, string & ShorelinesFile, vector<plg_t> & polygons, const plg_t & GridLimits,
                  grid_t & sgrid, SGfield_t<float> & fresolution, double wsize, double factor, double isocontour, bool filter, point2D_t point,
                  bool removePools, bool show_map, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  create a polygon from shoreline dataset
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  string TmpPoly;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create selection polygon for shorelines extraction/filtering

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  frame_t GridFrame=plg_spherical_minmax(&GridLimits,1);
  double dilatation=0.05;

  dilatation=1./(1.* min(sgrid.nx,sgrid.ny) );
  GridFrame.dilatation(dilatation);

/*------------------------------------------------------------------------------
  grid limits issue when covered area is large*/
//   plg_t ExtendedLimits=plg_t(GridFrame, PLG_SPHERICAL);
  plg_t ExtendedLimits;
  ExtendedLimits.duplicate(GridLimits);

  TmpPoly=(string) rootname +"-frame.plg";
  status=plg_save(TmpPoly.c_str(), PLG_FORMAT_SCAN, &ExtendedLimits, 1);
  
  if(ShorelinesFile=="") {
    ShorelinesFile="/home/data/shorelines/gshhs-2.2/gshhs_f.cst";
    printf("shorelines filename not given, use default : %s\n",ShorelinesFile.c_str());
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create model limits polygons

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  plg_filter_t plg_filter;
  bool tuned=false, enforcement=false, show_mesh=false, subdivide=false;
  if(debug) show_map=true;
  double lengthscale=0.0;
  plg_filter.init(lengthscale, fresolution, factor, wsize, subdivide, tuned, enforcement, show_map, show_mesh);
  
  status=plg_CreateLimits(ShorelinesFile, ExtendedLimits, rootname, plg_filter, point, isocontour, filter, removePools, debug);
  if(status) TRAP_ERR_EXIT(status,"plg_CreateLimits() error\n");
  
  ExtendedLimits.destroy();
  
  TmpPoly=(string) rootname +"-boundaries.plg";
  status=plg_load(TmpPoly.c_str(), PLG_FORMAT_SCAN, polygons);
  if(status) TRAP_ERR_EXIT(status,"plg_load(\""+TmpPoly+"\",,) error\n");
      
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int process_landmask(const string & rootname, string & exclusion_polygon, vector<plg_t> & polygons,
                       const grid_t & z_grid, const grid_t & f_grid, bool checks, float* & landmask, bool fix_isolated, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  //int target=2;
  
  printf("#################################################################\n");
  printf("compute landmask from model limits polygons \n");
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  process exclusion polygons if any
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(exclusion_polygon!="") {
    vector<plg_t> selection;
    status=plg_load(exclusion_polygon, PLG_FORMAT_SCAN, selection);
    vector<plg_t> extraction=plg_extract_out(polygons, selection, rootname, PLG_SPHERICAL, debug);
    plg_destroy_entries(polygons);
    polygons=extraction;
    string TmpPoly=(string) rootname +"-with-exclusion.plg";
//     printf("#################################################################\n");
//     printf("save reduced boundary polygons : %s\n",TmpPoly.c_str());
    status=plg_save(TmpPoly.c_str(), PLG_FORMAT_SCAN, polygons);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  compute landmask from model limits polygons
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  landmask=set_landmask02(z_grid,polygons,checks);
  plg_destroy_entries(polygons);
  
  if(fix_isolated==false)
    return 0;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  fix pinched/isolated cells
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  size_t count[2];
  
  count[0]=occurence(1.0f, landmask, z_grid.Hsize());
  
  float *topo=0;
  float topomask=NC_FILL_FLOAT;
  mesh_t mesh;
  int FixPinched=1;
  int verbose=1;
  char *flag=new char[z_grid.Hsize()];
  for(size_t m=0;m<z_grid.Hsize();m++) flag[m]=(landmask[m]==1);
  
  status=fe_ZGrid2Quadrangle(z_grid, f_grid, flag, topo, topomask, mesh, FixPinched, verbose);
  
  for(size_t m=0;m<z_grid.Hsize();m++) landmask[m]=-1;

  for(size_t m=0;m<mesh.nquadrangles;m++) {
    double t,p;
    int i,j;
    status=fe_position(mesh, mesh.quadrangles[m], &t, &p,0);
    status=map_index(z_grid,t,p,&i,&j);
    size_t n=j*z_grid.nx+i;
    landmask[n]=1;
    }
  
  count[1]=occurence(1.0f, landmask, z_grid.Hsize());
  
  printf("landmask checked for isolated/pinched cells:  initial=%d, final=%d \n", count[0], count[1]);
  
  printf("landmask sucessfully completed\n");
      
  return(status);
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
    "  %s [OPTIONS] notebookFile\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Set landmask and interpolate topography for a given grid.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  Show this help and exit.\n"
    "  -inv  switch off altitude to depth conversion\n"
    "  -average/-pointwise  switches between averaging/pointwise interpolations of the bathymetry. Default: pointwise.\n"
    "  -zmin  topography values will be constrained above this minimum. Default:-99999.\n"
    "  -format  fprintf format for the symphonie ASCII file, preferably with leading blank. Default: \" %%6.1f\".\n"
    "  -r  followed by the root name. Default : symphonie\n"
    "  -b  followed by the path to a bathymetry. Several bathymetries are allowed.\n"
    "  -g  followed by the path to the grid file. Disables -n.\n"
    "  -variable_t  followed by the grid variable name at T point\n"
    "  -n  followed by the path to the notebook file\n"
    "  --no_mask : disable mask creation. Disables -p\n"
    "  -p  followed by a path to a land mask polygons filename. Disables -s. If not given, compute one from the shorelines dataset (see below).\n"
    "  -s : followed by the path to a shorelines filename. Default : /home/data/shorelines/gshhs-2.2/gshhs_f.cst\n"
    "  --point  followed by the coordinates of a reference as \"[<lon>;<lat>]\"\n"
    "  --increment : followed by increment for grid frontier. Default : 10\n"
    "  --debug\n"
    "  --keep-pools : keep CREATED pools\n"
    "  --isocontour : followed by an isocontour value for the smoothed shoreline : -1 for land, 1 for ocean. Default : -0.4\n"
    );
  printf("\n"
    "FILE FORMATS :\n");
  plg_print_formats();
  /** \endcode */
}


extern ostringstream assertionCmd;
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  geo_t projection;
  float  mask=1.e+10f;
  int i,j,k,n,status;
  FILE *out;
  char *keyword;
  int flag=-1;
  char *rootname=NULL,*output=NULL,*input=NULL;
  char *notebook=NULL,*gridfile=NULL,*format=NULL;
  char *varname_t=NULL, *varname_f=NULL;
  string exclusion_polygon, PolygonsFile="";
  vector<string> bathymetry;
  grid_t cgrid,sgrid;

  float  *landmask,*topo;
  cdfgbl_t data_info,grid_info;
  cdfvar_t h_info;

  int inversion=1,average=0,increment=10;
  float zmin=-99999;
  double *dx, *dy, factor=1.000, wsize=0;
  
  int skip=1;
  
  bool do_landmask=true, do_bathymetry=true;
  char *vmask=0, *vtopo=0;
  poc_data_t<float>  meta_H_z, meta_H_f;

  string ShorelinesFile, TmpPoly, PointString;
  point2D_t point;
  double isocontour=-0.4;
  bool checks=false, filter=false, reduce=true, removePools=true;
  bool show_map=false;
  bool debug=false;
  
  poc_var_t z_gridVar, f_gridVar;
  grid_t z_grid,f_grid;
  
  fct_echo( argc, argv);
  
  if(argc<2) {
    print_help(argv[0]);
    exit(0);
    }
  
  n=1;
  while (n < argc) {
    keyword=argv[n];
    switch (keyword[0]) {
      case '-':
        if(strcmp(keyword,"--help")==0) {
          print_help(argv[0]);
          exit(0);
          }
        if(strcmp(keyword,"-inv")==0) {
          inversion=0;
          n++;
          continue;
          }
        if(strcmp(keyword,"--exclude")==0) {
          exclusion_polygon=argv[n+1];
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"--checks")==0) {
          checks=true;
          n++;
          continue;
          }
        if(strcmp(keyword,"--debug")==0) {
          debug=true;
          n++;
          continue;
          }
        if(strcmp(keyword,"--smooth")==0) {
          filter=true;
          n++;
          continue;
          }
        if(strcmp(keyword,"--keep-pools")==0) {
          removePools=false;
          n++;
          continue;
          }
        if(strcmp(keyword,"--increment")==0) {
          sscanf(argv[n+1],"%d",&increment);
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"-average")==0) {
          average=1;
          n++;
          continue;
          }
        if(strcmp(keyword,"-pointwise")==0) {
          average=0;
          n++;
          continue;
          }
        if(strcmp(keyword,"--show_map")==0){
          show_map=true;
          n++;
          continue;
          }
        if(strcmp(keyword,"-dx")==0) {
          sscanf(argv[n+1],"%lf",&wsize);
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"-zmin")==0) {
          sscanf(argv[n+1],"%f",&zmin);
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"-format")==0) {
          format=strdup(argv[n+1]);
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"-vmask")==0) {
          vmask=strdup(argv[n+1]);
          do_landmask=false;
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"--no_mask")==0) {
          do_landmask=false;
          n++;
          continue;
          }
        if(strcmp(keyword,"-vtopo")==0) {
          vtopo=strdup(argv[n+1]);
          do_bathymetry=false;
          n++;
          n++;
          continue;
          }
        else if(strncmp("--point",keyword)==0){
          PointString= argv[n+1];
          n++;
          n++;
          continue;
          }
        else if(strncmp("--isocontour",keyword)==0){
          sscanf(argv[n+1],"%lf",&isocontour);
          n++;
          n++;
          continue;
          }
        else if(strncmp("-variable_t",keyword)==0){
          varname_t= strdup(argv[n+1]);
          n++;
          n++;
          continue;
          }
        else if(strncmp("-variable_f",keyword)==0){
          varname_f= strdup(argv[n+1]);
          n++;
          n++;
          continue;
          }

        switch (keyword[1]) {
          case 'h' :
            print_help(argv[0]);
            exit(0);

          case 'r' :
            rootname= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'b' :
            bathymetry.push_back ((string) argv[n+1]);
            n++;
            n++;
            break;

          case 'f' :
            sscanf(argv[n+1],"%lf",&factor);
            n++;
            n++;
            break;

          case 'g' :
            gridfile= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'n' :
            notebook= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'p' :
            PolygonsFile= argv[n+1];
            n++;
            n++;
            break;

          case 's' :
            ShorelinesFile=argv[n+1];
            n++;
            n++;
            break;
        
          default:
            STDOUT_BASE_LINE("unknown option %s\n",keyword);
              exit(-1);
          }
        break;

      default:
        if(notebook==NULL) {
          notebook= strdup(keyword);
          n++;
          }
        else {
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
          }
        break;
      }
    }

  input=new char[1024];
  output=new char[1024];

  if(rootname==NULL) {
    printf("rootname for outputs not specified, use <symphonie>\n");
    rootname= strdup("symphonie");
    }
  else {
    printf("rootname for outputs : %s\n",rootname);
    }

  if(!format) format=strdup(" %6.1f");

  if(PointString!="") {
    status=plg_DecodePosition(PointString, point);
    if(status!=NC_NOERR) NC_TRAP_ERROR(wexit,status,1,"ocean mark prescription"+PointString+" error");
    }

  if(notebook && !gridfile) {
/*-----------------------------------------------------------------------------
    read notebook data and process grid */
    printf("#################################################################\n");
    printf("load notebook file: %s\n",notebook);
    status=load_notebook(notebook, &cgrid, &sgrid, &projection);
    if(status) TRAP_ERR_EXIT(status,"load_notebook(\"%s\",) error\n",notebook);
    printf("notebook successfuly processed\n",notebook);
    }

  if(gridfile) {
/*-----------------------------------------------------------------------------
    read an existing grid */
    if(varname_t==0) TRAP_ERR_EXIT(-1,"please provide a tracer variable name to identify working grid in %s (using command line option: -variable_t <variable name>)\n",gridfile);
    printf("#################################################################\n");
    printf("load tracer grid : file=%s variable=%s\n",gridfile, varname_t);
    status=poc_get_grid(gridfile,varname_t,&sgrid);
    if(status) TRAP_ERR_EXIT(status,"poc_get*grid*(\"%s\",\"%s\",) error\n",gridfile,varname_t);
    status=map_completegridaxis_2(&sgrid);
    printf("tracer grid successfuly processed\n");
    double reference_lat=(sgrid.ymax-sgrid.ymin)/2.0;
    double reference_lon=(sgrid.xmax-sgrid.xmin)/2.0;
    projPJ projection=assign_StereoOblique(reference_lat,reference_lon);
    cgrid=map_get_cartesian(projection,sgrid);
    pj_free(projection);
    skip=0;
    if(varname_f!=0) {
      printf("#################################################################\n");
      printf("load vorticity grid : file=%s variable=%s\n",gridfile, varname_f);
      status=poc_get_grid(gridfile,varname_f,&f_grid);
      if(status) TRAP_ERR_EXIT(status,"poc_get*grid*(\"%s\",\"%s\",) error\n",gridfile,varname_f);
      }
    }
  
  fflush_std();
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  compute grid resolution

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*-----------------------------------------------------------------------------
  compute grid resolution (cartesian)                                        */
  float *cresolution;
  
  status=map_cartesian_resolution(cgrid, &dx, &dy);
  if(status!=0) TRAP_ERR_EXIT(status,"map_cartesian_resolution((modeH=%d),,) error\n",cgrid.modeH);
  cresolution=new float[cgrid.Hsize()];
  for(int m=0;m<cgrid.Hsize();m++) {
    cresolution[m]=sqrt(dx[m]*dy[m]);
    }
  deletep(&dx);
  deletep(&dy);
  
/*-----------------------------------------------------------------------------
  compute grid resolution (spherical)                                        */
  float *sresolution;

  status=map_spherical_resolution(sgrid, &dx, &dy);
  if(status!=0) TRAP_ERR_EXIT(status,"map_spherical_resolution((modeH=%d),,) error\n",sgrid.modeH);
  sresolution=new float[cgrid.Hsize()];
  for(int m=0;m<cgrid.Hsize();m++) {
    sresolution[m]=sqrt(dx[m]*dy[m]);
    }
  deletep(&dx);
  deletep(&dy);
  
/*-----------------------------------------------------------------------------
  initialise resolution structured field                                     */
  SGfield_t<float> fresolution;
  fresolution.x=sresolution;
  fresolution.grid=&sgrid;
  fresolution.mask=-1;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    get landmask from model limits polygons
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  plg_t GridLimits=plg_GridLimits( sgrid,increment);
  TmpPoly=(string)rootname +"-grid.plg";
  status=plg_save(TmpPoly.c_str(), PLG_FORMAT_SCAN, &GridLimits, 1);
  
  if(do_landmask) {
    
    vector<plg_t> polygons;
    
    if(PolygonsFile!="") {
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
      load a model limits polygon file
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

      status=get_polygon(rootname, PolygonsFile, polygons, GridLimits, sgrid, fresolution,
                         wsize, factor, isocontour, filter, point, reduce, debug);
      }
    else {
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
      create a model limits polygon from shoreline dataset
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

      if(isinf(point.x) or isnan(point.x)) {
        TRAP_ERR_EXIT(-1,"\nplease provide a reference position (using command line option: --point \"[<lon>;<lat>]\")\n");
        }
      
      status=create_polygon(rootname, ShorelinesFile, polygons, GridLimits,  sgrid, fresolution,
                            wsize, factor, isocontour, filter, point,  removePools,  show_map, debug);
      }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    then process landmask
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
    
    if(removePools==true and f_grid.modeH==MODEH_UNSET)
      f_grid=map_vgrid(sgrid);
    
    status=process_landmask(rootname, exclusion_polygon, polygons, sgrid, f_grid, checks, landmask, removePools, debug);

    }
  else {
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    or load a landmask buffer from file
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    poc_data_t<float> meta_m;
    
    if(vmask!=0){
      status=meta_m.init(gridfile,vmask);
      NC_TRAP_ERROR(return,status,1,"poc_data_t::init(\"file=%s\",var=\"%s\") error",gridfile,vmask);
      status=meta_m.read_data(gridfile,0,0,1);
      }
    else{
      status=NC_ENOTVAR;
      }
    
    landmask=0;
    if(status==0) {
      swapValues(&landmask,&meta_m.data);
      }
    else {
      landmask=aset(sgrid.Hsize(),-1.f);
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    construct bathymetry
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(bathymetry.size()==0) TRAP_ERR_EXIT(-1,"\nplease provide at least one bathymetry file (using command line option: -b <bathymetry filename>)\n");
  printf("#################################################################\n");
  printf("load and interpolate bathymetry database (priority given in reverse input order)\n");
  fflush_std();
  status=topo_interpolate(bathymetry, sgrid, topo, mask, average);
  if(status) TRAP_ERR_EXIT(status,"topo_interpolate() error\n");
  printf("bathymetry interpolation sucessfully completed\n");

  for(int m=0;m<cgrid.Hsize();m++) {
    sresolution[m]/=1000.0;
    cresolution[m]/=1000.0;
    }

  sprintf(output,"%s.chk.spherical.nc",rootname);
  printf("#################################################################\n");
  printf("create netcdf file (spherical coordinates): %s\n",output);
  fflush_std();
  status=symio_createfile(output,(size_t) sgrid.nx,(size_t) sgrid.ny, sgrid);
  status=symio_writefile(output, landmask, topo, cresolution, sresolution);

  sprintf(output,"%s.chk.mercator.nc",rootname);
  printf("#################################################################\n");
  printf("create netcdf file (mercator coordinates): %s\n",output);
  fflush_std();
  status=symio_createfile(output,(size_t) cgrid.nx,(size_t) cgrid.ny, cgrid);
  status=symio_writefile(output, landmask, topo, cresolution, sresolution);
  
  deletep(&cresolution);
  deletep(&sresolution);

  if(inversion==1) {
    printf("#################################################################\n");
    printf("invert bathymetry sign\n");
    fflush_std();
    for(j=0;j<sgrid.ny;j++)
      for(i=0;i<sgrid.nx;i++) {
        k=j*sgrid.nx+i;
        if (topo[k]!=mask) {
/**------------------------------------------------------------------------------
          set ocean positive depth*/
          topo[k]=-topo[k];
          }
        }
      }

  printf("#################################################################\n");
  printf("check landmask/bathymetry consistency with zmin = %lf\n", zmin);
  fflush_std();
  for(j=0;j<sgrid.ny;j++)
    for(i=0;i<sgrid.nx;i++) {
      k=j*sgrid.nx+i;
      if (topo[k]!=mask) {
// /*------------------------------------------------------------------------------
//         force 0 elevation on land*/
//         if (landmask[k]==-1.) topo[k]=0.0;
// /*------------------------------------------------------------------------------
//         fix land point with ocean depth*/
//         if ((landmask[k]==-1.) && (topo[k]>0.0))  topo[k]=0.0;
/*------------------------------------------------------------------------------
        fix ocean point with land elevation*/
        if ((landmask[k]==1.) && (topo[k]<zmin))  topo[k]=zmin;
        }
      else {
         if (landmask[k]==1.) {
           printf("warning: oceanic point is masked (%lf %lf)\n",sgrid.x[k],sgrid.y[k]);
           }
         }
      }

  sprintf(output,"%s.chk.topo.asc",rootname);
  status= ascii_saver1 (output, cgrid, topo, mask,"%5.0f ",cgrid.nx);

  sprintf(output,"%s.chk.mask.asc",rootname);
  status= ascii_saver1 (output, cgrid, landmask, mask,"%5.0f ",cgrid.nx);

  sprintf(output,"%s.bathycote_in.dat",rootname);
  printf("#################################################################\n");
  printf("create symphonie ASCII file: %s\n",output);
  fflush_std();
  out=fopen(output,"w");

  for(i=skip;i<sgrid.nx-skip;i++) {
    for(j=skip;j<sgrid.ny-skip;j++) {
      k=sgrid.nx*j+i;
      switch ((int) landmask[k]) {
        case -1:
          flag=0;
          break;
        case 0:
          flag=0;
          break;
        case 1:
          flag=1;
          break;
        default:
          TRAP_ERR_EXIT(-1,"invalid mask value\n");
          break;
        }
      fprintf(out,"%1d",flag);
      }
    fprintf(out,"\n");
    }
  for(i=skip;i<sgrid.nx-skip;i++) {
    for(j=skip;j<sgrid.ny-skip;j++) {
      k=sgrid.nx*j+i;
      if (topo[k]==mask)  topo[k]=0.0;
      fprintf(out,format,topo[k]);
      }
    fprintf(out,"\n");
    }

  fclose(out);
  
  

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  create T-UGOm/sts compatible SG grid
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  {
  bool standard=false;
  
  printf("#################################################################\n");
  if (standard) {
    printf("landmask is on F grid : deriving T grid\n");
    z_grid=map_f2zgrid(sgrid);
    f_grid=sgrid;
    }
  else {
    printf("landmask is on T grid : deriving F grid\n");
    z_grid=sgrid;
    f_grid=map_vgrid(sgrid);
    }
  
  string filename=(string)rootname+"-symtools-comodo.nc";
  printf("create : %s\n",filename.c_str());
  fflush_std();

  f_gridVar=poc_save_grid(filename, f_grid, __FILE__,__LINE__,"f");
  
  z_gridVar=poc_save_grid(filename, z_grid, 0, 0, "t");
  
  delete[] topo;
  status=topo_interpolate(bathymetry, z_grid, topo, mask, average);
  if(status) TRAP_ERR_EXIT(status,"topo_interpolate() error\n");
  
  if(inversion==1) {
    for(j=0;j<z_grid.ny;j++)
      for(i=0;i<z_grid.nx;i++) {
        k=j*z_grid.nx+i;
        if (topo[k]!=mask) {
/*------------------------------------------------------------------------------
          set ocean positive depth*/
          topo[k]=-topo[k];
          }
        }
    }

  poc_var_t vtopo;
  
  vtopo=z_gridVar;
  vtopo.init("bathymetry_t",NC_FLOAT,"bathymetry_t","m",mask);
  vtopo << poc_att_t("axis","YX");
  status=poc_put_vara(filename,vtopo,0,topo,0);
  
  delete[] topo;
  status=topo_interpolate(bathymetry, f_grid, topo, mask, average);
  if(status) TRAP_ERR_EXIT(status,"topo_interpolate() error\n");
  
  if(inversion==1) {
    for(j=0;j<f_grid.ny;j++)
      for(i=0;i<f_grid.nx;i++) {
        k=j*f_grid.nx+i;
        if (topo[k]!=mask) {
/**------------------------------------------------------------------------------
          set ocean positive depth*/
          topo[k]=-topo[k];
          }
        }
    }
  
  vtopo=f_gridVar;
  vtopo.init("bathymetry_f",NC_FLOAT,"bathymetry_f","m",mask);
  poc_put_vara(filename,vtopo,0,topo,0);
  
  float *tlandmask=0;
  
  if(standard==true){
    printf("deriving T landmask from F landmask\n");
    
    tlandmask=new float[z_grid.Hsize()];
    
    for(j=0;j<z_grid.ny;j++) {
      for(i=0;i<z_grid.nx;i++) {
        k=j*z_grid.nx+i;
        
        int n1=j*f_grid.nx+i;
        int n2=(j+1)*f_grid.nx+i;
        int n3=(j+1)*f_grid.nx+i+1;
        int n4=j*f_grid.nx+i+1;
        
        if (landmask[n1]==1 and landmask[n2]==1 and landmask[n3]==1 and landmask[n4]==1) {
          tlandmask[k]=1;
          }
        else {
          tlandmask[k]=-1;
          }
        
        
        }
      }
    
    }
  else{
    tlandmask=landmask;
    
    printf("deriving F landmask from T landmask\n");
    
    landmask=aset(f_grid.Hsize(),-1.f);
    
    for(j=0;j<z_grid.ny;j++) {
      for(i=0;i<z_grid.nx;i++) {
        k=j*z_grid.nx+i;
        
        if(tlandmask[k]==-1)
          continue;
        
        int n1=j*f_grid.nx+i;
        int n2=(j+1)*f_grid.nx+i;
        int n3=(j+1)*f_grid.nx+i+1;
        int n4=j*f_grid.nx+i+1;
        
        landmask[n1]=1;
        landmask[n2]=1;
        landmask[n3]=1;
        landmask[n4]=1;
        }
      }
    
    }
  
  poc_var_t vmask;
  
  vmask=z_gridVar;
  vmask.init("landmask_t",NC_BYTE,"landmask_t","dimensionless",2);
  poc_put_vara(filename,vmask,0,tlandmask,0);
  
  vmask=f_gridVar;
  vmask.init("landmask_f",NC_BYTE,"landmask_f","dimensionless",2);
  poc_put_vara(filename,vmask,0,landmask,0);
  
//   poc_var_t vresolution=gridVar;
//   vresolution.init("resolution",NC_FLOAT,"resolution","km");
//   poc_put_vara(filename,vresolution,resolution,1);
  if (standard) {
    z_grid.free();
    }
  else {
    f_grid.free();
    }
  
  delete[] tlandmask;
  }
  
  cgrid.free();
  sgrid.free();
  delete[] topo;
  delete[] landmask;
  
  deletep(&output);
  deletep(&input);

  printf("#################################################################\n");
  STDOUT_BASE_LINE("symtools successfuly completed\n");
  
  return 0;
}
