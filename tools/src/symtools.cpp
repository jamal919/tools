
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
  int status;

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
  size_t size;
  bool repair=true;

  int inFormat=plg_find_format(ShorelinesFile.c_str());
  
  status=plg_load(ShorelinesFile.c_str(), inFormat, shorelines_base);
  if(status!=0) return(extraction);
  
  extraction=plg_extract_out(shorelines_base, selection, rootname, mode, debug);
  
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
    "  -v  followed by the grid variable name\n"
    "  -n  followed by the path to the notebook file\n"
    "  -p  followed by a path to a land mask polygons filename. Disables -s. If not given, compute one from the shorelines dataset (see below).\n"
    "  -s  followed by the path to a shorelines filename\n"
    "  --point  followed by the coordinates of a reference as \"[<lon>;<lat>]\"\n"
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
  int flag;
  char *rootname=NULL,*output=NULL,*input=NULL;
  char *notebook=NULL,*poly=NULL,*gridfile=NULL,*format=NULL,*varname=NULL;
  string exclusion_polygon;
  vector<string> bathymetry;
  grid_t cgrid,sgrid;

  float  *landmask,*topo;
  cdfgbl_t data_info,grid_info;
  cdfvar_t h_info;

  int inversion=1,average=0;
  float zmin=-99999;
  double *dx, *dy, factor=1.000, wsize=0;
  
  int skip=1;
  
  bool do_landmask=true, do_bathymetry=true;
  char *vmask=0, *vtopo=0;
  poc_data_t<float>  meta_H_z, meta_H_f;

  string ShorelinesFile, TmpPoly, PointString;
  point2D_t point;
  double isocontour=-0.4;
  bool checks=false, filter=false, reduce=true;
  bool show_map=false;
  bool debug=false;
  
  fct_echo( argc, argv);
  
  if(argc<2) {
    print_help(argv[0]);
    exit(0);
    }
    
  printf("\n\n\n WARNING : symtools will soon be deprecated, do rather use grid-config (check for some minor arguments changes)\n\n\n");
  
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
          sscanf(argv[n+1],"%f",&isocontour);
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
            poly= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 's' :
            ShorelinesFile=argv[n+1];
            n++;
            n++;
            break;
        
          case 'v' :
            varname= strdup(argv[n+1]);
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

  input=(char *)malloc(1024);
  output=(char *)malloc(1024);

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
    if(varname==0) TRAP_ERR_EXIT(status,"please provide a variable name to identify working grid in %s (using command line option: -v <variable name>)\n",gridfile);
    printf("#################################################################\n");
    printf("load grid file: %s\n",gridfile);
    status=poc_get_grid(gridfile,varname,&sgrid);
//     status=poc_getgrid2d(gridfile,varname,&sgrid);
    if(status) TRAP_ERR_EXIT(status,"poc_get*grid*(\"%s\",\"%s\",) error\n",gridfile,varname);
    printf("grid successfuly processed\n");
    double reference_lat=(sgrid.ymax-sgrid.ymin)/2.0;
    double reference_lon=(sgrid.xmax-sgrid.ymin)/2.0;
    projPJ projection=assign_StereoOblique(reference_lat,reference_lon);
    cgrid=map_get_cartesian(projection,sgrid);
    skip=0;
    }
  
/*-----------------------------------------------------------------------------
  compute grid resolution (cartesian) */
  float *cresolution;
  
  status=map_cartesian_resolution(cgrid, &dx, &dy);
  cresolution=new float[cgrid.Hsize()];
  for(int m=0;m<cgrid.Hsize();m++) {
    cresolution[m]=sqrt(dx[m]*dy[m]);
    }
  delete[] dx;
  delete[] dy;
  
/*-----------------------------------------------------------------------------
  compute grid resolution (spherical) */
  float *sresolution;

  status=map_spherical_resolution(sgrid, &dx, &dy);
  sresolution=new float[cgrid.Hsize()];
  for(int m=0;m<cgrid.Hsize();m++) {
    sresolution[m]=sqrt(dx[m]*dy[m]);
    }
  delete[] dx;
  delete[] dy;
  
/*-----------------------------------------------------------------------------
  initialise resolution structured field */
  SGfield_t<float> fresolution;
  fresolution.x=sresolution;
  fresolution.grid=&sgrid;
  fresolution.mask=-1;
  
  plg_t GridLimits=plg_GridLimits(sgrid, 1);
  TmpPoly=(string)rootname +"-grid.plg";
  status=plg_save(TmpPoly.c_str(), PLG_FORMAT_SCAN, &GridLimits, 1);
  
  if(do_landmask) {
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    get landmask from model limits polygons
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    vector<plg_t> polygons;
    
    if(poly!=NULL) {
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
      either load a dedicated polygon file
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

      sprintf(input,"%s",poly);
      printf("#################################################################\n");
      printf("load landmask polygons file: %s\n",poly);
      status=plg_load(input, PLG_FORMAT_SCAN, polygons);
      if(status) TRAP_ERR_EXIT(status,"plg_load(\"%s\",,) error\n",input);
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
        reduce polygon set to grid extent to minimize computational costs
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

      if(reduce) {

        frame_t GridFrame=plg_spherical_minmax(&GridLimits,1);
        double dilatation=0.05;
        dilatation=1./(1.* min(sgrid.nx,sgrid.ny) );
      
        GridFrame.dilatation(dilatation);

        plg_t ExtendedLimits=plg_t(GridFrame, PLG_INIT_SEPARATE);
      
        int target=0;

        vector<plg_t>dilated;
        dilated.push_back(ExtendedLimits);
        vector<plg_t> extraction=plg_extract(polygons, dilated, rootname, PLG_SPHERICAL, target, debug);
        plg_destroy_entries(polygons);
        polygons=extraction;
        TmpPoly=(string) rootname +"-boundaries-reduced.plg";
        printf("#################################################################\n");
        printf("save reduced boundary polygons : %s\n",TmpPoly.c_str());
        status=plg_save(TmpPoly.c_str(), PLG_FORMAT_SCAN, polygons);
        }
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
        filter polygon set
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

      if(filter) {
        vector<float> isovalues;
        for(int k=0;k<7;k++) isovalues.push_back(0.2-k*0.1);
            
        int target=2;
        double size=wsize, lengthscale=0.0;
        bool enforcement=false, tuned=true, show_mesh=false, subdivide=false;
        if(debug) show_map=true;
    
        plg_filter_t plg_filter;
        plg_filter.init(lengthscale, fresolution, factor, size, subdivide, tuned, enforcement, show_map, show_mesh);
    
        vector<plg_t> smoothed=plg_smooth(rootname, polygons, plg_filter, isovalues, isocontour, point, target, debug);
            
        plg_destroy_entries(polygons);
        polygons=smoothed;
        printf("#################################################################\n");
        STDOUT_BASE_LINE("remove swimming pools\n");
        status=plg_RemoveInlandLimits(polygons, point.x, point.y, (projPJ) 0, PLG_POINT_EXTERIOR, debug);
        TmpPoly=(string) rootname +"-boundaries-filtered.plg";
        printf("#################################################################\n");
        printf("save filtered boundary polygons : %s\n",TmpPoly.c_str());
        status=plg_save(TmpPoly.c_str(), PLG_FORMAT_SCAN, polygons);
        }
      }
    else {
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
      or create a polygon from shoreline dataset
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
      if(isinf(point.x) or isnan(point.x)) {
        TRAP_ERR_EXIT(-1,"\nplease provide a reference position (using command line option: --point \"[<lon>;<lat>]\")\n");
        }
        
      frame_t GridFrame=plg_spherical_minmax(&GridLimits,1);
      double dilatation=0.05;

      dilatation=1./(1.* min(sgrid.nx,sgrid.ny) );
      GridFrame.dilatation(dilatation);

      plg_t ExtendedLimits=plg_t(GridFrame, PLG_SPHERICAL);
      TmpPoly=(string) rootname +"-frame.plg";
      status=plg_save(TmpPoly.c_str(), PLG_FORMAT_SCAN, &ExtendedLimits, 1);
      
      bool enforcement=false, tuned=false, show_mesh=false, subdivide=false;
      if(debug) show_map=true;
      double lengthscale=0.0;

      plg_filter_t plg_filter;
      plg_filter.init(lengthscale, fresolution, factor, wsize, subdivide, tuned, enforcement, show_map, show_mesh);
       
      if(ShorelinesFile=="") {
        ShorelinesFile="/home/data/shorelines/gshhs-2.2/gshhs_f.cst";
        printf("shorelines filename not given, use default : %s\n",ShorelinesFile.c_str());
        }

      
      status=plg_CreateLimits(ShorelinesFile, ExtendedLimits, rootname, plg_filter, point, isocontour, filter, true, debug);
//       status=plg_CreateLimits(ShorelinesFile, ExtendedLimits, rootname, fresolution, factor, wsize, point, isocontour, enforcement, filter, show_map, debug);
      
      if(status) TRAP_ERR_EXIT(status,"plg_CreateLimits() error\n");
      TmpPoly=(string) rootname +"-boundaries.plg";
      status=plg_load(TmpPoly.c_str(), PLG_FORMAT_SCAN, polygons);
      if(status) TRAP_ERR_EXIT(status,"plg_load(\""+TmpPoly+"\",,) error\n");
      }

/*------------------------------------------------------------------------------
    exclusion polygons */
    if(exclusion_polygon!="") {
      vector<plg_t> selection;
      status=plg_load(exclusion_polygon, PLG_FORMAT_SCAN, selection);
      vector<plg_t> extraction=plg_extract_out(polygons, selection, rootname, PLG_SPHERICAL, debug);
      plg_destroy_entries(polygons);
      polygons=extraction;
      TmpPoly=(string) rootname +"-with-exclusion.plg";
//       printf("#################################################################\n");
//       printf("save reduced boundary polygons : %s\n",TmpPoly.c_str());
      status=plg_save(TmpPoly.c_str(), PLG_FORMAT_SCAN, polygons);
      }
/*------------------------------------------------------------------------------
    compute landmask from polygons */    
    landmask=set_landmask02(sgrid,polygons,checks);
    plg_destroy_entries(polygons);
    printf("land mask sucessfully completed\n");
//     extern int fe_ZGrid2Quadrangle(const grid_t & zgrid, const grid_t & fgrid, char  *landmask, float *topo, float topomask, mesh_t & mesh, int FixPinched, int verbose);

    }
  else {
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    or load a mask array form file
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    poc_data_t<float> meta_m;
    
    status=meta_m.init(gridfile,vmask);
    NC_TRAP_ERROR(return,status,1,"poc_data_t::init(\"file=%s\",var=\"%s\") error",gridfile,vmask);
    status=meta_m.read_data(gridfile,0,0,1);
    
    landmask=0;
    if(status==0) {
      swapValues(&landmask,&meta_m.data);
      }
    else {
      landmask=aset(sgrid.Hsize(),-1.f);
      }
    }
  
/*-----------------------------------------------------------------------------
  construct bathymetry */
  if(bathymetry.size()==0) TRAP_ERR_EXIT(-1,"\nplease provide at least one bathymetry file (using command line option: -b <bathymetry filename>)\n");
  printf("#################################################################\n");
  printf("load and interpolate bathymetry database (priority given in reverse input order)\n");
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
  status=symio_createfile(output,(size_t) sgrid.nx,(size_t) sgrid.ny, sgrid);
  status=symio_writefile(output, landmask, topo, cresolution, sresolution);

  sprintf(output,"%s.chk.mercator.nc",rootname);
  printf("#################################################################\n");
  printf("create netcdf file (mercator coordinates): %s\n",output);
  status=symio_createfile(output,(size_t) cgrid.nx,(size_t) cgrid.ny, cgrid);
  status=symio_writefile(output, landmask, topo, cresolution, sresolution);
  
  deletep(&cresolution);
  deletep(&sresolution);

  if(inversion==1) {
    printf("#################################################################\n");
    printf("invert bathymetry sign\n");
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
  poc_var_t zgridVar, fgridVar;
  grid_t zgrid,fgrid;
  
  bool standard=true;
  if (standard) {
    zgrid=map_f2zgrid(sgrid);
    fgrid=sgrid;
    }
  else {
    zgrid=sgrid;
    fgrid=map_vgrid(sgrid);
    }
  
  string filename=(string)rootname+"-symtools-comodo.nc";
  printf("#################################################################\n");
  printf("create : %s\n",filename.c_str());

  fgridVar=poc_save_grid(filename.c_str(), fgrid, __FILE__,__LINE__,"f");
  
  zgridVar=poc_save_grid(filename.c_str(), zgrid, 0, 0, "t");
  
  poc_var_t vtopo=zgridVar;
  poc_var_t vmask=zgridVar;
  
  delete[] topo;
  status=topo_interpolate(bathymetry, zgrid, topo, mask, average);
  if(status) TRAP_ERR_EXIT(status,"topo_interpolate() error\n");
  
  if(inversion==1) {
    for(j=0;j<zgrid.ny;j++)
      for(i=0;i<zgrid.nx;i++) {
        k=j*zgrid.nx+i;
        if (topo[k]!=mask) {
/*------------------------------------------------------------------------------
          set ocean positive depth*/
          topo[k]=-topo[k];
          }
        }
    }

  vtopo.init("bathymetry_t",NC_FLOAT,"bathymetry_t","m",mask);  
  
  vtopo<<poc_att_t("axis","YX");

  status=poc_put_vara(filename,vtopo,0,topo,0);
  
//   vmask.init("landmask_t",NC_FLOAT,"landmask","dimensionless",2);
//   poc_put_vara(filename,vmask,landmask,0);
  
  delete[] topo;
  status=topo_interpolate(bathymetry, fgrid, topo, mask, average);
  if(status) TRAP_ERR_EXIT(status,"topo_interpolate() error\n");
  
  if(inversion==1) {
    for(j=0;j<fgrid.ny;j++)
      for(i=0;i<fgrid.nx;i++) {
        k=j*fgrid.nx+i;
        if (topo[k]!=mask) {
/**------------------------------------------------------------------------------
          set ocean positive depth*/
          topo[k]=-topo[k];
          }
        }
    }
  
  vtopo=fgridVar;
  vtopo.init("bathymetry_f",NC_FLOAT,"bathymetry_f","m",mask);
  poc_put_vara(filename.c_str(),vtopo,0,topo,0);
  
  float *tlandmask=new float[zgrid.Hsize()];
  for(j=0;j<zgrid.ny;j++) {
    for(i=0;i<zgrid.nx;i++) {
      k=j*zgrid.nx+i;
      int n1=j*fgrid.nx+i;
      int n2=(j+1)*fgrid.nx+i;
      int n3=(j+1)*fgrid.nx+i+1;
      int n4=j*fgrid.nx+i+1;
      if (landmask[n1]==1 and landmask[n2]==1 and landmask[n3]==1 and landmask[n4]==1) {
        tlandmask[k]=1;
        }
      else {
        tlandmask[k]=-1;
        }
      }
    }
  
  vmask.init("landmask_t",NC_BYTE,"landmask_t","dimensionless",2);
  poc_put_vara(filename.c_str(),vmask,0,tlandmask,0);
  
  vmask=fgridVar;
  vmask.init("landmask_f",NC_BYTE,"landmask_f","dimensionless",2);
  poc_put_vara(filename.c_str(),vmask,0,landmask,0);
  
//   poc_var_t vresolution=gridVar;
//   vresolution.init("resolution",NC_FLOAT,"resolution","km");
//   poc_put_vara(filename,vresolution,resolution,1);
  if (standard) {
    zgrid.free();
    }
  else {
    fgrid.free();
    }
  
  }
  
  cgrid.free();
  sgrid.free();
  delete[] topo;
  delete[] landmask;

  printf("#################################################################\n");
  STDOUT_BASE_LINE("symtools successfuly completed\n");
  
  return 0;
}
