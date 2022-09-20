
/*******************************************************************************

  T-UGO tools, 2006-2016

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada
\author  Clement MAYET      LEGOS, Toulouse, France (PhD)

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Convert between different mesh formats.

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
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>


#include "tools-structures.h"

#include "geo.h"
#include "fe.h"
#include "map.h"
#include "polygones.h"
#include "functions.h"
#include "swap.h"
#include "poc-time.h"
#include "archive.h"


extern int fe_readQmesh_ascii(const string &polygonsFileName, mesh_t *mesh);


// a map containing mesh_file_format codes (integer) as defined in fe.h, indexed by mesh format names (string)
  map<string,int> mesh_format_map;


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
    " NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s -i input_filename -o output_filename [OPTIONS]\" \n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "   Convert between different mesh formats.\n");
  printf("\n"
   "OPTIONS\n"
    "  -h,--help  print this help and exit\n"
    +isFormat_help("  ")+
    "  -b  followed by path to boundary elements file.\n"
    "  -d  followed by path to .s2r bathymetry file. Only relevant if the output can contain the bathymetry.\n"
    "  --inv  bathymetry is positive downward. Default is negative downward. Will only be used to overwrite boundary codes.\n"
    "  -i  followed by path to the input mesh\n"
    "  -o  followed by path to the output mesh\n"
    "  -r  run 'fe_reducebw(mesh, 100);'\n"
    "  -c  create mode (Default)\n"
    "  -u  update mode, updates (lon,lat) variables in output file from input file\n"
   "\n"
   "FORMATS\n"
   "  Defaults are 'nei' if mesh is .nei or 'nc2d' if mesh is .nc .\n"
   "  Format names are: " + get_key_list(mesh_format_map) + " .\n"
   "In these previous names:\n"
   "- 'nei' stands for neighbour\n"
   "- 'nc' stands for NetCDF\n"
   "- 'gmsh' stands for gmesh (see its documentation)\n"
   "- 'ww' stands for WaveWatch III (see its documentation)\n"
   );
  
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_V2cartesian(mesh_t & mesh, const char *proj4)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
{
  int status;
  int ndata;
  double *x,*y;
  
  ndata=mesh.nvtxs;
  
  x=new double[ndata];
  y=new double[ndata];
  
  for(size_t n=0;n<mesh.nvtxs;n++) {
    x[n]=mesh.vertices[n].lon;
    y[n]=mesh.vertices[n].lat;
    }
    
  status=projection_to_geo (proj4, x, y, ndata);
  
  for(size_t n=0;n<mesh.nvtxs;n++) {
    mesh.vertices[n].lon=x[n];
    mesh.vertices[n].lat=y[n];
    }
  
  delete[] x;
  delete[] y;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_V2spherical(mesh_t & mesh, const char *proj4)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int ndata;
  double *x,*y;
  
  ndata=mesh.nvtxs;
  
  x=new double[ndata];
  y=new double[ndata];
  
  for(size_t n=0;n<mesh.nvtxs;n++) {
    x[n]=mesh.vertices[n].lon;
    y[n]=mesh.vertices[n].lat;
    }
    
  status=projection_to_geo (proj4, x, y, ndata);
  
  for(size_t n=0;n<mesh.nvtxs;n++) {
    mesh.vertices[n].lon=x[n];
    mesh.vertices[n].lat=y[n];
    }
  
  delete[] x;
  delete[] y;
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int parse_GridOptions(const string options, char* & gridfile, char* & maskfile, metagrid_t & meta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  example file=ORCA025-TRBB36s001_mesh_hgr.nc vlon_f=glamba

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  vector<string> tokens,keys,values;
  int k;
  string delimiter=" ";
  
  gridfile=0;
  
  delimiter=" ";
  tokens=string_split(options, delimiter);
  
  delimiter="=";
  for(k=0;k<tokens.size();k++) {
    vector<string> tmp=string_split(tokens[k], delimiter);
    keys.push_back(tmp[0]);
    values.push_back(tmp[1]);
    }
  
  for(k=0;k<keys.size();k++) {
    if(keys[k]=="file") {
      meta.gridfile=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="maskfile") {
      meta.maskfile=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="topofile") {
      meta.topofile=values[k];
      continue;
      }
    if(keys[k]=="vlon_t") {
      meta.z_gnames.vlon=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="vlat_t") {
      meta.z_gnames.vlat=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="vmask_t") {
      meta.z_gnames.vmask=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="vtopo_t") {
      meta.z_gnames.vtopo=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="vlon_f") {
      meta.f_gnames.vlon=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="vlat_f") {
      meta.f_gnames.vlat=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="vmask_f") {
      meta.f_gnames.vmask=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="vtopo_f") {
      meta.f_gnames.vtopo=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="projection") {
      meta.projection=strdup(values[k].c_str());
      continue;
      }
    printf("\n%s: keyword <%s> not reckognized \n\n", __func__, keys[k].c_str());
    TRAP_ERR_EXIT(-1,"exiting\n");
    }
  
  if(meta.maskfile=="") meta.maskfile=meta.gridfile.c_str();
  
  maskfile=strdup(meta.maskfile.c_str());
  gridfile=strdup(meta.gridfile.c_str());

//   if(meta.topofile=="") meta.topofile=gridfile;
    
  return(0);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,fmt_in=MESH_FILE_FORMAT_UNKNOWN, fmt_out=MESH_FILE_FORMAT_UNKNOWN;
  int n;
  int reduce=0, save_open_boundary=0, save_bathy=0;
  char *meshfile=NULL, *depthfile=NULL, *belfile=NULL, *output=NULL, *iproj4=NULL, *oproj4=NULL, *wproj4=NULL;
  mesh_t mesh;
  ostringstream report;
  int mode=0;
  double depth_conv=1;
  string input_format_name, output_format_name;
  int stopon_EdgeError, stopon_PinchError;
  char *comments[2];
  int element=FE_TRIANGLE;
  vector<plg_t> polygons;
  string grid_options="";
  double tag=NAN;
  bool allow_dry=true;

  fct_echo( argc, argv);

  mesh_format_map["nc2d"]            =MESH_FILE_FORMAT_NC2D;
  mesh_format_map["nc3d"]            =MESH_FILE_FORMAT_NC3D;
  mesh_format_map["nei"]             =MESH_FILE_FORMAT_TRIGRID;
  mesh_format_map["gmsh"]            =MESH_FILE_FORMAT_GMSH;
  mesh_format_map["gmsh_ww"]         =MESH_FILE_FORMAT_GMSH_WW;
  mesh_format_map["schism"]          =MESH_FILE_FORMAT_SCHISM;
  mesh_format_map["quadrangle"]      =MESH_FILE_FORMAT_QUADRANGLE_ASCII;
  mesh_format_map["telemac_ascii"]   =MESH_FILE_FORMAT_TELEMAC_ASCII;
  mesh_format_map["telemac_binary"]  =MESH_FILE_FORMAT_TELEMAC_BINARY;
  mesh_format_map["telemac_swapped"] =MESH_FILE_FORMAT_TELEMAC_SWAPPED;
  mesh_format_map["quoddy"]          =MESH_FILE_FORMAT_QUODDY;
  mesh_format_map["netcdf"]          =MESH_FILE_FORMAT_STRUCTURED_NETCDF;

  n=1;
  while (n < argc) {
    const char *keyword=argv[n];
    
    if(isFormat(argv,&n,&input_format_name,&output_format_name)){
      continue;
      }
      
    if(strcmp(argv[n],"-g")==0){
      grid_options=strdup(argv[n+1]);
      n++;
      n++;
      continue;
      }
   
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'i' :
/*-----------------------------------------------------------------------------
          input mesh filename*/
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'b' :
/*-----------------------------------------------------------------------------
          input bel file (boundary elements description)*/
          belfile= strdup(argv[n+1]);
          save_open_boundary=1;
          n++;
          n++;
          break;

        case 'd' :
/*-----------------------------------------------------------------------------
          input depth file (s2r)*/
          depthfile= strdup(argv[n+1]);
          save_bathy=1;
          n++;
          n++;
          break;

        case 'o' :
/*-----------------------------------------------------------------------------
          output mesh filename*/
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case '-' :
/*-----------------------------------------------------------------------------
          help */
          if(strcmp(keyword,"--help")==0){
            print_help(argv[0]);
            exit(0);
            }
          
          if(strcmp(keyword,"--quadrangle")==0){
            element=FE_QUADRANGLE;
            n++;
            continue;
            }
          
          if(strcmp(keyword,"--inv")==0){
/*-----------------------------------------------------------------------------
            inverse depth convention */
            depth_conv= -1;
            n++;
            break;
            }
          if(strcmp(keyword,"--dry=no")==0){
/*-----------------------------------------------------------------------------
             */
            allow_dry=false;
            n++;
            break;
            }
          if(strcmp(keyword,"--dry=yes")==0){
/*-----------------------------------------------------------------------------
             */
            allow_dry=true;
            n++;
            break;
            }
          if(strstr(keyword,"--wproj=")!=0){
            wproj4=strdup(keyword);
            n++;
            break;
            }
          if(strstr(keyword,"--iproj=")!=0){
            iproj4=strdup(keyword);
            n++;
            break;
            }
          if(strstr(keyword,"--oproj=")!=0){
            oproj4=strdup(keyword);
            n++;
            break;
            }
          break;

        case 'c' :
/*-----------------------------------------------------------------------------
          create mode*/
          mode=0;
          n++;
          break;

        case 'u' :
/*-----------------------------------------------------------------------------
          update mode*/
          mode=1;
          n++;
          break;

        case 'r' :
          reduce=1;
          n++;
          break;

        case 'h' :
/*-----------------------------------------------------------------------------
          help */
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
        printf("unknown option %s\n",keyword);
        print_help(argv[0]);
        exit(-1);
      }
    
    }

//  status= fe_readQmesh_ascii( meshfile, &mesh);

  if(meshfile==NULL and grid_options=="") {
    printf("no input mesh file specified; abort...\n");
    print_help(argv[0]);
    exit(-1);
    }

  if(output == NULL) {
    printf("no output mesh file specified; abort...\n");
    print_help(argv[0]);
    exit(-1);
    }
  
  if(output_format_name!="")
    fmt_out = mesh_format_map[output_format_name];
  
  if(input_format_name!="")
    fmt_in = mesh_format_map[input_format_name];
  
  if(fmt_in==MESH_FILE_FORMAT_UNKNOWN)
    fmt_in=fe_find_format(meshfile);
  
  if(fmt_out==MESH_FILE_FORMAT_UNKNOWN)
    fmt_out=fe_find_format(output);
 
  if(fmt_in==MESH_FILE_FORMAT_UNKNOWN or fmt_out==MESH_FILE_FORMAT_UNKNOWN) {
    if(fmt_in==MESH_FILE_FORMAT_UNKNOWN)
      printf("no input mesh format specified; abort...\n");
    if(fmt_out==MESH_FILE_FORMAT_UNKNOWN)
      printf("no output mesh format specified; abort...\n");
    print_help(argv[0]);
    exit(-1);
    }
 

/* *----------------------------------------------------------------------------
  MODULEF P2 ASCII format*/
//  status=fe_readmesh_MODULEFP2(meshfile,&mesh);
/* *----------------------------------------------------------------------------
  Quadrangle format*/
//  status= fe_readQmesh_ascii( meshfile, &mesh);

  printf("#################################################################\n");
  printf("load mesh from %s\n",meshfile);
  if(grid_options!="") {
    char *gridfile=0, *maskfile=0;
    metagrid_t meta;
    grid_t grid;
    mesh_t mesh;
    string rootname="mesh-format";
    bool debug=true;
    
    status= parse_GridOptions(grid_options, gridfile, maskfile, meta);
    if(status!=0) {
      TRAP_ERR_EXIT(status,"parse_GridOptions(\""+grid_options+"\",\"%s\",\"%s\",) error\n", gridfile, maskfile);
      }
    status=quadrangle_ImportStructured(gridfile, maskfile, mesh, meta, tag, rootname, wproj4, debug);
    if(status!=0) {
      TRAP_ERR_EXIT(status,"quadrangle_ImportStructured(\"%s\",\"%s\",,,%g,\""+rootname+"\",%d) error\n", gridfile, maskfile, tag, debug);
      }
    meshfile=strdup("mesh-format-quadrangle.nei");
    status=fe_savemesh(meshfile,MESH_FILE_FORMAT_TRIGRID, mesh);
    }
  else {
    status=fe_readmesh(meshfile,fmt_in, &mesh, iproj4);
    if(status!=0) TRAP_ERR_EXIT(status,"fe_readmesh(\"%s\",,,) error %d\n",meshfile,status);
    }
    
  if(fmt_in == MESH_FILE_FORMAT_TRIGRID) {
/*-----------------------------------------------------------------------------
    build element list from the original neighbour list */
    status=fe_list(&mesh, element);
    if(status!=0) TRAP_ERR_EXIT(status,"fe_list(,) error %d\n",status);
    }

  if(    fmt_in == MESH_FILE_FORMAT_GMSH
      or fmt_in == MESH_FILE_FORMAT_GMSH_WW
      or fmt_in == MESH_FILE_FORMAT_TELEMAC_ASCII
      or fmt_in == MESH_FILE_FORMAT_TELEMAC_BINARY
      or fmt_in == MESH_FILE_FORMAT_SCHISM ) {
/*-----------------------------------------------------------------------------
    gmesh and telemac formats do not contain a list of neighbours, we need to create it */
    status= fe_e2n(&mesh);
    if(status!=0) TRAP_ERR_EXIT(status,"fe_e2n() error %d\n",status);
    }
 
  status=fe_savemesh("check.nei", MESH_FILE_FORMAT_TRIGRID, mesh, oproj4);

  status=fe_connex(mesh);

  status= fe_edgetable(&mesh,0,0);
  if(status!=0) {
    STDERR_BASE_LINE("fe_edgetable(,0,0) error %d. Attempting to save to %s\n",status,output);
    status=fe_savemesh(output,fmt_out,mesh);
    wexit(-1);
    }
  
  printf("#################################################################\n");
  printf("build additional mesh tables/information\n");
  status= fe_vertex_element_tables(&mesh);
  if(status!=0) TRAP_ERR_EXIT(status,"fe_vertex_element_tables() error %d\n",status);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  builds a table containing geometric code (interior, extern boundary, island,
  etc...) from mesh geometry  
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  stopon_EdgeError=1;
  stopon_PinchError=0;
  status= fe_codetable1(&mesh, 0, stopon_EdgeError, stopon_PinchError);
  if(status!=0) {
    TRAP_ERR_EXIT(status,"fe_codetable1(,0,%d,%d) error %d\n", stopon_EdgeError, stopon_PinchError,status);
    }

  status= fe_bandwidth(&mesh);

  report << "#Number of elements   : " << mesh.ntriangles <<endl;
  report << "#Half-bandwidth       : " << mesh.hbw <<endl;
  report << "#Number of edges      : " << mesh.nedges <<endl;
  cout << report.str() << flush;

//  if(save_bathy==1 && fmt_in==MESH_FILE_FORMAT_TRIGRID) {
  if(save_bathy==1) {
    printf("#################################################################\n");
    printf("add LGP1 bathymetry to mesh vertex\n");
    float* depth_buffer=NULL;
    exitIfNull(
      depth_buffer = new float[mesh.nvtxs]
      );
    status=quoddy_loadr1(depthfile, mesh.nvtxs,depth_buffer); // read bathy (s2r) and write it in mesh_t
    for(int i=0; i<mesh.nvtxs; i++){
      mesh.vertices[i].h=depth_buffer[i];
      }
    delete [] depth_buffer;
    }

  if(save_bathy==1 && fmt_in==MESH_FILE_FORMAT_GMSH_WW) {
    printf("#################################################################\n");
    printf("add LGP1 bathymetry to mesh vertex\n");
    float* depth_buffer=NULL;
    exitIfNull(
      depth_buffer = new float[mesh.nvtxs]
      );

    for(int i=0; i<mesh.nvtxs; i++){
      depth_buffer[i]=mesh.vertices[i].h;
      }

    comments[0]=new char[1024];
    comments[1]=new char[1024];
    sprintf(comments[0]," bathymetry, unit in meters");
    sprintf(comments[1],"Created by mesh-format from %s",meshfile);
    status=quoddy_saver1(depthfile, mesh.nvtxs,depth_buffer,comments); // read bathy from gmsh meshfile and write it in s2r
    delete [] depth_buffer;
    }

  if(save_open_boundary==1) {
    status=fe_read_boundarycode(belfile,mesh,2) ; // read bel file and write it in mesh_t limits
    if(status!=0) TRAP_ERR_EXIT(status,"fe_read_boundarycode(\"%s\",,2) error %d\n",status);
    n=0;
    /** ***************************************************
     \bug 22 mars 2013 :  Clement MAYET :     some nodes already have code 5 (from nei mesh file),
      and there is conflict with belfile boundary code definition. Overwriting (lines below is necessary.
    ***************************************************** */

    for (int i=0;i<mesh.nvtxs;i++){
      mesh.vertices[i].code=-1;
      }
    for (int i=0;i<mesh.nedges;i++) {
      if(mesh.edges[i].code==5) {
        n++;
        //printf("open boundary N°%d  1 :   %d\n",n,mesh.edges[i].extremity[0]);
        //printf("open boundary N°%d  2 :   %d\n",n,mesh.edges[i].extremity[1]);
        int n1=mesh.edges[i].extremity[0];
        int n2=mesh.edges[i].extremity[1];
        if (depth_conv * mesh.vertices[n1].h <= 0) {
          mesh.vertices[n1].code=5;
          }
        else {
          printf("warning, boundary node IS DRY ! : vertex number  = %d   depth = %lf (depth factor = %lf)\n", n1, mesh.vertices[n1].h, depth_conv);
          if(allow_dry) mesh.vertices[n1].code=5;
          else printf("boundary node %d discarded\n", n1);
          }
        if (depth_conv * mesh.vertices[n2].h <= 0) {
          mesh.vertices[n2].code=5;
          }
        else {
          printf("warning, boundary node IS DRY ! : vertex number  = %d   depth = %lf (depth factor = %lf)\n", n2, mesh.vertices[n2].h, depth_conv);
          if(allow_dry) mesh.vertices[n2].code=5;
          else printf("boundary node %d discarded\n", n2);
          }
        }
      }
    }

  if(reduce==1) status= fe_reducebw(mesh, 100);
  
  status=fe_limits2poly(mesh, polygons, 0, true);
  status=plg_save("extracted.plg", PLG_FORMAT_SCAN, polygons);


  printf("#################################################################\n");
  printf("save mesh in %s\n",output);
  switch(mode) {
    case 0:
      status=fe_savemesh(output, fmt_out, mesh, oproj4);
      break;

    case 1:
      status=fe_updatemeshNC3D(output, mesh,  0);
      break;
    }

  STDOUT_BASE_LINE("end of mesh-format ... \n");
  exit(0);
}
