
/*******************************************************************************

  T-UGO tools, 2006-2019

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Generates a 2D mesh according to the given criteria file.

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

#include "tools-structures.h"

#include "functions.h"
#include "geo.h"
#include "fe.h"
#include "archive.h"
#include "map.h"
#include "topo.h"
#include "polygones.h"
#include "grd.h"
#include "list.h"
#include "maths.h"
#include "map.def"
#include "mass.h"
#include "netcdf-proto.h"
#include "xyz.h"
#include "statistic.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int addnodes(mesh_t & mesh, point_t *points, int npoints, double dmin, projPJ projection, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m, n, status;
  int ninteriors, ndata;
  int verbose=0;
  double *x=0,*y=0,*z=0;
  char *keep=0;
  vector<plg_t> polygons;
  mesh_t tmp;
  double xx,yy;
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  register interior nodes from mesh and add extra nodes
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ninteriors=0;
  for(n=0; n<mesh.nvtxs; n++) {
    if(mesh.vertices[n].code==0) ninteriors++;
    }
  
  ndata=mesh.nvtxs+npoints;
  
  x=new double[ndata];
  y=new double[ndata];
  z=new double[ndata];
  keep=new char[ndata];
  
  m=0;
  for(n=0; n<mesh.nvtxs; n++) {
    if(mesh.vertices[n].code==0) {
      x[m]=mesh.vertices[n].lon;
      y[m]=mesh.vertices[n].lat;
      m++;
      }
    }
  
  for(n=0; n<npoints; n++) {
    x[m]=points[n].t;
    y[m]=points[n].p;
    m++;
    }
  
  for (m=0;m<ndata;m++) {
    keep[m]=1;
    }
  
  printf("#################################################################\n");
  printf("check for duplicates dmin = %lf\n",dmin);
  int bkp=ndata;
  
  status=xyz_CheckDuplicate_cartesian(x, y, z, ndata, keep, dmin);
  status=xyz_reduce(x, y, z, keep, ndata);
  printf("initial count %d, after reduction %d\n", bkp, ndata);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  import boundaries from mesh
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_limits2poly(mesh,polygons,0,debug);
  mesh.destroy();

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  assembly node set
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_createnodes(polygons, x, y, 0, ndata, tmp);
  status=fe_savenodes("tmp.nod", NODE_FILE_FORMAT_TRIGRID, tmp);
  
  if(debug) {
    for (n=0;n<tmp.nvtxs;n++) {
      xx=tmp.vertices[n].lon;
      yy=tmp.vertices[n].lat;
      projection_to_geo(projection,&(tmp.vertices[n].lat),&(tmp.vertices[n].lon),xx,yy);
      }
    status=fe_savenodes("tmp.nod", NODE_FILE_FORMAT_TRIGRID, tmp);
    for (n=0;n<tmp.nvtxs;n++) {
      xx=tmp.vertices[n].lon;
      yy=tmp.vertices[n].lat;
      geo_to_projection(projection,yy,xx,&(tmp.vertices[n].lon),&(tmp.vertices[n].lat));
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create triangle elements and finalize mesh
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_triangulate(tmp, projection, mesh, verbose, debug);
  tmp.destroy();

  status=fe_edgetable(&mesh,0,0);
  status=fe_vertex_crosstables02(&mesh);
  status=fe_codetable2(&mesh,0,1,1);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int import(char *filename, char *proj4_options, plg_t *polygons, int npolygons, mesh_t & mesh, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  mesh_t tmp;
  int format;
  int m,n,ndata,bkp;
  double dmask;
  double *x=0,*y=0,*z=0;
  char *keep=0;
  range_t<double> lon_range,lat_range;
  string output=filename;
  frame_t plgframe;
  string header;
  projPJ *projection=0;
  bool duplicate=true;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  load random data
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  char *s=strstr(filename,".shp");
  
  if(s==NULL) {
    format=ASCII;
    }
  else{
    format=SHAPE;
    }
  
  switch (format) {
    case ASCII:
      status=xyz_loadraw (filename, header, proj4_options, x, y, z, &dmask, ndata, debug);
      break;
    case SHAPE:
      status=xyz_load_shp (filename, proj4_options, x, y, z, &dmask, ndata);
      break;
    }
  
  if(status!=0) TRAP_ERR_EXIT(status," error by loading %s\n",filename);
  
  if(polygons!=0) {
    printf("#################################################################\n");
    printf("interior nodeset reduction by mesh boundary polygons\n");
    bkp=ndata;
    status=xyz_PolygonSelection(polygons, npolygons, x, y, z, ndata, plgframe);
    printf("initial count %d, after reduction %d\n", bkp, ndata);
    }
  
/*------------------------------------------------------------------------------
  check duplicates */
  if(duplicate) {
    double dmin=1.0;
    printf("#################################################################\n");
    printf("check for duplicates dmin = %lf m\n",dmin);
    bkp=ndata;
    keep=new char[ndata];
    for (m=0;m<ndata;m++) {
      keep[m]=1;
      }
    status=xyz_CheckDuplicate_spherical(x, y, z, ndata, keep, dmin);
    status=xyz_reduce(x, y, z, keep, ndata);
    printf("initial count %d, after reduction %d\n", bkp, ndata);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  assembly node set
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_createnodes(polygons, npolygons, x, y, ndata,tmp);
  status=fe_savenodes("tmp.nod", NODE_FILE_FORMAT_TRIGRID, tmp);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create triangle elements and finalize mesh
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_triangulate(tmp, projection, mesh, 0, debug);
  tmp.destroy();

  status=fe_edgetable(&mesh,0,0);
  status=fe_vertex_crosstables02(&mesh);
  status=fe_codetable2(&mesh,0,1,1);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  back to spherical coordinates
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if( projection!=0) {
    for (n=0;n<mesh.nvtxs;n++) {
      double xx=mesh.vertices[n].lon;
      double yy=mesh.vertices[n].lat;
      projection_to_geo(projection,&(mesh.vertices[n].lat),&(mesh.vertices[n].lon),xx,yy);
      }
    }
  
  mesh.type=0;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int generate_nodes(criteria_t criteria, grid_t sgrid, float *density, float mask, plg_t *polygons, int npolygons, point_t* & points, int & npoints, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k;
  int ii,jj,kk;
  int count=0;
  char *interior;
  grid_t *grid;
  
  grid=&sgrid;
  
/*------------------------------------------------------------------------------
  flag interior cells*/
  grid->projection=sgrid.projection;
  interior=plg_test_grid(*grid,polygons,npolygons);

  count=0;
  for(j=0;j<grid->ny;j++) {
    for(i=0;i<grid->nx;i++) {
      k=grid->nx*j+i;
      if(interior[k]==1) count++;
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  eliminate some unnecessary points
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  count=0;
  for(j=2;j<grid->ny-2;j+=2) {
    for(i=2;i<grid->nx-2;i+=2) {
      k=grid->nx*j+i;
      if(density[k]==mask) continue;
      if(interior[k]!=1)   continue;
      if(density[k]>3.*grid->dx) {
        for(jj=j-2;jj<j+3;jj++) {
          for(ii=i-2;ii<i+3;ii++) {
            kk=grid->nx*jj+ii;
//            if(interior[kk]!=1)   goto next;
            if(density[kk]==mask) goto next;
            if(density[kk]<3.*grid->dx) {
              goto next;
              }
            }
          }
        for(jj=j-1;jj<j+2;jj++) {
          for(ii=i-1;ii<i+2;ii++) {
            kk=grid->nx*jj+ii;
            if(kk!=k) interior[kk]=-1;
            }
          }
next:
        count++;
        }
      }
    }

  count=0;
  for(j=0;j<grid->ny;j++) {
    for(i=0;i<grid->nx;i++) {
      k=grid->nx*j+i;
      if(interior[k]==1) count++;
      }
    }

  points=new point_t[count];

  count=0;
  for(j=0;j<grid->ny;j++) {
    for(i=0;i<grid->nx;i++) {
      k=grid->nx*j+i;
      if(interior[k]==1) {
        points[count].t=grid->x[k];
        points[count].p=grid->y[k];
        count++;
        }
      }
    }

  npoints=count;
  
  delete[] interior;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_adjust(criteria_t criteria, grid_t sgrid, float *density, float mask, plg_t *polygons, int npolygons, bool reshape, mesh_t & spherical, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,step,status;
  double x,y;
  grid_t *grid;
  mesh_t previous;
//   const double km2m=1e3;
//   double Lmin=min(criteria.minsize,criteria.shelf_minsize)*km2m;
  int previous_nvtxs=0;
  
  grid=&sgrid;
  
  if(verbose) printf(
    "#################################################################\n"
    "decimate node set according to mesh resolution\n\n");
  
  for (step=0;step<criteria.niterations;step++) {
/*------------------------------------------------------------------------------
    eliminate nodes*/
    if(verbose) printf("step %3d: nvtxs=%6d\n", step,spherical.nvtxs);

    spherical.type=1;

    status=fe_geometry(&spherical);

    status=fe_edgetable(&spherical,0,0);
    if(status!=0) {
      for (n=0;n<spherical.nvtxs;n++) {
        x=spherical.vertices[n].lon;
        y=spherical.vertices[n].lat;
        projection_to_geo(sgrid.projection,&(spherical.vertices[n].lat),&(spherical.vertices[n].lon),x,y);
        }
      status=fe_savemesh("mesh-abort.nei",MESH_FILE_FORMAT_TRIGRID, spherical);
      STDOUT_BASE_LINE("abandon mesh generation, current mesh in %s\n", "mesh-abort.nei");
      exit(-1);
      }
    status=fe_vertex_crosstables02(&spherical);
    status=fe_codetable2(&spherical,0,1,1);
    if(status!=0) {
      for (n=0;n<spherical.nvtxs;n++) {
        x=spherical.vertices[n].lon;
        y=spherical.vertices[n].lat;
        projection_to_geo(sgrid.projection,&(spherical.vertices[n].lat),&(spherical.vertices[n].lon),x,y);
        }
      status=fe_savemesh("mesh-abort.nei",MESH_FILE_FORMAT_TRIGRID, spherical);
      STDOUT_BASE_LINE("abandon mesh generation, current mesh in %s\n", "mesh-abort.nei");
      exit(-1);
      }

    if(reshape) status=fe_reshapeall(spherical,3);

    if(debug) {
      for (n=0;n<spherical.nvtxs;n++) {
        x=spherical.vertices[n].lon;
        y=spherical.vertices[n].lat;
        projection_to_geo(sgrid.projection,&(spherical.vertices[n].lat),&(spherical.vertices[n].lon),x,y);
        }
      char filename[1024];
      sprintf(filename,"tmp.%2.2d.nei",step);
      status=fe_savemesh(filename,MESH_FILE_FORMAT_TRIGRID, spherical);
      for (n=0;n<spherical.nvtxs;n++) {
        x=spherical.vertices[n].lon;
        y=spherical.vertices[n].lat;
        geo_to_projection(sgrid.projection,y,x,&(spherical.vertices[n].lon),&(spherical.vertices[n].lat));
        }
      }

/**-----------------------------------------------------------------------------
    to be optimized */
    status=fe_DecimateInterior(criteria, spherical, sgrid, density, mask, polygons, npolygons);
    
    status=fe_savenodes("tmp.nod", NODE_FILE_FORMAT_TRIGRID, spherical);
    spherical.destroy();
    
    status=fe_readnodes("tmp.nod", NODE_FILE_FORMAT_TRIGRID, &previous);

    status=system("rm -f tmp.nod");
      
    status=fe_triangulate(previous, 0, spherical, verbose, debug);
    previous.destroy();
    
    if(spherical.nvtxs==previous_nvtxs) {
      break;
      }
    previous_nvtxs=spherical.nvtxs;
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_nodit_new(criteria_t criteria, grid_t sgrid, float *density, float mask, plg_t *polygons, int npolygons, mesh_t & spherical, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  generate finite element mesh by creating interior nodes inside prescribed
  boundaries

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int n,status;
  int npoints;
  point_t *points;
  mesh_t *mesh;
  double x,y;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  generate prior interior nodes set

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=generate_nodes(criteria, sgrid, density, mask, polygons, npolygons, points, npoints, verbose, debug);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  add boundary nodes to interior nodes

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  mesh=fe_createnodes(polygons, npolygons, points, npoints);
  delete[] points;
  
  if(debug) {
    for (n=0;n<mesh->nvtxs;n++) {
      x=mesh->vertices[n].lon;
      y=mesh->vertices[n].lat;
      projection_to_geo(sgrid.projection,&(mesh->vertices[n].lat),&(mesh->vertices[n].lon),x,y);
      }
    status=fe_savenodes("tmp.nod", NODE_FILE_FORMAT_TRIGRID, *mesh);
    for (n=0;n<mesh->nvtxs;n++) {
      x=mesh->vertices[n].lon;
      y=mesh->vertices[n].lat;
      geo_to_projection(sgrid.projection,y,x,&(mesh->vertices[n].lon),&(mesh->vertices[n].lat));
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  generate initial mesh

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_triangulate(*mesh, sgrid.projection, spherical, 0, debug);
  mesh->destroy();

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  adjust (decimate) mesh interior nodes set to fit nodes density criteria

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_adjust(criteria, sgrid, density, mask, polygons, npolygons, true, spherical, verbose, debug);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
 finalize mesh

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_edgetable(&spherical,0,0);
  status=fe_vertex_crosstables02(&spherical);
  status=fe_codetable2(&spherical,0,1,1);
  
/*------------------------------------------------------------------------------
  finalize mesh*/
  for (n=0;n<spherical.nvtxs;n++) {
    x=spherical.vertices[n].lon;
    y=spherical.vertices[n].lat;
    projection_to_geo(sgrid.projection,&(spherical.vertices[n].lat),&(spherical.vertices[n].lon),x,y);
    }
  
  spherical.type=0;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_complete(const char *output, char *proj4_options, plg_t *polygons, int npolygons, criteria_t & criteria, mesh_t & mesh, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, status;
  int verbose=0;
//   criteria_t criteria;
  grid_t grid;
  float *density=0;
  float mask;
  int npoints;
  point_t *points;
  mesh_t work;
  double x,y;
  
//   double resolution=5.0;
  double dmin=0.9*criteria.minsize;
    
  status=initialise_meshgrid(output, criteria, &grid, &density, &mask, polygons, npolygons);
   
  for(size_t m=0;m<grid.Hsize();m++) density[m]=criteria.minsize*1000.;

//   status=fe_nodit_new(criteria, grid, density, mask, polygons, npolygons, mesh, 0, false);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  generate prior interior nodes set

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(mesh.type==0)
    for (n=0;n<mesh.nvtxs;n++) {
      x=mesh.vertices[n].lon;
      y=mesh.vertices[n].lat;
      geo_to_projection(grid.projection,y,x,&(mesh.vertices[n].lon),&(mesh.vertices[n].lat));
      }

  mesh.type=1;
      
  status=generate_nodes(criteria, grid, density, mask, polygons, npolygons, points, npoints, verbose, debug);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  add nodes to mesh

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=addnodes(mesh, points, npoints, dmin*1000., grid.projection, debug);
//   mesh.destroy();
// 
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//   
//   generate initial mesh
// 
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// 
//   status=fe_triangulate(work, grid.projection, mesh, 0, debug);
//   work.destroy();

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  adjust (decimate) mesh interior nodes set to fit nodes density criteria

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_adjust(criteria, grid, density, mask, polygons, npolygons, false, mesh, verbose, debug);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
 finalize mesh

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_edgetable(&mesh,0,0);
  status=fe_vertex_crosstables02(&mesh);
  status=fe_codetable2(&mesh,0,1,1);
  
/*------------------------------------------------------------------------------
  finalize mesh*/
  for (n=0;n<mesh.nvtxs;n++) {
    x=mesh.vertices[n].lon;
    y=mesh.vertices[n].lat;
    projection_to_geo(grid.projection,&(mesh.vertices[n].lat),&(mesh.vertices[n].lon),x,y);
    }
  
  mesh.type=0;
  
  return(0);
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
    "  %s [OPTIONS] criteria.crt\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Generates a 2D mesh according to the given criteria file.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  Show this help and exit.\n"
    "  --limits-only  produce everything (density and limits) but the mesh\n"
    "  --rootname  root name of the MANY intermediate output files. Default: anonymous\n"
    "  --debug  produce extra output files for debugging\n"
    "  -b  followed by path to bathymetry file. Can also be given in TOPOGRAPHY_FILE entry of criteria file.\n"
    "  -f  followed by path to mesh boundaries file. Disables -d and -p\n"
    "  -d  followed by path to genesis boundary descriptor. Disables -p\n"
    "  -p  followed by path to raw polygons for boundaries\n"
    "  -m  followed by path to a mesh to merge with\n"
    "  -o  followed by path to the output mesh. Default: <ROOTNAME>-mesh.nei\n"
    "  -s  followed by path to mesh density file\n"
    "  -z  followed by path to local mesh density prescription\n"
    );
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  char *keyword,*zone=NULL;
  char *input=NULL,*meshfile=NULL,*nodefile=0;
  char *depthfile=NULL,*output=NULL,*poly=NULL,*meshsize=NULL;
  char *boundaryfile=NULL,*descriptor=NULL;
  char neighfile[1024],additional[1024];
  mesh_t mesh;
  mesh_t external,splitted[2],final;
  grid_t grid;
  criteria_t criteria;
  plg_t *polygons;
  plg_array_t set;
  int format=PLG_FORMAT_BOUNDARIES,npolygons;
  float *density,mask;
  double x,y;
  int unsafe=0, mesh_build=1;
  bool debug=false, show_map=false, reshape=true;
  string rootname;
  string echofile;
  double resolution=5000.0;
  char *proj4_options=0;
  grid_t topogrid;
  float *topobase, topomask;

  plg_t test,target;

  if(argc<=1){
    print_help(argv[0]);
    exit(0);
    }
  
  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        if(strcmp(keyword,"--limits-only")==0) {
          mesh_build=0;
          n++;
          break;
          }
        if(strcmp(keyword,"--debug")==0) {
          debug=true;
          n++;
          break;
          }
        if(strcmp(keyword,"--no-reshape")==0) {
          reshape=false;
          n++;
          break;
          }
        if(strcmp(keyword,"--show_map")==0) {
          show_map=true;
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

        case 'b' :
          depthfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'd' :
          descriptor= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'f' :
          boundaryfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'n' :
          nodefile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          poly= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'r' :
          sscanf(argv[n+1],"%lf",&resolution);
          n++;
          n++;
          break;

        case 's' :
          meshsize= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        input= strdup(argv[n]);
        n++;
        break;
      }
    free(keyword);
    }
  
  if(rootname=="") rootname="anonymous";

  if(not show_map) {
    printf("\nWARNING : netcdf dump of mesh density not activated (as it used to be by default); please use --show_map to activate it)\n");
    }

  if(output==0) {
    string tmp=rootname+"-mesh.nei";
    output=strdup(tmp.c_str());
    printf("\nWARNING : no filename specified for output (-o [filename]); using default name (%s)...\n",output);
    }

  if(input!=0){
/*------------------------------------------------------------------------------
    load and initialize criteria structure*/
    printf("#################################################################\n");
    printf("\nload resolution criteria file : %s\n",input);
    status=fe_load_criteria(input, &criteria);
    if(status!=0) TRAP_ERR_EXIT(status,"fe_load_criteria(\"%s\" error,) error (%d: %s)\n",input,status,strerror(status));
    }
  else if(nodefile!=0) {
    criteria.uniform(resolution/1000.0);
    }
  else {
    printf("\n*** FATAL ERROR : no criteria file nor nodes file specified, abort ***\n\n");
    print_help(argv[0]);
    exit(-1);
    }
  
  if(depthfile!=0) {
    printf("#################################################################\n");
    printf("substitute criteria depth file (%s) with %s\n", criteria.regulardepth, depthfile);
    strcpy(criteria.regulardepth,depthfile);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  load mesh boundaries
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  if(boundaryfile!=0) {
/*------------------------------------------------------------------------------
    from mesh boundaries file (trigrid format)*/
    printf("load mesh boundary file : %s\n",boundaryfile);
    status=fe_loadboundaries(boundaryfile, format, &polygons, &npolygons);
    if(status!=0) TRAP_ERR_EXIT(status,"cannot load (or use) boundary file %s (%d: %s)\n",boundaryfile,status,strerror(status));
    }
  else if(descriptor!=0) {
/*------------------------------------------------------------------------------
    from genesis boundary descriptor*/
    printf("load boundary descriptor : %s\n",descriptor);
    status=plg_read_descriptor(descriptor, polygons, npolygons, true, debug);
    if(status!=0) TRAP_ERR_EXIT(status,"cannot load (or use) descriptor file %s (%d: %s)\n",descriptor,status,strerror(status));
    if(debug) {
      string tmp=rootname+"-mesh-boundaries.fr";
      status=plg_write_boundaries (tmp.c_str(), polygons, npolygons);
      }
    }
  else if(poly!=0) {
/*------------------------------------------------------------------------------
    from raw polygons*/
    printf("load mesh boundaries from raw polygons : %s\n",poly);
    status=plg_load_scan(poly, &polygons, &npolygons);
    if(status!=0) TRAP_ERR_EXIT(status,"cannot load (or use) polygons file %s (%d: %s)\n",poly,status,strerror(status));
    frame_t frame=plg_recale(polygons, npolygons);
    status=plg_setdirect(polygons, npolygons);
//     status=plg_write_boundaries ("mesh-generator.fr", polygons, npolygons);
    if(debug) {
      string tmp=rootname+"-mesh-boundaries.fr";
      status=plg_write_boundaries (tmp.c_str(), polygons, npolygons);
      }
    }
  else  {
    printf("*** no boundary information given ***\n");
    print_help(argv[0]);
    exit(-1);
    }
  
  if(debug) {
    string tmp=rootname+"-mesh-boundaries.plg";
    status=plg_save (tmp.c_str(), PLG_FORMAT_SCAN,  polygons, npolygons);
    }
  
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  import/create interior nodes
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(nodefile!=0) {
    status=import(nodefile, proj4_options, polygons, npolygons, mesh, debug);
    status=fe_complete("test.nc", proj4_options, polygons, npolygons, criteria, mesh, debug);
    reshape=false;
    }
  else {
    
/*------------------------------------------------------------------------------
  mesh size mapping*/
  criteria.resample_obsolete=( (criteria.resample_openlimits==1) && (criteria.resample_rigidlimits==1) );

  if(meshsize==0) {
    string tmp;
    printf("#################################################################\n");
    printf("compute mesh density\n");
    if(show_map) tmp=rootname+"-criteria.nc";
    else tmp="";
    status=fe_ComputeMeshsize(tmp.c_str(), criteria, &grid, &density, &mask, polygons, npolygons, 1, false);
    }
  else {
    printf("#################################################################\n");
    printf("load mesh density from %s\n",meshsize);
    status=defproj(&grid, polygons, npolygons);
    status=fe_ReloadMeshsize(meshsize, &grid, &density, &mask, polygons, npolygons);
    }

  printf("#################################################################\n");
  printf("check for isolated boundaries\n");
  unsafe=0;
  for(int s=1; s<  npolygons; s++) {
    char *position=plg_TestInterior(polygons[s].x,polygons[s].y,polygons[s].npt, polygons, 1, 0, debug);
    for(int k=0;k<polygons[s].npt;k++) {
      if(position[k]!=PLG_POINT_INTERIOR) {
        polygons[s].npt=0;
        unsafe++;
        printf("\n !!! boundary polygon %d found to be outside the domain limits (later discarded) \n",s);
        printf(" !!! Remember: when giving polygons, first one must be the main boundary \n\n\a",s);
        break;
        }
      }
    delete[] position;
    }
  

  if(criteria.resample_openlimits) {
    string tmp;
    printf("#################################################################\n");
    printf("open boundary adjustment\n");
/*-----------------------------------------------------------------------------------
    next step necessary because of projection/spherical issue */
    for(int k=0;k<20;k++) status=map_persistence(grid, density, mask, 0.);
    if(debug) {
      tmp=rootname+"-criteria-extended.nc";
      status=save_SG(tmp.c_str(), grid, density, mask, "density", "m", "density", 1);
      }
/*-----------------------------------------------------------------------------------
    re-sample polygons by following meshing criteria */
    for(int s=0; s<  npolygons; s++) {
      if(polygons[s].flag!=0) {
        plg_t q=plg_sample(polygons[s], grid, density, mask, criteria, 1, debug);
        polygons[s].duplicate(q);
        }
      }
    status=plg_spherical(grid.projection, polygons, npolygons);
/*------------------------------------------------------------------------------
    might be inaccurate for polar grids, to be thought again */
    status=plg_checkAutoSecant(polygons, npolygons, PLG_SPHERICAL);
    tmp=rootname+"-open-adjusted.plg";
    status=plg_save (tmp.c_str(), PLG_FORMAT_SCAN,  polygons, npolygons);
    grid.free();
    delete[] density;
    density=0;
    criteria.resample_obsolete=0;
    criteria.resample_openlimits=0;
    if(show_map) tmp=rootname+"-criteria-adjusted.nc";
    else tmp="";
/*------------------------------------------------------------------------------
    no needs to redo everythings, to be fixed */
    status=fe_ComputeMeshsize(tmp.c_str(), criteria, &grid, &density, &mask, polygons, npolygons, 1, false);
/*-----------------------------------------------------------------------------------
    next step necessary because of projection/spherical issue */
    for(int k=0;k<20;k++) status=map_persistence(grid, density, mask, 0.);
    }

  if(criteria.resample_rigidlimits) {
    string tmp;
    printf("#################################################################\n");
    printf("rigid boundary adjustment\n");
/*-----------------------------------------------------------------------------------
    re-sample polygons by following meshing criteria */
    for(int s=0; s<  npolygons; s++) {
      if(polygons[s].npt==0) {
        printf("empty polygon %d\n",s);
        continue;
        }
      if(polygons[s].flag!=0) {
        plg_t q=plg_resample_rigid_obsolete(polygons[s], grid.projection, criteria.shelf_minsize, 1);
        polygons[s].duplicate(q);
        }
      }
    status=plg_spherical(grid.projection, polygons, npolygons);
    status=plg_checkAutoSecant(polygons, npolygons, PLG_SPHERICAL);
    tmp=rootname+"-rigid-adjusted.plg";
    status=plg_save(tmp.c_str(), PLG_FORMAT_SCAN,  polygons, npolygons);
    status=plg_load_scan(tmp.c_str(), &polygons, &npolygons);
    grid.free();
    delete[] density;
    density=0;
    criteria.resample_obsolete=0;
    criteria.resample_rigidlimits=0;
    if(show_map) tmp=rootname+"-criteria-adjusted.nc";
    else tmp="";
/*------------------------------------------------------------------------------
    no needs to redo everythings, to be fixed */
    status=fe_ComputeMeshsize(tmp.c_str(), criteria, &grid, &density, &mask, polygons, npolygons, 1, false);
/*-----------------------------------------------------------------------------------
    next step necessary because of projection/spherical issue */
    for(int k=0;k<20;k++) status=map_persistence(grid, density, mask, 0.);
    }

/**-----------------------------------------------------------------------------
  echo criteria file */
  echofile=rootname+"-echo.crt";
  status= fe_save_criteria(echofile.c_str(), criteria);
  
/*-------------------------------------------------------------------------------------
  local tuning using extra criteria prescription file*/
  if(zone!=0) {
    string tmp;
    printf("#################################################################\n");
    printf("load local mesh density prescription %s\n",zone);
    criteria_t local;
    FILE *in;
    in=fopen(zone,"r");
    fscanf(in,"%s",neighfile);
    string s=(string) neighfile;
    set=plg_readneigh (s);
    int nitems=3;
    while(nitems==3) {
      nitems=fscanf(in,"%lf %lf %s",&x,&y,additional);
      if(nitems!=3) break;
      status=fe_load_criteria(additional, &local);
      if(show_map) tmp=rootname+"-criteria-local.nc";
      else tmp="";
      status=fe_ModifyMeshsize(tmp.c_str(), local, grid, density, mask, set.p, set.n);
      }
    fclose(in);
    }

  if(mesh_build==0) goto end;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  mesh generation
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("mesh generation\n");
  mesh=fe_nodit(criteria, grid, density, mask, polygons, npolygons, false);
  
  if(debug) {
    string tmp=rootname+"-mesh-no-reshape.nei";
    status= fe_savemesh(tmp.c_str(),MESH_FILE_FORMAT_TRIGRID, mesh);
    }
  }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  optionally merge with another mesh
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(meshfile!=0) {
    printf("#################################################################\n");
    printf("save generated mesh : %s\n", output);
    string raw=rootname+"-generated.nei";
    status= fe_savemesh(raw.c_str(),MESH_FILE_FORMAT_TRIGRID, mesh);
    printf("#################################################################\n");
    printf("merge meshes : adding %s\n", meshfile);
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&external);
    if(status!=0) TRAP_ERR_EXIT(status,"unable to read external mesh from %s (%d: %s)\n",meshfile,status,strerror(status));
    status=fe_list(&external);
    if(status!=0) TRAP_ERR_EXIT(status,"unable to build the element list from the original mesh\n");
    status= fe_edgetable(&external,0,0);
    status= fe_codetable2(&external,0,1,1);
    if(status!=0) TRAP_ERR_EXIT(status,"unable to rebuild the limits table and codes of the original mesh\n");
    splitted[0]=external;
    splitted[1]=mesh;
    final=fe_merge(splitted,0.1, (char *) 0);
    mesh=final;
    }
  
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  optionally reshape mesh
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(reshape) {
/**-----------------------------------------------------------------------------
    reshape mesh */
    printf("#################################################################\n");
    printf("reshape mesh\n");
    status= fe_reshapeall(mesh,3);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  interpolate depth
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("interpolate topography on mesh from %s\n", criteria.regulardepth);
  
  status=topo_loadfield(criteria.regulardepth, &topogrid, &topobase, &topomask, debug);
  if(status!=0){
    printf("WARNING not a workable topography !!! Setting to default %g m.\n",criteria.hmax);
    
    topogrid=get_zonegrid("30.");
    topogrid.print(cout,"topo:  ");
    
    topobase=aset(topogrid.Hsize(),(float)criteria.hmax);
    topomask=-1.;
    }
#pragma omp parallel for private(status)
  for (n=0;n<mesh.nvtxs;n++) {
    double x=mesh.vertices[n].lon;
    double y=mesh.vertices[n].lat;
    float z;
    status=map_interpolation(topogrid, topobase,topomask,x,y,&z);
    if (z!=topomask)
      mesh.vertices[n].h=-z;
    else {
      mesh.vertices[n].h=-20.;
      }
    }
  



/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  save final mesh
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("save final mesh : %s\n", output);
  status= fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID, mesh);
  
end:
  printf("mesh-generator sucessfully completed\n");

  STDOUT_BASE_LINE("end of mesh-generator ... \n");
  exit(0);
}
