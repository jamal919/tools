
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
*/
/*----------------------------------------------------------------------------*/

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
#include "xyz.h"
#include "geo.h"
#include "grd.h"
#include "topo.h"
#include "archive.h"
#include "netcdf-proto.h"
#include "poc-time.h"
#include "functions.h"
#include "poc-netcdf-data.hpp"

#include "polygons.hpp"
#include "exceptions.hpp"

int recursion_maxcount=2000;

class metafield_t {
private:
public:
  string variable;
  string format;
  string gridfile, maskfile, datafile;
  string type;                           // structured/unstructured           
  string frame;
  float tag;
  
  void destroy() {
    variable.clear();
    format.clear();
    gridfile.clear();
    maskfile.clear();
    datafile.clear();
    type.clear();
    frame.clear();
   }
};

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int parse_FieldOptions(const string & options, metafield_t * meta, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
  parse field input setting
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  vector<string> tokens,keys,values;
  int k;
  string delimiter=" ";
      
/*------------------------------------------------------------------------------
  first separate space-separated tokens */
  delimiter=" ";
  tokens=string_split(options, delimiter);
  
  if(debug) printf("tokens founds: %d\n",tokens.size());
  
  if(tokens.size()==1) {
    meta->datafile=options;
    return(0);
    }
  
/*------------------------------------------------------------------------------
  then identify key=value paires */
  delimiter="=";
  for(k=0;k<tokens.size();k++) {
    vector<string> tmp=string_split(tokens[k], delimiter);
    if(debug) printf("token %d, paire sizes: %d\n",k,tmp.size());
    keys.push_back(tmp[0]);
    if(tmp.size()==1) {
/*------------------------------------------------------------------------------
      old-fashioned input, no more supported */
      printf("wrong input setting %s; please set <file=%s format=NETCDF variable=xxx etc...\n",options.c_str(),options.c_str());
      TRAP_ERR_EXIT(-1, "parsing of %s failed\n", options.c_str());
      }
    values.push_back(tmp[1]);
    }
  
  for(k=0;k<keys.size();k++) {
    const string *keyk=&keys[k];
    char *valuek=poc_strdup(values[k].c_str());
    
    if(*keyk=="file") {
      meta->datafile=valuek;
      delete[]valuek;
      continue;
      }
    if(*keyk=="grid") {
      meta->gridfile=valuek;
      delete[]valuek;
      continue;
      }
    if(*keyk=="mask") {
      meta->maskfile=valuek;
      delete[]valuek;
      continue;
      }
    if(*keyk=="format") {
      meta->format=valuek;
      delete[]valuek;
      continue;
      }
    if(*keyk=="variable") {
      meta->variable=valuek;
      delete[]valuek;
      continue;
      }
    if(*keyk=="type") {
      meta->type=valuek;
      delete[]valuek;
      continue;
      }
    if(*keyk=="frame") {
      meta->frame=valuek;
      delete[]valuek;
      continue;
      }
    
    STDERR_BASE_LINE_FUNC("key \""+*keyk+"\" not recognised\n");
    delete[]valuek;
    }
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int store_topo(const char *rootname, mesh_t & mesh, frame_t & mapping, double grid_resolution, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  grid_t TOPOgrid;
  float *TOPOdata=0, TOPOmask=1.e+10;
  char output[1024];
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  create bathymetry structured grid 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("create bathymetry on structured grid: dx=%lf (%lf m) arc-sec dy=%lf (%lf m) arc-sec\n",
         grid_resolution*3600., grid_resolution*110000., grid_resolution*3600., grid_resolution*110000.);
    
  double dx=grid_resolution;
  double dy=grid_resolution;
  
  TOPOgrid=map_Rgrid(mapping, dx, dy, 0);
  map_printgrid(TOPOgrid);
  
  TOPOdata=new float[TOPOgrid.Hsize()];
  TOPOmask=1.e+10;
    
  int algorithm=0, verbose=0;
  
  status=fe_mapVdepth(mesh, TOPOgrid, TOPOdata, TOPOmask, algorithm, verbose);
  
//   status=checks(rootname,TOPOgrid,TOPOdata,TOPOmask,false);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  save data, to be replaced with topo_save. Example:
  
  status=topo_save(rootname, output, "netcdf", TOPOgrid, TOPOdata, TOPOmask, debug);
  
  topo_save will hande mirroring...

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=grd_mirror_r( TOPOgrid, TOPOgrid.nx, TOPOdata, TOPOmask);
  sprintf(output, "%s.grd",rootname);
  status=grd_save(output,TOPOgrid, TOPOgrid.nx, TOPOdata, TOPOmask);

  STDOUT_BASE_LINE("Finished with %s\n",__func__);
  
  TOPOgrid.free();
  delete[] TOPOdata;
  
  return(0);
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_delaunay_improve(const char *rootname, mesh_t & mesh, double ElementMaxSize, double ElementMaxRatio, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  post-processing (mostly for node's cloud stand-alone mesh)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  int k,m,n;
  size_t count;
  double d;
  char output[1024];
  

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  filter elements by size (to clean Delaunay mesh), pretty unsafe
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  int filter=(ElementMaxSize>0);
  if(filter==1) {
    printf("#################################################################\n");
    printf("remove oversized element connexion: maxsize=%lf m\n",ElementMaxSize);
    count=0;
    for(n=0;n<mesh.nvtxs;n++) {
      for(k=0;k<mesh.vertices[n].nngh;k++) {
        m=mesh.vertices[n].ngh[k];
        if(m<n) continue;
/*------------------------------------------------------------------------------
        warning : fe_distance now return distance in meters */
        d=fe_distance(mesh, n, m);
        if(d>ElementMaxSize) {
          status= fe_disconnectvertices(mesh, n,m);
          k=k-1;
          count++;
          }
        }
      }
    printf("remove oversized edges, discarded connections=%d\n", count);
    status=fe_cleanvertices(&mesh, false);
    if(count!=0) {
      if(mesh.triangles!=0) mesh.triangles->destroy();
      if(mesh.edges!=0)     mesh.edges->destroy();
      mesh.triangles=0;
      mesh.edges=0;
      }
    if(mesh.triangles==0) {
      status=fe_list(&mesh);
      status=fe_e2n(&mesh);
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  filter elements by size ratio (to clean Delaunay mesh), pretty unsafe
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  filter=(ElementMaxRatio!=0);
  if(filter==1) {
    printf("#################################################################\n");
    printf("remove oversized element connexion: maxsize=%lf m\n",ElementMaxSize);
    int nboundaries, verbose=0;
    count=0;
    for(n=0;n<mesh.nvtxs;n++) {
      double dmin=+1.e+10;
      double dmax=-1.e+10;
      vector<double> lengths;
      if(mesh.vertices[n].nngh<5) continue;
      for(k=0;k<mesh.vertices[n].nngh;k++) {
        m=mesh.vertices[n].ngh[k];
        d=fe_distance(mesh, n, m);
        lengths.push_back(d);
        dmin=min(dmin,d);
        dmax=max(dmax,d);
        }
      for(k=0;k<mesh.vertices[n].nngh;k++) {
        m=mesh.vertices[n].ngh[k];
        if(m<n) continue;
        d=fe_distance(mesh, n, m);
/*------------------------------------------------------------------------------
        following check should be limited to boundary edges */
        if(d>ElementMaxRatio*dmin) {
          status= fe_disconnectvertices(mesh, n,m);
          k=k-1;
          count++;
          }
        }
      }
    printf("remove unbalanced elements, discarded connections=%d\n", count);
    status=fe_cleanvertices(&mesh, false);
    if(count!=0) {
      if(mesh.triangles!=0) mesh.triangles->destroy();
      if(mesh.edges!=0)     mesh.edges->destroy();
      mesh.triangles=0;
      mesh.edges=0;
      }
    if(mesh.triangles==0) {
      status=fe_list(&mesh);
      status=fe_e2n(&mesh);
      }
    }

  printf("#################################################################\n");
  printf("construct face tables\n");
  
  mesh.nedges=0;
  for(n=0;n<mesh.nvtxs;n++) {
    for(k=0;k<mesh.vertices[n].nngh;k++) {
      if(n<mesh.vertices[n].ngh[k]) mesh.nedges++;
      }
    }
   
  mesh.faces=new face_t[mesh.nedges];
  mesh.nedges=0;
  for(n=0;n<mesh.nvtxs;n++) {
    for(k=0;k<mesh.vertices[n].nngh;k++) {
      if(n<mesh.vertices[n].ngh[k]) {
        m=mesh.nedges;
        mesh.faces[m].extremity[0]=n;
        mesh.faces[m].extremity[1]=mesh.vertices[n].ngh[k];
        mesh.nedges++;
        }
      }
    }
   
  if(debug) {
    sprintf(output, "%s-02.nei",rootname);
    status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID,mesh);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  improve vertex connection to favour isobaths
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("improve node connections\n");
  count=0;
  for(n=0;n<mesh.nedges;n++) {
    int n1,n2,q[4];
    double dh[2],angle;
    n1=mesh.faces[n].extremity[0];
    n2=mesh.faces[n].extremity[1];
/*------------------------------------------------------------------------------
    form direct order quadrangle*/
    status=fe_quadnodes(mesh, n1, n2, q);
    if(status!=0) continue;
/*------------------------------------------------------------------------------
    check convexity*/
    for(k=0;k<4;k++) {
      angle=fe_angle(mesh,q[k],q[(k+1)%4],q[(k+3)%4]);
      if(angle<M_PI/6.) {
        status=-1;
        }
      }
    if(status!=0) continue;
/*------------------------------------------------------------------------------
    check aspect ratio*/
    k=1;
    angle=fe_angle(mesh,q[k],q[(k+1)%4],q[(k+2)%4]);
    if(angle<M_PI/3.) continue;
    k=3;
    angle=fe_angle(mesh,q[k],q[(k+1)%4],q[(k+2)%4]);
    if(angle<M_PI/3.) continue;
/*------------------------------------------------------------------------------
    check elevation connexion*/
    dh[0]=fabs(mesh.vertices[q[0]].h-mesh.vertices[q[2]].h);
    dh[1]=fabs(mesh.vertices[q[1]].h-mesh.vertices[q[3]].h);
    if(dh[0]>dh[1]) {
      status=fe_disconnectvertices(mesh, n1, n2);
      if(status!=0) {
        printf("%s : trouble\n",__func__);
        }
      n1=q[1];
      n2=q[3];
      status=fe_connectvertices(mesh, n1, n2);
      if(status!=0) {
        printf("%s : trouble\n",__func__);
        }
      mesh.faces[n].extremity[0]=n1;
      mesh.faces[n].extremity[1]=n2;
      count++;
      }
    }

  delete [] mesh.faces;
  
  if(mesh.edges!=0) mesh.edges->destroy();
  
  status=fe_list(&mesh);
  
  }
  
    
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_delaunay_create(const char *rootname, mesh_t & mesh, mesh_t & delaunay, bool meshout, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int n;
  char output[1024];
  triangulateio in, out;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  initialize triangle output structure and create triangulation

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(debug) status=fe_savenodes("out.node", NODE_FILE_FORMAT_TRIANGLE, mesh);
  
  status= fe_triangulateio_init(in);
  status= fe_triangulateio_init(out);
  
  status= fe_node2triangle(mesh, &in);
  mesh.destroy();
  
  triangulate("z", &in, &out, (triangulateio *) NULL);

/*------------------------------------------------------------------------------
  output structure inherit from input structure's pointlist, to handle with care */
  out.pointlist = in.pointlist;

  printf("#################################################################\n");
  printf("convert triangle mesh into sirocco mesh\n");
  
  status=fe_triangle2mesh(out, &delaunay);
  
/*------------------------------------------------------------------------------
  output structure inherit from input structure's pointlist, to handle with care */
  in.pointlist = out.pointlist;
  fe_FreeTriangleIO(in);
  out.pointlist = 0;
  fe_FreeTriangleIO(out);
  
//   if(debug) {
//     sprintf(output, "%s-01.nei",rootname);
//     status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID,*finished);
//     }

  printf("#################################################################\n");
  printf("construct final topography mesh\n");
  delaunay.nnghm=0;
  for(n=0;n<delaunay.nvtxs;n++) {
    updatemax(&delaunay.nnghm,delaunay.vertices[n].nngh);
    }
  if(delaunay.edges!=0) delaunay.edges->destroy();
  
  status=fe_list(&delaunay);

//   if(meshout) {
//     sprintf(output, "%s.nei",rootname);
//     printf("save topography mesh in %s \n",output);
//     status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID,*finished);
//     }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int load_cloud(char *filename, const char *rootname, char *proj4_options, int DecimateNone, 
                 frame_t focus, const char *polygons, bool CheckDuplicates, 
                 double* & x, double* & y, double* & z, int & ndata, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int positive;
  int k,m,n,bkp,nn, format;
  size_t count;
  double dmask;
  char *keep;
  range_t<double> lon_range,lat_range;
  double d;
  char output[1024];
  frame_t plgframe;
  string header="";
   
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  load soundings cloud data 

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
      status=xyz_loadraw (filename, header, proj4_options, x, y, z, &dmask, ndata);
      break;
    case SHAPE:
      status=xyz_load_shp (filename, proj4_options, x, y, z, &dmask, ndata);
      break;
    }
    
  if(status!=0) {
     TRAP_ERR_EXIT(-1,"cannot load souding database, exiting\n");
    }
    

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  frame selection of depth soundings cloud
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/  

  if(focus.initialised()==1) {
    keep=new char[ndata];
    for (m=0;m<ndata;m++) {
      keep[m]=1;
      }
    bkp=ndata;
    for(n=0;n<ndata;n++) {
      if(focus.inside(x[n],y[n])!=1) keep[n]=0;
      }
    status=xyz_reduce(x, y, z, keep, ndata);
    printf("#################################################################\n");
    printf("dataset reduction by frame : \n");
    printf("initial count %d, after reduction %d\n", bkp, ndata);
    delete[] keep;
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  polygon selection of depth soundings cloud
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/  

  if(polygons!=0) {
    bkp=ndata;
    status=xyz_PolygonSelection(polygons, x, y, z, ndata, plgframe);
    printf("#################################################################\n");
    printf("dataset reduction by polygons : %s\n",polygons);
    printf("initial count %d, after reduction %d\n", bkp, ndata);
    }


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  depths correction (sign, vertical reference offset)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/  

//   if(scale!=1) {
//     for (m=0;m<ndata;m++) {
//       z[m]*=scale;
//       }
//     }
//   if(d_offset!=0.) {
//     for (m=0;m<ndata;m++) {
//       z[m]+=d_offset;
//       }
//     }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  convert "zero-hydro related" depths to "mean level related" depths
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/  

//   if(PBMA!="") {
//     int signus=1;
// //     if(negative) signus=+1;
// //     else         signus=-1;
//     status= xyz_hydro(PBMA.c_str(), signus, x, y, z, dmask, ndata, 5);
//     sprintf(output, "%s-meanlevel.xyz",rootname);
//     status=xyz_save (output, x, y, z, dmask, ndata);
//     }
  
  if(debug) {
    sprintf(output, "%s-dump.xyz",rootname);
    status=xyz_save (output, x, y, z, dmask, ndata);
    }
    
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_loadmap_UG(char *filename, const char *rootname, char *proj4_options, 
                        double resolution, double factor, double grid_resolution,
                        int DecimateNone, double ElementMaxSize, double ElementMaxRatio, 
                        frame_t focus, const char *e_poly, const char *polygons, frame_t mapping, bool CheckDuplicates, bool meshout, mesh_t & delaunay, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int positive;
  int k,m,n,ndata,bkp,nn, format;
  size_t count;
  double dmask;
  mesh_t mesh;
  char *keep;
  range_t<double> lon_range,lat_range;
  double *x,*y,*z,d;
  grid_t sTOPOgrid;
  grid_t TOPOgrid;
  float *TOPOdata=0, TOPOmask=1.e+10;
  char output[1024];
  frame_t plgframe;
  string header="";
   
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  load soundings cloud data 

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=load_cloud(filename, rootname, proj4_options, DecimateNone, 
                    focus, polygons, CheckDuplicates, x, y, z, ndata, debug);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  establish completion frame
  
  - if mapping frame not specified, addition frame is a limited dilatation of 
    data frame
  
  - if mapping frame specified, addition frame is a limited dilatation of 
    mapping frame
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/  

/*------------------------------------------------------------------------------
  compute data points geographic range */
  lon_range=range_t<double>(x,ndata);
  lat_range=range_t<double>(y,ndata);

  frame_t dataframe(lon_range,lat_range);
  
/*------------------------------------------------------------------------------
  this is an issue when data outside of mapping */
  if(mapping.xmin==INFINITY) mapping=dataframe;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  decimation : gather data by cells (clustering)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
//   keep=new char[ndata];
//   for (m=0;m<ndata;m++) {
//     keep[m]=0;
//     }

//   if(DecimateNone==0) {
//     printf("#################################################################\n");
//     printf("decimate node set\n");
//     double r=resolution;
//     if(r==0.0) r=grid_resolution/2.0;
//     count=clusterize(mapping, x, y, z, keep, ndata, r, factor);
//     }
//   else {
//     for (m=0;m<ndata;m++) {
//       keep[m]=1;
//       }
//     count=ndata;
// /*------------------------------------------------------------------------------
//     check duplicates */
//     if(CheckDuplicates) {
//       status=check_duplicate(x, y, z, ndata, keep, 0.001);
//       status=xyz_reduce(x, y, z, keep, ndata);
//       }
//    }
  
  printf("#################################################################\n");
  printf("generate unstructured mesh\n\n");
  
  mesh.nvtxs=ndata;
  mesh.vertices=new vertex_t[mesh.nvtxs];
  
  n=0;
  for(nn=0;nn<ndata;nn++) {
    if(keep[nn]==1) {
      mesh.vertices[n].lon=x[nn];
      mesh.vertices[n].lat=y[nn];
      mesh.vertices[n].h=-z[nn];
      n++;
      }
    }
  
  delete[] x;
  delete[] y;
  delete[] z;
  
  mesh.nvtxs=n;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  initialize triangle output structure and create delaunay triangulation

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
      
  status=fe_delaunay_create(rootname, mesh, delaunay, meshout, debug);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  mesh post-processing (mostly for node's cloud stand-alone mesh)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_delaunay_improve(rootname, delaunay, ElementMaxSize, ElementMaxRatio, debug);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  create bathymetry structured grid 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=store_topo(rootname, delaunay, mapping, grid_resolution, debug);
  
  return(status);
}

