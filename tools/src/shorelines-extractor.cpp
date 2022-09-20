
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
  
  if (meta->type=="") meta->type=(string) "structured";
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int set_quadrangles(const grid_t & zgrid, bool *cellmask, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  allocates and set all that is necessary and only that for contours
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i,j,k,l,m,n;
  size_t q_count,v_count,e_count;
  int *e_index, *v_index;
  
  int v_tab[4][2]={{0,0},{1,0},{1,1},{0,1}};
  int e_tab[4];
  
  grid_t vgrid=map_vgrid(zgrid);
  
  v_index=new int[vgrid.Hsize()];
  
  mesh.destroy();
  
  size_t size=occurence(true, cellmask, zgrid.Hsize());
  
  mesh.ntriangles  =0;
  mesh.nquadrangles=size;
  
  mesh.quadrangles=new quadrangle_t[mesh.nquadrangles];
  
  for(m=0;m<vgrid.Hsize();m++) v_index[m]=-1;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  register vertices
  
    - numerotation increase with column for a given line
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  q_count=0;
  v_count=0;
  for(j=0;j<zgrid.ny;j++) {
    for(i=0;i<zgrid.nx;i++) {
      m=zgrid.nx*j+i;
      if(not cellmask[m]) continue;
/*------------------------------------------------------------------------------
      register element */
      for(k=0;k<4;k++) {
        int ii=i+v_tab[k][0];
        int jj=j+v_tab[k][1];
        size_t mm=vgrid.nx*jj+ii;
        if(v_index[mm]==-1) {
/*------------------------------------------------------------------------------
          register vertex */
          v_index[mm]=v_count;
          v_count++;
          }
        mesh.quadrangles[q_count].vertex[k]=v_index[mm];
        }
      q_count++;
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  finalize vertices
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  mesh.nvtxs=v_count;
  mesh.vertices=new vertex_t[mesh.nvtxs];
  
  for(size_t j=0;j<vgrid.ny;j++) {
    for(size_t i=0;i<vgrid.nx;i++) {
      size_t m=vgrid.nx*j+i;
      if(v_index[m]==-1) continue;
      size_t  n=v_index[m];
//       mesh.vertices[n].lon=vgrid.x[m];
//       mesh.vertices[n].lat=vgrid.y[m];
      double x,y;
      vgrid.xy(i,j,x,y);
      mesh.vertices[n].lon=x;
      mesh.vertices[n].lat=y;
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  register edges:
  
    - numerotation increase with column for a given line
    - offset is number of horizontal edges, then offset is start of vertcal edges
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   int nrows=vgrid.ny-1;
//   int ncols=vgrid.nx-1;
  size_t offset=vgrid.ny*(vgrid.nx-1);
  
/*------------------------------------------------------------------------------
         horizontal edges      vertical edges    */
  size=vgrid.ny*(vgrid.nx-1)+(vgrid.ny-1)*vgrid.nx;
  e_index=new int[size];
  for(m=0;m<size;m++) e_index[m]=-1;
  
  e_count=0;
  q_count=0;
  for(j=0;j<zgrid.ny;j++) {
    for(i=0;i<zgrid.nx;i++) {
      m=zgrid.nx*j+i;
      if(not cellmask[m]) continue;
      e_tab[0]=m;               // lower edge
      e_tab[1]=offset+m+1;      // right edge
      e_tab[2]=m+vgrid.nx-1;    // upper edge
      e_tab[3]=offset+m;        // left  edge
      for(k=0;k<4;k++) {
        size_t mm=e_tab[k];
        if(e_index[mm]==-1) {
/*------------------------------------------------------------------------------
          register vertex */
          e_index[mm]=e_count;
          e_count++;
          }
        mesh.quadrangles[q_count].edges[k]=e_index[mm];
        }
      mesh.quadrangles[q_count].ancestro=m;
      q_count++;
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  finalize edges
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  mesh.nedges=e_count;
  mesh.edges=new edge_t[mesh.nedges];
  
  for(m=0;m<mesh.nquadrangles;m++) {
    for(k=0;k<4;k++) {
      l=(k+1) % 4;
      n=mesh.quadrangles[m].edges[k];
      int n1=mesh.quadrangles[m].vertex[k];
      int n2=mesh.quadrangles[m].vertex[l];
      mesh.edges[n].extremity[0]=n1;
      mesh.edges[n].extremity[1]=n2;
      mesh.vertices[n1].nngh++;
      mesh.vertices[n2].nngh++;
      mesh.edges[n].shared[mesh.edges[n].nshared]=m;
      mesh.edges[n].nshared++;
      }
    }
  
  for(n=0;n<mesh.nvtxs;n++) {
    mesh.vertices[n].nedges=mesh.vertices[n].nngh;
    mesh.vertices[n].ngh  =new int[mesh.vertices[n].nngh];
    mesh.vertices[n].edges=new int[mesh.vertices[n].nedges];
    mesh.vertices[n].nngh  =0;
    mesh.vertices[n].nedges=0;
    }
  
  for(n=0;n<mesh.nedges;n++) {
    int n1=mesh.edges[n].extremity[0];
    int n2=mesh.edges[n].extremity[1];
    mesh.vertices[n1].ngh[mesh.vertices[n1].nngh]=n2;
    mesh.vertices[n1].nngh++;
    mesh.vertices[n1].edges[mesh.vertices[n1].nedges]=n;
    mesh.vertices[n1].nedges++;
    mesh.vertices[n2].ngh[mesh.vertices[n2].nngh]=n1;
    mesh.vertices[n2].nngh++;
    mesh.vertices[n2].edges[mesh.vertices[n2].nedges]=n;
    mesh.vertices[n2].nedges++;
    }
  
  delete[] v_index;
  delete[] e_index;
  
  vgrid.free();
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float *reduced_01(grid_t & grid, float* z, float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  optimize contours detection by masking areas where no contours could be 
  expected
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  float *masked=new float[grid.Hsize()];
  
  for(size_t m=0;m<grid.Hsize();m++) masked[m]=z[m];

  int depth=3;
  int border=max(depth+2,3);
 
#pragma omp parallel for
  for(int j=border;j<grid.ny-border;j++) {
    for(int i=border;i<grid.nx-border;i++) {
      size_t k=grid.nx*j+i;
      float value=z[k];
      if( (value!=1) and  (value!=-1) ) continue;
      int count=0;
/*------------------------------------------------------------------------------
      for a given value, see if neighbour cells show identical values*/
      for(int jj=j-depth;jj<j+depth+1;jj++) {
        for(int ii=i-depth;ii<i+depth+1;ii++) {
          size_t kk=grid.nx*jj+ii;
          if(kk==k)       continue;
          if(z[kk]==mask) continue;
          if(z[kk]!=value) {
            count++;
            goto finish;
            }
          }
        }
    finish:
      if(count==0) {
        masked[k]=mask;
        }
      else {
        masked[k]=value;
        }
      }
    }
  
  return(masked);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float *reduced_02(grid_t & grid, float* z, float mask, float zmin, float zmax)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  optimize contours detection by masking areas where no contours could be 
  expected
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  float *masked=new float[grid.Hsize()];
  
  for(size_t m=0;m<grid.Hsize();m++) masked[m]=z[m];

  int depth=3;
  int border=max(depth+2,3);
 
#pragma omp parallel for
  for(int j=border;j<grid.ny-border;j++) {
    for(int i=border;i<grid.nx-border;i++) {
      size_t k=grid.nx*j+i;
      float value=z[k];
      bool superior=(value>zmax);
      bool inferior=(value<zmin);
      if( (value>=zmin) and  (value<=zmax) ) continue;
      int count=0;
/*------------------------------------------------------------------------------
      for a given value, see if neighbour cells show identical values*/
      if(inferior)
      for(int jj=j-depth;jj<j+depth+1;jj++) {
        for(int ii=i-depth;ii<i+depth+1;ii++) {
          size_t kk=grid.nx*jj+ii;
          if(kk==k)       continue;
          if(z[kk]==mask) continue;
          if(z[kk]>=zmin) {
            count++;
            goto finish;
            }
          }
        }
      if(superior)
      for(int jj=j-depth;jj<j+depth+1;jj++) {
        for(int ii=i-depth;ii<i+depth+1;ii++) {
          size_t kk=grid.nx*jj+ii;
          if(kk==k)       continue;
          if(z[kk]==mask) continue;
          if(z[kk]<=zmax) {
            count++;
            goto finish;
            }
          }
        }
    finish:
      if(count==0) {
        masked[k]=mask;
        }
      else {
        masked[k]=value;
        }
      }
    }
  
  return(masked);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int isocontours(const string & rootname, const string & ShorelinesFile, grid_t & grid, float *z, float mask, vector<float> & isovalues, bool smooth)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, verbose=0;
  mesh_t mesh;
  bool show_mesh=false, debug=false;
  float *UGarray;
  vector<plg_t> contours;

  string meshfile=rootname+"-quandrangle.nei";
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  mask cells that are unnecessary in isocontour computation
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  float zmin,zmax;
  quick_sort(isovalues);
    
  zmin=isovalues[0]-10;
  zmax=isovalues[isovalues.size()-1]+10;
  
  float *masked=reduced_02(grid, z, mask, zmin, zmax);
//   size_t size=occurence(mask, masked,grid.Hsize());
  
  bool *cellmask=new bool[grid.Hsize()];
  
#pragma omp parallel for
  for(size_t m=0;m<grid.Hsize(); m++) {
    if(masked[m]==mask) cellmask[m]=false;
    else cellmask[m]=true;
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  pre-process quadrangle mesh (faster)
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  status=set_quadrangles(grid, cellmask, mesh);
  
  UGarray=new float[mesh.nquadrangles];
  for(size_t m=0;m<mesh.nquadrangles;m++) {
    size_t mm=mesh.quadrangles[m].ancestro;
    UGarray[m]=masked[mm];
    }
  
  if(debug) status=fe_savemesh(meshfile.c_str(),MESH_FILE_FORMAT_TRIGRID, mesh);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  compute isocontour polygons
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  vector<plg_t> exhaustive;
  
  for(int k=0;k<isovalues.size();k++) {
    string output;
    std::stringstream tmp;
    tmp.width(4);
    tmp.precision(2);
    tmp.setf( std::ios::fixed, std:: ios::floatfield );
    tmp << fabs(isovalues[k]);
    
    plg_destroy_entries(contours);
    contours.clear();
    
    STDOUT_BASE_LINE("treating iso-contour=%lf\n",isovalues[k]);
    
    contours=plg_FindCountours(grid, masked, mask, isovalues[k], mesh, UGarray, show_mesh, debug, verbose);
      
    status=plg_cartesian((projPJ) 0, contours);
    
/*------------------------------------------------------------------------------
    low-level line smoothing, must be done prior to decimation */
    if(smooth) plg_smoothlines(contours);
    
/*------------------------------------------------------------------------------
    iso-contour detection can create a lot of aligned points, remove them */
    status=plg_decimate(contours, 0.0, debug);
    
    for(int s=0;s<contours.size();s++){
      double area=contours[s].area();
      if(debug) printf("polygon %d, isovalue=%f: area=%g\n", s, isovalues[k], area);
      }
    status=plg_add_entries(exhaustive, contours);
    if(isovalues[k]<0)
      output=rootname+"-contours-" + tmp.str() +".plg";
    else
      output=rootname+"-contours+" + tmp.str() +".plg";
    status=plg_save(output.c_str(), PLG_FORMAT_SCAN, contours);
    }
  
  status=plg_save(ShorelinesFile.c_str(), PLG_FORMAT_SCAN, exhaustive);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int surround(grid_t & grid, int i, int j, bool *passed, float *z, int *founds, int & nfounds)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count=0;
  
  size_t m=grid.nx*j+i;
  
  vector<int> tab;
  
  tab.push_back(-1);         // left  cell
  tab.push_back(grid.nx);    // upper cell
  tab.push_back(1);          // right cell
  tab.push_back(-grid.nx);   // lower cell
      
  for(int k=0;k<tab.size();k++) {
    size_t n=m+tab[k];
    if(n<0 or n> grid.Hsize()-1) continue;
    if(z[n]==z[m] and not passed[n]) {
      founds[nfounds]=n;
      nfounds++;
      passed[n]=true;
      count++;
      }
    }
  
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int scan_topo(string options, string rootname, string EchoFile, float zmin, float zmax, bool swap, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, verbose=0;
  grid_t grid;
  float *z, mask;
  string filename, varname;
  statistic_t s;
  
  metafield_t meta;
  
  parse_FieldOptions(options, &meta, debug);
  
  filename=meta.datafile;
  varname=meta.variable;
  
  if(meta.type=="structured") {
    if(varname=="") {
/*------------------------------------------------------------------------------
      use standard set of varnames for topography */
      status=topo_loadfield(filename.c_str(), &grid, &z, &mask, debug);
      }
    else {
      status=topo_loadfield(filename.c_str(), varname.c_str(), &grid, &z, &mask, debug);
      }
    }
  else if(meta.type=="xyz") {
    
    }
  else if(meta.type=="unstructured") {
    
    }
    
  if(status!=0) return(status);
  
  s=get_statistics(z, mask, grid.Hsize(),1);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  assume missing value to be driven by shoreline/model targetting (disputable choice)
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  float boarder;
  
  if(swap) boarder=+10000.0;
  else     boarder=-10000.0;
    
  for(size_t m=0;m<grid.Hsize(); m++) {
    if(z[m]==mask) z[m]=boarder;
    }
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  mask out of range data
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   for(size_t m=0;m<grid.Hsize(); m++) {
//     if(z[m]>zmax) z[m]=mask;
//     }
//  
//   for(size_t m=0;m<grid.Hsize(); m++) {
//     if(z[m]<zmin) z[m]=mask;
//     }
  float *masked=reduced_02(grid, z, mask, zmin, zmax);
    
  s=get_statistics(z, mask, grid.Hsize(),1);
  
  string out;
  status=topo_save(EchoFile, "netcdf", grid, masked, mask, debug);
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int scan_SRTM(string filename, string rootname, string EchoFile, float zmin, float zmax, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, verbose=0;
  grid_t grid;
  
  poc_data_t<float> elevation;
  poc_data_t<double> longitude, latitude;
  float *z, mask;
  bool *seed, *keep, *passed;
  size_t size, previous, k;
  statistic_t s;
    
//   string EchoFile=rootname+"-echo.nc";
  
  status=poc_get_grid(filename, "Band1", &grid, verbose, -1);
  if(status !=0) TRAP_ERR_EXIT(status,"poc_get_grid error\n");
 
  status=elevation.init(filename,"Band1",1);
  status=longitude.init(filename,"longitude",1);
  status=latitude .init(filename,"latitude",1);

  status=elevation.read_data(filename,-1,1);
  status=longitude.read_data(filename,-1,1);
  status=latitude .read_data(filename,-1,1);

  mask=elevation.mask;
  
  z=elevation.data;
  
  s=get_statistics(z, mask, grid.Hsize(),1);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  assume missing value to be ocean (disputable choice)
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(size_t m=0;m<grid.Hsize(); m++) {
    if(z[m]==mask) z[m]=0;
    }
 
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  mask out of range data
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(size_t m=0;m<grid.Hsize(); m++) {
    if(z[m]>zmax) z[m]=mask;
    }
 
  for(size_t m=0;m<grid.Hsize(); m++) {
    if(z[m]<zmin) z[m]=mask;
    }


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  address non-zero water levels
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  seed=new bool[grid.Hsize()];
  keep=new bool[grid.Hsize()];
  passed=new bool[grid.Hsize()];
  if(passed==0) TRAP_ERR_EXIT(-1,"memory allocation failure\n");

#pragma omp parallel for
  for(size_t m=0;m<grid.Hsize(); m++) {
    seed[m]  =false;
    keep[m]  =false;
    passed[m]=false;
    }
  
/*------------------------------------------------------------------------------
  flatness, preliminary selection (masking) */
#pragma omp parallel for
  for(int j=0;j<grid.ny; j++) {
    for(int i=0;i<grid.nx; i++) {
      size_t m=grid.nx*j+i;
      if(z[m]==mask) continue;
      int count=0;
      for(int jj=max(0,j-1);jj<min(grid.ny,j+2);jj++) {
        for(int ii=max(0,i-1);ii<min(grid.nx,i+2);ii++) {
          size_t n=jj*grid.nx+ii;
          if(z[n]==z[m]) {
            count++;
            }
          }
        }
/*------------------------------------------------------------------------------
      less 4 cells with similar values around, mask it */
      if(count<4) z[m]=mask;
      }
    }


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  seeds are starting cells to create clusters of identical value
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
/*------------------------------------------------------------------------------
  seek for seed (cell with all neigbours having identical values) */
  for(int j=0;j<grid.ny; j++) {
    for(int i=0;i<grid.nx; i++) {
      size_t m=grid.nx*j+i;
      if(z[m]==mask) continue;
      int count=0;
      for(int jj=max(0,j-1);jj<min(grid.ny,j+2);jj++) {
        for(int ii=max(0,i-1);ii<min(grid.nx,i+2);ii++) {
          size_t n=jj*grid.nx+ii;
          if(z[n]==z[m]) {
            count++;
            }
          }
        }
      if(count==9) seed[m]=true;
      if(seed[m]) {
        for(int jj=max(0,j-1);jj<min(grid.ny,j+2);jj++) {
          for(int ii=max(0,i-1);ii<min(grid.nx,i+2);ii++) {
            size_t n=jj*grid.nx+ii;
            keep[n]=true;
            }
          }
        }
      }
    }
  
  size=occurence(true, keep, grid.Hsize());
  printf("step #1 : kept=%d\n", size);
  
  size=0;
  previous=0;
  k=0;
  vector<int> tab;
  
  tab.push_back(-1);        // left  cell
  tab.push_back(grid.nx);   // upper cell
  tab.push_back(1);         // right cell
  tab.push_back(-grid.nx);  // lower cell
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  aggregate around seeds 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  do {
    previous=size;
    for(int j=0;j<grid.ny; j++) {
      for(int i=0;i<grid.nx; i++) {
        size_t m=grid.nx*j+i;
        if(not keep[m]) continue;
        if(seed[m])     continue;
        for(int k=0;k<tab.size();k++) {
          size_t n=m+tab[k];
          if(n<0 or n> grid.Hsize()-1) continue;
          if(z[n]==z[m]) {
            keep[n]=true;
            }
          }
        seed[m]=true;
        }
      }
    size=occurence(true, keep, grid.Hsize());
    if(debug) printf("iteration %d : %d value kept\n",k,size);
    k++;
    } while(size!=previous and k < 100);
    
  for(size_t m=0;m<grid.Hsize(); m++) {
    if(not keep[m]) z[m]=mask;
    }
  
  delete[] seed;
     
  size=occurence(true, keep, grid.Hsize());
  printf("step #2 : kept=%d\n", size);
     
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  select "big enough" cluster of identical values; mask others
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
/*------------------------------------------------------------------------------
  clustering */
  int *found=new int[grid.Hsize()];
  if(found==0) TRAP_ERR_EXIT(-1,"memory allocation failure\n");
  
  bool *done=new bool[grid.Hsize()];
  if(done==0) TRAP_ERR_EXIT(-1,"memory allocation failure\n");

  for(size_t m=0;m<grid.Hsize(); m++) done[m]=false;
  
  for(int j=0;j<grid.ny; j++) {
    for(int i=0;i<grid.nx; i++) {
      size_t m=grid.nx*j+i;
      if(not keep[m]) continue;
      if(passed[m])   continue;
      if(done[m])     continue;
      int count=0, nfounds=0;
      found[0]=m;
      nfounds++;
      int start=0;
      do {
        int n=found[start];
        int jj=n/grid.nx;
        int ii=n-jj*grid.nx;
        count=surround(grid, ii, jj, passed, z, found, nfounds);
        start++;
        } while(start!=nfounds);
      if(debug) printf("start %d : %d value passed (status=%d)\n",m,nfounds,status);
      if(nfounds<1000) {
/*------------------------------------------------------------------------------
        recursion limit not reached, cluster found, but too small */
        for(size_t n=0;n<nfounds; n++) keep[found[n]]=false;
        }
      else {
/*------------------------------------------------------------------------------
        recursion limit not reached, cluster found, big enough */
        for(size_t n=0;n<nfounds; n++) done[found[n]]=true;
        }
      }
    }
  
  size=occurence(true, keep, grid.Hsize());
  printf("step #3 : kept=%d\n", size);
     
  for(size_t m=0;m<grid.Hsize(); m++) {
    if(not keep[m]) z[m]=mask;
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  all dry cells are now masked, and ocean/river cells have 0 or other positive 
  values
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  delete[] keep;
  delete[] passed;
  delete[] done;
  
  status=remove(EchoFile.c_str());
  
  status=longitude.write_data(EchoFile);
  if(status!=0) TRAP_ERR_EXIT(-1,"wrinting %s failed\n", EchoFile.c_str());
  status=latitude.write_data (EchoFile);
  if(status!=0) TRAP_ERR_EXIT(-1,"wrinting %s failed\n", EchoFile.c_str());
    
  status=elevation.write_data(EchoFile);
  if(status!=0) TRAP_ERR_EXIT(-1,"wrinting %s failed\n", EchoFile.c_str());
    
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int extend_grid(grid_t & grid, float* & data, float boarder, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, verbose=0;
  resize_t resize;
  float *z;
  
  grid_t tmp[2];
  
  resize.init(grid, -1, 1, -1, 1);
  
  tmp[0]=map_vgrid(grid);
  tmp[1]=map_vgrid(tmp[0]);
  
  z=new float[tmp[1].Hsize()];
  for(size_t m=0;m<tmp[1].Hsize();m++) z[m]=boarder;
  
  for(int j=0;j<grid.ny;j++) {
    int jj=j-resize.start[1];
    for(int i=0;i<grid.nx;i++) {
      int ii=i-resize.start[0];
      size_t m =grid.Hindex(i,j);
      size_t mm=tmp[1].Hindex(ii,jj);
      z[mm]=data[m];
      }
    }
  
  tmp[0].free();
  grid.free();
  
  delete[] data;
  
  data=z;
  grid=tmp[1];
            
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int isoSRTM_shorelines(string filename, const char *varname, string rootname, string ShorelinesFile, 
                         vector<float> isovalues, bool closed, bool smooth, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, verbose=0;
  grid_t grid;
  
  poc_data_t<float> elevation;
  float *z, mask;
  
  status=poc_get_grid(filename, varname, &grid, verbose, -1);
  if(status !=0) TRAP_ERR_EXIT(status,"poc_get_grid error\n");
 
  status=elevation.init(filename, varname, 1);
  status=elevation.read_data(filename, -1, 1);

  mask=elevation.mask;
  
  if(closed) {
    float boarder=0.0;
    status=extend_grid(grid, elevation.data, boarder, debug);
    }
    
//  vector<plg_t>  plg_GridLimits(grid_t grid, float *mask, bool debug, int verbose);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  SRTM setting : initialize landmask (z buffer), water is -1, land is 1. After
                 SRTM pre-processing, land is cells with mask value
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
#pragma omp parallel for
  for(size_t m=0;m<grid.Hsize(); m++) {
    if(z[m]!=mask) z[m]=-1.0;
    }
  
#pragma omp parallel for
  for(size_t m=0;m<grid.Hsize(); m++) {
    if(z[m]==mask) z[m]=+1.0;
    }
 
  z=elevation.data;
  status=isocontours(rootname, ShorelinesFile, grid, z, mask, isovalues, smooth);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int isoTOPO_shorelines(const string & filename, const vector<string> varnames, const string & rootname, const string & ShorelinesFile, 
                         vector<float> & isovalues, bool closed, bool smooth, bool swap, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, verbose=0;
  grid_t grid;
  int k;
  poc_data_t<float> elevation;
  float *z, mask;
  
  for (k=0; k< varnames.size(); k++) {
    status=poc_get_grid(filename, varnames[k].c_str(), &grid, verbose, -1);
    if(status==0) break;
    }
  if(status !=0) TRAP_ERR_EXIT(status,"poc_get_grid error\n");
 
  status=elevation.init(filename, varnames[k].c_str(), 1);
  status=elevation.read_data(filename, -1, 1);

  mask=elevation.mask;
  
  if(closed) {
    float boarder;
    if(swap) boarder=+10000.0;
    else     boarder=-10000.0;
    status=extend_grid(grid, elevation.data, boarder, debug);
    }
    
//  vector<plg_t>  plg_GridLimits(grid_t grid, float *mask, bool debug, int verbose);

  z=elevation.data;
  status=isocontours(rootname, ShorelinesFile, grid, z, mask, isovalues, smooth);
  
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
    "  %s file1 [ file2 ... ] [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Convert between different coastline formats.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  Show this help and exit.\n"
    "  --debug  \n"
    "  --frame followed by \"[lon_min;lon_max;lat_min:lat_max]\". example: --frame \"[-10.0:0.0;42.5:50]\" !\n"
    "  --swap  swap from shoreline extraction type to model's limits extraction type\n"
    "  --topography followed by \"file<filename> variable=<varname>\" \n"
    "  -r  followed by the root name\n"
    "  -o  followed by the path of the output file. Default from the root name.\n"
    "  -O  followed by the output format\n"
    "  -p  followed by a path to a limiting polygone\n"
    "  -s  followed by the path to the shorelines. Default is /home/data/shorelines/gshhs-2.0/gshhs_f.cst\n"
    "  -v  followed by isocontour value\n"
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
  int n,nitems,status;
  string keyword,PolygonFile="", FrameString="", rootname="", zone="", output;
  string outFormatName="";
  string ShorelinesFile="", ElevationFile="", FilteredElevationFile="";
  string TopographyFile="";
  int inFormat,outFormat;
  vector<plg_t> extraction;
  frame_t frame;
  vector<plg_t> limits;
  point2D_t point;
  bool debug=false;
  int target=0;
  bool decimate=false, smooth=false;
  vector<float> isovalues;
  float tmp;
  double alignment=0.0;
  
  fct_echo( argc, argv);

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
          ShorelinesFile=argv[n+1];
          n++;
          n++;
          break;
        
        case 'r' :
          rootname=argv[n+1];
          n++;
          n++;
          break;
        
        case 'v' :
          nitems=sscanf(argv[n+1],"%f", &tmp);
//           if(nitems==0) 
          isovalues.push_back(tmp);
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
          else if(strncmp("--frame",keyword)==0){
            FrameString= argv[n+1];
            n++;
            n++;
            }
          else if(strncmp("--elevation",keyword)==0){
            ElevationFile= argv[n+1];
            n++;
            n++;
            }
          else if(strncmp("--filtered",keyword)==0){
            FilteredElevationFile= argv[n+1];
            n++;
            n++;
            }
          else if(strncmp("--topography",keyword)==0){
            TopographyFile= argv[n+1];
            n++;
            n++;
            }
          else if(strncmp("--swap",keyword)==0){
/*------------------------------------------------------------------------------
            reverse meaning of shorelines interior/exterior */
            target=1;
            n++;
            }
          else if(strncmp("--decimate",keyword)==0){
/*------------------------------------------------------------------------------
            decimate aligned points in shorelines lines */
            TRAP_ERR_EXIT(-1,"--decimate option deprecated, use instead --decimate-aligned followed by angle value\n");
            n++;
            n++;
            }
          else if(strncmp("--decimate-aligned",keyword)==0){
/*------------------------------------------------------------------------------
            decimate aligned points in shorelines lines */
            decimate=true;
            sscanf(argv[n+1],"%lf",&alignment);
            n++;
            }
          else if(strncmp("--smooth",keyword)==0){
/*------------------------------------------------------------------------------
            smooth staircase lines */
            smooth=true;
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

  if(rootname=="") rootname="anonymous";
  
  if(ShorelinesFile!="" and (ElevationFile!="" or FilteredElevationFile!="")) {
    fprintf(stderr,"shorelines file over-documented, use -s FILE or --elevation FILE but not both together\n");
    print_help(argv[0]);
    return -1;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  set extraction region from zone naming
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(zone!="") {
    status=plg_SetZone(zone.c_str(), frame, rootname, point);
    limits.push_back(plg_t(frame, PLG_INIT_SEPARATE));
    status=plg_cartesian((projPJ) 0, limits);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  set extraction region from frame string : --frame "[-10.0:0.0;42.5:50]"
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(FrameString!="") {
    status=plg_DecodeFrame(FrameString, frame);
    limits.push_back(plg_t(frame, PLG_INIT_SEPARATE));
    status=plg_cartesian((projPJ) 0, limits);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  set extraction region from polygons
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(PolygonFile!="") {
    inFormat=plg_find_format(PolygonFile);
    status=plg_load(PolygonFile, inFormat, limits);
    if(status!=NC_NOERR) NC_TRAP_ERROR(wexit,status,1,"plg_load(\""+PolygonFile+"\",%d,,) error",inFormat);
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  extract shoreline as SRTM iso-countour lines
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(ElevationFile!="" or FilteredElevationFile!="") {
    bool process=false, closed=true, smooth=false;
    float zmin=0.0,zmax=10.0;
    
    isovalues.clear();
    isovalues.push_back(0.0);
        
    if(FilteredElevationFile=="") {
      FilteredElevationFile=rootname+"-filtered.nc";
      process=true;
      }
    ShorelinesFile=rootname+"-shorelines.plg";

    if(process) status=scan_SRTM(ElevationFile, rootname, FilteredElevationFile, zmin, zmax, debug);
    char *varname="Band1";
    status=isoSRTM_shorelines(FilteredElevationFile, varname, rootname, ShorelinesFile, isovalues, closed, smooth, debug);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  extract shoreline as bathymetry iso-countour lines
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(TopographyFile!="") {
    bool process=false, closed=true;
    float zmin,zmax;
    string FilteredFile;
    bool swap=false;
    
    if(target==1) swap=true;
    
    if(isovalues.size()==0) isovalues.push_back(0.0);
    
    quick_sort(isovalues);
    
    zmin=isovalues[0]-10;
    zmax=isovalues[isovalues.size()-1]+10;
    
    FilteredFile=rootname+"-filtered.nc";
    ShorelinesFile=rootname+"-shorelines.plg";

    status=scan_topo(TopographyFile, rootname, FilteredFile, zmin, zmax, swap, debug);
    
//     char *varname="bathymetry";
    vector<string> varnames;
    varnames.push_back("bathymetry");
    varnames.push_back("Band1");
    status=isoTOPO_shorelines(FilteredFile, varnames, rootname, ShorelinesFile, isovalues, closed, smooth, swap, debug);
    }
  
  if(ShorelinesFile=="") {
    fprintf(stderr,"shorelines file not documented, use -s FILE or --elevation \"elevation=XXX filtered=XXX zmin=XXX zmax=XXX\" \n");
    print_help(argv[0]);
    return -1;
    }
  
  if(output=="") output=rootname+"-extraction.plg";
  

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  reduce data set to given region/sampling
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(limits.size()!=0) {
/*------------------------------------------------------------------------------
    extract shorelines dataset*/
    printf("#################################################################\n");
    printf("extract shorelines\n");

    extraction=plg_extract(ShorelinesFile, limits, rootname, PLG_SPHERICAL, target, debug);
    if(extraction.size()==0) goto error;

    if(decimate) for(int k=0;k<3;k++) status=plg_decimate(extraction, alignment, debug);
  
    if(outFormatName!="") {
      outFormat=plg_format_from_name(outFormatName);
      if(outFormat==PLG_FORMAT_UNKNOWN){
        plg_print_formats();
        TRAP_ERR_EXIT(-1,"exiting\n");
        }
      }
    else {
      outFormat=PLG_FORMAT_SCAN;
      }
    
    status=plg_save(output.c_str(), outFormat, extraction);
    }
  
  STDOUT_BASE_LINE(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(0);
error:
  STDOUT_BASE_LINE("error detected, quit ... \n");
  exit(-1);
}
