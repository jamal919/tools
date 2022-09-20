
/*******************************************************************************

  T-UGO tools, 2006-2013

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

\brief iso-lines functions
*/
/*----------------------------------------------------------------------------*/


#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <string.h>
#include <fstream>
#include <unistd.h> //unlink

#include "fe.h"
#include "fe-proto.h"
#include "functions.h"
#include "filter.h"

#include "polygones.h"

class seek_t {
private :
public :
  int first;
  bool do_interrupted;
  
  seek_t() {
    first=0;
    do_interrupted=true;
  }
};


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_FindCountoursStart_fast(const mesh_t & mesh, char *flag, seek_t & seek, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m,m1,n;
  
/*------------------------------------------------------------------------------
  get the first you can */
  for(m=seek.first;m<mesh.nedges && !flag[m];m++){}
  seek.first=m;
  
/*------------------------------------------------------------------------------
  go to the end */
  n=-1;
  do{
    m1=m;
    flag[m]=2;
    
    if(mesh.edges[m].extremity[0]==n)
      n=mesh.edges[m].extremity[1];
    else
      n=mesh.edges[m].extremity[0];
    
    vertex_t *vertex=&mesh.vertices[n];
    int k;
    for(k=0;k<vertex->nedges;k++){
      m=vertex->edges[k];
      if(m==m1)continue;
      if(flag[m]==1) break;
      }
    
    if(k==vertex->nedges)//line ends or wrap
      m=m1;
    
    }while(flag[m]!=2);
  
  return m;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_FindCountoursStart_old(const mesh_t & mesh, char *flag, seek_t & seek, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m,n;
  
/*------------------------------------------------------------------------------
  first try unclosed contours (i.e. intercepting grid limits, masked regions) */
  if(seek.do_interrupted) {
    for(n=seek.first;n<mesh.nvtxs;n++) {
      int count=0;
      for(int k=0;k<mesh.vertices[n].nedges;k++) {
        m=mesh.vertices[n].edges[k];
        if(flag[m]==1) count++;
        }
      if(count!=1) continue;
      for(int k=0;k<mesh.vertices[n].nedges;k++) {
        m=mesh.vertices[n].edges[k];
        if(flag[m]==1) {
          if(verbose)
            printf("field borderline starting point %d\n", m);
          seek.first=n;
          return(m);
          }
        }
      }
    }
  
  if(seek.do_interrupted) {
    seek.do_interrupted=false;
    seek.first=0;
    }
  
/*------------------------------------------------------------------------------
  then try closed contours */
  for(n=seek.first;n<mesh.nvtxs;n++) {
    int count=0;
    for(int k=0;k<mesh.vertices[n].nedges;k++) {
      m=mesh.vertices[n].edges[k];
      if(flag[m]==1) count++;
      }
    if(count!=2) continue;
    for(int k=0;k<mesh.vertices[n].nedges;k++) {
      m=mesh.vertices[n].edges[k];
      if(flag[m]==1) {
        if(verbose)
          printf("field inside line starting point %d\n", m);
        seek.first=n;
        return(m);
        }
      }
    }
  
  //if(verbose)
    printf("what ????\n");
  return(-1);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_FindCountoursNext_new(mesh_t & mesh, char *flag, int previous_vertex, int previous_edge)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  return next edge of vertex "previous" along contour */
{
  int status=0;
  int m,m1=-1;
  int count=0;
  int rotation=1;
  
  for(int k=0;k<mesh.vertices[previous_vertex].nedges;k++) {
    m=mesh.vertices[previous_vertex].edges[k];
    if(flag[m]){
      count++;
      if(m1<0)
        m1=m;
      }
    }
  if(count==2) {
    return(-2);
    }
  if(count==3) {
    status=fe_SortEdges(mesh, previous_vertex, previous_edge, rotation);
    for(int k=0;k<mesh.vertices[previous_vertex].nedges;k++) {
      m=mesh.vertices[previous_vertex].edges[k];
      if(flag[m]) return(m);
      }
    return(-3);
    }
  if(count!=1) return(-1);
  
  return(m1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_FindCountoursNext_old(mesh_t & mesh, char *flag, int previous_vertex, int previous_edge)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  return next edge of vertex "previous" along contour */
{
  int status=0;
  int m;
  int count=0;
  int rotation=1;
  
  vertex_t & vertex=mesh.vertices[previous_vertex];
  
  if(vertex.edges==0) TRAP_ERR_EXIT(-1,"vertices cross tables not initialised, abort code\n");
  
/*------------------------------------------------------------------------------
  count number of elligible edge */
  for(int k=0;k<vertex.nedges;k++) {
    m=vertex.edges[k];
    if(flag[m]) count++;
    }
  
  if(count==2) {
    return(-1);
    }
  
  if(count==3) {
    status=fe_SortEdges(mesh, previous_vertex, previous_edge, rotation);
    for(int k=0;k<vertex.nedges;k++) {
      m=vertex.edges[k];
      if(flag[m]) return(m);
      }
    return(-1);
    }
  
  if(count!=1) return(-1);
  
/*------------------------------------------------------------------------------
  one elligible edge only, return edge index */
  for(int k=0;k<vertex.nedges;k++) {
    m=vertex.edges[k];
    if(flag[m]) return(m);
    }
  
  return(status);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void setVtxAndEdge(mesh_t *mesh,const grid_t & vgrid,int j,int i,int m,int *k,int m2Offset,int l,int sharedOffset)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vertex_t *vertex=&mesh->vertices[m];
  int m2=m+m2Offset;
  
  vertex->ngh[*k]=m2;
  vertex->edges[*k]=l;
  
  if(m2Offset>0){
    edge_t *edge=&mesh->edges[l];
    int *nshared=&edge->nshared;
    
    edge->extremity[0]=m;
    edge->extremity[1]=m2;
    
    if(sharedOffset>0) TRAP_ERR_EXIT(ENOEXEC,"programming error\n");
    
    *nshared=0;
    m-=j;
    if(i<vgrid.nx-1 && j<vgrid.ny-1){
      edge->shared[*nshared]=m;
      (*nshared)++;
      }
    if(sharedOffset){
      edge->shared[*nshared]=m+sharedOffset;
      (*nshared)++;
      }
    
    }
  
  (*k)++;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void vgrid2edges_of_quadrangles(const grid_t & zgrid, const grid_t & vgrid, mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  allocates and set all that is necessary and only that for contours
------------------------------------------------------------------------------*/
{
  int sharedOffset;
  int i,j,k,l,m=0;
  
  mesh->destroy();
  /* mesh->vertices[].{edges[],lon,lat,ngh[]} */
  /* mesh->edges[].{extremity[],shared[],nshared} */
  mesh->nvtxs=vgrid.Hsize();
  mesh->nedges= (vgrid.nx-1)*vgrid.ny+ vgrid.nx*(vgrid.ny-1);
  mesh->ntriangles=0;
  mesh->nquadrangles=0;
  mesh->allocate();
  mesh->nquadrangles=zgrid.Hsize();
  
  for(j=0;j<vgrid.ny;j++)
    for(i=0;i<vgrid.nx;i++){
    vertex_t *vertex=&mesh->vertices[m];
    
    vgrid.xy(m,vertex->lon,vertex->lat);
    
    vertex->nngh=(i>0) + (i<vgrid.nx-1) + (j>0) + (j<vgrid.ny-1);
    vertex->nedges=vertex->nngh;
    vertex->allocate();
    
    l=j*(vgrid.nx-1+vgrid.nx);//index of the first horizontal edge
    
    k=0;
    
    if(i<vgrid.nx-1){
      if(j>0)
        sharedOffset=-zgrid.nx;
      else
        sharedOffset=0;
      
      setVtxAndEdge(mesh,vgrid,j,i,m,&k,
        1,
        l+i,
        sharedOffset);
      }
    
    if(j>0){
      setVtxAndEdge(mesh,vgrid,j,i,m,&k,
        -vgrid.nx,
        l-vgrid.nx+i,
        1);
      }
    
    if(i>0){
      setVtxAndEdge(mesh,vgrid,j,i,m,&k,
        -1,
        l+i-1,
        1);
      }
    
    if(j<vgrid.ny-1){
      setVtxAndEdge(mesh,vgrid,j,i,m,&k,
        vgrid.nx,
        l+vgrid.nx-1+i,
        -(i>0));
      }
    
    m++;
    }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t> plg_FindCountours(const grid_t & zgrid, float *zarray, float mask, const vector<float> & values, mesh_t & mesh, float* & UGarray, bool show_mesh, bool debug, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*--------------------------- IMPORTANT NOTE -----------------------------------
  PLEASE DO NOT WRITE (e.g. with printf,...) TO STANDARD OUTPUT IN THIS FUNCTION:
                    THIS FUNCTION IS USED IN POCViP !
------------------------------------------------------------------------------*/
  int status;
  vector<plg_t> contours;
  bool masked=false;
  struct timeval b4;
//   bool debug=(verbose>0);
  char *echo=0;
  
  if(mesh.nquadrangles==0) {
/*------------------------------------------------------------------------------
    at first use, build the quadrangle mesh */
    if(verbose>0) STDERR_BASE_LINE("%s: create quadrangle mesh from %dx%d zgrid\n",__func__,zgrid.nx,zgrid.ny);
    grid_t vgrid=map_vgrid(zgrid, zarray, mask,0);
    if(zgrid.projection){
      float *varray;
      char  *landmask=0;
      varray=map_extrapolate_z2v(zgrid, vgrid, zarray, mask);
      grid_t sgrid=map_get_spherical(vgrid.projection,vgrid);
      landmask=new char[sgrid.Hsize()];
      for(int m=0;m<sgrid.nx*sgrid.ny;m++) {
        if(varray[m]==mask) {
          landmask[m]=0;
          }
        else{
          landmask[m]=1;
          }
        }
      if(debug || show_mesh) {
        echo=poc_strdup("fe_FGrid2Quadrangle.nc");
        }
      status=fe_FGrid2Quadrangle(sgrid, vgrid, landmask, varray, mask, mesh, echo, debug);  // need  to be optimized
      deletep(&echo);
      delete[] varray;
      delete[] landmask;
      sgrid.free();
      if(status!=0) return(contours);
      masked=(mesh.nquadrangles!=zgrid.Hsize());
      if(masked) {
        sgrid=map_get_spherical(zgrid.projection,zgrid);
        UGarray=new float[mesh.nquadrangles];
/*------------------------------------------------------------------------------
           spherical coordinates interpolation will fail in polar regions */
//         for(int m=0;m<mesh.nquadrangles;m++) {
//           double t,p;
//           status=fe_position(mesh, mesh.quadrangles[m], &t, &p,0);
//           status=map_interpolation(sgrid,zarray,mask,t,p,&UGarray[m]);
//           if(status!=0) {
//             status=map_interpolation(sgrid,zarray,mask,t,p,&UGarray[m]);
//             }
//          }
/*------------------------------------------------------------------------------
        thus use cartesian coordinates interpolation */
        double *x=new double[mesh.nquadrangles];
        double *y=new double[mesh.nquadrangles];
        for(int m=0;m<mesh.nquadrangles;m++) {
          status=fe_position(mesh, mesh.quadrangles[m], &x[m], &y[m],0);
           }
        status=geo_to_projection (zgrid.proj4options, x, y, mesh.nquadrangles);
        for(int m=0;m<mesh.nquadrangles;m++) {
          status=map_interpolation(zgrid,zarray,mask, x[m], y[m],&UGarray[m]);
          if(status!=0) {
            status=map_interpolation(zgrid,zarray,mask,x[m], y[m],&UGarray[m]);
            }
          }
        delete[] x;
        delete[] y;
        sgrid.free();
        }
      else {
        UGarray=zarray;
        }
      }
    else{
      vgrid2edges_of_quadrangles(zgrid,vgrid,&mesh);
      }
    vgrid.projection=0;/* given by map_vgrid(zgrid,...) */
    vgrid.free();
    }
  
//   if(!zgrid.projection){
  if(UGarray==0){
    UGarray=zarray;
    }
  
  if(values.size()<=0) TRAP_ERR_RETURN(contours,verbose,"%d isovalues!\n",values.size());
  
  gettimeofday(&b4);
  
/*------------------------------------------------------------------------------
  flag eligible edges */
  int count=0,n;
  char *flag=aset(mesh.nedges,'\0');
  
#pragma omp parallel for
  for(n=0;n<mesh.nedges;n++) {
    edge_t *edge=&mesh.edges[n];
    int m1,m2;
    float am1,am2;
    bool c;
    int i;
    float value;
    
    if(edge->nshared<=1) continue;
    
    m1=edge->shared[0];
    m2=edge->shared[1];
    
    am1=UGarray[m1];
    if(am1==mask) {
      if(verbose>1) STDERR_BASE_LINE("%d %d\n",n,m1);
      continue;
      }
    
    am2=UGarray[m2];
    if(am2==mask) {
      if(verbose>1) STDERR_BASE_LINE("%d %d\n",n,m2);
      continue;
      }
    
    for(i=0;i<values.size();i++){
      value=values[i];
      c=( (am1>value) == (am2>value) );
      if(!c) break;
      }
    
    if(c) continue;
    
    flag[n]=1;
#pragma omp critical
      {
      count++;
      }
    }
  
  if(verbose>0) STDERR_BASE_LINE("%gs:%d edges flagged. creating polygones...\n",difftime(&b4),count);
  
/*------------------------------------------------------------------------------
  create polygones */
  int previous_edge;
  seek_t seek;
  double startT=0.,nextT=0.;
  
  while(count>0) {
    vector<int> nodes;
    struct timeval b4;
    gettimeofday(&b4);
    n=plg_FindCountoursStart_fast(mesh, flag, seek, verbose);
    //n=plg_FindCountoursStart_old(mesh, flag, seek, verbose);
    if(n<0) break;
    flag[n]=0;
    count--;
    nodes.push_back(mesh.edges[n].extremity[0]);
    nodes.push_back(mesh.edges[n].extremity[1]);
    previous_edge=n;
    startT+=difftime(&b4);
    do {
      int vertex=nodes.back();
      //n=plg_FindCountoursNext_new(mesh,flag,vertex,previous_edge);
      n=plg_FindCountoursNext_old(mesh,flag,vertex,previous_edge);
//       if(n<=0) break;
      if(n<0) break;
      flag[n]=0;
      count--;
      if(mesh.edges[n].extremity[0]==vertex)
        nodes.push_back(mesh.edges[n].extremity[1]);
      else
        nodes.push_back(mesh.edges[n].extremity[0]);
      previous_edge=n;
//       } while(n>0);
      } while(n>=0);
    nextT+=difftime(b4);
    plg_t polygon;
    polygon.init(nodes.size(),PLG_INIT_SEPARATE);
    for(int k=0;k<polygon.npt;k++) {
      n=nodes[k];
      polygon.t[k]=mesh.vertices[n].lon;
      polygon.p[k]=mesh.vertices[n].lat;
      }
    contours.push_back(polygon);
    }
  
  if(verbose>0) STDERR_BASE_LINE("%gs=%g+%g+:polygones done. cleaning up.\n",difftime(b4),startT,nextT);
  
  deletep(&flag);
  
//  status=plg_save("contours.plg", PLG_FORMAT_SCAN, contours);
  
  return(contours);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t> plg_FindCountours(const grid_t & zgrid, float*array, float mask, float value, mesh_t & mesh, float* & UGarray, bool show_mesh, bool debug, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<float> vs;
  
  vs.push_back(value);
  
  vector<plg_t> contours;
  
  contours=plg_FindCountours(zgrid, array, mask, vs, mesh, UGarray, show_mesh, debug, verbose);
  
  return contours;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_t plg_GridLimits(const grid_t & grid, int increment)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,n=0;
  double x,y;
  
  plg_t polygon;
  
  polygon.init(
    ((grid.nx-1+increment-1)/increment+1)*2+
    ((grid.ny-1+increment-1)/increment+1)*2-3,PLG_INIT_SEPARATE);
  
  j=0;
  for (i=0;i<grid.nx;i+=increment) {
    grid.xy(i,j,x,y);
    polygon.t[n]=x;
    polygon.p[n]=y;
    n++;
    }
  if(i-increment!=grid.nx-1) {
    i=grid.nx-1;
    grid.xy(i,j,x,y);
    polygon.t[n]=x;
    polygon.p[n]=y;
    n++;
    }
  
  i=grid.nx-1;
  for (j=increment;j<grid.ny;j+=increment) {
    grid.xy(i,j,x,y);
    polygon.t[n]=x;
    polygon.p[n]=y;
    n++;
    }
  if(j-increment!=grid.ny-1) {
    j=grid.ny-1;
    grid.xy(i,j,x,y);
    polygon.t[n]=x;
    polygon.p[n]=y;
    n++;
    }
  
  j=grid.ny-1;
  for (i=grid.nx-1-increment;i>=0;i-=increment) {
    grid.xy(i,j,x,y);
    polygon.t[n]=x;
    polygon.p[n]=y;
    n++;
    }
  if(i+increment!=0) {
    i=0;
    grid.xy(i,j,x,y);
    polygon.t[n]=x;
    polygon.p[n]=y;
    n++;
    }
  
  i=0;
  for (j=grid.ny-1-increment;j>=0;j-=increment) {
    grid.xy(i,j,x,y);
    polygon.t[n]=x;
    polygon.p[n]=y;
    n++;
    }
  if(j+increment!=0) {
    j=0;
    grid.xy(i,j,x,y);
    polygon.t[n]=x;
    polygon.p[n]=y;
    n++;
    }
  
  if(n!=polygon.npt) TRAP_ERR_EXIT(ENOEXEC,"coding error : polygon has %d points set when %d=ceil(%d/%d)*2+ceil(%d/%d)*2-3 points where expected\n",n,polygon.npt,grid.nx,increment,grid.ny,increment);
  
  for(j=n-1;j>=0;j--){
    if( isfinite(polygon.t[j]) and isfinite(polygon.t[j]) )
      continue;
    plg_delete_point(&polygon,j);
    }
  
  return(polygon);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t>  plg_GridLimits(grid_t grid, float *mask, bool debug, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<plg_t> polygons;
  
  float value=0.0, fakemask=1.e+10;
  mesh_t mesh;
  float *UGarray=0;
  
  float *buffer=aset<float>(grid.Hsize(), 1.0f);
  
  for(size_t m=0; m<grid.Hsize(); m++) {
    if(mask[m]==0.0) buffer[m]=-1;
    }
  polygons=plg_FindCountours(grid, buffer, fakemask, value, mesh, UGarray, false, debug, verbose);
  
  delete[] buffer;
  if(UGarray!=buffer) delete[] UGarray;
  
  mesh.destroy();
  
  return(polygons);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_RemoveImbrications(vector<plg_t> & p, double x, double y, int check)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int s,status;
  int count=0;
  int *imbrication=new int[p.size()];
  vector<int> continent;
  vector<point2D_t> points;
  
  for(s=0;s<p.size();s++) imbrication[s]=0;
  
/*------------------------------------------------------------------------------
  create an inside point for each polygon */
  for(int l=0;l<p.size();l++) {
    point2D_t point=p[l].InsidePoint();
    points.push_back(point);
    }
  
/*------------------------------------------------------------------------------
  detect imbricated polygons */

  #pragma omp parallel for
  for(int k=0;k<p.size();k++) {
    for(int l=0;l<p.size();l++) {
      if(l==k) continue;
      point2D_t & point=points[l];
      int flag=0;
      flag=plg_single(p[k],point.x,point.y,&flag);
      if(flag==PLG_POINT_INTERIOR) {
/*------------------------------------------------------------------------------
        inside point of polygon l is inside polygon k */
        #pragma omp critical(imbrication)
        {
        imbrication[l]++;
        }
        }
      }
    }

/*------------------------------------------------------------------------------
  remove imbricated polygons */
  for(int k=0;k<p.size();k++) {
    if(imbrication[k]!=0) {
      count++;
      p[k].npt=0;
      }
    }

  status=plg_compact_entries (p);
  
  delete[] imbrication;
  
  points.clear();
  
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int field_interpolation(const SGfield_t<float> & field, double x, double y, float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result=0;
  int64_t accel;
  
//  result=map_interpolation(*field.grid,field.x,field.mask,x,y,z);
  index_interpolation(*field.grid, x, y, &accel, field.x, field.mask, z, 0);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  SGfield_t<float> Smoothedlandmask(const string & filename, const vector<plg_t> & polygons,const plg_filter_t & plg_filter, int mode, float & value, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  lengthscale : filtering length scale
  
  dx: underlying working (cartesian) grid resolution
  
  mode:
  
    0 -
    1 -
    2 -
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int count, bkp, status;
  float mask=-2;
  float mean_bkp,mean;
  plg_array_t shorelines;
  vector<plg_t> limits;
  SGfield_t<float> field;
  float *r=0;
  
  double lengthscale=plg_filter.lengthscale;
  double dx=plg_filter.dx;
  double factor=plg_filter.factor;
  
  bool tuned=plg_filter.tuned;
  bool enforcement=plg_filter.enforcement;
    
  SGfield_t<float> resolution=plg_filter.resolution;

  shorelines=plg_vector2array(polygons);
  field.grid=new grid_t;
  grid_t & grid=*field.grid;
  
  debug=true;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  define cartesian working grid and compute landmask
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  build working grid (cartesian, structured), apply projection to shorelines */
  status=fe_defgrid(&grid, shorelines.p, shorelines.n, limits, dx);
  
  plg_destroy_entries(limits);
  
/*------------------------------------------------------------------------------
  set interior/exterior flag                                                  */
  printf("#------------------------------------------------------\n");
  printf("step 1: set landmask flag\n");
  float *landmask=set_landmask01(grid, shorelines.p, shorelines.n);
  
/*------------------------------------------------------------------------------
  in case where polygons are limits, not shorelines                           */
  if(mode!=0) {
    printf("boundary mode: use inverted landmask\n");
    for (size_t m=0; m<grid.Hsize();m++) landmask[m]=-landmask[m];
    }
  else {
    printf("shoreline mode: use native landmask\n");
    }
  
/*------------------------------------------------------------------------------
  delete pre-existing file if already existing                                */
  if(filename!="") remove(filename.c_str());

/*------------------------------------------------------------------------------
  save first landmask                                                         */
  if(filename!="") {
    status=save_SGXY(filename.c_str(), grid, landmask, mask, "landmask-0", "NONE", "landmask-0", 1, 1);
    }
  
/*------------------------------------------------------------------------------
  made to preserve mono-block penninsula, to be revised                       */
  if(enforcement) status=enforce_landmask(grid, shorelines.p, shorelines.n, landmask, 1.0, tuned);
  
/*------------------------------------------------------------------------------
  save final landmask                                                         */
  if(filename!="") {
    status=save_SGXY(filename.c_str(), grid, landmask, mask, "landmask", "NONE", "landmask", 0, 0);
    }
  
  float *flag=0;
//  flag=set_landmask03(grid, polygons);
  float *smoothed = new float[grid.Hsize()];
  float *delta    = new float[grid.Hsize()];

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  smooth land/sea flag
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#------------------------------------------------------\n");
  if(resolution.x==0) {
/*------------------------------------------------------------------------------
    fix smoothing lengthscale */
    printf("step 2: smooth landmask with fixed lengthscale=%lf mmap_smooth\n",lengthscale);
    status=map_smooth(LOESS, grid, landmask, mask, lengthscale, smoothed);
    }
  else {
/*------------------------------------------------------------------------------
    variable smoothing lengthscale */
    printf("step 2: smooth landmask with variable lengthscale (%lf x local resolution)\n", factor);
    r=new float[grid.Hsize()];
    set_grid_list(resolution.grid, 0);
    int nprocs=initialize_OPENMP(-1, 0);
    
    projPJ *refs = init_projection_parallel(grid.proj4options,nprocs,true,0);
    
#pragma omp parallel for if(nprocs>1)
    for(int m=0;m<grid.Hsize();m++) {
      double x, y, t, p;
      grid.xy(m,x,y);
//       projection_to_geo(grid.projection, &p, &t, x, y);
      projection_to_geo(refs[omp_get_thread_num()], &p, &t, x, y);
      status=field_interpolation(resolution, t, p, &r[m]);
      if(r[m]==resolution.mask) {
        r[m]=lengthscale;
        }
      else {
        r[m]=min(10.*lengthscale, r[m]*factor);
        }
      }
    
    deletep2D(&refs,nprocs,pj_free);
    
    status=map_smooth(grid, landmask, mask, r, smoothed);
    deletep(&resolution.grid->list);
    }
  
  if(filename!="") {
    status=save_SGXY(filename.c_str(), grid, smoothed, mask, "smoothed-0", "NONE", "smoothed-0", 0, 0);
    }
  
/*------------------------------------------------------------------------------
  made to preserve mono-block penninsula, to be revised                       */
  if(enforcement) status=enforce_landmask(grid, shorelines.p, shorelines.n, smoothed, 1, tuned);
  
  for(int m=0;m<grid.Hsize();m++) {
    if(landmask[m]!=mask) delta[m]=smoothed[m]-landmask[m];
    else delta[m]=mask;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  try to find "best" iso-value (identical sea/land ratio)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  bkp=0;
  mean_bkp=0;
  size_t unmasked=0;
  
  statistic_t sraw=get_statistics(landmask, mask, grid.Hsize(), 0);
  statistic_t ssmt=get_statistics(smoothed, mask, grid.Hsize(), 0);

/*------------------------------------------------------------------------------
  count initial number of land points                                         */
  for(int n=0;n<grid.Hsize();n++) {
    if(smoothed[n]!=mask) {
      if(landmask[n]>0.) {
        bkp++;
        }
      mean_bkp+=landmask[n];
      unmasked++;
      }
    }
  
redo:
  count=0;
  mean=0;
  unmasked=0;
/*------------------------------------------------------------------------------
  count actual number of land points if using isocontour value                */
  for(int n=0;n<grid.Hsize();n++) {
    if(smoothed[n]!=mask) {
      if(smoothed[n]>value) {
        count++;
        }
      mean+=landmask[n];
      unmasked++;
      }
    }
  if(count<bkp) {
    value-=0.01;
    goto redo;
    }
  
  mean/=count;
  mean_bkp/= bkp;
  
  printf("\n %s : best isovalue for smoothing=%f \n\n",__func__, value);
   
  if(filename!="") {
    status=save_SGXY(filename.c_str(), grid, smoothed, mask, "smoothed", "NONE", "smoothed", 0, 0);
    if(flag!=0) status=save_SGXY(filename.c_str(), grid, flag,  mask, "flag",  "NONE", "flag", 0,  0);
    if(r!=0)    status=save_SGXY(filename.c_str(), grid, r,     mask, "scale", "NONE", "scale", 0, 0);
    status=save_SGXY(filename.c_str(), grid, delta, mask, "delta", "NONE", "delta", 0);
    }
  
  delete[] landmask;
  if(flag!=0) delete[] flag;
  if(r!=0)    delete[] r;
  delete[] delta;

  field.x=smoothed;
  field.mask=mask;
  
  delete[] shorelines.p;
  
  return(field);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//   vector<plg_t> plg_smooth(string rootname, const vector<plg_t> & polygons, float scale, const SGfield_t<float> & resolution,
//                            double factor, double dx, vector<float> & isovalues, float target, point2D_t point, int mode,
//                            bool tuned, bool enforcement, bool show_map, bool show_mesh, bool debug)

  vector<plg_t> plg_smooth(const string & rootname, const vector<plg_t> & polygons, plg_filter_t plg_filter, vector<float> & isovalues, float target, point2D_t point, int mode, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 filter a shoreline
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  float mask=-10;
  string filename;
  SGfield_t<float> field;
  vector<plg_t> contours, null;
  float *UGarray;
  int verbose=(debug==true);
  float optimal_isovalue=0.2;
  
  frame_t frame=plg_spherical_minmax(polygons);
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 compute smoothed landsea mask (structured)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(plg_filter.save_map) filename=rootname+"-mapping.nc";
  else                    filename="";
  
  plg_filter.check();
  
  field=Smoothedlandmask(filename, polygons, plg_filter, mode, optimal_isovalue, debug);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  optimize contours detection by masking areas where no contours could be
  expected
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  grid_t & grid=*field.grid;
  float *smoothed=field.x;
  
  float *masked=new float[grid.Hsize()];
  for(int m=0;m<grid.Hsize();m++) masked[m]=smoothed[m];

  int border=max(10,10);
  int depth=9;
 
#pragma omp parallel for
  for(int j=border;j<grid.ny-border;j++) {
    for(int i=border;i<grid.nx-border;i++) {
      int k=grid.nx*j+i;
      float value=smoothed[k];
      if( (value!=1) &&  (value!=-1) ) continue;
      int count=0;
/*------------------------------------------------------------------------------
      for a given value, see if neighbour cells show identical values*/
      for(int jj=j-depth;jj<j+depth+1;jj++) {
        for(int ii=i-depth;ii<i+depth+1;ii++) {
          int kk=grid.nx*jj+ii;
          if(kk==k) continue;
          if(smoothed[kk]==mask) continue;
          if(smoothed[kk]!=value) {
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
  delete[] smoothed;
  smoothed=masked;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  compute  iso-contours
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#------------------------------------------------------\n");
  printf("step 3: compute iso_contours\n");
  mesh_t mesh;
  
  if(target==-1) isovalues.push_back(optimal_isovalue);
  else isovalues.push_back(target);
  
  for(int k=0;k<isovalues.size();k++) {
    
    plg_destroy_entries(contours);
    contours.clear();
    
    STDOUT_BASE_LINE("treating iso-contour=%lf\n",isovalues[k]);
    
    string output;
    
/*------------------------------------------------------------------------------
    compute landmask isocontour, cartesian coordinates in input, geocentric in output */
    contours=plg_FindCountours(grid, smoothed, mask, isovalues[k], mesh, UGarray, plg_filter.save_mesh, debug, verbose);
    if(debug) {
      asprintf(output="",rootname+"-raw_contours%+04.2f.plg",isovalues[k]);
      status=plg_save(output.c_str(), PLG_FORMAT_SCAN, contours);
      }
    
    if(contours.size()==0) goto cleanUp;
    
/*------------------------------------------------------------------------------
    -180°/180° issue */
    status=plg_cartesian((projPJ) 0, contours);
    status=plg_recale(contours, frame, PLG_SPHERICAL);
    
/*------------------------------------------------------------------------------
    slightly smooth isocontour */
    plg_smoothlines(contours);
    status=plg_cartesian((projPJ) 0, contours);
    
    if(debug) {
      asprintf(output="",rootname+"-unfiltered_contours%+04.2f.plg",isovalues[k]);
      status=plg_save(output.c_str(), PLG_FORMAT_SCAN, contours);
      }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    make use of reference point to select final set
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    switch(mode) {
      case 0:
/*------------------------------------------------------------------------------
        shorelines limits */
        status=plg_RemoveImbrications(contours, point.x, point.y, PLG_POINT_INTERIOR);
        break;
      case 1:
/*------------------------------------------------------------------------------
        model limits */
        status=plg_RemoveInlandLimits(contours, point.x, point.y, (projPJ) 0, PLG_POINT_EXTERIOR, debug);
        break;
      case 2:
/*------------------------------------------------------------------------------
        keep both */
        status=0;
        break;
      }
    
    if(status!=0) {
      return(null);
      }
    
    asprintf(output="",rootname+"-raw_contours%+04.2f.plg",isovalues[k]);
    status=plg_save(output.c_str(), PLG_FORMAT_SCAN, contours);
    }

cleanUp:
  mesh.destroy();
  if(UGarray!=smoothed) delete[] UGarray;
  delete[] smoothed;
  grid.free();
  delete field.grid;
  
  return(contours);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_CreateLimits(const string & ShorelinesFile, const plg_t & polygon, string rootname, plg_filter_t & plg_filter, point2D_t point, double isocontour, bool smooth, bool removePools, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
  
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  int target;
  vector<plg_t> limits, smoothed, boundaries, raw;
  vector<plg_t> extended, extraction;
  string output;
  range_t<double> r;
  vector<float>  isovalues;

  if(isinf(point.x) or isnan(point.x)) TRAP_ERR_RETURN(-1,1," no reference position given, abandon\n");
  
  if(rootname=="") rootname="symtools";
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  shorelines extraction

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  copy spherical coordinates on cartesian coordinates */
  limits.push_back(polygon);
  status=plg_cartesian((projPJ) 0, limits);

/*------------------------------------------------------------------------------
  extend shorelines extraction zone */
  extended=plg_dilatation_cartesian(limits, 1.05, debug); // MERCATOR issue

/*------------------------------------------------------------------------------
  extract shorelines */
  printf("#################################################################\n");
  printf("extract shorelines from %s\n", ShorelinesFile.c_str());
  fflush_std();
  target=0;
  extraction=plg_extract(ShorelinesFile, extended, rootname, PLG_SPHERICAL, target, debug);
  status=plg_RemoveImbrications(extraction, point.x, point.y, PLG_POINT_INTERIOR);
  
  if(extraction.size()<1) {
    status=plg_save("extended.plg", PLG_FORMAT_SCAN, extended);
    TRAP_ERR_RETURN(-1,1," %d shorelines extracted, abandon\n",extraction.size());
    }
  
  plg_destroy_entries(extended);
  
  output=rootname+"-extracted-shorelines.plg";
  printf("#################################################################\n");
  printf("store extracted shorelines in %s\n",output.c_str());
  fflush_std();
  status=plg_save(output.c_str(), PLG_FORMAT_SCAN, extraction);
  
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  shorelines filtering

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  plg_filter.check();
  
  double dx=plg_filter.dx;
  double factor=plg_filter.factor;
  double lengthscale=plg_filter.lengthscale;
  
  printf("#################################################################\n");
  printf("filter shorelines: working grid resolution=%gm, lengthscale=%lfm or factor=%lfm\n",dx,lengthscale,factor);
  fflush_std();
  isovalues.clear();
  for(int k=0;k<7;k++) isovalues.push_back(0.2-k*0.1);
  
  int mode;
  if(removePools)
    mode=0;
  else
    mode=2;
  smoothed=plg_smooth(rootname, extraction, plg_filter, isovalues, isocontour, point, mode, debug);
    
  if(smoothed.size()==0) return(-1);
  
  output=rootname+"-filtered-shorelines.plg";
  printf("#################################################################\n");
  printf("store filtered shorelines in %s\n",output.c_str());
  fflush_std();
  status=plg_save(output.c_str(), PLG_FORMAT_SCAN, smoothed);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create model polygons

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  extend filtered extraction zone */
  extended=plg_dilatation_cartesian(limits, 1.02, debug); // MERCATOR issue
  
  printf("#################################################################\n");
  printf("extract filtered boundaries\n");
  fflush_std();
  target=1;
  boundaries=plg_extract(smoothed, extended, rootname, PLG_SPHERICAL, target, debug);
  
  plg_destroy_entries(extended);
  
  if(removePools){
    printf("#################################################################\n");
    STDOUT_BASE_LINE("remove swimming pools\n");
    TRAP_ERR_EXIT(ENOEXEC,"fail safe\n");
    fflush_std();
    status=plg_RemoveInlandLimits(boundaries, point.x, point.y, (projPJ) 0, PLG_POINT_EXTERIOR, debug);
    //status=plg_RemoveImbrications(boundaries, point.x, point.y, PLG_POINT_INTERIOR);
    }
  
  output=rootname+"-boundaries.plg";
  printf("#################################################################\n");
  printf("store boundaries in %s\n",output.c_str());
  fflush_std();
  status=plg_save(output.c_str(), PLG_FORMAT_SCAN, boundaries);
  
  plg_destroy_entries(smoothed);
  plg_destroy_entries(boundaries);
  
  return(status);
}

