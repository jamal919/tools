
/*******************************************************************************

  T-UGO tools, 2006-2015

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Create model limits.

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
#include "topo.h"

#include "polygons.hpp"
#include "exceptions.hpp"

using namespace Polygons;// for PolygonSet

class assembly_t {
private :
public :
  paire_t connecters[2];
  double  distance[2];
  assembly_t() {
    for(int k=0;k<2;k++) distance[k]=NAN;
    }
  bool mono() {
    bool chk=false;
    if(connecters[0].value[0]!=-1 and connecters[1].value[0]!=-1) {
      chk=(connecters[0].value[0]==connecters[1].value[0]);
      }
    }
};

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void plg_joind_and_cut (plg_t p, int m, plg_t q, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//   int backup,s,status;
//   double t,p,x,y;
//   plg_point_t point;
//   plg_target_t paire;
//   int *list;
//
// /**----------------------------------------------------------------------------
//   Exclude primary point from search database*/
//   x=gPolygones.p[move_paire.line].x[move_paire.point];
//   y=gPolygones.p[move_paire.line].y[move_paire.point];
//
//   gPolygones.p[move_paire.line].x[move_paire.point]=1.e+10;
//   gPolygones.p[move_paire.line].y[move_paire.point]=1.e+10;
//
//   paire=plg_clicked_point (*mouse_x,*mouse_y);
//
// /**----------------------------------------------------------------------------
//   Restore primary point*/
//   gPolygones.p[move_paire.line].x[move_paire.point]=x;
//   gPolygones.p[move_paire.line].y[move_paire.point]=y;
//
//   status= plg_colocalize(gPolygones.p, paire.line, paire.point,move_paire.line, move_paire.point);
//   if(move_paire.line==paire.line) {
//     if(paire.point>move_paire.point) {
//       list=plg_cut_polygone (&(gPolygones.p), &(gPolygones.n), paire.line, paire.point);
//       move_paire.line=list[0];
//       delete [] list;
//       list=plg_cut_polygone (&(gPolygones.p), &(gPolygones.n), move_paire.line, move_paire.point);
//       delete [] list;
//       }
//     else {
//       list=plg_cut_polygone (&(gPolygones.p), &(gPolygones.n), move_paire.line, move_paire.point);
//       paire.line=list[0];
//       delete [] list;
//       list=plg_cut_polygone (&(gPolygones.p), &(gPolygones.n), paire.line, paire.point);
//       delete [] list;
//       }
//     }
//   else {
//     list=plg_cut_polygone (&(gPolygones.p), &(gPolygones.n), move_paire.line, move_paire.point);
//     delete [] list;
//     list=plg_cut_polygone (&(gPolygones.p), &(gPolygones.n), paire.line, paire.point);
//     delete [] list;
//     }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int resolution_smooth(grid_t grid, float *buffer, float mask, range_t<float> range, vector<plg_t> & polygons, float slopemax)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  project buffer on a cartesian grid, potentially set mask from polygons,
  then smooth it
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  int count;
  grid_t cgrid;
  bool created=false;
  vector<plg_t> limits;
  
  double dx=grid.dx*110000.;
  signed char*selected=0;
  float *ctopo;
    
  if (polygons.size()==0) {
    frame_t frame(grid.xmin,grid.xmax,grid.ymin,grid.ymax);
    plg_t p(frame, PLG_SPHERICAL);
    polygons.push_back(p);
    created=true;
    }
  
  status=fe_defgrid(&cgrid, polygons, limits, dx);
  if(created) polygons.clear();
  
  
  cgrid.modeH=0;
  
  ctopo=new float[cgrid.Hsize()];

  int nprocs=initialize_OPENMP(-1);

#pragma omp parallel for private(status) if(nprocs>1)
  for(int m=0;m<cgrid.Hsize();m++) {
    double x, y, t, p;
    cgrid.xy(m,x,y);
    projection_to_geo(cgrid.projection, &p, &t, x, y);
    status=map_interpolation(grid, buffer, mask, t, p, &ctopo[m]);
    }
  
  selected=new signed char[cgrid.Hsize()];
  for(int j=0;j<cgrid.ny;j++) {
    for(int i=0;i<cgrid.nx;i++) {
      int m=j*cgrid.nx+i;
      selected[m]=1;
      }
    }

  if(polygons.size()!=0) {
    count=topo_PolygonsSelection(cgrid, polygons, selected, false, 0);
    if(count<0) return(-1);
    }
  
  if (range.width()!=0) {
    status=topo_RangeSelection(cgrid, ctopo, mask, range, selected, 0);
    }
  
  for(int m=0;m<cgrid.Hsize();m++) {
    if(selected[m]==0) ctopo[m]=mask;
    }
  
  status=smooth_density_new(cgrid, ctopo, mask, slopemax);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  projection back to spherical buffer
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
#pragma omp parallel for private(status) if(nprocs>1)
  for(int m=0;m<grid.Hsize();m++) {
    double x, y, t, p;
    float z;
    grid.xy(m,t,p);
    geo_to_projection(cgrid.projection, p, t, &x, &y);
    status=map_interpolation(cgrid, ctopo, mask, x, y, &z);
    if(status==0) buffer[m]=z;
    }
  
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int CheckSecant(vector<plg_t> & p, const vector<plg_t> & q, bool repair, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  check p and q intersection
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status=0;
  int npoints, annoying;
  point2D_t **points=0;
  double **angles=0;
  double *aindexes=0, *bindexes=0;
  vector<int> rejected;
  
// #pragma omp parallel for private(status,npoints)
  for(size_t l=0;l<p.size();l++) {
    if(p[l].npt==0) continue;
    for(size_t k=0;k<q.size();k++) {
      if(q[k].npt==0) continue;
      annoying=0;
      npoints=plg_secantpoints(p[l], q[k],  points, angles, &aindexes, &bindexes);
      if(npoints!=0) {
        if(repair) {
          if(debug) printf("\npolygons %4d and %4d have %2d points of intersection\n",l,k,npoints);
          for(int n=0;n<npoints;n++) {
            int ll=floor(aindexes[n]+0.5);
            int kk=floor(bindexes[n]+0.5);
/*------------------------------------------------------------------------------
            just extremeties colocated, not an issue */
            if( ( (kk==0) or (kk==q[k].npt-1) ) and ( (ll==0) or (ll==p[l].npt-1) )) {
              if(debug) printf("intersection %d : actually colocated extremities (ok)\n",n);
              continue;
              }
            else {
              if(debug) printf("intersection %d : truely secant\n",n);
              }
            if(debug) printf("position: %lf %lf, t=%lf p=%lf (closed=%d)\n",aindexes[n],bindexes[n], p[l].t[ll], p[l].p[ll], p[l].closed());
            annoying++;
            }
          if(debug and annoying!=0) printf("annoying intersections:%d\n",annoying);
          }
        }
/*------------------------------------------------------------------------------
      ??? */
      if(annoying==2) {
//         printf("%s intersection patch: suppress polygon %d (closed=%d)\n",__func__, l, p[l].closed());
//         p[l].destroy(); // HERE !!! supress points instead!
        int target;
        if(aindexes[0]<aindexes[1]) {
          target=ceil (aindexes[0]);
          }
        else {
          target=floor (aindexes[0]);
          }
        printf("%s intersection patch: suppress point %d in polygon %d (closed=%d)\n",__func__, target, l, p[l].closed());
        status=plg_delete_point(&p[l],target);
        l--;
        break;
        }
      rejected.push_back(l);

      }
    }
    
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_cleanAutoSecant(vector<plg_t> & polygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int repaired, status;
  double x,y;
  
// /*------------------------------------------------------------------------------
//   check angles in a given polygon */
//   for(size_t p=0;p<polygons.size();p++) {
//     
//     if(polygons[p].npt==0) continue;
//     
//     
//     for(size_t k=0;k<polygons[p].npt-1;k++) {
//       double angle=polygons[p].angle(k);
//       if(angle<0) {
//         printf("negative angle: %lfE %lfN\n", polygons[p].x[k], polygons[p].y[k]);
//         }
//       }
//     }
  
  vector<int> *flags;
  
  flags=new vector<int>[polygons.size()];

/*------------------------------------------------------------------------------
  check intersection in a given polygon */
  for(size_t p=0;p<polygons.size();p++) {
    redo:
    repaired=0;
    if(polygons[p].npt==0) continue;
    
    for(size_t k=0;k<polygons[p].npt-1;k++) {
      line_t a=line_t(polygons[p],k);
      for(size_t l=k+2;l<polygons[p].npt-1;l++) {
        if(fabs(k-l)<2) continue;
        line_t b=line_t(polygons[p],l);
        int flag=plg_secantpoint(a,b,&x,&y);
        if(flag==PLG_LINES_SECANT) {
          printf("polygon/segment %d is auto-secant, lines are %d %d (flag %c %c)\n",p,k,l,polygons[p].flag[k],polygons[p].flag[l]);
          printf("intersection position: %lfE %lfN\n",x,y);
          if(polygons[p].flag[k]!=polygons[p].flag[l]) {
/*------------------------------------------------------------------------------
            rigid/open limits crossing, patch the rigid section */
            int lineT, lineM, target;
            if(polygons[p].flag[k]=='T') {
              lineT=k;
              lineM=l;
              }
            if(polygons[p].flag[l]=='T') {
              lineT=l;
              lineM=k;
              }
//             int km,kp;
//             km=(lineT-2+polygons[p].npt-1) % (polygons[p].npt-1);
//             kp=(lineT+2) % (polygons[p].npt-1);
//             if     (lineM==km) target=lineT;
//             else if(lineM==kp) target=lineT+1;
//             else               target=lineT;
            if(lineT>lineM) target=lineT;
            else            target=lineT+1;
            printf("%s intersection patch: suppress point %d in polygon %d, position %lfE %lfN\n",__func__, target,p,polygons[p].t[target],polygons[p].p[target]);
            status=plg_delete_point(&polygons[p],target);
            flags[p].push_back(target);
            repaired++;
            }
          }
        }
      }
    if(repaired!=0) goto redo;
    }

  status=0;
  for(size_t p=0;p<polygons.size();p++) {
    if(flags[p].size()!=0) {
      status=-1;
      }
    }
  
  delete [] flags;
 
  if (status!=0) return(-1);
  
  flags=new vector<int>[polygons.size()];

/*------------------------------------------------------------------------------
  check duplicate points polygon */
#pragma omp parallel for
  for(size_t p=0;p<polygons.size();p++) {
    
    if(polygons[p].npt==0) continue;
    
    for(size_t k=0;k<polygons[p].npt-1;k++) {
      for(size_t l=k+1;l<polygons[p].npt-1;l++) {
        double d=plg_distance(polygons[p], k, l, PLG_CARTESIAN);
        int flag=(d==0.0);
        if(flag==1) {
          printf("polygon %d has duplicated points, %d %d\n",p,k,l);
          printf("duplicated position: %lfE %lfN\n",polygons[p].t[k],polygons[p].p[k]);
          flags[p].push_back(k);
          if(l==k+1) status=plg_delete_point(&polygons[p], l);
          break;
          }
        }
      }
    }

  status=0;
  for(size_t p=0;p<polygons.size();p++) {
    if(flags[p].size()!=0) {
      status=-1;
      }
    }
  
  delete [] flags;
  
  if (status!=0) return(-1);
  
  flags=new vector<int>[polygons.size()];

/*------------------------------------------------------------------------------
  check intersection with a another polygon */
#pragma omp parallel for
  for(size_t p=0;p<polygons.size();p++) {
    
    if(polygons[p].npt==0) continue;
    
    for(size_t k=0;k<polygons[p].npt-1;k++) {
      line_t a=line_t(polygons[p],k);
      for(size_t q=p+1;q<polygons.size();q++) {
    
        if(polygons[q].npt==0) continue;
    
        for(size_t l=0;l<polygons[q].npt-1;l++) {
          line_t b=line_t(polygons[q],l);
          int flag=plg_secantpoint(a,b,&x,&y);
          if(flag==PLG_LINES_SECANT) {
            printf("polygon/segment %d is secant with %d, lines are %d %d\n",p,q,k,l);
            printf("intersection position: %lfE %lfN\n",x,y);
            flags[p].push_back(q);
//             if()
            break;
            }
          }
        }
      }
    }
  
  status=0;
  for(size_t p=0;p<polygons.size();p++) {
    if(flags[p].size()!=0) {
      status=-1;
      }
    }

  delete [] flags;
    
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t> plg_assembly_limits(vector<plg_t> & p, vector<plg_t> & q, projPJ projection, string & rootname, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  join p polygons (typically rigid limits) on q extremities (typically open limits)

  Special case: split bounary does not cross shorelines (1 open segemnt), to
                be further checked
                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  double d;
  int k,n;
  paire_t paire;
  plg_point_t point;
  int list;
  string filename;
  mesh_t mesh;
  bool repair, strict;
  vector<int> unclosed;
  vector<paire_t> connected;
  vector<plg_t> r;

  if(rootname=="") rootname="anonymous";
      
/*------------------------------------------------------------------------------
  insure mono-block connecting polygons */
  repair=true;
  printf("check constructed set\n");
  status=plg_CheckClosure(p, repair, 0, debug);
  printf("check imported set\n");
  status=plg_CheckClosure(q, repair, 0, debug);
    
/*------------------------------------------------------------------------------
  -180°/180° longitude issue */
  double t0=p[0].t[0];
  
  q[0].t[0]=degree_recale(q[0].t[0],t0);
  
  frame_t frame=plg_recale(q, t0);
  
  for(int k=0; k<p.size();k++)
    if(p[k].flag==0) p[k].SetFlag('T');
  
  for(int l=0; l<q.size();l++)
    if(q[l].flag==0) q[l].SetFlag('M');
  
  for(int k=0; k<p.size();k++) if(!p[k].closed()) unclosed.push_back(k);
//   printf("%d unclosed rigid polygons found\n",unclosed.size());
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  scan open limits polygons, identify and cut rigid polygon to insure common
  extremities
  
  step 1: following joining extremities distance, either
    - move rigid (possiblly internal) point to the open limit one (first/last)
    - add a new point in open limt at start/end
    
  This a recent change (06/10/2018), can be reset to previous behaviour by setting
  ratio to zero
    
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  double ratio=0.5;
  
  for(int l=0; l<q.size();l++) {
    double distance, length;
/*------------------------------------------------------------------------------
    first (open limit) q point*/
    point=plg_point_t(q[l],0);
    paire=plg_find_polygon(p, unclosed, point, d);
    k=paire.value[0];
    n=paire.value[1];
    if(n!=-1) {
      distance= plg_distance(p[k],n, q[l], 0, PLG_CARTESIAN);
      length=plg_minlength(p[k], PLG_CARTESIAN);
      printf("%s : %d %d first point at %lf (min=%lf)\n",__func__,k,l,distance,length);
      if(distance<length*ratio) {
        status=plg_move_point(p[k], n, point);
        }
      else {
        point=plg_point_t(p[k],n);
        status=plg_insert_point(q[l],0,point);
        }
      connected.push_back(paire);
      }
/*------------------------------------------------------------------------------
    last (open limit) q point*/
    point=plg_point_t(q[l],q[l].npt-1);
    paire=plg_find_polygon(p, unclosed, point,d);
    k=paire.value[0];
    n=paire.value[1];
    if(n!=-1) {
      distance= plg_distance(p[k],n, q[l], q[l].npt-1, PLG_CARTESIAN);
      length=plg_minlength(p[k], PLG_CARTESIAN);
      printf("%s : %d %d last point at %lf (min=%lf)\n",__func__,k,l,distance,length);
      if(distance<length*ratio) {
        status=plg_move_point(p[k], n, point);
        }
      else {
        point=plg_point_t(p[k],n);
        status=plg_insert_point(q[l],q[l].npt,point);
        }
      connected.push_back(paire);
      }
    }

  int *count=new int[p.size()];
  for(int l=0; l<p.size();l++) count[l]=0;
  
  for(int l=0; l<connected.size();l++) {
    count[connected[l].value[0]]++;
    }

  for(int l=0; l<connected.size();l++) {
    if(count[connected[l].value[0]]==2) continue;
    printf("single-point connection found with imported open limits : segment %d point %d\n",connected[l].value[0],connected[l].value[1]);
    }
  
  for(int l=0; l<p.size();l++) {
    if(count[l]!=1) continue;
    k=(l+1) % p.size();
    if(count[k]!=1) continue;
    printf("connection found %d %d\n", l, k);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  prepare p set for assembly
  
  pevious version failed if imported open limit is closed 
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(int l=0; l<q.size();l++) {
/*------------------------------------------------------------------------------
    first q point*/
    point=plg_point_t(q[l],0);
//     paire=plg_find_polygon(p, point, d);
    paire=plg_find_polygon(p, unclosed, point, d);
    k=paire.value[0];
    n=paire.value[1];
    if(n!=-1) {
      status=plg_move_point(p[k], n, point);
      list= plg_cut_polygone(p, k, n);
      }
/*------------------------------------------------------------------------------
    last q point*/
    point=plg_point_t(q[l],q[l].npt-1);
//     paire=plg_find_polygon(p, point,d);
    paire=plg_find_polygon(p, unclosed, point, d);
    k=paire.value[0];
    n=paire.value[1];
    if(n!=-1) {
      status=plg_move_point(p[k], n, point);
      list= plg_cut_polygone(p, k, n);
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  check and repair possible intersection between constructed and imported set
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(debug) {
    vector<plg_t> tmp;
    tmp=plg_duplicate(q);
    status=plg_merge(tmp,p);
    filename=rootname+"-merged-unfiltered.plg";
    status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, tmp);
    }
  
  repair=true;
  CheckSecant(p, q, repair, debug);
    
  status=plg_merge(q,p);
  
  if(debug) {
    filename=rootname+"-merged.plg";
    status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, q);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  assembly limit set by using mesh cycles detection
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  create degenerated mesh from polygons                                       */
  strict=true;
  bool do_edges=true;
  status=plg_export2nei(q, mesh, strict, do_edges, debug);
  
/*------------------------------------------------------------------------------
  then get elementary cycles                                                  */
  bool check_islands=true;
  int verbose=0;
  r=Nei2Polygons (mesh, check_islands, verbose, debug);
  
  if(debug) {
    filename=rootname+"-elementary-cycles.plg";
    status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, r);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  filter cycles
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(int s=0;s<r.size();s++){
    double area=r[s].area();
    if(area<=0.) {
      if(debug) printf("polygon %d eliminated, negative area=%lf\n",s, area);
      r.erase(r.begin()+s);
      s--;
      continue;
      }
    if(plg_isclosed(r[s])==0) {
      if(debug) printf("polygon %d eliminated, unclosed: area=%lf\n",s, area);
      r.erase(r.begin()+s);
      s--;
      }
    }
  
/*------------------------------------------------------------------------------
  -180°/180° longitude issue */
  frame=plg_recale(r, r[0].t[0]);
    
  status=plg_cartesian(projection, r);

  status=plg_cleanAutoSecant(r);
  
  if(debug) {
    filename=rootname+"-elementary-cycles-cleaned.plg";
    status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, r);
    }
  
  for(int s=0;s<r.size();s++){
    for(int k=0;k<r[s].npt-1;k++) {
      if(r[s].flag[k]=='O') {
        r[s].flag[k]='M';
        }
      }
    }
 
  return(r);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_CheckConnection(plg_t & pp, int end, plg_t & qq, int mode, int & e)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double distance;
  
  if(pp.npt==0) return(-1);
  if(qq.npt==0) return(-1);
  
  if(pp.closed()) return(-1);
  if(qq.closed()) return(-1);
  
  if(end==0) {
    distance= plg_distance(pp, 0, qq, 0, mode);
    if(distance==0.0) {
      e=0;
      return(0);
      }
    distance= plg_distance(pp, 0, qq, qq.npt-1, mode);
    if(distance==0.0) {
      e=1;
      return(0);
      }
    }
  else {
    distance= plg_distance(pp, pp.npt-1, qq, 0, mode);
    if(distance==0.0) {
      e=0;
      return(0);
      }
    distance= plg_distance(pp, pp.npt-1, qq, qq.npt-1, mode);
    if(distance==0.0) {
      e=1;
      return(0);
      }
    }

  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<int> plg_FindPath(vector<plg_t> & p, paire_t p1, paire_t p2, int mode, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int end, status;
  char *used;
  int current,flag=-1,count;
  vector<int> path;
  
  used=new char[p.size()];
  
  for (int k=0; k< p.size();k++) used[k]=0;
  
  current=p1.value[0];
  
  if(fabs(p[current].npt-1-p1.value[1]) > p1.value[1])
    end=0;
  else
    end=1;
  
  end=1-end;
  
  path.push_back(current);
  
  count=0;
  do {
    for (int k=0; k< p.size();k++) {
      count++;
      if(k==current) continue;
      if(p[k].closed()) continue;
      status=plg_CheckConnection(p[current], end, p[k], mode, flag);
      if(status==0) {
        current=k;
        path.push_back(current);
        end=1-flag;
        break;
        }
      }
    if(count > p.size()) {
      path.clear();
      break;
      }
    } while (current!=p2.value[0]);
    
  if(debug) printf("path size=%d\n",path.size());

  return(path);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_GetProxy(vector<plg_t> & p, vector<plg_t> & q, double dmax, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  paire_t paire[2];
  plg_point_t point;
  double d[2];
  vector<int> found;
  
/*------------------------------------------------------------------------------
  temporary neutralising the rigid boundaries */
  for(int s=0;s<q.size();s++) {
    if(q[s].flag[0]=='T') q[s].npt*=-1;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  identify corresponding limits in both sets
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(int l=0; l<p.size();l++) {
/*------------------------------------------------------------------------------
    first p point*/
    point=plg_point_t(p[l],0);
    paire[0]=plg_find_polygon(q, point, d[0]);
/*------------------------------------------------------------------------------
    last p point*/
    point=plg_point_t(p[l],p[l].npt-1);
    paire[1]=plg_find_polygon(q, point, d[1]);
    if(paire[0].value[0]==paire[1].value[0] and max(d[0],d[1])<dmax) {
      found.push_back(paire[0].value[0]);
      }
    }

/*------------------------------------------------------------------------------
  restore rigid boundaries */
  for(int s=0;s<q.size();s++) {
    if(q[s].flag[0]=='T') q[s].npt*=-1;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  add complementary open boundary
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(int s=0;s<q.size();s++) {
    if(q[s].flag[0]=='T') continue;
    int pos=vpos(s,found);
    if(pos!=-1) {
      printf("open segment %d added to imported open boundaries",s);
      plg_t tmp;
      tmp.duplicate(q[s]);
      tmp.SetFlag('M');
      p.push_back(tmp);
      }
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_GetCompanions(vector<plg_t> & p, vector<plg_t> & q, char flag, vector<assembly_t> & connections, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  plg_point_t point;
  double d[2];
  
/*------------------------------------------------------------------------------
  temporary neutralising boundaries with different flag*/
  for(int s=0;s<q.size();s++) {
    if(q[s].flag[0]!=flag) q[s].npt*=-1;
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  identify p set coresponding limits in q set
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  connections.clear();

  for(int l=0; l<p.size();l++) {
    assembly_t *tmp=new assembly_t;
/*------------------------------------------------------------------------------
    first p point*/
    point=plg_point_t(p[l],0);
    tmp->connecters[0]=plg_find_polygon(q, point, d[0]);
    tmp->distance[0]=d[0];
/*------------------------------------------------------------------------------
    last p point*/
    point=plg_point_t(p[l],p[l].npt-1);
    tmp->connecters[1]=plg_find_polygon(q, point, d[1]);
    tmp->distance[1]=d[1];
    connections.push_back(*tmp);
    }

/*------------------------------------------------------------------------------
  restore neutralized boundaries */
  for(int s=0;s<q.size();s++) {
    if(q[s].flag[0]!=flag) q[s].npt*=-1;
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_GetComplementary(vector<plg_t> & p, vector<plg_t> & q, projPJ projection, double dmax, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  treat the case where imported open boundaries are only part of full open set
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int n=-1, status;
  bool equalize=false;
  vector<assembly_t> connections, patching;
  vector<int> issues;
  
  status=plg_GetCompanions(p, q, 'T',  connections, debug);
  for(int l=0; l<p.size();l++) {
    if(connections[l].distance[0]>100.0 or connections[l].distance[1]>100.0) {
      int s,k;
      double d;
      s=connections[l].connecters[0].value[0];
      k=connections[l].connecters[0].value[1];
      d=connections[l].distance[0];
      printf("imported %d first point: connexion  segment=%4d point=%6d (%lf km)\n",l,s,k,d);
      s=connections[l].connecters[1].value[0];
      k=connections[l].connecters[1].value[1];
      d=connections[l].distance[1];
      printf("imported %d last point:  connexion  segment=%4d point=%6d (%lf km)\n",l,s,k,d);
      issues.push_back(l);
      }
    }
  
  if(issues.size()==0) return(0);
     
// /*------------------------------------------------------------------------------
//   neutralize rigid boundaries */
//   for(int s=0;s<q.size();s++) {
//     if(q[s].flag[0]=='T') q[s].npt*=-1;
//     }
    
  vector<plg_t> tmp;
  tmp=plg_subdivide(q,1000.0,equalize, PLG_CARTESIAN,debug);
    
  status=plg_spherical(projection, tmp);
  status=plg_GetCompanions(p, tmp, 'M',  patching, debug);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  identify coressponding limits in both sets
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(int ll=0; ll<issues.size();ll++) {
    int l=issues[ll];
    int s,k,r;
    double d, threshold=1.0e-04;
    if(not patching[l].mono()) continue;
/*------------------------------------------------------------------------------
    p point*/
    if(connections[l].distance[0]>100.0) {
      s=patching[l].connecters[0].value[0];
      k=patching[l].connecters[0].value[1];
      d=patching[l].distance[0];
      printf("cut %d first point: connexion  segment=%4d point=%6d (%lf km)\n",l,s,k,d);
      r=plg_cut_polygone(tmp, k, n);
      status=plg_decimate(tmp[s], threshold, debug);
      status=plg_decimate(tmp[r], threshold, debug);
      tmp[s].SetFlag('O');
      tmp[r].SetFlag('M');
      q[s].destroy();
      q.push_back(tmp[s]);
      q.push_back(tmp[r]);
      }
    if(connections[l].distance[1]>100.0) {
      s=patching[l].connecters[1].value[0];
      k=patching[l].connecters[1].value[1];
      d=patching[l].distance[1];
      printf("cut %d last point: connexion  segment=%4d point=%6d (%lf km)\n",l,s,k,d);
      r=plg_cut_polygone(tmp, s, k);
      status=plg_decimate(tmp[s], threshold, debug);
      status=plg_decimate(tmp[r], threshold, debug);
      tmp[s].SetFlag('O');
      tmp[r].SetFlag('M');
      q[s].destroy();
      q.push_back(tmp[s]);
      q.push_back(tmp[r]);
      }
    }
  
/*------------------------------------------------------------------------------
  remove empty polygons */
  status=plg_compact_entries(q);

  status=plg_save("debug.plg", PLG_FORMAT_SCAN, q);
 
// /*------------------------------------------------------------------------------
//   restore rigid boundaries */
//   for(int s=0;s<q.size();s++) {
//     if(q[s].flag[0]=='T') q[s].npt*=-1;
//     }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_CheckOpenConsistency(string rootname, const vector<plg_t> & boundaries, vector<plg_t> & imported, vector<plg_t> & selected, projPJ projection, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  parse imported open boundaries with the result of initial limits construction
 
  highly empirical, might fail in many conditions !!!
  
  initially built to have all expected open limits from imported set.
  
  patch #0 implemented to add missing open limits
  
  patch #1 implemented to deal with...
  
  patch #2 generalisation of #1 ?
  
  Special case: split bounary does not cross shorelines (1 open segemnt), to
                be further checked
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int  k,n,status;
  vector<plg_t> splitted, opened;
  string OpensName=rootname+"-open.plg";
  string ImportsName=rootname+"-imported-mofified.plg";
  paire_t paire;
  plg_point_t point;
  vector<paire_t> connected;
  vector<int> bridged;
  double d;
  int nissues=0;
  int verbose=0;
  
/*------------------------------------------------------------------------------
  split rigid and open segments in separate polygons */
  splitted=plg_split(boundaries, debug);
  
/*------------------------------------------------------------------------------
  identify open pieces and store for further assembly if required */
  for(int s=0;s<splitted.size();s++) {
    if(splitted[s].flag[0]=='M') {
      plg_t *tmp=new plg_t;
      tmp->duplicate(splitted[s]);
      opened.push_back(*tmp);
      }
    }
  if(debug) status=plg_save(OpensName.c_str(), PLG_FORMAT_SCAN, opened);

/*------------------------------------------------------------------------------
  recent extension */
  if(opened.size()>imported.size()) {
    printf("\nWARNING : actual (mesh-limits) opened limits %d, imported open limits %d\n", opened.size(), imported.size());
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    patch #0 : try to add necessary open boundaries to import
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    double dmax=1.0;
    printf("try patch #0-a (add some missing open limits from constructed set\n");
    status=plg_GetProxy(imported, splitted, dmax, debug);
    if(debug) status=plg_save(ImportsName.c_str(), PLG_FORMAT_SCAN, imported);
    }
  else if(opened.size()<imported.size()) {
    printf("\nWARNING : actual (mesh-limits) opened limits %d, imported open limits %d\n", opened.size(), imported.size());
    printf("this could make the assembled boundaries unworkable, check your selection polygon\n\n");
    if(debug) status=plg_save(OpensName.c_str(), PLG_FORMAT_SCAN, opened);
    return(-1);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    patch #0-bis : check if imported is only a sub-section of generated limits
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(opened.size()==imported.size()) {
    printf("try patch #0-b (add parts of open limits from constructed set\n");
    bool partial=false;
    for(int s=0;s<opened.size();s++) {
      double length[2], ratio;
      length[0]=imported[s].cartesian_size()/1000.;
      length[1]=opened[s].cartesian_size()/1000.;
      printf("lengths %d: imported=%fkm constructed=%fkm\n", s, length[0], length[1]);
      ratio=0.5*(length[0]+length[1])/sqrt(length[0]*length[1]);
      if(ratio<0.5 or ratio>1.5) partial=true;
      }
    if(partial) status=plg_GetComplementary(imported, splitted, projection, 10.0, debug);
    }
  
  
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  temporary neutralising the open boundaries
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
/*------------------------------------------------------------------------------
  neutralise open boundaries                                                  */
  for(int s=0;s<splitted.size();s++) {
    if(splitted[s].npt==0) continue;
    if(splitted[s].flag[0]=='M') {
      splitted[s].npt*=-1;
      }
    }
  
/*------------------------------------------------------------------------------
  neutralise closed rigid boundaries                                                  */
  for(int s=0;s<splitted.size();s++) {
    if(splitted[s].npt==0) continue;
    if(splitted[s].flag[0]=='T' and splitted[s].closed()) {
      splitted[s].npt*=-1;
      }
    }
  
  bool repair=true;
  status=plg_CheckClosure(splitted, repair, verbose, debug);
  status=plg_save("debug.plg", PLG_FORMAT_SCAN, splitted);
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  find connexion points between imported open boundaries and working rigid
  boundaries. First check if imported limits connect to rigid limits
  
  each imported limit should have 2 connections (first/last points)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  debug=true;
/*------------------------------------------------------------------------------
  scan imported limits connections with rigid limits                          */
  for(int l=0; l<imported.size();l++) {
/*------------------------------------------------------------------------------
    first q point*/
    point=plg_point_t(imported[l],0);
    paire=plg_find_polygon(splitted, point, d);
    k=paire.value[0];
    n=paire.value[1];
    if(n!=-1) {
      connected.push_back(paire);
      if(debug) printf("open %d first point: connexion found, segment=%4d point=%6d (%lf km)\n",l,paire.value[0],paire.value[1],d);
      }
/*------------------------------------------------------------------------------
    last q point*/
    point=plg_point_t(imported[l],imported[l].npt-1);
    paire=plg_find_polygon(splitted, point,d);
    k=paire.value[0];
    n=paire.value[1];
    if(n!=-1) {
      connected.push_back(paire);
      if(debug) printf("open %d  last point: connexion found, segment=%4d point=%6d (%lf km)\n",l,paire.value[0],paire.value[1],d);
      }
    }
  
/*------------------------------------------------------------------------------
  re-activate neutralized boundaries */
  for(int s=0;s<splitted.size();s++) {
    if(splitted[s].npt==0) continue;
//     if(splitted[s].flag[0]=='M') splitted[s].npt*=-1;
    if(splitted[s].npt<0) splitted[s].npt*=-1;
    }

/*------------------------------------------------------------------------------
  select all rigid boundaries */
  for(int s=0;s<splitted.size();s++) {
    if(splitted[s].npt==0) continue;
    if(splitted[s].flag[0]!='M') {
      plg_t *tmp=new plg_t;
      tmp->duplicate(splitted[s]);
      selected.push_back(*tmp);
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  count number of paires found for each rigid limits
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  int *count=new int[splitted.size()];
  for(int l=0; l<splitted.size();l++) count[l]=0;
  
  for(int l=0;l<connected.size();l++) {
/*------------------------------------------------------------------------------
    selected limit index */
    int k=connected[l].value[0];
    count[k]++;
    }

  nissues=0;
  for(int l=0;l<connected.size();l++) {
    if(count[connected[l].value[0]]==2) continue;
    if(debug) printf("hole found, segment=%4d point=%6d\n",connected[l].value[0],connected[l].value[1]);
    nissues++;
    }
  if(nissues==0) goto clear;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  patch #1 : try imported polygons as bridging open limit between to rigid ones
  
  symptom : a rigid limit has one connection only with the imported open limits.
            could be that a new open rigid has inserted inside an imported open
            limit (because of more detailed digitalizing)
            
            -----XooooooooooooooooooooooooooooooooooX-------------------
            
                2                                  1        1
                !                                  !
            ----XooooooooooooooooooooooooooooooooooX----XoooX---------
  
  see if 2nd next (i.e. next rigid) has also one point connection too
  
  confidence level quite low on this patch
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("try patch #1 (currently %d issues)\n", nissues);
  for(int l=0; l<splitted.size();l++) {
    if(count[l]!=1) continue;
    k=(l+2) % splitted.size();
    if(count[k]!=1) continue;
    int target=(l+1) % splitted.size();
    printf("bridge found between coastal segments %d %d: using opened segment %d (attemptive patch will apply)\n", l, k, target);
    for(int kk=0;kk<splitted[target].npt-1;kk++) splitted[target].flag[kk]='S';
    plg_t *tmp=new plg_t;
    tmp->duplicate(splitted[target]);
    selected.push_back(*tmp);
    count[l]=0;
    count[k]=0;
    }
  
  nissues=0;
  for(int l=0; l<splitted.size();l++) {
    if(count[l]==0) continue;
    if(count[l]==2) continue;
    printf("patch #1: residual problem at segment %d (%d)\n", l, count[l]);
    nissues++;
    }
  if(nissues==0) goto end;
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  patch #2 : idem but try a more complex path
  
  confidence level quite low on this patch
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("try patch #2 (currently %d issues)\n", nissues);
  for(int l=0; l<connected.size();l++) {
    if (count[connected[l].value[0]]!=1) continue;
    for(int k=0; k<connected.size();k++) {
      if (count[connected[k].value[0]]!=1) continue;
      if( connected[l].value[0] >= connected[k].value[0] ) continue;
      if(debug) printf("try bridge between coastal segments %d %d\n", connected[l].value[0], connected[k].value[0]);
      vector<int> path=plg_FindPath(splitted, connected[l], connected[k], PLG_SPHERICAL, debug);
      if(path.size()!=0) {
        printf("path found between %d %d (attemptive patch will apply)\n", connected[l].value[0], connected[k].value[0]);
        for (int ll=0 ; ll < path.size(); ll++) {
          int target=path[ll];
          if(splitted[target].flag[0]=='M') {
            for (int kk=0;kk<splitted[target].npt-1;kk++) splitted[target].flag[kk]='S';
            printf("open segment %d set as rigid segment\n", target);
            plg_t *tmp=new plg_t;
            tmp->duplicate(splitted[target]);
            selected.push_back(*tmp);
//             rigid.push_back(splitted[target]);
            count[connected[l].value[0]]=0;
            count[connected[k].value[0]]=0;
            }
          }
        }
      }
    }
  
  nissues=0;
  for(int l=0; l<splitted.size();l++) {
    if(count[l]==0) continue;
    if(count[l]==2) continue;
    printf("patch #2 residual problem at segment %d (%d)\n", l, count[l]);
    nissues++;
    }
  if(nissues==0) goto end;
  
// failed:
  status=plg_save("failed.plg", PLG_FORMAT_SCAN, imported);
  printf("\nrigid set is NOT safe for assembly with imported open limtis, assembly procedure might fail\n\n");
  for(int s=0;s<splitted.size();s++) splitted[s].destroy();
  splitted.clear();
  for(int s=0;s<opened.size();s++) opened[s].destroy();
  opened.clear();
  return(-1);

end:
  printf("\nrigid set is safe, all issues have been fixed\n\n");
  for(int s=0;s<splitted.size();s++) splitted[s].destroy();
  splitted.clear();
  for(int s=0;s<opened.size();s++) opened[s].destroy();
  opened.clear();
  return(1);
  
clear:
  printf("\nrigid set is safe, no issues found\n\n");
  for(int s=0;s<splitted.size();s++) splitted[s].destroy();
  splitted.clear();
  for(int s=0;s<opened.size();s++) opened[s].destroy();
  opened.clear();
  return(0);
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
    TRAP_ERR_EXIT(-1,"exiting\n");
    }

  return(fmt);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int check_limits(vector<plg_t> p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=plg_SortBySize(p,0);
  
  for(int k=1; k<p.size();k++) {
    point2D_t point=p[k].InsidePoint();
    int check=0;
    int flag=0;
    flag=plg_single(p[0], point.x, point.y, &flag, PLG_SPHERICAL, check);
    if(flag!=PLG_POINT_INTERIOR) {
      printf("polygon %d not include in polygon %d\n",k,0);
      return(-1);
      }
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int adjust_assembly(vector<plg_t> & p, const vector<plg_t> & circles, const vector<plg_t> & centres, const vector<double> & sizes, projPJ projection, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  orignally implemented to adapt fine rigid reolution to possibly larger open
  limits
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i, status;
  
  debug=true;
  
  if(debug) {
    printf("%s sizes:", __func__);
    for(int k=0;k<sizes.size(); k++) printf(" %lf",sizes[k]);
    printf("\n");
    }
  
//   status=plg_SortBySize(p,0);
  
  for(int k=0; k<p.size();k++) {
    if(p[k].flag[0]!='T') continue;
    int check=0;
    point2D_t point, test;
/*------------------------------------------------------------------------------
    check first rigid polygon point    */
    point=point2D_t(p[k].t[0], p[k].p[0]);
    for(int l=0; l<centres.size();l++) {
      int flag=0;
      flag=plg_single(centres[l], point.x, point.y, &flag, PLG_SPHERICAL, check);
      if(flag==PLG_POINT_INTERIOR) {
/*------------------------------------------------------------------------------
        should tune circles/sizes with respect to rigid segment resolution    */
        double L=plg_length(p[k], 0, 1, PLG_SPHERICAL)*1000.;
//         status=plg_circles(opened, projection, 1.5, circles, sizes);
        i=0;
        while(flag==PLG_POINT_INTERIOR) {
          if(i==p[k].npt) {i--;break;}
          test=point2D_t(p[k].t[i], p[k].p[i]);
          flag=0;
          flag=plg_single(circles[l], test.x, test.y, &flag, PLG_SPHERICAL, check);
          i++;
          }
        printf("polygon %d, points from 0 to %d included in circle %d\n",k,i,l);
/*------------------------------------------------------------------------------
        split p into two parts, first one to be re-sampled                    */
        i=min(i+1, p[k].npt-1);
        plg_t q1;
        q1.duplicate(p[k],0,i);
        plg_t q2;
        q2.duplicate(p[k],i,p[k].npt-1);
        status=plg_flip(q1);
        vector<double> prescribed;
        double radius=sizes[l];
        L=plg_length(q1, 0, 1, PLG_SPHERICAL)*1000.;
        if(L<radius) {
          prescribed.push_back(L);
          while(L<radius) {
//             L=L*1.333;
            L=L*1.25;
            prescribed.push_back(L);
            }
          }
        else {
          prescribed.push_back(L);
          while(L>radius) {
//             L=L/1.333;
            L=L/1.25;
            prescribed.push_back(L);
            }
          }
        int equalize=1;
        plg_t resampled=plg_resample(q1, prescribed, equalize, debug);
        resampled.SetFlag('T');
        status=plg_spherical(projection, &resampled, 1);
        status=plg_flip(resampled);
/*------------------------------------------------------------------------------
        projection got involved, accuracy issue may appear at extremities; fix that */
        plg_copy_point(&resampled, 0, q1, q1.npt-1);
        plg_copy_point(&resampled, resampled.npt-1, q1, 0);
        status=plg_concat(resampled,q2);
        p[k].destroy();
        p[k]=resampled;
        }
      }
/*------------------------------------------------------------------------------
    check last rigid polygon point    */
    point=point2D_t(p[k].t[p[k].npt-1], p[k].p[p[k].npt-1]);
    for(int l=0; l<centres.size();l++) {
      int flag=0;
      flag=plg_single(centres[l], point.x, point.y, &flag, PLG_SPHERICAL, check);
      if(flag==PLG_POINT_INTERIOR) {
        i=p[k].npt-1;
        while(flag==PLG_POINT_INTERIOR) {
          if(i==-1) {i++;break;}
          test=point2D_t(p[k].t[i], p[k].p[i]);
          flag=0;
          flag=plg_single(circles[l], test.x, test.y, &flag, PLG_SPHERICAL, check);
          i--;
          }
        printf("polygon %d, points from %d to %d included in circle %d\n",k,i,p[k].npt-1,l);
/*------------------------------------------------------------------------------
        split p into two parts, second one to be re-sampled                   */
        i=max(i-1, 0);
        plg_t q1;
        q1.duplicate(p[k],0,i);
        plg_t q2;
        q2.duplicate(p[k],i,p[k].npt-1);
        vector<double> prescribed;
        double L=plg_length(q2, 0, 1, PLG_SPHERICAL)*1000.;
        double radius=sizes[l];
        L=plg_length(q1, 0, 1, PLG_SPHERICAL)*1000.;
        if(L<radius) {
          prescribed.push_back(L);
          while(L<radius) {
//             L=L*1.333;
            L=L*1.25;
            prescribed.push_back(L);
            }
          }
        else {
          prescribed.push_back(L);
          while(L>radius) {
//             L=L/1.333;
            L=L/1.25;
            prescribed.push_back(L);
            }
          }
        int equalize=1;
        plg_t resampled=plg_resample(q2, prescribed, equalize, debug);
        resampled.SetFlag('T');
        status=plg_spherical(projection, &resampled, 1);
/*------------------------------------------------------------------------------
        projection got involved, accuracy issue may appear at extremities; fix that */
        plg_copy_point(&resampled, 0, q2, 0);
        plg_copy_point(&resampled, resampled.npt-1, q2, q2.npt-1);
        status=plg_concat(q1, resampled);
        p[k].destroy();
        p[k]=q1;
        }
      }
    }
  
//   for(int k=0; k<p.size();k++) {
//     int check=0;
//     point2D_t point=point2D_t(p[k].t[0], p[k].p[0]);
//     for(int l=0; l<circles.size();l++) {
//       int flag=0;
//       flag=plg_single(circles[l], point.x, point.y, &flag, PLG_SPHERICAL, check);
//       if(flag==PLG_POINT_INTERIOR) {
//         printf("polygon %d include in circle %d\n",k,l);
//         }
//       }
//     point=point2D_t(p[k].t[p[k].npt-1], p[k].p[p[k].npt-1]);
//     for(int l=0; l<circles.size();l++) {
//       int flag=0;
//       flag=plg_single(circles[l], point.x, point.y, &flag, PLG_SPHERICAL, check);
//       if(flag==PLG_POINT_INTERIOR) {
//         printf("polygon %d include in circle %d\n",k,l);
//         }
//       }
//     }

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
    "  %s [OPTIONS] input\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Create model limits.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  Show this help and exit.\n"
    "  -a  followed by the path of a pre-existing open limits file (aimed to work with mesh-split/mesh-assembly).\n"
    "  -b  followed by the path of a pre-processed limits file (aimed to operate on post-processing section only). Disables --extract, --smooth, and --limits\n"
    "  -c  followed by path to template criteria file\n"
    "  -p  followed by path to domain polygon\n"
    "  -o  followed by the path of the output file. Default from the input file and output format if both are given.\n"
    "  -s  followed by the path of the raw shoreline file\n"
    "  -z  followed by the zone. Will override --point\n"
    "  -I  followed by the input format. Default from the extension of the input file.\n"
    "  -O  followed by the output format. Default from the extension of the output file.\n"
    "  -r  followed by the background smoothing grid resolution (default is 1/4 of -l, i.e. 625m when -l and --sampling are at their defaults).\n"
    "  -l  followed by the smoothing length in m (default is 1/2 of --sampling, i.e. 2500m when --sampling is at its default).\n"
    "  --frame ...  delimiting rectangular frame\n"
    "  --extract  extract shorelines. Disabled by -b\n"
    "  --sampling  resolution of the ouput in KM. Default: 5km. Now obsolete, please use --sampling_length instead\n"
    "  --sampling_length  resolution of the ouput in METERS. Default: 5000m\n"
    "  --rootname  root name of the MANY intermediate output files. Default: anonymous\n"
    "  --point  followed by [lon;lat] of any ocean point inside the domain\n"
    "  --resolution  followed by the paths of polygons (.nei format) and values files concatenated by a space!\n"
    "      typical option argument: --resolution '../resolution.nei ../resolution.dat taht contains <lon lat value> setting'\n"
    "  --debug  enable debug outputs\n"
    "  --smooth  smooth shoreline. Disabled by -b\n"
    "  --limits  build model limits. Disabled by -b\n"
    "  --adjust  in case of assembly with pre-existing open limits, will adjust coastal resolution to fit neighbouring open limits\n"
    "  --lake when not open sea shorelines\n"
    "  --isocontour followed by targeted isovalue (optimal:-1, default:-0.4, max:<1)\n"
    "\n"
    "EXAMPLE\n"
    "  %s -p test.scan -s world-shoreline.shp --point '[-5;48]'\n",prog_name);
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nitems;
  int n,status;
  string keyword,input, PolygonFile, FrameString, PointString, rootname;
  string ShorelinesFile, BoundariesFile, MeshFile, output, zone;
  string OpenedFile, RigidFile, ResolutionString;
  string DescriptorFile, DescriptorRoot, CriteriaFile, BelFile;
  string inFormatName,outFormatName,TmpFile;
  string SampleFile;
  string zonesfilename;
  int inFormat,outFormat;
  vector<plg_t> extraction, limits, shorelines, sampled, smoothed, boundaries, opened, rigid;
  frame_t frame;
  point2D_t point;
  projPJ projection;
  bool debug=false;
  bool extract=false, smooth=true, shorelines2limits=true;
  bool LimitsSpecified=false;
  bool variable_resolution=false, tuned=true;
  bool enforcement=true, show_map=false, show_mesh=false, adjust=false, all=false;
  bool subdivide;
  double factor=-1;
  double lengthscale=0.f;
  double dx=0.f, mindx, sampling=-1.0;
  bool open_sea=true;
  float isocontour=0.0;
  criteria_t criteria;
  
  fct_echo( argc, argv);

  string path=argv[0];
  size_t pos = path.find_last_of("/");
  if( pos != string::npos ){
    path = path.substr(0,pos+1);
    }
  else {
    path="";
    }
  string generator;

  n=1;
  while (n < argc) {
    keyword=argv[n];
    if("-sampling"==keyword){
      printf("\ERROR: -sampling xxx is now deprecated sampling distance setting (given in km)\n");
      printf("please further use --sampling_length xxx, given in meters\n\n");
      print_help(argv[0]);
      exit(-1);
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case '-' :
          if(keyword=="--help"){
            print_help(argv[0]);
            exit(0);
            }
          else if("--extract"==keyword){
            extract=true;
            n++;
            }
          else if("--sampling"==keyword){
            printf("\nWARNING: --sampling xxx is obsolete sampling distance setting (given in km)\n");
            printf("please further use --sampling_length xxx, given in meters\n");
            sscanf(argv[n+1],"%lf",&sampling);
            sampling*=1000.0;
            n++;
            n++;
            }
          else if("--sampling_length"==keyword){
            sscanf(argv[n+1],"%lf",&sampling);
            n++;
            n++;
            }
          else if("--rootname"==keyword){
            rootname= argv[n+1];
            n++;
            n++;
            }
          else if("--point"==keyword){
            PointString= argv[n+1];
            n++;
            n++;
            }
          else if("--isocontour"==keyword){
            nitems=sscanf(argv[n+1],"%f",&isocontour);
            n++;
            n++;
            }
          else if("--factor"==keyword){
            sscanf(argv[n+1],"%f",&factor);
            n++;
            n++;
            }
          else if("--resolution"==keyword){
            ResolutionString= argv[n+1];
            variable_resolution=true;
            n++;
            n++;
            }
          else if("--debug"==keyword){
            debug=true;
            n++;
            }
          else if("--not-tuned"==keyword){
            tuned=false;
            n++;
            }
          else if("--no-smooth"==keyword){
            smooth=false;
            n++;
            }
          else if("--no-enforcement"==keyword){
            enforcement=false;
            n++;
            }
          else if("--show_map"==keyword){
            show_map=true;
            n++;
            }
          else if("--show_mesh"==keyword){
            show_mesh=true;
            n++;
            }
          else if("--limits"==keyword){
            shorelines2limits=false;
            n++;
            }
          else if("--adjust"==keyword){
            adjust=true;
            n++;
            }
           else if("--lake"==keyword){
            open_sea=false;
            n++;
            }
         else if("--all"==keyword){
            all=true;
            n++;
            }
          else {
            printf("unknown option "+keyword+"\n");
            print_help(argv[0]);
            exit(-1);
            }
          break;

        case 'a' :
          OpenedFile=argv[n+1];
          n++;
          n++;
          break;
        
        case 'm' :
          MeshFile=argv[n+1];
          n++;
          n++;
          break;
        
        case 'b' :
          BoundariesFile=argv[n+1];
          n++;
          n++;
          break;
        
        case 'c' :
          CriteriaFile=argv[n+1];
          n++;
          n++;
          break;
        
        case 'd' :
          DescriptorFile=argv[n+1];
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
        
        case 'h' :
          print_help(argv[0]);
          exit(0);
        default:
          printf("unknown option "+keyword+"\n");
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
        if(input=="") {
          input= argv[n];
          n++;
          }
        else {
          printf("unknown option "+keyword+"\n");
          print_help(argv[0]);
          exit(-1);
          }
        break;
      }
    }
  
  vector<float>  isovalues;
  SGfield_t<float> resolution;
  grid_t grid;
  
  if(rootname=="") rootname="anonymous";
  
  if(DescriptorFile!="") {
    bool check=false;
    status=plg_read_descriptor(DescriptorFile.c_str(), sampled, check, debug);
    projection=plg_DefineProjection(sampled, 0);
    goto finished;
    }
  
/*-----------------------------------------------------------------------------
  default coastal sampling when not given in program options */
  if(sampling==-1) {
    sampling=5000.0;
    printf("\nWARNING : coastal sampling not specified, use default %lf km; user specification available through: --sampling <value>\n\n", sampling);
    }
/*-----------------------------------------------------------------------------
  default smoothing ratio when not given in program options */
  if(factor==-1) {
    factor=0.5;
    printf("\nWARNING : (smoothing length / resolution) factor not specified, use default %lf; user specification available through: --factor <value>\n\n", factor);
    }
/*-----------------------------------------------------------------------------
  default smoothing scale when not given in program options */
  if(lengthscale==0.) {
    lengthscale=sampling*factor;
    printf("\nWARNING : shoreline smoothing length not specified, use default %lf m; user specification available through: -l <value>\n\n", lengthscale);
    }
  
/*-----------------------------------------------------------------------------
  default working grid resolution when not given in program options */
  if(dx==0.) dx=lengthscale/4.0;
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  model limits setup. It can be performed by giving a:
  
  - region name
  - frame
  - polygon
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  
  if(zone!="") {
    status=plg_SetZone(zone.c_str(), frame, rootname, point);
    if(status!=NC_NOERR) NC_TRAP_ERROR(wexit,status,1,"zone prescription"+zone+" error");
    limits.push_back(plg_t(frame, PLG_INIT_SEPARATE));
    status=plg_cartesian((projPJ) 0, limits);
    extract=true;
    LimitsSpecified=true;
    }
  
//   if(notebook!="") {
// /*-----------------------------------------------------------------------------
//     Read notebook data */
//     printf("#################################################################\n");
//     printf("load notebook file: %s\n",notebook);
//     status=load_notebook(notebook, &cgrid, &sgrid, &projection);
//     if(status==0) printf("notebook successfuly processed\n",notebook);
//     else  {
//       TRAP_ERR_EXIT(-1,"exiting\n");
//       }
//     plg_t GridLimits=plg_GridLimits(sgrid);
//     status=plg_save("grid.plg", PLG_FORMAT_SCAN, &GridLimits, 1);
//     }
  
  if(FrameString!="") {
    status=plg_DecodeFrame(FrameString, frame);
    if(status!=0) TRAP_ERR_EXIT(status," frame prescription"+FrameString+" error, exit\n");
    limits.push_back(plg_t(frame, PLG_INIT_SEPARATE));
    status=plg_cartesian((projPJ) 0, limits);
    extract=true;
    LimitsSpecified=true;
    }
  
  if(PolygonFile!="") {
    inFormat=plg_find_format(PolygonFile);
    status=plg_load(PolygonFile, inFormat, limits);
    if(status!=0) TRAP_ERR_EXIT(status," plg_load(\""+PolygonFile+"\", %d) error, exit\n",inFormat);
    extract=true;
    LimitsSpecified=true;
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  point gives "wet" indicator
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(PointString!="") {
    status=plg_DecodePosition(PointString, point);
    if(status!=0) TRAP_ERR_EXIT(status," ocean mark prescription"+PointString+" error, exit\n");
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  input is the working shoreline file; if not given, will be extracted from the
  shoreline database
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(input =="") {
    input=rootname+"-shoreline-extraction.plg";
    extract=true;
    }
  else {
    extract=false;
    }

  if(BoundariesFile!="") {
    printf("#################################################################\n");
    printf("reload boundary file: %s\n", BoundariesFile.c_str());
    fflush(stdout);
/*------------------------------------------------------------------------------
    re-loading directly a boundary file means loosing marine/terrestrial flag information*/
    inFormat=plg_find_format(BoundariesFile);
    status=plg_load(BoundariesFile, inFormat, boundaries);
    if(status!=NC_NOERR) NC_TRAP_ERROR(wexit,status,1,"plg_load(\""+input+"\",%d,,) error",inFormat);
    extract=false;
    LimitsSpecified=true;
    smooth=false;
    shorelines2limits=false;
    }
  
  if(!LimitsSpecified) {
    printf("*** model limits not specified ***\n");
    print_help(argv[0]);
    exit(-1);
    }
  if(output=="") output=rootname+"-boundaries.plg";
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  in case where we start from a raw shorelines set, make a first extraction to
  limit computational burden.
  
  it uses the model prescribed limits as a selection polygon
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(extract) {
/*-----------------------------------------------------------------------------
    default shorelines */
    if(ShorelinesFile=="") {
      ShorelinesFile="/home/data/shorelines/gshhs-2.2/gshhs_f.cst";
      printf("\nWARNING : using "+ShorelinesFile+" as default shoreline file (as none was defined with -s option)\n\n");
      }
    status=check_limits(limits);
    
    if(status!=0) {
      printf("limits for extraction is NOT consistent\n");
      goto error;
      }
/*------------------------------------------------------------------------------
    -180°/180° longitude issue*/
    frame_t frame=plg_recale(limits, limits[0].t[0]);
    
/*------------------------------------------------------------------------------
    extract shorelines dataset*/
    printf("#################################################################\n");
    printf("extract shorelines from %s\n", ShorelinesFile.c_str());
    fflush(stdout);
    vector<plg_t> dilated=plg_dilatation_cartesian(limits, 1.05, debug);
    string echo="";
    if(debug) {
      echo=rootname;
      }
    int target;
    if(open_sea) target=0;
    else target=1;
    extraction=plg_extract(ShorelinesFile, dilated, echo, PLG_SPHERICAL, target, debug);
    if(extraction.size()==0) {
      printf("extraction failed, abort\n");
      goto error;
      }
    printf("#################################################################\n");
    printf("store extracted shorelines in %s\n",input.c_str());
    fflush(stdout);
//     printf("to re-use change shorelines input with : -s %s\n",input.c_str());
    status=plg_save(input.c_str(), PLG_FORMAT_SCAN, extraction);
    inFormatName="SCAN";
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  variable resolution sampling can be perfomed by giving a polygon-set (nei
  format) and point value paires in an separate file
   
  typical option argument: --resolution "../resolution.nei ../resolution.dat"
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(variable_resolution) {
    printf("#################################################################\n");
    printf("load variable resolution from %s\n", ResolutionString.c_str());
    fflush(stdout);
    vector<string> arguments=string_split(ResolutionString, " ");
    frame_t frame=plg_spherical_minmax(limits);
//     status=map_set2Dgrid(&grid, frame, 100);
    status=map_set2Dgrid(&grid, frame, 400);
    resolution.x=new float[grid.Hsize()];
    const vector<float> v =Polygons::initialiseValue<double, float>(arguments[0], arguments[1], grid, sampling);
    for (size_t i = 0; i < v.size(); i++) resolution.x[i] = v[i];
    resolution.grid=&grid;
    map_printgrid(grid);
    float slopemax=0.2;
    range_t<float> range;
    status=resolution_smooth(grid, resolution.x, resolution.mask, range, limits, slopemax);
    float mask=-1;
    status=save_SGXY("resolution.nc", grid, resolution.x, mask, "resolution", "m", "resolution", 1);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  smoothing is used to filter sub-scale details (with respect to desired coastal
  sampling length) and aggregate packed micro-structures having a macro impact
  
  it acts on the shorelines polygons
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(smooth) {
    printf("#################################################################\n");
    printf("filter shorelines: working grid resolution=%lf m\n", dx);
    fflush(stdout);
/*------------------------------------------------------------------------------
    load new shorelines data set*/
    inFormat=format(input, inFormatName);
    status=plg_load(input, inFormat, shorelines);
    if(status!=0) goto error;
    status=plg_CheckDuplicated(shorelines, true, rootname, false);
    
/*------------------------------------------------------------------------------
    -180°/180° longitude issue*/
    frame_t frame=plg_recale(shorelines, limits[0].t[0]);
/*------------------------------------------------------------------------------
    define isovalues for smoothed contours*/
    isovalues.clear();
    for(int k=-2;k<5;k++) isovalues.push_back(0.-k*0.1);
    
/*------------------------------------------------------------------------------
    initialize filter parameters*/
    plg_filter_t plg_filter;
    subdivide=false;
    plg_filter.init(lengthscale, resolution, factor, dx, subdivide, tuned, enforcement, show_map, show_mesh);
    
    smoothed=plg_smooth(rootname, shorelines, plg_filter, isovalues, isocontour, point, 0, debug);
    
    string tmp=rootname+"-shoreline-smoothed" +".plg";
    status=plg_save(tmp.c_str(), PLG_FORMAT_SCAN, smoothed);
    }
  else {
/*------------------------------------------------------------------------------
    load new shorelines data set*/
    if(shorelines2limits) {
      inFormat=format(input, inFormatName);
      status=plg_load(input, inFormat, smoothed);
      if(status!=0) goto error;
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  create model boundary polygons from (potentially filtered) shorelines and
  selection polygon:
  
  -keep ocean sections unchanged (and flagged as 'M')
  
  -set shorelines along coastal sections (and flagged as 'T')
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(shorelines2limits) {
    printf("#################################################################\n");
    printf("build model limits\n");
    fflush(stdout);
/*------------------------------------------------------------------------------
    restore original extraction frame*/
    plg_destroy_entries(extraction);
    string echo="";
    if(debug) {
      echo=rootname+"-debug-limits";
      }
/*------------------------------------------------------------------------------
    -180°/180° longitude issue */
    frame_t frame=plg_recale(limits, limits[0].t[0]);
    
    degree_recale(&(point.x), frame.x_center());
    
    boundaries=plg_extract(smoothed, limits, echo, PLG_SPHERICAL, 1, debug);
    outFormat=format(output, outFormatName);
    status=plg_save(output.c_str(), outFormat, boundaries);
    
    status=plg_RemoveInlandLimits(boundaries, point.x, point.y, (projPJ) 0, PLG_POINT_EXTERIOR, debug);
    if(status==-1) {
      printf("ocean point (%lf;%lf) not included in model polygons\n", point.x, point.y);
      string FailName=rootname+"-failed.plg";
      status=plg_save(FailName.c_str(), PLG_FORMAT_SCAN, boundaries);
      goto error;
      }
    printf("#################################################################\n");
    printf("store model limits %s\n",output.c_str());
    fflush(stdout);
    outFormat=format(output, outFormatName);
    status=plg_save(output.c_str(), outFormat, boundaries);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  sample rigid sections (T flag) of boundary polygons at required resolution
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("sample shorelines (rigid limits) : default size=%lf meters\n", sampling);
  fflush(stdout);
  projection=plg_DefineProjection(boundaries, 0);
  status=plg_cartesian(projection, boundaries);
  for(int s=0;s<boundaries.size();s++) {
    if(boundaries[s].flag==0) {
      boundaries[s].SetFlag('T');
//       boundaries[s].flag=new char[boundaries[s].npt-1];
//       for(int kk=0;kk<boundaries[s].npt-1;kk++) boundaries[s].flag[kk]='T';
      }
    int equalize=1, verbose=0;
    if(debug) verbose=1;
    plg_t test;
    char mask='M';
    if(all) mask='?';
    if(!variable_resolution) {
      test=plg_resample_rigid_standard(boundaries[s], projection, sampling, mask, equalize, verbose);
      }
    else {
      test=plg_resample_rigid_variable(boundaries[s], projection, resolution, sampling, mask, equalize, verbose);
      }
/*------------------------------------------------------------------------------
    -180°/180° longitude issue */
    frame_t frame=plg_recale(&test, 1);
    if(test.npt>2) {
      status=plg_checkAutoSecant((plg_t *) &test, 1, PLG_SPHERICAL);
      if(status!=0) {
        printf("sample shorelines failed for segment size=%lf km\n", s, sampling);
        }
      else sampled.push_back(test);
      }
    }
 
  SampleFile=rootname+"-boundaries-resampled" +".plg";
  status=plg_save(SampleFile.c_str(), PLG_FORMAT_SCAN, sampled);
  
/*------------------------------------------------------------------------------
  create boundary descriptor file needed by mesh generator*/
  DescriptorRoot=rootname+"-boundaries";
  status=plg_write_descriptor(DescriptorRoot, sampled, debug);
    
finished:

  if(OpenedFile!="") {
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    assembly newly created model boundary polygons with prescribed open limits
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
    
    vector<plg_t> imported;
    vector<plg_t> selected;

    printf("#################################################################\n");
    printf("aggregate processed boundaries with prescribed open limits (%s)\n", OpenedFile.c_str());
    fflush(stdout);
    RigidFile=rootname+"-boundaries-rigid.scan";

    status=plg_load(OpenedFile, PLG_FORMAT_SCAN, imported);
    if(status!=0) goto error;

    status=plg_cartesian(projection, imported);
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    prepare resolution adjustment
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    vector<plg_t> circles, centres;
    vector<double> sizes, dum;
    status=plg_circles(imported, projection, 0.1, centres, dum);
    status=plg_circles(imported, projection, 1.5, circles, sizes);
    
/*------------------------------------------------------------------------------
    insure mono-block connecting polygons */
    bool repair=true;
    status=plg_CheckClosure(imported, repair, 0, debug);
    
/*------------------------------------------------------------------------------
    imported open limits special flag */
    for(int s=0;s<imported.size();s++) imported[s].SetFlag('X');
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    check consistency (i.e. ability to connect together) between actual open
    boundaries and imported ones
    
    - originally aimed for mesh upgrade purposes, i.e. imported open boundaries
      will be of similar position and number as the newly created ones
      
    - now extended to more general use
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    status=plg_CheckOpenConsistency(rootname, sampled, imported, selected, projection, debug);
    
//     if(status==0) {
// /*------------------------------------------------------------------------------
//       no problem detetected, nominal procedure apply */    
//       status=plg_load(RigidFile, PLG_FORMAT_SCAN, selected);
//       if(status!=0) goto error;
//       status=plg_cartesian(projection, selected);
//       }
//     else {
    
    if(status!=0) {
/*------------------------------------------------------------------------------
      problems detetected, patch applies to re-create a more convenient rigid set */    
      for(int s=0;s<selected.size();s++) {
        if(selected[s].flag==0) continue;
        if(selected[s].flag[0]!='S') continue;
/*------------------------------------------------------------------------------
        WARNING : S flag will turn into no resampling... evolution bug ? to be fixed? */
        int equalize=1, verbose=0;
        if(debug) verbose=1;
        plg_t test;
        char mask='M';
        if(!variable_resolution) {
          test=plg_resample_rigid_standard(selected[s], projection, sampling, mask, equalize, verbose);
          }
        else {
          test=plg_resample_rigid_variable(selected[s], projection, resolution, sampling, mask, equalize, verbose);
          }
        if(test.npt>2) {
          test.SetFlag('T');
          selected.push_back(test);
          selected[s].destroy();
          }
        else {
          selected[s].SetFlag('T');
          }
        }
      status=plg_cartesian(projection, selected);
      RigidFile=rootname+"-boundaries-selected.scan";
      status=plg_save(RigidFile.c_str(), PLG_FORMAT_SCAN, selected);
      }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    assembly newly constructed and imported limits
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    DescriptorRoot=rootname+"-assembled";
//     debug=true;
    sampled=plg_assembly_limits(selected, imported, projection, DescriptorRoot, debug);

/*------------------------------------------------------------------------------
    assembly might create undesirable external polygons, discard them         */
    status=plg_RemoveInlandLimits(sampled, point.x, point.y, projection, PLG_POINT_EXTERIOR, debug);
    
    SampleFile=rootname+"-assembled-boundaries" +".plg";
    status=plg_save(SampleFile.c_str(), PLG_FORMAT_SCAN, sampled);

/*------------------------------------------------------------------------------
    creates the model boundaries descriptor need by mesh-generator */
    status=plg_write_descriptor(DescriptorRoot, sampled, debug);
    
    if(adjust) {
      printf("#################################################################\n");
      printf("adjust coastal resolution to fit imported open limits\n");
      fflush(stdout);
      vector<plg_t> splitted=plg_split(sampled, debug);
      vector<plg_t> tmp;
      
      status=plg_merge(tmp, splitted);
      status=plg_merge(tmp, circles);
      status=plg_merge(tmp, centres);
      
      TmpFile=rootname+"-adjust-base.plg";
      status=plg_save(TmpFile.c_str(), PLG_FORMAT_SCAN, tmp);
      
      status=adjust_assembly(splitted, circles, centres, sizes, projection, debug);

      SampleFile=rootname+"-assembled-splitted-adjusted" +".plg";
      status=plg_save(SampleFile.c_str(), PLG_FORMAT_SCAN, splitted);
    
      status=plg_CheckClosure(splitted,true, 1, debug);
    
      status=plg_write_descriptor(DescriptorRoot, splitted, debug);
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  if a criteria template file given, produce a tuned one
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(CriteriaFile!="") {
    printf("#################################################################\n");
    printf("load resolution criteria file : %s\n",CriteriaFile.c_str());
    fflush(stdout);
    status=fe_load_criteria(CriteriaFile.c_str(), &criteria);
    if(status !=0) {
      goto error;
      }
    }
  else {
    printf("#################################################################\n");
    printf("initialise resolution criteria from COASTAL default\n");
    fflush(stdout);
    status=fe_InitCriteria("COASTAL",criteria);
    }
  
/*------------------------------------------------------------------------------
  adjust criteria parameters accordingly to mesh-limits tuning */
//   criteria.shelf_minsize=sampling*3.0/5.0/1000.0;
  criteria.shelf_minsize=sampling/1000.0;
  criteria.shelf_maxsize=criteria.shelf_minsize*5;
  criteria.cellsize=criteria.shelf_minsize/2.;
  
  mindx=criteria.cellsize;
  
  criteria.resample_openlimits=1;
    
  DescriptorFile=DescriptorRoot+".desc";
  strcpy(criteria.descriptor,DescriptorFile.c_str());
  
  if(variable_resolution) {
    vector<string> arguments=string_split(ResolutionString, " ");
    PolygonSetBase<double,float> p=Polygons::load<double,float>(arguments[0],arguments[1]);
    FILE *zones, *points;
    zonesfilename=rootname+"-zones.dat";
    int count=0;
    criteria_t local;
    zones=fopen(zonesfilename.c_str(),"w");
    points=fopen(arguments[1].c_str(),"r");
    fprintf(zones, "%s\n", arguments[0].c_str());
    for (PolygonSetBase<double,float>::iterator i = p.begin(); i != p.end();i++) {
      double x,y;
      char value[1024];
      PolygonBase<double, float> q = *i;
      local=criteria;
      local.shelf_minsize=q.z0()/1000.0;
      local.shelf_maxsize=local.shelf_minsize*10;
      local.cellsize=local.shelf_minsize/2.;
      mindx=min(mindx,local.cellsize);
      fscanf(points,"%lf %lf %s",&x,&y,value);
      sprintf(value,"%2.2d",count);
      CriteriaFile=rootname+"-"+value+".crt";
      status=fe_save_criteria(CriteriaFile.c_str(), local);
      fprintf(zones,"%9.4lf %9.4lf %s\n", x,y,CriteriaFile.c_str());
      count++;
      }
    fclose(zones);
    fclose(points);
    }
  
  criteria.cellsize=mindx;
 
  CriteriaFile=rootname+".crt";
  printf("#################################################################\n");
  printf("save resolution criteria: %s\n", CriteriaFile.c_str());
  fflush(stdout);
  status=fe_save_criteria(CriteriaFile.c_str(), criteria);
  if(status !=0) {
    goto error;
    }
  
/*------------------------------------------------------------------------------
  Create a BEL file for T-UGOm */
  BelFile=rootname+".bel";
  printf("#################################################################\n");
  printf("save BEL file: %s\n", BelFile.c_str());
  fflush(stdout);
  status=fe_write_boundarycode(BelFile.c_str(), sampled);
  sampled.clear();
  
  generator=path+"mesh-generator";
  printf("\nYou might now check the criteria file "+CriteriaFile+", then run the mesh generation command:\n\n");
  printf(generator+" -d "+DescriptorRoot+".desc "+CriteriaFile);
  if(variable_resolution) printf(" -z "+zonesfilename);
  printf(" --rootname "+rootname+" -b <bathymetry file>\n\n");
 
  STDOUT_BASE_LINE(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(0);
  
error:
  STDOUT_BASE_LINE("error detected, quit ... \n");
  exit(-1);
}
