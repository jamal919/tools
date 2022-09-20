
/*******************************************************************************

  T-UGO tools, 2006-2013

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief area and polygone selection functions
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
#include "geo.h"
#include "functions.h"

#include "polygons.hpp"
#include "exceptions.hpp"

using namespace Polygons; // for ElementaryCycles and Adjacency

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_DecodePosition(const string & input, point2D_t & point)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//  [lon;lat]
  int nitems;
  
  nitems=sscanf(input.c_str(),"[%lf;%lf]", &point.x, &point.y);
  if(nitems!=2) return (-1);
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_SetZone(const char *zone, frame_t & frame, string & rootname, point2D_t & point)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
/*------------------------------------------------------------------------------
  extraction frame*/

  status=-1;

  if(strcmp(zone,"korea")==0) {
    frame=frame_t(125.0, 128.0, 33.5, 36);   // Korea
    point.x=125.5;
    point.y= 35.0;
    rootname="korea";
    status=0;
    goto finish;
    }

  if(strcmp(zone,"taiwan")==0) {
    frame=frame_t(119, 120.5, 23, 24);       // Taïwan
    point.x=120.0;
    point.y= 23.5;
    rootname="taiwan";
    status=0;
    goto finish;
    }

  if(strcmp(zone,"north-china")==0) {
    frame=frame_t (113.0, 123.0, 26.0, 31.0);    // North Continental China
    point.x=122.5;
    point.y= 27.5;
    rootname="north-china";
    status=0;
    goto finish;
    }

  if(strcmp(zone,"central-china")==0) {
//    frame=frame_t(115.0, 120.5, 22.5, 26);   // Continental China
    frame=frame_t(112.0, 123.0, 21.0, 31.0);   // Continental China
    point.x=117.5;
    point.y= 23.0;
    rootname="central-china";
    status=0;
    goto finish;
    }

  if(strcmp(zone,"south-china")==0) {
    frame=frame_t(105.0, 115.0, 18, 24);     // South Continental China
    point.x=112.5;
    point.y= 21.0;
    rootname="south-china";
    status=0;
    goto finish;
    }

  if(strcmp(zone,"lofoten")==0) {
    frame=frame_t(4, 17, 60, 70);             // Lofoten
    point.x= 10.0;
    point.y= 65.0;
    rootname="lofoten";
    status=0;
    goto finish;
    }
  
//   frame_t frame(17, 27, 68, 72);             // Barents

finish:
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_DecodeFrame(const string input, frame_t & frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//  [tmin:tmax:dt; pmin:pmax:dp]
  int status;
  string s=input;
  string substring,longitude, latitude;
  size_t pointer;
  int count,nitems;
  
  pointer = s.find('[');
  if(pointer==string::npos) return (-1);
  pointer = s.find(';');
  if(pointer==string::npos) return (-1);
  pointer = s.find(']');
  if(pointer==string::npos) return (-1);
  
  pointer = s.find('[');
  s[pointer]=' ';
  pointer = s.find(']');
  s[pointer]=' ';
  
  pointer = s.find(';');
  s[pointer]=' ';
  longitude =s.substr(0,pointer);
  latitude  =s.substr(pointer+1);
  
  s=longitude;
  pointer = s.find(':');
  if(pointer==string::npos) return (-1);
  
  count=0;
  while(pointer!=string::npos) {
    s[pointer]=' ';
    pointer = s.find(':');
    count++;
    }
  switch(count) {
    case 1:
      nitems=sscanf(s.c_str(),"%lf %lf", &frame.xmin, &frame.xmax);
      if(nitems!=2) return (-1);
      break;
    case 2:
      nitems=sscanf(s.c_str(),"%lf %lf", &frame.xmin, &frame.xmax);
      if(nitems!=2) return (-1);
      break;
    default:
     return (-1);
    }
  
  s=latitude;
  pointer = s.find(':');
  if(pointer==string::npos) return (-1);
  
  count=0;
  while(pointer!=string::npos) {
    s[pointer]=' ';
    pointer = s.find(':');
    count++;
    }
  switch(count) {
    case 1:
      nitems=sscanf(s.c_str(),"%lf %lf", &frame.ymin, &frame.ymax);
      if(nitems!=2) return (-1);
      break;
    case 2:
      nitems=sscanf(s.c_str(),"%lf %lf", &frame.ymin, &frame.ymax);
      if(nitems!=2) return (-1);
      break;
    default:
     return (-1);
    }
  
  status=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_t  plg_convex(plg_t & p, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  plg_t q;
  
  bool closed=p.closed();
  
  q.duplicate(p);
  
  for(int s=0;s<q.npt;s++) {
    double angle=q.angle(s);
    if(angle<0) {
      status=plg_delete_point(&q,s);
      if((s==0) && closed) plg_copy_point(&q,q.npt-1,0);
      s--;
      if((s==q.npt-1) && closed) plg_copy_point(&q,0,q.npt-1);
      }
    if (debug) printf("point %d, angle %lf\n",s, angle);
    }
  
  for(int s=q.npt-1;s>=0;s--) {
    double angle=q.angle(s);
    if(angle<0) {
      status=plg_delete_point(&q,s);
      if((s==q.npt-1) && closed) plg_copy_point(&q,0,q.npt-1);
      s++;
      if((s==0) && closed) plg_copy_point(&q,q.npt-1,0);
      }
    if (debug) printf("point %d, angle %lf\n",s, angle);
    }
  
  return(q);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_dilatation(plg_t & p, double factor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  apply polygon dilatation with respect to center

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double x,y;
  
  if(p.npt==0) return(-1);
  
  x=y=0.0;
  for(int k=0;k<p.npt;k++) {
    x+=p.x[k];
    y+=p.y[k];
    }
  x/=(double) p.npt;
  y/=(double) p.npt;
  
  for(int k=0;k<p.npt;k++) {
    double dx=factor*(p.x[k]-x);
    double dy=factor*(p.y[k]-y);
    p.x[k]=x+dx;
    p.y[k]=y+dy;
    p.t[k]=x+dx;
    p.p[k]=y+dy;
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t>  plg_dilatation_cartesian(vector<plg_t> & p, double factor, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  dilatation of a given polygon
  
  use cartesian world, sub-samples original polygon
  
  calls plg_setdirect(p)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  vector<plg_t> q;
  
  projPJ projection=plg_DefineProjection(p, 0);
  status=plg_cartesian(projection, p);
  status=plg_setdirect(p);
  
  for(int s=0;s<p.size();s++) {
    plg_t tmp, r;
    r=plg_convex(p[s],debug);
    status=plg_dilatation(r, factor);
    q.push_back(r);
    }
  
  if(debug) {
    status=plg_spherical(projection, q);
    status=plg_save("dilatation-raw.plg",PLG_FORMAT_UNKNOWN, q);
    status=plg_cartesian(projection, q);
    }
  
  status=plg_checkSecant(p, q, debug);
  
  if(status!=0) {
    if(not debug){/* if not saved already */
      /* save */
      status=plg_spherical(projection, q);
      status=plg_save("dilatation-raw.plg",PLG_FORMAT_UNKNOWN, q);
      status=plg_cartesian(projection, q);
      }
    printf("%s: dilated polygon intersect original one, failure\n",__func__);
    q.clear();
    }
  
  status=plg_spherical(projection, q);
  
  if(debug) status=plg_save("dilatation-1.plg", PLG_FORMAT_SCAN, q);
  
  double radius=50000.0;
  int equalize=1;

  vector<plg_t> qq=plg_subdivide(q, radius, equalize, PLG_CARTESIAN, debug);
  status=plg_spherical(projection, qq);
  
  pj_free(projection);
  
  plg_destroy_entries(q);
  q.clear();
  
  q=qq;
  
  frame_t frame=plg_recale(q, p[0].t[0]);
  
  if(debug) status=plg_save("dilatation-2.plg", PLG_FORMAT_SCAN, q);
  
  status=plg_cartesian((projPJ) 0, q);
  status=plg_cartesian((projPJ) 0, p);
  
  return(q);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_FrameSelection(plg_array_t *array, frame_t frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------------
  Delete points, and polygones if no points left
------------------------------------------------------------------------------*/
  int i, j, intersection, status;
  double x,y;
  bool ok[2];
  int *list,cut;

/**----------------------------------------------------------------------------
  first detect polygones totally out of selection*/
  for(i=0;i<array->n;i++) {
    intersection=0;
    for(j=0;j<(array->p)[i].npt;j++) {
      x=(array->p)[i].x[j];
      y=(array->p)[i].y[j];
      ok[0]=frame.inside(x, y);
      if(ok[0]) {
        intersection++;
        }
      }
/**----------------------------------------------------------------------------
    why 2 ? */
//     if(intersection<2) (array->p)[i].npt=0;
    if(intersection<1) (array->p)[i].npt=0;
    }

  status= plg_compact_entries (&(array->p),&(array->n));

/**----------------------------------------------------------------------------
  then detect and cut polygones not totally out of selection*/
redo:
  cut=0;
  for(i=0;i<array->n;i++) {
//    if((array->p)[i].npt-1;j++) {
    for(j=0;j<(array->p)[i].npt-1;j++) {
      x=(array->p)[i].x[j];
      y=(array->p)[i].y[j];
      ok[0]=frame.inside(x, y);
      x=(array->p)[i].x[j+1];
      y=(array->p)[i].y[j+1];
      ok[1]=frame.inside(x, y);
      if(ok[0] && !ok[1] && (j!=(array->p)[i].npt-2)) {
        list=plg_cut_polygone(&(array->p),&(array->n), i, j+1);
        plg_swap_entries(array->p,list[0],array->n-2);
        plg_swap_entries(array->p,list[1],array->n-1);
//        updatemin(&i,array->n-3);
        delete[] list;
        cut++;
        break;
        }
      if( !ok[0] && ok[1] && (j!=0)) {
        list=plg_cut_polygone(&(array->p),&(array->n), i, j);
        plg_swap_entries(array->p,list[0],array->n-2);
        plg_swap_entries(array->p,list[1],array->n-1);
//        updatemin(&i,array->n-3);
        cut++;
        delete[] list;
        break;
        }
      }
    }

/**----------------------------------------------------------------------------
  re-detect polygones totally out of selection*/
  for(i=0;i<array->n;i++) {
    intersection=0;
    for(j=0;j<(array->p)[i].npt;j++) {
      x=(array->p)[i].x[j];
      y=(array->p)[i].y[j];
      ok[0]=frame.inside(x, y);
      if(ok[0]) {
        intersection++;
        }
      }
//    if(intersection<2) (array->p)[i].npt=0;
    if(intersection<1) (array->p)[i].npt=0;
    }

  status= plg_compact_entries (&(array->p),&(array->n));
  if(cut!=0) goto redo;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_FrameSelection(plg_t **polygones, int *npolygones, frame_t frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  plg_array_t array;
  
  array.p=*polygones;
  array.n=*npolygones;
  
  status=plg_FrameSelection(&array, frame);
  
  *polygones=array.p;
  *npolygones=array.n;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_FrameSelection(vector<plg_t> & polygons, frame_t frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  plg_array_t array=plg_vector2array(polygons);
    
  status=plg_FrameSelection(&array, frame);
  
  polygons=plg_array2vector(array);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_reduce_selection(vector<plg_t> & polygons, const vector<plg_t> & selection, int mode, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// keep from a set of polygones only the points within a selected area
/**
\param mode PLG_SPHERICAL or PLG_CARTESIAN
\return -1 if the selection is not closed, 0 otherwise
*/
/*----------------------------------------------------------------------------*/
{
  int i, j, status;
  frame_t frame;
  int inside,flag,check=0;
  int cut;
  int list;
  bool ok[2];
  double x,y;
  
  if(plg_isclosed(selection[0])==0) return(-1);
  
  frame= plg_spherical_minmax(selection);
  status=plg_QuickFrameReduction(polygons, frame, mode);
  
redo:

  cut=0;
  for(i=0;i<polygons.size();i++) {
    for(j=0;j<polygons[i].npt-1;j++) {
      x=polygons[i].x[j];
      y=polygons[i].y[j];
      inside=0;
      flag=plg_single(selection[0], x, y, &inside, mode, check);
//      ok[0]=(flag==PLG_POINT_INTERIOR);
      ok[0]=(inside==1);
      x=polygons[i].x[j+1];
      y=polygons[i].y[j+1];
      inside=0;
      flag=plg_single(selection[0], x, y, &inside, mode, check);
//      ok[1]=(flag==PLG_POINT_INTERIOR);
      ok[1]=(inside==1);
      if(ok[0] && !ok[1] && (j!=polygons[i].npt-2)) {
        list=plg_cut_polygone(polygons, i, j+1);
        plg_swap_entries(polygons,i,polygons.size()-2);
        plg_swap_entries(polygons,list,polygons.size()-1);
        cut++;
        break;
        }
      if( !ok[0] && ok[1] && (j!=0)) {
        list=plg_cut_polygone(polygons, i, j);
        plg_swap_entries(polygons,i,polygons.size()-2);
        plg_swap_entries(polygons,list,polygons.size()-1);
        cut++;
        break;
        }
      }
    }

  bool strict=true;
  status=plg_PolygonReduction(polygons, selection, strict, mode, 0, debug);
  if(cut!=0) goto redo;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_QuickFrameReduction(plg_array_t *array, frame_t frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------------
  Delete polygones if no points out of selection
------------------------------------------------------------------------------*/
  int i, j, intersection, status;
  double x,y;
  bool ok[2];

/**----------------------------------------------------------------------------
  first detect polygones totally out of selection*/
  for(i=0;i<array->n;i++) {
    intersection=0;
    for(j=0;j<(array->p)[i].npt;j++) {
      x=(array->p)[i].x[j];
      y=(array->p)[i].y[j];
      ok[0]=frame.inside(x, y);
      if(ok[0]) {
        intersection++;
        }
      }
    if(intersection==0) (array->p)[i].npt=0;
    }

  status=plg_compact_entries (&(array->p),&(array->n));

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_QuickFrameReduction(plg_t **polygones, int *npolygones, frame_t frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  plg_array_t array;
  
  array.p=*polygones;
  array.n=*npolygones;
  
  status=plg_QuickFrameReduction(&array, frame);
  
  *polygones=array.p;
  *npolygones=array.n;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_QuickFrameReduction(vector<plg_t> & polygons, frame_t frame, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int intersection, status;
  int i,j,pivot;
  double x,y;
  bool  ok;
  const size_t n=polygons.size();
  bool *keep=new bool[n];
  double center=frame.x_center();
  
  for(i=0;i<n;i++) {
    intersection=0;
    for(j=0;j<polygons[i].npt;j++) {
      x=polygons[i].x[j];
      y=polygons[i].y[j];
      if(mode==PLG_SPHERICAL)
        x=degree_recale(x, center);
      ok=frame.inside(x, y);
      if(ok) {
        intersection++;
        pivot=j;
        break;
        }
      }
    if(intersection==0) {
      polygons[i].destroy();
      keep[i]=false;
      }
    else {
      if(mode==PLG_SPHERICAL)
        polygons[i].x[pivot]=degree_recale(x, center);
      polygons[i].t[pivot]=polygons[i].x[pivot];
      for(j=pivot+1;j<polygons[i].npt;j++) {
        if(mode==PLG_SPHERICAL) {
          polygons[i].x[j]=degree_recale(polygons[i].x[j], polygons[i].x[j-1]);
          polygons[i].t[j]=polygons[i].x[j];
          }
        }
      for(j=pivot-1;j>=0;j--) {
        if(mode==PLG_SPHERICAL) {
          polygons[i].x[j]=degree_recale(polygons[i].x[j], polygons[i].x[j+1]);
          polygons[i].t[j]=polygons[i].x[j];
          }
        }
      keep[i]=true;
      }
    }
  
  vector<plg_t> p;
  for(i=0;i<n;i++) {
    if(keep[i]) {
      plg_t tmp;
      tmp.duplicate(polygons[i]);
      tmp.id=i;
      p.push_back(tmp);
      polygons[i].destroy();
      }
    }
  
  polygons.clear();
  polygons=p;
  
  delete[] keep;
  
  status=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_PolygonReduction(vector<plg_t> & polygons, const vector<plg_t> & selection, bool strict, int coordinates, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  filter polygons, keep interior set (i.e. selection-wise interior)
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int i, j, flag, intersection, status;
  int check=0;
  int inside;
  frame_t frame;
  double x,y;
  int size;
  
  size=polygons.size();

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  do first quick selection based on min/max frame
    
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  switch(coordinates) {
    case PLG_SPHERICAL:
      frame= plg_spherical_minmax(selection);
      break;
    case PLG_CARTESIAN:
      frame= plg_cartesian_minmax(selection);
      break;
    }
  status=plg_QuickFrameReduction(polygons, frame, coordinates);
  
  if(verbose==1) {
    printf("%s, 1st step : initial size=%d, after quick reduction=%d \n", __func__, size, polygons.size());
    }
  
  if(debug) {
    status=plg_save("PolygonReduction-1.plg", PLG_FORMAT_SCAN, polygons);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  then remove shorelines with no point inside selection polygons
    
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(i=0;i<polygons.size();i++) {
    intersection=0;
    for(j=0;j<polygons[i].npt;j++) {
      x=polygons[i].x[j];
      y=polygons[i].y[j];
      inside=0;
      flag=plg_single(selection[0], x, y, &inside, coordinates, check);
      switch (strict) {
        case true:
          if(flag==PLG_POINT_INTERIOR) {
            intersection++;
            }
          break;
        case false:
          if(flag!=PLG_POINT_EXTERIOR) {
            intersection++;
            }
          break;
        }
      }
    if(intersection==0) {
      polygons[i].destroy();
      polygons.erase(polygons.begin()+i);
      i--;
      }
    }
  
  if(verbose==1) {
    printf("%s, 2nd step : initial size=%d, after full reduction=%d \n", __func__, size, polygons.size());
    }
  
  if(debug) {
    status=plg_save("PolygonReduction-2.plg", PLG_FORMAT_SCAN, polygons);
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_CutAtIntersections(vector<plg_t> & p, vector<plg_t> & q, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  create a new point at intersection of any polygons in p and any another polygon
  in q, if any, and cut both polygons
  
  warning : critical function plg_secantpoint works on spherical coordinates
  
  p and q modified
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k,l,status;
  plg_point_t point;
  int vertex;
  double d;
  int list;
  bool p_truncated=false;
  bool q_truncated=false;
  bool debug=false;
  
  
  for(k=0;k<p.size();k++){
    if(debug) printf("line %d, size=%d\n",k,p[k].npt);
    for(l=0;l<q.size();l++){
      for(size_t kk=0;kk<p[k].npt-1;kk++) {
        line_t a=line_t(p[k],kk);
        for(size_t ll=0;ll<q[l].npt-1;ll++) {
          line_t b=line_t(q[l],ll);
          int flag=plg_secantpoint(a,b,point);
          switch(flag) {
            case PLG_LINES_SECANT:
/*------------------------------------------------------------------------------
              strictly secant lines, add new point and cut */
              status= plg_insert_point(p[k], kk+1, point);
              list= plg_cut_polygone(p, k, kk+1);
              status= plg_insert_point(q[l], ll+1, point);
              list= plg_cut_polygone(q, l, ll+1);
/*------------------------------------------------------------------------------
              p[k] now finished by truncation */
              goto finish_p;
              break;
            case PLG_LINES_SECANT_AT_EXTRIMITY:
              p_truncated=false;
              q_truncated=false;
              vertex= plg_find_point(p[k], PLG_SPHERICAL, point, &d);
              if(d!=0) {
                status= plg_insert_point(p[k], kk+1, point);
                vertex=kk+1;
                p_truncated=true;
                }
              list= plg_cut_polygone(p, k, vertex);
              vertex= plg_find_point(q[l], PLG_SPHERICAL, point, &d);
              if(d!=0) {
                status= plg_insert_point(q[l], ll+1, point);
                vertex=ll+1;
                q_truncated=true;
                }
              list= plg_cut_polygone(q, l, vertex);
              if(p_truncated) goto finish_p;
              if(q_truncated) goto finish_q;
              break;
            }
          }
finish_q:
        if(debug) printf("frame, size=%d\n",q.size());
        continue;
        }
finish_p:
    if(debug) printf("line %d, size=%d\n",k,p[k].npt);
    continue;
      }
// finish_p:
//     if(debug) printf("line %d, size=%d\n",k,p[k].npt);
//     continue;
    }
  
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_VectorSizeV(const vector<plg_t> p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,size;
  
  size=0;
  for(k=0;k<p.size();k++){
    size+=p[k].npt;
    }
  return(size);
}
  

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_VectorSizeL(const vector<plg_t> p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,size;
  
  size=0;
  for(k=0;k<p.size();k++){
    size+=p[k].npt-1;
    }
  return(size);
}
  

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_ImportVertex(vertex_t & vertex, plg_t p, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vertex.lon=p.t[target];
  vertex.lat=p.p[target];
  vertex.h=0;
  vertex.nngh=0;
  vertex.code=1;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_ExportVertex(vertex_t & vertex, plg_t p, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  p.t[target]=vertex.lon;
  p.p[target]=vertex.lat;
  p.x[target]=vertex.lon;
  p.y[target]=vertex.lat;
  
  return(0);
}
  

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_FindVertex(mesh_t & mesh, double t, double p, double & distance)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int found=-1,n;
  double dmin=1.e+10;
  double *d=new double[mesh.nvtxs];
  
  int nprocs=initialize_OPENMP(-1, 0);
#pragma omp parallel for if(nprocs>1)
  for(n=0;n<mesh.nvtxs;n++) {
//    d[n]=geo_distance(mesh.vertices[n].lon,mesh.vertices[n].lat,t,p);
    d[n]=geo_haversin_km(mesh.vertices[n].lon,mesh.vertices[n].lat,t,p);
    }
  
  for(n=0;n<mesh.nvtxs;n++) {
    if(d[n]<dmin) {
      found=n;
      dmin=d[n];
      }
    }
  delete[] d;
  distance=dmin;
  return(found);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_SeekCycles(mesh_t *mesh, Polygons::ElementaryCycles &elementaryCycles, int maxsize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  Make an element list from nodes' neighbours description.
  If consecutive neighbours of a node are connected, they
  form an element with the node
    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i,j,k,l,m,n,rotation=-1;
  int max_nb,nndes,*flag,*type,*size;
  int ntriangles,nquadrangles,npentacles,nelements;
  vector<int>  *tmp;
  int status;
  int nnpe[100];
  int nodes[maxsize-1];
  int cardinal[maxsize-1];

  if(mesh->type==-1) mesh->type=0;

  nndes=mesh->nvtxs;
  max_nb=mesh->nnghm;

  for(i=0; i<nndes; i++) {
    status=fe_order(*mesh,i,-1., 0.,rotation);
    }

  flag=new int[nndes];
  for(i=0; i<nndes; i++) flag[i]=0;

  for(i=0; i<nndes; i++) {
    if(mesh->vertices[i].nngh < 2 ) {
      printf("fe_SeekQuadrangles: pathetic node %d\n",i);
      flag[i]=1;
      }
    }

  tmp= new vector<int> [2*nndes];
  type=new int[2*nndes];
  size=new int[2*nndes];
  
/**-----------------------------------------------------------------------------
  this makes the algorithm unconvenient for polygone detection */
//  for(i=0; i<2*nndes; i++) tmp[i]=new int[maxsize-1];

  ntriangles   = 0;
  nquadrangles = 0;
  npentacles=0;
  nelements=0;

  for(l=0;l<maxsize;l++) cardinal[l]=0;

  for(n=0;n<nndes;n++) {
    if(flag[n]!=0) continue;
    status=fe_order(*mesh,n,-1., 0.,rotation);
    nodes[0]=n;
    for(k=0; k<mesh->vertices[n].nngh; k++) {
      nodes[1]=mesh->vertices[nodes[0]].ngh[k];
      if(nodes[1] <= nodes[0]) continue;
      status = fe_order2(*mesh, nodes[1], nodes[0], rotation);
      nodes[2]=mesh->vertices[nodes[1]].ngh[0];
      if(nodes[2] == nodes[1]) continue;
      if(nodes[2] <= nodes[0]) continue;
//      status = fe_order2(*mesh, nodes[2], nodes[1], rotation);
      for(l=3;l<maxsize;l++) {
        int previous=nodes[l-1];
        if(mesh->vertices[previous].nngh==2) {
          int first=mesh->vertices[previous].ngh[0];
          nodes[l]=( (first==nodes[l-2]) ? mesh->vertices[previous].ngh[1] : first);
          }
        else {
          status = fe_order2(*mesh, previous, nodes[l-2], rotation);
          nodes[l]=mesh->vertices[previous].ngh[0];
          }
        if(nodes[l] < nodes[0]) goto skip;
/* *----------------------------------------------------------------------------
        to avoid mis-formed cycles*/
//        if(pos(nodes[l],&(nodes[1]),l-1)!=-1) goto skip;
        if(nodes[l]==nodes[0]) {
          for(j=0;j<l;j++) tmp[nelements].push_back(nodes[j]);
          size[nelements]=l;
          cardinal[l]++;
          nelements++;
          goto skip;
          }
        }
skip:
      continue;
      }
    }

  for(m=0;m<nelements;m++) {
    switch (size[m]){
      case 3:
        type[m]=FE_TRIANGLE;
        ntriangles++;
        break;
      case 4:
        type[m]=FE_QUADRANGLE;
        nquadrangles++;
        break;
      case 5:
        type[m]=FE_PENTACLE;
        npentacles++;
        break;
      default:
        type[m]=FE_UNDEFINED;
        break;
      }
    }
  nnpe[FE_TRIANGLE]  =3;
  nnpe[FE_QUADRANGLE]=4;
  nnpe[FE_PENTACLE]  =5;

  for(i=0; i<nndes; i++) flag[i]=-1;

  for(m=0;m<nelements;m++) {
    Polygons::Adjacency cycle;
    for(k=0;k<size[m];k++) {
      flag[tmp[m][k]]=0;
      cycle.push_back(tmp[m][k]);
      }
    elementaryCycles.push_back(cycle);
    }

  for(i=0; i<nndes; i++) {
    if(flag[i]==-1) {
      STDOUT_BASE_LINE_FUNC("unused node %d\n",i);
      }
    }

  for(m=0;m<nelements;m++) {
    if(size[m]>5) {
      printf("oversized element %d: 1st node %d (%d points)\n",m,tmp[m][0], size[m]);
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_SeekCycles(mesh_t *mesh, vector<cycle_t> & ElementaryCycles, int maxsize, bool check_islands, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  Make an element list from nodes' neighbours description.
  
  If consecutive neighbours of a node are connected, they
  form an element with the node.
    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i,k,l,m,n,rotation=-1;
  int max_nb,nndes,*flag;
  int ntriangles,nquadrangles,npentacles,nelements;
  int status;
  int nnpe[100];

  if(mesh->type==-1) mesh->type=0;

  nndes=mesh->nvtxs;
  max_nb=mesh->nnghm;

  for(i=0; i<nndes; i++) {
    status=fe_order(*mesh,i,-1., 0.,rotation);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  search starting point only for nodes not flagged (i.e. set to 0)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  flag=new int[nndes];
  for(i=0; i<nndes; i++) flag[i]=0;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  check and remove for isolated/dead-ending nodes
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  for(i=0; i<nndes; i++) {
    if(mesh->vertices[i].nngh < 2 ) {
      if(debug) printf("%s : pathetic node %d\n",__func__,i);
      flag[i]=1;
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  if no islands expected or wanted, check and remove nodes which are not vertex
  (thus optimizing cycle search)
  
  26/02/2018:
  
  In that case, cycle search is fundamentally buggy, except if mesh has been
  created from 2 different polygon sets, first one being made out of
  monoblock segments (i.e. connections only with polygons from the second set).
  It was the case so far, so the bug has remained harmless.
  
  It has been detected in academic estuary developments
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(!check_islands) {
//     for(i=0; i<nndes; i++) {
    for(i=1; i<nndes; i++) {
      if(mesh->vertices[i].nngh < 3 ) {
        flag[i]=1;
        }
      }
    }
  
//   int flagged=occurence(0,flag,nndes);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  cycle: successive neighbours
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  ntriangles   = 0;
  nquadrangles = 0;
  npentacles   = 0;
  nelements    = 0;

  for(n=0;n<nndes;n++) {
    cycle_t cycle;
    if(flag[n]!=0) continue;
    status=fe_order(*mesh,n,-1., 0.,rotation);
    for(k=0; k<mesh->vertices[n].nngh; k++) {
      int n1,n2;
/*------------------------------------------------------------------------------
      second node, k-th neigbour of n */ 
      n1=mesh->vertices[n].ngh[k];
      if(n1 < n) {
        if(debug) printf("%s : eliminated seed %d %d \n",__func__,n,n1);
        continue;
        }
      if(n1 == n) {
        printf("%s : pathetic nodes %d %d \n",__func__,n,n1);
        continue;
        }
/*------------------------------------------------------------------------------
      third node, left-most neighbour with respect to [n,n1] */ 
      status = fe_order2(*mesh, n1, n, rotation);
      n2=mesh->vertices[n1].ngh[0];
      if(n2 == n) {
        printf("%s : pathetic element %d %d \n",__func__,n,n1);
        continue;
        }
      if(n2 == n1) {
        printf("%s : pathetic nodes %d %d %d\n",__func__,n,n1,n2);
        continue;
        }
      if(n2 < n) {
        if(debug) printf("%s : eliminated seed %d %d %d\n",__func__,n,n1,n2);
        continue;
        }
/*------------------------------------------------------------------------------
      initialize cycle */ 
      cycle.push_back(n);
      cycle.push_back(n1);
      cycle.push_back(n2);
      for(l=3;l<maxsize;l++) {
/*------------------------------------------------------------------------------
        find next node */ 
        int previous=cycle[l-1],next;
        if(mesh->vertices[previous].nngh==2) {
          int first=mesh->vertices[previous].ngh[0];
          next=( (first==cycle[l-2]) ? mesh->vertices[previous].ngh[1] : first);
          }
        else {
          status = fe_order2(*mesh, previous, cycle[l-2], rotation);
          next=mesh->vertices[previous].ngh[0];
          }
        if(next < cycle[0]) {
          if(debug) printf("%s : eliminated cycle %d, n0=%d next=%d size=%d\n",__func__,ElementaryCycles.size(), n,next, cycle.size());
          goto skip;
          }
/*------------------------------------------------------------------------------
        push node in cycle */ 
        cycle.push_back(next);
        if(next==n) {
          ElementaryCycles.push_back(cycle);
          nelements++;
          if(verbose) printf("%s : acquired cycle %d, n0=%d n1=%d n2=%d size=%d\n",__func__,ElementaryCycles.size(), n, n1, n2, cycle.size());
          goto skip;
          }
        }
skip:
      cycle.clear();
      continue;
      }
    }
  
  deletep(&flag);

  for(m=0;m<nelements;m++) {
    switch (ElementaryCycles[m].size()){
      case 3:
//         type[m]=FE_TRIANGLE;
        ntriangles++;
        break;
      case 4:
//         type[m]=FE_QUADRANGLE;
        nquadrangles++;
        break;
      case 5:
//         type[m]=FE_PENTACLE;
        npentacles++;
        break;
      default:
//         type[m]=FE_UNDEFINED;
        break;
      }
    }
  nnpe[FE_TRIANGLE]  =3;
  nnpe[FE_QUADRANGLE]=4;
  nnpe[FE_PENTACLE]  =5;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t> plg_ConvertCycles( const ElementaryCycles & elementaryCycles, const mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  vector<plg_t> set;
  int k;

/* *----------------------------------------------------------------------------
  convert LeBars polygon set into Lyard polygon set*/
  for (ElementaryCycles::const_iterator i = elementaryCycles.begin(); i != elementaryCycles.end();i++) {
    Adjacency p = *i;
    plg_t q;
    const size_t n = p.size()+1;
    q.npt=n;
    q.x=new double[q.npt];
    q.y=new double[q.npt];
    q.t=new double[q.npt];
    q.p=new double[q.npt];
    k=0;
    for (Adjacency::iterator j = p.begin(); j != p.end();j++) {
      const size_t v=*j;
      status=fe_ExportVertex(mesh.vertices[v], q, k);
//       double x=a.x();
//       q.t[k]=a.x();
//       q.p[k]=a.y();
//       q.x[k]=a.x();
//       q.y[k]=a.y();
      k++;
      }
    q.t[k]=q.t[0];
    q.p[k]=q.p[0];
    q.x[k]=q.t[0];
    q.y[k]=q.p[0];
    set.push_back(q);
    }
  return(set);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t> plg_ConvertCycles( const vector<cycle_t> & ElementaryCycles, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  convert LeBars polygons set into Lyard polygons set
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  vector<plg_t> set;
  int k,m;

  for (int i=0;i<ElementaryCycles.size();i++) {
    cycle_t p = ElementaryCycles[i];
    plg_t q;
    const size_t n = p.size();
    q.init(n,PLG_INIT_SEPARATE);
    for (int j = 0; j <p.size();j++) {
      const size_t v=p[j];
      status=fe_ExportVertex(mesh.vertices[v], q, j);
      }
    q.flag=new char[q.npt-1];
    for (int j = 0; j <p.size()-1;j++) {
      const size_t v1=p[j];
      const size_t v2=p[j+1];
      for(k=0;k<mesh.vertices[v1].nedges;k++) {
        bool found;
        m=mesh.vertices[v1].edges[k];
        found=( (mesh.edges[m].extremity[0]==v1 && mesh.edges[m].extremity[1]==v2) ||
                (mesh.edges[m].extremity[0]==v2 && mesh.edges[m].extremity[1]==v1) );
        if(found) break;
        }
      m=fe_isedge(mesh, v1, v2);
/*------------------------------------------------------------------------------
      inherit flag from edge code                                             */
      if(mesh.edges!=0) q.flag[j]=(char) mesh.edges[m].code;
      }
    set.push_back(q);
    }
  return(set);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t> Nei2Polygons(mesh_t & mesh, bool check_islands, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  convert mesh elementary cycles into polygons
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  vector<cycle_t> ElementaryCycles;
  vector<plg_t> p;
  int maxsize=mesh.nvtxs+1;
  
  printf("#build all elementary cycles from mesh (check islands=%d)\n", check_islands);
  status=fe_SeekCycles(&mesh, ElementaryCycles, maxsize, check_islands, verbose, debug);
  
  printf("#export %d cycles as closed polygons\n", ElementaryCycles.size());
  p=plg_ConvertCycles(ElementaryCycles, mesh);

  if(debug) {
    for(int s=0;s<p.size();s++){
      double area=p[s].area();
      printf("polygon %d: area=%lf\n",s, area);
      }
    }
  
  return(p);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_anomalies(mesh_t & mesh, vector<int> & list, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  plg_t p;
  
  if(list.size()==0) return(0);
  
  p.init(list.size(), PLG_INIT_SEPARATE);
  
  for(int k=0; k<list.size(); k++) {
    p.t[k]=mesh.vertices[list[k]].lon;
    p.p[k]=mesh.vertices[list[k]].lat;
    }
  status=plg_save("anomalies.plg", PLG_FORMAT_SCAN, &p, 1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_export2nei(vector<plg_t> & p, mesh_t & mesh, bool strict, bool do_edges, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  convert polygons into mesh, boundary type flag information preserved
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int nn,k,s,size,duplicate,n,nedges,status;
  double d;
  string filename="cycles.nei";
  bool first_is_closed=false;
    
  printf("#export polygons as unstructured mesh\n");

  size=plg_VectorSizeV(p);
  mesh.vertices=new vertex_t[size];
 
  nedges=0;
 
  if(do_edges) {
    size=plg_VectorSizeL(p);
    mesh.edges=new edge_t[size];
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  take all vertices of first polygon : warning, closed polygon special case
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  s=0;
  if(plg_isclosed(p[s])==1) {
    p[s].npt--;
    first_is_closed=true;
    }
  
  for(k=0;k<p[s].npt;k++) {
    status=fe_ImportVertex(mesh.vertices[k], p[s], k);
    }
  mesh.nvtxs=p[s].npt;
  for(k=1;k<p[s].npt;k++) {
    status=fe_connectvertices(mesh, k-1, k);
    if(do_edges) {
      mesh.edges[nedges].extremity[0]=k-1;
      mesh.edges[nedges].extremity[1]=k;
/*------------------------------------------------------------------------------
      export polygon flag as edge code                                        */      
      if(p[s].flag!=0) mesh.edges[nedges].code=(int) p[s].flag[k-1];
      nedges++;
      }
    }
  if(first_is_closed==1) {
    status=fe_connectvertices(mesh, 0, mesh.nvtxs-1);
    if(do_edges) {
      mesh.edges[nedges].extremity[0]=0;
      mesh.edges[nedges].extremity[1]=mesh.nvtxs-1;
/*------------------------------------------------------------------------------
      export polygon flag as edge code                                        */      
      if(p[s].flag!=0) mesh.edges[nedges].code=(int) p[s].flag[0];
      nedges++;
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  scan further polygons
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  n=mesh.nvtxs;
  for(s=1;s<p.size();s++){
    vector<int> added;
    n=mesh.nvtxs;
    for(k=0;k<p[s].npt;k++) {
/*------------------------------------------------------------------------------
      check if vertex already exist */
      duplicate=-1;
      if(strict) {
        duplicate=fe_FindVertex(mesh, p[s].t[k], p[s].p[k], d);
        }
      else {
        if( (k==0) or (k==p[s].npt-1) ) duplicate=fe_FindVertex(mesh, p[s].t[k], p[s].p[k], d);
        else d=1.e+10;
        }
      if(d!=0) {
/*------------------------------------------------------------------------------
        vertex does not exist yet, register and create it */
        status=fe_ImportVertex(mesh.vertices[n], p[s], k);
        added.push_back(n);
        n++;
        mesh.nvtxs=n;
        }
      else {
/*------------------------------------------------------------------------------
        vertex already exists, just register */
        added.push_back(duplicate);
        }
      }
    for(k=1;k<added.size();k++) {
      status=fe_connectvertices(mesh, added[k-1], added[k]);
      if(do_edges) {
        if(nedges==size) {
          printf("%s : nedges=%d size=%d polygon=%d (%d points)\n",__func__,nedges,size,s,p[s].npt);
          }
        mesh.edges[nedges].extremity[0]=added[k-1];
        mesh.edges[nedges].extremity[1]=added[k];
/*------------------------------------------------------------------------------
        export polygon flag as edge code                                      */      
        if(p[s].flag!=0) mesh.edges[nedges].code=(int) p[s].flag[k-1];
        nedges++;
        }
      }
    added.clear();
    mesh.nvtxs=n;
    }

  mesh.nvtxs=n;
  mesh.nedges=nedges;
  if(debug) status=fe_savemesh(filename.c_str(), MESH_FILE_FORMAT_TRIGRID, mesh);
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  check for "dead ends", i.e. vertices inherited from polygons that will not be
  connected at both extremities
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  vector<int> deadends;
  for(n=0;n<mesh.nvtxs;n++) {
    if(mesh.vertices[n].nngh==1) {
      deadends.push_back(n);
      }
    }
  status=plg_anomalies(mesh, deadends, debug);
  
  vector<int> disconnected;
  
  if(deadends.size()!=0) {
    if(debug) printf("warning : found %d deadends, patch applies\n", deadends.size());
    for(k=0;k<deadends.size();k++) {
      n=deadends[k];
      while(mesh.vertices[n].nngh==1) {
        int next=mesh.vertices[n].ngh[0];
        status=fe_disconnectvertex(mesh, n);
        disconnected.push_back(n);
        n=next;
        }
      }
    if(debug) printf("%d nodes disconnected\n", disconnected.size());
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  finalize vertices/edges tables
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(do_edges) {
    for(n=0;n<mesh.nvtxs;n++) {
      mesh.vertices[n].nedges=0;
      }
    for(n=0;n<mesh.nedges;n++) {
      nn=mesh.edges[n].extremity[0];
      mesh.vertices[nn].nedges++;
      nn=mesh.edges[n].extremity[1];
      mesh.vertices[nn].nedges++;
      }
    for(n=0;n<mesh.nvtxs;n++) {
      mesh.vertices[n].edges=new int[mesh.vertices[n].nedges];
      mesh.vertices[n].nedges=0;
      }
    for(n=0;n<mesh.nedges;n++) {
      nn=mesh.edges[n].extremity[0];
      mesh.vertices[nn].edges[mesh.vertices[nn].nedges]=n;
      mesh.vertices[nn].nedges++;
      nn=mesh.edges[n].extremity[1];
      mesh.vertices[nn].edges[mesh.vertices[nn].nedges]=n;
      mesh.vertices[nn].nedges++;
      }
    }
  
  if(debug) status=fe_savemesh(filename.c_str(), MESH_FILE_FORMAT_TRIGRID, mesh);
    
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_CutAtIntersectionsMultiples(const vector<plg_t> & p, vector<plg_t> & q, vector<plg_t> & islands, vector<plg_t> & external, vector<plg_t> & targeted, const vector<plg_t> & limits, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  create a new point at intersection of p and q polygons, if any, and cut polygons
  
  returns uncut, closed polygons in a separate list (islands)
  
  external contains cut polygons of p that will participate to external limits
  
  targetted contains uncut polygons of p that will participate to external limits
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k,status;
  
  printf("#cut polygons at frame limits: p=%d q=%d \n", p.size(), q.size());
  
  for(k=0;k<p.size();k++){
    vector<plg_t> pp;
    plg_t tmpp;
    tmpp.duplicate(p[k]);
    pp.push_back(tmpp);
    if(debug) status=plg_save("debug-selection-interim-0.plg", PLG_FORMAT_SCAN, pp);
    status=plg_CutAtIntersections(pp, q, PLG_SPHERICAL);
/*------------------------------------------------------------------------------
    if nothing was cut, polygon is inside frame, nothing else to do*/
    if(pp.size()==1 && plg_isclosed(pp[0])==1) {
      islands.push_back(pp[0]);
      continue;
      }
    if(debug) status=plg_save("debug-selection-interim-1.plg", PLG_FORMAT_SCAN, pp);
    
/*------------------------------------------------------------------------------
    polygon was cut, removes outside sections if limits given */
    if(limits.size()!=0) {
      status=plg_reduce_selection(pp, limits, PLG_SPHERICAL, debug);
      if(status!=0) return(-1);
      }
    
    if(debug) {
      status=plg_save("debug-selection-interim-2.plg", PLG_FORMAT_SCAN, pp);
      }
    
/*------------------------------------------------------------------------------
    enforce mono-block segments for external*/
    int target[2];
    for (target[0]=0; target[0]< pp.size();target[0]++) {
      if(pp[target[0]].npt==0) continue;
      for (target[1]=0; target[1]< pp.size();target[1]++) {
        if(target[0]==target[1]) continue;
        if(pp[target[1]].npt==0) continue;
        status=plg_concat(target,pp,PLG_SPHERICAL);
        }
      }
    if(debug) {
      status=plg_save("debug-selection-interim-3.plg", PLG_FORMAT_SCAN, pp);
      status=plg_save("debug-limits-1.plg", PLG_FORMAT_SCAN, q);
      }
    
    for(int kk=0;kk<pp.size();kk++) {
      if(pp[kk].npt>1)
        external.push_back(pp[kk]);
      else
        pp[kk].destroy();
      }
    targeted.push_back(p[k]);
    if(debug) status=plg_save("debug-external-interim-2.plg", PLG_FORMAT_SCAN, external);
    }
  
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void plg_smoothlines(plg_t & plg)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k;
  plg_t q;
  
  if(plg.npt<3) return;
  
  int isclosed=plg_isclosed(plg);
  q.init(plg.npt-!isclosed,PLG_INIT_SHARED);
  
  for(k=0;k<plg.npt-1;k++) {
    q.t[k]=0.5*(plg.t[k]+plg.t[k+1]);
    q.p[k]=0.5*(plg.p[k]+plg.p[k+1]);
    }
  
  if(isclosed){
    q.t[k]=q.t[0];
    q.p[k]=q.p[0];
    }
  
  q.setModeFromShared(PLG_INIT_SEPARATE);
  
  plg.destroy();
  plg=q;
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void plg_smoothlines(vector<plg_t> & p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,s;
  
  for(s=0;s<p.size();s++) {
    plg_t *plg=&p[s],q;
    if(plg->npt<3) continue;
    int isclosed=plg_isclosed(*plg);
    q.init(plg->npt-!isclosed,PLG_INIT_SHARED);
    for(k=0;k<plg->npt-1;k++) {
      q.t[k]=0.5*(plg->t[k]+plg->t[k+1]);
      q.p[k]=0.5*(plg->p[k]+plg->p[k+1]);
      }
    if(isclosed){
      q.t[k]=q.t[0];
      q.p[k]=q.p[0];
      }
    
    q.setModeFromShared(PLG_INIT_SEPARATE);
    
    plg->destroy();
    *plg=q;
    }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_RemoveInlandLimits(vector<plg_t> & p, double x, double y, projPJ projection, int position, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int c,s,status;
  int count=0;
  char *keep=new char[p.size()];
  vector<int> continent;
  string output;
    
  for(s=0;s<p.size();s++) keep[s]=1;
  
/*------------------------------------------------------------------------------
  find main boundary */
  for(s=0;s<p.size();s++) {
    int flag=0;
    flag=plg_single(p[s],x,y,&flag,PLG_SPHERICAL);
    if(flag==PLG_POINT_INTERIOR) {
      continent.push_back(s);
      break;
      }
    }
  
  if(continent.size()!=1) {
    return(-1);
    }
  c=continent[0];
  
  if(debug) {
    output="continent.plg";
    status=plg_save(output.c_str(), PLG_FORMAT_SCAN, &p[c], 1);
    }
  
/*------------------------------------------------------------------------------
  first remove inside/exterior polygons (position) */
  for(int k=1;k<p[c].npt;k++) {
    p[c].t[k]=degree_recale(p[c].t[k],p[c].t[k-1]);
    }
  
  double *xx=new double[p.size()];
  double *yy=new double[p.size()];
 
  for(s=0;s<p.size();s++) {
    xx[s]=p[s].t[0];
    yy[s]=p[s].p[0];
    point2D_t point=p[s].InsidePoint();
    xx[s]=point.x;
    yy[s]=point.y;
    geo_to_projection(projection, yy[s], xx[s], &xx[s], &yy[s]);
    }
  
  char *flag=plg_TestInterior(xx, yy, p.size(), &(p[c]), 1, 0, debug);
  
  for(s=p.size()-1;s>=0;s--) {
    if(s==c) continue;
    if(flag[s]==position) {
      keep[s]=0;
      }
    if(keep[s]==0) {
      count++;
      p.erase (p.begin()+s);
      if(s<c) c--;
      }
    }
 
/*------------------------------------------------------------------------------
  then remove imbricated polygons */
  for(int k=0;k<p.size();k++) {
    if(k==c) continue;
    for(int l=0;l<p.size();l++) {
      if(l==k) continue;
      if(l==c) continue;
      point2D_t point=p[l].InsidePoint(PLG_CARTESIAN);
      int flag=0;
      flag=plg_single(p[k],point.x,point.y,&flag, PLG_CARTESIAN, 0);
      if(flag==PLG_POINT_INTERIOR) {
        p.erase (p.begin()+l);
        if(l<c) c--;
        if(l<k) {
          k=-1;
          break;
          }
        else {
          l--;
          }
        }
      }
    }

  delete[] xx;
  delete[] yy;
  
  delete[] flag;
  delete[] keep;
  
  return (0);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int plg_CheckOverlapped(vector<plg_t> & polygons, bool repair, string rootname, bool verbose)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   int repaired=0;
//   
//   for(size_t p=0;p<polygons.size()-1;p++) {
//     frame_t frame=plg_spherical_minmax(&polygons[p], 1);
//     for(size_t q=p+1;p<polygons.size();q++) {
//       
//       }
//     }
//     
//   return (0);
// }
  

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_CheckDuplicated(vector<plg_t> & polygons, bool repair, string rootname, bool verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int repaired=0;
  
  vector<int> *flags=new vector<int>[polygons.size()];
  
// #pragma omp parallel for
  for(size_t p=0;p<polygons.size();p++) {
    
    if(polygons[p].npt==0) continue;
    bool touched=false;
    for(size_t k=0;k<polygons[p].npt-1;k++) {
//       for(size_t l=k+1;l<polygons[p].npt-1;l++) {
      for(size_t l=k+1;l<k+2;l++) {
        if(polygons[p].p[k]==polygons[p].p[l]) {
          double d=plg_distance(polygons[p], k, l, PLG_CARTESIAN);
          int flag=(d==0.0);
          if(flag==1) {
            if(verbose) {
              printf("plg_CheckDuplicated: polygon %d has duplicated points, %d %d\n",p,k,l);
              printf("duplicated position: %lfE %lfN\n",polygons[p].t[k],polygons[p].p[k]);
              }
            if(repair) {
              status=plg_delete_point(&(polygons[p]),l);
              repaired++;
              touched=true;
              if(k==polygons[p].npt-1) break;
              l--;
              }
            else {
              flags[p].push_back(k);
              break;
              }
            }
          }
        }
      }
    if(touched) {
      printf("polygon %4d, id=%4d, npt=%5d: polygon has duplicated points, ",p, polygons[p].id,polygons[p].npt);
      printf("first position: %9.3lfE %9.3lfN\n",polygons[p].t[0],polygons[p].p[0]);
      }
//     else {
//       printf("polygon %d ok\n",p);
//       }
    }
 
  if(repaired!=0) {
    string filename;
    if(rootname!="") filename=rootname+"-repaired.plg";
    else filename="repaired.plg";
    status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, polygons);
    printf("\n %s : repaired set of polygons saved in %s\n\n",__func__,filename.c_str());
    }
  
  status=0;
  for(size_t p=0;p<polygons.size();p++) {
    if(flags[p].size()!=0) {
      status=-1;
      }
    }
  
  delete[]flags;
  
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_CheckDuplicated(plg_t *polygons, int npolygons, bool repair, string rootname, bool verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  vector<plg_t> p;
  
  for(int k=0;k<npolygons;k++) {
    p.push_back(polygons[k]);
    }

  status=plg_CheckDuplicated(p, repair, rootname, verbose);
  
  for(int k=0;k<npolygons;k++) {
    polygons[k]=p[k];
    }
  
  
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_CheckClosure(vector<plg_t> & shorelines_base, bool repair, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, status;
  int step;
  vector<plg_t> p;
  vector<int> list;
  size_t size;
  
  
/*------------------------------------------------------------------------------
  scan closed/unclosed polygons */
  for(int k=0;k<shorelines_base.size();k++) {
    if(plg_isclosed(shorelines_base[k])==0) {
      if(verbose) {
        printf("%s : initial shorelines polygon %d not closed (npt=%d), ", __func__, k,shorelines_base[k].npt);
        n=0;
        printf("t=%9.3lf p=%9.3lf, ",shorelines_base[k].t[n],shorelines_base[k].p[n]);
        n=shorelines_base[k].npt-1;
        printf("t=%9.3lf p=%9.3lf\n",shorelines_base[k].t[n],shorelines_base[k].p[n]);
        }
      list.push_back(k);
      }
    }
  
  if(list.size()!=0) {
    if(verbose) printf("%s, initial shorelines: %d unclosed polygons\n", __func__, list.size());
    }
  else {
    return(0);
    }

  step=0;
  do {
    size=list.size();
    step++;
    if(debug) printf("step %d: %d unclosed polygons\n",step,size);
    if(list.size()>1 && repair) {
/*------------------------------------------------------------------------------
      try to assembly unclosed polygons (iterative process) */
      int merged;
      do {
        merged=0;
        for(int k=0;k<list.size()-1;k++) {
          if(shorelines_base[list[k]].npt==0) continue;
          for(int l=k+1;l<list.size();l++) {
            if(shorelines_base[list[l]].npt==0) continue;
            status=plg_concat(list[k],list[l], shorelines_base, PLG_SPHERICAL);
            if(status==0) {
              if(debug) printf("polygons %d and %d merged\n",list[k],list[l]);
              merged++;
              }
            if(plg_isclosed(shorelines_base[list[k]])==1) {
              break;
              }
            }
          }
/*------------------------------------------------------------------------------
        scan closed/unclosed polygons again */
        list.clear();
        for(int k=0;k<shorelines_base.size();k++) {
          if(shorelines_base[k].npt==0) continue;
          if(plg_isclosed(shorelines_base[k])==0) {
            list.push_back(k);
            }
          }
        if(list.size()<2) break;
        } while(merged!=0);
      }
    
//     if(verbose) printf("%s, merged segments : %d \n", __func__, merged);
    
/*------------------------------------------------------------------------------
    remove empty polygons */
    status=plg_compact_entries(shorelines_base);

/*------------------------------------------------------------------------------
    scan closed/unclosed polygons again */
    list.clear();
    for(int k=0;k<shorelines_base.size();k++) {
//       if(plg_isclosed(shorelines_base[k])==0) {
//         double d=plg_distance(shorelines_base[k], 0, shorelines_base[k].npt-1, PLG_SPHERICAL);
//         if(d<0.010) {
//           plg_point_t point;
//           point.t=shorelines_base[k].t[0];
//           point.p=shorelines_base[k].p[0];
//           point.x=shorelines_base[k].x[0];
//           point.y=shorelines_base[k].y[0];
//           status= plg_insert_point(shorelines_base[k],shorelines_base[k].npt, point);
//           }
//         }
      if(plg_isclosed(shorelines_base[k])==0) {
//         if(debug) printf("repaired shorelines polygon %d not closed (%d points)\n",k,shorelines_base[k].npt);
        list.push_back(k);
        }
      }
    } while(list.size()!=size);
  
  
/*------------------------------------------------------------------------------
  try to close unclosed polygons */
  do {
    size=list.size();
    step++;
    if(debug) printf("step %d: %d unclosed polygons\n",step,size);

    list.clear();
    for(int k=0;k<shorelines_base.size();k++) {
      if(plg_isclosed(shorelines_base[k])==0) {
        double d=plg_distance(shorelines_base[k], 0, shorelines_base[k].npt-1, PLG_SPHERICAL);
        if(d<0.010) {
          status=close_plg(&shorelines_base[k]);
//          if(debug) printf("repaired shorelines polygon %d not closed (%d points)\n",k,shorelines_base[k].npt);
          }
        }
      if(plg_isclosed(shorelines_base[k])==0) {
        list.push_back(k);
        }
      }
    } while(list.size()!=size);

  if(debug) status=plg_save("debug-repaired.plg", PLG_FORMAT_SCAN, shorelines_base);
   
  if(verbose) printf("%s, after repairing : %d unclosed polygons\n", __func__, list.size());
  
  if(list.size()!=0) {
    if(debug) {
      for(int k=0;k<list.size();k++) printf("%d ",list[k]);
      printf("\n");
      }
    status=-1;
    }
  else {
    status=0;
    }
  
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t> plg_extract(vector<plg_t> & shorelines_base, const vector<plg_t> & selection, const string & rootname, int coordinates, int target, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  extract the polygons from a coastline within a selected area

  rootname debug root file name. If empty and if debug==false then nothing is saved
  
  mode   : PLG_SPHERICAL or PLG_CARTESIAN workspace
  target : if 0 then return shorelines else return modelling polygone

  return the output selected by \c target

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  string filename;
  vector<plg_t> limits;
  int flag;
  vector<plg_t> islands, external, modeling, shorelines, empty;
 
  if(selection.size()==0) {
    return(empty);
    }
  
  filename=rootname+"-selection.plg";
  if(rootname!="") status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, selection);

  status=plg_CheckClosure(shorelines_base, false, 1, debug);
  if(status!=0)
    return (empty);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  sample shorelines database inside given selection polygon
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  bool strict=true;
  status=plg_PolygonReduction(shorelines_base, selection, strict, coordinates, 0, debug);
  if(debug) {
    filename=rootname+"-working-selection-0.plg";
    status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, shorelines_base);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  cut shorelines and selection polygon; targeted are original shorelines actually used
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  vector<plg_t> targeted;
  plg_t tmp;
  tmp.duplicate(selection[0]);
  
  limits.push_back(tmp);

  status=plg_CutAtIntersectionsMultiples(shorelines_base, limits, islands, external, targeted, selection, debug);
  
  if(status!=0) return(external);

  for(int k=0;k<targeted.size();k++) {
    if(plg_isclosed(targeted[k])==0) {
      printf("%s: polygon %d not closed\n",__func__,k);
      }
    }

  if(debug) {
    filename=rootname+"-targeted.plg";
    status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, targeted);
    filename=rootname+"-limits.plg";
    status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, limits);
    filename=rootname+"-islands.plg";
    status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, islands);
    filename=rootname+"-external.plg";
    status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, external);
    }
 
  for(int k=0;k<limits.size();k++)     limits[k].SetFlag('M');
  for(int k=0;k<islands.size();k++)   islands[k].SetFlag('T');
  for(int k=0;k<external.size();k++) external[k].SetFlag('T');

/*------------------------------------------------------------------------------
  save shorelines and frame pieces*/
  if(debug) {
    vector<plg_t> raw;
    status=plg_merge(raw, external);
    status=plg_merge(raw, limits);
    status=plg_merge(raw, islands);
    filename=rootname+"-shorelines-raw.plg";
    status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, raw);
    raw.clear();
    }
  
  if(target==2) {
    plg_merge(external, islands);
    return (external);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  assembly external shorelines and frame in closed polygons (elementary cycles)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=plg_merge(external, limits);
  if(debug) {
    filename=rootname+"-external.plg";
    status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, external);
    }
  
  vector<plg_t> r;
  if(external.size()!=1) {
    mesh_t mesh;
    bool strict=false, do_edges=true;
    status=plg_export2nei(external, mesh, strict, do_edges, debug);
    bool check_islands=false;
    int verbose=0;
    r=Nei2Polygons(mesh, check_islands, verbose, debug);
    if(debug) {
      filename=rootname+"-elementary-cycles.plg";
      status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, r);
      }
    mesh.destroy();
    }
  else {
    r=plg_duplicate(limits);
    }
  
  /* this will also destroy entries in limits*/
  plg_destroy_entries(external);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  select model limits ("wet region") external polygons, test based on original
  shorelines (assumes dry land is contained inside closed polygons)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  bool *destroyR=aset(r.size(),true);

  frame_t frame=plg_spherical_minmax(selection);
  status=plg_recale(targeted, frame, PLG_SPHERICAL);
  for(int k=0;k<r.size();k++) {
    if(r[k].area()<0) continue;
/*------------------------------------------------------------------------------
    infer a point inside polygon r[k] */
    point2D_t point=r[k].InsidePoint();
    bool keep=true;
    for(int l=0;l<targeted.size();l++) {
      int check=0;
      flag=0;
      flag=plg_single(targeted[l], point.x, point.y, &flag, PLG_SPHERICAL, check);
      if(flag==PLG_POINT_INTERIOR) {
/*------------------------------------------------------------------------------
        lie inside a shoreline polygon, hence "dry", discard */
        keep=false;
        break;
        }
      }
    if(keep) {
      modeling.push_back(r[k]);
      destroyR[k]=false;
      }
    }

  status=plg_merge(modeling, plg_duplicate(islands));
  if(rootname!="") {
    filename=rootname+"-model.plg";
    status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, modeling);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  select shorelines limits ("dry region") external polygons, test based on
  original shorelines
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(int k=0;k<r.size();k++) {
    if(r[k].area()<0) continue;
/*------------------------------------------------------------------------------
    infer a point inside polygon r[k] */
    point2D_t point=r[k].InsidePoint();
    if(debug) printf("%d %lf %lf\n", k, point.x, point.y);
    bool keep=false;
    for(int l=0;l<targeted.size();l++) {
      int check=0;
      flag=0;
      flag=plg_single(targeted[l], point.x, point.y, &flag, PLG_SPHERICAL, check);
      if(flag==PLG_POINT_INTERIOR) {
        keep=!keep;
        }
      }
    if(keep){
/*------------------------------------------------------------------------------
      odd inclusion in shoreline polygon,  ??? */
      shorelines.push_back(r[k]);
      destroyR[k]=false;
      }
    }
  
  for(int k=0;k<r.size();k++)
    if(destroyR[k])
      r[k].destroy();
  
  deletep(&destroyR);
  
//   if(rootname!="") {
//     filename=rootname+"-shorelines-no-islands.plg";
//     status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, shorelines);
//     }

  status=plg_merge(shorelines, islands);
  if(rootname!="") {
    filename=rootname+"-shorelines.plg";
    status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, shorelines);
    }
//   legend_t   *legends=NULL;
//   legends=lgd_import(x, y, buf, ndata, ncols);
//   status=lgd_save("topo-compare.lgd", legends, 1, NULL, NULL);
  
/*------------------------------------------------------------------------------
  shorelines contains the entries of islands */
  if(target==0) {
    plg_destroy_entries(modeling);
    return (shorelines);
    }
  else {
    plg_destroy_entries(shorelines);
    return (modeling);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_SortBySize(vector<plg_t> & p, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  size_t size=p.size();
  
  double *area=new double[p.size()];
  
  for(int k=0;k<size;k++) {
    area[k]=fabs(p[k].area());
    }
  
  for(int k=0;k<size-1;k++) {
    pos=maxpos(&(area[k]), size-k);
    if(pos!=0) {
//       double a=area[k];
//       double b=area[pos];
      plg_t tmp_p=p[k];
      p[k]=p[k+pos];
      p[k+pos]=tmp_p;
      double tmp_a=area[k];
      area[k]=area[k+pos];
      area[k+pos]=tmp_a;
      }
    }

  if(verbose==1) {
    for(int k=0;k<p.size();k++) {
      double a=fabs(p[k].area());
      printf("polygon %d: area=%lf\n",k,a);
      }
    }

  delete[] area;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t> plg_extract(const string & ShorelinesFile, const vector<plg_t> & selection, const string & rootname, int coordinates, int target, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  sample shorelines database inside given frame
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  vector<plg_t> shorelines_base;
  vector<plg_t> extraction;
  size_t size;
  string tmp;
  
//   debug=true;
 

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  load shorelines and perform sanity checks
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  int inFormat=plg_find_format(ShorelinesFile);
  status=plg_load(ShorelinesFile, inFormat, shorelines_base);
  if(status!=0) return(extraction);
  
  status=plg_CheckDuplicated(shorelines_base, true, rootname, debug);
 
  size=shorelines_base.size();

  status=plg_CheckClosure(shorelines_base, true, 0, debug);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  reduce shorelines database to optimize extraction
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  bool strict=true;
  status=plg_PolygonReduction(shorelines_base, selection, strict, coordinates, 1, true);
  if(shorelines_base.size()==0) return(extraction);
  
  if(debug) {
    string filename=rootname+"-base-shorelines-reduced.plg";
    status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, shorelines_base);
    }
  status=plg_CheckDuplicated(shorelines_base, true, rootname, debug);
  
  size=shorelines_base.size();

  status=plg_CheckClosure(shorelines_base, true, 1, debug);
  
  if(status!=0) {
    string filename=rootname+"-base-shorelines-reduced.plg";
    status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, shorelines_base);
    return (extraction);
    }
  
  if(shorelines_base.size()!=size) {
    string output;
    size_t pos=ShorelinesFile.rfind("/");
    if(pos!=string::npos) output="repaired-"+ShorelinesFile.substr(pos+1);
    else output="repaired-"+ShorelinesFile;
    status=plg_save(output.c_str(), inFormat, shorelines_base);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  extract shorelines database
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//  status=plg_SortBySize(shorelines_base);
  
  extraction=plg_extract(shorelines_base, selection, rootname, coordinates, target, debug);
  
  plg_destroy_entries(shorelines_base);
  
  return (extraction);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_QuickFrameDeletion(vector<plg_t> & polygons, vector<plg_t> & outside, vector<plg_t> & intersecting, frame_t frame, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  split polygon set into outside and possibly intersecting sets
  
  now based on polygon copy instead of addressing
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int intersection, status;
  int i,j;
  double x,y;
  bool ok;
  
  bool *keep=new bool[polygons.size()];
  double center=frame.x_center();
  
  for(i=0;i<polygons.size();i++) {
    intersection=0;
    for(j=0;j<polygons[i].npt;j++) {
      x=polygons[i].x[j];
      y=polygons[i].y[j];
      if(mode==PLG_SPHERICAL) x=degree_recale(x, center);
      ok=frame.inside(x, y);
      if(ok) {
        intersection++;
        break;
        }
      }
    if(intersection!=0) {
/*------------------------------------------------------------------------------
      at least one point inside frame, could be intersecting                  */
      status=plg_add_entry(intersecting, polygons[i]);
      keep[i]=false;
      }
    else {
/*------------------------------------------------------------------------------
      no point inside frame, could be outside; not bullet-proof test          */
      status=plg_add_entry(outside, polygons[i]);
      keep[i]=true;
      }
    }
  
//   for(i=0;i<polygons.size();i++) {
//     if(keep[i]) {
//       status=plg_add_entry(outside, polygons[i]);
//       }
//     }
    
  status=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_SplitPolygons(vector<plg_t> & polygons, const vector<plg_t> & selection, vector<plg_t> & outside, vector<plg_t> & inside, vector<plg_t> & intersecting, bool strict, int coordinates)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  split polygon set into inside, outside and intersecting sets
  
  now based on polygon copy instead of addressing
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i, j, flag, intersection, status;
  int check=0;
  int is_inside;
  frame_t frame;
  double x,y;
  vector<plg_t>  temporary;

/*------------------------------------------------------------------------------
  do first quick selection based on min/max frame                             */
  frame= plg_spherical_minmax(selection);
  status=plg_QuickFrameDeletion(polygons, outside, temporary, frame, coordinates);
  
/*------------------------------------------------------------------------------
  then remove shorelines with no point inside selection polygons              */
  double center=frame.x_center();
  for(i=0;i<temporary.size();i++) {
    intersection=0;
    for(j=0;j<temporary[i].npt;j++) {
      x=temporary[i].x[j];
      y=temporary[i].y[j];
      if(coordinates==PLG_SPHERICAL) x=degree_recale(x, center);
      is_inside=0;
      flag=plg_single(selection[0], x, y, &is_inside, coordinates, check);
      switch (strict) {
        case true:
          if(flag==PLG_POINT_INTERIOR) {
            intersection++;
            }
          break;
        case false:
          if(flag!=PLG_POINT_EXTERIOR) {
            intersection++;
            }
          break;
        }
      }
    if(intersection==temporary[i].npt) {
/*------------------------------------------------------------------------------
      all polygon points inside selection polygon, i.e. is inside             */
      status=plg_add_entry(inside, temporary[i]);
      }
    else if(intersection==0) {
/*------------------------------------------------------------------------------
      all polygon points outside selection polygon, i.e. is outside           */
      status=plg_add_entry(outside, temporary[i]);
      }
    else {
/*------------------------------------------------------------------------------
      mix of points inside and outside selection polygon, i.e. is intersecting */
      status=plg_add_entry(intersecting, temporary[i]);
      }
    }
  
  plg_destroy_entries(temporary);
  temporary.clear();
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_DeleteSelection(vector<plg_t> & polygons, const vector<plg_t> & selection, bool strict, int coordinates)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  wrapper
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  vector<plg_t> inside, intersecting,outside;
  
  status=plg_SplitPolygons(polygons, selection, outside, inside, intersecting, strict, coordinates);
  
  status=plg_destroy_entries(inside);
  status=plg_destroy_entries(intersecting);
  status=plg_destroy_entries(polygons);
  
  polygons=outside;
 
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t> plg_extract_out(vector<plg_t> & shorelines_base, const vector<plg_t> & selection, string rootname, int mode, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
  reduce shorelines database by removing lines inside given selection polygons
   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  vector<plg_t> extraction;
  size_t size;
  
  vector<plg_t> outside, inside, intersecting, limits, external;
  vector<plg_t> shorelines;
 
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
  
  cut base polygons intersecting exclusion polygon

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  printf("#cut intersecting polygons (%d) at selection limits\n", intersecting.size());
 
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
/*------------------------------------------------------------------------------
    enforce mono-block segments for external*/
    int target[2];
    for (target[0]=0; target[0]< pp.size();target[0]++) {
      if(pp[target[0]].npt==0) continue;
      for (target[1]=0; target[1]< pp.size();target[1]++) {
        if(target[0]==target[1]) continue;
        if(pp[target[1]].npt==0) continue;
        status=plg_concat(target,pp,PLG_SPHERICAL);
        }
      }
    for(int kk=0;kk<pp.size();kk++) {
      if(pp[kk].npt>1) external.push_back(pp[kk]);
      }
    pp.clear();
    }
  
  for(int k=0;k<limits.size();k++)     limits[k].SetFlag('M');
  for(int k=0;k<outside.size();k++)   outside[k].SetFlag('T');
  for(int k=0;k<external.size();k++) external[k].SetFlag('T');
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  seek all polygons that can be formed by external+limits

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
    
  status=plg_merge(external, limits);
  
  vector<plg_t> r;
  if(external.size()!=1) {
    mesh_t mesh;
    bool strict=false, do_edges=true;
    status=plg_export2nei(external, mesh, strict, do_edges, debug);
    bool check_islands=false;
    int verbose=0;
    r=Nei2Polygons(mesh, check_islands, verbose, debug);
    mesh.destroy();
    }
  else {
    r=plg_duplicate(limits);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  select shorelines limits ("dry region") external polygons, test based on
  original model limits
  
  (assumes wet is contained inside closed polygons)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  bool *destroyR=aset(r.size(),true);

  frame_t frame=plg_spherical_minmax(selection);
  status=plg_recale(intersecting, frame, PLG_SPHERICAL);
  
  for(int k=0;k<r.size();k++) {
    if(r[k].area()<0) continue;
/*------------------------------------------------------------------------------
    infer a point inside polygon r[k] */
    point2D_t point=r[k].InsidePoint();
    if(debug) printf("%d %lf %lf\n", k, point.x, point.y);
    bool keep=false;
    for(int l=0;l<intersecting.size();l++) {
      int check=0;
      int flag=0;
      flag=plg_single(intersecting[l], point.x, point.y, &flag, PLG_SPHERICAL, check);
      if(flag==PLG_POINT_INTERIOR) {
        keep=!keep;
        }
      }
    if(keep){
/*------------------------------------------------------------------------------
      odd inclusion in shoreline polygon, ??? */
      shorelines.push_back(r[k]);
      destroyR[k]=false;
      }
    }
  
  for(int k=0;k<r.size();k++)
    if(destroyR[k])
      r[k].destroy();
  
  
  status=plg_merge(shorelines, outside);
  
  return (shorelines);
}


