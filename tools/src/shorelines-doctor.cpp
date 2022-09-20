
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
#include "geo.h"
#include "archive.h"
#include "netcdf-proto.h"
#include "poc-time.h"
#include "functions.h"

#include "polygons.hpp"
#include "exceptions.hpp"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int defgrid(grid_t & grid, vector<plg_t> & polygons, double step, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  return a cartesian grid at the required resolution
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
  frame_t frame;
  projPJ projection;
  double ref_lat,ref_lon,x,y;
  char *pj_parameters=new char[512];

  frame= plg_recale(polygons);

  ref_lat=(frame.ymax+frame.ymin)/2;
  ref_lon=(frame.xmax+frame.xmin)/2;
  
/*------------------------------------------------------------------------------
  define stereographic oblique projection */
  projection=assign_StereoOblique(ref_lat,ref_lon,pj_parameters);
  
  geo_to_projection(projection, ref_lat, ref_lon, &x, &y);
    
  status=plg_cartesian(projection,polygons);
    
  frame= plg_cartesian_minmax(polygons);
    
  y=(frame.ymax+frame.ymin)/2;
  x=(frame.xmax+frame.xmin)/2;
    
  projection_to_geo(projection, &ref_lat, &ref_lon, x, y);

  printf("\ndefine background working grid: resolution=%lf m, projection=%s\n", step, pj_parameters);
  
  frame.dilatation(0.05);

  grid.xmin=frame.xmin;
  grid.xmax=frame.xmax;

  grid.ymin=frame.ymin;
  grid.ymax=frame.ymax;

  grid.dx=step;
  grid.dy=step;

  grid.xmin=NINT(grid.xmin/grid.dx)*grid.dx-grid.dx;
  grid.xmax=NINT(grid.xmax/grid.dx)*grid.dx+grid.dx;

  grid.ymin=NINT(grid.ymin/grid.dy)*grid.dy-grid.dy;
  grid.ymax=NINT(grid.ymax/grid.dy)*grid.dy+grid.dy;

  grid.nx=(int) NINT((grid.xmax-grid.xmin)/grid.dx)+1;
  grid.ny=(int) NINT((grid.ymax-grid.ymin)/grid.dy)+1;
  
  grid.modeH=0;

  status=map_completegridaxis(&grid,1);
  if(status!=0) return(status);

  grid.projection=projection;
  grid.proj4options=poc_strdup(pj_parameters);
  
  if(verbose==1) {
    map_printgrid(grid);
    printf("\n");
    }
    
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int FootPrint(vector<plg_t> & polygons, double step, string rootname, grid_t & grid, vector<size_t>* & pixels, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  pixelize
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
  vector<size_t> active;
  int *passed;
  int verbose=0;
  
  frame_t frame=plg_spherical_minmax(polygons);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  define the footprint master grid
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  status=defgrid(grid, polygons, step, verbose);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  detect polygons pixel
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  pixels=new vector<size_t> [polygons.size()];
  
  passed=new int[grid.Hsize()];
  
  for(size_t m=0; m<grid.Hsize(); m++) passed[m]=0;
    
  for(size_t p=0;p<polygons.size();p++) {
    for(size_t k=0;k<polygons[p].npt;k++) {
      double x=polygons[p].x[k];
      double y=polygons[p].y[k];
      int i,j;
      size_t m, pos;
      status=map_index(grid, x, y, &i, &j);
      m=j*grid.nx+i;
      passed[m]=1;
      pos=vpos(m,pixels[p]);
      if (pos==-1) {
        pixels[p].push_back(m);
        }
//       pos=vpos(m,active);
//       if (pos==-1) {
//         active.push_back(m);
//         }
      }
    if(debug) printf("%s : polygon %d #pixels=%d\n", __func__, p, pixels[p].size());
    }
  
  delete[] passed;
  
  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int vcommon(const vector<T> & p, const vector<T> & q, vector<T> & r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  for(size_t k=0;k<p.size(); k++) {
    size_t pos;
    pos=vpos(p[k],q);
    if (pos!=-1) {
      r.push_back(p[k]);
      }
    }
  
  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int vadd(vector<T> & p, T m)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  size_t pos;
  pos=vpos(m,p);
  if(pos==-1) {
    p.push_back(m);
    }
  
  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_reorganize(plg_t & p, plg_t q, vector<plg_t> & polygons, int p1, int p2, int q1, int q2, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
 
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
  
  vector<int> limitP;
  vector<int> limitQ;
  
  limitP.push_back(p1);
  limitP.push_back(p2);
  
  limitQ.push_back(q1);
  limitQ.push_back(q2);
  
  plg_t workP[2], workQ[2];
  bool closedP=plg_isclosed(p), closedQ=plg_isclosed(q);
  
  if(limitP[0]<limitP[1]) {
    if(limitP[0]!= 0)       workP[0].duplicate(p, 0, limitP[0]); 
    if(limitP[1]!= p.npt-1) workP[1].duplicate(p, limitP[1], p.npt-1); 
    }
  else  {
    if(limitP[0]!= p.npt-1) workP[0].duplicate(p, limitP[0], p.npt-1); 
    if(limitP[1]!= 0)       workP[1].duplicate(p, 0,limitP[1]); 
    }
  
  if(limitQ[0]<limitQ[1]) {
    if(not closedQ) {
      if(limitQ[0]!= 0)       workQ[0].duplicate(q, 0, limitQ[0]); 
      if(limitQ[1]!= q.npt-1) workQ[1].duplicate(q, limitQ[1], q.npt-1);
      }
    else {
      if(limitQ[0]!= 0)       workQ[0].duplicate(q, 0, limitQ[0]); 
      if(limitQ[1]!= q.npt-1) workQ[1].duplicate(q, limitQ[1], q.npt-1);
      }
    }
  else  {
    if(limitQ[0]!= q.npt-1) workQ[0].duplicate(q, limitQ[0], q.npt-1);
    if(limitQ[1]!= 0)       workQ[1].duplicate(q, 0, limitQ[1]);
    }
  
  if(workQ[0].npt!=0) {
    status=plg_concat(workP[0],workQ[0]);
    }
    
  if(workQ[1].npt!=0) {
    status=plg_concat(workP[1],workQ[1]);
    }
  
  polygons.push_back(workP[0]);
  polygons.push_back(workP[1]);
//   polygons[p].duplicate(workP[0]);
//   polygons[q].duplicate(workP[1]);

  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_CheckOverlapped(vector<plg_t> & polygons, bool join, bool split, bool overlapped, string rootname, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  detect and fix overlapping polygons
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
  int repaired=0;
  vector<size_t> *pixels;
  size_t count, bkp;
  double step;
  grid_t grid;
  int solved=0, unsolved=0;
  vector<plg_t> anomalies;
  vector<size_t> tmp;
  string s;

  bkp=polygons.size();
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  get polygons pixels
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  step=5.e+03;
  status=FootPrint(polygons, step, rootname, grid, pixels, debug);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  scan polygons
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  count=0;
  for(size_t p=0;p<polygons.size()-1;p++) {
    if(p>=bkp) break;
    if(polygons[p].npt==0) continue;
    for(size_t q=p+1;q<polygons.size();q++) {
      vector<size_t> common;
      vector<int> foundP, foundQ;
      if(q>=bkp) break;
      if(polygons[q].npt==0) continue;
/*------------------------------------------------------------------------------
      check if p and q have common pixels */
      status=vcommon(pixels[p],pixels[q],common);
      if(common.size()==0) continue;
      if(debug) printf("%s : polygons %d %d #pixels=%d\n", __func__, p, q, common.size());
      count++;
/*------------------------------------------------------------------------------
      check if p and q overlapp */
      for(int k=0;k<polygons[p].npt;k++) {
        size_t pos,posX,posY;
        int i,j;
        size_t m;
        status=map_index(grid, polygons[p].x[k], polygons[p].y[k], &i, &j);
        m=j*grid.nx+i;
        pos=vpos(m, common);
        if(pos==-1) continue;
        posX=vpos(polygons[p].x[k], polygons[q].x, polygons[q].npt);
        if(posX==-1) continue;
        posY=vpos(polygons[p].y[k], polygons[q].y, polygons[q].npt);
        if(posY!=posX) continue;
        foundP.push_back(k);
        foundQ.push_back(posX);
        if(debug) printf("%s : polygons %d %d k %d pos %d %d\n", __func__, p, q, k, posX, posY);
        }
      if(foundQ.size()==0) continue;
      if(foundQ.size()==1) {
        if( (foundQ[0]==0 or foundQ[0]==polygons[p].npt) and (foundP[0]==0 or foundP[0]==polygons[p].npt) and join) {
/*------------------------------------------------------------------------------
          1 common point only at extremities, concat p and q */
          printf("%s : polygons %d (%d) %d (%d) #connectable %d, %d %d \n", __func__, p, polygons[p].npt, q, polygons[q].npt, foundQ.size(), foundP[0], foundQ[0]);
          status=plg_concat(polygons[p],polygons[q]);
          polygons[q].npt=0;
          solved++;
          break;
          }
//         printf("%s : polygons %d (%d) %d (%d) #connectable %d, %d %d \n", __func__, p, polygons[p].npt, q, polygons[q].npt, foundQ.size(), foundP[0], foundQ[0]);
        unsolved++;
        status=vadd(tmp,p);
        status=vadd(tmp,q);
        continue;
        }
      if(foundQ.size()==2) {
        if(abs(foundP[0]-foundP[1])==1 and abs(foundQ[0]-foundQ[1])==1 and split) {
/*------------------------------------------------------------------------------
          channel connecting on "main" river */
          printf("%s : polygons %d (%d) %d (%d) #connectable %d, %d %d, %d %d \n", __func__, p, polygons[p].npt, q, polygons[q].npt, foundQ.size(), foundP[0], foundP[1], foundQ[0], foundQ[1]);
          status=plg_reorganize(polygons[p], polygons[q], polygons, foundP[0], foundP[1], foundQ[0], foundQ[1], debug);
          polygons[p].npt=0;
          polygons[q].npt=0;
          solved++;
          break;
          }
        unsolved++;
        status=vadd(tmp,p);
        status=vadd(tmp,q);
        continue;
        }
      if(foundQ.size()==3) {
        bool closedP=plg_isclosed(polygons[p]), closedQ=plg_isclosed(polygons[q]);
        bool process=(closedP and foundP[0]==0 and foundP[2]==polygons[p].npt-1) or (closedQ and foundQ[0]==0 and foundQ[2]==polygons[q].npt-1);
        if(abs(foundP[0]-foundP[1])==1 and abs(foundQ[0]-foundQ[1])==1 and process and split) {
/*------------------------------------------------------------------------------
          channel connecting on "main" river, special case */
          printf("%s : polygons %d (%d) %d (%d) #connectable %d, %d %d, %d %d \n", __func__, p, polygons[p].npt, q, polygons[q].npt, foundQ.size(), foundP[0], foundP[1], foundQ[0], foundQ[1]);
          status=plg_reorganize(polygons[p], polygons[q], polygons, foundP[0], foundP[1], foundQ[0], foundQ[1], debug);
          polygons[p].npt=0;
          polygons[q].npt=0;
          solved++;
          break;
          }
        if(abs(foundP[1]-foundP[2])==1 and abs(foundQ[1]-foundQ[2])==1 and process and split) {
/*------------------------------------------------------------------------------
          channel connecting on "main" river, special case */
          printf("%s : polygons %d (%d) %d (%d) #connectable %d, %d %d, %d %d \n", __func__, p, polygons[p].npt, q, polygons[q].npt, foundQ.size(), foundP[1], foundP[2], foundQ[1], foundQ[2]);
          status=plg_reorganize(polygons[p], polygons[q], polygons, foundP[1], foundP[2], foundQ[1], foundQ[2], debug);
          polygons[p].npt=0;
          polygons[q].npt=0;
          solved++;
          break;
          }
//         if(abs(foundP[0]-foundP[2])==2 and abs(foundQ[0]-foundQ[2])==2 and split) {
//           printf("%s : polygons %d (%d) %d (%d) #connectable %d, %d %d %d, %d %d %d \n", __func__, p, polygons[p].npt, q, polygons[q].npt, foundQ.size(), 
//                                                                                         foundP[0], foundP[1], foundP[2], foundQ[0], foundQ[1], foundQ[2]);
//           status=plg_reorganize(polygons[p], polygons[q], polygons, foundP[0], foundP[2], foundQ[0], foundQ[2], debug);
//           polygons[p].npt=0;
//           polygons[q].npt=0;
//           solved++;
//           break;
//           }
//         printf("%s : polygons %d (%d) %d (%d) #connectable %d, %d %d %d, %d %d %d \n", __func__, p, polygons[p].npt, q, polygons[q].npt, foundQ.size(), 
//                                                                                                  foundP[0], foundP[1], foundP[2], foundQ[0], foundQ[1], foundQ[2]);
        unsolved++;
        status=vadd(tmp,p);
        status=vadd(tmp,q);
        continue;
        }
      if(foundQ.size()==polygons[p].npt and foundQ.size()==polygons[q].npt and overlapped) {
/*------------------------------------------------------------------------------
        duplicated polygon */
        printf("%s : polygons %d (%d) %d (%d) #identical %d \n", __func__, p, polygons[p].npt, q, polygons[q].npt, foundQ.size());
        polygons[q].npt=0;
        solved++;
        continue;
        }
      if(foundQ.size()==polygons[p].npt xor foundQ.size()==polygons[q].npt and overlapped) {
/*------------------------------------------------------------------------------
        polygon repeating part of an existing one*/
        printf("%s : polygons %d (%d) %d (%d) #overlapped %d \n", __func__, p, polygons[p].npt, q, polygons[q].npt, foundQ.size());
        if(foundQ.size()==polygons[q].npt) polygons[q].npt=0;
        if(foundQ.size()==polygons[p].npt) {
          polygons[p].npt=0;
          break;
          }
        solved++;
        continue;
        }
        
//       if(foundQ.size()==polygons[p].npt xor foundQ.size()==polygons[q].npt and overlapped) {
// /*------------------------------------------------------------------------------
//         polygons repeating parts of each other*/
//         printf("%s : polygons %d (%d) %d (%d) #overlapped %d \n", __func__, p, polygons[p].npt, q, polygons[q].npt, foundQ.size());
//         if(foundQ.size()==polygons[q].npt) polygons[q].npt=0;
//         if(foundQ.size()==polygons[p].npt) {
//           polygons[p].npt=0;
//           break;
//           }
//         solved++;
//         continue;
//         }
        
//       bool closedP=plg_isclosed(polygons[p]), closedQ=plg_isclosed(polygons[q]);
//       bool process=(not closedP and not closedQ);
//       if(process and overlapped) {
//         int k=1;
//         while (abs(foundP[k]-foundP[k-1])==1 and abs(foundQ[k]-foundQ[k-1])==1) {
//           if(k==foundP.size()-1) break;
//           k++;
//           }
//         solved++;
//         continue;
//         }
//       printf("%s : polygons %d %d #overlapped %d \n", __func__, p, q, found.size());
      unsolved++;
      status=vadd(tmp,p);
      status=vadd(tmp,q);
      }
    }

  s=rootname+"-unsolved.plg";
  for(int k=0;k<tmp.size();k++) anomalies.push_back(polygons[tmp[k]]);
  status=plg_save(s.c_str(), PLG_FORMAT_SCAN, anomalies);
  
  for(size_t p=0;p<polygons.size();p++) {
    if(polygons[p].npt<=1) {
      polygons.erase (polygons.begin()+p);
      p--;
      }
    }
    
  for(int k=0;k<bkp;k++) pixels[k].clear();
  delete[]  pixels;
  grid.free();
  
  printf("%s : polygons size before=%d after=%d (solved=%d unsolved=%d)\n", __func__,bkp,polygons.size(), solved, unsolved); 
  
  return (solved);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_CheckSelfOverlapped(vector<plg_t> & polygons, bool repair, string rootname, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  detect and fix self-overlapping polygons
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
  int repaired=0;
  vector<size_t> *pixels;
  size_t count, bkp;
  double step;
  grid_t grid;
  int solved=0, unsolved=0;
  vector<plg_t> anomalies;
  plg_t work;
  
  bkp=polygons.size();
    
  count=0;
  for(size_t p=0;p<polygons.size();p++) {
    if(polygons[p].npt==0) continue;
    double x=polygons[p].x[0];
    double y=polygons[p].y[0];
    bool closed=plg_isclosed(polygons[p]);
    int npt;
    if(closed)
      npt=polygons[p].npt-2;
    else
      npt=polygons[p].npt-1;
    size_t posX,posY;
    posX=vpos(x, &(polygons[p].x[1]), npt);
    if(posX==-1) continue;
    posY=vpos(y, &(polygons[p].y[1]), npt);
    if(posY==-1) continue;
    if(posY!=posX) continue;
    printf("%s : polygons %d (%d) #overlapped %d \n", __func__, p, polygons[p].npt, posX);

/*------------------------------------------------------------------------------
    not safe, commented out */
//     work.duplicate(polygons[p], 0, (int) posX+1);
//     polygons.push_back(work);
//     polygons[p].npt=0;
//     solved++;
    }
  
  for(size_t p=0;p<polygons.size();p++) {
    if(polygons[p].npt<=1) {
      polygons.erase (polygons.begin()+p);
      p--;
      }
    }
    
  printf("%s : polygons size before=%d after=%d (solved=%d)\n", __func__,bkp,polygons.size(), solved); 
    
  return(solved);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_CheckOverlapped(vector<plg_t> & polygons, string rootname, bool repair, bool do_split, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  detect and fix overlapping polygons
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status, solved;
  bool join, split, overlapped;
  int count=0;
  
/*------------------------------------------------------------------------------
  check self-overlapping polygons */
//   count=0;
//   do {
//     solved=plg_CheckSelfOverlapped(polygons, repair, rootname, debug);
//     status=plg_save("tmp.plg", PLG_FORMAT_SCAN, polygons);
//     } while (solved!=0 and count < 10);
  
/*------------------------------------------------------------------------------
  check overlapping polygons, fix identicals */
  join=false; split=false; overlapped=true;
  solved=plg_CheckOverlapped(polygons, join, split, overlapped, rootname, debug);
  if(debug) status=plg_save("tmp.plg", PLG_FORMAT_SCAN, polygons);
  
/*------------------------------------------------------------------------------
  check overlapping polygons, fix connections */
  if(do_split) {
    join=false; split=true; overlapped=true;
    count=0;
    do {
      solved=plg_CheckOverlapped(polygons, join, split, overlapped, rootname, debug);
      if(debug) status=plg_save("tmp.plg", PLG_FORMAT_SCAN, polygons);
      count++;
      } while (solved!=0 and count < 10);
    }
    
  repair=true;
  status=plg_CheckClosure(polygons, repair, 1, debug);
  if(debug) status=plg_save("tmp.plg", PLG_FORMAT_SCAN, polygons);
  
/*------------------------------------------------------------------------------
  check overlapping polygons, fix concatable */
  join=true; split=true; overlapped=true;
  count=0;
  do {
    solved=plg_CheckOverlapped(polygons, join, split, overlapped, rootname, debug);
    if(debug) status=plg_save("tmp.plg", PLG_FORMAT_SCAN, polygons);
    count++;
    } while (solved!=0 and count < 20);
  
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
    "  merge a polygon set with an existing shoreline inside a polygon-defined extent \n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  Show this help and exit.\n"
    "  --debug  \n"
    "  --frame  followed by frame. WILL BE OVERWRITTEN TO [-10.0:0.0;42.5:50] !\n"
    "  -r  followed by the root name\n"
    "  -o  followed by the path of the output file. Default from the root name.\n"
    "  -O  followed by the output format\n"
    "  -p  followed by a path to a limiting polygone\n"
    "  -s  followed by the path to the shorelines. Default is /home/data/shorelines/gshhs-2.0/gshhs_f.cst\n"
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
  int n,status;
  string keyword,PolygonFile="", FrameString="", rootname="", zone="", output;
  string outFormatName="";
  string ShorelinesFile;
  int inFormat,outFormat;
  vector<plg_t> shorelines, tmp;
  frame_t frame;
  vector<plg_t> limits;
  point2D_t point;
  bool debug=false, repair;
  int size[2];
  vector<int> list[2];
  double rmax=500.0;
  bool do_join=true, do_overlapped=true, do_split=false;
  
  fct_echo( argc, argv);
  
  ShorelinesFile="";
  
  n=1;
  while (n < argc) {
    keyword=argv[n];
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

//         case 'p' :
//           PolygonFile=argv[n+1];
//           n++;
//           n++;
//           break;
        
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
        
//         case 'z' :
//           zone=argv[n+1];
//           n++;
//           n++;
//           break;
        
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
          else if(strncmp("--debug",keyword)==0){
            debug=true;
            n++;
            }
          else if(strncmp("--connections",keyword)==0){
            do_split=true;
            n++;
            }
//           else if(strncmp("--rmax",keyword)==0){
//             sscanf(argv[n+1],"%lf",&rmax);
//             n++;
//             n++;
//             }
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
  
  if(outFormatName!="") {
    outFormat=plg_format_from_name(outFormatName);
    if(outFormat==PLG_FORMAT_UNKNOWN){
      plg_print_formats();
      exit(-1);
      }
    }
  else {
    outFormat=PLG_FORMAT_SCAN;
    }
    
  if(output=="") output=rootname+"-doctored.plg";

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  detect format and open shoreline polygons
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(ShorelinesFile=="") {
    TRAP_ERR_EXIT(-1,"no shorelines file given, abort...\n");
    }
   
  inFormat=plg_find_format(ShorelinesFile);
  status=plg_load(ShorelinesFile, inFormat, shorelines);
  
  printf("%s : polygons size=%d format=%d\n", ShorelinesFile.c_str(), shorelines.size(), inFormat); 
    
//   if(zone!="") {
//     status=plg_SetZone(zone.c_str(), frame, rootname, point);
//     limits.push_back(plg_t(frame, PLG_INIT_SEPARATE));
//     status=plg_cartesian((projPJ) 0, limits);
//     }
// 
// //   FrameString= "[-10.0:0.0;42.5:50]";
//   if(FrameString!="") {
//     status=plg_DecodeFrame(FrameString, frame);
//     limits.push_back(plg_t(frame, PLG_INIT_SEPARATE));
//     status=plg_cartesian((projPJ) 0, limits);
//     }
// 
//   if(PolygonFile!="") {
//     inFormat=plg_find_format(PolygonFile);
//     status=plg_load(PolygonFile, inFormat, limits);
//     if(status!=0) TRAP_ERR_EXIT(status,"plg_load(\""+PolygonFile+"\",%d,,) error",inFormat);
//     }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  remove empty or one-point polygons
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  status=plg_CheckDuplicated(shorelines, true, rootname, false);

  for(size_t p=0;p<shorelines.size();p++) {
    if(shorelines[p].npt<=1) {
      shorelines.erase (shorelines.begin()+p);
      p--;
      }
    }
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  extract shorelines dataset in limits
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/



/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  check closure, perform concatenation when needed (repair=true)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  repair=true;
  status=plg_CheckOverlapped(shorelines, rootname, repair, do_split, debug);

//   status=plg_CheckOverlapped(shorelines, true, rootname, debug);
// 
//   status=plg_CheckClosure(shorelines, repair, 1, debug);
    
  printf("-----------------------> save %s\n",output.c_str());
  status=plg_save(output.c_str(), outFormat, shorelines);
  
  __OUT_BASE_LINE__(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(0);
error:
  __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}
