
/*******************************************************************************

  T-UGO tools, 2006-2018

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Assemble 2 meshes.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main <a href=#func-members>function</a> for how this works
and the print_help <a href=#func-members>function</a> for how to use this.
*/
/*----------------------------------------------------------------------------*/


#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "geo.h"
#include "fe.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "functions.h"


class split_t {
private :
public :
  mesh_t meshes[2];
  vector<plg_t> limits[2], opened[2];
  int *flags[2];
  bool *frontier[2];
  
  split_t() {
    flags[0]=0;
    flags[1]=0;
    frontier[0]=0;
    frontier[1]=0;
    }
};

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_concat(mesh_t meshes[2], mesh_t & final, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  concat 2 meshes
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status=0;
  string filename, rootname;
  size_t Voffset, Eoffset, Toffset;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  concat meshes
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  final.nvtxs=meshes[0].nvtxs+meshes[1].nvtxs;
  
  final.vertices=new vertex_t[final.nvtxs];
  
  final.ntriangles=meshes[0].ntriangles+meshes[1].ntriangles;
  
  final.triangles=new triangle_t[final.ntriangles];
  
  final.nedges=meshes[0].nedges+meshes[1].nedges;
  
  final.edges=new edge_t[final.nedges];
  
  Voffset=meshes[0].nvtxs;
  Eoffset=meshes[0].nedges;
  Toffset=meshes[0].ntriangles;

  for(size_t n=0; n<meshes[0].nvtxs; n++) {
    final.vertices[n]=meshes[0].vertices[n];
    }
   for(size_t n=0; n<meshes[1].nvtxs; n++) {
    final.vertices[n+Voffset]=meshes[1].vertices[n];
    for(int k=0; k<final.vertices[n+Voffset].nngh;k++) {
      final.vertices[n+Voffset].ngh[k]+=Voffset;
      }
    }
  
  for(size_t n=0; n<meshes[0].ntriangles; n++) {
    final.triangles[n]=meshes[0].triangles[n];
    }
   for(size_t n=0; n<meshes[1].ntriangles; n++) {
    final.triangles[n+Toffset]=meshes[1].triangles[n];
    for(int k=0; k<3;k++) {
      final.triangles[n+Toffset].vertex[k]+=Voffset;
      }
    }

  for(size_t n=0; n<meshes[0].nedges; n++) {
    final.edges[n]=meshes[0].edges[n];
    }
   for(size_t n=0; n<meshes[1].nedges; n++) {
    final.edges[n+Eoffset]=meshes[1].edges[n];
    for(int k=0; k<2;k++) {
      final.edges[n+Eoffset].extremity[k]+=Voffset;
      }
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_concat(mesh_t & m1, mesh_t & m2, mesh_t & final, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  concat 2 meshes
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
  string filename, rootname;
  
  mesh_t meshes[2];

  meshes[0]=m1;
  meshes[1]=m2;
  
  status=fe_concat(meshes, final, debug);
  
  
  return(status);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_merge02(mesh_t meshes[2], double dmax, char *dum, bool limited, mesh_t & final, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  concat 2 overlapping meshes
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
  bool exact=true;
  vector<plg_t> polygons[2];

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  get mesh limits as polygons
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_limits2poly(meshes[0],polygons[0],0,debug);
  status=fe_limits2poly(meshes[1],polygons[1],0,debug);

//   if(debug) status=plg_save("mesh-assembly-0.plg", PLG_FORMAT_SCAN, polygons[0]);
  
  status=plg_merge(polygons[0],polygons[1]);
  
  status=plg_save("mesh-assembly-merged.plg", PLG_FORMAT_SCAN, polygons[0]);
   
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  cut mesh 1 using mesh 0 polygons, creating internal and external sub-domains
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_MeshCut(meshes[1], polygons[0], exact, debug);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  merge mesh 0 and mesh 1 external sub-domains
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  final=fe_merge(meshes, dmax, (char *) 0, limited, debug);

  return(status);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_split(mesh_t & mesh, int *selected, int size, split_t & split , bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  selected: edge flag for interior/exterior splitting (1 is interior)
  
  size    : number of edges needed to fulfill criterion. if size > 0, then option=0
            else option=1
  
  option 0: select triangle as interior if at least "size" edges are eligible
  
    when size is 3, only triangles strictly included in the selection are taken as
    interior. when size is 1, boarder triangles are taken as interior
  
  option 1: select triangle as exterior if at least "size" edges are not eligible
  
    when size is 3, only triangles strictly excluded from the selection are taken as
    exterior. when size is 1, boarder triangles are taken as exterior
  
    non-symmetric condition
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int   l,m,m1,m2,i,n1,n2,n,status;
  int   count=0;
  int   nndes;
  int   ninternal,nexternal;
  bool  *keep=new bool[mesh.ntriangles];
  int   *used=NULL;
  mesh_t *work=NULL;
  string filename;
  int option;
  bool extras=true;

  work=split.meshes;

  printf ("number of nodes    (original mesh): %6d\n",mesh.nvtxs);
  printf ("number of elements (original mesh): %6d\n",mesh.ntriangles);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  screen edges

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ninternal=0;
  nexternal=0;

  if(size>0) option=0;
  else {
    option =1;
    size=-size;
    }

  for (m=0;m<mesh.ntriangles;m++) {
    keep[m]=false;
    count=0;
    switch (option) {
      case 0:
        for(i=0;i<3;i++) {
          n=mesh.triangles[m].edges[i];
          if(selected[n]==1) {
            count++;
            }
          }
        keep[m]=(count>=size);
        break;
      case 1:
        for(i=0;i<3;i++) {
          n=mesh.triangles[m].edges[i];
          if(selected[n]==0) {
            count++;
            }
          }
        keep[m]=not (count>=size);
        break;
      }
    switch(keep[m]) {
      case false:
        nexternal++;
        break;
      case true:
        ninternal++;
        break;
      }
    }

  work[0].triangles=new triangle_t [nexternal];
  work[1].triangles=new triangle_t [ninternal];

  work[0].ntriangles=nexternal;
  work[1].ntriangles=ninternal;

  work[0].type=SPHERICAL;
  work[1].type=SPHERICAL;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  screen elements

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ninternal=0;
  nexternal=0;
  for (m=0;m<mesh.ntriangles;m++) {
    switch(keep[m]) {
      case false:
        work[0].triangles[nexternal].ancestor=m;
        nexternal++;
        break;
      case true:
        work[1].triangles[ninternal].ancestor=m;
        ninternal++;
        break;
      }
    }

  used=new int[mesh.nvtxs];
  
  for(l=0;l<2;l++) {
    nndes=0;
    for (n=0;n<mesh.nvtxs;n++) used[n]=-1;
    for (m=0;m<work[l].ntriangles;m++) {
      for(i=0;i<3;i++) {
        n=mesh.triangles[work[l].triangles[m].ancestor].vertex[i];
        if(used[n]==-1) {
          used[n]=nndes;
          nndes++;
          }
        }
      }
    work[l].nvtxs = nndes;
    work[l].vertices=new vertex_t[work[l].nvtxs];
    for (n=0; n<work[l].nvtxs; n++) work[l].vertices[n].null_value();
    for (n=0;n<mesh.nvtxs;n++) {
       if(used[n]!=-1) {
         work[l].vertices[used[n]].lon=mesh.vertices[n].lon;
         work[l].vertices[used[n]].lat=mesh.vertices[n].lat;
         work[l].vertices[used[n]].h  =mesh.vertices[n].h;
         work[l].vertices[used[n]].code=0;
         work[l].vertices[used[n]].ancestor=n;
         }
       }
    for (m=0;m<work[l].ntriangles;m++) {
      for(i=0;i<3;i++) {
        n=mesh.triangles[work[l].triangles[m].ancestor].vertex[i];
        if(used[n]!=-1) {
          work[l].triangles[m].vertex[i]=used[n];
          }
        else {
          printf("error...\n");
          }
        }
      }
    printf ("number of nodes    (splitted mesh %d): %6d\n",l,work[l].nvtxs);
    printf ("number of elements (splitted mesh %d): %6d\n",l,work[l].ntriangles);

    status=fe_e2n(&(work[l]));
    if(debug) status=fe_savemesh("tmp-00.nei",MESH_FILE_FORMAT_TRIGRID,work[l]);

    status=fe_geometry(&(work[l]));
    status=fe_edgetable(&(work[l]),0,0);
    status=fe_vertex_crosstables02(&(work[l]));

    int RecycleCodes=0, SetLimits=1, verbose=0;
    status=fe_codetable(work[l], RecycleCodes, SetLimits, verbose);
    if(status!=0) {
      printf ("code table failed for mesh %d\n",l);
      return(0);
      }
    
/*------------------------------------------------------------------------------
    build edge's ancestor, use elements information*/
    for (m1=0;m1<work[l].ntriangles;m1++) {
      m2=work[l].triangles[m1].ancestor;
      for(i=0;i<3;i++) {
        n1=work[l].triangles[m1].edges[i];
        n2=mesh.triangles[m2].edges[i];
        work[l].edges[n1].ancestor=n2;
        }
      }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      
  frontier: edge flag, true if located along the interior/exterior splitting limit
            
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    if(extras) {
      split.frontier[l]=new bool[work[l].nedges];
      for (n=0;n<work[l].nedges;n++) {
        int nn=work[l].edges[n].ancestor;
        if((mesh.edges[nn].code==MESH_INTERIOR_EDGE) && (work[l].edges[n].code==1)) {
/*------------------------------------------------------------------------------
          internal/external limit */
          split.frontier[l][n]=true;
          }
        else {
          split.frontier[l][n]=false;
          }
        }
      }
    }

  if(debug) {
    filename="debug-external.nei";
    status=fe_savemesh(filename.c_str(),MESH_FILE_FORMAT_TRIGRID,work[0]);
    filename="debug-internal.nei";
    status=fe_savemesh(filename.c_str(),MESH_FILE_FORMAT_TRIGRID,work[1]);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create boundary polygon file for further mesh assembly
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(l=0;l<2;l++) {
    char *flag=new char[work[l].nedges];
    
    count=0;
    for(n=0;n<work[l].nedges;n++) {
      if(work[l].edges[n].code==0) flag[n]='I';
      else flag[n]='T';
      }
    
    for (n=0;n<work[l].nedges;n++) {
      n2=work[l].edges[n].ancestor;
      if((mesh.edges[n2].code==MESH_INTERIOR_EDGE) && (work[l].edges[n].code>=1)) {
/*------------------------------------------------------------------------------
        internal/external limit */
        flag[n]='M';
        work[l].edges[n].code=0;
        count++;
        }
      }
    
/*------------------------------------------------------------------------------
    export mesh limits as polygons, passing open/rigid flag */
    vector<plg_t> polygons, splitted;
    status=fe_limits2poly(work[l], polygons, flag, false);
    
/*------------------------------------------------------------------------------
    split polygons to separate open/rigid segments */
    splitted=plg_split(polygons, false);
    for(int s=0;s<splitted.size();s++) {
      if(splitted[s].flag[0]=='M') split.opened[l].push_back(splitted[s]);
      }
    
/*------------------------------------------------------------------------------
    open limits can be made of 2 joining polygons */
    if(split.opened[l].size()>1) {
      status=plg_concat(split.opened[l][0], split.opened[l][1]);
      }
    }
  
  delete[] used;

//   printf("split completed\n");

  return(0);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int find_companions(plg_t & p, plg_t & q, vector<int>* & connected, bool debug)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//  
//   find p points that attract q points
//   
//   for a given p point, connected return best q point if any
//  
// cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
// {
//   int n, status;
// 
//   vector<int>    *neighbors=new vector<int>   [p.npt];
//   vector<double> *distances=new vector<double>[p.npt];
//   
//   for(int k=0; k<q.npt; k++) {
//     double d;
//     int mode=PLG_SPHERICAL, n;
//     n=plg_find_point(p, mode, q.t[k], q.p[k], &d);
// /*------------------------------------------------------------------------------
//     n-th p's point is the closest point of k-th q's point */
//     neighbors[n].push_back(k);
//     distances[n].push_back(d);
//     }
//   
//   connected=new vector<int>[p.npt];
//   
//   for(int k=0; k<p.npt; k++) {
//     double dmin;
//     printf("%3d : (#neighbours=%2d)",k,neighbors[k].size());
//     if(neighbors[k].size()==0) {
//       printf("\n");
//       continue;
//       }
//     if(neighbors[k].size()==1) {
//       connected[k].push_back(neighbors[k][0]);
//       dmin=distances[k][0];
//       }
//     if(neighbors[k].size()>1) {
//       int pos=minpos(distances[k], &dmin);
//       connected[k].push_back(neighbors[k][pos]);
//       }
//     for(int i=0; i<connected[k].size(); i++) printf(" %3d (%9.3lf)",connected[k][i], dmin);
//     printf("\n");
//     }
//   
//   for(int k=0; k<p.npt; k++) neighbors[k].clear();
//   delete[] neighbors;
//   for(int k=0; k<p.npt; k++) distances[k].clear();
//   delete[] distances;
//   
//   return(0);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int fe_merge03(mesh_t meshes[2], string PolygonsFile, double dmax, char *dum, bool limited, mesh_t & final, bool debug)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//  
//   concat 2 overlapping meshes
//  
// cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
// {
//   int n, status;
//   vector<plg_t> polygons;
//   bool exact=false;
//   mesh_t work[2], *splitted;
//   split_t split[2];
//   int *selected, nselected;
//   int size=3;
//   string filename, rootname;
//   
//   if(PolygonsFile=="") TRAP_ERR_EXIT(-1,"%s : polygons file not set\n", __func__);
//   status=plg_load(PolygonsFile, PLG_FORMAT_UNKNOWN, polygons);
//   
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//  
//   split meshes
//   
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// 
//   for(int k=0;k<2;k++) {
//     exact=true;
//     if(exact) {
//       status=fe_presplit(meshes[k], polygons, 0.5, "split-adjusted");
//       if(status!=0) return(-1);
//       }
//     selected=new int[meshes[k].nedges];
//     for (n=0;n<meshes[k].nedges;n++) {
//       selected[n]=1;
//       }
//     int option=0;
//     nselected=fe_selectedges_01(meshes[k], polygons, selected, option, false);
//     if(nselected==0) return(0);
//     status=fe_split(meshes[k], selected, size, split[k], debug);
//     delete[] selected;
//     }
//     
//   plg_destroy_entries(polygons);
//   polygons.clear();
//   
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//  
//   save interior and exterior meshes
//   
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// 
//   if(rootname=="") rootname="assembly";
//   
//   if(debug) {
//     filename=rootname+"-external.nei";
//     status=fe_savemesh(filename.c_str(),MESH_FILE_FORMAT_TRIGRID,split[0].meshes[0]);
//   
//     filename=rootname+"-external.plg";
//     status=fe_SaveLimits(filename.c_str(),split[0].meshes[0]);
//   
//     filename=rootname+"-internal.nei";
//     status=fe_savemesh(filename.c_str(),MESH_FILE_FORMAT_TRIGRID,split[1].meshes[1]);
//   
//     filename=rootname+"-internal.plg";
//     status=fe_SaveLimits(filename.c_str(),split[1].meshes[1]);
//     }
//     
//   filename=rootname+"-external-opened.plg";
//   status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, split[0].opened[0]);
//   
//   filename=rootname+"-internal-opened.plg";
//   status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, split[1].opened[1]);
//   
//   for(int k=0; k<2; k++) {
//     printf("connection polygons, mesh %3d : (%2d %2d)\n",k,split[k].opened[0].size(),split[k].opened[1].size());
//     }
//   
//   plg_t tmp[2];
//   
//   tmp[0].duplicate(split[0].opened[0][0]);
//   polygons.push_back(tmp[0]);
//   
//   tmp[1].duplicate(split[1].opened[1][0]);
//   polygons.push_back(tmp[1]);
//   
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//  
//   concat mesh 0 external mesh and mesh 1 internal sub-meshes
//   
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// 
//   for(int k=0; k<polygons[0].npt; k++) {
//     double t,p;
//     int n;
//     t=polygons[0].t[k];
//     p=polygons[0].p[k];
//     n=fe_nearest_boundaryvertex(split[0].meshes[0], t, p, 0, 0);
//     split[0].meshes[0].vertices[n].code=-1;
//     }
//     
//   for(int k=0; k<polygons[0].npt; k++) {
//     double t,p;
//     int n;
//     t=polygons[1].t[k];
//     p=polygons[1].p[k];
//     n=fe_nearest_boundaryvertex(split[1].meshes[1], t, p, 0, 0);
//     split[1].meshes[1].vertices[n].code=-2;
//     }
//     
//   status=fe_concat(split[0].meshes[0], split[1].meshes[1], final, debug);
//   filename=rootname+"-raw.nei";
//   status=fe_savemesh(filename.c_str(),MESH_FILE_FORMAT_TRIGRID,final);
//   
//   vector<int> list[2];
//   
//   for(int k=0; k<polygons[0].npt; k++) {
//     double t,p;
//     int n;
//     t=polygons[0].t[k];
//     p=polygons[0].p[k];
//     n=fe_nearest_vertex(final, t, p, -1);
//     list[0].push_back(n);
//     }
//     
//   for(int k=0; k<polygons[1].npt; k++) {
//     double t,p;
//     int n;
//     t=polygons[1].t[k];
//     p=polygons[1].p[k];
//     n=fe_nearest_vertex(final, t, p, -2);
//     list[1].push_back(n);
//     }
//     
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//  
//   adjust mesh limits to prepare merging
//   
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// 
//   
// #if 0
//   int npt=max(polygons[0].npt,polygons[1].npt);
//   
//   vector<int> *paires=new vector<int>[npt];
//   vector<int> vertex[2];
//  
//   paires[0].push_back(0);
//   paires[0].push_back(0);
//   
//   n=fe_nearest_boundaryvertex(split[0].meshes[0], polygons[0].t[0], polygons[0].p[0], 0, 0);
//   vertex[0].push_back(n);
//   n=fe_nearest_boundaryvertex(split[1].meshes[1], polygons[1].t[0], polygons[1].p[0], 0, 0);
//   vertex[1].push_back(n);
//   
//   status=plg_delete_point(&polygons[0],0);
//   status=plg_delete_point(&polygons[1],0);
//   
//   for(int k=1; k<npt-1; k++) {
//     double t,p,d;
//     int mode=PLG_SPHERICAL, n1, n2;
//     t=polygons[1].t[0];
//     p=polygons[1].p[0];
//     n1=plg_find_point(polygons[0], mode, t, p, &d);
//     t=polygons[0].t[0];
//     p=polygons[0].p[0];
//     n2=plg_find_point(polygons[1], mode, t, p, &d);
// /*------------------------------------------------------------------------------
//      */
//     if(n1==0 and n2==0) {
//       paires[k].push_back(k);
//       paires[k].push_back(k);
//       
//       n=fe_nearest_boundaryvertex(split[0].meshes[0], polygons[0].t[0], polygons[0].p[0], 0, 0);
//       vertex[0].push_back(n);
//       n=fe_nearest_boundaryvertex(split[1].meshes[1], polygons[1].t[0], polygons[1].p[0], 0, 0);
//       vertex[1].push_back(n);
//       
//       status=plg_delete_point(&polygons[0],0);
//       status=plg_delete_point(&polygons[1],0);
//       continue;
//       }
//     if(n1>n2) {
//       plg_point_t point=plg_point_t(polygons[0],0);
//       status=plg_insert_point(polygons[1],0,point);
//       
//       n1=vertex[1][vertex[1].size()-1];
//       n2=fe_nearest_boundaryvertex(split[1].meshes[1], polygons[1].t[0], polygons[1].p[0], 0, 0);
//       
//       k--;
//       continue;
//       }
//     if(n2>n1) {
//       plg_point_t point=plg_point_t(polygons[1],0);
//       status=plg_insert_point(polygons[0],0,point);
//       k--;
//       continue;
//       }
//     }
//     
//   paires[0].push_back(npt-1);
//   paires[0].push_back(npt-1);
//   
//   status=plg_delete_point(&polygons[0],0);
//   status=plg_delete_point(&polygons[1],0);
// #endif
//     
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//  
//   loop on polygons[0] points, and get their best polygons[1] points attractors
//   
//   may find several polygons[0] points for one polygons[1] attractors, takes the
//   closest one
//   
//   some polygons[1] points may attract none and will be left alone
//   
//   non-symetric algorithm (sic)
//   
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// 
//   vector<int> *connected[2];
//   status=find_companions(polygons[1], polygons[0], connected[0], debug);
// 
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//  
//   merge companion vertices
//   
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
//   
//   for(int k=0; k<polygons[1].npt; k++) {
//     if(connected[0][k].size()!=1) continue;
//     double tt1,t1,p1;
//     double tt2,t2,p2;
//     double t,p;
//     int l,n1,n2,remove=0;
//     l=connected[0][k][0];
//     t1=polygons[0].t[l];
//     p1=polygons[0].p[l];
//     n1=fe_nearest_vertex(final, t1, p1, -1);
//     t2=polygons[1].t[k];
//     p2=polygons[1].p[k];
//     n2=fe_nearest_vertex(final, t2, p2, -2);
//     printf("connect %6d with %6d", n1,n2);
//     printf("\n");
// /*------------------------------------------------------------------------------
//     move vertices to their mid-position                                       */
//     tt1=degree_recale(t1, t2);
//     t=0.5*(tt1+t2);
//     p=0.5*(p1+p2);
//     status=fe_displacevertex(final, n1, t, p);
//     status=fe_displacevertex(final, n2, t, p);
// /*------------------------------------------------------------------------------
//     and merge vertices                                                        */
//     status=fe_mergevertex(&final, n1, n2, remove);
// //     final.vertices[n2].code=0;
// /*------------------------------------------------------------------------------
//     move polygons points                                                      */
//     tt2=degree_recale(t2, t1);
//     t=0.5*(t1+tt2);
//     p=0.5*(p1+p2);
//     polygons[0].t[l]=t;
//     polygons[0].p[l]=p;
//     tt1=degree_recale(t1, t2);
//     t=0.5*(tt1+t2);
//     p=0.5*(p1+p2);
//     polygons[1].t[k]=t;
//     polygons[1].p[k]=p;
//     }
//     
//   status=fe_cleanvertices(&final, false);
//   status=fe_list(&final);
//   status=fe_edgetable(&final,0,0);
// 
//   filename=rootname+"-raw-1.nei";
//   status=fe_savemesh(filename.c_str(),MESH_FILE_FORMAT_TRIGRID,final);
//   
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//  
//   create interlaced vertices
//   
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
//   
//   status=fe_Vertex_EdgesXtables(&final);
//   
//   for(int k=1; k<polygons[1].npt-1; k++) {
//     double dmin;
//     if(connected[0][k].size()!=0) {
//       continue;
//       }
//     if(connected[0][k-1].size()==0 or connected[0][k+1].size()==0) {
//       continue;
//       }
//     double tt1,t1,p1;
//     double tt2,t2,p2;
//     double t,p;
//     int l,m,n,n1,n2,v,remove=0;
//     l=connected[0][k-1][0];
//     t1=polygons[0].t[l];
//     p1=polygons[0].p[l];
//     n1=fe_nearest_vertex(final, t1, p1, -1);
//     l=connected[0][k+1][0];
//     t2=polygons[0].t[l];
//     p2=polygons[0].p[l];
//     n2=fe_nearest_vertex(final, t2, p2, -1);
//     n=fe_isedge(final, n1,n2);
//     if(n==-1) continue;
//     v=fe_AddVertexOnEdge(&final, n);
//     t2=polygons[1].t[k];
//     p2=polygons[1].p[k];
//     n2=fe_nearest_vertex(final, t2, p2, -2);
// /*------------------------------------------------------------------------------
//     merge vertices */
//     printf("connect %6d with %6d", v,n2);
//     status=fe_mergevertex(&final, v, n2, remove);
// //     l=connected[0][k-1][0];
// //     plg_point_t point=plg_point_t(polygons[1],0);
// //     status=plg_insert_point(polygons[0],0,point);
//     status=fe_Vertex_EdgesXtables(&final);
//     }
// 
//   status=fe_cleanvertices(&final, false);
//   status=fe_list(&final);
//   status=fe_edgetable(&final,0,0);
// 
//   filename=rootname+"-raw-2.nei";
//   status=fe_savemesh(filename.c_str(),MESH_FILE_FORMAT_TRIGRID,final);
//   
//   return(status);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int find_companions(mesh_t & mesh, vector<int> p,  vector<int> q, int code, vector<int>* & connected, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  find p points that attract q points
  
  for a given p point, connected return best q point if any
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  vector<int>    *neighbors=new vector<int>   [p.size()];
  vector<double> *distances=new vector<double>[p.size()];
  
  for(int k=0; k<q.size(); k++) {
    double tt,pp,d;
    int n;
    tt=mesh.vertices[q[k]].lon;
    pp=mesh.vertices[q[k]].lat;
    n=fe_nearest_vertex(mesh, tt, pp, code, d);
/*------------------------------------------------------------------------------
    n-th p's point is the closest point of k-th q's point */
    int pos=vpos(n,p);
    neighbors[pos].push_back(q[k]);
    distances[pos].push_back(d);
    }
  
  connected=new vector<int>[p.size()];
  
  for(int k=0; k<p.size(); k++) {
    double dmin;
    printf("%3d %6d: (#neighbours=%2d)",k,p[k],neighbors[k].size());
    if(neighbors[k].size()==0) {
      printf("\n");
      continue;
      }
    if(neighbors[k].size()==1) {
      connected[k].push_back(neighbors[k][0]);
      dmin=distances[k][0];
      }
    if(neighbors[k].size()>1) {
      int pos=minpos(distances[k], &dmin);
      connected[k].push_back(neighbors[k][pos]);
      }
    for(int i=0; i<connected[k].size(); i++) printf(" %3d (%9.3lf)",connected[k][i], dmin);
    printf("\n");
    }
  
  for(int k=0; k<p.size(); k++) neighbors[k].clear();
  delete[] neighbors;
  for(int k=0; k<p.size(); k++) distances[k].clear();
  delete[] distances;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_merge03(mesh_t meshes[2], string PolygonsFile, double dmax, char *dum, bool limited, mesh_t & final, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  concat 2 overlapping meshes
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int n, status;
  vector<plg_t> polygons;
  bool exact=false;
  split_t split[2];
  int *selected, nselected;
  int size=3;
  string filename, rootname;
  
  if(PolygonsFile=="") TRAP_ERR_EXIT(-1,"%s : polygons file not set\n", __func__);
  status=plg_load(PolygonsFile, PLG_FORMAT_UNKNOWN, polygons);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  split meshes
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(int k=0;k<2;k++) {
    exact=true;
    if(exact) {
      status=fe_presplit(meshes[k], polygons, 0.5, "split-adjusted");
      if(status!=0) return(-1);
      }
    status=fe_savemesh("adjusted.nei",MESH_FILE_FORMAT_TRIGRID,meshes[0]);
    selected=new int[meshes[k].nedges];
    for (n=0;n<meshes[k].nedges;n++) {
      selected[n]=1;
      }
    int option=0;
    nselected=fe_selectedges_01(meshes[k], polygons, selected, option, false);
    if(nselected==0) return(0);
//     debug=true;
    status=fe_split(meshes[k], selected, size, split[k], debug);
    delete[] selected;
    }
  
  plg_destroy_entries(polygons);
  polygons.clear();
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  save interior and exterior meshes
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(rootname=="") rootname="assembly";
  
  if(debug) {
    filename=rootname+"-external.nei";
    status=fe_savemesh(filename.c_str(),MESH_FILE_FORMAT_TRIGRID,split[0].meshes[0]);
  
    filename=rootname+"-external.plg";
    status=fe_SaveLimits(filename.c_str(),split[0].meshes[0]);
  
    filename=rootname+"-internal.nei";
    status=fe_savemesh(filename.c_str(),MESH_FILE_FORMAT_TRIGRID,split[1].meshes[1]);
  
    filename=rootname+"-internal.plg";
    status=fe_SaveLimits(filename.c_str(),split[1].meshes[1]);
    }
  
  filename=rootname+"-external-opened.plg";
  status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, split[0].opened[0]);
  
  filename=rootname+"-internal-opened.plg";
  status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, split[1].opened[1]);
  
  for(int k=0; k<2; k++) {
    printf("connection polygons, mesh %3d : (%2d %2d)\n",k,split[k].opened[0].size(),split[k].opened[1].size());
    }
  
  plg_t tmp[2];
  
  tmp[0].duplicate(split[0].opened[0][0]);
  polygons.push_back(tmp[0]);
  
  tmp[1].duplicate(split[1].opened[1][0]);
  polygons.push_back(tmp[1]);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  concat mesh 0 external mesh and mesh 1 internal sub-meshes
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(int k=0; k<polygons[0].npt; k++) {
    double t,p;
    int n;
    t=polygons[0].t[k];
    p=polygons[0].p[k];
    n=fe_nearest_boundaryvertex(split[0].meshes[0], t, p, 0, 0);
    split[0].meshes[0].vertices[n].code=-1;
    }
  
  for(int k=0; k<polygons[0].npt; k++) {
    double t,p;
    int n;
    t=polygons[1].t[k];
    p=polygons[1].p[k];
    n=fe_nearest_boundaryvertex(split[1].meshes[1], t, p, 0, 0);
    split[1].meshes[1].vertices[n].code=-2;
    }
  
  status=fe_concat(split[0].meshes[0], split[1].meshes[1], final, debug);
  filename=rootname+"-raw.nei";
  status=fe_savemesh(filename.c_str(),MESH_FILE_FORMAT_TRIGRID,final);
  
  vector<int> list[2];
  
  for(int k=0; k<polygons[0].npt; k++) {
    double t,p,d;
    int n;
    t=polygons[0].t[k];
    p=polygons[0].p[k];
    n=fe_nearest_vertex(final, t, p, -1, d);
    list[0].push_back(n);
    }
  
  for(int k=0; k<polygons[1].npt; k++) {
    double t,p,d;
    int n;
    t=polygons[1].t[k];
    p=polygons[1].p[k];
    n=fe_nearest_vertex(final, t, p, -2, d);
    list[1].push_back(n);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  loop on list[0] vertices, and get their best list[1] vertices attractors
  
  may find several list[0] vertices for one list[1] attractors, takes the
  closest one
  
  some list[1] vertices may attract none and will be left alone
  
  connected contails list[0] references
  
  non-symetric algorithm (sic)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  vector<int> *attracted[2];
  status=find_companions(final, list[1], list[0], -2, attracted[0], debug);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  merge companion vertices
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  vector<paire_t> mergers;
  
  for(int k=0; k<list[1].size(); k++) {
    if(attracted[0][k].size()!=1) continue;
    double tt1,t1,p1;
    double t2,p2;
    double t,p;
    int n1,n2;
    n1=list[1][k];
    n2=attracted[0][k][0];
    t1=final.vertices[n1].lon;
    p1=final.vertices[n1].lat;
    t2=final.vertices[n2].lon;
    p2=final.vertices[n2].lat;
/*------------------------------------------------------------------------------
    move vertices to their mid-position                                       */
    tt1=degree_recale(t1, t2);
    t=0.5*(tt1+t2);
    p=0.5*(p1+p2);
    status=fe_displacevertex(final, n1, t, p);
    status=fe_displacevertex(final, n2, t, p);
/*------------------------------------------------------------------------------
    register mergers                                                          */
    printf("connect %6d with %6d\n", n1,n2);
    paire_t paire(n1,n2);
    mergers.push_back(paire);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create interlaced vertices
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  int* companions[2];
  
  companions[0]=new int[list[0].size()];
  for(int k=0; k<list[0].size(); k++) companions[0][k]=-1;
  
  companions[1]=new int[list[1].size()];
  for(int k=0; k<list[1].size(); k++) companions[1][k]=-1;
      
  for(int k=0; k<list[1].size(); k++) {
    if(attracted[0][k].size()!=1) continue;
    int n=attracted[0][k][0];
    companions[1][k]=n;
    int pos=vpos(n, list[0]);
    companions[0][pos]=list[1][k];
    }
  
  printf("mesh #0 : %d unsolved vertices\n",occurence(-1, companions[0], list[0].size()));
  printf("mesh #1 : %d unsolved vertices\n",occurence(-1, companions[1], list[1].size()));
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create interlaced vertices
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  status=fe_Vertex_EdgesXtables(&final);
  
  for(int k=1; k<list[1].size()-1; k++) {
    if(companions[1][k]!=-1) {
      continue;
      }
    int n1, n2;
    n1=companions[1][k-1];
    n2=companions[1][k+1];
    if(n1==-1 or n2==-1) {
      continue;
      }
    int n, nn, v;
    n=fe_isedge(final, n1, n2);
    if(n==-1) continue;
    v=fe_AddVertexOnEdge(&final, n);
    status=fe_list(&final);
    status=fe_edgetable(&final,0,0);
    nn=list[1][k];
/*------------------------------------------------------------------------------
    register mergers                                                          */
    printf("create %6d (between %6d %6d) and connect with %6d\n", v, n1, n2, nn);
    paire_t paire(v,nn);
    mergers.push_back(paire);
    status=fe_Vertex_EdgesXtables(&final);
    }

  for(int k=1; k<list[0].size()-1; k++) {
    if(companions[0][k]!=-1) {
      continue;
      }
    int n1, n2;
    n1=companions[0][k-1];
    n2=companions[0][k+1];
    if(n1==-1 or n2==-1) {
      continue;
      }
    int n, nn, v;
    n=fe_isedge(final, n1, n2);
    if(n==-1) continue;
    v=fe_AddVertexOnEdge(&final, n);
    status=fe_list(&final);
    status=fe_edgetable(&final,0,0);
    nn=list[0][k];
/*------------------------------------------------------------------------------
    register mergers                                                          */
    printf("create %6d (between %6d %6d) and connect with %6d\n", v, n1, n2, nn);
    paire_t paire(v,nn);
    mergers.push_back(paire);
    status=fe_Vertex_EdgesXtables(&final);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  merge vertices
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  for(int k=0; k<mergers.size(); k++) {
    int remove=0;
    printf("merge %6d with %6d\n", mergers[k].value[0], mergers[k].value[1]);
    status=fe_mergevertex(&final, mergers[k].value[0], mergers[k].value[1], remove);
    for(int l=k+1; l<mergers.size(); l++) {
      if(mergers[l].value[0]==mergers[k].value[0]){
        mergers[l].value[0]=mergers[k].value[1];
        }
      }
    for(int l=k+1; l<mergers.size(); l++) {
      if(mergers[l].value[1]==mergers[k].value[0]){
        mergers[l].value[1]=mergers[k].value[1];
        }
      }
    }
  
  status=fe_cleanvertices(&final, false);
  status=fe_list(&final);
  status=fe_edgetable(&final,0,0);

  filename=rootname+"-raw-2.nei";
  status=fe_savemesh(filename.c_str(),MESH_FILE_FORMAT_TRIGRID,final);
  
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
    "  %s small.nei big.nei [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Assemble 2 meshes.\n"
    "For better performance, put the small one first\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help : Show this help and exit.\n"
    "  --debug : extra debugging output and save intermediate file merged.nei\n"
    "  --limited : limit merge to first boundary in first mesh\n"
    "  -d : followed by maximum merging distance in km. Default: 0.1km\n"
    "  -o : followed by path to assembled mesh. Default: mesh-assembly.nei\n"
    "  -p : followed by path to seperating polygon. Only for use with --overlapped\n"
    "  --raw : only contact meshes\n"
    "  --overlapped : merge exterior of first mesh with interior of second mesh, with respect to polygon given with -p\n"
    "  --no-reshape : do not reshape final mesh\n"
    );
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int k,n;
  char *keyword;
  char *meshfile[2]={NULL,NULL},*output=NULL;
  string PolygonsFile;
  mesh_t mesh[2],final;
  double dmax=0.1;
  int nmesh=0;
  bool limited=false, reshape=true, overlapped=false,raw=false;
  bool debug=false;

  fct_echo(argc,argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        if(strcmp(keyword,"--limited")==0) {
          limited=true;
          n++;
          break;
          }
        if(strcmp(keyword,"--debug")==0) {
          debug=true;
          n++;
          break;
          }
        if(strcmp(keyword,"--raw")==0) {
          raw=true;
          n++;
          break;
          }
        if(strcmp(keyword,"--overlapped")==0) {
          overlapped=true;
          n++;
          break;
          }
        if(strcmp(keyword,"--no-reshape")==0) {
          reshape=false;
          n++;
          continue;
          }
        if(strcmp(keyword,"--help")==0) {
          print_help(argv[0]);
          wexit(0);
          }
        switch (keyword[1]) {
        case 'h' :
          print_help(argv[0]);
          wexit(0);

        case 'd' :
          sscanf(argv[n+1],"%lf",&dmax);
          n++;
          n++;
          break;

        case 'p' :
          PolygonsFile=argv[n+1];
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          wexit(-1);
        }
        break;

      default:
        if(nmesh>=2){
          STDOUT_BASE_LINE("*** Only 2 meshes ***\n");
          print_help(argv[0]);
          wexit(-1);
          }
        
        meshfile[nmesh]= strdup(argv[n]);
        n++;
        nmesh++;
        break;
      }
    free(keyword);
    }

  if(nmesh!=2){
    STDOUT_BASE_LINE("*** Only %d mesh given when 2 are needed ***\n",nmesh);
    print_help(argv[0]);
    wexit(-1);
    }
  
  if(output==0) {
    output=strdup("mesh-assembly.nei");
    printf("no filename specified for output (-o [filename]); using default name (%s)...\n",output);
    }

  printf("#################################################################\n");
  printf("load meshes ant initialize related tables\n");
  for (k=0;k<2;k++) {
    if(meshfile[k] != NULL) {
      printf("load %d mesh : %s\n",k+1,meshfile[k]);
      status=fe_readmesh(meshfile[k],MESH_FILE_FORMAT_TRIGRID,&mesh[k]);
      if(status!=0) {
        printf("unable to read the original mesh in %s\n",meshfile[k]);
        TRAP_ERR_EXIT(status,"fe_readmesh() error\n");
        }
      status=fe_list(&mesh[k]);
      if(status!=0) {
        printf("unable to build the element list from the original mesh\n");
        TRAP_ERR_EXIT(status,"fe_list() error\n");
        }
      status= fe_edgetable(&mesh[k],0,0);
      int RecycleCodes=0/*, StopOn_EdgeError=1, StopOn_PinchError=0*/;
      int SetLimits=1;
//       status= fe_codetable2(&mesh[k], RecycleCodes, StopOn_EdgeError, StopOn_PinchError, 0);
      status=fe_codetable(mesh[k], RecycleCodes, SetLimits, 0);
      }
    else {
      printf("no %dth mesh filename specified; abort...\n",k);
      print_help(argv[0]);
      wexit(-1);
      }
    }

  if(raw) {
    printf("#################################################################\n");
    printf("concat meshes\n");
    status=fe_concat(mesh, final, debug);
    }
  else if(overlapped) {
    printf("#################################################################\n");
    printf("merge exterior of mesh %s with interior of mesh %s\n", meshfile[0],meshfile[1]);
    printf("interior/exterior defined with %s, maximum merging distance=%lf km\n",PolygonsFile.c_str(),dmax);
//     status=fe_merge02(mesh, dmax, (char *) 0, limited, final, debug);
    status=fe_merge03(mesh, PolygonsFile, dmax, (char *) 0, limited, final, debug);
    }
  else {
    printf("#################################################################\n");
    printf("merge meshes, maximum merging distance=%lf km\n",dmax);
    final=fe_merge(mesh, dmax, (char *) 0, limited, debug);
    }
  
  status= fe_savemesh("mesh-assembly-no-reshape.nei",MESH_FILE_FORMAT_TRIGRID, final);
  
  if(reshape) {
    printf("#################################################################\n");
    printf("reshape final mesh\n");
    status= fe_reshapeall(final,3);
    }
  
  printf("#################################################################\n");
  printf("save reshaped mesh : %s\n", output);
  status= fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID, final);

  STDOUT_BASE_LINE("end of mesh-assembly ... \n");
  exit(0);
}
