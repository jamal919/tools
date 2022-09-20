
/*******************************************************************************

  T-UGO tools, 2006-2015

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Refines a mesh.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main <a href=#func-members>function</a> for how this works
and the print_help <a href=#func-members>function</a> for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"
#include "constants.h"

#include "geo.h"
#include "fe.h"
#include "fe-proto.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "functions.h"
#include "statistic.h"

#include "zapper.h"     /*  rutin.h contains common utility routines  */


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int remove_duplicates(vector<int> & v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int /*k,*/n,nvalues,status=0;
  int *used;
  range_t<int> r=poc_minmax(v,-1);
  
//   printf("initial size=%d\n",v.size());
  
  nvalues=r.max+1;
  
  used=new int[nvalues];
  
  for(n=0;n<nvalues;n++)     used[n]=-1;

//   for(k=0;k<v.size();k++) {
  for(std::vector<int>::iterator s=v.begin(); s<v.end();s++) {
    int m=*s;
    if(used[m]==-1) {
      used[m]=1;
      }
    else {
      v.erase(s);
      }
    }

  delete[] used;
//   for(std::vector<int>::iterator s=v.begin(); s<v.end();s++) {
//     int m=*s;
//     for(std::vector<int>::iterator r=s+1;r<v.end();r++) {
//       int n=*r;
//       if(m==n) v.erase(r);
//       }
//     }
//   printf("final size=%d\n",v.size());
  
  return(status);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_submesh(mesh_t & mesh, vector<int> elements, const char *filename, mesh_t *work)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m, n, status;
  int s;
  
//  mesh_t *work=new mesh_t;
  
  vector<int> nodes;
  vector <vector<int> > partition;

  int *kept=new int[mesh.nvtxs];
  for(n=0;n<mesh.nvtxs;n++) kept[n]=-1;

  for(s=0;s<elements.size();s++) {
    m=elements[s];
    for(int k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
      if(kept[n]==-1) {
        kept[n]=nodes.size();
        nodes.push_back(n);
        }
      }
    }
  
  work->nvtxs=nodes.size();
  work->vertices=new vertex_t[work->nvtxs];
  
  work->ntriangles=elements.size();
  work->triangles=new triangle_t[work->ntriangles];
  
  work->type=mesh.type;
 
  for(s=0;s<elements.size();s++) {
    m=elements[s];
    for(int k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
      work->triangles[s].vertex[k]=kept[n];
      }
    work->triangles[s].ancestor=m;
    }
  
  for(s=0;s<nodes.size();s++) {
    n=nodes[s];
    work->vertices[s].lon=degree_recale(mesh.vertices[n].lon, 0.0);
    work->vertices[s].lat=mesh.vertices[n].lat;
    work->vertices[s].code=0;
    work->vertices[s].ancestor=n;
    }
  
  status=fe_e2n(work);

  status=fe_edgetable(work,0,0);
  
//   status=fe_vertex_element_tables(work);
//   
//   status=fe_GetConnex(*work, partition);
  
  if(filename!=0) status=fe_savemesh(filename,MESH_FILE_FORMAT_TRIGRID,*work);
  
  delete[] kept;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<int> check_pinched(mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, status;
  vector<int> nodes;
  
  status=fe_vertex_element_tables(&mesh);
  
/* *---------------------------------------------------------------------------
  screen pinched boundaries*/
  for(n=0;n<mesh.nvtxs;n++) {
    if(mesh.vertices[n].nngh-mesh.vertices[n].nelmts==2) {
      nodes.push_back(n);
      }
    }
  
  return(nodes);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<int> check_isolated(mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m, s, status;
  vector<int> elements;
  
  vector <vector<int> > partition;
  
/* *---------------------------------------------------------------------------
  screen island triangles*/
  status=fe_vertex_element_tables(&mesh);
  
  status=fe_GetConnex(mesh, partition);
  
  for(s=0;s<partition.size();s++) {
    if(partition[s].size()==1) {
      m=partition[s][0];
      elements.push_back(m);
      printf("orphan triangle %d\n",m);
      }
    partition[s].clear();
    }
    
  partition.clear();
  
  return(elements);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int split(mesh_t & mesh, mesh_t *splitted, vector<int> & target, int rank)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m, n, status;
  int s, smax;
  
  vector<int> nodes,interior,exterior,weird;
  
  int *kept=new int[mesh.nvtxs];

  for(n=0;n<mesh.nvtxs;n++) kept[n]=-1;
  
/**---------------------------------------------------------------------------
  reduce target to first one, insuring connex interior splitted mesh */ 
  nodes.push_back(target[0]);
  kept[target[0]]=0;
  
/**---------------------------------------------------------------------------
  extend vertices list up to rank degree */ 
  for(int k=0; k<rank;k++) {
    smax=nodes.size();
    for(s=0;s<smax;s++) {
      n=nodes[s];
      for(int l=0;l<mesh.vertices[n].nngh;l++) {
        m=mesh.vertices[n].ngh[l];
        if(kept[m]==-1) {
          nodes.push_back(m);
          kept[m]=k;
          }
        }
      }
    }
  delete[] kept;
  
/**---------------------------------------------------------------------------
  create corresponding element list */  
  kept=new int[mesh.ntriangles];
  for(m=0;m<mesh.ntriangles;m++) kept[m]=-1;
    
  if(mesh.vertices[0].nelmts==-1) {
    status=fe_vertex_element_tables(&mesh);
    }
  for(s=0;s<nodes.size();s++) {
    n=nodes[s];
    for(int l=0;l<mesh.vertices[n].nelmts;l++) {
      int m=mesh.vertices[n].elmts[l];
      if(kept[m]==-1) {
        interior.push_back(m);
        kept[m]=1;
        }
      }
    }
  status=fe_submesh(mesh, interior, "interior.nei", &(splitted[1]));
  
/**---------------------------------------------------------------------------
  take care of pinched limits issue */  
  weird=check_pinched(splitted[1]);
  if(weird.size()!=0) {
    for(s=0;s<weird.size();s++) {
      n=splitted[1].vertices[weird[s]].ancestor;
      for(int l=0;l<mesh.vertices[n].nelmts;l++) {
        int m=mesh.vertices[n].elmts[l];
        if(kept[m]==-1) {
          interior.push_back(m);
          kept[m]=1;
          }
        }
      }
    splitted[1].destroy();
    status=fe_submesh(mesh, interior,"interior.nei", &splitted[1]);
    }
  weird.clear();
  
  for(m=0;m<mesh.ntriangles;m++) if(kept[m]==-1) exterior.push_back(m);
  
  status=fe_submesh(mesh, exterior,"exterior.nei", &splitted[0]);
    
/**---------------------------------------------------------------------------
  take care of isolated element issue */  
  weird=check_isolated(splitted[0]);
  if(weird.size()!=0) {
    for(s=0;s<weird.size();s++) {
      int mm=weird[s];
      m=splitted[0].triangles[mm].ancestor;
      interior.push_back(m);
      kept[m]=1;
      }
    splitted[1].destroy();
    status=fe_submesh(mesh, interior,"interior.nei", &splitted[1]);
    exterior.clear();
    splitted[0].destroy();
    for(m=0;m<mesh.ntriangles;m++) if(kept[m]==-1) exterior.push_back(m);
    status=fe_submesh(mesh, exterior,"exterior.nei", &splitted[0]);
    }
  weird.clear();
  
  delete[] kept;

  return(status);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_ImportMeshsize(const mesh_t & parent, float *mainLGP0density, mesh_t & mesh, grid_t *grid, float **density, float *mask, plg_t *polygones, int npolygones, criteria_t & criteria)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int power,status;
  range_t<float> r;
  double step, dmax, ratio;
  vector<plg_t> limits;
  
  float *LGP0density=mesh_LGP0resolution(mesh);
  r=poc_minmax(LGP0density, mesh.ntriangles, 0.);
//  step=r.min/3.0;
  step=r.min/2.1;
  
  dmax=criteria.maxsize*1000.;
  ratio=dmax/step;
  
  power=ceil(log(ratio)/log(2.0)+0.5)-1;
  
  step=dmax/exp(power*log(2.0));
//  step=2000.0;
  
  delete[] LGP0density;
    
/* *---------------------------------------------------------------------------
  create cartesian grid to support mesh density */  
  status=fe_defgrid(grid, polygones, npolygones, limits, step);
  
  grid_t sgrid=map_get_spherical(grid->projection, *grid);
  
  *density=new float[grid->Hsize()];
  
//   int *elts=fe_scan_elements(mesh,sgrid,0,0);
//   status=fe_map(mesh, LGP0density, LGP0, sgrid, elts, *density, *mask);
  
  *mask=-1;
  int *elts=fe_scan_elements(parent,sgrid,0,1);
  status=fe_map(parent, mainLGP0density, LGP0, sgrid, elts, *density, *mask);
  
  
  for(int k=0;k<20;k++) status=map_persistence(sgrid, *density, *mask, 0.);

  criteria.minsize=criteria.shelf_minsize=r.min/1000.;
  criteria.cellsize=step/1000.;

  sgrid.free();
  delete[] elts;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t global(mesh_t & mesh, const mesh_t &parent, float *mainLGP0density, vector<int> & edges, int rank, double dmax, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, status;
  mesh_t final,refined;
  mesh_t splitted[2];
  vector<int> nodes;
  plg_t *polygones;
  int npolygones;
  criteria_t criteria;
  int  RecycleCodes, stopon_EdgeError, stopon_PinchError;
  frame_t frame;
  vector <vector<int> > partition;
  extern int set_density03(criteria_t criteria,grid_t grid, float *density, float mask,  plg_t *polygones, int npolygones);

  grid_t grid;
  float *density, mask;

  printf("#################################################################\n");
  printf("initial size=%d\n",mesh.ntriangles);
  
/**---------------------------------------------------------------------------
  create vertices list */  
  for(int s=0;s<edges.size();s++) {
    n=edges[s];
    nodes.push_back(mesh.edges[n].extremity[0]);
    nodes.push_back(mesh.edges[n].extremity[1]);
    }
    
  status=remove_duplicates(nodes);

  while(nodes.size()!=0) {
/* *---------------------------------------------------------------------------
    partition original mesh*/  
    status=split(mesh, splitted, nodes, rank);
    RecycleCodes=0; stopon_EdgeError=1; stopon_PinchError=0;
    status=fe_codetable1(&splitted[1], RecycleCodes, stopon_EdgeError, stopon_PinchError);
    
   status=fe_minmax(splitted[1], frame);
    
/* *---------------------------------------------------------------------------
    create boundary polygons from partitioned mesh limits*/  
    status=fe_limits2poly(splitted[1], &polygones, &npolygones, true);

    status=fe_geometry(&splitted[1]);
//    if(status!=0) return(-1);
    
/* *---------------------------------------------------------------------------
    get actual partitioned mesh resolution */  
    criteria.maxsize=dmax;
    criteria.niterations=100;
    status=fe_ImportMeshsize(parent, mainLGP0density, splitted[1], &grid, &density, &mask, polygones, npolygones, criteria);

    for(n=0; n< grid.Hsize();n++) if(density[n]!=mask) density[n]=min(density[n], (float) (dmax*1000.) );

    printf("\nmesh density, step #6 (model limits resolution consistency)\n");
    status= set_density03(criteria, grid, density, mask, polygones, npolygones);
    
/* *---------------------------------------------------------------------------
    generate refined mesh */  
    refined=fe_nodit(criteria, grid, density, mask, polygones, npolygones, debug);
    
    delete[] density;
    grid.free();
    
    for(int p=0; p<npolygones; p++) polygones[p].destroy();
    delete[] polygones;
    
    status=fe_codetable1(&refined,RecycleCodes, stopon_EdgeError, stopon_PinchError);
    
    nodes.clear();
    }

  splitted[1].destroy();
  
  if(debug) status=fe_savemesh("interior-refined.nei",MESH_FILE_FORMAT_TRIGRID,refined);
  
  status=fe_vertex_element_tables(&splitted[0]);
  
  status=fe_GetConnex(splitted[0], partition);
  
/* *---------------------------------------------------------------------------
  re-assembly mesh, depends upon number of connex partitions in splitted[0] */  
  switch(partition.size()) {
    case 1:
      stopon_PinchError=0;
      status=fe_codetable1(&splitted[0],RecycleCodes, stopon_EdgeError, stopon_PinchError);
      final=fe_merge(splitted[0], refined,(double) 0.1*dmax);
      refined.destroy();
      partition[0].clear();
      break;
      
    default:
      stopon_PinchError=0;
      for(int s=0; s<partition.size();s++) {
        mesh_t tmp;
        status=fe_submesh(splitted[0], partition[s], 0, &tmp);
//         status=fe_savemesh("part1.nei",MESH_FILE_FORMAT_TRIGRID,tmp);
//         status=fe_savemesh("part2.nei",MESH_FILE_FORMAT_TRIGRID,refined);
        status=fe_codetable1(&tmp,RecycleCodes, stopon_EdgeError, stopon_PinchError);
//         status=fe_savemesh("tmp.nei",MESH_FILE_FORMAT_TRIGRID,tmp);
        final=fe_merge(tmp, refined,(double) 0.1*dmax);
        refined.destroy();
        tmp.destroy();
        refined=final;
        partition[s].clear();
        }
      break;
    }

  partition.clear();
  splitted[0].destroy();

//  status=fe_geometry(&final);
//  status=mesh_resolution(final);
  
  stopon_PinchError=1;
  status=fe_codetable2(&final,RecycleCodes, stopon_EdgeError, stopon_PinchError);
  printf("final size=%d\n",final.ntriangles);
  
  return(final);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<int> fe_OverSized(mesh_t & mesh, int *frontier, double dmax)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   /*i,m,*/n/*,n1,n2*/;
  int   count=0;
  int   nelts,nndes,nedges;

  double /*t,p,*/d/*,t1,t2,p1,p2*/;

  vector<int> edges;

  nelts =mesh.ntriangles;
  nndes =mesh.nvtxs;
  nedges=mesh.nedges;

  printf ("check size criterion...\n");
  printf ("max size for refinement: %lf\n",dmax);

/*-----------------------------------------------------------------------------
  distance criterion */
  for (n=0;n<nedges;n++) {
    if(frontier[n]==1) continue;
    if(mesh.edges[n].code!=MESH_INTERIOR_EDGE)  continue;
    d=mesh.edges[n].L/1000.;
    if(d>1.25*dmax) {
      edges.push_back(n);
      count++;
      }
    }
  printf ("number of edges selected %d (total=%d)\n",count,nedges);
  
  return edges;
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
    
  status=geo_to_projection (proj4, x, y, ndata);
  
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

  mesh_t fe_RefineLimits(mesh_t & mesh, int *locked, double dmax, int equalize, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,status;
  int   count=0, neligible=0;
  vector<plg_t> polygons, resampled;
  double d, radius=dmax*1000.;
  vector<int> edges;
  mesh_t *refined, nodes;
  char *flag=0, *parameters;

  printf ("check edge size criterion...\n");
  printf ("max size for refinement: %lf\n",dmax);
  
  flag=new char[mesh.nedges];
  for (n=0;n<mesh.nedges;n++) flag[n]='I';

  for (int l=0; l<mesh.nlimits; l++) {
    for (int k=0; k<mesh.limits[l].nedges; k++) {
      n=mesh.limits[l].edges[k];
      flag[n]='M';
      }
    }

/*-----------------------------------------------------------------------------
  distance criterion */
  for (n=0;n<mesh.nedges;n++) {
    if(mesh.edges[n].code==MESH_INTERIOR_EDGE) continue;
    neligible++;
    if(locked[n]==1) continue;
    d=mesh.edges[n].L;
    if(d>1.25*radius) {
      edges.push_back(n);
      count++;
      }
    }
    
  printf ("number of boundary edges selected %d (total eligible=%d)\n", count, neligible);
  
  for (int k=0;k<edges.size();k++) {
    n=edges[k];
    flag[n]='T';
    }
  
  status=fe_limits2poly(mesh, polygons, flag, true);
  status=plg_save("extracted.plg", PLG_FORMAT_SCAN, polygons);

  parameters=new char[1024];
  projPJ projection=plg_DefineProjection(polygons, parameters);
  
  status=plg_cartesian(projection, polygons);
  
  for (int l=0; l<polygons.size(); l++) {
    vector<plg_t> p, q;
    plg_t r,s;
    p.push_back(polygons[l]);
/*------------------------------------------------------------------------------
    split polygons to separate open/rigid segments */
    q=plg_split(p, false);
    
    status=plg_save("splitted.plg", PLG_FORMAT_SCAN, q);

    for (int k=0; k<q.size(); k++) {
      if(q[k].flag[0]=='T') {
        r=plg_subdivide(q[k], radius, equalize, PLG_CARTESIAN,debug);
        }
      else r=q[k];
      status=plg_concat(s, r, PLG_CARTESIAN);
      if(status!=0) {
         printf ("concatenation issue\n");
        }
      }
    resampled.push_back(s);
    }
  
  status=plg_spherical(projection, resampled);
    
  point_t *points=new point_t[mesh.nvtxs];

  status=fe_V2cartesian(mesh, parameters);
  
  count=0;
  for(n=0;n<mesh.nvtxs;n++) {
    if(mesh.vertices[n].code==0) {
      points[count].t=mesh.vertices[n].lon;
      points[count].p=mesh.vertices[n].lat;
      count++;
      }
    }

/**-----------------------------------------------------------------------------
  add boundary nodes to interior nodes*/
  status=fe_createnodes(resampled, points, count, nodes);
  delete[] points;
  
  refined=fe_triangulate(nodes, 0, 0, debug);
  
  status=fe_V2spherical(*refined, parameters);
  status=fe_V2spherical(nodes, parameters);
  
  status=fe_geometry(refined);

  status=fe_edgetable(refined,0,0);

  int RecycleCodes=0; int stopon_EdgeError=1; int stopon_PinchError=0;
  status=fe_codetable1(refined, RecycleCodes, stopon_EdgeError, stopon_PinchError);
        
//   mesh->destroy();
  if(debug) status=fe_savenodes("fe_RefineLimits.nod" ,NODE_FILE_FORMAT_TRIGRID, nodes);
  if(debug) status=fe_savemesh ("fe_RefineLimits.nei" ,MESH_FILE_FORMAT_TRIGRID, *refined);

  return (*refined);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_selectcompound_01(mesh_t mesh, double dmax, double cmax, char *poly, int *selected)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  create selection list from :
  
  - max size (in meters)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int   i,m,n,n1,n2;
  int   count=0;
  int   nelts,nndes,nedges;
  int   inside[2];
  edge_t *edges;
  triangle_t *elt;
  double t,p,d,t1,t2,p1,p2;
  double H1,H2,c;
  plg_t *polygones=NULL;
  int npolygones=0;

  elt=mesh.triangles;

  nelts =mesh.ntriangles;
  nndes =mesh.nvtxs;
  nedges=mesh.nedges;
  edges =mesh.edges;
  printf ("check size criterion...\n");
  printf ("max size for refinement: %lf\n",dmax);

  if(dmax==0.0) {
    for (n=0;n<nedges;n++) {
      selected[n]=1;
      count++;
      }
    goto step2;
    }
/*-----------------------------------------------------------------------------
  distance criterion */
  for (n=0;n<nedges;n++) {
    n1=edges[n].extremity[0];
    n2=edges[n].extremity[1];
    t1=mesh.vertices[n1].lon;
    t2=mesh.vertices[n2].lon;
    p1=mesh.vertices[n1].lat;
    p2=mesh.vertices[n2].lat;
    d=geo_distance(t1,p1,t2,p2);
    if(d>dmax) {
      selected[n]=1;
      count++;
      }
    else {
      selected[n]=0;
      }
    }
  printf ("number of edges selected %d (/%d)\n",count,nedges);

 step2:
/*-----------------------------------------------------------------------------
  area criterion */
  printf ("check polygon criterion...\n");
  if(poly!=NULL) plg_load(poly, &polygones, &npolygones);
  else goto step3;
  for (n=0;n<nedges;n++) {
    if(selected[n]==0) continue;
    for(i=0;i<2;i++) {
      m=edges[n].extremity[i];
      t=mesh.vertices[m].lon;
      p=mesh.vertices[m].lat;
      inside[i]=plg_TestInterior(t,p,polygones,npolygones);
      }
    if((inside[0]==-1) || (inside[1]==-1)) {
      selected[n]=0;
      count--;
      }
    }
  printf ("number of edges selected %d (/%d)\n",count,nedges);

 step3:
/*-----------------------------------------------------------------------------
  area criterion */
  printf ("check dH/H criterion...\n");
  if(cmax==0.0) goto end;
  for (n=0;n<nedges;n++) {
    if(selected[n]==0) continue;
    n1=edges[n].extremity[0];
    n2=edges[n].extremity[1];
    H1=mesh.vertices[n1].h;
    H2=mesh.vertices[n2].h;
    c=fabs(H2-H1)/max(H1,H2);
    c=2.*fabs(H2-H1)/(H1+H2);
    if(c<cmax) {
      selected[n]=0;
      count--;
      }
    }
 end:
  printf ("number of edges selected %d (/%d)\n",count,nedges);
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_selectcompound_02(mesh_t mesh, double dmax, char *poly, int *selected)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,m,n,n1,n2,n3,channel;
  int   count=0;
  int   inside[2];
  edge_t *edges;
  double t,p,angle[3];
  plg_t *polygones=NULL;
  int npolygones=0;

  for (n=0;n<mesh.nedges;n++) {
    selected[n]=0;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  select triangles with all vertices being boundary vertex 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for (m=0;m<mesh.ntriangles;m++) {
    channel=1;
    for (i=0;i<3;i++) {
      n=mesh.triangles[m].vertex[i];
      if(mesh.vertices[n].code==0) {
        channel=0;
        break;
        }
      }
    if(channel==0) continue;
    for (i=0;i<3;i++) {
      n=mesh.triangles[m].edges[i];
      if(selected[n]==0) {
        selected[n]=1;
        count++;
        }
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  select bridging edges, i.e. connecting 2 boundary vertices belonging to 
  different limits
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for (n=0;n<mesh.nedges;n++) {
    if(selected[n]==1) continue;
    n1=mesh.edges[n].extremity[0];
    if(mesh.vertices[n1].code==0) continue;
    n2=mesh.edges[n].extremity[1];
    if(mesh.vertices[n2].code==0) continue;
    if(mesh.vertices[n1].code!=mesh.vertices[n2].code) {
/*-----------------------------------------------------------------------------
      bridging 2 different limits*/
      selected[n]=1;
      count++;
      continue;
      }
    if(mesh.edges[n].nshared==2) {
/*-----------------------------------------------------------------------------
      bridging the same limit*/
      selected[n]=1;
      count++;
      continue;
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  select small angles
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for (m=0;m<mesh.ntriangles;m++) {
    for (i=0;i<3;i++) {
      n1=mesh.triangles[m].vertex[i];
      n2=mesh.triangles[m].vertex[(i+2)%3];
      n3=mesh.triangles[m].vertex[(i+1)%3];
      angle[i]=fe_angle(mesh, n1,n2,n3);
      if(fabs(angle[i])<15.*d2r) {
        n=mesh.triangles[m].edges[(i+2)%3];
        if(selected[n]==0) {
          selected[n]=1;
          count++;
          }
        n=mesh.triangles[m].edges[(i+1)%3];
        if(selected[n]==0) {
          selected[n]=1;
          count++;
          }
        }
      }
    }

  goto end;

/*-----------------------------------------------------------------------------
  next step is no more used as it could modify the edges shared by the external
  and internal mesh after splitting, making merge ill-condittionned*/

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  exclude out of area criterion
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(poly!=NULL) plg_load(poly, &polygones, &npolygones);
  else goto end;
  for (n=0;n<mesh.nedges;n++) {
    if(selected[n]==0) continue;
    for(i=0;i<2;i++) {
      m=edges[n].extremity[i];
      t=mesh.vertices[m-1].lon;
      p=mesh.vertices[m-1].lat;
      inside[i]=plg_TestInterior(t,p,polygones,npolygones);
      }
    if((inside[0]==-1) || (inside[1]==-1)) {
      selected[n]=0;
      count--;
      }
    }
  
 end:
  printf ("number of edges selected %d\n",count);
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

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
    "  %s [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Refines a mesh.\n"
    "\n"
    "OPTIONS\n"
    "  -h,--help  Show this help and exit.\n"
    "  --debug  \n"
    "  -m  followed by meshfile\n"
    "  -d  followed by dmax\n"
    "  -e  followed by cmax\n"
    "  -c  enable channels\n"
    "  -o  followed by output. Default: refined-reshaped.nei\n"
    "  -p  followed by poly\n"
    );
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nndes,status;
  int n;
  int option,channels=0;
  char *meshfile=NULL,*output=NULL,*poly=NULL;
  mesh_t mesh,refined,internal,external,*splitted,*final;
  double dmax=0,cmax=0;
  int *selected=0, *targeted=0, *frontier=0, nselected;
  bool debug;
  int element=FE_TRIANGLE;
  int RecycleCodes, stopon_EdgeError, stopon_PinchError,SetLimits;
  int verbose=0;
  
  fct_echo(argc,argv);

  n=1;
  while (n < argc) {
    const char *keyword=argv[n];
    
    if( strcmp(keyword,"--help")==0 or strcmp(keyword,"-h")==0 ){
      print_help(argv[0]);
      exit(0);
      }
    
    if(strncmp("--debug",keyword)==0){
      debug=true;
      n++;
      continue;
      }
    if(strcmp(argv[n],"--quadrangle")==0){
      element=FE_QUADRANGLE;
      n++;
      continue;
      }
    
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'd' :
          sscanf(argv[n+1],"%lf",&dmax);
          n++;
          n++;
          break;

        case 'e' :
          sscanf(argv[n+1],"%lf",&cmax);
          n++;
          n++;
          break;

        case 'c' :
          channels=1;
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

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        n++;
        break;
      }
    
    }

  if(output==0) output=strdup("refined-reshaped.nei");

  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status!=0) TRAP_ERR_EXIT(status,"fe_readmesh(\"%s\",,) error (%d %s)\n",status,strerror(status));
    status=fe_list(&mesh, element);
    if(status!=0) TRAP_ERR_EXIT(status,"unable to build element list of input mesh: fe_list() error %d\n",status);
    }
  else {
    printf("*** Please specify mesh with -m ***\n");
    print_help(argv[0]);
    wexit(-1);
    }

  nndes=mesh.nvtxs;
  status=fe_edgetable(&mesh,0,0);
  if(status!=0) TRAP_ERR_EXIT(status,"unable to build edge list of input mesh: fe_edgetable() error %d\n",status);

/*   status=fe_reshape(mesh,3);  */

  option=0;

  if((dmax!=0) || (cmax!=0) || (poly!=NULL)) {
    option=1;
    if(channels!=0) option=2;
    }
     
//  option=3;
     
     
  switch (option) {
    case 0:
    {
      refined=fe_double(mesh, element);
      printf("double done...\n");
      status=fe_list(&refined, element);
      status=fe_edgetable(&refined,0,0);
      int RecycleCodes=0, stopon_EdgeError=1, stopon_PinchError=0;
//       status=fe_codetable2(&refined, RecycleCodes, stopon_EdgeError, stopon_PinchError);
      status=fe_codetable1(&refined, RecycleCodes, stopon_EdgeError, stopon_PinchError);
      status=fe_savemesh("refined.nei",MESH_FILE_FORMAT_TRIGRID,refined);
      if(element==FE_TRIANGLE) status=fe_reshapeall(refined,3);
      status=fe_savemesh("refined-reshaped.nei",MESH_FILE_FORMAT_TRIGRID,refined);
    }
      break;

    case 1:
      status=fe_codetable2(&mesh,0,1,0);
      if(status!=0) TRAP_ERR_EXIT(status,"unable to build limits table and codes of input mesh: fe_codetable2(,0,1,0) error %d\n",status);
/*       dmax=15.0; */
      selected=new int[mesh.nedges];
      nselected=fe_selectcompound_01( mesh, dmax, cmax, poly, selected);
      if(nselected==0) goto end;
/* *-----------------------------------------------------------------------------
      temporary unsafe mesh, neighbour information limited to boundary vertices*/
      refined=fe_refine( mesh,  selected, 1);
//      status=fe_savemesh("temporary.nei",MESH_FILE_FORMAT_TRIGRID,refined);
      RecycleCodes=1; stopon_EdgeError=1; stopon_PinchError=0;
      SetLimits=0;
//       status =fe_codetable2(&refined, RecycleCodes, stopon_EdgeError, stopon_PinchError);
//       status=fe_list(refined);
      final=fe_node2mesh(refined, debug);
      status =fe_codetable(refined, RecycleCodes, SetLimits, verbose);
      status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID,*final);
      break;

    case 2:
      status=fe_codetable2(&mesh,0,1,0);
      if(status!=0) TRAP_ERR_EXIT(status,"unable to build limits table and codes of input mesh: fe_codetable2(,0,1,0) error %d\n",status);
      selected=new int[mesh.nedges];
      targeted=new int[mesh.nedges];
/*-----------------------------------------------------------------------------
      select edges to be refined */
      nselected=fe_selectcompound_01(mesh, dmax, cmax, poly, selected);
      if(nselected==0) goto end;
      status=fe_selectcompound_02(mesh, dmax, poly, targeted);
/*-----------------------------------------------------------------------------
      extract submesh from selection to allow decent cartesian reshape */
      splitted=fe_split(mesh, selected, targeted, 0, 1, (string) "mesh-refine", debug);
      if(splitted==NULL) TRAP_ERR_EXIT(-2,"unable to split input mesh: fe_split(...) error\n");
/*-----------------------------------------------------------------------------
      refine selected edges; on return, refined mesh contains only nodes and
      boundary informations */
      refined =fe_refine(splitted[1], selected, 0);
/*-----------------------------------------------------------------------------
      RecycleCodes=1 in fe_codetable2, limited operation for node+boundary mesh*/
      RecycleCodes=1; stopon_EdgeError=1; stopon_PinchError=0;
      status  =fe_codetable2(&refined,RecycleCodes, stopon_EdgeError, stopon_PinchError);
      if(status!=0) TRAP_ERR_EXIT(status,"unable to build limits table and codes of refined mesh: fe_codetable2(...) error %d\n",status);
/*-----------------------------------------------------------------------------
      triangle the nodes and boundaries in refined mesh*/
      final=fe_node2mesh(refined, debug);
      if(final==NULL) TRAP_ERR_EXIT(-2,"fe_node2mesh(...) error\n");
      splitted[1]=*final;
      refined=fe_merge(splitted,(double) 0.1, (char *) 0);
      status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID,refined);
      break;

    case 3:
      vector<int> edges;
      bool partial=false;
      frontier=new int[mesh.nedges];
      for (n=0;n<mesh.nedges;n++) {
        frontier[n]=0;
        }
      if(poly!=0) {
        selected=new int[mesh.nedges];
        for (n=0;n<mesh.nedges;n++) {
          selected[n]=1;
          }
        status=fe_selectedges_01(mesh, poly, selected);
        RecycleCodes=0; stopon_EdgeError=1; stopon_PinchError=0;
        status=fe_codetable1(&mesh, RecycleCodes, stopon_EdgeError, stopon_PinchError);
        splitted=fe_split(mesh, selected, targeted, frontier, 1, (string) "mesh-refine", debug);
        mesh.destroy();
        mesh=splitted[1];
//         status=fe_geometry(&mesh);
//         if(status!=0) return(-1);
//         status=fe_edgetable(&mesh,0);
        delete[] selected;
        RecycleCodes=0; stopon_EdgeError=1; stopon_PinchError=0;
        status=fe_codetable1(&mesh, RecycleCodes, stopon_EdgeError, stopon_PinchError);
        status=fe_savemesh("work.nei",MESH_FILE_FORMAT_TRIGRID,mesh);
        partial=true;
        }
      int equalize=1;
      mesh_t work=fe_RefineLimits(mesh, frontier, dmax, equalize, true);
      status=fe_reshapeall(work,3);
      float *mainLGP0density=mesh_LGP0resolution(work);
      edges=fe_OverSized(work, frontier, dmax);
      mesh_t swap;
      int last=0;
      while (edges.size()!=0) {
        if(swap.ntriangles==0) {
          refined=global(work, work, mainLGP0density, edges, 30, dmax, debug);
          }
        else  {
          refined=global(swap, work, mainLGP0density, edges, 30, dmax, debug);
          swap.destroy();
          }
        swap=refined;
        edges.clear();
        edges=fe_OverSized(swap, frontier, dmax);
        status=fe_savemesh("refined.nei",MESH_FILE_FORMAT_TRIGRID,refined);
        if(last>refined.nvtxs) break;
        last=refined.nvtxs;
        }
      status=fe_reshapeall(refined,3);
      status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID,refined);
      mesh.destroy();
      refined.destroy();
      break;
    }

end: STDOUT_BASE_LINE("end of refine ... \n");
  exit(0);
}
