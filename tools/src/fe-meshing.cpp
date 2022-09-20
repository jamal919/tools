
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief mesh edition and diagnostics functions
*/
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>

#include "tools-structures.h"
#include "constants.h"

#include "geo.h"
#include "fe.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "archive.h"
#include "functions.h"
#include "matrix.h"

#include "zapper.h"     /*  rutin.h contains common utility routines  */


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_connex(mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,nn,start,first,low;
  int status=0;
  mesh_t work;
  int *tmp=NULL,*old=NULL;
  
  for(m = 0; m < mesh.ntriangles; m++) {
    int n1=mesh.triangles[m].vertex[0];
    int n2=mesh.triangles[m].vertex[1];
    int n3=mesh.triangles[m].vertex[2];
    bool isolated=( (mesh.vertices[n1].nngh==2) && (mesh.vertices[n2].nngh==2) && (mesh.vertices[n3].nngh==2) );
    if(isolated) printf("orphan triangle %d\n",m);
    }
  
  tmp      =new int[mesh.nvtxs];
  old      =new int[mesh.nvtxs];
/*------------------------------------------------------------------------------
  */
  start=0;
  for(n = 0; n < mesh.nvtxs; n++) tmp[n]=-1;
  
  work.nvtxs=0;
  tmp[start]=work.nvtxs;
  old[work.nvtxs]=start;
  work.nvtxs++;
  
  low=0;
  first=low;
  low=work.nvtxs;
  for(nn = first; nn < work.nvtxs; nn++) {
    n=old[nn];
    for(k=0;k<mesh.vertices[n].nngh;k++) {
      m=mesh.vertices[n].ngh[k];
      if(tmp[m]==-1) {
        tmp[m]=work.nvtxs;
        old[work.nvtxs]=m;
        work.nvtxs++;
        }
      }
    }
  
  if(work.nvtxs!=mesh.nvtxs) {
    printf("%s: mesh vertices=%d, connex=%d\n",__func__,mesh.nvtxs,work.nvtxs);
    for(n = 0; n < mesh.nvtxs; n++) {
      if(tmp[n]==-1) {
        printf("orphan vertex %d\n",n);
        }
      }
    printf("%s: mesh is not connex!!!\n",__func__);
    status=-1;
    }

  delete[] tmp;
  delete[] old;

  work.destroy();
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_GetConnex(mesh_t & mesh, vector <vector<int> > & partition)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  int start, sum=0;
  
  bool *used=new bool[mesh.ntriangles];
  
  bool incremented;
//  vector <vector<int> > partition;
  
  for(int m = 0; m < mesh.ntriangles; m++) {
    int n1=mesh.triangles[m].vertex[0];
    int n2=mesh.triangles[m].vertex[1];
    int n3=mesh.triangles[m].vertex[2];
    bool isolated=( (mesh.vertices[n1].nngh==2) && (mesh.vertices[n2].nngh==2) && (mesh.vertices[n3].nngh==2) );
    if(isolated) printf("orphan triangle %d\n",m);
    }
  
  for(int m = 0; m < mesh.ntriangles; m++) used[m]=false;
  start=0;
  
  do {
    vector<int> *elements=new vector<int>;
    elements->push_back(start);
    used[start]=true;
  
    incremented=false;
    int s=0;
    while(s<elements->size()) {
      int m=(*elements)[s];
      for(int k=0;k <3;k++) {
        int n=mesh.triangles[m].vertex[k];
        for(int l=0;l <mesh.vertices[n].nelmts;l++) {
          int mm=mesh.vertices[n].elmts[l];
          if(!used[mm]) {
            elements->push_back(mm);
            used[mm]=true;
            }
          }
        }
      s++;
      }

    partition.push_back(*elements);
    sum+=elements->size();
    
    printf("connex partition %d: %6d elements, cumulated=%6d (over %d)\n", partition.size(), elements->size(), sum, mesh.ntriangles);

    for(int m = 0; m < mesh.ntriangles; m++) {
      if(!used[m]) {
        start=m;
        incremented=true;
        break;
        }
      }
//     elements->clear();
    delete elements;
    } while(incremented);
    
  if(sum!=mesh.ntriangles) {
    printf("mismatch triangle %d\n",sum);
    return(-1);
    }
  
  delete[] used;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float *mesh_LGP0resolution(mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int i,m;
  
  float *resolution_LGP0=new float[mesh.ntriangles];
  
  status= discretisation_init(&mesh, LGP0);
 
  for(m=0;m<mesh.ntriangles;m++) {
    double L=1.e+10;
    for(i=0;i<3;i++) {
      L=min(L,mesh.triangles[m].l[i]);
      }
    L=0;
    for(i=0;i<3;i++) {
      L+=mesh.triangles[m].l[i]/3.0;
      }
    resolution_LGP0[m]=(float) L;
    }

  return(resolution_LGP0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mesh_resolution(mesh_t & mesh, float* & resolution_LGP0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int i,m;
  
  float mask=1.e+10;
  float *resolution_LGP1=new float[mesh.nvtxs];

  fe_integrale_init();
  
  status= discretisation_init(&mesh, LGP0);
  status= discretisation_init(&mesh, LGP1);
 
  resolution_LGP0=new float[mesh.ntriangles];
  
  for(m=0;m<mesh.ntriangles;m++) {
    double L=1.e+10;
    for(i=0;i<3;i++) {
      L=min(L,mesh.triangles[m].l[i]);
      }
    resolution_LGP0[m]=(float) L;
    }
  
  status=archiving_UGdummy2D("mesh-resolution.nc", mesh, "resolution_LGP0", "m", resolution_LGP0, mask, LGP0);

  status=fe_projection( mesh, resolution_LGP0, LGP0, resolution_LGP1, LGP1);
  
  status=archiving_UGdummy2D("mesh-resolution.nc", mesh, "resolution_LGP1", "m", resolution_LGP1, mask, LGP1);
  
  delete[] resolution_LGP1;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mesh_resolution(mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  float *resolution_LGP0=new float[mesh.ntriangles];

  status=mesh_resolution(mesh, resolution_LGP0);
  
  delete[] resolution_LGP0;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_limits2poly(const mesh_t & mesh, plg_t **polygones, int *npolygones, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,n,status;
  double lon;
  int marine=0, coastal=0;
  
  *npolygones=mesh.nlimits;
  *polygones=new plg_t[*npolygones];
  
  lon=mesh.vertices[mesh.limits[0].vertex[0]].lon;
  
  for (l=0; l<mesh.nlimits; l++) {
    (*polygones)[l].init(mesh.limits[l].nvertex+1,PLG_INIT_SHARED);
    for (k=0; k<mesh.limits[l].nvertex; k++) {
      n=mesh.limits[l].vertex[k];
      (*polygones)[l].t[k]=degree_recale(mesh.vertices[n].lon,lon);
      (*polygones)[l].p[k]=mesh.vertices[n].lat;
      }
    n=mesh.limits[l].vertex[0];
    (*polygones)[l].t[k]=degree_recale(mesh.vertices[n].lon,lon);
    (*polygones)[l].p[k]=mesh.vertices[n].lat;
    (*polygones)[l].flag=new char[mesh.limits[l].nedges];
    for (k=0; k<mesh.limits[l].nedges; k++) {
      n=mesh.limits[l].edges[k];
      if(mesh.edges[n].code==MESH_FLAGGED_EDGE) {
        (*polygones)[l].flag[k]='M';
        marine++;
        }
      else {
        (*polygones)[l].flag[k]='T';
        coastal++;
        }
      }
    }
  
  printf("%s : %d marine edges, %d coastal edges\n",__func__, marine, coastal);
  
  if(debug) status=plg_save("mesh-upgrade.plg", PLG_FORMAT_SCAN, (*polygones), *npolygones);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_limits2poly(const mesh_t & mesh, vector<plg_t> & polygons, char *flag, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,n,status;
  double lon;
  
  if(mesh.limits==0) {
    TRAP_ERR_EXIT(-1,"mesh limits not initiated");
    return (-1);
    }
  
  if(mesh.type==0) lon=mesh.vertices[mesh.limits[0].vertex[0]].lon;
  
  for (l=0; l<mesh.nlimits; l++) {
    plg_t p;
    p.init(mesh.limits[l].nvertex+1,PLG_INIT_SEPARATE);
    for (k=0; k<mesh.limits[l].nvertex; k++) {
      n=mesh.limits[l].vertex[k];
      p.t[k]=mesh.vertices[n].lon;
      if(mesh.type==0) p.t[k]=degree_recale(p.t[k],lon);
      p.p[k]=mesh.vertices[n].lat;
      }
    p.t[k]=p.t[0];
    p.p[k]=p.p[0];
    if(flag!=0) {
      p.flag=new char[mesh.limits[l].nedges];
      for (k=0; k<mesh.limits[l].nedges; k++) {
        n=mesh.limits[l].edges[k];
        p.flag[k]=flag[n];
        }
      }
    polygons.push_back(p);
    }
 
  status=plg_cartesian((projPJ) 0, polygons);
 
  if(debug) status=plg_save("limits2poly.plg", PLG_FORMAT_SCAN, polygons);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_CureOverconnections(mesh_t *mesh,int nghmax, int InteriorsOnly, vector<int> & watch, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  cleave vertices with 8 or more neighbours
  remove interior vertices with 3 neighbours

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int   n,status;
  int   count=0, mark, total;
  int stopon_EdgeError=1, stopon_PinchError=0;
  
  total=0;
  for (n=0;n<mesh->nvtxs;n++) {
    if(mesh->vertices[n].nngh>nghmax) total++;
    }
  printf("%d vertices to be cleaved (nvertices=%d)\n",total, mesh->nvtxs);
  
  mark=max(1,total/10);
  for (n=0;n<mesh->nvtxs;n++) {
    if( (InteriorsOnly==1) && (mesh->vertices[n].code!=0) ) continue;
    if(mesh->vertices[n].nngh>nghmax) {
      int n1,n2;
      if(watch.size()!=0) {
        size_t count=occurence(n,watch);
        debug=(count!=0);
        }
      status=fe_cleave(mesh,n,n1,n2,debug);
      if(debug) {
        watch.push_back(n1);
        watch.push_back(n2);
        }
      count++;
      if(count%mark==0) printf("%3d percent done...\n",(count/mark)*10);
      }
    }

  printf("%3d percent done...\n",100);

  for (n=0;n<mesh->nvtxs;n++) {
    if(mesh->vertices[n].nngh==3) {
      if(mesh->vertices[n].code==0)
        status=fe_removevertex(mesh,n);
        }
    }

/*------------------------------------------------------------------------------
  Remove unused vertices from mesh description*/
  status=fe_cleanvertices(mesh, watch, false);

  status=fe_list(mesh);
/*------------------------------------------------------------------------------
  rebuild edge table*/
  status=fe_edgetable(mesh,0,0);
  printf("#----------------------------------------------------------------\n");
  printf("reconstruct boundary codes and limits table\n");
  status=fe_codetable2(mesh,0, stopon_EdgeError, stopon_PinchError);
  if(status!=0) return(status);
  
  printf("#----------------------------------------------------------------\n");
  printf("reshape mesh\n");
  status=fe_reshapeall(*mesh,3);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_fix_RiverMouth(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *-----------------------------------------------------------------------

  remove sharp river mouth

------------------------------------------------------------------------*/
{
  int   i,m,n,n1,n2,n3,status;
  double d,angle[3],minangle,unsafe=10.0;
  int *flag;

  flag=new int[mesh->nedges];
  for (n=0;n<mesh->nedges;n++) flag[n]=0;

/*-----------------------------------------------------------------------------
  small angles*/
  for (m=0;m<mesh->ntriangles;m++) {
    minangle=1.e+10;
    for (i=0;i<3;i++) {
      n1=mesh->triangles[m].vertex[i];
      n2=mesh->triangles[m].vertex[(i+2)%3];
      n3=mesh->triangles[m].vertex[(i+1)%3];
      angle[i]=fe_angle(*mesh, n1,n2,n3);
      updatemin(&minangle,fabs(angle[i]));
      if(fabs(angle[i])<unsafe*d2r) {
        n=mesh->triangles[m].vertex[i];
        if((mesh->vertices[n].nngh==2)&&(mesh->vertices[n].code!=0)) {
          printf("element %d has minimum angle of %lf (alert set at %lf)\n",m,fabs(angle[i])*r2d,unsafe);
          printf("vertex  %d, #neighbours=%d\n",n,mesh->vertices[n].nngh);
 //         status=fe_removevertex(mesh,n);
          status=fe_disconnectvertex(*mesh, n1, n2);
          status=fe_disconnectvertex(*mesh, n1, n3);
          mesh->vertices[n].nngh=0;
          if(mesh->vertices[n2].code==mesh->vertices[n3].code) {
            flag[mesh->triangles[m].edges[i]]=1;
            }
          break;
          }
        }
      }
    }
  for (n=0;n<mesh->nedges;n++) {
    if(flag[n]==1) {
      n2=mesh->edges[n].extremity[0];
      n3=mesh->edges[n].extremity[1];
      d=fe_distance(*mesh, n2, n3);
      status=fe_mergevertex(mesh, n2, n3, 0);
      }
    }

/*------------------------------------------------------------------------------
  Remove unused vertices from mesh description*/
  status=fe_cleanvertices(mesh, false);

  status=fe_list(mesh);
  status=fe_edgetable(mesh,0,0);
  status=fe_codetable2(mesh,0,1,0);

  delete[] flag;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_clean03(mesh_t *mesh, double threshold)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  remove assles following sharp river mooth cleaning

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int   k,m,n,nn,target,status;
  double d,dmin;

//  printf("#suppress vertces with %d neighbours\n",targetted);
  status=fe_vertex_crosstables02(mesh);

  n=0;
  do {
/*------------------------------------------------------------------------------
    skip interior nodes*/
    if(mesh->vertices[n].code==0) {
      n++;
      continue;
      }
    dmin=+1.e+10;
    for(k=0;k<mesh->vertices[n].nngh;k++) {
      m=mesh->vertices[n].ngh[k];
      if(mesh->vertices[n].code!=mesh->vertices[m].code) continue;
      nn=mesh->vertices[n].edges[k];
/*------------------------------------------------------------------------------
      do not act on bridging edges*/
      if(mesh->edges[n].nshared==2) continue;
      d=fe_distance(*mesh, n, m);
      if(d<dmin) {
        dmin=d;
        target=m;
        }
      }
//      printf("merge %d and %d (nndes=%d)\n",n,target,mesh->nvtxs);
    if(dmin<threshold) {
      if(n<target) {
        status=fe_mergevertex(mesh, n, target, 0);
        }
      }
    n++;
    } while (n<mesh->nvtxs);

/*------------------------------------------------------------------------------
  Remove unused vertices from mesh description*/
  status=fe_cleanvertices(mesh, false);

  status=fe_list(mesh);
  status=fe_edgetable(mesh,0,0);
  status=fe_codetable2(mesh,0,1,0);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_selectedges_01(const mesh_t & mesh, vector<plg_t> & polygons, int *selected, bool polar)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  Detect edges in polygons
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int   i,k,m,n,status;
  int   count=0;
  int   inside[2];
  double t,p;
  projPJ proj=0;
  
/*-----------------------------------------------------------------------------
  assumes (???) lon-lat projection */
  t=polygons[0].t[0];
  for(k=0;k<polygons.size();k++) {
    if(!polygons[k].closed()) {
      printf("checking %dth polygon\n",k);
      TRAP_ERR_EXIT(-1, "selection polygon not closed\n");
      }
    for (n=0;n<polygons[k].npt;n++) {
      polygons[k].t[n]=degree_recale(polygons[k].t[n], t);
      polygons[k].x[n]=polygons[k].t[n];
      }
    }
  
  status=plg_checkAutoSecant(polygons, PLG_SPHERICAL);
  if(status==-1 or polar) {
    printf("check if polygons are polar polygons\n");
    status=CheckPolar(polygons, proj);
    }

  for (n=0;n<mesh.nedges;n++) {
    if(selected[n]==1) count++;
    }

  for (n=0;n<mesh.nedges;n++) {
    if(selected[n]==0) continue;
    for(i=0;i<2;i++) {
      m=mesh.edges[n].extremity[i];
      t=mesh.vertices[m].lon;
      p=mesh.vertices[m].lat;
//       for(k=0;k<npolygones;k++) {
//         for (m=0;m<polygones[k].npt;m++) {
//           polygones[k].x[m]=degree_recale(polygones[k].x[m], t);;
//           }
//         }
      t=degree_recale(t, polygons[0].t[0]);
      if(proj==0) {
        inside[i]=plg_TestInterior(t,p,polygons, PLG_SPHERICAL);
        }
      else {
        double x,y;
        geo_to_projection(proj, p,t, &x, &y);
        inside[i]=plg_TestInterior(x,y,polygons, PLG_CARTESIAN);
        }
      }
    if((inside[0]==PLG_POINT_EXTERIOR) or (inside[1]==PLG_POINT_EXTERIOR)) {
      selected[n]=0;
      count--;
      }
    }
  
//   polygons.clear();
  
  printf ("number of edges selected %d\n",count);
  return(count);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_selectedges_01(const mesh_t & mesh, const char *polygonfile, int *selected, bool polar)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  Detect edges in polygons
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int   count,status;
  vector<plg_t> polygons;
  
  if(polygonfile!=NULL) {
    status=plg_load( polygonfile, polygons);
    if(status!=0) return(-1);
    }
  else return(0);

  count=fe_selectedges_01(mesh, polygons, selected, polar);
  
  plg_destroy_entries(polygons);
  
  polygons.clear();
  
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_selectedges_01(const mesh_t & mesh, vector<plg_t> & polygons, int *selected, int option, bool polar)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  Detect edges in polygons
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int   i,k,m,n,status;
  int   count=0;
  int   inside[2];
  double t,p;
  projPJ proj=0;
  
/*-----------------------------------------------------------------------------
  assumes (???) lon-lat projection */
  t=polygons[0].t[0];
  for(k=0;k<polygons.size();k++) {
    if(!polygons[k].closed()) {
      printf("checking %dth polygon\n",k);
      TRAP_ERR_EXIT(-1, "selection polygon not closed\n");
      }
    for (n=0;n<polygons[k].npt;n++) {
      polygons[k].t[n]=degree_recale(polygons[k].t[n], t);
      polygons[k].x[n]=polygons[k].t[n];
      }
    }
  
  status=plg_checkAutoSecant(polygons, PLG_SPHERICAL);
  if(status==-1 or polar) {
    printf("check if polygons are polar polygons\n");
    status=CheckPolar(polygons, proj);
    }

  for (n=0;n<mesh.nedges;n++) {
    if(selected[n]==1) count++;
    }

  for (n=0;n<mesh.nedges;n++) {
    if(selected[n]==0) continue;
    for(i=0;i<2;i++) {
      m=mesh.edges[n].extremity[i];
      t=mesh.vertices[m].lon;
      p=mesh.vertices[m].lat;
//       for(k=0;k<npolygones;k++) {
//         for (m=0;m<polygones[k].npt;m++) {
//           polygones[k].x[m]=degree_recale(polygones[k].x[m], t);;
//           }
//         }
      t=degree_recale(t, polygons[0].t[0]);
      if(proj==0) {
        inside[i]=plg_TestInterior(t,p,polygons, PLG_SPHERICAL);
        }
      else {
        double x,y;
        geo_to_projection(proj, p,t, &x, &y);
        inside[i]=plg_TestInterior(x,y,polygons, PLG_CARTESIAN);
        }
      }
    if(option==0) {
      if((inside[0]==PLG_POINT_EXTERIOR) or (inside[1]==PLG_POINT_EXTERIOR)) {
        selected[n]=0;
        count--;
        }
      }
    else {
      if((inside[0]==PLG_POINT_EXTERIOR) and (inside[1]==PLG_POINT_EXTERIOR)) {
        selected[n]=0;
        count--;
        }
      }
    }
  
//   polygons.clear();
  
  printf ("number of edges selected %d\n",count);
  return(count);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_selectedges_02(const mesh_t & mesh, int *selected, bool initialise)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  Detect channel edges
-----------------------------------------------------------------------------*/
{
  int   i,m,n,channel;
  int   count=0;

  if(initialise) {
    for (n=0;n<mesh.nedges;n++) {
      selected[n]=0;
      }
    }

/*-----------------------------------------------------------------------------
  triangles with all nodes being boundaries*/
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

  printf ("number of channel edges %d\n",count);
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_selectedges_03(const mesh_t & mesh, int *selected, bool initialise)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  Detect bridging edges
-----------------------------------------------------------------------------*/
{
  int   i,m,n,n1,n2;
  int   count=0;
  vector <int> elements;

  printf ("#checking edges connecting 2 different boundaries (bridge)\n");
  if(initialise) {
    for (n=0;n<mesh.nedges;n++) {
      selected[n]=0;
      }
    }

/*-----------------------------------------------------------------------------
  bridging edges*/
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
  printf ("number of bridging edges %d\n",count);
 
  for(int loop=0;loop<5;loop++) {
  for (m=0;m<mesh.ntriangles;m++) {
    for (i=0;i<3;i++) {
      n=mesh.triangles[m].edges[i];
      if(selected[n]==1) elements.push_back(m);
      }
    }
  for (int k=0;k<elements.size();k++) {
    m=elements[k];
    for (i=0;i<3;i++) {
      n=mesh.triangles[m].edges[i];
      selected[n]=1;
      }
    }
  elements.clear();
  }
  count=0;
  for (n=0;n<mesh.nedges;n++) {
    if(selected[n]==1) count++;
    }
  printf ("number of bridging edges %d\n",count);
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_selectedges_04(const mesh_t & mesh, double dmax, int *selected)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  Detect wide edges
-----------------------------------------------------------------------------*/
{
  int   n,n1,n2;
  int   count=0;
  int   nelts,nndes,nedges;
  edge_t *edges;
  double d,t1,t2,p1,p2;

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
    return(count);
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


  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int fe_Echk_crossing(mesh_t mesh,int *targeted)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  Detect crossing edges
-----------------------------------------------------------------------------*/
{
  int m1,m2,n,nn;

  double *dx,*dy;
  double epsilon,error,nx,ny,x,x1,x2,y;
  double ddx,ddy,ddn;
  double alpha1,alpha2;
/*------------------------------------------------------------------------------
  */

  dx=new double[mesh.nedges];
  dy=new double[mesh.nedges];

  error=1.e-04;
  epsilon=1.e-12;

  for (n=0;n<mesh.nedges;n++) {
//    if(targeted[n]==0) continue;
    m1=mesh.edges[n].extremity[0];
    m2=mesh.edges[n].extremity[1];
    x1=mesh.vertices[m1].lon;
    x2=degree_recale(mesh.vertices[m2].lon,x1);
    dx[n]=x2-x1;
    dy[n]=mesh.vertices[m2].lat-mesh.vertices[m1].lat;
    }

  for (n=0;n<mesh.nedges;n++) {
    if(targeted[n]==0) continue;
    m1=mesh.edges[n].extremity[0];
    x=mesh.vertices[m1].lon;
    y=mesh.vertices[m1].lat;
    for (nn=n+1;nn<mesh.nedges;nn++) {
//      if(targeted[nn]==0) continue;
      m2=mesh.edges[nn].extremity[0];
/*------------------------------------------------------------------------------
      check if nn line cross n edge*/
      nx=-dy[nn];
      ny=+dx[nn];
//      ddx=degree_recale(mesh.vertices[m2].lon,x)-x;
      ddx=mesh.vertices[m2].lon-x;
      if(ddx > 180.0) ddx-=360.0;
      else if(ddx < -180.0) ddx+=360.0;
      ddy=mesh.vertices[m2].lat-y;
      ddn=dx[n]*nx+dy[n]*ny;
      if(ddn==0.0) continue; /*parallel*/
      alpha1=(ddx*nx+ddy*ny)/ddn;
      if ((alpha1 < 0.0) || (alpha1 > 1.0)) continue;
/*------------------------------------------------------------------------------
      check if n line cross nn edge*/
      nx=-dy[n];
      ny=+dx[n];
      alpha2=(ddx*nx+ddy*ny)/ddn;
      if ((alpha2 >= 0.0)&&(alpha2 < 1.0)) {
        if(m2==mesh.edges[n].extremity[0]) continue;
        if(m2==mesh.edges[n].extremity[1]) continue;
        m2=mesh.edges[nn].extremity[1];
        if(m2==mesh.edges[n].extremity[0]) continue;
        if(m2==mesh.edges[n].extremity[1]) continue;
        printf("cross edges: %d %d\n",n,nn);
        printf("%lf %lf \n %lf %lf\n\n",mesh.vertices[mesh.edges[n].extremity[0]].lon,mesh.vertices[mesh.edges[n].extremity[0]].lat,
                           mesh.vertices[mesh.edges[n].extremity[1]].lon,mesh.vertices[mesh.edges[n].extremity[1]].lat);
        printf("%lf %lf \n %lf %lf\n",mesh.vertices[mesh.edges[nn].extremity[0]].lon,mesh.vertices[mesh.edges[nn].extremity[0]].lat,
                           mesh.vertices[mesh.edges[nn].extremity[1]].lon,mesh.vertices[mesh.edges[nn].extremity[1]].lat);
        }
      }
    }

  free(dx);
  free(dy);

  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_Echk_angles(mesh_t & mesh, int *selected, double threshold, int mode, bool exterior_only, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  Detect flat or small angles (edge-wise)
------------------------------------------------------------------------------*/
{
  int   i,m,n,n1,n2,n3;
  int   count=0;
  double angle[3];
  
  for (n=0;n<mesh.nedges;n++) {
    selected[n]=0;
    }

/*-----------------------------------------------------------------------------
  small angles*/
  for (m=0;m<mesh.ntriangles;m++) {
    for (i=0;i<3;i++) {
      n1=mesh.triangles[m].vertex[i];
      n2=mesh.triangles[m].vertex[(i+2)%3];
      n3=mesh.triangles[m].vertex[(i+1)%3];
      angle[i]=fe_angle(mesh, n1,n2,n3);
      }
    for (i=0;i<3;i++) {
      switch(mode) {
        case 0:
          if(fabs(fabs(angle[i])-M_PI)<threshold*d2r) {
            n=mesh.triangles[m].edges[i];
            if(mesh.edges[n].code==MESH_INTERIOR_EDGE and exterior_only) continue;
            if(selected[n]==0) {
              selected[n]=1;
              count++;
              break;
              }
            }
          break;
        case 1:
          if(fabs(angle[i])<threshold*d2r) {
            n1=mesh.triangles[m].edges[(i+2)%3];
            n2=mesh.triangles[m].edges[(i+1)%3];
            if(mesh.edges[n1].code==MESH_INTERIOR_EDGE and mesh.edges[n2].code!=MESH_INTERIOR_EDGE)
              n=n2;
            else if(mesh.edges[n1].code!=MESH_INTERIOR_EDGE and mesh.edges[n2].code==MESH_INTERIOR_EDGE)
              n=n1;
            else continue;
            if(selected[n]==0) {
              selected[n]=1;
              count++;
              break;
              }
            }
          break;
        }
      }
    }
  printf ("number of small angle edges (threshold=%lf) : %d\n",threshold,count);
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_Tchk_angles(mesh_t & mesh, int *selected, double threshold, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  Detect small angles (vertex-wise)
-----------------------------------------------------------------------------*/
{
  int   i,m,n,n1,n2,n3;
  int   count=0;
  double angle[3],minangle,unsafe=5.0;

  if(selected==0) selected=new int [mesh.ntriangles];
  for (m=0;m<mesh.ntriangles;m++) {
    selected[m]=0;
    }

/*-----------------------------------------------------------------------------
  small angles*/
  for (m=0;m<mesh.ntriangles;m++) {
    minangle=1.e+10;
    for (i=0;i<3;i++) {
      n1=mesh.triangles[m].vertex[i];
      n2=mesh.triangles[m].vertex[(i+2)%3];
      n3=mesh.triangles[m].vertex[(i+1)%3];
      angle[i]=fe_angle(mesh, n1, n2, n3);
      updatemin(&minangle,fabs(angle[i]));
      if(fabs(angle[i])<threshold*d2r) {
        if(selected[m]==0) {
          selected[m]=1;
          count++;
//          break;
          }
        }
      if(fabs(angle[i])<unsafe*d2r) {
        printf("element %d has minimum angle of %lf (alert set at %lf)\n",m,fabs(angle[i])*r2d,unsafe);
        n=mesh.triangles[m].vertex[i];
        printf("vertex  %d, #neighbours=%d\n",n,mesh.vertices[n].nngh);
        }
      }
    if(minangle<unsafe*d2r) {
      printf("element %d has minimum angle of %lf (alert set at %lf)\n",m,minangle*r2d,unsafe);
      }
    }
  printf ("number of small angle elements %d\n",count);
  
  
  
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_AspectRatio(mesh_t mesh, int m)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  compute triangle aspect ratio
-----------------------------------------------------------------------------*/
{
  int   i;
  double a,S,f;
  
  S=mesh.triangles[m].TrueArea;
  a=0;
  for (i=0;i<3;i++) {
    a+=mesh.triangles[m].l[i];
    }
  a/=3.0;
  f=4.*S/(sqrt(3.0)*a*a);
  return(f);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_Tchk_AspectRatio(mesh_t & mesh, int *selected, double min_ratio,float *flags)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  Detect bad aspect ratio
-----------------------------------------------------------------------------*/
{
  int   i,m,n;
  int   count=0;
  double a,S,f;
  double rms=0,mean=0,alert_ratio=0.1;

  for (n=0;n<mesh.ntriangles;n++) {
    selected[n]=0;
    }

/*-------------------------------------------------------------------------------
  equilateral triangle surface = sqrt(3*a⁴/16)= sqrt(3)*a²/4

  l=perimeter/3=average side length
  aspect factor= 4*S/sqrt(3)*l²
  aspect factor is 1 for equilateral triangle

  flat element:
  a=c
  b=1.5*a

  s=7/4a
  l=7/6a

  S=a²sqrt(7/4 * (7/4-1)(7/4-1)(7/4-3/2))=a²sqrt(7/4 * 9/16 * 1/2)=a²sqrt(7/2)/8
  aspect factor= sqrt(7./2.)/2./sqrt(3.) / (7./6.)²



  flat element:
  a=c
  b=2*x*a

  s=(2+2x)a/2=a(1+x)
  l=2/3*a(1+x)

  S=sqrt(a(1+x)*ax*ax*a(1-x))=a²sqrt(x²(1-x²))=a²x*sqrt(1-x²)
  aspect factor= x*sqrt(1-x²)*4/sqrt(3)/(2/3(1+x))²=x*sqrt(1-x²)*9/sqrt(3)/(1+x)²
  
---------------------------------------------------------------------------*/

  for (m=0;m<mesh.ntriangles;m++) {
    S=mesh.triangles[m].TrueArea;
    a=0;
    for (i=0;i<3;i++) {
      a+=mesh.triangles[m].l[i];
      }
    a/=3.0;
    f=4.*S/(sqrt(3.0)*a*a);
    rms+=f*f;
    mean+=f;
    if(f<min_ratio) {
      if(selected[m]==0) {
        selected[m]=1;
        count++;
        }
      }
    if(f<alert_ratio) {
      printf("element %d has aspect ratio of %lf, smaller than %lf\n",m,f,alert_ratio);
      if(flags!=0) {
        for(i=0;i<3;i++) {
          n=mesh.triangles[m].vertex[i];
          flags[n]=1;
          }
        }
      }
    }

  mean/=(double) mesh.ntriangles;
  rms/=(double) mesh.ntriangles;
  rms=sqrt(rms-mean*mean);
  printf ("number of bad aspect ratio elements %d\n",count);
  printf ("aspect factor statistics : mean=%6.2f rms=%6.2f\n",mean,rms);
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_chkelement_02(mesh_t mesh, int *selected)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,m,n,channel;
  int   count=0;

  for (n=0;n<mesh.ntriangles;n++) {
    selected[n]=0;
    }

/*-----------------------------------------------------------------------------
  triangles with all nodes being boundaries*/
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
    if(selected[m]==0) {
      selected[m]=1;
      count++;
      }
    }

  printf ("number of channel elements %d\n",count);
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_chkvertex_consistency(mesh_t *mesh, int force)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  Detect dissymmetric (erroneous) neihgbour lists
-----------------------------------------------------------------------------*/
{
  int  m,n,k,l,status;
  int  count=0;
  const char *msg;

/*-------------------------------------------------------------------------------
  check reciprocity of node connexions*/
  for (n=0;n<mesh->nvtxs;n++) {
    for(k=0;k<mesh->vertices[n].nngh;k++) {
      m=mesh->vertices[n].ngh[k];
      if(m==n) {
        count++;
        printf ("#node %d has inconsistent neighbour list (self-connected): %d %d\n",n,k,m);
        if(force==1) {
          fe_disconnectvertex(*mesh, n, n);
          }
        }
      }
    }

  for (n=0;n<mesh->nvtxs;n++) {
    for(k=0;k<mesh->vertices[n].nngh;k++) {
      m=mesh->vertices[n].ngh[k];
      l=fe_anb(*mesh,n,m);
      if(l==-1) {
        count++;
        printf ("#node %d has inconsistent neighbour list (un-symmetrical connection): %d %d\n",n,k,m);
        if(force==1) {
          status=fe_insertvertex(mesh, n);
          }
        }
      }
    }

  if(count==0) {
    msg="ok";
    }
  else {
    msg="trouble";
    }
  if(count!=0) printf ("#nodes with inconsistent neighbour list=%d (%s)\n",count,msg);
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_chkvertex_00(mesh_t mesh, int *selected, int *histogram)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  Detect dissymmetric (erroneous) neihgbour lists
-----------------------------------------------------------------------------*/
{
  int   m,n,k,l;
  int   count=0;

  for (n=0;n<mesh.nvtxs;n++) {
    selected[n]=0;
    }

/*-------------------------------------------------------------------------------
  check reciprocity of node connexions*/
  for (n=0;n<mesh.nvtxs;n++) {
    for(k=0;k<mesh.vertices[n].nngh;k++) {
      m=mesh.vertices[n].ngh[k];
      l=fe_anb(mesh,n,m);
      if(l==-1) {
        if(selected[n]==0) {
          selected[n]=1;
          count++;
          printf ("#node %d has inconsistent neighbour list: %d %d\n",n,k,m);
          break;
          }
        }
      }
    }

  printf ("#nodes with inconsistent neighbour list=%d\n",count);
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_chkvertex_01(mesh_t mesh, int *selected, int *histogram)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  connections statistics
-----------------------------------------------------------------------------*/
{
  int   n,k;
  int   count=0;
  int   kmax,kmin;

  for (n=0;n<mesh.ntriangles;n++) {
    selected[n]=0;
    }

/*-----------------------------------------------------------------------------
  full node set*/
  kmin=10000;
  kmax=0;
  for (k=0;k<20;k++) {
    histogram[k]=0;
    }
  for (n=0;n<mesh.nvtxs;n++) {
    histogram[mesh.vertices[n].nngh]++;
    kmin=min(kmin,mesh.vertices[n].nngh);
    kmax=max(kmax,mesh.vertices[n].nngh);
    }

  printf ("Full vertex set histogram, count=%d\n",mesh.nvtxs);
  for (k=kmin;k<kmax;k++) {
    printf ("k=%3d: %6d %6.1f\n",k,histogram[k],100.*histogram[k]/mesh.nvtxs);
    }

/*-----------------------------------------------------------------------------
  interior node set*/
  kmin=10000;
  kmax=0;
  count=0;
  for (k=0;k<20;k++) {
    histogram[k]=0;
    }
  for (n=0;n<mesh.nvtxs;n++) {
    if(mesh.vertices[n].code==0) {
      histogram[mesh.vertices[n].nngh]++;
      kmin=min(kmin,mesh.vertices[n].nngh);
      kmax=max(kmax,mesh.vertices[n].nngh);
      count++;
      }
    }

  printf ("Interior vertex set histogram, count=%d\n",count);
  for (k=kmin;k<kmax;k++) {
    printf ("k=%3d: %6d %6.1f\n",k,histogram[k],100.*histogram[k]/count);
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_supress_diamond(mesh_t *mesh, int targetted, int nloopmax, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Suppress vertices with the tragetted number of neighbours by merging with
  their closes neighbour

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int   m,n,k,target,status;
  int   count=0,max_nngh;
  int   mark,total,loop,processed;
  double d,dmin;

  printf("#suppress vertices with %d neighbours (nloops=%d)\n",targetted,nloopmax);
  
  
  for(loop=0;loop<nloopmax;loop++){
    
    processed=0;
/*------------------------------------------------------------------------------
  first try to cleave neighbours*/
  total=0;
  for(n=0;n<mesh->nvtxs;n++) {
    if(mesh->vertices[n].code!=0) {
      n++;
      continue;
      }
    if(mesh->vertices[n].nngh==targetted) {
      total++;
      }
    }
  printf("loop %d: %d diamonds to be processed, try first cleaving method (nndes=%d)\n",loop,total,mesh->nvtxs);
  
  mark=max(total/10,1);
  
  count=0;
  n=0;
  do {
    if(mesh->vertices[n].code!=0) {
      n++;
      continue;
      }
    if(mesh->vertices[n].nngh==targetted) {
      max_nngh=0;
      for(k=0;k<mesh->vertices[n].nngh;k++) {
        m=mesh->vertices[n].ngh[k];
        if(mesh->vertices[m].code!=0) continue;
        if(max_nngh<mesh->vertices[m].nngh) {
          max_nngh=mesh->vertices[m].nngh;
          target=m;
          }
        }
      if(max_nngh>6) {
        status=fe_cleave(mesh, target,debug);
        if(status==0) {
          count++;
          if(count%mark==0) printf("%3d percent done...\n",(count/mark)*10);
          }
        }
      }
    n++;
    } while (n<mesh->nvtxs);
  printf("%d diamonds processed by cleaving neighbours (nvertices=%d)\n",count,mesh->nvtxs);
  processed+=count;
  if(debug) status=fe_savemesh("mesh-doctor.tmp.nei",MESH_FILE_FORMAT_TRIGRID,*mesh);

/*------------------------------------------------------------------------------
  then try to merge with the nearest neighbour*/
  total=0;
  for(n=0;n<mesh->nvtxs;n++) {
    if(mesh->vertices[n].code!=0) {
      n++;
      continue;
      }
    if(mesh->vertices[n].nngh==targetted) {
      total++;
      }
    }
  printf("loop %d: %d diamonds to be processed, try now merging method (nvertices=%d)\n",loop, total,mesh->nvtxs);
  
  mark=max(total/10,1);
  
  count=0;
  n=0;
  do {
    if(mesh->vertices[n].code!=0) {
      n++;
      continue;
      }
    if(mesh->vertices[n].nngh==targetted) {
      dmin=+1.e+10;
      for(k=0;k<mesh->vertices[n].nngh;k++) {
        m=mesh->vertices[n].ngh[k];
        d=fe_distance(*mesh, n, m);
        if(d<dmin) {
          dmin=d;
          target=m;
          }
        }
//      printf("merge %d and %d (nndes=%d)\n",n,target,mesh->nvtxs);
      status=fe_mergevertex(mesh, n, target, 0);
      if(status==0) {
        count++;
        if(count%mark==0) printf("%3d percent done...\n",(count/mark)*10);
        }
      }
    n++;
    } while (n<mesh->nvtxs);
    
  printf("%d diamonds processed by merging neighbours (nvertices=%d)\n",count,mesh->nvtxs);
  processed+=count;

  total=0;
  for(n=0;n<mesh->nvtxs;n++) {
    if(mesh->vertices[n].code!=0) {
      n++;
      continue;
      }
    if(mesh->vertices[n].nngh==targetted) {
      total++;
      }
    }
  if(total!=0) {
    if(loop==nloopmax-1)
      printf("loop %d: %d diamonds left after processing (need manual editing, or run mesh-doctor again) (nndes=%d)\n",loop, total,mesh->nvtxs);
    else
      printf("loop %d: %d diamonds left after processing (nndes=%d)\n",loop, total,mesh->nvtxs);
   }
  else {
    printf("loop %d: all diamonds removed\n",loop, total,mesh->nvtxs);
    }
  
/*------------------------------------------------------------------------------
  Remove unused vertices from mesh description*/
  status=fe_cleanvertices(mesh, false);
  if(debug) status=fe_savemesh("mesh-doctor.tmp.nei",MESH_FILE_FORMAT_TRIGRID,*mesh);

/*------------------------------------------------------------------------------
  Rebuild triangle structure */
  status=fe_list(mesh);
  status=fe_edgetable(mesh,0,0);
  printf("#----------------------------------------------------------------\n");
  printf("reconstruct boundary codes and limits table\n");
  status=fe_codetable2(mesh,0,1,0);

/*------------------------------------------------------------------------------
  Reshape mesh*/
  printf("#----------------------------------------------------------------\n");
  printf("reshape mesh\n");
  if(debug) status=fe_savemesh("mesh-doctor.tmp.nei",MESH_FILE_FORMAT_TRIGRID,*mesh);
  status=fe_reshapeall(*mesh,3);
  if(processed==0) break;
  }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_load_criteria(const char *filename, criteria_t *criteria, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in;
  char control,*keyword=0,*line,param[1024];
  int count=0;

  line=new char[1024];
  
  in=fopen(filename,"r");
  
  if(in==0) {
    count=errno;
    if(verbose>=0) printf("%s: cannot open %s (%d %s)\n",__func__,filename,count,strerror(count));
    return count;
    }

  while (true) {

    fgets(line,1024,in);
    count++;
    if(feof(in)) goto finished;
    
    if(strlen(line)==0) continue;
    
    control=line[0];
    
    if(keyword!=0) free(keyword);
//    keyword=strdup(&(line[1]));
// /**-----------------------------------------------------------------------------
//     remove end-of-line tails */
//     keyword[strlen(keyword)-1]='\0';
    keyword=(char *) malloc(strlen(line)*sizeof(char));
    sscanf(&line[1],"%s", keyword);
    
/**-----------------------------------------------------------------------------
    remove blank tails */
    for(size_t k=strlen(keyword)-1;k>0;k--) {
      if(keyword[k]!=' ') break;
      keyword[k]=0;
      }

    if(control != '#') continue;

    fgets(param,1024,in);
    count++;

/*------------------------------------------------------------------------------
    bathymetry range, min */
    if(strcmp(keyword, "OPEN_MIN_DEPTH")==0) {
      sscanf(param,"%lf", &(criteria->hmin));
      continue;
      }
    if(strcmp(keyword, "OPEN_MIN_DEPTH_2")==0) {
      sscanf(param,"%lf", &criteria->hmin2);
      continue;
      }
/*------------------------------------------------------------------------------
    bathymetry range, max */
    if(strcmp(keyword, "OPEN_MAX_DEPTH")==0) {
      sscanf(param,"%lf", &criteria->hmax);
      continue;
      }
    if(strcmp(keyword, "OPEN_MIN_ANGLE")==0) {
      sscanf(param,"%lf", &criteria->minangle);
      continue;
      }
/*----------------------------------------------------------------------------
    depth limit between open waters and shelf   */
    if(strcmp(keyword, "OPEN_SHELF_LIMIT")==0) {
      sscanf(param,"%lf", &criteria->open_shelf_limit);
      continue;
      }
    if(strcmp(keyword, "SHELF_SHORE_LIMIT")==0) {
      sscanf(param,"%lf", &criteria->shelf_shore_limit);
      continue;
      }

/*------------------------------------------------------------------------------
    minimum H'/H value */
    if(strcmp(keyword, "OPEN_MIN_RATIO")==0) {
      sscanf(param,"%lf", &criteria->minratio);
      continue;
      }
      
/*------------------------------------------------------------------------------
    minimum element size in deep ocean */
    if(strcmp(keyword, "OPEN_MIN_RADIUS")==0) {
      sscanf(param,"%lf", &criteria->minsize);
      continue;
      }
/*------------------------------------------------------------------------------
    maximum element size in deep ocean */
    if(strcmp(keyword, "OPEN_MAX_RADIUS")==0) {
      sscanf(param,"%lf", &criteria->maxsize);
      continue;
      }
/*------------------------------------------------------------------------------
    minimum element size in shelf seas */
    if(strcmp(keyword, "SHELF_MIN_RADIUS")==0) {
      sscanf(param,"%lf", &criteria->shelf_minsize);
      continue;
      }
/*------------------------------------------------------------------------------
    maximum element size in shelf seas */
    if(strcmp(keyword, "SHELF_MAX_RADIUS")==0) {
      sscanf(param,"%lf", &criteria->shelf_maxsize);
      continue;
      }
/*------------------------------------------------------------------------------
    minimum element size in shelf seas */
    if(strcmp(keyword, "SHORE_MIN_RADIUS")==0) {
      sscanf(param,"%lf", &criteria->shelf_minsize);
      continue;
      }
/*------------------------------------------------------------------------------
    maximum element size in shelf seas */
    if(strcmp(keyword, "SHORE_MAX_RADIUS")==0) {
      sscanf(param,"%lf", &criteria->shore_maxsize);
      continue;
      }
      
/*------------------------------------------------------------------------------
    background cgrid resolution (must be less than targeted mesh resolution) */
    if(strcmp(keyword, "OPEN_DX")==0) {
      sscanf(param,"%lf", &criteria->cellsize);
      continue;
      }
/*------------------------------------------------------------------------------
    gravity wavelength factor T(s)/15/1000 */
    if(strcmp(keyword, "SAMPLE_SCALE_1")==0) {
      sscanf(param,"%lf", &criteria->factor1);
      continue;
      }
/*------------------------------------------------------------------------------
    H'/H factor, ~0.4 */
    if(strcmp(keyword, "SAMPLE_SCALE_2")==0) {
      sscanf(param,"%lf", &criteria->factor2);
      continue;
      }
/*------------------------------------------------------------------------------
    CFL factor, not used yet */
    if(strcmp(keyword, "SAMPLE_SCALE_3")==0) {
      sscanf(param,"%lf", &criteria->factor3);
      continue;
      }
/*------------------------------------------------------------------------------
    surface wave factor */
    if(strcmp(keyword, "SAMPLE_SCALE_4")==0) {
      sscanf(param,"%lf", &criteria->factor4);
      continue;
      }

/*----------------------------------------------------------------------------
    surface wave period    */
    if(strcmp(keyword, "SURFACE_WAVE_PERIOD")==0) {
      sscanf(param,"%lf", &criteria->surface_wave_period);
      continue;
      }

// /*----------------------------------------------------------------------------
//     surface wave period    */
//     if(strcmp(keyword, "SURFACE_WAVE_CFL")==0) {
//       int nitems, dum;
//       nitems=sscanf(param,"%d", &dum);
//       if(nitems==1) {
//         criteria->surface_wave_cfl=(dum==1);;
//         }
//       else {
//         criteria->surface_wave_cfl=(strstr(param,"T")!=0);
//         }
//       continue;
//       }

/*----------------------------------------------------------------------------
    surface wave period    */
    if(strcmp(keyword, "SURFACE_WAVE_DT")==0) {
      sscanf(param,"%lf", &criteria->surface_wave_dt);
      continue;
      }

/*------------------------------------------------------------------------------
    1 * (gravity==true) + 2 * (H'/H==true) + 4 * (CFL==true) + 8 * (waves==true) */
    if(strcmp(keyword, "OPEN_SAMPLE_MODE")==0) {
      sscanf(param,"%d", &criteria->mode);
      continue;
      }
/*------------------------------------------------------------------------------
    allow/disallow boundary resampling accordingly with criteria */
    if(strcmp(keyword, "RE_SAMPLE")==0) {
      int nitems, dum;
      nitems=sscanf(param,"%d", &dum);
      if(nitems==1) {
        criteria->resample_obsolete=(dum==1);;
        }
      else {
        criteria->resample_obsolete=(strstr(param,"T")!=0);
        }
      continue;
      }
/*------------------------------------------------------------------------------
    allow/disallow boundary resampling accordingly with criteria */
    if(strcmp(keyword, "RE_SAMPLE_RIGID")==0) {
      int nitems, dum;
      nitems=sscanf(param,"%d", &dum);
      if(nitems==1) {
        criteria->resample_rigidlimits=(dum==1);;
        }
      else {
        criteria->resample_rigidlimits=(strstr(param,"T")!=0);
        }
      continue;
      }
/*------------------------------------------------------------------------------
    allow/disallow boundary resampling accordingly with criteria */
    if(strcmp(keyword, "RE_SAMPLE_OPEN")==0) {
      int nitems, dum;
      nitems=sscanf(param,"%d", &dum);
      if(nitems==1) {
        criteria->resample_openlimits=(dum==1);;
        }
      else {
        criteria->resample_openlimits=(strstr(param,"T")!=0);
        }
      continue;
      }
/*------------------------------------------------------------------------------
    allow/disallow boundary resampling accordingly with criteria */
    if(strcmp(keyword, "#SELECTIVE_OPEN_SAMPLING")==0) {
      int nitems, dum;
      nitems=sscanf(param,"%d", &dum);
      if(nitems==1) {
        criteria->selective_open_sampling=(dum==1);;
        }
      else {
        criteria->selective_open_sampling=(strstr(param,"T")!=0);
        }
      continue;
      }
/*------------------------------------------------------------------------------
    allow/disallow boundary resampling accordingly with criteria */
    if(strcmp(keyword, "#RE_SAMPLE_LOOSE_INDIVIDUALS")==0) {
      int nitems, dum;
      vector<string> tokens;
      tokens=string_split((string) param, " ");
      for(int k=0; k<tokens.size();k++) {
        nitems=sscanf(tokens[k].c_str(),"%d", &dum);
        if(nitems==1) {
          criteria->loose_openlimits.push_back(dum);
          }
        }
      tokens.clear();
      continue;
      }
/*------------------------------------------------------------------------------
    additional smoothing around islands */
    if(strcmp(keyword, "ISLAND_SHELF")==0) {
      int nitems,dum;
      nitems=sscanf(param,"%d", &dum);
      if(nitems==1) {
        criteria->shelf=(dum==1);
        }
      else {
        criteria->shelf=(strstr(param,"T")!=0);
        }
      continue;
      }
/*------------------------------------------------------------------------------
    resolution map smoothing */
    if(strcmp(keyword, "SMOOTH_CRITERION")==0) {
      int nitems,dum;
      nitems=sscanf(param,"%d", &dum);
      if(nitems==1) {
        criteria->smooth=(dum==1);;
        }
      else {
        criteria->smooth=(strstr(param,"T")!=0);
        }
      continue;
      }
/*------------------------------------------------------------------------------
    smoothing tuning*/
    if(strcmp(keyword, "SMOOTH_RATE")==0) {
      sscanf(param,"%lf", &criteria->maxrate);
      continue;
      }
/*------------------------------------------------------------------------------
    smoothing tuning*/
    if(strcmp(keyword, "SMOOTH_NITERATIONS")==0) {
      sscanf(param,"%d", &criteria->niterations);
      continue;
      }
    
/*------------------------------------------------------------------------------
    */
    if(strcmp(keyword, "LIMITS_CONSISTENCY_EXTENT")==0) {
      sscanf(param,"%lf", &criteria->limits_consistency_extent);
      continue;
      }
    if(strcmp(keyword, "LIMITS_CONSISTENCY_FACTOR")==0) {
      sscanf(param,"%lf", &criteria->limits_consistency_factor);
      continue;
      }
    
    if(strcmp(keyword, "INTERIOR_LINES")==0) {
      sscanf(param,"%s", &criteria->interior_lines);
      continue;
      }
    
      
/*------------------------------------------------------------------------------
    input data*/
    if(strcmp(keyword, "TOPOGRAPHY_FILE")==0) {
      sscanf(param,"%s", &criteria->regulardepth);
      continue;
      }
    if(strcmp(keyword, "BOUNDARY_FILE")==0) {
      sscanf(param,"%s", &criteria->boundary);
      continue;
      }
    if(strcmp(keyword, "RANDOM_DEPTH_FILE")==0) {
      sscanf(param,"%s", &criteria->randomdepth);
      continue;
      }
    if(strcmp(keyword, "DELAUNEY_DEPTH_FILE")==0) {
      sscanf(param,"%s", &criteria->delaunay);
      continue;
      }
    if(strcmp(keyword, "MESH_FILE_ROOTNAME")==0) {
      sscanf(param,"%s", &criteria->rootname);
      continue;
      }
    criteria->descriptor[0]='\0';

    if(strcmp(keyword, "DESCRIPTOR_FILE")==0) {
      sscanf(param,"%s", &criteria->descriptor);
//     CALL str_expandpath(criteria->descriptor,rstatus)
//     CALL strg_tocstring(criteria->descriptor)
      continue;
      }
    line[strlen(line)-1]='\0';
    if(verbose>=0) printf("%s: unreckognized keyword (%s) in line %d (%s)\n",__func__,keyword, count, line);
    }
  
  if(criteria->hmin2==0.0) criteria->hmin2=200.0;
  if(criteria->shelf_minsize==0.0) criteria->shelf_minsize=criteria->minsize;
  if(criteria->shelf_maxsize==0.0) criteria->shelf_maxsize=criteria->maxsize;
  
//  criteria.resample=0;
  if(criteria->maxrate==0) criteria->maxrate=0.75;

finished:
  delete[] line;
  fclose(in);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  string fe_criteria_2_formula(const criteria_t & criteria)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  ostringstream o,critList;
  
  /* NOTE: if(h[]<0){h[]=2} in set_density0{1,2} !!! */
  
  /* set_density00 */
  o << "OPEN_SHELF_LIMIT=" << criteria.open_shelf_limit << "\n";
  o << "shelf=h<OPEN_SHELF_LIMIT\n";
  
  o << "OPEN_MAX_RADIUS=" << criteria.maxsize << "\n";
  o << "SHELF_MAX_RADIUS=" << criteria.shelf_maxsize << "\n";
  o << "maxRadius=shelf*SHELF_MAX_RADIUS+(1-shelf)*OPEN_MAX_RADIUS\n";
  critList << "maxRadius";
  
  if(criteria.use(0)){
    /* set_density01 */
    o << "SAMPLE_SCALE_1=" << criteria.factor1 << "\n";
    o << "wlC=1e3*SAMPLE_SCALE_1* h^.5\n";
    critList << ",wlC";
    }
  
  if(criteria.use(1)){
    /* set_density02 */
    o << "OPEN_MIN_DEPTH_2=" << criteria.hmin2 << "\n";
    o << "OPEN_MIN_RATIO=" << criteria.minratio << "\n";
    o << "SAMPLE_SCALE_2=" << criteria.factor2 << "\n";
    o << "slope=gradient_modulus(h)\n";
    o << "slopeApplies=h>OPEN_MIN_DEPTH_2\n";
    o << "slopeC=slopeApplies*SAMPLE_SCALE_2*max(1e3*OPEN_MIN_RATIO,h/slope)+(1-slopeApplies)*maxRadius\n";
    critList << ",slopeC";
    }
  
  o << "roughRadius=min(" << critList.str() << ")\n";
  
//   if(criteria.use(3) || criteria.surface_wave_cfl){
  if(criteria.use(3)){
    o << "# !!! ocean waves not finished\n";
    }
  
  if(criteria.use(4)){
    o << "# !!! lamda4_15 criteria not finished\n";
    }

  if(criteria.smooth==1) {
    o << "SMOOTH_RATE=" << criteria.maxrate << "\n";
    o << "radius=smooth_density(roughRadius,SMOOTH_RATE)\n";
    }
  else{
    o << "radius=roughRadius\n";
    }
  
  return o.str();
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_save_criteria(const char *filename,criteria_t criteria)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in;

  in=fopen(filename,"w");
  
  if(in==0) {
    printf("fe_save_criteria, cannot open %s\n",filename);
    return(-1);
    }


  fprintf(in,"%s\n", "#OPEN_MIN_DEPTH");
  fprintf(in,"%lf\n", (criteria.hmin));

  fprintf(in,"%s\n", "#OPEN_MIN_DEPTH_2");
  fprintf(in,"%lf\n", criteria.hmin2);

  fprintf(in,"%s\n", "#OPEN_MAX_DEPTH");
  fprintf(in,"%lf\n", criteria.hmax);

  fprintf(in,"\n");

  fprintf(in,"%s\n", "#OPEN_MIN_ANGLE");
  fprintf(in,"%lf\n", criteria.minangle);

  fprintf(in,"%s\n", "#OPEN_MIN_RATIO");
  fprintf(in,"%lf\n", criteria.minratio);

  fprintf(in,"\n");

  fprintf(in,"%s\n", "#OPEN_MIN_RADIUS");
  fprintf(in,"%lf\n", criteria.minsize);

  fprintf(in,"%s\n", "#OPEN_MAX_RADIUS");
  fprintf(in,"%lf\n", criteria.maxsize);

  fprintf(in,"%s\n", "#SHELF_MIN_RADIUS");
  fprintf(in,"%lf\n", criteria.shelf_minsize);

  fprintf(in,"%s\n", "#SHELF_MAX_RADIUS");
  fprintf(in,"%lf\n", criteria.shelf_maxsize);

  fprintf(in,"%s\n", "#SHORE_MIN_RADIUS");
  fprintf(in,"%lf\n", criteria.shore_minsize);

  fprintf(in,"%s\n", "#SHORE_MAX_RADIUS");
  fprintf(in,"%lf\n", criteria.shore_maxsize);

  fprintf(in,"%s\n", "#OPEN_SHELF_LIMIT");
  fprintf(in,"%lf\n", criteria.open_shelf_limit);

  fprintf(in,"%s\n", "#SHELF_SHORE_LIMIT");
  fprintf(in,"%lf\n", criteria.shelf_shore_limit);

  fprintf(in,"\n");

  fprintf(in,"%s\n", "#OPEN_DX");
  fprintf(in,"%lf\n", criteria.cellsize);

  fprintf(in,"\n");

  fprintf(in,"%s\n", "#SAMPLE_SCALE_1");
  fprintf(in,"%lf\n", criteria.factor1);

  fprintf(in,"%s\n", "#SAMPLE_SCALE_2");
  fprintf(in,"%lf\n", criteria.factor2);

  fprintf(in,"%s\n", "#SAMPLE_SCALE_3");
  fprintf(in,"%lf\n", criteria.factor3);

  fprintf(in,"%s\n", "#SAMPLE_SCALE_4");
  fprintf(in,"%lf\n", criteria.factor4);

  fprintf(in,"\n");

  fprintf(in,"%s\n", "SURFACE_WAVE_PERIOD");
  fprintf(in,"%lf\n", criteria.surface_wave_period);
  
  fprintf(in,"%s\n", "SURFACE_WAVE_DT");
  fprintf(in,"%lf\n", criteria.surface_wave_dt);
  
  fprintf(in,"\n");

  fprintf(in,"%s\n", "#OPEN_SAMPLE_MODE");
  fprintf(in,"%d\n", criteria.mode);

  fprintf(in,"\n");

  fprintf(in,"%s\n", "#RE_SAMPLE");
  fprintf(in,"%d\n", criteria.resample_obsolete);

  fprintf(in,"%s\n", "#RE_SAMPLE_OPEN");
  fprintf(in,"%d\n", criteria.resample_openlimits);

  fprintf(in,"%s\n", "#RE_SAMPLE_RIGID");
  fprintf(in,"%d\n", criteria.resample_rigidlimits);

  fprintf(in,"%s\n", "#SELECTIVE_OPEN_SAMPLING");
  fprintf(in,"%d\n", criteria.selective_open_sampling);
  
  fprintf(in,"%s\n", "#RE_SAMPLE_LOOSE_INDIVIDUALS");
  for(int k=0; k<criteria.loose_openlimits.size();k++) fprintf(in,"%d", criteria.loose_openlimits[k]);
  fprintf(in,"\n");

  fprintf(in,"\n");

  fprintf(in,"%s\n", "#ISLAND_SHELF");
  fprintf(in,"%d\n", criteria.shelf);

  fprintf(in,"%s\n", "#SMOOTH_CRITERION");
  fprintf(in,"%d\n", criteria.smooth);

  fprintf(in,"%s\n", "#SMOOTH_RATE");
  fprintf(in,"%lf\n", criteria.maxrate);

  fprintf(in,"%s\n", "#SMOOTH_NITERATIONS");
  fprintf(in,"%d\n", criteria.niterations);

  fprintf(in,"\n");

  fprintf(in,"%s\n", "#LIMITS_CONSISTENCY_EXTENT");
  fprintf(in,"%lf\n", criteria.limits_consistency_extent);

  fprintf(in,"%s\n", "#LIMITS_CONSISTENCY_FACTOR");
  fprintf(in,"%lf\n", criteria.limits_consistency_factor);
  
  fprintf(in,"\n");
  
  fprintf(in,"%s\n", "#INTERIOR_LINES");
  fprintf(in,"%s\n", criteria.interior_lines);

  fprintf(in,"\n");

  fprintf(in,"%s\n", "#TOPOGRAPHY_FILE");
  fprintf(in,"%s\n", criteria.regulardepth);

  fprintf(in,"%s\n", "#BOUNDARY_FILE");
  fprintf(in,"%s\n", criteria.boundary);

  fprintf(in,"%s\n", "#RANDOM_DEPTH_FILE");
  fprintf(in,"%s\n", criteria.randomdepth);

  fprintf(in,"%s\n", "#DELAUNEY_DEPTH_FILE");
  fprintf(in,"%s\n", criteria.delaunay);

  fprintf(in,"%s\n", "#MESH_FILE_ROOTNAME");
  fprintf(in,"%s\n", criteria.rootname);

  fprintf(in,"%s\n", "#DESCRIPTOR_FILE");
  fprintf(in,"%s\n", criteria.descriptor);

  fclose(in);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_InitCriteria(const char *specification,criteria_t & criteria)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
   
  criteria.hmin=2.0;
  criteria.hmin2=500.0;
  criteria.hmax=2000.0;
 
  criteria.minratio=5.0;
  criteria.minangle=25.0;
 
  criteria.factor1=5.0;
  criteria.factor2=0.4;
  criteria.factor3=1.0;
  criteria.factor4=0.0;
 
  criteria.open_shelf_limit=450.0;
  
  criteria.mode=3;
  
  criteria.maxrate=0.75;
  criteria.smooth=1;
  criteria.niterations=1000;
  
  criteria.shelf=1;
  
  criteria.surface_wave_period=10.0;
  criteria.surface_wave_dt=10.0;
  
  criteria.limits_consistency_factor=1.5;
  
  strcpy(criteria.descriptor,    "UNSPECIFIED");
  strcpy(criteria.boundary,      "UNSPECIFIED");
  strcpy(criteria.regulardepth,  "UNSPECIFIED");
  strcpy(criteria.randomdepth,   "UNSPECIFIED");
  strcpy(criteria.delaunay,      "UNSPECIFIED");
  strcpy(criteria.rootname,      "UNSPECIFIED");
 
  
  if(strcmp(specification,"COASTAL")==0) {
    criteria.shelf_minsize=1.0;
    }
  else if(strcmp(specification,"REGIONAL")==0) {
    criteria.shelf_minsize=2.5;
    }
  
  criteria.minsize=5*criteria.shelf_minsize;
  
  criteria.shelf_maxsize=5*criteria.shelf_minsize;
  criteria.maxsize=5*criteria.minsize;
  
  criteria.cellsize=criteria.shelf_minsize/2.0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_loadboundaries(char *filename, int format, plg_t **polygones, int *npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,status;

  status=plg_load(filename, format, polygones, npolygones);
  for(k=0;k<*npolygones;k++) {
    status=close_plg(&(*polygones)[k]);
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int fe_mimicgrid_template(grid_t vgrid, mesh_t *mesh, float *v_landmask, float *v_topo, grid_t zgrid, T *z_landmask, const char *filename, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
 
  create a triangular mesh from a structured grid by cutting quadrangles into
  two triangles
  
  some fixes are needed:
  
    * use of tracer landmask instead of vorticity landmask
    
    * changes must be made consistent with use in symmic
 
------------------------------------------------------------------------------*/
{
  int    i,j,k,l,kk,ll,m,n;
  int    count,status;
  double x,y;
  T  z;
  int *incidence;
  int verbose=0;
  T mask=99999;

  incidence=new int[vgrid.nx*vgrid.ny];

  mesh->destroy();

#if 0
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  check for y-periodicity, remove first and last 3 columns
  
  (can't remember what for it was done, surface wave test case?)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  count=0;
  for(i=0;i<vgrid.nx;i++) {
    double ym,yn;
    j=0;
    m=vgrid.nx*j+i;
    j=vgrid.ny-3;
    n=vgrid.nx*j+i;
    x=  vgrid.x[m];
    ym= vgrid.y[m];
    yn= vgrid.y[n];
    if(ym==yn) count++;
    }

  for(i=0;i<vgrid.nx;i++) {
    j=0;
    m=vgrid.nx*j+i;
    v_landmask[m]=0;
    j=1;
    m=vgrid.nx*j+i;
    v_landmask[m]=0;
    j=2;
    m=vgrid.nx*j+i;
    v_landmask[m]=0;
    j=vgrid.ny-1;
    m=vgrid.nx*j+i;
    v_landmask[m]=0;
    j=vgrid.ny-2;
    m=vgrid.nx*j+i;
    v_landmask[m]=0;
    j=vgrid.ny-3;
    m=vgrid.nx*j+i;
    v_landmask[m]=0;
    }
    
#endif

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  minimize band width by numbering along the shortest side
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  mesh->nvtxs=0;
  if(vgrid.nx<vgrid.ny) {
    for(j=0;j<vgrid.ny;j++) {
      for(i=0;i<vgrid.nx;i++) {
        m=vgrid.nx*j+i;
        if(v_landmask[m]>0) {
          incidence[m]=mesh->nvtxs;
          mesh->nvtxs++;
          }
        }
      }
    }
  else {
    for(i=0;i<vgrid.nx;i++) {
      for(j=0;j<vgrid.ny;j++) {
        m=vgrid.nx*j+i;
        if(v_landmask[m]>0) {
          incidence[m]=mesh->nvtxs;
          mesh->nvtxs++;
          }
        }
      }
    }

  mesh->vertices= new vertex_t[mesh->nvtxs];
  mesh->nlayers=vgrid.nz-1;
  mesh->nlevels=vgrid.nz;

  mesh->nnghm=6;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  initialize vertices
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(j=0;j<vgrid.ny;j++) {
    for(i=0;i<vgrid.nx;i++) {
      m=vgrid.nx*j+i;
      if(v_landmask[m]>0.) {
        x= vgrid.x[m];
        y= vgrid.y[m];
        n=incidence[m];
        mesh->vertices[n].lon=x;
        mesh->vertices[n].lat=y;
        mesh->vertices[n].nngh=0;
        mesh->vertices[n].h=v_topo[m];
        mesh->vertices[n].code=0;
        mesh->vertices[n].ngh=new int[mesh->nnghm];
        for(k=0;k<mesh->nnghm;k++) mesh->vertices[n].ngh[k]=-1;
        if(vgrid.z!=0) {
          mesh->vertices[n].sigma  =new double[mesh->nlevels];
          mesh->vertices[n].zlevels=new double[mesh->nlevels];
          for (l=0;l<mesh->nlevels;l++) {
            mesh->vertices[n].zlevels[l]=vgrid.z[m+l*vgrid.ny*vgrid.nx];
            }
          for (l=0;l<mesh->nlevels;l++) {
            mesh->vertices[n].sigma[l]=-vgrid.z[m+l*vgrid.ny*vgrid.nx]/v_topo[m];
            }
          }
        }
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  initialize vertices connections
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(j=0;j<vgrid.ny;j++) {
    for(i=0;i<vgrid.nx;i++) {
      k=vgrid.nx*j+i;
      if(v_landmask[k]==0) continue;
      kk=incidence[k];
/*------------------------------------------------------------------------------
      righthand side neighbour*/
      if(i!=vgrid.nx-1) {
        l=vgrid.nx*j+i+1;
        ll=incidence[l];
        if(v_landmask[l]>0) {
          mesh->vertices[kk].ngh[mesh->vertices[kk].nngh]=ll;
          mesh->vertices[ll].ngh[mesh->vertices[ll].nngh]=kk;
          mesh->vertices[kk].nngh++;
          mesh->vertices[ll].nngh++;
          }
        }
/*------------------------------------------------------------------------------
      upper righthand side neighbour*/
      if((i!=vgrid.nx-1) and (j!=vgrid.ny-1)) {
        l=vgrid.nx*(j+1)+i+1;
        ll=incidence[l];
        if(v_landmask[l]>0.) {
          mesh->vertices[kk].ngh[mesh->vertices[kk].nngh]=ll;
          mesh->vertices[ll].ngh[mesh->vertices[ll].nngh]=kk;
          mesh->vertices[kk].nngh++;
          mesh->vertices[ll].nngh++;
          }
        }
/*------------------------------------------------------------------------------
      roof neighbour*/
      if(j!=vgrid.ny-1) {
        l=vgrid.nx*(j+1)+i;
        ll=incidence[l];
        if(v_landmask[l]>0.) {
          mesh->vertices[kk].ngh[mesh->vertices[kk].nngh]=ll;
          mesh->vertices[ll].ngh[mesh->vertices[ll].nngh]=kk;
          mesh->vertices[kk].nngh++;
          mesh->vertices[ll].nngh++;
          }
        }
      }
    }

/*------------------------------------------------------------------------------
  remove vertices having 1 neighbours or less */
  for(n=0;n<mesh->nvtxs;n++) {
    if(mesh->vertices[n].nngh<2) {
      status=fe_removevertex(mesh,n);
      n=-1;
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  finalize 2D mesh
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_list(mesh);
  if(debug) status=fe_savemesh("mimic-debug-01.nei",MESH_FILE_FORMAT_TRIGRID,*mesh);

  if( (vgrid.nx==zgrid.nx+1) && (vgrid.ny==zgrid.ny+1) ) goto final;
  
/*------------------------------------------------------------------------------
  remove undue connections, aimed for malformed grids (out of structured models)*/
  count=0;
  for(n=0;n<mesh->nvtxs;n++) {
/*-------------------------------------------------------------------------------
    test if edge centre is masked or not in the regular grid model*/
    for(k=0;k<mesh->vertices[n].nngh;k++) {
      m=mesh->vertices[n].ngh[k];
      x=0.5*(mesh->vertices[n].lon+mesh->vertices[m].lon);
      y=0.5*(mesh->vertices[n].lat+mesh->vertices[m].lat);
      status= map_interpolation(zgrid, z_landmask,mask,x,y,&z);
      if(z<0.3) {
        if(verbose==1) printf("disconnect vertices %d %d\n",n,m);
        status= fe_disconnectvertices(*mesh, n,m);
        k=k-1;
        }
      else if(z==mask) {
        if(verbose==1) printf("disconnect vertices %d %d\n",n,m);
        status= fe_disconnectvertices(*mesh, n,m);
        k=k-1;
        }
      else {
        count++;
        }
      }
    }

  status=fe_cleanvertices(mesh, false);
  
  status=fe_list(mesh);
  if(debug) status=fe_savemesh("mimic-debug.nei",MESH_FILE_FORMAT_TRIGRID,*mesh);

  for(n=0;n<mesh->nvtxs;n++) {
    delete[] mesh->vertices[n].ngh;
    }

  status=fe_e2n (mesh);
  status=fe_cleanvertices(mesh, false);
  
  status=fe_list(mesh);
  if(debug) status=fe_savemesh("mimic-debug.nei",MESH_FILE_FORMAT_TRIGRID,*mesh);

final:
  status= fe_connex(*mesh);
  if(status!=0) {
    status= fe_reducebw(*mesh, 1);
    if(status==-2) {
      mesh->destroy();
      status=fe_readmesh("mimic-reduced.nei",MESH_FILE_FORMAT_TRIGRID,mesh);
      }
    }
  status=fe_cleanvertices(mesh, false);
  status=fe_list(mesh);
  if(debug) status=fe_savemesh("mimic-debug-04.nei",MESH_FILE_FORMAT_TRIGRID,*mesh);
  
  status=fe_edgetable(mesh,0,0);
  status=fe_codetable2(mesh,0,1,0);
  status= fe_connex(*mesh);
  if(status!=0) {
    status= fe_reducebw(*mesh, 1);
    if(status==-2) {
      mesh->destroy();
      status=fe_readmesh("mimic-reduced.nei",MESH_FILE_FORMAT_TRIGRID,mesh);
      }
    }
  
  if(filename!=0)
    status=fe_savemesh(filename,MESH_FILE_FORMAT_TRIGRID,*mesh);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_mimicgrid(grid_t vgrid, mesh_t *mesh, float *v_landmask, float *v_topo, grid_t zgrid, float *z_landmask, const char *filename, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=fe_mimicgrid_template(vgrid, mesh, v_landmask, v_topo, zgrid, z_landmask, filename, debug);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_mimicgrid(grid_t vgrid, mesh_t *mesh, char *v_landmask, float *v_topo, grid_t zgrid, char *z_landmask, const char *filename, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  float *z_tmp=new float[zgrid.Hsize()];
  float *v_tmp=new float[vgrid.Hsize()];
  
  for(size_t n=0;n<zgrid.Hsize();n++) z_tmp[n]=z_landmask[n];
  for(size_t n=0;n<vgrid.Hsize();n++) v_tmp[n]=v_landmask[n];
  
  status=fe_mimicgrid(vgrid, mesh, v_tmp, v_topo, zgrid, z_tmp, filename, debug);
  
  delete[] z_tmp;
  delete[] v_tmp;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_Mesh2DElasticityT(mesh_t & mesh, double* & K, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
 *
 * aimed to provide a method for smooth mesh distorsion by computing an elasticity
 * matrix from the original mesh that can be used to recompute vertices positions
 * in the distorded mesh after boundary edges modification
 *
 * here it is coded for triangle meshes, easily adaptable to any other types of
 * meshes (quadrangles,...)
 *
 * Next step will be a displacement operator to fit mesh boundaries with a given
 * line or polygon.
 *
 * Boundary line-wise elasticity to be implemented:
 *
 *   - 1 fixed vertex per boundary
 *   - equilibrium condition for the others: k1 * L1 + k2 * L2 =0
 *
 *------------------------------------------------------------------------------
  
  Compute edge geometrical elasticity
  -----------------------------------
  
  Constrained problem:
  
    n interior vertices equations (to be strictly implemented)

    min potential energy :  E=(sum on neighbours) ki[(xn-xi)²+(yn-yi)²]
    
    It leads to 2 conditions:
           -> dE/dx: (sum on neighbours) ki xn = (sum) ki xi
           -> dE/dy: (sum on neighbours) ki yn = (sum) ki yi
  
    min departure from 1 :  D=(sum on edges) (ki-1)²
  
  
  Unconstrained problem (Lagrange mutipliers):
    
    min J = min E+D =(sum on edges) (ki-1)² + (sum on vertices) l2n   (sum on neighbours j) kj(xn-xj)
                                            + (sum on vertices) l2n+1 (sum on neighbours j) kj(yn-yj),
    where ki and ln are unknowns
    
  
    dJ/dki -> 2(ki-1) + (sum) l2n (xn-xi) (2 vertices n involved) + (sum) l2n+1 (yn-yi) (2 vertices n involved)
  
    dJ/dli -> (sum on neighbours j) kj(x-xj) = 0 (nngh edges involved)
  
    boundary conditions:
  
    ki = 1 for boundary edges
  
    li = 0 for boundary vertices
  
  
  Reference : https://fr.wikipedia.org/wiki/Multiplicateur_de_Lagrange
  
------------------------------------------------------------------------------*/
{
  int e,k,m,n,n1,n2,status;
  int row,col;
  size_t pos;
  hypermatrix_t R, P;
  double *rhs, *weight,Energy,*energyP1;
  double *x,*y;
  int count;
  
  size_t neq=2*mesh.nvtxs+mesh.nedges;
  
  int *vertex=new int[mesh.nvtxs];
  
  double mask=-1;
  
  if(K==0) K=new double[mesh.nedges];
  
  status=discretisation_init(&mesh,LGP1,0);
  
  x=new double[mesh.nvtxs];
  y=new double[mesh.nvtxs];
  
  weight=new double[mesh.nedges];
  for(n=0;n<mesh.nedges; n++) weight[n]=1.0;

  for(n=0;n<mesh.nvtxs; n++) {
    x[n]=mesh.vertices[n].lon;
    y[n]=mesh.vertices[n].lat;
    }
  
  for(n=0;n<mesh.nvtxs; n++) vertex[n]=-1;

redo:
/*------------------------------------------------------------------------------
  Rigidity matrix construction: unknowns are 2*nvertex multipliers (l) + nedges
  elasticity coefficient (k) */
  R.ordering=new ordering_t;
  
  R.ordering->nrows=neq;
  
  R.ordering->cardinal=new int[neq];
  R.ordering->pointer =new int[neq+1];
  
  for(n=0;n<mesh.nvtxs; n++) {
    if(mesh.vertices[n].code!=0) {
/*------------------------------------------------------------------------------
      boundary vertices, clamped condition (l=0) */
      R.ordering->cardinal[2*n]  =1;
      R.ordering->cardinal[2*n+1]=1;
      }
    else {
/*------------------------------------------------------------------------------
      interior vertices, #edges connected to vertex */
      R.ordering->cardinal[2*n]  =mesh.vertices[n].nngh;
      R.ordering->cardinal[2*n+1]=mesh.vertices[n].nngh;
      }
    }
  
  for(n=0;n<mesh.nedges; n++) {
    if(mesh.edges[n].code!=MESH_INTERIOR_EDGE) {
/*------------------------------------------------------------------------------
      boundary edges, clamped condition (k=1) */
      R.ordering->cardinal[n+2*mesh.nvtxs] = 1;
      }
    else {
/*------------------------------------------------------------------------------
      interior edges, 1 edge + 2 x 2 vertices */
      R.ordering->cardinal[n+2*mesh.nvtxs] = 5;
      }
    }
  
  R.ordering->pointer[0]=0;
  for(n=0;n<neq; n++) {
    R.ordering->pointer[n+1]=R.ordering->pointer[n]+R.ordering->cardinal[n];
    }
  
  size_t nnz=R.ordering->pointer[neq];
  R.ordering->incidence=new int[nnz];
  
/*------------------------------------------------------------------------------
  multipliers lines */
  for(n=0;n<mesh.nvtxs; n++) {
    if(mesh.vertices[n].code!=0) {
/*------------------------------------------------------------------------------
      boundary vertices, self-connection */
      pos=R.ordering->pointer[2*n];
      R.ordering->incidence[pos]=2*n;
      pos=R.ordering->pointer[2*n+1];
      R.ordering->incidence[pos]=2*n+1;
      }
    else {
/*------------------------------------------------------------------------------
      interior vertices, edges connected to vertex */
      int count=0;
      for(k=0;k<mesh.vertices[n].nngh;k++) {
        m=mesh.vertices[n].ngh[k];
        e=fe_isedge(mesh,n,m);
        pos=R.ordering->pointer[2*n];
        R.ordering->incidence[pos+count]=e+2*mesh.nvtxs;
        pos=R.ordering->pointer[2*n+1];
        R.ordering->incidence[pos+count]=e+2*mesh.nvtxs;
        count++;
        }
      }
    }
  
/*------------------------------------------------------------------------------
  elastic coefficients lines */
  for(n=0;n<mesh.nedges; n++) {
    row=n+2*mesh.nvtxs;
    if(mesh.edges[n].code!=MESH_INTERIOR_EDGE) {
/*------------------------------------------------------------------------------
      boundary edges, self-connection */
      pos=R.ordering->pointer[row];
      R.ordering->incidence[pos] = row;
      }
    else {
/*------------------------------------------------------------------------------
      interior edges, self-connection + 2 vertices x 2 (->x,y) */
      pos=R.ordering->pointer[row];
      R.ordering->incidence[pos] = row;
      R.ordering->incidence[pos+1]=2*mesh.edges[n].extremity[0];
      R.ordering->incidence[pos+2]=2*mesh.edges[n].extremity[1];
      R.ordering->incidence[pos+3]=2*mesh.edges[n].extremity[0]+1;
      R.ordering->incidence[pos+4]=2*mesh.edges[n].extremity[1]+1;
      }
    }
  
  status=matrix_reorder(R.ordering);
  
  R.allocate();
  rhs=new double[neq];
  
  R.assign(0.0);
  
  for(n=0;n<mesh.nvtxs; n++) {
    if(mesh.vertices[n].code!=0) {
      row=2*n;
      pos=R.ordering->pos(row,row);
      R.packed[pos] =1.0;
      rhs[row]=0;
      row=2*n+1;
      pos=R.ordering->pos(row,row);
      R.packed[pos] =1.0;
      rhs[row]=0;
      }
    else {
      for(k=0;k<mesh.vertices[n].nngh;k++) {
        m=mesh.vertices[n].ngh[k];
        e=fe_isedge(mesh,n,m);
        col=e+2*mesh.nvtxs;
        double dx=x[n]-x[m];
        row=2*n;
        pos=R.ordering->pos(row,col);
        R.packed[pos] =dx;
        rhs[row]=0.0;
        double dy=y[n]-y[m];
        row=2*n+1;
        pos=R.ordering->pos(row,col);
        R.packed[pos] =dy;
        rhs[row]=0.0;
        }
      }
    }
  
  for(n=0;n<mesh.nedges; n++) {
    row=n+2*mesh.nvtxs;
    if(mesh.edges[n].code!=MESH_INTERIOR_EDGE) {
/*------------------------------------------------------------------------------
      boundary edges */
      pos=R.ordering->pos(row,row);
      R.packed[pos] =1.0;
      rhs[row]=1.0;
      }
    else {
/*------------------------------------------------------------------------------
      potential energy minimisation */
      n1=mesh.edges[n].extremity[0];
      n2=mesh.edges[n].extremity[1];
      double dx=x[n2]-x[n1];
      col=2*mesh.edges[n].extremity[0];
      pos=R.ordering->pos(row,col);
      R.packed[pos] = -dx;
      col=2*mesh.edges[n].extremity[1];
      pos=R.ordering->pos(row,col);
      R.packed[pos] = +dx;
      double dy=y[n2]-y[n1];
      col=2*mesh.edges[n].extremity[0]+1;
      pos=R.ordering->pos(row,col);
      R.packed[pos] = -dy;
      col=2*mesh.edges[n].extremity[1]+1;
      pos=R.ordering->pos(row,col);
      R.packed[pos] = +dy;
      pos=R.ordering->pos(row,row);
/*------------------------------------------------------------------------------
      departure minimisation */
      R.packed[pos]=weight[n]*2.0;
      rhs[row]=weight[n]*2.0;
      }
    }
  
/*------------------------------------------------------------------------------
  solve rigidity problem */
  R.distributor=new distribution_t;
  R.distributor->context=SEQUENTIAL_COMPUTING;
      
  bool force_duplicate=false;
  int verbose=0;
  status=LinearSystem_initialise(&R, 0, 0, SOLVER_ID_UMFPACK, force_duplicate, verbose, debug);
  status=LinearSystem_solve(R, rhs);
  
/*------------------------------------------------------------------------------
  enforce positive rigidity by increasing weights in regularisation term */
  count=0;
  for(n=0;n<mesh.nedges; n++) {
    row=n+2*mesh.nvtxs;
    K[n]=rhs[row];
    if(K[n]<0.0) {
//       printf("edge %d: K=%lf\n",K[n]);
      weight[n]*=2.0;
      count++;
      }
    else if(K[n]<0.1) {
//       printf("edge %d: K=%lf\n",K[n]);
      weight[n]+=2*(1-K[n]);
      count++;
      }
    }
  
  
  status= archiving_UGdummy2D("test.nc", mesh, "K", "N/m", K, mask, NCP1);
  status= archiving_UGdummy2D("test.nc", mesh, "weight", "none", weight, mask, NCP1);
  
  energyP1=new double[mesh.nvtxs];
  
  for(n=0;n<mesh.nvtxs; n++) {
    energyP1[n]=0;
   if(mesh.vertices[n].code!=0) continue;
   for(k=0;k<mesh.vertices[n].nngh;k++) {
      m=mesh.vertices[n].ngh[k];
//       if(mesh.vertices[n].code!=0) {
      double dx=x[m]-x[n];
      double dy=y[m]-y[n];
      e=fe_isedge(mesh,n,m);
      energyP1[n]+=0.5*K[e]*(dx*dx+dy*dy)*1.e+10;
//       }
//     else {
//       
      }
    energyP1[n]/=mesh.vertices[n].mw;
    }
  status= archiving_UGdummy2D("test.nc", mesh, "energy", "Nm", energyP1, mask, LGP1);
  
  Energy=0;
  for(n=0;n<mesh.nedges; n++) {
    n1=mesh.edges[n].extremity[0];
    n2=mesh.edges[n].extremity[1];
    double dx=x[n2]-x[n1];
    double dy=y[n2]-y[n1];
    Energy+=K[n]*(dx*dx+dy*dy);
    }

  R.destroy();
  delete[] rhs;
  
  printf("#negative K edges=%d energy=%lf\n",count,Energy);
  
/*------------------------------------------------------------------------------
  if positive rigidity not reached, try again */
  if(count!=0) goto redo;
  
/*------------------------------------------------------------------------------
  now check systm by computing positions (nothing should move as boundary nodes kept unmoved) */
  neq=mesh.nvtxs;
  
  P.ordering=new ordering_t;
  
  P.ordering->nrows=neq;
  
  P.ordering->cardinal=new int[neq];
  P.ordering->pointer =new int[neq+1];
  
  for(n=0;n<mesh.nvtxs; n++) {
    if(mesh.vertices[n].code!=0) {
      P.ordering->cardinal[n]  =1;
      }
    else {
      P.ordering->cardinal[n]  =1+mesh.vertices[n].nngh;
      }
    }
  
  P.ordering->pointer[0]=0;
  for(n=0;n<neq; n++) {
    P.ordering->pointer[n+1]=P.ordering->pointer[n]+P.ordering->cardinal[n];
    }
  
  nnz=P.ordering->pointer[neq];
  P.ordering->incidence=new int[nnz];
  
  for(n=0;n<mesh.nvtxs; n++) {
    if(mesh.vertices[n].code!=0) {
      pos=P.ordering->pointer[n];
      P.ordering->incidence[pos]=n;
      }
    else {
      int count=0;
      pos=P.ordering->pointer[n];
      P.ordering->incidence[pos+count]=n;
      count++;
      for(k=0;k<mesh.vertices[n].nngh;k++) {
        m=mesh.vertices[n].ngh[k];
        P.ordering->incidence[pos+count]=m;
        count++;
        }
      }
    }
  
  status=matrix_reorder(P.ordering);
  
  P.allocate();
  rhs=new double[neq];
  
  P.assign(0.0);
  
//   x[0]-=0.02;
  
  for(n=0;n<mesh.nvtxs; n++) {
    if(mesh.vertices[n].code!=0) {
      row=n;
      pos=P.ordering->pos(row,row);
      P.packed[pos] =1.0;
      rhs[row]=x[n];
      }
    else {
      for(k=0;k<mesh.vertices[n].nngh;k++) {
        int m=mesh.vertices[n].ngh[k];
        int e=fe_isedge(mesh,n,m);
//         double dx=x[n]-x[m];
        row=n;
        pos=P.ordering->pos(row,row);
        P.packed[pos] += K[e];
        col=m;
        pos=P.ordering->pos(row,col);
        P.packed[pos] = -K[e];
        rhs[row]=0.0;
        }
      }
    }
  
/*------------------------------------------------------------------------------
  solve position problem */
  P.distributor=new distribution_t;
  P.distributor->context=SEQUENTIAL_COMPUTING;
      
  status=LinearSystem_initialise(&P, 0, 0, SOLVER_ID_UMFPACK, force_duplicate, verbose, debug);
  status=LinearSystem_solve(P, rhs);
  
/*------------------------------------------------------------------------------
  solve position problem */
  count=0;
  for(n=0;n<mesh.nvtxs; n++) {
    mesh.vertices[n].lon=rhs[n];
    double d=fabs(x[n]-rhs[n]);
    if(d>1.e-06) {
//       printf("edge %d: K=%lf\n",K[n]);
      count++;
      }
    }

  P.destroy();
  
  delete[] rhs;
  delete[] weight;
  
  delete[] x;
  delete[] y;
  
  return(0);
}
