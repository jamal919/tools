/*******************************************************************************

  T-UGO tools, 2006-2015

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief finite elements geometry edition functions
*/
/*----------------------------------------------------------------------------*/


#include <stdio.h>
#include <string.h>

#include "tools-structures.h"

#include "fe.def"
#include "fe.h"

#include "map.h"
#include "geo.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_ReallocateVertices(mesh_t *mesh, int nadd)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  re-allocate vertices array */
{
  int    first;
  vertex_t *dum=NULL;

  exitIfNull(
    dum=new vertex_t[mesh->nvtxs+nadd]
    );

  bcopy(mesh->vertices,dum,(mesh->nvtxs)*sizeof(vertex_t));
//  for(n=0;n<mesh->nvtxs;n++) mesh->vertices[n].destroy();
  delete[] mesh->vertices; /// HERE !!!
  mesh->vertices=dum;

  first=mesh->nvtxs;
  mesh->nvtxs+=nadd;

  return(first);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_setvertex(mesh_t *mesh, double lon, double lat, float z, int *ngh, int nngh, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  add vertex */
{
  int    k,n;

  n=target;
  
  updatemax(&mesh->nnghm,nngh);

  mesh->vertices[n].ngh=new int[nngh];

  mesh->vertices[n].nngh=nngh;
  for(k=0;k<nngh;k++) {
    mesh->vertices[n].ngh[k]=ngh[k];
    }

  mesh->vertices[n].lon=lon;
  mesh->vertices[n].lat=lat;
  mesh->vertices[n].h=z;

  return(n);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_addvertex(mesh_t *mesh, double lon, double lat, float z, int *ngh, int nngh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  add vertex

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int    k,n;

/*------------------------------------------------------------------------------
  resize vertices array */
  n=fe_ReallocateVertices(mesh, 1);
  
  mesh->vertices[n].nngh=nngh;
  
  if(nngh>0){
    updatemax(&mesh->nnghm,nngh);
    
    mesh->vertices[n].ngh=new int[nngh];
    
    for(k=0;k<nngh;k++) {
      mesh->vertices[n].ngh[k]=ngh[k];
      }
    }

  mesh->vertices[n].lon=lon;
  mesh->vertices[n].lat=lat;
  mesh->vertices[n].h=z;

  return(n);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_addvertex3D(mesh_t *mesh, double lon, double lat, double h, double *zlevel,int *ngh, int nngh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    l, n;

  n=fe_addvertex(mesh, lon, lat, (float) h, ngh, nngh);

  if(zlevel==0) goto finished;

  mesh->vertices[n].zlevels=new double[mesh->nlevels];
  mesh->vertices[n].sigma=new double[mesh->nlevels];
  for (l=0;l<mesh->nlevels;l++) {
    mesh->vertices[n].zlevels[l]=zlevel[l];
    mesh->vertices[n].sigma[l]=zlevel[l]/zlevel[0];
    }
finished:
  return(n);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_insertvertex(mesh_t *mesh, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Purpose : insert n in its neighbours' neighbour list (cross information 
            completion)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k,l,m;
  int *dum=NULL;

  for(k=0;k<mesh->vertices[n].nngh;k++) {
    m=mesh->vertices[n].ngh[k];
    if(m==-1) {
      printf("troubles at node %d, has -1  as %dth neighbour\n",m,k+1);
      }
    l=fe_anb(*mesh, n, m);
    if(l==-1) {
      mesh->vertices[m].nngh++;
      dum= new int[mesh->vertices[m].nngh];
      for(l=0;l<mesh->vertices[m].nngh-1;l++) dum[l]=mesh->vertices[m].ngh[l];
//      bcopy(mesh->vertices,dum,(mesh->vertices[m].nngh-1)*sizeof(int));
      delete[] mesh->vertices[m].ngh;
      mesh->vertices[m].ngh=dum;
      mesh->vertices[m].ngh[mesh->vertices[m].nngh-1]=n;
//      printf("add %d in %d list\n",n,m);
      }
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_disconnectvertex(mesh_t & mesh, int n1, int n2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Purpose : remove vertex n1 from n2's neighbour list

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k,l;

  k=fe_anb(mesh, n1, n2);
  if(k!=-1) {
    mesh.vertices[n2].nngh--;
    for(l=k;l<mesh.vertices[n2].nngh;l++) {
      mesh.vertices[n2].ngh[l]=mesh.vertices[n2].ngh[l+1];
      }
    mesh.vertices[n2].ngh[mesh.vertices[n2].nngh]=-1;
    return(0);
    }
  else
    return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_disconnectvertices(mesh_t & mesh, int n1, int n2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Purpose : remove vertex n1 in neighbours list of n2 and reverse

  it can be used to delete a mesh line

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int pos1,pos2,l;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  pos1 = position of n1 in neighbours list of n2
  pos2 = position of n2 in neighbours list of n1
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  pos1=fe_anb(mesh, n1, n2);
  pos2=fe_anb(mesh, n2, n1);

  if((pos1==-1) or (pos2==-1)) {
    printf("inconsistent neighbours: %d %d \n",n1,n2);
    return(-1);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  disconnect vertices
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  mesh.vertices[n2].nngh--;
  for(l=pos1;l<mesh.vertices[n2].nngh;l++) {
    mesh.vertices[n2].ngh[l]=mesh.vertices[n2].ngh[l+1];
    }
  mesh.vertices[n2].ngh[mesh.vertices[n2].nngh]=-1;

  mesh.vertices[n1].nngh--;
  for(l=pos2;l<mesh.vertices[n1].nngh;l++) {
    mesh.vertices[n1].ngh[l]=mesh.vertices[n1].ngh[l+1];
    }
  mesh.vertices[n1].ngh[mesh.vertices[n1].nngh]=-1;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  check mesh.nnghm. Pretty unsafe (and bugged) version until 05/03/2018
  
  replaced by a tedious but safer check (could be optimized ?)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   if(max(mesh.vertices[n2].nngh,mesh.vertices[n1].nngh)+1==mesh.nnghm){
//     int n;
//     for(n=0;n<mesh.nvtxs;n++){
//       if(mesh.vertices[n].nngh==mesh.nnghm)
//         break;
//       }
//     if(n==mesh.nvtxs)
// /*------------------------------------------------------------------------------
//       05/03/2018 : bug fix */
// //       mesh.nvtxs--;
//       mesh.nnghm--;
//     }

  mesh.nnghm=0;
  for(int n=0;n<mesh.nvtxs;n++) {
    mesh.nnghm=max(mesh.nnghm, mesh.vertices[n].nngh);
    }
    
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_substitute_vertex(mesh_t & mesh, int n1, int n2, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Purpose : substitute n1 with target in n2 neigbour list, and symetrically

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int pos1,pos2;

  pos1=fe_anb(mesh, n1, n2);
  pos2=fe_anb(mesh, n2, n1);

  if((pos1==-1) || (pos2==-1)) {
    printf("inconsistent neighbours: %d %d \n",n1,n2);
    return(-1);
    }
  mesh.vertices[n2].ngh[pos1]=target;
  mesh.vertices[n1].ngh[pos2]=target;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_disconnectvertex(mesh_t & mesh, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Purpose : fully disconnect vertex n

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int l,m,status;
  int nngh,*ngh;
  
  nngh=mesh.vertices[n].nngh;
  ngh=new int[nngh];
  
  for(l=0;l<mesh.vertices[n].nngh;l++) {
    ngh[l]=mesh.vertices[n].ngh[l];
    }

  for(l=0;l<nngh;l++) {
    m=ngh[l];
    status=fe_disconnectvertices(mesh, n, m);
    if(status==-1) {
      printf("%s : anomaly at node %d\n",__func__,n);
      }
    }
  
  delete[] ngh;
  
  if(mesh.vertices[n].nngh==0) {
    return(0);
    }
  else {
    return(-1);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_connectvertices(mesh_t & mesh, int n1, int n2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  add vertex n1 in n2 neighbour list and reprocically 
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int pos1,pos2,i,l,n;
  int *tmp;
  vertex_t *vertex;
  
  if(n1==-1 or n2==-1) return(-1);

  pos1=fe_anb(mesh, n1, n2);
  pos2=fe_anb(mesh, n2, n1);

  if((pos1!=-1) and (pos2!=-1)) return(0);
  if((pos1==-1) and (pos2!=-1)) return(-1);
  if((pos1!=-1) and (pos2==-1)) return(-1);
  
  for(i=0;i<2;i++) {
    if(i==0){
      n=n2;
      l=n1;
      }
    else{
      n=n1;
      l=n2;
      }
    vertex=&mesh.vertices[n]; 
    updatemax(&vertex->nngh,0);
    vertex->nngh++;
    updatemax(&mesh.nnghm,vertex->nngh);
    tmp=new int[vertex->nngh];
    if(vertex->nngh>1) valcpy(tmp,vertex->ngh,vertex->nngh-1);
    tmp[vertex->nngh-1]=l;
    if(vertex->ngh!=0) delete[] vertex->ngh;
    vertex->ngh=tmp;
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_removevertex(mesh_t *mesh, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  Remove target from vertices description

  Warning : element description is no more adequate !!!

-----------------------------------------------------------------------------*/
{
  int k,l,m;

/*------------------------------------------------------------------------------
  remove target from neighbours list*/
  for (m=0;m<mesh->nvtxs;m++) {
    k=fe_anb(*mesh, target, m);
    if(k!=-1) {
      vertex_t *vertex=&mesh->vertices[m];
      vertex->nngh--;
      for(l=k;l<vertex->nngh;l++) {
        vertex->ngh[l]=vertex->ngh[l+1];
        }
      vertex->ngh[vertex->nngh]=-1;
      }
    }

/*------------------------------------------------------------------------------
  remove target from vertex array*/
  mesh->nvtxs--;

  for (m=target;m<mesh->nvtxs;m++) {
    mesh->vertices[m]=mesh->vertices[m+1];
    }

/*------------------------------------------------------------------------------
  re-number entries in neighbours list*/
  mesh->nnghm=0;
  for (m=0;m<mesh->nvtxs;m++) {
    vertex_t *vertex=&mesh->vertices[m];
    updatemax(&mesh->nnghm,vertex->nngh);
    for(k=0;k<vertex->nngh;k++) {
      if(vertex->ngh[k]>target) {
        vertex->ngh[k]--;
        }
      }
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_remove_vertices(mesh_t & mesh, int *target, int ntargets)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  Remove target from vertices description

  Warning : element description is no more adequate !!!

-----------------------------------------------------------------------------*/
{
//  int nvtxs=mesh.nvtxs-ntargets;
  mesh_t work;

/*------------------------------------------------------------------------------
  remove target from neighbours list*/

//   for (m=0;m<mesh->nvtxs;m++) {
//     k=fe_anb(*mesh, target, m);
//     if(k!=-1) {
//       mesh->vertices[m].nngh--;
//       for(l=k;l<mesh->vertices[m].nngh;l++) {
//         mesh->vertices[m].ngh[l]=mesh->vertices[m].ngh[l+1];
//         }
//       mesh->vertices[m].ngh[mesh->vertices[m].nngh]=-1;
//       }
//     }
//
// /*------------------------------------------------------------------------------
//   remove target from vertex array*/
//   mesh->nvtxs--;
//
//   for (m=target;m<mesh->nvtxs;m++) {
//     mesh->vertices[m]=mesh->vertices[m+1];
//     }
//
// /*------------------------------------------------------------------------------
//   re-number entries in neighbours list*/
//   for (m=0;m<mesh->nvtxs;m++) {
//     for(k=0;k<mesh->vertices[m].nngh;k++) {
//       if(mesh->vertices[m].ngh[k]>target) {
//         mesh->vertices[m].ngh[k]--;
//         }
//       }
//     }
//
//   return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_displacevertex(mesh_t & mesh, int target, double lon, double lat)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Displace target

  Warning : element description (geometry) is no more adequate !!!

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k,l,n1,n2,n3,rotation=1,status;
  double angle;
  double lon_bkp,lat_bkp;
  double epsilon=1.e-10;
  
  vertex_t *vertex=&mesh.vertices[target];

  lon_bkp=vertex->lon;
  lat_bkp=vertex->lat;
  
  if(vertex->ngh!=0)
    status=fe_order2(mesh,target,vertex->ngh[0],rotation);
  
  vertex->lon=degree_recale(lon,0.);
  vertex->lat=lat;
  
  /* vertex is on boundary */
  if(vertex->code!=0)
    return 0;
  
  /* vertex is not of a trianlge */
  if(vertex->nelmts<=0)
    return 0;
  
/*-----------------------------------------------------------------------------
  check elements good shape */
  for(k=0;k<vertex->nngh;k++) {
    n1=target;
    n2=vertex->ngh[k];
    if(k==vertex->nngh-1) l=0;
    else l=k+1;
    n3=vertex->ngh[l];
//    angle=fe_angle_cartesian(mesh, n1,n2,n3);
    angle=fe_angle(mesh, n1,n2,n3);
/*-----------------------------------------------------------------------------
    negative angle means unsafe move for node n; if occures, cancel it*/
    if(angle  < epsilon) {
      vertex->lon=lon_bkp;
      vertex->lat=lat_bkp;
      return(-1);
      }
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_cleanvertices(mesh_t *mesh, vector<int> & watch, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Remove unused vertices from mesh description

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k,m,n,status;
  int *ancestor=NULL,*successor=NULL;
  int count=0;
  vertex_t *dum=NULL;

  int previous=mesh->nvtxs;
  
  ancestor =new int[mesh->nvtxs];
  successor=new int[mesh->nvtxs];
  for (m=0;m<mesh->nvtxs;m++) successor[m]=-1;

  for (m=0;m<mesh->nvtxs;m++) {
    k=vpos(-1,mesh->vertices[m].ngh,mesh->vertices[m].nngh);
    if(k!=-1) {
      TRAP_ERR_EXIT(-1,"");perror("fe_cleanvertices");
      }
//    if(mesh->vertices[m].nngh>0) {
    if(mesh->vertices[m].nngh>1) {
      ancestor[count]=m;
      successor[m]=count;
      count++;
      }
    else if(mesh->vertices[m].nngh==1) {
      if(debug) printf("node %d has %d neighbours\n",m,mesh->vertices[m].nngh);
      status=fe_disconnectvertices(*mesh, mesh->vertices[m].ngh[0],m);
      for (m=0;m<mesh->nvtxs;m++) successor[m]=-1;
      count=0;
      m=-1;
      }
    else if(mesh->vertices[m].nngh==0) {
      if(debug) printf("node %d has %d neighbours\n",m,mesh->vertices[m].nngh);
      }
    }

  for(int k=0;k<watch.size();k++) {
    n=watch[k];
    printf("node %d successor %d\n",n,successor[n]);
    }
  
  for(n=0;n<mesh->nvtxs;n++) {
    if(successor[n]==-1) mesh->vertices[n].destroy();
    }
  
  for (m=0;m<count;m++) {
    mesh->vertices[m]=mesh->vertices[ancestor[m]];
    }

  for (m=0;m<count;m++) {
    for(k=0;k<mesh->vertices[m].nngh;k++) {
      n=mesh->vertices[m].ngh[k];
      mesh->vertices[m].ngh[k]=successor[n];
      }
    }

  exitIfNull(
    dum=new vertex_t[count]
    );
  /// HERE !!!
  bcopy(mesh->vertices,dum,(count)*sizeof(vertex_t));
  
  delete[] mesh->vertices;
  mesh->vertices=dum;

  mesh->nvtxs=count;
  
  if(previous!=mesh->nvtxs) printf("#fe_cleanvertices: %d in entry, %d after processing\n",previous, mesh->nvtxs);
  
  if(previous!=mesh->nvtxs) {
    if(mesh->triangles!=0) mesh->triangles->destroy();
    if(mesh->edges!=0)     mesh->edges->destroy();
    mesh->triangles=0;
    mesh->edges=0;
    }
  
  delete[] ancestor;
  delete[] successor;

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_cleanvertices(mesh_t *mesh, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Remove unused vertices from mesh description

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  vector<int> watch;
  
  status=fe_cleanvertices(mesh, watch, debug);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_cleanisolated(mesh_t & mesh, int maxpercent)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,nn,start,first,low;
  int nmax;
  mesh_t work;
  int *idx=NULL,*tmp=NULL,*old=NULL;
  int connex=1;

  if(maxpercent==100) {
    nmax=mesh.nvtxs;
    }
  else {
    nmax=mesh.nvtxs*(maxpercent/100.);
    }
  
  idx      =new int[mesh.nvtxs];
  tmp      =new int[mesh.nvtxs];
  old      =new int[mesh.nvtxs];
/*------------------------------------------------------------------------------
  */
  start=-1;
  while(start<nmax) {
//next:
    start++;
    if (start==mesh.nvtxs) break;
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
          updatemin(&low,nn);
          }
        }
      }
    if(work.nvtxs!=mesh.nvtxs) {
//      printf("fe_cleanisolated: mesh is not connex\n");
      check_error(-3,"\nWARNING : mesh is not connex\n", __LINE__, __FILE__, 0);
      connex=0;
      }
    for(n = 0; n < mesh.nvtxs; n++) idx[n]=tmp[n];
    if(connex==1) break;
    if(connex==0) break;
    }

  if(connex==1) work.nvtxs=mesh.nvtxs;
  
  work.nnghm=mesh.nnghm;
  work.vertices=new vertex_t[work.nvtxs];
  for (n=0; n<work.nvtxs; n++) work.vertices[n].null_value();

  m=0;
  for(n=0;n<mesh.nvtxs;n++) {
    int m=idx[n];
    if(m==-1) continue;
    work.vertices[m].lon =mesh.vertices[n].lon;
    work.vertices[m].lat =mesh.vertices[n].lat;
    work.vertices[m].h   =mesh.vertices[n].h;
    work.vertices[m].code=0;
    work.vertices[m].nngh=mesh.vertices[n].nngh;
    work.vertices[m].ngh=new int[work.vertices[m].nngh];
    for(k=0;k<work.vertices[m].nngh;k++) {
      work.vertices[m].ngh[k]=idx[mesh.vertices[n].ngh[k]];
      }
    }

  delete[] idx;
  delete[] tmp;
  delete[] old;

  if(connex!=1) {
    printf("#fe_cleanisolated: #initial=%d #final=%d vertices (#isolated=%d)\n\n", mesh.nvtxs, work.nvtxs, mesh.nvtxs-work.nvtxs);
    }
  else {
    printf("#fe_cleanisolated, no action to be taken\n");
    }
  
  work.type=mesh.type;
  
  mesh.destroy();
  
  mesh=work;
  
  if(connex==1) return(0);
  else return(-2);

//error:
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_is_neighbour(mesh_t *mesh, int n , int m)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;

  for(i=0;i<mesh->vertices[n].nngh;i++)
    if(mesh->vertices[n].ngh[i]==m) return(0);

  return(-99);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_mergevertex(mesh_t *mesh, int n1 , int n2, int remove)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

 if n1 is a neighbour of n2 --> then n1 gives all his neighbours to n2 and
 is removed !!!

------------------------------------------------------------------------------*/
{
  int i,k,n,status;
  int count;
  int skip;
  int *neighbours=NULL;

/*------------------------------------------------------------------------------
  check if merge will not create mal-formed triangles*/
  status=fe_displacevertex(*mesh, n1, mesh->vertices[n2].lon, mesh->vertices[n2].lat);
  if(status!=0) return(status);

/*------------------------------------------------------------------------------
  copy n2 neighbour list*/
  neighbours=new int[mesh->vertices[n1].nngh+mesh->vertices[n2].nngh];

  count=0;
  for(i=0;i<mesh->vertices[n2].nngh;i++) {
    if(mesh->vertices[n2].ngh[i]==n1) continue;
    neighbours[count]=mesh->vertices[n2].ngh[i];
    count++;
    }

/*------------------------------------------------------------------------------
  substitute n1 with n2 in n1 neighbours' neighbour list*/
  for(i=0;i<mesh->vertices[n1].nngh;i++) {
    n=mesh->vertices[n1].ngh[i];
    if(n==n2) continue;
/*------------------------------------------------------------------------------
    remove n2 from n neighbour list ??? : thought to be wrong */
//     for(k=0;k<mesh->vertices[n].nngh;k++) {
//       status=fe_disconnectvertex(*mesh,n2,n);
//       }
/*------------------------------------------------------------------------------
    remove n2 from n neighbour list (to avoid to have it twice if common neighbour) */
    status=fe_disconnectvertex(*mesh,n2,n);
    for(k=0;k<mesh->vertices[n].nngh;k++) {
      if(mesh->vertices[n].ngh[k]==n1) {
        mesh->vertices[n].ngh[k]=n2;
        break;
        }
      }
    }

/*------------------------------------------------------------------------------
  rebuild n2 neighbour list*/
  for(i=0;i<mesh->vertices[n1].nngh;i++) {
    n=mesh->vertices[n1].ngh[i];
    if(n==n2) continue;
    skip=0;
/*------------------------------------------------------------------------------
    check if n already in neighbour list*/
    for(k=0;k<count;k++) {
      if(neighbours[k]==n) {
        skip=1;
        break;
        }
      }
/*------------------------------------------------------------------------------
    if n not already in neighbour list, add it*/
    if(skip==0) {
      neighbours[count]=n;
      count++;
      }
    }

  mesh->vertices[n1].nngh=0;
  mesh->vertices[n2].nngh=count;
  delete[] mesh->vertices[n2].ngh;
  mesh->vertices[n2].ngh=new int[mesh->vertices[n2].nngh];
  for(i=0;i<mesh->vertices[n2].nngh;i++) mesh->vertices[n2].ngh[i]=neighbours[i];

  if(remove==1) status=fe_removevertex(mesh,n1);

  delete[] neighbours;
  return(0);
}
