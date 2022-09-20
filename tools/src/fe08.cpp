
/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/


#include <stdio.h>
#include <string.h>


#include "tools-structures.h"

#include "geo.h"
#include "fe.h"
#include "matrix.h"
#include "zapper.h"     /*  rutin.h contains common utility routines  */


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int update_edge_table(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l, i, j, j1, j2, jm1, jp1, m, n1, n2, n, status;
  int nedges, chk = 0, count = 0, tmp = 0;
  int nelts, nndes;
  int interior = 0, boundary = 0, weird = 0;
  int ring[4] = { 0, 1, 2, 0 };
  int *ncells, ncellsmax = 0, **cells;
  edge_t    *edges;
  triangle_t *triangles, *ptr;

  nedges = 0;

  triangles = mesh->triangles;
  edges     = mesh->edges;

  nelts  = mesh->ntriangles;
  nndes  = mesh->nvtxs;
  nedges = mesh->nedges;

  for(n = 0; n < nedges; n++) {
    mesh->edges[n].nshared = 0;
    }

  for(m = 0; m < mesh->ntriangles; m++) {
    for(j = 0; j < 3; j++) {
      n = mesh->triangles[m].edges[j];
      if(edges[n].nshared > 1) {
        printf("found %d as weird edge, nshared=%d\n", count, edges[count].nshared);
        printf("elements: %d %d\n", edges[count].shared[0],edges[count].shared[1]);
        printf("vertices: %d %d\n", edges[count].extremity[0],edges[count].extremity[1]);
        }
      edges[n].shared[edges[n].nshared] = m;
      edges[n].nshared++;
      }
    }

  boundary = 0;
  interior = 0;
  for(n = 0; n < nedges; n++) {
    if(edges[n].nshared == 0) {
      weird++;
      }
    if(edges[n].nshared > 2) {
      weird++;
      printf("found %d as weird edge, nshared=%d\n", n, edges[n].nshared);
      printf("elements: %d %d %d\n", edges[n].shared[0],edges[n].shared[1],edges[n].shared[2]);
      printf("vertices: %d %d\n", edges[n].extremity[0],edges[n].extremity[1]);
      }
    if(edges[n].nshared == 1)
      boundary++;
    if(edges[n].nshared == 2)
      interior++;
    }

  if(weird!=0) {
    STDOUT_BASE_LINE_FUNC("\nfound %d total edges\n", nedges);
    report_edge_table( interior, boundary, weird,(interior * 2 + boundary) / 3, nelts);
    TRAP_ERR_EXIT(-1, "insane mesh, abort");
    }

/*------------------------------------------------------------------------------
  scan boundary edges for outward normal determination*/
  for(n = 0; n < mesh->nedges; n++) {
    if(mesh->edges[n].nshared == 1) {
      n1 = mesh->edges[n].extremity[0];
      n2 = mesh->edges[n].extremity[1];
      for(j1 = 0; j1 < 3; j1++)
        if(triangles[edges[n].shared[0]].vertex[j1] == n1)
          break;
      for(j2 = 0; j2 < 3; j2++)
        if(triangles[edges[n].shared[0]].vertex[j2] == n2)
          break;
      if((j2 - j1 == -1) || (j2 - j1 == 2)) {
/*-----------------------------------------------------------------------------
        put in direct rotation, so normal vector will be easier to compute*/
        edges[n].extremity[0] = n2;
        edges[n].extremity[1] = n1;
        }
      }
    }

/*-----------------------------------------------------------------------------
  associate appropriate edges to each element (kind of reverse table) */
  for(k = 0; k < nelts; k++) {
    for(j = 0; j < 3; j++) {
      triangles[k].edges[j] = -1;
      }
    }

  for(n = 0; n < nedges; n++) {
    for(k = 0; k < edges[n].nshared; k++) {
/*       ptr=edges[n].shared[k]; */
      ptr = &(triangles[edges[n].shared[k]]);
      n1 = edges[n].extremity[0];
      n2 = edges[n].extremity[1];
/*-----------------------------------------------------------------------------
      find j1 = reference index of n1 in element */
      for(j1 = 0; j1 < 3; j1++)
        if((*ptr).vertex[j1] == n1)
          break;
/*-----------------------------------------------------------------------------
      find j2 reference index of n2 in element */
      for(j2 = 0; j2 < 3; j2++)
        if((*ptr).vertex[j2] == n2)
          break;
/*------------------------------------------------------------------------------
      Initially we used:
        node 0-1 : side 0, but edge 2 (opposite to node 2)
        node 1-2 : side 1, but edge 0
        node 2-0 : side 2, but edge 1
      Now we use:
        side index in triangle have been changed to match edge index*/
      if((j2 - j1 == 1) || (j2 - j1 == -2)) {
/*------------------------------------------------------------------------------
        direct orientation, use j1 (smaller node index in element)
        to select edge index in triangle arrays nx,ny and l*/
        switch (j1) {
          case 0:
            (*ptr).edges[2] = n;
            edges[n].eindex[k] = 2;
            break;
          case 1:
            (*ptr).edges[0] = n;
            edges[n].eindex[k] = 0;
            break;
          case 2:
            (*ptr).edges[1] = n;
            edges[n].eindex[k] = 1;
            break;
          }
/*------------------------------------------------------------------------------
        store extremity's index in triangle, in direct order*/
        edges[n].vindex[k][0] = j1;
        edges[n].vindex[k][1] = j2;
        }
      else {
/*------------------------------------------------------------------------------
        reverse orientation, use j2 (smaller node index in element)
        to select edge index in triangle arrays nx,ny and l*/
        switch (j2) {
          case 0:
            (*ptr).edges[2] = n;
            edges[n].eindex[k] = 2;
            break;
          case 1:
            (*ptr).edges[0] = n;
            edges[n].eindex[k] = 0;
            break;
          case 2:
            (*ptr).edges[1] = n;
            edges[n].eindex[k] = 1;
            break;
          }
/*------------------------------------------------------------------------------
        store extremity's index in triangle, in direct order*/
        edges[n].vindex[k][0] = j2;
        edges[n].vindex[k][1] = j1;
        }
/*------------------------------------------------------------------------------
      side index in triangle have been changed to match edge index*/
      j=edges[n].eindex[k];
      edges[n].Tx =  (*ptr).ny[j];
      edges[n].Ty = -(*ptr).nx[j];
      edges[n].L  =  (*ptr).l[j];
      }
    }

/*------------------------------------------------------------------------------
  set edge neighbours list */
  for(n = 0; n < nedges; n++) {
    edges[n].ngh=new int [edges[n].nshared*2];
    edges[n].nngh=0;
    }

  for(n = 0; n < nedges; n++) {
    for(k=0;k<edges[n].nshared;k++) {
      m=edges[n].shared[k];
      j=edges[n].eindex[k];
      jm1=(j+2)%3;
      jp1=(j+1)%3;
      edges[n].ngh[edges[n].nngh] = triangles[m].edges[jm1];
      edges[n].nngh++;
      edges[n].ngh[edges[n].nngh] = triangles[m].edges[jp1];
      edges[n].nngh++;
      }
    }

/*-----------------------------------------------------------------------------
  cheks if previous step is correct */
  chk = 0;
  for(k = 0; k < mesh->ntriangles; k++) {
    for(j = 0; j < 3; j++) {
      if(mesh->triangles[k].edges[j] == -1) {
        chk++;
        STDOUT_BASE_LINE("found anomaly %d at element %d, edge %d\n", chk, k, j);
        }
      }
    }

/*------------------------------------------------------------------------------
  set edges' vertex neighbours cross table */
  status=fe_edge_crosstables01(mesh);

/*------------------------------------------------------------------------------
  set vertices' edge neighbours cross table */
  status=fe_vertex_crosstables02(mesh);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int update_vertextable(mesh_t * mesh, int **neigh, int nndes, int max_nb)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------

  - update the vertex neighbour list from external array
  - set cross tables to zero

----------------------------------------------------------------------*/
{
  int k, l, j, j1, j2, n1, n2, n;
  int nedges, chk = 0, count = 0, tmp = 0;
  int nelt, node;
  int interior = 0, boundary = 0, weird = 0;
  int ring[4] = { 0, 1, 2, 0 };
  int *ncells, ncellsmax = 0, **cells;
  vertex_t *vertices;
  triangle_t *elt, *ptr;

  mesh->nvtxs = nndes;
  mesh->nnghm = max_nb;
  vertices = new vertex_t[nndes];

  for(n = 0; n < nndes; n++) {
    vertices[n].ngh    = new int[max_nb];
    vertices[n].nngh   = 0;
    vertices[n].nedges = 0;
    vertices[n].elmts  = 0;
    for(j = 0; j < max_nb; j++) {
      node = neigh[n][j];
      if(node != -1)
        (vertices[n].nngh)++;
      vertices[n].ngh[j] = neigh[n][j];
      }
    }

  mesh->vertices = vertices;

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_vertex_crosstables01(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------
  build the element cross-tables for vertices
-----------------------------------------------------------------------------*/
{
  int nmax, status;
  int i, j, k, m, n;

/*------------------------------------------------------------------------------
  */
  nmax = 0;
  for(n = 0; n < mesh->nvtxs; n++)
    mesh->vertices[n].nelmts = 0;

  for(j = 0; j < mesh->ntriangles; j++) {
    for(i = 0; i < 3; i++) {
      n = mesh->triangles[j].vertex[i];
      mesh->vertices[n].nelmts++;
      updatemax(&nmax, mesh->vertices[n].nelmts);
      }
    }

  printf("fe_vertex_crosstables01: maximum elements/node = %d \n", nmax);

  for(n = 0; n < mesh->nvtxs; n++) {
    mesh->vertices[n].elmts = new int[mesh->vertices[n].nelmts];
    }

  for(n = 0; n < mesh->nvtxs; n++)
    mesh->vertices[n].nelmts = 0;

  for(j = 0; j < mesh->ntriangles; j++) {
    for(i = 0; i < 3; i++) {
      n = mesh->triangles[j].vertex[i];
      mesh->vertices[n].elmts[mesh->vertices[n].nelmts] = j;
      mesh->vertices[n].nelmts++;
      }
    }
  return (0);

error:
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_vertex_crosstables02(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  build the edges cross-tables for vertices
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int duplicated, nmax, status;
  int i, j, k, l, m, n;
  
  if(mesh->ntriangles==0)
    status=fe_list(mesh);
  if(mesh->ntriangles==0)
    return(-1);
  
/*------------------------------------------------------------------------------
  initialize vertices' edge neighbours list */
  for(n = 0; n < mesh->nvtxs; n++) {
    if(mesh->vertices[n].edges!=0) delete[] mesh->vertices[n].edges;
    mesh->vertices[n].edges=0 ;
    mesh->vertices[n].nedges=0;
    }

/*------------------------------------------------------------------------------
  edges directly connected to vertex*/
  for(n = 0; n < mesh->nedges; n++) {
    for (j=0;j<2;j++) {
      m=mesh->edges[n].extremity[j];
      mesh->vertices[m].nedges++;
      }
    }

/*------------------------------------------------------------------------------
  opposite edge*/
  for(n = 0; n < mesh->ntriangles; n++) {
    for (j=0;j<3;j++) {
      m=mesh->triangles[n].vertex[j];
      mesh->vertices[m].nedges+=1;
      }
    }

  for(n = 0; n < mesh->nvtxs; n++) {
    mesh->vertices[n].edges=new int [mesh->vertices[n].nedges];
    mesh->vertices[n].nedges=0;
    }

  for(n = 0; n < mesh->nedges; n++) {
    for (j=0;j<2;j++) {
      m=mesh->edges[n].extremity[j];
      duplicated=pos(n, mesh->vertices[m].edges,mesh->vertices[m].nedges);
      if(duplicated!=-1) {
        printf("anomaly, duplicated edge=%d, index=%d, vertex=%d\n",n,duplicated,m);
        }
      mesh->vertices[m].edges[mesh->vertices[m].nedges]=n;
      mesh->vertices[m].nedges++;
      }
    }

  for(l = 0; l < mesh->ntriangles; l++) {
    for (j=0;j<3;j++) {
      m=mesh->triangles[l].vertex[j];
      n=mesh->triangles[l].edges[j];
      duplicated=pos(n, mesh->vertices[m].edges,mesh->vertices[m].nedges);
      if(duplicated!=-1) {
        printf("anomaly, duplicated edge=%d, index=%d, vertex=%d\n",n,duplicated,m);
        }
      mesh->vertices[m].edges[mesh->vertices[m].nedges]=n;
      mesh->vertices[m].nedges+=1;
      }
    }
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_Vertex_EdgesXtables(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  
  build the geometrical edges cross-tables for vertices
  
  Compliant for T and Q meshes
  
------------------------------------------------------------------------------*/
{
  int duplicated, nmax, status;
  int i, j, k, l, m, n;
  
  if(mesh->vertices==0 and mesh->edges==0) return(-1);
  
/*------------------------------------------------------------------------------
  initialize vertices' edge neighbours list */
  for(n = 0; n < mesh->nvtxs; n++) {
    if(mesh->vertices[n].edges!=0) delete[] mesh->vertices[n].edges;
    mesh->vertices[n].edges=0 ;
    mesh->vertices[n].nedges=0;
    }

/*------------------------------------------------------------------------------
  edges directly connected to vertex*/
  for(n = 0; n < mesh->nedges; n++) {
    for (j=0;j<2;j++) {
      m=mesh->edges[n].extremity[j];
      mesh->vertices[m].nedges++;
      }
    }

  for(n = 0; n < mesh->nvtxs; n++) {
    mesh->vertices[n].edges=new int [mesh->vertices[n].nedges];
    mesh->vertices[n].nedges=0;
    }

  for(n = 0; n < mesh->nedges; n++) {
    for (j=0;j<2;j++) {
      m=mesh->edges[n].extremity[j];
      duplicated=pos(n, mesh->vertices[m].edges,mesh->vertices[m].nedges);
      if(duplicated!=-1) {
        printf("anomaly, duplicated edge=%d, index=%d, vertex=%d\n",n,duplicated,m);
        }
      mesh->vertices[m].edges[mesh->vertices[m].nedges]=n;
      mesh->vertices[m].nedges++;
      }
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_crosstables_T(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int k, l, i, j, m, n;
  int chk = 0, count = 0, tmp = 0;
  int ring[4] = { 0, 1, 2, 0 };

  if(mesh->edges[0].vertex!=0) return(0);

/*------------------------------------------------------------------------------
  set edges' vertex neighbours list */
  for(n = 0; n < mesh->nedges; n++) {
    mesh->edges[n].vertex=new int [mesh->edges[n].nshared+2];
    mesh->edges[n].nvertex=0;
    }

  for(n = 0; n < mesh->nedges; n++) {
    m=mesh->edges[n].shared[0];
    for (j=0;j<3;j++) {
      mesh->edges[n].vertex[j] = mesh->triangles[m].vertex[j];
      }
    mesh->edges[n].nvertex=3;
    if (mesh->edges[n].nshared==1) continue;
    m=mesh->edges[n].shared[1];
    j=3-(mesh->edges[n].vindex[1][0]+mesh->edges[n].vindex[1][1]);
    mesh->edges[n].vertex[3]=mesh->triangles[m].vertex[j];
    mesh->edges[n].nvertex=4;
    }

/*------------------------------------------------------------------------------
  set vertices' edge neighbours list */
  for(n = 0; n < mesh->nvtxs; n++) {
    mesh->vertices[n].edges =0;
    mesh->vertices[n].nedges=0;
    }

  for(n = 0; n < mesh->nedges; n++) {
    for (j=0;j<2;j++) {
      m=mesh->edges[n].extremity[j];
      mesh->vertices[m].nedges++;
      }
    }

  for(n = 0; n < mesh->ntriangles; n++) {
    for (j=0;j<3;j++) {
      m=mesh->triangles[n].vertex[j];
      mesh->vertices[m].nedges+=1;
      }
    }

  for(n = 0; n < mesh->nvtxs; n++) {
    mesh->vertices[n].edges=new int [mesh->vertices[n].nedges];
    mesh->vertices[n].nedges=0;
    }

  for(n = 0; n < mesh->nedges; n++) {
    for (j=0;j<2;j++) {
      m=mesh->edges[n].extremity[j];
      mesh->vertices[m].edges[mesh->vertices[m].nedges]=n;
      mesh->vertices[m].nedges++;
      }
    }

  for(n = 0; n < mesh->ntriangles; n++) {
    for (j=0;j<3;j++) {
      m=mesh->triangles[n].vertex[j];
      mesh->vertices[m].edges[mesh->vertices[m].nedges]=mesh->triangles[n].edges[j];
      mesh->vertices[m].nedges+=1;
      }
    }

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_crosstables_Q(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int k, l, i, j, m, n;
  int chk = 0, count = 0, tmp = 0;
  int ring[4] = { 0, 1, 2, 0 };

  if(mesh->edges[0].vertex!=0) return(0);

/**-----------------------------------------------------------------------------
  set edges' vertex neighbours list */
  for(n = 0; n < mesh->nedges; n++) {
    mesh->edges[n].vertex=new int [mesh->edges[n].nshared+2];
    mesh->edges[n].nvertex=0;
    }

  for(n = 0; n < mesh->nedges; n++) {
    m=mesh->edges[n].shared[0];
    for (j=0;j<3;j++) {
      mesh->edges[n].vertex[j] = mesh->quadrangles[m].vertex[j];
      }
    mesh->edges[n].nvertex=3;
    if (mesh->edges[n].nshared==1) continue;
    m=mesh->edges[n].shared[1];
    j=3-(mesh->edges[n].vindex[1][0]+mesh->edges[n].vindex[1][1]);
    mesh->edges[n].vertex[3]=mesh->quadrangles[m].vertex[j];
    mesh->edges[n].nvertex=4;
    }

/**-----------------------------------------------------------------------------
  set vertices' edge neighbours list */
  for(n = 0; n < mesh->nvtxs; n++) {
    mesh->vertices[n].edges =0;
    mesh->vertices[n].nedges=0;
    }

  for(n = 0; n < mesh->nedges; n++) {
    for (j=0;j<2;j++) {
      m=mesh->edges[n].extremity[j];
      mesh->vertices[m].nedges++;
      }
    }

  for(n = 0; n < mesh->ntriangles; n++) {
    for (j=0;j<3;j++) {
      m=mesh->triangles[n].vertex[j];
      mesh->vertices[m].nedges+=1;
      }
    }

  for(n = 0; n < mesh->nvtxs; n++) {
    mesh->vertices[n].edges=new int [mesh->vertices[n].nedges];
    mesh->vertices[n].nedges=0;
    }

  for(n = 0; n < mesh->nedges; n++) {
    for (j=0;j<2;j++) {
      m=mesh->edges[n].extremity[j];
      mesh->vertices[m].edges[mesh->vertices[m].nedges]=n;
      mesh->vertices[m].nedges++;
      }
    }

  for(n = 0; n < mesh->ntriangles; n++) {
    for (j=0;j<3;j++) {
      m=mesh->triangles[n].vertex[j];
      mesh->vertices[m].edges[mesh->vertices[m].nedges]=mesh->triangles[n].edges[j];
      mesh->vertices[m].nedges+=1;
      }
    }

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_crosstables(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  F. Lyard 20/02/2016: fe_crosstables needs to be checked and documented
  
  Q version looks funny; xtables seem to be numerical, not geometrical

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;

  if(mesh->ntriangles!=0) {
    status=fe_crosstables_T(mesh);
    }
  else {
    status=fe_crosstables_Q(mesh);
    }
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_edge_crosstables01(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int k, l, i, j, m, n;
/*------------------------------------------------------------------------------
  set edges' vertex neighbours cross table */
  for(n = 0; n < mesh->nedges; n++) {
    if(mesh->edges[n].vertex!=0) delete[] mesh->edges[n].vertex;
    mesh->edges[n].vertex=new int [mesh->edges[n].nshared+2];
    mesh->edges[n].nvertex=0;
    }

  for(n = 0; n < mesh->nedges; n++) {
    if (mesh->edges[n].nshared<1) continue;
    m=mesh->edges[n].shared[0];
    for (j=0;j<3;j++) {
      mesh->edges[n].vertex[j] = mesh->triangles[m].vertex[j];
      }
    mesh->edges[n].nvertex=3;
    if (mesh->edges[n].nshared==1) continue;
    m=mesh->edges[n].shared[1];
    j=3-(mesh->edges[n].vindex[1][0]+mesh->edges[n].vindex[1][1]);
    mesh->edges[n].vertex[3]=mesh->triangles[m].vertex[j];
    mesh->edges[n].nvertex=4;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_element_crosstables(mesh_t mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j,k,l,m,m1,m2,mm,mmm,n,pos;
  int ncols,*tmp,count,nngh;

  if(mesh.ntriangles<=0) return(-1);

  for(m = 0; m < mesh.ntriangles; m++) {
    mesh.triangles[m].nngh=0;
    for(k = 0; k <  3; k++) {
      mesh.triangles[m].ngh[k]=-1;
      mesh.triangles[m].opposite[k]=-1;
      }
    }

  for(n = 0; n < mesh.nedges; n++) {
    if(mesh.edges[n].nshared==2) {
      m1=mesh.edges[n].shared[0];
      m2=mesh.edges[n].shared[1];
      mesh.triangles[m1].ngh[mesh.triangles[m1].nngh]=m2;
      mesh.triangles[m2].ngh[mesh.triangles[m2].nngh]=m1;
      mesh.triangles[m1].nngh++;
      mesh.triangles[m2].nngh++;
      }
    }

  for(m = 0; m < mesh.ntriangles; m++) {
    for(j = 0; j <  3; j++) {
      n = mesh.triangles[m].edges[j];
      if(mesh.edges[n].nshared==1) {
        mesh.triangles[m].opposite[j]=-1;
        }
      else {
        if(mesh.edges[n].shared[0]==m) mesh.triangles[m].opposite[j]=mesh.edges[n].shared[1];
        else mesh.triangles[m].opposite[j]=mesh.edges[n].shared[0];
        }
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_element_surrounding(mesh_t mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,mm,n,status;
  int duplicated,*incidence;
  int ncols;

  ncols=mesh.nnghm*3;
/* *----------------------------------------------------------------------------
  bulid element surrounding list */
  incidence=new int[ncols];

  for(m=0;m<mesh.ntriangles;m++) {
    mesh.triangles[m].nsurrounding=0;
    for (i=0;i<3;i++) {
      n=mesh.triangles[m].vertex[i];
      for (k=0;k<mesh.vertices[n].nelmts;k++) {
        mm=mesh.vertices[n].elmts[k];
        duplicated=pos(mm,incidence,mesh.triangles[m].nsurrounding);
        if(duplicated==-1) {
          if(mesh.triangles[m].nsurrounding==ncols) TRAP_ERR_EXIT(-1, "array index overflow");
          incidence[mesh.triangles[m].nsurrounding]=mm;
          mesh.triangles[m].nsurrounding++;
          }
        }
      }
    mesh.triangles[m].surrounding=new int[mesh.triangles[m].nsurrounding];
    mesh.triangles[m].nsurrounding=0;
    for (i=0;i<3;i++) {
      n=mesh.triangles[m].vertex[i];
      for (k=0;k<mesh.vertices[n].nelmts;k++) {
        mm=mesh.vertices[n].elmts[k];
        duplicated=pos(mm,incidence,mesh.triangles[m].nsurrounding);
        if(duplicated==-1) {
          mesh.triangles[m].surrounding[mesh.triangles[m].nsurrounding]=mm;
          mesh.triangles[m].nsurrounding++;
          }
        }
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int build_vindex(mesh_t mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  build edge neighbours index (element-wise) table:

    - egde neighbours
    - connected vertices

-------------------------------------------------------------------------------*/
{
  int i,j,k,l,m,n,status;
  int n1,n2,n3,j1,j2,j3;
  triangle_t *ptr;

  for(n = 0; n < mesh.nedges; n++) {
    for(k = 0; k < mesh.edges[n].nshared; k++) {
/*       ptr=edges[n].shared[k]; */
      ptr = &(mesh.triangles[mesh.edges[n].shared[k]]);
      n1 = mesh.edges[n].extremity[0];
      n2 = mesh.edges[n].extremity[1];
/*-----------------------------------------------------------------------------
      find j1 = reference index of n1 in element */
      for(j1 = 0; j1 < 3; j1++)
        if((*ptr).vertex[j1] == n1)
          break;
/*-----------------------------------------------------------------------------
      find j2 reference index of n2 in element */
      for(j2 = 0; j2 < 3; j2++)
        if((*ptr).vertex[j2] == n2)
          break;
/*------------------------------------------------------------------------------
      Initially we used:
        node 0-1 : side 0, but edge 2 (opposite to node 2)
        node 1-2 : side 1, but edge 0
        node 2-0 : side 2, but edge 1
      Now we use:
        side index in triangle have been changed to match edge index*/
      if((j2 - j1 == 1) || (j2 - j1 == -2)) {
/*------------------------------------------------------------------------------
        direct orientation, use j1 (smaller node index in element)
        to select edge index in triangle arrays nx,ny and l*/
        switch (j1) {
          case 0:
            (*ptr).edges[2] = n;
            mesh.edges[n].eindex[k] = 2;
            break;
          case 1:
            (*ptr).edges[0] = n;
            mesh.edges[n].eindex[k] = 0;
            break;
          case 2:
            (*ptr).edges[1] = n;
            mesh.edges[n].eindex[k] = 1;
            break;
          }
/*------------------------------------------------------------------------------
        store extremity's index in triangle, in direct order*/
        mesh.edges[n].vindex[k][0] = j1;
        mesh.edges[n].vindex[k][1] = j2;
        }
      else {
/*------------------------------------------------------------------------------
        reverse orientation, use j2 (smaller node index in element)
        to select edge index in triangle arrays nx,ny and l*/
        switch (j2) {
          case 0:
            (*ptr).edges[2] = n;
            mesh.edges[n].eindex[k] = 2;
            break;
          case 1:
            (*ptr).edges[0] = n;
            mesh.edges[n].eindex[k] = 0;
            break;
          case 2:
            (*ptr).edges[1] = n;
            mesh.edges[n].eindex[k] = 1;
            break;
          }
/*------------------------------------------------------------------------------
        store extremity's index in triangle, in direct order*/
        mesh.edges[n].vindex[k][0] = j2;
        mesh.edges[n].vindex[k][1] = j1;
        }
/*------------------------------------------------------------------------------
      side index in triangle have been changed to match edge index*/
//       j=edges[n].eindex[k];
//       edges[n].Tx =  (*ptr).ny[j];
//       edges[n].Ty = -(*ptr).nx[j];
//       edges[n].L  =  (*ptr).l[j];
      }
    }

}

