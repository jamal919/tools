/*******************************************************************************

  T-UGO tools, 2006-2018

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief finite elements geometry edition functions
*/
/*----------------------------------------------------------------------------*/

#include <stdio.h>

#include "tools-structures.h"
#include "constants.h"

// #include "fe.def" // for MESH_FILE_FORMAT_TRIGRID
#include "fe.h"

#include "geo.h"
#include "zapper.h"


/*----------------------------------------------------------------------------*/
/// make sure fe_order2() is exported even though it is inline
/** whatever the compiler version. This might be related to
https://gcc.gnu.org/bugzilla/show_bug.cgi?id=55015#c13
KEEP THIS AT THE BEGIN OF THE FILE WHERE fe_order2 is defined!
SEE ALSO AT THE END OF THIS FILE. */
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define fe_order2 fe_order2_inline
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_minmax(const mesh_t & mesh,frame_t & frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  get mesh frame
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int  n;
  double t,p;

  if (mesh.vertices==0) return(EFAULT);

  frame.xmin=mesh.vertices[0].lon;
  frame.xmax=mesh.vertices[0].lon;
  frame.ymin=mesh.vertices[0].lat;
  frame.ymax=mesh.vertices[0].lat;

  for(n=1;n<mesh.nvtxs;n++) {
    t=mesh.vertices[n].lon;
    if(mesh.type==0) t=degree_recale(t,mesh.vertices[0].lon);
    p=mesh.vertices[n].lat;
    updatemin(&frame.xmin,t);
    updatemax(&frame.xmax,t);
    updatemin(&frame.ymin,p);
    updatemax(&frame.ymax,p);
    }
  
  
  if(mesh.triangles!=0){
/*------------------------------------------------------------------------------
    scan polar triangles */
    int n1,n2,n3;
    range_t<double> tR;
    
    for(n=0;n<mesh.ntriangles;n++) {
      const triangle_t *triangle=&mesh.triangles[n];
      
      n1=triangle->vertex[0];
      n2=triangle->vertex[1];
      n3=triangle->vertex[2];
      
      tR.init();
      tR<<mesh.vertices[n1].lon;
      tR<<mesh.vertices[n2].lon;
      tR<<mesh.vertices[n3].lon;
      
      const double dt=tR.d();
      
      if(dt<-180. or 180.<dt){
        /* triangle is polar */
        frame.xmin=-180.;
        frame.xmax=+180.;
        
        if(mesh.vertices[n1].lat>0)
          frame.ymax=+90.;
        else
          frame.ymin=-90.;
        
        }
      
      }
    
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  frame_t fe_ElementExtent(const mesh_t & mesh, const triangle_t & element)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k, n;
  double t,p,base;
  frame_t frame;
  bool spherical=(mesh.type==0);

  n=element.vertex[0];
  frame.xmin=mesh.vertices[n].lon;
  frame.xmax=mesh.vertices[n].lon;
  frame.ymin=mesh.vertices[n].lat;
  frame.ymax=mesh.vertices[n].lat;
  
  base=mesh.vertices[n].lon;

  for(k=1;k<3;k++) {
    n=element.vertex[k];
    t=mesh.vertices[n].lon;
    if(spherical) t=degree_recale(t,base);
    p=mesh.vertices[n].lat;
    updatemin(&frame.xmin,t);
    updatemax(&frame.xmax,t);
    updatemin(&frame.ymin,p);
    updatemax(&frame.ymax,p);
    }
  
  const double dt=frame.xmax-frame.xmin;
  if(dt<-180. or 180.<dt){
    /* triangle is polar */
    if(frame.ymax>0)
      frame.ymax= 90.;
    else
      frame.ymin=-90.;
    frame.xmin=-INFINITY;
    frame.xmax=+INFINITY;
    }

  return(frame);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  frame_t fe_ElementExtent(const mesh_t & mesh, const quadrangle_t & element)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k, n;
  double t,p,base;
  frame_t frame;
  bool spherical=(mesh.type==0);

  n=element.vertex[0];
  frame.xmin=mesh.vertices[n].lon;
  frame.xmax=mesh.vertices[n].lon;
  frame.ymin=mesh.vertices[n].lat;
  frame.ymax=mesh.vertices[n].lat;
  
  base=mesh.vertices[n].lon;

  for(k=1;k<4;k++) {
    n=element.vertex[k];
    t=mesh.vertices[n].lon;
    if(spherical) t=degree_recale(t,base);
    p=mesh.vertices[n].lat;
    updatemin(&frame.xmin,t);
    updatemax(&frame.xmax,t);
    updatemin(&frame.ymin,p);
    updatemax(&frame.ymax,p);
    }

  return(frame);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <class T> int fe_InExtent_template(const mesh_t & mesh, const T & element, double x, double y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int chk;
  frame_t frame;
  double base;
  bool spherical=(mesh.type==0);
  
  frame=fe_ElementExtent(mesh, element);
  
  if(spherical) {
    base=frame.x_center();
    x=degree_recale(x,base);
    }
  
  chk=frame.inside(x,y);
  
  return(chk);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_InExtent(const mesh_t & mesh, const triangle_t & element, double x, double y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int chk;
  
  chk=fe_InExtent_template(mesh, element, x, y);
  
  return(chk);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_InExtent(const mesh_t & mesh, const quadrangle_t & element, double x, double y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int chk;
  
  chk=fe_InExtent_template(mesh, element, x, y);
  
  return(chk);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int report_edge_table(int interior,int boundary,int weird,int presumed,int nelts)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  printf("found %d interior edges\n", interior);
  printf("found %d boundary edges\n", boundary);
  printf("found %d weird edges (should be zero) \n", weird);
  
  int difference;
  difference=presumed-nelts;
  printf("found presumedly %d elements, actually %d, difference: %d\n",presumed, nelts, difference);
  
  return difference;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_edge_table(mesh_t * mesh, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l, j, j1, j2, jm1, jp1, m, n1, n2, n;
  int nedges, chk = 0, count = 0;
  int nelts, nndes;
  int interior = 0, boundary = 0, weird = 0;
  int *ncells=NULL, **cells=NULL;
  edge_t *edges=NULL;
  triangle_t *elt=NULL;
  
  struct timeval b4;
  gettimeofday(&b4);
  
  nedges = 0;

  elt = mesh->triangles;
  nelts = mesh->ntriangles;
  nndes = mesh->nvtxs;

/*-----------------------------------------------------------------------------
  count edges*/
  for(n = 0; n < nndes; n++) {
    for(j = 0; j < mesh->vertices[n].nngh; j++) {
      n2 = mesh->vertices[n].ngh[j];
      if(n2 > n)
        nedges++;
      }
    }
  if(verbose) STDERR_BASE_LINE_FUNC("%gs\n",difftime(&b4));

/*-----------------------------------------------------------------------------
  count number of element connected to a given P1-node*/
  ncells = aset(nndes + 1,0);
  if(verbose) STDERR_BASE_LINE_FUNC("%gs\n",difftime(&b4));
  cells = new int *[nndes + 1];
  if(verbose) STDERR_BASE_LINE_FUNC("%gs\n",difftime(&b4));
  
  for(k = 0; k < nelts; k++) {
    for(j = 0; j < 3; j++) {
      n = elt[k].vertex[j];
      ncells[n]++;
      }
    }
  if(verbose) STDERR_BASE_LINE_FUNC("%gs\n",difftime(&b4));

  for(n = 0; n < nndes; n++) {
    if(ncells[n] == 0) STDERR_BASE_LINE_FUNC("unused node %d %d\n", n, mesh->vertices[n].nngh);
    cells[n] = new int[ncells[n]];
    ncells[n] = 0;
    }
  if(verbose) STDERR_BASE_LINE_FUNC("%gs\n",difftime(&b4));

/*-----------------------------------------------------------------------------
  build element connection table*/
  for(k = 0; k < nelts; k++) {
    for(j = 0; j < 3; j++) {
      n = elt[k].vertex[j];
      cells[n][ncells[n]] = k;
      ncells[n]++;
      }
    }
  if(verbose) STDERR_BASE_LINE_FUNC("%gs\n",difftime(&b4));

  edges = new edge_t[nedges];
  if(verbose) STDERR_BASE_LINE_FUNC("%gs\n",difftime(&b4));
  for(n = 0; n < nedges; n++) {
    edges[n].nshared = 0;
    }
  if(verbose) STDERR_BASE_LINE_FUNC("%gs\n",difftime(&b4));

  //#pragma omp parallel for private(n1,j,n2,k,l)
  for(n1 = 0; n1 < nndes; n1++) {
    const vertex_t *vertex=&mesh->vertices[n1];
    for(j = 0; j < vertex->nngh; j++) {
      n2 = vertex->ngh[j];
      if(n2 > n1){
/*-----------------------------------------------------------------------------
        create a new edge entry*/
        edge_t *edge;
        //#pragma omp critical(count)
        {
        edge=&edges[count];
        count++;
        }
      
        edge->extremity[0] = n1;
        edge->extremity[1] = n2;
/*-----------------------------------------------------------------------------
        scan elements shared by n1 and n2*/
        for(k = 0; k < ncells[n1]; k++)
          for(l = 0; l < ncells[n2]; l++)
            if(cells[n1][k] == cells[n2][l]) {
/*               edge->shared[edge->nshared]=&(elt[cells[n1][k]]); */
              if(edge->nshared > 1) {
                printf("found %d as weird edge, nshared=%d\n", count-1, edge->nshared);
                printf("elements: %d %d\n", edge->shared[0],edge->shared[1]);
                printf("vertices: %d %d\n", edge->extremity[0],edge->extremity[1]);
                }
              edge->shared[edge->nshared] = cells[n1][k];
              edge->nshared++;
              }
        }
      }
    }
  if(verbose) STDERR_BASE_LINE_FUNC("%gs\n",difftime(&b4));
  
  deletep(&ncells);
  deletep2D(&cells,nndes);
  if(verbose) STDERR_BASE_LINE_FUNC("%gs\n",difftime(&b4));

  boundary = 0;
  interior = 0;
  for(n = 0; n < nedges; n++) {
    edge_t *edge=&edges[n];
    if(edge->nshared == 0) {
      weird++;
      }
    if(edge->nshared > 2) {
      weird++;
      printf("found %d as weird edge, nshared=%d\n", n, edge->nshared);
      printf("elements: %d %d %d\n", edge->shared[0],edge->shared[1],edge->shared[2]);
      printf("vertices: %d %d\n", edge->extremity[0],edge->extremity[1]);
      }
    if(edge->nshared == 1)
      boundary++;
    if(edge->nshared == 2)
      interior++;
    }
  if(verbose) STDERR_BASE_LINE_FUNC("%gs\n",difftime(&b4));

  if(weird!=0) {
    STDOUT_BASE_LINE_FUNC("\nfound %d total edges\n", nedges);
    report_edge_table( interior, boundary, weird,(interior * 2 + boundary) / 3, nelts);
    }

/* *-----------------------------------------------------------------------------
  scan boundary edges for outward normal determination*/
  for(n = 0; n < nedges; n++) {
    edge_t *edge=&edges[n];
    if(edge->nshared == 1) {
      n1 = edge->extremity[0];
      n2 = edge->extremity[1];
      for(j1 = 0; j1 < 3; j1++)
        if(elt[edge->shared[0]].vertex[j1] == n1)
          break;
      for(j2 = 0; j2 < 3; j2++)
        if(elt[edge->shared[0]].vertex[j2] == n2)
          break;
      if((j2 - j1 == -1) || (j2 - j1 == 2)) {
/*-----------------------------------------------------------------------------
        put in direct rotation, so normal vector will be easier to compute*/
        edge->extremity[0] = n2;
        edge->extremity[1] = n1;
        }
      }
    }
  if(verbose) STDERR_BASE_LINE_FUNC("%gs\n",difftime(&b4));

/*-----------------------------------------------------------------------------
  associate appropriate edges to each element (kind of reverse table) */
  for(k = 0; k < nelts; k++) {
    aset(elt[k].edges,3,-1);
    }
  if(verbose) STDERR_BASE_LINE_FUNC("%gs\n",difftime(&b4));

  for(n = 0; n < nedges; n++) {
    edge_t *edge=&edges[n];
    for(k = 0; k < edge->nshared; k++) {
      triangle_t *triangle = &elt[edge->shared[k]];
      n1 = edge->extremity[0];
      n2 = edge->extremity[1];
/*-----------------------------------------------------------------------------
      find j1 = reference index of n1 in element */
      for(j1 = 0; j1 < 3; j1++)
        if(triangle->vertex[j1] == n1)
          break;
/*-----------------------------------------------------------------------------
      find j2 reference index of n2 in element */
      for(j2 = 0; j2 < 3; j2++)
        if(triangle->vertex[j2] == n2)
          break;
/* *-----------------------------------------------------------------------------
      Initially we used:
        node 0-1 : side 0, but edge 2 (opposite to node 2)
        node 1-2 : side 1, but edge 0
        node 2-0 : side 2, but edge 1
      Now we use:
        side index in triangle have been changed to match edge index*/
      if((j2 - j1 == 1) || (j2 - j1 == -2)) {
/* *-----------------------------------------------------------------------------
        direct orientation, use j1 (smaller node index in element)
        to select edge index in triangle arrays nx,ny and l*/
        switch (j1) {
          case 0:
            triangle->edges[2] = n;
            edge->eindex[k] = 2;
            break;
          case 1:
            triangle->edges[0] = n;
            edge->eindex[k] = 0;
            break;
          case 2:
            triangle->edges[1] = n;
            edge->eindex[k] = 1;
            break;
          }
/* *-----------------------------------------------------------------------------
        store extremity's index in triangle, in direct order*/
        edge->vindex[k][0] = j1;
        edge->vindex[k][1] = j2;
        }
      else {
/* *-----------------------------------------------------------------------------
        reverse orientation, use j2 (smaller node index in element)
        to select edge index in triangle arrays nx,ny and l*/
        switch (j2) {
          case 0:
            triangle->edges[2] = n;
            edge->eindex[k] = 2;
            break;
          case 1:
            triangle->edges[0] = n;
            edge->eindex[k] = 0;
            break;
          case 2:
            triangle->edges[1] = n;
            edge->eindex[k] = 1;
            break;
          }
/* *-----------------------------------------------------------------------------
        store extremity's index in triangle, in direct order*/
        edge->vindex[k][0] = j2;
        edge->vindex[k][1] = j1;
        }
/* *-----------------------------------------------------------------------------
      side index in triangle have been changed to match edge index*/
      j=edge->eindex[k];
      edge->Tx =  triangle->ny[j];
      edge->Ty = -triangle->nx[j];
      edge->L  =  triangle->l[j];
      }
    }

/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

    05/08/2008

    following block should be encapsulated in a routine as it is of
    general and repetitive interest

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/* *-----------------------------------------------------------------------------
  construct edge neighbours list */
  if(verbose) STDERR_BASE_LINE_FUNC("%gs\n",difftime(&b4));
  for(n = 0; n < nedges; n++) {
    edge_t *edge=&edges[n];
    edge->ngh=new int [edge->nshared*2];
    edge->nngh=0;
    }
  if(verbose) STDERR_BASE_LINE_FUNC("%gs\n",difftime(&b4));

  for(n = 0; n < nedges; n++) {
    edge_t *edge=&edges[n];
    for(k=0;k<edge->nshared;k++) {
      m=edge->shared[k];
      j=edge->eindex[k];
      jm1=(j+2)%3;
      jp1=(j+1)%3;
      edge->ngh[edge->nngh] = elt[m].edges[jm1];
      edge->nngh++;
      edge->ngh[edge->nngh] = elt[m].edges[jp1];
      edge->nngh++;
      }
    }
  if(verbose) STDERR_BASE_LINE_FUNC("%gs\n",difftime(&b4));

  mesh->nedges = nedges;
  deletep(&mesh->edges);
  mesh->edges = edges;

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
  if(verbose) STDERR_BASE_LINE_FUNC("%gs\n",difftime(&b4));

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_edgetable_T(mesh_t *mesh, int task, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l, j, j1, j2, n1, n2, n, status;
  int nedges, chk = 0, count = 0;
  int nelts, nndes;
  int interior = 0, boundary = 0, weird = 0;
  int *ncells, ncellsmax = 0, **cells;
  edge_t *edges;
  triangle_t *elts, *ptr;

  nedges = 0;
  int nprocs __attribute__((unused)) =initialize_OPENMP(-1,0);

  elts  = mesh->triangles;
  nelts = mesh->ntriangles;
  nndes = mesh->nvtxs;
  
  if(mesh->edges!=0){
    
    for(n=0;n<mesh->nedges;n++){
      mesh->edges[n].destroy();
      }
    
    deletep(&mesh->edges);
    }

/*-----------------------------------------------------------------------------
  count edges*/
  for(n = 0; n < nndes; n++) {
    for(j=0; j<mesh->vertices[n].nngh;j++) {
      n2 = mesh->vertices[n].ngh[j];
      if(n2 > n)
        nedges++;
      }
    }

/*-----------------------------------------------------------------------------
  count max number of element connected to a given P1-node*/
  ncells = new int  [nndes + 1];
  cells  = new int *[nndes + 1];

  for(n = 0; n < nndes; n++)
    ncells[n] = 0;

  for(k = 0; k < nelts; k++) {
    for(j = 0; j < 3; j++) {
      n = elts[k].vertex[j];
      ncells[n]++;
      }
    }

  for(n = 0; n < nndes; n++) {
    if(ncells[n] == 0) {
      printf("unused node %d (nngh=%d) : position lon=%lf lat=%lf\n", n, mesh->vertices[n].nngh,mesh->vertices[n].lon,mesh->vertices[n].lat);
      }
    }

  for(n = 0; n < nndes; n++)
    updatemax(&ncellsmax, ncells[n]);
  for(n = 0; n < nndes; n++)
    cells[n] = new int[ncellsmax];
  for(n = 0; n < nndes; n++)
    ncells[n] = 0;

/*-----------------------------------------------------------------------------
  build element connection table*/
  for(k = 0; k < nelts; k++) {
    for(j = 0; j < 3; j++) {
      n = elts[k].vertex[j];
      cells[n][ncells[n]] = k;
      ncells[n]++;
      }
    }

  edges = new edge_t[nedges];
  for(n = 0; n < nedges; n++) {
    edges[n].nshared = 0;
    }

  for(n1 = 0; n1 < nndes; n1++) {
    for(j=0; j<mesh->vertices[n1].nngh;j++) {
      n2 = mesh->vertices[n1].ngh[j];
      if(n2 > n1) {
/*-----------------------------------------------------------------------------
        create a new edge entry*/
        edges[count].extremity[0] = n1;
        edges[count].extremity[1] = n2;
/*-----------------------------------------------------------------------------
        scan elements shared by n1 and n2*/
        for(k = 0; k < ncells[n1]; k++)
          for(l = 0; l < ncells[n2]; l++)
            if(cells[n1][k] == cells[n2][l]) {
/*               edges[count].shared[edges[count].nshared]=&(elt[cells[n1][k]]); */
              if(edges[count].nshared > 1) {
                printf("found edge %d to be weird edge, nshared=%d\n", count, edges[count].nshared);
                printf("elements: %d %d\n", edges[count].shared[0],edges[count].shared[1]);
                printf("vertices: %d %d\n", edges[count].extremity[0],edges[count].extremity[1]);
                }
              edges[count].shared[edges[count].nshared] = cells[n1][k];
              edges[count].nshared++;
              }
        count++;
        }
      }
    }

  boundary = 0;
  interior = 0;
  for(n = 0; n < nedges; n++) {
//     if(edges[n].nshared == 0) {
//       weird++;
//       printf("found %d as weird edge, nshared=%d\n", n, edges[n].nshared);
//       }
    if(edges[n].nshared > 2) {
      weird++;
      printf("found %d as weird edge, nshared=%d\n", n, edges[n].nshared);
      printf("elements indices: %d %d %d\n", edges[n].shared[0],edges[n].shared[1],edges[n].shared[2]);
      printf("vertices indices: %d %d\n", edges[n].extremity[0],edges[n].extremity[1]);
      int n1=edges[n].extremity[0],n2=edges[n].extremity[1];
      printf("node position %6d lon=%9.3lf lat=%9.3lf\n",n1,mesh->vertices[n1].lon,mesh->vertices[n1].lat);
      printf("node position %6d lon=%9.3lf lat=%9.3lf\n",n2,mesh->vertices[n2].lon,mesh->vertices[n2].lat);
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

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  scan boundary edges for outward normal determination
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(n = 0; n < nedges; n++) {
    if(edges[n].nshared == 1) {
      n1 = edges[n].extremity[0];
      n2 = edges[n].extremity[1];
      for(j1 = 0; j1 < 3; j1++)
        if(elts[edges[n].shared[0]].vertex[j1] == n1)
          break;
      for(j2 = 0; j2 < 3; j2++)
        if(elts[edges[n].shared[0]].vertex[j2] == n2)
          break;
      if((j2 - j1 == -1) || (j2 - j1 == 2)) {
/*-----------------------------------------------------------------------------
        put in direct rotation, so normal vector will be easier to compute*/
        edges[n].extremity[0] = n2;
        edges[n].extremity[1] = n1;
        }
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  associate appropriate edges to each element (kind of reverse table) 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(k = 0; k < nelts; k++) {
    for(j = 0; j < 3; j++) {
      elts[k].edges[j] = -1;
      }
    }

  for(n = 0; n < nedges; n++) {
    for(k = 0; k < edges[n].nshared; k++) {
      ptr = &(elts[edges[n].shared[k]]);
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
/**-----------------------------------------------------------------------------
      Initially we used:
        node 0-1 : side 0, but edge 2 (opposite to node 2)
        node 1-2 : side 1, but edge 0
        node 2-0 : side 2, but edge 1
      Now we use:
        side index in triangle have been changed to match edge index*/
      if((j2 - j1 == 1) || (j2 - j1 == -2)) {
/**-----------------------------------------------------------------------------
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
/**-----------------------------------------------------------------------------
        store extremity's index in triangle, in direct order*/
        edges[n].vindex[k][0] = j1;
        edges[n].vindex[k][1] = j2;
        }
      else {
/**-----------------------------------------------------------------------------
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
/**-----------------------------------------------------------------------------
        store extremity's index in triangle, in direct order*/
        edges[n].vindex[k][0] = j2;
        edges[n].vindex[k][1] = j1;
        }
/**-----------------------------------------------------------------------------
      side index in triangel have been changed to match edge index*/
      j=edges[n].eindex[k];
      edges[n].Tx =  (*ptr).ny[j];
      edges[n].Ty = -(*ptr).nx[j];
      edges[n].L  =  (*ptr).l[j];
      }
    }

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

    05/08/2008

    following block should be encapsulated in a routine as it is of
    general and repetitive interest

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/**-----------------------------------------------------------------------------
  construct edge neighbours list */
// #pragma omp parallel for
  for(n = 0; n < nedges; n++) {
    if(edges[n].ngh!=0) delete[] edges[n].ngh;
    edges[n].ngh=new int [edges[n].nshared*2];
    edges[n].nngh=0;
    }

#pragma omp parallel for
  for(n = 0; n < nedges; n++) {
    for(int k=0;k<edges[n].nshared;k++) {
      int m=edges[n].shared[k];
      int j=edges[n].eindex[k];
      int jm1=(j+2)%3;
      int jp1=(j+1)%3;
      edges[n].ngh[edges[n].nngh] = elts[m].edges[jm1];
      edges[n].nngh++;
      edges[n].ngh[edges[n].nngh] = elts[m].edges[jp1];
      edges[n].nngh++;
      }
    }

  mesh->nedges = nedges;
  mesh->edges  = edges;

/*-----------------------------------------------------------------------------
  cheks if previous step is correct */
  if(debug) {
    chk = 0;
    for(k = 0; k < mesh->ntriangles; k++) {
      for(j = 0; j < 3; j++) {
        if(mesh->triangles[k].edges[j] == -1) {
          chk++;
          STDOUT_BASE_LINE("found anomaly %d at element %d, edge %d\n", chk, k, j);
          }
        }
      }
    }

  zaparr(ncells);
  for(n = 0; n < nndes; n++)
    zaparr(cells[n]);
  zaparr(cells);

/**-----------------------------------------------------------------------------
  set edges' vertex neighbours cross table */
  if(task==0) status=fe_edge_crosstables01(mesh);

/**-----------------------------------------------------------------------------
  set vertices' edge neighbours cross table */
  if(task==0) status=fe_vertex_crosstables02(mesh);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_edgetable(mesh_t *mesh, int task, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  switch(mesh->ntriangles) {
    case 0:
      status=fe_edgetable_Q(mesh, verbose);
      break;
    default:
      status=fe_edgetable_T(mesh, task, debug);
      break;
    }
  return (status);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t fe_doubleT(const mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,n1,n2,n,status;
  int   nedges;
  int   nelts,nndes;
  edge_t *edges=NULL;
  triangle_t *elt=NULL;
  double t1,t2,p1,p2,h1,h2;
  mesh_t work;

  elt=mesh.triangles;

  nelts =mesh.ntriangles;
  nndes =mesh.nvtxs;
  nedges=mesh.nedges;
  edges =mesh.edges;

  work.type=0;
/*-----------------------------------------------------------------------------
  final mesh : original vertices + 1 per edges vertex */
  work.vertices=new vertex_t[nedges+nndes];

  for (n=0;n<nndes;n++) {
    work.vertices[n].lon=mesh.vertices[n].lon;
    work.vertices[n].lat=mesh.vertices[n].lat;
    work.vertices[n].h=mesh.vertices[n].h;
    work.vertices[n].nngh=0;
    }

  for (n=0;n<nedges;n++) {
    n1=edges[n].extremity[0];
    n2=edges[n].extremity[1];
/*-----------------------------------------------------------------------------
    create a new node at the middle of the edge */
    t1=mesh.vertices[n1].lon;
    t2=mesh.vertices[n2].lon;
    if(mesh.type==0) t2=degree_recale(t2,t1);
    p1=mesh.vertices[n1].lat;
    p2=mesh.vertices[n2].lat;
    h1=mesh.vertices[n1].h;
    h2=mesh.vertices[n2].h;
    work.vertices[n+nndes].lon=0.5*(t1+t2);
    work.vertices[n+nndes].lat=0.5*(p1+p2);
    work.vertices[n+nndes].h=0.5*(h1+h2);
    work.vertices[n+nndes].nngh=0;
    }

 work.nvtxs=nedges+nndes;
 work.ntriangles=0;
 work.triangles= new triangle_t[4*nelts];

/*-----------------------------------------------------------------------------
  list old element to create new elements */
  for(k=0; k<nelts; k++) {
    work.triangles[work.ntriangles].vertex[0]=elt[k].vertex[0];
    work.triangles[work.ntriangles].vertex[1]=nndes+elt[k].edges[2];
    work.triangles[work.ntriangles].vertex[2]=nndes+elt[k].edges[1];
    status=fe_initaffine(&work,work.ntriangles);
    work.ntriangles++;
    work.triangles[work.ntriangles].vertex[0]=elt[k].vertex[1];
    work.triangles[work.ntriangles].vertex[1]=nndes+elt[k].edges[0];
    work.triangles[work.ntriangles].vertex[2]=nndes+elt[k].edges[2];
    status=fe_initaffine(&work,work.ntriangles);
    work.ntriangles++;
    work.triangles[work.ntriangles].vertex[0]=elt[k].vertex[2];
    work.triangles[work.ntriangles].vertex[1]=nndes+elt[k].edges[1];
    work.triangles[work.ntriangles].vertex[2]=nndes+elt[k].edges[0];
    status=fe_initaffine(&work,work.ntriangles);
    work.ntriangles++;
    work.triangles[work.ntriangles].vertex[0]=nndes+elt[k].edges[2];
    work.triangles[work.ntriangles].vertex[1]=nndes+elt[k].edges[0];
    work.triangles[work.ntriangles].vertex[2]=nndes+elt[k].edges[1];
    status=fe_initaffine(&work,work.ntriangles);
    work.ntriangles++;
    }

/*-----------------------------------------------------------------------------
  rebuild neighbour list */
  status=fe_e2n(&work);
//   work.nnghm=mesh.nnghm;
//
//   for (n=0; n<work.nvtxs; n++) {
//     work.vertices[n].ngh= (int *) malloc(work.nnghm*sizeof(int));
//     for(k=0;k<work.nnghm;k++) work.vertices[n].ngh[k]=-1;
//     work.vertices[n].nngh=0;
//     }
//
//   for(m=0;m<work.ntriangles;m++) {
//     for(i=0;i<3;i++) {
//       n=work.triangles[m].vertex[i];
//       for(j=0;j<3;j++) {
//         if(i==j) continue;
//         n1=work.triangles[m].vertex[j];
//         for(k=0;k<work.vertices[n].nngh;k++)
//           if(work.vertices[n].ngh[k]==n1) goto skip2;
//         if(work.vertices[n].nngh==mesh.nnghm) {
//        printf("alert\n");
// /*     return(NULL); */
//        }
//      work.vertices[n].nngh++;
//         work.vertices[n].ngh[k]=n1;
//      skip2:
//         count++;
//      }
//       }
//     }

  return(work);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t fe_doubleQ(const mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,n1,n2,n,status;
  int   nedges;
  int   nelts;
  edge_t *edges=NULL;
  quadrangle_t *elt=NULL;
  double t1,t2,p1,p2,h1,h2;
  mesh_t work;

  elt=mesh.quadrangles;

  nelts =mesh.nquadrangles;
//   nndes =mesh.nvtxs;
  nedges=mesh.nedges;
  edges =mesh.edges;

  work.type=0;
/*-----------------------------------------------------------------------------
  final mesh : original vertices + 1 per edges vertex */
  work.vertices=new vertex_t[mesh.nquadrangles+mesh.nedges+mesh.nvtxs];

  for (n=0;n<mesh.nvtxs;n++) {
    work.vertices[n].lon=mesh.vertices[n].lon;
    work.vertices[n].lat=mesh.vertices[n].lat;
    work.vertices[n].h=mesh.vertices[n].h;
    work.vertices[n].nngh=0;
    }

  for (n=0;n<nedges;n++) {
    n1=edges[n].extremity[0];
    n2=edges[n].extremity[1];
/*-----------------------------------------------------------------------------
    create a new node at the middle of the edge */
    t1=mesh.vertices[n1].lon;
    t2=mesh.vertices[n2].lon;
    if(mesh.type==0) t2=degree_recale(t2,t1);
    p1=mesh.vertices[n1].lat;
    p2=mesh.vertices[n2].lat;
    h1=mesh.vertices[n1].h;
    h2=mesh.vertices[n2].h;
    work.vertices[n+mesh.nvtxs].lon=0.5*(t1+t2);
    work.vertices[n+mesh.nvtxs].lat=0.5*(p1+p2);
    work.vertices[n+mesh.nvtxs].h=0.5*(h1+h2);
    work.vertices[n+mesh.nvtxs].nngh=0;
    }
  
  for(int m=0;m<mesh.nquadrangles;m++) {
    double t,p;
    status=fe_position(mesh, mesh.quadrangles[m], &t, &p,0);
//     for(i=0;i<4;i++) {
//       }
    work.vertices[m+mesh.nvtxs+mesh.nedges].lon=t;
    work.vertices[m+mesh.nvtxs+mesh.nedges].lat=p;
    work.vertices[m+mesh.nvtxs+mesh.nedges].h=0;
    work.vertices[m+mesh.nvtxs+mesh.nedges].nngh=0;
    }

  work.nvtxs=nedges+mesh.nvtxs+mesh.nquadrangles;
  work.nquadrangles=0;
  work.quadrangles= new quadrangle_t[4*nelts];

/*-----------------------------------------------------------------------------
  list old element to create new elements */
  for(k=0; k<nelts; k++) {
    work.quadrangles[work.nquadrangles].vertex[0]=elt[k].vertex[0];
    work.quadrangles[work.nquadrangles].vertex[1]=mesh.nvtxs+elt[k].edges[0];
    work.quadrangles[work.nquadrangles].vertex[2]=mesh.nvtxs+mesh.nedges+k;
    work.quadrangles[work.nquadrangles].vertex[3]=mesh.nvtxs+elt[k].edges[3];
//     status=fe_initaffine(&work,work.nquadrangles);
    work.nquadrangles++;
    work.quadrangles[work.nquadrangles].vertex[0]=elt[k].vertex[1];
    work.quadrangles[work.nquadrangles].vertex[1]=mesh.nvtxs+elt[k].edges[1];
    work.quadrangles[work.nquadrangles].vertex[2]=mesh.nvtxs+mesh.nedges+k;
    work.quadrangles[work.nquadrangles].vertex[3]=mesh.nvtxs+elt[k].edges[0];
//     status=fe_initaffine(&work,work.nquadrangles);
    work.nquadrangles++;
    work.quadrangles[work.nquadrangles].vertex[0]=elt[k].vertex[2];
    work.quadrangles[work.nquadrangles].vertex[1]=mesh.nvtxs+elt[k].edges[2];
    work.quadrangles[work.nquadrangles].vertex[2]=mesh.nvtxs+mesh.nedges+k;
    work.quadrangles[work.nquadrangles].vertex[3]=mesh.nvtxs+elt[k].edges[1];
//     status=fe_initaffine(&work,work.nquadrangles);
    work.nquadrangles++;
    work.quadrangles[work.nquadrangles].vertex[0]=elt[k].vertex[3];
    work.quadrangles[work.nquadrangles].vertex[1]=mesh.nvtxs+elt[k].edges[3];
    work.quadrangles[work.nquadrangles].vertex[2]=mesh.nvtxs+mesh.nedges+k;
    work.quadrangles[work.nquadrangles].vertex[3]=mesh.nvtxs+elt[k].edges[2];
//     status=fe_initaffine(&work,work.nquadrangles);
    work.nquadrangles++;
    }

/*-----------------------------------------------------------------------------
  rebuild neighbour list */
  status=fe_e2n(&work);

  return(work);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t fe_double(const mesh_t  &mesh, int ElementType)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//   int status;
  mesh_t work;
  
  switch(ElementType) {
    case FE_TRIANGLE:
      work=fe_doubleT(mesh);
      break;
    case FE_QUADRANGLE:
      work=fe_doubleQ(mesh);
//       for(size_t m=0;m<mesh->nquadrangles;m++) {
//         status=fe_initaffine_spherical(*mesh, &(mesh->quadrangles[m]), m);
//         }
      break;
    default:
      TRAP_ERR_EXIT(-1, "element type not reckognized\n");
      break;
    }
  return (work);
    
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_order(mesh_t & mesh, int n, double x, double y, int rotation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,m,m1,m2;
  double t1,t2,p1,p2,dummy;
  double *alpha;
  double beta,factor;
  
  vertex_t *vertex=&mesh.vertices[n];
  
  if(vertex->nngh<2) return(0);

  switch (rotation) {
    case -1:
      factor=-1.;
      break;

    case +1:
      factor=+1.;
      break;
    }

  beta=factor*atan2(y,x);
  
  alpha=new double[vertex->nngh];

  t1=vertex->lon;
  p1=vertex->lat;
  for(k=0;k<vertex->nngh;k++) {
    m=vertex->ngh[k];
    const vertex_t *neigh=&mesh.vertices[m];
    t2=neigh->lon;
    if(mesh.type==0) degree_recale(&t2,t1);
    p2=neigh->lat;
    alpha[k]=factor*atan2(p2-p1,t2-t1);
    if(alpha[k]<=beta) alpha[k]=alpha[k]+2.*M_PI;
    }

  for(k=1;k<vertex->nngh;k++) {
    if(alpha[k]<alpha[k-1]) {
      m1=vertex->ngh[k];
      m2=vertex->ngh[k-1];
      vertex->ngh[k]=m2;
      vertex->ngh[k-1]=m1;
      dummy=alpha[k-1];
      alpha[k-1]=alpha[k];
      alpha[k]=dummy;
      if(k>1) k=k-2;
      }
    }

  delete[] alpha;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

inline int fe_order2(mesh_t & mesh, int n, int from, int rotation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
   Re-order neighbours list of node n so that:

    - if orientation counter-clockwise (direct), first neighbours is right-most one
    - if orientation clockwise, first neighbours is left-most one

 SEE THE NOTES ABOUT
  - #fe_order2 AT THE START  OF THIS FILE
  - AND fe_order2() AT THE END OF THIS FILE
  BEFORE MOVING THIS FUNCTION!

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k, m, m1, m2;
  double t1, t2, p1, p2,  dummy;
  double *alpha;
  double beta, factor;

  switch (rotation) {
    case -1:
      factor = -1.;
      break;

    case +1:
      factor = +1.;
      break;
    }

  t1 = mesh.vertices[n].lon;
  p1 = mesh.vertices[n].lat;

  t2 = mesh.vertices[from].lon;
  if(mesh.type==0) t2 = degree_recale(t2, t1);
  p2 = mesh.vertices[from].lat;
/* *----------------------------------------------------------------------------
  reference direction (angle) */
  beta = factor * atan2(p2 - p1, t2 - t1);
  
  alpha=new double[mesh.vertices[n].nngh];

  for(k = 0; k < mesh.vertices[n].nngh; k++) {
    m = mesh.vertices[n].ngh[k];
    t2 = mesh.vertices[m].lon;
    if(mesh.type==0) t2 = degree_recale(t2, t1);
    p2 = mesh.vertices[m].lat;
    alpha[k] = factor * atan2(p2 - p1, t2 - t1);
    if(m == from)
      alpha[k] = alpha[k] + 2. * M_PI;
    else if(alpha[k] <= beta)
      alpha[k] = alpha[k] + 2. * M_PI;
    }

/* *----------------------------------------------------------------------------
  order neighbours in crescending angle order */
  for(k = 1; k < mesh.vertices[n].nngh; k++) {
    if(alpha[k] < alpha[k - 1]) {
      m1 = mesh.vertices[n].ngh[k];
      m2 = mesh.vertices[n].ngh[k - 1];
      mesh.vertices[n].ngh[k]     = m2;
      mesh.vertices[n].ngh[k - 1] = m1;
      dummy = alpha[k - 1];
      alpha[k - 1] = alpha[k];
      alpha[k] = dummy;
      if(k > 1) {
        k = k - 2;
        }
      }
    }
  
  delete[] alpha;
  
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_oppositeneighbour(mesh_t mesh, int n, int m0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    k,m,m1;
  double t1,t2,p1,p2,alpha[100],dummy=1.e+10;
  double beta;

  t1=mesh.vertices[n].lon;
  p1=mesh.vertices[n].lat;
  t2=mesh.vertices[m0].lon;
  p2=mesh.vertices[m0].lat;

  beta=atan2(p2-p1,t2-t1)+M_PI;

  t1=mesh.vertices[n].lon;
  p1=mesh.vertices[n].lat;
  for(k=0;k<mesh.vertices[n].nngh;k++) {
    m=mesh.vertices[n].ngh[k];
    t2=mesh.vertices[m].lon;
    p2=mesh.vertices[m].lat;
    alpha[k]=fmod(atan2(p2-p1,t2-t1)-beta,2.*M_PI);
    if(alpha[k]> M_PI) alpha[k]-=2.*M_PI;
    if(alpha[k]<-M_PI) alpha[k]+=2.*M_PI;
    if(fabs(alpha[k]) < dummy) {
      dummy=fabs(alpha[k]);
      m1=m;
      }
    }

  return(m1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_anb(const mesh_t & mesh, int a, int b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///  returns index of a in b's neighbour list; -1 if not found
/*----------------------------------------------------------------------------*/
{
  int    k,m;
  const vertex_t *vertex=&mesh.vertices[b];
  
  for(k=0;k<vertex->nngh;k++) {
    m=vertex->ngh[k];
    if(m==a) {
      return(k);
      }
    }
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_quadnodes(mesh_t & mesh, int a, int b, int q[4])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*-----------------------------------------------------------------------------
   */
{
  int k,m;
  int count,found[2];
  double angle;
  
  if(fe_anb(mesh,a,b)==-1) return(-1);
  
  count=0;
  for(k=0;k<mesh.vertices[b].nngh;k++) {
    m=mesh.vertices[b].ngh[k];
    if(fe_anb(mesh,a,m)!=-1) {
      found[count]=m;
      count++;
      if(count>2) return(-1);
      }
    }
  
/*-----------------------------------------------------------------------------
  [a,b] is an orphan edge**/
  if(count==0) return(-1);
  
/*-----------------------------------------------------------------------------
  [a,b] is an ??? edge**/
  if(count>2) return(-1);
  
/*-----------------------------------------------------------------------------
  [a,b] is a boundary edge*/
  if(count==1) return(-1);
  
  angle=fe_angle(mesh,a,b,found[0]);
  if(angle>0) {
    q[0]=a;
    q[1]=found[1];
    q[2]=b;
    q[3]=found[0];
    }
  else {
    q[0]=a;
    q[1]=found[0];
    q[2]=b;
    q[3]=found[1];
    }
  
//   if(q[3]==0) {
//     printf("troubles\n");
//     }
  
   return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_isedgeT(mesh_t & mesh, int a, int b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  find a and b common edge, not element versatile (triangle only) */
{
  int k, i, m, n, n1, n2;

  for(k = 0; k < mesh.vertices[b].nelmts; k++) {
    m = mesh.vertices[b].elmts[k];
    for(i=0;i<3;i++) {
      n=mesh.triangles[m].edges[i];
      n1=mesh.edges[n].extremity[0];
      n2=mesh.edges[n].extremity[1];
      if((a==n1) && (b==n2)) goto finished;
      if((b==n1) && (a==n2)) goto finished;
      }
    }
  //failure:
  return (-1);

  finished:
  return(n);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_isedge(mesh_t & mesh, int a, int b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  seek common edge for a and b vertices, element (triangle/quadrangle) versatile

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
  int k, n, n1, n2;

  if(mesh.vertices[b].edges==0) {
    printf("\nfe_isedge changed, now needs mesh.vertices[b].edges to be initialized before use  (fe_vertex_crosstables02)\n");
    status=fe_vertex_crosstables02(&mesh);
    if(status==-1) TRAP_ERR_EXIT(-1,"%s : fe_vertex_crosstables02 failure, exiting\n",__func__);
    }
  
  for(k = 0; k < mesh.vertices[b].nedges; k++) {
    n = mesh.vertices[b].edges[k];
    n1=mesh.edges[n].extremity[0];
    n2=mesh.edges[n].extremity[1];
    if((a==n1) and (b==n2)) goto finished;
    if((b==n1) and (a==n2)) goto finished;
    }
  //failure:
  return (-1);

  finished:
  return(n);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_distance(const mesh_t & mesh, int n, int m)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double d,t1,t2,p1,p2;

  t1=mesh.vertices[n].lon;
  p1=mesh.vertices[n].lat;
  t2=mesh.vertices[m].lon;
  p2=mesh.vertices[m].lat;

  d=geo_distance(t1,p1,t2,p2);

  return(d);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_distance(const mesh_t & mesh, int n, double t, double p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double d,t1,p1;

  t1=mesh.vertices[n].lon;
  p1=mesh.vertices[n].lat;

  d=geo_distance(t1,p1,t,p);

  return(d);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_distance(const discretisation_t & descriptor, int n, double t, double p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double d,t1,p1;

  t1=descriptor.nodes[n].lon;
  p1=descriptor.nodes[n].lat;

/*------------------------------------------------------------------------------
  warning : geo_distance now return distance in meters */
  d=geo_distance(t1,p1,t,p);

  return(d);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_angle(mesh_t & mesh, int m1, int m2, int m3) 

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  compute angles from longitude nad latitude, vulnerable */
{
  double alpha,t1,t2,t3,p1,p2,p3,p;
  double c=1.0,ctx,cpy,cty,cpx;
  double sine,cosine;

  t1=mesh.vertices[m1].lon;
  p1=mesh.vertices[m1].lat;
  t2=mesh.vertices[m2].lon;
  p2=mesh.vertices[m2].lat;
  if(mesh.type==0) t2=degree_recale(t2,t1);
  t3=mesh.vertices[m3].lon;
  if(mesh.type==0) t3=degree_recale(t3,t1);
  p3=mesh.vertices[m3].lat;

  p=(p1+p2+p3)/3.0;
  if(mesh.type==0) c=cos(p*d2r);

  ctx=  (t2-t1)*c;
  cty=  (t3-t1)*c;
  cpx=  p2-p1;
  cpy=  p3-p1;
  
  cosine=ctx*cty+cpx*cpy;
  sine  =ctx*cpy-cty*cpx;

  alpha=atan2(sine,cosine);

  return(alpha);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_angle_NoRecale(mesh_t & mesh, int m1, int m2, int m3)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  compute angles from presumedly cartesian coordinates that have been substituded
  to actual longitude and latitude in vertex objects, vulnerable */
{
  double alpha,t1,t2,t3,p1,p2,p3;
  double ctx,cpy,cty,cpx;
  double sine,cosine;

  t1=mesh.vertices[m1].lon;
  p1=mesh.vertices[m1].lat;
  t2=mesh.vertices[m2].lon;
  p2=mesh.vertices[m2].lat;
  t3=mesh.vertices[m3].lon;
  p3=mesh.vertices[m3].lat;

  ctx=  t2-t1;
  cty=  t3-t1;
  cpx=  p2-p1;
  cpy=  p3-p1;

  cosine=ctx*cty+cpx*cpy;
  sine  =ctx*cpy-cty*cpx;

  alpha=atan2(sine,cosine);

  return(alpha);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_angle_cartesian(mesh_t & mesh, int m1, int m2, int m3, double *x, double *y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  compute angles from cartesian coordinates given in x and y arrays */
{
  double alpha,t1,t2,t3,p1,p2,p3;
  double ctx,cpy,cty,cpx;
  double sine,cosine;

  t1=x[m1];
  p1=y[m1];
  t2=x[m2];
  p2=y[m2];
  t3=x[m3];
  p3=y[m3];

  ctx=  t2-t1;
  cty=  t3-t1;
  cpx=  p2-p1;
  cpy=  p3-p1;

  cosine=ctx*cty+cpx*cpy;
  sine  =ctx*cpy-cty*cpx;
//  sine=(ctx*cpy-cty*cpx)/sqrt((ctx*ctx+cpx*cpx)*(cty*cty+cpy*cpy));

  alpha=atan2(sine,cosine);

  return(alpha);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_angle_cartesian(mesh_t & mesh, int m1, int m2, int m3)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  compute angles in rads
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double alpha,t1,t2,t3,p1,p2,p3;
  
  t1=mesh.vertices[m1].lon;
  p1=mesh.vertices[m1].lat;
  t2=mesh.vertices[m2].lon;
  p2=mesh.vertices[m2].lat;
  t3=mesh.vertices[m3].lon;
  p3=mesh.vertices[m3].lat;
  
  alpha=geo_angle_radian(t1, p1, t2, p2, t3, p3);

  return(alpha);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_angle_cartesian_obsolete(mesh_t & mesh, int m1, int m2, int m3)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  compute angles in rads from cartesian coordinates locally computed 
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double alpha,t1,t2,t3,p1,p2,p3,x1,x2,x3,y1,y2,y3;
  double ctx,cpy,cty,cpx;
  double sine,cosine;
  
  t1=mesh.vertices[m1].lon;
  p1=mesh.vertices[m1].lat;
  t2=mesh.vertices[m2].lon;
  p2=mesh.vertices[m2].lat;
  t3=mesh.vertices[m3].lon;
  p3=mesh.vertices[m3].lat;
  
  projPJ projection=assign_StereoOblique(p1,t1);
  
  geo_to_projection(projection, p1, t1, &x1, &y1);
  geo_to_projection(projection, p2, t2, &x2, &y2);
  geo_to_projection(projection, p3, t3, &x3, &y3);

  ctx=  x2-x1;
  cty=  x3-x1;
  cpx=  y2-y1;
  cpy=  y3-y1;

  cosine=ctx*cty+cpx*cpy;
  sine  =ctx*cpy-cty*cpx;
  
  alpha =atan2(sine,cosine);
  
  if(alpha<0) {
    printf("%s : warning, negative angle\n", __func__);
    }
  
  free(projection);

  return(alpha);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_vertex_Etable(mesh_t & mesh, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  build the element lists for vertex n
------------------------------------------------------------------------------*/
/**
\param *mesh the mesh
\returns always 0
*/
{
  vector<int> elements;
  int k,m;

  mesh.vertices[n].nelmts=0;

#pragma omp parallel for
  for(m=0; m<mesh.ntriangles; m++) {
    for(int i=0; i<3; i++) {
      int nn =mesh.triangles[m].vertex[i];
      if(nn==n) {
#pragma omp critical(fe_vertex_Etable)
        {
        elements.push_back(m);
        }
        }
      }
    }
  
  mesh.vertices[n].nelmts=elements.size();

  if(mesh.vertices[n].elmts!=0) delete[] mesh.vertices[n].elmts;
  mesh.vertices[n].elmts=new int[mesh.vertices[n].nelmts];

  for(k=0;k<elements.size();k++) mesh.vertices[n].elmts[k]=elements[k];
  
  elements.clear();
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_vertex_element_tables(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// build the element lists for vertices
/**
\param *mesh the mesh
\returns always 0
*/
/*----------------------------------------------------------------------------*/
{
  int nmax;//maximum number of triangles per node
  int i,j,n;//indices for corners, triangles and nodes

/*----------------------------------------------------------------------------*/
  ///It counts the elements

  nmax=0;
  for(n=0; n<mesh->nvtxs; n++) mesh->vertices[n].nelmts=0;

  for(j=0; j<mesh->ntriangles; j++) {
    for(i=0; i<3; i++) {
      n =mesh->triangles[j].vertex[i];
      mesh->vertices[n].nelmts++;
      updatemax(&nmax,mesh->vertices[n].nelmts);
      }
    }

  ///in order to allocate the element lists.
  for(n=0; n<mesh->nvtxs; n++) {
    if(mesh->vertices[n].elmts!=0) delete[] mesh->vertices[n].elmts;
    mesh->vertices[n].elmts=new int[mesh->vertices[n].nelmts];
    }

  ///It then fills the element lists
  for(n=0; n<mesh->nvtxs; n++) mesh->vertices[n].nelmts=0;

  ///triangle by triangle, so <b>the lists are ordered</b>.
  for(j=0; j<mesh->ntriangles; j++) {
    for(i=0; i<3; i++) {
      n =mesh->triangles[j].vertex[i];
      mesh->vertices[n].elmts[mesh->vertices[n].nelmts]=j;
      mesh->vertices[n].nelmts++;
      }
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_setOBCflags(mesh_t & mesh, int n1, int n2, int flag)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,m,n;
  size_t count;
  
  m=0;
  count=0;
  for(l=0;l<mesh.nlimits;l++) {
    for(k=0;k<mesh.limits[l].nvertex;k++) {
      n=mesh.limits[l].vertex[k];
      if(n==n1) {
        do {
          n=mesh.limits[l].edges[k];
          mesh.edges[n].code=flag;
          count++;
          k=(k+1) % mesh.limits[l].nvertex;
          n=mesh.limits[l].vertex[k];
          } while(n!=n2);
        break;
        }
      }
    }

  printf("%d edges received computational flag=%d\n",count,flag);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* KEEP THIS AT THE END OF THE FILE WHERE fe_order2_inline IS DEFINED
KEEPING IN MIND THE CONTENT OF fe_order2 MACRO AT THE START OF THIS FILE! */
#undef fe_order2
int fe_order2(mesh_t & mesh, int n, int from, int rotation){
  return fe_order2_inline(mesh,n,from,rotation);}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
