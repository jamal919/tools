/*******************************************************************************

  T-UGO tools, 2006-2015

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

\brief finite elements geometry edition functions
*/
/*----------------------------------------------------------------------------*/

#include <stdio.h>

#include "tools-structures.h"

#include "fe.h"
#include "geo.h"
#include "zapper.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_SortEdges(mesh_t & mesh, int n, int edge, int rotation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Re-order neighbours list of node n so that:

    - if orientation counter-clockwise (direct), first neighbours is right-most one

    - if orientation clockwise, first neighbours is left-most one
    
  

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int k, m, m1, m2;
  double t1, t2, p1, p2,  dummy;
  double *alpha;
  double beta, factor;
  int from;

  switch (rotation) {
    case -1:
    case +1:
      break;
    default:
      TRAP_ERR_EXIT(ENOEXEC,"programming error : rotation is either +1 or -1, not %d\n",rotation);
    }
  factor=rotation;

  t1 = mesh.vertices[n].lon;
  p1 = mesh.vertices[n].lat;
  
  from=mesh.edges[edge].opposite_vertex(n);

  t2 = mesh.vertices[from].lon;
  if(mesh.type==0) t2 = geo_recale(t2, t1, (double) 180.0);
  p2 = mesh.vertices[from].lat;
  
/*------------------------------------------------------------------------------
  reference direction (angle) */
  beta = factor * atan2(p2 - p1, t2 - t1);
  
  alpha=new double[mesh.vertices[n].nedges];

  for(k = 0; k < mesh.vertices[n].nedges; k++) {
    int e=mesh.vertices[n].edges[k];
    m = mesh.edges[e].opposite_vertex(n);
    t2 = mesh.vertices[m].lon;
    if(mesh.type==0) t2 = geo_recale(t2, t1, (double) 180.0);
    p2 = mesh.vertices[m].lat;
    alpha[k] = factor * atan2(p2 - p1, t2 - t1);
    if(m == from)
      alpha[k] = alpha[k] + 2. * M_PI;
    else if(alpha[k] <= beta)
      alpha[k] = alpha[k] + 2. * M_PI;
    }

/*------------------------------------------------------------------------------
  order neighbours in crescending angle order */
  for(k = 1; k < mesh.vertices[n].nedges; k++) {
    if(alpha[k] < alpha[k - 1]) {
      int e1=mesh.vertices[n].edges[k];
      int e2=mesh.vertices[n].edges[k-1];
      m1 = mesh.edges[e1].opposite_vertex(n);
      m2 = mesh.edges[e2].opposite_vertex(n);
      mesh.vertices[n].ngh[k]     = m2;
      mesh.vertices[n].ngh[k - 1] = m1;
      mesh.vertices[n].edges[k]   = e2;
      mesh.vertices[n].edges[k-1] = e1;
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

  int fe_codetable1(mesh_t * mesh, int RecycleCodes, int stopon_EdgeError, int stopon_PinchError)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*-----------------------------------------------------------------------------

option :

0 : initialize of vertex codes (interior/exterior flag)  from edges table
1 : no initilization, assumes codes already set

-----------------------------------------------------------------------------*/
{
  int k, l, m, n1, n2, n, status;
  int nedges, total, code, rotation;
  int nelts, nndes;
  int ninteriors = 0, nboundaries = 0, nweirds = 0;
  edge_t *edges=NULL;
  double latmin;
  int *count=NULL, *first=NULL, nlmax = 10000;
  int verbose=0;
  vector<size_t> Vsequence;
  vector<int> targeted;

  count = new int[nlmax];
  first = new int[nlmax];

  nelts = mesh->ntriangles;
  nndes = mesh->nvtxs;
  nedges = mesh->nedges;
  edges = mesh->edges;

  switch (RecycleCodes) {
    case 0:
/*-----------------------------------------------------------------------------
      mesh is assumed to be fully constructed*/
      for(n = 0; n < mesh->nvtxs; n++) {
        mesh->vertices[n].code = 0;
        }
      for(n = 0; n < nedges; n++) {
        if(edges[n].nshared == 0) {
          printf("unused edge %d\n", n);
          nweirds++;
          }
        if(edges[n].nshared > 2) {
          printf("overused edge %d shared=%d\n", n, edges[n].nshared);
          nweirds++;
          }
        if(edges[n].nshared == 1)
          nboundaries++;
        if(edges[n].nshared == 2)
          ninteriors++;
        if(edges[n].nshared == 1) {
          n1 = edges[n].extremity[0];
          n2 = edges[n].extremity[1];
          mesh->vertices[n1].code = -1;
          mesh->vertices[n2].code = -1;
          }
        }
      if(verbose) {
        report_edge_table( ninteriors, nboundaries, nweirds, (ninteriors * 2 + nboundaries) / 3, nelts);
        }
      if((nweirds!=0) && (stopon_EdgeError==1)) {
        return(-1);
        }
      break;

    case 1:
/*-----------------------------------------------------------------------------
      mesh is assumed to be partially constructed*/
      for(n = 0; n < mesh->nvtxs; n++) {
        if(mesh->vertices[n].code != 0) {
          mesh->vertices[n].code = -1;
          nboundaries++;
          }
        }
      break;

    default:
      break;
    }
    
//    status=fe_savemesh("test.nei",MESH_FILE_FORMAT_TRIGRID,*mesh);

  for(n = 0; n < mesh->nvtxs; n++) {
    if(mesh->vertices[n].code == -1) targeted.push_back(n);
    }

  total = 0;
  code = 1;
  mesh->nlimits = 0;

  while(total < nboundaries) {
    mesh->nlimits++;
/*-----------------------------------------------------------------------------
    find the first point */
    latmin = -M_PI / 2.0;
    latmin = +1.e+10;
//    for(n = 0; n < mesh->nvtxs; n++) {
    for(int s = 0; s <targeted.size(); s++) {
      n=targeted[s];
      if(mesh->vertices[n].code == -1) {
        if(mesh->vertices[n].lat < latmin) {
          latmin = mesh->vertices[n].lat;
          first[code - 1] = n;
          }
        }
      }
    mesh->vertices[first[code - 1]].code = code;
/*-----------------------------------------------------------------------------
    initialize search */
    if(code == 1) {
/*-----------------------------------------------------------------------------
      scan in direct order (counter-clockwise) */
      rotation = +1;
      status = fe_order(*mesh, first[code - 1], (double) -1., (double) 0., rotation);
      }
    else {
/*-----------------------------------------------------------------------------
      scan in indirect order (clockwise) */
      rotation = -1;
      status = fe_order(*mesh, first[code - 1], (double) -1., (double) -1.e-06, rotation);
      k = 0;
      n = first[code - 1];
      m = mesh->vertices[n].ngh[k];
/*-----------------------------------------------------------------------------
      need to skip first neighbours if they are interior nodes */
      while(mesh->vertices[m].code != -1) {
        k++;
        if(k == mesh->vertices[n].nngh) {
          printf("cannot initiate boundary %d scanning (%d %d)\n", code, n, m);
          break;
          }
        m = mesh->vertices[n].ngh[k];
        }
/*-----------------------------------------------------------------------------
      m is first elligible node, then need to find the most righ-handed one */
      if(k < mesh->vertices[n].nngh-1) {
/*-----------------------------------------------------------------------------
        check if m and next neighbour of n are themselve neighbours */
        while(fe_anb(*mesh, m, mesh->vertices[n].ngh[k + 1]) != -1) {
          k++;
          if(k == mesh->vertices[n].nngh-1) {
/*-----------------------------------------------------------------------------
            we scanned all neighbours of n after m,   */
            printf("#potential trouble (pinched boundary or single triangle island?) at node %d\n", n);
            printf("position lon=%lf lat=%lf\n",mesh->vertices[n].lon,mesh->vertices[n].lat);
            if(stopon_PinchError==1) {
              return(-1);
              }
            break;
            }
          m = mesh->vertices[n].ngh[k];
          }
        }
      n = first[code - 1];
      status = fe_order2(*mesh, n,m, rotation);
      }
    total += 1;
    count[code - 1] = 1;
    n = first[code - 1];
    Vsequence.push_back(n);
    m = -1;
    while(m != first[code - 1]) {
      m = mesh->vertices[n].ngh[0];
      if(m == first[code - 1])
        break;
      if(mesh->vertices[m].code != -1) {
        printf("trouble (pinched boundary?) code %d, count= %d, node= %d %d\n", code, count[code - 1], n, m);
        printf("position lon=%lf lat=%lf\n",mesh->vertices[n].lon,mesh->vertices[n].lat);
        if(stopon_PinchError==1) {
          return(-1);
          }
        }
      status = fe_order2(*mesh, m, n, rotation);
      n = m;
      Vsequence.push_back(n);
      mesh->vertices[n].code = code;
      total++;
      count[code - 1]++;
      }
    if(verbose)
      printf("code %d: total %d count %d (expected total %d)\n", code, total, count[code-1], nboundaries);
    code++;
    }
    
  if(mesh->limits!=0) for(l = 0; l < mesh->nlimits; l++) mesh->limits[l].destroy();
  delete[] mesh->limits;
  
  mesh->limits = new limit_t[mesh->nlimits];
  total=0;
  
  for(l = 0; l < mesh->nlimits; l++) {
    mesh->limits[l].code = l + 1;
    mesh->limits[l].nvertex = count[l];
    mesh->limits[l].vertex = new int[mesh->limits[l].nvertex];
    for(k = 0; k < count[l]; k++) {
      m=Vsequence[total];
      mesh->limits[l].vertex[k] = m;
      total++;
      }
    mesh->limits[l].nedges = count[l];
    mesh->limits[l].edges  = new int[mesh->limits[l].nedges];
    mesh->limits[l].length = 0.0;
    for(k = 0; k < count[l]; k++) {
      n1=mesh->limits[l].vertex[k];
      n2=mesh->limits[l].vertex[(k+1) % mesh->limits[l].nvertex];
      n = fe_isedge(*mesh,n1,n2);
      mesh->limits[l].edges[k]= n;
      mesh->limits[l].length += mesh->edges[n].L;
      mesh->edges[n].code=l+1;
      }
    }

    
  delete[]count;
  delete[]first;
  
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_codetable_obsolete(mesh_t * mesh, int RecycleCodes, int StopOnPinched, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*-----------------------------------------------------------------------------

option :

0 : initialize of vertex codes (interior/exterior flag)  from edges table
1 : no initilization, assumes codes already set

-----------------------------------------------------------------------------*/
{
  int k, l, m, n1, n2, n, status;
  int PinchedDetetected=0,error=0;
  int nedges, total, code, rotation;
  int nelts, nndes;
  int ninteriors = 0, nboundaries = 0, nweirds = 0;
  edge_t *edges;
  double latmin, t1, t2, p1, p2;
  int *count, *first, nlmax = 10000;
  vector<size_t> Vsequence;
  
  count = new int[nlmax];
  first = new int[nlmax];

  nelts = mesh->ntriangles;
  nndes = mesh->nvtxs;
  nedges = mesh->nedges;
  edges = mesh->edges;

  switch (RecycleCodes) {
    case 0:
      for(n = 0; n < mesh->nvtxs; n++) {
        mesh->vertices[n].code = 0;
        }

      for(n = 0; n < nedges; n++) {
        if(edges[n].nshared == 0) {
          printf("unused edge %d\n", n);
          nweirds++;
          }
        if(edges[n].nshared > 2) {
          printf("overused edge %d %d\n", n, edges[n].nshared);
          nweirds++;
          }
        if(edges[n].nshared == 1)
          nboundaries++;
        if(edges[n].nshared == 2)
          ninteriors++;
        if(edges[n].nshared == 1) {
          n1 = edges[n].extremity[0];
          n2 = edges[n].extremity[1];
          mesh->vertices[n1].code = -1;
          mesh->vertices[n2].code = -1;
          }
        }
      if(verbose) {
        report_edge_table( ninteriors, nboundaries, nweirds,(ninteriors * 2 + nboundaries) / 3, nelts);
        }
      break;
      
    case 1:
      for(n = 0; n < mesh->nvtxs; n++) {
        if(mesh->vertices[n].code != 0) {
          mesh->vertices[n].code = -1;
          nboundaries++;
          }
        }
      break;


    default:
      break;
    }

  total = 0;
  code = 1;
  mesh->nlimits = 0;

  while(total < nboundaries) {
//    mesh->nlimits++;
/**-----------------------------------------------------------------------------
    find the first point */
    latmin = -M_PI / 2.0;
    latmin = +1.e+10;
    for(n = 0; n < mesh->nvtxs; n++) {
      if(mesh->vertices[n].code == -1) {
        if(mesh->vertices[n].lat < latmin) {
          latmin = mesh->vertices[n].lat;
          first[code - 1] = n;
          }
        }
      }
    mesh->vertices[first[code - 1]].code = code;
/**-----------------------------------------------------------------------------
    initialize search */
    if(code == 1) {
/*-----------------------------------------------------------------------------
      scan in direct order (counter-clockwise) */
      rotation = +1;
      status =  fe_order(*mesh, first[code - 1], (double) -1., (double) 0.,rotation);
      }
    else {
/*-----------------------------------------------------------------------------
      scan in indirect order (clockwise) + finding first follower of n */
      rotation = -1;
      status = fe_order(*mesh, first[code - 1], (double) -1., (double) -1.e-06, rotation);
      k = 0;
      n = first[code - 1];
      m = mesh->vertices[n].ngh[k];
/*-----------------------------------------------------------------------------
      need to skip first neighbours if they are interior nodes */
      while(mesh->vertices[m].code != -1) {
        k++;
        if(k == mesh->vertices[n].nngh) {
          printf("cannot find next boundary node after %d\n", m);
          break;
          }
        m = mesh->vertices[n].ngh[k];
        }
//      whilefe_anb(*mesh, m, mesh->vertices[n].ngh[k + 1]) != -1) {
      while((k < mesh->vertices[n].nngh-1) && (fe_anb(*mesh, m, mesh->vertices[n].ngh[k + 1])) != -1) {
        k++;
        if(k == mesh->vertices[n].nngh) {
          printf("cannot find next boundary node after %d\n", m);
          break;
          }
        m = mesh->vertices[n].ngh[k];
        }
      t1 = mesh->vertices[first[code - 1]].lon;
      p1 = mesh->vertices[first[code - 1]].lat;
      t2 = mesh->vertices[m].lon;
      t2 = geo_recale(t2, t1, (double) 180.0);
      p2 = mesh->vertices[m].lat;
      status = fe_order(*mesh, n, t2 - t1, p2 - p1, rotation);
      }
/**-----------------------------------------------------------------------------
    do boundary sequence search */
    total += 1;
    count[code - 1] = 1;
    n = first[code - 1];
    Vsequence.push_back(n);
    m = -1;
    while(m != first[code - 1]) {
      m = mesh->vertices[n].ngh[0];
      if(m == first[code - 1])
/*-----------------------------------------------------------------------------
        boundary completed */
        break;
      if(mesh->vertices[m].code != -1) {
/**-----------------------------------------------------------------------------
        so-called pinched boundary */
        printf("trouble ( pinched boundary) in boundary %d: count= %d, current=%d next=%d\n", code, count[code - 1], n, m);
        if(StopOnPinched==1) {
          return(-1);
          }
        PinchedDetetected=1;
        }
//       t1 = mesh->vertices[m].lon;
//       p1 = mesh->vertices[m].lat;
//       t2 = mesh->vertices[n].lon;
//       p2 = mesh->vertices[n].lat;
//       t2 = drecale(t2, t1, (double) 180.0);
//       status=fe_order(*mesh,m,t2-t1,p2-p1,rotation);
      status = fe_order2(*mesh, m, n, rotation);
      n = m;
      Vsequence.push_back(n);
      mesh->vertices[n].code = code;
      total++;
      count[code - 1]++;
      }
    if(verbose)
      printf("code %d: total %d count %d (expected total %d)\n", code, total, count[code - 1], nboundaries);
    code++;
    }

/**----------------------------------------------------------------------------
  mesh limits construction is not safe in pinched boundary error detected */
  if(PinchedDetetected==1) {
    error=-1;
    return(error);
    }
    
  mesh->limits->destroy();
  if(mesh->limits!=0) for(l = 0; l < mesh->nlimits; l++) mesh->limits[l].destroy();
  delete[] mesh->limits;
    
/** ***************************************************
 \bug Jan 23, 2012 :  Clement MAYET : nlimits was not initialized. To be securised
***************************************************** */
  mesh->nlimits=code-1;
  mesh->limits = new limit_t[mesh->nlimits];
  total=0;

/**----------------------------------------------------------------------------
  vulnerable to pinched boundaries, therefore troublesome in metis partition */
  for(l = 0; l < mesh->nlimits; l++) {
    mesh->limits[l].code = l + 1;
    mesh->limits[l].nvertex = count[l];
    mesh->limits[l].vertex = new int[mesh->limits[l].nvertex];
    mesh->limits[l].vertex[0] = first[l];
//     m = first[l];
//     for(k = 1; k < count[l]; k++) {
//       m = mesh->vertices[m].ngh[0];
//       mesh->limits[l].vertex[k] = m;
//       }
    for(k = 0; k < count[l]; k++) {
      m=Vsequence[total];
      mesh->limits[l].vertex[k] = m;
      total++;
      }
    mesh->limits[l].nedges = count[l];
    mesh->limits[l].edges  = new int[mesh->limits[l].nedges];
    mesh->limits[l].length = 0.0;
    for(k = 0; k < count[l]; k++) {
      n1=mesh->limits[l].vertex[k];
      n2=mesh->limits[l].vertex[(k+1) % mesh->limits[l].nvertex];
      n= fe_isedge(*mesh,n1,n2);
      if(n==-1) TRAP_ERR_EXIT(-1, "fe_isedge failed\n");
      mesh->limits[l].edges[k]=n;
      mesh->limits[l].length += mesh->edges[mesh->limits[l].edges[k]].L;
      }
    }

  delete[]count;
  delete[]first;
    
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_codetable2(mesh_t *mesh, int RecycleCodes, int stopon_EdgeError, int stopon_PinchError, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

option :

0 : initialize from edges table (sane mesh)
1 : no initilization (insane mesh), for temporary mesh building steps

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int   k,kk,l,m,n1,n2,n,status;
  int   nedges,total,code,rotation;
  int   nelts,nndes;
  int   ninteriors,nboundaries,nweirds;
//   int   verbose=0;
  edge_t *edges=NULL;
  int nlimits;
  double latmin;
  int *count=NULL, *first=NULL, nlmax=10000;

  exitIfNull(
    count=(int *) malloc(nlmax*sizeof(int))
    );
  exitIfNull(
    first=(int *) malloc(nlmax*sizeof(int))
    );

start:

  ninteriors=0;
  nboundaries=0;
  nweirds=0;

  nelts =mesh->ntriangles;
  nndes =mesh->nvtxs;
  nedges=mesh->nedges;
  edges =mesh->edges;
  
  bool fail=false;
  
  if(edges==0){
    fail=true;
    STDERR_BASE_LINE_FUNC("no edge table found: YOU MUST CALL fe_edgetable BEFORE %s\n",__func__);
    }

  if(mesh->triangles==0 and mesh->quadrangles==0){
    fail=true;
    STDERR_BASE_LINE_FUNC("no element table found: YOU MUST CALL fe_list (or something else ?) BEFORE %s\n",__func__);
    }
  
  if(fail) TRAP_ERR_EXIT(ENOEXEC,"%s:exiting\n",__func__);

  switch (RecycleCodes) {
    case 0:
/*------------------------------------------------------------------------------
      rebuild boundary codes informations from mesh geometry*/
      for (n=0; n<mesh->nvtxs; n++) {
        mesh->vertices[n].code=0;
        }
      for (n=0;n<nedges;n++) {
        bool weird=false;
        if (edges[n].nshared==0) {
          printf("unused edge %d: ",n);
          weird=true;
          }
        if (edges[n].nshared >2) {
          printf("overused edge %d nshared=%d\n",n,edges[n].nshared);
          for(k=0;k<edges[n].nshared;k++) printf("element %d=%d\n",k,edges[n].shared[k]);
          weird=true;
          }
        if (edges[n].nshared==1) nboundaries++;
        if (edges[n].nshared==2) ninteriors++;
        if (edges[n].nshared==1) {
          n1=edges[n].extremity[0];
          n2=edges[n].extremity[1];
          mesh->vertices[n1].code=-1;
          mesh->vertices[n2].code=-1;
          }
        if (weird){
          nweirds++;
          n1=edges[n].extremity[0];
          n2=edges[n].extremity[1];
          printf("nodes %d (lon %g lat %g) and %d (lon %g lat %g)\n",
            n1,mesh->vertices[n1].lon,mesh->vertices[n1].lat,
            n2,mesh->vertices[n2].lon,mesh->vertices[n2].lat);
          }
        }
      if(verbose==1) {
        report_edge_table(ninteriors,nboundaries,nweirds,(ninteriors*2+nboundaries)/3,nelts);
        }
      if(nweirds!=0) {
        report_edge_table(ninteriors,nboundaries,nweirds,(ninteriors*2+nboundaries)/3,nelts);
        goto error;
        }
      break;

    case 1:
/*------------------------------------------------------------------------------
      recycle existing boundary codes informations*/
      for (n=0; n<mesh->nvtxs; n++) {
        if(mesh->vertices[n].code!=0) {
          mesh->vertices[n].code=-1;
          nboundaries++;
          }
        }
      break;

    default:
      break;
    }

  total=0;
  code=1;
//   mesh->nlimits=0;
  nlimits=0;

  while (total<nboundaries) {
    nlimits++;
/*-----------------------------------------------------------------------------
    find the first point */
    latmin=1e+10;
    for (n=0; n<mesh->nvtxs; n++) {
      const vertex_t *vertex=&mesh->vertices[n];
      
      if(vertex->code==-1) {
        if(vertex->lat<latmin) {
          latmin=vertex->lat;
          first[code-1]=n;
          }
        }
      }
    mesh->vertices[first[code-1]].code=code;

    if(code==1) {
/*-----------------------------------------------------------------------------
      test whether really outer boundary */
      
      const vertex_t *vertex=&mesh->vertices[first[code-1]];
      
      rotation=+1;
      
      for(k=0;k<vertex->nngh;k++){
        m=vertex->ngh[k];
        if(mesh->vertices[m].lat<latmin){
          /* island */
          rotation=-1;
          break;
          }
        }
      
      }
    else {
/*-----------------------------------------------------------------------------
      inner boudary, i.e. island */
      rotation=-1;
      }
    
/*------------------------------------------------------------------------------
    initialize search */
    if(rotation>0) {
/*-----------------------------------------------------------------------------
      scan in direct order (counter-clockwise) */
      status=fe_order(*mesh,first[code-1],(double) -1., (double) 0.,rotation);
      }
    else {
/*-----------------------------------------------------------------------------
      scan in indirect order (clockwise) */
      status=fe_order(*mesh,first[code-1],(double) -1., (double) -1.e-06,rotation);
      n = first[code - 1];
      k=0;
      m=mesh->vertices[n].ngh[k];
      
/* *----------------------------------------------------------------------------
      need to skip first neighbours if interior nodes (should not happen ?!?)*/
      while(mesh->vertices[m].code!=-1) {
        k++;
        if(k==mesh->vertices[n].nngh) {
          printf ("trouble with node %d, cannot find boundary neighbour\n",m);
          goto error;
          }
        m=mesh->vertices[n].ngh[k];
        }
        
/*------------------------------------------------------------------------------
      need to detect pinched boundaries: */
      if(k < mesh->vertices[n].nngh-1) {
        while(fe_anb(*mesh, m, mesh->vertices[n].ngh[k + 1]) != -1) {
          k++;
          m = mesh->vertices[n].ngh[k];
          if(k == mesh->vertices[n].nngh-1) {
            printf("trouble (pinched boundary or non-connex mesh) at vertex %d\n", n);
            printf("position lon=%lf lat=%lf\n",mesh->vertices[n].lon,mesh->vertices[n].lat);
            if(stopon_PinchError==1) {
              return(-n-1);
              }
            break;
            }
          }
        }
      status = fe_order2(*mesh, n, m, rotation);
      }
    total+=1;
    count[code-1]=1;
    n=first[code-1];
    m=-1;
    while (m!=first[code-1]) {
      m=mesh->vertices[n].ngh[0];
/*------------------------------------------------------------------------------
      check if boundary is completed*/
      if(m==first[code-1]) break;
      if(mesh->vertices[m].code!=-1) {
/*------------------------------------------------------------------------------
        check if boundary is completed*/
        printf ("trouble with node %d, code=%d, code of first neighbour %d is not -1 but %d\n",n,code,m,mesh->vertices[m].code);
        printf("position lon=%lf lat=%lf\n",mesh->vertices[n].lon,mesh->vertices[n].lat);
        if(stopon_PinchError==1) {
          return(-n-1);
          }
        status=fe_cleavelimit(mesh, m);
        if(status!=0) goto error;
        else goto start;
        }
      status = fe_order2(*mesh, m, n, rotation);
      if(RecycleCodes==0) {
        for(k=0;k<mesh->vertices[m].nngh-1;k++) {
          if(fe_anb(*mesh, mesh->vertices[m].ngh[k], mesh->vertices[m].ngh[k+1]) == -1) {
            printf("#potential trouble (pinched boundary) at vertices current=%d next=%d\n", n, m);
            printf("position lon=%lf lat=%lf\n",mesh->vertices[n].lon,mesh->vertices[n].lat);
            for(kk=0;kk<mesh->vertices[m].nngh;kk++) printf("%dth neighbour is %d\n",kk+1,mesh->vertices[m].ngh[kk]);
            if(stopon_PinchError==1) {
              return(-m-1);
              }
            status=fe_cleavelimit(mesh, m);
            if(status!=0) {
              goto error;
              }
            else goto start;
            break;
            }
          }
        }
      n=m;
      mesh->vertices[n].code=code;
      total++;
      count[code-1]++;
      }
//    printf ("code %d: total %d count %d (expected total %d)\n",code, total,count[code-1],nboundaries);
    code++;
    }

  if(mesh->limits!=0) for(l = 0; l < mesh->nlimits; l++) mesh->limits[l].destroy();
  delete[] mesh->limits;
  
  mesh->nlimits=nlimits;
  exitIfNull(mesh->limits=new limit_t[mesh->nlimits]);

  for(l=0;l<mesh->nlimits;l++) {
    mesh->limits[l].code=l+1;
    mesh->limits[l].nvertex=count[l];
    exitIfNull(mesh->limits[l].vertex=new int[mesh->limits[l].nvertex]);
    mesh->limits[l].vertex[0]=first[l];
    m=first[l];
    for(k=1;k<count[l];k++) {
      m=mesh->vertices[m].ngh[0];
      mesh->limits[l].vertex[k]=m;
      }
    }

  if(RecycleCodes==1) goto end;

  if(mesh->vertices[0].nelmts==-1) {
    status=fe_vertex_element_tables(mesh);
    }
    
  for(l=0;l<mesh->nlimits;l++) {
    mesh->limits[l].nedges = count[l];
    mesh->limits[l].edges  = new int[mesh->limits[l].nedges];
    mesh->limits[l].length = 0.0;
    for(k = 0; k < count[l]; k++) {
      n1=mesh->limits[l].vertex[k];
      n2=mesh->limits[l].vertex[(k+1) % mesh->limits[l].nvertex];
      mesh->limits[l].edges[k]= fe_isedge(*mesh,n1,n2);
      mesh->limits[l].length += mesh->edges[mesh->limits[l].edges[k]].L;
      mesh->edges[mesh->limits[l].edges[k]].code=l+1;
      }
    }

end:
  free(count);
  free(first);

  return(0);
error:
  free(count);
  free(first);

  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_initVcodes_01(mesh_t & mesh, int & nboundaries, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  set geometric boundary codes, algorithm based on edges table
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int n1, n2, n, status;
  int total,code;
  int ninteriors,nweirds,nelements;
  
  ninteriors=0;
  nboundaries=0;
  nweirds=0;

  bool fail=false;
  
  if(mesh.edges==0){
    fail=true;
    STDERR_BASE_LINE_FUNC("no edge table found: YOU MUST CALL fe_edgetable BEFORE %s\n",__func__);
    }

  if(mesh.triangles==0 and mesh.quadrangles==0){
    fail=true;
    STDERR_BASE_LINE_FUNC("no element table found: YOU MUST CALL fe_list (or something else ?) BEFORE %s\n",__func__);
    }
  
  if(fail) TRAP_ERR_EXIT(ENOEXEC,"%s:exiting\n",__func__);
  
  if(mesh.nelements!=0) {
    nelements=mesh.nelements;
    }
  else if(mesh.ntriangles!=0) {
    nelements=mesh.ntriangles;
    }
  else if(mesh.nquadrangles!=0) {
    nelements=mesh.nquadrangles;
    }
    

/*------------------------------------------------------------------------------
  rebuild boundary codes informations from mesh geometry*/
  for (n=0; n<mesh.nvtxs; n++) {
    mesh.vertices[n].code=0;
    }
  for (n=0;n<mesh.nedges;n++) {
    bool weird=false;
    if (mesh.edges[n].nshared==0) {
      printf("unused edge %d: ",n);
      weird=true;
      }
    if (mesh.edges[n].nshared >2) {
      printf("overused edge %d nshared=%d\n",n,mesh.edges[n].nshared);
      for(int k=0;k<mesh.edges[n].nshared;k++) printf("element %d=%d\n",k,mesh.edges[n].shared[k]);
      weird=true;
      }
    if (mesh.edges[n].nshared==1) nboundaries++;
    if (mesh.edges[n].nshared==2) ninteriors++;
    if (mesh.edges[n].nshared==1) {
      n1=mesh.edges[n].extremity[0];
      n2=mesh.edges[n].extremity[1];
      mesh.vertices[n1].code=-1;
      mesh.vertices[n2].code=-1;
      }
    if (weird){
      nweirds++;
      n1=mesh.edges[n].extremity[0];
      n2=mesh.edges[n].extremity[1];
      printf("nodes %d (lon %g lat %g) and %d (lon %g lat %g)\n",
              n1,mesh.vertices[n1].lon,mesh.vertices[n1].lat,
              n2,mesh.vertices[n2].lon,mesh.vertices[n2].lat);
      }
    }
    
  if(verbose==1) {
    report_edge_table(ninteriors,nboundaries,nweirds,(ninteriors*2+nboundaries)/3,nelements);
    }
    
  if(nweirds!=0) {
    report_edge_table(ninteriors,nboundaries,nweirds,(ninteriors*2+nboundaries)/3,nelements);
    status=-1;
    }
  else status=0;

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_initVcodes_02(mesh_t & mesh, int & nboundaries, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  set geometric boundary codes, algorithm based on vertex table
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int m, n1, n2, n, status;
  int total,code;
  int ninteriors,nweirds;
  int *nelts;
  
  ninteriors=0;
  nboundaries=0;
  nweirds=0;
  
  nelts=new int[mesh.nvtxs];
  for(n=0;n<mesh.nvtxs;n++) nelts[n]=0;
  
  for(m=0;m<mesh.ntriangles;m++) {
    for (int k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
      nelts[n]++;
      }
    }

/*------------------------------------------------------------------------------
  rebuild boundary codes informations from mesh geometry*/
  for (n=0; n<mesh.nvtxs; n++) {
    mesh.vertices[n].code=0;
    if(nelts[n]!=mesh.vertices[n].nngh) {
      mesh.vertices[n].code=-1;
      nboundaries++;
      }
    }
    
  if(nweirds!=0) {
    status=-1;
    }
  else status=0;
  
  delete[] nelts;

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_initVcodes(mesh_t & mesh, int & nboundaries, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int n1, n2, n, status;
  int total,code;
  int ninteriors,nweirds;
  
  ninteriors=0;
  nboundaries=0;
  nweirds=0;

  bool fail=false;
  
  if(mesh.triangles==0 and mesh.quadrangles==0){
    fail=true;
    STDERR_BASE_LINE_FUNC("no element table found: YOU MUST CALL fe_list (or something else ?) BEFORE %s\n",__func__);
    }
  
  if(fail) TRAP_ERR_EXIT(ENOEXEC,"%s:exiting\n",__func__);
  
  if(mesh.edges==0){
    status=fe_initVcodes_02(mesh, nboundaries, verbose);
    }
  else {
    status=fe_initVcodes_01(mesh, nboundaries, verbose);
    }
  return(status);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_initEcodes(mesh_t & mesh, int & nboundaries, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int n1, n2, n, status;
  int total,code;
  int ninteriors,nweirds;
  
  ninteriors=0;
  nboundaries=0;
  nweirds=0;

  for (n=0; n<mesh.nedges; n++) {
    if(mesh.edges[n].nshared<=0) {
      nweirds++;
      }
    if(mesh.edges[n].nshared==1) {
      mesh.edges[n].code=MESH_UNDEFINED_EDGE;
      nboundaries++;
      }
    if(mesh.edges[n].nshared==2) {
      mesh.edges[n].code=MESH_INTERIOR_EDGE;
      ninteriors++;
      }
    if(mesh.edges[n].nshared>2) {
      nweirds++;
      }
    }
  
  if(nweirds==0) status=0;
  else status=-1;
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *get_Vboundarylist(mesh_t & mesh, int nboundaries, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  size_t m=0;
  int *list=0;
  
  if(nboundaries==0) return(list);
  
  list=new int[nboundaries];

  for (n=0; n<mesh.nvtxs; n++) {
    if(mesh.vertices[n].code==-1) {
      list[m]=n;
      m++;
      }
    }
    
  return(list);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *get_Eboundarylist(mesh_t & mesh, int nboundaries, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  size_t m=0;
  int *list=0;
  
  if(nboundaries==0) return(list);
  
  list=new int[nboundaries];

  for (n=0; n<mesh.nedges; n++) {
    if(mesh.edges[n].code==MESH_UNDEFINED_EDGE) {
      list[m]=n;
      m++;
      }
    }
    
  return(list);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int get_Estart(mesh_t & mesh, int *boundaries, int nboundaries, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  select a suitable edge to start mesh limits detection

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int n1,n2,nstart,n,status;
  int first=-1;
  double latmin=+1.0e+10;
    
  for (int k=0; k<nboundaries; k++) {
    n=boundaries[k];
    n1=mesh.edges[n].extremity[0];
    n2=mesh.edges[n].extremity[1];
    if(mesh.vertices[n1].lat<latmin) {
      first=n;
      latmin=mesh.vertices[n1].lat;
      nstart=n1;
      }
    if(mesh.vertices[n2].lat<latmin) {
      first=n;
      latmin=mesh.vertices[n2].lat;
      nstart=n2;
      }
    } 
  
  return(first);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_SetLimits(mesh_t & mesh, int *sequence, vector<int> & first, vector<int> & last, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int l, m, n, n1, n2, status;
  
  for(l = 0; l < mesh.nlimits; l++) {
    mesh.limits[l].destroy();
    }
  deletep(&mesh.limits);
  
  mesh.nlimits=first.size();
  mesh.limits = new limit_t[mesh.nlimits];

  for(l = 0; l < mesh.nlimits; l++) {
    mesh.limits[l].code = l + 1;
    mesh.limits[l].nvertex = last[l]-first[l]+1;
    deletep(&mesh.limits[l].vertex);
    mesh.limits[l].vertex = new int[mesh.limits[l].nvertex];
    for(int k = first[l]; k <= last[l]; k++) {
      m=sequence[k];
      mesh.limits[l].vertex[k-first[l]] = m;
      }
    mesh.limits[l].nedges = last[l]-first[l]+1;
    deletep(&mesh.limits[l].edges);
    mesh.limits[l].edges  = new int[mesh.limits[l].nedges];
    mesh.limits[l].length = 0.0;
    for(int k = first[l]; k <= last[l]; k++) {
      n1=mesh.limits[l].vertex[k-first[l]];
      n2=mesh.limits[l].vertex[(k+1-first[l]) % mesh.limits[l].nvertex];
      n= fe_isedge(mesh,n1,n2);
      if(n==-1) check_error(-1, "fe_isedge failed", __LINE__, __FILE__, 1);
      mesh.limits[l].edges[k-first[l]]=n;
      mesh.limits[l].length += mesh.edges[n].L;
      }
    }
  
  return(status);
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_SortEdges__(mesh_t & mesh, int n, int edge, int rotation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Re-order neighbours list of node n so that:

    - if orientation counter-clockwise (direct), first neighbours is right-most one

    - if orientation clockwise, first neighbours is left-most one
    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k, m, m1, m2;
  double t1, t2, p1, p2,  dummy;
  double *alpha;
  double beta, factor;
  int from;

  switch (rotation) {
    case -1:
    case +1:
      break;
    default:
      TRAP_ERR_EXIT(ENOEXEC,"programming error : rotation is either +1 or -1, not %d\n",rotation);
    }
  factor=rotation;

  t1 = mesh.vertices[n].lon;
  p1 = mesh.vertices[n].lat;
  
  from=mesh.edges[edge].opposite_vertex(n);

  t2 = mesh.vertices[from].lon;
  if(mesh.type==0) t2 = geo_recale(t2, t1, (double) 180.0);
  p2 = mesh.vertices[from].lat;
  
/*------------------------------------------------------------------------------
  reference direction (angle) */
  beta = factor * atan2(p2 - p1, t2 - t1);
  
  alpha=new double[mesh.vertices[n].nedges];

  for(k = 0; k < mesh.vertices[n].nedges; k++) {
    int e=mesh.vertices[n].edges[k];
    m = mesh.edges[e].opposite_vertex(n);
    t2 = mesh.vertices[m].lon;
    if(mesh.type==0) t2 = geo_recale(t2, t1, (double) 180.0);
    p2 = mesh.vertices[m].lat;
    alpha[k] = factor * atan2(p2 - p1, t2 - t1);
    if(m == from)
      alpha[k] = alpha[k] + 2. * M_PI;
    else if(alpha[k] <= beta)
      alpha[k] = alpha[k] + 2. * M_PI;
    }

/*------------------------------------------------------------------------------
  sort neighbours in ascending angle order */
//   for(k = 1; k < mesh.vertices[n].nedges; k++) {
  for(k = 1; k < mesh.vertices[n].nngh; k++) {
    if(alpha[k] < alpha[k - 1]) {
      int e1=mesh.vertices[n].edges[k];
      int e2=mesh.vertices[n].edges[k-1];
      if(k>=mesh.vertices[n].nngh) {
        printf("anomaly\n");
        }
      if(k<=0) {
        printf("anomaly\n");
        }
      mesh.vertices[n].edges[k]   = e2;
      mesh.vertices[n].edges[k-1] = e1;
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

  int fe_codetable(mesh_t & mesh, int RecycleCodes, int SetLimits, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  attempt for a robust boundary code parser
  
  based on edge topology
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int v1, v2, n,next,status;
  int count,rotation;
  int ninteriors,nboundaries,nweirds;
  int nlimits;
  double latmin;
  vector<int> Efirst,Elast;
  vector<int> Vfirst,Vlast;
  int *boundaries=0;
  int *Esequence=0, *Vsequence=0;
  bool fail=false;
  
  if(mesh.nvtxs==0) return(0);
  
  ninteriors=0;
  nboundaries=0;
  nweirds=0;
    
  if(mesh.triangles==0 and mesh.quadrangles==0){
    fail=true;
    STDERR_BASE_LINE_FUNC("no element table found: YOU MUST CALL fe_list (or something else ?) BEFORE %s\n",__func__);
    }
  
  if(fail) TRAP_ERR_EXIT(ENOEXEC,"%s:exiting\n",__func__);

  switch (RecycleCodes) {
    case 0:
/*------------------------------------------------------------------------------
      build boundary codes informations from mesh geometry*/
      status=fe_initVcodes(mesh, nboundaries, verbose);
      status=fe_initEcodes(mesh, nboundaries, verbose);
      if(status!=0) {
        goto error;
        }
      break;

    case 1:
/*------------------------------------------------------------------------------
      recycle existing boundary codes informations*/
      for (n=0; n<mesh.nedges; n++) {
        if(mesh.edges[n].code!=MESH_INTERIOR_EDGE) {
          mesh.edges[n].code=MESH_UNDEFINED_EDGE;
          nboundaries++;
          }
        }
      for (n=0; n<mesh.nvtxs; n++) {
        if(mesh.vertices[n].code!=0) {
          mesh.vertices[n].code=-1;
          }
        }
      break;

    default:
      break;
    }
    
  if(mesh.vertices[0].elmts!=0)
  for(n = 0; n < mesh.nvtxs; n++) {
    int misfit=mesh.vertices[n].nngh-mesh.vertices[n].nelmts;
    switch (misfit) {
      case 0:
        if(mesh.vertices[n].code!=0)  TRAP_ERR_EXIT(-1,"%s : anomaly #1 n=%d nngh=%d nelmts=%d\n",__func__,n,mesh.vertices[n].nngh,mesh.vertices[n].nelmts);
        break;
      case 1:
        if(mesh.vertices[n].code!=-1) TRAP_ERR_EXIT(-1,"%s : anomaly #2 n=%d nngh=%d nelmts=%d\n",__func__,n,mesh.vertices[n].nngh,mesh.vertices[n].nelmts);
        break;
      default:
//         if(verbose) printf("cpu %3d : pinched boundary n=%d nngh=%d nelmts=%d\n",gCPU_ID, n,mesh.vertices[n].nngh,mesh.vertices[n].nelmts);
        if(verbose) printf("%s : pinched boundary n=%d nngh=%d nelmts=%d\n",__func__, n,mesh.vertices[n].nngh,mesh.vertices[n].nelmts);
//         TRAP_ERR_EXIT(-1,"anomaly #3 n=%d nngh=%d nelmts=%d\n",n,mesh.vertices[n].nngh,mesh.vertices[n].nelmts);
        break;
      }
    }

  if(mesh.vertices[0].edges==0) {
    TRAP_ERR_EXIT(ENOEXEC,"%s: vertex edge's table not set, exiting. call fe_edgetable to fix issue\n",__func__);
    }

  count=0;
  nlimits=0;
  
  Esequence=new int[nboundaries];
  Vsequence=new int[nboundaries];
 
  while (count<nboundaries) {
    boundaries=get_Eboundarylist(mesh, nboundaries-count, verbose);
/*-----------------------------------------------------------------------------
    find the first point */
    n=get_Estart(mesh, boundaries, nboundaries-count, verbose);
    v1=mesh.edges[n].extremity[0];
    v2=mesh.edges[n].extremity[1];
    Efirst.push_back(count);
    Vfirst.push_back(count);
    Esequence[count]=n; Vsequence[count]=v1; count++;
    rotation = +1;
    status=fe_SortEdges__(mesh, v2, n,rotation);
    next=mesh.vertices[v2].edges[0];
    while(next!=Esequence[Efirst[nlimits]]) {
      n=next;
      if(mesh.edges[next].code!=-1) {
        TRAP_ERR_EXIT(-1,"%s: not boundary edge\n",__func__);
        }
      Esequence[count]=n; Vsequence[count]=v2; count++;
      if(count==nboundaries) {
        break;
        TRAP_ERR_EXIT(-1,"anomaly : n=%d next=%d v2=%d\n",n,next,v2);
        }
      v2=mesh.edges[n].extremity[1];
      status=fe_SortEdges__(mesh, v2, n,rotation);
      next=mesh.vertices[v2].edges[0];
      }
    Elast.push_back(count-1);
    Vlast.push_back(count-1);
    for(int k=Efirst[nlimits];k<=Elast[nlimits];k++) {
      n=Esequence[k];
      mesh.edges[n].code=nlimits+1;
      v1=mesh.edges[n].extremity[0];
      v2=mesh.edges[n].extremity[1];
      mesh.vertices[v1].code=nlimits+1;
      mesh.vertices[v2].code=nlimits+1;
      }
    nlimits++;
    deletep(&boundaries);
    }

  if(SetLimits==1) status=fe_SetLimits(mesh, Vsequence, Vfirst, Vlast, verbose);
  
//   status=fe_savemesh("processed.nei",MESH_FILE_FORMAT_TRIGRID, mesh);
  
  deletep(&Esequence);
  deletep(&Vsequence);
    
end:

  return(0);
  
error:

  return(-1);
}
