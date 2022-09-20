
/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2012

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France

***************************************************************************/

#include  "tugo-prototypes.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void treat_boundarycode_CQP0(mesh_t& mesh, specification_obsolete_t *spec, int code)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l;
  int m, mm, n, n1, n2;
  int nitems, found;
  int nedges, count = 0;
  int width;
  edge_t *edges;
  int *tmp;

  edges  = mesh.edges;
  nedges = mesh.nedges;

  tmp = new int[mesh.nquadrangles];
  for(n = 0; n < mesh.nquadrangles; n++)
    tmp[n] = -1;

  spec->nedges = 0;
  spec->nnodes = 0;

/**----------------------------------------------------------------------
  count all vertices and edges concerned with code*/
  for(n = 0; n < nedges; n++) {
    m = edges[n].shared[0];
    if(edges[n].code == code) {
      if(edges[n].nshared!=1) {
        printf("#boundary anomaly %d, edges %d, nodes %d\n", code, n, m);
        }
      spec->nedges++;
      if(tmp[m] == -1) {
        tmp[m] = spec->nnodes;
        spec->nnodes++;
        }
      else {
        printf("#boundary redundancy %d, edges %d, nodes %d\n", code, n, m);
        printf("lon=%lf lat=%lf\n",edges[n].lon/d2r,edges[n].lat/d2r);
        }
      }
    }

  spec->edges  = new int[spec->nedges];
  spec->nodes  = new int[spec->nnodes];
  spec->table  = new int[spec->nnodes];
  spec->r_edge = new int[spec->nnodes];
  spec->l_edge = new int[spec->nnodes];
  spec->r_node = new int[spec->nnodes];
  spec->l_node = new int[spec->nnodes];

  for(n = 0; n < spec->nedges; n++) {
    spec->edges[n] = -1;
    }

  for(n = 0; n < spec->nnodes; n++) {
    spec->r_node[n] = -1;
    spec->l_node[n] = -1;
    spec->r_edge[n] = -1;
    spec->l_edge[n] = -1;
    }

  spec->nedges = 0;
  spec->nnodes = 0;
  width = 0;
  for(n = 0; n < mesh.nquadrangles; n++)
    tmp[n] = -1;

  for(n = 0; n < nedges; n++) {
    m = edges[n].shared[0];
    if(edges[n].code == code) {
      if(tmp[m] == -1) {
        spec->nodes[spec->nnodes] = m;
        tmp[m] = spec->nnodes;
        spec->nnodes++;
/**----------------------------------------------------------------------------
        add neighbours elements*/
//         for(k=0;k<mesh.triangles[m].nngh;k++) {
//           mm=mesh.triangles[m].ngh[k];
//           if(tmp[mm] == -1) {
//             spec->nodes[spec->nnodes] = mm;
//             spec->nnodes++;
//             tmp[mm] = spec->nnodes;
//             }
//           }
        }
      spec->edges[spec->nedges] = n;
      spec->nedges++;
      }
    }

  zaparr(tmp);

  spec->hbw = width;

  printf("#boundary specification %d, edges %d, nodes %d\n", code, spec->nedges, spec->nnodes);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void treat_boundarycode_CQN1(mesh_t& mesh, specification_obsolete_t *spec, int code)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j, k, l;
  int m, n, n1, n2, n3;
  int nitems, found;
  int count = 0;
  int width;
  int *tmp;
  int eindex[4][2];
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check : 

  Note:

  20/08/2012 : for more generic codeing...

    to be added to element structure:
  
      - number of edge/elements
    
    to be added to discretisation structure:
  
      - numbering of nodes along edges
    

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  discretisation_t descriptor=mesh.CQN1descriptor;
  
  if(descriptor.nnodes==0) check_error(-1, "discretisation not initialized", __LINE__, __FILE__, 1);

  eindex[0][0]=0;
  eindex[0][1]=1;

  eindex[1][0]=1;
  eindex[1][1]=2;

  eindex[2][0]=2;
  eindex[2][1]=3;

  eindex[3][0]=3;
  eindex[3][1]=0;

  tmp = new int[descriptor.nnodes];
  for(n = 0; n < descriptor.nnodes; n++)
    tmp[n] = -1;

  spec->nedges = 0;
  spec->nnodes = 0;

/**----------------------------------------------------------------------
  count all vertices and edges concerned with code*/
  for(m = 0; m < mesh.nquadrangles; m++) {
    for(j=0;j<4;j++) {
      n=mesh.quadrangles[m].edges[j];
      if(mesh.edges[n].code == code) {
        spec->nedges++;
        if(tmp[n] == -1) {
          spec->nnodes++;
          tmp[n] = spec->nnodes;
          }
        }
      }
    }

  spec->edges  = new int[spec->nedges];
  spec->nodes  = new int[spec->nnodes];
  spec->table  = new int[spec->nnodes];
  spec->r_edge = new int[spec->nnodes];
  spec->l_edge = new int[spec->nnodes];
  spec->r_node = new int[spec->nnodes];
  spec->l_node = new int[spec->nnodes];

  for(n = 0; n < spec->nedges; n++) {
    spec->edges[n] = -1;
    }

  for(n = 0; n < spec->nnodes; n++) {
    spec->r_node[n] = -1;
    spec->l_node[n] = -1;
    spec->r_edge[n] = -1;
    spec->l_edge[n] = -1;
    }

  spec->nedges = 0;
  spec->nnodes = 0;
  width = 0;

  for(n = 0; n < descriptor.nnodes; n++)
    tmp[n] = -1;

  for(m = 0; m < mesh.nquadrangles; m++) {
    for(j=0;j<4;j++) {
      n=mesh.quadrangles[m].edges[j];
      if(mesh.edges[n].code == code) {
        spec->edges[spec->nedges] = n;
        spec->nedges++;
        if(tmp[n] == -1) {
          spec->nodes[spec->nnodes] = n;
          tmp[n] = spec->nnodes;
          spec->nnodes++;
          }
        }
      }
    }

//   for(n = 0; n < spec->nnodes; n++) {
//     width = MAX(width, abs(n- spec->r_node[n]));
//     width = MAX(width, abs(n- spec->l_node[n]));
//     }
//   spec->hbw = width;

  zaparr(tmp);

  printf("#boundary specification %d, edges %d, nodes %d\n", code, spec->nedges, spec->nnodes);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void treat_boundarycode(mesh_t& mesh, specification_obsolete_t *spec, int discretisation, int code)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l;
  int m, n, n1, n2;
  int nitems, found;
  int nedges, count = 0;
  int width;
  edge_t *edges;
  int *tmp;
  discretisation_t descriptor;

  spec->discretisation=discretisation;
  spec->code=code;

  switch (discretisation) {
    case CQP0:
      treat_boundarycode_CQP0(mesh, spec, code);
      break;

    case CQP1:
      treat_boundarycode_CQN1(mesh, spec, code);
      break;

    case CQN1:
      treat_boundarycode_CQN1(mesh, spec, code);
      break;

    default:
      check_error(-1, "illegal discretisation", __LINE__, __FILE__, 1);
      break;
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void process_boundarycode(mesh_t & mesh, int level, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------

  treat boundary codes for boundary condition processing

----------------------------------------------------------------------*/
{
  int k, l, e, status;
  int dum1, dum2, m, n, nd1, nd2, code, n1, n2;
  int nitems, found;
  int count = 0;
  char *msg;
  int *tmp;

/**----------------------------------------------------------------------
  read .bel file and gather mesh boundary description */
  fe_boundarycode(BelFile,mesh,level);

  for(n = 0; n < mesh.nedges; n++) {
    n1 = mesh.edges[n].extremity[0];
    n2 = mesh.edges[n].extremity[1];
    if((mesh.edges[n].nshared == 1) && (mesh.edges[n].code == MESH_UNDEFINED_NODE)) {
      printf("boundary error #1: edge=%d nodes= %d %d code=%d\n", n, n1, n2, mesh.edges[n].code);
      }
    if((mesh.edges[n].nshared == 2) && (mesh.edges[n].code != MESH_UNDEFINED_NODE)) {
      printf("boundary error #2: edge=%d nodes= %d %d code=%d\n", n, n1, n2, mesh.edges[n].code);
      }
    }

/**----------------------------------------------------------------------
  Elevation boundary conditions */
/*----------------------------------------------------------------------
  set open limits and river limits with similar code*/
  for(n = 0; n < mesh.nedges; n++) {
    if(mesh.edges[n].code == MESH_PERMEABLE_NODE) {
      mesh.edges[n].code = MESH_ELEVATION_NODE;
      }
    }
  treat_boundarycode(mesh, &gOpenBCs[level], z2D_discretisation, MESH_ELEVATION_NODE);
  printf("#open boundary, edges %d, nodes %d\n", gOpenBCs[level].nedges, gOpenBCs[level].nnodes);

  boundary_hbw[level]    = gOpenBCs[level].hbw;
  boundary_matrix[level] = new double[gOpenBCs[level].nnodes * (3 * boundary_hbw[level] + 1)];
  boundary_pivot[level]  = new int[gOpenBCs[level].nnodes];

/**----------------------------------------------------------------------
  Flux boundary conditions*/
/*----------------------------------------------------------------------
  set external land limits (1) and islands (2) with similar code*/
  for(n = 0; n < mesh.nedges; n++) {
    if(mesh.edges[n].code == MESH_ISLAND_NODE) {
      mesh.edges[n].code = MESH_LAND_NODE;
      }
    }

  treat_boundarycode(mesh, &gFluxBCs[level], u2D_discretisation, MESH_LAND_NODE);
  printf("#flux boundary, edges %d, nodes %d\n", gFluxBCs[level].nedges, gFluxBCs[level].nnodes);
  
/**----------------------------------------------------------------------
  Periodic boundary conditions*/
  /* STS: NO periodic */
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int boundary_normal_LGP0(specification_obsolete_t spec, mesh_t mesh, int k, double *nx, double *ny, double *L)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
  compute the normal vector, outward oriented, of the boundary at node n
  Its components are [nx,ny], length of segment is nn
----------------------------------------------------------------------*/
{
  int j,m,n,side[2],status;
  int count;

  m=spec.nodes[k];

  count=0;
  for(j=0;j<3;j++) {
    n=mesh.triangles[m].edges[j];
    if(mesh.edges[n].code==spec.code) {
      side[count]=j;
      count++;
      }
    }

  switch(count) {
    case 0:
      check_error(-1, "mal-formed solid boundary information", __LINE__, __FILE__, 1);
      break;
    case 1:
      j=side[0];
      *nx=mesh.triangles[m].nx[j];
      *ny=mesh.triangles[m].ny[j];
      *L =mesh.triangles[m].l [j];
      break;
    case 2:
      j=3-side[0]-side[1];
      *nx=mesh.triangles[m].nx[j];
      *ny=mesh.triangles[m].ny[j];
      *L =mesh.triangles[m].l [j];
      break;
    }

  status=0;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int boundary_normal_LGP1(specification_obsolete_t spec, mesh_t mesh, int n, double *nx, double *ny, double *nn)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
  compute the normal vector, outward oriented, of the boundary at node n
  Its components are [nx,ny], length of segment is nn
----------------------------------------------------------------------*/
{
  int status;
  double dX, dY, ds;
  double alpha, sinus, cosinus, size;
  double lon1, lon2, lat1, lat2;
  int k, n1, n2;

/**----------------------------------------------------------------------
  check right-hand side  neighbour*/
  k = spec.r_node[n];
  if(k == -1) {
    k = n;
    }
  n1 = spec.nodes[k];

/**----------------------------------------------------------------------
  check left-hand side  neighbour*/
  k = spec.l_node[n];
  if(k == -1) {
    k = n;
    }
  n2 = spec.nodes[k];

/**----------------------------------------------------------------------------
  coordinates in degrees, conversion needed */
  lon1 = mesh.vertices[n1].lon*d2r;
  lon2 = mesh.vertices[n2].lon*d2r;

  lat1 = mesh.vertices[n1].lat*d2r;
  lat2 = mesh.vertices[n2].lat*d2r;

  alpha = 0.5 * (lat1 + lat2);
/**----------------------------------------------------------------------------
  alpha is the mid latitude of node n1 and n2 */
  dX = R * cos(alpha) * (lon2 - lon1);
  dY = R * (lat2 - lat1);

  ds = sqrt(dX * dX + dY * dY);

  *nx = +dY / ds;
  *ny = -dX / ds;
  *nn = ds;
  status=0;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int boundary_normal_DGP1(specification_obsolete_t spec, mesh_t mesh, int n, double *nx, double *ny, double *L)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
  compute the normal vector, outward oriented, of the boundary at node n
  Its components are [nx,ny], length of segment is nn
----------------------------------------------------------------------*/
{
  int status;
  double dX, dY, ds;
  double alpha, sinus, cosinus, size;
  double lon1, lon2, lat1, lat2;
  int k, n1, n2;
  discretisation_t descriptor=mesh.DGP1descriptor;

/**----------------------------------------------------------------------
  check right-hand side  neighbour*/
  k = spec.r_node[n];
  if(k == -1) {
    k = n;
    }
  n1 = spec.nodes[k];

/**----------------------------------------------------------------------
  check left-hand side  neighbour*/
  k = spec.l_node[n];
  if(k == -1) {
    k = n;
    }
  n2 = spec.nodes[k];

/**----------------------------------------------------------------------------
  coordinates in degrees, conversion needed */
  lon1 = descriptor.nodes[n1].lon*d2r;
  lon2 = descriptor.nodes[n2].lon*d2r;

  lat1 = descriptor.nodes[n1].lat*d2r;
  lat2 = descriptor.nodes[n2].lat*d2r;

  alpha = 0.5 * (lat1 + lat2);
/**----------------------------------------------------------------------------
  alpha is the mid latitude of node n1 and n2 */
  dX = R * cos(alpha) * (lon2 - lon1);
  dY = R * (lat2 - lat1);

  ds = sqrt(dX * dX + dY * dY);

  *nx = +dY / ds;
  *ny = -dX / ds;
  *L  = ds;
  status=0;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int boundary_normal_NCP1(specification_obsolete_t spec, mesh_t mesh, int k, double *nx, double *ny, double *L)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
  compute the normal vector, outward oriented, of the boundary at node n
  Its components are [nx,ny], length of segment is L
----------------------------------------------------------------------*/
{
  int status;
  int n,node;
  discretisation_t descriptor=mesh.NCP1descriptor;

  node = spec.nodes[k];

/*----------------------------------------------------------------------
  capillo-tracté */ 
  n=mesh.NCP1descriptor.nodes[node].edges[0];

  *nx =  mesh.edges[n].Ty;
  *ny = -mesh.edges[n].Tx;
  *L  =  mesh.edges[n].L;
  status=0;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int boundary_normal_DNP1(specification_obsolete_t spec, mesh_t mesh, int k, double *nx, double *ny, double *L)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
  compute the normal vector, outward oriented, of the boundary at node n
  Its components are [nx,ny], length of segment is L
----------------------------------------------------------------------*/
{
  int status;
  int n,node;
  discretisation_t descriptor=mesh.DNP1descriptor;

  node = spec.nodes[k];

  n=mesh.DNP1descriptor.nodes[node].edges[0];

  *nx =  mesh.edges[n].Ty;
  *ny = -mesh.edges[n].Tx;
  *L  =  mesh.edges[n].L;
  status=0;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int boundary_normal_LGP2(specification_obsolete_t spec, mesh_t mesh, int n, double *nx, double *ny, double *nn)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
  compute the normal vector, outward oriented, of the boundary at node n
  Its components are [nx,ny], length of segment is nn
----------------------------------------------------------------------*/
{
  int status;
  double dX, dY, ds;
  double alpha, sinus, cosinus, size;
  double lon1, lon2, lat1, lat2;
  int k, n1, n2;
  node_t *P2nodes=mesh.LGP2descriptor.nodes;
  discretisation_t descriptor=mesh.LGP2descriptor;


/**----------------------------------------------------------------------
  check right-hand side  neighbour*/
  k = spec.r_node[n];
  if(k == -1) {
    k = n;
    }
  n1 = spec.nodes[k];

/**----------------------------------------------------------------------
  check left-hand side  neighbour*/
  k = spec.l_node[n];
  if(k == -1) {
    k = n;
    }
  n2 = spec.nodes[k];

/**----------------------------------------------------------------------------
  coordinates in degrees, conversion needed */
  lon1 = P2nodes[n1].lon*d2r;
  lon2 = P2nodes[n2].lon*d2r;

  lat1 = P2nodes[n1].lat*d2r;
  lat2 = P2nodes[n2].lat*d2r;

  alpha = 0.5 * (lat1 + lat2);
/**----------------------------------------------------------------------------
  alpha is the mid latitude of node n1 and n2 */
  dX = R * cos(alpha) * (lon2 - lon1);
  dY = R * (lat2 - lat1);

  ds = sqrt(dX * dX + dY * dY);

  *nx = +dY / ds;
  *ny = -dX / ds;
  *nn = ds;
  status=0;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int boundary_normal_CQP1(specification_obsolete_t spec, mesh_t mesh, int k, double *nx, double *ny, double *L)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
  compute the normal vector, outward oriented, of the boundary at node n
  Its components are [nx,ny], length of segment is nn
----------------------------------------------------------------------*/
{
  int status;
  int n,node;
  discretisation_t descriptor=mesh.CQP1descriptor;

  n = spec.nodes[k];

//  n=descriptor.nodes[node].edges[0];

  *nx =  mesh.edges[n].Ty;
  *ny = -mesh.edges[n].Tx;
  *L  =  mesh.edges[n].L;
  status=0;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int boundary_normal_CQN1(specification_obsolete_t spec, mesh_t mesh, int k, double *nx, double *ny, double *L)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
  compute the normal vector, outward oriented, of the boundary at node n
  Its components are [nx,ny], length of segment is nn
----------------------------------------------------------------------*/
{
  int status;
  int n,node;
  discretisation_t descriptor=mesh.CQN1descriptor;

  n = spec.nodes[k];

//  n=descriptor.nodes[node].edges[0];

  *nx =  mesh.edges[n].Ty;
  *ny = -mesh.edges[n].Tx;
  *L  =  mesh.edges[n].L;
  status=0;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int boundary_normal(specification_obsolete_t & spec, mesh_t & mesh, int node, double *nx, double *ny, double *L)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
  compute the normal vector, outward oriented, of the boundary at node n
  Its components are [nx,ny], length of segment is nn
----------------------------------------------------------------------*/
{
  int status;
  int discretisation=spec.discretisation;

  switch (discretisation) {
    case LGP0:
/*----------------------------------------------------------------------
      LGP0 discretisation */
      status=boundary_normal_LGP0(spec, mesh, node, nx, ny, L);
      break;

    case LGP1:
/*----------------------------------------------------------------------
      LGP1 discretisation */
      status=boundary_normal_LGP1(spec, mesh, node, nx, ny, L);
      break;

    case DGP1:
/*----------------------------------------------------------------------
      LGP1 discretisation */
      status=boundary_normal_DGP1(spec, mesh, node, nx, ny, L);
      break;

    case NCP1:
/*----------------------------------------------------------------------
      NCP1 discretisation */
      status=boundary_normal_NCP1(spec, mesh, node, nx, ny, L);
      break;

    case DNP1:
/*----------------------------------------------------------------------
      LGP1 discretisation */
      status=boundary_normal_DNP1(spec, mesh, node, nx, ny, L);
      break;

    case LGP2:
/*----------------------------------------------------------------------
      LGP2 discretisation */
      status=boundary_normal_LGP1(spec, mesh, node, nx, ny, L);  /// HERE !!!
      break;

    case CQP1:
/*----------------------------------------------------------------------
      CQP1 discretisation */
      status=boundary_normal_CQP1(spec, mesh, node, nx, ny, L);
      break;

    case CQN1:
/*----------------------------------------------------------------------
      CQP1 discretisation */
      status=boundary_normal_CQN1(spec, mesh, node, nx, ny, L);
      break;

    default:
      check_error(-1, "illegal discretisation", __LINE__, __FILE__, 1);
      break;
    }
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int boundary_normal(specification_obsolete_t & spec, mesh_t & mesh, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  if(spec.fluxes==0) {
    
    }
  status=boundary_normal(spec, mesh, n, &(spec.data[n].nx), &(spec.data[n].ny), &(spec.data[n].L));
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void build_boundary_conditions(mesh_t mesh, state2d_t *state2D, parameter_t data2D)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, n, n1, n2, n3, l;
  int Nboundary,status;
  double nx, ny, nn;

/**----------------------------------------------------------------------------
  Find boundary nodes and neighbours*/
  Nboundary = 0;
  for(l = 0; l < mesh.nlimits; l++) {
    Nboundary += mesh.limits[l].nvertex;
    }

  gFluxBCs[0].fluxes=new fluxspec_t[gFluxBCs[0].nnodes];
  gFluxBCs[0].data  =new dataspec_t[gFluxBCs[0].nnodes];
  for(i=0;i<gFluxBCs[0].nnodes;i++) {
    boundary_normal(gFluxBCs[0], mesh, i);
    }

  status=VelocityBCs_init(gFluxBCs[0]);
  
  printf("#Number of boundary vertices: %d \n", Nboundary);
  fprintf(echo, "#Number of boundary vertices: %d \n", Nboundary);

/*----------------------------------------------------------------------
  Write in buffer the boundary conditions ELEVATION*/
  if(gOpenBCs[0].nnodes != 0) {
    gOpenBCs[0].data  =new dataspec_t[gOpenBCs[0].nnodes];
    for(i = 0; i < gOpenBCs[0].nnodes; i++) {
      n=gOpenBCs[0].nodes[i];
      gOpenBCs[0].data[i].node = n;
      gOpenBCs[0].data[i].ancestor = i;
      switch (z2D_discretisation) {
        case LGP0:
          status=boundary_normal_LGP0(gOpenBCs[0], mesh, i, &nx, &ny, &nn);
          break;
        case LGP1:
          boundary_normal_LGP1(gOpenBCs[0], mesh, i, &nx, &ny, &nn);
          break;
        case LGP2:
          boundary_normal_LGP2(gOpenBCs[0], mesh, i, &nx, &ny, &nn);
          break;
        case CQP0:
          break;
        default:
          check_error(-1, "unknown discretisation", __LINE__, __FILE__, 1);
          break;
        }
      gOpenBCs[0].data[i].nx = nx;
      gOpenBCs[0].data[i].ny = ny;
      gOpenBCs[0].data[i].L  = nn;
      }

    if(tidal_bc) {
      gOpenBCs[0].tides=new tidespec_t[gOpenBCs[0].nnodes];
      for (n=0;n<gOpenBCs[0].nnodes;n++) {
	gOpenBCs[0].tides[n].node=gOpenBCs[0].nodes[n];
        gOpenBCs[0].tides[n].spectrum=&(WaveList);
        gOpenBCs[0].tides[n].ztide.allocate(WaveList.n);
        gOpenBCs[0].tides[n].utide.allocate(WaveList.n);
        gOpenBCs[0].tides[n].vtide.allocate(WaveList.n);
        }
      discretisation_t z_descriptor=get_descriptor(gFEmesh[0], z2D_discretisation);
      status=TidalBoundaryConditions(gFEmesh[0], z_descriptor, gOpenBCs[0]);
      }

/**----------------------------------------------------------------------
    initialization routine, flushing ok*/
    fflush(echo);

/*----------------------------------------------------------------------
    initialize z-code to -1 for open boundary conditions */
    for(i=0;i<gOpenBCs[0].nnodes;i++) {
      n=gOpenBCs[0].nodes[i];
      gdata2D[0].zCode[n] = FE_EXTERNAL_BOUNDARY_NODE;
      }
    }
}
