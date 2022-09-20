

/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

***************************************************************************/

#include "tools-structures.h"
#include "fe.h"
#include "matrix.h"
#include "zapper.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_LGP2xLGP2_0(mesh_t &mesh, discretisation_t &descriptor, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  LGP2 x LGP2, based on first level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int status;
  int i, m, n, n1, n2, n3, j, k, l, row, nrows, col, ncols;
  double A, hC, hc, temp1, temp2;
  double h;
  triangle_t *elt;
  int nndes, nelts;
  double *M;
  int *tmp;
  size_t count;

  nndes = descriptor.nnodes;
  nelts = mesh.ntriangles;

  nrows = mesh.nnghm + 1;
  ncols = 3*mesh.nnghm + 2;

  elt = mesh.triangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal = new int[nndes];
  ordering->pointer  = new int[nndes + 1];
  ordering->nrows    = nndes;
  tmp       = new int[nndes * ncols];

  for(n = 0; n < nndes * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  allocate arrays*/
  for(n = 0; n < nndes; n++) {
    tmp[n * ncols] = n;
    for(k = 0; k < descriptor.nodes[n].nnghbs; k++)
      tmp[n * ncols + k + 1] = descriptor.nodes[n].nghbs[k];
      }

  count = 0;
  for(n = 0; n < nndes; n++) {
    (ordering->pointer)[n]  = count;      /*pointer= column start index */
    (ordering->cardinal)[n] = descriptor.nodes[n].nnghbs + 1;
    count += (ordering->cardinal)[n];
    if((ordering->cardinal)[n]>=ncols) {
      check_error(-1, "matrix index overflow", __LINE__, __FILE__, 1);
      }
    }
  (ordering->pointer)[n] = count;

  ordering->incidence = new int[count];

  for(n = 0; n < nndes; n++) {
    for(k = 0; k < (ordering->cardinal)[n]; k++) {
      (ordering->incidence)[(ordering->pointer)[n] + k] = tmp[n * ncols + k];
      if((ordering->incidence)[(ordering->pointer)[n] + k] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }
    
  ordering->nrows     = nndes;

  status= matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,nndes);
  status=matrix_check_symmetry(ordering->incidence,ordering->pointer,ordering->cardinal,nndes);

  zaparr(tmp);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_LGP2xLGP2_1(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  LGP2 x LGP2 based on first and second level neighbouring relationships

  To be used with LGP1 velocities

  ----------------------------------------------------------------------*/
{
  int i, m, n, n1, n2, n3, j, k, l, status;
  int opposite[2][3], row, nrows, col, ncols,nndes;
  double A, hC, hc, temp1, temp2;
  double h;
  triangle_t *elt;
  double *M;
  int add,*tmp;
  int *incidence, *cardinal, *pointer;

  size_t count;

  discretisation_t descriptor=mesh.LGP2descriptor;
  node_t *P2nodes=mesh.LGP2descriptor.nodes;

  nndes = descriptor.nnodes;
  nrows = -1;
  ncols = 9*mesh.nnghm*mesh.nnghm + 1;

  elt = mesh.triangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  cardinal = new int[nndes];
  pointer  = new int[nndes + 1];
  tmp      = new int[nndes * ncols];

/*----------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < nndes * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < nndes; n++) {
    cardinal[n] = P2nodes[n].nnghbs  + 1;
    if(cardinal[n]>=ncols) {
      check_error(-1, "matrix index arrays overflow", __LINE__, __FILE__, 1);
      }
    tmp[n * ncols] = n;
/* *----------------------------------------------------------------------
    1st level of neighbours*/
    for(k = 0; k < P2nodes[n].nnghbs; k++) {
      tmp[n * ncols + k + 1] = P2nodes[n].nghbs[k];
      }
/* *----------------------------------------------------------------------
    2nd level of neighbours*/
    for(k = 0; k < P2nodes[n].nnghbs; k++) {
      m=P2nodes[n].nghbs[k];
      for(l = 0; l < P2nodes[m].nnghbs; l++) {
        if(array_index(&(tmp[n * ncols]),cardinal[n],P2nodes[m].nghbs[l]) == -1) {
          tmp[n * ncols + cardinal[n]] = P2nodes[m].nghbs[l];
          cardinal[n]++;
          }
        }
      }
    if(cardinal[n]>ncols) {
      check_error(-1, "matrix index arrays overflow", __LINE__, __FILE__, 1);
      }
    }

  count = 0;
  for(n = 0; n < nndes; n++) {
    pointer[n] = count;      /*pointer= column start index */
    count += cardinal[n];
    }

  pointer[nndes] = count;
  incidence = new int[count];
  for(n = 0; n < count; n++) {
    incidence[n]=-1;
    }

  for(n = 0; n < nndes; n++) {
    for(k = 0; k < cardinal[n]; k++) {
      incidence[pointer[n] + k] = tmp[n * ncols + k];
      }
    }
  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows     = nndes;

  status= matrix_reorder(ordering);

  status=matrix_check_consistency(*ordering);
  status=matrix_check_symmetry(*ordering);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_LGP2xLGP2_2(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  LGP2 x LGP2 based on first and second level neighbouring relationships

  To be used with NCP1 velocities

  ----------------------------------------------------------------------*/
{
  int status;
  int i, m, n, n1, n2, n3, j, k, l;
  int opposite[2][3], row, nrows, col, ncols,nndes;
  double A, hC, hc, temp1, temp2;
  double h;
  triangle_t *elt;
  double *M;
  int add,*tmp;
  int *incidence, *cardinal, *pointer;

  size_t count;

  discretisation_t descriptor=mesh.LGP2descriptor;
  node_t *P2nodes=mesh.LGP2descriptor.nodes;

  nndes = descriptor.nnodes;
  nrows = 2*mesh.nnghm + 1;
  ncols = 6*mesh.nnghm + 1;

  elt = mesh.triangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  cardinal = new int[nndes];
  pointer  = new int[nndes + 1];
  tmp      = new int[nndes * ncols];

/*----------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < nndes * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < nndes; n++) {
//    pointer[n] = count;      /*pointer= column start index */
    cardinal[n] = P2nodes[n].nnghbs  + 1;
    if(cardinal[n]>=ncols) {
      check_error(-1, "matrix index arrays overflow", __LINE__, __FILE__, 1);
      }
//    count += cardinal[n];
    tmp[n * ncols] = n;
    for(k = 0; k < P2nodes[n].nnghbs; k++) {
      tmp[n * ncols + k + 1] = P2nodes[n].nghbs[k];
      }
    }

/* *----------------------------------------------------------------------
  2nd level of neighbours*/
  for(n = 0; n < mesh.nedges; n++) {
    if( mesh.edges[n].nshared==1) continue;
    n1=mesh.edges[n].extremity[0];
    n2=mesh.edges[n].extremity[1];
    for(k = 0; k < mesh.edges[n].nshared; k++) {
      m=mesh.edges[n].shared[k];
      j=mesh.edges[n].eindex[k];
      switch (j) {
        case 0:
          opposite[k][0]=descriptor.NIbE[m][5];
          opposite[k][1]=descriptor.NIbE[m][0];
          opposite[k][2]=descriptor.NIbE[m][1];
          break;
        case 1:
          opposite[k][0]=descriptor.NIbE[m][1];
          opposite[k][1]=descriptor.NIbE[m][2];
          opposite[k][2]=descriptor.NIbE[m][3];
          break;
        case 2:
          opposite[k][0]=descriptor.NIbE[m][3];
          opposite[k][1]=descriptor.NIbE[m][4];
          opposite[k][2]=descriptor.NIbE[m][5];
          break;
        }
      }
    for(j = 0; j < 3; j++) {
      n1=opposite[0][j];
      for(i = 0; i < 3; i++) {
        n2=opposite[1][i];
/* *----------------------------------------------------------------------
        inner nodes with 4 neigbours duplicate connections, avoid it*/
        add=1;
        for(k = 0; k < cardinal[n1]; k++) {
          m=tmp[n1 * ncols + k];
          if(m==n2) {
            add=0;
            }
          }
        if(add==0) continue;
        k=cardinal[n1];
        if(k > ncols-1) {
          check_error(-1, "matrix index arrays overflow", __LINE__, __FILE__, 1);
          }
        tmp[n1 * ncols + k] = n2;
        cardinal[n1]++;
        count++;
        k=cardinal[n2];
        if(k > ncols-1) {
          check_error(-1, "matrix index arrays overflow", __LINE__, __FILE__, 1);
          }
        tmp[n2 * ncols + k] = n1;
        cardinal[n2]++;
        count++;
        }
      }
    }

  count = 0;
  for(n = 0; n < nndes; n++) {
    pointer[n] = count;      /*pointer= column start index */
    count += cardinal[n];
    }

  pointer[nndes] = count;
  incidence = new int[count];
  for(n = 0; n < count; n++) {
    incidence[n]=-1;
    }

  for(n = 0; n < nndes; n++) {
    for(k = 0; k < cardinal[n]; k++) {
      incidence[pointer[n] + k] = tmp[n * ncols + k];
      }
    }
  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows     = nndes;

  status=matrix_reorder(ordering);
  status=matrix_check_consistency(incidence,pointer,cardinal,nndes);
  status=matrix_check_symmetry(incidence,pointer,cardinal,nndes);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_LGP2xNCP1(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  P2 x ncP1 based on first level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int status;
  int i, m, n, n1, n2, n3, j, k, l;
  int opposite[2], row, nrows, col, ncols;
  double h;
  triangle_t *elt;
  int *tmp;
  int *incidence, *cardinal, *pointer;

  size_t count;

  discretisation_t descriptor=mesh.LGP2descriptor;
  node_t *P2nodes=mesh.LGP2descriptor.nodes;

  nrows = 2*mesh.nnghm;
  ncols = 2*mesh.nnghm;

  elt = mesh.triangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  cardinal = new int[descriptor.nnodes];
  pointer  = new int[descriptor.nnodes + 1];
  tmp      = new int[descriptor.nnodes * ncols];

/*----------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < descriptor.nnodes * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < descriptor.nnodes; n++) {
    pointer[n] = count;      /*pointer= column start index */
    cardinal[n] = P2nodes[n].nedges;
    if(cardinal[n]>=ncols) {
      check_error(-1, "matrix index arrays overflow", __LINE__, __FILE__, 1);
      }
    count += cardinal[n];
    for(k = 0; k <  P2nodes[n].nedges; k++) {
      tmp[n*ncols + k] = P2nodes[n].edges[k];
      }
    }
  pointer[n] = count;

  incidence = new int[count];

  for(n = 0; n < descriptor.nnodes; n++) {
    for(k = 0; k < cardinal[n]; k++) {
      incidence[pointer[n] + k] = tmp[n * ncols + k];
      }
    }
  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows     = descriptor.nnodes;

  status=matrix_reorder(ordering);
  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,descriptor.nnodes);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_LGP2xLGP1(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  LGP2 x LGP1 based on first level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int status;
  int i, m, n, n1, n2, n3, j, k, l;
  int opposite[2], row, nrows, col, ncols;
  double h;
  triangle_t *elt;
  int *tmp;
  int *incidence, *cardinal, *pointer;

  size_t count;

  discretisation_t descriptor=mesh.LGP2descriptor;
  node_t *P2nodes=mesh.LGP2descriptor.nodes;

  nrows = -1;
  ncols = 2*mesh.nnghm;

  elt = mesh.triangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  cardinal = new int[descriptor.nnodes];
  pointer  = new int[descriptor.nnodes + 1];
  tmp      = new int[descriptor.nnodes * ncols];

/*----------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < descriptor.nnodes * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < descriptor.nnodes; n++) {
    pointer[n] = count;      /*pointer= column start index */
    cardinal[n] = P2nodes[n].nvtces;
    if(cardinal[n]>=ncols) {
      check_error(-1, "matrix index arrays overflow", __LINE__, __FILE__, 1);
      }
    count += cardinal[n];
    for(k = 0; k <  P2nodes[n].nvtces; k++) {
      tmp[n*ncols + k] = P2nodes[n].vtces[k];
      }
    }
  pointer[n] = count;

  incidence = new int[count];

  for(n = 0; n < descriptor.nnodes; n++) {
    for(k = 0; k < cardinal[n]; k++) {
      incidence[pointer[n] + k] = tmp[n * ncols + k];
      }
    }
  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows     = descriptor.nnodes;

  status=matrix_reorder(ordering);
  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,descriptor.nnodes);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_LGP1xLGP2(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  LGP1 (rows) x LGP2 (columns) based on first level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int status;
  int i, m, n, nn, j, k, l;
  int opposite[2], *connexions,nconnexions, row, nrows, col, ncols;
  double A, hC, hc, temp1, temp2;
  double h;
  triangle_t *elt;
  double *M;
  int *tmp;
  int *incidence, *cardinal, *pointer;
  int duplicated;
  size_t count;
  discretisation_t descriptor=mesh.LGP2descriptor;

  nrows = -1;
  ncols =  3*(mesh.nnghm+1);

  elt = mesh.triangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  cardinal   = new int[mesh.nvtxs];
  pointer    = new int[mesh.nvtxs + 1];
  tmp        = new int[mesh.nvtxs * ncols];
  connexions = new int[ncols];

/*----------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < mesh.nvtxs * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < mesh.nvtxs; n++) {
    pointer[n]  = count;      /*pointer= column start index */
    nconnexions=0;
    for(k = 0; k < mesh.vertices[n].nelmts; k++) {
      m=mesh.vertices[n].elmts[k];
      for (j=0;j<descriptor.nnpe;j++) {
        nn = descriptor.NIbE[m][j];
        duplicated=pos(nn,connexions,nconnexions);
        if(duplicated==-1) {
          connexions[nconnexions] = nn;
          nconnexions++;
          if(nconnexions>=ncols) {
            check_error(-1, "matrix index arrays overflow", __LINE__, __FILE__, 1);
            }
          }
        }
      }
next:
    cardinal[n] = nconnexions;
    count += cardinal[n];
    for(k = 0; k < nconnexions; k++) {
      tmp[n*ncols + k] = connexions[k];
      }
    }
  pointer[n] = count;

  incidence = new int[count];

  for(n = 0; n < mesh.nvtxs; n++) {
    for(k = 0; k < cardinal[n]; k++) {
      incidence[pointer[n] + k] = tmp[n * ncols + k];
      }
    }

  zaparr(tmp);
  delete[] connexions;

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows     = mesh.nvtxs;

/* *----------------------------------------------------------------------
  re-ordering of the column index*/
  status=matrix_reorder(ordering);
/* *----------------------------------------------------------------------
  check matrix duplicated connections*/
  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,mesh.nvtxs);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_DGP1xLGP2(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  DGP1 (rows) x LGP2 (columns) based on first level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int status;
  int i, m, n, nn, j, k, l;
  int opposite[2], *connexions,nconnexions, row, nrows, col, ncols;
  double A, hC, hc, temp1, temp2;
  double h;
  triangle_t *elt;
  double *M;
  int *tmp;
  int *incidence, *cardinal, *pointer;
  int duplicated;
  size_t count;

  discretisation_t descriptor=mesh.LGP2descriptor;

  nrows = -1;
  ncols =  6;

  elt = mesh.triangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  cardinal   = new int[mesh.DGP1descriptor.nnodes];
  pointer    = new int[mesh.DGP1descriptor.nnodes + 1];
  tmp        = new int[mesh.DGP1descriptor.nnodes * ncols];

/*----------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < mesh.DGP1descriptor.nnodes * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k = 0; k < mesh.DGP1descriptor.nnpe; k++) {
      n=mesh.DGP1descriptor.NIbE[m][k];
      pointer[n]  = count;      /*pointer= column start index */
      cardinal[n] = 6;
      for (j=0;j<descriptor.nnpe;j++) {
        tmp[n*ncols + j] = descriptor.NIbE[m][j];
        }
      count += cardinal[n];
      }
    }
  pointer[mesh.DGP1descriptor.nnodes] = count;

  incidence = new int[count];

  for(n = 0; n < mesh.DGP1descriptor.nnodes; n++) {
    for(k = 0; k < cardinal[n]; k++) {
      incidence[pointer[n] + k] = tmp[n * ncols + k];
      }
    }

  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows     = mesh.DGP1descriptor.nnodes;

/* *----------------------------------------------------------------------
  re-ordering of the column index*/
  status=matrix_reorder(ordering);
/* *----------------------------------------------------------------------
  check matrix duplicated connections*/
  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,mesh.nvtxs);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_NCP1xLGP2(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  ncP1 (rows) x P1 (columns) based on first level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int status;
  int i, m, n, n1, n2, n3, j, k, l;
  int opposite[2], connexions[9],nconnexions, row, nrows, col, ncols;
  double A, hC, hc, temp1, temp2;
  double h;
  triangle_t *elt;
  double *M;
  int *tmp;
  int *incidence, *cardinal, *pointer;

  size_t count;

  discretisation_t descriptor=mesh.LGP2descriptor;
//  node_t *P2nodes=gP2nodes;

  nrows = 4;
  ncols = 9;

  elt = mesh.triangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  cardinal = new int[mesh.nedges];
  pointer  = new int[mesh.nedges + 1];
  tmp      = new int[mesh.nedges * ncols];

/*----------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < mesh.nedges * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < mesh.nedges; n++) {
    pointer[n]  = count;      /*pointer= column start index */

    m=mesh.edges[n].shared[0];
    for (j=0;j<descriptor.nnpe;j++) {
      connexions[j] = descriptor.NIbE[m][j];
      }
    nconnexions=descriptor.nnpe;
    if (mesh.edges[n].nshared==1) goto next;
    m=mesh.edges[n].shared[1];
    j=mesh.edges[n].eindex[1];
    switch (j) {
      case 0:
        connexions[6]=descriptor.NIbE[m][5];
        connexions[7]=descriptor.NIbE[m][0];
        connexions[8]=descriptor.NIbE[m][1];
        break;
      case 1:
        connexions[6]=descriptor.NIbE[m][1];
        connexions[7]=descriptor.NIbE[m][2];
        connexions[8]=descriptor.NIbE[m][3];
        break;
      case 2:
        connexions[6]=descriptor.NIbE[m][3];
        connexions[7]=descriptor.NIbE[m][4];
        connexions[8]=descriptor.NIbE[m][5];
        break;
      }
    nconnexions+=3;
next:
    cardinal[n] = nconnexions;
    count += cardinal[n];
    for(k = 0; k < nconnexions; k++) {
      tmp[n*ncols + k] = connexions[k];
      }
    }
  pointer[n] = count;

  incidence = new int[count];

  for(n = 0; n < mesh.nedges; n++) {
    for(k = 0; k < cardinal[n]; k++) {
      incidence[pointer[n] + k] = tmp[n * ncols + k];
      }
    }

  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows     = mesh.nedges;

/* *----------------------------------------------------------------------
  re-ordering of the column index*/
  status=matrix_reorder(ordering);
/* *----------------------------------------------------------------------
  check matrix duplicated connections*/
  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,mesh.nedges);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_LGP2xDGP1(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  LGP2 x DGP1 based on first level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int status;
  int i, m, n, n1, n2, n3, j, k, l;
  int opposite[2], row, nrows, col, ncols;
  double h;
  triangle_t *elt;
  int *tmp;
  int *incidence, *cardinal, *pointer;

  size_t count;

  discretisation_t descriptor=mesh.LGP2descriptor;
  node_t *P2nodes=mesh.LGP2descriptor.nodes;

  nrows = -1;
  ncols = 3*(mesh.nnghm+1);

  elt = mesh.triangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  cardinal = new int[descriptor.nnodes];
  pointer  = new int[descriptor.nnodes + 1];
  tmp      = new int[descriptor.nnodes * ncols];

/*----------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < descriptor.nnodes * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < descriptor.nnodes; n++) {
    pointer[n]  = count;      /*pointer= column start index */
    cardinal[n] = P2nodes[n].nelmts*3;
    if(cardinal[n]>=ncols) {
      check_error(-1, "matrix index arrays overflow", __LINE__, __FILE__, 1);
      }
    count += cardinal[n];
    for(k = 0; k <  P2nodes[n].nelmts; k++) {
      for(l = 0; l <  3; l++) {
        tmp[n*ncols + 3*k+l] = P2nodes[n].elmts[k]*3+l;
        }
      }
    }
  pointer[n] = count;

  incidence = new int[count];

  for(n = 0; n < descriptor.nnodes; n++) {
    for(k = 0; k < cardinal[n]; k++) {
      incidence[pointer[n] + k] = tmp[n * ncols + k];
      }
    }
  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows     = descriptor.nnodes;

  status= matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,descriptor.nnodes);

  return (0);
}

