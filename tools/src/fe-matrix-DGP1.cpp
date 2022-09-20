
/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2012

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

  int init_packed_DGP1xDGP1_0(mesh_t &mesh, discretisation_t &descriptor, ordering_t *ordering)
//  int init_packed_LGP2xLGP2_0(mesh_t mesh,  ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  DGP1 x DGP1, based on first level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int i, m, n, j, k, l, row, nrows, col, ncols, status;
  int nndes, nelts;
  double *M;
  int *tmp;
  size_t count;

  nndes = descriptor.nnodes;
  nelts = mesh.ntriangles;

  nrows = -1;
  ncols =  3;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal = new int[nndes];
  ordering->pointer  = new int[nndes + 1];
  ordering->nrows    = nndes;
  tmp                = new int[nndes * ncols];

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
    if((ordering->cardinal)[n]>ncols) {
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

  ordering->nrows = descriptor.nnodes;

  status= matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,nndes);
  status=matrix_check_symmetry(ordering->incidence,ordering->pointer,ordering->cardinal,nndes);

  zaparr(tmp);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_DNP1xDNP1_0(mesh_t &mesh, discretisation_t &descriptor, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometryic matrix symetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  DGP1 x DGP1, based on first level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int i, m, n, j, k, l, row, nrows, col, ncols, status;
  int nndes, nelts;
  double *M;
  int *tmp;
  size_t count;

  nndes = descriptor.nnodes;
  nelts = mesh.ntriangles;

  nrows = -1;
  ncols =  3;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal = new int[nndes];
  ordering->pointer  = new int[nndes + 1];
  ordering->nrows    = nndes;
  tmp                = new int[nndes * ncols];

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
    if((ordering->cardinal)[n]>ncols) {
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

  ordering->nrows = descriptor.nnodes;

  status= matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,nndes);
//  status=matrix_check_symmetry(ordering->incidence,ordering->pointer,ordering->cardinal,nndes);

  zaparr(tmp);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_DNP1xDNP1_1(mesh_t &mesh, discretisation_t &descriptor, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometryic matrix symetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  DNP1 x DNP1, based on cross-edge relationships

  ----------------------------------------------------------------------*/
{
  int i, m, n, n1, n2, j, k, l, row, nrows, col, ncols, status;
  int nndes, nelts;
  double *M;
  int *tmp;
  size_t count;

  nndes = descriptor.nnodes;
  nelts = mesh.ntriangles;

  nrows = -1;
  ncols =  2;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal = new int[nndes];
  ordering->pointer  = new int[nndes + 1];
  ordering->nrows    = nndes;
  tmp                = new int[nndes * ncols];

  for(n = 0; n < nndes * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  allocate arrays*/
  count = 0;
  for(n = 0; n < nndes; n++) {
    tmp[n * ncols] = n;
    (ordering->cardinal)[n] = 1;
    }

  for(n = 0; n < mesh.nedges; n++) {
    if(mesh.edges[n].nshared==1) continue;
    i=pos(n,mesh.triangles[mesh.edges[n].shared[0]].edges,3);
    j=pos(n,mesh.triangles[mesh.edges[n].shared[1]].edges,3);
    n1=descriptor.NIbE[mesh.edges[n].shared[0]][i];
    n2=descriptor.NIbE[mesh.edges[n].shared[1]][j];
    if ( (ordering->cardinal)[n1] == 1) {
      tmp[n1 * ncols + 1 ] = n2;
      (ordering->cardinal)[n1] ++;
      }
    if ( (ordering->cardinal)[n2] == 1) {
      tmp[n2 * ncols + 1 ] = n1;
      (ordering->cardinal)[n2] ++;
      }
    }

  count = 0;
  for(n = 0; n < nndes; n++) {
    (ordering->pointer)[n]  = count;
    count += (ordering->cardinal)[n];
    if((ordering->cardinal)[n]>ncols) {
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

  ordering->nrows = descriptor.nnodes;

  status= matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,nndes);
//  status=matrix_check_symmetry(ordering->incidence,ordering->pointer,ordering->cardinal,nndes);

  zaparr(tmp);

  return (0);
}

