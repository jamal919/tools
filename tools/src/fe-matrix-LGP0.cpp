
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

  int init_packed_LGP0xLGP0_0(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------

  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  LGP0 x LGP0, diagonal

  ----------------------------------------------------------------------*/
{
  int m;
  int nelts;

  nelts = mesh.ntriangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal  = new int[nelts];
  ordering->pointer   = new int[nelts + 1];
  ordering->incidence = new int[nelts];
  ordering->nrows=nelts;

/*----------------------------------------------------------------------
  */

  for(m = 0; m < nelts; m++) {
    (ordering->cardinal)[m]  = 1;
    (ordering->pointer)[m]   = m;
    (ordering->incidence)[m] = m;
    }
  (ordering->pointer)[m] = mesh.ntriangles;

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_LGP0xLGP0_1(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------

  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  lGP0 x LGP0, 1 level of neighbouring

  ----------------------------------------------------------------------*/
{
  int j,k,l,m,m1,m2,n;
  int ncols,*tmp,count;

  ncols = 4;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal = new int[mesh.ntriangles];
  ordering->pointer  = new int[mesh.ntriangles + 1];
  ordering->nrows=mesh.ntriangles;
  tmp                = new int[mesh.ntriangles * ncols];

/*----------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < mesh.ntriangles * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(m = 0; m < mesh.ntriangles; m++) {
    ordering->pointer[m] = count;      /*pointer= column start index */
    ordering->cardinal[m] = mesh.triangles[m].nngh;
    count += ordering->cardinal[m];
    for(k = 0; k < mesh.triangles[m].nngh; k++) {
      tmp[m*ncols + k] = mesh.triangles[m].ngh[k];
      }
    }
  ordering->pointer[n] = count;

  ordering->incidence = new int[count];

/*----------------------------------------------------------------------
  re-ordering of the column index*/
  for(m = 0; m < mesh.ntriangles; m++) {
  iterate:
    for(k = 0; k < ordering->cardinal[m]; k++) {
      n = tmp[m * ncols + k];
      for(l = ordering->cardinal[m] - 1; l > k; l--) {
        if(tmp[m * ncols + l] < n) {
          tmp[m * ncols + k] = tmp[m * ncols + l];
          tmp[m * ncols + l] = n;
          if(n == -1) {
            check_error(-1, "matrix index arrays overflow", __LINE__, __FILE__, 1);
            }
          goto iterate;
          }
        }
      }
    }

  for(n = 0; n < mesh.ntriangles; n++) {
    for(k = 0; k < ordering->cardinal[n]; k++) {
       /*incidence= row index */
      ordering->incidence[ordering->pointer[n] + k] = tmp[n * ncols + k];
      if(ordering->incidence[ordering->pointer[n] + k] == -1) {
        check_error(-1, "matrix index arrays overflow", __LINE__, __FILE__, 1);
        }
      }
    }

  zaparr(tmp);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_LGP0xLGP0_LGP1(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------

  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  LGP0 x LGP0, 2 level of neighbouring

  To be used with LGP1 velocities

  ----------------------------------------------------------------------*/
{
  int j,k,l,m,m1,m2,mm,mmm,n,pos,status;
  int ncols,*tmp,count,nngh;

  status= fe_element_surrounding(mesh);

  ncols = mesh.nnghm*mesh.nnghm+1;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal = new int[mesh.ntriangles];
  ordering->pointer  = new int[mesh.ntriangles + 1];
  ordering->nrows=mesh.ntriangles;
  tmp                = new int[mesh.ntriangles * ncols];

/*----------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < mesh.ntriangles * ncols; n++)
    tmp[n] = -1;

  count = 0;
  for(m = 0; m < mesh.ntriangles; m++) {
    ordering->pointer[m]  = count;      /*pointer= column start index */
//    ordering->cardinal[m] = mesh.triangles[m].nsurrounding+1;
    ordering->cardinal[m] = mesh.triangles[m].nsurrounding;
//    tmp[m*ncols] = m;
/*----------------------------------------------------------------------
    first level of neighbours*/
    for(k = 0; k < mesh.triangles[m].nsurrounding; k++) {
//      tmp[m*ncols + k +1] = mesh.triangles[m].surrounding[k];
      tmp[m*ncols + k] = mesh.triangles[m].surrounding[k];
      }
/*----------------------------------------------------------------------
    second level of neighbours*/
    for(k = 0; k < mesh.triangles[m].nsurrounding; k++) {
      mm = mesh.triangles[m].surrounding[k];
      for(l = 0; l < mesh.triangles[mm].nsurrounding; l++) {
        mmm=mesh.triangles[mm].surrounding[l];
        pos=array_index(&(tmp[m*ncols]), ordering->cardinal[m], mmm);
        if(pos==-1) {
          tmp[m*ncols + ordering->cardinal[m]] = mmm;
          ordering->cardinal[m]++;
          }
        }
      }
    count += ordering->cardinal[m];
    if(ordering->cardinal[m]>ncols) {
      check_error(-1, "matrix index overflow", __LINE__, __FILE__, 1);
      }
    }

  ordering->pointer[mesh.ntriangles] = count;

  ordering->incidence = new int[count];

/*----------------------------------------------------------------------
  re-ordering of the column index*/
  for(m = 0; m < mesh.ntriangles; m++) {
  iterate:
    for(k = 0; k < ordering->cardinal[m]; k++) {
      n = tmp[m * ncols + k];
      for(l = ordering->cardinal[m] - 1; l > k; l--) {
        if(tmp[m * ncols + l] < n) {
          tmp[m * ncols + k] = tmp[m * ncols + l];
          tmp[m * ncols + l] = n;
          if(n == -1) {
            check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
            }
          goto iterate;
          }
        }
      }
    }

  for(n = 0; n < mesh.ntriangles; n++) {
    for(k = 0; k < ordering->cardinal[n]; k++) {
/* *----------------------------------------------------------------------------
      incidence= row index */
      ordering->incidence[ordering->pointer[n] + k] = tmp[n * ncols + k];
      if(ordering->incidence[ordering->pointer[n] + k] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }

/* *-----------------------------------------------------------------------
  some testing*/
//   for(n = 0; n <count; n++) {
//     if(ordering->incidence[n] == -1) {
//       check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
//       }
//     }
/* *-----------------------------------------------------------------------
  symmetry testing*/
//   for(n = 0; n < mesh.ntriangles; n++) {
//     for(k = 0; k < ordering->cardinal[n]; k++) {
// /* *----------------------------------------------------------------------------
//       incidence= row index */
//       m=ordering->incidence[ordering->pointer[n] + k];
//       pos=array_index(&(tmp[m*ncols]), ordering->cardinal[m], n);
//       if(pos== -1) {
//         check_error(-1, "matrix symmetry anomaly", __LINE__, __FILE__, 1);
//         }
//       }
//     }
  zaparr(tmp);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,mesh.ntriangles);
  status=matrix_check_symmetry(ordering->incidence,ordering->pointer,ordering->cardinal,mesh.ntriangles);


  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_LGP0xLGP0_2(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------

  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  LGP0 x LGP0, 2 level of neighbouring

  To be used with NCP1 velocities

  ----------------------------------------------------------------------*/
{
  int j,k,l,m,m1,m2,mm,mmm,n,pos;
  int ncols,*tmp,count,nngh;

  ncols = 10;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal = new int[mesh.ntriangles];
  ordering->pointer  = new int[mesh.ntriangles + 1];
  ordering->nrows=mesh.ntriangles;
  tmp                = new int[mesh.ntriangles * ncols];

/*----------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < mesh.ntriangles * ncols; n++)
    tmp[n] = -1;

  count = 0;
  for(m = 0; m < mesh.ntriangles; m++) {
    ordering->pointer[m]  = count;      /*pointer= column start index */
    ordering->cardinal[m] = mesh.triangles[m].nngh+1;
    tmp[m*ncols] = m;
/*----------------------------------------------------------------------
    first level of neighbours*/
    for(k = 0; k < mesh.triangles[m].nngh; k++) {
      tmp[m*ncols + k +1] = mesh.triangles[m].ngh[k];
      }
/*----------------------------------------------------------------------
    second level of neighbours*/
    for(k = 0; k < mesh.triangles[m].nngh; k++) {
      mm = mesh.triangles[m].ngh[k];
      for(l = 0; l < mesh.triangles[mm].nngh; l++) {
        mmm=mesh.triangles[mm].ngh[l];
        pos=array_index(&(tmp[m*ncols]), ordering->cardinal[m], mmm);
        if(pos==-1) {
          tmp[m*ncols + ordering->cardinal[m]] = mmm;
          ordering->cardinal[m]++;
          }
        }
      }
    count += ordering->cardinal[m];
    if(ordering->cardinal[m]>ncols) {
      check_error(-1, "matrix index overflow", __LINE__, __FILE__, 1);
      }
    }

  ordering->pointer[mesh.ntriangles] = count;

  ordering->incidence = new int[count];

/*----------------------------------------------------------------------
  re-ordering of the column index*/
  for(m = 0; m < mesh.ntriangles; m++) {
  iterate:
    for(k = 0; k < ordering->cardinal[m]; k++) {
      n = tmp[m * ncols + k];
      for(l = ordering->cardinal[m] - 1; l > k; l--) {
        if(tmp[m * ncols + l] < n) {
          tmp[m * ncols + k] = tmp[m * ncols + l];
          tmp[m * ncols + l] = n;
          if(n == -1) {
            check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
            }
          goto iterate;
          }
        }
      }
    }

  for(n = 0; n < mesh.ntriangles; n++) {
    for(k = 0; k < ordering->cardinal[n]; k++) {
/* *----------------------------------------------------------------------------
      incidence= row index */
      ordering->incidence[ordering->pointer[n] + k] = tmp[n * ncols + k];
      if(ordering->incidence[ordering->pointer[n] + k] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }

/* *-----------------------------------------------------------------------
  some testing*/
  for(n = 0; n <count; n++) {
    if(ordering->incidence[n] == -1) {
      check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
      }
    }
/* *-----------------------------------------------------------------------
  symmetry testing*/
  for(n = 0; n < mesh.ntriangles; n++) {
    for(k = 0; k < ordering->cardinal[n]; k++) {
/* *----------------------------------------------------------------------------
      incidence= row index */
      m=ordering->incidence[ordering->pointer[n] + k];
      pos=array_index(&(tmp[m*ncols]), ordering->cardinal[m], n);
      if(pos== -1) {
        check_error(-1, "matrix symmetry anomaly", __LINE__, __FILE__, 1);
        }
      }
    }
  zaparr(tmp);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_LGP0xNCP1(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  LGP0 (rows) x ncP1 (column) based on first level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int i, m, n, n1, n2, n3, j, k, l;
  int opposite[2], row, col, ncols;
  double h;
  triangle_t *elt;
  int *tmp;
  int *incidence, *cardinal, *pointer;

  size_t count;

  ncols = 3;

  elt = mesh.triangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  cardinal = new int[mesh.ntriangles];
  pointer  = new int[mesh.ntriangles + 1];
  tmp      = new int[mesh.ntriangles * ncols];

/*----------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < mesh.ntriangles * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < mesh.ntriangles; n++) {
    pointer[n] = count;      /*pointer= column start index */
    cardinal[n] = 3;
    count += cardinal[n];
    for(k = 0; k <  3; k++) {
      tmp[n*ncols + k] = mesh.triangles[n].edges[k];
      }
    }
  pointer[n] = count;

  incidence = new int[count];

/*----------------------------------------------------------------------
  re-ordering of the column index*/
  for(n = 0; n < mesh.ntriangles; n++) {
  iterate:
    for(k = 0; k < cardinal[n]; k++) {
      m = tmp[n * ncols + k];
      for(l = cardinal[n] - 1; l > k; l--) {
//        if(tmp[n * ncols + l] == -1) continue;
        if(tmp[n * ncols + l] < m) {
          tmp[n * ncols + k] = tmp[n * ncols + l];
          tmp[n * ncols + l] = m;
          if(m == -1) {
            check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
            }
          goto iterate;
          }
        }
      }
    }

  for(n = 0; n < mesh.ntriangles; n++) {
    for(k = 0; k < cardinal[n]; k++) {
       /*incidence= row index */
      incidence[pointer[n] + k] = tmp[n * ncols + k];
      if(incidence[pointer[n] + k] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }

  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows=mesh.ntriangles;

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_LGP0xLGP1(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  LGP0 (rows) x LGP1 (column) based on first level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int i, m, n, n1, n2, n3, j, k, l;
  int opposite[2], row, col, ncols;
  double h;
  triangle_t *elt;
  int *tmp;
  int *incidence, *cardinal, *pointer;

  size_t count;

  ncols = 3;

  elt = mesh.triangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  cardinal = new int[mesh.ntriangles];
  pointer  = new int[mesh.ntriangles + 1];
  tmp      = new int[mesh.ntriangles * ncols];

/*----------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < mesh.ntriangles * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < mesh.ntriangles; n++) {
    pointer[n] = count;      /*pointer= column start index */
    cardinal[n] = 3;
    count += cardinal[n];
    for(k = 0; k <  3; k++) {
      tmp[n*ncols + k] = mesh.triangles[n].vertex[k];
      }
    }
  pointer[n] = count;

  incidence = new int[count];

/*----------------------------------------------------------------------
  re-ordering of the column index*/
  for(n = 0; n < mesh.ntriangles; n++) {
  iterate:
    for(k = 0; k < cardinal[n]; k++) {
      m = tmp[n * ncols + k];
      for(l = cardinal[n] - 1; l > k; l--) {
//        if(tmp[n * ncols + l] == -1) continue;
        if(tmp[n * ncols + l] < m) {
          tmp[n * ncols + k] = tmp[n * ncols + l];
          tmp[n * ncols + l] = m;
          if(m == -1) {
            check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
            }
          goto iterate;
          }
        }
      }
    }

  for(n = 0; n < mesh.ntriangles; n++) {
    for(k = 0; k < cardinal[n]; k++) {
       /*incidence= row index */
      incidence[pointer[n] + k] = tmp[n * ncols + k];
      if(incidence[pointer[n] + k] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }

  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows=mesh.ntriangles;

  int status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,mesh.ntriangles);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_LGP1xLGP0(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  LGP1 (rows) x LGP0 (columns) based on first level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int i, m, n, n1, n2, n3, j, k, l;
  int opposite[2], row, nrows, col, ncols;
  double A, hC, hc, temp1, temp2;
  double h;
  triangle_t *elt;
  double *M;
  int *tmp;
  int *incidence, *cardinal, *pointer;

  size_t count;

  nrows = -1;
  ncols = mesh.nnghm;

  elt = mesh.triangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  cardinal = new int[mesh.nvtxs];
  pointer  = new int[mesh.nvtxs + 1];
  tmp      = new int[mesh.nvtxs * ncols];

/*----------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < mesh.nvtxs * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < mesh.nvtxs; n++) {
    pointer[n]  = count;      /*pointer= column start index */
    cardinal[n] = mesh.vertices[n].nelmts;
    if(cardinal[n]>=ncols) {
      check_error(-1, "matrix index arrays overflow", __LINE__, __FILE__, 1);
      }
    count += cardinal[n];
    for(k = 0; k < mesh.vertices[n].nelmts; k++) {
      tmp[n*ncols + k] = mesh.vertices[n].elmts[k];
      }
    }
  pointer[n] = count;

  incidence = new int[count];

/*----------------------------------------------------------------------
  re-ordering of the column index*/
  for(n = 0; n < mesh.nvtxs; n++) {
  iterate:
    for(k = 0; k < cardinal[n]; k++) {
      m = tmp[n * ncols + k];
      for(l = cardinal[n] - 1; l > k; l--) {
        if(tmp[n * ncols + l] < m) {
          tmp[n * ncols + k] = tmp[n * ncols + l];
          tmp[n * ncols + l] = m;
          if(m == -1) {
            check_error(-1, "matrix index arrays overflow", __LINE__, __FILE__, 1);
            }
          goto iterate;
          }
        }
      }
    }

  for(n = 0; n < mesh.nvtxs; n++) {
    for(k = 0; k < cardinal[n]; k++) {
       /*incidence= row index */
      incidence[pointer[n] + k] = tmp[n * ncols + k];
      if(incidence[pointer[n] + k] == -1) {
        check_error(-1, "matrix index arrays overflow", __LINE__, __FILE__, 1);
        }
      }
    }

  zaparr(tmp);


  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows=mesh.nvtxs;

/* *----------------------------------------------------------------------
  check matrix duplicated connections*/
  int status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,mesh.nvtxs);

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_NCP1xLGP0(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  ncP1 (rows) x LGP0 (columns) based on first level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int i, m, n, n1, n2, n3, j, k, l;
  int opposite[2], row, nrows, col, ncols;
  double A, hC, hc, temp1, temp2;
  double h;
  triangle_t *elt;
  double *M;
  int *tmp;
  int *incidence, *cardinal, *pointer;

  size_t count;

  nrows = 2;
  ncols = 2;

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
    cardinal[n] = mesh.edges[n].nshared;
    count += cardinal[n];
    for(k = 0; k < mesh.edges[n].nshared; k++) {
      tmp[n*ncols + k] = mesh.edges[n].shared[k];
      }
    }
  pointer[n] = count;

  incidence = new int[count];

/*----------------------------------------------------------------------
  re-ordering of the column index*/
  for(n = 0; n < mesh.nedges; n++) {
  iterate:
    for(k = 0; k < cardinal[n]; k++) {
      m = tmp[n * ncols + k];
      for(l = cardinal[n] - 1; l > k; l--) {
        if(tmp[n * ncols + l] < m) {
          tmp[n * ncols + k] = tmp[n * ncols + l];
          tmp[n * ncols + l] = m;
          if(m == -1) {
            check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
            }
          goto iterate;
          }
        }
      }
    }

  for(n = 0; n < mesh.nedges; n++) {
    for(k = 0; k < cardinal[n]; k++) {
       /*incidence= row index */
      incidence[pointer[n] + k] = tmp[n * ncols + k];
      if(incidence[pointer[n] + k] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }

  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows=mesh.nedges;

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *MassPMatrix_P0(mesh_t mesh, int *incidence, int *cardinal, int *pointer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *----------------------------------------------------------------------

  Integrale version of P0 mass matrix

  Packed matrix where neighbours list as row index

----------------------------------------------------------------------*/
  int count, i, j, k, m, n, status;
  int gauss_n;
  int n1, n2;
  int ni, nj, row, column;
  double *M, p[3],q[3],c[3];
  double C, area;
  double beta[3];
  double gauss_x[16], gauss_y[16], gauss_w[16];

  count = pointer[mesh.ntriangles];

/*----------------------------------------------------------------------
  allocate packed matrix*/
  M = new double[count];

  for(n = 0; n < count; n++)
    M[n] = 0;

  for(m = 0; m < mesh.ntriangles; m++) {
    for(i = 0; i < 3; i++) {
      n = mesh.triangles[m].vertex[i];
      c[i]=mesh.vertices[n].c;
      }
    area = mesh.triangles[m].Area;
    M[m] += 2. * area * fe_integraleLGP1_2D(c);
    }
  return (M);
}

