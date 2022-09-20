
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
#include "zapper.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_LGP1xLGP1_0(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  P1 x P1, based on first level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int i, m, n, n1, n2, n3, j, k, l, row, nrows, col, ncols;
  double A, hC, hc, temp1, temp2;
  double h;
  triangle_t *elt;
  int nndes, nelts;
  double *M;
  int *tmp;

  size_t count;

  nndes = mesh.nvtxs;
  nelts = mesh.ntriangles;
  nrows = mesh.nnghm + 1;
  ncols = mesh.nnghm + 1;

  elt = mesh.triangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal = new int[nndes];
  ordering->pointer  = new int[nndes + 1];
  tmp       = new int[nndes * ncols];

  for(n = 0; n < nndes * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  allocate arrays*/
  for(n = 0; n < nndes; n++) {
    tmp[n * ncols] = n;
    for(k = 0; k < mesh.vertices[n].nngh; k++)
      tmp[n * ncols + k + 1] = mesh.vertices[n].ngh[k];
      }

  count = 0;
  for(n = 0; n < nndes; n++) {
    (ordering->pointer)[n] = count;      /*pointer= column start index */
    (ordering->cardinal)[n] = mesh.vertices[n].nngh + 1;
    count += (ordering->cardinal)[n];
    }
  (ordering->pointer)[n] = count;

  ordering->incidence = new int[count];

/*----------------------------------------------------------------------
  re-ordering of the column index*/
  for(n = 0; n < nndes; n++) {
  iterate:
    for(k = 0; k < (ordering->cardinal)[n]; k++) {
      m = tmp[n * ncols + k];
      for(l = (ordering->cardinal)[n] - 1; l > k; l--) {
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

  for(n = 0; n < nndes; n++) {
    for(k = 0; k < (ordering->cardinal)[n]; k++) {
       /*incidence= row index */
      (ordering->incidence)[(ordering->pointer)[n] + k] = tmp[n * ncols + k];
      if((ordering->incidence)[(ordering->pointer)[n] + k] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }

  ordering->nrows = nndes;
  zaparr(tmp);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int init_packed_LGP1xLGP1_3(mesh_t mesh, ordering_t *ordering, paire_t *distant_connexion, int ndistant_connexion)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------

  Compressed Row storage/Compressed Sparse Row

  P1 x P1, based on first level neighbouring relationships + extra connexions

----------------------------------------------------------------------*/
{
  int i, m, n, n1, n2, n3, j, k, l, row, nrows, col, ncols;
  int extra;
  triangle_t *elt;
  int nndes, nelts;
  int *tmp;
  int *nadditional;

  size_t count;

  nndes = mesh.nvtxs;
  nelts = mesh.ntriangles;
  
  nadditional=new int[mesh.nvtxs];
  for(n = 0; n < mesh.nvtxs; n++) {
    nadditional[n]=0;
    }

  for(i = 0; i < ndistant_connexion; i++) {
    nadditional[distant_connexion[i].value[0]]++;
    nadditional[distant_connexion[i].value[1]]++;
    }

  extra=0;
  for(n = 0; n < nndes; n++)
    extra = MAX(extra,nadditional[n]);

  nrows = mesh.nnghm + 1 + extra;
  ncols = mesh.nnghm + 1 + extra;

  elt = mesh.triangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal = new int[nndes];
  ordering->pointer  = new int[nndes + 1];
  ordering->nrows = nndes;
  tmp       = new int[nndes * ncols];

/*----------------------------------------------------------------------
  temporary incidence array with over-dimensionned size*/
  for(n = 0; n < nndes * ncols; n++)
    tmp[n] = -1;

  for(n = 0; n < nndes; n++) {
    tmp[n * ncols] = n;
    for(k = 0; k < mesh.vertices[n].nngh; k++) {
      tmp[n * ncols + k + 1] = mesh.vertices[n].ngh[k];
      }
/*
    for(l = 0; l < nadditional[n]; l++) {
      tmp[n * ncols + k + 1 + l] = additional[n][l];
      }
*/
    }
  
  for(n = 0; n < mesh.nvtxs; n++) {
    nadditional[n]=0;
    }
  
  for(i = 0; i < ndistant_connexion; i++) {
    n1=distant_connexion[i].value[0];
    n2=distant_connexion[i].value[1];
    tmp[n1 * ncols + mesh.vertices[n1].nngh + nadditional[n1] + 1] = n2;
    tmp[n2 * ncols + mesh.vertices[n2].nngh + nadditional[n2] + 1] = n1;
    nadditional[n1]++;
    nadditional[n2]++;
    }

/*----------------------------------------------------------------------
  build pointer and cardinal arrays*/
  count = 0;
  for(n = 0; n < nndes; n++) {
    ordering->pointer[n] = count;      /*pointer= column start index */
    ordering->cardinal[n] = mesh.vertices[n].nngh + nadditional[n] + 1;
    count += ordering->cardinal[n];
    }
  ordering->pointer[n] = count;

  ordering->incidence = new int[count];

/*----------------------------------------------------------------------
  re-ordering of the column index*/
  for(n = 0; n < nndes; n++) {
  iterate:
    for(k = 0; k < ordering->cardinal[n]; k++) {
      m = tmp[n * ncols + k];
      for(l = ordering->cardinal[n] - 1; l > k; l--) {
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

  for(n = 0; n < nndes; n++) {
    for(k = 0; k < ordering->cardinal[n]; k++) {
       /*incidence= row index */
      ordering->incidence[ordering->pointer[n] + k] = tmp[n * ncols + k];
      if(ordering->incidence[ordering->pointer[n] + k] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }

  zaparr(tmp);
  zaparr(nadditional);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_LGP1xLGP1_1(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  LGP1 x LGP1 based on first and second level neighbouring relationships

  To be used with LGP1 velocities

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

  nrows = 2*mesh.nnghm + 1;
  ncols = 4*mesh.nnghm + 1;

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

  count = 0;
  for(n = 0; n < mesh.nvtxs; n++) {
    pointer[n] = count;      /*pointer= column start index */
    cardinal[n] = mesh.vertices[n].nngh + 1;
    count += cardinal[n];
    tmp[n * ncols] = n;
/* *----------------------------------------------------------------------
    first level of neighbours*/
    for(k = 0; k < mesh.vertices[n].nngh; k++) {
      tmp[n * ncols + k + 1] = mesh.vertices[n].ngh[k];
      }
/* *----------------------------------------------------------------------
    2nd level of neighbours*/
    for(k = 0; k < mesh.vertices[n].nngh; k++) {
      m=mesh.vertices[n].ngh[k];
      for(l = 0; l < mesh.vertices[m].nngh; l++) {
        if(array_index(&(tmp[n * ncols]),cardinal[n],mesh.vertices[m].ngh[l]) == -1) {
          tmp[n * ncols + cardinal[n]] = mesh.vertices[m].ngh[l];
          cardinal[n]++;
          count++;
          }
        }
      }
    if(cardinal[n]>ncols) {
      check_error(-1, "matrix index arrays overflow", __LINE__, __FILE__, 1);
      }
    }


  count = 0;
  for(n = 0; n < mesh.nvtxs; n++) {
    pointer[n] = count;      /*pointer= column start index */
    count += cardinal[n];
    }



  pointer[mesh.nvtxs] = count;
  incidence = new int[count];

/*----------------------------------------------------------------------
  re-ordering of the column index*/
  for(n = 0; n < mesh.nvtxs; n++) {
  iterate:
    for(k = 0; k < cardinal[n]; k++) {
      m = tmp[n * ncols + k];
      for(l = cardinal[n] - 1; l > k; l--) {
        if(tmp[n * ncols + l] == -1) continue;
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

  for(n = 0; n < mesh.nvtxs; n++) {
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
  ordering->nrows=mesh.nvtxs;

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_LGP1xLGP1_2(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  P1 x P1 based on first and second level neighbouring relationships

  To be used with NCP1 velocities

  ----------------------------------------------------------------------*/
{
  int i, m, n, n1, n2, n3, j, k, l;
  int opposite[2], row, nrows, col, ncols;
  double A, hC, hc, temp1, temp2;
  double h;
  triangle_t *elt;
  double *M;
  int add,*tmp;
  int *incidence, *cardinal, *pointer;

  size_t count;

  nrows = 2*mesh.nnghm + 1;
  ncols = 2*mesh.nnghm + 1;

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
    pointer[n] = count;      /*pointer= column start index */
    cardinal[n] = mesh.vertices[n].nngh + 1;
    count += cardinal[n];
    tmp[n * ncols] = n;
    for(k = 0; k < mesh.vertices[n].nngh; k++) {
      tmp[n * ncols + k + 1] = mesh.vertices[n].ngh[k];
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
      j=3-(mesh.edges[n].vindex[k][0]+mesh.edges[n].vindex[k][1]);
      opposite[k]=mesh.triangles[m].vertex[j];
      }
    n1=opposite[0];
    n2=opposite[1];
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
      printf("alert\n");
      }
    tmp[n1 * ncols + k] = n2;
    cardinal[n1]++;
    count++;
    k=cardinal[n2];
    if(k > ncols-1) {
      printf("alert\n");
      }
    tmp[n2 * ncols + k] = n1;
    cardinal[n2]++;
    count++;
    }

  count = 0;
  for(n = 0; n < mesh.nvtxs; n++) {
    pointer[n] = count;      /*pointer= column start index */
    count += cardinal[n];
    }

  pointer[mesh.nvtxs] = count;
  incidence = new int[count];
  for(n = 0; n < count; n++) {
    incidence[n]=-1;
    }

/* *----------------------------------------------------------------------
  check matrix duplicated connections*/
//   for(n = 0; n < mesh.nvtxs; n++) {
//     for(k = 0; k < cardinal[n]-1; k++) {
//       m=tmp[n * ncols + k];
//       if(tmp[n * ncols + k+1]==m) {
//         printf("anomaly, duplicated node=%d, index=%d\n",n,k);
//         }
//       }
//     }

/*----------------------------------------------------------------------
  re-ordering of the column index*/
  for(n = 0; n < mesh.nvtxs; n++) {
  iterate:
    for(k = 0; k < cardinal[n]; k++) {
      m = tmp[n * ncols + k];
      for(l = cardinal[n] - 1; l > k; l--) {
        if(tmp[n * ncols + l] == -1) continue;
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

  for(n = 0; n < mesh.nvtxs; n++) {
    if(cardinal[n] == 0) {
      printf("anomaly, null cardinal for node=%d \n",n);
      }
    for(k = 0; k < cardinal[n]; k++) {
       /*incidence= row index */
      incidence[pointer[n] + k] = tmp[n * ncols + k];
      if(incidence[pointer[n] + k] == -1) {
        printf("anomaly, invalid connection at node=%, index=%d\n",n,k);
        }
      }
    }

//   for(n = 0; n < count; n++) {
//     if(incidence[n] == -1) {
//       printf("anomaly, invalid connection at count=%d\n",count);
//       }
//     }

/* *----------------------------------------------------------------------
  check matrix geometry*/
//   for(n = 0; n < mesh.nvtxs; n++) {
//     for(k = 0; k < cardinal[n]; k++) {
//       m=incidence[pointer[n] + k];
//       for(l = 0; l < cardinal[m]; l++) {
//         if(incidence[pointer[m] + l]==n) goto ok;
//         }
//       printf("anomaly, non-symmetric node=%, index=%d\n",n,k);
// ok:
//       count=0;
//       }
//     }
/* *----------------------------------------------------------------------
  check matrix duplicated connections*/
//   for(n = 0; n < mesh.nvtxs; n++) {
//     for(k = 0; k < cardinal[n]-1; k++) {
//       m=incidence[pointer[n] + k];
//       if(incidence[pointer[n] + k+1]==m) {
//         printf("anomaly, duplicated node=%d, index=%d, connected=%d\n",n,k,m);
//         }
//       }
//     }
  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows=mesh.nvtxs;

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_NCP1xLGP1(mesh_t mesh, ordering_t *ordering)

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
  int i, m, n, n1, n2, n3, j, k, l;
  int opposite[2], row, nrows, col, ncols;
  double A, hC, hc, temp1, temp2;
  double h;
  triangle_t *elt;
  double *M;
  int *tmp;
  int *incidence, *cardinal, *pointer;

  size_t count;

  nrows = 4;
  ncols = 4;

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
    cardinal[n] = mesh.edges[n].nvertex;
    count += cardinal[n];
    for(k = 0; k < mesh.edges[n].nvertex; k++) {
      tmp[n*ncols + k] = mesh.edges[n].vertex[k];
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

  int init_packed_LGP1xNCP1(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  P1 x ncP1 based on first level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int i, m, n, n1, n2, n3, j, k, l;
  int opposite[2], row, nrows, col, ncols;
  double h;
  triangle_t *elt;
  int *tmp;
  int *incidence, *cardinal, *pointer;

  size_t count;

  nrows = 2*mesh.nnghm;
  ncols = 2*mesh.nnghm;

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
    pointer[n] = count;      /*pointer= column start index */
    cardinal[n] = mesh.vertices[n].nedges;
    count += cardinal[n];
    for(k = 0; k <  mesh.vertices[n].nedges; k++) {
      tmp[n*ncols + k] = mesh.vertices[n].edges[k];
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

  for(n = 0; n < mesh.nvtxs; n++) {
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

  int QLP1_node(mesh_t mesh, int vertex, int level)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  ----------------------------------------------------------------------*/
{
  int m;

  m=vertex*mesh.nlevels+level;

  return(m);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_QLP1_3D(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  QLP1 x QLP1, based on first level prism's neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int i, m, n, nn, n1, n2, n3, j, k, l, row, nrows, col, ncols, pos;
  double A, hC, hc, temp1, temp2;
  double h;
  triangle_t *elt;
  int nelts;
  double *M;
  int *tmp;

  size_t count;

/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

  15/04/2008 :

  eacho node x level

  n = (node index) x nlevels + (level index)

  horizontal connection betwen nodes' neighbours
  vertical connection betwen 3 successive levels

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  nelts = mesh.ntriangles;

  nrows = mesh.nnghm + 1;
  ncols = (mesh.nnghm + 1)*3;

  elt = mesh.triangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal = new int[mesh.nvtxs*mesh.nlevels];
  ordering->pointer  = new int[mesh.nvtxs*mesh.nlevels + 1];
  ordering->nrows=mesh.nvtxs;
  tmp                = new int[mesh.nvtxs*mesh.nlevels * ncols];

  for(n = 0; n < mesh.nvtxs * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  allocate arrays*/
  for(n = 0; n < mesh.nvtxs; n++) {
    for(l = 0; l < mesh.nlevels; l++) {
      pos=QLP1_node(mesh,n,l);
      tmp[pos] = QLP1_node(mesh,n,l);
      for(k = 0; k < mesh.vertices[n].nngh; k++) {
        nn=mesh.vertices[n].ngh[k];
        tmp[pos + k + 1] = QLP1_node(mesh,nn,l);
        }
      ordering->cardinal[m] = mesh.vertices[n].nngh + 1;
      }
    for(l = 1; l < mesh.nlevels; l++) {
      m=QLP1_node(mesh,n,l);
      pos=m+ordering->cardinal[m];
      tmp[pos] = QLP1_node(mesh,n,l-1);
      for(k = 0; k < mesh.vertices[n].nngh; k++) {
        nn=mesh.vertices[n].ngh[k];
        tmp[pos + k + 1] = QLP1_node(mesh,nn,l-1);
        }
      ordering->cardinal[m] = mesh.vertices[n].nngh + 1;
      }
    for(l = 0; l < mesh.nlevels-1; l++) {
      m=QLP1_node(mesh,n,l+1);
      pos=m+ordering->cardinal[m];
      tmp[pos] = QLP1_node(mesh,n,l+1);
      for(k = 0; k < mesh.vertices[n].nngh; k++) {
        nn=mesh.vertices[n].ngh[k];
        tmp[pos + k + 1] = QLP1_node(mesh,nn,l+1);
        }
      ordering->cardinal[m] = mesh.vertices[n].nngh + 1;
      }
   }

  count = 0;
  for(n = 0; n < mesh.nvtxs*mesh.nlevels; n++) {
    ordering->pointer[n] = count;      /*pointer= column start index */
    count += (ordering->cardinal)[n];
    }
  ordering->pointer[n] = count;

  ordering->incidence = new int[count];

/*----------------------------------------------------------------------
  re-ordering of the column index*/
  for(n = 0; n < mesh.nvtxs*mesh.nlevels; n++) {
  iterate:
    for(k = 0; k < ordering->cardinal[n]; k++) {
      m = tmp[n * ncols + k];
      for(l = ordering->cardinal[n] - 1; l > k; l--) {
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

  for(n = 0; n < mesh.nvtxs*mesh.nlevels; n++) {
    for(k = 0; k < ordering->cardinal[n]; k++) {
       /*incidence= row index */
      ordering->incidence[ordering->pointer[n] + k] = tmp[n * ncols + k];
      if(ordering->incidence[(ordering->pointer)[n] + k] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }

  zaparr(tmp);

  return (0);
}

