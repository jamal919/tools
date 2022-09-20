
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

  int init_packed_NCP1xNCP1(mesh_t mesh, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  ncP1 (rows) x ncP1 (columns) based on first level neighbouring relationships

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

  nrows =-1;
  ncols = 5;

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
    cardinal[n] = mesh.edges[n].nngh+1;
    count += cardinal[n];
    tmp[n*ncols]=n;
    for(k = 0; k < mesh.edges[n].nngh; k++) {
      tmp[n*ncols + k +1] = mesh.edges[n].ngh[k];
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
  ordering->nrows     = mesh.nedges;

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_NCP1xNCP1_02(mesh_t mesh, ordering_t *ordering, int depth)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  ncP1 (rows) x ncP1 (columns) based on n levels neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int i, m, n, nn, j, k, l;
  int row, nrows, col, ncols,loop,first,last,additional;
  int *tmp;
  int *incidence, *cardinal, *pointer;

  size_t count;

  nrows = -1;
  ncols =  5;
  additional=4;
  for(loop=1;loop<depth;loop++) {
    ncols+=additional*2;
    additional=additional*2;
    }

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
    tmp[n*ncols] = n;
    cardinal[n] = 1;
    first=0;
    for(loop=0;loop<depth;loop++) {
      last=cardinal[n];
      for(k=first;k<last;k++) {
        m=tmp[n*ncols+k];
        for(l=0;l<mesh.edges[m].nngh;l++) {
          nn=mesh.edges[m].ngh[l];
 //         if(pos(nn,&(tmp[n*ncols]),cardinal[n])==-1) {
          if(array_index(&(tmp[n * ncols]),cardinal[n],nn) == -1) {
            if(cardinal[n]==ncols) {
              check_error(-1, "matrix index overflow", __LINE__, __FILE__, 1);
              }
            tmp[n*ncols+cardinal[n]] = nn;
            cardinal[n]++;
            }
          }
        }
      first=last;
      }
    count += cardinal[n];
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
  ordering->nrows     = mesh.nedges;

  int status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,mesh.nedges);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_NCP1xLGP0_02(mesh_t mesh, ordering_t *ordering,int depth)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  ncP1 (rows) x LGP0 (columns) based on 2 level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int i, m, mm, n, n1, n2, n3, j, k, l;
  int opposite[2], row, nrows, col, ncols, loop, first, last, additional;
  double A, hC, hc, temp1, temp2;
  double h;
  triangle_t *elt;
  double *M;
  int *tmp;
  int *incidence, *cardinal, *pointer;

  size_t count;

  nrows = 2;
  ncols = 6;
  additional=4;
  for(loop=1;loop<depth;loop++) {
    ncols+=additional*2;
    additional=additional*2;
    }

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
    for(k = 0; k < mesh.edges[n].nshared; k++) {
      tmp[n*ncols + k] = mesh.edges[n].shared[k];
      }
    first=0;
    for(loop=0;loop<depth;loop++) {
      last=cardinal[n];
      for(k = first; k < last; k++) {
        m=tmp[n*ncols + k];
        for(l = 0; l < mesh.triangles[m].nngh; l++) {
          mm=mesh.triangles[m].ngh[l];
 //         if(pos(mm,&(tmp[n*ncols]),cardinal[n])==-1) {
          if(array_index(&(tmp[n * ncols]),cardinal[n],mm) == -1) {
            tmp[n*ncols + cardinal[n]] = mm;
            if(cardinal[n]==ncols) {
              check_error(-1, "matrix index overflow", __LINE__, __FILE__, 1);
              }
            cardinal[n]++;
            }
          }
        }
      first=last;
      }
    count += cardinal[n];
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
  ordering->nrows     = mesh.nedges;

  int status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,mesh.nedges);


  return (0);
}

