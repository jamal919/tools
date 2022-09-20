
/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#include <stdio.h>
#include <string.h>
#include <time.h>

#include "tools-structures.h"

#include "fe.h"

//int solver_type;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int syspack_ordering(int **neighbours, int *card, int neq, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------

  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometric matrix symmetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  ----------------------------------------------------------------------*/
{
  int i, m, mm, first, n, n1, n2, n3, j, k, l, row, nrows, col, ncols;
  double h;
  double *M;

  size_t count;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal = new int[neq];
  ordering->pointer  = new int[neq + 1];

  count = 0;
  for(n = 0; n < neq; n++) {
    (ordering->pointer)[n] = count;      /*pointer= column start index */
    (ordering->cardinal)[n] = card[n];
    count += (ordering->cardinal)[n];
    }
  (ordering->pointer)[n] = count;

  ordering->incidence = new int[count];
  for(n = 0; n < neq; n++) {
    first=(ordering->pointer)[n];
    for(k = 0; k < card[n]; k++) {
      (ordering->incidence)[first + k]=neighbours[n][k];
      }
    }

/*----------------------------------------------------------------------
  re-ordering of the column index*/
  for(n = 0; n < neq; n++) {
    first=(ordering->pointer)[n];
  iterate:
    for(k = 0; k < (ordering->cardinal)[n]; k++) {
      m = (ordering->incidence)[first + k];
      for(l = (ordering->cardinal)[n] - 1; l > k; l--) {
        mm= (ordering->incidence)[first + l];
        if(mm < m) {
          (ordering->incidence)[first + k] = mm;
          (ordering->incidence)[first + l] = m;
          goto iterate;
          }
        }
      }
    }

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *syspack_matrix(double **in, int **neighbours, int *card, int neq, ordering_t ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count, i, j, k, l, m, n, status;
  int gauss_n;
  int ni, nj, row, column;
  double *M, p[3],q[3],r[3];

  count = ordering.pointer[neq];

/*----------------------------------------------------------------------
  allocate packed matrix*/
  M = new double[count];

  for(n = 0; n < count; n++)
    M[n] = 0;

  for(m = 0; m < neq; m++) {
    for(l = 0; l < card[m]; l++) {
      n=neighbours[m][l];
      for(k = 0; k < ordering.cardinal[m]; k++) {
        if(ordering.incidence[ordering.pointer[m] + k] == n) {
          row = k;
          break;
          }
        }
      if(k==ordering.cardinal[m]) {
        printf("troubles\n");
        }
      M[ordering.pointer[m] + row] = in[n][l];
      }
    }
  for(n = 0; n < count; n++){
    if(M[n] == 0) {
      printf("troubles\n");
      }
    }

  return (M);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int syspack_factorisation(int neq, double *packed,ordering_t ordering,triplet *Triplet, solver_t **Solver)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes

  Check :

  Notes

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  int i, lda, status;
  int *colstr, *rowind;
  char mtxtyp[2], cpivot, ordmthd;
  int outunt, msglvl = 0, error, nrhs, ldrhs, neqns, nz;
  float *rhs;

  char c;
  int *Tj;
  FILE *out;
  clock_t d1, d2;
  double duree1, duree2;
  int count, m, k, n;
  char **args;
  triplet *M1;
  solver_t *solveur1;

  // init triplet
  Triplet->comptype=CSC;
  Triplet->base=0;

  Triplet->stype=0;

  Triplet->nnz=ordering.pointer[neq];
  Triplet->nrow=neq;
  Triplet->ncol=neq;

  Triplet->i = ordering.pointer;
  Triplet->j = ordering.incidence;
  Triplet->x = packed;

  M1 = Triplet;

  status=0;
  solveur1=(solver_t *) init_solver_obsolete("SpDOMESTIC", M1->stype);

  status = factorize(solveur1,M1);
  if (status != 0) {
    __ERR_BASE_LINE__("exiting\n");exit(-1);
    }
  *Solver = solveur1;

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int syspack_solve(triplet *Triplet, solver_t *Solver, double *B)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, i, ndiff, ierr;
  int hbw,neqs;

  char cjob = 'N';
  void *Numeric;
  clock_t d1, d2;
  double duree1, duree2;
  /*Vec biter,xiter; */
  triplet *M1;
  solver_t *solveur1;

  M1 =  Triplet;
  solveur1=Solver;

  d1 = clock();
  ierr = solve(solveur1,B,1,false);
  d2 = clock();
  duree1 = d2 - d1;

  return(ierr);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int syspack_inverse(double *in, double *out, int neq, hypermatrix_t matrix)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n,status;

  for(n = 0; n < neq; n++)
    out[n] = in[n];

  status = syspack_solve(&(matrix.t),matrix.s,out);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int syspack_init(double **in,int **neighbours, int *card, int neq, hypermatrix_t *matrix)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n, n1, n2, status;

  matrix->ordering=new ordering_t;
  status=syspack_ordering(neighbours, card, neq,  matrix->ordering);
  matrix->packed=syspack_matrix(in, neighbours, card, neq, *(matrix->ordering));
  status=syspack_factorisation(neq, matrix->packed,*(matrix->ordering),&(matrix->t),&(matrix->s));
  return(status);
}


