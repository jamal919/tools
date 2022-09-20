
/*******************************************************************************

  T-UGO tools, 2006-2019

  Unstructured Ocean Grid initiative

*******************************************************************************/

#include "constants.h"

#include "tools-structures.h"
#include "matrix.h"
#include "solverlib.h"
#include "maths.h"


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  
  warning:
  
  matrix must follow fortran index order, i.e. fastest index being rox index
  
  
  mat[i,j]=mat[j*nrows+i]


xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_print(double *matrix, int nrow, int ncol)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j;
  
  for(i=0;i<nrow;i++) {
    for (j=0;j<ncol;j++) {
      printf("%lf ",matrix[j*nrow+i]);
      }
    printf("\n");
    }
  printf("\n");
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_inverse(double *matrix, int neq, double* & inverse, bool preserve, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, /*l,*/ lda, ldb, status;
  double *work;
  int neq2=square(neq);
  
  int pivot[neq];
  
  lda=neq;
  ldb=neq;
  
  if(preserve) {
    work=new double[neq2];
    for(k = 0; k < neq2; k++){
      work[k] = matrix[k];
      }
    }
  else work=matrix;

  for(k = 0; k < neq2; k++){
    if(isnan(matrix[k])==1) TRAP_ERR_EXIT(-1, "nan found in matrix, abort matrix_inverse\n");
    }
  
/*------------------------------------------------------------------------------
  warning : poc_getrf DOES modify the input matrix */ 
  status=poc_getrf(neq, work, pivot);
  if(status!=0) {
    printf("%s : matrix (%d,%d) is singular\n", __func__, neq, neq);
    status=matrix_print(work, neq, neq);
    TRAP_ERR_EXIT(-1, "factorization failed\n");
    }
  if (inverse==0) inverse=new double[neq2];
  
  aset(inverse,neq2,0.);

  for(k = 0; k < neq; k++){
    inverse[k*neq+k] = 1.0;
    }

  int nrhs=neq;
  
  status=poc_getrs(neq, nrhs, work, pivot, inverse);
  if(status!=0) TRAP_ERR_EXIT(-1, "solver failed\n");
  
//   double sum=0;
//   l=0;
//   for(k = 0; k < neq; k++){
//     sum+=matrix[k*neq+l]*inverse[l*neq+k];
//     }
//   sum=0;
//   for(k = 0; k < neq; k++){
//     sum+=matrix[l*neq+k]*inverse[k*neq+l];
//     }
    
  if(preserve) {
    delete[] work;
    }
  
  if(verbose>=1) status=matrix_print(inverse, neq, neq);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_transpose(double *matrix, int nrow, int ncol, double* & transpose, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l;
  
  if (transpose==0) transpose=new double[nrow*ncol];
  
  for(k = 0; k < nrow*ncol; k++){
    transpose[k] = 0.0;
    }

  for(l = 0; l < nrow; l++){
    for(k = 0; k < ncol; k++){
      transpose[l*ncol+k] = matrix[k*nrow+l];
      }
    }
 
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int matrix_transpose_template(T & A, T *tA, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccLinearSystem_solve
 
  compute a transpose matrix from a sparse, possibly distributed matrix with
  CSR ordering
  
  MPI-compliant if A is complete
  
  NOT COMPLETED
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  int asymmetric;
  ordering_t *CSCordering;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  asymmetric= matrix_check_symmetry(*A.ordering);
  if(asymmetric==0) {
/*------------------------------------------------------------------------------
    if ordering is symmetric, CSR and CSC ordering tables are the same*/
    status= matrix_CSR2CSC(&A);
    }
  else {
/*------------------------------------------------------------------------------
    if ordering is not symmetric, we need a separate CSC ordering*/
    CSCordering=ordering_CSR2CSC(A.ordering);
    status= matrix_CSR2CSC(&A,CSCordering);
    A.ordering->destroy();
    A.ordering=CSCordering;
    if(debug) status= matrix_check_consistency(*CSCordering);
    }
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_transpose(hyperzmatrix_t & A, hyperzmatrix_t *tA, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=matrix_transpose_template(A, tA, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_transpose(hypermatrix_t & A, hypermatrix_t *tA, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=matrix_transpose_template(A, tA, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename P, typename Q> int matrix_product_template(P *A, Q *B, int nrowA, int ncolA, int nrowB, int ncolB, Q* & product, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, k, l;
  
  if(ncolA!=nrowB) return(-1);
      
  int nrow=nrowA;
  int ncol=ncolB;
  
  if (product==0) product=new Q[nrow*ncol];
  
  for(k = 0; k < nrow*ncol; k++) {
    product[k] = 0.0;
    }

  for(l = 0; l < ncol; l++) {
    for(k = 0; k < nrow; k++) {
      for(i = 0; i < ncolA; i++) {
        product[l*nrow+k] += A[i*nrowA+k]*B[l*nrowB+i];
        }
      }
    }
 
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_product(double *A, double *B, int nrowA, int ncolA, int nrowB, int ncolB, double* & product, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=matrix_product_template(A, B, nrowA, ncolA, nrowB, ncolB, product, verbose);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_product(double *A, double *B, int ndim, double* & product, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l, status;
  double *tmp=0;
  
  if(product==A) {
    status=matrix_product(A, B, ndim, ndim, ndim, ndim, tmp, verbose);
    for(l = 0; l < ndim; l++){
      for(k = 0; k < ndim; k++){
        A[l*ndim+k]= tmp[l*ndim+k];
        }
      }
    delete[] tmp;
    }
  else if(product==B) {
    status=matrix_product(A, B, ndim, ndim, ndim, ndim, tmp, verbose);
    for(l = 0; l < ndim; l++){
      for(k = 0; k < ndim; k++){
        B[l*ndim+k]= tmp[l*ndim+k];
        }
      }
    delete[] tmp;
    }
  else {
    status=matrix_product(A, B, ndim, ndim, ndim, ndim, product, verbose);
    }
  
  if(status!=0) TRAP_ERR_EXIT(-1, "matrix product failed\n");
  
  if(verbose==1) status=matrix_print(product, ndim, ndim);
    
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename P, typename Q> int matrix_operation(P *A, Q *B, int nrowA, int ncolA, Q* & product, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=matrix_product_template(A, B, nrowA, ncolA, ncolA, 1, product, verbose);
 
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_operation(double *A, double *B, int nrowA, int ncolA, double* & product, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=matrix_product_template(A, B, nrowA, ncolA, ncolA, 1, product, verbose);
 
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_operation(double *A, complex<double> *B, int nrowA, int ncolA, complex<double>* & product, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=matrix_product_template(A, B, nrowA, ncolA, ncolA, 1, product, verbose);
 
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_least_square(double *matrix, int nrow, int ncol, double* & leastsquare, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  compute inv(tM x M) x tM
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k, status;
  int neq=ncol;
  double *transpose=0, *product=0, *inverse=0;
  
  if (leastsquare==0) leastsquare=new double[ncol*nrow];
  for(k = 0; k < ncol*nrow; k++){
    leastsquare[k] = 0.0;
    }
  
/**----------------------------------------------------------------------------
  compute tM x M (ncolxncol matrix)*/
  status=matrix_transpose(matrix, nrow, ncol,  transpose, verbose);
  
  status=matrix_product(transpose, matrix, ncol, nrow, nrow, ncol,  product, verbose);
  
/**----------------------------------------------------------------------------
  compute inverse (tM x M)  (ncolxncol matrix)*/
  status=matrix_inverse(product, neq,  inverse, false, verbose);
  
/**----------------------------------------------------------------------------
  compute inverse (tM x M) x tM  (ncol x nrow matrix)*/
  status=matrix_product(inverse, transpose, neq, neq, ncol, nrow,  leastsquare, verbose);
  
  if(verbose==1) status=matrix_print(leastsquare, ncol, nrow);

/**----------------------------------------------------------------------------
  verification */
  
  delete[] transpose;
  delete[] product;
  delete[] inverse;

  return(0);
}

