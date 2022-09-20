/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2017

  Unstructured Ocean Grid initiative

Developers:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
Contributors:
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France
  Frédéric Dupont    Canada environement, Canada

***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <time.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "solvers-functions.h"
#include "linear-gmres.h"

#if 1

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int vector_checkInf(double *A, int zdim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m;
  for(m = 0; m < zdim; m++) {
    double value=A[m];
    if(isinf(value)!=0) {
      TRAP_ERR_EXIT(-1, "nan encountered at n=%d\n",m);
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int vector_checkInf(complex<double> *A, int zdim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m;
  for(m = 0; m < zdim; m++) {
    double value=real(A[m]);
    if(isinf(value)!=0) {
      TRAP_ERR_EXIT(-1, "nan encountered at n=%d\n",m);
      }
    value=imag(A[m]);
    if(isinf(value)!=0) {
      TRAP_ERR_EXIT(-1, "nan encountered at n=%d\n",m);
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int vector_checkNaN(double *A, int zdim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m;
  for(m = 0; m < zdim; m++) {
    double value=A[m];
    if(isnan(value)!=0) {
      TRAP_ERR_EXIT(-1, "nan encountered at n=%d\n",m);
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int vector_checkNaN(float *A, int zdim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m;
  for(m = 0; m < zdim; m++) {
    float value=A[m];
    if(isnan(value)!=0) {
      TRAP_ERR_EXIT(-1, "nan encountered at n=%d\n",m);
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T, typename V> int matrix_operation_CSR_block(T& A, V  *x, V  *y, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  compute y+=Ax
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */ 
{
  int i,k,l,n,nn;
  int row2D, col2D;
  int rowsize, colsize;

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

  In layered ocean equation, some 3D matrices have a 2D ordering pointing
  to structured blocks:
  
  A: 3D operator, with 2D ordering (mxn) with m block lines (m*rowsize)
  
  x: vector with n block lines (n*colsize)
  
  y: vector with m block lines (m*rowsize)
  

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  rowsize=A.ordering->rowsize;
  colsize=A.ordering->colsize;
  
  if(init==1) {
    for(n = 0; n < A.neq*rowsize; n++) {
      y[n]=0;
      }
    }

  for(row2D = 0; row2D < A.neq; row2D++) {
    for(i=0; i < A.ordering->cardinal[row2D]; i++) {
      col2D=A.ordering->incidence[A.ordering->pointer[row2D]+i];
      for(k = 0; k < rowsize; k++) {
        for(l = 0; l < colsize; l++) {
//           nn=indexLGP0xLGP1_2D3D(A.ordering->pointer, rowsize, colsize, row2D, i, k, l);
          nn=index2D3D(A.ordering->pointer, A.ordering->cardinal, rowsize, colsize, row2D, i, k, l);
          y[row2D*rowsize+k]+=A.packed[nn]*x[col2D*colsize+l];
          }
        }
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T, typename V> int matrix_operation_COO(T& A, V  *x, V  *y, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  compute y+=Ax
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */ 
{
//   int i,j,n;
//   vptr_t *row;
//   int *col;
//   int base;
// 
//   if(init==1) {
//     for(n = 0; n < A.neq; n++) {
//       y[n]=0;
//       }
//     }
//   
//   row=A.ordering->pointer;
//   col=A.ordering->incidence;
//   
//   base=A.t.base;
//   
//   for(size_t m = 0; m < A.t.nnz; m++) { // GMRES IMPLEMENTATION, TO BE REWORKED
//     i=row[m]-base;
//     j=col[m]-base;
//     y[i]+=A.packed[m]*x[j];
//     }
    
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T, typename V> int matrix_operation_CSR(T& A, V  *x, V  *y, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  compute y+=Ax
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */ 
{
  int n;
  
#pragma omp parallel
  {
  size_t size,pos;
  int i,pni,m,n;
  
#pragma omp for
  for(n = 0; n < A.neq; n++) {
    
    if(init==1)
      y[n]=0;
    
    size=A.ordering->cardinal[n],
    pos=A.ordering->pointer[n];
    
    for(i=0; i < size ; i++) {
      m=A.ordering->incidence[pos];
      y[n]+=A.packed[pos]*x[m];
      pos++;
      }
    }
  
  }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T, typename V> int matrix_operation_CSC(T& A, V  *x, V  *y, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  compute y+=Ax
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */ 
{
  int n;
  int offset;
  
//   const int nprocs __attribute__((unused)) =gOPENMP_nCPUs;

  if(init==1) {
    for(n = 0; n < A.neq; n++) {
      y[n]=0;
      }
    }
/*------------------------------------------------------------------------------
  if solver has spoiled our nice C numbering*/
  if(A.ordering->incidence[0]!=0) A.ordering->numbering=1;
  switch(A.ordering->numbering) {
    case 0:
      offset=0;
      break;
    case 1:
      offset=1;
      break;
    default:
      TRAP_ERR_EXIT(-1, "ordering has no valid numbering");
      break;
    }
/*------------------------------------------------------------------------------
  following loop CAN NOT be parallelized !!! */
// #pragma omp parallel for if(nprocs>1)
  for(n = 0; n < A.neq; n++) { /** column */
    for(int i=0; i < A.ordering->cardinal[n]; i++) {
      size_t pos=A.ordering->pointer[n]-offset+i;
      int m=A.ordering->incidence[pos]-offset; /** row */
      y[m]+=A.packed[pos]*x[n];
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T, typename V> int matrix_operation_template(T& A, V  *x, V  *y, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  compute A x +=y
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */ 
{
  
  if(A.ordering==0) {
    TRAP_ERR_EXIT(-1, "matrix has no valid ordering");
    }
    
  int status,blocksize=A.ordering->_blocksize();

  if(A.neq<=0) {
    TRAP_ERR_EXIT(-1, "matrix has no valid neq (please check matrix initialization)");
    }

  switch(A.ordering->type) {
    
    case CSR:
      switch(blocksize) {
        case 1:
          status=matrix_operation_CSR(A, x, y, init);
          break;
        default:
          status=matrix_operation_CSR_block(A, x, y, init);
          break;
        }
      break;
      
    case CSC:
      status=matrix_operation_CSC(A, x, y, init);
      break;
      
    case COO:
      status=matrix_operation_COO(A, x, y, init);
      break;
      
    default:
      TRAP_ERR_EXIT(-1, "unknown ordering or ordering not implemented yet");
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_operation(hypermatrix_t& A, complex<double>  *x, complex<double>  *y, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  compute A x =y
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */ 
{
  int status;

  status=matrix_operation_template( A, x, y, init);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_operation(hyperzmatrix_t& A, complex<double>  *x, complex<double>  *y, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  compute A x =y
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */ 
{
  int status;

  status=matrix_operation_template( A, x, y, init);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_operation(hypermatrix_t& A, double  *x, double  *y, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  compute A x =y
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */ 
{
  int status;

  status=matrix_operation_template( A, x, y, init);
  return(status);
}



#endif
