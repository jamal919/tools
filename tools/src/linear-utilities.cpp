
/*******************************************************************************

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

*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <time.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "fe.h"

#include "parallel.h"

#include "periodic.h"

#include "zapper.h" //for zaparr()

int printOK;
double ITER=0.0;
int nbsys=0;

int  LinearSystem_setprint() {
  printOK = 1;
}
int  LinearSystem_unsetprint() {
  printOK = 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int free_ordering(ordering_t **ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(*ordering != NULL) {
    if((*ordering)->pointer!=NULL)   delete[] (*ordering)->pointer;
    (*ordering)->pointer=NULL;
    if((*ordering)->incidence!=NULL) delete[] (*ordering)->incidence;
    (*ordering)->incidence=NULL;
    if((*ordering)->cardinal!=NULL)  delete[] (*ordering)->cardinal;
    (*ordering)->cardinal=NULL;
    }
  *ordering = NULL;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_check_consistency(int *incidence,int *pointer,int *cardinal,int nnodes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n;
  int duplicated,error;
  int cmax;

  cmax=-1;
/**----------------------------------------------------------------------
  check matrix duplicated connections*/
  for(n = 0; n < nnodes; n++) {
/**----------------------------------------------------------------------
    node without connections*/
    if(cardinal[n] == 0) {
      printf("anomaly, null cardinal for node=%d \n",n);
      }
    updatemax(&cmax,cardinal[n]);
    for(k = 0; k < cardinal[n]-1; k++) {
      m=incidence[pointer[n] + k];
/**----------------------------------------------------------------------
      uninitialized connections*/
      if(m == -1) {
        printf("anomaly, uninitialized node=%d, index=%d, connected=%d\n",n,k,m);
        }
/**----------------------------------------------------------------------
      duplicated connections*/
      duplicated=vpos(m,&(incidence[pointer[n] + k +1]),cardinal[n]-k-1);
      if(duplicated!=-1) {
        printf("anomaly, duplicated node=%d, index=%d, connected=%d\n",n,k,m);
        }
      }
    }
  error=0;
  
#ifdef VERBOSE
  printf("matrix consistency check, nnodes=%7d, size=%8d, max cardinal=%3d\n",nnodes,pointer[nnodes],cmax);
#endif
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_print_ordering(const char *filename, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,row;
  FILE *out;

/**----------------------------------------------------------------------------
  MPI unsafe */
  out=fopen(filename,"w");
  
  for(row = 0; row < ordering->nrows; row++) {
    fprintf(out, "%6d :",row);
    for(k = 0; k < ordering->cardinal[row]; k++) {
      m = ordering->incidence[ordering->pointer[row]+k];
      fprintf(out, " %d",m);
      }
    fprintf(out, "\n");
    }
    
  fclose(out);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_print(const char *filename, hypermatrix_t & M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,row;
  FILE *out;
  ordering_t *ordering=M.ordering;
  double value;
  
/**----------------------------------------------------------------------------
  MPI unsafe */
  out=fopen(filename,"w");
  
  for(row = 0; row < ordering->nrows; row++) {
    fprintf(out, "%6d :",row);
    for(k = 0; k < ordering->cardinal[row]; k++) {
      m = ordering->incidence[ordering->pointer[row]+k];
//       fprintf(out, " %d",m);
//       }
//     fprintf(out, "\n");
//     for(k = 0; k < ordering->cardinal[row]; k++) {
      value = M.packed[ordering->pointer[row]+k];
      fprintf(out, "%6d %9.0lf",m, value);
      }
    fprintf(out, "\n");
    }
    
  fclose(out);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_check_consistency(ordering_t ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=matrix_check_consistency(ordering.incidence,ordering.pointer,ordering.cardinal,ordering.nrows);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_check_symmetry(int *incidence,int *pointer,int *cardinal,int nnodes, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n;
  int found,error;
  error=0;

/**----------------------------------------------------------------------
  check matrix duplicated connections*/
  for(n = 0; n < nnodes; n++) {
    for(k = 0; k < cardinal[n]-1; k++) {
      m=incidence[pointer[n] + k];
/**----------------------------------------------------------------------
      non-symetric connections*/
      found=vpos(n,&(incidence[pointer[m]]),cardinal[m]);
      if(found==-1) {
        if(debug) printf("anomaly, non-symetric node=%d, index=%d, connected=%d\n",n,k,m);
        error++;
        }
      }
    }

 return(error);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_check_symmetry(ordering_t ordering, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=matrix_check_symmetry(ordering.incidence,ordering.pointer,ordering.cardinal,ordering.nrows, debug);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_check_use(hyperzmatrix_t & A)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j,n,col,count;
  
  count=0;
  for(n = 0; n < A.neq(); n++) {
    for(j = A.ordering->pointer[n]; j < A.ordering->pointer[n + 1]; j++) {
      col=A.ordering->incidence[j];
      if(A.packed[j]==0.) {
        count++;
        }
      }
    }
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_reorder_01(ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,m,n,row;

  for(row = 0; row < ordering->nrows; row++) {
  iterate:
    for(k = 0; k < ordering->cardinal[row]; k++) {
      m = ordering->incidence[ordering->pointer[row]+k];
      for(l = ordering->cardinal[row] - 1; l > k; l--) {
        n = ordering->incidence[ordering->pointer[row]+l];
        if(n < m) {
          ordering->incidence[ordering->pointer[row]+k] = n;
          ordering->incidence[ordering->pointer[row]+l] = m;
          if(m == -1) {
            printf("anomaly\n");
            }
          goto iterate;
          }
        }
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_reorder(ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const int nprocs=omp_get_num_procs();
#ifdef VERBOSE
  printf("matrix_reorder : OPEN-MP optimzation, nprocs=%d\n",nprocs);
#endif

#pragma omp parallel for if(nprocs>1)
  for(size_t row = 0; row < ordering->nrows; row++) {
    size_t start, pos=ordering->pointer[row], nvalues=ordering->cardinal[row];
    start=0;
iterate:
    for(size_t k = start; k < nvalues; k++) {
      size_t m,n;
      m = ordering->incidence[pos+k];
      for(size_t l = nvalues - 1; l > k; l--) {
        n = ordering->incidence[pos+l];
        if(n < m) {
          ordering->incidence[pos+k] = n;
          ordering->incidence[pos+l] = m;
          if(m == -1) {
            printf("anomaly\n");
            }
          start=k;
          goto iterate;
          }
        }
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_hbw(ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,row;
  int delta,hbw;

  hbw=0;
  for(row = 0; row < ordering->nrows; row++) {
    for(k = 0; k < ordering->cardinal[row]; k++) {
      m = ordering->incidence[ordering->pointer[row]+k];
      delta=abs(m-row);
      hbw=max(delta,hbw);
      }
    }
  return(hbw);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_relative_pos(ordering_t *ordering, int row, int target, int StopOnError)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------

  find the index of column "target" in row line

-----------------------------------------------------------------------------*/
{
  int column;
  int offset=ordering->numbering;

  column=vpos(target,&(ordering->incidence[ordering->pointer[row]-offset]),ordering->cardinal[row]);
  if(column==-1) {
    printf("seeking node %d in %d matrix line (offset=%d, type=%d)\n", target, row, offset, ordering->type);
    if(StopOnError==1) check_error(-1, "ill-posed matrix", __LINE__, __FILE__, 1);
    }
  return(column);
}
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_relative_pos2(ordering_t *ordering, int row, int target, int start, int StopOnError)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  find the index of column "target" in row line

-----------------------------------------------------------------------------*/
{
  int column;
  int offset=ordering->numbering;

  column=vpos(target,&(ordering->incidence[ordering->pointer[row]-offset+start]),ordering->cardinal[row]-start);
  if(column==-1) {
    printf("\nseeking node %d in matrix line %d (offset=%d, type=%d)\n", target, row, offset, ordering->type);
    if(StopOnError==1) check_error(-1, "ill-posed matrix", __LINE__, __FILE__, 1);
    }
  return(column+start);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_absolute_pos(ordering_t *ordering, int row, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------

  find the index of column "target" in row line

-----------------------------------------------------------------------------*/
{
  int column;

  column=vpos(target,&(ordering->incidence[ordering->pointer[row]]),ordering->cardinal[row]);
  if(column==-1) {
    TRAP_ERR_EXIT(-1, "ill-posed matrix\n");
    }
  return(column+ordering->pointer[row]);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  size_t matrix_absolute_pos2(ordering_t *ordering, int row, int start, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  find the index of column "target" in row line

-----------------------------------------------------------------------------*/
{
  size_t column;
  int offset=ordering->numbering;

  column=start+vpos(target,&(ordering->incidence[ordering->pointer[row+start]]),ordering->cardinal[row]-start);
  if(column==-1) {
    printf("\nseeking node %d in matrix line %d (offset=%d, type=%d)\n", target, row, offset, ordering->type);
    TRAP_ERR_EXIT(-1, "ill-posed matrix");
    }
  return(column+ordering->pointer[row]);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_diagonal(hyperzmatrix_t *A)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j,n,col,row,status=0;
  int *diag;
  complex<double> sum;

  diag=new int[A->neq()];

  for(n = 0; n < A->neq(); n++) {
    sum=0;
    for(j = A->ordering->pointer[n]; j < A->ordering->pointer[n + 1]; j++) {
      col=A->ordering->incidence[j];
      row=n;
      if (col==row) {
//        printf("%d %lf\n",n,abs(A->packed[j]));
        diag[n]=j;
        }
      sum+=A->packed[j];
      }
    printf("%d %lf\n",n,abs(sum));
//    A->packed[diag[n]]-=sum;
    for(j = A->ordering->pointer[n]; j < A->ordering->pointer[n + 1]; j++) {
      col=A->ordering->incidence[j];
      row=n;
      if (col!=row) {
        A->packed[j]/=A->packed[diag[n]];
        }
      }
    A->packed[diag[n]]=1.0;
    }


  delete[] diag;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int matrix_IsDiagonal_template(T & A)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,col;
//  int *diag;
  complex<double> sum;

//  diag=new int[A.neq()];
  if(A.ordering->pointer[A.neq()]!=A.neq()) {
    A.IsDiagonal=false;
    return(0);
    }
    
  for(n = 0; n < A.neq(); n++) {
    if(A.ordering->cardinal[n]!=1) {
      A.IsDiagonal=false;
      return(0);
      }
    col=A.ordering->incidence[A.ordering->pointer[n]];
    if(col!=n) {
      A.IsDiagonal=false;
      return(0);
      }
    }

  A.IsDiagonal=true;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_IsDiagonal(hypermatrix_t & A)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  
  status=matrix_IsDiagonal_template(A);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_IsDiagonal(hyperzmatrix_t & A)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  
  status=matrix_IsDiagonal_template(A);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_check_diagonalNaN(hypermatrix_t A)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m,zdim=A.ordering->nrows;
  for(m = 0; m < A.ordering->pointer[zdim]; m++) {
    if(isnan(A.packed[m])!=0) {
      printf("nan encountered at n=%d\n",m);
      TRAP_ERR_EXIT(-1, "nan encountered\n");
      }
    if(A.packed[m]==0) {
      printf("zero encountered at n=%d\n",m);
      TRAP_ERR_EXIT(-1, "nan encountered\n");
      }
   }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_check_diagonalNaN(hyperzmatrix_t A)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m,zdim=A.ordering->nrows;
  int n,col;
  
  for(n = 0; n < A.neq(); n++) {
    col=matrix_relative_pos(A.ordering,n,n,1);
    if(abs(A.packed[A.ordering->pointer[n]+col])==0) {
      printf("zero encountered at n=%d\n",n);
      TRAP_ERR_EXIT(-1, "diagonal zero encountered\n");
      }
   }
  for(m = 0; m < A.ordering->pointer[zdim]; m++) {
    if(isnan(A.packed[m].real())!=0) {
      printf("nan encountered at n=%d\n",m);
      TRAP_ERR_EXIT(-1, "nan encountered\n");
      }
    if(isnan(A.packed[m].imag())!=0) {
      printf("nan encountered at n=%d\n",m);
      TRAP_ERR_EXIT(-1, "nan encountered\n");
      }
//     if(abs(A.packed[m])==0) {
//       printf("zero encountered at n=%d\n",m);
//       check_error(-1, "nan encountered", __LINE__, __FILE__, 1);
//       }
   }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_check_NaN(hypermatrix_t & A)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j,m,zdim;

  if(A.ordering==0) {
    TRAP_ERR_EXIT(-1, "no valid ordering\n");
    }

  zdim=A.ordering->nrows;

  for(m = 0; m < zdim; m++) {
    for(j = A.ordering->pointer[m]; j < A.ordering->pointer[m + 1]; j++) {
      if(isnan(A.packed[j])!=0) {
        printf("nan encountered at n=%d col=%d\n",m,j-A.ordering->pointer[m]);
        TRAP_ERR_EXIT(-1, "nan encountered\n");
        }
      }
   }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int vector_checkNaN(double *A, int zdim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m;
  for(m = 0; m < zdim; m++) {
    double value=A[m];
    if(isnan(value)!=0) {
      printf("nan encountered at n=%d\n",m);
      TRAP_ERR_EXIT(-1, "nan encountered\n");
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int vector_checkNaN(complex<double> *A, int zdim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m;
  for(m = 0; m < zdim; m++) {
    double value=real(A[m]);
    if(isnan(value)!=0) {
      printf("nan encountered at n=%d\n",m);
      TRAP_ERR_EXIT(-1, "nan encountered\n");
      }
    value=imag(A[m]);
    if(isnan(value)!=0) {
      printf("nan encountered at n=%d\n",m);
      TRAP_ERR_EXIT(-1, "nan encountered\n");
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_check_Zero(hypermatrix_t A)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j,m,zdim,count,global;

  if(A.ordering==0) {
    TRAP_ERR_EXIT(-1, "no valid ordering\n");
    }

  zdim=A.ordering->nrows;
  
  global=0;
  for(m = 0; m < zdim; m++) {
    count=0;
    for(j = A.ordering->pointer[m]; j < A.ordering->pointer[m + 1]; j++) {
//      if(A.packed[j]!=0) {
      if(fabs(A.packed[j]) < 1.e-10) {
        count++;
        global++;
        }
      }
    if(count==A.ordering->cardinal[m]) {
      printf("zeroed line at n=%d \n",m);
      TRAP_ERR_EXIT(-1, "singular line encountered\n");
      }
    for(j = A.ordering->pointer[m]; j < A.ordering->pointer[m + 1]; j++) {
//      if(A.packed[j]==0) {
      if(fabs(A.packed[j]) < 1.e-10) {
//        printf("zero encountered at n=%d col=%d\n",m,j-A.ordering->pointer[m]);
//        check_error(-1, "nan encountered", __LINE__, __FILE__, 1);
        }
      }
   }
 return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//  template<typename T> int matrix_CSR2CSC_template(T *A)
  int matrix_CSR2CSC(hyperzmatrix_t *A)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

    This works for geometrically symetric matrices (at ordering level)

    Principle: keep ordering unchanged, transpose matrix (not conjugate)

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

{
  int i,j,m,status=0;
  int nvalues,check,col,row;
//  typeof(A->packed) *tmp;
  complex<double> *tmp;

  if(A->ordering->type!=ORDER_CSR) {
    return(-1);
    }
  if(A->ordering->type==ORDER_CSC) return(0);

  nvalues=A->ordering->pointer[A->neq()];
//  tmp=new typeof(A->packed)[nvalues];
  tmp=new complex<double> [nvalues];
  for(row = 0; row < A->neq(); row++) {
    for(j = A->ordering->pointer[row]; j < A->ordering->pointer[row + 1]; j++) {
      col=A->ordering->incidence[j];
      i=matrix_relative_pos(A->ordering,col,row,1);
      check=matrix_relative_pos(A->ordering,row,col,1);
      m=A->ordering->pointer[col]+i;
/**----------------------------------------------------------------------------
      tmp[m] (=tmp[col,row]) =A[j] (=A[row,col]) */
      tmp[m]=A->packed[j];
      }
    }

  delete[] A->packed;

  A->packed=tmp;
  A->ordering->type=ORDER_CSC;

//  A->t.x = A->packed;
  return(status);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int matrix_CSR2CSC(hyperzmatrix_t *A)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   
//   status=matrix_CSR2CSC_template(A);
//   return(status);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int matrix_CSR2CSC(hypermatrix_t *A)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   
//   status=matrix_CSR2CSC_template(A);
//   return(status);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_CSR2CSC(hypermatrix_t *A)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

    This works for geometrically symetric matrices
    
    Actually it is equivalent to a transpose operation


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int i,j,m,status=0;
  int nvalues,check,col,row;
  double *tmp;

  if(A->ordering->type!=ORDER_CSR) {
    return(-1);
    }
  if(A->ordering->type==ORDER_CSC) return(0);

  nvalues=A->ordering->pointer[A->neq()];
  tmp=new double[nvalues];
  for(row = 0; row < A->neq(); row++) {
    for(j = A->ordering->pointer[row]; j < A->ordering->pointer[row + 1]; j++) {
      col=A->ordering->incidence[j];
      i=matrix_relative_pos(A->ordering,col,row,1);
      check=matrix_relative_pos(A->ordering,row,col,1);
      m=A->ordering->pointer[col]+i;
      tmp[m]=A->packed[j];
      }
    }

  delete[] A->packed;

  A->packed=tmp;
  A->ordering->type=ORDER_CSC;
  
  A->t.x = A->packed;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//  template <typename T> int matrix_CSR2CSC_template(hypermatrix_t *A, ordering_t *CSCordering)
  int matrix_CSR2CSC(hypermatrix_t *A, ordering_t *CSCordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,m,status=0;
  int nvalues,/*check,*/col,row;
  double *tmp;

  if(A->ordering->type!=ORDER_CSR) {
    return(-1);
    }
  if(A->ordering->type==ORDER_CSC) return(0);

  nvalues=A->ordering->pointer[A->neq()];
  tmp=new double[nvalues];
  for(row = 0; row < A->neq(); row++) {
    for(j = A->ordering->pointer[row]; j < A->ordering->pointer[row + 1]; j++) {
      col=A->ordering->incidence[j];
      i=matrix_relative_pos(CSCordering,col,row,1);
//      check=matrix_relative_pos(A->ordering,row,col);
      m=CSCordering->pointer[col]+i;
      tmp[m]=A->packed[j];
      }
    }

  delete[] A->packed;

  A->packed=tmp;
  A->ordering->type=ORDER_CSC;
  
  A->t.x = A->packed;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_CSR2CSC(hyperzmatrix_t *A, ordering_t *CSCordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,m,status=0;
  int nvalues,/*check,*/col,row;
  complex<double> *tmp;

  if(A->ordering->type!=ORDER_CSR) {
    return(-1);
    }
  if(A->ordering->type==ORDER_CSC) return(0);

  nvalues=A->ordering->pointer[A->neq()];
  tmp=new complex<double>[nvalues];
  for(row = 0; row < A->neq(); row++) {
    for(j = A->ordering->pointer[row]; j < A->ordering->pointer[row + 1]; j++) {
      col=A->ordering->incidence[j];
      i=matrix_relative_pos(CSCordering,col,row,1);
//      check=matrix_relative_pos(A->ordering,row,col);
      m=CSCordering->pointer[col]+i;
      tmp[m]=A->packed[j];
      }
    }

  delete[] A->packed;

  A->packed=tmp;
  A->ordering->type=ORDER_CSC;
  
  A->t.x = A->packed;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_CSR2COO(hypermatrix_t *M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

    This works for geometrically symetric matrices


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status = 0;
  int i, k;
  int nndes;
  int *pointer;
  int *incidence;
  triplet *M1;
  int nrow,nnz,idi;
  int *li,*co;
  double *val;
  
  if(M->t.comptype==COO) return(0);
  
  cout << "convert matrix from CSR to COO" << endl;
  nndes = M->t.nrow;

/**----------------------------------------------------------------------------
  init triplet */
  M->t.comptype=COO;
  M->t.base=1;

  M->t.stype=0;

  if(M->neq()==0) {
    TRAP_ERR_EXIT(-1, "neq field is zero\n");
    }
  M->t.nnz=M->ordering->pointer[M->neq()];
  M->t.nrow=nndes;
  M->t.ncol=nndes;

  M->t.i = M->ordering->pointer;
  M->t.j = M->ordering->incidence;
  M->t.x = M->packed;

  M1 = &(M->t);

  pointer   = M1->i;
  incidence = M1->j;

/**----------------------------------------------------------------------------
  unpack matrix to get nnz*/
  nnz = 0;
  nrow = 0;
  for (i=0; i <nndes; i++) {
    nrow++;
    nnz+=pointer[i + 1]-pointer[i];
    }
  M->t.nrow=nrow;
  M->t.ncol=nrow;
  M->t.nnz=nnz;

//  cout << "nnz=" << nnz << endl;

  li  = new int[nnz];
  co  = new int[nnz];
  val = new double[nnz];

/**----------------------------------------------------------------------------
  MUMPS needs COO one based matrix */
  M->t.comptype=COO;
 
  idi=0;
  for (i=0; i <nndes; i++) {
    for ( k=pointer[i];k < pointer[i+1]; k++) {
      li[idi] = i+OFFSET;
      co[idi] = incidence[k] + OFFSET;
      val[idi] = M->packed[k];
      idi ++;
      }
    }

//   free(M->t.j);
//   free(M->t.i);
//   free(M->t.x);
  delete[] M->t.j;
  delete[] M->t.i;
  delete[] M->t.x;

  M->t.i = li;
  M->t.j = co;
  M->t.x = val;
   
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  ordering_t *ordering_CSR2CSC(ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j,status=0;
  int nvalues,col,row;
  int *incidence, *cardinal, *pointer;
  ordering_t *CSCordering;
  int ncols;
  int last;

  if(ordering->type!=ORDER_CSR) {
    return(0);
    }
  if(ordering->type==ORDER_CSC) return(0);

  ncols=-1;
  for(row = 0; row < ordering->nrows; row++) {
    last=ordering->pointer[row]+ordering->cardinal[row]-1;
    updatemax(&ncols,ordering->incidence[last]);
    }
  ncols++;
  
  nvalues=ordering->pointer[ordering->nrows];
  
  incidence=new int[nvalues];
  cardinal =new int[ncols];
  pointer  =new int[ncols+1];
  for(row = 0; row < ncols; row++) {
    cardinal[row]=0;
    }
  
  for(row = 0; row < ordering->nrows; row++) {
    for(j = ordering->pointer[row]; j < ordering->pointer[row + 1]; j++) {
      col=ordering->incidence[j];
      cardinal[col]++;
      }
    }
  pointer[0]=0;
  for(row = 1; row < ncols+1; row++) {
    pointer[row]=pointer[row-1]+cardinal[row-1];
    }

  for(row = 0; row < ncols; row++) {
    cardinal[row]=0;
    }
  
  for(row = 0; row < ordering->nrows; row++) {
    for(j = ordering->pointer[row]; j < ordering->pointer[row + 1]; j++) {
      col=ordering->incidence[j];
      incidence[pointer[col]+cardinal[col]]=row;
      cardinal[col]++;
      }
    }

//   delete[] ordering->pointer;
//   delete[] ordering->cardinal;
//   delete[] ordering->incidence;
  
  CSCordering=new ordering_t;

  CSCordering->nrows=ncols;
  
  CSCordering->pointer=pointer;
  CSCordering->cardinal=cardinal;
  CSCordering->incidence=incidence;
  
  CSCordering->type=ORDER_CSC;
  
  status=matrix_reorder(CSCordering);
  
  return(CSCordering);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_import(hypermatrix_t A, hypermatrix_t B, hypermatrix_t C)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,m,n;
  int row, col, count;
  int duplicated;
  int *incidence;

  C.ordering->nrows=A.ordering->nrows;

  C.ordering->cardinal =new int[C.ordering->nrows];
  C.ordering->pointer  =new int[C.ordering->nrows+1];
  for(n = 0; n < C.ordering->nrows; n++) {
    C.ordering->cardinal[n]=0;
    }

  incidence=new int[C.ordering->nrows];
//  incidence=new int[B.ordering->nrows];
  for(row = 0; row < A.ordering->nrows; row++) {
/**-----------------------------------------------------------------------
    row is line index in matrix A and matrix C*/
    for(k=0; k < A.ordering->cardinal[row]; k++) {
/**-----------------------------------------------------------------------
      m is column index in matrix A and line index in matrix B*/
      m=A.ordering->incidence[A.ordering->pointer[row]+k];
      for(l=0; l < B.ordering->cardinal[m]; l++) {
/**-----------------------------------------------------------------------
        col is column index in matrix B and matrix C*/
        col=B.ordering->incidence[B.ordering->pointer[m]+l];
        duplicated=vpos(col,&(incidence[0]),C.ordering->cardinal[row]);
        if(duplicated==-1) {
          incidence[C.ordering->cardinal[row]]=col;
          C.ordering->cardinal[row]++;
          }
        }
      }
    }

  count=0;
  for(n = 0; n < C.ordering->nrows; n++) {
    C.ordering->pointer[n]=count;
    count+=C.ordering->cardinal[n];
    C.ordering->cardinal[n]=0;
    }
  C.ordering->pointer[n]=count;

  C.ordering->incidence=new int[count];
  for(n = 0; n < count; n++) {
    C.ordering->incidence[n]=-1;
    }

  for(row = 0; row < A.neq(); row++) {
/**-----------------------------------------------------------------------
    row is line index in matrix A and matrix C*/
    for(k=0; k < A.ordering->cardinal[row]; k++) {
/**-----------------------------------------------------------------------
      m is column index in matrix A and line index in matrix B*/
      m=A.ordering->incidence[A.ordering->pointer[row]+k];
      for(l=0; l < B.ordering->cardinal[m]; l++) {
/**-----------------------------------------------------------------------
        col is column index in matrix B and matrix C*/
        col=B.ordering->incidence[B.ordering->pointer[m]+l];
        duplicated=vpos(col,&(C.ordering->incidence[C.ordering->pointer[row]]),C.ordering->cardinal[row]);
        if(duplicated==-1) {
          C.ordering->incidence[C.ordering->pointer[row] + C.ordering->cardinal[row]]=col;
          C.ordering->cardinal[row]++;
          }
        }
      }
    }

  matrix_reorder(C.ordering);

  int status=matrix_check_consistency(C.ordering->incidence,C.ordering->pointer,C.ordering->cardinal,C.ordering->nrows);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ordering_import(ordering_t *A, ordering_t *B, ordering_t *C)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,m,n;
  int row, col, count, ncolmax;
  int duplicated;
  int *incidence;

  C->nrows=A->nrows;

  C->cardinal =new int[C->nrows];
  C->pointer  =new int[C->nrows+1];
  for(n = 0; n < C->nrows; n++) {
    C->cardinal[n]=0;
    }

  ncolmax=C->nrows;
  
  incidence=new int[ncolmax];
  for(row = 0; row < A->nrows; row++) {
/**-----------------------------------------------------------------------
    row is line index in matrix A and matrix C*/
    for(k=0; k < A->cardinal[row]; k++) {
/**-----------------------------------------------------------------------
      m is column index in matrix A and line index in matrix B*/
      m=A->incidence[A->pointer[row]+k];
      for(l=0; l < B->cardinal[m]; l++) {
/**-----------------------------------------------------------------------
        col is column index in matrix B and matrix C*/
        col=B->incidence[B->pointer[m]+l];
        duplicated=vpos(col,&(incidence[0]),C->cardinal[row]);
        if(duplicated==-1) {
          if(C->cardinal[row]>=ncolmax) {
            TRAP_ERR_EXIT(-1, "under-sized vector\n");
            }
          incidence[C->cardinal[row]]=col;
          C->cardinal[row]++;
          }
        }
      }
    }

  count=0;
  for(n = 0; n < C->nrows; n++) {
    C->pointer[n]=count;
    count+=C->cardinal[n];
    C->cardinal[n]=0;
    }
  C->pointer[n]=count;

  C->incidence=new int[count];
  for(n = 0; n < count; n++) {
    C->incidence[n]=-1;
    }

  for(row = 0; row < A->nrows; row++) {
/**-----------------------------------------------------------------------
    row is line index in matrix A and matrix C*/
    for(k=0; k < A->cardinal[row]; k++) {
/**-----------------------------------------------------------------------
      m is column index in matrix A and line index in matrix B*/
      m=A->incidence[A->pointer[row]+k];
      for(l=0; l < B->cardinal[m]; l++) {
/**-----------------------------------------------------------------------
        col is column index in matrix B and matrix C*/
        col=B->incidence[B->pointer[m]+l];
        duplicated=vpos(col,&(C->incidence[C->pointer[row]]),C->cardinal[row]);
        if(duplicated==-1) {
          C->incidence[C->pointer[row] + C->cardinal[row]]=col;
          C->cardinal[row]++;
          }
        }
      }
    }

  matrix_reorder(C);

  int status=matrix_check_consistency(C->incidence,C->pointer,C->cardinal,C->nrows);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_product(hypermatrix_t A, hypermatrix_t B, hypermatrix_t C, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,m,n;
  int row, col, pos;

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

//   for(row = 0; row < A.neq(); row++) {
//     for(k=0; k < A.ordering->cardinal[row]; k++) {
//       m=A.ordering->incidence[A.ordering->pointer[row]+k];
//       if(A.packed[m]==0) {
// //        check_error(-1, "ill-posed matrix_product", __LINE__, __FILE__, 1);
//         printf("%d %d %d\n",row,k,m);
//         }
//       }
//     }
//   for(row = 0; row < B.neq(); row++) {
//     for(k=0; k < B.ordering->cardinal[row]; k++) {
//       m=B.ordering->incidence[B.ordering->pointer[row]+k];
//       if(B.packed[m]==0) {
// //        check_error(-1, "ill-posed matrix_product", __LINE__, __FILE__, 1);
//         printf("%d %d %d\n",row,k,m);
//         }
//       }
//     }

  if(init==1) {
    for(n = 0; n < C.ordering->pointer[C.neq()]; n++) {
      C.packed[n]=0;
      }
    }

  for(row = 0; row < A.neq(); row++) {
/**-----------------------------------------------------------------------
    row is line index in matrix A and matrix C*/
    for(k=0; k < A.ordering->cardinal[row]; k++) {
/**-----------------------------------------------------------------------
      m is column index in matrix A and line index in matrix B*/
      m=A.ordering->incidence[A.ordering->pointer[row]+k];
      for(l=0; l < B.ordering->cardinal[m]; l++) {
/**-----------------------------------------------------------------------
        col is column index in matrix B and matrix C*/
        col=B.ordering->incidence[B.ordering->pointer[m]+l];
        pos=matrix_relative_pos(C.ordering,row,col,1);
//         for(i=0; i < C.ordering->cardinal[row]; i++) {
//           if(C.ordering->incidence[C.ordering->pointer[row]+i]==col) goto next;
//           }
//         check_error(-1, "ill-posed matrix_product", __LINE__, __FILE__, 1);
// next:
        C.packed[C.ordering->pointer[row]+pos]+=A.packed[A.ordering->pointer[row]+k]*B.packed[B.ordering->pointer[m]+l];
        }
      }
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_product(hypermatrix_t A, hyperzmatrix_t B, hyperzmatrix_t C, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,l,m,n;
  int row, col;
//   char filename[50];
//   FILE *out;

  

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  if(init==1) {
    for(n = 0; n < C.ordering->pointer[C.neq()]; n++) {
      C.packed[n]=0;
      }
    }
 
//   nbsys++;
//   //sprintf(filename,"matrixProductDomestic.%d.txt",nbsys);
//   //sprintf(filename,"matrixProductPastiX.%d.txt",nbsys);
//   sprintf(filename,"matrixProduct.%d.%d.txt",solvid,nbsys);
//   out=fopen(filename,"w");
//   for(row = 0; row < 10; row++) {
//   for(k=0; k < A.ordering->cardinal[row]; k++) {
//     m=A.ordering->incidence[A.ordering->pointer[row]+k];
//     for(l=0; l < B.ordering->cardinal[m]; l++) {
//       col=B.ordering->incidence[B.ordering->pointer[m]+l];
//       fprintf(out,"m=%d, col=%d  \n  ",m,col);
//       for(i=0; i < C.ordering->cardinal[row]; i++) {
// 	fprintf(out,"%d %d | ",C.ordering->incidence[C.ordering->pointer[row]+i], C.ordering->pointer[row] );
//       }
//       fprintf(out,"\n");
//     }
//   }
//   }
//   for(row = A.neq()-10; row < A.neq(); row++) {
//   for(k=0; k < A.ordering->cardinal[row]; k++) {
//     m=A.ordering->incidence[A.ordering->pointer[row]+k];
//     for(l=0; l < B.ordering->cardinal[m]; l++) {
//       col=B.ordering->incidence[B.ordering->pointer[m]+l];
//       fprintf(out,"m=%d, col=%d  \n  ",m,col);
//       for(i=0; i < C.ordering->cardinal[row]; i++) {
// 	fprintf(out,"%d %d | ",C.ordering->incidence[C.ordering->pointer[row]+i], C.ordering->pointer[row] );
//       }
//       fprintf(out,"\n");
//     }
//   }
//   }
//   fclose(out);
  

  for(row = 0; row < A.neq(); row++) {
/**-----------------------------------------------------------------------
    row is line index in matrix A and matrix C*/
    for(k=0; k < A.ordering->cardinal[row]; k++) {
/**-----------------------------------------------------------------------
      m is column index in matrix A and line index in matrix B*/
      m=A.ordering->incidence[A.ordering->pointer[row]+k];
      for(l=0; l < B.ordering->cardinal[m]; l++) {
/**-----------------------------------------------------------------------
        col is column index in matrix B and matrix C*/
        col=B.ordering->incidence[B.ordering->pointer[m]+l];

        for(i=0; i < C.ordering->cardinal[row]; i++) {
// 	printf("m=%d,col=%d,C=%d \n",m,col,C.ordering->incidence[C.ordering->pointer[row]+i]);	
          if(C.ordering->incidence[C.ordering->pointer[row]+i]==col) goto next;
          }
        TRAP_ERR_EXIT(-1, "ill-posed matrix_product\n");
next:
        C.packed[C.ordering->pointer[row]+i]+=A.packed[A.ordering->pointer[row]+k]*B.packed[B.ordering->pointer[m]+l];
        }
      }
    }
//  if(nbsys == 3) TRAP_ERR_EXIT(1,"exiting\n");
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_product_new(hypermatrix_t A, hyperzmatrix_t B, hyperzmatrix_t C, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,m,n;
  int kk,ll;
  int row, col, pos;

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

  same as above, improved coding

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  if(init==1) {
    for(n = 0; n < C.ordering->pointer[C.neq()]; n++) {
      C.packed[n]=0;
      }
    }
 
  for(row = 0; row < A.neq(); row++) {
/**-----------------------------------------------------------------------
    row is line index in matrix A and matrix C*/
    for(k=0; k < A.ordering->cardinal[row]; k++) {
/**-----------------------------------------------------------------------
      m is true (unpacked) column index in matrix A and line index in matrix B*/
      kk=A.ordering->pointer[row]+k;
      m=A.ordering->incidence[kk];
      for(l=0; l < B.ordering->cardinal[m]; l++) {
/**-----------------------------------------------------------------------
        col is column index in matrix B and matrix C*/
        ll=B.ordering->pointer[m]+l;
        col=B.ordering->incidence[B.ordering->pointer[m]+l];
/**-----------------------------------------------------------------------
        pos is absolute position in matrix C of true (row,col) coefficient*/
        pos=C.ordering->pos(row,col);
        if(pos==-1) {
          TRAP_ERR_EXIT(-1, "ill-posed matrix_product\n");
          }
        C.packed[pos]+=A.packed[kk]*B.packed[ll];
        }
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int matrix_product_2Dx3D_template(hypermatrix_t& A, T& B, T& C, int vdimA, int vdimB, int nlayers, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

  In layered ocean equation, some 3D operators are a repetition of 2D
  operators:
  
  A: 2D operator, with 2D ordering (repetitive)
  
  B: 3D operator, with 2D ordering
  
  C: 3D operator, with 3D ordering

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  if(init==1) {
    for(size_t n = 0; n < C.ordering->pointer[C.neq()]; n++) {
      C.packed[n]=0;
      }
    }

  const int nprocs=omp_get_num_procs();
  printf("matrix_product_2Dx3D_template : OPEN-MP optimzation, nprocs=%d\n",nprocs);

#pragma omp parallel for if(nprocs>1)
  for(int layer=0;layer<nlayers;layer++) {
    for(size_t row2D = 0; row2D < A.neq(); row2D++) {
/**-----------------------------------------------------------------------
      true (unpacked) 3D line index in A and C */
      size_t row3D=structured3DIndex(row2D, layer+1, A.neq(), vdimA);
      for(int i=0; i < A.ordering->cardinal[row2D]; i++) {
/*------------------------------------------------------------------------------
        m2D is true (unpacked) 2D column index in matrix A*/
/*------------------------------------------------------------------------------
        m2D is true (unpacked) 2D line index in matrix B*/
        size_t pos2D_A=A.ordering->pointer[row2D]+i;
        size_t m2D=A.ordering->incidence[pos2D_A];
        for(int j=0; j < B.ordering->cardinal[m2D]; j++) {
          size_t pos2D_B=B.ordering->pointer[m2D]+j;
/**------------------------------------------------------------------------
          col2D is 2D column index in matrix B*/
          size_t col2D=B.ordering->incidence[pos2D_B];
          for(int l=0; l < vdimA; l++) {
/**-----------------------------------------------------------------------
            pos is absolute position in matrix B of true (row2D,col2D) coefficient*/
            size_t pos3D_B=indexLGP0xLGP1_2D3D(B.ordering->pointer,vdimB,vdimA,m2D,j,layer,l);
/**-----------------------------------------------------------------------
            must be PerLayers ? ... unsafe ? */
            size_t col3D=structured3DIndex(col2D, l, A.neq(), vdimA);
            size_t pos=C.ordering->pos(row3D,col3D);
            if(pos==-1) {
              TRAP_ERR_EXIT(-1, "ill-posed matrix_product\n");
              }
            C.packed[pos]+=A.packed[pos2D_A]*B.packed[pos3D_B];
//            C.packed[pos]+=A.packed[pos2D_A];
//            C.packed[pos]+=B.packed[pos3D_B];
            }
          }
        }
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_product_2Dx3D(hypermatrix_t& A, hypermatrix_t& B, hypermatrix_t& C, int vdimA, int vdimB, int nlayers, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=matrix_product_2Dx3D_template(A, B, C, vdimA, vdimB, nlayers, init);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_product_2Dx3D(hypermatrix_t& A, hyperzmatrix_t& B, hyperzmatrix_t& C, int vdimA, int vdimB, int nlayers, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=matrix_product_2Dx3D_template(A, B, C, vdimA, vdimB, nlayers, init);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_axpy(hypermatrix_t & A, hypermatrix_t & B, double factor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**-----------------------------------------------------------------------
  compute A = A + factor *B       */
  int status=0;
  int k,m,n;
  int posB, offsetA, offsetB, pos;

  if(A.neq()!=B.neq()) {
    TRAP_ERR_EXIT(-1, "ill-posed matrix_axpy\n");
    }

/**-----------------------------------------------------------------------
  if solver has spoiled our nice C numbering*/
  if(A.ordering->incidence[0]!=0) A.ordering->numbering=1;
  switch(A.ordering->numbering) {
    case 0:
      offsetA=0;
      break;
    case 1:
      offsetA=1;
      break;
    default:
      TRAP_ERR_EXIT(-1, "ordering has no valid numbering\n");
      break;
    }
  if(B.ordering->incidence[0]!=0) B.ordering->numbering=1;
  switch(B.ordering->numbering) {
    case 0:
      offsetB=0;
      break;
    case 1:
      offsetB=1;
      printf("first pointer=%d, first incidence=%d, offsetB=%d\n", B.ordering->pointer[0], B.ordering->incidence[0], offsetB);
      break;
    default:
      TRAP_ERR_EXIT(-1, "ordering has no valid numbering\n");
      break;
    }

  for(n = 0; n < A.neq(); n++) {
    for(k=0; k < B.ordering->cardinal[n]; k++) {
      posB=B.ordering->pointer[n]-offsetB+k;          /// posB in C-numbering as array index operates with C-numbering
      m=B.ordering->incidence[posB]-offsetB+offsetA;  /// m in A matrix numbering (could be F- or C-numbering)
      if(m<0) {
        printf("seeking %dth node %d in %d matrix line (offsetA=%d, offsetB=%d)\n", k, m, n, offsetA, offsetB);
//        if(StopOnError==1) check_error(-1, "ill-posed matrix", __LINE__, __FILE__, 1);
        }
      pos=matrix_relative_pos(A.ordering,n,m,1);
      A.packed[A.ordering->pointer[n]+pos]+=factor*B.packed[posB];
      }
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_axpy(hyperzmatrix_t & A, hypermatrix_t & B, double factor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**-----------------------------------------------------------------------
  compute A = A + factor *B       */
  int i,k,m,n;

  if(A.neq()!=B.neq()) {
    TRAP_ERR_EXIT(-1, "ill-posed matrix_axpy\n");
    }

  for(n = 0; n < A.neq(); n++) {
    for(k=0; k < B.ordering->cardinal[n]; k++) {
      m=B.ordering->incidence[B.ordering->pointer[n]+k];
      for(i=0; i < A.ordering->cardinal[n]; i++) {
        if(A.ordering->incidence[A.ordering->pointer[n]+i]==m) goto next;
        }
      TRAP_ERR_EXIT(-1, "ill-posed matrix_axpy\n");
next:
      A.packed[A.ordering->pointer[n]+i]+=factor*B.packed[B.ordering->pointer[n]+k];
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_axpy_CSR(hyperzmatrix_t A, hypermatrix_t B, complex<double> factor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check : MANDATORY !!!

  Note:

  12/02/2013:

    It does work for CSR only. It will be accidentally ok if CSC on
    symmetric matrices

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
/**-----------------------------------------------------------------------
  compute A = A + factor *B */
  int k,m,n;
  int posB, pos;
  int offsetA,offsetB;

//  printf("#matrix_axpy %d %s : indexing type A=%d B=%d)\n",  __LINE__, __FILE__, A.ordering->type, B.ordering->type);

  if(A.neq()!=B.neq()) {
    TRAP_ERR_EXIT(-1, "ill-posed matrix_axpy\n");
    }

  if(A.ordering->incidence[0]!=0) A.ordering->numbering=1;
  switch(A.ordering->numbering) {
    case 0:
      offsetA=0;
      break;
    case 1:
      offsetA=1;
      break;
    default:
      TRAP_ERR_EXIT(-1, "ordering has no valid numbering\n");
      break;
    }
  
  if(B.ordering->incidence[0]!=0) B.ordering->numbering=1;
  switch(B.ordering->numbering) {
    case 0:
      offsetB=0;
      break;
    case 1:
      offsetB=1;
      printf("first pointer=%d, first incidence=%d, offsetB=%d\n", B.ordering->pointer[0], B.ordering->incidence[0], offsetB);
      break;
    default:
      TRAP_ERR_EXIT(-1, "ordering has no valid numbering\n");
      break;
    }

//   for(n = 0; n < A.neq(); n++) {
//     for(k=0; k < B.ordering->cardinal[n]; k++) {
//       m=B.ordering->incidence[B.ordering->pointer[n]+k];
//       for(i=0; i < A.ordering->cardinal[n]; i++) {
//         if(A.ordering->incidence[A.ordering->pointer[n]+i]==m) goto next;
//         }
//       check_error(-1, "ill-posed matrix_axpy", __LINE__, __FILE__, 1);
// next:
//       A.packed[A.ordering->pointer[n]+i]+=factor*B.packed[B.ordering->pointer[n]+k];
//       }
//     }

  for(n = 0; n < A.neq(); n++) {
    for(k=0; k < B.ordering->cardinal[n]; k++) {
      posB=B.ordering->pointer[n]-offsetB+k;
      m=B.ordering->incidence[posB]-offsetB+offsetA;  ///HERE
      if(m<0) {
        printf("seeking %dth node %d in %d matrix line (offsetA=%d, offsetB=%d posB=%d)\n", k, m, n, offsetA, offsetB, posB);
//        if(StopOnError==1) check_error(-1, "ill-posed matrix", __LINE__, __FILE__, 1);
        }
      pos=matrix_relative_pos(A.ordering,n,m,1);
      if(pos==-1) check_error(-1, "ill-posed matrix_axpy", __LINE__, __FILE__, 1);
      A.packed[A.ordering->pointer[n]+pos]+=factor*B.packed[B.ordering->pointer[n]-offsetB+k];
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename TA, typename TF> int matrix_axpy_template(TA & A, hypermatrix_t B, TF factor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**-----------------------------------------------------------------------
  compute A = A+ factor * B */
  int status,blocksize=A.ordering->_blocksize();
  
  if(A.ordering==0) {
    TRAP_ERR_EXIT(-1, "matrix has no valid ordering\n");
    }

  if(A.neq()<=0) {
    TRAP_ERR_EXIT(-1, "matrix has no valid neq (please check matrix initialization)\n");
    }

/**-----------------------------------------------------------------------
  control on A is not safe enough */
  switch(A.ordering->type) {
    case ORDER_CSR:
      switch(blocksize) {
        case 1:
          status=matrix_axpy_CSR(A, B, factor);
          break;
        default:
//          status=matrix_operation_CSR_block(A, x, y, init);
          break;
        }
      break;
      
    case ORDER_CSC:
      TRAP_ERR_EXIT(-1, "COO ordering not implemented yet\n");
      break;
      
    case ORDER_COO:
      TRAP_ERR_EXIT(-1, "COO ordering not implemented yet\n");
      break;
      
    default:
      TRAP_ERR_EXIT(-1, "unknown ordering or ordering not implemented yet\n");
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_axpy(hyperzmatrix_t & A, hypermatrix_t & B, complex<double> factor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=matrix_axpy_template(A, B, factor);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_axpy_3Dx2D(hyperzmatrix_t A, hypermatrix_t B, complex<double> factor, int shift, int nlayers)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

  In layered ocean equation, some 3D operators are a repetition of 2D
  operators:
  
  A: 3D operator, with 3D ordering
  
  B: 2D operator, with 2D ordering (repetitive)
  

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/**-----------------------------------------------------------------------
  compute A = A + factor *B       */
  int k;
  int row, col, layer;
  int row3D, col3D;
  int first,last;
  int offsetB;
  
  if(A.neq()!=B.neq()*nlayers) {
    TRAP_ERR_EXIT(-1, "ill-posed matrix_axpy\n");
    }

  if(shift>0) {
    first=0;
    last=nlayers-shift;
    }
  else if(shift<0) {
    first=-shift;
    last=nlayers;
    }
  else {
    first=0;
    last=nlayers;
    }
    
/**-----------------------------------------------------------------------
  if solver has spoiled our nice C numbering*/
  if(B.ordering->incidence[0]!=0) {                              // HERE !!!
    printf("\nmatrix_axpy_3Dx2D : warning, solver has probably spoiled our nice C numbering\n\n");
    B.ordering->numbering=1;
    }
  switch(B.ordering->numbering) {
    case 0:
      offsetB=0;
      break;
    case 1:
      offsetB=1;
      break;
    default:
      TRAP_ERR_EXIT(-1, "ordering has no valid numbering\n");
      break;
    }
    
  for(layer=first;layer<last;layer++) {
    for(row = 0; row < B.neq(); row++) {
      row3D=structured3DIndex(row, layer, B.neq(), nlayers);
      for(k=0; k < B.ordering->cardinal[row]; k++) {
        size_t posB=B.ordering->pointer[row]-offsetB+k;             // HERE !!!
        col=B.ordering->incidence[posB]-offsetB;                     // HERE !!!
        col3D=structured3DIndex(col, layer+shift, B.neq(), nlayers);
        size_t posA=A.ordering->pos(row3D,col3D);
        if(posA==-1) {
          posA=A.ordering->pos(row3D,col3D);
          TRAP_ERR_EXIT(-1, "matrix_axpy_3Dx2D: possibly ordering changed by solver (pastix does...)\n");
          }
//        A.packed[pos]+=factor*B.packed[B.ordering->pointer[row]+k];
        A.packed[posA]+=factor*B.packed[posB];
        }
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_operation(hypermatrix_t& A,const double *x, double *y, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**-----------------------------------------------------------------------
  compute A x = y */
  int i,m,n;
  int pos, offset;

  const int nprocs=omp_get_num_procs();

  if(init==1) {
#pragma omp parallel for if(nprocs>1)
    for(n = 0; n < A.neq(); n++) {
      y[n]=0;
      }
    }
    
/**-----------------------------------------------------------------------
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
      TRAP_ERR_EXIT(-1, "ordering has no valid numbering\n");
      break;
    }
    
#pragma omp parallel for private(i,m,pos) if(nprocs>1)
  for(n = 0; n < A.neq(); n++) {
    for(i=0; i < A.ordering->cardinal[n]; i++) {
      pos=A.ordering->pointer[n]-offset+i;
      m=A.ordering->incidence[pos]-offset;
      y[n]+=A.packed[pos]*x[m];
//       if(isnan(y[n])!=0) {
// 	printf("trouble\n");
//         }
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename M, typename V> int matrix_operation_CSR_block(M & A,const V *x, complex<double> *y, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**-----------------------------------------------------------------------
  compute A x =y */
  int i,k,l,n,nn;
  int row2D, col2D;
  int rowsize, colsize;

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

  In layered ocean equation, some 3D matrices have a 2D ordering pointing
  to structured blocks:
  
  A: 3D operator, with 2D ordering (mxn) with m block lines (m*rowsize)
  
  x: vector with n block lines (n*colsize)
  
  y: vector with m block lines (m*rowsize)
  

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  rowsize=A.ordering->rowsize;
  colsize=A.ordering->colsize;
  
  if(init==1) {
    for(n = 0; n < A.neq()*rowsize; n++) {
      y[n]=0;
      }
    }

  for(row2D = 0; row2D < A.neq(); row2D++) {
    for(i=0; i < A.ordering->cardinal[row2D]; i++) {
      col2D=A.ordering->incidence[A.ordering->pointer[row2D]+i];
      for(k = 0; k < rowsize; k++) {
        for(l = 0; l < colsize; l++) {
          nn=indexLGP0xLGP1_2D3D(A.ordering->pointer, rowsize, colsize, row2D, i, k, l);
          y[row2D*rowsize+k]+=A.packed[nn]*x[col2D*colsize+l];
          }
        }
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename M, typename V> int matrix_operation_CSR(M & A,const V *x, complex<double> *y, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**-----------------------------------------------------------------------
  compute A x = y */
  int i,m,n;
  
#define CHKNAN

  if(init==1) {
    for(n = 0; n < A.neq(); n++) {
      y[n]=0;
      }
    }
  for(n = 0; n < A.neq(); n++) {
    for(i=0; i < A.ordering->cardinal[n]; i++) {
      m=A.ordering->incidence[A.ordering->pointer[n]+i];
      y[n]+=A.packed[A.ordering->pointer[n]+i]*x[m];
#ifdef CHKNAN
      if(isnan(y[n].real())==1){
        printf("Nan %d %d\n",n,i);
        }
      if(isnan(y[n].imag())==1){
        printf("Nan %d %d\n",n,i);
        }
#endif
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename M, typename V> int matrix_operation_CSC(M & A,const V *x, complex<double> *y, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**-----------------------------------------------------------------------
  compute A x  =y */
  int i,m,n;

  if(init==1) {
    for(n = 0; n < A.neq(); n++) {
      y[n]=0;
      }
    }
  for(n = 0; n < A.neq(); n++) { /** column */
    for(i=0; i < A.ordering->cardinal[n]; i++) {
      m=A.ordering->incidence[A.ordering->pointer[n]+i]; /** row */
      y[m]+=A.packed[A.ordering->pointer[n]+i]*x[n];
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename M, typename V> int matrix_operation_template(M & A,const V *x, complex<double> *y, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**-----------------------------------------------------------------------
  compute A x =y */
  int status,blocksize=A.ordering->_blocksize();
  
  if(A.ordering==0) {
    TRAP_ERR_EXIT(-1, "matrix has no valid ordering\n");
    }

  if(A.neq()<=0) {
    TRAP_ERR_EXIT(-1, "matrix has no valid neq (please check matrix initialization)\n");
    }

  switch(A.ordering->type) {
    
    case ORDER_CSR:
      switch(blocksize) {
        case 1:
          status=matrix_operation_CSR(A, x, y, init);
          break;
        default:
          status=matrix_operation_CSR_block(A, x, y, init);
          break;
        }
      break;
      
    case ORDER_CSC:
      status=matrix_operation_CSC(A, x, y, init);
      break;
      
    case ORDER_COO:
      TRAP_ERR_EXIT(-1, "COO ordering not implemented yet\n");
      break;
      
    default:
      TRAP_ERR_EXIT(-1, "unknown ordering or ordering not implemented yet\n");
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_operation(hypermatrix_t& A,const complex<double> *x, complex<double> *y, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**-----------------------------------------------------------------------
  compute A x = y */
  int status;

  status=matrix_operation_template( A, x, y, init);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_operation(hyperzmatrix_t& A,const double *x, complex<double> *y, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**-----------------------------------------------------------------------
  compute A x = y */
  int status;

  status=matrix_operation_template( A, x, y, init);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_operation(hyperzmatrix_t& A,const complex<double> *x, complex<double> *y, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**-----------------------------------------------------------------------
  compute A x = y */
  int status;

  status=matrix_operation_template( A, x, y, init);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_operation(hypermatrix_t& A, complex<double>  *x, complex<double>  *y, int init)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**-----------------------------------------------------------------------
  compute A = A + factor *B*/
  int i,j,k,l,m,n;
  int row, col;

  if(init==1) {
    for(n = 0; n < A.neq(); n++) {
      y[n]=0;
      }
    }
  for(n = 0; n < A.neq(); n++) {
    for(i=0; i < A.ordering->cardinal[n]; i++) {
      m=A.ordering->incidence[A.ordering->pointer[n]+i];
      y[n]+=A.packed[A.ordering->pointer[n]+i]*x[m];
      }
    }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_shrink(hyperzmatrix_t & M, const vector<paire_t> & paires, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  remove zero coefficients in matrix to minor size
  
  MPI_Barrier(TUGO_COMM_WORLD) is MPI non-compliant if called for local solving

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  complex<double> *buffer, zero=complex<double>(0.0, 0.0);
  ordering_t *tmp=new ordering_t;
//   T tmp;
  bool spokesman=(gCPU_ID==gCPU_MASTER);
  
  if(verbose==1 and spokesman) printf("%s : matrix shrink\n",__func__);
  
  size_t n,nvalues=M.ordering->size();
  
  size_t nzeroes=M.nzero();
  
//   T.ordering=new ordering_t;
  
  tmp->nrows=M.ordering->nrows;
  tmp->ncols=M.ordering->ncols;
  tmp->type=M.ordering->type;
  
  tmp->cardinal=new int[M.ordering->nrows];
  tmp->pointer =new vptr_t[M.ordering->nrows+1];
  
  tmp->pointer[0]=0;
  for(size_t n=0;n<M.ordering->nrows;n++) {
    tmp->cardinal[n]=0;
    for(size_t pos=M.ordering->pointer[n];pos<M.ordering->pointer[n+1]; pos++) {
      if(M.packed[pos]!=zero) tmp->cardinal[n]++;
      }
    tmp->pointer[n+1]=tmp->pointer[n]+tmp->cardinal[n];
    }
  
  n=M.ordering->nrows;
  tmp->incidence =new int[tmp->pointer[n]];
  
  if(tmp->incidence==0) check_error(-1, "memory allocation failure", __LINE__, __FILE__, 1);
  
  for(size_t n=0;n<M.ordering->nrows;n++) {
    int k=0;
    for(size_t pos=M.ordering->pointer[n];pos<M.ordering->pointer[n+1]; pos++) {
      if(M.packed[pos]!=zero) {
        tmp->incidence[tmp->pointer[n]+k]=M.ordering->incidence[pos];
        k++;
        }
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  handle block multi-variate matrices
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(int k=0;k<M.ordering->hdim.size();k++) {
    int hdim=M.ordering->hdim[k];
    int vdim=M.ordering->vdim[k];
    if(tmp->start.size()==0) tmp->start.push_back(0);
    else {
      int v=tmp->start.size()-1;
      size_t offset=tmp->start[v]+tmp->hdim[v]*tmp->vdim[v];
      tmp->start.push_back(offset);
      }
    tmp->hdim.push_back(hdim);
    tmp->vdim.push_back(vdim);
    }
   
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  handle periodicity
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(paires.size()!=0) {
    vector<paire_t> paires3D;
    size_t offset=0;
  
    for(int k=0;k<M.ordering->hdim.size();k++) {
      int hdim=M.ordering->hdim[k];
      int vdim=M.ordering->vdim[k];
      offset=tmp->start[k];
//       if(tmp->start.size()==0) tmp->start.push_back(0);
//       else {
//         int v=tmp->start.size()-1;
//         size_t offset=tmp->start[v]+tmp->hdim[v]*tmp->vdim[v];
//         tmp->start.push_back(offset);
//         }
//       tmp->hdim.push_back(hdim);
//       tmp->vdim.push_back(vdim);
      for (int l=0; l< paires.size(); l++ ) {
        int n1=paires[l].value[0];
        int n2=paires[l].value[1];
        for (int j=0; j<vdim; j++ ) {
          paire_t *p=new paire_t;
          *p=paire_t(offset+n1*vdim+j,offset+n2*vdim+j);
          int m1=offset+structured3DIndex(n1, j, hdim, vdim);
          int m2=offset+structured3DIndex(n2, j, hdim, vdim);
          paires3D.push_back(*p);
          }
       }
      offset+=hdim*vdim;
      }
   
    status=periodic_ordering(tmp, paires3D, true);

    paires3D.clear();
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  diagnostics
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  double matrix_size;
  if(verbose==1) {
/*------------------------------------------------------------------------------
    MPI non-compliant if called for local solving*/
//     status=P_MPI_Barrier(TUGO_COMM_WORLD);
    size_t size1, size2;
    double sizeMo1, sizeMo2;
    size1=M.ordering->size();
    size2=tmp->size();
    sizeMo1=16*size1/1024./1024.;
    sizeMo2=16*size2/1024./1024.;
    printf("%s, cpu=%d : initial size of matrix %d (%.3lfMo), shrinked size %d (%.3lfMo), reduction ratio=%lf\n", __func__, gCPU_ID, size1, sizeMo1, size2, sizeMo2, sizeMo2/sizeMo1);
    matrix_size=16*tmp->size()/1024./1024.;
//     matrix_size=P_sum(matrix_size);
//     status=P_MPI_Barrier(TUGO_COMM_WORLD);
    }
  
  if(spokesman and verbose==1) {
    double ordering_size=matrix_size/4, system_size=matrix_size+ordering_size,duplication_size=matrix_size+2*ordering_size, total=system_size;
    printf("system memory use : matrix=%.3lfMo, ordering=%.3lfMo, system=%.3lfMo\n", matrix_size, ordering_size, system_size);
    duplication_size=matrix_size+2*ordering_size;
    total+=duplication_size;
    printf("overheads for duplication : matrix=%.3lfMo, ordering=%.3lfMo, system=%.3lfMo, total=%.3lfMo\n", matrix_size, 2*ordering_size, duplication_size, total);
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  set-up new packed array
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  buffer=new complex<double>[tmp->size()];
  
  if(buffer==0) printf("%s, cpu=%d: shrinked matrix allocation failure\n", __FUNCTION__, gCPU_ID);
  
  for(size_t m=0;m<tmp->size();m++) buffer[n]=0;
  
  for(size_t n=0;n<M.ordering->nrows;n++) {
    int k=0;
    size_t pos2=0;
    for(size_t pos=M.ordering->pointer[n];pos<M.ordering->pointer[n+1]; pos++) {
      if(M.packed[pos]!=zero) {
        size_t col=M.ordering->incidence[pos];
        pos2=matrix_relative_pos2(tmp, n, col, pos2, 1);
        buffer[tmp->pointer[n]+pos2]=M.packed[pos];
        k++;
        }
      }
    }
  
  delete[] M.packed;
  M.packed=buffer;
   
  M.ordering->destroy();
  *M.ordering=*tmp;
        
  if(debug) status= matrix_check_consistency(*tmp);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_shrink(hypermatrix_t & M, const vector<paire_t> & paires, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  remove zero coefficients in matrix to minor size
  
  MPI_Barrier(TUGO_COMM_WORLD) is MPI non-compliant if called for local solving

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  double *buffer, zero=0.0;
  ordering_t *tmp=new ordering_t;
//   T tmp;
  bool spokesman=(gCPU_ID==gCPU_MASTER);
  
  if(verbose==1 and spokesman) printf("%s : matrix shrink\n",__func__);
  
  size_t n,nvalues=M.ordering->size();
  
  size_t nzeroes=M.nzero();
  
//   T.ordering=new ordering_t;
  
  tmp->nrows=M.ordering->nrows;
  tmp->ncols=M.ordering->ncols;
  tmp->type=M.ordering->type;
  
  tmp->cardinal=new int[M.ordering->nrows];
  tmp->pointer =new vptr_t[M.ordering->nrows+1];
  
  tmp->pointer[0]=0;
  for(size_t n=0;n<M.ordering->nrows;n++) {
    tmp->cardinal[n]=0;
    for(size_t pos=M.ordering->pointer[n];pos<M.ordering->pointer[n+1]; pos++) {
      if(M.packed[pos]!=zero) tmp->cardinal[n]++;
      }
    tmp->pointer[n+1]=tmp->pointer[n]+tmp->cardinal[n];
    }
  
  n=M.ordering->nrows;
  tmp->incidence =new int[tmp->pointer[n]];
  
  if(tmp->incidence==0) check_error(-1, "memory allocation failure", __LINE__, __FILE__, 1);
  
  for(size_t n=0;n<M.ordering->nrows;n++) {
    int k=0;
    for(size_t pos=M.ordering->pointer[n];pos<M.ordering->pointer[n+1]; pos++) {
      if(M.packed[pos]!=zero) {
        tmp->incidence[tmp->pointer[n]+k]=M.ordering->incidence[pos];
        k++;
        }
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  handle block multi-variate matrices
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(int k=0;k<M.ordering->hdim.size();k++) {
    int hdim=M.ordering->hdim[k];
    int vdim=M.ordering->vdim[k];
    if(tmp->start.size()==0) tmp->start.push_back(0);
    else {
      int v=tmp->start.size()-1;
      size_t offset=tmp->start[v]+tmp->hdim[v]*tmp->vdim[v];
      tmp->start.push_back(offset);
      }
    tmp->hdim.push_back(hdim);
    tmp->vdim.push_back(vdim);
    }
   
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  handle periodicity
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  fragile 2D/3D detect */
  bool is3D=(M.ordering->vdim.size()!=0);
  
  if(paires.size()!=0) {
    vector<paire_t> paires3D;
    size_t offset=0;
  
    for(int k=0;k<M.ordering->hdim.size();k++) {
      int hdim=M.ordering->hdim[k];
      int vdim=M.ordering->vdim[k];
      offset=tmp->start[k];
//       if(tmp->start.size()==0) tmp->start.push_back(0);
//       else {
//         int v=tmp->start.size()-1;
//         size_t offset=tmp->start[v]+tmp->hdim[v]*tmp->vdim[v];
//         tmp->start.push_back(offset);
//         }
//       tmp->hdim.push_back(hdim);
//       tmp->vdim.push_back(vdim);
      for (int l=0; l< paires.size(); l++ ) {
        int n1=paires[l].value[0];
        int n2=paires[l].value[1];
        for (int j=0; j<vdim; j++ ) {
          paire_t *p=new paire_t;
          *p=paire_t(offset+n1*vdim+j,offset+n2*vdim+j);
          int m1=offset+structured3DIndex(n1, j, hdim, vdim);
          int m2=offset+structured3DIndex(n2, j, hdim, vdim);
          paires3D.push_back(*p);
          }
       }
      offset+=hdim*vdim;
      }
   
    if(is3D) {
      status=periodic_ordering(tmp, paires3D, true);
      }
    else {
      status=periodic_ordering(tmp, paires, true);
      }

    paires3D.clear();
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  diagnostics
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  double matrix_size;
  if(verbose==1) {
/*------------------------------------------------------------------------------
    MPI non-compliant if called for local solving*/
//     status=P_MPI_Barrier(TUGO_COMM_WORLD);
    size_t size1, size2;
    double sizeMo1, sizeMo2;
    size1=M.ordering->size();
    size2=tmp->size();
    sizeMo1=16*size1/1024./1024.;
    sizeMo2=16*size2/1024./1024.;
    printf("%s, cpu=%d : initial size of matrix %d (%.3lfMo), shrinked size %d (%.3lfMo ratio=%lf)\n", __func__, gCPU_ID, size1, sizeMo1, size2, sizeMo2, sizeMo2/sizeMo1);
    matrix_size=16*tmp->size()/1024./1024.;
//     matrix_size=P_sum(matrix_size);
//     status=P_MPI_Barrier(TUGO_COMM_WORLD);
    }
  
  if(spokesman and verbose==1) {
    double ordering_size=matrix_size/4, system_size=matrix_size+ordering_size,duplication_size=matrix_size+2*ordering_size, total=system_size;
    printf("system memory use : matrix=%.3lfMo, ordering=%.3lfMo, system=%.3lfMo\n", matrix_size, ordering_size, system_size);
    duplication_size=matrix_size+2*ordering_size;
    total+=duplication_size;
    printf("overheads for duplication : matrix=%.3lfMo, ordering=%.3lfMo, system=%.3lfMo, total=%.3lfMo\n", matrix_size, 2*ordering_size, duplication_size, total);
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  set-up new packed array
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  buffer=new double[tmp->size()];
  
  if(buffer==0) printf("%s, cpu=%d: shrinked matrix allocation failure\n", __FUNCTION__, gCPU_ID);
  
  for(size_t m=0;m<tmp->size();m++) buffer[m]=0;
  
  for(size_t n=0;n<M.ordering->nrows;n++) {
    int k=0;
    size_t pos2=0;
    for(size_t pos=M.ordering->pointer[n];pos<M.ordering->pointer[n+1]; pos++) {
      if(M.packed[pos]!=zero) {
        size_t col=M.ordering->incidence[pos];
        pos2=matrix_relative_pos2(tmp, n, col, pos2, 1);
        buffer[tmp->pointer[n]+pos2]=M.packed[pos];
        k++;
        }
      }
    }
  
  delete[] M.packed;
  M.packed=buffer;
   
  M.ordering->destroy();
  *M.ordering=*tmp;
        
  if(debug) status= matrix_check_consistency(*tmp);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_shrink(hyperzmatrix_t & M, paire_t *paires, int npaires, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  vector<paire_t> v;
  
  for (int l=0; l< npaires; l++ ) {
    int n1=paires[l].value[0];
    int n2=paires[l].value[1];
    paire_t *p=new paire_t;
    *p=paire_t(n1,n2);
    v.push_back(*p);
    }
  
  status=matrix_shrink(M, v, verbose, debug);
  
  v.clear();
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ordering_expand3D(ordering_t & ordering2D, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Expand a 2D ordering into a 3D ordering:
  
  * multiply number of rows by rowsize
  * connect each (expanded) 3D row with full (expanded) 3D column

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i, m, n, j, k, l, status;
  int row, col, nrows, colsize;
  int *tmp=0;
  int *incidence=0, *cardinal=0;
  vptr_t *pointer=0;
  bool debug=false;

  size_t count,size;
  
  if(ordering2D.ncols==-1) {
    TRAP_ERR_EXIT(-1,"matrix ncols not initialised\n");
    }
  
  range_t<int> r=ordering2D.stat();
  
  if(debug) {
    printf("%s : nrows=%d rowsize=%d rmax=%d\n",__FUNCTION__,ordering2D.nrows,ordering2D.rowsize,r.max);
    printf("%s : ncols=%d colsize=%d\n",__FUNCTION__,ordering2D.ncols,ordering2D.colsize,r);
    }
    
  colsize = ordering2D.colsize;
  nrows = ordering2D.nrows*ordering2D.rowsize;

/*------------------------------------------------------------------------------
  allocate arrays*/
  cardinal = new int[nrows];
  pointer  = new vptr_t[nrows + 1];
  
  size=nrows * colsize * r.max;
  tmp = new int[size];

/*------------------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < size; n++)
    tmp[n] = -1;
  
  if(debug) {
    printf("%s : allocate temporary array of %d\n",__FUNCTION__,size);
    }

/*------------------------------------------------------------------------------
  */
  count = 0;
  for(int m2D = 0; m2D < ordering2D.nrows; m2D++) {
    for(k = 0; k <  ordering2D.rowsize; k++) {
      m=m2D*ordering2D.rowsize+k;
      pointer [m] = count;
      cardinal[m] = ordering2D.cardinal[m2D]*colsize;
      count+=cardinal[m];
      int pos=0;
      for(int j = 0; j < ordering2D.cardinal[m2D]; j++) {
        for(int l = 0; l < ordering2D.colsize; l++) {
          size_t nn,nz2D,nz3D;
          nz2D=ordering2D.incidence[ordering2D.pointer[m2D]+j];
          nz3D=structured3DIndex(nz2D, l, ordering2D.ncols, ordering2D.colsize);
          tmp[m*colsize*r.max + pos] = nz3D;
          pos++;
          }
        }
      }
    }
    
  pointer[nrows] = count;
  
  if(debug) {
    printf("%s : array actual size %d\n",__FUNCTION__,count);
    }

  incidence = new int[count];

  for(n = 0; n < nrows; n++) {
    for(k = 0; k < cardinal[n]; k++) {
       /*incidence= row index */
      incidence[pointer[n] + k] = tmp[n * colsize * r.max + k];
      if(incidence[pointer[n] + k] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }

  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows     = nrows;
  ordering->ncols     = ordering2D.ncols*ordering2D.colsize;

//   status= matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,nrows);
  
  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ordering_expand3D(ordering_t & ordering2D, ordering_t *ordering, expand_t & connexions)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Expand a 2D ordering into a 3D ordering:
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i, m, n, j, k, l, status;
  int row, col, nrows, colsize;
  int *tmp=0;
  int *incidence=0, *cardinal=0;
  vptr_t *pointer=0;
  bool debug=false;

  size_t count,size;
  
  if(ordering2D.ncols==-1) {
    TRAP_ERR_EXIT(-1,"matrix ncols not initialised\n");
    }
  
  range_t<int> r=ordering2D.stat();
  
  if(debug) {
    printf("%s : nrows=%d rowsize=%d rmax=%d\n",__FUNCTION__,ordering2D.nrows,ordering2D.rowsize,r.max);
    printf("%s : ncols=%d colsize=%d\n",__FUNCTION__,ordering2D.ncols,ordering2D.colsize,r);
    }
    
  colsize = connexions.maxsize();
  nrows   = ordering2D.nrows*ordering2D.rowsize;

/*------------------------------------------------------------------------------
  allocate arrays*/
  cardinal = new int[nrows];
  pointer  = new vptr_t[nrows + 1];
  
  size=nrows * colsize * r.max;
  tmp = new int[size];

/*------------------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < size; n++)
    tmp[n] = -1;
  
  if(debug) {
    printf("%s : allocate temporary array of %d\n",__FUNCTION__,size);
    }

/*------------------------------------------------------------------------------
  */
  count = 0;
  for(int m2D = 0; m2D < ordering2D.nrows; m2D++) {
    for(k = 0; k <  ordering2D.rowsize; k++) {
      m=m2D*ordering2D.rowsize+k;
      pointer [m] = count;
      cardinal[m] = ordering2D.cardinal[m2D]*connexions.rows[k].incidence.size();
      count+=cardinal[m];
      int pos=0;
      for(int j = 0; j < ordering2D.cardinal[m2D]; j++) {
        for(int ll = 0; ll < connexions.rows[k].incidence.size(); ll++) {
          l=connexions.rows[k].incidence[ll];
          size_t nn,nz2D,nz3D;
          nz2D=ordering2D.incidence[ordering2D.pointer[m2D]+j];
          nz3D=structured3DIndex(nz2D, l, ordering2D.ncols, ordering2D.colsize);
          tmp[m*colsize*r.max + pos] = nz3D;
          pos++;
          }
        }
      }
    }
    
  pointer[nrows] = count;
  
  if(debug) {
    printf("%s : array actual size %d\n",__FUNCTION__,count);
    }

  incidence = new int[count];

  for(n = 0; n < nrows; n++) {
    for(k = 0; k < cardinal[n]; k++) {
       /*incidence= row index */
      incidence[pointer[n] + k] = tmp[n * colsize * r.max + k];
      if(incidence[pointer[n] + k] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }

  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows     = nrows;
  ordering->ncols     = ordering2D.ncols*ordering2D.colsize;

//   status= matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,nrows);
  
  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ordering_expand3D_01(ordering_t & ordering2D, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Expand a 2D ordering into a 3D ordering:
  
  * multiply number of rows by rowsize
  * connect each (expanded) 3D row with full (expanded) 3D column
  
  WARNING : vertical scheme and discretisation dependent

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i, m, n, j, k, l, status;
  int row, col, nrows, ncols;
  int *tmp;
  int *incidence=0, *cardinal=0;
  vptr_t *pointer=0;

  size_t count,size;
  
  range_t<int> r=ordering2D.stat();

  ncols = ordering2D.colsize;
  nrows = ordering2D.nrows*ordering2D.rowsize;

/*------------------------------------------------------------------------------
  allocate arrays*/
  cardinal = new int[nrows];
  pointer  = new vptr_t[nrows + 1];
  
  size=nrows * r.max;
  tmp = new int[size];

/*------------------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < size; n++)
    tmp[n] = -1;

/*------------------------------------------------------------------------------
  */
  count = 0;
  for(int m2D = 0; m2D < ordering2D.nrows; m2D++) {
    for(k = 0; k <  ordering2D.rowsize; k++) {
      m=m2D*ordering2D.rowsize+k;
      pointer [m] = count;
      cardinal[m] = ordering2D.cardinal[m2D];
      count+=cardinal[m];
      int pos=0;
      for(int j = 0; j < ordering2D.cardinal[m2D]; j++) {
        int l=k/2;
        size_t nz2D,nz3D;
        nz2D=ordering2D.incidence[ordering2D.pointer[m2D]+j];
        nz3D=structured3DIndex(nz2D, l, ordering2D.ncols, ordering2D.colsize);
        tmp[m*r.max + pos] = nz3D;
        pos++;
        }
      }
    }
    
  pointer[nrows] = count;

  incidence = new int[count];

  for(n = 0; n < nrows; n++) {
    for(k = 0; k < cardinal[n]; k++) {
       /*incidence= row index */
      incidence[pointer[n] + k] = tmp[n * r.max + k];
      if(incidence[pointer[n] + k] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }

  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows     = nrows;
  ordering->ncols     = ordering2D.ncols*ordering2D.colsize;

  status= matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,nrows);
  
  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ordering_expand3D_02(ordering_t & ordering2D, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Expand a 2D ordering into a 3D ordering:
  
  * multiply number of rows by rowsize
  * connect each (expanded) 3D row with full (expanded) 3D column
  
  WARNING : vertical scheme and discretisation dependent

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i, m, n, j, k, l, status;
  int row, col, nrows, ncols;
  int *tmp=0;
  int *incidence=0, *cardinal=0;
  vptr_t *pointer=0;

  size_t count,size;
  
  range_t<int> r=ordering2D.stat();

  ncols = ordering2D.colsize;
  nrows = ordering2D.nrows*ordering2D.rowsize;

/*------------------------------------------------------------------------------
  allocate arrays*/
  cardinal = new int[nrows];
  pointer  = new vptr_t[nrows + 1];
  
  size=nrows * r.max;
  tmp = new int[size];

/*------------------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < size; n++)
    tmp[n] = -1;

/*------------------------------------------------------------------------------
  construct tmp */
  count = 0;
  for(int m2D = 0; m2D < ordering2D.nrows; m2D++) {
    for(k = 0; k <  ordering2D.rowsize; k++) {
      m=m2D*ordering2D.rowsize+k;  // 3D row number 
      pointer [m] = count;
      cardinal[m] = ordering2D.cardinal[m2D];
      count+=cardinal[m];
      int pos=0;
      for(int j = 0; j < cardinal[m]; j++) {
        int l=max(k-1,0);
        size_t nz2D,nz3D;
        nz2D=ordering2D.incidence[ordering2D.pointer[m2D]+j];
        nz3D=structured3DIndex(nz2D, l, ordering2D.ncols, ordering2D.colsize);
        tmp[m*r.max + pos] = nz3D;
        pos++;
        }
      }
    }
  pointer[nrows] = count;

/*------------------------------------------------------------------------------
  allocate incidence and copy tmp */
  incidence = new int[count];

  for(n = 0; n < nrows; n++) {
    for(k = 0; k < cardinal[n]; k++) {
       /*incidence= row index */
      incidence[pointer[n] + k] = tmp[n * r.max + k];
      if(incidence[pointer[n] + k] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }

  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows     = nrows;
  ordering->ncols     = ordering2D.ncols*ordering2D.colsize;

//   status= matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,nrows);
  
  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_Hconcat_01(const hyperzmatrix_t & M1, const hyperzmatrix_t & M2, hyperzmatrix_t & M3)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  concat two matrices with same rows (or projection space) but cols corresponding
  to disjoint state vectors. 
  
  Concatenations will result in a extended state vector system
  
  No 3D expansion performed

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  
  if(M1.ordering->nrows!=M2.ordering->nrows) {
    check_error(-1,"inconsistent number of rows", __LINE__, __FILE__, 1);
    }

/*------------------------------------------------------------------------------
  at the moment, different colsize can be an issue for 3D expansion */
  if(M1.ordering->colsize!=M2.ordering->colsize) {
    check_error(-1,"inconsistent dimension for vertical expansion", __LINE__, __FILE__, 1);
    }

  size_t ncols=M1.ordering->ncols+M2.ordering->ncols;
  size_t nrows=M1.ordering->nrows;
  size_t size =M1.ordering->pointer[nrows]+M2.ordering->pointer[nrows];
  
  size_t offset=0;
  offset=M1.ordering->ncols*M1.ordering->colsize;
  
  vptr_t *pointer;
  
/*------------------------------------------------------------------------------
  first concat ordering */
  if(M3.ordering==0) {
  int *cardinal =new int[nrows];
  pointer  =new vptr_t[nrows+1];
  int *incidence=new int[size];
  
  for(size_t n=0;n<M1.ordering->nrows;n++) {
    cardinal[n]=M1.ordering->cardinal[n]+M2.ordering->cardinal[n];
    }
    
  pointer[0]=0;
  for(size_t n=0;n<M1.ordering->nrows;n++) {
    size_t count=pointer[n];
    for(int k=0;k<M1.ordering->cardinal[n];k++) {
      incidence[count]=M1.ordering->incidence[M1.ordering->pointer[n]+k];
      count++;
      }
    for(int k=0;k<M2.ordering->cardinal[n];k++) {
      incidence[count]=M2.ordering->incidence[M2.ordering->pointer[n]+k]+offset;
      count++;
      }
    pointer[n+1]=pointer[n]+cardinal[n];
    }
    
/*------------------------------------------------------------------------------
  keep all setting from M1/M2 */
  M3.ordering=new ordering_t;
  
  *(M3.ordering)=*(M1.ordering);
  
/*------------------------------------------------------------------------------
  except arrays */
  M3.ordering->pointer=pointer;
  M3.ordering->cardinal=cardinal;
  M3.ordering->incidence=incidence;
  M3.ordering->ncols=ncols;
//   M3.ordering->colsize=M1.ordering->colsize+M2.ordering->colsize;
  status=matrix_check_consistency(*M3.ordering);
  }
  else pointer=M3.ordering->pointer;
  
/*------------------------------------------------------------------------------
  then concat matrix */
  M3.allocate();
  M3.assign(0.0);
  
  for(size_t n=0;n<M1.ordering->nrows;n++) {
    size_t count=pointer[n];
    for(int k=0;k<M1.ordering->cardinal[n];k++) {
      M3.packed[count]=M1.packed[M1.ordering->pointer[n]+k];
      count++;
      }
    for(int k=0;k<M2.ordering->cardinal[n];k++) {
      M3.packed[count]=M2.packed[M2.ordering->pointer[n]+k];
      count++;
      }
    }
    
  
  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_Hconcat_02(hyperzmatrix_t & M1, hyperzmatrix_t & M2, hyperzmatrix_t & M3)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  concat two matrices with same rows (or projection space) but cols corresponding
  to disjoint state vectors. 
  
  Concatenations will result in a extended state vector system
  
  3D expansion performed

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  bool debug=false;
  
  if(M1.ordering->nrows!=M2.ordering->nrows) {
    check_error(-1,"inconsistent number of rows", __LINE__, __FILE__, 1);
    }
  
  if(debug) printf("%s : expand 1st matrix\n",__FUNCTION__);
  ordering_t *tmp1=new ordering_t;
  status=ordering_expand3D(*M1.ordering, tmp1);
  
  M1.ordering->destroy();
  *M1.ordering=*tmp1;

  if(debug) printf("%s : expand 2nd matrix\n",__FUNCTION__);
  ordering_t *tmp2=new ordering_t;
  status=ordering_expand3D(*M2.ordering, tmp2);

  M2.ordering->destroy();
  *M2.ordering=*tmp2;

  if(debug) printf("%s : concat matrix\n",__FUNCTION__);
  status=matrix_Hconcat_01(M1, M2, M3);
  
  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int matrix_Hconcat(hyperzmatrix_t & M1, hyperzmatrix_t & M2, hyperzmatrix_t & M3)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  bool expand=true;
  
  if(not expand) {
//   if(M1.ordering->colsize==M2.ordering->colsize) {
/*------------------------------------------------------------------------------
    no 3D expansion needed */
    status=matrix_Hconcat_01(M1, M2, M3);
    }
  else {
/*------------------------------------------------------------------------------
    3D expansion mandatory */
    status=matrix_Hconcat_02(M1, M2, M3);
    }
  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ordering_Vconcat(ordering_t *M1, ordering_t *M2, ordering_t *M3)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  expand and concat two ordering with similar cols but disjoint projection space 
  
  Concatenations will result in a extended state vector system
      
  if 2D nodes connected, assume corresponding expanded 3D nodes to be connected

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  size_t m;
  
  if(M1->ncols!=M2->ncols) return (-1);
  
  size_t ncols=M1->ncols;
/*------------------------------------------------------------------------------
  total number of 3D rows*/
  size_t nrows=M1->nrows*M1->rowsize+M2->nrows*M2->rowsize;
  
  size_t offset=0;
  
  int *cardinal =new int[nrows];
  vptr_t *pointer  =new vptr_t[nrows+1];
  
  pointer[0]=0;
  for(size_t n=0;n<M1->nrows;n++) {
    for(int k=0;k<M1->rowsize;k++) {
      m=n*M1->rowsize+k;
/*------------------------------------------------------------------------------
      */
//       cardinal[m]=M1->cardinal[n]*M1->rowsize+M2->cardinal[n]*M2->rowsize;
      cardinal[m]=M1->cardinal[n]*M1->rowsize;
      pointer[m+1]=pointer[m]+cardinal[m];
      }
    }
  
  offset=M1->nrows*M1->rowsize;
  for(size_t n=0;n<M2->nrows;n++) {
    for(int k=0;k<M2->rowsize;k++) {
      m=n*M2->rowsize+k+offset;
/*------------------------------------------------------------------------------
      */
//       cardinal[m]=M1->cardinal[n]*M1->rowsize+M2->cardinal[n]*M2->rowsize;
      cardinal[m]=M2->cardinal[n]*M2->rowsize;
      pointer[m+1]=pointer[m]+cardinal[m];
      }
    }
  
  size_t size=pointer[nrows];
  int *incidence=new int[size];
  
  offset=M1->nrows*M1->rowsize;
/*------------------------------------------------------------------------------
  2D rows loop */
  for(size_t n=0;n<M1->nrows;n++) {
/*------------------------------------------------------------------------------
    expand to 3D rows loop */
    for(int k=0;k<M1->rowsize;k++) {
      size_t m=n*M1->rowsize+k;
      size_t count=pointer[m];
/*------------------------------------------------------------------------------
      2D cols loop */
      for(int i=0;i<M1->cardinal[n];i++) {
        size_t n2D=M1->incidence[M1->pointer[n]+i];
/*------------------------------------------------------------------------------
        expand to 3D cols loop */
        for(int l=0;l<M1->rowsize;l++) {
          incidence[count]=n2D*M1->rowsize+l;
          count++;
          }
        }
//       for(int i=0;i<M2->cardinal[n];i++) {
//         size_t n2D=M2->incidence[M2->pointer[n]+i];
//         for(int l=0;l<M2->rowsize;l++) {
//           incidence[count]=n2D*M2->rowsize+l+offset;
//           count++;
//           }
//         }
      if(count-pointer[m]!=cardinal[m]) {
        check_error(-1,"ordering V-concatenation anomaly", __LINE__, __FILE__, 1);
        }
      }
    }
    
  offset=M1->nrows*M1->rowsize;
  for(size_t n=0;n<M2->nrows;n++) {
    for(int k=0;k<M2->rowsize;k++) {
      size_t m=n*M2->rowsize+k+offset;
      size_t count=pointer[m];
//       for(int i=0;i<M1->cardinal[n];i++) {
//         for(int l=0;l<M1->rowsize;l++) {
// 	  size_t n2D=M1->incidence[M1->pointer[n]+i];
//           incidence[count]=n2D*M1->rowsize+l;
//           count++;
//           }
//         }
      for(int i=0;i<M2->cardinal[n];i++) {
        for(int l=0;l<M2->rowsize;l++) {
	  size_t n2D=M2->incidence[M2->pointer[n]+i];
          incidence[count]=n2D*M2->rowsize+l;
          count++;
          }
        }
      if(count-pointer[m]!=cardinal[m]) {
        check_error(-1,"ordering V-concatenation anomaly", __LINE__, __FILE__, 1);
        }
      }
    }

/*------------------------------------------------------------------------------
  */
  M3->pointer=pointer;
  M3->cardinal=cardinal;
  M3->incidence=incidence;
  
  M3->nrows=nrows;
  M3->ncols=ncols;
    
  status=matrix_check_consistency(*M3);
//   status=matrix_check_symmetry(*M3);
  return (0);
}
