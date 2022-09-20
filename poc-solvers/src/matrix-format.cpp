
/**************************************************************************

  POC-SOLVERS interface, 2006-2012

  Part of the Unstructured Ocean Grid initiative and SYMPHONIE suite

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      PhD, LEGOS, Toulouse, France

***************************************************************************/

#include <string.h>
#include <omp.h>
#include <complex>
#include <vector>

#include "solvers-functions.h"
#include "solvers-interface.h"
#include "solverlib.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int matrix_relative_pos(packed_t<T> *packed, int row, int target, int StopOnError)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  find the index of column "target" in row line

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int column;
  int offset=packed->numbering;
  size_t cardinal;
  
  size_t pos=packed->pointer[row]-offset;

  cardinal=packed->pointer[row+1]-packed->pointer[row];
  
  column=vpos((int) target,&(packed->incidence[pos]),cardinal);
  
  if(column==-1) {
    printf("\nseeking node %d in matrix line %d (offset=%d, type=%d)\n", target, row, offset, packed->ordering);
    if(StopOnError==1) exit(-1);
    }
  return(column);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 template <typename T> int matrix_CSR2CSC_template(SOLVER_t<T> & solver, packed_t<T> & CSCpacked, int target_numbering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  int row, offset;
  T *tmp;
  size_t nvalues;
  
  if(solver.packed==0) return(-1);

  if(solver.packed->ordering!=CSR) {
    return(-1);
    }
    
  if(solver.packed->ordering==CSC) {
    return(0);
    }

  nvalues=solver.packed->pointer[solver.packed->nrows];
  tmp=new T[nvalues];
  
  offset=target_numbering-solver.packed->numbering;
  
// #pragma omp parallel for if(gOPENMP_nCPUs>1)
  for(row = 0; row < solver.packed->nrows; row++) {
    for(int j = solver.packed->pointer[row]-solver.packed->numbering; j < solver.packed->pointer[row + 1]-solver.packed->numbering; j++) {
      int col=solver.packed->incidence[j]+offset;
      int i=matrix_relative_pos(&CSCpacked,col,row+target_numbering,1);
      size_t m=solver.packed->pointer[col]+i+offset;
      tmp[m]=solver.packed->x[j];
      }
    }

  delete[] solver.packed->x;

  solver.packed->x=tmp;
  solver.packed->ordering=CSC;
  
  solver.packed->numbering=target_numbering;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int matrix_CSR2CSC(SOLVER_t<double> & solver, packed_t<double> CSCpacked, int target_numbering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=matrix_CSR2CSC_template(solver, CSCpacked, target_numbering);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int matrix_CSR2CSC(SOLVER_t< complex<double> > & solver, packed_t< complex<double> > CSCpacked, int target_numbering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=matrix_CSR2CSC_template(solver, CSCpacked, target_numbering);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int matrix_CSR2COO_template(SOLVER_t<T> & solver, triplet_t<T> COOtriplet, int target_numbering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status = 0;
  int i;
  size_t k;
  int nrows;
  size_t *pointer;
  int *incidence;
  size_t nnz;
  int offset;
  
  if(solver.packed->ordering!=CSR) return(-1);

  solver.COOtriplet->destroy();

  nrows = solver.packed->nrows;

/*------------------------------------------------------------------------------
  init triplet */
  solver.COOtriplet->ordering=COO_ROW_MAJOR;
  solver.COOtriplet->numbering=target_numbering;

  solver.COOtriplet->stype=0;

  solver.COOtriplet->nnz=solver.packed->pointer[nrows];
  solver.COOtriplet->nrows=nrows;
  solver.COOtriplet->ncols=nrows;

  pointer   = solver.packed->pointer;
  incidence = solver.packed->inidence;

  solver.COOtriplet->i = new int[nnz];
  solver.COOtriplet->j = new int[nnz];
  solver.COOtriplet->x = new double[nnz];

/*------------------------------------------------------------------------------
  MUMPS needs COO one based matrix */
  solver.COOtriplet->comptype=COO;
 
  offset=target_numbering-solver.packed->numbering;
  
  for (i=0; i <nrows; i++) {
    for (k=pointer[i]-solver.packed->numbering; k < pointer[i+1]-solver.packed->numbering; k++) {
      solver.COOtriplet->i[k] = i            + offset;
      solver.COOtriplet->j[k] = incidence[k] + offset;
      solver.COOtriplet->x[k] = solver.packed->x[k];
      }
    }
   
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class T> int packed2packed(SOLVER_t<T> & solver, int targeted_ordering, int targeted_numbering, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int offset;
  int status=-1;
      
  offset = targeted_numbering - solver.packed->numbering;
  
  switch(solver.packed->ordering) {
    case CSR:
      switch(targeted_ordering) {
        case CSR: 
          status=0;
          break;
        case CSC:
          printf("conversion from %d ordering toward %d ordering not yet implemented\n",solver.packed->ordering, targeted_ordering);
          break;
        default:
          printf("conversion from %d ordering toward %d ordering not yet implemented\n",solver.packed->ordering, targeted_ordering);
          break;
        }
      break;
      
    case CSC:
      switch(targeted_ordering) {
        case CSR: 
          printf("conversion from %d ordering toward %d ordering not yet implemented\n",solver.packed->ordering, targeted_ordering);
          status=0;
          break;
        case CSC:
          status=0;
          break;
        default:
          printf("conversion from %d ordering toward %d ordering not yet implemented\n",solver.packed->ordering, targeted_ordering);
          break;
        }
      break;
    default:
      printf("ordering %d not yet implemented\n",solver.packed->ordering);
      break;
    }
  
  return(status);
}
  

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class T> int packed2coo(SOLVER_t<T> & solver, int targeted_ordering, int targeted_numbering, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int offset;
  int status=-1;
      
  offset = targeted_numbering - solver.packed->numbering;
  
  switch(solver.packed->ordering) {
    case CSR:
      switch(targeted_ordering) {
        case CSR: 
          status=0;
          break;
        case CSC:
          printf("(double matrix) conversion from coo to csr not yet implemented\n");
          break;
        default:
          printf("ordering %d not yet implemented\n",targeted_ordering);
          break;
        }
      break;
      
    case CSC:
      switch(targeted_ordering) {
        case CSR: 
          status=0;
          break;
        case CSC:
          printf("(double matrix) conversion from coo to csr not yet implemented\n");
          status=0;
          break;
        default:
          printf("ordering %d not yet implemented\n",targeted_ordering);
          break;
        }
      break;
    default:
      printf("ordering %d not yet implemented\n",solver.packed->ordering);
      break;
    }
  
  return(status);
}
  

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class T> int convert_template(SOLVER_t<T> & solver, int targeted_packing, int targeted_ordering, int targeted_numbering, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  aimed to perform matrix modification to suite current solver requirements
  
  not operational

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int offset;
  int status=-1;
  
  switch(solver.format) {
    case MATRIX_NATIVE:
      offset = targeted_numbering - solver.COOtriplet->numbering;
      switch(targeted_packing) {
        case COO: 
          status=0;
          break;
        case PACKED:
          status=0;
//           printf("(double matrix) conversion from coo to csr not yet implemented\n");
          break;
        case BAND:
          printf("(double matrix) conversion from coo to bnd not yet implemented\n");
          break;
        default:
          printf("(double matrix) conversion from coo to %d not yet implemented\n",targeted_ordering);
          break;
        }
      break;
      
    case MATRIX_PACKED:
      offset = targeted_numbering - solver.packed->numbering;
      switch(targeted_packing) {
        case COO: 
          status=0;
          break;
        case PACKED:
          status=packed2packed(solver, targeted_ordering, targeted_numbering, verbose, debug);
          break;
        case BAND:
          printf("(double matrix) conversion from coo to bnd not yet implemented\n");
          break;
        default:
          printf("(double matrix) conversion from coo to %d not yet implemented\n",targeted_ordering);
          break;
        }
      break;
      
    case MATRIX_BAND:
      switch(targeted_packing) {
        case COO: 
          status=0;
          break;
        case PACKED:
          printf("(double matrix) conversion from coo to csr not yet implemented\n");
          break;
        case BAND:
          printf("(double matrix) conversion from coo to bnd not yet implemented\n");
          break;
        default:
          printf("(double matrix) conversion from coo to %d not yet implemented\n",targeted_ordering);
          break;
        }
      break;
      
    case MATRIX_USER:
      status=0;
      break;
      
    default:
      printf("conversion from %d not yet implemented\n",solver.format);
      break;
    }
    
//   offset = base - Mat->base;
//   
//   if (solver.packed->comptype == destform && offset==0) {
// /**-----------------------------------------------------------------------------
//     nothing to do, assumes it was properly done by the user                   */
//     if(verbose==1) printf("matrix format conversion from %d toward %d with offset=%d, nothing to be done\n",Mat->comptype,destform,offset);
//     return(0);
//     }
// 
//   if(verbose==1)
//     printf("matrix format conversion from %d toward %d with offset=%d (source=%d target=%d)\n",Mat->comptype,destform,offset,Mat->base, base);
  
//   switch(Mat->comptype) {
//     
//     case COO:
//       switch(destform) {
//         case CSC: 
//           coo2csc(Mat,offset);
//           Mat->comptype = destform;
//           status=0;
//           break;
//         case CSR: /* To csr */
//           printf("(double matrix) conversion from coo to csr not yet implemented\n");
//           break;
//         case BND: /* To bnd */
//           printf("(double matrix) conversion from coo to bnd not yet implemented\n");
//           break;
//         case COO: /* To csr */
//           printf("(double matrix) conversion from coo to coo not yet implemented\n");
//           break;
//         case COO_COL_ARRANGED:
//           printf("(double matrix) conversion from coo to coo (COL arranged) not yet implemented, offset=%d\n",offset);
//           break;
//         default:
//           printf("(double matrix) conversion from coo to %d not yet implemented\n",destform);
//           break;
//         }
//       break;
//       
//     case CSC:
//       switch(destform) {
//         case COO:
//           if (Mat->stype >  0) csc2coosym(Mat, offset);
//           if (Mat->stype == 0) status=csc2coo(Mat, offset);
//           if(status!=0) return(status);
//           Mat->comptype = destform;
//           status=0;
//           break;
//         case COO_COL_ARRANGED: /* To coo for pastix */
//           if(verbose) printf("complex matrix conversion from csc to coo\n");
// /*------------------------------------------------------------------------------
//           poc-solvers refurbishing : bug found */
// //           csr2cooz(Mat, offset);
//           csr2coo(Mat, offset);
//           Mat->comptype = destform;
//           status=0;
//           break;
//         case BND:
//           status = csc2bnd(Mat, offset);
//           Mat->comptype = destform;
//           status=0;
//           break;
//         default:
//           printf("(double matrix) conversion from csc to %d not yet implemented\n",destform);
//           break;
//         }
//       break;
//       
//     case CSR: /* *** From csr *** */
//         switch(destform) {   
//           case COO: /* To coo */
//             if(verbose) printf("complex matrix conversion from csr to coo\n");
//             if (Mat->stype >  0) csc2coosym(Mat, offset);
//             else if (Mat->stype == 0) csr2coo(Mat, offset);
//             Mat->comptype = destform;
//             status=0;
//             break;
//           default:
//             printf("(complex) matrix conversion from csc to %d not yet implemented\n",destform);
//             break;
//           }
//       break;
//       
//     case BND: /* *** From bnd *** */
//       printf("(double matrix) conversion from BND not yet implemented\n",destform);
//       break;
//       
//     default : /*Not yet*/
//       printf("(double matrix) matrix format not yet supported \n");
//       break;
//     }
  return(status);
}
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int convert(SOLVER_t<double> & solver, int targeted_packing, int targeted_ordering, int targeted_numbering, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=convert_template(solver, targeted_packing, targeted_ordering, targeted_numbering, verbose, debug);
  return(status);
}
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int convert(SOLVER_t< complex<double> > & solver, int targeted_packing, int targeted_ordering, int targeted_numbering, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=convert_template(solver, targeted_packing, targeted_ordering, targeted_numbering, verbose, debug);
  return(status);
}

