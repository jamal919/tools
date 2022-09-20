
/**************************************************************************

  POC-SOLVERS interface, 2006-2012

  Part of the Unstructured Ocean Grid initiative and SYMPHONIE suite

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      PhD, LEGOS, Toulouse, France

***************************************************************************/


#ifndef MATRIX_IO_H

#define MATRIX_IO_H

#ifdef HAVE_MPI
#include <mpi.h>
#endif


#include "cs.h"

#define MAXLINE 1024
#define EMPTY (-1)

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#ifndef MAX
#define BOOLEAN(x) ((x) ? TRUE : FALSE)
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif

typedef struct triplet_struct
{
  int  nrow;   /* the matrix is nrow-by-ncol */
  int  ncol;
  long nzmax;  /* maximum number of entries in the matrix */
  long nnz;    /* number of nonzeros in the matrix */

  int ml;      /* bnd format */
  int mu; 
  int lda;

  int stype ;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
   Describes what parts of the matrix are considered:

  0: matrix is "unsymmetric": use both upper and lower triangular parts
     (the matrix may actually be symmetric in pattern and value, but
     both parts are explicitly stored and used).  May be square or
     rectangular.
 >0: matrix is square and symmetric.  Entries in the lower triangular
     part.
 <0: matrix is square and symmetric.  Entries in the upper triangular part.

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  int comptype ;       /* compression type csc csr coo bnd*/
                       /* coo=0 csc=1 csr=2 lapack-band=3*/
  int base ;           /*  0 : rows,cols from 0 to n-1 */
                       /*  1 : rows,cols from 1 to n */

  /* for coo, csc, csr format*/
  size_t *pointer64 ;     /* [0..nzmax-1],                    */
  int *i ;              /* [0..nzmax-1], the row indices    */
  int *j ;              /* [0..nzmax-1], the column indices */
  
  double *x ;          /* size nzmax or 2*nzmax, if present */
  int *loc2glob;       /* local2global colmn numer( resp row in csr) */
  /* for bnd format */
  double *a;
  int *pivot;

  /* SpDOMESTIC */
  cs  *A;
  css *S;
  csn *N;

  /* MAphys */
  int *mydof;
  int *mysizeIntrf;
  int *myinterface;
  int *myidexVi;
  int *myptrindeVi;
  int *myindexintrf;
  
  int private_memory;

} triplet ;

typedef struct tripletz_struct
{
  int  nrow;    /* the matrix is nrow-by-ncol              */
  int  ncol;    /* the matrix is nrow-by-ncol              */
  long nzmax;   /* maximum number of entries in the matrix */
  long nnz;     /* number of nonzeros in the matrix        */

  int ml;       /* band format                             */
  int mu; 
  int lda;

  int stype ;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
   Describes what parts of the matrix are considered:

  0: matrix is "unsymmetric": use both upper and lower triangular parts
     (the matrix may actually be symmetric in pattern and value, but
     both parts are explicitly stored and used).  May be square or
     rectangular.
 >0: matrix is square and symmetric.  Entries in the lower triangular
     part.
 <0: matrix is square and symmetric.  Entries in the upper triangular part.

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  int comptype ;       /* compression type csc csr coo bnd */
                       /* coo=0 csc=1 csr=2 lapack-band=3  */
  int base ;           /*  0 : rows,cols from 0 to n-1     */
                       /*  1 : rows,cols from 1 to n       */

  /* for coo, csc, csr format*/
  size_t *pointer64 ;     /* [0..nzmax-1],                    */
  int *i ;              /* [0..nzmax-1], the row indices    */
  int *j ;              /* [0..nzmax-1], the column indices */
  
  complex<double> *x ;		/* size nzmax or 2*nzmax, if present  */
  /* for bnd format */
  complex<double> *a;
  int *loc2glob;        /* local2global colmn numer( resp row in csr) */
  int *pivot;

  /* SpDOMESTIC */
  cs_ci  *A ;
  cs_cis *S ;
  cs_cin *N ;

  /* MAphys */
  int *mydof;
  int *mysizeIntrf;
  int *myinterface;
  int *myidexVi;
  int *myptrindeVi;
  int *myindexintrf;
  
  int private_memory;

} tripletz ;
    
/** Double precision **/   
#ifdef __cplusplus
extern "C" {
#endif
triplet *read_triplet(char *file,int wanaOne_based);
int read_vect(char *file,long sz,double *Tx);
int coo_issym(triplet *Mat);
int convert(triplet *Mat,int destform, int base, int verbose);

int coo2csc    (triplet *Mat, int offset) ;
int csr2coo    (triplet *Mat, int offset) ;
int csc2coo    (triplet *Mat, int offset) ;
int csc2coosym (triplet *Mat, int offset) ;
int getbwd     (triplet *Mat, int *ml, int *mu) ;
int csc2bnd    (triplet *Mat, int offset) ;
int coo_issym  (triplet *Mat) ;

/** Complex precision **/   
tripletz *read_tripletz(char *file,int wanaOne_based);
int read_vectz(char *file,long sz,complex<double> *Tx);
int cooz_issym(tripletz *Mat);
int convertz(tripletz *Mat,int destform, int base, int verbose);

int coo2cscz    (tripletz *Mat, int offset) ;
int csc2cooz    (tripletz *Mat, int offset) ;
int csr2cooz    (tripletz *Mat, int offset) ;
int csc2coosymz (tripletz *Mat, int offset) ;
int getbwdz     (tripletz *Mat, int *ml, int *mu) ;
int csc2bndz    (tripletz *Mat, int offset) ;
int cooz_issym  (tripletz *Mat);

#if HAVE_MPI
int d_coo2glob(MPI_Comm world, int ncol,int nnz,int *Ap, int *Ai,double *Ax,int *ncolglob, int *nnzglob, int **Apglobin,int **Aiglobin,double **Axglobin);
int d_globcoo2csc(MPI_Comm world, int ncolglob, int nnzglob, int *Apglob,int *Aiglob,double *Axglob, int **Ap,int **Ai,double **Ax);
#endif

int z_coo2glob(int ncol,int nnz,int *Ap, int *Ai,complex<double> *Ax, int *ncolglob, int *nnzglob, int **Apglobin,int **Aiglobin,complex<double> **Axglobin);
int z_globcoo2csc(int ncolglob, int nnzglob, int *Apglob,int *Aiglob,complex<double> *Axglob,int **Ap,int **Ai,complex<double> **Ax);

extern  char *get_SolverName(int solver_id);

#ifdef __cplusplus
}
#endif

#endif
