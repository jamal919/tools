
/**************************************************************************

  POC-SOLVERS interface, 2006-2012

  Part of the Unstructured Ocean Grid initiative and SYMPHONIE suite

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      PhD, LEGOS, Toulouse, France

***************************************************************************/

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
 OLD_SOLVER_VERSION flag:
 
   if activated, use historical (1.0) code (still available in the c-code branch)
   not maintained, kept only for testing purposes, aimed to be removed
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#define OLD_SOLVER_VERSION 0

#include <stdio.h>
#include <stdbool.h>

#ifndef SOLVER_INTERFACE_H
#define SOLVER_INTERFACE_H

#define SOLVER_VERSION 2

#ifdef HAVE_MPI
#include <mpi.h>
#else
#define MPI_COMM_WORLD 0
#define MPI_Comm        size_t
#endif

#include "solvers-functions.h"

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  solvers include files
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#include "cs.h"

#ifdef MUMPS
#include "dmumps_c.h"
#include "zmumps_c.h" 
#endif

#ifdef __cplusplus
#include <complex>
#define complex_t complex<double>
#else
#include <complex.h>
#define complex_t complex<double>
#endif

#ifdef UMFPACK
#ifdef __cplusplus
extern "C" {
#endif
#include <umfpack.h>
#ifdef __cplusplus
}
#endif

#endif

#ifdef HYPRE
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#endif

#ifdef HIPS

#ifdef __cplusplus
extern "C" {
#endif
  
#include <dhips.h>
#include <zhips.h>
  
#ifdef __cplusplus
}
#endif

#endif

#ifdef PASTIX
#include <stdint.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif
  
#include <pastix.h>
#include <read_matrix.h>
#include <get_options.h>
  
#ifdef __cplusplus
}
#endif

#endif


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  following classes are made to ease solver data handling in a compact manner
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#ifdef LAPACK
/* see LAPACK's doxygen documentation on
  http://www.netlib.org/lapack/explore-html/index.html */
#ifdef __cplusplus
extern "C" {
#endif
extern void dgetrf_(const int *m,const int *n,double    *A,const int *ldA,int *pivot,int *status);
extern void zgetrf_(const int *m,const int *n,complex_t *A,const int *ldA,int *pivot,int *status);
extern void dgetrs_(const char *trans,const int *nA,const int *nrhs,const double    *A,const int *ldA,const int *pivot,double    *b,const int *ldb,int *status);
extern void zgetrs_(const char *trans,const int *nA,const int *nrhs,const complex_t *A,const int *ldA,const int *pivot,complex_t *b,const int *ldb,int *status);
#define ZLAPACK 1

extern void dgbtrf_(int *,int *,int *,int *, double *, int *,  int  *, int *);
extern void zgbtrf_(int *,int *,int *,int *, complex_t *, int *,  int  *, int *);
#ifdef __cplusplus
extern void dgbtrs_(char *,int *,int *,int *,int *, double *, int *,int *,  double *, int  *, int *);
extern void zgbtrs_(char *,int *,int *,int *,int *, complex_t *, int *,int *,  complex_t *, int  *, int *);
#else
extern void dgbtrs_(char *trans, int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, int *ipiv, double *b, int *ldb, int *info, int trans_len);
extern void zgbtrs_(char *trans, int *n, int *kl, int *ku, int *nrhs, complex_t *ab, int *ldab, int *ipiv, complex_t *b, int *ldb, int *info, int trans_len);
#endif

extern void dpbsv_ (char *, int *, int *, int *, double *, int *, double *, int *, int *);
extern int  dpotrf_(char *, int *, double *, int *, int *);

extern void dsterf_(int *, double *, double *, int *);
extern void dstevd_(char *, int *, double *, double *, double *, int *, double *, int *, int *, int *, int *);
extern void dsygvd_(int *, char *, char *, int *, double *, int *, double *, int *, double *, double *, int *, int *, int *, int *);

extern  void dgeev_(char *, char *, int *, double *, int *, double *, double *, double *, int *, double *, int *, double *, int *, int *);
extern  void dggev_(char *, char *, int *, double *, int *, double *, int *, double *, double *, double *, double *, int *, double *, int *, double *, int *, int *);

#define LAPACKF_ 1
#ifdef __cplusplus
}
#endif
#endif

#if HAVE_LIBGSLCBLAS == 1
#include <gsl/gsl_cblas.h>
#define CBLAS 1
#endif

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  CXSPARSE convenience class for c++-code version
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <class T> class cxsparse_cs_t { /* matrix in compressed-column or triplet form */
private :
public  :
    int nzmax ;     /* maximum number of entries */
    int m ;         /* number of rows */
    int n ;         /* number of columns */
    int *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
    int *i ;        /* row indices, size nzmax */
    T *x ;          /* numerical values, size nzmax */
    int nz ;        /* # of entries in triplet matrix, -1 for compressed-col */
};

template <class T> class  cxsparse_symbolic_t  /* symbolic Cholesky, LU, or QR analysis */
{
private :
public  :
    int *pinv ;     /* inverse row perm. for QR, fill red. perm for Chol */
    int *q ;        /* fill-reducing column permutation for LU and QR */
    int *parent ;   /* elimination tree for Cholesky and QR */
    int *cp ;       /* column pointers for Cholesky, row counts for QR */
    int *leftmost ; /* leftmost[i] = min(find(A(i,:))), for QR */
    int m2 ;        /* # of rows for QR, after adding fictitious rows */
    double lnz ;    /* # entries in L for LU or Cholesky; in V for QR */
    double unz ;    /* # entries in U for LU; in R for QR */
};

template <class T> class  cxsparse_numeric_t   /* numeric Cholesky, LU, or QR factorization */
{
private :
public  :
    cxsparse_cs_t<T> *L ;      /* L for LU and Cholesky, V for QR */
    cxsparse_cs_t<T> *U ;      /* U for LU, r for QR, not used for Cholesky */
    int *pinv ;     /* partial pivoting for LU */
    double *B ;     /* beta [0..n-1] for QR */
};

template <class T> class cxsparse_t { /* matrix in compressed-column or triplet form */
private :
public  :
    cxsparse_cs_t<T>       *matrix;
    cxsparse_symbolic_t<T> *symbolic;
    cxsparse_numeric_t<T>  *numeric;
};

#ifdef LAPACK

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  LAPACK convenience class for c++-code version
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <class T> class lapack_t {
private :
public  :
  int tmp;
};

#endif

#ifdef UMFPACK

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  UMFPACK historical convenience structure for c-code version, soon obsolete
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

typedef struct
{
  double Info    [UMFPACK_INFO];
  double Control [UMFPACK_CONTROL];
  void *Symbolic, *Numeric ; 
} UMFPACK_STRUC_C;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  UMFPACK convenience class for c++-code version
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <class T> class umfpack_t {
private :
public  :
  double Info    [UMFPACK_INFO];
  double Control [UMFPACK_CONTROL];
  void *Symbolic, *Numeric;
  double *matrix_real, *matrix_imag;
  umfpack_t () {
    Symbolic=Numeric=0;
    matrix_real=matrix_imag=0;
    }
};

#endif


#ifdef HYPRE
typedef struct
{   
   HYPRE_IJMatrix A;
   HYPRE_ParCSRMatrix parcsr_A;
   HYPRE_IJVector b;
   HYPRE_ParVector par_b;
   HYPRE_IJVector x;
   HYPRE_ParVector par_x;
   HYPRE_Solver hy_solver, hy_precond;
}  HYPRE_STRUC_C;
#endif

#ifdef HIPS

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  HIPS historical convenience structure for c-code version, soon obsolete
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

typedef struct
{
  INTS id_hips, idnbr, i, j;
  INTS *unknownlist;
  double *x, *rhsloc;
  INTS proc_id, n, ln;
  INTL *ia, nnz;
  INTS *ja;
  double *a;
  INTS domsize, nproc;
  INTS pbegin, pend;
  INTS ierr;
} dHIPS_STRUC_C;

typedef struct
{
  INTS id_hips, idnbr, i, j;
  INTS *unknownlist;
  complex_t *x, *rhsloc;
  INTS proc_id, n, ln;
  INTL *ia, nnz;
  INTS *ja;
  complex_t *a;
  INTS domsize, nproc;
  INTS pbegin, pend;
  INTS ierr;
} zHIPS_STRUC_C;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  HIPS convenience class for c++-code version
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <class T> class hips_t {
private :
public  :
  INTS id_hips, idnbr, i, j;
  INTS *unknownlist;
  T *x, *rhsloc;
  INTS proc_id, n, ln;
  INTL *ia, nnz;
  INTS *ja;
  T *a;
  INTS domsize, nproc;
  INTS pbegin, pend;
  INTS ierr;
  
  hips_t() {
    unknownlist=0;
    x=0;
    rhsloc=0;
    ia=0;
    ja=0;
    a=0;
    }
  
  void destroy() {
    ptr_delete(unknownlist);
    }
};

#endif

#ifdef MUMPS

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  MUMPS convenience class for c++-code version
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <class T> class mumps_t {
private :
public  :
    MUMPS_INT      sym, par, job;
    MUMPS_INT      comm_fortran;    /* Fortran communicator */
    MUMPS_INT      icntl[40];
    MUMPS_INT      keep[500];
    ZMUMPS_REAL    cntl[15];
    ZMUMPS_REAL    dkeep[230];
    MUMPS_INT8     keep8[150];
    MUMPS_INT      n;

    MUMPS_INT      nz_alloc; /* used in matlab interface to decide if we
                                free + malloc when we have large variation */

    /* Assembled entry */
    MUMPS_INT      nz;
    MUMPS_INT8     nnz;
    MUMPS_INT      *irn;
    MUMPS_INT      *jcn;
    T *a;

    /* Distributed entry */
    MUMPS_INT      nz_loc;
    MUMPS_INT8     nnz_loc;
    MUMPS_INT      *irn_loc;
    MUMPS_INT      *jcn_loc;
    T *a_loc;

    /* Element entry */
    MUMPS_INT      nelt;
    MUMPS_INT      *eltptr;
    MUMPS_INT      *eltvar;
    T *a_elt;

    /* Ordering, if given by user */
    MUMPS_INT      *perm_in;

    /* Orderings returned to user */
    MUMPS_INT      *sym_perm;    /* symmetric permutation */
    MUMPS_INT      *uns_perm;    /* column permutation */

    /* Scaling (inout but complicated) */
    ZMUMPS_REAL    *colsca;
    ZMUMPS_REAL    *rowsca;
    MUMPS_INT colsca_from_mumps;
    MUMPS_INT rowsca_from_mumps;

    /* RHS, solution, ouptput data and statistics */
    T *rhs, *redrhs, *rhs_sparse, *sol_loc;
    MUMPS_INT      *irhs_sparse, *irhs_ptr, *isol_loc;
    MUMPS_INT      nrhs, lrhs, lredrhs, nz_rhs, lsol_loc;
    MUMPS_INT      schur_mloc, schur_nloc, schur_lld;
    MUMPS_INT      mblock, nblock, nprow, npcol;
    MUMPS_INT      info[40],infog[40];
    ZMUMPS_REAL    rinfo[40], rinfog[40];

    /* Null space */
    MUMPS_INT      deficiency;
    MUMPS_INT      *pivnul_list;
    MUMPS_INT      *mapping;

    /* Schur */
    MUMPS_INT      size_schur;
    MUMPS_INT      *listvar_schur;
    T *schur;

    /* Internal parameters */
    MUMPS_INT      instance_number;
    T *wk_user;

    /* Version number: length=14 in FORTRAN + 1 for final \0 + 1 for alignment */
    char version_number[MUMPS_VERSION_MAX_LEN + 1 + 1];
    /* For out-of-core */
    char ooc_tmpdir[256];
    char ooc_prefix[64];
    /* To save the matrix in matrix market format */
    char write_problem[256];
    MUMPS_INT      lwk_user;
    /* For save/restore feature */
    char save_dir[256];
    char save_prefix[256];
  
    int init(){
      }
};

extern void mumps_c(mumps_t<double> *id);

extern void mumps_c(mumps_t< complex<double> > *id) ;

#endif

#ifdef PASTIX

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  PASTIX historical convenience structure for c-code version, soon obsolete
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

typedef struct
{
  pastix_data_t  *pastix_data ; /* Pointer to a storage structure needed by pastix           */
  pastix_int_t    ncol;         /* Size of the matrix                                        */
  pastix_int_t   *colptr      ; /* Indexes of first element of each column in row and values */
  pastix_int_t   *rows        ; /* Row of each element of the matrix                         */
  double *values      ; /* Value of each element of the matrix                       */
  double *rhs         ; /* right hand side                                           */
  double *rhssaved    ; /* right hand side (save)                                    */
  double *ax          ; /* A times X product                                         */ 
  pastix_int_t    iparm[IPARM_SIZE];  /* integer parameters for pastix                             */
  double          dparm[DPARM_SIZE];  /* floating parameters for pastix                            */
  pastix_int_t   *perm        ; /* Permutation tabular                                       */
  pastix_int_t   *invp        ; /* Reverse permutation tabular                               */
  char           *type        ; /* type of the matrix                                        */
  char           *rhstype     ; /* type of the right hand side                               */
#ifdef HAVE_MPI
  int             required;           /* MPI thread level required                                 */
  int             provided;           /* MPI thread level provided                                 */
  int             rank;
#endif
  driver_type_t  *driver_type;        /* Matrix driver(s) requested by user                        */
  char          **filename;           /* Filename(s) given by user                                 */
  int             nbmatrices;         /* Number of matrices given by user                          */
  int             nbthread;           /* Number of thread wanted by user                           */
  int             verbosemode;        /* Level of verbose mode (0, 1, 2)                           */
  int             ordering;           /* Ordering to use                                           */
  int             nbrhs;
  int             incomplete;         /* Indicate if we want to use incomplete factorisation       */
  int             level_of_fill;      /* Level of fill for incomplete factorisation                */
  int             amalgamation;       /* Level of amalgamation for Kass                            */
  int             ooc;                /* OOC limit (Mo/percent depending on compilation options)   */
  long            i,j;
  double norme1, norme2;
} PASTIX_STRUC_C;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  PASTIX convenience class for c++-code version
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <class T> class pastix_t {
private :
public  :
  pastix_data_t  *pastix_data ; /* Pointer to a storage structure needed by pastix           */
  pastix_int_t    ncol;         /* Size of the matrix                                        */
  pastix_int_t   *colptr      ; /* Indexes of first element of each column in row and values */
  pastix_int_t   *rows        ; /* Row of each element of the matrix                         */
  T *values      ; /* Value of each element of the matrix                       */
  T *rhs         ; /* right hand side                                           */
  T *rhssaved    ; /* right hand side (save)                                    */
  T *ax          ; /* A times X product                                         */ 
  pastix_int_t    iparm[IPARM_SIZE];  /* integer parameters for pastix                             */
  double          dparm[DPARM_SIZE];  /* floating parameters for pastix                            */
  pastix_int_t   *perm        ; /* Permutation tabular                                       */
  pastix_int_t   *invp        ; /* Reverse permutation tabular                               */
  char           *type        ; /* type of the matrix                                        */
  char           *rhstype     ; /* type of the right hand side                               */
#ifdef HAVE_MPI
  int             required;           /* MPI thread level required                                 */
  int             provided;           /* MPI thread level provided                                 */
  int             rank;
#endif
  driver_type_t  *driver_type;        /* Matrix driver(s) requested by user                        */
  char          **filename;           /* Filename(s) given by user                                 */
  int             nbmatrices;         /* Number of matrices given by user                          */
  int             nbthread;           /* Number of thread wanted by user                           */
  int             verbosemode;        /* Level of verbose mode (0, 1, 2)                           */
  int             ordering;           /* Ordering to use                                           */
  int             nbrhs;
  int             incomplete;         /* Indicate if we want to use incomplete factorisation       */
  int             level_of_fill;      /* Level of fill for incomplete factorisation                */
  int             amalgamation;       /* Level of amalgamation for Kass                            */
  int             ooc;                /* OOC limit (Mo/percent depending on compilation options)   */
  long            i,j;
  double norme1, norme2;
  
  void init() {
    pastix_data=0;
    ncol=0;
    colptr=rows=0;
    values=rhs=rhssaved=ax=0;
    perm=invp=0;
    }
    
   pastix_t<T> (){
     this->init();
   }
  
} ;
#endif  

#ifdef MAPHYS
#include "dmph_maphys_type_c.h"
#endif

#ifdef MAPHYS
typedef struct

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  historical convenience structure for c-code version, soon obsolete
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

{   
   DMPH_maphys_t_d mpahys_id;
   int n;
   int nnz;
   int *irn;
   int *jcn;
   double *a;
   double *rhs;
   double *sol;
  int *mydof;
  int *mysizeIntrf;
  int *myinterface;
  int *myidexVi;
  int *myptrindeVi;
  int *myindexintrf;
   
}  DMAPHYS_STRUC_C;

typedef struct
{   
   ZMPH_maphys_t_d mpahys_id;
   int n;
   int nnz;
   int *irn;
   int *jcn;
   complex_t *a;
   complex_t *rhs;
   complex_t *sol;   
   int *mydof;
   int *mysizeIntrf;
   int *myinterface;
   int *myidexVi;
   int *myptrindeVi;
   int *myindexintrf;
   
}  ZMAPHYS_STRUC_C; 

#endif

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  DIAGONAL convenience class for c++-code version
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <class T> class diagonal_t {
private :
public  :
int ndof;
};

#ifdef PETSC

#include "petscksp.h"
#include "petscmat.h" 
#include "petscvec.h" 

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  PETSC convenience class for c++-code version
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <class T> class petsc_t {
private :
public  :
  KSP ksp;
  Mat matrix;
  Vec RHS,x;
  
  void init() {
    ksp=0;
    matrix=0;
    RHS=0;
    x=0;
    }
    
   petsc_t<T> () {
     this->init();
   }
}; 

#endif



#endif
