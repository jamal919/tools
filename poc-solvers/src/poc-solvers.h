
/**************************************************************************

  POC-SOLVERS interface, 2006-2019

  Part of the SIROCCO national service (INSU, France)

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France

***************************************************************************/

#include <stdio.h>
#include <stdbool.h>

#ifndef POC_SOLVERS_H

#define POC_SOLVERS_H

#define POC_SOLVERS_VERSION_RELEASE 2
#define POC_SOLVERS_VERSION_MAJOR   1
#define POC_SOLVERS_VERSION_MINOR   0
#define POC_SOLVERS_LIB_2_0

#ifdef HAVE_MPI
#include <mpi.h>
#else
#define MPI_COMM_WORLD 0
#define MPI_Comm       size_t
#endif

#ifdef __cplusplus
#include <complex>
#define complex_t complex<double>
#else
#include <complex.h>
#define complex_t complex<double>
#endif

#include "solvers-functions.h"
#include "matrix_io.h"

#define POC_SOLVER_CXX 1

#define SEPARATOR_1 "\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
#define SEPARATOR_2 "\n-------------------------------------------------------------------------------\n"

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  matrix format
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#define COO 0
#define CSC 1
#define CSR 2
#define BND 3

#define COO_ROW_ARRANGED 0
#define COO_COL_ARRANGED 4

#define COO_ROW_MAJOR 0
#define COO_COL_MAJOR 4

#define COO_NAME "COO"
#define CSC_NAME "CSC"
#define CSR_NAME "CSR"
#define BND_NAME "BND"

#define PACKED 1
#define BAND   2

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  matrix format, revised
  
  CSR is MATRIX_PACKED, MATRIX_ROW_MAJOR, MATRIX_POINTER
  CSC is MATRIX_PACKED, MATRIX_COL_MAJOR, MATRIX_POINTER
  
  COO can be MATRIX_PACKED, MATRIX_ROW_MAJOR, MATRIX_TRIPLET (MUMPS, HIPS)
          or MATRIX_PACKED, MATRIX_COL_MAJOR, MATRIX_TRIPLET (PASTIX)
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
matrix value's kind of storage                                                */
#define MATRIX_NATIVE 0  // regular full matrix
#define MATRIX_PACKED 1  // regular packed matrix
#define MATRIX_BAND   2  // regular band matrix
#define MATRIX_USER  10  // user-defined matrix 

/*------------------------------------------------------------------------------
matrix value's kind of arrangement                                            */
#define MATRIX_ROW_MAJOR 0
#define MATRIX_COL_MAJOR 1

/*------------------------------------------------------------------------------
matrix value's kind of addressing                                             */
#define MATRIX_TRIPLET 0
#define MATRIX_POINTER 1
#define MATRIX_DIRECT  2

/*------------------------------------------------------------------------------
matrix value's kind of addressing                                             */
#define MATRIX_C_NUMBERING 0
#define MATRIX_F_NUMBERING 1



/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  some solvers get funny when handling diagnonal matrices (PASTIX and HIPS)
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#define SOLVER_ID_DIAGONAL   -99
#define SOLVER_ID_DOMESTIC     0
#define SOLVER_ID_LAPACK       1
#define SOLVER_ID_SUNPERF      2
#define SOLVER_ID_UMFPACK      3
#define SOLVER_ID_PETSC_GMRES  4
#define SOLVER_ID_MUMPS        5
#define SOLVER_ID_MUMPS_SYM    6
#define SOLVER_ID_SpDOMESTIC   7
#define SOLVER_ID_HYPRE        8
#define SOLVER_ID_HIPS         9
#define SOLVER_ID_PASTIX      10
#define SOLVER_ID_MAPHYS      11
#define SOLVER_ID_PASTIX_SEQUENTIAL 12

#define SOLVER_NAME_DIAGONAL      "DIAGONAL"
#define SOLVER_NAME_DOMESTIC      "DOMESTIC"
#define SOLVER_NAME_LAPACK        "LAPACK"
#define SOLVER_NAME_SUNPERF       "SUNPERF"
#define SOLVER_NAME_UMFPACK       "UMFPACK"
#define SOLVER_NAME_PETSC_GMRES   "PETSC_GMRES"
#define SOLVER_NAME_MUMPS         "MUMPS"
#define SOLVER_NAME_MUMPS_SYM     "MUMPS_SYM"
#define SOLVER_NAME_SpDOMESTIC    "SpDOMESTIC"
#define SOLVER_NAME_HYPRE         "HYPRE"
#define SOLVER_NAME_HIPS          "HIPS"
#define SOLVER_NAME_PASTIX        "PASTIX"
#define SOLVER_NAME_MAPHYS        "MAPHYS"

#define SOLVER_ID_TUGO_GMRES    100
#define SOLVER_NAME_TUGO_GMRES  "TUGO_GMRES"

/*------------------------------------------------------------------------------
  left for compatibility reasons; should be removed soon */
#define SOLVER_MATRIX_CENTRALIZED 0
#define SOLVER_MATRIX_DISTRIBUTED 1

#define SOLVER_RHS_CENTRALIZED 0
#define SOLVER_RHS_DISTRIBUTED 1

#define SOLVER_SOLUTION_CENTRALIZED 0
#define SOLVER_SOLUTION_DISTRIBUTED 1

/*----------------------------------------------------------------------------*/

// #define SOLVER_SEQUENTIAL -1
#define SOLVER_CENTRALIZED 0
#define SOLVER_DISTRIBUTED 1
#define SOLVER_OVERLAPPED  2


#ifdef HAVE_MPI
#define PASTIX_DISTRIBUTED 1
typedef struct {
  MPI_Comm communicator;
} mpi_handle_t;
#else
#define PASTIX_DISTRIBUTED 0
typedef struct {
  void *communicator;
} mpi_handle_t;
#endif

class solver_t {
private :
public  :
  char *name; 
  int  id;
/*------------------------------------------------------------------------------
  pointer to solvers specific anciliary data */
  void *parameters;  
  void *Mat;
  int RHS_distribution, solution_distribution, matrix_distribution;
/*------------------------------------------------------------------------------
  transpose flag added to ease pastix MPI distributed mode */
  int transpose_matrix;
  int format;
  MPI_Comm communicator;
  int rank;
};

template <class T> class triplet_t {
private :
public  :
  int nrows, ncols, ngdof;
  unsigned int *i ;     /* [0..nzmax-1], the row indices           */
  unsigned int *j ;     /* [0..nzmax-1], the column indices        */
  int ordering;
  int numbering;
  unsigned long nnzmax; /* maximum number of entries in the matrix */
  unsigned long nnz;    /* number of nonzeros in the matrix        */
  bool symmetric;
  T *x;
  
  void init() {
    nrows=0;
    ncols=0;
    ngdof=0;
    i=j=0;
    x=0;
    ordering=-1;
    numbering=-1;
    nnzmax=nnz=0;
    }
    
  triplet_t() {
    init();
    }
    
  void allocate() {
//     x=new T[nnz];
    i=new unsigned int[nnz];
    j=new unsigned int[nnz];
    }
    
  void destroy() {
    ptr_delete(i);
    ptr_delete(j);
    ptr_delete(x);
    init();
    }
};

template <class T> class band_t {
private :
public  :
  int nrows, ncols, ngdof;
  int ml;
  int mu; 
  int lda;
  int *pivot;
  T *x;
  
  band_t() {
    nrows=0;
    ncols=0;
    ngdof=0;
    ml=mu=lda=0;
    x=0;
    }
};

template <typename T> class packed_t {
private :
public  :
  int nrows, ncols, ngdof;
  int *incidence;
//   bool *solved;                /* number of non-zero values linked to solved columns (Petsc diagonal sub-matrix) */
  unsigned long *pointer;
  int *pointer32;
  int ordering;
  int numbering;
  unsigned long  nnzmax, nnz;
  int *loc2glob;
  bool symmetric;
  T *x;
  
  void init() {
    nrows=0;
    ncols=0;
    ngdof=0;
    incidence=0;
//     solved=0;
    pointer=0; 
    pointer32=0; 
    ordering=-1;
    numbering=-1;
    x=0;
    nnzmax=nnz=0;
    loc2glob=0;
    }
    
  packed_t() {
    init();
    }
    
  template <typename P> int cardinal(P* & c) {
    c=new P[nrows];
    for(size_t n=0;n<nrows;n++) {
      c[n]=(P) (pointer[n+1]-pointer[n]);
      }
    return(0);
    }
    
  template <typename P>  int set_nnz(P* & d_nnz, P* & o_nnz) {
//     if(solved==0) return(-1);
    d_nnz=new P[nrows];
    o_nnz=new P[nrows];
    for(size_t n=0; n<nrows; n++) {
      size_t count=0, cardinal=pointer[n+1]-pointer[n];
      for(size_t k=pointer[n]; k<pointer[n+1];k++) {
//         size_t m=incidence[k];
//         size_t pos=vpos((int) m,loc2glob,nrows);
        int m=incidence[k];
        int pos=vpos(m,loc2glob,nrows);
        if(pos==-1) {
//           printf("set_nnz: %d %d %d %d \n", n, loc2glob[n], m, k);
          }
        else {
//           if(not solved[pos]) {
//             printf("set_nnz anomaly: %d %d %d %d \n", n, loc2glob[n], m, k);
//             }
          count++;
          }
        }
      d_nnz[n]=count;
      o_nnz[n]=cardinal-count;
      }
    return(0);
    }
    
  void destroy() {
    ptr_delete(incidence);
//     ptr_delete(solved);
    ptr_delete(pointer);
    ptr_delete(pointer32);
    ptr_delete(loc2glob);
    ptr_delete(x);
    init();
    }
};

#define OPTIONS_SHOW_FACTORIZE_PERF  0
#define OPTIONS_SHOW_PRECOND_PERF    1
#define OPTIONS_SHOW_SOLVE_PERF      2
#define OPTIONS_SHOW_SOLVE_DIAGS     3

template <class T> class SOLVER_t {
private :
public  :
  char *name; 
  int  id;
/*------------------------------------------------------------------------------
  matrix storage informations */
  int format;
  int targeted_packing;          // packed, band or full
  int targeted_ordering;         // row-major or col-major
  int targeted_addressing;       // incidence pointer or row and col index
  int targeted_numbering;        // C (starting at 0) or Fortran (starting at 1)
/*------------------------------------------------------------------------------
  containers */
  triplet_t<T>   *COOtriplet;
  packed_t<T>    *packed;
  band_t<T>      *band;
  void *matrix;
/*------------------------------------------------------------------------------
  pointer to solvers specific anciliary data */
  void *parameters;  
  int RHS_distribution, solution_distribution, matrix_distribution;
//   distribution_t *native_distributor;
//   exchange_t *native_exchange;
/*------------------------------------------------------------------------------
  performance meter */
  void *perftimer;
/*------------------------------------------------------------------------------
  transpose flag added to ease pastix MPI distributed mode */
  int transpose_matrix;
/*------------------------------------------------------------------------------
  MPI data */
  MPI_Comm communicator;
  int rank, size;
  MPI_Comm world_communicator;
  int world_rank, world_size;
/*------------------------------------------------------------------------------
  specialised functions */
  int (*factorize)  (SOLVER_t<T> * solver, int verbose, bool debug);
  int (*solve)      (SOLVER_t<T> * solver, T *RHS, T *x, int transpose, int verbose, bool debug);
  int (*initialize) (SOLVER_t<T> & solver, int verbose, bool debug);
  int (*terminate)  (SOLVER_t<T> * solver, int verbose, bool debug);
  int (*error)      (int status, int verbose, bool debug);
/*------------------------------------------------------------------------------
  verbose */
  int options[100];
  
  void init() {
    name=0;
    id=-1;
    format=-1;
    targeted_packing=targeted_ordering=targeted_addressing=targeted_numbering=-1;
    COOtriplet=0;
    packed=0;
    band=0;
    matrix=0;
    parameters=0;
    RHS_distribution=solution_distribution=matrix_distribution=-1;
//     native_distributor=0;
//     native_exchange=0;
    perftimer=0;
    transpose_matrix=-1;
    communicator=0;
    rank=-1;
    factorize=0;
    solve=0;
    initialize=0;
    terminate=0;
    error=0;
    for(int k=0;k<100;k++) options[k]=-99;
    };
    
  SOLVER_t() {
    this->init();
    };
    
  void SetOption(int flag, int value) {
    options[flag]=value;
    }

  void destroy() {
    ptr_delete(parameters);
    ptr_delete(name);
    if(packed!=0) packed->destroy();
    if(COOtriplet!=0) COOtriplet->destroy();
    this->init();
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  solver interface compatibility with older versions
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  void d_import2packed(solver_t *slv) {
    int status;
    
    if(packed==0) packed=new packed_t<T>;
    
    triplet *Mat;
    Mat=(triplet *) slv->Mat;
    
    packed->incidence= (int *)  Mat->j;
    packed->x        = Mat->x;
  
    packed->nrows    = Mat->nrow;
    packed->nnz      = Mat->nnz;
    packed->loc2glob = Mat->loc2glob;

// //     packed->pointer  = (unsigned long *) Mat->i;
//     status=recast(Mat->i, packed->pointer, false, packed->nnz);
    
    packed->pointer32  = Mat->i;
    
    this->parameters=slv->parameters;
    
    this->transpose_matrix=slv->transpose_matrix;

    };
    
  void z_import2packed(solver_t *slv) {
    int status;
    
    if(packed==0) packed=new packed_t<T>;
    
    tripletz *Mat;
    Mat=(tripletz *) slv->Mat;
    
    packed->incidence= (int *)  Mat->j;
    packed->x        = Mat->x;
  
    packed->nrows    = Mat->nrow;
    packed->nnz      = Mat->nnz;
    packed->loc2glob = Mat->loc2glob;
 
// //     packed->pointer  = (unsigned long *) Mat->i;
//     status=recast(Mat->i, packed->pointer, false, packed->nnz);
    
    packed->pointer32  = Mat->i;
    
    this->parameters=slv->parameters;

    this->transpose_matrix=slv->transpose_matrix;
    
    };
    
  void d_export_packed(solver_t * slv) {
    int status;
    
    triplet *Mat;
    Mat=(triplet *) slv->Mat;
    
    if(packed->pointer!=0) {
      status=recast(packed->pointer, Mat->i, false, packed->nnz);
      }
    else Mat->i  = packed->pointer32; 
    Mat->j  = packed->incidence;
    Mat->x  = packed->x;
  
    Mat->nrow     = packed->nrows;
    Mat->nnz      = packed->nnz;
    Mat->loc2glob = packed->loc2glob;
    };
    
  void z_export_packed(solver_t * slv) {
    int status;
    
    tripletz *Mat;
    Mat=(tripletz *) slv->Mat;
    
    if(packed->pointer!=0) {
      status=recast(packed->pointer, Mat->i, false, packed->nnz);
      }
    else Mat->i  = packed->pointer32; 
    Mat->j  = packed->incidence;
    Mat->x  = packed->x;
  
    Mat->nrow     = packed->nrows;
    Mat->nnz      = packed->nnz;
    Mat->loc2glob = packed->loc2glob; 
    };
    
  void export_parameters(solver_t * slv) {
    slv->name=strdup(name); 
    slv->id=id;
    slv->format=format;
    slv->parameters=parameters;  
    slv->RHS_distribution=RHS_distribution;
    slv->solution_distribution=solution_distribution;
    slv->matrix_distribution=matrix_distribution;
    slv->transpose_matrix=transpose_matrix;
    slv->communicator=communicator;
    slv->rank=matrix_distribution;
    };
};

template <class T> class SOLVER_REG_t {
private :
public  :
  char *name; 
  int  id;
  int format, ordering;
  int (*factorize)  (SOLVER_t<T> * solver, int verbose, bool debug);
  int (*solve)      (SOLVER_t<T> * solver, T *, T *, int transpose, int verbose, bool debug);
  int (*initialize) (SOLVER_t<T> & solver, int verbose, bool debug);
  int (*terminate)  (SOLVER_t<T> * solver, int verbose, bool debug);
  int (*error)      (int status, int verbose, bool debug);
  SOLVER_REG_t() {
    name=0;
    id=-1;
    factorize=0;
    solve=0;
    initialize=0;
    terminate=0;
    error=0;
    };
    
  void destroy() {
    
    }
};
  
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
#ifdef __cplusplus
extern void dgbtrs_(char *,int *,int *,int *,int *, double *, int *,int *,  double *, int  *, int *);
#else
extern void dgbtrs_(char *trans, int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, int *ipiv, double *b, int *ldb, int *info, int trans_len);
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


/* Double precision routines */
#ifdef __cplusplus
extern "C" {
#define DEFAULT_ZERO =0
#else
#define DEFAULT_ZERO
#endif
  
int ishere(const char *solver_name); /* return 1/0  (solver is present/ not present) */
int get_SolverId(const char *solver_name, bool debug);
int LinearSystem_identify(const char *solver_name, bool debug);

int init_solverstruct(solver_t *solver);
int solver_terminate(solver_t *solveur);

void set_MPI_Init_thread_DONE(int);

solver_t *init_solver_obsolete(const char *solver_name,int typsym, int verbose DEFAULT_ZERO);
solver_t *init_solver(int solver_id, MPI_Comm world, int nprocs, int nunknowns, int verbose, bool debug);

int free_solver(solver_t *slv);

int factorize                   (solver_t *slv, triplet *Mat, int verbose DEFAULT_ZERO);
int factorize_mumps             (solver_t *slv, triplet *Mat, int verbose, bool debug);
int factorize_umfpack           (solver_t *slv, triplet *Mat,int verbose);
int factorize_lapack            (solver_t *slv, triplet *Mat);
int factorize_spdomestic        (solver_t *slv, triplet *Mat);
int factorize_hypre             (solver_t *slv, triplet *Mat);
int factorize_hips              (solver_t *slv, triplet *Mat, int verbose, bool debug);
int factorize_pastix_sequential (solver_t *slv, triplet *Mat, int verbose, bool debug);
int factorize_pastix_parallel   (solver_t *slv, triplet *Mat, int verbose, bool debug);
int factorize_pastix_parallel_distributed   (solver_t *slv, triplet *Mat, int verbose, bool debug);

int solve                       (solver_t *slv, double *RHS, int transp, int verbose DEFAULT_ZERO) ;
int solve_mumps                 (solver_t *slv, double *RHS, int transp, int verbose, bool debug) ;
int solve_umfpack               (solver_t *slv, double *RHS, int transp) ;
int solve_lapack                (solver_t *slv, double *RHS, int transp) ;
int solve_spdomestic            (solver_t *slv, double *RHS, int transp) ;
int solve_hypre                 (solver_t *slv, double *RHS, int transp) ;
int solve_hips                  (solver_t *slv, double *RHS, int transp, bool debug) ;
int solve_pastix                (solver_t *slv, double *RHS, int transp, bool debug) ;
int solve_pastix_parallel       (solver_t *slv, double *RHS, int transp, int verbose, bool debug) ;
int solve_pastix_sequential     (solver_t *slv, double *RHS, int transp) ;
#ifdef __cplusplus
}
#endif
/* Complexe precision routines */

#ifdef __cplusplus
extern "C" {
#endif
  
solver_t *initz_solver(char *solver_name,int typsym,int verbose DEFAULT_ZERO);

int freez_solver(solver_t *slv);

int factorizez                   (solver_t *slv, tripletz *Mat, int verbose DEFAULT_ZERO);
int factorizez_mumps             (solver_t *slv, tripletz *Mat, int verbose, bool debug);
int factorizez_umfpack           (solver_t *slv, tripletz *Mat, int verbose);
int factorizez_lapack            (solver_t *slv, tripletz *Mat);
int factorizez_spdomestic        (solver_t *slv, tripletz *Mat);
int factorizez_hips              (solver_t *slv, tripletz *Mat, int verbose, bool debug);
int factorizez_hips_complex      (solver_t *slv, tripletz *Mat, int verbose, bool debug);
int factorizez_pastix_sequential (solver_t *slv, tripletz *Mat, int verbose, bool debug);
int factorizez_pastix_parallel   (solver_t *slv, tripletz *Mat, int verbose, bool debug);
int factorizez_pastix_parallel_distributed   (solver_t *slv, tripletz *Mat, int verbose, bool debug);

int solvez                       (solver_t *slv, complex_t *RHS, int transp, int verbose DEFAULT_ZERO);
int solvez_mumps                 (solver_t *slv, complex_t *RHS, int transp, int verbose, bool debug) ;
int solvez_umfpack               (solver_t *slv, complex_t *RHS, int transp) ;
int solvez_lapack                (solver_t *slv, complex_t *RHS, int transp) ;
int solvez_spdomestic            (solver_t *slv, complex_t *RHS, int transp) ;
int solvez_hips                  (solver_t *slv, complex_t *RHS, int transp, int verbose, bool debug) ;
int solvez_hips_complex          (solver_t *slv, complex_t *RHS, int transp, int verbose, bool debug) ;
int solvez_pastix_sequential     (solver_t *slv, complex_t *RHS, int transp, int verbose, bool debug) ;
int solvez_pastix_parallel       (solver_t *slv, complex_t *RHS, int transp, int verbose, bool debug) ;
#ifdef __cplusplus
}
#endif

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  interface for internal use, c++ version

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

extern int mumps_factorize(SOLVER_t<double> *solver, int verbose, bool debug);
extern int mumps_factorize(SOLVER_t< complex<double> > *solver, int verbose, bool debug);

extern int mumps_solve(SOLVER_t<double> *solver, double *RHS, double *x,int transposed, int verbose, bool debug);
extern int mumps_solve(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x,int transposed, int verbose, bool debug);
 
extern int mumps_initialize(SOLVER_t<double> & solver, int verbose, bool debug);
extern int mumps_initialize(SOLVER_t< complex<double> > & solver, int verbose, bool debug);
  
extern int mumps_terminate(SOLVER_t<double> *solver, int verbose, bool debug);
extern int mumps_terminate(SOLVER_t< complex<double> > *solver, int verbose, bool debug);



extern int umfpack_factorize(SOLVER_t<double> *solver, int verbose, bool debug);
extern int umfpack_factorize(SOLVER_t< complex<double> > *solver, int verbose, bool debug);

extern int umfpack_solve(SOLVER_t<double> *solver, double *RHS, double *x,int transposed, int verbose, bool debug);
extern int umfpack_solve(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x,int transposed, int verbose, bool debug);

extern int umfpack_initialize(SOLVER_t<double> & solver, int verbose, bool debug);
extern int umfpack_initialize(SOLVER_t< complex<double> > & solver, int verbose, bool debug);
  
extern int umfpack_terminate(SOLVER_t<double> *solver, int verbose, bool debug);
extern int umfpack_terminate(SOLVER_t< complex<double> > *solver, int verbose, bool debug);



extern int hips_factorize(SOLVER_t<double> *solver, int verbose, bool debug);
extern int hips_factorize(SOLVER_t< complex<double> > *solver, int verbose, bool debug);

extern int hips_solve(SOLVER_t<double> *solver, double *RHS, double *x,int transposed, int verbose, bool debug);
extern int hips_solve(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x,int transposed, int verbose, bool debug);

extern int hips_initialize(SOLVER_t<double> & solver, int verbose, bool debug);
extern int hips_initialize(SOLVER_t< complex<double> > & solver, int verbose, bool debug);
  
extern int hips_terminate(SOLVER_t<double> *solver, int verbose, bool debug);
extern int hips_terminate(SOLVER_t< complex<double> > *solver, int verbose, bool debug);



extern int lapack_factorize(SOLVER_t<double> *solver, int verbose, bool debug);
extern int lapack_factorize(SOLVER_t< complex<double> > *solver, int verbose, bool debug);

extern int lapack_solve(SOLVER_t<double> *solver, double *RHS, double *x,int transposed, int verbose, bool debug);
extern int lapack_solve(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x,int transposed, int verbose, bool debug);

extern int lapack_initialize(SOLVER_t<double> & solver, int verbose, bool debug);
extern int lapack_initialize(SOLVER_t< complex<double> > & solver, int verbose, bool debug);
  
extern int lapack_terminate(SOLVER_t<double> *solver, int verbose, bool debug);
extern int lapack_terminate(SOLVER_t< complex<double> > *solver, int verbose, bool debug);



extern int spdomestic_factorize(SOLVER_t<double> *solver, int verbose, bool debug);
extern int spdomestic_factorize(SOLVER_t< complex<double> > *solver, int verbose, bool debug);
  
extern int spdomestic_solve(SOLVER_t<double> *solver, double *RHS, double *x,int transposed, int verbose, bool debug);
extern int spdomestic_solve(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x,int transposed, int verbose, bool debug);

extern int spdomestic_initialize(SOLVER_t<double> & solver, int verbose, bool debug);
extern int spdomestic_initialize(SOLVER_t< complex<double> > & solver, int verbose, bool debug);
  
extern int spdomestic_terminate(SOLVER_t<double> *solver, int verbose, bool debug);
extern int spdomestic_terminate(SOLVER_t< complex<double> > *solver, int verbose, bool debug);



extern int petsc_factorize(SOLVER_t<double> *solver, int verbose, bool debug);
extern int petsc_factorize(SOLVER_t< complex<double> > *solver, int verbose, bool debug);
  
extern int petsc_solve(SOLVER_t<double> *solver, double *RHS, double *x,int transposed, int verbose, bool debug);
extern int petsc_solve(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x,int transposed, int verbose, bool debug);

extern int petsc_initialize(SOLVER_t<double> & solver, int verbose, bool debug);
extern int petsc_initialize(SOLVER_t< complex<double> > & solver, int verbose, bool debug);
  
extern int petsc_terminate(SOLVER_t<double> *solver, int verbose, bool debug);
extern int petsc_terminate(SOLVER_t< complex<double> > *solver, int verbose, bool debug);



extern int pastix_factorize_parallel_distributed(SOLVER_t<double> *solver, int verbose, bool debug);
extern int pastix_factorize_parallel_distributed(SOLVER_t< complex<double> > *solver, int verbose, bool debug);

extern int pastix_factorize_sequential(SOLVER_t<double> *solver, int verbose, bool debug);
extern int pastix_factorize_sequential(SOLVER_t< complex<double> > *solver, int verbose, bool debug);

extern int pastix_solve_sequential(SOLVER_t<double> *solver, double *RHS, double *x,int transposed, int verbose, bool debug);
extern int pastix_solve_sequential(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x,int transposed, int verbose, bool debug);

extern int pastix_solve_parallel(SOLVER_t<double> *solver, double *RHS, double *x,int transposed, int verbose, bool debug);
extern int pastix_solve_parallel(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x,int transposed, int verbose, bool debug);

extern int pastix_initialize(SOLVER_t<double> & solver, int verbose, bool debug);
extern int pastix_initialize(SOLVER_t< complex<double> > & solver, int verbose, bool debug);
  
extern int pastix_terminate(SOLVER_t<double> *solver, int verbose, bool debug);
extern int pastix_terminate(SOLVER_t< complex<double> > *solver, int verbose, bool debug);



extern int gmres_factorize(SOLVER_t<double> *solver,            int verbose, bool debug);
extern int gmres_factorize(SOLVER_t< complex<double> > *solver, int verbose, bool debug);

extern int gmres_solve(SOLVER_t<double> *solver, double *RHS, double *x,                    int transposed, int verbose, bool debug);
extern int gmres_solve(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x,int transposed, int verbose, bool debug);

extern int gmres_initialize(SOLVER_t<double> & solver,            int verbose, bool debug);
extern int gmres_initialize(SOLVER_t< complex<double> > & solver, int verbose, bool debug);
  
extern int gmres_terminate(SOLVER_t<double> *solver,            int verbose, bool debug);
extern int gmres_terminate(SOLVER_t< complex<double> > *solver, int verbose, bool debug);


extern int diagonal_factorize(SOLVER_t<double> *solver,            int verbose, bool debug);
extern int diagonal_factorize(SOLVER_t< complex<double> > *solver, int verbose, bool debug);

extern int diagonal_solve(SOLVER_t<double> *solver, double *RHS, double *x,                    int transposed, int verbose, bool debug);
extern int diagonal_solve(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x,int transposed, int verbose, bool debug);

extern int diagonal_initialize(SOLVER_t<double> & solver,            int verbose, bool debug);
extern int diagonal_initialize(SOLVER_t< complex<double> > & solver, int verbose, bool debug);
  
extern int diagonal_terminate(SOLVER_t<double> *solver,            int verbose, bool debug);
extern int diagonal_terminate(SOLVER_t< complex<double> > *solver, int verbose, bool debug);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  dynamical solver registration

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

extern  int solver_register_dynamic_d(const char *name,
                                     int (*factorize)  (SOLVER_t<double> *, int, bool),
                                     int (*solve)      (SOLVER_t<double> *, double *, double *, int, int, bool),
                                     int (*initialize) (SOLVER_t<double> &, int, bool),
                                     int (*terminate)  (SOLVER_t<double> *, int, bool),
                                     int (*error)      (int, int, bool),
                                     int verbose, bool debug);

extern  int solver_register_dynamic_z(const char *name,
                                     int (*factorize)  (SOLVER_t< complex<double> > *, int, bool),
                                     int (*solve)      (SOLVER_t< complex<double> > *, complex<double> *, complex<double> *, int, int, bool),
                                     int (*initialize) (SOLVER_t< complex<double> > &, int, bool),
                                     int (*terminate)  (SOLVER_t< complex<double> > *, int, bool),
                                     int (*error)      (int, int, bool),
                                     int verbose, bool debug);

extern  int solver_initialize_dynamic(SOLVER_t<double> & solver,            int id, int verbose, bool debug);
extern  int solver_initialize_dynamic(SOLVER_t< complex<double> > & solver, int id, int verbose, bool debug);

extern  int solver_identify_dynamic(SOLVER_t<double> & solver,            int id, int verbose, bool debug);
extern  int solver_identify_dynamic(SOLVER_t< complex<double> > & solver, int id, int verbose, bool debug);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  interface for external use

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

extern int solver_factorize(SOLVER_t<double> *solver           , packed_t<double> & matrix,  int verbose, bool debug);
extern int solver_factorize(SOLVER_t< complex<double> > *solver, packed_t< complex<double> > & matrix, int verbose, bool debug);

extern int solver_solve(SOLVER_t<double> *solver, double *RHS, double *x,        int transposed, int verbose, bool debug);
extern int solver_solve(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x,int transposed, int verbose, bool debug);

extern int solver_init(int solver_id, MPI_Comm world, int nprocs, int nunknowns, SOLVER_t<double>* & solver,            int verbose, bool debug);
extern int solver_init(int solver_id, MPI_Comm world, int nprocs, int nunknowns, SOLVER_t< complex<double> >* & solver, int verbose, bool debug);

extern int solver_terminate(SOLVER_t<double> *solver, int verbose, bool debug);
extern int solver_terminate(SOLVER_t< complex<double> > *solver, int verbose, bool debug);

extern int solver_terminate(SOLVER_t<double> *solver);
extern int solver_terminate(SOLVER_t< complex<double> > *solver);

extern int convert(SOLVER_t<double> & solver,            int targeted_packing, int targeted_ordering, int targeted_numbering, int verbose, bool debug);
extern int convert(SOLVER_t< complex<double> > & solver, int targeted_packing, int targeted_ordering, int targeted_numbering, int verbose, bool debug);

#endif
