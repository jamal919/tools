
/*******************************************************************************

  poc_solvers, SIROCCO 2006-2018

  Unstructured Ocean Grid initiative

Developers:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
Contributors:
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France

*******************************************************************************/

#include <string.h>
#include <omp.h>
#include "poc-solvers.h"
#include "solvers-interface.h"

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  wrappers necessary to template use
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

cxsparse_cs_t<double> * cs_sqr(int order, const cxsparse_cs_t<double> *A, int qr) {
  cxsparse_cs_t<double> *S;
  S=(cxsparse_cs_t<double> *) cs_di_sqr(order, (const cs_di *) A, qr);
  return(S);
  }
  
cxsparse_cs_t< complex<double> > * cs_sqr(int order, const cxsparse_cs_t< complex<double> > *A, int qr) {
  cxsparse_cs_t< complex<double> > *S;
  S=(cxsparse_cs_t< complex<double> > *) cs_ci_sqr(order, (const cs_ci *) A, qr);
  return(S);
  }
  
cxsparse_cs_t<double> * cs_lu(const cxsparse_cs_t<double> *A, const cxsparse_symbolic_t<double> *S, double tol) {
  cxsparse_cs_t<double> *N;
  N=(cxsparse_cs_t<double> *) cs_di_lu((const cs_di *) A, (const cs_dis *) S, tol);
  return(N);
  }
    
cxsparse_cs_t< complex<double> > * cs_lu(const cxsparse_cs_t< complex<double> > *A, const cxsparse_symbolic_t< complex<double> > *S, double tol) {
  cxsparse_cs_t< complex<double> > *N;
  N=(cxsparse_cs_t< complex<double> > *) cs_ci_lu((const cs_ci *) A, (const cs_cis *) S, tol);
  return(N);
  }
  
int cs_lsolve(const cxsparse_cs_t<double> *L, double *x) {
  cs_di_lsolve ((cs_di *) L, x) ;
  }
  
int cs_lsolve(const cxsparse_cs_t< complex<double> > *L, complex<double> *x) {
  cs_ci_lsolve ((cs_ci *) L, x) ;
  }
  
int cs_ltsolve(const cxsparse_cs_t<double> *L, double *x) {
  cs_di_ltsolve ((cs_di *) L, x) ;
  }
  
int cs_ltsolve(const cxsparse_cs_t< complex<double> > *L, complex<double> *x) {
  cs_ci_ltsolve ((cs_ci *) L, x) ;
  }
  
int cs_usolve(const cxsparse_cs_t<double> *U, double *x) {
  cs_di_usolve ((cs_di *) U, x) ;
  }

int cs_usolve(const cxsparse_cs_t< complex<double> > *U, complex<double> *x) {
  cs_ci_usolve ((cs_ci *) U, x) ;
  }

int cs_utsolve(const cxsparse_cs_t<double> *U, double *x) {
  cs_di_utsolve ((cs_di *) U, x) ;
  }

int cs_utsolve(const cxsparse_cs_t< complex<double> > *U, complex<double> *x) {
  cs_ci_utsolve ((cs_ci *) U, x) ;
  }

void cs_ipvec_interface(const int *p, const double *b, double *x, int n) {
  cs_di_ipvec (p, b, x, n);
}

void cs_ipvec_interface(const int *p, const complex<double> *b, complex<double> *x, int n) {
  cs_ci_ipvec (p, b, x, n);
}

int cs_sfree(cxsparse_symbolic_t<double> *S) {
  cs_sfree(S) ;
  }

int cs_sfree(cxsparse_symbolic_t< complex<double> > *S) {
  cs_sfree(S) ;
  }

int cs_nfree(cxsparse_numeric_t<double> *N) {
  cs_nfree(N) ;
  }

int cs_nfree(cxsparse_numeric_t< complex<double> > *N) {
  cs_nfree(N) ;
  }


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  spdomestic interface
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int spdomestic_factorize_template(SOLVER_t<T> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0,n,order=3,bcl;
  double tol=1.0; /* unsymetrique usage*/
  cxsparse_t<T> *parameters;
  packed_t<T>  *packed;
  
  packed=solver->packed;
  
  cxsparse_cs_t<T>       *A;
  cxsparse_symbolic_t<T> *S ;
  cxsparse_numeric_t <T> *N ;
  
  parameters = (cxsparse_t<T> *) solver->parameters;
  
  A = parameters->matrix;
  S = parameters->symbolic;
  N = parameters->numeric;  
  
  A->nzmax = packed->nnz;
  A->m     = packed->nrows;
  A->n     = packed->ncols;
  A->p     = (int*) packed->pointer;
  A->i     = (int*) packed->incidence;
  A->x     = packed->x;
  A->nz    = -1;
  
  S = (cxsparse_symbolic_t<T> *) cs_sqr (order, A, 0) ;      /* ordering and symbolic analysis */

  N = (cxsparse_numeric_t <T> *) cs_lu (A, S, tol) ;         /* numeric LU factorization       */
  
  parameters->symbolic = S;
  parameters->numeric  = N;
  
  solver->packed = packed;
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int spdomestic_factorize(SOLVER_t<double> *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=spdomestic_factorize_template(solver, verbose, debug);
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int spdomestic_factorize(SOLVER_t< complex<double> > *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=spdomestic_factorize_template(solver, verbose, debug);
  return(status); 
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int spdomestic_solve_template(SOLVER_t<T> *solver, T *RHS, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0,n;
  cxsparse_t<T> *parameters;

  cxsparse_cs_t<T>       *A;
  cxsparse_symbolic_t<T> *S ;
  cxsparse_numeric_t <T> *N ;
  T *x;

  parameters = (cxsparse_t<T> *) solver->parameters;

  A = parameters->matrix;
  S = parameters->symbolic;
  N = parameters->numeric;
  n = A->n ;
  
  x = (T *) cs_ci_malloc (n, sizeof (T)) ;      /* get workspace */

  cs_ipvec_interface (N->pinv, RHS, x, n) ;            /* x = b(p) */
  
  if (transposed) {
    cs_utsolve (N->U, x) ;                   /* x = U\x */
    cs_ltsolve (N->L, x) ;                   /* x = L\x */
    }
  else {
    cs_lsolve (N->L, x) ;                    /* x = L\x */
    cs_usolve (N->U, x) ;                    /* x = U\x */
    }

  cs_ipvec_interface (S->q, x, RHS, n) ;               /* b(q) = x */
  cs_ci_free (x) ;

  return(status); 

}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int spdomestic_solve(SOLVER_t<double> *solver, double *RHS, double *x,int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=spdomestic_solve_template(solver, RHS, transposed, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int spdomestic_solve(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x,int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=spdomestic_solve_template(solver, RHS, transposed, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int spdomestic_initialize_template(SOLVER_t<T> & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  MPI_Comm world;
  int rank;
    
  cxsparse_t<T> *parameters=new cxsparse_t<T>;
  
  solver.parameters=parameters;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  initialisation might be moved to factorization
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  parameters->matrix=new cxsparse_cs_t<T>;
  parameters->symbolic=new cxsparse_symbolic_t<T>;
  parameters->numeric=new cxsparse_numeric_t<T>;
     
  solver.targeted_packing    = PACKED;
  solver.targeted_ordering   = CSC;
  solver.targeted_addressing = MATRIX_POINTER;
  solver.targeted_numbering  = 0;

  solver.initialize = spdomestic_initialize;
  solver.terminate  = spdomestic_terminate;
  solver.factorize  = spdomestic_factorize;
  solver.solve      = spdomestic_solve;
  
  solver.id=SOLVER_ID_SpDOMESTIC;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int spdomestic_initialize(SOLVER_t<double> & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=spdomestic_initialize_template(solver, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int spdomestic_initialize(SOLVER_t< complex<double> > & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=spdomestic_initialize_template(solver, verbose, debug);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int spdomestic_terminate_template(SOLVER_t<T> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  cxsparse_t<T> *parameters;

  cxsparse_cs_t<T>       *A;
  cxsparse_symbolic_t<T> *S ;
  cxsparse_numeric_t <T> *N ;
  T *x;

  parameters = (cxsparse_t<T> *) solver->parameters;
  
  A = parameters->matrix;
  S = parameters->symbolic;
  N = parameters->numeric;
  
  cs_sfree (S) ;
  cs_nfree (N) ;
      
  return(status);
  
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int spdomestic_terminate(SOLVER_t<double> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=spdomestic_terminate_template(solver, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int spdomestic_terminate(SOLVER_t< complex<double> > *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=spdomestic_terminate_template(solver, verbose, debug);
  
  return(status);
}



/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  solver interface compatibility with older versions
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#if OLD_SOLVER_VERSION
// if activated, use historical (1.0) code (still available in the c-code branch)
// not maintained, kept only for testing purposes, aimed to be removed
#else
// if not, provide historical (1.0) interface to new code (2.0)
// designed for smooth transition purposes, aimed to be removed

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorize_spdomestic(solver_t *slv, triplet *Mat)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t<double> solver;
  
  int verbose=0;
  bool debug=false;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);

  slv->Mat = Mat; 
      
  solver.d_import2packed(slv);

  status=spdomestic_factorize(&solver, verbose, debug);
  
  solver.d_export_packed(slv);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorizez_spdomestic(solver_t *slv, tripletz *Mat)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t< complex<double> > solver;
 
  int verbose=0;
  bool debug=false;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  slv->Mat = Mat;
  
  solver.z_import2packed(slv);

  status=spdomestic_factorize(&solver, verbose, debug);
  
  solver.z_export_packed(slv);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solve_spdomestic(solver_t *slv, double *RHS, int transposed)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int verbose=0;
  bool debug=false;
  
  SOLVER_t<double> solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  solver.d_import2packed(slv);

  status=spdomestic_solve(&solver, RHS, 0, transposed, verbose, debug);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solvez_spdomestic(solver_t *slv, complex<double> *RHS, int transposed)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int verbose=0;
  bool debug=false;
  
  SOLVER_t< complex<double> > solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  solver.z_import2packed(slv);

  status=spdomestic_solve(&solver, RHS, 0, transposed, verbose, debug);
  
  return(status); 
}
#endif

