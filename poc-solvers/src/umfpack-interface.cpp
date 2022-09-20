
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
#include <complex.h>

#include "poc-solvers.h"
#include "solvers-interface.h"

extern int nbsyslinearz2;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  wrappers necessary to template use
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#ifdef UMFPACK

int umfpack_symbolic(SOLVER_t<double> & solver, int n_row, int n_col, const int Ap[], const int Ai[], double Ax[], void **Symbolic, const double Control[UMFPACK_CONTROL], double Info[UMFPACK_INFO]) {
  int status;
  
  status = umfpack_di_symbolic (n_row, n_col, Ap, Ai, Ax, Symbolic, Control, Info) ; 
    
  return(status);
}

int umfpack_numeric(SOLVER_t<double> & solver, int n_row, int n_col, const int Ap[], const int Ai[], double Ax[], void *Symbolic, void **Numeric, const double Control[UMFPACK_CONTROL], double Info[UMFPACK_INFO]) {
  int status;

  status = umfpack_di_numeric (Ap, Ai, Ax, Symbolic, Numeric, Control, Info) ;
  
  return(status);
}

int umfpack_symbolic(SOLVER_t< complex<double> > & solver, int n_row, int n_col, const int Ap[], const int Ai[], const complex<double> Ax[], void **Symbolic, const double Control[UMFPACK_CONTROL], double Info[UMFPACK_INFO]) {
  int status;
  double  *Ar, *Az;
  size_t nnz=solver.packed->nnz;
  umfpack_t< complex<double> > *id;
  
  id = (umfpack_t< complex<double> > *) solver.parameters;
    
  Az = new double[nnz];
  Ar = new double[nnz];
  
  for(size_t bcl=0;bcl<nnz;bcl++) {
    Ar[bcl] = real(Ax[bcl]);
    Az[bcl] = imag(Ax[bcl]);
    }
    
  status = umfpack_zi_symbolic (n_row, n_col, Ap, Ai, Ar, Az, Symbolic, Control, Info) ; 
  
  id->matrix_real=Ar;
  id->matrix_imag=Az;
  
  return(status);
}

int umfpack_numeric(SOLVER_t< complex<double> > & solver, int n_row, int n_col, const int Ap[], const int Ai[], const complex<double> Ax[], void *Symbolic, void **Numeric, const double Control[UMFPACK_CONTROL], double Info[UMFPACK_INFO]) {
  int status;
  umfpack_t< complex<double> > *id;
  
  id = (umfpack_t< complex<double> > *) solver.parameters;

  status = umfpack_zi_numeric (Ap, Ai, id->matrix_real, id->matrix_imag, Symbolic, Numeric, Control, Info) ;
  
  return(status);
}

int umfpack_solve(SOLVER_t<double> & solver, int sys, const int Ap [ ], const int Ai [ ], const double Ax [ ], double RHS[ ], void *Numeric, const double Control [UMFPACK_CONTROL],  double Info [UMFPACK_INFO]) 
{
  int status;
  double  *tmp;
  int n=solver.packed->nrows;
  
  tmp = (double *) malloc (n * sizeof (double)) ;
  
  for(int i=0;i<n;i++) tmp[i]=RHS[i];

  status = umfpack_di_solve (sys, Ap, Ai, Ax, RHS, tmp, Numeric, Control, Info) ;

  free(tmp);
  
  return(status);
}

int umfpack_solve(SOLVER_t< complex<double> > & solver, int sys, const int Ap [ ], const int Ai [ ], const complex<double> Ax [ ], complex<double> RHS[ ], void *Numeric, const double Control [UMFPACK_CONTROL],  double Info [UMFPACK_INFO]) 
{
  int status;
  double  *tmpr, *tmpz, *rhsr, *rhsz;
  int n=solver.packed->nrows;
  umfpack_t< complex<double> > *id;
  
  id = (umfpack_t< complex<double> > *) solver.parameters;
  
  tmpr = new double[n];
  tmpz = new double[n];
  rhsr = new double[n];
  rhsz = new double[n];

  for(size_t bcl=0;bcl<n;bcl++) {
    rhsr[bcl] = real(RHS[bcl]);
    rhsz[bcl] = imag(RHS[bcl]);
    }
      
  status = umfpack_zi_solve (sys, Ap, Ai, id->matrix_real, id->matrix_imag, tmpr, tmpz, rhsr, rhsz, Numeric, Control, Info) ;
  
  for(size_t bcl=0;bcl<n;bcl++)  RHS[bcl] = tmpr[bcl] + tmpz[bcl] * _Complex_I; 

  ptr_delete(tmpr);
  ptr_delete(tmpz);
  ptr_delete(rhsr);
  ptr_delete(rhsz);
  
  return(status);
}

#endif

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  umfpack interface
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int umfpack_factorize_template(SOLVER_t<T> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
#ifdef UMFPACK
  umfpack_t<T> *parameters;
  char *solver_name;
  double *Info, *Control;
  T *Ax;
  int *Ap, *Ai,n,nnz,bcl;
  void *Symbolic, *Numeric ;
  packed_t<T>  *packed;
  
  solver_name=solver->name;
  if(verbose==1) printf("factorize UMFPACK ...\n");
  
  parameters=(umfpack_t<T> *) solver->parameters;
  
  packed=solver->packed;
  
  Symbolic=parameters->Symbolic;
  Numeric =parameters->Numeric;
  Info    =parameters->Info;
  Control =parameters->Control;
  
  n=packed->nrows;
  nnz=packed->nnz;
  
  if(solver->packed->pointer32==0) {
    bool duplicate=false;
    status=recast(solver->packed->pointer, solver->packed->pointer32, duplicate, n+1);
    packed->pointer32=solver->packed->pointer32;
    }
  
  Ap=packed->pointer32;
  Ai=packed->incidence;
  Ax=packed->x;
  
  status=umfpack_symbolic(*solver, n, n, Ap, Ai, Ax, &Symbolic, Control, Info) ; 
  if (status < 0) {
    printf("%s : status=%d\n",__func__,status);
    umfpack_zi_report_info (Control, Info) ;
    umfpack_zi_report_status (Control, status) ;
    return(-1);
    }
  parameters->Symbolic=Symbolic;
    
  status=umfpack_numeric(*solver, n, n, Ap, Ai, Ax, Symbolic, &Numeric, Control, Info) ; 
  if (status < 0) {
    printf("%s : status=%d\n",__func__,status);
    umfpack_zi_report_info (Control, Info) ;
    umfpack_zi_report_status (Control, status) ;
    return(-1);
    }
  parameters->Numeric=Numeric;
    
  return(status);
#else
  status=-1;
  printf("Please compile poc-solvers library with -DUMFPACK \n");
  return(status);
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int umfpack_factorize(SOLVER_t<double> *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=umfpack_factorize_template(solver, verbose, debug);
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int umfpack_factorize(SOLVER_t< complex<double> > *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=umfpack_factorize_template(solver, verbose, debug);
  return(status); 
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int umfpack_solve_template(SOLVER_t<T> *solver, T *RHS, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
#ifdef UMFPACK
  umfpack_t<T> *parameters;
  char *solver_name;
  double *Info, *Control, *tmpr, *tmpz;
  T *Ax;
  int *Ap, *Ai,n,nnz;
  void *Symbolic, *Numeric ;
  packed_t<T> *packed;
 
  status = 0;
  solver_name=solver->name;

  parameters = (umfpack_t<T> *) solver->parameters;
  
  Symbolic=parameters->Symbolic;
  Numeric=parameters->Numeric;
  Info=parameters->Info;
  Control=parameters->Control;
  
  packed=(packed_t<T> *) solver->packed;
  
  Ap=packed->pointer32;
  Ai=packed->incidence;
  Ax=packed->x;
  
  n=packed->nrows;
  nnz=packed->nnz;

  if (transposed==1) 
    status=umfpack_solve (*solver, UMFPACK_At, Ap, Ai, Ax, RHS, Numeric, Control, Info) ;
  else 
    status=umfpack_solve (*solver, UMFPACK_A, Ap,  Ai, Ax, RHS, Numeric, Control, Info) ;
    
  if (status < 0) {
    umfpack_zi_report_info (Control, Info) ;
    umfpack_zi_report_status (Control, status) ;
    }
#endif
    
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int umfpack_solve(SOLVER_t<double> *solver, double *RHS, double *x, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=umfpack_solve_template(solver, RHS, transposed, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int umfpack_solve(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=umfpack_solve_template(solver, RHS, transposed, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int umfpack_initialize_template(SOLVER_t<T> & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  MPI_Comm world;
  int rank;
  
#ifdef UMFPACK
  umfpack_t<T> *parameters=new umfpack_t<T>;
  solver.parameters = parameters;
  umfpack_di_defaults (parameters->Control);
/*     parameters->Control [UMFPACK_PRL] = 6 ; */
/*     parameters->Control [UMFPACK_PRL] = 5 ; */
/*     umfpack_di_report_control (parameters->Control) ; */

  solver.targeted_packing    = PACKED;
  solver.targeted_ordering   = CSC;
  solver.targeted_addressing = MATRIX_POINTER;
  solver.targeted_numbering  = 0;
  
  solver.initialize = umfpack_initialize;
  solver.terminate  = umfpack_terminate;
  solver.factorize  = umfpack_factorize;
  solver.solve      = umfpack_solve;
  
  solver.id=SOLVER_ID_UMFPACK;
  
  return(0);
#else 
    printf("Please compile poc-solvers library with -DUMFPACK \n");
    exit(-1);
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int umfpack_initialize(SOLVER_t<double> & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=umfpack_initialize_template(solver, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int umfpack_initialize(SOLVER_t< complex<double> > & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=umfpack_initialize_template(solver, verbose, debug);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int umfpack_terminate_template(SOLVER_t<T> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  MPI_Comm world;
  
#ifdef UMFPACK
  void *Symbolic, *Numeric;
  
  umfpack_t<T> *parameters = (umfpack_t<T> *) solver->parameters;
  
  Symbolic = parameters->Symbolic;
  Numeric  = parameters->Numeric;
  umfpack_di_free_symbolic (&Symbolic);
  umfpack_di_free_numeric  (&Numeric);
  
  ptr_delete(parameters->matrix_real);
  ptr_delete(parameters->matrix_imag);
  
  status=0;
  
  return(status);
#else
  return(status);
#endif 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int umfpack_terminate(SOLVER_t<double> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=umfpack_terminate_template(solver, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int umfpack_terminate(SOLVER_t< complex<double> > *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=umfpack_terminate_template(solver, verbose, debug);
  
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

  int factorize_umfpack(solver_t *slv, triplet *Mat, int verbose)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t<double> solver;

  bool debug=false;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);

  slv->Mat = Mat; 
      
  solver.d_import2packed(slv);

  status=umfpack_factorize(&solver, verbose, debug);
  
  solver.d_export_packed(slv);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorizez_umfpack(solver_t *slv, tripletz *Mat, int verbose)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t< complex<double> > solver;
 
  bool debug=false;

  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  slv->Mat = Mat;
  
  solver.z_import2packed(slv);

  status=umfpack_factorize(&solver, verbose, debug);
  
  solver.z_export_packed(slv);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solve_umfpack(solver_t *slv, double *RHS, int transposed)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int verbose=0;
  bool debug=false;
  
  SOLVER_t<double> solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  solver.d_import2packed(slv);

  status=umfpack_solve(&solver, RHS, 0, transposed, verbose, debug);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solvez_umfpack(solver_t *slv, complex<double> *RHS, int transposed)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t< complex<double> > solver;
 
  int verbose=0;
  bool debug=false;
  
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  solver.z_import2packed(slv);

  status=umfpack_solve(&solver, RHS, 0, transposed, verbose, debug);
  
  return(status); 
}

#endif
