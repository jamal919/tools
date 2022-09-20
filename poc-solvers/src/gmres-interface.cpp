
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
#include "linear-gmres.h"
#include "poc-solvers.h"
#include "solvers-interface.h"

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  gmres interface
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int gmres_factorize_template(SOLVER_t<T> *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
#ifdef HAVE_MPI
  gmres_t<T> *parameters;
  return(status);  
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gmres_factorize(SOLVER_t<double> *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=gmres_factorize_template(solver, verbose, debug);
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gmres_factorize(SOLVER_t< complex<double> > *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=gmres_factorize_template(solver, verbose, debug);
  return(status); 
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int gmres_solve_template(SOLVER_t<T> *solver, T *RHS, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
#ifdef HAVE_MPI
  gmres_t<T> *parameters;
  return(status); 

#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gmres_solve(SOLVER_t<double> *solver, double *RHS, double *x, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=gmres_solve_template(solver, RHS, transposed, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gmres_solve(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=gmres_solve_template(solver, RHS, transposed, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int gmres_initialize_template(SOLVER_t<T> & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  MPI_Comm world;
  int rank;
  
  gmres_t<T> *parameters=new gmres_t<T>;
  solver.parameters = parameters;

  solver.targeted_packing    = PACKED;
  solver.targeted_ordering   = CSR;
  solver.targeted_addressing = MATRIX_POINTER;
  solver.targeted_numbering  = 0;
  
  solver.initialize = gmres_initialize;
  solver.terminate  = gmres_terminate;
  solver.factorize  = gmres_factorize;
  solver.solve      = gmres_solve;
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gmres_initialize(SOLVER_t<double> & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=gmres_initialize_template(solver, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gmres_initialize(SOLVER_t< complex<double> > & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=gmres_initialize_template(solver, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int gmres_terminate_template(SOLVER_t<T> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  MPI_Comm world;
  
  gmres_t<T> *parameters = (gmres_t<T> *) solver->parameters;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gmres_terminate(SOLVER_t<double> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=gmres_terminate_template(solver, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gmres_terminate(SOLVER_t< complex<double> > *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=gmres_terminate_template(solver, verbose, debug);
  
  return(status);
}


#if OLD_SOLVER_VERSION
// if activated, use historical (1.0) code (still available in the c-code branch)
// not maintained, kept only for testing purposes, aimed to be removed
#else
// if not, provide historical (1.0) interface to new code (2.0)
// designed for smooth transition purposes, aimed to be removed

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorize_gmres(solver_t *slv, triplet *Mat, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t<double> solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);

  slv->Mat = Mat; 
      
  solver.d_import2packed(slv);

  status=gmres_factorize(&solver, verbose, debug);
  
  solver.d_export_packed(slv);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorizez_gmres(solver_t *slv, tripletz *Mat, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t< complex<double> > solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  slv->Mat = Mat;
  
  solver.z_import2packed(slv);

  status=gmres_factorize(&solver, verbose, debug);
  
  solver.z_export_packed(slv);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solve_gmres(solver_t *slv, double *RHS, int transposed, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int verbose=0;
  
  SOLVER_t<double> solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  solver.d_import2packed(slv);

  status=gmres_solve(&solver, RHS, 0, transposed, verbose, debug);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solvez_gmres(solver_t *slv, complex<double> *RHS, int transposed, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t< complex<double> > solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  solver.z_import2packed(slv);

  status=gmres_solve(&solver, RHS, 0, transposed, verbose, debug);
  
  return(status); 
}

#endif


