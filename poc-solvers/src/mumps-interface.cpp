
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

#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match MUMPS documentation */

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  see http://mumps.enseeiht.fr/doc/userguide_5.2.0.pdf
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  communicator strategy
  
  TODO : needs further testing and improvement to allow sequential solving in
         MPI T-UGOm
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#define USE_DEFAULT_COMMUNICATOR 1

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  wrappers necessary to template use
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#ifdef MUMPS

  void mumps_c(mumps_t<double> *id) {
    dmumps_c((DMUMPS_STRUC_C *) id);
  }

  void mumps_c(mumps_t< complex<double> > *id) {
    zmumps_c((ZMUMPS_STRUC_C *) id);
  }

#endif

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  mumps interface
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int mumps_factorize_template(SOLVER_t<T> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#ifdef MUMPS
  mumps_t<T> *parameters;
  int myid, ierr, nglob=0, count, rank, size;
  char *solver_name;
  triplet_t<T> *Mat;
  MPI_Comm communicator;
  
#if USE_DEFAULT_COMMUNICATOR
#else
/*------------------------------------------------------------------------------
  handle reduced cpu set case (small system do) */
  if(solver->rank==-1) return(0);
#endif
  
  communicator=solver->communicator;
  
  Mat=solver->COOtriplet;
  
#ifdef HAVE_MPI
  ierr = MPI_Comm_rank(communicator, &rank);
  ierr = MPI_Comm_size(communicator, &size);
#else
  ierr = 0;
#endif
  
  parameters = (mumps_t<T> *) solver->parameters;
  (*parameters).job=4;
  
/*------------------------------------------------------------------------------
  get global row number */
  count=1;
#ifdef HAVE_MPI
  ierr = MPI_Allreduce (  &(Mat->nrows), &nglob, count, MPI_INT, MPI_SUM, communicator );
  (*parameters).n = nglob;
#else
  (*parameters).n = Mat->nrows;
  ierr = 0;
#endif
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

ICNTL(1) is the output stream for error messages
ICNTL(2) is the output stream for diagnostic printing, statistics, and warning message
ICNTL(3) is the output stream for global information, collected on the host
ICNTL(4) is the level of printing for error, warning, and diagnostic messages
ICNTL(5) controls the matrix input format
ICNTL(6) permutes the matrix to a zero-free diagonal and/or scale the matrix
ICNTL(7) computes a symmetric permutation in case of sequential analysis
ICNTL(8) describes the scaling strategy
ICNTL(9) computes the solution using A or A T
ICNTL(10) applies the iterative refinement to the computed solution
ICNTL(11) computes statistics related to an error analysis
ICNTL(12) defines an ordering strategy for symmetric matrices
ICNTL(13) controls the parallelism of the root node
ICNTL(14) controls the percentage increase in the estimated working space
ICNTL(18) defines the strategy for the distributed input matrix
ICNTL(19) computes the Schur complement matrix
ICNTL(20) determines the format (dense or sparse) of the right-hand sides
ICNTL(21) determines the distribution (centralized or distributed) of the solution vectors
ICNTL(22) controls the in-core/out-of-core (OOC) factorization and solve
ICNTL(23) corresponds to the maximum size of the working memory in MegaBytes that MUMPS can allocate per working processor
ICNTL(24) controls the detection of “null pivot rows”
ICNTL(25) allows the computation of a solution of a deficient matrix and also of a null space basis
ICNTL(26) drives the solution phase if a Schur complement matrix
ICNTL(27) controls the blocking size for multiple right-hand sides
ICNTL(28) determines whether a sequential or parallel computation of the ordering is performed
ICNTL(29) defines the parallel ordering tool to be used to compute the fill-in reducing permutation
ICNTL(30) computes a user-specified set of entries in the inverse A−1 of the original matrix
ICNTL(31) indicates which factors may be discarded during the factorizatio
ICNTL(32) performs the forward elimination of the right-hand sides during the factorization
ICNTL(33) computes the determinant of the input matrix
ICNTL(35) controls the activation of Block Lock Rank (BLR) based factorization

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  (*parameters).ICNTL(18)=3;
  
/*------------------------------------------------------------------------------
  prior to 5.1.1 */
//     (*parameters).nz_loc = Mat->nnz;

/*------------------------------------------------------------------------------
  from 5.1.1 and later */
  (*parameters).nnz_loc = Mat->nnz;
  
  (*parameters).irn_loc= (int *) Mat->i;  
  (*parameters).jcn_loc= (int *) Mat->j; 
  (*parameters).a_loc  = (T *)   Mat->x;  

  (*parameters).ICNTL(1) =-1; 
  (*parameters).ICNTL(2) =-1; 
  (*parameters).ICNTL(3) =-1; 
  (*parameters).ICNTL(4) = 4;
  (*parameters).ICNTL(14)=20; //corresponds to the percentage increase in the estimated working space
//   (*parameters).ICNTL(14)=50; //corresponds to the percentage increase in the estimated working space
  
//   (*parameters).ICNTL(1)=6;
//   (*parameters).ICNTL(2)=0;
//   (*parameters).ICNTL(3)=6;
  
/*------------------------------------------------------------------------------
  dense RHS  */
  (*parameters).ICNTL(20)=0;
  
/*------------------------------------------------------------------------------
  centralized solution  */
  (*parameters).ICNTL(21)=0;

  mumps_c(parameters);
  
  ierr=(*parameters).info[0];
  
  if(ierr==-1)  printf("cpu %d: MUMPS factorisation failed, status=%d due to cpu %d\n",rank, ierr, parameters->info[1]);
  if(ierr==-9)  printf("cpu %d: MUMPS factorisation failed, status=%d due to memory deficiency %d\n",rank, ierr, parameters->info[1]);
  if(ierr==-10) printf("cpu %d: MUMPS factorisation failed, status=%d due to singular matrix\n",rank, ierr);
  if(ierr==-13) printf("cpu %d: MUMPS factorisation failed, status=%d due to memory deficience %d\n",rank, ierr, parameters->info[2]);
  if(ierr==0) {
    if(verbose==1) {
      printf("%s, cpu=%d: MEMORY USE, %d used by cpu, %d total Mo, status=%d\n",__FUNCTION__,rank, (*parameters).infog[20],(*parameters).infog[21], ierr);
      }
    if(verbose==1 && rank==0) {
      printf("MEMORY SCALING DIAG %d %d\n",size,(*parameters).infog[21]);
      }
    }
    
  return(ierr);
  
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mumps_factorize(SOLVER_t<double> *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=mumps_factorize_template(solver, verbose, debug);
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mumps_factorize(SOLVER_t< complex<double> > *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=mumps_factorize_template(solver, verbose, debug);
  return(status); 
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int mumps_solve_template(SOLVER_t<T> *solver, T *RHS, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#ifdef MUMPS
  mumps_t<T> *parameters;
  char *solver_name;
  int ierr, rank;
  double cput1, cput2, cput3;
  int status;
  MPI_Comm communicator;
  
#if USE_DEFAULT_COMMUNICATOR
#else
/*------------------------------------------------------------------------------
  handle reduced cpu set case (small system do) */
  if(solver->rank==-1) return(0);
#endif
  
  communicator=solver->communicator;
  
#ifdef HAVE_MPI
  cput1 = MPI_Wtime();
#endif
  
  ierr = MPI_Comm_rank(communicator, &rank);
  solver_name=solver->name;

  parameters = (mumps_t<T> *) solver->parameters;
  
  if (transposed == 1) {
    (*parameters).ICNTL(9)=2;
    }
    
  if (rank == 0) {
    (*parameters).rhs = (T *) RHS;
    }
  (*parameters).job=3;
  
  mumps_c(parameters);

#ifdef HAVE_MPI
  cput2 = MPI_Wtime(); 
  cput3 = cput2-cput1;
  if (rank == 0) {
    if(verbose==1) printf(" MUMPS resolution time =%lf\n",cput3);
    }
#endif

  ierr=(*parameters).info[0];
  if(ierr==-1) {
    printf("cpu %d: MUMPS solver failed, status=%d due to cpu %d\n",rank, ierr, parameters->info[1]);
    }
  if(ierr==-22) {
    printf("cpu %d: MUMPS solver failed, status=%d array fault %d\n",rank, ierr, parameters->info[1]);
    }
  return(ierr);
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mumps_solve(SOLVER_t<double> *solver, double *RHS, double *x, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=mumps_solve_template(solver, RHS, transposed, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mumps_solve(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=mumps_solve_template(solver, RHS, transposed, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int mumps_initialize_template(SOLVER_t<T> & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  MPI_Comm world;
  int rank;
  MPI_Comm communicator;
  
#if USE_DEFAULT_COMMUNICATOR
  if(solver.rank!=-1 and solver.communicator!=MPI_COMM_WORLD) {
    status=MPI_Comm_free(&solver.communicator);
    }
  solver.communicator=MPI_COMM_WORLD; // TESTING
#endif

  communicator=solver.communicator;
  
#ifdef MUMPS
  mumps_t<T> *parameters=new mumps_t<T>;
  
  solver.parameters=parameters;
  
  status = MPI_Comm_rank(communicator, &rank);
  
  parameters->job=JOB_INIT; 
  parameters->par=1; 
  parameters->sym=0;
#if USE_DEFAULT_COMMUNICATOR
  parameters->comm_fortran=USE_COMM_WORLD;
#else
  parameters->comm_fortran=(MUMPS_INT) MPI_Comm_c2f(communicator);
#endif
  
  mumps_c(parameters);
  
  if(verbose==1 && rank==0) printf("dmumps_par->instance_number %d \n",parameters->instance_number);
  if(verbose==1) printf("init MUMPS solver (cpu=%d), finish \n",rank);
    
  solver.targeted_packing    = PACKED;
  solver.targeted_ordering   = MATRIX_ROW_MAJOR;
  solver.targeted_addressing = MATRIX_TRIPLET;
  solver.targeted_numbering  = MATRIX_F_NUMBERING;
  
  solver.initialize = mumps_initialize;
  solver.terminate  = mumps_terminate;
  solver.factorize  = mumps_factorize;
  solver.solve      = mumps_solve;
  
  solver.id=SOLVER_ID_MUMPS;
  
  return(status);
#else 
    printf("Please compile poc-solvers library with -DMUMPS \n");
    exit(-1);
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mumps_initialize(SOLVER_t<double> & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=mumps_initialize_template(solver, verbose, debug);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mumps_initialize(SOLVER_t< complex<double> > & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=mumps_initialize_template(solver, verbose, debug);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int mumps_terminate_template(SOLVER_t<T> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  MPI_Comm world;
  
#ifdef MUMPS

  mumps_t<T> *parameters = (mumps_t<T> *) solver->parameters;
  (*parameters).job=JOB_END;
  mumps_c(parameters); 
  
  status=(*parameters).info[0];
  
/*------------------------------------------------------------------------------
  delete done in solver destroy() */
//   delete solver->parameters;
//   solver->parameters=0;
    
  return(status);
#else
  return(status);
#endif 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mumps_terminate(SOLVER_t<double> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=mumps_terminate_template(solver, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mumps_terminate(SOLVER_t< complex<double> > *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=mumps_terminate_template(solver, verbose, debug);
  
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

  int factorize_mumps(solver_t *slv, triplet *Mat, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t<double> solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);

  slv->Mat = Mat; 
      
  solver.d_import2packed(slv);

  status=mumps_factorize(&solver, verbose, debug);
  
  solver.d_export_packed(slv);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorizez_mumps(solver_t *slv, tripletz *Mat, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t< complex<double> > solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  slv->Mat = Mat;
  
  solver.z_import2packed(slv);

  status=mumps_factorize(&solver, verbose, debug);
  
  solver.z_export_packed(slv);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solve_mumps(solver_t *slv, double *RHS, int transposed, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  SOLVER_t<double> solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  solver.d_import2packed(slv);

  status=mumps_solve(&solver, RHS, 0, transposed, verbose, debug);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solvez_mumps(solver_t *slv, complex<double> *RHS, int transposed, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t< complex<double> > solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  solver.z_import2packed(slv);

  status=mumps_solve(&solver, RHS, 0, transposed, verbose, debug);
  
  return(status); 
}

#endif
