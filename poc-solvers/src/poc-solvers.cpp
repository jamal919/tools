
/**************************************************************************

  POC-SOLVERS interface, 2006-2019

  Part of the SIROCCO national service (INSU, France)

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France

***************************************************************************/

#define VERBOSE

#include <string.h>
#include <omp.h>
#include <complex.h>

#include "solvers-interface.h"
#include "linear-gmres.h"
#include "solverlib.h"
#include "mpi-reducers.h"
#include "perftimer.h"

extern int MPI_Init_thread_DONE;

static bool CheckSynchronicity=false;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  triplet : i,j,value (row-arranged) or j,i,value (column-arranged)
  
  CSR : value array, column index array, pointer array (optionally cardinal array)
  CSC : value array,    row index array, pointer array (optionally cardinal array)

             natural ordering    natural format    natural numbering  
              
  PASTIX   compressed columns     col-major CSC              FORTRAN
  UMFPACK  compressed columns     col-major CSC                    C
  CXSPARSE compressed columns     col-major CSC                    C
  MUMPS       compressed rows     row-major COO/CSR          FORTRAN
  HIPS        compressed rows     row-major COO/CSR        C/FORTRAN (HIPS_FORTRAN_NUMBERING switch)
  PETSC       compressed rows     row-major CSR                    C?
  MAPHYS      


  T-UGOm natural format is CSR (compressed rows ordering using pointers to 
  column indices for each row)
  
  Issues:
  
  - MPI distributed pastix kept in CSC, then transpose system solved
  
  - SCOTCH found to be sensitive to partitions contiguousity (in pastix)
 


xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            assignable communicator
              
  PASTIX                        YES
  UMFPACK        do not support MPI
  CXSPARSE       do not support MPI
  MUMPS                         YES
  HIPS                          YES
  PETSC                         ???
  MAPHYS                        ???      


xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  Some OpenMPI distributions do not support MPI_THREAD_MULTIPLE, that is wanted
  by default by SCOTCH when invoked by PASTIX. You can check that by launching:
  
  >ompi_info | grep -i thread
  
  In that case, you have 2 options :
   - change your MPI library for one supporting MPI_THREAD_MULTIPLE 
   - Rebuild Scotch without SCOTCH_PTHREAD and PaStiX with FORCE_NOSMP

 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  We have encountered a SMP issue with umfpack (requesting libopenblas_pthreads)
  libumfpack-5.7.6.so version, detected in August 2017

  issue : multi-threading disabled
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  MPI RHS/SOLUTION : distributed, centralized (on master proc), on demand (distributed/centralized)

             MPI             RHS      SOLUTION  
             
  PASTIX     yes     distributed   distributed
  UMFPACK     no               -             -
  CXSPARSE    no               -             -
  MUMPS      yes     centralized     on demand 
  HIPS       yes       on demand     on demand
  MAPHYS     yes       on demand     on demand

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class T> int factorize_template(SOLVER_t<T> *solver, packed_t<T> & matrix, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  factorize linear system matrix
  
  (optionally) operate first a format conversion to fit solver requirements
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  char *solver_name;
  int status;
  double default_precision=1.e-08;
  double cputime1,cputime2;
  double ellapsed, max, average;
  bool show_perf=(solver->options[OPTIONS_SHOW_FACTORIZE_PERF]==1);
  int target_format, target_ordering, target_numbering;
  
  bool track_call=false;

  solver_name=solver->name;
  
  if(track_call) printf("%s cpu=%d : ---> 1\n",__func__,solver->rank);
  if(CheckSynchronicity) status=mpi_synchronicity(solver->world_communicator, solver->world_rank, solver->world_size, __FUNCTION__ , __LINE__, false);
    
//   int proc_id,ncpus,ierr;
//   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
//   ierr = MPI_Comm_size(MPI_COMM_WORLD, &ncpus);
  
  if(solver==0) TRAP_ERR_EXIT(-1,"solver not initialized; terminate job");
  
  if(show_perf) cputime1 = MPI_Wtime();
  
  if(track_call) printf("%s cpu=%d : ---> 2\n",__func__,solver->rank);
  
/*------------------------------------------------------------------------------
  (future) format conversion; presently non-operating routine */
  status=convert(*solver, solver->targeted_packing, solver->targeted_ordering, solver->targeted_numbering, verbose, debug);
  if(status!=0) return(-1);
    
  if(track_call) printf("%s cpu=%d : ---> 3 status=%d\n",__func__,solver->rank,status);
  
/*------------------------------------------------------------------------------
  factorizing */
  if(solver->factorize==0) TRAP_ERR_EXIT(-1,"solver %s factorize function not initialized; terminate job",solver_name);
  status=solver->factorize(solver, verbose, debug);
    
  if(track_call) printf("%s cpu=%d : ---> 4 status=%d\n",__func__,solver->rank,status);
  
  if(show_perf) {
    cputime2 = MPI_Wtime();
    ellapsed=1000*(cputime2-cputime1); average=P_average(ellapsed); max=P_max(ellapsed);
    }
  
  if(solver->rank==0 and show_perf) {
    printf("%s : solver %d (neq=%d) solving time average=%lfms max=%lfms \n", __func__, solver->id, 0, average, max);
    }
    
  if(track_call) printf("%s cpu=%d : ---> 5\n",__func__,solver->rank);
  if(CheckSynchronicity) status=mpi_synchronicity(solver->world_communicator, solver->world_rank, solver->world_size, __FUNCTION__ , __LINE__, false);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_factorize(SOLVER_t<double> *solver, packed_t<double> & matrix, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=factorize_template(solver, matrix, verbose, debug);
  return(status);
}
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_factorize(SOLVER_t< complex<double> > *solver, packed_t< complex<double> > & matrix, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=factorize_template(solver, matrix, verbose, debug);
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class T> int solve_template(SOLVER_t<T> *solver, T *RHS, T *x, int transposed, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  char *solver_name;
  int status;
  double cputime1,cputime2;
  double ellapsed, max, average;
  const int ntics=solver->options[OPTIONS_SHOW_SOLVE_PERF];
  const bool show_perf=(ntics>0);
  perftimer_t *perftimer=(perftimer_t *) solver->perftimer;
  
  if(perftimer->ntics!=ntics) {
    perftimer->destroy();
    if(ntics>0) perftimer->init(ntics);
    }
  
  if(MPI_Init_thread_DONE==0) TRAP_ERR_EXIT(-1,"MPI not initialized\n");
  
//   int proc_id, ncpus, ierr;
//   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
//   ierr = MPI_Comm_size(MPI_COMM_WORLD, &ncpus);
   
//   if(show_perf) cputime1 = MPI_Wtime();
  if(show_perf) perftimer->acquire_start();

  if(solver==0) TRAP_ERR_EXIT(-1,"solver not initialized; terminate job");
  
  solver_name=solver->name;
     
  if(verbose) printf("%s real-system solver, transpose=%d \n",solver_name, transposed);
  
  switch (solver->id) {
    case SOLVER_ID_DIAGONAL:     /* no real solver used for diagonal matrices */
      break;

    case SOLVER_ID_DOMESTIC:
    case SOLVER_ID_SUNPERF:
    case SOLVER_ID_HYPRE:
      TRAP_ERR_EXIT(-1,"solver %s not implemented; terminate job", solver_name);
      break;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    external, static solvers
    
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    case SOLVER_ID_LAPACK:
    case SOLVER_ID_UMFPACK:  
    case SOLVER_ID_PETSC_GMRES:
    case SOLVER_ID_MUMPS: 
    case SOLVER_ID_MUMPS_SYM:
    case SOLVER_ID_SpDOMESTIC:
    case SOLVER_ID_HIPS:
    case SOLVER_ID_PASTIX:
    case SOLVER_ID_PASTIX_SEQUENTIAL:
      if(solver->solve==0) TRAP_ERR_EXIT(-1,"solver %s solve function not initialized; terminate job",solver_name);
      status=solver->solve(solver, RHS, x, transposed, verbose, debug);
      break;  

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    ad hoc, dynamically implemented solvers
    
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    default:
      if(solver->solve==0) TRAP_ERR_EXIT(-1,"solver %s solve function not initialized; terminate job",solver_name);
      status=solver->solve(solver, RHS, x, transposed, verbose, debug);
      break;
    }
    
//   if(show_perf) {
//     cputime2 = MPI_Wtime();
//     ellapsed=1000*(cputime2-cputime1); average=P_average(ellapsed); max=P_max(ellapsed);
//     }
//   
//   if(solver->rank==0 and show_perf) {
//     printf("%s : solver %d (neq=%d) solving time average=%lfms max=%lfms \n", __func__, solver->id, 0, average, max);
//     }

  if(show_perf) {
    bool ready;
    ready=perftimer->acquire_stop();
    if(ready and solver->rank==0) printf("%s : solver %d (neq=%d) solving time %s\n", __func__, solver->id, 0, perftimer->report.c_str());
    }
    
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_solve(SOLVER_t<double> *solver, double *RHS, double *x, int transposed, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=solve_template(solver, RHS, x, transposed, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_solve(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x, int transposed, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=solve_template(solver, RHS, x, transposed, verbose, debug);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int solver_init_template(int solver_id, MPI_Comm world, int nprocs, int nunknowns, SOLVER_t<T>* & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  solver initializing

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int ierr;
  int status;
    
  int  argc=0;
  char **argv=0;
  int  required;    /* MPI thread level required       */
  int  provided;    /* MPI thread level provided       */
  int  rank, size;
  perftimer_t *perftimer;
  
//   const int ntics=solver->options[OPTIONS_SHOW_SOLVE_PERF];
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  solver structure default initialisation
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  solver=new SOLVER_t<T>;
    
  char *solver_name=get_SolverName(solver_id);
  
  if(solver_name==0) TRAP_ERR_EXIT(-1,"solver %d not known",solver_id);
  
  solver->name = strdup(solver_name);
  
/*------------------------------------------------------------------------------
  default setting, possibly set to 1 in case of PASTIX solver */
  solver->transpose_matrix=0;
 
  solver->parameters =0;
  
  solver->RHS_distribution=SOLVER_CENTRALIZED;
  solver->solution_distribution=SOLVER_CENTRALIZED;
  
  solver->communicator=0;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  perftimer
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  perftimer=new perftimer_t;
//   if(ntics>0) perftimer->init(ntics);
  solver->perftimer=(void *) perftimer;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  minimal MPI initialisation for sequential use
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  switch(solver_id) {
    case SOLVER_ID_MUMPS:
    case SOLVER_ID_MUMPS_SYM:
    case SOLVER_ID_HIPS:
    case SOLVER_ID_PASTIX:
    case SOLVER_ID_PETSC_GMRES:
    case SOLVER_ID_MAPHYS:
    default:
      required = MPI_THREAD_MULTIPLE;
      provided = -1;
      if(MPI_Init_thread_DONE==0) {
        MPI_Init_thread(&argc, &argv, required, &provided);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
          switch (provided) {
            case MPI_THREAD_SINGLE:
              printf("MPI_Init_thread level = MPI_THREAD_SINGLE\n");
              break;
            case MPI_THREAD_FUNNELED:
              printf("MPI_Init_thread level = MPI_THREAD_FUNNELED\n");
              break;
            case MPI_THREAD_SERIALIZED:
              printf("MPI_Init_thread level = MPI_THREAD_SERIALIZED\n");
              break;
            case MPI_THREAD_MULTIPLE:
              printf("MPI_Init_thread level = MPI_THREAD_MULTIPLE\n");
              break;
            default:
              printf("MPI_Init_thread level = ???\n");
              break;
            }
          }
        MPI_Init_thread_DONE=1;
        }
      break;
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  dedicated MPI communicators initialisation
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(world==0) world=MPI_COMM_WORLD;
  
  ierr = MPI_Comm_rank(world, &rank);
  ierr = MPI_Comm_size(world, &size);
  
  solver->world_communicator=world;
  solver->world_rank=rank;
  solver->world_size=size;
  
  if(debug) printf("%s cpu %d (in calling group of %d) : start, id=%d ndof=%d\n",__FUNCTION__,rank, size, solver_id, nunknowns);  

  int color=0, key=0;
  MPI_Comm communicator=MPI_COMM_NULL;
  rank=-1;

/*------------------------------------------------------------------------------
  default value */
  if(nprocs==-1) nprocs=size;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  handle communicators (duplicate/split) :
  
  in some cases (such as T-UGOm subcycling, which sub-domains would be to small
  to have efficient parallel computing with too many processors), we may need 
  less processors than initially requested at run launch.
  
  see http://mpitutorial.com/tutorials/introduction-to-groups-and-communicators
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(size==nprocs) {
    ierr=MPI_Comm_dup(world, &communicator);
    if(ierr!=0) TRAP_ERR_EXIT(-1,"cpu %d : MPI_Comm_split failed", solver->world_rank);
    ierr=MPI_Comm_rank(communicator, &rank);
    }
  else {
    if(nunknowns==0) color=MPI_UNDEFINED; 
    ierr=MPI_Comm_split(world, color, key, &communicator);
    if(ierr!=0) TRAP_ERR_EXIT(-1,"cpu %d : MPI_Comm_dup failed", solver->world_rank);
    if(nunknowns!=0) ierr=MPI_Comm_rank(communicator, &rank);
    }
  
  solver->communicator=communicator;
  solver->rank=rank;
  solver->size=size;
  
//   printf("%s cpu %d (in calling group of %d) : communictor=%d (MPI_COM_WORLD=%d)\n",__FUNCTION__,rank, size,communicator, MPI_COMM_WORLD);
  
  switch(solver_id) {
    case SOLVER_ID_DIAGONAL:
      solver->initialize=diagonal_initialize;
      break;

    case SOLVER_ID_DOMESTIC:
    case SOLVER_ID_SUNPERF:
    case SOLVER_ID_HYPRE:
      TRAP_ERR_EXIT(-1,"solver %s not implemented", solver_name);
      break;

    case SOLVER_ID_LAPACK:
      solver->initialize=lapack_initialize;
      break;

    case SOLVER_ID_UMFPACK:  
      solver->initialize=umfpack_initialize;
      break;

    case SOLVER_ID_PETSC_GMRES:
      solver->initialize=petsc_initialize;
      break;  

    case SOLVER_ID_MUMPS: 
      solver->initialize=mumps_initialize;
      break;

    case SOLVER_ID_MUMPS_SYM:
      solver->initialize=mumps_initialize;
      break;

    case SOLVER_ID_SpDOMESTIC:
      solver->initialize=spdomestic_initialize;
      break;

    case SOLVER_ID_HIPS:
      solver->initialize=hips_initialize;
      break;  

    case SOLVER_ID_PASTIX:
      solver->initialize=pastix_initialize;
      break;  

    default:
      status=solver_initialize_dynamic(*solver, solver_id, verbose, debug);
      break;
    }
  
//   debug=true;
  if(debug) printf("%s : cpu=%3d (cpu=%3d), solver_id=%d name=%s\n",__func__, solver->rank, solver->world_rank, solver_id, solver_name);
    
  if(solver->initialize==0) TRAP_ERR_EXIT(-1,"cpu=%d : unknown/unimplemented solver %s, id=%d", solver->rank, solver_name, solver_id);
    
  if(CheckSynchronicity) status=mpi_synchronicity(solver->world_communicator, solver->world_rank, solver->world_size, __FUNCTION__ , __LINE__, false);
  
  status=solver->initialize(*solver, verbose, debug);
    
  if(CheckSynchronicity) status=mpi_synchronicity(solver->world_communicator, solver->world_rank, solver->world_size, __FUNCTION__ , __LINE__, false);
  
  free(solver_name);

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_init(int solver_id, MPI_Comm world, int nprocs, int nunknowns, SOLVER_t<double>* & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=solver_init_template(solver_id, world, nprocs, nunknowns, solver, verbose, debug);

  return(status);
}

 
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_init(int solver_id, MPI_Comm world, int nprocs, int nunknowns, SOLVER_t< complex<double> >* & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=solver_init_template(solver_id, world, nprocs, nunknowns, solver, verbose, debug);

  return(status);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_init(int solver_id, size_t world, int nprocs, int nunknowns, SOLVER_t<double>* & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  interface needed for non-MPI compilation of T-UGOm

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  
  status=solver_init_template(solver_id, (MPI_Comm) world, nprocs, nunknowns, solver, verbose, debug);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_init(int solver_id, size_t world, int nprocs, int nunknowns, SOLVER_t< complex<double> >* & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  interface needed for non-MPI compilation of T-UGOm

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  
  status=solver_init_template(solver_id, (MPI_Comm) world, nprocs, nunknowns, solver, verbose, debug);

  return(status);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   template <typename T> int solver_init_template(int solver_id, int nprocs, packed_t<T> & packed, packed_t<T> & P, MPI_Comm world, SOLVER_t<T>* & solver, int verbose, bool debug)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   
//   status=solver_init(solver_id, world, nprocs, packed.nrows, solver, verbose, debug);
//   
//   status=factorize(solver, verbose, debug);
//   
//   return(status);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int solver_init(int solver_id, int nprocs, packed_t<double> & packed, packed_t<double> & P, MPI_Comm world, SOLVER_t<double>* & solver, int verbose, bool debug)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   
//   status=solver_init(solver_id, world, nprocs, packed.nrows, solver, verbose, debug);
//   
//   status=factorize(solver, packed, verbose, debug);
//   
//   return(status);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int solver_init(int solver_id, int nprocs, packed_t< complex<double> > & packed, packed_t< complex<double> > & P, MPI_Comm world, SOLVER_t< complex<double> >* & solver, int verbose, bool debug)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   
//   status=solver_init(solver_id, world, nprocs, packed.nrows, solver, verbose, debug);
//   
//   status=factorize(solver, packed, verbose, debug);
//   
//   return(status);
// }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T>int solver_terminate_template(SOLVER_t<T> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0, ierr;
  MPI_Comm world;

  if(solver==0) return(0);

  switch (solver->id) {

    case SOLVER_ID_DOMESTIC:
    case SOLVER_ID_SUNPERF:
    case SOLVER_ID_HYPRE:
      TRAP_ERR_EXIT(-1,"solver %d not implemented", solver->id);
      break;

    case SOLVER_ID_LAPACK:
    case SOLVER_ID_DIAGONAL:
    case SOLVER_ID_UMFPACK:  
    case SOLVER_ID_PETSC_GMRES:
    case SOLVER_ID_MUMPS: 
    case SOLVER_ID_MUMPS_SYM:
    case SOLVER_ID_SpDOMESTIC:
    case SOLVER_ID_HIPS:
    case SOLVER_ID_PASTIX:
      break;  

    default:
      int id=solver_identify_dynamic(*solver, solver->id, verbose, debug);
      if(id==-1) TRAP_ERR_EXIT(-1,"solver %s not implemented", solver->id);
      break;
    }
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  solver-dependant terminate
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(solver->terminate==0) TRAP_ERR_EXIT(-1,"cpu=%d : unknown solver, id=%d", solver->rank, solver->id);
  
  status=solver->terminate(solver, verbose, debug);
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  communicator
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#if HAVE_MPI
  if(solver->rank!=-1) {
    if(solver->communicator==MPI_COMM_WORLD) {
      if(debug) printf("%s cpu %d (in calling group of %d) : alert communicator=%d (MPI_COM_WORLD=%d)\n",__FUNCTION__,solver->rank, solver->size,solver->communicator, MPI_COMM_WORLD);
      }
    else status=MPI_Comm_free(&solver->communicator);
    }
#endif

  if(status!=0) {
    printf("%s : solver =%d status=%d\n",__func__, solver->id, status);
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  perftimer
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(solver->perftimer!=0) {
    perftimer_t *perftimer=(perftimer_t *) solver->perftimer;
    perftimer->destroy();
    delete perftimer;
    solver->perftimer=0;
    }
    
  solver->destroy();

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  TODO : re-activate fater changing passing argument to SOLVER_t<T> * & solver

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
//   delete solver;
//   solver=0;

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_terminate(SOLVER_t<double> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=solver_terminate_template(solver, verbose, debug);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_terminate(SOLVER_t< complex<double> > *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=solver_terminate_template(solver, verbose, debug);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_terminate(SOLVER_t<double> *solver)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int verbose=0;
  bool debug=false;
  
  status=solver_terminate_template(solver, verbose, debug);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_terminate(SOLVER_t< complex<double> > *solver)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int verbose=0;
  bool debug=false;
  
  status=solver_terminate_template(solver, verbose, debug);

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

  int factorize(solver_t *slv, triplet *Mat, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t<double> solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);

  slv->Mat = Mat; 
      
  solver.d_import2packed(slv);

  status=solver_factorize(&solver, *solver.packed, verbose, debug);
  
  solver.d_export_packed(slv);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorizez(solver_t *slv, tripletz *Mat, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t< complex<double> > solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  slv->Mat = Mat;
  
  solver.z_import2packed(slv);

  status=solver_factorize(&solver, *solver.packed, verbose, debug);
  
  solver.z_export_packed(slv);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solve(solver_t *slv, double *RHS, int transposed)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int verbose=0;
  bool debug=false;
  
  SOLVER_t<double> solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  solver.d_import2packed(slv);

  status=solver_solve(&solver, RHS, 0, transposed, verbose, debug);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solvez(solver_t *slv, complex<double> *RHS, int transposed)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int verbose=0;
  bool debug=false;
  SOLVER_t< complex<double> > solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  solver.z_import2packed(slv);

  status=solver_solve(&solver, RHS, 0, transposed, verbose, debug);
  
  return(status); 
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  solver_t *init_solver(int solver_id, MPI_Comm world, int nprocs, int nunknowns, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  solver_t *slv;
  SOLVER_t<double> *solver;
//   HYPERMATRIX_t<double> matrix;
  
//   status=solver_init(solver_id, matrix, world, nprocs, nunknowns, solver, verbose, debug);
  status=solver_init(solver_id, world, nprocs, nunknowns, solver, verbose, debug);
  
  slv = (solver_t *) malloc(sizeof(solver_t));
  solver->export_parameters(slv);

  return(slv); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  solver_t *initz_solver(char *solver_name, int typsym, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  bool debug=false;
  MPI_Comm world=0;
  int nprocs=-1;
  int nunknowns;
  int solver_id;
  solver_t *slv;
  SOLVER_t< complex<double> > *solver;
//   HYPERMATRIX_t< complex<double> > matrix;
  
  solver_id=get_SolverId(solver_name, debug);

//   status=solver_init(solver_id, matrix, world, nprocs, nunknowns, solver, verbose, debug);
  status=solver_init(solver_id, world, nprocs, nunknowns, solver, verbose, debug);
    
  slv = (solver_t *) malloc(sizeof(solver_t));
  solver->export_parameters(slv);

  return(slv); 
}


#endif
