
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
// #include "functions.h"
#include "poc-solvers.h"
#include "solvers-interface.h"

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  wrappers necessary to template use
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#ifdef PASTIX

void dpastix(pastix_data_t **pastix_data, MPI_Comm pastix_comm,
             PASTIX_INT n, PASTIX_INT *colptr, PASTIX_INT *row,
             double *avals, PASTIX_INT * loc2glob, PASTIX_INT *perm, PASTIX_INT *invp,
             double *b, PASTIX_INT rhs, PASTIX_INT *iparm, double *dparm) {
  
  d_dpastix(pastix_data, pastix_comm, n, colptr, row, avals, loc2glob, perm, invp, b, rhs, iparm,  dparm);
}

void dpastix(pastix_data_t **pastix_data, MPI_Comm pastix_comm,
             PASTIX_INT n, PASTIX_INT *colptr, PASTIX_INT *row,
             complex<double> *avals, PASTIX_INT * loc2glob, PASTIX_INT *perm, PASTIX_INT *invp,
             complex<double> *b, PASTIX_INT rhs, PASTIX_INT *iparm, double *dparm) {
  
  z_dpastix(pastix_data, pastix_comm, n, colptr, row, avals, loc2glob, perm, invp, b, rhs, iparm,  dparm);
}


void pastix(pastix_data_t **pastix_data, MPI_Comm pastix_comm,
             PASTIX_INT n, PASTIX_INT *colptr, PASTIX_INT *row,
             double *avals, PASTIX_INT *perm, PASTIX_INT *invp,
             double *b, PASTIX_INT rhs, PASTIX_INT *iparm, double *dparm) {
  
d_pastix(pastix_data, pastix_comm, n, colptr, row, avals, perm,  invp, b, rhs,  iparm,  dparm);
}

void pastix(pastix_data_t **pastix_data, MPI_Comm pastix_comm,
             PASTIX_INT n, PASTIX_INT *colptr, PASTIX_INT *row,
             complex<double> *avals, PASTIX_INT *perm, PASTIX_INT *invp,
             complex<double> *b, PASTIX_INT rhs, PASTIX_INT *iparm, double *dparm) {
  
  z_pastix(pastix_data, pastix_comm, n, colptr, row, avals, perm,  invp, b, rhs,  iparm,  dparm);
}

PASTIX_INT pastix_checkMatrix(MPI_Comm pastix_comm, PASTIX_INT verb, PASTIX_INT flagsym, PASTIX_INT flagcor,
                       PASTIX_INT n, PASTIX_INT **colptr, PASTIX_INT **row, double **avals,
                       PASTIX_INT **loc2glob, PASTIX_INT dof) {
  
  d_pastix_checkMatrix(pastix_comm, verb, flagsym, flagcor,  n, colptr, row, avals, loc2glob, dof);
}

PASTIX_INT pastix_checkMatrix(MPI_Comm pastix_comm, PASTIX_INT verb, PASTIX_INT flagsym, PASTIX_INT flagcor,
                       PASTIX_INT n, PASTIX_INT **colptr, PASTIX_INT **row, complex<double> **avals,
                       PASTIX_INT **loc2glob, PASTIX_INT dof) {
  
  z_pastix_checkMatrix(pastix_comm, verb, flagsym, flagcor,  n, colptr, row, avals, loc2glob, dof);
}

#endif

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  pastix interface
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int pastix_factorize_sequential_template(SOLVER_t<T> *solver, int verbose, bool debug) {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/    
  int status;
#ifdef PASTIX
  pastix_t<T> *id = (pastix_t<T> *) solver->parameters;
  char *solver_name;

  T *Ax;
  pastix_int_t *Ap, *Ai;
  int  n,nnz,nbthread;
  T *rhs         = NULL;                /* right hand side                           */
  int               verbosemode;        /* Level of verbose mode (0, 1, 2)           */
  char             *type        = NULL; /* type of the matrix                        */
  int    flagsym=0;
  pastix_int_t    ncol;                 /* Size of the matrix                        */
  
  MPI_Comm communicator=solver->communicator;

  n=solver->packed->nrows;
  nnz=solver->packed->nnz;

  if (solver->packed->pointer32==0) status=recast(solver->packed->pointer, solver->packed->pointer32, false, n+1);

  Ap=solver->packed->pointer32;
  Ai=(pastix_int_t *) solver->packed->incidence;
  Ax=solver->packed->x;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  Initialize parameters to default values
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  verbosemode=0;    
  
  pastix_checkMatrix(communicator, verbosemode, API_SYM_NO, API_YES, n, &Ap, &Ai, &Ax, NULL,1);

/*------------------------------------------------------------------------------
  get (possibly changed) new addresses */
  
  solver->packed->pointer32=Ap;
  solver->packed->incidence=Ai;
  solver->packed->x=Ax;
  
  
  id->iparm[IPARM_MODIFY_PARAMETER] = API_NO;
  
  pastix(&(id->pastix_data), communicator, n, Ap, Ai, Ax, id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
  
  status=id->iparm[IPARM_ERROR_NUMBER];
  if (status < 0)  { 
    printf("%s : status=%d\n",__func__,status);  
    return(-1);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  Customize some parameters
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
#ifdef OMP_H
  nbthread = omp_get_max_threads(); 
  if(debug) printf("%s : using %d threads (omp_get_max_threads)...................... \n",__func__,nbthread);
//   nbthread = omp_get_num_threads();
//   if(debug) printf("%s : using %d threads (omp_get_num_threads)...................... \n",__func__,nbthread);
#else
  nbthread = 1;  
  if(debug) printf("%s : using only 1 thread!...................... \n",__func__,nbthread);
#endif
  
  verbosemode=0;    
  id->iparm[IPARM_SYM]           = API_SYM_NO;
  id->iparm[IPARM_FACTORIZATION] = API_FACT_LU;
  id->iparm[IPARM_THREAD_NBR]    = nbthread;
  id->iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
  id->iparm[IPARM_LEVEL_OF_FILL] = 0;
  id->iparm[IPARM_START_TASK]    = API_TASK_ORDERING;
  id->iparm[IPARM_END_TASK]      = API_TASK_NUMFACT;
  
//   id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-8;
  id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-12;
  
  id->iparm[IPARM_VERBOSE]       = verbosemode;

  id->perm = new pastix_int_t[n];
  id->invp = new pastix_int_t[n];
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  Factorisation using centralized pastix
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  pastix(&(id->pastix_data), communicator, n, Ap, Ai, Ax, id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
  
  status=id->iparm[IPARM_ERROR_NUMBER];
  if (status < 0)  { 
    printf("%s : status=%d\n",__func__,status);  
    return(-1);
    }
 
  status=0;        
  return(status);
#else
  status=-1;
  printf("%s (line %d) : please compile poc-solvers library with -DPASTIX \n", __func__,__LINE__);
  return(status);
#endif 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int pastix_factorize_parallel_distributed_template(SOLVER_t<T> *solver, int verbose, bool debug) {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  int k,status;
#ifdef PASTIX
#ifdef HAVE_MPI
  pastix_t<T> *id = (pastix_t<T> *) solver->parameters;
  char *solver_name;
  T *Ax;
  pastix_int_t *Ap, *Ai;
  pastix_int_t  n,nnz;
  int  nbthread;
  T *rhs = NULL;               /* right hand side                      */
  int    verbosemode;          /* Level of verbose mode (0, 1, 2)      */
  char   *type        = NULL;  /* type of the matrix                   */
  int    bcl, bcl2;
  int    flagsym=0;
  pastix_int_t    ncol;        /* Size of the matrix                   */
  int   mpid,sz;
  pastix_int_t *loc2glb;
  
  MPI_Comm communicator=solver->communicator;

  MPI_Barrier(communicator);
          
  n=solver->packed->nrows;
  nnz=solver->packed->nnz;

  if (solver->packed->pointer32==0) status=recast(solver->packed->pointer, solver->packed->pointer32, false, n+1);

  Ap=solver->packed->pointer32;
  Ai=solver->packed->incidence;
  Ax=solver->packed->x;
  
  loc2glb=new pastix_int_t[n];

  for(k=0;k<n;k++) loc2glb[k]=solver->packed->loc2glob[k];

  MPI_Barrier(communicator);

  MPI_Comm_rank(communicator, &mpid);
  MPI_Comm_size(communicator, &sz);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Initialize parameters to default values
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  verbosemode=verbose;  
  
  pastix_checkMatrix(communicator, verbosemode, API_SYM_NO, API_YES, n, &Ap, &Ai, &Ax, &loc2glb, 1);
  
/*------------------------------------------------------------------------------
  get (possibly changed) new addresses */

  solver->packed->pointer32=Ap;
  solver->packed->incidence=Ai;
  solver->packed->x=Ax;
  
  id->iparm[IPARM_MODIFY_PARAMETER] = API_NO;
  
  dpastix(&(id->pastix_data), communicator, n, Ap, Ai, Ax, loc2glb, id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
  
  status=id->iparm[IPARM_ERROR_NUMBER];
  if (status < 0)  { 
    printf("%s : status=%d\n",__func__,status);  
    return(-1);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Customize some parameters
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  debug=true;
#ifdef OMP_H
//   nbthread = omp_get_max_threads(); 
//   if(debug) printf("%s : using %d threads (omp_get_max_threads)...................... \n",__func__,nbthread);
  nbthread = omp_get_num_threads();
  if(debug) printf("%s : using %d threads (omp_get_num_threads)...................... \n",__func__,nbthread);
#else
  nbthread = 1;  
  if(debug) printf("%s : using only 1 thread!...................... \n",__func__,nbthread);
#endif
    
//   verbose=1; // HERE
  
  verbosemode=verbose;    
  id->iparm[IPARM_SYM]           = API_SYM_NO;
  id->iparm[IPARM_FACTORIZATION] = API_FACT_LU;
  id->iparm[IPARM_THREAD_NBR]    = nbthread;
  id->iparm[IPARM_MATRIX_VERIFICATION] = API_YES;
  id->iparm[IPARM_LEVEL_OF_FILL] = 0;
  id->iparm[IPARM_START_TASK] = API_TASK_ORDERING;
  id->iparm[IPARM_END_TASK]   = API_TASK_NUMFACT;
  
//   id->iparm[IPARM_END_TASK]   = API_TASK_ORDERING; // HERE

//   id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-8;
  id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-12;
  
  id->iparm[IPARM_VERBOSE]            = verbosemode;

  /*iparm[IPARM_END_TASK]   = API_TASK_REFINE;*/

  id->perm = new pastix_int_t[n];
  
//   id->invp = malloc(ncolglob*sizeof(pastix_int_t));
  
  MPI_Barrier(communicator);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  Factorisation using distributed pastix
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(verbose==1) printf("%s pastix factorisation start\n",__FUNCTION__);
  
  dpastix(&(id->pastix_data), communicator, n, Ap, Ai, Ax, loc2glb, id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
    
  status=id->iparm[IPARM_ERROR_NUMBER];
  
  if(verbose==1) printf("%s cpu=%d: pastix factorisation status=\n",__FUNCTION__, solver->rank,status);
  
  if (status < 0)  { 
    printf("%s : status=%d\n",__func__,status);  
    return(-1);
    }
 
  status=0;    // ???    
  return(status);
#endif 
#else
  status=-1;
  printf("%s (line %d) : please compile poc-solvers library with -DPASTIX \n", __func__,__LINE__);
  return(status);
#endif 
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int pastix_factorize_sequential(SOLVER_t<double> *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=pastix_factorize_sequential_template(solver, verbose, debug);
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int pastix_factorize_sequential(SOLVER_t< complex<double> > *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=pastix_factorize_sequential_template(solver, verbose, debug);
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int pastix_factorize_parallel_distributed(SOLVER_t<double> *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=pastix_factorize_parallel_distributed_template(solver, verbose, debug);
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int pastix_factorize_parallel_distributed(SOLVER_t< complex<double> > *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=pastix_factorize_parallel_distributed_template(solver, verbose, debug);
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int pastix_factorize(SOLVER_t<double> *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  bool sequential=(solver->packed->loc2glob==0);
  
  if(sequential) {
    status=pastix_factorize_sequential_template(solver, verbose, debug);
    }
  else {
    status=pastix_factorize_parallel_distributed_template(solver, verbose, debug);
    }
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int pastix_factorize(SOLVER_t< complex<double> > *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  bool sequential=(solver->packed->loc2glob==0);
  
  if(sequential) {
    status=pastix_factorize_sequential_template(solver, verbose, debug);
    }
  else {
    status=pastix_factorize_parallel_distributed_template(solver, verbose, debug);
    }
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int pastix_solve_sequential_template(SOLVER_t<T> *solver, T *RHS,int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
#ifdef PASTIX
  pastix_t<T> *id = (pastix_t<T> *) solver->parameters;
  char *solver_name;
  T *Ax;
  pastix_int_t *Ap, *Ai,n,nnz;
  int INCXY;
 
  MPI_Comm communicator=solver->communicator;

  status = 0;
  solver_name=solver->name;

  Ap=solver->packed->pointer32;
  Ai=solver->packed->incidence;
  Ax=solver->packed->x;
  
  n=solver->packed->nrows;
  nnz=solver->packed->nnz;

  if(verbose==1) printf("%s : n=%d nnz=%d\n",__FUNCTION__, n,nnz);
  
  id->iparm[IPARM_START_TASK] = API_TASK_SOLVE;
  id->iparm[IPARM_END_TASK]   = API_TASK_REFINE;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  in case where transpose matrix has been provided instead of direct one, then
  swap transposed solution flag
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(solver->transpose_matrix==1) transposed=1-transposed;
  
  id->iparm[IPARM_TRANSPOSE_SOLVE] = transposed;

  pastix(&(id->pastix_data), communicator,  n, Ap, Ai, Ax, id->perm, id->invp, RHS, 1, id->iparm, id->dparm);
  
  status=id->iparm[IPARM_ERROR_NUMBER];
  if (status < 0)  { 
    printf("%s : status=%d\n",__func__,status);  
    return(-1);
    }
 
  
  solver->communicator=communicator;
  
#endif
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int pastix_solve_parallel_template(SOLVER_t<T> *solver, T *RHS,int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int ierr=-1;
#ifdef PASTIX
#ifdef HAVE_MPI
  pastix_t<T> *id = (pastix_t<T> *) solver->parameters;
  char *solver_name;
  T *Ax;
  int *Ap, *Ai,n,nnz;
  tripletz *Mat;
  int INCXY;
  int rank,sz,i;
  int *loc2glb;
  
  MPI_Comm communicator=solver->communicator;

  ierr = MPI_Comm_rank(communicator, &rank);
  ierr = MPI_Comm_size(communicator, &sz);
  
  status = 0;
  solver_name=solver->name;

  Ap=solver->packed->pointer32;
  Ai=solver->packed->incidence;
  Ax=solver->packed->x;
  
  n=solver->packed->nrows;
  nnz=solver->packed->nnz;

  loc2glb=solver->packed->loc2glob;
//   loc2glb=new pastix_int_t[n];
//   for(int k=0;k<n;k++) loc2glb[k]=solver->packed->loc2glob[k];

  MPI_Barrier(communicator);
  if(verbose==1) printf("%s cpu=%3d : n=%d nnz=%d\n",__FUNCTION__,rank, n,nnz);
 
  id->iparm[IPARM_START_TASK] = API_TASK_SOLVE;
  id->iparm[IPARM_END_TASK]   = API_TASK_REFINE;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  in case where transpose matrix has been provided instead of direct one, then
  swap transposed solution flag
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(solver->transpose_matrix==1) transposed=1-transposed;
  
  id->iparm[IPARM_TRANSPOSE_SOLVE] = transposed;

  dpastix(&(id->pastix_data), communicator, n, Ap, Ai, Ax, loc2glb, id->perm, id->invp, RHS, 1, id->iparm, id->dparm);

  status=id->iparm[IPARM_ERROR_NUMBER];
  if (status < 0)  { 
    printf("%s : status=%d\n",__func__,status);  
    return(-1);
    }
 
  ierr=0;
  
#endif 
#endif 
  return(ierr);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int pastix_solve_sequential(SOLVER_t<double> *solver, double *RHS,int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=pastix_solve_sequential_template(solver, RHS, transposed, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int pastix_solve_sequential(SOLVER_t< complex<double> > *solver, complex<double> *RHS,int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=pastix_solve_sequential_template(solver, RHS, transposed, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int pastix_solve_parallel(SOLVER_t<double> *solver, double *RHS,int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=pastix_solve_parallel_template(solver, RHS, transposed, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int pastix_solve_parallel(SOLVER_t< complex<double> > *solver, complex<double> *RHS,int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=pastix_solve_parallel_template(solver, RHS, transposed, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int pastix_solve(SOLVER_t<double> *solver, double *RHS, double *x,int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  bool sequential=(solver->packed->loc2glob==0);
  
  if(sequential) {
    status=pastix_solve_sequential_template(solver, RHS, transposed, verbose, debug);
    }
  else {
    status=pastix_solve_parallel_template(solver, RHS, transposed, verbose, debug);
    }
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int pastix_solve(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x,int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  bool sequential=(solver->packed->loc2glob==0);
  
  if(sequential) {
    status=pastix_solve_sequential_template(solver, RHS, transposed, verbose, debug);
    }
  else {
    status=pastix_solve_parallel_template(solver, RHS, transposed, verbose, debug);
    }
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int pastix_initialize_template(SOLVER_t<T> & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  MPI_Comm world;
  int rank;
  
  solver.communicator=MPI_COMM_WORLD; // TESTING
  
#ifdef PASTIX
  
  pastix_t<T> *parameters=new pastix_t<T>;
  
  solver.parameters=parameters;
  
  solver.targeted_packing    = PACKED;
  solver.targeted_ordering   = CSC;
  solver.targeted_addressing = MATRIX_POINTER;
  solver.targeted_numbering  = 1;
  
  solver.initialize = pastix_initialize;
  solver.terminate  = pastix_terminate;
  solver.factorize  = pastix_factorize;
  solver.solve      = pastix_solve;
  
  solver.id=SOLVER_ID_PASTIX;
  
  return(0);
  
#else 
  printf("Please compile poc-solvers library with -DPASTIX \n");
  exit(-1);
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int pastix_initialize(SOLVER_t<double> & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=pastix_initialize_template(solver, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int pastix_initialize(SOLVER_t< complex<double> > & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=pastix_initialize_template(solver, verbose, debug);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int pastix_terminate_template(SOLVER_t<T> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  MPI_Comm communicator;
  
#ifdef PASTIX
  pastix_t<T> *parameters = (pastix_t<T> *) solver->parameters;
  
  int *Ap=(int *) solver->packed->pointer;
  int *Ai=(int *) solver->packed->incidence;
  T* Ax=solver->packed->x;
  T* rhs;
  int n;
  
  parameters->iparm[IPARM_START_TASK] = API_TASK_CLEAN;
  parameters->iparm[IPARM_END_TASK]   = API_TASK_CLEAN;
  
  communicator=solver->communicator;

  pastix(&(parameters->pastix_data), communicator, n, Ap, Ai, Ax,
         parameters->perm,  parameters->invp, rhs, 1, 
         parameters->iparm, parameters->dparm);
  
  ptr_delete(parameters->perm);
  ptr_delete(parameters->invp);
  
  delete parameters;
  solver->parameters=0;
  
  status=0;
  return(status);
#else
  return(status);
#endif 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int pastix_terminate(SOLVER_t<double> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=pastix_terminate_template(solver, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int pastix_terminate(SOLVER_t< complex<double> > *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=pastix_terminate_template(solver, verbose, debug);
  
  return(status);
}



/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  solver interface compatibility with older versions
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#if OLD_SOLVER_VERSION
#else

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorize_pastix_sequential(solver_t *slv, triplet *Mat, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t<double> solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);

  slv->Mat = Mat; 
      
  solver.d_import2packed(slv);
  
  status=pastix_factorize_sequential(&solver, verbose, debug);
  
  solver.d_export_packed(slv);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorizez_pastix_sequential(solver_t *slv, tripletz *Mat, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t< complex<double> > solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  slv->Mat = Mat;
  
  solver.z_import2packed(slv);

  status=pastix_factorize_sequential(&solver, verbose, debug);
  
  solver.z_export_packed(slv);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorize_pastix_parallel_distributed(solver_t *slv, triplet *Mat, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t<double> solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  slv->Mat = Mat;
  
  solver.d_import2packed(slv);

  status=pastix_factorize_parallel_distributed(&solver, verbose, debug);
  
  solver.d_export_packed(slv);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorizez_pastix_parallel_distributed(solver_t *slv, tripletz *Mat, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t< complex<double> > solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  slv->Mat = Mat;
  
  solver.z_import2packed(slv);

  status=pastix_factorize_parallel_distributed(&solver, verbose, debug);
  
  solver.z_export_packed(slv);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solve_pastix_sequential(solver_t *slv, double *RHS, int transposed)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int verbose=0;
  bool debug=false;
  
  SOLVER_t<double> solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  solver.d_import2packed(slv);

  status=pastix_solve_sequential(&solver, RHS, transposed, verbose, debug);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solvez_pastix_sequential(solver_t *slv, complex<double> *RHS, int transposed, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t< complex<double> > solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  solver.z_import2packed(slv);

  status=pastix_solve_sequential(&solver, RHS, transposed, verbose, debug);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solve_pastix_parallel(solver_t *slv, double *RHS, int transposed, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t<double> solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  solver.d_import2packed(slv);

  status=pastix_solve_parallel(&solver, RHS, transposed, verbose, debug);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solvez_pastix_parallel(solver_t *slv, complex<double> *RHS, int transposed, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t< complex<double> > solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  solver.z_import2packed(slv);

  status=pastix_solve_parallel(&solver, RHS, transposed, verbose, debug);
  
  return(status); 
}

#endif
