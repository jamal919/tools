
/**************************************************************************

  POC-SOLVERS interface, 2006-2012

  Part of the Unstructured Ocean Grid initiative and SYMPHONIE suite

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      PhD, LEGOS, Toulouse, France

***************************************************************************/

#define VERBOSE

#include <string.h>
#include <omp.h>

#include "solvers-interface.h"
#include "poc-solvers.h"

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match MUMPS documentation */

extern int MPI_Init_thread_DONE;
extern int nbsyslinearz2;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solve(solver_t *slv, double *RHS, int transposed, int verbose)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  char *solver_name;
  int rtn;
  double cputime1,cputime2;
  int proc_id, ncpus, ierr;
  bool debug=false;
  
#ifdef HAVE_MPI
   if(MPI_Init_thread_DONE==1) {
     ierr = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
     ierr = MPI_Comm_size(MPI_COMM_WORLD, &ncpus);
     cputime1 = MPI_Wtime();
     }
#endif
  
  solver_name=slv->name;
     
  if(verbose) printf("%s real-system solver, transpose=%d \n",solver_name, transposed);

  if (strcmp(solver_name,"MUMPS") == 0) { 
    rtn = solve_mumps(slv, RHS,transposed,verbose,false);
    }
  else if  (strcmp(solver_name,"UMFPACK") == 0) {
    rtn=solve_umfpack(slv, RHS,transposed);
    }
  else if  (strcmp(solver_name,"LAPACK") == 0) {
    rtn=solve_lapack(slv, RHS,transposed);
    }
  else if  (strcmp(solver_name,"SpDOMESTIC") == 0) {
    rtn=solve_spdomestic(slv, RHS,transposed);
    }
  else if  (strcmp(solver_name,"HIPS") == 0) {
    rtn=solve_hips(slv, RHS,transposed,debug);
    }
  else if  (strcmp(solver_name,"PASTIX") == 0) {
    if(ncpus>1) {  
    
#if PASTIX_DISTRIBUTED
      MPI_Barrier(MPI_COMM_WORLD);
      rtn=solve_pastix_parallel(slv, RHS,transposed, verbose, debug);
#else
      rtn=solve_pastix_parallel(slv, RHS,transposed, verbose, debug);
#endif
      }
    else {
      rtn=solve_pastix_sequential(slv, RHS,transposed);
      }     
    }
  else if  (strcmp(solver_name,"PASTIX-SEQUENTIAL") == 0) {
    rtn=solve_pastix_sequential(slv, RHS,transposed);
    }
  else {
    printf("solve, unknown solver : %s \n",solver_name);
    rtn=-1;
    }
     
#ifdef  HAVE_MPI
//   if(MPI_Init_thread_DONE==1 && verbose==1) {
  if(MPI_Init_thread_DONE==1) {
    cputime2 = MPI_Wtime();
    if (verbose==1 && proc_id == 0) {
      printf("--------------------------------------------\n");
      printf("%s: %s solving ellapsed time= %e\n",__FUNCTION__, solver_name, cputime2-cputime1);
      printf("--------------------------------------------\n");
      }
    }
#endif

  return(rtn);
}

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  depecrated 1.0 sources
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#if OLD_SOLVER_VERSION
// if activated, use historical (1.0) code (still available in the c-code branch)
// not maintained, kept only for testing purposes, aimed to be removed

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solve_mumps(solver_t *slv,double *RHS, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#ifdef MUMPS
  DMUMPS_STRUC_C *id;
  char *solver_name;
  int ierr,rank;
  double cput1, cput2, cput3;
  
#ifdef HAVE_MPI
  cput1 = MPI_Wtime();
#endif
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   solver_name=slv->name;
//   if (strcmp(solver_name,"MUMPS") == 0) {
 
  id = (DMUMPS_STRUC_C *) slv->parameters;
  if (transposed == 1) {
    (*id).ICNTL(9)=2;
    }
  if (rank == 0) {
    (*id).rhs = RHS;
    }
  (*id).job=3;
  dmumps_c(id);
  
#ifdef HAVE_MPI
  cput2 = MPI_Wtime(); 
  cput3 = cput2-cput1;
  if (rank == 0) {
    if(verbose==1) printf("  DMUMPS resolution time =%lf\n",cput3);
  }
#endif
  ierr=(*id).info[0];
  if(ierr==-1) {
    printf("cpu %d: DMUMPS solver failed, status=%d due to cpu %d\n",rank, ierr, id->info[1]);
    }
  if(ierr==-22) {
    printf("cpu %d: DMUMPS solver failed, status=%d array fault %d\n",rank, ierr, id->info[1]);
    }
  return(ierr);
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solve_umfpack(solver_t *slv,double *RHS,int transposed)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1,bcl;
#ifdef UMFPACK
  UMFPACK_STRUC_C *id;
  char *solver_name;
  double *Info, *Control, *Ax, *tmp;
  int *Ap, *Ai,n,nnz;
  void *Symbolic, *Numeric ;
  triplet *Mat;
  int INCXY;
 
  status = 0;
  solver_name=slv->name;

  id = (UMFPACK_STRUC_C *) slv->parameters;
  Symbolic=id->Symbolic;
  Numeric=id->Numeric;
  Info=id->Info;
  Control=id->Control;
  Mat=(triplet *) slv->Mat;
  Ap=Mat->i;
  Ai=Mat->j;
  Ax=Mat->x;
  n=Mat->nrow;
  nnz=Mat->nnz;
  tmp = (double *) malloc (n * sizeof (double)) ;
  INCXY=1;
  /* BLASC */
  /*dcopy(n,RHS,1,tmp,1);*/
  /*BLASF*/
  /*dcopy_(n,RHS,INCXY,tmp,INCXY);*/
  /*no blas Xcopy*/
  for (bcl=0;bcl<n;bcl++) tmp[bcl]=RHS[bcl];
  if (transposed == 1) status = umfpack_di_solve (UMFPACK_At, Ap, Ai, Ax, RHS, tmp,
    Numeric, Control, Info) ;
  else status = umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, RHS, tmp,
    Numeric, Control, Info) ;
  
  if (status < 0) {
    umfpack_di_report_info (Control, Info) ;
    umfpack_di_report_status (Control, status) ;
    }
    
  free(tmp);
#endif
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solve_lapack(solver_t *slv, double *RHS, int transposed)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#ifdef LAPACK
  int status,nrhs=1,trans_len;
  char cjob = 'N';
  triplet *Mat;

  Mat=(triplet *) slv->Mat;

  dgbtrs_(&cjob,  &Mat->nrow, &Mat->ml, &Mat->mu, &nrhs, 
          Mat->a, &Mat->lda,Mat->pivot, RHS, &Mat->nrow, &status);

  return(status);
 
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solve_spdomestic(solver_t *slv,double *RHS,int transposed)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0,bcl,n;
  triplet *Mat;
  cs *A;
  css *S ;
  csn *N ;
  double *x;

  Mat = (triplet *) slv->Mat;
  A = Mat->A;
  S = Mat->S;
  N = Mat->N;
  n = A->n ;
  x = (double *) cs_malloc (n, sizeof (double)) ;    /* get workspace */

  bcl = N->L->nz;

  cs_ipvec (N->pinv, RHS, x, n) ;            /* x = b(p) */
  if (transposed) {
    cs_utsolve (N->U, x) ;                   /* x = U\x  */
    cs_ltsolve (N->L, x) ;                   /* x = L\x  */
    }
  else {
    cs_lsolve (N->L, x) ;                    /* x = L\x  */
    cs_usolve (N->U, x) ;                    /* x = U\x  */
    }
  cs_ipvec (S->q, x, RHS, n) ;               /* b(q) = x */
  cs_free (x) ;
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solve_hips(solver_t *slv, double *RHS, int transposed, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#ifdef HIPS
#ifdef HAVE_MPI
  dHIPS_STRUC_C*id;
  char *solver_name;
  int proc_id, nprocs, myid=-1, ierr, status, n, nnz,bcl, nglob;
  double cput1, cput2, cput3;
  INTS id_hips;
  INTS *unknownlist;
  MPI_Comm world;
  
//   debug=true;
  
  if(slv->rank==-1) return(0);

  id = (dHIPS_STRUC_C*) slv->parameters; 
  
  id_hips=id->id_hips;

  world=slv->communicator;
  n=id->n;
  
  cput1 = MPI_Wtime();
  ierr  = MPI_Comm_rank(world, &proc_id);
  ierr  = MPI_Comm_size(world, &nprocs);
  
  if(debug) printf("%s cpu=%d (over %d) : id_hips=%d\n", __FUNCTION__, proc_id, nprocs, id_hips);

  nnz=id->nnz;
  unknownlist = id->unknownlist;
    
  /** TEst set RHS like x=1 */
/*   x  = (double *) malloc(sizeof(double)*nnz); */
/*   x2 = (double *) malloc(sizeof(double)*nnz); */
/*   for (bcl=0;bcl < nnz;bcl++) x[bcl]=1.0; */

/*   dHIPS_MatrixVectorProduct (id_hips, x, x2); */
/*   for (bcl=0;bcl < n;bcl++) RHS[bcl]=x2[bcl]; */

  if(debug) printf("%s cpu=%d (over %d) : set RHS\n", __FUNCTION__, proc_id, nprocs, id_hips);
  
//   int mode=SOLVER_CENTRALIZED; 
  int mode=slv->RHS_distribution; 
  switch(mode) {
    case SOLVER_CENTRALIZED:
      ierr = dHIPS_SetGlobalRHS(id_hips, RHS, 0, dHIPS_ASSEMBLY_OVW);
      dHIPS_ExitOnError(ierr);
      break;
    case SOLVER_DISTRIBUTED:
      ierr = dHIPS_SetRHS (id_hips, n, unknownlist, RHS,  dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_FOOL);
//       ierr = dHIPS_SetRHS (id_hips, n, unknownlist, RHS,  dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_RESPECT);
      dHIPS_ExitOnError(ierr);
      break;
    }

  cput1 = MPI_Wtime();
  
  if(debug) printf("%s cpu=%d (over %d) : get solution\n", __FUNCTION__, proc_id, nprocs, id_hips);
  
//   mode=SOLVER_CENTRALIZED;
  mode=slv->solution_distribution; 
  switch(mode) {
    case SOLVER_CENTRALIZED:
      ierr = dHIPS_GetGlobalSolution(id_hips,RHS,-1);
      dHIPS_ExitOnError(ierr);
      break;
    case SOLVER_DISTRIBUTED:
      ierr = dHIPS_GetSolution(id_hips,n,unknownlist, RHS, dHIPS_ASSEMBLY_FOOL);
      dHIPS_ExitOnError(ierr);
      break;      
    }
  
//   dHIPS_SetOptionINT (id_hips,dHIPS_DISABLE_PRECOND, 1);

  cput2 = MPI_Wtime();
  cput3 = cput2-cput1;
  
#ifdef VERBOSE  
  if (myid == 0) {
    printf("  HIPS Resolution time =%lf \n",cput3);
  }
#endif 

  if(debug) printf("%s cpu=%d (over %d) : done...\n", __FUNCTION__, proc_id, nprocs, id_hips);
  
  if(ierr==1) ierr=0;
  else {
    printf("%s cpu=%d (over %d) : failed...\n", __FUNCTION__, proc_id, nprocs, id_hips);
    }
  return(ierr);
#endif
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solvez_hips_complex(solver_t *slv, complex_t *RHS, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#ifdef HIPS
#ifdef HAVE_MPI
  dHIPS_STRUC_C*id;
  char *solver_name;
  int proc_id, nprocs, myid=-1, ierr, status, n, nnz,bcl, nglob;
  double cput1, cput2, cput3;
  INTS id_hips;
  INTS *unknownlist;
  MPI_Comm world;
  
//   debug=true;
  
  if(slv->rank==-1) return(0);

  id = (dHIPS_STRUC_C*) slv->parameters; 
  
  id_hips=id->id_hips;

  world=slv->communicator;
  n=id->n;
  
  cput1 = MPI_Wtime();
  ierr  = MPI_Comm_rank(world, &proc_id);
  ierr  = MPI_Comm_size(world, &nprocs);
  
  if(debug) printf("%s cpu=%d (over %d) : id_hips=%d\n", __FUNCTION__, proc_id, nprocs, id_hips);

  nnz=id->nnz;
  unknownlist = id->unknownlist;
    
  if(debug) printf("%s cpu=%d (over %d) : set RHS\n", __FUNCTION__, proc_id, nprocs, id_hips);
  
  int mode=slv->RHS_distribution; 
  switch(mode) {
    case SOLVER_CENTRALIZED:
      ierr = dHIPS_SetGlobalRHS(id_hips, (double *) RHS, 0, dHIPS_ASSEMBLY_OVW);
      dHIPS_ExitOnError(ierr);
      break;
    case SOLVER_DISTRIBUTED:
      ierr = dHIPS_SetRHS (id_hips, n, unknownlist, (double *) RHS,  dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_FOOL);
      dHIPS_ExitOnError(ierr);
      break;
    }

  cput1 = MPI_Wtime();
  
  if(debug) printf("%s cpu=%d (over %d) : get solution\n", __FUNCTION__, proc_id, nprocs, id_hips);
  
  mode=slv->solution_distribution; 
  switch(mode) {
    case SOLVER_CENTRALIZED:
      ierr = dHIPS_GetGlobalSolution(id_hips, (double *) RHS,-1);
      dHIPS_ExitOnError(ierr);
      break;
    case SOLVER_DISTRIBUTED:
      ierr = dHIPS_GetSolution(id_hips,n,unknownlist, (double *) RHS, dHIPS_ASSEMBLY_FOOL);
      dHIPS_ExitOnError(ierr);
      break;      
    }
  
//   dHIPS_SetOptionINT (id_hips,dHIPS_DISABLE_PRECOND, 1);

  cput2 = MPI_Wtime();
  cput3 = cput2-cput1;
  
#ifdef VERBOSE
  if (myid == 0) {
    printf("  HIPS Resolution time =%lf \n",cput3);
  }
#endif 

  if(debug) printf("%s cpu=%d (over %d) : done...\n", __FUNCTION__, proc_id, nprocs, id_hips);
  
  if(ierr==1) ierr=0;
  else {
    printf("%s cpu=%d (over %d) : failed...\n", __FUNCTION__, proc_id, nprocs, id_hips);
    }
  return(ierr);
#endif
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solve_pastix_sequential(solver_t *slv,double *RHS,int transposed)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,bcl,bcl2;
#ifdef PASTIX
  PASTIX_STRUC_C *id;
  char *solver_name;
  double *Ax;
  int *Ap, *Ai,n,nnz;
  triplet *Mat;
  int INCXY;
  int rank,sz;
  bool debug=false;
  
  status = 0;
  
#ifdef HAVE_MPI
   status = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   status = MPI_Comm_size(MPI_COMM_WORLD, &sz);
#endif  
  
  if(debug) printf("%s cpu=%3d : start...\n",__func__,rank);
  
  id = (PASTIX_STRUC_C *) slv->parameters;
  Mat=(triplet *) slv->Mat;
  Ap=Mat->i;
  Ai=Mat->j;
  Ax=Mat->x;
  n=Mat->nrow;
  nnz=Mat->nnz;

  id->iparm[IPARM_START_TASK] = API_TASK_SOLVE;
  id->iparm[IPARM_END_TASK]   = API_TASK_REFINE;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  in case where transpose matrix has been provided instead of direct one, then
  swap transposed solution flag
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(slv->transpose_matrix==1) transposed=1-transposed;
  
  id->iparm[IPARM_TRANSPOSE_SOLVE] = transposed;

  d_pastix(&(id->pastix_data), MPI_COMM_WORLD, n, Ap, Ai, Ax,
             id->perm, id->invp, RHS, 1, id->iparm, id->dparm);
  return(status);
#endif
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solve_pastix_parallel(solver_t *slv, double *RHS, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,bcl,bcl2;
  int ierr=-1;
#ifdef PASTIX
#ifdef HAVE_MPI
  PASTIX_STRUC_C *id;
  char *solver_name;
  double *Ax;
  int *Ap, *Ai,n,nnz;
  triplet *Mat;
  int INCXY;
  int rank,sz,i;
  int *loc2glb;
  
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &sz);
  
  status = 0;
  solver_name=slv->name;

  id = (PASTIX_STRUC_C *) slv->parameters;
  Mat=(triplet *) slv->Mat;
  Ap=Mat->i;
  Ai=Mat->j;
  Ax=Mat->x;
  n=Mat->nrow;
  nnz=Mat->nnz;
  
  loc2glb=Mat->loc2glob;
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(verbose==1) printf("%s cpu=%3d : n=%d nnz=%d\n",__FUNCTION__,rank, n,nnz);
 
  id->iparm[IPARM_START_TASK] = API_TASK_SOLVE;
  id->iparm[IPARM_END_TASK]   = API_TASK_REFINE;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  in case where transpose matrix has been provided instead of direct one, then
  swap transposed solution flag
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(slv->transpose_matrix==1) transposed=1-transposed;
  
  id->iparm[IPARM_TRANSPOSE_SOLVE] = transposed;

  d_dpastix(&(id->pastix_data), MPI_COMM_WORLD,
           n, Ap, Ai, Ax, loc2glb, id->perm, id->invp, RHS, 1, id->iparm, id->dparm);

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

#endif 
