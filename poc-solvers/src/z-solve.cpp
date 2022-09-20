
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

  int solvez(solver_t *slv, complex<double> *RHS, int transposed, int verbose)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  char *solver_name;
  int rtn;
  double cputime1, cputime2;
  int ierr, proc_id, ncpus;
  solver_name=slv->name;
  bool debug=false; 

#ifdef  HAVE_MPI
  if(MPI_Init_thread_DONE==1) {
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &ncpus);
    cputime1 = MPI_Wtime();
    }
#endif
   
  if(verbose) printf("%s complex-system solver, transpose=%d \n",solver_name, transposed);
   
  if (strcmp(solver_name,"MUMPS") == 0) { 
    rtn = solvez_mumps(slv, RHS,transposed, verbose, debug);
    }
  else if  (strcmp(solver_name,"UMFPACK") == 0) {
    rtn=solvez_umfpack(slv, RHS,transposed);
    }
  else if  (strcmp(solver_name,"LAPACK") == 0) {
    /*rtn=solvez_lapack(slv, RHS,transposed);*/
    printf("solver not implemented for complex matrices : %s \n",solver_name);
    }
  else if  (strcmp(solver_name,"SpDOMESTIC") == 0) {
    rtn=solvez_spdomestic(slv, RHS,transposed);
    }
  else if  (strcmp(solver_name,"HIPS") == 0) {
#ifdef __DEBUG_HIPS
    debug=true;
#endif
    rtn=solvez_hips(slv, RHS,transposed, verbose, debug);
//     rtn=solvez_hips_complex(slv, RHS,transposed, verbose, debug);
    }
  else if  (strcmp(solver_name,"PASTIX") == 0) {
    if(ncpus>1) {
    
#if PASTIX_DISTRIBUTED
      MPI_Barrier(MPI_COMM_WORLD);
      rtn=solvez_pastix_parallel(slv, RHS,transposed, verbose, debug);
#else
      rtn=solvez_pastix_sequential(slv, RHS,transposed, verbose, debug);
#endif    
      }
    else {
      rtn=solvez_pastix_sequential(slv, RHS,transposed, verbose, debug);
      }
    }
  else if  (strcmp(solver_name,"PASTIX-SEQUENTIAL") == 0) {
    rtn=solvez_pastix_sequential(slv, RHS,transposed, verbose, debug);
    }
  else {
    printf("solvez, unknown solver : %s \n",solver_name);
    }
     
#ifdef  HAVE_MPI
  if(MPI_Init_thread_DONE==1 && verbose==1) {
//   if(MPI_Init_thread_DONE==1) {
    cputime2 = MPI_Wtime();
    if (proc_id == 0 && verbose==1) {
      printf("--------------------------------------------\n");
      printf("%s: %s solving ellapsed time= %e\n", __FUNCTION__, solver_name, cputime2-cputime1);
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

  int solvez_umfpack(solver_t *slv,complex<double> *RHS,int transposed)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,bcl;
  status=-1;
#ifdef UMFPACK
  UMFPACK_STRUC_C *id;
  char *solver_name;
//  double *Info, *Control, *Ax, *tmpr, *tmpz;
  double *Info, *Control, *tmpr, *tmpz;
  complex<double> *Ax;
  complex<double> *chk;
  double *Az, *Ar, *rhsr, *rhsz;
  int *Ap, *Ai,n,nnz;
  void *Symbolic, *Numeric ;
  tripletz *Mat;
  int INCXY;
  int verbose=0;
 
  status = 0;
  solver_name=slv->name;
  if (strcmp(solver_name,"UMFPACK") == 0) {
    id = (UMFPACK_STRUC_C *) slv->parameters;
    Symbolic=id->Symbolic;
    Numeric=id->Numeric;
    Info=id->Info;
    Control=id->Control;
    Mat=(tripletz *) slv->Mat;
    Ap=Mat->i;
    Ai=Mat->j;
    Ax=Mat->x;
    if(verbose==1) printf("complex size = %d\n",sizeof(complex<double>));
    n=Mat->nrow;
    nnz=Mat->nnz;
    tmpr = (double *) malloc (n * sizeof (double)) ;
    tmpz = (double *) malloc (n * sizeof (double)) ;
    rhsr = (double *) malloc (n * sizeof (double)) ;
    rhsz = (double *) malloc (n * sizeof (double)) ;
    Az = (double *) malloc (nnz * sizeof (double)) ;
    Ar = (double *) malloc (nnz * sizeof (double)) ;
    for(bcl=0;bcl<nnz;bcl++) {
      Ar[bcl] = real(Ax[bcl]);
      Az[bcl] = imag(Ax[bcl]);
      }
    for(bcl=0;bcl<n;bcl++) {
      rhsr[bcl] = real(RHS[bcl]);
      rhsz[bcl] = imag(RHS[bcl]);
      }

    if (transposed == 1) 
      status = umfpack_zi_solve (UMFPACK_At, Ap, Ai, Ar, Az, tmpr, tmpz, rhsr, rhsz, Numeric, Control, Info) ;
    else 
      status = umfpack_zi_solve (UMFPACK_A, Ap,  Ai, Ar, Az, tmpr, tmpz, rhsr, rhsz, Numeric, Control, Info) ;
    
    for(bcl=0;bcl<n;bcl++)  RHS[bcl] = tmpr[bcl] + tmpz[bcl] * _Complex_I; 
    
    if (status < 0) {
      umfpack_zi_report_info (Control, Info) ;
      umfpack_zi_report_status (Control, status) ;
        /*error ("umfpack_di_solve failed") ;*/
      }
    free(tmpr);
    free(tmpz);
    free(rhsr);
    free(rhsz);
    free(Az);
    free(Ar);
    return(status);
  }
#endif
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solvez_spdomestic(solver_t *slv,complex<double> *RHS,int transposed)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0,bcl,n;
  tripletz *Mat;
  cs_ci *A;
  cs_cis *S ;
  cs_cin *N ;
  complex<double> *x;

  Mat = (tripletz *) slv->Mat;
  A = Mat->A;
  S = Mat->S;
  N = Mat->N;
  n = A->n ;
  x = (complex<double> *) cs_ci_malloc (n, sizeof (complex<double>)) ;    /* get workspace */

  bcl = N->L->nz;

  cs_ci_ipvec (N->pinv, RHS, x, n) ;            /* x = b(p) */
  if (transposed) {
    cs_ci_utsolve (N->U, x) ;                   /* x = U\x */
    cs_ci_ltsolve (N->L, x) ;                   /* x = L\x */
  }
  else {
    cs_ci_lsolve (N->L, x) ;                    /* x = L\x */
    cs_ci_usolve (N->U, x) ;                    /* x = U\x */
  }

  cs_ci_ipvec (S->q, x, RHS, n) ;               /* b(q) = x */
  cs_ci_free (x) ;
  

  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solvez_hips(solver_t *slv, complex<double>  *RHS, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#ifdef HIPS
#ifdef HAVE_MPI
  zHIPS_STRUC_C *id;
  char *solver_name;
  int proc_id, myid=-1, ierr, status, n, nnz,bcl, nglob;
  double cput1, cput2, cput3;
  INTS id_hips;
  INTS *unknownlist;
  char filename[50];
  FILE *out;
  complex<double> *x,*x2;
  
  cput1 = MPI_Wtime();
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

  if (debug && proc_id == 0) {printf("cpu=%d, %s starts... \n",proc_id,__FUNCTION__);}
  
  id = (zHIPS_STRUC_C *) slv->parameters;
  id_hips =id->id_hips;
  n=id->n;

  if (debug) printf("cpu=%d, %s resolution du system numero %d, n=%d \n",proc_id,__FUNCTION__,id_hips,n);

  nnz=id->nnz;
  unknownlist = id->unknownlist;
   
  /** TEst set RHS like x=1 */
//   x  = (complex<double> *) malloc(sizeof(complex<double>)*nnz);
//   x2 = (complex<double> *) malloc(sizeof(complex<double>)*nnz);
//   for (bcl=0;bcl < nnz;bcl++) x[bcl]=1.0; */
// 
//   zHIPS_MatrixVectorProduct (id_hips, x, x2);
//   for (bcl=0;bcl < n;bcl++) RHS[bcl]=x2[bcl];

  int mode=slv->RHS_distribution; 
  switch(mode) {
/*------------------------------------------------------------------------------
    version globale*/
    case SOLVER_CENTRALIZED:
      ierr = zHIPS_SetGlobalRHS(id_hips, RHS, 0,  zHIPS_ASSEMBLY_OVW);
      break;
    case SOLVER_DISTRIBUTED:
/*------------------------------------------------------------------------------
      version locale */
      ierr = zHIPS_SetRHS (id_hips, n, unknownlist, RHS,  zHIPS_ASSEMBLY_OVW, zHIPS_ASSEMBLY_OVW, zHIPS_ASSEMBLY_FOOL);
//       ierr = zHIPS_SetRHS (id_hips, n, unknownlist, RHS,  zHIPS_ASSEMBLY_OVW, zHIPS_ASSEMBLY_OVW, zHIPS_ASSEMBLY_RESPECT);
      break;
      
    }

  if (debug) printf("cpu=%d, %s zHIPS_SetRHS OK.... \n",proc_id,__FUNCTION__);

  cput1 = MPI_Wtime();
  
// /*------------------------------------------------------------------------------
//   version locale */
// //   ierr = zHIPS_GetSolution(id_hips,n,unknownlist, RHS, zHIPS_ASSEMBLY_FOOL);
// 
// /*------------------------------------------------------------------------------
//   version globale*/
//   ierr = zHIPS_GetGlobalSolution(id_hips,RHS,-1);
  
  mode=slv->solution_distribution; 
  switch(mode) {
    case SOLVER_CENTRALIZED:
      ierr = zHIPS_GetGlobalSolution(id_hips,RHS,-1);
      dHIPS_ExitOnError(ierr);
      break;
    case SOLVER_DISTRIBUTED:
      ierr = zHIPS_GetSolution(id_hips,n,unknownlist, RHS, zHIPS_ASSEMBLY_FOOL);
      dHIPS_ExitOnError(ierr);
      break;      
    }
  
//   zHIPS_SetOptionINT (id_hips,zHIPS_DISABLE_PRECOND, 1);

  cput2 = MPI_Wtime(); 
  cput3 = cput2-cput1;
  if (verbose == 1 && proc_id == 0) {
    printf("%s:  HIPS Resolution time =%lf \n",__FUNCTION__,cput3);
    }
    
//    ierr = MPI_Allreduce (  &n, &nglob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );   
//    x  = (complex<double> *) malloc(sizeof(complex<double>)*nglob);
//    ierr = zHIPS_GetGlobalSolution(id_hips,x,-1);
//    sprintf(filename,"solver_solve_hips.GlobalSolution.%d.txt",proc_id);
//    out=fopen(filename,"w");
//    for (bcl=0; bcl < nglob; bcl++) fprintf(out,"%d => %lf\n",bcl,x[bcl]);
//    fclose(out); 
//    free(x);
//    RHS = x;

  ierr=0;
  return(ierr);
#endif
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solvez_pastix_sequential(solver_t *slv,complex<double> *RHS,int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
#ifdef PASTIX
  PASTIX_STRUC_C *id;
  char *solver_name;
  complex<double> *Ax;
  int *Ap, *Ai,n,nnz;
  tripletz *Mat;
  int INCXY;
  char filerhs[50];
  int myid,sz,i,ierr;
  FILE *out;
 
  status = 0;
  solver_name=slv->name;

  id = (PASTIX_STRUC_C *) slv->parameters;
  Mat=(tripletz *) slv->Mat;
  Ap=Mat->i;
  Ai=Mat->j;
  Ax=Mat->x;
  n=Mat->nrow;
  nnz=Mat->nnz;

  if(verbose==1) printf("%s : n=%d nnz=%d\n",__FUNCTION__, n,nnz);
  
  id->iparm[IPARM_START_TASK] = API_TASK_SOLVE;
  id->iparm[IPARM_END_TASK]   = API_TASK_REFINE;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  in case where transpose matrix has been provided instead of direct one, then
  swap transposed solution flag
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(slv->transpose_matrix==1) transposed=1-transposed;
  
  id->iparm[IPARM_TRANSPOSE_SOLVE] = transposed;

  z_pastix(&(id->pastix_data), MPI_COMM_WORLD,
         n, Ap, Ai, Ax,
         id->perm, id->invp, RHS, 1, id->iparm, id->dparm);
      
#endif
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solvez_pastix_parallel(solver_t *slv, complex<double> *RHS, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,bcl,bcl2;
  int ierr=-1;
#ifdef PASTIX
#ifdef HAVE_MPI
  PASTIX_STRUC_C *id;
  char *solver_name;
  complex<double> *Ax;
  int *Ap, *Ai,n,nnz;
  tripletz *Mat;
  int INCXY;
  int rank,sz,i;
  int *loc2glb;
  
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &sz);
  
  status = 0;
  solver_name=slv->name;

  id = (PASTIX_STRUC_C *) slv->parameters;
  Mat=(tripletz *) slv->Mat;
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

  z_dpastix(&(id->pastix_data), MPI_COMM_WORLD,
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


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solvez_mumps(solver_t *slv, complex<double> *RHS, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#ifdef MUMPS
  ZMUMPS_STRUC_C *id;
  char *solver_name;
  int ierr, rank;
  double cput1, cput2, cput3;
  int status;
  
#ifdef HAVE_MPI
  cput1 = MPI_Wtime();
#endif
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  solver_name=slv->name;

  id = (ZMUMPS_STRUC_C *) slv->parameters;
  
  if (transposed == 1) {
    (*id).ICNTL(9)=2;
    }
    
  if (rank == 0) {
    (*id).rhs = (mumps_double_complex *) RHS;
    }
  (*id).job=3;
  
  zmumps_c(id);

#ifdef HAVE_MPI
  cput2 = MPI_Wtime(); 
  cput3 = cput2-cput1;
  if (rank == 0) {
    if(verbose==1) printf(" ZMUMPS resolution time =%lf\n",cput3);
    }
#endif
  ierr=(*id).info[0];
  if(ierr==-1) {
    printf("cpu %d: ZMUMPS solver failed, status=%d due to cpu %d\n",rank, ierr, id->info[1]);
    }
  if(ierr==-22) {
    printf("cpu %d: ZMUMPS solver failed, status=%d array fault %d\n",rank, ierr, id->info[1]);
    }
  return(ierr);
#endif
}

#endif

