
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
extern int HIPSSOLVERNBZ;
extern int nbsyslinearz2;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorizez(solver_t *slv,tripletz *Mat, int verbose)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
   char *solver_name;
   int rtn,base,format;
   int status;
   double default_precision=1.e-08;
   double cputime1, cputime2;
   int proc_id, ncpus, ierr;  
   bool debug =false;
 
#ifdef  HAVE_MPI
  if(MPI_Init_thread_DONE==1) {
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &ncpus);
    cputime1 = MPI_Wtime();
    }
#endif      
   
   solver_name=slv->name;
   
   if(verbose==1) printf("%s complex-matrix factorization : %s \n",__FUNCTION__,solver_name);

   if (strcmp(solver_name,"MUMPS") == 0) {
     base=1;
     format=COO;
     rtn = convertz(Mat, format, base, verbose);
     if(rtn!=0) return(-1);
     rtn = factorizez_mumps(slv, Mat, verbose, false);
     }
   else if  (strcmp(solver_name,"UMFPACK") == 0) {
     base=0;
     format=CSC;
     status=convertz(Mat,format,base, verbose);
     if(status!=0) return(-1);
     rtn=factorizez_umfpack(slv, Mat, verbose);
     }
   else if  (strcmp(solver_name,"LAPACK") == 0) {
     base=1;
     /*convertz(Mat,3,base);*/
     /*rtn=factorizez_lapack(slv, Mat);*/
     }
   else if  (strcmp(solver_name,"SpDOMESTIC") == 0) {
     base=0;
     format=CSC;
     status=convertz(Mat,format,base, verbose);
     if(status!=0) return(-1);
     rtn=factorizez_spdomestic(slv, Mat);
     }
  else if (strcmp(solver_name,"HIPS") == 0){
    base=1;
    format=COO;
    rtn = convertz(Mat, format, base, verbose);
    if(rtn!=0) return(-1);
#ifdef __DEBUG_HIPS
    debug=true;
#endif
    rtn=factorizez_hips(slv, Mat, verbose, debug);
//     rtn=factorizez_hips_complex(slv, Mat, verbose, debug);
    }
  else if (strcmp(solver_name,"PASTIX") == 0){
    if(ncpus>1) {
    
#if PASTIX_DISTRIBUTED
      base=1;
      verbose=0;
      format=CSC;
      rtn = convertz(Mat, format, base, verbose);
      rtn=factorizez_pastix_parallel_distributed(slv, Mat, verbose, debug);    
#else
      // Attention en mode parallel la matrix arrive en FORMAT COO
      // avec des indices globaux...
      // la redistibution, conversion, ...
      // Les jobs sont faits dans factorize_paxtix_parallel 
      base=1;
//       format=COO;
      format=COO_COL_ARRANGED;
      rtn = convertz(Mat, format, base, verbose);
      rtn=factorizez_pastix_parallel(slv, Mat, verbose, debug);
#endif
      }
    else {
      base=0; // this is wrong but was corrected through the call to pastix_checkMatrix
      format=CSC;
      status=convertz(Mat, format, base, verbose);
      rtn=factorizez_pastix_sequential(slv, Mat, verbose, debug);
      }
    }
  else if (strcmp(solver_name,"PASTIX-SEQUENTIAL") == 0){
    base=0; // this is wrong but was corrected through the call to pastix_checkMatrix
    format=CSC;
    status=convertz(Mat, format, base, verbose);
    if(status!=0) return(-1);
    rtn=factorizez_pastix_sequential(slv, Mat, verbose, debug);
    }
  else {
    printf("factorizez, unknown solver : %s \n",solver_name);
    rtn = -1;
    }
     
#ifdef  HAVE_MPI
   if(MPI_Init_thread_DONE==1 && verbose==1) {
     cputime2 = MPI_Wtime();
     if (verbose==1 && proc_id == 0) {
       printf("--------------------------------------------\n");
       printf("%s: %s factorization ellapsed time= %e\n",__FUNCTION__, solver_name, cputime2-cputime1);
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

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorizez_spdomestic(solver_t *slv, tripletz *Mat)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0,n,order=3,bcl;
  double tol=1.0; /* unsymetrique usage*/
  
  cs_ci *A = (cs_ci *) cs_calloc (1, sizeof (cs_ci)) ;
  cs_cis *S ;
  cs_cin *N ;
  
  Mat->A = A;
  A->nzmax = Mat->nnz;
  A->m =  Mat->nrow;
  A->n =  Mat->ncol;
  A->p =  Mat->i;
  A->i =  Mat->j;
  A->x =  Mat->x;
  A->nz=-1;
  n = A->n;

  S = (cs_cis *) cs_ci_sqr (order, A, 0) ;      /* ordering and symbolic analysis */

  N = (cs_cin *) cs_ci_lu (A, S, tol) ;         /* numeric LU factorization       */
  
  bcl = N->L->nz;
  Mat->S = S;
  Mat->N = N;
  slv->Mat = Mat;
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorizez_umfpack(solver_t *slv, tripletz *Mat, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
#ifdef UMFPACK
  UMFPACK_STRUC_C *id;
  char *solver_name;
  double *Info, *Control;
  complex<double> *Ax;
  double *Az, *Ar;
  int *Ap, *Ai,n,nnz,bcl;
  void *Symbolic, *Numeric ;
  
  solver_name=slv->name;
  if (strcmp(solver_name,"UMFPACK") == 0) {
    if(verbose==1) printf("factorize UMFPACK ...\n");
    id = (UMFPACK_STRUC_C *) slv->parameters;
    slv->Mat = Mat;
    Symbolic=id->Symbolic;
    Numeric=id->Numeric;
    Info=id->Info;
    Control=id->Control;
    Ap=Mat->i;
    Ai=Mat->j;
    Ax=Mat->x;
    n=Mat->nrow;
    nnz=Mat->nnz;
    
    Az = (double *) malloc (nnz * sizeof (double)) ;
    Ar = (double *) malloc (nnz * sizeof (double)) ;
    for(bcl=0;bcl<nnz;bcl++) {
      Ar[bcl] = real(Ax[bcl]);
      Az[bcl] = imag(Ax[bcl]);
      }

    status = umfpack_zi_symbolic (n, n, Ap, Ai, Ar, Az, &Symbolic, Control, Info) ; 

    if (status < 0) {
      printf("umfpack_zi_symbolic status=%d\n",status);
      umfpack_zi_report_info (Control, Info) ;
      umfpack_zi_report_status (Control, status) ;
      return(-1);
      /*error ("umfpack_di_symbolic failed") ;*/
      }

    status = umfpack_zi_numeric (Ap, Ai, Ar, Az, Symbolic, &Numeric, Control, Info) ;
    id->Symbolic=Symbolic;

    id->Numeric=Numeric; 
    if (status < 0) {
      printf("umfpack_zi_numeric status=%d\n",status);
      umfpack_zi_report_info (Control, Info) ;
      umfpack_zi_report_status (Control, status) ;
      return(-1);
      /*error ("umfpack_di_numeric failed") ;*/
      }
    free(Ar);
    free(Az);
    return(status);
  }
#else
  status=-1;
  printf("Please compile poc-solvers library with -DUMFPACK \n");
  return(status);
#endif
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorizez_mumps(solver_t *slv, tripletz *Mat, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#ifdef MUMPS
  ZMUMPS_STRUC_C *id;
  int myid, ierr, nglob=0, count, rank, size;
  char *solver_name;
  
#ifdef HAVE_MPI
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
  ierr = 0;
#endif
  
//   solver_name=slv->name;
//   if (strcmp(solver_name,"MUMPS") == 0) {
//     if(debug) printf("cpu=%d, factorizez MUMPS ...\n",myid);
  id = (ZMUMPS_STRUC_C *) slv->parameters;
  (*id).job=4;
  
  /* get global row number */
  count=1;
#ifdef HAVE_MPI
  ierr = MPI_Allreduce (  &(Mat->nrow), &nglob, count, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  (*id).n = nglob;
#else
  (*id).n = Mat->nrow;
  ierr = 0;
#endif
  
  (*id).ICNTL(18)=3;
/*------------------------------------------------------------------------------
    5.1.1 and later */
    (*id).nnz_loc = Mat->nnz;
//     (*id).nz_loc = Mat->nnz;
  (*id).irn_loc= Mat->i;  
  (*id).jcn_loc= Mat->j; 
  (*id).a_loc = (mumps_double_complex *) Mat->x;  

  (*id).ICNTL(1)=-1; 
  /*(*id).ICNTL(1)=6; */
  (*id).ICNTL(2)=-1; 
  /*(*id).ICNTL(2)=0; */
  (*id).ICNTL(3)=-1; 
  /*(*id).ICNTL(3)=6; */
  (*id).ICNTL(4)=4;
  (*id).ICNTL(14)=20;
  
/*------------------------------------------------------------------------------
  dense RHS  */
  (*id).ICNTL(20)=0;
/*------------------------------------------------------------------------------
  centralized solution  */
  (*id).ICNTL(21)=0;

  zmumps_c(id);
  
  ierr=(*id).info[0];
  if(ierr==-1)  printf("cpu %d: ZMUMPS factorisation failed, status=%d due to cpu %d\n",rank, ierr, id->info[1]);
  if(ierr==-9)  printf("cpu %d: ZMUMPS factorisation failed, status=%d due to memory deficiency %d\n",rank, ierr, id->info[1]);
  if(ierr==-10) printf("cpu %d: ZMUMPS factorisation failed, status=%d due to singular matrix\n",rank, ierr);
  if(ierr==-13) printf("cpu %d: DMUMPS factorisation failed, status=%d due to memory deficience %d\n",rank, ierr, id->info[2]);
  if(ierr==0) {
    if(verbose==1) {
      printf("%s, cpu=%d: MEMORY USE, %d used by cpu, %d total Mo, status=%d\n",__FUNCTION__,rank, (*id).infog[20],(*id).infog[21], ierr);
      }
    if(verbose==1 && rank==0) {
      printf("MEMORY SCALING DIAG %d %d\n",size,(*id).infog[21]);
      }
    }
  return(ierr);
  
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorizez_hips(solver_t *slv,tripletz *Mat, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,ierr;
#ifdef HIPS
#ifdef HAVE_MPI
  zHIPS_STRUC_C *id;
  INTS  sym_pattern, sym_matrix;
  INTS id_hips, idnbr, i, j;
  INTS *unknownlist, *hipsnodelist;
  COEFZ *x, *rhsloc;
  INTS proc_id, n, ln;
  INTL *ia, nnz;
  INTS *ja;
  COEFZ *a;
  INTS domsize, nproc;
  INTS pbegin, pend;
  INTS *mapptr,  *mapptr2,*mapp,*mapp2;
  INTS p;

  int nglob, nnzglob, tmp, bcl, idi, info;
  FILE *out;
  char filename[50];
  complex_t *x1,*x2;
  double tol;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  
  id = (zHIPS_STRUC_C *) slv->parameters; 
  
  /* nglob */
  tmp=Mat->nrow;
  
  if (debug && proc_id == 0) {printf("cpu=%d,  %s nloc=%d\n",proc_id,__FUNCTION__,tmp);}
  
  MPI_Allreduce( &tmp, &nglob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  
  if (debug && proc_id == 0) printf("cpu=%d,  %s nglob=%d\n",proc_id,__FUNCTION__,nglob);

//   idnbr = 100; /* total */
//   if (HIPSSOLVERNBZ == 0) {
//     ierr = zHIPS_Initialize(idnbr);
//     zHIPS_ExitOnError(ierr);
//     }
//   
//   id_hips = HIPSSOLVERNBZ; /** id_hips of the linear system **/
//   HIPSSOLVERNBZ++;
   
  id_hips = id->id_hips;
  
  status = MPI_Barrier(MPI_COMM_WORLD);
   
  if (debug && proc_id == 0) {
    printf("cpu=%d,  %s solver numero=%d\n",proc_id,__FUNCTION__,id_hips);
    }
    
  status = MPI_Barrier(MPI_COMM_WORLD);

  domsize = nglob/nproc/2;
//   domsize = nglob/(nproc+1);

  if (debug && proc_id == 0) {
    printf("cpu=%d,  %s domsize=%d\n",proc_id,__FUNCTION__,domsize);
    }
    
  status = MPI_Barrier(MPI_COMM_WORLD);

  zHIPS_SetOptionINT(id_hips, zHIPS_PARTITION_TYPE, 0); 
  zHIPS_SetOptionINT(id_hips, zHIPS_DOMSIZE, domsize);
  
  status = MPI_Barrier(MPI_COMM_WORLD);
  
#if 0
  zHIPS_SetDefaultOptions (id_hips, zHIPS_HYBRID);
//   zHIPS_SetDefaultOptions (id_hips, zHIPS_DIRECT);
#else
  zHIPS_SetDefaultOptions (id_hips, zHIPS_ITERATIVE);
  
  tol=1e-12;
  zHIPS_SetOptionREAL(id_hips,zHIPS_DROPTOL0,tol);
  tol=1e-12;
  zHIPS_SetOptionREAL(id_hips,zHIPS_DROPTOL1,tol);
//   tol=1e-12;
//   zHIPS_SetOptionREAL(id_hips,zHIPS_DROPTOLE,tol);
   
//   zHIPS_SetOptionINT(id_hips, zHIPS_ITMAX,150);
//   zHIPS_SetOptionINT(id_hips, zHIPS_ITMAX_SCHUR,150);
#endif
  
  status = MPI_Barrier(MPI_COMM_WORLD);
   
  if(verbose==1)  info=4;
  else            info=0;

  zHIPS_SetOptionINT (id_hips, zHIPS_VERBOSE, info);

  sym_matrix=0;
  zHIPS_SetOptionINT(id_hips, zHIPS_SYMMETRIC, sym_matrix);
  zHIPS_SetOptionINT(id_hips, zHIPS_GRAPH_SYM,0);
  
  /** C : numbering starts from 0 **/
  zHIPS_SetOptionINT(id_hips, zHIPS_FORTRAN_NUMBERING, 0);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  ENTER THE GRAPH
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  status = MPI_Barrier(MPI_COMM_WORLD);
  
  ierr = zHIPS_GraphBegin(id_hips,nglob, Mat->nnz);
  zHIPS_ExitOnError(ierr);
  for (bcl=0;bcl<Mat->nnz;bcl++) {
    ierr = zHIPS_GraphEdge(id_hips, Mat->i[bcl]-1, Mat->j[bcl]-1);
    }
  ierr = zHIPS_GraphEnd(id_hips);
  zHIPS_ExitOnError(ierr);
  status = MPI_Barrier(MPI_COMM_WORLD);
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   ENTER A USER PARTITION
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
/* Only the master processor (0 by default) needs to enter the partition
   collect inforamtion on proc 0 */
  mapptr = (INTS *)malloc(sizeof(INTS)*(nproc+1));
  mapptr2= (INTS *)malloc(sizeof(INTS)*(nproc+1));
  
  mapp =   (INTS *)malloc(sizeof(INTS)*(nglob+1));
  mapp2=   (INTS *)malloc(sizeof(INTS)*(nglob+1));

  for (bcl=0;bcl<nproc+1;bcl++)   mapptr2[bcl]=0;
  for (bcl=0;bcl<nglob+1;bcl++)   mapp2[bcl]  =0;
  
  for (bcl=0;bcl<nproc+1;bcl++)   mapptr[bcl] =0;
  for (bcl=0;bcl<nglob+1;bcl++)   mapp[bcl]   =0;
  
  for (bcl=proc_id+1;bcl<nproc+1;bcl++) mapptr2[bcl]=Mat->nrow;

  MPI_Allreduce( mapptr2, mapptr, nproc+1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  
  idi=0;
  mapp2[idi+mapptr[proc_id]] =  Mat->i[0]-1;
  for (bcl=1; bcl < Mat->nnz; bcl++) {
    if ( mapp2[idi+mapptr[proc_id]] != (Mat->i[bcl]-1) ) {
      idi++;
      mapp2[idi+mapptr[proc_id]] = Mat->i[bcl]-1;    
      }
    }
    
//   status = MPI_Barrier(MPI_COMM_WORLD);
//   
//   printf("cpu=%d ; ndofs %d %d\n", proc_id, idi+1,Mat->nrow);
  
  MPI_Allreduce( mapp2, mapp, nglob, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Set the unknownlist 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  unknownlist= (INTS *)malloc(sizeof(INTS)*(Mat->nrow));
  unknownlist[0]=Mat->i[0]-1;
  idi=0;
  for (bcl=1; bcl < Mat->nnz; bcl++) {
    if ( unknownlist[idi] != (Mat->i[bcl]-1) ) {
      idi++;
      unknownlist[idi] = Mat->i[bcl]-1;    
      }
    }
  
  status = MPI_Barrier(MPI_COMM_WORLD);
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Enter Partition 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status = MPI_Barrier(MPI_COMM_WORLD);
  if (proc_id ==0) {
    ierr = zHIPS_SetPartition(id_hips, nproc, mapptr, mapp);
    zHIPS_ExitOnError(ierr);
    }
  status = MPI_Barrier(MPI_COMM_WORLD);
    
  free(mapptr2);
  free(mapptr);
  free(mapp2);
  free(mapp);
  
/***************************************************/
/*            GET THE LOCAL UNKNOWN LIST           */
/***************************************************/
//   ierr = zHIPS_GetLocalUnknownNbr(id_hips, &ln);
//   zHIPS_ExitOnError(ierr);
//   printf("cpu=%d, %s nombre de noeuds interne %d\n",proc_id,__FUNCTION__,ln);
//   hipsnodelist = (INTS *)malloc(sizeof(INTS)*ln);
//   ierr = zHIPS_GetLocalUnknownList(id_hips, hipsnodelist);
//   zHIPS_ExitOnError(ierr);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  ENTER THE MATRIX COEFFICIENT
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/* The processor enter the rows pbegin to pend-1 (naive partition) 
   In the zHIPS_ASSEMBLY_FOOL mode any processor can enter any coefficient :
   nevertheless it is better to enter as much as possible coefficients
   that are in the local matrix of the processor */
  
    
  status = MPI_Barrier(MPI_COMM_WORLD);
  
  if (debug) {
    printf("cpu=%d,  %s now begin assembly\n",proc_id,__FUNCTION__);
    }
  status = MPI_Barrier(MPI_COMM_WORLD);
  
  ierr = zHIPS_AssemblyBegin(id_hips, Mat->nnz, zHIPS_ASSEMBLY_OVW, zHIPS_ASSEMBLY_OVW, zHIPS_ASSEMBLY_FOOL, sym_matrix);
  zHIPS_ExitOnError(ierr);
  
  for(bcl=0;bcl<Mat->nnz;bcl++) {
    ierr = zHIPS_AssemblySetValue(id_hips, Mat->i[bcl]-1, Mat->j[bcl]-1, Mat->x[bcl]);
    zHIPS_ExitOnError(ierr);
    }
  ierr = zHIPS_AssemblyEnd(id_hips);
  zHIPS_ExitOnError(ierr);
    
  status = MPI_Barrier(MPI_COMM_WORLD);
  
  if (debug) {
    printf("cpu=%d,  %s now set assembly done\n",proc_id,__FUNCTION__);
    }
  status = MPI_Barrier(MPI_COMM_WORLD);
  

  id->id_hips = id_hips;
  id->unknownlist = unknownlist;
  id->n = Mat->nrow;
  id->nnz = Mat->nnz;

  if (debug && proc_id == 0) printf("cpu=%d, sortie %s \n",proc_id,__FUNCTION__);

  status = MPI_Barrier(MPI_COMM_WORLD);

  ierr=0;
  status=0;
  
  return(status);
  
#endif 
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorizez_pastix_sequential(solver_t *slv, tripletz *Mat, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
#ifdef PASTIX
  PASTIX_STRUC_C *id;
  char *solver_name;
  complex<double> *Ax;
  pastix_int_t *Ap, *Ai;
  int  n,nnz,nbthread;
  complex<double> *rhs = NULL; /* right hand side                      */
  int    verbosemode;         /* Level of verbose mode (0, 1, 2)      */
  char   *type        = NULL; /* type of the matrix                   */
  int    bcl, bcl2;
  int    flagsym=0;
  pastix_int_t ncol;          /* Size of the matrix                   */
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
    
  id = (PASTIX_STRUC_C *) slv->parameters;
  slv->Mat = Mat; 
  Ap=Mat->i;
  Ai=Mat->j;
  Ax=Mat->x;
  n=Mat->nrow;
  nnz=Mat->nnz;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Initialize parameters to default values
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  verbosemode=verbose;  
  
  z_pastix_checkMatrix(MPI_COMM_WORLD, verbosemode,
                       API_SYM_NO, API_YES, n, &Ap, &Ai, &Ax, NULL,1);

  Mat->x=Ax;
  Mat->i=Ap;
  Mat->j=Ai;
  
  id->iparm[IPARM_MODIFY_PARAMETER] = API_NO;
  
  z_pastix(&(id->pastix_data), MPI_COMM_WORLD, n, Ap, Ai, Ax,
             id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Customize some parameters 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
#ifdef OMP_H
  nbthread = omp_get_max_threads(); 
  if(debug) printf("%s : using %d threads (omp_get_max_threads)...................... \n",__func__,nbthread);
//   nbthread = omp_get_num_threads();
//   printf("%s : using %d threads (omp_get_num_threads)...................... \n",__func__,nbthread);
#else
  nbthread = 1;  
  if(debug) printf("%s : using only 1 thread!...................... \n",__func__,nbthread);
#endif
  
  id->iparm[IPARM_SYM]           = API_SYM_NO;
  id->iparm[IPARM_FACTORIZATION] = API_FACT_LU;
  id->iparm[IPARM_THREAD_NBR]    = nbthread;
  id->iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
  id->iparm[IPARM_LEVEL_OF_FILL] = 0;
  id->iparm[IPARM_START_TASK]    = API_TASK_ORDERING;
  id->iparm[IPARM_END_TASK]      = API_TASK_NUMFACT;
//   id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-8;
  id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-12;
  id->iparm[IPARM_VERBOSE]            = verbosemode;
  /*iparm[IPARM_END_TASK]   = API_TASK_REFINE;*/

//   id->perm = malloc(n*sizeof(pastix_int_t));
//   id->invp = malloc(n*sizeof(pastix_int_t));
  id->perm = new pastix_int_t[n];
  id->invp = new pastix_int_t[n];
  

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Factorisation
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  z_pastix(&(id->pastix_data), MPI_COMM_WORLD,  n, Ap, Ai, Ax,
         id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
  
  status=id->iparm[IPARM_ERROR_NUMBER];
  
  if (status < 0)  {
    printf("%s factorization status=%d\n",__FUNCTION__,status);  
    return(-1);
    }
 
  if(verbose==1) printf("Factorizez pastix done (%s)...\n",__FUNCTION__);
  status=0;        
  return(status);
  
#else
  status=-1;
  printf("%s (line %d) : please compile poc-solvers library with -DPASTIX \n", __func__,__LINE__);
  return(status);
#endif 
  return(status); 
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorizez_pastix_parallel_distributed(solver_t *slv, tripletz *Mat, int verbose, bool debug) {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  int k,status;
#ifdef PASTIX
#ifdef HAVE_MPI
  PASTIX_STRUC_C *id;
  char *solver_name;
  complex<double> *Ax;
  pastix_int_t *Ap, *Ai;
  int  n,nnz,nbthread;
  complex<double> *rhs = NULL; /* right hand side                      */
  int    verbosemode;         /* Level of verbose mode (0, 1, 2)      */
  char   *type        = NULL; /* type of the matrix                   */
  int    bcl, bcl2;
  int    flagsym=0;
  pastix_int_t    ncol;       /* Size of the matrix                   */
  int   mpid,sz;
  pastix_int_t *loc2glb;

  MPI_Barrier(MPI_COMM_WORLD);
      
  id = (PASTIX_STRUC_C *) slv->parameters;
  
  slv->Mat = Mat;
  
  Ap=Mat->i;
  Ai=Mat->j;
  Ax=Mat->x;
  n=Mat->nrow;

  nnz=Mat->nnz;

//   loc2glb=malloc(n*sizeof(pastix_int_t));
  loc2glb=new pastix_int_t[n];

  for(k=0;k<n;k++) loc2glb[k]=Mat->loc2glob[k];

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Comm_rank(MPI_COMM_WORLD, &mpid);
  MPI_Comm_size(MPI_COMM_WORLD, &sz);
  
   
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Initialize parameters to default values
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  verbosemode=verbose;  
  
  z_pastix_checkMatrix(MPI_COMM_WORLD, verbosemode, API_SYM_NO, API_YES,
                      n, &Ap, &Ai, &Ax, &loc2glb, 1);
  
/*------------------------------------------------------------------------------
  get (possibly changed) new addresses */
  Mat->i=Ap;
  Mat->j=Ai;
  Mat->x=Ax;
  
  id->iparm[IPARM_MODIFY_PARAMETER] = API_NO;
  
  z_dpastix(&(id->pastix_data), MPI_COMM_WORLD, n, Ap, Ai, Ax, loc2glb,
              id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Customize some parameters
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
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
  
  MPI_Barrier(MPI_COMM_WORLD);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Factorisation
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(verbose==1) printf("%s pastix factorisation start\n",__FUNCTION__);
  
  z_dpastix(&(id->pastix_data), MPI_COMM_WORLD, n, Ap, Ai, Ax, loc2glb,
              id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
  
  if(verbose==1) printf("%s pastix factorisation end\n",__FUNCTION__);
  
  status=id->iparm[IPARM_ERROR_NUMBER];

  if (status < 0)  { 
    printf("z_pastix factorization status=%d\n",status);  
    return(-1);
    }
 
  status=0;        
  return(status);
#endif 
#else
  status=-1;
  printf("%s (line %d) : please compile poc-solvers library with -DPASTIX \n", __func__,__LINE__);
  return(status);
#endif 
  return(status); 
}

// #endif //moved at ebd

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorizez_pastix_parallel(solver_t *slv, tripletz *Mat, int verbose, bool debug) {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  int status;
#ifdef PASTIX
#ifdef HAVE_MPI
  PASTIX_STRUC_C *id;
  char *solver_name;
  complex<double> *Ax;
  pastix_int_t *Ap, *Ai;
  int  n,nnz,nbthread;
  complex<double> *rhs = NULL; /* right hand side                      */
  int    verbosemode;         /* Level of verbose mode (0, 1, 2)      */
  char   *type        = NULL; /* type of the matrix                   */
  int    bcl, bcl2;
  int    flagsym=0;
  pastix_int_t    ncol;       /* Size of the matrix                   */
  // Global variables
  int ncolglob, nnzglob;
  int *Apglob, *Aiglob, *colptr, *rows;
  complex<double> *Axglob, *values;
  char filemat[50];
  FILE *f_out2;
  int   mpid,sz;

  MPI_Barrier(MPI_COMM_WORLD);
      
  id = (PASTIX_STRUC_C *) slv->parameters;
  slv->Mat = Mat; 
  Ap=Mat->i;
  Ai=Mat->j;
  Ax=Mat->x;
  n=Mat->nrow;
  nnz=Mat->nnz;

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Comm_rank(MPI_COMM_WORLD, &mpid);
  MPI_Comm_size(MPI_COMM_WORLD, &sz);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Dispatch Matrix global matrix
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  z_coo2glob(n, nnz, Ap, Ai, Ax, &ncolglob, &nnzglob, &Apglob, &Aiglob, &Axglob);
  
  free(Ap);
  free(Ai);
  free(Ax);
  
  MPI_Barrier(MPI_COMM_WORLD);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  convert structure to CSC
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  z_globcoo2csc(ncolglob, nnzglob, Apglob, Aiglob, Axglob, &colptr, &rows, &values);
  
  free(Apglob);
  free(Aiglob);
  free(Axglob);
  
  MPI_Barrier(MPI_COMM_WORLD); 

  if(debug==1) {
/*------------------------------------------------------------------------------
    PRINT out */
    sprintf(filemat,"Solverz_MatrixGlobal_CSC.%d.%d.txt",nbsyslinearz2,mpid);
    MPI_Barrier(MPI_COMM_WORLD);
    printf("%d::: z_globcoo2csc filemat=%s \n",mpid, filemat);
    f_out2 = fopen(filemat,"w");
    fprintf(f_out2,"  %%%MatrixMarket matrix coordinate complex general \n");
    fprintf(f_out2,"%d  %d  %d 1 0 \n",ncolglob,ncolglob,nnzglob);
    for (bcl=0; bcl < ncolglob; bcl++) {
//       fprintf(f_out,"%d %d\n",bcl+1,colptr[bcl]); // Use Fortran Numbering 
      for (bcl2=colptr[bcl]-1;bcl2<colptr[bcl+1]-1;bcl2++) {
        fprintf(f_out2,"%d  %d  %e %e \n",bcl+1,rows[bcl2],values[bcl2]);
        }
      }
    fclose(f_out2); 
    }

//     free(Ap);
//     free(Ai);
//     free(Ax);
    
  Ap = colptr;
  Ai = rows;
  Ax = values;
  
  Mat->nrow = ncolglob;
  Mat->nnz  = nnzglob;
   
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Initialize parameters to default values
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  verbosemode=verbose;  
  
  z_pastix_checkMatrix(MPI_COMM_WORLD, verbosemode,
                     API_SYM_NO, API_YES,
                     ncolglob, &Ap, &Ai, &Ax, NULL,1);
  Mat->i=Ap;
  Mat->j=Ai;
  Mat->x=Ax;
  Mat->comptype=CSC;
  
  id->iparm[IPARM_MODIFY_PARAMETER] = API_NO;
  z_pastix(&(id->pastix_data), MPI_COMM_WORLD,
         ncolglob, Ap, Ai, Ax,
         id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Customize some parameters
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
#ifdef OMP_H
//   nbthread = omp_get_max_threads(); 
//   if(debug) printf("%s : using %d threads (omp_get_max_threads)...................... \n",__func__,nbthread);
  nbthread = omp_get_num_threads();
  if(debug) printf("%s : using %d threads (omp_get_num_threads)...................... \n",__func__,nbthread);
#else
  nbthread = 1;  
  if(debug) printf("%s : using only 1 thread!...................... \n",__func__,nbthread);
#endif
    
  verbosemode=verbose;
  
  id->iparm[IPARM_SYM]           = API_SYM_NO;
  id->iparm[IPARM_FACTORIZATION] = API_FACT_LU;
  id->iparm[IPARM_THREAD_NBR] = nbthread;
  id->iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
  id->iparm[IPARM_LEVEL_OF_FILL] = 0;
  id->iparm[IPARM_START_TASK] = API_TASK_ORDERING;
  id->iparm[IPARM_END_TASK]   = API_TASK_NUMFACT;
//   id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-8;
  id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-12;
  id->iparm[IPARM_VERBOSE]             = verbosemode;

  /*iparm[IPARM_END_TASK]   = API_TASK_REFINE;*/
//   id->perm = malloc(ncolglob*sizeof(pastix_int_t));
//   id->invp = malloc(ncolglob*sizeof(pastix_int_t) );
  id->perm = new pastix_int_t[ncolglob];
  id->invp = new pastix_int_t[ncolglob];
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Factorisation
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(verbose==1) printf("%s pastix factorisation start\n",__FUNCTION__);
  z_pastix(&(id->pastix_data), MPI_COMM_WORLD,
         ncolglob, Ap, Ai, Ax,
         id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
  if(verbose==1) printf("%s pastix factorisation end\n",__FUNCTION__);
  
  status=id->iparm[IPARM_ERROR_NUMBER];

  if (status < 0)  { 
    printf("z_pastix factorization status=%d\n",status);  
    return(-1);
    }
 
  status=0;        
  return(status);
#endif 
#else
  status=-1;
  printf("%s (line %d) : please compile poc-solvers library with -DPASTIX \n", __func__,__LINE__);
  return(status);
#endif 
  return(status); 
}

#endif //moved at ebd
