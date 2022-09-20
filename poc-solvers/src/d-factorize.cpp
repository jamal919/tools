
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
#include <complex.h>

#include "solvers-interface.h"
#include "solverlib.h"

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match MUMPS documentation */

extern int MPI_Init_thread_DONE;
extern int nbsyslinearz2;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorize(solver_t *slv, triplet *Mat, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  char *solver_name;
  int rtn,base,format;
  double default_precision=1.e-08;
  double cputime1,cputime2;
  int proc_id,ncpus,ierr;  
  bool debug=false; 

  solver_name=slv->name;
  
  slv->Mat=Mat;
  
#ifdef HAVE_MPI
   if(MPI_Init_thread_DONE==1) {
     ierr = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
     ierr = MPI_Comm_size(MPI_COMM_WORLD, &ncpus);
     cputime1 = MPI_Wtime();
     }
#endif
  
  if (strcmp(solver_name,"MUMPS") == 0) {
    base=1;
    format=COO;
    rtn = convert(Mat, format, base, verbose);
    if(rtn!=0) return(-1);
    rtn = factorize_mumps(slv, Mat, verbose, false);
    }
  else if  (strcmp(solver_name,"UMFPACK") == 0) {
    base=0;
    format=CSC;
    rtn = convert(Mat, format, base, verbose);
    if(rtn!=0) return(-1);
    rtn=factorize_umfpack(slv, Mat, verbose);
    }
  else if  (strcmp(solver_name,"HYPRE") == 0) {
    base=1;
    format=COO;
    rtn = convert(Mat, format, base, verbose);
    if(rtn!=0) return(-1);
    rtn=factorize_hypre(slv, Mat);
    }
  else if  (strcmp(solver_name,"HIPS") == 0) {
    base=1;
    format=COO;
    rtn = convert(Mat, format, base, verbose);
    if(rtn!=0) return(-1);
    if(MPI_Init_thread_DONE==0) {
      }
#ifdef __DEBUG_HIPS
    debug=true;
#endif
    rtn=factorize_hips(slv, Mat, verbose, debug);
    }
  else if  (strcmp(solver_name,"LAPACK") == 0) {
    base=1;
    format=BND; 
    rtn = convert(Mat, format, base, verbose);
    if(rtn!=0) return(-1);
    rtn=factorize_lapack(slv, Mat);
    }
  else if  (strcmp(solver_name,"SpDOMESTIC") == 0) {
    base=0;
    format=CSC;
    rtn = convert(Mat, format, base, verbose);
    if(rtn!=0) return(-1);
    rtn=factorize_spdomestic(slv, Mat);
    }
  else if  (strcmp(solver_name,"PASTIX") == 0) {
    if(ncpus>1) {
        
#if PASTIX_DISTRIBUTED
      base=1;
      verbose=0;
      format=CSC;
      rtn = convert(Mat, format, base, verbose);
      rtn=factorize_pastix_parallel_distributed(slv, Mat, verbose, debug);    
#else
    // Attention en mode parallel les matrix arrivent en FORMAT COO
    // avec des indices globaux...
    // la redistibution, conversion, ...
    // Les jobs sont faits dans factorize_paxtix_parallel
      base=1;
      format=COO_COL_ARRANGED;
      rtn = convert(Mat, format, base, verbose);
      if(debug) printf("call factorize_pastix_parallel .......;\n");
      rtn=factorize_pastix_parallel(slv, Mat, verbose, debug);
#endif
      }
    else {
      base=0; // this is wrong but was corrected through the call to pastix_checkMatrix
      format=CSC;
      rtn = convert(Mat, format, base, verbose);
      if(rtn!=0) return(-1);
      rtn=factorize_pastix_sequential(slv, Mat, verbose, debug);
      }
    }
  else if  (strcmp(solver_name,"PASTIX-SEQUENTIAL") == 0) {
    base=0; // this is wrong but was corrected through the call to pastix_checkMatrix
    format=CSC;
    rtn = convert(Mat, format, base, verbose);
    if(rtn!=0) return(-1);
    rtn=factorize_pastix_sequential(slv, Mat, verbose, debug);
    }
  else {
    printf("factorize, unknown solver : %s \n",solver_name);
    rtn = -1;
    }
    
#ifdef  HAVE_MPI
   if(MPI_Init_thread_DONE==1) {
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

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorize_mumps(solver_t *slv,triplet *Mat, int verbose, bool debug) {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#ifdef MUMPS
  DMUMPS_STRUC_C *id;
  int myid, ierr, nglob=0,count,rank,size;
  char *solver_name;
    
#ifdef HAVE_MPI
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
  ierr = 0;
#endif
    
  id = (DMUMPS_STRUC_C *) slv->parameters;
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
  (*id).a_loc =  Mat->x;  

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


  dmumps_c(id);
  ierr=(*id).info[0];
  if(ierr==-1)  printf("cpu %d: DMUMPS factorisation failed, status=%d due to cpu %d\n",rank, ierr, id->info[1]);
  if(ierr==-9)  printf("cpu %d: DMUMPS factorisation failed, status=%d due to memory deficience %d\n",rank, ierr, id->info[1]);
  if(ierr==-10) printf("cpu %d: DMUMPS factorisation failed, status=%d due to singular matrix\n",rank, ierr);
  if(ierr==-13) printf("cpu %d: DMUMPS factorisation failed, status=%d due to memory deficience %d\n",rank, ierr, id->info[2]);
  if(ierr==0) {
    if(verbose==1) {
      printf("%s : MEMORY USE, %d used by cpu, %d total Mo, status=%d\n",__FUNCTION__, (*id).infog[20],(*id).infog[21], ierr);
      }
    if(verbose==1 && rank==0) {
      printf("MEMORY SCALING DIAG %d %d\n",size,(*id).infog[21]);
      }
    }
  return(ierr);
  
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorize_umfpack(solver_t *slv, triplet *Mat, int verbose) {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  int status;
  
#ifdef UMFPACK
  UMFPACK_STRUC_C *id;
  char *solver_name;
  double *Info, *Control, *Ax;
  int *Ap, *Ai,n,nnz;
  void *Symbolic, *Numeric ;
  
  solver_name=slv->name;
  if (strcmp(solver_name,"UMFPACK") == 0) {
    if(verbose) printf("Factorize UMFPACK...\n");
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
    status = umfpack_di_symbolic (n, n, Ap, Ai, Ax, &Symbolic, Control, Info) ; 
    if (status != 0) {
      printf("Factorisation status=%d\n",status);
      umfpack_di_report_info (Control, Info) ;
      umfpack_di_report_status (Control, status) ;
        /*error ("umfpack_di_symbolic failed") ;*/
      }
    status = umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, Control, Info) ;
    id->Symbolic=Symbolic;
    id->Numeric=Numeric; 
    if (status != 0) {
      printf("Factorisation status=%d\n",status);
      umfpack_di_report_info (Control, Info) ;
      umfpack_di_report_status (Control, status) ;
        /*error ("umfpack_di_numeric failed") ;*/
      }
    return(status);
  }
#else
  status=-1;
  printf("Please compile poc-solvers library with -DUMFPACK \n");
  return(status);
#endif
  return(status); 
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorize_hips(solver_t *slv,triplet *Mat, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,ierr;
#ifdef HIPS
#ifdef HAVE_MPI
  dHIPS_STRUC_C *id;
  INTS  sym_pattern, sym_matrix;
  INTS id_hips, i, j;
  INTS *unknownlist, *hipsnodelist;
  COEFD *x, *rhsloc;
  INTS proc_id, n, ln;
  INTL *ia, nnz;
  INTS *ja;
  COEFD *a;
  INTS domsize, nprocs;
  INTS pbegin, pend;
  INTS *mapptr,  *mapptr2,*mapp,*mapp2;
  INTS p; 
  MPI_Comm world;
  int nglob, nnzglob, tmp, bcl, idi, info;
  double *x1,*x2, tol;
  bool talk;
  int root=0;
  
  id = (dHIPS_STRUC_C *) slv->parameters; 
  
  if(slv->rank==-1) return(0);
    
  id_hips=id->id_hips;
  world=slv->communicator;
  
  MPI_Comm_rank(world, &proc_id);
  MPI_Comm_size(world, &nprocs);
  
  status = MPI_Barrier(world);
  
//   debug=true;
  talk=(debug || (verbose==1 && proc_id==0));
  
  /* nglob */
  tmp=Mat->nrow; 
    
  MPI_Allreduce( &tmp, &nglob, 1, MPI_INT, MPI_SUM, world );
  
  if (talk) printf("%s cpu=%d (over %d) : id_hips=%d nloc=%d nglob=%d nnz=%d\n", __FUNCTION__, proc_id, nprocs, id_hips, tmp, nglob, Mat->nnz);
  
/*------------------------------------------------------------------------------
  parameter domsize is an argument of testHIPS.exe */
  domsize = nglob/nprocs/2;

  if (talk) printf("%s cpu=%d : domsize=%d\n",__FUNCTION__,proc_id,domsize);

  dHIPS_SetOptionINT(id_hips, dHIPS_PARTITION_TYPE, 0); 
  dHIPS_SetOptionINT(id_hips, dHIPS_DOMSIZE, domsize);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  OPTIONS
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   ierr= dHIPS_SetDefaultOptions (id_hips, dHIPS_HYBRID); 
  ierr = dHIPS_SetDefaultOptions (id_hips, dHIPS_ITERATIVE); 
  tol=1e-12;
  dHIPS_SetOptionREAL(id_hips,dHIPS_DROPTOL0,tol);
  dHIPS_SetOptionREAL(id_hips,dHIPS_DROPTOL1,tol);
  
  if(verbose==1)
    info=4;
  else
    info=0;
  dHIPS_SetOptionINT (id_hips, dHIPS_VERBOSE, info);
  
//    dHIPS_SetOptionREAL(id_hips, INTS number, REAL value);
  sym_matrix=0;
  dHIPS_SetOptionINT (id_hips, dHIPS_SYMMETRIC, sym_matrix);
/*------------------------------------------------------------------------------
  C : numbering starts from 0 but here starts from 1 */
  dHIPS_SetOptionINT(id_hips, dHIPS_FORTRAN_NUMBERING, 0);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  ENTER THE GRAPH
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  ierr = dHIPS_GraphBegin(id_hips,nglob, Mat->nnz);
  dHIPS_ExitOnError(ierr);
  
  for (bcl=0;bcl<Mat->nnz;bcl++) {
    ierr = dHIPS_GraphEdge(id_hips, Mat->i[bcl]-1, Mat->j[bcl]-1);
    }
  ierr = dHIPS_GraphEnd(id_hips);
  dHIPS_ExitOnError(ierr);
  
  status = MPI_Barrier(world);
  
  if (talk) printf("%s cpu=%d : GraphEnd OK...\n",__FUNCTION__,proc_id);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   ENTER A USER PARTITION
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  mapptr = (INTS *)malloc(sizeof(INTS)*(nprocs+1));
  mapptr2= (INTS *)malloc(sizeof(INTS)*(nprocs+1));
  
  mapp =   (INTS *)malloc(sizeof(INTS)*nglob);
  mapp2=   (INTS *)malloc(sizeof(INTS)*nglob);

  for (bcl=0;bcl<nprocs;bcl++)  mapptr2[bcl]=0;
  for (bcl=0;bcl<nglob;bcl++)   mapp2[bcl]  =0;
  for (bcl=0;bcl<nprocs;bcl++)  mapptr[bcl] =0;
  for (bcl=0;bcl<nglob;bcl++)   mapp[bcl]   =0;
  
  for (bcl=proc_id+1;bcl<nprocs+1;bcl++) mapptr2[bcl]=Mat->nrow;

  MPI_Allreduce(mapptr2, mapptr, nprocs+1, MPI_INT, MPI_SUM, world );
 
/*------------------------------------------------------------------------------
  handle empty partition case */
  if(Mat->nnz!=0) {
    idi=0;
    mapp2[idi+mapptr[proc_id]] =  Mat->i[0]-1;
    for (bcl=1; bcl < Mat->nnz; bcl++) {
      if ( mapp2[idi+mapptr[proc_id]] != (Mat->i[bcl]-1) ) {
        idi++;
        mapp2[idi+mapptr[proc_id]] = Mat->i[bcl]-1;
        }
      }
    }
    
/*------------------------------------------------------------------------------
  Only the master processor (0 by default) needs to enter the partition */
//   MPI_Allreduce( mapp2, mapp, nglob, MPI_INT, MPI_SUM, world );
  MPI_Reduce( mapp2, mapp, nglob, MPI_INT, MPI_SUM, root, world );
            
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Set the unknownlist 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  handle empty partition case */
  if(Mat->nnz!=0) {
    unknownlist= (INTS *)malloc(sizeof(INTS)*(Mat->nrow));
    unknownlist[0]=Mat->i[0]-1;
    idi=0;
    for (bcl=1; bcl < Mat->nnz; bcl++) {
      if ( unknownlist[idi] != (Mat->i[bcl]-1) ) {
        idi++;
        unknownlist[idi] = Mat->i[bcl]-1;
        }
      }
    }
  else unknownlist=(INTS *) malloc(sizeof(INTS)*(1));

  if (talk) printf("%s cpu %d : ndofs=%d (over %d)\n",__FUNCTION__,proc_id,idi+1,Mat->nrow);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Enter Partition 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (proc_id ==0) {
    ierr = dHIPS_SetPartition(id_hips, nprocs, mapptr, mapp);
    dHIPS_ExitOnError(ierr);
    }
  
  status = MPI_Barrier(world);
  if (talk) printf("%s cpu %d : dHIPS_SetPartition done\n",__FUNCTION__,proc_id);
  
  free(mapptr2);
  free(mapptr);
  free(mapp2);
  free(mapp);
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  GET THE LOCAL UNKNOWN LIST
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   ierr = dHIPS_GetLocalUnknownNbr(id_hips, &ln);
//   dHIPS_ExitOnError(ierr);
//   hipsnodelist = (INTS *)malloc(sizeof(INTS)*ln);
//   ierr = dHIPS_GetLocalUnknownList(id_hips, hipsnodelist);
//   dHIPS_ExitOnError(ierr);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  ENTER THE MATRIX COEFFICIENT
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /** The processor enter the rows pbegin to pend-1 (naive partition) **/
  /** In the dHIPS_ASSEMBLY_FOOL mode any processor can enter any coefficient :
      nevertheless it is better to enter as much as possible coefficients
      that are in the local matrix of the processor **/

  ierr = dHIPS_AssemblyBegin(id_hips, Mat->nnz, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_FOOL, sym_matrix);
  dHIPS_ExitOnError(ierr);
  
  if (talk) printf("%s cpu %d : AssemblyBegin OK...\n",__FUNCTION__,proc_id);

  for(bcl=0;bcl<Mat->nnz;bcl++) {
    ierr = dHIPS_AssemblySetValue(id_hips, Mat->i[bcl]-1, Mat->j[bcl]-1, Mat->x[bcl]);
    dHIPS_ExitOnError(ierr);
    }

  ierr = dHIPS_AssemblyEnd(id_hips);
  dHIPS_ExitOnError(ierr);
   
  id->unknownlist = unknownlist;
  id->n   = Mat->nrow;
  id->nnz = Mat->nnz;

  if (talk) printf("%s cpu %d : done...\n",__FUNCTION__,proc_id);

  status = MPI_Barrier(world);
  
  ierr=0;
  status=0;
  return(status);
#endif 
#endif
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int factorizez_hips_complex(solver_t *slv,tripletz *Mat, int verbose, bool debug)
//   
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status,ierr;
// #ifdef HIPS
// #ifdef HAVE_MPI
// /*------------------------------------------------------------------------------
//           poc-solvers refurbishing : bug found */
// //   dHIPS_STRUC_C *id;
//   zHIPS_STRUC_C *id;
//   INTS  sym_pattern, sym_matrix;
//   INTS id_hips, i, j;
//   INTS *unknownlist, *hipsnodelist;
//   INTS proc_id, n, ln;
//   INTL *ia, nnz;
//   INTS *ja;
//   INTS domsize, nprocs;
//   INTS pbegin, pend;
//   INTS *mapptr,  *mapptr2,*mapp,*mapp2;
//   INTS p; 
//   MPI_Comm world;
//   int nglob, nnzglob, tmp, bcl, idi, info;
//   double *x1,*x2, tol;
//   bool talk;
//   int root=0;
//   
//   id = (zHIPS_STRUC_C *) slv->parameters; 
//   
//   if(slv->rank==-1) return(0);
//     
//   id_hips=id->id_hips;
//   world=slv->communicator;
//   
//   MPI_Comm_rank(world, &proc_id);
//   MPI_Comm_size(world, &nprocs);
//   
//   status = MPI_Barrier(world);
//   
//   talk=(debug || (verbose==1 && proc_id==0));
//   
// /*------------------------------------------------------------------------------
//   2x ndofs */
//   tmp=Mat->nrow*2; 
//   nnz=4*Mat->nnz;
//   
//   MPI_Allreduce( &tmp, &nglob, 1, MPI_INT, MPI_SUM, world );
//   
//   if (talk) printf("%s cpu=%d (over %d) : id_hips=%d nloc=%d nglob=%d nnz=%d\n", __FUNCTION__, proc_id, nprocs, id_hips, tmp, nglob, nnz);
//   
// /*------------------------------------------------------------------------------
//   parameter domsize is an argument of testHIPS.exe */
//   domsize = nglob/nprocs/2;
// 
//   if (talk) printf("%s cpu=%d : domsize=%d\n",__FUNCTION__,proc_id,domsize);
// 
//   dHIPS_SetOptionINT(id_hips, dHIPS_PARTITION_TYPE, 0); 
//   dHIPS_SetOptionINT(id_hips, dHIPS_DOMSIZE, domsize);
// 
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//   OPTIONS
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// 
//   ierr = dHIPS_SetDefaultOptions (id_hips, dHIPS_ITERATIVE); 
//   tol=1e-18;
//   dHIPS_SetOptionREAL(id_hips,dHIPS_DROPTOL0,tol);
//   dHIPS_SetOptionREAL(id_hips,dHIPS_DROPTOL1,tol);
//   
//   if(verbose==1)
//     info=4;
//   else
//     info=0;
//   dHIPS_SetOptionINT (id_hips, dHIPS_VERBOSE, info);
//   
//   sym_matrix=0;
//   dHIPS_SetOptionINT (id_hips, dHIPS_SYMMETRIC, sym_matrix);
// /*------------------------------------------------------------------------------
//   C : numbering starts from 0 but here starts from 1 */
//   dHIPS_SetOptionINT(id_hips, dHIPS_FORTRAN_NUMBERING, 0);
//   
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//   ENTER THE GRAPH
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
//   
//   ierr = dHIPS_GraphBegin(id_hips,nglob, nnz);
//   dHIPS_ExitOnError(ierr);
//   
//   for (bcl=0;bcl<Mat->nnz;bcl++) {
//     for(i=0;i<2;i++) {     
//       for(j=0;j<2;j++) {
//         ierr = dHIPS_GraphEdge(id_hips, 2*(Mat->i[bcl]-1)+i, 2*(Mat->j[bcl]-1)+j);
//         }
//       }
//     }
//   ierr = dHIPS_GraphEnd(id_hips);
//   dHIPS_ExitOnError(ierr);
//   
//   status = MPI_Barrier(world);
//   
//   if (talk) printf("%s cpu=%d : GraphEnd OK...\n",__FUNCTION__,proc_id);
// 
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//    ENTER A USER PARTITION
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
//   
//   mapptr = (INTS *)malloc(sizeof(INTS)*(nprocs+1));
//   mapptr2= (INTS *)malloc(sizeof(INTS)*(nprocs+1));
//   
//   mapp =   (INTS *)malloc(sizeof(INTS)*nglob);
//   mapp2=   (INTS *)malloc(sizeof(INTS)*nglob);
// 
//   for (bcl=0;bcl<nprocs;bcl++)  mapptr2[bcl]=0;
//   for (bcl=0;bcl<nglob;bcl++)   mapp2[bcl]  =0;
//   for (bcl=0;bcl<nprocs;bcl++)  mapptr[bcl] =0;
//   for (bcl=0;bcl<nglob;bcl++)   mapp[bcl]   =0;
//   
//   for (bcl=proc_id+1;bcl<nprocs+1;bcl++) mapptr2[bcl]=2*Mat->nrow;
// 
//   MPI_Allreduce(mapptr2, mapptr, nprocs+1, MPI_INT, MPI_SUM, world );
//  
// /*------------------------------------------------------------------------------
//   handle empty partition case */
//   if(Mat->nnz!=0) {
//     int running;
//     idi=0;
//     mapp2[idi+mapptr[proc_id]] =  2*(Mat->i[0]-1);
//     running=Mat->i[0]-1;
//     idi++;
//     mapp2[idi+mapptr[proc_id]] =  2*(Mat->i[0]-1)+1;
//     for (bcl=1; bcl < Mat->nnz; bcl++) {
//       if ( running != Mat->i[bcl]-1) {
//         idi++;
//         mapp2[idi+mapptr[proc_id]] = 2*(Mat->i[bcl]-1);
//         running=Mat->i[bcl]-1;
//         idi++;
//         mapp2[idi+mapptr[proc_id]] =  2*(Mat->i[bcl]-1)+1;
//         }
//       }
//     }
//     
// /*------------------------------------------------------------------------------
//   Only the master processor (0 by default) needs to enter the partition */
//   MPI_Reduce( mapp2, mapp, nglob, MPI_INT, MPI_SUM, root, world );
//             
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//   Set the unknownlist 
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// 
// /*------------------------------------------------------------------------------
//   handle empty partition case */
//   if(Mat->nnz!=0) {
//     int running;
//     unknownlist= (INTS *)malloc(sizeof(INTS)*(2*Mat->nrow));
//     unknownlist[0]=2*(Mat->i[0]-1);
//     idi=0;
//     running=Mat->i[0]-1;
//     idi++;
//     for (bcl=1; bcl < Mat->nnz; bcl++) {
//       if ( running != Mat->i[bcl]-1 ) {
//         idi++;
//         unknownlist[idi] = 2*(Mat->i[bcl]-1);
//         running=Mat->i[bcl]-1;
//         idi++;
//         unknownlist[idi] = 2*(Mat->i[bcl]-1)+1;
//         }
//       }
//     }
//   else unknownlist=(INTS *) malloc(sizeof(INTS)*(1));
// 
//   if (talk) printf("%s cpu %d : ndofs=%d (over %d)\n",__FUNCTION__,proc_id,idi+1,2*Mat->nrow);
// 
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//   Enter Partition 
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// 
//   if (proc_id ==0) {
//     ierr = dHIPS_SetPartition(id_hips, nprocs, mapptr, mapp);
//     dHIPS_ExitOnError(ierr);
//     }
//   
//   status = MPI_Barrier(world);
//   if (talk) printf("%s cpu %d : dHIPS_SetPartition done\n",__FUNCTION__,proc_id);
//   
//   free(mapptr2);
//   free(mapptr);
//   free(mapp2);
//   free(mapp);
//      
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//   ENTER THE MATRIX COEFFICIENT
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// 
//   ierr = dHIPS_AssemblyBegin(id_hips, nnz, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_FOOL, sym_matrix);
//   dHIPS_ExitOnError(ierr);
//   
//   if (talk) printf("%s cpu %d : AssemblyBegin OK...\n",__FUNCTION__,proc_id);
// 
//   for(bcl=0;bcl<Mat->nnz;bcl++) {
//     int i=Mat->i[bcl]-1;
//     int j=Mat->j[bcl]-1;
//     double value;
//     value=real(Mat->x[bcl]);
//     ierr = dHIPS_AssemblySetValue(id_hips, 2*i,   2*j,   value);
//     ierr = dHIPS_AssemblySetValue(id_hips, 2*i+1, 2*j+1, value);
//     value=imag(Mat->x[bcl]);
//     ierr = dHIPS_AssemblySetValue(id_hips, 2*i,   2*j+1, -value);
//     ierr = dHIPS_AssemblySetValue(id_hips, 2*i+1, 2*j,    value);
//     dHIPS_ExitOnError(ierr);
//     }
// 
//   ierr = dHIPS_AssemblyEnd(id_hips);
//   dHIPS_ExitOnError(ierr);
//    
//   id->unknownlist = unknownlist;
//   id->n   = 2*Mat->nrow;
//   id->nnz = 4*Mat->nnz;
// 
//   if (talk) printf("%s cpu %d : done...\n",__FUNCTION__,proc_id);
// 
//   status = MPI_Barrier(world);
//   
//   ierr=0;
//   status=0;
//   return(status);
// #endif 
// #endif
// }



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorize_lapack(solver_t *slv, triplet *Mat) {
    
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

#ifdef LAPACK
  int status; 

  slv->Mat=Mat;
  Mat->pivot = (int *) malloc (Mat->nrow * sizeof (int)) ;
/*   printf("Appel a dgbtrf Mat->lda=%d \n",Mat->lda); */
/*   printf("               Mat->ml=%d, Mat->mu=%d \n", */
/*          Mat->ml, Mat->mu); */
/*   for (bcl=0; bcl < Mat->nrow; bcl++) { */
/*     printf(" Mat->pivot[%d]=%d \n",bcl,Mat->pivot[bcl]); */
/*   } */
  dgbtrf_(         &Mat->nrow, &Mat->nrow,& Mat->ml, &Mat->mu, 
                   Mat->a , &Mat->lda, Mat->pivot,&status);
  /*  printf("Apres dgbtrf \n");*/
/*   for (bcl=0; bcl < Mat->nrow; bcl++) { */
/*     printf(" Mat->pivot[%d]=%d \n",bcl,Mat->pivot[bcl]); */
/*   } */
  /*DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )*/
  return(status);
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorize_spdomestic(solver_t *slv, triplet *Mat) {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  int status=0,n,order=3,bcl;
  double tol=1.0; /* unsymetrique usage*/
  cs *A = (cs *) cs_calloc (1, sizeof (cs)) ;    /* allocate the cs struct */  ; 
  css *S ;
  csn *N ;
  
  Mat->A = A;
  A->nzmax = Mat->nnz;
  A->m =  Mat->nrow;
  A->n =  Mat->ncol;
  A->p =  Mat->i;
  A->i =  Mat->j;
  A->x =  Mat->x;
  A->nz=-1;
  n = A->n;
  S = cs_sqr (order, A, 0) ;                    /* ordering and symbolic analysis */
  N = cs_lu (A, S, tol) ;                    /* numeric LU factorization */
  bcl = N->L->nz;
  Mat->S = S;
  Mat->N = N;
  slv->Mat = Mat;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorize_pastix_sequential(solver_t *slv,triplet *Mat, int verbose, bool debug) {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/    
  int status;
#ifdef PASTIX
  PASTIX_STRUC_C *id;
  char *solver_name;

  double *Ax;
  pastix_int_t *Ap, *Ai;
  int  n,nnz,nbthread;
  double *rhs         = NULL;           /* right hand side                           */
  int               verbosemode;        /* Level of verbose mode (0, 1, 2)           */
  char             *type        = NULL; /* type of the matrix                        */
  int     bcl, bcl2;
  int    flagsym=0;
  pastix_int_t    ncol;                 /* Size of the matrix                        */
//   int verbose=0;
  
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

    verbosemode=0;    
    
    d_pastix_checkMatrix(MPI_COMM_WORLD, verbosemode,
                     API_SYM_NO, API_YES,
                     n, &Ap, &Ai, &Ax, NULL,1);

  
    id->iparm[IPARM_MODIFY_PARAMETER] = API_NO;
    d_pastix(&(id->pastix_data), MPI_COMM_WORLD,
         n, Ap, Ai, Ax,
         id->perm, id->invp, rhs, 1, id->iparm, id->dparm);

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
  id->iparm[IPARM_THREAD_NBR] = nbthread;
  id->iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
  id->iparm[IPARM_LEVEL_OF_FILL] = 0;
  id->iparm[IPARM_START_TASK] = API_TASK_ORDERING;
  id->iparm[IPARM_END_TASK]   = API_TASK_NUMFACT;
//   id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-8;
  id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-12;
  id->iparm[IPARM_VERBOSE]             = verbosemode;

  /*iparm[IPARM_END_TASK]   = API_TASK_REFINE;*/
//   id->perm = malloc(n*sizeof(pastix_int_t));
//   id->invp = malloc(n*sizeof(pastix_int_t));
  id->perm = new pastix_int_t[n];
  id->invp = new pastix_int_t[n];
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Factorisation
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  d_pastix(&(id->pastix_data), MPI_COMM_WORLD,
         n, Ap, Ai, Ax,
         id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
  
  status=id->iparm[IPARM_ERROR_NUMBER];

  if (status < 0)  { printf("Factorisation status=%d\n",status);  }
 
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

  int factorize_pastix_parallel_distributed(solver_t *slv, triplet *Mat, int verbose, bool debug) {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  int k,status;
#ifdef PASTIX
#ifdef HAVE_MPI
  PASTIX_STRUC_C *id;
  char *solver_name;
  double *Ax;
  pastix_int_t *Ap, *Ai;
  int  n,nnz,nbthread;
  double *rhs = NULL; /* right hand side                      */
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
  
  d_pastix_checkMatrix(MPI_COMM_WORLD, verbosemode, API_SYM_NO, API_YES,
                      n, &Ap, &Ai, &Ax, &loc2glb, 1);
  
/*------------------------------------------------------------------------------
  get (possibly changed) new addresses */
  Mat->i=Ap;
  Mat->j=Ai;
  Mat->x=Ax;
  
  id->iparm[IPARM_MODIFY_PARAMETER] = API_NO;
  
  d_dpastix(&(id->pastix_data), MPI_COMM_WORLD, n, Ap, Ai, Ax, loc2glb,
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
//   id->perm = malloc(n*sizeof(pastix_int_t));
  id->perm = new pastix_int_t[n];

//   id->invp = malloc(ncolglob*sizeof(pastix_int_t));
  
  MPI_Barrier(MPI_COMM_WORLD);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Factorisation
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(verbose==1) printf("%s pastix factorisation start\n",__FUNCTION__);
  
  d_dpastix(&(id->pastix_data), MPI_COMM_WORLD, n, Ap, Ai, Ax, loc2glb,
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

// #endif // moved at end

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int factorize_pastix_parallel(solver_t *slv, triplet *Mat, int verbose, bool debug) {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  int status;
#ifdef PASTIX
#ifdef HAVE_MPI
  PASTIX_STRUC_C *id;
  char *solver_name;
  double *Ax;
  pastix_int_t *Ap, *Ai;
  int  n,nnz,nbthread;
  double *rhs         = NULL; /* right hand side                     */
  int    verbosemode;         /* Level of verbose mode (0, 1, 2)     */
  char   *type        = NULL; /* type of the matrix                  */
  int    bcl, bcl2;
  int    flagsym=0;
  pastix_int_t    ncol;       /* Size of the matrix                  */
  // Global variables
  int ncolglob, nnzglob;
  int *Apglob, *Aiglob, *colptr, *rows;
  double *Axglob, *values;
  int             mpid,sz;
  MPI_Comm world;
  
  if(slv->rank==-1) return(0);

  solver_name=slv->name;
    
  id = (PASTIX_STRUC_C *) slv->parameters;
  slv->Mat = Mat; 
  Ap=Mat->i;
  Ai=Mat->j;
  Ax=Mat->x;
  n=Mat->nrow;
  nnz=Mat->nnz;

  world=slv->communicator;
//   world=MPI_COMM_WORLD;

//   MPI_Comm_rank(MPI_COMM_WORLD, &mpid);
//   MPI_Comm_size(MPI_COMM_WORLD, &sz);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Dispatch Matrix global matrix
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  d_coo2glob(world, n,nnz,Ap,Ai,Ax, &ncolglob, &nnzglob, &Apglob, &Aiglob, &Axglob);
  
  free(Ap);
  free(Ai);
  free(Ax);
   
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  convert structure to CSC
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  d_globcoo2csc(world, ncolglob, nnzglob, Apglob , Aiglob, Axglob , &colptr, &rows, &values);
  
  free(Apglob);
  free(Aiglob);
  free(Axglob);
   
  Ap = colptr;
  Ai = rows;
  Ax = values;
   
  Mat->nrow = ncolglob;
  Mat->nnz  = nnzglob;
   
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Initialize parameters to default values
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  verbosemode=verbose;  
    
  d_pastix_checkMatrix(world, verbosemode, API_SYM_NO, API_YES, ncolglob, &Ap, &Ai, &Ax, NULL,1);

  Mat->i=Ap;
  Mat->j=Ai;
  Mat->x=Ax;
  Mat->comptype=CSC;
    
  id->iparm[IPARM_MODIFY_PARAMETER] = API_NO;              /* ?                      */
  
  d_pastix(&(id->pastix_data), world,  ncolglob, Ap, Ai, Ax, id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
    
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
  id->iparm[IPARM_SYM]                 = API_SYM_NO;         /* non-symmetric matrix   */
  id->iparm[IPARM_FACTORIZATION]       = API_FACT_LU;        /* LU factorization       */
  id->iparm[IPARM_THREAD_NBR]          = nbthread;           /* number of threads      */
  id->iparm[IPARM_MATRIX_VERIFICATION] = API_NO;             /* no matrix verification */
  id->iparm[IPARM_LEVEL_OF_FILL]       = 0;                  /* ?                      */
  id->iparm[IPARM_START_TASK]          = API_TASK_ORDERING;  /* first task             */
  id->iparm[IPARM_END_TASK]            = API_TASK_NUMFACT;   /* last task              */
//   id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-8;
  id->dparm[DPARM_EPSILON_REFINEMENT]  = 1e-12;              /* ?                      */
  id->iparm[IPARM_VERBOSE]             = verbosemode;

  /*iparm[IPARM_END_TASK]   = API_TASK_REFINE;*/

//   id->perm = malloc(ncolglob*sizeof(pastix_int_t));
//   id->invp = malloc(ncolglob*sizeof(pastix_int_t));
  id->perm = new pastix_int_t[ncolglob];
  id->invp = new pastix_int_t[ncolglob];
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Factorisation
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(verbose==1) printf("%s pastix factorisation start\n",__FUNCTION__);
  
  d_pastix(&(id->pastix_data), world, ncolglob, Ap, Ai, Ax, id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
  
  if(verbose==1) printf("%s pastix factorisation end\n",__FUNCTION__);
  
  status=id->iparm[IPARM_ERROR_NUMBER];

  if (status < 0)  {
    printf("d_pastix factorization status=%d\n",status);  
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

#endif 

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  depecrated 1.0 sources with no replcament in 2.0 code
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorize_hypre(solver_t *slv,triplet *Mat) {
    
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  int status;
#ifdef HYPRE
#ifdef HAVE_MPI
  HYPRE_STRUC_C *id_hypre; 
  HYPRE_IJMatrix A;
  HYPRE_ParCSRMatrix parcsr_A;
  HYPRE_IJVector b;
  HYPRE_ParVector par_b;
  HYPRE_IJVector x;
  HYPRE_ParVector par_x;
  HYPRE_Solver hy_solver, hy_precond;
  HYPRE_STRUC_C *id;
  int ilower,iupper,bcl1,bcl2,ptr;
  char *solver_name;
  char filename[50];
  int myid, ierr, ncolmax, ncol, oldlig;
  FILE *out;
  int lig, col;
  double val;
 

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  sprintf(filename,"solver_fatorize_hyre.1.%d.txt",myid);
  solver_name=slv->name;
  id = slv->parameters;
  printf("%d, n=%d,nnz=%d \n",myid,Mat->nrow,Mat->nnz);
  printf("%d, ligdeb=%d, ligfin=%d \n",myid,Mat->i[0],Mat->i[Mat->nrow-1]);
  out=fopen(filename,"w");
  ilower=Mat->nrow;
  iupper=0;
  ncolmax=0;
  ncol=1;
  oldlig=Mat->i[0];
  bcl1=0;
  fprintf(out, "%d %d %lf \n",Mat->i[bcl1],Mat->j[bcl1],Mat->x[bcl1]); 
  for(bcl1 = 1; bcl1 < Mat->nnz; bcl1++) {
    fprintf(out, "%d %d %lf \n",Mat->i[bcl1],Mat->j[bcl1],Mat->x[bcl1]); 
    if (ilower > Mat->i[bcl1]) ilower=Mat->i[bcl1];
    if (iupper < Mat->i[bcl1]) iupper=Mat->i[bcl1];
    if (Mat->i[bcl1] == oldlig) {
      ncol=ncol+1;
    }
    else {
      oldlig=Mat->i[bcl1];
      if (ncol > ncolmax) ncolmax=ncol;
      ncol=1;
    }
  }
  fclose(out);
  printf("%d, ilower=%d, iupper=%d ncolmax=%d \n",myid,ilower,iupper,ncolmax);
  /* Create the matrix.
      Note that this is a square matrix, so we indicate the row partition
      size twice (since number of rows = number of cols) */
   HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
   status = MPI_Barrier(MPI_COMM_WORLD);
   printf(" IJMatrixCreate OK.....\n");
   /* Choose a parallel csr format storage (see the User's Manual) */
   HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
   /* Initialize before setting coefficients */
   HYPRE_IJMatrixInitialize(A);
   status = MPI_Barrier(MPI_COMM_WORLD);
   printf(" HYPRE_IJMatrixInitialize OK.....\n");

  sprintf(filename,"solver_fatorize_hyre.2.%d.txt",myid);
  out=fopen(filename,"w");   
  for(bcl1 = 0; bcl1 < Mat->nnz; bcl1++) {
    lig =  Mat->i[bcl1]-1;
    col =  Mat->j[bcl1]-1;
    val =  Mat->x[bcl1];
    ncol=1;
    fprintf(out, "%d %d %lf \n",lig,col,val); 
    /*HYPRE_IJMatrixSetValues(A, 1, &ncol, &lig, &col, &val);*/

  }
  fclose(out);
   status = MPI_Barrier(MPI_COMM_WORLD);
   exit(-23);

   /* Assemble after setting the coefficients */
   HYPRE_IJMatrixAssemble(A);
   /* Get the parcsr matrix object to use */
   HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);
   id->A=A;
   id->parcsr_A=parcsr_A;
   /*---pour test */
   HYPRE_IJMatrixPrint(id->A, "IJ.out.A");
   status = MPI_Barrier(MPI_COMM_WORLD);
   exit(-23);
#endif
#endif 
  return(status); 
}
