
/*******************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2012

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France

*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <time.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "fe.h"
#include "matrix.h"

extern int matrix_diagonal(hyperzmatrix_t *A);

int nbsyslinearz=0;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LinearSystem_factorisation_sequential(hyperzmatrix_t *M, int solver_id)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes

  Check :

  Notes

    stype=0: unsymetric,
    stype>0: matrix is square and symmetric. Entries in the lower triangular
    stype<0: matrix is square and symmetric. Entries in the upper triangular

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  int i, lda, status, ierr;
  int *colstr, *rowind;
  char mtxtyp[2], cpivot, ordmthd;
  int outunt, msglvl = 0, error, nrhs, ldrhs, neqns, nz;
  float *rhs;
  double handle[150];
  int hbw, neq;
  int verbose;

  char c;
  int *Tj;
  FILE *out;
  clock_t d1, d2;
  double duree1, duree2;
  int count, m, k, n;
  char **args;
  tripletz *M1;
  solver_t *solveur1;
  
  // init triplet
  M->t.comptype=CSC;
  M->t.base=0;

  M->t.stype=0;

  if(M->neq()==0) {
    check_error(-1, "neq field is zero", __LINE__, __FILE__, 1);
    }

  M->t.nnz=M->ordering->pointer[M->neq()];
  M->t.nrow=M->neq();
  M->t.ncol=M->neq();

  M->t.i = M->ordering->pointer;
  M->t.j = M->ordering->incidence;
  M->t.x = M->packed;

  M1 = &(M->t);

  switch (solver_id) {
/*------------------------------------------------------------------------------
    DOMESTIC*/
    case SOLVER_ID_DOMESTIC:
      printf("built-in solver\n");
//       status = band_unpack02(M->packed, M->ordering->incidence, M->ordering->cardinal, M->neq, M->hbw,&(M->band));
//       status = d_LUDV(M->band, M->neq, M->hbw);
      check_error(-1, "complex solver not implemented", __LINE__, __FILE__, 1);
      break;

/*------------------------------------------------------------------------------
    LAPACKF_*/
    case SOLVER_ID_LAPACK:
      printf("lapack band matrix solver\n");
#   ifdef LAPACKF_
      hbw = M->hbw;
      if(hbw==-1) {
        M->ordering->gethbw();
        hbw   = M->ordering->hbw;
        }
      if(hbw==-1) {
        check_error(-1, "hbw (half-bandwidth) field is zero, check matrix construction routine...", __LINE__, __FILE__, 1);
        }
      lda = 3 * hbw + 1;
      M->pivot = new int[M->neq()];
      status = band_unpack01(M->packed, M->ordering->incidence, M->ordering->cardinal, M->neq, hbw,&(M->band));
      zgbtrf_(&M->neq, &M->neq, &hbw, &hbw, M->band, &lda, M->pivot, &status);
      if(status!=0) {
        check_error(-1, "zgbtrf (double complex band matrix) factorization failed", __LINE__, __FILE__, 1);
        }
#   else
      check_error(-1, "please compile with -DLAPACKF_", __LINE__, __FILE__, 1);
#   endif
      break;

/*------------------------------------------------------------------------------
    SUNPERF*/
    case SOLVER_ID_SUNPERF:
      printf("sunperf packed matrix solver\n");
      check_error(-1, "complex solver not implemented", __LINE__, __FILE__, 1);

/**----------------------------------------------------------------------
    UMFPACK*/
    case SOLVER_ID_UMFPACK:
      ierr=0;
      solveur1=(solver_t *) init_solver_obsolete("UMFPACK", M1->stype);
      ierr = factorizez(solveur1,M1,verbose);
      if (ierr != 0) {
        cout << "UMFPACK factorize error:" << ierr << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
      M->s = solveur1;
      status=ierr;
      break;

/*------------------------------------------------------------------------------
    MUMPS ASSYMETRIC*/
    case SOLVER_ID_MUMPS:

#ifdef MUMPS
      check_error(-1, "complex solver not implemented", __LINE__, __FILE__, 1);
#   else
      check_error(-1, "please compile with -DMUMPS", __LINE__, __FILE__, 1);
#   endif
      break;

/*------------------------------------------------------------------------------
    MUMPS SYMETRIC*/
    case SOLVER_ID_MUMPS_SYM:
#ifdef MUMPS
      check_error(-1, "complex solver not implemented", __LINE__, __FILE__, 1);
#   else
      check_error(-1, "please compile with -DMUMPS", __LINE__, __FILE__, 1);
#   endif
      break;

/**----------------------------------------------------------------------
    SpDOMESTIC*/
    case SOLVER_ID_SpDOMESTIC:
      ierr=0;
      solveur1=(solver_t *) init_solver_obsolete("SpDOMESTIC", M1->stype);
      ierr = factorizez(solveur1,M1,verbose);
      if (ierr != 0) {
        cout << "SpDOMESTIC factorize error:" << ierr << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
      M->s = solveur1;
      status=ierr;
      break;

/**----------------------------------------------------------------------
    HIPS*/
    case SOLVER_ID_HIPS:
      ierr=0;
      solveur1=(solver_t *) init_solver_obsolete("HIPS", M1->stype);
      ierr = factorizez(solveur1,M1,verbose);
      if (ierr != 0) {
        cout << "HIPS factorize error:" << ierr << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
      M->s = solveur1;
      status=ierr;
      break;

/**----------------------------------------------------------------------
    PASTIX*/
    case SOLVER_ID_PASTIX:
      ierr=0;
      solveur1=(solver_t *) init_solver_obsolete("PASTIX", M1->stype);
      ierr = factorizez(solveur1,M1,verbose);
      if (ierr != 0) {
        cout << "PASTIX factorize error:" << ierr << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
      M->s = solveur1;
      status=ierr;
      break;

    default:
      check_error(-1, "solver not implemented", __LINE__, __FILE__, 1);
      break;
    }

/**----------------------------------------------------------------------------
  store solver id */

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes

  Check : MANDATORY

  Notes

    If the following lines trigger a compiling fault, you must either:
    
      1) upgrade poc-solvers library to version 2 and successors and
         update T-UGOm Makefile consequently
      
      2) keep your usual poc-solvers version and add the "id" field
         in solver structure define in solverlib.h

         typedef struct solver_struct
         {
           char *name;
           int id;               <---------------------- here
           void *parameters;
           void *Mat;
         } solver;

         Re-compile library after change was made, the solverlib.h will
         be copied into the include directory.
         You might need to comment "#include mip.h" line in the
         include/solverlib.h as it is needed for library test program but
         redundant (and incomplete) for T-UGOm
         
    1st option is recommended, option 2 to be used as the last bullet
    
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  if(M->s==0) M->s=new solver_t;
  M->s->id=solver_id;

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LinearSystem_solve_sequential(hyperzmatrix_t & M, complex <double> *rhs)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, i, ndiff, ierr, solver_id;
  int hbw,neqs;

  int job, nrhs = 1,bandwidth;
  char cjob = 'N';
  double *x, *tmp, maxdiff, maxdiff2, rnorm;

  void *Numeric;
  clock_t d1, d2;
  double duree1, duree2;
  tripletz *M1;
  solver_t *solveur1;

  int *lign,*colo;
  double *val;
  char filerhs[30], filemat[30];
  int transposed;
  char filename[50];
  FILE *out;

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes

  Check : MANDATORY

  Notes

    If the following lines trigger a compiling fault, you must either:
    
      1) upgrade poc-solvers library to version 2 and successors and
         update T-UGOm Makefile consequently
      
      2) keep your usual poc-solvers version and add the "id" field
         in solver structure define in solverlib.h

         typedef struct solver_struct
         {
           char *name;
           int id;               <---------------------- here
           void *parameters;
           void *Mat;
         } solver;

         Re-compile library after change was made, the solverlib.h will
         be copied into the include directory.
         You might need to comment "#include mip.h" line in the
         include/solverlib.h as it is needed for library test program but
         redundant (and incomplete) for T-UGOm
    
    1st option is recommended, option 2 to be used as the last bullet
    
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  if(M.s==0)
    check_error(-1, "solver structure not properly initialized", __LINE__, __FILE__, 1);
  solver_id=M.s->id;
  
  M1 =  &(M.t);
  solveur1=M.s;

  switch (solver_id) {
    case SOLVER_ID_DOMESTIC:
/**----------------------------------------------------------------------------
      DOMESTIC */
      check_error(-1, "solver not implemented", __LINE__, __FILE__, 1);
      break;

    case SOLVER_ID_LAPACK:
/**----------------------------------------------------------------------------
      LAPACK */
      hbw   = M.hbw;
      hbw = M.hbw;
      if(hbw==-1) {
        M.ordering->gethbw();
        hbw   = M.ordering->hbw;
        }
      if(hbw==-1) {
        check_error(-1, "hbw (half-bandwidth) field is zero, check matrix construction routine...", __LINE__, __FILE__, 1);
        }
      neqs   = M.t.nrow;
      bandwidth = 3 * hbw + 1;
#ifdef LAPACKF_
      zgbtrs_(&cjob, &neqs, &hbw, &hbw, &nrhs, M.band, &bandwidth, M.pivot, rhs, &neqs, &status);
      if(status!=0) {
        check_error(-1, "zgbtrs (double complex band matrix) solver failed", __LINE__, __FILE__, 1);
        }
#endif
      break;

    case SOLVER_ID_SUNPERF:
/**----------------------------------------------------------------------------
      SUNPERF */
      check_error(-1, "solver not implemented", __LINE__, __FILE__, 1);
      break;

    case SOLVER_ID_UMFPACK:
/**----------------------------------------------------------------------------
      UMFPACK */
      d1 = clock();
      switch(M.ordering->type) {
        case ORDER_CSR:
          transposed=1;
          check_error(-1, "unstable mode?", __LINE__, __FILE__, 0);
          break;
        case ORDER_CSC:
          transposed=0;
          break;
        default:
          check_error(-1, "invalid ordering type", __LINE__, __FILE__, 1);
        }
//      transposed=0;
      ierr = solvez(solveur1,rhs,transposed,0);
      d2 = clock();
      duree1 = d2 - d1;
      break;

    case SOLVER_ID_MUMPS:
/**----------------------------------------------------------------------------
      MUMPS */
      check_error(-1, "solver not implemented", __LINE__, __FILE__, 1);
      break;

    case SOLVER_ID_MUMPS_SYM:
/**----------------------------------------------------------------------------
      MUMPS-SYM */
      check_error(-1, "solver not implemented", __LINE__, __FILE__, 1);
      break;

    case SOLVER_ID_SpDOMESTIC:
/**----------------------------------------------------------------------------
      SpDOMESTIC */
      d1 = clock();
      switch(M.ordering->type) {
        case ORDER_CSR:
          transposed=1;
          check_error(-1, "unstable mode?", __LINE__, __FILE__, 0);
          break;
        case ORDER_CSC:
          transposed=0;
          break;
        default:
          check_error(-1, "invalid ordering type", __LINE__, __FILE__, 1);
        }
      ierr = solvez(solveur1,rhs,transposed,0);
      d2 = clock();
      duree1 = d2 - d1;
//      printf("SpDOMESTIC Solve: duree=%e\n",duree1/CLOCKS_PER_SEC);
      printf("SpDOMESTIC solving, elapsed time=%.1e(s)\n",duree1/CLOCKS_PER_SEC);
      break;


    case SOLVER_ID_PASTIX:
/**----------------------------------------------------------------------------
      PASTIX */
      d1 = clock();
      switch(M.ordering->type) {
        case ORDER_CSR:
          transposed=1;
          printf("No transpose resoltution for PASTIX yet \n");
          check_error(-1, "unstable mode?", __LINE__, __FILE__, 0);
          break;
        case ORDER_CSC:
          transposed=0;
          break;
        default:
          check_error(-1, "invalid ordering type", __LINE__, __FILE__, 1);
        }
      transposed=0;
      ierr = solvez(solveur1,rhs,transposed,0);
      d2 = clock();
      duree1 = d2 - d1;
//      printf("PASTIX Solvez: duree=%e\n",duree1/CLOCKS_PER_SEC);
      printf("PASTIX solving, elapsed time=%.1e(s)\n",duree1/CLOCKS_PER_SEC);
      break;


    default:
      check_error(-1, "solver not implemented", __LINE__, __FILE__, 1);
    }

//   nbsyslinearz++;
//   sprintf(filename,"RHSz.%d.%d.txt",solver_id,nbsyslinearz);
//   out=fopen(filename,"w");
//   for (i=0; i <neqs; i++)   fprintf(out,"      %d -> %e+i%e \n",i,real(rhs[i]),imag(rhs[i]));
//   fclose(out);

  return(ierr);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int LinearSystem_solve_parallel(hyperzmatrix_t & M, complex<double> *B)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, i, m, n, ierr, solver_id;
  int neqs;
  double d1, d2;
  double duree1, duree2;
  tripletz *M1;
  solver_t *solveur1;
  int transposed;
  complex<double> *Bglob,*Bglobtmp;
  int nglob,rank, nndes, nloc, root,bcl;
  int allocated=0;
  int gmin,gmax,gsize,proc_id;
  int HIPSRESONB;
  char filename[50];
  FILE *out;
  
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes

  Check : MANDATORY

  Notes

    If the following lines trigger a compiling fault, you must either:
    
      1) upgrade poc-solvers library to version 2 and successors and
         update T-UGOm Makefile consequently
      
      2) keep your usual poc-solvers version and add the "id" field
         in solver structure define in solverlib.h

         typedef struct solver_struct
           {
           char *name;
           int id;               <---------------------- here
           void *parameters;
           void *Mat;
           } solver;

         Re-compile library after change was made, the solverlib.h will
         be copied into the include directory.
         You might need to comment "#include mip.h" line in the
         include/solverlib.h as it is needed for library test program but
         redundant (and incomplete) for T-UGOm
    
    1st option is recommended, option 2 to be used as the last bullet
    
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
#ifdef HAVE_MPI
   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
#endif

  if(M.s==0)
    check_error(-1, "solver structure not properly initialized", __LINE__, __FILE__, 1);
  solver_id=M.s->id;

  M1 =  &(M.t);
  solveur1=M.s;
  nndes = M.distributor->ndof;
  nloc  = M.distributor->nsolved;
  nglob = M.distributor->ngdof;

  gmin  = M.distributor->gmin;
  gmax  = M.distributor->gmax;
  gsize = M.distributor->gsize;

  switch (solver_id) {
    
    case SOLVER_ID_MUMPS:
/**----------------------------------------------------------------------------
      MUMPS */
#ifdef HAVE_MPI
/**----------------------------------------------------------------------------
      Centralised RHS on host 0 */
  
//   status = MPI_Barrier(MPI_COMM_WORLD);
//   printf("cpu %d: nglob=%d nloc=%d nndes=%d\n",gCPU_ID,nglob,nloc,nndes);
//   status = MPI_Barrier(MPI_COMM_WORLD);

/**----------------------------------------------------------------------------
      initialise exchange vector from local values (solved nodes) */
      Bglob    = M.mpi_global;
      Bglobtmp = M.mpi_tmp;
      for (i=0; i <nglob; i++) {
        Bglob[i]   =0.0;
        }
      for (i=0; i <M.d->nsolved; i++) {
        n = M.d->solved[i];
        m = M.d->gsolved[i];
        Bglobtmp[m] = B[n];
        }
      for (i=0; i <M.d->nsolved; i++) {
        m = M.d->gsolved[i];
        }

/**----------------------------------------------------------------------------
      global vector construction: ok if domain are disjoint (summation) */
      status = MPI_Allreduce ( Bglobtmp, Bglob, nglob, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
//      status = MPI_Allreduce ( &(Bglobtmp[gmin]), &(Bglob[gmin]), gsize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#else
      allocated=0;
      Bglob = B;
#endif
      
#ifdef HAVE_MPI
      d1 = MPI_Wtime();
#endif
/**----------------------------------------------------------------------------
      call solver_t */
      ierr = solvez(solveur1,Bglob,0, 0);
#ifdef HAVE_MPI
      d2 = MPI_Wtime();
#endif
      duree1 = d2 - d1;
      printf("MUMPS parallel solver, ellapsed=%e\n",duree1);
#ifdef HAVE_MPI
/**----------------------------------------------------------------------------
      broadcast distributed solution*/
      root=0;
      status=  MPI_Bcast( Bglob, nglob, MPI_DOUBLE, root, MPI_COMM_WORLD);
/**----------------------------------------------------------------------------
      update local solution vector (all noddes)*/
      for (n=0; n <nndes; n++) {
        m = M.d->gindex[n];
        B[n] = Bglob[m];
        }
#endif
      break;

    case SOLVER_ID_HIPS:
/**----------------------------------------------------------------------------
      HIPS */
#ifdef HAVE_MPI
      Bglob    = M.mpi_global;
      Bglobtmp = M.mpi_tmp;
      for (i=0; i <nglob; i++) Bglob[i]   =0.0;
   
      for (i=0; i <M.d->nsolved; i++) {
        n = M.d->solved[i];
        m = M.d->gsolved[i];
        Bglobtmp[m] = B[n];
        }
/**----------------------------------------------------------------------------
      global vector construction: ok if domain are disjoint (summation) */
      status = MPI_Allreduce ( Bglobtmp, Bglob, nglob, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

      ierr = solvez(solveur1,Bglob,0, 0);
  
/**----------------------------------------------------------------------------
      update local solution vector (all noddes)*/
      for (n=0; n <nndes; n++) {
        m = M.d->gindex[n];
        B[n] = Bglob[m];
        }
#endif
      break;

/**----------------------------------------------------------------------------
    OTHERS */
    default:
      printf("matrix %s\n",M.name.c_str());
      check_error(-1, "solver not implemented for parallel computing", __LINE__, __FILE__, 1);
      break;
    }
    
#ifdef HAVE_MPI
  status = MPI_Barrier(MPI_COMM_WORLD);
#endif
 
  return(ierr);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LinearSystem_solve(hyperzmatrix_t & M, complex<double> *B)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  extern int print_vector(char* file, mesh_t mesh, int discretisation, double *val);
  
  if(M.distributor==0) {
#ifdef HAVE_MPI
    status = MPI_Barrier(MPI_COMM_WORLD);
#endif
    check_error(-1, "WARNING !!! node distribution structure not initialized, to be fixed", __LINE__, __FILE__, 0);
    M.distributor=new distribution_t;
    M.distributor->context=SEQUENTIAL_COMPUTING;
#ifdef HAVE_MPI
    status = MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
    
/**----------------------------------------------------------------------------
  special treatment for diagonal matrices */
  if(M.IsDiagonal) {
/**----------------------------------------------------------------------------
    WARNING : 2013-06-25, inversion cancelled to avoid matrix modification*/
    for(size_t n=0; n< M.neq(); n++) {
//      B[n]*=M.packed[n];
      B[n]/=M.packed[n];
      }
    return(0);
    }
  
  switch(M.distributor->context) {
    case SEQUENTIAL_COMPUTING:
//      printf("call LinearSystem_solve_sequential, context=%d\n",M.d->context);
      status= LinearSystem_solve_sequential(M, B);
      break;
    case PARALLEL_COMPUTING:
//      printf("call LinearSystem_solve_parallel_new, context=%d\n",M.d->context);
      status= LinearSystem_solve_parallel(M, B);
      break;
    default:
      printf("LinearSystem_solve, context=%d\n",M.distributor->context);
      check_error(-1, "invalid solver context", __LINE__, __FILE__, 1);
      break;
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int LinearSystem_factorisation_parallel(hyperzmatrix_t *M, int solver_id)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    stype=0: unsymetric
    stype>0: matrix is square and symmetric. Entries in the lower triangular
    stype<0: matrix is square and symmetric. En

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status, ierr;
  int verbose=0;
  int error;
  int neq;
  tripletz *triplet;
  solver_t *solveur1;
  
  if(M->distributor==0) {
    check_error(-1, "node distribution structure not initialized", __LINE__, __FILE__, 1);
    }

  triplet = &(M->t);
  
  printf("(complex) LinearSystem_factorisation_parallel, solver id=%d\n",solver_id);
  
  switch (solver_id) {
/**----------------------------------------------------------------------------
    DOMESTIC*/
    case SOLVER_ID_DOMESTIC:
      check_error(-1, "not parallelized", __LINE__, __FILE__, 1);
      break;
  
/**----------------------------------------------------------------------------
    MUMPS ASSYMETRIC*/
    case SOLVER_ID_MUMPS:
      ierr=0;
      solveur1=(solver_t *) init_solver_obsolete("MUMPS", triplet->stype);
      ierr = factorizez(solveur1,triplet,verbose);
      if (ierr != 0) {
        cout << "MUMPS factorizez error:" << ierr << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
//      printf(" Symetrie = %d \n",coo_issym(triplet));
      M->s  = solveur1;
      status=ierr;
      break;
  
/**----------------------------------------------------------------------------
    MUMPS SYMETRIC*/
    case SOLVER_ID_MUMPS_SYM:
      ierr=0;
      triplet->stype=1; /*assume symetrie*/
      solveur1=(solver_t *) init_solver_obsolete("MUMPS", triplet->stype);
      ierr = factorizez(solveur1,triplet,verbose);
      if (ierr != 0) {
        cout << "MUMPS factorizez error:" << ierr << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
//      printf(" Symetrie = %d \n",coo_issym(triplet));
      M->s = solveur1;
      status=ierr;
      break;
  
/**----------------------------------------------------------------------------
    SpDOMESTIC*/
    case SOLVER_ID_SpDOMESTIC:
      check_error(-1, "not parallelized", __LINE__, __FILE__, 1);
      break;
      
/**----------------------------------------------------------------------------
    HYPRE */
    case SOLVER_ID_HYPRE:
#ifdef HYPRE
      ierr=0;
      solveur1=(solver_t *) init_solver_obsolete("HYPRE", triplet->stype);
      ierr = factorizez(solveur1,triplet,verbose);
      if (ierr != 0) {
        cout << "HYPRE factorizez error:" << ierr << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
      printf(" Symetrie = %d \n",coo_issym(triplet));
      M->s  = solveur1;
#   else
      check_error(-1, "please compile with -DHYPRE", __LINE__, __FILE__, 1);
#   endif
      status=ierr;
      break;
      
/**----------------------------------------------------------------------------
    HIPS */
    case SOLVER_ID_HIPS:
      ierr=0;
      solveur1=(solver_t *) init_solver_obsolete("HIPS", triplet->stype);
      printf("(complex) LinearSystem_factorisation_parallel, initialisation done\n");
      ierr = factorizez(solveur1,triplet,verbose);
      if (ierr != 0) {
        cout << "HIPS factorizez error:" << ierr << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
      M->s  = solveur1;
      status=ierr;
      break;
   
/**----------------------------------------------------------------------------
    PASTIX */
    case SOLVER_ID_PASTIX:
      ierr=0;
      solveur1=(solver_t *) init_solver_obsolete("PASTIX", triplet->stype);
      ierr = factorizez(solveur1,triplet,verbose);
      if (ierr != 0) {
        cout << "PASTIX factorizez error:" << ierr << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
      M->s  = solveur1;
      status=ierr;
      break;
   
/**----------------------------------------------------------------------------
    */
    default:
      check_error(-1, "solver not implemented", __LINE__, __FILE__, 1);
    }
   
/**----------------------------------------------------------------------------
  store solver id */
  if(M->s==0) M->s=new solver_t;
  M->s->id=solver_id;
  printf(" LinearSystem_factorisation (parallel) status: %d\n",ierr);

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int LinearSystem_factorisation(hyperzmatrix_t *M, int solver_id)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, lda, status, ierr;
  if(M->distributor==0) {
#ifdef HAVE_MPI
    status = MPI_Barrier (  MPI_COMM_WORLD);
#endif
    check_error(-1, "node distribution structure not initialized", __LINE__, __FILE__, 1);
#ifdef HAVE_MPI
    status = MPI_Barrier (  MPI_COMM_WORLD);
#endif
    }
/**----------------------------------------------------------------------------
  special treatment for diagonal matrices */
  status=matrix_IsDiagonal(*M);
  if(M->IsDiagonal) {
    for(size_t n=0; n< M->neq(); n++) {
/**----------------------------------------------------------------------------
      WARNING : 2013-06-25, inversion cancelled to avoid matrix modification*/
      if(M->packed[n]!=0.) {
//        M->packed[n]=1./M->packed[n];
        }
      else {
        return(-1);
        }
      }
    return(0);
    }
    
    
  switch(M->distributor->context) {
    case SEQUENTIAL_COMPUTING:
      status=LinearSystem_factorisation_sequential(M,  solver_id);
      break;
    case PARALLEL_COMPUTING:
      status=LinearSystem_factorisation_parallel(M,  solver_id);
      M->mpi_global = new complex<double>[M->distributor->ngdof];
      M->mpi_tmp    = new complex<double>[M->distributor->ngdof];
      for(i=0;i<M->distributor->ngdof;i++) M->mpi_tmp[i]=0;
      break;
    default:
      printf("LinearSystem_solve, context=%d\n",M->distributor->context);
      check_error(-1, "invalid solver context", __LINE__, __FILE__, 1);
      break;
    }
//  M->ordering->type=M->t.comptype;
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LinearSystem_initialise(hyperzmatrix_t *A, int *clamped, int nclamped, int solver_id)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int j, k, n;
  int reorder;
  int verbose=0;
  ordering_t *CSCordering;

/**----------------------------------------------------------------------
  set the matrix row to 0---0 1 0---0 at listed nodes*/
  for(k = 0; k < nclamped; k++) {
    n=clamped[k];
    for(j = A->ordering->pointer[n]; j < A->ordering->pointer[n + 1]; j++) {
      if(A->ordering->incidence[j] == n) {
        A->packed[j] = diagonal;
        }
      else {
        A->packed[j] = 0.0;
        }
      }
    }
//  status= matrix_diagonal(A);

  if(verbose) printf("Complex-valued linear system factorisation...\n");

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

  24/10/2010:

    KILLER : CSC expect column sparse matrix, we provide row sparse matrix
             up to now, this was overcome by using the transpose solver, it
             does work as our matrices are geometrically symmetric

             It does NOT work with complex-valued matrices as in addition to
             transposition, a conjugate operation is done

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  switch (solver_id) {
    case SOLVER_ID_DOMESTIC:
    case SOLVER_ID_LAPACK:
    case SOLVER_ID_SUNPERF:
    case SOLVER_ID_HIPS:
      break;
    default:
/**----------------------------------------------------------------------
      conversion to CSR ordering */
//       status=matrix_check_symmetry(*A->ordering);
//       status=matrix_CSR2CSC(A);
      status= matrix_check_symmetry(*A->ordering);
      if(status==0) {
/**----------------------------------------------------------------------
        if ordering is symmetric, CSR and CSC ordering are the same*/
        reorder=0;
        status= matrix_CSR2CSC(A);
        }
      else {
/**----------------------------------------------------------------------
        if ordering is not symmetric, we need a separate CSC ordering*/
        reorder=1;
//        status=matrix_check_Zero(*A);
        CSCordering=ordering_CSR2CSC(A->ordering);
        status= matrix_CSR2CSC(A,CSCordering);
        A->ordering->destroy();
        A->ordering=CSCordering;
        status= matrix_check_consistency(*CSCordering);
//        status=matrix_check_Zero(*A);
        }
    }

  if(A->distributor==0) {
#ifdef HAVE_MPI
    status = MPI_Barrier (  MPI_COMM_WORLD);
#endif
    check_error(-1, "node distribution structure not initialized", __LINE__, __FILE__, 1);
#ifdef HAVE_MPI
    status = MPI_Barrier (  MPI_COMM_WORLD);
#endif
    }

/**----------------------------------------------------------------------
  initialise triplet, temporary patch for PARALLEL runs */
  switch(A->distributor->context) {
    case SEQUENTIAL_COMPUTING:
      status=LinearSystem_InitTriplet(A,false,0,false);
      break;
    case PARALLEL_COMPUTING:
/**----------------------------------------------------------------------
      duplicate arrays, else sparse matrix operations garbled */
      status=LinearSystem_InitTriplet(A,true,0,false);
      break;
    default:
      printf("LinearSystem_factorization, context=%d\n",A->distributor->context);
      check_error(-1, "invalid solver context", __LINE__, __FILE__, 1);
      break;
    }
 
/**----------------------------------------------------------------------
  factorize matrix*/
  status = LinearSystem_factorisation(A, solver_id);

  if(verbose) printf("Complex-valued linear system factorisation status: %d\n", status);

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LinearSystem_terminate(hyperzmatrix_t *A)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int j, k, n;
  int verbose=0;
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

  30/12/2011: yet another killer

    KILLER : some solvers change the C numbering into F numbering
    this is the case of PASTIX...

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  if(A->s==0) {
    check_error(-1, "LinearSystem_terminate warning: try to free uniniatlised matrix", __LINE__, __FILE__, 0);
    return (-1);
    }

//  switch (gLinearSolverID) {
  switch (A->s->id) {
    case SOLVER_ID_DOMESTIC:
    case SOLVER_ID_LAPACK:
    case SOLVER_ID_SUNPERF:
      delete[] A->band;
      delete[] A->pivot;
      break;
    default:
      status=freez_solver(A->s);
    }

  if(A->ordering->incidence[0]==1) {
/**----------------------------------------------------------------------------
    change F-numbering (base 1) into C-numbering (base 0)*/
//    check_error(-1, "warning: numbering in ordering structure was changed by solver (pastix does...)", __LINE__, __FILE__, 0);
/**----------------------------------------------------------------------------
    decrement positions (block start)*/
    for (size_t n=0;n<A->ordering->nrows+1;n++) {
      A->ordering->pointer[n]--;
      }
/**----------------------------------------------------------------------------
    decrement node indices*/
    for (size_t n=0;n<A->ordering->pointer[A->ordering->nrows];n++) {
      A->ordering->incidence[n]--;
      }
    }
  A->ordering->type=ORDER_CSR;
  return (status);
}

