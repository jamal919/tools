
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

#ifdef PETSC
#include "petsc.h"
#include "petscksp.h"
#include "petscmat.h" 
#include "petscvec.h" 
#endif

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  communicator strategy
  
  warning : USE_DEFAULT_COMMUNICATOR=0 will hang on nuwa
  
  TODO : needs further testing and improvement to allow sequential solving in
         MPI T-UGOm
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#define USE_DEFAULT_COMMUNICATOR 1

int MyPetscInitialized=0;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  wrappers necessary to template use
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#ifdef PETSC

//   PetscErrorCode MatSetValues(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],const PetscScalar v[],InsertMode addv)
//   {
//   PetscErrorCode ierr;
//   ierr=MatSetValues(mat, m, idxm, n,idxn, v, addv);
//   return(ierr);
//   }

#endif

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  petsc interface
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int petsc_factorize_template(SOLVER_t<T> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  PETSC (personal) interpretation :
 
  * global row/col numbers are NOT relative to caller native partitioning, but
    to PETSC ranges
   
  * PETSC_COMM_WORLD needs to be set at any access to avoid confusion
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
    
#ifdef PETSC
  int count, m, k, n, neq, size, rank;
  PetscErrorCode ierr;
  PetscInt Istart, Iend;
  MPI_Comm communicator;
      
  PetscInt d_nz=0, o_nz=0;
  PetscInt *d_nnz, *o_nnz;
  
  communicator=solver->communicator;
  
  petsc_t<T> *parameters = (petsc_t<T> *) solver->parameters;

#if USE_DEFAULT_COMMUNICATOR
#else
/*------------------------------------------------------------------------------
  handle reduced cpu set case (small system do) */
  if(solver->rank==-1) return(0);
#endif
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create MATRIX in Petsc format :init 

  PetscErrorCode  MatCreateSeqAIJ(MPI_Comm comm,PetscInt m,PetscInt n,PetscInt nz,const PetscInt nnz[],Mat *A)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
/*------------------------------------------------------------------------------
  vector of number of non-zero value per row */
  PetscInt *nnz;
  status=solver->packed->cardinal(nnz);
  
  neq=solver->packed->nrows;
  
//   MPI_Comm_size(communicator,&size);
//   MPI_Comm_rank(communicator,&rank);
  
  size=solver->size;
  rank=solver->rank;
  
  if(debug) printf("%s cpu=%d : neq=%d\n", __func__, rank,neq);
  
  switch (size) {
    case 1:
      ierr = MatCreateSeqAIJ(communicator, neq, neq, 0, nnz, &parameters->matrix);
      if(debug and ierr!=0) printf("%s cpu=%d : MatCreateSeqAIJ status=%d \n", __func__, rank, ierr);
      break;
    
    default:
      ierr=MatCreate(communicator, &parameters->matrix);
      if(debug and ierr!=0) printf("%s cpu=%d : MatCreate status=%d \n", __func__, rank, ierr);
      ierr=MatSetType(parameters->matrix, MATMPIAIJ);
      if(debug and ierr!=0) printf("%s cpu=%d : MatSetType status=%d \n", __func__, rank, ierr);
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        
      m n  number of local rows (or PETSC_DECIDE), number of local columns (or PETSC_DECIDE)
      M N  number of global rows (or PETSC_DECIDE), number of global columns (or PETSC_DECIDE)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
      PetscInt m, n, M, N;
      m=solver->packed->nrows;
      n=solver->packed->ncols;
      M=solver->packed->ngdof;
      N=solver->packed->ngdof;
      
//       ierr=MatSetSizes(parameters->matrix, m, n, M, N);
      
//       ierr=MatSetSizes(parameters->matrix, m, n, PETSC_DETERMINE , PETSC_DETERMINE );
      ierr=MatSetSizes(parameters->matrix, m, m, PETSC_DETERMINE , PETSC_DETERMINE );
//       ierr=MatGetSize(parameters->matrix, &M, &N);
//       if(debug) printf("%s cpu=%d : MatSetSizes m=%d n=%d M=%d N=%d\n", __func__, rank, m, n, M, N);
      
//       ierr=MatSetSizes(parameters->matrix, PETSC_DECIDE, PETSC_DECIDE, M, N); // PETSC_DECIDE will likely not fit tugo partitioning
      
      if(debug and ierr!=0) printf("%s cpu=%d : MatSetSizes m=%d n=%d M=%d N=%d\n", __func__, rank, m, n, M, N);
      if(debug and ierr!=0) printf("%s cpu=%d : MatSetSizes status=%d \n", __func__, rank, ierr);
      CHKERRQ(ierr);
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        
      d_nz  - number of nonzeros per row in DIAGONAL portion of local submatrix 
              (same value is used for all local rows)
      d_nnz - array containing the number of nonzeros in the various rows of 
              the DIAGONAL portion of the local submatrix (possibly different 
              for each row) or NULL (PETSC_NULL_INTEGER in Fortran), if d_nz 
              is used to specify the nonzero structure. The size of this array 
              is equal to the number of local rows, i.e 'm'. For matrices that 
              will be factored, you must leave room for (and set) the diagonal 
              entry even if it is zero.
      o_nz -  number of nonzeros per row in the OFF-DIAGONAL portion of local 
              submatrix (same value is used for all local rows).
      o_nnz - array containing the number of nonzeros in the various rows of 
              the OFF-DIAGONAL portion of the local submatrix (possibly different
              for each row) or NULL (PETSC_NULL_INTEGER in Fortran), if o_nz 
              is used to specify the nonzero structure. The size of this array
              is equal to the number of local rows, i.e 'm'.
              
      If the *_nnz parameter is given then the *_nz parameter is ignored 

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

      if(solver->packed->nrows==0) {
        ierr=MatMPIAIJSetPreallocation(parameters->matrix, 0, 0, 0, 0);
        }
      else {        
        status=solver->packed->set_nnz(d_nnz, o_nnz);
        if(status!=0) TRAP_ERR_EXIT(-1,"set_nnz failed\n");
//         ierr=MatSeqAIJSetPreallocation(parameters->matrix, d_nz, nnz);
        ierr=MatMPIAIJSetPreallocation(parameters->matrix, d_nz, d_nnz, o_nz, o_nnz);
        }
      if(debug and ierr!=0) printf("%s cpu=%d : MatMPIAIJSetPreallocation status=%d \n", __func__, rank, ierr);
      CHKERRQ(ierr);
//       MatSetOption(parameters->matrix,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_FALSE);
      break;
    }
  
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  PetscErrorCode MatGetOwnershipRange(Mat mat,PetscInt *m,PetscInt *n)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ierr = MatGetOwnershipRange(parameters->matrix, &Istart, &Iend);
  if(debug and ierr!=0) printf("%s cpu=%d : MatGetOwnershipRange status=%d Istart=%d Iend=%d \n", __func__, rank, ierr, Istart, Iend);
  CHKERRQ(ierr);
  
  MPI_Barrier(communicator);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  now built glob2glob (i.e. native partitioning to Petsc partitioning) global 
  row/col numbers correspondance
  
  each cpu can infer Petsc partitioning global row/col numbers for locally solved rows:
  
  Petsc glob = local + Istart
  
  Then cpus must exchange those numbers for non-locally solved columns
  
  only for local numbers in incidence
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#if 0
  PetscInt *petsc_global=new PetscInt[solver->packed->ncols];
  
  if(solver->packed->loc2glob!=0) {
//     exchange_t & mpi_exchange=solver->native_distributor->mpi_exchange;
    if(debug) printf("%s cpu=%d : exchange=%lu \n", __func__, rank, solver->native_exchange); 
    MPI_Barrier(communicator);
    exchange_t & mpi_exchange=*solver->native_exchange;
    if(debug) status=communications_summary(mpi_exchange, debug);
    MPI_Barrier(communicator);
    for(m = 0; m < neq; m++) {
      int row;
      row=solver->packed->loc2glob[m];  // native global
      petsc_global[m]=Istart+m;         // petsc global
      }
    if(debug) printf("%s cpu=%d : call exchange_atom \n", __func__, rank);
    status=exchange_atom(mpi_exchange, rank, size, petsc_global, petsc_global, false);
    if(debug) printf("%s cpu=%d : exchange_atom status=%d \n", __func__, rank, status);
/*------------------------------------------------------------------------------
    modify loc2glob; safe ? */
    for(m = 0; m < neq; m++) {
      int row;
      row=solver->packed->loc2glob[m];              // native global
      solver->packed->loc2glob[m]=petsc_global[m];  // petsc global
      }
    }
#endif
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  now built glob2glob (i.e. native partitioning to Petsc partitioning) global 
  row/col numbers correspondance
  
  each cpu can infer Petsc partitioning global row/col numbers for locally solved rows:
  
  Petsc glob = local + Istart
  
  Then cpus must exchange those numbers for non-locally solved columns
  
  only for globam numbers in incidence
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  PetscInt *petsc_global=0, *tmp=0;
  
  if(solver->packed->loc2glob!=0) {
    petsc_global=new PetscInt[solver->packed->ngdof];
    tmp=new PetscInt[solver->packed->ngdof];
    for(m = 0; m < solver->packed->ngdof ; m++) {
      petsc_global[m]=0;
      tmp[m]=0;
      }
    for(m = 0; m < neq; m++) {
      int row;
      row=solver->packed->loc2glob[m];  // native global
      tmp[row]=Istart+m;                // petsc global
      }
    status = MPI_Allreduce(tmp, petsc_global, solver->packed->ngdof, MPI_INT, MPI_SUM, communicator ); 
/*------------------------------------------------------------------------------
    modify loc2glob; safe ? */
    for(m = 0; m < neq; m++) {
      int row;
      row=solver->packed->loc2glob[m];                // native global
      solver->packed->loc2glob[m]=petsc_global[row];  // petsc global
      }
    delete[] tmp;
    }

  MPI_Barrier(communicator);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  insert values
  
  row-per-row insertion, more optimal

  PetscErrorCode MatSetValues(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],const PetscScalar v[],InsertMode addv)

  MatSetValues() uses 0-based row and column numbers in Fortran as well as in C.
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  size_t offset=0;
  for(m = 0; m < neq; m++) {
    PetscInt row;
    PetscInt cardinal=solver->packed->pointer[m + 1]-solver->packed->pointer[m];
    if(solver->packed->loc2glob==0) {
      row=m;
      ierr = MatSetValues(parameters->matrix, 1, &row, cardinal, &(solver->packed->incidence[offset]), &(solver->packed->x[offset]), INSERT_VALUES);
      if(debug and ierr!=0) {
        printf("%s cpu=%d : MatSetValues status=%d m=%d row=%d cardinal=%d d=%d o=%d, ", __func__, rank,  ierr, m, row, cardinal, d_nnz[m], o_nnz[m]);
        for(size_t k=solver->packed->pointer[m];k<solver->packed->pointer[m+1];k++) printf(" %d",solver->packed->incidence[k]);
        printf("\n");
        CHKERRQ( ierr);
        }
      }
    else {
      row=Istart+m;  // global for PETCS-MPI
/*------------------------------------------------------------------------------
      build petsc global indices */
      PetscInt *indices= new PetscInt[cardinal];
      int count=0;
      for(size_t k=solver->packed->pointer[m]; k<solver->packed->pointer[m + 1];k++) {
        int col=solver->packed->incidence[k];
        indices[count]=petsc_global[col];
        count++;
        }
      ierr = MatSetValues(parameters->matrix, 1, &row, cardinal, indices, &(solver->packed->x[offset]), INSERT_VALUES);
      if(debug and ierr!=0) {
        printf("%s cpu=%d : MatSetValues status=%d m=%d row=%d cardinal=%d d=%d o=%d, ", __func__, rank,  ierr, m, row, cardinal, d_nnz[m], o_nnz[m]);
        for(size_t k=0;k<cardinal;k++) printf(" %d",indices[k]);
        printf("\n");
        CHKERRQ( ierr);
        }
      delete[] indices;
      }
    offset+=cardinal;
    }
  if(debug) printf("%s cpu=%d : MatSetValues successfuly passed\n", __func__, rank);
  
  ptr_delete(petsc_global);
  
  ierr = MatAssemblyBegin(parameters->matrix, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  
  ierr = MatAssemblyEnd(parameters->matrix, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  
  if(debug) PetscPrintf(communicator,"matrix initialization done\n");

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create KSP
  
  PetscErrorCode  KSPCreate(MPI_Comm comm,KSP *inksp)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ierr = KSPCreate(communicator, &parameters->ksp);
  if(debug and ierr!=0) printf("KSPCreate status=%d \n", ierr);
  CHKERRQ(ierr);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  Set operators. Here the matrix that defines the linear system
  also serves as the preconditioning matrix.
  
  PetscErrorCode  KSPSetOperators(KSP ksp,Mat Amat,Mat Pmat)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ierr = KSPSetOperators(parameters->ksp, parameters->matrix, parameters->matrix);
  if(debug and ierr!=0) printf("KSPSetOperators status=%d \n", ierr);
  CHKERRQ(ierr);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  GMRES, conjugate gradient etc...

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   ierr = KSPSetType(parameters->ksp, KSPBCGS);
  ierr = KSPSetType(parameters->ksp, KSPGMRES); 
  CHKERRQ(ierr);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  pre-conditioner

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  PC pc;
  PCType type;
  
  ierr=KSPGetPC(parameters->ksp, &pc);
  switch (size) {
    case 1:
      ierr=PCSetType(pc, "ilu");
      break;
    
    default:
      ierr=PCSetType(pc, "bjacobi");
      break;
    }
     
  ierr=PCGetType(pc, &type);
  if(rank==0) printf("%s cpu=%d : PC type=%s\n", __func__, rank, type);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  int KSPSetTolerances(KSP ksp,PetscReal rtol,PetscReal atol,PetscReal dtol,int maxits)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  size_t ngdof=solver->packed->ngdof;
  double relative_tolerance= 1.e-2 / ((ngdof + 1) * (ngdof + 1));
  double absolute_tolerance= 1.e-3;
  ierr = KSPSetTolerances(parameters->ksp, relative_tolerance, absolute_tolerance, PETSC_DEFAULT, PETSC_DEFAULT);
  if(debug and ierr!=0) printf("%s cpu=%d : KSPSetTolerances status=%d\n", __func__, rank, ierr);
  CHKERRQ(ierr);
  
  /*ierr = KSPSetMonitor(parameters->ksp,KSPDefaultMonitor,PETSC_NULL, PETSC_NULL); */
  status=ierr;
   
  if(debug) PetscPrintf(communicator,"solver factorizing done\n");
  
  return(status);
#else
  status=-1;
  printf("Please compile poc-solvers library with -DPETSC \n");
//       TRAP_ERR_EXIT(-1, "please compile with -DPETSC\n");
  return(status);
#endif
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int petsc_factorize(SOLVER_t<double> *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=petsc_factorize_template(solver, verbose, debug);
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int petsc_factorize(SOLVER_t< complex<double> > *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
//   commented until versatile double/complex<double> implementation done
//   status=petsc_factorize_template(solver, verbose, debug); //LEFT OVER
  return(status); 
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int petsc_solve_template(SOLVER_t<T> *solver, T *RHS, T *x, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  distributed RHS/solution version
  
  make use of:
  
  PetscErrorCode  VecSetSizes(Vec v, PetscInt n, PetscInt N)

  PetscErrorCode  VecGetValues(Vec x,PetscInt ni,const PetscInt ix[],PetscScalar y[])
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  debug=false;
  status=-1;
  
#ifdef PETSC
//   Vec xiter, biter;
  PetscBool flg;
  PetscScalar v;
  PetscErrorCode ierr;
  PetscInt neq;
  MPI_Comm communicator;
  
#if USE_DEFAULT_COMMUNICATOR
#else
/*------------------------------------------------------------------------------
  handle reduced cpu set case (small system do) */
  if(solver->rank==-1) return(0);
#endif
  
  communicator=solver->communicator;
  
  petsc_t<T> *parameters = (petsc_t<T> *) solver->parameters;

  Vec & biter=parameters->RHS;
  Vec & xiter=parameters->x;
  
  neq   = solver->packed->nrows;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create RHS vector
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(biter==0) {
    ierr = VecCreate(communicator, &biter);
    CHKERRQ(ierr);
    ierr = VecSetSizes(biter, neq, PETSC_DECIDE);
    CHKERRQ(ierr);
    }
    
  switch (solver->size) {
    case 1:
      VecSetType(biter, VECSEQ);
      for(PetscInt i = 0; i < neq; i++) {
        v = RHS[i];
        VecSetValue(biter, i, v, INSERT_VALUES);
        }
      break;
    default:
      ierr = VecSetType(biter, VECMPI);
      ierr = VecSetValues(biter, neq, solver->packed->loc2glob, RHS, INSERT_VALUES);
      break;
    }      

  ierr = VecAssemblyBegin(biter);
  CHKERRQ(ierr); if(ierr!=0) return(-1);
  ierr = VecAssemblyEnd(biter);
  CHKERRQ(ierr); if(ierr!=0) return(-1);

//   commented until versatile double/complex<double> implementation done
//   int VecValid(Vec v,PetscTruth *flg)
//   VecValid(xiter, &flg); // LEFT OVER
//   if(!flg)
//     ierr = VecDuplicate(biter, &xiter);
//   CHKERRQ(ierr);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create solution vector
  
  if x==0 create xiter and duplicate biter at 1st solve. Then keep last solution
  as new prior one
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#define REINITIALIZE 0
  
  if(xiter==0) {
    if(x==0) {
/*------------------------------------------------------------------------------
      no prior given, allocate vector and use RHS values */
      ierr = VecDuplicate(biter, &xiter);
      }
    else {
/*------------------------------------------------------------------------------
      prior given, allocate vector */
      ierr = VecCreate(communicator, &xiter);
      CHKERRQ(ierr);
      ierr = VecSetSizes(xiter, neq, PETSC_DECIDE);
      CHKERRQ(ierr);
      }
    }
    
  if(x!=0) {
/*------------------------------------------------------------------------------
    prior given, set values */
    ierr = VecSetType(xiter, VECMPI);
    ierr = VecSetValues(xiter, neq, solver->packed->loc2glob, x, INSERT_VALUES);
    ierr = VecAssemblyBegin(xiter);
    CHKERRQ(ierr); if(ierr!=0) return(-1);
    ierr = VecAssemblyEnd(xiter);
    CHKERRQ(ierr); if(ierr!=0) return(-1);
    }
  else {
#if REINITIALIZE
    ierr = VecSetType(xiter, VECMPI);
    ierr = VecSetValues(xiter, neq, solver->packed->loc2glob, RHS, INSERT_VALUES);
    ierr = VecAssemblyBegin(xiter);
    CHKERRQ(ierr); if(ierr!=0) return(-1);
    ierr = VecAssemblyEnd(xiter);
    CHKERRQ(ierr); if(ierr!=0) return(-1);
#else
/*------------------------------------------------------------------------------
    no prior given, keep last solution as prior */
#endif
    }
    
  ierr = KSPSolve(parameters->ksp, biter, xiter);
  CHKERRQ(ierr); if(ierr!=0) return(-1);
          
  if(debug and ierr!=0) printf("KSPSolve ok ierr=%d\n", ierr);
    
  switch (solver->size) {
    case 1:
      for(PetscInt i = 0; i < neq; i++) {
        ierr=VecGetValues(xiter, 1, &i, &v);
        CHKERRQ(ierr);
        RHS[i] = v;
        }
      break;
    default:
      ierr=VecGetValues(xiter, neq, solver->packed->loc2glob, RHS);
      CHKERRQ(ierr);
      break;
    }
    
//   if(x==0) {    
//     ierr = VecDestroy(&xiter);
//     CHKERRQ(ierr);
//     xiter=0;
//     }
    
  if(debug) PetscPrintf(communicator,"petsc solve done\n");

  return(0);
#endif
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int petsc_solve(SOLVER_t<double> *solver, double *RHS, double *x, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=petsc_solve_template(solver, RHS, (double *) 0, transposed, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int petsc_solve(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
//   commented until versatile double/complex<double> implementation done
//   status=petsc_solve_template(solver, RHS, 0, transposed, verbose, debug); //LEFT OVER
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int petsc_initialize_template(SOLVER_t<T> & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  MPI_Comm communicator;
  int rank;
  
#ifdef PETSC
    
  int ierr, size=0, argc=0;
  char **args;
  const char help[]="";
  
#if USE_DEFAULT_COMMUNICATOR
  solver.communicator=MPI_COMM_WORLD; // TESTING
#endif

  communicator=solver.communicator;
 
  PETSC_COMM_WORLD=communicator;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  PetscErrorCode  PetscInitialize(int *argc,char ***args,const char file[],const char help[])
  
  it is not clear if PETSC_COMM_WORLD should be limited to solver communicator, as
  simultaneous use of various Petsc solver invokation would use possibly different
  communicator
  
  I inferred that it is actually dealt with in matrix and ksp creation (factorization)
  and keeping MPI_COMM_WORLD as default communicator for initialization is ok
  
  warning : USE_DEFAULT_COMMUNICATOR=0 will hang on nuwa
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   ierr=PetscInitialize(&argc,&args,(char*)0,help);
//   CHKERRQ(ierr);
//   if (ierr) return ierr;
  
  PETSC_COMM_WORLD=MPI_COMM_WORLD;
  if(MyPetscInitialized==0) {
    ierr=PetscInitialize(&argc,&args,(char*)0,help);
    CHKERRQ(ierr);
    if (ierr) return ierr;
    MyPetscInitialized=1;
    }
  
//   if(solver.rank!=-1) {
//     communicator=solver.communicator;
//     MPI_Comm_size(communicator,&size);
//     }
    
  petsc_t<T> *parameters=new petsc_t<T>;
  solver.parameters = parameters;
  
  solver.targeted_packing    = PACKED;
  solver.targeted_ordering   = CSR;
  solver.targeted_addressing = MATRIX_POINTER;
  solver.targeted_numbering  = 0;
  
  solver.initialize = petsc_initialize;
  solver.terminate  = petsc_terminate;
  solver.factorize  = petsc_factorize;
  solver.solve      = petsc_solve;
  
  solver.id=SOLVER_ID_PETSC_GMRES;

  if(debug) printf("%s cpu=%d (in calling group of %d) : communicator=%d (MPI_COM_WORLD=%d)\n",__FUNCTION__,solver.rank, solver.size,solver.communicator, MPI_COMM_WORLD);
  
  return(0);
#else 
    printf("Please compile poc-solvers library with -DPETSC \n");
    exit(-1);
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int petsc_initialize(SOLVER_t<double> & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=petsc_initialize_template(solver, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int petsc_initialize(SOLVER_t< complex<double> > & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=petsc_initialize_template(solver, verbose, debug);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int petsc_terminate_template(SOLVER_t<T> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  MPI_Comm communicator;
  
#ifdef PETSC
  
#if USE_DEFAULT_COMMUNICATOR
#else
/*------------------------------------------------------------------------------
  handle reduced cpu set case (small system do) */
  if(solver->rank==-1) return(0);
#endif
  
  petsc_t<T> *parameters = (petsc_t<T> *) solver->parameters;
  PetscErrorCode ierr;
  
  ierr = MatDestroy(&parameters->matrix);
  CHKERRQ(ierr);
  ierr = KSPDestroy(&parameters->ksp);
  CHKERRQ(ierr);
  
  ierr = VecDestroy(&parameters->x);
  CHKERRQ(ierr);
  ierr = VecDestroy(&parameters->RHS);
  CHKERRQ(ierr);
  
  status=0;
  
  return(status);
#else
  return(status);
#endif 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int petsc_terminate(SOLVER_t<double> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=petsc_terminate_template(solver, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int petsc_terminate(SOLVER_t< complex<double> > *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=petsc_terminate_template(solver, verbose, debug);
  
  return(status);
}


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  solver interface compatibility with older versions
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#if OLD_SOLVER_VERSION
#else


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorize_petsc(solver_t *slv, triplet *Mat, int verbose)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t<double> solver;

  bool debug=false;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);

  slv->Mat = Mat; 
      
  solver.d_import2packed(slv);

  status=petsc_factorize(&solver, verbose, debug);
  
  solver.d_export_packed(slv);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorizez_petsc(solver_t *slv, tripletz *Mat, int verbose)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t< complex<double> > solver;
 
  bool debug=false;

  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  slv->Mat = Mat;
  
  solver.z_import2packed(slv);

  status=petsc_factorize(&solver, verbose, debug);
  
  solver.z_export_packed(slv);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solve_petsc(solver_t *slv, double *RHS, int transposed)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int verbose=0;
  bool debug=false;
  
  SOLVER_t<double> solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  solver.d_import2packed(slv);

  status=petsc_solve(&solver, RHS, 0, transposed, verbose, debug);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solvez_petsc(solver_t *slv, complex<double> *RHS, int transposed)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t< complex<double> > solver;
 
  int verbose=0;
  bool debug=false;
  
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  solver.z_import2packed(slv);

  status=petsc_solve(&solver, RHS, 0, transposed, verbose, debug);
  
  return(status); 
}

#endif
