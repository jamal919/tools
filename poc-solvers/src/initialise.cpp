
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

/**-----------------------------------------------------------------------------
  added because of PASTIX/HIPS... initilisation */
  int MPI_Init_thread_DONE=0;
  
  int HIPSSOLVERNB =0;
  int HIPSSOLVERNBZ=0;
  
  static int *DHIPSSOLVER=0, *ZHIPSSOLVER=0;

  int nbsyslinearz2=0;
  
  
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  triplet : i,j,value (row-arranged) or j,i,value (column-arranged)
  
  CSR : value array, column index array, pointer array (optionally cardinal array)
  CSC : value array,    row index array, pointer array (optionally cardinal array)

             natural ordering    natural format    natural numbering  
              
  PASTIX   compressed columns               CSC              FORTRAN
  UMFPACK  compressed columns               CSC                    C
  CXSPARSE compressed columns               CSC                    C
  MUMPS       compressed rows     row-major COO              FORTRAN
  HIPS        compressed rows     row-major COO            C/FORTRAN (HIPS_FORTRAN_NUMBERING switch)
  MAPHYS      


  T-UGOm natural format is CSR (compressed rows ordering using pointers to 
  column indices for each row)
  
  Issues:
  
  - MPI distributed pastix kept in CSC, then transpose system solved
  
  - SCOTCH found to be sensitive to partitions contiguousity (in pastix)
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Some OpenMPI distributions do not support MPI_THREAD_MULTIPLE, that is wanted
  by default by SCOTCH when invoked by PASTIX. You can check that by launching:
  
  >ompi_info | grep -i thread
  
  In that case, you have 2 options :
   - change your MPI library for one supporting MPI_THREAD_MULTIPLE 
   - Rebuild Scotch without SCOTCH_PTHREAD and PaStiX with FORCE_NOSMP

 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  We have encountered a SMP issue with umfpack (requesting libopenblas_pthreads)
  libumfpack-5.7.6.so version, detected in August 2017

  issue : multi-threading disabled
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  MPI RHS/SOLUTION : distributed, centralized (on master proc), on demand (distributed/centralized)

             MPI             RHS      SOLUTION  
             
  PASTIX     yes     distributed   distributed
  UMFPACK     no               -             -
  CXSPARSE    no               -             -
  MUMPS      yes     centralized     on demand 
  HIPS       yes       on demand     on demand
  MAPHYS     yes       on demand     on demand

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void  set_MPI_Init_thread_DONE(int done) 
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  MPI_Init_thread_DONE=done;
}
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ishere(const char *solver_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{  
  if (strcmp(solver_name,"MUMPS") == 0) {
#ifdef MUMPS
    return(1);
#else
    return(0);
#endif
    }
  else if (strcmp(solver_name,"UMFPACK") == 0){
#ifdef UMFPACK
    return(1);
#else
    return(0);
#endif
    }
  else if (strcmp(solver_name,"HIPS") == 0){
#ifdef HIPS
    return(1);
#else
    return(0);
#endif
    }
  else if (strcmp(solver_name,"LAPACK") == 0){
#ifdef LAPACK
    return(1);
#else
    return(0);
#endif
    }
  else if (strcmp(solver_name,"HIPS") == 0){
#ifdef HIPS
    return(1);
#else
    return(0);
#endif
    }
  else if (strcmp(solver_name,"PASTIX") == 0){
#ifdef PASTIX
    return(1);
#else
    return(0);
#endif
    }
  else if (strcmp(solver_name,"MAPHYS") == 0){
#ifdef MAPHYS
    return(1);
#else
    return(0);
#endif
    }
  else if (strcmp(solver_name,"SpDOMESTIC") == 0){
    return(1);
    }
  else {
    return(0);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int get_SolverId(const char *solver_name, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int id;
  
  id = -1;
  
  if(strcmp(solver_name, "DOMESTIC")          == 0) 
    id = SOLVER_ID_DOMESTIC;
  else if(strcmp(solver_name, "LAPACK")       == 0)
    id = SOLVER_ID_LAPACK;
  else if(strcmp(solver_name, "SUNPERF")      == 0)
    id = SOLVER_ID_SUNPERF;
  else if(strcmp(solver_name, "UMFPACK")      == 0)
    id = SOLVER_ID_UMFPACK;
  else if(strcmp(solver_name, "PETSC_GMRES")  == 0)
    id = SOLVER_ID_PETSC_GMRES;
  else if(strcmp(solver_name, "MUMPS")        == 0)
    id = SOLVER_ID_MUMPS;
  else if(strcmp(solver_name, "MUMPS_SYM")    == 0)
    id = SOLVER_ID_MUMPS_SYM;
  else if(strcmp(solver_name, "SpDOMESTIC")   == 0)
    id = SOLVER_ID_SpDOMESTIC;
  else if(strcmp(solver_name, "PASTIX")       == 0)
    id = SOLVER_ID_PASTIX;
  else if(strcmp(solver_name, "HIPS")         == 0)
    id = SOLVER_ID_HIPS;
  else if(strcmp(solver_name, "TEST1")        == 0)
    id = 99;
  else if(strcmp(solver_name, "TEST2")        == 0)
    id = 98;
/*------------------------------------------------------------------------------
  temporary, until dynamic solver handling fully implemented */
  else if(strcmp(solver_name, "TUGO_GMRES")   == 0)
    id = SOLVER_ID_TUGO_GMRES;
  
  if(debug) {
    printf("%s : solver name=%s id=%d\n", __func__, solver_name, id);
    }

  return id;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LinearSystem_identify(const char *solver_name, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int id;
  
  id=get_SolverId(solver_name, debug);
  
  return id;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *get_ordering_name(int ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  char *name;

  switch (ordering) {
    case CSR:
      name=strdup(CSR_NAME);
      break;

    case CSC:
      name=strdup(CSC_NAME);
      break;

    case COO:
      name=strdup(COO_NAME);
      break;

    default:
      name=strdup("unknown");
      break;
    }
    
  return(name);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *get_SolverName(int solver_id)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  char *name;

  switch (solver_id) {
    case SOLVER_ID_DIAGONAL:
      name=strdup(SOLVER_NAME_DIAGONAL);
      break;

    case SOLVER_ID_DOMESTIC:
      name=strdup(SOLVER_NAME_DOMESTIC);
      break;

    case SOLVER_ID_LAPACK:
      name=strdup(SOLVER_NAME_LAPACK);
      break;

    case SOLVER_ID_SUNPERF:
      name=strdup(SOLVER_NAME_SUNPERF);
      break;

    case SOLVER_ID_UMFPACK:
      name=strdup(SOLVER_NAME_UMFPACK);
      break;

    case SOLVER_ID_PETSC_GMRES:
      name=strdup(SOLVER_NAME_PETSC_GMRES);
      break;

/*----------------------------------------------------------------------
    MUMPS ASSYMETRIC*/
    case SOLVER_ID_MUMPS:
      name=strdup(SOLVER_NAME_MUMPS);
      break;
 
/*----------------------------------------------------------------------
    MUMPS SYMETRIC*/
    case SOLVER_ID_MUMPS_SYM:
      name=strdup(SOLVER_NAME_MUMPS_SYM);
      break;

    case SOLVER_ID_HIPS:
      name=strdup(SOLVER_NAME_HIPS);
      break;

    case SOLVER_ID_HYPRE:
      name=strdup(SOLVER_NAME_HYPRE);
      break;

    case SOLVER_ID_SpDOMESTIC: 
      name=strdup(SOLVER_NAME_SpDOMESTIC);
      break;

    case SOLVER_ID_PASTIX:
      name=strdup(SOLVER_NAME_PASTIX);
      break;

    case SOLVER_ID_MAPHYS:
      name=strdup(SOLVER_NAME_MAPHYS);
      break;

/*------------------------------------------------------------------------------
    temporary, until dynamic solver handling fully implemented */
    case SOLVER_ID_TUGO_GMRES:
      name=strdup(SOLVER_NAME_TUGO_GMRES);
      break;

    default:
      return(0);
      break;
    }
    
  return(name);
}


/*************************************************************************/
/*======================= FIN INFOS ==================================== */
/*=======================================================================*/
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int DHIPSinitialise_IDs()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n;
  if(DHIPSSOLVER==0) {
    DHIPSSOLVER=(int *)  malloc(100*sizeof(int));
    for(n=0;n<100; n++) DHIPSSOLVER[n]=-1;
    }
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ZHIPSinitialise_IDs()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n;
  if(ZHIPSSOLVER==0) {
    ZHIPSSOLVER=(int *) malloc(100*sizeof(int));
    for(n=0;n<100; n++) ZHIPSSOLVER[n]=-1;
    }
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int DHIPSget_IDs()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int id=0;
  while (DHIPSSOLVER[id]!=-1) id++;
  DHIPSSOLVER[id]=id;
  return(id);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int DHIPSfree_IDs(int id)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  DHIPSSOLVER[id]=-1;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ZHIPSget_IDs()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int id=0;
  while (ZHIPSSOLVER[id]!=-1) id++;
  ZHIPSSOLVER[id]=id;
  
  return(id);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ZHIPSfree_IDs(int id)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  ZHIPSSOLVER[id]=-1;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
 
  int init_solverstruct(solver_t *solver)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
   solver->name=0;
   solver->parameters=0;
   solver->Mat=0;
   solver->RHS_distribution=solver->solution_distribution=solver->matrix_distribution=-1;
   solver->transpose_matrix=-1;
   solver->format=-1;
   solver->communicator=0;
   solver->rank=0;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  solver_t *init_solver_obsolete(const char *solver_name,int typsym, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  solver_t *slv;
  int ierr,myid;
  int status;
 
#ifdef MUMPS 
  DMUMPS_STRUC_C *par_mumps;
#endif
#ifdef UMFPACK
  UMFPACK_STRUC_C *par_umf; 
  double *Control;
#endif
#ifdef HYPRE
  HYPRE_STRUC_C *par_hypre; 
#endif
#ifdef HIPS
  dHIPS_STRUC_C *par_hips; 
#endif
#ifdef PASTIX
  PASTIX_STRUC_C *par_pastix; 
#endif
#ifdef MAPHYS
  MAPHYS_STRUC_C *id_maphys; 
  MAPHYS_IJMatrix A;
  MAPHYS_ParCSRMatrix parcsr_A;
  MAPHYS_IJVector b;
  MAPHYS_ParVector par_b;
  MAPHYS_IJVector x;
  MAPHYS_ParVector par_x;
  MAPHYS_Solver maphys_solver;
#endif
  
  int argc=0;
  char **argv=0;
  int  required;    /* MPI thread level required       */
  int  provided;    /* MPI thread level provided       */
  int  rank;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  solver structure default initialisation
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  slv = (solver_t *) malloc(sizeof(solver_t));
  slv->name = strdup(solver_name);
  
/*------------------------------------------------------------------------------
  default setting, possibly set to 1 in case of PASTIX solver */
  slv->transpose_matrix=0;
 
  slv->parameters =0;
  
  slv->RHS_distribution=SOLVER_CENTRALIZED;
  slv->solution_distribution=SOLVER_CENTRALIZED;
  
  slv->communicator=MPI_COMM_WORLD;
  slv->rank=-1;
  

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT MUMPS SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  if (strcmp(solver_name,"MUMPS") == 0) {
#ifdef MUMPS
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
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    slv->parameters = (DMUMPS_STRUC_C *) calloc(1,sizeof(DMUMPS_STRUC_C));
    par_mumps = (DMUMPS_STRUC_C *) slv->parameters;
    /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
    /*(*par_mumps).job=JOB_INIT; (*par_mumps).par=1; (*par_mumps).sym=typsym;*/
    par_mumps->job=JOB_INIT; par_mumps->par=1; par_mumps->sym=0;
    par_mumps->comm_fortran=USE_COMM_WORLD;
    dmumps_c(par_mumps);
    status=par_mumps->info[0];
    if(status==-1)  printf("cpu %d: DMUMPS factorisation failed, status=%d due to cpu %d\n",rank, status, par_mumps->info[1]);
    if(verbose==1 && myid==0) printf("dmumps_par->instance_number %d \n",par_mumps->instance_number);
    if(verbose==1) printf("init MUMPS solver (cpu=%d), finish \n",myid);
#else 
    printf("Please compile poc-solvers library with -DMUMPS \n");
    exit(-1);
#endif
    }


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT UMFPACK SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"UMFPACK") == 0){
#ifdef UMFPACK
    slv->parameters = (UMFPACK_STRUC_C *) malloc(sizeof(UMFPACK_STRUC_C));
    par_umf = (UMFPACK_STRUC_C *) slv->parameters;
    Control=par_umf->Control;
    umfpack_di_defaults (par_umf->Control);
/*     Control [UMFPACK_PRL] = 6 ; */
/*     Control [UMFPACK_PRL] = 5 ; */
/*     umfpack_di_report_control (Control) ; */
#else 
    printf("Please compile poc-solvers library with -DUMFPACK \n");
    exit(-1);
#endif
  }  

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT HYPRE SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"HYPRE") == 0){
#ifdef HYPRE
    slv->parameters = (HYPRE_STRUC_C *) malloc(sizeof(HYPRE_STRUC_C));
    par_hypre = slv->parameters;
#else 
    printf("Please compile poc-solvers library with -DHYPRE \n");
    exit(-1);
#endif
  }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT HIPS SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"HIPS") == 0){
#ifdef HIPS
    slv->parameters = (dHIPS_STRUC_C *) malloc(sizeof(dHIPS_STRUC_C));
    status=DHIPSinitialise_IDs();
    par_hips = (dHIPS_STRUC_C *) slv->parameters;
    par_hips->unknownlist=0;
    required = MPI_THREAD_MULTIPLE;
    provided = -1;
    if(MPI_Init_thread_DONE==0) {
/*------------------------------------------------------------------------------
      solver used in a non-parallel context, provide minimal MPI initialisation */
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
    if (HIPSSOLVERNB == 0) {
      int idnbrmax = 100;
      ierr = dHIPS_Initialize(idnbrmax);
      dHIPS_ExitOnError(ierr);
      }
//     par_hips->id_hips = HIPSSOLVERNB; /** id_hips of the linear system **/
    par_hips->id_hips = DHIPSget_IDs(); /** id_hips of the linear system **/
    HIPSSOLVERNB++;
    ierr = dHIPS_SetCommunicator(par_hips->id_hips, MPI_COMM_WORLD);
    if(ierr!=1) {
      dHIPS_ExitOnError(ierr);
      }
    slv->communicator=MPI_COMM_WORLD;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    slv->rank=rank;
    slv->RHS_distribution=SOLVER_CENTRALIZED;
    slv->solution_distribution=SOLVER_CENTRALIZED;
#else
    printf("Please compile poc-solvers library with -DHIPS \n");
    exit(-1);
#endif
  }

  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT PASTIX SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"PASTIX") == 0){
#ifdef PASTIX
    slv->parameters = (PASTIX_STRUC_C *) calloc(1,sizeof(PASTIX_STRUC_C));
    par_pastix = (PASTIX_STRUC_C *) slv->parameters;
#ifdef HAVE_MPI
    if(verbose) printf("init_solver: PASTIX parallel mpi version \n");
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
#endif
#else 
    printf("%s (line %d) : please compile poc-solvers library with -DPASTIX \n", __func__,__LINE__);
    exit(-1);
#endif
    } 
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT PETSC SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"PETSC") == 0){
    printf("PETSC pas encore implimente\n");
    exit(-1);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT LAPACK SOLVER 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"LAPACK") == 0){
    slv->parameters = NULL;
#ifdef LAPACK
#else
    printf("Please compile poc-solvers library with -DLAPACK \n");
    exit(-1);
#endif
    }
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT SpDOMESTIC SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"SpDOMESTIC") == 0){
    slv->parameters = NULL;   
    }

  else {
    printf("Connait pas ce solveur : %s \n",solver_name);
    slv = NULL;
    exit(-1);
    }
    
//   free(solver_name);

  return(slv);
}

#if OLD_SOLVER_VERSION
// if activated, use historical (1.0) code (still available in the c-code branch)
// not maintained, kept only for testing purposes, aimed to be removed

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  solver_t *init_solver(int solver_id, MPI_Comm world, int nprocs, int nunknowns, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  solver_t *slv;
  int ierr;
  int status;
  
#ifdef MUMPS
  DMUMPS_STRUC_C *par_mumps;
#endif
#ifdef UMFPACK
  UMFPACK_STRUC_C *par_umf; 
  double *Control;
#endif
#ifdef HYPRE
  HYPRE_STRUC_C *par_hypre; 
#endif
#ifdef HIPS
  dHIPS_STRUC_C *par_hips; 
#endif
#ifdef PASTIX
  PASTIX_STRUC_C *par_pastix; 
#endif
  
  int  argc=0;
  char **argv=0;
  int  required;    /* MPI thread level required       */
  int  provided;    /* MPI thread level provided       */
  int  rank, size;
  

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  solver structure default initialisation
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  slv = (solver_t *) malloc(sizeof(solver_t));
  
  char *solver_name=get_SolverName(solver_id);
  slv->name = strdup(solver_name);
  
/*------------------------------------------------------------------------------
  default setting, possibly set to 1 in case of PASTIX solver */
  slv->transpose_matrix=0;
 
  slv->parameters =0;
  
  slv->RHS_distribution=SOLVER_CENTRALIZED;
  slv->solution_distribution=SOLVER_CENTRALIZED;
  
  slv->communicator=0;
  
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  minimal MPI initialisation for sequential use
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#ifdef HAVE_MPI
  switch(solver_id) {
    case SOLVER_ID_MUMPS:
    case SOLVER_ID_MUMPS_SYM:
    case SOLVER_ID_HIPS:
    case SOLVER_ID_PASTIX:
    case SOLVER_ID_MAPHYS:
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
  
  ierr = MPI_Comm_rank(world, &rank);
  ierr = MPI_Comm_size(world, &size);
  
  if(debug) printf("%s cpu %d (in calling group of %d) : start\n",__FUNCTION__,rank, size);  

  int color=0, key=0;
  MPI_Comm newcomm;
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
    ierr = MPI_Comm_dup(world, &newcomm);
    ierr = MPI_Comm_rank(newcomm, &rank);
    }
  else {
    if(nunknowns==0) color=MPI_UNDEFINED;
    ierr = MPI_Comm_split(world, color, key, &newcomm);
    if(nunknowns!=0) ierr = MPI_Comm_rank(world, &rank);
    }
  
  slv->communicator=newcomm;
  slv->rank=rank;
#else
  slv->rank=0;
#endif
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT MUMPS SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  if (strcmp(solver_name,"MUMPS") == 0) {
#ifdef MUMPS
    slv->parameters = (DMUMPS_STRUC_C *) calloc(1,sizeof(DMUMPS_STRUC_C));
    par_mumps = (DMUMPS_STRUC_C *) slv->parameters;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
    /*(*par_mumps).job=JOB_INIT; (*par_mumps).par=1; (*par_mumps).sym=typsym;*/
    par_mumps->job=JOB_INIT; 
    par_mumps->par=1; 
    par_mumps->sym=0;
    par_mumps->comm_fortran=USE_COMM_WORLD;
    dmumps_c(par_mumps);
    if(verbose==1 && rank==0) printf("dmumps_par->instance_number %d \n",par_mumps->instance_number);
    if(verbose==1) printf("init MUMPS solver (cpu=%d), finish \n",rank);
#else 
    printf("Please compile poc-solvers library with -DMUMPS \n");
    exit(-1);
#endif
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT UMFPACK SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"UMFPACK") == 0) {
#ifdef UMFPACK
    slv->parameters = (UMFPACK_STRUC_C *) malloc(sizeof(UMFPACK_STRUC_C));
    par_umf = (UMFPACK_STRUC_C *) slv->parameters;
    Control=par_umf->Control;
    umfpack_di_defaults (par_umf->Control);
/*     Control [UMFPACK_PRL] = 6 ; */
/*     Control [UMFPACK_PRL] = 5 ; */
/*     umfpack_di_report_control (Control) ; */
#else 
    printf("Please compile poc-solvers library with -DUMFPACK \n");
    exit(-1);
#endif
    }  

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT HYPRE SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"HYPRE") == 0) {
#ifdef HYPRE
    slv->parameters = (HYPRE_STRUC_C *) malloc(sizeof(HYPRE_STRUC_C));
#else 
    printf("Please compile poc-solvers library with -DHYPRE \n");
    exit(-1);
#endif
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT HIPS SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"HIPS") == 0) {
#ifdef HIPS
    slv->parameters = (dHIPS_STRUC_C *) malloc(sizeof(dHIPS_STRUC_C));
    status=DHIPSinitialise_IDs();
    par_hips = (dHIPS_STRUC_C *) slv->parameters;
    par_hips->unknownlist=0;
    if (HIPSSOLVERNB == 0) {
      int idnbrmax = 100;
      ierr = dHIPS_Initialize(idnbrmax);
      dHIPS_ExitOnError(ierr);
      }
    par_hips->id_hips = DHIPSget_IDs(); /** id_hips of the linear system **/
    HIPSSOLVERNB++;
    if(rank!=-1) {
      ierr = dHIPS_SetCommunicator(par_hips->id_hips, newcomm);
      dHIPS_ExitOnError(ierr);
      }
#else
    printf("Please compile poc-solvers library with -DHIPS \n");
    exit(-1);
#endif
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT PASTIX SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"PASTIX") == 0) {
#ifdef PASTIX
    slv->parameters = (PASTIX_STRUC_C *) calloc(1,sizeof(PASTIX_STRUC_C));
    par_pastix = (PASTIX_STRUC_C *) slv->parameters;
#else 
    printf("%s (line %d) : please compile poc-solvers library with -DPASTIX \n", __func__,__LINE__);
    exit(-1);
#endif
    } 
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT PETSC SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"PETSC") == 0){
    slv->parameters = NULL;   
    printf("PETSC pas encore implimente\n");
    exit(-1);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT LAPACK SOLVER 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"LAPACK") == 0){
    slv->parameters = NULL;
#ifdef LAPACK
#else
    printf("Please compile poc-solvers library with -DLAPACK \n");
    exit(-1);
#endif
    }
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT SpDOMESTIC SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"SpDOMESTIC") == 0){
    slv->parameters = NULL;   
    }

  else {
    printf("Connait pas ce solveur : %s \n",solver_name);
    slv = NULL;
    exit(-1);
    }
    
  free(solver_name);

  return(slv);
}

#endif

/*=======================================================================*/
/*                 Solver Terminate                                      */
/*************************************************************************/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int free_solver(solver_t *slv)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
#ifdef MUMPS
  DMUMPS_STRUC_C *id_mumps;
#endif
#ifdef UMFPACK
  UMFPACK_STRUC_C *id_umf;
#endif
#ifdef PASTIX
  PASTIX_STRUC_C *id_pastix;
  pastix_data_t  *pastix_data = NULL;
  pastix_int_t    *iparm;  /* integer parameters for pastix                             */
  double          *dparm;  /* floating parameters for pastix                            */
  double *Ax;
  int *Ap, *Ai,n,nnz,nbthread;
  pastix_int_t   *perm        = NULL; /* Permutation tabular                             */
  pastix_int_t   *invp        = NULL; /* Reverse permutation tabular                     */
  double *rhs         = NULL; /* right hand side                                         */
#endif  /*double *Control;*/
#ifdef HIPS
  dHIPS_STRUC_C* parameters;
  INTS id_hips;
  int ierr;
#endif
  char *solver_name;
  void *Symbolic, *Numeric ;
  triplet *Mat;
  FILE *out;
  MPI_Comm world;

  solver_name=slv->name;
  Mat = (triplet *) slv->Mat;

  if (solver_name != NULL) {
/* ---------------- FREE MUMPS SOLVER -------------------- */
#ifdef MUMPS
    if (strcmp(solver_name,"MUMPS") == 0) {
      id_mumps = (DMUMPS_STRUC_C *) slv->parameters;
      (*id_mumps).job=JOB_END;
      dmumps_c(id_mumps);
//       (*id_mumps).irn_loc=0;
//       (*id_mumps).jcn_loc=0;
//       (*id_mumps).a_loc=0;
      Mat->i=0;
      Mat->j=0;
      Mat->x=0;
      }
#endif

/* ---------------- FREE UMFPACK SOLVER -------------------- */
#ifdef UMFPACK
    if (strcmp(solver_name,"UMFPACK") == 0) {
      id_umf = (UMFPACK_STRUC_C *) slv->parameters;
      Symbolic=id_umf->Symbolic;
      Numeric=id_umf->Numeric;
      umfpack_di_free_symbolic (&Symbolic) ;
      umfpack_di_free_numeric (&Numeric) ;
      }
#endif

/* ---------------- FREE LAPACK SOLVER -------------------- */
#ifdef LAPACK
    if (strcmp(solver_name,"LAPACK") == 0) {
      free(Mat->a);
      free(Mat->pivot);
      }
#endif

/* ---------------- FREE PASTIX SOLVER -------------------- */
#ifdef PASTIX
    if (strcmp(solver_name,"PASTIX") == 0) {
      id_pastix = (PASTIX_STRUC_C *) slv->parameters;
      Ap=Mat->i;
      Ai=Mat->j;
      Ax=Mat->x;
      
      id_pastix->iparm[IPARM_START_TASK] = API_TASK_CLEAN;
      id_pastix->iparm[IPARM_END_TASK]   = API_TASK_CLEAN;
      
      world=slv->communicator;

      d_pastix(&(id_pastix->pastix_data), world,
         n, Ap, Ai, Ax,
         id_pastix->perm, id_pastix->invp, rhs, 1, 
         id_pastix->iparm, id_pastix->dparm);

      free(id_pastix->perm);
      free(id_pastix->invp);
      if(Mat->private_memory==1) {
        free(Ax);
        free(Ap);
        free(Ai);
        }
      }
#endif

#ifdef HIPS
    if (strcmp(solver_name,"HIPS") == 0) {
      /************************************************/
      /* Free zHIPS internal structure for problem id  */
      /************************************************/
      parameters  = (dHIPS_STRUC_C*) slv->parameters;
      id_hips =parameters->id_hips;
      ierr = dHIPS_Clean(id_hips);
      dHIPS_ExitOnError(ierr);
      DHIPSfree_IDs(id_hips);
/*------------------------------------------------------------------------------
      apparently hips communicator has been freed here */
      slv->rank=-1;
      HIPSSOLVERNB--;
      
      /************************************************/
      /* Free zHIPS internal structures               */
      /* (to be done once only)                       */
      /************************************************/
      if(HIPSSOLVERNB==0) {
        ierr = dHIPS_Finalize();
        dHIPS_ExitOnError(ierr);
        }
      if(parameters->unknownlist!=0) free(parameters->unknownlist);
      }
#endif
    
/* ---------------- FREE SpDOMESTIC SOLVER -------------------- */
    if (strcmp(solver_name,"SpDOMESTIC") == 0) {
      cs_sfree (Mat->S) ;
      cs_nfree (Mat->N) ;
      }
    
    free(slv->parameters);
    free(slv->name);
    /*free(slv);*/
#if HAVE_MPI
    if(slv->rank!=-1) status=MPI_Comm_free(&slv->communicator);
//     if(slv->communicator!=MPI_COMM_WORLD) status=MPI_Comm_free(&slv->communicator);
//     slv->communicator=0;
#endif
    slv->parameters = NULL;
    slv->name = NULL;
    }
  else {
    printf(" free already done with this solver \n");
    }
  
  return(status);
}

/*************************************************************************/
/*=======================================================================*/
/*======================FIN TERMINATE ===================================*/
/*=======================================================================*/



/*************************************************************************/
/*************************************************************************/
/*=============== Complex Double precision routines =============================*/
/*************************************************************************/
/*************************************************************************/

#if OLD_SOLVER_VERSION
// if activated, use historical (1.0) code (still available in the c-code branch)
// not maintained, kept only for testing purposes, aimed to be removed

solver_t *initz_solver(char *solver_name, int typsym, int verbose)
{
  solver_t *slv;
  int ierr,myid,idnbr;
  int status;
  
#ifdef MUMPS 
  ZMUMPS_STRUC_C *par_mumps;
  ZMUMPS_STRUC_C *id;
#endif

#ifdef UMFPACK
  UMFPACK_STRUC_C *id_umf; 
  double *Control;
#endif
  
#ifdef HIPS
/*------------------------------------------------------------------------------
          poc-solvers refurbishing : bug found */
//   dHIPS_STRUC_C *par_hips; 
  zHIPS_STRUC_C *par_hips; 
#endif
  
#ifdef PASTIX
  PASTIX_STRUC_C *id_pastix; 
#endif
  
  int argc;
  char **argv;
  int             required;           /* MPI thread level required                                 */
  int             provided;           /* MPI thread level provided                                 */
  int             rank;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  solver structure default initialisation
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  slv = (solver_t *) malloc(sizeof(solver_t));
  slv->name = strdup(solver_name);
 
/*------------------------------------------------------------------------------
  default setting, possibly set to 1 in case of PASTIX solver */
  slv->transpose_matrix=0;
 
  slv->parameters =0;

  slv->RHS_distribution=SOLVER_CENTRALIZED;
  slv->solution_distribution=SOLVER_CENTRALIZED;

  slv->communicator=MPI_COMM_WORLD;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT MUMPS SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  if (strcmp(solver_name,"MUMPS") == 0) {
#ifdef MUMPS
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
    ierr = MPI_Barrier(MPI_COMM_WORLD);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    slv->parameters = (ZMUMPS_STRUC_C *) calloc(1,sizeof(ZMUMPS_STRUC_C));
    id =  (ZMUMPS_STRUC_C *) slv->parameters;
    /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
    /*(*id_mumps).job=JOB_INIT; (*id_mumps).par=1; (*id_mumps).sym=typsym;*/
    id->job=JOB_INIT; id->par=1; id->sym=0;
    id->comm_fortran=USE_COMM_WORLD;
    id->ICNTL(4)=4;
    /*(*id_mumps).comm_fortran=USE_COMM_WORLD;*/
    /*dmumps_c(id_mumps);*/
    ierr = MPI_Barrier(MPI_COMM_WORLD);
    zmumps_c(id);
    if(verbose==1 && myid==0) printf("dmumps_par->instance_number %d \n",id->instance_number);
#else 
    printf("Please compile poc-solvers library with -DMUMPS \n");
    exit(-1);
#endif
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT UMFPACK SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 else if (strcmp(solver_name,"UMFPACK") == 0){
#ifdef UMFPACK
    slv->parameters = (UMFPACK_STRUC_C *) malloc(sizeof(UMFPACK_STRUC_C));
    id_umf = (UMFPACK_STRUC_C *) slv->parameters;
    Control=id_umf->Control;
    umfpack_zi_defaults (id_umf->Control) ;
/*     Control [UMFPACK_PRL] = 6 ; */
/*     Control [UMFPACK_PRL] = 5 ; */
/*     umfpack_di_report_control (Control) ; */
#else 
    printf("Please compile poc-solvers library with -DUMFPACK \n");
    exit(-1);
#endif
   }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT HIPS SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"HIPS") == 0){
#ifdef HIPS
    slv->parameters = (zHIPS_STRUC_C *) malloc(sizeof(zHIPS_STRUC_C));
    par_hips = (zHIPS_STRUC_C *) slv->parameters;
    status=ZHIPSinitialise_IDs();
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
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    slv->rank=rank;
    idnbr = 100; /* total */
    if (HIPSSOLVERNBZ == 0) {
      ierr = zHIPS_Initialize(idnbr);
      zHIPS_ExitOnError(ierr);
      }
    HIPSSOLVERNBZ++;    
//     id_hips = HIPSSOLVERNBZ; /** id_hips of the linear system **/
    par_hips->id_hips = ZHIPSget_IDs(); /** id_hips of the linear system **/
#else 
    printf("Please compile poc-solvers library with -DHIPS  \n");
    exit(-1);
#endif
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT PASTX SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"PASTIX") == 0){
#ifdef PASTIX
#if HAVE_MPI
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
#endif
    slv->parameters = (PASTIX_STRUC_C *) calloc(1,sizeof(PASTIX_STRUC_C));
    id_pastix = (PASTIX_STRUC_C *) slv->parameters;
#else
    printf("%s (line %d) : please compile poc-solvers library with -DPASTIX \n", __func__,__LINE__);
    exit(-1);
#endif
  } 

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT PETSC SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"PETSC") == 0){
    printf("PETSC pas encore implimente\n");
    exit(-1);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT LAPACK SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"LAPACK") == 0){
    slv->parameters = NULL;
#ifdef LAPACK
#else
    printf("Please compile poc-solvers library with -DLAPACK \n");
    exit(-1);
#endif
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INIT SpDOMESTIC SOLVER
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  else if (strcmp(solver_name,"SpDOMESTIC") == 0){
    slv->parameters = NULL;   
    }

  else {
    printf("Connait pas ce solveur : %s \n",solver_name);
    slv = NULL;
    exit(-1);
    }

  return(slv);
}

#endif

/*=======================================================================*/
/*                 Solver Terminate                                      */
/*************************************************************************/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int freez_solver(solver_t *slv)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  /*mumps_solver mumpsslv;*/
  /*int len,ierr,myid;*/
  int status=0;
#ifdef MUMPS
  ZMUMPS_STRUC_C *id_mumps;
#endif
#ifdef UMFPACK
  UMFPACK_STRUC_C *id_umf;
#endif
#ifdef PASTIX
  PASTIX_STRUC_C *id_pastix;
  complex<double> *Ax=0, *rhs=0;
  int *Ap=0, *Ai=0,n,nnz;
#endif
  /*double *Control;*/
#ifdef HIPS
  zHIPS_STRUC_C *id;
  INTS id_hips;
  FILE *out;
  int ierr;
#endif
  char *solver_name=0;
  void *Symbolic, *Numeric ;
  tripletz *Mat;

  if(slv!=0) {
    solver_name=slv->name;
    Mat = (tripletz *) slv->Mat;
    }

  if (solver_name != NULL) {
//  status=initialisation_terminate();  

/* ---------------- FREE MUMPS SOLVER -------------------- */
#ifdef MUMPS
    if (strcmp(solver_name,"MUMPS") == 0) {
      id_mumps = (ZMUMPS_STRUC_C *) slv->parameters;
      (*id_mumps).job=JOB_END;
      zmumps_c(id_mumps);
      }
#endif

/* ---------------- FREE UMFPACK SOLVER -------------------- */
#ifdef UMFPACK
    if (strcmp(solver_name,"UMFPACK") == 0) {
      id_umf = (UMFPACK_STRUC_C *) slv->parameters;
      Symbolic=id_umf->Symbolic;
      Numeric=id_umf->Numeric;
      umfpack_zi_free_symbolic (&Symbolic) ;
      umfpack_zi_free_numeric (&Numeric) ;
      }
#endif

/* ---------------- FREE LAPACK SOLVER -------------------- */
#ifdef LAPACK
    if (strcmp(solver_name,"LAPACK") == 0) {
      free(Mat->a);
      free(Mat->pivot);
      }
#endif
    
/* ---------------- FREE SpDOMESTIC SOLVER -------------------- */
    if (strcmp(solver_name,"SpDOMESTIC") == 0) {
      cs_ci_sfree (Mat->S) ;
      cs_ci_nfree (Mat->N) ;
      }

/* ---------------- FREE PASTIX SOLVER -------------------- */
#ifdef PASTIX
    if (strcmp(solver_name,"PASTIX") == 0) {
      Ap=Mat->i;
      Ai=Mat->j;
      Ax=Mat->x;
      id_pastix = (PASTIX_STRUC_C *) slv->parameters;
      id_pastix->iparm[IPARM_START_TASK] = API_TASK_CLEAN;
      id_pastix->iparm[IPARM_END_TASK]   = API_TASK_CLEAN;

      z_pastix(&(id_pastix->pastix_data), MPI_COMM_WORLD,
         n, Ap, Ai, Ax,
         id_pastix->perm, id_pastix->invp, rhs, 1, 
         id_pastix->iparm, id_pastix->dparm);
      free(id_pastix->perm);
      free(id_pastix->invp);
      }
#endif

/* ---------------- FREE HIPS SOLVER -------------------- */
#ifdef HIPS
    if (strcmp(solver_name,"HIPS") == 0) {
      /************************************************/
      /* Free zHIPS internal structure for problem id  */
      /************************************************/
      id = (zHIPS_STRUC_C *) slv->parameters;
      id_hips =id->id_hips;
      ierr = zHIPS_Clean(id_hips);
      zHIPS_ExitOnError(ierr);
      ZHIPSfree_IDs(id_hips);
      HIPSSOLVERNBZ--;
      
      /**********************************/
      /* Free zHIPS internal structures  */
      /**********************************/
      if(HIPSSOLVERNBZ==0) {
        ierr = zHIPS_Finalize();
        zHIPS_ExitOnError(ierr);
        }
      }
#endif

//    if(Mat->private_memory==1) {
//       free(Ax);
//       free(Ap);
//       free(Ai);
//       }
    free(slv->parameters);
    free(slv->name);
    /*free(slv);*/
    slv->parameters = NULL;
    slv->name = NULL;
    /*slv = NULL;*/
    }
  else {
    printf(" free already done with this solver \n");
    }
    
//   slv=0;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_terminate(solver_t *solveur)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0, ierr;

  if(solveur==0) return(0);

  switch (solveur->id) {
    case SOLVER_ID_DIAGONAL:     /* no real solver used for diagonal matrices */
      break;

    case SOLVER_ID_DOMESTIC:     /* DOMESTIC */
      break;

    case SOLVER_ID_LAPACK:       /* LAPACK */
      break;

    case SOLVER_ID_SUNPERF:      /* SUNPERF */
      break;

    case SOLVER_ID_UMFPACK:      /* UMFPACK */
      status=free_solver(solveur);
      break;

    case SOLVER_ID_PETSC_GMRES:        /* ITERATIVE SOLVER */
      status=free_solver(solveur);
      break;

    case SOLVER_ID_MUMPS:        /* MUMPS */
      status=free_solver(solveur);
      break;

    case SOLVER_ID_MUMPS_SYM:    /* MUMPS-SYM */
      status=free_solver(solveur);
      break;

    case SOLVER_ID_SpDOMESTIC:    /* SpDOMESTIC */
      status=free_solver(solveur);
      break;

//    case 8:                      /* HYPRE */
//      break;
    
    case SOLVER_ID_HIPS:           /* HIPS */
      status=free_solver(solveur);
      break;  

    case SOLVER_ID_PASTIX:         /* PASTIX */
      status=free_solver(solveur);
      break;  

    default:
      exit(-1);
      break;
    }

  if(status!=0) {
    printf("solver =%d\n",solveur->id);
    printf("solver_terminate status=%d ...\n",status);
    }

  return (status);
}


/*************************************************************************/
/*=======================================================================*/
/*======================FIN TERMINATE ===================================*/
/*=======================================================================*/


