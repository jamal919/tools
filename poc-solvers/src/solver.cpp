
// /**************************************************************************
// 
//   POC-SOLVERS interface, 2006-2012
// 
//   Part of the Unstructured Ocean Grid initiative and SYMPHONIE suite
// 
// Contributors:
// 
//   Florent Lyard      LEGOS/CNRS, Toulouse, France
//   Cyril Nguyen       LA/CNRS,    Toulouse, France
//   Damien Allain      LEGOS/CNRS, Toulouse, France
//   Yoann Le Bars      PhD, LEGOS, Toulouse, France
// 
// ***************************************************************************/
// 
// // #define __DEBUG_HIPS
// 
// #define VERBOSE
// 
// #include <string.h>
// #include <omp.h>
// 
// #include "solvers-interface.h"
// #include "poc-solvers.h"
// 
// #define JOB_INIT -1
// #define JOB_END -2
// #define USE_COMM_WORLD -987654
// #define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match MUMPS documentation */
// 
// /*------------------------------------------------------------------------------
//   temporary patch for compatibility */
// // #define solver solver_t
// 
// /**-----------------------------------------------------------------------------
//   added because of PASTIXinitilisation */
//   int MPI_Init_thread_DONE=0;
//   
//   int HIPSSOLVERNB=0;
//   int HIPSSOLVERNBZ=0;
// 
//   int nbsyslinearz2=0;
// 
// /*=======================================================================*/
// /*                            INFOS                                      */
// /*************************************************************************/
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// void  set_MPI_Init_thread_DONE(int done) {
//   MPI_Init_thread_DONE=done;
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// int ishere(char *solver_name) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   
//   if (strcmp(solver_name,"MUMPS") == 0) {
// #ifdef MUMPS
//     return(1);
// #else
//     return(0);
// #endif
//     }
//   else if (strcmp(solver_name,"UMFPACK") == 0){
// #ifdef UMFPACK
//     return(1);
// #else
//     return(0);
// #endif
//     }
//   else if (strcmp(solver_name,"HIPS") == 0){
// #ifdef HIPS
//     return(1);
// #else
//     return(0);
// #endif
//     }
//   else if (strcmp(solver_name,"LAPACK") == 0){
// #ifdef LAPACK
//     return(1);
// #else
//     return(0);
// #endif
//     }
//   else if (strcmp(solver_name,"HIPS") == 0){
// #ifdef HIPS
//     return(1);
// #else
//     return(0);
// #endif
//     }
//   else if (strcmp(solver_name,"PASTIX") == 0){
// #ifdef PASTIX
//     return(1);
// #else
//     return(0);
// #endif
//     }
//   else if (strcmp(solver_name,"MAPHYS") == 0){
// #ifdef MAPHYS
//     return(1);
// #else
//     return(0);
// #endif
//     }
//   else if (strcmp(solver_name,"SpDOMESTIC") == 0){
//     return(1);
//     }
//   else {
//     return(0);
//     }
//       
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int LinearSystem_identify(char *solver_name)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int id;
//   id = -1;
//   if(strcmp(solver_name, "DOMESTIC")   == 0) 
//     id = SOLVER_ID_DOMESTIC;
//   if(strcmp(solver_name, "LAPACK")     == 0)
//     id = SOLVER_ID_LAPACK;
//   if(strcmp(solver_name, "SUNPERF")    == 0)
//     id = SOLVER_ID_SUNPERF;
//   if(strcmp(solver_name, "UMFPACK")    == 0)
//     id = SOLVER_ID_UMFPACK;
//   if(strcmp(solver_name, "GMRES")      == 0)
//     id = SOLVER_ID_GMRES;
//   if(strcmp(solver_name, "MUMPS")      == 0)
//     id = SOLVER_ID_MUMPS;
//   if(strcmp(solver_name, "MUMPS_SYM")  == 0)
//     id = SOLVER_ID_MUMPS_SYM;
//   if(strcmp(solver_name, "SpDOMESTIC") == 0)
//     id = SOLVER_ID_SpDOMESTIC;
//   if(strcmp(solver_name, "PASTIX") == 0)
//     id = SOLVER_ID_PASTIX;
//   if(strcmp(solver_name, "HIPS")       == 0)
//     id = SOLVER_ID_HIPS;
//   if(strcmp(solver_name, "MAPHYS") == 0)
//     id = SOLVER_ID_MAPHYS;
//   if(strcmp(solver_name, "TEST1")      == 0)
//     id = 99;
//   if(strcmp(solver_name, "TEST2")      == 0)
//     id = 98;
// 
//   return id;
// }
// 
// 
// /*************************************************************************/
// /*======================= FIN INFOS ==================================== */
// /*=======================================================================*/
//   
// 
// /*=======================================================================*/
// /*                            INIT                                       */
// /*************************************************************************/
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// solver *init_solver(char *solver_name,int typsym, int verbose)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// {
//   solver *slv;
//   /*mumps_solver mumpsslv;*/
//   int ierr,myid; 
//  
// #ifdef MUMPS 
//   DMUMPS_STRUC_C *id_mumps;
//   DMUMPS_STRUC_C *id;
//  
// #endif
// #ifdef UMFPACK
//   UMFPACK_STRUC_C *id_umf; 
//   double *Control;
// #endif
// #ifdef HYPRE
//   HYPRE_STRUC_C *id_hypre; 
//   HYPRE_IJMatrix A;
//   HYPRE_ParCSRMatrix parcsr_A;
//   HYPRE_IJVector b;
//   HYPRE_ParVector par_b;
//   HYPRE_IJVector x;
//   HYPRE_ParVector par_x;
//   HYPRE_Solver hy_solver, hy_precond;
// #endif
// #ifdef HIPS
//   INTS *id_hips; 
//   INTS idnbr, i, j;
//   INTS *unknownlist;
//   COEFD *x, *rhsloc;
//   INTS proc_id, n, ln;
//   INTL *ia, nnz;
//   INTS *ja;
//   COEFD *a;
//   INTS domsize, nproc;
//   INTS pbegin, pend;
// #endif
// #ifdef PASTIX
//   PASTIX_STRUC_C *id_pastix; 
// #endif
// #ifdef MAPHYS
//   MAPHYS_STRUC_C *id_maphys; 
//   MAPHYS_IJMatrix A;
//   MAPHYS_ParCSRMatrix parcsr_A;
//   MAPHYS_IJVector b;
//   MAPHYS_ParVector par_b;
//   MAPHYS_IJVector x;
//   MAPHYS_ParVector par_x;
//   MAPHYS_Solver maphys_solver;
// #endif
//   int argc=0;
//   char **argv=0;
//   int             required;           /* MPI thread level required                                 */
//   int             provided;           /* MPI thread level provided                                 */
//   int             rank;
//   /*slv = (solver *) calloc(1,sizeof(solver));*/
//   slv = (solver *) malloc(sizeof(solver));
//   slv->name = strdup(solver_name);
// 
//   /*-------------------------------------------------------------*/
//   /* *********            INIT VARIOUS SOLVERS             *******/
// 
// 
// 
// 
//   /* ---------------- INIT MUMPS SOLVER -------------------- */
//   if (strcmp(solver_name,"MUMPS") == 0) {
// #ifdef MUMPS
//     required = MPI_THREAD_MULTIPLE;
//     provided = -1;
//     if(MPI_Init_thread_DONE==0) {
//       MPI_Init_thread(&argc, &argv, required, &provided);
//       MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//       if (rank == 0) {
//         switch (provided) {
//           case MPI_THREAD_SINGLE:
//             printf("MPI_Init_thread level = MPI_THREAD_SINGLE\n");
//             break;
//           case MPI_THREAD_FUNNELED:
//             printf("MPI_Init_thread level = MPI_THREAD_FUNNELED\n");
//             break;
//           case MPI_THREAD_SERIALIZED:
//             printf("MPI_Init_thread level = MPI_THREAD_SERIALIZED\n");
//             break;
//           case MPI_THREAD_MULTIPLE:
//             printf("MPI_Init_thread level = MPI_THREAD_MULTIPLE\n");
//             break;
//           default:
//             printf("MPI_Init_thread level = ???\n");
//             break;
//           }
//         }
//       MPI_Init_thread_DONE=1;
//       }
// //     else {
// //       printf("%s : %s initialisation ; MPI already initialised \n",__FUNCTION__, solver_name);
// //       }
//     ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
// #if 0
//     ierr = MPI_Barrier(MPI_COMM_WORLD);
//     if(verbose==1) printf("init MUMPS solver (cpu=%d), start \n",myid);
// #endif
//     /*slv->parameters = (DMUMPS_STRUC_C *) malloc(sizeof(DMUMPS_STRUC_C));*/
//     slv->parameters = (DMUMPS_STRUC_C *) calloc(1,sizeof(DMUMPS_STRUC_C));
// //    id = (DMUMPS_STRUC_C *) calloc(1,sizeof(DMUMPS_STRUC_C*));
//     id =  slv->parameters;
//     /*id_mumps = slv->parameters;*/
//     /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
//     /*(*id_mumps).job=JOB_INIT; (*id_mumps).par=1; (*id_mumps).sym=typsym;*/
//     id->job=JOB_INIT; id->par=1; id->sym=0;
//     id->comm_fortran=USE_COMM_WORLD;
//     /*(*id_mumps).comm_fortran=USE_COMM_WORLD;*/
//     /*dmumps_c(id_mumps);*/
// #if 0
//     ierr = MPI_Barrier(MPI_COMM_WORLD);
// #endif
//     dmumps_c(id);
//     if(verbose==1 && myid==0) printf("dmumps_par->instance_number %d \n",id->instance_number);
//     /*slv->parameters = id_mumps;*/
//     /*slv->parameters = id;*/
//     /* ---------------- FIN INIT MUMPS SOLVER -------------------- */
//     if(verbose==1) printf("init MUMPS solver (cpu=%d), finish \n",myid);
// #else 
//     printf("Please compile poc-solvers library with -DMUMPS \n");
//     exit(-1);
// #endif
//     }
// 
// 
// 
//   /* --------------- INIT UMFPACK SOLVER ------------------ */
//   else if (strcmp(solver_name,"UMFPACK") == 0){
// #ifdef UMFPACK
//     slv->parameters = (UMFPACK_STRUC_C *) malloc(sizeof(UMFPACK_STRUC_C));
//     id_umf = slv->parameters;
//     Control=id_umf->Control;
//     umfpack_di_defaults (id_umf->Control);
// /*     Control [UMFPACK_PRL] = 6 ; */
// /*     Control [UMFPACK_PRL] = 5 ; */
// /*     umfpack_di_report_control (Control) ; */
// /*    printf("Fin init umf\n");*/
//     /* --------------- FIN INIT UMFPACK SOLVER ------------------ */
// #else 
//     printf("Please compile poc-solvers library with -DUMFPACK \n");
//     exit(-1);
// #endif
//   }  
// 
//   /* --------------- INIT HYPRE SOLVER ------------------ */
//   else if (strcmp(solver_name,"HYPRE") == 0){
// #ifdef HYPRE
//     slv->parameters = (HYPRE_STRUC_C *) malloc(sizeof(HYPRE_STRUC_C));
//     id_hypre = slv->parameters;
//   /* --------------- FIN INIT HYPRE SOLVER ------------------ */
// #else 
//     printf("Please compile poc-solvers library with -DHYPRE \n");
//     exit(-1);
// #endif
//   }
// 
//   /* --------------- INIT HIPS SOLVER ------------------ */
//   else if (strcmp(solver_name,"HIPS") == 0){
// #ifdef HIPS
//     slv->parameters = (dHIPS_STRUC_C *) malloc(sizeof(dHIPS_STRUC_C));
// //     id_hips = slv->parameters;
//     required = MPI_THREAD_MULTIPLE;
//     provided = -1;
//     if(MPI_Init_thread_DONE==0) {
//       MPI_Init_thread(&argc, &argv, required, &provided);
//       MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//       if (rank == 0) {
//         switch (provided) {
//           case MPI_THREAD_SINGLE:
//             printf("MPI_Init_thread level (HIPS) = MPI_THREAD_SINGLE\n");
//             break;
//           case MPI_THREAD_FUNNELED:
//             printf("MPI_Init_thread level (HIPS) = MPI_THREAD_FUNNELED\n");
//             break;
//           case MPI_THREAD_SERIALIZED:
//             printf("MPI_Init_thread level (HIPS) = MPI_THREAD_SERIALIZED\n");
//             break;
//           case MPI_THREAD_MULTIPLE:
//             printf("MPI_Init_thread level (HIPS) = MPI_THREAD_MULTIPLE\n");
//             break;
//           default:
//             printf("MPI_Init_thread level (HIPS) = ???\n");
//             break;
//           }
//         }
//       MPI_Init_thread_DONE=1;
//       }
//     /* --------------- FIN INIT HYPRE SOLVER ------------------ */
// #else
//     printf("Please compile poc-solvers library with -DHIPS \n");
//     exit(-1);
// #endif
//   }
// 
//   
//   /* --------------- INIT PASTX SOLVER ------------------ */
//   else if (strcmp(solver_name,"PASTIX") == 0){
// #ifdef PASTIX
//     slv->parameters = (PASTIX_STRUC_C *) calloc(1,sizeof(PASTIX_STRUC_C));
//     /*id_pastix = (PASTIX_STRUC_C *) calloc(1,sizeof(PASTIX_STRUC_C));*/
//     
//     id_pastix = slv->parameters;
// #ifdef HAVE_MPI
//     if(verbose) printf("init_solver: PASTIX parallel mpi version \n");
//     required = MPI_THREAD_MULTIPLE;
//     provided = -1;
//     if(MPI_Init_thread_DONE==0) {
//       MPI_Init_thread(&argc, &argv, required, &provided);
//       MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//       if (rank == 0) {
//         switch (provided) {
//           case MPI_THREAD_SINGLE:
//             printf("MPI_Init_thread level = MPI_THREAD_SINGLE\n");
//             break;
//           case MPI_THREAD_FUNNELED:
//             printf("MPI_Init_thread level = MPI_THREAD_FUNNELED\n");
//             break;
//           case MPI_THREAD_SERIALIZED:
//             printf("MPI_Init_thread level = MPI_THREAD_SERIALIZED\n");
//             break;
//           case MPI_THREAD_MULTIPLE:
//             printf("MPI_Init_thread level = MPI_THREAD_MULTIPLE\n");
//             break;
//           default:
//             printf("MPI_Init_thread level = ???\n");
//             break;
//           }
//         }
//       MPI_Init_thread_DONE=1;
//       }
// //     else {
// //       printf("%s : %s initialisation ; MPI already initialised \n",__FUNCTION__, solver_name);
// //       }
// #endif
//     /* --------------- FIN INIT PASTIX SOLVER ------------------ */
// #else 
//     printf("Please compile poc-solvers library with -DPASTIX \n");
//     exit(-1);
// #endif
//     } 
//   
//   /* --------------- INIT PETSC SOLVER ------------------ */
//   else if (strcmp(solver_name,"PETSC") == 0){
//     printf("PETSC pas encore implimente\n");
//     exit(-1);
//     }
// 
//   /* --------------- INIT LAPACK SOLVER ------------------ */
//   else if (strcmp(solver_name,"LAPACK") == 0){
//     slv->parameters = NULL;
//     /* --------------- FIN INIT LAPACK SOLVER ------------------ */
// #ifdef LAPACK
// #else
//     printf("Please compile poc-solvers library with -DLAPACK \n");
//     exit(-1);
// #endif
//     }
//  
//   /* --------------- INIT SpDOMESTIC SOLVER ------------------ */
//   else if (strcmp(solver_name,"SpDOMESTIC") == 0){
//     slv->parameters = NULL;   
//     /* --------------- FIN INIT LAPACK SOLVER ------------------ */
//     /*exit(-1);*/
//     }
// 
//   else {
//     printf("Connait pas ce solveur : %s \n",solver_name);
//     slv = NULL;
//     exit(-1);
//     }
// 
//   return(slv);
// }
// /*************************************************************************/
// /*======================= FIN INIT  ==================================== */
// /*=======================================================================*/
// 
// 
// /*************************************************************************/
// /*************************************************************************/
// /*=============== Double precision routines =============================*/
// /*************************************************************************/
// /*************************************************************************/
// 
// /*=======================================================================*/
// /*=======================================================================*/
// /*                               FACTORIZE                               */
// /*************************************************************************/
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int factorize(solver *slv,triplet *Mat, int verbose) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   char *solver_name;
//   int rtn,base,format;
//   double default_precision=1.e-08;
//   double cputime1,cputime2;
//   int proc_id,ierr;  
//   bool debug=false; 
// 
//   solver_name=slv->name;
//   
//   slv->Mat=Mat;
//   
// #ifdef HAVE_MPI
//    if(MPI_Init_thread_DONE==1) {
//      ierr = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
//      cputime1 = MPI_Wtime();
//      }
// #endif
//   
//   if (strcmp(solver_name,"MUMPS") == 0) {
//     base=1;
//     format=COO;
//     rtn = convert(Mat, format, base, verbose);
//     if(rtn!=0) return(-1);
//     rtn = factorize_mumps(slv, Mat, verbose, false);
//     }
//   else if  (strcmp(solver_name,"UMFPACK") == 0) {
//     base=0;
//     format=CSC;
//     rtn = convert(Mat, format, base, verbose);
//     if(rtn!=0) return(-1);
//     rtn=factorize_umfpack(slv, Mat, verbose);
//     }
//   else if  (strcmp(solver_name,"HYPRE") == 0) {
//     base=1;
//     format=COO;
//     rtn = convert(Mat, format, base, verbose);
//     if(rtn!=0) return(-1);
//     rtn=factorize_hypre(slv, Mat);
//     }
//   else if  (strcmp(solver_name,"HIPS") == 0) {
//     base=1;
//     format=COO;
//     rtn = convert(Mat, format, base, verbose);
//     if(rtn!=0) return(-1);
//     if(MPI_Init_thread_DONE==0) {
//       }
// #ifdef __DEBUG_HIPS
//     debug=true;
// #endif
//     rtn=factorize_hips(slv, Mat, verbose, debug);
//     }
//   else if  (strcmp(solver_name,"LAPACK") == 0) {
//     base=1;
//     format=BND; 
//     rtn = convert(Mat, format, base, verbose);
//     if(rtn!=0) return(-1);
//     rtn=factorize_lapack(slv, Mat);
//     }
//   else if  (strcmp(solver_name,"SpDOMESTIC") == 0) {
//     base=0;
//     format=CSC;
//     rtn = convert(Mat, format, base, verbose);
//     if(rtn!=0) return(-1);
//     rtn=factorize_spdomestic(slv, Mat);
//     }
//   else if  (strcmp(solver_name,"PASTIX") == 0) {
// //#ifdef  HAVE_MPI_HERE
// #ifdef  HAVE_MPI
//     // Attention en mode parallel les matrix arrivent en FORMAT COO
//     // avec des indices globaux...
//     // la redistibution, conversion, ...
//     // Les jobs sont faits dans factorize_paxtix_parallel
//     base=1;
// //     format=COO_COL;
//     format=COO_COL_ARRANGED;
//     rtn = convert(Mat, format, base, verbose);
//     if(debug) printf("call factorize_pastix_parallel .......;\n");
//     rtn=factorize_pastix_parallel(slv, Mat, verbose);
// #else    
//     base=0;
//     format=CSC;
//     rtn = convert(Mat, format, base, verbose);
//     if(rtn!=0) return(-1);
//     rtn=factorize_pastix_sequential(slv, Mat);
// #endif
//     }
//   else {
//     printf("factorize, unknown solver : %s \n",solver_name);
//     rtn = -1;
//     }
//     
// #ifdef  HAVE_MPI
//    if(MPI_Init_thread_DONE==1) {
//      cputime2 = MPI_Wtime();
//      if (verbose==1 && proc_id == 0) {
//        printf("--------------------------------------------\n");
//        printf("%s: %s factorization ellapsed time= %e\n",__FUNCTION__, solver_name, cputime2-cputime1);
//        printf("--------------------------------------------\n"); 
//        }
//      }
// #endif
// 
//    return(rtn);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int factorize_mumps(solver *slv,triplet *Mat, int verbose, bool debug) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// #ifdef MUMPS
//   DMUMPS_STRUC_C *id;
//   int myid, ierr, nglob=0,count,rank,size;
//   char *solver_name;
//     
// #ifdef HAVE_MPI
//   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
// #else
//   ierr = 0;
// #endif
//   
//   solver_name=slv->name;
//   if (strcmp(solver_name,"MUMPS") == 0) {
//     if(debug) printf("cpu=%d, factorize MUMPS ...\n",myid);
//     id = slv->parameters;
//     (*id).job=4;
//     
//     /* get global row number */
//     count=1;
// #ifdef HAVE_MPI
//     ierr = MPI_Allreduce (  &(Mat->nrow), &nglob, count, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
//     (*id).n = nglob;
// #else
//     (*id).n = Mat->nrow;
//     ierr = 0;
// #endif
//   
// /**-----------------------------------------------------------------------------
//     triplets format */
//     (*id).ICNTL(18)=3;
//     (*id).nz_loc = Mat->nnz;
//     (*id).irn_loc= Mat->i;  
//     (*id).jcn_loc= Mat->j; 
//     (*id).a_loc =  Mat->x;  
// 
// /**-----------------------------------------------------------------------------
//     output streams */
//     (*id).ICNTL(1)=-1; 
//     /*(*id).ICNTL(1)=6; */
//     (*id).ICNTL(2)=-1; 
//     /*(*id).ICNTL(2)=0; */
//     (*id).ICNTL(3)=-1; 
//     /*(*id).ICNTL(3)=6; */
// /**-----------------------------------------------------------------------------
//     level of printing */
//     (*id).ICNTL(4)=4;
// /**-----------------------------------------------------------------------------
//     percentage increase in estimated workspace */
//     (*id).ICNTL(14)=100;
// /**-----------------------------------------------------------------------------
//     dense/sparse RHS */
//     (*id).ICNTL(20)=0;
// /**-----------------------------------------------------------------------------
//     centralized/distributed solution vector */
//     (*id).ICNTL(21)=0;
//     }
//   else {
//     }
//   dmumps_c(id);
//   ierr=(*id).info[0];
//   if(ierr==-1)  printf("cpu %d: DMUMPS factorisation failed, status=%d due to cpu %d\n",rank, ierr, id->info[1]);
//   if(ierr==-9)  printf("cpu %d: DMUMPS factorisation failed, status=%d due to memory deficience %d\n",rank, ierr, id->info[1]);
//   if(ierr==-10) printf("cpu %d: DMUMPS factorisation failed, status=%d due to singular matrix\n",rank, ierr);
//   if(ierr==-13) printf("cpu %d: DMUMPS factorisation failed, status=%d due to memory deficience %d\n",rank, ierr, id->info[2]);
//   if(ierr==0) {
//     if(verbose==1) {
//       printf("%s : MEMORY USE, %d used by cpu, %d total Mo, status=%d\n",__FUNCTION__, (*id).infog[20],(*id).infog[21], ierr);
//       }
//     if(verbose==1 && rank==0) {
//       printf("MEMORY SCALING DIAG %d %d\n",size,(*id).infog[21]);
//       }
//     }
//   return(ierr);
//   
// #endif
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int factorizez_mumps(solver *slv, tripletz *Mat, int verbose, bool debug) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// #ifdef MUMPS
//   ZMUMPS_STRUC_C *id;
//   int myid, ierr, nglob=0, count, rank, size;
//   char *solver_name;
//   
// #ifdef HAVE_MPI
//   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
// #else
//   ierr = 0;
// #endif
//   
//   solver_name=slv->name;
//   if (strcmp(solver_name,"MUMPS") == 0) {
//     if(debug) printf("cpu=%d, factorizez MUMPS ...\n",myid);
//     id = slv->parameters;
//     (*id).job=4;
//     
//     /* get global row number */
//     count=1;
// #ifdef HAVE_MPI
//     ierr = MPI_Allreduce (  &(Mat->nrow), &nglob, count, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
//     (*id).n = nglob;
// #else
//     (*id).n = Mat->nrow;
//     ierr = 0;
// #endif
//   
//     (*id).ICNTL(18)=3;
// /*------------------------------------------------------------------------------
//     5.1.1 and later */
//     (*id).nnz_loc = Mat->nnz;
// //     (*id).nz_loc = Mat->nnz;
//     (*id).irn_loc= Mat->i;  
//     (*id).jcn_loc= Mat->j; 
//     (*id).a_loc =  Mat->x;  
// 
//     (*id).ICNTL(1)=-1; 
//     /*(*id).ICNTL(1)=6; */
//     (*id).ICNTL(2)=-1; 
//     /*(*id).ICNTL(2)=0; */
//     (*id).ICNTL(3)=-1; 
//     /*(*id).ICNTL(3)=6; */
//     (*id).ICNTL(4)=4;
//     (*id).ICNTL(14)=100;
//     /*...... Centtralized RHS, SOl ... ***/
//     (*id).ICNTL(20)=0;
//     (*id).ICNTL(21)=0;
//     }
//   else {
//     }
//   zmumps_c(id);
//   
//   ierr=(*id).info[0];
//   if(ierr==-1)  printf("cpu %d: ZMUMPS factorisation failed, status=%d due to cpu %d\n",rank, ierr, id->info[1]);
//   if(ierr==-9)  printf("cpu %d: ZMUMPS factorisation failed, status=%d due to memory deficiency %d\n",rank, ierr, id->info[1]);
//   if(ierr==-10) printf("cpu %d: ZMUMPS factorisation failed, status=%d due to singular matrix\n",rank, ierr);
//   if(ierr==-13) printf("cpu %d: DMUMPS factorisation failed, status=%d due to memory deficience %d\n",rank, ierr, id->info[2]);
//   if(ierr==0) {
//     if(verbose==1) {
//       printf("%s, cpu=%d: MEMORY USE, %d used by cpu, %d total Mo, status=%d\n",__FUNCTION__,rank, (*id).infog[20],(*id).infog[21], ierr);
//       }
//     if(verbose==1 && rank==0) {
//       printf("MEMORY SCALING DIAG %d %d\n",size,(*id).infog[21]);
//       }
//     }
//   return(ierr);
//   
// #endif
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int factorize_mumps_V0(solver *slv,triplet *Mat) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// #ifdef MUMPS
// #ifdef HAVE_MPI
//   DMUMPS_STRUC_C *id;
//   int myid, ierr, nglob,count;
//   char *solver_name;
//   
//   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
// 
//   solver_name=slv->name;
//   if (strcmp(solver_name,"MUMPS") == 0) {
//     id = slv->parameters;
//     /*printf(" MUMPS parameter sym =%d myid=%d \n",(*id).sym, myid);*/
//     /*printf("dmumps_par->instance_number %d \n",id->instance_number);*/
//     /*printf("nnz=%d, nzmax=%d  \n",Mat->nnz,Mat->nzmax);*/
//     (*id).job=4;
//     
//     /* get global row nuber */
//     count=1;
//     ierr = MPI_Allreduce (  &(Mat->nrow), &nglob, count, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
//   
//     (*id).n = nglob; 
//     (*id).ICNTL(18)=3;
//     (*id).nz_loc = Mat->nnz;
//     (*id).irn_loc= Mat->i;  
//     (*id).jcn_loc= Mat->j; 
//     (*id).a_loc =  Mat->x;  
// 
//     (*id).ICNTL(1)=-1; 
//     /*(*id).ICNTL(1)=6; */
//     (*id).ICNTL(2)=-1; 
//     /*(*id).ICNTL(2)=0; */
//     (*id).ICNTL(3)=-1; 
//     /*(*id).ICNTL(3)=6; */
//     (*id).ICNTL(4)=4;
// 
//     /*...... Centtralized RHS, SOl ... ***/
//     (*id).ICNTL(20)=0;
//     (*id).ICNTL(21)=0;
//     }
//   else {
//     }
//   dmumps_c(id);
//   printf(" MEM = %d, %d Mo\n",(*id).infog[20],(*id).infog[21]);
//   printf(" ierr=%d \n",ierr);
//   return(ierr);
//   
// #endif
// #endif
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int factorize_umfpack(solver *slv, triplet *Mat, int verbose) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   int status;
//   
// #ifdef UMFPACK
//   UMFPACK_STRUC_C *id;
//   char *solver_name;
//   double *Info, *Control, *Ax;
//   int *Ap, *Ai,n,nnz;
//   void *Symbolic, *Numeric ;
//   
//   solver_name=slv->name;
//   if (strcmp(solver_name,"UMFPACK") == 0) {
//     if(verbose) printf("Factorize UMFPACK...\n");
//     id = slv->parameters;
//     slv->Mat = Mat;
//     Symbolic=id->Symbolic;
//     Numeric=id->Numeric;
//     Info=id->Info;
//     Control=id->Control;
//     Ap=Mat->i;
//     Ai=Mat->j;
//     Ax=Mat->x;
//     n=Mat->nrow;
//     nnz=Mat->nnz;
//     status = umfpack_di_symbolic (n, n, Ap, Ai, Ax, &Symbolic, Control, Info) ; 
//     if (status != 0) {
//       printf("Factorisation status=%d\n",status);
//       umfpack_di_report_info (Control, Info) ;
//       umfpack_di_report_status (Control, status) ;
//         /*error ("umfpack_di_symbolic failed") ;*/
//       }
//     status = umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, Control, Info) ;
//     id->Symbolic=Symbolic;
//     id->Numeric=Numeric; 
//     if (status != 0) {
//       printf("Factorisation status=%d\n",status);
//       umfpack_di_report_info (Control, Info) ;
//       umfpack_di_report_status (Control, status) ;
//         /*error ("umfpack_di_numeric failed") ;*/
//       }
//     return(status);
//   }
// #else
//   status=-1;
//   printf("Please compile poc-solvers library with -DUMFPACK \n");
//   return(status);
// #endif
//   return(status); 
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int factorize_hypre(solver *slv,triplet *Mat) {
//     
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   int status;
// #ifdef HYPRE
// #ifdef HAVE_MPI
//   HYPRE_STRUC_C *id_hypre; 
//   HYPRE_IJMatrix A;
//   HYPRE_ParCSRMatrix parcsr_A;
//   HYPRE_IJVector b;
//   HYPRE_ParVector par_b;
//   HYPRE_IJVector x;
//   HYPRE_ParVector par_x;
//   HYPRE_Solver hy_solver, hy_precond;
//   HYPRE_STRUC_C *id;
//   int ilower,iupper,bcl1,bcl2,ptr;
//   char *solver_name;
//   char filename[50];
//   int myid, ierr, ncolmax, ncol, oldlig;
//   FILE *out;
//   int lig, col;
//   double val;
//  
// 
//   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//   sprintf(filename,"solver_fatorize_hyre.1.%d.txt",myid);
//   solver_name=slv->name;
//   id = slv->parameters;
//   printf("%d, n=%d,nnz=%d \n",myid,Mat->nrow,Mat->nnz);
//   printf("%d, ligdeb=%d, ligfin=%d \n",myid,Mat->i[0],Mat->i[Mat->nrow-1]);
//   out=fopen(filename,"w");
//   ilower=Mat->nrow;
//   iupper=0;
//   ncolmax=0;
//   ncol=1;
//   oldlig=Mat->i[0];
//   bcl1=0;
//   fprintf(out, "%d %d %lf \n",Mat->i[bcl1],Mat->j[bcl1],Mat->x[bcl1]); 
//   for(bcl1 = 1; bcl1 < Mat->nnz; bcl1++) {
//     fprintf(out, "%d %d %lf \n",Mat->i[bcl1],Mat->j[bcl1],Mat->x[bcl1]); 
//     if (ilower > Mat->i[bcl1]) ilower=Mat->i[bcl1];
//     if (iupper < Mat->i[bcl1]) iupper=Mat->i[bcl1];
//     if (Mat->i[bcl1] == oldlig) {
//       ncol=ncol+1;
//     }
//     else {
//       oldlig=Mat->i[bcl1];
//       if (ncol > ncolmax) ncolmax=ncol;
//       ncol=1;
//     }
//   }
//   fclose(out);
//   printf("%d, ilower=%d, iupper=%d ncolmax=%d \n",myid,ilower,iupper,ncolmax);
//   /* Create the matrix.
//       Note that this is a square matrix, so we indicate the row partition
//       size twice (since number of rows = number of cols) */
//    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
//    status = MPI_Barrier(MPI_COMM_WORLD);
//    printf(" IJMatrixCreate OK.....\n");
//    /* Choose a parallel csr format storage (see the User's Manual) */
//    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
//    /* Initialize before setting coefficients */
//    HYPRE_IJMatrixInitialize(A);
//    status = MPI_Barrier(MPI_COMM_WORLD);
//    printf(" HYPRE_IJMatrixInitialize OK.....\n");
// 
//   sprintf(filename,"solver_fatorize_hyre.2.%d.txt",myid);
//   out=fopen(filename,"w");   
//   for(bcl1 = 0; bcl1 < Mat->nnz; bcl1++) {
//     lig =  Mat->i[bcl1]-1;
//     col =  Mat->j[bcl1]-1;
//     val =  Mat->x[bcl1];
//     ncol=1;
//     fprintf(out, "%d %d %lf \n",lig,col,val); 
//     /*HYPRE_IJMatrixSetValues(A, 1, &ncol, &lig, &col, &val);*/
// 
//   }
//   fclose(out);
//    status = MPI_Barrier(MPI_COMM_WORLD);
//    exit(-23);
// 
//    /* Assemble after setting the coefficients */
//    HYPRE_IJMatrixAssemble(A);
//    /* Get the parcsr matrix object to use */
//    HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);
//    id->A=A;
//    id->parcsr_A=parcsr_A;
//    /*---pour test */
//    HYPRE_IJMatrixPrint(id->A, "IJ.out.A");
//    status = MPI_Barrier(MPI_COMM_WORLD);
//    exit(-23);
// #endif
// #endif 
//   return(status); 
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int factorize_hips(solver *slv,triplet *Mat, int verbose, bool debug)
//   
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status,ierr;
// #ifdef HIPS
// #ifdef HAVE_MPI
//   dHIPS_STRUC_C *id;
//   INTS  sym_pattern, sym_matrix;
//   INTS id_hips, idnbr, i, j;
//   INTS *unknownlist, *hipsnodelist;
//   COEFD *x, *rhsloc;
//   INTS proc_id, n, ln;
//   INTL *ia, nnz;
//   INTS *ja;
//   COEFD *a;
//   INTS domsize, nproc;
//   INTS pbegin, pend;
//   INTS *mapptr,  *mapptr2,*mapp,*mapp2;
//   INTS p;
// 
//   int nglob, nnzglob, tmp, bcl, idi, info;
//   FILE *out;
//   char filename[50];
//   double *x1,*x2, tol;
//   
//   MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
//   MPI_Comm_size(MPI_COMM_WORLD, &nproc);
//   
//   if (debug && proc_id == 0) {printf("cpu=%d, Entree %s \n",proc_id,__FUNCTION__);}
// 
//   id = slv->parameters; 
//   
//   /* nglob */
//   tmp=Mat->nrow; 
//   
//   if (debug) printf("cpu=%d,  %s nloc=%d\n",proc_id,__FUNCTION__,tmp); 
//   
//   MPI_Allreduce( &tmp, &nglob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
//   
//   if (debug) printf(" %s : cpu=%d, id_hips=%d nglob=%d\n",__FUNCTION__, proc_id, HIPSSOLVERNB, nglob);
// 
//   idnbr = 100; /* total */
//   if (HIPSSOLVERNB == 0) {
// //     printf(" factorize_hips : dHIPS_Initialize \n");
//     ierr = dHIPS_Initialize(idnbr);
// //     printf(" factorize_hips : dHIPS_Initialize ierr=%d\n",ierr);
//     dHIPS_ExitOnError(ierr);
//     }
//   id_hips = HIPSSOLVERNB; /** id_hips of the linear system **/
//   HIPSSOLVERNB++;
//   
//   /** parameter domsize is an argument of testHIPS.ex **/
// 
//   if (debug && proc_id == 0) {printf("cpu=%d,  %s solver numero=%d\n",proc_id,__FUNCTION__,id_hips);}
// 
//   domsize = nglob/nproc/2;
// 
//   if (debug && proc_id == 0) {printf("cpu=%d,  %s domsize=%d\n",proc_id,__FUNCTION__,domsize);}
// 
//   dHIPS_SetOptionINT(id_hips, dHIPS_PARTITION_TYPE, 0); 
//   dHIPS_SetOptionINT(id_hips, dHIPS_DOMSIZE, domsize);
// 
// /*------------------------------------------------------------------------------
//   OPTIONS */
// //   ierr= dHIPS_SetDefaultOptions (id_hips, dHIPS_HYBRID);
//   ierr = dHIPS_SetDefaultOptions (id_hips, dHIPS_ITERATIVE); 
//   tol=1e-12;
//   dHIPS_SetOptionREAL(id_hips,dHIPS_DROPTOL0,tol);
//   dHIPS_SetOptionREAL(id_hips,dHIPS_DROPTOL1,tol);
//   
//   if(verbose==1)
//     info=4;
//   else
//     info=0;
//   dHIPS_SetOptionINT (id_hips, dHIPS_VERBOSE, info);
//   
// //    dHIPS_SetOptionREAL(id_hips, INTS number, REAL value);
//   sym_matrix=0;
//   dHIPS_SetOptionINT (id_hips, dHIPS_SYMMETRIC, sym_matrix);
//   /** C : numbering starts from 0 **/
//   /* but here starts from 1 **/
//   dHIPS_SetOptionINT(id_hips, dHIPS_FORTRAN_NUMBERING, 0);
//   /***************************************************/
//   /*            ENTER THE GRAPH                      */
//   /***************************************************/
//   ierr = dHIPS_GraphBegin(id_hips,nglob, Mat->nnz);
//   dHIPS_ExitOnError(ierr);
//   
//   for (bcl=0;bcl<Mat->nnz;bcl++) {
//     ierr = dHIPS_GraphEdge(id_hips, Mat->i[bcl]-1, Mat->j[bcl]-1);
//     }
//   ierr = dHIPS_GraphEnd(id_hips);
//   dHIPS_ExitOnError(ierr);
//   
//   status = MPI_Barrier(MPI_COMM_WORLD);
//   
//   if (debug && proc_id == 0) {printf("cpu=%d,  %s GraphEnd OK...\n",proc_id,__FUNCTION__);}
// 
//    /***************************************************/
//   /*            ENTER A USER PARTITION               */
//   /***************************************************/
//   /** Only the master processor (0 by default) needs to enter the partition **/
//   /* collect information on proc 0 */
//   mapptr = (INTS *)malloc(sizeof(INTS)*(nproc+1));
//   mapptr2= (INTS *)malloc(sizeof(INTS)*(nproc+1));
//   mapp = (INTS *)malloc(sizeof(INTS)*nglob);
//   mapp2= (INTS *)malloc(sizeof(INTS)*nglob);
// 
//   for (bcl=0;bcl<nproc;bcl++)   mapptr2[bcl]=0;
//   for (bcl=0;bcl<nglob;bcl++)   mapp2[bcl]=0;
//   for (bcl=0;bcl<nproc;bcl++)   mapptr[bcl]=0;
//   for (bcl=0;bcl<nglob;bcl++)   mapp[bcl]=0; 
//   for (bcl=proc_id+1;bcl<nproc+1;bcl++) mapptr2[bcl]=Mat->nrow;
// 
//   MPI_Allreduce( mapptr2, mapptr, nproc+1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
//   
// // #ifdef __DEBUG_HIPS
// //       for (bcl=0;bcl<nproc+1;bcl++) printf("...%d, ",mapptr[bcl]);
// //       printf("\n");
// // #endif
//  
//   idi=0;
//   mapp2[idi+mapptr[proc_id]] =  Mat->i[0]-1;
//   for (bcl=1; bcl < Mat->nnz; bcl++) {
//     if ( mapp2[idi+mapptr[proc_id]] != (Mat->i[bcl]-1) ) {
//       idi++;
//       /*printf("bcl=%d, idi=%d, [idi+mapptr[proc_id]=%d nglob=%d \n",bcl,idi,idi+mapptr[proc_id], nglob);*/
//       mapp2[idi+mapptr[proc_id]] = Mat->i[bcl]-1;    
//       }
//     }
//     
// //       sprintf(filename,"solver_fatorize_hips.1.%d.txt",proc_id);
// //       out=fopen(filename,"w");
// //       for (bcl=0;bcl<nglob;bcl++)   fprintf(out, "%d \n",mapp2[bcl]);
// //       fclose(out);
// 
//   MPI_Allreduce( mapp2, mapp, nglob, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
//       
// //       sprintf(filename,"solver_fatorize_hips.2.%d.txt",proc_id);
// //       out=fopen(filename,"w");
// //       for (bcl=0;bcl<nglob;bcl++)   fprintf(out, "%d \n",mapp[bcl]);
// //       fclose(out);
//       
//   /* Set the unknownlist */
//   unknownlist= (INTS *)malloc(sizeof(INTS)*(Mat->nrow));
//   unknownlist[0]=Mat->i[0]-1;
//   idi=0;
//   for (bcl=1; bcl < Mat->nnz; bcl++) {
//     if ( unknownlist[idi] != (Mat->i[bcl]-1) ) {
//       idi++;
//       unknownlist[idi] = Mat->i[bcl]-1;    
//       }
//     }
//       
// //       sprintf(filename,"solver_fatorize_hips.unknonwlist.%d.txt",proc_id);
// //       out=fopen(filename,"w");
// //       for (bcl=0; bcl < idi; bcl++)  fprintf(out, "%d \n",unknownlist[bcl]);
// //       fclose(out);
// 
//    /* Enter Partition */
//   if (proc_id ==0) {
//     ierr = dHIPS_SetPartition(id_hips, nproc, mapptr, mapp);
//     dHIPS_ExitOnError(ierr);
//     }
//   
//   status = MPI_Barrier(MPI_COMM_WORLD);
// //       printf("%d,  factorize_hips SetPartition OK...\n",proc_id);
//   
//   free(mapptr2);
//   free(mapptr);
//   free(mapp2);
//   free(mapp);
//       
//   /***************************************************/
//   /*            GET THE LOCAL UNKNOWN LIST           */
//   /***************************************************/
// //   ierr = dHIPS_GetLocalUnknownNbr(id_hips, &ln);
// //   dHIPS_ExitOnError(ierr);
// //   printf("%d, factorize_hips nombre de noeuds interne %d\n",proc_id,ln);
// //   hipsnodelist = (INTS *)malloc(sizeof(INTS)*ln);
// //   ierr = dHIPS_GetLocalUnknownList(id_hips, hipsnodelist);
// //   dHIPS_ExitOnError(ierr);
//   
//   /***************************************************/
//   /*          ENTER THE MATRIX COEFFICIENT           */
//   /***************************************************/
//   /** The processor enter the rows pbegin to pend-1 (naive partition) **/
//   /** In the dHIPS_ASSEMBLY_FOOL mode any processor can enter any coefficient :
//       nevertheless it is better to enter as much as possible coefficients
//       that are in the local matrix of the processor **/
//   
// //   sprintf(filename,"solver_fatorize_hips.Matrice.%d.txt",proc_id);
// //   out=fopen(filename,"w");
// //   for(bcl=0;bcl<Mat->nnz;bcl++) {
// //     fprintf(out, "%d : %d,%d,%lf  \n",bcl,Mat->i[bcl], Mat->j[bcl], Mat->x[bcl]);
// //     }
// //   fclose(out);
// 
//   ierr = dHIPS_AssemblyBegin(id_hips, Mat->nnz, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_FOOL, sym_matrix);
//   dHIPS_ExitOnError(ierr);
//   
//   if (debug && proc_id == 0) {printf("cpu %d,  %s AssemblyBegin OK...\n",proc_id,__FUNCTION__);}
// 
//   for(bcl=0;bcl<Mat->nnz;bcl++) {
//     ierr = dHIPS_AssemblySetValue(id_hips, Mat->i[bcl]-1, Mat->j[bcl]-1, Mat->x[bcl]);
//     dHIPS_ExitOnError(ierr);
//     }
//     
//   ierr = dHIPS_AssemblyEnd(id_hips);
//   dHIPS_ExitOnError(ierr);
// 
// //   sprintf(filename,"solver_fatorize_hips.4.%d.txt",proc_id);
// //   out=fopen(filename,"w");
// //   fprintf(out, "%d  %d \n",idi, Mat->nrow);
// //   for (bcl=0; bcl < Mat->nrow; bcl++) fprintf(out, "%d => %d \n",bcl,unknownlist[bcl]);
// //   fclose(out);
//     
//   id->id_hips = id_hips;
//   id->unknownlist = unknownlist;
//   id->n = Mat->nrow;
//   id->nnz = Mat->nnz;
// 
//   if (debug && proc_id == 0) {printf("cpu %d, Sortie %s \n",proc_id,__FUNCTION__);}
// 
//   status = MPI_Barrier(MPI_COMM_WORLD);
//   
//   /*exit(-23);*/
// 
//   /* une resolution pour test travail usr le local nodelist*/
// //   x1 = (double *) malloc(sizeof(double)*Mat->nrow);
// //   x2 = (double *) malloc(sizeof(double)*Mat->nrow);
// //   for (bcl=0;bcl < Mat->nrow;bcl++) x1[bcl]=x2[bcl]=1.0;
// //   
// //   dHIPS_MatrixVectorProduct (id_hips, Mat->x, x1);
// //   
// //   ierr = dHIPS_SetRHS(id_hips, Mat->nrow, unknownlist, x1, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_FOOL);
// //   dHIPS_ExitOnError(ierr);
// //   ierr = dHIPS_GetSolution(id_hips,Mat->nrow,unknownlist, x2, dHIPS_ASSEMBLY_FOOL);
// 
// //   //sprintf(filename,"solver_fatorize_hips.RESO1.%d.txt",proc_id);
// //   //out=fopen(filename,"w");
// //   //fprintf(out, "%d  %d \n",idi, Mat->nrow);
// //   //for (bcl=0; bcl < Mat->nrow; bcl++) fprintf(out, "%d => %d \n",bcl,x2[bcl]);
// //   //fclose(out);
// 
//   ierr=0;
//   status=0;
//   return(status);
// #endif 
// #endif
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int factorize_hips_sequential(solver *slv,triplet *Mat) {
//     
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int status,ierr;
// #ifdef HIPS
//   dHIPS_STRUC_C *id;
//   INTS  sym_pattern, sym_matrix;
//   INTS id_hips, idnbr, i, j;
//   INTS *unknownlist, *hipsnodelist;
//   COEFD *x, *rhsloc;
//   INTS proc_id, n, ln;
//   INTL *ia, nnz;
//   INTS *ja;
//   COEFD *a;
//   INTS domsize, nproc;
//   INTS pbegin, pend;
//   INTS *mapptr,  *mapptr2,*mapp,*mapp2;
//   INTS p;
// 
//   int nglob, nnzglob, tmp, bcl, idi, info;
//   FILE *out;
//   char filename[50];
//   double *x1,*x2, tol;
//   
// 
//   /* nglob */
//   nglob=Mat->nrow;  
//   
//   idnbr = 100; /* total */
//   if (HIPSSOLVERNB == 0) {
//     ierr = dHIPS_Initialize(idnbr);
//     dHIPS_ExitOnError(ierr);
//     }
//   id_hips = HIPSSOLVERNB; /** id_hips of the linear system **/
//   HIPSSOLVERNB++;
//   
//   /** parameter domsize is an argument of testHIPS.ex **/
// #ifdef __DEBUG_HIPS
//   printf("%d,  factorize_hips solver numero=%d\n",proc_id,id_hips);
// #endif
//   domsize = nglob/nproc/2;
// #ifdef __DEBUG_HIPS
//   printf("%d,  factorize_hips domsize=%d\n",proc_id,domsize);
// #endif
//    dHIPS_SetOptionINT(id_hips, dHIPS_PARTITION_TYPE, 0); 
//    dHIPS_SetOptionINT(id_hips, dHIPS_DOMSIZE, domsize); 
// /*   /\* OPTIONS*\/ */
// //  dHIPS_SetDefaultOptions (id_hips, dHIPS_HYBRID);
//   dHIPS_SetDefaultOptions (id_hips, dHIPS_ITERATIVE);
//   tol=1e-12;
//   dHIPS_SetOptionREAL(id_hips,dHIPS_DROPTOL0,tol);
//   dHIPS_SetOptionREAL(id_hips,dHIPS_DROPTOL1,tol);
// #ifdef VERBOSE  
//   info=4;
// #else
//   infor=0
// #endif  
//   dHIPS_SetOptionINT (id_hips, dHIPS_VERBOSE, info);
//   /*  dHIPS_SetOptionREAL(id_hips, INTS number, REAL value);*/
//   sym_matrix=0;
//   dHIPS_SetOptionINT (id_hips, dHIPS_SYMMETRIC, sym_matrix);
//   /** C : numbering starts from 0 **/
//   /* but here starts from 1 **/
//   dHIPS_SetOptionINT(id_hips, dHIPS_FORTRAN_NUMBERING, 0);
//   
//   /***************************************************/
//   /*            ENTER THE GRAPH                      */
//   /***************************************************/
//   ierr = dHIPS_GraphBegin(id_hips,nglob, Mat->nnz);
//   dHIPS_ExitOnError(ierr);
//   for (bcl=0;bcl<Mat->nnz;bcl++) {
//     ierr = dHIPS_GraphEdge(id_hips, Mat->i[bcl]-1, Mat->j[bcl]-1);
//   }
//   ierr = dHIPS_GraphEnd(id_hips);
//   dHIPS_ExitOnError(ierr);
//   status = MPI_Barrier(MPI_COMM_WORLD);
//   
// #ifdef __DEBUG_HIPS
//   printf("%d,  factorize_hips GraphEnd OK...\n",proc_id);
// #endif
//   
//   /***************************************************/
//   /*            ENTER A USER PARTITION               */
//   /***************************************************/
//   /*** MARCHE PAS POUR LE MOMENT ***/
//   /** Only the master processor (0 by default) needs to enter the partition **/
//   /* collect inforamtion on proc 0 */
//       mapptr = (INTS *)malloc(sizeof(INTS)*(nproc+1));
//       mapptr2= (INTS *)malloc(sizeof(INTS)*(nproc+1));
//       mapp = (INTS *)malloc(sizeof(INTS)*nglob);
//       mapp2= (INTS *)malloc(sizeof(INTS)*nglob);
// 
//       for (bcl=0;bcl<nproc;bcl++)   mapptr2[bcl]=0;
//       for (bcl=0;bcl<nglob;bcl++)   mapp2[bcl]=0;
//       for (bcl=0;bcl<nproc;bcl++)   mapptr[bcl]=0;
//       for (bcl=0;bcl<nglob;bcl++)   mapp[bcl]=0;
//       for (bcl=proc_id+1;bcl<nproc+1;bcl++) mapptr2[bcl]=Mat->nrow;
// 
// ///      MPI_Allreduce( mapptr2, mapptr, nproc+1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
//       
// // #ifdef __DEBUG_HIPS
// //       for (bcl=0;bcl<nproc+1;bcl++) printf("...%d, ",mapptr[bcl]);
// //       printf("\n");
// // #endif
//  
//       idi=0;
//       mapp2[idi+mapptr[proc_id]] =  Mat->i[0]-1;
//       for (bcl=1; bcl < Mat->nnz; bcl++) {
//         if ( mapp2[idi+mapptr[proc_id]] != (Mat->i[bcl]-1) ) {
//           idi++;
//           /*printf("bcl=%d, idi=%d, [idi+mapptr[proc_id]=%d nglob=%d \n",bcl,idi,idi+mapptr[proc_id], nglob);*/
//           mapp2[idi+mapptr[proc_id]] = Mat->i[bcl]-1;        
//         }
//       }
// ///      MPI_Allreduce( mapp2, mapp, nglob, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
//       /* Set the unknownlist */
//       unknownlist= (INTS *)malloc(sizeof(INTS)*(Mat->nrow));
//       unknownlist[0]=Mat->i[0]-1;
//       idi=0;
//       for (bcl=1; bcl < Mat->nnz; bcl++) {
//         if ( unknownlist[idi] != (Mat->i[bcl]-1) ) {
//           idi++;
//           unknownlist[idi] = Mat->i[bcl]-1;        
//         }
//       }
//        /* Enter Partition */
//       if (proc_id ==0) {
//         ierr = dHIPS_SetPartition(id_hips, nproc, mapptr, mapp);
//         dHIPS_ExitOnError(ierr);
//       }
//       
// 
//       free(mapptr2);
//       free(mapptr);
//       free(mapp2);
//       free(mapp);
// 
//   ierr = dHIPS_AssemblyBegin(id_hips, Mat->nnz, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_FOOL, sym_matrix);
//   dHIPS_ExitOnError(ierr);
//   
// // #ifdef __DEBUG_HIPS
// //   printf("%d,  factorize_hips AssemblyBegin OK...\n",proc_id);
// // #endif    
//   sprintf(filename,"solver_fatorize_hips.4.%d.txt",proc_id);
//   for(bcl=0;bcl<Mat->nnz;bcl++) {
//     ierr = dHIPS_AssemblySetValue(id_hips, Mat->i[bcl]-1, Mat->j[bcl]-1, Mat->x[bcl]);
//     dHIPS_ExitOnError(ierr);
//   }
//   ierr = dHIPS_AssemblyEnd(id_hips);
//   dHIPS_ExitOnError(ierr);
//     
//   id->id_hips = id_hips;
//   id->unknownlist = unknownlist;
//   id->n = Mat->nrow;
//   id->nnz = Mat->nnz;
// 
// // #ifdef __DEBUG_HIPS
// //   printf("%d, Sortie factorize_hips \n",proc_id);
// // #endif    
// //   status = MPI_Barrier(MPI_COMM_WORLD);
//   /*exit(-23);*/
// 
//   /* une resolution pour test travail usr le local nodelist*/
//   x1 = (double *) malloc(sizeof(double)*Mat->nrow);
//   x2 = (double *) malloc(sizeof(double)*Mat->nrow);
//   for (bcl=0;bcl < Mat->nrow;bcl++) x1[bcl]=1.0;
//   /*dHIPS_MatrixVectorProduct (id_hips, x, x2);*/
//   ierr = dHIPS_SetRHS(id_hips, Mat->nrow, unknownlist, x1, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_FOOL);
//   dHIPS_ExitOnError(ierr);
//   ierr = dHIPS_GetSolution(id_hips,Mat->nrow,unknownlist, x2, dHIPS_ASSEMBLY_FOOL);
// 
//   ierr=0;
//   status=0;
//   return(status);
// #endif
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int factorize_lapack(solver *slv,triplet *Mat) {
//     
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// #ifdef LAPACK
//   int status; 
// 
//   slv->Mat=Mat;
//   Mat->pivot = (int *) malloc (Mat->nrow * sizeof (int)) ;
// /*   printf("Appel a dgbtrf Mat->lda=%d \n",Mat->lda); */
// /*   printf("               Mat->ml=%d, Mat->mu=%d \n", */
// /*          Mat->ml, Mat->mu); */
// /*   for (bcl=0; bcl < Mat->nrow; bcl++) { */
// /*     printf(" Mat->pivot[%d]=%d \n",bcl,Mat->pivot[bcl]); */
// /*   } */
//   dgbtrf_(         &Mat->nrow, &Mat->nrow,& Mat->ml, &Mat->mu, 
//                    Mat->a , &Mat->lda, Mat->pivot,&status);
//   /*  printf("Apres dgbtrf \n");*/
// /*   for (bcl=0; bcl < Mat->nrow; bcl++) { */
// /*     printf(" Mat->pivot[%d]=%d \n",bcl,Mat->pivot[bcl]); */
// /*   } */
//   /*DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )*/
//   return(status);
// #endif
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int factorize_spdomestic(solver *slv,triplet *Mat) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   int status=0,n,order=3,bcl;
//   double tol=1.0; /* unsymetrique usage*/
//   cs *A = cs_calloc (1, sizeof (cs)) ;    /* allocate the cs struct */  ; 
//   css *S ;
//   csn *N ;
//   
//   Mat->A = A;
//   A->nzmax = Mat->nnz;
//   A->m =  Mat->nrow;
//   A->n =  Mat->ncol;
//   A->p =  Mat->i;
//   A->i =  Mat->j;
//   A->x =  Mat->x;
//   A->nz=-1;
//   n = A->n;
//   S = cs_sqr (order, A, 0) ;                    /* ordering and symbolic analysis */
//   N = cs_lu (A, S, tol) ;                    /* numeric LU factorization */
//   bcl = N->L->nz;
//   Mat->S = S;
//   Mat->N = N;
//   slv->Mat = Mat;
//   return(status);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int factorize_pastix(solver *slv,triplet *Mat) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/    
//   int status;
// #ifdef PASTIX
//   PASTIX_STRUC_C *id;
//   char *solver_name;
//   /*pastix_data_t  *pastix_data = NULL;*/
//   double *Ax;
//   pastix_int_t *Ap, *Ai;
//   int  n,nnz,nbthread;
//   double *rhs         = NULL; /* right hand side                                           */
//   int               verbosemode;        /* Level of verbose mode (0, 1, 2)                           */
//   char             *type        = NULL; /* type of the matrix                                        */
//   int      bcl, bcl2;
//   int    flagsym=0;
//   pastix_int_t    ncol;               /* Size of the matrix                                        */
//   int verbose=0;
//   
//   solver_name=slv->name;
//   if (strcmp(solver_name,"PASTIX") == 0) {
//     if(verbose==1) printf("factorize PASTIX ...\n");
//     id = slv->parameters;
//     slv->Mat = Mat; 
//     Ap=Mat->i;
//     Ai=Mat->j;
//     Ax=Mat->x;
//     n=Mat->nrow;
//     nnz=Mat->nnz;
// /*
//     pastix_data = id->pastix_data;
//     perm = id->perm;
//     invp = id->invp;
//     iparm = id->iparm;
//     dparm = id->dparm;
// */
//   
//     /*pastix_data=id->pastix_data;*/
// 
//   /*******************************************/
//   /* Initialize parameters to default values */
//   /*******************************************/
// 
//     verbosemode=0;    
//     
// //    printf(" check matrix ncol=%d \n",n);
//     d_pastix_checkMatrix(MPI_COMM_WORLD, verbosemode,
//                      API_SYM_NO, API_YES,
//                      n, &Ap, &Ai, &Ax, NULL,1);
// 
//   
//     id->iparm[IPARM_MODIFY_PARAMETER] = API_NO;
//     d_pastix(&(id->pastix_data), MPI_COMM_WORLD,
//          n, Ap, Ai, Ax,
//          id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
//   /*******************************************/
//   /*       Customize some parameters         */
//   /*******************************************/
// #ifdef OMP_H
//   nbthread = omp_get_max_threads(); 
//   printf("#using %d threads (omp_get_max_threads)...................... \n",nbthread);
// #else
//   nbthread = 1;  
//   printf("#using only 1 thread!...................... \n",nbthread);
// #endif
//   verbosemode=0;    
//   id->iparm[IPARM_SYM]           = API_SYM_NO;
//   id->iparm[IPARM_FACTORIZATION] = API_FACT_LU;
//   id->iparm[IPARM_THREAD_NBR] = nbthread;
//   id->iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
//   id->iparm[IPARM_LEVEL_OF_FILL] = 0;
//   id->iparm[IPARM_START_TASK] = API_TASK_ORDERING;
//   id->iparm[IPARM_END_TASK]   = API_TASK_NUMFACT;
// //   id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-8;
//   id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-12;
//   id->iparm[IPARM_VERBOSE]             = verbosemode;
// 
//   /*iparm[IPARM_END_TASK]   = API_TASK_REFINE;*/
// //  printf("Ici 1 ....\n");
//   id->perm = malloc(n*sizeof(pastix_int_t));
//   id->invp = malloc(n*sizeof(pastix_int_t));
//   /*rhs = malloc(n*sizeof(double));*/
// //  printf("Ici 2 ....\n");
//   
// /*Factorisation.........*/ 
//   d_pastix(&(id->pastix_data), MPI_COMM_WORLD,
//          n, Ap, Ai, Ax,
//          id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
// //  printf("Facto OK...\n");
//   
//   status=id->iparm[IPARM_ERROR_NUMBER];
// 
//   if (status < 0)  { printf("Factorisation status=%d\n",status);  }
//  
// 
// /*   Mat->i = Ap; */
// /*   Mat->j = Ai; */
// /*   Mat->x = Ax; */
// /*   slv->Mat = Mat; */
// 
// /*  /\*TEST SOLVE *\/ */
// /*   Mat=slv->Mat; */
// /*     Ap=Mat->i; */
// /*     Ai=Mat->j; */
// /*     Ax=Mat->x; */
// /*     n=Mat->nrow; */
// /*     nnz=Mat->nnz; */
//  
// /*    printf("Première resolution......\n"); */
// /*   id->iparm[IPARM_START_TASK] = API_TASK_SOLVE; */
// /*   id->iparm[IPARM_END_TASK]   = API_TASK_REFINE; */
// /*   d_pastix(&(id->pastix_data), MPI_COMM_WORLD, */
// /*          n, Ap, Ai, Ax, */
// /*          id->perm, id->invp, rhs, 1, id->iparm, id->dparm);  */
// /*   id->iparm[IPARM_START_TASK] = API_TASK_SOLVE; */
// /*   id->iparm[IPARM_END_TASK]   = API_TASK_REFINE; */
// /*   d_pastix(&(id->pastix_data), MPI_COMM_WORLD, */
// /*          n, Ap, Ai, Ax, */
// /*          id->perm, id->invp, rhs, 1, id->iparm, id->dparm);  */
// 
//  
//   status=0;        
//   return(status);
//   }
// #else
//   status=-1;
//   printf("Please compile poc-solvers library with -DPASTIX \n");
//   return(status);
// #endif 
//   return(status); 
// }
// 
// // /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// // 
// //  int factorize_pastix_parallel_new(solver *slv,triplet *Mat, int verbose) {
// // 
// // /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// //   int status;
// // #ifdef PASTIX
// // #ifdef HAVE_MPI
// //   PASTIX_STRUC_C *id;
// //   char *solver_name;
// //   double *Ax;
// //   pastix_int_t *Ap, *Ai;
// //   int  n,nnz,nbthread;
// //   double *rhs         = NULL; /* right hand side                     */
// //   int    verbosemode;         /* Level of verbose mode (0, 1, 2)     */
// //   char   *type        = NULL; /* type of the matrix                  */
// //   int    bcl, bcl2;
// //   int    flagsym=0;
// //   pastix_int_t    ncol;       /* Size of the matrix                  */
// //   // Global variables
// //   int ncolglob, nnzglob;
// //   int *Apglob, *Aiglob, *colptr, *rows;
// //   double *Axglob, *values;
// //   char filemat[50];
// //   FILE *f_out;
// //   int             mpid,sz;
// // 
// //   solver_name=slv->name;
// //   
// // //   if (strcmp(solver_name,"PASTIX") == 0) {
// //   
// //   id = slv->parameters;
// //   slv->Mat = Mat; 
// //   Ap=Mat->i;
// //   Ai=Mat->j;
// //   Ax=Mat->x;
// //   n=Mat->nrow;
// //   nnz=Mat->nnz;
// // 
// //   MPI_Comm_rank(MPI_COMM_WORLD, &mpid);
// //   MPI_Comm_size(MPI_COMM_WORLD, &sz);
// // 
// // // /*------------------------------------------------------------------------------
// // //   Dispatch Matrix global matrix */
// // //   d_coo2glob(n,nnz,Ap,Ai,Ax, &ncolglob, &nnzglob, &Apglob, &Aiglob, &Axglob);
// // //   
// // //   free(Ap);
// // //   free(Ai);
// // //   free(Ax);
// //    
// // /*------------------------------------------------------------------------------
// //   convert structure to CSC */
// //   d_globcoo2csc(ncolglob,nnzglob,Ap , Ai, Ax , &colptr, &rows, &values);
// //   
// //   free(Apglob);
// //   free(Aiglob);
// //   free(Axglob);
// //    
// //   Ap = colptr;
// //   Ai = rows;
// //   Ax = values;
// //    
// //   Mat->nrow = ncolglob;
// //   Mat->nnz = nnzglob;
// //    
// //   /*******************************************/
// //   /* Initialize parameters to default values */
// //   /*******************************************/
// // 
// //   verbosemode=verbose;  
// //   
// //   d_pastix_checkMatrix(MPI_COMM_WORLD, verbosemode,
// //                      API_SYM_NO, API_YES, ncolglob, &Ap, &Ai, &Ax, NULL,1);
// // 
// //   Mat->i=Ap;
// //   Mat->j=Ai;
// //   Mat->x=Ax;
// //   Mat->comptype=CSC;
// //     
// //   id->iparm[IPARM_MODIFY_PARAMETER] = API_NO;
// //   d_pastix(&(id->pastix_data), MPI_COMM_WORLD,
// //            ncolglob, Ap, Ai, Ax, id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
// //   
// //   /*******************************************/
// //   /*       Customize some parameters         */
// //   /*******************************************/
// // #ifdef OMP_H
// //   nbthread = omp_get_max_threads(); 
// //   printf("#using %d threads (omp_get_max_threads)...................... \n",nbthread);
// // #else
// //   nbthread = 1;  
// //   printf("#using only 1 thread!...................... \n",nbthread);
// // #endif
// //   verbosemode=verbose;    
// //   id->iparm[IPARM_SYM]           = API_SYM_NO;
// //   id->iparm[IPARM_FACTORIZATION] = API_FACT_LU;
// //   id->iparm[IPARM_THREAD_NBR] = nbthread;
// //   id->iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
// //   id->iparm[IPARM_LEVEL_OF_FILL] = 0;
// //   id->iparm[IPARM_START_TASK] = API_TASK_ORDERING;
// //   id->iparm[IPARM_END_TASK]   = API_TASK_NUMFACT;
// // //   id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-8;
// //   id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-12;
// //   id->iparm[IPARM_VERBOSE]             = verbosemode;
// // 
// //   /*iparm[IPARM_END_TASK]   = API_TASK_REFINE;*/
// // 
// //   id->perm = malloc(ncolglob*sizeof(pastix_int_t));
// //   id->invp = malloc(ncolglob*sizeof(pastix_int_t));
// //   
// // /*Factorisation.........*/ 
// //   if(verbose==1) printf("%s pastix factorisation start\n",__FUNCTION__);
// //   
// //   d_pastix(&(id->pastix_data), MPI_COMM_WORLD,
// //            ncolglob, Ap, Ai, Ax, id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
// //   
// //   if(verbose==1) printf("%s pastix factorisation end\n",__FUNCTION__);
// //   
// //   status=id->iparm[IPARM_ERROR_NUMBER];
// // 
// //   if (status < 0)  {
// //     printf("d_pastix factorization status=%d\n",status);  
// //     return(-1);
// //     }
// //  
// //   status=0;        
// //   return(status);
// // 
// // #endif 
// // #else
// //   status=-1;
// //   printf("Please compile poc-solvers library with -DPASTIX \n");
// //   return(status);
// // #endif 
// //   return(status); 
// // }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//  int factorize_pastix_parallel(solver *slv,triplet *Mat, int verbose) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   int status;
// #ifdef PASTIX
// #ifdef HAVE_MPI
//   PASTIX_STRUC_C *id;
//   char *solver_name;
//   double *Ax;
//   pastix_int_t *Ap, *Ai;
//   int  n,nnz,nbthread;
//   double *rhs         = NULL; /* right hand side                     */
//   int    verbosemode;         /* Level of verbose mode (0, 1, 2)     */
//   char   *type        = NULL; /* type of the matrix                  */
//   int    bcl, bcl2;
//   int    flagsym=0;
//   pastix_int_t    ncol;       /* Size of the matrix                  */
//   // Global variables
//   int ncolglob, nnzglob;
//   int *Apglob, *Aiglob, *colptr, *rows;
//   double *Axglob, *values;
//   char filemat[50];
//   FILE *f_out;
//   int             mpid,sz;
// 
//   solver_name=slv->name;
//   
// //   if (strcmp(solver_name,"PASTIX") == 0) {
//   
//   id = slv->parameters;
//   slv->Mat = Mat; 
//   Ap=Mat->i;
//   Ai=Mat->j;
//   Ax=Mat->x;
//   n=Mat->nrow;
//   nnz=Mat->nnz;
// 
//   MPI_Comm_rank(MPI_COMM_WORLD, &mpid);
//   MPI_Comm_size(MPI_COMM_WORLD, &sz);
// 
// /*------------------------------------------------------------------------------
//   Dispatch Matrix global matrix */
//   d_coo2glob(n,nnz,Ap,Ai,Ax, &ncolglob, &nnzglob, &Apglob, &Aiglob, &Axglob);
//   
//   free(Ap);
//   free(Ai);
//   free(Ax);
//    
// /*------------------------------------------------------------------------------
//   convert structure to CSC */
//   d_globcoo2csc(ncolglob,nnzglob,Apglob , Aiglob, Axglob , &colptr, &rows, &values);
//   
//   free(Apglob);
//   free(Aiglob);
//   free(Axglob);
//    
//   Ap = colptr;
//   Ai = rows;
//   Ax = values;
//    
//   Mat->nrow = ncolglob;
//   Mat->nnz = nnzglob;
//    
//   /*******************************************/
//   /* Initialize parameters to default values */
//   /*******************************************/
// 
//   verbosemode=verbose;  
//   
//   d_pastix_checkMatrix(MPI_COMM_WORLD, verbosemode,
//                      API_SYM_NO, API_YES, ncolglob, &Ap, &Ai, &Ax, NULL,1);
// 
//   Mat->i=Ap;
//   Mat->j=Ai;
//   Mat->x=Ax;
//   Mat->comptype=CSC;
//     
//   id->iparm[IPARM_MODIFY_PARAMETER] = API_NO;
//   d_pastix(&(id->pastix_data), MPI_COMM_WORLD,
//            ncolglob, Ap, Ai, Ax, id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
//   
//   /*******************************************/
//   /*       Customize some parameters         */
//   /*******************************************/
// #ifdef OMP_H
//   nbthread = omp_get_max_threads(); 
//   printf("#using %d threads (omp_get_max_threads)...................... \n",nbthread);
// #else
//   nbthread = 1;  
//   printf("#using only 1 thread!...................... \n",nbthread);
// #endif
//   verbosemode=verbose;    
//   id->iparm[IPARM_SYM]           = API_SYM_NO;
//   id->iparm[IPARM_FACTORIZATION] = API_FACT_LU;
//   id->iparm[IPARM_THREAD_NBR] = nbthread;
//   id->iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
//   id->iparm[IPARM_LEVEL_OF_FILL] = 0;
//   id->iparm[IPARM_START_TASK] = API_TASK_ORDERING;
//   id->iparm[IPARM_END_TASK]   = API_TASK_NUMFACT;
// //   id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-8;
//   id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-12;
//   id->iparm[IPARM_VERBOSE]             = verbosemode;
// 
//   /*iparm[IPARM_END_TASK]   = API_TASK_REFINE;*/
// 
// //   id->perm = malloc(ncolglob*sizeof(pastix_int_t));
// //   id->invp = malloc(ncolglob*sizeof(pastix_int_t));
//   id->perm = new pastix_int_t[ncolglob];
//   id->invp = new pastix_int_t[ncolglob];
//   
// /*Factorisation.........*/ 
//   if(verbose==1) printf("%s pastix factorisation start\n",__FUNCTION__);
//   
//   d_pastix(&(id->pastix_data), MPI_COMM_WORLD,
//            ncolglob, Ap, Ai, Ax, id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
//   
//   if(verbose==1) printf("%s pastix factorisation end\n",__FUNCTION__);
//   
//   status=id->iparm[IPARM_ERROR_NUMBER];
// 
//   if (status < 0)  {
//     printf("d_pastix factorization status=%d\n",status);  
//     return(-1);
//     }
//  
//   status=0;        
//   return(status);
// 
// #endif 
// #else
//   status=-1;
//   printf("Please compile poc-solvers library with -DPASTIX \n");
//   return(status);
// #endif 
//   return(status); 
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// /*************************************************************************/
// /*************************************************************************/
// /*=======================================================================*/
// /*====================== FIN FACTORISE ==================================*/
// /*=======================================================================*/
// 
// 
// 
// /*=======================================================================*/
// /*                              SOLVE                                    */
// /*************************************************************************/
// /*=======================================================================*/
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int solve(solver *slv,double *RHS,int transp, int verbose) {
//   
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   char *solver_name;
//   int rtn;
//   double cputime1,cputime2;
//   int proc_id, ierr;
//   bool debug=false;
//   
// #ifdef  HAVE_MPI
//   if(MPI_Init_thread_DONE==1) {
//     ierr = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
//     cputime1 = MPI_Wtime();
//     }
// #endif
//   
//   solver_name=slv->name;
//      
//   if(verbose) printf("%s real-system solver, transpose=%d \n",solver_name, transp);
// 
//   if (strcmp(solver_name,"MUMPS") == 0) { 
//     rtn = solve_mumps(slv, RHS,transp,verbose,false);
//     }
//   else if  (strcmp(solver_name,"UMFPACK") == 0) {
//     rtn=solve_umfpack(slv, RHS,transp);
//     }
//   else if  (strcmp(solver_name,"LAPACK") == 0) {
//     rtn=solve_lapack(slv, RHS,transp);
//     }
//   else if  (strcmp(solver_name,"SpDOMESTIC") == 0) {
//     rtn=solve_spdomestic(slv, RHS,transp);
//     }
//   else if  (strcmp(solver_name,"HIPS") == 0) {
//     rtn=solve_hips(slv, RHS,transp);
//     }
//   else if  (strcmp(solver_name,"PASTIX") == 0) {
// #ifdef HAVE_MPI    
//     MPI_Barrier(MPI_COMM_WORLD);
//     rtn=solve_pastix_parallel(slv, RHS,transp);
// #else
//     rtn=solve_pastix(slv, RHS,transp);
// #endif     
//     }
//   else if  (strcmp(solver_name,"PASTIX-SEQUENTIAL") == 0) {
//     rtn=solve_pastix(slv, RHS,transp);
//     }
//   else {
//     printf("solve, unknown solver : %s \n",solver_name);
//     rtn=-1;
//     }
//      
// #ifdef  HAVE_MPI
// //   if(MPI_Init_thread_DONE==1 && verbose==1) {
//   if(MPI_Init_thread_DONE==1) {
//     cputime2 = MPI_Wtime();
//     if (verbose==1 && proc_id == 0) {
//       printf("--------------------------------------------\n");
//       printf("%s: %s solving ellapsed time= %e\n",__FUNCTION__, solver_name, cputime2-cputime1);
//       printf("--------------------------------------------\n");
//       }
//     }
// #endif
// 
//   return(rtn);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int solve_mumps(solver *slv,double *RHS, int transp, int verbose, bool debug) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// #ifdef MUMPS
//   DMUMPS_STRUC_C *id;
//   char *solver_name;
//   int ierr,rank;
//   double cput1, cput2, cput3;
//   
// #ifdef HAVE_MPI
//   cput1 = MPI_Wtime();
// #endif
//   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   solver_name=slv->name;
//   if (strcmp(solver_name,"MUMPS") == 0) {
//     id = slv->parameters;
//     if (transp == 1) {
//       (*id).ICNTL(9)=2;
//       }
//     if (rank == 0) {
//       /*printf("Solve MUMPS ...\n");*/
//       (*id).rhs = RHS;
//       }
//     (*id).job=3;
//     dmumps_c(id);
//     }
// #ifdef HAVE_MPI
//   cput2 = MPI_Wtime(); 
//   cput3 = cput2-cput1;
//   if (rank == 0) {
//     if(verbose==1) printf("  DMUMPS resolution time =%lf\n",cput3);
//     }
// #endif
//   ierr=(*id).info[0];
//   if(ierr==-1) {
//     printf("cpu %d: DMUMPS solver failed, status=%d due to cpu %d\n",rank, ierr, id->info[1]);
//     }
//   if(ierr==-22) {
//     printf("cpu %d: DMUMPS solver failed, status=%d array fault %d\n",rank, ierr, id->info[1]);
//     }
//   return(ierr);
// #endif
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// int solve_umfpack(solver *slv,double *RHS,int transp) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   int status,bcl;
// #ifdef UMFPACK
//   UMFPACK_STRUC_C *id;
//   char *solver_name;
//   double *Info, *Control, *Ax, *tmp;
//   int *Ap, *Ai,n,nnz;
//   void *Symbolic, *Numeric ;
//   triplet *Mat;
//   int INCXY;
//  
//   status = 0;
//   solver_name=slv->name;
//   if (strcmp(solver_name,"UMFPACK") == 0) {
//     id = slv->parameters;
//     Symbolic=id->Symbolic;
//     Numeric=id->Numeric;
//     Info=id->Info;
//     Control=id->Control;
//     Mat=slv->Mat;
//     Ap=Mat->i;
//     Ai=Mat->j;
//     Ax=Mat->x;
//     n=Mat->nrow;
//     nnz=Mat->nnz;
//     tmp = (double *) malloc (n * sizeof (double)) ;
//     INCXY=1;
//     /* BLASC */
//     /*dcopy(n,RHS,1,tmp,1);*/
//     /*BLASF*/
//     /*dcopy_(n,RHS,INCXY,tmp,INCXY);*/
//     /*no blas Xcopy*/
//     for (bcl=0;bcl<n;bcl++) tmp[bcl]=RHS[bcl];
//     if (transp == 1) status = umfpack_di_solve (UMFPACK_At, Ap, Ai, Ax, RHS, tmp,
//       Numeric, Control, Info) ;
//     else status = umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, RHS, tmp,
//       Numeric, Control, Info) ;
//     
//     if (status < 0) {
//       umfpack_di_report_info (Control, Info) ;
//       umfpack_di_report_status (Control, status) ;
//       /*error ("umfpack_di_solve failed") ;*/
//       }
//     free(tmp);
//     return(status);
//   }
// #endif
//   return(status);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// int solve_lapack(solver *slv,double *RHS,int transp) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// #ifdef LAPACK
//   int status,nrhs=1,trans_len;
//   char cjob = 'N';
//   /*  double *B;*/
//   triplet *Mat;
// 
//   Mat=slv->Mat;
//   /*  B = (double *) malloc (Mat->nrow * sizeof (double)) ;*/
//   /*printf("Appel a dgbtrs_ \n");*/
// /*   for (bcl=0; bcl < Mat->nrow; bcl++) { */
// /*     printf(" Mat->pivot[%d]=%d \n",bcl,Mat->pivot[bcl]); */
// /*   } */
// /*   for (bcl=0; bcl < Mat->nrow; bcl++)  */
// /*     printf("RHS[%d]=%e \n",bcl,RHS[bcl]); */
//   dgbtrs_(&cjob,  &Mat->nrow, &Mat->ml, &Mat->mu, &nrhs, 
//           Mat->a, &Mat->lda,Mat->pivot, RHS, &Mat->nrow, &status,trans_len);
//   /*for (bcl=0;bcl<Mat->nrow;bcl++) RHS[bcl]=B[bcl];*/
//   /*printf("Appel a dgbtrs_ status=%d\n",status);*/
//   /*for (bcl=0; bcl < Mat->nrow; bcl++) 
//     printf("RHS[%d]=%e \n",bcl,RHS[bcl]);*/
//   /*  free(B);*/
//   return(status);
//  
// #endif
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// int solve_spdomestic(solver *slv,double *RHS,int transp) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   int status=0,bcl,n;
//   triplet *Mat;
//   cs *A;
//   css *S ;
//   csn *N ;
//   double *x;
// 
//   Mat = slv->Mat;
//   A = Mat->A;
//   S = Mat->S;
//   N = Mat->N;
//   n = A->n ;
//   x = cs_malloc (n, sizeof (double)) ;    /* get workspace */
// 
// 
//   /*printf("\n");*/
//   /*printf("....... SOLVE .............\n");*/
//   bcl = N->L->nz;
//   /*  for(bcl=0;bcl<A->n;bcl++) printf("RHS[%d]=%e\n",bcl,RHS[bcl]);*/
// 
//   cs_ipvec (N->pinv, RHS, x, n) ;            /* x = b(p) */
//   if (transp) {
//     cs_utsolve (N->U, x) ;                    /* x = U\x */
//     cs_ltsolve (N->L, x) ;                    /* x = L\x */
//   }
//   else {
//     cs_lsolve (N->L, x) ;                    /* x = L\x */
//     cs_usolve (N->U, x) ;                    /* x = U\x */
//   }
//   cs_ipvec (S->q, x, RHS, n) ;                   /* b(q) = x */
//   cs_free (x) ;
//   
// /*   printf("A->P : "); */
// /*   for (bcl=0;bcl<A->n;bcl++) printf("%d ",A->p[bcl]); */
// /*   printf("\n"); */
// /*   for (bcl=0;bcl<A->nzmax;bcl++) { */
// /*     printf(" A : i=%d, x=%e \n",A->i[bcl],A->x[bcl]); */
// /*   } */
// /*   printf("L->P : "); */
// /*   for (bcl=0;bcl<N->L->n;bcl++) printf("%d ",N->L->p[bcl]); */
// /*   printf("\n"); */
// /*   for (bcl=0;bcl<N->L->nzmax;bcl++) { */
// /*     printf(" L : i=%d, x=%e \n",N->L->i[bcl],N->L->x[bcl]);  */
// /*   } */
// /*   printf("A->U : "); */
// /*   for (bcl=0;bcl<N->U->n;bcl++) printf("%d ",N->U->p[bcl]); */
// /*   printf("\n"); */
// /*   for (bcl=0;bcl<N->U->nzmax;bcl++) { */
// /*     printf(" U : i=%d, x=%e \n",N->U->i[bcl],N->U->x[bcl]);  */
// /*   } */
// /*   for(bcl=0;bcl<A->n;bcl++) printf("RHS[%d]=%e\n",bcl,RHS[bcl]); */
// 
//   return(status); 
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int solve_hips(solver *slv,double *RHS, int transp, bool debug) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// #ifdef HIPS
// #ifdef HAVE_MPI
//   dHIPS_STRUC_C*id;
//   char *solver_name;
//   int proc_id, myid=-1, ierr, status, n, nnz,bcl, nglob;
//   double cput1, cput2, cput3;
//   INTS id_hips;
//   INTS *unknownlist;
//   char filename[50];
//   FILE *out;
//   double *x,*x2;
//   
// 
//   cput1 = MPI_Wtime();
//   ierr  = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
// 
// // #ifdef __DEBUG_HIPS
// //   printf("%d, Entreee solve_hips \n",proc_id);
// // #endif
//   
//   id = slv->parameters;
//   id_hips =id->id_hips;
//   n=id->n;
// 
// // #ifdef __DEBUG_HIPS
// //   printf("%d, solve_hips resolution du system numero %d, n=%d \n",proc_id,id_hips,n);
// // #endif
// 
//   nnz=id->nnz;
//   unknownlist = id->unknownlist;
//   
// //   sprintf(filename,"solver_solve_hips.1.%d.txt",proc_id);
// //   out=fopen(filename,"w");   
// //   for (bcl=0; bcl < n; bcl++) fprintf(out, "%d => %d \n",bcl,unknownlist[bcl]);
// //   fclose(out);
//   
//   /** TEst set RHS like x=1 */
// /*   x  = (double *) malloc(sizeof(double)*nnz); */
// /*   x2 = (double *) malloc(sizeof(double)*nnz); */
// /*   for (bcl=0;bcl < nnz;bcl++) x[bcl]=1.0; */
// 
// /*   dHIPS_MatrixVectorProduct (id_hips, x, x2); */
// /*   for (bcl=0;bcl < n;bcl++) RHS[bcl]=x2[bcl]; */
// 
//   int mode=1;
//   switch(mode) {
//     case 0:
// /*------------------------------------------------------------------------------
//       version locale */
// //       ierr = dHIPS_SetRHS (id_hips, n, unknownlist, RHS,  dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_FOOL);
//       ierr = dHIPS_SetRHS (id_hips, n, unknownlist, RHS,  dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_RESPECT);
//       break;
//       
// /*------------------------------------------------------------------------------
//     version globale*/
//     case 1:
//       ierr = dHIPS_SetGlobalRHS(id_hips,RHS,0,  dHIPS_ASSEMBLY_OVW);
//       break;
//     }
// 
// // #ifdef __DEBUG_HIPS
// //   printf("%d, solve_hips dHIPS_SetRHS OK.... \n",proc_id);
// // #endif
// 
// //   sprintf(filename,"solver_solve_hips.RHS.2.%d.txt",proc_id);
// //   out=fopen(filename,"w");   
// //   for (bcl=0; bcl < n; bcl++) fprintf(out, "%d => %d, %lf\n",bcl,unknownlist[bcl]-1,RHS[bcl]);
// //   fclose(out);
// 
//   cput1 = MPI_Wtime();
//   
//   switch(mode) {
//     case 0:
// /*------------------------------------------------------------------------------
//       version locale */
//       ierr = dHIPS_GetSolution(id_hips,n,unknownlist, RHS, dHIPS_ASSEMBLY_FOOL);
//       break;
//       
// /*------------------------------------------------------------------------------
//     version globale*/
//     case 1:
//       ierr = dHIPS_GetGlobalSolution(id_hips,RHS,-1);
//       break;
//     }
//   
// //   sprintf(filename,"solver_solve_hips.SOL.2.%d.txt",proc_id);
// //   out=fopen(filename,"w");   
// //   for (bcl=0; bcl < n; bcl++) fprintf(out, "%d => %d, %lf\n",bcl,unknownlist[bcl],RHS[bcl]);
// //   fclose(out);
// 
// //   dHIPS_SetOptionINT (id_hips,dHIPS_DISABLE_PRECOND, 1);
// 
//   cput2 = MPI_Wtime(); 
//   cput3 = cput2-cput1;
// #ifdef VERBOSE  
//   if (myid == 0) {
//     printf("  HIPS Resolution time =%lf \n",cput3);
//   }
// #endif  
// /*    ierr = MPI_Allreduce (  &n, &nglob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );   */
// /*    x  = (double *) malloc(sizeof(double)*nglob); */
// /*    ierr = dHIPS_GetGlobalSolution(id_hips,x,-1); */
// /*    sprintf(filename,"solver_solve_hips.GlobalSolution.%d.txt",proc_id); */
// /*    out=fopen(filename,"w"); */
// /*    for (bcl=0; bcl < nglob; bcl++) fprintf(out,"%d => %lf\n",bcl,x[bcl]); */
// /*    fclose(out);  */
// /*    free(x); */
//    /*RHS = x; */
//    ierr=0;
//   return(ierr);
// #endif
// #endif
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// int solve_pastix(solver *slv,double *RHS,int transp) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   int status,bcl,bcl2;
// #ifdef PASTIX
//   PASTIX_STRUC_C *id;
//   char *solver_name;
//   double *Ax;
//   int *Ap, *Ai,n,nnz;
//   triplet *Mat;
//   int INCXY;
//  
//   
//   
//   status = 0;
//   solver_name=slv->name;
//   if (strcmp(solver_name,"PASTIX") == 0) {
//     id = slv->parameters;
//     Mat=slv->Mat;
//     Ap=Mat->i;
//     Ai=Mat->j;
//     Ax=Mat->x;
//     n=Mat->nrow;
//     nnz=Mat->nnz;
// 
//    
//   id->iparm[IPARM_START_TASK] = API_TASK_SOLVE;
//   id->iparm[IPARM_END_TASK]   = API_TASK_REFINE;
// 
//   d_pastix(&(id->pastix_data), MPI_COMM_WORLD,
//          n, Ap, Ai, Ax,
//          id->perm, id->invp, RHS, 1, id->iparm, id->dparm);
//     return(status);
//   }
// #endif
//   return(status);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// int solve_pastix_parallel(solver *slv,double *RHS,int transp) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   int status,bcl,bcl2;
// #ifdef PASTIX
// #ifdef HAVE_MPI
//   PASTIX_STRUC_C *id;
//   char *solver_name;
//   double *Ax;
//   int *Ap, *Ai,n,nnz;
//   triplet *Mat;
//   int INCXY;
//   char filerhs[50];
//   int myid,sz,i,ierr;
//   FILE *out;
//   
// #ifdef HAVE_MPI
//      ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//      ierr = MPI_Comm_size(MPI_COMM_WORLD, &sz);
// #endif  
//   
//   status = 0;
//   solver_name=slv->name;
//   if (strcmp(solver_name,"PASTIX") == 0) {
//     id = slv->parameters;
//     Mat=slv->Mat;
//     Ap=Mat->i;
//     Ai=Mat->j;
//     Ax=Mat->x;
//     n=Mat->nrow;
//     nnz=Mat->nnz;
// 
// //    sprintf(filerhs,"RHS.%d.%d.txt",sz,myid);
// //    out=fopen(filerhs,"w");
// //    for (i=0; i <n; i++)   fprintf(out,"      %d -> %e \n",i,RHS[i]);
// //    fclose(out);
// //    MPI_Barrier(MPI_COMM_WORLD);
// //    printf("solve_pastix_parallel n=%d nnz=%d\n",n,nnz);
//    
//   id->iparm[IPARM_START_TASK] = API_TASK_SOLVE;
//   id->iparm[IPARM_END_TASK]   = API_TASK_REFINE;
// 
//   d_pastix(&(id->pastix_data), MPI_COMM_WORLD,
//          n, Ap, Ai, Ax,
//          id->perm, id->invp, RHS, 1, id->iparm, id->dparm);
// //    sprintf(filerhs,"SOL.%d.%d.txt",sz,myid);
// //    out=fopen(filerhs,"w");
// //    for (i=0; i <n; i++)   fprintf(out,"      %d -> %e \n",i,RHS[i]);
// //    fclose(out);
// 
//    return(status);
//   }
// #endif
// #endif
//   return(status);
// }
// /*=======================================================================*/
// /*====================  FIN SOLVE =======================================*/
// /*=======================================================================*/
// 
// /*=======================================================================*/
// /*                 Solver Terminate                                      */
// /*************************************************************************/
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// int free_solver(solver *slv)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   /*mumps_solver mumpsslv;*/
//   /*int len,ierr,myid;*/
//   int status=0;
// #ifdef MUMPS
//   DMUMPS_STRUC_C *id_mumps;
// #endif
// #ifdef UMFPACK
//   UMFPACK_STRUC_C *id_umf;
// #endif
// #ifdef PASTIX
//   PASTIX_STRUC_C *id_pastix;
//   pastix_data_t  *pastix_data = NULL;
//   pastix_int_t    *iparm;  /* integer parameters for pastix                             */
//   double          *dparm;  /* floating parameters for pastix                            */
//   double *Ax;
//   int *Ap, *Ai,n,nnz,nbthread;
//   pastix_int_t   *perm        = NULL; /* Permutation tabular                             */
//   pastix_int_t   *invp        = NULL; /* Reverse permutation tabular                     */
//   double *rhs         = NULL; /* right hand side                                         */
// #endif  /*double *Control;*/
// #ifdef HIPS
//   dHIPS_STRUC_C*id;
//   INTS id_hips;
//   int ierr;
// #endif
//   char *solver_name;
//   void *Symbolic, *Numeric ;
//   triplet *Mat;
//   FILE *out;
// 
//   solver_name=slv->name;
// //  printf("free solver : %s \n",solver_name);
//   Mat = slv->Mat;
// 
//   if (solver_name != NULL) {
// /* ---------------- FREE MUMPS SOLVER -------------------- */
// #ifdef MUMPS
//     if (strcmp(solver_name,"MUMPS") == 0) {
//       id_mumps = slv->parameters;
//       (*id_mumps).job=JOB_END;
//       dmumps_c(id_mumps);
// //       (*id_mumps).irn_loc=0;
// //       (*id_mumps).jcn_loc=0;
// //       (*id_mumps).a_loc=0;
//       Mat->i=0;
//       Mat->j=0;
//       Mat->x=0;
// /* ---------------- FREE MUMPS SOLVER -------------------- */
//       }
// #endif
// 
// /* ---------------- FREE UMFPACK SOLVER -------------------- */
// #ifdef UMFPACK
//     if (strcmp(solver_name,"UMFPACK") == 0) {
//       id_umf = slv->parameters;
//       Symbolic=id_umf->Symbolic;
//       Numeric=id_umf->Numeric;
//       umfpack_di_free_symbolic (&Symbolic) ;
//       umfpack_di_free_numeric (&Numeric) ;
// /* ---------------- FREE UMFPACK SOLVER -------------------- */
//       }
// #endif
// 
// /* ---------------- FREE LAPACK SOLVER -------------------- */
// #ifdef LAPACK
//     if (strcmp(solver_name,"LAPACK") == 0) {
//       free(Mat->a);
//       free(Mat->pivot);
// /* ---------------- FREE LAPACK SOLVER -------------------- */
//       }
// #endif
// 
// /* ---------------- FREE PASTIX SOLVER -------------------- */
// #ifdef PASTIX
//     if (strcmp(solver_name,"PASTIX") == 0) {
//       id_pastix = slv->parameters;
//       Ap=Mat->i;
//       Ai=Mat->j;
//       Ax=Mat->x;
//       id_pastix->iparm[IPARM_START_TASK] = API_TASK_CLEAN;
//       id_pastix->iparm[IPARM_END_TASK]   = API_TASK_CLEAN;
// 
//       d_pastix(&(id_pastix->pastix_data), MPI_COMM_WORLD,
//          n, Ap, Ai, Ax,
//          id_pastix->perm, id_pastix->invp, rhs, 1, 
//          id_pastix->iparm, id_pastix->dparm);
//       free(id_pastix->perm);
//       free(id_pastix->invp);
//       if(Mat->private_memory==1) {
//         free(Ax);
//         free(Ap);
//         free(Ai);
//         }
// /* ---------------- FREE PASTIX SOLVER -------------------- */
//       }
// #endif
// #ifdef HIPS
//     if (strcmp(solver_name,"HIPS") == 0) {
//       /************************************************/
//       /* Free zHIPS internal structure for problem id */
//       /************************************************/
//       id = slv->parameters;
//       id_hips =id->id_hips;
//       ierr = dHIPS_Clean(id_hips);
//       dHIPS_ExitOnError(ierr);
//       HIPSSOLVERNB--;
// 
//       /************************************************/
//       /* Free zHIPS internal structures               */
//       /* (to be done once only)                       */
//       /************************************************/
//       if(HIPSSOLVERNB==0) {
//         ierr = dHIPS_Finalize();
//         dHIPS_ExitOnError(ierr);
//         }
//       
//       }
// #endif
//     
// /* ---------------- FREE SpDOMESTIC SOLVER -------------------- */
//     if (strcmp(solver_name,"SpDOMESTIC") == 0) {
//       cs_sfree (Mat->S) ;
//       cs_nfree (Mat->N) ;
// /* ---------------- FREE SpDOMESTIC SOLVER -------------------- */
//       }
//     
//     free(slv->parameters);
//     free(slv->name);
//     /*free(slv);*/
//     slv->parameters = NULL;
//     slv->name = NULL;
//     }
//   else {
//     printf(" free already done with this solver \n");
//     }
//   
//   return(status);
// }
// 
// /*************************************************************************/
// /*=======================================================================*/
// /*======================FIN TERMINATE ===================================*/
// /*=======================================================================*/
// 
// 
// 
// /*************************************************************************/
// /*************************************************************************/
// /*=============== Complex Double precision routines =============================*/
// /*************************************************************************/
// /*************************************************************************/
// /*=======================================================================*/
// /*=======================================================================*/
// /*                            INIT                                       */
// /*************************************************************************/
// solver *initz_solver(char *solver_name,int typsym, int verbose DEFAULT_ZERO)
// {
//   solver *slv;
//   /*mumps_solver mumpsslv;*/
//   int ierr,myid;
//  
// #ifdef MUMPS 
//   ZMUMPS_STRUC_C *id_mumps;
//   ZMUMPS_STRUC_C *id;
// #endif
//   
// #ifdef UMFPACK
//   UMFPACK_STRUC_C *id_umf; 
//   double *Control;
// #endif
//   
// #ifdef HIPS
//   INTS id_hips, idnbr, i, j;
//   INTS *unknownlist;
//   COEFZ *x, *rhsloc;
//   INTS proc_id, n, ln;
//   INTL *ia, nnz;
//   INTS *ja;
//   COEFZ *a;
//   INTS domsize, nproc;
//   INTS pbegin, pend;
// #endif
//   
// #ifdef PASTIX
//   PASTIX_STRUC_C *id_pastix; 
//   int argc;
//   int **argv;
//   int             required;           /* MPI thread level required                                 */
//   int             provided;           /* MPI thread level provided                                 */
//   int             rank;
// #endif
//   
//   /*slv = (solver *) calloc(1,sizeof(solver));*/
//   slv = (solver *) malloc(sizeof(solver));
//   slv->name = strdup(solver_name);
// 
//   /*-------------------------------------------------------------*/
//   /* *********            INIT VARIOUS SOLVERS             *******/
// 
// /** ---------------- INIT MUMPS SOLVER -------------------- */
//   /* ---------------- INIT MUMPS SOLVER -------------------- */
//   if (strcmp(solver_name,"MUMPS") == 0) {
// #ifdef MUMPS
//     required = MPI_THREAD_MULTIPLE;
//     provided = -1;
//     if(MPI_Init_thread_DONE==0) {
//       MPI_Init_thread(&argc, &argv, required, &provided);
//       MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//       if (rank == 0) {
//         switch (provided) {
//           case MPI_THREAD_SINGLE:
//             printf("MPI_Init_thread level = MPI_THREAD_SINGLE\n");
//             break;
//           case MPI_THREAD_FUNNELED:
//             printf("MPI_Init_thread level = MPI_THREAD_FUNNELED\n");
//             break;
//           case MPI_THREAD_SERIALIZED:
//             printf("MPI_Init_thread level = MPI_THREAD_SERIALIZED\n");
//             break;
//           case MPI_THREAD_MULTIPLE:
//             printf("MPI_Init_thread level = MPI_THREAD_MULTIPLE\n");
//             break;
//           default:
//             printf("MPI_Init_thread level = ???\n");
//             break;
//           }
//         }
//       MPI_Init_thread_DONE=1;
//       }
// //     else {
// //       printf("%s : %s initialisation ; MPI already initialised \n",__FUNCTION__, solver_name);
// //       }
//     ierr = MPI_Barrier(MPI_COMM_WORLD);
//     ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
// #if 0
//     if(verbose==1) printf("init MUMPS solver (cpu=%d), start \n",myid);
// #endif
//     /*slv->parameters = (DMUMPS_STRUC_C *) malloc(sizeof(DMUMPS_STRUC_C));*/
//     slv->parameters = (ZMUMPS_STRUC_C *) calloc(1,sizeof(ZMUMPS_STRUC_C));
// //    id = (DMUMPS_STRUC_C *) calloc(1,sizeof(DMUMPS_STRUC_C*));
//     id =  slv->parameters;
//     /*id_mumps = slv->parameters;*/
//     /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
//     /*(*id_mumps).job=JOB_INIT; (*id_mumps).par=1; (*id_mumps).sym=typsym;*/
//     id->job=JOB_INIT; id->par=1; id->sym=0;
//     id->comm_fortran=USE_COMM_WORLD;
//     id->ICNTL(4)=4;
//     /*(*id_mumps).comm_fortran=USE_COMM_WORLD;*/
//     /*dmumps_c(id_mumps);*/
//     ierr = MPI_Barrier(MPI_COMM_WORLD);
//     zmumps_c(id);
//     if(verbose==1 && myid==0) printf("dmumps_par->instance_number %d \n",id->instance_number);
//     /*slv->parameters = id_mumps;*/
//     /*slv->parameters = id;*/
// #if 0
//     /* ---------------- FIN INIT MUMPS SOLVER -------------------- */
//     if(verbose==1) printf("init MUMPS solver, finish \n");
// #endif
// #else 
//     printf("Please compile poc-solvers library with -DMUMPS \n");
//     exit(-1);
// #endif
//     }
// /** --------------- INIT UMFPACK SOLVER ------------------ */
//   else if (strcmp(solver_name,"UMFPACK") == 0){
// #ifdef UMFPACK
//     slv->parameters = (UMFPACK_STRUC_C *) malloc(sizeof(UMFPACK_STRUC_C));
//     id_umf = slv->parameters;
//     Control=id_umf->Control;
//     umfpack_zi_defaults (id_umf->Control) ;
// /*     Control [UMFPACK_PRL] = 6 ; */
// /*     Control [UMFPACK_PRL] = 5 ; */
// /*     umfpack_di_report_control (Control) ; */
// /*    printf("Fin init umf\n");*/
// /* --------------- FIN INIT UMFPACK SOLVER ------------------ */
// #else 
//     printf("Please compile poc-solvers library with -DUMFPACK \n");
//     exit(-1);
// #endif
//    }
// 
// /** --------------- INIT HIPS SOLVER ------------------ */
//   else if (strcmp(solver_name,"HIPS") == 0){
// #ifdef HIPS
//     required = MPI_THREAD_MULTIPLE;
//     provided = -1;
//     if(MPI_Init_thread_DONE==0) {
//       MPI_Init_thread(&argc, &argv, required, &provided);
//       MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//       if (rank == 0) {
//         switch (provided) {
//           case MPI_THREAD_SINGLE:
//             printf("MPI_Init_thread level = MPI_THREAD_SINGLE\n");
//             break;
//           case MPI_THREAD_FUNNELED:
//             printf("MPI_Init_thread level = MPI_THREAD_FUNNELED\n");
//             break;
//           case MPI_THREAD_SERIALIZED:
//             printf("MPI_Init_thread level = MPI_THREAD_SERIALIZED\n");
//             break;
//           case MPI_THREAD_MULTIPLE:
//             printf("MPI_Init_thread level = MPI_THREAD_MULTIPLE\n");
//             break;
//           default:
//             printf("MPI_Init_thread level = ???\n");
//             break;
//           }
//         }
//       MPI_Init_thread_DONE=1;
//       }
// //     else {
// //       printf("%s : %s initialisation ; MPI already initialised \n",__FUNCTION__, solver_name);
// //       }
//     slv->parameters = (zHIPS_STRUC_C *) malloc(sizeof(zHIPS_STRUC_C));
// /* --------------- FIN INIT HIPS SOLVER ------------------ */
// #else 
//     printf("Please compile poc-solvers library with -DHIPS  \n");
//     exit(-1);
// #endif
//     }
// 
//   
//   /* --------------- INIT PASTX SOLVER ------------------ */
//   else if (strcmp(solver_name,"PASTIX") == 0){
// #ifdef PASTIX
//     required = MPI_THREAD_MULTIPLE;
//     provided = -1;
//     if(MPI_Init_thread_DONE==0) {
//       MPI_Init_thread(&argc, &argv, required, &provided);
//       MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//       if (rank == 0) {
//         switch (provided) {
//           case MPI_THREAD_SINGLE:
//             printf("MPI_Init_thread level = MPI_THREAD_SINGLE\n");
//             break;
//           case MPI_THREAD_FUNNELED:
//             printf("MPI_Init_thread level = MPI_THREAD_FUNNELED\n");
//             break;
//           case MPI_THREAD_SERIALIZED:
//             printf("MPI_Init_thread level = MPI_THREAD_SERIALIZED\n");
//             break;
//           case MPI_THREAD_MULTIPLE:
//             printf("MPI_Init_thread level = MPI_THREAD_MULTIPLE\n");
//             break;
//           default:
//             printf("MPI_Init_thread level = ???\n");
//             break;
//           }
//         }
//       MPI_Init_thread_DONE=1;
//       }
// //     else {
// //       printf("%s : %s initialisation ; MPI already initialised \n",__FUNCTION__, solver_name);
// //       }
//     slv->parameters = (PASTIX_STRUC_C *) calloc(1,sizeof(PASTIX_STRUC_C));
//     /*id_pastix = (PASTIX_STRUC_C *) calloc(1,sizeof(PASTIX_STRUC_C));*/
//     
//     id_pastix = slv->parameters;
//     /* --------------- FIN INIT PASTIX SOLVER ------------------ */
// #else 
//     printf("Please compile poc-solvers library with -DPASTIX \n");
//     exit(-1);
// #endif
//   } 
// 
// /** --------------- INIT PETSC SOLVER ------------------ */
//   else if (strcmp(solver_name,"PETSC") == 0){
//     printf("PETSC pas encore implimente\n");
//     exit(-1);
// /* --------------- FIN INIT PETSC SOLVER ------------------ */
//     }
// 
// /** --------------- INIT LAPACK SOLVER ------------------ */
//   else if (strcmp(solver_name,"LAPACK") == 0){
//     slv->parameters = NULL;
// /* --------------- FIN INIT LAPACK SOLVER ------------------ */
// #ifdef LAPACK
// #else
//     printf("Please compile poc-solvers library with -DLAPACK \n");
//     exit(-1);
// #endif
//     }
// 
// /** --------------- INIT SpDOMESTIC SOLVER ------------------ */
//   else if (strcmp(solver_name,"SpDOMESTIC") == 0){
//     slv->parameters = NULL;   
// /* --------------- FIN INIT LAPACK SOLVER ------------------ */
//     /*exit(-1);*/
//     }
// 
//   else {
//     printf("Connait pas ce solveur : %s \n",solver_name);
//     slv = NULL;
//     exit(-1);
//     }
// 
//   return(slv);
// }
// 
// /*************************************************************************//*
// /*                               FACTORIZE                               */
// /*************************************************************************/
// 
// /*************************************************************************/
// 
//   int factorizez(solver *slv,tripletz *Mat, int verbose) {
//   
// /*************************************************************************/
// 
//   char *solver_name;
//   int rtn,base,format;
//   int status;
//   double default_precision=1.e-08;
//   double cputime1, cputime2;
//   int proc_id, ierr;  
//   bool debug =false;
//  
// #ifdef  HAVE_MPI
//   if(MPI_Init_thread_DONE==1) {
//     ierr = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
//     cputime1 = MPI_Wtime();
//     }
// #endif     
//    
//   solver_name=slv->name;
//   
//   if(verbose==1) printf("%s complex-matrix factorization : %s \n",__FUNCTION__,solver_name);
// 
//   if (strcmp(solver_name,"MUMPS") == 0) {
//     base=1;
//     format=COO;
//     rtn = convertz(Mat, format, base, verbose);
//     if(rtn!=0) return(-1);
//     rtn = factorizez_mumps(slv, Mat, verbose, false);
//     }
//   else if  (strcmp(solver_name,"UMFPACK") == 0) {
//     base=0;
//     format=CSC;
//     status=convertz(Mat,format,base, verbose);
//     if(status!=0) return(-1);
//     rtn=factorizez_umfpack(slv, Mat, verbose);
//     }
//   else if  (strcmp(solver_name,"LAPACK") == 0) {
//     base=1;
//     /*convertz(Mat,3,base);*/
//     /*rtn=factorizez_lapack(slv, Mat);*/
//     }
//   else if  (strcmp(solver_name,"SpDOMESTIC") == 0) {
//     base=0;
//     format=CSC;
//     status=convertz(Mat,format,base, verbose);
//     if(status!=0) return(-1);
//     rtn=factorizez_spdomestic(slv, Mat);
//     }
//   else if (strcmp(solver_name,"HIPS") == 0){
//     base=1;
//     format=COO;
// //     format=COO_COL_ARRANGED;
//     rtn = convertz(Mat, format, base, verbose);
//     if(rtn!=0) return(-1);
// #ifdef __DEBUG_HIPS
//     debug=true;
// #endif
//     rtn=factorizez_hips(slv, Mat, verbose, debug);
//     }
//   else if (strcmp(solver_name,"PASTIX") == 0){
// // #ifdef  HAVE_MPI_HERE
// #ifdef  HAVE_MPI
//     // Attention en mode parallel la matrix arrive en FORMAT COO
//     // avec des indices globaux...
//     // la redistibution, conversion, ...
//     // Les jobs sont faits dans factorize_paxtix_parallel 
//     base=1;
// //     format=COO;
//     format=COO_COL_ARRANGED;
//     rtn = convertz(Mat, format, base, verbose);
//     rtn=factorizez_pastix_parallel(slv, Mat, verbose, debug);
// #else    
//      base=1; 
//      status=convertz(Mat,1,0, verbose);
//      if(status!=0) return(-1);
//      rtn=factorizez_pastix(slv, Mat, verbose);
// #endif     
//      }
//   else if (strcmp(solver_name,"PASTIX-SEQUENTIAL") == 0){
//      base=1; 
//      status=convertz(Mat,1,0, verbose);
//      if(status!=0) return(-1);
//      rtn=factorizez_pastix(slv, Mat, verbose);
//      }
//    else {
//      printf("factorizez, unknown solver : %s \n",solver_name);
//      rtn = -1;
//      }
//      
// #ifdef  HAVE_MPI
//   if(MPI_Init_thread_DONE==1 && verbose==1) {
//     cputime2 = MPI_Wtime();
//     if (verbose==1 && proc_id == 0) {
//       printf("--------------------------------------------\n");
//       printf("%s: %s factorization ellapsed time= %e\n",__FUNCTION__, solver_name, cputime2-cputime1);
//       printf("--------------------------------------------\n"); 
//       }
//     }
// #endif
//    return(rtn);
//    
// }
// 
// /*************************************************************************/
// 
// int factorizez_spdomestic(solver *slv,tripletz *Mat) {
//   
// /*************************************************************************/
//   int status=0,n,order=3,bcl;
//   complex<double> tol=1.0; /* unsymetrique usage*/
//   cs_ci *A = cs_calloc (1, sizeof (cs_ci)) ;    /* allocate the cs struct */  ; 
//   cs_cis *S ;
//   cs_cin *N ;
//   
//   Mat->A = A;
//   A->nzmax = Mat->nnz;
//   A->m =  Mat->nrow;
//   A->n =  Mat->ncol;
//   A->p =  Mat->i;
//   A->i =  Mat->j;
//   A->x =  Mat->x;
//   A->nz=-1;
//   n = A->n;
// //  printf("factorizez_spdomestic : compute symbolics \n");
//   S = cs_ci_sqr (order, A, 0) ;                    /* ordering and symbolic analysis */
// //  printf("factorizez_spdomestic : LU factorization \n");
//   N = cs_ci_lu (A, S, tol) ;                    /* numeric LU factorization */
// //  printf("factorizez_spdomestic : done \n");
//   bcl = N->L->nz;
//   Mat->S = S;
//   Mat->N = N;
//   slv->Mat = Mat;
//   return(status);
// }
// 
// /*************************************************************************/
// 
// int factorizez_umfpack(solver *slv,tripletz *Mat, int verbose)
// 
// /*************************************************************************/
// {
//   int status;
// #ifdef UMFPACK
//   UMFPACK_STRUC_C *id;
//   char *solver_name;
//   double *Info, *Control;
//   complex<double> *Ax;
//   double *Az, *Ar;
//   int *Ap, *Ai,n,nnz,bcl;
//   void *Symbolic, *Numeric ;
//   
//   solver_name=slv->name;
//   if (strcmp(solver_name,"UMFPACK") == 0) {
//     if(verbose==1) printf("factorize UMFPACK ...\n");
//     id = slv->parameters;
//     slv->Mat = Mat;
//     Symbolic=id->Symbolic;
//     Numeric=id->Numeric;
//     Info=id->Info;
//     Control=id->Control;
//     Ap=Mat->i;
//     Ai=Mat->j;
//     Ax=Mat->x;
//     n=Mat->nrow;
//     nnz=Mat->nnz;
//     
//     Az = (double *) malloc (nnz * sizeof (double)) ;
//     Ar = (double *) malloc (nnz * sizeof (double)) ;
//     for(bcl=0;bcl<nnz;bcl++) {
//       Ar[bcl] = creal(Ax[bcl]);
//       Az[bcl] = cimag(Ax[bcl]);
//       }
// 
//     status = umfpack_zi_symbolic (n, n, Ap, Ai, Ar, Az, &Symbolic, Control, Info) ; 
// 
//     if (status < 0) {
//       printf("umfpack_zi_symbolic status=%d\n",status);
//       umfpack_zi_report_info (Control, Info) ;
//       umfpack_zi_report_status (Control, status) ;
//       return(-1);
//       /*error ("umfpack_di_symbolic failed") ;*/
//       }
// 
//     status = umfpack_zi_numeric (Ap, Ai, Ar, Az, Symbolic, &Numeric, Control, Info) ;
//     id->Symbolic=Symbolic;
// 
//     id->Numeric=Numeric; 
//     if (status < 0) {
//       printf("umfpack_zi_numeric status=%d\n",status);
//       umfpack_zi_report_info (Control, Info) ;
//       umfpack_zi_report_status (Control, status) ;
//       return(-1);
//       /*error ("umfpack_di_numeric failed") ;*/
//       }
//     free(Ar);
//     free(Az);
//     return(status);
//   }
// #else
//   status=-1;
//   printf("Please compile poc-solvers library with -DUMFPACK \n");
//   return(status);
// #endif
//   return(status); 
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int factorizez_hips_cyril(solver *slv,tripletz *Mat, int verbose, bool debug)
//   
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status,ierr;
// #ifdef HIPS
// #ifdef HAVE_MPI
//   zHIPS_STRUC_C *id;
//   INTS  sym_pattern, sym_matrix;
//   INTS id_hips, idnbr, i, j;
//   INTS *unknownlist, *hipsnodelist;
//   COEFZ *x, *rhsloc;
//   INTS proc_id, n, ln;
//   INTL *ia, nnz;
//   INTS *ja;
//   COEFZ *a;
//   INTS domsize, nproc;
//   INTS pbegin, pend;
//   INTS *mapptr,  *mapptr2,*mapp,*mapp2;
//   INTS p;
// 
//   int nglob, nnzglob, tmp, bcl, idi, info;
//   FILE *out;
//   char filename[50];
//   complex_t *x1,*x2, tol;
//   
//   MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
//   MPI_Comm_size(MPI_COMM_WORLD, &nproc);
//   
//   id = slv->parameters; 
//   
//   /* nglob */
//   tmp=Mat->nrow;
//   
//   if (debug && proc_id == 0) {printf("cpu=%d,  %s nloc=%d\n",proc_id,__FUNCTION__,tmp);}
//   
//   MPI_Allreduce( &tmp, &nglob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
//   
//   if (debug && proc_id == 0) printf("cpu=%d,  %s nglob=%d\n",proc_id,__FUNCTION__,nglob);
// 
//   idnbr = 100; /* total */
//   if (HIPSSOLVERNBZ == 0) {
//     ierr = zHIPS_Initialize(idnbr);
//     zHIPS_ExitOnError(ierr);
//     }
//   id_hips = HIPSSOLVERNBZ; /** id_hips of the linear system **/
//   HIPSSOLVERNBZ++;
//   
//   /** parameter domsize is an argument of testHIPS.ex **/
//   if (debug) {
//     printf("cpu=%d,  %s solver numero=%d\n",proc_id,__FUNCTION__,id_hips);
//     }
// 
//   domsize = nglob/nproc/2;
// 
//   if (debug && proc_id == 0) {printf("cpu=%d,  %s domsize=%d\n",proc_id,__FUNCTION__,domsize);}
// 
// //   zHIPS_SetDefaultOptions (id_hips, zHIPS_HYBRID);
//   zHIPS_SetDefaultOptions (id_hips, zHIPS_ITERATIVE);
//   
//   zHIPS_SetOptionINT(id_hips, zHIPS_PARTITION_TYPE, 0); 
//   zHIPS_SetOptionINT(id_hips, zHIPS_DOMSIZE, domsize); 
// /*    OPTIONS */
// 
//   tol=1e-12;
//   zHIPS_SetOptionREAL(id_hips,zHIPS_DROPTOL0,tol);
//   zHIPS_SetOptionREAL(id_hips,zHIPS_DROPTOL1,tol);
//   
//   if(verbose==1)  info=4;
//   else            info=0;
// 
//   zHIPS_SetOptionINT (id_hips, zHIPS_VERBOSE, info);
// //     zHIPS_SetOptionREAL(id_hips, INTS number, REAL value);
//   sym_matrix=0;
//   zHIPS_SetOptionINT (id_hips, zHIPS_SYMMETRIC, sym_matrix);
//   
//   /** C : numbering starts from 0 **/
//   /* but here starts from 1 **/
//   
//   zHIPS_SetOptionINT(id_hips, zHIPS_FORTRAN_NUMBERING, 0);
//   
// /***************************************************/
// /*            ENTER THE GRAPH                      */
// /***************************************************/
//   ierr = zHIPS_GraphBegin(id_hips,nglob, Mat->nnz);
//   zHIPS_ExitOnError(ierr);
//   for (bcl=0;bcl<Mat->nnz;bcl++) {
//     ierr = zHIPS_GraphEdge(id_hips, Mat->i[bcl]-1, Mat->j[bcl]-1);
//     }
//   ierr = zHIPS_GraphEnd(id_hips);
//   zHIPS_ExitOnError(ierr);
//   status = MPI_Barrier(MPI_COMM_WORLD);
//   
// // #ifdef __DEBUG_HIPS
// //   printf("%d,  factorize_hips GraphEnd OK...\n",proc_id);
// // #endif
//   
// /***************************************************/
// /*            ENTER A USER PARTITION               */
// /***************************************************/
// /* Only the master processor (0 by default) needs to enter the partition
//    collect inforamtion on proc 0 */
//   mapptr = (INTS *)malloc(sizeof(INTS)*(nproc+1));
//   mapptr2= (INTS *)malloc(sizeof(INTS)*(nproc+1));
//   
//   mapp =   (INTS *)malloc(sizeof(INTS)*nglob);
//   mapp2=   (INTS *)malloc(sizeof(INTS)*nglob);
// 
//   for (bcl=0;bcl<nproc+1;bcl++)   mapptr2[bcl]=0;
//   for (bcl=0;bcl<nglob;bcl++)   mapp2[bcl]  =0;
//   for (bcl=0;bcl<nproc+1;bcl++)   mapptr[bcl] =0;
//   for (bcl=0;bcl<nglob;bcl++)   mapp[bcl]   =0;
//   
//   for (bcl=proc_id+1;bcl<nproc+1;bcl++) mapptr2[bcl]=Mat->nrow;
// 
//   MPI_Allreduce( mapptr2, mapptr, nproc+1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
//       
// // #ifdef __DEBUG_HIPS
// //       for (bcl=0;bcl<nproc+1;bcl++) printf("...%d, ",mapptr[bcl]);
// //       printf("\n");
// // #endif
//  
//   idi=0;
//   mapp2[idi+mapptr[proc_id]] =  Mat->i[0]-1;
//   for (bcl=1; bcl < Mat->nnz; bcl++) {
//     if ( mapp2[idi+mapptr[proc_id]] != (Mat->i[bcl]-1) ) {
//       idi++;
//       /*printf("bcl=%d, idi=%d, [idi+mapptr[proc_id]=%d nglob=%d \n",bcl,idi,idi+mapptr[proc_id], nglob);*/
//       mapp2[idi+mapptr[proc_id]] = Mat->i[bcl]-1;    
//       }
//     }
//     
// //       sprintf(filename,"solver_fatorize_hips.1.%d.txt",proc_id);
// //       out=fopen(filename,"w");
// //       for (bcl=0;bcl<nglob;bcl++)   fprintf(out, "%d \n",mapp2[bcl]);
// //       fclose(out);
// 
//   MPI_Allreduce( mapp2, mapp, nglob, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
//   
// //       sprintf(filename,"solver_fatorize_hips.2.%d.txt",proc_id);
// //       out=fopen(filename,"w");
// //       for (bcl=0;bcl<nglob;bcl++)   fprintf(out, "%d \n",mapp[bcl]);
// //       fclose(out);
// 
//   /* Set the unknownlist */
//   unknownlist= (INTS *)malloc(sizeof(INTS)*(Mat->nrow));
//   unknownlist[0]=Mat->i[0]-1;
//   idi=0;
//   for (bcl=1; bcl < Mat->nnz; bcl++) {
//     if ( unknownlist[idi] != (Mat->i[bcl]-1) ) {
//       idi++;
//       unknownlist[idi] = Mat->i[bcl]-1;    
//       }
//     }
//   
// //       sprintf(filename,"solver_fatorize_hips.unknonwlist.%d.txt",proc_id);
// //       out=fopen(filename,"w");
// //       for (bcl=0; bcl < idi; bcl++)  fprintf(out, "%d \n",unknownlist[bcl]);
// //       fclose(out);
// 
//    /* Enter Partition */
//   status = MPI_Barrier(MPI_COMM_WORLD);
//   if (proc_id ==0) {
//     ierr = zHIPS_SetPartition(id_hips, nproc, mapptr, mapp);
//     zHIPS_ExitOnError(ierr);
//     }
//       
// //       status = MPI_Barrier(MPI_COMM_WORLD);
// //       printf("cpu=%d,  %s SetPartition OK...\n",proc_id,__FUNCTION__);
// 
//   free(mapptr2);
//   free(mapptr);
//   free(mapp2);
//   free(mapp);
//   
// /***************************************************/
// /*            GET THE LOCAL UNKNOWN LIST           */
// /***************************************************/
// //   ierr = zHIPS_GetLocalUnknownNbr(id_hips, &ln);
// //   zHIPS_ExitOnError(ierr);
// //   printf("cpu=%d, %s nombre de noeuds interne %d\n",proc_id,__FUNCTION__,ln);
// //   hipsnodelist = (INTS *)malloc(sizeof(INTS)*ln);
// //   ierr = zHIPS_GetLocalUnknownList(id_hips, hipsnodelist);
// //   zHIPS_ExitOnError(ierr);
//   
// /***************************************************/
// /*          ENTER THE MATRIX COEFFICIENT           */
// /***************************************************/
// /* The processor enter the rows pbegin to pend-1 (naive partition) 
//    In the zHIPS_ASSEMBLY_FOOL mode any processor can enter any coefficient :
//    nevertheless it is better to enter as much as possible coefficients
//    that are in the local matrix of the processor */
//   
// //   sprintf(filename,"solver_fatorize_hips.Matrice.%d.txt",proc_id);
// //   out=fopen(filename,"w");
// //   for(bcl=0;bcl<Mat->nnz;bcl++) {
// //     fprintf(out, "%d : %d,%d,%lf  \n",bcl,Mat->i[bcl], Mat->j[bcl], Mat->x[bcl]);
// //     }
// //   fclose(out);
// 
//   ierr = zHIPS_AssemblyBegin(id_hips, Mat->nnz, zHIPS_ASSEMBLY_OVW, zHIPS_ASSEMBLY_OVW, zHIPS_ASSEMBLY_FOOL, sym_matrix);
//   zHIPS_ExitOnError(ierr);
//   
// // #ifdef __DEBUG_HIPS
// //   printf("%d,  factorize_hips AssemblyBegin OK...\n",proc_id);
// // #endif    
// 
//   for(bcl=0;bcl<Mat->nnz;bcl++) {
//     ierr = zHIPS_AssemblySetValue(id_hips, Mat->i[bcl]-1, Mat->j[bcl]-1, Mat->x[bcl]);
//     zHIPS_ExitOnError(ierr);
//     }
//   ierr = zHIPS_AssemblyEnd(id_hips);
//   zHIPS_ExitOnError(ierr);
// 
// //   sprintf(filename,"solver_fatorize_hips.4.%d.txt",proc_id);
// //   out=fopen(filename,"w");
// //   fprintf(out, "%d  %d \n",idi, Mat->nrow);
// //   for (bcl=0; bcl < Mat->nrow; bcl++) fprintf(out, "%d => %d \n",bcl,unknownlist[bcl]);
// //   fclose(out);
//     
//   id->id_hips = id_hips;
//   id->unknownlist = unknownlist;
//   id->n = Mat->nrow;
//   id->nnz = Mat->nnz;
// 
//   if (debug && proc_id == 0) printf("cpu=%d, sortie %s \n",proc_id,__FUNCTION__);
// 
//   status = MPI_Barrier(MPI_COMM_WORLD);
// 
// //   if(debug) {
// //     /* une resolution pour test travail usr le local nodelist*/
// //     x1 = (complex_t *) malloc(sizeof(complex_t)*Mat->nrow);
// //     x2 = (complex_t *) malloc(sizeof(complex_t)*Mat->nrow);
// //     for (bcl=0;bcl < Mat->nrow;bcl++) x1[bcl]=x2[bcl]=1.0;
// //   
// //     ierr = zHIPS_MatrixVectorProduct (id_hips, Mat->x, x1);
// //   
// //     ierr = zHIPS_SetRHS(id_hips, Mat->nrow, unknownlist, x1, zHIPS_ASSEMBLY_OVW, zHIPS_ASSEMBLY_OVW, zHIPS_ASSEMBLY_FOOL);
// //     zHIPS_ExitOnError(ierr);
// //     ierr = zHIPS_GetSolution(id_hips,Mat->nrow,unknownlist, x2, zHIPS_ASSEMBLY_FOOL);
// //     }
//     
// //   sprintf(filename,"solver_fatorize_hips.RESO1.%d.txt",proc_id);
// //   out=fopen(filename,"w");
// //   fprintf(out, "%d  %d \n",idi, Mat->nrow);
// //   for (bcl=0; bcl < Mat->nrow; bcl++) fprintf(out, "%d => %d \n",bcl,x2[bcl]);
// //   fclose(out);
//   
//   ierr=0;
//   status=0;
//   
//   return(status);
//   
// #endif 
// #endif
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int factorizez_hips(solver *slv,tripletz *Mat, int verbose, bool debug)
//     
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status,ierr;
// #ifdef HIPS
// #ifdef HAVE_MPI
//   zHIPS_STRUC_C *id;
//   INTS  sym_pattern, sym_matrix;
//   INTS id_hips, idnbr, i, j;
//   INTS *unknownlist, *hipsnodelist;
//   COEFZ *x, *rhsloc;
//   INTS proc_id, n, ln;
//   INTL *ia, nnz;
//   INTS *ja;
//   COEFZ *a;
//   INTS domsize, nproc;
//   INTS pbegin, pend;
//   INTS *mapptr,  *mapptr2,*mapp,*mapp2;
//   INTS p;
// 
//   int nglob, nnzglob, tmp, bcl, idi, info;
//   FILE *out;
//   char filename[50];
//   complex_t *x1,*x2, tol;
//   
//   MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
//   MPI_Comm_size(MPI_COMM_WORLD, &nproc);
//   
//   id = slv->parameters; 
//   
//   /* nglob */
//   tmp=Mat->nrow;
//   
//   if (debug && proc_id == 0) {printf("cpu=%d,  %s nloc=%d\n",proc_id,__FUNCTION__,tmp);}
//   
//   MPI_Allreduce( &tmp, &nglob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
//   
//   if (debug && proc_id == 0) printf("cpu=%d,  %s nglob=%d\n",proc_id,__FUNCTION__,nglob);
// 
//   idnbr = 100; /* total */
//   if (HIPSSOLVERNBZ == 0) {
//     ierr = zHIPS_Initialize(idnbr);
//     zHIPS_ExitOnError(ierr);
//     }
//   id_hips = HIPSSOLVERNBZ; /** id_hips of the linear system **/
//   HIPSSOLVERNBZ++;
//   
//   status = MPI_Barrier(MPI_COMM_WORLD);
//    
//   if (debug && proc_id == 0) {
//     printf("cpu=%d,  %s solver numero=%d\n",proc_id,__FUNCTION__,id_hips);
//     }
//     
//   status = MPI_Barrier(MPI_COMM_WORLD);
// 
//   domsize = nglob/nproc/2;
// //   domsize = nglob/(nproc+1);
// 
//   if (debug && proc_id == 0) {
//     printf("cpu=%d,  %s domsize=%d\n",proc_id,__FUNCTION__,domsize);
//     }
//     
//   status = MPI_Barrier(MPI_COMM_WORLD);
// 
//   zHIPS_SetOptionINT(id_hips, zHIPS_PARTITION_TYPE, 0); 
//   zHIPS_SetOptionINT(id_hips, zHIPS_DOMSIZE, domsize);
//   
//   status = MPI_Barrier(MPI_COMM_WORLD);
//   
// //   zHIPS_SetDefaultOptions (id_hips, zHIPS_HYBRID);
//   zHIPS_SetDefaultOptions (id_hips, zHIPS_ITERATIVE);
// //   zHIPS_SetDefaultOptions (id_hips, zHIPS_DIRECT);
//   
//   status = MPI_Barrier(MPI_COMM_WORLD);
//   
//   tol=1e-12;
//   zHIPS_SetOptionREAL(id_hips,zHIPS_DROPTOL0,tol);
//   tol=1e-12;
//   zHIPS_SetOptionREAL(id_hips,zHIPS_DROPTOL1,tol);
// //   tol=1e-15;
// //   zHIPS_SetOptionREAL(id_hips,zHIPS_DROPTOLE,tol);
//    
// //   zHIPS_SetOptionINT(id_hips, zHIPS_ITMAX,150);
// //   zHIPS_SetOptionINT(id_hips, zHIPS_ITMAX_SCHUR,150);
//    
//   if(verbose==1)  info=4;
//   else            info=0;
// 
//   zHIPS_SetOptionINT (id_hips, zHIPS_VERBOSE, info);
// 
//   sym_matrix=0;
//   zHIPS_SetOptionINT (id_hips, zHIPS_SYMMETRIC, sym_matrix);
//   zHIPS_SetOptionINT(id_hips,zHIPS_GRAPH_SYM,0);
//   
//   /** C : numbering starts from 0 **/
//   zHIPS_SetOptionINT(id_hips, zHIPS_FORTRAN_NUMBERING, 0);
//   
// /***************************************************/
// /*            ENTER THE GRAPH                      */
// /***************************************************/
//   status = MPI_Barrier(MPI_COMM_WORLD);
//   
//   ierr = zHIPS_GraphBegin(id_hips,nglob, Mat->nnz);
//   zHIPS_ExitOnError(ierr);
//   for (bcl=0;bcl<Mat->nnz;bcl++) {
//     ierr = zHIPS_GraphEdge(id_hips, Mat->i[bcl]-1, Mat->j[bcl]-1);
//     }
//   ierr = zHIPS_GraphEnd(id_hips);
//   zHIPS_ExitOnError(ierr);
//   status = MPI_Barrier(MPI_COMM_WORLD);
//     
// /***************************************************/
// /*            ENTER A USER PARTITION               */
// /***************************************************/
// /* Only the master processor (0 by default) needs to enter the partition
//    collect inforamtion on proc 0 */
//   mapptr = (INTS *)malloc(sizeof(INTS)*(nproc+1));
//   mapptr2= (INTS *)malloc(sizeof(INTS)*(nproc+1));
//   
//   mapp =   (INTS *)malloc(sizeof(INTS)*(nglob+1));
//   mapp2=   (INTS *)malloc(sizeof(INTS)*(nglob+1));
// 
//   for (bcl=0;bcl<nproc+1;bcl++)   mapptr2[bcl]=0;
//   for (bcl=0;bcl<nglob+1;bcl++)   mapp2[bcl]  =0;
//   
//   for (bcl=0;bcl<nproc+1;bcl++)   mapptr[bcl] =0;
//   for (bcl=0;bcl<nglob+1;bcl++)   mapp[bcl]   =0;
//   
//   for (bcl=proc_id+1;bcl<nproc+1;bcl++) mapptr2[bcl]=Mat->nrow;
// 
//   MPI_Allreduce( mapptr2, mapptr, nproc+1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
//   
//   idi=0;
//   mapp2[idi+mapptr[proc_id]] =  Mat->i[0]-1;
//   for (bcl=1; bcl < Mat->nnz; bcl++) {
//     if ( mapp2[idi+mapptr[proc_id]] != (Mat->i[bcl]-1) ) {
//       idi++;
//       mapp2[idi+mapptr[proc_id]] = Mat->i[bcl]-1;    
//       }
//     }
//     
// //   status = MPI_Barrier(MPI_COMM_WORLD);
// //   
// //   printf("cpu=%d ; ndofs %d %d\n", proc_id, idi+1,Mat->nrow);
//   
//   MPI_Allreduce( mapp2, mapp, nglob, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
//   
//   /* Set the unknownlist */
//   unknownlist= (INTS *)malloc(sizeof(INTS)*(Mat->nrow));
//   unknownlist[0]=Mat->i[0]-1;
//   idi=0;
//   for (bcl=1; bcl < Mat->nnz; bcl++) {
//     if ( unknownlist[idi] != (Mat->i[bcl]-1) ) {
//       idi++;
//       unknownlist[idi] = Mat->i[bcl]-1;    
//       }
//     }
//   
//   status = MPI_Barrier(MPI_COMM_WORLD);
//  
//    /* Enter Partition */
//   status = MPI_Barrier(MPI_COMM_WORLD);
//   if (proc_id ==0) {
//     ierr = zHIPS_SetPartition(id_hips, nproc, mapptr, mapp);
//     zHIPS_ExitOnError(ierr);
//     }
//   status = MPI_Barrier(MPI_COMM_WORLD);
//     
//   free(mapptr2);
//   free(mapptr);
//   free(mapp2);
//   free(mapp);
//   
// /***************************************************/
// /*            GET THE LOCAL UNKNOWN LIST           */
// /***************************************************/
// //   ierr = zHIPS_GetLocalUnknownNbr(id_hips, &ln);
// //   zHIPS_ExitOnError(ierr);
// //   printf("cpu=%d, %s nombre de noeuds interne %d\n",proc_id,__FUNCTION__,ln);
// //   hipsnodelist = (INTS *)malloc(sizeof(INTS)*ln);
// //   ierr = zHIPS_GetLocalUnknownList(id_hips, hipsnodelist);
// //   zHIPS_ExitOnError(ierr);
//   
// /***************************************************/
// /*          ENTER THE MATRIX COEFFICIENT           */
// /***************************************************/
// /* The processor enter the rows pbegin to pend-1 (naive partition) 
//    In the zHIPS_ASSEMBLY_FOOL mode any processor can enter any coefficient :
//    nevertheless it is better to enter as much as possible coefficients
//    that are in the local matrix of the processor */
//   
//     
//   status = MPI_Barrier(MPI_COMM_WORLD);
//   
//   if (debug) {
//     printf("cpu=%d,  %s now begin assembly\n",proc_id,__FUNCTION__);
//     }
//   status = MPI_Barrier(MPI_COMM_WORLD);
//   
//   ierr = zHIPS_AssemblyBegin(id_hips, Mat->nnz, zHIPS_ASSEMBLY_OVW, zHIPS_ASSEMBLY_OVW, zHIPS_ASSEMBLY_FOOL, sym_matrix);
//   zHIPS_ExitOnError(ierr);
//   
//   for(bcl=0;bcl<Mat->nnz;bcl++) {
//     ierr = zHIPS_AssemblySetValue(id_hips, Mat->i[bcl]-1, Mat->j[bcl]-1, Mat->x[bcl]);
//     zHIPS_ExitOnError(ierr);
//     }
//   ierr = zHIPS_AssemblyEnd(id_hips);
//   zHIPS_ExitOnError(ierr);
//     
//   status = MPI_Barrier(MPI_COMM_WORLD);
//   
//   if (debug) {
//     printf("cpu=%d,  %s now set assembly done\n",proc_id,__FUNCTION__);
//     }
//   status = MPI_Barrier(MPI_COMM_WORLD);
//   
// 
//   id->id_hips = id_hips;
//   id->unknownlist = unknownlist;
//   id->n = Mat->nrow;
//   id->nnz = Mat->nnz;
// 
//   if (debug && proc_id == 0) printf("cpu=%d, sortie %s \n",proc_id,__FUNCTION__);
// 
//   status = MPI_Barrier(MPI_COMM_WORLD);
// 
//   ierr=0;
//   status=0;
//   
//   return(status);
//   
// #endif 
// #endif
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int factorizez_pastix(solver *slv,tripletz *Mat, int verbose)
//   
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
// #ifdef PASTIX
//   PASTIX_STRUC_C *id;
//   char *solver_name;
//   complex<double> *Ax;
//   pastix_int_t *Ap, *Ai;
//   int  n,nnz,nbthread;
//   complex<double> *rhs = NULL; /* right hand side                      */
//   int    verbosemode;         /* Level of verbose mode (0, 1, 2)      */
//   char   *type        = NULL; /* type of the matrix                   */
//   int    bcl, bcl2;
//   int    flagsym=0;
//   pastix_int_t ncol;          /* Size of the matrix                   */
//  
//   if(verbose==1) printf("enter %s\n",__FUNCTION__);
//   
//   solver_name=slv->name;
//   
// //   if (strcmp(solver_name,"PASTIX") == 0) {
//   id = slv->parameters;
//   slv->Mat = Mat; 
//   Ap=Mat->i;
//   Ai=Mat->j;
//   Ax=Mat->x;
//   n=Mat->nrow;
//   nnz=Mat->nnz;
// 
// //   pastix_data = id->pastix_data;
// //   perm = id->perm;
// //   invp = id->invp;
// //   iparm = id->iparm;
// //   dparm = id->dparm;
// 
//   
// //   pastix_data=id->pastix_data;
// 
// /*******************************************/
// /* Initialize parameters to default values */
// /*******************************************/
// 
//   verbosemode=verbose;  
//   
//   z_pastix_checkMatrix(MPI_COMM_WORLD, verbosemode,
//                        API_SYM_NO, API_YES, n, &Ap, &Ai, &Ax, NULL,1);
// 
//   Mat->x=Ax;
//   Mat->i=Ap;
//   Mat->j=Ai;
//   
//   id->iparm[IPARM_MODIFY_PARAMETER] = API_NO;
//   z_pastix(&(id->pastix_data), MPI_COMM_WORLD,
//          n, Ap, Ai, Ax, id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
//   
// /*******************************************/
// /*       Customize some parameters         */
// /*******************************************/
// #ifdef OMP_H
//   nbthread = omp_get_max_threads(); 
//   if(verbose==1) printf("#using %d threads (omp_get_max_threads)...................... \n",nbthread);
// #else
//   nbthread = 1;  
//   if(verbose==1) printf("#using only 1 thread!...................... \n",nbthread);
// #endif
//   id->iparm[IPARM_SYM]           = API_SYM_NO;
//   id->iparm[IPARM_FACTORIZATION] = API_FACT_LU;
//   id->iparm[IPARM_THREAD_NBR] = nbthread;
//   id->iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
//   id->iparm[IPARM_LEVEL_OF_FILL] = 0;
//   id->iparm[IPARM_START_TASK] = API_TASK_ORDERING;
//   id->iparm[IPARM_END_TASK]   = API_TASK_NUMFACT;
// //   id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-8;
//   id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-12;
//   id->iparm[IPARM_VERBOSE]             = verbosemode;
//   /*iparm[IPARM_END_TASK]   = API_TASK_REFINE;*/
// 
//   id->perm = malloc(n*sizeof(pastix_int_t));
//   id->invp = malloc(n*sizeof(pastix_int_t));
//   
// /*Factorisation.........*/ 
//   z_pastix(&(id->pastix_data), MPI_COMM_WORLD,
//          n, Ap, Ai, Ax,
//          id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
//   
//   if (status < 0)  {
//     printf("%s factorization status=%d\n",__FUNCTION__,status);  
//     return(-1);
//     }
// 
// /*   Mat->i = Ap; */
// /*   Mat->j = Ai; */
// /*   Mat->x = Ax; */
// /*   slv->Mat = Mat; */
// 
// /*  /\*TEST SOLVE *\/ */
// /*   Mat=slv->Mat; */
// /*     Ap=Mat->i; */
// /*     Ai=Mat->j; */
// /*     Ax=Mat->x; */
// /*     n=Mat->nrow; */
// /*     nnz=Mat->nnz; */
//  
// /*    printf("Première resolution......\n"); */
// /*   id->iparm[IPARM_START_TASK] = API_TASK_SOLVE; */
// /*   id->iparm[IPARM_END_TASK]   = API_TASK_REFINE; */
// /*   printf("iparm dans facto\n"); */
// /*   /\*for (bcl=0;bcl<IPARM_SIZE;bcl++) printf("iparm[%d] = %d, ",bcl,iparm[bcl]);*\/ */
// /*   printf(" solve dans facto 1\n"); */
// /*   d_pastix(&(id->pastix_data), MPI_COMM_WORLD, */
// /*          n, Ap, Ai, Ax, */
// /*          id->perm, id->invp, rhs, 1, id->iparm, id->dparm);  */
// /*   printf(" solve dans facto 2\n"); */
// /*   id->iparm[IPARM_START_TASK] = API_TASK_SOLVE; */
// /*   id->iparm[IPARM_END_TASK]   = API_TASK_REFINE; */
// /*   d_pastix(&(id->pastix_data), MPI_COMM_WORLD, */
// /*          n, Ap, Ai, Ax, */
// /*          id->perm, id->invp, rhs, 1, id->iparm, id->dparm);  */
// 
// /*   /\*printf("Adresse pastix_data %d \n",&pastix_data);*\/ */
// /*   printf("Adresse pastix_data %d \n",&(id->pastix_data)); */
// /*   printf("Adresse id %d \n",id); */
// /*   printf("Adresse solver-parameter %d \n",slv->parameters); */
//  
//   if(verbose==1) printf("Factorizez pastix done...\n");
//   status=0;        
//   return(status);
//   
// #else
//   status=-1;
//   printf("Please compile poc-solvers library with -DPASTIX \n");
//   return(status);
// #endif 
//   return(status); 
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int factorizez_pastix_parallel(solver *slv,tripletz *Mat, int verbose, bool debug) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   int status;
// #ifdef PASTIX
// #ifdef HAVE_MPI
//   PASTIX_STRUC_C *id;
//   char *solver_name;
//   complex<double> *Ax;
//   pastix_int_t *Ap, *Ai;
//   int  n,nnz,nbthread;
//   complex<double> *rhs = NULL; /* right hand side                      */
//   int    verbosemode;         /* Level of verbose mode (0, 1, 2)      */
//   char   *type        = NULL; /* type of the matrix                   */
//   int    bcl, bcl2;
//   int    flagsym=0;
//   pastix_int_t    ncol;       /* Size of the matrix                   */
//   // Global variables
//   int ncolglob, nnzglob;
//   int *Apglob, *Aiglob, *colptr, *rows;
//   complex<double> *Axglob, *values;
//   char filemat[50];
//   FILE *f_out2;
//   int   mpid,sz;
// 
//   MPI_Barrier(MPI_COMM_WORLD);
// 
//   if(verbose) printf("enter %s\n",__FUNCTION__);
//   
//   solver_name=slv->name;
//   
// //   if (strcmp(solver_name,"PASTIX") == 0) {
//     
//   id = slv->parameters;
//   slv->Mat = Mat; 
//   Ap=Mat->i;
//   Ai=Mat->j;
//   Ax=Mat->x;
//   n=Mat->nrow;
//   nnz=Mat->nnz;
// 
//   MPI_Barrier(MPI_COMM_WORLD);
// 
//   MPI_Comm_rank(MPI_COMM_WORLD, &mpid);
//   MPI_Comm_size(MPI_COMM_WORLD, &sz);
//   
// /*------------------------------------------------------------------------------
//   Dispatch Matrix global matrix */
//   z_coo2glob(n, nnz, Ap, Ai, Ax, &ncolglob, &nnzglob, &Apglob, &Aiglob, &Axglob);
//   
//   free(Ap);
//   free(Ai);
//   free(Ax);
//   
//   MPI_Barrier(MPI_COMM_WORLD);
//   
// /*------------------------------------------------------------------------------
//   convert structure to CSC */
//   z_globcoo2csc(ncolglob, nnzglob, Apglob, Aiglob, Axglob, &colptr, &rows, &values);
//   
//   free(Apglob);
//   free(Aiglob);
//   free(Axglob);
//   
//   MPI_Barrier(MPI_COMM_WORLD);
// 
//   if(debug==1) {
// /*------------------------------------------------------------------------------
//     PRINT out */
//     sprintf(filemat,"Solverz_MatrixGlobal_CSC.%d.%d.txt",nbsyslinearz2,mpid);
//     MPI_Barrier(MPI_COMM_WORLD);
//     printf("%d::: z_globcoo2csc filemat=%s \n",mpid, filemat);
//     f_out2 = fopen(filemat,"w");
//     fprintf(f_out2,"  %%%MatrixMarket matrix coordinate complex general \n");
//     fprintf(f_out2,"%d  %d  %d 1 0 \n",ncolglob,ncolglob,nnzglob);
//     for (bcl=0; bcl < ncolglob; bcl++) {
// //       fprintf(f_out,"%d %d\n",bcl+1,colptr[bcl]); // Use Fortran Numbering 
//       for (bcl2=colptr[bcl]-1;bcl2<colptr[bcl+1]-1;bcl2++) {
//         fprintf(f_out2,"%d  %d  %e %e \n",bcl+1,rows[bcl2],values[bcl2]);
//         }
//       }
//     fclose(f_out2); 
//     }
// 
// //     free(Ap);
// //     free(Ai);
// //     free(Ax);
//     
//   Ap = colptr;
//   Ai = rows;
//   Ax = values;
//   
//   Mat->nrow = ncolglob;
//   Mat->nnz  = nnzglob;
//    
// /*******************************************/
// /* Initialize parameters to default values */
// /*******************************************/
// 
//   verbosemode=verbose;  
//   
//   z_pastix_checkMatrix(MPI_COMM_WORLD, verbosemode,
//                      API_SYM_NO, API_YES,
//                      ncolglob, &Ap, &Ai, &Ax, NULL,1);
//   Mat->i=Ap;
//   Mat->j=Ai;
//   Mat->x=Ax;
//   Mat->comptype=CSC;
//   
//   id->iparm[IPARM_MODIFY_PARAMETER] = API_NO;
//   z_pastix(&(id->pastix_data), MPI_COMM_WORLD,
//          ncolglob, Ap, Ai, Ax,
//          id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
// /*******************************************/
// /*       Customize some parameters         */
// /*******************************************/
// #ifdef OMP_H
//   nbthread = omp_get_max_threads(); 
//   if(verbose) printf("#using %d threads (omp_get_max_threads)...................... \n",nbthread);
// #else
//   nbthread = 1;  
//   if(verbose) printf("#using only 1 thread!...................... \n",nbthread);
// #endif
//   verbosemode=verbose;    
//   id->iparm[IPARM_SYM]           = API_SYM_NO;
//   id->iparm[IPARM_FACTORIZATION] = API_FACT_LU;
//   id->iparm[IPARM_THREAD_NBR] = nbthread;
//   id->iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
//   id->iparm[IPARM_LEVEL_OF_FILL] = 0;
//   id->iparm[IPARM_START_TASK] = API_TASK_ORDERING;
//   id->iparm[IPARM_END_TASK]   = API_TASK_NUMFACT;
// //   id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-8;
//   id->dparm[DPARM_EPSILON_REFINEMENT] = 1e-12;
//   id->iparm[IPARM_VERBOSE]             = verbosemode;
// 
//   /*iparm[IPARM_END_TASK]   = API_TASK_REFINE;*/
// 
// //   id->perm = malloc(ncolglob*sizeof(pastix_int_t));
// //   id->invp = malloc(ncolglob*sizeof(pastix_int_t));
//   id->perm = new pastix_int_t[ncolglob];
//   id->invp = new pastix_int_t[ncolglob];
//   
//   /*rhs = malloc(n*sizeof(double));*/
//   
// /*Factorisation.........*/ 
//   if(verbose==1) printf("%s pastix factorisation start\n",__FUNCTION__);
//   z_pastix(&(id->pastix_data), MPI_COMM_WORLD,
//          ncolglob, Ap, Ai, Ax,
//          id->perm, id->invp, rhs, 1, id->iparm, id->dparm);
//   if(verbose==1) printf("%s pastix factorisation end\n",__FUNCTION__);
//   
//   status=id->iparm[IPARM_ERROR_NUMBER];
// 
//   if (status < 0)  { 
//     printf("z_pastix factorization status=%d\n",status);  
//     return(-1);
//     }
//  
//   status=0;        
//   return(status);
// #endif 
// #else
//   status=-1;
//   printf("Please compile poc-solvers library with -DPASTIX \n");
//   return(status);
// #endif 
//   return(status); 
// }
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// 
// /*************************************************************************/
// /*=======================================================================*/
// /*====================== FIN FACTORISE ==================================*/
// /*=======================================================================*/
// 
// 
// 
// /*=======================================================================*/
// /*                              SOLVE                                    */
// /*************************************************************************/
// /*=======================================================================*/
// /*=======================================================================*/
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int solvez(solver *slv,complex<double> *RHS, int transp, int verbose) {
//   
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   char *solver_name;
//   int rtn;
//   double cputime1, cputime2;
//   int ierr, proc_id;
//   solver_name=slv->name;
//   bool debug=false; 
// 
// #ifdef  HAVE_MPI
//   if(MPI_Init_thread_DONE==1) {
//     ierr = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
//     cputime1 = MPI_Wtime();
//     }
// #endif
//    
//   if(verbose) printf("%s complex-system solver, transpose=%d \n",solver_name, transp);
//    
//   if (strcmp(solver_name,"MUMPS") == 0) { 
//     rtn = solvez_mumps(slv, RHS,transp, verbose, debug);
// //     printf("solver not implemented for complex matrices : %s \n",solver_name);
//     }
//   else if  (strcmp(solver_name,"UMFPACK") == 0) {
//     rtn=solvez_umfpack(slv, RHS,transp);
//     }
//   else if  (strcmp(solver_name,"LAPACK") == 0) {
//     /*rtn=solvez_lapack(slv, RHS,transp);*/
//     printf("solver not implemented for complex matrices : %s \n",solver_name);
//     }
//   else if  (strcmp(solver_name,"SpDOMESTIC") == 0) {
//     rtn=solvez_spdomestic(slv, RHS,transp);
//     }
//   else if  (strcmp(solver_name,"HIPS") == 0) {
// #ifdef __DEBUG_HIPS
//     debug=true;
// #endif
//     rtn=solvez_hips(slv, RHS,transp, verbose, debug);
//     }
//   else if  (strcmp(solver_name,"PASTIX") == 0) {
// #ifdef HAVE_MPI    
//     MPI_Barrier(MPI_COMM_WORLD);
//     rtn=solvez_pastix_parallel(slv, RHS,transp, verbose, debug);
// #else
//     rtn=solvez_pastix_sequential(slv, RHS,transp, verbose, debug);
// #endif    
//     }
//   else if  (strcmp(solver_name,"PASTIX-SEQUENTIAL") == 0) {
//     rtn=solvez_pastix_sequential(slv, RHS,transp, verbose, debug);
//     }
//   else {
//     printf("solvez, unknown solver : %s \n",solver_name);
//     }
//      
// #ifdef  HAVE_MPI
//   if(MPI_Init_thread_DONE==1 && verbose==1) {
// //   if(MPI_Init_thread_DONE==1) {
//     cputime2 = MPI_Wtime();
//     if (proc_id == 0 && verbose==1) {
//       printf("--------------------------------------------\n");
//       printf("%s: %s solving ellapsed time= %e\n", __FUNCTION__, solver_name, cputime2-cputime1);
//       printf("--------------------------------------------\n"); 
//       }
//     }
// #endif
// 
//    return(rtn);
//   
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// int solvez_umfpack(solver *slv,complex<double> *RHS,int transp) {
//   
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   int status,bcl;
//   status=-1;
// #ifdef UMFPACK
//   UMFPACK_STRUC_C *id;
//   char *solver_name;
// //  double *Info, *Control, *Ax, *tmpr, *tmpz;
//   double *Info, *Control, *tmpr, *tmpz;
//   complex<double> *Ax;
//   complex<double> *chk;
//   double *Az, *Ar, *rhsr, *rhsz;
//   int *Ap, *Ai,n,nnz;
//   void *Symbolic, *Numeric ;
//   tripletz *Mat;
//   int INCXY;
//   int verbose=0;
//  
//   status = 0;
//   solver_name=slv->name;
//   if (strcmp(solver_name,"UMFPACK") == 0) {
//     id = slv->parameters;
//     Symbolic=id->Symbolic;
//     Numeric=id->Numeric;
//     Info=id->Info;
//     Control=id->Control;
//     Mat=slv->Mat;
//     Ap=Mat->i;
//     Ai=Mat->j;
//     Ax=Mat->x;
//     if(verbose==1) printf("complex size = %d\n",sizeof(complex<double>));
//     n=Mat->nrow;
//     nnz=Mat->nnz;
//     tmpr = (double *) malloc (n * sizeof (double)) ;
//     tmpz = (double *) malloc (n * sizeof (double)) ;
//     rhsr = (double *) malloc (n * sizeof (double)) ;
//     rhsz = (double *) malloc (n * sizeof (double)) ;
//     Az = (double *) malloc (nnz * sizeof (double)) ;
//     Ar = (double *) malloc (nnz * sizeof (double)) ;
//     for(bcl=0;bcl<nnz;bcl++) {
//       Ar[bcl] = creal(Ax[bcl]);
//       Az[bcl] = cimag(Ax[bcl]);
//       }
//     for(bcl=0;bcl<n;bcl++) {
//       rhsr[bcl] = creal(RHS[bcl]);
//       rhsz[bcl] = cimag(RHS[bcl]);
//       }
// 
//     if (transp == 1) 
//       status = umfpack_zi_solve (UMFPACK_At, Ap, Ai, Ar, Az, tmpr, tmpz, rhsr, rhsz, Numeric, Control, Info) ;
//     else 
// //      status = umfpack_zi_solve (UMFPACK_A, Ap,  Ai, Ax, Az, tmpr, tmpz, rhsr, rhsz, Numeric, Control, Info) ;
//       status = umfpack_zi_solve (UMFPACK_A, Ap,  Ai, Ar, Az, tmpr, tmpz, rhsr, rhsz, Numeric, Control, Info) ;
//     
//     for(bcl=0;bcl<n;bcl++)  RHS[bcl] = tmpr[bcl] + tmpz[bcl] * _Complex_I; 
//     
//     if (status < 0) {
//       umfpack_zi_report_info (Control, Info) ;
//       umfpack_zi_report_status (Control, status) ;
//         /*error ("umfpack_di_solve failed") ;*/
//       }
//     free(tmpr);
//     free(tmpz);
//     free(rhsr);
//     free(rhsz);
//     free(Az);
//     free(Ar);
//     return(status);
//   }
// #endif
//   return(status);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// int solvez_spdomestic(solver *slv,complex<double> *RHS,int transp) {
//   
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   int status=0,bcl,n;
//   tripletz *Mat;
//   cs_ci *A;
//   cs_cis *S ;
//   cs_cin *N ;
//   complex<double> *x;
// 
//   Mat = slv->Mat;
//   A = Mat->A;
//   S = Mat->S;
//   N = Mat->N;
//   n = A->n ;
//   x = cs_ci_malloc (n, sizeof (complex<double>)) ;    /* get workspace */
// 
// 
//   /*printf("\n");*/
//   /*printf("....... SOLVE .............\n");*/
//   bcl = N->L->nz;
//   /*  for(bcl=0;bcl<A->n;bcl++) printf("RHS[%d]=%e\n",bcl,RHS[bcl]);*/
// 
//   cs_ci_ipvec (N->pinv, RHS, x, n) ;            /* x = b(p) */
//   if (transp) {
//     cs_ci_utsolve (N->U, x) ;                    /* x = U\x */
//     cs_ci_ltsolve (N->L, x) ;                    /* x = L\x */
//   }
//   else {
// //    printf("solvez_spdomestic : x = L\\x  \n");
//     cs_ci_lsolve (N->L, x) ;                    /* x = L\x */
// //    printf("solvez_spdomestic : x = U\\x  \n");
//     cs_ci_usolve (N->U, x) ;                    /* x = U\x */
//   }
// //  printf("solvez_spdomestic : b(q) = x  \n");
//   cs_ci_ipvec (S->q, x, RHS, n) ;                   /* b(q) = x */
// //  printf("solvez_spdomestic : free\n");
//   cs_ci_free (x) ;
// //  printf("solvez_spdomestic : done...  \n");   
//   
// 
//   return(status); 
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int solvez_hips(solver *slv, complex<double>  *RHS, int transp, int verbose, bool debug)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
// #ifdef HIPS
// #ifdef HAVE_MPI
//   zHIPS_STRUC_C*id;
//   char *solver_name;
//   int proc_id, myid=-1, ierr, status, n, nnz,bcl, nglob;
//   double cput1, cput2, cput3;
//   INTS id_hips;
//   INTS *unknownlist;
//   char filename[50];
//   FILE *out;
//   complex<double> *x,*x2;
//   
//   cput1 = MPI_Wtime();
//   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
// 
//   if (debug && proc_id == 0) {printf("cpu=%d, %s starts... \n",proc_id,__FUNCTION__);}
//   
//   id = slv->parameters;
//   id_hips =id->id_hips;
//   n=id->n;
// 
//   if (debug) printf("cpu=%d, %s resolution du system numero %d, n=%d \n",proc_id,__FUNCTION__,id_hips,n);
// 
//   nnz=id->nnz;
//   unknownlist = id->unknownlist;
//   
// //   sprintf(filename,"solver_solve_hips.1.%d.txt",proc_id);
// //   out=fopen(filename,"w");   
// //   for (bcl=0; bcl < n; bcl++) fprintf(out, "%d => %d \n",bcl,unknownlist[bcl]);
// //   fclose(out);
//   
//   /** TEst set RHS like x=1 */
// //   x  = (complex<double> *) malloc(sizeof(complex<double>)*nnz);
// //   x2 = (complex<double> *) malloc(sizeof(complex<double>)*nnz);
// //   for (bcl=0;bcl < nnz;bcl++) x[bcl]=1.0; */
// // 
// //   zHIPS_MatrixVectorProduct (id_hips, x, x2);
// //   for (bcl=0;bcl < n;bcl++) RHS[bcl]=x2[bcl];
// 
//   int mode=1;
//   switch(mode) {
//     case 0:
// /*------------------------------------------------------------------------------
//       version locale */
// //       ierr = zHIPS_SetRHS (id_hips, n, unknownlist, RHS,  zHIPS_ASSEMBLY_OVW, zHIPS_ASSEMBLY_OVW, zHIPS_ASSEMBLY_FOOL);
//       ierr = zHIPS_SetRHS (id_hips, n, unknownlist, RHS,  zHIPS_ASSEMBLY_OVW, zHIPS_ASSEMBLY_OVW, zHIPS_ASSEMBLY_RESPECT);
//       break;
//       
// /*------------------------------------------------------------------------------
//     version globale*/
//     case 1:
//       ierr = zHIPS_SetGlobalRHS(id_hips, RHS, 0,  zHIPS_ASSEMBLY_OVW);
//       break;
//     }
// 
//   if (debug) printf("cpu=%d, %s zHIPS_SetRHS OK.... \n",proc_id,__FUNCTION__);
// 
// //   sprintf(filename,"solver_solve_hips.RHS.2.%d.txt",proc_id);
// //   out=fopen(filename,"w");   
// //   for (bcl=0; bcl < n; bcl++) fprintf(out, "%d => %d, %lf\n",bcl,unknownlist[bcl]-1,RHS[bcl]);
// //   fclose(out);
// 
//   cput1 = MPI_Wtime();
//   
// /*------------------------------------------------------------------------------
//   version locale */
// //   ierr = zHIPS_GetSolution(id_hips,n,unknownlist, RHS, zHIPS_ASSEMBLY_FOOL);
// 
// /*------------------------------------------------------------------------------
//   version globale*/
//   ierr = zHIPS_GetGlobalSolution(id_hips,RHS,-1);
//   
// //   sprintf(filename,"solver_solve_hips.SOL.2.%d.txt",proc_id);
// //   out=fopen(filename,"w");   
// //   for (bcl=0; bcl < n; bcl++) fprintf(out, "%d => %d, %lf\n",bcl,unknownlist[bcl],RHS[bcl]);
// //   fclose(out);
// 
// //   zHIPS_SetOptionINT (id_hips,zHIPS_DISABLE_PRECOND, 1);
// 
//   cput2 = MPI_Wtime(); 
//   cput3 = cput2-cput1;
//   if (verbose == 1 && proc_id == 0) {
//     printf("%s:  HIPS Resolution time =%lf \n",__FUNCTION__,cput3);
//     }
//     
// //    ierr = MPI_Allreduce (  &n, &nglob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );   
// //    x  = (complex<double> *) malloc(sizeof(complex<double>)*nglob);
// //    ierr = zHIPS_GetGlobalSolution(id_hips,x,-1);
// //    sprintf(filename,"solver_solve_hips.GlobalSolution.%d.txt",proc_id);
// //    out=fopen(filename,"w");
// //    for (bcl=0; bcl < nglob; bcl++) fprintf(out,"%d => %lf\n",bcl,x[bcl]);
// //    fclose(out); 
// //    free(x);
// //    RHS = x;
// 
//   ierr=0;
//   return(ierr);
// #endif
// #endif
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int solvez_pastix_sequential(solver *slv,complex<double> *RHS,int transp, int verbose, bool debug) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   int status,bcl,bcl2;
// #ifdef PASTIX
//   PASTIX_STRUC_C *id;
//   char *solver_name;
//   complex<double> *Ax;
//   int *Ap, *Ai,n,nnz;
//   tripletz *Mat;
//   int INCXY;
//   char filerhs[50];
//   int myid,sz,i,ierr;
//   FILE *out;
//  
//   status = 0;
//   solver_name=slv->name;
//   if (strcmp(solver_name,"PASTIX") == 0) {
//     id = slv->parameters;
//     Mat=slv->Mat;
//     Ap=Mat->i;
//     Ai=Mat->j;
//     Ax=Mat->x;
//     n=Mat->nrow;
//     nnz=Mat->nnz;
// 
// //    sprintf(filerhs,"RHS.%d.%d.txt",sz,myid);
// //    out=fopen(filerhs,"w");
// //    for (i=0; i <n; i++)   fprintf(out,"      %d -> %e+i%e \n",i,RHS[i]);
// //    fclose(out);
// //    MPI_Barrier(MPI_COMM_WORLD);
//    printf("solve_pastix_parallel n=%d nnz=%d\n",n,nnz);
//    
//   id->iparm[IPARM_START_TASK] = API_TASK_SOLVE;
//   id->iparm[IPARM_END_TASK]   = API_TASK_REFINE;
// 
//   z_pastix(&(id->pastix_data), MPI_COMM_WORLD,
//          n, Ap, Ai, Ax,
//          id->perm, id->invp, RHS, 1, id->iparm, id->dparm);
//    
// //    sprintf(filerhs,"SOL.%d.%d.txt",sz,myid);
// //    out=fopen(filerhs,"w");
// //    for (i=0; i <n; i++)   fprintf(out,"      %d -> %e+i%e \n",i,RHS[i]);
// //    fclose(out);
// //    
// //    exit(-23);
//    
//     return(status);
//   }
// #endif
//   return(status);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int solvez_pastix_parallel(solver *slv,complex<double> *RHS,int transp, int verbose, bool debug) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   int status,bcl,bcl2;
//   int ierr=-1;
// #ifdef PASTIX
// #ifdef HAVE_MPI
//   PASTIX_STRUC_C *id;
//   char *solver_name;
//   complex<double> *Ax;
//   int *Ap, *Ai,n,nnz;
//   tripletz *Mat;
//   int INCXY;
//   char filerhs[50];
//   int myid,sz,i;
//   FILE *out;
// 
// #ifdef HAVE_MPI
//   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//   ierr = MPI_Comm_size(MPI_COMM_WORLD, &sz);
// #endif
//   
//   status = 0;
//   solver_name=slv->name;
// //   if (strcmp(solver_name,"PASTIX") == 0) {
//   id = slv->parameters;
//   Mat=slv->Mat;
//   Ap=Mat->i;
//   Ai=Mat->j;
//   Ax=Mat->x;
//   n=Mat->nrow;
//   nnz=Mat->nnz;
// 
// //  sprintf(filerhs,"RHSz.%d.%d.txt",sz,myid);
// //  out=fopen(filerhs,"w");
// //  for (i=0; i <n; i++)   fprintf(out,"    %d -> %e+i%e \n",i,RHS[i]);
// //  fclose(out);
//   
//   MPI_Barrier(MPI_COMM_WORLD);
//   if(verbose==1) printf("cpu=%d : %s n=%d nnz=%d\n",myid,__FUNCTION__, n,nnz);
//  
//   id->iparm[IPARM_START_TASK] = API_TASK_SOLVE;
//   id->iparm[IPARM_END_TASK]   = API_TASK_REFINE;
// 
//   z_pastix(&(id->pastix_data), MPI_COMM_WORLD,
//            n, Ap, Ai, Ax, id->perm, id->invp, RHS, 1, id->iparm, id->dparm);
// 
// //    sprintf(filerhs,"SOLz.%d.%d.txt",sz,myid);
// //    out=fopen(filerhs,"w");
// //    for (i=0; i <n; i++)   fprintf(out,"      %d -> %e+i%e \n",i,RHS[i]);
// //    fclose(out);
//   
//   ierr=0;
//   
// #endif 
// #endif 
//   return(ierr);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int solvez_mumps(solver *slv,complex<double> *RHS,int transp, int verbose, bool debug) {
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// #ifdef MUMPS
//   ZMUMPS_STRUC_C *id;
//   char *solver_name;
//   int ierr, rank;
//   double cput1, cput2, cput3;
//   int status;
//   
// #ifdef HAVE_MPI
//   cput1 = MPI_Wtime();
// #endif
//   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   solver_name=slv->name;
//   if (strcmp(solver_name,"MUMPS") == 0) {
//     id = slv->parameters;
//     if (transp == 1) {
//       (*id).ICNTL(9)=2;
//     }
//     if (rank == 0) {
//       /*printf("Solve MUMPS ...\n");*/
//       (*id).rhs = RHS;
//     }
//     (*id).job=3;
//     zmumps_c(id);
//   }
// #ifdef HAVE_MPI
//   cput2 = MPI_Wtime(); 
//   cput3 = cput2-cput1;
//   if (rank == 0) {
//     if(verbose==1) printf(" ZMUMPS resolution time =%lf\n",cput3);
//   }
// #endif
//   ierr=(*id).info[0];
//   if(ierr==-1) {
//     printf("cpu %d: ZMUMPS solver failed, status=%d due to cpu %d\n",rank, ierr, id->info[1]);
//     }
//   if(ierr==-22) {
//     printf("cpu %d: ZMUMPS solver failed, status=%d array fault %d\n",rank, ierr, id->info[1]);
//     }
//   return(ierr);
// #endif
// }
// 
// /*=======================================================================*/
// /*=======================================================================*/
// /*====================  FIN SOLVE =======================================*/
// /*=======================================================================*/
// 
// 
// 
// 
// 
//   
// 
// 
// 
// /*=======================================================================*/
// /*                 Solver Terminate                                      */
// /*************************************************************************/
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// int freez_solver(solver *slv)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   /*mumps_solver mumpsslv;*/
//   /*int len,ierr,myid;*/
//   int status=0;
// #ifdef MUMPS
//   ZMUMPS_STRUC_C *id_mumps;
// #endif
// #ifdef UMFPACK
//   UMFPACK_STRUC_C *id_umf;
// #endif
// #ifdef PASTIX
//   PASTIX_STRUC_C *id_pastix;
//   complex<double> *Ax=0, *rhs=0;
//   int *Ap=0, *Ai=0,n,nnz;
// //   tripletz *Mat;
// //   Mat=slv->Mat;
//   
// #endif
//   /*double *Control;*/
// #ifdef HIPS
//   zHIPS_STRUC_C*id;
//   INTS id_hips;
//   FILE *out;
//   int ierr;
// #endif
//   char *solver_name=0;
//   void *Symbolic, *Numeric ;
//   tripletz *Mat;
// 
//   if(slv!=0) {
//     solver_name=slv->name;
//     Mat = slv->Mat;
//     }
// 
//   if (solver_name != NULL) {
// //  status=initialisation_terminate();  
// 
// /* ---------------- FREE MUMPS SOLVER -------------------- */
// #ifdef MUMPS
//     if (strcmp(solver_name,"MUMPS") == 0) {
//       id_mumps = slv->parameters;
//       (*id_mumps).job=JOB_END;
//       zmumps_c(id_mumps);
// /* ---------------- FREE MUMPS SOLVER -------------------- */
//       }
// #endif
// 
// /* ---------------- FREE UMFPACK SOLVER -------------------- */
// #ifdef UMFPACK
//     if (strcmp(solver_name,"UMFPACK") == 0) {
//       id_umf = slv->parameters;
//       Symbolic=id_umf->Symbolic;
//       Numeric=id_umf->Numeric;
//       umfpack_zi_free_symbolic (&Symbolic) ;
//       umfpack_zi_free_numeric (&Numeric) ;
// /* ---------------- FREE UMFPACK SOLVER -------------------- */
//       }
// #endif
// 
// #ifdef LAPACK
//     if (strcmp(solver_name,"LAPACK") == 0) {
//       free(Mat->a);
//       free(Mat->pivot);
//       /* ---------------- FREE LAPACK SOLVER -------------------- */
//       }
// #endif
//     
// /* ---------------- FREE SpDOMESTIC SOLVER -------------------- */
//     if (strcmp(solver_name,"SpDOMESTIC") == 0) {
//       cs_ci_sfree (Mat->S) ;
//       cs_ci_nfree (Mat->N) ;
// /* ---------------- FREE SpDOMESTIC SOLVER -------------------- */
//       }
// 
// /* ---------------- FREE PASTIX SOLVER -------------------- */
// #ifdef PASTIX
//     if (strcmp(solver_name,"PASTIX") == 0) {
//       Ap=Mat->i;
//       Ai=Mat->j;
//       Ax=Mat->x;
//       id_pastix = slv->parameters;
//       id_pastix->iparm[IPARM_START_TASK] = API_TASK_CLEAN;
//       id_pastix->iparm[IPARM_END_TASK]   = API_TASK_CLEAN;
// 
//       z_pastix(&(id_pastix->pastix_data), MPI_COMM_WORLD,
//          n, Ap, Ai, Ax,
//          id_pastix->perm, id_pastix->invp, rhs, 1, 
//          id_pastix->iparm, id_pastix->dparm);
//       free(id_pastix->perm);
//       free(id_pastix->invp);
// /* ---------------- FREE PASTIX SOLVER -------------------- */
//       }
// #endif
// /* ---------------- FREE HIPS SOLVER -------------------- */
// #ifdef HIPS
//     if (strcmp(solver_name,"HIPS") == 0) {
//   /************************************************/
//   /* Free zHIPS internal structure for problem id  */
//   /************************************************/
//       id = slv->parameters;
//       id_hips =id->id_hips;
//       ierr = zHIPS_Clean(id_hips);
//       zHIPS_ExitOnError(ierr);
//       HIPSSOLVERNBZ--;
//       
//   /**********************************/
//   /* Free zHIPS internal structures  */
//   /**********************************/
//       if(HIPSSOLVERNBZ==0) {
//         ierr = zHIPS_Finalize();
//         zHIPS_ExitOnError(ierr);
//         }
//       
//       }
// #endif
// /* ---------------- FREE HIPS SOLVER -------------------- */
// 
// //    if(Mat->private_memory==1) {
// //       free(Ax);
// //       free(Ap);
// //       free(Ai);
// //       }
//     free(slv->parameters);
//     free(slv->name);
//     /*free(slv);*/
//     slv->parameters = NULL;
//     slv->name = NULL;
//     /*slv = NULL;*/
//     }
//   else {
//     printf(" free already done with this solver \n");
//     }
//     
// //   slv=0;
//   
//   return(status);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int solver_terminate(solver *solveur)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status=0, ierr;
// 
//   if(solveur==0) return(0);
// 
// //  switch (gLinearSolverID) { /// HERE!!!
//   switch (solveur->id) {
//     case SOLVER_ID_DOMESTIC:     /* DOMESTIC */
//       break;
// 
//     case SOLVER_ID_LAPACK:       /* LAPACK */
//       break;
// 
//     case SOLVER_ID_SUNPERF:      /* SUNPERF */
//       break;
// 
//     case SOLVER_ID_UMFPACK:      /* UMFPACK */
//       status=free_solver(solveur);
//       break;
// 
//     case SOLVER_ID_GMRES:        /* ITERATIV SOLVER */
// //      printf("iterativ packed matrix solver\n");
//       status=free_solver(solveur);
//       break;
// 
//     case SOLVER_ID_MUMPS:        /* MUMPS */
//       status=free_solver(solveur);
//       break;
// 
//     case SOLVER_ID_MUMPS_SYM:    /* MUMPS-SYM */
//       status=free_solver(solveur);
//       break;
// 
//     case SOLVER_ID_SpDOMESTIC:    /* SpDOMESTIC */
//       status=free_solver(solveur);
//       break;
// 
// //    case 8:                      /* HYPRE */
// //      break;
//     
//     case SOLVER_ID_HIPS:           /* HIPS */
//       status=free_solver(solveur);
//       break;  
// 
//     case SOLVER_ID_PASTIX:         /* PASTIX */
//       status=free_solver(solveur);
//       break;  
// 
//     default:
//       exit(-1);
//       break;
//     }
// 
//   if(status!=0) {
//     printf("solver =%d\n",solveur->id);
//     printf("solver_terminate status=%d ...\n",status);
//     }
// 
//   return (status);
// }
// 
// 
// /*************************************************************************/
// /*=======================================================================*/
// /*======================FIN TERMINATE ===================================*/
// /*=======================================================================*/
// 
// 
// 
// /*=======================================================================*/
// /*======================UTILITAIRES MODE PARALLEL =======================*/
// /*=======================================================================*/
// 
// #ifdef HAVE_MPI
// 
// int d_coo2glob(int ncol,int nnz,int *Ap, int *Ai,double *Ax,
//                         int *ncolglob, int *nnzglob, int **Apglobin,int **Aiglobin,double **Axglobin) {  
//       // Ap, Ai, Ax in COO FORMAT    
//       int mpid,sz;
//       int nc,nn;
//       int *Aptmpg;
//       int *Aitmpg;
//       double *Axtmpg;
//       int bcl;
//       int *debtab;
//       int *debtab2; 
//       char filemat[50];
//       FILE *f_out;
//       int *Apglob;
//       int *Aiglob;
//       double *Axglob;
//       int debug=0;
//         
//       MPI_Comm_rank(MPI_COMM_WORLD, &mpid);
//       MPI_Comm_size(MPI_COMM_WORLD, &sz);
//       if(debug==1) printf("%d::: Entree  cooredispatch2csc\n",mpid );
//       if(debug==1) printf("%d::: Entree  cooredi... ncol=%d nnz=%d \n",mpid,ncol,nnz );
//       
// //      sprintf(filemat,"Matrixlocal.%d.txt",mpid);
// //      f_out = fopen(filemat,"w");
// //      for (bcl=0; bcl < nnz; bcl++) {
// //         fprintf(f_out,"%d %d %lf\n",Ap[bcl],Ai[bcl],Ax[bcl]);
// //      }
// //      fclose(f_out); 
// 
//      // Global Values
//       nc=ncol;
//       nn=nnz;
//       if(debug==1) printf("%d::: Entree  cooredi2... ncol=%d nnz=%d \n",mpid,nc,nn );
//       MPI_Allreduce(&nc, ncolglob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
//       MPI_Allreduce(&nn, nnzglob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
//       MPI_Barrier(MPI_COMM_WORLD);
//       if(debug==1) printf("%d::: ncolglob=%d nnzglob=%d \n",mpid,*ncolglob,*nnzglob );
//       nn = *nnzglob;
//       Aptmpg=(int *)malloc(sizeof(int)*nn);
//       Apglob=(int *)malloc(sizeof(int)*nn);
// //    MPI_Barrier(MPI_COMM_WORLD);
// //       printf("%d::: Alloc1 passed \n",mpid );
//       Aitmpg=(int *)malloc(sizeof(int)*nn);
//       Aiglob=(int *)malloc(sizeof(int)*nn);
// //    MPI_Barrier(MPI_COMM_WORLD);
// //       printf("%d::: Alloc2 passed \n",mpid );
//       Axtmpg=(double *)malloc(sizeof(double)*nn);
//       Axglob=(double *)malloc(sizeof(double)*nn);
// //     MPI_Barrier(MPI_COMM_WORLD);
// //       printf("%d::: Alloc2 passed \n",mpid );
//       memset(Aptmpg, 0, nn*sizeof(int));
//       memset(Aitmpg, 0, nn*sizeof(int));
//       memset(Axtmpg, 0.0, nn*sizeof(double));
// //    MPI_Barrier(MPI_COMM_WORLD);
// //       printf("%d::: memset passed \n",mpid );
//       
//       // Each proc fill in his part of globla Matrix
//       debtab =(int *)malloc(sizeof(int)*(sz+1));
//       memset(debtab, 0, (sz+1)*sizeof(int));
//       debtab2 =(int *)malloc(sizeof(int)*(sz+1));
//       memset(debtab2, 0, (sz+1)*sizeof(int));
//       for (bcl=mpid+1;bcl<(sz+1);bcl++)  {
//          debtab2[bcl]=nnz;
//       }
//       MPI_Allreduce(debtab2, debtab, sz+1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
//       free(debtab2);
// //      for (bcl=0; bcl<(sz+1);bcl++) printf("   %d::: debtab[%d]=%d \n",mpid,bcl,debtab[bcl]);
// //       printf("%d::: debtab=%d \n",mpid,debtab[mpid] );
// //       printf("%d::: fintab=%d nnzglob=%d\n",mpid,debtab[mpid]+nnz-1,*nnzglob );
// //       
// //       printf("%d::: Ap[1]=%d .........\n",mpid,Ap[1]);
// //       printf("%d::: Ap[%d]=%d .........\n",mpid,nnz-1,Ap[nnz-1]);
// //       printf("%d::: Ai[1]=%d .........\n",mpid,Ai[1]);
// //       printf("%d::: Ai[%d]=%d .........\n",mpid,nnz-1,Ai[nnz-1]);
//      
//      for (bcl=0; bcl < nnz; bcl++) {
//           Aptmpg[debtab[mpid]+bcl] = Ap[bcl];
//            Aitmpg[debtab[mpid]+bcl] = Ai[bcl];
//             Axtmpg[debtab[mpid]+bcl] = Ax[bcl];
//      }
//      MPI_Allreduce(Aptmpg, Apglob, nn, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
//      MPI_Allreduce(Aitmpg, Aiglob, nn, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
//      MPI_Allreduce(Axtmpg, Axglob, nn, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
//      free(Aptmpg);
//      free(Aitmpg);
//      free(Axtmpg);
//      
//      *Apglobin = Apglob;
//      *Aiglobin = Aiglob;
//      *Axglobin = Axglob;
//     
// //      for (bcl=0; bcl < 10; bcl++) {
// //       printf("%d::: %d => Apglob[bcl]=%d .........\n",mpid,bcl,Apglob[bcl]);
// //      }
// //      for (bcl=nn-10; bcl < nn; bcl++) {
// //       printf("%d::: %d => Apglob[bcl]=%d .........\n",mpid,bcl,Apglob[bcl]);
// //      }
//      
// //      sprintf(filemat,"MatrixGlobal.%d.txt",mpid);
// //      f_out = fopen(filemat,"w");
// //      for (bcl=0; bcl < nn; bcl++) {
// //         fprintf(f_out,"%d %d %lf\n",Apglob[bcl],Aiglob[bcl],Axglob[bcl]);
// //      }
// //      fclose(f_out); 
//      
//      return EXIT_SUCCESS;
// }
// 
// int d_globcoo2csc(int ncolglob, int nnzglob, int *Apglob,int *Aiglob,double *Axglob,
//       int **Ap,int **Ai,double **Ax) {
//       /* Thanks to Yousef Saad SPARSKIT */
//       int mpid,sz;
//       int *Aptmpg;
//       int *Aitmpg;
//       double *Axtmpg;
//       int ncol, nnz, k, i, j, k0, iad, bcl, bcl2;
//       double x;
//       char filemat[50];
//       FILE *f_out;
//       int debug=0;
//      
//       MPI_Comm_rank(MPI_COMM_WORLD, &mpid);
//       MPI_Comm_size(MPI_COMM_WORLD, &sz);
//       ncol=ncolglob;
//       nnz=nnzglob;
//       
//       Aptmpg=(int *)malloc(sizeof(int   )*(ncol+1));
//       Aitmpg=(int *)malloc(sizeof(int   )*nnz);
//       Axtmpg=(double *)malloc(sizeof(double)*nnz);
//     
// //       for(bcl=0; bcl<ncol+1;bcl++) Aptmpg[bcl]=0;
//       memset(Aptmpg, 0, (ncol+1)*sizeof(int));
//    MPI_Barrier(MPI_COMM_WORLD);
//       if(debug==1) printf("%d::: d_globcoo2csc malloc passed \n",mpid );
// // !
// // !  Determine the row lengths.
// // !
//       for(k = 0; k < nnz; k++) {
//           //printf("%d : k=%d, nnz=%d ncol=%d Apglob[k]=%d \n",mpid,k,nnz,ncol,Apglob[k]);
//           Aptmpg[Apglob[k]] ++;
//       }
//    MPI_Barrier(MPI_COMM_WORLD);
//        if(debug==1) printf("%d::: d_globcoo2csc row length passed \n",mpid );
// // !
// // !  The starting position of each row.
// // !
//   k = 0;
//   for (j=0; j<ncol+1; j++) {
//     k0 = Aptmpg[j];
//     Aptmpg[j] = k;
//     k = k+k0;
//   }
//    MPI_Barrier(MPI_COMM_WORLD);
//        if(debug==1) printf("%d::: d_globcoo2cs  starting position passed \n",mpid );
// // !
// // !  Go through the structure once more.  Fill in output matrix.
// // !
//   for (k=0;k<nnz;k++) {
//    i =  Apglob[k];
//    j =  Aiglob[k]; // Fortran Numbering
//    x =  Axglob[k];
//    iad = Aptmpg[i];
//    Axtmpg[iad] = x;
//    Aitmpg[iad] = j; // Fortran Numbering
//    Aptmpg[i] = iad+1;
//   }
//    MPI_Barrier(MPI_COMM_WORLD);
//        if(debug==1) printf("%d::: d_globcoo2csc fil matrix passed \n",mpid );
// 
// // !
// // !  Shift back IAO.
// // !
// //   for (j=ncol; j > 0; j--) {
// //     Aptmpg[j] = Aptmpg[j-1];
// //   }
//    for (j=0; j<ncol+1;j++) Aptmpg[j]++;
//        
//   *Ap=Aptmpg;
//   *Ai=Aitmpg;
//   *Ax=Axtmpg;
// 
//   return EXIT_SUCCESS;
// }
// 
// int z_coo2glob(int ncol,int nnz,int *Ap, int *Ai,complex<double> *Ax,
//                         int *ncolglob, int *nnzglob, int **Apglobin,int **Aiglobin,complex<double> **Axglobin) {  
//       // Ap, Ai, Ax in COO FORMAT    
//   int mpid,sz;
//   int nc,nn;
//   int *Aptmpg;
//   int *Aitmpg;
//   complex<double> *Axtmpg;
//   int bcl;
//   int *debtab;
//   int *debtab2; 
//   char filemat[50];
//   FILE *f_out;
//   int *Apglob;
//   int *Aiglob;
//   complex<double> *Axglob;
//   int debug=0;
//     
//   MPI_Comm_rank(MPI_COMM_WORLD, &mpid);
//   MPI_Comm_size(MPI_COMM_WORLD, &sz);
//   if(debug==1) printf("%d::: Entree  Zcooredispatch2csc\n",mpid );
//   if(debug==1) printf("%d::: Entree  Zcooredi... ncol=%d nnz=%d \n",mpid,ncol,nnz );
//   
// /*------------------------------------------------------------------------------
//   evaluate global values for number of columns and non-zero values */
//   nc=ncol;
//   nn=nnz;
//   if(debug==1) printf("%d::: Entree  Zcooredi2... ncol=%d nnz=%d \n",mpid,nc,nn );
//   
//   MPI_Allreduce(&nc, ncolglob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
//   MPI_Allreduce(&nn, nnzglob,  1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
//   MPI_Barrier(MPI_COMM_WORLD);
//   
//   if(debug==1) printf("%d::: Zncolglob=%d Znnzglob=%d \n",mpid,*ncolglob,*nnzglob );
//   
// /*------------------------------------------------------------------------------
//   allocate global arrays */
//   nn = *nnzglob;
//   Aptmpg=(int *)malloc(sizeof(int)*nn);
//   Apglob=(int *)malloc(sizeof(int)*nn);
//   MPI_Barrier(MPI_COMM_WORLD);
//   
//   if(debug==1) printf("%d::: Alloc1 passed \n",mpid );
//   Aitmpg=(int *)malloc(sizeof(int)*nn);
//   Aiglob=(int *)malloc(sizeof(int)*nn);
//   MPI_Barrier(MPI_COMM_WORLD);
//   
//   if(debug==1) printf("%d::: Alloc2 passed \n",mpid );
//   Axtmpg=(complex<double>*)malloc(sizeof(complex<double>)*nn);
//   Axglob=(complex<double>*)malloc(sizeof(complex<double>)*nn);
//   MPI_Barrier(MPI_COMM_WORLD);
//   
//   if(debug==1) printf("%d::: Alloc2 passed \n",mpid );
//   
// /*------------------------------------------------------------------------------
//   initialize global arrays */
//   memset(Aptmpg, 0, nn*sizeof(int));
//   memset(Aitmpg, 0, nn*sizeof(int));
//   memset(Axtmpg, (0.0, 0.0), nn*sizeof(complex<double>));
//   MPI_Barrier(MPI_COMM_WORLD);
//   
//   if(debug==1) printf("%d::: memset passed \n",mpid );
//   
// /*------------------------------------------------------------------------------
//   each proc fills in his part of global matrices */
//   debtab =(int *)malloc(sizeof(int)*(sz+1));
//   memset(debtab, 0, (sz+1)*sizeof(int));
//   debtab2 =(int *)malloc(sizeof(int)*(sz+1));
//   memset(debtab2, 0, (sz+1)*sizeof(int));
//   for (bcl=mpid+1;bcl<(sz+1);bcl++)  {
//     debtab2[bcl]=nnz;
//     }
//   MPI_Allreduce(debtab2, debtab, sz+1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD); // HERE !!!
//   
//   free(debtab2);
//   
// //  for (bcl=0; bcl<(sz+1);bcl++) printf("   %d::: debtab[%d]=%d \n",mpid,bcl,debtab[bcl]);
// //   printf("%d::: debtab=%d \n",mpid,debtab[mpid] );
// //   printf("%d::: fintab=%d nnzglob=%d\n",mpid,debtab[mpid]+nnz-1,*nnzglob );
// //   
// //   printf("%d::: Ap[1]=%d .........\n",mpid,Ap[1]);
// //   printf("%d::: Ap[%d]=%d .........\n",mpid,nnz-1,Ap[nnz-1]);
// //   printf("%d::: Ai[1]=%d .........\n",mpid,Ai[1]);
// //   printf("%d::: Ai[%d]=%d .........\n",mpid,nnz-1,Ai[nnz-1]);
//      
//   for (bcl=0; bcl < nnz; bcl++) {
//     Aptmpg[debtab[mpid]+bcl] = Ap[bcl];
//     Aitmpg[debtab[mpid]+bcl] = Ai[bcl];
//     Axtmpg[debtab[mpid]+bcl] = Ax[bcl];
//     }
//   MPI_Allreduce(Aptmpg, Apglob, nn, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
//   MPI_Allreduce(Aitmpg, Aiglob, nn, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
//   MPI_Barrier(MPI_COMM_WORLD);
//    
//   if(debug==1) printf("%d::: Reduce Integer passed \n",mpid );
//    
//   MPI_Allreduce(Axtmpg, Axglob, nn, MPI_DOUBLE_COMPLEX , MPI_SUM, MPI_COMM_WORLD);
//   MPI_Barrier(MPI_COMM_WORLD);
//   
//   if(debug==1) printf("%d::: Reduce Copmplex passed \n",mpid );
//    
//   free(Aptmpg);
//   free(Aitmpg);
//   free(Axtmpg);
//   
//   *Apglobin = Apglob;
//   *Aiglobin = Aiglob;
//   *Axglobin = Axglob;
//     
// //   for (bcl=0; bcl < 10; bcl++) {
// //    printf("%d::: %d => Apglob[bcl]=%d .........\n",mpid,bcl,Apglob[bcl]);
// //   }
// //   for (bcl=nn-10; bcl < nn; bcl++) {
// //    printf("%d::: %d => Apglob[bcl]=%d .........\n",mpid,bcl,Apglob[bcl]);
// //   }
//   
// //   nbsyslinearz2++;
// //   sprintf(filemat,"Z_MatrixGlobal.%d.%d.txt",nbsyslinearz2,mpid);
// //   f_out = fopen(filemat,"w");
// //   for (bcl=0; bcl < nn; bcl++) {
// //         fprintf(f_out,"%d %d  %g %g\n",Apglob[bcl],Aiglob[bcl],Axglob[bcl]);
// //   }
// //   fclose(f_out); 
//   
//   return EXIT_SUCCESS;
// }
// 
// int z_globcoo2csc(int ncolglob, int nnzglob, int *Apglob,int *Aiglob,complex<double> *Axglob,
//   int **Ap,int **Ai,complex<double> **Ax) {
//   /* Thanks to Yousef Saad SPARSKIT */
//   int mpid,sz;
//   int *Aptmpg;
//   int *Aitmpg;
//   complex<double> *Axtmpg;
//   int ncol, nnz, k, i, j, k0, iad, bcl, bcl2;
//   complex<double> x;
//   char filemat[50];
//   FILE *f_out;
//   int debug=0;
//      
//   MPI_Comm_rank(MPI_COMM_WORLD, &mpid);
//   MPI_Comm_size(MPI_COMM_WORLD, &sz);
//   ncol=ncolglob;
//   nnz=nnzglob;
//   
//   Aptmpg=(int *)malloc(sizeof(int   )*(ncol+1));
//   Aitmpg=(int *)malloc(sizeof(int   )*nnz);
//   Axtmpg=(complex<double> *)malloc(sizeof(complex<double>)*nnz);
//     
// //   for(bcl=0; bcl<ncol+1;bcl++) Aptmpg[bcl]=0;
//   memset(Aptmpg, 0, (ncol+1)*sizeof(int));
//   MPI_Barrier(MPI_COMM_WORLD);
//       
//   if(debug==1) printf("%d::: Z_d_globcoo2csc malloc passed \n",mpid );
// 
// /*------------------------------------------------------------------------------
//    Determine the row lengths, use Aptmpg as cardinal counter */
//   for(k = 0; k < nnz; k++) {
// /*------------------------------------------------------------------------------
//     Apglob[k] is column index */
//     Aptmpg[Apglob[k]] ++;
//     }
//   MPI_Barrier(MPI_COMM_WORLD);
//        
//   if(debug==1) printf("%d::: Z_d_globcoo2csc row length passed \n",mpid );
// 
// /*------------------------------------------------------------------------------
//   The starting position of each row, set Aptmpg as pointer */
//   k = 0;
//   for (j=0; j<ncol+1; j++) {
//     k0 = Aptmpg[j];
//     Aptmpg[j] = k;
//     k = k+k0;
//     }
//    
//   MPI_Barrier(MPI_COMM_WORLD);
//        
//   if(debug==1) printf("%d::: Z_d_globcoo2cs  starting position passed \n",mpid );
// 
// /*------------------------------------------------------------------------------
//   Go through the structure once more.  Fill in output matrix */
//   for (k=0;k<nnz;k++) {
//     i =  Apglob[k];
//     j =  Aiglob[k]; // Fortran Numbering
//     x =  Axglob[k];
//     iad = Aptmpg[i];
//     Axtmpg[iad] = x;
//     Aitmpg[iad] = j; // Fortran Numbering
//     Aptmpg[i] = iad+1;
//     }
//    
//   MPI_Barrier(MPI_COMM_WORLD);
//   
//   if(debug==1) printf("%d::: Z_d_globcoo2csc fil matrix passed \n",mpid );
// 
// // !
// // !  Shift back IAO.
// // !
// //   for (j=ncol; j > 0; j--) {
// //     Aptmpg[j] = Aptmpg[j-1];
// //   }
//    for (j=0; j<ncol+1;j++) Aptmpg[j]++;
//        
//   *Ap=Aptmpg;
//   *Ai=Aitmpg;
//   *Ax=Axtmpg;
// 
//   return EXIT_SUCCESS;
// }
// 
// #endif
