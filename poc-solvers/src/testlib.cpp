
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "solvers-interface.h"
#include "poc-solvers.h"

#define onebased 1
#define zerobased 0

/* -------------------------------------------------------------------------- */
/* triplet form of the matrix.  The triplets can be in any order. */
/* -------------------------------------------------------------------------- */


int main(int argc, char ** argv)
{
  if( argc>1 and strcmp(argv[1],"--version")==0 ) {
    printf("%d.%d\n",POC_SOLVERS_VERSION_RELEASE,POC_SOLVERS_VERSION_MAJOR);
    return 0;
    }
  
  triplet *M1;
  double *RHS;
  int bcl,bcl2,ierr,status,nz1,myid,nbr;
  char *solv_name;
  solver_t *solveur1,*solveur2,*solveur3,*solveur4,*solveur5;
  FILE *f_out;
  char *fileMat, *fileRHS;
  int *Ap, *Ai, *Atmp;
  double *Ax;
  clock_t t1,t2,t3,tsolve1,tsolve2;
  double duree;
  static int    *Arow ;
  static int    *Acol ;
  static double *Aval ;
  char           *type        = NULL; /* type of the matrix                                        */
  char           *rhstype     = NULL; /* type of the right hand side                               */
  int ncol;
#ifdef PASTIX
  driver_type_t  driver_type;
  printf("Nombre de threads = %d \n",omp_get_max_threads());
#endif
  
#ifdef HAVE_MPI
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif
  if(argc > 1) {
    fileMat=argv[1];
    fileRHS=argv[2];
  }
  
  
  if (myid == 0) {
    /* Lecture de la matrice et RHS au format coo */
    M1 = (triplet *) read_triplet(fileMat,onebased);
    M1->comptype=COO;
    M1->base=onebased;
    printf("nnz = %d stype=%d nrow=%d \n",M1->nnz,M1->stype,M1->nrow);
    RHS = (double *)malloc(sizeof(double)*M1->nrow);
    printf("sizeof(RHS) = %d  ---- 1--- \n",sizeof(RHS));
    for (bcl=0  ; bcl < 5; bcl++) RHS[bcl]=bcl;
    read_vect(fileRHS,M1->nrow,RHS);
    printf("Lecture Ok \n");
    printf("M1 compression = %d\n",  M1->comptype);
    for (bcl=0  ; bcl < 5; bcl++) {
        printf("RHS[%d]=%e ----\n",bcl,RHS[bcl]);
        printf("%d, i,j,v=%d,%d,%e \n",bcl,M1->i[bcl],M1->j[bcl],M1->x[bcl]);
      }
  }
  
  printf("RESOLUTION \n\n");
  

    if ( ishere("PASTIX")) {
/*
    M1 = (triplet *) read_triplet(fileMat,onebased);
    M1->comptype=COO;
    M1->base=onebased;
*/
/*Avec Pastix Read car ci dessu problème de transpose ou format...*/  
    M1 = (triplet *)  malloc(sizeof(triplet));
    M1->comptype=CSC;
    M1->base=onebased;
#ifdef PASTIX
    driver_type = MM;

    d_read_matrix(fileMat, &ncol, &Ap, &Ai, &Ax,
              &RHS, &type, &rhstype, driver_type, MPI_COMM_WORLD);
#endif
    RHS = (double *)malloc(sizeof(double)*ncol);
    
    printf("*****************************************\n");
    /*read_vect(fileRHS,M1->nrow,RHS);*/
    for (bcl=0  ; bcl < ncol; bcl++) RHS[bcl] = 1.0;
    t1=clock();
    /*----------------------------------------------*/
    /* Resolution avec PASTIX */
    printf("******* Resolution avec PASTIX  ***********************\n");
    solveur5=(solver_t *) init_solver_obsolete("PASTIX", M1->stype,0);
    printf("Ok 2\n");
    printf("main resolution avec %s \n",solveur5->name);
    ierr = 0;
    M1->i = Ap;
    M1->j = Ai;
    M1->x = Ax;
    M1->nrow = ncol;
    /*factorize*/
    for (bcl=0  ; bcl < 5; bcl++) {
        printf("RHS[%d]=%e ----\n",RHS[bcl]);
        printf("%d, i,j,v=%d,%d,%e \n",bcl,M1->i[bcl],M1->j[bcl],M1->x[bcl]);
      }
/*Bon pour format CSC ici coo
    bcl = M1->nrow+1;
    bcl2 = M1->i[bcl];
    printf("%d, i,j,v=%d,%d,%e \n",bcl2,bcl,M1->j[bcl2],M1->x[bcl2]);
    printf("Check ...\n");
    printf("ncol=%d\n", M1->nrow);
    for (bcl=0  ; bcl < M1->nrow; bcl++) {
        printf("bcl=%d, %d...........\n",bcl,M1->i[bcl]);
         for (bcl2=  M1->i[bcl]; bcl2 <  M1->i[bcl+1]; bcl2++) {
             printf("%d, %d, %e \n",bcl,M1->j[bcl2],M1->x[bcl2]);
         }
     }
   
*/
    printf("MAtrice dans testlib 1\n");
    for(bcl=0; bcl < 10; bcl++) {
      printf("Ap[%d]=%d Ai[...]= ",bcl,Ap[bcl]);
      for(bcl2=M1->i[bcl];bcl2<M1->i[bcl+1];bcl2++) {
          printf("%d ",M1->j[bcl2]);
      }
      printf("\n");
     }
     printf("\n");
     
    ierr = factorize(solveur5,M1,1); /*....................*/
    
    printf("MAtrice dans testlib 2\n");
    for(bcl=0; bcl < 10; bcl++) {
      printf("Ap[%d]=%d Ai[...]= ",bcl,Ap[bcl]);
      for(bcl2=M1->i[bcl];bcl2<M1->i[bcl+1];bcl2++) {
          printf("%d ",M1->j[bcl2]);
      }
      printf("\n");
     }
     printf("\n");
    
    /*solve*/
    ierr = solve(solveur5,RHS,0,0);
    printf("solve status=%d \n",ierr);
    /*print solution */
    f_out = fopen("solutionPASTIX.txt","w");
    for (bcl=0;bcl < M1->nrow;bcl++) {
      fprintf(f_out,"%e\n",RHS[bcl]);
    }
    fclose(f_out);
    t2=clock();
    duree = (t2-t1)  *1000. / CLOCKS_PER_SEC;
    printf("Duree resoltion avec PASTIX  : %e (ms)\n",duree);
    free_solver(solveur5);
    free(M1);
    printf("*****************************************\n");
  }
  else {
    printf("PASTIX SOLVER not present \n");
  }
 
  
 

#ifdef HAVE_MPI
ierr = MPI_Finalize();
#endif
}
