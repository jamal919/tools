
#include <stdio.h>
#include <malloc.h>
#include <time.h>

#include "solvers-interface.h"
#include "poc-solvers.h"
// #include "revision.h"

#define onebased 1
#define zerobased 0

/* -------------------------------------------------------------------------- */
/* tripletz form of the matrix.  The tripletzs can be in any order. */
/* -------------------------------------------------------------------------- */


int main(int argc, char ** argv)
{
  tripletz *M1;
  complex<double> *RHS;
  int bcl,ierr,status,nz1,myid,nbr;
  char *solv_name;
  solver_t *solveur1,*solveur2,*solveur3,*solveur4;
  FILE *f_out;
  char *fileMat, *fileRHS;
  int *Ap, *Ai;
  complex<double> *Ax;
  clock_t t1,t2,t3,tsolve1,tsolve2;
  double duree;
  static int    *Arow ;
  static int    *Acol ;
  static complex<double> *Aval ;
  int chk_double=0;

//   printf("poc-solver testlibz program, " REVISION "\n");

  if(argc > 1) {
    fileMat=argv[1];
    fileRHS=argv[2];
    }

  if(chk_double==1) {
    /* Lecture de la matrice et RHS au format coo */
    M1 = (tripletz *) read_triplet(fileMat,onebased);
    M1->comptype=COO;
    M1->base=onebased;
    printf("nnz = %d stype=%d nrow=%d \n",M1->nnz,M1->stype,M1->nrow);
    RHS = (complex<double> *)malloc(sizeof(complex<double>)*M1->nrow); 
    printf("sizeof(RHS) = %d  ---- 1--- \n",sizeof(RHS));
    for (bcl=0  ; bcl < 5; bcl++) RHS[bcl]=bcl;
    read_vectz(fileRHS,M1->nrow,RHS);
    printf("Lecture Ok \n");
    printf("M1 compression = %d\n",  M1->comptype);
    }

  if(ishere("SpDOMESTIC")) {
    M1 = (tripletz *) read_tripletz(fileMat,onebased);
    M1->comptype=COO; 
    M1->base=onebased;
    printf("*****************************************\n");
    read_vectz(fileRHS,M1->nrow,RHS);
    t1=clock();
    /*----------------------------------------------*/
    /* Resolution avec SpDOMESTIC */
    solveur4=(solver_t *) init_solver_obsolete("SpDOMESTIC", M1->stype,0);
    /*----------------------------------------------*/
    if (solveur4 == NULL) {
      exit(-1);
      }
    printf("Ok 2\n");
    printf("main resolution avec %s \n",solveur4->name);
    ierr = 0;
    /*factorize*/
    printf("Compression =%d, base=%d \n",M1->comptype,M1->base);
    ierr = factorizez(solveur4,M1,0);
    printf("factorize retour status=%d\n",ierr);
    printf("Compression =%d, base=%d \n",M1->comptype,M1->base);
    /*solve*/
    printf("SpDOMESTIC solve...\n");
    ierr = solvez(solveur4,RHS,1,0); 
    printf("solve status=%d \n",ierr); 
     /*print solution */ 
    f_out = fopen("solutionSpDOMESTIC.txt","w"); 
    for (bcl=0;bcl < M1->nrow;bcl++) {
      fprintf(f_out,"%e\n",RHS[bcl]);
      }
    fclose(f_out); 
    t2=clock();
    duree = (t2-t1)  *1000. / CLOCKS_PER_SEC;
    printf("main resolution avec %s finie...\n",solveur4->name);
    printf("Duree resoltion avec SpDOMESTIC  : %e (ms)\n",duree);
    freez_solver(solveur4);
    }
  else {
    printf("SpDOMESTIC SOLVER not present \n");
  }

}
