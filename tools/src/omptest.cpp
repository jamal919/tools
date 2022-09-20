#include <stdio.h>
#include <omp.h>

/*
/home/softs/poc-solvers/src/initialise.c (rev 42332e5716af) says:
  
  We have encountered a SMP issue with umfpack (requesting libopenblas_pthreads)
  libumfpack-5.7.6.so version, detected in August 2017
  
  issue : multi-threading disabled
*/

int main(void){
  printf("testing omp_get_num_procs()=%d , omp_get_max_threads()=%d\n",
    omp_get_num_procs(), omp_get_max_threads() );
  }
