
/*******************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2017

  Unstructured Ocean Grid initiative

Developers:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
Contributors:
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France
  Frédéric Dupont    Canada environement, Canada

*******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex>
#include <cmath>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace std;
  
static size_t call_id=0;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mpi_synchronicity(MPI_Comm communicator, int CPUid, int nCPUs, const char *caller, int line, bool disciplined)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  verify that CPUs are all synchronized at collectve tasks
  
  TODO : communicator dependance to be implemented

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status=0;
//   bool spokesman=(gCPU_ID==gCPU_MASTER);
  
//   printf("%s : caller=%s line=%d cpu=%d \n", __func__, caller, line, CPUid);
  
  call_id++;

#ifdef HAVE_MPI
  int check[nCPUs];
  bool first=true, ok=true;
  
  if(disciplined) {
/*------------------------------------------------------------------------------
    force output to be ordered */
    status = MPI_Barrier(communicator);
    fflush(stdout);
    status = MPI_Barrier(communicator);
    for(int p=0;p<nCPUs;p++) {
      if(CPUid==p) {
        printf("%s id=%d : caller=%s line=%d cpu=%3d\n", __func__, call_id, caller, line, CPUid);
        }
      status = MPI_Barrier(communicator);
      fflush(stdout);
      status = MPI_Barrier(communicator);
      }
    }
  else {
    printf("%s id=%d : caller=%s line=%d cpu=%3d\n", __func__, call_id, caller, line, CPUid);
    }
  
  status = MPI_Barrier(communicator);
  fflush(stdout);
  status = MPI_Barrier(communicator);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  transmit line number
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status = MPI_Allgather(&line, 1, MPI_INT, check, 1, MPI_INT, communicator);
  first=true;
  for(int p=0;p<nCPUs;p++) {
    if(check[p]!=check[0]) {
      if(first) printf(" not ok \n");
      first=false; ok=false;
      printf("%s cpu=%3d : lines differ %d %d \n", __func__, p, check[p], check[0]);
      fflush(stdout);
      }
    }
  status = MPI_Barrier(communicator);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  transmit call id
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status = MPI_Allgather(&call_id, 1, MPI_INT, check, 1, MPI_INT, communicator);
  first=true;
  for(int p=0;p<nCPUs;p++) {
    if(check[p]!=check[0]) {
      if(first) printf(" not ok \n");
      first=false;
//       ok=false;
      printf("%s cpu=%3d : call ids differ %d %d \n", __func__, p, check[p], check[0]);
      fflush(stdout);
      }
    }
  
  status = MPI_Barrier(communicator);
  fflush(stdout);
  status = MPI_Barrier(communicator);
  if(ok and CPUid==0) printf("%s id=%d : caller=%s line=%d ok\n", __func__, call_id, caller, line);
  status = MPI_Barrier(communicator);
  fflush(stdout);
  status = MPI_Barrier(communicator);
#endif
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mpi_synchronicity(int CPUid, int nCPUs, const char *caller, int line, bool disciplined)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status=0;
  
#ifdef HAVE_MPI
  status=mpi_synchronicity(MPI_COMM_WORLD, CPUid, nCPUs, caller, line, disciplined);
#endif
  
  return(status);
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int P_MPI_Barrier(MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status=0;
//   
// #ifdef HAVE_MPI
//   status = MPI_Barrier(communicator);
// #endif
//   return(status);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int P_MPI_Barrier(MPI_Comm communicator, bool check)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status=0;
//   
// #ifdef HAVE_MPI
//   status = MPI_Barrier(communicator);
// #endif
//   return(status);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int P_MPI_BarrierFlush(MPI_Comm communicator, FILE *stream)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status=0;
//   
// #ifdef HAVE_MPI
//   status=MPI_Barrier(communicator);
//   status=fflush(stream);
//   status=MPI_Barrier(communicator);
// #else
//   status=fflush(stream);
// #endif
//   
//   return(status);
// }


