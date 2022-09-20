
/**************************************************************************

  POC-SOLVERS interface, 2006-2019

  Part of the SIROCCO national service (INSU, France)

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France

***************************************************************************/

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 MPI communication function and classes
 
 to be strictly synchronized with T-UGOm sources (native development sources)
 
 not used for the moment in poc-solvers

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <time.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "solvers-functions.h"
#include "linear-gmres.h"
#include "mpi-communications.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int exchange_atom_template(const exchange_t & exchange, int CPUid, int nCPUs, MPI_Datatype MPItype, T **to_send, T **to_receive, size_t blocksize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  synchronize (i.e. export/import) values of partition-wise local vectors with
  identical sizes between CPUs
  
  Use exchange_t class for communication handling
  
  MPI_Barrier use not wanted here 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k,l,m,n,p; 
  int status=0, error=0; 
  T **bufferR, **bufferS;
  int tagZcount, count;
  int nSend, nRecv, tagz;
  MPI_Comm communicator=exchange.communicator;

#ifdef HAVE_MPI
  MPI_Request Srequest[nCPUs],   Rrequest[nCPUs];
  MPI_Status  Smpistatus[nCPUs], Rmpistatus[nCPUs];
#endif
  
  bufferS=new T*[nCPUs];

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  send values from CPUid to the other CPUs
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nSend = 0;
  for (p=0;p<nCPUs;p++) {
    if(exchange.cpu_export.count[p]==0) continue;
    bufferS[p]=new T[exchange.cpu_export.count[p]*blocksize];
    T* pointer=bufferS[p];
    for(k=0;k<exchange.cpu_export.count[p];k++) {
      n=exchange.cpu_export.distributed[p][k];
      memmove(pointer, to_send[n], blocksize*sizeof(T));
      pointer+=blocksize;
      }
    tagz = CPUid*(nCPUs-1) + p;
    count=exchange.cpu_export.count[p];
#ifdef HAVE_MPI
    status=MPI_Issend(bufferS[p], blocksize*count, MPItype, p, tagz, communicator, &Srequest[nSend]);
#endif
    nSend ++;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  receive values from the other CPUs to CPUid
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nRecv = 0;
  bufferR=new T*[nCPUs];

  for (p=0;p<nCPUs;p++) {
    if(exchange.cpu_import.count[p] == 0) continue;
    bufferR[p] = new T [exchange.cpu_import.count[p]*blocksize];
    tagz = p*(nCPUs-1) + CPUid;
    count=exchange.cpu_import.count[p];
/**----------------------------------------------------------------------------
    receive list of nodes exported from the other CPUs to CPUid*/
#ifdef HAVE_MPI
    status = MPI_Irecv(bufferR[p], blocksize*count, MPItype, p, tagz, communicator, &Rrequest[nRecv]);
#endif
    nRecv ++;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  wait exchange completion
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#ifdef HAVE_MPI
  status = MPI_Waitall(nSend, Srequest, Smpistatus); 
  status = MPI_Waitall(nRecv, Rrequest, Rmpistatus); 
#endif

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  update buffer and release memory
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for (p=0; p<nCPUs; p++) {
    if(exchange.cpu_import.count[p] == 0) continue;
    T* pointer=bufferR[p];
    for(k=0;k<exchange.cpu_import.count[p];k++) {
      n=exchange.cpu_import.distributed[p][k];
//       memmove(&(to_receive[n]), &(bufferR[p]), blocksize*sizeof(T));
      memmove(to_receive[n], pointer, blocksize*sizeof(T));
      pointer+=blocksize;
      }
    delete[] bufferR[p];
    } 
  delete[] bufferR;
  
  for (p=0;p<nCPUs;p++) {
    if(exchange.cpu_export.count[p]==0) continue;
    delete[] bufferS[p];
    }
  delete[] bufferS;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int exchange_atom(const exchange_t & exchange, int CPUid, int nCPUs, char **to_send, char **to_receive, size_t blocksize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=exchange_atom_template(exchange, CPUid, nCPUs, MPI_CHAR, to_send, to_receive, blocksize);
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int exchange_atom(const exchange_t & exchange, int CPUid, int nCPUs, complex<double> **to_send, complex<double> **to_receive, size_t blocksize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=exchange_atom_template(exchange, CPUid, nCPUs, MPI_DOUBLE_COMPLEX, to_send, to_receive, blocksize);
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int exchange_atom_template(const exchange_t & exchange, int CPUid, int nCPUs, MPI_Datatype MPItype, T *to_send, T *to_receive, size_t blocksize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  synchronize (i.e. export/import) values of a partition-wise local vector between
  CPUs
  
  Use exchange_t class for communication handling
  
  MPI_Barrier use not wanted here 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k,l,m,n,p; 
  int status=0, error=0; 
  T **bufferR, **bufferS;
  int tagZcount, count;
  int nSend, nRecv, tagz;
  MPI_Comm communicator=exchange.communicator;

#ifdef HAVE_MPI
  MPI_Request Srequest[nCPUs], Rrequest[nCPUs];
  MPI_Status  Smpistatus[nCPUs], Rmpistatus[nCPUs];
  MPI_Request *requestZcount, *requestZnode;
#endif
  
  bufferS=new T*[nCPUs];

/*-----------------------------------------------------------------------------
  send list of nodes exported from CPUid to the other CPUs */
  nSend = 0;
  for (p=0;p<nCPUs;p++) {
    if(exchange.cpu_export.count[p]==0) continue;
    bufferS[p]=new T[exchange.cpu_export.count[p]*blocksize];
    T* target=bufferS[p];
    T* source;
    for(k=0;k<exchange.cpu_export.count[p];k++) {
      n=exchange.cpu_export.distributed[p][k];
      source=to_send+n*blocksize;
      memmove(target, source, blocksize*sizeof(T));
      target+=blocksize;
      }
/**----------------------------------------------------------------------------
    What is tagz for ? */
    tagz = CPUid*(nCPUs-1) + p;
    count=exchange.cpu_export.count[p];
//     printf("%s cpu=%3d : send to p=%d count=%d\n",__func__,CPUid,p,blocksize*count);
#ifdef HAVE_MPI
    status=MPI_Issend(bufferS[p], blocksize*count, MPItype, p, tagz, communicator, &Srequest[nSend]);
#endif
    nSend ++;
    }

  nRecv = 0;
  bufferR=new T*[nCPUs];

  for (p=0;p<nCPUs;p++) {
    if(exchange.cpu_import.count[p] == 0) continue;
    bufferR[p] = new T [exchange.cpu_import.count[p]*blocksize];
    tagz = p*(nCPUs-1) + CPUid;
    count=exchange.cpu_import.count[p];
/**----------------------------------------------------------------------------
    receive list of nodes exported from the other CPUs to CPUid*/
//     printf("%s cpu=%3d : receive from p=%d count=%d\n",__func__,CPUid,p,blocksize*count);
#ifdef HAVE_MPI
    status = MPI_Irecv(bufferR[p], blocksize*count, MPItype, p, tagz, communicator, &Rrequest[nRecv]);
#endif
    nRecv ++;
    }

/**----------------------------------------------------------------------------
  wait exchange completion*/
#ifdef HAVE_MPI
  status = MPI_Waitall(nSend, Srequest, Smpistatus); 
  status = MPI_Waitall(nRecv, Rrequest, Rmpistatus); 
#endif

  for (p=0; p<nCPUs; p++) {
    if(exchange.cpu_import.count[p] == 0) continue;
    T* target;
    T* source=bufferR[p];
    for(k=0;k<exchange.cpu_import.count[p];k++) {
      n=exchange.cpu_import.distributed[p][k];
      target=to_receive+n*blocksize;
      memmove(target, source, blocksize*sizeof(T));
      source+=blocksize;
      }
    delete[] bufferR[p];
    } 
  delete[] bufferR;
  
  for (p=0;p<nCPUs;p++) {
    if(exchange.cpu_export.count[p]==0) continue;
    delete[] bufferS[p];
    }
  delete[] bufferS;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int exchange_atom(const exchange_t & exchange, int CPUid, int nCPUs, char *to_send, char *to_receive, size_t blocksize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=exchange_atom_template(exchange, CPUid, nCPUs, MPI_CHAR, to_send, to_receive, blocksize);
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int exchange_atom_template(const exchange_t & exchange, int CPUid, int nCPUs, MPI_Datatype MPItype, T *to_send, T *to_receive, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  synchronize (i.e. export/import) values of a partition-wise local vector between
  CPUs
  
  Use exchange_t class for communication handling
  
  based on distributed index (in array addressing)
  
  MPI_Barrier use not wanted here 

-------------------------------------------------------------------------------*/
{
  int k,l,m,n,p; 
  int status=0, error=0;
  T **bufferR, **bufferS;
  int tagZcount, count;
  int nSend, nRecv, tagz;
  MPI_Comm communicator=exchange.communicator;

#ifdef HAVE_MPI
  MPI_Request Srequest[nCPUs], Rrequest[nCPUs];
  MPI_Status  Smpistatus[nCPUs], Rmpistatus[nCPUs]; 
  MPI_Request *requestZcount, *requestZnode;
#endif
  
/*------------------------------------------------------------------------------
  null partition issue */
//   status=synchronicity(CPUid, gnCPUs, __FUNCTION__, __LINE__,false);
  
  status=exchange.check();
// /**----------------------------------------------------------------------------
//   sequential mode */
//   if(status==-1) return(0);

  bufferS=new T*[nCPUs];
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  send list of nodes exported from CPUid to the other CPUs 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nSend = 0;
  if(exchange.cpu_export.count!=0)
  for (p=0;p<nCPUs;p++) {
    if(exchange.cpu_export.count[p]==0) continue;
    bufferS[p]=new T[exchange.cpu_export.count[p]];
    for(k=0;k<exchange.cpu_export.count[p];k++) {
      n=exchange.cpu_export.distributed[p][k];
      bufferS[p][k]=to_send[n];
      if(isnan(to_send[n])) {
        printf("%s cpu=%3d : sent to p=%d : n=%d k=%d to_send[n]=%lf\n",__func__, CPUid,p,n, k, to_send[n]);
        error=-1;
        }
      }
/**----------------------------------------------------------------------------
    What is tagz for ? */
    tagz = CPUid*(nCPUs-1) + p;
    count=exchange.cpu_export.count[p];
#ifdef HAVE_MPI
    status=MPI_Issend(bufferS[p], count, MPItype, p, tagz, communicator, &Srequest[nSend]);
#endif
    nSend ++;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  receive value imported from the other CPUs 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nRecv = 0;
  bufferR=new T*[nCPUs];

  if(exchange.cpu_import.count!=0)
  for (p=0;p<nCPUs;p++) {
    if(exchange.cpu_import.count[p] == 0) continue;
    bufferR[p] = new T [exchange.cpu_import.count[p]];
    tagz = p*(nCPUs-1) + CPUid;
    count=exchange.cpu_import.count[p];
/**----------------------------------------------------------------------------
    receive list of nodes exported from the other CPUs to CPUid*/
#ifdef HAVE_MPI
    status = MPI_Irecv(bufferR[p], count, MPItype, p, tagz, communicator, &Rrequest[nRecv]);
#endif
    nRecv ++;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  wait exchange completion 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#ifdef HAVE_MPI
  bool check=false;
  if(check) printf("%s cpu=%d : waiting send completion (%d)\n", __func__, CPUid,nSend);
  if(nSend!=0) status = MPI_Waitall(nSend, Srequest, Smpistatus); 
  if(check) printf("%s cpu=%d : waiting recv completion (%d)\n", __func__, CPUid,nRecv);
  if(nRecv!=0) status = MPI_Waitall(nRecv, Rrequest, Rmpistatus); 
  if(check) printf("%s cpu=%d : transaction done\n", __func__, CPUid);
#endif

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  finalize receive transaction 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(exchange.cpu_import.count!=0)
  for (p=0; p<nCPUs; p++) {
    if(exchange.cpu_import.count[p] == 0) continue;
    if(to_receive==0) TRAP_ERR_EXIT(-1,"%s cpu=%3d : receiving buffer is null\n",__func__,CPUid);
    for(k=0;k<exchange.cpu_import.count[p];k++) {
      n=exchange.cpu_import.distributed[p][k];
      m=exchange.cpu_import.centralized[p][k];
//       if(debug) {
      if(debug) {
        T error=(to_receive[n]-bufferR[p][k])/(to_receive[n]+bufferR[p][k]);
        if(abs(error)>1.e-7) {
          printf("%s cpu=%3d : received from p=%d : n=%5d (centralized=%5d) k=%4d original=%lf received=%lf\n",__func__, CPUid, p, n, m, k, to_receive[n], bufferR[p][k]);
          }
        }
      to_receive[n]=bufferR[p][k];
      if(isnan(to_receive[n])) {
        printf("%s cpu=%3d : received from p=%d : n=%d k=%d to_receive[n]=%lf\n",__func__,CPUid,p,n, k, to_receive[n]);
        error=-1;
        }
      }
    delete[] bufferR[p];
    } 
  deletep(&bufferR);
  
  if(exchange.cpu_export.count!=0)
  for (p=0;p<nCPUs;p++) {
    if(exchange.cpu_export.count[p]==0) continue;
    delete[] bufferS[p];
    }
  deletep(&bufferS);
  
  status=error;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int exchange_atom(const exchange_t & exchange, int CPUid, int nCPUs, char *to_send, char *to_receive, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=exchange_atom_template(exchange, CPUid, nCPUs, MPI_CHAR, to_send, to_receive, debug);
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int exchange_atom(const exchange_t & exchange, int CPUid, int nCPUs, int *to_send, int *to_receive, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=exchange_atom_template(exchange, CPUid, nCPUs, MPI_INT, to_send, to_receive, debug);
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int exchange_atom(const exchange_t & exchange, int CPUid, int nCPUs, float *to_send, float *to_receive, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=exchange_atom_template(exchange, CPUid, nCPUs, MPI_FLOAT, to_send, to_receive, debug);
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int exchange_atom(const exchange_t & exchange, int CPUid, int nCPUs, double *to_send, double *to_receive, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=exchange_atom_template(exchange, CPUid, nCPUs, MPI_DOUBLE, to_send, to_receive, debug);
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int exchange_atom(const exchange_t & exchange, int CPUid, int nCPUs, complex<double> *to_send, complex<double> *to_receive, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=exchange_atom_template(exchange, CPUid, nCPUs, MPI_DOUBLE_COMPLEX, to_send, to_receive, debug);
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *cpu_table(distribution_t & distributor, const int CPUid, MPI_Comm communicator, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  returns a centralized vector (of global dimension ngdof) containing managing 
  processor id.
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int n, status;
  
  if(distributor.ngdof==-1) TRAP_ERR_EXIT(-1,"%s cpu=%3d : error ngdof=%d\n", __func__, CPUid, distributor.ngdof);
  if(distributor.ngdof==0)  TRAP_ERR_EXIT(-1,"%s cpu=%3d : error ngdof=%d\n", __func__, CPUid, distributor.ngdof);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  set current cpu contribution
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  int *buffer    =new int[distributor.ngdof];
  int *processors=new int[distributor.ngdof];

  for(n=0;n<distributor.ngdof;n++) buffer[n]=0;
  
  for(n=0;n<distributor.nsolved;n++) {
    int centralized=distributor.gsolved[n];
    buffer[centralized]=CPUid;
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  gather cpu community contribution
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#ifdef HAVE_MPI
  status = MPI_Allreduce (buffer, processors, distributor.ngdof, MPI_INT, MPI_SUM, communicator); 
#endif 
  
  delete[] buffer;
  
  return(processors);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int communication_set_distributed(distribution_t & distributor, communication_t & target, const int CPUid, const int nCPUs, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,p,n;
  
  int *centralized=new int[distributor.ndof];
  for(n=0;n<distributor.ndof;n++) centralized[n]=distributor.gindex[n];

  size_t size=occurence(-1, centralized, distributor.ndof);
  
  if(size!=0) TRAP_ERR_EXIT(-1,"%s cpu=%3d : globalized index not properly set (nnodes=%d)\n", __func__, CPUid, distributor.ndof);
    
  for (p=0;p<nCPUs;p++) {
    if(target.count[p]==0) continue;
    target.distributed[p]=new int[target.count[p]];
    for(k=0;k<target.count[p];k++) {
      n=target.centralized[p][k];
      int pos=vpos(n,centralized,distributor.ndof);
      if(pos==-1) TRAP_ERR_EXIT(-1,"%s cpu=%3d : error k=%d n=%d (nnodes=%d)\n", __func__, CPUid, k,n,distributor.ndof);
      target.distributed[p][k]=pos;
      }
    }
    
  delete[] centralized;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int communication_dualize(communication_t & source, communication_t & target, const int CPUid, const int nCPUs, MPI_Comm communicator, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  exchange information : send exports and receive imports  or reverse
    
  comments made for first case
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k, p, n, status;
  int *buffer, *bufferS[nCPUs], **bufferR,count;
  int tagZcount, tagZnode;
  int tnbSend[nCPUs][nCPUs],tnbRecv[nCPUs*nCPUs];
  int nSend, nRecv, tagz;
  
#ifdef HAVE_MPI
  MPI_Request Srequest[nCPUs],   Rrequest[nCPUs];
  MPI_Status  Smpistatus[nCPUs], Rmpistatus[nCPUs];
#endif

  status=P_MPI_Barrier(communicator);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  first send the number of exported nodes per CPU 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for (p=0;p<nCPUs;p++) {
    tnbSend[CPUid][p] = source.count[p];
    }

#ifdef HAVE_MPI
  status = MPI_Allgather(tnbSend[CPUid], nCPUs, MPI_INTEGER,  /* send buf,count,type */
                         tnbRecv,        nCPUs, MPI_INTEGER,  /* recv buf,count,type */
                         communicator);                       /* comm,flag           */
#endif
/**-----------------------------------------------------------------------------
  count export by p to q is count import from p to q */
  for (p=0;p<nCPUs;p++) {
    target.count[p]=tnbRecv[p*nCPUs+CPUid];
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  then send list of nodes (centralized index) buffers exported from CPUid to 
  the other CPUs
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nSend = 0;
  for (p=0;p<nCPUs;p++) {
    if(source.count[p]==0) continue;
    if(source.count[p]<0) TRAP_ERR_EXIT(-1,"cpu =%d send count=%d\n",CPUid,source.count[p]);
    bufferS[p]=new int[source.count[p]];
    for(k=0;k<source.count[p];k++) {
      n=source.centralized[p][k];
      bufferS[p][k]=n;
      }
/*------------------------------------------------------------------------------
    send list of nodes exported to p */
    count=source.count[p];
/**----------------------------------------------------------------------------
    */
    tagz = CPUid*(nCPUs-1) + p;
#ifdef HAVE_MPI
    status=MPI_Issend(bufferS[p], source.count[p], MPI_INTEGER, p, tagz, communicator, &Srequest[nSend]);
#endif
    nSend ++;
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  receive buffer as import tables
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nRecv = 0;
  bufferR=new int* [nCPUs];
  for (p=0;p<nCPUs;p++) {
    if(target.count[p] == 0) continue;
    bufferR[p] = new int [target.count[p]];
    tagz = p*(nCPUs-1) + CPUid;
/**----------------------------------------------------------------------------
    receive list of nodes exported from the other CPUs to CPUid*/
#ifdef HAVE_MPI
    status = MPI_Irecv(bufferR[p], target.count[p], MPI_INTEGER, p, tagz, communicator, &Rrequest[nRecv]);
#endif
    nRecv ++;
    }  

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  wait exchange completion
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#ifdef HAVE_MPI
  status = MPI_Waitall(nSend, Srequest, Smpistatus); 
  status = MPI_Waitall(nRecv, Rrequest, Rmpistatus); 
#endif

  for (p=0; p<nCPUs; p++) {
    if(target.count[p] == 0) continue;
    target.centralized[p]=new int[target.count[p]];
    for(k=0;k<target.count[p];k++) {
      target.centralized[p][k]=bufferR[p][k];
      }
    delete[] bufferR[p];
    }
  delete[] bufferR;

  for (p=0;p<nCPUs;p++) {
    if(source.count[p]==0) continue;
    delete[] bufferS[p];
    }
    
  status=P_MPI_Barrier(communicator);

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int communication_SetImportTables(distribution_t & distributor, int CPUid, MPI_Comm communicator, int **LocalTable, int **UniversalTable, int *Count, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  centralizer based import table setup
  
  root processor (that hold centralizer) parses and scatter informations to other 
  processors
  
  needed for mpi exchange setup in general case

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int j,k,l,p,status;
  size_t i,n,ngdof=distributor.ngdof, ndof=distributor.ndof;
  int *gprocessor=0, *processor=0;
  vector<int> solving_processors;
  
  centralizer_t centralizer=distributor.centralizer;
  
  if(debug) {
    printf("\n%s: initialize export table betwen CPU %d and other CPUs \n", __func__, CPUid);
    status=P_MPI_Barrier(communicator);
    }
    
#if 0
// de-activated 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  root processor build global dof's processor table
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(CPUid==gCPU_MASTER) {
    gprocessor=new int[ngdof];
    for(n=0;n<ngdof;n++) gprocessor[n]=-1;
/*------------------------------------------------------------------------------
    setup dof's handling processor table */    
    for(p=0;p<centralizer.nCPUs;p++) {
      for(i=0;i<centralizer.SolvedCount[p];i++) {
        n=centralizer.UniversalSolvedTable[p][i];
        gprocessor[n]=p;
        }
      }
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  then scatter to each processor
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  processor=new int[ngdof];
    
  status=vector_RootScatter(distributor, processor, gprocessor, CPUid, gnCPUs, gCPU_MASTER, communicator);
    
  if(CPUid==gCPU_MASTER) delete[] gprocessor;

#else
  processor=cpu_table(distributor, CPUid, communicator, debug);
#endif
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  build import tables
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  get list of processors from which local cpu needs to import */    
  for(n=0;n<ndof;n++) {
    p=processor[n];
    if(p==CPUid) continue;
    int k=vpos(processor[n], solving_processors);
    if(k==-1) solving_processors.push_back(processor[n]);
    }

/*------------------------------------------------------------------------------
  should have been allocated prior to call */    
  for(i=0;i<centralizer.nCPUs;i++) {
    Count[i]=0;
    LocalTable[i]    =0;
    UniversalTable[i]=0;
    }
    
  for(i=0;i<solving_processors.size();i++) {
    p=solving_processors[i];
    Count[p]=occurence(p, processor, ndof);
    LocalTable[p]    =new int[Count[p]];
    UniversalTable[p]=new int[Count[p]];
    int m=0;
    for(n=0;n<ndof;n++) {
      if(processor[n]==p) {
        LocalTable[p][m]=n;
        UniversalTable[p][m]=distributor.gindex[n];
        m++;
        }
      }
    }
  
  if(debug) {
    for(i=0;i<solving_processors.size();i++) {
      p=solving_processors[i];
      printf("\n%s: CPU %d imports %d values from CPU %d (over %d %d) \n", __func__, CPUid, Count[p], p, solving_processors.size(), centralizer.nCPUs);
      }
    status=P_MPI_Barrier(communicator);
    }
      
  delete[] processor;
  
  status=P_MPI_Barrier(communicator);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int communications_init(distribution_t distributor, int CPUid, int nCPUs, MPI_Comm communicator, exchange_t & exchange, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  contsruct import/export tables for parallel computing
  
  >>>>>> distributor dependent

  A processor (A) needs to know:

    1- which processors (B,C,...) to address to get values
       * must provide index (in B) of node n (in A)
       * count

    2- which processor to address to send values

  Information should be symetrical, so only step 1 or 2 needs to be 
  actually computed, then MPI exhanges might be used to complete the tables

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  init local/universal import tables
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  exchange.init(CPUid, nCPUs, communicator);

  status=communication_SetImportTables(distributor, CPUid, exchange.communicator, exchange.cpu_import.distributed, exchange.cpu_import.centralized, exchange.cpu_import.count, debug);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  exchange information : send imports and receive exports
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=communication_dualize(exchange.cpu_import, exchange.cpu_export, CPUid, nCPUs, exchange.communicator, debug);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  complete export information (build local importation tables from universal tables)
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=communication_set_distributed(distributor, exchange.cpu_export, CPUid, nCPUs, debug);

  if(debug) status=communications_summary(exchange, debug);

  status=P_MPI_Barrier(communicator);

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int communications_summary(const exchange_t & exchange, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  exchange diagnostic (summary)
  
  MPI_Barrier not to be used in this routine
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k,l,m,n,p,q;
  int status;
  FILE *summary;
  char filename[64];
  int discretisation;
  
  size_t nimported=exchange.cpu_import.undisciplined_size();
  size_t nexported=exchange.cpu_export.undisciplined_size();

  if(nimported==0 and nexported==0) {
    if(debug) printf("%s, cpu=%3d : no import and no export \n", __func__, exchange.CPUid);
    return(0);
    }
  
/**-----------------------------------------------------------------------------
  MPI safe */
  discretisation=exchange.discretisation;
//   sprintf(filename,"summary-%s.%d.%d.txt",discretisation_name(discretisation),exchange.nCPUs,exchange.CPUid);
  sprintf(filename,"summary-%s.%d.%d.txt","generic",exchange.nCPUs,exchange.CPUid);

  summary=fopen(filename,"w");

  p=exchange.CPUid;

  fprintf(summary,SEPARATOR_1);
  for(q=0;q<exchange.nCPUs;q++) {
//     if(p==q) continue;
    fprintf(summary,SEPARATOR_1);
    fprintf(summary,"Export tables from CPU %d to CPU %d: count = %d \n",p,q,exchange.cpu_export.count[q]);
    if(debug) {
      fprintf(summary,"local    : ");
      for(k=0;k<exchange.cpu_export.count[q];k++) {
        fprintf(summary,"%3d ",exchange.cpu_export.distributed[q][k]);
        }
      fprintf(summary,"\n");
      fprintf(summary,"universal : ");
      for(k=0;k<exchange.cpu_export.count[q];k++) {
        fprintf(summary,"%3d ",exchange.cpu_export.centralized[q][k]);
        }
      fprintf(summary,"\n");
      }
    }

  fprintf(summary,SEPARATOR_1);
  for(q=0;q<exchange.nCPUs;q++) {
//     if(p==q) continue;
    fprintf(summary,SEPARATOR_1);
    fprintf(summary,"Import tables to CPU %d from CPU %d : count = %d \n",p,q,exchange.cpu_import.count[q]);
    if(debug) {
      fprintf(summary,"local    : ");
      for(k=0;k<exchange.cpu_import.count[q];k++) {
        fprintf(summary,"%3d ",exchange.cpu_import.distributed[q][k]);
        }
      fprintf(summary,"\n");
      fprintf(summary,"universal : ");
      for(k=0;k<exchange.cpu_import.count[q];k++) {
        fprintf(summary,"%3d ",exchange.cpu_import.centralized[q][k]);
        }
      fprintf(summary,"\n");
      }
    }

//   fprintf(summary,SEPARATOR_1);
//   fprintf(summary, "CPU %d: z-total =%d u-total =%d\n",p, exchange.Count[p], exchange.Count[p]);
//   fprintf(summary, "CPU %d: z-solved=%d u-solved=%d\n",p, exchange.SolvedCount[p], exchange.SolvedCount[p]);

  fclose(summary);

}




