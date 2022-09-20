
// /*******************************************************************************
// 
//   T-UGOm hydrodynamic ocean model, 2006-2017
// 
//   Unstructured Ocean Grid initiative
// 
// *******************************************************************************/
// /** \file
// 
// \brief MPI reducer functions
// */
// /*----------------------------------------------------------------------------*/
// 
// #include <mpi.h>
// #include <complex>
// 
// #include "mpi-communications.h"
// 
// using namespace std;

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   template <typename T> T P_sum_template(T value, MPI_Datatype MPItype, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
// #if PARALLEL && HAVE_MPI
//   T global;
//   int status __attribute__((unused));
//   global=0;
// /**-----------------------------------------------------------------------------
//   transmit local cardinal and sum for each domain */
//   status = MPI_Allreduce ( &value, &global, 1, MPItype, MPI_SUM, communicator);
//   return(global);
// #else
//   return(value);
// #endif
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   template <typename T> T P_norm_template(T value, MPI_Datatype MPItype, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
// #if PARALLEL && HAVE_MPI
//   T global;
//   int status __attribute__((unused));
//   global=0;
//   T norm=value*value;
// /**-----------------------------------------------------------------------------
//   transmit local norm and sum for each domain, then take square root */
//   status = MPI_Allreduce ( &norm, &global, 1, MPItype, MPI_SUM, communicator);
//   global=sqrt(global);
//   return(global);
// #else
//   return(value);
// #endif
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   template <typename T> T P_max_template(T value, MPI_Datatype MPItype, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
// #if PARALLEL && HAVE_MPI
//   T global=-INFINITY;
//   const bool debug=false;
//   
//   int status __attribute__((unused));
//   
//   if(debug) printf("%s cpu %d: value=%g\n" ,__func__, gCPU_ID,value );fflush(stdout);
//   status=MPI_Allreduce (&value, &global, 1, MPItype,  MPI_MAX, communicator);
//   if(debug) printf("%s cpu %d: global=%g\n",__func__, gCPU_ID,global);fflush(stdout);
//   
//   return(global);
// #else
//   return(value);
// #endif
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// template <typename T> T P_min_template(T value, MPI_Datatype MPItype, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
// #if PARALLEL && HAVE_MPI
//   T global;
//   int status __attribute__((unused));
//   status=MPI_Allreduce (&value, &global, 1, MPItype,  MPI_MIN, communicator);
//   return(global);
// #else
//   return(value);
// #endif
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int P_sum(int value, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int global;
//   global=P_sum_template(value, MPI_INT, communicator);
//   return(global);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   size_t P_sum(size_t value, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   size_t global;
//   global=P_sum_template(value, MPI_UNSIGNED_LONG, communicator);
//   return(global);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   double P_norm(double value, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   double global;
//   global=P_norm_template(value, MPI_DOUBLE, communicator);
//   return(global);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   double P_sum(double value, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   double global;
//   global=P_sum_template(value, MPI_DOUBLE, communicator);
//   return(global);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   double P_max(double value, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   double global;
//   global=P_max_template(value, MPI_DOUBLE, communicator);
//   return(global);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   double P_min(double value, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   double global;
//   global=P_min_template(value, MPI_DOUBLE, communicator);
//   return(global);
// }
// 
// 
// // /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// // 
// //   double P_norm(complex<double> value)
// // 
// // /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// // {
// //   double global;
// //   global=P_norm_template(value, MPI_DOUBLE_COMPLEX);
// //   return(global);
// // }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   complex<double> P_sum(complex<double> value, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   complex<double> global;
//   global=P_sum_template(value, MPI_DOUBLE_COMPLEX, communicator);
//   return(global);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int P_maxloc(double & value, int & n, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
// #if PARALLEL && HAVE_MPI
//   int status __attribute__((unused));
//   
//   struct {
//     double x;
//     int n;
//     } in, out;
// 
//   in.x=value;
//   in.n=n;
//   status=MPI_Allreduce (&in, &out, 1, MPI_DOUBLE_INT,  MPI_MAXLOC, communicator);
//   value=out.x;
//   n=out.n;
//   return(0);
// #else
//   return(0);
// #endif
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int P_minloc(double & value, int & n, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
// #if PARALLEL && HAVE_MPI
//   int status __attribute__((unused));
//   
//   struct {
//     double x;
//     int n;
//     } in, out;
// 
//   in.x=value;
//   in.n=n;
//   status=MPI_Allreduce (&in, &out, 1, MPI_DOUBLE_INT,  MPI_MINLOC, communicator);
//   value=out.x;
//   n=out.n;
//   return(0);
// #else
//   return(0);
// #endif
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   template <typename T> int vector_Allreduce_template(const distribution_t & distributor,const T *v, T* & global, MPI_Datatype MPItype, int cpu, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// 
//   feed global vectors construction from partioned ones and returns global vector
//   
//   compliant wiht non-MPI and MPI context
//   
//   allows for the user to provide global buffer (thus avoiding countless 
//   allocation and de-allocation in linear solvers routine)
// 
// cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
// {
//   T *tmp;
//   int nloc,nglob,status,rank, i, m, n;
//   int root;
//   int allocated=0;
//   bool debug=false;
// 
// #ifdef HAVE_MPI
// 
//   nloc  = distributor.nsolved;
//   nglob = distributor.ngdof;
//   
//   if(nglob==-1) {
//     TRAP_ERR_EXIT(-1,"%s cpu=%d : distributor not initialized, local=%d, total=%d\n",__func__,cpu,nloc,nglob);
//     }
//     
//   if(global==0) {
//     global    = new T[nglob];
//     if(debug) printf("%s cpu=%d : memory allocation, size=%d\n",__func__,cpu,nglob);
//     if(global==0) {
//       TRAP_ERR_EXIT(-1,"%s cpu=%d : memory allocation failed, size=%d\n",__func__,cpu,nglob);
//       }
//     allocated=1;
//     }
//   else if(debug) printf("%s cpu=%d : memory already allocated, size=%d\n",__func__,cpu,nglob);
//   
//   tmp = new T[nglob];
//   if(tmp==0) {
//     TRAP_ERR_EXIT(-1,"%s cpu=%d : memory allocation failed, size=%d\n",__func__,cpu,nglob);
//     }
// 
//   for (n=0;n<nglob;n++) {
//     tmp[n]    = 0.0;
//     global[n] = 0.0;
//     }
// 
// /**----------------------------------------------------------------------------
//   construction of local values ... */
//   for (i=0; i <distributor.nsolved; i++) {
//     n=distributor.solved[i];
//     m = distributor.gsolved[i];
//     tmp[m] = v[n];
//     }
// /**----------------------------------------------------------------------------
//   global vector construction: ok if domain are disjoint (summation) */
//   status = MPI_Allreduce (tmp, global, nglob, MPItype, MPI_SUM, communicator ); 
//   delete[] tmp;
// #else
//   global=(T *) v;
// #endif
//   
//   return (allocated);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int vector_Allreduce(const distribution_t & distributor,const int *v, int* & global, int cpu, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int allocated=0;
//   
//   allocated=vector_Allreduce_template(distributor, v, global, MPI_INT, cpu, communicator); 
//   
//   return(allocated);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int vector_Allreduce(const distribution_t & distributor,const double *v, double* & global, int cpu, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int allocated=0;
//   
//   allocated=vector_Allreduce_template(distributor, v, global, MPI_DOUBLE, cpu, communicator); 
//   
//   return(allocated);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int vector_Allreduce(const distribution_t & distributor,const complex<double> *v, complex<double>* & global, int cpu, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int allocated=0;
//   
//   allocated=vector_Allreduce_template(distributor, v, global, MPI_DOUBLE_COMPLEX, cpu, communicator); 
//   
//   return(allocated);
// }

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   template <typename T> int vector_Allreduce_template(const int *loc2glob,const T *v, T* & global, MPI_Datatype MPItype, int cpu, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// 
//   feed global vectors construction from partioned ones and returns global vector
//   
//   compliant wiht non-MPI and MPI context
//   
//   allows for the user to provide global buffer (thus avoiding countless 
//   allocation and de-allocation in linear solvers routine)
// 
// cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
// {
//   T *tmp;
//   int nloc,nglob,status,rank, i, m, n;
//   int root;
//   int allocated=0;
//   bool debug=false;
// 
// #ifdef HAVE_MPI
// 
//   nloc  = distributor.nsolved;
//   nglob = distributor.ngdof;
//   
//   if(nglob==-1) {
//     TRAP_ERR_EXIT(-1,"%s cpu=%d : distributor not initialized, local=%d, total=%d\n",__func__,cpu,nloc,nglob);
//     }
//     
//   if(global==0) {
//     global    = new T[nglob];
//     if(debug) printf("%s cpu=%d : memory allocation, size=%d\n",__func__,cpu,nglob);
//     if(global==0) {
//       TRAP_ERR_EXIT(-1,"%s cpu=%d : memory allocation failed, size=%d\n",__func__,cpu,nglob);
//       }
//     allocated=1;
//     }
//   else if(debug) printf("%s cpu=%d : memory already allocated, size=%d\n",__func__,cpu,nglob);
//   
//   tmp = new T[nglob];
//   if(tmp==0) {
//     TRAP_ERR_EXIT(-1,"%s cpu=%d : memory allocation failed, size=%d\n",__func__,cpu,nglob);
//     }
// 
//   for (n=0;n<nglob;n++) {
//     tmp[n]    = 0.0;
//     global[n] = 0.0;
//     }
// 
// /**----------------------------------------------------------------------------
//   construction of local values ... */
//   for (i=0; i <distributor.nsolved; i++) {
//     n=distributor.solved[i];
//     m = distributor.gsolved[i];
//     tmp[m] = v[n];
//     }
// /**----------------------------------------------------------------------------
//   global vector construction: ok if domain are disjoint (summation) */
//   status = MPI_Allreduce (tmp, global, nglob, MPItype, MPI_SUM, communicator ); 
//   delete[] tmp;
// #else
//   global=(T *) v;
// #endif
//   
//   return (allocated);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int vector_Allreduce(const int *loc2glob,const int *v, int* & global, int cpu, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int allocated=0;
//   
//   allocated=vector_Allreduce_template(distributor, v, global, MPI_INT, cpu, communicator); 
//   
//   return(allocated);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int vector_Allreduce(const int *loc2glob,const double *v, double* & global, int cpu, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int allocated=0;
//   
//   allocated=vector_Allreduce_template(distributor, v, global, MPI_DOUBLE, cpu, communicator); 
//   
//   return(allocated);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int vector_Allreduce(const int *loc2glob,const complex<double> *v, complex<double>* & global, int cpu, MPI_Comm communicator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int allocated=0;
//   
//   allocated=vector_Allreduce_template(distributor, v, global, MPI_DOUBLE_COMPLEX, cpu, communicator); 
//   
//   return(allocated);
// }

