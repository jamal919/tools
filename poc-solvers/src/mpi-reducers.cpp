
/*******************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/

#include <mpi.h>

#include <complex>
using namespace std;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> T P_sum_template(T value, MPI_Datatype MPItype)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#if PARALLEL && HAVE_MPI
  T global;
  int status __attribute__((unused));
  global=0;
/**-----------------------------------------------------------------------------
  transmit local cardinal and sum for each domain */
  status = MPI_Allreduce ( &value, &global, 1, MPItype, MPI_SUM, MPI_COMM_WORLD );
  return(global);
#else
  return(value);
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> T P_average_template(MPI_Comm communicator, T value, MPI_Datatype MPItype)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#if PARALLEL && HAVE_MPI
  T global;
  int status __attribute__((unused));
  int size;
  status=MPI_Comm_size(communicator,&size);
  global=0;
/**-----------------------------------------------------------------------------
  transmit local cardinal and sum for each domain */
  status = MPI_Allreduce ( &value, &global, 1, MPItype, MPI_SUM, communicator );
  global/=(double) size;
  return(global);
#else
  return(value);
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> T P_norm_template(MPI_Comm communicator, T value, MPI_Datatype MPItype)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#if PARALLEL && HAVE_MPI
  T global;
  int status __attribute__((unused));
  global=0;
  T norm=value*value;
/**-----------------------------------------------------------------------------
  transmit local norm and sum for each domain, then take square root */
  status = MPI_Allreduce ( &norm, &global, 1, MPItype, MPI_SUM, communicator);
  global=sqrt(global);
  return(global);
#else
  return(value);
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> T P_max_template(MPI_Comm communicator, T value, MPI_Datatype MPItype)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#if PARALLEL && HAVE_MPI
  T global=-INFINITY;
  const bool debug=false;
  
  int status __attribute__((unused));
  
  if(debug) printf("%s cpu %d: value=%g\n" ,__func__, gCPU_ID,value );fflush(stdout);
  
  status=MPI_Allreduce (&value, &global, 1, MPItype,  MPI_MAX, communicator);
  
  if(debug) printf("%s cpu %d: global=%g\n",__func__, gCPU_ID,global);fflush(stdout);
  
  return(global);
#else
  return(value);
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> T P_min_template(T value, MPI_Datatype MPItype)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#if PARALLEL && HAVE_MPI
  T global;
  int status __attribute__((unused));
  status=MPI_Allreduce (&value, &global, 1, MPItype,  MPI_MIN, MPI_COMM_WORLD);
  return(global);
#else
  return(value);
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int P_sum(int value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int global;
  global=P_sum_template(value, MPI_INT);
  return(global);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  size_t P_sum(size_t value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  size_t global;
  global=P_sum_template(value, MPI_UNSIGNED_LONG);
  return(global);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double P_sum(double value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double global;
  global=P_sum_template(value, MPI_DOUBLE);
  return(global);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double P_average(double value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double global;
  global=P_average_template(MPI_COMM_WORLD, value, MPI_DOUBLE);
  return(global);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double P_norm(double value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double global;
  global=P_norm_template(MPI_COMM_WORLD, value, MPI_DOUBLE);
  return(global);
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   double P_norm(complex<double> value)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   double global;
//   global=P_norm_template(value, MPI_DOUBLE_COMPLEX);
//   return(global);
// }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> P_sum(complex<double> value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> global;
  global=P_sum_template(value, MPI_DOUBLE_COMPLEX);
  return(global);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double P_max(double value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double global;
  global=P_max_template(MPI_COMM_WORLD, value, MPI_DOUBLE);
  return(global);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double P_min(double value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double global;
  global=P_min_template(value, MPI_DOUBLE);
  return(global);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int P_maxloc(double & value, int & n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#if PARALLEL && HAVE_MPI
  int status __attribute__((unused));
  
  struct {
    double x;
    int n;
    } in, out;

  in.x=value;
  in.n=n;
  status=MPI_Allreduce (&in, &out, 1, MPI_DOUBLE_INT,  MPI_MAXLOC, MPI_COMM_WORLD );
  value=out.x;
  n=out.n;
  return(0);
#else
  return(0);
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int P_minloc(double & value, int & n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#if PARALLEL && HAVE_MPI
  int status __attribute__((unused));
  
  struct {
    double x;
    int n;
    } in, out;

  in.x=value;
  in.n=n;
  status=MPI_Allreduce (&in, &out, 1, MPI_DOUBLE_INT,  MPI_MINLOC, MPI_COMM_WORLD );
  value=out.x;
  n=out.n;
  return(0);
#else
  return(0);
#endif
}


