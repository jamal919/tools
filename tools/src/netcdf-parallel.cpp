
/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2012

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France

***************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <cmath>
#include <netcdf.h>

#ifdef PARALLEL
#ifdef HAVE_MPI
#include "mpi.h"
#endif
#endif


#include "tools-structures.h"

#include "fe.def"

#include "poc-netcdf.def"
#include "netcdf-proto.h"
#include "map.h"
#include "fe.h"
#include "parallel.h"

// static int gCPU_ID=0,gCPU_MASTER=0;

static char *rootname=".",*output_path=".";
static int gArchiveExtent=SEQUENTIAL_OUTPUT;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_seqput_UG2D(const char *filename, const mesh_t & mesh, int id, const int *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n, N, recordsize;
  int status;
  size_t nitems;
  size_t start[2];
  size_t count[2],size[2];
  int ncid;
  cdfvar_t variable;
  
  status = nc_open(filename, NC_WRITE, &ncid);
  if(status != NC_NOERR) goto error;

  status = cdf_varinfo(ncid,id, &variable ,0);

  for(k=0;k<variable.ndim;k++) {
    size[k]=variable.dim[k].length;
    }

  start[0] = 0;
  count[0] = size[0];

  if(variable.type!=NC_INT)
    check_error(-1,"wrong buffer type", __LINE__, __FILE__, 1);

  status = nc_put_vara_int(ncid, id, start, count, z);
  if(status != NC_NOERR) goto error;

  status = nc_close(ncid);
  if(status != NC_NOERR) goto error;

  variable.destroy();
  return (0);

error:
  nc_check_error(status, __LINE__, __FILE__);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_seqput_UG2D(const char *filename, const mesh_t & mesh, int id, const float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n, N, recordsize;
  int status;
  size_t nitems;
  size_t start[2];
  size_t count[2],size[2];
  int ncid;
  cdfvar_t variable;
  
  status = nc_open(filename, NC_WRITE, &ncid);
  if(status != NC_NOERR) goto error;

//  status=init_ncvariable(&variable);
  status = cdf_varinfo(ncid,id, &variable ,0);

  if(variable.ndim!=1) {
    printf("poc_seqput_UG2D: inconsistent call for %s, number of dimensions=%d, should be 1\n", variable.name, variable.ndim);
    }
    
  for(k=0;k<variable.ndim;k++) {
    size[k]=variable.dim[k].length;
    }

  start[0] = 0;
  count[0] = size[0];

  if(variable.type!=NC_FLOAT)
    check_error(-1,"wrong buffer type", __LINE__, __FILE__, 1);

  status = nc_put_vara_float(ncid, id, start, count, z);
  if(status != NC_NOERR) goto error;

  status = nc_close(ncid);
  if(status != NC_NOERR) goto error;

  variable.destroy();
  return (0);

error:
  nc_check_error(status, __LINE__, __FILE__);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_seqput_UG2D(const char *filename, const mesh_t & mesh, int id, const double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n, N, recordsize;
  int status;
  size_t nitems;
  size_t start[5];
  size_t count[5],size[5];
  int ncid;
  cdfvar_t variable;
  
  status = nc_open(filename, NC_WRITE, &ncid);
  if(status != NC_NOERR) goto error;

  status = cdf_varinfo(ncid,id, &variable ,0);
  if(status != NC_NOERR) goto error;
  
  if(variable.ndim!=1) {
    printf("poc_seqput_UG2D: inconsistent call for %s, number of dimensions=%d, should be 1\n", variable.name, variable.ndim);
    }

  for(k=0;k<variable.ndim;k++) {
    size[k]=variable.dim[k].length;
    }

  start[0] = 0;
  count[0] = size[0];

  switch(variable.type) {
    case NC_DOUBLE:
      break;
    default:
      check_error(-1,"buffer type andvariable type do not match", __LINE__, __FILE__, 0);
      break;
    }
  // if(variable.type!=NC_DOUBLE)
  // check_error(-1,"wrong buffer type", __LINE__, __FILE__, 1);

  status = nc_put_vara_double(ncid, id, start, count, z);
//  status = nc_put_vara(ncid, id, start, count, z);
  if(status != NC_NOERR) goto error;

  status = nc_close(ncid);
  if(status != NC_NOERR) goto error;

  variable.destroy();
  return (0);

error:
  if(status!=NC_NOERR) {
    fprintf(stderr, "poc_seqput_UG2D error, file %s, id %d name=%s\n",filename,id, variable.name);
    nc_check_error(status, __LINE__, __FILE__);
    }
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int poc_put_UG2D(const char *filename, const mesh_t & mesh, int id, const int *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated, status;
  mesh_t global;
  int *gbuffer=0;
  int discretisation;
  extern int h_discretisation(const char *filename, int id);

/**----------------------------------------------------------------------------
  read discretisation so global output vector can be formed from local ones */
//  discretisation= h_discretisation (filename, id);
  
  switch(gArchiveExtent) {
    case SEQUENTIAL_OUTPUT:
      status=poc_seqput_UG2D(filename, mesh, id, buffer);
      break;
    case GLOBALIZED_OUTPUT:
      global=*(mesh.origin);
//      allocated=OutputVector(buffer, &gbuffer, mesh, discretisation, gCPU_ID, 0);
      if(gCPU_ID==gCPU_MASTER) status=poc_seqput_UG2D(filename, mesh, id, gbuffer);
      if(allocated) delete[] gbuffer;
      break;
    case PARTITIONED_OUTPUT:
      status=poc_seqput_UG2D(filename, mesh, id, buffer);
      break;
    default:
      check_error(-1, "illicit output mode", __LINE__, __FILE__, 1);
    }
  
#ifdef HAVE_MPI
  status = MPI_Barrier(MPI_COMM_WORLD);
#endif
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int poc_put_UG2D(const char *filename, const mesh_t & mesh, int id, const float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated, status;
  mesh_t global;
  float *gbuffer=0;
  int discretisation;
  extern int h_discretisation(const char *filename, int id);
  
//  discretisation= h_discretisation (filename, id);
 
#ifdef HAVE_MPI
  status = MPI_Barrier(MPI_COMM_WORLD);
#endif
  switch(gArchiveExtent) {
    case SEQUENTIAL_OUTPUT:
      status=poc_seqput_UG2D(filename, mesh, id, buffer);
      break;
    case GLOBALIZED_OUTPUT:
      global=*(mesh.origin);
//      allocated=OutputVector(buffer, &gbuffer, mesh, discretisation, gCPU_ID, 0);
#ifdef HAVE_MPI
  status = MPI_Barrier(MPI_COMM_WORLD);
#endif
      if(gCPU_ID==gCPU_MASTER) status=poc_seqput_UG2D(filename, mesh, id, gbuffer);
      if(allocated) delete[] gbuffer;
      break;
    case PARTITIONED_OUTPUT:
      status=poc_seqput_UG2D(filename, mesh, id, buffer);
      break;
    default:
      check_error(-1, "illicit output mode", __LINE__, __FILE__, 1);
    }
  
#ifdef HAVE_MPI
  status = MPI_Barrier(MPI_COMM_WORLD);
#endif
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int poc_put_UG2D(const char *filename, const mesh_t & mesh, int id, const double *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated, status;
  mesh_t global;
  double *gbuffer=0;
  int discretisation;
  extern int h_discretisation(const char *filename, int id);
  
  discretisation=  (filename, id);

  switch(gArchiveExtent) {
    case SEQUENTIAL_OUTPUT:
      status=poc_seqput_UG2D(filename, mesh, id, buffer);
      break;
    case GLOBALIZED_OUTPUT:
      global=*(mesh.origin);
//      allocated=OutputVector(buffer, &gbuffer, mesh, discretisation, gCPU_ID, 0);
      if(gCPU_ID==gCPU_MASTER) status=poc_seqput_UG2D(filename, mesh, id, gbuffer);
      else status=0;
      if(allocated) delete[] gbuffer;
      break;
    case PARTITIONED_OUTPUT:
      status=poc_seqput_UG2D(filename, mesh, id, buffer);
      break;
    default:
      check_error(-1, "illicit output mode", __LINE__, __FILE__, 1);
    }
  
  return(status);
  }









/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_seqput_UG3D(const char *filename, const mesh_t & mesh, int frame, int id, const int *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n, N, recordsize;
  int status;
  size_t nitems;
  size_t start[2];
  size_t count[2],size[2];
  int ncid;
  cdfvar_t variable;
  
  status = nc_open(filename, NC_WRITE, &ncid);
  if(status != NC_NOERR) goto error;

  status = cdf_varinfo(ncid,id, &variable ,0);

  for(k=0;k<variable.ndim;k++) {
    size[k]=variable.dim[k].length;
    }

  start[0] = frame;
  start[1] = 0;

  count[0] = 1;
  count[1] = size[1];

// #ifdef DEBUG_PARALLEL
//   printf("poc_seqput_UG3D, frame=%d, id=%d, address=%h\n", frame, id, z);
// #endif
  
  if(variable.type!=NC_INT)
    check_error(-1,"wrong buffer type", __LINE__, __FILE__, 1);

  status = nc_put_vara_int(ncid, id, start, count, z);
  if(status != NC_NOERR) goto error;

  status = nc_close(ncid);
  if(status != NC_NOERR) goto error;

  variable.destroy();
  return (0);

error:
  nc_check_error(status, __LINE__, __FILE__);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_seqput_UG3D(const char *filename, const mesh_t & mesh, int frame, int id, const float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n, N, recordsize;
  int status;
  size_t nitems;
  size_t start[2];
  size_t count[2],size[2];
  int ncid;
  cdfvar_t variable;
  
  status = nc_open(filename, NC_WRITE, &ncid);
  if(status != NC_NOERR) goto error;

  status = cdf_varinfo(ncid,id, &variable ,0);

  for(k=0;k<variable.ndim;k++) {
    size[k]=variable.dim[k].length;
    }

  start[0] = frame;
  start[1] = 0;

  count[0] = 1;
  count[1] = size[1];

  if(variable.type!=NC_FLOAT)
    check_error(-1,"wrong buffer type", __LINE__, __FILE__, 1);

  status = nc_put_vara_float(ncid, id, start, count, z);
  if(status != NC_NOERR) goto error;

  status = nc_close(ncid);
  if(status != NC_NOERR) goto error;

  variable.destroy();
  return (0);

error:
  nc_check_error(status, __LINE__, __FILE__);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_seqput_UG3D(const char *filename, const mesh_t & mesh, int frame, int id, const double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n, N, recordsize;
  int status;
  size_t nitems;
  size_t start[2];
  size_t count[2],size[2];
  int ncid;
  cdfvar_t variable;
  
// #ifdef DEBUG_PARALLEL
//   printf("poc_seqput_UG3D (double), frame=%d, id=%d, name=%s, address=%d\n", frame, id, z);
// #endif

  status = nc_open(filename, NC_WRITE, &ncid);
  if(status != NC_NOERR) goto error;

  status = cdf_varinfo(ncid,id, &variable ,0);
  if(status != NC_NOERR) goto error;
  
#ifdef DEBUG_PARALLEL
  printf("poc_seqput_UG3D (double), frame=%d, id=%d, name=%s, address=%d\n", frame, id, variable.name, z);
#endif

  for(k=0;k<variable.ndim;k++) {
    size[k]=variable.dim[k].length;
    }

  start[0] = frame;
  start[1] = 0;

  count[0] = 1;
  count[1] = size[1];

  if(variable.type!=NC_DOUBLE)
    check_error(-1,"wrong buffer type", __LINE__, __FILE__, 1);

  status = nc_put_vara_double(ncid, id, start, count, z);
  if(status != NC_NOERR) goto error;

  status = nc_close(ncid);
  if(status != NC_NOERR) goto error;

  variable.destroy();
  return (0);

error:
  if(status!=NC_NOERR) {
    fprintf(stderr, "poc_seqput_UG3D error, file %s, id %d name=%s\n",filename,id, variable.name);
    nc_check_error(status, __LINE__, __FILE__);
    }
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_seqput_UG3D(const char *filename, mesh_t & mesh, int id, double **z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n, N, recordsize;
  int status;
  size_t nitems;
  size_t start[2];
  size_t count[2],size[2];
  int ncid;
  cdfvar_t variable;

  status = nc_open(filename, NC_WRITE, &ncid);
  if(status != NC_NOERR) goto error;

  status = cdf_varinfo(ncid,id, &variable ,0);
  if(status != NC_NOERR) goto error;
  
  for(k=0;k<variable.ndim;k++) {
    size[k]=variable.dim[k].length;
    }

  start[0] = 0;
  start[1] = 0;

  count[0] = 1;
  count[1] = size[1];

  if(variable.type!=NC_DOUBLE)
    check_error(-1,"wrong buffer type", __LINE__, __FILE__, 1);

  for(k=0;k<size[0];k++) {
    start[0] = k;
    status = nc_put_vara_double(ncid, id, start, count, z[k]);
    if(status != NC_NOERR) goto error;
    }

  status = nc_close(ncid);
  if(status != NC_NOERR) goto error;

  variable.destroy();
  return (0);

error:
  if(status!=NC_NOERR) {
    fprintf(stderr, "poc_seqput_UG3D error, file %s, id %d name=%s\n",filename,id, variable.name);
    nc_check_error(status, __LINE__, __FILE__);
    }
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int poc_put_UG3D(const char *filename, const mesh_t & mesh, int frame, int id, const int *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated, status;
  mesh_t global;
  int *gbuffer=0;
  int discretisation;
  extern int h_discretisation(const char *filename, int id);

/**----------------------------------------------------------------------------
  read discretisation so global output vector can be formed from local ones */
  discretisation= h_discretisation (filename, id);
  
#ifdef PARALLEL
//  printf("poc_put_UG3D (int), frame=%d, id=%d, address=%h\n", frame, id, buffer);
#endif

  switch(gArchiveExtent) {
    case SEQUENTIAL_OUTPUT:
      status=poc_seqput_UG3D(filename, mesh, frame, id, buffer);
      break;
    case GLOBALIZED_OUTPUT:
      global=*(mesh.origin);
//      allocated=OutputVector(buffer, &gbuffer, mesh, discretisation, gCPU_ID, 0);
      if(gCPU_ID==gCPU_MASTER) status=poc_seqput_UG3D(filename, mesh, frame, id, gbuffer);
      if(allocated) delete[] gbuffer;
      break;
    case PARTITIONED_OUTPUT:
      status=poc_seqput_UG3D(filename, mesh, frame, id, buffer);
      break;
    default:
      check_error(-1, "illicit output mode", __LINE__, __FILE__, 1);
    }
  
#ifdef HAVE_MPI
  status = MPI_Barrier(MPI_COMM_WORLD);
#endif
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int poc_put_UG3D(const char *filename, const mesh_t & mesh, int frame, int id, const float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated, status;
  mesh_t global;
  float *gbuffer=0;
  int discretisation;
  extern int h_discretisation(const char *filename, int id);
  
  discretisation= h_discretisation (filename, id);

#ifdef PARALLEL
//  printf("poc_put_UG3D (float), frame=%d, id=%d, address=%d\n", frame, id, buffer);
#endif
 
#ifdef HAVE_MPI
  status = MPI_Barrier(MPI_COMM_WORLD);
#endif
  switch(gArchiveExtent) {
    case SEQUENTIAL_OUTPUT:
      status=poc_seqput_UG3D(filename, mesh, frame, id, buffer);
      break;
    case GLOBALIZED_OUTPUT:
      global=*(mesh.origin);
//      allocated=OutputVector(buffer, &gbuffer, mesh, discretisation, gCPU_ID, 0);
#ifdef HAVE_MPI
  status = MPI_Barrier(MPI_COMM_WORLD);
#endif
      if(gCPU_ID==gCPU_MASTER) status=poc_seqput_UG3D(filename, mesh, frame, id, gbuffer);
      if(allocated) delete[] gbuffer;
      break;
    case PARTITIONED_OUTPUT:
      status=poc_seqput_UG3D(filename, mesh, frame, id, buffer);
      break;
    default:
      check_error(-1, "illicit output mode", __LINE__, __FILE__, 1);
    }
  
#ifdef HAVE_MPI
  status = MPI_Barrier(MPI_COMM_WORLD);
#endif
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int poc_put_UG3D(const char *filename, const mesh_t & mesh, int frame, int id, const double *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated, status;
  mesh_t global;
  double *gbuffer=0;
  int discretisation;
  extern int h_discretisation(const char *filename, int id);
  
  discretisation= h_discretisation (filename, id);

// #ifdef HAVE_MPI
// #ifdef DEBUG_PARALLEL
//   printf("poc_put_UG3D (double), cpu=%d, frame=%d, id=%d, address=%d\n", gCPU_ID, frame, id, buffer);
// #endif
// #endif

  switch(gArchiveExtent) {
    case SEQUENTIAL_OUTPUT:
      status=poc_seqput_UG3D(filename, mesh, frame, id, buffer);
      break;
    case GLOBALIZED_OUTPUT:
      global=*(mesh.origin);
// #ifdef HAVE_MPI
// #ifdef DEBUG_PARALLEL
//      printf("poc_put_UG3D (double), cpu=%d, frame=%d, id=%d, discretisation=%d, address=%d\n", gCPU_ID, frame, id, discretisation, buffer);
//      status = MPI_Barrier(MPI_COMM_WORLD);
// #endif
// #endif
#ifdef HAVE_MPI
      status = MPI_Barrier(MPI_COMM_WORLD);
#endif
//      allocated=OutputVector(buffer, &gbuffer, mesh, discretisation, gCPU_ID, 0);
      
// #ifdef HAVE_MPI
// #ifdef DEBUG_PARALLEL
//   printf("poc_put_UG3D (double), cpu=%d, frame=%d, id=%d, discretisation=%d, address=%d\n", gCPU_ID, frame, id, discretisation, buffer);
//   status = MPI_Barrier(MPI_COMM_WORLD);
// #endif
// #endif
#ifdef HAVE_MPI
      status = MPI_Barrier(MPI_COMM_WORLD);
#endif
      if(gCPU_ID==gCPU_MASTER) status=poc_seqput_UG3D(filename, mesh, frame, id, gbuffer);
      else status=0;
#ifdef HAVE_MPI
      status = MPI_Barrier(MPI_COMM_WORLD);
#endif
      if(allocated) delete[] gbuffer;
      break;
    case PARTITIONED_OUTPUT:
      status=poc_seqput_UG3D(filename, mesh, frame, id, buffer);
      break;
    default:
      check_error(-1, "illicit output mode", __LINE__, __FILE__, 1);
    }
  
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int poc_put_UG3D(const char *filename, mesh_t & mesh, int id, double **buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated, status;
  mesh_t global;
  double *gbuffer=0;
  int discretisation;
  extern int h_discretisation(const char *filename, int id);
  
  discretisation= h_discretisation (filename, id);

  switch(gArchiveExtent) {
    case SEQUENTIAL_OUTPUT:
      status=poc_seqput_UG3D(filename, mesh, id, buffer);
      break;
    case GLOBALIZED_OUTPUT:
//       global=*(mesh.origin);
// #ifdef HAVE_MPI
//       status = MPI_Barrier(MPI_COMM_WORLD);
// #endif
//       allocated=OutputVector(buffer, &gbuffer, mesh, discretisation, gCPU_ID, 0);
// #ifdef HAVE_MPI
//       status = MPI_Barrier(MPI_COMM_WORLD);
// #endif
//       if(gCPU_ID==gCPU_MASTER) status=poc_seqput_UG3D(filename, mesh, id, gbuffer);
//       else status=0;
// #ifdef HAVE_MPI
//       status = MPI_Barrier(MPI_COMM_WORLD);
// #endif
//       if(allocated) delete[] gbuffer;
      break;
    case PARTITIONED_OUTPUT:
      status=poc_seqput_UG3D(filename, mesh, id, buffer);
      break;
    default:
      check_error(-1, "illicit output mode", __LINE__, __FILE__, 1);
    }
  
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_seqput_UG4D(const char *filename, mesh_t & mesh, int frame, int id, double **z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n, N, recordsize;
  int status;
  size_t nitems;
  size_t start[3];
  size_t count[3],size[3];
  int ncid;
  cdfvar_t variable;
  
  status = nc_open(filename, NC_WRITE, &ncid);
  if(status != NC_NOERR) goto error;

  status = cdf_varinfo(ncid,id, &variable ,0);

  for(k=0;k<variable.ndim;k++) {
    size[k]=variable.dim[k].length;
    }

  start[0] = frame;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = 1;
  count[2] = size[2];

  if(variable.type!=NC_DOUBLE)
    check_error(-1,"wrong buffer type", __LINE__, __FILE__, 1);

  for(k=0;k<size[1];k++) {
    start[1] = k;
    status = nc_put_vara_double(ncid, id, start, count, z[k]);
    if(status != NC_NOERR) goto error;
    }

  status = nc_close(ncid);
  if(status != NC_NOERR) goto error;


  variable.destroy();
  return (0);

error:
  nc_check_error(status, __LINE__, __FILE__);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_seqput_UG4D(const char *filename, mesh_t & mesh, int frame, int level, int id, double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n, N, recordsize;
  int status;
  size_t nitems;
  size_t start[3];
  size_t count[3],size[3];
  int ncid;
  cdfvar_t variable;
  
  status = nc_open(filename, NC_WRITE, &ncid);
  if(status != NC_NOERR) goto error;

  status = cdf_varinfo(ncid,id, &variable ,0);

  for(k=0;k<variable.ndim;k++) {
    size[k]=variable.dim[k].length;
    }

  start[0] = frame;
  start[1] = level;
  start[2] = 0;

  count[0] = 1;
  count[1] = 1;
  count[2] = size[2];
  
  if(variable.type!=NC_DOUBLE)
    check_error(-1,"wrong buffer type", __LINE__, __FILE__, 1);

  status = nc_put_vara_double(ncid, id, start, count, z);
  if(status != NC_NOERR) goto error;

  status = nc_close(ncid);
  if(status != NC_NOERR) goto error;

  variable.destroy();
  return (0);

error:
  nc_check_error(status, __LINE__, __FILE__);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_seqput_UG4D(const char *filename, mesh_t & mesh, int frame, int level, int id, float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n, N, recordsize;
  int status;
  size_t nitems;
  size_t start[3];
  size_t count[3],size[3];
  int ncid;
  cdfvar_t variable;
  
  status = nc_open(filename, NC_WRITE, &ncid);
  if(status != NC_NOERR) goto error;

  status = cdf_varinfo(ncid,id, &variable ,0);

  for(k=0;k<variable.ndim;k++) {
    size[k]=variable.dim[k].length;
    }

  start[0] = frame;
  start[1] = level;
  start[2] = 0;

  count[0] = 1;
  count[1] = 1;
  count[2] = size[2];
  
  if(variable.type!=NC_DOUBLE)
    check_error(-1,"wrong buffer type", __LINE__, __FILE__, 1);

  status = nc_put_vara_float(ncid, id, start, count, z);
  if(status != NC_NOERR) goto error;

  status = nc_close(ncid);
  if(status != NC_NOERR) goto error;

  variable.destroy();
  return (0);

error:
  nc_check_error(status, __LINE__, __FILE__);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int poc_put_UG4D(const char *filename, mesh_t & mesh, int frame, int id, double **buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated, status;
  mesh_t global;
  double *gbuffer=0;
  int discretisation;
  extern int h_discretisation(const char *filename, int id);
  
  discretisation= h_discretisation (filename, id);

#ifdef DEBUG_PARALLEL
  printf("poc_put_UG4D (double), cpu=%d, frame=%d, id=%d, address=%d\n", gCPU_ID, frame, id, buffer);
#endif

  switch(gArchiveExtent) {
    case SEQUENTIAL_OUTPUT:
      status=poc_seqput_UG4D(filename, mesh, frame, id, buffer);
      break;
    case GLOBALIZED_OUTPUT:
      global=*(mesh.origin);
      if(gCPU_ID==gCPU_MASTER) status=poc_seqput_UG4D(filename, mesh, frame, id, buffer);
      if(allocated) delete[] gbuffer;
      break;
    case PARTITIONED_OUTPUT:
      status=poc_seqput_UG4D(filename, mesh, frame, id, buffer);
      break;
    default:
      check_error(-1, "illicit output mode", __LINE__, __FILE__, 1);
    }
    return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int poc_put_UG4D(const char *filename, mesh_t & mesh, int frame, int level, int id, float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated, status;
  mesh_t global;
  float *gbuffer=0;
  int discretisation;
  extern int h_discretisation(const char *filename, int id);
  
  discretisation= h_discretisation (filename, id);

#ifdef DEBUG_PARALLEL
  printf("poc_put_UG4D (double), cpu=%d, frame=%d, id=%d, address=%d\n", gCPU_ID, frame, id, buffer);
#endif

  switch(gArchiveExtent) {
    case SEQUENTIAL_OUTPUT:
      status=poc_seqput_UG4D(filename, mesh, frame, level, id, buffer);
      break;
    case GLOBALIZED_OUTPUT:
      global=*(mesh.origin);
//      allocated=OutputVector(buffer, &gbuffer, mesh, discretisation, gCPU_ID, 0);
      if(gCPU_ID==gCPU_MASTER) status=poc_seqput_UG4D(filename, mesh, frame, level, id, gbuffer);
      if(allocated) delete[] gbuffer;
      break;
    case PARTITIONED_OUTPUT:
      status=poc_seqput_UG4D(filename, mesh, frame, level, id, buffer);
      break;
    default:
      check_error(-1, "illicit output mode", __LINE__, __FILE__, 1);
    }
    return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int poc_put_UG4D(const char *filename, mesh_t & mesh, int frame, int level, int id, double *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated, status;
  mesh_t global;
  double *gbuffer=0;
  int discretisation;
  extern int h_discretisation(const char *filename, int id);
  
  discretisation= h_discretisation (filename, id);

#ifdef DEBUG_PARALLEL
  printf("poc_put_UG4D (double), cpu=%d, frame=%d, id=%d, address=%d\n", gCPU_ID, frame, id, buffer);
#endif

  switch(gArchiveExtent) {
    case SEQUENTIAL_OUTPUT:
      status=poc_seqput_UG4D(filename, mesh, frame, level, id, buffer);
      break;
    case GLOBALIZED_OUTPUT:
      global=*(mesh.origin);
//      allocated=OutputVector(buffer, &gbuffer, mesh, discretisation, gCPU_ID, 0);
      if(gCPU_ID==gCPU_MASTER) status=poc_seqput_UG4D(filename, mesh, frame, level, id, gbuffer);
      if(allocated) delete[] gbuffer;
      break;
    case PARTITIONED_OUTPUT:
      status=poc_seqput_UG4D(filename, mesh, frame, level, id, buffer);
      break;
    default:
      check_error(-1, "illicit output mode", __LINE__, __FILE__, 1);
    }
    return(status);
  }


/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check : MANDATORY !!!

  Note:

  19/02/2010

    the following in not really necessary

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_seqwrite_xyt(const char *filename, grid_t grid, int frame, int id,float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int ncid;     /* netCDF id */
  size_t start[3];
  size_t count[3];
  int status = 0;

  start[0] = frame;
  start[1] = 0;
  start[2] = 0;
  count[0] = 1;
  count[1] = grid.ny;
  count[2] = grid.nx;

  if(id==0) {
    printf("netcdf id is zero, generally used by non generic variables: alert!!!\n");
    }

  status = nc_open(filename,(int) NC_WRITE, &ncid);
  if(status != 0){
    nc_check_error(status, __LINE__, __FILE__);
    goto error;
    }

  status = nc_put_vara_float(ncid, id, start, count, buffer);
  if(status != 0) {
    nc_check_error(status, __LINE__, __FILE__);
    int chk = nc_close(ncid);
    goto error;
    }

  status = nc_close(ncid);
  if(status != 0) {
    nc_check_error(status, __LINE__, __FILE__);
    goto error;
    }

error:
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_seqwrite_xyt(const char *filename, grid_t grid, int frame, int id,double *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int ncid;     /* netCDF id */
  size_t start[3];
  size_t count[3];
  int status = 0;

  start[0] = frame;
  start[1] = 0;
  start[2] = 0;
  count[0] = 1;
  count[1] = grid.ny;
  count[2] = grid.nx;

  if(id==0) {
    printf("netcdf id is zero, generally used by non generic variables: alert!!!\n");
    }

  status = nc_open(filename,(int) NC_WRITE, &ncid);
  if(status != 0){
    nc_check_error(status, __LINE__, __FILE__);
    goto error;
    }

  status = nc_put_vara_double(ncid, id, start, count, buffer);
  if(status != 0) {
    int close_status = nc_close(ncid);
    nc_check_error(status, __LINE__, __FILE__);
    goto error;
    }

  status = nc_close(ncid);
  if(status != 0) {
    nc_check_error(status, __LINE__, __FILE__);
    goto error;
    }

error:
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_write_xyt(const char *filename, grid_t grid, int frame, int id, float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated=0, status;
  mesh_t global;
  float *gbuffer=0;

// #ifdef DEBUG_PARALLEL
//   printf("poc_write_xyt (float), frame=%d, id=%d, address=%d\n", frame, id, buffer);
// #endif
 
  switch(gArchiveExtent) {
    case SEQUENTIAL_OUTPUT:
      status=poc_seqwrite_xyt(filename, grid, frame, id, buffer);
      break;
    case GLOBALIZED_OUTPUT:
// #ifdef HAVE_MPI
//      status = MPI_Barrier(MPI_COMM_WORLD);
// #endif
      if(gCPU_ID==gCPU_MASTER) status=poc_seqwrite_xyt(filename, grid, frame, id, buffer);
      break;
    case PARTITIONED_OUTPUT:
      status=poc_seqwrite_xyt(filename, grid, frame, id, buffer);
      break;
    default:
      check_error(-1, "illicit output mode", __LINE__, __FILE__, 1);
    }
  
// #ifdef HAVE_MPI
//   status = MPI_Barrier(MPI_COMM_WORLD);
// #endif
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_write_xyt(const char *filename, grid_t grid, int frame, int id, double *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated=0, status;
  mesh_t global;
  float *gbuffer=0;

// #ifdef DEBUG_PARALLEL
//   printf("poc_write_xyt (float), frame=%d, id=%d, address=%d\n", frame, id, buffer);
// #endif
 
  switch(gArchiveExtent) {
    case SEQUENTIAL_OUTPUT:
      status=poc_seqwrite_xyt(filename, grid, frame, id, buffer);
      break;
    case GLOBALIZED_OUTPUT:
// #ifdef HAVE_MPI
//      status = MPI_Barrier(MPI_COMM_WORLD);
// #endif
      if(gCPU_ID==gCPU_MASTER) status=poc_seqwrite_xyt(filename, grid, frame, id, buffer);
      break;
    case PARTITIONED_OUTPUT:
      status=poc_seqwrite_xyt(filename, grid, frame, id, buffer);
      break;
    default:
      check_error(-1, "illicit output mode", __LINE__, __FILE__, 1);
    }
  
// #ifdef HAVE_MPI
//   status = MPI_Barrier(MPI_COMM_WORLD);
// #endif
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_UG3D(const char *filename, mesh_t & mesh, int frame, int id, double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int discretisation,allocated;
  int gdim;
  extern int h_discretisation(const char *filename, int id);
  double *gbuffer;
  mesh_t global;
  
  switch(gArchiveExtent) {
    case SEQUENTIAL_OUTPUT:
      status=poc_seqget_UG3D(filename, mesh, frame, id, z);
      break;
    case GLOBALIZED_OUTPUT:
/**----------------------------------------------------------------------------
      read discretisation so global output vector can be formed from local ones */
      discretisation= h_discretisation (filename, id);
      global=*(mesh.origin);
//      status=vector_dimension(mesh, discretisation, &gdim, CONTEXT_GLOBAL);
      gbuffer=new double[gdim];
      status=poc_seqget_UG3D(filename, mesh, frame, id, gbuffer);
//      status=VectorGlobal2local(mesh, discretisation, z, gbuffer);
      delete[] gbuffer;
      break;
    case PARTITIONED_OUTPUT:
      status=poc_seqget_UG3D(filename, mesh, frame, id, z);
      break;
    default:
      check_error(-1, "illicit output mode", __LINE__, __FILE__, 1);
    }
  
  return (status);
}

