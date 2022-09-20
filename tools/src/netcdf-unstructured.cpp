
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief File creation and variable routine, UG (some obsolete)
*/
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"
#include "constants.h"

#include "poc-netcdf.def"
#include "netcdf-proto.h"
#include "fe.def"

/* Variables externes --------------------------------------------------------*/

/* Variables internes --------------------------------------------------------*/

#ifdef _add_

#endif

#ifdef _add__

#endif

extern date_t time_getdate(double, date_t, char, double *);


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   int poc_put_UG2D(const char *filename, mesh_t & mesh, int id, float *z)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int k,n, N, recordsize;
//   int status;
//   size_t nitems;
//   size_t start[2];
//   size_t count[2],size[2];
//   int ncid;
//   cdfvar_t variable;
//   status = nc_open(filename, NC_WRITE, &ncid);
//   if(status != NC_NOERR) goto error;
//
//
//   status = cdf_varinfo(ncid,id, &variable);
//
//   for(k=0;k<variable.ndim;k++) {
//     size[k]=variable.dim[k].length;
//     }
//
//   start[0] = 0;
//   count[0] = size[0];
//
//   status = nc_put_vara_float(ncid, id, start, count, z);
//   if(status != NC_NOERR) goto error;
//
//   status = nc_close(ncid);
//   if(status != NC_NOERR) goto error;
//
// //  status= free_ncvariable(&variable);
//   variable.destroy();
//   return (0);
//
// error:
// //  nc_check_error(status, __LINE__, __FILE__);
//   return (status);
// }


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   int poc_put_UG3D(const char *filename, mesh_t & mesh, int frame, int id, double *z)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int k,n, N, recordsize;
//   int status;
//   size_t nitems;
//   size_t start[2];
//   size_t count[2],size[2];
//   int ncid;
//   cdfvar_t variable;
//   status = nc_open(filename, NC_WRITE, &ncid);
//   if(status != NC_NOERR) goto error;
//
//   status = cdf_varinfo(ncid,id, &variable );
//
//   for(k=0;k<variable.ndim;k++) {
//     size[k]=variable.dim[k].length;
//     }
//
//   start[0] = frame;
//   start[1] = 0;
//
//   count[0] = 1;
//   count[1] = size[1];
//
//   status = nc_put_vara_double(ncid, id, start, count, z);
//   if(status != NC_NOERR) goto error;
//
//   status = nc_close(ncid);
//   if(status != NC_NOERR) goto error;
//
//   variable.destroy();
//   return (0);
//
// error:
//   nc_check_error(status, __LINE__, __FILE__);
//   return (status);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_att(int ncid, int id , const char* name, short *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_get_att_short(ncid, id, name, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_att(int ncid, int id , const char* name, int *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_get_att_int(ncid, id, name, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_att(int ncid, int id , const char* name, float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_get_att_float(ncid, id, name, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_att(int ncid, int id , const char* name, double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=nc_get_att_double(ncid, id, name, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_seqget_UG3D(const char *filename, mesh_t & mesh, int frame, int id, double *z)

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

//
  status = cdf_varinfo(ncid,id, &variable ,0);

  for(k=0;k<variable.ndim;k++) {
    size[k]=variable.dim[k].length;
    }

  start[0] = frame;
  start[1] = 0;

  count[0] = 1;
  count[1] = size[1];

  status = nc_get_vara_double(ncid, id, start, count, z);
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

template<typename T> int poc_get_UG3D_template(int ncid, int frame, int level, int id, T *z, int *n=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k;
  int status;
  size_t *start,*count,*size;
  cdfvar_t variable;

  status = cdf_varinfo(ncid,id, &variable );

  start=new size_t[variable.ndim];
  count=new size_t[variable.ndim];
  size =new size_t[variable.ndim];
  
  for(k=0;k<variable.ndim;k++) {
    size[k]=variable.dim[k].length;
    }

  if(frame==LASTFRAME)
    frame=size[0]-1;
  
  switch (variable.ndim) {
    case 1:
      frame==NOFRAME;
      start[0] = 0;
      count[0] = size[0];
      if(n)
        *n=size[0];
      break;
    case 2:
      if(frame==NOFRAME) {
        start[0] = 0;
        count[0] = size[0];
        }
      else {
        start[0] = frame;
        count[0] = 1;
        }
      start[1] = 0;
      count[1] = size[1];
      if(n)
        *n=size[1];
      break;
    case 3:
      if(level==LASTLEVEL)
        level=size[1]-1;
      start[0] = frame;
      start[1] = level;
      start[2] = 0;
      count[0] = 1;
      count[1] = 1;
      count[2] = size[2];
      if(n)
        *n=size[2];
      break;
    }

  status = poc_get_vara(ncid, id, start, count, z);
  if(status != NC_NOERR) goto error;

  delete[] start;
  delete[] count;
  delete[] size;
  
  return (0);

error:
  nc_check_error(status, __LINE__, __FILE__);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_UG3D(int ncid, int frame, int id, int *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
///
/**
\note NetCDF automatically carries out the conversion between the recorded type and the required type!
*/
{
  int level=-1;
  return poc_get_UG3D_template(ncid, frame, level, id, z);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_UG3D(int ncid, int frame, int id, float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int level=-1;
  return poc_get_UG3D_template(ncid, frame, level, id, z);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_UG3D(int ncid, int frame, int id, double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int level=-1;
  return poc_get_UG3D_template(ncid, frame, level, id, z);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_get_UG3D_template(int ncid, int frame, const char *varname, T *z, int *n=NULL, int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, id;
  int level=-1;

  status = nc_inq_varid(ncid, varname, &id);
  if(status != NC_NOERR) goto error;
  
  status=poc_get_UG3D_template(ncid, frame, level, id, z, n);
  return (status);

error:
  if(verbose)NC_CHKERR_LINE_FILE(status,"nc_inq_varid(,\"%s\",) error",varname);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_UG3D(int ncid, int frame, const char *varname, int *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_get_UG3D_template(ncid, frame, varname, z, NULL, verbose);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_UG3D(int ncid, int frame, const char *varname, float *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_get_UG3D_template(ncid, frame, varname, z, NULL, verbose);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_UG3D(int ncid, int frame, const char *varname, double *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_get_UG3D_template(ncid, frame, varname, z, NULL, verbose);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_UG3D(int ncid, int frame, const char *varname, dcomplex *z, double  *aBuf, double *GBuf)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int aLen,GLen,i,status;
  char *bufName;
  asprintf(&bufName,"a_%s",varname);
  status=poc_get_UG3D_template(ncid, frame, bufName, aBuf, &aLen);
  free(bufName);
  if(status)return status;
  asprintf(&bufName,"G_%s",varname);
  status=poc_get_UG3D_template(ncid, frame, bufName, GBuf, &GLen);
  free(bufName);
  if(status)return status;
  if(aLen!=GLen)return NC_EVARSIZE;
  for(i=0;i<aLen;i++){
    z[i]=polar(aBuf[i],-GBuf[i]*d2r);
    }
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_get_UG3D_template(const char *filename, int frame, const char *varname, T *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, ncid,nc_status;

  nc_status = nc_open(filename, 0, &ncid);
  if(nc_status != NC_NOERR) goto error;
  
  status=poc_get_UG3D_template(ncid, frame, varname, z, NULL, verbose);
  
  nc_status=nc_close(ncid);
  if(nc_status != NC_NOERR) goto error;
  
  return (status);
  
error:
  verbose=1;
  if(verbose)NC_CHKERR_LINE_FILE(nc_status,"nc_open(\"%s\",,) error",filename);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_UG3D(const char *filename, int frame, const char *varname, int *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_get_UG3D_template(filename, frame, varname, z, verbose);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_UG3D(const char *filename, int frame, const char *varname, float *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_get_UG3D_template(filename, frame, varname, z, verbose);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_UG3D(const char *filename, int frame, const char *varname, double *z, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_get_UG3D_template(filename, frame, varname, z, verbose);
  
  return result;
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   template <typename T> int poc_get_UG3D_template(int ncid, int frame, int id, T *z,int *length)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int k,n, N, recordsize;
//   int status;
//   size_t nitems;
//   size_t start[2];
//   size_t count[2],size[2];
//   cdfvar_t variable;
//
//   status=init_ncvariable(&variable);
//   status = cdf_varinfo(ncid,id, &variable );
//
//   for(k=0;k<variable.ndim;k++) {
//     size[k]=variable.dim[k].length;
//     }
//
//   start[0] = frame;
//   start[1] = 0;
//
//   count[0] = 1;
//   count[1] = size[1];
//
//   status = nc_get_vara(ncid, id, start, count, z);
//   if(status != NC_NOERR) goto error;
//
//   status= free_ncvariable(&variable);
//
//   *length=size[1];
//
//   return (0);
//
// error:
//   nc_check_error(status, __LINE__, __FILE__);
//   return (status);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int poc_get_UG3D_template(const char *filename, int frame, int level, int id, T *z, int *length)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int ncid;
  int status, local;
  
  status = nc_open(filename, NC_NOWRITE, &ncid);
  if(status != NC_NOERR) goto error;

  status=poc_get_UG3D_template(ncid, frame, level, id, z, length);
  if(status != NC_NOERR) {
    local = nc_close(ncid);
    goto error;
    }

  status = nc_close(ncid);
  if(status != NC_NOERR) goto error;

  return (0);

error:
  nc_check_error(status, __LINE__, __FILE__);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T, typename Z> int poc_get_UG3D_template(T input, int frame, int level, const char *amplitude, const char *phase, complex<Z> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int dim,aLen,GLen,i,status,id;
  char   *bufName;
  Z *aBuf,*GBuf;
//  double *aBuf,*GBuf,radian;
  size_t count[2],size[2];
  cdfvar_t variable;
  cdfgbl_t global;
  
  id = cdf_varid(input, amplitude);
  if(id == -1) {
    printf("amplitude variable %s not found\n", amplitude);
    return (-1);
    }

  status = cdf_varinfo(input,id, &variable );
  status = cdf_globalinfo(input,&global,0);

  switch (variable.ndim) {
    case 1:
      dim=variable.dim[0].length;
      frame=NOFRAME;
      break;
    case 2:
      dim=variable.dim[1].length;
      break;
    case 3:
      dim=variable.dim[2].length;
      break;
    }
  aBuf=new Z[dim];
  GBuf=new Z[dim];
//   aBuf=new double[dim]; ///HERE
//   GBuf=new double[dim]; ///HERE
  
  status=poc_get_UG3D_template(input, frame, level, id, aBuf, &aLen);
  if(status)return status;

  id = cdf_varid(input, phase);
  if(id == -1) {
    printf("phase variable %s not found\n",phase);
    return (-1);
    }
    
  status=poc_get_UG3D_template(input, frame, level, id, GBuf, &GLen);
  if(status)return status;

  if(aLen!=GLen) {
    printf("amplitude variable size (%d) does not fir phase-lage variable size (%d)\n",aLen,GLen);
    return NC_EVARSIZE;
    }

#pragma omp parallel for
  for(i=0;i<aLen;i++){
    Z radian=-GBuf[i]*d2r;
    z[i]=(complex<Z> )polar(aBuf[i],radian);
    }
  delete[] aBuf;
  delete[] GBuf;
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_UG3D(int ncid, int frame, const char *amplitude, const char *phase, complex<double> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int level =-1;
  
  status = poc_get_UG3D_template(ncid, frame, level, amplitude, phase, z);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_UG3D(int ncid, int frame, const char *amplitude, const char *phase, complex<float> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int level =-1;
  
  status = poc_get_UG3D_template(ncid, frame, level, amplitude, phase, z);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_UG3D(int ncid, int frame, int level, const char *amplitude, const char *phase, complex<double> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status = poc_get_UG3D_template(ncid, frame, level, amplitude, phase, z);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_UG3D(int ncid, int frame, int level, const char *amplitude, const char *phase, complex<float> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status = poc_get_UG3D_template(ncid, frame, level, amplitude, phase, z);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_UG3D(const char *filename, int frame, const char *amplitude, const char *phase, complex<double> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int level =-1;
  
  status = poc_get_UG3D_template(filename, frame, level, amplitude, phase, z);
  if(status!=0) {
    printf("read error in %s\n",filename);
    }
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_UG3D(const char *filename, int frame, const char *amplitude, const char *phase, complex<float> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int level =-1;
  
  status = poc_get_UG3D_template(filename, frame, level, amplitude, phase, z);
  if(status!=0) {
    printf("poc_get_UG3D read error from %s (amplitude=%s, phase=%s, iteration=%d)\n",filename,amplitude, phase, frame);
    }
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_UG3D(const char *filename, int frame, int level, const char *amplitude, const char *phase, complex<double> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status = poc_get_UG3D_template(filename, frame, level, amplitude, phase, z);
  if(status!=0) {
    printf("read error in %s\n",filename);
    }
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_UG3D(const char *filename, int frame, int level, const char *amplitude, const char *phase, complex<float> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status = poc_get_UG3D_template(filename, frame, level, amplitude, phase, z);
  if(status!=0) {
    printf("poc_get_UG3D read error from %s (amplitude=%s, phase=%s, iteration=%d)\n",filename,amplitude, phase, frame);
    }
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int h_discretisation(const char *filename, int id)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, discretisation, status;
  cdfvar_t var;
  cdfgbl_t global;
  char *dimname;

  status=cdf_globalinfo(filename,&global,0);
  status=cdf_varinfo(filename,id,&var);
  
//   for(k=0;k<var.ndim;k++) {
//     status = nc_inq_dimlen(ncid, variable.dimids[k], &(size[k]));
//     status=nc_inq_dim(ncid,variable.dimids[k],name,&lengthhp);
//     }

  dimname=global.dimension[var.dim[var.ndim-1].id].name;

  if(strcmp(dimname,"M")==0) {
    discretisation= LGP0;
    }
  else if(strcmp(dimname,"NQUADRANGLES")==0) {
    discretisation= LGP0;
    }
  else if(strcmp(dimname,"N")==0) {
    discretisation= LGP1;
    }
  else if(strcmp(dimname,DGP1_standards[NNODES_DIMNAME])==0) {
    discretisation= DGP1;
    }
  else if(strcmp(dimname,"E")==0) {
    discretisation= NCP1;
    }
  else if(strcmp(dimname,DNP1_standards[NNODES_DIMNAME])==0) {
    discretisation= DNP1;
    }
  else if(strcmp(dimname,LGP2_standards[NNODES_DIMNAME])==0) {
    discretisation= LGP2;
    }
  else if(strcmp(dimname,CQP0_standards[NNODES_DIMNAME])==0) {
    discretisation= CQP0;
    }
  else if(strcmp(dimname,CQN1_standards[NNODES_DIMNAME])==0) {
    discretisation= CQN1;
    }
  else {
    TRAP_ERR_EXIT(-1, "discretisation not implemented");
    }

//  printf("h_discretisation, filename %s, id %d, var %s, dimname %s, discretisation %d\n",filename,id,var.name,dimname,discretisation);

  var.destroy();
  global.destroy();

  return(discretisation);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_UG4D(const char *filename, mesh_t mesh, int frame, int id, double **z)

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

  status = cdf_varinfo(ncid,id, &variable );

  for(k=0;k<variable.ndim;k++) {
    size[k]=variable.dim[k].length;
    }

  start[0] = frame;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = 1;
  count[2] = size[2];


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

  int poc_put_UG4D(const char *filename, mesh_t mesh, int frame, int id, float **z)

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

  status = cdf_varinfo(ncid,id, &variable );

  for(k=0;k<variable.ndim;k++) {
    size[k]=variable.dim[k].length;
    }

  start[0] = frame;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = 1;
  count[2] = size[2];


  for(k=0;k<size[1];k++) {
    start[1] = k;
    status = nc_put_vara_float(ncid, id, start, count, z[k]);
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

  int poc_put_UG4D(const char *filename, mesh_t mesh, int frame, int id, float *z)

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


  status = cdf_varinfo(ncid,id, &variable );

  for(k=0;k<variable.ndim;k++) {
    size[k]=variable.dim[k].length;
    }

  start[0] = frame;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = size[1];
  count[2] = size[2];


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

  int poc_get_UG4D(const char *filename, mesh_t & mesh, int frame, int id, double **z)

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


  for(k=0;k<size[1];k++) {
    start[1] = k;
    status = nc_get_vara_double(ncid, id, start, count, z[k]);
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

  int poc_get_UG4D(int ncid, mesh_t & mesh, int frame, int id, double **z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n, N, recordsize;
  int status;
  size_t nitems;
  size_t start[3];
  size_t count[3],size[3];
  cdfvar_t variable;

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


  for(k=0;k<size[1];k++) {
    start[1] = k;
    status = nc_get_vara_double(ncid, id, start, count, z[k]);
    if(status != NC_NOERR) goto error;
    }


  variable.destroy();
  return (0);

error:
  nc_check_error(status, __LINE__, __FILE__);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_UG4D(int ncid, int frame, const char *AA, const char *GG, complex<double> **z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l, n, N, recordsize;
  int status;
  size_t nitems;
  size_t start[3];
  size_t count[3],size[3];
  cdfvar_t aVar, GVar;
  int id_A, id_G;
  double *a,*G;
  const complex< double > J=(0.,1.);
  
  status = cdf_varinfo(ncid, AA, &aVar ,0);
  status = cdf_varinfo(ncid, GG, &GVar ,0);

  for(k=0;k<aVar.ndim;k++) {
    size[k]=aVar.dim[k].length;
    }

  start[0] = frame;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = 1;
  count[2] = size[2];

  a=new double[size[2]];
  G=new double[size[2]];

  for(k=0;k<size[1];k++) {
    start[1] = k;
    status = nc_get_vara_double(ncid, aVar.id, start, count, a);
    if(status != NC_NOERR) goto error;
    status = nc_get_vara_double(ncid, GVar.id, start, count, G);
    if(status != NC_NOERR) goto error;
//    for(l=0;l<size[2];l++) z[k][l]+=tmp[l];
    for(l=0;l<size[2];l++) z[l][k]=polar<double>(a[l],-G[l]*d2r);
    }

  delete[] a;
  delete[] G;
  
  aVar.destroy();
  GVar.destroy();
  return (0);

error:
  nc_check_error(status, __LINE__, __FILE__);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_UG4D(const char *filename, int frame, const char *A, const char *G, complex<double> **z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int ncid, status;
  
  status = nc_open(filename, NC_WRITE, &ncid);
  if(status != NC_NOERR) goto error;

  status=poc_get_UG4D(ncid, frame, A, G, z);
  
  status = nc_close(ncid);
  if(status != NC_NOERR) goto error;

  return (status);
error:
  nc_check_error(status, __LINE__, __FILE__);
  return (status);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_getLast_template(int ncid, int frame, int frequency, int id, T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k;
  int status;
  size_t *start,*count,*size;
  cdfvar_t variable;

  status = cdf_varinfo(ncid,id, &variable );

  start=new size_t[variable.ndim];
  count=new size_t[variable.ndim];
  size =new size_t[variable.ndim];
  
  for(k=0;k<variable.ndim;k++) {
    size[k]=variable.dim[k].length;
    start[k] = 0;
    count[k] = size[k];
    }

  if(frame<0)
    frame+=size[0];
  start[0] = frame;
  start[1] = frequency;

  count[0] = 1;
  count[1] = 1;

  status = poc_get_vara(ncid, id, start, count, z);
  if(status != NC_NOERR) goto error;

  delete[] start;
  delete[] count;
  delete[] size;
  
  return (0);

error:
  nc_check_error(status, __LINE__, __FILE__);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_putLast_template(int ncid, int frame, int frequency, int id, T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k;
  int status;
  size_t *start,*count,*size;
  cdfvar_t variable;

  status = cdf_varinfo(ncid,id, &variable );

  start=new size_t[variable.ndim];
  count=new size_t[variable.ndim];
  size =new size_t[variable.ndim];
  
  for(k=0;k<variable.ndim;k++) {
    size[k]=variable.dim[k].length;
    start[k] = 0;
    count[k] = size[k];
    }

  if(frame<0)
    frame+=size[0];
  start[0] = frame;
  start[1] = frequency;

  count[0] = 1;
  count[1] = 1;

  status = poc_put_vara(ncid, id, start, count, z);
  if(status != NC_NOERR) goto error;

  delete[] start;
  delete[] count;
  delete[] size;
  
  return (0);

error:
  nc_check_error(status, __LINE__, __FILE__);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_putLast_template(const char *filename, int frame, int frequency, int id, T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int ncid, status;
  status=nc_open(filename,NC_WRITE,&ncid);
  if(status != NC_NOERR)
    goto error;

  status=poc_putLast_template(ncid, frame, frequency, id, z);

  status = nc_close(ncid);
  if(status != NC_NOERR)
    goto error;
  
  return(status);
error:
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_getLast(int ncid, int frame, int frequency, int id, short *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_getLast_template(ncid, frame, frequency, id, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_getLast(int ncid, int frame, int frequency, int id, int *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_getLast_template(ncid, frame, frequency, id, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_getLast(int ncid, int frame, int frequency, int id, float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_getLast_template(ncid, frame, frequency, id, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_getLast(int ncid, int frame, int frequency, int id, double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_getLast_template(ncid, frame, frequency, id, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_putLast(int ncid, int frame, int frequency, int id, short *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_putLast_template(ncid, frame, frequency, id, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_putLast(int ncid, int frame, int frequency, int id, int *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_putLast_template(ncid, frame, frequency, id, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_putLast(int ncid, int frame, int frequency, int id, float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_putLast_template(ncid, frame, frequency, id, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_putLast(int ncid, int frame, int frequency, int id, double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_putLast_template(ncid, frame, frequency, id, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_putLast(const char *filename, int frame, int frequency, int id, short *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_putLast_template(filename, frame, frequency, id, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_putLast(const char *filename, int frame, int frequency, int id, int *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_putLast_template(filename, frame, frequency, id, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_putLast(const char *filename, int frame, int frequency, int id, float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_putLast_template(filename, frame, frequency, id, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_putLast(const char *filename, int frame, int frequency, int id, double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_putLast_template(filename, frame, frequency, id, z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  cdfvar_t poc_shortvariable_nt(char *name, float mask, char *units,
                              float scale, float offset, char *standardname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cdfvar_t variable;
  int k, l, status;
  short smask;
  int ncid;

  smask = 32768;


  variable.id = -1;
  variable.type = NC_SHORT;

  variable.initdim(2);

  variable.natt = 9;
  variable.att = new cdfatt_t[variable.natt];

  variable.name = new char[strlen(name)+1];
  strcpy(variable.name,name);

  k = 0;
  variable.dim[k].name = strdup("T");
  k++;
  variable.dim[k].name = strdup("N");

  k = 0;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("units");
  variable.att[k].data = strdup(units);
  variable.att[k].length = strlen(variable.att[k].data);
/*
  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);
*/
  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("standard_name");
  variable.att[k].data = strdup(standardname);
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("short_name");
  variable.att[k].data = strdup(name);
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name= strdup("axis");
  variable.att[k].name = strdup("content");
  variable.att[k].data = strdup("TN");
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("associate");
  variable.att[k].data = strdup("time lat lon");
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_SHORT;
  variable.att[k].name = strdup("missing_value");
  variable.att[k].data = new char[sizeof(short)];
  for(l = 0; l < sizeof(short); l++)
    variable.att[k].data[l] = (char) *((char *) &smask + l);
  variable.att[k].length = sizeof(short);

  k++;
  variable.att[k].type = NC_SHORT;
  variable.att[k].name = strdup("_FillValue");
  variable.att[k].data = new char[sizeof(short)];
  for(l = 0; l < sizeof(short); l++)
    variable.att[k].data[l] = (char) *((char *) &smask + l);
  variable.att[k].length = sizeof(short);

  k++;
  variable.att[k].type = NC_FLOAT;
  variable.att[k].name = strdup("scale_factor");
  variable.att[k].data = new char[sizeof(float)];
  for(l = 0; l < sizeof(float); l++)
    variable.att[k].data[l] = (char) *((char *) &scale + l);
  bcopy((char *) &scale, variable.att[k].data, sizeof(float));
  variable.att[k].length = sizeof(float);

  k++;
  variable.att[k].type = NC_FLOAT;
  variable.att[k].name = strdup("offset");
  variable.att[k].data = new char[sizeof(float)];
  for(l = 0; l < sizeof(float); l++)
    variable.att[k].data[l] = (char) *((char *) &offset + l);
  variable.att[k].length = sizeof(float);

  return (variable);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  cdfvar_t poc_floatvariable_nt(char *name, float mask, char *units,
                              float scale, float offset, char *standardname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  cdfvar_t variable;
  int k, l, status;
  int ncid;



  variable.id = -1;
  variable.type = NC_FLOAT;

  variable.initdim(2);

  variable.natt = 9;
  variable.att = new cdfatt_t[variable.natt];

  variable.name = new char[strlen(name)+1];
  strcpy(variable.name,name);

  k = 0;
  variable.dim[k].name = strdup("T");
  k++;
  variable.dim[k].name = strdup("N");

  k = 0;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("units");
  variable.att[k].data = strdup(units);
  variable.att[k].length = strlen(variable.att[k].data);
/*
  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);
*/
  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("standard_name");
  variable.att[k].data = strdup(standardname);
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("short_name");
  variable.att[k].data = strdup(name);
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name= strdup("axis");
  variable.att[k].name = strdup("content");
  variable.att[k].data = strdup("TN");
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("associate");
  variable.att[k].data = strdup("time lat lon");
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_FLOAT;
  variable.att[k].name = strdup("missing_value");
  variable.att[k].data = new char[sizeof(float)];
  for(l = 0; l < sizeof(float); l++)
    variable.att[k].data[l] = (char) *((char *) &mask + l);
  variable.att[k].length = sizeof(float);

  k++;
  variable.att[k].type = NC_FLOAT;
  variable.att[k].name = strdup("_FillValue");
  variable.att[k].data = new char[sizeof(float)];
  for(l = 0; l < sizeof(float); l++)
    variable.att[k].data[l] = (char) *((char *) &mask + l);
  variable.att[k].length = sizeof(float);

  k++;
  variable.att[k].type = NC_FLOAT;
  variable.att[k].name = strdup("scale_factor");
  variable.att[k].data = new char[sizeof(float)];
  for(l = 0; l < sizeof(float); l++)
    variable.att[k].data[l] = (char) *((char *) &scale + l);
  bcopy((char *) &scale, variable.att[k].data, sizeof(float));
  variable.att[k].length = sizeof(float);

  k++;
  variable.att[k].type = NC_FLOAT;
  variable.att[k].name = strdup("offset");
  variable.att[k].data = new char[sizeof(float)];
  for(l = 0; l < sizeof(float); l++)
    variable.att[k].data[l] = (char) *((char *) &offset + l);
  variable.att[k].length = sizeof(float);

  return (variable);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> cdfvar_t poc_variable_UG2D_template(const char *name, T mask, const char *units, T scale, T offset,const char *standardname, const char *dim0, nc_type type)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cdfvar_t variable;
  int k, l, status;
  char *dim[2];

  dim[0]=strdup(dim0);

  variable.id   = -1;
  variable.type = type;

  variable.initdim(1);

  variable.natt = 9;
  variable.att = new cdfatt_t[variable.natt];

  variable.name = new char[strlen(name)+1];
  strcpy(variable.name,name);

  k = 0;
  variable.dim[k].name = poc_strdup(dim[k]);

  k = 0;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = poc_strdup("units");
  variable.att[k].data = poc_strdup(units);
  variable.att[k].length = strlen(variable.att[k].data);
/*
  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);
*/
  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = poc_strdup("standard_name");
  variable.att[k].data = poc_strdup(standardname);
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = poc_strdup("short_name");
  variable.att[k].data = poc_strdup(name);
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = poc_strdup("content");
  variable.att[k].data = poc_strdup("N");
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = poc_strdup("associate");
  variable.att[k].data = poc_strdup("lat lon");
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = type;
  variable.att[k].name = poc_strdup("missing_value");
  variable.att[k].data = new char[sizeof(T)];
  for(l = 0; l < sizeof(T); l++)
    variable.att[k].data[l] = (char) *((char *) &mask + l);
  variable.att[k].length = sizeof(T);

  k++;
  variable.att[k].type = type;
  variable.att[k].name = poc_strdup("_FillValue");
  variable.att[k].data = new char[sizeof(T)];
  for(l = 0; l < sizeof(T); l++)
    variable.att[k].data[l] = (char) *((char *) &mask + l);
  variable.att[k].length = sizeof(T);

  k++;
  variable.att[k].type = type;
  variable.att[k].name = poc_strdup("scale_factor");
  variable.att[k].data = new char[sizeof(T)];
  bcopy((char *) &scale, variable.att[k].data, sizeof(T));
  variable.att[k].length = sizeof(T);

  k++;
  variable.att[k].type = type;
  variable.att[k].name = poc_strdup("offset");
  variable.att[k].data = new char[sizeof(T)];
  bcopy((char *) &offset, variable.att[k].data, sizeof(T));
  variable.att[k].length = sizeof(T);

  for(k=0;k<1;k++) free(dim[k]);

  return (variable);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  cdfvar_t poc_variable_UG2D(const char *name, double mask, const char *units,
                             double scale, double offset,const char *standardname, const char *dim0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cdfvar_t variable;
  variable=poc_variable_UG2D_template(name, mask,  units, scale,  offset,  standardname,   dim0, NC_DOUBLE);
  return (variable);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  cdfvar_t poc_variable_UG2D(const char *name, float mask, const char *units,
                             float scale, float offset,const char *standardname, const char *dim0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cdfvar_t variable;
  variable << poc_variable_UG2D_template(name, mask,  units, scale,  offset,  standardname,   dim0, NC_FLOAT);
  return (variable);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   cdfvar_t poc_variable_UG2D(const char *name, float mask, const char *units,
//                              float scale, float offset,const char *standardname, const char *dim0)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   cdfvar_t variable;
//   int k, l, status;
//   char *dim[2];
//
//   dim[0]=strdup(dim0);
//
//
//
//   variable.id = -1;
//   variable.type = NC_FLOAT;
//
//   variable.initdim(1);
//
//   variable.natt = 9;
//   variable.att = new cdfatt_t[variable.natt];
//
//   variable.name = new char[strlen(name)+1];
//   strcpy(variable.name,name);
//
//   k = 0;
//   variable.dim[k].name = strdup(dim[k]);
//
//   k = 0;
//   variable.att[k].type = NC_CHAR;
//   variable.att[k].name = strdup("units");
//   variable.att[k].data = strdup(units);
//   variable.att[k].length = strlen(variable.att[k].data);
// /*
//   k++;
//   variable.att[k].type=NC_CHAR;
//   variable.att[k].name=strdup("long_name");
//   variable.att[k].data=strdup(longname);
//   variable.att[k].length=strlen( variable.att[k].data);
// */
//   k++;
//   variable.att[k].type = NC_CHAR;
//   variable.att[k].name = strdup("standard_name");
//   variable.att[k].data = strdup(standardname);
//   variable.att[k].length = strlen(variable.att[k].data);
//
//   k++;
//   variable.att[k].type = NC_CHAR;
//   variable.att[k].name = strdup("short_name");
//   variable.att[k].data = strdup(name);
//   variable.att[k].length = strlen(variable.att[k].data);
//
//   k++;
//   variable.att[k].type = NC_CHAR;
//   variable.att[k].name = strdup("content");
//   variable.att[k].data = strdup("N");
//   variable.att[k].length = strlen(variable.att[k].data);
//
//   k++;
//   variable.att[k].type = NC_CHAR;
//   variable.att[k].name = strdup("associate");
//   variable.att[k].data = strdup("lat lon");
//   variable.att[k].length = strlen(variable.att[k].data);
//
//   k++;
//   variable.att[k].type = NC_FLOAT;
//   variable.att[k].name = strdup("missing_value");
//   variable.att[k].data = new char[sizeof(float)];
//   for(l = 0; l < sizeof(float); l++)
//     variable.att[k].data[l] = (char) *((char *) &mask + l);
//   variable.att[k].length = sizeof(float);
//
//   k++;
//   variable.att[k].type = NC_FLOAT;
//   variable.att[k].name = strdup("_FillValue");
//   variable.att[k].data = new char[sizeof(float)];
//   for(l = 0; l < sizeof(float); l++)
//     variable.att[k].data[l] = (char) *((char *) &mask + l);
//   variable.att[k].length = sizeof(float);
//
//   k++;
//   variable.att[k].type = NC_FLOAT;
//   variable.att[k].name = strdup("scale_factor");
//   variable.att[k].data = new char[sizeof(float)];
//   bcopy((char *) &scale, variable.att[k].data, sizeof(float));
//   variable.att[k].length = sizeof(float);
//
//   k++;
//   variable.att[k].type = NC_FLOAT;
//   variable.att[k].name = strdup("offset");
//   variable.att[k].data = new char[sizeof(float)];
//   bcopy((char *) &offset, variable.att[k].data, sizeof(float));
//   variable.att[k].length = sizeof(float);
//
//   for(k=0;k<1;k++) free(dim[k]);
//
//   return (variable);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  cdfvar_t poc_variable_UG3D(const char *name, int mask,const char *units, int scale, int offset,const char *standardname, const char *dim0,const char *dim1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cdfvar_t variable;
  int k, l, status;
  char *dim[2];

  dim[0]=strdup(dim0);
  dim[1]=strdup(dim1);



  variable.id = -1;
  variable.type = NC_INT;

  variable.initdim(2);

  variable.natt = 9;
  variable.att = new cdfatt_t[variable.natt];

  variable.name = new char[strlen(name)+1];
  strcpy(variable.name,name);

  k = 0;
  variable.dim[k].name = (char *) poc_strdup(dim[k]);
  k++;
  variable.dim[k].name = (char *) poc_strdup(dim[k]);

  k = 0;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = (char *) poc_strdup("units");
  variable.att[k].data = (char *) poc_strdup(units);
  variable.att[k].length = strlen(variable.att[k].data);
/*
  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=(char *) poc_strdup("long_name");
  variable.att[k].data=(char *) poc_strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);
*/
  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = (char *) poc_strdup("standard_name");
  variable.att[k].data = (char *) poc_strdup(standardname);
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = (char *) poc_strdup("short_name");
  variable.att[k].data = (char *) poc_strdup(name);
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = (char *) poc_strdup("content");
  variable.att[k].data = (char *) poc_strdup("TN");
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = (char *) poc_strdup("associate");
  variable.att[k].data = (char *) poc_strdup("time level lat lon");
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_INT;
  variable.att[k].name = (char *) poc_strdup("missing_value");
  variable.att[k].data = new char[sizeof(double)];
  for(l = 0; l < sizeof(double); l++)
    variable.att[k].data[l] = (char) *((char *) &mask + l);
  variable.att[k].length = sizeof(double);

  k++;
  variable.att[k].type = NC_INT;
  variable.att[k].name = (char *) poc_strdup("_FillValue");
  variable.att[k].data = new char[sizeof(double)];
  for(l = 0; l < sizeof(double); l++)
    variable.att[k].data[l] = (char) *((char *) &mask + l);
  variable.att[k].length = sizeof(double);

  k++;
  variable.att[k].type = NC_INT;
  variable.att[k].name = (char *) poc_strdup("scale_factor");
  variable.att[k].data = new char[sizeof(double)];
  bcopy((char *) &scale, variable.att[k].data, sizeof(double));
  variable.att[k].length = sizeof(double);

  k++;
  variable.att[k].type = NC_INT;
  variable.att[k].name = (char *) poc_strdup("offset");
  variable.att[k].data = new char[sizeof(double)];
  bcopy((char *) &offset, variable.att[k].data, sizeof(double));
  variable.att[k].length = sizeof(double);

  for(k=0;k<2;k++) free(dim[k]);

  return (variable);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  cdfvar_t poc_variable_UG3D(const char *name, float mask, const char *units, float scale, float offset,const char *standardname,  const char *dim0, const char *dim1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cdfvar_t variable;
  int k, l, status;
  char *dim[2];

  dim[0]=strdup(dim0);
  dim[1]=strdup(dim1);



  variable.id = -1;
  variable.type = NC_FLOAT;

  variable.initdim(2);

  variable.natt = 9;
  variable.att = new cdfatt_t[variable.natt];

  variable.name = new char[strlen(name)+1];
  strcpy(variable.name,name);

  k = 0;
  variable.dim[k].name = (char *) poc_strdup(dim[k]);
  k++;
  variable.dim[k].name = (char *) poc_strdup(dim[k]);

  k = 0;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = (char *) poc_strdup("units");
  variable.att[k].data = (char *) poc_strdup(units);
  variable.att[k].length = strlen(variable.att[k].data);
/*
  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=(char *) poc_strdup("long_name");
  variable.att[k].data=(char *) poc_strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);
*/
  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = (char *) poc_strdup("standard_name");
  variable.att[k].data = (char *) poc_strdup(standardname);
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = (char *) poc_strdup("short_name");
  variable.att[k].data = (char *) poc_strdup(name);
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = (char *) poc_strdup("content");
  variable.att[k].data = (char *) poc_strdup("TN");
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = (char *) poc_strdup("associate");
  variable.att[k].data = (char *) poc_strdup("time lat lon");
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_FLOAT;
  variable.att[k].name = (char *) poc_strdup("missing_value");
  variable.att[k].data = new char[sizeof(float)];
  for(l = 0; l < sizeof(float); l++)
    variable.att[k].data[l] = (char) *((char *) &mask + l);
  variable.att[k].length = sizeof(float);

  k++;
  variable.att[k].type = NC_FLOAT;
  variable.att[k].name = (char *) poc_strdup("_FillValue");
  variable.att[k].data = new char[sizeof(float)];
  for(l = 0; l < sizeof(float); l++)
    variable.att[k].data[l] = (char) *((char *) &mask + l);
  variable.att[k].length = sizeof(float);

  k++;
  variable.att[k].type = NC_FLOAT;
  variable.att[k].name = (char *) poc_strdup("scale_factor");
  variable.att[k].data = new char[sizeof(float)];
  bcopy((char *) &scale, variable.att[k].data, sizeof(float));
  variable.att[k].length = sizeof(float);

  k++;
  variable.att[k].type = NC_FLOAT;
  variable.att[k].name = (char *) poc_strdup("offset");
  variable.att[k].data = new char[sizeof(float)];
  bcopy((char *) &offset, variable.att[k].data, sizeof(float));
  variable.att[k].length = sizeof(float);

  for(k=0;k<2;k++) free(dim[k]);

  return (variable);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  cdfvar_t poc_variable_UG3D(const char *name, double mask, const char *units,
                           double scale, double offset,const char *standardname,
                           const char *dim0, const char *dim1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cdfvar_t variable;
  int k, l, status;
  char *dim[2];

  dim[0]=strdup(dim0);
  dim[1]=strdup(dim1);



  variable.id = -1;
  variable.type = NC_DOUBLE;

  variable.initdim(2);

  variable.natt = 9;
  variable.att = new cdfatt_t[variable.natt];

  variable.name = new char[strlen(name)+1];
  strcpy(variable.name,name);

  k = 0;
  variable.dim[k].name = strdup(dim[k]);
  k++;
  variable.dim[k].name = strdup(dim[k]);

  k = 0;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("units");
  variable.att[k].data = strdup(units);
  variable.att[k].length = strlen(variable.att[k].data);
/*
  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);
*/
  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("standard_name");
  variable.att[k].data = strdup(standardname);
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("short_name");
  variable.att[k].data = strdup(name);
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name= strdup("axis");
  variable.att[k].name = strdup("content");
  variable.att[k].data = strdup("TN");
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("associate");
  variable.att[k].data = strdup("time level lat lon");
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_DOUBLE;
  variable.att[k].name = strdup("missing_value");
  variable.att[k].data = new char[sizeof(double)];
  for(l = 0; l < sizeof(double); l++)
    variable.att[k].data[l] = (char) *((char *) &mask + l);
  variable.att[k].length = sizeof(double);

  k++;
  variable.att[k].type = NC_DOUBLE;
  variable.att[k].name = strdup("_FillValue");
  variable.att[k].data = new char[sizeof(double)];
  for(l = 0; l < sizeof(double); l++)
    variable.att[k].data[l] = (char) *((char *) &mask + l);
  variable.att[k].length = sizeof(double);

  k++;
  variable.att[k].type = NC_DOUBLE;
  variable.att[k].name = strdup("scale_factor");
  variable.att[k].data = new char[sizeof(double)];
  bcopy((char *) &scale, variable.att[k].data, sizeof(double));
  variable.att[k].length = sizeof(double);

  k++;
  variable.att[k].type = NC_DOUBLE;
  variable.att[k].name = strdup("offset");
  variable.att[k].data = new char[sizeof(double)];
  bcopy((char *) &offset, variable.att[k].data, sizeof(double));
  variable.att[k].length = sizeof(double);

  for(k=0;k<2;k++) free(dim[k]);

  return (variable);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> cdfvar_t poc_variable_UG3D_template(const char *name, T mask, const char *units, T scale, T offset,const char *standardname, const char *dim0, int h_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cdfvar_t var;
  char *dim1;

  switch(h_discretisation) {

    case LGP0:
      dim1=strdup("M");
      break;

    case LGP1:
      dim1=strdup("N");
      break;

    case DGP1:
      dim1=strdup(DGP1_standards[NNODES_DIMNAME]);
      break;

    case NCP1:
      dim1=strdup("E");
      break;

    case DNP1:
      dim1=strdup(DNP1_standards[NNODES_DIMNAME]);
      break;

    case LGP2:
      dim1=strdup(LGP2_standards[NNODES_DIMNAME]);
      break;

    default:
     TRAP_ERR_EXIT(-1, "discretisation not implemented");
     break;
    }

  var=poc_variable_UG3D(name, mask, units, scale, offset, standardname, dim0, (const char *) dim1);

  free(dim1);

  return(var);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  cdfvar_t poc_variable_UG3D(const char *name, int mask, const char *units,
                           int scale, int offset,const char *standardname,
                           const char *dim0, int h_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cdfvar_t var;
  var=poc_variable_UG3D_template(name, mask, units, scale, offset, standardname, dim0, h_discretisation);
  return(var);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  cdfvar_t poc_variable_UG3D(const char *name, float mask, const char *units,
                           float scale, float offset,const char *standardname,
                           const char *dim0, int h_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cdfvar_t var;
  var=poc_variable_UG3D_template(name, mask, units, scale, offset, standardname, dim0, h_discretisation);
  return(var);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  cdfvar_t poc_variable_UG3D(const char *name, double mask, const char *units,
                           double scale, double offset,const char *standardname,
                           const char *dim0, int h_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cdfvar_t var;
  var=poc_variable_UG3D_template(name, mask, units, scale, offset, standardname, dim0, h_discretisation);
  return(var);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  cdfvar_t poc_variable_UG4D(const char *name, double mask, const char *units,
                           double scale, double offset,const char *standardname,
                           const char *dim0, const char *dim1, const char *dim2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cdfvar_t variable;
  int k, l, status;
  char *dim[3];

  dim[0]=strdup(dim0);
  dim[1]=strdup(dim1);
  dim[2]=strdup(dim2);



  variable.id = -1;
  variable.type = NC_DOUBLE;

  variable.initdim(3);

  variable.natt = 9;
  variable.att = new cdfatt_t[variable.natt];

  variable.name = new char[strlen(name)+1];
  strcpy(variable.name,name);

  k = 0;
  variable.dim[k].name = strdup(dim[k]);
  k++;
  variable.dim[k].name = strdup(dim[k]);
  k++;
  variable.dim[k].name = strdup(dim[k]);

  k = 0;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("units");
  variable.att[k].data = strdup(units);
  variable.att[k].length = strlen(variable.att[k].data);
/*
  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);
*/
  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("standard_name");
  variable.att[k].data = strdup(standardname);
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("short_name");
  variable.att[k].data = strdup(name);
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name= strdup("axis");
  variable.att[k].name = strdup("content");
  switch (dim[1][0]) {
    case 'N':
      variable.att[k].data = strdup("TNZ");
      break;
    case 'E':
      variable.att[k].data = strdup("TEZ");
      break;
    }
//  variable.att[k].data = strdup("TNZ");
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("associate");
  switch (dim[2][0]) {
    case 'L':
      variable.att[k].data = strdup("time level lat lon");
      break;
    case 'K':
      variable.att[k].data = strdup("time layer lat lon");
      break;
    }
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_DOUBLE;
  variable.att[k].name = strdup("missing_value");
  variable.att[k].data = new char[sizeof(double)];
  for(l = 0; l < sizeof(double); l++)
    variable.att[k].data[l] = (char) *((char *) &mask + l);
  variable.att[k].length = sizeof(double);

  k++;
  variable.att[k].type = NC_DOUBLE;
  variable.att[k].name = strdup("_FillValue");
  variable.att[k].data = new char[sizeof(double)];
  for(l = 0; l < sizeof(double); l++)
    variable.att[k].data[l] = (char) *((char *) &mask + l);
  variable.att[k].length = sizeof(double);

  k++;
  variable.att[k].type = NC_DOUBLE;
  variable.att[k].name = strdup("scale_factor");
  variable.att[k].data = new char[sizeof(double)];
  bcopy((char *) &scale, variable.att[k].data, sizeof(double));
  variable.att[k].length = sizeof(double);

  k++;
  variable.att[k].type = NC_DOUBLE;
  variable.att[k].name = strdup("offset");
  variable.att[k].data = new char[sizeof(double)];
  bcopy((char *) &offset, variable.att[k].data, sizeof(double));
  variable.att[k].length = sizeof(double);

  for(k=0;k<3;k++) free(dim[k]);

  return (variable);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  cdfvar_t poc_variable_UG4D(const char *name, float mask, const char *units,
                           float scale, float offset,const char *standardname,
                           const char *dim0, const char *dim1, const char *dim2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cdfvar_t variable;
  int k, l, status;
  char *dim[3];

  dim[0]=strdup(dim0);
  dim[1]=strdup(dim1);
  dim[2]=strdup(dim2);



  variable.id = -1;
  variable.type = NC_FLOAT;

  variable.initdim(3);

  variable.natt = 9;
  variable.att = new cdfatt_t[variable.natt];

  variable.name = new char[strlen(name)+1];
  strcpy(variable.name,name);

  k = 0;
  variable.dim[k].name = strdup(dim[k]);
  k++;
  variable.dim[k].name = strdup(dim[k]);
  k++;
  variable.dim[k].name = strdup(dim[k]);

  k = 0;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("units");
  variable.att[k].data = strdup(units);
  variable.att[k].length = strlen(variable.att[k].data);
/*
  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);
*/
  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("standard_name");
  variable.att[k].data = strdup(standardname);
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("short_name");
  variable.att[k].data = strdup(name);
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name= strdup("axis");
  variable.att[k].name = strdup("content");
  switch (dim[1][0]) {
    case 'K':
      variable.att[k].data = strdup("TZN");
      break;
    case 'N':
      variable.att[k].data = strdup("TNZ");
      break;
    case 'E':
      variable.att[k].data = strdup("TEZ");
      break;
    }
//  variable.att[k].data = strdup("TNZ");
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_CHAR;
  variable.att[k].name = strdup("associate");
  switch (dim[2][0]) {
    case 'L':
      variable.att[k].data = strdup("time level lat lon");
      break;
    case 'K':
      variable.att[k].data = strdup("time layer lat lon");
      break;
    case 'N':
      variable.att[k].data = strdup("time layer lat lon");
      break;
    }
  variable.att[k].length = strlen(variable.att[k].data);

  k++;
  variable.att[k].type = NC_FLOAT;
  variable.att[k].name = strdup("missing_value");
  variable.att[k].data = new char[sizeof(float)];
  for(l = 0; l < sizeof(float); l++)
    variable.att[k].data[l] = (char) *((char *) &mask + l);
  variable.att[k].length = sizeof(float);

  k++;
  variable.att[k].type = NC_FLOAT;
  variable.att[k].name = strdup("_FillValue");
  variable.att[k].data = new char[sizeof(float)];
  for(l = 0; l < sizeof(float); l++)
    variable.att[k].data[l] = (char) *((char *) &mask + l);
  variable.att[k].length = sizeof(float);

  k++;
  variable.att[k].type = NC_FLOAT;
  variable.att[k].name = strdup("scale_factor");
  variable.att[k].data = new char[sizeof(float)];
  bcopy((char *) &scale, variable.att[k].data, sizeof(float));
  variable.att[k].length = sizeof(float);

  k++;
  variable.att[k].type = NC_FLOAT;
  variable.att[k].name = strdup("offset");
  variable.att[k].data = new char[sizeof(float)];
  bcopy((char *) &offset, variable.att[k].data, sizeof(float));
  variable.att[k].length = sizeof(float);

  for(k=0;k<3;k++) free(dim[k]);

  return (variable);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  cdfvar_t poc_variable_UG4D(const char *name, double mask, const char *units,
                           double scale, double offset,const char *standardname,
                           const char *dim0, int h_discretisation, int v_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cdfvar_t var;
  char *dim1, *dim2;

  switch(h_discretisation) {

    case LGP0:
      dim1=strdup("M");
      break;

    case LGP1:
      dim1=strdup("N");
      break;

    case DGP1:
      dim1=strdup(DGP1_standards[NNODES_DIMNAME]);
      break;

    case NCP1:
      dim1=strdup("E");
      break;

    case DNP1:
      dim1=strdup(DNP1_standards[NNODES_DIMNAME]);
      break;

    case LGP2:
      dim1=strdup(LGP2_standards[NNODES_DIMNAME]);
      break;

    default:
      TRAP_ERR_EXIT(-1, "discretisation not implemented");
      break;
    }

  switch(v_discretisation) {

    case LAYERS:
      dim2=strdup("K");
      break;

    case LEVELS:
      dim2=strdup("L");
      break;

    default:
      TRAP_ERR_EXIT(-1, "discretisation not implemented");
      break;
    }
    
  var=poc_variable_UG4D(name, mask, units, scale, offset, standardname, dim0, (const char *) dim2, (const char *) dim1);

  free(dim1);
  free(dim2);

  return(var);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_read2d(int ncid,int varid, int layer, int frame, float **buffer, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,axis_id,n,status;
  int verbose=0;
  unsigned int size;
  float scale,offset;
  char *axis;
  int understood=0,pos;
  double mean,rms;
  size_t start[10],count[10],stride[10],imap[10];
  UGdecoded_t decoded;

  cdfgbl_t g_info;
  cdfvar_t v_info;

  date_t origine;
  double *time;
  int nframes;

  status= cdf_globalinfo(ncid,&g_info,verbose);
  status= cdf_varinfo(ncid,varid,&v_info);

parse:
  status=poc_UGdecode_axis( v_info, g_info, &decoded);
   if(status==0) {
    axis=strdup(decoded.axis);
    if(strcmp(axis,"N")==0) {
      start[0]=0;
      count[0]=v_info.dim[0].length;
      size=count[0];
      understood=1;
      }
    if(strcmp(axis,"TN")==0) {
      start[0]=frame;
      start[1]=0;
      count[0]=1;
      count[1]=v_info.dim[1].length;
      size=count[1];
      understood=1;
      }
    if(strcmp(axis,"NN")==0) {
      start[0]=0;
      start[1]=0;
      count[0]=v_info.dim[0].length;
      count[1]=v_info.dim[1].length;
      size=count[0]*count[1];
      understood=1;
      }
    if(strcmp(axis,"ZN")==0) {
      start[0]=layer;
      start[1]=0;
      count[0]=1;
      count[1]=v_info.dim[1].length;
      size=count[1];
      understood=1;
      }
    if(strcmp(axis,"NZ")==0) {
      start[0]=0;
      start[1]=layer;
      count[0]=v_info.dim[0].length;
      count[1]=1;
      size=count[0];
      understood=1;
      }
    if(strcmp(axis,"MP")==0) {
      start[0]=0;
      start[1]=0;
      count[0]=v_info.dim[0].length;
      count[1]=v_info.dim[1].length;
      size=count[0]*count[1];
      understood=1;
      }
    if(strcmp(axis,"TZN")==0) {
      start[0]=frame;
      start[1]=layer;
      start[2]=0;
      count[0]=1;
      count[1]=1;
      count[2]=v_info.dim[2].length;
      size=count[2];
      understood=1;
      }
    if(strcmp(axis,"TNZ")==0) {
      start[0]=frame;
      start[1]=0;
      start[2]=layer;
      count[0]=1;
      count[1]=v_info.dim[1].length;
      count[2]=1;
      size=count[1];
      understood=1;
      }
    free(axis);
    }
  else {
    status=-1;
    goto error;
    }

  if(understood==0) {
    status=-1;
    printf("unstructured array arrangement is not understood: %s\n",axis);
    goto error;
    }

  status= poc_UGdecode_discretisation(v_info, g_info, &decoded);
  if (status!=0) {
    printf("unstructured array discretisation is not reckognized\n");
    goto error;
    }

  if(*buffer!=0) delete[] *buffer;
  *buffer=new float[size];
/*------------------------------------------------------------------------------
  read  variable frame in netcdf file */
  status=nc_get_vara_float(ncid,varid,start,count,(*buffer));
//  status=ncvargetg(ncid, varid, start, count, stride, NULL, buffer);

  status=nc_get_att_float(ncid,varid,"scale_factor",&scale);
  if(status != NC_NOERR) scale=1.0;

  status=nc_get_att_float(ncid,varid,"add_offset",&offset);
  if(status != NC_NOERR) offset=0.0;

  status=nc_get_att_float(ncid,varid,"_FillValue",mask);
  if(status != NC_NOERR) *mask=1.0e+10;

  if((scale !=1.0) || (offset !=0.0)) {
    printf("scale and offset correction applies: scale=%lf offset=%lf\n", scale, offset);
    for(m=0;m<size;m++)
      if((*buffer)[m]!=*mask)  (*buffer)[m]=(*buffer)[m]*scale+offset;
    }

  n=0;
  mean=0;
  rms=0;
  for(m=0;m<size;m++)
    if((*buffer)[m]!=*mask) {
      mean+=(*buffer)[m];
      rms+=(*buffer)[m]*(*buffer)[m];
      n++;
      }

  if(n!=0) {
    mean/=(double) n;
    rms=sqrt(rms/(double)n-mean*mean);
    }
  else {
    status=-1;
    goto error;
    }

  status=0;
  return(status);

error:
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_read3d(int ncid,int varid, int frame, float **data, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,m,axis_id,n,status;
  int verbose=0;
  unsigned int size,hsize,vsize;
  float scale,offset,*tmp,*buffer;
  char *axis;
  int understood=0,rearrange=0;
  double mean,rms;
  size_t start[10],count[10];

  cdfgbl_t g_info;
  cdfvar_t v_info;
  UGdecoded_t decoded;

  date_t origine;
  double *time;
  int nframes;

  status= cdf_globalinfo(ncid,&g_info,verbose);
  if(status!=0) goto error;

  status= cdf_varinfo(ncid,varid,&v_info);
  if(status!=0) goto error;

/*------------------------------------------------------------------------------
  Find dimensions order: use "AXIS" attribute if available */
  axis_id=-1;
  for(n=0;n<v_info.natt;n++) {
    if(strcmp(v_info.att[n].name,"axis")==0) {
      axis_id=n;
      break;
      }
/* *------------------------------------------------------------------------
    attribute name changed to comply with CF standard*/
    if(strcmp(v_info.att[n].name,"content")==0) {
      axis_id=n;
      break;
      }
    }

  status=poc_UGdecode_axis( v_info, g_info, &decoded);
  if(status==0) {
    axis=strdup(decoded.axis);
    if(strcmp(axis,"TN")==0) {
      understood=0;
      }
    if(strcmp(axis,"ZN")==0) {
      start[0]=0;
      start[1]=0;
      count[0]=v_info.dim[0].length;
      count[1]=v_info.dim[1].length;
      size=count[1]*count[0];
      understood=1;
//      details->arrangement=PER_LAYER;
      }
    if(strcmp(axis,"TZN")==0) {
      start[0]=frame;
      start[1]=0;
      start[2]=0;
      count[0]=1;
      count[1]=v_info.dim[1].length;
      count[2]=v_info.dim[2].length;
      size=count[2]*count[1];
      understood=1;
//      details->arrangement=PER_LAYER;
      }
    if(strcmp(axis,"TZM")==0) {
      start[0]=frame;
      start[1]=0;
      start[2]=0;
      count[0]=1;
      count[1]=v_info.dim[1].length;
      count[2]=v_info.dim[2].length;
      size=count[2]*count[1];
      understood=1;
//      details->arrangement=PER_LAYER;
      }
    if(strcmp(axis,"TNZ")==0) {
      start[0]=frame;
      start[1]=0;
      start[2]=0;
      count[0]=1;
      count[1]=v_info.dim[1].length;
      count[2]=v_info.dim[2].length;
      size=count[2]*count[1];
      understood=1;
      hsize=count[1];
      vsize=count[2];
      rearrange=1;
//      details->arrangement=PER_COLUMN;
      }
    if(strcmp(axis,"TMZ")==0) {
      start[0]=frame;
      start[1]=0;
      start[2]=0;
      count[0]=1;
      count[1]=v_info.dim[1].length;
      count[2]=v_info.dim[2].length;
      size=count[2]*count[1];
      understood=1;
//      details->arrangement=PER_COLUMN;
      }
    free(axis);
    }
  else {
    status=-1;
    goto error;
    }

  if(understood==0) {
    status=-1;
    goto error;
    }

  status= poc_UGdecode_discretisation(v_info, g_info, &decoded);
  if (status!=0) goto error;

/*------------------------------------------------------------------------------
  read  variable frame in netcdf file */
  if(*data==0) buffer=new float[size];
  else buffer=*data;
  status=nc_get_vara_float(ncid,varid,start,count,buffer);

  status=nc_get_att_float(ncid,varid,"scale_factor",&scale);
  if(status != NC_NOERR) scale=1.0;

  status=nc_get_att_float(ncid,varid,"add_offset",&offset);
  if(status != NC_NOERR) offset=0.0;

  status=nc_get_att_float(ncid,varid,"_FillValue",mask);
  if(status != NC_NOERR) *mask=1.0e+10;

  if((scale !=1.0) || (offset !=0.0))
    for(m=0;m<size;m++)
      if(buffer[m]!=*mask)  buffer[m]=buffer[m]*scale+offset;

  if(rearrange==1) {
/* *------------------------------------------------------------------
    change NZ to ZN order*/
    tmp=(float *) malloc(size*sizeof(float));
    for(l=0;l<hsize;l++) {
      for(k=0;k<vsize;k++) {
        n=vsize*l+k;
        m=hsize*k+l;
        tmp[m]=buffer[n];
        }
      }
    for(l=0;l<size;l++) buffer[l]=tmp[l];
    free(tmp);
    }

  status=0;
  *data=buffer;
  return(status);

error:
  *data=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_read3d(const char *filename, int varid, int frame, float **data, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,close_status,ncid;
  
  status=nc_open(filename,NC_NOWRITE,&ncid);
  if(status!=NC_NOERR)return status;
  
  status=fe_read3d(ncid, varid, frame, data, mask);
  
  close_status=nc_close(ncid);
  
  if(status==NC_NOERR)
    status=close_status;

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_extract(int ncid, int varid, size_t * start, size_t * count, float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, m, status;
  int verbose = 0;
  unsigned int size;
  float scale, offset;

  cdfvar_t v_info;
  cdfgbl_t g_info;

  status = cdf_globalinfo(ncid, &g_info, verbose);

  status = cdf_varinfo(ncid, varid, &v_info);

/*------------------------------------------------------------------------------
  read field in netcdf file */
  status = nc_get_vara_float(ncid, varid, start, count, buffer);

  size = 1;
  for(k = 0; k < v_info.ndim; k++) {
    size *= count[k];
  }

  status = nc_get_att_float(ncid, varid, "scale_factor", &scale);
  if(status != NC_NOERR)
    scale = 1.0;

  status = nc_get_att_float(ncid, varid, "add_offset", &offset);
  if(status != NC_NOERR)
    offset = 0.0;

  if((scale != 1.0) || (offset != 0.0))
    for(m = 0; m < size; m++)
      buffer[m] = buffer[m] * scale + offset;

  return (status);

error:
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_getbuffer(int ncid,int varid, float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,status;
  int verbose=0;
  unsigned int size;
  float scale,offset;

  cdfgbl_t g_info;
  cdfvar_t v_info;

  status= cdf_globalinfo(ncid,&g_info,verbose);

  status= cdf_varinfo(ncid,varid,&v_info);

/*------------------------------------------------------------------------------
  read full variable field in netcdf file */
  status=nc_get_var_float(ncid,varid,buffer);

  size=1;
  for(k=0;k<v_info.ndim;k++) {
    size*=v_info.dim[k].length;
    }

  status=nc_get_att_float(ncid,varid,"scale_factor",&scale);
  if(status != NC_NOERR) scale=1.0;

  status=nc_get_att_float(ncid,varid,"add_offset",&offset);
  if(status != NC_NOERR) offset=0.0;

  if((scale !=1.0) || (offset !=0.0))
    for(m=0;m<size;m++) buffer[m]=buffer[m]*scale+offset;

  return(status);

error:
  return(status);

}

/*-----------------------------------------------------------------------------*/
/**
loads all the values of a variable
\date created 29 Apr 2011
\author Damien Allain
*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_getbuffer(int ncid,const char*varname, float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, id;

  status = nc_inq_varid(ncid, varname, &id);
  if(status != NC_NOERR) goto error;
  
  return poc_getbuffer(ncid, id, buffer);

error:
  nc_check_error(status, __LINE__, __FILE__);
  return (status);
}

