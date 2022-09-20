
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

\brief variable and grid loading poc-netcdf definitions

Old note : Grid loading routines, mostly obsolete
*/
/*------------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "tools-structures.h"
#include "tools-define.h"
#include "netcdf-proto.h"
#include "map.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_getgrid2d (const char* filename,cdfgbl_t global, cdfvar_t info, cdfvar_t x_info, cdfvar_t y_info, grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int i,j,k,status, ncid;
  size_t lengthhp;
  int var,vx=-1,vy=-1;
  nc_type xtypep;
  char *text=NULL;

  double *ddum=NULL,*vector=NULL,factor;
  decoded_t decoded;
  cdfvar_t vx_info,vy_info;
  cdfgbl_t grid_info,*grid_or_global;

  bool debug=false;
  
/*------------------------------------------------------------------------------
  check variable "associate" attribute */
  var=info.id;
  if(info.ndim <2) NC_TRAP_ERROR(return,NC_EBADDIM,1,"%s has %d dimensions which is <2",info.name,info.ndim);

  status= cdf_globalinfo (filename,&grid_info,0);
  if(status != NC_NOERR) NC_TRAP_ERROR(return,status,1,"cdf_globalinfo(\"%s\",,) error",filename);

  grid->circular=0; /*default*/

  for(i=0;i<2;i++){
    //try first grid file, then data file (where time should be)
    grid_or_global= (i==0? &grid_info:&global);
    status=poc_decode_associates(info,*grid_or_global, &decoded, i==0, debug);
    if(decoded.vx==-1 || decoded.vy==-1)
      continue;
    vx=decoded.vx;
    vy=decoded.vy;
    vx_info=grid_or_global->variable[vx];
    vy_info=grid_or_global->variable[vy];
    break;
    }

  if( (vx==-1) || (vy==-1) ) {
    status=-1;
    NC_TRAP_ERROR(return,status,1,"cannot identify grid variables : vx=%d vy=%d\n", vx, vy);
    }

  vx=vx_info.id;
  vy=vy_info.id;

/*------------------------------------------------------------------------------
  open file for reading */
  status=nc_open(filename,0,&ncid);
  if(status != NC_NOERR) NC_TRAP_ERROR(return,status,1,"nc_open(\"%s\",0,) error",filename);

// /*------------------------------------------------------------------------------
//   assume x dim is last and y dim is last-1 */
//   if(strcasestr(decoded.axis,"YX")!=0) {
//     grid->nx=info.dim[info.ndim-1].length;
//     grid->ny=info.dim[info.ndim-2].length;
//     }
//   else {
//     grid->nx=info.dim[info.ndim-2].length;
//     grid->ny=info.dim[info.ndim-1].length;
//     }
    
  grid->nx=decoded.xlen;
  grid->ny=decoded.ylen;
  grid->nz=1;

  ddum=new double [grid->nx*grid->ny];

/*------------------------------------------------------------------------------
  get x coordinate data */
  var=vx;
  switch (vx_info.ndim) {
    case 2:
      status=nc_get_var_double(ncid,var,ddum);
      if(status != NC_NOERR) {NC_CHKERR_LINE_FILE(status,"nc_get_var_double() error");goto error_after_open;}
      grid->modeH=2;
      break;

    case 1:
      exitIfNull(vector=new double[grid->nx]);
      status=nc_get_var_double(ncid,var,vector);
      if(status != NC_NOERR) {NC_CHKERR_LINE_FILE(status,"nc_get_var_double() error");goto error_after_open;}
      for(j=0;j<grid->ny;j++)
        for(i=0;i<grid->nx;i++) ddum[grid->nx*j+i]=vector[i];
      delete[] vector;
      grid->modeH=2;
      break;

    default:
      status=NC_EDIMSIZE;
      NC_CHKERR_LINE_FILE(status,"dimension of x=%d ???",vx_info.ndim);
      goto error_after_open;
      break;
    }

  grid->xmin=ddum[0];
  grid->xmax=ddum[grid->nx-1];

  factor=1.;
  text=vx_info.getattr("units");
  if(text!=NULL && strcmp(text,"radian_east")==0)factor=180./M_PI;
  //if(strcmp(vx_info.getattr("units"),"degree_east")==0)factor=1.;

  text=vx_info.getattr("topology");
  if(text!=NULL && strcmp(text,"circular")==0)grid->circular=1;

  grid->x=new double [grid->nx*grid->ny];
//  grid->x=(double *) malloc(grid->nx*grid->ny*sizeof(double));
  for(k=0;k<grid->nx*grid->ny;k++) grid->x[k]=ddum[k]*factor;
  for(k=0;k<grid->nx*grid->ny;k++) updatemin(&grid->xmin,grid->x[k]);
  for(k=0;k<grid->nx*grid->ny;k++) updatemax(&grid->xmax,grid->x[k]);
  grid->dx=(grid->xmax-grid->xmin)/(grid->nx-1);


/*------------------------------------------------------------------------------
  get y coordinate data */
  var=vy;
  switch (vy_info.ndim) {
    case 2:
      status=nc_get_var_double(ncid,var,ddum);
      if(status != NC_NOERR) {NC_CHKERR_LINE_FILE(status,"nc_get_var_double() error");goto error_after_open;}
      grid->modeH=2;
      break;

    case 1:
      exitIfNull(
        vector=new double[grid->ny]
        );
      status=nc_get_var_double(ncid,var,vector);
      if(status != NC_NOERR) {NC_CHKERR_LINE_FILE(status,"nc_get_var_double() error");goto error_after_open;}
      for(j=0;j<grid->ny;j++)
        for(i=0;i<grid->nx;i++) ddum[grid->nx*j+i]=vector[j];
      delete[] vector;
      grid->modeH=2;
      break;

    default:
      status=NC_EDIMSIZE;
      NC_CHKERR_LINE_FILE(status,"dimension of y=%d ???\n",vy_info.ndim);
      goto error_after_open;
      break;
    }

  grid->ymin=ddum[0];
  grid->ymax=ddum[grid->nx*grid->ny-1];

  factor=1.;
  #warning REDO AS ABOVE
  status=nc_inq_att(ncid,var,"units",&xtypep,&lengthhp);
  if(status == NC_NOERR) {
    exitIfNull(
      text=(char *) malloc(lengthhp+1)
      );
    status=nc_get_att_text(ncid,var,"units",text);
    text[lengthhp]='\0';
    if(strcmp(text,"radian_north")==0) factor=180./M_PI;
    if(strcmp(text,"degree_north")==0) factor=1.;
    free(text);
    }

  grid->y=new double [grid->nx*grid->ny];
//  grid->y=(double *) malloc(grid->nx*grid->ny*sizeof(double));
  for(k=0;k<grid->nx*grid->ny;k++) grid->y[k]=ddum[k]*factor;
  for(k=0;k<grid->nx*grid->ny;k++) updatemin(&grid->ymin,grid->y[k]);
  for(k=0;k<grid->nx*grid->ny;k++) updatemax(&grid->ymax,grid->y[k]);
  grid->dy=(grid->ymax-grid->ymin)/(grid->ny-1);

  delete[] ddum;

  grid->zmin=0.;
  grid->zmax=0.;
  grid->z=NULL;
  grid->nz=1;
  grid->dz=0.0;

  grid->modeH= 2;
  grid->modeV=-1;


/*------------------------------------------------------------------------------
  default */
  grid->circular=0;
  grid->overlapped=0;
  grid->connex=1;
 /*
  mapc_printgrid3d(*grid);
 */

  grid->connex=mapc_checkconnexity(*grid);
  if(grid->x[grid->nx-1]-grid->x[0]>360.0) {
     grid->circular=1;
     grid->overlapped=1;
     }

  if(grid->xmax-grid->xmin>359.0) {
     grid->circular=1;
     grid->overlapped=1;
     }

error_after_open:
  status=nc_close(ncid);
  if(status != NC_NOERR) NC_CHKERR_LINE_FILE(status,"nc_close() error on %s",filename);

/* *----------------------------------------------------------------------------
  patch: commented to avoid pb when time and space grid not in the same file */
//  status= poc_gettime(filename,  info,  global, &(grid->origine), &(grid->time), &(grid->nt));

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_getgrid2d (const char* filename,cdfgbl_t global, cdfvar_t v_info, grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  cdfvar_t x_info, y_info;
  
  status=poc_getgrid2d (filename, global,  v_info, x_info, y_info, grid);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_getgrid2d(const char *input, const char *variable, grid_t *grid,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    v,status;
  cdfgbl_t info;
  cdfvar_t h_info;

  status=cdf_globalinfo(input,&info,verbose);
  if(verbose){
    for (v=0;v<info.nvarsp;v++) {
      printf("variable %3d: name %s, type %d,ndim %d \n",v,(info.variable[v]).name,info.variable[v].type,info.variable[v].ndim);
      }
    }
  
  status=cdf_varinfo(input,variable,&h_info,verbose);
  status=poc_getgrid2d(input, info, h_info, grid);

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_getgrid2d(const char *input, const char *variable, const char *xname, const char *yname, grid_t *grid, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    v,status;
  cdfgbl_t info;
  cdfvar_t h_info;

  status=cdf_globalinfo(input,&info,verbose);
  if(verbose){
    for (v=0;v<info.nvarsp;v++) {
      printf("variable %3d: name %s, type %d,ndim %d \n",v,(info.variable[v]).name,info.variable[v].type,info.variable[v].ndim);
      }
    }
  
  status=cdf_varinfo(input,variable,&h_info,verbose);
  status=poc_getgrid2d(input, info, h_info, grid);

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int poc_getvar2d_template(const char* filename, int varid, int frame, T *buf, T *mask, cdfvar_t info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  In the following, we assume a standard field structure with typical TZYX
  arrangement.
------------------------------------------------------------------------*/
  {
  int status, ncid;
  size_t *start,*count;
  T scale=1.0,offset=0.0;
  int verbose=0;
  cdfgbl_t global;
  decoded_t decoded;

  if (verbose) status=cdf_info(filename);

  status=cdf_globalinfo(filename,&global,verbose);

  status=nc_open(filename,0,&ncid);
  if (status !=NC_NOERR) {return(status);}
/*------------------------------------------------------------------------------
  inquire start and count*/
  start=new size_t[info.ndim];
  count=new size_t[info.ndim];

  status= poc_decode_axis(info,global,&decoded);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"poc_decode_axis() error");

  start[decoded.xdim]=0;
  count[decoded.xdim]=decoded.xlen;

  start[decoded.ydim]=0;
  count[decoded.ydim]=decoded.ylen;

  if(decoded.zdim!=-1){
    start[decoded.zdim]=0;
    count[decoded.zdim]=1;
    }

  if(decoded.tdim!=-1){
    if(frame==-1) frame=decoded.tlen-1;
    start[decoded.tdim]=frame;
    count[decoded.tdim]=1;
    }

  status=poc_get_vara(ncid,varid,start,count,buf);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"poc_get_vara() error");

  status= poc_decode_mask(info, &decoded);
  if(status != NC_NOERR) {
    printf(";mask value not documented, 0 used by default.");
//     scale=1;
//     offset=0;
    status=0;
    }
//   else {
    scale=decoded.scale;
    offset=decoded.offset;
//     }
/*------------------------------------------------------------------------------
  pass mask */
  *mask=decoded.spec;
  
/*------------------------------------------------------------------------------
  aply scale and offset */
  if(scale!=1.0 or offset!=0) {
    for(size_t m=0;m<decoded.xlen*decoded.ylen;m++) {
      bool masked;
      if(isnan(*mask)) masked=isnan(buf[m]);
      else masked=(buf[m]==*mask);
      if(masked) {
        buf[m]=buf[m]*scale+offset;
//         if(isnan(buf[m])) {
//           printf("poc_getvar2d_template:Nan value out of mask in %s\n", filename);
//           buf[m]=*mask;
//           }
        }
      }
    }

//   if((isinf(*mask)!=0) || (isinf(*mask)!=0)) {
//     for(m=0;m<decoded.xlen*decoded.ylen;m++) {
//       if(buf[m]==*mask) {
// 	buf[m]=1.e10;
//         }
//       }
//     *mask=1.e10;
//     }
  bool swap_needed;
  switch (info.ndim){
    case 2:
      swap_needed=strcmp(decoded.axis,"XY")==0;
      break;
    case 3:
      swap_needed=strcmp(decoded.axis,"TXY")==0;
      break;
    }
  if(swap_needed) {
/*------------------------------------------------------------------------------
    X,Y,Z dimension are in the wrong order */
    float *tmp=new float[decoded.xlen*decoded.ylen];
    int i,j;
    size_t k1,k2;
    for(j=0;j<decoded.ylen;j++)
      for(i=0;i<decoded.xlen;i++) {
        k1=j*decoded.xlen+i;
        k2=i*decoded.ylen+j;
        tmp[k1]=buf[k2];
        }
    for(j=0;j<decoded.ylen;j++)
      for(i=0;i<decoded.xlen;i++) {
        k1=j*decoded.xlen+i;
        buf[k1]=tmp[k1];
        }
    delete[] tmp;
    }

  if(status != NC_NOERR) goto error;
  status=nc_close(ncid);
  status=0;
  return(status);

error:
  nc_check_error(status,__LINE__,__FILE__);
  status=nc_close(ncid);
  status=-1;
  printf("loading variable %d failed\n",varid);
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_getvar2d  (const char* filename, int varid, int frame, signed char *buf, signed char *mask, cdfvar_t info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_getvar2d_template(filename, varid, frame, buf, mask, info);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_getvar2d  (const char* filename, int varid, int frame, short *buf, short *mask, cdfvar_t info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_getvar2d_template(filename, varid, frame, buf, mask, info);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_getvar2d  (const char* filename, int varid, int frame, float *buf, float *mask, cdfvar_t info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_getvar2d_template(filename, varid, frame, buf, mask, info);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_getvar2d  (const char* filename, int varid, int frame, double *buf, double *mask, cdfvar_t info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=poc_getvar2d_template(filename, varid, frame, buf, mask, info);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_getvar3d (const char* filename, int varid, int frame, float *buf, float *mask, cdfvar_t info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  In the following, we assume a standard field structure with typical TZYX
  arrangement.
------------------------------------------------------------------------*/
  {
  int m,status, ncid;
  size_t *start,*count;
  float missing=1e+10;
  int verbose=0;
  cdfgbl_t global;
  decoded_t decoded;
  
  if(varid==-1) goto abort;

  if (verbose) status=cdf_info(filename);

  status=cdf_globalinfo(filename,&global,verbose);

  status=nc_open(filename,0,&ncid);
  if (status !=NC_NOERR) {
    status=-1;
    return(status);
    }
/*------------------------------------------------------------------------------
  inquire start and count*/
  start=new size_t[info.ndim];
  count=new size_t[info.ndim];

  status= poc_decode_axis(info,global,&decoded);

  start[decoded.xdim]=0;
  count[decoded.xdim]=decoded.xlen;

  start[decoded.ydim]=0;
  count[decoded.ydim]=decoded.ylen;

  if(decoded.zdim!=-1){
    start[decoded.zdim]=0;
    count[decoded.zdim]=decoded.zlen;
    }

  if(decoded.tdim!=-1){
    start[decoded.tdim]=frame;
    count[decoded.tdim]=1;
    }

  status=nc_get_vara_float(ncid,varid,start,count,buf);
  if(status != NC_NOERR) goto error;

  status= poc_decode_mask(info, &decoded);
/*------------------------------------------------------------------------------
  pass mask */
  if(status == NC_NOERR) {
    *mask=decoded.spec;
    }
  else {
    status = NC_NOERR;
    }
/*------------------------------------------------------------------------------
  aply scale and offset */
  for(m=0;m<decoded.xlen*decoded.ylen*decoded.zlen;m++) {
    if(isnan(buf[m])) {
//      printf("poc_getvar3d : Nan value out of mask in %s\n", filename);
      *mask=missing;
      buf[m]=missing;
      }
    if(buf[m]!=*mask) buf[m]=buf[m]*decoded.scale+decoded.offset;
    }

  if(status != NC_NOERR) goto error;
  status=nc_close(ncid);
  status=0;
  return(status);

error:
  nc_check_error(status,__LINE__,__FILE__);
  status=nc_close(ncid);
abort:
  status=-1;
  printf("loading variable %d failed\n",varid);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int poc_getgrid3d (const char* filename, cdfgbl_t g_info, cdfvar_t v_info, grid_t *grid, int *vdepth)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int i,j,k,/*k1,k2,*/m,n,status, ncid;
  size_t /*lengthhp,*/*dimlgth=NULL;
  int var,/*vx,vy,*/vz,vmask=-1,vsigma=-1,axis_id;
  int *ndim=NULL,**dimids=NULL;
//   nc_type xtypep;
  decoded_t decoded,vx_decoded,vy_decoded,vz_decoded;
  cdfvar_t vx_info,vy_info,vz_info;
  cdfgbl_t grid_info;

  double *ddum=NULL,*vector=NULL,/*tmp,*/factor,*zdum=NULL,z0;
//   int swap_needed=0;

  int verbose=0;
  bool debug=false;

/*------------------------------------------------------------------------------
  open file for reading */
  status=nc_open(filename,0,&ncid);
  NC_TRAP_ERROR(return,status,1,"nc_open(\"%s\",0,) error",filename);

  status= cdf_globalinfo(ncid, &grid_info,verbose);

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

  status=poc_decode_axis(v_info, g_info, &decoded);
/* *------------------------------------------------------------------------
  try first grid file, then data file (where time should be)*/
  status=poc_decode_associates(v_info, grid_info, &decoded,1, debug);
  status=poc_decode_associates(v_info, g_info, &decoded,0, debug);

  if(status!=0) goto error;

  vmask=-1;

  printf("%s grid parsing ----> variable x: %d, variable y: %d variable z: %d variable sigma: %d\n",__func__,decoded.vx,decoded.vy,decoded.vz,vsigma);

//  if( (vx==-1) || (vy==-1) || ((vz==-1) && (vsigma==-1)) ) goto error;
  if( (decoded.vx==-1) || (decoded.vy==-1)) goto error;

/*------------------------------------------------------------------------------
  get dimension */
  grid->nx=decoded.xlen;
  grid->ny=decoded.ylen;
  grid->nz=decoded.zlen;

  exitIfNull(
    ddum=(double *) malloc(grid->nx*grid->ny*sizeof(double))
    );

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  get x coordinate data 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  var=decoded.vx;
  status= cdf_varinfo(ncid,var,&vx_info);
  var=decoded.vy;
  status= cdf_varinfo(ncid,var,&vy_info);
  
  status=poc_decode_axis(vx_info, grid_info, &vx_decoded);
  status=poc_decode_axis(vy_info, grid_info, &vy_decoded);

  var=decoded.vx;
  switch (vx_info.ndim) {
    case 1:
/*------------------------------------------------------------------------------
      read vector and reconstruct x 2D-array*/
      exitIfNull(
        vector=(double *) malloc(grid->nx*sizeof(double))
        );
      status=nc_get_var_double(ncid,var,vector);
      if(status != NC_NOERR) goto error;
      for(j=0;j<grid->ny;j++)
        for(i=0;i<grid->nx;i++) ddum[grid->nx*j+i]=vector[i];
      free(vector);
      grid->modeH=2;
      break;

    case 2:
/*------------------------------------------------------------------------------
      read x 2D-array*/
      if(vx_decoded.xlen!=grid->nx) {
        printf("inconsistent grid, nx dimension\n");
        goto error;
        }
      if(vx_decoded.ylen!=grid->ny) {
        printf("inconsistent grid, ny dimension\n");
        goto error;
        }
      status=nc_get_var_double(ncid,var,ddum);
      if(status != NC_NOERR) goto error;
      break;

    default:
/*------------------------------------------------------------------------------
    structural anomaly */
      goto error;
      break;
    }

  grid->xmin=ddum[0];
  grid->xmax=ddum[grid->nx-1];

  factor=1.;
  status= poc_decode_units(vx_info,&vx_decoded,&factor);

  exitIfNull(
    grid->x=new double[grid->nx*grid->ny]
    );
  for(k=0;k<grid->nx*grid->ny;k++) grid->x[k]=ddum[k]*factor;
  for(k=0;k<grid->nx*grid->ny;k++) updatemin(&grid->xmin,grid->x[k]);
  for(k=0;k<grid->nx*grid->ny;k++) updatemax(&grid->xmax,grid->x[k]);
  grid->dx=(grid->xmax-grid->xmin)/(grid->nx-1);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  get y coordinate data 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  var=decoded.vy;
  switch (vy_info.ndim) {
    case 1:
/*------------------------------------------------------------------------------
      read vector and reconstruct y 2D-array*/
      exitIfNull(
        vector=(double *) malloc(grid->ny*sizeof(double))
        );
      status=nc_get_var_double(ncid,var,vector);
      if(status != NC_NOERR) goto error;
      for(j=0;j<grid->ny;j++)
        for(i=0;i<grid->nx;i++) ddum[grid->nx*j+i]=vector[j];
      free(vector);
      grid->modeH=2;
      break;

    case 2:
/*------------------------------------------------------------------------------
      read y 2D-array*/
      if(vy_decoded.xlen!=grid->nx) {
        printf("inconsistent grid, nx dimension\n");
        goto error;
        }
      if(vy_decoded.ylen!=grid->ny) {
        printf("inconsistent grid, ny dimension\n");
        goto error;
        }
    status=nc_get_var_double(ncid,var,ddum);
    if(status != NC_NOERR) goto error;
    grid->modeH=2;
    break;

    default:
/*------------------------------------------------------------------------------
      structural anomaly */
      goto error;
      break;
    }

  grid->ymin=ddum[0];
  grid->ymax=ddum[grid->nx*grid->ny-1];

  factor=1.;
  status= poc_decode_units(vy_info,&vy_decoded,&factor);

//   status=nc_inq_att(ncid,var,"topology",&xtypep,&lengthhp);
//   if(status == NC_NOERR) {
//     text=(char *) malloc(lengthhp+1);
//     status=nc_get_att_text(ncid,var,"topology",text);
//     text[lengthhp]='\0';
//     if(strcmp(text,"circular")==0) grid->circular=1;
//     free(text);
//     }

  exitIfNull(
    grid->y=new double[grid->nx*grid->ny]
    );
  for(k=0;k<grid->nx*grid->ny;k++) grid->y[k]=ddum[k]*factor;
  for(k=0;k<grid->nx*grid->ny;k++) updatemin(&grid->ymin,grid->y[k]);
  for(k=0;k<grid->nx*grid->ny;k++) updatemax(&grid->ymax,grid->y[k]);
  grid->dy=(grid->ymax-grid->ymin)/(grid->ny-1);

  free(ddum);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  get z coordinate data in case z_dim ok 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(decoded.vz!=-1) {
    var=decoded.vz;
    status= cdf_varinfo(ncid,var,&vz_info);

    exitIfNull(
      grid->z=new double[grid->nx*grid->ny*grid->nz]
      );

    status= poc_decode_mask(vz_info,&vz_decoded);
    if(status ==0) grid->zmask=vz_decoded.spec;
    else grid->zmask=0.0;

//    printf("cdf_loadvar_gd_3d : warning, zmask modification\n");
    status= cdf_varinfo(ncid,var,&vz_info);

    switch (vz_info.ndim) {
      case 1:
/*------------------------------------------------------------------------------
      uniform z levels */
        exitIfNull(
          ddum=new double[grid->nz]
          );
        status=nc_get_var_double(ncid,var,ddum);
        if(status != NC_NOERR) goto error;
        for(k=0;k<grid->nz;k++)
          for(m=0;m<grid->nx*grid->ny;m++)
            grid->z[k*grid->nx*grid->ny+m]=ddum[k];
        free(ddum);
        grid->modeV=1;
        break;

      case 2:
/*------------------------------------------------------------------------------
      structural anomaly */
/*       goto error; */
      break;

      case 3:
/*------------------------------------------------------------------------------
        x,y dependant z levels */
        status=nc_get_var_double(ncid,var,grid->z);
        if(status != NC_NOERR) goto error;
        grid->modeV=3;
        break;

      default:
/*------------------------------------------------------------------------------
      structural anomaly */
        goto error;
        break;
      }
    grid->zmin=1.e10;
    grid->zmax=1.e-10;

    for(k=0;k<grid->nx*grid->ny*grid->nz;k++) updatemin(&grid->zmin,grid->z[k]);
    for(k=0;k<grid->nx*grid->ny*grid->nz;k++) updatemax(&grid->zmax,grid->z[k]);
    grid->dz=(grid->zmax-grid->zmin)/(grid->nz-1);
    }

  if(decoded.vz!=-1)
  if(decoded.vz==decoded.vt) {
    z0=grid->z[0];
    for(k=0;k<grid->nz;k++) {
      for(m=0;m<grid->nx*grid->ny;m++) grid->z[k*grid->nx*grid->ny+m]-=z0;
      }
    for(k=0;k<grid->nz;k++) {
      for(m=0;m<grid->nx*grid->ny;m++) grid->z[k*grid->nx*grid->ny+m]/=3600.0;
      }
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  get s coordinate data in case s_dim ok 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(vsigma!=-1) {
    *vdepth=vz; /*detected bathymetry*/
    var=vsigma;
    grid->z=NULL;
    switch (ndim[var]) {
      case 1:
/*------------------------------------------------------------------------------
        assume uniform sigma levels */
        exitIfNull(
          grid->sigma=(double *) malloc(grid->nz*sizeof(double))
          );
        status=nc_get_var_double(ncid,vsigma,grid->sigma);
        if(status != NC_NOERR) goto error;
        grid->modeV=2;
/*------------------------------------------------------------------------------
        load bathymetry*/
        var=vz;
        if(dimlgth[dimids[var][0]]*dimlgth[dimids[var][1]]==grid->nx*grid->ny) {
          exitIfNull(
            grid->z=(double *) malloc(grid->nx*grid->ny*grid->nz*sizeof(double))
            );
          exitIfNull(
            zdum=(double *) malloc(dimlgth[dimids[var][0]]*dimlgth[dimids[var][1]]*sizeof(double))
            );
          status=nc_get_var_double(ncid,var,zdum);
/*------------------------------------------------------------------------------
          compute sigma level depth*/
          for(k=0;k<grid->nz;k++)
            for(m=0;m<grid->nx*grid->ny;m++) grid->z[k*grid->nx*grid->ny+m]=zdum[m]*grid->sigma[k];
          free(zdum);
          grid->zmin=1.e10;
          grid->zmax=1.e-10;

          for(k=0;k<grid->nx*grid->ny*grid->nz;k++) updatemin(&grid->zmin,grid->z[k]);
          for(k=0;k<grid->nx*grid->ny*grid->nz;k++) updatemax(&grid->zmax,grid->z[k]);
          grid->dz=(grid->zmax-grid->zmin)/(grid->nz-1);
          }
        break;

      case 2:
/*------------------------------------------------------------------------------
        structural anomaly */
        goto error;
        break;

      case 3:
/*------------------------------------------------------------------------------
        x,y dependant s levels */
        exitIfNull(
          grid->sigma=(double *) malloc(grid->nx*grid->ny*grid->nz*sizeof(double))
          );
        status=nc_get_var_double(ncid,var,grid->sigma);
        if(status != NC_NOERR) goto error;
        grid->modeV=3;
        break;

      default:
/*------------------------------------------------------------------------------
        structural anomaly */
        goto error;
        break;
      }
    }
//
//   if(vmask!=-1) {
//     zdum=(double *) malloc(dimlgth[dimids[vmask][0]]*dimlgth[dimids[vmask][1]]*sizeof(double));
//     grid->mask=(signed char *) malloc(grid->nx*grid->ny);
//     status=nc_get_var_double(ncid,vmask,zdum);
//     if(status != NC_NOERR) goto error;
//     for(m=0;m<grid->nx*grid->ny;m++) grid->mask[m]=floor(zdum[m]+0.5);
//     free(zdum);
//     }
//
//   if(swap_needed) {
//     for(j=0;j<grid->ny;j++)
//       for(i=0;i<grid->nx;i++) {
// 	k1=j*grid->nx+i;
// 	k2=i*grid->ny+j;
//         tmp=grid->x[k1];
//         grid->x[k1]=grid->x[k2];
//         grid->x[k2]=tmp;
//         tmp=grid->y[k1];
//         grid->y[k1]=grid->y[k2];
//         grid->y[k2]=tmp;
// 	}
//     }

  grid->modeH=2;

/*------------------------------------------------------------------------------
  default */
  grid->circular=0;
  grid->overlapped=0;
  grid->connex=1;

  grid->connex=mapc_checkconnexity(*grid);
  if(grid->x[grid->nx-1]-grid->x[0]>360.0){
     grid->circular=1;
     grid->overlapped=1;
     }

  if(grid->xmax-grid->xmin>359.0) {
     grid->circular=1;
     grid->overlapped=1;
     }

  status=nc_close(ncid);
  NC_TRAP_ERROR(return,status,1,"nc_close() error with %s",filename);

  return(status);

error:
  nc_check_error(status,__LINE__,__FILE__);
  if(dimlgth != NULL) free(dimlgth);
  status=nc_close(ncid);
  nc_check_error(status,__LINE__,__FILE__);
  status=-1;
  printf("loading 3D grid for variable %d failed\n",v_info.id);
  return(status);

}

