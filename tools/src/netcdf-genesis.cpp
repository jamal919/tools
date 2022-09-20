
/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief grid and variable loading routines (for compatibility with GENESIS, soon obsoletes)
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#ifdef IEEE
#include <ieeefp.h>
#endif

#include "tools-structures.h"
#include "tools-define.h"
#include "poc-netcdf.def"
#include "netcdf-proto.h"

#include "poc-time.h"
#include "map.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int cdf_loadvar_r1_3d (const char* filename, int varid, int frame, grid_t grid, float *buf, float *mask, variable_t *info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int option,status, ncid;
  int axis_id,n;
  int v,vmask=-1,vtime=-1;
  int i,j,k,k1,k2,m,fmt,swap_needed=0,d,missing_found=0;
  size_t *start,*count;
  float fill,missing=1.e+10,scale=1.0,offset=0.0,*tmp;
  int tdim,nframes;
  double *time,factor;
  int verbose=0,chk=0;
  cdfgbl_t g_info;
  cdfvar_t v_info;
  decoded_t decoded;
  date_t origine;

  info->name=NULL;
  info->standard_name=NULL;
  info->long_name=NULL;

  status=nc_open(filename,0,&ncid);
  if (status !=NC_NOERR) {
    status=-1;
    return(status);
    }

  status= cdf_globalinfo(ncid,&g_info,verbose);
  if(status!=0) goto error;

  status= cdf_varinfo(ncid,varid,&v_info);
  if(status!=0) goto error;

/*------------------------------------------------------------------------
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

  if(status!=0) goto error;

  vmask=-1;

  start=new size_t[v_info.ndim];
  count=new size_t[v_info.ndim];

  start[decoded.xdim]=0;
  count[decoded.xdim]=decoded.xlen;

  start[decoded.ydim]=0;
  count[decoded.ydim]=decoded.ylen;

  start[decoded.zdim]=0;
  count[decoded.zdim]=decoded.zlen;

  if(decoded.tdim!=-1){
/* *----------------------------------------------------------------------------
    warning: frame index start at 1, netcdf at 0*/
    start[decoded.tdim]=frame-1;
    count[decoded.tdim]=1;
    }

  status=nc_get_vara_float(ncid,varid,start,count,buf);
  if(status != NC_NOERR) goto error;

/*------------------------------------------------------------------------
  patch: force missing value from mask if exists */
  status= poc_decode_mask(v_info, &decoded);
  *mask=decoded.spec;
  scale=decoded.scale;
  offset=decoded.offset;

/*------------------------------------------------------------------------
  */
  if(isnan(*mask)) {
    for(k=0;k<grid.nz;k++) {
      for(j=0;j<grid.ny;j++) {
        for(i=0;i<grid.nx;i++) {
          m=k*grid.ny*grid.nx+j*grid.nx+i;
          if(isnan(buf[m])) {
            buf[m]=1.e+10;
            }
          }
        }
      }
    *mask=1.e+10;
    }
    
/*------------------------------------------------------------------------
  aply scale and offset */
  if((scale!=1.) && (offset!=0))
    for(k=0;k<grid.nz;k++)
      for(j=0;j<grid.ny;j++)
        for(i=0;i<grid.nx;i++) {
          m=k*grid.ny*grid.nx+j*grid.nx+i;
          if(isnan(buf[m])) {
            printf("Nan value out of mask\n");
            buf[m]=*mask;
            }
          if(buf[m]!=*mask) buf[m]=buf[m]*scale+offset;
          }

/*------------------------------------------------------------------------
  figure out time of frame */
  status= poc_gettime(ncid, v_info, g_info, &origine, &time, &nframes);
  if(status==0) {
/*------------------------------------------------------------------------
    time origin */
//    info->origin=origine;
    info->time=time[frame-1];
    delete[] time;
    }
  else
    info->time=0;

  status=nc_close(ncid);
  nc_check_error(status,__LINE__,__FILE__);

  info->name   =strdup(v_info.name);
  info->origin =strdup(filename);
  info->ndim   =v_info.ndim;
  info->id     =varid;

  status= poc_decode_names(v_info, &decoded);
  status= poc_decode_names(g_info, &decoded);

  info->standard_name=NULL;
  info->long_name=NULL;

  if(decoded.standard_name==0) decoded.standard_name=strdup("no_standard_name");
  if(decoded.long_name==0)    decoded.long_name=strdup("no_long_name");

  if(strlen(decoded.standard_name) != 0)
    info->standard_name=strdup(decoded.standard_name);
  if(strlen(decoded.long_name) != 0) {
    info->long_name=strdup(decoded.long_name);
    }
  else {
    info->long_name=strdup("no name");
    }

  if(info->standard_name==NULL) info->standard_name=strdup(info->long_name);
  for(k=0;k<strlen(info->standard_name);k++)
    if(info->standard_name[k]==' ') info->standard_name[k]='_';

  info->units=NULL;
  if(strlen(decoded.units) != 0)
    info->units=strdup(decoded.units);
  if(info->units==NULL) info->units=strdup("no units");

  swap_needed=strcmp(decoded.axis,"YX")==0;

  if(swap_needed) {
/*------------------------------------------------------------------------
    X,Y,Z dimension are in the wrong order */
    tmp=new float[grid.nx*grid.ny];
    for(j=0;j<grid.ny;j++)
      for(i=0;i<grid.nx;i++) {
        k1=j*grid.nx+i;
        k2=i*grid.ny+j;
        tmp[k1]=buf[k2];
        }
    for(j=0;j<grid.ny;j++)
      for(i=0;i<grid.nx;i++) {
        k1=j*grid.nx+i;
        buf[k1]=tmp[k1];
        }
    delete[] tmp;
    }

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

extern "C" int cdf_loadvar_gd_3d (char* filename,cdfgbl_t g_info,cdfvar_t v_info, grid_t *grid,int *vdepth)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  FILE *in;
  int i,j,k,k1,k2,m,n,option,status, ncid;
  size_t lengthhp,*dimlgth=NULL,index[1];
  float time;
  int d, dim,v,var,att,vx,vy,vz,vt,vmask=-1,vsigma=-1,axis_id;
  int *ndim=NULL,**dimids=NULL;
  nc_type *vartype,xtypep;
  char *associate=NULL,*maskname=NULL,*axis=NULL;
  decoded_t decoded,vx_decoded,vy_decoded,vz_decoded;
  size_t *start,*count;
  cdfvar_t vx_info,vy_info,vz_info;
  cdfgbl_t grid_info;

  char bdum;
  int idum;
  float fdum;
  double *ddum=NULL,*vector=NULL,tmp,factor,*zdum,*sdum,z0;
  int swap_needed=0;
  int x_dim=-1,y_dim=-1,z_dim=-1,s_dim=-1,t_dim=-1,p,check;

  int verbose=0;
  bool debug=false;


/*------------------------------------------------------------------------
  open file for reading */
  status=nc_open(filename,0,&ncid);
  if(status != NC_NOERR) goto open_error;

  status= cdf_globalinfo(ncid, &grid_info,verbose);

/*------------------------------------------------------------------------
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
  status=poc_decode_associates(v_info, grid_info, &decoded,1,debug);
  status=poc_decode_associates(v_info, g_info, &decoded,0,debug);

  if(status!=0) goto error;

  vmask=-1;

  printf("---->variable x: %d, variable y: %d variable z: %d variable sigma: %d\n",decoded.vx,decoded.vy,decoded.vz,vsigma);

//  if( (vx==-1) || (vy==-1) || ((vz==-1) && (vsigma==-1)) ) goto error;
  if( (decoded.vx==-1) || (decoded.vy==-1)) goto error;

/*------------------------------------------------------------------------
  get dimension */
  grid->nx=decoded.xlen;
  grid->ny=decoded.ylen;
  grid->nz=decoded.zlen;

  ddum=(double *) malloc(grid->nx*grid->ny*sizeof(double));

/*------------------------------------------------------------------------
  get x coordinate data */
  var=decoded.vx;
  status= cdf_varinfo(ncid,var,&vx_info);
  status=poc_decode_mask(vx_info, &vx_decoded);
  switch (vx_info.ndim) {
    case 1:
/*------------------------------------------------------------------------
      read vector and reconstruct x 2D-array*/
      vector=(double *) malloc(grid->nx*sizeof(double));
      status=nc_get_var_double(ncid,var,vector);
      if(status != NC_NOERR) goto error;
      for(j=0;j<grid->ny;j++)
        for(i=0;i<grid->nx;i++) ddum[grid->nx*j+i]=vector[i];
      free(vector);
      grid->modeH=2;
      break;

    case 2:
/*------------------------------------------------------------------------
      read x 2D-array*/
      status=nc_get_var_double(ncid,var,ddum);
      if(status != NC_NOERR) goto error;
      break;

    default:
/*------------------------------------------------------------------------
    structural anomaly */
      goto error;
      break;
    }

  if(vx_decoded.offset!=0) {
    for(n=0;n<decoded.ylen*decoded.xlen;n++) ddum[n]+=vx_decoded.offset;
    }
  if(vx_decoded.scale!=1) {
    for(n=0;n<decoded.ylen*decoded.xlen;n++) ddum[n]*=vx_decoded.scale;
    }

  grid->xmin=ddum[0];
  grid->xmax=ddum[grid->nx-1];

  factor=1.;
  status= poc_decode_units(vx_info,&vx_decoded,&factor);

  grid->x=(double *) malloc(grid->nx*grid->ny*sizeof(double));
  for(k=0;k<grid->nx*grid->ny;k++) grid->x[k]=ddum[k]*factor;
  for(k=0;k<grid->nx*grid->ny;k++) grid->xmin=MIN(grid->xmin,grid->x[k]);
  for(k=0;k<grid->nx*grid->ny;k++) grid->xmax=MAX(grid->xmax,grid->x[k]);
  grid->dx=(grid->xmax-grid->xmin)/(grid->nx-1);

/*------------------------------------------------------------------------
  get y coordinate data */
  var=decoded.vy;
  status= cdf_varinfo(ncid,var,&vy_info);
  status=poc_decode_mask(vy_info, &vy_decoded);
  switch (vy_info.ndim) {
    case 1:
/*------------------------------------------------------------------------
    read vector and reconstruct y 2D-array*/
    vector=(double *) malloc(grid->ny*sizeof(double));
    status=nc_get_var_double(ncid,var,vector);
    if(status != NC_NOERR) goto error;
    for(j=0;j<grid->ny;j++)
      for(i=0;i<grid->nx;i++) ddum[grid->nx*j+i]=vector[j];
    free(vector);
    grid->modeH=2;
    break;

    case 2:
/*------------------------------------------------------------------------
    read y 2D-array*/
    status=nc_get_var_double(ncid,var,ddum);
    if(status != NC_NOERR) goto error;
    grid->modeH=2;
    break;

    default:
/*------------------------------------------------------------------------
    structural anomaly */
    goto error;
    break;
    }

  if(vy_decoded.offset!=0) {
    for(n=0;n<decoded.ylen*decoded.xlen;n++) ddum[n]+=vy_decoded.offset;
    }
  if(vy_decoded.scale!=1) {
    for(n=0;n<decoded.ylen*decoded.xlen;n++) ddum[n]*=vy_decoded.scale;
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


  grid->y=(double *) malloc(grid->nx*grid->ny*sizeof(double));
  for(k=0;k<grid->nx*grid->ny;k++) grid->y[k]=ddum[k]*factor;
  for(k=0;k<grid->nx*grid->ny;k++) grid->ymin=MIN(grid->ymin,grid->y[k]);
  for(k=0;k<grid->nx*grid->ny;k++) grid->ymax=MAX(grid->ymax,grid->y[k]);
  grid->dy=(grid->ymax-grid->ymin)/(grid->ny-1);

  free(ddum);

/*------------------------------------------------------------------------
  get z coordinate data in case z_dim ok*/
  if(decoded.vz!=-1) {
    var=decoded.vz;
    status= cdf_varinfo(ncid,var,&vz_info);

    grid->z=(double *) malloc(grid->nx*grid->ny*grid->nz*sizeof(double));

    status= poc_decode_mask(vz_info,&vz_decoded);
    if(status ==0) grid->zmask=vz_decoded.spec;
    else grid->zmask=0.0;

    printf("cdf_loadvar_gd_3d : warning, zmask modification");
    status= cdf_varinfo(ncid,var,&vz_info);

    switch (vz_info.ndim) {
      case 1:
/*------------------------------------------------------------------------
      uniform z levels */
        ddum=(double *) malloc(grid->nz*sizeof(double));
        status=nc_get_var_double(ncid,var,ddum);
        if(status != NC_NOERR) goto error;
        if(vz_decoded.offset!=0) {
          for(n=0;n<decoded.zlen;n++) ddum[n]+=vz_decoded.offset;
          }
        if(vz_decoded.scale!=1) {
          for(n=0;n<decoded.zlen;n++) ddum[n]*=vz_decoded.scale;
          }
        for(k=0;k<grid->nz;k++)
          for(m=0;m<grid->nx*grid->ny;m++) grid->z[k*grid->nx*grid->ny+m]=ddum[k];
        free(ddum);
        grid->modeV=1;
        break;

      case 2:
/*------------------------------------------------------------------------
      structural anomaly */
/*       goto error; */
      break;

      case 3:
/*------------------------------------------------------------------------
        x,y dependant z levels */
        status=nc_get_var_double(ncid,var,grid->z);
        if(status != NC_NOERR) goto error;
        if(vz_decoded.offset!=0) {
          for(n=0;n<decoded.zlen*decoded.ylen*decoded.xlen;n++) ddum[n]+=vz_decoded.offset;
          }
        if(vz_decoded.scale!=1) {
          for(n=0;n<decoded.zlen*decoded.ylen*decoded.xlen;n++) ddum[n]*=vz_decoded.scale;
          }
        grid->modeV=3;
        break;

      default:
/*------------------------------------------------------------------------
      structural anomaly */
        goto error;
        break;
      }
    grid->zmin=1.e+10;
    grid->zmax=1.e-10;

    for(k=0;k<grid->nx*grid->ny*grid->nz;k++) grid->zmin=MIN(grid->zmin,grid->z[k]);
    for(k=0;k<grid->nx*grid->ny*grid->nz;k++) grid->zmax=MAX(grid->zmax,grid->z[k]);
    grid->dz=(grid->zmax-grid->zmin)/(grid->nz-1);
    }

  if(vz==vt) {
    z0=grid->z[0];
    for(k=0;k<grid->nz;k++) {
      for(m=0;m<grid->nx*grid->ny;m++) grid->z[k*grid->nx*grid->ny+m]-=z0;
      }
    for(k=0;k<grid->nz;k++) {
      for(m=0;m<grid->nx*grid->ny;m++) grid->z[k*grid->nx*grid->ny+m]/=3600.0;
      }
    }
/*------------------------------------------------------------------------
  get s coordinate data in case s_dim ok*/
  if(vsigma!=-1) {
    *vdepth=vz; /*detected bathymetry*/
    var=vsigma;
    grid->z=NULL;
    switch (ndim[var]) {
      case 1:
/*------------------------------------------------------------------------
        assume uniform sigma levels */
        grid->sigma=(double *) malloc(grid->nz*sizeof(double));
        status=nc_get_var_double(ncid,vsigma,grid->sigma);
        if(status != NC_NOERR) goto error;
        grid->modeV=3;
/*------------------------------------------------------------------------
    load bathymetry*/
        var=vz;
        if(dimlgth[dimids[var][0]]*dimlgth[dimids[var][1]]==grid->nx*grid->ny) {
          grid->z=(double *) malloc(grid->nx*grid->ny*grid->nz*sizeof(double));
          zdum=(double *) malloc(dimlgth[dimids[var][0]]*dimlgth[dimids[var][1]]*sizeof(double));
          status=nc_get_var_double(ncid,var,zdum);
/*------------------------------------------------------------------------
          compute sigma level depth*/
          for(k=0;k<grid->nz;k++)
            for(m=0;m<grid->nx*grid->ny;m++) grid->z[k*grid->nx*grid->ny+m]=zdum[m]*grid->sigma[k];
          free(zdum);
          grid->zmin=1.e+10;
          grid->zmax=1.e-10;

          for(k=0;k<grid->nx*grid->ny*grid->nz;k++) grid->zmin=MIN(grid->zmin,grid->z[k]);
          for(k=0;k<grid->nx*grid->ny*grid->nz;k++) grid->zmax=MAX(grid->zmax,grid->z[k]);
          grid->dz=(grid->zmax-grid->zmin)/(grid->nz-1);
          }
        break;

      case 2:
/*------------------------------------------------------------------------
        structural anomaly */
        goto error;
        break;

      case 3:
/*------------------------------------------------------------------------
        x,y dependant s levels */
        grid->sigma=(double *) malloc(grid->nx*grid->ny*grid->nz*sizeof(double));
        status=nc_get_var_double(ncid,var,grid->sigma);
        if(status != NC_NOERR) goto error;
        grid->modeV=2;
        break;

      default:
/*------------------------------------------------------------------------
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

/*------------------------------------------------------------------------
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
  nc_check_error(status,__LINE__,__FILE__);
  if(status != NC_NOERR) goto close_error;

  status=0;
  return(status);

error:
  nc_check_error(status,__LINE__,__FILE__);
  if(dimlgth != NULL) free(dimlgth);
  status=nc_close(ncid);
  nc_check_error(status,__LINE__,__FILE__);
  status=-1;
  printf("loading 3D grid for variable %d failed\n",v_info.id);
  return(status);

no_go:
  nc_check_error(status,__LINE__,__FILE__);
  status=nc_close(ncid);
  nc_check_error(status,__LINE__,__FILE__);
  status=1; /* not netcdf file*/
  return(status);

unknown_format:
  status=nc_close(ncid);
  status=2; /* netcdf file, but not knwon format*/
  return(status);

close_error:
  nc_check_error(status,__LINE__,__FILE__);
  status=3; /* not a file?*/
  return(status);

open_error:
  nc_check_error(status,__LINE__,__FILE__);
  status=3; /* not a file?*/
  return(status);

  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int cdf_loadvargrid_3d (const char* filename,int vdata, grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int status;
  int verbose=0,vdepth;
  cdfgbl_t g_info;
  cdfvar_t v_info;
  
  status= cdf_globalinfo(filename, &g_info,verbose);
  status= cdf_varinfo(filename, vdata, &v_info);
  status= cdf_loadvar_gd_3d ((char*) filename, g_info, v_info, grid,&vdepth);
  
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int cdf_loadvar_gd_2d (const char* filename, cdfgbl_t global, cdfvar_t info, grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in;
  int i,j,k,k1,k2,m,option,status, ncid;
  size_t lengthhp,*dimlgth=NULL,index[1];
  size_t start[3], count[3];
  float time;
  int *nattsp;
  int d, dim,v,var,att,vx,vy,vz,vt,vmask,n;
  nc_type *vartype,xtypep;
  char *text;

  double *ddum=NULL,*vector=NULL,tmp,factor;
  int swap_needed=0;
  decoded_t decoded,vx_decoded,vy_decoded;
  cdfvar_t vx_info,vy_info;
  int verbose=0;
  date_t origin;
  cdfgbl_t grid_info;
  
  bool debug=false;
  
  vx=-1;
  vy=-1;

  decoded.vx=-1;
  decoded.vy=-1;

  status= cdf_globalinfo(filename, &grid_info,verbose);

/*------------------------------------------------------------------------
  check variable "associate" attribute */
  var=info.id;
  if(info.ndim <2) goto error;

  grid->circular=0; /*default*/
//  status= poc_decode_associates(info, global, &decoded);
/* *------------------------------------------------------------------------
  try first grid file, then data file (where time should be)*/
  status=poc_decode_associates(info,grid_info, &decoded,0,debug);
  status=((decoded.vx!=-1) && (decoded.vy!=-1));
  if(status==1) {
    vx=decoded.vx;
    vy=decoded.vy;
    status=cdf_varinfo(filename,grid_info.variable[vx].name,&vx_info,verbose);
    status=cdf_varinfo(filename,grid_info.variable[vy].name,&vy_info,verbose);
    goto next;
    }

  status=poc_decode_associates(info, global, &decoded,0,debug);
  status=((decoded.vx!=-1) && (decoded.vy!=-1));
  if(status==1) {
    vx=decoded.vx;
    vy=decoded.vy;
    status=cdf_varinfo(filename,global.variable[vx].name,&vx_info,verbose);
    status=cdf_varinfo(filename,global.variable[vy].name,&vy_info,verbose);
    }

next:

  if( (vx==-1) || (vy==-1) ) goto error;


/*------------------------------------------------------------------------
  open file for reading */
  status=nc_open(filename,0,&ncid);
  if(status != NC_NOERR) goto open_error;

  grid->nz=1;

  status=poc_decode_associates(vx_info,grid_info, &vx_decoded,1, debug);
  status=poc_decode_associates(vy_info,grid_info, &vy_decoded,1, debug);

  status=poc_decode_mask(vx_info, &vx_decoded);
  status=poc_decode_mask(vy_info, &vy_decoded);

  if((vx_info.ndim==2)&&(vx_info.ndim==2)) {
    if(vx_decoded.xlen!=vy_decoded.xlen) {
      printf("x and y have different size, abort...\n");
      goto error;
      }
    if(vx_decoded.ylen!=vy_decoded.ylen) {
      printf("x and y have different size, abort...\n");
      goto error;
      }
    }
/* *------------------------------------------------------------------------
  get dimension: decode from variable is not adequate (why ?)*/
  status= poc_decode_axis(info, global, &decoded);
//   grid->nx=vx_decoded.xlen;
//   grid->ny=vx_decoded.ylen;

  if(vx_decoded.xlen!=decoded.xlen) {
    printf("variable and grid have different x-size, abort...\n");
    goto error;
    }
  if(vy_decoded.ylen!=decoded.ylen) {
    printf("variable and grid have different y-size, abort...\n");
    goto error;
    }

  grid->nx=decoded.xlen;
  grid->ny=decoded.ylen;
  grid->nz=1;

  ddum=(double *) malloc(grid->nx*grid->ny*sizeof(double));

/*------------------------------------------------------------------------
  get x coordinate data */
  var=vx;
  switch (vx_info.ndim) {
    case 3:
      count[0]=1;
      count[1]=decoded.ylen;
      count[2]=decoded.xlen;
      start[0]=0;
      start[1]=0;
      start[2]=0;
      status=nc_get_vara_double(ncid,var,start,count,ddum);
      if(status != NC_NOERR) goto error;
      grid->modeH=2;
      break;

    case 2:
      status=nc_get_var_double(ncid,var,ddum);
      if(status != NC_NOERR) goto error;
      grid->modeH=2;
      break;

    case 1:
      vector=(double *) malloc(grid->nx*sizeof(double));
      status=nc_get_var_double(ncid,var,vector);
      if(status != NC_NOERR) goto error;
      for(j=0;j<grid->ny;j++)
        for(i=0;i<grid->nx;i++) ddum[grid->nx*j+i]=vector[i];
      free(vector);
      grid->modeH=2;
      break;

    default:
      printf("dimension of x=%d ???\n",vx_info.ndim);
      goto error;
      break;
    }

  if(vx_decoded.offset!=0) {
    for(n=0;n<decoded.ylen*decoded.xlen;n++) ddum[n]+=vx_decoded.offset;
    }
  if(vx_decoded.scale!=1) {
    for(n=0;n<decoded.ylen*decoded.xlen;n++) ddum[n]*=vx_decoded.scale;
    }

  grid->xmin=ddum[0];
  grid->xmax=ddum[grid->nx-1];

  factor=1.;

  status=nc_inq_att(ncid,var,"units",&xtypep,&lengthhp);
  if(status == NC_NOERR) {
    text=(char *) malloc(lengthhp+1);
    status=nc_get_att_text(ncid,var,"units",text);
    text[lengthhp]='\0';
    if(strcmp(text,"radian_east") ==0) factor=180./M_PI;;
    if(strcmp(text,"degree_east")==0) factor=1.;
    if(strcmp(text,"degree_east") ==0) factor=1.;
    free(text);
    }

  status=nc_inq_att(ncid,var,"topology",&xtypep,&lengthhp);
  if(status == NC_NOERR) {
    text=(char *) malloc(lengthhp+1);
    status=nc_get_att_text(ncid,var,"topology",text);
    text[lengthhp]='\0';
    if(strcmp(text,"circular")==0) grid->circular=1;
    free(text);
    }

  grid->x=(double *) malloc(grid->nx*grid->ny*sizeof(double));
  for(k=0;k<grid->nx*grid->ny;k++) grid->x[k]=ddum[k]*factor;
  for(k=0;k<grid->nx*grid->ny;k++) grid->xmin=MIN(grid->xmin,grid->x[k]);
  for(k=0;k<grid->nx*grid->ny;k++) grid->xmax=MAX(grid->xmax,grid->x[k]);
  grid->dx=(grid->xmax-grid->xmin)/(grid->nx-1);


/*------------------------------------------------------------------------
  get y coordinate data */
  var=vy;
  switch (vy_info.ndim) {
    case 3:
      count[0]=1;
      count[1]=decoded.ylen;
      count[2]=decoded.xlen;
      start[0]=0;
      start[1]=0;
      start[2]=0;
      status=nc_get_vara_double(ncid,var,start,count,ddum);
      if(status != NC_NOERR) goto error;
      grid->modeH=2;
      break;

    case 2:
    status=nc_get_var_double(ncid,var,ddum);
    if(status != NC_NOERR) goto error;
    grid->modeH=2;
    break;

    case 1:
    vector=(double *) malloc(grid->ny*sizeof(double));
    status=nc_get_var_double(ncid,var,vector);
    if(status != NC_NOERR) goto error;
    for(j=0;j<grid->ny;j++)
      for(i=0;i<grid->nx;i++) ddum[grid->nx*j+i]=vector[j];
    free(vector);
    grid->modeH=2;
    break;

    default:
    printf("dimension of y=%d ???\n",vy_info.ndim);
    goto error;
    break;
    }

  if(vy_decoded.offset!=0) {
    for(n=0;n<decoded.ylen*decoded.xlen;n++) ddum[n]+=vy_decoded.offset;
    }
  if(vy_decoded.scale!=1) {
    for(n=0;n<decoded.ylen*decoded.xlen;n++) ddum[n]*=vy_decoded.scale;
    }
    
  grid->ymin=ddum[0];
  grid->ymax=ddum[grid->nx*grid->ny-1];

  factor=1.;

  status=nc_inq_att(ncid,var,"units",&xtypep,&lengthhp);
  if(status == NC_NOERR) {
    text=(char *) malloc(lengthhp+1);
    status=nc_get_att_text(ncid,var,"units",text);
    text[lengthhp]='\0';
    if(strcmp(text,"radian_north")==0) factor=180./M_PI;
    if(strcmp(text,"degree_north")==0) factor=1.;
    if(strcmp(text,"degree_north")==0) factor=1.;
    free(text);
    }

  grid->y=(double *) malloc(grid->nx*grid->ny*sizeof(double));
  for(k=0;k<grid->nx*grid->ny;k++) grid->y[k]=ddum[k]*factor;
  for(k=0;k<grid->nx*grid->ny;k++) grid->ymin=MIN(grid->ymin,grid->y[k]);
  for(k=0;k<grid->nx*grid->ny;k++) grid->ymax=MAX(grid->ymax,grid->y[k]);
  grid->dy=(grid->ymax-grid->ymin)/(grid->ny-1);

  free(ddum);

  grid->zmin=0.;
  grid->zmax=0.;
  grid->z=NULL;
  grid->nz=1;
  grid->dz=0.0;

  grid->modeH= 2;
  grid->modeV=-1;


/*------------------------------------------------------------------------
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

  status=nc_close(ncid);
  nc_check_error(status,__LINE__,__FILE__);
  if(status != NC_NOERR) goto close_error;

//  status= poc_gettime   (filename,  info,  global, &origin, &(grid->time), &(grid->nt));

  status=0;
  return(status);

  error:
//   check_err(status,__LINE__,__FILE__);
  if(dimlgth != NULL) free(dimlgth);
  status=nc_close(ncid);
//   check_err(status,__LINE__,__FILE__);
  status=-1;
  printf("loading 2D grid for variable %d failed\n",info.id);
  return(status);

  no_go:
//   check_err(status,__LINE__,__FILE__);
  status=nc_close(ncid);
//   check_err(status,__LINE__,__FILE__);
  status=1; /* not netcdf file*/
  return(status);

  unknown_format:
  status=nc_close(ncid);
  status=2; /* netcdf file, but not knwon format*/
  return(status);

  close_error:
//   check_err(status,__LINE__,__FILE__);
  status=3; /* not a file?*/
  return(status);

  open_error:
//   check_err(status,__LINE__,__FILE__);
  status=3; /* not a file?*/
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cdf_loadvargrid (const char* filename, int vdata, grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int verbose=0;
  cdfgbl_t g_info;
  cdfvar_t v_info;
  
  status= cdf_globalinfo(filename, &g_info,verbose);
  status= cdf_varinfo(filename, vdata, &v_info);
  status= cdf_loadvar_gd_2d (filename, g_info, v_info, grid);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cdf_loadvargrid_2d (const char* filename,int vdata, grid_t *grid,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  cdfgbl_t g_info;
  cdfvar_t v_info;
  
  status= cdf_globalinfo(filename, &g_info,verbose);
  status= cdf_varinfo(filename, vdata, &v_info);
  status= cdf_loadvar_gd_2d (filename, g_info, v_info, grid);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cdf_loadvar_r1_2d (const char* filename, int varid, int level, int frame, grid_t grid, float *buf, float *mask,variable_t *info, decoded_t & decoded)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int option,status, ncid;
  int axis_id,n;
  int v,vmask=-1,vtime=-1;
  int i,j,k,k1,k2,m,fmt,swap_needed=0,d,missing_found=0;
  char name[1024],*long_name=NULL,*associate=NULL,*maskname=NULL;
  size_t *start,*count;
  float fill,missing=1.e+10,scale=1.0,offset=0.0,*tmp;
  int tdim,nframes;
  double *time,factor;
  int verbose=0,chk=0;
  cdfgbl_t g_info;
  cdfvar_t v_info;
//  decoded_t decoded;
  date_t origine;

  info->name=NULL;
  info->standard_name=NULL;
  info->long_name=NULL;

  status=nc_open(filename,0,&ncid);
  if (status !=NC_NOERR) {
    status=-1;
    return(status);
    }

  status= cdf_globalinfo(ncid,&g_info,verbose);
  if(status!=0) goto error;

  status= cdf_varinfo(ncid,varid,&v_info,verbose);
  if(status!=0) goto error;

/*------------------------------------------------------------------------
  Find dimensions order: use "AXIS" attribute if available */
  axis_id=-1;
  for(n=0;n<v_info.natt;n++) {
    if(strcmp(v_info.att[n].name,"axis")==0) {
      axis_id=n;
      break;
      }
/**------------------------------------------------------------------------
    attribute name changed to comply with CF standard*/
    if(strcmp(v_info.att[n].name,"content")==0) {
      axis_id=n;
      break;
      }
    }

  status=poc_decode_axis(v_info, g_info, &decoded);

  if(status!=0) goto error;

  if(decoded.xlen!=grid.nx) {
    printf("variable and axis have different X-size (%d %d), abort...\n",decoded.xlen,grid.nx);
    goto error;
    }
  if(decoded.ylen!=grid.ny) {
    printf("variable and axis have different Y-size (%d %d), abort...\n",decoded.ylen,grid.ny);
    goto error;
    }

  vmask=-1;

  start=new size_t[v_info.ndim];
  count=new size_t[v_info.ndim];

  start[decoded.xdim]=0;
  count[decoded.xdim]=decoded.xlen;

  start[decoded.ydim]=0;
  count[decoded.ydim]=decoded.ylen;

  if(decoded.zdim!=-1){
    start[decoded.zdim]=level-1;
    count[decoded.zdim]=1;
    }

  if(decoded.tdim!=-1){
/**----------------------------------------------------------------------------
    warning: frame index start at 1, netcdf at 0*/
//    start[decoded.tdim]=frame-1;
    start[decoded.tdim]=frame;
    count[decoded.tdim]=1;
    }

  if(decoded.fdim!=-1){
/**----------------------------------------------------------------------------
    warning: frame index start at 1, netcdf at 0*/
//    start[decoded.fdim]=level-1;
    start[decoded.fdim]=level;
    count[decoded.fdim]=1;
    }

  status=nc_get_vara_float(ncid,varid,start,count,buf);
  if(status != NC_NOERR) goto error;

/*------------------------------------------------------------------------
  patch: force missing value from mask if exists */
  status= poc_decode_mask(v_info, &decoded);
  *mask=decoded.spec;
  scale=decoded.scale;
  offset=decoded.offset;

/*------------------------------------------------------------------------
  */
  if( (isnan(*mask)) || (isinf(*mask)) ){
    for(k=0;k<grid.nz;k++) {
      for(j=0;j<grid.ny;j++) {
        for(i=0;i<grid.nx;i++) {
          m=k*grid.ny*grid.nx+j*grid.nx+i;
          if((isnan(buf[m])) || (isinf(buf[m])) ) {
            buf[m]=1.e+10;
            }
          }
        }
      }
    *mask=1.e+10;
    }

/*------------------------------------------------------------------------
  aply scale and offset */
  if((scale!=1.) || (offset!=0))
    for(k=0;k<grid.nz;k++)
      for(j=0;j<grid.ny;j++)
        for(i=0;i<grid.nx;i++) {
          m=k*grid.ny*grid.nx+j*grid.nx+i;
          if(isnan(buf[m])) {
            printf("Nan value out of mask\n");
            }
          if(buf[m]!=*mask) {
            buf[m]=buf[m]*scale+offset;
            }
          else {
            buf[m]=*mask;
            }
          }

/*------------------------------------------------------------------------
  figure out time of frame */
  status= poc_gettime(ncid, v_info, g_info, &origine, &time, &nframes);
  if(status==0) {
/*------------------------------------------------------------------------
    time origin */
//    info->origin=origine;
    info->time=time[frame-1];
    delete[] time;
    }
  else
    info->time=0;

  status=nc_close(ncid);
  nc_check_error(status,__LINE__,__FILE__);

  info->name   =strdup(v_info.name);
  info->origin =strdup(filename);
  info->ndim   =v_info.ndim;
  info->id     =varid;

  status= poc_decode_names(v_info, &decoded);
  status= poc_decode_names(g_info, &decoded);

  info->standard_name=NULL;
  info->long_name=NULL;
  info->production=NULL;

  if(decoded.standard_name==0) decoded.standard_name=strdup("no_standard_name");
  if(decoded.long_name==0)     decoded.long_name    =strdup("no_long_name");
  if(decoded.production==0)    decoded.production   =strdup("");

  if(strlen(decoded.standard_name) != 0)
    info->standard_name=strdup(decoded.standard_name);
  if(strlen(decoded.long_name) != 0) {
    info->long_name=strdup(decoded.long_name);
    }
  else {
    info->long_name=strdup(info->name);
    }
  
  if(strlen(decoded.production) != 0)
    info->production=strdup(decoded.production);

  if(info->standard_name==NULL) info->standard_name=strdup(info->long_name);
  for(k=0;k<strlen(info->standard_name);k++)
    if(info->standard_name[k]==' ') info->standard_name[k]='_';

  info->units=NULL;
  if(decoded.units != 0)
    info->units=strdup(decoded.units);
  if(info->units==NULL) info->units=strdup("no units");

  switch (v_info.ndim){
    case 2:
      swap_needed=strcmp(decoded.axis,"XY")==0;
      break;
    case 3:
      swap_needed=strcmp(decoded.axis,"TXY")==0;
      break;
    }
  if(swap_needed) {
/*------------------------------------------------------------------------
    X,Y,Z dimension are in the wrong order */
    tmp=new float[grid.nx*grid.ny];
    for(j=0;j<grid.ny;j++)
      for(i=0;i<grid.nx;i++) {
        k1=j*grid.nx+i;
        k2=i*grid.ny+j;
        tmp[k1]=buf[k2];
        }
    for(j=0;j<grid.ny;j++)
      for(i=0;i<grid.nx;i++) {
        k1=j*grid.nx+i;
        buf[k1]=tmp[k1];
        }
    delete[] tmp;
    }

  status=0;
  return(status);

error:
//  check_err(status,__LINE__,__FILE__);
  status=nc_close(ncid);
  status=-1;
  printf("loading variable %d failed\n",varid);
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int cdf_loadvar_r1_2d (const char* filename, int v, int k, int t, grid_t grid, int n, float *buf, float *mask,variable_t *info)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  decoded_t decoded;
  int status;
  
  status=cdf_loadvar_r1_2d (filename,  v,  k,  t, grid, buf, mask,info,decoded);
  return(status);
  
}
