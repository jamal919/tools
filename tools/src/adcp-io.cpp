
/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Cyril Nguen        LA, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#include "config.h"

#include <stdio.h>
#include <string.h>

#include "tools-structures.h"

#include "tides.h"
#include "fe.h"
#include "archive.h"
#include "netcdf-proto.h"
#include "poc-time.h"
#include "mgr.h"
#include "functions.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int sadcp_info(char *filename, double **time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int ncid,status;
  int time_id,reference_id;
  int vdim,tdim;
  int u_id,v_id,w_id;
  cdfvar_t u_info,v_info,w_info,reference_info;
  cdfgbl_t global;
  char *s;
  date_t *origin;
  int hour,minute,seconds;

  origin=new date_t;

  int nlevels, nsamples;

  status= cdf_globalinfo(filename,&global,0);

  status = nc_open(filename, NC_NOWRITE, &ncid);

  status = nc_inq_varid(ncid, "JULD", &time_id);
  status = nc_inq_varid(ncid, "REFERENCE_DATE_TIME", &reference_id);

  status = nc_inq_varid(ncid, "UVEL_ADCP", &u_id);
  status = nc_inq_varid(ncid, "VVEL_ADCP", &v_id);
  status = nc_inq_varid(ncid, "WVEL_ADCP", &w_id);

  vdim=cdf_identify_dimension(global,"N_LEVEL");
  tdim=cdf_identify_dimension(global,"N_DATE_TIME");

  status= cdf_varinfo(ncid, u_id, &u_info);

  s=new char[32];
  status=nc_get_var_text(ncid,reference_id,s);
  nc_check_error(status,__LINE__,__FILE__);

  sscanf(s,"%4d%2d%2d%2d%2d%2d",&(origin->year),&(origin->month),&(origin->day),&hour,&minute,&seconds);
  origin->second=hour*3600.+minute*60.+seconds;

  nlevels=global.dimension[vdim].length;
  nsamples=global.dimension[tdim].length;

  *time=new double[nsamples];
  status=nc_get_var_double(ncid,time_id,*time);
  nc_check_error(status,__LINE__,__FILE__);

//  initial = time_getcnesdate(first, 's');
  status = nc_close(ncid);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int sadcp_read(char *filename, float **u, float **v, int level)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int ncid,status;
  int time_id,reference_id;
  int vdim,tdim;
  int u_id,v_id,w_id;
  cdfvar_t u_info,v_info,w_info,reference_info;
  cdfgbl_t global;
  char *s;
  date_t *origin;
  int hour,minute,seconds;
  size_t count[2],start[2];
  ptrdiff_t stride[2];

  origin=new date_t;

  int nlevels, nsamples;

  status= cdf_globalinfo(filename,&global,0);

  status = nc_open(filename, NC_NOWRITE, &ncid);

  status = nc_inq_varid(ncid, "JULD", &time_id);
  status = nc_inq_varid(ncid, "REFERENCE_DATE_TIME", &reference_id);

  status = nc_inq_varid(ncid, "UVEL_ADCP", &u_id);
  status = nc_inq_varid(ncid, "VVEL_ADCP", &v_id);
  status = nc_inq_varid(ncid, "WVEL_ADCP", &w_id);

  vdim=cdf_identify_dimension(global,"N_LEVEL");
  tdim=cdf_identify_dimension(global,"N_DATE_TIME");

  nlevels=global.dimension[vdim].length;
  nsamples=global.dimension[tdim].length;

  status= cdf_varinfo(ncid, u_id, &u_info);

  *u=new float[nsamples];
  *v=new float[nsamples];
//  *w=new float[nsamples];

  count[0]=nsamples;
  count[1]=1;

  start[0]=0;
  start[1]=level;

  stride[0]=1;
  stride[1]=nlevels;

  *u[0]=1;

  status=nc_get_vars_float(ncid,u_id,start,count,stride,*u);
  nc_check_error(status,__LINE__,__FILE__);

  status=nc_get_vars_float(ncid,v_id,start,count,stride,*v);
  nc_check_error(status,__LINE__,__FILE__);

//  initial = time_getcnesdate(first, 's');
  status = nc_close(ncid);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int sadcp_readall(char *filename, float **u, float **v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int ncid,status;
  int time_id,reference_id;
  int vdim,tdim;
  int u_id,v_id,w_id;
  cdfvar_t u_info,v_info,w_info,reference_info;
  cdfgbl_t global;
  char *s;
  date_t *origin;
  int hour,minute,seconds;
  size_t count[2],start[2];
  ptrdiff_t stride[2];

  origin=new date_t;

  int nlevels, nsamples;

  status= cdf_globalinfo(filename,&global,0);

  status = nc_open(filename, NC_NOWRITE, &ncid);

  status = nc_inq_varid(ncid, "JULD", &time_id);
  status = nc_inq_varid(ncid, "REFERENCE_DATE_TIME", &reference_id);

  status = nc_inq_varid(ncid, "UVEL_ADCP", &u_id);
  status = nc_inq_varid(ncid, "VVEL_ADCP", &v_id);
  status = nc_inq_varid(ncid, "WVEL_ADCP", &w_id);

  vdim=cdf_identify_dimension(global,"N_LEVEL");
  tdim=cdf_identify_dimension(global,"N_DATE_TIME");

  nlevels=global.dimension[vdim].length;
  nsamples=global.dimension[tdim].length;

  status= cdf_varinfo(ncid, u_id, &u_info);

  *u=new float[nsamples*nlevels];
  *v=new float[nsamples*nlevels];
//  *w=new float[nsamples];

  count[0]=nsamples;
  count[1]=nlevels;

  start[0]=0;
  start[1]=0;

  status=nc_get_vara_float(ncid,u_id,start,count,*u);
  nc_check_error(status,__LINE__,__FILE__);

  status=nc_get_vara_float(ncid,v_id,start,count,*v);
  nc_check_error(status,__LINE__,__FILE__);

//  initial = time_getcnesdate(first, 's');
  status=nc_get_var_float(ncid,u_id,*u);
  nc_check_error(status,__LINE__,__FILE__);

  status = nc_close(ncid);

}

