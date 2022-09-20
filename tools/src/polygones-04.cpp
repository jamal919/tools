

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

*/

#include "poc-netcdf.hpp"
#include "gshhs.h"
#include "polygones.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_write_netcdf (const char *filename, const plg_t *polygones, int nplg)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,verbose=0;
  int ncid,lon_id,lat_id,start_id,count_id,segments_id;
  int k,n,ns,nsmax,npts;
  double lon, lat;
  cdfgbl_t global;
  poc_dim_t dimL, dimP;
  int dimidL,dimidP;
  size_t start[1], count[1];
  poc_var_t vlon,vlat,vstart,vcount;
//  size_t *buffer;
  int *buffer;

  npts=0;
  for(n=0;n<nplg;n++) {
    npts+=polygones[n].npt;
    }
   
  status=nc_create(filename, NC_CLOBBER, &ncid);
  
  dimL.len=nplg;
  dimL.name="nlines";
  
  status=poc_def_dim(ncid,(const poc_dim_t) dimL,&dimidL);
  
  dimP.len=npts;
  dimP.name="npoints";
  
  status=poc_def_dim(ncid,(const poc_dim_t) dimP,&dimidP);
  
  vlon.name="lon";
  vlon.type=NC_DOUBLE;
  vlon.dimensions<<dimP;
  
  vlat.name="lat";
  vlat.type=NC_DOUBLE;
  vlat.dimensions<<dimP;
  
  status=poc_def_var(ncid,vlon,&lon_id);
  status=poc_def_var(ncid,vlat,&lat_id);
  
  vstart.name="start";
  vstart.type=NC_INT;
  vstart.dimensions<<dimL;

  vcount.name="count";
  vcount.type=NC_INT;
  vcount.dimensions<<dimL;
  
  status=poc_def_var(ncid,vstart,&start_id);
  status=poc_def_var(ncid,vcount,&count_id);
  
//   status=nc_inq_varid(ncid,"lon",&lon_id);
//   status=nc_inq_varid(ncid,"lat",&lat_id);
//   status=nc_inq_varid(ncid,"segments",&segments_id);
  status=nc_close(ncid);
  
  status=nc_open(filename, NC_WRITE, &ncid);
  
  start[0]=0;
  for(ns=0;ns<nplg;ns++) {
    count[0]=polygones[ns].npt;
    nc_put_vara_double(ncid, lon_id, start, count,polygones[ns].t);
    nc_put_vara_double(ncid, lat_id, start, count,polygones[ns].p);
    start[0]+=count[0];
    }

  start[0]=0;
  for(ns=0;ns<nplg;ns++) {
    count[0]=polygones[ns].npt;
    nc_put_vara_double(ncid, lon_id, start, count,polygones[ns].t);
    nc_put_vara_double(ncid, lat_id, start, count,polygones[ns].p);
    start[0]+=count[0];
    }
    
  buffer=new int[nplg];
//  buffer=new size_t[nplg];
  buffer[0]=0;
  for(ns=1;ns<nplg;ns++) {
    buffer[ns]=buffer[ns-1]+polygones[ns-1].npt;
    }
  nc_put_var_int(ncid, start_id, buffer);
//  nc_put_var_ulonglong(ncid, start_id, buffer);
  
  for(ns=0;ns<nplg;ns++) {
    buffer[ns]=polygones[ns].npt;
    }
  nc_put_var_int(ncid, count_id, buffer);
//  nc_put_var_ulonglong(ncid, count_id, buffer);
  
  status=nc_close(ncid);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_read_netcdf (const char *filename, plg_t *polygones, int nplg)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,verbose=0;
  int ncid,lon_id,lat_id,start_id,count_id;
  int k,ns,nsmax,npts,nlines;
  double lon, lat;
  size_t start[1], count[1];
  int *buffer;
  
  status=nc_open(filename, NC_NOWRITE, &ncid);
  
  status=nc_inq_varid(ncid,"lon",&lon_id);
  status=nc_inq_varid(ncid,"lat",&lat_id);
  
  status=nc_inq_varid(ncid,"start",&start_id);
  status=nc_inq_varid(ncid,"count",&count_id);
  
  buffer=new int[nplg];
  nc_get_var_int(ncid, count_id, buffer);
  for(ns=0;ns<nplg;ns++) {
    polygones[ns].npt=buffer[ns];
    }
  deletep(&buffer);
  
  start[0]=0;
  for(ns=0;ns<nplg;ns++) {
    polygones[ns].t= new double[polygones[ns].npt];
    polygones[ns].p= new double[polygones[ns].npt];
    polygones[ns].x= new double[polygones[ns].npt];
    polygones[ns].y= new double[polygones[ns].npt];
    count[0]=polygones[ns].npt;
    nc_get_vara_double(ncid, lon_id, start, count,polygones[ns].t);
    nc_get_vara_double(ncid, lat_id, start, count,polygones[ns].p);
    nc_get_vara_double(ncid, lon_id, start, count,polygones[ns].x);
    nc_get_vara_double(ncid, lat_id, start, count,polygones[ns].y);
    start[0]+=count[0];
    }

  status=nc_close(ncid);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_inquire_netcdf (const char *filename, int *ns)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,verbose=0;
  int ncid,lon_id,lat_id,start_id,count_id;
  cdfvar_t v;
  
  status=nc_open(filename, NC_NOWRITE, &ncid);
  if(status!=0) return(-1);
  
  status=nc_inq_varid(ncid,"lon",&lon_id);
  if(status!=0) return(-1);
  status=nc_inq_varid(ncid,"lat",&lat_id);
  if(status!=0) return(-1);
  
  status=nc_inq_varid(ncid,"start",&start_id);
  if(status!=0) return(-1);
  status=nc_inq_varid(ncid,"count",&count_id);
  if(status!=0) return(-1);
  
  status= cdf_varinfo(ncid,start_id, &v,verbose);
 
  *ns=v.dim[0].length;
  status=nc_close(ncid);
  
  v.destroy();

  return(0);
}

