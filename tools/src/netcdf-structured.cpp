
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

\brief data input/output poc-netcdf definitions
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>

#include "tools-structures.h"
#include "constants.h"
#include "netcdf-proto.h"
#include "poc-netcdf.def"


void delete_cdfvar_p(cdfvar_t **var)
{
  if(*var==0) return;
  (*var)->destroy();
  delete[]*var;
  *var=0;
}


/*----------------------------------------------------------------------------*/


int create_ncfile2d(char *filename,size_t x_len, size_t y_len, grid_t grid,int options[create_ncfile2d_OPTIONS])


/*----------------------------------------------------------------------------*/
/* *------------------------------------------------------------------------
  attribute name "axis" changed in "content" to comply with CF standard*/
{
   int  ncid;	/* netCDF id */

   /* dimension ids */
   int x_dim;
   int y_dim;
   int time_counter_dim;

   /* dimension lengths */
/*
   size_t x_len = 100;
   size_t y_len = 100;
*/
   size_t time_counter_len = NC_UNLIMITED;

   /* variable ids */

   int lon_id;
   int lat_id;
   int time_id;
   int h_id;
   int u_id;
   int v_id;
   int p_id;
   int ibd_id;

   /* variable shapes */
   int lon_dims[RANK_lon];
   int lat_dims[RANK_lat];
   int time_dims[RANK_time];
   int h_dims[RANK_h];
   int u_dims[RANK_u];
   int v_dims[RANK_v];
   int p_dims[RANK_p];
   int ibd_dims[RANK_ibd];

   /* attribute vectors */
   double lon_valid_min[1];
   double lon_valid_max[1];
   double lat_valid_min[1];
   double lat_valid_max[1];
   float h_missing_value[1];
   float h__FillValue[1];
   float h_scale_factor[1];
   float h_add_offset[1];
   float u_missing_value[1];
   float u__FillValue[1];
   float u_scale_factor[1];
   float u_add_offset[1];
   float v_missing_value[1];
   float v__FillValue[1];
   float v_scale_factor[1];
   float v_add_offset[1];
   float p_missing_value[1];
   float p__FillValue[1];
   float p_scale_factor[1];
   float p_add_offset[1];
   float ibd_missing_value[1];
   float ibd__FillValue[1];
   float ibd_scale_factor[1];
   float ibd_add_offset[1];

   /* enter define mode */
   int stat = nc_create(filename, NC_CLOBBER, &ncid);
   nc_check_error(stat,__LINE__,__FILE__);

   /* define dimensions */
   stat = nc_def_dim(ncid, "x", x_len, &x_dim);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_def_dim(ncid, "y", y_len, &y_dim);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_def_dim(ncid, "time_counter", time_counter_len, &time_counter_dim);
   nc_check_error(stat,__LINE__,__FILE__);

   /* define variables */

   lon_dims[0] = y_dim;
   lon_dims[1] = x_dim;
   stat = nc_def_var(ncid, "lon", NC_DOUBLE, RANK_lon, lon_dims, &lon_id);
   nc_check_error(stat,__LINE__,__FILE__);

   lat_dims[0] = y_dim;
   lat_dims[1] = x_dim;
   stat = nc_def_var(ncid, "lat", NC_DOUBLE, RANK_lat, lat_dims, &lat_id);
   nc_check_error(stat,__LINE__,__FILE__);

   time_dims[0] = time_counter_dim;
   stat = nc_def_var(ncid, "time", NC_DOUBLE, RANK_time, time_dims, &time_id);
   nc_check_error(stat,__LINE__,__FILE__);

   if(options[ELEVATION]==1) {
     h_dims[0] = time_counter_dim;
     h_dims[1] = y_dim;
     h_dims[2] = x_dim;
     stat = nc_def_var(ncid, "h", NC_FLOAT, RANK_h, h_dims, &h_id);
     nc_check_error(stat,__LINE__,__FILE__);
     }

   if(options[CURRENTS]==1) {
     u_dims[0] = time_counter_dim;
     u_dims[1] = y_dim;
     u_dims[2] = x_dim;
     stat = nc_def_var(ncid, "u", NC_FLOAT, RANK_u, u_dims, &u_id);
     nc_check_error(stat,__LINE__,__FILE__);

     v_dims[0] = time_counter_dim;
     v_dims[1] = y_dim;
     v_dims[2] = x_dim;
     stat = nc_def_var(ncid, "v", NC_FLOAT, RANK_v, v_dims, &v_id);
     nc_check_error(stat,__LINE__,__FILE__);
     }

   if(options[PRESSURE]==1) {
     p_dims[0] = time_counter_dim;
     p_dims[1] = y_dim;
     p_dims[2] = x_dim;
     stat = nc_def_var(ncid, "p", NC_FLOAT, RANK_p, p_dims, &p_id);
     nc_check_error(stat,__LINE__,__FILE__);
     }

   if(options[IBD]==1) {
     ibd_dims[0] = time_counter_dim;
     ibd_dims[1] = y_dim;
     ibd_dims[2] = x_dim;
     stat = nc_def_var(ncid, "ibd", NC_FLOAT, RANK_ibd, ibd_dims, &ibd_id);
     nc_check_error(stat,__LINE__,__FILE__);
     }

   /* assign attributes */
   stat = nc_put_att_text(ncid, lon_id, "units", 12, "degree_east");
   nc_check_error(stat,__LINE__,__FILE__);
   lon_valid_min[0] = -180;
   stat = nc_put_att_double(ncid, lon_id, "valid_min", NC_DOUBLE, 1, lon_valid_min);
   nc_check_error(stat,__LINE__,__FILE__);
   lon_valid_max[0] = 180;
   stat = nc_put_att_double(ncid, lon_id, "valid_max", NC_DOUBLE, 1, lon_valid_max);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lon_id, "long_name", 9, "longitude");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lon_id, "standard_name", 9, "longitude");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lon_id, "nav_model", 12, "Default grid");
   nc_check_error(stat,__LINE__,__FILE__);

   stat = nc_put_att_text(ncid, lat_id, "units", 13, "degree_north");
   nc_check_error(stat,__LINE__,__FILE__);
   lat_valid_min[0] = -90;
   stat = nc_put_att_double(ncid, lat_id, "valid_min", NC_DOUBLE, 1, lat_valid_min);
   nc_check_error(stat,__LINE__,__FILE__);
   lat_valid_max[0] = 90;
   stat = nc_put_att_double(ncid, lat_id, "valid_max", NC_DOUBLE, 1, lat_valid_max);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lat_id, "long_name", 8, "latitude");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lat_id, "standard_name", 8, "latitude");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lat_id, "nav_model", 12, "Default grid");
   nc_check_error(stat,__LINE__,__FILE__);

   stat = nc_put_att_text(ncid, time_id, "units", 7, "seconds");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, time_id, "calendar", 9, "gregorian");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, time_id, "title", 4, "Time");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, time_id, "long_name", 34, "Time axis elapsed from time_origin");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, time_id, "time_origin", 21, " 1950-JAN-01 00:00:00");
   nc_check_error(stat,__LINE__,__FILE__);

   if(options[ELEVATION]==1) {
     stat = nc_put_att_text(ncid, h_id, "units", 1, "m");
     nc_check_error(stat,__LINE__,__FILE__);
     h_missing_value[0] = 1e+35;
     stat = nc_put_att_float(ncid, h_id, "missing_value", NC_FLOAT, 1, h_missing_value);
     nc_check_error(stat,__LINE__,__FILE__);
     h__FillValue[0] = 1e+35;
     stat = nc_put_att_float(ncid, h_id, "_FillValue", NC_FLOAT, 1, h__FillValue);
     nc_check_error(stat,__LINE__,__FILE__);
     h_scale_factor[0] = 1;
     stat = nc_put_att_float(ncid, h_id, "scale_factor", NC_FLOAT, 1, h_scale_factor);
     nc_check_error(stat,__LINE__,__FILE__);
     h_add_offset[0] = 0;
     stat = nc_put_att_float(ncid, h_id, "add_offset", NC_FLOAT, 1, h_add_offset);
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, h_id, "long_name", 19, "sea level elevation");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, h_id, "standard_name", 21, "sea_surface_elevation");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, h_id, "float_name", 1, "h");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, h_id, "content", 3, "TYX");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, h_id, "associate", 12, "time lat lon");
     nc_check_error(stat,__LINE__,__FILE__);
     }

   if(options[CURRENTS]==1) {
     stat = nc_put_att_text(ncid, u_id, "units", 5, "m s-1");
     nc_check_error(stat,__LINE__,__FILE__);
     u_missing_value[0] = 1e+35;
     stat = nc_put_att_float(ncid, u_id, "missing_value", NC_FLOAT, 1, u_missing_value);
     nc_check_error(stat,__LINE__,__FILE__);
     u__FillValue[0] = 1e+35;
     stat = nc_put_att_float(ncid, u_id, "_FillValue", NC_FLOAT, 1, u__FillValue);
     nc_check_error(stat,__LINE__,__FILE__);
     u_scale_factor[0] = 1;
     stat = nc_put_att_float(ncid, u_id, "scale_factor", NC_FLOAT, 1, u_scale_factor);
     nc_check_error(stat,__LINE__,__FILE__);
     u_add_offset[0] = 0;
     stat = nc_put_att_float(ncid, u_id, "add_offset", NC_FLOAT, 1, u_add_offset);
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, u_id, "long_name", 28, "ocean barotropic E-component");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, u_id, "standard_name", 39, "b arotropic_eastward_sea_water_velocity");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, u_id, "float_name", 1, "u");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, u_id, "content", 3, "TYX");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, u_id, "associate", 12, "time lat lon");
     nc_check_error(stat,__LINE__,__FILE__);

     stat = nc_put_att_text(ncid, v_id, "units", 5, "m s-1");
     nc_check_error(stat,__LINE__,__FILE__);
     v_missing_value[0] = 1e+35;
     stat = nc_put_att_float(ncid, v_id, "missing_value", NC_FLOAT, 1, v_missing_value);
     nc_check_error(stat,__LINE__,__FILE__);
     v__FillValue[0] = 1e+35;
     stat = nc_put_att_float(ncid, v_id, "_FillValue", NC_FLOAT, 1, v__FillValue);
     nc_check_error(stat,__LINE__,__FILE__);
     v_scale_factor[0] = 1;
     stat = nc_put_att_float(ncid, v_id, "scale_factor", NC_FLOAT, 1, v_scale_factor);
     nc_check_error(stat,__LINE__,__FILE__);
     v_add_offset[0] = 0;
     stat = nc_put_att_float(ncid, v_id, "add_offset", NC_FLOAT, 1, v_add_offset);
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, v_id, "long_name", 28, "ocean barotropic N-component");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, v_id, "standard_name", 39, "b arotropic_eastward_sea_water_velocity");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, v_id, "float_name", 1, "v");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, v_id, "content", 3, "TYX");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, v_id, "associate", 12, "time lat lon");
     nc_check_error(stat,__LINE__,__FILE__);
     }

   if(options[PRESSURE]==1) {
     stat = nc_put_att_text(ncid, p_id, "units", 2, "Pa");
     nc_check_error(stat,__LINE__,__FILE__);
     p_missing_value[0] = 1e+35;
     stat = nc_put_att_float(ncid, p_id, "missing_value", NC_FLOAT, 1, p_missing_value);
     nc_check_error(stat,__LINE__,__FILE__);
     p__FillValue[0] = 1e+35;
     stat = nc_put_att_float(ncid, p_id, "_FillValue", NC_FLOAT, 1, p__FillValue);
     nc_check_error(stat,__LINE__,__FILE__);
     p_scale_factor[0] = 100;
     stat = nc_put_att_float(ncid, p_id, "scale_factor", NC_FLOAT, 1, p_scale_factor);
     nc_check_error(stat,__LINE__,__FILE__);
     p_add_offset[0] = 0;
     stat = nc_put_att_float(ncid, p_id, "add_offset", NC_FLOAT, 1, p_add_offset);
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, p_id, "long_name", 22, "sea level air pressure");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, p_id, "standard_name", 25, "air_pressure_at_sea_level");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, p_id, "float_name", 1, "P");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, p_id, "content", 3, "TYX");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, p_id, "associate", 12, "time lat lon");
     }

   if(options[IBD]==1) {
     stat = nc_put_att_text(ncid, ibd_id, "units", 2, "cm");
     nc_check_error(stat,__LINE__,__FILE__);
     ibd_missing_value[0] = 1e+35;
     stat = nc_put_att_float(ncid, ibd_id, "missing_value", NC_FLOAT, 1, ibd_missing_value);
     nc_check_error(stat,__LINE__,__FILE__);
     ibd__FillValue[0] = 1e+35;
     stat = nc_put_att_float(ncid, ibd_id, "_FillValue", NC_FLOAT, 1, ibd__FillValue);
     nc_check_error(stat,__LINE__,__FILE__);
     ibd_scale_factor[0] = 1;
     stat = nc_put_att_float(ncid, ibd_id, "scale_factor", NC_FLOAT, 1, ibd_scale_factor);
     nc_check_error(stat,__LINE__,__FILE__);
     ibd_add_offset[0] = 0;
     stat = nc_put_att_float(ncid, ibd_id, "add_offset", NC_FLOAT, 1, ibd_add_offset);
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, ibd_id, "long_name", 28, "inverted varometer departure");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, ibd_id, "standard_name", 31, "inverted_barometer_departure");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, ibd_id, "float_name", 1, "ibd");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, ibd_id, "content", 3, "TYX");
     nc_check_error(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, ibd_id, "associate", 12, "time lat lon");
     }

   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 5, "CF1.0");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "grid_type", 7, "regular");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "file_name", 7, "toto.nc");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "production", 32, "Florent Lyard produced this file");
   nc_check_error(stat,__LINE__,__FILE__);

   /* leave define mode */
   stat = nc_enddef (ncid);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_close(ncid);
   nc_check_error(stat,__LINE__,__FILE__);

   stat=nc_open(filename,NC_WRITE,&ncid);
   stat=nc_put_var_double(ncid,lon_id,grid.x);
   stat=nc_put_var_double(ncid,lat_id,grid.y);
   stat = nc_close(ncid);
   return 0;
}

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

int poc_write_xy(const char *filename, grid_t grid, int id, short *buffer)

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
{
   int  ncid;

   size_t start[2];
   size_t count[2];
   int status=0;

   start[0]=0;
   start[1]=0;
   count[0]=grid.ny;
   count[1]=grid.nx;

   status=nc_open(filename,NC_WRITE,&ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status=nc_put_vara_short(ncid,id,start,count,buffer);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status = nc_close(ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
 error:
  return(status);
}

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

int poc_write_xy(const char *filename, grid_t grid, int id, int *buffer)

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
{
   int  ncid;

   size_t start[2];
   size_t count[2];
   int status=0;

   start[0]=0;
   start[1]=0;
   count[0]=grid.ny;
   count[1]=grid.nx;

   status=nc_open(filename,NC_WRITE,&ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status=nc_put_vara_int(ncid,id,start,count,buffer);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status = nc_close(ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
 error:
  return(status);
}

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

int poc_write_xy(const char *filename, grid_t grid, int id, float *buffer)

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
{
   int  ncid;

   size_t start[2];
   size_t count[2];
   int status=0;

   start[0]=0;
   start[1]=0;
   count[0]=grid.ny;
   count[1]=grid.nx;

   status=nc_open(filename,NC_WRITE,&ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status=nc_put_vara_float(ncid,id,start,count,buffer);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status = nc_close(ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
 error:
  return(status);
}

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

int poc_write_xy(const char *filename, grid_t grid, int id, double *buffer)

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
{
   int  ncid;

   size_t start[2];
   size_t count[2];
   int status=0;

   start[0]=0;
   start[1]=0;
   count[0]=grid.ny;
   count[1]=grid.nx;

   status=nc_open(filename,NC_WRITE,&ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status=nc_put_vara_double(ncid,id,start,count,buffer);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status = nc_close(ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
 error:
  return(status);
}

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  int poc_write_xyz(const char *filename, grid_t grid, int id, float *buffer)

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
{
   int  ncid;

   size_t start[4];
   size_t count[4];
   int status=0;

   start[0]=0;
   start[1]=0;
   start[2]=0;
   count[0]=grid.nz;
   count[1]=grid.ny;
   count[2]=grid.nx;

   status=nc_open(filename,NC_WRITE,&ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status=nc_put_vara_float(ncid,id,start,count,buffer);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status = nc_close(ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
 error:
  return(status);
}

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  int poc_write_xyz(const char *filename, grid_t grid, int id, double *buffer)

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
{
   int  ncid;

   size_t start[4];
   size_t count[4];
   int status=0;

   start[0]=0;
   start[1]=0;
   start[2]=0;
   count[0]=grid.nz;
   count[1]=grid.ny;
   count[2]=grid.nx;

   status=nc_open(filename,NC_WRITE,&ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status=nc_put_vara_double(ncid,id,start,count,buffer);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status = nc_close(ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
 error:
  return(status);
}

// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
//
//   int poc_write_xyt(const char *filename, grid_t grid, int frame, int id, float *buffer)
//
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// {
//    int  ncid;
//
//    size_t start[3];
//    size_t count[3];
//    int status=0;
//
//    start[0]=frame;
//    start[1]=0;
//    start[2]=0;
//    count[0]=1;
//    count[1]=grid.ny;
//    count[2]=grid.nx;
//
//    status=nc_open(filename,NC_WRITE,&ncid);
//    nc_check_error(status,__LINE__,__FILE__);
//    if(status !=0) goto error;
//    status=nc_put_vara_float(ncid,id,start,count,buffer);
//    nc_check_error(status,__LINE__,__FILE__);
//    if(status !=0) goto error;
//    status = nc_close(ncid);
//    nc_check_error(status,__LINE__,__FILE__);
//    if(status !=0) goto error;
//  error:
//   return(status);
// }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  int poc_write_xyzt(const char *filename, grid_t grid, int frame, int id, float *buffer)

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
{
   int  ncid;

   size_t start[4];
   size_t count[4];
   int status=0;

   start[0]=frame;
   start[1]=0;
   start[2]=0;
   start[3]=0;
   count[0]=1;
   count[1]=grid.nz;
   count[2]=grid.ny;
   count[3]=grid.nx;

   status=nc_open(filename,NC_WRITE,&ncid);
   if(status !=0) goto error;
   
   status=nc_put_vara_float(ncid,id,start,count,buffer);
   if(status !=0) goto error;
   
   status = nc_close(ncid);
   if(status !=0) goto error;
   
 error:
   if(status!=NC_NOERR) {
     nc_check_error(status,__LINE__,__FILE__);
     }
   return(status);
}

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

int poc_write_xyzt(const char *filename, grid_t grid, int frame, int id, double *buffer)

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
{
   int  ncid;

   size_t start[4];
   size_t count[4];
   int status=0;

   start[0]=frame;
   start[1]=0;
   start[2]=0;
   start[3]=0;
   count[0]=1;
   count[1]=grid.nz;
   count[2]=grid.ny;
   count[3]=grid.nx;

   status=nc_open(filename,NC_WRITE,&ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status=nc_put_vara_double(ncid,id,start,count,buffer);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status = nc_close(ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
 error:
  return(status);
}


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <typename type> int poc_write_template(const char *filename, cdfvar_t var, size_t *start,size_t *count,const type *buffer)

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
{
  int m,ncid;     /* netCDF id */
  int status=0;
  cdfvar_t info;
  type scale,offset,mask,*tmp=0;
  size_t size;
  
  status=nc_open(filename,NC_WRITE,&ncid);
  nc_check_error(status,__LINE__,__FILE__);
  if(status !=0) goto error_0;

  cdf_varinfo(ncid, var.id, &info);
  poc_decode_mask(info,&scale,&offset,&mask);

  if((scale != 1.0) || (offset != 0.0)) {
    size=1;
    for(m = 0; m < var.ndim; m++)
      size*=count[m];
    tmp=new type[size];
    for(m = 0; m < size; m++)
      if(buffer[m] != mask) {
        tmp[m] = (buffer[m]-offset)/scale;
        }
      else {
        tmp[m] = mask;
        }
    }

  status=poc_put_vara(ncid,var.id,start,count,buffer);
  nc_check_error(status,__LINE__,__FILE__,"poc_put_vara error on %s",var.name);
  if(status !=0) goto error;
  status = nc_close(ncid);
  nc_check_error(status,__LINE__,__FILE__);
  if(status !=0) goto error;

  if(tmp!=0) delete[] tmp;
  return(status);

 error:
  status = nc_close(ncid);
  status=-1;
 error_0:
  return(status);
}

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxXXX*/

int poc_write(const char *filename, cdfvar_t var, size_t *start,size_t *count,char *buffer)
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
{
  int result;
  
  result=poc_write_template(filename, var, start, count, buffer);
  
  return result;
}

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxXXX*/

int poc_write(const char *filename, cdfvar_t var, size_t *start,size_t *count,short *buffer)
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
{
  int result;
  
  result=poc_write_template(filename, var, start, count, buffer);
  
  return result;
}

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxXXX*/

int poc_write(const char *filename, cdfvar_t var, size_t *start,size_t *count,float *buffer)
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
{
  int result;
  
  result=poc_write_template(filename, var, start, count, buffer);
  
  return result;
}

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxXXX*/

int poc_write(const char *filename, cdfvar_t var, size_t *start,size_t *count,double *buffer)
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
{
  int result;
  
  result=poc_write_template(filename, var, start, count, buffer);
  
  return result;
}

#if 0 //int poc_put_vara(int ncid, int varid, const size_t *start, const size_t *count,const fcomplex *out) NOT DECLARED
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxXXX*/

int poc_write(const char *filename, cdfvar_t var, size_t *start,size_t *count,fcomplex *cbuffer)
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
{
  int result;
  
  result=poc_write_template(filename, var, start, count, buffer);
  
  return result;
}
#endif

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

int write_ncfile2dframe(char *filename, grid_t grid, int frame, int id, float *buffer)

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
{
   int  ncid;

   size_t start[3];
   size_t count[3];
   int status=0;

   start[0]=frame;
   start[1]=0;
   start[2]=0;
   count[0]=1;
   count[1]=grid.ny;
   count[2]=grid.nx;

   status=nc_open(filename,NC_WRITE,&ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status=nc_put_vara_float(ncid,id,start,count,buffer);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status = nc_close(ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
 error:
  return(status);
}

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

int write_ncfile2d(char *filename, grid_t grid, int ndim,int *tstart, int id, float *buffer)

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
{
   int  ncid;

   size_t start[2+ndim];
   size_t count[2+ndim];
   int status=0;
   int bcl;

   for(bcl=0;bcl < ndim ;bcl++) {
     start[bcl] = tstart[bcl];
     count[bcl] = 1;
     }
   start[bcl] = 0;
   count[bcl]=grid.ny;
   bcl++;
   start[bcl] = 0;
   count[bcl]=grid.nx;

   status=nc_open(filename,NC_WRITE,&ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status=nc_put_vara_float(ncid,id,start,count,buffer);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status = nc_close(ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
 error:
  return(status);
}

/*------------------------------------------------------------------------

  Nom         :  write_ncfile2d

  --------------------------------------------------------------------------*/
int write_ncfile2d_obsolete(char *filename, grid_t grid, int frame, int id, float *buffer)
{
   int  ncid;			/* netCDF id */

   size_t start[3];
   size_t count[3];
   int status=0;

   start[0]=frame;
   start[1]=0;
   start[2]=0;
   count[0]=1;
   count[1]=grid.ny;
   count[2]=grid.nx;

   status=nc_open(filename,NC_WRITE,&ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status=nc_put_vara_float(ncid,id,start,count,buffer);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status = nc_close(ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
 error:
  return(status);
}

/*------------------------------------------------------------------------

  Nom         :  write_ncfile3d

  --------------------------------------------------------------------------*/
int write_ncfile3d(char *filename, grid_t grid, int ndim,int *tstart, int id, float *buffer)
{
   int  ncid;			/* netCDF id */

   size_t start[3+ndim];
   size_t count[3+ndim];
   int status=0;
   int bcl;

   for(bcl=0;bcl < ndim ;bcl++) {
     start[bcl] = tstart[bcl];
     count[bcl] = 1;
     }
   start[bcl] = 0;
   count[bcl]=grid.nz;
   bcl ++;
   start[bcl] = 0;
   count[bcl]=grid.ny;
   bcl++;
   start[bcl] = 0;
   count[bcl]=grid.nx;

   status=nc_open(filename,NC_WRITE,&ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status=nc_put_vara_float(ncid,id,start,count,buffer);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status = nc_close(ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
 error:
  return(status);
}

/*------------------------------------------------------------------------

  Nom         :  write_ncfile3d

  --------------------------------------------------------------------------*/
int write_ncfile3d_obsolete(char *filename, grid_t grid, int frame, int id, float *buffer)
{
   int  ncid;			/* netCDF id */

   size_t start[4];
   size_t count[4];
   int status=0;

   start[0]=frame;
   start[1]=0;
   start[2]=0;
   start[3]=0;
   count[0]=1;
   count[1]=grid.nz;
   count[2]=grid.ny;
   count[3]=grid.nx;

   status=nc_open(filename,NC_WRITE,&ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status=nc_put_vara_float(ncid,id,start,count,buffer);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
   status = nc_close(ncid);
   nc_check_error(status,__LINE__,__FILE__);
   if(status !=0) goto error;
 error:
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int save_SG(const char *output,grid_t grid,float *buffer,float mask,const char *varname,const char *unit,const char *standard,int create)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,n,status;
  fcomplex z;
  cdfvar_t variable;
  pocgrd_t ncgrid;
  float *rbufx=NULL,*rbufy=NULL,*topo=NULL, rmask;

  
  printf("#################################################################\n");
  printf("write %s in output file : %s\n",varname, output);
  
  switch(create) {
    case -1:
      status=poc_sphericalgrid_xy(output,&ncgrid);
      if(status==0) break;
    case 1:
      status= poc_createfile(output);
      grid.nz=1;
      grid.z=NULL;
      status=poc_sphericalgrid_xy(output,"",grid,&ncgrid);
      break;
    case 0:
      status=poc_sphericalgrid_xy(output,&ncgrid);
      break;
    }

  if(status!=0) return(-1);
  
  poc_standardvariable_xy(&variable,varname,mask,unit,1., 0.,standard,standard,standard,ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  grid, variable.id,buffer);
  variable.destroy();
  
  return(0);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int save_SG(const char *output,grid_t grid,fcomplex *tide,fcomplex cmask,const char **varname,const char **unit,const char **standard,int create)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,n,status;
  fcomplex z;
  cdfvar_t variable;
  pocgrd_t ncgrid;
  float *rbufx=NULL,*rbufy=NULL,*topo=NULL, rmask;

  printf("#################################################################\n");
  printf("convert complex to amplitude, phase lag, units scale=%f\n",1.0);
  rmask=cmask.real();
  
  exitIfNull(
    rbufx=new float[grid.nx*grid.ny]
    );
  exitIfNull(
    rbufy=new float[grid.nx*grid.ny]
    );
    
  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i++) {
      n=i+grid.nx*j;
      z=tide[n];
      if(z!=cmask){
        rbufx[n]=abs(z);
        rbufy[n]=-arg(z)*r2d;
        }
      else {
        rbufx[n]=rmask;
        rbufy[n]=rmask;
        }
      }
    }
  
  printf("#################################################################\n");
  printf("write (%s,%s) in output file : %s\n",varname[0],varname[1], output);
  if(create==1) {
    status= poc_createfile(output);
    grid.nz=1;
    grid.z=NULL;
    status=poc_sphericalgrid_xy(output,"",grid,&ncgrid);
    }
  else {
    status=poc_sphericalgrid_xy(output,&ncgrid);
    }
    
  poc_standardvariable_xy(&variable,varname[0],rmask,unit[0],1., 0.,standard[0],standard[0],standard[0],ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  grid, variable.id,rbufx);
  variable.destroy();

  poc_standardvariable_xy(&variable,varname[1],rmask,unit[1],1., 0.,standard[1],standard[1],standard[1],ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  grid, variable.id,rbufy);
  variable.destroy();

//   if(topofile!=NULL) {
//     poc_standardvariable_xy(&variable, "bathymetry",1e+11,"m",1., 0.,"bathymetry","bathymetry","bathymetry",ncgrid);
//     status=create_ncvariable(output, &variable);
//     status=poc_write_xy(output,  topogrid, variable.id,topo);
//     variable.destroy();
//     }
  delete[] rbufx;
  delete[] rbufy;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int save_SGXY_template(const char *output, grid_t grid, T *buffer, T mask, const char *varname,const char *unit, const char *standard, int create, int write_grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  fcomplex z;
  cdfvar_t variable;
  pocgrd_t ncgrid;

//   printf("#################################################################\n");
  printf("#----------------------------------------------------------------\n");
  printf("writing %s to %s\n",varname,output);
  if(create==1) {
    status= poc_createfile(output);
    }
  if(write_grid==1) {
    grid.nz=1;
    grid.z=NULL;
    if(grid.modeH==0) status=map_completegridaxis(&grid,1);
    status=poc_sphericalgrid_xy(output,"",grid,&ncgrid);
    }
  else {
    status=poc_sphericalgrid_xy(output,"",grid,&ncgrid);
    }

  poc_standardvariable_xy(&variable,varname,mask,unit,1., 0.,standard,standard,standard,ncgrid);
  ncgrid.destroy();
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  grid, variable.id,buffer);
  variable.destroy();
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_SGXY(const char *output, const grid_t & grid, int *buffer, int mask, const char *varname,const char *unit, const char *standard, int create, int write_grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=save_SGXY_template(output,grid,buffer,mask,varname,unit,standard,create,write_grid);
  
  return status;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_SGXY(const char *output, const grid_t & grid, float *buffer, float mask, const char *varname,const char *unit, const char *standard, int create, int write_grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=save_SGXY_template(output,grid,buffer,mask,varname,unit,standard,create,write_grid);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_SGXY(const char *output, const grid_t & grid, double *buffer, double mask, const char *varname,const char *unit, const char *standard, int create, int write_grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=save_SGXY_template(output,grid,buffer,mask,varname,unit,standard,create,write_grid);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int read_SGatlas(const char *filename,const char **varnames,grid_t & grid,complex<float> * tide,complex<float> *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, n, nerrors, status;
  int kk,ll,mm;
  const int nnext=24;
  float *buf[2],spec[2];
  fcomplex zz;
  double x, y;
  size_t count[3], start[3];
  cdfgbl_t global;
  int id;
  variable_t varinfo;

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = 1;
  count[2] = 1;

  status= cdf_globalinfo(filename,&global,1);
  if(status!=NC_NOERR)return status;
  id=cdf_identify(global,varnames[0]);
  if(id<0){
    status=NC_ENOTVAR;/* Variable not found */
    nc_check_error(status,__LINE__,__FILE__,"error with %s in %s",varnames[0],filename);
    return status;
    }

  status=cdf_loadvargrid_2d (filename,id, &grid);
  if(status !=0) goto error;

  exitIfNull(buf[0] =new float   [grid.nx*grid.ny]);
  exitIfNull(buf[1] =new float   [grid.nx*grid.ny]);
  exitIfNull(tide   =new fcomplex[grid.nx*grid.ny]);

/*-----------------------------------------------------------------------
  load netcdf variable */
  id=cdf_identify(global,varnames[0]);
  status= cdf_loadvar_r1_2d (filename, id, 0, 0, grid, grid.nx, buf[0], &spec[0] ,&varinfo);
  varinfo.reset();

  id=cdf_identify(global,varnames[1]);
  if(id<0){
    status=NC_ENOTVAR;/* Variable not found */
    nc_check_error(status,__LINE__,__FILE__,"error with %s in %s",varnames[1],filename);
    return status;
    }
  status= cdf_loadvar_r1_2d (filename, id, 0, 0, grid, grid.nx, buf[1], &spec[1] ,&varinfo);
  varinfo.reset();
  global.destroy();

  *mask=fcomplex(9999.,9999.);
  for (j=0;j<grid.ny;j++)
    for (i=0;i<grid.nx;i++) {
      n=i+grid.nx*j;
      if((buf[0][n]!=spec[0]) && (buf[1][n]!=spec[1])) {
        tide[n]=polar<float>(buf[0][n],-buf[1][n]*d2r);
        }
      else {
        tide[n]=*mask;
        }
      }
  delete[]buf[0];
  delete[]buf[1];



  return (status);

error:
  map_printgrid(grid);
  printf("interpolation error: node %d lon=%lf lat=%lf \n", n, x, y);
  return (status);
}
