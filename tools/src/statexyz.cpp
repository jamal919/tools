
/*******************************************************************************

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

*******************************************************************************/


#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "tools-structures.h"

#include "fe.h"
#include "map.h"
#include "netcdf-proto.h"


/*----------------------------------------------------------------------------*/
/**
\sa create_ncfile_xy() and create_ncfile_xz()
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int create_ncfile_xyz(char *filename,size_t x_len, size_t y_len, size_t z_len, grid_t grid,float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *------------------------------------------------------------------------
  attribute name "axis" changed in "content" to comply with CF standard*/
{
   int  ncid; /* netCDF id */

   /* dimension ids */
   int x_dim;
   int y_dim;
   int z_dim;
   int time_counter_dim;

   /* dimension lengths */
/*
   size_t x_len = 100;
   size_t y_len = 100;
   size_t z_len = 100;
*/
   size_t time_counter_len = NC_UNLIMITED;

   /* variable ids */
   int lon_id;
   int lat_id;
   int depth_id;
   int level_depth_id;
   int time_id;
   int h_id;
   int UMEAN_id;
   int VMEAN_id;
   int T_id;
   int S_id;
   int U_id;
   int V_id;
   int W_id;
   int rho_id;

   /* rank (number of dimensions) for each variable */
#  define RANK_lon 2
#  define RANK_lat 2
#  define RANK_depth 2
#  define RANK_level_depth 3
#  define RANK_time 1
#  define RANK_h 3
#  define RANK_UMEAN 3
#  define RANK_VMEAN 3
#  define RANK_T 4
#  define RANK_S 4
#  define RANK_U 4
#  define RANK_V 4
#  define RANK_W 4
#  define RANK_rho 4

   /* variable shapes */
   int lon_dims[RANK_lon];
   int lat_dims[RANK_lat];
   int depth_dims[RANK_depth];
   int level_depth_dims[RANK_level_depth];
   int time_dims[RANK_time];
   int h_dims[RANK_h];
   int UMEAN_dims[RANK_UMEAN];
   int VMEAN_dims[RANK_VMEAN];
   int T_dims[RANK_T];
   int S_dims[RANK_S];
   int U_dims[RANK_U];
   int V_dims[RANK_V];
   int W_dims[RANK_W];
   int rho_dims[RANK_rho];

   /* attribute vectors */
   float depth_missing_value[1];
   float depth__FillValue[1];
   float depth_scale_factor[1];
   float depth_add_offset[1];
   float level_depth_missing_value[1];
   float level_depth__FillValue[1];
   float level_depth_scale_factor[1];
   float level_depth_add_offset[1];
   float h_missing_value[1];
   float h__FillValue[1];
   float h_scale_factor[1];
   float h_add_offset[1];
   float UMEAN_missing_value[1];
   float UMEAN__FillValue[1];
   float UMEAN_scale_factor[1];
   float UMEAN_add_offset[1];
   float VMEAN_missing_value[1];
   float VMEAN__FillVMEANalue[1];
   float VMEAN_scale_factor[1];
   float VMEAN_add_offset[1];
   float T_missing_value[1];
   float T__FillValue[1];
   float T_scale_factor[1];
   float T_add_offset[1];
   float S_missing_value[1];
   float S__FillValue[1];
   float S_scale_factor[1];
   float S_add_offset[1];
   float U_missing_value[1];
   float U__FillValue[1];
   float U_scale_factor[1];
   float U_add_offset[1];
   float V_missing_value[1];
   float V__FillValue[1];
   float V_scale_factor[1];
   float V_add_offset[1];
   float W_missing_value[1];
   float W__FillWalue[1];
   float W_scale_factor[1];
   float W_add_offset[1];
   float rho_missing_value[1];
   float rho__FillValue[1];
   float rho_scale_factor[1];
   float rho_add_offset[1];

   /* enter define mode */
   int stat = nc_create(filename, NC_CLOBBER, &ncid);
   nc_check_error(stat,__LINE__,__FILE__);

   /* define dimensions */
   stat = nc_def_dim(ncid, "x", x_len, &x_dim);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_def_dim(ncid, "y", y_len, &y_dim);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_def_dim(ncid, "z", z_len, &z_dim);
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

   depth_dims[0] = y_dim;
   depth_dims[1] = x_dim;
   stat = nc_def_var(ncid, "depth", NC_FLOAT, RANK_depth, depth_dims, &depth_id);
   nc_check_error(stat,__LINE__,__FILE__);

   level_depth_dims[0] = z_dim;
   level_depth_dims[1] = y_dim;
   level_depth_dims[2] = x_dim;
   stat = nc_def_var(ncid, "level_depth", NC_FLOAT, RANK_level_depth, level_depth_dims, &level_depth_id);
   nc_check_error(stat,__LINE__,__FILE__);

   time_dims[0] = time_counter_dim;
   stat = nc_def_var(ncid, "time", NC_DOUBLE, RANK_time, time_dims, &time_id);
   nc_check_error(stat,__LINE__,__FILE__);

   h_dims[0] = time_counter_dim;
   h_dims[1] = y_dim;
   h_dims[2] = x_dim;
   stat = nc_def_var(ncid, "h", NC_FLOAT, RANK_h, h_dims, &h_id);
   nc_check_error(stat,__LINE__,__FILE__);

   UMEAN_dims[0] = time_counter_dim;
   UMEAN_dims[1] = y_dim;
   UMEAN_dims[2] = x_dim;
   stat = nc_def_var(ncid, "UMEAN", NC_FLOAT, RANK_UMEAN, UMEAN_dims, &UMEAN_id);
   nc_check_error(stat,__LINE__,__FILE__);

   VMEAN_dims[0] = time_counter_dim;
   VMEAN_dims[1] = y_dim;
   VMEAN_dims[2] = x_dim;
   stat = nc_def_var(ncid, "VMEAN", NC_FLOAT, RANK_VMEAN, VMEAN_dims, &VMEAN_id);
   nc_check_error(stat,__LINE__,__FILE__);

   T_dims[0] = time_counter_dim;
   T_dims[1] = z_dim;
   T_dims[2] = y_dim;
   T_dims[3] = x_dim;
   stat = nc_def_var(ncid, "T", NC_FLOAT, RANK_T, T_dims, &T_id);
   nc_check_error(stat,__LINE__,__FILE__);

   S_dims[0] = time_counter_dim;
   S_dims[1] = z_dim;
   S_dims[2] = y_dim;
   S_dims[3] = x_dim;
   stat = nc_def_var(ncid, "S", NC_FLOAT, RANK_S, S_dims, &S_id);
   nc_check_error(stat,__LINE__,__FILE__);

   U_dims[0] = time_counter_dim;
   U_dims[1] = z_dim;
   U_dims[2] = y_dim;
   U_dims[3] = x_dim;
   stat = nc_def_var(ncid, "U", NC_FLOAT, RANK_U, U_dims, &U_id);
   nc_check_error(stat,__LINE__,__FILE__);

   V_dims[0] = time_counter_dim;
   V_dims[1] = z_dim;
   V_dims[2] = y_dim;
   V_dims[3] = x_dim;
   stat = nc_def_var(ncid, "V", NC_FLOAT, RANK_V, V_dims, &V_id);
   nc_check_error(stat,__LINE__,__FILE__);

   W_dims[0] = time_counter_dim;
   W_dims[1] = z_dim;
   W_dims[2] = y_dim;
   W_dims[3] = x_dim;
   stat = nc_def_var(ncid, "W", NC_FLOAT, RANK_W, W_dims, &W_id);
   nc_check_error(stat,__LINE__,__FILE__);

   rho_dims[0] = time_counter_dim;
   rho_dims[1] = z_dim;
   rho_dims[2] = y_dim;
   rho_dims[3] = x_dim;
   stat = nc_def_var(ncid, "rho", NC_FLOAT, RANK_rho, rho_dims, &rho_id);
   nc_check_error(stat,__LINE__,__FILE__);

   /* assign attributes */
   stat = nc_put_att_text(ncid, lon_id, "units", 12, "degrees_east");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lon_id, "long_name", 9, "longitude");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lon_id, "standard_name", 9, "longitude");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lon_id, "nav_model", 12, "Default grid");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lat_id, "units", 13, "degrees_north");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lat_id, "long_name", 8, "latitude");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lat_id, "standard_name", 8, "latitude");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lat_id, "nav_model", 12, "Default grid");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, depth_id, "units", 1, "m");
   nc_check_error(stat,__LINE__,__FILE__);
   depth_missing_value[0] = mask;
   stat = nc_put_att_float(ncid, depth_id, "missing_value", NC_FLOAT, 1, depth_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   depth__FillValue[0] = mask;
   stat = nc_put_att_float(ncid, depth_id, "_FillValue", NC_FLOAT, 1, depth__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   depth_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, depth_id, "scale_factor", NC_FLOAT, 1, depth_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   depth_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, depth_id, "add_offset", NC_FLOAT, 1, depth_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, depth_id, "long_name", 10, "bathymetry");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, depth_id, "standard_name", 10, "bathymetry");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, depth_id, "content", 3, "TYX");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, depth_id, "associate", 7, "lat lon");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, level_depth_id, "units", 1, "m");
   nc_check_error(stat,__LINE__,__FILE__);
   level_depth_missing_value[0] = mask;
   stat = nc_put_att_float(ncid, level_depth_id, "missing_value", NC_FLOAT, 1, level_depth_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   level_depth__FillValue[0] = mask;
   stat = nc_put_att_float(ncid, level_depth_id, "_FillValue", NC_FLOAT, 1, level_depth__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   level_depth_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, level_depth_id, "scale_factor", NC_FLOAT, 1, level_depth_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   level_depth_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, level_depth_id, "add_offset", NC_FLOAT, 1, level_depth_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, level_depth_id, "long_name", 7, "Z level");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, level_depth_id, "standard_name", 7, "Z_level");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, level_depth_id, "content", 3, "ZYX");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, level_depth_id, "associate", 19, "level_depth lat lon");
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
   stat = nc_put_att_text(ncid, h_id, "units", 1, "m");
   nc_check_error(stat,__LINE__,__FILE__);
   h_missing_value[0] = mask;
   stat = nc_put_att_float(ncid, h_id, "missing_value", NC_FLOAT, 1, h_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   h__FillValue[0] = mask;
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
   stat = nc_put_att_text(ncid, h_id, "content", 3, "TYX");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, h_id, "associate", 12, "time lat lon");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, UMEAN_id, "units", 5, "m s-1");
   nc_check_error(stat,__LINE__,__FILE__);
   UMEAN_missing_value[0] = mask;
   stat = nc_put_att_float(ncid, UMEAN_id, "missing_value", NC_FLOAT, 1, UMEAN_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   UMEAN__FillValue[0] = mask;
   stat = nc_put_att_float(ncid, UMEAN_id, "_FillValue", NC_FLOAT, 1, UMEAN__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   UMEAN_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, UMEAN_id, "scale_factor", NC_FLOAT, 1, UMEAN_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   UMEAN_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, UMEAN_id, "add_offset", NC_FLOAT, 1, UMEAN_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, UMEAN_id, "long_name", 37, "ocean barotropic velocity E-component");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, UMEAN_id, "standard_name", 27, "eastward_sea_water_velocity");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, UMEAN_id, "content", 3, "TYX");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, UMEAN_id, "associate", 12, "time lat lon");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, VMEAN_id, "units", 5, "m s-1");
   nc_check_error(stat,__LINE__,__FILE__);
   VMEAN_missing_value[0] = mask;
   stat = nc_put_att_float(ncid, VMEAN_id, "missing_value", NC_FLOAT, 1, VMEAN_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   VMEAN__FillVMEANalue[0] = mask;
   stat = nc_put_att_float(ncid, VMEAN_id, "_FillVMEANalue", NC_FLOAT, 1, VMEAN__FillVMEANalue);
   nc_check_error(stat,__LINE__,__FILE__);
   VMEAN_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, VMEAN_id, "scale_factor", NC_FLOAT, 1, VMEAN_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   VMEAN_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, VMEAN_id, "add_offset", NC_FLOAT, 1, VMEAN_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, VMEAN_id, "long_name", 37, "ocean barotropic velocity N-component");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, VMEAN_id, "standard_name", 28, "northward_sea_water_velocity");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, VMEAN_id, "content", 3, "TYX");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, VMEAN_id, "associate", 12, "time lat lon");
   nc_check_error(stat,__LINE__,__FILE__);
   T_missing_value[0] = mask;
   stat = nc_put_att_float(ncid, T_id, "missing_value", NC_FLOAT, 1, T_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   T__FillValue[0] = mask;
   stat = nc_put_att_float(ncid, T_id, "_FillValue", NC_FLOAT, 1, T__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   T_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, T_id, "scale_factor", NC_FLOAT, 1, T_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   T_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, T_id, "add_offset", NC_FLOAT, 1, T_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, T_id, "standard_name", 21, "sea_water_temperature");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, T_id, "long_name", 21, "sea water temperature");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, T_id, "units", 2, "\260C");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, T_id, "content", 4, "TZYX");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, T_id, "associate", 24, "time level_depth lat lon");
   nc_check_error(stat,__LINE__,__FILE__);
   S_missing_value[0] = mask;
   stat = nc_put_att_float(ncid, S_id, "missing_value", NC_FLOAT, 1, S_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   S__FillValue[0] = mask;
   stat = nc_put_att_float(ncid, S_id, "_FillValue", NC_FLOAT, 1, S__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   S_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, S_id, "scale_factor", NC_FLOAT, 1, S_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   S_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, S_id, "add_offset", NC_FLOAT, 1, S_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, S_id, "standard_name", 18, "sea_water_salinity");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, S_id, "long_name", 18, "sea water salinity");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, S_id, "units", 6, "PSU   ");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, S_id, "content", 4, "TZYX");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, S_id, "associate", 24, "time level_depth lat lon");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, U_id, "units", 5, "m s-1");
   nc_check_error(stat,__LINE__,__FILE__);
   U_missing_value[0] = mask;
   stat = nc_put_att_float(ncid, U_id, "missing_value", NC_FLOAT, 1, U_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   U__FillValue[0] = mask;
   stat = nc_put_att_float(ncid, U_id, "_FillValue", NC_FLOAT, 1, U__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   U_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, U_id, "scale_factor", NC_FLOAT, 1, U_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   U_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, U_id, "add_offset", NC_FLOAT, 1, U_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, U_id, "long_name", 26, "ocean velocity E-component");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, U_id, "standard_name", 27, "eastward_sea_water_velocity");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, U_id, "content", 4, "TZYX");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, U_id, "associate", 24, "time level_depth lat lon");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, V_id, "units", 5, "m s-1");
   nc_check_error(stat,__LINE__,__FILE__);
   V_missing_value[0] = mask;
   stat = nc_put_att_float(ncid, V_id, "missing_value", NC_FLOAT, 1, V_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   V__FillValue[0] = mask;
   stat = nc_put_att_float(ncid, V_id, "_FillValue", NC_FLOAT, 1, V__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   V_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, V_id, "scale_factor", NC_FLOAT, 1, V_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   V_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, V_id, "add_offset", NC_FLOAT, 1, V_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, V_id, "long_name", 26, "ocean velocity N-component");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, V_id, "standard_name", 28, "northward_sea_water_velocity");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, V_id, "content", 4, "TZYX");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, V_id, "associate", 24, "time level_depth lat lon");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, W_id, "units", 5, "m s-1");
   nc_check_error(stat,__LINE__,__FILE__);
   W_missing_value[0] = mask;
   stat = nc_put_att_float(ncid, W_id, "missing_value", NC_FLOAT, 1, W_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   W__FillWalue[0] = mask;
   stat = nc_put_att_float(ncid, W_id, "_FillWalue", NC_FLOAT, 1, W__FillWalue);
   nc_check_error(stat,__LINE__,__FILE__);
   W_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, W_id, "scale_factor", NC_FLOAT, 1, W_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   W_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, W_id, "add_offset", NC_FLOAT, 1, W_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, W_id, "long_name", 26, "ocean velocity Z-component");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, W_id, "standard_name", 25, "upward_sea_water_velocity");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, W_id, "content", 4, "TZYX");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, W_id, "associate", 24, "time level_depth lat lon");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, rho_id, "units", 3, "PSI");
   nc_check_error(stat,__LINE__,__FILE__);
   rho_missing_value[0] = mask;
   stat = nc_put_att_float(ncid, rho_id, "missing_value", NC_FLOAT, 1, rho_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   rho__FillValue[0] = mask;
   stat = nc_put_att_float(ncid, rho_id, "_FillValue", NC_FLOAT, 1, rho__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   rho_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, rho_id, "scale_factor", NC_FLOAT, 1, rho_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   rho_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, rho_id, "add_offset", NC_FLOAT, 1, rho_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, rho_id, "long_name", 13, "water density");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, rho_id, "standard_name", 13, "water_density");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, rho_id, "content", 4, "TZYX");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, rho_id, "associate", 24, "time level_depth lat lon");
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
   return 0;


}

/*-----------------------------------------------------------------------------*/


int nc_writeXYZ(int ncid, grid_t grid, int var, float *z)


/*-----------------------------------------------------------------------------*/

{

/*-----------------------------------------------------------------------------*/

  int  status;
  size_t start[3];
  size_t count[3];


  start[0]=0;
  start[1]=0;
  start[2]=0;

  count[0]=grid.nz;
  count[1]=grid.ny;
  count[2]=grid.nx;
  
  status=nc_put_vara_float(ncid,var,start,count,z);
  if (status != NC_NOERR) goto error;
    
  return (0);

 error:
  nc_check_error(status,__LINE__,__FILE__);
  return(status);
}

/*-----------------------------------------------------------------------------*/


int nc_writeXYZT(int ncid, grid_t grid, int var, int frame, float *z)


/*-----------------------------------------------------------------------------*/

{

/*-----------------------------------------------------------------------------*/

  int  status;
  size_t start[4];
  size_t count[4];


  start[0]=frame;
  start[1]=0;
  start[2]=0;
  start[3]=0;

  count[0]=1;
  count[1]=grid.nz;
  count[2]=grid.ny;
  count[3]=grid.nx;
  
  status=nc_put_vara_float(ncid,var,start,count,z);
  if (status != NC_NOERR) goto error;
    
  return (0);

 error:
  nc_check_error(status,__LINE__,__FILE__);
  return(status);
}

/*-----------------------------------------------------------------------------*/


int nc_writeXYT(int ncid, grid_t grid, int var, int wave, float *z)


/*-----------------------------------------------------------------------------*/

{

/*-----------------------------------------------------------------------------*/

  int  status;
  size_t start[3];
  size_t count[3];


  start[0]=wave;
  start[1]=0;
  start[2]=0;

  count[0]=1;
  count[1]=grid.ny;
  count[2]=grid.nx;
  
  status=nc_put_vara_float(ncid,var,start,count,z);
  if (status != NC_NOERR) goto error;
    
  return (0);

 error:
  nc_check_error(status,__LINE__,__FILE__);
  return(status);
}

/*-----------------------------------------------------------------------------*/


int nc_writeXY(int ncid, grid_t grid, int var, float *z)


/*-----------------------------------------------------------------------------*/

{

/*-----------------------------------------------------------------------------*/

  int  status;
  size_t start[3];
  size_t count[3];


  start[0]=0;
  start[1]=0;

  count[0]=grid.ny;
  count[1]=grid.nx;
  
  status=nc_put_vara_float(ncid,var,start,count,z);
  if (status != NC_NOERR) goto error;
  return (0);

 error:
  nc_check_error(status,__LINE__,__FILE__);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t get_zonegrid_fatal(char *zone)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  grid_t grid;
  int zone_initialised=0;

/* NEA grid */
  if(strcmp(zone,"NEA")==0) {
  map_set2Dgrid(&grid,-20.0,+29.5,+15.0,+65.0,1./10.,1./10.);
  zone_initialised=1;
  }

/* NEA_hr grid haute resolution */
  if(strcmp(zone,"NEA_hr")==0) {
  map_set2Dgrid(&grid,-20.0,+29.5,+15.0,+65.0,1./50.,1./50.);
  zone_initialised=1;
  }

/* Manche_hr grid haute resolution */
  if(strcmp(zone,"Manche_hr")==0) {
  map_set2Dgrid(&grid,-6.0,+6.0,+48.0,+53.0,1./100.,1./100.);
  zone_initialised=1;
  }

 /* NEA_shom grid distribution shom  haute resolution */
  if(strcmp(zone,"NEA_shom")==0) {
  map_set2Dgrid(&grid,-20.0,+40.0,+15.0,+55.0,1./30.,1./30.);
  zone_initialised=1;
  }
      
/* Iroise grid */
  if(strcmp(zone,"iroise")==0) {
  map_set2Dgrid(&grid,-6.0,+47.25,-3.5,+49.0,1./60.,1./60.);
  zone_initialised=1;
  }
    
/* global grid */
  if(strcmp(zone,"global")==0) {
  map_set2Dgrid(&grid,0.0,-80.0,+360.0,+80.0,1./4.,1./4.);
  zone_initialised=1;
  }
  
/* global grid */
  if(strcmp(zone,"global-loren")==0) {
  map_set2Dgrid(&grid,0.0,-90.0,+360.0,+90.0,1.,1.);
  zone_initialised=1;
  }
  
    
/* mesdsea grid, normal resolution */
  if(strcmp(zone,"medsea")==0) {
  map_set2Dgrid(&grid,-10.0,27.5,40.0,47.5,1./10.,1./10.);
  zone_initialised=1;
  }

/* alobran grid, normal resolution */
  if(strcmp(zone,"alboran")==0) {
  map_set2Dgrid(&grid,-9.0,32,-2+1./300,38+1./300,1./60.,1./60.);
  zone_initialised=1;
  }
  
/* alobran grid, normal resolution */
  if(strcmp(zone,"strait-of-sicily")==0) {
  map_set2Dgrid(&grid,-9.0,32,-2+1./300,38+1./300,1./60.,1./60.);
  zone_initialised=1;
  }
 
  if(zone_initialised!=1) {
    STDOUT_BASE_LINE("no valid region specified (%s); abort...\n",zone);
    exit(-1);
    }
  grid.nz  = 1;

  return(grid);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mog3d_load(char *name, mesh_t mesh, int nlyrs, grid_t grid, int ncid, int *elts, float mask, int frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  
  int   k,m,n,nitems,var,status;
  FILE *in;
  float *rbufN,*rbufC,*rbuf2d,*rbuf3d;
  float  *buffer,*buffer4;
  double *buffer8,seconds;
  double time;
  int nndes,  nelts, offset;
  int date[10];
  grid_t grid2d;

  double *sigma;
  double *depth;
  float  *zlevels;

  in=fopen(name,"r");
  if(in == NULL) {
    printf("unable to open %s\n",name);
    return (-1);
    }

  nndes=mesh.nvtxs;
  nelts=mesh.ntriangles;

  grid2d= map_getgrid2d(grid);

  rbufN=(float *)malloc(nndes*sizeof(float));
  rbufC=(float *)malloc(nelts*sizeof(float));

  buffer4=(float *)malloc(nlyrs*nelts*sizeof(float));
  buffer8=(double *)malloc(nlyrs*nelts*sizeof(double));
  buffer=buffer4;

  rbuf2d=(float *)malloc(grid.nx*grid.ny*sizeof(float));
  rbuf3d=(float *)malloc(grid.nx*grid.ny*grid.nz*sizeof(float));

/*------------------------------------------------------------------------------
  special buffers */
  sigma=(double *)malloc(nlyrs*nelts*sizeof(double));
  //zlevels=(double *)malloc(nlyrs*nelts*sizeof(double)); Certainement un bug Thierry
  zlevels=(float *)malloc(nlyrs*nelts*sizeof(float));
  depth=(double *)malloc(nelts*sizeof(double));

/*------------------------------------------------------------------------------
  month,day,year,seconds */
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(&date,3,sizeof(int),in);
  nitems=fread(&seconds,1,sizeof(double),in);
  nitems=fread(&time,1,sizeof(double),in);
  nitems=fread(&offset,sizeof(int),1,in);

/*------------------------------------------------------------------------------
  nndes,nelts,nlyrs */
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(&m,sizeof(int),1,in);
  nitems=fread(&m,sizeof(int),1,in);
  nitems=fread(&nlyrs,sizeof(int),1,in);
  nitems=fread(&offset,sizeof(int),1,in);

  nlyrs=11;
/*------------------------------------------------------------------------------
  U,V,W centroid */
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nlyrs*nelts,in);
  if(nitems!=nlyrs*nelts) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nlyrs*nelts,in);
  if(nitems!=nlyrs*nelts) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nlyrs*nelts,in);
  if(nitems!=nlyrs*nelts) goto error;
  if(offset!=nlyrs*nelts*sizeof(double)) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

/*------------------------------------------------------------------------------
  turbulence centroid */
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nlyrs*nelts,in);
  if(nitems!=nlyrs*nelts) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nlyrs*nelts,in);
  if(nitems!=nlyrs*nelts) goto error;
  nitems=fread(&offset,sizeof(int),1,in);
 
/*------------------------------------------------------------------------------
  S noeuds */
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nlyrs*nndes,in);
  if(nitems!=nlyrs*nndes) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  for(k=0;k<nlyrs-1;k++) {
    for(n=0;n<nndes;n++) rbufN[n]=buffer8[n*nlyrs+k];
    status=fe_map(mesh,rbufN,grid2d,elts,rbuf2d,mask);
    if(status!=0) goto error;
    for(n=0;n<grid.nx*grid.ny;n++) rbuf3d[n+k*grid.nx*grid.ny]=rbuf2d[n];
    }
  for(n=0;n<grid.nx*grid.ny;n++) rbuf3d[n+k*grid.nx*grid.ny]=mask;
  status=nc_inq_varid(ncid,"S",&var);
  if(status!=0) goto error;
  status=nc_writeXYZT(ncid, grid, var, frame, rbuf3d);
  if(status!=0) goto error;

/*------------------------------------------------------------------------------
  T noeuds */
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nlyrs*nndes,in);
  if(nitems!=nlyrs*nndes) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  for(k=0;k<nlyrs-1;k++) {
    for(n=0;n<nndes;n++) rbufN[n]=buffer8[n*nlyrs+k];
    status=fe_map(mesh,rbufN,grid2d,elts,rbuf2d,mask);
    if(status!=0) goto error;
    for(n=0;n<grid.nx*grid.ny;n++) rbuf3d[n+k*grid.nx*grid.ny]=rbuf2d[n];
    }
  for(n=0;n<grid.nx*grid.ny;n++) rbuf3d[n+k*grid.nx*grid.ny]=mask;
  status=nc_inq_varid(ncid,"T",&var);
  if(status!=0) goto error;
  status=nc_writeXYZT(ncid, grid, var, frame, rbuf3d);
  if(status!=0) goto error;

/*------------------------------------------------------------------------------
  Rho noeuds */
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nlyrs*nndes,in);
  if(nitems!=nlyrs*nndes) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  for(k=0;k<nlyrs-1;k++) {
    for(n=0;n<nndes;n++) rbufN[n]=buffer8[n*nlyrs+k];
    status=fe_map(mesh,rbufN,grid2d,elts,rbuf2d,mask);
    if(status!=0) goto error;
    for(n=0;n<grid.nx*grid.ny;n++) rbuf3d[n+k*grid.nx*grid.ny]=rbuf2d[n];
    }
  status=nc_inq_varid(ncid,"rho",&var);
  if(status!=0) goto error;
  status=nc_writeXYZT(ncid, grid, var, frame, rbuf3d);
  if(status!=0) goto error;

/*------------------------------------------------------------------------------
  T,S,Rho background noeuds */
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nlyrs*nndes,in);
  if(nitems!=nlyrs*nndes) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nlyrs*nndes,in);
  if(nitems!=nlyrs*nndes) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nlyrs*nndes,in);
  if(nitems!=nlyrs*nndes) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

/*------------------------------------------------------------------------------
  diffusivit� centroid */
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nlyrs*nelts,in);
  if(nitems!=nlyrs*nelts) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nlyrs*nelts,in);
  if(nitems!=nlyrs*nelts) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nlyrs*nelts,in);
  if(nitems!=nlyrs*nelts) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

/*------------------------------------------------------------------------------
  U,V barotrope centroid */
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nelts,in);
  if(nitems!=nelts) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nelts,in);
  if(nitems!=nelts) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

/*------------------------------------------------------------------------------
  surface elevation t, t-1 centroid */
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nelts,in);
  if(nitems!=nelts) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nelts,in);
  if(nitems!=nelts) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

/*------------------------------------------------------------------------------
  profondeur moyenne,instantann�e t, t-1 centroid */
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nelts,in);
  if(nitems!=nelts) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nelts,in);
  if(nitems!=nelts) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nelts,in);
  if(nitems!=nelts) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

/*------------------------------------------------------------------------------
  vitesse vertical (sigma) aux noeuds*/
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nlyrs*nndes,in);
  if(nitems!=nlyrs*nndes) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

/*------------------------------------------------------------------------------
  surface elevation t, t-1 aux noeuds */
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nndes,in);
  if(nitems!=nndes) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  for(n=0;n<nndes;n++) rbufN[n]=buffer8[n];
  status=fe_map(mesh,rbufN,grid2d,elts,rbuf2d,mask);
  if(status!=0) goto error;
  status=nc_inq_varid(ncid,"h",&var);
  if(status!=0) goto error;
  status=nc_writeXYT(ncid, grid2d, var, frame, rbuf2d);
  if(status!=0) goto error;

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nndes,in);
  if(nitems!=nndes) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

/*------------------------------------------------------------------------------
  profondeur moyenne,instantann�e t, t-1 aux noeuds */
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(depth,sizeof(double),nndes,in);
  if(nitems!=nndes) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  for(n=0;n<nndes;n++) rbufN[n]=depth[n];
  status=fe_map(mesh,rbufN,grid2d,elts,rbuf2d,mask);
  if(status!=0) goto error;
  status=nc_inq_varid(ncid,"depth",&var);
  if(status!=0) goto error;
  status=nc_writeXY(ncid, grid2d, var, rbuf2d);
  if(status!=0) goto error;

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nndes,in);
  if(nitems!=nndes) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(double),nndes,in);
  if(nitems!=nndes) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

/*------------------------------------------------------------------------------
  drying t, t-1 aux noeuds */
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(int),nndes,in);
  if(nitems!=nndes) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(int),nndes,in);
  if(nitems!=nndes) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

/*------------------------------------------------------------------------------
  drying t, t-1 centroids */
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(rbufC,sizeof(int),nelts,in);
  if(nitems!=nelts) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(rbufC,sizeof(int),nelts,in);
  if(nitems!=nelts) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(rbufC,sizeof(int),nelts,in);
  if(nitems!=nelts) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

/*------------------------------------------------------------------------------
  sigma levels [0:-1] noeuds extremites et centre */
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(sigma,sizeof(double),nndes*nlyrs,in);
  if(nitems!=nndes*nlyrs) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

/*------------------------------------------------------------------------------
  compute zlevels from sigma and depth */
/*
  for(n=0;n<nndes;n++)
    for(l=0;l<nlyrs;l++)
      zlevels[n*nlyrs+l]=sigma[n*nlyrs+l]*depth[n];
*/
  for(k=0;k<nlyrs;k++) {
    for(n=0;n<nndes;n++) rbufN[n]=sigma[n*nlyrs+k]*depth[n];
    status=fe_map(mesh,rbufN,grid2d,elts,rbuf2d,mask);
    if(status!=0) goto error;
    for(n=0;n<grid.nx*grid.ny;n++) rbuf3d[n+k*grid.nx*grid.ny]=rbuf2d[n];
    }
  status=nc_inq_varid(ncid,"level_depth",&var);
  if(status!=0) goto error;
  status=nc_writeXYZ(ncid, grid, var, rbuf3d);
  if(status!=0) goto error;

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(int),nndes*nlyrs,in);
  if(nitems!=nndes*nlyrs) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

/*------------------------------------------------------------------------------
  sigma levels [0:-1] centroids  extremites et centre */
  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(int),nelts*nlyrs,in);
  if(nitems!=nelts*nlyrs) goto error;
  nitems=fread(&offset,sizeof(int),1,in);

  nitems=fread(&offset,sizeof(int),1,in);
  nitems=fread(buffer8,sizeof(int),nelts*nlyrs,in);
  if(nitems!=nelts*nlyrs) goto error;
  nitems=fread(&offset,sizeof(int),1,in);


 //finished:
  free(rbufN);
  free(rbufC);
  free(buffer4);
  free(buffer8);
  free(rbuf2d);
  free(rbuf3d);
  free(sigma);

  fclose(in);
  return (0);

 error:
  free(rbufN);
  free(rbufC);
  free(buffer4);
  free(buffer8);
  free(rbuf2d);
  free(rbuf3d);

  fclose(in);
  return (-1);

  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  rmask;
  float  spec=-9999;
  int nndes,status,*elts;
  int n;
  char *keyword,*zone=NULL,*sectionfile=NULL,*rootname=NULL;
  char *meshfile=NULL,*output=NULL,*path=NULL,*nodefile=NULL;
  grid_t grid;
  grid_t   grid2d;
  mesh_t mesh;
  int ncid,iterative=0,count;
  int nlyrs=43;

  fct_echo(argc,argv);
 
  n=1;
  while (n < argc)
    {
    keyword=strdup(argv[n]);
    switch (keyword[0])
      {
      case '-':
        switch (keyword[1])
        {
        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'n' :
          nodefile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          sscanf(argv[n+1],"%f",&spec);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          path= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'r' :
          rootname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
          sectionfile= strdup(argv[n]);
          n++;
          break;
        break;
      }
      free(keyword);
    }

  rmask=spec;
 
 if(path ==NULL) path=strdup(".");

 if(sectionfile==NULL) {
   iterative=1;
   sectionfile=(char *)malloc(1024);
   }
 
 if(meshfile != NULL) {
    status=fe_readmesh_TGL(meshfile,nodefile,&mesh);
    if(status !=0) goto error;
/*
    status=fe_list(&mesh);
    if(status !=0) goto error;
*/
    }
 else
   {
   printf("no mesh file specified; abort...\n");
   goto error;
   }

  mesh.type=1;
  nndes=mesh.nvtxs;
  //fe_initaffine(&mesh);

  if(zone !=NULL) grid=get_zonegrid_fatal(zone);
  else            grid=get_grid_d(mesh);
 
  status=map3d_completegridaxis(&grid);

  grid2d= map_getgrid2d(grid);

  elts=fe_scan_elements(mesh,grid2d,mesh.type);
  if(elts==NULL) goto error;

  output=(char *)malloc(1024);

/*   sprintf(output,"/home/ocean/lyard/%s.nc",rootname); */
  sprintf(output,"%s.nc",rootname);

  rmask=-9999.;

  nlyrs=11;
  grid.nz=nlyrs;

  status=create_ncfile_xyz(output, grid.nx, grid.ny, grid.nz, grid,rmask);
  if(status !=0) goto error;

  status=nc_open(output,NC_WRITE,&ncid);
  if(status !=0) goto error;

/**----------------------------------------------------------------------------
  obsolete call, will be suppressed */
  status=nc_writeaxis_3d(ncid,grid);
  if(status !=0) goto error;

  count=0;
  do
    {
    count++;
    sprintf(sectionfile,"%s.%d",rootname,count);
/*     sprintf(sectionfile,"%s",rootname); */
    printf("reading %s\n",sectionfile);
    status=mog3d_load(sectionfile, mesh, nlyrs, grid,  ncid,  elts,  rmask, count-1);
    if(status!=0) {
      printf("cannot read %s\n",sectionfile);
      goto finished;
      }
/*     else goto finished; */
    } while(status==0);
  
 finished:
  status = nc_close(ncid);

  printf("output nc file: %s\n",output);
  free(elts);
  TRAP_ERR_EXIT(0,"exiting\n");
error:
 STDOUT_BASE_LINE("error detected, quit ... \n");
  exit(-1);
}
