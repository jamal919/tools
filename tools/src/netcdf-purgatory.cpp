
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

\brief variable and grid loading poc-netcdf definitions

Old note : variables routines
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "tools-structures.h"

#include "netcdf-proto.h"
#include "poc-netcdf.def"
#include "map.h"
#include "poc-time.h"
#include "zapper.h"


/** Purgatoire */






int create_header(char *output, grid_t grid, spectrum_t spectrum)
/* *------------------------------------------------------------------------
  attribute name "axis" changed in "content" to comply with CF standard*/
 {			

   int  ncid;			/* netCDF id */

   /* dimension ids */
   int X_dim;
   int Y_dim;
   int W_dim;
   int L_dim;

   /* dimension lengths */
   size_t X_len = grid.nx;
   size_t Y_len = grid.ny;
   size_t W_len = spectrum.n;
   size_t L_len = 10;

   /* variable ids */
   int lon_id;
   int lat_id;
   int spectrum_id;
   int Ha_id;
   int Hg_id;

   /* variable shapes */
   int lon_dims[RANK_lon];
   int lat_dims[RANK_lat];
   int spectrum_dims[RANK_spectrum];
   int Ha_dims[RANK_Ha];
   int Hg_dims[RANK_Hg];

   /* attribute vectors */
   double lon_valid_min[1];
   double lon_valid_max[1];
   double lat_valid_min[1];
   double lat_valid_max[1];
   float Ha__FillValue[1];
   float Ha_missing_value[1];
   double Ha_scale_factor[1];
   double Ha_add_offset[1];
   float Hg__FillValue[1];
   float Hg_missing_value[1];
   double Hg_scale_factor[1];
   double Hg_add_offset[1];

   /* enter define mode */
   int stat = nc_create(output, NC_CLOBBER, &ncid);
   nc_check_error(stat,__LINE__,__FILE__);

   /* define dimensions */
   stat = nc_def_dim(ncid, "X", X_len, &X_dim);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_def_dim(ncid, "Y", Y_len, &Y_dim);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_def_dim(ncid, "W", W_len, &W_dim);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_def_dim(ncid, "L", L_len, &L_dim);
   nc_check_error(stat,__LINE__,__FILE__);

   /* define variables */

   lon_dims[0] = Y_dim;
   lon_dims[1] = X_dim;
   stat = nc_def_var(ncid, "lon", NC_DOUBLE, RANK_lon, lon_dims, &lon_id);
   nc_check_error(stat,__LINE__,__FILE__);

   lat_dims[0] = Y_dim;
   lat_dims[1] = X_dim;
   stat = nc_def_var(ncid, "lat", NC_DOUBLE, RANK_lat, lat_dims, &lat_id);
   nc_check_error(stat,__LINE__,__FILE__);

   spectrum_dims[0] = W_dim;
   spectrum_dims[1] = L_dim;
   stat = nc_def_var(ncid, "spectrum", NC_CHAR, RANK_spectrum, spectrum_dims, &spectrum_id);
   nc_check_error(stat,__LINE__,__FILE__);

   Ha_dims[0] = W_dim;
   Ha_dims[1] = Y_dim;
   Ha_dims[2] = X_dim;
   stat = nc_def_var(ncid, "Ha", NC_FLOAT, RANK_Ha, Ha_dims, &Ha_id);
   nc_check_error(stat,__LINE__,__FILE__);

   Hg_dims[0] = W_dim;
   Hg_dims[1] = Y_dim;
   Hg_dims[2] = X_dim;
   stat = nc_def_var(ncid, "Hg", NC_FLOAT, RANK_Hg, Hg_dims, &Hg_id);
   nc_check_error(stat,__LINE__,__FILE__);

   /* assign attributes */
   stat = nc_put_att_text(ncid, lon_id, "units", 12, "degree_east");
   nc_check_error(stat,__LINE__,__FILE__);
   lon_valid_min[0] = -180;
   stat = nc_put_att_double(ncid, lon_id, "valid_min", NC_DOUBLE, 1, lon_valid_min);
   nc_check_error(stat,__LINE__,__FILE__);
   lon_valid_max[0] = 180;
   stat = nc_put_att_double(ncid, lon_id, "valid_max", NC_DOUBLE, 1, lon_valid_max);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lon_id, "long_name", 9, "Longitude");
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
   stat = nc_put_att_text(ncid, lat_id, "long_name", 8, "Latitude");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lat_id, "standard_name", 8, "latitude");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lat_id, "nav_model", 12, "Default grid");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, spectrum_id, "long_name", 14, "tidal_spectrum");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Ha_id, "long_name", 19, "tidal_LSA_amplitude");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Ha_id, "short_name", 2, "Ha");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Ha_id, "units", 1, "m");
   nc_check_error(stat,__LINE__,__FILE__);
   Ha__FillValue[0] = -9999;
   stat = nc_put_att_float(ncid, Ha_id, "_FillValue", NC_FLOAT, 1, Ha__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   Ha_missing_value[0] = -9999;
   stat = nc_put_att_float(ncid, Ha_id, "missing_value", NC_FLOAT, 1, Ha_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   Ha_scale_factor[0] = 1;
   stat = nc_put_att_double(ncid, Ha_id, "scale_factor", NC_DOUBLE, 1, Ha_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   Ha_add_offset[0] = 0;
   stat = nc_put_att_double(ncid, Ha_id, "add_offset", NC_DOUBLE, 1, Ha_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Hg_id, "long_name", 19, "tidal_LSA_phase_lag");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Hg_id, "short_name", 2, "Hg");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Hg_id, "units", 7, "degree");
   nc_check_error(stat,__LINE__,__FILE__);
   Hg__FillValue[0] = -9999;
   stat = nc_put_att_float(ncid, Hg_id, "_FillValue", NC_FLOAT, 1, Hg__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   Hg_missing_value[0] = -9999;
   stat = nc_put_att_float(ncid, Hg_id, "missing_value", NC_FLOAT, 1, Hg_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   Hg_scale_factor[0] = 1;
   stat = nc_put_att_double(ncid, Hg_id, "scale_factor", NC_DOUBLE, 1, Hg_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   Hg_add_offset[0] = 0;
   stat = nc_put_att_double(ncid, Hg_id, "add_offset", NC_DOUBLE, 1, Hg_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, "CF 1.0");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "Topic", 13, "tidal loading");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "file_name", 6, "XXX.nc");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "production", 8, "POC 2004");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "history", 1, " ");
   nc_check_error(stat,__LINE__,__FILE__);

   /* leave define mode */
   stat = nc_enddef (ncid);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_close(ncid);
   nc_check_error(stat,__LINE__,__FILE__);
   return 0;

}


int create_headershort(char *output, grid_t grid, spectrum_t spectrum)
 {			

   int  ncid;			/* netCDF id */

   /* dimension ids */
   int X_dim;
   int Y_dim;
   int W_dim;
   int L_dim;

   /* dimension lengths */
   size_t X_len = grid.nx;
   size_t Y_len = grid.ny;
   size_t W_len = NC_UNLIMITED;
   size_t L_len = 10;

   /* variable ids */
   int lon_id;
   int lat_id;
   int spectrum_id;
   int depth_id;
   int Ha_id;
   int Hg_id;
   int Ua_id;
   int Ug_id;
   int Va_id;
   int Vg_id;

   /* variable shapes */
   int lon_dims[RANK_lon];
   int lat_dims[RANK_lat];
   int spectrum_dims[RANK_spectrum];
   int depth_dims[RANK_depth];
   int Ha_dims[RANK_Ha];
   int Hg_dims[RANK_Hg];
   int Ua_dims[RANK_Ua];
   int Ug_dims[RANK_Ug];
   int Va_dims[RANK_Va];
   int Vg_dims[RANK_Vg];

   /* attribute vectors */
   double lon_valid_min[1];
   double lon_valid_max[1];
   double lat_valid_min[1];
   double lat_valid_max[1];
   short depth_missing_value[1];
   short depth__FillValue[1];
   double depth_scale_factor[1];
   double depth_add_offset[1];
   short Ha__FillValue[1];
   short Ha_missing_value[1];
   double Ha_scale_factor[1];
   double Ha_add_offset[1];
   short Hg__FillValue[1];
   short Hg_missing_value[1];
   double Hg_scale_factor[1];
   double Hg_add_offset[1];
   short Ua__FillValue[1];
   short Ua_missing_value[1];
   double Ua_scale_factor[1];
   double Ua_add_offset[1];
   short Ug__FillValue[1];
   short Ug_missing_value[1];
   double Ug_scale_factor[1];
   double Ug_add_offset[1];
   short Va__FillValue[1];
   short Va_missing_value[1];
   double Va_scale_factor[1];
   double Va_add_offset[1];
   short Vg__FillValue[1];
   short Vg_missing_value[1];
   double Vg_scale_factor[1];
   double Vg_add_offset[1];

   /* enter define mode */
   int stat = nc_create(output, NC_CLOBBER, &ncid);
   nc_check_error(stat,__LINE__,__FILE__);

   /* define dimensions */
   stat = nc_def_dim(ncid, "X", X_len, &X_dim);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_def_dim(ncid, "Y", Y_len, &Y_dim);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_def_dim(ncid, "W", W_len, &W_dim);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_def_dim(ncid, "L", L_len, &L_dim);
   nc_check_error(stat,__LINE__,__FILE__);

   /* define variables */

   lon_dims[0] = Y_dim;
   lon_dims[1] = X_dim;
   stat = nc_def_var(ncid, "lon", NC_DOUBLE, RANK_lon, lon_dims, &lon_id);
   nc_check_error(stat,__LINE__,__FILE__);

   lat_dims[0] = Y_dim;
   lat_dims[1] = X_dim;
   stat = nc_def_var(ncid, "lat", NC_DOUBLE, RANK_lat, lat_dims, &lat_id);
   nc_check_error(stat,__LINE__,__FILE__);

   spectrum_dims[0] = W_dim;
   spectrum_dims[1] = L_dim;
   stat = nc_def_var(ncid, "spectrum", NC_CHAR, RANK_spectrum, spectrum_dims, &spectrum_id);
   nc_check_error(stat,__LINE__,__FILE__);

   depth_dims[0] = Y_dim;
   depth_dims[1] = X_dim;
   stat = nc_def_var(ncid, "depth", NC_SHORT, RANK_depth, depth_dims, &depth_id);
   nc_check_error(stat,__LINE__,__FILE__);

   Ha_dims[0] = W_dim;
   Ha_dims[1] = Y_dim;
   Ha_dims[2] = X_dim;
   stat = nc_def_var(ncid, "Ha", NC_SHORT, RANK_Ha, Ha_dims, &Ha_id);
   nc_check_error(stat,__LINE__,__FILE__);

   Hg_dims[0] = W_dim;
   Hg_dims[1] = Y_dim;
   Hg_dims[2] = X_dim;
   stat = nc_def_var(ncid, "Hg", NC_SHORT, RANK_Hg, Hg_dims, &Hg_id);
   nc_check_error(stat,__LINE__,__FILE__);

   Ua_dims[0] = W_dim;
   Ua_dims[1] = Y_dim;
   Ua_dims[2] = X_dim;
   stat = nc_def_var(ncid, "Ua", NC_SHORT, RANK_Ua, Ua_dims, &Ua_id);
   nc_check_error(stat,__LINE__,__FILE__);

   Ug_dims[0] = W_dim;
   Ug_dims[1] = Y_dim;
   Ug_dims[2] = X_dim;
   stat = nc_def_var(ncid, "Ug", NC_SHORT, RANK_Ug, Ug_dims, &Ug_id);
   nc_check_error(stat,__LINE__,__FILE__);

   Va_dims[0] = W_dim;
   Va_dims[1] = Y_dim;
   Va_dims[2] = X_dim;
   stat = nc_def_var(ncid, "Va", NC_SHORT, RANK_Va, Va_dims, &Va_id);
   nc_check_error(stat,__LINE__,__FILE__);

   Vg_dims[0] = W_dim;
   Vg_dims[1] = Y_dim;
   Vg_dims[2] = X_dim;
   stat = nc_def_var(ncid, "Vg", NC_SHORT, RANK_Vg, Vg_dims, &Vg_id);
   nc_check_error(stat,__LINE__,__FILE__);

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

   stat = nc_put_att_text(ncid, spectrum_id, "long_name", 14, "tidal_spectrum");
   nc_check_error(stat,__LINE__,__FILE__);

   stat = nc_put_att_text(ncid, depth_id, "long_name", 20, "model_positive_depth");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, depth_id, "short_name", 5, "depth");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, depth_id, "units", 1, "m");
   nc_check_error(stat,__LINE__,__FILE__);
   depth_missing_value[0] = -32768;
   stat = nc_put_att_short(ncid, depth_id, "missing_value", NC_SHORT, 1, depth_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   depth__FillValue[0] = -32768;
   stat = nc_put_att_short(ncid, depth_id, "_FillValue", NC_SHORT, 1, depth__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   depth_scale_factor[0] = 1;
   stat = nc_put_att_double(ncid, depth_id, "scale_factor", NC_DOUBLE, 1, depth_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   depth_add_offset[0] = 0;
   stat = nc_put_att_double(ncid, depth_id, "add_offset", NC_DOUBLE, 1, depth_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, depth_id, "associate", 7, "lat lon");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, depth_id, "content", 2, "YX");

   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Ha_id, "long_name", 25, "tidal_elevation_amplitude");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Ha_id, "short_name", 2, "Ha");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Ha_id, "units", 1, "m");
   nc_check_error(stat,__LINE__,__FILE__);
   Ha__FillValue[0] = -32768;
   stat = nc_put_att_short(ncid, Ha_id, "_FillValue", NC_SHORT, 1, Ha__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   Ha_missing_value[0] = -32768;
   stat = nc_put_att_short(ncid, Ha_id, "missing_value", NC_SHORT, 1, Ha_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   Ha_scale_factor[0] = 0.001;
   stat = nc_put_att_double(ncid, Ha_id, "scale_factor", NC_DOUBLE, 1, Ha_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   Ha_add_offset[0] = 0;
   stat = nc_put_att_double(ncid, Ha_id, "add_offset", NC_DOUBLE, 1, Ha_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Ha_id, "associate", 16, "spectrum lat lon");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Ha_id, "content", 3, "WYX");
   nc_check_error(stat,__LINE__,__FILE__);

   stat = nc_put_att_text(ncid, Hg_id, "long_name", 25, "tidal_elevation_phase_lag");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Hg_id, "short_name", 2, "Hg");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Hg_id, "units", 7, "degree");
   nc_check_error(stat,__LINE__,__FILE__);
   Hg__FillValue[0] = -32768;
   stat = nc_put_att_short(ncid, Hg_id, "_FillValue", NC_SHORT, 1, Hg__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   Hg_missing_value[0] = -32768;
   stat = nc_put_att_short(ncid, Hg_id, "missing_value", NC_SHORT, 1, Hg_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   Hg_scale_factor[0] = 0.1;
   stat = nc_put_att_double(ncid, Hg_id, "scale_factor", NC_DOUBLE, 1, Hg_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   Hg_add_offset[0] = 0;
   stat = nc_put_att_double(ncid, Hg_id, "add_offset", NC_DOUBLE, 1, Hg_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Hg_id, "associate", 16, "spectrum lat lon");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Hg_id, "content", 3, "WYX");
   nc_check_error(stat,__LINE__,__FILE__);
 
   stat = nc_put_att_text(ncid, Ua_id, "long_name", 32, "tidal_eastward_current_amplitude");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Ua_id, "short_name", 2, "Ua");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Ua_id, "units", 1, "m");
   nc_check_error(stat,__LINE__,__FILE__);
   Ua__FillValue[0] = -32768;
   stat = nc_put_att_short(ncid, Ua_id, "_FillValue", NC_SHORT, 1, Ua__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   Ua_missing_value[0] = -32768;
   stat = nc_put_att_short(ncid, Ua_id, "missing_value", NC_SHORT, 1, Ua_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   Ua_scale_factor[0] = 0.0001;
   stat = nc_put_att_double(ncid, Ua_id, "scale_factor", NC_DOUBLE, 1, Ua_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   Ua_add_offset[0] = 0;
   stat = nc_put_att_double(ncid, Ua_id, "add_offset", NC_DOUBLE, 1, Ua_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Ua_id, "associate", 16, "spectrum lat lon");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Ua_id, "content", 3, "WYX");
   nc_check_error(stat,__LINE__,__FILE__);

   stat = nc_put_att_text(ncid, Ug_id, "long_name", 32, "tidal_eastward_current_phase_lag");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Ug_id, "short_name", 2, "Ug");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Ug_id, "units", 7, "degree");
   nc_check_error(stat,__LINE__,__FILE__);
   Ug__FillValue[0] = -32768;
   stat = nc_put_att_short(ncid, Ug_id, "_FillValue", NC_SHORT, 1, Ug__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   Ug_missing_value[0] = -32768;
   stat = nc_put_att_short(ncid, Ug_id, "missing_value", NC_SHORT, 1, Ug_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   Ug_scale_factor[0] = 0.1;
   stat = nc_put_att_double(ncid, Ug_id, "scale_factor", NC_DOUBLE, 1, Ug_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   Ug_add_offset[0] = 0;
   stat = nc_put_att_double(ncid, Ug_id, "add_offset", NC_DOUBLE, 1, Ug_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Ug_id, "associate", 16, "spectrum lat lon");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Ug_id, "content", 3, "WYX");
   nc_check_error(stat,__LINE__,__FILE__);

   stat = nc_put_att_text(ncid, Va_id, "long_name", 34, "tidal_northtward_current_amplitude");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Va_id, "short_name", 2, "Va");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Va_id, "units", 1, "m");
   nc_check_error(stat,__LINE__,__FILE__);
   Va__FillValue[0] = -32768;
   stat = nc_put_att_short(ncid, Va_id, "_FillValue", NC_SHORT, 1, Va__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   Va_missing_value[0] = -32768;
   stat = nc_put_att_short(ncid, Va_id, "missing_value", NC_SHORT, 1, Va_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   Va_scale_factor[0] = 0.0001;
   stat = nc_put_att_double(ncid, Va_id, "scale_factor", NC_DOUBLE, 1, Va_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   Va_add_offset[0] = 0;
   stat = nc_put_att_double(ncid, Va_id, "add_offset", NC_DOUBLE, 1, Va_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Va_id, "associate", 16, "spectrum lat lon");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Va_id, "content", 3, "WYX");
   nc_check_error(stat,__LINE__,__FILE__);

   stat = nc_put_att_text(ncid, Vg_id, "long_name", 33, "tidal_northward_current_phase_lag");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Vg_id, "short_name", 2, "Vg");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Vg_id, "units", 7, "degree");
   nc_check_error(stat,__LINE__,__FILE__);
   Vg__FillValue[0] = -32768;
   stat = nc_put_att_short(ncid, Vg_id, "_FillValue", NC_SHORT, 1, Vg__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   Vg_missing_value[0] = -32768;
   stat = nc_put_att_short(ncid, Vg_id, "missing_value", NC_SHORT, 1, Vg_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   Vg_scale_factor[0] = 0.1;
   stat = nc_put_att_double(ncid, Vg_id, "scale_factor", NC_DOUBLE, 1, Vg_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   Vg_add_offset[0] = 0;
   stat = nc_put_att_double(ncid, Vg_id, "add_offset", NC_DOUBLE, 1, Vg_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Vg_id, "associate", 16, "spectrum lat lon");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Vg_id, "content", 3, "WYX");
   nc_check_error(stat,__LINE__,__FILE__);

   stat = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, "CF 1.0");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "Topic", 16, "barotropic tides");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "file_name", 6, "XXX.nc");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "production", 8, "POC 2004");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "history", 1, " ");
   nc_check_error(stat,__LINE__,__FILE__);

   /* leave define mode */
   stat = nc_enddef (ncid);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_close(ncid);
   nc_check_error(stat,__LINE__,__FILE__);
   return 0;
}

/*-----------------------------------------------------------------------------*/

int nc_writeaxis(int ncid, grid_t grid)

/*-----------------------------------------------------------------------------*/
{
  int  n,i,j,k,N,recordsize;
  int  status;
  size_t nitems;
  double *buffer;
  int varid;

  status=nc_inq_varid(ncid,"lon",&varid);
  status=nc_put_var_double(ncid,varid,grid.x);
  if (status != NC_NOERR) goto error;

  status=nc_inq_varid(ncid,"lat",&varid);
  status=nc_put_var_double(ncid,varid,grid.y);
  if (status != NC_NOERR) goto error;

  return (0);

 error:
  return(status);
}

/*-----------------------------------------------------------------------------*/

int nc_writeaxis_3d(int ncid, grid_t grid)

/*-----------------------------------------------------------------------------*/
{
  int  n,i,j,k,N,recordsize;
  int  status;
  size_t nitems;
  double *buffer;
  size_t start[3];
  size_t count[3];

  status=nc_put_var_double(ncid,0,grid.x);
  if (status != NC_NOERR) goto error;

  status=nc_put_var_double(ncid,1,grid.y);
  if (status != NC_NOERR) goto error;

  return (0);

 error:
  return(status);
}

/*-----------------------------------------------------------------------------*/

int nc_writespectrum(int ncid, spectrum_t spectrum)

/*-----------------------------------------------------------------------------*/
{
  int  k,wave;
  int  status;
  size_t nitems;
  size_t start[3];
  size_t count[3];
  char buffer[10];

  for(wave=0;wave<spectrum.n;wave++) {
    start[0]=wave;
    start[1]=0;
    count[0]=1;
    count[1]=10;
    strcpy(buffer,spectrum.waves[wave].name);
    status=nc_put_vara_text(ncid,2,start,count,buffer);
    if (status != NC_NOERR) goto error;
    }
  return (0);

error:
  nc_check_error(status,__LINE__,__FILE__);
  return(status);
}

/*-----------------------------------------------------------------------------*/

int nc_write_depth_2delete(int ncid, grid_t grid, int var, float *z, float mask)

/*-----------------------------------------------------------------------------*/
{
  int  n,N,recordsize;
  int  status;
  size_t nitems;
  short *buffer=NULL;
  size_t start[3];
  size_t count[3];
  float scale,offset;
  short missing;

  start[0]=0;
  start[1]=0;

  count[0]=grid.ny;
  count[1]=grid.nx;

  exitIfNull(
    buffer=(short *) malloc(grid.ny*grid.nx*sizeof(short))
    );

  status=nc_get_att_float(ncid,var,"scale_factor",&scale);
  if(status!=0) scale=1.0;

  status=nc_get_att_float(ncid,var,"add_offset",&offset);
  if(status!=0) offset=0.0;

  status=nc_get_att_short(ncid,var,"_FillValue",&missing);

  for(n=0;n<grid.ny*grid.nx;n++)
    if(z[n]!=mask) buffer[n]=(short)( (z[n]-offset)/scale );
    else           buffer[n]=missing;

  status=nc_put_vara_short(ncid,var,start,count,buffer);
  if (status != NC_NOERR) goto error;

  free(buffer);
  return (0);

 error:
  nc_check_error(status,__LINE__,__FILE__);
  return(status);
}

/*-----------------------------------------------------------------------------*/

int nc_write_r1(int ncid, grid_t grid, int var, int wave, float *z, float mask)

/*-----------------------------------------------------------------------------*/
{
  int  n,N,recordsize;
  int  status;
  size_t nitems;
  short *buffer=NULL;
  size_t start[3];
  size_t count[3];
  double scale,offset,zmax=0.0,range;
  short missing;



  start[0]=wave;
  start[1]=0;
  start[2]=0;

  count[0]=1;
  count[1]=grid.ny;
  count[2]=grid.nx;

  exitIfNull(
    buffer=(short *) malloc(grid.ny*grid.nx*sizeof(short))
    );

  status=nc_get_att_double(ncid,var,"scale_factor",&scale);
  if(status!=0) scale=1.0;

  for(n=0;n<grid.ny*grid.nx;n++)
    if(z[n]!=mask) zmax=MAX(zmax,fabs(z[n]));

  range=ceil(log(zmax/32500.)/log(10.));
  scale=exp(range*log(10.));
  status=nc_put_att_double(ncid,var,"scale_factor",NC_DOUBLE,1,&scale);

  status=nc_get_att_double(ncid,var,"add_offset",&offset);
  if(status!=0) offset=0.0;

  status=nc_get_att_short(ncid,var,"_FillValue",&missing);

  for(n=0;n<grid.ny*grid.nx;n++)
    if(z[n]!=mask) buffer[n]=(short)( floor((z[n]-offset)/scale+0.5) );
    else           buffer[n]=missing;

  status=nc_put_vara_short(ncid,var,start,count,buffer);
  if (status != NC_NOERR) goto error;

  free(buffer);
  return (0);

 error:
  nc_check_error(status,__LINE__,__FILE__);
  return(status);
}

/*-----------------------------------------------------------------------------*/

int nc_write_c1_2delete(int ncid, grid_t grid, int var, int wave, fcomplex *z)

/*-----------------------------------------------------------------------------*/
{
  int  n,N;
  int  status;
  size_t nitems;
  float *buffer=NULL;
  size_t start[3];
  size_t count[3];

  N=grid.nx*grid.ny;

  exitIfNull(
    buffer=(float *) malloc(N*sizeof(float))
    );

  for (n=0;n<N;n++) buffer[n]=real(z[n]);

  start[0]=wave;
  start[1]=0;
  start[2]=0;

  count[0]=1;
  count[1]=grid.ny;
  count[2]=grid.nx;

  status=nc_put_vara_float(ncid,var,start,count,buffer);
  if (status != NC_NOERR) goto error;

  free(buffer);

  return (0);

 error:
  return(status);
}

/*----------------------------------------------------------------------------*/
/// Lists the variables of an opened file
/**
\date last reviewed 28 Apr 2011
\author Damien Allain

\param ncid netCDF ID
\param *nvar number of variables found
\param ***name names of variables.
The array and its elements should be passed to free(3) to release the allocated storage when it is no longer needed.

\bug Why are so many arrays (like lengthhp) filled when the info is not used later ?
*/
/*----------------------------------------------------------------------------*/
int cdf_getvariable (int ncid,int *nvar, char ***name)
  {

  FILE *in=NULL;
  int option,status;//, ncid;
  size_t *lengthhp=NULL,index[1];
  float time;
  int ndimsp,nvarsp,ngattsp,unlimdimidp;
  int dim,var,vdata,vx,vy,vz,vt,att,nattsp;
  int *ndim=NULL,**dimids=NULL;
  nc_type xtypep;
  char **dimname=NULL,**varname=NULL;

  char bdum;
  short sdum;
  int idum;
  float fdum;
  double *ddum=NULL;

  /*for some reason, it was thought that one may only want the list of variables in a file and be done with it. I think one wants the names of the variables to actually read the file.
  status=nc_open(filename,0,&ncid);
  if(status != NC_NOERR) goto error;
  */

  status=nc_inq(ncid,&ndimsp,&nvarsp,&ngattsp,&unlimdimidp);
  if(status != NC_NOERR) goto no_go;

  printf("ncid %d ,ndimsp %d,nvarsp %d,ngattsp %d,unlimdimidp %d \n",
          ncid,ndimsp,nvarsp,ngattsp,unlimdimidp);

  exitIfNull(
    dimname=(char **) malloc(ndimsp*sizeof(char *))
    );
  exitIfNull(
    lengthhp=(size_t *) malloc(ndimsp*sizeof(size_t))
    );
  for (dim=0;dim<ndimsp;dim++) {
    exitIfNull(
      dimname[dim]=new char[NC_MAX_NAME+1]
      );
    status=nc_inq_dim(ncid,dim,dimname[dim],&lengthhp[dim]);
    if(status != NC_NOERR) goto error;
    printf("dimension %d,name %s, length %d \n",dim,dimname[dim],lengthhp[dim]);
    }

  exitIfNull(
    varname=(char **) malloc(nvarsp*sizeof(char *))
    );
  exitIfNull(
    dimids=(int **) malloc(nvarsp*sizeof(int *))
    );
  exitIfNull(
    ndim=(int *) malloc(nvarsp*sizeof(int))
    );

  for (var=0;var<nvarsp;var++) {
    exitIfNull(
      varname[var]=new char[NC_MAX_NAME+1]
      );
    status=nc_inq_varndims(ncid,var,&ndim[var]);
    if(status != NC_NOERR) goto error;
    exitIfNull(
      dimids[var]=(int *) malloc(ndim[var]*sizeof(int))
      );
    status=nc_inq_var(ncid,var,varname[var],&xtypep,&ndim[var],dimids[var],&nattsp);
    if(status != NC_NOERR) goto error;
    printf("var %d, name %s, type %d,ndim %d= ",
            var,varname[var],xtypep,ndim[var]);
    for(dim=0;dim< ndim[var];dim++) printf(" %d",dimids[var][dim]);
    printf(", natt %d\n", nattsp);
    }

  //status=nc_close(ncid);
  if(status != NC_NOERR) goto error;
  *name=varname;
  *nvar=nvarsp;
  return(status);

error:
  if(lengthhp != NULL) free(lengthhp);
  //status=nc_close(ncid);
  status=-1;
  return(status);

no_go:
  //status=nc_close(ncid);
  status=1;
  return(status);

unknown_format:
  //status=nc_close(ncid);
  status=2;
  return(status);

  }

  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_getshort_nt(const char *filename, mesh_t mesh, int frame, int id, float *z, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, N, recordsize;
  int status;
  size_t nitems;
  short *buffer=NULL;
  size_t start[3];
  size_t count[3];
  float scale, offset;
  short missing;
  double tmp;
  int ncid;

  start[0] = frame;
  start[1] = 0;

  count[0] = 1;
  count[1] = mesh.nvtxs;

  buffer = new short[mesh.nvtxs];

  status = nc_open(filename, NC_WRITE, &ncid);
  if(status != NC_NOERR)
    goto error;

  status = nc_get_att_float(ncid, id, "scale_factor", &scale);
  if(status != 0)
    scale = 1.0;

  status = nc_get_att_float(ncid, id, "add_offset", &offset);
  if(status != 0)
    offset = 0.0;

  status = nc_get_att_short(ncid, id, "_FillValue", &missing);

  status = nc_get_vara_short(ncid, id, start, count, buffer);
  if(status != NC_NOERR)
    goto error;

  *mask = 1.e+10;

  for(n = 0; n < mesh.nvtxs; n++)
    if(buffer[n] != missing) {
      z[n] = buffer[n] * scale + offset;
    } else
      z[n] = *mask;

  zaparr(buffer);

  status = nc_close(ncid);
  if(status != NC_NOERR)
    goto error;

  return (0);

error:
  nc_check_error(status, __LINE__, __FILE__);
  return (status);
}


  