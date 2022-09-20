
/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Fr�d�ric Dupont    Universit� de Laval � Qu�bec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#include <config.h>

#include <stdio.h>
#include <string.h>

#include "tools-structures.h"

#include "fe.def"
#include "fe.h"

#include "map.h"
#include "geo.h"
#include "swap.h"
#include "maths.h"
#include "poc-time.h"
#include "netcdf-proto.h"
#include "rutin.h"

#include "polygones.h"
#include "zapper.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <class Element> int fe_savemeshNC2D_template(const char *filename, Element & mesh, int option)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l, m, n;

  int ncid, stat;       /* netCDF id */

  /* dimension ids */
  int TRIANGLE_dim;
  int QUADRANGLE_dim;
  int N_dim;
  int E_dim;
  int P3_dim;
  int P4_dim;
  int T_dim;

  /* variable ids */
  int time_id;
  int lon_id;
  int lat_id;
  int triangles_id;
  int quadrangles_id;
  int elements_id;
  int bathymetry_id;
  int sigma_id;
  int depth_id;
  int surfelev_id;
  int ubar_id;
  int vbar_id;
  int U_id;
  int V_id;
  int W_id;
  int T_id;
  int S_id;
  int RHO_id;
  int TKE_id;
  int KH_id;
  /* dimension lengths */
  size_t TRIANGLE_len   = mesh.ntriangles;
  size_t QUADRANGLE_len = mesh.nquadrangles;
  size_t N_len = mesh.nvtxs;
  size_t E_len = mesh.nedges;
  size_t P3_len = 3;
  size_t P4_len = 4;
  size_t T_len = NC_UNLIMITED;

  /* rank (number of dimensions) for each variable */
#define RANK_lon 1
#define RANK_lat 1
#define RANK_element 2
#define RANK_bathymetry 1

  /* variable shapes */
  int lon_dims[RANK_lon];
  int lat_dims[RANK_lat];
  int element_dims[RANK_element];
  int bathymetry_dims[RANK_bathymetry];

  /* attribute vectors */
  double lon_valid_min[1];
  double lon_valid_max[1];
  double lat_valid_min[1];
  double lat_valid_max[1];
  double bathymetry_scale_factor[1];
  double bathymetry_add_offset[1];
  double sigma_scale_factor[1];
  double sigma_add_offset[1];
  double depth_scale_factor[1];
  double depth_add_offset[1];

  double *buffer;
  int *nv;
  char text[1024];

/*------------------------------------------------------------------------
  enter define mode */
  stat = nc_create(filename, NC_CLOBBER, &ncid);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------
  define dimensions */
  if(TRIANGLE_len!=0) {
    stat = nc_def_dim(ncid, "M", TRIANGLE_len, &TRIANGLE_dim);
    nc_check_error(stat, __LINE__, __FILE__);
    }
  stat = nc_def_dim(ncid, "NQUADRANGLES", QUADRANGLE_len, &QUADRANGLE_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "N", N_len, &N_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "E", E_len, &E_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "P3", P3_len, &P3_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "P4", P4_len, &P4_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "T", T_len, &T_dim);
  nc_check_error(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, "CF 1.0");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, NC_GLOBAL, "grid_type", 2, "UG");
  nc_check_error(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, NC_GLOBAL, "file_name", strlen(filename), filename);
  nc_check_error(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, NC_GLOBAL, "history", 1, " ");
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------
  leave define mode */
  if(option == 0) {
    stat = nc_enddef(ncid);
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_close(ncid);
    nc_check_error(stat, __LINE__, __FILE__);
    return (0);
    }

/*------------------------------------------------------------------------
   define other variables */

  lon_dims[0] = N_dim;
  stat = nc_def_var(ncid, "lon", NC_DOUBLE, RANK_lon, lon_dims, &lon_id);
  nc_check_error(stat, __LINE__, __FILE__);

  lat_dims[0] = N_dim;
  stat = nc_def_var(ncid, "lat", NC_DOUBLE, RANK_lat, lat_dims, &lat_id);
  nc_check_error(stat, __LINE__, __FILE__);

  if(TRIANGLE_len!=0) {
    element_dims[0] = TRIANGLE_dim;
    element_dims[1] = P3_dim;
    stat = nc_def_var(ncid, "triangles", NC_INT, RANK_element, element_dims, &triangles_id);
    nc_check_error(stat, __LINE__, __FILE__);
    }

  element_dims[0] = QUADRANGLE_dim;
  element_dims[1] = P4_dim;
  stat = nc_def_var(ncid, "quadrangles", NC_INT, RANK_element, element_dims, &quadrangles_id);
  nc_check_error(stat, __LINE__, __FILE__);

  /* assign attributes */
  stat = nc_put_att_text(ncid, lon_id, "units", 12, "degrees_east");
  nc_check_error(stat, __LINE__, __FILE__);
  lon_valid_min[0] = -180;
  stat = nc_put_att_double(ncid, lon_id, "valid_min", NC_DOUBLE, 1,lon_valid_min);
  nc_check_error(stat, __LINE__, __FILE__);

  lon_valid_max[0] = 180;
  stat = nc_put_att_double(ncid, lon_id, "valid_max", NC_DOUBLE, 1, lon_valid_max);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lon_id, "long_name", 9, "longitude");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lon_id, "standard_name", 9, "longitude");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lon_id, "subgrid", 5, "point");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lon_id, "content", 1, "N");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "units", 13, "degrees_north");
  nc_check_error(stat, __LINE__, __FILE__);

  lat_valid_min[0] = -90;
  stat = nc_put_att_double(ncid, lat_id, "valid_min", NC_DOUBLE, 1,lat_valid_min);
  nc_check_error(stat, __LINE__, __FILE__);
  lat_valid_max[0] = 90;
  stat = nc_put_att_double(ncid, lat_id, "valid_max", NC_DOUBLE, 1,lat_valid_max);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "long_name", 8, "latitude");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "standard_name", 8, "latitude");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "subgrid", 5, "point");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "content", 1, "N");
  nc_check_error(stat, __LINE__, __FILE__);

  if(TRIANGLE_len!=0) {
    elements_id=triangles_id;
    stat = nc_put_att_text(ncid, elements_id, "long_name", 20,"triangles_connectivity");
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid, elements_id, "standard_name", 9, "triangles");
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid, elements_id, "subgrid", 4, "cell");
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid, elements_id, "content", 2, "MP");
    nc_check_error(stat, __LINE__, __FILE__);
    }
    
  elements_id=quadrangles_id;
  stat = nc_put_att_text(ncid, elements_id, "long_name", 20,"quadrangles_connectivity");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, elements_id, "standard_name", 11, "quadrangles");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, elements_id, "subgrid", 4, "cell");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, elements_id, "content", 2, "MP");
  nc_check_error(stat, __LINE__, __FILE__);


  bathymetry_dims[0] = N_dim;
  stat = nc_def_var(ncid, "bathymetry", NC_FLOAT, RANK_bathymetry, bathymetry_dims, &bathymetry_id);
  nc_check_error(stat, __LINE__, __FILE__);

  strcpy(text,"model_positive_bathymetry_node");
  stat = nc_put_att_text(ncid, bathymetry_id, "long_name", strlen(text),text);
  nc_check_error(stat, __LINE__, __FILE__);

  strcpy(text,"bathymetry");
  stat = nc_put_att_text(ncid, bathymetry_id, "short_name", strlen(text),text);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "units", 1, "m");
  nc_check_error(stat, __LINE__, __FILE__);
  bathymetry_scale_factor[0] = 1;
  stat = nc_put_att_double(ncid, bathymetry_id, "scale_factor", NC_DOUBLE, 1,bathymetry_scale_factor);
  nc_check_error(stat, __LINE__, __FILE__);
  bathymetry_add_offset[0] = 0;
  stat = nc_put_att_double(ncid, bathymetry_id, "add_offset", NC_DOUBLE, 1,bathymetry_add_offset);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "associate", 7, "lon lat");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "subgrid", 5, "point");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "content", 1, "N");
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------
   leave define mode */
  stat = nc_enddef(ncid);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_close(ncid);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------
  write time independent fields in netcdf file */
  stat = nc_open(filename, NC_WRITE, &ncid);

  buffer = new double[N_len];

/*------------------------------------------------------------------
  longitude*/
  for(n = 0; n < N_len; n++)
    buffer[n] = mesh.vertices[n].lon;
  stat = nc_put_var_double(ncid, lon_id, buffer);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------
  latitude*/
  for(n = 0; n < N_len; n++)
    buffer[n] = mesh.vertices[n].lat;
  stat = nc_put_var_double(ncid, lat_id, buffer);
  nc_check_error(stat, __LINE__, __FILE__);

/* *------------------------------------------------------------------
  incidence for triangles*/
  if(TRIANGLE_len!=0) {
    elements_id=triangles_id;
    nv = new int[TRIANGLE_len * P3_len];
    for(m = 0; m < TRIANGLE_len; m++) {
      for(k = 0; k < P3_len; k++)
        nv[m * P3_len + k] = mesh.triangles[m].vertex[k];
      }
    stat = nc_put_var_int(ncid, elements_id, nv);
    nc_check_error(stat, __LINE__, __FILE__);
    zaparr(nv);
    }

/* *------------------------------------------------------------------
  incidence for quadrangles*/
  elements_id=quadrangles_id;
  nv = new int[QUADRANGLE_len * P4_len];
  for(m = 0; m < QUADRANGLE_len; m++) {
    for(k = 0; k < P4_len; k++)
      nv[m * P4_len + k] = mesh.quadrangles[m].vertex[k];
    }
  stat = nc_put_var_int(ncid, elements_id, nv);
  nc_check_error(stat, __LINE__, __FILE__);
  zaparr(nv);

/*------------------------------------------------------------------
  bathymétrie*/
  for(n = 0; n < N_len; n++)
    buffer[n] = mesh.vertices[n].h;
  stat = nc_put_var_double(ncid, bathymetry_id, buffer);
  nc_check_error(stat, __LINE__, __FILE__);
  zaparr(buffer);

  stat = nc_close(ncid);
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <class Element> int fe_savemeshNC3D_template(const char *filename, gmesh_t<Element> & mesh, int option)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *------------------------------------------------------------------------
  attribute name "axis" changed in "content" to comply with CF standard*/
{
  int k, l, m, n;

  int ncid, stat;       /* netCDF id */

  /* dimension ids */
  int M_dim;
  int N_dim;
  int E_dim;
  int P_dim;
  int K_dim;
  int L_dim;
  int T_dim;

  /* variable ids */
  int time_id;
  int lon_id;
  int lat_id;
  int bcode_id;
  int edges_id;
  int element_id;
  int bathymetry_id;
  int sigma_id;
  int depth_id;
  int surfelev_id;
  int ubar_id;
  int vbar_id;
  int U_id;
  int V_id;
  int W_id;
  int T_id;
  int S_id;
  int RHO_id;
  int TKE_id;
  int KH_id;
  /* dimension lengths */
  size_t M_len = mesh.nelts;
  size_t N_len = mesh.nvtxs;
  size_t E_len = mesh.nedges;
  size_t P_len = 4;
  size_t K_len = mesh.nlayers;
  size_t L_len = mesh.nlevels;
  size_t T_len = NC_UNLIMITED;

  /* rank (number of dimensions) for each variable */
#define RANK_time 1
#define RANK_lon 1
#define RANK_lat 1
#define RANK_element 2
#define RANK_edges 2
#define RANK_bathymetry 1
#define RANK_sigma 2
#define RANK_depth 2

  /* variable shapes */
  int time_dims[RANK_time];
  int lon_dims[RANK_lon];
  int lat_dims[RANK_lat];
  int element_dims[RANK_element];
  int edges_dims[RANK_element];
  int bathymetry_dims[RANK_bathymetry];
  int sigma_dims[RANK_sigma];
  int depth_dims[RANK_depth];

  /* attribute vectors */
  double lon_valid_min[1];
  double lon_valid_max[1];
  double lat_valid_min[1];
  double lat_valid_max[1];
  double bathymetry_scale_factor[1];
  double bathymetry_add_offset[1];
  double sigma_scale_factor[1];
  double sigma_add_offset[1];
  double depth_scale_factor[1];
  double depth_add_offset[1];

  double *buffer;
  int *nv;
  char text[1024];

/*------------------------------------------------------------------------
   enter define mode */
  stat = nc_create(filename, NC_CLOBBER, &ncid);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------
   define dimensions */
  stat = nc_def_dim(ncid, "M", M_len, &M_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "N", N_len, &N_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "E", E_len, &E_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "P", P_len, &P_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "K", K_len, &K_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "L", L_len, &L_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "T", T_len, &T_dim);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------
   define time variable */

  time_dims[0] = T_dim;
  stat = nc_def_var(ncid, "time", NC_DOUBLE, RANK_time, time_dims, &time_id);
  nc_check_error(stat, __LINE__, __FILE__);

  bathymetry_dims[0] = N_dim;
  stat = nc_def_var(ncid, "bathymetry", NC_FLOAT, RANK_bathymetry,
                    bathymetry_dims, &bathymetry_id);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------
   /* assign minimum attributes */
  stat = nc_put_att_text(ncid, time_id, "units", 7, "seconds");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, time_id, "calendar", 9, "gregorian");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, time_id, "long_name", 41, "Time elasped in seconds since time_origin");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, time_id, "title", 4, "Time");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, time_id, "time_origin", 20, "1950-JAN-01 00:00:00");
  nc_check_error(stat, __LINE__, __FILE__);

  strcpy(text,"model_positive_bathymetry_node");
  stat = nc_put_att_text(ncid, bathymetry_id, "long_name", strlen(text),text);
  nc_check_error(stat, __LINE__, __FILE__);
  
  strcpy(text,"bathymetry");
  stat = nc_put_att_text(ncid, bathymetry_id, "short_name", strlen(text),text);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "units", 1, "m");
  nc_check_error(stat, __LINE__, __FILE__);
  bathymetry_scale_factor[0] = 1;
  stat = nc_put_att_double(ncid, bathymetry_id, "scale_factor", NC_DOUBLE, 1,bathymetry_scale_factor);
  nc_check_error(stat, __LINE__, __FILE__);
  bathymetry_add_offset[0] = 0;
  stat = nc_put_att_double(ncid, bathymetry_id, "add_offset", NC_DOUBLE, 1,bathymetry_add_offset);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "associate", 7, "lon lat");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "subgrid", 5, "point");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "content", 1, "N");
  nc_check_error(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, "CF 1.0");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, NC_GLOBAL, "grid_type", 2, "UG");
  nc_check_error(stat, __LINE__, __FILE__);
  
  strcpy(text,"T-UGOm archive");
  stat = nc_put_att_text(ncid, NC_GLOBAL, "Topic", strlen(text),text);
  nc_check_error(stat, __LINE__, __FILE__);
  
  stat = nc_put_att_text(ncid, NC_GLOBAL, "file_name", strlen(filename), filename);
  nc_check_error(stat, __LINE__, __FILE__);
  
  strcpy(text,"F. Lyard/POC 2007");
  stat = nc_put_att_text(ncid, NC_GLOBAL, "production", strlen(text),text);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, NC_GLOBAL, "history", 1, " ");
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------
   leave define mode */
  if(option == 0) {
    stat = nc_enddef(ncid);
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_close(ncid);
    nc_check_error(stat, __LINE__, __FILE__);
    return (0);
  }

/*------------------------------------------------------------------------
   define other variables */

  lon_dims[0] = N_dim;
  stat = nc_def_var(ncid, "lon", NC_DOUBLE, RANK_lon, lon_dims, &lon_id);
  nc_check_error(stat, __LINE__, __FILE__);

  lat_dims[0] = N_dim;
  stat = nc_def_var(ncid, "lat", NC_DOUBLE, RANK_lat, lat_dims, &lat_id);
  nc_check_error(stat, __LINE__, __FILE__);

  element_dims[0] = M_dim;
  element_dims[1] = P_dim;
  stat = nc_def_var(ncid, "element", NC_INT, RANK_element, element_dims,&element_id);
  nc_check_error(stat, __LINE__, __FILE__);

  edges_dims[0] = M_dim;
  edges_dims[1] = P_dim;
  stat = nc_def_var(ncid, "edges", NC_INT, RANK_edges, edges_dims,&edges_id);
  nc_check_error(stat, __LINE__, __FILE__);

  sigma_dims[0] = N_dim;
  sigma_dims[1] = L_dim;
  stat = nc_def_var(ncid, "sigma", NC_FLOAT, RANK_sigma, sigma_dims,&sigma_id);
  nc_check_error(stat, __LINE__, __FILE__);

  depth_dims[0] = N_dim;
  depth_dims[1] = L_dim;
  stat = nc_def_var(ncid, "depth", NC_FLOAT, RANK_depth, depth_dims,
                    &depth_id);
  nc_check_error(stat, __LINE__, __FILE__);

  /* assign attributes */
  stat = nc_put_att_text(ncid, lon_id, "units", 12, "degrees_east");
  nc_check_error(stat, __LINE__, __FILE__);
  lon_valid_min[0] = -180;
  stat = nc_put_att_double(ncid, lon_id, "valid_min", NC_DOUBLE, 1,
                           lon_valid_min);
  nc_check_error(stat, __LINE__, __FILE__);

  lon_valid_max[0] = 180;
  stat = nc_put_att_double(ncid, lon_id, "valid_max", NC_DOUBLE, 1,
                           lon_valid_max);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lon_id, "long_name", 9, "longitude");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lon_id, "standard_name", 9, "longitude");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lon_id, "subgrid", 5, "point");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lon_id, "content", 1, "N");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "units", 13, "degrees_north");
  nc_check_error(stat, __LINE__, __FILE__);

  lat_valid_min[0] = -90;
  stat = nc_put_att_double(ncid, lat_id, "valid_min", NC_DOUBLE, 1,lat_valid_min);
  nc_check_error(stat, __LINE__, __FILE__);
  lat_valid_max[0] = 90;
  stat = nc_put_att_double(ncid, lat_id, "valid_max", NC_DOUBLE, 1,lat_valid_max);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "long_name", 8, "latitude");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "standard_name", 8, "latitude");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "subgrid", 5, "point");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "content", 1, "N");
  nc_check_error(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, element_id, "long_name", 20, "element_connectivity");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, element_id, "standard_name", 7, "element");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, element_id, "subgrid", 4, "cell");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, element_id, "content", 2, "MP");
  nc_check_error(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, edges_id, "long_name", 20, "element_edges");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, edges_id, "standard_name", 7, "edges");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, edges_id, "subgrid", 4, "cell");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, edges_id, "content", 2, "MP");
  nc_check_error(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, bcode_id, "long_name", strlen("boundary_code"), "boundary_code");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bcode_id, "standard_name", strlen("boundary_code"), "boundary_code");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bcode_id, "subgrid", 4, "cell");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bcode_id, "content", 1, "N");
  nc_check_error(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, sigma_id, "long_name", 37,"generalized sigma coordinate at nodes");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, sigma_id, "standard_name", 35,"general_ocean_sigma_coordinate_node");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, sigma_id, "short_name", 5, "sigma");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, sigma_id, "units", 1, "1");
  nc_check_error(stat, __LINE__, __FILE__);
  sigma_scale_factor[0] = 1;
  stat = nc_put_att_double(ncid, sigma_id, "scale_factor", NC_DOUBLE, 1,
                           sigma_scale_factor);
  nc_check_error(stat, __LINE__, __FILE__);
  sigma_add_offset[0] = 0;
  stat = nc_put_att_double(ncid, sigma_id, "add_offset", NC_DOUBLE, 1,
                           sigma_add_offset);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, sigma_id, "associate", 7, "lon lat");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, sigma_id, "subgrid", 5, "point");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, sigma_id, "content", 2, "NZ");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, sigma_id, "positive", 2, "up");
  nc_check_error(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, depth_id, "long_name", 37, "generalized depth coordinate at nodes");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, depth_id, "standard_name", 35, "general_ocean_depth_coordinate_node");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, depth_id, "short_name", 5, "depth");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, depth_id, "units", 1, "1");
  nc_check_error(stat, __LINE__, __FILE__);
  depth_scale_factor[0] = 1;
  stat = nc_put_att_double(ncid, depth_id, "scale_factor", NC_DOUBLE, 1,depth_scale_factor);
  nc_check_error(stat, __LINE__, __FILE__);
  depth_add_offset[0] = 0;
  stat = nc_put_att_double(ncid, depth_id, "add_offset", NC_DOUBLE, 1, depth_add_offset);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, depth_id, "associate", 7, "lon lat");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, depth_id, "subgrid", 5, "point");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, depth_id, "content", 2, "NZ");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, depth_id, "positive", 2, "up");
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------
   leave define mode */
  stat = nc_enddef(ncid);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_close(ncid);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------
  write time independent fields in netcdf file */
  stat = nc_open(filename, NC_WRITE, &ncid);

  buffer = new double[N_len];

/*------------------------------------------------------------------
  longitude*/
  for(n = 0; n < N_len; n++)
    buffer[n] = mesh.vertices[n].lon;
  stat = nc_put_var_double(ncid, lon_id, buffer);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------
  latitude*/
  for(n = 0; n < N_len; n++)
    buffer[n] = mesh.vertices[n].lat;
  stat = nc_put_var_double(ncid, lat_id, buffer);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------
  incidence*/
  nv = new int[M_len * P_len];
  for(m = 0; m < M_len; m++) {
    for(k = 0; k < P_len; k++)
      nv[m * P_len + k] = mesh.elements[m].vertex[k];
    }
  stat = nc_put_var_int(ncid, element_id, nv);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------
  edges*/
  for(m = 0; m < M_len; m++) {
    for(k = 0; k < P_len; k++)
      nv[m * P_len + k] = mesh.elements[m].edges[k];
    }
  stat = nc_put_var_int(ncid, edges_id, nv);
  nc_check_error(stat, __LINE__, __FILE__);
  zaparr(nv);

/*------------------------------------------------------------------
  boundary codes*/
  nv = new int[N_len];
  for(m = 0; m < N_len; m++) {
    nv[m] = mesh.vertices[m].code;
    }
  stat = nc_put_var_int(ncid, bcode_id, nv);
  nc_check_error(stat, __LINE__, __FILE__);
  zaparr(nv);

/*------------------------------------------------------------------
  bathymétrie*/
  for(n = 0; n < N_len; n++)
    buffer[n] = mesh.vertices[n].h;
  stat = nc_put_var_double(ncid, bathymetry_id, buffer);
  nc_check_error(stat, __LINE__, __FILE__);
  zaparr(buffer);

  buffer = new double[N_len*L_len];
/*------------------------------------------------------------------
  initial sigma levels*/
  for(n = 0; n < N_len; n++)
    for(l = 0; l < L_len; l++)
      buffer[n*L_len+l] = mesh.vertices[n].sigma[l];
  stat = nc_put_var_double(ncid, sigma_id, buffer);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------
  initial depth levels*/
  for(n = 0; n < N_len; n++) {
    if(mesh.vertices[n].zlevels==0) {
      for(l = 0; l < L_len; l++)
        buffer[n*L_len+l] = -mesh.vertices[n].sigma[l]*mesh.vertices[n].h;
      }
    else {
      for(l = 0; l < L_len; l++)
        buffer[n*L_len+l] = mesh.vertices[n].zlevels[l];
      }
    }
  stat = nc_put_var_double(ncid, depth_id, buffer);
  nc_check_error(stat, __LINE__, __FILE__);
  zaparr(buffer);

  stat = nc_close(ncid);
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_savemeshNC2D_new(const char *filename, mesh_t & mesh, int option)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_savemeshNC2D_template(filename, mesh, option);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_savemeshNC3D_new(const char *filename, gmesh_t<quadrangle_t> & mesh, int option)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_savemeshNC3D_template(filename, mesh, option);
  }

