
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
#include <cmath>
#include <netcdf.h>

#define pi M_PI

#ifdef PARALLEL
#ifdef HAVE_MPI
#include "mpi.h"
#endif
#endif

#include "tools-structures.h"

#include "constants.h"
#include "fe.h"
#include "map.h"
#include "poc-time.h"
#include "netcdf-proto.h"
#include "poc-netcdf.def"
#include "zapper.h"

static int first_access=1;
static int *incidence = NULL;
static char *rootname=0, *localname, *zone;
static grid_t ncgrid;
static double local_regular_sampling;
static int targets[1000];

static int lon_id, lat_id;
static int time_id;
static int h_id,z_id,ubar_id,vbar_id,Hubar_id,Hvbar_id;
static int zmin_id,zmax_id,MKE2D_id;
static int Ex2D_id,Ey2D_id,Kh2D_id;
static int Ax2D_id,Ay2D_id;
static int divHu2D_id,dhdt2D_id,omega2D_id,w2D_id;

static int Cd_id,z0_id,He_H_id;
static int lmts_id;

static int uwind_id,vwind_id,pa_id,wsx_id,wsy_id;

static int u_id,v_id,w_id,omega_id,omicro_id,dhdt_id,dsdt_id;
static int ibd_id, qkr_id;
static int T_id, S_id, RHO_id, p_id;

static int px_id,py_id;
static int qx_id,qy_id;
static int Ex_id,Ey_id,Kh_id;
static int Exmean3D_id,Eymean3D_id;

static int Ax_id,Ay_id;
static int Axmean3D_id,Aymean3D_id;

static int ubar3D_id,vbar3D_id;

static int sxx_id, sxy_id, syy_id;

static int divSxx_u_id, divSxx_v_id, divdivSxx_id;

static int tracer2D_id;

static int cycles_id;

#define SEA_FLOOR_DEPTH                 0
#define SEA_SURFACE_ELEVATION           1
#define BAROTROPIC_CURRENTS             2
#define BAROTROPIC_TRANSPORT            3
#define INVERSE_BAROMETER_DEPARTURE     4
#define OCEAN_BOTTOM_MOTION             5
#define BAROTROPIC_MOMENTUM_DIFFUSION   6
#define BAROTROPIC_MOMENTUM_ADVECTION   7
#define BAROTROPIC_TRANSPORT_DIVERGENCE 8


#define BAROCLINIC_CURRENTS              10
#define SEA_WATER_TEMPERATURE            11
#define SEA_WATER_SALINITY               12
#define SEA_WATER_DENSITY                13
#define BAROCLINIC_PRESSURE              14
#define PRESSURE_GRADIENT                15
#define MOMENTUM_DIFFUSION               16
#define MOMENTUM_DIFFUSION_COEFFICIENT   17
#define MOMENTUM_DIFFUSION_COEFFICIENT2D 18
#define LATERAL_MIXING                   19
#define MOMENTUM_ADVECTION               20

#define BOTTOM_FRICTION_COEFFICIENT      30
#define BOTTOM_RUGOSITY                  31
#define RELATIVE_BOTTOM_LAYER            32

#define SEA_SURFACE_ELEVATION_EXTREMA    40
#define MEAN_BAROTROPIC_KINETIC_ENERGY   41

#define MEAN_BAROCLINIC_CURRENTS               50
#define MEAN_BAROCLINIC_MOMENTUM_DIFFUSION     51
#define MEAN_BAROCLINIC_MOMENTUM_ADVECTION     52

#define WIND_STRESS                      60
#define ATMOSPHERIC_PRESSURE             61

#define RADIATION_STRESS		70
#define DIV_SXX				71
#define DIVDIV_SXX			72

#define TRACER_2D                       100

#define CYCLES                          200

extern char *sgetnewdate(date_t reference, double t);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_savemesh3d(const char *filename, mesh_t & mesh, int option, bool edges, bool levels, bool code, bool bathymetry)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l, m, n;

  int ncid, stat;       /* netCDF id */

  /* dimension ids */
  int M_dim;
  int Q_dim;
  int N_dim;
  int E_dim;
  int P_dim;
  int P2_dim;
  int K_dim;
  int L_dim;
  int L2D_dim;
  int T_dim;

  /* variable ids */
  int time_id;
  int lon_id;
  int lat_id;
  int edges_id;
  int CLGP2nodes_id;
  int element_id;
  int LGP1code_id, zCode_id, uCode_id;
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
  size_t M_len = mesh.ntriangles;
  size_t Q_len = mesh.nquadrangles;
  size_t N_len = mesh.nvtxs;
  size_t E_len = mesh.nedges;
//  size_t P_len = 3;
  size_t P_len = -1;
  size_t P2_len = 6;
  size_t K_len = mesh.nlayers;
  size_t L_len = mesh.nlevels;
  size_t L2D_len = 2;
  size_t T_len = NC_UNLIMITED;

  /* rank (number of dimensions) for each variable */
#define RANK_time 1
#define RANK_lon 1
#define RANK_lat 1
#define RANK_element 2
#define RANK_edges 2
#define RANK_connectivity 2
#define RANK_bathymetry 1
#define RANK_sigma 2
#define RANK_depth 2

  /* variable shapes */
  int time_dims[RANK_time];
  int lon_dims[RANK_lon];
  int lat_dims[RANK_lat];
  int element_dims[RANK_element];
  int edges_dims[RANK_edges];
  int CLGP2nodes_dims[RANK_connectivity];
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

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  enter define mode
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  stat = nc_create(filename, NC_CLOBBER, &ncid);
  nc_check_error(stat, __LINE__, __FILE__);
  if(stat!=0) {
    printf("cannot create %s\n",filename);
    return(-1);
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  define dimensions 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(M_len!=0){
    P_len=3;
    stat = nc_def_dim(ncid, "M", M_len, &M_dim);
    nc_check_error(stat, __LINE__, __FILE__);
    }
  if(Q_len!=0){
    P_len=4;
    stat = nc_def_dim(ncid, "Q", Q_len, &Q_dim);
    nc_check_error(stat, __LINE__, __FILE__);
    }
  
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
  stat = nc_def_dim(ncid, "L2D", L2D_len, &L2D_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "T", T_len, &T_dim);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------
  define time variable */
  time_dims[0] = T_dim;
  stat = nc_def_var(ncid, "time", NC_DOUBLE, RANK_time, time_dims, &time_id);
  nc_check_error(stat, __LINE__, __FILE__);

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

/*------------------------------------------------------------------------
  define bathymetry variable */
  if(bathymetry) {
    bathymetry_dims[0] = N_dim;
    stat = nc_def_var(ncid, "bathymetry", NC_FLOAT, RANK_bathymetry,  bathymetry_dims, &bathymetry_id);
    nc_check_error(stat, __LINE__, __FILE__);
    }
    
  stat = nc_put_att_text(ncid, bathymetry_id, "long_name", 25, "model_positive_bathymetry_node");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "short_name", 6, "bathymetryv");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "units", 1, "m");
  nc_check_error(stat, __LINE__, __FILE__);
  bathymetry_scale_factor[0] = 1;
  stat = nc_put_att_double(ncid, bathymetry_id, "scale_factor", NC_DOUBLE, 1, bathymetry_scale_factor);
  nc_check_error(stat, __LINE__, __FILE__);
  bathymetry_add_offset[0] = 0;
  stat = nc_put_att_double(ncid, bathymetry_id, "add_offset", NC_DOUBLE, 1, bathymetry_add_offset);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "associate", 7, "lon lat");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "subgrid", 5, "point");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "content", 1, "N");
  nc_check_error(stat, __LINE__, __FILE__);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  global attributes
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  stat = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, "CF 1.0");
  nc_check_error(stat, __LINE__, __FILE__);
  
  strcpy(text,"SIROCCO");
  stat = nc_put_att_text(ncid, NC_GLOBAL, "Extensions", strlen(text),text);
  nc_check_error(stat, __LINE__, __FILE__);
  
  stat = nc_put_att_text(ncid, NC_GLOBAL, "grid_type", 2, "UG");
  nc_check_error(stat, __LINE__, __FILE__);

  strcpy(text,"T-UGOm archive/restart");
  stat = nc_put_att_text(ncid, NC_GLOBAL, "Topic", strlen(text),text);
  nc_check_error(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, NC_GLOBAL, "file_name", strlen(filename), filename);
  nc_check_error(stat, __LINE__, __FILE__);

  sprintf(text,"made by %s around line %d of " __FILE__ " in T-UGOm 2011",__FUNCTION__,__LINE__);
  stat = nc_put_att_text(ncid, NC_GLOBAL, "production", strlen(text),text);
  nc_check_error(stat, __LINE__, __FILE__);
  
  stat = nc_put_att_text(ncid, NC_GLOBAL, "history", 1, " ");
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------
  leave define mode and return */
  if(option == 0) {
    stat = nc_enddef(ncid);
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_close(ncid);
    nc_check_error(stat, __LINE__, __FILE__);
    return (0);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  define other variables
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  lon_dims[0] = N_dim;
  stat = nc_def_var(ncid, "lon", NC_DOUBLE, RANK_lon, lon_dims, &lon_id);
  nc_check_error(stat, __LINE__, __FILE__);

  lat_dims[0] = N_dim;
  stat = nc_def_var(ncid, "lat", NC_DOUBLE, RANK_lat, lat_dims, &lat_id);
  nc_check_error(stat, __LINE__, __FILE__);

  if(M_len!=0) {
    element_dims[0] = M_dim;
    element_dims[1] = P_dim;
    stat = nc_def_var(ncid, "triangles", NC_INT, RANK_element, element_dims, &element_id);
    nc_check_error(stat, __LINE__, __FILE__);
    if(edges) {
      edges_dims[0] = M_dim;
      edges_dims[1] = P_dim;
      stat = nc_def_var(ncid, "edges", NC_INT, RANK_edges, edges_dims,&edges_id);
      nc_check_error(stat, __LINE__, __FILE__);
      }
    }

  if(Q_len!=0) {
    element_dims[0] = Q_dim;
    element_dims[1] = P_dim;
    stat = nc_def_var(ncid, "quadrangles", NC_INT, RANK_element, element_dims, &element_id);
    nc_check_error(stat, __LINE__, __FILE__);
    if(edges) {
      edges_dims[0] = Q_dim;
      edges_dims[1] = P_dim;
      stat = nc_def_var(ncid, "edges", NC_INT, RANK_edges, edges_dims,&edges_id);
      nc_check_error(stat, __LINE__, __FILE__);
      }
    }

  if(code) {
    lat_dims[0] = N_dim;
    stat = nc_def_var(ncid, "LGP1code", NC_DOUBLE, RANK_lat, lat_dims, &LGP1code_id);
    nc_check_error(stat, __LINE__, __FILE__);
    }
    
  if(levels) {
    depth_dims[0] = N_dim;
    depth_dims[1] = L_dim;
    stat = nc_def_var(ncid, "depth", NC_FLOAT, RANK_depth, depth_dims, &depth_id);
    nc_check_error(stat, __LINE__, __FILE__);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  assign attributes 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

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

  stat = nc_put_att_text(ncid, element_id, "long_name", 20, "element_connectivity");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, element_id, "standard_name", 7, "element");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, element_id, "subgrid", 4, "cell");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, element_id, "content", 2, "MP");
  nc_check_error(stat, __LINE__, __FILE__);

  if(edges) {
    stat = nc_put_att_text(ncid, edges_id, "long_name", 13, "element_edges");
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid, edges_id, "standard_name", 5, "edges");
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid, edges_id, "subgrid", 4, "cell");
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid, edges_id, "content", 2, "MP");
    nc_check_error(stat, __LINE__, __FILE__);
    }
    
  if(code) {
    stat = nc_put_att_text(ncid, LGP1code_id, "long_name", strlen("LGP1_boundary_code"), "LGP1_boundary_code");
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid, LGP1code_id, "standard_name", strlen("LGP1_boundary_code"), "LGP1_boundary_code");
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid, LGP1code_id, "subgrid", 4, "cell");
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid, LGP1code_id, "content", 1, "N");
    nc_check_error(stat, __LINE__, __FILE__);
    }

  if(mesh.vertices[0].zlevels!=NULL and levels) {
    stat = nc_put_att_text(ncid, depth_id, "long_name", 37, "generalized depth coordinate at nodes");
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid, depth_id, "standard_name", 35, "general_ocean_depth_coordinate_node");
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid, depth_id, "short_name", 5, "depth");
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid, depth_id, "units", 1, "1");
    nc_check_error(stat, __LINE__, __FILE__);
    depth_scale_factor[0] = 1;
    stat = nc_put_att_double(ncid, depth_id, "scale_factor", NC_DOUBLE, 1, depth_scale_factor);
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
    }
  
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

  for(n = 0; n < N_len; n++)
    buffer[n] = mesh.vertices[n].lon;
  stat = nc_put_var_double(ncid, lon_id, buffer);
  nc_check_error(stat, __LINE__, __FILE__);

  for(n = 0; n < N_len; n++)
    buffer[n] = mesh.vertices[n].lat;
  stat = nc_put_var_double(ncid, lat_id, buffer);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------
  vertex and edges incidence*/
  if(M_len!=0) {
    nv = new int[M_len * P_len];
    for(m = 0; m < M_len; m++) {
      for(k = 0; k < P_len; k++)
        nv[m * P_len + k] = mesh.triangles[m].vertex[k];
      }
    stat = nc_put_var_int(ncid, element_id, nv);
    nc_check_error(stat, __LINE__, __FILE__);
    if(edges) {
      for(m = 0; m < M_len; m++) {
        for(k = 0; k < P_len; k++)
          nv[m * P_len + k] = mesh.triangles[m].edges[k];
        }
      stat = nc_put_var_int(ncid, edges_id, nv);
      nc_check_error(stat, __LINE__, __FILE__);
      }
    zaparr(nv);
    }
    
  if(Q_len!=0) {
    nv = new int[Q_len * P_len];
    for(m = 0; m < Q_len; m++) {
      for(k = 0; k < P_len; k++)
        nv[m * P_len + k] = mesh.quadrangles[m].vertex[k];
        }
    stat = nc_put_var_int(ncid, element_id, nv);
    nc_check_error(stat, __LINE__, __FILE__);
    if(edges) {
      for(m = 0; m < Q_len; m++) {
        for(k = 0; k < P_len; k++)
          nv[m * P_len + k] = mesh.quadrangles[m].edges[k];
        }
      stat = nc_put_var_int(ncid, edges_id, nv);
      nc_check_error(stat, __LINE__, __FILE__);
      }
    zaparr(nv);
    }
    
/*------------------------------------------------------------------
  boundary codes*/
  if(code) {
    nv = new int[N_len];
    for(m = 0; m < N_len; m++) {
      nv[m] = mesh.vertices[m].code;
      }
    stat = nc_put_var_int(ncid, LGP1code_id, nv);
    nc_check_error(stat, __LINE__, __FILE__);
    zaparr(nv);
    }

/*------------------------------------------------------------------
  bathymétrie*/
  if(bathymetry) {
    for(n = 0; n < N_len; n++)
      buffer[n] = mesh.vertices[n].h;
    stat = nc_put_var_double(ncid, bathymetry_id, buffer);
    nc_check_error(stat, __LINE__, __FILE__);
    }

  zaparr(buffer);


  if(mesh.vertices[0].zlevels!=NULL and levels) {
    buffer = new double[N_len*L_len];
/*------------------------------------------------------------------
    initial depth levels*/
    for(n=0;n<N_len;n++) {
      for(l=0;l<L_len;l++) {
        buffer[n*L_len+l]=mesh.vertices[n].zlevels[l];
        }
      }
    stat=nc_put_var_double(ncid,depth_id,buffer);
    nc_check_error(stat,__LINE__,__FILE__);
    zaparr(buffer);
    }

  stat = nc_close(ncid);
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_savemesh3d(const char *filename, mesh_t & mesh, int option)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  bool levels=true, code=true, bathymetry=true, edges=true;
  
  status=fe_savemesh3d(filename, mesh, option, edges, levels, code, bathymetry);
  
  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int archiving_UGdummy3D(const char *localname, mesh_t & mesh,const char *name, double **buffer, int frame, int h_discretisation, int v_discretisation, int l_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int i,k,l,m,n;
  double *nvalues;
  FILE *in;
  double mask=1.e+10;
  int count,ncid,variable_id;
  int file_exist;
  cdfvar_t variable;
  double scale,offset;
  char filename[1024];
  char *output_path="./";
  
  if(rootname==0) rootname=strdup(output_path);
  char *s=strstr( (char *) localname, (char *) rootname);
  if(s==0) {
    sprintf(filename, "%s/%s", rootname, localname);
    }
  else {
//     printf("\n warning, output path found in filename, this is not a good practice\n");
//     printf("output path :%s \n",rootname);
//     printf("filename    :%s \n\n",localname);
    sprintf(filename, "%s", localname);
    }
/*------------------------------------------------------------------------
  see if file exist*/
  in = fopen(filename, "r");
  if(in == NULL) {
    file_exist=0;
    }
  else {
    file_exist=1;
    }

  if(file_exist == 0) {
    if(mesh.nlayers==0) mesh.nlayers=1;
    if(mesh.nlevels==0) mesh.nlevels=2;
    status = fe_savemesh3d(filename, mesh, 1);
    }
/**----------------------------------------------------------------------------
  temporary, for plotting in xscan (to be fixed in xscan) */
  status =fe_savediscretisation(filename, mesh, mesh.LGP1descriptor);

  switch (h_discretisation) {
    case LGP0:
      status =fe_savediscretisation(filename, mesh, mesh.LGP0descriptor);
      break;
    case LGP1:
//      status =fe_savediscretisation(filename, mesh, mesh.LGP1descriptor);
      break;
    case DGP1:
      status =fe_savediscretisation(filename, mesh, mesh.DGP1descriptor);
      break;
    case NCP1:
      break;
    case DNP1:
      status =fe_savediscretisation(filename, mesh, mesh.DNP1descriptor);
      break;
    case LGP2:
      status =fe_savediscretisation(filename, mesh, mesh.LGP2descriptor);
      break;
      }

/**------------------------------------------------------------------------

  2D section - 2D section - 2D section - 2D section - 2D section - 2D section

-------------------------------------------------------------------------*/

  count=frame;

  scale=1.0;
  offset=0.;

  variable =poc_variable_UG4D( name, mask, "m", scale, offset, name,"T",h_discretisation,v_discretisation);
//  variable =poc_variable_UG4D( name, mask, "m", scale, offset, name,"T",h_discretisation,v_discretisation,l_discretisation);
  switch (l_discretisation) {
    case LGP0:
      variable.att[4].data=(char *) poc_strdup("time z-LGP0 lat lon");
      break;
    case LGP1:
      variable.att[4].data=(char *) poc_strdup("time z-LGP1 lat lon");
      break;
    case DGP1:
      variable.att[4].data=(char *) poc_strdup("time z-DGP1 lat lon");
      break;
    case NCP1:
      variable.att[4].data=(char *) poc_strdup("time z-NCP1 lat lon");
      break;
    case DNP1:
      variable.att[4].data=(char *) poc_strdup("time z-DNP1 lat lon");
      break;
    case LGP2:
      variable.att[4].data=(char *) poc_strdup("time z-LGP2 lat lon");
      break;
      }

/**----------------------------------------------------------------------------
  check if variable already exists */
  variable.id= cdf_identify((const char *) filename, name);
  if(variable.id==-1) {
    status = cdf_createvariable((char *) filename, &(variable));
    }
  variable_id = variable.id;

  variable.destroy();
  
  status = poc_put_UG4D((char *) filename, mesh, count, variable_id, buffer);

terminate:
  return (0);

error:
  printf("archiving_UGarchive failed ..\n");
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int archiving_UGdummy3D(const char *localname, mesh_t & mesh, const char *name, double **buffer, int h_discretisation, int v_discretisation, int l_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int k,n,nnodes,nlevels;
  double **tmp;
  
  discretisation_t descriptor=get_descriptor(mesh,h_discretisation);
  nnodes=descriptor.nnodes;
  
  switch (v_discretisation) {
    case LAYERS:
      nlevels=mesh.nlayers;
      break;
    case LEVELS:
      nlevels=mesh.nlevels;
      break;
    default:
      check_error(-1, "discretisation not implemented in archive procedure", __LINE__, __FILE__, 1);
      break;
    }
/**----------------------------------------------------------------------------
  re-organize buffer order to save layer slices */
  tmp=new double*[nlevels];

  for(k=0;k<nlevels;k++) {
    tmp[k]=new double[nnodes];
    for (n=0;n<nnodes;n++) {
      tmp[k][n]=buffer[n][k];
      }
    }
  status=archiving_UGdummy3D(localname, mesh, name, tmp, (int) 0, h_discretisation, v_discretisation, l_discretisation);
  for(k=0;k<nlevels;k++) {
     if (tmp[k] != NULL) {
         delete[] tmp[k];
         tmp[k]=NULL;
         }
     }
  if (tmp != NULL){
   delete[] tmp;
   tmp=NULL;
  }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int archiving_UGdummy3D(const char *localname, mesh_t & mesh, const char *name1,const char *name2, complex <double> **buffer, int h_discretisation, int v_discretisation, int l_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int k,n,nnodes,nlevels;
  double **tmp;

  switch (h_discretisation) {
    case LGP2:
      nnodes=mesh.LGP2descriptor.nnodes;
      break;

    case LGP1:
      nnodes=mesh.LGP1descriptor.nnodes;
      break;

    case DGP1:
      nnodes=mesh.DGP1descriptor.nnodes;
      break;

    case LGP0:
      nnodes=mesh.LGP0descriptor.nnodes;
      break;

    case NCP1:
      nnodes=mesh.nedges;
      break;

    default:
      check_error(-1, "discretisation not implemented in archive procedure", __LINE__, __FILE__, 1);
      break;
    }
  switch (v_discretisation) {
    case LAYERS:
      nlevels=mesh.nlayers;
      break;
    case LEVELS:
      nlevels=mesh.nlevels;
      break;
    default:
      check_error(-1, "discretisation not implemented in archive procedure", __LINE__, __FILE__, 1);
      break;
    }

/**----------------------------------------------------------------------------
  re-arrange by layers */
  tmp=new double*[nlevels];
  for(k=0;k<nlevels;k++) {
    tmp[k]=new double[nnodes];
    for (n=0;n<nnodes;n++) {
      tmp[k][n]=abs(buffer[n][k]);
      }
    }
  status=archiving_UGdummy3D(localname, mesh, name1, tmp, (int) 0, h_discretisation, v_discretisation, l_discretisation);

  for(k=0;k<nlevels;k++) {
//    tmp[n]=new double[nnodes];
    for (n=0;n<nnodes;n++) {
      tmp[k][n]=-arg(buffer[n][k])*180./M_PI;
      }
    }
  status=archiving_UGdummy3D(localname, mesh, name2, tmp, (int) 0, h_discretisation, v_discretisation, l_discretisation);

  for(k=0;k<nlevels;k++) {
    delete[] tmp[k];
    }
  delete[] tmp;
  return (status);
}
