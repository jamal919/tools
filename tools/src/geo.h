
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief projection, distance and wrap-around function delarations
*/
/*----------------------------------------------------------------------------*/

#if GEO_H == 0
#define GEO_H 1

#include <proj_api.h> /* for projPJ */

#include "functions.h" /* for toponyms_t */


typedef struct {
  int    type;
  double scale;
  double lon,lat;
  double radius;
  double unit;
//  double pole_t,pole_p;
  double x_offset,y_offset;
  projPJ reference;
} geo_t;


#define NORTH_ATLANTIC     11
#define HUDSON_BAY              111
#define HUDSON_STRAIT           1111
#define HUDSON_FOX              1112
#define HUDSON_SOUTH            1113
#define NORTH_EAST_ATLANTIC     112

#define TROPICAL_ATLANTIC  12
#define SOUTH_ATLANTIC     13

#define MEDITERRANEAN      20

#define TROPICAL_INDIAN    32
#define SOUTH_INDIAN       33

#define NORTH_PACIFIC      41
#define YELLOW_SEA         411
#define OKHOTSK_SEA        412
#define BERING_SEA         413
#define JAPAN_SEA          414

#define TROPICAL_PACIFIC   42
#define CHINA_SEA          421
#define TONKIN_GULF        4211

#define SOUTH_PACIFIC      43

#define ARCTIC             50

extern int geo_mercator_directe(const geo_t & projection, double lon,double lat,double *x,double *y);
extern int geo_mercator_inverse(const geo_t & projection, double *lon,double *lat,double x,double y);
extern geo_t geo_mercator_init(double lon, double lat, double radius);

extern double geo_angle(double t0,double p0,double t1,double p1,double t2,double p2);

extern double geo_distance(double t1,double p1,double t2,double p2, char unit='m');
extern  float geo_distance(float  t1,float  p1,float  t2,float  p2, char unit='m');

extern double geo_km(double t1,double p1,double t2,double p2);
extern  float geo_km(float t1,float p1,float t2,float p2);

extern float  geo_haversin(const float& t1, const float& p1, const float& t2, const float& p2);
extern double geo_haversin(const double& t1,const double& p1,const double& t2,const double& p2);
extern float  geo_haversin_km(const float& t1, const float& p1, const float& t2, const float& p2);
extern double geo_haversin_km(const double& t1,const double& p1,const double& t2,const double& p2);

extern double geo_haversinR(const double& t1,const double& p1,const double& t2,const double& p2);
extern double geo_haversinR_km(const double& t1,const double& p1,const double& t2,const double& p2);

extern int geo_init_regions(toponyms_t *toponyms);

extern double gravitational_constant_01(double s,double c);
extern double gravitational_constant_01rad(double lat);
extern double gravitational_constant_01deg(double lat);


extern float  geo_recale(float  value, float  center, float  range=M_PI);
extern double geo_recale(double value, double center, double range=M_PI);
extern float  degree_recale(float value, float center);
extern double degree_recale(double value, double center);
extern void degree_recale(double *value, double center);
extern double radian_recale(double value, double center=0.);
extern void radian_recale(double *value, double center=0.);

#define GEO_STR_PROJ_NONE  "NONE"
#define GEO_STR_PROJ_MERC  "Mercator"
#define GEO_STR_PROJ_STEN  "Stereo nord"
#define GEO_STR_PROJ_STES  "Stereo sud"
#define GEO_STR_PROJ_LAMB  "Lambert conforme"
#define GEO_STR_PROJ_PTAN  "Plan tangent"
#define GEO_STR_PROJ_CTAN  "Cone tangent"
#define GEO_STR_PROJ_GARL  "Garlic"
#define GEO_STR_PROJ_TRAN  "Transverse Mercator"
#define GEO_STR_PROJ_MERC2 "Mercator  2"

#define GEO_PROJECTION_NONE  0
#define GEO_PROJECTION_MERC  1
#define GEO_PROJECTION_STEN  2
#define GEO_PROJECTION_STES  3
#define GEO_PROJECTION_LAMB  4
#define GEO_PROJECTION_PTAN  5
#define GEO_PROJECTION_CTAN  6
#define GEO_PROJECTION_GARL  7
#define GEO_PROJECTION_TRAN  8
#define GEO_PROJECTION_MERC2 9

#define   PROJ_TYPE_NONE   -1
#define   PROJ_TYPE_MERCA   0
#define   PROJ_TYPE_CC      1
#define   PROJ_TYPE_MILL    2
#define   PROJ_TYPE_EQC     3
#define   PROJ_TYPE_CASS    4
#define   PROJ_TYPE_SINU    5
#define   PROJ_TYPE_MOLL    6
#define   PROJ_TYPE_ROBIN   7
#define   PROJ_TYPE_COLLG   8
#define   PROJ_TYPE_LCC     9
#define   PROJ_TYPE_PCONIC 10
#define   PROJ_TYPE_STN    11
#define   PROJ_TYPE_STS    12
#define   PROJ_TYPE_ORTHO  13
#define   PROJ_TYPE_AUGUST 14
#define   PROJ_TYPE_LAGRAN 15

#define   PROJ_STRG_NONE   "NONE"
#define   PROJ_STRG_MERCA  "MERCATOR"
#define   PROJ_STRG_CC     "UNNAMED"
#define   PROJ_STRG_MILL   "UNNAMED"
#define   PROJ_STRG_EQC    "UNNAMED"
#define   PROJ_STRG_CASS   "UNNAMED"
#define   PROJ_STRG_SINU   "UNNAMED"
#define   PROJ_STRG_MOLL   "UNNAMED"
#define   PROJ_STRG_ROBIN  "UNNAMED"
#define   PROJ_STRG_COLLG  "UNNAMED"
#define   PROJ_STRG_LCC    "UNNAMED"
#define   PROJ_STRG_PCONIC "UNNAMED"
#define   PROJ_STRG_STN    "UNNAMED"
#define   PROJ_STRG_STS    "UNNAMED"
#define   PROJ_STRG_ORTHO  "UNNAMED"
#define   PROJ_STRG_AUGUST "UNNAMED"
#define   PROJ_STRG_LAGRAN "UNNAMED"

extern projPJ init_projection(const char *proj4_options,bool doStrstr,int verbose=1);
extern projPJ *init_projection_parallel(const char *proj4_options,int ncpu,bool doStrstr,int verbose=1);
extern void free_threadSafe_projection(projPJ proj);
extern void pj_fwd(projPJ proj,double lon,double lat,double *x,double *y);

extern void geo_to_projection(projPJ ref, double lat,double lon,double *x,double *y);
extern  int geo_to_projection(const char *proj4_options, double *x, double *y, int  ndata, projPJ *projection=0);
extern  int geo_to_projection(const char *proj4_options, double *lon, double *lat, double *x, double *y, int  ndata);

extern void pj_inv(projPJ proj,double x,double y,double *lon,double *lat);

extern int projection_to_geo(projPJ  ref, double *lat,double *lon,double x,double y);
extern int projection_to_geo(projPJ proj, double *x, double *y, int  ndata);

extern int projection_to_geo(const char *proj4_options, double *x, double *y, int  ndata, int verbose=0);
extern int projection_to_geo(const char *proj4_options, double *x, double *y, double *t, double *p, int  ndata);

extern projPJ assign_projection(int pclass, double *lat, double *lon, char *pj_parameters=0);

extern projPJ assign_StereoOblique(double, double, char* pj_parameters=0);
extern projPJ assign_StereoNorth(double, double, double, char* pj_parameters=0);

#include "maths.h"
#include "constants.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  inline double geo_angle_radian(double t0,double p0,double t1,double p1,double t2,double p2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  geo angle in radians from lon and lat in degrees

  t0 longitude in degrees
  p0 latitude  in degrees
  t1 longitude in degrees
  p1 latitude  in degrees
  t1 longitude in degrees
  p1 latitude  in degrees

  return angle in radians
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double alpha,scalar;
  vector3_t o,u,v,uv;
  
/*------------------------------------------------------------------------------
  returns unit vectors */
  o=math_polar2cartesian(t0*d2r,p0*d2r);
  u=math_polar2cartesian(t1*d2r,p1*d2r);
  v=math_polar2cartesian(t2*d2r,p2*d2r);
 
/*------------------------------------------------------------------------------
  operator * is scalar product */
  u-=o*(o*u);
  v-=o*(o*v);
  
  scalar=u*v;
  uv=u&v;
  
  alpha=atan2(hypot(uv),scalar);
  if(uv*o<0)
    alpha*=-1.;
  
  return alpha;
}


#endif
