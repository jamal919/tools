
/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

*******************************************************************************/

#include "maths.h" // for cartesian<->polar
#include "constants.h"
#include "geo.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> T geo_recale_template(T value, T center, T range)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///ensures value is within center +/- range
/*----------------------------------------------------------------------------*/
{
  T x=value;
  const T _2range=range*2;
  
  if(not isfinite(value) or not isfinite(range))
    return value;
  
  if(not isfinite(center)){
    if(isnan(center))
      return value;
    else
      return center;
    }
  
  if(value==center)
    return value;
  
  if(center+range==center || center-range==center)
    TRAP_ERR_RETURN(NAN,1,"infinitely small range: %s(%g,%g,%g)\n",__func__,value,center,range);
  
  while (x > center+range){
    x-=_2range;
    if(x==value)
      TRAP_ERR_RETURN(NAN,1,"infinite loop: %s(%.9g,%g,%g)\n",__func__,value,center,range);
    }
  while (x <  center-range){
    x+=_2range;
    if(x==value)
      TRAP_ERR_RETURN(NAN,1,"infinite loop: %s(%.9g,%g,%g)\n",__func__,value,center,range);
    }
  
  return x;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float geo_recale(float value, float center, float range)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float x;
  
  x=geo_recale_template(value,center,range);
  
  return x;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double geo_recale(double value, double center, double range)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double x;
  
  x=geo_recale_template(value,center,range);
  
  return x;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double degree_recale(double value, double center)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double x;
  
  x=geo_recale/*avoid perlReplace.pl*/(value,center,180.);
  
  return x;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void degree_recale(double *value, double center)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  *value=degree_recale(*value,center);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float degree_recale(float value, float center)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float x;
  
  x=geo_recale/*avoid perlReplace.pl*/(value,center,180.f);
  
  return x;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double radian_recale(double value, double center)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double x;
  
  x=geo_recale/*avoid perlReplace.pl*/(value,center,M_PI);
  
  return x;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void radian_recale(double *value, double center)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  *value=geo_recale/*avoid perlReplace.pl*/(*value,center,M_PI);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int geo_mercator_directe(const geo_t & projection, double lon,double lat,double *x,double *y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double t,c1,c2;

//   lon+=projection.pole_t;
//   lat+=90.-projection.pole_p;

  t=tan(0.5*lat*d2r);

  c1=(lon-projection.lon)*d2r;
  c2=log((1+t)/(1-t));

  *x=projection.scale*c1-projection.x_offset;
  *y=projection.scale*c2-projection.y_offset;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int geo_mercator_inverse(const geo_t & projection, double *lon,double *lat,double x,double y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double c1,c2;

  c1=(x+projection.x_offset)/projection.scale;
  c2=exp((y+projection.y_offset)/projection.scale);

  *lon=c1/d2r+projection.lon;
  *lat=atan2(c2-1,c2+1)*2.0/d2r;

//   *lon-=projection.pole_t;
//   *lat-=90.-projection.pole_p;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

geo_t geo_mercator_init(double lon,double lat, double radius)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double x,y;
  geo_t projection;
  int status;

  projection.radius=radius;
  projection.scale=cos(lat*d2r)*radius;

  projection.x_offset=0;
  projection.y_offset=0;

  projection.lon=lon;
  projection.lat=lat;

//   projection.pole_t=pole_t;
//   projection.pole_p=pole_p;

  status=geo_mercator_directe(projection,  lon, lat, &x, &y);

  projection.x_offset=x;
  projection.y_offset=y;

  return(projection);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  geo_t geo_init(int ptype)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  geo_t projection;

  switch (ptype) {
    case GEO_PROJECTION_MERC:
//      geo_mercator_init(double lon,double lat, double radius)
      break;
    }
  return(projection);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double geo_angle(double t0,double p0,double t1,double p1,double t2,double p2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  geo angle in degrees from lon and lat in degrees

  return angle in degrees
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double alpha;
  
  alpha=geo_angle_radian(t0, p0, t1, p1, t2, p2);
  alpha*=r2d;
  
  return alpha;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double geo_distance_d2r(double t1,double p1,double t2,double p2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// geo distance in radian from lon and lat in degrees
/**
\param t1 longitude in degrees
\param p1 latitude in degrees
\param t1 longitude in degrees
\param p1 latitude in degrees
\return angle in radians
*/
/*----------------------------------------------------------------------------*/
{
  double dt,dp;
  double a1,a2,b1,b2,pds,angle;
  
  dt=t2-t1;
  if(dt== 360.0) dt=0.0;
  if(dt==-360.0) dt=0.0;
  dp=p2-p1;

  if ((dt == 0.0) && (dp == 0.0)) return(0.);

  a1=t1*d2r;
  b1=p1*d2r;
  a2=t2*d2r;
  b2=p2*d2r;
  
  vector3_t u,v;
  u=math_polar2cartesian(a1,b1);
  v=math_polar2cartesian(a2,b2);
  
  pds=u*v;
  
  if(pds >= 1.)
    return 0.;
  
  angle=acos(pds);
  
  return angle;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> T geo_distance_template(T t1, T p1, T t2, T p2, char unit)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// returns spherical arc length in km from lon and lat in degrees of 2 positions
/**
\param t1 1st position longitude in degrees
\param p1 1st position latitude in degrees
\param t1 2nd position longitude in degrees
\param p1 2nd position latitude in degrees
\param unit m for meters (default); k for km; d for degrees
\return distance in requested unit, meters by default
*/
/*----------------------------------------------------------------------------*/
{
  T angle,ro;
  
  angle=geo_distance_d2r(t1,p1,t2,p2);
  
  if(angle==0.)
    return 0.;
  
  switch(unit){
  case 'm':
    ro=angle*MeanEarthRadius;
    break;
  case 'k':
    ro=angle*MeanEarthRadius*1e-3;
    break;
  case 'd':
    ro=angle*r2d;
    break;
  case 'r':
    ro=angle;
    break;
  default:
    TRAP_ERR_EXIT(ENOEXEC,"programming error: unit '%c' not recognised\n",unit);
    }
  
  return ro;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double geo_distance(double t1,double p1,double t2,double p2, char unit)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=geo_distance_template(t1,p1,t2,p2,unit);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float geo_distance(float t1,float p1,float t2,float p2, char unit)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float result;
  
  result=geo_distance_template(t1,p1,t2,p2,unit);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double geo_km(double t1,double p1,double t2,double p2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=geo_distance_template(t1,p1,t2,p2,'k');
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float geo_km(float t1,float p1,float t2,float p2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=geo_distance_template(t1,p1,t2,p2,'k');
  
  return result;
}


//     /// @brief The usual PI/180 constant
//     static const double DEG_TO_RAD = 0.017453292519943295769236907684886;
//     /// @brief Earth's quatratic mean radius for WGS-84
//     static const double EARTH_RADIUS_IN_METERS = 6372797.560856;
//
//     /** @brief Computes the arc, in radian, between two WGS-84 positions.
//       *
//       * The result is equal to <code>Distance(from,to)/EARTH_RADIUS_IN_METERS</code>
//       *    <code>= 2*asin(sqrt(h(d/EARTH_RADIUS_IN_METERS )))</code>
//       *
//       * where:<ul>
//       *    <li>d is the distance in meters between 'from' and 'to' positions.</li>
//       *    <li>h is the haversine function: <code>h(x)=sin²(x/2)</code></li>
//       * </ul>
//       *
//       * The haversine formula gives:
//       *    <code>h(d/R) = h(from.lat-to.lat)+h(from.lon-to.lon)+cos(from.lat)*cos(to.lat)</code>
//       *
//       * @sa http://en.wikipedia.org/wiki/Law_of_haversines
//       */
//     double ArcInRadians(const Position& from, const Position& to) {
//         double latitudeArc  = (from.lat - to.lat) * DEG_TO_RAD;
//         double longitudeArc = (from.lon - to.lon) * DEG_TO_RAD;
//         double latitudeH = sin(latitudeArc * 0.5);
//         latitudeH *= latitudeH;
//         double lontitudeH = sin(longitudeArc * 0.5);
//         lontitudeH *= lontitudeH;
//         double tmp = cos(from.lat*DEG_TO_RAD) * cos(to.lat*DEG_TO_RAD);
//         return 2.0 * asin(sqrt(latitudeH + tmp*lontitudeH));
//     }
//
//     /** @brief Computes the distance, in meters, between two WGS-84 positions.
//       *
//       * The result is equal to <code>EARTH_RADIUS_IN_METERS*ArcInRadians(from,to)</code>
//       *
//       * @sa ArcInRadians
//       */
//     double DistanceInMeters(const Position& from, const Position& to) {
//         return EARTH_RADIUS_IN_METERS*ArcInRadians(from, to);
//     }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> inline T geo_haversinR_template(const T& t1,const T& p1,const T& t2,const T& p2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T dt,dp;
  T c1,c2,s1,s2,pds,angle;
  T ro;
/// http://en.wikipedia.org/wiki/Haversine_formula
/// d = 2 r arcsin{ sqrt[ sin² (0.5*(lat2-lat1) ) + cos(lat1)cos(lat2)sin²(0.5*(lon2-lon1)) ] }

  dt=t2-t1;
  dp=p2-p1;

  if ((dt == 0.0) && (dp == 0.0)) return(0.);

  c1=cos(p1);
  c2=cos(p2);

  s1=sin(dt/2.0);
  s2=sin(dp/2.0);

  pds=s2*s2+c1*c2*s1*s1;

  if(pds < 1.0) {
    angle=2.0*asin(sqrt(pds));
    ro=(T)(MeanEarthRadius*angle);
    return(ro);
    }
  else return(0.);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double geo_haversinR_km(const double& t1,const double& p1,const double& t2,const double& p2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  return(geo_haversinR_template(t1,p1,t2,p2)/1000.0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double geo_haversinR(const double& t1,const double& p1,const double& t2,const double& p2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  return(geo_haversinR_template(t1,p1,t2,p2));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> inline T geo_haversin_template(const T& t1,const T& p1,const T& t2,const T& p2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  http://en.wikipedia.org/wiki/Haversine_formula
  
  d = 2 r arcsin{ sqrt[ sin² (0.5*(lat2-lat1) ) + cos(lat1)cos(lat2)sin²(0.5*(lon2-lon1)) ] }

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  T dt,dp;
  T c1,c2,s1,s2,pds,angle;
  T ro;

  dt=t2-t1;
  dp=p2-p1;

  if ((dt == 0.0) && (dp == 0.0)) return(0.);

  c1=cos(DEG_TO_RAD*p1);
  c2=cos(DEG_TO_RAD*p2);

  s1=sin(DEG_TO_RAD*dt/2.0);
  s2=sin(DEG_TO_RAD*dp/2.0);

  pds=s2*s2+c1*c2*s1*s1;

  if(pds <= 1.0) {
    angle=2.0*asin(sqrt(pds));
    ro=(T)(MeanEarthRadius*angle);
    return(ro);
    }
  else {
    angle=M_PI;
    ro=(T)(MeanEarthRadius*angle);
    return(ro);
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float geo_haversin_km(const float& t1, const float& p1, const float& t2, const float& p2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float ro=geo_haversin_template(t1,p1,t2,p2)/1000.0;
  return(ro);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double geo_haversin_km(const double& t1,const double& p1,const double& t2,const double& p2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double ro=geo_haversin_template(t1,p1,t2,p2)/1000.0;
  return(ro);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float geo_haversin(const float& t1, const float& p1, const float& t2, const float& p2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float ro=geo_haversin_template(t1,p1,t2,p2);
  return(ro);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double geo_haversin(const double& t1,const double& p1,const double& t2,const double& p2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double ro=geo_haversin_template(t1,p1,t2,p2);
  return(ro);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int geo_GreatCircle(double t1, double p1, double t2, double p2, double lamda, double &t3, double &p3)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  matrix3x3_t R;
  vector3_t axis,u,v,w;
  double radius=1,dum,angle,chk;

  u=math_polar2cartesian(t1*d2r,p1*d2r,radius);
  v=math_polar2cartesian(t2*d2r,p2*d2r,radius);
  
  chk=hypot(u);
  chk=hypot(v);

  angle=geo_haversinR(t1,p1,t2,p2);

  axis=u & v;
  
  chk=hypot(axis);
  
  axis/=chk;

  R=math_rotation3D_init( axis, angle*lamda);

  w=math_rotation3D('D', R,  u);
  chk=hypot(w);
  
  status=math_cartesian2polar(w,&t3,&p3,&dum);
  
  t3/=d2r;
  p3/=d2r;
  
  return (0.);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int geo_init_regions(toponyms_t *toponyms)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  (*toponyms)["NORTH_ATLANTIC"]   =NORTH_ATLANTIC;
  (*toponyms)["TROPICAL_ATLANTIC"]=TROPICAL_ATLANTIC;
  (*toponyms)["SOUTH_ATLANTIC"]   =SOUTH_ATLANTIC;
  
  (*toponyms)["HUDSON_BAY"]       =HUDSON_BAY;
  (*toponyms)["HUDSON_STRAIT"]    =HUDSON_STRAIT;
  (*toponyms)["HUDSON_FOX"]       =HUDSON_FOX;
  (*toponyms)["HUDSON_SOUTH"]     =HUDSON_SOUTH;
  
  (*toponyms)["NORTH_PACIFIC"]    =NORTH_PACIFIC;
  (*toponyms)["TROPICAL_PACIFIC"] =TROPICAL_PACIFIC;
  (*toponyms)["SOUTH_PACIFIC"]    =SOUTH_PACIFIC;
  
  (*toponyms)["TROPICAL_INDIAN"]  =TROPICAL_INDIAN;
  (*toponyms)["SOUTH_INDIAN"]     =SOUTH_INDIAN;
  
  (*toponyms)["MEDITERRANEAN"]    =MEDITERRANEAN;
  
  (*toponyms)["ARCTIC"]           =ARCTIC;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double gravitational_constant_01(double s,double c)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// latitude-dependant gravity of Earth
/**
\param s sinus of latitude
\param c cosinus of latitude
\returns gravity
*/
/*----------------------------------------------------------------------------*/
{
  double g;
  
  ///from from http://en.wikipedia.org/wiki/Gravity_of_Earth#Mathematical_models
  g=9.780327*(1+0.0053024*s*s-2*0.0000058*s*c*s*c);
  
  return g;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double gravitational_constant_01rad(double lat)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// latitude-dependant gravity of Earth
/**
\param lat latitude in radians
\returns gravity
*/
/*----------------------------------------------------------------------------*/
{
  const double
    s=sin(lat),
    c=cos(lat);
  double g;
  
  ///call gravitational_constant_01(double,double)
  g=gravitational_constant_01(s,c);
  
  return g;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double gravitational_constant_01deg(double lat)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// latitude-dependant gravity of Earth
/**
\param lat latitude in degrees
\returns gravity
*/
/*----------------------------------------------------------------------------*/
{
  const double latr=lat*d2r;
  double g;
  
  ///call gravitational_constant_01rad(double)
  g=gravitational_constant_01rad(latr);
  
  return g;
}
