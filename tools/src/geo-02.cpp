
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief proj lib wrappers
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <string.h>

#include "constants.h"
#include "tools-structures.h"

#include "geo.h"
#include "functions.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  string check_proj4_options(const char *proj4_options,bool doStrstr,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const char *s, *r, *chk;
  string tmp=proj4_options;
  bool ok=true;
  
/*------------------------------------------------------------------------------
  try to filter out old (and non-standrad) way of prescribing prj4 options. 
  It is now required to build options in a strict proj4 fashion.
------------------------------------------------------------------------------*/
  
  if(doStrstr) {
    s=strstr(proj4_options,"--proj=");
    if(s!=0) {
      printf("invoking command did not properly handled --proj option (%s), please contact sirocco developpers\n",proj4_options);
      TRAP_ERR_EXIT(ENOEXEC,"*** libproj PARAMETER UPDATE HAS EXPIRED ON %d ***\n");
      }
    s=strstr(proj4_options,"+proj=");
    r=strstr(proj4_options,"proj=");
    if(s==0 and r!=0) {
      printf("please provide +proj= instead of proj= next time (%s)\n",proj4_options);
      tmp="+"+tmp;
      s=strdup(tmp.c_str());
      }
    else if(s==0 and r==0) {
      printf("invoking command did not set +proj (%s) \n",proj4_options);
      TRAP_ERR_EXIT(ENOEXEC,"\n*** improper proj4 options, abort ***\n\n");
      }
    }
  else
    s=proj4_options;
  
  int n;
  string str=s;
  
  /* remove double spaces */
  do{
    n=replace(&str,"  "," ");
    }while(n>0);
  
  bool correctedSomething=false;
  
  /* replace all " " by " +" */
  for(n=0;;n+=2){
    
    /* find " " but not " +" */
    while(1){
      n=str.find(' ',n);
      if(n<0) break;
      if(str.substr(n,2)!=" +")break;
      n++;
      }
    
    if(n<0) break;
    
    str.replace(n,1," +");
    correctedSomething=true;
    }
  if(str.substr(0,1)!="+"){
    str.insert(0,"+");
    correctedSomething=true;
    }
  
  if(correctedSomething){
#undef REMOVAL_DATE
#define REMOVAL_DATE 20150125
    if(expire(REMOVAL_DATE))
    TRAP_ERR_EXIT(ENOEXEC,"\n*** improper proj4 options, missing + key (found %s, expected %s), abort ***\n\n",s,str.c_str());
    }
  
  if(verbose)
    printf("argument passed to libproj initialisation : "+str+"\n");
  
  return str;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  projPJ init_projection(const char *proj4_options, bool doStrstr, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  projPJ ref;
  string str=check_proj4_options(proj4_options,doStrstr,verbose);
  
  ref = pj_init_plus(str.c_str());
  if (!ref) __TRAP_ERR_EXIT__(1,"Projection initialization failed\n");
  
  return ref;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  projPJ *init_projection_parallel(const char *proj4_options,int ncpu,bool doStrstr,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  initialise projection buffers for OPENMP optimisation
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  string str=check_proj4_options(proj4_options,doStrstr,verbose);
  
  projPJ *refs=new projPJ[ncpu];
  
  if(verbose==1) printf("%s : initialise projection with paramaters=%s\n",__func__,str.c_str());
  
  for(int k=0;k<ncpu;k++){
    projPJ &ref=refs[k];
    ref=pj_init_plus_ctx(pj_ctx_alloc(),str.c_str());
//     if (!ref && verbose) __TRAP_ERR_EXIT__(1,"Projection initialization failed\n"); /// HERE !!!
    if (!ref) __TRAP_ERR_EXIT__(1,"Projection initialization failed\n");
    }
  
  return refs;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void free_threadSafe_projection(projPJ proj)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  pj_ctx_free(pj_get_ctx(proj));
  pj_free(proj);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void pj_fwd(projPJ proj,double lon,double lat,double *x,double *y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///overload
/**
\param lon longitude in degrees
\param lat latitude in degrees
*/
/*----------------------------------------------------------------------------*/
{
  projUV idata;
  projUV odata;
//   LP lp;
//   struct FACTORS *fac;
//   double h;
//   int status;

  idata.u=lon*DEG_TO_RAD;
  idata.v=lat*DEG_TO_RAD;

//  status= proj_factors( lp, ref, h, fac);

  odata=pj_fwd(idata,proj);
  
  *x=odata.u;
  *y=odata.v;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void geo_to_projection(projPJ ref, double lat,double lon,double *x,double *y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(ref==0) {
    *x=lon;
    *y=lat;
    }
  else {
    pj_fwd(ref,lon,lat,x,y);
    if( (*x==HUGE_VAL) || (*y==HUGE_VAL) ) {
      *x=-9999;
      *y=-9999;
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int geo_to_projection(const char *proj4_options, double *x, double *y, int  ndata, projPJ *projection)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  double t,p;

  status=0;
  if(projection!=0) *projection=0;

  if(proj4_options!=0) {
    
    int nprocs=initialize_OPENMP(-1, 0);
    
    projPJ *refs = init_projection_parallel(proj4_options,nprocs,true,0);
    if(projection!=0) *projection=refs[0];
#pragma omp parallel for private(t,p) if(nprocs>1)
    for(size_t n=0;n<ndata;n++) {
      t=x[n];
      p=y[n];
      geo_to_projection(refs[omp_get_thread_num()], p, t, &x[n], &y[n]);
      }
    
    deletep2D(&refs,nprocs,free_threadSafe_projection);
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int geo_to_projection(const char *proj4_options, double *lon, double *lat, double *x, double *y, int  ndata)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  double t,p;

  status=0;

  if(proj4_options!=0) {
    
    int nprocs=initialize_OPENMP(-1, 0);
    
    projPJ *refs = init_projection_parallel(proj4_options,nprocs,true,0);
    
#pragma omp parallel for private(t,p) if(nprocs>1)
    for(size_t n=0;n<ndata;n++) {
      t=lon[n];
      p=lat[n];
      geo_to_projection(refs[omp_get_thread_num()], p, t, &x[n], &y[n]);
      }
    
    deletep2D(&refs,nprocs,free_threadSafe_projection);
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void pj_inv(projPJ proj,double x,double y,double *lon,double *lat)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///overload
/**
\param *lon longitude in degrees
\param *lat latitude in degrees
*/
/*----------------------------------------------------------------------------*/
{
  projUV odata;
  projUV idata;

  idata.u=x;
  idata.v=y;

  odata=pj_inv(idata,proj);
  
  *lon=odata.u*r2d;
  *lat=odata.v*r2d;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_to_geo(projPJ ref, double *lat,double *lon,double x,double y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  
  if(ref==0) {
    *lon=x;
    *lat=y;
    }
  else {
    pj_inv(ref,x,y,lon,lat);
    if( (*lon==HUGE_VAL) || (*lat==HUGE_VAL) ) {
      *lon=-9999;
      *lat=-9999;
      status=-1;
      }
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_to_geo(const char *proj4_options, double *x, double *y, int  ndata, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  double t,p;

  status=0;

  if(proj4_options!=0) {
    
    int nprocs=initialize_OPENMP(-1, 0);
    
    projPJ *refs = init_projection_parallel(proj4_options,nprocs,true,verbose);
 
#ifdef EMODNET2015 // yet another shameful EMODNET distribution format
    double xx=0,yy=0;
    projection_to_geo(refs[0],&p,&t,xx,yy);
    
    xx=20;
    yy=20;
    projection_to_geo(refs[0],&p,&t,xx,yy);
#pragma omp parallel for private(t,p) if(nprocs>1)
    for(size_t n=0;n<ndata;n++) {
      x[n]+=-16.2520833333333;
      y[n]+=+39.9979166666667;
      }
#else    
#pragma omp parallel for private(t,p) if(nprocs>1)
    for(size_t n=0;n<ndata;n++) {
      projection_to_geo(refs[omp_get_thread_num()],&p,&t,x[n],y[n]);
      x[n]=t;
      y[n]=p;
      }
#endif

    deletep2D(&refs,nprocs,free_threadSafe_projection);
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_to_geo(projPJ proj, double *x, double *y, int  ndata)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int chk, status;
  double t,p;

  status=0;

  if(proj!=0) {
    
    for(size_t n=0;n<ndata;n++) {
      chk=projection_to_geo(proj, &p, &t, x[n], y[n]);
      x[n]=t;
      y[n]=p;
      if(chk==-1) status=-1;
      }
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_to_geo (const char *proj4_options, double *x, double *y, double *t, double *p, int  ndata)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=0;

  if(proj4_options!=0) {
    
    int nprocs=initialize_OPENMP(-1, 0);
    
    projPJ *refs = init_projection_parallel(proj4_options,nprocs,true,0);
    
#pragma omp parallel if(nprocs>1)
    for(size_t n=0;n<ndata;n++) {
      const double xx=x[n];
      const double yy=y[n];
      projection_to_geo(refs[omp_get_thread_num()],&(p[n]),&(t[n]),xx,yy);
//        if(t[n]==0) {
//          printf("alert %d %lf %lf %lf %lf %d %d\n",n, p[n],t[n],x[n],y[n], nprocs, omp_get_thread_num());
//          }
      }
    
    deletep2D(&refs,nprocs,free_threadSafe_projection);
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void projection_test()
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  projPJ ref;
  const char *parms = "+proj=merc +ellps=GRS80 +lon_0=0E";

  if ( ! (ref = pj_init_plus(parms)) ) {
    __ERR_BASE_LINE__( "Projection initialization failed\n"); exit(1);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  projPJ assign_StereoOblique(double lat0, double lon0, char *pj_parameters)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  WARNING: this is not stereo oblique!!! (still probably relevant in most cases)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  projPJ ref;
//  char *ellipsoid={"ellps=GRS80"};
  const char *ellipsoid="ellps=WGS84";
  char *parms;

//  char *parms[] = { "proj=utm", "ellps=WGS84", "zone=3"};
  
  if(lat0==0.) lat0=0.001; /* because of a bug in stere proj */
  
  asprintf(&parms,"+proj=stere +%s +lat_0=%lf +lon_0=%lf",ellipsoid,lat0,lon0);
  ref = pj_init_plus(parms);
  
  if(pj_parameters!=0) strcpy(pj_parameters,parms);
  
  free(parms);

  return(ref);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  projPJ assign_StereoNorth(double lat0, double lat_ts, double lon0, char *pj_parameters)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  projPJ ref;
//  char *ellipsoid={"ellps=GRS80"};
  const char *ellipsoid="ellps=WGS84";
  char *parms;

  asprintf(&parms,"+proj=stere +%s +lat_ts=%lf +lat_0=%lf +lon_0=%lf",ellipsoid,lat0,lat_ts,lon0);
  ref = pj_init_plus(parms);
  
  if(pj_parameters!=0) strcpy(pj_parameters,parms);
  
  free(parms);

  return(ref);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  projPJ assign_projection(int pclass, double *lat, double *lon, char *pj_parameters)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  projPJ ref;

//  char *ellipsoid={"ellps=GRS80"};
  const char *ellipsoid="ellps=WGS84";
  char *parms;

  switch(pclass) {
    case 0/**/ :
      asprintf(&parms,"+proj=merc +%s +lon_0=%lfE +lat_0=%lfE +to_meter=%lfE",ellipsoid,lon[0],lat[0],1./cos(lat[0]+M_PI/180.0));
      break;
    case 1/**/ :
      asprintf(&parms,"+proj=cc +%s +lon_0=0E",ellipsoid);
      break;
    case 2/**/ :
      asprintf(&parms,"+proj=mill +%s +lon_0=0E",ellipsoid);
      break;
    case 3/**/ :
      asprintf(&parms,"+proj=eqc +%s +lon_0=0E",ellipsoid);
      break;
    case 4/*PROJ_CASS*/ :
      asprintf(&parms,"+proj=cass +%s +lon_0=0E",ellipsoid);
      break;
    case 5/*PROJ_SINU*/ :
      asprintf(&parms,"+proj=sinu +%s +lon_0=0E",ellipsoid);
      break;
    case 6/*PROJ_MOLL*/ :
      asprintf(&parms,"+proj=moll +%s +lon_0=0E",ellipsoid);
      break;
    case 7/*PROJ_ROBIN*/ :
      asprintf(&parms,"+proj=robin +%s +lon_0=0E",ellipsoid);
      break;
    case 8/*PROJ_COLLG*/ :
      asprintf(&parms,"+proj=collg +%s +lon_0=0E",ellipsoid);
      break;
    case 9/*PROJ_LCC*/ :
      asprintf(&parms,"+proj=lcc +%s +lon_0=0E",ellipsoid);
      break;
    case 10/*PROJ_PCONIC*/ :
      asprintf(&parms,"+proj=pconic +%s +lon_0=0E",ellipsoid);
      break;
    case 11/*PROJ_STN*/ :
      asprintf(&parms,"+proj=ortho +%s +lon_0=0E","lat_0=-90",ellipsoid);
      break;
    case 12/*PROJ_STS*/ :
      asprintf(&parms,"+proj=ortho +%s +lon_0=0E","lat_0=90",ellipsoid);
      break;
    case 13/*PROJ_ORTHO*/ :
      asprintf(&parms,"+proj=ortho +%s +lon_0=0E",ellipsoid);
      break;
    case 14/*PROJ_AUGUST*/ :
      asprintf(&parms,"+proj=august +%s +lon_0=0E",ellipsoid);
      break;
    case 15/*PROJ_LAGRAN*/ :
      asprintf(&parms,"+proj=lagrng +%s +lon_0=0E",ellipsoid);
      break;
    }
  ref = pj_init_plus(parms);
  if(pj_parameters!=0) strcpy(pj_parameters,parms);
  free(parms);

  return(ref);
}
