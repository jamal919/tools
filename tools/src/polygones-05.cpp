
/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

*******************************************************************************/

#include <config.h>

#include <stdio.h>
#include <string>

#include "tools-structures.h"

#include "polygons.hpp"
#include "exceptions.hpp"

#include "polygones.h"
#include "geo.h"
#include "list.h"
#include "maths.h"

#include "functions.h"

using namespace Polygons;// for PolygonSet

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_cartesian_ellipse(double x0, double y0, double a, double b, range_t<double> angles, int npt, plg_t & plg)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n;
  double *x,*y;

  plg.npt=npt;
  
  x=new double[plg.npt];
  y=new double[plg.npt];

  plg.x=new double[plg.npt];
  plg.y=new double[plg.npt];
  
  double delta=angles.width();
  
  for(n=0;n<plg.npt;n++) {
    double alpha=angles.min+n*delta/(plg.npt-1.);
    x[n]=x0+a*cos(alpha);
    y[n]=y0+b*sin(alpha);
    }
  
  for(n=0;n<plg.npt;n++) {
    plg.x[n]=x[n];
    plg.y[n]=y[n];
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_create_ellipse(plg_array_t *plg, double a, double b, double t0, double p0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  int utm_zone;
  double *x,*y,x0,y0;
  char *proj4_options;
  const char *hemisphere;
  projPJ proj;
  
//  L=2.*M_PI*radius;

  plg->n=1;
  plg->p=new plg_t[plg->n];

  plg->p[0].npt=100;
  
  x=new double[plg->p[0].npt];
  y=new double[plg->p[0].npt];

  plg->p[0].t=new double[plg->p[0].npt];
  plg->p[0].p=new double[plg->p[0].npt];
  
  if(p0<0) hemisphere="+south";
  else hemisphere="";
  
  utm_zone=(t0+180)/6+1;
  
  asprintf(&proj4_options,"+proj=utm +zone=%d +ellps=WGS84 +datum=WGS84 +units=m +no_defs %s",utm_zone,hemisphere);

  proj = pj_init_plus(proj4_options);
  if (!proj) TRAP_ERR_EXIT(1,"Projection initialization failed\n");
  
  geo_to_projection(proj, p0, t0, &x0, &y0);
  pj_free(proj);
  
  for(n=0;n<plg->p[0].npt;n++) {
    double alpha=n*2*M_PI/(plg->p[0].npt-1.);
    x[n]=x0+a*cos(alpha);
    y[n]=y0+b*sin(alpha);
    }
  
  status=projection_to_geo (proj4_options, x, y, plg->p[0].npt);
  
  free(proj4_options);

  for(n=0;n<plg->p[0].npt;n++) {
    plg->p[0].t[n]=x[n];
    plg->p[0].p[n]=y[n];
    }

  status=plg_save("ellipse.plg",PLG_FORMAT_SCAN,plg->p,1);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_cartesian_circle(double x0, double y0, double radius, range_t<double> angles, int npt, plg_t & plg)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n;
  double *x,*y;
  double L;
  
  L=2.*M_PI*radius;

  plg.npt=npt;
  
  x=new double[plg.npt];
  y=new double[plg.npt];

  plg.x=new double[plg.npt];
  plg.y=new double[plg.npt];
  
  double delta=angles.width();
  
  for(n=0;n<plg.npt;n++) {
    double alpha=angles.min+n*delta/(plg.npt-1.);
    x[n]=x0+radius*cos(alpha);
    y[n]=y0+radius*sin(alpha);
    }
  
  for(n=0;n<plg.npt;n++) {
    plg.x[n]=x[n];
    plg.y[n]=y[n];
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_cartesian_circle(double t0, double p0, projPJ proj, double radius, range_t<double> angles, int npt, plg_t & plg)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  double *x,*y,x0,y0;
  double L;
  
  L=2.*M_PI*radius;

  plg.npt=npt;
  
  x=new double[plg.npt];
  y=new double[plg.npt];

  plg.t=new double[plg.npt];
  plg.p=new double[plg.npt];
    
  geo_to_projection(proj, p0, t0, &x0, &y0);
  
  double delta=angles.width();
  
  for(n=0;n<plg.npt;n++) {
    double alpha=n*delta/(plg.npt-1.);
    x[n]=x0+radius*cos(alpha);
    y[n]=y0+radius*sin(alpha);
    }
  
  status=projection_to_geo (proj, x, y, plg.npt);
    
  for(n=0;n<plg.npt;n++) {
    plg.t[n]=x[n];
    plg.p[n]=y[n];
    }

//   status=plg_save("circle.plg",PLG_FORMAT_SCAN, &plg,1);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_cartesian_circle(double t0, double p0, projPJ proj, double radius, plg_t & plg)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int npt, status;
  range_t<double> angles;
  
  npt=100;
  
  angles.min=0;
  angles.max=2*M_PI;
    
  status=plg_cartesian_circle(t0, p0, proj, radius, angles, npt, plg);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_circles(const plg_t & p, projPJ proj, double factor, vector<plg_t> & q, vector<double> & sizes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,status;
  double L, radius;
  
  plg_t *q1=new plg_t;
  k=0;
  L=p.cartesian_size(k);
  radius=factor*L;
  status=plg_cartesian_circle(p.t[k], p.p[k], proj, radius, *q1);
  q.push_back(*q1);
  sizes.push_back(L);
  
  plg_t *q2=new plg_t;
  k=p.npt-1;
  L=p.cartesian_size(k-1);
  radius=factor*L;
  status=plg_cartesian_circle(p.t[k], p.p[k], proj, radius, *q2);
  q.push_back(*q2);
  sizes.push_back(L);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_circles(const vector<plg_t> & p, projPJ proj, double factor, vector<plg_t> & q, vector<double> & sizes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  bool debug=false;
  
  q.clear();
  sizes.clear();
  
  for(n=0;n<p.size();n++) {
    status=plg_circles(p[n], proj, factor, q, sizes);
    }

  if(debug) status=plg_save("extremity-circles.plg",PLG_FORMAT_SCAN, q);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_t plg_resample(const plg_t & input, int no)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const int ni=input.npt;
  int i,m=-1;
  vector3_t *v,vi;
  double *d,di;
  plg_t o;
  
  v=new vector3_t[ni];
  d=new double   [ni];
  
  for(i=0;i<ni;i++){
    v[i]=math_polar2cartesian(input.t[i]*d2r,input.p[i]*d2r);
    
    if(i==0)
      di=0.;
    else
      di+=acos(v[i]*v[i-1]);
    
    d[i]=di;
    }
  
  const double k=(no-1)/di;
  
  for(i=0;i<ni;i++){
    d[i]*=k;
    }
  
  o.init(no);
  double *oti,*opi;
  
  for(i=0;i<no;i++){
    map_interpolate1D(v,d,ni,i,&vi,&m);
    
    oti=&o.t[i];
    opi=&o.p[i];
    math_cartesian2polar(vi,oti,opi);
    *oti*=r2d;
    *opi*=r2d;
    }
  
  return o;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_t plg_resample(const plg_t & target, double radius, int equalize, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Resample an existing polygon following a fix radius, possibly tuned
  if equalizing is enabled
  
  resampling along target lines, but pre-existing nodes not preserved
  
  if segment changes, first interception is unique

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  int i,j,n;
  vector2D_t u,v;
  vector2D_t w1,w2,w;
  vector2D_t e1,e2;
  point2D_t p,q;
  point2D_t r,o;
  double a,r1,r2,length,check,misfit,previous_misfit=0.5;
  vector<point2D_t *> extraction;
  point2D_t *dum;
  plg_t out;
  point2D_t first,previous;

  double c=0.75;
  int direction=1;
  
begin:
  first.x=target.x[0];
  first.y=target.y[0];

  length=0;

  r=first;
  previous=r;
  
/*------------------------------------------------------------------------------
  register first point*/
  dum=new point2D_t;
  *dum=r;
  extraction.push_back(dum);
  
  for(j=0;j<target.npt-1;j++) {
    p.x = target.x[j];
    p.y = target.y[j];
    q.x = target.x[j+1];
    q.y = target.y[j+1];
    w1=vector2D_t(r,p);
    w2=vector2D_t(r,q);
/*------------------------------------------------------------------------------
    check if this (new) segment is intercepted */
    r1=hypot(w1);
    r2=hypot(w2);
    if((r1-radius)*(r2-radius)>0) continue;
/*------------------------------------------------------------------------------
    first point on this segment*/
    status=math_orthobase(p, q, &e1, &e2);
    u=math_vector_coordinates01(w2, e1, e2);
    o=q+(-u.x)*e1;
    a=sqrt(radius*radius-u.y*u.y);
    r=o+a*e1;
/*------------------------------------------------------------------------------
    register new point*/
    dum=new point2D_t;
    *dum=r;
    extraction.push_back(dum);

    w=vector2D_t(previous,r);
    check=hypot(w);
    previous=r;
    
    length+=radius;
/*------------------------------------------------------------------------------
    further points on same segment*/
    w=vector2D_t(r,q);
    while(hypot(w)>radius) {
      r=r+radius*e1;
/*------------------------------------------------------------------------------
      register new point*/
      dum=new point2D_t;
      *dum=r;
      extraction.push_back(dum);
      
      w=vector2D_t(previous,r);
      check=hypot(w);
      previous=r;
      
      length+=radius;
      w=vector2D_t(r,q);
      }
    if(debug) printf("radius=%lf km, length=%lf km, line=%d\n",radius/1000.,length/1000.0,  j);
    }

/*------------------------------------------------------------------------------
  add last point */
  w=vector2D_t(r,q);
  if (hypot(w)!=0) {
    dum=new point2D_t;
    *dum=q;
    extraction.push_back(dum);
    length+=hypot(w);
    }

  n=extraction.size();
  
  misfit=hypot(w)/radius;
  
  double tolerance=0.25;
  
  tolerance=0.25*n/1000.0;
  updatemin(&tolerance,0.25);
  
/*------------------------------------------------------------------------------
  equalize sampling*/
  if((equalize==1) && (misfit>tolerance) && (misfit<1.0-tolerance)) {
/*------------------------------------------------------------------------------
    choose for larger or smaller radius according to last iteration*/
    if(previous_misfit<misfit) direction=direction-1;
    if(direction=0) {
/*------------------------------------------------------------------------------
      take a slightly larger radius*/
      radius=c*radius+(1.-c)*length/(n-2);
      }
    else {
/*------------------------------------------------------------------------------
      take a slightly smaller radius*/
      radius=c*radius+(1.-c)*length/(n-1);
      }
    for(size_t k=0;k<extraction.size();k++) delete extraction[k];
    extraction.clear();
    previous_misfit=misfit;
    goto begin;
    }

  n=extraction.size();
  if( (equalize==1) &&  (misfit<0.25) ) n=n-1;
  
  if(n==0) return(out);
  out.init(n,PLG_INIT_XY);

  for(i=0;i<n;i++) {
    out.x[i]=extraction[i]->x;
    out.y[i]=extraction[i]->y;
    }
  
  if(equalize==1) {
    out.x[n-1]=q.x;
    out.y[n-1]=q.y;
    }
  
  for(size_t k=0;k<extraction.size();k++) delete extraction[k];
  extraction.clear();
    
  return(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_duplicate(const vector<plg_t> & p, vector<plg_t> & q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  
  for(int k=0; k<p.size();k++) {
    plg_t *tmp=new plg_t;
    tmp->duplicate(p[k]);
    q.push_back(*tmp);
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_t plg_resample(const plg_t & target, const vector<double> & sizes, int equalize, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Resample an existing polygon following a given radius array, possibly tuned
  if equalizing is enabled
  
  resampling along target lines, but pre-existing nodes not preserved
  
  if segment changes, first interception is unique

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  int i,j,n;
  vector2D_t u,v;
  vector2D_t w1,w2,w;
  vector2D_t e1,e2;
  point2D_t p,q;
  point2D_t r,o;
  double a,r1,r2,length,check,misfit,previous_misfit=0.5;
  vector<point2D_t *> extraction;
  point2D_t *dum;
  plg_t out;
  point2D_t first,previous;
  vector<double> actual_radius;
  double c=0.95, radius;
  int direction=1;
  
  vector<double> prescribed;
  for(int k=0; k<sizes.size();k++) prescribed.push_back(sizes[k]);
  
begin:
  first.x=target.x[0];
  first.y=target.y[0];

  length=0;

  r=first;
  previous=r;
  
/*------------------------------------------------------------------------------
  register first point*/
  dum=new point2D_t;
  *dum=r;
  extraction.push_back(dum);
  
  radius=prescribed[min(extraction.size()-1,prescribed.size()-1)];
  actual_radius.push_back(radius);
  
  for(j=0;j<target.npt-1;j++) {
    p.x = target.x[j];
    p.y = target.y[j];
    q.x = target.x[j+1];
    q.y = target.y[j+1];
    w1=vector2D_t(r,p);
    w2=vector2D_t(r,q);
/*------------------------------------------------------------------------------
    check if this (new) segment is intercepted */
    r1=hypot(w1);
    r2=hypot(w2);
    if((r1-radius)*(r2-radius)>0) continue;
/*------------------------------------------------------------------------------
    first point on this segment*/
    status=math_orthobase(p, q, &e1, &e2);
    u=math_vector_coordinates01(w2, e1, e2);
    o=q+(-u.x)*e1;
    a=sqrt(radius*radius-u.y*u.y);
    r=o+a*e1;
/*------------------------------------------------------------------------------
    register new point*/
    dum=new point2D_t;
    *dum=r;
    extraction.push_back(dum);

    w=vector2D_t(previous,r);
    check=hypot(w);
    previous=r;
    
    length+=radius;
    if(debug) printf("radius=%lf km, length=%lf km, line=%d, point=%d\n",radius/1000.,length/1000.0,  j, extraction.size());
    
/*------------------------------------------------------------------------------
    further points on same segment*/
    radius=prescribed[min(extraction.size()-1,prescribed.size()-1)];
    actual_radius.push_back(radius);
    w=vector2D_t(r,q);
    while(hypot(w)>radius) {
      r=r+radius*e1;
/*------------------------------------------------------------------------------
      register new point*/
      dum=new point2D_t;
      *dum=r;
      extraction.push_back(dum);
      
      w=vector2D_t(previous,r);
      check=hypot(w);
      previous=r;
      
      length+=radius;
      w=vector2D_t(r,q);
      if(debug) printf("radius=%lf km, length=%lf km, line=%d, point=%d\n",radius/1000.,length/1000.0,  j, extraction.size());
      radius=prescribed[min(extraction.size()-1,prescribed.size()-1)];
      actual_radius.push_back(radius);
      }
    }

/*------------------------------------------------------------------------------
  add last point */
  w=vector2D_t(r,q);
  if (hypot(w)!=0) {
    dum=new point2D_t;
    *dum=q;
    extraction.push_back(dum);
    length+=hypot(w);
    }

  n=extraction.size();

  for(int k=actual_radius.size();k<prescribed.size();k++) actual_radius.push_back(prescribed[k]);
  misfit=hypot(w)/radius;
  
  c=1.-0.5/n;

/*------------------------------------------------------------------------------
  equalize sampling*/
  if((equalize==1) && (misfit>0.15) && (misfit<0.85)) {
/*------------------------------------------------------------------------------
    choose for larger or smaller radius according to last iteration*/
    if(previous_misfit<misfit) direction=1-direction;
    if(direction=0) {
/*------------------------------------------------------------------------------
      take a slightly larger radius*/
      c=1.+0.5*misfit/n;
      for(int k=0;k<actual_radius.size();k++) {
        actual_radius[k]=actual_radius[k]*c;
        }
      }
    else {
/*------------------------------------------------------------------------------
      take a slightly smaller radius*/
      c=1.-0.5*misfit/n;
      for(int k=0;k<actual_radius.size();k++) {
        actual_radius[k]=actual_radius[k]*c;
        }
      }
    for(size_t k=0;k<extraction.size();k++) delete extraction[k];
    extraction.clear();
    previous_misfit=misfit;
    prescribed.clear();
    for(int k=0;k<actual_radius.size();k++) prescribed.push_back(actual_radius[k]);
    actual_radius.clear();
    goto begin;
    }

  n=extraction.size();
  if( (equalize==1) &&  (misfit<0.25) ) n=n-1;
  
  if(n==0) return(out);
  out.init(n,PLG_INIT_XY);

  for(i=0;i<n;i++) {
    out.x[i]=extraction[i]->x;
    out.y[i]=extraction[i]->y;
    }
  
  if(equalize==1) {
    out.x[n-1]=q.x;
    out.y[n-1]=q.y;
    }
  
  for(size_t k=0;k<extraction.size();k++) delete extraction[k];
  extraction.clear();
    
  return(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t> plg_resample(const vector<plg_t> & target, double radius, int equalize, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k;
  vector<plg_t> out;
  
  for(k=0;k<target.size();k++) {
    plg_t *q=new plg_t;
    *q=plg_resample(target[k], radius, equalize, debug);
    out.push_back(*q);
    }
  
  return(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_t plg_sample_regular_cartesian(plg_t target, double radius, int equalize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  polygon resampling, cartesian version
  
  preserve first and last point only

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  int i,n;
  vector2D_t u,v,w;
  vector2D_t e1,e2;
  point2D_t p,q,r;
  double length;
  vector<point2D_t> extraction;
  point2D_t *dum;
  plg_t out;
  point2D_t first,last;

  first.x=target.x[0];
  first.y=target.y[0];
  last.x =target.x[target.npt-1];
  last.y =target.y[target.npt-1];

  status=math_orthobase(first, last, &e1, &e2);
  u=vector2D_t(first,last);

  n=NINT(hypot(u)/radius)+1;

  if(equalize==1) {
    radius=hypot(u)/double (n-1);
    }
  else {
    n=n-1;
    }
  length=0;

  r=first;
  dum=new point2D_t;
  *dum=r;
  extraction.push_back(*dum);
  delete dum;

  for(i=0;i<n-1;i++) {
    r=r+radius*e1;
/*------------------------------------------------------------------------------
    register new point*/
    dum=new point2D_t;
    *dum=r;
    extraction.push_back(*dum);
    delete dum;
    length+=radius;
    }

  if(equalize==0) {
    w=vector2D_t(r,last);
    dum=new point2D_t;
    *dum=last;
    extraction.push_back(*dum);
    delete dum;
    length+=hypot(w);
    }

  if(extraction.size()==1) extraction.push_back(last);
  n=extraction.size();
  
  out.init(n,PLG_INIT_SEPARATE);

  for(i=0;i<extraction.size();i++) {
    out.x[i]=extraction[i].x;
    out.y[i]=extraction[i].y;
    }

  extraction.clear();
  
  out.x[out.npt-1]=last.x;
  out.y[out.npt-1]=last.y;
  
  return(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_t plg_sample_regular_loxodromic(plg_t target, double radius, int equalize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  polygon resampling, cartesian version
  
  preserve first and last point only

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  int i,n;
  vector2D_t u,v,w;
  vector2D_t e1,e2;
  point2D_t p,q,r;
  double length;
  vector<point2D_t> extraction;
  point2D_t *dum;
  plg_t out;
  point2D_t first,last;

  first.x=target.t[0];
  first.y=target.p[0];
  last.x =target.t[target.npt-1];
  last.y =target.p[target.npt-1];

  status=math_orthobase(first, last, &e1, &e2);
  u=vector2D_t(first,last);

  n=NINT(hypot(u)/radius)+1;

  if(equalize==1) {
    radius=hypot(u)/double (n-1);
    }
  else {
    n=n-1;
    }
  length=0;

  r=first;
  dum=new point2D_t;
  *dum=r;
  extraction.push_back(*dum);
  delete dum;

  for(i=0;i<n-1;i++) {
    r=r+radius*e1;
/*------------------------------------------------------------------------------
    register new point*/
    dum=new point2D_t;
    *dum=r;
    extraction.push_back(*dum);
    length+=radius;
    }

  if(equalize==0) {
    w=vector2D_t(r,last);
    dum=new point2D_t;
    *dum=last;
    extraction.push_back(*dum);
    delete dum;
    length+=hypot(w);
    }

  if(extraction.size()==1) extraction.push_back(last);
  n=extraction.size();
  
  out.init(n,PLG_INIT_SEPARATE);

  for(i=0;i<extraction.size();i++) {
    out.t[i]=extraction[i].x;
    out.p[i]=extraction[i].y;
    }

  extraction.clear();
  
  out.t[out.npt-1]=last.x;
  out.p[out.npt-1]=last.y;
  
  return(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_t plg_subdivide(plg_t target, double radius, int equalize, int style, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  subdivide the polygon's segments at radius resolution
  
  assume (x,y) to be cartesian coordinates (m)
  
  projection/orthodromic/great circle option?

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int coordinates, status;
  plg_t out;
  
  for(int j=0;j<target.npt-1;j++) {
    plg_t p,q;
    switch (style) {
      case PLG_LOXODROMIC:
        p=plg_t(target.t[j],target.p[j],target.t[j+1],target.p[j+1]);
        q=plg_sample_regular_loxodromic(p, radius, equalize);
        coordinates=PLG_SPHERICAL;
        break;
      case PLG_CARTESIAN:
        p=plg_t(target.x[j],target.y[j],target.x[j+1],target.y[j+1]);
        q=plg_sample_regular_cartesian(p, radius, equalize);
        coordinates=PLG_CARTESIAN;
        break;
      }
//           if(target.flag!=0) p.SetFlag(target.flag[j]);

    if(target.flag!=0) q.SetFlag(target.flag[j]);
    status=plg_concat(out, q, coordinates);
    if(status==-1) {
      out.destroy();
      TRAP_ERR_EXIT(-1, " concat failed for sub-segment %d\n",j);
      }
    q.destroy();
    p.destroy();
    }
  return(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t> plg_subdivide(const vector<plg_t> & target, double radius, int equalize, int style, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<plg_t> out;
  
#pragma omp parallel for
  for(int k=0;k<target.size();k++) {
    plg_t q;
    q=plg_subdivide(target[k], radius, equalize, style, debug);
    #pragma omp critical(push)
    {
    out.push_back(q);
    }
    }
  
  return(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  projPJ plg_DefineProjection(plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  frame_t frame;
  projPJ projection;
  double ref_lat[2], ref_lon[2];

  frame= plg_spherical_minmax(polygones, npolygones);

  ref_lat[0]=(frame.ymax+frame.ymin)/2;
  ref_lon[0]=(frame.xmax+frame.xmin)/2;
  
/*------------------------------------------------------------------------------
  define Oblique Stereographic projection */
  projection=assign_StereoOblique(ref_lat[0],ref_lon[0]);

  status=plg_cartesian(projection,polygones, npolygones);

  frame= plg_cartesian_minmax(polygones, npolygones);

  return(projection);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  projPJ plg_DefineProjection(vector<plg_t> & polygons, char *parameters)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  frame_t frame;
  projPJ projection;
  double ref_lat[2], ref_lon[2];

  frame=plg_spherical_minmax(polygons);

  ref_lat[0]=(frame.ymax+frame.ymin)/2;
  ref_lon[0]=(frame.xmax+frame.xmin)/2;
  
/*------------------------------------------------------------------------------
  define Oblique Stereographic projection */
  projection=assign_StereoOblique(ref_lat[0],ref_lon[0], parameters);

  status=plg_cartesian(projection,polygons);

  frame= plg_cartesian_minmax(polygons);

  return(projection);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_t plg_sample_shorelines(double t1, double p1, double t2, double p2, plg_t *shorelines, int nshorelines, double radius, int equalize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  plg_t out;
  double x1,x2,y1,y2;
  double L1, L2, L;
  plg_target_t plg_target[2];
  int m1,m2,dum;
  plg_t target;
  vector2D_t u1,u2;
  
  int closed,island,reversed;
  
  projPJ projection=plg_DefineProjection(shorelines, nshorelines);
  
  status=plg_cartesian(projection, shorelines, nshorelines);
  
  geo_to_projection(projection, p1, t1, &x1, &y1);
  geo_to_projection(projection, p2, t2, &x2, &y2);
    
  plg_target[0].line=plg_NearestPolygon(shorelines, nshorelines, x1, y1, m1, dum, u1);
  plg_target[1].line=plg_NearestPolygon(shorelines, nshorelines, x2, y2, m2, dum, u2);
  
  if(plg_target[0].line!=plg_target[1].line) {
/*------------------------------------------------------------------------------
    to be further fixed */
    STDOUT_BASE_LINE("trouble\n");
    exit(-1);
    }
  
  size_t s=plg_target[0].line;
  
  closed=plg_isclosed(shorelines[s]);
  
  island=( (t1==t2) && (p1==p2) );
  
  if( (island==1) && (closed==0) ) {
    STDOUT_BASE_LINE("trouble\n");
    exit(-1);
    }
  
  if(island==0) {
    if(closed==1) {
      L1=plg_length(shorelines[s], m1, m2, CARTESIAN);
      L2=plg_length(shorelines[s], m2, m1, CARTESIAN);
      if(L2<L1) {
/*------------------------------------------------------------------------------
        choose the shortes shoreline section */
        int dum=m1;
        m1=m2;
        m2=dum;
        reversed=1;
        }
      else {
        reversed=0;
        }
      }
    target.duplicate(shorelines[s], m1, m2);
    if(reversed==0) {
      target.x[0]=x1;
      target.y[0]=y1;
      target.x[target.npt-1]=x2;
      target.y[target.npt-1]=y2;
      }
    else {
      target.x[0]=x2;
      target.y[0]=y2;
      target.x[target.npt-1]=x1;
      target.y[target.npt-1]=y1;
      }
    }
  else {
    target.duplicate(shorelines[s]);
    }
  
  L=plg_length(target,CARTESIAN);
  
  updatemin(&radius,L/10.0);
  
  out=plg_resample(target, radius, equalize);
  
  status=plg_spherical(projection, &out, 1);

  target.destroy();
  
  return(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t> plg_ConvertPolygonSet(PolygonSet polygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  vector<plg_t> set;
  int k;

/*------------------------------------------------------------------------------
  convert LeBars polygon set into Lyard polygon set*/
  for (PolygonSet::iterator i = polygons.begin(); i != polygons.end();i++) {
    Polygon p = *i;
    plg_t q;
    const size_t n = p.size()+1;
    q.npt=n;
    q.x=new double[q.npt];
    q.y=new double[q.npt];
    q.t=new double[q.npt];
    q.p=new double[q.npt];
    k=0;
    for (Polygon::iterator j = p.begin(); j != p.end();j++) {
      const Point a=*j;
      q.t[k]=a.x();
      q.p[k]=a.y();
      q.x[k]=a.x();
      q.y[k]=a.y();
      k++;
      }
    q.t[k]=q.t[0];
    q.p[k]=q.p[0];
    set.push_back(q);
    }
  return(set);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_array_t plg_readneigh (const string & polygonsFileName) throw ()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
using namespace Polygons;

  plg_array_t set;
  int k,m;

  try {
    /* Read the polygons in file polygonsFileName */
    int maxsize=-1;
    PolygonSet polygons = searchPolygons<double,double>(GraphBase<double> (polygonsFileName),maxsize);
/*------------------------------------------------------------------------------
    convert LeBars polygon set into Lyard polygon set*/
    set.n=polygons.size();
    set.p=new plg_t[set.n];

    m=0;
    for (PolygonSet::iterator i = polygons.begin(); i != polygons.end();i++) {
      Polygon p = *i;
      const size_t n = p.size()+1;
      set.p[m].npt=n;
      set.p[m].x=new double[set.p[m].npt];
      set.p[m].y=new double[set.p[m].npt];
      set.p[m].t=new double[set.p[m].npt];
      set.p[m].p=new double[set.p[m].npt];
      k=0;
      for (Polygon::iterator j = p.begin(); j != p.end();j++) {
        const Point a=*j;
        set.p[m].t[k]=a.x();
        set.p[m].p[k]=a.y();
        k++;
        }
      set.p[m].t[k]=set.p[m].t[0];
      set.p[m].p[k]=set.p[m].p[0];
      m++;
      }
    }
  
/*------------------------------------------------------------------------------
  Part of the code executed if a problem occur during file reading */
  catch (TugoExceptions::ReadError file) {
    TRAP_ERR_EXIT(-2, "An error occured when reading the file \"%s\".",file.fileName().c_str());
    }
  
/*------------------------------------------------------------------------------
  Part of the code executed if a reference point is on the boundary of a polygon */
  catch (TugoExceptions::OnBoundaryPoint point) {
    TRAP_ERR_EXIT(-3,
            "The reference point n°%l, which coordinates are (%d, %d) is on the boundary of a polygon.",
            point.ind(), point.x(), point.y());
    }

/*------------------------------------------------------------------------------
  Part of the code executed if trying to use a negative or null value for z0 */
  catch (TugoExceptions::Negative) {
    TRAP_ERR_EXIT(-4, "Trying to use a negative value for z0.");
    }
    
  return(set);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t> plg_readneigh(const string & polygonsFileName, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_array_t r=plg_readneigh (polygonsFileName);
  
  vector<plg_t> s=plg_array2vector(r);
  
  return(s);
  
}
