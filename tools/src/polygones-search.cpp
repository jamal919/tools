
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

\brief plg_t search functions
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>

#include "polygones.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int CheckPolar(vector<plg_t> polygons, projPJ &proj)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  Patch for circum-polar polygons
------------------------------------------------------------------------------*/
{
  int   k,n,status;
  double p;
  bool north;
  double lon, lat, x, y;
//   projPJ proj;
  
  for(k=0;k<polygons.size();k++) {
    if(!polygons[k].closed()) {
      TRAP_ERR_EXIT(-1, "polygons not closed");
      return(-1);
      }
    }
  p=polygons[0].p[0];
  
  if(p>0) north=true;
  else north=false;
  
  for(k=0;k<polygons.size();k++) {
    for (n=0;n<polygons[k].npt;n++) {
      if(north && polygons[k].p[n] <0) {
        TRAP_ERR_EXIT(-1, "polygons are ambiguous");
        return(-1);
        }
      if(!north && polygons[k].p[n] >0) {
        TRAP_ERR_EXIT(-1, "polygons are ambiguous");
        return(-1);
        }
      }
    }
  
  if(north) {
    proj=assign_StereoOblique(+90.0, 0.0);
    lon= 0.0;
    lat=90.0;
    }
  else  {
    proj=assign_StereoOblique(-90.0, 0.0);
    lon=  0.0;
    lat=-90.0;
    }
 
  geo_to_projection(proj, lat, lon, &x, &y);
  status=plg_cartesian(proj, polygons);
     
  int in=0;
  int inside=plg_single( polygons[0], x, y, &in, PLG_CARTESIAN);
  
  status=plg_checkAutoSecant(polygons, PLG_CARTESIAN);
  if(status==-1 and inside==PLG_POINT_EXTERIOR) {
    printf("polygones are NOT polar and auto-secant, this is an issue\n");
    return(-1);
    }
  else if(status==0 and inside!=PLG_POINT_EXTERIOR){
    printf("Ok, polygones are polar\n");
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int CheckPolar(plg_t & p, projPJ &proj)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  Patch for circum-polar polygons
------------------------------------------------------------------------------*/
{
  int status;

  vector<plg_t> polygons;
  
  polygons.push_back(p);
  
  status=CheckPolar(polygons, proj);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> frame_t plg_spherical_minmax_core_template(const T & polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  frame_t frame;
  range_t<double> rx,ry;
  
  for(i=0;i<npolygones;i++) {
    rx.init(polygones[i].t, polygones[i].npt);
    ry.init(polygones[i].p, polygones[i].npt);
    if(i==0) {
      frame.xmin=rx.min;
      frame.xmax=rx.max;
      frame.ymin=ry.min;
      frame.ymax=ry.max;
      }
    else {
      updatemin(&frame.xmin,rx.min);
      updatemax(&frame.xmax,rx.max);
      updatemin(&frame.ymin,ry.min);
      updatemax(&frame.ymax,ry.max);
      }
    }

  return(frame);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  frame_t plg_spherical_minmax(const plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  frame_t frame;
  
  if(polygones==0) return(frame);

  frame=plg_spherical_minmax_core_template(polygones,npolygones);

  return(frame);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  frame_t plg_spherical_minmax(const vector<plg_t> & polygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  frame_t frame;
  
  frame=plg_spherical_minmax_core_template(polygons,polygons.size());

  return(frame);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> frame_t plg_cartesian_minmax_core_template(const T & polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  frame_t frame;
  range_t<double> rx,ry;
  
  for(i=0;i<npolygones;i++) {
    const plg_t *plg=&polygones[i];
    
    rx.init(plg->x, plg->npt);
    ry.init(plg->y, plg->npt);
    if(i==0) {
      frame.xmin=rx.min;
      frame.xmax=rx.max;
      frame.ymin=ry.min;
      frame.ymax=ry.max;
      }
    else {
      updatemin(&frame.xmin,rx.min);
      updatemax(&frame.xmax,rx.max);
      updatemin(&frame.ymin,ry.min);
      updatemax(&frame.ymax,ry.max);
      }
    }

  return(frame);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  frame_t plg_cartesian_minmax(const plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  frame_t frame;
  
  if(polygones==0) return(frame);

  frame=plg_cartesian_minmax_core_template(polygones,npolygones);

  return(frame);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  frame_t plg_cartesian_minmax(const vector<plg_t> & polygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  frame_t frame;
  
  frame=plg_cartesian_minmax_core_template(polygons,polygons.size());

  return(frame);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

frame_t plg_cartesian_minmax(const plg_t & polygone)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  frame_t frame;
  
  frame=plg_cartesian_minmax(&polygone,1);

  return(frame);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_identify(int id, plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;

  for(i=0;i<npolygones;i++) {
    if(polygones[i].id==id) return(i);
    }
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_isclosed(const plg_t & polygon)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int closed, m1, m2;
  double distance;
  
  if(polygon.npt<1)
    return 1;
  
  m1=0;
  m2=polygon.npt-1;
  distance=geo_distance(polygon.t[m1],polygon.p[m1],polygon.t[m2],polygon.p[m2]);
  
  if(distance==0.0)
    closed=1;
  else {
    closed=0;
    }
  
  return(closed);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int close_plg(plg_t *polygon)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(plg_isclosed(*polygon)==1)
    return 0;
  
  int status;
  
  plg_point_t point(*polygon,0);
  
  status= plg_insert_point(*polygon,polygon->npt, point);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double plg_howfar(plg_t *polygones, int s, int m, double t, double p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double distance;
  
  distance=geo_distance(polygones[s].t[m],polygones[s].p[m],t,p);
  
  return(distance);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_find_point_cartesian(plg_t & polygon, double x, double y, double *d)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j,m;
  double distance,dmin,dx,dy;

  m=-1;
  dmin=1.e+10;

  for(j=0;j<polygon.npt;j++) {
    dx=polygon.x[j]-x;
    dy=polygon.y[j]-y;
    distance=dx*dx+dy*dy;
    if(distance<dmin) {
      m=j;
      dmin=distance;
      }
    }

  *d=sqrt(dmin);
  return(m);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_find_point_spherical(plg_t & polygon, double t, double p, double *d)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  warning: returns kilometeric distances; to be re-worked
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int j,m;
  double *distance,dmin;

  m=-1;
  dmin=1.e+10;

  distance=new double[polygon.npt];
  
  int nprocs=initialize_OPENMP(-1, 0);
#pragma omp parallel for if(nprocs>1)
  for(j=0;j<polygon.npt;j++) {
    distance[j]=geo_haversin_km(polygon.t[j],polygon.p[j],t,p);
    }

  for(j=0;j<polygon.npt;j++) {
    if(distance[j]<dmin) {
      m=j;
      dmin=distance[j];
      }
    }

  delete[] distance;
  
  *d=dmin;
  return(m);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_find_point(plg_t & polygon, int mode, double t, double p, double *d)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m;

  switch(mode) {
    case PLG_SPHERICAL:
      m=plg_find_point_spherical(polygon, t, p, d);
      break;
    case PLG_CARTESIAN:
      m=plg_find_point_cartesian(polygon, t, p, d);
      break;
    default:
      TRAP_ERR_EXIT(ENOEXEC,"programming error\n");
    }
  return(m);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_find_point(plg_t & polygon, int mode, plg_point_t & point, double *d)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m;

  m=plg_find_point(polygon, mode, point.t, point.p, d);
  return(m);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  paire_t plg_find_polygon(vector<plg_t> & polygons, plg_point_t point, double & dmin)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double d;
  paire_t p;
  
  dmin=1.e+10;
  
  p.value[0]=-1;
  p.value[1]=-1;

  for(int k=0; k<polygons.size();k++) {
//     if(plg_isclosed(polygons[k])==1) continue; // HERE!!!
    if(polygons[k].npt<=0) continue;
    int n=plg_find_point(polygons[k], PLG_SPHERICAL, point, &d);
    if(d<dmin) {
      p.value[0]=k;
      p.value[1]=n;
      dmin=d;
      }
    }
  
  return(p);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  paire_t plg_find_polygon(vector<plg_t> & polygons, vector<int> & list, plg_point_t point, double & dmin)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double d;
  paire_t p;
  
  dmin=1.e+10;
  
  p.value[0]=-1;
  p.value[1]=-1;

  for(int kk=0; kk<list.size();kk++) {
    int k=list[kk];
//     if(plg_isclosed(polygons[k])==1) continue; // HERE!!!
    int n=plg_find_point(polygons[k], PLG_SPHERICAL, point, &d);
    if(d<dmin) {
      p.value[0]=k;
      p.value[1]=n;
      dmin=d;
      }
    }
  
  return(p);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_NearestSegment(plg_t & polygon, double x, double y, int & m1, int & m2, double & dmin, vector2D_t & u)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j;
  double x1,x2,y1,y2,dx,dy,dx12,dy12;
  double d,ro;
  double utx,uty,unx,uny,xt,xn;
  bool test;
  
  dmin=1.e+10;

  for(j=0;j<polygon.npt;j++) {
    x1 = polygon.x[j];
    y1 = polygon.y[j];
    dx = x-x1;
    dy = y-y1;
    d=sqrt(dx*dx+dy*dy);
    if (d < dmin) {
      m1=j;
      m2=j;
      dmin=d;
      u.x=dx;
      u.y=dy;
      }
    }
  
  for(j=0;j<polygon.npt-1;j++) {
    x1 = polygon.x[j];
    y1 = polygon.y[j];
    x2 = polygon.x[j+1];
    y2 = polygon.y[j+1];
    if(fabs(polygon.t[j+1]-polygon.t[j]) > 180.0) goto next;
    dx12= x2-x1;
    dy12= y2-y1;
    dx  = x-x1;
    dy  = y-y1;
    ro=sqrt(dx12*dx12+dy12*dy12);
    if(ro == 0.0) {
      goto next;
      }
    utx=dx12/ro;
    uty=dy12/ro;
    xt=dx*utx+dy*uty;
    test=(xt*(xt-ro) <= 0.0);
    if (test) {
      unx=-uty;
      uny=+utx;
      xn=dx*unx+dy*uny;
      d=fabs(xn);
      if (d < dmin) {
        m1=j;
        m2=j+1;
        dmin=d;
        u.x=xn*unx;
        u.y=xn*uny;
        }
      }
  next:
    x1=x2;
    y1=y2;
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_NearestPolygon(plg_t *polygons, int npolygons, double x, double y, int & m1, int & m2, vector2D_t & u)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,status;
  int s=-1;
  double d,dmin;
  int mm1,mm2;
  vector2D_t uu;
  dmin=1.e+10;

  for(i=0;i<npolygons;i++) {
    status=plg_NearestSegment(polygons[i], x, y, mm1, mm2, d, uu);
    if(d<dmin) {
      s=i;
      m1=mm1;
      m2=mm2;
      u=uu;
      dmin=d;
      }
    }
  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_NearestPolygon(vector<plg_t> & polygons, double x, double y, int & m1, int & m2, vector2D_t & u)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,status;
  int s=-1;
  double d,dmin;
  int mm1,mm2;
  vector2D_t uu;
  dmin=1.e+10;

  for(i=0;i<polygons.size();i++) {
    status=plg_NearestSegment(polygons[i], x, y, mm1, mm2, d, uu);
    if(d<dmin) {
      s=i;
      m1=mm1;
      m2=mm2;
      u=uu;
      dmin=d;
      }
    }
  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_NearestPolygon(vector<plg_t> & polygons, vector<int> & exclusion, double x, double y, int & m1, int & m2, vector2D_t & u)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,pos,status;
  int s=-1;
  double d,dmin;
  int mm1,mm2;
  vector2D_t uu;
  dmin=1.e+10;

  for(i=0;i<polygons.size();i++) {
    pos=vpos(i,exclusion);
    if(pos!=-1) continue;
    status=plg_NearestSegment(polygons[i], x, y, mm1, mm2, d, uu);
    if(d<dmin) {
      s=i;
      m1=mm1;
      m2=mm2;
      u=uu;
      dmin=d;
      }
    }
  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double plg_distance_spherical (const plg_t & plg_a, int vertex_a, const plg_t & plg_b, int vertex_b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double t1,t2,p1,p2;
  double distance;
  int m;

  m=vertex_a;
  t1=plg_a.t[m];
  p1=plg_a.p[m];

  m=vertex_b;
  t2=plg_b.t[m];
  p2=plg_b.p[m];

//  distance=geo_distance(t1,p1,t2,p2);
  distance=geo_haversin_km(t1,p1,t2,p2);

  return(distance);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double plg_distance_cartesian (const plg_t & plg_a, int vertex_a, const plg_t & plg_b, int vertex_b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double x1,x2,y1,y2,dx,dy;
  double distance;
  int m;

  m=vertex_a;
  x1=plg_a.x[m];
  y1=plg_a.y[m];

  m=vertex_b;
  x2=plg_b.x[m];
  y2=plg_b.y[m];
  
  dx=x2-x1;
  dy=y2-y1;

  distance=sqrt(dx*dx+dy*dy);

  return(distance);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double plg_distance(const plg_t & plg_a, int vertex_a, const plg_t & plg_b, int vertex_b, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double distance;

  switch (mode) {
    case 0:
      distance=plg_distance_spherical(plg_a, vertex_a, plg_b, vertex_b);
      break;
    case 1:
      distance=plg_distance_cartesian(plg_a, vertex_a, plg_b, vertex_b);
      break;
    default:
      TRAP_ERR_EXIT(ENOEXEC,"programming error\n");
    }
  
  return(distance);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double plg_distance(plg_t plg_a, int vertex_a, int vertex_b, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double distance;

  distance=plg_distance_spherical(plg_a, vertex_a, plg_a, vertex_b);
  
  return(distance);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double plg_distance(const plg_t & plg_a, const plg_t & plg_b, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double d,distance=1.e+10;
  int k, l;

  switch (mode) {
    case 0:
      for(l=0;l<plg_b.npt;l++) {
        for(k=0;k<plg_a.npt;k++) {
          d=plg_distance_spherical(plg_a, k, plg_b, l);
          updatemin(&distance,d);
          }
        }
      break;
    case 1:
      for(l=0;l<plg_b.npt;l++) {
        for(k=0;k<plg_a.npt;k++) {
          d=plg_distance_cartesian(plg_a, k, plg_b, l);
          updatemin(&distance,d);
          }
        }
      break;
    default:
      TRAP_ERR_EXIT(ENOEXEC,"programming error\n");
    }
  
  return(distance);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double plg_length(plg_t plg_a, int vertex_a, int vertex_b, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double length=0;
  int reversed=0;
  
  if(vertex_a<vertex_b) {
    reversed=0;
    }
  else{
    reversed=1;
    int dum=vertex_b;
    vertex_b=vertex_a;
    vertex_a=dum;
    }
  

  switch (mode) {
    case 0:
      for(size_t k=vertex_a;k<vertex_b;k++) {
        length+=plg_distance_spherical(plg_a, k, plg_a, k+1);
        }
      break;
    case 1:
      for(size_t k=vertex_a;k<vertex_b;k++) {
        length+=plg_distance_cartesian(plg_a, k, plg_a, k+1);
        }
      break;
    }
  
  if(reversed==1) {
    length=plg_length(plg_a, mode)-length;
    }
  
  return(length);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double plg_length(plg_t plg_a, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double length;

  if(plg_a.npt<2) return(0.0);
  
  length=plg_length(plg_a, 0, plg_a.npt-1, mode);
  
  return(length);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double plg_minlength(plg_t plg_a, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double length, minlength=+INFINITY;

  if(plg_a.npt<2) return(0.0);
  
  for(int k=0;k<plg_a.npt;k++) {
    length=plg_length(plg_a, k, k+1, mode);
    updatemin(&minlength, length);
    }
  
  return(minlength);
}


