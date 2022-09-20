
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

\brief plg_t edition (add, remove, copy, ...) functions
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>

#include "tools-structures.h"

#include "polygones.h"
#include "geo.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void plg_copy_point(double *x, double *y, double *t, double *p, const plg_t & src, int srcPoint)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  *t=src.t[srcPoint];
  *p=src.p[srcPoint];
  *x=src.x[srcPoint];
  *y=src.y[srcPoint];
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_move_point(plg_t & dest, int destPoint, double x, double y, double t, double p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(dest.x!=0 && dest.y!=0){
    dest.x[destPoint]=x;
    dest.y[destPoint]=y;
    }
  
  if(isnan(t) && isnan(p) && dest.x==dest.t && dest.y==dest.p)
    return PLG_INIT_SHARED;
  
  dest.t[destPoint]=t;
  dest.p[destPoint]=p;
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void plg_copy_point(plg_t *dest, int destPoint, const plg_t & src, int srcPoint)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
  patch to avoid allocation mismatch, could be confusing */

  if(dest->t!=0 and src.t!=0) {
    dest->t[destPoint]=src.t[srcPoint];
    dest->p[destPoint]=src.p[srcPoint];
    }
  if(dest->x!=0 and src.x!=0) {
    dest->x[destPoint]=src.x[srcPoint];
    dest->y[destPoint]=src.y[srcPoint];
    }
  if(dest->z!=0 and src.z!=0) {
    dest->z[destPoint]=src.z[srcPoint];
    }
//   if(dest->flag!=0 and src.flag!=0 ) {
//     dest->flag[destPoint]=src.flag[srcPoint];
//     }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void plg_copy_point(plg_t *polygone, int dest, int src)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_copy_point(polygone,dest,*polygone,src);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_delete_point(plg_t *polygon, int point)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k;

  for(k=point;k<polygon->npt-1;k++) {
    plg_copy_point(polygon,k,k+1);
    }
  if(polygon->flag!=0) {
    for(k=point;k<polygon->npt-2;k++) {
      polygon->flag[k]=polygon->flag[k+1];
      }
    }
  polygon->npt--;

  return(-0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_point_t plg_newpoint_cartesian(plg_t *polygones, plg_interval_t interval)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double x1,y1,x2,y2;
  plg_point_t  point;

  x1=polygones[interval.line].x[interval.point[0]];
  y1=polygones[interval.line].y[interval.point[0]];

  x2=polygones[interval.line].x[interval.point[1]];
  y2=polygones[interval.line].y[interval.point[1]];

  point.x=0.5*(x1+x2);
  point.y=0.5*(y1+y2);

//  status=geo_retroction(&point.t,&point.p,&point.x,&point.y);

  return(point);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_insert_point(plg_t & polygon, int position, const plg_point_t & point)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  insert point in polygon at postion
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k,m;
  plg_t buf;
  
  buf.init(polygon.npt+1,PLG_INIT_SEPARATE);

  for(k=0;k<position;k++) {
    plg_copy_point(&buf,k,polygon,k);
    }
  for(k=position+1;k<buf.npt;k++) {
    plg_copy_point(&buf,k,polygon,k-1);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  flag, position :
  
    * at extremities, reproduce first/last flag
    * in interior, reproduce the segment flag
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(polygon.flag!=0) {
    buf.flag=new char[buf.npt-1];

/*------------------------------------------------------------------------------
    */
    if (position==0) {
      buf.flag[0]=polygon.flag[0];
      for(k=1;k<buf.npt-1;k++) {
        buf.flag[k]=polygon.flag[k-1];
        }
      }
    else if (position==polygon.npt) {
      for(k=0;k<buf.npt-2;k++) {
        buf.flag[k]=polygon.flag[k];
        }
      buf.flag[buf.npt-2]=polygon.flag[polygon.npt-2];
      }
    else {
      for(k=0;k<position;k++) {
        buf.flag[k]=polygon.flag[k];
        }
      buf.flag[position]=polygon.flag[position-1];
      for(k=position+1;k<buf.npt-1;k++) {
        buf.flag[k]=polygon.flag[k-1];
        }
      }
    }
  
  polygon.destroy();

  polygon=buf;

  plg_move_point(polygon,position,point);
  

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_insert_point(plg_t *polygones, int target, int position, const plg_point_t & point)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=plg_insert_point(polygones[target], position, point);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_colocalize(plg_t *polygones, plg_target_t target[2])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int s1, m1, s2, m2;

  s1=target[0].line;
  s2=target[1].line;
  m1=target[0].point;
  m2=target[1].point;
  
  plg_copy_point(&polygones[s2],m2,polygones[s1],m1);

  return(-0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_colocalize(plg_t *polygones, int s1, int m1, int s2, int m2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  plg_target_t target[2];

  target[0].line=s1;
  target[1].line=s2;
  target[0].point=m1;
  target[1].point=m2;
  status=plg_colocalize(polygones, target);
  return(-0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_colocalize(plg_t & p, int m, plg_t & q, int n, int option)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  double t1, p1, t2, p2, tmid, pmid;
  double x1, y1, x2, y2;
  
  t1=p.t[m];
  t2=q.t[m];
  
  p1=p.p[m];
  p2=q.p[m];
  
  tmid=0.5*(t1+geo_recale(t2, t1, 180.0));
  
  p.t[m]=tmid;
  
  tmid=0.5*(geo_recale(t1, t2, 180.0)+t2);
  
  q.t[m]=tmid;
  
  pmid=0.5*(p1+p2);
  
  p.p[m]=pmid;
  q.p[m]=pmid;
  
  return(-0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_move_point(plg_t & polygon, int position, const plg_point_t & point)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=plg_move_point(polygon,position, point.x, point.y, point.t, point.p);

  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_move_point(plg_t *polygones, int target, int position, const plg_point_t & point)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=plg_move_point(polygones[target], position , point);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void plg_distord (plg_t *polygones,plg_target_t target, double xx, double yy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m,k,s;
  double r,d,x,y,t,p,dx,dy;

  s=target.line;
  m=target.point;

  dx=xx-polygones[s].x[m];
  dy=yy-polygones[s].y[m];

  for(k=0;k<polygones[s].npt;k++) {
    plg_copy_point(&x,&y,&t,&p,polygones[s],k);
    d=fabs((float) m - (float) k);
    if( plg_isclosed(polygones[s])==1) {
      d=min(d,polygones[s].npt-1-d);
      }
    if(d>2.0) {
      r=0.;
      continue;
      }
    else {
      r=fabs(cos(d/2.0));
      }
    polygones[s].x[k]+=dx*r;
    polygones[s].y[k]+=dy*r;
//    ok=geo_retroction(&t,&p,&(polygones[s].x[k]),&(polygones[s].y[k]));
    polygones[s].t[k]=t;
    polygones[s].p[k]=p;
    }

// IF(scan_magneticbehaviour==ScanN_MAGNETIC_RULERS) THEN
//   !---------------------------------------------------------------------
//   ! Force x and y to coincide with the rulers grid
//   !---------------------------------------------------------------------
//   ok=map_pointingrid(axes_rulersgrid,x,y,k,l,rstatus)
//   x=map_gridx(axes_rulersgrid,k,l,rstatus)
//   y=map_gridy(axes_rulersgrid,k,l,rstatus)
// END IF

// dx=x-scan_x(mouse_s,mouse_p)
// dy=y-scan_y(mouse_s,mouse_p)

// IF (MOUSE_MOVEd) THEN
//   DO i=1,scan_npt(mouse_s)
//     IF(scan_closed(mouse_s,rstatus)) THEN
//       d=min(ABS(i-mouse_p),scan_npt(mouse_s)-ABS(i-mouse_p))
//     ELSE
//       d=ABS(i-mouse_p)
//     ENDIF
//     IF(d.LE.2) THEN
//       r=ABS(COS(d/2.0));
//     ELSE
//       r=0.
//       CYCLE
//     END IF
//     scan_x(mouse_s,i)=scan_x(mouse_s,i)+r*dx
//     scan_y(mouse_s,i)=scan_y(mouse_s,i)+r*dy
//     xx=scan_x(mouse_s,i)
//     yy=scan_y(mouse_s,i)
//     ok=geo_retroction (tt,pp,xx,yy)
//     scan_t(mouse_s,i)=tt
//     scan_p(mouse_s,i)=pp
//   END DO
//   scan_modified=1
//   MOUSE_MOVEd = .false.
// END IF
//
// CALL main_setmousecursor (MOUSE_CURSOR_CROSS)
//
// CALL mousef_pop_call (MOUSE_MOVE)
//
// CALL scan_unselect (0)
//
//   main_rebuilddefaultforeground ();
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_compact_entries (plg_t **polygones, int *npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  remove empty polygons from initial set */
{
  int k,kk,n,empty;
  plg_t *tmp;

  empty=0;
  for(k=0;k<*npolygones;k++) {
    if((*polygones)[k].npt==0) {
      empty++;
      }
    }

  n=*npolygones-empty;

  tmp=new plg_t[n];

  kk=0;
  for(k=0;k<*npolygones;k++) {
    if((*polygones)[k].npt!=0) {
      tmp[kk]=(*polygones)[k];
      kk++;
      }
    }

  delete[] (*polygones);

  *polygones=tmp;
  *npolygones=n;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_compact_entries (vector<plg_t> & polygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  remove empty polygons from initial set */
{
  int k,empty;

  empty=0;
  for(k=0;k<polygons.size();k++) {
    if(polygons[k].npt==0) {
      empty++;
      }
    }

  if(empty==0) return(0);
  
  for(k=polygons.size()-1;k>=0;k--) {
    if(polygons[k].npt==0) {
      polygons[k].destroy();
      polygons.erase(polygons.begin()+k);
      }
    }
    
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int *plg_increase_entries (plg_t **polygones, int *npolygones, int additional)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n,s;
  plg_t *tmp;
  int *list;

  list=new int[additional];

  n=*npolygones+additional;

  s=0;
  for(k=*npolygones;k<n;k++) {
    list[s]=k;
    s++;
    }


  tmp=new plg_t[n];
  for(k=0;k<*npolygones;k++) {
    tmp[k]=(*polygones)[k];
    }

  delete[] (*polygones);

  *polygones=tmp;
  *npolygones=n;
  return(list);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int *plg_increase_reuse_entries (plg_t **polygones, int *npolygones, int additional)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n,s;
  plg_t *tmp;
  int *list;

  if(additional<1) {
    return((int *) 0);
    }

  list=new int[additional];
  if(list==NULL) {
    return((int *) 0);
    }

  s=0;
  for(k=0;k<*npolygones;k++) {
    if((*polygones)[k].npt==0) {
      list[s]=k;
      s++;
      additional--;
      }
    if(additional==0) return(list);
    }

  n=*npolygones+additional;

  for(k=*npolygones;k<n;k++) {
    list[s]=k;
    s++;
    }

  tmp=new plg_t[n];
  for(k=0;k<*npolygones;k++) {
    tmp[k]=(*polygones)[k];
    }

  delete[] (*polygones);

  *polygones=tmp;
  *npolygones=n;
  return(list);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_add(vector<plg_t> *array,const plg_t & plg)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// add a polygone, increasing the size only if necessary
/**
\return the index of the added polygone
*/
/*----------------------------------------------------------------------------*/
{
  int i;
  const int n=array->size();
  
  /* searching an empty polygone */
  for(i=0;i<n && (*array)[i].npt;i++);
  
  if(i<n)
    (*array)[i]=plg;
  else
    array->push_back(plg);

  return i;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_add(vector<plg_t> *array,double x, double y, double t, double p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// add an isolated point, increasing the size only if necessary
/**
\return the index of the added polygone
*/
/*----------------------------------------------------------------------------*/
{
  plg_t plg;
  int i;
  plg_init_mode mode;
  
  if(isnan(t) && isnan(p))
    mode=PLG_INIT_SHARED;
  else
    mode=PLG_INIT_SEPARATE;
  
  plg.init(1,mode);
  
  plg_move_point(plg,0,x,y,t,p);

  i=plg_add(array,plg);

  return i;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

extern "C" int plg_add_(plg_array_t *array,int *npt, double *x, double *y, double *t, double *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n;
  int *list;
  plg_t **polygones;

  polygones=&(array->p);

  list=plg_increase_reuse_entries (polygones,&(array->n),1);

  if(list==(int *) 0) {
    return(-1);
    }
  n=list[0];

  (*polygones)[n].init(*npt,PLG_INIT_SEPARATE);

  for(k=0;k<*npt;k++) {
    plg_move_point((*polygones)[n],k,x[k],y[k],t[k],p[k]);
    }

  delete [] list;

  return(n);
}
