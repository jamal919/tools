
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

\brief plg_t processing (join, cut, etc...) functions
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>

#include "tools-structures.h"

#include "polygones.h"
#include "geo.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> void swapval(T *v1, T *v2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T tmp;
  tmp=*v1;
  *v1=*v2;
  *v2=tmp;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_flip(plg_t polygon)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k;

  for(k=0;k<polygon.npt/2;k++) {
    swapval(&(polygon.t[k]),&(polygon.t[polygon.npt-k-1]));
    swapval(&(polygon.p[k]),&(polygon.p[polygon.npt-k-1]));
    swapval(&(polygon.x[k]),&(polygon.x[polygon.npt-k-1]));
    swapval(&(polygon.y[k]),&(polygon.y[polygon.npt-k-1]));
    if(polygon.z!=0) swapval(&(polygon.z[k]),&(polygon.z[polygon.npt-k-1]));
    }
  if (polygon.flag!=0) {
    for(k=0;k<(polygon.npt-1)/2;k++) {
      swapval(&(polygon.flag[k]),&(polygon.flag[polygon.npt-1-k-1]));
      }
    }

  return(-0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_flip(plg_t *polygones, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=plg_flip(polygones[target]);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_merge (vector<plg_t> & p, vector<plg_t> q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  add polygons of q to p list; be aware that polygons are not duplicated

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  for(size_t k=0;k<q.size();k++) p.push_back(q[k]);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_concat (plg_t & pp, plg_t & qq, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  concat two existing polygones, trying all possible extremities

  mode given to plg_distance() :
    0:plg_distance_spherical()
    1:plg_distance_cartesian()

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k,kk,status;
  double distance;
  char *flag=0;
  double epsilon=1.e-10;
//   double *z=0;
  
  if(pp.npt==0) {
    pp.duplicate(qq);
    return(0);
    }
  if(qq.npt==0) return(-1);

  distance= plg_distance(pp, 0, qq, 0, mode);

  if(distance<epsilon) {
    status=plg_flip(pp);
    goto finish;
    }

  distance= plg_distance(pp, pp.npt-1, qq, 0, mode);

  if(distance<epsilon) {
    goto finish;
    }

  distance= plg_distance(pp, 0, qq, qq.npt-1, mode);

  if(distance<epsilon) {
    status=plg_flip(pp);
    status=plg_flip(qq);
    goto finish;
    }

  distance= plg_distance(pp, pp.npt-1, qq, qq.npt-1, mode);

  if(distance<epsilon) {
    status=plg_flip(qq);
    goto finish;
    }

  return(-1);

finish:
  plg_t buf;
  
  buf.init(pp.npt+qq.npt-1, PLG_INIT_SEPARATE);
  
  if(pp.flag!=0 and qq.flag!=0) {
    flag=new char[buf.npt-1];
    }
  
  if(pp.z!=0 and qq.z!=0) {
    buf.z=new double[buf.npt];
    }
  
  k=0;
  for(kk=0;kk<pp.npt;kk++) {
    plg_copy_point(&buf,k,pp,k);
    k++;
    }
  for(kk=1;kk<qq.npt;kk++) {
    plg_copy_point(&buf,k,qq,kk);
    k++;
    }
  
  k=0;
  for(kk=0;kk<pp.npt-1;kk++) {
    if(pp.flag!=0 and qq.flag!=0) flag[k]=pp.flag[k];
    k++;
    }
  for(kk=0;kk<qq.npt-1;kk++) {
    if(pp.flag!=0 and qq.flag!=0) flag[k]=qq.flag[kk];
    k++;
    }

  pp.destroy();
  qq.destroy();

  pp=buf;
  pp.flag=flag;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_concat (int target[2], plg_t *polygones, int npolygones, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=plg_concat(polygones[target[0]], polygones[target[1]], mode);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_concat (int target[2], vector<plg_t> & polygones, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=plg_concat(polygones[target[0]], polygones[target[1]], mode);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_concat (int p1, int p2, vector<plg_t> & polygones, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=plg_concat(polygones[p1], polygones[p2], mode);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_concat (plg_target_t target[2], plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// concat two existing polygones, with given extremities
/*----------------------------------------------------------------------------*/
{
  int k,kk,m,s,status;
  double t1,t2,p1,p2;
  double distance;

  s=target[0].line;
  m=target[0].point;
  t1=polygones[s].t[m];
  p1=polygones[s].p[m];

  s=target[1].line;
  m=target[1].point;
  t2=polygones[s].t[m];
  p2=polygones[s].p[m];

  distance=geo_distance(t1,p1,t2,p2);

  if(distance==0.0) {
    if(target[0].point==0) status=plg_flip(polygones,target[0].line);
    if(target[1].point!=0) status=plg_flip(polygones,target[1].line);
    goto finish;
    }

  return(-1);

finish:
  plg_t buf;
  buf.init(polygones[target[0].line].npt+polygones[target[1].line].npt-1,PLG_INIT_SEPARATE);
  k=0;
  for(kk=0;kk<polygones[target[0].line].npt;kk++) {
    plg_copy_point(&buf,k,polygones[target[0].line],k);
    k++;
    }
  for(kk=1;kk<polygones[target[1].line].npt;kk++) {
    plg_copy_point(&buf,k,polygones[target[1].line],kk);
    k++;
    }

  polygones[target[0].line].destroy();
  polygones[target[1].line].destroy();

  polygones[target[0].line]=buf;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_join (plg_t & p, plg_t & q, double join_radius)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  int k,nk,kmin,m[4][2];
  double t1,t2,p1,p2;
  double dmin,d[4];

  k=0;

  m[k][0]=p.npt-1;
  m[k][1]=0;
  k++;

  m[k][0]=0;
  m[k][1]=q.npt-1;
  k++;

  m[k][0]=0;
  m[k][1]=0;
  k++;

  m[k][0]=p.npt-1;
  m[k][1]=q.npt-1;
  k++;

  if(&p==&q) nk=2;
  else nk=4;

  dmin=1.e+10;
  for(k=0;k<nk;k++) {
    t1=p.t[m[k][0]];
    p1=p.p[m[k][0]];
    t2=q.t[m[k][1]];
    p2=q.p[m[k][1]];
    d[k]=geo_distance(t1,p1,t2,p2);
    if(d[k]<dmin) {
      kmin=k;
      dmin=d[k];
      }
    }

  if(dmin<join_radius) {
    k=kmin;
    plg_copy_point(&p,m[k][0], q,m[k][1]);
    d[k]=1.e+10;
    status=0;
    }

  dmin=1.e+10;
  for(k=0;k<nk;k++) {
    if(d[k]<dmin) {
      kmin=k;
      dmin=d[k];
      }
    }

  if(dmin<join_radius) {
    k=kmin;
    plg_copy_point(&p,m[k][0], q,m[k][1]);
    d[k]=1.e+10;
    status=0;
    }
  
  return(status);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_join (int target[2], plg_t *polygones, int npolygones, double join_radius)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=plg_join(polygones[target[0]], polygones[target[1]], join_radius);
  
  return(status);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_join (int target[2], vector<plg_t> polygones, double join_radius)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=plg_join(polygones[target[0]], polygones[target[1]], join_radius);
  
  return(status);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_swap_entries (plg_t *polygones, int m, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_t tmp;

  tmp=polygones[m];
  polygones[m]=polygones[n];
  polygones[n]=tmp;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_swap_entries (vector<plg_t> polygones, int m, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_t tmp;

  tmp=polygones[m];
  polygones[m]=polygones[n];
  polygones[n]=tmp;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *plg_cut_polygone (plg_t **polygones, int *npolygones, int target, int point)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n1,n2;
  int *list;
  int re_use=1;

  if(target<0 or target>=*npolygones)
    TRAP_ERR_EXIT(-1,"target=% when *npolygones=%d\n",target,*npolygones);
  
  if(re_use==1) {
    list=plg_increase_reuse_entries (polygones,npolygones,2);
    }
  else {
    list=plg_increase_entries (polygones,npolygones,2);
    }

  if(list==(int *) 0) {
    return(list);
    }
  n1=list[0];
  n2=list[1];

  (*polygones)[n1].init(point+1,PLG_INIT_SEPARATE);
  (*polygones)[n2].init((*polygones)[target].npt-point,PLG_INIT_SEPARATE);

  for(k=0;k<(*polygones)[n1].npt;k++) {
    plg_copy_point(&(*polygones)[n1],k,
                    (*polygones)[target],k);
    }

  for(k=0;k<(*polygones)[n2].npt;k++) {
    plg_copy_point(&(*polygones)[n2],k,
                    (*polygones)[target],k+point);
    }

  (*polygones)[target].destroy();

  return(list);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_cut_polygone(vector<plg_t> & polygons, int target, int point, int force)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  cut a polygon, increasing the size only if necessary

  point the cloned point

  return the index of the new polygon or -1 if point is out of polygones[target] range

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k;
  plg_t p1,p2;
  int list=-1;
  
  if( not force){
    if(point==0) {
      return(list);
      }
    if(point==polygons[target].npt-1) {
      return(list);
      }
    }

  p1.init(point+1,PLG_INIT_SEPARATE);
  p2.init(polygons[target].npt-point,PLG_INIT_SEPARATE);
  
  for(k=0;k<p1.npt;k++) {
    plg_copy_point(&p1,k,polygons[target],k);
    }

  for(k=0;k<p2.npt;k++) {
    plg_copy_point(&p2,k,polygons[target],k+point);
    }

  if(polygons[target].flag!=0) {
    p1.SetFlag(-1);
    p2.SetFlag(-1);
    for(k=0;k<p1.npt-1;k++) p1.flag[k]=polygons[target].flag[k];
    for(k=0;k<p2.npt-1;k++) p2.flag[k]=polygons[target].flag[k+point];
    }

//   plg_delete_polygone(polygones[target]);
//   polygones.erase (polygones.begin()+target-1);
//
//   list[0]=polygones.size();
//   polygones.push_back(p1);
//
//   list[1]=polygones.size();
//   polygones.push_back(p2);
  
  polygons[target].destroy();
  polygons[target]=p1;
  
  list=plg_add(&polygons,p2);
  
  return(list);
}
