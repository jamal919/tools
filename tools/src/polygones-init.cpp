
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

\brief plg_t initialisation, including projection
*/
/*----------------------------------------------------------------------------*/

#include <stdio.h>

#include "constants.h"
#include "polygones.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_orientation (const plg_t & plg_polygones, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double alpha;
  double c,angle,pdv,pds,dx1,dx2,dy1,dy2,dn1,dn2;
  int    size,k,km,kp;
  int    orientation;

  size=plg_polygones.npt-1;

  angle=0;
  for(k=0;k<size;k++) {
    km=(k+size-1) % size;
    kp=(k+size+1) % size;
    dx1=plg_polygones.x[k]-plg_polygones.x[km];
    dy1=plg_polygones.y[k]-plg_polygones.y[km];
    dn1=sqrt(dx1*dx1+dy1*dy1);
    dx2=plg_polygones.x[kp]-plg_polygones.x[k];
    dy2=plg_polygones.y[kp]-plg_polygones.y[k];
    dn2=sqrt(dx2*dx2+dy2*dy2);
    dx1/=dn1;
    dy1/=dn1;
    dx2/=dn2;
    dy2/=dn2;
    pdv=dx1*dy2-dx2*dy1;
    pds=dx1*dx2+dy1*dy2;
    if(pdv>0.0) {
      alpha= acos(pds);
      }
    else {
      alpha=-acos(pds);
      }
/*------------------------------------------------------------------------------
    it return alpha in [-180,180] interval*/
    alpha=alpha*r2d;
    angle+=alpha;
//    printf("km=%d,k=%d,kp=%d, alpha=%lf angle=%lf\n",km,k,kp,alpha,angle);
    }
  if(angle > 0.0) {
    orientation=+1;
    }
  else {
    orientation=-1;
    }

  if(debug) printf("angle=%lf, orientation=%d\n",angle,orientation);
  return(orientation);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_setdirect(plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,status,orientation;
  bool debug=false;
  
  for(i=0;i<npolygones;i++) {
    orientation=plg_orientation (polygones[i], debug);
    if(orientation<0) {
      printf("set polygon %d to direct rotation order\n",i);
      status=plg_flip(polygones,i);
      orientation=plg_orientation (polygones[i], debug);
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_setdirect(vector<plg_t> & polygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,status,orientation;
  bool debug=false;
  
  for(i=0;i<polygons.size();i++) {
    orientation=plg_orientation (polygons[i], debug);
    if(orientation<0) {
      printf("set polygon %d to direct rotation order\n",i);
      status=plg_flip(polygons[i]);
      orientation=plg_orientation (polygons[i], debug);
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void plg_t::init(int n,const plg_init_mode shared)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  bool doXY=false,doTP=false;
  
  destroy();
  NULLPointers();
  
  switch(shared){
    case PLG_INIT_SEPARATE:
      doXY=true;
      doTP=true;
      break;
    case PLG_INIT_SHARED:
    case PLG_INIT_XY:
      doXY=true;
      break;
    case PLG_INIT_TP:
      doTP=true;
      break;
    }
  if(n!=0) {
    if(doXY){
      x=new double[n];
      y=new double[n];
      }
    if(doTP){
      t=new double[n];
      p=new double[n];
      }
    else if(shared==PLG_INIT_SHARED){
      t=x;
      p=y;
      }
    }
  
  npt=n;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void plg_deletep(plg_t **polygones,const int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  for(int i=0;i<npolygones;i++){
    (*polygones)[i].destroy();
    }
  
  deletep(polygones);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_destroy_entries(vector<plg_t> & polygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k;

  for(k=polygons.size()-1;k>=0;k--) {
    polygons[k].destroy();
    }

  polygons.clear();
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_add_entry(vector<plg_t> & polygons, plg_t p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  plg_t *tmp=new plg_t;
  tmp->duplicate(p);
  polygons.push_back(*tmp);

  return 0;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_add_entries(vector<plg_t> & polygons, plg_t *p, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  for(int i=0;i<n;i++){
    plg_t *tmp=new plg_t;
    tmp->duplicate(p[i]);
    polygons.push_back(*tmp);
    }
  
  return 0;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_add_entries(vector<plg_t> & polygons, vector<plg_t> & p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  for(int i=0;i<p.size();i++){
    plg_t *tmp=new plg_t;
    tmp->duplicate(p[i]);
    polygons.push_back(*tmp);
    }
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t> plg_duplicate(const vector<plg_t> & p0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  vector<plg_t> p1;
  
  for(i=0;i<p0.size();i++){
    plg_t tmp;
    tmp.duplicate(p0[i]);
    p1.push_back(tmp);
    }
  
  return p1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t>  plg_array2vector(plg_array_t polygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<plg_t> p;
  
  copy(&p,polygones.p,polygones.n);
  
  return(p);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t>  plg_array2vector(plg_t *polygons, int npolygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<plg_t> p;
  
  copy(&p,polygons, npolygons);
  
  return(p);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_array_t plg_vector2array(vector<plg_t> polygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_array_t array(polygones.size());
  
  for(int k=0;k<polygones.size();k++) array.p[k]=polygones[k];
  
  return(array);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_vector2array(vector<plg_t> polygons, plg_t* & p, int & np)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  p=new plg_t[polygons.size()];
  np=polygons.size();
  
  for(int k=0;k<polygons.size();k++) p[k].duplicate(polygons[k]);
    
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  frame_t plg_recale(vector<plg_t> & polygons, double center)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j;
  frame_t frame;
  double t,t0;
  
//   printf("warning, plg_recale changed!!!\n");

  if(isnan(center)) t0=polygons[0].t[0];
  else t0=center;
  
  for(i=0;i<polygons.size();i++) {
//     for(j=0;j<polygones[i].npt;j++) {
//       t=polygones[i].t[j];
//       polygones[i].t[j]=geo_recale(t,t0,(double) 180.0);
//       }
    j=0;
    t=polygons[i].t[j];
    polygons[i].t[j]=geo_recale(t,t0,(double) 180.0);
    for(j=1;j<polygons[i].npt;j++) {
      t=polygons[i].t[j];
      polygons[i].t[j]=geo_recale(t,polygons[i].t[j-1],(double) 180.0);
      }
    }
    
  frame=plg_spherical_minmax(polygons);
  return(frame);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  frame_t plg_recale(plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j;
  frame_t frame;
  double t,t0;
  
//  printf("warning, plg_recale changed!!!\n");

  t0=polygones[0].t[0];
  for(i=0;i<npolygones;i++) {
//     for(j=0;j<polygones[i].npt;j++) {
//       t=polygones[i].t[j];
//       polygones[i].t[j]=geo_recale(t,t0,(double) 180.0);
//       }
    j=0;
    t=polygones[i].t[j];
    polygones[i].t[j]=geo_recale(t,t0,(double) 180.0);
    for(j=1;j<polygones[i].npt;j++) {
      t=polygones[i].t[j];
      polygones[i].t[j]=geo_recale(t,polygones[i].t[j-1],(double) 180.0);
      }
    }
    
  frame=plg_spherical_minmax(polygones, npolygones);
  return(frame);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_recale(vector<plg_t> & polygons, frame_t frame, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,pivot,intersection;
  double t,t0,x,y;
  double center=frame.x_center();
  bool ok;
  
  for(i=0;i<polygons.size();i++) {
/*------------------------------------------------------------------------------
    first find a point inside frame */    
    intersection=0;
    for(j=0;j<polygons[i].npt;j++) {
      x=polygons[i].x[j];
      y=polygons[i].y[j];
      if(mode==PLG_SPHERICAL) x=geo_recale(x, center, 180.0);
      ok=frame.inside(x, y);
      if(ok) {
        intersection++;
        pivot=j;
        break;
        }
      }
    if(intersection!=0) {
/*------------------------------------------------------------------------------
      then use as reference for geo_recale */    
      polygons[i].x[pivot]=geo_recale(x, center, 180.0);
      polygons[i].t[pivot]=polygons[i].x[pivot];
      for(j=pivot+1;j<polygons[i].npt;j++) {
        polygons[i].x[j]=geo_recale(polygons[i].x[j], polygons[i].x[j-1], 180.0);
        polygons[i].t[j]=polygons[i].x[j];
        }
      for(j=pivot-1;j>=0;j--) {
        polygons[i].x[j]=geo_recale(polygons[i].x[j], polygons[i].x[j+1], 180.0);
        polygons[i].t[j]=polygons[i].x[j];
        }
      }
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_recale(plg_t *polygones, int npolygones, frame_t frame, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  vector<plg_t> tmp;
  
  copy(&tmp,polygones, npolygones);
  
  status=plg_recale(tmp, frame, mode);
  
  tmp.clear();
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_spherical(projPJ ref, plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j;

  for(i=0;i<npolygones;i++) {
    if(polygones[i].t==0) {
      polygones[i].t=new double[polygones[i].npt];
      }
    if(polygones[i].p==0) {
      polygones[i].p=new double[polygones[i].npt];
      }
    for(j=0;j<polygones[i].npt;j++) {
      projection_to_geo(ref, &polygones[i].p[j],&polygones[i].t[j],(polygones[i].x[j]),(polygones[i].y[j]));
      }
    }
    
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void plg_degree_recale(vector<plg_t> *polygons, double center)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// recale the coastline
/*----------------------------------------------------------------------------*/
{
  for(int i=0;i<polygons->size();i++){
    const plg_t *plg=&(*polygons)[i];
    double *t=plg->t;
    for(int k=0;k<plg->npt;k++)
      degree_recale(&t[k],center);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_spherical(projPJ ref, vector<plg_t> & polygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j;

  for(i=0;i<polygons.size();i++) {
    if(polygons[i].t==0) {
      polygons[i].t=new double[polygons[i].npt];
      }
    if(polygons[i].p==0) {
      polygons[i].p=new double[polygons[i].npt];
      }
    for(j=0;j<polygons[i].npt;j++) {
      projection_to_geo(ref, &polygons[i].p[j],&polygons[i].t[j],(polygons[i].x[j]),(polygons[i].y[j]));
      }
    }
    
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> void plg_cartesian_core_template(projPJ ref, T & polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  compute cartesian coordinates

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i,j;

  for(i=0;i<npolygones;i++) {
    plg_t *plg=&polygones[i];
    
    if(plg->x==0 or plg->x==plg->t) {
      plg->x=new double[plg->npt];
//      printf("warning, x array not allocated; apply plan B\n");
      }
    if(plg->y==0 or plg->y==plg->p) {
      plg->y=new double[plg->npt];
//      printf("warning, y array not allocated; apply plan B\n");
      }
    for(j=0;j<plg->npt;j++) {
      double t=degree_recale(plg->t[j],0.0);
      geo_to_projection(ref, plg->p[j],t,&(plg->x[j]),&(plg->y[j]));
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_cartesian(projPJ ref, plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_cartesian_core_template(ref,polygones,npolygones);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_cartesian(projPJ ref, vector<plg_t> & polygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_cartesian_core_template(ref,polygones,polygones.size());
  
  return(0);
}
