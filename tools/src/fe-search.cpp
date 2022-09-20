
/*******************************************************************************

  T-UGO tools, 2006-2013

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief element search functions, often for interpolation initialisation
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits>

#include "tools-structures.h"
#include "constants.h"

#include "functions.h"

#include "fe.def"
#include "fe.h"

#include "map.h"
#include "rutin.h"
#include "geo.h"
#include "maths.h"

int fe_lastelement=0;
int fe_nlistmax=1000;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_nearest_vertex(const mesh_t & mesh,double x,double y,int *except,int nexcept)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,n,p;
  double dmin,dx,dy,d,xx;

  dmin = 1.e+10;

  for (i = 0;i<mesh.nvtxs;i++) {
    for( n=0;n<nexcept;n++)
      if(i == except[n]) goto skip;
/*    dx = mesh.vertices[i].lon-x;*/
    xx=mesh.vertices[i].lon;
    xx=degree_recale(xx,x);
    dx = xx-x;
    dy = mesh.vertices[i].lat-y;
    d  = dx*dx+dy*dy;
    if (d < dmin) {
      dmin   = d;
      p = i;
      }
    skip:
    continue;
    }

  return(p);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_nearest_boundaryvertex(const mesh_t & mesh,double x,double y,int *except,int nexcept)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,n,p;
  double dmin,dx,dy,d,xx;

  dmin = 1.e+10;

  for (i = 0;i<mesh.nvtxs;i++) {
    if(mesh.vertices[i].code<=0) continue;
    for( n=0;n<nexcept;n++)
      if(i == except[n]) goto skip;
/*    dx = mesh.vertices[i].lon-x;*/
    xx=mesh.vertices[i].lon;
    xx=degree_recale(xx,x);
    dx = xx-x;
    dy = mesh.vertices[i].lat-y;
    d  = dx*dx+dy*dy;
    if (d < dmin) {
      dmin   = d;
      p = i;
      }
    skip:
    continue;
    }

  return(p);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_nearest_vertex(const mesh_t & mesh, double x, double y, int code, double & d)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,m;
  double dmin,dx,dy,xx;
  double t,p;

  dmin = 1.e+10;

  for (i = 0;i<mesh.nvtxs;i++) {
    if(mesh.vertices[i].code!=code)
      continue;
    
    xx=mesh.vertices[i].lon;
    xx=degree_recale(xx,x);
    dx = xx-x;
    dy = mesh.vertices[i].lat-y;
    d  = dx*dx+dy*dy;
    if (d < dmin) {
      dmin   = d;
      m = i;
      }
    }

  t=mesh.vertices[m].lon;
  p=mesh.vertices[m].lat;
  
  d=geo_haversin(x,y,t,p);
  
  return(m);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_nearest_node(const discretisation_t & descriptor,double x,double y,int *except,int nexcept)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int p;
  double dmin;

  dmin = 1.e+10;
  p=-1;

#pragma omp for
  for (int i = 0;i<descriptor.nnodes;i++) {
    double dx,dy,d,xx;
    for(int n=0;n<nexcept;n++)
      if(i == except[n]) goto skip;
    xx=descriptor.nodes[i].lon;
    xx=degree_recale(xx,x);
    dx = xx-x;
    dy =descriptor.nodes[i].lat-y;
    d  = dx*dx+dy*dy;
    if (d < dmin) {
      #pragma omp critical(update)
        {
        dmin = d;
        p    = i;
        }
      }
    skip:
    continue;
    }
  return(p);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_whichnearest(const mesh_t & mesh,double t,double p,double *tnew,double *pnew)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,n1,n2,n3,found;
  double t1,t2,t3,p1,p2,p3,pinf,psup,tinf,tsup,d;
  double d1,d2,d3,dmin,dmax,tbar,pbar;
  double tt,pp;

  dmin=1.e+10;
  found=0;

  for (i=0;i<mesh.ntriangles;i++) {
  /*  for (i=max(fe_lastelement-250,0);i<min(fe_lastelement+250,mesh.ntriangles);i++)*/
    n1=mesh.triangles[i].vertex[0];
    n2=mesh.triangles[i].vertex[1];
    n3=mesh.triangles[i].vertex[2];
    t1=mesh.vertices[n1].lon;
    t2=mesh.vertices[n2].lon;
    t3=mesh.vertices[n3].lon;
    t2=degree_recale(t2,t1);
    t3=degree_recale(t3,t1);
    p1=mesh.vertices[n1].lat;
    p2=mesh.vertices[n2].lat;
    p3=mesh.vertices[n3].lat;
    tinf=min(t1,t2);
    updatemin(&tinf,t3);
    tsup=max(t1,t2);
    updatemax(&tsup,t3);
    pinf=min(p1,p2);
    updatemin(&pinf,p3);
    psup=max(p1,p2);
    updatemax(&psup,p3);
    t=degree_recale(t,0.5*(tinf+tsup));
    tbar=(t1+t2+t3)/3.0;
    pbar=(p1+p2+p3)/3.0;
    d=fabs(t-tbar)*cos(pbar*d2r)+fabs(p-pbar);
    if (d < dmin) {
      found=i;
      dmin=d;
      }
    }
  if(found==0) return(-1);

  n1=mesh.triangles[found].vertex[0];
  n2=mesh.triangles[found].vertex[1];
  n3=mesh.triangles[found].vertex[2];
  t1=mesh.vertices[n1].lon;
  t2=mesh.vertices[n2].lon;
  t3=mesh.vertices[n3].lon;
  t2=degree_recale(t2,t1);
  t3=degree_recale(t3,t1);
  p1=mesh.vertices[n1].lat;
  p2=mesh.vertices[n2].lat;
  p3=mesh.vertices[n3].lat;
  tinf=min(t1,t2);
  updatemin(&tinf,t3);
  tsup=max(t1,t2);
  updatemax(&tsup,t3);
  pinf=min(p1,p2);
  updatemin(&pinf,p3);
  psup=max(p1,p2);
  updatemax(&psup,p3);
  t=degree_recale(t,0.5*(tinf+tsup));
  tbar=(t1+t2+t3)/3.0;
  pbar=(p1+p2+p3)/3.0;

  d1=geo_km(t,p,t1,p1);
  d2=geo_km(t,p,t2,p2);
  d3=geo_km(t,p,t3,p3);

  dmax=max(d1,d2);
  updatemax(&dmax,d3);
  
/*
  printf("***************\n");
  printf("found=%d\n",found);
  printf("t,p=%f,%f\n",t,p);
  printf("t1,t2,t3=%f,%f,%f\n",t1,t2,t3);
  printf("p1,p2,p3=%f,%f,%f\n",p1,p2,p3);
  printf("d1,d2,d3=%f,%f,%f\n",d1,d2,d3);
  printf("dmax=%f\n",dmax);
*/
 
  if(dmax > 100.) return(-1);

  if(d1 == dmax) {
    tt=d3/(d2+d3)*t2+d2/(d2+d3)*t3;
    pp=d3/(d2+d3)*p2+d2/(d2+d3)*p3;
    *tnew=0.1*t1+0.9*tt;
    *pnew=0.1*p1+0.9*pp;
    }
  else if (d2 == dmax) {
    tt=d3/(d1+d3)*t1+d1/(d1+d3)*t3;
    pp=d3/(d1+d3)*p1+d1/(d1+d3)*p3;
    *tnew=0.1*t2+0.9*tt;
    *pnew=0.1*p2+0.9*pp;
    }
  else {
    tt=d2/(d1+d2)*t1+d1/(d1+d2)*t2;
    pp=d2/(d1+d2)*p1+d1/(d1+d2)*p2;
    *tnew=0.1*t3+0.9*tt;
    *pnew=0.1*p3+0.9*pp;
    }

  d=geo_km(t,p,*tnew,*pnew);
/*    printf("%f %f %f %f  %f %f %f %f %f %f\n",d1,d2,d3,d,t,p,tt,pp,*tnew,*pnew);  */
/*    printf("%f %f %f %f %f %f\n",t1,t2,t3,p1,p2,p3);  */
  return(found);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_testelement_v1(const mesh_t & mesh,int i,float t,float p,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
     sous-programme donnant le signe du produit vectoriel de deux
     vecteurs de r*r; ces vecteurs ont pour origine ou extremite
     l'un des sommets d'un triangle et un point dont on veut
     determiner s'il appartient ou non au triangle.
     
     en sortie, le parametre status permettra de resoudre la
     question: status=0 si le point m appartient au triangle
     status=1 sinon
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int n1,n2,n3;
  float t1,t2,t3,p1,p2,p3;
  double dt1,dp1,dt2,dp2,pdt12,dt3,dp3;
  double pdt23,pdt31,dror,d12,d23,d31,epsilon;

  epsilon=-1.e-04;

  n1=mesh.triangles[i].vertex[0];
  n2=mesh.triangles[i].vertex[1];
  n3=mesh.triangles[i].vertex[2];

  t1=mesh.vertices[n1].lon;
  t2=mesh.vertices[n2].lon;
  t3=mesh.vertices[n3].lon;

  if(mode==0) {
    t1=degree_recale(t1,t);
    t2=degree_recale(t2,t1);
    t3=degree_recale(t3,t1);
    }
  p1=mesh.vertices[n1].lat;
  p2=mesh.vertices[n2].lat;
  p3=mesh.vertices[n3].lat;

/*--------------------- produit vectoriel m1,m2 ------------------------*/

  dt1=(double) t1 - (double) t;
  dp1=(double) p1 - (double) p;
  dt2=(double) t2 - (double) t;
  dp2=(double) p2 - (double) p;

  pdt12=dt1*dp2-dt2*dp1;

/*-------------------------- produit m2,m3 ------------------------------*/

  dt3=(double) t3 - (double) t;
  dp3=(double) p3 - (double) p;

  pdt23=dt2*dp3-dt3*dp2;

/*-------------------------- produit m3,m1 ------------------------------*/

  pdt31=dt3*dp1-dt1*dp3;

  dror=pdt12+pdt23+pdt31;

  if (dror == 0.) {
    return(-1);
    }

  d12=pdt12/dror;

  if (d12 < epsilon) {
   return(0);
   }

  d23=pdt23/dror;

  if (d23 < epsilon) {
   return(0);
   }

  d31=pdt31/dror;
  if (d31 < epsilon) {
   return(0);
   }

/*   printf("%lf %lf %lf %lf %lf %lf %f %f\n",t1,t2,t3,p1,p2,p3,t,p); */
/*   printf("%lf %lf %lf %lf %lf %lf\n",dt1,dt2,dt3,dp1,dp2,dp3); */
/*   printf("%lf %lf %lf %lf %d %d %d %d\n",dror,d12,d23,d31,n1,n2,n3,i); */
  return(1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_testelement_v2(const mesh_t & mesh,const triangle_t & q,double t,double p,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Checks whether a point is inside a triangle.

  mesh unstructured grid
  q triangle
  t longitude
  p latitude
  mode if 0: degree_recale() the longitues
  
  returns
   2 if the triangle is polar and the point is closer to the pole than the corner furthest from the pole,
  -1 if the triangle is flat,
   1 if the point is in the triangle,
   0 if not.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int n1,n2,n3;
  double t1,t2,t3,p1,p2,p3;
  double dt1,dp1,dt2,dp2,pdt12,dt3,dp3;
  double pdt23,pdt31,dror,d12,d23,d31,epsilon;

  epsilon=-1.e-04;

/*------------------------------------------------------------------------------
  gets the coordinates of the corners of the triangle */
  n1=q.vertex[0];
  n2=q.vertex[1];
  n3=q.vertex[2];

  t1=mesh.vertices[n1].lon;
  t2=mesh.vertices[n2].lon;
  t3=mesh.vertices[n3].lon;

  p1=mesh.vertices[n1].lat;
  p2=mesh.vertices[n2].lat;
  p3=mesh.vertices[n3].lat;
  
#define NORTH_POLE_EXP

  if(mode==0){
    /* degree_recale() the longitudes */
    t1=degree_recale(t1,t);
    t2=degree_recale(t2,t1);
    t3=degree_recale(t3,t1);
#ifdef NORTH_POLE_EXP
/*------------------------------------------------------------------------------
    north pole longitude is ill-defined (singularity) */
    if(p1==90) {
      t3 = degree_recale(t3, t2);
      t1=0.5*(t2+t3);
//       printf("north pole patch applies\n");
      }
    if(p2==90) {
      t3 = degree_recale(t3, t1);
      t2=0.5*(t1+t3);
//       printf("north pole patch applies\n");
      }
    if(p3==90) {
      t3 = degree_recale(t2, t1);
      t3=0.5*(t1+t2);
//       printf("north pole patch applies\n");
      }
#endif
    
    if(fabs(t3-t2)>180.){
/*------------------------------------------------------------------------------
      triangle is polar */
      if(p1>0){
        /* triangle is North of the Equator */
        if( p>=p1 or
            p>=p2 or
            p>=p3 )
          return 2;
        }
      else{
        /* triangle is South of the Equator */
        if( p<=p1 or
            p<=p2 or
            p<=p3 )
          return 2;
        }
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  It calculates the vectorial products of the 3 angles made by the point and 
  each pair of corners.
  If it is 0, the point is on the side joining these corners. 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*--------------------- produit vectoriel m1,m2 ------------------------*/

  dt1=(double) t1 - (double) t;
  dp1=(double) p1 - (double) p;
  dt2=(double) t2 - (double) t;
  dp2=(double) p2 - (double) p;
  
  pdt12=dt1*dp2-dt2*dp1;

/*-------------------------- produit m2,m3 ------------------------------*/

  dt3=(double) t3 - (double) t;
  dp3=(double) p3 - (double) p;

  pdt23=dt2*dp3-dt3*dp2;

/*-------------------------- produit m3,m1 ------------------------------*/

  pdt31=dt3*dp1-dt1*dp3;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  It calculates the sum of these vectorial products whose absolute value is the
  area of the triangle and whose sign is + or - whether the triangle is, 
  respectively, counter-clockwise or clockwise.
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  dror=pdt12+pdt23+pdt31;

/*------------------------------------------------------------------------------
  If the triangle has an area of 0 */
  if (dror == 0.) {
    return(-1);
    }
    
/*------------------------------------------------------------------------------
  multiplying by pow(dror,-1.) rather than dividing by dror makes things a lot faster */
  dror=pow(dror,-1);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  If any of the vectorial product has a sign that is different from their sum</h1>
  if means the point is outside of the triangle, so it returns 0
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  d12=pdt12*dror;

  if (d12 < epsilon) {
   return(0);
   }

  d23=pdt23*dror;

  if (d23 < epsilon) {
   return(0);
   }

  d31=pdt31*dror;
  if (d31 < epsilon) {
   return(0);
   }
 
/*------------------------------------------------------------------------------
  Otherwise it means the point is in the triangle, so it returns 1 */

  return(1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_testelement_v2(const mesh_t & mesh,const quadrangle_t & q,double t,double p,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Checks whether a point is inside a rectangle.
/**
\param mesh unstructured grid
\param q rectangle
\param t longitude
\param p latitude
\param mode if 0: degree_recale() the longitues
\returns 1 if the point is in the rectangle, 0 if not.
\sa fe_testelement_v2(mesh_t,triangle_t,double,double,int)
*/
/*----------------------------------------------------------------------------*/
{
  int n1,n2,n3,n4;
  double t1,t2,t3,t4,p1,p2,p3,p4;
  vector2D_t u,v;

  const double epsilon=-1.e-09;

/*------------------------------------------------------------------------------
  gets the coordinates of the corners of the triangle */
  n1=q.vertex[0];
  n2=q.vertex[1];
  n3=q.vertex[2];
  n4=q.vertex[3];

  t1=mesh.vertices[n1].lon;
  t2=mesh.vertices[n2].lon;
  t3=mesh.vertices[n3].lon;
  t4=mesh.vertices[n4].lon;

  p1=mesh.vertices[n1].lat;
  p2=mesh.vertices[n2].lat;
  p3=mesh.vertices[n3].lat;
  p4=mesh.vertices[n4].lat;

  if(mode==0){
    /* degree_recale() the longitudes */
    t1=degree_recale(t1,t);
    t2=degree_recale(t2,t1);
    t3=degree_recale(t3,t1);
    t4=degree_recale(t4,t1);
    }

/*--------------------- produit vectoriel m1,m2 ------------------------*/
  u=vector2D_t(t,p,t1,p1);
  v=vector2D_t(t,p,t2,p2);
  if (u*v < epsilon) {
   return(0);
   }
  u=vector2D_t(t,p,t2,p2);
  v=vector2D_t(t,p,t3,p3);
  if (u*v < epsilon) {
   return(0);
   }
  u=vector2D_t(t,p,t3,p3);
  v=vector2D_t(t,p,t4,p4);
  if (u*v < epsilon) {
   return(0);
   }
  u=vector2D_t(t,p,t4,p4);
  v=vector2D_t(t,p,t1,p1);
  if (u*v < epsilon) {
   return(0);
   }

  return(1);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_testelement_v3(const mesh_t & mesh,const triangle_t & q, double t, double p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Checks whether a point is inside a triangle.

  returns
   1 if the point is in the triangle,
   0 if not.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int n1,n2,n3;
  double t1,t2,t3,p1,p2,p3;
  const double epsilon=-1.e-09;
  double alpha,scalar;
  vector3_t o,u,v,w;
  
/*------------------------------------------------------------------------------
  gets the coordinates of the corners of the triangle */
  n1=q.vertex[0];
  n2=q.vertex[1];
  n3=q.vertex[2];

  t1=mesh.vertices[n1].lon;
  t2=mesh.vertices[n2].lon;
  t3=mesh.vertices[n3].lon;

  p1=mesh.vertices[n1].lat;
  p2=mesh.vertices[n2].lat;
  p3=mesh.vertices[n3].lat;
  
  double angle;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  It calculates angle between point and successive triangle vertices
  
  if all angles positive, point is inside
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  returns unit vectors */
  o=math_polar2cartesian(t*d2r, p*d2r);
  u=math_polar2cartesian(t1*d2r,p1*d2r);
  v=math_polar2cartesian(t2*d2r,p2*d2r);
  
/*------------------------------------------------------------------------------
  operator * is scalar product */
  u-=o*(o*u);
  v-=o*(o*v);
  
  angle=u^v;
  if(angle<0) return(0);

  angle=v^w;
  if(angle<0) return(0);

  w=math_polar2cartesian(t3*d2r,p3*d2r);
  w-=o*(o*w);
  
  angle=w^u;
  if(angle<0) return(0);

//   angle=geo_angle_radian(t, p, t1, p1, t2, p2);
//   if(angle<0) return(0);
// 
//   angle=geo_angle_radian(t, p, t2, p2, t3, p3);
//   if(angle<0) return(0);
// 
//   angle=geo_angle_radian(t, p, t3, p3, t1, p1);
//   if(angle<0) return(0);

  return(1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_testelement_v3(const mesh_t & mesh, const triangle_t & q, int m, double *t, double *p, int *elements, size_t *index, int ndata)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Checks whether a point is inside a triangle.

  returns
   1 if the point is in the triangle,
   0 if not.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int n1,n2,n3;
  double t1,t2,t3,p1,p2,p3;
  const double epsilon=-1.e-09;
  double alpha,scalar;
  vector3_t o,u,v,w;
  
/*------------------------------------------------------------------------------
  gets the coordinates of the corners of the triangle */
  n1=q.vertex[0];
  n2=q.vertex[1];
  n3=q.vertex[2];

  t1=mesh.vertices[n1].lon;
  t2=mesh.vertices[n2].lon;
  t3=mesh.vertices[n3].lon;

  p1=mesh.vertices[n1].lat;
  p2=mesh.vertices[n2].lat;
  p3=mesh.vertices[n3].lat;
  
  double angle;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  It calculates angle between point and successive triangle vertices
  
  if all angles positive, point is inside
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  returns unit vectors */
  u=math_polar2cartesian(t1*d2r,p1*d2r);
  v=math_polar2cartesian(t2*d2r,p2*d2r);
  w=math_polar2cartesian(t3*d2r,p3*d2r);
  
  for(size_t nn=0;nn<ndata;nn++) {
    size_t n;
    vector3_t uu,vv,ww;
    if(index!=0) n=index[nn];
    else n=nn;
    if(elements[n]!=-1) continue;
    o=math_polar2cartesian(t[n]*d2r, p[n]*d2r);
/*------------------------------------------------------------------------------
    operator * is scalar product */
    uu=u-o*(o*u);
    vv=v-o*(o*v);
  
    angle=uu^vv;
    if(angle<0) continue;

    ww=w-o*(o*w);
  
    angle=vv^ww;
    if(angle<0) continue;

    angle=ww^uu;
    if(angle<0) continue;
    
    elements[n]=m;
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_testelement_v3(const mesh_t & mesh,const quadrangle_t & q, double t, double p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Checks whether a point is inside a triangle.

  returns
   1 if the point is in the triangle,
   0 if not.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int n1,n2,n3,n4;
  double t1,t2,t3,t4,p1,p2,p3,p4;
  const double epsilon=-1.e-09;

/*------------------------------------------------------------------------------
  gets the coordinates of the corners of the triangle */
  n1=q.vertex[0];
  n2=q.vertex[1];
  n3=q.vertex[2];
  n4=q.vertex[2];

  t1=mesh.vertices[n1].lon;
  t2=mesh.vertices[n2].lon;
  t3=mesh.vertices[n3].lon;
  t4=mesh.vertices[n4].lon;

  p1=mesh.vertices[n1].lat;
  p2=mesh.vertices[n2].lat;
  p3=mesh.vertices[n3].lat;
  p4=mesh.vertices[n4].lat;
  
  double angle;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  It calculates angle between point and successive triangle vertices
  
  if all angles positive, point is inside
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  angle=geo_angle_radian(t, p, t1, p1, t2, p2);
  if(angle<0) return(0);

  angle=geo_angle_radian(t, p, t2, p2, t3, p3);
  if(angle<0) return(0);

  angle=geo_angle_radian(t, p, t3, p3, t4, p4);
  if(angle<0) return(0);

  angle=geo_angle_radian(t, p, t4, p4, t1, p1);
  if(angle<0) return(0);

  return(1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_whichelement_v2(const mesh_t & mesh,double t,double p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k[10],actual,next[10],m;
  int depth=0;
  int mode=0; /*Thierry le 03/06/2004*/

/*------------------------------------------------------------------------------
  multi-thread unsafe, to re-worked */
  actual=fe_lastelement;

  if (fe_testelement_v2(mesh,mesh.triangles[actual],t,p,mode) == 1) {
    fe_lastelement=actual;
    return(actual);
    }

  for (depth=0,k[depth]=0;k[depth]<mesh.elt_nngh[actual];k[depth]++) {
    next[depth]=mesh.elt_nghl[actual][k[depth]];
    if (fe_testelement_v2(mesh,mesh.triangles[next[depth]],t,p,mode) == 1) {
      fe_lastelement=next[depth];
      return(next[depth]);
      }
    }
  for (depth=0,k[depth]=0;k[depth]<mesh.elt_nngh[actual];k[depth]++) {
    next[depth]=mesh.elt_nghl[actual][k[depth]];
    for (depth=1,k[depth]=0;k[depth]<mesh.elt_nngh[next[depth-1]];k[depth]++) {
      next[depth]=mesh.elt_nghl[next[depth-1]][k[depth]];
      if (fe_testelement_v2(mesh,mesh.triangles[next[depth]],t,p,mode) == 1) {
        fe_lastelement=next[depth];
        return(next[depth]);
        }
      }
    for (depth=1,k[depth]=0;k[depth]<mesh.elt_nngh[next[depth-1]];k[depth]++) {
      next[depth]=mesh.elt_nghl[next[depth-1]][k[depth]];
      for (depth=2,k[depth]=0;k[depth]<mesh.elt_nngh[next[depth-1]];k[depth]++) {
        next[depth]=mesh.elt_nghl[next[depth-1]][k[depth]];
        if (fe_testelement_v2(mesh,mesh.triangles[next[depth]],t,p,mode) == 1) {
          fe_lastelement=next[depth];
          return(next[depth]);
          }
        }
      depth=1;
      }
    depth=0;
    }

  for (m=fe_lastelement;m<mesh.ntriangles;m++) {
    if (fe_testelement_v2(mesh,mesh.triangles[m],t,p,mode) == 1) {
      fe_lastelement=m;
      return(m);
      break;
      }
    }

  for (m=fe_lastelement;m>=0;m--) {
    if (fe_testelement_v2(mesh,mesh.triangles[m],t,p,mode) == 1) {
      fe_lastelement=m;
      return(m);
      break;
      }
    }

  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_whichelement(const mesh_t & mesh,double t,double p, int hint)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  return the element number containing the point lon,lat=(t,p)
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int m,n1,n2,n3;
  int inside;
  double t1,t2,t3,p1,p2,p3,pinf,psup,tinf,tsup;
  double c1,c2;
  char  eligible;
  int mode=0;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  if hint (element) elligible, try it first
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(hint>-1 and hint < mesh.ntriangles) {
    m=hint;
    inside=fe_testelement_v2(mesh,mesh.triangles[m],t,p,mode);
    if (inside == 1) {
      return(m);
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  then scan all elements list
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  multi-thread unsafe, to re-worked */
//   for (i=fe_lastelement;i<mesh.ntriangles+fe_lastelement;i++) {
//     m=i%mesh.ntriangles;
  for (m=0;m<mesh.ntriangles;m++) {
    n1=mesh.triangles[m].vertex[0];
    n2=mesh.triangles[m].vertex[1];
    n3=mesh.triangles[m].vertex[2];
    t1=mesh.vertices[n1].lon;
    t2=mesh.vertices[n2].lon;
    t3=mesh.vertices[n3].lon;
    t2=degree_recale(t2,t1);
    t3=degree_recale(t3,t1);
    p1=mesh.vertices[n1].lat;
    p2=mesh.vertices[n2].lat;
    p3=mesh.vertices[n3].lat;
    tinf=min(t1,t2);
    updatemin(&tinf,t3);
    tsup=max(t1,t2);
    updatemax(&tsup,t3);
    pinf=min(p1,p2);
    updatemin(&pinf,p3);
    psup=max(p1,p2);
    updatemax(&psup,p3);
    t=degree_recale(t,0.5*(tinf+tsup));
    c1=(t-tinf)*(t-tsup);
    c2=(p-pinf)*(p-psup);
    eligible=((c1 <= 0.0) && (c2 <= 0.0));
    if (!eligible) continue;
    inside=fe_testelement_v2(mesh,mesh.triangles[m],t,p,mode);
    if (inside == 1) {
//      fe_lastelement=m;
      return(m);
      break;
      }
    }

  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_whichelement(const mesh_t & mesh,double t,double p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=fe_whichelement(mesh,t,p,-1);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class T> int fe_whichelement_inlist_template(const mesh_t & mesh,const T & list,int nlisted,double t,double p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,i,n1,n2,n3;
  double t1,t2,t3,p1,p2,p3,pinf,psup,tinf,tsup;
  double c1,c2;
  char  eligible;
  int mode=0; /*Thierry le 03/06/2004*/

  i=fe_lastelement;
  n1=mesh.triangles[i].vertex[0];
  n2=mesh.triangles[i].vertex[1];
  n3=mesh.triangles[i].vertex[2];
  t1=mesh.vertices[n1].lon;
  t2=mesh.vertices[n2].lon;
  t3=mesh.vertices[n3].lon;
  t2=degree_recale(t2,t1);
  t3=degree_recale(t3,t1);
  p1=mesh.vertices[n1].lat;
  p2=mesh.vertices[n2].lat;
  p3=mesh.vertices[n3].lat;
  tinf=min(t1,t2);
  updatemin(&tinf,t3);
  tsup=max(t1,t2);
  updatemax(&tsup,t3);
  pinf=min(p1,p2);
  updatemin(&pinf,p3);
  psup=max(p1,p2);
  updatemax(&psup,p3);
  t=degree_recale(t,0.5*(tinf+tsup));
  c1=(t-tinf)*(t-tsup);
  c2=(p-pinf)*(p-psup);
  eligible=((c1 <= 0.0) && (c2 <= 0.0));
  if (fe_testelement_v1(mesh,i,t,p,mode) == 1) {
    fe_lastelement=i;
    return(i);
    }
/*  printf("#fe_whichelement_inlist %d %d \n",fe_lastelement,nlisted);*/
  for (k=0;k<nlisted;k++) {
    i=list[k];
    n1=mesh.triangles[i].vertex[0];
    n2=mesh.triangles[i].vertex[1];
    n3=mesh.triangles[i].vertex[2];
    t1=mesh.vertices[n1].lon;
    t2=mesh.vertices[n2].lon;
    t3=mesh.vertices[n3].lon;
    t2=degree_recale(t2,t1);
    t3=degree_recale(t3,t1);
    p1=mesh.vertices[n1].lat;
    p2=mesh.vertices[n2].lat;
    p3=mesh.vertices[n3].lat;
    tinf=min(t1,t2);
    updatemin(&tinf,t3);
    tsup=max(t1,t2);
    updatemax(&tsup,t3);
    pinf=min(p1,p2);
    updatemin(&pinf,p3);
    psup=max(p1,p2);
    updatemax(&psup,p3);
    t=degree_recale(t,0.5*(tinf+tsup));
    c1=(t-tinf)*(t-tsup);
    c2=(p-pinf)*(p-psup);
    eligible=((c1 <= 0.0) && (c2 <= 0.0));
    if (!eligible) continue;
    if (fe_testelement_v1(mesh,i,t,p,mode) == 1) {
      fe_lastelement=i;
      return(i);
      }
    }
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_whichelement_inlist(const mesh_t & mesh,const int *list,int nlisted,double t,double p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=fe_whichelement_inlist_template(mesh,list,nlisted,t,p);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_whichelement_inlist(const mesh_t & mesh,const vector<int> & list,double t,double p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=fe_whichelement_inlist_template(mesh,list,list.size(),t,p);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_whichelement_inlist(const mesh_t & mesh,const vector<int> & list,int nlisted,double t,double p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
///only for fe_beta_inlist_template()
{
  int status;
  
  status=fe_whichelement_inlist_template(mesh,list,list.size(),t,p);
  
  return status;
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// void list_t::init(int gnx,int gny)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// ///initialise list for fe_whichelement_createlist()
// {
//   int k;
//   
//   nx=gnx;
//   ny=gny;
//   
//   elements  =new vector<int> *[nx];
//   locks     =new omp_lock_t  *[nx];
//   exitIfNull(elements);
//     
//   for (k=0;k<nx;k++) {
//     elements[k]=new vector<int> [ny];
//     exitIfNull(elements[k]);
//     locks[k]   =new omp_lock_t  [ny];
//     exitIfNull(locks[k]);
//     for(int l=0;l<ny;l++){
//       elements[k][l].reserve(4);
//       omp_init_lock(&locks[k][l]);
//       }
//     }
//    
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void fe_AllocateList(const mesh_t & mesh,grid_t *grid,list_t *list)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///initialise list for fe_whichelement_createlist()
/*----------------------------------------------------------------------------*/
{
  frame_t frame;

  fe_minmax(mesh, frame);
  
  grid->free();
  list->~list_t();
  
  grid->modeH=0;
  
  grid->xmin=frame.xmin;
  grid->ymin=frame.ymin;
  grid->xmax=frame.xmax;
  grid->ymax=frame.ymax;

  grid->xmin-=0.1;
  grid->ymin-=0.1;
  grid->xmax+=0.1;
  grid->ymax+=0.1;

  updatemin(&grid->xmax,grid->xmin+360.);
  
  const int n=16;//minimum is from 12 to 20
  grid->nx=36*n;
  grid->ny=18*n;

  grid->dx=(grid->xmax-grid->xmin)/grid->nx;
  grid->dy=(grid->ymax-grid->ymin)/grid->ny;
  
  list->init(grid->nx,grid->ny);
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// list_t::~list_t()
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /*----------------------------------------------------------------------------*/
// ///free what list_t::init() allocates
// /*----------------------------------------------------------------------------*/
// {
//   deletep2D(&elements,nx);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void fe_PushInList(list_t *list,int k1,int k2,int l1,int l2,int m)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l;
  
  for (k=k1;k<=k2;k++) {
    vector<int> *elementsk=list->elements[k];
    omp_lock_t  *locksk   =list->locks   [k];
    
    for (l=l1;l<=l2;l++) {
      vector<int> *elementskl=&elementsk[l];
      omp_lock_t  *lockskl   =&locksk   [l];
      
      omp_set_lock(lockskl);
      elementskl->push_back(m);
      omp_unset_lock(lockskl);
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void fe_CreateTriangleList(const grid_t & grid,const mesh_t & mesh,list_t *list)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  //#pragma omp parallel for
  for (int i=0;i<mesh.ntriangles;i++) {
    const triangle_t *triangle=&mesh.triangles[i];
    int    n1,n2,n3;
    double t1,t2,t3,p1,p2,p3;
    range_t<double> pR,tR;
    
    n1=triangle->vertex[0];
    n2=triangle->vertex[1];
    n3=triangle->vertex[2];
    
    t1=mesh.vertices[n1].lon;
    t2=mesh.vertices[n2].lon;
    t3=mesh.vertices[n3].lon;
    t2=degree_recale(t2,t1);
    t3=degree_recale(t3,t1);
    
    p1=mesh.vertices[n1].lat;
    p2=mesh.vertices[n2].lat;
    p3=mesh.vertices[n3].lat;
    
    tR<<t1;
    tR<<t2;
    tR<<t3;
    
    pR<<p1;
    pR<<p2;
    pR<<p3;
    
    const double dt=tR.d();
    if(dt<-180. or 180.<dt){
      /* triangle is polar */
      if(pR.max>0){
        pR<< 90.;
        }
      else{
        pR<<-90.;
        }
      tR.init(grid.xmin,grid.xmax);
      }
    
    fe_UpdateList(grid,list,pR,tR,i);
    }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void fe_CreateQuadrangleList(const grid_t & grid,const mesh_t & mesh,list_t *list)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  #pragma omp parallel for
  for (int i=0;i<mesh.nquadrangles;i++) {
    int    n1,n2,n3,n4;
    double t1,t2,t3,t4,p1,p2,p3,p4;
    range_t<double> pR,tR;
    
    n1=mesh.quadrangles[i].vertex[0];
    n2=mesh.quadrangles[i].vertex[1];
    n3=mesh.quadrangles[i].vertex[2];
    n4=mesh.quadrangles[i].vertex[3];
    
    t1=mesh.vertices[n1].lon;
    t2=mesh.vertices[n2].lon;
    t3=mesh.vertices[n3].lon;
    t4=mesh.vertices[n4].lon;
    t2=degree_recale(t2,t1);
    t3=degree_recale(t3,t1);
    t4=degree_recale(t4,t1);
   
    p1=mesh.vertices[n1].lat;
    p2=mesh.vertices[n2].lat;
    p3=mesh.vertices[n3].lat;
    p4=mesh.vertices[n4].lat;
   
    tR<<t1;
    tR<<t2;
    tR<<t3;
    tR<<t4;
    
    pR<<p1;
    pR<<p2;
    pR<<p3;
    pR<<p4;
    
    fe_UpdateList(grid,list,pR,tR,i);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_UpdateList(const grid_t & grid,list_t *list,const range_t<double> &pR,range_t<double> tR, int i, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k1,k2,l1,l2;
  
  tR.min=map_recale(grid,tR.min);
  tR.max=map_recale(grid,tR.max);
  
  l1=int( max(floor((pR.min-grid.ymin)/grid.dy),0.) );
  l2=int( min(floor((pR.max-grid.ymin)/grid.dy),grid.ny-1.) );

  if(tR.min>tR.max) {
    k1=0;
    k2=int( min(1+floor((tR.max-grid.xmin)/grid.dx),grid.nx-1.) );
    if(verbose>0) STDERR_BASE_LINE("fe_PushInList(,%d,%d,%d,%d,%d)\n",k1,  k2, l1, l2, i);
    fe_PushInList(list, k1,  k2, l1, l2, i);
    k1=int( max(floor((tR.min-grid.xmin)/grid.dx)-1,0.) );
    k2=grid.nx-1;
    if(verbose>0) STDERR_BASE_LINE("fe_PushInList(,%d,%d,%d,%d,%d)\n",k1,  k2, l1, l2, i);
    fe_PushInList(list, k1,  k2, l1, l2, i);
    }
  else {
    k1=int( max(floor((tR.min-grid.xmin)/grid.dx),0.) );
    k2=int( min(floor((tR.max-grid.xmin)/grid.dx),grid.nx-1.) );
    if(verbose>0) STDERR_BASE_LINE("fe_PushInList(,%d,%d,%d,%d,%d)\n",k1,  k2, l1, l2, i);
    fe_PushInList(list, k1,  k2, l1, l2, i);
    }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_Allocate_and_CreateList(const mesh_t & mesh,grid_t *grid, list_t *triangles, list_t *quadrangles)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  create a list of elligible elements for grid nodes.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  
  if(mesh.ntriangles!=0) {
    fe_AllocateList(mesh,grid,triangles);
    fe_CreateTriangleList(*grid,mesh,triangles);
    }
  
  if(mesh.nquadrangles!=0) {
    if(quadrangles==0)TRAP_ERR_RETURN(,1,"WARNING : no quadrangle list provided\n");
    fe_AllocateList(mesh,grid,quadrangles);
    fe_CreateQuadrangleList(*grid,mesh,quadrangles);
    }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename Element> int grid_limits(const mesh_t & mesh,const grid_t & grid,const Element & q,const int nnpe, irange_t & i_range, irange_t & j_range,int mode,int align)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n,status;
  double tmid;
  range_t<double> lon_range, lat_range;
  int nx[nnpe],ny[nnpe];
  double t[nnpe],p[nnpe];
  int margin=2;

  for(k=0;k<nnpe;k++) {
    n=q.vertex[k];
    t[k]=mesh.vertices[n].lon;
    p[k]=mesh.vertices[n].lat;
    }

/*------------------------------------------------------------------------------
   */
  if(mode==0) {
    for(k=1;k<nnpe;k++) {
      t[k]=degree_recale(t[k],t[0]);
      }
    }

  lon_range=range_t<double>(t,nnpe);
  lat_range=range_t<double>(p,nnpe);

  if(mode==0) {
    if(align==0) {
      tmid=map_recale(grid,lon_range.min);
      }
    else {
      tmid=map_recale(grid,lon_range.max);
      }
    for(k=0;k<nnpe;k++) {
      n=q.vertex[k];
      t[k]=degree_recale(mesh.vertices[n].lon,tmid);
      }
    lon_range=range_t<double>(t,nnpe);
    lat_range=range_t<double>(p,nnpe);
    }
  
  for(k=0;k<nnpe;k++) {
    status=map_index(grid,t[k],p[k],&nx[k],&ny[k]);
    if(status!=0) {
      status=map_index(grid,t[k],p[k],&nx[k],&ny[k]);
      }
    }

  i_range=irange_t(nx,nnpe);
  j_range=irange_t(ny,nnpe);
  
  i_range.adjust(margin, 0, grid.nx-1);
  j_range.adjust(margin, 0, grid.ny-1);
  
  if(i_range.min>i_range.max or j_range.min>j_range.max)
    return(-1);
  else
    return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename E> int fe_detect_region(const mesh_t &mesh,const E & q,const grid_t & grid, int *buffer,
                                            size_t & completed,const irange_t & i_range,const irange_t & j_range,int mode,int elt,int nprocs)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Tells which nodes of a structured grid are inside a particular element of 
  an unstructured grid.

  mesh         : unstructured grid
  q            : element
  grid         : structured grid
  *buffer      : array[grid.ny*grid.nx] of element indexes, or -1 if outside or undone
  i_range      : range of the element. See limits()
  j_range      : range of the element. See limits()
  mode         : When 0, make fe_testelement_v2() call geo_recale(), as it should
                 for spherical grids.
  elt          : element index
  nprocs If >1 : parallelise whith OpenMP.
 
 returns 0 or -1 if the element is a flat triangle.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  const size_t all=grid.ny*grid.nx;
  int keepLooping=1,flatTriangle=0,chk;

#pragma omp parallel for if(nprocs>1) reduction(or:flatTriangle)
  for (size_t l=j_range.min;l<j_range.max;l++) {
/*------------------------------------------------------------------------------
    See http://www.thinkingparallel.com/2007/06/29/breaking-out-of-loops-in-openmp/ */
    #pragma omp flush(keepLooping)
    if(!keepLooping)
      continue;
    size_t m;
    int status;
    double t,p;
    size_t completed0=0;
    for(size_t k=i_range.min;k<=i_range.max;k++) {
      m=k+grid.nx*l;
      grid.xy(k, l, t, p);
      if (buffer[m] == -1) {
        status=fe_testelement_v2(mesh,q,t,p,mode);
//         chk=fe_testelement_v3(mesh,q,t,p);
//         if(chk!=status) {
//           chk=fe_testelement_v3(mesh,q,t,p);
//           printf("%s alert : %d \n",__func__,elt);
//           }
        if(status==1) {
          buffer[m]=elt;
          completed0++;
          }
        else if (status != 0) {
          printf("%s alert : %d \n",__func__,elt);
          keepLooping=0;
          flatTriangle|=1;
          break;
          }
        }
      }
#pragma omp critical(completedIsCritical)
    {
    completed+=completed0;
    }
    if(completed == all or flatTriangle){
      keepLooping=0;
#pragma omp flush(keepLooping)
      }
    }

  if(flatTriangle)
    return -1;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_detect_region(const mesh_t &mesh, const triangle_t & q, int element, const grid_t & grid, int *buffer,
                       size_t & completed,const irange_t & i_range,const irange_t & j_range, int nprocs)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Tells which nodes of a structured grid are inside a particular element of 
  an unstructured grid.
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int n1,n2,n3;
  double t1,t2,t3,p1,p2,p3;
  const double epsilon=-1.e-09;
  vector3_t o,u,v,w;

  int flatTriangle=0,chk;
  size_t count=0;

  if(flatTriangle)
    return -1;
  
/*------------------------------------------------------------------------------
  gets the coordinates of the corners of the triangle */
  n1=q.vertex[0];
  n2=q.vertex[1];
  n3=q.vertex[2];

  t1=mesh.vertices[n1].lon;
  t2=mesh.vertices[n2].lon;
  t3=mesh.vertices[n3].lon;

  p1=mesh.vertices[n1].lat;
  p2=mesh.vertices[n2].lat;
  p3=mesh.vertices[n3].lat;
  
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  It calculates angle between point and successive triangle vertices
  
  if all angles positive, point is inside
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  returns unit vectors */
  u=math_polar2cartesian(t1*d2r,p1*d2r);
  v=math_polar2cartesian(t2*d2r,p2*d2r);
  w=math_polar2cartesian(t3*d2r,p3*d2r);
  
#pragma omp parallel for if(nprocs>1) reduction(+:count)
  for (size_t l=j_range.min;l<j_range.max;l++) {
    size_t m;
    double t,p,angle;
    for(size_t k=i_range.min;k<=i_range.max;k++) {
      m=k+grid.nx*l;
      if(buffer[m]!=-1) continue;
      grid.xy(k, l, t, p);
      vector3_t uu,vv,ww;
      o=math_polar2cartesian(t*d2r, p*d2r);
/*------------------------------------------------------------------------------
      operator * is scalar product */
      uu=u-o*(o*u);
      vv=v-o*(o*v);
  
      angle=uu^vv;
      if(angle<0) continue;

      ww=w-o*(o*w);
  
      angle=vv^ww;
      if(angle<0) continue;

      angle=ww^uu;
      if(angle<0) continue;
    
      buffer[m]=element;
      count++;
      }
    }
  
  completed+=count;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_detectregular(const mesh_t & mesh,const grid_t & grid,int *buffer,float *fullness,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  Tells which particular element of an unstructured grid contains which nodes
  of a structured grid

  Efficient when triangle are large compared to grid density
  
  mesh      : unstructured grid
  grid      : structured grid
  *buffer   : array[grid.ny*grid.nx] of element indexes, or -1 if outside or undone
  *fullness : percentage of nodes of the structured grid that are contained by
              an element of the unstructured grid
  mode      : See fe_detect_region() : When 0, make fe_testelement_v2() 
              call geo_recale(), as it should for sperical grids.

 returns always 0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  size_t completed,alldone;

  for(size_t j=0;j<grid.ny;j++)
    for (size_t i=0;i<grid.nx;i++) {
      buffer[i+grid.nx*j]=-1;
      }

  completed=0;
  alldone=grid.Hsize();
  
  int nprocs=initialize_OPENMP(-1);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  For all elements it calls:
    - limits() to get the range of the element
    - then fe_detect_region() 
    
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  size_t milestone=mesh.ntriangles/20;
  
  for (size_t i=0;i<mesh.ntriangles;i++) {
    if(i % milestone ==0) printf("%s : %2d%s done\n", __func__, 100*i/(20*milestone),"%");
    irange_t i_range,j_range;
    status=grid_limits(mesh, grid, mesh.triangles[i], 3, i_range, j_range, mode,0);
    if(status==0)
      status= fe_detect_region(mesh, mesh.triangles[i], grid, buffer, completed, i_range, j_range, mode,i, nprocs);
    if(completed == alldone) goto end;

    if(mode!=0) continue;

    status=grid_limits(mesh, grid, mesh.triangles[i], 3, i_range, j_range, mode, 1);
    if(status==0)
      status= fe_detect_region(mesh, mesh.triangles[i], grid, buffer, completed, i_range, j_range, mode,i, nprocs);
    if(completed == alldone) goto end;
    }

  for (size_t i=0;i<mesh.nquadrangles;i++) {
    irange_t i_range,j_range;
    status=grid_limits(mesh, grid, mesh.quadrangles[i], 4, i_range, j_range, mode,0);
    if(status==0)
      status= fe_detect_region(mesh, mesh.quadrangles[i], grid, buffer, completed, i_range, j_range, mode,i, nprocs);
    if(completed == alldone) goto end;

    if(mode!=0) continue;

    status=grid_limits(mesh, grid, mesh.quadrangles[i], 4, i_range, j_range, mode, 1);
    if(status==0)
      status= fe_detect_region(mesh, mesh.quadrangles[i], grid, buffer, completed, i_range, j_range, mode,i, nprocs);
    if(completed == alldone) goto end;
    }

 end:

  *fullness=100.0*completed/alldone;

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_detectregular_experimental(const mesh_t & mesh,const grid_t & grid,int *buffer,float *fullness,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  Tells which particular element of an unstructured grid contains which nodes
  of a structured grid

  Efficient when triangle are large compared to grid density
  
  mesh      : unstructured grid
  grid      : structured grid
  *buffer   : array[grid.ny*grid.nx] of element indexes, or -1 if outside or undone
  *fullness : percentage of nodes of the structured grid that are contained by
              an element of the unstructured grid
  mode      : See fe_detect_region() : When 0, make fe_testelement_v2() 
              call geo_recale(), as it should for sperical grids.

 returns always 0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  size_t completed,alldone,count=0;

  for(size_t j=0;j<grid.ny;j++)
    for (size_t i=0;i<grid.nx;i++) {
      buffer[i+grid.nx*j]=-1;
      }

  completed=0;
  alldone=grid.Hsize();
  
  int nprocs=max(1, omp_get_num_procs()/2); /// HERE !!!
  
  nprocs=initialize_OPENMP(nprocs, 1);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  For all elements it calls:
  
    - grid_limits() to get the grid index range for the element
    - then call fe_detect_region()
    
  cpu load balance to be optimized
    
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  size_t milestone;
  
  milestone=max(1,mesh.ntriangles/20/nprocs);

// #pragma omp parallel for if(nprocs>1) reduction(+:completed) firstprivate (count)
  int nnprocs=1;  
  milestone=max(1,mesh.ntriangles/20);
  for (size_t i=0;i<mesh.ntriangles;i++) {
    if(count % milestone ==0) {
      if(nnprocs!=1) {
        int p=omp_get_thread_num();
        if(p==0) printf("%s : thread %d %2d%s done\n", __func__, p, 100*count/(20*milestone),"%");
//         printf("%s : thread %d %2d%s done\n", __func__, p, 100*count/(20*milestone),"%");
        }
      else {
        printf("%s : %2d%s done\n", __func__, 100*count/(20*milestone),"%");
        }
      }
    count++;
    irange_t i_range,j_range;
    status=grid_limits(mesh, grid, mesh.triangles[i], 3, i_range, j_range, mode,0);
    if(status==0)
      status=fe_detect_region(mesh, mesh.triangles[i], i, grid, buffer, completed, i_range, j_range, nprocs);
//     if(completed == alldone) goto end;
    }

  for (size_t i=0;i<mesh.nquadrangles;i++) {
    irange_t i_range,j_range;
    status=grid_limits(mesh, grid, mesh.quadrangles[i], 4, i_range, j_range, mode,0);
    if(status==0)
      status=fe_detect_region(mesh, mesh.quadrangles[i], grid, buffer, completed, i_range, j_range, mode,i, nprocs);
//     if(completed == alldone) goto end;

    if(mode!=0) continue;

    status=grid_limits(mesh, grid, mesh.quadrangles[i], 4, i_range, j_range, mode, 1);
    if(status==0)
      status=fe_detect_region(mesh, mesh.quadrangles[i], grid, buffer, completed, i_range, j_range, mode,i, nprocs);
//     if(completed == alldone) goto end;
    }

 end:

  *fullness=100.0*completed/alldone;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_detectregular_safe(const mesh_t & mesh,const grid_t & grid,int *buffer,float *fullness,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j;
  float completed,alldone;
  irange_t i_range,j_range;

  completed=0;
  alldone=grid.ny*grid.nx;

  int nprocs=initialize_OPENMP(-1);
#pragma omp parallel for if(nprocs>1)
  for(j=0;j<grid.ny;j++) {
    int m=0;
    for (int i=0;i<grid.nx;i++) {
      int n,status;
      if(m>=mesh.ntriangles) m=0;
      double t,p;
      n=i+grid.nx*j;
      if(buffer[n]!=-1) continue;
      t=grid.x[n];
      p=grid.y[n];
      status=fe_testelement_v2(mesh,mesh.triangles[m],t,p,mode);
      if(status==1) {
        buffer[n]=m;
//         completed++;
        }
      else buffer[n]=-1;
      for (m=0; m<mesh.ntriangles;m++) {
        status=fe_testelement_v2(mesh,mesh.triangles[m],t,p,mode);
        if(status==1) {
          buffer[n]=m;
          break;
          }
        else if (status != 0) {
          status=fe_testelement_v2(mesh,mesh.triangles[m],t,p,mode);
          printf("#Warning: %d \n",m);
          break;
          }
        }
      }
    }

  for(j=0;j<grid.ny;j++) {
    for (int i=0;i<grid.nx;i++) {
      int n;
      n=i+grid.nx*j;
      if(buffer[n]!=-1) completed++;
      }
    }
  *fullness=100.0*completed/(grid.nx*grid.ny);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class T, class E> int fe_FindElement(const mesh_t & mesh, const E * & elements, const T & list, double t,double p, int & last)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,i,status,polar=-1;
  int mode=0;
  
  if(last>=0){
    i=last;
    status=fe_InExtent(mesh, elements[i], t, p);
    if(status!=0){
      status=fe_testelement_v2(mesh,elements[i],t,p,mode);
      if(status == 1) {
        return(i);
        }
      }
    }

  for (k=0;k<list.size();k++) {
    i=list[k];
    status=fe_InExtent(mesh, elements[i], t, p);
    if (status==0) continue;
    status=fe_testelement_v2(mesh,elements[i],t,p,mode);
    if (status == 1) {
      last=i;
      return(i);
      }
    if (status == 2) {
      polar=i;
      }
    }
  
  if(polar!=-1) last=polar;
  
  return(polar);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

const vector<int> *findInList(const list_t & list,const grid_t & grid,double dx,double dy,double t,double p,double *x,double *y,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  return the list of indexes; if outside: NULL

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int    k,l;
  
  if(list.elements==0)
    return NULL;
  
  *x=map_recale(grid,t);
  *y=p;
  
//   if(*x<grid.xmin || *x>grid.xmax ||  *y<grid.ymin || *y>grid.ymax) {
  if(*x<grid.xmin-dx or *x>grid.xmax+dx or  *y<grid.ymin-dy or *y>grid.ymax+dy) {
    /* outside */
    TRAP_ERR_RETURN(NULL,verbose-1,"(%g;%g) out of ([%g;%g];[%g;%g])\n",*x,*y,grid.xmin,grid.xmax,grid.ymin,grid.ymax);
    }
    
  k=int( floor((*x-grid.xmin)/dx) );
  l=int( floor((*y-grid.ymin)/dy) );
  
  if(k<0) {
    if(verbose) STDERR_BASE_LINE_FUNC("[%d,%d]\n",k,l);
    k=0;
    }
  
  if(l<0) {
    if(verbose) STDERR_BASE_LINE_FUNC("[%d,%d]\n",k,l);
    l=0;
    }
  
   if(k>=grid.nx) {
    if(verbose) STDERR_BASE_LINE_FUNC("[%d,%d]\n",k,l);
    k=grid.nx-1;
    }
  
  if(l>=grid.ny) {
    if(verbose) STDERR_BASE_LINE_FUNC("[%d,%d]\n",k,l);
    l=grid.ny-1;
    }
  
  if(verbose>0) STDERR_BASE_LINE_FUNC("[%d,%d]\n",k,l);
  
  return &list.elements[k][l];
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int fe_detect_template(const mesh_t & mesh, const T* elements, const list_t & list,const grid_t & grid, double t, double p, int & last)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  returns if any: the index of the element containing the point; if none: -1

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double x,y;
  int element;
  const vector<int> *targets;
  
  targets=findInList(list,grid,grid.dx,grid.dy,t,p,&x,&y);
  if(targets==0) {
    return -1;
    }
  if(targets->size()==0) {
    return -1;
    }
//   element=fe_whichelement_inlist(mesh,*eis,x,y);
  element=fe_FindElement(mesh, elements, *targets,x,y,last);
  
  return element;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_detect(const mesh_t & mesh, const triangle_t * elements, const list_t & list,const grid_t & grid, double t, double p, int & last)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int found;
  
  found=fe_detect_template(mesh,elements,list,grid,t,p,last);

  return found;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_detect(const mesh_t & mesh, const quadrangle_t * elements, const list_t & list,const grid_t & grid, double t, double p, int & last)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int found;
  
  found=fe_detect_template(mesh,elements,list,grid,t,p,last);

  return(found);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int *fe_detect_template(const mesh_t & mesh, const T* elements, const list_t & list,const grid_t & grid,const double *t,const double *p,int npositions)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n;
  int *found;
  
  found=new int[npositions];

#define DEBUG_fe_detect 0
#if !DEBUG_fe_detect
  int nprocs=initialize_OPENMP(-1, 0);
#pragma omp parallel if(nprocs>1)
#endif
  {
  int last=-1;
#if !DEBUG_fe_detect
#pragma omp for
#endif
  for(n=0; n<npositions; n++) {
    found[n]=fe_detect(mesh, elements, list, grid, t[n], p[n], last);
    }
  }

  return(found);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *fe_detect(const mesh_t & mesh, const triangle_t * elements, const list_t & list,const grid_t & grid,const double *t,const double *p,int npositions)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int *found;
  
  found=fe_detect_template(mesh,elements,list,grid,t,p,npositions);

  return(found);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *fe_detect(const mesh_t & mesh, const quadrangle_t * elements, const list_t & list,const grid_t & grid,const double *t,const double *p,int npositions)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int *found;
  
  found=fe_detect_template(mesh,elements,list,grid,t,p,npositions);

  return(found);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *fe_detect(const mesh_t & mesh, const double *t,const double *p, int npositions)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int *elements=0;
  list_t triangles, quadrangles;
  grid_t grid;

/*------------------------------------------------------------------------------
  create a list of elligible elements for grid nodes */
  fe_Allocate_and_CreateList(mesh,&grid,&triangles,&quadrangles);
  
/*------------------------------------------------------------------------------
  finalise elements detection by seeking in grid nodes' lists */
  if(triangles.elements!=0 and quadrangles.elements!=0)
    TRAP_ERR_EXIT(ENOEXEC,"%s not coded yet for mixed triangle and quadrangle lists\n",__func__);
  
  if(triangles.elements!=0)   elements=fe_detect(mesh, mesh.triangles,triangles,grid,t,p,npositions);
  
  if(quadrangles.elements!=0) elements=fe_detect(mesh, mesh.quadrangles,quadrangles,grid,t,p,npositions);
  
  return elements;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *fe_scan_elements(const mesh_t & mesh,const grid_t & grid,int mode,int option)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  Tells which particular element of an unstructured grid contains which nodes of
  a structured grid
  
  option 0 is standard option, effivcient when triangle lerge compred to grid cells
  
  option 1 is list-based; more RAM's demanding (~10Mb on a 500k mesh) and >100 times quicker.
  
  option 2 is a patch for horrible grids such as ORCA; we need to get a quick
  and safe procedure
  
  mesh unstructured grid
  grid structured grid
  mode given to fe_detectregular(), then fe_detect_region() : When 0 (default), make fe_testelement_v2() call geo_recale(), as it should for sperical grids.
  option When 0 (default), uses fe_detectregular() and takes forever. When 1, uses fe_detect() and only takes ages.

  returns an array[grid.ny*grid.nx] of integers that are equal to the element indexes or -1 if outside

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k, status,*buffer;
  float fullness;
  grid_t tmp;
  double *t,*p;
  
  const int n=grid.Hsize();
  if(n==0)
    return NULL;

  switch(grid.modeH) {
    case 0:
    case 1:
      tmp.duplicate(grid);
      status=map_completegridaxis(&tmp,2);
      t=tmp.x;
      p=tmp.y;
      break;
    default:
      t=grid.x;
      p=grid.y;
      break;
    }
  
  switch(option) {
    case 0:
/*------------------------------------------------------------------------------
      standard algorithm */
      exitIfNull(
        buffer=new int[n]
        );
      for(k=0;k<n;k++) buffer[k]=-1;
      status=fe_detectregular(mesh,grid,buffer,&fullness,mode);
//       status=fe_detectregular_safe(mesh,grid,buffer,&fullness,mode);
      break;
  
    case 1:
/*------------------------------------------------------------------------------
      list-based algorithm */
      buffer=fe_detect(mesh, t, p, n);
      break;
      
    case 2:
/*------------------------------------------------------------------------------
      ORCA stuff */
      exitIfNull(
        buffer=new int[n]
        );
      for(k=0;k<n;k++) buffer[k]=-1;
      status=fe_detectregular(mesh,grid,buffer,&fullness,mode);
      status=fe_detectregular_safe(mesh,grid,buffer,&fullness,mode);
      break;
  
    default:
      buffer=0;
      break;
    }

  tmp.free();

  return(buffer);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int *fe_scan_elements3d(mesh_t mesh,grid_t grid3d, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, status,*buffer;
  float fullness;
  grid_t grid; // a commenter thierry le 26/10/2005

//  buffer=(int *) malloc(grid.nx*grid.ny*sizeof(int));  thierry le 26/10/2005
  buffer=(int *) malloc(grid3d.nx*grid3d.ny*sizeof(int));
  if(buffer==NULL) goto end;

  for(k=0;k<grid3d.nx*grid3d.ny;k++) buffer[k]=0;
//  for(k=0;k<grid.nx*grid.ny;k++) buffer[k]=0;  thierry le 26/10/2005

//  status=fe_detectregular(mesh,grid3d,buffer,&fullness,mode);  A  faire thierry le 26/10/2005
  status=fe_detectregular(mesh,grid,buffer,&fullness,mode);
  if(status !=0) {
    free(buffer);
    buffer=NULL;
    goto end;
    }

end:
  return(buffer);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double rufdist(double lon1,double lat1,double lon2,double lat2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///ruf distance
/*----------------------------------------------------------------------------*/
{
  double r,lat=.5*(lat1+lat2);
  
  r=manhattan(cos(lat)*(lon1-lon2),lat1-lat2);
  
  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double rufdist(const vertex_t & vertex,double lon,double lat)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///ruf distance
/*----------------------------------------------------------------------------*/
{
  double r;
  
  r=rufdist(vertex.lon,vertex.lat,lon,lat);
  
  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fe_nearest_element(const mesh_t & mesh,double lon,double lat,int *m,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Gives element index
/**
\param lon longitude in degrees
\param lat latitude in degrees
\param[in,out] *m index of the closest element. Used as prior solution.
\returns 0 when inside, -d when outside, with d the ruf distance
*/
/*----------------------------------------------------------------------------*/
{
  const int nt=mesh.ntriangles;
  int n;//<node index
  const triangle_t *t;
  const vertex_t *v;
  int i,nmin,status;
  double l,lmin;
  
  if(*m<0 || *m>nt){
    *m=0;
    }
  
  t=&mesh.triangles[*m];
  status=fe_testelement_v2(mesh,*t,lon,lat);
  if(status==1)return 0.;
  
  lmin=+INFINITY;
  for(i=0;i<3;i++){
    n=t->vertex[i];
    l=rufdist(mesh.vertices[n],lon,lat);
    if(lmin>l){
      lmin=l;
      nmin=n;
      }
    }
  
  do{
    v=&mesh.vertices[nmin];
    for(i=0;i<v->nngh;i++){
      n=v->ngh[i];
      l=rufdist(mesh.vertices[n],lon,lat);
      if(lmin>l){
        lmin=l;
        nmin=n;
        }
      }
    }while(nmin!=n);
  
  lmin=+INFINITY;
  for(i=0;i<v->nelmts;i++){
    n=v->elmts[i];
    status=fe_testelement_v2(mesh,mesh.triangles[n],lon,lat);
    if(status==1){
      *m=n;
      return 0.;
      }
    }
  *m=n;
  
  switch(status){
  case 0:
    return -lmin;
  case -1:
    return NAN;
  default:
    TRAP_ERR_EXIT(ENOEXEC,"fe_testelement_v2() bad error %d\n",status);
    }
}
