
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

\brief
*/
/*----------------------------------------------------------------------------*/

#include "constants.h"

#include "maths.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  matrix2x2_t math_rotation2D_init(double angle)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  matrix2x2_t matrix;
  
  angle*=d2r;
  
  double s=sin(angle);
  double c=cos(angle);
  
  matrix.c[0][0]= c;
  matrix.c[0][1]= s;
  matrix.c[1][0]=-s;
  matrix.c[1][1]= c;
  
  return(matrix);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector2D_t math_rotation2D(const matrix2x2_t & R,const vector2D_t & u)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  vector2D_t v;

  v.x=0;
  v.y=0;

  i=0;
  v.x=R.c[0][i]*u.x+R.c[1][i]*u.y;
  i=1;
  v.y=R.c[0][i]*u.x+R.c[1][i]*u.y;

  return(v);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector3_t math_rotation3D_direct(const matrix3x3_t & R,const vector3_t & u)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  vector3_t v;

  v.x=0;
  v.y=0;
  v.z=0;

  i=0;
  v.x=R.c[0][i]*u.x+R.c[1][i]*u.y+R.c[2][i]*u.z;
  i=1;
  v.y=R.c[0][i]*u.x+R.c[1][i]*u.y+R.c[2][i]*u.z;
  i=2;
  v.z=R.c[0][i]*u.x+R.c[1][i]*u.y+R.c[2][i]*u.z;

  return(v);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector3_t math_rotation3D_inverse(const matrix3x3_t & R,const vector3_t & u)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  vector3_t v;

  v.x=0;
  v.y=0;
  v.z=0;

  i=0;
  v.x=R.c[i][0]*u.x+R.c[i][1]*u.y+R.c[i][2]*u.z;
  i=1;
  v.y=R.c[i][0]*u.x+R.c[i][1]*u.y+R.c[i][2]*u.z;
  i=2;
  v.z=R.c[i][0]*u.x+R.c[i][1]*u.y+R.c[i][2]*u.z;

  return(v);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector3_t math_rotation3D(char operation, const matrix3x3_t & R,const vector3_t & u)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector3_t v;

  switch (operation) {
    case 'D':
      v=math_rotation3D_direct(R, u);
      break;

    case 'I':
      v=math_rotation3D_inverse(R, u);
      break;
    }

  return(v);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int math_rotation3Ddeg(char operation,const matrix3x3_t & R1,const matrix3x3_t & R2, double *x,double *y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector3_t u,v,w;
  int status;
  
  *x*=d2r;
  *y*=d2r;
  u=math_polar2cartesian(*x,*y);
  v=math_rotation3D(operation, R1, u);
  w=math_rotation3D(operation, R2, v);
  status=math_cartesian2polar(w,x,y);
  *x/=d2r;
  *y/=d2r;
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int math_rotation3Ddeg(char operation,const matrix3x3_t & R, double *x,double *y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector3_t u,v;
  int status;
  
  *x*=d2r;
  *y*=d2r;
  u=math_polar2cartesian(*x,*y);
  v=math_rotation3D(operation, R, u);
  status=math_cartesian2polar(v,x,y);
  *x/=d2r;
  *y/=d2r;
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector3_t math_polar2cartesian(double lon, double lat, double radius)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** Converts 3D polar coordinates to vector3_t
\param lon in radians
\param lat in radians
\param radius Default:1.0
\returns vector3_t
*/
/*----------------------------------------------------------------------------*/
{
  vector3_t v;
  const double cl=cos(lat);
  
  v.x=cos(lon)*cl;
  v.y=sin(lon)*cl;
  v.z=sin(lat);
  
  if(radius!=1.)v*=radius;

  return(v);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int math_cartesian2polar(const vector3_t & v,double *lon, double *lat, double *radius)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** Converts vector3_t to 3D polar coordinates
\param v
\param *lon in radians
\param *lat in radians
\param *radius
\returns always 0
*/
/*----------------------------------------------------------------------------*/
{
  double _radius;
  
  _radius=hypot(v);
  *lon=atan2(v.y,v.x);
  *lat=asin(v.z/_radius);
  
  if(radius)
    *radius=_radius;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  matrix3x3_t math_rotation3D_init(const vector3_t & axis, double angle)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *-------------------------------------------------------------------------------------------
  http://fr.wikipedia.org/wiki/Rotation_vectorielle */
  matrix3x3_t M, I, A, B;
  int i,j;


  M.c[0][0]=0;  M.c[1][0]=0;  M.c[2][0]=0;
  M.c[0][1]=0;  M.c[1][1]=0;  M.c[2][1]=0;
  M.c[0][2]=0;  M.c[1][2]=0;  M.c[2][2]=0;

  I.c[0][0]=1;  I.c[1][0]=0;  I.c[2][0]=0;
  I.c[0][1]=0;  I.c[1][1]=1;  I.c[2][1]=0;
  I.c[0][2]=0;  I.c[1][2]=0;  I.c[2][2]=1;

  A.c[0][0]=axis.x*axis.x;  A.c[1][0]=axis.y*axis.x;  A.c[2][0]=axis.z*axis.x;
  A.c[0][1]=axis.x*axis.y;  A.c[1][1]=axis.y*axis.y;  A.c[2][1]=axis.z*axis.y;
  A.c[0][2]=axis.x*axis.z;  A.c[1][2]=axis.y*axis.z;  A.c[2][2]=axis.z*axis.z;

  B.c[0][0]=  0;       B.c[1][0]= -axis.z;  B.c[2][0]= +axis.y;
  B.c[0][1]= +axis.z;  B.c[1][1]=  0;       B.c[2][1]= -axis.x;
  B.c[0][2]= -axis.y;  B.c[1][2]= +axis.x;  B.c[2][2]=  0;

//  M=cos(angle)*I+(1.-cos(angle))*A+sin(angle)*B;
  for(j=0;j<3;j++) {
    for(i=0;i<3;i++) {
      M.c[j][i]=cos(angle)*I.c[j][i]+(1.-cos(angle))*A.c[j][i]+sin(angle)*B.c[j][i];
      }
    }

  return(M);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  matrix3x3_t math_rotation3D_init(char a, double angle)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  matrix3x3_t M;
  vector3_t axis;
  
  switch(a){
  case 'x':case 'X':
    axis.x=1.;
    break;
  case 'y':case 'Y':
    axis.y=1.;
    break;
  case 'z':case 'Z':
    axis.z=1.;
    break;
  default:
    TRAP_ERR_EXIT(ENOEXEC,"No axis %c.\n",a);
    }
  
  M=math_rotation3D_init(axis,angle);

  return M;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int math_orthobase(const point2D_t & p, const point2D_t & q, vector2D_t *u, vector2D_t *v, double *normp)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///give vector basis
/**
\param *normp if not NULL, set to the distance between \c p and \c q
\return 0 on success or -1 if \c *normp is 0.
*/
/*----------------------------------------------------------------------------*/
{
  double norm;

  u->x=q.x-p.x;
  u->y=q.y-p.y;

  norm=hypot(*u);

  if(norm==0.) return(-1);

  *u*=1./norm;

  v->x=-u->y;
  v->y= u->x;
  
  if(normp)
    *normp=norm;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector2D_t math_vector_coordinates01(const vector2D_t & r,const vector2D_t & u,const vector2D_t & v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  coordinates of r in u,v base
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  vector2D_t w;

/*------------------------------------------------------------------------------
  true if u,v is orthonormal */
  w.x=r.x*u.x+r.y*u.y;
  w.y=r.x*v.x+r.y*v.y;

  return(w);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector2D_t math_vector_coordinates02(const vector2D_t & r,const vector2D_t & u,const vector2D_t & v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  coordinates of r in u,v base 
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double det;
  vector2D_t w;

  det=u*v;
  w.x=r*v/det;
  w.y=u*r/det;

  return(w);
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   template <typename T> T math_vector_coordinates03(T a,const vector2D_t & u, T b,const vector2D_t & v)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /* *----------------------------------------------------------------------------
//   coordinates of vector with a projection on u and b rpojection on v
//   u.w = ux * wx + uy * wy = a
//   v.w = vx * wx + vy * wy = b
//   ux uy  wx
//   vx vy  wy
//  -----------------------------------------------------------------------------*/
// {
//   int status;
//   double norm,det;
//   
//   vector2D_t uu(u.x,v.x);
//   vector2D_t vv(u.y,v.y);
// 
//   vector2D_t r(a,b),w;
// 
//   det=u*v;
//   w.x=r*v/det;
//   w.y=u*r/det;
// 
//   return(w);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector2D_t math_vector_coordinates03(double a,const vector2D_t & u, double b,const vector2D_t & v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  coordinates of vector with a projection on u and b projection on v
  u.w = ux * wx + uy * wy = a
  v.w = vx * wx + vy * wy = b
  ux uy  wx
  vx vy  wy
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double det;
  
  vector2D_t uu(u.x,v.x);
  vector2D_t vv(u.y,v.y);

  vector2D_t r(a,b),w;

  det=u*v;
  w.x=r*v/det;
  w.y=u*r/det;

  return(w);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  cvector2D_t math_vector_coordinates03(complex<double> a,const vector2D_t & u, complex<double> b,const vector2D_t & v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  coordinates of vector with a projection on u and b rpojection on v
  u.w = ux * wx + uy * wy = a
  v.w = vx * wx + vy * wy = b
  ux uy  wx
  vx vy  wy
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double det;
  
  vector2D_t uu(u.x,v.x);
  vector2D_t vv(u.y,v.y);

  cvector2D_t r(a,b),w;

  det=uu*vv;
  w.x=r*vv/det;
  w.y=uu*r/det;

  return(w);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double hypot(double x,double y,double z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  3D distance
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double m;

  m=square(x)+square(y)+square(z);
  m=sqrt(m);

  return m;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double manhattan(double x,double y,double z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  3D manhattan distance

  note : not faster than square if not inline

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double m;

  m=fabs(x)+fabs(y)+fabs(z);

  return m;
}

