
/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief map interpolation definitions
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>

#include "tools-structures.h"

#include "maths.h"
#include "geo.h"
#include "constants.h"
#include "map.def"
#include "map.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float map_getf_xyt(const grid_t &grid, float *buf1, float *buf2, double t1, double t2, float mask,
                                float x, float y, double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,k0,l0;
  float w,r1,r2,s,xx,yy,teta;
  float z1,z2,z;
  int   n;

  n=grid.nx;

  teta=(t-t1)/(t2-t1);
/*
  if(teta < 0.0) {
    z=mask;
    return (z);
  }

  if(teta > 1.0) {
    z=mask;
    return (z);
  }
*/
  if(x < grid.xmin) {
    z=mask;
    return (z);
    }
    
  if (x > grid.xmax) {
    z=mask;
    return (z);
    }
 
  if(y < grid.ymin) {
    z=mask;
    return (z);
    }
 
  if(y > grid.ymax) {
    z=mask;
    return (z);
    }

  k0= (int) ((x-grid.xmin)/grid.dx);
  l0= (int) ((y-grid.ymin)/grid.dy);
 
  z1=0.0;
  z2=0.0;
  w=0.0;

  for( k=k0; k <= k0+1; k++) {
    xx=grid.xmin+k*grid.dx;
    if (k < 0)       break;
    if (k > grid.nx-1) break;
    for( l=l0; l <= l0+1; l++) {
      yy=grid.ymin+l*grid.dy;
      if(l < 0)       break;
      if(l > grid.ny-1) break;
      s=(grid.dx-fabs(x-xx))*(grid.dy-fabs(y-yy));
      r1=buf1[k+n*l];
      r2=buf2[k+n*l];
      if (r1 != mask) {
        z1 +=r1*s;
        z2 +=r2*s;
        w=w+s;
        }
      }
    }

if (w > 0.0) {
  z1 /=w;
  z2 /=w;
  z=(1.-teta)*z1+teta*z2;
  }
else
  z=mask;

return(z);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T,typename O> int map_interpolation00(const grid_t & grid,T *buf,T mask,double x,double y,O *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,k0,l0,status=0;
  double xx,yy;
  double w,s;
  T r;

  if(x < grid.xmin) {
    status=-1;
    goto error;
    }

  if (x > grid.xmax) {
    status=-1;
    goto error;
    }

  if(y < grid.ymin) {
    status=-1;
    goto error;
    }

  if(y > grid.ymax+0.5*grid.dy) {
    status=-1;
    goto error;
    }

  if((grid.xmax-grid.xmin)/grid.dx==grid.nx) {
/**----------------------------------------------------------------------------
    cell grid */
    k0= (int) ((x-grid.xmin)/grid.dx);
    l0= (int) ((y-grid.ymin)/grid.dy);
    *z=buf[k0+grid.nx*l0];;
    return(status);
    }
    
  k0= (int) ((x-grid.xmin)/grid.dx);
  l0= (int) ((y-grid.ymin)/grid.dy);

  *z=0.0;
  w=0.0;

  for( k=k0; k <= k0+1; k++) {
    xx=grid.xmin+k*grid.dx;
    if (k < 0)         continue;
    if (k > grid.nx-1) continue;
    for( l=l0; l <= l0+1; l++) {
      yy=grid.ymin+l*grid.dy;
      if(l < 0)         continue;
      if(l > grid.ny-1) continue;
      s=(grid.dx-fabs(x-xx))*(grid.dy-fabs(y-yy));
      r=buf[k+grid.nx*l];
      if (r != mask) {
        *z=*z+r*(O) s;
        w=w+s;
        }
      }
    }

  if (w > 0.0) *z /=w;
  else {
    *z=mask;
    }

  return(0);

 error:
  *z=mask;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T,typename O> int map_interpolation01(const grid_t & grid,T *buf,T mask,double x,double y,O *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,k0,l0,status;
//  double w,r,s,xx,yy;
  double w,s,xx,yy;
  T r;
  int attempt,kprev,lprev;
  double ux,uy,un,dx,dy;

  if(x < grid.xmin-0.5*grid.dy) {
    status=-1;
    goto error;
    }

  if (x > grid.xmax+0.5*grid.dy) {
    status=-1;
    goto error;
    }

  if(y < grid.ymin-0.5*grid.dy) {
    status=-1;
    goto error;
    }

  if(y > grid.ymax+0.5*grid.dy) {
    status=-1;
    goto error;
    }

  k0=grid.nx/2;
  l0=grid.ny/2;

  for(attempt=1;attempt<10;attempt++) {
    kprev=k0;
    lprev=l0;

    ux= grid.x[k0+1]-grid.x[k0];
    un=sqrt(ux*ux);
    ux=ux/un;
    k0= k0+int(floor((x-grid.x[k0])*ux/un));
    if(k0<0) {
      k0=0;
      }
    if(k0>=grid.nx-1) {
      k0=grid.nx-2;
      }
    uy= grid.y[l0+1]-grid.y[l0];
    un=sqrt(uy*uy);
    uy=uy/un;
    l0= l0+int( floor((y-grid.y[l0])*uy/un));
    if(l0<0) {
      l0=0;
      }
    if(l0>=grid.ny-1) {
      l0=grid.ny-2;
      }
    if((kprev==k0) && (lprev==l0)) break;
    }

  *z=0.0;
  w=0.0;

  dx=fabs(grid.x[k0+1]-grid.x[k0]);
  dy=fabs(grid.y[l0+1]-grid.y[l0]);
  for( k=k0; k <= k0+1; k++) {
    xx=grid.x[k];
    if (k < 0)         break;
    if (k > grid.nx-1) break;
    for( l=l0; l <= l0+1; l++) {
      yy=grid.y[l];
      if(l < 0)         break;
      if(l > grid.ny-1) break;
      s=(dx-fabs(x-xx))*(dy-fabs(y-yy));
      size_t m=grid.Hindex(k,l);
      r=buf[m];
      if (r != mask) {
        *z=*z+r*(O)s;
        w=w+s;
        }
      }
    }

  if (w > 0.0) *z /=w;
  else {
    *z=mask;
    return(-1);
    }

  return(0);

 error:
  *z=mask;
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T,typename O> int map_interpolation02(const grid_t & grid,T *dum,T mask,double x,double y,O *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k0,l0;
  double w,area[4],dx,dy;
  T zz;
  T value[4];
  double ux,uy,un,vx,vy;
  double xx[4],yy[4];
  int i,j,m,count;
  size_t m0;
  int attempt,kprev,lprev,status=-1;
  int redo_k,redo_l;
  double a,b,d;
  double epsilon=1.e-02;
  const int nnext=25;
  int nextk[nnext]={0, -1,+1, 0, 0,-1,+1,-1,+1,-2,+2, 0, 0,-2,+2,-2,+2,-1,-1,+1,+1,-2,-2,+2,+2};
  int nextl[nnext]={0,  0, 0,-1,+1,-1,+1,+1,-1, 0, 0,-2,+2,-1,-1,+1,+1,-2,+2,-2,+2,-2,+2,-2,+2};
  
  if(x < grid.xmin-epsilon*grid.dx) {
    *z=mask;
    return (-1);
    }

  if(x > grid.xmax+epsilon*grid.dx) {
    *z=mask;
    return (-1);
    }

  if(y < grid.ymin-epsilon*grid.dy) {
    *z=mask;
    return (-1);
    }

  if(y > grid.ymax+epsilon*grid.dy) {
    *z=mask;
    return (-1);
    }

  k0=grid.nx/2;
  l0=grid.ny/2;

  for(attempt=1;attempt<20;attempt++) {
    kprev=k0;
    lprev=l0;
    m0=grid.Hindex(k0,l0);
    ux= grid.x[m0+1]-grid.x[m0];
    uy= grid.y[m0+1]-grid.y[m0];
    un=sqrt(ux*ux+uy*uy);
    ux=ux/un;
    uy=uy/un;
    // k0= k0+int( ((x-grid.x[m0])*ux+(y-grid.y[m0])*uy)/un ); LR: 01/09/06
    k0= k0+int( floor( ((x-grid.x[m0])*ux+(y-grid.y[m0])*uy)/un) );
    if(k0<0) {
      k0=0;
      }
    if(k0>=grid.nx-1) {
      k0=grid.nx-2;
      }
    ux= grid.x[grid.nx+m0]-grid.x[m0];
    uy= grid.y[grid.nx+m0]-grid.y[m0];
    un=sqrt(ux*ux+uy*uy);
    ux=ux/un;
    uy=uy/un;
    // l0= l0+int( ((x-grid.x[m0])*ux+(y-grid.y[m0])*uy)/un ); LR: 01/09/06
    l0= l0+int( floor( ((x-grid.x[m0])*ux+(y-grid.y[m0])*uy)/un) );
    if(l0<0) {
      l0=0;
      }
    if(l0>=grid.ny-1) {
      l0=grid.ny-2;
      }
    if((kprev==k0) && (lprev==l0)) break;
    }

  redo_k=1;
  redo_l=1;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  we have a starting position

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

 redo:

  m0=grid.Hindex(k0,l0);

  ux= grid.x[m0+1]-grid.x[m0];
  uy= grid.y[m0+1]-grid.y[m0];

  vx= grid.x[grid.nx+m0]-grid.x[m0];
  vy= grid.y[grid.nx+m0]-grid.y[m0];

  dx=x-grid.x[m0];
  dy=y-grid.y[m0];

  d=ux*vy-uy*vx;

  a=(dx*vy-dy*vx)/d;
  b=(ux*dy-uy*dx)/d;
  
  for(int mm=0;mm<9;mm++) {
    int kk=k0+nextk[mm];
    int ll=l0+nextl[mm];
    if(kk<0) continue;
    if(ll<0) continue;
    if(kk>grid.nx-1) continue;
    if(ll>grid.ny-1) continue;
    m0=grid.Hindex(kk,ll);
    ux= grid.x[m0+1]-grid.x[m0];
    uy= grid.y[m0+1]-grid.y[m0];

    vx= grid.x[grid.nx+m0]-grid.x[m0];
    vy= grid.y[grid.nx+m0]-grid.y[m0];

    dx=x-grid.x[m0];
    dy=y-grid.y[m0];

    d=ux*vy-uy*vx;

    a=(dx*vy-dy*vx)/d;
    b=(ux*dy-uy*dx)/d;
    if(a > -1.0 and a < 2.0 and b>-1.0 and b< 2.0) {
      k0=kk;
      l0=ll;
      break;
      }
    }

  if((a < -1.0) || (a > 2.0)) {
    *z=mask;
    return (status);
    }

  if((b < -1.0) || (b > 2.0)) {
    *z=mask;
    return(status);
    }

  if(a < -epsilon) {
    if((k0>0) && (redo_k==1)) {
      k0=k0-1;
      redo_k=0;
      goto redo;
      }
    else {
      *z=mask;
      return (status);
      }
    }

  if(a > 1.+epsilon) {
    if((k0<grid.nx-2) && (redo_k==1)) {
     k0=k0+1;
      redo_k=0;
      goto redo;
      }
    else {
      *z=mask;
      return(status);
      }
    }

  if(b < -epsilon) {
    if((l0>0) && (redo_l==1)) {
      l0=l0-1;
      redo_l=0;
      goto redo;
      }
    else {
      *z=mask;
      return(status);
      }
    }

  if(b > 1.+epsilon) {
    if((l0<grid.ny-2) && (redo_l==1)) {
      l0=l0+1;
      redo_l=0;
      goto redo;
      }
    else {
      *z=mask;
      return(status);
      }
    }

  if( (k0<0) || (k0>=grid.nx-1) || (l0<0) || (l0>=grid.ny-1)) {
    *z=mask;
    return(status);
    }

  status=0;

/*-----------------------------------------------------------------------
  Interpolation (bilinear)*/

  count=0;

  xx[0]=1.-a;
  xx[1]=a;

  yy[0]=1.-b;
  yy[1]=b;

  for(i=0;i<2;i++) {
    for(j=0;j<2;j++) {
      size_t mm=grid.Hindex(k0+i,l0+j);
      zz=dum[mm];
      if (zz != mask) {
        value[count]=zz;
        area[count]=xx[i]*yy[j];
        count++;
        }
      }
    }

  w=0.0;
  for(m=0;m<count;m++) {
    w=w+area[m];
    }

  zz=0;
  if (w > 0.1) {
    for(m=0;m<count;m++) zz+=value[m]*(T)(area[m]/w); /// HERE !!!
    *z=zz;
    return(0);
    }
  else {
    *z=mask;
    return(-1);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_coeffcients02(const grid_t &grid, double x, double y, double *area, int *node, int *count)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k0,l0;
  double w,dx,dy;
  double ux,uy,un,vx,vy;
  double xx[4],yy[4];
  int i,j,m,m0;
  int attempt,kprev,lprev,status=-1;
  int redo_k,redo_l;
  double a,b,d;

  *count=0;

  if(x < grid.xmin) {
    return (status);
    }
    
  if (x > grid.xmax) {
    return (status);
    }
 
  if(y < grid.ymin) {
    return (status);
    }
 
  if(y > grid.ymax) {
    return (status);
    }

  k0=grid.nx/2;
  l0=grid.ny/2;

  for(attempt=1;attempt<10;attempt++) {
    kprev=k0;
    lprev=l0;
    m0=grid.nx*l0+k0;
    ux= grid.x[m0+1]-grid.x[m0];
    uy= grid.y[m0+1]-grid.y[m0];
    un=sqrt(ux*ux+uy*uy);
    ux=ux/un;
    uy=uy/un;
    // k0= k0+int( ((x-grid.x[m0])*ux+(y-grid.y[m0])*uy)/un ); LR: 01/09/06
    k0= k0+int( floor( ((x-grid.x[m0])*ux+(y-grid.y[m0])*uy)/un) );
    if(k0<0) {
      k0=0;
      }
    if(k0>=grid.nx-1) {
      k0=grid.nx-2;
      }
    m0=grid.nx*l0+k0;
    ux= grid.x[grid.nx+m0]-grid.x[m0];
    uy= grid.y[grid.nx+m0]-grid.y[m0];
    un=sqrt(ux*ux+uy*uy);
    ux=ux/un;
    uy=uy/un;

    // l0= l0+int( ((x-grid.x[m0])*ux+(y-grid.y[m0])*uy)/un ); LR: 01/09/06
    l0= l0+int( floor( ((x-grid.x[m0])*ux+(y-grid.y[m0])*uy)/un) );
    if(l0<0) {
      l0=0;
      }
    if(l0>=grid.ny-1) {
      l0=grid.ny-2;
      }
    if((kprev==k0) && (lprev==l0)) break;
    }

  redo_k=1;
  redo_l=1;

 redo:

  m0=grid.nx*l0+k0;

  ux= grid.x[m0+1]-grid.x[m0];
  uy= grid.y[m0+1]-grid.y[m0];

  vx= grid.x[grid.nx+m0]-grid.x[m0];
  vy= grid.y[grid.nx+m0]-grid.y[m0];

  dx=x-grid.x[m0];
  dy=y-grid.y[m0];

  d=ux*vy-uy*vx;

  a=(dx*vy-dy*vx)/d;
  b=(ux*dy-uy*dx)/d;

  if((a < -1.0) || (a > 2.0)) {
    return (status);
    }

  if((b < -1.0) || (b > 2.0)) {
    return(status);
    }

//   if(a < -1.e-05)
//     {
//     if((k0>0) && (redo_k==1))
//       {
//       k0=k0-1;
//       redo_k=0;
//       goto redo;
//       }
//     else
//       {
//       return (status);
//       }
//     }
  if(a < -1.e-05) {
    if((k0>0) && (redo_k==1)) {
      k0=k0-1;
      redo_k=0;
      goto redo;
      }
    else {
      return (status);
      }
    }

  if(a > 1.) {
    if((k0<grid.nx-1) && (redo_k==1)) {
      k0=k0+1;
      redo_k=0;
      goto redo;
      }
    else
      {
      return(status);
      }
    }

  if(b < -1.e-05) {
    if((l0>0) && (redo_l==1)) {
      l0=l0-1;
      redo_l=0;
      goto redo;
      }
    else
      {
      return(status);
      }
    }


  if(b > 1.) {
    if((l0<grid.ny-1) && (redo_l==1)) {
      l0=l0+1;
      redo_l=0;
      goto redo;
      }
    else
      {
      return(status);
      }
    }

  if( (k0<0) || (k0>=grid.nx-1) || (l0<0) || (l0>=grid.ny-1)) {
    return(status);
    }

  status=0;

/*-----------------------------------------------------------------------
  Interpolation (bilinear)*/

  xx[0]=1.-a;
  xx[1]=a;
    
  yy[0]=1.-b;
  yy[1]=b;

  *count=0;
  for(i=0;i<2;i++) {
    for(j=0;j<2;j++) {
      node[*count]=i+k0+(j+l0)*grid.nx;
      area[*count]=xx[i]*yy[j];
      (*count)++;
      }
    }

  w=0.0;
  for(m=0;m<*count;m++) {
    w=w+area[m];
    }

  if (w > 0.0) {
    for(m=0;m<*count;m++) area[m]/=w;
    return(0);
    }
  else return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_interpolation3d(const grid_t &grid,double x, double y, double z, float *buffer)
/* THANKFULLY UNUSED !!! */
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_index00(const grid_t &grid, double x, double y, int *ktrue, int *ltrue)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k0,l0;
  int status;
  int x_out=0,y_out=0;

  if(x < grid.xmin) {
    x_out=1;
    x=grid.xmin-1;
    }

  if (x > grid.xmax) {
    x_out=1;
    x=grid.xmax+1;
    }

  if(y < grid.ymin) {
    y_out=1;
    y=grid.ymin-1;
    }

  if(y > grid.ymax) {
    y_out=1;
    y=grid.ymax+1;
    }

  k0=int((x-grid.xmin)/grid.dx);
  l0=int((y-grid.ymin)/grid.dy);

  *ktrue=k0;
  *ltrue=l0;

  if((x_out==1) || (y_out==1)) status=-1;
  else                           status= 0;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_index02(const grid_t &grid, double x, double y, int *ktrue, int *ltrue)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
///Gives indexes of a 2D grid
/**
\returns 0 on success, -1 on error, and, err... 1 when I don't know! (see bug note below)

\date 2011-09-26 Damien Allain : started a review
\bug 2011-09-26 Damien Allain : map_index02() IS BADLY UNDOCUMENTED
*/
{
  int k,l,k0,l0;
  double area,area_max=-1.,dx,dy;
  double ux,uy,un,vx,vy,vn;
  double xx[4],yy[4];
  int i,j,m0;
  int attempt,kprev,lprev;
  int redo_k,redo_l;
  double a,b,d;
  double epsilon=1.e01,umax;
  double accuracy=5.e-02;
  bool try_again=false;
  bool twisted=false;

  dx= 0.51*(grid.xmax-grid.xmin)/grid.nx;
  dy= 0.51*(grid.ymax-grid.ymin)/grid.ny;

  *ktrue=0;
  *ltrue=0;

  if(x < grid.xmin) {
    x=grid.xmin;
    twisted=true;
    *ktrue=0;
    }

  if(x > grid.xmax) {
    x=grid.xmax;
    twisted=true;
    *ktrue=grid.nx-1;
    }

  if(y < grid.ymin) {
    y=grid.ymin;
    twisted=true;
    *ltrue=0;
    }

  if(y > grid.ymax) {
    y=grid.ymax;
    twisted=true;
    *ltrue=grid.ny-1;
    }

//   if((x < grid.xmin-dx) || (x > grid.xmax+dx) || (y < grid.ymin-dy) || (y > grid.ymax+dy)) {
//     return -1;
//     }

  if( ( (x < grid.xmin) || (x > grid.xmax) ) && ( (y < grid.ymin) || (y > grid.ymax) ) ) {
    return -1;
    }

/*------------------------------------------------------------------------------
  ruf search */
  k0=grid.nx/2;
  l0=grid.ny/2;

  umax=1.e+10;

  for(double s=0;s<0.9999;s+=1/27.) {
    for(double r=0;r<0.9999;r+=1/27.) {
      k=grid.nx*r;
      l=grid.ny*s;
      m0=grid.nx*l+k;
      ux= x-grid.x[m0];
      uy= y-grid.y[m0];
      un=ux*ux+uy*uy;
      if(un<umax) {
        k0=k;
        l0=l;
        umax=un;
        }
      }
    }

  for(attempt=1;attempt<20;attempt++) {
    kprev=k0;
    lprev=l0;
/* *----------------------------------------------------------------------------
    try point m0-(k0,l0]*/
    m0=grid.nx*l0+k0;
    ux= grid.x[m0+1]-grid.x[m0];
    uy= grid.y[m0+1]-grid.y[m0];
    un=sqrt(ux*ux+uy*uy);
    ux=ux/un;
    uy=uy/un;
    k0= k0+floor(((x-grid.x[m0])*ux+(y-grid.y[m0])*uy)/un);
    if(k0<0) {
      k0=0;
      }
    if(k0>=grid.nx-1) {
      k0=grid.nx-2;
      }

//    m0=grid.nx*l0+k0;
    vx= grid.x[grid.nx+m0]-grid.x[m0];
    vy= grid.y[grid.nx+m0]-grid.y[m0];
    vn=sqrt(vx*vx+vy*vy);
    vx=vx/vn;
    vy=vy/vn;

    l0= l0+floor(((x-grid.x[m0])*vx+(y-grid.y[m0])*vy)/vn);
    if(l0<0) {
      l0=0;
      }
    if(l0>=grid.ny-1) {
      l0=grid.ny-2;
      }
    if((kprev==k0) && (lprev==l0)) goto next;
    }

next:

  redo_k=0;
  redo_l=0;

redo:

  m0=grid.nx*l0+k0;

  ux= grid.x[m0+1]-grid.x[m0];
  uy= grid.y[m0+1]-grid.y[m0];

  vx= grid.x[grid.nx+m0]-grid.x[m0];
  vy= grid.y[grid.nx+m0]-grid.y[m0];

  dx=x-grid.x[m0];
  dy=y-grid.y[m0];

  d=ux*vy-uy*vx;

  a=(dx*vy-dy*vx)/d;
  b=(ux*dy-uy*dx)/d;
  
  try_again=false;

/* *-----------------------------------------------------------------------------
  definetely out of cell, thus grid */
  if((a < -1.0-epsilon) || (a > 2.0+epsilon)) {
    return -1;
    }
  if((b < -1.0-epsilon) || (b > 2.0+epsilon)) {
    return -1;
    }

/*-----------------------------------------------------------------------------
  slightly out of cell, give a second chance */
  if(a < 0.-accuracy) {
    if((k0>0) && (redo_k<10)) {
      k0=k0-1;
      redo_k++;
//       goto redo;
      try_again=true;
      }
    else {
      return 1;
      }
    }

  if(a > 1.+accuracy) {
    if((k0<grid.nx-2) && (redo_k<10)) {
     k0=k0+1;
     redo_k++;
//      goto redo;
     try_again=true;
     }
   else {
     return 1;
     }
   }

/*-----------------------------------------------------------------------------
  slightly out of cell, give a second chance */
  if(b < 0.-accuracy) {
    if((l0>0) && (redo_l<10)) {
      l0=l0-1;
      redo_l++;
//       goto redo;
      try_again=true;
      }
    else {
      return 1;
      }
    }


  if(b > 1.+accuracy) {
    if((l0<grid.ny-2) && (redo_l<10)) {
      l0=l0+1;
      redo_l++;
//       goto redo;
      try_again=true;
      }
    else {
      *ktrue=k0;
      *ltrue=l0+1;
      return 1;
      }
    }
  if(try_again) goto redo;

  if( (k0<0) || (k0>=grid.nx-1) || (l0<0) || (l0>=grid.ny-1)) {
    return 1;
    }

/*-----------------------------------------------------------------------
  */

  xx[0]=1.-a;
  xx[1]=a;

  yy[0]=1.-b;
  yy[1]=b;

  for(i=0;i<2;i++) {
    for(j=0;j<2;j++) {
      area=xx[i]*yy[j];
      if(area>area_max) {
        *ktrue=k0+i;
        *ltrue=l0+j;
        area_max=area;
        }
      }
    }

 if(twisted) 
   return(-1);
 else 
   return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_index(const grid_t &grid, double x, double y, int *ktrue, int *ltrue)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------
  return nearest grid point index */
{
  int status;

  switch (grid.modeH) {
    case 0:
      status=map_index00( grid,  x,  y, ktrue, ltrue);
      break;
    case 1:
      status=map_index00( grid,  x,  y, ktrue, ltrue);
      break;
    case 2:
      status=map_index02( grid,  x,  y, ktrue, ltrue);
      break;
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T,typename O> int map_interpolation_template(const grid_t & grid,T *dum,T mask,double x,double y,O *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  if(grid.projection==0)
  if(x < grid.xmin || x > grid.xmax){
    x=geo_recale(x,(grid.xmin+grid.xmax)/2,180.0);
    }

  switch(grid.modeH) {
    case 0:
      status=map_interpolation00(grid, dum, mask, x, y, z);
      break;

    case 1:
      status=map_interpolation01(grid, dum, mask, x, y, z);
      break;

    case 2:
      status=map_interpolation02(grid, dum, mask, x, y, z);
      break;

    default:
      status=-1;
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_interpolation(const grid_t & grid,double *dum,double mask,double x,double y,double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=map_interpolation_template(grid,dum,mask,x,y,z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_interpolation(const grid_t & grid,float *dum,float mask,double x,double y,float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=map_interpolation_template(grid,dum,mask,x,y,z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_interpolation(const grid_t & grid,signed char *dum,signed char mask,double x,double y,float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=map_interpolation_template(grid,dum,mask,x,y,z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_interpolation(const grid_t & grid,short *dum,short mask,double x,double y,float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=map_interpolation_template(grid,dum,mask,x,y,z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_interpolation(const grid_t & grid,complex<float> *dum,complex<float>  mask,double x,double y,complex<float> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=map_interpolation_template(grid,dum,mask,x,y,z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_interpolation(const grid_t & grid,double *dum,double mask,double x,double y,complex<double> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=map_interpolation_template(grid,dum,mask,x,y,z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_interpolation(const grid_t & grid,complex<double>  *dum,complex<double>  mask,double x,double y,complex<double> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=map_interpolation_template(grid,dum,mask,x,y,z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_nearestvalue00_template(const grid_t &grid, int n, T *buf, T mask, double x, double y, T *z, int maxdepth)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l, k0, l0;
  double w, wmax=1.e+10, s, xx, yy;
  int increment=0;
  T r;

  if(x < grid.xmin - grid.dx) { /*!!! FHL : allow extrapolation */
    *z = mask;
    return (-1);
    }

  if(x > grid.xmax + grid.dx) {
    *z = mask;
    return (-1);
    }

  if(y < grid.ymin - grid.dy) {
    *z = mask;
    return (-1);
    }

  if(y > grid.ymax + grid.dy) {
    *z = mask;
    return (-1);
    }

  k0 = (int) ((x - grid.xmin) / grid.dx);
  l0 = (int) ((y - grid.ymin) / grid.dy);

  *z = 0.0;
  w = 0.0;

process:

  for(k = k0-increment; k <= k0+increment + 1; k++) {
    xx = grid.xmin + k * grid.dx;
    if(k < 0)
      break;
    if(k > grid.nx - 1)
      break;
    for(l = l0-increment; l <= l0+increment + 1; l++) {
      yy = grid.ymin + l * grid.dy;
      if(l < 0)
        break;
      if(l > grid.ny - 1)
        break;
      s = (x-xx)*(x-xx)+(y - yy)*(y - yy);
//       if(fabs(s)<1.e-06) {
//         s=1.e-06;
//         }
      r = buf[k + n * l];
      if(r != mask) {
        if(s < wmax) {
          *z = r;
          wmax = s ;
          }
        }
      }
    }

  if((wmax == 1.e+10) &&(increment<maxdepth)) {
    increment++;
    goto process;
    }

  if(wmax < 1.e+10) {
    return (0);
    }
  else {
    *z = mask;
    return (-1);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_nearestvalue00(const grid_t &grid, int n, unsigned short *buf, unsigned short mask, double x, double y, unsigned short *z, int maxdepth)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status= map_nearestvalue00_template(grid, n, buf, mask, x, y, z, maxdepth);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_nearestvalue00(const grid_t &grid, int n, float *buf, float mask, double x, double y, float *z, int maxdepth)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status= map_nearestvalue00_template(grid, n, buf, mask, x, y, z, maxdepth);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_ClosestVertex(const grid_t &grid, double x, double y, int & m0, double & distance, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,k0,l0;
  double dx,dy;
  double ux,uy,un,vx,vy/*,vn*/;
  int i,j,m;
//   int attempt,kprev,lprev;
  double a,b,d;
  double umax,dmax;
  double dr=1./27.,ds=1./27.;

  dx= 0.51*(grid.xmax-grid.xmin)/grid.nx;
  dy= 0.51*(grid.ymax-grid.ymin)/grid.ny;

  if((x < grid.xmin-dx) || (x > grid.xmax+dx) || (y < grid.ymin-dy) || (y > grid.ymax+dy)) {
    return -1;
    }

/* *----------------------------------------------------------------------------
  initialization*/
  k0=grid.nx/2;
  l0=grid.ny/2;

/* *----------------------------------------------------------------------------
  initialization*/
  umax=1.e+10;

  for(double s=0;s<0.9999;s+=ds) {
    for(double r=0;r<0.9999;r+=dr) {
      k=grid.nx*r;
      l=grid.ny*s;
      m0=grid.nx*l+k;
      ux= x-grid.x[m0];
      uy= y-grid.y[m0];
      un=ux*ux+uy*uy;
      if(un<umax) {
        k0=k;
        l0=l;
        umax=un;
        }
      }
    }

//   for(attempt=1;attempt<20;attempt++) {
//     kprev=k0;
//     lprev=l0;
// /* *----------------------------------------------------------------------------
//     try point m0-(k0,l0]*/
//     m0=grid.nx*l0+k0;
//     ux= grid.x[m0+1]-grid.x[m0];
//     uy= grid.y[m0+1]-grid.y[m0];
//     un=sqrt(ux*ux+uy*uy);
//     ux=ux/un;
//     uy=uy/un;
//     k0= k0+floor(((x-grid.x[m0])*ux+(y-grid.y[m0])*uy)/un);
//     if(k0<0) {
//       k0=0;
//       }
//     if(k0>=grid.nx-1) {
//       k0=grid.nx-2;
//       }
// 
//     vx= grid.x[grid.nx+m0]-grid.x[m0];
//     vy= grid.y[grid.nx+m0]-grid.y[m0];
//     vn=sqrt(vx*vx+vy*vy);
//     vx=vx/vn;
//     vy=vy/vn;
// 
//     l0= l0+floor(((x-grid.x[m0])*vx+(y-grid.y[m0])*vy)/vn);
//     if(l0<0) {
//       l0=0;
//       }
//     if(l0>=grid.ny-1) {
//       l0=grid.ny-2;
//       }
//     if((kprev==k0) && (lprev==l0)) break;
//     }

  m0=grid.nx*l0+k0; 
  
//   if(debug) {
//     printf("rough seek: %d %d %d %lf\n",k0,l0,m0,umax);
//     }
  
  int  depth_k=grid.nx*dr;
  int  depth_l=grid.ny*ds;
  bool try_again=true;

redo:  
  umax=1.e+10;
  double c=cos(grid.y[m0]*M_PI/180.);
  for(j=max(0,l0-depth_l);j<min(grid.ny,l0+depth_l+1);j++) {
    for(i=max(0,k0-depth_k);i<min(grid.nx,k0+depth_k+1);i++) {
      m=grid.nx*j+i;
      ux= (x-grid.x[m])*c;
      uy= y-grid.y[m];
      un=ux*ux+uy*uy;
      if(un<umax) {
        k=i;
        l=j;
        umax=un;
        }
      }
    }
    
  k0=k;
  l0=l;
  dmax=1.e+10;
  for(j=max(0,l0-3);j<min(grid.ny,l0+1+3);j++) {
    for(i=max(0,k0-3);i<min(grid.nx,k0+3+1);i++) {
      m=grid.nx*j+i;
      un=geo_haversin_km(x,y,grid.x[m],grid.y[m]);
      if(un<dmax) {
        k=i;
        l=j;
        dmax=un;
        }
      }
    }
    
  m0=grid.nx*l+k;
  
  distance=geo_haversin_km(x,y,grid.x[m0],grid.y[m0]);
  
  if( (distance>0.01) && try_again){
    try_again=false;
    k0=k;
    l0=l;
    goto redo;
    }
 
  if(debug) {
    printf("sharp seek: k=%d l=%d m=%d u=%lf d=%lf\n",k,l,m0,sqrt(umax),distance);
    printf("lon= %lf   x[k,l]= %lf   x[k+1,l]= %lf\n",x,grid.x[m0],grid.x[m0+1]);
    printf("lat= %lf   y[k,l]= %lf   y[k,l+1]= %lf\n",y,grid.y[m0],grid.y[m0+grid.nx]);
    }
  else
    return 0; /* THE LINES BELOW THIS ONE ARE USELESS AND BUGGY! */
  
  ux= grid.x[m0+1]-grid.x[m0];
  uy= grid.y[m0+1]-grid.y[m0];

  if(grid.nx+m0>=grid.Hsize()) TRAP_ERR_EXIT(ENOEXEC,"programming error\n");
  vx= grid.x[grid.nx+m0]-grid.x[m0];
  vy= grid.y[grid.nx+m0]-grid.y[m0];

  dx=x-grid.x[m0];
  dy=y-grid.y[m0];

  d=ux*vy-uy*vx;

  a=(dx*vy-dy*vx)/d;
  b=(ux*dy-uy*dx)/d;
  
  return 0;
}

