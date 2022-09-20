#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"
#include "map.h"

#define nstat 12

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_earthquake00(grid_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

{
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;
  float *earthquake;

  earthquake=(float *) malloc(grid.nx*grid.ny*sizeof(float));

  Lx=grid.xmax-grid.xmin;
  Ly=grid.ymax-grid.ymin;

  for(j=0;j<grid.ny;j++)
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      x=grid.x[k];
      y=grid.y[k];
      s=x*2.*M_PI/Lx;
      r=y*2.*M_PI/Ly;
      earthquake[k]=8.1*sin(r)*exp(-0.37*r)*sin(s/2.);
      }

  return(earthquake);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_earthquake01(grid_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;
  float *earthquake;

  earthquake=(float *) malloc(grid.nx*grid.ny*sizeof(float));

  Lx=(grid.nx-1)*grid.dx;
  Ly=(grid.ny-1)*grid.dy;

  for(j=0;j<grid.ny;j++)
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      x=i*grid.dx;
      y=j*grid.dy;
      s=x*2.*M_PI/Lx;
      r=y*2.*M_PI/Ly;
      earthquake[k]=8.1*sin(r)*exp(-0.37*r)*sin(s/2.);
      }

  return(earthquake);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_earthquake02(grid_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;
  float *earthquake;

  earthquake=(float *) malloc(grid.nx*grid.ny*sizeof(float));

  Lx=(grid.nx-1)*grid.dx;
  Ly=(grid.ny-1)*grid.dy;

  for(j=0;j<grid.ny;j++)
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      x=i*grid.dx;
      y=j*grid.dy;
      s=x*2.*M_PI/Lx;
      r=y*2.*M_PI/Ly;
      earthquake[k]=(8.1*sin(r)*exp(-0.37*r)*(sin(s/0.3)+sin(s/3.)+sin(s/2.)));
      }

  return(earthquake);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_earthquake03(grid_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;
  float *earthquake;

  earthquake=(float *) malloc(grid.nx*grid.ny*sizeof(float));

  Lx=(grid.nx-1)*grid.dx;
  Ly=(grid.ny-1)*grid.dy;

  for(j=0;j<grid.ny;j++)
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      x=i*grid.dx;
      y=j*grid.dy;
      s=x*2.*M_PI/Lx;
      r=y*2.*M_PI/Ly;
      earthquake[k]=(8.1*sin(r)*exp(-0.37*r)*fabs(cos(2*s/2.)));
      }

  return(earthquake);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_earthquake04(grid_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;
  float *earthquake;

  earthquake=(float *) malloc(grid.nx*grid.ny*sizeof(float));

  Lx=(grid.nx-1)*grid.dx;
  Ly=(grid.ny-1)*grid.dy;

  for(j=0;j<grid.ny;j++)
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      x=i*grid.dx;
      y=j*grid.dy;
      s=x*2.*M_PI/Lx;
      r=y*2.*M_PI/Ly;
      earthquake[k]=8.1*sin(r)*exp(-0.37*r)*sin(s/3.)+fabs(sin(2*s/2.));
      }

  return(earthquake);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_earthquake05(grid_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;
  float *earthquake;

  earthquake=(float *) malloc(grid.nx*grid.ny*sizeof(float));

  Lx=(grid.nx-1)*grid.dx;
  Ly=(grid.ny-1)*grid.dy;

  for(j=0;j<grid.ny;j++)
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      x=i*grid.dx;
      y=j*grid.dy;
      s=(x*2.*M_PI/Lx)-Lx;
      r=y*2.*M_PI/Ly;
      earthquake[k]=(6.1*sin(r)*exp(-0.37*r)*sin(s/3.)+(sin(3*s/2.)));
      }

  return(earthquake);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_earthquake06(grid_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;
  float *earthquake;

  earthquake=(float *) malloc(grid.nx*grid.ny*sizeof(float));

  Lx=(grid.nx-1)*grid.dx;
  Ly=(grid.ny-1)*grid.dy;

  for(j=0;j<grid.ny;j++)
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      x=i*grid.dx;
      y=j*grid.dy;
      s=(x*2.*M_PI/Lx)-1.5*Lx;
      r=(y*2.*M_PI/Ly);
      earthquake[k]=(8.1*sin(r)*exp(-0.37*r)*sin(s/3.)+fabs(sin(3*s/2.)))*sin((s+1.5*Lx)/2.);
      }

  return(earthquake);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_earthquake07(grid_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;
  float *earthquake;

  earthquake=(float *) malloc(grid.nx*grid.ny*sizeof(float));

  Lx=(grid.nx-1)*grid.dx;
  Ly=(grid.ny-1)*grid.dy;

  for(j=0;j<grid.ny;j++)
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      x=i*grid.dx;
      y=j*grid.dy;
      s=(x*2.*M_PI/Lx)-1.5*Lx;
      r=y*2.*M_PI/Ly;
      earthquake[k]=8.1*sin(r)*exp(-0.37*r)*fabs(cos(2*s/2.))+(sin(3*s/2.));
      }

  return(earthquake);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_earthquake08(grid_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;
  float *earthquake;

  earthquake=(float *) malloc(grid.nx*grid.ny*sizeof(float));

  Lx=(grid.nx-1)*grid.dx;
  Ly=(grid.ny-1)*grid.dy;

  for(j=0;j<grid.ny;j++)
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      x=i*grid.dx;
      y=j*grid.dy;
      s=(x*2.*M_PI/Lx)-1.5*Lx;
      r=y*2.*M_PI/Ly;
      earthquake[k]=8.1*sin(r)*exp(-0.37*r)*fabs(cos(2*s/2.))+fabs(sin(3*s/2.));
      }

  return(earthquake);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_earthquake09(grid_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;
  float *earthquake;

  earthquake=(float *) malloc(grid.nx*grid.ny*sizeof(float));

  Lx=(grid.nx-1)*grid.dx;
  Ly=(grid.ny-1)*grid.dy;

  for(j=0;j<grid.ny;j++)
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      x=i*grid.dx;
      y=j*grid.dy;
      s=(x*2.*M_PI/Lx)-1.5*Lx;
      r=y*2.*M_PI/Ly;
      earthquake[k]=8.1*sin(r)*exp(-0.37*r)*sin(s/3.)+fabs(sin(3*s/2.));
      }

  return(earthquake);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_earthquake10(grid_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;
  float *earthquake;

  earthquake=(float *) malloc(grid.nx*grid.ny*sizeof(float));

  Lx=(grid.nx-1)*grid.dx;
  Ly=(grid.ny-1)*grid.dy;

  for(j=0;j<grid.ny;j++)
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      x=i*grid.dx;
      y=j*grid.dy;
      s=(x*2.*M_PI/Lx)-1.5*Lx;
      r=y*2.*M_PI/Ly;
      earthquake[k]=(8.1*sin(r)*exp(-0.37*r)*fabs(cos(2*s/2.))+fabs(sin(3*s/2.)))*fabs(sin((s+1.5*Lx)/2.));
      }

  return(earthquake);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_earthquake11(grid_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;
  float *earthquake;

  earthquake=(float *) malloc(grid.nx*grid.ny*sizeof(float));

  Lx=(grid.nx-1)*grid.dx;
  Ly=(grid.ny-1)*grid.dy;

  for(j=0;j<grid.ny;j++)
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      x=i*grid.dx;
      y=j*grid.dy;
      s=x*2.*M_PI/Lx;
      r=y*2.*M_PI/Ly;
      earthquake[k]=(9.6*sin(r)*exp(-0.4*r)*fabs(cos(2*s/2.)))*sin(s/2.)*sin(s/2.);
      }

  return(earthquake);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_earthquake12(grid_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;
  float *earthquake;

  earthquake=(float *) malloc(grid.nx*grid.ny*sizeof(float));

  Lx=(grid.nx-1)*grid.dx;
  Ly=(grid.ny-1)*grid.dy;

  for(j=0;j<grid.ny;j++)
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      x=i*grid.dx;
      y=j*grid.dy;
      s=(x*2.*M_PI/Lx)-Lx;
      r=y*2.*M_PI/Ly;
      earthquake[k]=(8.1*sin(r)*exp(-0.37*r)*sin(s/3.)+(sin(3*s/2.)))*sin((s+Lx)/2.);
      }

  return(earthquake);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_earthquake14(grid_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;
  float *earthquake;

  earthquake=(float *) malloc(grid.nx*grid.ny*sizeof(float));

  Lx=(grid.nx-1)*grid.dx;
  Ly=(grid.ny-1)*grid.dy;

  for(j=0;j<grid.ny;j++)
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      x=i*grid.dx;
      y=j*grid.dy;
      s=x*2.*M_PI/Lx;
      r=y*2.*M_PI/Ly;
      earthquake[k]=(8.1*sin(r)*exp(-0.37*r)*(sin(s/0.3)+sin(s/3.)+sin(s/2.)))*sin(s/2.);
      }

  return(earthquake);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_earthquake15(grid_t grid)

/* r(x), r(0)=2, r(Lx/2)=1,r(Lx)=2 */
/*(cos(x*M_PI/Lx)+1)*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;
  float *earthquake;

  earthquake=(float *) malloc(grid.nx*grid.ny*sizeof(float));

  Lx=(grid.nx-1)*grid.dx;
  Ly=(grid.ny-1)*grid.dy;

  for(j=0;j<grid.ny;j++)
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      x=i*grid.dx;
      y=j*grid.dy;
      s=x*2.*M_PI/Lx;
      r=(y*2.*M_PI/Ly)*(fabs(cos(x*M_PI/Lx))+1);
      earthquake[k]=(1*sin(r)*exp(-0.37*r))*(3*sin(s/2.)+2*sin(s*5/2.));
      }
 
  return(earthquake);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_earthquake16(grid_t grid)

/* r(x), r(0)=2, r(Lx/2)=1,r(Lx)=2 */
/*((4*x*x/(Lx*Lx))-(4*x/Lx)+2)*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;
  float *earthquake;

  earthquake=(float *) malloc(grid.nx*grid.ny*sizeof(float));

  Lx=(grid.nx-1)*grid.dx;
  Ly=(grid.ny-1)*grid.dy;

  for(j=0;j<grid.ny;j++)
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      x=i*grid.dx;
      y=j*grid.dy;
      s=x*2.*M_PI/Lx;
      r=(y*2.*M_PI/Ly)/*((4*x*x/(Lx*Lx))-(4*x/Lx)+2)*/;
      earthquake[k]=(sin(r)*exp(-0.37*r))*(3*sin(s/2.)+2*sin(s*5/2.));
      }
 /*       earthquake[k]=(1.5*sin(r)*exp(-0.37*r))*(3*sin(s/2.)+2*sin(s*5/2.)); */
/*       } */
  return(earthquake);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_earthquake_generic(grid_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float *earthquake;

  earthquake=set_earthquake05(grid);

  return(earthquake);
}

