
/**************************************************************************

  T-UGO tools, 2006-2013

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

\brief mesh area integration and mass matrices functions and element position functions
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <string.h>

#include "tools-structures.h"
#include "constants.h"

#include "fe.def"
#include "fe.h"

#include "map.h"
#include "geo.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_LGP0base(double x, double y, double *beta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  beta[0] = 1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_LGP1base(double x, double y, double *beta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  beta[0] = 1. - x - y;
  beta[1] = x;
  beta[2] = y;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_LGP1prime(double x, double y, double *beta_x,double *beta_y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  beta_x[0] = -1;
  beta_x[1] =  1;
  beta_x[2] =  0;

  beta_y[0] = -1;
  beta_y[1] =  0;
  beta_y[2] =  1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_NCP1base(double x, double y, double *beta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  beta[0] = 2 * x + 2 * y - 1;
  beta[1] = 1 - 2 * x;
  beta[2] = 1 - 2 * y;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_NCP1prime(double x, double y, double *beta_x,double *beta_y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  beta_x[0] =  2;
  beta_x[1] = -2;
  beta_x[2] =  0;

  beta_y[0] =  2;
  beta_y[1] =  0;
  beta_y[2] = -2;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void fe_LGP2base(double x, double y, double *beta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  beta[0] =  2. * (x + y - 0.5) * (x + y - 1.);   /*  2x²+2y²+4xy-3x-3y+1 */
  beta[1] = -4. * x * (x + y - 1.);               /* -4x²    -4xy+4x      */
  beta[2] =  2. * x * (x - 0.5);                  /*  2x²        - x      */
  beta[3] =  4. * x * y;                          /*         +4xy         */
  beta[4] =  2. * y * (y - 0.5);                  /*     +2y²       - y   */
  beta[5] = -4. * y * (x + y - 1.);               /*     -4y²-4xy   +4y   */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void fe_LGP2prime(double x, double y, double *beta_x,double *beta_y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  beta_x[0] =   4*x+4*y-3;   /*  4x+4y-3 */
  beta_x[1] =  -8*x-4*y+4;   /* -8x-4y+4 */
  beta_x[2] =   4*x    -1;   /*  4x   -1 */
  beta_x[3] =      +4*y  ;   /*    +4y   */
  beta_x[4] =           0;   /*	       0 */
  beta_x[5] =      -4*y  ;   /*    -4y   */

  beta_y[0] =   4*x+4*y-3;   /*  4x+4y-3 */
  beta_y[1] =  -4*x      ;   /* -4x      */
  beta_y[2] =           0;   /*        0 */
  beta_y[3] =   4*x      ;   /*  4x      */
  beta_y[4] =      +4*y-1;   /*    +4y-1 */
  beta_y[5] =  -4*x-8*y+4;   /* -4x-8y+4 */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_LGP3base(double x, double y, double *beta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float c1, c2;
  float d1, d2, d3, d4, d5, d6, d7, d8, d9;

  c1 = 1. / 3.;
  c2 = 2. / 3.;

  d1 = x;
  d2 = (x - c1) / c1;
  d3 = (x - c2) / c2;
  d4 = y;
  d5 = (y - c1) / c1;
  d6 = (y - c2) / c2;
  d7 = (x + y - c1) / c1;
  d8 = (x + y - c2) / c2;
  d9 = x + y - 1;

  beta[0] = -d7 * d8 * d9;
  beta[1] = d1 * d8 * d9 / (c1 * c1);
  beta[2] = -d1 * d2 * d9 / (c1 * c2);
  beta[3] = d1 * d2 * d3;
  beta[4] = d1 * d2 * d4 / (c1 * c2);
  beta[5] = d1 * d4 * d5 / (c1 * c2);
  beta[6] = d4 * d5 * d6;
  beta[7] = -d4 * d5 * d9 / (c1 * c2);
  beta[8] = d4 * d8 * d9 / (c1 * c1);
  beta[9] = -d1 * d4 * d9 / (c1 * c1 * c1);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_LGP4base(double x, double y, double *beta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double c1, c2, c3, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12;

  c1 = 0.25;
  c2 = 0.5;
  c3 = 0.75;

  d1 = x;
  d2 = (x - c1) / c1;
  d3 = (x - c2) / c2;
  d4 = (x - c3) / c3;
  d5 = y;
  d6 = (y - c1) / c1;
  d7 = (y - c2) / c2;
  d8 = (y - c3) / c3;
  d9 = (x + y - c1) / c1;
  d10 = (x + y - c2) / c2;
  d11 = (x + y - c3) / c3;
  d12 = x + y - 1;

  beta[0] = d9 * d10 * d11 * d12;
  beta[1] = -d1 * d10 * d11 * d12 / (c1 * c1);
  beta[2] = d1 * d2 * d11 * d12 / (c1 * c1 / c3);
  beta[3] = -d1 * d2 * d3 * d12 / (c1 * c3);
  beta[4] = d1 * d2 * d3 * d4;
  beta[5] = d1 * d2 * d3 * d5 / (c1 * c3);
  beta[6] = d1 * d2 * d5 * d6 / (c2 * c2);
  beta[7] = d1 * d5 * d6 * d7 / (c1 * c3);
  beta[8] = d5 * d6 * d7 * d8;
  beta[9] = -d5 * d6 * d7 * d12 / (c1 * c3);
  beta[10] = d5 * d6 * d11 * d12 / (c1 * c1 / c3);
  beta[11] = -d5 * d10 * d11 * d12 / (c1 * c1);
  beta[12] = d1 * d5 * d11 * d12 / (c1 * c2 * c1 * c1 / c3);
  beta[13] = -d1 * d2 * d5 * d12 / (c1 * c1 * c2);
  beta[14] = -d1 * d5 * d6 * d12 / (c1 * c1 * c2);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_QP0base(double x, double y, double *beta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double c1, c2, c3;

  beta[0] = 1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_QP1base(double x, double y, double *beta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double c1, c2, c3;

  beta[0] = (1. - x) * (1. - y);
  beta[1] = x * (1. - y);
  beta[2] = x * y;
  beta[3] = (1. - x) * y;
}




/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_BL3Dbase(double x, double y, double z, double *beta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double c1, c2, c3;

  beta[0] = (1. - x - y) * (1. - z);
  beta[1] = x * (1. - z);
  beta[2] = y * (1. - z);
  beta[3] = (1. - x - y) * z;
  beta[4] = x * z;
  beta[5] = y * z;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_BL3Dprime(double x, double y, double z, double *prime[3])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double c1, c2, c3;

  prime[0][0] = -(1. - z);
  prime[0][1] = -(1. - z);
  prime[0][2] = -(1. - x - y);

  prime[1][0] = (1. - z);
  prime[1][1] = 0.0;
  prime[1][2] = -x;

  prime[2][0] = 0.0;
  prime[2][1] = (1. - z);
  prime[2][2] = -y;

  prime[3][0] = -z;
  prime[3][1] = -z;
  prime[3][2] = (1. - x - y);

  prime[4][0] = z;
  prime[4][1] = 0.0;
  prime[4][2] = x;

  prime[5][0] = 0.0;
  prime[5][1] = z;
  prime[5][2] = y;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_LGPbase(const mesh_t & mesh,int discretisation,double x,double y,double *beta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, status, nnpe;
  status = 0;
  extern void fe_Q0base(double x, double y, double *beta);
  extern void fe_Q1base(double x, double y, double *beta);

  switch (discretisation) {
    case CQP0:
      nnpe=1;
      fe_QP0base(x, y, beta);
      break;
      
    case CQP1:
      nnpe=1;
      fe_QP1base(x, y, beta);
      break;
      
    case LGP0:
      nnpe=1;
      fe_LGP0base(x, y, beta);
      break;
      
    case LGP1:
    case DGP1:
      nnpe=3;
      fe_LGP1base(x, y, beta);
      break;

    case NCP1:
    case DNP1:
      nnpe=3;
      fe_NCP1base(x, y, beta);
      break;

    case LGP2:
    case DGP2:
      nnpe=6;
      fe_LGP2base(x, y, beta);
      break;

    default:
      status = -1;
      break;
  }
//       for(i = 0; i < nnpe; i++)
//         if((beta[i] - 1.) * beta[i] > 1.e-02) //complicated way of ensuring that 0<=beta[i]<=1
//           status = -1;
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_barycentre(mesh_t mesh,int i,double *t,double *p,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n1,n2,n3,status;
  double t1,t2,t3,p1,p2,p3;

  n1=mesh.triangles[i].vertex[0];
  n2=mesh.triangles[i].vertex[1];
  n3=mesh.triangles[i].vertex[2];

  t1=mesh.vertices[n1].lon;
  t2=mesh.vertices[n2].lon;
  t3=mesh.vertices[n3].lon;

  if(mode==0) {
    t1=geo_recale(t1,t1,(double) 180.0);
    t2=geo_recale(t2,t1,(double) 180.0);
    t3=geo_recale(t3,t1,(double) 180.0);
    }
  p1=mesh.vertices[n1].lat;
  p2=mesh.vertices[n2].lat;
  p3=mesh.vertices[n3].lat;

  *t=(t1+t2+t3)/3;
  *p=(p1+p2+p3)/3;

  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_integrale01(mesh_t mesh, float *buffer, double *sum)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  compute meah domain integrale of LGP1 buffer

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int	 i,j,k,m,status;
  double  C;
  double  z,a;
  double  beta[16][3];
  double mean,variance,rms,A;
  double tmp[10],budget;
  double delta1,delta2,c1,c2;
  double omega, factor=1.e+06*3600.0;
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int gauss_n=3;
  int n, verbose=0;
/*----------------------------------------------------------------------
  initialise numerical integrals*/
  status=gauss_init(gauss_n,gauss_x,gauss_y,gauss_w);
  for(k=0; k<gauss_n; k++) {
    fe_LGP1base(gauss_x[k],gauss_y[k],beta[k]);
    }

  for(j=0;j<2;j++) sum[j]=0.0;

  for(m=0; m<mesh.ntriangles; m++) {
/*----------------------------------------------------------------------
    convert A (element area) in km**2*/
    A=mesh.triangles[m].Area/(1.0e+6);
/*----------------------------------------------------------------------
    spatial integration*/
    tmp[0]=0;
    tmp[1]=0;
    for(k=0; k<gauss_n; k++) {
      z=0;
      a=0;
      for(j=0;j<3;j++) {
        n=mesh.triangles[m].vertex[j];
        C=cos(mesh.vertices[n].lat*d2r);
        z+=buffer[n]*beta[k][j]*C;
        a+=beta[k][j]*C;
        }
      tmp[0]+=gauss_w[k]*z;
      tmp[1]+=gauss_w[k]*a;
      }
    for(j=0;j<2;j++) sum[j]+=2.*A*tmp[j];
    }

  if(verbose==1) {
    printf("integrale (user units x km**2)--: %lf\n",sum[0]);
    printf("surface (10^6 km**2)------------: %lf\n",sum[1]/1.e+06);
    printf("mean (user units)= %lf\n",sum[0]/sum[1]);
   }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_integrale02(mesh_t mesh, float *buffer, double *sum, double latmax)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  compute meah domain integrale of LGP1 buffer within latitude range

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int	 i,j,k,m,status;
  double  C;
  double  z,a;
  double  beta[3];
  double mean,variance,rms,A;
  double tmp[10],budget;
  double delta1,delta2,c1,c2;
  double omega, factor=1.e+06*3600.0;
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int gauss_n=3;
  int n,verbose=0;
/*----------------------------------------------------------------------
  initialise numerical integrals*/
  status=gauss_init(gauss_n,gauss_x,gauss_y,gauss_w);

  for(j=0;j<2;j++) sum[j]=0.0;

  for(m=0; m<mesh.ntriangles; m++) {
/*----------------------------------------------------------------------
    convert A (element area) in km**2*/
    A=mesh.triangles[m].Area/(1.0e+6);
/*----------------------------------------------------------------------
    spatial integration*/
    tmp[0]=0;
    tmp[1]=0;
    for(j=0;j<3;j++) {
      n=mesh.triangles[m].vertex[j];
      if(fabs(mesh.vertices[n].lat)> latmax) goto skip;
      }
    for(k=0; k<gauss_n; k++) {
      fe_LGP1base(gauss_x[k],gauss_y[k],beta);
      z=0;
      a=0;
      for(j=0;j<3;j++) {
        n=mesh.triangles[m].vertex[j];
        C=cos(mesh.vertices[n].lat*d2r);
        z+=buffer[n]*beta[j]*C;
        a+=beta[j]*C;
        }
      tmp[0]+=gauss_w[k]*z;
      tmp[1]+=gauss_w[k]*a;
      }
    skip:
    for(j=0;j<2;j++) sum[j]+=2*A*tmp[j];
    }


  if(verbose==1) {
    printf("integrale (user units x km**2)--: %lf\n",sum[0]);
    printf("surface (10^6 km**2)------------: %lf\n",sum[1]/1.e+06);
    printf("mean (user units)= %lf\n",sum[0]/sum[1]);
    }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_massmatrix01(mesh_t mesh, double **A, double*b, int ngauss)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*----------------------------------------------------------------------------*/
  int	 i,j,k,m,status;
  double  C;
  double  beta[3];
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int gauss_n;
  int ni,nj,neq;
  int *pivot;
/*----------------------------------------------------------------------------*/

  gauss_n=ngauss;
  status=gauss_init(gauss_n,gauss_x,gauss_y,gauss_w);

  for(ni=0; ni<mesh.nvtxs; ni++)
    for(nj=0; nj<mesh.nvtxs; nj++)
       A[ni][nj]=0;

  for(m=0; m<mesh.ntriangles; m++) {
/*----------------------------------------------------------------------

    spatial integration:


    int{ beta-n * sum (a-m x  beta-m)}= int{ beta-n * buffer }

    over a given triangle:

    sum  C*wk * beta-n(k)*sum[a-m x  beta-m(k)] = sum C*wk*beta-n(k)*buffer(k)
     k

----------------------------------------------------------------------*/
    for(k=0; k<gauss_n; k++) {
      fe_LGP1base(gauss_x[k],gauss_y[k],beta);
      for(i=0;i<3;i++) {
        ni=mesh.triangles[m].vertex[i];
        C=cos(mesh.vertices[ni].lat*d2r);
        for(j=0;j<3;j++) {
          nj=mesh.triangles[m].vertex[j];
          A[ni][nj]+=beta[i]*beta[j]*gauss_w[k]*C;
          }
        }
/*
      b[n1]+=gauss_w[k]*buffer[m];
*/
      }
    }
  neq=mesh.nvtxs;
  /* needs specific factorisation*/

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_massmatrix02(mesh_t mesh, double *A, int *pivot, int hbw)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*----------------------------------------------------------------------------*/
  int	 i,j,k,m,status;
  int    gauss_n;
  int    n,neq,bandmax,bandwidth,ml,mu;
  int    ni,nj,row,column;
  double p[3],q[3];
  double C,area;
  double beta[3];
  double gauss_x[16], gauss_y[16], gauss_w[16];
/*----------------------------------------------------------------------------*/

  ml=hbw;
  mu=hbw;

  bandmax=2*ml+mu+1;
  bandwidth=ml+mu+1;

  printf("Enter fe_massmatrix02 \n");

  for(n=0; n<mesh.nvtxs*bandmax; n++) A[n]=0;

  for(m=1; m<=mesh.ntriangles; m++) {
    area=mesh.triangles[m].TrueArea;
    for(i=0;i<3;i++) {
      ni=mesh.triangles[m].vertex[i]-1;
      for(k=0; k<3; k++) q[k]=0;
      q[i]=1.;
      for(j=0;j<3;j++) {
        nj=mesh.triangles[m].vertex[j]-1;
        column=nj;
        row=bandwidth-1+ni-nj;
        for(k=0; k<3; k++) p[k]=0;
        p[j]=1.;
        A[column*bandmax+row]+=2.*area*integraleP1xP1(p,q);
        }
      }
    }

  neq=mesh.nvtxs;
/*
#ifdef LINPACKF_
  dgbfa_ (A,&bandmax,&neq,&hbw,&hbw,pivot,&status);
#elif LINPACKF
  dgbfa (A,&bandmax,&neq,&hbw,&hbw,pivot,&status);
#elif LINPACKC
  dgbfa (A,bandmax,neq,hbw,hbw,pivot,&status);
#elif  LAPACKC
  dgbtrf(bandmax,neq,hbw,hbw,A,bandmax,pivot,&status);
#endif
*/
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_massmatrix03(mesh_t mesh, double *A, int *pivot, int hbw, int ngauss)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*----------------------------------------------------------------------------*/
  int	 i,j,k,m,status;
  int    gauss_n;
  int    n,neq,bandmax,bandwidth,ml,mu;
  int    ni,nj,row,column;
  double C,area;
  double beta[3];
  double gauss_x[16], gauss_y[16], gauss_w[16];
/*----------------------------------------------------------------------------*/

  gauss_n=ngauss;
  status=gauss_init(gauss_n,gauss_x,gauss_y,gauss_w);

  ml=hbw;
  mu=hbw;

  bandmax=2*ml+mu+1;
  bandwidth=ml+mu+1;

  for(n=0; n<mesh.nvtxs*bandmax; n++) A[n]=0;

  for(m=1; m<=mesh.ntriangles; m++) {
    area=mesh.triangles[m].Area;
    for(k=0; k<gauss_n; k++) {
      fe_LGP1base(gauss_x[k],gauss_y[k],beta);
      for(i=0;i<3;i++) {
        ni=mesh.triangles[m].vertex[i]-1;
        C=cos(mesh.vertices[ni].lat*d2r);
        for(j=0;j<3;j++) {
          nj=mesh.triangles[m].vertex[j]-1;
          column=nj;
          row=bandwidth-1+ni-nj;
          A[column*bandmax+row]+=beta[i]*beta[j]*gauss_w[k]*C*area;
          }
        }
      }
    }

  neq=mesh.nvtxs;

/*
#ifdef LINPACKF_
  dgbfa_ (A,&bandmax,&neq,&hbw,&hbw,pivot,&status);
#elif LINPACKF
  dgbfa (A,&bandmax,&neq,&hbw,&hbw,pivot,&status);
#elif LINPACKC
  dgbfa (A,bandmax,neq,hbw,hbw,pivot,&status);
#elif  LAPACKC
  dgbtrf(bandmax,neq,hbw,hbw,A,bandmax,pivot,&status);
#endif
*/
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_position(const mesh_t &mesh, const vertex_t & vertex, double *t,double *p,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*---------------------------------------------------------------------
  should be used to return computational node position (element barycentre, edge
  mid point, vertex position etc....)
!---------------------------------------------------------------------*/
{
  int n1,n2,n3,status;

  *t=vertex.lon;
  *p=vertex.lat;

  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_position(const mesh_t &mesh, const triangle_t & triangle, double *t,double *p,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

/*---------------------------------------------------------------------
  should be used to return computational node position (element barycentre, edge
  mid point, vertex position etc....)
!---------------------------------------------------------------------*/
{
  int n1,n2,n3,status;
  double t1,t2,t3,p1,p2,p3;

  n1=triangle.vertex[0];
  n2=triangle.vertex[1];
  n3=triangle.vertex[2];

  t1=mesh.vertices[n1].lon;
  t2=mesh.vertices[n2].lon;
  t3=mesh.vertices[n3].lon;

  if(mode==0) {
    t1=geo_recale(t1,t1,(double) 180.0);
    t2=geo_recale(t2,t1,(double) 180.0);
    t3=geo_recale(t3,t1,(double) 180.0);
    }
  p1=mesh.vertices[n1].lat;
  p2=mesh.vertices[n2].lat;
  p3=mesh.vertices[n3].lat;

  *t=(t1+t2+t3)/3;
  *p=(p1+p2+p3)/3;
  
  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_position(const mesh_t &mesh, const quadrangle_t & quadrangle, double *t,double *p,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*---------------------------------------------------------------------
  should be used to return computational node position (element barycentre, edge
  mid point, vertex position etc....)
!---------------------------------------------------------------------*/
{
  int n1,n2,n3,n4,status;
  double t1,t2,t3,t4,p1,p2,p3,p4;

  n1=quadrangle.vertex[0];
  n2=quadrangle.vertex[1];
  n3=quadrangle.vertex[2];
  n4=quadrangle.vertex[3];

  t1=mesh.vertices[n1].lon;
  t2=mesh.vertices[n2].lon;
  t3=mesh.vertices[n3].lon;
  t4=mesh.vertices[n4].lon;

  if(mode==0) {
    t1=geo_recale(t1,t1,(double) 180.0);
    t2=geo_recale(t2,t1,(double) 180.0);
    t3=geo_recale(t3,t1,(double) 180.0);
    t4=geo_recale(t4,t1,(double) 180.0);
    }
  p1=mesh.vertices[n1].lat;
  p2=mesh.vertices[n2].lat;
  p3=mesh.vertices[n3].lat;
  p4=mesh.vertices[n4].lat;

  *t=(t1+t2+t3+t4)/4;
  *p=(p1+p2+p3+p4)/4;
  
  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_position(const mesh_t &mesh, int m, double *t,double *p,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*---------------------------------------------------------------------
  should be used to return computational node position (element barycentre, edge
  mid point, vertex position etc....)
!---------------------------------------------------------------------*/
{
  int status;
  if(mesh.ntriangles!=0) {
    status=fe_position(mesh, mesh.triangles[m], t, p, mode);
    }
  else  {
    status=fe_position(mesh, mesh.quadrangles[m], t, p, mode);
    }
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_position(const mesh_t &mesh, const edge_t & edge, double *t,double *p,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

/*---------------------------------------------------------------------
  should be used to return computational node position (element barycentre, edge
  mid point, vertex position etc....)
!---------------------------------------------------------------------*/
{
  int n1,n2,status;
  double t1,t2,p1,p2;

  n1=edge.extremity[0];
  n2=edge.extremity[1];

  t1=mesh.vertices[n1].lon;
  t2=mesh.vertices[n2].lon;

  if(mode==0) {
    t2=geo_recale(t2,t1,(double) 180.0);
    }
  p1=mesh.vertices[n1].lat;
  p2=mesh.vertices[n2].lat;

  *t=(t1+t2)/2.0;
  *p=(p1+p2)/2.0;

  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_position3D(const mesh_t &mesh, triangle_t triangle,int level, double *t,double *p, double *z,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*---------------------------------------------------------------------
  should be used to return computational node position (element barycentre, edge
  mid point, vertex position etc....)
!---------------------------------------------------------------------*/
{
  int n1,n2,n3,status;
  double t1,t2,t3,p1,p2,p3;
  double z1,z2,z3;

  n1=triangle.vertex[0];
  n2=triangle.vertex[1];
  n3=triangle.vertex[2];

  t1=mesh.vertices[n1].lon;
  t2=mesh.vertices[n2].lon;
  t3=mesh.vertices[n3].lon;

  if(mode==0) {
    t1=geo_recale(t1,t1,(double) 180.0);
    t2=geo_recale(t2,t1,(double) 180.0);
    t3=geo_recale(t3,t1,(double) 180.0);
    }
  p1=mesh.vertices[n1].lat;
  p2=mesh.vertices[n2].lat;
  p3=mesh.vertices[n3].lat;

  z1=mesh.vertices[n1].zlevels[level];
  z2=mesh.vertices[n2].zlevels[level];
  z3=mesh.vertices[n3].zlevels[level];

  *t=(t1+t2+t3)/3;
  *p=(p1+p2+p3)/3;
  *z=(z1+z2+z3)/3;

  return(0);
  }

