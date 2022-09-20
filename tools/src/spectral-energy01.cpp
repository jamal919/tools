
/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/**
\author  Florent Lyard      LEGOS/CNRS, Toulouse, France
  E-mail: florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada
*/

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>

#include "archive.h"

#include "poc-netcdf.hpp"
#include "solverlib.h"

#include "tides.h"
#include "fe.h"
#include "fe-integration.h"
#include "map.h"
#include "mass.h"

#include "filter.h"
#include "statistic.h"
#include "functions.h"
#include "maths.h"
#include "zapper.h"

#include "matrix.h"
#include "constants.h"

using namespace std;

// static double P_g=9.81,rho_water=1027.0,R=6400.*1000.;
static int idx = 0;

extern double **grad_opx,**grad_opy;

extern matrix2x2_t  *DragMatrix;
extern zmatrix2x2_t *FricMatrix;

extern void spectral_friction01(mesh_t, parameter_t, int, int);
extern void spectral_friction02(mesh_t, parameter_t, int, int);
extern void spectral_friction03(spectrum_t spectrum,double duration, parameter_t data, int nndes, int target, double **buffer);
extern void spectral_viscosity(spectrum_t ,double ,mesh_t mesh,parameter_t , int , double **, double **);

extern int fe_LGP1_beta_gradient(mesh_t *mesh, int n, int j, double *C, double *D);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void smooth_vector(mesh_t mesh, double *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j, k, l, m, n, nndes;
  int *cnt;
  double *tmp;
  nndes = mesh.nvtxs;

  cnt = new int[nndes];
  tmp = new double[nndes];

  for(n = 0; n < nndes; n++) {
    tmp[n] = mesh.vertices[n].mw*buffer[n];
    cnt[n] = mesh.vertices[n].mw;
    }

  for(n = 0; n < nndes; n++)
    for(k = 0; k < mesh.vertices[n].nngh; k++) {
      m = mesh.vertices[n].ngh[k];
      cnt[n]+=mesh.vertices[m].mw;
      }

  for(n = 0; n < nndes; n++) {
    for(k = 0; k < mesh.vertices[n].nngh; k++) {
      m = mesh.vertices[n].ngh[k];
      tmp[n] += mesh.vertices[m].mw*buffer[m];
      }
    buffer[n] = tmp[n]/cnt[n];
    }
  delete[] cnt;
  delete[] tmp;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<class T> int fe_LGP1gradient(mesh_t mesh, T *z, int m, T *gx, T *gy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Notes:

  derivation-related coefficient. Definition:

  1/R ds/dlon = dsdx / 2A
  1/R ds/dlat = dsdy / 2A

  dsdx= ds /dlon * 2A/R
  dsdy= ds /dlat * 2A/R

  for LGP1 variable:

  d beta[j] /d lon = triangle.dxdt * beta_x[j] + triangle.dydt * beta_y[j]
                   = +R DQ[j]/2A
  d beta[j] /d lat = triangle.dxdp * beta_x[j] + triangle.dydp * beta_y[j]
                   = -R DP[j]/2A

  dsdx=  sum( DQ[j] * s(j) )
  dsdy= -sum( DP[j] * s(j) )

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int i, n, p, status = 0;
  double x = 0., y = 0., C = 0.;
  T dzdx = 0., dzdy = 0.;
  //double beta_x[3], beta_y[3];//unused

/*   fe_initaffine(mesh,m); */
//  fe_LGP1prime(x, y, beta_x, beta_y);

//   dzdx = 0.;
//   dzdy = 0.;
//   for(i = 0; i < 3; i++) {
//     n = mesh.triangles[m].vertex[i];
// /*----------------------------------------------------------------------
//     Warning : assume latitude in radians*/
//     C += cos(mesh.vertices[n].lat *d2r) / (double) 3.;
//     dzdx += z[n] * beta_x[i];
//     dzdy += z[n] * beta_y[i];
//     }
//   *gx = (dzdx * mesh.triangles[m].dxdt + dzdy * mesh.triangles[m].dydt) / (R * C);
//   *gy = (dzdx * mesh.triangles[m].dxdp + dzdy * mesh.triangles[m].dydp) / R;

  dzdx = 0.;
  dzdy = 0.;
  for(i = 0; i < 3; i++) {
    n = mesh.triangles[m].vertex[i];
    C += cos(mesh.vertices[n].lat *d2r) / (double) 3.;
    dzdx += z[n] * mesh.triangles[m].DQ[i];
    dzdy -= z[n] * mesh.triangles[m].DP[i];
    }
  *gx = dzdx / mesh.triangles[m].Area / (double) 2.0 /C;
  *gy = dzdy / mesh.triangles[m].Area / (double) 2.0;

  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<class T> void fe_gradient02(mesh_t mesh,T *z, T *dzdx, T *dzdy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *------------------------------------------------------------------------
  It returns the gradient of z. Units are identical to input units of z per
  m
------------------------------------------------------------------------*/
{
  int i, j, k, l, m;
  int n;
  T gradx, grady;

  for(n = 0; n < mesh.nvtxs; n++) {
    dzdx[n] = 0;
    dzdy[n] = 0;
    }

  for(n = 0; n < mesh.nvtxs; n++) {
/*----------------------------------------------------------------------
    gradient computation*/
    gradx = grad_opx[n][0] * z[n];
    grady = grad_opy[n][0] * z[n];
    for(k = 0; k < mesh.vertices[n].nngh; k++) {
      m = mesh.vertices[n].ngh[k];
      if(m != -1) {
        gradx += grad_opx[n][k + 1] * z[m];
        grady += grad_opy[n][k + 1] * z[m];
        }
      }
    dzdx[n] = gradx;
    dzdy[n] = grady;
    }
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *energy_divergence(mesh_t mesh, double *buffer_x, double *buffer_y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *-----------------------------------------------------------------------

  Compute spectral energy divergence by flux approach over triangles equation

 -----------------------------------------------------------------------*/
{
  int frame,status;
  int i, i0, i1, i2, j, k, l, m, m1, m2, n,n0, n1, n2;
  triangle_t triangle;

  double x,y,U,V;
  double lateral[3];
  double pseudo,area;
  double mean;
  double c[3],p[3],q[3],r[3],s[3];
  double UU[3],VV[3],dlon,dlat;

  double *rhs,*div,w;

/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

  19/08/2009 :

    1) compute divergence integral over the mesh triangles by integrating
       fluxes on triangle's sides

    2) divide by total triangle area to obtain a mean value

    3) retrieve LGP1 discretisation by solving inverse problem

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */


/*-----------------------------------------------------------------------
  allocate arrays*/
  rhs=new double[mesh.nvtxs];
  div=new double[mesh.ntriangles];

  for(m = 0; m < mesh.ntriangles; m++) {
    triangle = mesh.triangles[m];
    area   = triangle.TrueArea;
    pseudo = triangle.Area;
/*-----------------------------------------------------------------------
    compute lateral flux, exact algorithm*/
    n0=triangle.vertex[0];
    n1=triangle.vertex[1];
    n2=triangle.vertex[2];
    UU[0]=buffer_x[n0];
    UU[1]=buffer_x[n1];
    UU[2]=buffer_x[n2];
    VV[0]=buffer_y[n0];
    VV[1]=buffer_y[n1];
    VV[2]=buffer_y[n2];
    div[m]=0;
    for(j = 0; j < 3; j++) {
      n=triangle.vertex[j];
      //c[j]=mesh.vertices[n].c;//value not used
      }
    for(j = 0; j < 3; j++) {
      n=triangle.edges[j];
/*-----------------------------------------------------------------------
      edge 0 : vertex 1 and 2 etc...*/
      i0=(j+1) %3;
      i1=(j+2) %3;
      n0 = mesh.triangles[m].vertex[i0];
      n1 = mesh.triangles[m].vertex[i1];
      p[0]= UU[i0];
      p[1]= UU[i1];
      q[0]= VV[i0];
      q[1]= VV[i1];
//       r[0]= state2D[present].H[n0];
//       r[1]= state2D[present].H[n1];
      s[0]= cos(mesh.vertices[n0].lat*d2r);
      s[1]= cos(mesh.vertices[n1].lat*d2r);
      dlon=(mesh.vertices[n1].lon-mesh.vertices[n0].lon)*d2r;
      dlat=(mesh.vertices[n1].lat-mesh.vertices[n0].lat)*d2r;
//       sum  =R*dlat*fe_integraleLGP1xLGP1_1D(p,r);
//       sum -=R*dlon*fe_integraleLGP1xLGP1xLGP1_1D(q,r,s);
      div[m] +=R*dlat*fe_integraleLGP1_1D(p);
      div[m] -=R*dlon*fe_integraleLGP1xLGP1_1D(q,s);
      }
    div[m]/=area;
    }

  for(n = 0; n < mesh.nvtxs; n++) {
    rhs[n]  = 0.0;
    }

/* *-----------------------------------------------------------------------
  Nodale quadrature */
  for(m = 0; m < mesh.ntriangles; m++) {
    pseudo = mesh.triangles[m].Area;
    for(j = 0; j < 3; j++) {
      n=mesh.triangles[m].vertex[j];
//      rhs[n]  += 3.*div[m]*mesh.vertices[n].mw/pseudo;
      rhs[n]  += div[m]*pseudo/3./mesh.vertices[n].mw;
      }
    }
 
/* *-----------------------------------------------------------------------
  Integrale scalar product */
//  for(m = 0; m < mesh.ntriangles; m++) {
//    pseudo = mesh.triangles[m].Area;
//     for(j = 0; j < 3; j++) {
//       n=triangle.vertex[j];
//       rhs[n]  += 2.0*pseudo*div[m]*fe_integraleLGP1xLGP1_2D(c,j);
//       }
//    }
//  status = MassMatrix_solve(mesh,&(mesh.matrix.t),mesh.matrix.s,rhs);

/*-----------------------------------------------------------------------
  deallocate arrays*/
  delete[] div;
  return(rhs);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  Compute astronomical potential gradient constants (complex)
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void tidal_potential_harmonics_gradient(tidal_wave w, mesh_t mesh, parameter_t P1_data,
                                  dcomplex *P[3])
{
  int n, k;
  double a, V, dV, V0, C, S, x, y;
  double k2 = 0.3, h2 = 0.6; // slope and gravity modification
  double omega;

/* *----------------------------------------------------------------------------
  warning: convert centimeters to meters*/
  a = w.Ap * (1 + k2 - h2);
  switch (w.nT) {
    case (0):

/*######################### Long period tide #########################

   potential/g = A*(1/2 -3/2 sin^2(Y)*cos(w*t+V0)

   dP/dx=  0
   dP/dy= -3*A*cos(Y)sin(Y)*cos(w*t+V0)

----------------------------------------------------------------------*/

      for(n = 0; n < mesh.nvtxs; n++) {
        dV = 0.0;
        C = P1_data.C[n];
        S = P1_data.S[n];
        P[0][n] = polar( a * (0.5 - 1.5 * S * S), dV);
        P[1][n] = 0;
        P[2][n] = polar( -3 * a * C * S, dV);
      }
      break;

    case (1):

/*########################### Diurnal tide ###########################

   potential/g = A*sin(2Y)*cos(w*t+V0+X)
   
   dP/dx=  -2*A*cos(Y)*sin(Y)*sin(w*t+V0+X)
   dP/dy=   2*A*cos(2Y)*cos(w*t+V0+X) = 2*A*[cos^2(Y)-sin^2(Y)]*cos(w*t+V0+X)

----------------------------------------------------------------------*/

      for(n = 0; n < mesh.nvtxs; n++) {
        dV = P1_data.lon[n];
        C = P1_data.C[n];
        S = P1_data.S[n];
        P[0][n] = polar( 2 * a * S * C, dV);
        P[1][n] = polar( -2 * a * C * S, dV - M_PI / 2);
        P[2][n] = polar( 2 * a * (C * C - S * S), dV);
        }
      break;

    case (2):

/*######################### Semi-diurnal tide #########################

   potential/g = A*cos^2(Y)*cos(w*t+V0+2*X)

   dP/dx=  -2*A*cos^2(Y)*sin(w*t+V0+2*X)
   dP/dy=  -2*A*cos(Y)*sin(Y)*cos(w*t+V0+2*X)

----------------------------------------------------------------------*/

      for(n = 0; n < mesh.nvtxs; n++) {
        dV = 2 * P1_data.lon[n];
        C = P1_data.C[n];
        S = P1_data.S[n];
        P[0][n] = polar( a * C * C, dV);
        P[1][n] = polar( -2 * a * C * C, dV - M_PI / 2);
        P[2][n] = polar( -2 * a * C * S, dV);
        }
      break;

/*####################### non-astronomical tide #######################
   potential/g = 0

----------------------------------------------------------------------*/
    default:
      for(n = 0; n < mesh.nvtxs; n++) {
        P[0][n] = 0;
        P[1][n] = 0;
        P[2][n] = 0;
        }
      break;

    }

/*######################### negative amplitudes #########################
   just invert phase
   not usefull for complex
----------------------------------------------------------------------*/
//   for(n = 0; n < mesh.nvtxs; n++) {
//     for(k = 0; k < 3; k++) {
//       if(Pa[k][n] < 0) {
//         Pa[k][n] = -Pa[k][n];
//         Pg[k][n] += M_PI;
//         }
//       }
//     }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void spectral_flux02(mesh_t mesh,parameter_t data, int w)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
  Compute the spectral energy fluxes at the open boundaries*/
{
  double Nx, Ny, L, gPSI, rhoPSI;
  double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0;
  double h, F, Fx, Fy;
  double a1, a2, G1, G2, a3, G3;
  double p[2],q[2],r[2],s[2];
  dcomplex u,u1,u2,v,v1,v2,z,z1,z2;
  int i,k,n,n1,n2,node;

  gPSI = P_g;
  rhoPSI = rho_water;

  for(n = 0; n < mesh.nedges; n++) {
    if(mesh.edges[n].code != MESH_INTERIOR_EDGE) {
/*-----------------------------------------------------------------------
      inward normal vector*/
      Nx = -mesh.edges[n].Ty;
      Ny =  mesh.edges[n].Tx;
      L =   mesh.edges[n].L;
      for(k=0;k<2;k++) {
        node=mesh.edges[n].extremity[k];
        p[k]=data.h[node];
        a1 = data.Htide[node].a[w];
        G1 = data.Htide[node].G[w];
        a2 = data.Utide[node].a[w];
        G2 = data.Utide[node].G[w];
        a3 = data.Vtide[node].a[w];
        G3 = data.Vtide[node].G[w];
        z1=a1*cos(G1);
        u1=a2*cos(G2);
        v1=a3*cos(G3);
        z2=a1*sin(G1);
        u2=a2*sin(G2);
        v2=a3*sin(G3);
        q[k]= 0.5 * a1 * a2 * cos(G2 - G1);
        r[k]= 0.5 * a1 * a3 * cos(G3 - G1);
        }
/*----------------------------------------------------------------------
      integration along the triangle side*/
      Fx = fe_integraleLGP1xLGP1_1D (p,q);
      Fy = fe_integraleLGP1xLGP1_1D (p,r);
      F  = gPSI * rhoPSI * L * (Fx * Nx + Fy * Ny);
/*----------------------------------------------------------------------
      Potential energy flux (energy budget form)*/
      sum2 += F;
      }
    }

  printf("Potential energy flux : %12.6f gW\n",sum2/1.e+09);

  for(n = 0; n < mesh.nedges; n++) {//for all edges
    if(mesh.edges[n].code == 5) {
/*-----------------------------------------------------------------------
      inward normal vector*/
      Nx = -mesh.edges[n].Ty;
      Ny =  mesh.edges[n].Tx;
      L =   mesh.edges[n].L;
      for(k=0;k<2;k++) {
        node=mesh.edges[n].extremity[k];
        p[k]=data.h[node];
        a1 = data.Htide[node].a[w];
        G1 = data.Htide[node].G[w];
        a2 = data.Utide[node].a[w];
        G2 = data.Utide[node].G[w];
        a3 = data.Vtide[node].a[w];
        G3 = data.Vtide[node].G[w];
        z1=a1*cos(G1);
        u1=a2*cos(G2);
        v1=a3*cos(G3);
        z2=a1*sin(G1);
        u2=a2*sin(G2);
        v2=a3*sin(G3);
        q[k]= a1;
        r[k]= 0.5 * a2 * cos(G2 - G1);
        s[k]= 0.5 * a3 * cos(G3 - G1);
        }
/*----------------------------------------------------------------------
      integration along the triangle side*/
      Fx = fe_integraleLGP1xLGP1xLGP1_1D (p,q,r);
      Fy = fe_integraleLGP1xLGP1xLGP1_1D (p,q,s);
      F  = gPSI * rhoPSI * L * (Fx * Nx + Fy * Ny);
/*----------------------------------------------------------------------
      Potential energy flux (energy budget form)*/
      if(F>0) {
        sum3 += F;
        }
      else {
        sum4 += F;
        }
      }
    }

  printf("Potential energy flux : %12.6f %12.6f %12.6f gW\n",sum3/1.e+09,sum4/1.e+09,(sum3+sum4)/1.e+09);

  for(n = 0; n < mesh.nedges; n++) {//for all edges
    if(mesh.edges[n].code == 5) {
      n1 = mesh.edges[n].extremity[0];
      n2 = mesh.edges[n].extremity[1];

/*-----------------------------------------------------------------------
      inward normal vector*/
      Nx = -mesh.edges[n].Ty;
      Ny =  mesh.edges[n].Tx;
      L =   mesh.edges[n].L;

      a1 = data.Htide[n1].a[w];
      G1 = data.Htide[n1].G[w];
      a2 = data.Utide[n1].a[w];
      G2 = data.Utide[n1].G[w];
      a3 = data.Vtide[n1].a[w];
      G3 = data.Vtide[n1].G[w];

      z1=a1*dcomplex(cos(G1),-sin(G1));
      u1=a2*dcomplex(cos(G2),-sin(G2));
      v1=a3*dcomplex(cos(G3),-sin(G3));

      a1 = data.Htide[n2].a[w];
      G1 = data.Htide[n2].G[w];
      a2 = data.Utide[n2].a[w];
      G2 = data.Utide[n2].G[w];
      a3 = data.Vtide[n2].a[w];
      G3 = data.Vtide[n2].G[w];

      z2=a1*dcomplex(cos(G1),-sin(G1));
      u2=a2*dcomplex(cos(G2),-sin(G2));
      v2=a3*dcomplex(cos(G3),-sin(G3));

      z=0.5*(z1+z2);
      u=0.5*(u1+u2);
      v=0.5*(v1+v2);

      h = 0.5*(data.h[n1]+data.h[n2]);

      a1 =  abs(z);
      G1 = -arg(z);
      a2 =  abs(u);
      G2 = -arg(u);
      a3 =  abs(v);
      G3 = -arg(v);
//  double fe_integraleLGP1xLGP1xLGP1_1D(double *p, double *q, double *r)

      Fx = 0.5 * h * gPSI * rhoPSI *a1 * a2 * cos(G2 - G1);
      Fy = 0.5 * h * gPSI * rhoPSI *a1 * a3 * cos(G3 - G1);

      F = L * (Fx * Nx + Fy * Ny);

/*----------------------------------------------------------------------
      Potential energy flux (energy budget form)*/
      sum1 += F;
      }
    }
  printf("Potential energy flux : %12.6f gW\n",sum1/1.e+09);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void spectral_flux03(mesh_t mesh,parameter_t data, int w, int ww)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
  Compute the spectral energy fluxes at the open boundaries*/
{
  double Nx, Ny, L, gPSI, rhoPSI;
  double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0;
  double h, F, Fx, Fy;
  double a1, a2, G1, G2, a3, G3, a[2],G[2];
  double az1,Gz1,au1,Gu1,av1,Gv1,az2,Gz2,au2,Gu2,av2,Gv2;
  dcomplex u[2],u1,u2,u3,u4,v[2],v1,v2,v3,v4,z[2],z1,z2,z3,z4;
  int i,k,n,n1,n2;

  gPSI = P_g;
  rhoPSI = rho_water;

  for(n = 0; n < mesh.nedges; n++) {//for all edges
    if(mesh.edges[n].code == 5) {
      n1 = mesh.edges[n].extremity[0];
      n2 = mesh.edges[n].extremity[1];

/*-----------------------------------------------------------------------
     outward normal vector*/
      Nx = -mesh.edges[n].Ty;
      Ny =  mesh.edges[n].Tx;
      L =   mesh.edges[n].L;

      a1 = data.Htide[n1].a[w];
      G1 = data.Htide[n1].G[w];
      a2 = data.Utide[n1].a[w];
      G2 = data.Utide[n1].G[w];
      a3 = data.Vtide[n1].a[w];
      G3 = data.Vtide[n1].G[w];

      z1=a1*dcomplex(cos(G1),-sin(G1));
      u1=a2*dcomplex(cos(G2),-sin(G2));
      v1=a3*dcomplex(cos(G3),-sin(G3));

      a1 = data.Htide[n2].a[w];
      G1 = data.Htide[n2].G[w];
      a2 = data.Utide[n2].a[w];
      G2 = data.Utide[n2].G[w];
      a3 = data.Vtide[n2].a[w];
      G3 = data.Vtide[n2].G[w];

      z2=a1*dcomplex(cos(G1),-sin(G1));
      u2=a2*dcomplex(cos(G2),-sin(G2));
      v2=a3*dcomplex(cos(G3),-sin(G3));

      z[0]=0.5*(z1+z2);
      u[0]=0.5*(u1+u2);
      v[0]=0.5*(v1+v2);

      a1 = data.Htide[n1].a[ww];
      G1 = data.Htide[n1].G[ww];
      a2 = data.Utide[n1].a[ww];
      G2 = data.Utide[n1].G[ww];
      a3 = data.Vtide[n1].a[ww];
      G3 = data.Vtide[n1].G[ww];
 
      z1=a1*dcomplex(cos(G1),-sin(G1));
      u1=a2*dcomplex(cos(G2),-sin(G2));
      v1=a3*dcomplex(cos(G3),-sin(G3));

      a1 = data.Htide[n2].a[ww];
      G1 = data.Htide[n2].G[ww];
      a2 = data.Utide[n2].a[ww];
      G2 = data.Utide[n2].G[ww];
      a3 = data.Vtide[n2].a[ww];
      G3 = data.Vtide[n2].G[ww];

      z2=a1*dcomplex(cos(G1),-sin(G1));
      u2=a2*dcomplex(cos(G2),-sin(G2));
      v2=a3*dcomplex(cos(G3),-sin(G3));

      z[1]=0.5*(z1+z2);
      u[1]=0.5*(u1+u2);
      v[1]=0.5*(v1+v2);

      h = 0.5*(data.h[n1]+data.h[n2]);

      az1 =  abs(z[0]);
      Gz1 = -arg(z[0]);
      au1 =  abs(u[0]);
      Gu1 = -arg(u[0]);
      av1 =  abs(v[0]);
      Gv1 = -arg(v[0]);

      az2 =  abs(z[1]);
      Gz2 = -arg(z[1]);
      au2 =  abs(u[1]);
      Gu2 = -arg(u[1]);
      av2 =  abs(v[1]);
      Gv2 = -arg(v[1]);
/*----------------------------------------------------------------------
     Non linear energy flux */
      Fx = 0.25 * gPSI * rhoPSI *az1 * az2 * au1 * cos(Gz1 - Gz2 + Gu1);
      Fy = 0.25 * gPSI * rhoPSI *az1 * az2 * av1 * cos(Gz1 - Gz2 + Gv1);
      F = L * (Fx * Nx + Fy * Ny);
/*----------------------------------------------------------------------
     Potential energy flux (energy budget form)*/
      sum1 += F;
/*----------------------------------------------------------------------
     Non linear kinetic energy flux */
      Fx = 0.25 * h * rhoPSI *au1 * au2 * au1 * cos(Gu1 - Gu2 + Gu1);
      Fy = 0.25 * h * rhoPSI *av1 * av2 * av1 * cos(Gv1 - Gv2 + Gv1);

      F = L * (Fx * Nx + Fy * Ny);
/*----------------------------------------------------------------------
      Potential energy flux (energy budget form)*/
      sum2 += F;
    }
  }
  printf("Non-linear potential energy flux : %8.1f mW\n",sum1/1.e+06);
  printf("Non-linear kinetic energy flux   : %8.1f mW\n",sum2/1.e+06);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double flux_integrale(mesh_t mesh, double **buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#warning flux_integrale is a duplicate
/*----------------------------------------------------------------------
  Compute the spectral energy fluxes at the open boundaries*/
{
  double Nx, Ny, L, gPSI, rhoPSI;
  double sum1 = 0;
  double h, F, Fx, Fy;
  double zx,zy;
  int i,k,n,n1,n2;

  gPSI = P_g;
  rhoPSI = rho_water;

  for(n = 0; n < mesh.nedges; n++) {//for all edges
    if(mesh.edges[n].code == 5) {
      n1 = mesh.edges[n].extremity[0];
      n2 = mesh.edges[n].extremity[1];

/*-----------------------------------------------------------------------
      outward normal vector*/
      Nx = -mesh.edges[n].Ty;
      Ny =  mesh.edges[n].Tx;
      L =   mesh.edges[n].L;

      zx=0.5*(buffer[n1][0]+buffer[n2][0]);
      zy=0.5*(buffer[n1][1]+buffer[n2][1]);

      F = L * rhoPSI * gPSI *(zx * Nx + zy * Ny);

/*----------------------------------------------------------------------
      Potential energy flux (energy budget form)*/
      sum1 += F;
    }
  }
  return(sum1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void spectral_flux_budget(mesh_t mesh,spectrum_t spectrum,double duration,parameter_t data)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double c,R1,R2,R3;
  double a1,G1,a2,G2,a3,G3;
  double *residuals,*taux,*tauy,tauxm,tauym;
  int k,m,n,nframe,w,status;
  int maxstep=1;
  double *u,*v,*H,*ax,*Gx,*ay,*Gy;
  harmonic_t harmonic;
  date_t start;
  double dt,*time;
  tidal_wave wave;
  double f,Vu,V, V0, omega, norm;
  double cs[3],sn[3];
  int nodal=0;
  double ***buffer,sum;

  start.year=1990;
  start.month=1;
  start.day=1;
  start.second=0;

  dt=3600.0;
  nframe=(int) NINT(duration/dt);
 /*-----------------------------------------------------------------------*/
  astro_angles_t astro_angles;
  init_argument(&astro_angles,start);
/*
  duration=3600.*24*15;
  dt=2.*M_PI/(spectrum.waves[0].omega*d2r/3600.)/24.0;
  nframe=(int) NINT(duration/dt);
*/
  buffer=new double **[spectrum.n];
  for(k=0;k<spectrum.n;k++) {
    buffer[k]=new double * [mesh.nvtxs];
    for(n=0;n<mesh.nvtxs;n++)  buffer[k][n]=new double [2];
    }

  time = new double[nframe];
  u    = new double[nframe];
  v    = new double[nframe];
  H    = new double[nframe];
  taux = new double[nframe];
  tauy = new double[nframe];
  residuals = new double[nframe];

  ax= new double [spectrum.n];
  Gx= new double [spectrum.n];
  ay= new double [spectrum.n];
  Gy= new double [spectrum.n];

  for(k=0;k<nframe;k++) time[k]=(double)k*dt;

  status=harmonic_init01(nframe, start, time, spectrum, &harmonic,0,astro_angles);

  for(n=0;n<mesh.nvtxs;n++) {
    for(m=0;m<nframe;m++) {
      u[m]=0;
      v[m]=0;
      H[m]=data.h[n];
      }
    for(w = 0; w < spectrum.n; w++) {
      a1=data.Utide[n].a[w];
      G1=data.Utide[n].G[w];
      a2=data.Vtide[n].a[w];
      G2=data.Vtide[n].G[w];
      a3=data.Htide[n].a[w];
      G3=data.Htide[n].G[w];
      cs[0]=a1*cos(G1);
      cs[1]=a2*cos(G2);
      cs[2]=a3*cos(G3);
      sn[0]=a1*sin(G1);
      sn[1]=a2*sin(G2);
      sn[2]=a3*sin(G3);
      for(m=0; m<nframe; m++) {
        u[m] += cs[0] * harmonic.cs[w][m] + sn[0] * harmonic.sn[w][m];
        v[m] += cs[1] * harmonic.cs[w][m] + sn[1] * harmonic.sn[w][m];
        H[m] += cs[2] * harmonic.cs[w][m] + sn[2] * harmonic.sn[w][m];
        }
      }

    for(m=0;m<nframe;m++) {
      taux[m]=u[m]*H[m];
      tauy[m]=v[m]*H[m];
      }
    status=harmonic_compute(taux,residuals, nframe, harmonic, spectrum, maxstep,ax,Gx);
    status=harmonic_compute(tauy,residuals, nframe, harmonic, spectrum, maxstep,ay,Gy);
    for(w = 0; w < spectrum.n; w++) {
      a1=data.Htide[n].a[w];
      G1=data.Htide[n].G[w];
      buffer[w][n][0]=0.5*ax[w]*a1*cos(Gx[w]-G1);
      buffer[w][n][1]=0.5*ay[w]*a1*cos(Gy[w]-G1);
      }
    }
  for(w = 0; w < spectrum.n; w++) {
    sum= flux_integrale(mesh, &(buffer[w][0]));
    printf("%s total potential energy flux : %12.4f gW\n",spectrum.waves[w].name,sum/1.e+09);
    }

  for(n=0;n<mesh.nvtxs;n++) {
    for(m=0;m<nframe;m++) {
      u[m]=0;
      v[m]=0;
      H[m]=0;
      }
    for(w = 0; w < spectrum.n; w++) {
      a1=data.Utide[n].a[w];
      G1=data.Utide[n].G[w];
      a2=data.Vtide[n].a[w];
      G2=data.Vtide[n].G[w];
      a3=data.Htide[n].a[w];
      G3=data.Htide[n].G[w];
      cs[0]=a1*cos(G1);
      cs[1]=a2*cos(G2);
      cs[2]=a3*cos(G3);
      sn[0]=a1*sin(G1);
      sn[1]=a2*sin(G2);
      sn[2]=a3*sin(G3);
      for(m=0; m<nframe; m++) {
        u[m] += cs[0] * harmonic.cs[w][m] + sn[0] * harmonic.sn[w][m];
        v[m] += cs[1] * harmonic.cs[w][m] + sn[1] * harmonic.sn[w][m];
        H[m] += cs[2] * harmonic.cs[w][m] + sn[2] * harmonic.sn[w][m];
        }
      }

    for(m=0;m<nframe;m++) {
      taux[m]=u[m]*H[m];
      tauy[m]=v[m]*H[m];
      }
    status=harmonic_compute(taux,residuals, nframe, harmonic, spectrum, maxstep,ax,Gx);
    status=harmonic_compute(tauy,residuals, nframe, harmonic, spectrum, maxstep,ay,Gy);
    for(w = 0; w < spectrum.n; w++) {
      a1=data.Htide[n].a[w];
      G1=data.Htide[n].G[w];
      buffer[w][n][0]=0.5*ax[w]*a1*cos(Gx[w]-G1);
      buffer[w][n][1]=0.5*ay[w]*a1*cos(Gy[w]-G1);
      }
    }
  for(w = 0; w < spectrum.n; w++) {
    sum= flux_integrale(mesh, &(buffer[w][0]));
    printf("%s partial potential energy flux : %12.4f gW\n",spectrum.waves[w].name,sum/1.e+09);
    }

  delete[] time;

  delete[] u;
  delete[] v;
  delete[] H;

  delete[] taux;
  delete[] tauy;
  delete[] residuals;

  delete[] ax;
  delete[] ay;
  delete[] Gx;
  delete[] Gy;
}



/*----------------------------------------------------------------------------*/
/// Gives the spectral energy budget.
/**
\date last reviewed 23 Jun 2011

\todo clean up

\param mesh grid
\param AnalysisList list of tidal frequencies to do the analysis for
\param *P1_action forcing on each node
\param *output_path path to output directory
\param *path path to directory of netcdf input
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void spectral_energy01(mesh_t & mesh,spectrum_t & AnalysisList, parameter_t & data2D, action_t *P1_action, const char *output_path, const char *path=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, status;					//Indexes of corners or wave, sum or cornered triangles, anything. Status
  int l, m, n, w;					//Indexes of edges, triangles, nodes and waves.
  double A, C, D, gPSI, rhoPSI, hpg_2;			//, U;Area, cos(latitude) or latitude gradient, longitude gradient, g, volumetric mass, buffer
  double uu, a2, a3;					//total speed and its amplitudes
  double zr;						//, zi;real and imaginary parts of pressure
  double row[10][3];
  double beta[3];
  double tmp[10], sum[10];				//, budget;buffer, integration
  double h, th,tLa,tLo;					//, H, dH, dHdT;node depth, corners' triangle mean: depth, lat, lon. 3 unused
  double ch[3],cc[3];					// corners depth and cosine
  double dum,tmpx, tmpy;				//, c1, c2;buffers
  dcomplex zx, zy;					//field gradient components
  double *Pa, *PG;					//pressure phase and amplitude
  int nMin;						//index of minimum
  dcomplex *P, *Px, *Py;				//complex pressure and its gradients
  dcomplex u, v, taux,tauy;				//complex speed components and its related friction
  double *uMax;						//maximum speed
  double G2, G3;					//speed components' phases
  double *buf, mask=1.e+10;				//buffer
  dcomplex *L, *Lx, *Ly;				//complex secondary potential and its LSA gradient
  double *RoW[10];
  double *meanKE,*meanWD;
  int kZ0,kM2,kM4;					//node index, wave frequency indicator
  tidal_wave wave;
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int gauss_n = 7;
  char filename[1024], ncfilename[1024];
  int haveP0;
  FILE *out;

  dcomplex *Ph[3];					//complex tidal potential harmonics
  double *HDFrow[100],*BFRrow[100],*ADVrow[100],alpha=0.1,rate,*gratio,*gvalue,*z0;
  double r;						//, c, eps;//Cd, 2 unused
  double tau,duration;					//wave period/frequency stuff
  dcomplex jw, jw_3;
  double *local[10];					//spatial integration
  
  int ncid;
  
  char *comment[2];//first 2 lines of data file
  comment[0]=new char[256];
  comment[1]=new char[256];

/*----------------------------------------------------------------------
  Compute the energy budget over the domain*/
  duration=15*24*3600.;
  for(i = 0; i < AnalysisList.n; i++) {
    for(j = i + 1; j < AnalysisList.n; j++) {
      tau = AnalysisList.waves[j].omega - AnalysisList.waves[i].omega;
      if(tau == 0.) {
        printf("wave: %10s %10s, separation: %9.3f days \n", AnalysisList.waves[i].name, AnalysisList.waves[j].name, tau);
        }
      else {
        tau = (1. / tau) * 360. *3600.;
        duration=MAX(duration,tau);
        }
      }
    }

  //This variable replacement is somewhat usefull for code parsing.
  gPSI = P_g;
  rhoPSI = rho_water;

  for(k = 0; k < 3; k++) {
    Ph[k] = new dcomplex[mesh.nvtxs];
    }

  for(k = 0; k < 10; k++) {
    RoW[k]   = new double[mesh.nvtxs];
    local[k] = new double[mesh.ntriangles];
    }

  buf = new double[mesh.nedges];//edges are the most numerous

  P = new dcomplex[mesh.nvtxs];
  Px = new dcomplex[mesh.nvtxs];
  Py = new dcomplex[mesh.nvtxs];
  Lx = new dcomplex[mesh.nvtxs];
  Ly = new dcomplex[mesh.nvtxs];
  L = new dcomplex[mesh.nvtxs];

  Pa = new double[mesh.nvtxs];
  PG = new double[mesh.nvtxs];

  uMax = new double[mesh.nvtxs];

  gratio = new double[mesh.nvtxs];
  gvalue = new double[mesh.nvtxs];

  z0     = new double[mesh.nvtxs];

  meanKE = new double[mesh.nvtxs];
  meanWD = new double[mesh.nvtxs];

  for(n = 0; n < mesh .nvtxs; n++) {// for all nodes
    gratio[n] = 0.0;
    gvalue[n] = 0.0;
    z0[n] = 2.e-03;
    meanKE[n] = 0.0;
    meanWD[n] = 0.0;
    }

  FricMatrix = new zmatrix2x2_t[mesh.nvtxs];//initialised in spectral_friction01 and spectral_friction02
  gauss_init(gauss_n, gauss_x, gauss_y, gauss_w);

  kZ0=kM2=kM4=-1;
  for(w = 0; w < AnalysisList.n; w++) {
    wave = AnalysisList.waves[w];
    if(strcmp("Z0",wave.name)==0) kZ0=w;
    if(strcmp("M2",wave.name)==0) kM2=w;
    if(strcmp("M4",wave.name)==0) kM4=w;
    ADVrow[w] = new double[mesh.nvtxs];
    HDFrow[w] = new double[mesh.nvtxs];
    BFRrow[w] = new double[mesh.nvtxs];
    }

/*----------------------------------------------------------------------
  maximum speeds */
  for(w = 0; w < AnalysisList.n; w++) {
    wave = AnalysisList.waves[w];
    if(kM2!=w)//if not M2
      continue;//skip
    for(n = 0; n < mesh.nvtxs; n++) {// for all nodes
      a2 = data2D.Utide[n].a[w];
      G2 = data2D.Utide[n].G[w];
      a3 = data2D.Vtide[n].a[w];
      G3 = data2D.Vtide[n].G[w];
      a2*=a2*.5;
      a3*=a3*.5;
      G2-=G3;
      uMax[n]=sqrt(a2+a3+sqrt(a2*a2+a3*a3+2*a2*a3*cos(2*G2)));
      }
    }
/*----------------------------------------------------------------------
  compute non-linear contributions*/
/* *----------------------------------------------------------------------
  disconnected*/
//  spectral_flux(mesh,AnalysisList,duration,state2D,data2D);

#define NONLINEAR
#ifdef NONLINEAR
  spectral_viscosity (AnalysisList,duration,mesh,data2D,w,HDFrow,ADVrow);
  spectral_friction03(AnalysisList,duration,data2D,mesh.nvtxs, w, BFRrow);
#endif

//  spectral_flux03(mesh, data2D, kM2,kM4);

  for(w = 0; w < AnalysisList.n; w++) {
    wave = AnalysisList.waves[w];
/*----------------------------------------------------------------------
    compute linearized friction coeffcients*/
    if(w==kM2){
/*----------------------------------------------------------------------
      dominant wave*/
      spectral_friction01(mesh, data2D,mesh.nvtxs, kM2);
      }
    else {
/*----------------------------------------------------------------------
      secondary wave*/
      if(kM2!=-1) {
        spectral_friction02(mesh, data2D,mesh.nvtxs, kM2);
        }
      else {
        spectral_friction01(mesh, data2D,mesh.nvtxs, w);
        }
      }
    spectral_flux02(mesh, data2D, w);
/* *----------------------------------------------------------------------------
    T-UGO specific*/
//    spectral_flux01( mesh, w, gZspec[0], gNobc_Z[0]);

/* *----------------------------------------------------------------------
    astronomical forcing*/
    tidal_potential_harmonics_gradient(wave, mesh, data2D, Ph);

/* *----------------------------------------------------------------------
    pressure gradient*/
    for(n = 0; n < mesh.nvtxs; n++) {// for all nodes
      P[n]=data2D.Htide[n].z[w]; // pressure <- tidal elevation
      PG[n] = data2D.Htide[n].G[w]; // phase
      Pa[n] = data2D.Htide[n].a[w]; // amplitude
      }
    fe_gradient02(mesh,P, Px, Py);

/* *----------------------------------------------------------------------
    LSA gradient*/
    for(n = 0; n < mesh.nvtxs; n++) {// for all nodes
      L[n] = P1_action[n].LSA[w]; // Loading and Self Attraction (LSA)
      }
    fe_gradient02(mesh,L, Lx, Ly);

    for(j = 0; j < 10; j++)
      sum[j] = 0.0;

    for(m = 0; m < mesh.ntriangles; m++) {
/*----------------------------------------------------------------------
      A (triangle area) in m**2*/
      A = mesh.triangles[m].Area;
      
      for(j = 0; j < 10; j++)
        tmp[j] = 0;
      for(i = 0; i < 3; i++) {
        n = mesh.triangles[m].vertex[i];
/*----------------------------------------------------------------------
        depth in m*/
        h = data2D.h[n]; // mean depth
        C = data2D.C[n]; // cosine of the latitude
        hpg_2=0.5 * h * rhoPSI * gPSI; //for code brevity
/*----------------------------------------------------------------------
        amplitude in m/s*/
        u = data2D.Utide[n].z[w];
        a2 = data2D.Utide[n].a[w];
        G2 = data2D.Utide[n].G[w];
        v = data2D.Vtide[n].z[w];
        a3 = data2D.Vtide[n].a[w];
        G3 = data2D.Vtide[n].G[w];

/* *----------------------------------------------------------------------
        astronomical potential rate of work*/
        zx = Ph[1][n];
        zy = Ph[2][n];
        //NOTE: operator%(complex,complex) : scalar product of 2 vectors
        row[0][i] = hpg_2 * ( u % zx / C + v % zy ) * C / R;
        RoW[0][n] = hpg_2 * ( u % zx / C + v % zy ) / R;
/* *----------------------------------------------------------------------
        secondary potential rate of work*/
        zx = Lx[n];
        zy = Ly[n];
        row[1][i] = hpg_2 * ( u % zx + v % zy ) * C;
        RoW[1][n] = hpg_2 * ( u % zx + v % zy );
/* *----------------------------------------------------------------------
        wave drag rate of work*/
        tmpx = DragMatrix[n].c[0][0] * a2 * a2 +
               DragMatrix[n].c[0][1] * a2 * a3 * cos(G3 - G2);
        tmpy = DragMatrix[n].c[1][0] * a2 * a3 * cos(G3 - G2) +
               DragMatrix[n].c[1][1] * a3 * a3;
        row[2][i] = 0.5 * h * rhoPSI * (tmpx + tmpy) * C;
        RoW[2][n] = 0.5 * h * rhoPSI * (tmpx + tmpy);
/* *----------------------------------------------------------------------
        pressure gradient rate of work*/
        zx = Px[n];
        zy = Py[n];
        row[3][i] = -hpg_2 * ( u % zx + v % zy ) * C;
        RoW[3][n] = -hpg_2 * ( u % zx + v % zy );
        
        status = fe_LGP1gradient(mesh, P, m, &zx, &zy);
        row[5][i] = -hpg_2 * ( u % zx + v % zy ) * C;
        RoW[5][n] = -hpg_2 * ( u % zx + v % zy );
/* *----------------------------------------------------------------------
        bottom friction rate of work, quasi-linearised algorithm*/
        taux = FricMatrix[n].c[0][0]*u + FricMatrix[n].c[1][0]*v;
        tauy = FricMatrix[n].c[0][1]*u + FricMatrix[n].c[1][1]*v;
        tmpx = taux % u * 0.5;
        tmpy = tauy % v * 0.5;
        row[6][i] = -rhoPSI * data2D.Cd[n] * (tmpx + tmpy) * C;
        RoW[6][n] = -rhoPSI * data2D.Cd[n] * (tmpx + tmpy);
/*----------------------------------------------------------------------
        bottom friction rate of work*/
/*
        if(h > 1000.) {
          eps = 200.0 / h;
          c = 1. / (eps + alpha * (1. - eps));
          r = c * c * (Cm + alpha * alpha * data2D[n].vfc);
          }
        else {
          r=data2D[n].vfc;
          }
*/
        r=data2D.Cd[n];
/*test change*/
//        row[7][i] = -rhoPSI* h * r * BFRrow[w][n] * C;
//        RoW[7][n] = -rhoPSI* h * r * BFRrow[w][n];
/* *----------------------------------------------------------------------
        bottom friction rate of work, spectral analysis algorithm*/
        row[7][i] = -rhoPSI * r * BFRrow[w][n] * C;
        RoW[7][n] = -rhoPSI * r * BFRrow[w][n];
/* *----------------------------------------------------------------------
        temporary patch !!! */
//         row[7][i] = row[6][i];
//         RoW[7][n] = RoW[6][n];

/* *----------------------------------------------------------------------
        horizontal momentum diffusion rate of work*/
        row[8][i] = rhoPSI * h * HDFrow[w][n] * C;
        RoW[8][n] = rhoPSI * h * HDFrow[w][n];
/* *----------------------------------------------------------------------
        energy advection*/
        row[4][i] = rhoPSI * h * ADVrow[w][n] * C;
        RoW[4][n] = rhoPSI * h * ADVrow[w][n];
        }

/*----------------------------------------------------------------------
      spatial integration*/
      for(k = 0; k < gauss_n; k++) {
        fe_LGP1base(gauss_x[k], gauss_y[k], beta);
        for(j = 0; j < 10; j++) {
          zr = 0.;
          for(i = 0; i < 3; i++) {
            zr += row[j][i] * beta[i];
            }
          tmp[j] += gauss_w[k] * zr;
          }
        }
      for(j = 0; j < 10; j++) {
        sum[j]    += 2 * A * tmp[j];
        local[j][m] =2 * A * tmp[j];
        }
      }   /* end of triangle loop */

    printf("wave %s\n", wave.name);
    printf("astronom. pot.  rate of work (Gw): %f\n", sum[0] / 1.e+09);
    printf("secondary pot.  rate of work (Gw): %f\n", sum[1] / 1.e+09);
    printf("wave drag       rate of work (Gw): %f\n", sum[2] / 1.e+09);
    printf("pressure forces rate of work (Gw): %f\n", sum[3] / 1.e+09);
    printf("pressure forces rate of work (Gw): %f\n", sum[5] / 1.e+09);
    printf("energy advection rate        (Gw): %f\n", sum[4] / 1.e+09);
    printf("bottom friction rate of work (Gw): %f\n", sum[6] / 1.e+09);
    printf("bottom friction rate of work (Gw): %f\n", sum[7] / 1.e+09);
    printf("horizontal dif. rate of work (Gw): %f\n", sum[8] / 1.e+09);
    printf("Tidal forcing rate of work   (Gw): %f\n", (sum[0]+sum[1]+sum[3]-sum[4]) / 1.e+09);
    printf("Dissipation rate of work     (Gw): %f\n", (sum[2]+sum[7]+sum[8]) / 1.e+09);
    printf("Energy balance closure       (Gw): %f\n", (sum[0]+sum[1]+sum[2]+sum[3]-sum[4]+sum[7]+sum[8]) / 1.e+09);

    sprintf(ncfilename, "%s/%s.Cd.nc" ,path,  wave.name);
    fprintf(stderr,"line %d of %s: writing to %s\n",__LINE__,__FILE__,ncfilename);
    
    sprintf(filename, "%s/%s.row-pot.%2.2d.s2r",output_path,  wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], " (w/m²)");
    status=quoddy_saver1(filename, mesh.nvtxs, RoW[0],comment);
    
    sprintf(filename, "%s/%s.row-lsa.%2.2d.s2r", output_path, wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], " (w/m²)");
    status=quoddy_saver1(filename, mesh.nvtxs, RoW[1],comment);

    sprintf(filename, "%s/%s.row-wd.%2.2d.s2r", output_path, wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], " (w/m²)");
    status=quoddy_saver1((const char *) filename, mesh.nvtxs, RoW[2],comment);

    //pressure gradient
    sprintf(filename, "%s/%s.prs-grd.%2.2d.v2c", output_path, wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], " (m/m ?)");
    status=quoddy_savec2((const char *) filename, mesh.nvtxs, Px,Py,comment);
    //RoW
    sprintf(filename, "%s/%s.row-prs.%2.2d.s2r", output_path, wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], " (w/m²)");
    status=quoddy_saver1((const char *) filename, mesh.nvtxs, RoW[3],comment);

/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Notes:

  29/08/2009:
    retrieve friction parameters from energy budget
    map of major energy terms to do asap

    pressure + ... RoW=friction RoW=-rho Cd |u| u�

    Cd = [0.4/(log(H/z0)-1)]�

    z0=h/exp(1+O.4/sqrt(Cd))

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/* *----------------------------------------------------------------------
    Real Cd */
    for(n = 0; n < mesh .nvtxs; n++) {
      buf[n]    = data2D.Cd[n];
      }
    archiving_UGdummy2D(ncfilename, mesh, "Cd", "none", buf, mask, LGP1);

/* *----------------------------------------------------------------------
    Cd inversion, rough estimate of dissipation work*/
    for(n = 0; n < mesh .nvtxs; n++) {
      a2 = data2D.Utide[n].a[w];
      a3 = data2D.Vtide[n].a[w];
      RoW[9][n] = RoW[3][n]+RoW[0][n]+RoW[1][n]+RoW[2][n]-RoW[4][n]+RoW[8][n];
    #if 1
      uu = sqrt(a2*a2+a3*a3);
      buf[n]    = RoW[9][n]/(uu*uu*uu*rhoPSI);
    #else
      uu = a2*a2+a3*a3;
      buf[n]    = RoW[9][n]/(uMax[n]*uu*.5*rhoPSI);
    #endif
      }
    sprintf(filename, "%s/%s.prs-equivalent-Cd-a.%2.2d.s2r", output_path, wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], "equivalent Cd using rough estimate of dissipation work");
    status=quoddy_saver1(filename, mesh.nvtxs, buf,comment);
    archiving_UGdummy2D(ncfilename, mesh, "Cd_a", "none", buf, mask, LGP1);

/* *----------------------------------------------------------------------
    z0 inversion*/
    for(n = 0; n < mesh .nvtxs; n++) {
      uu=MAX(0.0,buf[n]);
      uu=sqrt(uu);
      uu=exp(1.+0.4/uu);
      buf[n] = data2D.h[n]/uu;
      }
    sprintf(filename, "%s/%s.prs-equivalent-z0-a.%2.2d.s2r", output_path, wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], "equivalent z0");
    status=quoddy_saver1(filename, mesh.nvtxs, buf,comment);

/* *----------------------------------------------------------------------
    z0 inversion, should be identical to z0 for hydrodynamical solutions*/
    for(n = 0; n < mesh .nvtxs; n++) {
      uu=MAX(0.0,data2D.Cd[n]);
      uu=sqrt(uu);
      uu=exp(1.+0.4/uu);
      buf[n] = data2D.h[n]/uu;
      }
    sprintf(filename, "%s/%s.prs-equivalent-z0-true.%2.2d.s2r", output_path, wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], "equivalent z0");
    status=quoddy_saver1(filename, mesh.nvtxs, buf,comment);

/* *----------------------------------------------------------------------
    Cd inversion, linearised estimate of dissipation work*/
    for(n = 0; n < mesh .nvtxs; n++) {
      dum= -RoW[6][n]/data2D.Cd[n];
      buf[n] = RoW[9][n]/dum;
      }
    sprintf(filename, "%s/%s.prs-equivalent-Cd-b.%2.2d.s2r", output_path, wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], "equivalent Cd using linearised estimate of dissipation work");
    status=quoddy_saver1(filename, mesh.nvtxs, buf,comment);
    archiving_UGdummy2D(ncfilename, mesh, "Cd_b", "none", buf, mask, LGP1);

/* *----------------------------------------------------------------------
    z0 inversion*/
    for(n = 0; n < mesh .nvtxs; n++) {
      uu=MAX(0.0,buf[n]);
      uu=sqrt(uu);
      uu=exp(1.+0.4/uu);
      buf[n] = data2D.h[n]/uu;
      }

    sprintf(filename, "%s/%s.prs-equivalent-z0-b.%2.2d.s2r", output_path, wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], "equivalent z0");
    status=quoddy_saver1(filename, mesh.nvtxs, buf,comment);

/* *----------------------------------------------------------------------
    Cd inversion, analysed estimate of dissipation work*/
    for(n = 0; n < mesh .nvtxs; n++) {
      dum= -RoW[7][n]/data2D.Cd[n];
      buf[n] = RoW[9][n]/dum;
      }
    sprintf(filename, "%s/%s.prs-equivalent-Cd-c.%2.2d.s2r", output_path, wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], "equivalent Cd using analysed estimate of dissipation work");
    status=quoddy_saver1(filename, mesh.nvtxs, buf,comment);
    archiving_UGdummy2D(ncfilename, mesh, "Cd_c", "none", buf, mask, LGP1);

/* *----------------------------------------------------------------------
    z0 improvement*/
    for(n = 0; n < mesh .nvtxs; n++) {
      uu=MAX(0.0,buf[n]);
      uu=sqrt(uu);
      uu=exp(1.+0.4/uu);
      uu = data2D.h[n]/uu;
      rate=100*fabs(RoW[9][n])/(fabs(RoW[0][n])+fabs(RoW[1][n])+fabs(RoW[2][n])+fabs(RoW[3][n])+fabs(RoW[4][n])+fabs(RoW[8][n]));
      if((rate>gratio[n])&& (rate>75)&& (fabs(RoW[9][n])>1.e-03) && (data2D.h[n]<200.)) {
        z0[n]=uu;
        }
      gratio[n]=MAX(gratio[n],rate);
      gvalue[n]=MAX(gvalue[n],fabs(RoW[9][n]));
      if((rate>60.0) && (fabs(RoW[9][n])>1.e-03)) {
        buf[n]=uu;
        }
      else {
        buf[n]=0.002;
        }
      if(data2D.h[n]>200.) buf[n]=0.002;
      }
    smooth_vector( mesh, buf);
    smooth_vector( mesh, buf);
    smooth_vector( mesh, buf);

    sprintf(filename, "%s/%s.prs-equivalent-z0-c.%2.2d.s2r", output_path, wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], "equivalent z0");
    status=quoddy_saver1(filename, mesh.nvtxs, buf,comment);

/* *----------------------------------------------------------------------
    Cd reconstruction*/
    for(n = 0; n < mesh .nvtxs; n++) {
      uu=data2D.h[n]/buf[n];
      uu=0.4/(log(uu)-1);
      uu=uu*uu;
      buf[n] = uu;
      }

    sprintf(filename, "%s/%s.prs-equivalent-Cd-d.%2.2d.s2r", output_path, wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], "equivalent Cd using analysed estimate of dissipation work");
    status=quoddy_saver1(filename, mesh.nvtxs, buf,comment);
    archiving_UGdummy2D(ncfilename, mesh, "Cd_d", "none", buf, mask, LGP1);
    for(n = 0; n < mesh .nvtxs; n++) {
      buf[n] *= data2D.h[n];
      }
    archiving_UGdummy2D(ncfilename, mesh, "hCd_d", "m", buf, mask, LGP1);

    for(n = 0; n < mesh .nvtxs; n++) {
      buf[n] = 100*fabs(RoW[9][n])/(fabs(RoW[0][n])+fabs(RoW[1][n])+fabs(RoW[2][n])+fabs(RoW[3][n])+fabs(RoW[4][n])+fabs(RoW[8][n]));
      }

    sprintf(filename, "%s/%s.absolute-budget.%2.2d.s2r", output_path, wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], "absolute budget");
    status=quoddy_saver1(filename, mesh.nvtxs, buf,comment);

/* *----------------------------------------------------------------------
    bottom friction*/
    sprintf(filename, "%s/%s.row-bfr-linearised.%2.2d.s2r", output_path, wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], " (w/m²)");
    status=quoddy_saver1(filename, mesh.nvtxs, RoW[6],comment);
    archiving_UGdummy2D(ncfilename, mesh, "bfRoW_lin", "w", RoW[6], mask, LGP1);

    sprintf(filename, "%s/%s.row-bfr-analysed.%2.2d.s2r", output_path, wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], " (w/m²)");
    status=quoddy_saver1(filename, mesh.nvtxs, RoW[7],comment);
    archiving_UGdummy2D(ncfilename, mesh, "bfRoW_ana", "w", RoW[7], mask, LGP1);

    for(n = 0; n < mesh .nvtxs; n++) {
      buf[n] = -RoW[9][n];
      }
    sprintf(filename, "%s/%s.row-bfr-deduced.%2.2d.s2r", output_path, wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], " (w/m²)");
    status=quoddy_saver1(filename, mesh.nvtxs, buf,comment);
    archiving_UGdummy2D(ncfilename, mesh, "bfRoW_ded", "w", buf, mask, LGP1);

/* *----------------------------------------------------------------------
    horizontal diffusion*/
    sprintf(filename, "%s/%s.row-hdf.%2.2d.s2r", output_path, wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], " (w/m²)");
    status=quoddy_saver1(filename, mesh.nvtxs, RoW[8],comment);

    sprintf(filename, "%s/%s.row-adv.%2.2d.s2r", output_path, wave.name, idx);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], " (w/m²)");
    status=quoddy_saver1(filename, mesh.nvtxs, RoW[4],comment);

    sprintf(filename, "%s/%s.budget.%2.2d.s2r", output_path, wave.name, idx);
    out = fopen(filename, "w");
    fprintf(out, "%s\n", "off-line");
    fprintf(out, " (w/m²)\n");
    for(n = 0; n < mesh .nvtxs; n++) {
      for(j=0;j<mesh .vertices[n].nelmts;j++) {
        m=mesh .vertices[n].elmts[j];
        zr=local[0][m]+local[1][m]+local[2][m]+local[3][m]-local[4][m]+local[7][m]+local[8][m];
        }
      fprintf(out, "%d %lf\n", n + 1, zr/mesh .vertices[n].mw);
      }
    fclose(out);

/* *----------------------------------------------------------------------
    Total kinetic energy*/
    for(n = 0; n < mesh.nvtxs; n++) {
      a2 = data2D.Utide[n].a[w];
      a3 = data2D.Vtide[n].a[w];
      meanKE[n]+=(a2*a2+a3*a3);
      }

/* *----------------------------------------------------------------------
    Mean wavedrag energy*/
    for(n = 0; n < mesh.nvtxs; n++) {
      meanWD[n]+=RoW[2][n];
      }

    }

  sprintf(filename, "%s/gratio.%2.2d.s2r", output_path, idx);
  sprintf(comment[0], "%s", "off-line");
  sprintf(comment[1], "ratio");
  status=quoddy_saver1(filename, mesh.nvtxs, gratio,comment);

  sprintf(filename, "%s/gvalue.%2.2d.s2r", output_path, idx);
  sprintf(comment[0], "%s", "off-line");
  sprintf(comment[1], "value");
  status=quoddy_saver1(filename, mesh.nvtxs, gvalue,comment);

  sprintf(filename, "%s/totalKE.s2r", output_path, idx);
  sprintf(comment[0], "%s", "off-line");
  sprintf(comment[1], "total kinetic energy");
  status=quoddy_saver1(filename, mesh.nvtxs, meanKE,comment);

  sprintf(filename, "%s/totalWD.s2r", output_path, idx);
  sprintf(comment[0], "%s", "off-line");
  sprintf(comment[1], "total kinetic energy");
  status=quoddy_saver1(filename, mesh.nvtxs, meanWD,comment);

  for(n = 0; n < mesh .nvtxs; n++) {
    uu=data2D.h[n]/z0[n];
    uu=0.4/(log(uu)-1);
    uu=uu*uu;
    buf[n] = uu;
    }

  sprintf(filename, "%s/equivalent-Cd-a.%2.2d.s2r", output_path, idx);
  sprintf(comment[0], "%s", "off-line");
  sprintf(comment[1], "equivalent Cd using analysed estimate of dissipation work");
  status=quoddy_saver1(filename, mesh.nvtxs, buf,comment);

  sprintf(filename, "%s/equivalent-z0-a.%2.2d.s2r", output_path, idx);
  sprintf(comment[0], "%s", "off-line");
  sprintf(comment[1], "equivalent z0 using analysed estimate of dissipation work");
  status=quoddy_saver1(filename, mesh.nvtxs, z0,comment);

  smooth_vector( mesh, z0);
  smooth_vector( mesh, z0);
  smooth_vector( mesh, z0);

  for(n = 0; n < mesh .nvtxs; n++) {
    uu=data2D.h[n]/z0[n];
    uu=0.4/(log(uu)-1);
    uu=uu*uu;
    buf[n] = uu;
    }

  sprintf(filename, "%s/equivalent-Cd-b.%2.2d.s2r", output_path, idx);
  sprintf(comment[0], "%s", "off-line");
  sprintf(comment[1], "equivalent Cd using analysed estimate of dissipation work");
  status=quoddy_saver1(filename, mesh.nvtxs, buf,comment);

  sprintf(filename, "%s/equivalent-z0-b.%2.2d.s2r", output_path, idx);
  sprintf(comment[0], "%s", "off-line");
  sprintf(comment[1], "equivalent z0 using analysed estimate of dissipation work");
  status=quoddy_saver1(filename, mesh.nvtxs, z0,comment);

  for(k = 0; k < 10; k++)
    delete[] RoW[k];
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void spectral_transport(mesh_t & mesh, spectrum_t & AnalysisList, parameter_t & data2D, const char *output_path, char *extension)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, m, w, status;
  double uu, u, v, A, U, gPSI, rhoPSI, C;
  double a1, a2, G1, G2, z;
  double *meanx, *meany,*divergence;
  double h, H, dH;
  double factor = 1.e+06 * 3600.0;
  int n;
  char filename[1024];
  FILE *out;
  tidal_wave wave;
  char *comment[2];

/* *----------------------------------------------------------------------
  Compute and save spectral mass and energy transport*/

  meanx = new double[mesh .nvtxs];
  meany = new double[mesh .nvtxs];

  comment[0]=new char[256];
  comment[1]=new char[256];

  gPSI = P_g;
  rhoPSI = rho_water;

  for(k = 0; k < AnalysisList.n; k++) {
    wave = AnalysisList.waves[k];
    if(wave.omega != 0.0) {
/* *----------------------------------------------------------------------
     usual tidal harmonic transport*/
      for(n = 0; n < mesh .nvtxs; n++) {
        a1 = data2D.Htide[n].a[k];
        G1 = data2D.Htide[n].G[k];
        a2 = data2D.Utide[n].a[k];
        G2 = data2D.Utide[n].G[k];
        meanx[n] = 0.5 * a1 * a2 * cos(G2 - G1);
        a2 = data2D.Vtide[n].a[k];
        G2 = data2D.Vtide[n].G[k];
        meany[n] = 0.5 * a1 * a2 * cos(G2 - G1);
        }
      }
    else {
/* *----------------------------------------------------------------------
     tidal residual transport*/
      for(n = 0; n < mesh .nvtxs; n++) {
        a1 = data2D.Htide[n].a[k];
        G1 = data2D.Htide[n].G[k];
        a2 = data2D.Utide[n].a[k];
        G2 = data2D.Utide[n].G[k];
        meanx[n] = a2 * cos(G2) * (a1 * cos(G1) + data2D.h[n]);
        a2 = data2D.Vtide[n].a[k];
        G2 = data2D.Vtide[n].G[k];
        meany[n] = a2 * cos(G2) * (a1 * cos(G1) + data2D.h[n]);
        }
      }

    sprintf(filename, "%s/%s.mass-flux.%s.v2r", output_path, wave.name, extension);
    sprintf(comment[0], "%s", "off-line");
    sprintf(comment[1], "mean spectral transport (m^3/s/m)");
    status=quoddy_saver2(filename, mesh.nvtxs, meanx, meany,comment);
    }

  for(k = 0; k < AnalysisList.n; k++) {
    wave = AnalysisList.waves[k];
    if(wave.omega != 0.0) {
/*----------------------------------------------------------------------
      usual tidal harmonic energy flux, corrected 18/08/2005*/
      for(n = 0; n < mesh .nvtxs; n++) {
        a1 = data2D.Htide[n].a[k];
        G1 = data2D.Htide[n].G[k];
        a2 = data2D.Utide[n].a[k];
        G2 = data2D.Utide[n].G[k];
        meanx[n] = 9.81 * rhoPSI * 0.5 * a1 * a2 * cos(G2 - G1) * data2D.h[n];
        a2 = data2D.Vtide[n].a[k];
        G2 = data2D.Vtide[n].G[k];
        meany[n] = 9.81 * rhoPSI * 0.5 * a1 * a2 * cos(G2 - G1) * data2D.h[n];
        }
      sprintf(filename, "%s/%s.energy-flux.%s.v2r", output_path, wave.name, extension);
      sprintf(comment[0], "%s", "off-line");
      sprintf(comment[1], "mean spectral energy flux (W/m)");
      status=quoddy_saver2(filename, mesh.nvtxs, meanx, meany,comment);
//       for(n = 0; n < mesh .nvtxs; n++) {
//         meanx[n] = 1.0;
//         meany[n] = 0.0;
//         meanx[n] = 0.0;
//         meany[n] = 1.0;
//         meanx[n] = R*mesh.vertices[n].lon*d2r;
//         meany[n] = 0.0;
//         }
      divergence=energy_divergence( mesh, meanx, meany);
      sprintf(filename, "%s/%s.energy-divergence.%s.s2r", output_path, wave.name, extension);
      sprintf(comment[0], "%s", "off-line");
      sprintf(comment[1], "mean spectral energy divergence (W/m²)");
      status=quoddy_saver1(filename, mesh.nvtxs, divergence,comment);
/*----------------------------------------------------------------------
     div(E)=rho Cd |u| u²*/
      for(n = 0; n < mesh .nvtxs; n++) {
        a1 = data2D.Utide[n].a[k];
        a2 = data2D.Vtide[n].a[k];
        uu = sqrt(a1*a1+a2*a2);
        meanx[n] = divergence[n]/uu/uu/uu/rhoPSI;
        }
      sprintf(filename, "%s/%s.flux-equivalent-Cd.%s.s2r", output_path, wave.name, extension);
      sprintf(comment[0], "%s", "off-line");
      sprintf(comment[1], "equivalent Cd");
      status=quoddy_saver1(filename, mesh.nvtxs, meanx,comment);
      }
    }

  delete[] meanx;
  delete[] meany;
  delete[] comment[0];
  delete[] comment[1];
}
