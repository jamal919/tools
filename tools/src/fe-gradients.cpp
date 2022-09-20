
/*******************************************************************************

  T-UGOm hydrodynamic ocean model, 2006

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

*******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <cmath>

#include "tools-structures.h"
#include "fe.h"
#include "geo.h"
#include "map.h"
#include "constants.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_gradientLGP1(mesh_t mesh, double *z, int m, double *gx, double *gy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  compute local (element) gradient of a LGP1 filed
-----------------------------------------------------------------------------*/
{
  int i, n, p, status = 0;
  double x = 0., y = 0., C = 0.;
  double dzdx = 0., dzdy = 0.;
  double beta_x[3], beta_y[3];

  dzdx = 0.;
  dzdy = 0.;
  for(i = 0; i < 3; i++) {
    n = mesh.triangles[m].vertex[i];
    C += cos(mesh.vertices[n].lat * M_PI / 180.0) / (double) 3.;
    dzdx += z[n] * mesh.triangles[m].DQ[i];
    dzdy -= z[n] * mesh.triangles[m].DP[i];
    }
  *gx = dzdx / mesh.triangles[m].Area / (double) 2.0 /C;
  *gy = dzdy / mesh.triangles[m].Area / (double) 2.0;

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int fe_gradientLGP1_template(mesh_t mesh, discretisation_t descriptor, T *z, int m, T *gx, T *gy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  compute local (element) gradient of a LGP1 field
-----------------------------------------------------------------------------*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
  double C = 0.;
  T dzdx = 0., dzdy = 0.;

  dzdx = 0.;
  dzdy = 0.;
  for(i = 0; i < 3; i++) {
    n  = descriptor.NIbE[m][i];
    C += descriptor.nodes[n].c / (double) 3.;
    dzdx += z[n] *(T)  mesh.triangles[m].DQ[i];
    dzdy -= z[n] *(T)  mesh.triangles[m].DP[i];
    }
  gx[0] = dzdx / (T) mesh.triangles[m].Area / (T) 2.0 /(T) C;
  gy[0] = dzdy / (T) mesh.triangles[m].Area / (T) 2.0;

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_gradientNCP1(mesh_t mesh, double *z, int m, double *gx,double *gy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
  returns the true gradient of NCP1 buffer*/
  int i, n, p, status = 0;
  double x = 0., y = 0., C = 0.;
  double dzdx = 0., dzdy = 0.;
  double gamma_x[3], gamma_y[3];

/*------------------------------------------------------------------------------
  returns gradient in z units/meters*/
  C = 0.;
  dzdx = 0.;
  dzdy = 0.;
  for(i = 0; i < 3; i++) {
    n = mesh.triangles[m].edges[i];
    C += mesh.edges[n].c /3.0;
    dzdx -= z[n] * mesh.triangles[m].DQ[i];
    dzdy += z[n] * mesh.triangles[m].DP[i];
    }
  *gx = dzdx / mesh.triangles[m].Area /C;
  *gy = dzdy / mesh.triangles[m].Area;

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int fe_gradientNCP1_template(mesh_t mesh, discretisation_t descriptor, T *z, int m, T *gx, T *gy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  derivation-related coefficient. Definition:

  1/R ds/dlon = dsdx / 2A
  1/R ds/dlat = dsdy / 2A

  dsdx= ds /dlon * 2A/R
  dsdy= ds /dlat * 2A/R

  for NCP1 variable:

  d beta[j] /d lon = triangle.dxdt * beta_x[j] + triangle.dydt * beta_y[j]
                   = -R DQ[j]/A
  d beta[j] /d lat = triangle.dxdp * beta_x[j] + triangle.dydp * beta_y[j]
                   = +R DP[j]/A

  dsdx= -2*sum( DQ[j] * s(j) )
  dsdy= +2*sum( DP[j] * s(j) )

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*------------------------------------------------------------------------------
  returns the true gradient of NCP1 buffer*/
  int i, n, p, status = 0;
  double C = 0.;
  T dzdx = 0., dzdy = 0.;

/*------------------------------------------------------------------------------
  returns gradient in z units/meters*/
  C = 0.;
  dzdx = 0.;
  dzdy = 0.;
  for(i = 0; i < 3; i++) {
    n  = descriptor.NIbE[m][i];
    C += descriptor.nodes[n].c / (double) 3.;
    dzdx -= z[n] * (T) mesh.triangles[m].DQ[i];
    dzdy += z[n] *(T)  mesh.triangles[m].DP[i];
    }
  gx[0] = dzdx / (T) mesh.triangles[m].Area /(T) C;
  gy[0] = dzdy / (T) mesh.triangles[m].Area;

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int fe_gradientLGP2_template(mesh_t mesh, discretisation_t descriptor, T *z, int m, T *gx, T *gy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  compute local (element) gradient of a LGP2 field
-----------------------------------------------------------------------------*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  derivation-related coefficient. Definition:

  for LGP1 variable:

  d beta[j] /d lon = triangle.dxdt * beta_x[j] + triangle.dydt * beta_y[j]

  d beta[j] /d lat = triangle.dxdp * beta_x[j] + triangle.dydp * beta_y[j]

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int i, k, n, p, status = 0;
//  double C = 0.;
  T C = 0.;
  T dzdx = 0., dzdy = 0.;
  double dxdt,dydt,dxdp,dydp,factor;
  double x,y,beta_x[6][6],beta_y[6][6];
  
//  factor=2.0*mesh.triangles[m].Area/R;
  factor=1.0;
  dxdt=mesh.triangles[m].dxdt*factor;
  dydt=mesh.triangles[m].dydt*factor;
  dxdp=mesh.triangles[m].dxdp*factor;
  dydp=mesh.triangles[m].dydp*factor;
  
//   x=0;y=0;k=0;
//   fe_LGP2prime( x, y, beta_x[k], beta_y[k]);
//   x=1;y=0;k=1;
//   fe_LGP2prime( x, y, beta_x[k], beta_y[k]);
//   x=0;y=1;k=2;
//   fe_LGP2prime( x, y, beta_x[k], beta_y[k]);

  x=0;   y=0;  k=0;
  fe_LGP2prime( x, y, beta_x[k], beta_y[k]);
  x=0.5; y=0;  k=1;
  fe_LGP2prime( x, y, beta_x[k], beta_y[k]);
  x=1;   y=0;  k=2;
  fe_LGP2prime( x, y, beta_x[k], beta_y[k]);
  x=0.5; y=0.5; k=3;
  fe_LGP2prime( x, y, beta_x[k], beta_y[k]);
  x=0;   y=1;   k=4;
  fe_LGP2prime( x, y, beta_x[k], beta_y[k]);
  x=0;   y=0.5; k=5;
  fe_LGP2prime( x, y, beta_x[k], beta_y[k]);

  for(k=0;k<6;k++) {
    dzdx = 0.;
    dzdy = 0.;
//    C=0;
    C=mesh.triangles[m].cosinus;
    for(i = 0; i < descriptor.nnpe; i++) {
      n  = descriptor.NIbE[m][i];
//      C += descriptor.nodes[n].c / (double) 6.;
      dzdx += z[n] * (T) (dxdt*beta_x[k][i]+dydt*beta_y[k][i]);
      dzdy += z[n] * (T) (dxdp*beta_x[k][i]+dydp*beta_y[k][i]);
      }
    gx[k] = dzdx / (T) R /C;
    gy[k] = dzdy / (T) R;
    }
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int fe_gradient_template(mesh_t mesh, discretisation_t descriptor, T *z, int m, T *gx, T *gy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  switch (descriptor.type) {
    case LGP1:
    case DGP1:
      status=fe_gradientLGP1_template(mesh,  descriptor, z, m,  gx, gy);
      break;
    case DNP1:
    case NCP1:
      status=fe_gradientNCP1_template(mesh,  descriptor,  z, m, gx, gy);
      break;
    case LGP2:
//    case DGP2:
      status=fe_gradientLGP2_template(mesh,  descriptor, z, m,  gx, gy);
      break;
    default:
      break;
    }
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_gradient(mesh_t mesh, discretisation_t descriptor, float *z, int m, float *gx, float *gy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_gradient_template(mesh,  descriptor,  z, m, gx, gy);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_gradient(mesh_t mesh, discretisation_t descriptor, double *z, int m, double *gx, double *gy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_gradient_template(mesh,  descriptor,  z, m, gx, gy);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_gradient(mesh_t mesh, discretisation_t descriptor, complex<float> *z, int m, complex<float> *gx, complex<float> *gy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_gradient_template(mesh,  descriptor,  z, m, gx, gy);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_gradient(mesh_t mesh, discretisation_t descriptor, complex<double> *z, int m, complex<double> *gx, complex<double> *gy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_gradient_template(mesh,  descriptor,  z, m, gx, gy);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_P1NCgradient_test(mesh_t mesh, double *z, int m, double *gx,double *gy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
  returns the true gradient of NCP1 buffer*/
  int i, n, p, status = 0;
  double x = 0., y = 0., C = 0.;
  double dzdx = 0., dzdy = 0.;
  double gamma_x[3], gamma_y[3];

  fe_NCP1prime(x, y, gamma_x, gamma_y);

  for(i = 0; i < 3; i++) {
    n = mesh.triangles[m].edges[i];
//    C+=cos(mesh.edges[n].lat*M_PI/180.0)/3.;
/*------------------------------------------------------------------------------
    Warning : assume latitude in radians*/
    C += cos(mesh.edges[n].lat) / 3.;
    dzdx += z[n] * gamma_x[i];
    dzdy += z[n] * gamma_y[i];
    }

/*------------------------------------------------------------------------------
  returns grdient in z units/meters*/
  *gx = (dzdx * mesh.triangles[m].dxdt + dzdy * mesh.triangles[m].dydt) / (R * C);
  *gy = (dzdx * mesh.triangles[m].dxdp + dzdy * mesh.triangles[m].dydp) / R;

  C = 0.;
  dzdx = 0.;
  dzdy = 0.;
  for(i = 0; i < 3; i++) {
    n = mesh.triangles[m].edges[i];
    C += cos(mesh.edges[n].lat) / 3.;
    dzdx -= z[n] * mesh.triangles[m].DQ[i];
    dzdy += z[n] * mesh.triangles[m].DP[i];
    }
  *gx = dzdx / mesh.triangles[m].Area /C;
  *gy = dzdy / mesh.triangles[m].Area;

  return (status);
}

