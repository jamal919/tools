
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
// #include "fe-integration.h"
#include "map.h"
#include "mass.h"

#include "cefmo.h"
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

const char *tugo_version="tools", *gOutputPath=".";

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//  template <typename T> T integrale_time_2(complex<T> z1, complex<T>  z2, T (*sintegration) (T*, T*))
  template <typename T> static T integrale_time_2(complex<T> z1, complex<T>  z2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Time average over 1 period
-----------------------------------------------------------------------------*/
{
  T d;
  d = 0.5 * (real(z1) * real(z2) + imag(z1) * imag(z2));
  return (d);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  complex<double> local_integraleLGP1xLGP0xLGP0_2D(double *p, complex<double> *q, complex<double> *r)
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n;
  complex<double> sum;

  sum = fe_integraleLGP1_2D(p);
  sum*=q[0]*r[0];
  return (sum);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  complex<double> local_integraleLGP0xLGP0xLGP0_2D(double *p, complex<double> *q, complex<double> *r)
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n;
  complex<double> sum;

  sum=0.5*p[0]*q[0]*r[0];
  return (sum);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  complex<double> local_integraleLGP1xNCP1xLGP0_2D(double *p, complex<double> *q, complex<double> *r)
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n;
  complex<double> sum;

  sum = fe_integraleLGP1xNCP1_2D(p,q)*r[0];
  return (sum);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  complex<double> local_integraleLGP2xLGP2xLGP1_2D(double *p, complex<double> *q, complex<double> *r)
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n;
  complex<double> sum;

  sum = fe_integraleLGP2xLGP2xLGP1_2D(p,q,r);
  return (sum);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_budget(FILE *out, const char *region, char *wavename, double *sum)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, m, status;
  double budget;
    
  fprintf(out,"%s wave %s\n", region, wavename);
  fprintf(out,"astronom. pot.  rate of work (Gw): %f\n", sum[0] / 1.e+09);
  fprintf(out,"secundary pot.  rate of work (Gw): %f\n", sum[1] / 1.e+09);
  fprintf(out,"wave drag       rate of work (Gw): %f\n", sum[2] / 1.e+09);
  fprintf(out,"pressure forces rate of work (Gw): %f\n", sum[5] / 1.e+09);
  fprintf(out,"bottom friction rate of work (Gw): %f\n", sum[6] / 1.e+09);
  budget=sum[0]+sum[1]+sum[2]+sum[5]-sum[4]+sum[6]+sum[8];
  fprintf(out,"Energy balance closure       (Gw): %f\n", budget / 1.e+09);
  fprintf(out,"Mean Kinetic Energy          (TJ): %f\n", sum[10] / 1.e+12);
  fprintf(out,"Mean Potential Energy        (TJ): %f\n", sum[11] / 1.e+12);
  fprintf(out,"Q factor           (adimensional): %f\n", sum[12]);
  
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int SpArchive_Get2D(const char *filename, mesh_t & mesh, int frame, cefmo_t & cefmo, tide2D_t & state)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------

  Spectral solution archive

 -----------------------------------------------------------------------*/
{
  int i,j,k,l,m,n,status;
  int ncid,varid;
  int zdim,udim;
  int z_discretisation, u_discretisation;
  int wave_id;
  double *buffer[2];
  char varname_a[1024], varname_G[1024];

  status=paire_discretisation_id(cefmo.paire, &z_discretisation, &u_discretisation);
  if(status!=0) return(-1);
  
  const char *UNAME=discretisation_name(u_discretisation);
  const char *ZNAME=discretisation_name(z_discretisation);

/**----------------------------------------------------------------------------
  elevation */
  sprintf(varname_a,"a_eta_%s", ZNAME);
  sprintf(varname_G,"G_eta_%s", ZNAME);
  status=poc_get_UG3D(filename, frame, (const char *) varname_a, (const char *) varname_G, state.z);
  if(status!=0) goto error;
  
/**----------------------------------------------------------------------------
  currents */
  sprintf(varname_a,"a_u_%s", UNAME);
  sprintf(varname_G,"G_u_%s", UNAME);
  status=poc_get_UG3D(filename, frame, (const char *) varname_a, (const char *) varname_G, state.u);
  if(status!=0) goto error;
  sprintf(varname_a,"a_v_%s", UNAME);
  sprintf(varname_G,"G_v_%s", UNAME);
  status=poc_get_UG3D(filename, frame, (const char *) varname_a, (const char *) varname_G, state.v);
  if(status!=0) goto error;

  return(0);

error:
  printf("%s read error: v-amplitude=%s v-phaselag=%s frame=%d\n",__FUNCTION__, varname_a, varname_G, frame);
  return(-1);

}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int SpArchive_GetFull2D(const char *filename, mesh_t & mesh, int frame, cefmo_t & cefmo, tide2D_t & state)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------

  Spectral solution archive

 -----------------------------------------------------------------------*/
{
  int i,j,k,l,m,n,status;
  int ncid,varid;
  int zdim,udim;
  int z_discretisation, u_discretisation;
  int wave_id;
  double *buffer[2];
  char varname[1024], varname_a[1024], varname_G[1024];

  status=paire_discretisation_id(cefmo.paire, &z_discretisation, &u_discretisation);
  if(status!=0) return(-1);
  
  const char *UNAME=discretisation_name(u_discretisation);
  const char *ZNAME=discretisation_name(z_discretisation);

/**----------------------------------------------------------------------------
  elevation */
  sprintf(varname_a,"a_eta_%s", ZNAME);
  sprintf(varname_G,"G_eta_%s", ZNAME);
  status=poc_get_UG3D(filename, frame, (const char *) varname_a, (const char *) varname_G, state.z);
  if(status!=0) goto error;
  
/**----------------------------------------------------------------------------
  potential */
  sprintf(varname_a,"a_POT_%s", ZNAME);
  sprintf(varname_G,"G_POT_%s", ZNAME);
  status=poc_get_UG3D(filename, frame, (const char *) varname_a, (const char *) varname_G, state.z);
  if(status!=0) goto error;
  
/**----------------------------------------------------------------------------
  LSA */
  sprintf(varname_a,"a_LSA_%s", ZNAME);
  sprintf(varname_G,"G_LSA_%s", ZNAME);
  status=poc_get_UG3D(filename, frame, (const char *) varname_a, (const char *) varname_G, state.z);
  if(status!=0) goto error;
  
/**----------------------------------------------------------------------------
  depth */
  sprintf(varname,"h_%s", ZNAME);
  status=poc_get_UG3D(filename, frame, (const char *) varname, state.h_z);
  if(status!=0) goto error;
  
/**----------------------------------------------------------------------------
  currents */
  sprintf(varname_a,"a_u_%s", UNAME);
  sprintf(varname_G,"G_u_%s", UNAME);
  status=poc_get_UG3D(filename, frame, (const char *) varname_a, (const char *) varname_G, state.u);
  if(status!=0) goto error;
  sprintf(varname_a,"a_v_%s", UNAME);
  sprintf(varname_G,"G_v_%s", UNAME);
  status=poc_get_UG3D(filename, frame, (const char *) varname_a, (const char *) varname_G, state.v);
  if(status!=0) goto error;

  return(0);

error:
  printf("%s read error: v-amplitude=%s v-phaselag=%s frame=%d\n",__FUNCTION__, varname_a, varname_G, frame);
  return(-1);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int spectral_energy_T(mesh_t & mesh, tidal_wave wave, const cefmo_t & cefmo, const parameter_t & data,int iteration, double *RoW[10], double *MKE, double *MPE, double *Q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  F. Lyard 21/02/2016: embarassing bug found on energy quantics !!!
  
  Energy density fields (RoW, w/m²) were 1/2 the values they should be, due to 
  triangle integration that needs to be multiply by a factor 2 (as unitary triangle
  on which numerical integrale are performed has an area of 1/2.); 
  global/regional budget were ok as they got multiplied in summation 
  (i.e. 2*A*C = 2.0 x (triangle) area)
  
  Bug fixed
  
------------------------------------------------------------------------------*/
{
  int i, j, k, m, n, w, status;
  double A, gPSI, rhoPSI, C;
  double zz;
  complex<double> uu,vv;
  double sum[100], budget;
//  double *RoW[10];
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int gauss_n = 7;
  int idx = 0;
  FILE *out;
  char *sdate;
  char *comments[2];
  char ascii_file[256];
  double mask=9999.;
  int fmt=MESH_FAMILY_QUODDY;

  complex<double> (*fe_integraleZxZ_2D)    (complex<double> *, complex<double> *);
  complex<double> (*fe_integraleZxUxGZ_2D) (double*, complex<double> *, complex<double> *);
  complex<double> (*fe_integraleZxUxU_2D)  (double*, complex<double> *, complex<double> *);
  complex<double> (*fe_integraleUxUxU_2D)  (double*, complex<double> *, complex<double> *);
  
  int zdim, udim;
  int z_discretisation, u_discretisation, g_discretisation;
  
  discretisation_t z_descriptor;
  discretisation_t u_descriptor;
  discretisation_t g_descriptor;

  size_t   Z_NNPE, U_NNPE, GZ_NNPE;
  complex<double> *px,*py;
  complex<double> *taux,*tauy;
  complex<double> *z, *u, *v, *zstar, *ustar, *vstar;
  
  double mesh_area=0.0,omega=wave.omega*dph2rps;
  double T=2*M_PI/omega;
  
  double *h,*Cd;
  double sumx,sumy;
  
  status=paire_discretisation_id(cefmo.paire, &z_discretisation,&u_discretisation);
  g_discretisation=gradient_natural_discretisation(z_discretisation);

  z_descriptor= get_descriptor(mesh, z_discretisation);
  u_descriptor= get_descriptor(mesh, u_discretisation);
  g_descriptor= get_descriptor(mesh, g_discretisation);
  
  status=paire_dimension(mesh,cefmo.paire,&zdim,&udim);
  
  Z_NNPE  = z_descriptor.nnpe;
  U_NNPE  = u_descriptor.nnpe;
  GZ_NNPE = g_descriptor.nnpe;

  h=new double[Z_NNPE];

  Cd=new double[U_NNPE];

  z=new complex<double>[Z_NNPE];
  u=new complex<double>[U_NNPE];
  v=new complex<double>[U_NNPE];
  
  zstar=new complex<double>[Z_NNPE];
  ustar=new complex<double>[U_NNPE];
  vstar=new complex<double>[U_NNPE];
  
  taux=new complex<double>[U_NNPE];
  tauy=new complex<double>[U_NNPE];
  
  px=new complex<double>[GZ_NNPE];
  py=new complex<double>[GZ_NNPE];
  
  gPSI   = P_g;
  rhoPSI = rho_water;
  
  switch (cefmo.paire) {
    
    case LGP0xLGP1:
      fe_integraleZxZ_2D    = fe_integraleLGP1xLGP1_2D;
      fe_integraleZxUxGZ_2D = local_integraleLGP1xLGP0xLGP0_2D;
      fe_integraleZxUxU_2D  = local_integraleLGP0xLGP0xLGP0_2D;
      fe_integraleUxUxU_2D  = local_integraleLGP0xLGP0xLGP0_2D;
      break;
      
    case LGP1xLGP1:
      fe_integraleZxUxGZ_2D = fe_integraleLGP1xLGP1xLGP0_2D;
      break;
      
    case NCP1xLGP1:
      fe_integraleZxZ_2D    = fe_integraleLGP1xLGP1_2D;
      fe_integraleZxUxGZ_2D = local_integraleLGP1xNCP1xLGP0_2D;
      fe_integraleZxUxU_2D  = fe_integraleLGP1xNCP1xNCP1_2D;
      fe_integraleUxUxU_2D  = fe_integraleNCP1xNCP1xNCP1_2D;
      break;
      
    case DGP1xLGP2:
      fe_integraleZxZ_2D    = fe_integraleLGP1xLGP1_2D;
      fe_integraleZxUxGZ_2D = fe_integraleLGP2xLGP1xLGP1_2D;
      fe_integraleZxUxU_2D  = fe_integraleLGP2xLGP1xLGP1_2D;
      fe_integraleUxUxU_2D  = fe_integraleLGP1xLGP1xLGP1_2D;
      break;
      
    case DNP1xLGP2:
      fe_integraleZxZ_2D    = fe_integraleLGP1xLGP1_2D;
      fe_integraleZxUxGZ_2D = fe_integraleLGP2xLGP1xLGP1_2D;
      fe_integraleZxUxU_2D  = fe_integraleLGP2xLGP1xLGP1_2D;
      fe_integraleUxUxU_2D  = fe_integraleLGP1xLGP1xLGP1_2D;
      break;
      
    case LGP2xLGP2:
      fe_integraleZxZ_2D    = fe_integraleLGP2xLGP2_2D;
      fe_integraleZxUxGZ_2D = fe_integraleLGP2xLGP2xLGP1_2D;
      fe_integraleZxUxU_2D  = fe_integraleLGP2xLGP2xLGP2_2D;
      fe_integraleUxUxU_2D  = fe_integraleLGP2xLGP2xLGP2_2D;
      break;
      
    case CQN1xCQP0:
//       fe_integraleZxZ_2D    = fe_integraleLGP2xLGP2_2D;
//       fe_integraleZxUxGZ_2D = fe_integraleLGP2xLGP2xLGP1_2D;
//       fe_integraleZxUxU_2D  = fe_integraleLGP2xLGP2xLGP2_2D;
//       fe_integraleUxUxU_2D  = fe_integraleLGP2xLGP2xLGP2_2D;
      return(-1);
      break;
      
    default:
      check_error(-1, "illicit discretisation paire", __LINE__, __FILE__, 1);
      break;
    }

  gauss_init(gauss_n, gauss_x, gauss_y, gauss_w);

/**----------------------------------------------------------------------
    to complete */
//     spectral_flux02(mesh, data, w);
//     spectral_flux01( mesh, w, gZspec[0], gNobc_Z[0]);

  for(j = 0; j < 100; j++)
    sum[j] = 0.0;

/**----------------------------------------------------------------------
  compute energy estimates on triangles*/
  for(m = 0; m < mesh.ntriangles; m++) {
/*----------------------------------------------------------------------
    convert A (element area) in m**2*/
    A = mesh.triangles[m].Area;
    C = mesh.triangles[m].cosinus;
    double area=A*C;
    for(k = 0; k < 10; k++) {
      RoW[k][m] = 0;
      }
    for(i = 0; i < z_descriptor.nnpe; i++) {
      n=z_descriptor.NIbE[m][i];
/*----------------------------------------------------------------------
      convert dpth in m*/
      h[i]=cefmo.state.h_z[n];
      z[i]=cefmo.state.z[n];
      zstar[i]=conj(cefmo.state.z[n]);
      }
    for(i = 0; i < u_descriptor.nnpe; i++) {
      n=u_descriptor.NIbE[m][i];
      u[i]=cefmo.state.u[n];
      v[i]=cefmo.state.v[n];
/**----------------------------------------------------------------------
      use conjugate of u,v to optimize time integration over 1 period*/
      ustar[i]=conj(cefmo.state.u[n]);
      vstar[i]=conj(cefmo.state.v[n]);
      }
/**----------------------------------------------------------------------
    energy fluxes*/
//       sumx=real(fe_integraleLGP1xLGP1xLGP1_2D(h,z,u))*C;
//       sumy=real(fe_integraleLGP1xLGP1xLGP1_2D(h,z,v));
/**----------------------------------------------------------------------
    astronomical potential gradient one its natural discretisation*/
    status = fe_gradient(mesh, z_descriptor, cefmo.state.potential, m, px, py);
    sumx=real(fe_integraleZxUxGZ_2D(h,ustar,px));
    sumy=real(fe_integraleZxUxGZ_2D(h,vstar,py));
/**----------------------------------------------------------------------
    o.5 is for time-averaging, 2.0 is for space-averaging*/
    RoW[0][m] = 0.5 * rhoPSI * gPSI * (sumx+sumy) *2.0;
    sum[0]   += area*RoW[0][m];
/**----------------------------------------------------------------------
    LSA potential gradient one its natural discretisation*/
    status = fe_gradient(mesh, z_descriptor, cefmo.state.LSA, m, px, py);
    sumx=real(fe_integraleZxUxGZ_2D(h,ustar,px));
    sumy=real(fe_integraleZxUxGZ_2D(h,vstar,py));
    RoW[1][m] = 0.5 * rhoPSI * gPSI * (sumx+sumy) *2.0;
    sum[1]   += area*RoW[1][m];
/**----------------------------------------------------------------------
    wave drag rate of work*/
    for(i = 0; i < u_descriptor.nnpe; i++) {
      n=u_descriptor.NIbE[m][i];
      uu=cefmo.state.u[n];
      vv=cefmo.state.v[n];
      taux[i] = cefmo.DragMatrix[n].c[0][0]*uu + cefmo.DragMatrix[n].c[0][1]*vv;
      tauy[i] = cefmo.DragMatrix[n].c[1][0]*uu + cefmo.DragMatrix[n].c[1][1]*vv;
      }
    sumx=real(fe_integraleZxUxU_2D(h,ustar,taux));
    sumy=real(fe_integraleZxUxU_2D(h,vstar,tauy));
    RoW[2][m] = 0.5 * rhoPSI * (sumx+sumy) *2.0;
    sum[2]   += area*RoW[2][m];
/**----------------------------------------------------------------------
    elevation gradient one its natural discretisation*/
    status = fe_gradient(mesh, z_descriptor, cefmo.state.z, m, px, py);
/**----------------------------------------------------------------------
    pressure gradient rate of work*/
    sumx=real(fe_integraleZxUxGZ_2D(h,ustar,px));
    sumy=real(fe_integraleZxUxGZ_2D(h,vstar,py));
    RoW[5][m] = -0.5 * rhoPSI * gPSI * (sumx+sumy) *2.0;
    sum[5]   += area*RoW[5][m];
// /**----------------------------------------------------------------------
//     Linearized bottom friction rate of work */ /// HERE !!!
//     for(i = 0; i < u_descriptor.nnpe; i++) {
//       n=u_descriptor.NIbE[m][i];
//       uu=cefmo.state.u[n];
//       vv=cefmo.state.v[n];
//       Cd[i]=data.Cd[n];
// //       taux[i] = cefmo.FrictionMatrix[n].c[0][0]*uu + cefmo.FrictionMatrix[n].c[1][0]*vv;
// //       tauy[i] = cefmo.FrictionMatrix[n].c[0][1]*uu + cefmo.FrictionMatrix[n].c[1][1]*vv;
//       complex<double> r1, r2, r3, r4;
//       double r=data.rlinear[n];
//       double u0=data.u0[n];
// //      friction_coefficient(cefmo.FrictionMatrix[n], Cd, r, u0, h, CdRatio[n], &r1, &r2, &r3, &r4);
//       friction_coefficient(cefmo.FrictionMatrix[n], (double) data.Cd[n], r, u0, 1.0, 1.0, &r1, &r2, &r3, &r4);
//       Cd[i]=1.0;
//       taux[i] = r1*uu + r2*vv;
//       tauy[i] = r3*uu + r4*vv;
//       }
//     sumx=real(fe_integraleUxUxU_2D(Cd,ustar,taux));
//     sumy=real(fe_integraleUxUxU_2D(Cd,vstar,tauy));
//     RoW[6][m] = -0.5 * rhoPSI * (sumx+sumy) *2.0;
//     sum[6]   += area*RoW[6][m];

/**----------------------------------------------------------------------
    Mean kinetic energy */
    sumx=real(fe_integraleZxUxU_2D(h,ustar,u));
    sumy=real(fe_integraleZxUxU_2D(h,vstar,v));
    MKE[m] = 0.5 * rhoPSI * (sumx+sumy) *2.0;
    sum[10]   += area*MKE[m];
/**----------------------------------------------------------------------
    Mean potential energy */
    sumx=real(fe_integraleZxZ_2D(zstar,z));
    MPE[m] = 0.5 * rhoPSI* gPSI * sumx *2.0;
    sum[11]   += area*MPE[m];
/**----------------------------------------------------------------------
    Q  (quality) factor 
    definition can be found in : Dynamics of the long-period tides
    CARL WUNSCH, DALE B. HAIDVOGEL, MOHAMED ISKANDARANI and R. HUGHES */
    Q[m] = 2*M_PI*(MKE[m]+MPE[m])/(RoW[2][m]+RoW[6][m])/T;
    Q[m] = MIN(fabs(Q[m]),100.);
    sum[12]   += area*Q[m];
    mesh_area+=area;
    }

  sum[12]/=mesh_area;
  
  print_budget(stdout,"whole domain", wave.name, sum);
  
  sprintf(ascii_file, "%s/%s.energy-budget.log",gOutputPath, wave.name);
  if(iteration==1) {
    out=fopen(ascii_file,"w");
    }
  else {
    out=fopen(ascii_file,"a");
    }
    
  print_budget(out,"whole domain", wave.name, sum);
  fclose(out);
  
  comments[0]=new char[128];
  comments[1]=new char[128];

  double *zLGP1=new double[mesh.nvtxs];

  status=projection_LGP1(RoW[0],zLGP1,mesh,LGP0);
  sprintf(ascii_file, "%s/%s.row-pot.%2.2d.LGP1.s2r",gOutputPath, wave.name, iteration);
  sprintf(comments[0], "spectral solver, %s", tugo_version);
  sprintf(comments[1], "astronomic potential rate of work (w/m²)\n");
  status=fe_ascii_saver1(ascii_file, mesh, zLGP1, mask, LGP1, fmt, comments);

  status=projection_LGP1(RoW[1],zLGP1,mesh,LGP0);
  sprintf(ascii_file, "%s/%s.row-lsa.%2.2d.LGP1.s2r",gOutputPath, wave.name, iteration);
  sprintf(comments[0], "spectral solver, %s", tugo_version);
  sprintf(comments[1], "LSA potential rate of work (w/m²)\n");
  status=fe_ascii_saver1(ascii_file, mesh, zLGP1, mask, LGP1, fmt, comments);

  status=projection_LGP1(RoW[5],zLGP1,mesh,LGP0);
  sprintf(ascii_file, "%s/%s.row-prs.%2.2d.LGP1.s2r",gOutputPath, wave.name, iteration);
  sprintf(comments[0], "spectral solver, %s", tugo_version);
  sprintf(comments[1], "pressure forces rate of work (w/m²)\n");
  status=fe_ascii_saver1(ascii_file, mesh, zLGP1, mask, LGP1, fmt, comments);

  status=projection_LGP1(RoW[2],zLGP1,mesh,LGP0);
  sprintf(ascii_file, "%s/%s.row-wd.%2.2d.LGP1.s2r",gOutputPath, wave.name, iteration);
  sprintf(comments[0], "spectral solver, %s", tugo_version);
  sprintf(comments[1], "wave drag rate of work (w/m²)\n");
  status=fe_ascii_saver1(ascii_file, mesh, zLGP1, mask, LGP1, fmt, comments);

  status=projection_LGP1(RoW[6],zLGP1,mesh,LGP0);  
  sprintf(ascii_file, "%s/%s.row-bfr.%2.2d.LGP1.s2r",gOutputPath, wave.name, iteration);
  sprintf(comments[0], "spectral solver, %s", tugo_version);
  sprintf(comments[1], "bottom friction rate of work (w/m²)\n");
  status=fe_ascii_saver1(ascii_file, mesh, zLGP1, mask, LGP1, fmt, comments);

  delete[] zLGP1;
  
  delete[] z; delete[] u; delete[] v; delete[] zstar; delete[] ustar; delete[] vstar;
  delete[] px; delete[] py; delete[] taux; delete[] tauy;
  
  return(0);  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int spectral_EnergyBudget(const char *energy, const char *polygonsfile, spectrum_t spectrum, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,status;
  int verbose=0;
  mesh_t mesh;
  string varname[20];
  poc_var_t var;
  double *UGbufR[20], mask;
  int nvars;
  double RoW[20];
  vector<string> waves;
  vector<plg_t> polygons;
  double *x,*y;
  int *zone;
  
  for(k=0;k<20;k++) UGbufR[k]=0;
  
//   for (k=0;k<spectrum.n;k++) {
//     char *tmp;
//     status=tide_decode_atlasname(path, convention, spectrum.waves[k].name, 0, &tmp);
//     input.push_back(tmp);
//     }
  
  status=plg_load(polygonsfile, PLG_FORMAT_UNKNOWN, polygons);
  
  printf("#################################################################\n");
  printf("load mesh geometry from %s\n",energy);
  status=fe_readgeometry(energy, &mesh);
  if(status) NC_TRAP_ERROR(return,status,1,"fe_readgeometry(\"%s\",) error",energy);
  
  x=new double[mesh.ntriangles];
  y=new double[mesh.ntriangles];
  for(size_t m=0;m<mesh.ntriangles;m++) {
    fe_position(mesh, mesh.triangles[m], &x[m], &y[m],0);
    }
  
  zone=new int[mesh.ntriangles];
  for(size_t m=0;m<mesh.ntriangles;m++) zone[m]=-1;
  
  for(int p=0;p<polygons.size();p++) {
    char *inside=plg_TestInterior(x,y,mesh.ntriangles,polygons,verbose,debug);
    for(size_t m=0;m<mesh.ntriangles;m++) {
      if(inside[m]==1) zone[m]=p;
      }
    delete[] inside;
    }
  
  varname[0]="pressureRoW";
  varname[1]="BFrictionRoW";
  varname[2]="WaveDragRoW";
  varname[3]="AstroP_RoW";
  varname[4]="LSA_RoW";

  nvars=5;

//   for(int v=0;v<waves.size();v++) {
  for(int v=0;v<spectrum.n;v++) {
    for(k=0;k<20;k++) RoW[k]=0;
    for(k=0;k<nvars;k++) {
      status=poc_inq_var(energy, varname[k], &var);
      status=poc_decode_mask(var, 0, 0, &mask);
      if(UGbufR[k]==0) UGbufR[k]=new double[mesh.ntriangles];
      status=poc_get_UG3D(energy, -1, varname[k].c_str(), UGbufR[k]);
      RoW[k]=0;
      for(size_t m=0;m<mesh.ntriangles;m++) {
//         double area=mesh.triangles[m].TrueArea;
        double area=mesh.triangles[m].Area*mesh.triangles[m].cosinus;
//         RoW[k]+=2.0 * UGbufR[k][m]*area; //HERE!!!
        RoW[k]+=UGbufR[k][m]*area;
        } 
      }
    printf("entire mesh domain budget\n");
    printf("pressure forces rate of work (Gw): %f\n", RoW[0] / 1.e+09);
    printf("bottom friction rate of work (Gw): %f\n", RoW[1] / 1.e+09);
    printf("wave drag       rate of work (Gw): %f\n", RoW[2] / 1.e+09);
    printf("astronom. pot.  rate of work (Gw): %f\n", RoW[3] / 1.e+09);
    printf("LSA pot.        rate of work (Gw): %f\n", RoW[4] / 1.e+09);
    printf("Energy balance closure       (Gw): %f\n", (RoW[0]+RoW[1]+RoW[2]+RoW[3]+RoW[4]) / 1.e+09);
    for(int p=0;p<polygons.size();p++) {
      for(k=0;k<20;k++) RoW[k]=0;
      for(k=0;k<nvars;k++) {
        RoW[k]=0;
        for(size_t m=0;m<mesh.ntriangles;m++) {
         if(zone[m]!=p) continue;
          double area=mesh.triangles[m].TrueArea;
//           RoW[k]+=2.0 * UGbufR[k][m]*area; //HERE!!!
          RoW[k]+=UGbufR[k][m]*area;
          } 
        } 
      printf("region %d domain budget\n",p);
      printf("pressure forces rate of work (Gw): %f\n", RoW[0] / 1.e+09);
      printf("bottom friction rate of work (Gw): %f\n", RoW[1] / 1.e+09);
      printf("wave drag       rate of work (Gw): %f\n", RoW[2] / 1.e+09);
      printf("astronom. pot.  rate of work (Gw): %f\n", RoW[3] / 1.e+09);
      printf("LSA pot.        rate of work (Gw): %f\n", RoW[4] / 1.e+09);
      printf("Energy balance closure       (Gw): %f\n", (RoW[0]+RoW[1]+RoW[2]+RoW[3]+RoW[4]) / 1.e+09);
      }
    }
//     printf("wave %s\n", wave.name);
//     printf("astronom. pot.  rate of work (Gw): %f\n", sum[0] / 1.e+09);
//     printf("secondary pot.  rate of work (Gw): %f\n", sum[1] / 1.e+09);
//     printf("wave drag       rate of work (Gw): %f\n", sum[2] / 1.e+09);
//     printf("pressure forces rate of work (Gw): %f\n", sum[3] / 1.e+09);
//     printf("pressure forces rate of work (Gw): %f\n", sum[5] / 1.e+09);
//     printf("energy advection rate        (Gw): %f\n", sum[4] / 1.e+09);
//     printf("bottom friction rate of work (Gw): %f\n", sum[6] / 1.e+09);
//     printf("bottom friction rate of work (Gw): %f\n", sum[7] / 1.e+09);
//     printf("horizontal dif. rate of work (Gw): %f\n", sum[8] / 1.e+09);
//     printf("Tidal forcing rate of work   (Gw): %f\n", (sum[0]+sum[1]+sum[3]-sum[4]) / 1.e+09);
//     printf("Dissipation rate of work     (Gw): %f\n", (sum[2]+sum[7]+sum[8]) / 1.e+09);
//     printf("Energy balance closure       (Gw): %f\n", (sum[0]+sum[1]+sum[2]+sum[3]-sum[4]+sum[7]+sum[8]) / 1.e+09);
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int spectral_EnergyComputation(const char *meshfile, const char *discretisation, const char *path, const char *output_path, spectrum_t & spectrum, spectrum_t & dominants)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  
  cefmo_t cefmo;
  mesh_t mesh;
  tidal_wave wave;
  parameter_t data;
  int iteration;
  double *RoW[10], *MKE, *MPE, *Q;
  int nelements=0, element_discretisation=-1;
  int zdim,udim,paire;
  int z_discretisation, u_discretisation;
  zmatrix2x2_t *FrictionMatrix[2][2];
  char filename[1024];
  int frame=-1;
  tide2D_t state;
  
  if(discretisation==0) {
    printf("computational paire not specified, abort...\n");
    goto error;
    }
  else if(strcmp(discretisation,"LGP0xLGP1")==0) {
    paire=LGP0xLGP1;
    }
  else if(strcmp(discretisation,"DGP1xLGP2")==0) {
    paire=DGP1xLGP2;
    }
  else if(strcmp(discretisation,"DNP1xLGP2")==0) {
    paire=DNP1xLGP2;
    }
  else if(strcmp(discretisation,"Q1xQ0")==0) {
    paire=CQP1xCQP0;
    }
  else {
    goto error;
    }
  
  sprintf(filename,"%s/%s.spectral.nc",path, spectrum.waves[0].name);
  printf("#################################################################\n");
  printf("load mesh geometry from %s\n",filename);
  status=fe_readgeometry(filename, &mesh);
  if(status!=0) NC_TRAP_ERROR(return,status,1,"fe_readgeometry(\"%s\",) error",filename);
  
  
  switch(mesh.nquadrangles) {
    case 0:
      if(mesh.LGP0descriptor.nnodes==0) {
        status= discretisation_init(&mesh, LGP0, SEQUENTIAL_COMPUTING);
        }
      nelements=mesh.ntriangles;
      element_discretisation=LGP0;
      break;
    default:
      if(mesh.CQP0descriptor.nnodes==0) {
        status= discretisation_init(&mesh, CQP0, SEQUENTIAL_COMPUTING);
        }
      break;
      nelements=mesh.nquadrangles;
      element_discretisation=CQP0;
    }     
  
  for(int k = 0; k < 10; k++) {
    RoW[k] = new double[nelements]; 
    }
  MKE = new double[nelements]; 
  MPE = new double[nelements]; 
  Q   = new double[nelements]; 
  
/*------------------------------------------------------------------------------
  botttom friction is an issue, as it has to handled consistently with the
  original T-UGOm computation. Let it be handled prior to any thing else
------------------------------------------------------------------------------*/

  cefmo.paire=paire;
  for(int k = 0; k < dominants.n; k++) {
    cefmo.dominants[k]=dominants.waves[k];
    }

  status=paire_discretisation_id(paire, &z_discretisation,&u_discretisation); 
  
  status=discretisation_init(&mesh,z_discretisation,0);
  status=discretisation_init(&mesh,u_discretisation,0);

  status=paire_dimension(mesh,cefmo.paire,&zdim,&udim);
 
//   status=bottom_friction();
  
  state.z=new complex<double> [zdim];
  state.u=new complex<double> [udim];
  state.v=new complex<double> [udim];
  
  for(int k = 0; k < dominants.n; k++) {
    sprintf(filename,"%s/%s.spectral.nc",path, cefmo.dominants[0].name);
    printf("#################################################################\n");
    printf("load dominant wave state vector %s\n",filename);
    status=SpArchive_Get2D(filename, mesh, frame, cefmo, state);
    if(status!=0) NC_TRAP_ERROR(return,status,1,"SpArchive_Get2D(\"%s\",) error",filename);
    FrictionMatrix[k][0] = spectral_friction01(mesh, state.u, state.v, udim);
    FrictionMatrix[k][1] = spectral_friction02(mesh, state.u, state.v, udim);
    }
  
  
  state.h_z=new double[zdim];
  state.potential=new complex<double> [zdim];
  state.LSA=new complex<double> [zdim];

  cefmo.state=state;

  for(int k = 0; k < spectrum.n; k++) {
    sprintf(filename,"%s/%s.spectral.nc",path, spectrum.waves[k].name);
    status=SpArchive_GetFull2D(filename, mesh, frame, cefmo, cefmo.state);
    status=spectral_energy_T(mesh, spectrum.waves[k], cefmo, data, iteration, RoW, MKE, MPE, Q);
    }

  return (status);
  
  
  
error:
  return (-1);
}
