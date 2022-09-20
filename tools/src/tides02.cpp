
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

\brief tidal energy flux definitions UNUSED !!!
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <string.h>
#include <complex>
#include <cmath>


#include "tools-structures.h"

#include "geo.h"
#include "fe.h"

using namespace std;

extern double Energyflux, EpFlux, EkFlux, EppFlux, EkpFlux;
extern double MassP, Mass, MassFlux;
extern double wp, Wp, DivE, DivHu1, DivHu2, DivHu3;
extern double DailyBudget;
extern double Ep_previous, Ek_previous;
extern double Ep, Ek, Wf, Wv, Wtp, Wtl, Wws;
extern double Wpa, Cr, energy, work, Epp;
extern double Ek1, Ek2;

/*float g=9.81, ro=1027.0;*/

double P_g=9.81,rho_water=1027.0;

extern matrix2x2_t  *gDragMatrix;
extern zmatrix2x2_t *gFricMatrix;
/*----------------------------------------------------------------------
 
 Energy averaged over the water column
 
 Ep = 1/2 ro x g x H x dH  potential energy
 Ek = 1/2 ro x H x u**2    kinetic energy
 Ei = 1/2 ro x g x H x h   intrisinc energy
 
----------------------------------------------------------------------*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float integrale_time_2(fcomplex z1, fcomplex z2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float d;
/*************************************************************************
  Time average over 1 period
**************************************************************************/

/* FHL 04/06/2006 COMPLEX
  d = 0.5 * (z1.r * z2.r + z1.i * z2.i);
*/
  d = 0.5 * (real(z1) * real(z2) + imag(z1) * imag(z2));
  return (d);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float integrale_time_4(fcomplex z1, fcomplex z2, fcomplex z3, fcomplex z4)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float d;
/*
d=z1*CONJG(z2)*CONJG(z3)*z4 +z1*CONJG(z2)*z3*CONJG(z4)
      +z1*z2*CONJG(z3)*CONJG(z4) +CONJG(z1)*CONJG(z2)*z3*z4
      +CONJG(z1)*z2*CONJG(z3)*z4 +CONJG(z1)*z2*z3*CONJG(z4);

d=d/16.0; */

  return (d);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void spectralenergy02(mesh_t mesh, spectrum_t WaveList, parameter_t *gdata2D, action_t **P1_action)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, m, w;
  float u, v, A, U, gPSI, rhoPSI, C;
  float a1, a2, G1, G2, z, zr, zi;
  double row[10][3];
  double beta[3];
  double kappa1, kappa2;
  double tmp[10], sum[10], budget;
  double h, H, dH, dHdT;
  double delta1, delta2, c1, c2;
  double omega, factor = 1.e+06 * 3600.0;
  int n;
  tidal_wave wave;
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int gauss_n = 3;

/*----------------------------------------------------------------------
  Compute the energy budget over the CLOSE domain*/

  for(n = 0; n < mesh.nvtxs; n++) {
    P1_action[1][n].tpG = new float[WaveList.n];
    P1_action[2][n].tpG = new float[WaveList.n];
    }

//WARNING fudged P1_state to fit function prototype
  parameter_t P1_state;
/* *----------------------------------------------------------------------------
  to be replaced*/
//  tidal_potential_constants(mesh, P1_state, P1_action);

/* *----------------------------------------------------------------------------
  to be replaced*/
  gPSI = P_g;
  rhoPSI = rho_water;

  gauss_init(gauss_n, gauss_x, gauss_y, gauss_w);

  for(w = 0; w < WaveList.n; w++) {
    wave = WaveList.wave[w];
    omega = wave.omega * dph2rps;
/*-------------------------------------------------------
     u=gdata2D[0].Utide[n].a[k]*cos(gdata2D[0].Utide[n].G[k]);
     v=gdata2D[0].Vtide[n].a[k]*cos(gdata2D[0].Vtide[n].G[k]);
--------------------------------------------------------*/

    for(j = 0; j < 2; j++)
      sum[j] = 0.0;
    for(m = 0; m < mesh.ntriangles; m++) {
/*----------------------------------------------------------------------
      convert A (element area) in m**2*/
      A = mesh.triangles[m].Area / (1.0e+4);
      for(j = 0; j < 2; j++)
        tmp[j] = 0;
      for(i = 0; i < 3; i++) {
        n = mesh.triangles[m].vertex[i];
/*----------------------------------------------------------------------
        convert dpth in m*/
        h = gdata2D[0].h[n];
        C = gdata2D[0].C[n];
/*----------------------------------------------------------------------
        convert amplitude in m*/
        a1 = gdata2D[0].Htide[n].a[w];
        G1 = gdata2D[0].Htide[n].G[w];
/*----------------------------------------------------------------------
        astronomical potential rate of work*/
        a2 = P1_action[1][n].tpA[w];
        G2 = P1_action[1][n].tpA[w];
        row[0][i] =
              -0.5 * h * rhoPSI * gPSI * omega * a1 * a2 * sin(G2 - G1) * C;
/*----------------------------------------------------------------------
        secundary potential rate of work*/
/* FHL 04/06/2006 COMPLEX
        zr = P1_action[2][n].tlc[w].r;
        zi = P1_action[2][n].tlc[w].i;
        a2 = sqrt(zr * zr + zi * zi);
        G2 = -atan2(zi, zr);
*/

       a2 =  abs(P1_action[2][n].tlc[w]);
       G2 = -arg(P1_action[2][n].tlc[w]);
       row[1][i] =
              -0.5 * h * rhoPSI * gPSI * omega * a1 * a2 * sin(G2 - G1) * C;
      }
/*----------------------------------------------------------------------
      spatial integration*/
      for(k = 0; k < gauss_n; k++) {
        fe_LGP1base(gauss_x[k], gauss_y[k], beta);
        for(j = 0; j < 2; j++) {
          z = 0.;
          for(i = 0; i < 3; i++) {
            z += row[j][i] * beta[i];
          }
          tmp[j] += gauss_w[k] * z;
        }
      }
      for(j = 0; j < 2; j++)
        sum[j] += A * tmp[j];
    }   /* end of element loop */
    printf("wave %s\n", wave.name);
    printf("astronom. pot. rate of work (Gw): %f\n", sum[0] / 1.e+09);
    printf("secundary pot. rate of work (Gw): %f\n", sum[1] / 1.e+09);
  }     /* end of wave loop */

  for(n = 0; n < mesh.nvtxs; n++) {
    delete[] P1_action[1][n].tpA;
    delete[] P1_action[2][n].tpG;
  }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double flux_integrale(mesh_t mesh,parameter_t data, double **buffer)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
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

  for(n = 0; n < mesh.nedges; n++) {
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

void spectral_flux(mesh_t mesh,spectrum_t spectrum,double duration,state2d_t state, parameter_t data)

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
/*
  duration=3600.*24*15;
  dt=2.*M_PI/(spectrum.wave[0].omega*d2r/3600.)/24.0;
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

  status=harmonic_init01(nframe, start, time, spectrum, &harmonic,0);

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
    sum= flux_integrale(mesh,data, &(buffer[w][0]));
    printf("%s total potential energy flux : %12.4f gW\n",spectrum.wave[w].name,sum/1.e+09);
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
    sum= flux_integrale(mesh,data, &(buffer[w][0]));
    printf("%s partial potential energy flux : %12.4f gW\n",spectrum.wave[w].name,sum/1.e+09);
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


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void spectral_flux01(mesh_t mesh, int w, z_spec *zspec, int nobc)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
  Compute the spectral energy fluxes at the open boundaries*/
{
  double Nx, Ny, L, gPSI, rhoPSI;
  double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0;
  double h, u, v, F;
  double a1, a2, G1, G2, a3, G3;
  int i,k,n;

  gPSI = P_g;
  rhoPSI = rho_water;

  for(i = 0; i < nobc; i++) {
    n = zspec[i].node;
/*----------------------------------------------------------------------
    normal to bpundary. Outward?*/
    Nx = (double) zspec[i].cosine;
    Ny = (double) zspec[i].sine;
    L =  (double) zspec[i].size / 2.0;

    a1 = gdata2D[0].Htide[n].a[w];
    G1 = gdata2D[0].Htide[n].G[w];
    a2 = gdata2D[0].Utide[n].a[w];
    G2 = gdata2D[0].Utide[n].G[w];
    a3 = gdata2D[0].Vtide[n].a[w];
    G3 = gdata2D[0].Vtide[n].G[w];

    h = gdata2D[0].h[n];

    u = 0.5 * h * rhoPSI * gPSI *a1 * a2 * cos(G2 - G1);
    v = 0.5 * h * rhoPSI * gPSI *a1 * a3 * cos(G3 - G1);

    F = L * (u * Nx + v * Ny);

/*----------------------------------------------------------------------
    Potential energy flux (energy budget form)*/
    sum1 += F;
    }
  printf("Potential energy flux : %8.1f gW\n",sum1/1.0e+09);
}


 /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void TotalEnergy_RH5(FILE * out, double t, mesh_t mesh)
 /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n;
  float u, v, A, U, gPSI, rhoPSI, C;
  float row[10][3], beta[3], z;
  double tmp[10], sum1, sum2;
  double kappa1, kappa2;
  double budget;
  double h, H, dH, dHdT, T;
  double delta1, delta2, c1, c2, omega;
  double Mystery = 0.0, factor = 1.e+06 * 3600.0;
  double Aire_Totale;
  fcomplex cplx1, cplx2, cplx3;

/*----------------------------------------------------------------------
  compute the total energy of the model
  for the RH5 frequency (T=5 days) ...
----------------------------------------------------------------------*/

  gPSI = P_g;
  rhoPSI = rho_water;

  m = 0;
/*  while(strcmp(AnalysisList.wave[m].name, "RH5") != 0)*/
  while(strcmp(AnalysisList.wave[m].name, "M2") != 0)
    m++;

  printf("#total energy calculation for %s\n", AnalysisList.wave[m].name);
  T = 3600. * 360. / AnalysisList.wave[m].omega;        /* periode de l'onde
                                                                 en secondes */
  printf("#periode de %s = %lf \n", T);

  sum1 = 0.0;
  sum2 = 0.0;
  Aire_Totale = 0.0;

  for(k = 0; k < mesh.ntriangles; k++) {
/*----------------------------------------------------------------------
    convert A (element area) in m**2
----------------------------------------------------------------------*/
    A = mesh.triangles[k].Area / (1.0e+4);

    for(j = 0; j < 2; j++)
      tmp[j] = 0;

    for(i = 0; i < 3; i++) {
      n = mesh.triangles[k].vertex[i];
      H = gdata2D[0].Htide[n].a[m];        /* instantaneous harmonic elevation */
      h = gdata2D[0].h[n];    /* mean depth */
      dH = H - h;
      C = gdata2D[0].C[n];

      u = gdata2D[0].Utide[n].a[m];
      v = gdata2D[0].Vtide[n].a[m];

      c1 = rhoPSI * gPSI;
      c2 = rhoPSI * h;

/*----------------------------------------------------------------------
  Potential energy
----------------------------------------------------------------------*/
/* FHL 04/06/2006 COMPLEX
      cplx1.r = H * cos(gdata2D[0].Htide[n].G[m]);
      cplx1.i = H * sin(gdata2D[0].Htide[n].G[m]);
*/
      cplx1 =  fcomplex(H *cos(gdata2D[0].Htide[n].G[m]),H *sin(gdata2D[0].Htide[n].G[m]));
      tmp[0] += 0.5 * c1 * integrale_time_2(cplx1, cplx1);

/*----------------------------------------------------------------------
  Kinetic energy
----------------------------------------------------------------------*/
/* FHL 04/06/2006 COMPLEX
      cplx2.r = u * cos(gdata2D[0].Utide[n].G[m]);
      cplx2.i = u * sin(gdata2D[0].Utide[n].G[m]);
      cplx3.r = v * cos(gdata2D[0].Vtide[n].G[m]);
      cplx3.i = v * sin(gdata2D[0].Vtide[n].G[m]);
*/
      cplx2 = u * fcomplex(cos(gdata2D[0].Utide[n].G[m]),sin(gdata2D[0].Utide[n].G[m]));
      cplx3 = v * fcomplex(cos(gdata2D[0].Vtide[n].G[m]),sin(gdata2D[0].Vtide[n].G[m]));
      tmp[1] +=
            0.5 * c2 * (integrale_time_2(cplx2, cplx2) +
                        integrale_time_2(cplx3, cplx3));
    }

    tmp[0] *= C * A / 3.0;
    tmp[1] *= C * A / 3.0;

    sum1 += tmp[0];
    sum2 += tmp[1];
    Aire_Totale += C * A;
  }

  fprintf(out, "#RH5  potential energy : %lf\t %f\n", t, sum1 / (T * 1.e+09)); /* en Gw  */
  fprintf(out, "#RH5  kinetic   energy : %lf\t %f\n", t, sum2 / (T * 1.e+09)); /* en Gw  */
  fprintf(out, "#RH5  total     energy : %lf\t %f\n", t, (sum1 + sum2) / (T * 1.e+09));    /* en Gw  */
  fprintf(out, "#Aire Totale           : %lf km*km \t", Aire_Totale / 1.e+06);  /* en km */
}

 /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void Qfactor(FILE * out, double t, mesh_t mesh)
 /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, l, k, m, mm;
  float u, v, A, U, UM2, gPSI, rhoPSI, C;
  float row[10][3], beta[3], z;
  double tmp[10], sum[10];
  double kappa, T;
  double budget, ETot, DTot, qfactor, Td1, Td2;
  double h, H, dH, dHdT;
  double delta1, delta2, c1, c2;
  double Mystery = 0.0, factor = 1.e+06 * 3600.0;
  int n;
  fcomplex cplx1, cplx2, cplx3, cplx4, cplx5;
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------
  compute the Quality factor of the model
  through the total energy over the domain and
  the damping rate of the model
  energie integrated on the period T=1 hour or 5 days -> a confirmer ?...
----------------------------------------------------------------------*/
  gPSI = P_g;
  rhoPSI = rho_water;
  kappa = 2 * pi / (1500.);

  m = 0;
  while(strcmp(AnalysisList.wave[m].name, "RH5") != 0)
    m++;

  T = 3600. * 360. / AnalysisList.wave[m].omega;
                                         /* periode de l'onde en secondes */
  printf("#periode de %s = %lf \n", T);

  mm = 0;
  while(strcmp(AnalysisList.wave[mm].name, "M2") != 0)
    mm++;

  for(j = 0; j < 4; j++)
    sum[j] = 0.0;

  for(k = 0; k < mesh.ntriangles; k++) {
/*----------------------------------------------------------------------
    convert A (element area) in m**2
----------------------------------------------------------------------*/
    A = mesh.triangles[k].Area / (1.0e+4);

    for(j = 0; j < 4; j++)
      tmp[j] = 0;

    for(i = 0; i < 3; i++) {
      n = mesh.triangles[k].vertex[i];
      H = gdata2D[0].Htide[n].a[m];
      h = gdata2D[0].h[n];
      dH = H - h;
      C = gdata2D[0].C[n];

      u = gdata2D[0].Utide[n].a[mm];
      v = gdata2D[0].Vtide[n].a[mm];
      UM2 = sqrt(u * u + v * v + 1.e-30);

      u = gdata2D[0].Utide[n].a[m];
      v = gdata2D[0].Vtide[n].a[m];
      U = sqrt(u * u + v * v + 1.e-30);

      c1 = rhoPSI * gPSI;
      c2 = rhoPSI * h;

/*----------------------------------------------------------------------
  Potential energy
----------------------------------------------------------------------*/
/* FHL 04/06/2006 COMPLEX
      cplx1.r = H * cos(gdata2D[0].Htide[n].G[m]);
      cplx1.i = H * sin(gdata2D[0].Htide[n].G[m]);
*/
      cplx1 =  fcomplex(H *cos(gdata2D[0].Htide[n].G[m]),H *sin(gdata2D[0].Htide[n].G[m]));
      tmp[0] = 0.5 * c1 * integrale_time_2(cplx1, cplx1);

/*----------------------------------------------------------------------
  Kinetic energy
----------------------------------------------------------------------*/
/* FHL 04/06/2006 COMPLEX
      cplx2.r = u * cos(gdata2D[0].Utide[n].G[m]);
      cplx2.i = u * sin(gdata2D[0].Utide[n].G[m]);
      cplx3.r = v * cos(gdata2D[0].Vtide[n].G[m]);
      cplx3.i = v * sin(gdata2D[0].Vtide[n].G[m]);
*/
      cplx2 = u * fcomplex(cos(gdata2D[0].Utide[n].G[m]),sin(gdata2D[0].Utide[n].G[m]));
      cplx3 = v * fcomplex(cos(gdata2D[0].Vtide[n].G[m]),sin(gdata2D[0].Vtide[n].G[m]));
      tmp[1] =
            0.5 * c2 * (integrale_time_2(cplx2, cplx2) +
                        integrale_time_2(cplx3, cplx3));

/*----------------------------------------------------------------------
  damping rate for bottom friction dissipation (friction work)
  Considering the Depth dependent background velocity effect: U0 at 100 m
----------------------------------------------------------------------*/
/*      tmp[2] = -rhoPSI*Cd*(U+U0*sqrt(100./H))*U*U*CheckInterval/H;*/
      tmp[2] = -1.*rhoPSI * gdata2D[0].Cd[n] * (UM2 + U0 * sqrt(1. / H)) *
              (integrale_time_2(cplx2, cplx2) +integrale_time_2(cplx3, cplx3)) / H;

/*----------------------------------------------------------------------
  damping rate for internal waves genaration process
----------------------------------------------------------------------*/
      tmp[3] = (-1. / kappa) * rhoPSI * gstate2D[0][1].N[n] *
               (gdata2D[0].dhdx[n] * gdata2D[0].dhdx[n] *
               integrale_time_2(cplx2, cplx2) +
             gdata2D[0].dhdy[n] * gdata2D[0].dhdy[n] *
             integrale_time_2(cplx3, cplx3) +
             2. * gdata2D[0].dhdx[n] * gdata2D[0].dhdy[n] *
             integrale_time_2(cplx2, cplx3));

/*----------------------------------------------------------------------
  damping rate for internal viscous dissipation (viscosity work)
----------------------------------------------------------------------*/
      cplx4 = fcomplex(0.,0.);
      cplx5 = fcomplex(0.,0.);
      tmp[4] = c2 * (integrale_time_2(cplx2, cplx4) + integrale_time_2(cplx3, cplx5));
    }

/*----------------------------------------------------------------------
      spatial integration
----------------------------------------------------------------------*/
    tmp[0] *= C * A / 3.0;
    tmp[1] *= C * A / 3.0;
    tmp[2] *= C * A / 3.0;
    tmp[3] *= C * A / 3.0;
    tmp[4] *= C * A / 3.0;

    sum[0] += tmp[0];
    sum[1] += tmp[1];
    sum[2] += tmp[2];
    sum[3] += tmp[3];
    sum[4] += tmp[3];
  }
  ETot = sum[0] + sum[1];
  DTot = sum[2] + sum[3] + sum[4];

/*-------------------------------------------------------------------------
   T=5 jours pour RH5
   qfactor=pi2*ETot/(DTot*T); -> formule donn� par Wunsch 97, Jackson 75 ...
--------------------------------------------------------------------------*/
  qfactor = pi2 * ETot / (DTot * T);
  Td1 = qfactor / (reference_spectrum.wave[m].omega * d2r / 3600.);
  Td2 = ETot / (DTot / (5. * 24. * 3600.));
  fprintf(out, "Qfactor= %lf;\t Dissipation Time Scale= %lf, %lf\n", qfactor,
          Td1, Td2);
}
