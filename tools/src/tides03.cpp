
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

\brief tidal momentum definitions
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>

#include "tools-structures.h"
#include "tides.h"

extern int fft (double *x, int n);
matrix2x2_t  *DragMatrix;
zmatrix2x2_t *FricMatrix;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void accelere(double rn0,double rn1,double *rn2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*****************************************************************
 Cette routine accelere la convergence des calculs de RDAMP
**************************************************************** */
  double a1,a2, quot;

  a1=rn1-rn0;
  a2=*rn2-rn1;
  quot=a2-a1;
  if(fabs(quot)>1.e-10) *rn2=rn0-a1*a1/quot;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void spectral_friction01(mesh_t mesh, parameter_t data, int nndes,int wave)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/****************************************************************
 Cette routine traduit le processus iteratif, decrit dans le
 #TidesHarmonicConstituentsModeling_REF, et mis en oeuvre pour le
 calcul de l'onde dominante. Le calcul de l'acceleration est fait
 sur les coefficients de frottement ro et roprim
*****************************************************************/
  double c=1.0,R1,R2;
  double a,G;
  int n;
  dcomplex u,v,j=dcomplex(0.,1.);

  for(n=0;n<mesh.nvtxs;n++) {
    a=data.Utide[n].a[wave];
    G=data.Utide[n].G[wave];
    u=a * dcomplex(cos(G),sin(-G));
    a=data.Vtide[n].a[wave];
    G=data.Vtide[n].G[wave];
    v=a * dcomplex(cos(G),sin(-G));
    frodom(u,v,&R1,&R2);
    R1+=0.05;

/*----------------------------------------------------------------------
    c[j][i]=matrix's coefficent row i, column j

    r0 r1 u
    r2 r3 v
----------------------------------------------------------------------*/
    FricMatrix[n].c[0][0]=R1;
    FricMatrix[n].c[1][0]= j*R2;
    FricMatrix[n].c[0][1]=-j*R2;
    FricMatrix[n].c[1][1]=R1;
    }
/*
accelere(solve->buffer.pg[1].tri[lpg],solve->buffer.pg[3].tri[lpg],&solve->rdamp.pg[1].tri[lpg]);
accelere(solve->buffer.pg[2].tri[lpg],solve->buffer.pg[4].tri[lpg],&solve->rdamp.pg[2].tri[lpg]);
*/
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void spectral_friction02(mesh_t mesh, parameter_t data, int nndes, int wave)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/****************************************************************
 Cette routine traduit le processus iteratif, decrit dans le
 #TidesHarmonicConstituentsModeling_REF, et mis en oeuvre pour le
 calcul de l'onde secondaire.
*****************************************************************/
  double c=1.0,R1,R2,R3;
  double a,G;
  int n;
  dcomplex u,v;

  for(n=0;n<mesh.nvtxs;n++) {
    a=data.Utide[n].a[wave];
    G=data.Utide[n].G[wave];
    u=a * dcomplex(cos(G),sin(-G));
    a=data.Vtide[n].a[wave];
    G=data.Vtide[n].G[wave];
    v=a * dcomplex(cos(G),sin(-G));
    frosec(u,v,&R1,&R2,&R3);
    R1*=c;
    R2*=c;
    R3*=c;
/*
  r0(k)= (rdamp(3,l)+rdamp(4,l))*c/h
  r1(k)= rdamp(5,l)*c/h
  r2(k)= r1(k)
  r3(k)= (rdamp(3,l)-rdamp(4,l))*c/h
*/
    FricMatrix[n].c[0][0]=(R1+R2);
    FricMatrix[n].c[1][0]=R3;
    FricMatrix[n].c[0][1]=R3;
    FricMatrix[n].c[1][1]=(R1-R2);
    }
}


//#if DEFINED_global_astro_angles
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void spectral_friction03(spectrum_t spectrum,double duration, parameter_t data, int nndes, int target, double **buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double c,R1,R2,R3;
  double a1,G1,a2,G2,a3,G3;
  double *residuals,*taux,*tauy,*KE,tauxm,tauym;
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
 /*-----------------------------------------------------------------------*/
  astro_angles_t astro_angles;
  init_argument(&astro_angles,start);

  time = new double[nframe];
  u    = new double[nframe];
  v    = new double[nframe];
  H    = new double[nframe];
//  KE   = new double[nframe];
  taux = new double[nframe];
  tauy = new double[nframe];
  residuals = new double[nframe];

  ax= new double [spectrum.n];
  Gx= new double [spectrum.n];
  ay= new double [spectrum.n];
  Gy= new double [spectrum.n];

  for(k=0;k<nframe;k++) time[k]=(double)k*dt;

/* *----------------------------------------------------------------------------
  to be replaced*/
  status=harmonic_init01(nframe, start, time, spectrum, &harmonic,0, astro_angles);

  for(n=0;n<nndes;n++) {
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
      c=sqrt(u[m]*u[m]+v[m]*v[m]) + data.u0[n];
/*test change*/
//      taux[m]=c*u[m]/H[m];
//      tauy[m]=c*v[m]/H[m];
      taux[m]=c*u[m];
      tauy[m]=c*v[m];
      }
//    status=harmonic_compute(u,residuals, nframe, harmonic, spectrum, maxstep,ax,Gx);
//    status=harmonic_compute(v,residuals, nframe, harmonic, spectrum, maxstep,ay,Gy);
/* *----------------------------------------------------------------------------
  to be replaced*/
    status=harmonic_compute(taux,residuals, nframe, harmonic, spectrum, maxstep,ax,Gx);
/* *----------------------------------------------------------------------------
  to be replaced*/
    status=harmonic_compute(tauy,residuals, nframe, harmonic, spectrum, maxstep,ay,Gy);
//    if(data[n].Ua[0] > .2) status=fft(taux, nframe);
    for(w = 0; w < spectrum.n; w++) {
      a1=data.Utide[n].a[w];
      G1=data.Utide[n].G[w];
      a2=data.Vtide[n].a[w];
      G2=data.Vtide[n].G[w];
      buffer[w][n]=0.5*(ax[w]*a1*cos(Gx[w]-G1)+ay[w]*a2*cos(Gy[w]-G2));
      }

/*----------------------------------------------------------------------
    empirical verification*/
/*
    for(w = 0; w < spectrum.n; w++) {
      a1=data.Utide[n].a[w];
      G1=data.Utide[n].G[w];
      a2=data.Vtide[n].a[w];
      G2=data.Vtide[n].G[w];
      cs[0]=a1*cos(G1);
      cs[1]=a2*cos(G2);
      sn[0]=a1*sin(G1);
      sn[1]=a2*sin(G2);
      tauxm=0.;
      tauym=0.;
      for(m=0; m<nframe; m++) {
        u[m] = cs[0] * harmonic.cs[w][m] + sn[0] * harmonic.sn[w][m];
        v[m] = cs[1] * harmonic.cs[w][m] + sn[1] * harmonic.sn[w][m];
        tauxm+=u[m]*taux[m];
        tauym+=v[m]*tauy[m];
        }
      tauxm/=nframe;
      tauym/=nframe;
//      buffer[w][n]=0.5*(tauxm+tauym);
      }
*/
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
//#endif

double SMAGR=0.28;
double minKh=1.0;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void smagorinsky_P1NQ(state2d_t P1_state, parameter_t P1_data,mesh_t mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n;
  int n1, n2, n3;
  double D1, D2, AH, mean, rms, A, factor;
  double ux, uy, vx, vy, *tmp;

/*----------------------------------------------------------------------
  FACE VALUE (QUODDY 4 ???)
  double SMAGR  =   0.28;

  d/dx= DQ/2A (m^-1)
  d/dy=-DP/2A

  D1 and D2 are s^-1
  AH and Asmag are m^2/s
----------------------------------------------------------------------*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Notes:

  derivation-related coefficient. Definition:

  1/R ds/dlon = dsdx / 2A
  1/R ds/dlat = dsdy / 2A

  for P1 variable:

  1/R d beta /d lon =  DQ[j]/2A
  1/R d beta /d lat = -DP[j]/2A

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */


  double damp;

  tmp=new double[mesh.ntriangles];

  for(n = 0; n < mesh.nvtxs; n++) {
    P1_state.Ah[n] = 0.0;
    P1_state.ux[n] = 0.0;
    P1_state.uy[n] = 0.0;
    P1_state.vx[n] = 0.0;
    P1_state.vy[n] = 0.0;
    P1_state.lmts[n] = 0.0;
    }

#pragma omp parallel
#pragma shared (elt,P1_state,P1_data,nelt,SMAGR,smb1)
#pragma local (l,n1,n2,n3,AH)
#pragma local (D1,D2,A)

/*BRACKET FOR PARALLEL COMPUTATION*/
  {
#pragma pfor iterate (l=1;nelt;1)
/*----------------------------------------------------------------------
    Use of variational approach in the nodal quadrature frame*/
    for(l = 0; l < mesh.ntriangles; l++) {
      A = mesh.triangles[l].Area;
      n1 = mesh.triangles[l].vertex[0];
      n2 = mesh.triangles[l].vertex[1];
      n3 = mesh.triangles[l].vertex[2];
      ux =  (P1_state.u[n1] * mesh.triangles[l].DQ[0] / P1_data.C[n1]
           + P1_state.u[n2] * mesh.triangles[l].DQ[1] / P1_data.C[n2]
           + P1_state.u[n3] * mesh.triangles[l].DQ[2] / P1_data.C[n3]) / (2. * A);
      uy = -(P1_state.u[n1] * mesh.triangles[l].DP[0]
           + P1_state.u[n2] * mesh.triangles[l].DP[1]
           + P1_state.u[n3] * mesh.triangles[l].DP[2]) / (2. * A);
      vx =  (P1_state.v[n1] * mesh.triangles[l].DQ[0] / P1_data.C[n1]
           + P1_state.v[n2] * mesh.triangles[l].DQ[1] / P1_data.C[n2]
           + P1_state.v[n3] * mesh.triangles[l].DQ[2] / P1_data.C[n3]) / (2. * A);
      vy = -(P1_state.v[n1] * mesh.triangles[l].DP[0]
           + P1_state.v[n2] * mesh.triangles[l].DP[1]
           + P1_state.v[n3] * mesh.triangles[l].DP[2]) / (2. * A);

/*----------------------------------------------------------------------
      grad u/x -grad v/y*/
      D1 = ux - vy;
/*----------------------------------------------------------------------
      grad v/x +grad u/y*/
      D2 = vx + uy;
      AH = SMAGR * 2.0 * A * sqrt(D1 * D1 + D2 * D2);
      tmp[l] = AH * A / 3;
      for(j = 0; j < 3; j++) {
        n = mesh.triangles[l].vertex[j];
        P1_state.ux[n] += ux * A / 3;
        P1_state.uy[n] += uy * A / 3;
        P1_state.vx[n] += vx * A / 3;
        P1_state.vy[n] += vy * A / 3;
        P1_state.lmts[n] +=  sqrt(D1 * D1 + D2 * D2) * A / 3;
        }
      }
#pragma synchronize
/*BRACKET FOR PARALLEL COMPUTATION*/
    }

  for(l = 0; l < mesh.ntriangles; l++) {
    n1 = mesh.triangles[l].vertex[0];
    n2 = mesh.triangles[l].vertex[1];
    n3 = mesh.triangles[l].vertex[2];
    P1_state.Ah[n1] += tmp[l];
    P1_state.Ah[n2] += tmp[l];
    P1_state.Ah[n3] += tmp[l];
    }

  for(n = 0; n < mesh.nvtxs; n++) {
    factor=1./P1_data.lmm[n];
    P1_state.ux[n] *= factor;
    P1_state.uy[n] *= factor;
    P1_state.vx[n] *= factor;
    P1_state.vy[n] *= factor;
    P1_state.lmts[n] *= factor;
    P1_state.Ah[n] *= factor;
    if(P1_state.Ah[n] < minKh)
      P1_state.Ah[n] = minKh;
    }


  delete[] tmp;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void Hdiffusion2D_LGP1_CLX(mesh_t mesh, state2d_t  P1_state, parameter_t  P1_data)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *-----------------------------------------------------------------------
  Compution of the horizontal diffusion contribution

  LGP1, continuous Galerkin

-----------------------------------------------------------------------*/
  int i, j, l, n, k, m;
  int cnt, *list;
  double CC, CC1, uDQ, uDP, vDQ, vDP;
  double u, v, A;
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Notes:

  01/09/2008:

    We use the non-conservative formulation:

    div Kh grad u

    < div Kh grad u, beta>   = <div beta Kh grad u, 1> - <Kh grad u,grad beta>
                             = <beta Kh grad u,n>      - <Kh grad u,grad beta>


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  for(n = 0; n < mesh.nvtxs; n++) {
    P1_state.Ex[n] = 0.0;
    P1_state.Ey[n] = 0.0;
    }

  for(l = 0; l < mesh.ntriangles; l++) {
    CC = 0.0;
    CC1 = 0.0;
    uDQ = 0.0;
    uDP = 0.0;
    vDQ = 0.0;
    vDP = 0.0;
    A = mesh.triangles[l].Area * (double) 12.;
    for(j = 0; j < 3; j++) {
      n = mesh.triangles[l].vertex[j];
      CC += P1_data.C[n];
      CC1 += 1. / P1_data.C[n];
      u = (double) P1_state.u[n];
      v = (double) P1_state.v[n];
      uDQ += u * mesh.triangles[l].DQ[j];  /* = 2*A dudt */
      uDP += u * mesh.triangles[l].DP[j];
      vDQ += v * mesh.triangles[l].DQ[j];
      vDP += v * mesh.triangles[l].DP[j];
      }
    for(j = 0; j < 3; j++) {
      n = mesh.triangles[l].vertex[j];
      P1_state.Ex[n] += (mesh.triangles[l].DQ[j] * uDQ * CC1 + mesh.triangles[l].DP[j] * uDP * CC) / A;
      P1_state.Ey[n] += (mesh.triangles[l].DQ[j] * vDQ * CC1 + mesh.triangles[l].DP[j] * vDP * CC) / A;
      }
    }

/*----------------------------------------------------------------------
  lmm*C= sum over elts of (AC/3)*/
  for(n = 0; n < mesh.nvtxs; n++) {
    P1_state.Ex[n] *= -P1_state.Ah[n] / (P1_data.lmm[n] * P1_data.C[n]);
    P1_state.Ey[n] *= -P1_state.Ah[n] / (P1_data.lmm[n] * P1_data.C[n]);
    }
}


//#if DEFINED_global_astro_angles
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void spectral_viscosity(spectrum_t spectrum,double duration,mesh_t mesh,parameter_t LGP1_data, int target, double **diffusivity, double **advection)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double c=1.0,R1,R2,R3;
  double a1,G1,a2,G2;
  double zr,zi;
  int k,m,n,nframe,w,status;
  int maxstep=1;
  double *u,*v,*ax,*Gx,*ay,*Gy;
  harmonic_t harmonic;
  date_t start;
  double dt,*time;
  tidal_wave wave;
  double f,Vu,V, V0, omega, norm, advx, advy;
  double **ucs,**usn,**vcs,**vsn,**b[4];
  int nodal=0,nndes,neq,nrhs=1;
  char *task;
  state2d_t state;

  start.year=1950;
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
  dt=2.*M_PI/(spectrum.wave[0].omega*d2r/3600.)/24.0;
  nframe=(int) NINT(duration/dt);
*/

  nndes=mesh.nvtxs;

  time = new double[nframe];
  for(k=0;k<nframe;k++) time[k]=k*dt;
/* *----------------------------------------------------------------------------
  to be replaced*/
  status=harmonic_init01(nframe, start, time, spectrum, &harmonic,1,astro_angles);

  b[0]=new double * [nndes];
  b[1]=new double * [nndes];
  b[2]=new double * [nndes];
  b[3]=new double * [nndes];

  ucs=new double * [nndes];
  usn=new double * [nndes];
  vcs=new double * [nndes];
  vsn=new double * [nndes];

  for(n=0;n<nndes;n++) {
    b[0][n]=new double [harmonic.neq];
    b[1][n]=new double [harmonic.neq];
    b[2][n]=new double [harmonic.neq];
    b[3][n]=new double [harmonic.neq];
    for(k=0;k<harmonic.neq;k++) {
      b[0][n][k]=0;
      b[1][n][k]=0;
      b[2][n][k]=0;
      b[3][n][k]=0;
      }
    ucs[n]=new double [spectrum.n];
    usn[n]=new double [spectrum.n];
    vcs[n]=new double [spectrum.n];
    vsn[n]=new double [spectrum.n];
    }
  neq=harmonic.neq;

  ax= new double [spectrum.n];
  Gx= new double [spectrum.n];
  ay= new double [spectrum.n];
  Gy= new double [spectrum.n];

  for(n=0;n<nndes;n++) {
    for(w = 0; w < spectrum.n; w++) {
      a1=LGP1_data.Utide[n].a[w];
      G1=LGP1_data.Utide[n].G[w];
      a2=LGP1_data.Vtide[n].a[w];
      G2=LGP1_data.Vtide[n].G[w];
      ucs[n][w]=a1*cos(G1);
      vcs[n][w]=a2*cos(G2);
      usn[n][w]=a1*sin(G1);
      vsn[n][w]=a2*sin(G2);
      }
    }

  state.u  =new double[nndes];
  state.v  =new double[nndes];
  state.ux =new double[nndes];
  state.uy =new double[nndes];
  state.vx =new double[nndes];
  state.vy =new double[nndes];
  state.Ah =new double[nndes];
  state.Ex =new double[nndes];
  state.Ey =new double[nndes];
  state.lmts =new double[nndes];


  for(m=0;m<nframe;m++) {
    for(n=0;n<nndes;n++) {
      state.u[n]=0;
      state.v[n]=0;
      for(w = 0; w < spectrum.n; w++) {
        state.u[n]+=ucs[n][w] * harmonic.cs[w][m]+usn[n][w] * harmonic.sn[w][m];
        state.v[n]+=vcs[n][w] * harmonic.cs[w][m]+vsn[n][w] * harmonic.sn[w][m];
        }
      }
/* *----------------------------------------------------------------------------
    to be replaced*/
    smagorinsky_P1NQ(state, LGP1_data, mesh);
    Hdiffusion2D_LGP1_CLX(mesh, state, LGP1_data);
    for(n=0;n<nndes;n++) {
      for(w = 0; w < spectrum.n; w++) {
         b[0][n][2 * w]     += state.Ex[n] * harmonic.cs[w][m];
         b[0][n][2 * w + 1] += state.Ex[n] * harmonic.sn[w][m];
         b[1][n][2 * w]     += state.Ey[n] * harmonic.cs[w][m];
         b[1][n][2 * w + 1] += state.Ey[n] * harmonic.sn[w][m];
         advx=state.u[n]*state.ux[n]+state.v[n]*state.uy[n];
         b[2][n][2 * w]     += advx * harmonic.cs[w][m];
         b[2][n][2 * w + 1] += advx * harmonic.sn[w][m];
         advy=state.u[n]*state.vx[n]+state.v[n]*state.vy[n];
         b[3][n][2 * w]     += advy * harmonic.cs[w][m];
         b[3][n][2 * w + 1] += advy * harmonic.sn[w][m];
         }
       }
    }
  
  task=new char[2];
  for(n=0;n<nndes;n++) {
    for(k=0;k<4;k++) {
      status=poc_getrs(neq, nrhs, harmonic.A, harmonic.pivot, b[k][n]);
      }
    for(w = 0; w < spectrum.n; w++) {
      zr = b[0][n][2 * w];
      zi = b[0][n][2 * w + 1];
      ax[w] = sqrt(zr * zr + zi * zi);
      Gx[w] = atan2(zi, zr);
      zr = b[1][n][2 * w];
      zi = b[1][n][2 * w + 1];
      ay[w] = sqrt(zr * zr + zi * zi);
      Gy[w] = atan2(zi, zr);
      }
    for(w = 0; w < spectrum.n; w++) {
      a1=LGP1_data.Utide[n].a[w];
      G1=LGP1_data.Utide[n].G[w];
      a2=LGP1_data.Vtide[n].a[w];
      G2=LGP1_data.Vtide[n].G[w];
      diffusivity[w][n]=0.5*(ax[w]*a1*cos(Gx[w]-G1)+ay[w]*a2*cos(Gy[w]-G2));
      }
    for(w = 0; w < spectrum.n; w++) {
      zr = b[2][n][2 * w];
      zi = b[2][n][2 * w + 1];
      ax[w] = sqrt(zr * zr + zi * zi);
      Gx[w] = atan2(zi, zr);
      zr = b[3][n][2 * w];
      zi = b[3][n][2 * w + 1];
      ay[w] = sqrt(zr * zr + zi * zi);
      Gy[w] = atan2(zi, zr);
      }
    for(w = 0; w < spectrum.n; w++) {
      a1=LGP1_data.Utide[n].a[w];
      G1=LGP1_data.Utide[n].G[w];
      a2=LGP1_data.Vtide[n].a[w];
      G2=LGP1_data.Vtide[n].G[w];
      advection[w][n]=0.5*(ax[w]*a1*cos(Gx[w]-G1)+ay[w]*a2*cos(Gy[w]-G2));
      }
    }

  for(n=0;n<nndes;n++) {
    delete[] b[0][n];
    delete[] b[1][n];
    delete[] b[2][n];
    delete[] b[3][n];
    delete[] ucs[n];
    delete[] usn[n];
    delete[] vcs[n];
    delete[] vsn[n];
    }

  delete[] b[0];
  delete[] b[1];
  delete[] b[2];
  delete[] b[3];

  delete[] ucs;
  delete[] usn;
  delete[] vcs;
  delete[] vsn;

  delete[] ax;
  delete[] ay;
  delete[] Gx;
  delete[] Gy;

//  state.destroy();

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void spectral_advection(spectrum_t spectrum,double duration,state2d_t state, parameter_t data, int nndes, int target, double **buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/****************************************************************
*****************************************************************/
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

  start.year=1990;
  start.month=1;
  start.day=1;
  start.second=0;

  nframe=365*24;
  dt=1800.0;
  nframe=(int) NINT(duration/dt);

  duration=3600.*24*15;
  dt=2.*M_PI/(spectrum.waves[0].omega*d2r/3600.)/24.0;
  nframe=(int) NINT(duration/dt);

//  spectrum.n=4;

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

/* *----------------------------------------------------------------------------
  to be replaced*/
//  status=harmonic_init01(nframe, start, time, spectrum, &harmonic,0);

  for(n=0;n<nndes;n++) {
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
      c=sqrt(u[m]*u[m]+v[m]*v[m]) + data.u0[n];
      taux[m]=c*u[m]/H[m];
      tauy[m]=c*v[m]/H[m];
      }
//    status=harmonic_compute(u,residuals, nframe, harmonic, spectrum, maxstep,ax,Gx);
//    status=harmonic_compute(v,residuals, nframe, harmonic, spectrum, maxstep,ay,Gy);
/* *----------------------------------------------------------------------------
  to be replaced*/
//    status=harmonic_compute(taux,residuals, nframe, harmonic, spectrum, maxstep,ax,Gx);
/* *----------------------------------------------------------------------------
  to be replaced*/
//    status=harmonic_compute(tauy,residuals, nframe, harmonic, spectrum, maxstep,ay,Gy);
//    if(data[n].Ua[0] > .2) status=fft(taux, nframe);
    for(w = 0; w < spectrum.n; w++) {
      a1=data.Utide[n].a[w];
      G1=data.Utide[n].G[w];
      a2=data.Vtide[n].a[w];
      G2=data.Vtide[n].G[w];
      buffer[w][n]=0.5*(ax[w]*a1*cos(Gx[w]-G1)+ay[w]*a2*cos(Gy[w]-G2));
      }
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



