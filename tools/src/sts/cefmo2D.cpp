
/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2012

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France

***************************************************************************/

#include "tugo-prototypes.h"
#include "dry.h"
#include "ice.h"
#include "tmatrix.h"
#include "cefmo.h"
#include "spectrum.h"
#include "matrix.h"
#include "poc-netcdf.hpp"
#include "parallel.h"

extern int  LinearSystem_terminate(hyperzmatrix_t *A);
extern int initialise_compound(compound_t **compound, int & ncompound, spectrum_t WaveList);

extern void SpectralMomentum3D_UBASExZBASE(mesh_t mesh, cefmo_t *cefmo, tidal_wave wave, parameter_t data2D, tide2D_t state2D /*, tide3D_t state3D*/);

//#define ROUTINE_TRACK
/**----------------------------------------------------------------------------
runtime matrices, to be further reconsidered */
// #define M1 cefmo.M1
// #define M2 cefmo.M2

// #define G1 cefmo.G1
// #define G2 cefmo.G2

#define D1 cefmo.D1
#define D2 cefmo.D2

/**----------------------------------------------------------------------------
 to be further reconsidered */
#define IMPLICIT_ZEROFLUX 
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 \todo do it for 3D too... */

#define Z_NNPE cefmo.z_descriptor->nnpe

double *CdRatio;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void cefmo2D_implicitFBCs01(complex<double> **SpMu, complex<double> **SpMv)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,n;
  double Tx,Ty;
  complex <double> Axx,Ayy,Axy,Ayx,Aij,det;
  
#ifdef IMPLICIT_ZEROFLUX
/**----------------------------------------------------------------------------
  set boundary condition in momentum matrix*/
  for(i = 0; i < gSpectralRigidBCs.nnodes; i++) {
    n= gSpectralRigidBCs.nodes[i];
    Tx = -gSpectralRigidBCs.data[i].ny;
    Ty =  gSpectralRigidBCs.data[i].nx;
    if(fabs(Tx)>fabs(Ty)) {
/**----------------------------------------------------------------------------
      tangent projection of dynamics equation*/
      Axx = Tx * SpMu[n][0] + Ty * SpMv[n][0];
      Ayy = Tx * SpMu[n][1] + Ty * SpMv[n][1];
      SpMu[n][0] =  Axx;
      SpMu[n][1] =  Ayy;
/**----------------------------------------------------------------------------
      normal flux condition*/
      SpMv[n][0] = -Ty*diagonal;
      SpMv[n][1] =  Tx*diagonal;
      }
    else {
/**----------------------------------------------------------------------------
      tangent projection of dynamics equation*/ /// HERE !!!
      Axx = Tx * SpMu[n][0] + Ty * SpMv[n][0];
      Ayy = Tx * SpMu[n][1] + Ty * SpMv[n][1];
      SpMu[n][0] =  Axx;
      SpMu[n][1] =  Ayy;
/**----------------------------------------------------------------------------
      normal flux condition*/
      SpMv[n][0] = -Ty*diagonal;
      SpMv[n][1] =  Tx*diagonal;
      }
    }
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void cefmo2D_implicitFBCs02(cefmo_t & cefmo)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,n;
  double Tx,Ty;
  complex <double> Axx,Ayy,Axy,Ayx,Aij,det;
  
#ifdef IMPLICIT_ZEROFLUX
/**----------------------------------------------------------------------------
  modify inverse matrix to include rhs projection*/
  for(i = 0; i < gSpectralRigidBCs.nnodes; i++) {
    n= gSpectralRigidBCs.nodes[i];
    Tx = -gSpectralRigidBCs.data[i].ny;
    Ty =  gSpectralRigidBCs.data[i].nx;
    Axx=cefmo.InvSpMu[n][0];
    Axy=cefmo.InvSpMu[n][1];
    Ayx=cefmo.InvSpMv[n][0];
    Ayy=cefmo.InvSpMv[n][1];
    cefmo.InvSpMu[n][0] =  Axx*Tx;
    cefmo.InvSpMu[n][1] =  Axx*Ty;
    cefmo.InvSpMv[n][0] =  Ayx*Tx;
    cefmo.InvSpMv[n][1] =  Ayx*Ty;
    }
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void SpMomentum2D_RHS(mesh_t & mesh, cefmo_t & cefmo, parameter_t & data, double theta, complex <double> *rhsU, complex <double> *rhsV)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Compute the spectral momentum rhs

  Transport version

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i, j, m, n;
  int status;
  int udim,zdim;
  complex <double> *buf;
  complex <double> *prx,*pry;
  complex <double> sumx,sumy,mean;
  double dxdt,dydt,dxdp,dydp,factor;
  double C, Lmm, H;
  double S;
  complex <double> tmp,dzdx,dzdy;
  actionZ_t action;
  double pseudo;
  double h[Z_NNPE];
  complex <double> p[Z_NNPE],z[Z_NNPE];

  discretisation_t & z_descriptor=*cefmo.z_descriptor;
  discretisation_t & u_descriptor=*cefmo.u_descriptor;

  status=paire_dimension(mesh,cefmo.paire,&zdim,&udim);
  
/**-----------------------------------------------------------------------
  compute gradients of :
    -elevation
    -astronomic potential
    -atmospheric pressure
    -loading/self-attraction and 
    -ocean bottom vertical displacement 
-----------------------------------------------------------------------*/
  buf = new complex<double>[zdim];

  for(n = 0; n < zdim; n++) {
    buf[n] = theta*cefmo.state.z[n];
    }

  if(LSA_forcing) {
    for(n = 0; n < zdim; n++) {
      buf[n] -= cefmo.state.LSA[n];
      }
    }

  if(astronomic_forcing==1) {
    for(n = 0; n < zdim; n++) {
      buf[n] -= cefmo.state.potential[n];
      }
    }

/**-----------------------------------------------------------------------
  compute pressure gradients*/
  action=cefmo.action;

//   for(n=0;n<udim;n++) {
//     action.prx[n][0] = 0.0;
//     action.pry[n][0] = 0.0;
//     }
// 
//   for(m = 0; m < mesh.ntriangles; m++) {
//     C=mesh.triangles[m].cosinus;
//     mean=0.0;
//     for(j = 0; j < z_descriptor.nnpe; j++) {
//       n=z_descriptor.NIbE[m][j];
//       mean += buf[n]/(double) Z_NNPE;
//       }
//     for(j = 0; j < z_descriptor.nnpe; j++) {
//       n=z_descriptor.NIbE[m][j];
//       h[j]=cefmo.state.h_z[n];
//       p[j]=buf[n]-mean;
//       }
//     dzdx = 0.;
//     dzdy = 0.;
//     for(j = 0; j < 3; j++) {
//       n=z_descriptor.NIbE[m][j];
//       dzdx += p[j] * mesh.triangles[m].DQ[j];
//       dzdy -= p[j] * mesh.triangles[m].DP[j];
//       }
//     sumx=dzdx*fe_integraleLGP1_2D(h);
//     sumy=dzdy*fe_integraleLGP1_2D(h);
// /**-----------------------------------------------------------------------
//     integration kept unchanged...*/
//     action.prx[m][0] +=   P_g *sumx;
//     action.pry[m][0] += C*P_g *sumy;
//     }
  
  prx = new complex<double>[udim];
  pry = new complex<double>[udim];
  
  STDERR_BASE_LINE("buf=========> %g,%g\n",buf[3].real(),buf[3].imag());
  STDERR_BASE_LINE("cefmo.G1=========> %g\n",cefmo.G1.packed[0]);
  status=matrix_operation(cefmo.G1,buf,prx,1);
  STDERR_BASE_LINE("prx=========> %g,%g\n",prx[3].real(),prx[3].imag());
  status=matrix_operation(cefmo.G2,buf,pry,1);
  for(m = 0; m < udim; m++) {
    action.prx[m][0]=prx[m];
    action.pry[m][0]=pry[m];
    }
  
  zaparr(buf);
  
  delete[] prx;
  delete[] pry;

  for(n = 0; n < udim; n++) {
//     C   = (double) u_descriptor.nodes[n].c;
//     S   = (double) u_descriptor.nodes[n].s;
//     H   = cefmo.state.h_u[n];
    Lmm = cefmo.u_LumpedMW[n];

    tmp = 0;

/*#################### Pressure gradient contribution ##################
----------------------------------------------------------------------*/
    tmp -=       action.prx[n][0];

    tmp += Lmm * action.Fx[n];

    rhsU[n] = tmp;

    tmp = 0;

/*#################### Pressure gradient contribution ##################
----------------------------------------------------------------------*/
    tmp -=      action.pry[n][0];

    tmp += Lmm * action.Fy[n];

    rhsV[n] = tmp;
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void SpMomentum2D_Initialise(mesh_t & mesh, cefmo_t & cefmo, parameter_t & data, tidal_wave wave)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

    G1 and G2 are gradient operators

    M1 and M2 are momentum time gradient operators

    A(u,v)=(Fu,Fv)

    a b     d -b
    c d    -c  a  * (ad-bc)^-1

    Axx -Aij    Ayy  Aij
    Aij  Ayy   -Aij  Axx / (Axx²+Aij²)  NB: Axx=Ayy

  24/09/2009 :

    Implicit boundary contions need to be added by modifying equations:

    Normal flux condition:

      u * nx + v * ny =0

    Projection of 2D momentum equation along the tangential direction:

      A(u,v).t=(Fu,Fv).t

      (Axx * u  + Axy * v) * tx + (Ayx *u + Ayy * v) * ty   = Fu * tx +Fv *ty

      (Axx * tx + Ayx * ty) * u + (Axy * tx + Ayy * ty) * v = ...

                     | Axx * tx + Ayx * ty     Axy * tx + Ayy * ty |
     modified  A =   |                                             |
                     |                  nx                      ny |

    Projection :

      A u = P F

            |  tx     ty |
      P =   |            |
            |   0     0  |

            |  Bxx * tx         Bxx * ty  |
      B P = |                             |
            |  Byx * tx         Byx * ty  |

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int i, j, k, l, m, n, row;
  int status;
  int udim,zdim;
  double S, Lmm, Cd, h;
  complex <double> Axx,Ayy,Axy,Ayx,Aij,det;
  double Tx,Ty,C;
  ordering_t ordering;
  int *pointer,*incidence,*cardinal;
  zmatrix2x2_t *FrictionMatrix;
  matrix2x2_t  *DragMatrix;
  complex<double> J=complex<double>(0.,1.);
  double f;

  discretisation_t & u_descriptor=*cefmo.u_descriptor;
  discretisation_t & z_descriptor=*cefmo.z_descriptor;

#ifdef ROUTINE_TRACK
  printf("entering %d in %s\n",__LINE__, __FILE__);
#endif
  
  status=paire_dimension(mesh,cefmo.paire,&zdim,&udim);
  
  if(cefmo.M1.ordering==0) {
/**----------------------------------------------------------------------------
    to be done only once */ 
    cefmo.M1.ordering= new ordering_t;

    status=fe_crosstables(&mesh);
    status=cefmo.packed_UBASExZBASE(mesh, u_descriptor, cefmo.M1.ordering);

    /* STS: NO periodic */
    
    cefmo.M2.ordering=cefmo.M1.ordering;

    pointer   = cefmo.M1.ordering->pointer;
    incidence = cefmo.M1.ordering->incidence;
    cardinal  = cefmo.M1.ordering->cardinal;

    cefmo.M1.packed=new complex<double> [pointer[udim]];
    cefmo.M2.packed=new complex<double> [pointer[udim]];

    for(n = 0; n < pointer[udim]; n++) {
      cefmo.M1.packed[n]=0;
      }
    for(n = 0; n < pointer[udim]; n++) {
      cefmo.M2.packed[n]=0;
      }

/**----------------------------------------------------------------------------
    */ 
    cefmo.G1.packed=new double[pointer[udim]];
    cefmo.G2.packed=new double[pointer[udim]];

    cefmo.G1.ordering=cefmo.M1.ordering;
    cefmo.G2.ordering=cefmo.M1.ordering;

    for(n = 0; n < pointer[udim]; n++) {
      cefmo.G1.packed[n]=0;
      }
    for(n = 0; n < pointer[udim]; n++) {
      cefmo.G2.packed[n]=0;
      }

    cefmo.InvSpMu = tmatrix(cefmo.InvSpMu, complex <double> (1.e+10,0.), udim, 2);
    cefmo.InvSpMv = tmatrix(cefmo.InvSpMv, complex <double> (1.e+10,0.), udim, 2);

/**-----------------------------------------------------------------------
    compute pressure gradients matrices G1 & G2*/
    status=cefmo.SpGradient2D( mesh, cefmo, data);
/**-----------------------------------------------------------------------
    periodic mesh add-on*/
    if(mesh.periodic==1) {
      /* STS: NO periodic */
      }
   
    status=matrix_print("matrixG1-3.dat",cefmo.G1);
  STDERR_BASE_LINE("cefmo.G1=========> %g\n",cefmo.G1.packed[0]);
    status=matrix_print("matrixG2-3.dat",cefmo.G2);

    }
    
  pointer   = cefmo.M1.ordering->pointer;
  incidence = cefmo.M1.ordering->incidence;
  cardinal  = cefmo.M1.ordering->cardinal;
    
  FrictionMatrix = cefmo.FrictionMatrix;
  DragMatrix     = cefmo.DragMatrix;

  cefmo.SpMu = tmatrix(cefmo.SpMu, complex <double> (1.e+10,0.), udim, 2);
  cefmo.SpMv = tmatrix(cefmo.SpMv, complex <double> (1.e+10,0.), udim, 2);
  STDERR_BASE_LINE("cefmo.SpMv=========> %g\n",cefmo.SpMv[0]->real());

/**-----------------------------------------------------------------------
  compute direct spectral momentum matrix*/
  for(n = 0; n < udim; n++) {
    C   = (double) u_descriptor.nodes[n].c;
    S   = (double) u_descriptor.nodes[n].s;
    h   = cefmo.state.h_u[n];
    Lmm = cefmo.u_LumpedMW[n];
    f = S * two_Omega;

    Cd=data.Cd[n];
    complex<double> r1, r2, r3, r4;
    double r=data.rlinear[n];
    friction_coefficient(FrictionMatrix[n], Cd, r, h, &r1, &r2, &r3, &r4);
    cefmo.SpMu[n][0] = Lmm * C * (wave.omega*J*dph2rps + r1 - DragMatrix[n].c[0][0]);
    cefmo.SpMu[n][1] = Lmm * C * (-f                   + r2 - DragMatrix[n].c[0][1]);
    cefmo.SpMv[n][0] = Lmm * C * (+f                   + r3 - DragMatrix[n].c[1][0]);
    cefmo.SpMv[n][1] = Lmm * C * (wave.omega*J*dph2rps + r4 - DragMatrix[n].c[1][1]);
    }
  STDERR_BASE_LINE("cefmo.SpMv=========> %g\n",cefmo.SpMv[0]->real());

/**-----------------------------------------------------------------------
  periodic mesh add-on*/
//   if(mesh.periodic==1 && cefmo.paire!=0) {
  if(mesh.periodic==1) {
    /* STS: NO periodic */
    }

/**-----------------------------------------------------------------------
  implicit rigid BC, first step*/
  cefmo2D_implicitFBCs01(cefmo.SpMu, cefmo.SpMv);
  STDERR_BASE_LINE("cefmo.SpMv=========> %g\n",cefmo.SpMv[0]->real());
  
/**----------------------------------------------------------------------------
  compute inverse momentum matrix*/
  for(n = 0; n < udim; n++) {
    Axx=cefmo.SpMu[n][0];
    Axy=cefmo.SpMu[n][1];
    Ayx=cefmo.SpMv[n][0];
    Ayy=cefmo.SpMv[n][1];
    det=Axx*Ayy-Axy*Ayx;
    cefmo.InvSpMu[n][0] =  Ayy/det;
    cefmo.InvSpMu[n][1] = -Axy/det;
    cefmo.InvSpMv[n][0] = -Ayx/det;
    cefmo.InvSpMv[n][1] =  Axx/det;
    }

/**-----------------------------------------------------------------------
  implicit rigid BC, second step*/
  cefmo2D_implicitFBCs02(cefmo);

/*-----------------------------------------------------------------------------
  compute M1 & M2 as product of A and G1 & G2: hu=[AG1]z*/
  for(n = 0; n < udim; n++) {
    for(k = 0; k < cardinal[n]; k++) {
      const double g1=cefmo.G1.packed[pointer[n]+k];
      const double g2=cefmo.G2.packed[pointer[n]+k];
      cefmo.M1.packed[pointer[n]+k]=g1*cefmo.InvSpMu[n][0]+g2*cefmo.InvSpMu[n][1];
      cefmo.M2.packed[pointer[n]+k]=g1*cefmo.InvSpMv[n][0]+g2*cefmo.InvSpMv[n][1];
      }
    }

}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int SpWE2D_Initialise(mesh_t & mesh,cefmo_t & cefmo, parameter_t & data, hyperzmatrix_t & elevation, tidal_wave wave, int initialise)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------

  Semi-implicit continuity equation, Walter-Casulli type of

  Transport version

 -----------------------------------------------------------------------*/
{
  int status;
  int i, m, n;
  int udim,zdim;
  complex<double> J=complex<double>(0.,1.);

  hypermatrix_t  massmatrix;

  discretisation_t & u_descriptor=*cefmo.u_descriptor;
  discretisation_t & z_descriptor=*cefmo.z_descriptor;
  
  status=paire_dimension(mesh,cefmo.paire,&zdim,&udim);

  if(initialise==1) {
/**-----------------------------------------------------------------------
    packed matrix for ZBASE x UBASE divergence*/
    D1.ordering=new ordering_t;
    status=cefmo.packed_ZBASExUBASE(mesh, z_descriptor, D1.ordering);
    
    /* STS: NO periodic */

    D2.ordering=D1.ordering;

/*-----------------------------------------------------------------------
    allocate arrays*/
    D1.packed=new double[D1.ordering->pointer[zdim]];
    D2.packed=new double[D2.ordering->pointer[zdim]];

/*-----------------------------------------------------------------------
    initialize arrays*/
    for(n = 0; n < D1.ordering->pointer[zdim]; n++) {
      D1.packed[n]=0;
      }
    for(n = 0; n < D1.ordering->pointer[zdim]; n++) {
      D2.packed[n]=0;
      }

    status=cefmo.SpDivergence2D( mesh, cefmo, data);

/**-----------------------------------------------------------------------
    WE square packed matrix for ZBASE discretisation*/
    elevation.ordering=new ordering_t;

//     status= cefmo.packed_ZBASExZBASE(mesh, z_descriptor, elevation.ordering);
    status= ordering_import(D1.ordering,cefmo.M1.ordering,elevation.ordering);

    elevation.distributor=&z_descriptor.distributor;
    elevation.packed=new complex<double> [elevation.ordering->pointer[zdim]];
    }
  else {
//    elevation.ordering->destroy();
//    status= ordering_import(D1.ordering,cefmo.M1.ordering,elevation.ordering);e
    if(elevation.packed==0) elevation.packed=new complex<double> [elevation.ordering->pointer[zdim]];
    }
  
/**-----------------------------------------------------------------------
  initialise with D*A(-1)*G = D1*M1+D2*M2 */
  status= matrix_product(D1,cefmo.M1,elevation,1);
  status= matrix_product(D2,cefmo.M2,elevation,0);

  for(m = 0; m < elevation.ordering->pointer[zdim]; m++) {
    elevation.packed[m]=-elevation.packed[m];
    }

// #define ROUTINE_TRACK

/**-----------------------------------------------------------------------
  add mass matrix*/
#ifdef ROUTINE_TRACK
  printf("add massmatrix at %d in %s\n",__LINE__, __FILE__);
#endif
  if(mesh.periodic==0)
    massmatrix=z_descriptor.massmatrix;
  else {
    massmatrix.ordering=new ordering_t;
    status=initialise_MassMatrix(mesh, z_descriptor, &massmatrix, 0);
    }
  status = matrix_axpy(elevation,massmatrix,J*wave.omega*dph2rps);
#ifdef ROUTINE_TRACK
  printf("elevation matrix done at %d in %s\n",__LINE__, __FILE__);
#endif
  
  elevation.discretisation=z_descriptor.type;

#ifdef ROUTINE_TRACK
  printf("leaving at %d in %s\n",__LINE__, __FILE__);
#endif
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void SpWE2D_RHS(mesh_t & mesh, cefmo_t & cefmo, parameter_t & data, complex<double>  *rhs)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------

  Spectral elevation solver

 -----------------------------------------------------------------------*/
{
  int status;
  int i, j, k, l, m, n;
  int zdim,udim;
  complex<double> J=complex<double>(0.,1.);
  complex<double>  *Fz,*Fu,*Fv;
  double theta;

  status=paire_dimension(mesh,cefmo.paire,&zdim,&udim);
/*-----------------------------------------------------------------------
  allocate arrays*/
  Fu =new complex<double> [udim];
  Fv =new complex<double> [udim];

/**-----------------------------------------------------------------------

  Compute the partial transport U'=hu': theta is weight on z
    At this stage, state[2].U,state[2].V is partial transport

-----------------------------------------------------------------------*/
  theta=0.0;
  STDERR_BASE_LINE("pot=========> %g,%g\n",cefmo.state.potential[3].real(),cefmo.state.potential[3].imag());
  STDERR_BASE_LINE("cefmo.G1=========> %g\n",cefmo.G1.packed[0]);
  SpMomentum2D_RHS   (mesh, cefmo, data,theta,Fu,Fv);
  STDERR_BASE_LINE("Fu=========> %g,%g\n",Fu[3].real(),Fu[3].imag());
  SpMomentum2D_solve (mesh, cefmo, Fu,Fv);

/**-----------------------------------------------------------------------
  compute -Hu'.grad(beta)*/
  STDERR_BASE_LINE("u=========> %g,%g\n",cefmo.state.u[3].real(),cefmo.state.u[3].imag());
  status=matrix_operation(D1,cefmo.state.u,rhs,1);
  status=matrix_operation(D2,cefmo.state.v,rhs,0);

/**-----------------------------------------------------------------------
  apply -1 factor (i.e. move to rhs)*/
  for(n = 0; n < zdim; n++) {
    rhs[n]=-rhs[n];
    }

/*-----------------------------------------------------------------------
  deallocate arrays*/
  zaparr(Fu);
  zaparr(Fv);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void SpMomentum2D_Solver(mesh_t & mesh, cefmo_t & cefmo, parameter_t & data, tidal_wave wave)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------

  Spectral momentum solver

 -----------------------------------------------------------------------*/
{
  int status;
  int m,n;
  int zdim,udim;
  complex<double>  *Fu,*Fv;
  double theta;
  extern void SpMomentum2D_OBCPatch(mesh_t & mesh, cefmo_t & cefmo, parameter_t & data, tidal_wave wave, complex<double>  *Fu, complex<double>  *Fv);

  status=paire_dimension(mesh,cefmo.paire,&zdim,&udim);
/*-----------------------------------------------------------------------
  allocate arrays*/
  Fu =new complex<double> [udim];
  Fv =new complex<double> [udim];

/**-----------------------------------------------------------------------
  Compute spectral transport (h*u, h*v) */
  theta=1.0;
  
  SpMomentum2D_RHS   (mesh, cefmo, data,theta,Fu,Fv);

//   status=check_periodic2D(mesh, cefmo.u_descriptor->paires, Fu);
//   status=check_periodic2D(mesh, cefmo.u_descriptor->paires, Fv);
  
  SpMomentum2D_solve (mesh, cefmo, Fu,Fv);

/**-----------------------------------------------------------------------
  patch to adjust flux divergence at open boundary nodes */
  if(cefmo.z_descriptor->type==CQP0)
    SpMomentum2D_OBCPatch(mesh, cefmo, data, wave, Fu, Fv);
  
/**-----------------------------------------------------------------------
  convert h*u into u */
  for(n = 0; n < udim; n++) {
    cefmo.state.u[n]=cefmo.state.u[n]/cefmo.state.h_u[n];
    cefmo.state.v[n]=cefmo.state.v[n]/cefmo.state.h_u[n];
    }

/*-----------------------------------------------------------------------
  deallocate arrays*/
  zaparr(Fu);
  zaparr(Fv);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cefmo2D_compute(mesh_t & mesh, parameter_t & data, cefmo_t & cefmo, tidal_wave wave, int wave_id, double *N, complex<double> *rhs, complex<double> *divergence, int iteration)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int i,j,k,l,m,n,status;
  int zdim,udim;
  int z_discretisation, u_discretisation;
  int initialized;
  tide2D_t M2state;
  char *comment[2],filename[1024],varname[1024];
  double *zbuffer,*u_buffer;
  bool force_duplicate=false;
  int verbose=0;
  bool debug=false;

  status=paire_dimension(mesh,cefmo.paire,&zdim,&udim);
  status=paire_discretisation_id(cefmo.paire, &z_discretisation, &u_discretisation);
  
  const char *UNAME=discretisation_name(u_discretisation);
  const char *ZNAME=discretisation_name(z_discretisation);

  discretisation_t & u_descriptor=*cefmo.u_descriptor;
  discretisation_t & z_descriptor=*cefmo.z_descriptor;

/**-----------------------------------------------------------------------
  Solve for astronomical constituents*/
  STDOUT_BASE_LINE("wave %s, iteration %d\n",wave.name,iteration);

//   sprintf(filename, "%s.spectral.nc",wave.name);
//   
//   sprintf(varname,"h_eta_%s", UNAME);
//   status=archiving_UGdummy2D(filename, mesh, varname, cefmo.state.h_u,(int) u_discretisation);
//   sprintf(varname,"h_eta_%s", ZNAME);
//   status=archiving_UGdummy2D(filename, mesh, varname, cefmo.state.h_z,(int) z_discretisation);

/*------------------------------------------------------------------------
  internal wave drag*/
  cefmo.DragMatrix = SpectralDrag_init(u_descriptor, wave, N,cefmo.state.h_u,data, gwdHmin, gwdHmax, udim);

/**-----------------------------------------------------------------------------
  initialize tidal astronomical potential*/
  if(astronomic_forcing==1) {
    tidal_equilibrium_harmonics(wave, data, zdim, cefmo.state.potential);
    STDERR_BASE_LINE("pot=========> %g,%g\n",cefmo.state.potential[3].real(),cefmo.state.potential[3].imag());
    }

/**-----------------------------------------------------------------------------
  initialize tidal loading and self-attraction*/
  if(LSA_forcing) {
    for(n = 0; n < zdim; n++) {
      cefmo.state.LSA[n]=data.LSA[n].z[wave_id];
      }
    }

  initialized=(cefmo.SpMatrix.ordering!=0);

  switch (initialized) {
    case 0:
/**-----------------------------------------------------------------------------
      initialize momentum inverse matrix*/
      SpMomentum2D_Initialise(mesh, cefmo, data, wave);
/**-----------------------------------------------------------------------------
      initialize the spectral elevation matrix*/
      cefmo.SpMatrix.ordering=new ordering_t;
      status=SpWE2D_Initialise( mesh, cefmo, data, cefmo.SpMatrix, wave,1);
      break;
    default:
/**-----------------------------------------------------------------------------
      initialize momentum inverse matrix*/
      SpMomentum2D_Initialise(mesh, cefmo, data, wave);
/**-----------------------------------------------------------------------------
      initialize the spectral elevation matrix*/
      status=SpWE2D_Initialise( mesh, cefmo, data, cefmo.SpMatrix, wave,0);
      break;
    }

/*-----------------------------------------------------------------------------
  factorize the spectral elevation matrix*/
//   status=LinearSystem_initialise(&cefmo.SpMatrix, cefmo.clamped, cefmo.nclamped, gSpectralSolver, force_duplicate, verbose, debug);
  status=LinearSystem_initialise(&cefmo.SpMatrix, cefmo.clamped, cefmo.nclamped, gSpectralSolver);
  if(status!=0) {
    check_error(-1, "solver initialisation/factorisation failed", __LINE__, __FILE__, 1);
    }
    
/*-----------------------------------------------------------------------------
  initialize the spectral elevation RHS*/
  STDERR_BASE_LINE("cefmo.G1=========> %g\n",cefmo.G1.packed[0]);
  SpWE2D_RHS(mesh, cefmo,  data, rhs);
  if(divergence!=0) {
    for(n = 0;n < zdim; n++) {
      rhs[n]-=divergence[n];
      }
    }

/*-----------------------------------------------------------------------------
  assign boundary conditions to spectral elevation RHS*/
  for(k = 0; k < gSpectralOpenBCs.nnodes; k++) {
    n=gSpectralOpenBCs.nodes[k];
    rhs[n]=diagonal*zconvert(gSpectralOpenBCs.tides[k].ztide,wave_id);
    }
  
  cefmo.SpMatrix.distributor=&(z_descriptor.distributor);

/*-----------------------------------------------------------------------------
  solve spectral elevation system*/
  status=LinearSystem_solve(cefmo.SpMatrix,rhs);
  for(n = 0; n < zdim; n++) {
    cefmo.state.z[n]=rhs[n];
//     if(rhs[n]==0.) {
//       printf("elevation spectral solution is zero for n=%d (%d %d %d)\n",n, cefmo.SpMatrix.d->gindex[n], z_descriptor.nodes[n].rank, z_descriptor.nodes[n].globalized);
//       check_error(-1, "abort...", __LINE__, __FILE__, 1);
//       }
    if(isnan(real(rhs[n]))!=0) {
      printf("elevation spectral solution is nan for n=%d\n",n);
      check_error(-1, "abort...", __LINE__, __FILE__, 1);
      }
    }

/**-----------------------------------------------------------------------------
  solve spectral velocity system*/
  SpMomentum2D_Solver( mesh,  cefmo,  data, wave);

//   status=check_periodic2D(mesh, u_descriptor.paires, cefmo.state.u);
//   status=check_periodic2D(mesh, u_descriptor.paires, cefmo.state.v);

/**-----------------------------------------------------------------------------
  archive spectral solution*/
  SpArchive_Put2D(mesh, cefmo, data, wave, iteration);

/**-----------------------------------------------------------------------------
  energy budgets*/
  if(gSpectralEnergy == 1) {
    /* STS: NO energy budgets */
    }
/*-----------------------------------------------------------------------------
  free solver*/
  status= LinearSystem_terminate(&(cefmo.SpMatrix));

/*-----------------------------------------------------------------------------
  free wave drag matrix*/
  delete[] cefmo.DragMatrix;
  cefmo.DragMatrix=0;

  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cefmo2D(mesh_t mesh, parameter_t data, cefmo_t cefmo, atlas2D_t *atlas,int already_initialised)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------

  Implicit wave equation

 -----------------------------------------------------------------------*/
{
  int i,j,k,l,m,n,status;
  int zdim,udim;
  int wave_id;
  int K1_id,M2_id;
  int iteration,count;
  int max_iteration;
  tidal_wave wave;
  complex <double> *rhs;
  double *N;
  double *zbuffer,*u_buffer;
  double *h;
  tide2D_t M2state,K1state;
  char *comment[2],filename[1024],msg[128];
  int z_discretisation, u_discretisation;
  complex<double> *divergence=0;
  double duration=24*20*3600.;
  complex<double> **bufx, **bufy;
  complex<double> mask=1.e+10;

  pair<map<const char *,tide2D_t>::iterator,bool> atlas_status;
  int ncompound;

  spectrum_t solved(WaveList.n);
  spectrum_t prediction(100), analysis(100);
  int seq_zdim,seq_udim;
  
  status=paire_dimension(mesh,cefmo.paire,&zdim,&udim);
  status=paire_dimension(mesh,gSequentialPaire2D ,&seq_zdim,&seq_udim);
  status=paire_discretisation_id(cefmo.paire, &z_discretisation, &u_discretisation);

  discretisation_t & u_descriptor=*cefmo.u_descriptor;
  discretisation_t & z_descriptor=*cefmo.z_descriptor;

  max_iteration=tugo_cfg->tides->Spectral2DMaxIteration.numerical_value();

  status= cefmo_allocate(mesh, cefmo, already_initialised);

/**-----------------------------------------------------------------------
  Brünt-Vassala frequency needed for internal tide drag*/
  N=new double[udim];
  u_buffer=recast(gstate2D[0][0].N,seq_udim);
  status=fe_projection(mesh, u_buffer, u2D_discretisation, N, u_discretisation);
  delete[] u_buffer;
  
/**-----------------------------------------------------------------------
  topography and topography gradient needed for internal tide drag*/
  status=cefmo_topography(mesh, cefmo.state, data, cefmo.paire);
  
/**-----------------------------------------------------------------------
  For academic cases */
  status=coriolis_setup(z_descriptor);
  status=coriolis_setup(u_descriptor);

/**-----------------------------------------------------------------------
  import time-stepping parameters */
  status=fe_projection(mesh, gdata2D[0].Cd, u2D_discretisation, data.Cd, u_discretisation);
  status=fe_projection(mesh, gdata2D[0].z0, u2D_discretisation, data.z0, u_discretisation);
  status=fe_projection(mesh, gdata2D[0].u0, u2D_discretisation, data.u0, u_discretisation);
  status=fe_projection(mesh, gdata2D[0].FrictionRatio,u2D_discretisation, data.FrictionRatio, u_discretisation);
  
  status=SpectralFriction_init(u_descriptor, wave, data, udim);
  
  CdRatio=new double[u_descriptor.nnodes];
  for(n = 0; n < udim; n++) {
    CdRatio[n]=data.FrictionRatio[n];
    }
  status=friction_ratio(mesh, u_descriptor, data.Cd, CdRatio);

  rhs=new complex <double>[zdim];

/**-----------------------------------------------------------------------
  identify dominant wave*/
  K1_id=WaveList.wave_index(cefmo.dominants[1]);
  if(K1_id==-1) {
    sprintf(msg,"%s not included in tidal spectrum",cefmo.dominants[1].name);
    check_error(-1, msg, __LINE__, __FILE__, 0);
    cefmo.ndominants=1;
    }
  else {
    cefmo.ndominants=2;
    }

/**-----------------------------------------------------------------------
  identify dominant wave*/
  M2_id=WaveList.wave_index(cefmo.dominants[0]);
  if(M2_id==-1) {
    sprintf(msg,"%s not included in tidal spectrum",cefmo.dominants[0].name);
    check_error(-1, msg, __LINE__, __FILE__, 1);
    }

  M2state.z  =new complex <double>[zdim];
  M2state.u  =new complex <double>[udim];
  M2state.v  =new complex <double>[udim];

  K1state.z  =new complex <double>[zdim];
  K1state.u  =new complex <double>[udim];
  K1state.v  =new complex <double>[udim];

  for(n = 0; n < zdim; n++) {
    M2state.z[n]=0;
    K1state.z[n]=0;
    }
  for(n = 0; n < udim; n++) {
    M2state.u[n]=complex <double>(1.0,0.0);
    M2state.v[n]=complex <double>(0.0,0.0);
    K1state.u[n]=complex <double>(0.5,0.0);
    K1state.v[n]=complex <double>(0.0,0.0);
    if (already_initialised==1){
      M2state.u[n]=cefmo.state.u[n];
      M2state.v[n]=cefmo.state.v[n];
      }
    }

/*-----------------------------------------------------------------------------
  initialise open boundary condition structures*/
  treat_boundarycode(mesh, &gSpectralOpenBCs,  z_discretisation, MESH_ELEVATION_NODE);
  treat_boundarycode(mesh, &gSpectralRigidBCs, u_discretisation, MESH_LAND_NODE);
  
  gSpectralRigidBCs.fluxes=new fluxspec_t[gSpectralRigidBCs.nnodes];
  gSpectralRigidBCs.data  =new dataspec_t[gSpectralRigidBCs.nnodes];
  
  for(i=0;i<gSpectralRigidBCs.nnodes;i++) {
    boundary_normal(gSpectralRigidBCs, mesh, i);
    }
    
  if(tidal_bc) {
    if(gSpectralOpenBCs.nnodes!=0) {
      gSpectralOpenBCs.tides=new tidespec_t[gSpectralOpenBCs.nnodes];
      gSpectralOpenBCs.data=new dataspec_t[gSpectralOpenBCs.nnodes];
      }
    for(k = 0; k < gSpectralOpenBCs.nnodes; k++) {
      gSpectralOpenBCs.tides[k].ztide.a = new float[WaveList.n];
      gSpectralOpenBCs.tides[k].ztide.G = new float[WaveList.n];
      gSpectralOpenBCs.tides[k].ztide.s = 0;
      gSpectralOpenBCs.tides[k].utide.a = new float[WaveList.n];
      gSpectralOpenBCs.tides[k].utide.G = new float[WaveList.n];
      gSpectralOpenBCs.tides[k].utide.s = 0;
      gSpectralOpenBCs.tides[k].vtide.a = new float[WaveList.n];
      gSpectralOpenBCs.tides[k].vtide.G = new float[WaveList.n];
      }
    status= TidalBoundaryConditions( mesh, z_descriptor, gSpectralOpenBCs);
    }

/**----------------------------------------------------------------------------
  initialise open boundary condition array*/
  if(gSpectralOpenBCs.nnodes!=0) {
    cefmo.clamped=new int[gSpectralOpenBCs.nnodes];
    }
  cefmo.nclamped=gSpectralOpenBCs.nnodes;
  for(k = 0; k < gSpectralOpenBCs.nnodes; k++) {
    cefmo.clamped[k]=gSpectralOpenBCs.nodes[k];
    }

/**-----------------------------------------------------------------------
  look for a prior M2 and K1 already computed*/
  string prior=tugo_cfg->tides->PriorSolution.face_value();
  if(prior!="NONE") {
    int frame=9;
    sprintf(filename,"%s/%s.spectral.nc",prior.c_str(), cefmo.dominants[0].name);
    status=SpArchive_Get2D(filename, frame, cefmo, M2state);
    if(status!=0) {
      sprintf(msg,"error by loading 1st dominant wave solution from %s",filename);
      check_error(-1, msg, __LINE__, __FILE__, 1);
      }
    sprintf(filename,"%s/%s.spectral.nc",prior.c_str(), cefmo.dominants[1].name);
    status=SpArchive_Get2D(filename, frame, cefmo, K1state);
    if(status!=0) {
      sprintf(msg,"error by loading 2nd dominant wave solution from %s",filename);
      check_error(-1, msg, __LINE__, __FILE__, 1);
      }
    }

/**----------------------------------------------------------------------------
  iterate dominant wave solver*/
  for(iteration=0;iteration<max_iteration;iteration++) {
    wave_id=M2_id;
    wave=WaveList.waves[wave_id];
/*------------------------------------------------------------------------
    linearized friction coefficients*/
    cefmo.FrictionMatrix =  spectral_friction01(mesh, M2state.u, M2state.v, udim);
    if(cefmo.ndominants==2) spectral_friction02(mesh, cefmo.FrictionMatrix,  K1state.u, K1state.v, udim);
    status=cefmo2D_compute(mesh, data, cefmo, WaveList.waves[wave_id], wave_id, N,rhs,  divergence, iteration);
    for(n = 0; n < udim; n++) {
      M2state.u[n]=cefmo.state.u[n];
      M2state.v[n]=cefmo.state.v[n];
      }
    for(n = 0; n < zdim; n++) {
      M2state.z[n]=cefmo.state.z[n];
      }
/*-----------------------------------------------------------------------------
    Friction will be re-initialised to profit iteration improvments*/
    delete[] cefmo.FrictionMatrix;
    cefmo.FrictionMatrix=0;

    if(cefmo.ndominants==1) continue;

    wave_id=K1_id;
    wave=WaveList.waves[wave_id];
/*------------------------------------------------------------------------
    linearized friction coefficients*/
    cefmo.FrictionMatrix = spectral_friction01(mesh, K1state.u, K1state.v, udim);
    spectral_friction02(mesh, cefmo.FrictionMatrix,  M2state.u, M2state.v, udim);
    status=cefmo2D_compute(mesh, data, cefmo, WaveList.waves[wave_id], wave_id, N, rhs,  divergence, iteration);
    for(n = 0; n < udim; n++) {
      K1state.u[n]=cefmo.state.u[n];
      K1state.v[n]=cefmo.state.v[n];
      }
    for(n = 0; n < zdim; n++) {
      K1state.z[n]=cefmo.state.z[n];
      }
/*-----------------------------------------------------------------------------
    Friction will be re-initialised to profit iteration improvments*/
    delete[] cefmo.FrictionMatrix;
    cefmo.FrictionMatrix=0;
    }
    
/**-----------------------------------------------------------------------
  assume M2 and K1 already computed*/
  if(max_iteration<=0) {
    int frame=-max_iteration-1;
    sprintf(filename,"%s/%s.spectral.nc",gOutputPath,cefmo.dominants[0].name);
    status=SpArchive_Get2D(filename, frame, cefmo, M2state);
    if(status!=0) {
      sprintf(msg,"error by loading 1st dominant wave solution",filename);
      check_error(-1, msg, __LINE__, __FILE__, 1);
      }
    sprintf(filename,"%s/%s.spectral.nc",gOutputPath,cefmo.dominants[1].name);
    status=SpArchive_Get2D(filename, frame, cefmo, K1state);
    if(status!=0) {
      sprintf(msg,"error by loading 2nd dominant wave solution",filename);
      check_error(-1, msg, __LINE__, __FILE__, 1);
      }
    }

/*------------------------------------------------------------------------
  linearized friction coefficients (same for all but dominant tides)*/
  cefmo.FrictionMatrix = spectral_friction02(mesh, M2state.u, M2state.v, udim);
  if(cefmo.ndominants==2) spectral_friction02(mesh, cefmo.FrictionMatrix,  K1state.u, K1state.v, udim);

/**-----------------------------------------------------------------------
  register solution in atlas*/
  atlas_status=atlas->insert(std::pair<const char *,tide2D_t>(WaveList.waves[M2_id].name,M2state));
  if(cefmo.ndominants==2) atlas_status=atlas->insert(std::pair<const char *,tide2D_t>(WaveList.waves[K1_id].name,K1state));

  iteration=0;

/**-----------------------------------------------------------------------
  Solve for astronomical constituents*/
  for(wave_id=0;wave_id<WaveList.n;wave_id++) {
/*------------------------------------------------------------------------
    M2 already done*/
    if(wave_id==M2_id) continue;
/*------------------------------------------------------------------------
    K1 already done*/
    if(wave_id==K1_id) continue;
/*------------------------------------------------------------------------
    do not treat compound tide yet*/
    if(WaveList.waves[wave_id].Ap==0.) continue;

    status=cefmo2D_compute(mesh, data, cefmo, WaveList.waves[wave_id], wave_id, N,rhs,  divergence, iteration);
/**-----------------------------------------------------------------------
    register solution in atlas*/
    tide2D_t local;
    local.z  =new complex <double>[zdim];
    local.u  =new complex <double>[udim];
    local.v  =new complex <double>[udim];
    for(n = 0; n < udim; n++) {
      local.u[n]=cefmo.state.u[n];
      local.v[n]=cefmo.state.v[n];
      }
    for(n = 0; n < zdim; n++) {
      local.z[n]=cefmo.state.z[n];
      }
    atlas_status=atlas->insert(std::pair<const char *,tide2D_t>(WaveList.waves[wave_id].name,local));
//    local.destroy();
    }
    
#if 0
  max_iteration=tugo_cfg->tides->CompoundMaxIteration.numerical_value();

/**-----------------------------------------------------------------------
  parse compound constituents*/
  compound_t *compound[100];
  ncompound=0;
  
  status=initialise_compound(compound, ncompound, WaveList);  
  
  for(k=0;k<ncompound;k++) {
/**-----------------------------------------------------------------------
    register solution in atlas */
    tide2D_t local;
    local.z  =new complex <double>[zdim];
    local.u  =new complex <double>[udim];
    local.v  =new complex <double>[udim];
    for(n = 0; n < udim; n++) {
      local.u[n]=0;
      local.v[n]=0;
      }
    for(n = 0; n < zdim; n++) {
      local.z[n]=0;
      }
    wave_id=WaveList.wave_index(compound[k]->wave.name);
    atlas_status=atlas->insert(std::pair<const char *,tide2D_t>(WaveList.waves[wave_id].name,local));
    status=prediction.add( WaveList.waves[wave_id]);
    status=analysis.add(WaveList.waves[wave_id]);
    }

  for(size_t round=0;round<max_iteration;round++){
  for(size_t kk=0;kk<ncompound;kk++) {
    tidal_wave *generating;
    int ngenerating;
    int *signs;
    tide2D_t state[10];

    divergence=new complex<double>[zdim];
    for(n = 0;n < zdim; n++) {
      divergence[n]=0.0;
      }

/*------------------------------------------------------------------------
    initialize additional non-linear contributions*/
    for(m=0;m<udim;m++) {
      cefmo.action.Fx[m] =0.;
      cefmo.action.Fy[m] =0.;
      }

    wave_id=WaveList.wave_index(compound[kk]->wave.name);
    wave=WaveList.waves[wave_id];
    status=decode_compound(compound[kk]->formula[0], &generating, &signs, &ngenerating);

    status=prediction.add(wave);
    status=analysis.add(wave);
    
    status=analysis.add(wZ0);

    for(i=0;i<ngenerating;i++) {
      status=prediction.add(generating[i]);
      status=analysis.add(generating[i]);
      }
    if(WaveList.wave_index("S2")!=-1) {
      status=prediction.add(wS2);
      status=analysis.add(wS2);
      }
    if(WaveList.wave_index("N2")!=-1) {
      status=prediction.add(wN2);
      status=analysis.add(wN2);
      }
    if(WaveList.wave_index("K2")!=-1) {
      status=prediction.add(wK2);
      status=analysis.add(wK2);
      }
    if(WaveList.wave_index("K1")!=-1) {
      status=prediction.add(wK1);
      status=analysis.add(wK1);
      }
    if(WaveList.wave_index("O1")!=-1) {
      status=prediction.add(wO1);
      status=analysis.add(wO1);
      }
    if(WaveList.wave_index("Q1")!=-1) {
      status=prediction.add(wQ1);
      status=analysis.add(wQ1);
      }

    switch (ngenerating) {
      case 0:
//        check_error(-1, "no generating waves defined for compound tide", __LINE__, __FILE__, 1);
        break;
      case 1:
        state[0]=(*atlas)[generating[0].name];
        state[1]=(*atlas)[generating[0].name];
/**-----------------------------------------------------------------------
        non-linear divergence in continuity equation*/
        status= cefmo.SpDivergence2D_RHS(mesh, state[0].z, state[0].u, state[0].v, zdim, divergence, signs[0]);
/**-----------------------------------------------------------------------
        non-linear advection in momentum equation*/
        status= cefmo.SpAdvection2D(mesh, cefmo.state.h_u, state[0].u, state[0].v, state[0].u, state[0].v, cefmo.action.Fx, cefmo.action.Fy);
        break;
      case 2:
        state[0]=(*atlas)[generating[0].name];
        state[1]=(*atlas)[generating[1].name];
        if(state[0].z==0) break;
        if(state[1].z==0) break;
/**-----------------------------------------------------------------------
        non-linear divergence in continuity equation*/
        status= cefmo.SpDivergence2D_RHS(mesh, state[0].z, state[1].u, state[1].v, zdim, divergence, signs[1]);
        status= cefmo.SpDivergence2D_RHS(mesh, state[1].z, state[0].u, state[0].v, zdim, divergence, signs[1]);
/**-----------------------------------------------------------------------
        non-linear advection in momentum equation*/
        status= cefmo.SpAdvection2D(mesh, cefmo.state.h_u, state[0].u, state[0].v, state[1].u, state[1].v, cefmo.action.Fx, cefmo.action.Fy);
        status= cefmo.SpAdvection2D(mesh, cefmo.state.h_u, state[1].u, state[1].v, state[0].u, state[0].v, cefmo.action.Fx, cefmo.action.Fy);
        break;
      default:
        check_error(-1, "not implemented yet", __LINE__, __FILE__, 1);
        break;
      }

/**-----------------------------------------------------------------------
    non-linear friction in momentum equation*/
    bufx=new complex<double> *[analysis.n];
    bufy=new complex<double> *[analysis.n];
    for(k=0;k<analysis.n;k++) {
      bufx[k]=new complex<double> [udim];
      bufy[k]=new complex<double> [udim];
      for(m=0;m<udim;m++) {
        bufx[k][m]=0;
        bufy[k][m]=0;
        }
      }
/**----------------------------------------------------------------------------
    remove current compound wave to avoid redundancy with linearized friction */
    status=prediction.remove(wave);
/**----------------------------------------------------------------------------
    compute total friction generating forces at all analysis spectrum frequencies */
    cefmo_frictionRHS(prediction, analysis, duration, *atlas, data, udim, bufx, bufy);
/**----------------------------------------------------------------------------
    then put it back for further compound wave computation */
    status=prediction.add(wave);

/**-----------------------------------------------------------------------
    some debugging outputs*/
    if(gArchivesLevel>=CEFMO_ARCHIVES_STANDARD) {
      sprintf(filename, "%s.spectral.nc",wave.name);
      status=archiving_UGdummy2D(filename, mesh, "AdvectionRHSx_a", "AdvectionRHSx_G", "todo", cefmo.action.Fx, mask, round, u_discretisation);
      status=archiving_UGdummy2D(filename, mesh, "AdvectionRHSy_a", "AdvectionRHSy_G", "todo", cefmo.action.Fy, mask, round, u_discretisation);
      status=archiving_UGdummy2D(filename, mesh, "Divergence_a",    "Divergence_G",    "todo", divergence,      mask, round, z_discretisation);
      }
/**----------------------------------------------------------------------------
    get friction generating forces at compound wave frequency*/
    k=analysis.wave_index(wave);
    for(m=0;m<udim;m++) {
      cefmo.action.Fx[m] -=((double) data.Cd[m])*bufx[k][m];
      cefmo.action.Fy[m] -=((double) data.Cd[m])*bufy[k][m];
      }

/**-----------------------------------------------------------------------
    some debugging outputs*/
    if(gArchivesLevel>=CEFMO_ARCHIVES_STANDARD) {
      sprintf(filename, "%s.spectral.nc",wave.name);
      status=archiving_UGdummy2D(filename, mesh, "FrictionRHSx_a", "FrictionRHSx_G", "todo", bufx[k],    mask, round, u_discretisation);
      status=archiving_UGdummy2D(filename, mesh, "FrictionRHSy_a", "FrictionRHSy_G", "todo", bufy[k],    mask, round, u_discretisation);
      }
    for(k=0;k<analysis.n;k++) {
      delete[] bufx[k];
      delete[] bufy[k];
      }
    delete[] bufx;
    delete[] bufy;

    iteration=round;
    status=cefmo2D_compute(mesh, data, cefmo, WaveList.waves[wave_id], wave_id, N, rhs, divergence, iteration);

/**-----------------------------------------------------------------------
    update solution in atlas*/
    tide2D_t local=(*atlas)[WaveList.waves[wave_id].name];
    for(n = 0; n < udim; n++) {
      local.u[n]=cefmo.state.u[n];
      local.v[n]=cefmo.state.v[n];
      }
    for(n = 0; n < zdim; n++) {
      local.z[n]=cefmo.state.z[n];
      }

    delete[] divergence;
    }
  }
  /* STS: NO compound */
#endif

  delete[] rhs;
  delete[] N;

  delete[] cefmo.FrictionMatrix;
  cefmo.FrictionMatrix=0;

  status= SpWE2D_Terminate(mesh, cefmo);

  return(0);
}
