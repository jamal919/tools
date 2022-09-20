
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

\brief tidal friction definitions
*/
/*----------------------------------------------------------------------------*/

#include "tools-define.h"
#include "tides.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void frodom(complex<double> v1,complex<double> v2,double *ro,double *roprim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

/*****************************************************************
 Cette routine permet le calcul des coefficients de frottement
 ro  roprim, pour l'onde dominante M2, traduisant la quasi-
 linearite de la loi de Chezy.
 Ce calcul est fait en determinant les coefficients G00 et G02 a
 partir des integrales elliptiques E et F
**************************************************************** */
{
  double j,j1,m;
  double carj,racj,arac,prac,gml,g00,g02,qj1,vm;
  double epsilon,eps_fj1,unmfj1;
  double vpr,vqr,vpi,vqi,r2pc2;
  double fj1,aa,umv;

  vpr=real(v1);
  vpi=imag(v1);
  vqr=real(v2);
  vqi=imag(v2);

  eps_fj1=1.0E-10;

  r2pc2=vpr*vpr+vpi*vpi+vqr*vqr+vqi*vqi;
  umv=vpi*vqr-vpr*vqi;
  fj1=1.e+00-4.e+00*(umv)*(umv)/(r2pc2*r2pc2);
  aa=sqrt(r2pc2);

  if (fj1 < 0.0) {
    printf("FJ1 : \n",fj1);
    if (fabs(fj1) < eps_fj1) {
      fj1=0;
      printf("FJ1 modifie : \n",fj1);
      }
    }

  j1=sqrt(fj1);

  epsilon=1.0E-06;

/*-----------------------------------------------------------------------
  Calcul des integrales dans le cas ou 0.75<J1<1.*/

  if((j1 > 0.75)&&(j1 < 1.-epsilon)) {
    carj=2.*j1/(1.+j1);
    j=sqrt(carj);
    racj=sqrt(j);
    arac=1.-racj;
    prac=1.+racj;
    m=2.*(arac/prac)*log(2.*prac/arac);
    gml=sqrt(2./(2.-carj))/M_PI;
    g00=(gml/2.)*(prac*prac+((3.+j)*m));
    g02=(gml/3.)*(prac*prac*(2.-carj) -(2.+6.*j+3.*carj+j*carj)*m)/carj;
    }

/*-----------------------------------------------------------------------
  Calcul des integrales dans le cas ou 0.<J1<=0.75 */
  else if((j1 > 0.)&&(j1 <= 0.75)) {
    qj1=fj1*fj1;
    g00=1.-0.0625*fj1-0.01465*qj1;
    g02=(0.5+0.046875*fj1+0.0170898*qj1)*j1;
    }
/*-----------------------------------------------------------------------
  Calcul de ro et roprim pour 0.0<J1<1.*/

  if((j1 > 0.)&&(j1 < 1.-epsilon)) {
      vm=aa/sqrt(2.);
      g02=g02/(2.*j1);
      *ro=vm*(g00+g02);
      unmfj1=1-fj1;
      if(unmfj1<0.0) {
        printf("1-FJ1 : \n",unmfj1);
        if(fabs(unmfj1)<eps_fj1) {
          unmfj1=0.;
          printf("1-FJ1 modifie : \n",unmfj1);
          }
        }
      *roprim=vm*g02*sqrt(1.-fj1);
    }
  else
    {
/*-----------------------------------------------------------------------
    Calcul de ro et roprim pour J1=0. */
    if(j1==0.) {
      *ro=aa*0.88388;
      *roprim=aa*0.176777;
      }
/*-----------------------------------------------------------------------
    Calcul de ro et roprim pour J1=1.  ro=8/(3*M_PI)*aa */
    if(j1>=1.-epsilon) {
      *ro=aa*(8/(3*M_PI));
      *roprim=0.;
      }
    }
/*-----------------------------------------------------------------------
  Test sur umv pour le signe de roprim */
  if(umv<0.)*roprim=-*roprim;

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void frosec(complex<double> v1,complex<double> v2,double *r,double *rprim,double *rsecon)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

/****************************************************************
! Cette routine permet le calcul des coefficients de frottement
! R, RPRIM et RSECON, pour une onde secondaire du spectre
! Ce calcul est fait en determinant les coefficients G00 et G02 a
! partir des integrales elliptiques E et F
!****************************************************************/
{
double  j,j1,m;
double  vpr,vqr,vpi,vqi;
double  coef,r2,c2,r2mc2,r2pc2,aa,umv,upv,fj1;
double  carj,racj,arac,prac,prac2,gml,g00,g02j1,dmj2,qj1;
double  epsilon;

coef=3.0/(2.*sqrt(2.));

vpr=real(v1)+0.001; /// HERE !!!
vpi=imag(v1);
vqr=real(v2);
vqi=imag(v2);

r2=vpr*vpr+vpi*vpi;
c2=vqr*vqr+vqi*vqi;

r2pc2=r2+c2;
r2mc2=r2-c2;
aa=sqrt(r2pc2);
umv=vpr*vqi-vpi*vqr;
upv=vpr*vqr+vpi*vqi;
fj1=1.-4.*umv*umv/(r2pc2*r2pc2);
j1=sqrt(fj1);

epsilon=1.0E-06;

//----------------- Calcul de G00 et G02/j1 pour 0.75<J1<1. -------------

if((j1 > 0.75) && (j1 < 1.-epsilon)) {
  j=sqrt(2.*j1/(1.+j1));
  carj=j*j;
  dmj2=2.-carj;
  racj=sqrt(j);
  arac=1.-racj;
  prac=1.+racj;
  prac2=prac*prac;
  m=2.*(arac/prac)*log(2.*prac/arac);
  gml=sqrt(1./(2*dmj2))/M_PI;
  g00=gml*(prac2+(3.+j)*m);
  gml=sqrt(2.*dmj2)/(3.*M_PI*carj*carj);
  g02j1=gml*(prac2*dmj2-(2.+6.*j+3.*carj+j*carj)*m);
  }
  
//-------------- Calcul de G00 et G02/j1 pour 0.<J1<=0.75 ---------------
  
else if((j1 > 0.) && (j1 <= 0.75)) {
  qj1=fj1*fj1;
  g00=1.e+00-0.0625E+00 *fj1-0.01465E+00 *qj1;
  g02j1=0.5E+00 +0.046875E+00 *fj1+0.0170898E+00 *qj1;
  }

//----------------- Calcul de G00 et G02/j1 pour J1=0. ------------------
  
else if(j1 == 0.) {
  g00=1.e+00;
  g02j1=0.5E+00;
  }
  
//----------------- Calcul de G00 et G02/j1 pour J1=1. ------------------
  else if (j1 >= 1.-epsilon) {
  g00=0.9003163E+00;
  g02j1=0.6002109E+00;
  }

/*-----------------------------------------------------------------------
! Calculs des coefficients de frottements R,R' et R" + termes correctifs
!-----------------------------------------------------------------------*/

*r     =g00*aa*coef;
*rprim =g02j1*0.5*r2mc2*coef/aa;
*rsecon=g02j1*upv*coef/aa;

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void spectral_friction01(zmatrix2x2_t *fric,complex<double> u,complex<double> v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/****************************************************************
 Cette routine traduit le processus iteratif, decrit dans le
 #TidesHarmonicConstituentsModeling_REF, et mis en oeuvre pour le
 calcul de l'onde dominante. Le calcul de l'acceleration est fait
 sur les coefficients de frottement ro et roprim
*****************************************************************/
  double R1,R2;
  const complex<double> j=dcomplex(0.,1.);
  
  frodom(u,v,&R1,&R2);
  
  fric->c[0][0]=    R1;//r
  fric->c[1][0]=  j*R2;//r'
  fric->c[0][1]= -j*R2;//r''
  fric->c[1][1]=    R1;//r'''
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void spectral_friction02(zmatrix2x2_t *fric,complex<double> u,complex<double> v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/****************************************************************
 Cette routine traduit le processus iteratif, decrit dans le
 #TidesHarmonicConstituentsModeling_REF, et mis en oeuvre pour le
 calcul des ondes secondaires.
*****************************************************************/
  double R1,R2,R3;
  
  frosec(u,v,&R1,&R2,&R3);

  fric->c[0][0]=(R1+R2);
  fric->c[1][0]=R3;
  fric->c[0][1]=R3;
  fric->c[1][1]=(R1-R2);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

zmatrix2x2_t *spectral_friction01(mesh_t & mesh, complex <double> *u, complex <double> *v, int nnodes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/****************************************************************
 Cette routine traduit le processus iteratif, decrit dans le
 J.P.O (LE PROVOST, ROUGIER, PONCET), et mis en oeuvre pour le
 calcul de l'onde dominante. Le calcul de l'acceleration est fait
 sur les coefficients de frottement ro et roprim
*****************************************************************/
  double c=1.0,R1,R2;
  double a,G;
  int n,status;
  int udim,zdim;
  complex<double> j=complex<double>(0.,1.);
  zmatrix2x2_t *FrictionMatrix;
  
  extern void frodom(complex<double> v1,complex<double> v2,double *ro,double *roprim);

//  status=paire_dimension(mesh,&zdim,&udim);

  udim=nnodes;

  FrictionMatrix = new zmatrix2x2_t[udim];

  for(n=0;n<udim;n++) {
    frodom(u[n],v[n],&R1,&R2);
    
//     if(isnormal(R1)==0) {
//       printf("troubles\n");
//       }
//     if(isnormal(R2)==0) {
//       printf("troubles\n");
//       }
    
//     switch (counter%3) {
//       case 0:
//         RO[0][n]=R1;
//         ROPRIM[0][n]=R2;
//         break;
//       case 1:
//         RO[1][n]=R1;
//         ROPRIM[1][n]=R2;
//         break;
//       case 2:
//         accelere(RO[0][n],RO[1][n],&R1);
//         accelere(ROPRIM[0][n],ROPRIM[1][n],&R2);
//         break;
//       }
/**----------------------------------------------------------------------
    c[j][i]=matrix's coefficent row i, column j

    FRT u = r0 r1 u
    FRT v = r2 r3 v
----------------------------------------------------------------------*/
    FrictionMatrix[n].c[0][0]=    R1;
    FrictionMatrix[n].c[1][0]=  j*R2;
    FrictionMatrix[n].c[0][1]= -j*R2;
    FrictionMatrix[n].c[1][1]=    R1;
    }

  return(FrictionMatrix);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

zmatrix2x2_t *spectral_friction02(mesh_t & mesh, complex <double> *u, complex <double> *v, int nnodes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/****************************************************************
 Cette routine traduit le processus iteratif, decrit dans le
 J.P.O (LE PROVOST, ROUGIER, PONCET), et mis en oeuvre pour le
 calcul de l'onde secondaire.
*****************************************************************/
  double c=1.0,R1,R2,R3;
  double a,G;
  int n;
  zmatrix2x2_t *FrictionMatrix;

  FrictionMatrix=new zmatrix2x2_t[nnodes];

  for(n=0;n<nnodes;n++) {
    frosec(u[n],v[n],&R1,&R2,&R3);
    R1*=c;
    R2*=c;
    R3*=c;
/*
  r0(k)= (rdamp(3,l)+rdamp(4,l))*c/h
  r1(k)= rdamp(5,l)*c/h
  r2(k)= r1(k)
  r3(k)= (rdamp(3,l)-rdamp(4,l))*c/h
*/
    FrictionMatrix[n].c[0][0]=(R1+R2);
    FrictionMatrix[n].c[1][0]=R3;
    FrictionMatrix[n].c[0][1]=R3;
    FrictionMatrix[n].c[1][1]=(R1-R2);
    }

  return(FrictionMatrix);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void spectral_friction02(mesh_t mesh, zmatrix2x2_t *FrictionMatrix, complex <double> *u, complex <double> *v, int nnodes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double c=1.0,R1,R2,R3;
  double a,G;
  int n;

  for(n=0;n<nnodes;n++) {
    frosec(u[n],v[n],&R1,&R2,&R3);
    R1*=c;
    R2*=c;
    R3*=c;
    FrictionMatrix[n].c[0][0]+=(R1+R2);
    FrictionMatrix[n].c[1][0]+=R3;
    FrictionMatrix[n].c[0][1]+=R3;
    FrictionMatrix[n].c[1][1]+=(R1-R2);
    }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  matrix2x2_t *SpectralDrag_init(discretisation_t descriptor, tidal_wave wave, double *N, double *h, parameter_t data, double wdHmin, double wdHmax, int nndes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------

  Internal wave drag, CEFMO type of

 -----------------------------------------------------------------------*/
{
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

  18/09/2010: dependency upon wave currents still to be done
//   u0(k)=SQRT(ABS(vx_g(k,itri))**2+ABS(vy_g(k,itri))**2)
//   lamda(k)=(2./3.14)*u0(k)/freq
//   tau=1.-f(k)*f(k)/freq2


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  int n;
  vector2D_t grd;
  double u, v, w, kappa, lambda, omega2;
  double factor, lamda, hramp;
  double c, s, r, S,f, k_b, k_0, alpha, tau;
  double t,p;
//  double hmin = 300, hmax = 500;
  matrix2x2_t *matrix;

  alpha = 2.5;
/**----------------------------------------------------------------------------
  Slope change at 50 km, k= 2 * pi/50000 */
  k_b = 6.28 / 50000.;
  k_0 = 6.28 / 9000.;

  tau=2./M_PI/wave.omega;
  omega2=wave.omega*wave.omega;


  matrix = new matrix2x2_t[nndes];

  w=wave.omega*dph2rps;
//  lamda=150000.*sqrt(w/two_Omega);
  lamda=150000;
//  lamda=150000.*29./wave.omega;
  
  for(n = 0; n < nndes; n++) {
    matrix[n].c[0][0] = 0.0;
    matrix[n].c[0][1] = 0.0;
    matrix[n].c[1][0] = 0.0;
    matrix[n].c[1][1] = 0.0;
    S   = (double) descriptor.nodes[n].s;
//     t=descriptor.nodes[n].lon*180./M_PI;
//     p=descriptor.nodes[n].lat*180./M_PI;
//     t=geo_recale(t,-180.,180.);
//     if( (t>-70.) && (t<-55.) && (p>-60.) && (p<-50.) ) { /// HERE
//       data.Nbar[n]=5.e-04;
//       }
    f = S * two_Omega;
    c=1.-f*f/(w*w);
    
    c=max(0.2,c);
    
    if(c<0.) continue;
    c=sqrt(c);
    
    if(h[n] > wdHmin) {
/*----------------------------------------------------------------------
      linear ramp [0:1] between hmin and hmax */
      hramp = min(1., (h[n] - wdHmin) / (wdHmax - wdHmin));
/*----------------------------------------------------------------------
      sinusoidale ramp between hmin and hmax */
      hramp = 0.5 * (1. - cos(hramp * 3.14));
      
//      double cz=min(1., 500./h[n]);
      
      grd.x = data.dhdx[n];
      grd.y = data.dhdy[n];
//       lamda=data.celerity[n]*2.*M_PI/w;
//       lamda=150000.;
      kappa = 2 * M_PI / lamda;
      factor = c * hramp * data.wdc1[n] * N[n] / kappa;
//      factor = c * hramp * data.wdc1[n] * data.Nbar[n] / kappa;
/*----------------------------------------------------------------------
      divide by h to achieve average force */
      if(h[n]==0) {
        printf("SpectralDrag_init : trouble\n");
        }
      matrix[n].c[0][0] += -factor * grd.x * grd.x / h[n];
      matrix[n].c[0][1] += -factor * grd.x * grd.y / h[n];
      matrix[n].c[1][0] += -factor * grd.x * grd.y / h[n];
      matrix[n].c[1][1] += -factor * grd.y * grd.y / h[n];
      }
    }

  return (matrix);
}
