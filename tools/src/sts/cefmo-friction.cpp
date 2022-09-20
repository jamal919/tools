
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

/**----------------------------------------------------------------------------

  Tidal spectral solver (cefmo type)

  solves the quasi-linearized complex equations

  Available discretisation : LGP0xLGP1 and DLGP1xLGP2

  Mostly identical to CEFMO, except for:

    -numerical integration (exact here, Gauss quadrature in CEFMO)

    -bathymetry discretisation (computational nodes here, Gauss points in CEFMO)

    -M2 and K1 velocity mix (for friction coefficient computation) not implemented yet

-----------------------------------------------------------------------------*/
#include "tugo-prototypes.h"
#include <string>
#include "dry.h"
#include "ice.h"
#include "tmatrix.h"
#include "cefmo.h"
#include "tides.h"
#include "tides.def"
#include "geo.h"
#include "keywords.hpp"
#include "poc-netcdf.hpp"


extern int fft (const char *, double *x, double dt, int n);

/**----------------------------------------------------------------------------
for dominant wave iteration acceleration*/
int counter=0;

extern atlas2D_t cefmo_atlas2D;
extern atlas3D_t cefmo_atlas3D;

int fft_count=0;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int initialise_compound(compound_t **compound, int & ncompound, spectrum_t WaveList)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int wave_id, status;
  
/**-----------------------------------------------------------------------
  parse compound constituents*/
  
  ncompound=0;
  
  for(wave_id=0;wave_id<WaveList.n;wave_id++) {
//    spectrum_t generating;
    tidal_wave *generating;
    int ngenerating;
    int *signs;
    if(wave_id==WaveList.wave_index("M4")) {
      compound[ncompound]=new compound_t;
      compound[ncompound]->wave=wM4;
      compound[ncompound]->id=wave_id;
      compound[ncompound]->formula.push_back((string) "M2+M2");
//      compound[ncompound]->formula.push_back((string) "M6-M2");
      ncompound++;
      }
//     if(wave_id==WaveList.wave_index("M2")) {
//       compound[ncompound]=new compound_t;
//       compound[ncompound]->wave=wM2;
//       compound[ncompound]->formula.push_back((string) "M4-M2");
//       ncompound++;
//       }
    if(wave_id==WaveList.wave_index("M6")) {
      compound[ncompound]=new compound_t;
      compound[ncompound]->wave=wM6;
      compound[ncompound]->id=wave_id;
      compound[ncompound]->formula.push_back((string) "M2+M4");
      ncompound++;
      }
    if(wave_id==WaveList.wave_index("M8")) {
      compound[ncompound]=new compound_t;
      compound[ncompound]->wave=wM8;
      compound[ncompound]->id=wave_id;
//      compound[ncompound]->formula.push_back((string) "M4+M4");
      compound[ncompound]->formula.push_back((string) "M2+M6");
      ncompound++;
      }
    if(wave_id==WaveList.wave_index("M3")) {
      compound[ncompound]=new compound_t;
      compound[ncompound]->wave=wM3;
      compound[ncompound]->formula.push_back((string) "");
      compound[ncompound]->id=wave_id;
      status=decode_compound(compound[ncompound]->formula[0], &generating, &signs, &ngenerating);
      ncompound++;
      }
    if(wave_id==WaveList.wave_index("MS4")) {
      compound[ncompound]=new compound_t;
      compound[ncompound]->wave=wMS4;
      compound[ncompound]->id=wave_id;
      compound[ncompound]->formula.push_back((string) "M2+S2");
      ncompound++;
      }
    if(wave_id==WaveList.wave_index("MN4")) {
      compound[ncompound]=new compound_t;
      compound[ncompound]->wave=wMN4;
      compound[ncompound]->id=wave_id;
      compound[ncompound]->formula.push_back((string) "M2+N2");
      ncompound++;
      }
    if(wave_id==WaveList.wave_index("S4")) {
      compound[ncompound]=new compound_t;
      compound[ncompound]->wave=wS4;
      compound[ncompound]->id=wave_id;
      compound[ncompound]->formula.push_back((string) "S2+S2");
      ncompound++;
      }
    if(wave_id==WaveList.wave_index("N4")) {
      compound[ncompound]=new compound_t;
      compound[ncompound]->wave=wN4;
      compound[ncompound]->id=wave_id;
      compound[ncompound]->formula.push_back((string) "N2+N2");
      ncompound++;
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int decode_compound(string s, spectrum_t *generating, int **signs)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------------
  */
  size_t first, pos,plus,minus;
  int count;
  date_t origin;
  spectrum_t reference=initialize_tide(&gAstronomicAngles,origin);
  int k;

  plus =s.find("+");
  minus=s.find("-");
  pos=min(plus,minus);

  count=0;
  while(pos!=string::npos) {
    count++;
    plus =s.find("+",pos+1);
    minus=s.find("-",pos+1);
    pos=min(plus,minus);
    }

  *generating=spectrum_t(count+1);

  plus =s.find("+");
  minus=s.find("-");
  first=0;
  pos=min(plus,minus);

  count=0;
  while(pos!=string::npos) {
    string name=s.substr(first,pos-first);
    k=generating->add(reference,(const char *) name.c_str(),0);
    first=pos+1;
    count++;
    plus =s.find("+",pos+1);
    minus=s.find("-",pos+1);
    pos=min(plus,minus);
    }
  string name=s.substr(first,s.size());
  k=generating->add(reference,(const char *) name.c_str(),0);
  
  reference.destroy();
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int decode_compound(string s, tidal_wave **generating, int **signs, int *ngenerating)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------------
  */
  size_t first, pos,plus,minus;
  int count;
  date_t origin;
  spectrum_t reference=initialize_tide(&gAstronomicAngles,origin);
  int k;

  plus =s.find("+");
  minus=s.find("-");
  pos=min(plus,minus);
  
  if(pos==string::npos) {
    *ngenerating=0;
    return(0);
    }

  count=0;
  while(pos!=string::npos) {
    count++;
    plus =s.find("+",pos+1);
    minus=s.find("-",pos+1);
    pos=min(plus,minus);
    }

  *ngenerating=count+1;
  *generating=new tidal_wave[count+1];

  *signs=new int[count+1];
  (*signs)[0]=1;

  plus =s.find("+");
  minus=s.find("-");
  first=0;
  pos=min(plus,minus);

  count=0;
  while(pos!=string::npos) {
    if(plus<minus) (*signs)[count+1]=1;
    else  (*signs)[count+1]=-1;
    string name=s.substr(first,pos-first);
    (*generating)[count]=reference.wave((const char *) name.c_str());
    first=pos+1;
    count++;
    plus =s.find("+",pos+1);
    minus=s.find("-",pos+1);
    pos=min(plus,minus);
    }
  string name=s.substr(first,s.size());
  (*generating)[count]=reference.wave((const char *) name.c_str());
  
  switch (*ngenerating) {
    case 2:
      if(strcmp(((*generating)[0]).name,((*generating)[1].name))==0) (*ngenerating)--;
      break;
    }
  reference.destroy();
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int friction_ratio(mesh_t & mesh, discretisation_t & u_descriptor, float *Cd, double *ratio)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  using namespace Keywords;
  int k, m, n, status;
  double  *ratioLGP1,c,z;
  
  ratioLGP1=new double[mesh.nvtxs];

  status=fe_projection(mesh, ratio, u_descriptor.type, ratioLGP1, gP1discretisation);

  for(int depth=0;depth<50;depth++) {
    for(n = 0; n < mesh.nvtxs; n++) {
      c=z=0.0;
      for(k=0;k<mesh.vertices[n].nngh;k++) {
        m=mesh.vertices[n].ngh[k];
        z+=ratioLGP1[m];
        c++;
        }
      ratioLGP1[n]=z/c;
      }
    }
   
  status=fe_projection(mesh, ratioLGP1, gP1discretisation, ratio, u_descriptor.type);

  delete[] ratioLGP1;
  
  printf("smoothing done\n");
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void friction_coefficient(zmatrix2x2_t FrictionMatrix, double Cd, double r, double u0, double h, double ratio, complex<double> *r1, complex<double> *r2, complex<double> *r3, complex<double> *r4)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  using namespace Keywords;
  int status;
  
  double R=r/h;
//  double R=0.8*r/h+gR_default/50;
  
  double alpha=0.0;
//  alpha=-M_PI/6.*min(h/200.,1.);
//  alpha=-M_PI/6.;

  double cs=cos(alpha),sn=sin(alpha);
  complex< double > shift=complex< double >(cs,sn);
  
  switch(gBottomFrictionType) {
    case FRICTION_LINEAR:
      *r1=r/h;
      *r2=0.;
      *r3=0.;
      *r4=r/h;
//       *r1=  cs*R;
//       *r2= -sn*R;
//       *r3=  sn*R;
//       *r4=  cs*R;
//       *r1=gR_default/h;
//       *r2=0.;
//       *r3=0.;
//       *r4=gR_default/h;
      break;
      
    case FRICTION_QUADRATIC:
/**----------------------------------------------------------------------
      Quadratic friction with Cd formally prescribed*/
    case FRICTION_KARMAN:
/**----------------------------------------------------------------------
      Quadratic friction with Cd deduced from log velocity profile assumption*/
//      u0=gU0;
//      u0=min(u0*h/200,gU0);
      
      *r1=(1.-ratio)*(gR_default/h)+ratio*Cd*(FrictionMatrix.c[0][0]+u0)/h;
      *r2=ratio*Cd*(FrictionMatrix.c[1][0])/h;
      *r3=ratio*Cd*(FrictionMatrix.c[0][1])/h;
      *r4=(1.-ratio)*(gR_default/h)+ratio*Cd*(FrictionMatrix.c[1][1]+u0)/h;
      
//       *r1*=shift;
//       *r2*=shift;
//       *r3*=shift;
//       *r4*=shift;

//       *r1=(1.-ratio)*(gR_default/h)+ratio*Cd*(FrictionMatrix.c[0][0])/h;
//       *r2=ratio*Cd*(FrictionMatrix.c[1][0])/h;
//       *r3=ratio*Cd*(FrictionMatrix.c[0][1])/h;
//       *r4=(1.-ratio)*(gR_default/h)+ratio*Cd*(FrictionMatrix.c[1][1])/h;
//
//       *r1*=shift;
//       *r2*=shift;
//       *r3*=shift;
//       *r4*=shift;
//
//       *r1+=ratio*Cd*(u0)/h;
//       *r4+=ratio*Cd*(u0)/h;
      
      

//       u0=gU0*(1.-tanh(h/50.));
//       *r1=(1.-ratio)*(gR_default/h)+ratio*Cd*(FrictionMatrix.c[0][0]+u0)/h;
//       *r2=ratio*Cd*(FrictionMatrix.c[1][0]+u0)/h;
//       *r3=ratio*Cd*(FrictionMatrix.c[0][1]+u0)/h;
//       *r4=(1.-ratio)*(gR_default/h)+ratio*Cd*(FrictionMatrix.c[1][1]+u0)/h;

//       *r1=Cd*(FrictionMatrix.c[0][0]+gU0)/h; /// HERE!!!
//       *r2=Cd*FrictionMatrix.c[1][0]/h;
//       *r3=Cd*FrictionMatrix.c[0][1]/h;
//       *r4=Cd*(FrictionMatrix.c[1][1]+gU0)/h; /// HERE!!!
      break;
    default:
      check_error(-1, "illicit friction mode", __LINE__, __FILE__, 1);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void friction_coefficient(zmatrix2x2_t FrictionMatrix, double Cd, double r, double h, complex<double> *r1, complex<double> *r2, complex<double> *r3, complex<double> *r4)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double ratio=1.0,u0=gU0;
  friction_coefficient(FrictionMatrix, Cd, r, u0, h, ratio, r1, r2, r3, r4);
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// zmatrix2x2_t *spectral_friction01(mesh_t & mesh, complex <double> *u, complex <double> *v, int nnodes)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
// /****************************************************************
//  Cette routine traduit le processus iteratif, decrit dans le
//  J.P.O (LE PROVOST, ROUGIER, PONCET), et mis en oeuvre pour le
//  calcul de l'onde dominante. Le calcul de l'acceleration est fait
//  sur les coefficients de frottement ro et roprim
// *****************************************************************/
//   double c=1.0,R1,R2;
//   double a,G;
//   int n,status;
//   int udim,zdim;
//   complex<double> j=complex<double>(0.,1.);
//   zmatrix2x2_t *FrictionMatrix;
//   
//   extern void frodom(complex<double> v1,complex<double> v2,double *ro,double *roprim);
// 
// //  status=paire_dimension(mesh,&zdim,&udim);
// 
//   udim=nnodes;
// 
//   FrictionMatrix = new zmatrix2x2_t[udim];
// 
//   for(n=0;n<udim;n++) {
//     frodom(u[n],v[n],&R1,&R2);
//     
// //     if(isnormal(R1)==0) {
// //       printf("troubles\n");
// //       }
// //     if(isnormal(R2)==0) {
// //       printf("troubles\n");
// //       }
//     
// //     switch (counter%3) {
// //       case 0:
// //         RO[0][n]=R1;
// //         ROPRIM[0][n]=R2;
// //         break;
// //       case 1:
// //         RO[1][n]=R1;
// //         ROPRIM[1][n]=R2;
// //         break;
// //       case 2:
// //         accelere(RO[0][n],RO[1][n],&R1);
// //         accelere(ROPRIM[0][n],ROPRIM[1][n],&R2);
// //         break;
// //       }
// /**----------------------------------------------------------------------
//     c[j][i]=matrix's coefficent row i, column j
// 
//     FRT u = r0 r1 u
//     FRT v = r2 r3 v
// ----------------------------------------------------------------------*/
//     FrictionMatrix[n].c[0][0]=    R1;
//     FrictionMatrix[n].c[1][0]=  j*R2;
//     FrictionMatrix[n].c[0][1]= -j*R2;
//     FrictionMatrix[n].c[1][1]=    R1;
//     }
// 
//   counter++;
//   return(FrictionMatrix);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// zmatrix2x2_t *spectral_friction02(mesh_t & mesh, complex <double> *u, complex <double> *v, int nnodes)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
// /****************************************************************
//  Cette routine traduit le processus iteratif, decrit dans le
//  J.P.O (LE PROVOST, ROUGIER, PONCET), et mis en oeuvre pour le
//  calcul de l'onde secondaire.
// *****************************************************************/
//   double c=1.0,R1,R2,R3;
//   double a,G;
//   int n;
//   zmatrix2x2_t *FrictionMatrix;
// 
//   FrictionMatrix=new zmatrix2x2_t[nnodes];
// 
//   for(n=0;n<nnodes;n++) {
//     frosec(u[n],v[n],&R1,&R2,&R3);
//     R1*=c;
//     R2*=c;
//     R3*=c;
// /*
//   r0(k)= (rdamp(3,l)+rdamp(4,l))*c/h
//   r1(k)= rdamp(5,l)*c/h
//   r2(k)= r1(k)
//   r3(k)= (rdamp(3,l)-rdamp(4,l))*c/h
// */
//     FrictionMatrix[n].c[0][0]=(R1+R2);
//     FrictionMatrix[n].c[1][0]=R3;
//     FrictionMatrix[n].c[0][1]=R3;
//     FrictionMatrix[n].c[1][1]=(R1-R2);
//     }
// 
//   return(FrictionMatrix);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   void spectral_friction02(mesh_t mesh, zmatrix2x2_t *FrictionMatrix, complex <double> *u, complex <double> *v, int nnodes)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   double c=1.0,R1,R2,R3;
//   double a,G;
//   int n;
// 
//   for(n=0;n<nnodes;n++) {
//     frosec(u[n],v[n],&R1,&R2,&R3);
//     R1*=c;
//     R2*=c;
//     R3*=c;
//     FrictionMatrix[n].c[0][0]+=(R1+R2);
//     FrictionMatrix[n].c[1][0]+=R3;
//     FrictionMatrix[n].c[0][1]+=R3;
//     FrictionMatrix[n].c[1][1]+=(R1-R2);
//     }
// }
// 

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int SpectralFriction_init(discretisation_t descriptor, tidal_wave wave,parameter_t data, int nndes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,frame=0,status;
  using namespace Keywords;
  char filename[1024];
  double *buffer;
  
  switch(gBottomFrictionType) {
    case FRICTION_LINEAR:
      for(n=0;n<nndes;n++) data.rlinear[n]=gR_default; /// HERE
//       sprintf(filename,"%s/%s.spectral.nc",gOutputPath,"M2");
//       buffer=new double[nndes];
//       status=poc_get_UG3D(filename, frame, "Vmodule", buffer);
//       for(n=0;n<nndes;n++) data.r[n]=buffer[n]*2.0e-03;
//       delete[] buffer;
      break;
      
    case FRICTION_QUADRATIC:
    case FRICTION_KARMAN:
      break;
    default:
      check_error(-1, "illicit friction mode", __LINE__, __FILE__, 1);
      break;
    }
  return(0);
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   matrix2x2_t *SpectralDrag_init(discretisation_t descriptor, tidal_wave wave, double *N, double *h, parameter_t data, int nndes)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /**-----------------------------------------------------------------------
// 
//   Internal wave drag, CEFMO type of
// 
//  -----------------------------------------------------------------------*/
// {
// /**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// 
//   Development notes:
// 
//   Check :
// 
//   Note:
// 
//   18/09/2010: dependency upon wave currents still to be done
// //   u0(k)=SQRT(ABS(vx_g(k,itri))**2+ABS(vy_g(k,itri))**2)
// //   lamda(k)=(2./3.14)*u0(k)/freq
// //   tau=1.-f(k)*f(k)/freq2
// 
// 
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
//   int n;
//   vector2D_t grd;
//   double u, v, w, kappa, lambda, omega2;
//   double factor, lamda, hramp;
//   double c, s, r, S,f, k_b, k_0, alpha, tau;
//   double t,p;
// //  double hmin = 300, hmax = 500;
//   matrix2x2_t *matrix;
// 
//   alpha = 2.5;
// /**----------------------------------------------------------------------------
//   Slope change at 50 km, k= 2 * pi/50000 */
//   k_b = 6.28 / 50000.;
//   k_0 = 6.28 / 9000.;
// 
//   tau=2./M_PI/wave.omega;
//   omega2=wave.omega*wave.omega;
// 
// 
//   matrix = new matrix2x2_t[nndes];
// 
//   w=wave.omega*dph2rps;
// //  lamda=150000.*sqrt(w/two_Omega);
//   lamda=150000;
// //  lamda=150000.*29./wave.omega;
//   
//   for(n = 0; n < nndes; n++) {
//     matrix[n].c[0][0] = 0.0;
//     matrix[n].c[0][1] = 0.0;
//     matrix[n].c[1][0] = 0.0;
//     matrix[n].c[1][1] = 0.0;
//     S   = (double) descriptor.nodes[n].s;
// //     t=descriptor.nodes[n].lon*180./M_PI;
// //     p=descriptor.nodes[n].lat*180./M_PI;
// //     t=geo_recale(t,-180.,180.);
// //     if( (t>-70.) && (t<-55.) && (p>-60.) && (p<-50.) ) { /// HERE
// //       data.Nbar[n]=5.e-04;
// //       }
//     f = S * two_Omega;
//     c=1.-f*f/(w*w);
//     
//     c=max(0.2,c);
//     
//     if(c<0.) continue;
//     c=sqrt(c);
//     
//     if(h[n] > wdHmin) {
// /*----------------------------------------------------------------------
//       linear ramp [0:1] between hmin and hmax */
//       hramp = min(1., (h[n] - wdHmin) / (wdHmax - wdHmin));
// /*----------------------------------------------------------------------
//       sinusoidale ramp between hmin and hmax */
//       hramp = 0.5 * (1. - cos(hramp * 3.14));
//       
// //      double cz=min(1., 500./h[n]);
//       
//       grd.x = data.dhdx[n];
//       grd.y = data.dhdy[n];
// //       lamda=data.celerity[n]*2.*M_PI/w;
// //       lamda=150000.;
//       kappa = 2 * M_PI / lamda;
//       factor = c * hramp * data.wdc1[n] * N[n] / kappa;
// //      factor = c * hramp * data.wdc1[n] * data.Nbar[n] / kappa;
// /*----------------------------------------------------------------------
//       divide by h to achieve average force */
//       if(h[n]==0) {
//         printf("SpectralDrag_init : trouble\n");
//         }
//       matrix[n].c[0][0] += -factor * grd.x * grd.x / h[n];
//       matrix[n].c[0][1] += -factor * grd.x * grd.y / h[n];
//       matrix[n].c[1][0] += -factor * grd.x * grd.y / h[n];
//       matrix[n].c[1][1] += -factor * grd.y * grd.y / h[n];
//       }
//     }
// 
//   return (matrix);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void cefmo_frictionRHS_SEQ(spectrum_t prediction, spectrum_t analysis, double duration,
                       atlas2D_t atlas, parameter_t data, int nndes, complex<double> **bufx, complex<double> **bufy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**-----------------------------------------------------------------------
  Direct harmonic analysis of friction for non-linear RHS
-------------------------------------------------------------------------*/
  double c;
  double *residuals,tauxm,tauym;
  double *taux,*tauy;
  int k,m,n,nframe,w,status;
  int maxstep=1;
  double *u,*v,*ax,*Gx,*ay,*Gy;
  harmonic_t h1,h2;
  date_t start;
  double dt,*time;
  int nodal=0;
  char output[64];

  start.year=1990;
  start.month=1;
  start.day=1;
  start.second=0;

  dt=1800.0;
  status=harmonic_check(analysis, duration, dt);
  if(status!=0) {
    check_error(-1, "duration not suitable", __LINE__, __FILE__, 1);
    }

  nframe=(int) NINT(duration/dt);

  time = new double[nframe];
  u    = new double[nframe];
  v    = new double[nframe];
  taux = new double[nframe];
  tauy = new double[nframe];
  residuals = new double[nframe];

  ax= new double [analysis.n];
  Gx= new double [analysis.n];
  ay= new double [analysis.n];
  Gy= new double [analysis.n];

  for(k=0;k<nframe;k++) time[k]=(double)k*dt;

//  init_argument(start);
  status=harmonic_init01(nframe, start, time, analysis, &h1,0,gAstronomicAngles);// tools update
  status=harmonic_coefficients(nframe, start, time, prediction, &h2,0);

  for(n=0;n<nndes;n++) {
//     double *u,*v,*H,*ax,*Gx,*ay,*Gy;
//     double *taux,*tauy;
//     u    = new double[nframe];
//     v    = new double[nframe];
//     H    = new double[nframe];
//     taux = new double[nframe];
//     tauy = new double[nframe];
//     ax= new double [analysis.n];
//     Gx= new double [analysis.n];
//     ay= new double [analysis.n];
//     Gy= new double [analysis.n];
/**----------------------------------------------------------------------------
    create time serie (tidal prediction) */
    for(m=0;m<nframe;m++) {
      u[m]=0;
      v[m]=0;
      }
    for(w = 0; w < prediction.n; w++) {
      double cs[3],sn[3];
      double *wcs, *wsn;
      tide2D_t state=atlas[prediction.waves[w].name];
      cs[0]=  state.u[n].real();
      cs[1]=  state.v[n].real();
      cs[2]=  state.z[n].real();
      sn[0]= -state.u[n].imag();
      sn[1]= -state.v[n].imag();
      sn[2]= -state.z[n].imag();
      wcs=h2.cs[w];
      wsn=h2.sn[w];
      for(m=0; m<nframe; m++) {
        u[m] += cs[0] * wcs[m] + sn[0] * wsn[m];
        v[m] += cs[1] * wcs[m] + sn[1] * wsn[m];
        }
      }

/**----------------------------------------------------------------------------
    depth-integrated bottom friction */
    for(m=0;m<nframe;m++) {
//      c=sqrt(u[m]*u[m]+v[m]*v[m]) + data.u0[n];
      c=sqrt(u[m]*u[m]+v[m]*v[m]);
      taux[m]=c*u[m];
      tauy[m]=c*v[m];
      }
/**----------------------------------------------------------------------------
    analyze time serie */
//  status=harmonic_compute(taux,residuals, nframe, h1, analysis, maxstep, ax, Gx);
    status=harmonic_compute(taux, NULL, nframe, h1, analysis, maxstep, ax, Gx);
//  status=harmonic_compute(tauy, residuals, nframe, h1, analysis, maxstep, ay, Gy);
    status=harmonic_compute(tauy, NULL, nframe, h1, analysis, maxstep, ay, Gy);
/**----------------------------------------------------------------------------
    store constants */
    for(w = 0; w < analysis.n; w++) {
      bufx[w][n]=polar(ax[w],-Gx[w]);
      bufy[w][n]=polar(ay[w],-Gy[w]);
      }
//     delete[] u;
//     delete[] v;
//     delete[] H;
//     delete[] taux;
//     delete[] tauy;
//     delete[] ax;
//     delete[] Gx;
//     delete[] ay;
//     delete[] Gy;
    }

  zaparr(time);
  zaparr(residuals);

  delete[] u;
  delete[] v;
  delete[] taux;
  delete[] tauy;
  delete[] ax;
  delete[] Gx;
  delete[] ay;
  delete[] Gy;
  
  h1.destroy();
  h2.destroy();
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// void cefmo_frictionRHS(spectrum_t prediction, spectrum_t analysis, double duration,
//                        atlas2D_t atlas, parameter_t data, int nndes, complex<double> **bufx, complex<double> **bufy)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   double sampling=1800.;
//   
// //  if(duration==0.0)
//     harmonic_optimize(analysis, duration, sampling, &duration);
//   
//   switch(gOPENMP_nCPUs) {
//     case 1:
//       cefmo_frictionRHS_SEQ(prediction, analysis, duration, atlas, data, nndes, bufx, bufy);
//       break;
//     default:
//       cefmo_frictionRHS_OMP(prediction, analysis, duration, atlas, data, nndes, bufx, bufy);
//       break;
//     }
// }

