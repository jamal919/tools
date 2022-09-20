#ifndef ICE_H
#define ICE_H

#ifdef MAIN_SOURCE
#define PREFIX 
#else
#define PREFIX extern
#endif

/**         for ICE         **/

typedef struct
  {
  double u,v; /* ice velocity*/
  double a,h,hsnow,t[3]; /* concentration, volume, temperature */
  double p; /* internal ice pressure */
  double tau; /* frictional coefficient with water */
  } ice_t;


PREFIX int use_ice;
PREFIX double hinit, ainit, hsnoi;
PREFIX double dtice,tice;
PREFIX double rhoice;
PREFIX double thetaice;
PREFIX double uimin;
PREFIX complex<double> CCwi;
PREFIX double pstar, Cstrength, zetamax, etamax, Delta0, alpha_mr, eccent_VP,
              ellip2_VP, phi_cavitating, sinphi_cavitating;
PREFIX double ee_EVP;
PREFIX int istep_EVP;

PREFIX ice_t *Feice;
PREFIX char icefile[255];
PREFIX int ice_mom, ice_rheo;
PREFIX double **stressice;

#endif /* ICE_H */
