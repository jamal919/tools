


#ifndef STURM_LIOUVILLE_H

#define STURM_LIOUVILLE_H

  extern int compute_Wmode_v1(int nlevels, double *Z, double *density,double *modes);
  extern int compute_Wmode_v2(int nlevels, double *z, double *rho, double g, double *N, double *c, double **modes, int verbose);
  extern int compute_Wmode_v3(int nlevels, double *z, double *rho, double g, double *N, double *c, double **modes, bool rigid_lid, int verbose);
  extern int compute_Wmode_v4(int nlevels, double *z, double *rho, double g, double *N, double *c, double **modes, bool rigid_lid, int verbose);
  extern int compute_Wmode_v5(int nlevels, double *z, double *rho, double g, double *N, double *c, double **modes, int & nmodes, bool rigid_lid, int verbose);
  extern int compute_Wmode_v6(int nlevels, double *z, double *rho, double g, double *N, double *c, double **modes, int & nmodes, bool rigid_lid, int verbose);

  extern int compute_Pmode(int nlevels, double *z, double *rho, double g, double *N, double *c, double **modes, int & nmodes, bool rigid_lid, int verbose);
  
  extern int normalize_Wmodes(int nlevels, double *z, double **modes, double mask, int nmodes, double *R, double *N, int verbose);
  extern int normalize_Pmodes(int nlevels, double *z, double **modes, double mask, int nmodes, double *R, int verbose);
  extern int normalize_Umodes(int nlevels, double *z, double **modes, double mask, int nmodes, double *R, int verbose);

  extern int Umodes_decomposition_v1(int nlevels, double *z, double **modes, int nmodes, double *rho, complex<double> *u, complex<double> *decomposition, int verbose);
  
  extern int Umodes_decomposition_v2(int nlevels, double **modes, int nmodes, double          *u, double          *decomposition, bool debug, int verbose);
  extern int Umodes_decomposition_v2(int nlevels, double **modes, int nmodes, complex<float>  *u, complex<float>  *decomposition, bool debug, int verbose);
  extern int Umodes_decomposition_v2(int nlevels, double **modes, int nmodes, complex<double> *u, complex<double> *decomposition, bool debug, int verbose);
  
  extern int Wmodes_decomposition(int nlevels, double *z, double **modes, int nmodes, double *rho, double *N, complex<double> *u, complex<double> *decomposition, int verbose);
  extern int Wmodes_decomposition(int nlevels, double *z, double **modes, int nmodes, double *rho, complex<double> *u, complex<double> *decomposition, int verbose);

  extern int Wmodes_ProjectionMatrix_v1(int nlevels, double *z, double **modes, int nmodes, double *rho, double *N, double *coefficients, int verbose);
  extern int Wmodes_ProjectionMatrix_v2(int nlevels, double **modes, int nmodes, double *coefficients, bool debug, int verbose);
  extern int Wmodes_unresolved(int nlevels, double *z, double **modes, int nmodes, double *rho, double *N, double *coefficients, double *inverse, bool debug, int verbose);
  
  extern double ComputeModesNumerics(double omega, double c, double s, double dx);
 
#endif