//functions de fourier.c
extern void fct_dxytopolar_c(double u,double v,double *ro,double *teta);
extern void psd_r1_c(double dx,const double *z, double mask, int n, double *puissance,
                double *phase,double *f,int *nf,int option,int *status);
extern void fourier_r1_c(double, double *, double, int, double *, double *,int *,int *);
extern vector<double> spectrum_psd_average(const spectrum_t& s, const double *residuals, const double *times,
                                           const int nmes, const double Trepet, bool apodization);
