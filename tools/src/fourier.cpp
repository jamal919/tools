#include <config.h>

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "tools-structures.h"
#include "functions.h"
#include "rutin.h"
#include "fourier.h"
#include "tides.h"

#include "tools-define.h"

#if HAVE_LIBGSL == 1

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_math.h>
 
#elif HAVE_LIBFFTPACK == 1

//dfft sont des routines de fftpack ....
extern "C" {
extern void dffti_(int *, double *);
extern void dfftf_(int *, double *,double *);
}

#endif

#include <gsl/gsl_statistics_double.h>


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void psd_r1_c(double dx,const double *z, double mask, int n, double *puissance,
                double *phase,double *f,int *nf,int option,int *status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 Calcul de la densite spectrale de puissance

------------------------------------------------------------------------
 Inputs :
  N    = nombre d'observations
  Z    = tableau de dimension n contenant la sequence a transformer
  dX   = echantillonnage de la serie
  option=0 : Spectre de Puissance (unite de Z)^2\n
  option=1 : Densite spectrale de Puissance (unite de Z)^2/[cycle/(unite de x)]\n
  option=2 : Densite spectrale de Puissance en VARIANCE CONSERVEE (unite de Z)^2\n
------------------------------------------------------------------------
 Outputs :
  TF   = tableau de dimension nobs contenant les coeff. de la transformee
  Puissance= tableau de dimension (nobs+1)/2 des coefficients de la dsp
  Phase    = tableau de dimension (nobs+1)/2 des coefficients de la phase
------------------------------------------------------------------------

   On stocke la frequence 0 !!!

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int i;
  double *amplitude=NULL,*tf=NULL,tobs,a,b;

  if(n%2==0) {
    *nf=n/2;
    }
  else {
    *nf=(n+1)/2;
    }
  exitIfNull(
    amplitude=new double[*nf]
    );
  
  tf=copy(z,n);

/*
  if(option==0)
    printf("Spectre de Puissance (unite de Z)^2\n");
  else if(option==1)
    printf("Densite spectrale de Puissance (unite de Z)^2/[cycle/(unite de x)]\n");
  else if(option==2)
    printf("Densite spectrale de Puissance en VARIANCE CONSERVEE (unite de Z)^2\n");
*/

  FILE *ctrlF=0;
  const int minnTest=8;
  if(false and n>minnTest){
    /* initialise FT test */
    n=minnTest;
    
    const int nf=n/2+1;
    double rs[nf],is[nf];
    int j;
    double t,tj;
    
    for(j=0;j<nf;j++){
      rs[j]=2*j;
      is[j]=2*j+1;
      }
    rs[0]=1e-3;
    
    vector<double> Cs,Ss;
    
    for(i=0;i<n;i++){
      double *zi=&tf[i];
      t=-i*M_PIx2/n;
      
      *zi=0.;
      for(j=0;j<nf;j++){
        tj=t*j;
        *zi+=rs[j]*cos(tj)+is[j]*sin(tj);
        Cs.push_back(cos(tj));
        Ss.push_back(sin(tj));
        }
      }
    
    ctrlF=fopen("psd_r1_c.py","w");
    ASSERT_ARRAY(ctrlF,"rs",rs,nf,1);
    ASSERT_ARRAY(ctrlF,"Is",is,nf,1);
    ASSERT_ARRAY(ctrlF,"Cs",Cs,Cs.size(),1);
    fprintf(ctrlF,"Cs=matrix(Cs.reshape([%d,%d]))\n",n,nf);
    ASSERT_ARRAY(ctrlF,"Ss",Ss,Ss.size(),1);
    fprintf(ctrlF,"Ss=matrix(Ss.reshape([%d,%d]))\n",n,nf);
    ASSERT_ARRAY(ctrlF,"z",tf,n,1);
    }

/*------------------------------------------------------------------------------
  Calcul de la transformee de Fourier de la serie TF */
#if HAVE_LIBGSL == 1
  {
  //see info:/gsl-ref/Mixed-radix FFT routines for real data
  gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc (n);
  gsl_fft_real_wavetable *wata = gsl_fft_real_wavetable_alloc (n);
  
  *status=gsl_fft_real_transform (tf, 1, n, wata, work);
  
  gsl_fft_real_wavetable_free (wata);
  gsl_fft_real_workspace_free (work);
  }
#elif HAVE_LIBFFTPACK == 1
  wsave=new double[2*n+15];
  if(wsave == NULL) {
    STDOUT_BASE_LINE("error in allocating memory for wsave \n");
    exit(99);
    }
  dffti_(&n,         wsave);
  dfftf_(&n, &tf[0], wsave);
  delete[] wsave;
#else
  STDOUT_BASE_LINE("Neither GSL nor FFTPACK available, abort... \n");
  exit(-1);
#endif
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  nf   : nombre de frequences dans le spectre (0 incluse)
  Tobs : longueur d'observation

  Attention, Tobs=N*dX !!! (et non pas (N-1)*dX)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  tobs=n*dx;

/*------------------------------------------------------------------------------
  Calcul de la position des estimees spectrales */
  for(i=0;i<*nf;i++)
    f[i]=(double)(i)/tobs;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  NORMALISATION effectuee par les routines utilisees:

  Les subroutines multiplient les coefficients de Fourier par N.
  De plus, un facteur 1/2 intervient pour tenir compte du fait que
  la serie est decomposee en cos,sin et non exponentielle complexe conjuguee
  seulement pour ces valeurs qui sont decomposees en cos,sin donc pas pour
  la fréquences 0, ni celle de Nyquist dans le cas de N multiple de 2.

  2/T x integrale= 2/N x serie

  Avec cette normalisation, on retourne les amplitudes physiques

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  tf[0]/=n;
  int imax;
  if(n&1)
    imax=n;
  else{
    imax=n-1;
    tf[imax]/=n;
    }
  for(i=1;i<imax;i++)
    tf[i]*=2./n;

  if(ctrlF!=0){
    ASSERT_ARRAY(ctrlF,"tf",tf,n,1);
    fclose(ctrlF);
    TRAP_ERR_EXIT(ENOEXEC,"testing %s\n",__func__);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  Calcul de la puissance et de la phase associees a chaque harmonique:
  
  Amplitude          = sqrt(a^2+b^2)
  Amplitude efficace = Amplitude/sqrt(2)
  Puissance          = Amplitude^2
  
  ==> Puissance efficace = (Amplitude efficace)² = Amplitude²/2
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  1ere frequence <=> frequence 0 */
  amplitude[0] = tf[0];
  puissance[0] = amplitude[0]*amplitude[0];
  phase[0]     = 0.;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  Spectre CLASSIQUE: puissance en (unite de Z)^2/[cycle/(unite de x)]

  Pour la frequence zero, on a une puissance infinie car f(0)=0, donc
  pbm de visu ==> on impose Puissance(0)=0. pour pouvoir faire la visu.

  Du fait de la normalisation de DTFR_R1, on a:

  OneSided(i)=G(i)=2*S(i)=[a(i)**2+b(i)**2]*T/2

  Spectre en VARIANCE CONSERVEE = puissance en (unite de x)^2

  Si l'abscisse est en log, on fait le changement de variable: u=ln(f)
  d'ou du=1/f df

  OneSided(i)=G(i)*f(i)*ln(10)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(i=1;i<*nf;i++) {
    a=tf[2*(i-1)+1];
    if(i<=*nf-1 or n&1)
      b=tf[2*(i-1)+2];
    else
      b=0;
    fct_xy2polar(a, b, amplitude[i], phase[i]);
    if(option == 0)
      puissance[i]=(a*a+b*b);
    else if(option == 1)
      puissance[i]=(a*a+b*b)*tobs/2.0;
    else if (option == 2)
      puissance[i]=f[i]*log(10.)*(a*a+b*b)*tobs/2.0;
    }

  if(0){
    ctrlF=fopen("psd_r1_c.py","a");
    ASSERT_ARRAY(ctrlF,"z",z,n,1);
    ASSERT_ARRAY(ctrlF,"option",&option,1,1);
    ASSERT_ARRAY(ctrlF,"power",puissance,*nf,1);
    ASSERT_ARRAY(ctrlF,"phase",phase,*nf,1);
    ASSERT_ARRAY(ctrlF,"f",f,*nf,1);
    ASSERT_ARRAY(ctrlF,"status",status,1,1);
    fclose(ctrlF);
    TRAP_ERR_EXIT(ENOEXEC,"testing\n");
    }
  
  delete[]amplitude;
  delete[]tf;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fourier_r1_c(double dx, double *z, double mask, int n, double *tf,
                double *f,int *nf,int *status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

 Calcul de la transformee de fourier de z,
 retourne la tableau des coefficients (fcomplexes) de cette transformee.

------------------------------------------------------------------------
 Inputs :
  N    = nombre d'observations
  Z    = tableau de dimension n contenant la sequence a transformer
  dX   = echantillonnage de la serie
------------------------------------------------------------------------
 Outputs :
  TF   = tableau de dimension nobs contenant les coeff. de la transformee
------------------------------------------------------------------------

   On stocke la frequence 0 !!!

-----------------------------------------------------------------------*/
{
  int i,m;
  double xm,tobs;

  *nf=1+n/2;
  
/*  tf =(double *) malloc(n*sizeof(double)); */

  
 
/*------------------------------------------------------------------------------
     Calcul de la valeur moyenne de la serie z
------------------------------------------------------------------------*/

  xm=0.0;
  m=0;
  for(i=0;i<n;i++) {
  if(z[i] != mask) {
    xm=xm+z[i];
    m=m+1;
    }
  }
  
  xm=xm/m;
/*  printf("mean=%lf\n",xm); */

/*------------------------------------------------------------------------------
    On retranche cette valeur moyenne a la serie z
-----------------------------------------------------------------------*/

  for(i=0;i<n;i++) {
  if(z[i] != mask)
    tf[i]=z[i]-xm;
  else
    tf[i]=0.0;
  }

/*------------------------------------------------------------------------------
    Calcul de la transformee de Fourier de la serie TF
----------------------------------------------------------------------*/

#if HAVE_LIBGSL == 1
  {
  //see info:/gsl-ref/Mixed-radix FFT routines for real data
  gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc (n);
  gsl_fft_real_wavetable *wata = gsl_fft_real_wavetable_alloc (n);
  
  gsl_fft_real_transform (tf, 1, n, wata, work);
  
  gsl_fft_real_wavetable_free (wata);
  gsl_fft_real_workspace_free (work);
  }
#elif HAVE_LIBFFTPACK == 1
//with the fortran lib fftpack
  wsave=(double *)malloc((2*n+15)*sizeof(double));
  if(wsave == NULL) {
    STDOUT_BASE_LINE("error in allocating memory for wsave \n");
    exit(99);
    }
 dffti_(&n,         wsave);
 dfftf_(&n, &tf[0], wsave);

 free(wsave);
#else
  STDOUT_BASE_LINE("Neither GSL nor FFTPACK available, abort... \n");
  exit(-1);
#endif

/*------------------------------------------------------------------------------

     Nf   : nombre de frequences dans le spectre (0 incluse)
     Tobs : LONGueur d'observation

        Attention, Tobs=N*dX !!! (et non pas (N-1)*dX)

----------------------------------------------------------------------*/

  tobs=n*dx;

/*------------------------------------------------------------------------------
 Calcul de la position des estimees spectrales
----------------------------------------------------------------------*/

  for(i=0;i<*nf;i++)
    f[i]=(double)(i)/tobs;

/* *-----------------------------------------------------------------------

 NORMALISATION effectuee par les routines utilisees:

 Les subroutines multiplient les coefficients de Fourier par N.
 De plus, un facteur 1/2 intervient pour tenir compte du fait que
 la serie est decomposee en cos,sin et non exponentielle fcomplexe conjuguee

 2/T x integrale= 2/N x serie

 Avec cette normalisation, on retourne les amplitudes physiques

----------------------------------------------------------------------*/
  for(i=0;i<n;i++)
    tf[i]=2.0*tf[i]/n;

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <class T> size_t value_padding(vector<T>& data, vector<T>& times, const T Trepet, const T value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const int nmes=data.size();
  if( Trepet<=0. or nmes<=0 or nmes!=times.size() ){
    return -1;
    }
  const T t0=times[0];
  // define new length
  const int padded_n = (times.back() - t0) / Trepet + 1;
  
/*------------------------------------------------------------------------------
  build extended time vector
  unsafe if original time serie is not time-sorted */
  vector<T> padded_t(padded_n);
  vector<int> data_index(nmes);
  
  size_t n,m;
  for(n = 0; n < padded_n; n++){
    padded_t[n]=times[0] + n * Trepet;
    }
  for(m = 0; m < nmes; m++) {
    const T *t=&times[m];
    n=(*t-times[0])/Trepet;
    padded_t[n] = *t;
    data_index[m]=n;
    }
  
/*------------------------------------------------------------------------------
  build value-padded dataset */
  vector<double> padded_d(padded_n, value);
  for(m = 0; m < nmes; m++){
    padded_d[data_index[m]] = data[m];
    }
  
/*------------------------------------------------------------------------------
  return values */
  data.swap(padded_d);
  times.swap(padded_t);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <class T> size_t apodisation(const string & apodize_function, size_t len, vector<T>& w)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  size_t a = len / 2;
  double sW = 0;
  
/*------------------------------------------------------------------------------
  original apodization function from kstpsd.cpp*/
  if(apodize_function=="original"){
    for (size_t i = 0; i < len; i++){
      w.push_back(1. - cos(2 * M_PI * i / len));
      sW += w.at(i) * w.at(i);
      }
    sW = sqrt(sW / len);
    for(typename std::vector<T>::iterator p = w.begin(); p != w.end(); p++){
      *p /= sW;
      }
    return(0);
    }

/*------------------------------------------------------------------------------
  Bartlett function from kstpsd.cpp*/
  if(apodize_function=="bartlett"){
    for(size_t i = 0; i < len; i++) {
      double x = i - a;
      w.push_back(1 - fabs(x) / a);
      }
    return(0);
    }

/*------------------------------------------------------------------------------
  Blackman function from kstpsd.cpp*/
  if(apodize_function=="blackman"){
    for(size_t i = 0; i < len; i++){
      double x = i - a;
      w.push_back(0.42 + 0.5 * cos(M_PI * x / a) + 0.08 * cos(2 * M_PI * x / a));
      }
    return(0);
    }

/*------------------------------------------------------------------------------
  Connes function from kstpsd.cpp*/
  if(apodize_function=="connes"){
    for (size_t i = 0; i < len; i++){
      double x = i - a;
      w.push_back(gsl_pow_2(1 - (x * x) / (a * a)));
      }
    return(0);
    }

/*------------------------------------------------------------------------------
  cosine function from kstpsd.cpp*/
  if(apodize_function=="cosine"){
    for (size_t i = 0; i < len; i++) {
      double x = i - a;
      w.push_back(cos((M_PI * x) / (2 * a)));
      }
    return(0);
    }

  /*------------------------------------------------------------------------------
  Gaussian function */
  if(apodize_function=="gaussian"){
    double gaussianSigma = 0.4;
    for (size_t i = 0; i < len; i++) {
      double x = i - a;
      w.push_back(exp((-1 * x * x) / (2 * gaussianSigma * gaussianSigma)));
      }
    return(0);
    }

/*------------------------------------------------------------------------------
  Hamming function */
  if(apodize_function=="hamming"){
    for (size_t i = 0; i < len; i++) {
      double x = i - a;
      w.push_back(0.54 + 0.46 * cos(M_PI * x / a));
      }
    return(0);
    }


/*------------------------------------------------------------------------------
  Hann function */
  if(apodize_function=="hann"){
    for (size_t i = 0; i < len; i++){
      double x = i - a;
      w.push_back(gsl_pow_2(cos((M_PI * x) / (2 * a))));
      }
    return(0);
    }


/*------------------------------------------------------------------------------
  Welch function */
  if(apodize_function=="welch"){
    for (size_t i = 0; i < len; i++) {
      double x = i - a;
      w.push_back(1 - (x * x) / (a * a));
      }
    return(0);
    }
 
/*------------------------------------------------------------------------------
  rectangular function  */
  if(apodize_function=="rectangular"){
    for (size_t i = 0; i < len; i++){
      w.push_back(1.);
      }
    return(0);
    }
  else {
/*------------------------------------------------------------------------------
    default: original function  */
    for (size_t i = 0; i < len; i++){
      w.push_back(1. - cos(2 * M_PI * i / len));
      sW += w.at(i) * w.at(i);
      }
    sW = sqrt(sW / len);
    for(typename std::vector<T>::iterator p = w.begin(); p != w.end(); p++){
      *p /= sW;
      }
    return(1);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<double> spectrum_psd_average(const spectrum_t& s, const double *residuals, const double *times,
                                      const int nmes, double Trepet, bool apodization)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  compute error estimates based on non tidal ocean signal's
  contamination on tidal peaks
  
  inspired by Ponchaut et al. (2001 http://dx.doi.org/10.1175/1520-0426(2001)018<0077:AAOTTS>2.0.CO;2)

  NOTES
  times and Trepet in days
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
/*------------------------------------------------------------------------------
  local copy */
  const double
    mean = gsl_stats_mean(residuals, 1, nmes);
  vector<double> local_residuals(nmes);
  vector<double> local_times(nmes);
  int status = -1;
  int i;

/*------------------------------------------------------------------------------
  remove mean residuals */
  for(i = 0; i < nmes; i++){
    local_residuals[i]=residuals[i] - mean;
    local_times[i]=times[i];
    }

/*------------------------------------------------------------------------------
  zero-padding */
  status=value_padding(local_residuals, local_times, Trepet, 0.);
  if( status != 0){
    gmerror((string)"value_padding() failed");
    }
  
/*------------------------------------------------------------------------------
  re-scaling */

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
Zero-padding adds 0s to the series, so the series must be rescaled by:
  - sqrt(N1/N0) to compensate for the 0s added
  - sqrt(N1/N0) to compensate for the extra length that will diminish the amplitudes of the noise
and the product of these factors gives N1/N0.
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  size_t n = local_times.size();
  const double
    rescale=(double)n/nmes;
  for(i=0;i<n;i++)
    local_residuals[i]*=rescale;
  
/*------------------------------------------------------------------------------
  apodization */
  if(apodization == true){
    vector<double> w;
    
    apodisation("original", n, w);
    
    if(! w.empty()){
      for(size_t m = 0; m < w.size(); m++){
        local_residuals.at(m) *= w[m];
        }
      }
    
    }

/*------------------------------------------------------------------------------
  real FFT of the input signal  */

  double *res = new double[n];
  for(i = 0; i < n; i++){
    res[i] = local_residuals[i];
    }
  
  double *p = new double[n]; // unit will be unit of res
  double *phase = new double[n];
  double *f = new double[n]; // unit will be 1/(unit of Trepet)
  int nf = 0;
  status = -1;
  
/*------------------------------------------------------------------------------
  compute power spectrum */
  int option=0;
  psd_r1_c (Trepet, res, 99.9999, n, p, phase, f, &nf, option, &status);
  delete [] res;
  
/*------------------------------------------------------------------------------
  RMS for [ freq-K*df ; freq+K*df ] */
  const double
    trecord = times[nmes - 1] - times[0];  // time serie duration
  /* K is at least one per year and
  at least 7 so that there are 7*2+1=15 samples
  so 15*2=30 real and imaginary parts taken into account */
  const int
    K = MIN( ceil(trecord/365.25*2) , 7);
  
  double *running_rms = new double[nf];
  double sum_square,mean_square;
  int i0,ii,i1,nfi;
  for(i = 0; i < nf ; i++){
    i0=min(   0, i-K);
    i1=max(nf-1, i+K);
    nfi=i1-i0+1;
    sum_square=0.;
    for(ii=i0;ii<=i1;ii++)
      sum_square+=p[ii];
    mean_square=sum_square/(2*nfi-1);/* RMS on the population */
    running_rms[i] = sqrt(mean_square);
    }
  
/*------------------------------------------------------------------------------
  Noise level estimation */
  vector<double> error(s.n, -1.);
  double omega_cpd;
  const double alias=1./Trepet;
  for(size_t k = 0; k < s.n; k++){
    omega_cpd=s.waves[k].omega*dph2cpd;
    alias_frequency(alias,&omega_cpd);
    i=round(omega_cpd*trecord);/* same as dividing by df=1-trecord */
    error[k]=running_rms[i];
    }
  
  delete [] f;
  delete [] running_rms;
  
  return(error);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

vector<double> spectrum_psd_average_obsolete(const spectrum_t& s, const double *residuals, const double *times,
                                    const int nmes, double Trepet, bool apodization)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* -----------------------------------------------------------------------------
  compute error estimates based on non tidal ocean signal's
  contamination on tidal peaks
  
  translated from aktarus sources,
  inspired by Ponchaut et al., 1999

  validated for hourly predictions
  validated for altimetry predictions

  NOTES
  times and Trepet in days
  
-----------------------------------------------------------------------------*/
/* ----------------------------------------------------
  local copy */
  double mean = gsl_stats_mean(residuals, 1, nmes);
  vector<double> local_residuals;
  vector<double> local_times;
  int status = -1;

  
  for(size_t i = 0; i < nmes; i++){
    local_residuals.push_back(residuals[i] - mean); // remove mean residuals
    local_times.push_back(times[i]);
    }

/*------------------------------------------------------------------------------
  zero-padding */
  status=value_padding(local_residuals, local_times, Trepet, 0.);
  if( status != 0){
    gmerror((string)"value_padding() failed");
    }
  size_t n = local_times.size();

/*------------------------------------------------------------------------------
  apodization */
  if(apodization == true){
    vector<double> w;
    apodisation("original", n, w);
    if(! w.empty()){
      for(size_t m = 0; m < w.size(); m++){
        local_residuals.at(m) *= w[m];
        }
      }
    w.clear();
    }

/*------------------------------------------------------------------------------
  real FFT of the input signal  */
  if( fmod(n, 2.0) != 0){
    n--;
    }

  double *res = new double[n];
  for(size_t i = 0; i < n; i++){
    res[i] = local_residuals[i];
    }
 
#if defined(DEBUG)
  FILE *debugFile = NULL;
  debugFile=fopen("debug.time_series.dat", "w");
  if(debugFile==0){
      TRAP_ERR_EXIT(errno,"file opening issue : debug.time_series.dat (%d %s)\n",errno,strerror(errno));
    }
  for(size_t i = 0; i < n; i++){
    fprintf(debugFile, "%12.6f %9.6f %9.6f %9.6f\n", local_times[i], residuals[i], local_residuals[i], res[i]);
    }
  fclose(debugFile);
#endif


  double *p = new double[n]; // unit will be unit of res
  double *phase = new double[n];
  double *f = new double[n]; // unit will be 1/(unit of Trepet)
  int nf = 0;
  status = -1;
  
/*------------------------------------------------------------------------------
  compute power spectrum */
  int option=0;
  psd_r1_c (Trepet, res, 99.9999, n, p, phase, f, &nf, option, &status);// LR, comment: return status is non applicable
/*------------------------------------------------------------------------------
  convert power to amplitude */
  double *a = new double[n];
  for (size_t i = 0; i < nf ; i++){
    a[i] = sqrt(p[i]);
    }
  delete [] p;
  delete [] phase;
  
  /* ----------------------------------------------------
     internal checks */

#if defined(DEBUG)
  double max = 0 ;
  size_t index_max = 0;
  for (size_t i = 0; i < nf ; i++){
    if (max < a[i]){
      max = a[i];
      index_max = i;
       }
    }
  double Tmax = 1. / f[index_max];
  cout << "#DBG max= " << max
        << " at index= " << index_max
        << " (T= " << Tmax << " d)" << endl;

  double min = 999999.;
  size_t index_min = 0;
  for (size_t i = 1; i < nf ; i++){
    if (min > a[i]){
      min = a[i];
      index_min = i;
       }
    }

  double Tmin = 1. / f[index_min];
  cout << "#DBG min= " << min
        << " at index= " << index_min
        << " (T= " << Tmin << " d)" << endl;

  FILE *fftFile = NULL;
  fftFile=fopen("debug.fft", "w");
  if(fftFile==0) {
      TRAP_ERR_EXIT(errno,"file opening issue : debug.fft (%d %s)\n",errno,strerror(errno));
    }
  for(size_t i = 1; i < nf; i++){
    fprintf(fftFile, "%12.6f %9.4f\n", 1. / f[i], a[i]);
    }
  fclose(fftFile);
#endif

  local_residuals.clear();
  local_times.clear();
  delete [] res;


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  series of length T, sampling dT, N samples
  
  T=(N-1)xdT
  
  Highest frequency: f=1/(2xdT) Niquyst
  
  Lowest frequency:  f=1/T
  
  Fourier resolution: 1/T
  
  
  Averaging energy spectrum for f-K*df/2 < freq < f+K*df/2
  
  Presumely equivalent to: average of analysis over T/K
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
/*------------------------------------------------------------------------------
  Averaging for freq-df/2 < freq <freq+df/2 */
  double trecord = (times[nmes - 1] - times[0]) /  365.25;  // time serie duration in year
  int df = 1 + floor(trecord) / 0.5;                        // df frequencies ( equiv. 0.5 years) average

  df=min(df,5);
  
  double *mean_fi = new double[nf - df];
  for(size_t i = 0; i < (nf - df); i++){
    mean_fi[i] = gsl_stats_mean(a + i, 1, df);
    }

#if defined(DEBUG)
  fftFile = fopen("debug.meanfft", "w");
  if(fftFile==0) {
    TRAP_ERR_EXIT(errno,"file opening issue : debug.meanfft (%d %s)\n",errno,strerror(errno));
    }
  for(size_t i = 0; i != nf - df; ++i){
    fprintf(fftFile, "%12.6f %9.6f\n", 1. / f[i + df / 2], mean_fi[i]);
    }
  fclose(fftFile);
#endif

/*------------------------------------------------------------------------------
  Noise level estimation */
  double dph2cpd = 360. / 24.;
  std::vector<double> error(s.n, -1.);
  double fe = 1. / Trepet;
  
  for(size_t k = 0; k < s.n; k++){
    bool found = false;
#if defined(DEBUG)
    printf("%10s, aliased at %9.4f days, ", s.waves[k].name, dph2d / s.waves[k].aliased);
#endif
/*------------------------------------------------------------------------------
    interpolate averaged residual amplitude */
    for(size_t i = 0; i < (nf - df) - 1; i++){
      double r= s.waves[k].aliased / dph2cpd - f[i + df / 2];
      if( fabs(r) < (fe / n) ){                                        // freq in 1/day
        error[k] = mean_fi[i] +r*(mean_fi[i + 1] - mean_fi[i]) / (fe / n);
#if defined(DEBUG)
        printf("[%9.6f %9.6f], ", 1. / f[i + df / 2 + 1],  1. / f[i + df / 2]);
        printf("interpolate averaged amplitudes %9.6f %9.6f => %9.6f\n",
               mean_fi[i], mean_fi[i + 1], error[k]);
#endif
        found = true;
        break;
      }
    }

/*------------------------------------------------------------------------------
    backup solution: interpolate raw residual amplitude */
    if(found != true){
      for(size_t i = 0; i < nf - 1; ++i){
        double r= s.waves[k].aliased / dph2cpd - f[i];
        if( fabs(r) < (fe / n) ){                                       // freq in 1/day
          error[k] = a[i]+r*(a[i + 1] - a[i]) / (fe / n);
#if defined(DEBUG)
          printf("[%9.6f %9.6f], ",  1. / f[i + 1],  1. / f[i]);
          printf("interpolate raw amplitudes %9.6f %9.6f => %9.6f\n",
                 a[i], a[i + 1], error[k]);
#endif
          found = true;
          break;
        }
      }
    }
  
    // backup solution: assign raw residual amplitude
    if(found != true){
       for(size_t i = 0; i < nf - 1; ++i){
          double r= s.waves[k].aliased / dph2cpd - f[i];
          if( fabs(r) > (fe / n) ){ // freq in 1/day
            error[k] = a[i]; /// LR, BUG error could be mean value if i = 0, so 0 as average is removed, when s.waves[k].aliased / dph2d is close to Nyquist period
            found = true;
            break;
          }
       }
    }

  }
  delete [] mean_fi;

#if defined(DEBUG)
  fftFile = fopen("debug.avgfft", "w");
  if (fftFile == NULL) {
    TRAP_ERR_EXIT(errno,"file opening issue : debug.avgfft (%d %s)\n",errno,strerror(errno));
    }
  for(size_t i = 0; i < s.n; i++){
    fprintf(fftFile, "%10s %9.4f %9.6f\n", s.waves[i].name, dph2d / s.waves[i].aliased, error[i]);
    }
  fclose(fftFile);
#endif

  delete [] f;
  delete [] a;

  return(error);

}

