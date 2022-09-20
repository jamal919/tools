// #include "tools-structures.h"

#ifndef _STATISTIC_H
#define _STATISTIC_H

#include "poc-array-stats.hpp"

class statistic_t{
private :
public :
  double mean,median,mae;
  double std,rms;
  double cf,pof,nof,pof2x,nof2x;
  double min,max;
  
  statistic_t() {
    mean=median=mae=0;
    std=rms=0;
    cf=pof=nof=pof2x=nof2x;
    min=max=0;
    }
  
  double variance() {
    double v;
    v=std*std;
    return(v);
    }
  
  char *print() {
    char *s=new char[256];
    sprintf(s,"mean=%f std=%f",mean,std);
    return(s);
    }
  
};

class statistic_c_t{
private :
public :
  complex<double> mean, min,max;
  double std;
  
  statistic_c_t() {
    mean=0;
    std=0;
    min=max=0;
    }
  
  double variance() {
    double v;
    v=std*std;
    return(v);
    }
  
};

typedef struct {
  double trend,bias;
} regression_t;

class covariance_t {
private :
public :
  statistic_t s1,s2;
  double correlation, covariance;
  int n;
  
  covariance_t() {
    correlation=0;
    n=0;
    }
  
};

extern int print_stats(int n,int count,const complex<double> & mean,double variance,char unit='d');

extern statistic_t OutliersFrequency(double *h, double mask, int n, double limit);

extern statistic_t get_statistics(const int *h, int mask, size_t n, int verbose);

extern statistic_t get_statistics(const double *h, double mask, int n);
extern statistic_t get_statistics(const double *h, double mask, int n, int verbose);
extern statistic_t get_statistics(const double *h, double mask, range_t<size_t> range, int verbose);

extern statistic_t get_statistics(const vector<float> &h, float mask=1.e+10, int verbose=1);
extern statistic_t get_statistics(const vector<double> &h, double mask=1.e+10, int verbose=1);

extern statistic_t get_statistics(const float *h, float mask, int n);
extern statistic_t get_statistics(const float *h, float mask, int *list, int n, int verbose);
extern statistic_t get_statistics(const float *h, float mask, int n, int verbose);

extern statistic_t get_statistics(const double *h, int n);
extern statistic_t get_statistics(const float *h, int n);

extern statistic_c_t cget_statistics(complex<float> *h, complex<float> mask, int n);

// extern statistic_t sget_statistics01(float *h, float mask, int n, int verbose);

// template <typename type> statistic_t get_statistics(type *h, type mask, int n);
// template <typename type> statistic_t get_statistics(type *h, int n);

extern statistic_t sget_filteredstatistics(float *h, float mask, int n, float threshold);
extern statistic_t sget_filteredstatistics2(float *h1,float *h2, float mask, int n, float threshold);

extern regression_t get_trend(float *h, float *t, float mask, int n);

extern covariance_t get_covariance(float  *h1, float  *h2, float  mask, int n, int verbose=0);
extern covariance_t get_covariance(double *h1, double *h2, double mask, int n, int verbose=0);

extern double *get_autocorrelation(double *h,int n,int range,int test);
extern double *get_maxcorrelation(double *h1,double *h2,int n,int range, int test, int & max_positive, int & max_negative, int verbose=1);

extern void poc_minmax(double *v,int n,double mask,double *min,double *max);
extern void poc_minmax(float *v,int n,float mask,float *min,float *max);

extern range_t<double> poc_minmax(double *x, int nvalues, double mask);
extern range_t<float>  poc_minmax(float *x, int nvalues, float mask);
extern range_t<int>    poc_minmax(int *x, int nvalues, int mask);

extern range_t<int>    poc_minmax(const vector<int> & x, int mask);
extern range_t<double> poc_minmax(const vector<double> & x, double mask);

extern double get_timesampling(const double *t, double mask, int n, int verbose=0);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int median_shell (T *a, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int i,j,k;
  T aa;
  
  k=0;

  while(k<n) {
    k=3*k+1;
    }
    
  while(k>1) {
    k= k/3;
    if (k<0) return(0);
  
    for(i= k+1;i<n;i++) {
      aa= a[i];
      j= i;
_30:
      if (a[j-k]>=aa) goto _40;
      a[j]= a[j-k];
      j= j-k;
      if (j>k) goto _30 ;
_40:
      a[j]= aa;
      }
    }
    
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int median_shell (vector<T>  & a)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n=a.size();
  int i,j,k;
  T aa;
  
  k=0;

  while(k<n) {
    k=3*k+1;
    }
    
  while(k>1) {
    k= k/3;
    if (k<0) return(0);
  
    for(i= k+1;i<n;i++) {
      aa= a[i];
      j= i;
_30:
      if (a[j-k]>=aa) goto _40;
      a[j]= a[j-k];
      j= j-k;
      if (j>k) goto _30 ;
_40:
      a[j]= aa;
      }
    }
    
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> T median_value (const vector<T>  & a)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T median;
  int pos, status;
  vector<T>  b=a;
  
  if(a.size()==1) return(a[0]);
  
  quick_sort(b);
  
  if(a.size()%2==0)
    pos=a.size()/2;
  else
    pos=(a.size()-1)/2;
    
  median=b[pos];
  
  return(median);
}

#endif
