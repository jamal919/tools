
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief data series statistics
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"
#include "constants.h"

#include "map.h"

extern void fourier(char *save, double *h, double mask, int n, double dt);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<class A,typename T> void poc_minmax_template(A v,int n,T mask,T *min,T *max)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///get min and max of NetCDF values
/**
\param v values
\param n number of values
\param mask
\param *min
\param *max
*/
/*----------------------------------------------------------------------------*/
{
  int i;
  
  *min= (T) +INFINITY;
  *max= (T) -INFINITY;
  
  T *vi;
  
  for(i=0;i<n;i++){
    vi=&v[i];
    if(*vi==mask || isnan(*vi))
      continue;
    updatemin(min,*vi);
    updatemax(max,*vi);
    }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void poc_minmax(double *v,int n,double mask,double *min,double *max)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_minmax_template(v,n,mask,min,max);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void poc_minmax(float *v,int n,float mask,float *min,float *max)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_minmax_template(v,n,mask,min,max);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void poc_minmax(vector<double> v,double mask,double *min,double *max)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_minmax_template(v,v.size(),mask,min,max);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void poc_minmax(vector<float> v,float mask,float *min,float *max)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_minmax_template(v,v.size(),mask,min,max);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  range_t<double> poc_minmax(double *x, int nvalues, double mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  range_t<double> r;

  poc_minmax_template(x, nvalues, mask, &r.min, &r.max);

  return(r);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  range_t<float> poc_minmax(float *x, int nvalues, float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  range_t<float> r;

  poc_minmax_template(x, nvalues, mask, &r.min, &r.max);

  return(r);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  range_t<int> poc_minmax(int *x, int nvalues, int mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  range_t<int> r;

  poc_minmax_template(x, nvalues, mask, &r.min, &r.max);

  return(r);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  range_t<double> poc_minmax(const vector<double> & x, double mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  range_t<double> r;

  poc_minmax_template(x, x.size(), mask, &r.min, &r.max);

  return(r);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  range_t<int> poc_minmax(const vector<int> & x, int mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  range_t<int> r;

  poc_minmax_template(x, x.size(), mask, &r.min, &r.max);

  return(r);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int print_stats(int n,int count,double mean=NAN,double variance=NAN)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int chars;
  
  chars=printf("#samples: %6d %6d",n,count);
  if(count>0)
    chars+=printf(" mean: %8.4lf; variance: %8.4lf; deviation: %8.4lf",
      mean,variance,sqrt(variance));
  chars+=printf("\n");
  
  return chars;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int print_stats(int n,int count,const complex<double> & mean,double variance,char unit)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int chars;
  double argMean=arg(mean);
  
  switch(unit){
    case 'r': break;
    case 'd':
      argMean*=r2d;
      break;
    default:
      TRAP_ERR_EXIT(ENOEXEC,"coding error: unit='%c'\n",unit);
    }
  
  chars=printf("#samples: %6d/%6d; mean: %8.4lf %8.2lf; variance: %8.4lf; deviation: %8.4lf\n",
    n, count, abs(mean), argMean, variance, sqrt(variance));
  
  return chars;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

statistic_t examin_differences(double *h1, double *h2, double mask, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
-----------------------------------------------------------------------*/
{
  int  k,count;
  double mean,variance,z;
  double dmin,dmax;
  int kmin,kmax;
  double *h=NULL;
  statistic_t s;

  count=0;
  variance=0;
  mean=0;

  dmin=+1.0e+10;
  dmax=-1.0e+10;

  h=new double[n];
  for(k=0;k<n;k++) {
    if(h1[k] != mask) {
      h[k]=h1[k]-h2[k];
      mean=mean+h[k];
      if(fabs(h[k])<dmin) {
        dmin=fabs(h[k]);
        kmin=k;
        }
      if(fabs(h[k])>dmax) {
        dmax=fabs(h[k]);
        kmax=k;
        }
      count++;
      }
    else {
      h[k]=mask;
      }
    }

  mean=mean/count;

  for(k=0;k<n;k++) {
    if(h1[k] != mask) {
      h[k]-=mean;
      }
    else {
      h[k]=mask;
      }
    }

  for(k=0;k<n;k++) {
    if(h1[k] != mask) {
      z=h[k];
      variance=variance+z*z;
      }
    }

  mean=mean/count;
  variance=variance/count;
  variance=variance-mean*mean;
  print_stats(n,count,mean,variance);
  printf("#min/max: %d %g  %d %g\n",kmin,dmin,kmax,dmax);

  s.mean=mean;
  s.std=sqrt(variance);

  return(s);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> statistic_t OutliersFrequency_template(T *h, T mask, int n, T limit, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  statistic_t s;
  int  k;
  double count;
  double z;
  double pof,nof,pof2x,nof2x;

  pof=nof=pof2x=nof2x=0;
  count=0;
  for(k=0;k<n;k++) {
    if(h[k] != mask) {
      z=h[k];
      if(z>  limit) pof++;
      if(z< -limit) nof++;
      if(z>  2.*limit) pof2x++;
      if(z< -2.*limit) nof2x++;
      count++;
      }
    }
  if(count==0) {
    s.cf=mask;
    s.nof=mask;
    s.pof=mask;
    s.nof2x=mask;
    s.pof2x=mask;
    return(s);
    }
  
  s.pof=100.*pof/count;
  s.nof=100.*nof/count;
  s.pof2x=100.*pof2x/count;
  s.nof2x=100.*nof2x/count;
  s.cf=100.*(count-pof-nof)/count;

  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  statistic_t OutliersFrequency(double *h, double mask, int n, double limit)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  statistic_t s=OutliersFrequency_template(h, mask, n, limit, 0);
  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> statistic_t get_statistics_template(const T *h, T mask, int n, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k;
  double count;
  double mean,variance,mae,z;
  statistic_t s;

  count=0;
  variance=0;
  mean=mae=0;

  for(k=0;k<n;k++) {
    if(h[k] != mask) {
      z=h[k];
      if(isnan(z)!=0) continue;
      mean+=z;
      mae+=fabs(z);
      variance+=z*z;
      if(isnan(mean)!=0) {
        printf("anomaly in get_statistics_template\n");
        }
      count++;
      }
    }
  
  if(count==0) {
    s.min=mask;
    s.max=mask;
    s.mean=mask;
    s.std=mask;
    print_stats( n, count, mean, variance);
    return(s);
    }
  
  s.min=+INFINITY;
  s.max=-INFINITY;

  mean=mean/count;
  mae/=count;
  variance=variance/count;
  variance=variance-mean*mean;
  
  variance=0;
  for(k=0;k<n;k++) {
    if(h[k] != mask) {
      z=h[k]-mean;
      s.min=min((T) s.min,h[k]);
      s.max=max((T) s.max,h[k]);
      variance=variance+z*z;
      }
    }
  variance=variance/count;
    
  if(verbose==1)
    print_stats( n, count, mean, variance);

  s.mean=mean;
  s.mae=mae;
  s.std=sqrt(variance);

  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> statistic_t get_statistics_template(const vector<T> & h, T mask, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k;
  double count;
  double mean,variance,mae,z;
  statistic_t s;
  
  int n=h.size();

  count=0;
  variance=0;
  mean=mae=0;

  for(k=0;k<n;k++) {
    if(h[k] != mask) {
      z=h[k];
      mean=mean+z;
      mae+=fabs(z);
      variance=variance+z*z;
      count++;
      }
    }
  
  if(count==0) {
    s.min=mask;
    s.max=mask;
    s.mean=mask;
    s.std=mask;
    return(s);
    }
  
  s.min=+INFINITY;
  s.max=-INFINITY;

  mean=mean/count;
  mae/=count;
  variance=variance/count;
  variance=variance-mean*mean;
  
  variance=0;
  for(k=0;k<n;k++) {
    if(h[k] != mask) {
      z=h[k]-mean;
      s.min=min(s.min,(double) h[k]);
      s.max=max(s.max,(double) h[k]);
      variance=variance+z*z;
      }
    }
  variance=variance/count;
    
  if(verbose==1)
    print_stats( n, count, mean, variance);

  s.mean=mean;
  s.mae=mae;
  s.std=sqrt(variance);

  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> statistic_t get_statistics_template(const T *h, T mask, int *list, int n, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k,kk;
  double count;
  double mean,variance,z;
  statistic_t s;

  count=0;
  variance=0;
  mean=0;

  for(kk=0;kk<n;kk++) {
    k=list[kk];
    if(h[k] != mask) {
      z=h[k];
      mean=mean+z;
      variance=variance+z*z;
      count++;
      }
    }
  
  if(count==0) {
    s.min=mask;
    s.max=mask;
    s.mean=mask;
    s.std=mask;
    print_stats( n, count);
    return(s);
    }
  
  mean=mean/count;
  variance=variance/count;
  variance=variance-mean*mean;
  
  variance=0;
  for(kk=0;kk<n;kk++) {
    k=list[kk];
    if(h[k] != mask) {
      z=h[k]-mean;
      variance=variance+z*z;
      }
    }
  variance=variance/count;
    
  if(verbose==1)
    print_stats( n, count, mean, variance);

  s.mean=mean;
  s.std=sqrt(variance);

  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  statistic_t get_statistics(const int *h, int mask, size_t n, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  statistic_t s=get_statistics_template(h, mask, n, verbose);
  return(s);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  statistic_t get_statistics(const double *h, double mask, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  statistic_t s=get_statistics_template(h, mask, n, 1);
  return(s);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  statistic_t get_statistics(const vector<double> & h, double mask, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  statistic_t s=get_statistics_template(h, mask, verbose);
  return(s);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  statistic_t get_statistics(const vector<float> & h, float mask, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  statistic_t s=get_statistics_template(h, mask, verbose);
  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  statistic_t get_statistics(const double *h, double mask, int n, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  statistic_t s=get_statistics_template(h, mask, n, verbose);
  return(s);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  statistic_t get_statistics(const double *h, double mask, range_t<size_t> range, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  size_t m,n;
  m=range.min;
  n=range.width();
  
  statistic_t s=get_statistics_template(&h[m], mask, n, 1);
  return(s);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  statistic_t get_statistics(const float *h, float mask, int n, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  statistic_t s=get_statistics_template(h, mask, n, verbose);
  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  statistic_t get_statistics(const float *h, float mask, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  statistic_t s=get_statistics_template(h, mask, n, 1);
  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  statistic_t get_statistics(const float *h, float mask, int *list, int n, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  statistic_t s=get_statistics_template(h,  mask, list, n, verbose);
  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

statistic_c_t cget_statistics(fcomplex *h, fcomplex mask, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k,count;
  double variance;
  complex<double> mean,z;
  statistic_c_t s;

  count=0;
  mean=0;
  for(k=0;k<n;k++) {
    if(h[k] != mask) {
      count++;
      }
    }

  variance=0;
  for(k=0;k<n;k++) {
    if(h[k] != mask) {
      z=(complex<double>) h[k]-mean;
      variance=variance+norm(z)/2;
      }
    }

  variance=variance/count;

  print_stats(n, count, mean ,variance);

  mean=0;
  count=0;
  for(k=0;k<n;k++) {
    if(h[k] != mask) {
      z=(complex<double>) h[k];
      mean=mean+z;
      count++;
      }
    }

  mean=mean/complex<double> (count);

  variance=0;
  for(k=0;k<n;k++) {
    if(h[k] != mask) {
      z=(complex<double>) h[k]-mean;
      variance=variance+norm(z)/2;
      }
    }

  variance=variance/count;

  print_stats(n, count, mean ,variance);

  s.mean=mean;
  s.std=sqrt(variance);

  return(s);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  statistic_c_t get_statistics(complex<double> *h, complex<double> mask, int n, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k;
  double count;
  double variance,v;
  
  complex<double> mean,z, mae;
  statistic_c_t s;
  
  count=0;
  variance=0;
  mean=mae=0;

  for(k=0;k<n;k++) {
    if(h[k] != mask) {
      z=h[k];
      mean=mean+z;
      mae+=abs(z);
      count++;
      }
    }
  
  if(count==0) {
    s.min=mask;
    s.max=mask;
    s.mean=mask;
    s.std=1./0.;
    return(s);
    }
  
  mean=mean/count;
  mae/=count;
    
  variance=0;
  for(k=0;k<n;k++) {
    if(h[k] != mask) {
      v=abs(h[k]-mean);
      variance=variance+v*v;
      }
    }
  variance=variance/count;
    
  if(verbose==1)
    printf("#samples: total=%6d valid=%6d; mean: %g; variance: %g; deviation: %g\n",
          n,(int) count,abs(mean),variance,sqrt(variance));

  s.mean=mean;
//   s.mae=mae;
  s.std=sqrt(variance);

  return(s);
}



// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// statistic_t sget_statistics01(float *h, float mask, int n, int verbose)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /*------------------------------------------------------------------------------
// -----------------------------------------------------------------------*/
// {
//   int  k,count;
//   double mean,variance,z;
//   statistic_t s;
// 
//   count=0;
//   variance=0;
//   mean=0;
// 
//   for(k=0;k<n;k++) {
//     if(h[k] != mask) {
//       z=h[k];
//       mean=mean+z;
//       variance=variance+z*z;
//       count++;
//       }
//     }
//   mean=mean/count;
//   variance=variance/count;
//   variance=variance-mean*mean;
//   if(verbose)
//     printf("##: %6d; mean: %8.4lf; var: %8.4lf; std: %8.4lf\n",
//           count,mean,variance,sqrt(variance));
// 
//   s.mean=mean;
//   s.std=sqrt(variance);
// 
//   return(s);
// 
// }



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename type> statistic_t get_statistics_template(const type *h, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *-----------------------------------------------------------------------
  no mask value
-----------------------------------------------------------------------*/
{
  unsigned int  k,count;
  double mean,variance,z;
  statistic_t s;

  count=0;
  variance=0;
  mean=0;

  for(k=0;k<n;k++) {
    z=h[k];
    mean=mean+z;
    variance=variance+z*z;
    count++;
    }
  mean=mean/count;
  variance=variance/count;
  variance=variance-mean*mean;
  print_stats( n, count, mean, variance);

  s.mean=mean;
  s.std=sqrt(variance);

  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  statistic_t get_statistics(const double *h, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  statistic_t s=get_statistics_template(h,  n);
  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  statistic_t get_statistics(const float *h, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  statistic_t s=get_statistics_template(h,  n);
  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

statistic_t sget_filteredstatistics(float *h, float mask, int n, float threshold)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k;
  double variance;
  float *z=NULL;
  statistic_t s;

  s=get_statistics(h,  mask, n);

  variance= s.std*s.std;

  exitIfNull(
    z=(float *) malloc(n*sizeof(float))
    );
  for(k=0;k<n;k++) {
    if(h[k] != mask) {
      z[k]=(h[k]-s.mean);
      if(z[k]*z[k] > threshold*variance ) z[k]=mask;
      }
    else  z[k]=mask;
    }

  s=get_statistics(z,  mask, n);

  free(z);

  return(s);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

statistic_t sget_filteredstatistics2(float *h1,float *h2, float mask, int n, float threshold)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k;
  double variance1;
  double variance2;
  float *z1=NULL,*z2=NULL;
  statistic_t s1,s2;

  s1=get_statistics(h1,  mask, n);
  s2=get_statistics(h2,  mask, n);

  variance1= s1.std*s1.std;
  variance2= s2.std*s2.std;

  exitIfNull(
    z1=(float *) malloc(n*sizeof(float))
    );
  exitIfNull(
    z2=(float *) malloc(n*sizeof(float))
    );
  for(k=0;k<n;k++) {
    if(h1[k] != mask) {
      z1[k]=(h1[k]-s1.mean);
      if(z1[k]*z1[k] > threshold*variance1 ) z1[k]=mask;
      }
    else  z1[k]=mask;
    if(h2[k] != mask) {
      z2[k]=(h2[k]-s2.mean);
      if(z2[k]*z2[k] > threshold*variance2 ) z2[k]=mask;
      }
    else  z2[k]=mask;
    }

  for(k=0;k<n;k++) {
    if((z1[k] != mask) && (z2[k] != mask)) {
      continue;
      }
    else {
      z1[k]=mask;
      z2[k]=mask;
      }
    }

  s1=get_statistics(z1,  mask, n);
  s2=get_statistics(z2,  mask, n);

  free(z1);
  free(z2);

  return(s1);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  regression_t get_trend(float *h, float *t, float mask, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k,count;
  double covariance,z;
  statistic_t s1,s2;
  regression_t r;

  s1=get_statistics(h, mask, n,0);
  s2=get_statistics(t, mask, n,0);

  count=0;
  covariance=0;

  for(k=0;k<n;k++) {
    if(h[k] != mask) {
      z=(h[k]-s1.mean)*(t[k]-s2.mean);
      covariance=covariance+z;
      count++;
      }
    }

  covariance=covariance/count;

  r.trend=covariance/(s2.std*s2.std);
  r.bias=s1.mean-r.trend*s2.mean;

  printf("-------------------------------\n");
  printf("trend= %6.2lf   bias= %6.2lf \n",r.trend,r.bias);

  return(r);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> covariance_t get_covariance_template(T *h1, T *h2, T mask, int n, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k,count;
  double covariance,correlation;
  T z;
  statistic_t s1,s2;
  covariance_t c;
  
  s1=get_statistics(h1, mask, n, 0);
  s2=get_statistics(h2, mask, n, 0);

  count=0;
  covariance=0;

  for(k=0;k<n;k++) {
    if((h1[k] != mask)  && (h2[k] != mask)) {
      z=(h1[k]-s1.mean)*(h2[k]-s2.mean);
      covariance=covariance+z;
      count++;
      }
    }

  if(count==0) return(c);
  
  covariance=covariance/count;
  c.covariance=covariance;

  correlation=covariance/(s2.std*s1.std);

  c.correlation=correlation;
  c.s1=s1;
  c.s2=s2;
  
  c.n=count;

  if(verbose==1) {
    printf("-------------------------------\n");
    printf("covariance= %6.2lf   correlation= %6.2lf \n",covariance,c.correlation);
    }
  
  return(c);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  covariance_t get_covariance(float *h1, float *h2, float mask, int n, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  covariance_t c;
  c=get_covariance_template(h1, h2, mask, n, verbose);
  return(c);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  covariance_t get_covariance(double *h1, double *h2, double mask, int n, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  covariance_t c;
  c=get_covariance_template(h1, h2, mask, n, verbose);
  return(c);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *get_autocorrelation(double *h,int n,int range, int test)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k,m,count;
  double covariance,*correlation,z,v;
//   double d;
  double *tmp=NULL;
  statistic_t s;
  covariance_t c;

  s=get_statistics(h, n);
  v=s.std*s.std;

  correlation=new double[range];

  for(m=0;m<range;m++) {
/* *-----------------------------------------------------------------------
    in theory, should be treated as a call to covariance function...*/
    count=0;
    covariance=0;
    for(k=0;k<n-m;k++) {
//      z=(h[k]-s.mean)*(h[k+m]-s.mean)/v;
      z=(h[k]-s.mean)*(h[k+m]-s.mean);
      covariance+=z;
      count++;
      }
    covariance=covariance/count;
    correlation[m]=covariance/(s.std*s.std);
    }

//  fourier((char *) 0, correlation, correlation, 1.e+10, range, 1.0);

//  s=get_statistics(correlation, m);
/* *-----------------------------------------------------------------------
  sucessive data are incremented in cycles:

  we seek peaks in correlation (local maximum)

----------------------------------------------------------------------*/
  for(m=1;m<range-1;m++) {
//    d=fabs(correlation[m]-s.mean)/s.std;
    if((correlation[m]>correlation[m-1]) && (correlation[m]>correlation[m+1])) {
      printf("%d %lf\n",m,correlation[m]);
      if(test%m==0) {
        printf("%d %lf\n",m,correlation[m]);
        return(correlation);
        }
      }
    }

/* *-----------------------------------------------------------------------
  sucessive data are incremented by cycles

  we seek slow decrease in correlation

----------------------------------------------------------------------*/
  tmp=new double[n];
  for(k=0;k<n-1;k++) {
    tmp[k]=h[k+1]-h[k];
    }
  tmp[n-1]=h[0]-h[n-1];

  s=get_statistics(tmp, n);

  for(m=0;m<range;m++) {
/* *-----------------------------------------------------------------------
    in theory, should be treated as a call to covariance function...*/
    count=0;
    covariance=0;
    for(k=0;k<n-m-1;k++) {
//      z=(h[k]-s.mean)*(h[k+m]-s.mean)/v;
      z=(tmp[k]-s.mean)*(tmp[k+m]-s.mean);
      covariance+=z;
      count++;
      }
    covariance=covariance/count;
    correlation[m]=covariance/(s.std*s.std);
    }

  for(m=1;m<range-1;m++) {
//    d=fabs(correlation[m]-s.mean)/s.std;
    if((correlation[m]>correlation[m-1]) && (correlation[m]>correlation[m+1])) {
      printf("%d %lf\n",m,correlation[m]);
      if(test%m==0) {
        printf("%d %lf\n",m,correlation[m]);
        return(correlation);
        }
      }
    }

  delete[] tmp;
  for(m=1;m<range-1;m++) {
    }

  return(correlation);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *get_maxcorrelation(double *h1,double *h2,int n,int range, int test, int & max_positive, int & max_negative, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k,m,count;
  double covariance,*correlation,z;
//   double *tmp=NULL;
  statistic_t s1,s2;
  covariance_t c;
  double maxp=-1.e+10,maxn=1.e+10;
  double mask=1.e+10;

  s1=get_statistics(h1, mask, n, verbose);
  s2=get_statistics(h2, mask, n, verbose);

  correlation=new double[range];

  for(m=0;m<range;m++) {
/* *-----------------------------------------------------------------------
    in theory, should be treated as a call to covariance function...*/
    count=0;
    covariance=0;
    for(k=0;k<n-m;k++) {
      z=(h1[k]-s1.mean)*(h2[k+m]-s2.mean);
      covariance+=z;
      count++;
      }
    covariance=covariance/count;
    correlation[m]=covariance/(s1.std*s2.std);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  sucessive data are incremented in cycles:

  we seek peaks in correlation (local maximum)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(m=1;m<range-1;m++) {
    if((correlation[m]>correlation[m-1]) && (correlation[m]>correlation[m+1])) {
      if(verbose==1) printf("%d %lf\n",m,correlation[m]);
      if(maxp<correlation[m]){
        maxp=correlation[m];
        max_positive=m;
        }
      }
    if((correlation[m]<correlation[m-1]) && (correlation[m]<correlation[m+1])) {
      if(verbose==1) printf("%d %lf\n",m,correlation[m]);
      if(maxn>correlation[m]){
        maxn=correlation[m];
        max_negative=m;
        }
      }
    }

// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
// 
//   sucessive data are incremented by cycles
// 
//   we seek slow decrease in correlation
// 
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// 
//   tmp=new double[n];
//   for(k=0;k<n-1;k++) {
//     tmp[k]=h[k+1]-h[k];
//     }
//   tmp[n-1]=h[0]-h[n-1];
// 
//   s=get_statistics(tmp, n);
// 
//   for(m=0;m<range;m++) {
// /* *-----------------------------------------------------------------------
//     in theory, should be treated as a call to covariance function...*/
//     count=0;
//     covariance=0;
//     for(k=0;k<n-m-1;k++) {
//       z=(tmp[k]-s.mean)*(tmp[k+m]-s.mean);
//       covariance+=z;
//       count++;
//       }
//     covariance=covariance/count;
//     correlation[m]=covariance/(s.std*s.std);
//     }
// 
//   for(m=1;m<range-1;m++) {
// //    d=fabs(correlation[m]-s.mean)/s.std;
//     if((correlation[m]>correlation[m-1]) && (correlation[m]>correlation[m+1])) {
//       printf("%d %lf\n",m,correlation[m]);
//       if(test%m==0) {
//         printf("%d %lf\n",m,correlation[m]);
//         return(correlation);
//         }
//       }
//     }
// 
//   delete[] tmp;
//   for(m=1;m<range-1;m++) {
//     }

  return(correlation);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double get_timesampling(const double *t, double mask, int n, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double sampling;
  statistic_t s;
  vector<double> continuous;
  size_t niterations=0;
  double factor=3;
  
  if(n<=0) TRAP_ERR_RETURN(NAN,verbose>0,"%d<=0 samples\n",n);
  
  double *delta=aset(n-1,mask);
  
  for(size_t k=0;k<n-1;k++) {
    delta[k]=mask;
    if(t[k]==mask)   continue;
    if(t[k+1]==mask) continue;
    delta[k]=t[k+1]-t[k];
    }
  
  s=get_statistics(delta, mask, n-1, verbose);
  
  sampling=s.mean;
  
  while(s.std!=0.0) {
    int remasked=0;
    if(verbose==1) printf("%s:masking all dt out of (mean=%g)+/-(std=%g)*(factor=%g): ",__func__,s.mean,s.std,factor);
    for(size_t k=0;k<n-1;k++) {
      if(delta[k]==mask)   continue;
      if(fabs(delta[k]-s.mean) > factor*s.std) {
        delta[k]=mask;
        remasked++;
        }
      }
    
    if(verbose==1) printf("%d remasked: ",remasked);
    s=get_statistics(delta, mask, n-1, verbose);
    if(verbose==1) printf("dt in [%g;%g]\n",s.min,s.max);
    sampling=s.mean;
//     if(s.std==0.)
    if(s.max-s.min<1e-6*sampling)
      break;
    if(remasked==0) factor*=0.9;
    niterations++;
    if(niterations>1000000) {
      printf("get_timesampling : could not converge after %d iterations, abort\n",niterations);
      sampling=0.0;
      return(sampling);
      }
    }
    
  for(size_t k=0;k<n-1;k++) {
    if(t[k]==mask)   continue;
    if(t[k+1]==mask) continue;
    delta[k]=t[k+1]-t[k];
    if(fabs(delta[k]) > 1.25*s.mean) {
      delta[k]=mask;
      }
    }
  
  s=get_statistics(delta, mask, n-1, verbose);
  
  
  delete[] delta;
  return(sampling);
}
