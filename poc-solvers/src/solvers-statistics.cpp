
/**************************************************************************

  POC-SOLVERS interface, 2006-2019

  Part of the SIROCCO national service (INSU, France)

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France

***************************************************************************/

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 imported from T-UGOm sources (native development sources)
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#include "solvers-statistics.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> statistic_t<T> compute_statistics_template(const T *h, T mask, int n, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k;
  double count;
  double mean,variance,mae,z;
  statistic_t<T> s;

  count=0;
  variance=0;
  mean=mae=0;

  for(k=0;k<n;k++) {
    if(h[k] != mask) {
      z=h[k];
      if(isnan(z)!=0) continue;
      mean=mean+z;
      mae+=fabs(z);
      variance=variance+z*z;
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
    
//   if(verbose==1) print_stats( n, count, mean, variance);

  s.mean=mean;
  s.mae=mae;
  s.std=sqrt(variance);

  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  statistic_t<double> compute_statistics(const double *h, double mask, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  statistic_t<double> s=compute_statistics_template(h, mask, n, 1);
  return(s);
}

