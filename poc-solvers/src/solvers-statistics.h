
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

#ifndef _SOLVER_STATISTIC_H
#define _SOLVER_STATISTIC_H

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <complex>
#include <vector>

using namespace std;


template <typename T> class statistic_t{
private :
public :
  T mean,median,mae;
  double std,rms;
//   double cf,pof,nof,pof2x,nof2x;
  T min,max;
  
//   statistic_t<T> compute_statistics(const T *h, T mask, int n);
  
  statistic_t() {
    mean=median=mae=0;
    std=rms=0;
//     cf=pof=nof=pof2x=nof2x;
    min=max=0;
    }
      
  double variance() {
    double v;
    v=std*std;
    return(v);
    }
  
};

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


extern statistic_t<double> compute_statistics(const double *h, double mask, int n);


#endif
