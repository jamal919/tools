
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

#ifndef SOLVERS_FUNCTIONS_H

#define SOLVERS_FUNCTIONS_H

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <string.h>
#include <omp.h>
#include <complex>
#include <vector>
#include <cmath>

#include "poc-assertions.h"
// #include "solvers-statistics.h"

extern void ptr_delete_void(void *&p,void (*F)(void *));

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> void ptr_delete(T* & p,void (*F)(void *))

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  ptr_delete_void((void* &) p, F);
}
 
// template <typename T> void ptr_delete(T* & p, void (*F)(T *)=NULL);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> void ptr_delete(T* & p,void (*F)(T *)=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  delete a pointer and set it to NULL
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  if(p==0) return;
  
  if(F)
    F(p);
  else
    delete[]p;
  p=0;
}

#ifndef _FUNCTIONS_H

template <typename T> class range_t {
private :
public :
  T min,max;
  void init() {
    min= (T) +INFINITY;
    max= (T) -INFINITY;
    }
  range_t() {
    init();
    }
  void init(T min0,T max0) {
    min=min0;
    max=max0;
    }
  range_t(T min0,T max0) {
    init(min0,max0);
    }
  inline range_t<T> & operator <<(const T x){
    updatemin(&min,x);
    updatemax(&max,x);
    return *this;
    }
  void constrain(T *x) const {
    if(min>*x)
      *x=min;
    if(max<*x)
      *x=max;
    }
  void read(const T *p,size_t n) {
    size_t k;
    for(k=0;k<n;k++) {
      *this<<p[k];
      }
    }
  void read(const vector<T> p) {
    size_t k;
    for(k=0;k<p.size();k++) {
      *this<<p[k];
      }
    }
  void read(const T *p, T mask, size_t n) {
    size_t k;
    for(k=0;k<n;k++) {
      if(p[k]!=mask) {
        *this<<p[k];
        }
      }
    }
  void init(const T *p,size_t n) {
    init();
    read(p,n);
    }
  void init(const vector<T> p) {
    init();
    read(p);
    }
  void init(const T *p,T mask, size_t n) {
    init();
    read(p,mask,n);
    }
  range_t(const T *p, T mask, int n) {
    init(p,mask,n);
    }
  range_t(const T *p,int n) {
    init(p,n);
    }
  range_t(const vector<T> p) {
    init(p);
    }
  void narrow(T min0,T max0) {
    updatemax(&min,min0);
    updatemin(&max,max0);
    }
  void init(T center, T radius, T lmin, T lmax) {
    min=center-radius;
    max=center+radius;
    narrow(lmin,lmax);
    }
  bool has(const T x) const {
    return min<=x && x<=max;
    }
  bool in(const T xMin,const T xMax) const {
    return xMin<=min && min<=xMax && xMin<=max && max<=xMax;
    }
  T d() const {
    return max-min;
    }
  void adjust(T margin, T lmin, T lmax) {
    min -= margin;
    max += margin;
    narrow(lmin,lmax);
    }
  void print(const char *comment) {
    if(comment!=0) 
      printf("%s : min=%lf max =%lf\n",comment,min,max);
    else
      printf("min=%lf max =%lf\n",min,max);
    }
};


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename in, typename out, typename card> int recast(in *&x, out *&y, bool duplicate, card nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  if(nvalues==0) return(0);
  
  if(x==0) return(0);

//   if(y==0) y=new out[nvalues];
  y=new out[nvalues];
  
  if(y==0)  return(0);
  
  for(card n=0;n<nvalues;n++) {
    y[n]=(out)x[n];
    }
    
  if(not duplicate) {
    delete[] x;
    x=0;
    }
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename in, typename out> void recast(const in *x, out* &y, int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;

  if(nvalues==0) return;
  
  if(x==0) return;

  if(y==0) y=new out[nvalues];
  
  for(n=0;n<nvalues;n++) {
    y[n]=(out)x[n];
    }
}


extern void deletep_void(void **p,void (*F)(void *));

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> void deletep(T **p,void (*F)(void *))

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  deletep_void((void**)p,F);
  }
 
// template <typename T> void deletep(T **p,void (*F)(T *)=NULL);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> void deletep(T **p,void (*F)(T *)=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  delete a pointer and set it to NULL
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  if(*p==0) return;
  
  if(F)
    F(*p);
  else
    delete[]*p;
  *p=0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> void deletep2D(T ***p,int n,void (*F)(void *)=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// delete a double pointer and set it to \c NULL
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  if(!*p) return;
  
  if(n<=0) {
    printf("dimension n=%d\n",n);
    TRAP_ERR_EXIT(-1,"set the n parameter to the right value in this function call!\n");
    }
  
  if(F==0) {
    for(int j=0; j<n; j++) {
      delete[] (*p)[j];
      }
    delete[] (*p);
    }
  else {
    for(int j=0; j<n; j++) {
      deletep(&(*p)[j],F);
      }
    deletep(p);
    }

}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename type> int vpos(type y, type *x, int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  pos,n;

  pos=-1;
  for(n=0;n<nvalues;n++) {
    if(x[n]==y) {
      pos=n;
      return(pos);
      }
    }

  return(pos);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename type> int vpos(type y, type *x, size_t nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  pos,n;

  pos=-1;
  for(n=0;n<nvalues;n++) {
    if(x[n]==y) {
      pos=n;
      return(pos);
      }
    }

  return(pos);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename type> int vpos(type y, vector<type> x)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  pos,n;

  pos=-1;
  for(n=0;n<x.size();n++) {
    if(x[n]==y) {
      pos=n;
      return(pos);
      }
    }

  return(pos);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename type> int vpos(type y, vector<type> x, int first)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  pos,n;

  pos=-1;
  for(n=first;n<x.size();n++) {
    if(x[n]==y) {
      pos=n;
      return(pos);
      }
    }

  return(pos);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename type> int vpos(type y, vector<type> x, int first, bool condition)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  pos,n;

  pos=-1;
  for(n=first;n<x.size();n++) {
    bool chk=(x[n]==y);
    if(chk==condition) {
      pos=n;
      return(pos);
      }
    }

  return(pos);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> void quick_sort (T *a, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  
  http://www.academia.edu/1976253/Selection_of_Best_Sorting_Algorithm
  
  http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#C
  
------------------------------------------------------------------------------*/
{
    int i, j;
    T p, t;
    if (n < 2)
        return;
    p = a[n / 2];
    for (i = 0, j = n - 1;; i++, j--) {
        while (a[i] < p)
            i++;
        while (p < a[j])
            j--;
        if (i >= j)
            break;
        t = a[i];
        a[i] = a[j];
        a[j] = t;
    }
    quick_sort(a, i);
    quick_sort(a + i, n - i);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> void quick_sort (vector<T> & a)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  
  http://www.academia.edu/1976253/Selection_of_Best_Sorting_Algorithm
  
  http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#C
  
------------------------------------------------------------------------------*/
{
  size_t i;
    
  T *b=new T[a.size()];
  for(i=0;i<a.size();i++) b[i]=a[i];
   
  quick_sort(b, a.size());
  
  for(i=0;i<a.size();i++) a[i]=b[i];
  
  delete[]b;
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename type> int occurence(type y, vector<type> x)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  count,n;

  count=0;
  for(n=0;n<x.size();n++) {
    if(x[n]==y) {
      count++;
      }
    }

  return(count);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename type> int occurence(type y, type *x, int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  count,n;

  count=0;
  for(n=0;n<nvalues;n++) {
    if(x[n]==y) {
      count++;
      }
    }

  return(count);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename type> range_t<type> minmax(type *x, type mask, int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  n,m;
  type t;
  range_t<type> r;

  for(m=0;m<nvalues;m++) {
    if(x[m]==mask) continue;
    r.min=x[m];
    r.max=x[m];
    break;
    }

  for(n=m+1;n<nvalues;n++) {
    t=x[n];
    if(t==mask) continue;
    r.min=MIN(t,r.min);
    r.max=MAX(t,r.max);
    }

  return(r);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename type> range_t<type> minmax(type *x, int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  n;
  type t;
  range_t<type> r;

  r.min=x[0];
  r.max=x[0];

  for(n=1;n<nvalues;n++) {
    t=x[n];
    r.min=min(t,r.min);
    r.max=max(t,r.max);
    }

  return(r);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> bool isnan(const complex<T>& z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  bool result;
  
  result=isnan(real(z)) || isnan(imag(z));
  
  return result;
}


#endif

#endif
