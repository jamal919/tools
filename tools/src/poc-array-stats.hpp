
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief declarations of usefull array statistics functions and classes
*/
/*----------------------------------------------------------------------------*/

#if POC_ARRAY_STATS_HPP == 0
#define POC_ARRAY_STATS_HPP 1

// #include <math.h>
#include <maths.h>

#include "functions.h" //for swapval()
#include "poc-array.hpp"

#include "tools-define.h" //for MIN and MAX


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
min/max and position functions and classes
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename To,typename Ti> bool updatemin(To *m,const Ti & x)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const bool update=*m>x;
  
  if(update)
    *m=x;
  
  return update;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename To,typename Ti> bool updatemax(To *m,const Ti & x)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const bool update=*m<x;
  
  if(update)
    *m=x;
  
  return update;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

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
  T width() {
    T w=(T) 0;
    if(min== (T) +INFINITY and max== (T) -INFINITY) return(w);
    w=max-min;
    return(w);
    }
  void print(const char *comment) {
    if(comment!=0) 
      printf("%s : min=%lf max =%lf\n",comment,min,max);
    else
      printf("min=%lf max =%lf\n",min,max);
    }
};

typedef range_t<int> irange_t;

extern int isfinite(const range_t<double> &r);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename type> range_t<type> minmax(type *x, type mask, int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
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

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename type> range_t<type> minmax(type *x, int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  n;
  type t;
  range_t<type> r;

  r.min=x[0];
  r.max=x[0];

  for(n=1;n<nvalues;n++) {
    t=x[n];
    r.min=MIN(t,r.min);
    r.max=MAX(t,r.max);
    }

  return(r);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> T amin(T *x,int n,int *imin=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///minimum of an array
/**
\param n
\param x array[n] of T
\param *imin optional. Set to index of minimum.
\returns the minimum value, or anything if n<1.
*/
/*----------------------------------------------------------------------------*/
{
  int i,imin_=-1;
  
#if GCC_VERSION >= 40500  /* Test for GCC >=4.5.0 */
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wuninitialized"
#endif
  T xmin;
  
  if(n<1)return xmin;
#if GCC_VERSION >= 40500  /* Test for GCC >=4.5.0 */
  #pragma GCC diagnostic pop
#endif
  
  xmin=x[0];
  for(i=0;i<n;i++) {
    if(updatemin(&xmin,x[i]))
      imin_=i;
    }
  
  if(imin)
    *imin=imin_;

  return xmin;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename type> int array_index(const type *x,const int nvalues,const type & y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**
array_index() and pos() are NOT identical : the argument order is different
*/
{
  int  pos,n;

  pos=-1;
  for(n=0;n<nvalues;n++) {
    if(x[n]==y) {
      pos=n;
      break;
      }
    }

  return(pos);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename type,typename array> int maxpos_template(const array &x, int nvalues, type *xmax=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief return index of max value in x vector
\param x array[nvalues]
\param nvalues
\param *xmax if non \c NULL, will be set to the maximum
\return index of the maximum
*/
/*----------------------------------------------------------------------------*/
{
  int pos,n;
  type xmax_;
  
  if(xmax==NULL)
    xmax=&xmax_;

  *xmax=NAN;
  pos=-1;
  for(n=0;n<nvalues;n++) {
    if((isnan(*xmax) && !isnan(x[n])) || x[n]>*xmax) {
      *xmax=x[n];
      pos=n;
      }
    }

  return(pos);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename type> int maxpos(type *x, int nvalues, type *xmax=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=maxpos_template(x,nvalues,xmax);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename type> int maxpos(vector<type> x, type *xmax=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=maxpos_template(x,x.size(),xmax);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename type,typename array> int minpos_template(const array &x, int nvalues, type *xmin=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief return index of min value in x vector
\param x array[nvalues]
\param nvalues
\param *xmin if non \c NULL, will be set to the minimum
\return index of the minimum
*/
/*----------------------------------------------------------------------------*/
{
  int pos,n;
  type xmin_;
  
  if(xmin==NULL)
    xmin=&xmin_;

  *xmin=NAN;
  pos=-1;
  for(n=0;n<nvalues;n++) {
    if((isnan(*xmin) && !isnan(x[n])) || x[n]<*xmin) {
      *xmin=x[n];
      pos=n;
      }
    }

  return(pos);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename type> int minpos(type *x, int nvalues, type *xmin=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=minpos_template(x,nvalues,xmin);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename type> int minpos(vector<type> x, type *xmin=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=minpos_template(x,x.size(),xmin);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename type,typename array> int vpos_template(const type y, const array &x, int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
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


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename type> int vpos(const type y, const type *x, int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  
  pos=vpos_template(y,x,nvalues);
  
  return(pos);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename type> int vpos(const type y, const std::vector<type> &x)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  
  pos=vpos_template(y,x,x.size());
  
  return(pos);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename type> int vpush(type y, vector<type> & x)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;

  pos=vpos_template(y,x,x.size());
  
  if(pos==-1) {
    x.push_back(y);
    pos=x.size()-1;
    }

  return(pos);
}



/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
sort functions
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> size_t *sort(const T *array, int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l;
  T *x=new T[nvalues];
  size_t *pos=new size_t[nvalues];
  
  valcpy(x,array,nvalues);
  for(k=0;k<nvalues;k++) pos[k]=k;
  
/**-----------------------------------------------------------------------------
  re-arrange x array in ascending order*/
  for(k=0;k<nvalues;k++) {
    l=minpos(&(x[k]),nvalues-k)+k;
    if(k!=l) {
      swapval(x[k],x[l]);
      swapval(pos[k],pos[l]);
      }
    }
  
  delete[] x;
  
  return(pos);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> size_t *sort(vector<T> array)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  size_t *pos;
  
  pos=sort(&array[0],array.size());
  
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
{
  size_t i;
  
  T *b=new T[a.size()];
  for(i=0;i<a.size();i++) b[i]=a[i];
  
  quick_sort(b, a.size());
  
  for(i=0;i<a.size();i++) a[i]=b[i];
  
  delete[]b;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename cardinal, typename type, typename array> cardinal occurence_template(type y, const array x, cardinal nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cardinal  count,n;
  bool nan_mask=isnan(y);
  
  count=0;
  
  switch(nan_mask) {
    case 0:
      for(n=0;n<nvalues;n++) {
        if(x[n]==y) {
          count++;
          }
        }
      break;
    case 1:
      for(n=0;n<nvalues;n++) {
        if(isnan(x[n])) {
          count++;
          }
        }
      break;
    }    

  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename cardinal, typename type> cardinal occurence(type y, const type *x, cardinal nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cardinal count;
  
  count=occurence_template(y,x,nvalues);

  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename type> size_t occurence(type y, const vector<type> & x)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  size_t count;
  
  count=occurence_template(y,x,x.size());

  return(count);
}

#endif
