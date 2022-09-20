
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief declarations of usefull array convertion and allocation functions
*/
/*----------------------------------------------------------------------------*/

#if POC_ARRAY_HPP == 0
#define POC_ARRAY_HPP 1

#include <string.h>

#include <vector>
#include <errno.h>

#include "poc-assertions.h"

using namespace std;


/*----------------------------------------------------------------------------*/
/// gives a copy of an array
/**
by allocating a new array and initialising the values
\param s source
\param n number of elements
\return the array
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> T *copy(const T *s,int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T *r;//returned value
  
  r=new T[n];
  
  memcpy(r,s,n*sizeof(T));
  
  return r;
}

extern void deletep_void_void(void **p,void (*F)(void *));

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> void deletep_void(T **p,void (*F)(void *))

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  deletep_void_void((void**)p,F);
}

template <typename T> void deletep(T **p,void (*F)(T *)=NULL);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename Td,class A> void valcpy(Td * dest,const A & src, int size)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// copy values of an array to another, converting if needed
/*----------------------------------------------------------------------------*/
{
  exitIfNull(dest);
  
  for(size_t k=0;k<size;k++)
    dest[k]=(Td)src[k];
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename Td,typename Ts> void poc_copy(Td * & v,const Ts *src, int size)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// copy an array to another, RE-allocating the later
/*----------------------------------------------------------------------------*/
{
  deletep(&v);
  
  if(!src) return;
  
  v=new Td[size];
  
  valcpy(v,src,size);
}


/*----------------------------------------------------------------------------*/
/// gives a copy of an array
/**
by allocating a new array and initialising the values
\param s source
\return the array
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> T *copy(const vector<T> &s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T *r;//returned value
  const int n=s.size();
  
  r=new T[n];
  
  for(int i=0;i<n;i++)
    r[i]=s[i];
  
  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> void copy(vector<T> *v,const T* src, int size)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// array to vector
/*----------------------------------------------------------------------------*/
{
  v->resize(size);
  
  for(int i=0;i<size;i++)
    (*v)[i]=src[i];
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename in> double *recast(in *x, int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  double *buffer=0;
  
  if(nvalues==0) return(0);
  
  if(x==0) return(buffer);

  buffer=new double[nvalues];
  for(n=0;n<nvalues;n++) {
    buffer[n]=x[n];
    }

  return(buffer);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename in, typename out> void recast(const in *x, out* &y, int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n;

  if(nvalues==0) return;
  
  if(x==0) return;

  if(y==0) y=new out[nvalues];
  
  for(n=0;n<nvalues;n++) {
    y[n]=(out)x[n];
    }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> void aset(T *a,size_t n,const T & v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// sets all the values of an array to a value
/**
\param a array
\param n number of elements
\param v value
*/
/*----------------------------------------------------------------------------*/
{
  int i;
  for(i=0;i<n;i++){
    a[i]=v;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<class A,typename T> void aset(A *a,const T & v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// sets all the values of an array to a value
/**
\param a pointer to array object
\param v value
*/
/*----------------------------------------------------------------------------*/
{
  const int n=a->size();
  
  for(int i=0;i<n;i++){
    (*a)[i]=v;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> T *aset(size_t n,const T & v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// allocate an array and sets all its values
/**
\param n number of elements
\param v value
*/
/*----------------------------------------------------------------------------*/
{
  T *r;
  r=new T[n];
  
  if(r!=0)
    aset(r,n,v);
  
  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> void deletep(T **p,void (*F)(T *))

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// delete a pointer and set it to \c NULL
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  if(*p==0) return;
  
  if(F)
    F(*p);
  else
    delete[]*p;
  *p=0;
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> void newp2D(T ***p,int m,int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// allocated a 2D pointer
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  *p=new T*[m];
  
  for(int j=0; j<m; j++){
    (*p)[j]=new T[n];
    }
  /** \endcode */
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
  
  if(n<0)
    TRAP_ERR_EXIT(ENOEXEC,"set the n parameter to the right value in this function call!\n");
  
  for(int j=0; j<n; j++){
    if(F!=0)
      deletep_void(&(*p)[j],F);
    else
      deletep(&(*p)[j]);
    }
  
  deletep(p);
  /** \endcode */
}


#endif
