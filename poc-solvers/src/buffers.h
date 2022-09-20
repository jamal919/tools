

/*******************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2012

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France

*******************************************************************************/

#ifndef BUFFERS_H
#define BUFFERS_H

#include "solvers-functions.h"

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  aimed to provide flexible, resizeable memory buffers  
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template<typename T> class buffers_t {
private :
public :
/*------------------------------------------------------------------------------
  geometry parameters*/
  T *x;    
  size_t n;
  
  void init() {
    x=0;
    n=0;
    }
  
  int init(size_t required) {
    this->destroy();
    if(required==0)  return(0);
    if(required==-1) return(0);
    x=new T[required];
    if(x==0) {
      n=0;
      return(-1);
      }
    n=required;
    return(0);
    }
  
  int init(size_t required, T value) {
    int status;
    status=this->init(required);
    if(status!=0) return(-1);
    for(size_t j=0;j<n;j++) x[j]=value;
    return(0);
    }
  
  buffers_t() {
    init();
    }
    
  void destroy() {
    ptr_delete(x);
    n=0;
    }
  
  int check(size_t required) {
    int status=0;
    if(required>n) {
      this->destroy();
      status=this->init(required);
      }
    return(status);
    }
  
//   ~buffers_t() {
//     this->destroy();
//     deletep(&this->x);
//     }
};

class buffers2D_t {
private :
public :
/*------------------------------------------------------------------------------
  geometry parameters*/
  double **x;    
  size_t nrows;
  size_t ncols;
  
  void init() {
    x=0;
    nrows=0;
    ncols=0;
    }
  
  int init(size_t required_rows, size_t required_columns) {
    this->destroy();
    x=new double*[required_rows];
    if(x==0) {
      nrows=0;
      return(-1);
      }
    for(size_t j=0;j<required_rows;j++) {
      x[j]=new double[required_columns];
      if(x[j]==0) {
        nrows=0;
        ncols=0;
        return(-1);
        }
      }
    nrows=required_rows;
    ncols=required_columns;
    return(0);
    }
  
  int init(size_t required_rows, size_t required_columns, double value) {
    int status;
    status=this->init(required_rows, required_columns);
    if(status!=0) return(-1);
    for(size_t j=0;j<nrows;j++) {
      for(size_t i=0;i<ncols;i++) {
        x[j][i]=value;
        }
      }
    return(0);
    }
  
  buffers2D_t() {
    init();
    }
    
  void destroy() {
    if(x!=0) for(size_t j=0;j<nrows;j++) ptr_delete(x[j]);
    ptr_delete(x);
    nrows=0;
    ncols=0;
    }
  
  int check(size_t required_rows, size_t required_columns) {
    int status=0;
    if(required_rows>nrows or required_columns>ncols) {
      this->destroy();
      status=this->init(required_rows, required_columns);
      }
    return(status);
    }
  
  ~buffers2D_t() {
    this->destroy();
    }
};

#endif



