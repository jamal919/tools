
#ifndef _IMBRICATION_H
#define _IMBRICATION_H

#include "fe.def"
#include "fe-proto.h"

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  imbrication level 0 class
 
  basic_t provides interpolation information nucleus needed to perform a field
  interpolation at a given position

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

class basic_t{
private :
public :
  int element;                   /* element index                                         */
  mesh_t  *mesh;                 /* pointer to mesh                                       */
  int *nodes, nnodes;            /* nodes index and number of numerical nodes in elements */
  double *beta;                  /* interpolation weights                                 */
  double lon,lat,x,y;            /* positions (spherical, cartesian)                      */
  int discretisation;            /* element discretisation flag                           */
  discretisation_t *descriptor;  /* element discretisation structure                      */
  
  basic_t() {
    element=-1;
    lon=NAN;
    lat=NAN;
    discretisation=LGP1;
    nnodes=3;
    nodes=new int[nnodes];
    beta=new double[nnodes];
    for(int k=0;k<nnodes;k++) {
      nodes[k]=-1;
      }
    }
    
  basic_t(int n) {
    element=-1;
    lon=NAN;
    lat=NAN;
    discretisation=LGP1;
    nnodes=n;
    nodes=new int[nnodes];
    beta=new double[nnodes];
    for(int k=0;k<nnodes;k++) {
      nodes[k]=-1;
      }
    }

  void init(discretisation_t *d) {
    element=-1;
    lon=NAN;
    lat=NAN;
    discretisation=d->type;
    nnodes=d->nnpe;
    nodes=new int[nnodes];
    beta=new double[nnodes];
    for(int k=0;k<nnodes;k++) {
      nodes[k]=-1;
      }
    descriptor=d;
    }
    
  basic_t(discretisation_t *d) {
    init(d);
    }
    
  int check(double xx, double yy) {
    int status;
//     if(xx==lon and yy==lat) return(0);
    int hint=element;
    status = fe_beta(*mesh, xx, yy, hint, &element, nodes, beta);
    lon=xx;
    lat=yy;
    return(status);
    }
    
  template <typename T> T interpolate(T *buffer) {
    T z=0;
    for(int k=0;k<nnodes;k++) {
      int n=nodes[k];
      z+=buffer[n]*beta[k];
      }
    return(z);
    }
    
  template <typename T> T interpolate(T *buffer, T mask, T substitute) {
    T z=0;
    if(element==-1) {
      return(substitute);
      }
    for(int k=0;k<nnodes;k++) {
      int n=nodes[k];
      T tmp=buffer[n];
      if(tmp==mask) tmp=substitute;
      z+=tmp*beta[k];
      }
    return(z);
    }
    
  ~basic_t() {
    element=-1;
    discretisation=-1;
    nnodes=0;
    delete[] nodes;
    delete[] beta;
    }
};

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  imbrication level 1 class
 
  nesting_t provides interpolation information needed to perform a field
  interpolation at a given array of positions

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <typename T> class nesting_t {
private :
public :
  basic_t *basics;               /* array of interpolation weights and nodes index  */
  mesh_t  *mesh;                 /* pointer to mesh                                 */
  T* buffer;                     /* temporary buffer                                */
  int nnodes;                    /* number of nodes in array                        */
  discretisation_t *descriptor;  /* discretisation descriptor of external UG field  */
     
  int allocate(int s, discretisation_t *d, mesh_t *m) {
    mesh=0;
    nnodes=s;
    descriptor=d;
    basics=new basic_t[nnodes];
    mesh=m;
    for(int n=0;n<nnodes;n++) {
      basics[n].nodes =new int[d->nnpe];
      basics[n].beta  =new double[d->nnpe];
      basics[n].mesh  =m;
      }
    }
    
  nesting_t() {
    mesh=0;
    descriptor=0;
    nnodes=0;
    }
    
  nesting_t(int s, discretisation_t *d, mesh_t *m) {
    mesh=0;
    nnodes=s;
    descriptor=d;
    basics=new basic_t[nnodes];
    mesh=m;
    for(int n=0;n<nnodes;n++) {
      basics[n].nodes =new int[d->nnpe];
      basics[n].beta  =new double[d->nnpe];
      }
    }
      
  int check(vector<double> x, vector<double> y) {
    int status;
    for(int n=0;n<nnodes;n++) {
      status=basics[n].check(x[n],y[n]);
      }    
    }
    
  ~nesting_t() {
    mesh=0;
    descriptor=0;
    nnodes=0;
    }
};

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  imbrication level 2 class
  
  hold interpolation material (nesting) and value at positions
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <class T> class imbrication_t {
private :
public :
  nesting_t<double> *nesting;                     /* cells and related nodes for positions  */
  T *x;                                           /* field values at positions              */
  char * units;                                   /*                                        */
  T mask;                                         /*                                        */
  double time,next;                               /*                                        */
  int nvalues ;                                   /*                                        */
  int discretisation;                             /*                                        */
  int initialized;                                /*                                        */
  int volatile_data;                              /*                                        */
  int (*pre_processing)(void);                    /*                                        */
  int (*post_processing)(imbrication_t<T>);       /*                                        */
  int (*processor) (imbrication_t<T> &, double);  /*                                        */
  datastream_t< T > *stream;
  
  void init() {
    initialized=0;
    nesting=0;
    x=0;
    units=0;
    nvalues=0;
    volatile_data=0;
    pre_processing=0;
    post_processing=0;
    processor=0;
    stream=0;
    time=-1.e+10;
    next=-1.e+10;
    }
    
  imbrication_t() { 
    init();
    }
    
  int locate(vector<double> x, vector<double> y) {
    int status;
    status=nesting->check(x,y);
    return(status);
    }

  int check(double t) {
    int status;
    if(processor==0) return(-1);
    if(t > time) status=processor(*this, t);
    return(status);
    }

  int acquire(double t) {
    int status;
    if(processor==0) return(-1);
    status=processor(*this, t);
    if(post_processing!=0) post_processing(*this);
    return(status);
    }

  int allocate(int n) {
    nvalues=n;
    x=new T[nvalues];
//    nesting=new nesting_t<double>[nvalues];
    return(0);
    }

  imbrication_t(int size, nesting_t<double> *n, datastream_t < T > *s , int (*p) (imbrication_t<T> &, double)  ) {
    init();
    this->allocate(size);
    nesting=n;
    stream=s;
    processor=p;
    }
    
};
extern int stream_imbricationxUG(imbrication_t<double> & Ufield, double t);
extern int stream_imbricationxUG_02(imbrication_t<double> & Ufield, double t);

#endif
