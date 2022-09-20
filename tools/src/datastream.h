
/*******************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/

#ifndef DATASTREAM_H
#define DATASTREAM_H

#include <limits>
#include "functions.h"
#include "poc-netcdf.hpp"
#include "netcdf-classes.h"
#include "map-classes.h"

#define LINEAR    0
#define QUADRATIC 1

extern void print_stream_DecodeName_help(int cpu=-1);
extern int stream_DecodeName(double t,const string & filedir,const string & file_convention, char **filename, double *dt=0);
extern int stream_DecodeName(double t, int cpu,const string & filedir,const string & file_convention, char **filename, double *dt=0);


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  basic container
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <class T> class datastream_t {
private :
public :
  T *dum;
  void *initiator;
  double __time,__next;              /*                                  */
  void* auxiliary;                   /*                                  */
  
  virtual int check (double t) {
    return (0);
    }
  virtual double time () {
    return (0);
    }

};


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  stream url, network access
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <class T> class urlstream_t: public  datastream_t <T>  {
private :
public :
  string url;
  string username, password;
  string protocol;
  int (*processor) (urlstream_t<T> &, double);   /*                                  */
  datastream_t<T> *stream;                       /*                                  */
  int id;
  int (*identify)(const char *, const char *, int *, int);
  T dum;
  
  urlstream_t() {
    url="";
    identify=0;
    id=-1;
    }
  
  };

  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  stream file, names for file and id for variable
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

extern void print_stream_download_and_cleanup_help();

template <class T> class filestream_t;

extern void stream_download_and_cleanup(filestream_t<float > * filestream, double t, double dt);
extern void stream_download_and_cleanup(filestream_t<double> * filestream, double t, double dt);
extern void stream_download_and_cleanup(filestream_t<short > * filestream, double t, double dt);

template <class T> class filestream_t: public  datastream_t <T>  {
private :
public :
  string filename;
  string path,convention,varname;
  poc_global_t glob;
  string time_format, separator;
  int (*processor) (filestream_t<T> &, double);  /* could be used for wget/rsync download? */
  datastream_t<T> *stream;                       /* could be used for wget/rsync download? */
  int id;
  int (*identify)(const char *, const char *, int *, int);
  T dum;
  
  filestream_t() {
    filename=path=convention="";
    path="";
    convention="";
    varname="";
    identify=0;
    id=-1;
    }
  
  filestream_t(const char *p, const char *c, const char *v, int (*func)(const char *, const char *, int *, int)) {
    path=(string) p;
    convention=(string) c;
    if(v!=0) varname=(string) v;
    identify=func;
    id=-1;
    }
  
  int get_id(int verbose=0) {
    int status;
    int varid;
    
/**----------------------------------------------------------------------------
    identify variable in file */
    if(identify!=0) {
      status=this->identify(this->filename.c_str(),varname.c_str(),&varid,verbose);
      this->id=varid;
      }
    else {
      status=poc_inq(filename,&glob,verbose);
      if(status!=0)
        return status;
      id=glob.variables.find(varname);
      if(id>=0)
        status=0;
      else
        status=NC_ENOTVAR;
      }
    
    return(status);
    }
  
  int check (double t) {
/*------------------------------------------------------------------------------
    check will simply update stream value (filename) */
    char *tmp;
    int status;
    double dt;
    
/**----------------------------------------------------------------------------
    construct filename as a combination of convention and time */
    const char *tmpdir=getenv("TMPDIR");
    
    if(tmpdir==0)
      status = stream_DecodeName(t,this->path, this->convention, &tmp);
    else{
      status = stream_DecodeName(t,tmpdir    , this->convention, &tmp, &dt);
      stream_download_and_cleanup(this,t,dt);
      }
    
    if(this->filename!=tmp){
      this->filename=tmp;
      status=this->get_id();
      }
    else
      status=0;
    
    delete[]tmp;
    
    return status;
    }
  
  int finalize (double t,int verbose=0) {
/*------------------------------------------------------------------------------
    finalize will update stream value (filename) and inquire variable id */
    int status;
    
    status = check(t);
    
/**----------------------------------------------------------------------------
    identify variable in file */
    status = get_id(verbose);
    
    return(status);
    }
  
  int finalize (double t, const char *v) {
/*------------------------------------------------------------------------------
    wrapper */
    int status;
    varname=(string) v;
    status=finalize (t);
    return(status);
    }
};


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  stream for data handling, structured
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <class T> class SGfield_t : public  datastream_t <T> {
private :
public :
  T *x, mask;                                 /* field values at nodes            */
  char * units;                               /* field units                      */
  double time,next;                           /*                                  */
  grid_t *grid;                               /* pointer to related grid          */
//   grib_gaussian_t *gaussian;                  /* pointer to related grid          */
  int gtype;                                  /*                                  */
  int volatile_data;                          /*                                  */
  int (*pre_processing)(void);                /*                                  */
  int (*post_processing)(SGfield_t<T> &);     /*                                  */
  int (*processor) (SGfield_t<T> &, double);  /*                                  */
  datastream_t<T> *stream;                    /*                                  */
  
  void destroy() {
    deletep(&x);
    }
  
  void init() {
    x=0;
    units=0;
    time=-INFINITY;
    grid=0;
//     gaussian=0;
    gtype=0;
    volatile_data=0;
    pre_processing=0;
    post_processing=0;
    processor=0;
    stream=0;
    }
  
  SGfield_t() {
    this->init();
    }
  
  int check(double t) {
/*------------------------------------------------------------------------------
    check will conditionnally update stream value(s) */
    int status=0;
    if(processor==0)
      return(-1);
/*------------------------------------------------------------------------------
    originally meant for ascending time request, do it mor general */
//     if(t > time) status=processor(*this, t);
    if(t != time)
      status=processor(*this, t);
    return status;
    }

  int acquire(double t) {
/*------------------------------------------------------------------------------
    check will UN-conditionnally update stream value(s) */
    int status;
    status=processor(*this, t);
    if(post_processing!=0) post_processing(*this);
    return(status);
    }

  size_t nvalues() {
    const size_t
      s=grid->Hsize();
    return(s);
    }
  
  int allocate() {
    if(grid!=0) {
      const size_t
        s=grid->Hsize();
      x=new T[s];
      return(s*sizeof(T));
      }
// // // //     else {
// // // //       x=new T[gaussian->nvalues];
// // // //       return(gaussian->nvalues*sizeof(T));
// // // //       }
    }
  
  int allocated() {
    return(x!=0);
    }

  int finalize() {
    int status=0;
    if(post_processing!=0) post_processing(*this);
    return(status);
    }

  void init(filestream_t<T> *filestream, grid_t *g, int (*proc) (SGfield_t<T> &, double) ) {
    this->init();
    this->stream=filestream;
    this->grid=g;
    this->processor=proc;
    this->allocate();
    }
  
  SGfield_t(filestream_t<T> *filestream, grid_t *g, int (*proc) (SGfield_t<T> &, double) ) {
    this->init(filestream,g,proc);
    }
  
};


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  stream for data handling, unstructured
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <class T> class UGfield_t: public  datastream_t <T> {
private :
public :
  T *x, mask;                                 /* field values at nodes             */
  char * units;                               /*                                   */
  double time,next;                           /*                                   */
  mesh_t *mesh;                               /* pointer to related mesh           */
  discretisation_t *descriptor;               /* pointer to related descriptor     */
  list_t *list;                               /* pointer to related list           */
  grid_t *listGrid;                           /* pointer to related grid           */
  int volatile_data;                          /*                                   */
  int (*pre_processing)(void);                /*                                   */
  int (*post_processing)(UGfield_t<T> &);     /*                                   */
  int (*processor) (UGfield_t<T> &, double);  /* processor is linked to stream !!! */
  datastream_t< T > *stream;                  /* stream is linked to processor !!! */
  
  void destroy() {
    deletep(&x);
    }
  
  void init() {
    x=0;
    units=0;
    time=-INFINITY;
    next=-INFINITY;
    mesh=0;
    descriptor=0;
    volatile_data=0;
    pre_processing=0;
    post_processing=0;
    processor=0;
    stream=0;
    }
  
  UGfield_t() {
    this->init();
    }
  
//   int connect( int (*) (UGfield_t<T>, double) processor) {
//     UGfield_t<T> f;
//     return(0);
//     }

  int check(double t) {
/*------------------------------------------------------------------------------
    check will conditionnally update stream value(s) */
    int status __attribute__((unused))=0;
    if(processor==0) return(-1);
/*------------------------------------------------------------------------------
    originally meant for ascending time request, do it mor general */
//     if(t > time) status=processor(*this, t);
    if(t < time or t > next)
      status=processor(*this, t);
    return status;
    }

  int acquire(double t) {
/*------------------------------------------------------------------------------
    acquire will UN-conditionnally update stream value(s) */
    int status;
    if(processor==0) return(-1);
    status=processor(*this, t);
    if(post_processing!=0) post_processing(*this);
    return(status);
    }

  int nvalues() {
    return(descriptor->nnodes);
    }

  int allocate() {
    x=new T[descriptor->nnodes];
    return(descriptor->nnodes*sizeof(T));
    }

  int allocated() {
    return(x!=0);
    }

  int finalize() {
    int status=0;
    if(post_processing!=0) post_processing(*this);
    return(status);
    }
  
  void init(filestream_t<T> *filestream, mesh_t *m, discretisation_t *d, int (*proc) (UGfield_t<T> &, double)) {
    this->init();
    this->stream=filestream;
    this->mesh=m;
    if(d->nnodes<=0) TRAP_ERR_EXIT(ENOEXEC,"descriptor not initialised. You MUST call discretisation_init()!\n");
    this->descriptor=d;
    this->processor=proc;
    this->allocate();
    }
  
  UGfield_t(mesh_t *m, discretisation_t *d) {
    this->init(0,m,d,0);
    }
  
  UGfield_t(filestream_t<T> *filestream, mesh_t *m, discretisation_t *d, int (*proc) (UGfield_t<T> &, double)) {
    this->init(filestream,m,d,proc);
    }
  
  UGfield_t(const char *p, const char *c, const char *v, discretisation_t *d, int (*proc) (UGfield_t<T> &, double),
                                                              int (*identify)(const char *, const char *, int *, int)) {
    filestream_t<T> *filestream=new filestream_t<T>(p,c,v,identify);
    this->init(filestream,0,d,proc);
    }
  
  size_t address(int target) const {
    size_t s=(size_t) x;
    return(s);
    }

};


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  stream data, pointwise time serie
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <class T> class LOCALstream_t {
private :
public :
  T *x, mask;                                     /* set of timed values              */
  char * units;                                   /*                                  */
  double time,next;                               /*                                  */
  int nseries;                                    /* number of values in set          */
  int volatile_data;                              /*                                  */
  int (*pre_processing)(void);                    /*                                  */
  int (*post_processing)(LOCALstream_t<T> &);     /*                                  */
  int (*processor) (LOCALstream_t<T> &, double);  /*                                  */
  datastream_t< T > *stream;
  
  void destroy() {
    deletep(&x);
    }
  
  LOCALstream_t() {
    x=0;
    units=0;
    nseries=0;
    volatile_data=0;
    pre_processing=0;
    post_processing=0;
    processor=0;
    stream=0;
    time=-1.e+10;
    next=-1.e+10;
    }
  
//   int connect( int (*) (UGfield_t<T>, double) processor) {
//     UGfield_t<T> f;
//     return(0);
//     }

  int check(double t) {
/*------------------------------------------------------------------------------
    check will conditionnally update stream value(s) */
    int status;
    if(processor==0) return(-1);
    if(t > time) status=processor(*this, t);
    return(status);
    }

  int acquire(double t) {
/*------------------------------------------------------------------------------
    check will UN-conditionnally update stream value(s) */
    int status;
    if(processor==0) return(-1);
    status=processor(*this, t);
    if(post_processing!=0) post_processing(*this);
    return(0);
    }

  int nvalues() {
    return(nseries);
    }

  int allocate() {
    x=new T[nseries];
    return(nseries*sizeof(T));
    }

  int allocated() {
    return(x!=0);
    }

  int finalize() {
    int status;
    if(post_processing!=0) post_processing(*this);
    return(status);
    }

};


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  stream data, multi-variate
  
  aimed to deal with data that need conjunction of more than 1 field (typically
  wind interpolation by range and direction)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <class T> class combofield_t: public  datastream_t <T> {
private :
public :
  T  *fields;                                    /* field values at nodes            */
  double time,next;                              /*                                  */
  size_t nfields;                                /*                                  */
  int (*processor) (combofield_t<T> &, double);  /*                                  */
  datastream_t< T > *stream;
  
  void destroy() {
    for(size_t n=0;n<nfields;n++) fields[n].destroy();
    }
  
  combofield_t() {
    fields=0;
    nfields=0;
    processor=0;
    }
  
  combofield_t(size_t n) {
    nfields=n;
    fields=new T[nfields];
    processor=0;
    stream=0;
    }
  
  int check(double t) {
/*------------------------------------------------------------------------------
    check will conditionnally update stream value(s) */
    int status __attribute__((unused));
    for(size_t n=0;n<nfields;n++) status=fields[n].check(t);
/*-----------------------------------------------------------------------------
    it can be that combo update is limited to fields update*/
    if(processor!=0) status=processor(*this, t);
    return(0);
    }

  int acquire(double t) {
/*------------------------------------------------------------------------------
    check will UN-conditionnally update stream value(s) */
    int status=0;
    for(size_t n=0;n<nfields;n++) {
      if(fields[n].processor!=0){
        status=fields[n].acquire(t);
        if(status!=0)
          return(-1);
        }
      this->time=fields[n].time;
      this->next=fields[n].next;
      }
/*-----------------------------------------------------------------------------
    it can be that combo update is limited to fields update*/
    if(processor!=0) status=processor(*this, t);
    return(status);
    }

  int allocate() {
    fields=new T[nfields];
    return(nfields);
    }

  int allocate(size_t n) {
    nfields=n;
    fields=new T[nfields];
    return(nfields);
    }

  int allocated() {
    int count=0;
    if(fields==0) return(0);
    for(int n=0;n<nfields;n++) {
      if(fields[n].x!=0) count++;
      }
    return(count==nfields);
    }

  size_t address(int target) {
    size_t s=(size_t) &(fields[target]);
    return(s);
    }
  
  int nvalues() {
    int n=fields[0].nvalues();
    return(n);
    }


};


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  stream data, with temporal sequencing (i.e. sequence of class T)
  
  aimed to deal with time interpolation of data available at different frames
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <class T> class sequence_t: public  datastream_t <T> {
private :
  
public :
  T *frames;
  int nframes, order;
  int id;
  int initialised;
  int (*pre_processing)(void);                    /*                                  */
  int (*post_processing)(LOCALstream_t<T> &);     /*                                  */
  int (*processor) (LOCALstream_t<T> &, double);  /*                                  */
  datastream_t< T > *stream;
 
  sequence_t() {
    id=-1;
    order=-1;
    nframes=0;
    frames=0;
    initialised=0;
    pre_processing=0;
    post_processing=0;
    processor=0;
    }

  int init(double t) {
    int status=0;
    double target=t;
    for(int k=0;k<nframes;k++) {
      status=frames[k].check(target);
      if(status!=0) return status;
      target=frames[k].next;
      }
    return(status);
    }
  
  const int initialized() const {
    int status=this->initialised;
    return(status);
    }
  
  void rotate_frames() {
    int i;
    T tmp;
    
    tmp=frames[0];
    
    for(i=1;i<nframes;i++)
      frames[i-1]=frames[i];
    
    frames[nframes-1]=tmp;
    }

  int check(double t, int verbose=0) {
    int status=0;
//    if(frames[0].allocated()==0) return (0);
    if(frames==0) return (0);
    double next;
  
    switch (order) {
      case LINEAR:
        if(t < frames[0].time or t > frames[1].next) {
/*------------------------------------------------------------------------------
          t out of current sequence, re-init sequence */
          if(verbose>0) printf("t=%lf out of current sequence [%lf:%lf], re-init sequence\n",t,frames[0].time,frames[1].next);
          status=init(t);
          }
        if(status==0 and t > frames[1].time) {
          next=frames[1].next;
          rotate_frames();
          status=frames[1].acquire(next);
          }
        break;
      case QUADRATIC:
        if(t < 0.5*(frames[0].time+frames[1].time) or t > 0.5*(frames[2].time+frames[2].next)) {
/*------------------------------------------------------------------------------
          t out of current sequence, re-init sequence */
          status=init(t);
          }
        if(status==0 and t > 0.5*(frames[1].time+frames[2].time)) {
          next=frames[2].next;
          rotate_frames();
          status=frames[2].acquire(next);
          }
        break;
      default:
        TRAP_ERR_EXIT(-1, "illegal interpolation order\n\n");
        break;
      }
    
    if(status==0) initialised=1;
    
    return(status);
    }
  
  double get_r(double t) {
    /* THIS FUNCTION MUST BE CALLED BEFORE GETTING THE FRAME POINTERS
      BECAUSE check() UPDATES THEM !!! */
    int status;
    double r;
    
    status=check(t);
    if(status!=0) TRAP_ERR_EXIT(-1,"t=%g: datastream check failed",t);
    
    r=(frames[1].time-t)/(frames[1].time-frames[0].time);
    if(
        ( order==LINEAR and (r<0 or 1<r) ) or
        ( order==QUADRATIC and (r<-1 or 1<r) )
        )
      TRAP_ERR_EXIT(-1, "order=%d and r=%g: invalid interpolation\n",order,r);
    
    return r;
    }
  
/*------------------------------------------------------------------------------
 
    Linear interpolation :
    
    z=a*r+b*(1-r)  r E [0;+1]

    a=z0
    b=z1

------------------------------------------------------------------------------*/
  template <typename C> C interpolate_L(double r, C z0, C z1, C mask) {
    C zz;
    const double a=z0, b=z1;
    
    if(z0==mask or z1==mask) return(mask);
    
    zz = a * r + b * (1-r);
    
    return (zz);
    }

/*------------------------------------------------------------------------------
 
    Quadratic interpolation :
    
    z=a*r*r+b*r+c  r E [-1;+1]

      z=a*r*r+b*r+c  ->  z0=a-b+c    z1=c    z2=a+b+c

      c=z1    z2-z0=2b    z2+z0=2a+2c=2a+2z1

    a=(z2+z0)/2-z1    b=(z2-z0)/2    c=z1

------------------------------------------------------------------------------*/
  template <typename C> C interpolate_Q(double r, C z0, C z1, C z2, C mask) {
    double r2;
    C zz;
    double a, b, c;
    
    if(z0==mask or z1==mask or z2==mask) return(mask);
    
    r=-r;
    
    r2 = r * r;
    
    a = 0.5 * (z2 + z0) - z1;
    b = 0.5 * (z2 - z0);
    c = z1;
    zz = a * r2 + b * r + c;
    
    return (zz);
    }

/*------------------------------------------------------------------------------
  
  wrapper of above
  
------------------------------------------------------------------------------*/
  template <typename C> C interpolate_order(double r, C *a, C mask) {
    C z;
    
    switch (order) {
      case LINEAR:
        z=interpolate_L(r, a[0], a[1], mask);
        break;
      case QUADRATIC:
        z=interpolate_Q(r, a[0], a[1], a[2], mask);
        break;
      }
    
    return z;
    }

/*------------------------------------------------------------------------------
  
  wrapper of above
  
------------------------------------------------------------------------------*/
  template <typename C> C interpolate_order_n(double r, C *a, C *b, C *c, int n, C mask) {
    C z;
    
    switch (order) {
      case LINEAR:
        z=interpolate_L(r, a[n], b[n], mask);
        break;
      case QUADRATIC:
        z=interpolate_Q(r, a[n], b[n], c[n], mask);
        break;
      }
    
    return z;
    }

/*------------------------------------------------------------------------------
  
  looped version of above
  
------------------------------------------------------------------------------*/
  template <typename C> void interpolate_order_all(double r, C *a, C *b, C *c, C *z, C mask) {
    int n,nvalues;
    
    nvalues=frames[1].nvalues();
    
    switch (order) {
      case LINEAR:
        for(n=0;n<nvalues;n++)
          z[n]=interpolate_L(r, a[n], b[n], mask);
        break;
      case QUADRATIC:
        for(n=0;n<nvalues;n++)
          z[n]=interpolate_Q(r, a[n], b[n], c[n], mask);
        break;
      }
    
    }

/*------------------------------------------------------------------------------
  
  time interpolation of all fields in sequenced (multivariate) class
  
------------------------------------------------------------------------------*/
  template <typename C> int interpolate(double t, int target, C *z) {
    int status;
    double r;
    
    r=get_r(t);
    
    C mask=frames[0].mask;
    
    C *a=(C *) frames[0].address(target);
    C *b=(C *) frames[1].address(target);
    C *c=0;
    
    if(order>=QUADRATIC)
      c=(C *) frames[2].address(target);
    
    interpolate_order_all(r,a,b,c,z,mask);
    
    return(0);
    }
  
/*------------------------------------------------------------------------------
  
  time interpolation of nth field in sequenced (multivariate) class
  
------------------------------------------------------------------------------*/
  template <typename C> int interpolate(double t, int target, int n, C & z) {
    int status;
    double r;
    
    r=get_r(t);
    
    C mask=frames[0].fields[0].mask;
    
    C *a=(C *) frames[0].address(target);
    C *b=(C *) frames[1].address(target);
    C *c=0;
    
    if(order>=QUADRATIC)
      c=(C *) frames[2].address(target);
    
    z=interpolate_order_n(r,a,b,c,n,mask);
    
    return z;
    }

/*------------------------------------------------------------------------------
  
  time interpolation of nth field in sequenced (monovariate) class
  
------------------------------------------------------------------------------*/
  template <typename C> double interpolate_mask_template(double t, int n, C mask) {
    double r,z;
    
    r=get_r(t);
    
    C *a=frames[0].x;
    C *b=frames[1].x;
    C *c=0;
    
    if(order>=QUADRATIC)
      c=frames[2].x;
    
    z=interpolate_order_n(r,a,b,c,n,frames[0].mask);
    
    return(z);
    }
  
  double interpolate(double t, int n) {
    double z;
    
    z=interpolate_mask_template(t,n,frames[0].mask);
    
    return(z);
    }
  
/*------------------------------------------------------------------------------
  
  time interpolation of all fields in sequenced (monovariate) class
  
------------------------------------------------------------------------------*/
  template <typename C> int interpolate(double t, C *z) {
    double r;
    
    r=get_r(t);
    
    C mask=frames[0].mask;
    
    C *a=frames[0].x;
    C *b=frames[1].x;
    C *c=0;
    
    if(order>=QUADRATIC)
      c=frames[2].x;
    
    interpolate_order_all(r,a,b,c,z,mask);
    
    return(0);
    }
  
  int allocate(size_t n) {
    nframes=n;
    frames=new T[nframes];
    switch (n) {
      case 2:
        order=LINEAR;
        break;
      case 3:
        order=QUADRATIC;
        break;

      }
    return(nframes);
    }
};


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  stream data, generic wrapper?
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template < class C > class fieldstream_t: public  datastream_t <float> {
private :
public :
  C data;
  int id;
  double time() {
    return(data.time);
    }
  int check(double time) {
    return(0);
    }
  double update(double time) {
    return(data.time);
    }
  
  fieldstream_t() {
    id=-1;
    }
  
  fieldstream_t(C data) {
    this->data=data;
    }
};

typedef combofield_t< UGfield_t<float> >      Ucombo_t;
typedef combofield_t< SGfield_t<float> >      Scombo_t;
typedef combofield_t< LOCALstream_t<float> >  Tcombo_t;


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  stream data, atmopsheric stream
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <class T> class atmosphere_t {
private :
public :
  UGfield_t<T> p[3],u10m[3],v10m[3];           /* field values at nodes            */
  UGfield_t<T> uStar[3],vStar[3];              /* field values at nodes            */
  UGfield_t<T> uRange[3],uDirection[3];        /* field values at nodes            */
  UGfield_t<T> wsx[3],wsy[3];                  /* field values at nodes            */
  sequence_t<Ucombo_t> winds;                  /* field values at nodes            */
  int time_interpolation, nframes;
  int initialized;
  
  atmosphere_t() {
    initialized=0;
    time_interpolation=-1;
    }
};


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  stream data, surface waves
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

template <class T> class swave_t {
private :
public :
  UGfield_t<T> p[3],u10m[3],v10m[3];          /* field values at nodes            */
  sequence_t<Ucombo_t> winds;                 /* field values at nodes            */
  int time_interpolation, nframes;
  int initialized;
  
  swave_t() {
    initialized=0;
    time_interpolation=-1;
    }
};


extern int stream_UGxSG(UGfield_t<float> & Ufield, double t);

extern int stream_SGxNETCDF(SGfield_t<float>  & field, double t);
extern int stream_SGxNETCDF(SGfield_t<short>  & field, double t);
extern int stream_UGxNETCDF(UGfield_t<double> & field, double t);
extern int stream_UGxNETCDF(UGfield_t<float> & field, double t);

#define USE_UG_SEQUENCE 1
#if USE_UG_SEQUENCE
extern int filestream2UGfield_processor(UGfield_t<float> & UGfield,double time);
extern sequence_t< UGfield_t<float> >  initialise_UG_sequence(double time,
                                              filestream_t<float> *filestream,
                                              int nframes,mesh_t *mesh, discretisation_t *descriptor,
                                              int (* processor) (UGfield_t<float> & field, double t) );
extern sequence_t< UGfield_t<float> >  initialise_UG_sequence(double time,
                                              const string & path, const string & convention, const string & name,
                                              int nframes,mesh_t *mesh, discretisation_t *descriptor,
                                              int (* processor) (UGfield_t<float> & field, double t),
                                              int (*identify)(const char *, const char *, int *, int) );
#endif

extern  int polar_conversion(Ucombo_t & combo, double t);

#endif /* DATASTREAM_H */
