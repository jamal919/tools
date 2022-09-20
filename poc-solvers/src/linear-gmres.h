
/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2017

  Unstructured Ocean Grid initiative

Developers:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
Contributors:
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France
  Frédéric Dupont    Canada environement, Canada

***************************************************************************/

#ifndef LINEAR_GMRES_H
#define LINEAR_GMRES_H

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <complex>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "buffers.h"
#include "solvers-functions.h"
#include "mpi-communications.h"
#include "poc-solvers.h"

#define VPTR_64 1

#if VPTR_64
// typedef unsigned long vptr_t;
typedef size_t vptr_t;
#else
// typedef unsigned int vptr_t;
typedef int vptr_t;
#endif

class ordering_t {
private :
public :
  int *incidence;      /* array[pointer[nrows]] of column indexes of non-0 elements */
  int *cardinal;       /* array[nrows] of count of non-0 elements */
  vptr_t *pointer;     /* array[nrows+1] of indexes of first non-0 element concatenated with the total count of non-0 elements */
  
  int *rows;           /* if set (is it ever ?), array[<nrows] of the indexes of the rows with at least 1 non-0 element */
  int ncols, nrows, max_cardinal;
  int hbw;
  int shared;
  bool duplicated;
  int type,numbering;
  size_t blocksize, rowsize, colsize, offset;   /* for structured 3D extension of 2D ordering (i.e. 3D matrix with 2D ordering) */
  vector<size_t> hdim, vdim, start;             /* for multi-variate linear systems */
  
//  ordering_t& operator= (const ordering_t& m);   // Opérateur d'affectation
//  ordering_t& operator= (const ordering_t& m) {
//     ordering_t dummy;
//     dummy=m;
//     dummy.shared=1;
//     return (dummy);
//     }

  void push(ordering_t **o) {
    this->shared++;
    *o=this;
    }

  void pop() {
    this->shared--;
    }

  size_t _blocksize() {
    size_t tmp;
    tmp=rowsize*colsize;
    blocksize=tmp;
    return(tmp);
    }

  size_t size() {
/*------------------------------------------------------------------------------
    returns sparse matrix size (to be multiplied by matrix type's size to get
    size in bytes) */ 
    size_t offset,tmp;
    bool debug=false;
    
    switch(this->numbering) {
      case 0:
        offset=0;
        break;
      case 1:
        offset=1;
        break;
      default:
//         STDOUT_BASE_LINE("illegal numbering\n", this->numbering);
        exit(-1);
        break;
      }
    if(pointer!=0) {
/**----------------------------------------------------------------------------
      when numbering changed to F-numbering, not only incidence (i.e. row or
      column) numbers are incremented by 1, but also pointer as it is used to
      indicates the offset (position in array) in F-numbering
-----------------------------------------------------------------------------*/
      tmp=(pointer[nrows]-offset)*this->_blocksize();
      if(debug) printf("bloksize=%lu (%lux%lu), nblocks=%lu\n", this->_blocksize(), rowsize, colsize, pointer[nrows]-offset);
      }
    else {
      tmp=-1;
      printf("ordering pointer not initialised\n");
      }
    return(tmp);
    }

  int pos(int row, int col) {
/**----------------------------------------------------------------------------
    seek for col position in row incidence array, standard packed matrix */
    size_t k, cpos;
    int shift=0,scol;
    bool fortran=(numbering==1);
    if(fortran) {
      shift=1;
      }
    switch(this->type) {
      case COO:
        break;
      case CSR:
      case CSC:
//         int *tab=&(incidence[pointer[row]-shift]);
//         for(k=0;k<cardinal[row];k++) {
//           if(tab[k]==col+shift) {
//             return(pointer[row]-shift+k);
//             }
//           }
//         break;
        cpos=pointer[row]-shift;
        scol=col+shift;
        int *tab=&(incidence[cpos]);
        for(k=0;k<cardinal[row];k++) {
          if(tab[k]==scol) {
            return(cpos+k);
            }
          }
        break;
      }
    printf("%s : could not find pos for row=%d col=%d\n",__FUNCTION__,row,col);
    return(-1);
    }

  int *Istart(int row) {
/**----------------------------------------------------------------------------
    seek for col position in row incidence array, standard packed matrix */
    int shift=0;
    if(numbering==1) {
      shift=1;
      }
    switch(this->type) {
      case COO:
        break;
      case CSR:
      case CSC:
        int *tab=&(incidence[pointer[row]-shift]);
        return(tab);
        break;
      }
    printf("%s : could not find start for row=%d\n",__FUNCTION__,row);
    return(0);
    }

  int pos(int row, int col, int k, int l) {
/**----------------------------------------------------------------------------
    seek for col position in row incidence array, block-wise packed matrix*/
    size_t i;
    int *tab=&(incidence[pointer[row]]);
    for(i=0;i<cardinal[row];i++) {
      if(tab[i]==col) {
        return((pointer[row]+i)*blocksize+k*colsize+l);
        }
      }
    printf("%s : could not find pos for row=%d col=%d\n",__FUNCTION__,row,col);
    return(-1);
    }

  int col(int row, int k) {
/**----------------------------------------------------------------------------
    return k-th col index in row incidence array, standard packed matrix*/
    
    if(k<cardinal[row]) {
      return(incidence[pointer[row]+k]);
      }
    else {
      return(-1);
      }
    }

  ordering_t() {
    reset();
    }

  void reset() {
    incidence=NULL;
    cardinal=NULL;
    pointer=NULL;
    rows=0;
    ncols=-1;
    nrows=-1;
    max_cardinal=-1;
    duplicated=0;
    shared=0;
    hbw=-1;
    type=CSR;
    numbering=0;
    blocksize=rowsize=colsize=1;
    offset=0;
    hdim.clear(); vdim.clear(); start.clear();
    }

  void allocate(int nrows0) {
/*------------------------------------------------------------------------------
    when nrows is known, allocation of cardinal and pointer can be done  */
    nrows    = nrows0;
    cardinal = new int[nrows];
    pointer  = new vptr_t[nrows + 1];
    for(size_t n=0;n<nrows;n++) cardinal[n]=0;
    for(size_t n=0;n<nrows;n++) pointer[n]=-1;
    pointer[nrows]=0;
    }

  void allocate_finalize() {
/*------------------------------------------------------------------------------
    when cardinal is known, allocation/initialisation can be finalized  */ 
    pointer[0]=0;
    for(size_t n=0;n<nrows;n++) {
      pointer[n+1]=pointer[n]+cardinal[n];
      }
    size_t nnz=pointer[nrows];
    incidence =new int[nnz];
    }

  void duplicate(const ordering_t & source) {
    (*this)=source;
    int size=pointer[nrows];
    incidence= new int[size];
    cardinal = new int[nrows];
    pointer  = new vptr_t[nrows+1];
    memcpy(incidence,source.incidence,size*sizeof(int));
    memcpy(cardinal,source.cardinal,nrows*sizeof(int));
    memcpy(pointer,source.pointer,(nrows+1)*sizeof(vptr_t));
    }

  void destroy() {
    if(incidence!=NULL) delete[] incidence;
    if(cardinal!=NULL)  delete[] cardinal;
    if(pointer!=NULL)   delete[] pointer;
    if(rows!=NULL)      delete[] rows;
//     hdim.clear();
//     vdim.clear();
    (*this).reset();
    }


  void gethbw() {
    int k,m,row;
    int delta,hbw;
    hbw=-1;
    for(row = 0; row < this->nrows; row++) {
      for(k = 0; k < this->cardinal[row]; k++) {
        m = this->incidence[this->pointer[row]+k];
        delta=abs(m-row);
        hbw=max(delta,hbw);
        }
      }
    this->hbw=hbw;
    }
    
  range_t<int> stat() {
    range_t<int> r;
    r=minmax<int>(this->cardinal,this->nrows);
    return(r);
    }
};


template <typename T> class HYPERMATRIX_t{
private :
public :
  ordering_t *ordering;
  int ordering_owner;
  T *packed;
  T *normalisation;
  SOLVER_t<T> *s;
  distribution_t *distributor;
  exchange_t     *mpi_exchange;
  int type, neq, hbw;
  string name;
  buffers_t<double> mpi_buffers[2];
  void *periodic;
  int  discretisation;
  bool IsDiagonal;

  HYPERMATRIX_t() {
    ordering =0;
    packed   =0;
    normalisation=0;
    s=        0;
    distributor  =0;
    mpi_exchange =0;
    type     =0;
    neq      =-1;
    hbw      =-1;
    ordering_owner=1;
    name="";
    discretisation=-1;
    IsDiagonal=0;
    }

  int allocate(void) {
    int status=0;
    ptr_delete(this->packed);
    size_t N=this->ordering->size();
    this->packed=new double[N];
    if(this->packed==0) status=-1;
    return(status);
    }

  void assign(double value) {
    int n,N;
    if(this->ordering!=0) {
      N=this->ordering->size();
      for(n=0;n<N;n++) this->packed[n]=value;
      }
    }

  int allocate(double value) {
    int status=0;
    status=this->allocate();
    if(status!=0) return(-1);
    this->assign(value);
    return(0);
    }
    
  size_t nzero() {
    size_t count;
    double zero=0.0;
    count=occurence(zero, packed, ordering->size());
    return(count);
    }

  void destroy() {
    bool debug=false;
    int verbose=0;
    name.clear();
    
/**----------------------------------------------------------------------------
    might have been unintialized or freed just before */
    if(this->s!=0){
      int status __attribute__((unused));
      status=solver_terminate(this->s, verbose, debug);
/**----------------------------------------------------------------------------
      might be dangerous */
      free(this->s);/* You MUST use free because the poc-solvers are coded in C. */
      this->s=0;
      }
/**----------------------------------------------------------------------------
    KILLER: adresses could have been passed to solver, and potentially free'd
    This the case when using MUMPS solver
    to be fixed
-----------------------------------------------------------------------------*/
//     if((t.nnz!=0) && (t.x==0)) {
//       if(this->ordering!=0) ordering->incidence=0;
//       if(this->ordering!=0) ordering->pointer=0;
//       packed=0;
//       }
    if(this->ordering!=0) {
      if(this->ordering_owner==1) {
        (*this->ordering).destroy();
	delete this->ordering;
        }
      this->ordering=0;
      }
    if(packed!=0) {
      delete[] packed;
        ///  Probleme avec ce delete à la toute fin du programe ???
      packed   =0;
      }
    if(normalisation!=0) {
      delete[] normalisation;
      normalisation=0;
      }
    neq=-1;
    hbw=-1;
    distributor  =0;
    mpi_exchange =0;
    mpi_buffers[0].destroy();
    mpi_buffers[1].destroy();
    }

};

typedef  HYPERMATRIX_t<double> hypermatrix_t;
typedef  HYPERMATRIX_t< complex<double> > hyperzmatrix_t;

#define index2D3D(pointer,cardinal,K,L,nU,pos,k,l)  ((pointer[nU]*K+k*cardinal[nU]+pos)*L+l)

extern int vector_checkInf(double *A, int zdim);
extern int vector_checkInf(complex<double> *A, int zdim);

extern int vector_checkNaN(float  *A, int zdim);
extern int vector_checkNaN(double *A, int zdim);
extern int vector_checkNaN(complex<double> *A, int zdim);

extern int matrix_operation(hypermatrix_t&,  double *,           double *, int );
extern int matrix_operation(hyperzmatrix_t&, complex<double>  *, complex<double>  *, int );
extern int matrix_operation(hypermatrix_t&,  complex<double>  *, complex<double>  *, int );


template<typename T0> class vector_template_t {
private:
  bool aliased;
public:
  T0 *x;
  size_t size;

  void init(size_t n){
    size=n;
    x=new T0[size];
    aliased=false;
    }
      
  void init(size_t n, T0 c){
    this->init(n);
    for(size_t m=0; m<size;m++) x[m]=c;
    }
    
  void alias(T0 *p,size_t n){
    size=n;
    x=p;
    aliased=true;
    }
    
  vector_template_t(){
    x=0;
    size=0;
    aliased=false;
    }
      
  vector_template_t(size_t n){
    this->init(n);
    }
      
  void duplicate(T0 *p,size_t n){
    this->init(n);
    for(size_t m=0; m<size;m++) x[m]=p[m];
    }
    
  vector_template_t (T0 *p,size_t n) {
    this->alias(p, n);
    }
    
  template<typename T> inline vector_template_t& operator *= (const T c) {
    for(size_t m=0; m<size;m++) {
      x[m] *=  c;
      }
    return *this;
    }
  
  template<typename T> inline vector_template_t& operator /= (const T c) {
    for(size_t m=0; m<size;m++) {
      x[m] /=  c;
      }
    return *this;
    }
  
  inline vector_template_t & operator-=(const vector_template_t & a){
    for(size_t m=0; m<size;m++) {
      x[m] -= a.x[m];
      }
    return *this;
    }
  
  inline vector_template_t & operator+=(const vector_template_t & a){
    for(size_t m=0; m<size;m++) {
      x[m] += a.x[m];
      }
    return *this;
    }
  
  void clear(){
    if(this->x!=0 and not this->aliased) delete[] this->x;
    x=0;
    size=0;
    aliased=false;
    }
      
  ~vector_template_t(){
    if(this->x!=0 and not this->aliased) delete[] this->x;
    }
      
  };

// class dvector_t : public vector_template_t<double> {
// public:  
// };
  
// class zvector_t : public vector_template_t< complex<double> >{
// public:  
// };

typedef vector_template_t< double >          dvector_t;
typedef vector_template_t< complex<double> > zvector_t;
 
  inline double operator*(const dvector_t & p, const dvector_t & q) {
    double s=0;
    for(size_t m=0; m<p.size;m++) {
      s+=p.x[m]*q.x[m];
      }
    return s;
    }
  
  inline complex<double> operator*(const zvector_t & p, const zvector_t & q) {
    complex<double> s=0;
    for(size_t m=0; m<p.size;m++) {
      s+=p.x[m]*conj(q.x[m]);
      }
    return s;
    }
  
  inline double norm(const dvector_t & p) {
    double s=sqrt(p*p);
    return s;
    }
    
  inline double norm(const zvector_t & p) {
    complex<double> z=0.0;
    double s;
    for(size_t m=0; m<p.size;m++) {
      z+=p.x[m]*conj(p.x[m]);
      }
    s=sqrt(abs(z));
    return s;
    }
  
template <typename T> class gmres_t {
private:
public:
  int nitermax, nunknowns, niter;
  vector_template_t<T> r;         // residual vector
  vector_template_t<T> cs,sn;     // Given rotation coeffcients
  vector_template_t<T> g;         // pseudo-residual vector
  dvector_t e1;                   // residual norm vector
  T *y;
  T *H,*Q;
  bool initialised;
  
  gmres_t() {
    H=Q=0;
    initialised=false;
    }
    
  void init(int m, int n) {
    nitermax =m;
    nunknowns=n;
    r.init(n,0.0);
    sn.init(m,0.0);
    cs.init(m,0.0);
    g.init(m,0.0);
    e1.init(n,0.0);
    y=new T[m+1];
    Q=new T[n*m];
//     for(size_t k=0; k<m*n; k++) Q[k]=0.0;
    H=new T[(m+1)*m];
//     for(int k=0; k<(m+1)*m; k++) H[k]=0.0;
    initialised=true;
    }
    
  void clear(){
    sn.clear();
    cs.clear();
    g.clear();
    e1.clear();
    ptr_delete(y);
    ptr_delete(Q);
    ptr_delete(H);
    initialised=false;
    }

};

//   extern int gmres_solver(mesh_t & mesh, hypermatrix_t & A, hypermatrix_t & P, double *rhs, double *x, double *e, int max_iterations, double threshold, gmres_t & gmres);
  extern int gmres_solver(hypermatrix_t & A, hypermatrix_t & P, double *rhs, double *x, double *e, int max_iterations, double threshold, gmres_t<double> & gmres);

#endif
