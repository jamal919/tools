
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief sparse matrix operation declaration
*/
/*----------------------------------------------------------------------------*/

#if MATRIX_H == 0
#define MATRIX_H 1

#include <complex>

#include <poc-solvers.h> 

#include "poc-array.hpp" /* for deletep() */
#include "poc-array-stats.hpp" /* for class range_t */

#define fcomplex complex<float>
typedef complex<double> dcomplex;

#define ORDER_UNDEF -1
#define ORDER_COO   COO
#define ORDER_CSC   CSC
#define ORDER_CSR   CSR
#define ORDER_BND   BND

#define ORDER_UNDEF_NAME "UNDEFINED"
#define ORDER_COO_NAME   COO_NAME
#define ORDER_CSC_NAME   CSC_NAME
#define ORDER_CSR_NAME   CSR_NAME
#define ORDER_BND_NAME   BND_NAME

#if DEV_SOLVER_ACTIVATE
#define VPTR_64 1
#else
#define VPTR_64 0
#endif


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
  int *incidence;  /* array[pointer[nrows]] of column indexes of non-0 elements */
  int *cardinal;   /* array[nrows] of count of non-0 elements */
  int *pointer;    /* array[nrows+1] of indexes of first non-0 element concatenated with the total count of non-0 elements */
  int nrows, ncols, max_cardinal;
  int hbw;
  int shared;
  int type,numbering;
  size_t blocksize, rowsize, colsize;   /* introduced for structured 3D extension of 2D ordering */
  
  int *rows;       /* if set (is it ever ?), array[<nrows] of the indexes of the rows with at least 1 non-0 element */
  bool duplicated;
  size_t offset;   /* for structured 3D extension of 2D ordering (i.e. 3D matrix with 2D ordering) */
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

  size_t size(int verbose=0) {
    size_t offset,tmp;
    switch(this->numbering) {
      case 0:
        offset=0;
        break;
      case 1:
        offset=1;
        break;
      default:
        STDOUT_BASE_LINE("illegal numbering\n", this->numbering);
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
      if(verbose==1) printf("bloksize=%lu (%lux%lu), nblocks=%lu\n", this->_blocksize(), rowsize, colsize, pointer[nrows]-offset);
      }
    else {
      tmp=-1;
      printf("ordering pointer not initialised\n");
      }
    return(tmp);
    }

  int pos(int row, int col) {
    size_t k;
    int *tab=&(incidence[pointer[row]]);
    for(k=0;k<cardinal[row];k++) {
      if(tab[k]==col) {
        return(pointer[row]+k);
        }
      }
    printf("%s : could not find pos for row=%d col=%d\n",__FUNCTION__,row,col);
    return(-1);
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

  ordering_t() {
    reset();
    }

  void reset() {
    incidence=NULL;
    cardinal=NULL;
    pointer=NULL;
    nrows=-1;
    max_cardinal=-1;
    shared=0;
    hbw=-1;
    type=ORDER_CSR;
    numbering=0;
    blocksize=rowsize=colsize=1;
    }

  void destroy() {
    if(incidence!=NULL) delete[] incidence;
    if(cardinal!=NULL)  delete[] cardinal;
    if(pointer!=NULL)   delete[] pointer;
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

class distribution_t{
private :
public :
  size_t nsolved,ndof,ngdof;        /* locally solved nodes, number of nodes in partitioned mesh, total number of nodes */
  size_t *solved,*gsolved,*rsolved; /* list of solved nodes (local index), list of solved nodes (global index)          */
  size_t *gindex,*rindex;           /* index of local nodes in global numbering                                         */
  size_t gmin,gmax,gsize;
  int context;
  int reorder;

  distribution_t() {
    nsolved =  0;
    ndof    = -1;
    ngdof   = -1;
    solved  =  0;
    gsolved =  0;
    rsolved =  0;
    gindex  =  0;
    rindex  =  0;
    context = -1;
    reorder =  0;
    gmin    = -1;
    gmax    = -1;
    gsize   = -1;
    }
  
  void destroy() {
    nsolved =  0;
    ndof    = -1;
    ngdof   = -1;
    deletep(&solved);
    deletep(&gsolved);
    deletep(&rsolved);
    deletep(&gindex);
    deletep(&rindex);
    context = -1;
    reorder =  0;
    gmin    = -1;
    gmax    = -1;
    gsize   = -1;
    }
};

class hypermatrix_t{
private :
public :
  ordering_t *ordering;
  int shared_ordering;
  double *packed;
  double *band;
  int    *pivot;
  double *info, *control;
  triplet t;
  solver_t *s;
  distribution_t *distributor;
  int type, hbw;
  string name;
  double *mpi_global,*mpi_tmp;
  int discretisation;
  bool IsDiagonal;
  const vector<paire_t> * periodic;  // to be reworked

  hypermatrix_t() {
    ordering =0;
    packed   =0;
    band     =0;
    pivot    =0;
    info     =0;
    control  =0;
    s=        0;
    distributor=        0;
    type     =0;
    hbw      =-1;
    shared_ordering=1;
    name=(string) "";
    t.nrow=t.ncol=t.nzmax=t.nnz=0;
    t.mu=t.ml=t.lda=-1;
    t.stype=0;
    t.comptype=-1;
    t.base=0;
    t.i=t.j=0;
    t.x=0;
    t.a=0;
    t.pivot=0;
    t.A=0;
    t.S=0;
    t.N=0;
    mpi_global  = mpi_tmp  =0; /*buffer for parallel computing */
    discretisation=-1;
    IsDiagonal=0;
    }

  int allocate(int verbose=0) {
    deletep(&this->packed);
    size_t N=this->ordering->size(verbose);
    this->packed=new double[N];
    return(0);
    }

  void assign(double value, int verbose=0) {
    int n,N;
    if(this->ordering!=0) {
      N=this->ordering->size(verbose);
      for(n=0;n<N;n++) this->packed[n]=value;
      }
    }

  int allocate(double value, int verbose=0) {
    this->allocate(verbose);
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
    
    int status;
    
/**----------------------------------------------------------------------------
    might have been unintialized or freed just before */
    if(this->s!=0) status=solver_terminate(this->s);
    this->s=0;
/**----------------------------------------------------------------------------
    KILLER: adresses could have been passed to solver, and potentially freed
    This the case when using MUMPS solver
    to be fixed
-----------------------------------------------------------------------------*/
    if((t.nnz!=0) && (t.x==0)) {
      if(this->ordering!=0) ordering->incidence=0;
      if(this->ordering!=0) ordering->pointer=0;
      packed=0;
      }
    if(this->ordering!=0) {
      if((*this->ordering).shared==0) {
        (*this->ordering).destroy();
        }
      delete this->ordering;
      this->ordering=0;
      }
    if(packed!=0) {
      delete[] packed;
        ///  Probleme avec ce delete à la toute fin du programe ???
      packed   =0;
      }
    if(band!=0) {
      delete[] band;
      band=0;
      }
    if(pivot!=0) {
      delete[] pivot;
      pivot=0;
      }
    hbw=-1;
    deletep(&mpi_global);
    deletep(&mpi_tmp);
/**----------------------------------------------------------------------------
    to be fixed */
//     if(d!=0) {
//       delete d;
//       d=0;
//       }
    }

  int neq(){
    if(ordering==NULL)
      return -1;
    return ordering->nrows;
    }

  int *pneq(){
    if(ordering==NULL)
      return NULL;
    return &ordering->nrows;
    }

};

class hyperzmatrix_t{
private :
public :
  ordering_t *ordering;
  int shared_ordering;
  complex <double> *packed;
  double *info, *control;
  tripletz t;
  solver_t   *s;
  int type, hbw;
  complex <double>  *band;
  int    *pivot;
  distribution_t *distributor;
  string name;
  complex<double> *mpi_global,*mpi_tmp;
  int discretisation;
  bool IsDiagonal;
  const vector<paire_t> * periodic;  // to be reworked

  hyperzmatrix_t() {
    ordering =0;
    packed   =0;
    info     =0;
    control  =0;
    s=        0;
    distributor=        0;
    type     =0;
    hbw      =-1;
    band     =0;
    pivot    =0;
    shared_ordering=1;
    t.nrow=t.ncol=t.nzmax=t.nnz=0;
    t.mu=t.ml=t.lda=-1;
    t.stype=0;
    t.comptype=-1;
    t.base=0;
    t.i=t.j=0;
    t.x=0;
    t.a=0;
    t.pivot=0;
    t.A=0;
    t.S=0;
    t.N=0;
    mpi_global  = mpi_tmp  =0; /*buffer for parallel computing */
    discretisation=-1;
    IsDiagonal=0;
    }

  int allocate(void) {
    deletep(&this->packed);
    size_t N=this->ordering->size();
    this->packed=new complex<double>[N];
    return(0);
    }

  void assign(double value) {
    int n,N;
    if(this->ordering!=0) {
      N=this->ordering->size();
      for(n=0;n<N;n++) this->packed[n]=value;
      }
    }

  int allocate(double value) {
    this->allocate();
    this->assign(value);
    return(0);
    }

  size_t nzero() {
    size_t count;
    complex<double> zero=complex<double>(0.0, 0.0);
    count=occurence(zero, packed, ordering->size());
    return(count);
    }

  void destroy() {
    if(this->ordering!=0) {
      if((*this->ordering).shared==0) {
        (*this->ordering).destroy();
        }
      }
    if(packed!=0) {
      delete[] packed;
      packed=0;
      }
    if(band!=0) {
      delete[] band;
      band=0;
      }
    if(pivot!=0) {
      delete[] pivot;
      pivot=0;
      }
    hbw=-1;
    deletep(&mpi_global);
    deletep(&mpi_tmp);
    }

  int neq(){
    if(ordering==NULL)
      return -1;
    return ordering->nrows;
    }

  int *pneq(){
    if(ordering==NULL)
      return NULL;
    return &ordering->nrows;
    }

};

class incidence_t {
private :
public :
  vector<int> incidence;
  
  incidence_t() {
    }
};

class expand_t {
private :
public :
  vector<incidence_t> rows;
  
  void init() {
    }
    
  expand_t() {
    }
    
  int maxsize() {
    vector<int> tmp;
    for(int k=0; k<rows.size(); k++) {
      tmp.push_back(rows[k].incidence.size());
      }
    range_t<int> r(tmp);
    return(r.max);
    }
};

extern const double diagonal;

extern int sgmatrix_symmetry(float *, int);
extern int dgmatrix_symmetry(double *, int);

extern int gmatrix_ubw(float *, int, int);
extern int gmatrix_lbw(float *, int, int);

extern int dgmatrix_ubw(double *A, int, int);
extern int dgmatrix_lbw(double *A, int, int);

extern int dgmatrix_transpose(double *, double *, int, int);
extern int gmatrix_transpose(float *, float *, int, int);

extern float *gmatrix_product(float *, float *, int, int, int);
extern double *dmatrix_product(double *, double *, int, int, int);
extern double *dmatrix_product02(double *, double *, int, int, int);

extern float *bmatrix_product(float *, float *, int, int, int, int, int);

extern int matrix_transpose(hyperzmatrix_t *A);

extern int matrix_import(hypermatrix_t A, hypermatrix_t B, hypermatrix_t C);
extern int matrix_import(hyperzmatrix_t A, hyperzmatrix_t B, hyperzmatrix_t C);

extern int matrix_product(hypermatrix_t A, hypermatrix_t B, hypermatrix_t C, int init=0);
extern int matrix_product(hypermatrix_t A, hyperzmatrix_t B, hyperzmatrix_t C, int init=0);
extern int matrix_product(hyperzmatrix_t A, hypermatrix_t B, hyperzmatrix_t C, int init=0);
extern int matrix_product(hyperzmatrix_t A, hyperzmatrix_t B, hyperzmatrix_t C, int init=0);
extern int matrix_product(hypermatrix_t A, hypermatrix_t B, hypermatrix_t *C);
extern int matrix_product(hyperzmatrix_t A, hyperzmatrix_t B, hyperzmatrix_t *C);

extern int matrix_axpy(hypermatrix_t  &, hypermatrix_t &, double );
extern int matrix_axpy(hyperzmatrix_t &, hypermatrix_t &, double );
extern int matrix_axpy(hyperzmatrix_t &, hypermatrix_t &, complex<double>);

extern int matrix_operation(hypermatrix_t &A,const double* x, double* y, int init = 1 );
extern int matrix_operation(hypermatrix_t &A,const complex<double> *x, complex<double> *y, int init=1);
extern int matrix_operation(hyperzmatrix_t &A,const double* x, complex<double>* y, int init = 1 );
extern int matrix_operation(hyperzmatrix_t &A,const complex<double>* x, complex<double>* y, int init = 1 );

extern int vector_checkNaN(double *A, int zdim);
extern int vector_checkNaN(complex<double> *A, int zdim);

extern ordering_t *copy(const ordering_t *s);
extern hyperzmatrix_t copy(hyperzmatrix_t s);

extern int sparsePPM_(hypermatrix_t A,char *path,int d=10);
extern int sparsePPM_(hyperzmatrix_t A,char *path,int d=10);

extern int ordering_expand3D(ordering_t & ordering2D, ordering_t *ordering);
extern int ordering_expand3D(ordering_t & ordering2D, ordering_t *ordering, expand_t & connexions);
extern int ordering_expand3D_01(ordering_t & ordering2D, ordering_t *ordering);
extern int ordering_expand3D_02(ordering_t & ordering2D, ordering_t *ordering);
extern int matrix_Hconcat_01(const hyperzmatrix_t & M1, const hyperzmatrix_t & M2, hyperzmatrix_t & M3);
extern int matrix_Hconcat_02(hyperzmatrix_t & M1, hyperzmatrix_t & M2, hyperzmatrix_t & M3);
extern int matrix_Hconcat(hyperzmatrix_t & M1, hyperzmatrix_t & M2, hyperzmatrix_t & M3);
extern int ordering_Vconcat(ordering_t *M1, ordering_t *M2, ordering_t *M3);

/*----------------------------------------------------------------------------*/
/// shows the image of a matrix
/**
by saving a PPM file.
Calls safely sparsePPM_().
\date created 9 Jun 2011
\author Damien Allain

\param A the matrix
\param path the path of the PPM file
\param d the incrementation. Default: 10
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define sparsePPM(x,...) sparsePPM_(x,__FILE__ ":" TOSTRING(__LINE__) ":" #x ".ppm", ##__VA_ARGS__)
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

extern int free_ordering(ordering_t **ordering);

extern int matrix_check_consistency(int *incidence,int *pointer,int *cardinal,int nnodes);
extern int matrix_check_consistency(ordering_t ordering);

extern int matrix_check_symmetry(int *incidence,int *pointer,int *cardinal,int nnodes, bool debug=false);
extern int matrix_check_symmetry(ordering_t ordering, bool debug=false);

extern int matrix_IsDiagonal(hypermatrix_t & A);
extern int matrix_IsDiagonal(hyperzmatrix_t & A);

extern int matrix_reorder(ordering_t *ordering);

extern int matrix_absolute_pos(ordering_t *ordering, int row, int target);
extern size_t matrix_absolute_pos2(ordering_t *ordering, int row, int start, int target);

extern int matrix_relative_pos(ordering_t *ordering, int row, int target, int StopOnError=1);
extern int matrix_relative_pos2(ordering_t *ordering, int row, int target, int start, int StopOnError);

extern int matrix_CSR2CSC(hypermatrix_t *);
extern int matrix_CSR2CSC(hypermatrix_t *A, ordering_t *CSCordering);

extern int matrix_CSR2CSC(hyperzmatrix_t *A);
extern int matrix_CSR2CSC(hyperzmatrix_t *A, ordering_t *CSCordering);

extern ordering_t *ordering_CSR2CSC(ordering_t *ordering);

extern int matrix_check_Zero(hypermatrix_t A);

extern int ordering_import(ordering_t *A, ordering_t *B, ordering_t *C);

extern int LinearSystem_factorisation(hypermatrix_t *M);
extern int LinearSystem_factorisation(hyperzmatrix_t *M);

extern int LinearSystem_solve(hypermatrix_t &M, double *B);
extern int LinearSystem_solve(hyperzmatrix_t &M, complex <double> *rhs);

extern int LinearSystem_initialise(hypermatrix_t *,  int *, int,  int, bool force_duplicate, int verbose, bool debug);
extern int LinearSystem_initialise(hyperzmatrix_t *, int *, int,  int id=SOLVER_ID_UMFPACK);

extern int LinearSystem_InitTriplet(hypermatrix_t *M, bool duplicate, int verbose, bool debug);
extern int LinearSystem_InitTriplet(hyperzmatrix_t *M, bool duplicate, int verbose, bool debug);
extern int LinearSystem_identify(const char *solver_name);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename type> int pos(const type & y,const type *x,const int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief return index of given value in x vector

\date 2011-09-16 Damien ALLAIN : reviewed and documented

\param y given value
\param *x vector
\param nvalues number of elements in vector
\returns the index of the first value in x that equals y or -1 if none found

array_index() and pos() are NOT identical : the argument order is different
*/
/*----------------------------------------------------------------------------*/
{
  int  n;

  for(n=0;n<nvalues;n++) {
    if(x[n]==y) {
      return n;
      }
    }

  return -1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename type> int dupos(type *x, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief return index of duplicate value in x vector

\date 2011-09-19 Damien ALLAIN : creation

\param *x vector
\param n number of elements in vector
\returns the index of the first duplicate value in x or -1 if none found
*/
/*----------------------------------------------------------------------------*/
{
  int i,j;

  for(i=0;i<n-1;i++) {
    for(j=i+1;j<n;j++) {
      if(x[i]==x[j]) {
        return i;
        }
      }
    }

  return -1;
}

extern int matrix_print(const char *filename, hypermatrix_t & M);
extern int matrix_print_ordering(const char *filename, ordering_t *ordering);

extern int matrix_shrink(hypermatrix_t  & M, const vector<paire_t> & paires, int verbose, bool debug);
extern int matrix_shrink(hyperzmatrix_t & M, const vector<paire_t> & paires, int verbose, bool debug);

extern int matrix_shrink(hyperzmatrix_t & M, paire_t *paires, int npaires, int verbose, bool debug);

#endif
