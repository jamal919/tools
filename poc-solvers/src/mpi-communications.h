
/**************************************************************************

  POC-SOLVERS interface, 2006-2019

  Part of the SIROCCO national service (INSU, France)

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France

***************************************************************************/

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 MPI communication function and classes
 
 to be strictly synchronized with T-UGOm sources (native development sources)
 
 not used for the moment in poc-solvers

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


#ifndef MPI_COMMUNICATIONS_H
#define MPI_COMMUNICATIONS_H

#include <vector>
#include <complex>

#ifdef HAVE_MPI
#include <mpi.h>
#else
#define ompi_datatype_t size_t
#define MPI_Datatype    size_t
#define MPI_Comm        size_t
#define MPI_CHAR    0
#define MPI_INT     0
#define MPI_FLOAT   0
#define MPI_DOUBLE  0
#define MPI_LONG    0
#define MPI_DOUBLE_COMPLEX  0
#define MPI_UNSIGNED_LONG 0
#endif

#include "solvers-functions.h"


class eblock_t {
private :
public :
  int id;
  vector<int> elements;     /* list of elements in block             */
  vector<int> partitions;   /* list of partitions concerned by block */
  
  eblock_t() {
    id=-1;
    }
  void destroy() {
    elements.clear();
    }
};

class connexion_t {
private :
public :
  int **incidence, *cardinal;
  int nitems;
  int *CPUs, nCPUs;
  int maxdepth;

  connexion_t() {
    incidence=0;
    cardinal=0;
    nitems=-1;
    CPUs=0;
    nCPUs=-1;
    }

  void destroy() {
    size_t n;
    if(incidence!=0) {
      for(n=0;n<nitems;n++) {
        deletep(&(incidence[n]));
        }
      delete[] incidence;
      }
    deletep(&cardinal);
    deletep(&CPUs);
    nitems=-1;
    nCPUs=-1;
    }
};


#define NATIVE_SCOPE   0
#define LIMITED_SCOPE  1

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  
  distributed : index in distributed discretisation
  centralized : index in centralized discretisation
  
  when limiting scope (i.e. exchange vector is limited to a subset of discretisation)
  
  distributed : index in limited distributed discretisation (subset)
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

class communication_t {
private :
public :
  int **distributed, **centralized, *count;
  int *incidence,nincidence;
  int nCPUs, CPUid;
  int scope;
  
  void init(int id, int n) {
    this->destroy();
    distributed=new int*[n];
    centralized=new int*[n];
    count=new int[n];
    nCPUs=n;
    CPUid=id;
    for(int k=0;k<nCPUs;k++) {
      distributed[k]=0;
      centralized[k]=0;
      count[k]=0;
      }
    scope=NATIVE_SCOPE;
    }

  void init() {
    distributed=0;centralized=0;count=0;
    incidence=0;
    nincidence=-1;
    nCPUs=-1, CPUid=-1;
    scope=-1;
    }

  communication_t() {
    init();
    }

  communication_t(int id, int n) {
    init(id,n);
    }

  int undisciplined_size() const {
    int size=0;
    for(int k=0; k<nCPUs;k++) {
//       printf("%s: %d %d %d\n",__func__,k,count[k],nCPUs);
      size+=count[k];
      }
    return(size);
    }

  int duplicate(const communication_t c) {
    this->init(c.CPUid, c.nCPUs);
    for(int k=0;k<nCPUs;k++) {
      memmove(distributed[k], c.distributed[k], c.count[k]*sizeof(int));
      memmove(centralized[k], c.centralized[k], c.count[k]*sizeof(int));
      }
    memmove(count, c.count, c.nCPUs);
    nincidence=c.nincidence;
    if(nincidence!=0) {
      memmove(incidence, c.incidence, c.nCPUs*sizeof(int));
      }
    scope=c.scope;
    return(0);
    }

  size_t set_incidence() {
    size_t size;
    int pos, n;
    bool debug=false;
    vector<int> used;
    for(int k=0; k<nCPUs;k++) {
      for(int l=0; l<count[k];l++) {
        n=distributed[k][l];
        pos=vpos(n,used);
        if(pos==-1) used.push_back(n);
        }
      }
    size=used.size();
    incidence=new int[size];
    for(size_t m=0;m<size;m++) {
      incidence[m]=used[m];
      if (debug) printf("incidence %d %d\n",m,incidence[m]);
      }
    used.clear();
    nincidence=size;
    return(size);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
  replace numerical node index with position in numerical node incidence list
  
  It is used to reduce buffer size in communications when a specific initialisation 
  is needed (i.e. buffers not already existing but to be specially created for
  the exchange). buffer is then set only for incidence nodes.
  
  typically meteo forcing provided by top level toward subcycling levels, as time
  in subcycling not identical to top level

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  size_t reduce_scope() {
    size_t size;
    int pos;
    if(scope!=NATIVE_SCOPE) return(-1);
    size=set_incidence();
    for(int p=0;p<nCPUs;p++) {
      for(int k=0;k<count[p];k++) {
        int n=distributed[p][k];
        pos=vpos(n, incidence, nincidence);
        distributed[p][k]=pos;
        }
      }
    scope=LIMITED_SCOPE;
    return(size);
    }

  int *undisciplined_list() {
    int *buffer=new int[this->undisciplined_size()];
    int n=0;
    for(int k=0; k<nCPUs;k++) {
      for(int l=0; l<count[k];l++) {
        buffer[n]=distributed[k][l];
        n++;
        }
      }
    return(buffer);
    }

  int list(vector<int> & v) {
    v.clear();
    for(int k=0; k<nCPUs;k++) {
      for(int l=0; l<count[k];l++) {
        size_t pos=vpos(distributed[k][l],v);
        if(pos==-1) v.push_back(distributed[k][l]);
        }
      }
    quick_sort(v);
    return(0);
    }

  int check() const {
    int status=0;
    if(count==0) return(-1);
    return(status);
    }

  void destroy() {
    for(int k=0; k<nCPUs;k++) {
      deletep(&distributed[k]);
      deletep(&centralized[k]);
      }
    deletep(&distributed);
    deletep(&centralized);
    deletep(&count);
    deletep(&incidence);
    init();
    }
};

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

class exchange_t {
private :
public :
  communication_t cpu_import, cpu_export;
  int *CPUs, nCPUs, CPUid;
  MPI_Comm communicator;
  int discretisation;
  int scope;
  string label;
  
  void init(int id, int n, MPI_Comm c) {
    cpu_import.init(id, n);
    cpu_export.init(id, n);
    communicator=c;
/*------------------------------------------------------------------------------
    not used and probably memory leak false positive in valgrind */
//     CPUs=new int[n];
    nCPUs=n;
    CPUid=id;
    label="";
    scope=-1;
    }

  int check() const {
    int status;
    status=cpu_export.check();
    if(status==-1) return(status);
    status=cpu_import.check();
    if(status==-1) return(status);
    return(status);
    }

  int duplicate(const exchange_t e) {
    int status;
    status=cpu_export.duplicate(e.cpu_export);
    if(status==-1) return(status);
    status=cpu_import.duplicate(e.cpu_import);
    if(status==-1) return(status);
    deletep(&CPUs);
    CPUs=new int[e.nCPUs];
    memmove(CPUs, e.CPUs, e.nCPUs);
    nCPUs=e.nCPUs;
    CPUid=e.CPUid;
    communicator=e.communicator;
    discretisation=e.discretisation;
    scope=e.scope;
    label=e.label;
    return(status);
    }

  exchange_t() {
    CPUs=0;nCPUs=-1;CPUid=-1;communicator=0;discretisation=-1;scope=-1;label="";
    }

//   exchange_t(int id, int n) {
//     init(id,n);
//     }
  
  exchange_t(int id, int n, MPI_Comm c) {
    init(id,n,c);
    }
  
  void destroy() {
    cpu_import.destroy();
    cpu_export.destroy();
    if(CPUs!=0) printf("%s : CPUS array allocated (%x) \n",__func__,CPUs);
    deletep(&CPUs);
    CPUs=0;nCPUs=-1;CPUid=-1;discretisation=-1;
    }
};

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

class centralizer_t {
private :
public :
  size_t **SolvedTable,**UniversalSolvedTable, *SolvedCount;
  size_t **UniversalGlobalTable, *GlobalCount;
  int *CPUs, nCPUs, CPUid;

  void init() {
    SolvedTable=0; UniversalSolvedTable=0; SolvedCount=0;
    UniversalGlobalTable=0; GlobalCount=0;
    CPUs=0; nCPUs=-1, CPUid=-1;
    }

  centralizer_t() {
    this->init();
    }

  int allocate(int id, int n) {
    SolvedTable=new size_t*[n];
    UniversalSolvedTable=new size_t*[n];
    SolvedCount=new size_t[n];
    UniversalGlobalTable=new size_t*[n];
    GlobalCount=new size_t[n];
    for(int k=0;k<nCPUs;k++) {
      SolvedTable[k]=0;
      UniversalSolvedTable[k]=0;
      UniversalGlobalTable[k]=0;
      }
    CPUs=new int[n];nCPUs=n;
    CPUid=id;
    }

  void destroy() {
    deletep(&CPUs);
    deletep(&SolvedCount);
    deletep(&GlobalCount);
    for(int k=0;k<nCPUs;k++) {
      deletep(&SolvedTable[k]);
      deletep(&UniversalSolvedTable[k]);
      deletep(&UniversalGlobalTable[k]);
      }
    this->init();
    }
};


/*------------------------------------------------------------------------------

  parallel exchange class evolution, aimed for a unique dof

------------------------------------------------------------------------------*/
class MaPHyS_t {
private :
public :
  int CPUid, nCPUs;     /* partition ID (cpu), number of partitions (cpus)     */
  int ndofs, nnodes;    /* number of DoFs in partition, number of shared nodes */
  int *nodes;           /* shared nodes array, global numbering                */
  int npartitions;      /* number of connected partitions                      */
  int *partitions;      /* partitions array                                    */
  int *start;           /* incidence pointers array                            */
  int *incidence;       /* incidence array                                     */
  int *locals, nlocals; /* re-numbering (starting by locals), number of locals */

  void init() {
    nnodes=0;
    nodes=0;
    npartitions=0;
    partitions=0;
    start=0;
    incidence=0;
    locals=0;
    nlocals=0;
    }

  void init(size_t n, size_t p) {
    nnodes=n;
    nodes=new int[n];
    npartitions=p;
    partitions=new int[p];
    start=new int[p+1];
    incidence=0;
    }

  MaPHyS_t() {
    this->init();
    }

  void destroy() {
    deletep(&nodes);
    deletep(&partitions);
    deletep(&start);
    deletep(&incidence);
    deletep(&locals);
    this->init();
    }
};


class distribution_t{
private :
public :
  size_t nsolved,ndof,ngdof;           /* locally solved dofs, number of dofs in partitioned mesh, total number of dofs  */
  size_t *solved,*gsolved;             /* list of solved dofs (local index), list of solved dofs (global index)          */
  size_t *gindex;                      /* index of local dofs in global numbering                                        */
  size_t gmin,gmax,gsize;
  MaPHyS_t MaPHyS;                     /* MaPHyS structure                                                               */
  exchange_t mpi_exchange;             /* for fully distributed MPI exchanges                                            */
  centralizer_t centralizer;           /* for root gathering and scaterring                                              */
  int context;                         /* sequential or parallel                                                         */
  int reorder;
  int discretisation;

  distribution_t() {
    nsolved =  0;
    ndof    = -1;
    ngdof   = -1;
    solved  =  0;
    gsolved =  0;
    gindex  =  0;
    context = -1;
    reorder =  0;
    discretisation =  -1;
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
    deletep(&gindex);
    context = -1;
    reorder =  0;
    gmin    = -1;
    gmax    = -1;
    gsize   = -1;
    }
};

#define CONTEXT_LOCAL  0
#define CONTEXT_GLOBAL 1

extern int exchange_atom(const exchange_t & exchange, int CPUid, int nCPUs, char **to_send, char **to_receive, size_t blocksize);
extern int exchange_atom(const exchange_t & exchange, int CPUid, int nCPUs, complex<double> **to_send, complex<double> **to_receive, size_t blocksize);

extern int exchange_atom(const exchange_t & exchange, int CPUid, int nCPUs, char *to_send,  char *to_receive, size_t blocksize);

extern int exchange_atom(const exchange_t & exchange, int CPUid, int nCPUs, char   *to_send, char   *to_receive, bool debug);
extern int exchange_atom(const exchange_t & exchange, int CPUid, int nCPUs, int    *to_send, int    *to_receive, bool debug);
extern int exchange_atom(const exchange_t & exchange, int CPUid, int nCPUs, float  *to_send, float  *to_receive, bool debug);
extern int exchange_atom(const exchange_t & exchange, int CPUid, int nCPUs, double *to_send, double *to_receive, bool debug);
extern int exchange_atom(const exchange_t & exchange, int CPUid, int nCPUs, complex<double> *to_send, complex<double> *to_receive, bool debug);

extern int communication_dualize(communication_t & source, communication_t & target, const int CPUid, const int nCPUs, MPI_Comm communicator, bool debug);
extern int communication_set_distributed(distribution_t & distributor, communication_t & target, const int CPUid, const int nCPUs, bool debug);
extern int *cpu_table(distribution_t & distributor, const int CPUid, MPI_Comm communicator, bool debug);
extern int communications_summary(const exchange_t & exchange, bool debug);

extern int communications_init(distribution_t distributor, int CPUid, int nCPUs, MPI_Comm communicator, exchange_t & exchange, int verbose, bool debug);

extern int mpi_synchronicity(int CPUid, int nCPUs, const char *, int line, bool disciplined);
extern int mpi_synchronicity(MPI_Comm communicator, int CPUid, int nCPUs, const char *caller, int line, bool disciplined);


extern  int allocated_somewhere(size_t buffer, MPI_Comm communicator);
extern  int check_equal(int value, MPI_Comm communicator);

/*------------------------------------------------------------------------------
  protected MPI calls */
extern int P_MPI_Barrier(MPI_Comm communicator);
extern int P_MPI_BarrierFlush(MPI_Comm communicator, FILE *stream);

#endif
