
/**************************************************************************

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

***************************************************************************/

#ifdef PARALLEL
//#include "mpi.h"

// static int emulate_MPI_Allgather (void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,  MPI_Comm comm )
// {
//   int p,q,n,status;
//   int *in,*out;
//
//   for(p=0;p<gnCPUs;p++) {
//
//     }
//   printf("emulate_MPI_Allgather\n");
//
// }

#endif

#define SEQUENTIAL_OUTPUT  0
#define GLOBALIZED_OUTPUT  1
#define PARTITIONED_OUTPUT 2

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
    deletep2D(&incidence,nitems);
    deletep(&cardinal);
    deletep(&CPUs);
    nitems=-1;
    nCPUs=-1;
    }
};

class exchange_t {
private :
public :
  int **zSolvedTable,**zUniversalSolvedTable, *zSolvedCount;
  int **uSolvedTable,**uUniversalSolvedTable, *uSolvedCount;
  int **zExportTable,**zUniversalExportTable, *zExportCount;
  int **uExportTable,**uUniversalExportTable, *uExportCount;
  int **zImportTable,**zUniversalImportTable, *zImportCount;
  int **uImportTable,**uUniversalImportTable, *uImportCount;
  int *CPUs, nCPUs, CPUid;

  exchange_t() {
    zSolvedTable=0;zUniversalSolvedTable=0;zSolvedCount=0;
    uSolvedTable=0;uUniversalSolvedTable=0;uSolvedCount=0;
    zExportTable=0;zUniversalExportTable=0;zExportCount=0;
    uExportTable=0;uUniversalExportTable=0;uExportCount=0;
    zImportTable=0;zUniversalImportTable=0;zImportCount=0;
    uImportTable=0;uUniversalImportTable=0;uImportCount=0;
    CPUs=0;nCPUs=-1, CPUid=-1;
    }

  exchange_t(int id, int n) {
    zSolvedTable=new int*[n];zUniversalSolvedTable=new int*[n];zSolvedCount=new int[n];
    uSolvedTable=new int*[n];uUniversalSolvedTable=new int*[n];uSolvedCount=new int[n];
    zExportTable=new int*[n];zUniversalExportTable=new int*[n];zExportCount=new int[n];
    uExportTable=new int*[n];uUniversalExportTable=new int*[n];uExportCount=new int[n];
    zImportTable=new int*[n];zUniversalImportTable=new int*[n];zImportCount=new int[n];
    uImportTable=new int*[n];uUniversalImportTable=new int*[n];uImportCount=new int[n];
    CPUs=new int[n];nCPUs=n;
    CPUid=id;
    }

  void destroy() {
    }

};

class partition_t {
private :
public :
  connexion_t *eConnexion;
  mesh_t      *mesh;
  int         *CPUs, nCPUs;
  int maxdepth;

  partition_t() {
    }

  void destroy() {
    }
};

#define TAG_Z_VALUE   100
#define TAG_Z_NODE    101
#define TAG_Z_COUNT   102

#define TAG_U_VALUE   200
#define TAG_U_NODE    201
#define TAG_U_COUNT   202

#define CONTEXT_LOCAL  0
#define CONTEXT_GLOBAL 1

extern connexion_t fe_partition(mesh_t mesh, int, mesh_t *submesh, int target, int npartition);
extern exchange_t communications_init(mesh_t mesh, connexion_t eConnexions, int CPUid, int nCPUs, int verbose);

extern  mesh_t     *fe_partition_init(mesh_t mesh, int maxdepth, int npartition);
extern  partition_t fe_partition_AllCPUs(mesh_t mesh, int maxdepth, mesh_t *, int, int npartition);
extern  exchange_t  emulate_communications_init(mesh_t mesh, connexion_t eConnexions, int CPUid, int nCPUs);
extern  int communications_summary(exchange_t exchange);
extern  int MatrixLocal2Global( hypermatrix_t *M, mesh_t *mesh);
extern  int MatrixLocal2Global( hypermatrix_t *M, distribution_t distributor);
 
extern  connexion_t fe_partition_eConnexions(mesh_t mesh, mesh_t *splitted,int npartition);

extern int export_init_LGP1xLGP1(mesh_t, connexion_t, int, int **, int **, int *, int **, int **, int *);

extern int export_init_LGP1(mesh_t, connexion_t, int, int **, int **, int *);
extern int export_init_GENERIC(mesh_t, connexion_t, int, int, int **, int **, int *);

extern connexion_t partition_zConnexions_LGP1(mesh_t mesh, connexion_t eConnexions);

extern  connexion_t fe_eLocalConnexion(mesh_t mesh, connexion_t general, mesh_t partitioned, int target);

extern int OutputVector(int    *local, int    **global, mesh_t & mesh, int discretisation, int cpu, int target);
extern int OutputVector(float  *local, float  **global, mesh_t & mesh, int discretisation, int cpu, int target);
extern int OutputVector(double *local, double **global, mesh_t & mesh, int discretisation, int cpu, int target);

extern int OutputVector(complex<float>  *local, complex<float>  **global, mesh_t & mesh, int discretisation, int cpu, int target);
extern int OutputVector(complex<double> *local, complex<double> **global, mesh_t & mesh, int discretisation, int cpu, int target);

extern int VectorLocal2Global(discretisation_t descriptor, double *local, double **global, int cpu, int target);
extern int VectorLocal2Global(mesh_t mesh, int discretisation, double *local, double **global, int cpu, int target);
extern int VectorGlobal2local(mesh_t mesh, int discretisation, double *local, double *global);
extern int VectorGlobal2local(mesh_t mesh, int discretisation, int *local, int *global);
