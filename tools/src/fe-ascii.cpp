
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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex>
#include <cmath>

#ifdef HAVE_MPI
//#ifdef PARALLEL
#include "mpi.h"
#endif
using namespace std;

#include "tools-define.h"
#include "tools-structures.h"

#include "fe.h"
#include "archive.h"

#include "parallel.h"

#define NCOMMENTS 2

// static int gCPU_ID, gCPU_MASTER;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int OutputVector(double *local, double **global, mesh_t & mesh, int discretisation, int cpu, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated=0;
  
//  allocated=OutputVector_template(local, global, MPI_DOUBLE, mesh, discretisation, cpu, target);
  
  return(allocated);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int OutputVector(complex<double> *local, complex<double> **global, mesh_t & mesh, int discretisation, int cpu, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated=0;
  
//  allocated=OutputVector_template(local, global, mesh, discretisation, cpu, target);
  
  return(allocated);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int OutputVector(complex<float> *local, complex<float> **global, mesh_t & mesh, int discretisation, int cpu, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated=0;
  
//  allocated=OutputVector_template(local, global, mesh, discretisation, cpu, target);
  
  return(allocated);
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   mesh_t output_mesh(mesh_t local, int gArchiveExtent)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   mesh_t mesh;
//   switch(gArchiveExtent) {
//     case SEQUENTIAL_OUTPUT:
//       return(local);
//       break;
//     case GLOBALIZED_OUTPUT:
//       return(*(local.origin));
//       break;
//     case PARTITIONED_OUTPUT:
//       return(local);
//       break;
//     default:
//       TRAP_ERR_EXIT(-1, "illicit output mode\n");
//     }
//
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int copy_LGP0(mesh_t mesh, float *buffer, float *tmp)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m,mm;

  for (m=0;m<mesh.ntriangles;m++) {
    mm=mesh.triangles[m].origin;
    buffer[m]=tmp[mm];
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int copy_LGP1(mesh_t mesh, float *buffer, float *tmp)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m,mm;

  for (m=0;m<mesh.nvtxs;m++) {
    mm=mesh.vertices[m].origin;
    buffer[m]=tmp[mm];
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int copy_NCP1(mesh_t mesh, float *buffer, float *tmp)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m,mm;

  for (m=0;m<mesh.nedges;m++) {
    mm=mesh.edges[m].origin;
    buffer[m]=tmp[mm];
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int parallel_loadr1(const char *name, mesh_t mesh, float *buffer, float *mask, int fmt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,n,start=0,dum,nitems,nndes;
  char  line[1000];
  float *tmp;
  FILE *in;

  switch(fmt) {
    case LGP0:
      nndes=mesh.origin->ntriangles;
      break;
    case LGP1:
      nndes=mesh.origin->nvtxs;
      break;
    case NCP1:
      nndes=mesh.origin->nedges;
      break;
    }

  tmp=new float[nndes];

  in=fopen(name,"r");
  if(in == NULL) {
    printf("unable to open %s\n",name);
    return (-1);
    }

  n=999;
  for(k=0;k<NCOMMENTS;k++) {
    fgets(line,n,in);
    }

  nitems=sscanf(line,"%d %f",&dum,&(tmp[0]));

  if((nitems ==2) && (dum==1)) start=1;

  for(n=start;n<nndes;n++) {
    nitems=fscanf(in,"%d %f",&dum,&(tmp[n]));
    if(nitems !=2) goto error;
    }

  fclose(in);

  switch(fmt) {
    case LGP0:
      copy_LGP0(mesh, buffer, tmp);
      break;
    case LGP1:
      copy_LGP1(mesh, buffer, tmp);
      break;
    case NCP1:
      copy_NCP1(mesh, buffer, tmp);
      break;
    }

  delete[] tmp;
  return (0);

 error:
  fclose(in);
  delete[] tmp;
  return (-1);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int parallel_loadr2(const char *name, mesh_t mesh, float *bufx, float *bufy, float *mask, int fmt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,n,start=0,dum,nitems,nndes;
  char  line[1000];
  float *tmpx,*tmpy;
  FILE *in;

  switch(fmt) {
    case LGP0:
      nndes=mesh.origin->ntriangles;
      break;
    case LGP1:
      nndes=mesh.origin->nvtxs;
      break;
    case NCP1:
      nndes=mesh.origin->nedges;
      break;
    }

  tmpx=new float[nndes];
  tmpy=new float[nndes];

  in=fopen(name,"r");
  if(in == NULL) {
    printf("unable to open %s\n",name);
    return (-1);
    }

  n=999;
  for(k=0;k<NCOMMENTS;k++) {
    fgets(line,n,in);
    }

  nitems=sscanf(line,"%d %f %f",&dum,&tmpx[0],&tmpy[0]);

  if((nitems ==3) && (dum==1)) start=1;

  for(n=start;n<nndes;n++) {
    nitems=fscanf(in,"%d %f %f",&dum,&tmpx[n],&tmpy[n]);
    if(nitems !=3) goto error;
    }

  fclose(in);

  switch(fmt) {
    case LGP0:
      copy_LGP0(mesh, bufx, tmpx);
      copy_LGP0(mesh, bufy, tmpy);
      break;
    case LGP1:
      copy_LGP1(mesh, bufx, tmpx);
      copy_LGP1(mesh, bufy, tmpy);
      break;
    case NCP1:
      copy_NCP1(mesh, bufx, tmpx);
      copy_NCP1(mesh, bufy, tmpy);
      break;
    }

  delete[] tmpx;
  delete[] tmpy;
  return (0);

 error:
  fclose(in);
  delete[] tmpx;
  delete[] tmpy;
  return (-1);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_ascii_loadr1(const char *filename, mesh_t mesh, float *buffer, float *mask, int fmt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  if(mesh.origin==0) {
    discretisation_t descriptor=get_descriptor(mesh,fmt);
    int nnodes=descriptor.nnodes;
    if (nnodes==0) return(-1);
    status=quoddy_loadr1(filename, nnodes, buffer);
    }
  else {
    status=parallel_loadr1(filename, mesh, buffer, mask, fmt);
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_ascii_loadr2(const char *filename, mesh_t mesh, float *bufx, float *bufy, float *mask, int fmt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  if(mesh.origin==0) {
    status=quoddy_loadr2(filename, mesh.nvtxs, bufx, bufy);
    }
  else {
    status=parallel_loadr2(filename, mesh, bufx, bufy, mask, fmt);
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int fe_ascii_saver1_template(const char *filename, mesh_t & mesh, T *buffer, T mask, int discretisation, int fmt, char **comments, int gArchiveExtent)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated, status;
  mesh_t global;
  T *gbuffer=0;
 
  switch(gArchiveExtent) {
    case SEQUENTIAL_OUTPUT:
      status=quoddy_saver1(filename, mesh.nvtxs, buffer, comments);
      break;
    case GLOBALIZED_OUTPUT:
      global=*(mesh.origin);
//      allocated=OutputVector(buffer, &gbuffer, mesh,discretisation,  gCPU_ID, 0);
      if(gCPU_ID==gCPU_MASTER) status=quoddy_saver1(filename, global.nvtxs, gbuffer, comments);
      if(allocated) delete[] gbuffer;
      break;
    case PARTITIONED_OUTPUT:
      status=quoddy_saver1(filename, mesh.nvtxs, buffer, comments);
      break;
    default:
      TRAP_ERR_EXIT(-1, "illicit output mode\n");
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_ascii_saver1(const char *filename, mesh_t & mesh, double *buffer, double mask, int discretisation, int fmt, char **comments, int gArchiveExtent)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_ascii_saver1_template(filename, mesh, buffer, mask, discretisation, fmt, comments,gArchiveExtent);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_ascii_saver1(const char *filename, mesh_t & mesh, float *buffer, float mask, int discretisation, int fmt, char **comments, int gArchiveExtent)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_ascii_saver1_template(filename, mesh, buffer, mask, discretisation, fmt, comments,gArchiveExtent);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_ascii_saver2(const char *filename, mesh_t & mesh, double *bufx, double *bufy, double mask, int discretisation, int fmt, char **comments, int gArchiveExtent)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated, status;
  mesh_t global;
  double *gbufx=0,*gbufy=0;
  int nnodes;
  discretisation_t descriptor;
  
  descriptor=get_descriptor(mesh,  discretisation);
  
  switch(gArchiveExtent) {
    case SEQUENTIAL_OUTPUT:
      nnodes=descriptor.nnodes;
      status=quoddy_saver2(filename, nnodes, bufx, bufy, comments);
      break;
    case GLOBALIZED_OUTPUT:
      global=*(mesh.origin);
      allocated=OutputVector(bufx, &gbufx, mesh,discretisation,  gCPU_ID, 0);
      allocated=OutputVector(bufy, &gbufy, mesh,discretisation,  gCPU_ID, 0);
      nnodes=descriptor.distributor.ngdof;
      if(gCPU_ID==gCPU_MASTER) status=quoddy_saver2(filename, nnodes, gbufx, gbufy, comments);
      if(allocated) delete[] gbufx;
      if(allocated) delete[] gbufy;
      break;
    case PARTITIONED_OUTPUT:
      nnodes=descriptor.nnodes;
      status=quoddy_saver2(filename, nnodes, bufx, bufy, comments);
      break;
    default:
      TRAP_ERR_EXIT(-1, "illicit output mode\n");
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int fe_ascii_savec1_template(const char *filename, mesh_t & mesh,T *buf, T mask, int discretisation, int fmt, char **comments, int gArchiveExtent)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated, status;
  mesh_t global;
  T *gbuf=0;
  int nnodes;
  discretisation_t descriptor;
  
  descriptor=get_descriptor(mesh,  discretisation);
  
  switch(gArchiveExtent) {
    case SEQUENTIAL_OUTPUT:
      nnodes=descriptor.nnodes;
      status=quoddy_savec1(filename, nnodes, buf, comments);
      break;
    case GLOBALIZED_OUTPUT:
//      global=*(mesh.origin);
      nnodes=descriptor.distributor.ngdof;
      if(nnodes==-1) TRAP_ERR_EXIT(-1, "discretisation distributor not initialized\n");
      allocated=OutputVector(buf, &gbuf, mesh, discretisation,  gCPU_ID, 0);
      if(gCPU_ID==gCPU_MASTER) status=quoddy_savec1(filename, nnodes, gbuf, comments);
      if(allocated) delete[] gbuf;
      break;
    case PARTITIONED_OUTPUT:
      nnodes=descriptor.nnodes;
      status=quoddy_savec1(filename, nnodes, buf, comments);
      break;
    default:
      TRAP_ERR_EXIT(-1, "illicit output mode\n");
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_ascii_savec1(const char *filename, mesh_t & mesh, complex<double> *buf, complex<double>mask, int discretisation, int fmt, char **comments, int gArchiveExtent)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=fe_ascii_savec1_template(filename, mesh, buf, mask, discretisation, fmt, comments, gArchiveExtent);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_ascii_savec1(const char *filename, mesh_t & mesh, complex<float> *buf, complex<float>mask, int discretisation, int fmt, char **comments, int gArchiveExtent)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=fe_ascii_savec1_template(filename, mesh, buf, mask, discretisation, fmt, comments, gArchiveExtent);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_ascii_savec2(const char *filename, mesh_t & mesh, complex<double> *bufx, complex<double> *bufy, complex<double>mask, int discretisation, int fmt, char **comments, int gArchiveExtent)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int allocated, status;
  mesh_t global;
  complex<double> *gbufx=0,*gbufy=0;
  int nnodes;
  discretisation_t descriptor;

  descriptor=get_descriptor(mesh,  discretisation);

  switch(gArchiveExtent) {
    case SEQUENTIAL_OUTPUT:
      nnodes=descriptor.nnodes;
      status=quoddy_savec2(filename, nnodes, bufx, bufy, comments);
      break;
    case GLOBALIZED_OUTPUT:
      global=*(mesh.origin);
      allocated=OutputVector(bufx, &gbufx, mesh,discretisation,  gCPU_ID, 0);
      allocated=OutputVector(bufy, &gbufy, mesh,discretisation,  gCPU_ID, 0);
      if(gCPU_ID==gCPU_MASTER) status=quoddy_savec2(filename, global.nvtxs, gbufx, gbufy, comments);
      if(allocated) delete[] gbufx;
      if(allocated) delete[] gbufy;
      break;
    case PARTITIONED_OUTPUT:
      status=quoddy_savec2(filename, nnodes, bufx, bufy, comments);
      break;
    default:
      TRAP_ERR_EXIT(-1, "illicit output mode\n");
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void read_boundarycode_obsolete(const char *filename, mesh_t & mesh, int level)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------

  construct boundary codes from BEL file

  Input:
   -  mesh.vertices cross tables needed
   -  mesh.triangles cross tables needed


----------------------------------------------------------------------*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int k, l, e;
  int dum1, dum2, n, nd1, nd2, code, n1, n2;
  int nitems, found;
  int nedges, count = 0;
  char char_dum[500],*msg;
  FILE *in;
  edge_t *edges;

  nedges = 0;

  edges  = mesh.edges;
  nedges = mesh.nedges;

/**----------------------------------------------------------------------
  initialize edges' code to UNDEFINED VALUE */
  for(n = 0; n < nedges; n++)
    edges[n].code = MESH_UNDEFINED_EDGE;

  in = fopen(filename, "r");
  if(in==0) {
    msg=new char[1024];
    TRAP_ERR_EXIT(-1,"read_boundarycode failed (file opening: %s)",filename);
    }

  fgets(char_dum, n, in);
  fgets(char_dum, n, in);

  count = 0;

  while(!feof(in)) {
    nitems = fscanf(in, "%d %d %d %d %d", &dum1, &nd1, &nd2, &dum2, &code);
    if(nitems != 5)
      goto done;
/**----------------------------------------------------------------------
    C-array from fortan numbering*/
    nd1--;
    nd2--;
    count++;
/**----------------------------------------------------------------------
    find corresponding triangle edge */
    found = 0;
    for(k = 0; k < mesh.vertices[nd1].nelmts; k++) {
      e = mesh.vertices[nd1].elmts[k];
      for(l = 0; l < 3; l++) {
        n = mesh.triangles[e].edges[l];
        n1 = edges[n].extremity[0];
        n2 = edges[n].extremity[1];
        if((n1 == nd1) && (n2 == nd2)) {
          edges[n].code = code;
          found = 1;
          break;
          }
        if((n1 == nd2) && (n2 == nd1)) {
          edges[n].code = code;
          found = 1;
          break;
          }
        }
      }
    if(found != 1) {
      msg=new char[256];
      TRAP_ERR_EXIT(-1,"boundary error at edge %d, nodes %d %d, code %d --- edge not found\n", dum1, nd1, nd2, code);
      }
    }

done:

  fclose(in);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void read_boundarycode(const char *filename, mesh_t & mesh, int level)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------

  construct boundary codes from BEL file

  Input:
   -  mesh.vertices cross tables needed
   -  mesh.triangles cross tables needed


----------------------------------------------------------------------*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int dum1, dum2, n, nd1, nd2, code;
  int nitems, found;
  int nedges, count = 0;
  char char_dum[500],*msg;
  FILE *in;
  edge_t *edges;

  nedges = 0;

  edges  = mesh.edges;
  nedges = mesh.nedges;

/**----------------------------------------------------------------------
  initialize edges' code to UNDEFINED VALUE */
  for(n = 0; n < nedges; n++)
    edges[n].code = MESH_UNDEFINED_EDGE;

  in = fopen(filename, "r");
  if(in==0) {
    msg=new char[1024];
    TRAP_ERR_EXIT(-1,"read_boundarycode failed (file opening: %s)",filename);
    }

  fgets(char_dum, n, in);
  fgets(char_dum, n, in);

  count = 0;

  while(!feof(in)) {
    nitems = fscanf(in, "%d %d %d %d %d", &dum1, &nd1, &nd2, &dum2, &code);
    if(nitems != 5)
      goto done;
/**----------------------------------------------------------------------
    C-array from fortan numbering*/
    nd1--;
    nd2--;
    count++;
/**----------------------------------------------------------------------
    find corresponding triangle edge */
    found = 0;
    n=fe_isedge(mesh, nd1,nd2);
    if(n!=-1) {
      edges[n].code=code;
      found=1;
      }
    if(found != 1) {
      msg=new char[256];
      TRAP_ERR_EXIT(-1,"boundary error at edge %d, nodes %d %d, code %d --- edge not found\n", dum1, nd1, nd2, code);
      }
    }

done:

  fclose(in);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_boundarycode(const char *filename, mesh_t & mesh, int level)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,n1,n2/*,status*/;

  if(mesh.origin==0) {
    fe_read_boundarycode(filename, mesh, level);
    }
  else {
    fe_read_boundarycode(filename, *mesh.origin, level);
    for(n=0;n<mesh.nedges;n++) {
      mesh.edges[n].code=mesh.origin->edges[mesh.edges[n].origin].code;
      }
    for(n = 0; n < mesh.nedges; n++) {
      n1 = mesh.edges[n].extremity[0];
      n2 = mesh.edges[n].extremity[1];
      if((mesh.edges[n].nshared == 1) && (mesh.edges[n].code == MESH_UNDEFINED_EDGE)) {
        mesh.edges[n].code = MESH_PARTITION_NODE;
        }
      }
    }
//  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void communications_receive(mesh_t mesh, int **zTable, int **uTable, connexion_t connexions)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//   int k,l,m,n,status;
//   int count,zdim,udim;
//   int *cardinal;
//   double *buffer;
//   int CPUid;
//
// // MPI_Send(buffer,count,MPI_DOUBLE,dest,tag,comm);
//   status=paire_dimension(gFEmesh[0],&zdim,&udim);
//
//   for (CPUid=0;CPUid<connexions.nCPUs;CPUid++) {
// //  status=MPI_Receive(buffer,count,MPI_DOUBLE,connexions.CPUs[CPUid],TAG_ELEVATION,MPI_COMM_WORLD);
//     for (l=0;l<cardinal[CPUid];l++) {
//       n=zTable[CPUid][l];
//       gstate2D[0][2].H[n]=buffer[l];
//       }
// //  status=MPI_Receive(buffer,count,MPI_DOUBLE,connexions.CPUs[CPUid],TAG_VELOCITIES,MPI_COMM_WORLD);
//     for (l=0;l<cardinal[CPUid];l++) {
//       n=uTable[CPUid][l];
//       gstate2D[0][2].u[n]=buffer[l];
//       }
//     }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int import_finalize_LGP1xLGP1(mesh_t mesh, exchange_t exchange)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Build importation tables from UniversalImportTables

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int k,p,n,nn;
  int *origin;

#ifdef HAVE_MPI
  int status;
  
  status = MPI_Barrier (  MPI_COMM_WORLD);
//  cout << "import_finalize_LGP1xLGP1 : exchange.nCPUs=" << exchange.nCPUs << endl ;
//  cout << "import_finalize_LGP1xLGP1 : exchange.rank=" << exchange.CPUid << endl ;
//  status = MPI_Barrier (  MPI_COMM_WORLD);
  for (p=0;p<exchange.nCPUs;p++) {
    if(p==exchange.CPUid) continue;
    cout << "import_finalize_LGP1xLGP1, cpu=" << exchange.CPUid << " : rank " << p << ", inport count="<< exchange.zImportCount[p] << endl;
    }
  status = MPI_Barrier (  MPI_COMM_WORLD);
#endif

  origin=new int[mesh.nvtxs];
  for (n=0;n<mesh.nvtxs;n++) {
    origin[n]=mesh.vertices[n].origin;
    }

  for (p=0;p<exchange.nCPUs;p++) {
    if(p==exchange.CPUid) continue;
    exchange.zImportTable[p] = new int[exchange.zImportCount[p]];
    for (k=0;k<exchange.zImportCount[p];k++) {
/**----------------------------------------------------------------------------
      nn is node index in  master mesh*/
      nn=exchange.zUniversalImportTable[p][k];
/**----------------------------------------------------------------------------
      find out local node corresponding to master mesh */
      n=vpos(nn,origin,mesh.nvtxs);
      exchange.zImportTable[p][k]=n;
      // exchange.zImportTable[p][k]=mesh.origin->vertices[n].child;  // MODIF CN 2011-02-11
      // exchange.zImportTable[p][k]=n;  // MODIF CN 2011-02-11 -----> faux identique a dessus
//      cout << "n=" << n << endl;
      }
/**----------------------------------------------------------------------------
    LGP1xLGP1, u tables are identical to z tables */
    exchange.uImportTable[p]=exchange.zImportTable[p];
    }

  delete[] origin;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int import_finalize_GENERIC(mesh_t mesh, int discretisation, exchange_t exchange, char target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Build importation tables from UniversalImportTables

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int k,p,n,nn,zdim;
  int *origin;
  discretisation_t descriptor=get_descriptor(mesh, discretisation);

  zdim=descriptor.nnodes;

#ifdef HAVE_MPI
  int status;
  
  status = MPI_Barrier(MPI_COMM_WORLD);
  for (p=0;p<exchange.nCPUs;p++) {
    if(p==exchange.CPUid) continue;
    cout << "import_finalize_GENERIC, cpu=" << exchange.CPUid << " : rank " << p << ", inport count="<< exchange.zImportCount[p] << endl;
    }
  status = MPI_Barrier(MPI_COMM_WORLD);
#endif

  origin=new int[zdim];
  for (n=0;n<zdim;n++) {
    origin[n]=descriptor.nodes[n].globalized;
    }

  switch(target) {
    case 'z':
  for (p=0;p<exchange.nCPUs;p++) {
    if(p==exchange.CPUid) continue;
    exchange.zImportTable[p] = new int[exchange.zImportCount[p]];
    for (k=0;k<exchange.zImportCount[p];k++) {
/**----------------------------------------------------------------------------
      nn is node index in  master mesh*/
      nn=exchange.zUniversalImportTable[p][k];
/**----------------------------------------------------------------------------
      find out local node corresponding to master mesh */
      n=vpos(nn,origin,zdim);
      exchange.zImportTable[p][k]=n;
      }
    }
      break;
    case 'u':
  for (p=0;p<exchange.nCPUs;p++) {
    if(p==exchange.CPUid) continue;
    exchange.uImportTable[p] = new int[exchange.uImportCount[p]];
    for (k=0;k<exchange.uImportCount[p];k++) {
      nn=exchange.uUniversalImportTable[p][k];
      n=vpos(nn,origin,zdim);
      exchange.uImportTable[p][k]=n;
      }
    }
      break;
    default:
      TRAP_ERR_EXIT(-1, "illegal paire\n");
      break;
    }

  delete[] origin;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solved_init_LGP1xLGP1(mesh_t mesh, exchange_t exchange)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Build solved tables from zImportTable

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int k,p,n;
  int *imported=new int[mesh.nvtxs];
  int count;
//   char filesolved[14];
//   int nCPUs,CPU_ID;
//   FILE *out;

/**----------------------------------------------------------------------------
  count solved nodes (i.e. not imported) */
  count=mesh.nvtxs;
  for (n=0;n<mesh.nvtxs;n++) {
    imported[n]=0;
    }
  for (p=0;p<exchange.nCPUs;p++) {
    if(p==exchange.CPUid) continue;
    for (k=0;k<exchange.zImportCount[p];k++) {
      n=exchange.zImportTable[p][k];
      imported[n]=1;
      count--;
      }
    }

  p=exchange.CPUid;

  if(count==0) {
    printf( "processor %d solves zero nodes\n",p);
    TRAP_ERR_EXIT(-1, "parellel tables initialisation failure\n");
    }

  exchange.zSolvedTable[p]=new int[count];
  exchange.zUniversalSolvedTable[p]=new int[count];
  exchange.zSolvedCount[p]=count;

#ifdef CHECK_PARALLEL
  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  cout << "nsolved = " << count << endl;
  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
#endif
  
  count=0;
  for (n=0;n<mesh.nvtxs;n++) {
    if(imported[n]==0) {
/**----------------------------------------------------------------------------
      n is solved */
      exchange.zSolvedTable[p][count]=n;
      exchange.zUniversalSolvedTable[p][count]=mesh.vertices[n].origin;
      count++;
      }
    }

/**----------------------------------------------------------------------------
  LGP1xLGP1, u tables are identical to z tables */
  exchange.uSolvedTable[p]=exchange.zSolvedTable[p];
  exchange.uUniversalSolvedTable[p]=exchange.zUniversalSolvedTable[p];
  exchange.uSolvedCount[p]=exchange.zSolvedCount[p];

//  sprintf(filesolved,"solved.%d.%d.txt",nCPUs,CPU_ID);
//  cout << filesolved << "%%%%%%%%%%%%%%%%" << endl;

  //out = fopen(filesolved, "w");
  //fprintf(out, "%d  \n",exchange.uSolvedCount[p] );

//   printf( "%d  \n",exchange.zSolvedCount[p] );
//   for (n=0; n<exchange.uSolvedCount[p];n++) {
//     //  printf(out, "%d , universal num: %d\n", exchange.uSolvedTable[p][n],exchange.uUniversalSolvedTable[p][n]);
//     printf("%d , universal num: %d,%d\n", exchange.zSolvedTable[p][n],exchange.zUniversalSolvedTable[p][n],
// 	   exchange.zUniversalSolvedTable[p][n]+1);
//     }

  //fclose(out);
//#ifdef HAVE_MPI
//      MPI_Barrier (  MPI_COMM_WORLD);
//      cout << "Fin ecriture fichier...";
//      cout << endl;
//#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solved_init_GENERIC(mesh_t mesh, exchange_t exchange, int paire)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Build solved tables from zImportTable
  
  Actually this is redundant with distribution_t structure that comes
  attached to discretisation_t...
  
  \todo To be rationalized...

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int k,p,n;
  int *imported=new int[mesh.nvtxs];
  int count;

  return(0);
  
/**----------------------------------------------------------------------------
  count solved nodes (i.e. not imported) */
  count=mesh.nvtxs;
  for (n=0;n<mesh.nvtxs;n++) {
    imported[n]=0;
    }
  for (p=0;p<exchange.nCPUs;p++) {
    if(p==exchange.CPUid) continue;
    for (k=0;k<exchange.zImportCount[p];k++) {
      n=exchange.zImportTable[p][k];
      imported[n]=1;
      count--;
      }
    }

  p=exchange.CPUid;

  if(count==0) {
    printf( "processor %d solves zero nodes\n",p);
    TRAP_ERR_EXIT(-1, "parellel tables initialisation failure\n");
    }

  exchange.zSolvedTable[p]=new int[count];
  exchange.zUniversalSolvedTable[p]=new int[count];
  exchange.zSolvedCount[p]=count;

#ifdef CHECK_PARALLEL
  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  cout << "nsolved = " << count << endl;
  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
#endif
  
  count=0;
  for (n=0;n<mesh.nvtxs;n++) {
    if(imported[n]==0) {
/**----------------------------------------------------------------------------
      n is solved */
      exchange.zSolvedTable[p][count]=n;
      exchange.zUniversalSolvedTable[p][count]=mesh.vertices[n].origin;
      count++;
      }
    }

/**----------------------------------------------------------------------------
  LGP1xLGP1, u tables are identical to z tables */
  exchange.uSolvedTable[p]=exchange.zSolvedTable[p];
  exchange.uUniversalSolvedTable[p]=exchange.zUniversalSolvedTable[p];
  exchange.uSolvedCount[p]=exchange.zSolvedCount[p];

}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   exchange_t communications_init(mesh_t mesh, connexion_t eConnexions, int CPUid, int nCPUs, int gSequentialPaire2D, int verbose)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//   contsruct import/export tables for parallel computing
//
//   mesh is partioned mesh treated by cpu CPUid
//   computational discretisation dependent
//
//   A processor (A) needs to know:
//
//     1- which processors (B,C,...) to address to get values
//        * must provide index (in B) of node n (in A)
//        * count
//
//     2- which processor to address to send values
//
//   Information should be symetrical, so only step 1 or 2 needs to be
//   actually computed, then MPI exhanges might be used to complete the tables
//
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
// {
//   int k,l,m,n,p;
//   int status;
//   int *buffer, *bufferS, **bufferR,count;
//   int tagZcount, tagZnode;
//   exchange_t exchange(CPUid, nCPUs);
//   int tnbSend[nCPUs][nCPUs];
//   int nSend, nRecv, tagz;
//
// #ifdef HAVE_MPI
//   MPI_Request request[nCPUs];
//   MPI_Status mpistatus[nCPUs];
//   MPI_Request *requestZcount, *requestZnode;
// #endif
//
//
// //   printf(SEPARATOR_1);
//    printf("\n initialize communication betwen CPU %d and other CPUs\n", CPUid);
//
// /**----------------------------------------------------------------------------
//   step 1 : build export tables */
//   switch(gSequentialPaire2D) {
//     case LGP1xLGP1:
// #ifdef HAVE_MPI
//       status = MPI_Barrier (MPI_COMM_WORLD);
// #endif
//       status=export_init_LGP1xLGP1(mesh, eConnexions, CPUid, exchange.zExportTable, exchange.zUniversalExportTable, exchange.zExportCount,
//                                                              exchange.uExportTable, exchange.uUniversalExportTable, exchange.uExportCount);
// #ifdef HAVE_MPI
//       status = MPI_Barrier (MPI_COMM_WORLD);
// #endif
//       break;
//
//     case NCP1xLGP1:
//       status=export_init_LGP1(mesh, eConnexions, CPUid, exchange.zExportTable, exchange.zUniversalExportTable, exchange.zExportCount);
//       status=export_init_GENERIC(mesh, eConnexions, CPUid, NCP1, exchange.uExportTable, exchange.uUniversalExportTable, exchange.uExportCount);
//
// //      status=export_init_LGP1(mesh, eConnexions, CPUid, exchange.zExportTable, exchange.zUniversalExportTable, exchange.zExportCount);
// //      status=export_init_LGP1(mesh, eConnexions, CPUid, exchange.uExportTable, exchange.uUniversalExportTable, exchange.uExportCount);
//
// /**----------------------------------------------------------------------------
//       something buggy in export_init, use the wrong one for testing*/
// //       status=export_init_LGP1xLGP1(mesh, eConnexions, CPUid, exchange.zExportTable, exchange.zUniversalExportTable, exchange.zExportCount,
// //                                                              exchange.uExportTable, exchange.uUniversalExportTable, exchange.uExportCount);
//       break;
//
//     case NCP1xLGP0:
//       TRAP_ERR_EXIT(-1, "paire not yet implemented for parallel computing\n");
//       break;
//
//     case NCP1xLGP2:
//       status=export_init_GENERIC(mesh, eConnexions, CPUid, LGP2, exchange.zExportTable, exchange.zUniversalExportTable, exchange.zExportCount);
//       status=export_init_GENERIC(mesh, eConnexions, CPUid, NCP1, exchange.uExportTable, exchange.uUniversalExportTable, exchange.uExportCount);
//       break;
//
//     case DNP1xLGP2:
//       status=export_init_GENERIC(mesh, eConnexions, CPUid, LGP2, exchange.zExportTable, exchange.zUniversalExportTable, exchange.zExportCount);
//       status=export_init_GENERIC(mesh, eConnexions, CPUid, DNP1, exchange.uExportTable, exchange.uUniversalExportTable, exchange.uExportCount);
//       break;
//
//     default:
//       TRAP_ERR_EXIT(-1, "illegal paire\n");
//       break;
//     }
//
// /**----------------------------------------------------------------------------
//   send information on Z exchanges*/
//
// /*-----------------------------------------------------------------------------
//   first send the number of exported nodes per CPU */
//   for (p=0;p<eConnexions.nCPUs;p++) {
//     tnbSend[CPUid][p] = exchange.zExportCount[p];
//     }
//
// #ifdef HAVE_MPI
//   status = MPI_Allgather(tnbSend[CPUid], nCPUs,MPI_INTEGER,  /* send buf,count,type */
//                          tnbSend, nCPUs,MPI_INTEGER,         /* recv buf,count,type */
//                          MPI_COMM_WORLD);                    /* comm,flag           */
// #endif
//   for (p=0;p<eConnexions.nCPUs;p++) {
//     exchange.zImportCount[p]=  tnbSend[p][CPUid];
//     }
//
// /*-----------------------------------------------------------------------------
//   send list of nodes exported from CPUid to the other CPUs */
//   nSend = 0;
//   for (p=0;p<eConnexions.nCPUs;p++) {
//     // Send/recv buff
//     if(exchange.zExportCount[p]==0) continue;
//     bufferS=new int[exchange.zExportCount[p]];
//     if(verbose) cout << "Send Nb =" << exchange.zExportCount[p] << " vers " << p << endl;
//     for(k=0;k<exchange.zExportCount[p];k++) {
//       n=exchange.zUniversalExportTable[p][k];
//       bufferS[k]=n;
//       }
// /*-----------------------------------------------------------------------------
//     send list of nodes exported to p */
//     count=exchange.zExportCount[p];
// /**----------------------------------------------------------------------------
//     ??? */
//     tagz = CPUid*(eConnexions.nCPUs-1) + p;
// #ifdef HAVE_MPI
//     status=MPI_Issend(bufferS, exchange.zExportCount[p], MPI_INTEGER,
// 		      p,  tagz, MPI_COMM_WORLD, &request[nSend]);
// #endif
//     nSend ++;
//     delete[] bufferS;
//     }
//
// /*-----------------------------------------------------------------------------
//   step 2 : deduce import tables*/
//
//   nRecv = 0;
//   bufferR=new int* [eConnexions.nCPUs];
//
//   for (p=0;p<eConnexions.nCPUs;p++) {
// /**----------------------------------------------------------------------------
//     ??? */
//     exchange.zImportCount[p]=  tnbSend[p][CPUid];
//     if(tnbSend[p][CPUid] == 0) continue;
//     bufferR[p] = new int [exchange.zImportCount[p]];
//     tagz = p*(eConnexions.nCPUs-1) + CPUid;
// /**----------------------------------------------------------------------------
//     receive list of nodes exported from the other CPUs to CPUid*/
// #ifdef HAVE_MPI
//     status = MPI_Irecv(bufferR[p], exchange.zImportCount[p], MPI_INTEGER, p, tagz, MPI_COMM_WORLD, &request[nRecv]);
// #endif
//     nRecv ++;
//     }
//
// /**----------------------------------------------------------------------------
//   wait exchange completion*/
// #ifdef HAVE_MPI
//   status = MPI_Waitall(nRecv, &request[0], &mpistatus[0]);
// #endif
//
//   for (p=0; p<eConnexions.nCPUs; p++) {
//     if(exchange.zImportCount[p] == 0) continue;
//     exchange.zUniversalImportTable[p]=new int[exchange.zImportCount[p]];
//     for(k=0;k<exchange.zImportCount[p];k++) {
//       exchange.zUniversalImportTable[p][k]=bufferR[p][k];
//       }
//     delete[] bufferR[p];
//     }
//
//   delete[] bufferR;
//
// /**----------------------------------------------------------------------------
//   send information on U exchanges*/
//
//   for (p=0;p<eConnexions.nCPUs;p++) {
//     tnbSend[CPUid][p] = exchange.uExportCount[p];
//     }
//
// #ifdef HAVE_MPI
//   status = MPI_Allgather(tnbSend[CPUid], nCPUs, MPI_INTEGER,     /* send buf,count,type */
//                          tnbSend,        nCPUs, MPI_INTEGER,     /* recv buf,count,type */
//                          MPI_COMM_WORLD);                        /* comm,flag           */
// #endif
//   for (p=0;p<eConnexions.nCPUs;p++) {
//     exchange.uImportCount[p]=  tnbSend[p][CPUid];
//     }
//   // Send
//   nSend = 0;
//   for (p=0;p<eConnexions.nCPUs;p++) {
//     // Send/recv buff
//     if(exchange.uExportCount[p]==0) continue;
//     bufferS=new int[exchange.uExportCount[p]];
//     for(k=0;k<exchange.uExportCount[p];k++) {
//       n=exchange.uUniversalExportTable[p][k];
//       bufferS[k]=n;
//       }
//     count=exchange.uExportCount[p];
//     tagz = CPUid*(eConnexions.nCPUs-1) + p;
// #ifdef HAVE_MPI
//     status=MPI_Issend(&bufferS[0], exchange.uExportCount[p]   , MPI_INTEGER,
// 		     p,  tagz, MPI_COMM_WORLD, &request[nSend]);
// #endif
//     nSend ++;
//     delete[] bufferS;
//     }
//   nRecv = 0;
//   bufferR=new int* [eConnexions.nCPUs];
//   // RECV
//   for (p=0;p<eConnexions.nCPUs;p++) {
//     exchange.uImportCount[p]=  tnbSend[p][CPUid];
//     if(tnbSend[p][CPUid] == 0) continue;
//     bufferR[p] = new int [exchange.uImportCount[p]];
//     tagz = p*(eConnexions.nCPUs-1) + CPUid;
// #ifdef HAVE_MPI
//     status = MPI_Irecv(&bufferR[p][0], exchange.uImportCount[p], MPI_INTEGER,
//     		       p,  tagz, MPI_COMM_WORLD, &request[nRecv]);
// #endif
//     nRecv ++;
//     }
//   //Wait
// #ifdef HAVE_MPI
//   status = MPI_Waitall(nRecv, &request[0], &mpistatus[0]);
//   status = MPI_Barrier (  MPI_COMM_WORLD);
// #endif
//   if(verbose) cout << CPUid << " Salle d'attente "<< nRecv << endl;
//   for (p=0; p<eConnexions.nCPUs; p++) {
//     if(verbose) cout << CPUid << " Nelts recu = " << exchange.uImportCount[p] << endl;
//     if(exchange.uImportCount[p] == 0) continue;
//     exchange.uUniversalImportTable[p]=new int[exchange.uImportCount[p]];
//     if(verbose) cout << CPUid << " liste de recu : " ;
//     for(k=0;k<exchange.uImportCount[p];k++) {
//       if(verbose) cout << bufferR[p][k] << ",";
//       exchange.uUniversalImportTable[p][k]=bufferR[p][k];
//       }
//     delete[] bufferR[p];
//     if(verbose) cout << endl;
//     }
//   delete[] bufferR;
//
// #ifdef HAVE_MPI
//   status = MPI_Barrier (  MPI_COMM_WORLD);
// #ifdef DEBUG_PARALLEL
//   cout << "synchronize before completing import information " << endl;
//   status = MPI_Barrier (  MPI_COMM_WORLD);
// #endif
// #endif
//
// /**----------------------------------------------------------------------------
//   complete import information */
//   switch(gSequentialPaire2D) {
//     case LGP1xLGP1:
//       status=import_finalize_LGP1xLGP1(mesh, exchange);
//       break;
//
//     case NCP1xLGP1:
//       status=import_finalize_GENERIC(mesh, LGP1, exchange, 'z');
//       status=import_finalize_GENERIC(mesh, NCP1, exchange, 'u');
//       break;
//
//     case NCP1xLGP2:
//       status=import_finalize_GENERIC(mesh, LGP2, exchange, 'z');
//       status=import_finalize_GENERIC(mesh, NCP1, exchange, 'u');
//       break;
//
//     case DNP1xLGP2:
//       status=import_finalize_GENERIC(mesh, LGP2, exchange, 'z');
//       status=import_finalize_GENERIC(mesh, DNP1, exchange, 'u');
//       break;
//
//     default:
//       TRAP_ERR_EXIT(-1, "illegal paire\n");
//       break;
//     }
//
// #ifdef HAVE_MPI
//   status = MPI_Barrier (  MPI_COMM_WORLD);
// #ifdef DEBUG_PARALLEL
//   cout << "synchronize after completing import information " << endl;
//   status = MPI_Barrier (  MPI_COMM_WORLD);
// #endif
// #endif
//
// /**----------------------------------------------------------------------------
//   complete solved information */
//   switch(gSequentialPaire2D) {
//     case LGP1xLGP1:
//       status=solved_init_LGP1xLGP1(mesh, exchange);
//       break;
//
//     case NCP1xLGP1:
//       status=solved_init_GENERIC(mesh, exchange,gSequentialPaire2D);
//       break;
//
//     case NCP1xLGP2:
//       status=solved_init_GENERIC(mesh, exchange,gSequentialPaire2D);
//       break;
//
//     case DNP1xLGP2:
//       status=solved_init_GENERIC(mesh, exchange,gSequentialPaire2D);
//       break;
//
//     default:
//       TRAP_ERR_EXIT(-1, "illegal paire\n");
//       break;
//     }
//
// #ifdef HAVE_MPI
//   status = MPI_Barrier (  MPI_COMM_WORLD);
// #ifdef DEBUG_PARALLEL
//   cout << "STOP apres appel complete export information " << endl;
//   status = MPI_Barrier (  MPI_COMM_WORLD);
// #endif
// #endif
//
// #ifdef DEBUG_PARALLEL
//   status=communications_summary(exchange);
// #endif
// #ifdef HAVE_MPI
//   status = MPI_Barrier (  MPI_COMM_WORLD);
// #endif

// next:
//   return(exchange);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int communications_summary(exchange_t exchange)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int k,p,q;
  FILE *summary;
  char filename[64];

/**----------------------------------------------------------------------------
  MPI safe */
  sprintf(filename,"summary.%d.%d.txt",exchange.nCPUs,exchange.CPUid);
  summary=fopen(filename,"w");

  p=exchange.CPUid;

//  fprintf(summary,SEPARATOR_1);
  for(q=0;q<exchange.nCPUs;q++) {
    if(p==q) continue;
//     fprintf(summary,SEPARATOR_1);
    fprintf(summary,"zExport tables from CPU %d to CPU %d: count = %d \n",p,q,exchange.zExportCount[q]);
    fprintf(summary,"local    : ");
    for(k=0;k<exchange.zExportCount[q];k++) {
      fprintf(summary,"%3d ",exchange.zExportTable[q][k]);
      }
    fprintf(summary,"\n");
    fprintf(summary,"universal : ");
    for(k=0;k<exchange.zExportCount[q];k++) {
      fprintf(summary,"%3d ",exchange.zUniversalExportTable[q][k]);
      }
    fprintf(summary,"\n");
    }

//  fprintf(summary,SEPARATOR_1);
  for(q=0;q<exchange.nCPUs;q++) {
    if(p==q) continue;
    fprintf(summary,"zImport tables from CPU %d to CPU %d: count = %d \n",p,q,exchange.zImportCount[q]);
    fprintf(summary,"local    : ");
    for(k=0;k<exchange.zImportCount[q];k++) {
      fprintf(summary,"%3d ",exchange.zImportTable[q][k]);
      }
    fprintf(summary,"\n");
    fprintf(summary,"universal : ");
    for(k=0;k<exchange.zImportCount[q];k++) {
      fprintf(summary,"%3d ",exchange.zUniversalImportTable[q][k]);
      }
    fprintf(summary,"\n");
    }

  fclose(summary);

}



