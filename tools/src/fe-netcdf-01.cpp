
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief finite element poc-netcdf definitions
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "fe.def"
#include "fe.h"

#include "map.h"
#include "rutin.h"
#include "geo.h"
#include "poc-time.h"
#include "netcdf-proto.h"
#include "poc-netcdf.def"
#include "functions.h"
#include "zapper.h"

#define RANK_connectivity 2
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_readdiscretisation_TUGO(const char *filename,mesh_t *mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///read TUGO mesh from NetCDF file given the discretisation
/**
\param filename
\param *mesh mesh initialised with fe_readgeometry()
\param discretisation
\returns NC_NOERR on success or the NetCDF error status on failure
*/
/*----------------------------------------------------------------------------*/
{
  int n,k,status, ncid;
  int varid;
  int *v;
  int nmax=-1;
  bool completed=false;

  status=nc_open(filename,0,&ncid);
  if(status != NC_NOERR) goto error;

  switch(discretisation) {
    case LGP0:
/* *----------------------------------------------------------------------------
      T-UGOM, LGP0 nodes discretisation*/
      status=nc_inq_varid (ncid,"C-LGP0",&varid);
      if(status==NC_NOERR) {
        v=(int *) malloc(mesh->ntriangles*sizeof(int));
        status=nc_get_var_int   (ncid,varid,v);
        mesh->LGP0descriptor.NIbE=new int*[mesh->ntriangles];
        for (n=0;n<mesh->ntriangles;n++) {
          mesh->LGP0descriptor.NIbE[n]=new int[1];
          for(k=0;k<1;k++) {
            mesh->LGP0descriptor.NIbE[n][k]=v[n+k];
            }
          }
        free(v);
        mesh->LGP0descriptor.nnodes=mesh->ntriangles;
        mesh->LGP0descriptor.nnpe=1;
        mesh->LGP1descriptor.type=LGP0;
        mesh->LGP1descriptor.nodes=new node_t[mesh->LGP0descriptor.nnodes];
        completed=true;
        }
      break;

    case LGP1:
/* *----------------------------------------------------------------------------
      T-UGOM, LGP1 nodes discretisation*/
      status=nc_inq_varid (ncid,"C-LGP1",&varid);
      if(status==NC_NOERR) {
        v=(int *) malloc(mesh->ntriangles*3*sizeof(int));
        status=nc_get_var_int   (ncid,varid,v);
        mesh->LGP1descriptor.NIbE=new int*[mesh->ntriangles];
        for (n=0;n<mesh->ntriangles;n++) {
          mesh->LGP1descriptor.NIbE[n]=new int[3];
          for(k=0;k<3;k++) {
            mesh->LGP1descriptor.NIbE[n][k]=v[n*3+k];
            }
          }
        free(v);
        mesh->LGP1descriptor.nelmts=mesh->ntriangles;
        mesh->LGP1descriptor.nnodes=mesh->nvtxs;
        mesh->LGP1descriptor.nnpe=3;
        mesh->LGP1descriptor.type=LGP1;
        mesh->LGP1descriptor.nodes=new node_t[mesh->LGP1descriptor.nnodes];
        for(size_t m = 0; m < mesh->ntriangles; m++) {
          int nn;
          for(k=0;k<3;k++) {
            nn=mesh->triangles[m].vertex[k];
            n=mesh->LGP1descriptor.NIbE[m][k];
            mesh->LGP1descriptor.nodes[n].lon=mesh->vertices[nn].lon;
            mesh->LGP1descriptor.nodes[n].lat=mesh->vertices[nn].lat;
            }
          }
        completed=true;
        }
      break;

    case DGP1:
/* *----------------------------------------------------------------------------
      T-UGOM, DGP1 nodes discretisation*/
      status=nc_inq_varid (ncid,"C-DGP1",&varid);
      if(status==NC_NOERR) {
        v=(int *) malloc(mesh->ntriangles*3*sizeof(int));
        status=nc_get_var_int   (ncid,varid,v);
        mesh->DGP1descriptor.NIbE=new int*[mesh->ntriangles];
        for (n=0;n<mesh->ntriangles;n++) {
          mesh->DGP1descriptor.NIbE[n]=new int[3];
          for(k=0;k<3;k++) {
            mesh->DGP1descriptor.NIbE[n][k]=v[n*3+k];
            }
          }
        free(v);
        mesh->DGP1descriptor.nelmts=mesh->ntriangles;
        mesh->DGP1descriptor.nnodes=mesh->ntriangles*3;
        mesh->DGP1descriptor.nnpe=3;
        completed=true;
        }
      break;

    case NCP1:
/* *----------------------------------------------------------------------------
      T-UGOM, NCP1 nodes discretisation*/
      status=nc_inq_varid (ncid,"C-NCP1",&varid);
      if(status==NC_NOERR) {
        v=(int *) malloc(mesh->ntriangles*3*sizeof(int));
        status=nc_get_var_int   (ncid,varid,v);
        mesh->NCP1descriptor.NIbE=new int*[mesh->ntriangles];
        for (n=0;n<mesh->ntriangles;n++) {
          mesh->NCP1descriptor.NIbE[n]=new int[3];
          for(k=0;k<3;k++) {
            mesh->NCP1descriptor.NIbE[n][k]=v[n*3+k];
            nmax=max(nmax,mesh->NCP1descriptor.NIbE[n][k]);
            }
          }
        free(v);
        mesh->NCP1descriptor.nelmts=mesh->ntriangles;
        mesh->NCP1descriptor.nnodes=nmax+1;
        mesh->NCP1descriptor.nnpe=3;
        mesh->NCP1descriptor.type=NCP1;

/* *----------------------------------------------------------------------------
        create NCP1 nodes NIbE list (for each element)*/
        mesh->NCP1descriptor.nodes=new node_t[mesh->NCP1descriptor.nnodes];
        for(size_t m = 0; m < mesh->ntriangles; m++) {
          int nn;
          for(k=0;k<3;k++) {
            nn=mesh->triangles[m].edges[k];
            n=mesh->NCP1descriptor.NIbE[m][2*k+1];
            mesh->NCP1descriptor.nodes[n].lon=mesh->edges[nn].lon;
            mesh->NCP1descriptor.nodes[n].lat=mesh->edges[nn].lat;
            }
          }
        completed=true;
        }
      break;

    case DNP1:
/* *----------------------------------------------------------------------------
      T-UGOM, DNP1 nodes discretisation*/
      status=nc_inq_varid (ncid,"C-DNP1",&varid);
      if(status==NC_NOERR) {
        v=new int[mesh->ntriangles*3];
        status=nc_get_var_int   (ncid,varid,v);
        mesh->DNP1descriptor.NIbE=new int*[mesh->ntriangles];
        for (n=0;n<mesh->ntriangles;n++) {
          mesh->DNP1descriptor.NIbE[n]=new int[3];
          for(k=0;k<3;k++) {
            mesh->DNP1descriptor.NIbE[n][k]=v[n*3+k];
            nmax=max(nmax,mesh->DNP1descriptor.NIbE[n][k]);
            }
          }
        delete[] v;
        mesh->DNP1descriptor.nelmts=mesh->ntriangles;
        mesh->DNP1descriptor.nnodes=nmax+1;
        mesh->DNP1descriptor.nnpe=3;
        mesh->DNP1descriptor.type=DNP1;

/* *----------------------------------------------------------------------------
        create DNP1 nodes NIbE list (for each element)*/
        mesh->DNP1descriptor.nodes=new node_t[mesh->DNP1descriptor.nnodes];
        for(size_t m = 0; m < mesh->ntriangles; m++) {
          int nn;
          for(k=0;k<3;k++) {
            nn=mesh->triangles[m].edges[k];
            n=mesh->DNP1descriptor.NIbE[m][k];
            mesh->DNP1descriptor.nodes[n].lon=mesh->edges[nn].lon;
            mesh->DNP1descriptor.nodes[n].lat=mesh->edges[nn].lat;
            }
          }
        completed=true;
        }
      break;

    case LGP2:
/* *----------------------------------------------------------------------------
      T-UGOM, LGP2 nodes discretisation*/
      status=nc_inq_varid (ncid,"C-LGP2",&varid);
      if(status==NC_NOERR) {
        v=(int *) malloc(mesh->ntriangles*6*sizeof(int));
        status=nc_get_var_int   (ncid,varid,v);
        mesh->LGP2descriptor.NIbE=new int*[mesh->ntriangles];
        for (n=0;n<mesh->ntriangles;n++) {
          mesh->LGP2descriptor.NIbE[n]=new int[6];
          for(k=0;k<6;k++) {
            mesh->LGP2descriptor.NIbE[n][k]=v[n*6+k];
            nmax=max(nmax,mesh->LGP2descriptor.NIbE[n][k]);
            }
          }
        free(v);
        mesh->LGP2descriptor.nelmts=mesh->ntriangles;
        mesh->LGP2descriptor.nnodes=nmax+1;
        mesh->LGP2descriptor.nnpe=6;
        mesh->LGP2descriptor.type=LGP2;

/* *----------------------------------------------------------------------------
        create LGP2 nodes NIbE list (for each element)*/
        mesh->LGP2descriptor.nodes=new node_t[mesh->LGP2descriptor.nnodes];
        for(size_t m = 0; m < mesh->ntriangles; m++) {
          int nn;
          for(k=0;k<3;k++) {
            nn=mesh->triangles[m].vertex[k];
            n=mesh->LGP2descriptor.NIbE[m][2*k];
            mesh->LGP2descriptor.nodes[n].lon=mesh->vertices[nn].lon;
            mesh->LGP2descriptor.nodes[n].lat=mesh->vertices[nn].lat;
            }
          for(k=0;k<3;k++) {
            nn=mesh->triangles[m].edges[k];
            n=mesh->LGP2descriptor.NIbE[m][2*k+1];
            mesh->LGP2descriptor.nodes[n].lon=mesh->edges[nn].lon;
            mesh->LGP2descriptor.nodes[n].lat=mesh->edges[nn].lat;
            }
          }
        completed=true;
        }
      break;
    }
    
  status = nc_close(ncid);
  
  if(completed) status= 0;
  else          status=-1;
  
  return(status);

error:
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_readdiscretisation_WWW3(const char *filename,mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, ncid;

  status=nc_open(filename,0,&ncid);
  if(status != NC_NOERR) goto error;
  status = nc_close(ncid);
  return(status);

error:
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_readdiscretisation(const char *filename,mesh_t *mesh, int fmt, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///read mesh from NetCDF file given the discretisation
/**
\param filename
\param *mesh
\param fmt grid format: 0. Returns NC_EINVAL if anything else.
\param discretisation
\returns NC_NOERR on success or the NetCDF error status on failure
*/
/*----------------------------------------------------------------------------*/
{
  int status=NC_EINVAL;
  
  switch(fmt) {
    case 0:
      ///If the format is 0, calls fe_readdiscretisation_TUGO()
      status=fe_readdiscretisation_TUGO(filename,mesh,discretisation);
      break;
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_e2edges (mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  build the neigh list from the element list*/
{
  int i,m,n;
  int n1,n2;

  for(n=0;n<mesh->nedges;n++) mesh->edges[n].nshared=0;

  for(m=0;m<mesh->ntriangles;m++)
    for(i=0;i<3;i++) {
      n=mesh->triangles[m].edges[i];
      mesh->edges[n].nshared+=1;
      n1=mesh->triangles[m].vertex[(i+1)%3];
      n2=mesh->triangles[m].vertex[(i+2)%3];
      mesh->edges[n].extremity[0]=n1;
      mesh->edges[n].extremity[1]=n2;
      mesh->edges[n].lon=0.5*(mesh->vertices[n1].lon+mesh->vertices[n2].lon);
      mesh->edges[n].lat=0.5*(mesh->vertices[n1].lat+mesh->vertices[n2].lat);
      }

  for(m=0;m<mesh->nquadrangles;m++)
    for(i=0;i<4;i++) {
      n=mesh->quadrangles[m].edges[i];
      mesh->edges[n].nshared+=1;
      n1=mesh->quadrangles[m].vertex[(i+1)%3];
      n2=mesh->quadrangles[m].vertex[(i+2)%3];
      mesh->edges[n].extremity[0]=n1;
      mesh->edges[n].extremity[1]=n2;
      mesh->edges[n].lon=0.5*(mesh->vertices[n1].lon+mesh->vertices[n2].lon);
      mesh->edges[n].lat=0.5*(mesh->vertices[n1].lat+mesh->vertices[n2].lat);
      }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> void fe_set_triangles_template(mesh_t *mesh,const T *corners,int index0=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///set triangles from corners buffer
/*----------------------------------------------------------------------------*/
{
  int n,k;//<triangles and corners indexes
  
  if(!mesh->triangles)
    mesh->triangles=new triangle_t[mesh->ntriangles];
  
  for (n=0;n<mesh->ntriangles;n++) {
    for(k=0;k<3;k++) {
      mesh->triangles[n].vertex[k]=corners[k+3*n]-index0;
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void fe_set_triangles(mesh_t *mesh,const int *corners,int index0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  fe_set_triangles_template(mesh,corners,index0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void fe_set_triangles(mesh_t *mesh,const double *corners,int index0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  fe_set_triangles_template(mesh,corners,index0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> void fe_set_quadrangles_template(mesh_t *mesh,const T *corners,int index0=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///set quadrangles from corners buffer
/*----------------------------------------------------------------------------*/
{
  int n,k;//<quadrangles and corners indexes
  
  if(!mesh->quadrangles) 
    mesh->quadrangles=new quadrangle_t[mesh->nquadrangles];
  
  for (n=0;n<mesh->nquadrangles;n++) {
    for(k=0;k<4;k++) {
      mesh->quadrangles[n].vertex[k]=corners[k+4*n]-index0;
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void fe_set_quadrangles(mesh_t *mesh,const int *corners,int index0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  fe_set_quadrangles_template(mesh,corners,index0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void fe_set_quadrangles(mesh_t *mesh,const double *corners,int index0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  fe_set_quadrangles_template(mesh,corners,index0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_readgeometry(const char *filename,mesh_t *mesh,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///read mesh from NetCDF file
/**
\param filename
\param *mesh
\param verbose default:1
\returns NC_NOERR on success or the NetCDF error status on failure
*/
/*----------------------------------------------------------------------------*/
{
  int n,k,status, ncid;
  int varid,dim;
  double *x=0,*y=0,*z=0,*h=0;
  int *v;
  variable_t details;
  cdfgbl_t fileinfo;

  
  status=nc_open(filename,0,&ncid);
  if(status) NC_TRAP_ERROR(return,status,1,"nc_open(\"%s\",0,) error",filename);
  
/* *------------------------------------------------------------------------
  inquire file informations */
  status= cdf_globalinfo(ncid, &fileinfo, 0);

/* *----------------------------------------------------------------------------
  T-UGOM*/

  dim=cdf_identify_dimension(fileinfo,"N");
  if(dim!=-1) mesh->nvtxs=fileinfo.dimension[dim].length;

  dim=cdf_identify_dimension(fileinfo,"M");
  if(dim!=-1) mesh->ntriangles=fileinfo.dimension[dim].length;

  dim=cdf_identify_dimension(fileinfo,"NQUADRANGLES");
  if(dim!=-1) mesh->nquadrangles=fileinfo.dimension[dim].length;

  dim=cdf_identify_dimension(fileinfo,"Q");
  if(dim!=-1) mesh->nquadrangles=fileinfo.dimension[dim].length;

  dim=cdf_identify_dimension(fileinfo,"E");
  if(dim!=-1) mesh->nedges=fileinfo.dimension[dim].length;

  dim=cdf_identify_dimension(fileinfo,"K");
  if(dim!=-1) mesh->nlayers=fileinfo.dimension[dim].length;

  dim=cdf_identify_dimension(fileinfo,"L");
  if(dim!=-1) mesh->nlevels=fileinfo.dimension[dim].length;

/* *----------------------------------------------------------------------------
  WW3*/
  dim=cdf_identify_dimension(fileinfo,"node");
  if(dim!=-1) mesh->nvtxs=fileinfo.dimension[dim].length;

  dim=cdf_identify_dimension(fileinfo,"element");
  if(dim!=-1) mesh->ntriangles=fileinfo.dimension[dim].length;

/* *----------------------------------------------------------------------------
  FVCOM*/
  dim=cdf_identify_dimension(fileinfo,"node");
  if(dim!=-1) mesh->nvtxs=fileinfo.dimension[dim].length;

  dim=cdf_identify_dimension(fileinfo,"nele");
  if(dim!=-1) mesh->ntriangles=fileinfo.dimension[dim].length;

  dim=cdf_identify_dimension(fileinfo,"E");
  if(dim!=-1) mesh->nedges=fileinfo.dimension[dim].length;

  dim=cdf_identify_dimension(fileinfo,"siglay");
  if(dim!=-1) mesh->nlayers=fileinfo.dimension[dim].length;

  dim=cdf_identify_dimension(fileinfo,"siglev");
  if(dim!=-1) mesh->nlevels=fileinfo.dimension[dim].length;
  
  if(mesh->nvtxs<=0)NC_TRAP_ERROR(return,NC_EBADDIM,verbose,"Could not find number of vertices.");

  mesh->vertices = new vertex_t[mesh->nvtxs];
  mesh->vertices[0].elmts=NULL;
  // if(mesh->nedges!=0)       mesh->edges       = new edge_t[mesh->nedges];
  if(mesh->ntriangles!=0)   mesh->triangles   = new triangle_t[mesh->ntriangles];
  if(mesh->nquadrangles!=0) mesh->quadrangles = new quadrangle_t[mesh->nquadrangles];

  x=new double[mesh->nvtxs];
  y=new double[mesh->nvtxs];
  h=new double[mesh->nvtxs];
  z=0;
  if(mesh->nlevels!=0)  z=new double[mesh->nvtxs*mesh->nlevels];

/* *----------------------------------------------------------------------------
  T-UGOM geometry*/
  status=nc_inq_varid (ncid,"lon",&varid);
  if(status!=0) {
/* *----------------------------------------------------------------------------
    WW3*/
    status=nc_inq_varid (ncid,"longitude",&varid);
    }
  status=nc_get_var_double(ncid,varid,x);

  status=nc_inq_varid (ncid,"lat",&varid);
  if(status!=0) {
/* *----------------------------------------------------------------------------
    WW3*/
    status=nc_inq_varid (ncid,"latitude",&varid);
    }
  status=nc_get_var_double(ncid,varid,y);

  mesh->type=SPHERICAL;

/* *----------------------------------------------------------------------------
  T-UGOM bathymetry*/
  status=nc_inq_varid     (ncid,"bathymetry",&varid);
  if(status==NC_NOERR) {
    status=nc_get_var_double(ncid,varid,h);
    }

/* *----------------------------------------------------------------------------
  FVCOM bathymetry*/
  status=nc_inq_varid     (ncid,"h",&varid);
  if(status==NC_NOERR) {
    status=nc_get_var_double(ncid,varid,h);
    }

  for (n=0;n<mesh->nvtxs;n++) {
    mesh->vertices[n].lon=x[n];
    mesh->vertices[n].lat=y[n];
    mesh->vertices[n].h=h[n];
    }

  if(x!=0) delete[] x;
  if(y!=0) delete[] y;
  if(z!=0) delete[] z;
  
  delete[] h;

/* *----------------------------------------------------------------------------
  T-UGOM, former format*/
  status=nc_inq_varid     (ncid,"element",&varid);
  if(status==NC_NOERR) {
    v=new int[mesh->ntriangles*3];
    status=nc_get_var_int   (ncid,varid,v);
    fe_set_triangles(mesh,v);
    delete[] v;
    }

/* *----------------------------------------------------------------------------
  T-UGOM, new format*/
  status=nc_inq_varid     (ncid,"triangles",&varid);
  if(status==NC_NOERR) {
    v=new int[mesh->ntriangles*3];
    status=nc_get_var_int   (ncid,varid,v);
    fe_set_triangles(mesh,v);
    delete[] v;
    }

  status=nc_inq_varid     (ncid,"quadrangles",&varid);
  if(status==NC_NOERR) {
    v=new int[mesh->nquadrangles*4];
    status=nc_get_var_int   (ncid,varid,v);
    fe_set_quadrangles(mesh,v);
    delete[] v;
    }

/* *----------------------------------------------------------------------------
  WW3*/
  status=nc_inq_varid     (ncid,"tri",&varid);
  if(status==NC_NOERR) {
    v=new int[mesh->ntriangles*3];
    status=nc_get_var_int   (ncid,varid,v);
    fe_set_triangles(mesh,v,1);
    delete[] v;
    }

/* *----------------------------------------------------------------------------
  FVCOM*/
  status=nc_inq_varid     (ncid,"nv",&varid);
  if(status==NC_NOERR) {
    v=new int[mesh->ntriangles*3];
    status=nc_get_var_int   (ncid,varid,v);
    for (n=0;n<mesh->ntriangles;n++) {
      for(k=0;k<3;k++) {
        mesh->triangles[n].vertex[k]=v[n+mesh->ntriangles*k]-1;
        }
      }
    delete[] v;
    }

  status = nc_inq_varid(ncid, "edges", &varid);
  if(status==NC_NOERR) {
    v=new int[mesh->ntriangles*3+mesh->nquadrangles*4];
    status = nc_get_var_int(ncid, varid, v);
    for(n = 0; n < mesh->ntriangles; n++) {
      for(k = 0; k < 3; k++) {
        mesh->triangles[n].edges[k] = v[n * 3 + k];
        }
      }
    for(n = 0; n < mesh->nquadrangles; n++) {
      for(k = 0; k < 4; k++) {
        mesh->quadrangles[n].edges[k] = v[mesh->ntriangles*3 + n * 4 + k];
        }
      }
    delete[] v;
    }
  else {
    if(verbose) STDOUT_BASE_LINE("warning: edge table not found in %s\n",filename);
    mesh->nedges=0;
    deletep(&mesh->edges);
    }
  
  fe_init_from_elements(mesh);

  fileinfo.destroy();

  status = nc_close(ncid);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void fe_init_from_elements(mesh_t *mesh,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///init mesh from its triangles
/**
\param *mesh
*/
/*----------------------------------------------------------------------------*/
{
  int status,n;
  bool debug=false;
  struct timeval b4;
  gettimeofday(&b4);
  
/*------------------------------------------------------------------------------
  elements to neighbours */
  if(verbose>0) STDERR_BASE_LINE(" fe_e2n...\n");
  status=fe_e2n (mesh,max(0,verbose-1));
  
  if(mesh->nedges==0 || mesh->edges==0) {
/*------------------------------------------------------------------------------
    edge table need to be built from scratch */
    if(verbose>0) STDERR_BASE_LINE("%gs: init_edge_table...\n",difftime(&b4));
//     status=init_edge_table(mesh,max(0,verbose-1));
    status=fe_edgetable(mesh,1,verbose,debug);
    }
  
  if(verbose>0) STDERR_BASE_LINE("%gs: fe_e2edges...\n",difftime(&b4));
  status=fe_e2edges (mesh);
  
/*------------------------------------------------------------------------------
  initialise affine constants */
  if(verbose>0) STDERR_BASE_LINE("%gs: fe_initaffine...\n",difftime(&b4));
  #pragma omp parallel for private(n,status)
  for(n=0; n<mesh->ntriangles; n++) {
    status=fe_initaffine(mesh,n);
    if(mesh->triangles[n].Area <=0.0)
      #pragma omp critical(consoleOutput)
      {
      printf("element %d has negative area (cw order)\n",n);
      }
    }
  if(verbose>0) STDERR_BASE_LINE("%gs: fe_initaffine_spherical(,quadrangle_t &,)...\n",difftime(&b4));
  #pragma omp parallel for private(n,status)
  for(n=0; n<mesh->nquadrangles; n++) {
    status=fe_initaffine_spherical(*mesh,&mesh->quadrangles[n],n);
    if(mesh->quadrangles[n].Area <=0.0)
      #pragma omp critical(consoleOutput)
      {
      printf("element %d has negative area (cw order)\n",n);
      }
    }
  if(verbose>0) STDERR_BASE_LINE("%gs.\n",difftime(b4));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_readmesh3d(const char *filename,mesh_t *mesh, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, k, status, ncid;
  size_t *dimlgth;
  int ndimsp,nvarsp,ngattsp,unlimdimidp,dim,*nattsp,**dimids,*ndim;
  nc_type *vartype;
  char **dimname,**varname;
  int var,varid;
  double *x=0,*y=0,*z=0,*h=0,H;
  float *sigma=0,*elevation=0,mask;
  int *v;
  variable_t details;
/*----------------------------------------------------------------------------*/
/* read time dependent fields in netcdf file */


  status=nc_open(filename,0,&ncid);
  if(status != NC_NOERR) goto error;

  status=nc_inq(ncid,&ndimsp,&nvarsp,&ngattsp,&unlimdimidp);
  if(status != NC_NOERR) goto error;
  printf("ncid %d ,ndimsp %d,nvarsp %d,ngattsp %d,unlimdimidp %d \n",
          ncid,ndimsp,nvarsp,ngattsp,unlimdimidp);

  varname = new char *[nvarsp];
  nattsp  = new int[nvarsp];
  vartype = new nc_type[nvarsp];
  dimids  = new int *[nvarsp];
  ndim    = new int[nvarsp];
  dimname = new char *[ndimsp];
  dimlgth = new size_t[ndimsp];

/*------------------------------------------------------------------------------
  inquire dimensions name and length*/
  for(dim = 0; dim < ndimsp; dim++) {
    dimname[dim] = new char[1024];
    status = nc_inq_dim(ncid, dim, dimname[dim], &dimlgth[dim]);
    if(status != NC_NOERR) goto error;
    if (verbose) printf("dimension %d,name %s, length %d \n",dim,dimname[dim],dimlgth[dim]);
    }

/*------------------------------------------------------------------------------
  inquire variables*/
  for (var=0;var<nvarsp;var++) {
    varname[var] = new char[1024];
    status = nc_inq_varndims(ncid, var, &ndim[var]);
    if(status != NC_NOERR)
      goto error;
    dimids[var] = new int[ndim[var]];
    status = nc_inq_var(ncid, var, varname[var], &vartype[var], &ndim[var],
                     dimids[var], &nattsp[var]);
    }

  mesh->destroy();

/* *----------------------------------------------------------------------------
  T-UGOM*/
  for(dim=0;dim<ndimsp;dim++) {
    if(strcmp(dimname[dim],"N")==0) {
      mesh->nvtxs=dimlgth[dim];
      continue;
      }
    if(strcmp(dimname[dim],"M")==0) {
      mesh->ntriangles=dimlgth[dim];
      continue;
      }
    if(strcmp(dimname[dim],"NQUADRANGLES")==0) {
      mesh->nquadrangles=dimlgth[dim];
      continue;
      }
    if(strcmp(dimname[dim],"E")==0) {
      mesh->nedges=dimlgth[dim];
      continue;
      }
    if(strcmp(dimname[dim],"K")==0) {
      mesh->nlayers=dimlgth[dim];
      continue;
      }
    if(strcmp(dimname[dim],"L")==0) {
      mesh->nlevels=dimlgth[dim];
      continue;
      }
    }
/* *----------------------------------------------------------------------------
  WW3*/
  for(dim=0;dim<ndimsp;dim++) {
    if(strcmp(dimname[dim],"node")==0) {
      mesh->nvtxs=dimlgth[dim];
      continue;
      }
    if(strcmp(dimname[dim],"level")==0) {
      mesh->nlayers=dimlgth[dim];
      continue;
      }
    if(strcmp(dimname[dim],"element")==0) {
      mesh->ntriangles=dimlgth[dim];
      continue;
      }
    }

/* *----------------------------------------------------------------------------
  FVCOM*/
  for(dim=0;dim<ndimsp;dim++) {
    if(strcmp(dimname[dim],"node")==0) {
      mesh->nvtxs=dimlgth[dim];
      continue;
      }
    if(strcmp(dimname[dim],"nele")==0) {
      mesh->ntriangles=dimlgth[dim];
      continue;
      }
    if(strcmp(dimname[dim],"E")==0) {
      mesh->nedges=dimlgth[dim];
      continue;
      }
    if(strcmp(dimname[dim],"siglay")==0) {
      mesh->nlayers=dimlgth[dim];
      continue;
      }
    if(strcmp(dimname[dim],"siglev")==0) {
      mesh->nlevels=dimlgth[dim];
      continue;
      }
    }

  mesh->vertices = new vertex_t[mesh->nvtxs];
  if(mesh->nedges!=0)       mesh->edges       = new edge_t[mesh->nedges];
  if(mesh->ntriangles!=0)   mesh->triangles   = new triangle_t[mesh->ntriangles];
  if(mesh->nquadrangles!=0) mesh->quadrangles = new quadrangle_t[mesh->nquadrangles];

  x = new double[mesh->nvtxs];
  y = new double[mesh->nvtxs];
  h = new double[mesh->nvtxs];
  z=0;
  if(mesh->nlevels!=0)  z=new double[mesh->nvtxs*mesh->nlevels];
/* *----------------------------------------------------------------------------
  T-UGOM geometry*/
  status=nc_inq_varid (ncid,"lon",&varid);
  if(status!=0) {
/* *----------------------------------------------------------------------------
    WW3*/
    status=nc_inq_varid (ncid,"longitude",&varid);
    }
  status=nc_get_var_double(ncid,varid,x);

  status=nc_inq_varid (ncid,"lat",&varid);
  if(status!=0) {
/* *----------------------------------------------------------------------------
    WW3*/
    status=nc_inq_varid (ncid,"latitude",&varid);
    }
  status=nc_get_var_double(ncid,varid,y);

//  mesh->type=SPHERICAL;
    mesh->type=0;

/* *----------------------------------------------------------------------------
  T-UGOM bathymetry*/
  status=nc_inq_varid     (ncid,"bathymetry",&varid);
  if(status==NC_NOERR) {
    status=nc_get_var_double(ncid,varid,h);
    }

/* *----------------------------------------------------------------------------
  FVCOM bathymetry*/
  status=nc_inq_varid     (ncid,"h",&varid);
  if(status==NC_NOERR) {
    status=nc_get_var_double(ncid,varid,h);
//     for (n=0;n<mesh->nvtxs;n++) {
//       x[n]/=10000.;
//       y[n]/=10000.;
//       }
    }

/* *----------------------------------------------------------------------------
  T-UGOM levels*/
//  status=nc_inq_varid     (ncid,"depth",&varid);
  status=nc_inq_varid     (ncid,"z-LGP1",&varid);
  if(status!=NC_NOERR) status=nc_inq_varid(ncid,"depth",&varid);
  if(status==NC_NOERR) {
    status=nc_get_var_double(ncid,varid,z);
    }
  status=nc_inq_varid     (ncid,"elevation",&varid);
   if(status==NC_NOERR) {
     status=fe_read2d(ncid, varid, 0, 0, &elevation, &mask);
    }
  else {
    elevation=new float[mesh->nvtxs];
    for (n=0;n<mesh->nvtxs;n++) elevation[n]=0;
    }
  status=nc_inq_varid     (ncid,"H",&varid);
  if(status==NC_NOERR) {
    status=fe_read2d(ncid, varid, 0, 0, &elevation, &mask);
    for (n=0;n<mesh->nvtxs;n++) {
      elevation[n]-=h[n];
      }
    }
/* *----------------------------------------------------------------------------
  FVCOM levels*/
  status=nc_inq_varid     (ncid,"siglev",&varid);
  if(status==NC_NOERR) {
    sigma=new float[mesh->nlevels];
    status=nc_get_var_float(ncid,varid,sigma);
    for (n=0;n<mesh->nvtxs;n++) {
      for(k=0;k<mesh->nlevels;k++) {
        z[n*mesh->nlevels+k]=h[n]*sigma[k];
        }
      }
    }
//   status=nc_inq_varid     (ncid,"zeta",&varid);
//   if(status==NC_NOERR) {
//     status=fe_read2d(ncid, varid, 0, frame, &elevation, &mask,&details);
//     }

  if(elevation==0) {
    elevation=new float[mesh->nvtxs];
    for (n=0;n<mesh->nvtxs;n++) elevation[n]=1.;
    }
  for (n=0;n<mesh->nvtxs;n++) {
    mesh->vertices[n].lon=x[n];
    mesh->vertices[n].lat=y[n];
    mesh->vertices[n].h=h[n];
    mesh->vertices[n].sigma  =(double *) malloc(mesh->nlevels*sizeof(double));
    mesh->vertices[n].zlevels=(double *) malloc(mesh->nlevels*sizeof(double));
    if(mesh->nlevels>2) {
      for(k=0;k<mesh->nlevels;k++) {
//        mesh->vertices[n].zlevels[k]=z[n*mesh->nlevels+k];
        mesh->vertices[n].zlevels[k]=z[k*mesh->nvtxs+n];
        if(isnan(mesh->vertices[n].zlevels[k])==1) {
          printf("nan value at node= %d, level= %d\n",n,k);
          }
        }
      if(sigma==0) {
        H=-mesh->vertices[n].zlevels[mesh->nlevels-1]+mesh->vertices[n].zlevels[0];
/* *----------------------------------------------------------------------------
        http://www.ipp.mpg.de/~rfs/comas/Helsinki/helsinki04/CompScience/csep/csep1.phy.ornl.gov/om/node30.html
        s=(z-elevation)/(h+elevation)*/
        for(k=0;k<mesh->nlevels;k++) {
//          mesh->vertices[n].sigma[k]=(mesh->vertices[n].zlevels[k]-mesh->vertices[n].zlevels[0])/H;
          mesh->vertices[n].sigma[k]=(mesh->vertices[n].zlevels[k]-elevation[n])/(H+elevation[n]);
          if(isnan(mesh->vertices[n].sigma[k])==1) {
            printf("nan value at node= %d, level= %d\n",n,k);
            }
          }
        }
      else {
        for(k=0;k<mesh->nlevels;k++) {
          mesh->vertices[n].sigma[k]=sigma[k];
          }
        }
/* *----------------------------------------------------------------------------
        http://www.ipp.mpg.de/~rfs/comas/Helsinki/helsinki04/CompScience/csep/csep1.phy.ornl.gov/om/node30.html
        s=(z-elevation)/(h+elevation)*/
//      H=-mesh->vertices[n].zlevels[mesh->nlevels-1]+elevation[n];
      H=-mesh->vertices[n].zlevels[mesh->nlevels-1]+mesh->vertices[n].zlevels[0];
      for(k=0;k<mesh->nlevels;k++) {
        mesh->vertices[n].zlevels[k]=mesh->vertices[n].sigma[k]*(H+elevation[n])+elevation[n];
        if(isnan(mesh->vertices[n].zlevels[k])==1) {
          printf("nan value at node= %d, level= %d\n",n,k);
          }
        }
      }
    }

//   free(x);
//   free(y);
  delete[] x;
  delete[] y;
  if(z!=0) delete[] z;
//  free(h);
  delete[] h;

/* *----------------------------------------------------------------------------
  T-UGOM, former format*/
  status=nc_inq_varid     (ncid,"element",&varid);
  if(status==NC_NOERR) {
    v=new int[mesh->ntriangles*3];
    status=nc_get_var_int   (ncid,varid,v);
    fe_set_triangles(mesh,v);
    delete[] v;
    }

/* *----------------------------------------------------------------------------
  T-UGOM, new format*/
  status=nc_inq_varid     (ncid,"triangles",&varid);
  if(status==NC_NOERR) {
    v=new int[mesh->ntriangles*3];
    status=nc_get_var_int   (ncid,varid,v);
    fe_set_triangles(mesh,v);
    delete[] v;
    }

  status=nc_inq_varid     (ncid,"quadrangles",&varid);
  if(status==NC_NOERR) {
    v=new int[mesh->nquadrangles*4];
    status=nc_get_var_int   (ncid,varid,v);
    fe_set_quadrangles(mesh,v);
    delete[] v;
    }

/* *----------------------------------------------------------------------------
  WW3*/
  status=nc_inq_varid     (ncid,"tri",&varid);
  if(status==NC_NOERR) {
    v=new int[mesh->ntriangles*3];
    status=nc_get_var_int   (ncid,varid,v);
    fe_set_triangles(mesh,v,1);
    delete[] v;
    }

/* *----------------------------------------------------------------------------
  FVCOM*/
  status=nc_inq_varid     (ncid,"nv",&varid);
  if(status==NC_NOERR) {
    v=new int[mesh->ntriangles*3];
    status=nc_get_var_int   (ncid,varid,v);
    for (n=0;n<mesh->ntriangles;n++) {
      for(k=0;k<3;k++) {
        mesh->triangles[n].vertex[k]=v[n+mesh->ntriangles*k]-1;
        }
      }
    delete[] v;
    }

  status = nc_inq_varid(ncid, "edges", &varid);
  if(status==NC_NOERR) {
    v=new int[mesh->ntriangles*3+mesh->nquadrangles*4];
    status = nc_get_var_int(ncid, varid, v);
    for(n = 0; n < mesh->ntriangles; n++) {
      for(k = 0; k < 3; k++) {
        mesh->triangles[n].edges[k] = v[n * 3 + k];
        }
      }
    for(n = 0; n < mesh->nquadrangles; n++) {
      for(k = 0; k < 4; k++) {
        mesh->quadrangles[n].edges[k] = v[mesh->ntriangles*3 + n * 4 + k];
        }
      }
    delete[] v;
    }

/* *----------------------------------------------------------------------------
  T-UGOM, LGP2 nodes discretisation*/
  status=nc_inq_varid (ncid,"C-LGP2",&varid);
  if(status==NC_NOERR) {
    v=(int *) malloc(mesh->ntriangles*6*sizeof(int));
    status=nc_get_var_int   (ncid,varid,v);
    mesh->LGP2descriptor.NIbE=new int*[mesh->ntriangles];
    for (n=0;n<mesh->ntriangles;n++) {
      mesh->LGP2descriptor.NIbE[n]=new int[6];
      for(k=0;k<6;k++) {
        mesh->LGP2descriptor.NIbE[n][k]=v[n*6+k];
        }
      }
    free(v);
    }

  status=fe_e2n (mesh);

  if(mesh->nedges==0) {
/* *----------------------------------------------------------------------------
    edge table need to be built from scratch*/
    status=fe_edgetable(mesh,0,0);
    //status=build_edgetable(mesh);
    printf("mesh->nedges==0\n");
    }
  else {
/* *----------------------------------------------------------------------------
    use file's edge table*/
    //status=fe_e2edges (mesh);
    }

  status = nc_close(ncid);
  return(status);

error:
  return(status);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_savemeshNC3D(const char *filename, mesh_t & mesh, int option)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *------------------------------------------------------------------------
  attribute name "axis" changed in "content" to comply with CF standard*/
{
  int k, m, n;

  int ncid, stat;       /* netCDF id */

  /* dimension ids */
  int M_dim;
  int N_dim;
  int E_dim;
  int P_dim;
  int K_dim;
  int L_dim;
  int T_dim;

  /* variable ids */
  int time_id;
  int lon_id;
  int lat_id;
  int bcode_id;
  int edges_id;
  int element_id;
  int bathymetry_id;
  /* dimension lengths */
  size_t M_len = mesh.ntriangles;
  size_t N_len = mesh.nvtxs;
  size_t E_len = mesh.nedges;
  size_t P_len = 3;
  size_t K_len = mesh.nlayers;
  size_t L_len = mesh.nlevels;
  size_t T_len = NC_UNLIMITED;

  /* rank (number of dimensions) for each variable */
#define RANK_time 1
#define RANK_lon 1
#define RANK_lat 1
#define RANK_element 2
#define RANK_edges 2
#define RANK_bathymetry 1
#define RANK_sigma 2
#define RANK_depth 2

  /* variable shapes */
  int time_dims[RANK_time];
  int lon_dims[RANK_lon];
  int lat_dims[RANK_lat];
  int element_dims[RANK_element];
  int edges_dims[RANK_element];
  int bathymetry_dims[RANK_bathymetry];

  /* attribute vectors */
  double lon_valid_min[1];
  double lon_valid_max[1];
  double lat_valid_min[1];
  double lat_valid_max[1];
  double bathymetry_scale_factor[1];
  double bathymetry_add_offset[1];

  double *buffer;
  int *nv;
  char text[1024];

/*------------------------------------------------------------------------------
   enter define mode */
  stat = nc_create(filename, NC_CLOBBER, &ncid);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------------
   define dimensions */
  stat = nc_def_dim(ncid, "M", M_len, &M_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "N", N_len, &N_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "E", E_len, &E_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "P", P_len, &P_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "K", K_len, &K_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "L", L_len, &L_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "T", T_len, &T_dim);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------------
   define time variable */

  time_dims[0] = T_dim;
  stat = nc_def_var(ncid, "time", NC_DOUBLE, RANK_time, time_dims, &time_id);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------------
   /* assign minimum attributes */
  stat = nc_put_att_text(ncid, time_id, "units", 7, "seconds");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, time_id, "calendar", 9, "gregorian");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, time_id, "long_name", 41, "Time elasped in seconds since time_origin");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, time_id, "title", 4, "Time");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, time_id, "time_origin", 20, "1950-JAN-01 00:00:00");
  nc_check_error(stat, __LINE__, __FILE__);

  bathymetry_dims[0] = N_dim;
  stat = nc_def_var(ncid, "bathymetry", NC_FLOAT, RANK_bathymetry, bathymetry_dims, &bathymetry_id);
  nc_check_error(stat, __LINE__, __FILE__);

  strcpy(text,"model_positive_bathymetry_node");
  stat = nc_put_att_text(ncid, bathymetry_id, "long_name", strlen(text),text);
  nc_check_error(stat, __LINE__, __FILE__);
  
  strcpy(text,"bathymetry");
  stat = nc_put_att_text(ncid, bathymetry_id, "short_name", strlen(text),text);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "units", 1, "m");
  nc_check_error(stat, __LINE__, __FILE__);
  bathymetry_scale_factor[0] = 1;
  stat = nc_put_att_double(ncid, bathymetry_id, "scale_factor", NC_DOUBLE, 1,bathymetry_scale_factor);
  nc_check_error(stat, __LINE__, __FILE__);
  bathymetry_add_offset[0] = 0;
  stat = nc_put_att_double(ncid, bathymetry_id, "add_offset", NC_DOUBLE, 1,bathymetry_add_offset);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "associate", 7, "lon lat");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "subgrid", 5, "point");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "content", 1, "N");
  nc_check_error(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, "CF 1.0");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, NC_GLOBAL, "grid_type", 2, "UG");
  nc_check_error(stat, __LINE__, __FILE__);
  
  strcpy(text,"T-UGOm archive");
  stat = nc_put_att_text(ncid, NC_GLOBAL, "Topic", strlen(text),text);
  nc_check_error(stat, __LINE__, __FILE__);
  
  stat = nc_put_att_text(ncid, NC_GLOBAL, "file_name", strlen(filename), filename);
  nc_check_error(stat, __LINE__, __FILE__);
  
  strcpy(text,"F. Lyard/POC 2007");
  stat = nc_put_att_text(ncid, NC_GLOBAL, "production", strlen(text),text);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, NC_GLOBAL, "history", 1, " ");
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------------
   leave define mode */
  if(option == 0) {
    stat = nc_enddef(ncid);
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_close(ncid);
    nc_check_error(stat, __LINE__, __FILE__);
    return (0);
    }

/*------------------------------------------------------------------------------
   define other variables */

  lon_dims[0] = N_dim;
  stat = nc_def_var(ncid, "lon", NC_DOUBLE, RANK_lon, lon_dims, &lon_id);
  nc_check_error(stat, __LINE__, __FILE__);

  lat_dims[0] = N_dim;
  stat = nc_def_var(ncid, "lat", NC_DOUBLE, RANK_lat, lat_dims, &lat_id);
  nc_check_error(stat, __LINE__, __FILE__);

  element_dims[0] = M_dim;
  element_dims[1] = P_dim;
  stat = nc_def_var(ncid, "element", NC_INT, RANK_element, element_dims,&element_id);
  nc_check_error(stat, __LINE__, __FILE__);

  edges_dims[0] = M_dim;
  edges_dims[1] = P_dim;
  stat = nc_def_var(ncid, "edges", NC_INT, RANK_edges, edges_dims,&edges_id);
  nc_check_error(stat, __LINE__, __FILE__);

//   sigma_dims[0] = N_dim;
//   sigma_dims[1] = L_dim;
//   stat = nc_def_var(ncid, "sigma", NC_FLOAT, RANK_sigma, sigma_dims,&sigma_id);
//   nc_check_error(stat, __LINE__, __FILE__);
//
//   depth_dims[0] = N_dim;
//   depth_dims[1] = L_dim;
//   stat = nc_def_var(ncid, "depth", NC_FLOAT, RANK_depth, depth_dims, &depth_id);
//   nc_check_error(stat, __LINE__, __FILE__);

  /* assign attributes */
  stat = nc_put_att_text(ncid, lon_id, "units", 12, "degrees_east");
  nc_check_error(stat, __LINE__, __FILE__);
  lon_valid_min[0] = -180;
  stat = nc_put_att_double(ncid, lon_id, "valid_min", NC_DOUBLE, 1, lon_valid_min);
  nc_check_error(stat, __LINE__, __FILE__);

  lon_valid_max[0] = 180;
  stat = nc_put_att_double(ncid, lon_id, "valid_max", NC_DOUBLE, 1, lon_valid_max);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lon_id, "long_name", 9, "longitude");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lon_id, "standard_name", 9, "longitude");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lon_id, "subgrid", 5, "point");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lon_id, "content", 1, "N");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "units", 13, "degrees_north");
  nc_check_error(stat, __LINE__, __FILE__);

  lat_valid_min[0] = -90;
  stat = nc_put_att_double(ncid, lat_id, "valid_min", NC_DOUBLE, 1,lat_valid_min);
  nc_check_error(stat, __LINE__, __FILE__);
  lat_valid_max[0] = 90;
  stat = nc_put_att_double(ncid, lat_id, "valid_max", NC_DOUBLE, 1,lat_valid_max);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "long_name", 8, "latitude");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "standard_name", 8, "latitude");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "subgrid", 5, "point");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "content", 1, "N");
  nc_check_error(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, element_id, "long_name", 20, "element_connectivity");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, element_id, "standard_name", 7, "element");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, element_id, "subgrid", 4, "cell");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, element_id, "content", 2, "MP");
  nc_check_error(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, edges_id, "long_name", 20, "element_edges");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, edges_id, "standard_name", 7, "edges");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, edges_id, "subgrid", 4, "cell");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, edges_id, "content", 2, "MP");
  nc_check_error(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, bcode_id, "long_name", strlen("boundary_code"), "boundary_code");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bcode_id, "standard_name", strlen("boundary_code"), "boundary_code");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bcode_id, "subgrid", 4, "cell");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bcode_id, "content", 1, "N");
  nc_check_error(stat, __LINE__, __FILE__);

//   stat = nc_put_att_text(ncid, sigma_id, "long_name", 37,"generalized sigma coordinate at nodes");
//   nc_check_error(stat, __LINE__, __FILE__);
//   stat = nc_put_att_text(ncid, sigma_id, "standard_name", 35,"general_ocean_sigma_coordinate_node");
//   nc_check_error(stat, __LINE__, __FILE__);
//   stat = nc_put_att_text(ncid, sigma_id, "short_name", 5, "sigma");
//   nc_check_error(stat, __LINE__, __FILE__);
//   stat = nc_put_att_text(ncid, sigma_id, "units", 1, "1");
//   nc_check_error(stat, __LINE__, __FILE__);
//   sigma_scale_factor[0] = 1;
//   stat = nc_put_att_double(ncid, sigma_id, "scale_factor", NC_DOUBLE, 1, sigma_scale_factor);
//   nc_check_error(stat, __LINE__, __FILE__);
//   sigma_add_offset[0] = 0;
//   stat = nc_put_att_double(ncid, sigma_id, "add_offset", NC_DOUBLE, 1, sigma_add_offset);
//   nc_check_error(stat, __LINE__, __FILE__);
//   stat = nc_put_att_text(ncid, sigma_id, "associate", 7, "lon lat");
//   nc_check_error(stat, __LINE__, __FILE__);
//   stat = nc_put_att_text(ncid, sigma_id, "subgrid", 5, "point");
//   nc_check_error(stat, __LINE__, __FILE__);
//   stat = nc_put_att_text(ncid, sigma_id, "content", 2, "NZ");
//   nc_check_error(stat, __LINE__, __FILE__);
//   stat = nc_put_att_text(ncid, sigma_id, "positive", 2, "up");
//   nc_check_error(stat, __LINE__, __FILE__);

//   stat = nc_put_att_text(ncid, depth_id, "long_name", 37, "generalized depth coordinate at nodes");
//   nc_check_error(stat, __LINE__, __FILE__);
//   stat = nc_put_att_text(ncid, depth_id, "standard_name", 35, "general_ocean_depth_coordinate_node");
//   nc_check_error(stat, __LINE__, __FILE__);
//   stat = nc_put_att_text(ncid, depth_id, "short_name", 5, "depth");
//   nc_check_error(stat, __LINE__, __FILE__);
//   stat = nc_put_att_text(ncid, depth_id, "units", 1, "1");
//   nc_check_error(stat, __LINE__, __FILE__);
//   depth_scale_factor[0] = 1;
//   stat = nc_put_att_double(ncid, depth_id, "scale_factor", NC_DOUBLE, 1,depth_scale_factor);
//   nc_check_error(stat, __LINE__, __FILE__);
//   depth_add_offset[0] = 0;
//   stat = nc_put_att_double(ncid, depth_id, "add_offset", NC_DOUBLE, 1, depth_add_offset);
//   nc_check_error(stat, __LINE__, __FILE__);
//   stat = nc_put_att_text(ncid, depth_id, "associate", 7, "lon lat");
//   nc_check_error(stat, __LINE__, __FILE__);
//   stat = nc_put_att_text(ncid, depth_id, "subgrid", 5, "point");
//   nc_check_error(stat, __LINE__, __FILE__);
//   stat = nc_put_att_text(ncid, depth_id, "content", 2, "NZ");
//   nc_check_error(stat, __LINE__, __FILE__);
//   stat = nc_put_att_text(ncid, depth_id, "positive", 2, "up");
//   nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------------
   leave define mode */
  stat = nc_enddef(ncid);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_close(ncid);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------------
  write time independent fields in netcdf file */
  stat = nc_open(filename, NC_WRITE, &ncid);

  buffer = new double[N_len];

/*------------------------------------------------------------------------------
  longitude*/
  for(n = 0; n < N_len; n++)
    buffer[n] = mesh.vertices[n].lon;
  stat = nc_put_var_double(ncid, lon_id, buffer);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------------
  latitude*/
  for(n = 0; n < N_len; n++)
    buffer[n] = mesh.vertices[n].lat;
  stat = nc_put_var_double(ncid, lat_id, buffer);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------------
  incidence*/
  nv = new int[M_len * P_len];
  for(m = 0; m < M_len; m++) {
    for(k = 0; k < P_len; k++)
      nv[m * P_len + k] = mesh.triangles[m].vertex[k];
    }
  stat = nc_put_var_int(ncid, element_id, nv);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------------
  edges*/
  for(m = 0; m < M_len; m++) {
    for(k = 0; k < P_len; k++)
      nv[m * P_len + k] = mesh.triangles[m].edges[k];
    }
  stat = nc_put_var_int(ncid, edges_id, nv);
  nc_check_error(stat, __LINE__, __FILE__);
  zaparr(nv);

/*------------------------------------------------------------------------------
  boundary codes*/
  nv = new int[N_len];
  for(m = 0; m < N_len; m++) {
    nv[m] = mesh.vertices[m].code;
    }
  stat = nc_put_var_int(ncid, bcode_id, nv);
  nc_check_error(stat, __LINE__, __FILE__);
  zaparr(nv);

/*------------------------------------------------------------------------------
  bathymétrie*/
  for(n = 0; n < N_len; n++)
    buffer[n] = mesh.vertices[n].h;
  stat = nc_put_var_double(ncid, bathymetry_id, buffer);
  nc_check_error(stat, __LINE__, __FILE__);
  zaparr(buffer);

  buffer = new double[N_len*L_len];
/*------------------------------------------------------------------------------
  initial sigma levels*/
//   for(n = 0; n < N_len; n++)
//     for(l = 0; l < L_len; l++)
//       buffer[n*L_len+l] = mesh.vertices[n].sigma[l];
//   stat = nc_put_var_double(ncid, sigma_id, buffer);
//   nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------------
  initial depth levels*/
//   for(n = 0; n < N_len; n++) {
//     if(mesh.vertices[n].zlevels==0) {
//       for(l = 0; l < L_len; l++)
//         buffer[n*L_len+l] = -mesh.vertices[n].sigma[l]*mesh.vertices[n].h;
//       }
//     else {
//       for(l = 0; l < L_len; l++)
//         buffer[n*L_len+l] = mesh.vertices[n].zlevels[l];
//       }
//     }
//   stat = nc_put_var_double(ncid, depth_id, buffer);
//   nc_check_error(stat, __LINE__, __FILE__);
  zaparr(buffer);

  stat = nc_close(ncid);
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_updatemeshNC3D(const char *filename, mesh_t mesh, int option)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n;

  int ncid, stat;       /* netCDF id */

  /* variable ids */
  int lon_id;
  int lat_id;
  /* dimension lengths */
  size_t N_len = mesh.nvtxs;

  /* rank (number of dimensions) for each variable */

  double *buffer;
  cdfvar_t info;
  int status;

/*------------------------------------------------------------------------------
  write time independent fields in netcdf file */
  stat = nc_open(filename, NC_WRITE, &ncid);

  buffer = new double[N_len];

/*------------------------------------------------------------------------------
  longitude*/
  status = nc_inq_varid(ncid, "lon", &lon_id);
//  status= cdf_varinfo( ncid, lon_id, &info);
  for(n = 0; n < N_len; n++)
    buffer[n] = mesh.vertices[n].lon;
  stat = nc_put_var_double(ncid, lon_id, buffer);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------------
  latitude*/
  status = nc_inq_varid(ncid, "lat", &lat_id);
  for(n = 0; n < N_len; n++)
    buffer[n] = mesh.vertices[n].lat;
  stat = nc_put_var_double(ncid, lat_id, buffer);
  nc_check_error(stat, __LINE__, __FILE__);


  stat = nc_close(ncid);
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_savemeshNC2D(const char *filename, mesh_t & mesh, int option)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, m, n;

  int ncid, stat;       /* netCDF id */

  /* dimension ids */
  int M_dim;
  int N_dim;
  int E_dim;
  int P_dim;
  int T_dim;

  /* variable ids */
  int lon_id;
  int lat_id;
  int element_id;
  int bathymetry_id;
  /* dimension lengths */
  size_t M_len = mesh.ntriangles;
  size_t N_len = mesh.nvtxs;
  size_t E_len = mesh.nedges;
  size_t P_len = 3;
  size_t T_len = NC_UNLIMITED;

  /* rank (number of dimensions) for each variable */
#define RANK_lon 1
#define RANK_lat 1
#define RANK_element 2
#define RANK_bathymetry 1

  /* variable shapes */
  int lon_dims[RANK_lon];
  int lat_dims[RANK_lat];
  int element_dims[RANK_element];
  int bathymetry_dims[RANK_bathymetry];

  /* attribute vectors */
  double lon_valid_min[1];
  double lon_valid_max[1];
  double lat_valid_min[1];
  double lat_valid_max[1];
  double bathymetry_scale_factor[1];
  double bathymetry_add_offset[1];

  double *buffer;
  int *nv;
  char text[1024];

/*------------------------------------------------------------------------------
   enter define mode */
  stat = nc_create(filename, NC_CLOBBER, &ncid);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------------
   define dimensions */
  stat = nc_def_dim(ncid, "M", M_len, &M_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "N", N_len, &N_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "E", E_len, &E_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "P", P_len, &P_dim);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "T", T_len, &T_dim);
  nc_check_error(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, "CF 1.0");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, NC_GLOBAL, "grid_type", 2, "UG");
  nc_check_error(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, NC_GLOBAL, "file_name", strlen(filename), filename);
  nc_check_error(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, NC_GLOBAL, "history", 1, " ");
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------------
   leave define mode */
  if(option == 0) {
    stat = nc_enddef(ncid);
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_close(ncid);
    nc_check_error(stat, __LINE__, __FILE__);
    return (0);
    }

/*------------------------------------------------------------------------------
   define other variables */

  lon_dims[0] = N_dim;
  stat = nc_def_var(ncid, "lon", NC_DOUBLE, RANK_lon, lon_dims, &lon_id);
  nc_check_error(stat, __LINE__, __FILE__);

  lat_dims[0] = N_dim;
  stat = nc_def_var(ncid, "lat", NC_DOUBLE, RANK_lat, lat_dims, &lat_id);
  nc_check_error(stat, __LINE__, __FILE__);

  element_dims[0] = M_dim;
  element_dims[1] = P_dim;
  stat = nc_def_var(ncid, "element", NC_INT, RANK_element, element_dims,&element_id);
  nc_check_error(stat, __LINE__, __FILE__);

  /* assign attributes */
  stat = nc_put_att_text(ncid, lon_id, "units", 12, "degrees_east");
  nc_check_error(stat, __LINE__, __FILE__);
  lon_valid_min[0] = -180;
  stat = nc_put_att_double(ncid, lon_id, "valid_min", NC_DOUBLE, 1,lon_valid_min);
  nc_check_error(stat, __LINE__, __FILE__);

  lon_valid_max[0] = 180;
  stat = nc_put_att_double(ncid, lon_id, "valid_max", NC_DOUBLE, 1, lon_valid_max);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lon_id, "long_name", 9, "longitude");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lon_id, "standard_name", 9, "longitude");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lon_id, "subgrid", 5, "point");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lon_id, "content", 1, "N");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "units", 13, "degrees_north");
  nc_check_error(stat, __LINE__, __FILE__);

  lat_valid_min[0] = -90;
  stat = nc_put_att_double(ncid, lat_id, "valid_min", NC_DOUBLE, 1,lat_valid_min);
  nc_check_error(stat, __LINE__, __FILE__);
  lat_valid_max[0] = 90;
  stat = nc_put_att_double(ncid, lat_id, "valid_max", NC_DOUBLE, 1,lat_valid_max);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "long_name", 8, "latitude");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "standard_name", 8, "latitude");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "subgrid", 5, "point");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, lat_id, "content", 1, "N");
  nc_check_error(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, element_id, "long_name", 20,"element_connectivity");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, element_id, "standard_name", 7, "element");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, element_id, "subgrid", 4, "cell");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, element_id, "content", 2, "MP");
  nc_check_error(stat, __LINE__, __FILE__);

  bathymetry_dims[0] = N_dim;
  stat = nc_def_var(ncid, "bathymetry", NC_FLOAT, RANK_bathymetry, bathymetry_dims, &bathymetry_id);
  nc_check_error(stat, __LINE__, __FILE__);

  strcpy(text,"model_positive_bathymetry_node");
  stat = nc_put_att_text(ncid, bathymetry_id, "long_name", strlen(text),text);
  nc_check_error(stat, __LINE__, __FILE__);

  strcpy(text,"bathymetry");
  stat = nc_put_att_text(ncid, bathymetry_id, "short_name", strlen(text),text);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "units", 1, "m");
  nc_check_error(stat, __LINE__, __FILE__);
  bathymetry_scale_factor[0] = 1;
  stat = nc_put_att_double(ncid, bathymetry_id, "scale_factor", NC_DOUBLE, 1,bathymetry_scale_factor);
  nc_check_error(stat, __LINE__, __FILE__);
  bathymetry_add_offset[0] = 0;
  stat = nc_put_att_double(ncid, bathymetry_id, "add_offset", NC_DOUBLE, 1,bathymetry_add_offset);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "associate", 7, "lon lat");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "subgrid", 5, "point");
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, bathymetry_id, "content", 1, "N");
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------------
   leave define mode */
  stat = nc_enddef(ncid);
  nc_check_error(stat, __LINE__, __FILE__);
  stat = nc_close(ncid);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------------
  write time independent fields in netcdf file */
  stat = nc_open(filename, NC_WRITE, &ncid);

  buffer = new double[N_len];

/*------------------------------------------------------------------------------
  longitude*/
  for(n = 0; n < N_len; n++)
    buffer[n] = mesh.vertices[n].lon;
  stat = nc_put_var_double(ncid, lon_id, buffer);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------------
  latitude*/
  for(n = 0; n < N_len; n++)
    buffer[n] = mesh.vertices[n].lat;
  stat = nc_put_var_double(ncid, lat_id, buffer);
  nc_check_error(stat, __LINE__, __FILE__);

/*------------------------------------------------------------------------------
  incidence*/
  nv = new int[M_len * P_len];
  for(m = 0; m < M_len; m++) {
    for(k = 0; k < P_len; k++)
      nv[m * P_len + k] = mesh.triangles[m].vertex[k];
    }
  stat = nc_put_var_int(ncid, element_id, nv);
  nc_check_error(stat, __LINE__, __FILE__);
  zaparr(nv);

/*------------------------------------------------------------------------------
  bathymétrie*/
  for(n = 0; n < N_len; n++)
    buffer[n] = mesh.vertices[n].h;
  stat = nc_put_var_double(ncid, bathymetry_id, buffer);
  nc_check_error(stat, __LINE__, __FILE__);
  zaparr(buffer);

  stat = nc_close(ncid);
  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_CheckDimensions(const char *filename,const mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  char **standards;
  int k, l, m, n, status;
  int ncid, stat;

  /* dimension ids */
  int M_dim;
  int Q_dim;
  int N_dim;
  int E_dim;
  int P_dim;
  int NNPE_dim;

  int K_dim;
  int L_dim;
  int L2D_dim;
  int T_dim;
  int nodes_dim;

  int nodes_id,levels_id;

  /* dimension lengths */
  size_t M_len = mesh.ntriangles;
  size_t Q_len = mesh.nquadrangles;
  size_t N_len = mesh.nvtxs;
  size_t E_len = mesh.nedges;
  size_t P_len  = 3;

  size_t K_len = mesh.nlayers;
  size_t L_len = mesh.nlevels;
  size_t L2D_len = 2;
  size_t T_len = NC_UNLIMITED;
  
  int    Element_dim;
  size_t Element_len;
  
  bool do_connectivity=false, do_levels=false;
  
  /* rank (number of dimensions) for each variable */
#define RANK_NODES 2

  /* variable shapes */
  int nodes_dims[RANK_NODES];

  int *nv;
  double *z;
  char text[1024];
  
/**----------------------------------------------------------------------
  */
	  
  stat = nc_open(filename, NC_WRITE, &ncid);
  if(stat!=NC_NOERR) {
    nc_check_error(stat, __LINE__, __FILE__);
    printf("\n\n fe_savediscretisation file=%s cpu=%d \n\n", filename, 0);
    check_error(stat, "fe_savediscretisation, nc_open failed, return bad status", __LINE__, __FILE__, 0);
    return (-1);
    }

#define dimLenAndMeshMismatchMsg "\n*** %d %s IN %s AND %d IN mesh. ***\n*** ARCHIVE WILL BE FIXED UP CONSISTENTLY ***\n"

  if(M_len!=0) {
    stat = nc_inq_dimid(ncid, "M", &M_dim);
    nc_check_error(stat, __LINE__, __FILE__);
    check_error(stat, "fe_savediscretisation failed, stop", __LINE__, __FILE__, 1);;
    
    size_t lenInFile;
    stat = nc_inq_dimlen(ncid, M_dim, &lenInFile);
    nc_check_error(stat, __LINE__, __FILE__);
    check_error(stat, "fe_savediscretisation failed, stop", __LINE__, __FILE__, 1);;
    if(lenInFile!=M_len) {
      TRAP_ERR_RETURN(NC_EDIMSIZE, 1, dimLenAndMeshMismatchMsg, lenInFile,"TRIANGLES", filename, M_len);
      }
    }
    
  if(Q_len!=0) {
    stat = nc_inq_dimid(ncid, "Q", &Q_dim);
    nc_check_error(stat, __LINE__, __FILE__);
    check_error(stat, "fe_savediscretisation failed, stop", __LINE__, __FILE__, 1);;

    size_t lenInFile;
    stat = nc_inq_dimlen(ncid, Q_dim, &lenInFile);
    nc_check_error(stat, __LINE__, __FILE__);
    check_error(stat, "fe_savediscretisation failed, stop", __LINE__, __FILE__, 1);;
    if(lenInFile!=Q_len) {
      TRAP_ERR_RETURN(NC_EDIMSIZE, 1, dimLenAndMeshMismatchMsg, lenInFile, "QUADRANGLES", filename, Q_len);
      }
    }
    
  if(N_len!=0) {
    stat = nc_inq_dimid(ncid, "N", &N_dim);
    nc_check_error(stat, __LINE__, __FILE__);
    check_error(stat, "fe_savediscretisation failed, stop", __LINE__, __FILE__, 1);;

    size_t lenInFile;
    stat = nc_inq_dimlen(ncid, N_dim, &lenInFile);
    nc_check_error(stat, __LINE__, __FILE__);
    check_error(stat, "fe_savediscretisation failed, stop", __LINE__, __FILE__, 1);;
    if(lenInFile!=N_len) {
      TRAP_ERR_RETURN(NC_EDIMSIZE, 1, dimLenAndMeshMismatchMsg, lenInFile, "VERTICES", filename, N_len);
      }
    }
    
  stat = nc_close(ncid);
  nc_check_error(stat, __LINE__, __FILE__);
  
  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_savediscretisation(const char *filename, mesh_t & mesh, discretisation_t & discretisation, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  char **standards;
  int k, m, status;
  int ncid, stat;

  /* dimension ids */
  int M_dim;
  int Q_dim;
  int NNPE_dim;

  int L_dim;
  int nodes_dim;

  int nodes_id,levels_id;

  /* dimension lengths */
  size_t M_len = mesh.ntriangles;
  size_t Q_len = mesh.nquadrangles;
  size_t NNPE_len = discretisation.nnpe;

  size_t L_len = mesh.nlevels;

  size_t nodes_len = discretisation.nnodes;
  
  int    Element_dim;
  size_t Element_len;
  
  /* rank (number of dimensions) for each variable */
#define RANK_NODES 2

  /* variable shapes */
  int nodes_dims[RANK_NODES];

  int *nv;
  double *z;

/**----------------------------------------------------------------------
  */
  switch(discretisation.type) {
    case LGP0:
      standards=(char **)LGP0_standards;
      break;

    case LGP1:
      standards=(char **)LGP1_standards;
      break;

    case DGP1:
      standards=(char **)DGP1_standards;
      break;

    case NCP1:
      check_error(-1, "discretisation not implemented (edge-implicit)", __LINE__, __FILE__, 1);
      break;

    case DNP1:
      standards=(char **)DNP1_standards;
      break;

    case LGP2:
      standards=(char **)LGP2_standards;
      break;

    case CQP0:
      standards=(char **)CQP0_standards;
      break;

    case CQN1:
      standards=(char **)CQN1_standards;
      break;

    default:
      TRAP_ERR_EXIT(-1,"discretisation %d not yet implemented\n",discretisation.type);
      break;
    }

/**------------------------------------------------------------------------
  enter define mode */
  stat = nc_open(filename, NC_WRITE, &ncid);
  nc_check_error(stat, __LINE__, __FILE__);

  stat=nc_redef(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

  if(M_len!=0) {
    stat = nc_inq_dimid(ncid, "M", &M_dim);
    nc_check_error(stat, __LINE__, __FILE__);
    check_error(stat, "fe_savediscretisation failed, stop", __LINE__, __FILE__, 1);;

    stat = nc_inq_dimlen(ncid, M_dim, &M_len);
    nc_check_error(stat, __LINE__, __FILE__);
    check_error(stat, "fe_savediscretisation failed, stop", __LINE__, __FILE__, 1);;
    }
    
  if(Q_len!=0) {
    stat = nc_inq_dimid(ncid, "Q", &Q_dim);
    nc_check_error(stat, __LINE__, __FILE__);
    check_error(stat, "fe_savediscretisation failed, stop", __LINE__, __FILE__, 1);;

    stat = nc_inq_dimlen(ncid, Q_dim, &Q_len);
    nc_check_error(stat, __LINE__, __FILE__);
    check_error(stat, "fe_savediscretisation failed, stop", __LINE__, __FILE__, 1);;
    }
    
/**------------------------------------------------------------------------
  define dimensions */
  stat = nc_inq_dimid(ncid, standards[NNPE_DIMNAME], &NNPE_dim);
  if(stat!=0) {
    stat = nc_def_dim(ncid, standards[NNPE_DIMNAME], NNPE_len, &NNPE_dim);
    printf("discretisation %s, NNPE_len=%d\n",standards[NNPE_DIMNAME],NNPE_len);
    nc_check_error(stat, __LINE__, __FILE__);
    check_error(stat, "fe_savediscretisation failed, stop", __LINE__, __FILE__, 1);;
    }

  stat = nc_inq_dimid(ncid, standards[NNODES_DIMNAME], &nodes_dim);
  if(stat!=0) {
    stat = nc_def_dim(ncid, standards[NNODES_DIMNAME], nodes_len, &nodes_dim);
    nc_check_error(stat, __LINE__, __FILE__);
    check_error(stat, "fe_savediscretisation failed, stop", __LINE__, __FILE__, 1);;
    }

  stat = nc_inq_dimid(ncid, "L", &L_dim);
//   if(stat!=0) {
//     stat = nc_def_dim(ncid, standards[NNODES_DIMNAME], nodes_len, &nodes_dim);
//     nc_check_error(stat, __LINE__, __FILE__);
//     }

  if(mesh.ntriangles!=0) {
    Element_dim=M_dim;
    Element_len=M_len;
    }
  if(mesh.quadrangles!=0) {
    Element_dim=Q_dim;
    Element_len=Q_len;
    }
    
/**----------------------------------------------------------------------------
  define connectivity variable*/
  status = nc_inq_varid(ncid,standards[CONNECTIVITY_VARNAME], &nodes_id);
  if(status!=0) {
    nodes_dims[0] = Element_dim;
    nodes_dims[1] = NNPE_dim;
    stat = nc_def_var(ncid, standards[CONNECTIVITY_VARNAME], NC_INT, RANK_connectivity, nodes_dims, &nodes_id);
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid, nodes_id, "long_name",     strlen(standards[CONNECTIVITY_LNGNAME]), standards[CONNECTIVITY_LNGNAME]);
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid, nodes_id, "standard_name", strlen(standards[CONNECTIVITY_STDNAME]), standards[CONNECTIVITY_STDNAME]);
    nc_check_error(stat, __LINE__, __FILE__);
    }

/**----------------------------------------------------------------------------
  define levels variable*/
  status = nc_inq_varid(ncid,standards[LEVELS_VARNAME], &levels_id);
  if(status!=0) {
    nodes_dims[0] = L_dim;
    nodes_dims[1] = nodes_dim;
    stat = nc_def_var(ncid, standards[LEVELS_VARNAME], NC_DOUBLE, RANK_connectivity, nodes_dims,&levels_id);
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid, levels_id, "long_name",     strlen(standards[LEVELS_LNGNAME]), standards[LEVELS_LNGNAME]);
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid, levels_id, "standard_name", strlen(standards[LEVELS_STDNAME]), standards[LEVELS_STDNAME]);
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid, levels_id, "content", 2, "ZN");
    nc_check_error(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid, levels_id, "associates", 13, "layer lat lon");
    nc_check_error(stat, __LINE__, __FILE__);
    }

/**------------------------------------------------------------------
  leave define mode */
  stat = nc_enddef(ncid);
  nc_check_error(stat, __LINE__, __FILE__);

/**----------------------------------------------------------------------------
  write connectivity table*/
  nv = new int[Element_len * NNPE_len];
  for(m = 0; m < Element_len; m++) {
    for(k = 0; k < NNPE_len; k++)
      nv[m * NNPE_len + k] = discretisation.NIbE[m][k];
      }
  stat = nc_put_var_int(ncid, nodes_id, nv);
  nc_check_error(stat, __LINE__, __FILE__);
  delete[] nv;

/**----------------------------------------------------------------------------
  write levels table: might be done somewhere else */
  if(discretisation.nodes!=0) {
  if(discretisation.nodes[0].zlevels!=0) {
    z = new double[nodes_len * L_len];
    for(k = 0; k < L_len; k++) {
      for(m = 0; m < nodes_len; m++) {
        z[k * nodes_len + m] = discretisation.nodes[m].zlevels[k];
        }
      }
    stat = nc_put_var_double(ncid, levels_id, z);
    nc_check_error(stat, __LINE__, __FILE__);
    delete[] z;
    }
  }
  stat = nc_close(ncid);
  nc_check_error(stat, __LINE__, __FILE__);
  return (0);
}

