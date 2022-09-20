
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
#include "config.h"

#include "tugo-prototypes.h"
#include "tmatrix.h"
#include "fe.h"
#include "geo.h"
#include <cmath>
#include <list>

#include "constants.h"
#include "functions.h"
#include "poc-assertions.h"
#include "matrix.h"
#include "map.h"
#include "poc-netcdf.hpp"


using namespace std;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_MassWeights_T(mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,i,m;
  int n1,n2,n3,n4;
  
  for(i = 0; i < mesh.nedges; i++) {
    mesh.edges[i].mw = 0.0;
    for(k = 0; k < mesh.edges[i].nshared; k++) {
      m = mesh.edges[i].shared[k];
      mesh.edges[i].mw += mesh.edges[i].c * mesh.triangles[m].Area;
      }
    mesh.edges[i].mw /= 3.0;
    }
    
  for(i = 0; i < mesh.nvtxs; i++)
    mesh.vertices[i].mw = 0;

  for(l = 0; l < mesh.ntriangles; l++) {
    n1 = mesh.triangles[l].vertex[0];
    n2 = mesh.triangles[l].vertex[1];
    n3 = mesh.triangles[l].vertex[2];
    mesh.vertices[n1].mw += mesh.triangles[l].Area / 3.;
    mesh.vertices[n2].mw += mesh.triangles[l].Area / 3.;
    mesh.vertices[n3].mw += mesh.triangles[l].Area / 3.;
    }
  return(0);

}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_MassWeights_Q(mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,i,m;
  int n1,n2,n3,n4;
  
  for(i = 0; i < mesh.nedges; i++) {
    mesh.edges[i].mw = 0.0;
    for(k = 0; k < mesh.edges[i].nshared; k++) {
      m = mesh.edges[i].shared[k];
      mesh.edges[i].mw += mesh.edges[i].c * mesh.quadrangles[m].Area;
      }
    mesh.edges[i].mw /= 4.0;
    }
    
  for(i = 0; i < mesh.nvtxs; i++)
    mesh.vertices[i].mw = 0;

  for(l = 0; l < mesh.nquadrangles; l++) {
    n1 = mesh.quadrangles[l].vertex[0];
    n2 = mesh.quadrangles[l].vertex[1];
    n3 = mesh.quadrangles[l].vertex[2];
    n4 = mesh.quadrangles[l].vertex[2];
    mesh.vertices[n1].mw += mesh.quadrangles[l].Area / 4.;
    mesh.vertices[n2].mw += mesh.quadrangles[l].Area / 4.;
    mesh.vertices[n3].mw += mesh.quadrangles[l].Area / 4.;
    mesh.vertices[n3].mw += mesh.quadrangles[l].Area / 4.;
    }
  return(0);

}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_MassWeights(mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  switch (mesh.ntriangles) {
    case 0:
      status=fe_MassWeights_Q(mesh);
      break;
    default:
      status=fe_MassWeights_T(mesh);
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void sidelength02(mesh_t &mesh, int n1, int n2, double *sinus, double *cosinus,double *size, double ref)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
  return the length and !outward! normal orientation of the element side
  [n1,n2]  !!!
  length given in meters
----------------------------------------------------------------------*/
{
  double dX, dY, ds;
  double alpha;
  double lon1, lon2;
  double lat1, lat2;
  double dlon,dlat;

  lon2 = mesh.vertices[n2].lon;
  lon1 = mesh.vertices[n1].lon;
  lon2 = degree_recale(lon2, lon1);
  lat2 = mesh.vertices[n2].lat;
  lat1 = mesh.vertices[n1].lat;

  alpha = (double) 0.5 * (lat1+lat2) * d2r;

  dlon=(lon2-lon1)*d2r;
  dlat=(lat2-lat1)*d2r;

  dX = R * cos(alpha) * dlon;
  dY = R * dlat;

  ds = sqrt(dX * dX + dY * dY);

/*----------------------------------------------------------------------
  outward normal is clockwise rotated */
  *cosinus = +dY / ds;
  *sinus   = -dX / ds;
  *size    = ds;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void sidelength_T(mesh_t *mesh, int n1, int n2, int m, int s, double ref)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  return the length and !outward! normal orientation of the element side
  [n1,n2]; length given in meters
----------------------------------------------------------------------*/
{
  sidelength02(*mesh, n1, n2, &mesh->triangles[m].ny[s],
                              &mesh->triangles[m].nx[s],
                              &mesh->triangles[m].l[s], ref);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void sidelength_Q(mesh_t *mesh, int n1, int n2, int m, int s, double ref)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  return the length and !outward! normal orientation of the element side
  [n1,n2]; length given in meters
----------------------------------------------------------------------*/
{
  sidelength02(*mesh, n1, n2, &mesh->quadrangles[m].ny[s],
                              &mesh->quadrangles[m].nx[s],
                              &mesh->quadrangles[m].l[s], ref);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int read_mesh(const char *type,const char *rootname, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, l, n, status;
  int idum, ndum;
  char BoundaryFile[256], NodeFile[256], EleFile[256], BatFile[256],
        NeighbourFile[256];
  char output[256],msg[256];
  float dummy, fdum;
  int *localmax, k, kb, check;
  poc_var_t var;
  string vname;
  char   *vtopo=0, *vmask=0;
  metagrid_t meta;
  double tag=NAN;


  if(!strncmp(type, "nei", 4)) {
/*------------------------------------------------------------------------------
    input from Falconer Henry style neighbour file (TRIGRID format)*/
    l = strlen(rootname);
    strncpy(NeighbourFile, rootname, l + 1);
    strcat(NeighbourFile, "ll.nei");
    fprintf(echo, "Trying neighbour file:- %s\n", NeighbourFile);

    status = fe_readmesh(NeighbourFile, MESH_FILE_FORMAT_TRIGRID, &mesh);
    if(status != 0) {
      strcpy(NeighbourFile, GridName);
      status = fe_readmesh(NeighbourFile, MESH_FILE_FORMAT_TRIGRID, &mesh);
      }

    if(status != 0) {
      sprintf(msg,"fe_readmesh error: %s",NeighbourFile);
      check_error(status,msg, __LINE__, __FILE__, 1);
      }
      
   switch(gMeshFileType) {
     case FE_TRIANGLE:
       status = fe_list(&mesh );
       break;
       
     case FE_QUADRANGLE:
       status = fe_list(&mesh, FE_QUADRANGLE );
       break;
       
     default:
       check_error(status, "illegal mesh type", __LINE__, __FILE__, 1);
       break;
     }
    if(status != 0)
      check_error(status, "fe_list error", __LINE__, __FILE__, 1);
    }
  else if(!strncmp(type, "nod", 3)) {
/*------------------------------------------------------------------------------
    input from QUODDY style*/
    check_error(-1, "*.nod format not implemented", __LINE__, __FILE__, 1);
    } 
  else if(!strncmp(type, "nc", 2)) {
/*------------------------------------------------------------------------------
    netcdf format*/
   switch(gMeshFileType) {
     case FE_TRIANGLE:
       strcpy(NeighbourFile, rootname);
       fprintf(echo, " Using neighbour file:- %s\n", NeighbourFile);
       status = fe_readmesh3d(NeighbourFile, &mesh, 0);
       if(status != 0) check_error(status, "fe_readmesh3d error", __LINE__, __FILE__, 1);
       status = fe_list(&mesh);
       if(status != 0) check_error(status, "fe_list error", __LINE__, __FILE__, 1);
       break;
       
     case FE_QUADRANGLE:
       check_error(-1, "*.nc format not implemented", __LINE__, __FILE__, 1);
       break;
       
     case FE_STRUCTURED:
       strcpy(NeighbourFile, rootname);
/*------------------------------------------------------------------------------
       try setup key*/
       if(tugo_cfg->model->MeshMeta.face_value()!="NONE")
        parse_GridOptions(tugo_cfg->model->MeshMeta.face_value(), meta);
/*------------------------------------------------------------------------------
       guess*/
       else {
         vname="longitude";
         status=poc_inq_var( NeighbourFile, vname,&var,0);
         if(status==0) meta.z_gnames.vlon=strdup(vname.c_str());
         if(status!=0) {
           vname="lon";
           status=poc_inq_var( NeighbourFile, vname,&var,0);
           if(status==0) meta.z_gnames.vlon=strdup(vname.c_str());
           }
         vname="longitude_f";
         status=poc_inq_var( NeighbourFile, vname,&var,0);
         if(status==0) meta.f_gnames.vlon=strdup(vname.c_str());
         if(status!=0) {
           vname="lon_f";
           status=poc_inq_var( NeighbourFile, vname,&var,0);
           if(status==0) meta.f_gnames.vlon=strdup(vname.c_str());
           }
/*------------------------------------------------------------------------------
         vorticity nodes should be privileged for bathymetry precription*/
         vname="bathymetry_f";
         status=poc_inq_var( NeighbourFile, vname,&var,0);
         if(status==0) meta.z_gnames.vtopo=strdup(vname.c_str());
/*------------------------------------------------------------------------------
         if not the case, then look for pressure node bahymetry*/
         if(status!=0) {
           vname="bathymetry_t";
           status=poc_inq_var( NeighbourFile, vname,&var,0);
           if(status==0) meta.z_gnames.vtopo=strdup(vname.c_str());
           }
         if(status!=0) {
           vname="bathymetry";
           status=poc_inq_var( NeighbourFile, vname,&var,0);
           if(status==0) meta.z_gnames.vtopo=strdup(vname.c_str());
           }
         if(status!=0) {
           vname="H0";
           status=poc_inq_var( NeighbourFile, vname,&var,0);
           if(status==0) meta.z_gnames.vtopo=strdup(vname.c_str());
           }
         vname="landmask";
         status=poc_inq_var( NeighbourFile,vname,&var,0);
         if(status==0) meta.z_gnames.vmask=strdup(vname.c_str());
         }
//       if(meta.z_gnames.vmask==0) tag=0.;
       status=quadrangle_ImportStructured(NeighbourFile, mesh, meta, gCgrid, tag);
       if(status != 0) check_error(status, "quadrangle_ImportStructured error, abort", __LINE__, __FILE__, 1);
       break;
       
     default:
       check_error(status, "illegal mesh type", __LINE__, __FILE__, 1);
       break;
     }

    } 
  else {
    check_error(-1, "unrecognized mesh file extension", __LINE__,__FILE__, 1);
    }

/*----------------------------------------------------------------------
  */
  for(i = 0; i < mesh.nvtxs; i++) {
    mesh.vertices[i].c=cos(mesh.vertices[i].lat*d2r);
    mesh.vertices[i].s=sin(mesh.vertices[i].lat*d2r);
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int geometry(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------

  Partial (2D) geometry initialisation

-----------------------------------------------------------------------------*/
{
  int i, ii, i1, i2, i3, j, k, l, m, n, target;
  int n1, n2, n3, n4;
  int InitializeCodes, StopOnPinched, status;
  double ref,lon1, lon2, lon3, lat1, lat2, lat3;
  float C, hC, hc;
  float temp1, temp2;
  double S, dP, dQ, dZdx, dZdy;
  float h;
  discretisation_t z_descriptor,u_descriptor;

// #ifdef HAVE_MPI
// //#ifdef PARALLEL
//   int context=PARALLEL_COMPUTING;
// #else
//   int context=SEQUENTIAL_COMPUTING;
// #endif
  int context=gLinearSolverMode;

/*----------------------------------------------------------------------
  geometry initialisation */
  for(l = 0; l < mesh->ntriangles; l++) {
    status = fe_initaffine(mesh, l);
    n1 = mesh->triangles[l].vertex[0];
    n2 = mesh->triangles[l].vertex[1];
    n3 = mesh->triangles[l].vertex[2];
/**----------------------------------------------------------------------
    must be consistent with edge's index*/
    target=2;
    sidelength_T(mesh, n1, n2, l, target, mesh->triangles[l].cosinus);
    target=0;
    sidelength_T(mesh, n2, n3, l, target, mesh->triangles[l].cosinus);
    target=1;
    sidelength_T(mesh, n3, n1, l, target, mesh->triangles[l].cosinus);
    }

  for(l = 0; l < mesh->nquadrangles; l++) {
//    status = fe_initaffine(mesh, l); /// HERE !!!
    n1 = mesh->quadrangles[l].vertex[0];
    n2 = mesh->quadrangles[l].vertex[1];
    n3 = mesh->quadrangles[l].vertex[2];
    n4 = mesh->quadrangles[l].vertex[3];
/**----------------------------------------------------------------------
    must be consistent with edge's index*/
    target=0;
    sidelength_Q(mesh, n1, n2, l, target, mesh->quadrangles[l].cosinus);
    target=1;
    sidelength_Q(mesh, n2, n3, l, target, mesh->quadrangles[l].cosinus);
    target=2;
    sidelength_Q(mesh, n3, n4, l, target, mesh->quadrangles[l].cosinus);
    target=3;
    sidelength_Q(mesh, n4, n1, l, target, mesh->quadrangles[l].cosinus);
    }

  /* STS: main_area NOT used... */

/*----------------------------------------------------------------------
  building table of element connected to vertices*/
  status = fe_vertex_element_tables(mesh);

  if(gParallel==0) {
/**----------------------------------------------------------------------
    building edge connectivity*/
    status = fe_edgetable(mesh);
    }
  else {
/**----------------------------------------------------------------------
    update edge tables*/
    status = update_edge_table(mesh);
    }

  for(i = 0; i < mesh->nedges; i++) {
    n1 = mesh->edges[i].extremity[0];
    n2 = mesh->edges[i].extremity[1];
    ref=mesh->vertices[n1].lon;
    mesh->edges[i].lon = 0.5 * (mesh->vertices[n1].lon + geo_recale(mesh->vertices[n2].lon,ref,180.0))*d2r;
    mesh->edges[i].lat = 0.5 * (mesh->vertices[n1].lat + mesh->vertices[n2].lat)*d2r;
    mesh->edges[i].c = cos(mesh->edges[i].lat);
    mesh->edges[i].s = sin(mesh->edges[i].lat);
    }

/*----------------------------------------------------------------------
  building element connectivity*/
  status = fe_element_crosstables(*mesh);

/*----------------------------------------------------------------------
  building mesh limits table*/
  InitializeCodes = 1;
  StopOnPinched   = 0;
  status = fe_codetable_obsolete(mesh, InitializeCodes, StopOnPinched, 0);

/**----------------------------------------------------------------------
  read boundary numerical type from *.bel file*/
  status=discretisation_init(mesh, z2D_discretisation, context);
  status=discretisation_init(mesh, u2D_discretisation, context); /// HERE !!!
  process_boundarycode(*mesh, 0, false);
  
/**----------------------------------------------------------------------
  insert periodic mesh elements (if any)*/
  /* STS: NO periodic */

//   for(i = 0; i < mesh->nedges; i++) {
//     mesh->edges[i].mw = 0.0;
//     for(k = 0; k < mesh->edges[i].nshared; k++) {
//       m = mesh->edges[i].shared[k];
//       mesh->edges[i].mw += mesh->edges[i].c * mesh->triangles[m].Area;
//       }
//     mesh->edges[i].mw /= 3.0;
//     }
  status= fe_MassWeights(*mesh);
  
/**----------------------------------------------------------------------
  initialize computational structures*/
  status=discretisation_init(mesh, z2D_discretisation, context);
  status=discretisation_init(mesh, u2D_discretisation, context);
  status=discretisation_init(mesh, g2D_discretisation, context);

/**----------------------------------------------------------------------
  initialize computational structures for LGP1 outputs*/
  switch (mesh->ntriangles) {
    case 0:
      status=discretisation_init(mesh, CQP1, context);
      break;
    default:
      status=discretisation_init(mesh, LGP1, context);
      break;
    }

  status=paire_discretisation(*mesh, gSequentialPaire2D, &z_descriptor, &u_descriptor);
  for(n=0;n<z_descriptor.nnodes;n++) {
    z_descriptor.nodes[n].ancestro=n;
    }
  for(n=0;n<u_descriptor.nnodes;n++) {
    u_descriptor.nodes[n].ancestro=n;
    }

  mesh->level  = 0;
  mesh->units  = 0;      /*degrees */
  
  return(0);

}
