
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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <cmath>
#include <netcdf.h>

#include "config.h"
#include "tugo-prototypes.h"
#include "constants.h"
#include "fe.h"
#include "map.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int parameters_allocate2D(mesh_t mesh, int paire, parameter_t *data)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int udim;
  int zdim;

  status=paire_dimension(mesh,paire, &zdim,&udim);
  if(zdim==0)
      check_error(-1, "zdim uninitialized", __LINE__, __FILE__, 1);
  if(udim==0)
      check_error(-1, "udim uninitialized", __LINE__, __FILE__, 1);

  data->lon  = new double[zdim];
  data->lat  = new double[zdim];
  data->h    = new double[zdim];
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Notes:

  11/12/2010:

    Forcing terms: should be given at u nodes?

  17/02/2012:

    allocation is wrong, but did not trigger any error when checking
    memory use... strange...
  
  was
  
  data->dhdx = new double[zdim];
  data->dhdy = new double[zdim];
  
  should be 

  data->dhdx = new double[udim];
  data->dhdy = new double[udim];
  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

//   data->dhdx = new double[zdim];
//   data->dhdy = new double[zdim];
  data->dhdx = new double[udim];
  data->dhdy = new double[udim];

  data->C    = new double[zdim];
  data->S    = new double[zdim];

  data->lmm  = new float[zdim];

  data->zAncestor     = new int[zdim];
  data->zParent       = new int[zdim];
  data->zCode         = new int[zdim];

/**----------------------------------------------------------------------------
  stability control*/
//   data->countdown = new int[zdim];
//   data->unstable  = new int[zdim];
  data->countdown = new int[mesh.ntriangles];
  data->unstable  = new int[mesh.ntriangles];

  data->sponge = new int[zdim];

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Notes:

  30/10/2009:

    the strategy for tidal analysis discretisation stiil to be fixed:
    - LGP1 nodes ?
    - natural elevation and velocity nodes ?

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  data->Htide    = new hconstant_t[zdim];
  data->Utide    = new hconstant_t[zdim];
  data->Vtide    = new hconstant_t[zdim];

  data->LSA      = new hconstant_t[zdim];

  data->u0   = new double[udim];

  data->z0   = new float[udim];
  data->Cd   = new float[udim];
  data->rlinear    = new float[udim];
  data->FrictionRatio = new float[udim];
  data->vwr      = new float[udim];
  data->wdc1     = new float[udim];
  data->celerity = new float[udim];
  data->Nbar     = new float[udim];

  data->uParent  = new int[udim];
  data->uCode    = new int[udim];

  if(udim==zdim) {
    data->u0__ = data->u0;
    data->colocation=1;
    }
  else {
    data->u0__ = new double[zdim];
    data->colocation=0;
    }

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  parameter_t parameters_initialise2D(mesh_t mesh,int paire)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,n1,n2,n3,nLGP1,nNCP1;
  int nlayers,status;
  int udim, zdim;
  int z_discretisation,u_discretisation;
  double h, shallowest, *buffer;
  parameter_t local;

  status=paire_dimension(mesh,paire,&zdim,&udim);
  status=paire_discretisation_id(paire, &z_discretisation,&u_discretisation);

  status=parameters_allocate2D( mesh, paire, &local);

/**-------------------------------------------------------------------------
  By default, parameters are given at element vertices*/
  switch (z_discretisation) {
    case LGP0:
      status=projection_LGP1xLGP0(gP1data2D.h,   local.h,mesh);
      status=projection_LGP1xLGP0(gP1data2D.dhdx,local.dhdx,mesh);
      status=projection_LGP1xLGP0(gP1data2D.dhdy,local.dhdy,mesh);
/**-------------------------------------------------------------------------
      projection fails if field has a periodic nature*/
//       buffer=new double[mesh.nvtxs];
//       for(i = 0; i < mesh.nvtxs; i++)  buffer[i] = mesh.vertices[i].lon*d2r;
//       status=projection_LGP1xLGP0(buffer,   local.lon,mesh);
//       for(i = 0; i < mesh.nvtxs; i++)  buffer[i] = mesh.vertices[i].lat*d2r;
//       status=projection_LGP1xLGP0(buffer,   local.lat,mesh);
//       delete[] buffer;
      for(i = 0; i < mesh.LGP0descriptor.nnodes; i++) {
        local.lon[i] = mesh.LGP0descriptor.nodes[i].lon;
        local.lat[i] = mesh.LGP0descriptor.nodes[i].lat;
        local.lon[i] *= d2r;
        local.lat[i] *= d2r;
        }
      break;

    case LGP1:
      for(i = 0; i < mesh.nvtxs; i++) {
        local.lon[i] = mesh.vertices[i].lon;
        local.lat[i] = mesh.vertices[i].lat;
        local.lon[i] *= d2r;
        local.lat[i] *= d2r;
        local.h[i]     = mesh.vertices[i].h;
        local.dhdx[i]  = gP1data2D.dhdx[i];
        local.dhdy[i]  = gP1data2D.dhdy[i];
        }
      break;

    case LGP2:
/**-------------------------------------------------------------------------
      not necessary if later done by cefmo_topography*/
/**-------------------------------------------------------------------------
      need to handle more properly in use out of cefmo context*/   /// HERE !!!
//       status=projection_LGP1xLGP2(gP1data2D.h,   local.h,mesh);
//       status=projection_LGP1xLGP2(gP1data2D.dhdx,local.dhdx,mesh);
//       status=projection_LGP1xLGP2(gP1data2D.dhdy,local.dhdy,mesh);
/**-------------------------------------------------------------------------
      projection fails if field has a periodic nature*/
//       buffer=new double[mesh.nvtxs];
//       for(i = 0; i < mesh.nvtxs; i++)  buffer[i] = mesh.vertices[i].lon*d2r;
//       status=projection_LGP1xLGP2(buffer,   local.lon,mesh);
//       for(i = 0; i < mesh.nvtxs; i++)  buffer[i] = mesh.vertices[i].lat*d2r;
//       status=projection_LGP1xLGP2(buffer,   local.lat,mesh);
//       delete[] buffer;
      for(i = 0; i < mesh.LGP2descriptor.nnodes; i++) {
        local.lon[i] = mesh.LGP2descriptor.nodes[i].lon;
        local.lat[i] = mesh.LGP2descriptor.nodes[i].lat;
        local.lon[i] *= d2r;
        local.lat[i] *= d2r;
        }
      break;

    case CQP0:
      for(i = 0; i < mesh.CQP0descriptor.nnodes; i++) {
        local.lon[i] = mesh.CQP0descriptor.nodes[i].lon;
        local.lat[i] = mesh.CQP0descriptor.nodes[i].lat;
        local.lon[i] *= d2r;
        local.lat[i] *= d2r;
        }
      break;

    default:
      check_error(-1, "discretisation not implemented", __LINE__, __FILE__, 1);
      break;
    }

  for(i = 0; i < udim; i++) {
    local.uParent[i] =  i;
    local.uCode[i]   = FE_UNDEFINED_NODE;
    }

  for(i = 0; i < zdim; i++) {
    local.zParent[i]   = i;
    local.zAncestor[i] = i;
    local.C[i] = cos(local.lat[i]);
    local.S[i] = sin(local.lat[i]);
    local.lmm[i]   = 0.0;
    local.zCode[i] = 0;
    }
  for(i = 0; i < mesh.ntriangles; i++) {
    local.countdown[i] =  1;
    local.unstable[i]  =  0;
    }

  for(i = 0; i < zdim; i++) {
    local.Htide[i].a = NULL;
    local.Htide[i].G = NULL;
    local.Utide[i].a = NULL;
    local.Utide[i].G = NULL;
    local.Vtide[i].a = NULL;
    local.Vtide[i].G = NULL;
    }

  switch (z_discretisation) {

    case LGP1:
      for(i = 0; i < mesh.nvtxs; i++)
        local.lmm[i]=mesh.vertices[i].mw;
      break;

    default:
      for(i = 0; i < zdim; i++)
        local.lmm[i]=1./0.;
      break;
    }

  switch (u_discretisation) {

    case LGP0:

//       status=projection_LGP0(gP1data2D.z0,  local.z0,   mesh, u_discretisation);
//       status=projection_LGP0(gP1data2D.Cd,  local.Cd,   mesh, u_discretisation);
//       status=projection_LGP0(gP1data2D.vwr, local.vwr,  mesh, u_discretisation);
//       status=projection_LGP0(gP1data2D.wdc1,local.wdc1, mesh, u_discretisation);
//       status=projection_LGP0(gP1data2D.u0,  local.u0,   mesh, u_discretisation);
      status=fe_projection(mesh,gP1data2D.z0,   LGP1, local.z0,   u_discretisation);
      status=fe_projection(mesh,gP1data2D.Cd,   LGP1, local.Cd,   u_discretisation);
      status=fe_projection(mesh,gP1data2D.vwr,  LGP1, local.vwr,  u_discretisation);
      status=fe_projection(mesh,gP1data2D.wdc1, LGP1, local.wdc1, u_discretisation);
      status=fe_projection(mesh,gP1data2D.u0,   LGP1, local.u0,   u_discretisation);
      break;

    case LGP1:
      for(n = 0; n < mesh.nvtxs; n++) {
        local.z0[n]   = gP1data2D.z0[n];
        local.Cd[n]   = gP1data2D.Cd[n];
        local.vwr[n]  = gP1data2D.vwr[n];
        local.wdc1[n] = gP1data2D.wdc1[n];
        local.u0[n]   = gP1data2D.u0[n];
        }
      break;

    case NCP1:

      for(i = 0; i < mesh.nedges; i++) {
        n1 = mesh.edges[i].extremity[0];
        n2 = mesh.edges[i].extremity[1];
        local.z0[i]   = 0.5 * (gP1data2D.z0[n1]   + gP1data2D.z0[n2]);
        local.Cd[i]   = 0.5 * (gP1data2D.Cd[n1]   + gP1data2D.Cd[n2]);
        local.vwr[i]  = 0.5 * (gP1data2D.vwr[n1]  + gP1data2D.vwr[n2]);
        local.wdc1[i] = 0.5 * (gP1data2D.wdc1[n1] + gP1data2D.wdc1[n2]);
        local.u0[i]   = 0.5 * (gP1data2D.u0[n1]   + gP1data2D.u0[n2]);
        }
      break;

    case DGP1:

      for(m = 0; m < mesh.ntriangles; m++) {
        for(k=0;k<mesh.DGP1descriptor.nnpe;k++) {
          n=mesh.DGP1descriptor.NIbE[m][k];
          nLGP1=mesh.triangles[m].vertex[k];
          local.z0[n]   = gP1data2D.z0[nLGP1];
          local.Cd[n]   = gP1data2D.Cd[nLGP1];
          local.vwr[n]  = gP1data2D.vwr[nLGP1];
          local.wdc1[n] = gP1data2D.wdc1[nLGP1];
          local.u0[n]   = gP1data2D.u0[nLGP1];
          }
        }
      break;

    case DNP1:

      for(m = 0; m < mesh.ntriangles; m++) {
        for(k=0;k<mesh.DNP1descriptor.nnpe;k++) {
          n=mesh.DNP1descriptor.NIbE[m][k];
          nNCP1=mesh.triangles[m].edges[k];
          n1 = mesh.edges[nNCP1].extremity[0];
          n2 = mesh.edges[nNCP1].extremity[1];
          local.z0[n]   = 0.5 * (gP1data2D.z0[n1]   + gP1data2D.z0[n2]);
          local.Cd[n]   = 0.5 * (gP1data2D.Cd[n1]   + gP1data2D.Cd[n2]);
          local.vwr[n]  = 0.5 * (gP1data2D.vwr[n1]  + gP1data2D.vwr[n2]);
          local.wdc1[n] = 0.5 * (gP1data2D.wdc1[n1] + gP1data2D.wdc1[n2]);
          local.u0[n]   = 0.5 * (gP1data2D.u0[n1]   + gP1data2D.u0[n2]);
          }
        }
      break;

    case IPG:
      status=fe_projection(mesh,gP1data2D.z0,   LGP1, local.z0,   u_discretisation);
      status=fe_projection(mesh,gP1data2D.Cd,   LGP1, local.Cd,   u_discretisation);
      status=fe_projection(mesh,gP1data2D.vwr,  LGP1, local.vwr,  u_discretisation);
      status=fe_projection(mesh,gP1data2D.wdc1, LGP1, local.wdc1, u_discretisation);
      status=fe_projection(mesh,gP1data2D.u0,   LGP1, local.u0,   u_discretisation);
      break;

    case LGP2:

      status=fe_projection(mesh,gP1data2D.z0,   LGP1, local.z0,   u_discretisation);
      status=fe_projection(mesh,gP1data2D.Cd,   LGP1, local.Cd,   u_discretisation);
      status=fe_projection(mesh,gP1data2D.vwr,  LGP1, local.vwr,  u_discretisation);
      status=fe_projection(mesh,gP1data2D.wdc1, LGP1, local.wdc1, u_discretisation);
      status=fe_projection(mesh,gP1data2D.u0,   LGP1, local.u0,   u_discretisation);
      break;

    case CQP1:

      status=fe_projection(mesh,gP1data2D.z0,   CQP1, local.z0,   u_discretisation);
      status=fe_projection(mesh,gP1data2D.Cd,   CQP1, local.Cd,   u_discretisation);
      status=fe_projection(mesh,gP1data2D.vwr,  CQP1, local.vwr,  u_discretisation);
      status=fe_projection(mesh,gP1data2D.wdc1, CQP1, local.wdc1, u_discretisation);
      status=fe_projection(mesh,gP1data2D.u0,   CQP1, local.u0,   u_discretisation);
      break;

    case CQN1:

      status=fe_projection(mesh,gP1data2D.z0,   CQP1, local.z0,   u_discretisation);
      status=fe_projection(mesh,gP1data2D.Cd,   CQP1, local.Cd,   u_discretisation);
      status=fe_projection(mesh,gP1data2D.vwr,  CQP1, local.vwr,  u_discretisation);
      status=fe_projection(mesh,gP1data2D.wdc1, CQP1, local.wdc1, u_discretisation);
      status=fe_projection(mesh,gP1data2D.u0,   CQP1, local.u0,   u_discretisation);
      break;

    default:
      check_error(-1, "discretisation not implemented", __LINE__, __FILE__, 1);
      break;
    }
  status=fe_projection(mesh, gP1data2D.celerity, LGP1, local.celerity, u_discretisation);
  status=fe_projection(mesh, gP1data2D.Nbar,     LGP1, local.Nbar,     u_discretisation);

  status=0;
  return(local);
}

