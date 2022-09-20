
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
#include <string.h>
#include <unistd.h>
#include <cmath>
#include <netcdf.h>

#include "tugo-prototypes.h"
#include "constants.h"
#include "fe.h"
#include "map.h"
#include "parallel.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gravitational_constant_01(discretisation_t & descriptor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  http://en.wikipedia.org/wiki/Gravity_of_Earth
-----------------------------------------------------------------------------*/
{
  int n, status;
  printf("use latitude-dependant gravitational constant\n");
  for(n=0;n<descriptor.nnodes;n++) {
    const double s=descriptor.nodes[n].s;
    const double c=descriptor.nodes[n].c;
    descriptor.nodes[n].g=9.780327*(1+0.0053024*s*s-2*0.0000058*s*c*s*c);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gravitational_constant_00(discretisation_t & descriptor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  http://en.wikipedia.org/wiki/Standard_gravity
-----------------------------------------------------------------------------*/
{
  int n, status;
  printf("use uniform gravitational constant\n");
  for(n=0;n<descriptor.nnodes;n++) {
    descriptor.nodes[n].g=P_g;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gravitational_constant(discretisation_t & descriptor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  if(tugo_cfg->physics->gravity_mode.face_value() == "CONSTANT") {
    gGravityMode=0;
    }
  else if (tugo_cfg->physics->gravity_mode.face_value() == "LATITUDE-DEPENDANT") {
    gGravityMode=1;
    }
  else {
    check_error(-1, "unreckognized option for gravity mode", __LINE__, __FILE__, 1);
    }
  
  switch(gGravityMode) {
    case 0:
      status=gravitational_constant_00(descriptor);
      break;
    case 1:
      status=gravitational_constant_01(descriptor);
      break;
    }
    
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int coriolis_setup(discretisation_t & descriptor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int n;
  double latitude=tugo_cfg->physics->coriolis_latitude.numerical_value<double>();
  double dtr=M_PI/180.;
  
  if (tugo_cfg->physics->coriolis_mode.face_value() == "LATITUDE-DEPENDANT") {
    gCoriolisMode=0;
    }
  else if (tugo_cfg->physics->coriolis_mode.face_value() == "CONSTANT") {
    gCoriolisMode=1;
    }
  else if (tugo_cfg->physics->coriolis_mode.face_value() == "BETA-PLAN") {
    gCoriolisMode=1;
    }
  else {
    check_error(-1, "unreckognized option for gravity mode", __LINE__, __FILE__, 1);
    }
  
  switch(gCoriolisMode) {
    case 0:
      for(n=0;n<descriptor.nnodes;n++) {
        descriptor.nodes[n].s=sin(descriptor.nodes[n].lat*dtr);
        }
      break;
    case 1:
      for(n=0;n<descriptor.nnodes;n++) {
        descriptor.nodes[n].s=sin(latitude*dtr);
        }
      break;
    case 2:
      for(n=0;n<descriptor.nnodes;n++) {
        descriptor.nodes[n].s=sin(latitude*dtr);
        }
      break;
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int coriolis_setup(parameter_t & data, int nnodes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int n;
  double latitude=tugo_cfg->physics->coriolis_latitude.numerical_value<double>();
  double dtr=M_PI/180.;
  
  if (tugo_cfg->physics->coriolis_mode.face_value() == "LATITUDE-DEPENDANT") {
    gCoriolisMode=0;
    }
  else if (tugo_cfg->physics->coriolis_mode.face_value() == "CONSTANT") {
    gCoriolisMode=1;
    }
  else if (tugo_cfg->physics->coriolis_mode.face_value() == "BETA-PLAN") {
    gCoriolisMode=1;
    }
  else {
    check_error(-1, "unreckognized option for gravity mode", __LINE__, __FILE__, 1);
    }
  
  switch(gCoriolisMode) {
    case 0:
      for(n=0;n<nnodes;n++) {
        data.S[n]=sin(data.lat[n]);
        }
      break;
    case 1:
      for(n=0;n<nnodes;n++) {
        data.S[n]=sin(latitude*dtr);
        }
      break;
    case 2:
      for(n=0;n<nnodes;n++) {
        data.S[n]=sin(latitude*dtr);
        }
      break;
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int allocate_action2D(mesh_t & mesh, int level)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,nndes,status;
  int udim,tdim=0;;
  int zdim;

  status=paire_dimension(mesh,&zdim,&udim);
  if(zdim==0)
      check_error(-1, "zdim uninitialized", __LINE__, __FILE__, 1);
  if(udim==0)
      check_error(-1, "udim uninitialized", __LINE__, __FILE__, 1);
  
/**----------------------------------------------------------------------
  create generic forcing vector*/
  gUseReformedAction=1;
  for (k=0;k<gNTF;k++) gAction2D[level][k].allocate(zdim,udim,tdim);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int allocate_data2D(mesh_t mesh,parameter_t *data)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int udim;
  int zdim;

  status=paire_dimension(mesh,&zdim,&udim);
  if(zdim==0)
      check_error(-1, "zdim uninitialized", __LINE__, __FILE__, 1);
  if(udim==0)
      check_error(-1, "udim uninitialized", __LINE__, __FILE__, 1);

  data->lon  = new double[zdim];
  data->lat  = new double[zdim];
  data->h    = new double[zdim];
  
/**----------------------------------------------------------------------------
  preparing new forcing coding (toward true state vector discretisation)*/
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
  data->rlinear = new float[udim];
  data->FrictionRatio   = new float[udim];
  data->vwr  = new float[udim];
  data->wdc1     = new float[udim];
  data->Nbar     = new float[udim];
  data->celerity = new float[udim];

  data->uParent  = new int[udim];
  data->uCode    = new int[udim];

  if(udim==zdim) {
    data->h__  = data->h;
    data->u0__ = data->u0;
    data->colocation=1;
    }
  else {
    data->h__  = new double[udim];
    data->u0__ = new double[zdim];
    data->colocation=0;
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int allocate_state2D(mesh_t mesh,state2d_t *state)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int udim,zdim,g_udim;
  int i;

  status=paire_dimension(mesh,&zdim, &udim, &g_udim);

//   g_udim=max(mesh.ntriangles,udim);
//   g_udim=max(mesh.nedges,g_udim);

  state->H = new double[zdim];
  state->z = new double[zdim];

#ifdef DRY
  state->wet = new int[zdim];
#endif      

  state->P  = new double[zdim];
  state->pb = new double[zdim];

  state->u = new double[udim];
  state->v = new double[udim];

  state->U = new double[udim];
  state->V = new double[udim];

  state->Kh = new double[g_udim];

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Notes:

  22/04/2011:

    allocation changed to accommodate various way of computing Smagorinski
    coeffficient, with possible side effects on advection

    To be reviewed....

\author Damien Allain
\date reviewed 12 May 2011 : bug fixed

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  

  switch(AD_discretisation) {

    case NODE_WISE:
//       state->ux = new double[udim];
//       state->uy = new double[udim];
//       state->vx = new double[udim];
//       state->vy = new double[udim];
      state->ux = new double[g_udim];
      state->uy = new double[g_udim];
      state->vx = new double[g_udim];
      state->vy = new double[g_udim];
      state->Ah = new double[udim];
//      state->lmts = new double[udim];
      state->lmts = new double[g_udim];
      break;

    case ELEMENT_WISE:
      state->ux = new double[g_udim];
      state->uy = new double[g_udim];
      state->vx = new double[g_udim];
      state->vy = new double[g_udim];
      state->Ah = new double[udim];
      state->lmts = new double[g_udim];
      break;

    default:
      check_error(-1, "illicit advection-diffusion discretisation mode", __LINE__, __FILE__, 1);
      break;
    }

  state->Ax = new double[udim];
  state->Ay = new double[udim];

  state->Ex = new double[udim];
  state->Ey = new double[udim];

  state->rtau = new double[udim];
  state->rx   = new double[udim];
  state->ry   = new double[udim];
  state->taux = new double[udim];
  state->tauy = new double[udim];

  state->delta3Dx = new double[udim];
  state->delta3Dy = new double[udim];

  state->vic = new float[udim];
  state->vie = new float[udim];
  state->N   = new float[udim];

  if(udim==zdim) {
    state->H__ = state->H;
    state->z__ = state->z;
    state->u__ = state->u;
    state->v__ = state->v;
    state->U__ = state->U;
    state->V__ = state->V;
    state->Ex__   = state->Ex;
    state->Ey__   = state->Ey;
    state->taux__ = state->taux;
    state->tauy__ = state->tauy;
    state->Ah__   = state->Ah;
    state->colocation=1;
    }
  else {
    state->H__ = new double[udim];
    state->z__ = new double[udim];
    state->u__ = new double[zdim];
    state->v__ = new double[zdim];
    state->U__ = new double[zdim];
    state->V__ = new double[zdim];
    state->Ex__   = new double[zdim];
    state->Ey__   = new double[zdim];
    state->taux__ = new double[zdim];
    state->tauy__ = new double[zdim];
    state->Ah__   = new double[zdim];
    state->colocation=0;
    }

  if (tugo_cfg->dynamics2D->BarotropicW.face_value() == (string) "TRUE") {
    state->wres = new double[zdim];
    state->w    = new double*[zdim];
    for(i = 0; i < zdim; i++) {
      state->w[i]= new double[2];
      }
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int allocate_UGO2D(mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,n1,n2,n3;
  int nlayers,mode,status;
  int udim;
  int zdim;

  status=paire_dimension(mesh,&zdim,&udim);

/**-------------------------------------------------------------------------
  state vectors allocation*/
  status=allocate_data2D(mesh,&(gdata2D[0]));
  
/**-------------------------------------------------------------------------
  state vectors allocation*/
  for(k = 0; k < 3; k++) {
    status=allocate_state2D(mesh,&(gstate2D[0][k]));
    }
  
  status=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int initialise_UGO2D(mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,n1,n2,n3;
  int nlayers,status;
  int udim, zdim, g_udim;
  double h, shallowest, *buffer;

//  status=paire_dimension(mesh,&zdim,&udim);
  status=paire_dimension(mesh,&zdim,&udim, &g_udim);

/**-------------------------------------------------------------------------
  By default, parameters are given at element vertices*/
  switch (z2D_discretisation) {

    case CQP0:
      status=projection_CQP1xCQP0(mesh, gP1data2D.h, gdata2D[0].h);
      for(i = 0; i < mesh.CQP0descriptor.nnodes; i++) {
        gdata2D[0].lon[i] = mesh.CQP0descriptor.nodes[i].lon;
        gdata2D[0].lat[i] = mesh.CQP0descriptor.nodes[i].lat;
        gdata2D[0].lon[i] *= d2r;
        gdata2D[0].lat[i] *= d2r;
        }
      break;

    default:
      check_error(-1, "discretisation not implemented", __LINE__, __FILE__, 1);
      break;
    }

  for(i = 0; i < udim; i++) {
    gdata2D[0].uParent[i] =  i;
    gdata2D[0].uCode[i]   = FE_UNDEFINED_NODE;
    }

  for(i = 0; i < zdim; i++) {
    gdata2D[0].zParent[i]   = i;
    gdata2D[0].zAncestor[i] = i;
    for(k = 0; k < 3; k++) {
      gstate2D[0][k].H[i] = gdata2D[0].h[i];
      gstate2D[0][k].z[i]  = 0.0;
      gstate2D[0][k].P[i]  = 0.0;
      gstate2D[0][k].pb[i] = 0.0;
#ifdef DRY
      gstate2D[0][k].wet[i] = -1;
#endif      
      }
    gdata2D[0].C[i] = cos(gdata2D[0].lat[i]);
    gdata2D[0].S[i] = sin(gdata2D[0].lat[i]);
    gdata2D[0].lmm[i]   = 0.0;
    gdata2D[0].zCode[i] = 0;
    }

/**-----------------------------------------------------------------------------
  overload gdata2D[0].S initialisation, for academic tests */
  status=coriolis_setup(gdata2D[0], zdim);
  
  for(i = 0; i < mesh.ntriangles; i++) {
    gdata2D[0].countdown[i] =  1;
    gdata2D[0].unstable[i]  =  0;
    }

  for(i = 0; i < udim; i++) {
    for(k = 0; k < 3; k++) {
      gstate2D[0][k].u[i] = 0.0;
      gstate2D[0][k].v[i] = 0.0;
      gstate2D[0][k].U[i] = 0.0;
      gstate2D[0][k].V[i] = 0.0;
      gstate2D[0][k].Ax[i] = 0.0;
      gstate2D[0][k].Ay[i] = 0.0;
      gstate2D[0][k].Ex[i] = 0.0;
      gstate2D[0][k].Ey[i] = 0.0;
/**----------------------------------------------------------------------------
  Ah is implicitely u node difusion coefficient, Kh is u-gradient node value */
      gstate2D[0][k].Ah[i] = 0.0;
//      gstate2D[0][k].H__[i] = 0.0;
      gstate2D[0][k].taux[i] = 0.0;
      gstate2D[0][k].tauy[i] = 0.0;
      gstate2D[0][k].delta3Dx[i] = 0.0;
      gstate2D[0][k].delta3Dy[i] = 0.0;
//      gstate2D[0][k].lmts[i] = 0.0;
      }
    }
/**----------------------------------------------------------------------------
  Ah is implicitely u node difusion coefficient, Kh is u-gradient node value */
  for(i = 0; i < g_udim; i++) {
    for(k = 0; k < 3; k++) {
      gstate2D[0][k].Kh[i] = 0.0;
      }
    }

  switch (AD_discretisation) {
    case ELEMENT_WISE:
      for(i = 0; i < gFEmesh[0].ntriangles; i++) gstate2D[0][1].lmts[i]=0.0;
      for(i = 0; i < gFEmesh[0].ntriangles; i++) gstate2D[0][1].Kh[i]=0.0;
      break;
    case NODE_WISE:
      for(i = 0; i < udim; i++) gstate2D[0][1].lmts[i]=0.0;
      break;
    }

  if(gstate2D[0][0].colocation==0) {
    for(k = 0; k < 3; k++) {
      for(i = 0; i < zdim; i++) {
        gstate2D[0][k].u__[i] = 0.0;
        gstate2D[0][k].v__[i] = 0.0;
        }
      for(i = 0; i < udim; i++) {
        gstate2D[0][k].H__[i] = 0.0;
        }
      }
    }

  if (tugo_cfg->dynamics2D->BarotropicW.face_value() == (string) "TRUE") {
    for(i = 0; i < zdim; i++) {
      for(k = 0; k < 3; k++) {
        gstate2D[0][k].wres[i] = 0.0;
        gstate2D[0][k].w[i][0] = 0.0;
        gstate2D[0][k].w[i][1] = 0.0;
        }
      }
    }

  for(i = 0; i < zdim; i++) {
    gdata2D[0].Htide[i].a = NULL;
    gdata2D[0].Htide[i].G = NULL;
    gdata2D[0].Utide[i].a = NULL;
    gdata2D[0].Utide[i].G = NULL;
    gdata2D[0].Vtide[i].a = NULL;
    gdata2D[0].Vtide[i].G = NULL;
    }

  for(i = 0; i < gFEmesh[0].nvtxs; i++)
    gFEmesh[0].vertices[i].mw = 0;

  for(l = 0; l < gFEmesh[0].ntriangles; l++) {
    n1 = gFEmesh[0].triangles[l].vertex[0];
    n2 = gFEmesh[0].triangles[l].vertex[1];
    n3 = gFEmesh[0].triangles[l].vertex[2];
    gFEmesh[0].vertices[n1].mw += gFEmesh[0].triangles[l].Area / 3.;
    gFEmesh[0].vertices[n2].mw += gFEmesh[0].triangles[l].Area / 3.;
    gFEmesh[0].vertices[n3].mw += gFEmesh[0].triangles[l].Area / 3.;
    }

  switch (z2D_discretisation) {

    case LGP1:
      for(i = 0; i < gFEmesh[0].nvtxs; i++)
        gdata2D[0].lmm[i]=gFEmesh[0].vertices[i].mw;
      break;

    default:
      for(i = 0; i < zdim; i++)
        gdata2D[0].lmm[i]=1./0.;
      break;
    }
    
  switch (u2D_discretisation) {
    case LGP0:
//      gaction2D_obsolete  = 0;
      status=fe_projection(mesh, gP1data2D.z0,   LGP1, gdata2D[0].z0,   u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.Cd,   LGP1, gdata2D[0].Cd,   u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.rlinear,    LGP1, gdata2D[0].rlinear,    u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.FrictionRatio,   LGP1, gdata2D[0].FrictionRatio,   u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.vwr,  LGP1, gdata2D[0].vwr,  u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.wdc1, LGP1, gdata2D[0].wdc1, u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.u0,   LGP1, gdata2D[0].u0,   u2D_discretisation);
      for(k=0;k<3;k++) {
        status=fe_projection(mesh, gP1state2D.vic, LGP1, gstate2D[0][k].vic, u2D_discretisation);
        status=fe_projection(mesh, gP1state2D.vie, LGP1, gstate2D[0][k].vie, u2D_discretisation);
        status=fe_projection(mesh, gP1state2D.N,   LGP1, gstate2D[0][k].N,   u2D_discretisation);
        }
      break;

    case LGP1:
      for(n = 0; n < gFEmesh[0].nvtxs; n++) {
        gdata2D[0].z0[n]   = gP1data2D.z0[n];
        gdata2D[0].Cd[n]   = gP1data2D.Cd[n];
        gdata2D[0].rlinear[n]    = gP1data2D.rlinear[n];
        gdata2D[0].FrictionRatio[n]   = gP1data2D.FrictionRatio[n];
        gdata2D[0].vwr[n]  = gP1data2D.vwr[n];
        gdata2D[0].wdc1[n] = gP1data2D.wdc1[n];
        gdata2D[0].u0[n]   = gP1data2D.u0[n];
//        for(k=0;k<3;k++) gstate2D[0][k].u0[n]  = gP1state2D.u0[n];
        for(k=0;k<3;k++) gstate2D[0][k].vic[n] = gP1state2D.vic[n];
        for(k=0;k<3;k++) gstate2D[0][k].vie[n] = gP1state2D.vie[n];
        for(k=0;k<3;k++) gstate2D[0][k].N[n]   = gP1state2D.N[n];
        }
      break;

    case DGP1:
//      gaction2D_obsolete  = 0;
      status=fe_projection(mesh, gP1data2D.z0,   LGP1, gdata2D[0].z0,   u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.Cd,   LGP1, gdata2D[0].Cd,   u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.rlinear,    LGP1, gdata2D[0].rlinear,    u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.FrictionRatio,   LGP1, gdata2D[0].FrictionRatio,   u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.vwr,  LGP1, gdata2D[0].vwr,  u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.wdc1, LGP1, gdata2D[0].wdc1, u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.u0,   LGP1, gdata2D[0].u0,   u2D_discretisation);
      for(k=0;k<3;k++) {
        status=fe_projection(mesh, gP1state2D.vic, LGP1, gstate2D[0][k].vic, u2D_discretisation);
        status=fe_projection(mesh, gP1state2D.vie, LGP1, gstate2D[0][k].vie, u2D_discretisation);
        status=fe_projection(mesh, gP1state2D.N,   LGP1, gstate2D[0][k].N,   u2D_discretisation);
        }
      break;

    case DNP1:
//      gaction2D_obsolete  = 0;
      status=fe_projection(mesh, gP1data2D.z0,   LGP1, gdata2D[0].z0,   u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.Cd,   LGP1, gdata2D[0].Cd,   u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.rlinear,    LGP1, gdata2D[0].rlinear,    u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.FrictionRatio,   LGP1, gdata2D[0].FrictionRatio,   u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.vwr,  LGP1, gdata2D[0].vwr,  u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.wdc1, LGP1, gdata2D[0].wdc1, u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.u0,   LGP1, gdata2D[0].u0,   u2D_discretisation);
      for(k=0;k<3;k++) {
        status=fe_projection(mesh, gP1state2D.vic, LGP1, gstate2D[0][k].vic, u2D_discretisation);
        status=fe_projection(mesh, gP1state2D.vie, LGP1, gstate2D[0][k].vie, u2D_discretisation);
        status=fe_projection(mesh, gP1state2D.N,   LGP1, gstate2D[0][k].N,   u2D_discretisation);
        }
      break;

    case NCP1:
      for(i = 0; i < gFEmesh[0].nedges; i++) {
        n1 = gFEmesh[0].edges[i].extremity[0];
        n2 = gFEmesh[0].edges[i].extremity[1];
//        gdata2D[0].u0[i] =    0.5 * (gP1data2D.u0[n1] + gP1data2D.u0[n2]);
        gdata2D[0].z0[i]   = 0.5 * (gP1data2D.z0[n1]   + gP1data2D.z0[n2]);
        gdata2D[0].Cd[i]   = 0.5 * (gP1data2D.Cd[n1]   + gP1data2D.Cd[n2]);
        gdata2D[0].rlinear[i]    = 0.5 * (gP1data2D.rlinear[n1]    + gP1data2D.rlinear[n2]);
        gdata2D[0].FrictionRatio[i]   = 0.5 * (gP1data2D.FrictionRatio[n1]   + gP1data2D.FrictionRatio[n2]);
        gdata2D[0].vwr[i]  = 0.5 * (gP1data2D.vwr[n1]  + gP1data2D.vwr[n2]);
        gdata2D[0].wdc1[i] = 0.5 * (gP1data2D.wdc1[n1] + gP1data2D.wdc1[n2]);
        gdata2D[0].u0[i]   = 0.5 * (gP1data2D.u0[n1]   + gP1data2D.u0[n2]);
//        for(k=0;k<3;k++) gstate2D[0][k].u0[i]  =0.5 * (gP1state2D.u0[n1]  + gP1state2D.u0[n2]);
        for(k=0;k<3;k++) gstate2D[0][k].vic[i] =0.5 * (gP1state2D.vic[n1] + gP1state2D.vic[n2]);
        for(k=0;k<3;k++) gstate2D[0][k].vie[i] =0.5 * (gP1state2D.vie[n1] + gP1state2D.vie[n2]);
        for(k=0;k<3;k++) gstate2D[0][k].N[i]   =0.5 * (gP1state2D.N[n1]   + gP1state2D.N[n2]);
        }
      break;

    case CQP1:
      for(n = 0; n < gFEmesh[0].nvtxs; n++) {
        gdata2D[0].z0[n]   = gP1data2D.z0[n];
        gdata2D[0].Cd[n]   = gP1data2D.Cd[n];
        gdata2D[0].rlinear[n]    = gP1data2D.rlinear[n];
        gdata2D[0].FrictionRatio[n]   = gP1data2D.FrictionRatio[n];
        gdata2D[0].vwr[n]  = gP1data2D.vwr[n];
        gdata2D[0].wdc1[n] = gP1data2D.wdc1[n];
        gdata2D[0].u0[n]   = gP1data2D.u0[n];
//        for(k=0;k<3;k++) gstate2D[0][k].u0[n]  = gP1state2D.u0[n];
        for(k=0;k<3;k++) gstate2D[0][k].vic[n] = gP1state2D.vic[n];
        for(k=0;k<3;k++) gstate2D[0][k].vie[n] = gP1state2D.vie[n];
        for(k=0;k<3;k++) gstate2D[0][k].N[n]   = gP1state2D.N[n];
        }
      break;

    case CQN1:
      status=fe_projection(mesh, gP1data2D.z0,   CQP1, gdata2D[0].z0,   u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.Cd,   CQP1, gdata2D[0].Cd,   u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.rlinear,    CQP1, gdata2D[0].rlinear,    u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.FrictionRatio,   CQP1, gdata2D[0].FrictionRatio,   u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.vwr,  CQP1, gdata2D[0].vwr,  u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.wdc1, CQP1, gdata2D[0].wdc1, u2D_discretisation);
      status=fe_projection(mesh, gP1data2D.u0,   CQP1, gdata2D[0].u0,   u2D_discretisation);
      for(k=0;k<3;k++) {
        status=fe_projection(mesh, gP1state2D.vic, CQP1, gstate2D[0][k].vic, u2D_discretisation);
        status=fe_projection(mesh, gP1state2D.vie, CQP1, gstate2D[0][k].vie, u2D_discretisation);
        status=fe_projection(mesh, gP1state2D.N,   CQP1, gstate2D[0][k].N,   u2D_discretisation);
        }
      break;
    default:
      check_error(-1, "discretisation not implemented", __LINE__, __FILE__, 1);
      break;
    }
    

  status=fe_projection(mesh, gP1data2D.celerity, LGP1, gdata2D[0].celerity, u2D_discretisation);

  status=fe_projection(mesh, gP1data2D.dhdx, LGP1, gdata2D[0].dhdx, u2D_discretisation);
  status=fe_projection(mesh, gP1data2D.dhdy, LGP1, gdata2D[0].dhdy, u2D_discretisation);
  
  status=fe_projection(mesh, gdata2D[0].h, z2D_discretisation, gdata2D[0].h__, u2D_discretisation);

  for(k=0;k<3;k++) {
    update_uNodes_elevation(gFEmesh[0],gstate2D[0][k], gdata2D[0]);
    }

#ifdef DRY
/*------------------------------------------------------------------------------
  Initialization of the elevation */
  int max_nb   = mesh.nnghm;
  for(i = 0; i < mesh.nvtxs; i++) {
    shallowest = 0.0;
    h = gdata2D[0].h[i];
/*------------------------------------------------------------------------------
    shallowest is upward oriented variable, always positive (bounded to zero).
    If it reaches 0, it sticks with it.
    if all surrounding nodes are dry, (positive) altitude of the lowest one
    0 elsewhere
    */
    shallowest = max(shallowest, -0.01 - h);
    for(j = 0; j < max_nb; j++) {
      int ii = gFEmesh[0].vertices[i].ngh[j];
      if(ii > -1) {
        h = gdata2D[0].h[ii];
        shallowest = min(shallowest, max(0.0, -0.01 - h));
        }
      }
    gstate2D[0][2].H[i] = shallowest + gdata2D[0].h[i];
    gstate2D[0][1].H[i] = shallowest + gdata2D[0].h[i];
    gstate2D[0][0].H[i] = shallowest + gdata2D[0].h[i];
    }
#endif

/**----------------------------------------------------------------------------
  generic forcing initialization*/
  if(gUseReformedAction) {
  for(n = 0; n < udim; n++) {
    for(k = 1; k < 3; k++) {
      gAction2D[0][k].prx[n] = 0.0;
      gAction2D[0][k].pry[n] = 0.0;
      gAction2D[0][k].tpx[n] = 0.0;
      gAction2D[0][k].tpy[n] = 0.0;
      gAction2D[0][k].wdx[n] = 0.0;
      gAction2D[0][k].wdy[n] = 0.0;
      gAction2D[0][k].wsx[n] = 0.0;
      gAction2D[0][k].wsy[n] = 0.0;
      gAction2D[0][k].RSdvg_x[n] = 0.0;
      gAction2D[0][k].RSdvg_y[n] = 0.0;
      }
    }
  for(n = 0; n < zdim; n++) {
    for(k = 1; k < 3; k++) {
      gAction2D[0][k].Pa[n]  = 0.0;
      gAction2D[0][k].LSA[n] = 0.0;
      }
    }
  }
  status=0;
  return(status);
}

