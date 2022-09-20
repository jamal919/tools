
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

/**----------------------------------------------------------------------------

  Tidal spectral solver (cefmo type)

  solves the quasi-linearized complex equations

  Available discretisation : LGP0xLGP1 and DLGP1xLGP2

  Mostly identical to CEFMO, except for:

    -numerical integration (exact here, Gauss quadrature in CEFMO)

    -bathymetry discretisation (computational nodes here, Gauss points in CEFMO)

    -M2 and K1 velocity mix (for friction coefficient computation) not implemented yet

-----------------------------------------------------------------------------*/
#include "tugo-prototypes.h"
#include <string>
#include <unistd.h>
#include "cefmo.h"
#include "matrix.h"

extern int gravitational_constant(discretisation_t & descriptor);

/**----------------------------------------------------------------------------
for dominant wave iteration acceleration*/
#define ACCELERE
//double *RO[2]={0,0},*ROPRIM[2]={0,0};
//int counter=0;

atlas2D_t cefmo_atlas2D;
atlas3D_t cefmo_atlas3D;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int cefmo_allocate(mesh_t & mesh, cefmo_t & cefmo, int already_initialised)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  int z_discretisation, u_discretisation;
  int zdim,udim;
  
  status=paire_dimension(mesh,cefmo.paire,&zdim,&udim);
  status=paire_discretisation_id(cefmo.paire, &z_discretisation, &u_discretisation);

/**----------------------------------------------------------------------------
  bathymetry at elevation and velocity nodes */
  cefmo.state.h_z=new double[zdim];
  cefmo.state.h_u=new double[udim];

/**----------------------------------------------------------------------------
  complex state vector */
  cefmo.state.u=new complex <double>[udim];
  cefmo.state.v=new complex <double>[udim];
  cefmo.state.z=new complex <double>[zdim];

/**----------------------------------------------------------------------------
  complex tidal potential */
  cefmo.state.LSA=new complex <double>[zdim];
  cefmo.state.potential=new complex <double>[zdim];

  mesh.nlayers=1;
  mesh.nlevels=mesh.nlayers+1;

  if (already_initialised==0){
    cefmo.action.prx=new complex<double>*[udim];
    cefmo.action.pry=new complex<double>*[udim];
    cefmo.action.Fx=new complex<double>[udim];
    cefmo.action.Fy=new complex<double>[udim];
    for(n = 0; n < udim; n++) {
      cefmo.action.prx[n]=new complex<double>[mesh.nlayers];
      cefmo.action.pry[n]=new complex<double>[mesh.nlayers];
      cefmo.action.Fx[n]=0.0;
      cefmo.action.Fy[n]=0.0;
      }
    }
  
  for(n = 0; n < zdim; n++) {
    cefmo.state.z[n]=0.0;
    cefmo.state.LSA[n]=0.0;
    cefmo.state.potential[n]=0.0;
    }

/**-----------------------------------------------------------------------
  start with strong currenst to initiate iterative process*/
  for(n = 0; n < udim; n++) {
    cefmo.state.u[n]=complex <double>(0.5,0.0);
    cefmo.state.v[n]=complex <double>(0.5,0.0);
    }
  
  gPMLattenuation_x=0;
  gPMLattenuation_y=0;
  gBufferAttenuation=0;

  return(0);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int SpArchive_Get(const char *filename, mesh_t mesh, cefmo_t cefmo, tide2D_t state)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /**-----------------------------------------------------------------------
// 
//   Spectral solution archive
// 
//  -----------------------------------------------------------------------*/
// {
//   int i,j,k,l,m,n,status;
//   int ncid,varid,frame=0;
//   int zdim,udim;
//   int z_discretisation, u_discretisation;
//   int wave_id;
//   double *buffer[2];
//   char varname[1024];
// 
//   status=paire_dimension(mesh,cefmo.paire,&zdim,&udim);
//   status=paire_discretisation_id(cefmo.paire, &z_discretisation, &u_discretisation);
//   
//   const char *UNAME=discretisation_name(u_discretisation);
//   const char *ZNAME=discretisation_name(z_discretisation);
// 
// 
//   buffer[0]=new double [zdim];
//   buffer[1]=new double [zdim];
// 
//   status = nc_open(filename, NC_NOWRITE, &ncid);
//   if(status != 0) goto error;
// 
// /**----------------------------------------------------------------------------
//   elevation */
//   status = nc_inq_varid(ncid, "a_eta_LGP1", &varid);
//   if(status != 0)  goto error;
//   status = poc_get_UG3D(ncid, mesh, frame, varid, buffer[0]);
//   if(status != 0)  goto error;
//   status = nc_inq_varid(ncid, "G_eta_LGP1", &varid);
//   if(status != 0)  goto error;
//   status = poc_get_UG3D(ncid, mesh, frame, varid, buffer[0]);
//   if(status != 0)  goto error;
//   for(n = 0; n < zdim; n++) {
//     state.z[n]=polar(buffer[0][n],-buffer[1][n]*M_PI/180.0);
//     }
// 
// /**----------------------------------------------------------------------------
//   currents */
//   status = nc_inq_varid(ncid, "a_u_LGP0", &varid);
//   if(status != 0)  goto error;
//   status = poc_get_UG3D(ncid, mesh, frame, varid, buffer[0]);
//   if(status != 0)  goto error;
//   status = nc_inq_varid(ncid, "G_u_LGP0", &varid);
//   if(status != 0)  goto error;
//   status = poc_get_UG3D(ncid, mesh, frame, varid, buffer[0]);
//   if(status != 0)  goto error;
//   for(n = 0; n < udim; n++) {
//     state.u[n]=polar(buffer[0][n],-buffer[1][n]*M_PI/180.0);
//     }
//   status = nc_inq_varid(ncid, "a_v_LGP0", &varid);
//   if(status != 0)  goto error;
//   status = poc_get_UG3D(ncid, mesh, frame, varid, buffer[0]);
//   if(status != 0)  goto error;
//   status = nc_inq_varid(ncid, "G_v_LGP0", &varid);
//   if(status != 0)  goto error;
//   status = poc_get_UG3D(ncid, mesh, frame, varid, buffer[0]);
//   if(status != 0)  goto error;
//   for(n = 0; n < udim; n++) {
//     state.v[n]=polar(buffer[0][n],-buffer[1][n]*M_PI/180.0);
//     }
// 
//   return(0);
// 
// error:
//   return(-1);
// 
// }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int SpArchive_Get2D(const char *filename, int frame, cefmo_t cefmo, tide2D_t state)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------

  Spectral solution archive

 -----------------------------------------------------------------------*/
{
  int i,j,k,l,m,n,status;
  int ncid,varid;
  int zdim,udim;
  int z_discretisation, u_discretisation;
  int wave_id;
  double *buffer[2];
  char varname_a[1024], varname_G[1024];

  status=paire_discretisation_id(cefmo.paire, &z_discretisation, &u_discretisation);
  if(status!=0) return(-1);
  
  const char *UNAME=discretisation_name(u_discretisation);
  const char *ZNAME=discretisation_name(z_discretisation);

/**----------------------------------------------------------------------------
  elevation */
  sprintf(varname_a,"a_eta_%s", ZNAME);
  sprintf(varname_G,"G_eta_%s", ZNAME);
  status=poc_get_UG3D(filename, frame, (const char *) varname_a, (const char *) varname_G, state.z);
  if(status!=0) goto error;
  
/**----------------------------------------------------------------------------
  currents */
  sprintf(varname_a,"a_u_%s", UNAME);
  sprintf(varname_G,"G_u_%s", UNAME);
  status=poc_get_UG3D(filename, frame, (const char *) varname_a, (const char *) varname_G, state.u);
  if(status!=0) goto error;
  sprintf(varname_a,"a_v_%s", UNAME);
  sprintf(varname_G,"G_v_%s", UNAME);
  status=poc_get_UG3D(filename, frame, (const char *) varname_a, (const char *) varname_G, state.v);
  if(status!=0) goto error;

  return(0);

error:
  return(-1);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void SpArchive_Put2D(const mesh_t & mesh,const cefmo_t & cefmo,const parameter_t & data,const tidal_wave & wave,int iteration)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------

  Spectral solution archive

 -----------------------------------------------------------------------*/
{
  int i,j,k,l,m,n,status;
  int zdim,udim;
  int z_discretisation, u_discretisation;
//   int archives_level;
//   int archives_maxframe=-1;
  int frame;
  int wave_id;
  complex <double> *zLGP1,*uLGP1,*vLGP1;
  double *Zbuffer,*Ubuffer;
  char *comment[2],filename[1024],varname[1024];
  double mask=1.e+10;
  bool proceed;
  
  extern int archiving_SGtides2D(const mesh_t & mesh,const cefmo_t & cefmo,const tidal_wave & wave,int iteration,const Cgrid_t & Cgrid);

  status=paire_dimension(mesh,cefmo.paire,&zdim,&udim);
  status=paire_discretisation_id(cefmo.paire, &z_discretisation, &u_discretisation);
  
  const char *UNAME=discretisation_name(u_discretisation);
  const char *ZNAME=discretisation_name(z_discretisation);

  comment[0]=new char[128];
  comment[1]=new char[128];
  
/**-----------------------------------------------------------------------------
  Structured tidal constants archive */
  if(gHarmonic_StructuredArchive==1) {
    status=archiving_SGtides2D(mesh, cefmo, wave, iteration, gCgrid);
    }

//   if(tugo_cfg->tides->SpectralArchiveLevel.face_value()==(string) "MINIMAL") {
//     gArchivesLevel=CEFMO_ARCHIVES_MINIMAL;
//     archives_maxframe=0;
//     }
//   else if(tugo_cfg->tides->SpectralArchiveLevel.face_value()==(string) "STANDARD") {
//     gArchivesLevel=CEFMO_ARCHIVES_STANDARD;
//     }
//   else if(tugo_cfg->tides->SpectralArchiveLevel.face_value()==(string) "DEBUGGING") {
//     gArchivesLevel=CEFMO_ARCHIVES_DEBUGGING;
//     }
//   else if(tugo_cfg->tides->SpectralArchiveLevel.face_value()==(string) "EXPLORATION") {
//     gArchivesLevel=CEFMO_ARCHIVES_EXPLORATION;
//     }
//   else {
//     check_error(-1, "illegal setting spectral solution archive content", __LINE__, __FILE__, 1);
//     }

  proceed=(gArchivesLevel==CEFMO_ARCHIVES_STANDARD) || (gArchivesLevel==CEFMO_ARCHIVES_DEBUGGING);
  if(proceed) {
    sprintf(filename, "%s/%s.uv.%2.2d.v2c",gOutputPath, wave.name, iteration);
    sprintf(comment[0], "spectral solver, %s", tugo_version);
    sprintf(comment[1], "%s velocity (m/s)\n",UNAME);
    status=quoddy_savec2((const char*) filename, udim, cefmo.state.u, cefmo.state.v, (char **) comment);

    sprintf(filename, "%s/%s.ele.%2.2d.s2c",gOutputPath, wave.name, iteration);
    sprintf(comment[0], "spectral solver, %s", tugo_version);
    sprintf(comment[1], "%s elevation (m)\n", ZNAME);
    status=quoddy_savec1((const char*)filename, zdim, cefmo.state.z, (char **) comment);
#if 0
    zLGP1=new complex <double>[mesh.nvtxs];
    status=projection_LGP1(cefmo.state.z,zLGP1,mesh,z_discretisation);

    sprintf(filename, "%s/%s.ele.%2.2d.LGP1.s2c",gOutputPath, wave.name, iteration);
    sprintf(comment[0], "spectral solver, %s", tugo_version);
    sprintf(comment[1], "LGP1 elevation (m)\n");
    status=quoddy_savec1((const char*)filename, mesh.nvtxs, zLGP1, (char **) comment);

    uLGP1=new complex <double>[mesh.nvtxs];
    vLGP1=new complex <double>[mesh.nvtxs];
    status=projection_LGP1(cefmo.state.u,uLGP1,mesh,u_discretisation);
    status=projection_LGP1(cefmo.state.v,vLGP1,mesh,u_discretisation);

    sprintf(filename, "%s/%s.uv.%2.2d.LGP1.v2c",gOutputPath, wave.name, iteration);
    sprintf(comment[0], "spectral solver, %s", tugo_version);
    sprintf(comment[1], "LGP1 velocitiy (m/s)\n");
    status=quoddy_savec2((const char*)filename, mesh.nvtxs, uLGP1, vLGP1, (char **) comment);

    delete[] zLGP1;
    delete[] uLGP1;
    delete[] vLGP1;

#endif

  }
  
/**----------------------------------------------------------------------------
  output path will be added in archiving_UGdummy2D routine*/
  sprintf(filename, "%s.spectral.nc", wave.name);

  if(gArchivesMaxFrame!=-1) frame=MIN(gArchivesMaxFrame, iteration);
  else frame=iteration;
  
  proceed=(gArchivesLevel==CEFMO_ARCHIVES_STANDARD) || (gArchivesLevel==CEFMO_ARCHIVES_DEBUGGING);
  if( (iteration==0) && (proceed)) {
    sprintf(varname,"h_%s", UNAME);
    status=archiving_UGdummy2D(filename, mesh, varname, "m", cefmo.state.h_u, mask, NOFRAME, u_discretisation);
    sprintf(varname,"dhdx_%s", UNAME);
    status=archiving_UGdummy2D(filename, mesh, varname, "none", data.dhdx, mask, NOFRAME, u_discretisation);
    sprintf(varname,"dhdy_%s", UNAME);
    status=archiving_UGdummy2D(filename, mesh, varname, "none", data.dhdy, mask,  NOFRAME, u_discretisation);
    sprintf(varname,"h_%s", ZNAME);
    status=archiving_UGdummy2D(filename, mesh, varname, "m", cefmo.state.h_z, mask, NOFRAME, z_discretisation);
    status=archiving_UGdummy2D(filename, mesh, "Z0", "m",    data.z0, mask, u_discretisation);
    status=archiving_UGdummy2D(filename, mesh, "u0", "m/s",  data.u0, mask, u_discretisation);
//    status=archiving_UGdummy2D(filename, mesh, "Cd", "none", data.Cd, mask, frame, u_discretisation);
    status=archiving_UGdummy2D(filename, mesh, "Cd", "none", data.Cd, mask, u_discretisation);
    }
  
  Ubuffer=new double [udim];
  for(n = 0; n < udim; n++) {
    Ubuffer[n]=cefmo.u_descriptor->nodes[n].g;
    }
  sprintf(varname,"g_%s", UNAME);
  status=archiving_UGdummy2D(filename, mesh, varname, "ms-2", Ubuffer, mask, NOFRAME, u_discretisation);
  delete[] Ubuffer;
   
  Zbuffer=new double [zdim];

/**----------------------------------------------------------------------------
  elevation */
  for(n = 0; n < zdim; n++) {
    Zbuffer[n]=abs(cefmo.state.z[n]);
    }
  sprintf(varname,"a_eta_%s", ZNAME);
  status=archiving_UGdummy2D(filename, mesh, varname, "m", Zbuffer, mask, frame,z_discretisation);
  for(n = 0; n < zdim; n++) {
    Zbuffer[n]=-arg(cefmo.state.z[n])*180/M_PI;
    }
  sprintf(varname,"G_eta_%s", ZNAME);
  status=archiving_UGdummy2D(filename, mesh, varname, "degrees", Zbuffer, mask, frame,z_discretisation);

  proceed=(gArchivesLevel==CEFMO_ARCHIVES_STANDARD) || (gArchivesLevel==CEFMO_ARCHIVES_DEBUGGING);
  if(proceed) {
    frame=NOFRAME;
/**----------------------------------------------------------------------------
    potential */
    for(n = 0; n < zdim; n++) {
      Zbuffer[n]=abs(cefmo.state.potential[n]);
      }
    sprintf(varname,"a_POT_%s", ZNAME);
    status=archiving_UGdummy2D(filename, mesh, varname, "m", Zbuffer, mask, frame,z_discretisation);
    for(n = 0; n < zdim; n++) {
      Zbuffer[n]=-arg(cefmo.state.potential[n])*180/M_PI;
      }
    sprintf(varname,"G_POT_%s", ZNAME);
    status=archiving_UGdummy2D(filename, mesh, varname, "degrees", Zbuffer, mask, frame,z_discretisation);

/**----------------------------------------------------------------------------
    loading/self-attraction */
    for(n = 0; n < zdim; n++) {
      Zbuffer[n]=abs(cefmo.state.LSA[n]);
      }
    sprintf(varname,"a_LSA_%s", ZNAME);
    status=archiving_UGdummy2D(filename, mesh, varname, "m", Zbuffer, mask,frame,z_discretisation);
    for(n = 0; n < zdim; n++) {
      Zbuffer[n]=-arg(cefmo.state.LSA[n])*180/M_PI;
      }
    sprintf(varname,"G_LSA_%s", ZNAME);
    status=archiving_UGdummy2D(filename, mesh, varname, "degrees", Zbuffer, mask,frame,z_discretisation);
/**----------------------------------------------------------------------------
    Cd ratio */
//    status=archiving_UGdummy2D(filename, mesh, "Cd_ratio", CdRatio,frame,u_discretisation);
    }
    
/**----------------------------------------------------------------------------
  currents */
  if(gArchivesMaxFrame!=-1) frame=MIN(gArchivesMaxFrame, iteration);
  else frame=iteration;
  
  proceed=(gArchivesLevel==CEFMO_ARCHIVES_MINIMAL) || (gArchivesLevel==CEFMO_ARCHIVES_STANDARD) || (gArchivesLevel==CEFMO_ARCHIVES_DEBUGGING);
  if(proceed) {
    Ubuffer=new double [udim];

    for(n = 0; n < udim; n++) {
      Ubuffer[n]=abs(cefmo.state.u[n]);
      }
    sprintf(varname,"a_u_%s", UNAME);
    status=archiving_UGdummy2D(filename, mesh, varname, "m/s", Ubuffer, mask, frame,u_discretisation);
    for(n = 0; n < udim; n++) {
      Ubuffer[n]=-arg(cefmo.state.u[n])*180/M_PI;
      }
    sprintf(varname,"G_u_%s", UNAME);
    status=archiving_UGdummy2D(filename, mesh, varname, "degrees", Ubuffer, mask, frame,u_discretisation);
    for(n = 0; n < udim; n++) {
      Ubuffer[n]=abs(cefmo.state.v[n]);
      }
    sprintf(varname,"a_v_%s", UNAME);
    status=archiving_UGdummy2D(filename, mesh, varname, "m/s", Ubuffer, mask, frame,u_discretisation);
    for(n = 0; n < udim; n++) {
      Ubuffer[n]=-arg(cefmo.state.v[n])*180/M_PI;
      }
    sprintf(varname,"G_v_%s", UNAME);
    status=archiving_UGdummy2D(filename, mesh, varname, "degrees", Ubuffer, mask, frame,u_discretisation);
    delete[] Ubuffer;
    }
    
/**-----------------------------------------------------------------------
  ustar */
  proceed=(gArchivesLevel==CEFMO_ARCHIVES_STANDARD) || (gArchivesLevel==CEFMO_ARCHIVES_DEBUGGING);
  if(proceed) {
    complex<double> *ustar=new complex<double>[udim];
    complex<double> *vstar=new complex<double>[udim];
    complex<double> cmask=complex<double>(mask,mask);
    for(n = 0; n < udim; n++) {
      complex<double> uu,vv;
      double a, b, pol, dir, time;
      double CdRoot=sqrt(data.Cd[n]);
      uu=cefmo.FrictionMatrix[n].c[0][0]*cefmo.state.u[n]+cefmo.FrictionMatrix[n].c[1][0]*cefmo.state.v[n];
      vv=cefmo.FrictionMatrix[n].c[0][1]*cefmo.state.u[n]+cefmo.FrictionMatrix[n].c[1][1]*cefmo.state.v[n];
//      ellipse_parameter(uu, vv, &a, &b, &pol, &dir, &time);
      ustar[n]=uu*CdRoot;
      vstar[n]=vv*CdRoot;
      }
  
    status=archiving_UGdummy2D(filename, mesh, "a_uStar", "G_uStar", "m/s", ustar, cmask, frame, u_discretisation);
    status=archiving_UGdummy2D(filename, mesh, "a_vStar", "G_vStar", "m/s", vstar, cmask, frame, u_discretisation);
    delete[] ustar;
    delete[] vstar;
    }
  
  delete[] Zbuffer;

  delete[] comment[0];
  delete[] comment[1];
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void SpMomentum2D_solve(mesh_t & mesh, cefmo_t & model, complex<double> *rhsU,complex<double> *rhsV)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i, j, k, l, m, n;
  int udim,zdim;
  int status;

  status=paire_dimension(mesh,model.paire,&zdim,&udim);

  for(i = 0; i < udim; i++) {
    model.state.u[i] = rhsU[i]*model.InvSpMu[i][0] + rhsV[i]*model.InvSpMu[i][1];
    model.state.v[i] = rhsU[i]*model.InvSpMv[i][0] + rhsV[i]*model.InvSpMv[i][1];
    }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int SpWE2D_Terminate(mesh_t mesh,cefmo_t model)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------

  delete arrays and structures

 -----------------------------------------------------------------------*/
{
  int n,status;
  int udim,zdim;

/**-----------------------------------------------------------------------
  */
  status=paire_dimension(mesh,model.paire,&zdim,&udim);

  if(model.InvSpMu!=0) {
    for(n=0;n<udim;n++) {
      delete[] model.InvSpMu[n];
      }
    delete[] model.InvSpMu;
    }
  if(model.InvSpMv!=0) {
    for(n=0;n<udim;n++) {
      delete[] model.InvSpMv[n];
      }
    delete[] model.InvSpMv;
    }
  model.InvSpMu=0;
  model.InvSpMv=0;

  free_ordering(&(model.M1.ordering));

  delete[] model.M1.packed;
  delete[] model.M2.packed;

  delete[] model.G1.packed;
  delete[] model.G2.packed;

  model.M1.packed=0;
  model.M2.packed=0;

  model.G1.packed=0;
  model.G2.packed=0;

  free_ordering(&(model.D1.ordering));

  delete[] model.D1.packed;
  delete[] model.D2.packed;

  model.D1.packed=0;
  model.D2.packed=0;

  delete[] model.SpMatrix.packed;
  model.SpMatrix.packed=0;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> zconvert(hconstant_t x, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------------
  polar to complex conversion (tidal convention on phase)*/
  double c=1.0,R1,R2;
  complex<double> z;
  double a,G;

  a=(double) x.a[k];
  G=(double) x.G[k];

  z=polar(a,-G);
  return(z);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<float> cconvert(hconstant_t x, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------------
  polar to complex conversion (tidal convention on phase)*/
  double c=1.0,R1,R2;
  complex<float> z;
  double a,G;

  a=(float) x.a[k];
  G=(float) x.G[k];

  z=polar(a,-G);
  return(z);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void tidal_equilibrium_harmonics(tidal_wave w, parameter_t data,int nnodes, complex<double> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, k;
  double a, G, V, dV, V0, C, S, x, y;
  complex< double > equilibrium;
  double k2 = love_k2, h2 = love_h2;
 
  double omega;
  
  complex< double > anelastic;
  
  switch(gAnelasticBodyTide) {
    case 0:
      anelastic=polar(1.0,0.0);
    case 1:
      anelastic=polar(1.0,0.2*M_PI/180.);
    }

/*----------------------------------------------------------------------
  Compute astronomical potential constants (amplitude and phase lag)
----------------------------------------------------------------------*/

  a=w.Ap;
  
/*------------------------------------------------------------------------------
  designed for non tidal frequency simulation*/
  if(a==-1) {
    for(n = 0; n < nnodes; n++) {
      z[n]=0.0;
      }
    return;
    }
  
  equilibrium = a * (1. + k2*anelastic - h2*anelastic);
  switch (w.nT) {
    case (0):

/*######################### Long period tide #########################

   potential/g = A*(1/2 -3/2 sin^2(Y)*cos(w*t+V0)

   dP/dx=  0
   dP/dy= -3*A*cos(Y)sin(Y)*cos(w*t+V0)

----------------------------------------------------------------------*/

      for(n = 0; n < nnodes; n++) {
        dV = 0.0;
        C = data.C[n];
        S = data.S[n];
        a = (0.5 - 1.5 * S * S);
        G = -dV;
        z[n]=equilibrium * polar(a,-G);
        }
      break;

    case (1):

/*########################### Diurnal tide ###########################

   potential/g = A*sin(2Y)*cos(w*t+V0+X)

   dP/dx=  -2*A*cos(Y)*sin(Y)*sin(w*t+V0+X)
   dP/dy=   2*A*cos(2Y)*cos(w*t+V0+X) = 2*A*[cos^2(Y)-sin^2(Y)]*cos(w*t+V0+X)

----------------------------------------------------------------------*/

      for(n = 0; n < nnodes; n++) {
        dV = data.lon[n];
        C = data.C[n];
        S = data.S[n];
        a = 2 * S * C;
        G = -dV;
        z[n]=equilibrium * polar(a,-G);
        }
      break;

    case (2):

/*######################### Semi-diurnal tide #########################

   potential/g = A*cos^2(Y)*cos(w*t+V0+2*X)

   dP/dx=  -2*A*cos^2(Y)*sin(w*t+V0+2*X)
   dP/dy=  -2*A*cos(Y)*sin(Y)*cos(w*t+V0+2*X)

----------------------------------------------------------------------*/

      for(n = 0; n < nnodes; n++) {
        dV = 2 * data.lon[n];
        C = data.C[n];
        S = data.S[n];
        a = C * C;
        G = -dV;
        z[n]=equilibrium * polar(a,-G);
        }
      break;

/*####################### non-astronomical tide #######################
   potential/g = 0

----------------------------------------------------------------------*/
    default:
      for(n = 0; n < nnodes; n++) {
        a = 0;
        G = 0;
        z[n]=polar(a,-G);
        }
      break;

  }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cefmo_topography_T(mesh_t & mesh, tide2D_t & state, parameter_t & data, int paire)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  int z_discretisation, u_discretisation;
  int udim,zdim;
  bool specific;
  string s="";

  cout << "\n-----------------------------------------------------------------------------" << endl;
  cout << "ocean bottom topography setup" << endl;
  status=paire_discretisation_id(paire, &z_discretisation, &u_discretisation);
  status=paire_dimension(mesh,paire,&zdim,&udim);
  
/**-----------------------------------------------------------------------
  topography */
  specific=(tugo_cfg->tides->SpecificData.face_value()==(string) "TRUE");
  if(specific) {
    s=tugo_cfg->tides->TopoFile.face_value();
    if(s==(string) "NONE") {
      check_error(-1, "discretisation-specific database ON, but topography not set", __LINE__, __FILE__, 0);
      }
    }
  else {
    s=(string) "NONE";
    }
    
  if(s!=(string) "NONE") {
    printf("cefmo_topography, loading %s discretisation-specific topography: %s\n",discretisation_name(u_discretisation), s.c_str());
    status=quoddy_loadr1((const char*) s.c_str(), udim, state.h_u);
    if(status!=0) {
      check_error(-1, "error in topography file reading", __LINE__, __FILE__, 1);
      }
    for(n = 0; n < udim; n++) {
      state.h_u[n]=MAX(state.h_u[n],10.0);
      }
    status=fe_projection(mesh, state.h_u, u_discretisation, state.h_z, z_discretisation);
    }
  else {
    printf("cefmo_topography, interpolating LGP1 topography\n");
    status=fe_projection(mesh, gP1data2D.h,  LGP1,  state.h_z, z_discretisation);
    status=fe_projection(mesh, gP1data2D.h,  LGP1,  state.h_u, u_discretisation);
    }

/**-----------------------------------------------------------------------
  topography gradient needed for internal tide drag*/
  
  if(specific) {
    s=tugo_cfg->tides->SlopeFile.face_value();
    if(s==(string) "NONE") {
      check_error(-1, "discretisation-specific database ON, but topography slope not set", __LINE__, __FILE__, 0);
      }
    }
  else {
    s=(string) "NONE";
    }
  if(s!=(string) "NONE") {
    printf("cefmo_topography, loading %s discretisation-specific topography slope: %s\n",discretisation_name(u_discretisation),s.c_str());
    status=quoddy_loadr2((const char*) s.c_str(), udim, data.dhdx, data.dhdy);
    if(status!=0) {
      check_error(-1, "error in topography file reading", __LINE__, __FILE__, 1);
      }

/**-----------------------------------------------------------------------------
    slope is given in m/km, must be divide by 1000. to be really dimensionless*/
    for(n = 0; n < udim; n++) {
      data.dhdx[n] = data.dhdx[n]/1000.;
      data.dhdy[n] = data.dhdy[n]/1000.;
      }
    }
  else {
    printf("cefmo_topography, interpolating LGP1 topography slope\n");
//     status=smoother2D_LGP1xLGP0 (mesh, gP1data2D.dhdx);/// HERE 
    status=fe_projection(mesh, gP1data2D.dhdx, LGP1, data.dhdx, u_discretisation);
//     status=smoother2D_LGP1xLGP0 (mesh, gP1data2D.dhdy);/// HERE 
    status=fe_projection(mesh, gP1data2D.dhdy, LGP1, data.dhdy, u_discretisation);
//    exit(-1);
    }
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cefmo_topography_Q(mesh_t & mesh, tide2D_t & state, parameter_t & data, int paire)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  int z_discretisation, u_discretisation;
  int udim,zdim;
  bool specific;
  string s="";

  cout << "\n-----------------------------------------------------------------------------" << endl;
  cout << "ocean bottom topography setup" << endl;
  status=paire_discretisation_id(paire, &z_discretisation, &u_discretisation);
  status=paire_dimension(mesh,paire,&zdim,&udim);
  
/**-----------------------------------------------------------------------
  topography */
  specific=(tugo_cfg->tides->SpecificData.face_value()==(string) "TRUE");
  if(specific) {
    s=tugo_cfg->tides->TopoFile.face_value();
    if(s==(string) "NONE") {
      check_error(-1, "discretisation-specific database ON, but topography not set", __LINE__, __FILE__, 0);
      }
    }
  else {
    s=(string) "NONE";
    }
    
  if(s!=(string) "NONE") {
    printf("cefmo_topography, loading %s discretisation-specific topography: %s\n",discretisation_name(u_discretisation), s.c_str());
    status=quoddy_loadr1((const char*) s.c_str(), udim, state.h_u);
    if(status!=0) {
      check_error(-1, "error in topography file reading", __LINE__, __FILE__, 1);
      }
    for(n = 0; n < udim; n++) {
      state.h_u[n]=MAX(state.h_u[n],10.0);
      }
    status=fe_projection(mesh, state.h_u, u_discretisation, state.h_z, z_discretisation);
    }
  else {
    printf("cefmo_topography, interpolating CQP1 topography\n");
    status=fe_projection(mesh, gP1data2D.h,  CQP1,  state.h_z, z_discretisation);
    status=fe_projection(mesh, gP1data2D.h,  CQP1,  state.h_u, u_discretisation);
    }

/**-----------------------------------------------------------------------
  topography gradient needed for internal tide drag*/
  
  if(specific) {
    s=tugo_cfg->tides->SlopeFile.face_value();
    if(s==(string) "NONE") {
      check_error(-1, "discretisation-specific database ON, but topography slope not set", __LINE__, __FILE__, 0);
      }
    }
  else {
    s=(string) "NONE";
    }
  if(s!=(string) "NONE") {
    printf("cefmo_topography, loading %s discretisation-specific topography slope: %s\n",discretisation_name(u_discretisation),s.c_str());
    status=quoddy_loadr2((const char*) s.c_str(), udim, data.dhdx, data.dhdy);
    if(status!=0) {
      check_error(-1, "error in topography file reading", __LINE__, __FILE__, 1);
      }
/**-----------------------------------------------------------------------------
    slope is given in m/km, must be divide by 1000. to be really dimensionless*/
    for(n = 0; n < udim; n++) {
      data.dhdx[n] = data.dhdx[n]/1000.;
      data.dhdy[n] = data.dhdy[n]/1000.;
      }
    }
  else {
    printf("cefmo_topography, interpolating CQP1 topography slope\n");
    status=fe_projection(mesh, gP1data2D.dhdx, CQP1, data.dhdx, u_discretisation);
    status=fe_projection(mesh, gP1data2D.dhdy, CQP1, data.dhdy, u_discretisation);
    }
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cefmo_topography(mesh_t & mesh, tide2D_t & state, parameter_t & data, int paire)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  switch (mesh.nature()) {
    case FE_TRIANGLE:
      status=cefmo_topography_T(mesh,  state,  data, paire);
      break;
    case FE_QUADRANGLE:
      status=cefmo_topography_Q(mesh,  state,  data, paire);
      break;
    }
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cefmo_solver(mesh_t *mesh, parameter_t data)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------

  Spectral tidal solver

 -----------------------------------------------------------------------*/
{
  int i,j,k,l,m,n,status;
  cefmo_t cefmo;
  extern parameter_t parameters_initialise2D(mesh_t,int);
  int zdim,udim;
  parameter_t SpParameters;
  int z_discretisation, u_discretisation;
  int recycle=1;

  cout << "\n-----------------------------------------------------------------------------" << endl;
  cout << "Spectral tidal solver: start" << endl;

  string PaireName=tugo_cfg->tides->SpectralPaire.face_value();
  cefmo.paire=discretisation_PaireIdFromName(PaireName.c_str());
    
/**-----------------------------------------------------------------------------
  identify dominant waves */
  {
  char *name;
  int id;
  name=strdup(tugo_cfg->tides->DominantWave_1.face_value().c_str());
  cefmo.dominants[0]=WaveList.wave((const char *) name);
  free(name);
  name=strdup(tugo_cfg->tides->DominantWave_2.face_value().c_str());
  cefmo.dominants[1]=WaveList.wave((const char *) name);
  free(name);
  cefmo.ndominants=0;
  id=WaveList.wave_index(cefmo.dominants[0]);
  if(id!=-1) cefmo.ndominants=1;
  id=WaveList.wave_index(cefmo.dominants[1]);
  if(id!=-1) cefmo.ndominants=2;
  }
  
/**-----------------------------------------------------------------------------
  init paire and corresponding discretisation descriptors */
  status=paire_discretisation_id(cefmo.paire, &z_discretisation, &u_discretisation);
  
#ifdef PARALLEL
//#ifdef HAVE_MPI
  int context=PARALLEL_COMPUTING;
#else
  int context=SEQUENTIAL_COMPUTING;
#endif

  status= discretisation_init(mesh, u_discretisation, context);
  status= discretisation_init(mesh, z_discretisation, context);
  
  status=paire_dimension(*mesh,cefmo.paire,&zdim,&udim);
  
/**-----------------------------------------------------------------------------
  allocate memory for and initialize model parameters */
  SpParameters=parameters_initialise2D(*mesh,cefmo.paire);
  
  delete[] SpParameters.dhdx;
  delete[] SpParameters.dhdy;
  SpParameters.dhdx=new double[udim];
  SpParameters.dhdy=new double[udim];
  
/**-----------------------------------------------------------------------------
  initialize loading/self-attraction */
  if(LSA_forcing)
    {/* STS: NO LSA */}
  
  cefmo.u_descriptor=get_descriptor_address(*mesh, u_discretisation);
  cefmo.z_descriptor=get_descriptor_address(*mesh, z_discretisation);
  
/**-----------------------------------------------------------------------------
  initialize nodes gravitational constant */
  status=gravitational_constant(*cefmo.u_descriptor);
  status=gravitational_constant(*cefmo.z_descriptor);

  if(tugo_cfg->tides->SpectralArchiveLevel.face_value()==(string) "MINIMAL") {
    gArchivesLevel=CEFMO_ARCHIVES_MINIMAL;
    gArchivesMaxFrame=0;
    }
  else if(tugo_cfg->tides->SpectralArchiveLevel.face_value()==(string) "STANDARD") {
    gArchivesLevel=CEFMO_ARCHIVES_STANDARD;
    }
  else if(tugo_cfg->tides->SpectralArchiveLevel.face_value()==(string) "DEBUGGING") {
    gArchivesLevel=CEFMO_ARCHIVES_DEBUGGING;
    }
  else if(tugo_cfg->tides->SpectralArchiveLevel.face_value()==(string) "EXPLORATION") {
    gArchivesLevel=CEFMO_ARCHIVES_EXPLORATION;
    }
  else {
    check_error(-1, "illegal setting spectral solution archive content", __LINE__, __FILE__, 1);
    }
  
  switch(cefmo.paire) {
    case CQN1xCQP0:
/**----------------------------------------------------------------------------
      two modes to be prepared, 1 FV the other FE */
/**----------------------------------------------------------------------------
      U node x Z node connexion (for pressure gradient) */
      cefmo.packed_UBASExZBASE=init_packed_CQN1xCQP0;
/**----------------------------------------------------------------------------
      Z node x Z node connexion (wave equation), now implicit */
      cefmo.packed_ZBASExZBASE=init_packed_CQP0xCQP0_01;
/**----------------------------------------------------------------------------
      Z node x U node connexion (for transport divergence) */
      cefmo.packed_ZBASExUBASE=init_packed_CQP0xCQN1;
/**----------------------------------------------------------------------------
      pressure gradient matrix */
      cefmo.SpGradient2D=SpGradient2D_CQN1xCQP0;
/**----------------------------------------------------------------------------
      transport divergence matrix */
      cefmo.SpDivergence2D=SpDivergence2D_CQN1xCQP0;
/**----------------------------------------------------------------------------
      non-linear divergence rhs */
      cefmo.SpDivergence2D_RHS=WE2D_CQN1xCQP0_SpDivergence;
/**----------------------------------------------------------------------------
      non-linear advection rhs */
      cefmo.SpAdvection2D=momentum2d_CQN1xCQP0_SpAdvection;
      
//       status= cefmo2D(*mesh, SpParameters, cefmo, &cefmo_atlas2D, 0);
      if(cefmo.ndominants!=0) {
        status= cefmo2D(*mesh, SpParameters, cefmo, &cefmo_atlas2D, 0);
        }
      else {
        /* STS: NO experimental */
        }
      if(Spectral3DRun!=UNDEFINED_MODE) {
        check_error(-1, "implementation not finalized...", __LINE__, __FILE__, 1);
//        status= cefmo3D_LGP0xLGP1(*mesh, SpParameters, cefmo, &cefmo_atlas3D, cefmo_atlas2D);
        }
      break;
      
    default:
      check_error(-1, "illegal spectral sover discretisation...", __LINE__, __FILE__, 1);
      break;
    }

  if(status!=0) {
    check_error(-1, "spectral sover failed...", __LINE__, __FILE__, 1);
    }
  cout << "Spectral tidal solver: OK" << endl;
  
  /*STS END OF PROGRAMME*/
  exit(0);

  return(status);
}

