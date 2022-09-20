
/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/**
\author  Florent Lyard      LEGOS/CNRS, Toulouse, France
  E-mail: florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada
*/

#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>

#include "archive.h"

#include "poc-netcdf.hpp"
#include "solverlib.h"

#include "tides.h"
#include "fe.h"
#include "fe-integration.h"
#include "map.h"
#include "mass.h"
#include "cefmo.h"

#include "filter.h"
#include "statistic.h"
#include "functions.h"
#include "maths.h"
#include "zapper.h"

#include "matrix.h"

#include "constants.h"

using namespace std;


double **grad_opx,**grad_opy;
int idx = 0;

extern matrix2x2_t  *DragMatrix;
extern zmatrix2x2_t *FricMatrix;
extern void spectral_friction01(mesh_t, parameter_t, int, int);
extern void spectral_friction02(mesh_t, parameter_t, int, int);
extern void spectral_friction03(spectrum_t spectrum,double duration, parameter_t data, int nndes, int target, double **buffer);
extern void spectral_viscosity(spectrum_t ,double ,mesh_t mesh,parameter_t , int , double **, double **);  


extern void spectral_transport(mesh_t & mesh, spectrum_t & AnalysisList, parameter_t & data2D, const char *output_path, char *extension);
extern void spectral_energy01(mesh_t & mesh,spectrum_t & AnalysisList, parameter_t & data2D, action_t *P1_action, const char *output_path, const char *path=NULL);
extern int  spectral_EnergyBudget(const char *energy, const char *polygonsfile, spectrum_t spectrum, bool debug);
extern void spectral_depth(mesh_t mesh,spectrum_t AnalysisList, parameter_t data2D, char *output_path, char *path);

extern  int spectral_EnergyComputation(const char *meshfile, const char *discretisation, const char * path, const char *output_path, spectrum_t & spectrum, spectrum_t & dominants);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void gradient_operator(mesh_t mesh,parameter_t data2D)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n;
  int dim;
  int nndes;
  int status;
  float C, hC, hc;
  double S, dP, dQ, dZdx, dZdy;

  nndes = mesh.nvtxs;
  grad_opx = new double*[nndes];
  grad_opy = new double*[nndes];
  for(n = 0; n < nndes; n++) {
    grad_opx[n] = new double[mesh.vertices[n].nngh+1];
    grad_opy[n] = new double[mesh.vertices[n].nngh+1];
    }
/*----------------------------------------------------------------------
  gradient operator initialisation  */
  for(n = 0; n < nndes; n++) {
    for(i = 0; i <= mesh.vertices[n].nngh; i++) {
      grad_opx[n][i] = 0;
      grad_opy[n][i] = 0;
      }
    }

/*-----------------------------------------------------------------------
  compute contribution of triangles weighted by triangle area */
  for(l = 0; l < mesh.ntriangles; l++) {
    for(j = 0; j < 3; j++) {
      n = mesh.triangles[l].vertex[j];
      C = data2D.C[n];
      dP = mesh.triangles[l].DP[j] / 2.0e0; //  no /A, == multiply then divide
      dQ = mesh.triangles[l].DQ[j] / 2.0e0;
      dZdx = +dQ / C;
      dZdy = -dP;
      for(i = 0; i < 3; i++) {
        m = mesh.triangles[l].vertex[i];
        S = 3. * data2D.lmm[m];   /* total surface weight */
        if(m == n) {
          grad_opx[m][0] += dZdx / S;
          grad_opy[m][0] += dZdy / S;
          }
        else {
          for(k = 0; k < mesh.vertices[m].nngh; k++) {
            if(mesh.vertices[m].ngh[k] == n) {
              grad_opx[m][k + 1] += dZdx / S;
              grad_opy[m][k + 1] += dZdy / S;
              break;
              }
            }
          }
        }
      }
    }
}




/*----------------------------------------------------------------------------*/
/// Gives the gradient of LGP1 Beta for node n at cornered triangle j
/**
\date reviewed 28 Apr 2011
\author Damien Allain

\param mesh grid
\param n node index
\param j cornered triangle index
\param *C longitude gradient component
\param *D latitude gradient component
\returns index of the triangle
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_LGP1_beta_gradient(mesh_t *mesh, int n, int j, double *C, double *D)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Development notes on B(eta) gradient components C and D :

C x0 + D y0 = B0=1
C x1 + D y1 = B1=0  =>  C (x0-x1) + D (y0-y1) = B0-B1=1
C x2 + D y2 = B2=0      C (x0-x2) + D (y0-y2) = B0-B2=1

From http://en.wikipedia.org/wiki/Cramer%27s_rule :
            Det= (x0-x1)(y0-y2) - (y0-y1)(x0-x2)
     =>     C=((y0-y2)-(y0-y1))/Det  = (y1-y2)/Det
            D=((x0-x1)-(x0-x2))/Det  = -(x1-x2)/Det

From fe01.cpp:fe_initaffine_spherical, this gives :
  Det= TrueArea*2 = Area*cos(lat)*2
  C=-DQ0/(TrueArea*2)
  D=-DP0*cos(lat)/(Area*cos(lat)*2)=-DP0/(Area*2)

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  int m, i;//indexes of triangle and its corners
  triangle_t *tp;//pointer to triangle
  m=mesh->vertices[n].elmts[j];
  tp = &mesh->triangles[m];
  for(i = 0; i < 3; i++) {// for all corners
    if(tp->vertex[i]==n){
      *C=tp->DQ[i];
      *D=-tp->DP[i];
      return m;
      }
    }
  fprintf(stderr,"line %d of %s: WRONG MESH\n",__LINE__,__FILE__);
  return -1;//should not reach here
}


// /*----------------------------------------------------------------------------*/
// /// Gives the gradient of NCP1 Beta for node n at cornered triangle j
// /**
// \date created 28 Apr 2011
// \author Damien Allain
// 
// \param mesh grid
// \param n node index
// \param j cornered triangle index
// \param *C longitude gradient component
// \param *D latitude gradient component
// \returns index of the triangle
// */
// /*----------------------------------------------------------------------------*/
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int fe_NCP1_beta_gradient(mesh_t *mesh, int n, int j, double *C, double *D)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
// /* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// 
// Development notes
// 
// See fe_LGP1_beta_gradient notes.
// 
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
//   int m, i;//indexes of triangle and its corners
//   triangle_t *tp;//pointer to triangle
//   m=mesh->vertices[n].elmts[j];
//   tp = &mesh->triangles[m];
//   for(i = 0; i < 3; i++) {// for all corners
//     if(tp->vertex[i]==n){
//       *C=-2*tp->DQ[i];
//       *D=2*tp->DP[i];
//       return m;
//       }
//     }
//   fprintf(stderr,"line %d of %s: WRONG MESH\n",__LINE__,__FILE__);
//   return -1;//should not reach here
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int decode_loadingname(char *wave, char **filename,const char *loading_directory,const char *loading_convention)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
  Atmospheric pressure and wind from AM fields; load the field at or
  just before t
 ----------------------------------------------------------------------*/
{
  int k, i, nt, nframe, status, hour;
  char dummy[256], *pointer;
  char *tmp,*p2;
  FILE *out;

  size_t *start, *count;
  cdfvar_t info;

/*----------------------------------------------------------------------
  build the atlas file name*/
  (*filename) = new char[strlen(loading_directory) + 256];
  sprintf((*filename), "%s/%s", loading_directory, loading_convention);

/*----------------------------------------------------------------------
  check existence name*/
  out = fopen(*filename, "r");
  status = (out == NULL);

  switch (status) {
    case 0:
/*----------------------------------------------------------------------
      file exists, do nothing more*/
      fclose(out);
      break;

    default:
/*----------------------------------------------------------------------
      use format information*/
      pointer = strstr((*filename), "WAVE");
      if(pointer != NULL) {
        tmp=strdup(*filename);
        p2 = strstr(tmp, "WAVE");
        sprintf(dummy, "%s", wave);
        strcpy(pointer, dummy);
        strcat(pointer, p2+4);
        status=0;
        }
      pointer = strstr((*filename), "wave");
      if(pointer != NULL) {
        tmp=strdup(*filename);
        p2 = strstr(tmp, "wave");
        sprintf(dummy, "%s", wave);
        int k;
        for(k=0;k<strlen(dummy);k++) {
          dummy[k]=tolower(dummy[k]);
          }
        strcpy(pointer, dummy);
        strcat(pointer, p2+4);
        status=0;
        }
    }

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int read_loading02(spectrum_t WaveList, parameter_t data2D, int nndes, action_t *P1_action, int *found)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
 Loading and self-attraction from diagnostic fields
 Read the amplitude (meter) and phase lag, convert into
 real and imaginary part
----------------------------------------------------------------------*/
{
  int i, k, m, n, items, status;
  float zr, zi, a, G;
  float *amp, *pha, *zx, *zy, factor, mask;
  fcomplex *z,cmask;
  float dummy;
  double x, y, dum;
  char *filename=0;
  size_t count[3], start[3], var;
  size_t length;
  int ncid, amp_id, pha_id;
  variable_t varinfo;

  grid_t LSAgrid;

//  filename = new char[1024];

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = 1;
  count[2] = 1;

  for(k = 0; k < WaveList.n; k++) {
//    status=decode_loadingname(WaveList.waves[k].name,&filename,"/home/softs/data/loading/FES2004/","WAVE.load.nc");
    status=decode_loadingname(WaveList.waves[k].name,&filename,"/home/softs/data/loading/FES2004/","WAVE.load.nc");
    printf("trying %s\n",filename);
    status = cdf_loadvargrid_2d(filename, 1, &LSAgrid);
    if(status != 0) {
      delete[] filename;
      continue;
      }
    count[1] = LSAgrid.ny;
    count[2] = LSAgrid.nx;
    amp = new float[LSAgrid.nx * LSAgrid.ny];
    pha = new float[LSAgrid.nx * LSAgrid.ny];

    status = nc_open(filename, NC_NOWRITE, &ncid);
    status = nc_inq_varid(ncid, "Ha", &amp_id);
    if(status != 0)
      amp_id = -1;
    status = nc_inq_varid(ncid, "Hg", &pha_id);
    if(status != 0)
      pha_id = -1;
    nc_close(ncid);

    printf("file %s, amplitude id:%d phase id:%d\n",filename, amp_id, pha_id);
//    status = cdf_loadvar_r1(filename, amp_id, start, count, amp, &mask, units);
    status= cdf_loadvar_r1_2d (filename, amp_id, 0, 0, LSAgrid, LSAgrid.nx, amp, &mask ,&varinfo);
    if(strcmp(varinfo.units, "cm") == 0)
      factor = 1.0e-02;
    else if(strcmp(varinfo.units, "m") == 0)
      factor = 1.0;
    else {
      printf("no units given in %s, assume loading in meters\n", filename);
      factor = 1.0;
      }

//    status = cdf_loadvar_r1(filename, pha_id, start, count, pha, &mask,units);
    status= cdf_loadvar_r1_2d (filename, pha_id, 0, 0, LSAgrid, LSAgrid.nx, pha, &mask ,&varinfo);
    if(status != 0)  goto error;

    cmask=fcomplex(mask,mask);
    z = new fcomplex[LSAgrid.nx * LSAgrid.ny];
    for(m = 0; m < LSAgrid.ny * LSAgrid.nx; m++) {
      if((amp[m] != mask) && (pha[m] != mask)) {
        z[m] = polar<float>(factor * amp[m],-pha[m] *d2r);
        }
      else {
        z[m] = cmask;
        }
      }
    for(n = 0; n < nndes; n++) {
      x = data2D.lon[n] *r2d;
      y = data2D.lat[n] *r2d;
      if(x < LSAgrid.xmin - LSAgrid.dx / 2.)
        x = x + 360.0;
      if(x > LSAgrid.xmax + LSAgrid.dx / 2.)
        x = x - 360.0;
      status = map_interpolation(LSAgrid, z, cmask, x, y, &(data2D.LSA[n].z[k]));
      if(status != 0) {
        goto error;
        }
      }
    found[k] = 1;
    zaparr(amp);
    zaparr(pha);
    zaparr(z);
    zaparr(LSAgrid.x);
    zaparr(LSAgrid.y);
    delete[] filename;
    }

  status=0;
  return (status);

error:
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  matrix2x2_t *drag_init_obsolete(parameter_t P1_data, int nndes, float wdHmin , float wdHmax )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#warning drag_init_obsolete is a mutant
{
  int n;
  vector2D_t grd;
  float u, v, w, mean_u, tau[4], N = 2.0e-03, kappa = 2 * M_PI / 1000.;
  float factor, hramp;
  float s, r, f, k_b, k_0, alpha;
//  float hmin = 300, hmax = 500;
//  float wdHmin = 300, wdHmax = 500;
  float wdc1=20,wdc2=0;
  double h;
  matrix2x2_t *matrix;

  alpha = 2.5;
  /*  Slope change at 50 km, k=2pi/50000 */
  k_b = 6.28 / 50000.;
  k_0 = 6.28 / 9000.;

  matrix = new matrix2x2_t[nndes];

  for(n = 0; n < nndes; n++) {
    h = P1_data.h[n];
//    N = P1_state[1].N[n];       /*uninitialized at first time step */
    matrix[n].c[0][0] = 0.0;
    matrix[n].c[0][1] = 0.0;
    matrix[n].c[1][0] = 0.0;
    matrix[n].c[1][1] = 0.0;
    if(h > wdHmin) {
/*----------------------------------------------------------------------
      linear ramp [0:1] between hmin and hmax */
      hramp = MIN(1., (h - wdHmin) / (wdHmax - wdHmin));
/*----------------------------------------------------------------------
      sinusoidale ramp between hmin and hmax */
      hramp = 0.5 * (1. - cos(hramp * M_PI));
      grd.x = P1_data.dhdx[n];
      grd.y = P1_data.dhdy[n];
      kappa = 2 * M_PI / (1500.);
      factor = hramp * P1_data.wdc1[n] * N / kappa;
/*----------------------------------------------------------------------
      divide by h to achieve average force, FL 11/04/2001 */
      matrix[n].c[0][0] += -factor * grd.x * grd.x / h;
      matrix[n].c[0][1] += -factor * grd.x * grd.y / h;
      matrix[n].c[1][0] += -factor * grd.x * grd.y / h;
      matrix[n].c[1][1] += -factor * grd.y * grd.y / h;
      }
    if(wdc2 == 0.0)
      continue;
    if(h > 1500.) {
      r = P1_data.vwr[n];
      f = fabs(2 * P1_data.S[n] * two_Omega);
/*----------------------------------------------------------------------
      kappa depends on tidal excursion, problem to initialise*/
/*
      h=P1_data.h[n];
      u=P1_state[1].u[n];
      v=P1_state[1].v[n];
      w=MAX(sqrt(u*u+v*v+1.e-10),1.0);
*/
      w = 5.0;
      kappa = f / w;
      if(kappa > k_b) {
        s = r * exp((1 - alpha) * log(kappa / k_0));
        }
      else {
        s = r * exp((1 - alpha) * log(k_b / k_0));
        s = s * exp(-0.1 * log(kappa / k_b));
        }
      factor = wdc2 * 2 * N * s * kappa / h;
      matrix[n].c[0][0] += -factor;
      matrix[n].c[1][1] += -factor;
      }
    }
  return (matrix);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int spectral_EnergyComputation_obsolete(const char *meshfile, const char *depthfile, const char *slopefile, const char *belfile, const char *Cdfile, const char *path, const char *output_path, spectrum_t & spectrum)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  char *rootname=NULL,hfile[256]="\0",ufile[256]="\0";
  mesh_t mesh;
  complex<float> *z,*u,*v;
  float *h;
//   spectrum_t spectrum;
  state2d_t state;
  action_t *P1_action;
  parameter_t data;
  parameter_t data2D;
  tidal_wave wave;
  int *found;
  bool debug=false;

/*------------------------------------------------------------------------------
  read mesh*/
  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status!=0) {
      printf("unable to read the original mesh in %s\n",meshfile);
      goto error;
      }
    status=fe_list(&mesh);
    if(status!=0) {
      printf("unable to build the element list from the original mesh\n");
      goto error;
      }
    }
 else {
   printf("no mesh file specified; abort...\n");
   goto error;
   }

  status= fe_edgetable(&mesh,0,0);
  if(status!=0) {
    printf("unable to build the edge list from the original mesh\n");
    goto error;
    }
  
  status= fe_codetable2(&mesh,0,1,0);
  if(status!=0) {
    printf("unable to rebuild the limits table and codes of the original mesh\n");
    goto error;
    }

/*------------------------------------------------------------------------------
  persian gulf experiment specific, commented */    
  status=fe_read_boundarycode(belfile,mesh,0);
  if(status!=0) {
    printf("line %d of %s: cannot read %s\n",__LINE__,__FILE__,belfile);
    goto error;
    }
    
/* *----------------------------------------------------------------------------
  to be replaced*/

  h=new float[mesh.nvtxs];
  z=new complex<float>[mesh.nvtxs];
  u=new complex<float>[mesh.nvtxs];
  v=new complex<float>[mesh.nvtxs];

  status=quoddy_loadr1(depthfile, mesh.nvtxs, h);
  if(status!=0) {
    printf("line %d of %s: cannot read %s\n",__LINE__,__FILE__,depthfile);
    goto error;
    }
    
#if 0
  for(n = 0; n < mesh.nvtxs; n++) {
    h[n] = MAX(10.,h[n]);
    }
#endif

  data2D.lon=new double[mesh.nvtxs];
  data2D.lat=new double[mesh.nvtxs];
  data2D.C  =new double[mesh.nvtxs];
  data2D.S  =new double[mesh.nvtxs];
  data2D.h  =new double[mesh.nvtxs];
  data2D.lmm=new float[mesh.nvtxs];

  data2D.dhdx  =new double[mesh.nvtxs];
  data2D.dhdy  =new double[mesh.nvtxs];
  data2D.Cd    =new float[mesh.nvtxs];
  data2D.u0    =new double[mesh.nvtxs];
  data2D.wdc1  =new float[mesh.nvtxs];

  status = quoddy_loadr2(slopefile, mesh.nvtxs, data2D.dhdx, data2D.dhdy);
  for(n = 0; n < mesh.nvtxs; n++) {
    data2D.dhdx[n] /= 1000.;
    data2D.dhdy[n] /= 1000.;
    }

  data2D.Htide=new hconstant_t[mesh.nvtxs];
  data2D.Utide=new hconstant_t[mesh.nvtxs];
  data2D.Vtide=new hconstant_t[mesh.nvtxs];

  P1_action=new action_t[mesh.nvtxs];

  for(n = 0; n < mesh .nvtxs; n++) {
    data2D.Htide[n].z=new complex<float>[spectrum.n];
    data2D.Htide[n].size=spectrum.n;
    data2D.Utide[n].z=new complex<float>[spectrum.n];
    data2D.Utide[n].size=spectrum.n;
    data2D.Vtide[n].z=new complex<float>[spectrum.n];
    data2D.Vtide[n].size=spectrum.n;
    data2D.LSA[n].z=new complex<float>[spectrum.n];
    }

  for(n = 0; n < mesh .nvtxs; n++) {
    data2D.lon[n]= mesh.vertices[n].lon*d2r;
    data2D.lat[n]= mesh.vertices[n].lat*d2r;
    data2D.C[n]= cos(data2D.lat[n]);
    data2D.S[n]= sin(data2D.lat[n]);
    data2D.h[n]= h[n];
    mesh.vertices[n].h= h[n];//for netcdf
/* *----------------------------------------------------------------------
    temporary patch*/
    data2D.Cd[n]  = 2.5e-03;
    data2D.u0[n]  = 0.05;
    data2D.wdc1[n]= 15.;
    data2D.LSA[n].z[0]=0;
    }
  
  for (int w=0;w<spectrum.n;w++) {
      sprintf(hfile,"%s/%s.ele.%s.s2c",path,spectrum.waves[w].name,rootname);
      sprintf(ufile,"%s/%s.uv.%s.v2c",path, spectrum.waves[w].name,rootname);
      status=quoddy_loadc1(hfile, mesh.nvtxs, z);
      if(status!=0) {
        printf("line %d of %s: cannot read %s : removing %s\n",__LINE__,__FILE__,hfile,spectrum.waves[w].name);
        removewave(&spectrum,w);continue;//rather than
        goto error;
        }
      status=quoddy_loadc2(ufile, mesh.nvtxs, u,v);
      if(status!=0) {
        printf("line %d of %s: cannot read %s : removing %s\n",__LINE__,__FILE__,ufile,spectrum.waves[w].name);
        removewave(&spectrum,w);continue;//rather than
        goto error;
        }
//       Cdfile=strdup("../data/friction-coefficient.s2r");
      status=quoddy_loadr1(Cdfile, mesh.nvtxs, data2D.Cd);
        if(status!=0) {
        printf("warning : cannot read  %s, use default Cd value \n",Cdfile);
        }
      
    for(n = 0; n < mesh .nvtxs; n++) {
      data2D.Htide[n].z[w]=z[n];
      data2D.Utide[n].z[w]=u[n];
      data2D.Vtide[n].z[w]=v[n];
/* *----------------------------------------------------------------------
      temporary patch*/
      P1_action[n].LSA[w]=0;
      }
    }
  for(n = 0; n < mesh .nvtxs; n++) {
    // initialise a (amplitude) and G (phase) in data2D.{H,U,V}tide[n]
    data2D.Htide[n].set_polar();
    data2D.Utide[n].set_polar();
    data2D.Vtide[n].set_polar();
    }

  for(int l = 0; l < mesh.ntriangles; l++) {
    int n1 = mesh.triangles[l].vertex[0];
    int n2 = mesh.triangles[l].vertex[1];
    int n3 = mesh.triangles[l].vertex[2];
    data2D.lmm[n1] += mesh.triangles[l].Area / 3.;
    data2D.lmm[n2] += mesh.triangles[l].Area / 3.;
    data2D.lmm[n3] += mesh.triangles[l].Area / 3.;
    }

  found=new int[spectrum.n];
  status= read_loading02( spectrum,  data2D, mesh.nvtxs, P1_action, found);

  gradient_operator(mesh, data2D);
  DragMatrix = drag_init_obsolete(data2D,mesh.nvtxs,200.,300.);

  fe_integrale_init();

  spectral_transport(mesh, spectrum, data2D, output_path, rootname);

/*------------------------------------------------------------------------------
  persian gulf experiment specific, commented */    
//   sscanf(rootname,"%d",&idx);
// #define showVarVal(f,x) printf(#x "=" #f,x);
// #define showLine printf("line %d of %s: ",__LINE__,__FILE__);
//   showLine showVarVal(%d;,mesh.nvtxs) showVarVal(%d;,mesh.ntriangles) showVarVal(%d;,mesh.nedges) printf("\n");
// 
//   spectral_depth( mesh,spectrum, data2D, output_path, path);
//   printf("line %d of %s: *** SKIPPING spectral_energy01 ***\n",__LINE__,__FILE__);
//   goto end;
  
  spectral_energy01( mesh,spectrum, data2D, P1_action, output_path, path);

  return (status);
error:
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  char *input=0,*keyword,*path=NULL,*output_path=NULL,*s;
  char *meshfile=NULL,*belfile=NULL,*depthfile=NULL,*slopefile=NULL,*Cdfile=NULL,*zonefile=NULL;
  char *output=NULL,*wavename=NULL;
  char *rootname=NULL, *discretisation=0;
  string InputFormat="";
  spectrum_t spectrum(100), dominants(100);
  tidal_wave wave;
  bool compute_energy=true, compute_budget=false;
  bool debug=false;

  fct_echo( argc, argv);

//   spectrum.n=0;
//   spectrum.nmax=100;
//   spectrum.waves=new tidal_wave[spectrum.nmax];

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"--compute-energy=yes")==0) {
      compute_energy=true;
      n++;
      continue;
      }
    if(strcmp(keyword,"--compute-energy=no")==0) {
      compute_energy=false;
      n++;
      continue;
      }
    if(strcmp(keyword,"--compute-budget=yes")==0) {
      compute_budget=true;
      n++;
      continue;
      }
    if(strcmp(keyword,"--compute-budget=no")==0) {
      compute_budget=false;
      n++;
      continue;
      }
    if(strcmp(keyword,"-discretisation")==0) {
      discretisation= strdup(argv[n+1]);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-iF")==0) {
      InputFormat=argv[n+1];
      n++;
      n++;
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 'b' :
          belfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'c' :
          Cdfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'd' :
          depthfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          path= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output_path= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'r' :
          rootname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          slopefile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'w' :
          wavename= strdup(argv[n+1]);
          n++;
          n++;
          wave=wave_from_name(wavename);
          addwave(&spectrum, wave);
          break;

        case 'z' :
          zonefile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        if(input==NULL) {
          input=strdup(argv[n]);
          n++;
          }
        else {
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
          }
        break;
      }
      free(keyword);
    }
  
//   for (k=0;k<spectrum.n;k++) {
//     char *tmp;
//     status=tide_decode_atlasname(path, convention, spectrum.waves[k].name, 0, &tmp);
//     input.push_back(tmp);
//     }

  if(path==NULL) path=strdup(".");
  if(output_path==NULL) output_path=strdup(".");
 
  if(compute_energy) {
    if(InputFormat=="ascii") {
      status=spectral_EnergyComputation_obsolete(meshfile, depthfile, slopefile, belfile, Cdfile, path,  output_path, spectrum);
      }
    else {
      wave=wave_from_name("M2");
      addwave(&dominants, wave);
      status=spectral_EnergyComputation(meshfile, discretisation, path,  output_path, spectrum, dominants);
      }
    }
    
  if(compute_budget) {
    status=spectral_EnergyBudget(input, zonefile, spectrum, debug);
    }
  

error:

end:
  __OUT_BASE_LINE__(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(0);
}
