
/*******************************************************************************

  T-UGO tools, 2006-2016

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada
\author  Yves Soufflet      LEGOS, Toulouse, France

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Concatenate model grid and output files.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>

#include <errno.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "fe.h"
#include "poc-time.h"
#include "list.h"
#include "map.h"
#undef CARTESIAN /* not used here */
#include "map.def"
#include "topo.h"
#include "grd.h"
#include "geo.h"
#undef CARTESIAN /* not used here */
#include "cefmo.h"
#include "functions.h"
#include "tides.h"
#include "netcdf-proto.h"
#include "poc-netcdf-data.hpp"
#include "sturm-liouville.h"

#include "nemo-api.h"

extern double *get_layer_thickness(const string & rootname,const grid_t & grid,int verbose);

extern int VerticalModes_SnapshotDecomposition(const string *input, const string & filename, const string *variables, int maxmodes, bool debug);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void compute_atlas_names(const string & wave,const string & varname,string *filename,string *varnameA,string *varnameG)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  *filename=wave+"-"+varname+"-atlas.nc";
  *varnameA=varname+"_a";
  *varnameG=varname+"_G";
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int compute_critical_slope(grid_t & grid, double omega, float *N3D, double mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
 rough version */
/**----------------------------------------------------------------------------
 
  slope of internal wave characteristics :
 
    alpha = sqrt[(w²-f²)/(N²-w²)] aspect ratio (Gill, p. 259)
    
    actually alpha = cotan(phi); phi angle with vertical; at first order alpha=angle with horizontal
    
    typically alpha = 1.4e-04/2.e-03 = 7.0e-02 (70 m / km, 14km for 1000m)
              alpha = 1.4e-04/1.e-02 = 1.4e-02 (14 m / km, 70km for 1000m)
    
    gamma = grad(h) / alpha
  
  regimes:
 
    sub-critical   : grad(h) < alpha (gamma < 1), propagation is onshore and offshore
    critical       : grad(h) = alpha (gamma = 1), resonant
    super-critical : grad(h) > alpha (gamma > 1), propagation is offshore only
 
-----------------------------------------------------------------------------*/
{
  int status;
  size_t i,j,k,m;
  grid_t topogrid;
  float *topo,*topo_x,*topo_y, topomask;
  float *gamma,*N;
  float *topo_N;
  char *filename;
  size_t size;
  bool debug=false;
  
  filename=strdup("/media/8d0f6fd3-80b8-4696-b618-6f74efb78e3a/models/crozet/bathymetry/merge-cea+ifremer-v4/gridone.grd");
  
  printf("#################################################################\n");
  printf("load bathymetry database\n");
//  status=map_loadfield((const char *) filename, (const char *) 0, &topogrid, &topo, &topomask);
  status=topo_loadfield((const char *) filename, &topogrid, &topo, &topomask, debug);

  printf("#################################################################\n");
  printf("interpolate depth's gradients\n");
  size=topogrid.Hsize();
  topo_x=new float[size];
  topo_y=new float[size];
  status= map_gradient(topogrid, topogrid.nx, topo, topomask, GEOCENTRIC, topo_x, topo_y);

  topo_N=new float[size];
  gamma=new float[size];
  
  N=new float[grid.Hsize()];
  
  for(m=0; m<grid.Hsize(); m++) {
    N[m]=mask;
    for(k=grid.nz-1;k>0;k--) {
      size_t m3D=k*grid.Hsize()+m;
      if(N3D[m3D]!=0) {
        N[m]=N3D[m3D];
        break;
        }
      }
    }
  
/**----------------------------------------------------------------------------
  element-wise computation (such as energy diags) */
  for (j=0; j<topogrid.ny; j++) {
    for (i=0; i<topogrid.nx; i++) {
      m=topogrid.nx*j+i;
      double t,p;
      topogrid.xy(i,j,t,p);
      status=map_interpolation(grid,N,mask,t,p,&topo_N[m]);
      if(topo_N[m]==mask) {
        gamma[m]=topomask;
        continue;
        }
/**----------------------------------------------------------------------------
      compute gamma */
      float s=sqrt(topo_x[m]*topo_x[m]+topo_y[m]*topo_y[m]);
      float alpha=2*omega/topo_N[m];
      gamma[m]=s/alpha;
      if(gamma[m]==0) {
        gamma[m]=topomask;
        continue;
        }
      }
    }
  status=map_completegridaxis(&topogrid,1);
  status=save_SGXY("critical.nc", topogrid, gamma, topomask, "gamma", "none", "gamma", 1);
  status=save_SGXY("critical.nc", topogrid, topo_N, topomask, "N", "none", "N", 0);
  
  delete[] topo;
  delete[] topo_x;
  delete[] topo_y;
  delete[] topo_N;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int VerticalModes_normalize(const string & rootname,int maxmodes, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,n,status;
  int verbose=0;
  grid_t grid, wgrid;
  double *modes=0, *z, mask;
  int nkept, nmodes;
  string filename, varname;
  poc_data_t<double> scalarData;
  
  nkept=5;
  nmodes=maxmodes;
  
/*------------------------------------------------------------------------------
  */
  filename=rootname+".vertical-modes.nc";
  
  varname="Wmodes";
  status=poc_get_grid(filename, varname, &wgrid, verbose, -1);
  z=new double[wgrid.nz];
  status=scalarData.init(filename,varname);
  for(i=0;i<nmodes;i++) {
    status=scalarData.read_data(filename,i);
    nmodes=scalarData.nframes;
    if(modes==0) {
      modes=new double[wgrid.nz];
      }
    mask=scalarData.mask;
    for(n=0;n<wgrid.Hsize();n++) {
      for(k=0;k<wgrid.nz;k++) {
        size_t m=k*wgrid.Hsize()+n;
        modes[k]=scalarData.data[m];
        z[k]=wgrid.z[m];
        }
      status=normalize_Wmodes(wgrid.nz, z, &modes, mask, 1, 0, 0, verbose);
      for(k=0;k<wgrid.nz;k++) {
        size_t m=k*wgrid.Hsize()+n;
        scalarData.data[m]=modes[k];
        }
      }
    status=scalarData.write_data(filename,i);
    }
  
  delete[] modes;
  modes=0;
  
  scalarData.destroy_data();
//   grid.free();
  
  varname="Pmodes";
  status=poc_get_grid(filename, varname, &grid, verbose, -1);
//   z=new double[grid.nz];
  status=scalarData.init(filename,varname);
  for(i=0;i<nmodes;i++) {
    status=scalarData.read_data(filename,i);
    nmodes=scalarData.nframes;
    if(modes==0) {
      modes=new double[grid.nz];
      }
    mask=scalarData.mask;
    for(n=0;n<grid.Hsize();n++) {
      for(k=0;k<grid.nz;k++) {
        size_t m=k*grid.Hsize()+n;
        modes[k]=scalarData.data[m];
        }
      for(k=0;k<wgrid.nz;k++) {
        size_t m=k*wgrid.Hsize()+n;
        z[k]=wgrid.z[m];
        }
      status=normalize_Pmodes(wgrid.nz, z, &modes, mask, 1, 0, verbose);
      for(k=0;k<grid.nz;k++) {
        size_t m=k*grid.Hsize()+n;
        scalarData.data[m]=modes[k];
        }
      }
    status=scalarData.write_data(filename,i);
    }
  
  delete[] modes;
  modes=0;
  
  scalarData.destroy_data();
  grid.free();
  
  varname="Umodes";
  
  status=poc_get_grid(filename, varname, &grid, verbose, -1);
//   z=new double[grid.nz];
  status=scalarData.init(filename,varname);
  for(i=0;i<nmodes;i++) {
    status=scalarData.read_data(filename,i);
    nmodes=scalarData.nframes;
    if(modes==0) {
      modes=new double[grid.nz];
      }
    mask=scalarData.mask;
    for(n=0;n<grid.Hsize();n++) {
      for(k=0;k<grid.nz;k++) {
        size_t m=k*grid.Hsize()+n;
        modes[k]=scalarData.data[m];
        z[k]=grid.z[m];
        }
      for(k=0;k<wgrid.nz;k++) {
        size_t m=k*wgrid.Hsize()+n;
        z[k]=wgrid.z[m];
        }
      status=normalize_Umodes(wgrid.nz, z, &modes, mask, 1, 0, verbose);
      for(k=0;k<grid.nz;k++) {
        size_t m=k*grid.Hsize()+n;
        scalarData.data[m]=modes[k];
        }
      }
    status=scalarData.write_data(filename,i);
    }
  
  delete[] modes;
  modes=0;
  
  scalarData.destroy_data();
  grid.free();
  
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int compute_vertical_mode(const char *bathymetry, int nRequestedProcs, grid_t & t_grid, grid_t & w_grid, float *T, float *S, float *R, float mask,const char *rootname, int maxmodes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int j,k,l,m;
  char filename[1024];
  float **celerity, **Wmodes, *Nbar=0, *N=0;
  double *z;
  int *nmodes;
  grid_t topogrid;
  float *topo=0,topomask;
  double f,omega = 7.292e-5;
  bool rigid_lid=false;
  
  if(isfinite(ctrlLon) and isfinite(ctrlLat)){
    set_grid_list(&t_grid);
    index_interpolation(t_grid,ctrlLon,ctrlLat,&ctrlM,t_grid.x,(double)NAN,(double*)0,0);
    t_grid.ij(ctrlM,&ctrlI,&ctrlJ);
    }
  if(ctrlM>=0){
    ctrlF=fopen("ctrlF.py","w");
    ASSERT_ARRAY(ctrlF,"m",&ctrlM,1,1);
    }
  
/**-----------------------------------------------------------------------------
  get additional bathymetry if needed (basically for observational climatology database)*/
//   if(bathymetry==0) bathymetry=strdup("/home/softs/genesis/data/topography/gebco/gridone.grd");
  if(bathymetry!=NULL) {
    status=grd_loadgrid(bathymetry,&topogrid);
    if(status !=0) {
      STDOUT_BASE_LINE("cannot load bathymetry file=%s\n",bathymetry);
      exit(-1);
      }
    exitIfNull(
      topo=new float[topogrid.nx*topogrid.ny]
      );
    topogrid.modeH=0;
    status=  grd_loadr1(bathymetry,topogrid,topo,&topomask);
    }
  else {
    topo=0;
    }
  
  sprintf(filename,"%s.verification.nc",rootname);
  status=poc_def_att(filename,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
  
  if(R!=0){
    status=save_SGXYZ(filename, t_grid, R, mask, "RHO", "kg/m3", "density", 1, 1, "", "");
    if(t_grid.modeV==1)
      ASSERT_ARRAY(ctrlF,"z",&t_grid.z[ctrlM],t_grid.nz,1);
    else if(t_grid.modeV==3)
      ASSERT_ARRAY(ctrlF,"z",&t_grid.z[ctrlM],t_grid.nz,t_grid.Hsize());
    ASSERT_ARRAY(ctrlF,"R",&R[ctrlM],t_grid.nz,t_grid.Hsize());
    }
  
  int nprocs=initialize_OPENMP(nRequestedProcs);
  
/**----------------------------------------------------------------------------
  compute Brünt-Vaïsala frequency */
  if(R==0) TRAP_ERR_EXIT(-1,"MUST have density to continue...\n");
  Nbar=new float[t_grid.Hsize()];
  N=new float [(t_grid.nz+1)*t_grid.Hsize()];
  z=new double[(t_grid.nz+1)*t_grid.Hsize()];
  
  fflush(stdout);
  
#pragma omp parallel for private(m,status) if(nprocs>1)
  for(j=0;j<t_grid.ny;j++) {
    for(int i=0;i<t_grid.nx;i++) {
//       int m=j*t_grid.nx+i;
      status=BruntVassala(t_grid, w_grid, T, S, R, mask, topogrid, topo, topomask, i, j, N, z);
      }
    }
  
  fflush(stdout);
  
  grid_t wgrid=t_grid;
  wgrid.nz++;
  wgrid.z=z;
  status=save_SGXYZ(filename, wgrid, N, mask, "N", "s-1", "N", 0, 1, "", "w");
  
//   status=compute_critical_slope(grid, omega, N, (float) mask);
  
  if(maxmodes==-1) maxmodes=wgrid.nz;
  
/**----------------------------------------------------------------------------
  compute vertical modes */
  celerity=new float *[maxmodes];
  for(k=0;k<maxmodes;k++) {
    celerity[k]=new float[t_grid.ny*t_grid.nx];
    }
  
  Wmodes=new float *[maxmodes];
  for(k=0;k<maxmodes;k++) {
    Wmodes[k]=new float[wgrid.nz*t_grid.ny*t_grid.nx];
    }
  for(k=0;k<maxmodes;k++) {
    for(j=0;j<t_grid.ny;j++) {
      for(int i=0;i<t_grid.nx;i++) {
        m=j*t_grid.nx+i;
        celerity[k][m]=mask;
        }
      }
    for(m=0;m<wgrid.size();m++) {
      Wmodes[k][m]=mask;
      }
    }
  
  nmodes=aset(wgrid.Hsize(),(int) mask);
  
  bool wascending;
  double wfactor;
  
  status=check_vertical_direction(wgrid,&wascending,&wfactor);
  
#pragma omp parallel for if(nprocs>1)
  for(j=0;j<t_grid.ny;j++) {
    int status;
    for(int i=0;i<t_grid.nx;i++) {
//       size_t m=j*t_grid.nx+i;
      status=SG_Wmodes_v2(wgrid, wascending, wfactor, T, S, R, mask, topogrid, topo, topomask, i, j, Nbar, celerity, Wmodes, nmodes, maxmodes, rigid_lid);
      if(status!=0) {
        if(status!=-1) printf("vertical modes failed: i=%d j=%d status=%d\n",i,j,status);
        }
      }
    }

  sprintf(filename,"%s.vertical-modes.nc",rootname);
  status=poc_def_att(filename,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
  int range[2]={0,maxmodes-1};

/*------------------------------------------------------------------------------
  3D W-mode profiles */
  status=save_SGXYZT(filename, wgrid, Wmodes, (float) mask, "Wmodes",   "dimensionless",   "Wmodes", 1, 1, "w", "w",range);
  status=save_SGXY((const char*) filename, (const grid_t) wgrid, nmodes, (int)   mask, "nWmodes",   "dimensionless",   "nWmodes", 0, 0);
  
/*------------------------------------------------------------------------------
  3D P-mode and U-mode profiles */
  float *Pmodes=new float[t_grid.size()];
  float *Umodes=new float[t_grid.size()];
  for(k=0;k<maxmodes;k++) {
    for(j=0;j<t_grid.ny;j++) {
      for(int i=0;i<t_grid.nx;i++) {
        m=j*t_grid.nx+i;
        for(l=0;l<t_grid.nz;l++) {
          int n =l*t_grid.nx*t_grid.ny+m;
          int n1=l*t_grid.nx*t_grid.ny+m;
          int n2=(l+1)*t_grid.nx*t_grid.ny+m;
          double dH=z[n2]-z[n1];
          if( (Wmodes[k][n2] != mask) and (Wmodes[k][n1] != mask) ) {
            Pmodes[n]=celerity[k][m]*celerity[k][m]*R[n]*(Wmodes[k][n2]-Wmodes[k][n1])/dH;
            Umodes[n]=Pmodes[n]/9.81/R[n];
            }
          else {
            Pmodes[n]=mask;
            Umodes[n]=mask;
            }
          }
        }
      }
    int frame=k;
    int create_file=0;
    int create_grid=(k==0);
    status=save_SGXYZT(filename, t_grid, Pmodes, (float) mask, "Pmodes", "dimensionless", "Pmodes", create_file, create_grid, "", "",frame);
    create_grid=0;
    status=save_SGXYZT(filename, t_grid, Umodes, (float) mask, "Umodes", "dimensionless", "Umodes", create_file, create_grid, "", "",frame);
    }
  delete[] Pmodes;
  delete[] Umodes;
  
  status=save_SGXY((const char*) filename, (const grid_t) wgrid, nmodes, (int)   mask, "nWmodes",   "dimensionless",   "nWmodes", 0, 0);

/*------------------------------------------------------------------------------
  2D mode celerity */
  status=save_SGXYT(filename, t_grid, celerity, (float) mask, "celerity", "m/s", "celerity", 0, 1, range);

/*------------------------------------------------------------------------------
  2D M2 wavelength */
  for(k=0;k<maxmodes;k++) {
    for(j=0;j<t_grid.ny;j++) {
      for(int i=0;i<t_grid.nx;i++) {
        m=j*t_grid.nx+i;
        if(celerity[k][m]!=mask) celerity[k][m]*=12.5*3.6;
        }
      }
    }
  
  status=save_SGXYT(filename, t_grid, celerity, (float) mask, "M2_lamda", "km", "M2_lamda", 0, 1, range);
  
/*------------------------------------------------------------------------------
  2D Rossby radius */
  for(k=0;k<maxmodes;k++) {
    for(j=0;j<t_grid.ny;j++) {
      for(int i=0;i<t_grid.nx;i++) {
        m=j*t_grid.nx+i;
        if(celerity[k][m]!=mask) {
          celerity[k][m]/=12.5*3.6;
          double x,y;
          t_grid.xy(i,j,x,y);
          if(y>5 || y<(-5)) {
            f=fabs(2.*omega*sin(y*M_PI/180.));
            celerity[k][m]/=f*1000;
            }
          else {
            f=4.*omega*cos(y*M_PI/180.)/6371;
            celerity[k][m]=sqrt(celerity[k][m]/(f*1000));
            }
          }
        }
      }
    }
  
  status=save_SGXYT(filename, t_grid, celerity, (float) mask, "rossby_radius", "km", "rossby_radius", 0, 1, range);

  status=save_SGXY(filename, t_grid, Nbar, (float) mask, "Nbar_1", "s^-1", "Nbar_1", 0);
  
  delete[] Nbar;
  delete[] N;
  
  for(k=0;k<maxmodes;k++) {
    delete[] celerity[k];
    }
  delete[] celerity;
  
  for(k=0;k<maxmodes;k++) {
    delete[] Wmodes[k];
    }
  delete[] Wmodes;
  
  delete[] topo;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int levitus(const string *input,const char *rootname,const char *bathymetry, int nRequestedProcs, int maxmodes, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int i,j,k,m;

  char filename[1024]="/home/data/climatology/levitus/data-1994.nc";
  grid_t grid, w_grid;
  float *R=0, *T=0, *S=0;
  float mask;

  if(input[RHO_ID]!="") strcpy(filename,input[RHO_ID].c_str());

  status= map_loadfield3D(filename, filename, "potdens", grid, R, mask, verbose);
  if(status !=0) {
    STDOUT_BASE_LINE("cannot load climatology file=%s\n",filename);
    exit(-1);
    }
  for(k=0;k<grid.nz;k++) {
    for(j=0;j<grid.ny;j++) {
      for(i=0;i<grid.nx;i++) {
        m=k*grid.ny*grid.nx+j*grid.nx+i;
        if(R[m]!=mask) {
          R[m]=R[m]*1000.+1000.;
          }
        }
      }
    }
  
//   int maxmodes=5;
  status= compute_vertical_mode(bathymetry, nRequestedProcs, grid, w_grid, T, S, R, mask, rootname, maxmodes);
  
  delete[] R;
  
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int WOA_annual(const string *input,const string *variables,const char *rootname,const char *bathymetry, int nRequestedProcs, int maxmodes, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int j,k,m;
  char filename[1024];
  grid_t grid, w_grid;
  float *R,*T,*S, mask;
  grid_t topogrid;
 
  strcpy(filename,input[TEM_ID].c_str());
  
  status= map_loadfield3D(filename, filename, variables[0].c_str(), grid, T, mask, verbose);
  if(status !=0) {
    STDOUT_BASE_LINE("cannot load %s climatology file=%s\n", "temperature", filename);
    exit(-1);
    }
  
  strcpy(filename,input[SAL_ID].c_str());
  
  status= map_loadfield3D(filename, filename, variables[1].c_str(), grid, S, mask, verbose);
  if(status !=0) {
    STDOUT_BASE_LINE("cannot load %s climatology file=%s\n", "salinity", filename);
    exit(-1);
    }
  
  R=new float[grid.nz*grid.ny*grid.nx];
  for(k=0;k<grid.nz;k++) {
    for(j=0;j<grid.ny;j++) {
      for(int i=0;i<grid.nx;i++) {
        m=k*grid.ny*grid.nx+j*grid.nx+i;
        if(T[m]!=mask) {
          R[m]=water_density(T[m],S[m],grid.z[m]);
// 	   if(R[m]==0.) {
// 	     printf("%d %d %d\n",i,j,k);
//              }
          }
        else {
          R[m]=mask;
          }
        }
      }
    }
  
//   int maxmodes=5;
  status= compute_vertical_mode(bathymetry, nRequestedProcs, grid, w_grid, T, S, R, mask, rootname, maxmodes);
  
  delete[] R;
  
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int WOA_monthly(const string *input,const string *variables,const char *rootname, int frame,const char *bathymetry, int nRequestedProcs, int maxmodes, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int j,k,m;
  char filename[1024];
  grid_t grid, w_grid;
  float *R,*T,*S, mask;
  grid_t topogrid;
 
  strcpy(filename,input[TEM_ID].c_str());
  
  status= map_loadfield3D(filename, filename, variables[0].c_str(), frame, grid, T, mask, verbose);
  if(status !=0) {
    STDOUT_BASE_LINE("cannot load %s climatology file=%s\n", "temperature", filename);
    exit(-1);
    }
  
  strcpy(filename,input[SAL_ID].c_str());
  
  status= map_loadfield3D(filename, filename, variables[1].c_str(), frame, grid, S, mask, verbose);
  if(status !=0) {
    STDOUT_BASE_LINE("cannot load %s climatology file=%s\n", "salinity", filename);
    exit(-1);
    }
  
  R=new float[grid.nz*grid.ny*grid.nx];
  for(k=0;k<grid.nz;k++) {
    for(j=0;j<grid.ny;j++) {
      for(int i=0;i<grid.nx;i++) {
        m=k*grid.ny*grid.nx+j*grid.nx+i;
        if(T[m]!=mask) {
          R[m]=water_density(T[m],S[m],grid.z[m]);
// 	   if(R[m]==0.) {
// 	     printf("%d %d %d\n",i,j,k);
//              }
          }
        else {
          R[m]=mask;
          }
        }
      }
    }
  
//   int maxmodes=5;
  status= compute_vertical_mode(bathymetry, nRequestedProcs, grid, w_grid, T, S, R, mask, rootname, maxmodes);
  
  delete[] R;
  
  
  return(0);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int get_w_grid(const string & filename,const grid_t & grid, grid_t & w_grid, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=poc_get_grid(filename,"depth_w", &w_grid, verbose, 0);
  NC_CHKERR_BASE_LINE(status,"poc_get_grid(\""+filename+"\",\"depth_w\",,,0) error");
  if(status!=0){
    bool increasing;
    double factor;
    
    status=check_vertical_direction(grid,&increasing,&factor);
    if(status!=0) TRAP_ERR_EXIT(status,"check_vertical_direction() error\n");
    
    STDERR_BASE_LINE("MAKING w GRID from increasing=%d, factor=%g\n",increasing,factor);
    int i,k,m,mw;
    const int
      n=grid.Hsize();
    
    w_grid=grid;
    w_grid.nz++;
    w_grid.z=aset(w_grid.size(),w_grid.zmask);
    
    poc_data_t<double> thickness;
    status=thickness.init(filename,"h",verbose);
    
    if(increasing){
      m=0;
      mw=m;
      }
    else{
      m=(grid.nz-1)*n;
      mw=m+n;
      }
    
    for(i=0;i<n;i++,m++,mw++){
      const double
        zm=grid.z[m],
        thicknessm=thickness.data[m];
      
      if(zm==grid.zmask or thicknessm==thickness.mask){
        w_grid.z[mw]=w_grid.zmask;
        continue;
        }
      
      w_grid.z[mw]=factor*zm-thicknessm*.5;
      }
    
    for(k=0;k<grid.nz;k++){
      m=k*n;
      mw=m;
      if(increasing) mw+=n;
      
      for(i=0;i<n;i++,m++,mw++){
        const double
          zm=grid.z[m],
          thicknessm=thickness.data[m];
        
        if(zm==grid.zmask or thicknessm==thickness.mask){
          w_grid.z[mw]=w_grid.zmask;
          continue;
          }
        
        w_grid.z[mw]=factor*zm+thicknessm*.5;
        }
      }
    
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int SYMPHONIE_spectral(const string *input,const char *rootname,const string varnames[5],const char *bathymetry, int nRequestedProcs, const char *units, int maxmodes, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,m,n;
  string filename;
  grid_t grid, w_grid;
  float *R=0,*T=0,*S=0, mask;
  
/**----------------------------------------------------------------------------
  get density */
  filename=input[RHO_ID];
  poc_data_t<float> density;
  status=density.init(filename,varnames[3],verbose);
  NC_TRAP_ERROR(wexit,status,1,"poc_data_t::init(\""+filename+"\",\""+varnames[3]+"\",) error");
  status=poc_get_grid(filename,density.info, &grid, verbose, 0);
  NC_TRAP_ERROR(wexit,status,1,"poc_get_grid(\""+filename+"\",(\"rhp\"),,,0) error");
  mask=density.mask;
  swapval(R,density.data);
  
  n=grid.size();
  for(m=0;m<n;m++){
    float *Rm=&R[m];
    if(*Rm==mask) continue;
    *Rm+=1000.f;
    }
  
  get_w_grid(filename,grid,w_grid,verbose);
  
  status=compute_vertical_mode(bathymetry, nRequestedProcs, grid, w_grid, T, S, R, mask, rootname, maxmodes);
  
  delete[] R;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ECCO_all(const string *input,const char *rootname,const char *bathymetry, int nRequestedProcs, int maxmodes, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int j,k,m;
  char *filename;
  grid_t grid, w_grid;
  float *R,*T,*S, mask;
  grid_t topogrid;
  
  filename=new char[1024];
 
  strcpy(filename,input[TEM_ID].c_str());
  
  status= map_loadfield3D(filename, filename, "THETA", grid, T, mask, verbose);
  if(status !=0) {
    STDOUT_BASE_LINE("cannot load %s climatology file=%s\n", "temperature", filename);
    exit(-1);
    }
  
  strcpy(filename,input[SAL_ID].c_str());
  
  status= map_loadfield3D(filename, filename, "SALT", grid, S, mask, verbose);
  if(status !=0) {
    STDOUT_BASE_LINE("cannot load %s climatology file=%s\n", "salinity", filename);
    exit(-1);
    }
  
  for(k=0;k<grid.nz;k++) {
    for(j=0;j<grid.ny;j++) {
      for(int i=0;i<grid.nx;i++) {
        m=k*grid.ny*grid.nx+j*grid.nx+i;
        if(T[m]<-1000) {
          T[m]=mask;
          }
        if(S[m]<-1000) {
          S[m]=mask;
          }
        }
      }
    }
  
  R=new float[grid.nz*grid.ny*grid.nx];
  for(k=0;k<grid.nz;k++) {
    for(j=0;j<grid.ny;j++) {
      for(int i=0;i<grid.nx;i++) {
        m=k*grid.ny*grid.nx+j*grid.nx+i;
        if(T[m]!=mask) {
          R[m]=water_density(T[m],S[m],grid.z[m]);
// 	   if(R[m]==0.) {
// 	     printf("%d %d %d\n",i,j,k);
//              }
          }
        else {
          R[m]=mask;
          }
        }
      }
    }
  
//   int maxmodes=5;
  status= compute_vertical_mode(bathymetry, nRequestedProcs, grid, w_grid, T, S, R, mask, rootname, maxmodes);
  
  delete[] R;
  
  delete[] filename;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int SYMPHONIE(const string *input,const char *rootname,const char *bathymetry, int nRequestedProcs, const char *units, int maxmodes, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int j,k,m;
  char datafile[1024], gridfile[1024];
  grid_t t_grid, w_grid;
  float *R=0, *T=0, *S=0, *depth_w=0;
  float mask;
  grid_t topogrid;
  int count=0;

  strcpy(datafile,input[TEM_ID].c_str());
  
  if(input[G_ID]=="") strcpy(gridfile,input[TEM_ID].c_str());
  else strcpy(gridfile,input[G_ID].c_str());
  
  status= map_loadfield3D(datafile, gridfile, "tem", t_grid, T, mask, verbose);
  if(status !=0) {
    STDOUT_BASE_LINE("cannot load %s climatology file=%s\n", "temperature", datafile);
    exit(-1);
    }
  
  strcpy(datafile,input[SAL_ID].c_str());
  if(input[G_ID]=="") strcpy(gridfile,input[SAL_ID].c_str());
  else strcpy(gridfile,input[G_ID].c_str());
  
  status= map_loadfield3D(datafile, gridfile, "sal", t_grid, S, mask, verbose);
  if(status !=0) {
    STDOUT_BASE_LINE("cannot load %s climatology file=%s\n", "salinity", datafile);
    exit(-1);
    }

  status= map_loadfield3D(gridfile, gridfile, "depth_w", w_grid, depth_w, mask, verbose);
  if(status !=0) {
    STDOUT_BASE_LINE("cannot load %s climatology file=%s\n", "depth_w", datafile);
    exit(-1);
    }

  if(strcmp(units,"K")==0) {
    for(k=0;k<t_grid.nz;k++) {
      for(j=0;j<t_grid.ny;j++) {
        for(int i=0;i<t_grid.nx;i++) {
          m=k*t_grid.ny*t_grid.nx+j*t_grid.nx+i;
          if(T[m]!=mask) {
            T[m]-=273.15;
            }
          }
        }
      }
    }
  
  R=new float[t_grid.size()];
  for(k=0;k<t_grid.nz;k++) {
    for(j=0;j<t_grid.ny;j++) {
      for(int i=0;i<t_grid.nx;i++) {
        m=k*t_grid.ny*t_grid.nx+j*t_grid.nx+i;
        if(T[m]!=mask) {
          if(T[m]<0.0) {
            if(verbose) printf("T anomaly at i=%4d j=%4d k=%4d : T=%f S=%f \n",i,j,k,T[m],S[m]);
            count++;
            T[m]=4.;
            }
          R[m]=water_density(T[m],S[m],t_grid.z[m]);
          if(R[m]>1050.) {
            if(verbose) printf("R anomaly at i=%d j=%d k=%d : T=%f S=%f R=%f \n",i,j,k,T[m],S[m],R[m]);
            }
          }
        else {
          R[m]=mask;
          }
        }
      }
    }
  
  poc_data_t<float> difv;
  
  
  size_t nvalids=occurence(mask, T, t_grid.size());
  if(count>0) printf("warning: %d T anomaly ((negative temperature) found (over %d)\n",count,nvalids);

  range_t<float> range=range_t<float> (R,mask,t_grid.size());
  
  status= compute_vertical_mode(bathymetry, nRequestedProcs, t_grid, w_grid, T, S, R, mask, rootname, maxmodes);

  delete[] R;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int compute_vertical_mode_template(mesh_t mesh, T state, int paire)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//   int status;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int compute_vertical_mode(mesh_t mesh, ugo_state_t state, int paire)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=compute_vertical_mode_template(mesh, state, paire);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int compute_vertical_mode(mesh_t mesh, tide3D_t state, int paire)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=compute_vertical_mode_template(mesh, state, paire);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : argv[0]
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    "NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Computes vertical modes, spectral decomposition and batrotropic/batroclinic energy fluxes from harmonic.\n"
    "\n"
    "OPTIONS\n"
    "  -h,--help  show this help and exit\n"
    "  --debug  set debug mode\n"
    "  --compute-modes=[yes|no] whether of not to (re)compute vertical modes (for NEMO and SYMPHONIE models only). Default: yes. Tip: they are stored in a file named <rootname>.vertical-modes.nc\n"
    "  --decomposition=[yes|no] whether of not to (re)compute spectral decomposition (for NEMO and SYMPHONIE models only). Default: yes. Tip: it is stored in a file named modal-decomposition.nc\n"
    "  --recomposition=[yes|no] whether of not to (re)compute spectral recomposition (for NEMO models only). Default: yes. Tip: it is stored in a file named modal-recomposition.nc\n"
    "  -nmodes  followed by the number of vertical modes to compute. Default: 10\n"
    "  -nrecomposed  followed by . Default: 10\n"
    "  -nprocs  followed by the number of (logical) cores to use for parallelisation, in particular for vertical mode computation. Default: use all\n"
    "  -unstructured  followed by the discretisation of the unstructured grid\n"
    "  -tfile FILE : get temperature from FILE\n"
    "  -sfile FILE : get salinity from FILE\n"
    "  -rfile FILE : get POTENTIAL density from FILE\n"
    "  -ufile FILE : get zonal velocity from FILE\n"
    "  -gfile FILE : get grid from FILE\n"
    "  -tvar VAR : use VAR as temperature variable name (for WOA and NEMO models only)\n"
    "  -svar VAR : use VAR as salinity variable name (for WOA and NEMO models only)\n"
    "  -rvar VAR : use VAR as POTENTIAL density variable name (for NEMO models only)\n"
    "  --verbose  set verbose mode\n"
    "  --equation_of_state  UNUSED\n"
    "  -b  followed by the bathymetry\n"
    "  -m  followed by the model name, one of: NEMO SYMPHONIE SYMPHONIE-spectral ORCA ORCA12 ORCA12-monthly GLORYS WOA2005 WOA2009 ECCO LEVITUS\n"
    "  -r  followed by the root name. Default: same as the model name (see above)\n"
    "  -v  followed by one list of space-separated variable names (for NEMO and SYMPHONIE models only). Default: \"vovecrtz vozocrtx vomecrty Pbc\"\n"
    "  -w  followed by the tidal wave name (for NEMO and SYMPHONIE-spectral models only). Default: M2\n"
    "\n"
    "BUGS\n"
    "  Many file or variable names are model-dependant! This means that reading this documentation is likely not enough to know how to use this program...\n"
    ); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// function that will be ran when the executable is started
/** See the print_help <a href=#func-members>function</a> for its use.
*/
/*----------------------------------------------------------------------------*/
{
  int i,n,status;
  char *rootname=NULL;
  const char *keyword=NULL;
  const char *equation_of_state=NULL;
  date_t date;
  string input[10], variables[10], grid_variables[10], model;
  
  const char *bathymetry=NULL;
  int   structured=1;

  char *discretisation;
  int nRequestedProcs=-1, maxmodes=10, nrecomposed=10;
  int verbose=0;
  
  bool compute_modes=true, decomposition=true, recomposition=true;
  string wave="M2";
  string varnames[5];
  bool debug=false;
  
  keyword=getenv("LON");
  if(keyword!=0)
    status=sscanf(keyword,"%lg",&ctrlLon);
  keyword=getenv("LAT");
  if(keyword!=0)
    status=sscanf(keyword,"%lg",&ctrlLat);
  
  fct_echo( argc, argv);
  
  varnames[0]="vovecrtz";
  varnames[1]="vozocrtx";
  varnames[2]="vomecrty";
  varnames[3]="Pbc";

  n=1;
  while (n < argc) {
    keyword=argv[n];
    switch (keyword[0]) {
      case '-':
        if( strcmp(keyword,"--help")==0 or strcmp(keyword,"-h")==0 ){
          print_help(argv[0]);
          exit(0);
          }
        if(strcmp(keyword,"--debug")==0) {
          debug=true;
          n++;
          break;
          }
        if(strcmp(keyword,"--compute-modes=yes")==0) {
          compute_modes=true;
          n++;
          break;
          }
        if(strcmp(keyword,"--compute-modes=no")==0) {
          compute_modes=false;
          n++;
          break;
          }
        if(strcmp(keyword,"--decomposition=yes")==0) {
          decomposition=true;
          n++;
          break;
          }
        if(strcmp(keyword,"--decomposition=no")==0) {
          decomposition=false;
          n++;
          break;
          }
        if(strcmp(keyword,"--recomposition=yes")==0) {
          recomposition=true;
          n++;
          break;
          }
        if(strcmp(keyword,"--recomposition=no")==0) {
          recomposition=false;
          n++;
          break;
          }
        if(strcmp(keyword,"-nmodes")==0) {
          sscanf(argv[n+1],"%d",&maxmodes);
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"-nrecomposed")==0) {
          sscanf(argv[n+1],"%d",&nrecomposed);
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"-nprocs")==0) {
          sscanf(argv[n+1],"%d",&nRequestedProcs);
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"-unstructured")==0) {
          discretisation= argv[n+1];
          structured=0;
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"-tfile")==0) {
          input[TEM_ID]= argv[n+1];
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"-sfile")==0) {
          input[SAL_ID]= argv[n+1];
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"-rfile")==0) {
          input[RHO_ID]= argv[n+1];
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"-ufile")==0) {
          input[U_ID]= argv[n+1];
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"-gfile")==0) {
          input[G_ID]= argv[n+1];
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"-mfile")==0) {
          input[M_ID]= argv[n+1];
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"-tvar")==0) {
          variables[TEM_ID]= argv[n+1];
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"-svar")==0) {
          variables[SAL_ID]= argv[n+1];
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"-rvar")==0) {
          variables[RHO_ID]= argv[n+1];
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"-u_var")==0) {
          variables[U_ID]= argv[n+1];
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"-v_var")==0) {
          variables[V_ID]= argv[n+1];
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"-w_var")==0) {
          variables[W_ID]= argv[n+1];
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"-ssh_var")==0) {
          variables[SSH_ID]= argv[n+1];
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"--verbose")==0) {
          verbose=1;
          n++;
          break;
          }
        if(strcmp(keyword,"--equation_of_state")==0) {
          equation_of_state= argv[n+1];
          n++;
          n++;
          break;
          }
        switch (keyword[1]) {

/*------------------------------------------------------------------------------
        bathymetry database*/
        case 'b' :
          bathymetry= argv[n+1];
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        naming convention for tidal atlases*/
        case 'm' :
          model= argv[n+1];
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        output file name*/
        case 'r' :
          rootname= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        harmonic varnames*/
        case 'v' :
          {
          string tmp=argv[n+1];
          vector<string> tokens=string_split(tmp," ");
          for(int k=0;k<min((size_t) 5,tokens.size());k++) {
            varnames[k]=tokens[k];
            }
          }
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        tidal wave*/
        case 'w' :
          wave=argv[n+1];
          n++;
          n++;
          break;

        default:
          printf("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
          break;
        }
        break;

      default:
        printf("unknown option %s\n",keyword);
        print_help(argv[0]);
        exit(-1);
        break;
      }
    }
  
  
/*-----------------------------------------------------------------------------
  argument checks */
  
  status=0;
  
  if(model==""){
    printf("*** Please specify the model name with the -m option ***\n");
    status=-1;
    }

#if 0
  if(==0){
    printf("*** Please specify the  with the - option ***\n");
    status=-1;
    }
#endif
  
  if(status!=0){
    print_help(argv[0]);
    exit(-1);
    }
  
/*----------------------------------------------------------------------------*/
  
  if(rootname==0) rootname=strdup(model.c_str());
  
  water_density_check(0);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  vertical modes, 3 steps processing:
  
  1- compute vertical modes from climatologic density (u/v,w,P modes)
  2- decompose/recompose harmonic solutions
  3- compute energy budget
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(model=="NEMO-snapshot") {
    string ModesFile=(string) rootname+".vertical-modes.nc";
/*------------------------------------------------------------------------------
    mask filename to be pushed in options */
    string maskFile="mesh_mask.nc";
    if(input[M_ID]!="") maskFile=input[M_ID];
    string e3t_var="e3t_0",e3u_var="e3u_0",e3v_var="e3v_0",e3w_var="e3w_0";
    string filenames[4];
    status=0;
    if(compute_modes) {
      status=NEMO(input, maskFile, grid_variables, variables, rootname, bathymetry, nRequestedProcs, "C", maxmodes, verbose);
//       status=NEMO_FV(input, variable, rootname, bathymetry, nRequestedProcs, "C", maxmodes,verbose);
      status=VerticalModes_normalize(rootname, maxmodes, debug);
      }
    if(input[U_ID]!="") status=VerticalModes_SnapshotDecomposition(input, "NEMO.vertical-modes.nc", variables, maxmodes, debug);
    }
  else if(model=="NEMO") {
    string ModesFile=(string) rootname+".vertical-modes.nc";
/*------------------------------------------------------------------------------
    mask filename to be pushed in options */
    const string maskFile="mesh_mask.nc";
    string e3t_var="e3t_0",e3u_var="e3u_0",e3v_var="e3v_0",e3w_var="e3w_0";
    string filenames[4];
    status=0;
    if(compute_modes) {
      status=NEMO(input, maskFile, grid_variables, variables, rootname, bathymetry, nRequestedProcs, "C", maxmodes,verbose);
//       status=NEMO_FV(input, variable, rootname, bathymetry, nRequestedProcs, "C", maxmodes,verbose);
      status=VerticalModes_normalize(rootname, maxmodes, debug);
      }
    
    for(i=0;i<4;i++){
      string *filenamei=&filenames[i];
      if(*filenamei==""){
        *filenamei=varnames[i];
        if(i<4)
          *filenamei+="-atlas";
        *filenamei+=".nc";
        }
      }
//     filenames[0]="vovecrtz-atlas.nc";
//     filenames[1]="vozocrtx-atlas.nc";
//     filenames[2]="vomecrty-atlas.nc";
//     filenames[3]="Pbc.nc";
    
    if(decomposition) status=VerticalModes_SpectralDecomposition_NEMO(rootname, maskFile, maxmodes, wave, filenames, varnames, nrecomposed, debug, recomposition);
    if(status==0) {
      status=VerticalModes_SpectralDiagnostics(maskFile, ModesFile, maxmodes, wave, nrecomposed, verbose, debug);
      status=VerticalModes_SpectralEnergy(maskFile, ModesFile, maxmodes, wave, nrecomposed, verbose, debug);
      }
    }
  else if(model=="SYMPHONIE") {
    if(compute_modes) status=SYMPHONIE(input, rootname, bathymetry, nRequestedProcs, "", maxmodes, verbose);
    if(input[U_ID]!="") status=VerticalModes_SnapshotDecomposition(input, "SYMPHONIE.vertical-modes.nc", variables, maxmodes, debug);
    }
  else if(model=="SYMPHONIE-spectral") {
    string ModesFile=(string) rootname+".vertical-modes.nc";
    status=0;
    if(compute_modes) {
      if(compute_modes) status=SYMPHONIE_spectral(input, rootname, varnames, bathymetry, nRequestedProcs, "", maxmodes, verbose);
      status=VerticalModes_normalize(rootname, maxmodes, debug);
      }
    if(decomposition) status=VerticalModes_SpectralDecomposition(rootname, maxmodes, wave, varnames, nrecomposed, debug);
    const string maskFile="mesh_mask.nc";
    if(status==0) status=VerticalModes_SpectralEnergy(maskFile, ModesFile, maxmodes, wave, nrecomposed, verbose, debug);
    }
  else if(model=="ORCA") {
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "C", maxmodes,verbose);
    }
  else if(model=="ORCA12") {
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "C", maxmodes,verbose);
    }
  else if(model=="ORCA12-monthly") {
    
//     ORCA12.L46-MAL105b_y1998-2007mMM_TS.nc
    
    rootname=strdup("ORCA12-01");
    input[SAL_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m01_TS.nc");
    input[TEM_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m01_TS.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "C", maxmodes, verbose);
    rootname=strdup("ORCA12-02");
    input[SAL_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m02_TS.nc");
    input[TEM_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m02_TS.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "C", maxmodes, verbose);
    rootname=strdup("ORCA12-03");
    input[SAL_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m03_TS.nc");
    input[TEM_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m03_TS.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "C", maxmodes, verbose);
    rootname=strdup("ORCA12-04");
    input[SAL_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m04_TS.nc");
    input[TEM_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m04_TS.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "C", maxmodes, verbose);
    rootname=strdup("ORCA12-05");
    input[SAL_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m05_TS.nc");
    input[TEM_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m05_TS.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "C", maxmodes, verbose);
    rootname=strdup("ORCA12-06");
    input[SAL_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m06_TS.nc");
    input[TEM_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m06_TS.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "C", maxmodes, verbose);
    rootname=strdup("ORCA12-07");
    input[SAL_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m07_TS.nc");
    input[TEM_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m07_TS.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "C", maxmodes, verbose);
    rootname=strdup("ORCA12-08");
    input[SAL_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m08_TS.nc");
    input[TEM_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m08_TS.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "C", maxmodes, verbose);
    rootname=strdup("ORCA12-09");
    input[SAL_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m09_TS.nc");
    input[TEM_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m09_TS.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "C", maxmodes, verbose);
    rootname=strdup("ORCA12-10");
    input[SAL_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m10_TS.nc");
    input[TEM_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m10_TS.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "C", maxmodes, verbose);
    rootname=strdup("ORCA12-11");
    input[SAL_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m11_TS.nc");
    input[TEM_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m11_TS.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "C", maxmodes, verbose);
    rootname=strdup("ORCA12-12");
    input[SAL_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m12_TS.nc");
    input[TEM_ID]=strdup("./ORCA12.L46-MAL105b_y1998-2007m12_TS.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "C", maxmodes, verbose);
    }
  else if(model=="GLORYS") {
    /*status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "K", maxmodes, verbose);
    rootname=strdup("GLORYS-01");
    input[SAL_ID]=strdup("./GLORYS_y2001-2009_Sa_MONTH_01.nc");
    input[TEM_ID]=strdup("./GLORYS_y2001-2009_T_MONTH_01.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "K", maxmodes, verbose);
    rootname=strdup("GLORYS-02");
    input[SAL_ID]=strdup("./GLORYS_y2001-2009_Sa_MONTH_02.nc");
    input[TEM_ID]=strdup("./GLORYS_y2001-2009_T_MONTH_02.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "K", maxmodes, verbose);
    rootname=strdup("GLORYS-03");
    input[SAL_ID]=strdup("./GLORYS_y2001-2009_Sa_MONTH_03.nc");
    input[TEM_ID]=strdup("./GLORYS_y2001-2009_T_MONTH_03.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "K", maxmodes, verbose);
    rootname=strdup("GLORYS-04");
    input[SAL_ID]=strdup("./GLORYS_y2001-2009_Sa_MONTH_04.nc");
    input[TEM_ID]=strdup("./GLORYS_y2001-2009_T_MONTH_04.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "K", maxmodes, verbose);
    rootname=strdup("GLORYS-05");
    input[SAL_ID]=strdup("./GLORYS_y2001-2009_Sa_MONTH_05.nc");
    input[TEM_ID]=strdup("./GLORYS_y2001-2009_T_MONTH_05.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "K", maxmodes, verbose);
    rootname=strdup("GLORYS-06");
    input[SAL_ID]=strdup("./GLORYS_y2001-2009_Sa_MONTH_06.nc");
    input[TEM_ID]=strdup("./GLORYS_y2001-2009_T_MONTH_06.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "K", maxmodes, verbose);
    rootname=strdup("GLORYS-07");
    input[SAL_ID]=strdup("./GLORYS_y2001-2009_Sa_MONTH_07.nc");
    input[TEM_ID]=strdup("./GLORYS_y2001-2009_T_MONTH_07.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "K", maxmodes, verbose);
    rootname=strdup("GLORYS-08");
    input[SAL_ID]=strdup("./GLORYS_y2001-2009_Sa_MONTH_08.nc");
    input[TEM_ID]=strdup("./GLORYS_y2001-2009_T_MONTH_08.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "K", maxmodes, verbose);*/
    rootname=strdup("GLORYS-09");
    input[SAL_ID]=strdup("./GLORYS_y2001-2009_Sa_MONTH_09.nc");
    input[TEM_ID]=strdup("./GLORYS_y2001-2009_T_MONTH_09.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "K", maxmodes, verbose);
    rootname=strdup("GLORYS-10");
    input[SAL_ID]=strdup("./GLORYS_y2001-2009_Sa_MONTH_10.nc");
    input[TEM_ID]=strdup("./GLORYS_y2001-2009_T_MONTH_10.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "K", maxmodes, verbose);
    rootname=strdup("GLORYS-11");
    input[SAL_ID]=strdup("./GLORYS_y2001-2009_Sa_MONTH_11.nc");
    input[TEM_ID]=strdup("./GLORYS_y2001-2009_T_MONTH_11.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "K", maxmodes, verbose);
    rootname=strdup("GLORYS-12");
    input[SAL_ID]=strdup("./GLORYS_y2001-2009_Sa_MONTH_12.nc");
    input[TEM_ID]=strdup("./GLORYS_y2001-2009_T_MONTH_12.nc");
    status=ORCA_4(input, rootname, bathymetry, nRequestedProcs, "K", maxmodes, verbose);
    }
  else if(model=="WOA2005") {
//     if(input[TEM_ID]=="") input[TEM_ID]="/home/softs/data/climatology/WOA2005/temp_ann.cdf";
//     if(input[SAL_ID]=="") input[SAL_ID]="/home/softs/data/climatology/WOA2005/salt_ann.cdf";
    if(input[TEM_ID]=="") input[TEM_ID]="./temp_ann.cdf";
    if(input[SAL_ID]=="") input[SAL_ID]="./salt_ann.cdf";
    if(variables[0]=="") variables[0]="temperature";
    if(variables[1]=="") variables[1]="salinity";
//    status=WOA_annual(input, variable, model, bathymetry, nRequestedProcs, maxmodes);
    for (int frame=0; frame <12; frame++) {
      char tmp[1024];
     // sprintf(tmp,"%s-%2.2d",rootname,frame+1);
//       if(input[TEM_ID]=="") input[TEM_ID]="/home/softs/data/climatology/WOA2005/temp_month.cdf";
//       if(input[SAL_ID]=="") input[SAL_ID]="/home/softs/data/climatology/WOA2005/salt_month.cdf";
      if(input[TEM_ID]=="") input[TEM_ID]="./temp_month.cdf";
      if(input[SAL_ID]=="") input[SAL_ID]="./salt_month.cdf";
      status=WOA_monthly(input, variables, tmp, frame, bathymetry, nRequestedProcs, maxmodes, verbose);
      }
    }
  else if(model=="WOA2009") {
//     input[TEM_ID]=strdup("/home/softs/data/climatology/WOA2009/temperature_annual_1deg.nc");
//     input[SAL_ID]=strdup("/home/softs/data/climatology/WOA2009/salinity_annual_1deg.nc");
    if(input[TEM_ID]=="") input[TEM_ID]="./temperature_annual_1deg.nc";
    if(input[SAL_ID]=="") input[SAL_ID]="./salinity_annual_1deg.nc";
    if(variables[0]=="") variables[0]="t_an";
    if(variables[1]=="") variables[1]="s_an";
    status=WOA_annual(input, variables, rootname, bathymetry, nRequestedProcs, maxmodes,verbose);
    for (int frame=0; frame <12; frame++) {
      char tmp[1024];
      sprintf(tmp,"%s-%2.2d",rootname,frame+1);
//       input[TEM_ID]=strdup("/home/softs/data/climatology/WOA2009/temperature_monthly_1deg.nc");
//       input[SAL_ID]=strdup("/home/softs/data/climatology/WOA2009/salinity_monthly_1deg.nc");
      if(input[TEM_ID]=="") input[TEM_ID]="./temperature_monthly_1deg.nc";
      if(input[SAL_ID]=="") input[SAL_ID]="./salinity_monthly_1deg.nc";
      if(variables[0]=="") variables[0]="t_an";
      if(variables[1]=="") variables[1]="s_an";
      status=WOA_monthly(input, variables, tmp, frame, bathymetry, nRequestedProcs, maxmodes,verbose);
      }
    }
  else if(model=="ECCO") {
//     input[TEM_ID]=strdup("/home/softs/data/climatology/WOA2009/temperature_annual_1deg.nc");
//     input[SAL_ID]=strdup("/home/softs/data/climatology/WOA2009/salinity_annual_1deg.nc");
    //if(input[TEM_ID]=="") input[TEM_ID]="./temperature_annual_1deg.nc";
    //if(input[SAL_ID]=="") input[SAL_ID]="./salinity_annual_1deg.nc";
    status=ECCO_all(input, rootname, bathymetry, nRequestedProcs, maxmodes,verbose);
    /*rootname=strdup("ECCO-01");
    input[SAL_ID]=strdup("./ECCO_SA_MONTH_01.nc");
    input[TEM_ID]=strdup("./ECCO_T_MONTH_01.nc");
    status=ECCO_all(input, variable, rootname, bathymetry, nRequestedProcs, maxmodes);
    rootname=strdup("ECCO-02");
    input[SAL_ID]=strdup("./ECCO_SA_MONTH_02.nc");
    input[TEM_ID]=strdup("./ECCO_T_MONTH_02.nc");
    status=ECCO_all(input, variable, rootname, bathymetry, nRequestedProcs, maxmodes);
    rootname=strdup("ECCO-03");
    input[SAL_ID]=strdup("./ECCO_SA_MONTH_03.nc");
    input[TEM_ID]=strdup("./ECCO_T_MONTH_03.nc");
    status=ECCO_all(input, variable, rootname, bathymetry, nRequestedProcs, maxmodes);
    rootname=strdup("ECCO-04");
    input[SAL_ID]=strdup("./ECCO_SA_MONTH_04.nc");
    input[TEM_ID]=strdup("./ECCO_T_MONTH_04.nc");
    status=ECCO_all(input, variable, rootname, bathymetry, nRequestedProcs, maxmodes);
    rootname=strdup("ECCO-05");
    input[SAL_ID]=strdup("./ECCO_SA_MONTH_05.nc");
    input[TEM_ID]=strdup("./ECCO_T_MONTH_05.nc");
    status=ECCO_all(input, variable, rootname, bathymetry, nRequestedProcs, maxmodes);
    rootname=strdup("ECCO-06");
    input[SAL_ID]=strdup("./ECCO_SA_MONTH_06.nc");
    input[TEM_ID]=strdup("./ECCO_T_MONTH_06.nc");
    status=ECCO_all(input, variable, rootname, bathymetry, nRequestedProcs, maxmodes);
    rootname=strdup("ECCO-07");
    input[SAL_ID]=strdup("./ECCO_SA_MONTH_07.nc");
    input[TEM_ID]=strdup("./ECCO_T_MONTH_07.nc");
    status=ECCO_all(input, variable, rootname, bathymetry, nRequestedProcs, maxmodes);
    rootname=strdup("ECCO-08");
    input[SAL_ID]=strdup("./ECCO_SA_MONTH_08.nc");
    input[TEM_ID]=strdup("./ECCO_T_MONTH_08.nc");
    status=ECCO_all(input, variable, rootname, bathymetry, nRequestedProcs, maxmodes);
    rootname=strdup("ECCO-09");
    input[SAL_ID]=strdup("./ECCO_SA_MONTH_09.nc");
    input[TEM_ID]=strdup("./ECCO_T_MONTH_09.nc");
    status=ECCO_all(input, variable, rootname, bathymetry, nRequestedProcs, maxmodes);
    rootname=strdup("ECCO-10");
    input[SAL_ID]=strdup("./ECCO_SA_MONTH_10.nc");
    input[TEM_ID]=strdup("./ECCO_T_MONTH_10.nc");
    status=ECCO_all(input, variable, rootname, bathymetry, nRequestedProcs, maxmodes);
    rootname=strdup("ECCO-11");
    input[SAL_ID]=strdup("./ECCO_SA_MONTH_11.nc");
    input[TEM_ID]=strdup("./ECCO_T_MONTH_11.nc");
    status=ECCO_all(input, variable, rootname, bathymetry, nRequestedProcs, maxmodes);
    rootname=strdup("ECCO-12");
    input[SAL_ID]=strdup("./ECCO_SA_MONTH_12.nc");
    input[TEM_ID]=strdup("./ECCO_T_MONTH_12.nc");
    status=ECCO_all(input, variable, rootname, bathymetry, nRequestedProcs, maxmodes);*/
//     for (int frame=0; frame <12; frame++) {
//       char tmp[1024];
//       sprintf(input[TEM_ID],"%s_%2.2d.nc","./ECCO_SA_MONTH",frame+1);
//       sprintf(input[SAL_ID],"%s_%2.2d.nc","./ECCO_T_MONTH",frame+1);
//       status=ECCO_monthly(input, variable, tmp, frame, bathymetry, nRequestedProcs, maxmodes);
//       }
    }
  else if(model=="LEVITUS") {
    status=levitus(input, rootname, bathymetry, nRequestedProcs, maxmodes, verbose);
    }
  else{
    printf("*** please specify model with -m ***\n");
    print_help(argv[0]);
    wexit(-1);
    }
  
  TRAP_ERR_EXIT(0,"%s -computation complete ^^^^^^^^^\n",argv[0]);
}
