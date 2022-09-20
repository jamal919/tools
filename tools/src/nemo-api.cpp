
/*******************************************************************************
T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative
  
E-mail: florent.lyard@legos.obs-mip.fr
*******************************************************************************/
/**
\file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada
\author  Yves Soufflet      LEGOS, Toulouse, France
**/


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

extern  void compute_atlas_names(const string & wave,const string & varname,string *filename,string *varnameA,string *varnameG);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int NEMO_u2p(const grid_t ugrid, bool *u_mask, double *e3u, const grid_t tgrid, bool *t_mask, double *e3t, complex<double> *u, complex<double> **out, complex<double> mask, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k;
  complex<double> u1,u2;
  
  for(int j=0;j<tgrid.ny;j++) {
    for(i=0;i<tgrid.nx;i++) {
      size_t n=j*tgrid.nx+i;
      for(k=0;k<tgrid.nz;k++) {
        out[n][k]=mask;
        }
      }
    }
 
  for(int j=0;j<tgrid.ny;j++) {
    for(i=2;i<tgrid.nx-2;i++) {
      size_t n=j*tgrid.nx+i;
      for(k=0;k<tgrid.nz;k++) {
/*------------------------------------------------------------------------------
        t-node ndex */
        size_t m=k*tgrid.Hsize()+j*tgrid.nx+i;
        if(t_mask[m]==0) continue;
/*------------------------------------------------------------------------------
        u-node index */
        size_t m1=k*ugrid.Hsize()+j*ugrid.nx+i-1;
        size_t m2=k*ugrid.Hsize()+j*ugrid.nx+i;
        if(u_mask[m1]==0 and t_mask[m]==1) u1=0;
        else u1=u[m1];
        if(u_mask[m2]==0 and t_mask[m]==1) u2=0;
        else u2=u[m2];
//         if(u[m1]!=umask and u[m2]!=umask) {
          out[n][k]=0.5*(e3u[m1]*u[m1]+e3u[m2]*u[m2])/e3t[m];
//           }
//         else  {
/*------------------------------------------------------------------------------
          */
//           out[n][k]=mask;
//           }
        }
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int VerticalModes_SpectralDecomposition_NEMO(const string & rootname, const string maskFile, int maxmodes, string & wave, string filenames[4], string varnames[4], int nrecomposed, bool debug, bool recomposition)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Interface to decompose harmonic fields into modal contributions
  
  Presently hard-coded for NEMO-derived files, badly-formed structured
  ouputs are just a pain in the ... as usual
  
    v f v f v f
    t u t u t u   last T line unusable as v (above) masked
    v f v f v f
    t u t u t u   first T line unusable as v start after
    
  
  u,v and t have same dimensions: 2 lines and columns too much in t

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i,k,n,status;
  grid_t grid,ugrid,vgrid,wgrid;
  double **Umodes=0, **Pmodes=0, **Wmodes=0, mask;
  int MinUmodes;
  complex<double> **u=0, **v=0, **w=0, **p=0, umask, *buffer=0, cmask, *eta=0;
  complex<double> **u_decomposition=0, **v_decomposition=0, **p_decomposition=0, **w_decomposition=0;
  complex<double> *ubarotropic=0, *ubaroclinic=0, *vbarotropic=0, *vbaroclinic=0, *wbarotropic=0, *wbaroclinic=0, *pbarotropic=0, *pbaroclinic=0;
  string filename,varnameA,varnameG,varname,output_decomposition,output_recomposition;
  int nlayers,nlevels;
  int verbose=0;
  poc_data_t<float> scalarData;
  poc_data_t<int8_t> t_mask,u_mask,v_mask,w_mask;
  int range[2]={0,-1};
  int create_file=1, create_grid=1;
  int nkept;
  poc_data_t<double> e3w,e3t,e3u,e3v;/*< for partial steps */
  double g=9.81, rho0=1000.0;
   
  nkept=min(maxmodes, nrecomposed);
  
  create_file=1;
  output_decomposition=wave+"-modal-decomposition.nc";
  output_recomposition=wave+"-modal-recomposition.nc";

//   goto pmodes;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  process vertical velocity modal decomposition
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
/*------------------------------------------------------------------------------
  hard coded */
  
  status=w_mask.init(maskFile,"tmask");
  status=w_mask.read_data(maskFile,-1,1);

  filename=wave+"-"+filenames[0];
  
  varnameA=varnames[0]+"_a";
  varnameG=varnames[0]+"_G";
  
  status=poc_get_grid(filename, varnameA, &grid, verbose, -1);
  if(status!=0) {
/*------------------------------------------------------------------------------
    get tracer grid from pressure field; not safe*/
    filename=wave+"-Pbc.nc";
    varnameA=varnames[3]+"_a";
    varnameG=varnames[3]+"_G";
    status=poc_get_grid(filename, varnameA, &grid, verbose, -1);
    if(status==0) goto umodes;
    else goto finish;
    }
    
  if(grid.ny==1) {
    grid.y[0]=49.8;
    if(grid.modeH==1) status=map_completegridaxis(&grid,2);
    }

  nlevels=grid.nz+1;
  
  buffer=new complex<double>[grid.size()];
  
  status=poc_get_cvara(filename, varnameA, varnameG,-1, buffer);
  
  umask=NC_FILL_COMPLEX;
  w=new complex<double>* [grid.Hsize()];
  for(n=0;n<grid.Hsize();n++) {
    w[n]=new complex<double>[nlevels];
    for(k=0;k<nlevels;k++) {
      w[n][k]=umask;
      }
    }
  
  for(n=0;n<grid.Hsize();n++) {
    for(k=0;k<nlevels-1;k++) {
      size_t m=k*grid.Hsize()+n;
      if(w_mask.data[m]==0) continue;
      if(buffer[m]!=umask) {
        w[n][k]=buffer[m];
        }
      else  {
        w[n][k]=umask;
        }
      }
    }
  delete[] buffer;
  
/*------------------------------------------------------------------------------
  add 0 at bottom level */
  for(n=0;n<grid.Hsize();n++) {
    if(w[n][0]==umask) continue;
    for(k=0;k<nlevels-1;k++) {
//       size_t m=k*grid.Hsize()+n;
      if(w[n][k]==umask) {
        w[n][k]=0.0;
        break;
        }
      }
    }
  
  filename=rootname+".vertical-modes.nc";
  varname="Wmodes";
  status=VerticalModes_SpectralDecomposition(filename, varname, maxmodes, Wmodes, mask, w, umask, w_decomposition, MinUmodes, debug);
//   grid.free();
 
  grid.nz++;
  
/*------------------------------------------------------------------------------
  store modal decomposition coefficients */
  range[0]=0;
  range[1]=MinUmodes-1;
  create_grid=1;
  status=save_SGXYT_C(output_decomposition.c_str(), grid, w_decomposition, umask, "Wc_a", "Wc_G", "dimensionless", "w_coeffcients", create_file, create_grid, range);

/*------------------------------------------------------------------------------
  re-compose vertical velocity modes */
  if(recomposition){
    wbarotropic=new complex<double>[grid.size()];
    for(size_t m=0;m<grid.size();m++) wbarotropic[m]=umask;
    wbaroclinic=new complex<double>[grid.size()];
    for(size_t m=0;m<grid.size();m++) wbaroclinic[m]=umask;
    
    for(i=0;i<nkept;i++) {
      for(n=0;n<grid.Hsize();n++) {
        for(int k=0;k<nlevels;k++) {
          size_t m=k*grid.Hsize()+n;
          wbaroclinic[m]=umask;
          if(Wmodes[i][m]==mask) continue;
          if(w_decomposition[i][n]!=umask) wbaroclinic[m]=Wmodes[i][m]*w_decomposition[i][n];
          }
        }
      create_file=(i==0);
      create_grid=(i==0);
      int frame=i;
      status=save_SGXYZT_C(output_recomposition.c_str(), grid, wbaroclinic, umask, "w_modal_a","w_modal_G","m/s", "w_modal", create_file, create_grid,"w",frame);
      }
    
/*------------------------------------------------------------------------------
    re-compose vertical velocity barotropic and baroclinic modes */
    i=0;
    for(n=0;n<grid.Hsize();n++) {
      for(int k=0;k<nlevels;k++) {
        size_t m=k*grid.Hsize()+n;
        wbarotropic[m]=umask;
        if(Wmodes[i][m]==mask) continue;
        if(w_decomposition[i][n]!=umask) wbarotropic[m]=Wmodes[i][m]*w_decomposition[i][n];
        }
      }
    create_file=0;
    create_grid=0;
    status=save_SGXYZ_C(output_recomposition.c_str(), grid, wbarotropic, umask, "w_barotropic_a","w_barotropic_G","m/s", "w_barotropic", create_file, create_grid,"","w");
    
    i=1;
    for(n=0;n<grid.Hsize();n++) {
      for(int k=0;k<nlevels;k++) {
        size_t m=k*grid.Hsize()+n;
        wbaroclinic[m]=umask;
        if(Wmodes[i][m]==mask) continue;
        if(w_decomposition[i][n]!=umask) wbaroclinic[m]=Wmodes[i][m]*w_decomposition[i][n];
        }
      }
    for(i=2;i<MinUmodes-1;i++) {
      for(n=0;n<grid.Hsize();n++) {
        for(int k=0;k<nlevels;k++) {
          size_t m=k*grid.Hsize()+n;
          if(wbaroclinic[m]==umask) continue;
          if(Wmodes[i][m]==mask) continue;
          if(w_decomposition[i][n]!=umask) wbaroclinic[m]+=Wmodes[i][m]*w_decomposition[i][n];
          }
        }
      }
    create_file=0;
    create_grid=0;
    status=save_SGXYZ_C(output_recomposition.c_str(), grid, wbaroclinic, umask, "w_baroclinic_a","w_baroclinic_G","m/s", "w_baroclinic", create_file, create_grid,"","w");
    }
  
  grid.nz--;
  create_file=0;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  process horizontal velocity modal decomposition
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
umodes:
  nlayers=grid.nz;
  
  umask=NC_FILL_COMPLEX;
  u=new complex<double>* [grid.Hsize()];
  v=new complex<double>* [grid.Hsize()];
  for(n=0;n<grid.Hsize();n++) {
    u[n]=new complex<double>[nlayers];
    v[n]=new complex<double>[nlayers];
    for(k=0;k<nlayers;k++) {
      u[n][k]=umask;
      v[n][k]=umask;
      }
    }
  
  status=t_mask.init(maskFile,"tmask");
  status=u_mask.init(maskFile,"umask");
  status=v_mask.init(maskFile,"vmask");
  
  status=t_mask.read_data(maskFile,-1,1);
  status=u_mask.read_data(maskFile,-1,1);
  status=v_mask.read_data(maskFile,-1,1);
  
  status=e3t.init(maskFile,"e3t");
  status=e3u.init(maskFile,"e3u");
  status=e3v.init(maskFile,"e3v");
  
  status=e3t.read_data(maskFile,-1,1);
  status=e3u.read_data(maskFile,-1,1);
  status=e3v.read_data(maskFile,-1,1);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  zonal velocity, load and interpolate at tracer nodes
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
//   filename=wave+"-"+varnames[1]+"-atlas.nc";
  if(filenames[1]=="") goto pmodes;
  
  filename=wave+"-"+filenames[1];

//   varnameA="vozocrtx_a";
//   varnameG="vozocrtx_G";
  varnameA=varnames[1]+"_a";
  varnameG=varnames[1]+"_G";
  
  status=poc_get_grid(filename, varnameA, &ugrid, verbose, -1);
  if(status!=0) NC_TRAP_ERROR(wexit,status,1,"poc_get_grid error with "+varnameA+" in "+filename+"\n");

  if(grid.ny==1) grid.y[0]=49.8;
  
  buffer=new complex<double>[ugrid.size()];
  
  status=poc_get_cvara(filename,varnameA, varnameG,-1, buffer);
  
/*------------------------------------------------------------------------------
  zonal velocity, interpolate at tracer nodes */
//   for(int j=2;j<grid.ny-2;j++) {
  for(int j=0;j<grid.ny;j++) {
    for(i=2;i<grid.nx-2;i++) {
      size_t n=j*grid.nx+i;
      for(k=0;k<nlayers;k++) {
/*------------------------------------------------------------------------------
        t-node index */
        size_t m=k*grid.Hsize()+j*grid.nx+i;
        if(t_mask.data[m]==0) continue;
/*------------------------------------------------------------------------------
        u-node index */
        size_t m1=k*ugrid.Hsize()+j*ugrid.nx+i-1;
        size_t m2=k*ugrid.Hsize()+j*ugrid.nx+i;
        if(u_mask.data[m1]==0 and t_mask.data[m]==1) buffer[m1]=0;
        if(u_mask.data[m2]==0 and t_mask.data[m]==1) buffer[m2]=0;
        if(buffer[m1]!=umask and buffer[m2]!=umask) {
/*------------------------------------------------------------------------------
          too rough, dismissed */
//           u[n][k]=0.5*(buffer[m1]+buffer[m2]);
/*------------------------------------------------------------------------------
          Rachid hint:
          Ut(i,j,k) =0.5*( e3u(i,j,k)*u(i,j,k) +  e3u(i+1,j,k)*u(i+1,j,k) ) / e3t(i,j,k) */
          u[n][k]=0.5*(e3u.data[m1]*buffer[m1]+e3u.data[m2]*buffer[m2])/e3t.data[m];
          }
        else  {
/*------------------------------------------------------------------------------
          ??? */
//           u[n][k]=umask;
          u[n][k]=0.0;
          }
        }
      }
    }

  delete[] buffer;
  ugrid.free();
  

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  meridional velocity, load and interpolate at tracer nodes
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
//   filename=wave+"-"+varnames[2]+"-atlas.nc";
  if(filenames[2]=="") goto pmodes;
  
  filename=wave+"-"+filenames[2];
  
//   varnameA="vomecrty_a";
//   varnameG="vomecrty_G";
  varnameA=varnames[2]+"_a";
  varnameG=varnames[2]+"_G";
  
  status=poc_get_grid(filename, varnameA, &vgrid, verbose, -1);
  if(status!=0) NC_TRAP_ERROR(wexit,status,1,"poc_get_grid error with "+varnameA+" in "+filename+"\n");

  if(grid.ny==1) grid.y[0]=49.8;
  
  buffer=new complex<double>[vgrid.size()];
  
  status=poc_get_cvara(filename,varnameA, varnameG,-1, buffer);
  
/*------------------------------------------------------------------------------
  zonal velocity, interpolate at tracer nodes */
  if(grid.ny!=1) {
/*------------------------------------------------------------------------------
    regular grid */
    for(int j=2;j<grid.ny-2;j++) {
      for(int i=2;i<grid.nx-2;i++) {
        size_t n=j*grid.nx+i;
        for(k=0;k<nlayers;k++) {
          size_t m=k*grid.Hsize()+j*grid.nx+i;
          if(t_mask.data[m]==0) continue;
          size_t m1=k*vgrid.Hsize()+(j-1)*vgrid.nx+i;
          size_t m2=k*vgrid.Hsize()+j*vgrid.nx+i;
          if(v_mask.data[m1]==0 and t_mask.data[m]==1) buffer[m1]=0;
          if(v_mask.data[m2]==0 and t_mask.data[m]==1) buffer[m2]=0;
          if(buffer[m1]!=umask and buffer[m2]!=umask) {
/*------------------------------------------------------------------------------
            too rough, dismissed */
//             v[n][k]=0.5*(buffer[m1]+buffer[m2]);
/*------------------------------------------------------------------------------
            Rachid hint:
            Ut(i,j,k) =0.5*( e3u(i,jk)*u(i,j,k) +  e3u(i+1,j,k)*u(i+1,j,k) ) / e3t(i,j,k) */
            v[n][k]=0.5*(e3v.data[m1]*buffer[m1]+e3v.data[m2]*buffer[m2])/e3t.data[m];
            }
          else  {
            v[n][k]=0;
//             v[n][k]=umask;
            }
          }
        }
      }
    }
  else {
/*------------------------------------------------------------------------------
    y-periodic grid (COMODO IW test case) */
    int j=0;
    for(int i=2;i<grid.nx-2;i++) {
      size_t n=j*grid.nx+i;
      for(k=0;k<nlayers;k++) {
        size_t m=k*grid.Hsize()+j*grid.nx+i;
        if(t_mask.data[m]==0) continue;
        size_t m1=k*vgrid.Hsize()+j*vgrid.nx+i;
        size_t m2=k*vgrid.Hsize()+j*vgrid.nx+i;
        if(v_mask.data[m1]==0 and t_mask.data[m]==1) buffer[m1]=0;
        if(v_mask.data[m2]==0 and t_mask.data[m]==1) buffer[m2]=0;
        if(buffer[m1]!=umask and buffer[m2]!=umask) {
          v[n][k]=0.5*(buffer[m1]+buffer[m2]);
          }
        else  {
          v[n][k]=0;
//           v[n][k]=umask;
          }
        }
      }
    }
  delete[] buffer;
  vgrid.free();
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  compute u,v modal decomposition
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  filename=rootname+".vertical-modes.nc";
  varname="Umodes";
  
  status=VerticalModes_SpectralDecomposition(filename, varname, maxmodes, Umodes, mask, u, umask, u_decomposition, MinUmodes, debug);
  status=VerticalModes_SpectralDecomposition(filename, varname, maxmodes, Umodes, mask, v, umask, v_decomposition, MinUmodes, debug);
  
  create_grid=1;
  range[0]=0;
  range[1]=MinUmodes-1;
/*------------------------------------------------------------------------------
  modal decomposition coefficients */
  status=save_SGXYT_C(output_decomposition.c_str(), grid, u_decomposition, umask, "Uc_a", "Uc_G", "dimensionless", "u_coeffcients", create_file, create_grid, range);
  create_file=0;
  create_grid=1;
  status=save_SGXYT_C(output_decomposition.c_str(), grid, v_decomposition, umask, "Vc_a", "Vc_G", "dimensionless", "u_coeffcients", create_file, create_grid, range);
  
if(recomposition){
  ubarotropic=new complex<double>[grid.size()];
  for(size_t m=0;m<grid.size();m++) ubarotropic[m]=umask;
  ubaroclinic=new complex<double>[grid.size()];
  for(size_t m=0;m<grid.size();m++) ubaroclinic[m]=umask;
  
  vbarotropic=new complex<double>[grid.size()];
  for(size_t m=0;m<grid.size();m++) vbarotropic[m]=umask;
  vbaroclinic=new complex<double>[grid.size()];
  for(size_t m=0;m<grid.size();m++) vbaroclinic[m]=umask;
  
//   create_file=0;
//   create_grid=1;
// /*------------------------------------------------------------------------------
//   modal fields */
//   status=save_SGXYZ_C(output_decomposition.c_str(), grid, ubarotropic, umask, "Ubt_a", "Ubt_G", "m/s", "u_barotropic", create_file, create_grid);
//   create_file=0;
//   status=save_SGXYZ_C(output_decomposition.c_str(), grid, ubaroclinic, umask, "Ubc_a", "Ubc_G", "m/s", "u_baroclinic", create_file, create_grid);
//   status=save_SGXYZ_C(output_decomposition.c_str(), grid, vbarotropic, umask, "Vbt_a", "Vbt_G", "m/s", "v_barotropic", create_file, create_grid);
//   status=save_SGXYZ_C(output_decomposition.c_str(), grid, vbaroclinic, umask, "Vbc_a", "Vbc_G", "m/s", "v_baroclinic", create_file, create_grid);
  
  for(i=0;i<nkept;i++) {
    for(n=0;n<grid.Hsize();n++) {
      for(int k=0;k<nlayers;k++) {
        size_t m=k*grid.Hsize()+n;
        ubaroclinic[m]=umask;
        vbaroclinic[m]=umask;
        if(Umodes[i][m]==mask) continue;
        if(u_decomposition[i][n]!=umask) ubaroclinic[m]=Umodes[i][m]*u_decomposition[i][n];
        if(v_decomposition[i][n]!=umask) vbaroclinic[m]=Umodes[i][m]*v_decomposition[i][n];
        }
      }
    create_file=0;
    create_grid=(i==0);
    int frame=i;
    status=save_SGXYZT_C(output_recomposition.c_str(), grid, ubaroclinic, umask, "u_modal_a","u_modal_G","m/s", "u_modal", create_file, create_grid,"",frame);
    create_file=0;
    create_grid=0;
    status=save_SGXYZT_C(output_recomposition.c_str(), grid, vbaroclinic, umask, "v_modal_a","v_modal_G","m/s", "v_modal", create_file, create_grid,"",frame);
    }

/*------------------------------------------------------------------------------
  re-compose horizontal velocity barotropic and baroclinic modes */
  i=0;
  for(n=0;n<grid.Hsize();n++) {
    for(int k=0;k<nlayers;k++) {
      size_t m=k*grid.Hsize()+n;
      ubarotropic[m]=umask;
      vbarotropic[m]=umask;
      if(Umodes[i][m]==mask) continue;
      if(u_decomposition[i][n]!=umask) ubarotropic[m]=Umodes[i][m]*u_decomposition[i][n];
      if(v_decomposition[i][n]!=umask) vbarotropic[m]=Umodes[i][m]*v_decomposition[i][n];
      }
    }
  create_file=0;
  create_grid=0;
  status=save_SGXYZ_C(output_recomposition.c_str(), grid, ubarotropic, umask, "u_barotropic_a","u_barotropic_G","m/s", "u_barotropic", create_file, create_grid,"","");
  status=save_SGXYZ_C(output_recomposition.c_str(), grid, vbarotropic, umask, "v_barotropic_a","v_barotropic_G","m/s", "v_barotropic", create_file, create_grid,"","");
  
  i=1;
  for(n=0;n<grid.Hsize();n++) {
    for(int k=0;k<nlayers;k++) {
      size_t m=k*grid.Hsize()+n;
      ubaroclinic[m]=umask;
      vbaroclinic[m]=umask;
      if(Umodes[i][m]==mask) continue;
      if(u_decomposition[i][n]!=umask) ubaroclinic[m]=Umodes[i][m]*u_decomposition[i][n];
      if(v_decomposition[i][n]!=umask) vbaroclinic[m]=Umodes[i][m]*v_decomposition[i][n];
      }
    }
  for(i=2;i<MinUmodes-1;i++) {
    for(n=0;n<grid.Hsize();n++) {
      for(int k=0;k<nlayers;k++) {
        size_t m=k*grid.Hsize()+n;
        if(ubaroclinic[m]==umask) continue;
        if(Umodes[i][m]==mask) continue;
        if(u_decomposition[i][n]!=umask) ubaroclinic[m]+=Umodes[i][m]*u_decomposition[i][n];
        if(v_decomposition[i][n]!=umask) vbaroclinic[m]+=Umodes[i][m]*v_decomposition[i][n];
        }
      }
    }
  create_file=0;
  create_grid=0;
  status=save_SGXYZ_C(output_recomposition.c_str(), grid, ubaroclinic, umask, "u_baroclinic_a","u_baroclinic_G","m/s", "u_baroclinic", create_file, create_grid,"","");
  status=save_SGXYZ_C(output_recomposition.c_str(), grid, vbaroclinic, umask, "v_baroclinic_a","v_baroclinic_G","m/s", "v_baroclinic", create_file, create_grid,"","");

  create_file=0;
  }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  process pressure modal decomposition
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
pmodes:
  filename=wave+"-"+filenames[3];
  
  varnameA=varnames[3]+"_a";
  varnameG=varnames[3]+"_G";
  
//   status=poc_get_grid(filename, varnameA, &grid, verbose, -1);
//   nlayers=grid.nz;
//   if(t_mask.data==0) {
//     status=t_mask.init(maskFile,"tmask");
//     status=t_mask.read_data(maskFile,-1,1);
//     }
  
  buffer=new complex<double>[grid.size()];
  
  umask=NC_FILL_COMPLEX;
  p=new complex<double>* [grid.Hsize()];
  for(n=0;n<grid.Hsize();n++) {
    p[n]=new complex<double>[nlayers];
    for(k=0;k<nlayers;k++) {
      p[n][k]=umask;
      }
    }
  
  status=poc_get_cvara(filename,varnameA, varnameG,-1, buffer);
  if(status!=0) goto finish;
  
#ifdef SSH_CORRECTION
  complex<double> *ssh;
  ssh=new complex<double>[grid.Hsize()];
  compute_atlas_names(wave,varnames[4],&filename,&varnameA,&varnameG);
  status=poc_get_cvara(filename,varnameA, varnameG,0, ssh);
  NC_TRAP_ERROR(wexit,status,1,"poc_get_cvara(\""+filename+"\",\""+varnameA+"\",\""+varnameG+"\",-1,) error");
  status=e3w.init(maskFile,"e3w");
  status=e3w.read_data(maskFile,-1,1);
#endif
  
//   status=e3t.init(maskFile,"e3t");
//   status=e3t.read_data(maskFile,-1,1);
  
  umask=NC_FILL_COMPLEX;
    
  for(n=0;n<grid.Hsize();n++) {
    if(t_mask.data[n]==0) continue;
#ifdef SSH_CORRECTION
    double depth=0;
    for(k=0;k<nlayers;k++) {
      size_t m=k*grid.Hsize()+n;
      if(t_mask.data[m]==0) continue;
      depth+=e3t.data[m];
      }
    double immersion=e3t.data[n]/2.;
    double theta=(depth-immersion)/depth;
    p[n][0]=buffer[n]+g*rho0*ssh[n]*theta;
    for(k=1;k<nlayers;k++) {
      size_t m=k*grid.Hsize()+n;
      if(t_mask.data[m]==0) continue;
      immersion+=e3w.data[m];
      theta=(depth-immersion)/depth;
      if(buffer[m]!=umask) {
//         p[n][k]=buffer[m];
        p[n][k]=buffer[m]+g*rho0*ssh[n]*theta;
        }
      else  {
        p[n][k]=umask;
        }
      }
#else
    for(k=0;k<nlayers;k++) {
      size_t m=k*grid.Hsize()+n;
      if(t_mask.data[m]==0) continue;
      if(buffer[m]!=umask) {
        p[n][k]=buffer[m];
        }
      else  {
        p[n][k]=umask;
        }
      }
#endif
    }
  delete[] buffer;
  
  pbarotropic=new complex<double>[grid.size()];
  for(size_t m=0;m<grid.size();m++) pbarotropic[m]=umask;
  pbaroclinic=new complex<double>[grid.size()];
  for(size_t m=0;m<grid.size();m++) pbaroclinic[m]=umask;
  
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  compute p modal decomposition
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  filename=rootname+".vertical-modes.nc";
  varname="Pmodes";
  
  status=VerticalModes_SpectralDecomposition(filename, varname, maxmodes, Pmodes, mask, p, umask, p_decomposition, MinUmodes, debug);
  
  create_file=0;
  create_grid=1;
  range[0]=0;
  range[1]=MinUmodes-1;
/*------------------------------------------------------------------------------
  modal decomposition coefficients */
  status=save_SGXYT_C(output_decomposition.c_str(), grid, p_decomposition, umask, "Pc_a", "Pc_G", "dimensionless", "p_coeffcients", create_file, create_grid, range);

//   create_file=0;
//   create_grid=1;
// /*------------------------------------------------------------------------------
//   modal fields */
//   status=save_SGXYZ_C(output_decomposition.c_str(), grid, pbarotropic, umask, "Pbt_a", "Pbt_G", "N/m²", "p_barotropic", create_file, create_grid);
//   create_file=0;
//   create_grid=1;
//   status=save_SGXYZ_C(output_decomposition.c_str(), grid, pbaroclinic, umask, "Pbc_a", "Pbc_G", "N/m²", "p_baroclinic", create_file, create_grid);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  compute p modal recomposition
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(recomposition){
    for(i=0;i<nkept;i++) {
      for(n=0;n<grid.Hsize();n++) {
        for(int k=0;k<nlayers;k++) {
          size_t m=k*grid.Hsize()+n;
          pbaroclinic[m]=umask;
          if(Pmodes[i][m]==mask) continue;
          if(p_decomposition[i][n]!=umask) pbaroclinic[m]=Pmodes[i][m]*p_decomposition[i][n];
          }
        }
      create_file=0;
      create_grid=0;
      int frame=i;
      status=save_SGXYZT_C(output_recomposition.c_str(), grid, pbaroclinic, umask, "p_modal_a","p_modal_G","N/m²", "p_modal", create_file, create_grid,"",frame);
      }
  
/*------------------------------------------------------------------------------
    re-compose pressure barotropic and baroclinic modes */
    i=0;
    for(n=0;n<grid.Hsize();n++) {
      for(int k=0;k<nlayers;k++) {
        size_t m=k*grid.Hsize()+n;
        pbarotropic[m]=umask;
        if(Pmodes[i][m]==mask) continue;
        if(p_decomposition[i][n]!=umask) pbarotropic[m]=Pmodes[i][m]*p_decomposition[i][n];
        }
      }
    create_file=0;
    create_grid=0;
    status=save_SGXYZ_C(output_recomposition.c_str(), grid, pbarotropic, umask, "p_barotropic_a","p_barotropic_G","m/s", "p_barotropic", create_file, create_grid,"","");
  
    eta=new complex<double> [grid.Hsize()];
    
    for(n=0;n<grid.Hsize();n++) {
      k=0;
      size_t m=k*grid.Hsize()+n;
      if(pbarotropic[m]!=umask) eta[n]=pbarotropic[m]/(double) 1020/ 9.81;
      else eta[n]=umask;
      }
    create_file=0;
    create_grid=0;
    status=save_SGXY_C(output_recomposition.c_str(), grid, eta, umask, "eta_barotropic_a","eta_barotropic_G","m/s", "p_barotropic", create_file, create_grid);
    
    i=1;
    for(n=0;n<grid.Hsize();n++) {
      for(int k=0;k<nlayers;k++) {
        size_t m=k*grid.Hsize()+n;
        pbaroclinic[m]=umask;
        if(Pmodes[i][m]==mask) continue;
        if(p_decomposition[i][n]!=umask) pbaroclinic[m]=Pmodes[i][m]*p_decomposition[i][n];
        }
      }
    for(i=2;i<MinUmodes-1;i++) {
      for(n=0;n<grid.Hsize();n++) {
        for(int k=0;k<nlayers;k++) {
          size_t m=k*grid.Hsize()+n;
          if(pbaroclinic[m]==umask) continue;
          if(Pmodes[i][m]==mask) continue;
          if(p_decomposition[i][n]!=umask) pbaroclinic[m]+=Pmodes[i][m]*p_decomposition[i][n];
          }
        }
      }
    create_file=0;
    create_grid=0;
    status=save_SGXYZ_C(output_recomposition.c_str(), grid, pbaroclinic, umask, "p_baroclinic_a","p_baroclinic_G","m/s", "p_baroclinic", create_file, create_grid,"","");
    for(n=0;n<grid.Hsize();n++) {
      k=0;
      size_t m=k*grid.Hsize()+n;
      if(pbaroclinic[m]!=umask) eta[n]=pbaroclinic[m]/(double) 1020/ 9.81;
      else eta[n]=umask;
      }
    create_file=0;
    create_grid=0;
    status=save_SGXY_C(output_recomposition.c_str(), grid, eta, umask, "eta_baroclinic_a","eta_baroclinic_G","m/s", "p_barotropic", create_file, create_grid);
    }

finish:

  status=SpectralEnergy_RAW(maskFile, grid, u, v, w, p, wave, verbose, debug);
  
  for(n=0;n<grid.Hsize();n++) {
    delete[] u[n];
    delete[] v[n];
    delete[] w[n];
    delete[] p[n];
    }
    
  delete[] u;
  delete[] v;
  delete[] w;
  delete[] p;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ORCA_4(const string *input,const char *rootname,const char *bathymetry, int nRequestedProcs, const char *units, int maxmodes, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int j,k,m;
  char filename[1024];
  grid_t grid, w_grid;
  float *R,*T,*S, mask;
  grid_t topogrid;

  if(input[TEM_ID]=="") sprintf(filename,"/home/data/climatology/ORCA/orca-4/ORCA025.L75-MJM95_y2000-2009_gridT.nc");
  else strcpy(filename,input[TEM_ID].c_str());
  
  status= map_loadfield3D(filename, filename, "votemper", grid, T, mask, verbose);
  if(status !=0) {
    STDOUT_BASE_LINE("cannot load %s climatology file=%s\n", "temperature", filename);
    exit(-1);
    }
  
  if(input[SAL_ID]=="") sprintf(filename,"/home/data/climatology/ORCA/orca-4/ORCA025.L75-MJM95_y2000-2009_gridT.nc");
  else strcpy(filename,input[SAL_ID].c_str());
  
  status= map_loadfield3D(filename, filename, "vosaline", grid, S, mask, verbose);
  if(status !=0) {
    STDOUT_BASE_LINE("cannot load %s climatology file=%s\n", "salinity", filename);
    exit(-1);
    }

    
  if(strcmp(units,"K")==0) {
    for(k=0;k<grid.nz;k++) {
      for(j=0;j<grid.ny;j++) {
        for(int i=0;i<grid.nx;i++) {
           m=k*grid.ny*grid.nx+j*grid.nx+i;
           if(T[m]!=mask) {
             T[m]-=273.15;
             }
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
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int NEMO(const string *input, const string maskFile, const string *grid_variables, const string *variables, const char *rootname,const char *bathymetry, int nRequestedProcs, const char *units, int maxmodes, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int j,k,m;
  string filename,varname;
  poc_var_t info;
  grid_t grid, w_grid;
  float *R=0,*T=0,*S=0, mask;
  double *z=0;
  grid_t topogrid;
  
/**----------------------------------------------------------------------------
  get additional bathymetry */
//   if(bathymetry==0) bathymetry=strdup("/home/softs/genesis/data/topography/gebco/gridone.grd");
//   if(bathymetry!=NULL) {
//     status=grd_loadgrid(bathymetry,&topogrid);
//     if(status !=0) {
//       STDOUT_BASE_LINE("cannot load bathymetry file=%s\n",bathymetry);
//       TRAP_ERR_EXIT(-1,"exiting\n");
//       }
//     exitIfNull(
//       topo=new float[topogrid.nx*topogrid.ny]
//       );
//     topogrid.modeH=0;
//     status=  grd_loadr1(bathymetry,topogrid,topo,&topomask);
//     }
//   else {
//     topo=0;
//     }
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  read temperature (if any)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(input[TEM_ID]!="") {
    filename=input[TEM_ID];
    varname=variables[TEM_ID];
    if(varname=="") TRAP_ERR_EXIT(-1,"please provide a temperature variable name using -tvar XXX option\n");
    status=poc_inq_var(filename,varname,&info);
    NC_TRAP_ERROR(wexit,status,1,"poc_inq_var(\""+filename+"\",\""+varname+"\",) error");
    status=poc_get_grid(filename, info, &grid, verbose, 0);
    NC_TRAP_ERROR(wexit,status,1,"poc_get_grid(\""+filename+"\",(\""+varname+"\"),,,) error");
    T=new float[grid.size()];
    status=poc_get_vara(filename,info,0,T);
    NC_TRAP_ERROR(wexit,status,1,"poc_get_vara(\""+filename+"\",(\""+varname+"\"),,0,) error");
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  read salinity (if any)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(input[SAL_ID]!="") {
    filename=input[SAL_ID];
    varname=variables[SAL_ID];
    if(varname=="") TRAP_ERR_EXIT(-1,"please provide a salinity variable name using -svar XXX option\n");
    status=poc_inq_var(filename,varname,&info);
    NC_TRAP_ERROR(wexit,status,1,"poc_inq_var(\""+filename+"\",\""+varname+"\",) error");
    status=poc_get_grid(filename, info, &grid, verbose, 0);
    NC_TRAP_ERROR(wexit,status,1,"poc_get_grid(\""+filename+"\",(\""+varname+"\"),,,) error");
    T=new float[grid.size()];
    status=poc_get_vara(filename,info,0,T);
    NC_TRAP_ERROR(wexit,status,1,"poc_get_vara(\""+filename+"\",(\""+varname+"\"),,0,) error");
    }

//   if(strcmp(units,"K")==0) {
//     for(k=0;k<grid.nz;k++) {
//       for(j=0;j<grid.ny;j++) {
//         for(int i=0;i<grid.nx;i++) {
//            m=k*grid.ny*grid.nx+j*grid.nx+i;
//            if(T[m]!=mask) {
//              T[m]-=273.15;
//              }
//            }
//          }
//        }
//      }
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  read potential density (if any)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  size_t count=0;
  if(input[RHO_ID]!="") {
    poc_data_t<float> density;
    
    filename=input[RHO_ID];
    varname=variables[RHO_ID];
    if(varname=="") TRAP_ERR_EXIT(-1,"please provide a potential density variable name using -rvar XXX option\n");

    status=poc_get_grid(filename, varname, &grid, verbose, 0);
    status=density.init(filename,varname);
    NC_TRAP_ERROR(wexit,status,1,"cannot load density climatology\npoc_data_t::init(\""+filename+"\",\""+varname+"\") error\n");
    status=density.read_data(filename,0);
    
    mask=density.mask;
    swapval(R,density.data);
    
    if(mask==0.){
      printf("mask value of "+varname+" is 0, which is error-prone: changing to NC_FILL_FLOAT=%g.\n",NC_FILL_FLOAT);
      for(m=0;m<density.length;m++){
        float *Rm=&R[m];
        if(*Rm==mask)
          *Rm=NC_FILL_FLOAT;
        }
      mask=NC_FILL_FLOAT;
      }
    
    k=grid.nz-1;
    count=occurence(mask, &R[k*grid.ny*grid.nx], grid.Hsize());
    
/*------------------------------------------------------------------------------
    NEMO mask is often badly documented, replace zero value with mask to avoid further assles */
    if(count==0) {
//       mask=0.0;
//       count=occurence(mask, &R[k*grid.ny*grid.nx], grid.Hsize());
      for(size_t k=0;k<grid.size();k++) if(R[k]==0.0) R[k]=mask;
      count=occurence(mask, &R[k*grid.ny*grid.nx], grid.Hsize());
      }
/*------------------------------------------------------------------------------
    check range of density, add 1000 if reduced density */
    range_t<float> r=poc_minmax(R, grid.size(),mask);
    if(r.max<100.0) for(m=0;m<density.length;m++) if(R[m]!=mask) R[m]+=1000.;
    if(count==grid.Hsize()) printf(varname+" : last layer entirely masked (NEMO does)\n");
    }
  
  if(grid.modeV==1) {
    z=new double[grid.nz*grid.ny*grid.nx];
    for(k=0;k<grid.nz;k++) {
      for(j=0;j<grid.ny;j++) {
        for(int i=0;i<grid.nx;i++) {
           m=k*grid.ny*grid.nx+j*grid.nx+i;
           z[m]=grid.z[k];
           }
         }
       }
    delete[] grid.z;
    grid.z=z;
    grid.modeV=3;
    }
  
  if(grid.ny==1) grid.y[0]=49.8;
  
/*------------------------------------------------------------------------------
  mask filename to be pushed in options */
/*------------------------------------------------------------------------------
  possible issue: we assume tmask and rho are on same grid */
  poc_data_t<int8_t> tmask,umask;
  status=tmask.init(maskFile,"tmask",1);
  status=tmask.read_data(maskFile,-1,1);
  
  k=grid.nz-1;
  count=occurence((int8_t) 0, &tmask.data[k*grid.ny*grid.nx], grid.Hsize());
  if(count==grid.Hsize()) printf("tmask : last layer entirely masked (NEMO does)\n");
  
//   count=occurence((int8_t) 0, tmask.data, grid.size());
//   count=occurence(mask, R, grid.size());

#if 0
/*------------------------------------------------------------------------------
  check mask and _FillValue consistency */
  if(R!=0 and tmask.data!=0) {
    count=0;
    const int gridHSize=grid.Hsize();
    for(j=0;j<grid.ny;j++) {
      for(int i=0;i<grid.nx;i++) {
        for(k=0;k<grid.nz;k++) {
          m=k*gridHSize+j*grid.nx+i;
          if(tmask.data[m]==0 and R[m]!=mask) {
            if(k>0){
              printf("%s: inconsistent mask at i=%d j=%d k=%d R[m]=%g R[m-HSize]=%g\n",__func__,i,j,k,R[m],R[m-gridHSize]);
              count++;
              }
            break;
            }
          }
        }
      }
    }
#endif

  if(R!=0){
    count=occurence(mask, R, grid.size());
    for(size_t m=0;m<grid.size();m++) {
      if(tmask.data[m]==0) R[m]=mask;
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  e3t is layer thickness
  e3w is distance between 2 tracer nodes
  
  zw[0]=0
  zw[k+1]=zw[k]+e3t[k]
  
  zt[k]=zw[k]+0.5*e3w[k] (i.e. levels are set a tracer nodes mid-position)
    
  zt[k+1]-zt[k]=e3w[k+1]
               =zw[k+1]+0.5*e3w[k+1]-zw[k]-0.5*e3w[k]
               =e3t[k]+0.5*(e3w[k+1]-e3w[k])
  
  so we should verify e3w[k+1]=e3t[k]+0.5*(e3w[k+1]-e3w[k]), i.e. e3t[k]=0.5*(e3w[k+1]+e3w[k])
  
  This condition hold for layers valid at k and k+1
  
  Partial step issue (k+1):
  
  zt[k+1]-zt[k]=e3w[k+1] holds
  
  zt[k+1]=zw[k+1]+0.5*e3w[k+1] DOES NOT hold
  
  zt[k+1]<zw[k+1]+0.5*e3w[k+1]
  
  e3w[k+1]<e3t[k]+0.5*(e3w[k+1]-e3w[k]) i.e. e3t[k]>0.5*(e3w[k+1]+e3w[k])
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
/*------------------------------------------------------------------------------
  vertical grid variables to be pushed in options */

  poc_data_t<double> e3t,e3w;
  
  status=e3t.init(maskFile,"e3t");
//   status=e3t.init(maskFile,"e3t_0");
  status=e3t.read_data(maskFile,-1);
  NC_TRAP_ERROR(wexit,status,1,"e3t reading failed\n");
  
  status=e3w.init(maskFile,"e3w");
//   status=e3w.init(maskFile,"e3w_0");
  status=e3w.read_data(maskFile,-1);
  NC_TRAP_ERROR(wexit,status,1,"e3w reading failed\n");
  
  for(size_t m=0;m<grid.size();m++) {
//     e3t.data[m]*=tmask.data[m];
    if(tmask.data[m]==0) e3t.data[m]=mask;
    }
  for(j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      for(k=0;k<grid.nz;k++) {
        size_t m=k*grid.ny*grid.nx+j*grid.nx+i;
        if(tmask.data[m]==0) e3w.data[m]=mask;
        }
      }
    }
  
//   /* compatibility check beween e3t and e3w */
//   for(j=0;j<grid.ny;j++) {
//     for(int i=0;i<grid.nx;i++) {
//       size_t n=j*grid.nx+i;
//       for(k=0;k<grid.nz-1;k++) {
// //         size_t m0=(k-1)*grid.ny*grid.nx+n;
//         size_t m1=k*grid.ny*grid.nx+n;
//         if(tmask.data[m1]==0) continue;
//         size_t m2=(k+1)*grid.ny*grid.nx+n;
//         if(tmask.data[m2]==0) continue;
//         double v=e3t.data[m1]-0.5*e3w.data[m1]-0.5*e3w.data[m2];
//         if(fabs(v)>1.0) {
//           v=0;
//           }
//         }
//       }
//     }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  recomputing T grid from the integration of e3w
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  aset(grid.z,grid.size(),grid.zmask);
  
  for(j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      size_t n=j*grid.nx+i;
      double depth=0.0;
      k=0;
      m=k*grid.ny*grid.nx+j*grid.nx+i;
      if(e3w.data[m]!=mask) grid.z[m]=-e3w.data[m]/2.0;;
      for(k=1;k<grid.nz;k++) {
        m=k*grid.ny*grid.nx+n;
        if(e3w.data[m]==mask) continue;
        depth+=e3w.data[m];
        grid.z[m]=-depth;
        }
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  recomputing W grid from the integration of e3t
  implemented because W depth is not always available from mesh_mask.nc
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  w_grid=grid;
  w_grid.nz++;
  w_grid.z=aset(w_grid.size(),w_grid.zmask);
  
  for(j=0;j<w_grid.ny;j++) {
    for(int i=0;i<w_grid.nx;i++) {
      size_t n=j*grid.nx+i;
      double depth=0.0;
      k=0;
      m=k*w_grid.ny*w_grid.nx+j*w_grid.nx+i;
      if(e3t.data[m]!=mask) w_grid.z[m]=0;
      for(k=1;k<w_grid.nz;k++) {
        size_t m0=(k-1)*grid.ny*grid.nx+n;
        if(e3t.data[m0]==mask) continue;
        m=k*w_grid.ny*w_grid.nx+j*w_grid.nx+i;
        depth+=e3t.data[m0];
        w_grid.z[m]=-depth;
        }
      }
    }


//   R=new float[grid.nz*grid.ny*grid.nx];
//   for(k=0;k<grid.nz;k++) {
//     for(j=0;j<grid.ny;j++) {
//       for(int i=0;i<grid.nx;i++) {
//          m=k*grid.ny*grid.nx+j*grid.nx+i;
//          x=map_grid_x (grid, i, j);
//          y=map_grid_y (grid, i, j);
//          x=map_recale(topogrid,x);
//          status=map_interpolation(topogrid, topo,topomask,x,y,&h);
//          if(T[m]!=mask) {
//            R[m]=water_density(T[m],S[m],grid.z[m]);
//            }
//          else {
//            R[m]=mask;
//            }
//          if(grid.z[m]>-h) {
//            T[m]=mask;
//            S[m]=mask;
//            R[m]=mask;
//            }
//          }
//        }
//      }

//   delete[] topo;
  topogrid.free();
  
  status=compute_vertical_mode(bathymetry, nRequestedProcs, grid, w_grid, T, S, R, mask, rootname, maxmodes);

  deletep(&R);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int NEMO_FV(const string *input,const string maskFile, const string *variables,const char *rootname,const char *bathymetry, int nRequestedProcs, const char *units, int maxmodes, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int j,k,m;
  string filename,varname;
  poc_var_t info;
  grid_t grid, w_grid;
  float *R=0,*T=0,*S=0, mask;
  double *z=0;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  read potential density (if any)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  size_t count=0;
  if(input[RHO_ID]!="") {
    poc_data_t<float> density;
    
    filename=input[RHO_ID];
    varname=variables[RHO_ID];
    if(varname=="") TRAP_ERR_EXIT(-1,"please provide a potential density variable name using -rvar XXX option\n");

    status=poc_get_grid(filename, varname, &grid, verbose, 0);
    status=density.init(filename,varname);
    NC_TRAP_ERROR(wexit,status,1,"cannot load density climatology\npoc_data_t::init(\""+filename+"\",\""+varname+"\") error\n");
    status=density.read_data(filename,0);
    NC_TRAP_ERROR(wexit,status,1,"cannot load density climatology\npoc_data_t(,\""+varname+"\")::read_data(\""+filename+"\",) error\n");
    
    mask=density.mask;
    swapval(R,density.data);
    
    if(mask==0.){
      printf("mask value of "+varname+" is 0, which is error-prone: changing to NC_FILL_FLOAT=%g.\n",NC_FILL_FLOAT);
      for(m=0;m<density.length;m++){
        float *Rm=&R[m];
        if(*Rm==mask)
          *Rm=NC_FILL_FLOAT;
        }
      mask=NC_FILL_FLOAT;
      }
    
    k=grid.nz-1;
    count=occurence(mask, &R[k*grid.ny*grid.nx], grid.Hsize());
/*------------------------------------------------------------------------------
    check range of density, add 1000 if reduced density */
    range_t<float> r=poc_minmax(R, grid.size(),mask);
    if(r.max<100.0) for(m=0;m<density.length;m++) if(R[m]!=mask) R[m]+=1000.;
    if(count==grid.Hsize()) printf(varname+" : last layer entirely masked (NEMO does)\n");
    }
  
  if(grid.modeV==1) {
    z=new double[grid.nz*grid.ny*grid.nx];
    for(k=0;k<grid.nz;k++) {
      for(j=0;j<grid.ny;j++) {
        for(int i=0;i<grid.nx;i++) {
           m=k*grid.ny*grid.nx+j*grid.nx+i;
           z[m]=grid.z[k];
           }
         }
       }
    delete[] grid.z;
    grid.z=z;
    grid.modeV=3;
    }
  
  if(grid.ny==1) grid.y[0]=49.8;
  
/*------------------------------------------------------------------------------
  possible issue: we assume tmask and rho are on same grid */
  poc_data_t<int8_t> umask;
  status=umask.init(maskFile,"umask",1);
  status=umask.read_data(maskFile,-1,1);
  
  k=grid.nz-1;
  count=occurence((int8_t) 0, &umask.data[k*grid.ny*grid.nx], grid.Hsize());
  if(count==grid.Hsize()) printf("umask : last layer entirely masked (NEMO does)\n");
  
//   count=occurence((int8_t) 0, tmask.data, grid.size());
//   count=occurence(mask, R, grid.size());


  if(R!=0){
    count=occurence(mask, R, grid.size());
    for(size_t m=0;m<grid.size();m++) {
      if(umask.data[m]==0) R[m]=mask;
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  e3t is layer thickness
  e3w is distance between 2 tracer nodes
  
  zw[0]=0
  zw[k+1]=zw[k]+e3t[k]
  
  zt[k]=zw[k]+0.5*e3w[k] (i.e. levels are set a tracer nodes mid-position)
    
  zt[k+1]-zt[k]=e3w[k+1]
               =zw[k+1]+0.5*e3w[k+1]-zw[k]-0.5*e3w[k]
               =e3t[k]+0.5*(e3w[k+1]-e3w[k])
  
  so we should verify :
  
  e3w[k+1]=e3t[k]+0.5*(e3w[k+1]-e3w[k]), i.e. e3t[k]=0.5*(e3w[k+1]+e3w[k])
  
  This condition hold for layers valid at k and k+1
  
  
  Partial step issue (k+1):
  =========================
  
  * zt[k+1]-zt[k]=e3w[k+1] holds
  
  * zt[k+1]=zw[k+1]+0.5*e3w[k+1] DOES NOT hold
  
      zt[k+1]<zw[k+1]+0.5*e3w[k+1]
  
      e3w[k+1]<e3t[k]+0.5*(e3w[k+1]-e3w[k]) i.e. e3t[k]>0.5*(e3w[k+1]+e3w[k])
      
  
  Reverse construction (layer-immersion from level-immersion):
  
  zt[0]=0.5*zw[1]
  
  zt[k]=zt[k-1]-2.0*(zt[k-1]-zw[k])=zw[k]-(zt[k-1]-zw[k])
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  
  poc_data_t<double> e3u;
  status=e3u.init(maskFile,"e3u");
  NC_TRAP_ERROR(wexit,status,1,"e3u init failed");
  status=e3u.read_data(maskFile,-1);
  NC_TRAP_ERROR(wexit,status,1,"e3u reading failed");
  
  for(size_t m=0;m<grid.size();m++) {
    if(umask.data[m]==0) e3u.data[m]=mask;
    }
  
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  recomputing leve-grid immersion from the integration of e3u
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  w_grid=grid;
  w_grid.nz++;
  w_grid.z=aset(w_grid.size(),w_grid.zmask);
  
  for(j=0;j<w_grid.ny;j++) {
    for(int i=0;i<w_grid.nx;i++) {
      size_t m,n=j*w_grid.nx+i;
      double depth=0.0;
      k=0;
      m=k*w_grid.Hsize()+n;
//       if(e3u.data[m]!=mask) w_grid.z[m]=0;
      if(umask.data[m]!=0) w_grid.z[m]=0;
      else continue;
      for(k=1;k<w_grid.nz;k++) {
        depth+=e3u.data[m];
        m=k*w_grid.Hsize()+n;
        w_grid.z[m]=-depth;
        }
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  recomputing layer-grid immersion from the integration of W grid
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  aset(grid.z,grid.size(),grid.zmask);
  
  for(j=0;j<grid.ny;j++) {
    for(int i=0;i<grid.nx;i++) {
      size_t m,m1,m2,n=j*grid.nx+i;
      double e3w;
      k=0;
      m=k*grid.Hsize()+n;
      m2=(k+1)*grid.Hsize()+n;
      if(e3u.data[m]!=mask) grid.z[m]=w_grid.z[m2]/2;
      else continue;
      for(k=1;k<grid.nz;k++) {
        m=k*grid.ny*grid.nx+n;
        m1=(k-1)*grid.Hsize()+n;
        if(e3u.data[m]==mask) continue;
        e3w=2.0*(grid.z[m1]-w_grid.z[m]);
        grid.z[m]=w_grid.z[m]-e3w/2.0;
        }
      }
    }
  
  status=compute_vertical_mode(bathymetry, nRequestedProcs, grid, w_grid, T, S, R, mask, rootname, maxmodes);

  deletep(&R);
  
  return(0);
}

