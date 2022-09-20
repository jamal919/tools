
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

extern  double * get_layer_thickness(const string & rootname,const grid_t & grid,int verbose);

extern  int VerticalModes_SnapshotDecomposition(const string *input, const string & filename, const string *variables, int maxmodes, bool debug);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double integrale_time_2(complex<double> z1, complex<double> z2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double d;
/*******************************************************************************
  Time average over 1 period
*******************************************************************************/

  d = 0.5 * (real(z1) * real(z2) + imag(z1) * imag(z2));
  return (d);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int SpectralEnergy_RAW(const string maskFile, grid_t & grid, complex<double> **u, complex<double> **v, complex<double> **w, complex<double> **p, string & wave, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nlayers, status;
  complex<double> umask;
  string filename, output_energy;
  int  create_file=1, create_grid=1;
  int frame=0;
  poc_data_t<double> scalarData;
  double *uflux,*vflux,dmask=NC_FILL_DOUBLE;
    
  output_energy=wave+"-raw-energy.nc";
  
  umask=NC_FILL_COMPLEX;
  
//   const string maskFile="mesh_mask.nc";
  
  poc_data_t<double> e3t, e3w;
  status=e3t.init(maskFile,"e3t");
  status=e3t.read_data(maskFile,0);
  
//   if(status!=0){
//     NC_CHKERR_BASE_LINE(status,"poc_data_t::read_data(\"+maskFile+\",\"e3t\",0,) error");
//     filename=ModesFile;
//     const string rootname=filename.substr(0,filename.find('.'));
//     STDERR_BASE_LINE("computing it from depths_w in "+filename+"\n");
//     e3t.data=get_layer_thickness(rootname,grid,verbose);
//     }
  
  uflux=new double[grid.Hsize()];
  vflux=new double[grid.Hsize()];
  
  nlayers=grid.nz;
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  modal depth-integrated energy fluxes

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    for(size_t n=0;n<grid.Hsize();n++) {
      uflux[n]=0;
      vflux[n]=0;
      double weight=0;
      for(int k=0;k<nlayers;k++) {
        size_t m=k*grid.Hsize()+n;
        if(p[n][k]!=umask and u[n][k]!=umask and v[n][k]!=umask) {
          uflux[n]+=e3t.data[m]*integrale_time_2(p[n][k], u[n][k]);
          vflux[n]+=e3t.data[m]*integrale_time_2(p[n][k], v[n][k]);
          weight+=e3t.data[m];
          }
        }
      if(weight==0) {
        uflux[n]=dmask;
        vflux[n]=dmask;
        }
      }
      
    create_file=1;
    create_grid=1;
    
    status=save_SGXY(output_energy.c_str(), grid, uflux, dmask, "rawflux_x","N/m²", "depth_integrated_energy_flux", create_file, create_grid);
    
    create_file=0;
    create_grid=0;
    
    status=save_SGXY(output_energy.c_str(), grid, vflux, dmask, "rawflux_y","N/m²", "depth_integrated_energy_flux", create_file, create_grid);
    
    double *divergence=new double[grid.Hsize()];
    status=map_divergence(grid, uflux, vflux, dmask, GEOCENTRIC, divergence);
    status=save_SGXY(output_energy.c_str(), grid, divergence, dmask, "divergence","w/m²", "depth_integrated_energy_flux_divergence", create_file, create_grid);
    delete[] divergence;
    

  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int VerticalModes_SpectralEnergy(const string maskFile, string ModesFile, int maxmodes, string & wave, int nrecomposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,n,nmasked,status;
  grid_t grid;
  complex<double> **p_decomposition=0, **u_decomposition=0, **v_decomposition=0, umask;
  double **Pmodes=0, **Umodes=0, mask;
  complex<double>  *pbaroclinic=0, *ubaroclinic=0, *vbaroclinic=0;
  int nlayers, nkept, MinUmodes;
  int create_file, create_grid;
  double *uflux,*vflux,dmask=NC_FILL_DOUBLE;
  double *uflux3D,*vflux3D;
  string filename, varname, varnameA, varnameG, output_energy;
  poc_data_t<double> scalarData;
  
  nkept=min(maxmodes, nrecomposed);
  
  output_energy=wave+"-modal-energy.nc";
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  reload vertical modes

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  filename=ModesFile;
  varname="Pmodes";
  
  status=poc_get_grid(filename, varname, &grid, verbose, -1);
  
  nlayers=grid.nz;
  
  MinUmodes=nlayers;
  
  status=scalarData.init(filename,varname);
  for(i=0;i<MinUmodes;i++) {
    status=scalarData.read_data(filename,i);
    MinUmodes=min(scalarData.nframes, maxmodes);
    if(Pmodes==0) {
      Pmodes=new double*[MinUmodes];
      for(k=0;k<MinUmodes;k++) Pmodes[k]=0;
      }
    if(Pmodes[i]==0) Pmodes[i]=new double[grid.size()];
    for(n=0;n<grid.size();n++) Pmodes[i][n]=scalarData.data[n];
    mask=scalarData.mask;
    nmasked=occurence<int>(mask, Pmodes[i], grid.size());
    }
  scalarData.destroy_data();
  
  varname="Umodes";
  
  status=scalarData.init(filename,varname);
  for(i=0;i<MinUmodes;i++) {
    status=scalarData.read_data(filename,i);
    MinUmodes=min(scalarData.nframes, maxmodes);
    if(Umodes==0) {
      Umodes=new double*[MinUmodes];
      for(k=0;k<MinUmodes;k++) Umodes[k]=0;
      }
    if(Umodes[i]==0) Umodes[i]=new double[grid.size()];
    for(n=0;n<grid.size();n++) Umodes[i][n]=scalarData.data[n];
    mask=scalarData.mask;
    nmasked=occurence<int>(mask, Umodes[i], grid.size());
    }
  scalarData.destroy_data();
  

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  reload modal decomposition

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  p_decomposition=new complex<double>*[MinUmodes];
  for(i=0;i<MinUmodes;i++) {
    p_decomposition[i]=new complex<double>[grid.Hsize()];
    for(n=0;n<grid.Hsize();n++) p_decomposition[i][n]=umask;
    }
  u_decomposition=new complex<double>*[MinUmodes];
  for(i=0;i<MinUmodes;i++) {
    u_decomposition[i]=new complex<double>[grid.Hsize()];
    for(n=0;n<grid.Hsize();n++) u_decomposition[i][n]=umask;
    }
  v_decomposition=new complex<double>*[MinUmodes];
  for(i=0;i<MinUmodes;i++) {
    v_decomposition[i]=new complex<double>[grid.Hsize()];
    for(n=0;n<grid.Hsize();n++) v_decomposition[i][n]=umask;
    }
  
  filename=wave+"-modal-decomposition.nc";
  
  complex<double> *buffer=new complex<double>[grid.size()];
  
  varnameA="Uc_a";
  varnameG="Uc_G";
  for(i=0;i<MinUmodes;i++) {
    status=poc_get_cvara(filename,varnameA, varnameG, i, buffer);
    for(n=0;n<grid.Hsize();n++) u_decomposition[i][n]=buffer[n];
    }
  
  varnameA="Vc_a";
  varnameG="Vc_G";
  for(i=0;i<MinUmodes;i++) {
    status=poc_get_cvara(filename,varnameA, varnameG, i, buffer);
    for(n=0;n<grid.Hsize();n++) v_decomposition[i][n]=buffer[n];
    }
  
  varnameA="Pc_a";
  varnameG="Pc_G";
  for(i=0;i<MinUmodes;i++) {
    status=poc_get_cvara(filename,varnameA, varnameG, i, buffer);
    for(n=0;n<grid.Hsize();n++) p_decomposition[i][n]=buffer[n];
    }
  
  umask=NC_FILL_COMPLEX;
  dmask=mask;
  
  poc_data_t<double> e3t, e3w;
  status=e3t.init(maskFile,"e3t");
  if(status!=0){
    NC_CHKERR_BASE_LINE(status,"poc_data_t::read_data(\"mesh_mask.nc\",\"e3t\",0,) error");
    filename=ModesFile;
    const string rootname=filename.substr(0,filename.find('.'));
    STDERR_BASE_LINE("computing it from depths_w in "+filename+"\n");
    e3t.data=get_layer_thickness(rootname,grid,verbose);
    }
  else
    status=e3t.read_data(maskFile,0);
  
//   poc_data_t<float> hdepw;
//   filename=maskFile;
//   status=hdepw.init(filename,"hdepw",0);
//   if(status!=0) do{
//     NC_CHKERR_BASE_LINE(status,"poc_data_t::init(\""+filename+"\",\"hdepw\",0,) error");
//     
//     STDERR_BASE_LINE("getting it from gdepw in "+filename+"\n");
//     status=hdepw.init(filename,"gdepw",1);
//     if(status==0) break;
//     
//     filename="grid.nc";
//     STDERR_BASE_LINE("getting it from h_w in "+filename+"\n");
//     status=hdepw.init(filename,"h_w",1);
//     if(status==0) break;
//     
//     STDERR_BASE_LINE("no, getting it from depp in "+filename+"\n");
//     status=hdepw.init(filename,"depp",1);
//     
//     NC_TRAP_ERROR(wexit,status,1,"poc_data_t::init(\""+filename+"\",(\"hdepw\" or \"gdepw\" or \"h_w\" or \"depp\"),0,) error");
//     }while(0);
//   status=hdepw.read_data(maskFile,0,0);
//   NC_TRAP_ERROR(wexit,status,1,"poc_data_t::read_data(\""+filename+"\",\""+hdepw.info.name+"\",0,) error");
  
  uflux=new double[grid.Hsize()];
  vflux=new double[grid.Hsize()];
  
  uflux3D=new double[grid.Hsize()*nlayers];
  vflux3D=new double[grid.Hsize()*nlayers];
  
  for(size_t m=0;m<grid.Hsize();m++) uflux[m]=dmask;
  for(size_t m=0;m<grid.Hsize();m++) vflux[m]=dmask;
  
  pbaroclinic=new complex<double>[grid.size()];
  for(size_t m=0;m<grid.size();m++) pbaroclinic[m]=umask;
  ubaroclinic=new complex<double>[grid.size()];
  for(size_t m=0;m<grid.size();m++) ubaroclinic[m]=umask;
  vbaroclinic=new complex<double>[grid.size()];
  for(size_t m=0;m<grid.size();m++) vbaroclinic[m]=umask;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  modal depth-integrated energy fluxes

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(i=0;i<nkept;i++) {
    for(n=0;n<grid.Hsize();n++) {
      for(k=0;k<nlayers;k++) {
        size_t m=k*grid.Hsize()+n;
        pbaroclinic[m]=umask;
        ubaroclinic[m]=umask;
        vbaroclinic[m]=umask;
        if(Pmodes[i][m]==mask or p_decomposition[i][n]==umask) continue;
        if(Umodes[i][m]==mask or u_decomposition[i][n]==umask) continue;
        if(Umodes[i][m]==mask or v_decomposition[i][n]==umask) continue;
        pbaroclinic[m]=Pmodes[i][m]*p_decomposition[i][n];
        ubaroclinic[m]=Umodes[i][m]*u_decomposition[i][n];
        vbaroclinic[m]=Umodes[i][m]*v_decomposition[i][n];
        }
      uflux[n]=0;
      vflux[n]=0;
      double weight=0;
      for(k=0;k<nlayers;k++) {
        size_t m=k*grid.Hsize()+n;
        if(pbaroclinic[m]!=umask) {
          uflux[n]+=e3t.data[m]*integrale_time_2(pbaroclinic[m], ubaroclinic[m]);
          vflux[n]+=e3t.data[m]*integrale_time_2(pbaroclinic[m], vbaroclinic[m]);
          weight+=e3t.data[m];
          }
        }
      if(weight==0) {
        uflux[n]=dmask;
        vflux[n]=dmask;
        }
      }
    create_file=(i==0);
    create_grid=(i==0);
    int frame=i;
    status=save_SGXYT(output_energy.c_str(), grid, uflux, dmask, "Fmodal_x","N/m²", "Fmodal_x", create_file, create_grid, frame);
    create_file=0;
    create_grid=0;
    status=save_SGXYT(output_energy.c_str(), grid, vflux, dmask, "Fmodal_y","N/m²", "Fmodal_y", create_file, create_grid, frame);
    if(i==0) {
      status=save_SGXY(output_energy.c_str(), grid, uflux, dmask, "Fbarotropic_x","N/m²", "Fbarotropic_x", create_file, create_grid);
      status=save_SGXY(output_energy.c_str(), grid, vflux, dmask, "Fbarotropic_y","N/m²", "Fbarotropic_y", create_file, create_grid);
      double *divergence=new double[grid.Hsize()];
      status=map_divergence(grid, uflux, vflux, dmask, GEOCENTRIC, divergence);
      status=save_SGXY(output_energy.c_str(), grid, divergence, dmask, "Dbarotropic","w/m²", "Dbarotropic", create_file, create_grid);
      delete[] divergence;
      }
    }
    

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  baroclinic depth-integrated energy fluxes

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(n=0;n<grid.Hsize();n++) {
    uflux[n]=0;
    vflux[n]=0;
    for(k=0;k<nlayers;k++) {
      size_t m=k*grid.Hsize()+n;
      uflux3D[m]=0.0;
      vflux3D[m]=0.0;
      }
    for(i=1;i<nkept;i++) {
      for(k=0;k<nlayers;k++) {
        size_t m=k*grid.Hsize()+n;
        pbaroclinic[m]=umask;
        ubaroclinic[m]=umask;
        vbaroclinic[m]=umask;
        if(Pmodes[i][m]==mask or p_decomposition[i][n]==umask) continue;
        if(Umodes[i][m]==mask or u_decomposition[i][n]==umask) continue;
        if(Umodes[i][m]==mask or v_decomposition[i][n]==umask) continue;
        pbaroclinic[m]=Pmodes[i][m]*p_decomposition[i][n];
        ubaroclinic[m]=Umodes[i][m]*u_decomposition[i][n];
        vbaroclinic[m]=Umodes[i][m]*v_decomposition[i][n];
        }
      double weight=0;
      for(k=0;k<nlayers;k++) {
        size_t m=k*grid.Hsize()+n;
        if(pbaroclinic[m]!=umask) {
          uflux[n]+=e3t.data[m]*integrale_time_2(pbaroclinic[m], ubaroclinic[m]);
          vflux[n]+=e3t.data[m]*integrale_time_2(pbaroclinic[m], vbaroclinic[m]);
          uflux3D[m]+=integrale_time_2(pbaroclinic[m], ubaroclinic[m]);
          vflux3D[m]+=integrale_time_2(pbaroclinic[m], vbaroclinic[m]);
          weight+=e3t.data[m];
          }
        }
      if(weight==0) {
        uflux[n]=dmask;
        vflux[n]=dmask;
        }
      }
    }
  status=save_SGXY(output_energy.c_str(), grid, uflux, dmask, "Fbaroclinic_x","N/m²", "Fbaroclinic_x", create_file, create_grid);
  status=save_SGXY(output_energy.c_str(), grid, vflux, dmask, "Fbaroclinic_y","N/m²", "Fbaroclinic_y", create_file, create_grid);
  
  create_grid=1;
  
  status=save_SGXYZ(output_energy.c_str(), grid, uflux3D, dmask, "Fbaroclinic3D_x","N/m²", "Fbaroclinic3D_x", create_file, create_grid,"","");
  
  create_grid=0;
  
  status=save_SGXYZ(output_energy.c_str(), grid, vflux3D, dmask, "Fbaroclinic3D_y","N/m²", "Fbaroclinic3D_y", create_file, create_grid,"","");
  
  double *divergence=new double[grid.Hsize()];
  status=map_divergence(grid, uflux, vflux, dmask, GEOCENTRIC, divergence);
  status=save_SGXY(output_energy.c_str(), grid, divergence, dmask, "Dbaroclinic","w/m²", "Dbaroclinic", create_file, create_grid);
  delete[] divergence;
  
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  internal tide production

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  float *hx=new float[grid.Hsize()];
  float *hy=new float[grid.Hsize()];
  float *h=new float[grid.Hsize()];
  
#if 0
  for(n=0;n<grid.Hsize();n++) {
    size_t m=n;
    k=0;
    while(Pmodes[0][m]!=mask and k<nlayers-1) {
      k++;
      m=k*grid.Hsize()+n;
      }
    if(Pmodes[0][m]==mask) k--;
    m=k*grid.Hsize()+n;
    h[n]=hdepw.data[m];
    }
#endif
  
  status=e3w.init(maskFile,"e3w");
  status=e3w.read_data(maskFile,0);
  for(n=0;n<grid.Hsize();n++) {
    h[n]=0.0;
    size_t m=n;
    k=0;
    while(Pmodes[0][m]!=mask and k<nlayers) {
      h[n]+=e3w.data[m];
      k++;
      m=k*grid.Hsize()+n;
      }
    }
  
  status=map_gradient(grid, grid.nx, h, mask,  GEOCENTRIC, hx, hy);
  
  double *production=new double[grid.Hsize()];

  for(n=0;n<grid.Hsize();n++) {
    production[n]=dmask;
    complex<double>  pbottom=0;
    k=0;
    size_t m=n;
    i=0;
    while(Pmodes[i][m]!=mask and k<nlayers-1) {
      k++;
      m=k*grid.Hsize()+n;
      }
    if(Pmodes[i][m]==mask) k--;
    if(k<0) continue;
    m=k*grid.Hsize()+n;
    for(i=1;i<nkept;i++) {
      if(Pmodes[i][m]==mask or p_decomposition[i][n]==umask) continue;
      pbottom+=Pmodes[i][m]*p_decomposition[i][n];
      }
    i=0;
    if(Umodes[i][m]==mask or u_decomposition[i][n]==umask) continue;
    if(Umodes[i][m]==mask or v_decomposition[i][n]==umask) continue;
    complex<double> uu=Umodes[i][m]*u_decomposition[i][n];
    complex<double> vv=Umodes[i][m]*v_decomposition[i][n];
    production[n] = hx[n]*integrale_time_2(pbottom, uu);
    production[n]+= hy[n]*integrale_time_2(pbottom, vv);
//     if(hx[n]!=0) {
//       printf("%d %lf %lf (%lf %lf) \n",n,h[n],hx[n],h[n-1],h[n+1]);
//       }
    if(isnan(production[n])==1) {
      production[n]=dmask;
      }
    }
  status=save_SGXY(output_energy.c_str(), grid, production, dmask, "production","w/m²", "production", create_file, create_grid);
  delete[] production;
  
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  additional diagnostics (wavelength)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  complex<double> *laplacian=0;
  int ref=0;
  float *wavelength=new float[grid.Hsize()], rmask=1.e+10;
    
//   for(int i=0;i<MinUmodes;i++) {
  status=map_laplacian(grid, u_decomposition[1], umask, ref, laplacian);
  for(size_t m=0;m<grid.Hsize();m++) {
    wavelength[m]=rmask;
    if(laplacian[m]!=umask and u_decomposition[1][m]!=0.0) {
      complex<double> ratio=laplacian[m]/u_decomposition[1][m];
      double k2=-real(ratio);
      if(k2>0.0) wavelength[m]=2.*M_PI/sqrt(k2);
      else {
        if(verbose==1) printf("%s anomaly\n",__func__);
        }
      }
    }
  create_file=1;
  status=save_SGXY("modal-diagnostics.nc", grid, wavelength, rmask, "uL", "m", "uL", create_file);
//     }

  status=map_laplacian(grid, p_decomposition[1], umask, ref, laplacian);
  for(size_t m=0;m<grid.Hsize();m++) {
    wavelength[m]=rmask;
    if(laplacian[m]!=umask and p_decomposition[1][m]!=0.0) {
      complex<double> ratio=laplacian[m]/p_decomposition[1][m];
      double k2=-real(ratio);
      if(k2>0.0) wavelength[m]=2.*M_PI/sqrt(k2);
      else {
        if(verbose==1) printf("%s anomaly\n",__func__);
        }
      }
    }
  create_file=0;
  status=save_SGXY("modal-diagnostics.nc", grid, wavelength, rmask, "pL", "m", "pL", create_file);
  
    
  deletep2D(&p_decomposition,MinUmodes);
  deletep2D(&u_decomposition,MinUmodes);
  deletep2D(&v_decomposition,MinUmodes);

/*------------------------------------------------------------------------------
  only for verification, commented to keep files small */ 
//   for(n=0;n<grid.Hsize();n++) {
//     for(k=0;k<nlayers;k++) {
//       size_t m=k*grid.Hsize()+n;
//       if(u[n][k]!=umask) pbarotropic[m]=u[n][k]/9.81/1019.2;
//       else pbarotropic[m]=umask;
//       }
//     }
//   create_file=0;
//   =1;
//   status=save_SGXYZ_C("modal-decomposition.nc", grid, pbarotropic, umask, "P_a", "P_G", "m", "u_coeffcients", create_file, create_grid);
 
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int VerticalModes_SpectralDiagnostics(const string maskFile, string ModesFile, int maxmodes, string & wave, int nrecomposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,n,nmasked,status;
  grid_t grid;
  complex<double> **p_decomposition=0, **u_decomposition=0, **v_decomposition=0, umask;
  double **Pmodes=0, **Umodes=0, mask;
  complex<double>  *pbaroclinic=0, *ubaroclinic=0, *vbaroclinic=0;
  int nlayers, nkept, MinUmodes;
  int create_file, create_grid;
  double *uflux,*vflux,dmask=NC_FILL_DOUBLE;
  double *uflux3D,*vflux3D;
  string filename, varname, varnameA, varnameG;
  poc_data_t<double> scalarData;
  
  nkept=min(maxmodes, nrecomposed);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  reload vertical modes

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  filename=ModesFile;
  varname="Pmodes";
  
  status=poc_get_grid(filename, varname, &grid, verbose, -1);
  
  nlayers=grid.nz;
  
  MinUmodes=nlayers;
  
  status=scalarData.init(filename,varname);
  for(i=0;i<MinUmodes;i++) {
    status=scalarData.read_data(filename,i);
    MinUmodes=min(scalarData.nframes, maxmodes);
    if(Pmodes==0) {
      Pmodes=new double*[MinUmodes];
      for(k=0;k<MinUmodes;k++) Pmodes[k]=0;
      }
    if(Pmodes[i]==0) Pmodes[i]=new double[grid.size()];
    for(n=0;n<grid.size();n++) Pmodes[i][n]=scalarData.data[n];
    mask=scalarData.mask;
    nmasked=occurence<int>(mask, Pmodes[i], grid.size());
    }
  scalarData.destroy_data();
  
/*------------------------------------------------------------------------------
  reload modal decomposition */
  p_decomposition=new complex<double>*[MinUmodes];
  for(i=0;i<MinUmodes;i++) {
    p_decomposition[i]=new complex<double>[grid.Hsize()];
    for(n=0;n<grid.Hsize();n++) p_decomposition[i][n]=umask;
    }
  u_decomposition=new complex<double>*[MinUmodes];
  for(i=0;i<MinUmodes;i++) {
    u_decomposition[i]=new complex<double>[grid.Hsize()];
    for(n=0;n<grid.Hsize();n++) u_decomposition[i][n]=umask;
    }
  v_decomposition=new complex<double>*[MinUmodes];
  for(i=0;i<MinUmodes;i++) {
    v_decomposition[i]=new complex<double>[grid.Hsize()];
    for(n=0;n<grid.Hsize();n++) v_decomposition[i][n]=umask;
    }
  
//   filename="modal-decomposition.nc";
  filename=wave+"-modal-decomposition.nc";
  
  complex<double> *buffer=new complex<double>[grid.size()];
  
  varnameA="Uc_a";
  varnameG="Uc_G";
  for(i=0;i<MinUmodes;i++) {
    status=poc_get_cvara(filename,varnameA, varnameG, i, buffer);
    for(n=0;n<grid.Hsize();n++) u_decomposition[i][n]=buffer[n];
    }
  
  varnameA="Vc_a";
  varnameG="Vc_G";
  for(i=0;i<MinUmodes;i++) {
    status=poc_get_cvara(filename,varnameA, varnameG, i, buffer);
    for(n=0;n<grid.Hsize();n++) v_decomposition[i][n]=buffer[n];
    }
  
  varnameA="Pc_a";
  varnameG="Pc_G";
  for(i=0;i<MinUmodes;i++) {
    status=poc_get_cvara(filename,varnameA, varnameG, i, buffer);
    for(n=0;n<grid.Hsize();n++) p_decomposition[i][n]=buffer[n];
    }
  
  umask=NC_FILL_COMPLEX;
  dmask=mask;
  
  poc_data_t<double> e3t, e3w;
  status=e3t.init(maskFile,"e3t",0);
  if(status!=0){
    NC_CHKERR_BASE_LINE(status,"poc_data_t::init(\"mesh_mask.nc\",\"e3t\") error");
    filename=ModesFile;
    const string rootname=filename.substr(0,filename.find('.'));
    STDERR_BASE_LINE("computing it from depths_w in "+filename+"\n");
    e3t.data=get_layer_thickness(rootname,grid,verbose);
    }
  else
    status=e3t.read_data(maskFile,0,0);
  
  poc_data_t<float> hdepw;
  filename=maskFile;
  status=hdepw.init(filename,"hdepw",0);
  if(status!=0) do{
    NC_CHKERR_BASE_LINE(status,"poc_data_t::init(\""+filename+"\",\"hdepw\",0,) error");
    
    STDERR_BASE_LINE("getting it from gdepw in "+filename+"\n");
    status=hdepw.init(filename,"gdepw",1);
    if(status==0) break;
    
    filename="grid.nc";
    STDERR_BASE_LINE("getting it from h_w in "+filename+"\n");
    status=hdepw.init(filename,"h_w",1);
    if(status==0) break;
    
    STDERR_BASE_LINE("no, getting it from depp in "+filename+"\n");
    status=hdepw.init(filename,"depp",1);
    
    NC_TRAP_ERROR(wexit,status,1,"poc_data_t::init(\""+filename+"\",(\"hdepw\" or \"gdepw\" or \"h_w\" or \"depp\"),0,) error");
    }while(0);
  status=hdepw.read_data(maskFile,0,0);
  NC_TRAP_ERROR(wexit,status,1,"poc_data_t::read_data(\""+filename+"\",\""+hdepw.info.name+"\",0,) error");
  
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  additional diagnostics (wavelength)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  complex<double> *laplacian=0, *laplacian_u=0, *laplacian_v=0;
  int ref=0;
  float *wavelength=new float[grid.Hsize()], rmask=1.e+10;
  
//   for(size_t m=0;m<grid.Hsize();m++) if(u_decomposition[1][m]!=umask) u_decomposition[1][m]/=abs(u_decomposition[1][m]);
    
//   for(int i=0;i<MinUmodes;i++) {
  status=map_laplacian(grid, u_decomposition[1], umask, ref, laplacian_u);
  for(size_t m=0;m<grid.Hsize();m++) {
    wavelength[m]=rmask;
    if(laplacian_u[m]!=umask and u_decomposition[1][m]!=0.0) {
      complex<double> ratio;
      double k2;
      ratio=laplacian_u[m]/u_decomposition[1][m];
      k2=-real(ratio);
      if(k2>0.0) wavelength[m]=2.*M_PI/sqrt(k2);
      else wavelength[m]=-1;
//       else {
//         if(verbose==1) printf("%s anomaly\n",__func__);
//         }
      }
    }
  create_file=1;
  status=save_SGXY("modal-diagnostics.nc", grid, wavelength, rmask, "uL", "m", "uL", create_file);
  
  status=map_laplacian(grid, v_decomposition[1], umask, ref, laplacian_v);
  for(size_t m=0;m<grid.Hsize();m++) {
    wavelength[m]=rmask;
    if(laplacian_v[m]!=umask and v_decomposition[1][m]!=0.0) {
      complex<double> ratio;
      double k2;
      ratio=laplacian_v[m]/v_decomposition[1][m];
      k2=-real(ratio);
      if(k2>0.0) wavelength[m]=2.*M_PI/sqrt(k2);
      else wavelength[m]=-1;
//       else {
//         if(verbose==1) printf("%s anomaly\n",__func__);
//         }
      }
    }
  create_file=0;
  status=save_SGXY("modal-diagnostics.nc", grid, wavelength, rmask, "vL", "m", "vL", create_file);
//     }

  deletep2D(&u_decomposition,MinUmodes);
  deletep2D(&v_decomposition,MinUmodes);
  deletep(&laplacian_u);
  deletep(&laplacian_v);

  status=map_laplacian(grid, p_decomposition[1], umask, ref, laplacian);
  for(size_t m=0;m<grid.Hsize();m++) {
    wavelength[m]=rmask;
    if(laplacian[m]!=umask and p_decomposition[1][m]!=0.0) {
      complex<double> ratio=laplacian[m]/p_decomposition[1][m];
      double k2=-real(ratio);
      if(k2>0.0) wavelength[m]=2.*M_PI/sqrt(k2);
      else {
        if(verbose==1) printf("%s anomaly\n",__func__);
        }
      }
    }
  create_file=0;
  status=save_SGXY("modal-diagnostics.nc", grid, wavelength, rmask, "pL", "m", "pL", create_file);
  
    
  deletep2D(&p_decomposition,MinUmodes);
  deletep(&laplacian);
  
/*------------------------------------------------------------------------------
  only for verification, commented to keep files small */ 
//   for(n=0;n<grid.Hsize();n++) {
//     for(k=0;k<nlayers;k++) {
//       size_t m=k*grid.Hsize()+n;
//       if(u[n][k]!=umask) pbarotropic[m]=u[n][k]/9.81/1019.2;
//       else pbarotropic[m]=umask;
//       }
//     }
//   create_file=0;
//   =1;
//   status=save_SGXYZ_C("modal-decomposition.nc", grid, pbarotropic, umask, "P_a", "P_G", "m", "u_coeffcients", create_file, create_grid);
 
  return(0);
}


