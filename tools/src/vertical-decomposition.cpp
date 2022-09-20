
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

extern  void compute_atlas_names(const string & wave,const string & varname,string *filename,string *varnameA,string *varnameG);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double * get_layer_thickness(const string & rootname,const grid_t & grid,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,n;
  grid_t wgrid;
  const string
    filename=rootname+".vertical-modes.nc";
  
  status=poc_get_grid(filename, "Wmodes", &wgrid, verbose, 0);
  NC_TRAP_ERROR(wexit,status,1,"poc_get_grid(\""+filename+"\",\"Wmodes\",...) error");
  const int
    s=grid.size(),
    gridHSize=grid.Hsize();
  
  if(wgrid.Hsize()!=gridHSize) TRAP_ERR_EXIT(ENOEXEC,"not coded for that\n");
  double *layer_thickness=0;
  swapval(layer_thickness,wgrid.z);
  wgrid.free();
  for(n=0;n<s;n++){
    layer_thickness[n]=fabs(layer_thickness[n]-layer_thickness[n+gridHSize]);
    }
  
  return(layer_thickness);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int VerticalModes_SpectralDecomposition(const string & filename,const string & varname, int maxmodes,
                                          double ** &modes, double & mask, complex<double> **p, complex<double> umask, complex<double> ** &p_decomposition, int & MinUmodes, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,/*k,*/n,status, verbose=0;
  int nmasked, nprocs;
  grid_t grid;
  const string rootname=filename.substr(0,filename.find('.'));
  
  poc_data_t<double> scalarData;
  
/*------------------------------------------------------------------------------
   */
//   filename="NEMO.vertical-modes.nc";
//   varname="Pmodes";
  
  status=poc_get_grid(filename, varname, &grid, verbose, -1);
  
  const int
    nlayers=grid.nz,
    gridHSize=grid.Hsize();
  
  double *layer_thickness;
  layer_thickness=get_layer_thickness(rootname,grid,verbose);
  
  MinUmodes=nlayers;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  load pre-computed vertical modes
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  status=scalarData.init(filename,varname);
  
  for(i=0;i<MinUmodes;i++) {
    status=scalarData.read_data(filename,i);
    MinUmodes=min(scalarData.nframes, maxmodes);
    if(modes==0) {
      modes=aset(MinUmodes,(double*)0);
      }
    if(modes[i]==0) modes[i]=new double[grid.size()];
    for(n=0;n<grid.size();n++) modes[i][n]=scalarData.data[n];
    mask=scalarData.mask;
    nmasked=occurence<int>(mask, modes[i], grid.size());
    }
  
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  compute modal decomposition
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  p_decomposition=new complex<double>*[MinUmodes];
  for(i=0;i<MinUmodes;i++) {
    p_decomposition[i]=aset(gridHSize,umask);
    }
  
  int nRequestedProcs=-1;
  nprocs=initialize_OPENMP(nRequestedProcs);
  
// #pragma omp parallel for private(status) if(nprocs>1)
//   for(int jj=2;jj<grid.ny-2;jj++) {
  for(int jj=0;jj<grid.ny;jj++) {
    complex<double> *decomposition=new complex<double>[nlayers];
    double **lmodes=new double*[MinUmodes];
    for(int i=0;i<MinUmodes;i++) {
      lmodes[i]=new double[nlayers];
      }
    for(int ii=2;ii<grid.nx-2;ii++) {
      size_t n=jj*grid.nx+ii;
      int ntruelayers=0;
      int ntruemodes=0;
      int count;
      bool skip=false;
      string skipReason;
      for(int k=0;k<nlayers;k++) {
        if(p[n][k]==umask) continue;
        size_t m=k*gridHSize+n;
        if(layer_thickness[m]<=0.) continue;
        ntruelayers++;
        }
      for(int i=0;i<MinUmodes;i++) {
        double std=0;
        count=0;
        for(int k=0;k<nlayers;k++) {
          size_t m=k*gridHSize+n;
          lmodes[i][k]=modes[i][m];
          if(lmodes[i][k]!=mask) {
            std+=square(lmodes[i][k]);
            count++;
            }
          }
        if(count!=0 and std!=0.0) {
/*------------------------------------------------------------------------------
          count effective number of valid */        
          ntruemodes++;
/*------------------------------------------------------------------------------
          check added because of tracer and velocity grids mismatches */        
          if(count!=ntruelayers) {
            skip=true;
            asprintf(skipReason,"[%d]%d!=%d\n",i,count,ntruelayers);
            }
          std=sqrt(std/count);
/*------------------------------------------------------------------------------
          assume normalisation already done */        
//           for(int k=0;k<nlayers;k++) {
//             size_t m=k*gridHSize+n;
//             if(lmodes[i][k]!=mask) lmodes[i][k]/=std;
//             if(modes[i][m]!=mask) modes[i][m]/=std;
//             }
          }
        else {
          if(i==0){
            skip=true;
            asprintf(skipReason,"[%d] count=%d , std=%g\n",i,count,std);
            }
          }
        }
/*------------------------------------------------------------------------------
      check added because of tracer and velocity grids mismatches */
      int nmodes;
      if(MinUmodes>ntruelayers) {
        nmodes=min(ntruemodes,ntruelayers);
        }
      else {
        nmodes=min(ntruemodes,MinUmodes);
        }
      if(skip){
//         STDERR_BASE_LINE_FUNC("[%d,%d] skipped:\n"+skipReason+"\n",ii,jj);
        continue;
        }
      status=Umodes_decomposition_v2(ntruelayers, lmodes, nmodes, p[n], decomposition, debug, 0);
      count=0;
      for(int i=0;i<MinUmodes;i++) {
        p_decomposition[i][n]=decomposition[i];
        }
      }
    deletep2D(&lmodes,MinUmodes);
    delete[] decomposition;
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int VerticalModes_SpectralDecomposition(const string & rootname, int maxmodes, string & wave,const string varnames[5], int nrecomposed, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  Interface to decompose harmonic fields into modal contributions
  
  Presently hard-coded for NEMO-derived files, badly-formed structured
  ouputs are just a pain in the ... as usual
  
    v f v f v f
    t u t u t u   last T line unusable as v (above) masked
    v f v f v f
    t u t u t u   first T line unusable as v start after
    
  
  u,v and t have same dimensions: 2 lines and columns too much in t

------------------------------------------------------------------------------*/
{
  int i,k,n,status;
  grid_t grid,ugrid,vgrid;
  double **Umodes=0, **Pmodes=0, mask;
  int MinUmodes;
  complex<double> **u=0, **v=0, umask, *buffer=0, cmask;
  complex<double> **u_decomposition=0, **v_decomposition=0, **p_decomposition=0;
  complex<double> *ubarotropic=0, *ubaroclinic=0, *vbarotropic=0, *vbaroclinic=0, *pbarotropic=0, *pbaroclinic=0;
  string filename,varnameA,varnameG,varname;
  int verbose=0;
  int range[2]={0,-1};
  int create_file=1, create_grid=1;
  int nkept;
  
  nkept=min(maxmodes, nrecomposed);
  
  create_file=1;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  horizontal velocity modal decomposition
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
/*------------------------------------------------------------------------------
  load horizontal velocity vertical modes */
  
  filename=wave+"-"+varnames[3]+"-atlas.nc";
  status=poc_get_grid(filename,varnames[3]+"_a",&grid,-1,0);
  if(status!=0){
    NC_CHKERR_BASE_LINE(status,"poc_get_grid(\""+filename+"\",\""+varnames[3]+"_a\",...) error");
    STDERR_BASE_LINE_FUNC("getting grid from z atlas\n");
    
    varname="z";
    varnameA=varname+"_a";
    varnameG=varname+"_G";
    filename=wave+"-"+varname+"-atlas.nc";
    
    status=poc_get_grid(filename,varnameA,&grid,verbose,0);
    NC_TRAP_ERROR(wexit,status,1,"poc_get_grid(\""+filename+"\",\""+varnameA+"\",...) error");
    }
  
  const int
    nlayers=grid.nz,
    gridHSize=grid.Hsize();
  
  umask=NC_FILL_COMPLEX;
  u=new complex<double>* [gridHSize];
  v=new complex<double>* [gridHSize];
  for(n=0;n<gridHSize;n++) {
    u[n]=aset(nlayers,umask);
    v[n]=aset(nlayers,umask);
    }
  
/*------------------------------------------------------------------------------
  zonal velocity, load and interpolate at tracer nodes */
  filename=wave+"-"+varnames[1]+"-atlas.nc";
  varnameA=varnames[1]+"_a";
  varnameG=varnames[1]+"_G";
  
  status=poc_get_grid(filename, varnameA, &ugrid, verbose, -1);
  if(status!=0) goto pmodes;
  
  buffer=new complex<double>[ugrid.size()];
  
  status=poc_get_cvara(filename,varnameA, varnameG,-1, buffer);
  
//   for(int j=2;j<grid.ny-2;j++) {
  for(int j=0;j<grid.ny;j++) {
    for(i=2;i<grid.nx-2;i++) {
      size_t n=j*grid.nx+i;
      for(k=0;k<nlayers;k++) {
/*------------------------------------------------------------------------------
        u-node ndex */
        size_t m1=k*ugrid.Hsize()+j*ugrid.nx+i-1;
        size_t m2=k*ugrid.Hsize()+j*ugrid.nx+i;
        if(buffer[m1]!=umask and buffer[m2]!=umask) {
          u[n][k]=0.5*(buffer[m1]+buffer[m2]);
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
  
/*------------------------------------------------------------------------------
  meridian velocity, load and interpolate at tracer nodes */
  filename=wave+"-"+varnames[2]+"-atlas.nc";
  varnameA=varnames[2]+"_a";
  varnameG=varnames[2]+"_G";
  
  status=poc_get_grid(filename, varnameA, &vgrid, verbose, -1);
  if(status!=0) goto pmodes;
  
  buffer=new complex<double>[vgrid.size()];
  
  status=poc_get_cvara(filename,varnameA, varnameG,-1, buffer);
  
  if(grid.ny!=1) {
/*------------------------------------------------------------------------------
    regular grid */
    for(int j=1;j<grid.ny;j++) {
      for(int i=2;i<grid.nx-2;i++) {
        size_t n=j*grid.nx+i;
        for(k=0;k<nlayers;k++) {
          size_t m1=k*vgrid.Hsize()+(j-1)*vgrid.nx+i;
          size_t m2=k*vgrid.Hsize()+j*vgrid.nx+i;
          if(buffer[m1]!=umask and buffer[m2]!=umask) {
            v[n][k]=0.5*(buffer[m1]+buffer[m2]);
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
        size_t m1=k*vgrid.Hsize()+j*vgrid.nx+i;
        size_t m2=k*vgrid.Hsize()+j*vgrid.nx+i;
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
  
  filename=rootname+".vertical-modes.nc";
  varname="Umodes";
  
  status=VerticalModes_SpectralDecomposition(filename, varname, maxmodes, Umodes, mask, u, umask, u_decomposition, MinUmodes, debug);
  status=VerticalModes_SpectralDecomposition(filename, varname, maxmodes, Umodes, mask, v, umask, v_decomposition, MinUmodes, debug);
  
  create_grid=1;
  range[0]=0;
  range[1]=MinUmodes-1;
/*------------------------------------------------------------------------------
  modal decomposition coefficients */
  status=save_SGXYT_C("modal-decomposition.nc", grid, u_decomposition, umask, "Uc_a", "Uc_G", "dimensionless", "u_coeffcients", create_file, create_grid, range);
  create_file=0;
  create_grid=1;
  status=save_SGXYT_C("modal-decomposition.nc", grid, v_decomposition, umask, "Vc_a", "Vc_G", "dimensionless", "u_coeffcients", create_file, create_grid, range);
  
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
//   status=save_SGXYZ_C("modal-decomposition.nc", grid, ubarotropic, umask, "Ubt_a", "Ubt_G", "m/s", "u_barotropic", create_file, create_grid);
//   create_file=0;
//   status=save_SGXYZ_C("modal-decomposition.nc", grid, ubaroclinic, umask, "Ubc_a", "Ubc_G", "m/s", "u_baroclinic", create_file, create_grid);
//   status=save_SGXYZ_C("modal-decomposition.nc", grid, vbarotropic, umask, "Vbt_a", "Vbt_G", "m/s", "v_barotropic", create_file, create_grid);
//   status=save_SGXYZ_C("modal-decomposition.nc", grid, vbaroclinic, umask, "Vbc_a", "Vbc_G", "m/s", "v_baroclinic", create_file, create_grid);
  
  for(i=0;i<nkept;i++) {
    for(n=0;n<gridHSize;n++) {
      for(int k=0;k<nlayers;k++) {
        size_t m=k*gridHSize+n;
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
    status=save_SGXYZT_C("modal-recomposition.nc", grid, ubaroclinic, umask, "u_modal_a","u_modal_G","m/s", "u_modal", create_file, create_grid,"",frame);
    create_file=0;
    create_grid=0;
    status=save_SGXYZT_C("modal-recomposition.nc", grid, vbaroclinic, umask, "v_modal_a","v_modal_G","m/s", "v_modal", create_file, create_grid,"",frame);
    }

/*------------------------------------------------------------------------------
  re-compose horizontal velocity barotropic and baroclinic modes */
  i=0;
  for(n=0;n<gridHSize;n++) {
    for(int k=0;k<nlayers;k++) {
      size_t m=k*gridHSize+n;
      ubarotropic[m]=umask;
      vbarotropic[m]=umask;
      if(Umodes[i][m]==mask) continue;
      if(u_decomposition[i][n]!=umask) ubarotropic[m]=Umodes[i][m]*u_decomposition[i][n];
      if(v_decomposition[i][n]!=umask) vbarotropic[m]=Umodes[i][m]*v_decomposition[i][n];
      }
    }
  create_file=0;
  create_grid=0;
  status=save_SGXYZ_C("modal-recomposition.nc", grid, ubarotropic, umask, "u_barotropic_a","u_barotropic_G","m/s", "u_barotropic", create_file, create_grid,"","");
  status=save_SGXYZ_C("modal-recomposition.nc", grid, vbarotropic, umask, "v_barotropic_a","v_barotropic_G","m/s", "v_barotropic", create_file, create_grid,"","");
  
  i=1;
  for(n=0;n<gridHSize;n++) {
    for(int k=0;k<nlayers;k++) {
      size_t m=k*gridHSize+n;
      ubaroclinic[m]=umask;
      vbaroclinic[m]=umask;
      if(Umodes[i][m]==mask) continue;
      if(u_decomposition[i][n]!=umask) ubaroclinic[m]=Umodes[i][m]*u_decomposition[i][n];
      if(v_decomposition[i][n]!=umask) vbaroclinic[m]=Umodes[i][m]*v_decomposition[i][n];
      }
    }
  for(i=2;i<MinUmodes-1;i++) {
    for(n=0;n<gridHSize;n++) {
      for(int k=0;k<nlayers;k++) {
        size_t m=k*gridHSize+n;
        if(ubaroclinic[m]==umask) continue;
        if(Umodes[i][m]==mask) continue;
        if(u_decomposition[i][n]!=umask) ubaroclinic[m]+=Umodes[i][m]*u_decomposition[i][n];
        if(v_decomposition[i][n]!=umask) vbaroclinic[m]+=Umodes[i][m]*v_decomposition[i][n];
        }
      }
    }
  create_file=0;
  create_grid=0;
  status=save_SGXYZ_C("modal-recomposition.nc", grid, ubaroclinic, umask, "u_baroclinic_a","u_baroclinic_G","m/s", "u_baroclinic", create_file, create_grid,"","");
  status=save_SGXYZ_C("modal-recomposition.nc", grid, vbaroclinic, umask, "v_baroclinic_a","v_baroclinic_G","m/s", "v_baroclinic", create_file, create_grid,"","");

  create_file=0;


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  pressure modal decomposition
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
pmodes:
  
  int k0,kmax,dk;
  bool increasing;
  double factor;
  const double rho0=999.2,g=9.81;
  
  status=check_vertical_direction(grid,&increasing,&factor);
  if(status!=0) TRAP_ERR_EXIT(status,"check_vertical_direction() error\n");
  if(1 or verbose) STDOUT_BASE_LINE("increasing=%d,factor=%g\n",(int)increasing,factor);
  
  fflush(stdout);
  
  if(increasing){
    k0=0;
    kmax=grid.nz;
    dk=1;
    }
  else{
    k0=grid.nz-1;
    kmax=-1;
    dk=-1;
    }
  
  umask=NC_FILL_COMPLEX;
  
  compute_atlas_names(wave,varnames[3],&filename,&varnameA,&varnameG);
  buffer=new complex<double>[grid.size()];
  status=poc_get_cvara(filename,varnameA, varnameG,0, buffer, -1);
  
  if(status!=0){
/*------------------------------------------------------------------------------
    isopycnal/Lagrangian */
    NC_CHKERR_BASE_LINE(status,"poc_get_cvara(\""+filename+"\",\""+varnameA+"\",\""+varnameG+"\",0,) error");
    STDERR_BASE_LINE("computing isopycnal/Lagrangian pressure.\n");
    
    compute_atlas_names(wave,"z",&filename,&varnameA,&varnameG);
    status=poc_get_cvara(filename,varnameA, varnameG,0, buffer);
    NC_TRAP_ERROR(wexit,status,1,"poc_get_cvara(\""+filename+"\",\""+varnameA+"\",\""+varnameG+"\",0,) error");
    
    poc_data_t<double> sigma;
    status=sigma.init("grid.nc","sigma");
    NC_TRAP_ERROR(wexit,status,1,"poc_data_t::init(\"grid.nc\",\"sigma\") error");
    
    poc_cdata_t layer_thickness;
    compute_atlas_names(wave,"h",&filename,&varnameA,&varnameG);
    status=layer_thickness.init(filename,varnameA,varnameG);
    NC_TRAP_ERROR(wexit,status,1,"poc_cdata_t::init(\""+filename+"\",\""+varnameA+"\",\""+varnameG+"\",0,) error");
    
    for(n=0;n<gridHSize;n++){
      complex<double> mass(0.,0.),p;/* integrated mass */
      
      for(k=k0;k!=kmax;k+=dk) {
        size_t m=k*gridHSize+n;
        
        const complex<double>
          *bufferm=&buffer[m],/* z' */
          *layer_thicknessm=&layer_thickness.data[m];
        const double
          sigmam=sigma.data[m],
          rhom=sigma.data[m]+1e3;
        
        if(*bufferm!=umask and sigmam!=sigma.mask and *layer_thicknessm!=umask) {
          
          p=(mass-(*bufferm)*rhom)*g;
          
          u[n][k]=p;
          
          mass+=(*layer_thicknessm)*rhom;
          }
        else {
          u[n][k]=umask;
          }
        
        if(n==448){
          printf("[%d][%d]%11gm+%11gm,%11g° %11gm,%11g° %11gkg/m3 %11gkg,%11g° %11gP,%11g°\n",n,k,grid.z[m],abs(*bufferm),arg(*bufferm),abs(*layer_thicknessm),arg(*layer_thicknessm),rhom,abs(mass),arg(mass),abs(p),arg(p));
          if(*bufferm==umask and k==k0)
            TRAP_ERR_EXIT(ENOEXEC,"wrong test at %d\n",n);
          }
        }
      }
    
    }
  else{
/*------------------------------------------------------------------------------
    sigma/Eulerian */
    
    complex<double> *ssh;
    ssh=new complex<double>[gridHSize];
    compute_atlas_names(wave,varnames[4],&filename,&varnameA,&varnameG);
    status=poc_get_cvara(filename,varnameA, varnameG,0, ssh);
    NC_TRAP_ERROR(wexit,status,1,"poc_get_cvara(\""+filename+"\",\""+varnameA+"\",\""+varnameG+"\",-1,) error");
    
    double *layer_thickness;
    layer_thickness=get_layer_thickness(rootname,grid,verbose);
    
    for(n=0;n<gridHSize;n++) {
      complex<double> p(0.,0.); /* pressure */
      double dh0,dh1;           /* layer thicknesses */
      complex<double> lm0,lm1;  /* layer masses */
      for(k=k0;k!=kmax;k+=dk) {
        size_t m=k*gridHSize+n;
        complex<double> *bufferm=&buffer[m];/* rho' */
        
        dh0=dh1;
        dh1=layer_thickness[m];
        
        if(*bufferm!=umask) {
          if(k==k0){
            lm1=dh1* *bufferm;
            p=g*(rho0*ssh[n]+0.5*lm1);
            }
          else {
            lm0=lm1;
            lm1=dh1* *bufferm;
            p+=g*0.5*(lm0+lm1);
            }
          
          if(dh1<0.) TRAP_ERR_EXIT(ENOEXEC,"[%d][%d] dh1=%g<=0\n",n,k,dh1);
          
          u[n][k]=p;
          }
        else {
          u[n][k]=umask;
          }
        
        if(n==448){
          printf("[%d][%d]%11gm %11gm %11gkg/m3,%11g° %11gkg,%11g° %11gP,%11g°\n",n,k,grid.z[m],dh1,abs(*bufferm),arg(*bufferm),abs(lm1),arg(lm1),abs(p),arg(p));
          if(*bufferm==umask and k==k0)
            TRAP_ERR_EXIT(ENOEXEC,"wrong test at %d\n",n);
          }
        }
      }
    
    delete[] layer_thickness;
    delete[] ssh;
    }
  fflush(stdout);
  delete[] buffer;
  
  pbarotropic=new complex<double>[grid.size()];
  for(size_t m=0;m<grid.size();m++) pbarotropic[m]=umask;
  pbaroclinic=new complex<double>[grid.size()];
  for(size_t m=0;m<grid.size();m++) pbaroclinic[m]=umask;
  
  
  filename=rootname+".vertical-modes.nc";
  varname="Pmodes";
  
  status=VerticalModes_SpectralDecomposition(filename, varname, maxmodes, Pmodes, mask, u, umask, p_decomposition, MinUmodes, debug);
    
  create_file=0;
  create_grid=1;
  range[0]=0;
  range[1]=MinUmodes-1;
/*------------------------------------------------------------------------------
  modal decomposition coefficients */
  status=save_SGXYT_C("modal-decomposition.nc", grid, p_decomposition, umask, "Pc_a", "Pc_G", "dimensionless", "u_coeffcients", create_file, create_grid, range);

//   create_file=0;
//   create_grid=1;
// /*------------------------------------------------------------------------------
//   modal fields */
//   status=save_SGXYZ_C("modal-decomposition.nc", grid, pbarotropic, umask, "Pbt_a", "Pbt_G", "N/m²", "p_barotropic", create_file, create_grid);
//   create_file=0;
//   create_grid=1;
//   status=save_SGXYZ_C("modal-decomposition.nc", grid, pbaroclinic, umask, "Pbc_a", "Pbc_G", "N/m²", "p_baroclinic", create_file, create_grid);
  
  for(i=0;i<nkept;i++) {
    for(n=0;n<gridHSize;n++) {
      for(int k=0;k<nlayers;k++) {
        size_t m=k*gridHSize+n;
        pbaroclinic[m]=umask;
        if(Pmodes[i][m]==mask) continue;
        if(p_decomposition[i][n]!=umask) pbaroclinic[m]=Pmodes[i][m]*p_decomposition[i][n];
        }
      }
    create_file=0;
    create_grid=0;
    int frame=i;
    status=save_SGXYZT_C("modal-recomposition.nc", grid, pbaroclinic, umask, "p_modal_a","p_modal_G","N/m²", "p_modal", create_file, create_grid,"",frame);
    }
  
/*------------------------------------------------------------------------------
  re-compose pressure barotropic and baroclinic modes */
  i=0;
  for(n=0;n<gridHSize;n++) {
    for(int k=0;k<nlayers;k++) {
      size_t m=k*gridHSize+n;
      pbarotropic[m]=umask;
      if(Pmodes[i][m]==mask) continue;
      if(p_decomposition[i][n]!=umask) pbarotropic[m]=Pmodes[i][m]*p_decomposition[i][n];
      }
    }
  create_file=0;
  create_grid=0;
  status=save_SGXYZ_C("modal-recomposition.nc", grid, pbarotropic, umask, "p_barotropic_a","p_barotropic_G","m/s", "p_barotropic", create_file, create_grid,"","");
    
  i=1;
  for(n=0;n<gridHSize;n++) {
    for(int k=0;k<nlayers;k++) {
      size_t m=k*gridHSize+n;
      pbaroclinic[m]=umask;
      if(Pmodes[i][m]==mask) continue;
      if(p_decomposition[i][n]!=umask) pbaroclinic[m]=Pmodes[i][m]*p_decomposition[i][n];
      }
    }
  for(i=2;i<MinUmodes-1;i++) {
    for(n=0;n<gridHSize;n++) {
      for(int k=0;k<nlayers;k++) {
        size_t m=k*gridHSize+n;
        if(pbaroclinic[m]==umask) continue;
        if(Pmodes[i][m]==mask) continue;
        if(p_decomposition[i][n]!=umask) pbaroclinic[m]+=Pmodes[i][m]*p_decomposition[i][n];
        }
      }
    }
  create_file=0;
  create_grid=0;
  status=save_SGXYZ_C("modal-recomposition.nc", grid, pbaroclinic, umask, "p_baroclinic_a","p_baroclinic_G","m/s", "p_baroclinic", create_file, create_grid,"","");

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int VerticalModes_SnapshotDecomposition(const string *input, const string & filename, const string *variables, int maxmodes, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Interface to decompose temporal fields into modal contributions

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i,k,n,status;
  grid_t grid, ugrid, vgrid;
  double **modes=0, **Umodes=0, **Pmodes=0, mask;
  int MinUmodes;
  double *decomposition=0, **u=0, **v=0, **p=0, umask, *buffer=0;
  double **u_decomposition=0, **v_decomposition=0, **p_decomposition=0;
  typedef float TYPE;
  TYPE **u_modal=0, **v_modal=0, **p_modal=0, **ssh_modal, TYPEmask;
  float *dz_t;
  string varname;
  int nlayers;
  int verbose=0;
  poc_data_t<float> scalarData;
  string gridfile, ufile;
  float *P=0, *R=0, *T=0, *S=0, *ssh=0;
  
  string uname, vname, wname, rname, tname, sname, ssh_name, dz_name, EP_name;
  
#if 1  
   uname="uoce"; vname="voce"; wname="woce"; rname="rhop"; tname=""; sname=""; ssh_name="ssh"; dz_name="e3t"; EP_name="eulpres";
#else
   uname="u"; vname="v"; wname="woce"; rname="rhop"; tname="tem"; sname="sal", ssh_name="ssh_ib"; dz_name="dz_t"; EP_name="";
#endif
   
  gridfile=input[G_ID];
  ufile=input[U_ID];
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  3D currents 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  varname="Umodes";
  
  status=poc_get_grid(filename, varname, &grid, verbose, -1);
  
  MinUmodes=grid.nt;
  nlayers=grid.nz;
  
  MinUmodes=nlayers;
  
  int nmasked;
  status=scalarData.init(filename,varname);
  for(i=0;i<MinUmodes;i++) {
    status=scalarData.read_data(filename,i);
    MinUmodes=min(scalarData.nframes, maxmodes);
    if(Umodes==0) Umodes=new double*[MinUmodes];
    Umodes[i]=new double[grid.size()];
    for(n=0;n<grid.size();n++) Umodes[i][n]=scalarData.data[n];
    mask=scalarData.mask;
    nmasked=occurence<int>(mask, Umodes[i], grid.size());
    }
  
  modes=new double*[MinUmodes];
  for(i=0;i<MinUmodes;i++) {
    modes[i]=new double[nlayers];
    }

  poc_var_t vU,vV,vT,vS,vR,vEP,vSSH,vdz;

/*------------------------------------------------------------------------------
  identify variable in data file, get grid in grid file */
  status=poc_inq_var(ufile, uname, &vU, verbose);
  status=poc_get_grid(gridfile, vU, &ugrid, verbose, -1);

  status=poc_inq_var(ufile, vname, &vV, verbose);
  status=poc_get_grid(gridfile, vV, &vgrid, verbose, -1);
  
     
/*------------------------------------------------------------------------------
  identify layer thickness in grid file */
  dz_t=new float[grid.size()];
  status=poc_inq_var(ufile, dz_name, &vdz, verbose);
  status=poc_get_vara(gridfile, vdz, 0, dz_t);

  status=poc_inq_var(ufile, ssh_name, &vSSH, verbose);
  ssh=new float[grid.Hsize()];
  status=poc_get_vara(ufile, vSSH, 0, ssh, 0);
  
/*------------------------------------------------------------------------------
  compute density, first step toward pressure */
  
  R=new float[grid.size()];
  
  if(tname!="" and sname!="") {
    
    status=poc_inq_var(ufile, tname, &vT, verbose);
    status=poc_inq_var(ufile, sname, &vS, verbose);
    
    T=new float[grid.size()];
    status=poc_get_var(ufile, tname, T);
  
    S=new float[grid.size()];
    status=poc_get_var(ufile, sname, S);
  
    for(k=0;k<grid.nz;k++) {
      for(int j=0;j<grid.ny;j++) {
        for(int i=0;i<grid.nx;i++) {
          size_t m=k*grid.ny*grid.nx+j*grid.nx+i;
          if(T[m]!=mask) {
            if(T[m]<0.0) {
              if(verbose) printf("T anomaly at i=%4d j=%4d k=%4d : T=%f S=%f \n",i,j,k,T[m],S[m]);
              T[m]=4.;
              }
            R[m]=water_density(T[m],S[m],grid.z[m]);
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
    }
  else {
    status=poc_inq_var(ufile, rname, &vR, verbose);
    status=poc_get_vara(ufile, vR, 0, R, 0);
    }
  
/*------------------------------------------------------------------------------
  vertical modes need pressure perturbation (hence density perturbation) */
  P=new float[grid.size()];
  for(size_t m=0; m<grid.size(); m++) {
    if(R[m]!=mask) P[m]=0.0;
    else P[m]=mask;
    }
  
/*------------------------------------------------------------------------------
  brutal global average */
  statistic_t s=get_statistics(R,mask,grid.size(),1);
  
/*------------------------------------------------------------------------------
  depth-dependant average at reference depths; perfectible */
  vector<double> z;
  double *Rref, *count;
  
// #define SMART

#ifdef SMART  
  for(double zz=0;zz>-100.;zz-=5.) {
    z.push_back(zz);
    }
  for(double zz=-100;zz>-500.;zz-=10.) {
    z.push_back(zz);
    }
  for(double zz=-500;zz>-1000.;zz-=20.) {
    z.push_back(zz);
    }
  for(double zz=-1000;zz>-3000.;zz-=50.) {
    z.push_back(zz);
    }
  
  const size_t zn=z.size();
  
  Rref=aset(zn,0.);
  count=aset(zn,0.);
  
  for(k=0;k<grid.nz;k++) {
    for(int j=0;j<grid.ny;j++) {
      for(int i=0;i<grid.nx;i++) {
        size_t m=k*grid.ny*grid.nx+j*grid.nx+i;
        if(R[m]!=mask) {
          double zz=grid.z[m];
          int pos=vpos(zz, z);
          if(pos<0 or zn<=pos) TRAP_ERR_EXIT(ENOEXEC,"programming error : wrong index %d/%u\n",pos,zn);
          Rref[pos]+=R[m]-1000.0;
          count[pos]+=1;
          if(pos!=0) pos--;
          Rref[pos]+=R[m]-1000.0;
          count[pos]+=1;
          }
        }
      }
    }
  for(k=0;k<z.size();k++) {
    Rref[k]/=count[k];
    Rref[k]+=1000.;
    }

/*------------------------------------------------------------------------------
  sort density in ascending order */
  int nswap;
  do {
    nswap=0;
    for(k=0;k<z.size()-1;k++) {
      if(Rref[k]>Rref[k+1]) {
        double tmp=Rref[k];
        Rref[k]=Rref[k+1];
        Rref[k+1]=tmp;
        nswap++;
        }
      }
    } while(nswap!=0);

  FILE *out=fopen("profile.dat","w");
  fprintf(out,"%d\n", z.size());
  for(k=0;k<z.size();k++) {
    double N2;
    if(k==0 or k==z.size()-1) N2=0;
    else N2=sqrt(-9.81*(Rref[k+1]-Rref[k-1])/(z[k+1]-z[k-1])/Rref[k]);
    fprintf(out,"%lf %lf %lf\n", z[k], Rref[k], N2);
    }
  fclose(out);
#endif
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  compute perturbation pressure
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// #define SMART
  k=grid.nz-1;
  for(size_t n=0; n<grid.Hsize(); n++) {
    size_t m=k*grid.Hsize()+n;
    if(R[m]!=mask and ssh[n]!=mask) {
      double mean=s.mean;
#ifdef SMART
      double zz=grid.z[m];
      int pos=vpos(zz, z);
      mean=Rref[pos];
#endif
      double RR=R[m]-mean;
      P[m]=9.81*(R[m]*ssh[n]+RR*dz_t[m]/2.0);
      }
    else {
      P[m]=mask;
      }
    }
  
  for(k=grid.nz-2;k>=0;k--) {
    for(int j=0;j<grid.ny;j++) {
      for(int i=0;i<grid.nx;i++) {
        size_t mup=(k+1)*grid.ny*grid.nx+j*grid.nx+i;
        size_t m=k*grid.ny*grid.nx+j*grid.nx+i;
        if(R[m]!=mask) {
          double mean=s.mean;
 #ifdef SMART
         double zz=grid.z[m];
          int pos=vpos(zz, z);
          mean=0.5*Rref[pos];
          zz=grid.z[mup];
          pos=vpos(zz, z);
          mean+=0.5*Rref[pos];
#endif
          P[m]=P[mup]+9.81*((R[mup]-mean)*dz_t[mup]/2.0+(R[m]-mean)*dz_t[m]/2.0);
          }
        else {
          P[m]=mask;
          }
        }
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  snapshot modal decomposition, input arrays processing
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  umask=-9999.;
  
  u=new double* [grid.Hsize()];
  v=new double* [grid.Hsize()];
  p=new double* [grid.Hsize()];
  for(n=0;n<grid.Hsize();n++) {
    u[n]=new double[nlayers];
    v[n]=new double[nlayers];
    p[n]=new double[nlayers];
    for(k=0;k<nlayers;k++) {
      u[n][k]=umask;
      v[n][k]=umask;
      p[n][k]=umask;
      }
    }
  
  buffer=new double[ugrid.size()];
  status=poc_get_vara(ufile, vU, 0, buffer, 0);
  
  for(int j=1;j<grid.ny-1;j++) {
    for(i=1;i<grid.nx-1;i++) {
      size_t n=j*grid.nx+i;
      for(k=0;k<nlayers;k++) {
        size_t m1=k*ugrid.Hsize()+j*ugrid.nx+i-1;
        size_t m2=k*ugrid.Hsize()+j*ugrid.nx+i;
        if(buffer[m1]!=umask and buffer[m2]!=umask) {
          u[n][k]=0.5*(buffer[m1]+buffer[m2]);
          if(u[n][k]==0) u[n][k]=umask;
          }
        else  {
          u[n][k]=umask;
          }
        }
      }
    }
  delete[] buffer;
  
  buffer=new double[vgrid.size()];
  status=poc_get_vara(ufile, vV, 0, buffer, 0);
  
  for(int j=1;j<grid.ny-1;j++) {
    for(i=1;i<grid.nx-1;i++) {
      size_t n=j*grid.nx+i;
      for(k=0;k<nlayers;k++) {
        size_t m1=k*vgrid.Hsize()+(j-1)*vgrid.nx+i;
        size_t m2=k*vgrid.Hsize()+j*vgrid.nx+i;
        if(buffer[m1]!=umask and buffer[m2]!=umask) {
          v[n][k]=0.5*(buffer[m1]+buffer[m2]);
          if(v[n][k]==0) v[n][k]=umask;
          }
        else  {
          v[n][k]=umask;
          }
        }
      }
    }
  delete[] buffer;
  
  for(int j=1;j<grid.ny-1;j++) {
    for(i=1;i<grid.nx-1;i++) {
      size_t n=j*grid.nx+i;
      for(k=0;k<nlayers;k++) {
        size_t m=k*grid.Hsize()+j*grid.nx+i;
        p[n][k]=P[m];
        }
      }
    }
    
  decomposition=new double[nlayers];
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  3D currents decomposition
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  u_decomposition=new double*[MinUmodes];
  v_decomposition=new double*[MinUmodes];
  p_decomposition=new double*[MinUmodes];
  for(i=0;i<MinUmodes;i++) {
    u_decomposition[i]=aset(grid.Hsize(),umask);
    v_decomposition[i]=aset(grid.Hsize(),umask);
    p_decomposition[i]=aset(grid.Hsize(),umask);
    }
  
  TYPEmask=umask;
  
  int nkept=5;
  
  nkept=MinUmodes;
  
  u_modal=new TYPE *[nkept];
  for(i=0;i<nkept;i++) {
    u_modal[i]=new TYPE[grid.size()];
    for(size_t m=0;m<grid.size();m++) u_modal[i][m]=TYPEmask;
    }
  
  v_modal=new TYPE *[nkept];
  for(i=0;i<nkept;i++) {
    v_modal[i]=new TYPE[grid.size()];
    for(size_t m=0;m<grid.size();m++) v_modal[i][m]=TYPEmask;
    }
  
  for(n=0;n<grid.Hsize();n++) {
    int ntruelayers=0;
    int ntruemodes=MinUmodes;
    double count;
    bool skip=false;
    for(k=0;k<nlayers;k++) {
      if(u[n][k]!=umask) ntruelayers++;
      }
/*------------------------------------------------------------------------------
    euclidian normalization of modes (so we can compare coefficients) */
    ntruemodes=0;       
    for(i=0;i<MinUmodes;i++) {
      double std=0;
      count=0;
      for(k=0;k<nlayers;k++) {
        size_t m=k*grid.Hsize()+n;
        modes[i][k]=Umodes[i][m];
        if(modes[i][k]!=mask) {
          std+=modes[i][k]*modes[i][k]*dz_t[m];
//           count+=dz_t[m];
          count++;
          }
        }
      if(count) ntruemodes++;
//       if(count!=0 and std!=0.0) {
// /*------------------------------------------------------------------------------
//         count effective number of valid */        
//         ntruemodes++;
//         std=sqrt(std/count);
//         for(k=0;k<nlayers;k++) {
//           if(modes[i][k]!=mask) modes[i][k]/=std;
//           }
//         }
//       else {
//         if(i==0) skip=true;
//         }
      }
/*------------------------------------------------------------------------------
    check added because of tracer and velocity grids mismatches */
    int nmodes;
    if(MinUmodes>ntruelayers) {
      nmodes=min(ntruemodes,ntruelayers);
      }
    else {
      nmodes=min(ntruemodes,MinUmodes);
      }
    if(skip or nmodes==0) continue;
    status=Umodes_decomposition_v2(ntruelayers, modes, nmodes, u[n], decomposition, debug, 0);
    for(i=0;i<MinUmodes;i++) {
      u_decomposition[i][n]=decomposition[i];
      }
    status=Umodes_decomposition_v2(ntruelayers, modes, nmodes, v[n], decomposition, debug, 0);
    for(i=0;i<MinUmodes;i++) {
      v_decomposition[i][n]=decomposition[i];
      }
    for(k=0;k<nlayers;k++) {
      size_t m=k*grid.Hsize()+n;
      for(i=0;i<min(nmodes,nkept-1);i++) {
        if(Umodes[i][m]!=mask) {
          u_modal[i][m]=modes[i][k]*u_decomposition[i][n];
          v_modal[i][m]=modes[i][k]*v_decomposition[i][n];
          }
        }
      int baroclinic=nkept-1;
      u_modal[baroclinic][m]=0;
      v_modal[baroclinic][m]=0;
      for(i=1;i<nmodes;i++) {
        if(Umodes[i][m]!=mask) {
          u_modal[baroclinic][m]+=modes[i][k]*u_decomposition[i][n];
          v_modal[baroclinic][m]+=modes[i][k]*v_decomposition[i][n];
          }
        }
      }
    }
  
  int range[2]={0,MinUmodes-1};
  
  grid.time=new double[MinUmodes];
  grid.nt=MinUmodes;
  
  int create_file=1;
  int create_grid=1;
  for(int k=0;k<MinUmodes;k++) grid.time[k]=k;
  
  status=save_SGXYT("modal-decomposition.nc", grid, u_decomposition, umask, "Uc", "dimensionless", "u_coeffcients", create_file, create_grid, range);
  create_file=0;
  create_grid=0;
  status=save_SGXYT("modal-decomposition.nc", grid, v_decomposition, umask, "Vc", "dimensionless", "v_coeffcients", create_file, create_grid, range);
  
  for(n=0;n<grid.Hsize();n++) {
    for(i=0;i<MinUmodes;i++) {
      if(u_decomposition[i][n]!=umask and v_decomposition[i][n]!=umask) {
        double uu=u_decomposition[i][n];
        double vv=v_decomposition[i][n];
        u_decomposition[i][n]=uu*uu+vv*vv;
        }
      else {
        u_decomposition[i][n]=umask;
        }
      }
    }
  status=save_SGXYT("modal-decomposition.nc", grid, u_decomposition, umask, "Ek", "dimensionless", "Ek", create_file, create_grid, range);

  for(i=0;i<nkept;i++) {
    create_file=0;
    create_grid=(i==0);
    range[0]=i;
    range[1]=i;
    status=save_SGXYZT("modal-decomposition.nc", grid, u_modal, TYPEmask, "u_modal", "m/s", "u_baroclinic", create_file, create_grid,"","",range);
    create_grid=0;
    status=save_SGXYZT("modal-decomposition.nc", grid, v_modal, TYPEmask, "v_modal", "m/s", "v_barotropic", create_file, create_grid,"","",range);
    }
  
  deletep2D(&Umodes,MinUmodes);
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  pressure decomposition
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  varname="Pmodes";
  
//   status=poc_get_grid(filename, varname, &grid, verbose, -1);
  
  MinUmodes=grid.nt;
  nlayers=grid.nz;
  
  MinUmodes=nlayers;
  
  status=scalarData.init(filename,varname);
  for(i=0;i<MinUmodes;i++) {
    status=scalarData.read_data(filename,i);
    MinUmodes=min(scalarData.nframes, maxmodes);
    if(Pmodes==0) Pmodes=new double*[MinUmodes];
    Pmodes[i]=new double[grid.size()];
    for(n=0;n<grid.size();n++) Pmodes[i][n]=scalarData.data[n];
    mask=scalarData.mask;
    nmasked=occurence<int>(mask, Pmodes[i], grid.size());
    }
  
  p_modal=new TYPE *[nkept];
  for(i=0;i<nkept;i++) {
    p_modal[i]=new TYPE[grid.size()];
    for(size_t m=0;m<grid.size();m++) p_modal[i][m]=TYPEmask;
    }
  
  ssh_modal=new TYPE *[nkept];
  for(i=0;i<nkept;i++) {
    ssh_modal[i]=new TYPE[grid.Hsize()];
    for(size_t m=0;m<grid.Hsize();m++) ssh_modal[i][m]=TYPEmask;
    }
  
  for(n=0;n<grid.Hsize();n++) {
    int ntruelayers=0;
    int ntruemodes=0;
    double count;
    bool skip=false;
    for(k=0;k<nlayers;k++) {
//       if(p[n][k]==0) p[n][k]=umask;
      if(p[n][k]!=umask) ntruelayers++;
      }
/*------------------------------------------------------------------------------
    euclidian normalization of modes (so we can compare coefficients) */        
    for(i=0;i<MinUmodes;i++) {
      double std=0;
      count=0;
      for(k=0;k<nlayers;k++) {
        size_t m=k*grid.Hsize()+n;
        modes[i][k]=Pmodes[i][m];
        if(modes[i][k]!=mask) {
          std+=modes[i][k]*modes[i][k]*dz_t[m];
          count+=dz_t[m];
          }
        }
      if(count!=0 and std!=0.0) {
/*------------------------------------------------------------------------------
        count effective number of valid */        
        ntruemodes++;
        std=sqrt(std/count);
        for(k=0;k<nlayers;k++) {
          if(modes[i][k]!=mask) modes[i][k]/=std;
          }
        }
      else {
        if(i==0) skip=true;
        }
      }
/*------------------------------------------------------------------------------
    check added because of tracer and velocity grids mismatches */
    int nmodes;
    if(MinUmodes>ntruelayers) {
      nmodes=min(ntruemodes,ntruelayers);
      }
    else {
      nmodes=min(ntruemodes,MinUmodes);
      }
    if(skip) continue;
    status=Umodes_decomposition_v2(ntruelayers, modes, nmodes, p[n], decomposition, debug, 0);
    for(i=0;i<MinUmodes;i++) {
      p_decomposition[i][n]=decomposition[i];
      }
    for(k=0;k<nlayers;k++) {
      size_t m=k*grid.Hsize()+n;
      for(i=0;i<min(nmodes,nkept-1);i++) {
        if(Pmodes[i][m]!=mask) {
          p_modal[i][m]=modes[i][k]*p_decomposition[i][n];
          }
        }
      int ibc=nkept-1;
      p_modal[ibc][m]=0;
      for(i=1;i<nmodes;i++) {
        if(Pmodes[i][m]!=mask) {
          p_modal[ibc][m]+=modes[i][k]*p_decomposition[i][n];
          }
        }
      }
    for(i=0;i<MinUmodes;i++) {
      k=nlayers-1;
      size_t m=k*grid.Hsize()+n;
      if(p_modal[i][m]!=mask) {
        ssh_modal[i][n]=p_modal[i][m]/Rref[0]/9.81;
        }
      else{
        ssh_modal[i][n]=mask;
        }
      }
    }
  
  status=save_SGXYZ("modal-decomposition.nc", grid, P, umask, "p", "dimensionless", "p", create_file, create_grid,"","");
  
  range[0]=0;
  range[1]=MinUmodes-1;
  status=save_SGXYT("modal-decomposition.nc", grid, p_decomposition, umask, "Pc", "dimensionless", "p_coeffcients", create_file, create_grid, range);
  
  range[0]=0;
  range[1]=MinUmodes-1;
  status=save_SGXYT("modal-decomposition.nc", grid, ssh_modal, umask, "ssh_modal", "m", "ssh_modal", create_file, create_grid, range);
  
  for(i=0;i<nkept;i++) {
    create_file=0;
    create_grid=(i==0);
    range[0]=i;
    range[1]=i;
    status=save_SGXYZT("modal-decomposition.nc", grid, p_modal, TYPEmask, "p_modal", "N/m²", "p_baroclinic", create_file, create_grid,"","",range);
    }

  deletep2D(&u,grid.Hsize());
  deletep2D(&v,grid.Hsize());
  deletep2D(&p,grid.Hsize());
  
  delete[] decomposition;
 
  deletep2D(&u_decomposition,MinUmodes);
  deletep2D(&v_decomposition,MinUmodes);
  deletep2D(&p_decomposition,MinUmodes);
  
  double *uflux,*vflux,dmask=NC_FILL_DOUBLE;
  uflux=new double[grid.Hsize()];
  vflux=new double[grid.Hsize()];
  
  for(n=0;n<grid.Hsize();n++) {
    uflux[n]=0;
    vflux[n]=0;
    double weight=0;
    for(k=0;k<nlayers;k++) {
      size_t m=k*grid.Hsize()+n;
      if(p_modal[0][m]!=TYPEmask and u_modal[0][m]!=TYPEmask) {
        uflux[n]+=dz_t[m]*p_modal[0][m]*u_modal[0][m];
        vflux[n]+=dz_t[m]*p_modal[0][m]*v_modal[0][m];
        weight+=dz_t[m];
        }
      }
    if(weight==0) {
      uflux[n]=dmask;
      vflux[n]=dmask;
      }
    }
  create_file=0;
  create_grid=0;
  status=save_SGXY("modal-decomposition.nc", grid, uflux, dmask, "uFbt_x", "N/s", "barotropic_energy_flux", 0);
  status=save_SGXY("modal-decomposition.nc", grid, vflux, dmask, "uFbt_y", "N/s", "barotropic_energy_flux", 0);
      
  int baroclinic=nkept-1;
  
  for(n=0;n<grid.Hsize();n++) {
    uflux[n]=0;
    vflux[n]=0;
    double weight=0;
    for(k=0;k<nlayers;k++) {
      size_t m=k*grid.Hsize()+n;
      if(p_modal[baroclinic][m]!=TYPEmask and u_modal[baroclinic][m]!=TYPEmask) {
        uflux[n]+=dz_t[m]*p_modal[baroclinic][m]*u_modal[baroclinic][m];
        vflux[n]+=dz_t[m]*p_modal[baroclinic][m]*v_modal[baroclinic][m];
        weight+=dz_t[m];
        }
      }
    if(weight==0) {
      uflux[n]=dmask;
      vflux[n]=dmask;
      }
    }
  create_file=0;
  create_grid=0;
  status=save_SGXY("modal-decomposition.nc", grid, uflux, dmask, "uFbc_x", "N/s", "baroclinic_energy_flux", 0);
  status=save_SGXY("modal-decomposition.nc", grid, vflux, dmask, "uFbc_y", "N/s", "baroclinic_energy_flux", 0);
  
  deletep2D(&u_modal,nkept);
  deletep2D(&v_modal,nkept);
  deletep2D(&p_modal,nkept);
  
  poc_var_t vUmodal,vVmodal;

  return(0);
}

