
/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

#define MAIN_SOURCE


#include <stdio.h>
#include <stdlib.h>
#include <string.h>


 
#include "tools-structures.h"

#include "functions.h"
#include "map.h"
#include "xyz.h"
#include "fe.h"
#include "netcdf-proto.h"
#include "archive.h"
#include "geo.h"
#include "sym-io.h"
#include "tides.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double x,y;

  float  spec=-9999;
  int nndes,status;
  int i,k,n;
  int *z_elts,*u_elts,*v_elts;
  float *z_mask,*u_mask,*v_mask;

  char *keyword,*notebook=NULL;
  char *meshfile=NULL, *depthfile=NULL,*output=NULL,*path=NULL;
  char *zgridfile=NULL,*ugridfile=NULL,*vgridfile=NULL;
  char hfile[1024],ufile[1024],vname[1024],LSAfile[1024];
  char *wave[1024];
  grid_t zgrid,ugrid,vgrid,LSAGRID;
  grid_t cartesian_zgrid,cartesian_ugrid,cartesian_vgrid;
  float *rbuf;
  float *rbufx,*rbufy;
  fcomplex *cbuf;
  float  *rbuffer,rmask;
  fcomplex *cbuffer,cmask;
  fcomplex *cbufferx,*cbuffery;
  mesh_t mesh;
  int nwave=0,analysis;
  spectrum_t spectrum;
  string cmd;
  cdfvar_t variable;
  pocgrd_t z_ncgrid, u_ncgrid, v_ncgrid;
  geo_t projection;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
//         case 'g' :
//           gridfile= strdup(argv[n+1]);
//           n++;
//           n++;
//           break;²

        case 'z' :
          zgridfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'u' :
          ugridfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'v' :
          vgridfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'n' :
          notebook= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'd' :
          depthfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'a' :
          sscanf(argv[n+1],"%d",&analysis);
          n++;
          n++;
          break;

        case 's' :
          sscanf(argv[n+1],"%f",&spec);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          path= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
          wave[nwave]= strdup(argv[n]);
//        printf("input wave=%s\n",wave[nwave]);
          spectrum.n=nwave+1;
          nwave++;
          n++;
        break;
      }
      free(keyword);
    }

  spectrum.waves=new tidal_wave[spectrum.n];
  for (k=0;k<nwave;k++) {
    strcpy(spectrum.waves[k].name,wave[k]);
    }

  rmask=spec;
  cmask=fcomplex(spec,spec);

//  if((gridfile == NULL) && (zgridfile == NULL)) goto error;
  
  status=0;
  
  if(output==NULL){
    status=1;
    printf("*** please specify  with option  ***\n");
    }
  
  if(notebook==NULL){
    status=1;
    printf("*** please specify notebook file with option -n ***\n");
    }
  
  if(meshfile==NULL){
    status=1;
    printf("*** please specify mesh file with option -m ***\n");
    }
  
  if(status!=0)
    TRAP_ERR_EXIT(-1,"sorry, no more help\n");
  
  if(path ==NULL) path=strdup(".");

/* *-----------------------------------------------------------------------------
  Read notebook data */
  if(notebook!=0) {
//    status=load_notebook_obsolete(notebook,&cartesian_zgrid,&projection);
    status=load_notebook(notebook, &cartesian_zgrid, &zgrid, &projection);
    if(status!=0) {
      STDOUT_BASE_LINE("unable to load notebook file: %s\n",notebook);
      exit(-1);
      }
    cartesian_ugrid=map_z2ugrid(cartesian_zgrid);
    cartesian_vgrid=map_z2vgrid(cartesian_zgrid);
//    zgrid=get_spherical(projection,cartesian_zgrid);
    ugrid=map_get_spherical(projection,cartesian_ugrid);
    vgrid=map_get_spherical(projection,cartesian_vgrid);
    printf("%s (notebook file) processed\n",notebook);
    }

/* *-----------------------------------------------------------------------------
  Read FE mesh */
 if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status !=0) TRAP_ERR_EXIT(status,"fe_readmesh(\"%s\",) error\n",meshfile);
    status=fe_list(&mesh);
    if(status !=0) TRAP_ERR_EXIT(status,"fe_list() error\n");
    }

  nndes=mesh.nvtxs;
  //fe_initaffine(&mesh);

//   status=nc_writespectrum(ncid, spectrum);
//   if(status !=0) TRAP_ERR_EXIT(status,"nc_writespectrum() error\n");

  z_elts=fe_scan_elements(mesh,zgrid,0);
  if(z_elts==NULL) TRAP_ERR_EXIT(status,"fe_scan_elements() error\n");
  zgrid.nz=1;
  zgrid.z=NULL;
  z_mask=new float[zgrid.nx*zgrid.ny];

  u_elts=fe_scan_elements(mesh,ugrid,0);
  if(u_elts==NULL) TRAP_ERR_EXIT(status,"fe_scan_elements() error\n");
  ugrid.nz=1;
  ugrid.z=NULL;
  u_mask=new float[ugrid.nx*ugrid.ny];

  v_elts=fe_scan_elements(mesh,vgrid,0);
  if(v_elts==NULL) TRAP_ERR_EXIT(status,"fe_scan_elements() error\n");
  vgrid.nz=1;
  vgrid.z=NULL;
  v_mask=new float[vgrid.nx*vgrid.ny];

  atlas_grid_or_mesh_t gm;
  
  for (k=0;k<nwave;k++) {

    sprintf(output,"%s.nc",wave[k]);
    sprintf(LSAfile,"/home/softs/data/loading/FES2004/%s.load.nc",wave[k]);

/* *------------------------------------------------------------------------
    grids*/
    status= poc_createfile(output);

    status=poc_sphericalgrid_xy(output,"ZHL",zgrid,&z_ncgrid);
    status=poc_sphericalgrid_xy(output,"XHL",ugrid,&u_ncgrid);
    status=poc_sphericalgrid_xy(output,"YHL",vgrid,&v_ncgrid);

    rbuffer =new float[mesh.nvtxs];
    cbuffer =new fcomplex[mesh.nvtxs];
    cbufferx=new fcomplex[mesh.nvtxs];
    cbuffery=new fcomplex[mesh.nvtxs];

/* *------------------------------------------------------------------------
    depths*/
    if(depthfile!=NULL) {
      status=quoddy_loadr1(depthfile, nndes, rbuffer);
      if(status!=0) TRAP_ERR_EXIT(status,"quoddy_loadr1(\"%s\") error\n",depthfile);
      rbuf=new float[zgrid.nx*zgrid.ny];
      status=fe_map(mesh,rbuffer,zgrid,z_elts,rbuf,rmask);
      if(status!=0) TRAP_ERR_EXIT(status,"fe_map() error\n");
      poc_standardvariable_xy(&variable,"depth",rmask,"m",1., 0.,"model_positive_depth","model positive depth","depth",z_ncgrid);
      status=create_ncvariable(output, &variable);
      status=poc_write_xy(output, zgrid, variable.id,rbuf);
        variable.destroy();
      delete[] rbuf;
      }

/* *------------------------------------------------------------------------
    LSA*/
    rbufx =new float[zgrid.nx*zgrid.ny];
    rbufy =new float[zgrid.nx*zgrid.ny];
    status=tide_atlas2positions(LSAfile,"Ha","Hg",zgrid.x,zgrid.y,zgrid.nx*zgrid.ny,rbufx,rbufy,rmask,&gm);
//    sprintf(vname,"%s_Ha",wave[k]);
    sprintf(vname,"LSAa");
    poc_standardvariable_xy(&variable,vname,rmask,"m",1., 0.,"tidal_LSA_amplitude","tidal LSA amplitude",vname,z_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output, zgrid, variable.id,rbufx);
    variable.destroy();

//    sprintf(vname,"%s_Hg",wave[k]);
    sprintf(vname,"LSAg");
    poc_standardvariable_xy(&variable,vname,rmask,"degree",1., 0.,"tidal_LSA_phase_lag","tidal LSA phase lag",vname,z_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output, zgrid, variable.id,rbufy);
    variable.destroy();

    delete[] rbufx;
    delete[] rbufy;

/* *------------------------------------------------------------------------
    tidal elevation*/
    sprintf(hfile,"%s/%s.ele.%2.2d.s2c",path,wave[k],analysis);
    printf("gridit: treating %s wave from harmonic analysis #%d \n",wave[k],analysis);
    status=quoddy_loadc1(hfile, nndes, cbuffer);
    if(status!=0) {
      TRAP_ERR_EXIT(status,"cannot read %s (%d %s)\n",hfile,status,strerror(status));
      }

    cbuf=new fcomplex[zgrid.nx*zgrid.ny];
    status=fe_map(mesh,cbuffer,zgrid,z_elts,cbuf,cmask);

    rbufx =new float[zgrid.nx*zgrid.ny];
    rbufy =new float[zgrid.nx*zgrid.ny];
    for (i=0;i<zgrid.nx*zgrid.ny;i++) {
      x=real(cbuf[i]);
      y=imag(cbuf[i]);
      if(x!=rmask) {
/*------------------------------------------------------------------------------
        convert fcomplex h in amplitude (meters) and phase lag (degrees)*/
        rbufx[i]=sqrt(x*x+y*y);
        rbufy[i]=atan2(-y,x)*r2d;
        }
      else {
        rbufx[i]=rmask;
        rbufy[i]=rmask;
        if(z_mask[i]==1) {
          printf("problem\n");
          }
        }
      }

//    sprintf(vname,"%s_Ha",wave[k]);
    sprintf(vname,"Ha");
    poc_standardvariable_xy(&variable,vname,rmask,"m",1., 0.,"tidal_elevation_amplitude","tidal elevation amplitude",vname,z_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output, zgrid, variable.id,rbufx);
//    status=free_ncvariable(&variable);
    variable.destroy();

//    sprintf(vname,"%s_Hg",wave[k]);
    sprintf(vname,"Hg");
    poc_standardvariable_xy(&variable,vname,rmask,"degree",1., 0.,"tidal_elevation_phase_lag","tidal elevation phase lag",vname,z_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output, zgrid, variable.id,rbufy);
//    status=free_ncvariable(&variable);
    variable.destroy();

    delete[] cbuf;
    delete[] rbufx;
    delete[] rbufy;

/* *------------------------------------------------------------------------
    tidal currents*/
    sprintf(ufile,"%s/%s.uv.%2.2d.v2c",path,wave[k],analysis);
    status=quoddy_loadc2(ufile, nndes, cbufferx, cbuffery);
    if(status!=0) {
      TRAP_ERR_EXIT(status,"cannot read %s (%d %s)\n",ufile,status,strerror(status));
      }

    cbuf=new fcomplex[ugrid.nx*ugrid.ny];
    status=fe_map(mesh,cbufferx,ugrid,u_elts,cbuf,cmask);
    if(status!=0) TRAP_ERR_EXIT(status,"fe_map() error\n");

    rbufx =new float[ugrid.nx*ugrid.ny];
    rbufy =new float[ugrid.nx*ugrid.ny];
    for (i=0;i<ugrid.nx*ugrid.ny;i++) {
      x=real(cbuf[i]);
      y=imag(cbuf[i]);
      if(x!=rmask) {
/*------------------------------------------------------------------------------
        convert fcomplex u in amplitude (meters) and phase lag (degrees) */
        rbufx[i]=sqrt(x*x+y*y);
        rbufy[i]=atan2(-y,x)*r2d;
        }
      else {
        rbufx[i]=rmask;
        rbufy[i]=rmask;
        if(u_mask[i]==1) {
          printf("problem\n");
          }
        }
      }

//    sprintf(vname,"%s_Ua",wave[k]);
    sprintf(vname,"Ua");
    poc_standardvariable_xy(&variable,vname,rmask,"m/s",1., 0.,"tidal_eastward_current_amplitude","tidal eastward current amplitude",vname,u_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output,  ugrid, variable.id,rbufx);
//    status=free_ncvariable(&variable);
    variable.destroy();

//    sprintf(vname,"%s_Ug",wave[k]);
    sprintf(vname,"Ug");
    poc_standardvariable_xy(&variable,vname,rmask,"degree",1., 0.,"tidal_eastward_current_phase_lag","tidal eastward current phase lag",vname,u_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output,  ugrid, variable.id,rbufy);
//    status=free_ncvariable(&variable);
    variable.destroy();

    delete[] cbuf;
    delete[] rbufx;
    delete[] rbufy;

    cbuf=new fcomplex[vgrid.nx*vgrid.ny];
    status=fe_map(mesh,cbuffery,vgrid,v_elts,cbuf,cmask);
    if(status!=0) TRAP_ERR_EXIT(status,"fe_map() error\n");

    rbufx =new float[vgrid.nx*vgrid.ny];
    rbufy =new float[vgrid.nx*vgrid.ny];
    for (i=0;i<vgrid.nx*vgrid.ny;i++) {
      x=real(cbuf[i]);
      y=imag(cbuf[i]);
      if(x!=rmask) {
/*------------------------------------------------------------------------------
        convert fcomplex v in amplitude (meters) and phase lag (degrees) */
        rbufx[i]=sqrt(x*x+y*y);
        rbufy[i]=atan2(-y,x)*r2d;
        }
      else {
        rbufx[i]=rmask;
        rbufy[i]=rmask;
        if(v_mask[i]==1) {
          printf("problem\n");
          }
        }
      }

//    sprintf(vname,"%s_Va",wave[k]);
    sprintf(vname,"Va");
    poc_standardvariable_xy(&variable,vname,rmask,"m/s",1., 0.,"tidal_northward_current_amplitude","tidal north current amplitude",vname,v_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output,  vgrid, variable.id,rbufx);
//    status=free_ncvariable(&variable);
    variable.destroy();

//    sprintf(vname,"%s_Vg",wave[k]);
    sprintf(vname,"Vg");
    poc_standardvariable_xy(&variable,vname,rmask,"degree",1., 0.,"tidal_northward_current_phase_lag","tidal north current phase lag",vname,v_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output,  vgrid, variable.id,rbufy);
//    status=free_ncvariable(&variable);
    variable.destroy();

    delete[] cbuf;
    delete[] rbufx;
    delete[] rbufy;

    }

  delete[] z_elts;
  delete[] u_elts;
  delete[] v_elts;

  delete[] rbuffer;
  delete[] cbuffer;
  delete[] cbufferx;
  delete[] cbuffery;

  TRAP_ERR_EXIT(0,"exiting\n");
}
