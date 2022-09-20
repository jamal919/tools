
/*******************************************************************************

  T-UGO tools, 2006-2013

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Converts unstructured to structured tidal atlases

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "version-macros.def" //for VERSION and REVISION

#include "tools-structures.h"

#include "functions.h"
#include "map.h"
#include "xyz.h"
#include "fe.h"
#include "poc-netcdf-data.hpp"
#include "archive.h"
#include "geo.h"
#include "sym-io.h"
#include "tides.h"
#include "cefmo.h"

extern   grid_t get_zgrid(grid_t zgrid);
extern   grid_t get_zgrid(grid_t zgrid);

string cmd;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int gridit_NETCDF(char *varnames[6],const grid_t grid[3], int ngrids, char **wave, int nwave, char *atlas_directory, char *atlas_convention, int iteration, char *rootname, char *paire, int persistence)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  extern toponyms_t paire_map(void);
  struct timeval before;
  
  int i,k, status,id;
  char amplitude[256],phase[256];
  char *input,output[1024];
  tide2D_t state;
  atlas2D_t atlas;
  mesh_t mesh;
  discretisation_t *descriptor;
  int *elts[3]={NULL,NULL,NULL};
  fcomplex *UGbuf[2], mask=1.e+10, *SGbuf[2];
  float *buf=0;

  poc_var_t gvar[3],av,gv;
  
  toponyms_t Pmap=paire_map();
  
  int z_discretisation=UNSET, u_discretisation=UNSET;
  char *euv;
  int ugrid,vgrid;

  status=tide_decode_atlasname(atlas_directory,atlas_convention,wave[0], 0,&input);
  
  printf("#################################################################\n");
  printf("load mesh geometry from %s\n",input);
  if(varnames[0]!=0) {
    status=poc_inq_var(input,varnames[0],&av,1);
    status=poc_get_mesh(input,av,&mesh,1,-1, &z_discretisation);
    if(status) NC_TRAP_ERROR(return,status,1,"poc_get_mesh(\"%s\",(\""+av.name+"\"),..) error",input);
    if(varnames[2]!=0){
      status=poc_inq_var(input,varnames[2],&av,1);
      u_discretisation=poc_get_discretisation(av,&mesh);
      }
    }
  else{
    status=fe_readgeometry(input, &mesh);
    if(status!=0) NC_TRAP_ERROR(return,status,1,"fe_readgeometry(\"%s\",) error",input);
    
    status=paire_discretisation_id(Pmap[paire], &z_discretisation, &u_discretisation);
    if(status!=0) TRAP_ERR_RETURN(status,1,"paire_discretisation_id() error\n");
    }
  
  switch(ngrids) {
    case 1:
      euv=0;
      ugrid=0;
      vgrid=0;
      break;
    case 3:
      euv=strdup("EUV");
      ugrid=1;
      vgrid=2;
      break;
    }

/**----------------------------------------------------------------------------
  detect elements corresponding to grid nodes */
  printf("scanning elements (may take ages): ",input);
  fflush(stdout);
  gettimeofday(&before);
  
  for(i=0;i<ngrids;i++){
    if(euv!=0){
      printf("%c, ",euv[i]); fflush(stdout);
      }
    elts[i]=fe_scan_elements(mesh,grid[i],0,1);
    if(elts[i]==NULL) TRAP_ERR_RETURN(-1,1,"fe_scan_elements() error\n");
    }
  
  printf(" Took %gs.\n",__FUNCTION__,difftime(before));

/**----------------------------------------------------------------------------
  save structured grid(s) */
  for (k=0;k<nwave;k++) {
    sprintf(output,"%s.%s.nc",wave[k],rootname);
    printf("creating file %s\n", output);
/*------------------------------------------------------------------------------
    grids*/
    for(i=0;i<ngrids;i++){
      string loc;
      if(euv!=0)
        loc=string(1,euv[i]);
      int create=(i==0);
      status=poc_save_grid(output,&gvar[i],grid[i],create,1,loc);
      }
    status=poc_def_att(output,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
    status=poc_def_att(output,poc_att_t("history",cmd));
    }

/**----------------------------------------------------------------------------
  first do elevation */
  id=z_discretisation;
  
  printf("#################################################################\n");
  printf("treating elevation, load mesh %s discretisation from %s\n",discretisation_name(z_discretisation),input);
  status=discretisation_init(&mesh,id,0);
  if(status) NC_TRAP_ERROR(return,status,1,"discretisation_init(,%d,0) error",id);
  
  descriptor=get_descriptor_address(mesh,id);
  if(descriptor->nnodes==0) TRAP_ERR_RETURN(-1,1,"get_descriptor() error\n");
  
  UGbuf[0]=new fcomplex [descriptor->nnodes];
  int nxy=grid[0].Hsize();
  SGbuf[0]=new fcomplex[nxy];
  if(varnames[0]==0) {
    sprintf(amplitude,"a_eta_%s",discretisation_name(z_discretisation));
    sprintf(phase,    "G_eta_%s",discretisation_name(z_discretisation));
    }
  else {
    sprintf(amplitude,"%s",varnames[0]);
    if(varnames[1])
      sprintf(phase,"%s",varnames[1]);
    else
      phase[0]='\0';
    }
  if(phase[0]==0)
    buf=new float[max(descriptor->nnodes,nxy)];

  for (k=0;k<nwave;k++) {
    status=tide_decode_atlasname(atlas_directory,atlas_convention,wave[k], 0,&input);
    if(phase[0]!=0)
      status=poc_get_UG3D(input, iteration, amplitude, phase, UGbuf[0]);
    else{
      float inputMask;
      status=poc_inq_var(input, amplitude, &av);
      status=poc_decode_mask(av, 0, 0, &inputMask);
      
      status=poc_get_UG3D(input, iteration, amplitude, buf);
      for(i=0;i<descriptor->nnodes;i++){
        float *bufi=&buf[i];
        if(*bufi==inputMask)
          *bufi=NAN;
        UGbuf[0][i]=*bufi;
        }
      }
    if(status!=0) {
      printf("elevation not done\n");
      continue;
//      return(-1);
      }
    gettimeofday(&before);
    status=fe_map(mesh, UGbuf[0], id, grid[0], elts[0], SGbuf[0], mask);
    if(persistence==1) for(int iteration=0; iteration<1; iteration++) status=map_persistence(grid[0], SGbuf[0], mask, 0);
    STDERR_BASE_LINE_FUNC("fe_map() took %gs\n",difftime(before));
    sprintf(output,"%s.%s.nc",wave[k],rootname);
#if 0
    status=tides_savec1(output, "Ha", "Hg", "elevation", "m", grid[0], SGbuf[0], mask, z_ncgrid);
#else
    av=gv=gvar[0];
    if(phase[0]){
      av.init("Ha",NC_FLOAT,comodo_standard_name("",wave[k]),"m",real(mask));
      gv.init("Hg",NC_FLOAT,comodo_standard_name("",wave[k]),"degrees",real(mask));
      poc_put_cvara(output,av,gv,0,SGbuf[0]);
      }
    else{
      av.init(amplitude,NC_FLOAT,"","",real(mask));
      av<<poc_att_t("production","many attributes not coded yet at " __LINE_FILE_PACKAGE_REVISION);
      for(i=0;i<nxy;i++){
        float *bufi=&buf[i];
        *bufi=real(SGbuf[0][i]);
        if(isnan(*bufi))
          *bufi=real(mask);
        }
      poc_put_vara(output,av,0,buf);
      }
#endif
    }
  
  deletep(&buf);
  delete[] SGbuf[0];
  delete[] UGbuf[0];
  descriptor->destroy();
  
  if(varnames[2]==0) return 0;
  
/**----------------------------------------------------------------------------
  then do velocities */
  id=u_discretisation;
  
  status=tide_decode_atlasname(atlas_directory,atlas_convention,wave[0], 0,&input);
  printf("#################################################################\n");
  printf("treating currents, load mesh %s discretisation from %s\n",discretisation_name(u_discretisation),input);
  status=discretisation_init(&mesh,id,0);
  if(status) NC_TRAP_ERROR(return,status,1,"discretisation_init(,%d,0) error",id);
  
  descriptor=get_descriptor_address(mesh,id);
  if(descriptor->nnodes==0) TRAP_ERR_RETURN(-1,1,"get_descriptor() error\n");
  
  UGbuf[0]=new fcomplex [descriptor->nnodes];
  UGbuf[1]=new fcomplex [descriptor->nnodes];
  
  for (k=0;k<nwave;k++) {
    status=tide_decode_atlasname(atlas_directory,atlas_convention,wave[k], 0,&input);
    if(varnames[2]==0) {
      sprintf(amplitude,"a_u_%s",discretisation_name(u_discretisation));
      sprintf(phase,    "G_u_%s",discretisation_name(u_discretisation));
      }
    else {
      sprintf(amplitude,"%s",varnames[2]);
      sprintf(phase,    "%s",varnames[3]);
      }
    printf("loading u %s %s\n",amplitude, phase);
    status=poc_get_UG3D(input, iteration, amplitude, phase, UGbuf[0]);
    if(status!=0) {
      printf("currents not done\n");
      continue;
      }
    if(varnames[4]==0) {
      sprintf(amplitude,"a_v_%s",discretisation_name(u_discretisation));
      sprintf(phase,    "G_v_%s",discretisation_name(u_discretisation));
      }
    else {
      sprintf(amplitude,"%s",varnames[4]);
      sprintf(phase,    "%s",varnames[5]);
      }
    printf("loading v %s %s\n",amplitude, phase);
    status=poc_get_UG3D(input, iteration, amplitude, phase, UGbuf[1]);
    if(status) NC_TRAP_ERROR(return,status,1,"poc_get_UG3D() error");
    SGbuf[0]=new fcomplex[grid[ugrid].Hsize()];
    SGbuf[1]=new fcomplex[grid[vgrid].Hsize()];
    sprintf(output,"%s.%s.nc",wave[k],rootname);
    status=fe_map(mesh, UGbuf[0], id, grid[ugrid], elts[ugrid], SGbuf[0], mask);
#if 0
    status=tides_savec1(output, "Ua", "Ug", "currents", "m/s", grid[1], SGbuf[0], mask, z_ncgrid);
#else
    av=gv=gvar[ugrid];
    av.init("Ua",NC_FLOAT,comodo_standard_name("",wave[k]),"m",real(mask));
    gv.init("Ug",NC_FLOAT,comodo_standard_name("",wave[k]),"degrees",real(mask));
    status=poc_put_cvara(output,av,gv,0,SGbuf[0]);
#endif
    status=fe_map(mesh, UGbuf[1], id, grid[vgrid], elts[vgrid], SGbuf[1], mask);
#if 0
    status=tides_savec1(output, "Va", "Vg", "currents", "m/s", grid[2], SGbuf[1], mask, z_ncgrid);
#else
    av=gv=gvar[vgrid];
    av.init("Va",NC_FLOAT,comodo_standard_name("",wave[k]),"m",real(mask));
    gv.init("Vg",NC_FLOAT,comodo_standard_name("",wave[k]),"degrees",real(mask));
    status=poc_put_cvara(output,av,gv,0,SGbuf[1]);
#endif
    delete[] SGbuf[0];
    delete[] SGbuf[1];
    }

  delete[] UGbuf[0];
  delete[] UGbuf[1];
  descriptor->destroy();
  
  Pmap.clear();
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gridit_ASCII(char *path, char *meshfile, char *depthfile, grid_t zgrid, char **wave, int nwave, char *atlas_convention, int analysis, char *rootname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double x,y;

  int nndes,status;
  int i,k;
  int *z_elts=NULL,*u_elts=NULL,*v_elts=NULL;
  //float *z_mask,*u_mask,*v_mask;
  char output[1024];
  char hfile[1024],ufile[1024],vname[1024],LSAfile[1024];
  grid_t LSAGRID;
  float    *rbuf=NULL;
  float    *rbufx=NULL,*rbufy=NULL;
  fcomplex *cbuf=NULL;
  float    *rbuffer=NULL,rmask=1.e+10;
  fcomplex *cbuffer=NULL,cmask=complex<float>(rmask,rmask);
  fcomplex *cbufferx=NULL,*cbuffery=NULL;
  mesh_t mesh;

  cdfvar_t variable;
  pocgrd_t z_ncgrid, u_ncgrid, v_ncgrid;

  frame_t frame;

/*----------------------------------------------------------------------------
  Read FE mesh */
  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status !=0) goto error;
    status=fe_list(&mesh);
    if(status !=0) goto error;
    }
  else {
    printf("no mesh file specified; abort...\n");
    goto error;
    }

  status=fe_minmax(mesh, frame);

  nndes=mesh.nvtxs;
//  status=fe_initaffine(&mesh);

  z_elts=fe_scan_elements(mesh,zgrid,0);
  if(z_elts==NULL) goto error;

  u_elts=z_elts;
  v_elts=z_elts;

  zgrid.nz=1;
  zgrid.z=NULL;
  //z_mask=new float[zgrid.nx*zgrid.ny];

  //u_mask=z_mask;
  //v_mask=z_mask;

  for (k=0;k<nwave;k++) {

    sprintf(output,"%s.%s.nc",wave[k],rootname);
    printf("%s\n", output);
    sprintf(LSAfile,"/home/data/loading/FES2004/%s.load.nc",wave[k]);

/*------------------------------------------------------------------------------
    grids*/
    status= poc_createfile(output);

    status=poc_sphericalgrid_xy(output,"",zgrid,&z_ncgrid);

    u_ncgrid=z_ncgrid;
    v_ncgrid=z_ncgrid;

    rbuffer =new float[mesh.nvtxs];
    cbuffer =new fcomplex[mesh.nvtxs];
    cbufferx=new fcomplex[mesh.nvtxs];
    cbuffery=new fcomplex[mesh.nvtxs];

/*------------------------------------------------------------------------------
    depths*/
    if(depthfile!=NULL) {
      status=quoddy_loadr1(depthfile, nndes, rbuffer);
      if(status!=0) goto error;
      rbuf=new float[zgrid.nx*zgrid.ny];
      status=fe_map(mesh,rbuffer,zgrid,z_elts,rbuf,rmask);
      if(status!=0) goto error;
      poc_standardvariable_xy(&variable,"depth",rmask,"m",1., 0.,"model_positive_depth","model positive depth","depth",z_ncgrid);
      status=create_ncvariable(output, &variable);
      status=poc_write_xy(output, zgrid, variable.id,rbuf);
      variable.destroy();
      delete[] rbuf;
      }

/*------------------------------------------------------------------------------
    LSA*/
    sprintf(hfile,"%s/%s.load.s2c",path,wave[k],analysis);
    printf("gridit: treating %s LSA \n",wave[k]);
    status=quoddy_loadc1(hfile, nndes, cbuffer);
    if(status!=0) {
      printf("cannot read %s\n",hfile);
      goto elevation;
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
  /*    if(z_mask[i]==1) {
          printf("problem\n");
          } */
        }
      }
    sprintf(vname,"LSAa");
    poc_standardvariable_xy(&variable,vname,rmask,"m",1., 0.,"tidal_LSA_amplitude","tidal LSA amplitude",vname,z_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output, zgrid, variable.id,rbufx);
    variable.destroy();

    sprintf(vname,"LSAg");
    poc_standardvariable_xy(&variable,vname,rmask,"degree",1., 0.,"tidal_LSA_phase_lag","tidal LSA phase lag",vname,z_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output, zgrid, variable.id,rbufy);
    variable.destroy();

    delete[] cbuf;
    delete[] rbufx;
    delete[] rbufy;

elevation:
/*------------------------------------------------------------------------------
    tidal elevation*/
    if(atlas_convention==0) {
      sprintf(hfile,"%s/%s.ele.%2.2d.s2c",path,wave[k],analysis);
      }
    else {
//       char *root;
//       status=tide_decode_atlasname(path,atlas_convention,wave[k],&root);
//       sprintf(hfile,"%s.s2c",root);
      sprintf(hfile,"%s/%s.ele.%s.s2c",path,wave[k],atlas_convention);
      }
    printf("gridit: treating %s wave elevation from harmonic analysis #%d \n",wave[k],analysis);
    status=quoddy_loadc1(hfile, nndes, cbuffer);
    if(status!=0) {
      printf("cannot read %s\n",hfile);
      goto error;
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
   /*   if(z_mask[i]==1) {
          printf("problem\n");
          } */
        }
      }

//    sprintf(vname,"%s_Ha",wave[k]);
    sprintf(vname,"Ha");
    poc_standardvariable_xy(&variable,vname,rmask,"m",1., 0.,"tidal_elevation_amplitude","tidal elevation amplitude",vname,z_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output, zgrid, variable.id,rbufx);
    variable.destroy();

//    sprintf(vname,"%s_Hg",wave[k]);
    sprintf(vname,"Hg");
    poc_standardvariable_xy(&variable,vname,rmask,"degree",1., 0.,"tidal_elevation_phase_lag","tidal elevation phase lag",vname,z_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output, zgrid, variable.id,rbufy);
    variable.destroy();

    for (i=0;i<zgrid.nx*zgrid.ny;i++) {
      rbufx[i]=z_elts[i];
      }

    sprintf(vname,"triangle");
    poc_standardvariable_xy(&variable,vname,rmask,"none",1., 0.,"element","element",vname,z_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output, zgrid, variable.id,rbufx);
    variable.destroy();

    delete[] cbuf;
    delete[] rbufx;
    delete[] rbufy;

/*------------------------------------------------------------------------------
    tidal currents*/
    if(atlas_convention==0) {
      sprintf(ufile,"%s/%s.uv.%2.2d.v2c",path,wave[k],analysis);
      }
    else {
      sprintf(ufile,"%s/%s.uv.%s.v2c",path,wave[k],atlas_convention);
      }
    printf("gridit: treating %s wave currents from harmonic analysis #%d \n",wave[k],analysis);
    status=quoddy_loadc2(ufile, nndes, cbufferx, cbuffery);
    if(status!=0) {
      printf("cannot read %s\n",ufile);
      goto error;
      }

    cbuf=new fcomplex[zgrid.nx*zgrid.ny];
    status=fe_map(mesh,cbufferx,zgrid,u_elts,cbuf,cmask);
    if(status!=0) goto error;

    rbufx =new float[zgrid.nx*zgrid.ny];
    rbufy =new float[zgrid.nx*zgrid.ny];
    for (i=0;i<zgrid.nx*zgrid.ny;i++) {
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
   /*   if(u_mask[i]==1) {
        printf("problem\n");
         } */
        }
      }

//    sprintf(vname,"%s_Ua",wave[k]);
    sprintf(vname,"Ua");
    poc_standardvariable_xy(&variable,vname,rmask,"m/s",1., 0.,"tidal_eastward_current_amplitude","tidal eastward current amplitude",vname,u_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output,  zgrid, variable.id,rbufx);
    variable.destroy();


//    sprintf(vname,"%s_Ug",wave[k]);
    sprintf(vname,"Ug");
    poc_standardvariable_xy(&variable,vname,rmask,"degree",1., 0.,"tidal_eastward_current_phase_lag","tidal eastward current phase lag",vname,u_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output,  zgrid, variable.id,rbufy);
    variable.destroy();

    delete[] cbuf;
    delete[] rbufx;
    delete[] rbufy;

    cbuf=new fcomplex[zgrid.nx*zgrid.ny];
    status=fe_map(mesh,cbuffery,zgrid,v_elts,cbuf,cmask);
    if(status!=0) goto error;

    rbufx =new float[zgrid.nx*zgrid.ny];
    rbufy =new float[zgrid.nx*zgrid.ny];
    for (i=0;i<zgrid.nx*zgrid.ny;i++) {
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
    /*  if(v_mask[i]==1) {
          printf("problem\n");
          }*/
        }
      }

//    sprintf(vname,"%s_Va",wave[k]);
    sprintf(vname,"Va");
    poc_standardvariable_xy(&variable,vname,rmask,"m/s",1., 0.,"tidal_northward_current_amplitude","tidal northward current amplitude",vname,v_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output,  zgrid, variable.id,rbufx);
    variable.destroy();

//    sprintf(vname,"%s_Vg",wave[k]);
    sprintf(vname,"Vg");
    poc_standardvariable_xy(&variable,vname,rmask,"degree",1., 0.,"tidal_northward_current_phase_lag","tidal nortwardh current phase lag",vname,v_ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output,  zgrid, variable.id,rbufy);
    variable.destroy();

    delete[] cbuf;
    delete[] rbufx;
    delete[] rbufy;

    }

  delete[] z_elts;
//   delete[] u_elts;
//   delete[] v_elts;

  delete[] rbuffer;
  delete[] cbuffer;
  delete[] cbufferx;
  delete[] cbuffery;

  return(0);
error:
  printf("gridit: error detected, quit ... \n");
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gridit_energy_ASCII(char *path, char *meshfile, char *depthfile, grid_t zgrid, char **wave, int nwave, char *atlas_convention, int analysis, char *rootname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nndes,status;
  int k;
  int *z_elts=NULL;
  char output[1024];
  char hfile[1024];
  float    *rbuf=NULL;
  float    *rbuffer=NULL,rmask=1.e+37;
  mesh_t mesh;
  char varname[64],content[64];

  frame_t frame;

/*----------------------------------------------------------------------------
  Read FE mesh */
  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status !=0) goto error;
    status=fe_list(&mesh);
    if(status !=0) goto error;
    }
  else {
    printf("no mesh file specified; abort...\n");
    goto error;
    }

  status=fe_minmax(mesh, frame);

  nndes=mesh.nvtxs;
  //fe_initaffine(mesh);

  z_elts=fe_scan_elements(mesh,zgrid,0);
  if(z_elts==NULL) goto error;

  zgrid.nz=1;
  zgrid.z=NULL;

  for (k=0;k<nwave;k++) {
//    sprintf(output,"%s.%s.nc",wave[k],rootname);
    sprintf(output,"%s.nc",rootname);
    printf("%s\n", output);

/*------------------------------------------------------------------------------
    grids*/
    rbuffer =new float[mesh.nvtxs];

/*------------------------------------------------------------------------------
    tidal elevation*/
    if(atlas_convention==0) {
      sprintf(hfile,"%s/%s.row-wd.%2.2d.LGP1.s2r",path,wave[k],analysis);
      }
    else {
      sprintf(hfile,"%s/%s.row-wd.%2.2d.LGP1.s2r",path,wave[k],analysis);
//      sprintf(hfile,"%s/%s.ele.%s.s2c",path,wave[k],atlas_convention);
      }
//    status=tide_decode_atlasname(path,atlas_convention,wave[k],&root);
    status=quoddy_loadr1(hfile, nndes, rbuffer,content);
    if(status!=0) goto error;
    
    rbuf=new float[zgrid.nx*zgrid.ny];
    status=fe_map(mesh,rbuffer,zgrid,z_elts,rbuf,rmask);
    if(status!=0) goto error;
    
    sprintf(varname,"WD_RoW_%s",wave[k]);
    content[strlen(content)-1]=0;
    status=save_SG(output, zgrid, rbuf, rmask, varname, "w/m²", content, -1);
    delete[] rbuf;
    }

  delete[] z_elts;

  delete[] rbuffer;

  return(0);
error:
  printf("gridit: error detected, quit ... \n");
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

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
    " %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s input [OPTIONS] wave1 [wave2 ...]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Converts unstructured to structured tidal atlases.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  show this help and exit\n"
    );
  print_zone_arg_help();
  printf(
    "  -discretisation  followed by discretisation pair of netcdf input atlas. Used for default variable names. Only used if variable names are not given.\n"
    "  -c  followed by netcdf input atlas file name convention. See below.\n"
    "  -p  followed by directory containing the atlas files. Default: .\n"
    "  -v  followed by an empty-terminated list of up to 6 variable names for the amplitudes and phases of e, u and v. Defaults computed from discretisation as TUGO would.\n"
    "  -a  followed by frame index. Default: -1 (last).\n"
    "  -g  followed by grid file AND by an empty-terminated list of up to 3 gridded variable names for e, u and v\n"
    "  -n  followed by notebook file\n"
    "  -m  followed by mesh file\n"
    "  -z  followed by zone name\n"
    "  -o  followed by some output root name. Default: generic.\n"
    "  -f  followed by input format. Default: ascii. If anything else: netcdf.\n"
    "\n"
    "EXAMPLES :\n"
    "  gridit-cdf -z global -f nc -c WAVE.spectral.nc -o FES2012 M2 -discretisation DNP1xLGP2\n"
    "  gridit-cdf -dx 1m -dy 1m -f nc -c WAVE.spectral.energy.nc -v WaveDragRow '' -o energy-RG M2\n"
    "\n"
    "SEE ALSO: tides-converter\n"
    );
  
  print_tide_decode_atlasname_help();
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int k,n;

  char *keyword=NULL,*notebook=NULL;
  char *meshfile=NULL, *depthfile=NULL,*rootname=NULL,*path=NULL, *gridfile=NULL, *griddedvar[3]={NULL,NULL,NULL}, *zone=NULL;
  char *atlas_convention=0,*discretisation=0,*format=0;
  char *wave[1024];
  grid_t grid[3];
  grid_t cartesian_zgrid;
  mesh_t mesh;
  int nwave=0,analysis=-1;
  geo_t projection;
  frame_t frame,prescribed;
  double dx=NAN,dy=NAN;
  int persistence=0;
  const int nvars=6;
  char *varnames[nvars]={0,0,0,0,0,0};
  int defaultVars=1;

  cmd=fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(read_zone_arg(keyword,argv[n+1],&prescribed,&dx,&dy)) {
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"--persistence")==0) {
      persistence=1;
      n++;
      continue;
      }
    if(strcmp(keyword,"--help")==0) {
      print_help(argv[0]);
      exit(0);
      }
    switch (keyword[0]) {
      case '-':
        if(strcmp(keyword,"-discretisation")==0) {
          discretisation= strdup(argv[n+1]);
          n++;
          n++;
          break;
          }
        switch (keyword[1]) {
/*------------------------------------------------------------------------------
        naming convention for tidal atlases*/
        case 'c' :
          atlas_convention= strdup(argv[n+1]);
          n++;
          n++;
          break;
          
        case 'g' :
          n++;
          gridfile= strdup(argv[n]);
          n++;
/* ----------------------------------------------------------------------------
          empty-terminated list of up to 3 variable names */
          for(k=0;k<3;k++){
            if(argv[n+k][0]==0){
              k++;
              break;
              }
            griddedvar[k]= strdup(argv[n+k]);
            }
          n+=k;
          break;

        case 'h' :
          print_help(argv[0]);
          exit(0);

        case 'z' :
          zone= strdup(argv[n+1]);
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

        case 'o' :
          rootname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'f' :
          format= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          path= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'v' :
          if(varnames[0]){
            if(defaultVars)
              defaultVars=0;
            else{
              __FILE_LINE__(stdout,"multiple use of -v : the last one overrides the previous.");
              for(k=0;k<nvars;k++){
                deletep((void**)(&varnames[k]),free);
                }
              }
            }
/* ----------------------------------------------------------------------------
          empty-terminated list of up to 6 variable names */
          for(k=0;k<nvars;k++){
            if(argv[n+1+k][0]==0){
              k++;
              break;
              }
            varnames[k]=strdup(argv[n+1+k]);
            }
          n++;
          n+=k;
          break;
        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
          wave[nwave]= strdup(argv[n]);
//        printf("input wave=%s\n",wave[nwave]);
          nwave++;
          n++;
        break;
      }
      free(keyword);
    }

  if(format!=0 && strcmp(format,"ascii")!=0 && varnames[0]==0 && discretisation==0) {
    fprintf(stderr,"*** discretisation missing even though input format is NetCDF and variable name is missing ***\n");
    print_help(argv[0]);
    exit(-1);
    }

  if(nwave==0) {
    fprintf(stderr,"*** empty tidal wave list ***\n");
    print_help(argv[0]);
    exit(-1);
    }

  if(rootname ==NULL) rootname=strdup("generic");
  if(path ==NULL) path=strdup(".");
  if(format==0) format=strdup("ascii");

/*----------------------------------------------------------------------------
  Read FE mesh */
  if(meshfile != NULL) {
    printf("#################################################################\n");
    printf("load ascii mesh from file : %s\n",meshfile);
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status) TRAP_ERR_EXIT(status,"fe_readmesh() error\n");
    status=fe_list(&mesh);
    if(status) TRAP_ERR_EXIT(status,"fe_list() error\n");
    }
  else {
    if(strcmp(format,"ascii")==0) {
      fprintf(stderr,"*** no mesh file specified ***\n");
      print_help(argv[0]);
      exit(-1);
      }
    }


/*----------------------------------------------------------------------------
  Regular grid prescription: priority sequence*/
  if(notebook!=0) {
/*----------------------------------------------------------------------------
    Read notebook data */
    printf("#################################################################\n");
    printf("load grid from notebook file : %s\n",notebook);
    status=load_notebook(notebook, &cartesian_zgrid, &grid[0], &projection);
    printf("%s (notebook file) processed\n",notebook);
    if(status!=0) {
      STDOUT_BASE_LINE("unable to load notebook file: %s\n",notebook);
      exit(-1);
      }
    }
  else if(gridfile != NULL) {
/*----------------------------------------------------------------------------
    Read grid file */
    printf("#################################################################\n");
    printf("load regular grid from %s for",gridfile);
    for(k=0;k<3;k++){
      if(griddedvar[k]==0) break;
      printf(" %s",griddedvar[k]);fflush(stdout);
      status=poc_get_grid(gridfile,griddedvar[k],&grid[k],0);
      if(status) NC_TRAP_ERROR(wexit,status,1,"poc_get_grid(\"%s\",\"%s\",) error",gridfile,griddedvar[k]);
      if(!grid[k].nx || !grid[k].ny) TRAP_ERR_EXIT(status,"empty grid\n");
      switch(k){
        case 0:printf(",");break;
        case 1:printf(" and");break;
        case 2:printf("\n");break;}
      }
    }
  else if(zone != NULL) {
/*----------------------------------------------------------------------------
    Apply zone definition by name*/
    printf("#################################################################\n");
    printf("define grid from zone definition : %s\n",zone);
    grid[0]=get_zonegrid(zone);
/**----------------------------------------------------------------------------
    allow for lighter grid storage */
//    status=map_completegridaxis(&grid[0],2);
    status=map_completegridaxis(&grid[0],1);
    }
  else {
    if(!mesh.vertices){
      status=tide_decode_atlasname(path,atlas_convention,wave[0], 0,&meshfile);
      printf("#################################################################\n");
      printf("load mesh geometry from %s\n",meshfile);
      status=fe_readgeometry(meshfile, &mesh);
      if(status) NC_TRAP_ERROR(return,status,1,"fe_readgeometry(\"%s\",) error",meshfile);
      }
/*----------------------------------------------------------------------------
    Apply frame definition*/
    printf("#################################################################\n");
    printf("define grid from arguments: ");fflush(stdout);
    
    status=fe_minmax(mesh, frame);
    if(status)TRAP_ERR_EXIT(status,"fe_minmax() error\n");
    mesh.destroy();
    
/*----------------------------------------------------------------------------
    default frame (mesh limits)*/
    grid[0].xmin=frame.xmin;
    grid[0].ymin=frame.ymin;
    grid[0].xmax=frame.xmax;
    grid[0].ymax=frame.ymax;
    
/*----------------------------------------------------------------------------
    apply prescribed limits and resolution*/
    apply_zone_arg(&grid[0],prescribed,dx,dy);
    grid[0].modeH=0;
    
    grid[0].brief_print();
    
    status=map_completegridaxis_2(&grid[0]);
    }

  if(strcmp(format,"ascii")==0) {
    printf("#################################################################\n");
    printf("treat ascii files\n");
    
    status=gridit_ASCII(path, meshfile, depthfile, grid[0], wave, nwave, atlas_convention, analysis, rootname);
    if(status) TRAP_ERR_EXIT(status,"gridit_ASCII() error\n");
/*------------------------------------------------------------------------------
    temporary patch for Ariane */
//     status=gridit_energy_ASCII(path, meshfile, depthfile, grid[0], wave, nwave, atlas_convention, analysis, rootname);
//     if(status) TRAP_ERR_EXIT(status,"gridit_energy_ASCII() error\n");
    }
  else {
/*----------------------------------------------------------------------------
    single or staggered grids*/
    int ngrids=1;
    for(k=1;k<3;k++){
      if(grid[k].x==0) continue;
      grid[k]=grid[0];
      ngrids++;
      }
    printf("#################################################################\n");
    printf("treat netcdf files\n");
    status=gridit_NETCDF(varnames, grid, ngrids, wave, nwave, path, atlas_convention, analysis, rootname, discretisation,persistence);
    if(status) TRAP_ERR_EXIT(status,"gridit_NETCDF() error\n");
    }
  
  return 0;
}
