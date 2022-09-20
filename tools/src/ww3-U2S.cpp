
/***************************************************************************
T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative
  
E-mail: florent.lyard@legos.obs-mip.fr

\file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada
\author  Clement MAYET      LEGOS, Toulouse, France (PhD)

VERSION :
\date Jul 20, 2011 : Clement MAYET : add the possibility to write a NC_SHORT variable
\date Jul 20, 2011 : Clement MAYET : bugfix : apply variable scale factor from netcdf attribute
\date  21/06/2011 : Clement MAYET : first working version
\date  10/06/2011 : Florent LYARD : first pieces

\brief WaveWatch III Unstructured netcdf output converter (unstructured to structured grid)

\todo 21/06/2011 : appropriate output file name, argument reading, help function, possibility to interpolate on any output grid

***************************************************************************/


#define MAIN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"

#include "tools-define.h"
#include "tools-structures.h"

#include "poc-netcdf.def"

#include "fe.h"
#include "map.h"
#include "geo.h"
#include "sym-io.h"
#include "polygones.h"
#include "grd.h"
#include "netcdf-proto.h"
#include "poc-time.h"

#define T_VAR     0
#define S_VAR     1
#define RHO_VAR   2
#define U_VAR     3
#define V_VAR     4
#define W_VAR     5
#define UBAR_VAR  6
#define VBAR_VAR  7


class index_t{
  private:
  public:
  int i,j,k;

  index_t() {
    i=j=k=-1;
  }
};




/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int colocate(mesh_t mesh, grid_t grid, index_t *indices)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k(-1),l(-1),n(0),status(-1);
  double x,y;
  
  for(n=0;n<mesh.nvtxs;n++) {
    x = mesh.vertices[n].lon;
    y = mesh.vertices[n].lat;
    status=map_index(grid,x,y,&k,&l);
    if(status==0){
      indices[n].i=k;
      indices[n].j=l;
      }
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int duplicate_data(const char *input, const char *output, cdfvar_t variable, T proto)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int ncid,status;
  
  T *z;

  z=new T[variable.size()];
  status=nc_open(input,NC_NOWRITE,&ncid);
  if(status != NC_NOERR)
    goto error;

  status=poc_get_var(ncid, variable.id,z);

  status = nc_close(ncid);
  if(status != NC_NOERR)
    goto error;
  
  status=create_ncvariable(output, &variable);

  status=nc_open(output,NC_WRITE,&ncid);
  if(status != NC_NOERR)
    goto error;

  status=poc_put_var(ncid, variable.id, z);

  status = nc_close(ncid);
  if(status != NC_NOERR)
    goto error;
  
  return(status);
error:
  return(-1);
}

/**
\brief Regrid from .... to .....

\todo 21/06/2011 : Clement MAYET : MEMORY ALLOCATION of zSM, add the possibility to interpolate on any output grid
@param
@return
*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int regrid(const char *input, grid_t meshgrid, grid_t grid, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int f,time,m,n,varNum,status,verbose=0;
  int ncid;
  cdfgbl_t g_info;
  cdfvar_t vU_info,vS_info;
  mesh_t mesh;
  float *zU(0);		// variable from Unstructured grid
  float *zSM(0);	// variable on Structured Vorticity grid
  float *zSV(0);	// variable on Structured Mass grid
  index_t *indices(0);
  grid_t vgrid;
  UGdecoded_t UGdecoded;
  decoded_t decoded;
  string output;
  pocgrd_t ncgrid;
  double factor;


  output=input;
  output="SG_"+output;
  
/*----------------------------------------------------------------------------
  get unstructured netcdf file details*/
  status=cdf_globalinfo(input,&g_info,verbose);
  
/*----------------------------------------------------------------------------
  read mesh */

  status=fe_readmesh3d(input,&mesh,verbose);
  
  switch(target) {
    case 0:
      vgrid = meshgrid;
      break;
    case 1:
/*----------------------------------------------------------------------------
  create vorticity grid from mass grid */
      vgrid = map_vgrid(meshgrid);
      break;
    }
    
  indices = new index_t[mesh.nvtxs];
  zU	  = new float[mesh.nvtxs];
  zSV	  = new float[vgrid.nx*vgrid.ny];
  
  
  status  = colocate(mesh, vgrid, indices);
  
/*----------------------------------------------------------------------------
  create netcdf outpput file */
  status = poc_createfile(output.c_str());
  
/*----------------------------------------------------------------------------
  add grid */
  int fdim=cdf_identify_dimension(g_info,"f");
  if(fdim==-1) {
    status=poc_sphericalgrid_xyzt(output.c_str(),"","",vgrid,&ncgrid);
    }
  else {
    status=poc_sphericalgrid_xyzwt(output.c_str(),"",(int) g_info.dimension[fdim].length, vgrid,&ncgrid);
    }

  status=nc_open(input,0, &ncid);

  for(varNum=0;varNum<g_info.nvarsp;varNum++) {
/*----------------------------------------------------------------------------
    get unstructured netcdf variable details*/
    status=cdf_varinfo(input, varNum, &vU_info);
    status=poc_UGdecode_axis(vU_info,g_info,&UGdecoded);
    status=poc_decode_mask(vU_info,&decoded);
    status=poc_decode_names(vU_info,&decoded);
    status=poc_decode_units(vU_info,&decoded,&factor);
    if ( strcmp(vU_info.name,"f")== 0 ) {
      status=duplicate_data(input, output.c_str(), vU_info, (float) 1.);
      continue;
      }
    if ( strcmp(vU_info.name,"time")== 0 ) {
      status=duplicate_data(input, output.c_str(), vU_info, (double) 1.);
      continue;
      }
    // work only on the non-grid variables
    if ( (strcmp(vU_info.name,"longitude")== 0) || (strcmp(vU_info.name,"latitude")== 0) ) {
      continue;
      }
    if ( (strcmp(vU_info.name,"tri")== 0) || (strcmp(vU_info.name,"MAPSTA")== 0) ) {
      continue;
      }
    if(strcmp(UGdecoded.axis,"TN")==0) {
/*----------------------------------------------------------------------------
      create structured netcdf variable */
      switch (vU_info.type) {
        case NC_SHORT:
          status=poc_standardvariable_xyt(vU_info.name, short(decoded.spec) ,decoded.units,decoded.scale, decoded.offset,decoded.standard_name,decoded.long_name,vU_info.name,ncgrid,vS_info);
          break;
        case NC_FLOAT:
          status=poc_standardvariable_xyt(vU_info.name,decoded.spec,decoded.units,decoded.scale, decoded.offset,decoded.standard_name,decoded.long_name,vU_info.name,ncgrid,vS_info);
          break;
        default:
          status=poc_standardvariable_xyt(vU_info.name,decoded.spec,decoded.units,decoded.scale, decoded.offset,decoded.standard_name,decoded.long_name,vU_info.name,ncgrid,vS_info);
          break;
        }
      status=create_ncvariable(output.c_str(), &vS_info);
/*----------------------------------------------------------------------------
      put UG solution on SG grid */
      for(time=0;time<UGdecoded.tlen;time++) {
        for (int i=0;i<vgrid.nx*vgrid.ny;i++) zSV[i]= decoded.spec; // initialization
        status=poc_get_UG3D(ncid,time,varNum,zU);
        for(n=0;n<mesh.nvtxs;n++) {
          m=indices[n].j*(vgrid.nx)+indices[n].i;
          zSV[m]=zU[n];
          }
        status=poc_write_xyt(output.c_str(),  vgrid, time, vS_info.id,zSV);
        }
      if(status==0) continue;
      }
    else if(strcmp(UGdecoded.axis,"TFN")==0) {
      switch (vU_info.type) {
        case NC_SHORT:
          vS_info=poc_standardvariable_xywt(vU_info.name, short(decoded.spec) ,decoded.units,decoded.scale, decoded.offset,decoded.standard_name,decoded.long_name,vU_info.name,ncgrid);
          break;
        case NC_FLOAT:
          vS_info=poc_standardvariable_xywt(vU_info.name,decoded.spec,decoded.units,decoded.scale, decoded.offset,decoded.standard_name,decoded.long_name,vU_info.name,ncgrid);
          break;
        default:
//          vS_info=poc_standardvariable_xyt(vU_info.name,decoded.spec,decoded.units,decoded.scale, decoded.offset,decoded.standard_name,decoded.long_name,vU_info.name,ncgrid);
          break;
        }
      status=create_ncvariable(output.c_str(), &vS_info);
/*----------------------------------------------------------------------------
      put UG solution on SG grid */
      for(time=0;time<UGdecoded.tlen;time++) {
        for(f=0;f<UGdecoded.flen;f++) {
          for (int i=0;i<vgrid.nx*vgrid.ny;i++) zSV[i]= decoded.spec; // initialization
          status=poc_getLast(ncid,time,f,varNum,zU);
          for(n=0;n<mesh.nvtxs;n++) {
            m=indices[n].j*(vgrid.nx)+indices[n].i;
            zSV[m]=zU[n];
            }
          status=poc_putLast(output.c_str(), time, f, vS_info.id, zSV);
          }
        if(status==0) continue;
        }
      }
    }
    
// cdfvar_t poc_standardvariable_xyw(char *name, float mask,char *units,float scale,float offset, char *standardname,char *longname,char *shortname,pocgrd_t grid)

  delete indices;
  indices=0;
  delete zU;
  zU=0;
  delete zSV;
  zSV=0;
  
  status=nc_close(ncid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  geo_t projection;
  int n,status;
  int target;
  char *keyword;
  grid_t cgrid,sgrid;
  grid_t cmeshgrid,smeshgrid;
  char *rootname=NULL,*output=NULL,*input=NULL,*discretisation=NULL,*bathymetry=NULL;
  char *notebook=NULL,*meshbook=NULL;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'd' :
          discretisation= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'r' :
          rootname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          meshbook= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'n' :
          notebook= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
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
          input= strdup(argv[n]);
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

  if(rootname==NULL) rootname= strdup("quaker");

  if(discretisation==0) {
    printf("target grid not passed in arguments, use T-GRID as a default\n");
    discretisation=strdup("T-GRID");
    }
  if(strcmp(discretisation,"T-GRID")==0) {
    target=0;
    }
  else if(strcmp(discretisation,"F-GRID")==0) {
    target=1;
    }
    
/*-----------------------------------------------------------------------------
  build structured grid from notebook */
  if(notebook!=NULL) {
    status=load_notebook(notebook, &cgrid, &sgrid, &projection);
    printf("%s (mesh book file) processed\n",notebook);
    }
/*-----------------------------------------------------------------------------
  build structured grid from notebook */
  if(meshbook!=NULL) {
    status=load_notebook(meshbook, &cmeshgrid, &smeshgrid, &projection);
    printf("%s (mesh book file) processed\n",meshbook);
    }
  else {
    smeshgrid=sgrid;
    cmeshgrid=cgrid;
    }
  status=regrid(rootname, smeshgrid, sgrid, target);  // spherical sgrid=spherical grid  (lon/lat) ( cgrid = cartesian grid (km) )
    
  __ERR_BASE_LINE__("exiting\n");exit(0);

error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);

}
/*
      variable[0] =poc_variable_UG4D("u", mask, "m/s", scale, offset,"eastward_velocity","T","E","K");
      variable[0] =poc_variable_UG4D("v", mask, "m/s", scale, offset,"northward_velocity","T","E","K");
      variable[0] =poc_variable_UG4D("w", mask, "m/s", scale, offset,"vertical_velocity","T","M","L");
      variable[0] =poc_variable_UG4D("U", mask, "kg/m^2s", scale, offset,"eastward_momentum","T","E","K");
      variable[0] =poc_variable_UG4D("V", mask, "kg/m^2s", scale, offset,"northward_momentum","T","E","K");
      variable[0] =poc_variable_UG4D("W", mask, "kg/m^2s", scale, offset,"vertical_momentum","T","M","L");
      variable[0] =poc_variable_UG4D("T", mask, "C", scale, offset,"potential_sea_water_temperature","T","M","K");
      variable[0] =poc_variable_UG4D("S", mask, "PSU", scale, offset,"sea_water_sanility","T","M","K");
      variable[0] =poc_variable_UG3D("sfd", dmask, "m", dscale, doffset,"sea_floor_deformation");
      variable[0] =poc_variable_UG3D("elevation", dmask, "m", dscale, doffset,"sea_surface_elevation");
      variable[0] =poc_floatvariable_nt("ubar", mask, "m/s", scale, offset,"mean_eastward_velocity");
      variable[0] =poc_floatvariable_nt("vbar", mask, "m/s", scale, offset,"mean_northward_velocity");
      variable[0] =poc_variable_UG4D("u", mask, "m/s", scale, offset,"eastward_velocity","T","E","K");
      variable[0] =poc_variable_UG4D("v", mask, "m/s", scale, offset,"northward_velocity","T","E","K");
      variable[0] =poc_variable_UG4D("w", mask, "m/s", scale, offset,"vertical_velocity","T","M","L");
      variable[0] =poc_variable_UG4D("U", mask, "kg/m^2s", scale, offset,"eastward_momentum","T","E","K");
      variable[0] =poc_variable_UG4D("V", mask, "kg/m^2s", scale, offset,"northward_momentum","T","E","K");
      variable[0] =poc_variable_UG4D("W", mask, "kg/m^2s", scale, offset,"vertical_momentum","T","M","L");
      variable[0] =poc_variable_UG4D("T", mask, "C", scale, offset,"potential_sea_water_temperature","T","M","K");
      variable[0] =poc_variable_UG4D("S", mask, "PSU", scale, offset,"sea_water_sanility","T","M","K");
*/
