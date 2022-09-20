
/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

***************************************************************************/

#include <stdio.h>
#include <string.h>
#include <cmath>
#include <netcdf.h>
#include <stdlib.h>

#define pi M_PI

#ifdef PARALLEL
#ifdef HAVE_MPI
#include "mpi.h"
#endif
#endif

#include "tugo-prototypes.h"
#include "constants.h"
#include "fe.h"
#include "map.h"
#include "poc-time.h"
#include "netcdf-proto.h"
#include "poc-netcdf-data.hpp"

static int *incidence = NULL;
static char *rootname, *zone;
static grid_t ncgrid;
static double local_regular_sampling;
static int targets[100];

static int h_id,Ha_id,Hg_id,Ua_id,Ug_id,Va_id,Vg_id;
static int LSAa_id,LSAg_id;

#define SEA_FLOOR_DEPTH                 0
#define SEA_SURFACE_ELEVATION           1
#define BAROTROPIC_CURRENTS             2
#define TIDAL_LSA_POTENTIAL             3

extern char *sgetnewdate(date_t reference, double t);
extern tugo_cfg_C *tugo_cfg;
extern double *ctracers[3];
extern   int poc_write_xyt(const char *filename, grid_t grid, int frame, int id, float *buffer);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int archiving_SGtides2D(const mesh_t & mesh,const cefmo_t & cefmo,const tidal_wave & wave,int iteration,const Cgrid_t & Cgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  /*NOTE : static variables are not destroyed when the function returns
  and are not re-initialised when the function is called again,
  a bit like global variables */
  static int nu,nv;
#if 1
  #define T double
  #define NC_T NC_DOUBLE
  #define NC_FILL_T NC_FILL_DOUBLE
#else
  #define T float
  #define NC_T NC_FLOAT
  #define NC_FILL_T NC_FILL_FLOAT
#endif
  static complex<T> *zbuffer=NULL,*ubuffer=NULL,*vbuffer=NULL;
  static const T spec=NC_FILL_T;
  static complex<T> const mask=spec;
  
  int status;
  
  string const waveName=wave.name;
  string const pSNs[2]={"amplitude","phase"};
  
  string const SNS="at_"+waveName+"_frequency_due_to_non_equilibrium_ocean_tide";
  
//  string oSNs[3]={"sea_surface_height_above_geoid","sea_water_x_velocity","sea_water_y_velocity"};
  string oSNs[5]={"sea_surface_height_above_geoid","sea_water_eastward_velocity","sea_water_northward_velocity","astronomical_potential","loading_self-attraction_potential"};

  string const filename=(string) gOutputPath+ "/" + waveName + ".spectral-CGRID.nc";
  
  poc_var_t ZgridVar,UgridVar,VgridVar,FgridVar;
  poc_var_t aV,gV,aVu,gVu,aVv,gVv;
  
  status=poc_save_grid(filename, &ZgridVar, *(Cgrid.z_grid),1,1);
  if(status)__NC_CHKERR_LINE_FILE__(status,"poc_save_grid() error");
  
  status=poc_def_att(filename,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
  
  status=poc_save_grid(filename, &UgridVar, *(Cgrid.u_grid),0,1,"u");
  status=poc_def_att(filename,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
  status=poc_save_grid(filename, &VgridVar, *(Cgrid.v_grid),0,1,"v");

  /* is this really necessary ??? */
  ZgridVar<<poc_att_t(_FillValue,spec);
  UgridVar<<poc_att_t(_FillValue,spec);
  VgridVar<<poc_att_t(_FillValue,spec);

/*------------------------------------------------------------------------------
  elevation */
  aV=ZgridVar;
  aV.init(waveName+"_elevation_a",NC_T,oSNs[0]+"_"+pSNs[0]+"_"+SNS,"m",spec);
  gV=ZgridVar;
  gV.init(waveName+"_elevation_G",NC_T,oSNs[0]+"_"+pSNs[1]+"_"+SNS,"degrees",spec);

  zbuffer=new complex<T> [Cgrid.z_grid->Hsize()];
  
  for(int m=0;m<Cgrid.z_grid->Hsize();m++){
    int n=Cgrid.z_incidence[m];
    if(n<0){
      zbuffer[m]=mask;
      continue;
      }
    zbuffer[m]=cefmo.state.z[n];
    }
  
  status=poc_put_cvara(filename,aV,gV,0,zbuffer);
  
  delete[] zbuffer;
  
/*------------------------------------------------------------------------------
  astronomical potential */
  aV=ZgridVar;
  aV.init(waveName+"_POT_a",NC_T,oSNs[3]+"_"+pSNs[0],"m",spec);
  gV=ZgridVar;
  gV.init(waveName+"_POT_G",NC_T,oSNs[3]+"_"+pSNs[1],"degrees",spec);

  zbuffer=new complex<T> [Cgrid.z_grid->Hsize()];
  
  for(int m=0;m<Cgrid.z_grid->Hsize();m++){
    int n=Cgrid.z_incidence[m];
    if(n<0)
      zbuffer[m]=mask;
    else
      zbuffer[m]=cefmo.state.potential[n];
    }
  
  status=poc_put_cvara(filename,aV,gV,0,zbuffer);
  
  delete[] zbuffer;
  
/*------------------------------------------------------------------------------
  LSA potential */
  aV=ZgridVar;
  aV.init(waveName+"_LSA_a",NC_T,oSNs[4]+"_"+pSNs[0],"m",spec);
  gV=ZgridVar;
  gV.init(waveName+"_LSA_G",NC_T,oSNs[4]+"_"+pSNs[1],"degrees",spec);

  zbuffer=new complex<T> [Cgrid.z_grid->Hsize()];
  
  for(int m=0;m<Cgrid.z_grid->Hsize();m++){
    int n=Cgrid.z_incidence[m];
    if(n<0)
      zbuffer[m]=mask;
    else
      zbuffer[m]=cefmo.state.LSA[n];
    }
  
  status=poc_put_cvara(filename,aV,gV,0,zbuffer);
  
  delete[] zbuffer;
  
/*------------------------------------------------------------------------------
  tidal currents, u nodes */
  aVu=UgridVar;
  aVu.init(waveName+"_u_a_u",NC_T,oSNs[1]+"_"+pSNs[0]+"_"+SNS,"m s-1",spec);
  gVu=UgridVar;
  gVu.init(waveName+"_u_G_u",NC_T,oSNs[1]+"_"+pSNs[1]+"_"+SNS,"degrees",spec);
  aVv=UgridVar;
  aVv.init(waveName+"_v_a_u",NC_T,oSNs[2]+"_"+pSNs[0]+"_"+SNS,"m s-1",spec);
  gVv=UgridVar;
  gVv.init(waveName+"_v_G_u",NC_T,oSNs[2]+"_"+pSNs[1]+"_"+SNS,"degrees",spec);
  
  ubuffer=new complex<T> [Cgrid.u_grid->Hsize()];
  vbuffer=new complex<T> [Cgrid.u_grid->Hsize()];
  
  for(int m=0;m<Cgrid.u_grid->Hsize();m++){
    int n=Cgrid.u_incidence[m];
    if(n<0) {
      ubuffer[m]=mask;
      vbuffer[m]=mask;
      }
    else {
      ubuffer[m]=cefmo.state.u[n];
      vbuffer[m]=cefmo.state.v[n];
      }
    }
 
  status=poc_put_cvara(filename,aVu,gVu,0,ubuffer);
  status=poc_put_cvara(filename,aVv,gVv,0,vbuffer);
  
  delete[] ubuffer;
  delete[] vbuffer;

/*------------------------------------------------------------------------------
  tidal currents, v nodes  */
  aVu=VgridVar;
  aVu.init(waveName+"_u_a_v",NC_T,oSNs[1]+"_"+pSNs[0]+"_"+SNS,"m s-1",spec);
  gVu=VgridVar;
  gVu.init(waveName+"_u_G_v",NC_T,oSNs[1]+"_"+pSNs[1]+"_"+SNS,"degrees",spec);
  aVv=VgridVar;
  aVv.init(waveName+"_v_a_v",NC_T,oSNs[2]+"_"+pSNs[0]+"_"+SNS,"m s-1",spec);
  gVv=VgridVar;
  gVv.init(waveName+"_v_G_v",NC_T,oSNs[2]+"_"+pSNs[1]+"_"+SNS,"degrees",spec);
  
  ubuffer=new complex<T> [Cgrid.v_grid->Hsize()];
  vbuffer=new complex<T> [Cgrid.v_grid->Hsize()];
  
  for(int m=0;m<Cgrid.v_grid->Hsize();m++){
    int n=Cgrid.v_incidence[m];
    if(n<0) {
      ubuffer[m]=mask;
      vbuffer[m]=mask;
      }
    else {
      ubuffer[m]=cefmo.state.u[n];
      vbuffer[m]=cefmo.state.v[n];
      }
    }
 
  status=poc_put_cvara(filename,aVu,gVu,0,ubuffer);
  status=poc_put_cvara(filename,aVv,gVv,0,vbuffer);
  
  delete[] ubuffer;
  delete[] vbuffer;

  return(status);
}
