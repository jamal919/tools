
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

#include "tugo-prototypes.h"
#include "constants.h"
#include "keywords.hpp"
#include "drag.hpp"
#include "map.h"
#include "grd.h"
#include "geo.h"

#include "polygons.hpp"

extern float tide_prediction(double, spectrum_t, hconstant_t, int);
extern  int fe_ascii_loadr1(const char *, mesh_t, float *, float *, int);



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int initializeZ0_ByRegions (const char *ValuesFilename, const mesh_t &mesh, float *z0, double defaultZ0) 

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  int nitems,flag;
  FILE *in;
  double t,p,value;
  char toponym[64];
  grid_t grid;
  float *code, mask, *region, z;
  size_t count;
  
  const char *filename=strdup(tugo_cfg->misc->ToponymsFile.face_value().c_str());
  
  pair<map<const char *,int>::iterator,bool> toponym_status;
  toponyms_t *toponyms;
  
  std::map<int, double> z0map;
  std::map<int, double>:: iterator it;
  pair<map<int, double>::iterator,bool> z0_status;
  
  status=grd_loadgrid (filename,&grid);
  if(status!=0) {
    check_error(-1, "error in toponyms file reading", __LINE__, __FILE__, 1);
    }
  code=new float[grid.nx*grid.ny];
  status=grd_loadr1 (filename,  grid, code, &mask);
  if(status!=0) {
    check_error(-1, "error in toponyms file reading", __LINE__, __FILE__, 1);
    }
//  status=grd_mirror_r (grid,grid.nx,code,mask); /// HERE !!!
  
  region=new float[mesh.nvtxs];
  
  for(n=0;n<mesh.nvtxs;n++) {
    t=mesh.vertices[n].lon;
    p=mesh.vertices[n].lat;
    t=map_recale(grid,t);
    status=map_nearestvalue00(grid, grid.nx,code,mask,t,p,&z);
    if(z==mask) {
      region[n]=1100;
      }
    else if(z==3000) {
      region[n]=1100;
      }
    else if(z==4000) {
      region[n]=1100;
      }
    else {
      region[n]=z;
      }
    }
    
  toponyms=new toponyms_t;
  status=geo_init_regions(toponyms);

  
  in=fopen(ValuesFilename,"r");
  if(in==0) {
    check_error(-1, "error in Z0 regional value file reading", __LINE__, __FILE__, 1);
    }
  while (!feof(in)) {
    nitems=fscanf(in,"%s %lf",toponym,&value);
    if(nitems!=2) break;
    z0_status=z0map.insert(std::pair<int,double>((*toponyms)[toponym],value));
    }
  fclose(in);

/**----------------------------------------------------------------------------
  initialisation */
  for(n=0;n<mesh.nvtxs;n++) {
    t=mesh.vertices[n].lon;
    p=mesh.vertices[n].lat;
    z0[n] = defaultZ0;
    }
/**----------------------------------------------------------------------------
  2nd level */
  for(n=0;n<mesh.nvtxs;n++) {
    t=mesh.vertices[n].lon;
    p=mesh.vertices[n].lat;
    flag=((int) NINT(region[n])) / 100;
//    it=z0map.find(flag);
    count=z0map.count(flag);
/**----------------------------------------------------------------------------
    check if flag registered in z0map */
    if(count>0) {
      z0[n] = z0map[flag];
      }
    if(z0[n]==0.) {
      printf("initializeZ0_ByRegions : trouble\n");
      }
    }
/**----------------------------------------------------------------------------
  3rd level */
  for(n=0;n<mesh.nvtxs;n++) {
    t=mesh.vertices[n].lon;
    p=mesh.vertices[n].lat;
    flag=((int) NINT(region[n])) / 10;
//    it=z0map.find(flag);
    count=z0map.count(flag);
/**----------------------------------------------------------------------------
    check if flag registered in z0map */
    if(count>0) {
      z0[n] = z0map[flag];
      }
    if(z0[n]==0.) {
      printf("initializeZ0_ByRegions : trouble\n");
      }
    }
  delete[] region;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void initialize_bottomfriction2D(mesh_t & mesh, parameter_t & data)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  using namespace Keywords;
  int n,n1,n2,status,zdim,udim;
  float *Cd,*z0;
  float *Cdmin;
  char filename[1024], rootname[1024];
  int  verbose=(gCPU_ID==gCPU_MASTER);
  
  if(verbose) {
    printf(SEPARATOR_1);
    printf("\ninitialize 2D friction coefficients\n");
    }
    
  status=paire_dimension(mesh,&zdim,&udim);
  
/**----------------------------------------------------------------------
  input data firts set on mesh vertices, might be revised*/
  udim=mesh.nvtxs;
  for(n = 0; n < udim; n++) {
    data.z0[n] = -1;
    data.Cd[n] = -1;
    data.FrictionRatio[n] = 1;
    }

  switch(gBottomFrictionType) {
    case FRICTION_QUADRATIC:
/**----------------------------------------------------------------------
      Quadratic friction with Cd formally prescribed*/
      Cd=new float[udim];
      if(!Cd) check_error(-1, "allocation failure", __LINE__, __FILE__, 1);
      if (tugo_cfg->dissipation->QFCPolygons.face_value() != "NONE")
/*----------------------------------------------------------------------
        variable Cd coefficient from polygon-wise distribution*/
        {/* STS: NO space-variable setting */}
      else if(strcmp(FrictionFile, "NONE") != 0) {
/*----------------------------------------------------------------------
        variable Cd coefficient from node-wise distribution*/
        strcpy(filename,FrictionFile);
        status = load_s2r(filename, Cd, udim, 0);
        if(status != 0) check_error(-1, "bottom friction coefficient file input failure", __LINE__, __FILE__, 1);
        }
      else {
/*----------------------------------------------------------------------
        uniform bottom friction coefficient */
        if(verbose) printf("use uniform bottom friction coefficient %f \n", gCd_default);
        for(n = 0; n < udim; n++) {
          Cd[n] = gCd_default;
          }
        }
      for(n = 0; n < udim; n++) {
        data.Cd[n] = Cd[n];
        }
      delete[] Cd;
      break;
      
    case FRICTION_LINEAR:
      for(n = 0; n < udim; n++) {
        data.rlinear[n] = gR_default;
        }
      
    case FRICTION_KARMAN:
/**----------------------------------------------------------------------
      Quadratic friction with Cd deduced from log velocity profile assumption*/
      z0=new float[udim];
      if (gz0_default <= 0) gz0_default = z0_def;
      printf(SEPARATOR_2);

      if (tugo_cfg->dissipation->RugosityPolygons.face_value() != "NONE") {
/*----------------------------------------------------------------------
        variable z0 coefficient from polygon-wise distribution*/
        /* STS: NO space-variable setting */
        }
      else if (tugo_cfg->dissipation->RugosityValuesbyRegions.face_value() != "NONE") {
/*----------------------------------------------------------------------
        variable z0 coefficient from region-wise distribution*/
        if(verbose) printf("initialize z0 from polygons\n");
        const char * Filename   = tugo_cfg->dissipation->RugosityValuesbyRegions.face_value().c_str();
        initializeZ0_ByRegions(Filename, gFEmesh[0], z0, gz0_default);
        }
      else {
/*----------------------------------------------------------------------
        uniform z0 (default) coefficient */
        if(verbose) printf("initialize z0 from default value\n");
        for (n=0; n < udim; n++) {
          z0[n] = gz0_default;
          }
        }

      Cdmin=new float[udim];
      if (tugo_cfg->dissipation->QFCMPolygons.face_value() != "NONE")
/*----------------------------------------------------------------------
        variable z0 coefficient from polygon-wise distribution*/
        {/* STS: NO space-variable setting */}
      else {
/*----------------------------------------------------------------------
        uniform z0 coefficient */
        for (n=0; n < udim; n++) {
          Cdmin[n] = gCd_minimum;
          }
        }

      for(n = 0; n < udim; n++) {
        data.z0[n] = z0[n];
        }
      delete[] z0;
      karman_Cd( data, data.h, udim, Cdmin);
      delete[] Cdmin;
      break;
      
    default:
      check_error(-1, "Invalid bottom friction type", __LINE__, __FILE__, 1) ;
      break;
    }
    
  if(gCPU_ID==gCPU_MASTER) {
    printf(SEPARATOR_2);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void initialize_u0 (mesh_t &mesh, parameter_t data2D, int fmt, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,nndes,status;
  float *u0,mask;
  char filename[1024], rootname[1024];
  
  verbose=(gCPU_ID==gCPU_MASTER);

  if(verbose) {
    printf(SEPARATOR_1);
    printf("\ninitialize background velocity\n");
    }

  rootname[0] = 0;

  switch(fmt) {
    case LGP0:
      nndes=mesh.ntriangles;
      break;
    case LGP1:
      nndes=mesh.nvtxs;
      break;
    case NCP1:
      nndes=mesh.nedges;
      break;
    }

/*----------------------------------------------------------------------
  variable background velocity */
  if (tugo_cfg->dissipation->BKGPolygons.face_value() != "NONE") {
/*----------------------------------------------------------------------
    variable u0 from polygon-wise distribution*/
    /* STS: NO space-variable setting */
    }
  else if(strcmp(BackgroundFile, "NONE") != 0) {
/*----------------------------------------------------------------------
    variable u0 from node-wise distribution*/
    /* STS: NO space-variable setting */
    }
  else {
    if(verbose) printf("use uniform background velocity %f m/s \n", gU0);
    for(n = 0; n < nndes; n++) {
      data2D.u0[n] = gU0;
      }
    }
/**----------------------------------------------------------------------
    Depth dependent background velocity effect: U0 at 100 m
    could be revived through options
#ifndef DRY
    state[i].tau = sqrt(u*u+v*v +1.e-30)+u0*sqrt(100./H);
#endif
*/

}
