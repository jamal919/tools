
/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#define MAIN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"

#include "tools-define.h"
#include "tools-structures.h"
#include "functions.h"
#include "tides.h"

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

extern  int cefmo(mesh_t *mesh, parameter_t data);

#warning check_error REDEFINED TO nc_check_error
#define check_error nc_check_error

#define T_GRID 0
#define U_GRID 1
#define V_GRID 2
#define F_GRID 3
#define W_GRID 4

#define COMODO    0
#define SYMTOOLS  1
#define SYMPHONIE 2
#define ROMS      3

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int initialise_UGO3D_NCP1_QLP0_QLP1(mesh_t mesh, ugo_state_t *state)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status = -1;
  int i, k, l, m, n, n1, n2;
  int type;
  double T, S;
  double R, r;
  double t, p;

/*-----------------------------------------------------------------------------
  NC-P1 nodes*/
  state->u     = new double*[mesh.nedges];
  state->v     = new double*[mesh.nedges];
  state->umean = new double [mesh.nedges];
  state->vmean = new double [mesh.nedges];
  state->U     = new double*[mesh.nedges];
  state->V     = new double*[mesh.nedges];
  for(n = 0; n < mesh.nedges; n++) {
    state->u[n]   = new double[mesh.nlayers];
    state->v[n]   = new double[mesh.nlayers];
    state->U[n]   = new double[mesh.nlayers];
    state->V[n]   = new double[mesh.nlayers];
    }

/*-----------------------------------------------------------------------------
  P0 nodes*/
  state->elevation = new double[mesh.ntriangles];
  state->T     = new double*[mesh.ntriangles];
  state->S     = new double*[mesh.ntriangles];
  state->rho   = new double*[mesh.ntriangles];
  for(n = 0; n < mesh.ntriangles; n++) {
    state->T[n]     = new double[mesh.nlayers];
    state->S[n]     = new double[mesh.nlayers];
    state->rho[n]   = new double[mesh.nlayers];
    }
  return (status);

error:
  status = -1;
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int initialise_UGO3D_NCP1_QLP1_QLP1(mesh_t mesh, ugo_state_t *state)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status = -1;
  int i, k, l, m, n, n1, n2;
  int type;
  double T, S;
  double R, r;
  double t, p;

/*-----------------------------------------------------------------------------
  NC-P1 nodes*/
  state->u     = new double*[mesh.nedges];
  state->v     = new double*[mesh.nedges];
  state->umean = new double [mesh.nedges];
  state->vmean = new double [mesh.nedges];
  state->U     = new double*[mesh.nedges];
  state->V     = new double*[mesh.nedges];
  for(n = 0; n < mesh.nedges; n++) {
    state->u[n]   = new double[mesh.nlayers];
    state->v[n]   = new double[mesh.nlayers];
    state->U[n]   = new double[mesh.nlayers];
    state->V[n]   = new double[mesh.nlayers];
    }
/*-----------------------------------------------------------------------------
  P1 nodes*/
  state->elevation = new double[mesh.nvtxs];
  state->T     = new double*[mesh.nvtxs];
  state->S     = new double*[mesh.nvtxs];
  state->rho   = new double*[mesh.nvtxs];
  for(n = 0; n < mesh.ntriangles; n++) {
    state->T[n]     = new double[mesh.nlevels];
    state->S[n]     = new double[mesh.nlevels];
    state->rho[n]   = new double[mesh.nlevels];
    }

  return (status);

error:
  status = -1;
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int initialise_UGO3D_NCP1_QLP1_LGP1(mesh_t mesh, ugo_state_t *state)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status = -1;
  int i, k, l, m, n, n1, n2;
  int type;
  double T, S;
  double R, r;
  double t, p;

/*-----------------------------------------------------------------------------
  NC-P1 nodes*/
  state->u     = new double*[mesh.nedges];
  state->v     = new double*[mesh.nedges];
  state->umean = new double [mesh.nedges];
  state->vmean = new double [mesh.nedges];
  state->U     = new double*[mesh.nedges];
  state->V     = new double*[mesh.nedges];
  for(n = 0; n < mesh.nedges; n++) {
    state->u[n]   = new double[mesh.nlayers];
    state->v[n]   = new double[mesh.nlayers];
    state->U[n]   = new double[mesh.nlayers];
    state->V[n]   = new double[mesh.nlayers];
    }
    
/*-----------------------------------------------------------------------------
  P1 nodes*/
  state->elevation = new double[mesh.nvtxs];
  state->T     = new double*[mesh.nvtxs];
  state->S     = new double*[mesh.nvtxs];
  state->rho   = new double*[mesh.nvtxs];
  for(n = 0; n < mesh.ntriangles; n++) {
    state->T[n]     = new double[mesh.nlayers];
    state->S[n]     = new double[mesh.nlayers];
    state->rho[n]   = new double[mesh.nlayers];
    }

  return (status);

error:
  status = -1;
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int initialise_UGO3D(mesh_t mesh, ugo_state_t *state, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n;
  int nlayers,status;


  switch(mode) {
    case 0:
/*-----------------------------------------------------------------------
      u NCP1, w QLP0, (T,S) LGP0*/
//      status=initialise_UGO3D_NCP1_QLP0_LGP0(mesh, state);
      break;

    case 1:
/*-----------------------------------------------------------------------
      u NCP1, w QLP0, (T,S) LGP1*/
//      check_error(-1, "paire not implemented yet", __LINE__, __FILE__, 1);
      break;

    case 2:
/*-----------------------------------------------------------------------
      u NCP1, w QLP0, (T,S) QLP1*/
//      status=initialise_UGO3D_NCP1_QLP0_QLP1(mesh, state);
      break;

    case NCP1xQLP1xLGP0:
/*-----------------------------------------------------------------------
      u NCP1, w QLP1, (T,S) LGP0*/
//      status=initialise_UGO3D_NCP1_QLP1_LGP0(mesh, state);
      break;

    case NCP1xQLP1xLGP1:
/*-----------------------------------------------------------------------
      u NCP1, w QLP1, (T,S) LGP1*/
      status=initialise_UGO3D_NCP1_QLP1_LGP1(mesh, state);
      break;


    case 6:
/*-----------------------------------------------------------------------
      u NCP1, w QLP0, (T,S) QLP1*/
//      status=initialise_UGO3D_NCP1_QLP1_QLP1(mesh, state);
      break;

    default:
//      check_error(-1, "unknown paire", __LINE__, __FILE__, 1);
      break;
    }


  return(status);
}




/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#define rsize 1000

  geo_t projection;
  int i,j,k,l,m,n,v,count,status;

  char *keyword;
  cdfgbl_t data_info,grid_info;
  cdfvar_t T_info,S_info,RHO_info;
  cdfvar_t u_info,v_info,w_info;
  cdfvar_t U_info,V_info,W_info;
  cdfvar_t SSE_info,SSF_info;
  cdfvar_t UBAR_info,VBAR_info;
  cdfvar_t depth_z_info;
  cdfvar_t mask_info,topo_info,f_topo_info;
//  cdfvar_t info[50];
  int tracer_discretisation=LGP1,frame=0,verbose;

  char *rootname=NULL,*output=NULL,*input=NULL,*discretisation=NULL,*bathymetry=NULL;
  char *symgridfile=NULL,*notebook=NULL,*meshbook=NULL,*poly=NULL;
  char *source=NULL;

  grid_t cmeshgrid,smeshgrid;
  grid_t grid[10],topogrid3d,extended,vgrid;
  mesh_t mesh;
  float  *z_topo=0,*z_landmask=0,*buffer[10],spec[10];
  float  *v_topo=0,*v_landmask=0;
  double scale,offset,x,y;
  float z;
  ugo_state_t state;
  float *z_landmask3D;
  
  plg_t *polygones=NULL;
  int npolygones=0;

  pocgrd_t ncgrid;
  cdfvar_t variable,variable_T,variable_S,variable_SSE;
  int vdepth;
  bool debug=false;
  
  int source_id, target;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        if(strcmp(keyword,"--debug")==0) {
          debug=true;
          n++;
          continue;
          }
        switch (keyword[1]) {
        case 'r' :
          rootname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'd' :
          discretisation= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'b' :
          bathymetry= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'g' :
          symgridfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          meshbook= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          poly= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          if(strcmp(keyword,"-source")==0) {
            source= strdup(argv[n+1]);
            n++;
            n++;
            break;
            }
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

  if(source==0) {
    printf("source format not passed in arguments, use symtools as a default\n");
    source=strdup("symtools");
    }
  if(strcmp(source,"comodo")==0) {
    source_id=COMODO;
    }
  else if(strcmp(source,"symtools")==0) {
    source_id=SYMTOOLS;
    }
  else if(strcmp(source,"symphonie")==0) {
    source_id=SYMPHONIE;
    }
  else if(strcmp(source,"roms")==0) {
    source_id=ROMS;
    }
    
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

  parameter_t data;
//  status= cefmo(&mesh, data);

  if(output==NULL) output=strdup("3Dmesh.nc");

  if(rootname==NULL) rootname= strdup("quaker");

/*-----------------------------------------------------------------------------
  build local FE mesh from notebook_grid */
  if(meshbook!=NULL) {
    printf("#################################################################\n");
    printf("load notebook file: %s\n",notebook);
    status=load_notebook(notebook, &cmeshgrid, &smeshgrid, &projection);
    printf("%s (mesh book file) processed\n",meshbook);
    }

  printf("#################################################################\n");
  printf("check netcdf input file: %s\n", input);
  verbose=0;
  status=cdf_globalinfo(input,&data_info,verbose);
  for (v=0;v<data_info.nvarsp;v++) {
    printf("variable %3d: name %s, type %d,ndim %d \n",v,(data_info.variable[v]).name,data_info.variable[v].type,data_info.variable[v].ndim);
    }
  status=cdf_globalinfo(symgridfile,&grid_info,verbose);

//  status=varinfo(input,"T",&T_info,1);
  status=cdf_varinfo(input,"tem",&T_info,1);
//  status=cdf_varinfo(input,"S",&S_info,1);
  status=cdf_varinfo(input,"sal",&S_info,1);
//  status=cdf_varinfo(input,"RHO",&RHO_info,1);
//   status=cdf_varinfo(input,"U",&U_info,1);
//   status=cdf_varinfo(input,"V",&V_info,1);
  status=cdf_varinfo(input,"vlx",&U_info,1);
  status=cdf_varinfo(input,"vly",&V_info,1);
  status=cdf_varinfo(input,"sse",&SSE_info,1);
  status=cdf_varinfo(input,"ssf",&SSF_info,1);
//   status=cdf_varinfo(input,"W",&W_info,1);
//   status=cdf_varinfo(input,"UBAR",&UBAR_info,1);
//   status=cdf_varinfo(input,"VBAR",&VBAR_info,1);


  if(T_info.id==-1) {
    printf("#################################################################\n");
    printf("temperature variable not found, look for tracer bathymetry: %s\n", input);
    switch (source_id) {
      case COMODO:
//        status=;
        break;
      case SYMTOOLS:
        status=cdf_varinfo(input,"topo",    &topo_info,1);
        status=cdf_varinfo(input,"landmask",&mask_info,1);
        status= poc_getgrid2d (input, grid_info, topo_info, &grid[T_GRID]);
        status= poc_getgrid2d (input, grid_info, mask_info, &grid[T_GRID]);
        break;
      case ROMS :
        status=cdf_varinfo(input,"h",       &topo_info,1);
        status=cdf_varinfo(input,"mask_rho",&mask_info,1);
/**----------------------------------------------------------------------------
        h is 2D variable */
        status= poc_getgrid2d (input, grid_info, topo_info, &grid[T_GRID]);
/**----------------------------------------------------------------------------
        mask_rho is 2D variable */
//        status= poc_getgrid3d (input, grid_info, mask_info, &grid[T_GRID],&vdepth);
        break;
      case SYMPHONIE:
        status=cdf_varinfo(input,"depth_w",&topo_info,1);
//        status=cdf_varinfo(input,"h_w",&topo_info,1);
        status=cdf_varinfo(input,"mask_t",&mask_info,1);
        status= poc_getgrid3d (input, grid_info, topo_info, &grid[W_GRID],&vdepth);
        status= poc_getgrid3d (input, grid_info, mask_info, &grid[T_GRID],&vdepth);
        break;
      }
    }
  else
    status= poc_getgrid3d (input, grid_info, T_info, &grid[T_GRID],&vdepth);

  if(status !=0) goto error;

  switch (source_id) {
    case SYMTOOLS:
      z_topo     = new float[grid[T_GRID].Hsize()];
      z_landmask = new float[grid[T_GRID].Hsize()];
      status= poc_getvar2d (input, topo_info.id, frame, z_topo,     &spec[0], topo_info);
      status= poc_getvar2d (input, mask_info.id, frame, z_landmask, &spec[1], mask_info);
      break;
    case SYMPHONIE :
      z_topo     = new float[grid[W_GRID].Hsize()];
      z_landmask = new float[grid[T_GRID].Hsize()];
      for(j=0;j<grid[W_GRID].ny;j++) {
        for(i=0;i<grid[W_GRID].nx;i++) {
          k=j*grid[W_GRID].nx+i;
/* *----------------------------------------------------------------------------
          use 1st (deepest) level for bathymetry */
          z_topo[k]=-grid[W_GRID].z[k];
          }
        }
//      status= poc_getvar3d (input, topo_info.id, frame, z_topo,     &spec[0] ,topo_info);
      status= poc_getvar2d (input, mask_info.id, frame, z_landmask, &spec[1], mask_info);
      break;
    case ROMS:
      z_topo       = new float[grid[T_GRID].Hsize()];
      status= poc_getvar2d (input, topo_info.id, frame, z_topo,     &spec[0], topo_info);
      z_landmask   = new float[grid[T_GRID].Hsize()];
      status= poc_getvar2d (input, mask_info.id, frame, z_landmask, &spec[1], mask_info);
//       z_landmask3D = new float[grid[T_GRID].nx*grid[T_GRID].ny*grid[T_GRID].nz];
//       status= poc_getvar3d (input, mask_info.id, frame, z_landmask3D, &spec[0] ,mask_info);
//       for(j=0;j<grid[T_GRID].ny;j++) {
//         for(i=0;i<grid[T_GRID].nx;i++) {
//           k=j*grid[T_GRID].nx+i;
// /* *----------------------------------------------------------------------------
//           use last (shallowest) level for landmask */
//           z_landmask[k]=z_landmask3D[k+grid[T_GRID].nx*grid[T_GRID].ny*(grid[T_GRID].nz-1)];
//           if(z_landmask[k]==-1) {
//             z_landmask[k]=0.;
//             }
//           if(z_landmask[k]==spec[0]) {
//             z_landmask[k]=0.;
//             }
// //           if(z_landmask[k]==0) {
// //             z_topo[k]=spec[0];
// //             }
//           }
//         }
      break;
    }
    
  goto next;

/**  temporary
-----------------------------------------------------------------------------*/
  
  for(j=0;j<grid[0].ny;j++)
    for(i=0;i<grid[0].nx;i++) {
/*-----------------------------------------------------------------------------
      levels downwards*/
      k=grid[0].ny*grid[0].nx*(grid[0].nz-1)+j*grid[0].nx+i;
/*-----------------------------------------------------------------------------
      levels upwards*/
      k=j*grid[0].nx+i;
      z_topo[k]=grid[0].z[k];
      }

  if(poly!=NULL) {
/*-----------------------------------------------------------------------------
    get landmask from polygons */
    sprintf(input,"%s.plg",poly);
    printf("#################################################################\n");
    printf("load polygons file: %s\n",poly);
    if(poly!=NULL) {
      sprintf(input,"%s.plg",poly);
      status=plg_load_scan(input, &polygones, &npolygones);
      if(status !=0) {__ERR_BASE_LINE__("exiting\n");exit(-1);}
      }
    else __ERR_BASE_LINE__("exiting\n");exit(-1);
//    landmask=set_landmask02(smeshgrid,polygones,npolygones);
    printf("land mask sucessfully completed\n");
    }
  else {
/*-----------------------------------------------------------------------------
    get landmask from variable mask */
    buffer[0]=new float[grid[0].nx*grid[0].ny*grid[0].nz];
    status= poc_getvar3d (input, T_info.id, frame, buffer[0], &spec[0] ,T_info);
    z_landmask=new float[grid[0].nx*grid[0].ny];
    for(j=0;j<grid[0].ny;j++) {
      for(i=0;i<grid[0].nx;i++) {
        k=j*grid[0].nx+i;
        if(buffer[0][k]==spec[0]) {
          z_landmask[k]=0.;
          }
        else {
          z_landmask[k]=1.;
          }
        }
      }
    }

next:
/*-----------------------------------------------------------------------------
  assume spherical grid*/
  projection.type=0;
 
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  build local FE mesh from regular grid
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  switch(target) {
    case 0:
/*------------------------------------------------------------------------------
      create mesh grid from tracer grid*/
/*------------------------------------------------------------------------------
      extend tracer grid*/
//       printf("#################################################################\n");
//       printf("extend tracer grid\n");
//       extended= map_extend_zgrid(grid[T_GRID],z_topo);
      grid[T_GRID].z=new double[grid[T_GRID].nx*grid[T_GRID].ny*grid[T_GRID].nz];
      for(j=0;j<grid[T_GRID].ny;j++) {
        for(i=0;i<grid[T_GRID].nx;i++) {
          k=j*grid[T_GRID].nx+i;
          grid[T_GRID].z[k]=z_topo[k];
          }
        }
      printf("\n#################################################################\n");
      printf("create mesh grid from tracer grid\n");
      status= fe_mimicgrid(grid[T_GRID], &mesh, z_landmask, z_topo, grid[T_GRID], z_landmask, "symmic.nei", debug);
      break;
    case 1:
/*------------------------------------------------------------------------------
      create vorticity grid from tracer grid*/
      status=cdf_varinfo(input,"h_f",&f_topo_info,1);
      status= poc_getgrid2d (input, grid_info, f_topo_info, &grid[F_GRID]);
      v_topo     = new float[grid[F_GRID].Hsize()];
      v_landmask = new float[grid[F_GRID].Hsize()];
      status= poc_getvar2d (input, f_topo_info.id, frame, v_topo, &spec[1], f_topo_info);
      status=cdf_varinfo(input,"mask2d_f",&mask_info,1);
      status= poc_getvar2d (input, mask_info.id, frame, v_landmask, &spec[1], mask_info);
      vgrid=grid[F_GRID];
//       printf("\n#################################################################\n");
//       printf("create voticity grid from tracer grid\n");
//       vgrid      = map_vgrid(grid[T_GRID], z_topo, spec[0],0);
//       v_topo     = map_extrapolate_z2v(grid[T_GRID], vgrid, z_topo,     spec[0]);
//       v_landmask = map_extrapolate_z2v(grid[T_GRID], vgrid, z_landmask, spec[1]);
      printf("#################################################################\n");
      printf("create mesh grid from vorticity grid\n");
      debug=true;
      status= fe_mimicgrid(vgrid, &mesh, v_landmask, v_topo, grid[T_GRID], z_landmask, "symmic.nei", debug);
      break;
    }

  switch (source_id) {
    case SYMTOOLS:
      __OUT_BASE_LINE__("only 2D conversion available, symmic stops here...\n");
      exit(0);
      break;
    case SYMPHONIE :
      break;
    case ROMS:
      __OUT_BASE_LINE__("(at the moment) only 2D conversion available, symmic stops here...\n");
      exit(0);
      break;
    }

  
  status= fe_savemeshNC3D("3Dmesh.nc",mesh, 1);

  scale  = 1.0;
  offset = 0.0;
  count  = 0;

  status=initialise_UGO3D(mesh,&state, NCP1xQLP1xLGP1);
/* *----------------------------------------------------------------------------
  T,S and RHO 3D variables */

  for(k=0;k<3;k++) buffer[k]=new float[grid[0].nx*grid[0].ny*grid[0].nz];

  status= poc_getvar3d (input, T_info.id,   frame, buffer[0], &spec[0] ,T_info);
  if(status!=0) goto error;
  status= poc_getvar3d (input, S_info.id,   frame, buffer[1], &spec[1] ,S_info);
  if(status!=0) goto error;
//  status= poc_getvar3d (input, RHO_info.id, frame, grid[0], buffer[2], &spec[2] ,RHO_info);

  buffer[0]=map_duplicate3D(grid[0],buffer[0],extended);
  buffer[1]=map_duplicate3D(grid[0],buffer[1],extended);

  grid[0]=extended;

  switch (tracer_discretisation) {
    case LGP0:
      for (n=0;n<mesh.ntriangles;n++) {
        status=fe_barycentre(mesh,n,&x,&y,0);
        for (k=0;k<mesh.nlayers;k++) {
          status= map_interpolation( grid[0], &(buffer[0][k*grid[0].nx*grid[0].ny]), spec[0],x,y,&z);
          state.T[n][k]  = z;
          status= map_interpolation( grid[0], &(buffer[1][k*grid[0].nx*grid[0].ny]), spec[1],x,y,&z);
          state.S[n][k]  = z;
//       status= map_interpolation( grid[0], &(buffer[2][k*grid[0].nx*grid[0].ny]), spec[2],x,y,&z);
//       state.rho[n][k]  = z;
          }
        }
      variable_T = poc_variable_UG4D("T", (double) spec[0], "C", scale, offset,"potential_sea_water_temperature","T","M","K");
      variable_S = poc_variable_UG4D("S", (double) spec[1], "PSU", scale, offset,"sea_water_sanility","T","M","K");
      break;

    case LGP1:
      for (n=0;n<mesh.nvtxs;n++) {
        x=mesh.vertices[n].lon;
        y=mesh.vertices[n].lat;
        for (k=0;k<mesh.nlayers;k++) {
          status= map_interpolation( grid[0], &(buffer[0][k*grid[0].nx*grid[0].ny]), spec[0],x,y,&z);
          state.T[n][k]  = z;
          status= map_interpolation( grid[0], &(buffer[1][k*grid[0].nx*grid[0].ny]), spec[1],x,y,&z);
          state.S[n][k]  = z;
//       status= map_interpolation( grid[0], &(buffer[2][k*grid[0].nx*grid[0].ny]), spec[2],x,y,&z);
//       state.rho[n][k]  = z;
          }
        }
      variable_T = poc_variable_UG4D("T", (double) spec[0], "C", scale, offset,"potential_sea_water_temperature","T","N","K");
      status   = cdf_createvariable(output, &(variable_T));
      variable_S = poc_variable_UG4D("S", (double) spec[1], "PSU", scale, offset,"sea_water_sanility","T","N","K");
      status   = cdf_createvariable(output, &(variable_S));
      break;
    }
  status   = poc_put_UG4D(output, mesh, count, variable_T.id, state.T);
  variable_T.destroy();

  status   = poc_put_UG4D(output, mesh, count, variable_S.id, state.S);
  variable_S.destroy();
//   variable = poc_variable_UG4D("RHO", (double) spec[2], "kg/m^3", scale, offset,"sea_water_density","T","M","K");
//   status   = cdf_createvariable(output, &(variable));
//   status   = poc_put_UG4D(output, mesh, count, variable.id, state.rho);
//   status   = free_ncvariable(&variable);
//   for(n = 0; n < mesh.ntriangles; n++) {
//     delete (state.T[n]);
//     delete (state.S[n]);
//     delete (state.rho[n]);
//     }
//   delete (state.T);
//   delete (state.S);
//   delete (state.rho);


/* *----------------------------------------------------------------------------
  u,v 3D variables */

  status= poc_getgrid3d (input, grid_info, U_info, &grid[1],&vdepth);
  status= poc_getvar3d (input, U_info.id, frame, buffer[0], &spec[0] ,U_info);

  status= poc_getgrid3d (input, grid_info, V_info, &grid[2],&vdepth);
  status= poc_getvar3d (input, V_info.id, frame, buffer[1], &spec[1] ,V_info);

  grid[3]= map_duplicategrid(grid[1]);
  grid[4]= map_duplicategrid(grid[2]);

  buffer[0]=map_duplicate3D(grid[1],buffer[0],grid[3]);
  buffer[1]=map_duplicate3D(grid[2],buffer[1],grid[4]);

  for (n=0;n<mesh.nedges;n++) {
    status=fe_position(mesh,mesh.edges[n],&x,&y,0);
    for (k=0;k<mesh.nlayers;k++) {
      status= map_interpolation( grid[1], &(buffer[0][k*grid[0].nx*grid[0].ny]), spec[0],x,y,&z);
      state.u[n][k]  = z;
      status= map_interpolation( grid[2], &(buffer[1][k*grid[0].nx*grid[0].ny]), spec[1],x,y,&z);
      state.v[n][k]  = z;
      }
    }
  variable = poc_variable_UG4D("u", (double) spec[0], "m/s", scale, offset,"eastward_velocity","T","E","K");
  status   = cdf_createvariable(output, &(variable));
  status   = poc_put_UG4D(output, mesh, count, variable.id, state.u);
  variable.destroy();

  variable = poc_variable_UG4D("v", (double) spec[1], "m/s", scale, offset,"northward_velocity","T","E","K");
  status   = cdf_createvariable(output, &(variable));
  status   = poc_put_UG4D(output, mesh, count, variable.id, state.v);
  variable.destroy();

  for(k=0;k<3;k++) delete[] buffer[k];

/* *----------------------------------------------------------------------------
  elevation 2D variables */

  status= poc_getgrid2d (input, grid_info, SSE_info, &grid[0]);

  for(k=0;k<3;k++) buffer[k]=new float[grid[0].nx*grid[0].ny];

  status= poc_getvar2d (input, SSE_info.id, frame, buffer[0], &(spec[0]) ,SSE_info);

  buffer[0]=map_duplicate2D(grid[0],buffer[0],extended);

  grid[0]=extended;

  switch (tracer_discretisation) {
    case LGP0:
      for (n=0;n<mesh.ntriangles;n++) {
        status=fe_barycentre(mesh,n,&x,&y,0);
        status= map_interpolation( grid[0], buffer[0], spec[0],x,y,&z);
        state.T[n][k]  = z;
        }
      variable_SSE = poc_variable_UG3D("sse", (double) spec[0], "m", scale, offset,"sea_surface_elevation","T","M");
      break;

    case LGP1:
      for (n=0;n<mesh.nvtxs;n++) {
        x=mesh.vertices[n].lon;
        y=mesh.vertices[n].lat;
        status= map_interpolation( grid[0], buffer[0], spec[0],x,y,&z);
        state.elevation[n]  = z;
        }
      variable_SSE = poc_variable_UG3D("sse", (double) spec[0], "m", scale, offset,"sea_surface_elevation","T","N");
      status   = cdf_createvariable(output, &(variable_SSE));
      break;
    }
  status   = poc_put_UG3D(output, mesh, count, variable_SSE.id, state.elevation);
  variable_SSE.destroy();

  for(k=0;k<3;k++) delete[] buffer[k];

/* *----------------------------------------------------------------------------
  mean uv 2D variables */

  status= poc_getvar3d (input, UBAR_info.id, frame, buffer[0], &spec[0] ,UBAR_info);
  status= poc_getvar3d (input, VBAR_info.id, frame, buffer[1], &spec[1] ,VBAR_info);

  for (n=0;n<mesh.nedges;n++) {
    status=fe_position(mesh,mesh.edges[n],&x,&y,0);
    status= map_interpolation( grid[0], &(buffer[0][k*grid[0].nx*grid[0].ny]), spec[0],x,y,&z);
    state.umean[n] = z;
    status= map_interpolation( grid[0], &(buffer[1][k*grid[0].nx*grid[0].ny]), spec[1],x,y,&z);
    state.vmean[n] = z;
    }

  variable = poc_variable_UG3D("ubar", (double) spec[0], "m/s", scale, offset,"eastward_mean_velocity","T","E");
  status   = cdf_createvariable(output, &(variable));
  status   = poc_put_UG3D(output, mesh, count, variable.id, state.umean);
  variable.destroy();

  variable = poc_variable_UG3D("vbar", (double) spec[1], "m/s", scale, offset,"northward_mean_velocity","T","E");
  status   = cdf_createvariable(output, &(variable));
  status   = poc_put_UG3D(output, mesh, count, variable.id, state.vmean);
  variable.destroy();

/*
  variable = poc_variable_UG4D("U", mask, "kg/m^2s", scale, offset,"eastward_momentum","T","E","K");
  status   = cdf_createvariable(output, &(variable));
  status   = poc_put_UG4D(output, mesh, count, variable.id, state.U);
  status   = free_ncvariable(&variable);

  variable = poc_variable_UG4D("V", mask, "kg/m^2s", scale, offset,"northward_momentum","T","E","K");
  status   = cdf_createvariable(output, &(variable));
  status   = poc_put_UG4D(output, mesh, count, variable.id, state.V);
  status   = free_ncvariable(&variable);

*/

/*-----------------------------------------------------------------------------
  QL-P0 nodes*/
  state.w     = new double*[mesh.ntriangles];
  state.W     = new double*[mesh.ntriangles];
  for(n = 0; n < mesh.ntriangles; n++) {
    state.w[n]   = new double[mesh.nlevels];
    state.W[n]   = new double[mesh.nlevels];
    }

//  status= load_field3d (input, w_info.id, frame, grid[0], buffer[0], &spec[0] ,w_info);
//  status= load_field3d (input, W_info.id, frame, grid[0], buffer[1], &spec[1] ,W_info);

  for (n=0;n<mesh.ntriangles;n++) {
    status=fe_barycentre(mesh,n,&x,&y,0);
    for (k=0;k<mesh.nlevels;k++) {
      state.w[n][k]  = 0.0;
      state.W[n][k]  = 0.0;
      }
    }

  variable = poc_variable_UG4D("w", (double) spec[2], "m/s", scale, offset,"vertical_velocity","T","M","L");
  status   = cdf_createvariable(output, &(variable));
  status   = poc_put_UG4D(output, mesh, count, variable.id, state.w);
  variable.destroy();
/*

  variable = poc_variable_UG4D("W", mask, "kg/m^2s", scale, offset,"vertical_momentum","T","M","L");
  status   = cdf_createvariable(output, &(variable));
  status   = poc_put_UG4D(output, mesh, count, variable.id, state.W);
  status   = free_ncvariable(&variable);
*/
/*-----------------------------------------------------------------------------
  P1 nodes*/
  state.elevation = new double[mesh.nvtxs];

//  status= load_field3d (input, w_info.id, frame, grid[0], buffer[0], &spec[0] ,w_info);

  for (n=0;n<mesh.nvtxs;n++) {
    status=fe_position(mesh,mesh.vertices[n],&x,&y,0);
    state.elevation[n]=0.0;
    state.elevation[n] = -mesh.vertices[n].sigma[mesh.nlevels-1]*mesh.vertices[n].h;
    }

  variable = poc_variable_UG3D("elevation", (double) spec[0], "m", scale, offset,"sea_surface_elevation","T","N");
  status   = cdf_createvariable(output, &(variable));
  status   = poc_put_UG3D(output, mesh, count, variable.id, state.elevation);
  variable.destroy();

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
