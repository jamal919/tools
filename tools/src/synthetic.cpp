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

#include "fe.h"
#include "map.h"
#include "geo.h"
#include "sym-io.h"
#include "polygones.h"
#include "grd.h"
#include "netcdf-proto.h"

#define nstat 12
extern tidal_wave wSa;
extern float *set_earthquake_generic(grid_t );
extern int fe_setsigma03(mesh_t * mesh, int nlayers, int type);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_mimicgrid02(grid_t cgrid, int *incidence, mesh_t *mesh, geo_t projection, float *landmask, float *topo)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k,l,kk,ll,m,n,nn;
  int    nitems,type,status;
  int    nz,type1,type2;
  FILE   *in,*out;
  double lon,lat,angle,radius,x,y;
  double tx,ty,nx,ny,x_shift,y_shift;
  char   line[1024],nodefile[1024];

  mesh->destroy();

/*------------------------------------------------------------------------
  minimize band width by numbering along the shortest side*/
  mesh->nvtxs=0;
  if(cgrid.nx<cgrid.ny) {
    for(j=0;j<cgrid.ny;j++) {
      for(i=0;i<cgrid.nx;i++) {
        k=cgrid.nx*j+i;
        if(landmask[k]==1.) {
          incidence[k]=mesh->nvtxs;
          mesh->nvtxs++;
          }
        }
      }
    }
  else {
    for(i=0;i<cgrid.nx;i++) {
      for(j=0;j<cgrid.ny;j++) {
        k=cgrid.nx*j+i;
        if(landmask[k]==1.) {
          incidence[k]=mesh->nvtxs;
          mesh->nvtxs++;
          }
        }
      }
    }

  mesh->vertices= new vertex_t[mesh->nvtxs];
  mesh->nlayers=cgrid.nz-1;
  mesh->nlevels=cgrid.nz;

  mesh->nnghm=6;

  for(j=0;j<cgrid.ny;j++) {
    for(i=0;i<cgrid.nx;i++) {
      m=cgrid.nx*j+i;
      if(landmask[m]==1.) {
        x= cgrid.x[m];
        y= cgrid.y[m];
        n=incidence[m];
        status=geo_mercator_inverse(projection,&(mesh->vertices[n].lon),&(mesh->vertices[n].lat),x,y);
//        mesh->vertices[n].lon=x;
//        mesh->vertices[n].lat=y;
        mesh->vertices[n].nngh=0;
        mesh->vertices[n].h=-topo[m];
        mesh->vertices[n].code=0;
        mesh->vertices[n].ngh=(int *) malloc(mesh->nnghm*sizeof(int));
        for(k=0;k<mesh->nnghm;k++) mesh->vertices[n].ngh[k]=-1;
/*------------------------------------------------------------------------
        set vertical layers position*/
        if(cgrid.z!=NULL) {
          mesh->vertices[n].zlevels=new double[mesh->nlevels];
          for (l=0;l<mesh->nlevels;l++) {
            mesh->vertices[n].zlevels[l]=cgrid.z[m+l*cgrid.ny*cgrid.nx]/topo[m];
            }
          }
        }
      }
    }

  if(cgrid.z==NULL) status=fe_setsigma03(mesh, mesh->nlayers,type);

  for(j=0;j<cgrid.ny;j++) {
    for(i=0;i<cgrid.nx;i++) {
      k=cgrid.nx*j+i;
      if(landmask[k]!=1.) continue;
      kk=incidence[k];
/*------------------------------------------------------------------------
      righthand side neighbour*/
      if(i!=cgrid.nx-1) {
        l=cgrid.nx*j+i+1;
        ll=incidence[l];
        if(landmask[l]==1.) {
          mesh->vertices[kk].ngh[mesh->vertices[kk].nngh]=ll;
          mesh->vertices[ll].ngh[mesh->vertices[ll].nngh]=kk;
          mesh->vertices[kk].nngh++;
          mesh->vertices[ll].nngh++;
          }
        }
/*------------------------------------------------------------------------
      upper righthand side neighbour*/
      if((i!=cgrid.nx-1)&&(j!=cgrid.ny-1)) {
        l=cgrid.nx*(j+1)+i+1;
        ll=incidence[l];
        if(landmask[l]==1.) {
          mesh->vertices[kk].ngh[mesh->vertices[kk].nngh]=ll;
          mesh->vertices[ll].ngh[mesh->vertices[ll].nngh]=kk;
          mesh->vertices[kk].nngh++;
          mesh->vertices[ll].nngh++;
          }
        }
/*------------------------------------------------------------------------
      roof neighbour*/
      if(j!=cgrid.ny-1) {
        l=cgrid.nx*(j+1)+i;
        ll=incidence[l];
        if(landmask[l]==1.) {
          mesh->vertices[kk].ngh[mesh->vertices[kk].nngh]=ll;
          mesh->vertices[ll].ngh[mesh->vertices[ll].nngh]=kk;
          mesh->vertices[kk].nngh++;
          mesh->vertices[ll].nngh++;
          }
        }
      }
    }

  for(n=0;n<mesh->nvtxs;n++) {
    if(mesh->vertices[n].nngh<2) {
      status=fe_removevertex(mesh,n);
      n=-1;
      }
    }

  status=fe_list(mesh);
  status=fe_savemesh("test.nei",MESH_FILE_FORMAT_TRIGRID,*mesh);

  for(n=0;n<mesh->nvtxs;n++) {
    free(mesh->vertices[n].ngh);
    }
  status= fe_e2n (mesh);
  status=fe_list(mesh);
 
  status=fe_edgetable(mesh,0,0);
  status=fe_codetable2(mesh,0,1,0);
  status=fe_savemesh("regular.nei",MESH_FILE_FORMAT_TRIGRID,*mesh);
  
/*-----------------------------------------------------------------------------
  3D part */
  status= fe_savemeshNC3D("3Dmesh.nc",*mesh, 1);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_mimicgrid01(grid_t cgrid, mesh_t *mesh, geo_t projection, float *landmask, float *topo)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k,l,kk,ll,n=1024,nn;
  int    nitems,status;
  int    nz,type1,type2;
  FILE   *in,*out;
  double lon,lat,angle,radius,x,y;
  double tx,ty,nx,ny,x_shift,y_shift;
  char   line[1024],nodefile[1024];
  int *incidence;

  mesh->destroy();

  incidence=(int *)malloc(cgrid.nx*cgrid.ny*sizeof(int));

  mesh->nvtxs=0;
  if(cgrid.nx<cgrid.ny) {
    for(j=0;j<cgrid.ny;j++) {
      for(i=0;i<cgrid.nx;i++) {
        k=cgrid.nx*j+i;
        if(landmask[k]==1.) {
          incidence[k]=mesh->nvtxs;
          mesh->nvtxs++;
          }
        }
      }
    }
  else {
    for(i=0;i<cgrid.nx;i++) {
      for(j=0;j<cgrid.ny;j++) {
        k=cgrid.nx*j+i;
        if(landmask[k]==1.) {
          incidence[k]=mesh->nvtxs;
          mesh->nvtxs++;
          }
        }
      }
    }

  mesh->vertices=new vertex_t[mesh->nvtxs];

  mesh->nnghm=6;
  for(j=0;j<cgrid.ny;j++) {
    for(i=0;i<cgrid.nx;i++) {
      k=cgrid.nx*j+i;
      if(landmask[k]==1.) {
        x= cgrid.x[k];
        y= cgrid.y[k];
        kk=incidence[k];
        status=geo_mercator_inverse(projection,&(mesh->vertices[kk].lon),&(mesh->vertices[kk].lat),x,y);
        mesh->vertices[kk].nngh=0;
        mesh->vertices[kk].h=-topo[k];
        mesh->vertices[kk].code=0;
        mesh->vertices[kk].ngh=(int *) malloc(mesh->nnghm*sizeof(int));
        for(nn=0;nn<mesh->nnghm;nn++) mesh->vertices[kk].ngh[nn]=-1;
        }
      }
    }

  for(j=0;j<cgrid.ny;j++) {
    for(i=0;i<cgrid.nx;i++) {
      k=cgrid.nx*j+i;
      if(landmask[k]!=1.) continue;
      kk=incidence[k];
/*------------------------------------------------------------------------
      righthand side neighbour*/
      if(i!=cgrid.nx-1) {
        l=cgrid.nx*j+i+1;
        ll=incidence[l];
        if(landmask[l]==1.) {
          mesh->vertices[kk].ngh[mesh->vertices[kk].nngh]=ll;
          mesh->vertices[ll].ngh[mesh->vertices[ll].nngh]=kk;
          mesh->vertices[kk].nngh++;
          mesh->vertices[ll].nngh++;
          }
        }
/*------------------------------------------------------------------------
      upper righthand side neighbour*/
      if((i!=cgrid.nx-1)&&(j!=cgrid.ny-1)) {
        l=cgrid.nx*(j+1)+i+1;
        ll=incidence[l];
        if(landmask[l]==1.) {
          mesh->vertices[kk].ngh[mesh->vertices[kk].nngh]=ll;
          mesh->vertices[ll].ngh[mesh->vertices[ll].nngh]=kk;
          mesh->vertices[kk].nngh++;
          mesh->vertices[ll].nngh++;
          }
        }
/*------------------------------------------------------------------------
      roof neighbour*/
      if(j!=cgrid.ny-1) {
        l=cgrid.nx*(j+1)+i;
        ll=incidence[l];
        if(landmask[l]==1.) {
          mesh->vertices[kk].ngh[mesh->vertices[kk].nngh]=ll;
          mesh->vertices[ll].ngh[mesh->vertices[ll].nngh]=kk;
          mesh->vertices[kk].nngh++;
          mesh->vertices[ll].nngh++;
          }
        }
      }
    }

  for(n=0;n<mesh->nvtxs;n++) {
    if(mesh->vertices[n].nngh<2) {
      status=fe_removevertex(mesh,n);
      n=-1;
      }
    }

  status=fe_list(mesh);
  status=fe_savemesh("test.nei",MESH_FILE_FORMAT_TRIGRID,*mesh);

  for(n=0;n<mesh->nvtxs;n++) {
    free(mesh->vertices[n].ngh);
    }
  status= fe_e2n (mesh);
  status=fe_list(mesh);
 
  status=fe_edgetable(mesh,0,0);
  status=fe_codetable2(mesh,0,1,0);
  status=fe_savemesh("regular.nei",MESH_FILE_FORMAT_TRIGRID,*mesh);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tsunami(char *notebook,char *meshbook,char *poly,char *bathymetry,char *rootname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  z;
  geo_t projection;

  float  topo_mask=1.e+10,mask=1.e+10;
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;

  size_t size;
  FILE *in,*out;
  char *keyword,*zone;
  int flag;

  char *output=NULL,*input=NULL,*pathname=NULL;
  char *bathy=NULL;

  grid_t cgrid,sgrid,cmeshgrid,smeshgrid;
  grid_t grid,topogrid;
  mesh_t mesh;
  float  *landmask,*topo, *earthquake, *topobase,*tmp,*motion;
  float *dmdx,*dmdy,*dmdn;

  plg_t *polygones=NULL;
  int npolygones=0;

  pocgrd_t ncgrid;
  cdfvar_t variable;
  int dimlgth[4]={10,20,30,40};

  input =(char *) malloc(1024);
  output=(char *) malloc(1024);

  if(rootname==NULL) rootname= strdup("quaker");

/*-----------------------------------------------------------------------------
  Read notebook data */
  if(notebook==NULL) notebook= strdup("notebook_grid");

  status=load_notebook(notebook, &cgrid, &sgrid, &projection);
  printf("%s (notebook file) processed\n",notebook);

/*-----------------------------------------------------------------------------
  build local FE mesh from notebook_grid */
  if(meshbook==NULL) meshbook= strdup("notebook_mesh");

  status=load_notebook(notebook, &cmeshgrid, &smeshgrid, &projection);
  printf("%s (mesh book file) processed\n",meshbook);

/*-----------------------------------------------------------------------------
  get landmask from polygons */
  if(poly!=NULL) {
    sprintf(input,"%s.plg",poly);
    if(poly!=NULL) {
      sprintf(input,"%s.plg",poly);
      status=plg_load_scan(input, &polygones, &npolygones);
      if(status !=0) {__ERR_BASE_LINE__("exiting\n");exit(-1);}
      }
    else __ERR_BASE_LINE__("exiting\n");exit(-1);
    landmask=set_landmask02(smeshgrid,polygones,npolygones);
    printf("land mask sucessfully completed\n");
    }
  else {
    landmask=(float *) malloc(smeshgrid.nx*smeshgrid.ny*sizeof(float));
    for(j=0;j<sgrid.ny;j++)
      for(i=0;i<sgrid.nx;i++) {
        k=j*sgrid.nx+i;
        landmask[k]=1.;
        }
    }

/*-----------------------------------------------------------------------------
  set the earthquake characteristic*/
  earthquake=set_earthquake_generic(cgrid);

/*-----------------------------------------------------------------------------
  read topo database */
  sprintf(input,bathymetry);

  status=grd_loadgrid(input,&topogrid);
  if(status !=0) {
    __OUT_BASE_LINE__("cannot load bathymetry file=%f\n",input);
    exit(-1);
    }

  topogrid.xmin=floor(smeshgrid.xmin/topogrid.dx-1.0)*topogrid.dx;
  topogrid.xmax=floor(smeshgrid.xmax/topogrid.dx+1.0)*topogrid.dx;
  topogrid.ymin=floor(smeshgrid.ymin/topogrid.dy-1.0)*topogrid.dy;
  topogrid.ymax=floor(smeshgrid.ymax/topogrid.dy+1.0)*topogrid.dy;

  topogrid.nx=(int)( ((topogrid.xmax-topogrid.xmin)/topogrid.dx)+1 );
  topogrid.ny=(int)( ((topogrid.ymax-topogrid.ymin)/topogrid.dx)+1 );

  tmp=(float *) malloc(topogrid.nx*topogrid.ny*sizeof(float));
  topogrid.modeH=0;
  status=  grd_extract(input,topogrid,topogrid.nx,tmp, &mask);
  status=  grd_mirror_r(topogrid,topogrid.nx,tmp,topo_mask);

  topobase=tmp;

  //status=map_completegridaxis_2(&topogrid);
  status=map3d_completegridaxis(&topogrid);
/*------------------------------------------------------------------------------
  Interpolate and store topo on symphonie grid*/

  mask=-99999.;
  topo  =(float *) malloc(smeshgrid.nx*smeshgrid.ny*sizeof(float));
  motion=(float *) malloc(smeshgrid.nx*smeshgrid.ny*sizeof(float));

//  status=copy_grid_to_grid1d(topogrid,&topogrid1d);

  for(j=0;j<smeshgrid.ny;j++)
    for(i=0;i<smeshgrid.nx;i++) {
      k=j*smeshgrid.nx+i;
      x=smeshgrid.x[k];
      y=smeshgrid.y[k];
      status=map_interpolation(topogrid, topobase,topo_mask,x,y,&z);
      if (z!=topo_mask)
        topo[k]=z;
      else
        topo[k]=mask;
      status=map_interpolation(sgrid, earthquake,mask,x,y,&z);
      if (z!=mask)
        motion[k]=z;
      else
        motion[k]=0;
      }
  printf("bathymetry sucessfully completed\n");

/*-----------------------------------------------------------------------------
  build local FE mesh from notebook_grid */
  status=fe_mimicgrid01(cmeshgrid,&mesh,projection,landmask,topo);
  dmdx=(float *) malloc(smeshgrid.nx*smeshgrid.ny*sizeof(float));
  dmdy=(float *) malloc(smeshgrid.nx*smeshgrid.ny*sizeof(float));
  dmdn=(float *) malloc(smeshgrid.nx*smeshgrid.ny*sizeof(float));
  status= map_gradient(smeshgrid,smeshgrid.nx, motion, mask,0, dmdx,dmdy);


/*-----------------------------------------------------------------------------
  save netcdf file */
  sprintf(output,"%s.spherical.nc",rootname);

  status= poc_createfile(output);

  grid= map_getgrid3d(smeshgrid);
  status=poc_sphericalgrid_xyzt(output,"","",grid,&ncgrid);

  poc_standardvariable_xy(&variable,
    "deformation",mask,"m",1., 0.,"sea_floor_deformation","sea floor deformation","sfd",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  smeshgrid, variable.id,motion);
  variable.destroy();

  poc_standardvariable_xy(&variable,
    "dmdx",mask,"m",1., 0.,"sea_floor_deformation_x_gradient","sea floor deformation x-gradient","dmdx",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  smeshgrid, variable.id,dmdx);
  variable.destroy();

  poc_standardvariable_xy(&variable,
    "dmdy",mask,"m",1., 0.,"sea_floor_deformation_y_gradient","sea floor deformation y-gradient","dmdy",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  smeshgrid, variable.id,dmdy);
  variable.destroy();

  poc_standardvariable_xy(&variable,
    "bathymetry",mask,"m",1., 0.,"bathymetry","bathymetry","bathymetry",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  smeshgrid, variable.id,topo);
  variable.destroy();

  free(topo);
  free(landmask);

  __ERR_BASE_LINE__("exiting\n");exit(0);

 error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int symphonie(char *notebook,char *meshbook,char *poly,char *bathymetry,char *rootname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  z;
  geo_t projection;

  float  topo_mask=1.e+10,mask=1.e+10;
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;
  int nitems,*incidence;

  size_t size;
  FILE *in,*out;
  char *keyword,*zone;
//  int flag;

  char *output=NULL,*input=NULL,*pathname=NULL;
  char *bathy=NULL;

  grid_t   cgrid,sgrid,cmeshgrid,smeshgrid,topogrid1d;
  grid_t grid,topogrid;
  mesh_t mesh;
  float  *landmask,*topo, *earthquake, *topobase,*tmp,*motion;
  float *dmdx,*dmdy,*dmdn;

  plg_t *polygones=NULL;
  int npolygones=0;

  pocgrd_t ncgrid;
  cdfvar_t variable;
  int dimlgth[4]={10,20,30,40};
  char flag;

  input =(char *) malloc(1024);
  output=(char *) malloc(1024);

  if(rootname==NULL) rootname= strdup("quaker");

/*-----------------------------------------------------------------------------
  build local FE mesh from notebook_grid */
  if(meshbook==NULL) meshbook= strdup("notebook_mesh");

  status=load_notebook(notebook, &cmeshgrid, &smeshgrid, &projection);
  printf("%s (mesh book file) processed\n",meshbook);

  landmask=(float *) malloc(smeshgrid.nx*smeshgrid.ny*sizeof(float));
  for(j=0;j<smeshgrid.ny;j++)
    for(i=0;i<smeshgrid.nx;i++) {
      k=j*smeshgrid.nx+i;
      landmask[k]=0;
      }

  topo=(float *) malloc(smeshgrid.nx*smeshgrid.ny*sizeof(float));
  for(j=0;j<smeshgrid.ny;j++)
    for(i=0;i<smeshgrid.nx;i++) {
      k=j*smeshgrid.nx+i;
      topo[k]=-1.;
      }

/*-----------------------------------------------------------------------------
  read topo database */
  sprintf(output,"bathycote_in.dat");
  out=fopen(output,"r");

  for(i=1;i<smeshgrid.nx-1;i++) {
    for(j=1;j<smeshgrid.ny-1;j++) {
      k=smeshgrid.nx*j+i;
      nitems=fscanf(out,"%c",&flag);
//      nitems=fscanf(out,"%1f",&landmask[k]);
//      printf("%d %d %d %d %d\n",i,j,k,nitems,flag);
      if(nitems==0) {
        printf("%d %d %d\n",k,nitems,landmask[k]);
        }
      switch (flag) {
        case '0':
          landmask[k]=0;
          break;
        case '1':
          landmask[k]=1;
          break;
        default:
          j--;
          break;
        }
      }
    }
  for(i=1;i<smeshgrid.nx-1;i++) {
    for(j=1;j<smeshgrid.ny-1;j++) {
      k=smeshgrid.nx*j+i;
      fscanf(out,"%f",&topo[k]);
      }
    }

  fclose(out);

/*------------------------------------------------------------------------------
  Interpolate and store topo on symphonie grid*/

/*-----------------------------------------------------------------------------
  build local FE mesh from notebook_grid */
  incidence=new int[smeshgrid.nx*smeshgrid.ny];
  status=fe_mimicgrid02(cmeshgrid,incidence,&mesh,projection,landmask,topo);
  free(topo);
  free(landmask);

  __ERR_BASE_LINE__("exiting\n");exit(0);

 error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  z;
  geo_t projection;

  float  topo_mask=1.e+10,mask=1.e+10;
  double  x,y;
  double  r,s,Lx,Ly;

  int i,j,k,l,m,n,status;

  size_t size;
  FILE *in,*out;
  char *keyword,*zone;
  int flag;

  char *rootname=NULL,*output=NULL,*input=NULL,*pathname=NULL,*bathymetry=NULL;
  char *bathy=NULL,*notebook=NULL,*meshbook=NULL,*poly=NULL;

  grid_t   cgrid,sgrid,cmeshgrid,smeshgrid,topogrid1d;
  grid_t grid,topogrid;
  mesh_t mesh;
  float  *landmask,*topo, *earthquake, *topobase,*tmp,*motion;
  float *dmdx,*dmdy,*dmdn;

  plg_t *polygones=NULL;
  int npolygones=0;

  pocgrd_t ncgrid;
  cdfvar_t variable;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'r' :
          rootname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'd' :
          pathname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'b' :
          bathymetry= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'n' :
          notebook= strdup(argv[n+1]);
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

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        if(input==NULL) {
          notebook= strdup(argv[n]);
          n++;
          }
        else
          {
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
          }
        break;
      }
      free(keyword);
    }


//  tsunami(notebook,meshbook,poly,bathymetry,rootname);
  status=symphonie(notebook,meshbook,poly,bathymetry,rootname);
  __ERR_BASE_LINE__("exiting\n");exit(0);

 error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
  
}
