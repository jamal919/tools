/**************************************************************************

  Nonlinear finite element time stepping model

  David Greenberg    Bedford Institute of Oceanography   October 1991

***************************************************************************/

#define MAIN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "tools-structures.h"

#include "geo.h"
#include "fe.h"
#include "archive.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "functions.h"

#include "map.def"

#include "zapper.h"     /*  rutin.h contains common utility routines  */

//extern  int dMassMatrix(mesh_t * mesh, double **M, ordering_t *ordering, int sproduct2D);
extern  void fe_LGP0_to_LGP1(double *in, double *out, mesh_t mesh);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_topo01(char *input, char *output, mesh_t mesh, grid_t topogrid, float *topo, float topomask, char *dataset)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,k,l,m,n,status;
  double t,p;
  float h;
  float z,*depths=NULL,*slope_x=NULL,*slope_y=NULL;
  float *topo_x=NULL,*topo_y=NULL;
  size_t size;
  char *comments[2];

  comments[0]=new char[1024];
  comments[1]=new char[1024];

  depths=new float[mesh.nvtxs];

//  status=quoddy_loadr1(input, mesh.nvtxs, depths);

  if(input!=NULL) {
    status=quoddy_loadr1(input, mesh.nvtxs, depths);
    if(status!=0) return(-1);
    }
  else {
    for (n=0;n<mesh.nvtxs;n++) {
      depths[n]=mesh.vertices[n].h;
      }
    }

  for (n=0;n<mesh.nvtxs;n++) {
    t=mesh.vertices[n].lon;
    p=mesh.vertices[n].lat;
    t=map_recale(topogrid,t);
/* *----------------------------------------------------------------------------
    minimu elevation is positive, model depth also !!! */
    status=map_interpolation(topogrid,topo,topomask,t,p,&z);
    if(depths[n]==0.0) {
      if(z!=topomask) {
        depths[n]=MAX(2.,1.5*z);
        }
      else {
        depths[n]=2;
        }
      }
//    depths[n]=-depths[n];
    }

  sprintf(comments[0],"LGP1 bathymetry, unit in meters");
  sprintf(comments[1],"Created by mesh-zmin, direct interpolation from %s",dataset);

  if(output==NULL) output=strdup("topo-0.s2r");
  status=quoddy_saver1(output, mesh.nvtxs, depths,comments);

  delete[] depths;

  delete[] comments[0];
  delete[] comments[1];

  printf("fe_topo01 completed\n");

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_select_nodes(mesh_t mesh, char *poly, int *selected)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,m,n,n1,n2,status;
  int   count=0;
  int   inside;
  edge_t *edges=NULL;
  triangle_t *elt=NULL,*ptr=NULL;
  double t,p,d,t1,t2,p1,p2;
  double H1,H2,c;
  plg_t *polygones=NULL;
  int npolygones=0;

/*-----------------------------------------------------------------------------
  */
  printf ("check polygon criterion...\n");
  if(poly!=NULL) plg_load_scan(poly, &polygones, &npolygones);
  else return(0);
  for (n=0;n<mesh.nvtxs;n++) {
    t=mesh.vertices[n].lon;
    p=mesh.vertices[n].lat;
    inside=plg_TestInterior(t,p,polygones,npolygones);
    if(inside==-1) {
      selected[n]=0;
      }
    else {
      selected[n]=1;
      count++;
      }
    }

  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_topo02(char *input,char *output, char *poly,mesh_t mesh,float zmin)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,k,l,m,n,status;
  double t,p;
  float h;
  float z,*depths=NULL,*slope_x=NULL,*slope_y=NULL;
  float *topo_x=NULL,*topo_y=NULL;
  size_t size;
  char *comments[2];
  int *selected=NULL,count;

  comments[0]=new char[1024];
  comments[1]=new char[1024];

  depths=new float[mesh.nvtxs];
  selected=new int[mesh.nvtxs];
//bathymetry/base-10/
  if(input!=NULL) {
    status=quoddy_loadr1(input, mesh.nvtxs, depths);
    if(status!=0) return(-1);
    }
  else {
    for (n=0;n<mesh.nvtxs;n++) {
      depths[n]=mesh.vertices[n].h;
      }
    }

  count=fe_select_nodes(mesh, poly, selected);

  for (n=0;n<mesh.nvtxs;n++) {
    t=mesh.vertices[n].lon;
    p=mesh.vertices[n].lat;
    if(depths[n]<zmin) {
      if(selected[n]==1) {
        depths[n]=zmin;
        }
      }
    }

  sprintf(comments[0],"LGP1 bathymetry, unit in meters");
  sprintf(comments[1],"Created by mesh-zmin, zmin=%f prescribed in %s",zmin,poly);

  if(output==NULL) output=strdup("topo-0.s2r");
  status=quoddy_saver1(output, mesh.nvtxs, depths,comments);

  delete[] depths;
  delete[] selected;

  delete[] comments[0];
  delete[] comments[1];

  printf("fe_topo02 completed\n");

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nndes,status,*elts=NULL,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int option,channels=0;
  FILE *file=NULL;
  FILE *out=NULL;
  char *keyword=NULL,*zone=NULL,*gridfile=NULL,*input=NULL;
  char *meshfile=NULL,*depthfile=NULL,*output=NULL,*poly=NULL;
  mesh_t mesh;
  double dmax=0,cmax=0;
  int *selected=NULL,*targeted=NULL;
  grid_t topogrid;
  float *topo=NULL,topomask,h,zmin;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'b' :
          depthfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          poly= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          if(strcmp(keyword,"-hmin")==0) {
            sscanf(argv[n+1],"%f",&zmin);
            n++;
            n++;
            continue;
            }

          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        input= strdup(argv[n]);
        n++;
        break;
      }
    free(keyword);
    }

/* *----------------------------------------------------------------------
  load and initialize mesh structure*/
  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status!=0) {
      printf("unable to read the original mesh in %s\n",meshfile);
      goto error;
      }
    status=fe_list(&mesh);
    if(status!=0) {
      printf("unable to build the element list from the original mesh\n");
      goto error;
      }
    }
 else {
   printf("no mesh file specified; abort...\n");
   goto error;
   }

  status= fe_edgetable(&mesh,0,0);
  if(status!=0) {
    printf("unable to build the edge list from the original mesh\n");
    goto error;
    }

/* *----------------------------------------------------------------------
  load zmin grid*/
  if(depthfile!=NULL){
    status=grd_loadgrid(depthfile,&topogrid);
    if(status !=0) {
      __OUT_BASE_LINE__("cannot load zmin file=%f\n",depthfile);
      exit(-1);
      }
/* *----------------------------------------------------------------------
    load zmin data*/

    exitIfNull(
      topo=(float *) malloc(topogrid.nx*topogrid.ny*sizeof(float))
      );
    topogrid.modeH=0;
    status=  grd_loadr1(depthfile,topogrid,topo,&topomask);
    status= fe_topo01(input, output, mesh, topogrid, topo, topomask, depthfile);
    }

  if(poly!=NULL){
    status= fe_topo02(input, output, poly, mesh, zmin);
    }

  free(topo);

end: __OUT_BASE_LINE__("end of mesh-zmin ... \n");
  exit(0);

error:
  __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}
