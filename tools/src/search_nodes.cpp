#define MAIN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"
#include "constants.h"

#include "polygones.h"
#include "fe.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_selectnodes(mesh_t mesh, char *poly, int *selected)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,m,n,n1,n2,status;
  int   count=0;
  int   nelts,nndes,nedges;
  int   inside[2];
  edge_t *edges;
  triangle_t *elt,*ptr;
  double t,p,d,t1,t2,p1,p2;
  plg_t *polygones=NULL;
  int npolygones=0;

  elt=mesh.triangles;

  nelts =mesh.ntriangles;
  nndes =mesh.nvtxs;
  nedges=mesh.nedges;
  edges =mesh.edges;

/*-----------------------------------------------------------------------------
  area criterion */
  if(poly!=NULL) status=plg_load_scan(poly, &polygones, &npolygones);
  else goto end;

  for (n=0;n<nndes;n++) {
    t=mesh.vertices[n].lon;
    p=mesh.vertices[n].lat;
    selected[n]=plg_TestInterior(t,p,polygones,npolygones);
    if(selected[n]==1) count++;
    }

 end:
  printf ("number of nodes selected %d\n",count);
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double lat,lon;
  int rep;
  int nndes,status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int option;
  FILE *file;
  FILE *out;
  char *keyword,*zone=NULL,*gridfile=NULL;
  char *meshfile=NULL,*depthfile=NULL,*output=NULL,*poly=NULL;
  mesh_t mesh;
  double dmax=0;
  int *selected;

  fct_echo(argc,argv);
 
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
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
          n++;
        break;
      }
      free(keyword);
    }

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

 selected=(int *)malloc(mesh.nvtxs*sizeof(int));

 status=fe_selectnodes( mesh, poly, selected);

 error:
 __ERR_BASE_LINE__("exiting\n");exit(-1);
}
