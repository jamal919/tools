
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
 

#include "tools-structures.h"

#include "geo.h"
#include "fe.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "parallel.h"

#include "zapper.h"     /*  rutin.h contains common utility routines  */


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t *fe_partition_01_(mesh_t mesh,int *partition,int npartition)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Partitioned meshes are directly taken from metis element list.

  Not suitable for parallel computing (no meshes overlapping).

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int   k,l,m,i,j,n,p,status;
  int   chk=0,count=0,nndes=0;
  int   ninternal,nexternal,keep;
  int   *used;
  mesh_t *work;
  char output[1024];
  bool verbose=0;

  if(verbose) printf ("number of nodes    (original mesh): %6d\n",mesh.nvtxs);
  if(verbose) printf ("number of elements (original mesh): %6d\n",mesh.ntriangles);

  work=new mesh_t[npartition];

  for(p=0;p<npartition;p++) {
    ninternal=0;
    nexternal=0;

    for (m=0;m<mesh.ntriangles;m++) {
      keep=(partition[m]==p);
      switch(keep) {
        case 0:
          nexternal++;
          break;
        case 1:
          ninternal++;
          break;
        }
      }

    work[p].triangles=new triangle_t [ninternal];

    work[p].ntriangles=ninternal;

    ninternal=0;
    nexternal=0;

/* *-----------------------------------------------------------------------------
    screen elements*/
    for (m=0;m<mesh.ntriangles;m++) {
      keep=(partition[m]==p);
      switch(keep) {
        case 0:
          nexternal++;
          break;
        case 1:
          work[p].triangles[ninternal].origin=m;
          ninternal++;
          break;
        }
      }

/* *-----------------------------------------------------------------------------
    count vertices*/
    used=new int[mesh.nvtxs];
    nndes=0;
    for (n=0;n<mesh.nvtxs;n++) used[n]=-1;
    for (m=0;m<work[p].ntriangles;m++) {
      for(i=0;i<3;i++) {
        n=mesh.triangles[work[p].triangles[m].origin].vertex[i];
        if(used[n]==-1) {
          used[n]=nndes;
          nndes++;
          }
        }
      }
/* *-----------------------------------------------------------------------------
    build vertices and elements table*/
    work[p].nvtxs = nndes;
    work[p].vertices=new vertex_t[work[p].nvtxs];
    for (n=0; n<work[p].nvtxs; n++) work[p].vertices[n].null_value();
    for (n=0;n<mesh.nvtxs;n++) {
      if(used[n]!=-1) {
        work[p].vertices[used[n]].lon=mesh.vertices[n].lon;
        work[p].vertices[used[n]].lat=mesh.vertices[n].lat;
        work[p].vertices[used[n]].code=0;
        }
      }
    for (m=0;m<work[p].ntriangles;m++) {
      for(i=0;i<3;i++) {
        n=mesh.triangles[work[p].triangles[m].origin].vertex[i];
        if(used[n]!=-1) {
          work[p].triangles[m].vertex[i]=used[n];
          }
        else {
          printf("triangle %d vertex %d is not registred\n",m,i);
          }
        }
      }
    if(verbose) printf ("number of nodes    (splitted mesh %d): %6d\n",p,work[p].nvtxs);
    if(verbose) printf ("number of elements (splitted mesh %d): %6d\n",p,work[p].ntriangles);
//    work[p].degree = 1;
    status= fe_e2n(&(work[p]));
//    status= build_edgetable(&work[p],0);
//    status= fe_codetable2(&work[p],0);
    if(verbose) {
      sprintf(output,"%s-raw.%2.2d.nei","internal",p);
      status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID,work[p]);
      }
    }

  delete[] used;

  if(verbose) printf("split completed\n");

  return(work);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nndes,status,*elts=NULL,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int option;
  FILE *file=NULL;
  FILE *out=NULL;
  char *keyword=NULL,*zone=NULL,*gridfile=NULL;
  char *meshfile=NULL,*depthfile=NULL,*output=NULL,*poly=NULL;
  char cmdline[1024];
  mesh_t mesh,refined,*splitted=NULL,*final=NULL;
  int  *partition=NULL, npartition;

  fct_echo(argc,argv);
 
  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

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

  status= fe_codetable2(&mesh,0,1,0);
  if(status!=0) {
    printf("unable to rebuild the limits table and codes of the original mesh\n");
    goto error;
    }

/*-----------------------------------------------------------------------------
  extract submesh from selection to allow decent cartesian reshape */
  status=fe_savemesh("mesh.metis",MESH_FILE_FORMAT_METIS,mesh);

  npartition=4;
  sprintf(cmdline,"%s %s %d","/home/softs/metis-4.0/partnmesh","mesh.metis",npartition);
  status=system(cmdline);

//   partition=fe_read_epartition("mesh.metis.epart.4",mesh);
//   splitted=fe_partition_01_(mesh,partition,npartition);

  partition=fe_read_epartition("mesh.metis.epart.4",mesh);
  splitted=fe_partition_03(mesh,partition,npartition);

  __ERR_BASE_LINE__("exiting\n");exit(0);

end: __OUT_BASE_LINE__("end of refine ... \n");
  exit(0);

error:
  __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}
